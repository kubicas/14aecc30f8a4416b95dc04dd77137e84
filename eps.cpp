#define EPS
#include "eps/eps.h"
#include <fstream>
#include <stdexcept>
#include <sstream>
#include <vector>
#include <limits>
#include <algorithm>
#include <stack>
#include <regex>
#include <set>
#include <iostream>

namespace // anonymous
{

void calculate_arc(
    eps::point_t const C,
    eps::vect_t rx,
    eps::vect_t& ry,
    eps::point_t& B,
    eps::point_t& E,
    eps::transformation_t& transformation,
    bool& is_ellipse,
    float const epsilon)
{
    float ax = abs(rx);
    float ay = abs(ry);
    if (abs(ax - ay) <= epsilon)
    {
        is_ellipse = false;
        return;
    }
    rx *= 1.0f / ax;
    transformation.m_r.m_x.m_x = +rx.m_x; transformation.m_r.m_x.m_y = +rx.m_y;
    transformation.m_r.m_y.m_x = -rx.m_y; transformation.m_r.m_y.m_y = +rx.m_x;
    transformation.m_r.m_x *= ax;
    transformation.m_r.m_y *= ay;
    transformation.m_t = C;
    eps::transformation_t transformation_inv = ~transformation;
    B *= transformation_inv;
    E *= transformation_inv;
    ry *= abs(E);
    is_ellipse = true;
}

void bezier_points_of_tangency(
    eps::point_t a, eps::point_t ai, eps::point_t bi, eps::point_t b, float epsilon,
    float& minx, float& maxx, float& miny, float& maxy)
{
    float t1 = 0;
    float t2 = 0;
    int dummy;
    eps::abc_formula(-a.m_x + 3 * ai.m_x - 3 * bi.m_x + b.m_x, 2 * a.m_x - 4 * ai.m_x + 2 * bi.m_x, ai.m_x - a.m_x, epsilon, dummy, t1, t2);
    t1 = eps::clip(t1, 0.f, 1.f);
    t2 = eps::clip(t2, 0.f, 1.f);
    float v1 = eps::bezier(a.m_x, ai.m_x, bi.m_x, b.m_x, t1);
    float v2 = eps::bezier(a.m_x, ai.m_x, bi.m_x, b.m_x, t2);
    minx = std::min(minx, std::min(v1, v2));
    maxx = std::max(maxx, std::max(v1, v2));
    t1 = 0;
    t2 = 0;
    eps::abc_formula(-a.m_y + 3 * ai.m_y - 3 * bi.m_y + b.m_y, 2 * a.m_y - 4 * ai.m_y + 2 * bi.m_y, ai.m_y - a.m_y, epsilon, dummy, t1, t2);
    t1 = eps::clip(t1, 0.f, 1.f);
    t2 = eps::clip(t2, 0.f, 1.f);
    v1 = eps::bezier(a.m_y, ai.m_y, bi.m_y, b.m_y, t1);
    v2 = eps::bezier(a.m_y, ai.m_y, bi.m_y, b.m_y, t2);
    miny = std::min(miny, std::min(v1, v2));
    maxy = std::max(maxy, std::max(v1, v2));
}

void ellipse_points_of_tangency(
    eps::point_t a, eps::point_t c, eps::vect_t rx, eps::vect_t ry, float epsilon,
    float& minx, float& maxx, float& miny, float& maxy)
{
    float xx2 = rx.m_x*rx.m_x;
    float xy2 = rx.m_y*rx.m_y;
    float yx2 = ry.m_x*ry.m_x;
    float yy2 = ry.m_y*ry.m_y;
    float ax2 = xx2 + xy2;
    float ay2 = yx2 + yy2;
    float dx = std::sqrt((ax2*xx2 + ay2*xy2) / ax2);
    float dy = std::sqrt((ax2*yx2 + ay2*yy2) / ay2);
    maxx = std::max(maxx, c.m_x + dx);
    maxy = std::max(maxy, c.m_y + dy);
    minx = std::min(minx, c.m_x - dx);
    miny = std::min(miny, c.m_y - dy);
}

void arc_points_of_tangency(
    eps::point_t a, eps::point_t c, eps::vect_t rx, eps::vect_t ry, eps::point_t b, float epsilon,
    float& minx, float& maxx, float& miny, float& maxy)
{
    float xx2 = rx.m_x*rx.m_x;
    float xy2 = rx.m_y*rx.m_y;
    float yx2 = ry.m_x*ry.m_x;
    float yy2 = ry.m_y*ry.m_y;
    float ax2 = xx2 + xy2;
    float ay2 = yx2 + yy2;
    float dx = std::sqrt((ax2*xx2 + ay2 * xy2) / ax2);
    float dy = std::sqrt((ax2*yx2 + ay2 * yy2) / ay2);
    bool positive = (rx.m_x * ry.m_y - rx.m_y * ry.m_x) >= 0;
    eps::vect_t ra = a - c;
    eps::vect_t rb = b - c;
    float aa = std::atan2(ra.m_y, ra.m_x);
    float ab = std::atan2(rb.m_y, rb.m_x);
    if (eps::in_arc(0.f, aa, ab, positive))
    {
        maxx = std::max(maxx, c.m_x + dx);
    }
    if (eps::in_arc(static_cast<float>(0.5*eps::pi), aa, ab, positive))
    {
        maxy = std::max(maxy, c.m_y + dy);
    }
    if (eps::in_arc(static_cast<float>(eps::pi), aa, ab, positive))
    {
        minx = std::min(minx, c.m_x - dx);
    }
    if (eps::in_arc(static_cast<float>(1.5*eps::pi), aa, ab, positive))
    {
        miny = std::min(miny, c.m_y - dy);
    }
    if (eps::in_arc(static_cast<float>(2 * eps::pi), aa, ab, positive))
    {
        maxx = std::max(maxx, c.m_x + dx);
    }
}

}; // anonymous

namespace eps
{

std::ostream& operator<<(std::ostream& o, eps::vect_t v)
{
    return o << v.m_x << ' ' << v.m_y;
}

std::ostream& operator<<(std::ostream& o, eps::point_t v)
{
    return o << v.m_x << ' ' << v.m_y;
}

std::ostream& operator<<(std::ostream& o, eps::area_t a)
{
    return o << a.m_min << ' ' << a.m_max;
}

std::ostream& operator<<(std::ostream& o, eps::rotation_t r)
{
    return o << r.m_x << ' ' << r.m_y;
}

std::ostream& operator<<(std::ostream& o, eps::transformation_t t)
{
    return o << "[ " << t.m_r << ' ' << t.m_t << " ]";
}

area_t null_bounding_box()
{
    return area_t(
        point_t(
            std::numeric_limits<float>::max(),
            std::numeric_limits<float>::max()),
        point_t(
            std::numeric_limits<float>::min(),
            std::numeric_limits<float>::min()));
}

struct lineending_none_t
    : lineending_t
{
    void draw_procedure(std::ostream&) const override {}
    void draw(std::ostream&, eps::point_t, float, float,
        float, float, float) const override {}
} l_lineending_none;
lineending_t* lineending_none() { return &l_lineending_none; }

struct linestyle_none_t
    : linestyle_t
{
    void draw(std::ostream& stream) const override
    {
        stream << "[] 0 setdash\n";
    }
} l_linestyle_none;
linestyle_t* linestyle_none() { return &l_linestyle_none; }

float get_epsilon(graphicsstate_t& graphicsstate)
{
    return graphicsstate.epsilon();
}

void new_path(std::ostream& stream)
{
    stream << "newpath\n";
}

void closepath(std::ostream& stream)
{
    stream << "closepath\n";
}

void stroke(std::ostream& stream, graphicsstate_t& graphicsstate, iproperties_t const& properties)
{
    if (properties.linewidth() != graphicsstate.linewidth())
    {
        graphicsstate.setlinewidth(properties.linewidth());
        stream << graphicsstate.linewidth() << " setlinewidth\n";
    }
    if ((properties.linercolor() != graphicsstate.linercolor()) ||
        (properties.linegcolor() != graphicsstate.linegcolor()) ||
        (properties.linebcolor() != graphicsstate.linebcolor()))
    {
        graphicsstate.setlinergbcolor(
            properties.linercolor(), properties.linegcolor(), properties.linebcolor());
        if ((graphicsstate.linercolor() == graphicsstate.linegcolor()) && (graphicsstate.linegcolor() == graphicsstate.linebcolor()))
        {
            stream << graphicsstate.linercolor() << " setgray\n";
        }
        else
        {
            stream << graphicsstate.linercolor() << ' ' << graphicsstate.linegcolor() << ' ' << graphicsstate.linebcolor() << ' ' << " setrgbcolor\n";
        }
        graphicsstate.setfillrgbcolor(properties.linercolor(), properties.linegcolor(), properties.linebcolor());
    }
    if (properties.linestyle() != graphicsstate.linestyle())
    {
        graphicsstate.setlinestyle(properties.linestyle());
        graphicsstate.linestyle()->draw(stream);
    }
    if (properties.linecap() != graphicsstate.linecap())
    {
        graphicsstate.setlinecap(properties.linecap());
        stream << static_cast<int>(graphicsstate.linecap()) << " setlinecap\n";
    }
    if (properties.linejoin() != graphicsstate.linejoin())
    {
        graphicsstate.setlinejoin(properties.linejoin());
        stream << static_cast<int>(graphicsstate.linejoin()) << " setlinejoin\n";
    }
    if (properties.miterlimit() != graphicsstate.miterlimit())
    {
        graphicsstate.setmiterlimit(properties.miterlimit());
        stream << graphicsstate.miterlimit() << " setmiterlimit\n";
    }
    stream << "stroke\n";
}

void fill(std::ostream& stream, graphicsstate_t& graphicsstate, iproperties_t const& properties, bool and_stroke) // clears moveto data!!
{
    if ((properties.fillrcolor() != graphicsstate.fillrcolor()) ||
        (properties.fillgcolor() != graphicsstate.fillgcolor()) ||
        (properties.fillbcolor() != graphicsstate.fillbcolor()))
    {
        graphicsstate.setfillrgbcolor(
            properties.fillrcolor(), properties.fillgcolor(), properties.fillbcolor());
        if ((graphicsstate.fillrcolor() == graphicsstate.fillgcolor()) && (graphicsstate.fillgcolor() == graphicsstate.fillbcolor()))
        {
            stream << graphicsstate.fillrcolor() << " setgray\n";
        }
        else
        { 
            stream << graphicsstate.fillrcolor() << ' ' << graphicsstate.fillgcolor() << ' ' << graphicsstate.fillbcolor() << ' ' << " setrgbcolor\n";
        }
        graphicsstate.setlinergbcolor(properties.fillrcolor(), properties.fillgcolor(), properties.fillbcolor());
    }
    if (and_stroke)
    {
        stream << "gsave\n";
        stream << "fill\n";
        stream << "grestore\n";
        stroke(stream, graphicsstate, properties);
    }
    else
    {
        stream << "fill\n";
    }
}

void gsave(std::ostream& stream)
{
    stream << "gsave\n";
}

void grestore(std::ostream& stream)
{
    stream << "grestore\n";
}

void initgraphics(std::ostream& stream, graphicsstate_t& graphicsstate)
{
    stream << "initgraphics\n";
    graphicsstate = graphicsstate_t();
}

void moveto(std::ostream& stream, point_t p)
{
    stream << p << " moveto\n";
}

void rmoveto(std::ostream& stream, vect_t v)
{
    stream << v << " rmoveto\n";
}

void lineto(std::ostream& stream, point_t p)
{
    stream << p << " lineto\n";
}

void rlineto(std::ostream& stream, vect_t v)
{
    stream << v << " rlineto\n";
}

void curveto(std::ostream& stream, eps::point_t tangent1, eps::point_t tangent2, eps::point_t end)
{
    stream << tangent1 << ' ' << tangent2 << ' ' << end << " curveto\n";
}

void rcurveto(std::ostream& stream, eps::vect_t tangent1, eps::vect_t tangent2, eps::vect_t end)
{
    stream << tangent1 << ' ' << tangent2 << ' ' << end << " rcurveto\n";
}

void arc(std::ostream& stream, eps::point_t center, float radius, float begin_angle, float end_angle)
{
    stream << center << ' ' << radius << ' ' << begin_angle << ' ' << end_angle << " arc\n";
}

void arcn(std::ostream& stream, eps::point_t center, float radius, float begin_angle, float end_angle)
{
    stream << center << ' ' << radius << ' ' << begin_angle << ' ' << end_angle << " arcn\n";
}

void arct(std::ostream& stream, eps::point_t tangent, eps::point_t end, float radius)
{
    stream << tangent << ' ' << end << ' ' << radius << " arct\n";
}

void show(std::ostream& stream, std::string const& text)
{
    std::string t(text);
	t = std::regex_replace(t, std::regex("\\("), "\\(");
	t = std::regex_replace(t, std::regex("\\)"), "\\)");
    stream << "(" << t << ") show\n";
}

// see http://texdoc.net/texmf-dist/doc/latex/psfrag/pfgguide.pdf
void showlatex(std::ostream& stream, std::string const& text, text_ref_t text_ref, float scale, float rotate/*deg*/)
{
    std::string t(text);
    t = std::regex_replace(t, std::regex("\\("), "\\(");
    t = std::regex_replace(t, std::regex("\\)"), "\\)");
    t = std::regex_replace(t, std::regex("\\\\"), "\\\\");
    static char const* to_text_ref[static_cast<int>(text_ref_t::number_of_text_refs)] =
    { "bl][Bl", "bc][Bl", "br][Bl", "cl][Bl", "cc][Bl", "cr][Bl", "tl][Bl", "tc][Bl", "tr][Bl", "Bl][Bl", "Bc][Bl", "Br][Bl" };
    stream << "(\\\\tex[" << to_text_ref[static_cast<int>(text_ref)] << "][" << scale << "][" << rotate << "]{" << t << "}) show\n";
}

void clip(std::ostream& stream)
{
    stream << "clip\n";
}

void pushmatrix(std::ostream& stream)
{
    stream << "matrix currentmatrix\n";
}

void concatmatrix(std::ostream& stream, eps::transformation_t const& transformation)
{
    stream << transformation << " concat\n";
}

void popmatrix(std::ostream& stream)
{
    stream << "setmatrix\n";
}

void scale(std::ostream& stream, float x, float y)
{
    stream << x << ' ' << y << " scale\n";
}

void rotate(std::ostream& stream, float angle)
{
    stream << angle << " rotate\n";
}

void concat(std::ostream& stream, transformation_t t)
{
    stream << t << " concat\n";
}

void begin_procedure(std::ostream& stream, std::string const& name, std::vector<char const*> l)
{
    stream << "/" << name << " {\n";
    for (std::vector<char const*>::const_reverse_iterator it = l.crbegin(); it != l.crend(); ++it)
    {
        stream << "/" << *it << " exch def\n";
    }
    gsave(stream);
    stream << "matrix currentmatrix currentlinecap currentlinejoin currentmiterlimit initgraphics setmiterlimit setlinejoin setlinecap setmatrix\n";
}

void end_procedure(std::ostream& stream)
{
    eps::grestore(stream);
    stream << "} bind def\n";
}

void call_procedure(std::ostream& stream, std::string const& name)
{
    stream << name << "\n";
}

void draw_ellipse(std::ostream& stream, point_t a, point_t c, vect_t rx, vect_t ry, float epsilon)
{
    transformation_t transformation;
    bool is_ellipse;
    point_t b(a);
    calculate_arc(c, rx, ry, a, b, transformation, is_ellipse, epsilon);
    if (is_ellipse)
    {
        pushmatrix(stream);
        concatmatrix(stream, transformation);
        float a1 = to_deg(std::atan2(a.m_y, a.m_x));
        float a2 = a1 + 360;
        arc(stream, point_t(0.f, 0.f), 1.f, a1, a2);
        popmatrix(stream);
    }
    else
    {
        float r = abs(ry);
        vect_t vB = a - c;
        vect_t vE = b - c;
        float a1 = to_deg(std::atan2(vB.m_y, vB.m_x));
        float a2 = a1 + 360;
        arc(stream, c, r, a1, a2);
    }
}

void draw_arc(std::ostream& stream, point_t a, point_t c, vect_t rx, vect_t ry, point_t b, float epsilon)
{
    transformation_t transformation;
    bool is_ellipse;
    calculate_arc(c, rx, ry, a, b, transformation, is_ellipse, epsilon);
    bool positive = (rx.m_x * ry.m_y - rx.m_y * ry.m_x) >= 0;
    if (is_ellipse)
    {
        pushmatrix(stream);
        concatmatrix(stream, transformation);
        float a1 = to_deg(std::atan2(a.m_y, a.m_x));
        float a2 = to_deg(std::atan2(b.m_y, b.m_x));
        if (positive)
        {
            arc(stream, point_t(0.f, 0.f), 1.f, a1, a2);
        }
        else
        {
            arcn(stream, point_t(0.f, 0.f), 1.f, a1, a2);
        }
        popmatrix(stream);
    }
    else
    {
        float r = abs(ry);
        vect_t vB = a - c;
        vect_t vE = b - c;
        float a1 = to_deg(std::atan2(vB.m_y, vB.m_x));
        float a2 = to_deg(std::atan2(vE.m_y, vE.m_x));
        if (positive)
        {
            arc(stream, c, r, a1, a2);
        }
        else
        {
            arcn(stream, c, r, a1, a2);
        }
    }
}

area_t bezier_bounding_box(point_t a, point_t ai, point_t bi, point_t b, float epsilon)
{
    float minx = std::min(a.m_x, b.m_x);
    float maxx = std::max(a.m_x, b.m_x);
    float miny = std::min(a.m_y, b.m_y);
    float maxy = std::max(a.m_y, b.m_y);
    bezier_points_of_tangency(a, ai, bi, b, epsilon, minx, maxx, miny, maxy);
    return area_t(point_t(minx, miny), point_t(maxx, maxy));
}

area_t ellipse_bounding_box(point_t a, point_t c, vect_t rx, vect_t ry, float epsilon)
{
    eps::area_t area = null_bounding_box();
    ellipse_points_of_tangency(a, c, rx, ry, epsilon, area.m_min.m_x, area.m_max.m_x, area.m_min.m_y, area.m_max.m_y);
    return area;
}

area_t arc_bounding_box(point_t a, point_t c, vect_t rx, vect_t ry, point_t b, float epsilon)
{
    float minx = std::min(a.m_x, b.m_x);
    float maxx = std::max(a.m_x, b.m_x);
    float miny = std::min(a.m_y, b.m_y);
    float maxy = std::max(a.m_y, b.m_y);
    arc_points_of_tangency(a, c, rx, ry, b, epsilon, minx, maxx, miny, maxy);
    return area_t(point_t(minx, miny), point_t(maxx, maxy));
}

std::string single_line_it(std::string const& in)
{
    std::stringstream ss;
    ss  << R"(\tiny\textit{)"
        << in
        << '}';
    return ss.str();
}

std::string single_line_listing(std::string const& in, char delim_char)
{
    std::stringstream ss;
    ss << R"(\lstinline[basicstyle=\tiny])"
        << delim_char
        << in
        << delim_char;
    return ss.str();
}

std::string multi_line_listing(float width, std::string const& in, char delim_char)
{
    std::stringstream ss;
    ss  << R"(\raisebox{\dimexpr\depth-\fontchardp\font`y}{\parbox{)"
        << width
        << R"(bp}{\centering \linespread{0.3}\selectfont\lstinline[basicstyle=\tiny])"
        << delim_char
        << in
        << delim_char
        << R"(}})";
        return ss.str();
}

/////////////////////////////////////////////

using properties_mem_mgr_t = std::set<eps::properties_override_t>;
static properties_mem_mgr_t properties_mem_mgr;
properties_mem_mgr_t& get_mem_mgr()
{
    return properties_mem_mgr;
}

shape_t::shape_t(iproperties_t const& parent_properties)
    : m_parent_properties(parent_properties)
    , m_pproperties_override(nullptr)
{}
shape_t::shape_t(shape_t const& other)
    : m_parent_properties(other.m_parent_properties)
    , m_pproperties_override(other.m_pproperties_override)
{
    inc_ref();
}
shape_t::~shape_t()
{
    dec_ref();
}

#define SET_FUNCTION_1(TYPE, NAME) \
void shape_t::set##NAME(TYPE NAME) \
{ \
    properties_override_t properties_override; \
    get(properties_override); \
    properties_override.set##NAME(NAME); \
    add(properties_override); \
}

#define SET_FUNCTION_1B(TYPE, NAME, NAMEB) \
void shape_t::set##NAME(TYPE NAME) \
{ \
    properties_override_t properties_override; \
    get(properties_override); \
    properties_override.set##NAMEB(NAME, NAME, NAME); \
    add(properties_override); \
}

#define SET_FUNCTION_3(TYPE, NAME, NAME1, NAME2, NAME3) \
void shape_t::set##NAME(TYPE NAME1, TYPE NAME2, TYPE NAME3) \
{ \
    properties_override_t properties_override; \
    get(properties_override); \
    properties_override.set##NAME(NAME1, NAME2, NAME3); \
    add(properties_override); \
}
SET_FUNCTION_1(float, linewidth)
SET_FUNCTION_1B(float, linegray, linergbcolor)
SET_FUNCTION_3(float, linergbcolor, linercolor, linegcolor, linebcolor)
SET_FUNCTION_1B(float, fillgray, fillrgbcolor)
SET_FUNCTION_3(float, fillrgbcolor, fillrcolor, fillgcolor, fillbcolor)
SET_FUNCTION_1(cap_t, linecap)
SET_FUNCTION_1(join_t, linejoin)
SET_FUNCTION_1(float, miterlimit)
SET_FUNCTION_1(float, epsilon)
SET_FUNCTION_1(lineending_t*, lineend)
SET_FUNCTION_1(lineending_t*, linebegin)
SET_FUNCTION_1(linestyle_t*, linestyle)
#undef SET_FUNCTION_1
#undef SET_FUNCTION_1B
#undef SET_FUNCTION_3

#define GET_FUNCTION_1(TYPE, NAME) \
TYPE shape_t::NAME() const \
{ \
    return m_pproperties_override ? \
    m_pproperties_override->NAME(m_parent_properties.NAME()) : \
    m_parent_properties.NAME(); \
}
GET_FUNCTION_1(float, linewidth)
GET_FUNCTION_1(float, linercolor)
GET_FUNCTION_1(float, linegcolor)
GET_FUNCTION_1(float, linebcolor)
GET_FUNCTION_1(float, fillrcolor)
GET_FUNCTION_1(float, fillgcolor)
GET_FUNCTION_1(float, fillbcolor)
GET_FUNCTION_1(cap_t, linecap)
GET_FUNCTION_1(join_t, linejoin)
GET_FUNCTION_1(float, miterlimit)
GET_FUNCTION_1(float, epsilon)
GET_FUNCTION_1(lineending_t*, lineend)
GET_FUNCTION_1(lineending_t*, linebegin)
GET_FUNCTION_1(linestyle_t*, linestyle)
#undef GET_FUNCTION_1

void shape_t::inc_ref()
{
    if (!m_pproperties_override)
    {
        return;
    }
    ++m_pproperties_override->m_ref_count;
}
void shape_t::dec_ref()
{
    if (!m_pproperties_override)
    {
        return;
    }
    if (m_pproperties_override->m_ref_count > 0)
    {
        --m_pproperties_override->m_ref_count;
    }
    if (!m_pproperties_override->m_ref_count)
    {
        get_mem_mgr().erase(*m_pproperties_override);
    }
    m_pproperties_override = nullptr;
}
void shape_t::get(properties_override_t& properties_override)
{
    if (!m_pproperties_override)
    {
        return;
    }
    properties_override = *m_pproperties_override;
    dec_ref();
}
void shape_t::add(properties_override_t& properties_override)
{
    std::pair<properties_mem_mgr_t::iterator, bool> ret =
        get_mem_mgr().insert(properties_override);
    m_pproperties_override = &*(ret.first);
    inc_ref();
}


void group_t::add(std::unique_ptr<shape_t>&& o)
{
    m_shapes.push_back(std::move(o));
}


area_t group_t::bounding_box(float epsilon)
{
    eps::area_t area = null_bounding_box();
    for (std::unique_ptr<shape_t> &i : m_shapes)
    {
        eps::area_t shape_area = i->bounding_box(epsilon);
        min_bounding_box(area.m_min, shape_area.m_min);
        max_bounding_box(area.m_max, shape_area.m_max);
    }
    return area;
}

void group_t::draw(std::ostream& stream, eps::graphicsstate_t& graphicsstate) const
{
    for (std::unique_ptr<shape_t> const &i : m_shapes)
    {
        i->draw(stream, graphicsstate);
    }
}

void group_t::apply(transformation_t const& t, bool excluding_text)
{
    for (std::unique_ptr<shape_t> &i : m_shapes)
    {
        i->apply(t, excluding_text);
    }
}

struct canvas_impl_t
    : public eps::canvas_t
{
public:
    canvas_impl_t(graphicsstate_t const& root_properties, std::string const& filename)
        : eps::canvas_t(root_properties)
        , m_ofs(filename, std::ofstream::out)
    {
        if (!m_ofs.is_open())
        {
            THROW(std::runtime_error, "E0001", << "Cannot open'" << filename << "'");
        }
    }
    void draw() override
    {
        eps::graphicsstate_t graphicsstate;
        area_t area = bounding_box(graphicsstate.epsilon());
        area.m_min.m_x = std::floor(area.m_min.m_x);
        area.m_min.m_y = std::floor(area.m_min.m_y);
        area.m_max.m_x = std::ceil(area.m_max.m_x);
        area.m_max.m_y = std::ceil(area.m_max.m_y);
        m_ofs << "%!PS-Adobe-3.0\n" << "%%BoundingBox: " << area << std::endl;
        m_ofs << "/Times-Roman 10 selectfont\n"; // select one font so that psfrag works
        for (eps::properties_override_t const& p : properties_mem_mgr)
        {
            if (p.lineend(nullptr) && p.lineend(nullptr)->m_used)
            {
                p.lineend(nullptr)->m_used = false;
                p.lineend(nullptr)->draw_procedure(m_ofs);
            }
            if (p.linebegin(nullptr) && p.linebegin(nullptr)->m_used)
            {
                p.linebegin(nullptr)->m_used = false;
                p.linebegin(nullptr)->draw_procedure(m_ofs);
            }
        }
        group_t::draw(m_ofs, graphicsstate);
    }
    std::ofstream m_ofs;
};

static graphicsstate_t const root_properties;

std::unique_ptr<canvas_t> create_canvas(
    std::string const& filename)
{
    return std::make_unique<eps::canvas_impl_t>(root_properties, filename);
}

EPS_API void handle_exception()
{
    try
    {
        throw;
    }
    catch (std::exception& e)
    {
        std::cerr << "error " << e.what() << std::endl;
    }
    catch (...)
    {
        std::cerr << "error non-standard exception" << std::endl;
    }
}

}; // namespace eps