#define EPS
#include "eps/eps_basic_shapes.h"

#include <iterator>

namespace // anonymous
{

class beginpoint_t
    : public eps::section_t
{
public:
    beginpoint_t(int a) : m_a(a) {}
    int m_a;
};

class linesection_t
    : public eps::section_t
{
public:
    linesection_t(int b) : m_b(b) {}
    int m_b;
};

class beziersection_t
    : public eps::section_t
{
public:
    beziersection_t(int a, int ai, int bi, int b) : m_a(a), m_ai(ai), m_bi(bi), m_b(b) {}
    int m_a;
    int m_ai;
    int m_bi;
    int m_b;
};

class arcsection_t
    : public eps::section_t
{
public:
    arcsection_t(int a, int c, int ax, int ay, int b) : m_a(a), m_c(c), m_ax(ax), m_ay(ay), m_b(b) {}
    int m_a;
    int m_c;
    int m_ax;
    int m_ay;
    int m_b;
};

class closepath_t
    : public eps::section_t
{};

}; // namespace anonymous

namespace eps
{

path_t::path_t(iproperties_t const& parent_properties)
    : dynamic_shape_t(parent_properties)
    , m_fill(false)
{}

path_t::path_t(path_t const& rhs)
    : dynamic_shape_t(rhs)
    , m_fill(rhs.m_fill)
{
    m_sections.reserve(rhs.m_sections.size());
    std::transform(rhs.m_sections.begin(), rhs.m_sections.end(), std::back_inserter(m_sections), 
        [](const std::unique_ptr<section_t>& section) -> std::unique_ptr<section_t>
    {
        if (beginpoint_t* s = dynamic_cast<beginpoint_t*>(&*section))
        {
            return std::make_unique<beginpoint_t>(*s);
        }
        else if (linesection_t* s = dynamic_cast<linesection_t*>(&*section))
        {
            return std::make_unique<linesection_t>(*s);
        }
        else if (beziersection_t* s = dynamic_cast<beziersection_t*>(&*section))
        {
            return std::make_unique<beziersection_t>(*s);
        }
        else if (arcsection_t* s = dynamic_cast<arcsection_t*>(&*section))
        {
            return std::make_unique<arcsection_t>(*s);
        }
        else if (closepath_t* s = dynamic_cast<closepath_t*>(&*section))
        {
            return std::make_unique<closepath_t>(*s);
        }
        else
        {
            THROW(std::logic_error, "E0101", << "Unsupported closepath_t type");
        }
    });
}

void path_t::draw(std::ostream& stream, graphicsstate_t& graphicsstate) const
{
    new_path(stream);
    for (std::unique_ptr<section_t> const& s : m_sections)
    {
        if (beginpoint_t const* p = dynamic_cast<beginpoint_t const*>(s.get()))
        {
            eps::moveto(stream, m_[p->m_a]);
        }
        else if (linesection_t const* p = dynamic_cast<linesection_t const*>(s.get()))
        {
            eps::lineto(stream, m_[p->m_b]);
        }
        else if (beziersection_t const* p = dynamic_cast<beziersection_t const*>(s.get()))
        {
            eps::curveto(stream, m_[p->m_ai], m_[p->m_bi], m_[p->m_b]);
        }
        else if (arcsection_t const* p = dynamic_cast<arcsection_t const*>(s.get()))
        {
            eps::draw_arc(stream, m_[p->m_a], m_[p->m_c], m_[p->m_ax] - m_[p->m_c], m_[p->m_ay] - m_[p->m_c], m_[p->m_b], get_epsilon(graphicsstate));
        }
        else if (closepath_t const* p = dynamic_cast<closepath_t const*>(s.get()))
        {
            eps::closepath(stream);
        }
        else
        {
            THROW(std::logic_error, "E0101", << "Unsupported closepath_t type");
        }
    }
    if (m_fill)
    {
        eps::fill(stream, graphicsstate, *this, true);
    }
    else
    {
        eps::stroke(stream, graphicsstate, *this);
    }
}

area_t path_t::bounding_box(float epsilon)
{
    area_t area = null_bounding_box();
    for (std::unique_ptr<section_t> const& s : m_sections)
    {
        if (beginpoint_t const* p = dynamic_cast<beginpoint_t const*>(s.get()))
        {
            min_bounding_box(area.m_min, m_[p->m_a]);
            max_bounding_box(area.m_max, m_[p->m_a]);
        }
        else if (linesection_t const* p = dynamic_cast<linesection_t const*>(s.get()))
        {
            min_bounding_box(area.m_min, m_[p->m_b]);
            max_bounding_box(area.m_max, m_[p->m_b]);
        }
        else if (beziersection_t const* p = dynamic_cast<beziersection_t const*>(s.get()))
        {
            area_t bb = bezier_bounding_box(m_[p->m_a], m_[p->m_ai], m_[p->m_bi], m_[p->m_b], epsilon);
            min_bounding_box(area.m_min, bb.m_min);
            max_bounding_box(area.m_max, bb.m_max);
        }
        else if (arcsection_t const* p = dynamic_cast<arcsection_t const*>(s.get()))
        {
            area_t bb = arc_bounding_box(m_[p->m_a], m_[p->m_c], m_[p->m_ax] - m_[p->m_c], m_[p->m_ay] - m_[p->m_c], m_[p->m_b], epsilon);
            min_bounding_box(area.m_min, bb.m_min);
            max_bounding_box(area.m_max, bb.m_max);
        }
    }
    return area;
}

void path_t::moveto(point_t p)
{
    m_.emplace_back(p);
    m_sections.emplace_back(std::make_unique<beginpoint_t>(static_cast<int>(m_.size()) - 1));
}

void path_t::lineto(point_t p)
{
    m_.emplace_back(p);
    m_sections.emplace_back(std::make_unique<linesection_t>(static_cast<int>(m_.size())-1));
}

void path_t::curveto(point_t tangent1, point_t tangent2, point_t end)
{
    m_.emplace_back(tangent1);
    m_.emplace_back(tangent2);
    m_.emplace_back(end);
    m_sections.emplace_back(std::make_unique<beziersection_t>(
		static_cast<int>(m_.size()) - 4, 
		static_cast<int>(m_.size()) - 3, 
		static_cast<int>(m_.size()) - 2, 
		static_cast<int>(m_.size()) - 1));
}

void path_t::arcto(point_t center, point_t x_ax, point_t y_ax, point_t end)
{
    m_.emplace_back(center);
    m_.emplace_back(x_ax);
    m_.emplace_back(y_ax);
    m_.emplace_back(end);
    m_sections.emplace_back(std::make_unique<arcsection_t>(
		static_cast<int>(m_.size()) - 5, 
		static_cast<int>(m_.size()) - 4, 
		static_cast<int>(m_.size()) - 3, 
		static_cast<int>(m_.size()) - 2, 
		static_cast<int>(m_.size()) - 1));
}

void path_t::closepath()
{
    m_sections.emplace_back(std::make_unique<closepath_t>());
}

}; // namespace eps