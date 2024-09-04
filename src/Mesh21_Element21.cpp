/**
 * @file Mesh21_Element21.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Mesh21/ELement21.hpp implementation.
 * @date 2024-07-25
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Ivo.hpp>

namespace ivo {

    // Constructors.
    
    /**
     * @brief Full constructor.
     * 
     * @param base Base.
     * @param height Height.
     * @param p Space degree.
     * @param q Time degree.
     */
    Element21::Element21(const Polygon21 &base, const Real &height, const Natural &p, const Natural &q): _base{base}, _height{height}, _p{p}, _q{q} {
        #ifndef NDEBUG // Integrity check.
        assert(height > ___zero);
        #endif
    }

    /**
     * @brief Low-degree constructor.
     * 
     * @param base Base.
     * @param height Height.
     */
    Element21::Element21(const Polygon21 &base, const Real &height): _base{base}, _height{height}, _p{1}, _q{1} {
        #ifndef NDEBUG // Integrity check.
        assert(height > ___zero);
        #endif
    }

    /**
     * @brief Copy constructor.
     * 
     * @param element Element.
     */
    Element21::Element21(const Element21 &element): _base{element._base}, _height{element._height}, _p{element._p}, _q{element._q} {}

    // Methods.

    /**
     * @brief Element's dofs.
     * 
     * @return Natural 
     */
    Natural Element21::dofs() const { return (this->_q + 1) * (this->_p + 1) * (this->_p + 2) / 2; } // [?]
    
    /**
     * @brief Element's bottom base.
     * 
     * @return Polygon21 
     */
    Polygon21 Element21::b_base() const { return this->_base; }

    /**
     * @brief Element's top base.
     * 
     * @return Polygon21 
     */
    Polygon21 Element21::t_base() const {
        std::vector<Point21> t_points;

        for(const auto &b_point: this->_base.points())
            t_points.emplace_back(b_point + this->_height * 1.0_t);

        return Polygon21{t_points};
    }

    /**
     * @brief Element's bottom edges.
     * 
     * @return std::vector<Edge21> 
     */
    std::vector<Edge21> Element21::b_edges() const { return this->_base.edges(); }

    /**
     * @brief Element's top edges.
     * 
     * @return std::vector<Edge21> 
     */
    std::vector<Edge21> Element21::t_edges() const {
        std::vector<Edge21> t_edges;

        for(const auto &edge: this->_base.edges())
            t_edges.emplace_back(edge(0) + this->_height * 1.0_t, edge(1) + this->_height * 1.0_t);

        return t_edges;
    }

    /**
     * @brief Element's faces.
     * 
     * @return std::vector<Polygon21> 
     */
    std::vector<Polygon21> Element21::faces() const {
        std::vector<Polygon21> faces;

        for(const auto &edge: this->_base.edges())
            faces.emplace_back(Polygon21{{edge(0), edge(1), edge(1) + this->_height * 1.0_t, edge(0) + this->_height * 1.0_t}});

        return faces;
    }

    /**
     * @brief Element's time interval.
     * 
     * @return std::array<Real, 2> 
     */
    std::array<Real, 2> Element21::interval() const {
        Real start = this->_base.points()[0](2);
        return {start, start + this->_height};
    }

    // Output.
    
    /**
     * @brief Element output. CSV format.
     * 
     * @param ost 
     * @param element Element.
     * @return std::ostream& 
     */
    std::ostream &operator <<(std::ostream &ost, const Element21 &element) {
        for(const auto &point: element._base.points())
            ost << point(0) << "," << point(1) << "," << point(2) << ",";

        return ost << element._height << std::flush;
    }

}