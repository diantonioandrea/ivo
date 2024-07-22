/**
 * @file Geometry_Line21.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Geometry/Line21.hpp implementation.
 * @date 2024-07-22
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Ivo.hpp>

namespace ivo {

    // Constructors and copy operators.

    /**
     * @brief Points constructor.
     * 
     * @param p Point.
     * @param q Point.
     */
    Line21::Line21(const Point21 &p, const Point21 &q): _a{q(0) - p(0)}, _b{q(1) - p(1)}, _c{q(2) - p(2)}, _x0{p(0)}, _y0{p(1)}, _t0{p(2)} {
        #ifndef NDEBUG
        assert(p != q);
        #endif
    }

    /**
     * @brief Edge constructor.
     * 
     * @param edge Edge.
     */
    Line21::Line21(const Edge21 &edge): Line21{edge(0), edge(1)} {}

    /**
     * @brief Copy constructor.
     * 
     * @param line Line.
     */
    Line21::Line21(const Line21 &line): _a{line._a}, _b{line._b}, _c{line._c}, _x0{line._x0}, _y0{line._y0}, _t0{line._t0} {}

    /**
     * @brief Copy operator.
     * 
     * @param line Line.
     * @return Line21& 
     */
    Line21 &Line21::operator =(const Line21 &line) {
        this->_a = line._a;
        this->_b = line._b;
        this->_c = line._c;
        this->_x0 = line._x0;
        this->_y0 = line._y0;
        this->_t0 = line._t0;
        return *this;
    }

    // Subscript operator, legacy scalar access (C++23).

    #if __cplusplus > 202002L

    /**
     * @brief Scalar reference access, legacy.
     * 
     * @param j Coordinate index.
     * @param k Parameter index.
     * @return Real& 
     */
    Real &Line21::operator [](const Natural &j, const Natural &k) {
        #ifndef NDEBUG // Integrity check.
        assert(j <= 2);
        assert(k <= 1);
        #endif

        if(j == 0) {
            if(k == 0)
                return this->_a;

            return this->_x0;
        }

        if(j == 1) {
            if(k == 0)
                return this->_b;

            return this->_y0;
        }

        if(k == 0)
            return this->_c;

        return this->_t0;
    }

    #endif

    // Call operator, subscript behaviour.

    /**
     * @brief Scalar access.
     * 
     * @param j Coordinate index.
     * @param k Parameter index.
     * @return Real 
     */
    Real Line21::operator ()(const Natural &j, const Natural &k) const {
        #ifndef NDEBUG // Integrity check.
        assert(j <= 2);
        assert(k <= 1);
        #endif

        if(j == 0) {
            if(k == 0)
                return this->_a;

            return this->_x0;
        }

        if(j == 1) {
            if(k == 0)
                return this->_b;

            return this->_y0;
        }

        if(k == 0)
            return this->_c;

        return this->_t0;
    }

    /**
     * @brief Scalar insert.
     * 
     * @param j Coordinate index.
     * @param k Parameter index.
     * @param parameter Parameter.
     */
    void Line21::operator ()(const Natural &j, const Natural &k, const Real &parameter) {
        #ifndef NDEBUG // Integrity check.
        assert(j <= 2);
        assert(k <= 1);
        #endif

        if(j == 0) {
            if(k == 0)
                this->_a = parameter;

            this->_x0 = parameter;
        }

        if(j == 1) {
            if(k == 0)
                this->_b = parameter;

            this->_y0 = parameter;
        }

        if(k == 0)
            this->_c = parameter;

        this->_t0 = parameter;
    }

    // Call operator.

    /**
     * @brief Line evaluation.
     * 
     * @param s Parameter.
     * @return Point21 
     */
    Point21 Line21::operator ()(const Real &s) const { return Point21{this->_a * s + this->_x0, this->_b * s + this->_y0, this->_c * s + this->_t0}; }

    // Output.

    /**
     * @brief Line output.
     * 
     * @param ost 
     * @param line Line.
     * @return std::ostream& 
     */
    std::ostream &operator <<(std::ostream &ost, const Line21 &line) {
        ost << "x: " << line._a << "s + " << line._x0 << std::endl;
        ost << "y: " << line._b << "s + " << line._y0 << std::endl;
        return ost << "t: " << line._c << "s + " << line._t0 << std::flush;
    }

}