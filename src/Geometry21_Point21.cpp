/**
 * @file Geometry21_Point21.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Geometry/Point21.hpp implementation.
 * @date 2024-07-22
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Ivo.hpp>

namespace ivo {

    /**
     * @brief Zero constructor.
     * 
     */
    Point21::Point21(): _x{0.0L}, _y{0.0L}, _t{0.0L} {}

    /**
     * @brief Space constructor.
     * 
     * @param x Point's x.
     * @param y Point's y.
     */
    Point21::Point21(const Real &x, const Real &y): _x{x}, _y{y}, _t{0.0L} {}

    /**
     * @brief Space-time constructor.
     * 
     * @param x Point's x.
     * @param y Point's y.
     * @param t Point's t.
     */
    Point21::Point21(const Real &x, const Real &y, const Real &t): _x{x}, _y{y}, _t{t} {}

    /**
     * @brief Copy constructor.
     * 
     * @param point Point.
     */
    Point21::Point21(const Point21 &point): _x{point._x}, _y{point._y}, _t{point._t} {}

    /**
     * @brief Copy operator.
     * 
     * @param point Point.
     * @return Point21& 
     */
    Point21 &Point21::operator =(const Point21 &point) {
        this->_x = point._x;
        this->_y = point._y;
        this->_t = point._t;
        return *this;
    }

    // Comparison.

    /**
     * @brief Point == point.
     * 
     * @param point Point.
     * @return true 
     * @return false 
     */
    bool Point21::operator ==(const Point21 &point) const {
        return distance(*this, point) <= GEOMETRICAL_ZERO;
    }

    /**
     * @brief Point != point.
     * 
     * @param point Point.
     * @return true 
     * @return false 
     */
    bool Point21::operator !=(const Point21 &point) const {
        return distance(*this, point) > GEOMETRICAL_ZERO;
    }

    // Subscript operator, legacy scalar access (C++23).

    #if __cplusplus > 202002L

    /**
     * @brief Scalar reference access, legacy.
     * 
     * @param j Index.
     * @return Real& 
     */
    Real &Point21::operator [](const Natural &j) {
        #ifndef NDEBUG // Integrity check.
        assert(j <= 2);
        #endif

        if(j == 0)
            return this->_x;

        if(j == 1)
            return this->_y;
        
        return this->_t;
    }

    #endif

    // Call operator, subscript behaviour.

    /**
     * @brief Scalar access.
     * 
     * @param j Index.
     * @return Real 
     */
    Real Point21::operator ()(const Natural &j) const {
        #ifndef NDEBUG // Integrity check.
        assert(j <= 2);
        #endif

        if(j == 0)
            return this->_x;

        if(j == 1)
            return this->_y;
        
        return this->_t;
    }

    /**
     * @brief Scalar insert.
     * 
     * @param j Index.
     * @param scalar Scalar.
     */
    void Point21::operator ()(const Natural &j, const Real &scalar) {
        #ifndef NDEBUG // Integrity check.
        assert(j <= 2);
        #endif

        if(j == 0)
            this->_x = scalar;

        if(j == 1)
            this->_y = scalar;
        
        this->_t = scalar;
    }
    
    // Operations.

    /**
     * @brief +Point.
     * 
     * @return Point21 
     */
    Point21 Point21::operator +() const { return *this; }

    /**
     * @brief -Point.
     * 
     * @return Point21 
     */
    Point21 Point21::operator -() const { return Point21{-this->_x, -this->_y, -this->_t}; }

    /**
     * @brief Point + scalar.
     * 
     * @param scalar Scalar.
     * @return Point21 
     */
    Point21 Point21::operator +(const Real &scalar) const { return Point21{this->_x + scalar, this->_y + scalar, this->_t + scalar}; }

    /**
     * @brief Scalar + point.
     * 
     * @param scalar Scalar.
     * @param point Point.
     * @return Point21 
     */
    Point21 operator +(const Real &scalar, const Point21 &point) { return Point21{scalar + point._x, scalar + point._y, scalar + point._t}; }

    /**
     * @brief Point += scalar.
     * 
     * @param scalar Scalar.
     * @return Point21& 
     */
    Point21 &Point21::operator +=(const Real &scalar) {
        this->_x += scalar;
        this->_y += scalar;
        this->_t += scalar;
        return *this;
    }

    /**
     * @brief Point - scalar.
     * 
     * @param scalar Scalar.
     * @return Point21 
     */
    Point21 Point21::operator -(const Real &scalar) const { return Point21{this->_x - scalar, this->_y - scalar, this->_t - scalar}; }

    /**
     * @brief Scalar - point.
     * 
     * @param scalar Scalar.
     * @param point Point.
     * @return Point21 
     */
    Point21 operator -(const Real &scalar, const Point21 &point) { return Point21{scalar - point._x, scalar - point._y, scalar - point._t}; }

    /**
     * @brief Point -= scalar.
     * 
     * @param scalar Scalar.
     * @return Point21& 
     */
    Point21 &Point21::operator -=(const Real &scalar) {
        this->_x -= scalar;
        this->_y -= scalar;
        this->_t -= scalar;
        return *this;
    }

    /**
     * @brief Point * scalar.
     * 
     * @param scalar Scalar.
     * @return Point21 
     */
    Point21 Point21::operator *(const Real &scalar) const { return Point21{this->_x * scalar, this->_y * scalar, this->_t * scalar}; }

    /**
     * @brief Scalar * point.
     * 
     * @param scalar Scalar.
     * @param point Point.
     * @return Point21 
     */
    Point21 operator *(const Real &scalar, const Point21 &point) { return Point21{scalar * point._x, scalar * point._y, scalar * point._t}; }

    /**
     * @brief Point *= scalar.
     * 
     * @param scalar Scalar.
     * @return Point21& 
     */
    Point21 &Point21::operator *=(const Real &scalar) {
        this->_x *= scalar;
        this->_y *= scalar;
        this->_t *= scalar;
        return *this;
    }

    /**
     * @brief Point / scalar.
     * 
     * @param scalar Scalar.
     * @return Point21 
     */
    Point21 Point21::operator /(const Real &scalar) const { return Point21{this->_x / scalar, this->_y / scalar, this->_t / scalar}; }

    /**
     * @brief Scalar / point.
     * 
     * @param scalar Scalar.
     * @param point Point.
     * @return Point21 
     */
    Point21 operator /(const Real &scalar, const Point21 &point) { return Point21{scalar / point._x, scalar / point._y, scalar / point._t}; }

    /**
     * @brief Point /= scalar.
     * 
     * @param scalar Scalar.
     * @return Point21& 
     */
    Point21 &Point21::operator /=(const Real &scalar) {
        this->_x /= scalar;
        this->_y /= scalar;
        this->_t /= scalar;
        return *this;
    }

    /**
     * @brief Point + point.
     * 
     * @param point Point.
     * @return Point21 
     */
    Point21 Point21::operator +(const Point21 &point) const { return Point21{this->_x + point._x, this->_y + point._y, this->_t + point._t}; }

    /**
     * @brief Point += point.
     * 
     * @param point Point.
     * @return Point21& 
     */
    Point21 &Point21::operator +=(const Point21 &point) {
        this->_x += point._x;
        this->_y += point._y;
        this->_t += point._t;
        return *this;
    }

    /**
     * @brief Point - point.
     * 
     * @param point Point.
     * @return Point21 
     */
    Point21 Point21::operator -(const Point21 &point) const { return Point21{this->_x - point._x, this->_y - point._y, this->_t - point._t}; }

    /**
     * @brief Point -= point.
     * 
     * @param point Point.
     * @return Point21& 
     */
    Point21 &Point21::operator -=(const Point21 &point) {
        this->_x -= point._x;
        this->_y -= point._y;
        this->_t -= point._t;
        return *this;
    }

    // Output.

    /**
     * @brief Point output.
     * 
     * @param ost 
     * @param point Point.
     * @return std::ostream& 
     */
    std::ostream &operator <<(std::ostream &ost, const Point21 &point) {
        return ost << "(" << point._x << ", " << point._y << "; " << point._t << ")" << std::flush;
    }

    // Literals.

    /**
     * @brief x-only point.
     * 
     * @param x Coordinate.
     * @return Point21 
     */
    Point21 operator ""_x(const Real x) { return Point21{x, 0.0L, 0.0L}; }

    /**
     * @brief y-only point.
     * 
     * @param y Coordinate.
     * @return Point21 
     */
    Point21 operator ""_y(const Real y) { return Point21{0.0L, y, 0.0L}; }

    /**
     * @brief t-only point.
     * 
     * @param t Coordinate.
     * @return Point21 
     */
    Point21 operator ""_t(const Real t) { return Point21{0.0L, 0.0L, t}; }

}