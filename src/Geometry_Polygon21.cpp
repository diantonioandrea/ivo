/**
 * @file Geometry_Polygon21.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Geometry/Polygon21.hpp implementation.
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
     * @param points 
     */
    Polygon21::Polygon21(const std::vector<Point21> &points) {
        #ifndef NDEBUG // Integrity check.
        assert(points.size() > 2);
        for(Natural j = 0; j < points.size(); ++j)
            for(Natural k = j + 1; k < points.size(); ++k)
                assert(points[j] != points[k]);
        #endif

        this->_points.resize(points.size());
        std::copy(points.begin(), points.end(), this->_points.begin());
    }

    /**
     * @brief Copy constructor.
     * 
     * @param polygon Polygon.
     */
    Polygon21::Polygon21(const Polygon21 &polygon) {
        this->_points.resize(polygon._points.size());
        std::copy(polygon._points.begin(), polygon._points.end(), this->_points.begin());
    }

    /**
     * @brief Copy operator.
     * 
     * @param polygon Polygon.
     * @return Polygon21& 
     */
    Polygon21 &Polygon21::operator =(const Polygon21 &polygon) {
        this->_points.resize(polygon._points.size());
        std::copy(polygon._points.begin(), polygon._points.end(), this->_points.begin());
        return *this;
    }

    // Subscript operator, legacy scalar access (C++23).

    #if __cplusplus > 202002L

    /**
     * @brief Scalar reference access, legacy.
     * 
     * @param j Index.
     * @return Point21& 
     */
    Point21 &Polygon21::operator [](const Natural &j) {
        #ifndef NDEBUG // Integrity check.
        assert(j < this->_points.size());
        #endif

        return this->_points[j];
    }

    #endif

    // Call operator, subscript behaviour.

    /**
     * @brief Scalar access.
     * 
     * @param j Index.
     * @return Point21 
     */
    Point21 Polygon21::operator ()(const Natural &j) const {
        #ifndef NDEBUG // Integrity check.
        assert(j < this->_points.size());
        #endif

        return this->_points[j];
    }

    /**
     * @brief Scalar insert.
     * 
     * @param j Index.
     * @param point Point.
     */
    void Polygon21::operator ()(const Natural &j, const Point21 &point) {
        #ifndef NDEBUG // Integrity check.
        assert(j < this->_points.size());
        #endif

        this->_points[j] = point;
    }

    // Output.
    
    /**
     * @brief Polygon output.
     * 
     * @param ost 
     * @param polygon Polygon.
     * @return std::ostream& 
     */
    std::ostream &operator <<(std::ostream &ost, const Polygon21 &polygon) {
        ost << "{";

        for(Natural j = 0; j < polygon._points.size() - 1; ++j)
            ost << polygon._points[j] << ", ";
        
        return ost << polygon._points[polygon._points.size() - 1] << "}" << std::flush;
    }

}