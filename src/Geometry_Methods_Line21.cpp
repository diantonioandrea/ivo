/**
 * @file Geometry_Methods_Line21.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Geometry/Methods/Line21.hpp implementation.
 * @date 2024-07-22
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Ivo.hpp>

namespace ivo {

    // Line methods.

    /**
     * @brief Normal(edge, point).
     * 
     * @param edge Edge.
     * @param point Point.
     * @return std::optional<Line21> 
     */
    std::optional<Line21> normal(const Edge21 &edge, const Point21 &point) {
        Line21 line = normal(Line21{edge}, point);

        if(intersections(line, edge).has_value())
            return line;

        return {};
    }

    /**
     * @brief Normal(line, point).
     * 
     * @param line Line.
     * @param point Point.
     * @return Line21 
     */
    Line21 normal(const Line21 &line, const Point21 &point) {

        // [!]

    }

    /**
     * @brief Bisector2(p, q).
     * Space bisector.
     * 
     * @param p Point.
     * @param q Point.
     * @return Line21 
     */
    Line21 bisector2(const Point21 &p, const Point21 &q) {
        #ifndef NDEBUG
        assert(p != q);
        assert(std::abs(p(2) - q(2)) < GEOMETRICAL_ZERO);
        #endif

        // [!]
    }

    // Intersections.

    /**
     * @brief Intersections(line, edge).
     * 
     * @param line Line.
     * @param edge Edge.
     * @return std::optional<Point21> 
     */
    std::optional<Point21> intersections(const Line21 &line, const Edge21 &edge) {
        std::optional<Point21> intersection = intersections(line, Line21{edge});

        if(intersection.has_value())
            if(contains(edge, intersection.value()))
                return intersection.value();

        return {};
    }

    std::optional<Point21> intersections(const Line21 &, const Line21 &);

    // Distances.

    /**
     * @brief Distance(line, point).
     * 
     * @param line Line.
     * @param point Point.
     * @return Real 
     */
    Real distance(const Line21 &line, const Point21 &point) {

        // [!]

    }

    /**
     * @brief Distance(line, edge).
     * 
     * @param line Line.
     * @param edge Edge.
     * @return Real 
     */
    Real distance(const Line21 &line, const Edge21 &edge) {
        
        // [!]

    }

    /**
     * @brief distance(line, line).
     * 
     * @param r Line.
     * @param s Line.
     * @return Real 
     */
    Real distance(const Line21 &r, const Line21 &s) {

        // [!]

    }

    // Containment.

    /**
     * @brief Contains(line, point).
     * 
     * @param line Line.
     * @param point Point.
     * @return true 
     * @return false 
     */
    bool contains(const Line21 &line, const Point21 &point) {
        for(Natural j = 0; j < 3; ++j) {
            if(std::abs(line(j, 0)) > GEOMETRICAL_ZERO) {
                Real s = (point(j) - line(j, 1)) / line(j, 0);

                for(Natural k = 0; k < 3; ++k) {
                    if(j == k)
                        continue;

                    if(std::abs(point(k) - s * line(k, 0) - line(k, 1)) > GEOMETRICAL_ZERO)
                        return false;
                }
            }
        }

        return true;
    }

    /**
     * @brief Contains(line, edge).
     * 
     * @param line Line.
     * @param edge Edge.
     * @return true 
     * @return false 
     */
    bool contains(const Line21 &line, const Edge21 &edge) { return contains(line, edge(0)) && contains(line, edge(1)); }

}