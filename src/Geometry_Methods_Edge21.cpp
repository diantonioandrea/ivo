/**
 * @file Geometry_Methods_Edge21.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Geometry/Methods/Edge21.hpp implementation.
 * @date 2024-07-22
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Ivo.hpp>

namespace ivo {

    // Intersections.

    /**
     * @brief Intersections(edge, edge).
     * 
     * @param ab Edge.
     * @param cd Edge.
     * @return std::optional<Point21> 
     */
    std::optional<Point21> intersections(const Edge21 &ab, const Edge21 &cd) {
        std::optional<Point21> intersection = intersections(Line21{ab}, Line21{cd});

        if(intersection.has_value())
            if(contains(ab, intersection.value()) && contains(cd, intersection.value()))
                return intersection.value();

        return {};
    }

    // Distances.

    /**
     * @brief Distance(edge, point).
     * Code-wise redundant with respect to distance(line, point).
     * 
     * @param edge Edge.
     * @param point Point.
     * @return Real 
     */
    Real distance(const Edge21 &edge, const Point21 &point) {

        // Naming.
        Line21 line{edge};

        // Point vector.
        Vector<Real> p{{point(0), point(1), point(2)}};

        // Line parameters and reference point.
        Vector<Real> rv{{line(0, 0), line(1, 0), line(2, 0)}};
        Vector<Real> p0{{line(0, 1), line(1, 1), line(2, 1)}};

        // Line's parameter.
        Real t = dot(rv, p0 - p) / dot(rv, rv);

        if(contains(edge, line(t)))
            return distance(line(t), point);

        // Distances.
        std::vector<Real> distances{distance(edge(0), point), distance(edge(1), point)};
        return *std::min_element(distances.begin(), distances.end());

    }

    // Containment.

    /**
     * @brief Contains(edge, point).
     * Evaluated by distances.
     * 
     * @param ab Edge.
     * @param p Point.
     * @return true 
     * @return false 
     */
    bool contains(const Edge21 &ab, const Point21 &p) {
        if((ab(0) == p) || (ab(1) == p))
            return true;
            
        Real ap = distance(ab(0), p);
        Real bp = distance(ab(1), p);

        return std::abs(distance(ab(0), ab(1)) - (ap + bp)) <= GEOMETRICAL_ZERO;
    }

    /**
     * @brief Contains(edge, edge).
     * 
     * @param ab Edge.
     * @param cd Edge.
     * @return true 
     * @return false 
     */
    bool contains(const Edge21 &ab, const Edge21 &cd) {
        return contains(ab, cd(0)) && contains(ab, cd(1));
    }

}