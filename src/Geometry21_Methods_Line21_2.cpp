/**
 * @file Geometry21_Methods_Line21_2.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Geometry/Methods/Line21_2.hpp implementation.
 * @date 2024-07-23
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Ivo.hpp>

namespace ivo {

    // Line methods.

    /**
     * @brief Edge's bisector.
     * Space only.
     * 
     * @param edge Edge.
     * @return Line21 
     */
    Line21 bisector2(const Edge21 &edge) {
        #ifndef NDEBUG // Integrity check.
        assert(spatial(edge));
        #endif

        // Midpoint (p0) and difference.
        Point21 p0 = (edge(0) + edge(1)) / 2.0L;
        Point21 difference = edge(1) - edge(0);

        // Edge vector.
        Vector<Real> e_vector{{difference(0), difference(1)}};

        // Normal vector.
        Vector<Real> n_vector{{-e_vector(1), e_vector(0)}};
        n_vector /= norm(n_vector);

        // Shifted midpoint.
        Point21 p1{p0(0) + n_vector(0), p0(1) + n_vector(1), p0(2)};

        return Line21{p0, p1};
    }

    /**
     * @brief Two points' bisector.
     * 
     * @param p Point.
     * @param q Point.
     * @return Line21 
     */
    Line21 bisector2(const Point21 &p, const Point21 &q) { return bisector2(Edge21{p, q}); }

}