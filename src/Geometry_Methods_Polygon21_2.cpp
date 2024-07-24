/**
 * @file Geometry_Methods_Polygon21_2.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Geometry/Methods/Polygon21_2.hpp implementation.
 * @date 2024-07-23
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Ivo.hpp>

namespace ivo {

    // Polygon methods.

    /**
     * @brief Polygon's bounding box.
     * 
     * @param polygon Polygon.
     * @return std::array<Point21, 2> 
     */
    std::array<Point21, 2> box(const Polygon21 &polygon) {
        #ifndef NDEBUG // Integrity check.
        assert(spatial(polygon));
        #endif

        // Points.
        std::vector<Point21> points = polygon.points();

        // Coordinates.
        Real min_x = points[0](0), min_y = points[0](1);
        Real max_x = points[0](0), max_y = points[0](1);
        Real t = points[0](2);

        for(const auto &point: polygon.points()) {
            min_x = (point(0) < min_x) ? point(0) : min_x;
            min_y = (point(1) < min_y) ? point(1) : min_y;
            max_x = (point(0) > max_x) ? point(0) : max_x;
            max_y = (point(1) > max_y) ? point(1) : max_y;
        }

        return {Point21{min_x, min_y, t}, Point21{max_x, max_y, t}};
    }

    // Containment.

    /**
     * @brief Contains(polygon, point).
     * Space only.
     * 
     * @param polygon Polygon.
     * @param point Point.
     * @return true 
     * @return false 
     */
    bool contains2(const Polygon21 &polygon, const Point21 &point) {
        #ifndef NDEBUG // Integrity check.
        assert(spatial(polygon));
        #endif

        // First check.
        std::vector<Point21> points = polygon.points();
        if(std::abs(point(2) - points[0](2)) > GEOMETRICAL_ZERO)
            return false;

        // Second check.
        Line21 line{Point21{point(0), 0.0L, 0.0L}, Point21{point(0) + 1.0L, 0.0L, 0.0L}};
        Natural counter = 0;

        for(const auto &intersection: intersections(line, polygon))
            if(intersection(0) >= point(0))
                ++counter;

        return static_cast<bool>(counter % 2);
    }

}