/**
 * @file Geometry21_Methods_Polygon21_2.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Geometry/Methods/Polygon21_2.hpp implementation.
 * @date 2024-07-23
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Ivo.hpp>

namespace ivo {

    // (Meshing) polygon methods.

    /**
     * @brief Polygon's reduction with respect to a line and a point.
     * Space only.
     * 
     * @param polygon 
     * @param point 
     * @param line 
     * @return Polygon21 
     */
    Polygon21 reduce2(const Polygon21 &polygon, const Line21 &line, const Point21 &point) {
        std::vector<Point21> p_points = polygon.points(); // Polygon's points.
        std::vector<Point21> i_points = intersections(line, polygon); // Intersections.

        #ifndef NDEBUG // Integrity check.
        assert(spatial(polygon));
        assert(spatial(line));
        assert(contains2(polygon, point));
        assert(std::abs(line(2, 1) - p_points[0](2)) <= GEOMETRICAL_ZERO);
        #endif
        
        // Not enough intersections.
        if(i_points.size() <= 1)
            return polygon;
        
        std::array<Point21, 2> points = {i_points[0], i_points[1]}; // New points.

        // More than two intersections (Shouldn't happen).
        if(i_points.size() > 2)
            for(Natural j = 0; j < i_points.size(); ++j)
                for(Natural k = 0; k < i_points.size(); ++k)
                    if(distance(i_points[j], i_points[k]) > distance(i_points[0], i_points[1]))
                        points = {i_points[j], i_points[k]};

        // Indices.
        std::vector<Edge21> p_edges = polygon.edges(); // Polygon's edges.
        std::array<Natural, 2> indices = {0, 1}; // New points' indices.

        for(Natural j = 0; j < p_edges.size(); ++j) {
            Edge21 p_edge = p_edges[j];

            if(contains(p_edge, points[0]) && contains(p_edge, points[1])) // Shouldn't happen.
                return polygon;

            if(contains(p_edge, points[0]) && (p_edge(1) != points[0]))
                indices[0] = j;
            
            if(contains(p_edge, points[1]) && (p_edge(1) != points[1]))
                indices[1] = j;
        }

        if(indices[0] > indices[1]) {
            std::swap(indices[0], indices[1]);
            std::swap(points[0], points[1]);
        }
        
        // Creates two new polygons A, B.
        std::vector<Point21> a_points, b_points;

        for(Natural j = 0; j < p_points.size(); ++j) {
            if((j <= indices[0]) || (j > indices[1]))
                a_points.emplace_back(p_points[j]);
            
            if(j == indices[0]) {
                if(points[0] != *--a_points.end()) // Avoids duplicates.
                    a_points.emplace_back(points[0]);

                b_points.emplace_back(points[0]);
            }

            if((j > indices[0]) && (j <= indices[1]))
                if(p_points[j] != *--b_points.end()) // Avoids duplicates.
                    b_points.emplace_back(p_points[j]);

            if(j == indices[1]) {
                a_points.emplace_back(points[1]);
                b_points.emplace_back(points[1]);
            }
        }

        // Polygon's choice.
        Polygon21 A{a_points}, B{b_points};

        if(contains2(A, point))
            return A;

        #ifndef NDEBUG // Integrity check.
        assert(contains2(B, point));
        #endif

        return B;
    }

    // Polygon methods.

    /**
     * @brief Polygon's bounding box.
     * Space only.
     * 
     * @param polygon Polygon.
     * @return std::array<Point21, 2> 
     */
    std::array<Point21, 2> box2(const Polygon21 &polygon) {
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