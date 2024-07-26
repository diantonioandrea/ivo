/**
 * @file Geometry21_Methods_Polygon21.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Geometry/Methods/Polygon21.hpp implementation.
 * @date 2024-07-23
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Ivo.hpp>

namespace ivo {

    // Polygon methods.

    /**
     * @brief Area of a polygon.
     * 
     * @param polygon Polygon.
     * @return Real 
     */
    Real area(const Polygon21 &polygon) {
        std::vector<Point21> points = polygon.points();
        points.emplace_back(points[0]);
        Vector<Real> _area{3};

        for(Natural j = 1; j < points.size(); ++j)
            _area += cross(Vector<Real>{points[j]}, Vector<Real>{points[j - 1]});

        return norm(0.5 * _area);
    }

    /**
     * @brief Centre of a polygon.
     * 
     * @param polygon 
     * @return Point21 
     */
    Point21 centre(const Polygon21 &polygon) {
        std::vector<Point21> points = polygon.points();
        Point21 _centre;

        for(Natural j = 0; j < points.size(); ++j)
            _centre += points[j];

        return _centre / static_cast<Real>(points.size());
    }

    /**
     * @brief Centroid of a polygon.
     * 
     * @param polygon Polygon.
     * @return Point21 
     */
    Point21 centroid(const Polygon21 &polygon) {
        if(polygon.points().size() == 3)
            return centre(polygon);

        // Average on triangles.
        std::vector<Polygon21> triangles;
        Point21 _centre{centre(polygon)};
        Point21 _centroid;

        for(const auto &edge: polygon.edges())
            triangles.emplace_back(Polygon21{{edge(0), edge(1), _centre}});

        for(const auto &triangle: triangles)
            _centroid += area(triangle) * centre(triangle);

        return _centroid / area(polygon);
    }

    // Checks.

    /**
     * @brief Spatial(polygon).
     * 
     * @param polygon Polygon.
     * @return true 
     * @return false 
     */
    bool spatial(const Polygon21 &polygon) {
        std::vector<Point21> points = polygon.points();

        for(Natural j = 1; j < points.size(); ++j)
            if(std::abs(points[j](2) - points[0](2)) > GEOMETRICAL_ZERO)
                return false;

        return true;
    }

}