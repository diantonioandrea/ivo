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

    /**
     * @brief Generates random points in a polygon
     * 
     * @param polygon Polygon.
     * @param number Number.
     * @return std::vector<Point21> 
     */
    std::vector<Point21> random2(const Polygon21 &polygon, const Natural &number) {

        // Seeding.
        std::srand(std::time(nullptr));

        #ifndef NDEBUG // integrity check.
        assert(number > 0);
        assert(spatial(polygon));
        #endif

        // Random points.
        std::vector<Point21> points(number, Point21{});

        // Boundaries.
        auto [min_xy, max_xy] = box2(polygon);
        Real t = min_xy(2);

        // Generation.
        for(Natural j = 0; j < number; ++j) {
            bool generation = true;

            while(generation) {
                Real x = min_xy(0) + (max_xy(0) - min_xy(0)) * (static_cast<Real>(std::rand()) / static_cast<Real>(RAND_MAX));
                Real y = min_xy(1) + (max_xy(1) - min_xy(1)) * (static_cast<Real>(std::rand()) / static_cast<Real>(RAND_MAX));
                Point21 point{x, y, t};

                if(!contains2(polygon, point))
                    continue;
                
                generation = false;
                for(Natural k = 0; k < j; ++k)
                    if(point == points[k]) {
                        generation = true;
                        break;
                    }

                if(!generation)
                    points[j] = point;
            }
        }

        return points;
    }

    /**
     * @brief Returns the Voronoi diagram of a given set of points inside a polygon.
     * 
     * @param polygon Polygon.
     * @param points Points.
     * @return std::vector<Polygon21> 
     */
    std::vector<Polygon21> voronoi2(const Polygon21 &polygon, const std::vector<Point21> &points) {
        #ifndef NDEBUG // Integrity check.
        assert(spatial(polygon));
        #endif

        // Cells.
        std::vector<Polygon21> cells;

        for(Natural j = 0; j < points.size(); ++j) {
            Polygon21 cell{polygon};

            #ifndef NDEBUG // Integrity check.
            assert(contains2(polygon, points[j]));
            #endif

            // Reductions.
            for(Natural k = 0; k < points.size(); ++k) {
                if(j == k)
                    continue;
                    
                cell = reduce2(cell, bisector2(points[j], points[k]), points[j]);
            }

            cells.emplace_back(cell);
        }

        return cells;
    }

    /**
     * @brief Returns the Voronoi diagram of a random set of points inside a polygon.
     * 
     * @param polygon Polygon.
     * @param number Number.
     * @return std::vector<Polygon21> 
     */
    std::vector<Polygon21> voronoi2(const Polygon21 &polygon, const Natural &number) {
        return voronoi2(polygon, random2(polygon, number));
    }

    // (Meshing) diagram postprocessing.

    /**
     * @brief Lloyd's relaxation process.
     * 
     * @param polygon Polygon.
     * @param diagram Polygons.
     */
    void lloyd2(const Polygon21 &polygon, std::vector<Polygon21> &diagram) {

        // Lloyd's steps.
        Natural steps = 16 + diagram.size();

        // Centroids.
        std::vector<Point21> centroids;

        for(const auto &cell: diagram)
            centroids.emplace_back(centroid(cell));

        // Diagram update, first step.
        diagram = voronoi2(polygon, centroids);

        for(Natural j = 1; j < steps; ++j) {

            // Centroids and residual update.
            Real residual = 0.0L;
            for(Natural k = 0; k < diagram.size(); ++k) {
                residual += distance(centroids[k], centroid(diagram[k]));
                centroids[k] = centroid(diagram[k]);
            }

            if(residual <= LLOYD_STOP * diagram.size())
                return;

            // Diagram update.
            diagram = voronoi2(polygon, centroids);
        }
    }

    /**
     * @brief Diagram collapser.
     * 
     * @param polygon Polygon.
     * @param diagram Polygons.
     */
    void collapse2(const Polygon21 &polygon, std::vector<Polygon21> &diagram) {
        for(Natural j = 0; j < diagram.size(); ++j) {

            // j-th cell.
            Polygon21 cell_j = diagram[j];
            std::vector<Edge21> edges_j = cell_j.edges();

            // j-th size.
            Real size_j = 0.0L;

            for(const auto &edge_j: edges_j)
                if(edge_j.size() > size_j)
                    size_j = edge_j.size();

            for(Natural ej = 0; ej < edges_j.size(); ++ej) {
                if(edges_j[ej].size() > COLLAPSE * size_j)
                    continue;

                // Edge to be collapsed and midpoint.
                Edge21 edge_j = edges_j[ej];
                Point21 midpoint = (edge_j(0) + edge_j(1)) / 2.0L;

                for(Natural k = 0; k < diagram.size(); ++k) {
                    if(j == k)
                        continue;
                    
                    // Collapse flag.
                    bool collapse = false;
                    
                    // k-th cell.
                    Polygon21 cell_k = diagram[k];
                    std::vector<Edge21> edges_k = cell_k.edges();
                    std::vector<Point21> points_k = cell_k.points();

                    // Looks for common edges.
                    for(Natural ek = 0; ek < edges_k.size(); ++ek) {
                        if(edges_k[ek] == edge_j) {
                            if(ek < edges_k.size() - 1) {
                                points_k.erase(points_k.begin() + ek + 1);
                                points_k[ek] = midpoint;
                            } else {
                                points_k.pop_back();
                                points_k[0] = midpoint;
                            }

                            diagram[k] = Polygon21(points_k);
                            collapse = true;
                            break;
                        }
                    }

                    if(collapse)
                        continue;

                    // Looks for common points.
                    for(Natural pk = 0; pk < points_k.size(); ++pk)
                        if((points_k[pk] == edge_j(0)) || (points_k[pk] == edge_j(1))) {
                            diagram[k](pk, midpoint);
                            break;
                        }
                }
            }
        }
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
        Line21 line{point, point + 1.0_x};
        Natural counter = 0;

        for(const auto &intersection: intersections(line, polygon))
            if(intersection(0) >= point(0))
                ++counter;

        return static_cast<bool>(counter % 2);
    }

}