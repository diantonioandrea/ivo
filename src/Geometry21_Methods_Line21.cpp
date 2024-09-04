/**
 * @file Geometry21_Methods_Line21.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Geometry/Methods/Line21.hpp implementation.
 * @date 2024-07-22
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Ivo.hpp>

namespace ivo {

    // Intersections.
    
    /**
     * @brief Intersections(line, polygon).
     * 
     * @param line 
     * @param polygon
     * @return std::set<Point21> 
     */
    std::vector<Point21> intersections(const Line21 &line, const Polygon21 &polygon) {
        std::vector<Point21> points;
        
        for(const auto &edge: polygon.edges()) {
            std::optional<Point21> intersection = intersections(line, edge);

            if(!intersection.has_value())
                continue;

            // Unique intersections only.
            bool flag = true;
            for(const auto &point: points)
                if(point == intersection)
                    flag = false;

            if(flag)
                points.emplace_back(intersection.value());
        }

        return points;
    }

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

    /**
     * @brief Intersections(line, line).
     * 
     * @param r Line.
     * @param s Line.
     * @return std::optional<Point21> 
     */
    std::optional<Point21> intersections(const Line21 &r, const Line21 &s) {
        if(distance(r, s) > constants::geometrical_zero)
            return {};

        Natural rj = 0, sj = 0;
        for(Natural j = 0; j < 3; ++j) {
            if(std::abs(r(j, 0)) > constants::geometrical_zero)
                for(Natural k = 0; k < 3; ++k) {
                    if(j == k)
                        continue;
                    
                    if((std::abs(s(k, 0)) > constants::geometrical_zero) && (std::abs(s(k, 0) * r(j, 0) - s(j, 0) * r(k, 0)) > constants::geometrical_zero)) {
                        rj = j;
                        sj = k;
                        break;
                    }
                }

            if(rj != sj)
                break;
        }

        if(rj == sj) // [?]
            return {};

        // r's parameter.
        Real t; 

        if(std::abs(s(rj, 0)) > constants::geometrical_zero)
            t = (s(sj, 0) * (s(rj, 1) - r(rj, 1)) - s(rj, 0) * (s(sj, 1) - r(sj, 1))) / (s(sj, 0) * r(rj, 0) - s(rj, 0) * r(sj, 0));
        else
            t = (s(rj, 1) - r(rj, 1)) / r(rj, 0);

        return r(t);
    }

    // Distances.

    /**
     * @brief Distance(line, point).
     * 
     * @param line Line.
     * @param point Point.
     * @return Real 
     */
    Real distance(const Line21 &line, const Point21 &point) {

        // Point vector.
        Vector<Real> p{{point(0), point(1), point(2)}};

        // Line parameters and reference point.
        Vector<Real> rv{{line(0, 0), line(1, 0), line(2, 0)}};
        Vector<Real> p0{{line(0, 1), line(1, 1), line(2, 1)}};

        // Line's parameter.
        Real t = dot(rv, p - p0) / dot(rv, rv);

        return distance(line(t), point);
    }

    /**
     * @brief Distance(line, edge).
     * Code-wise redundant with respect to distance(line, line).
     * 
     * @param line Line.
     * @param edge Edge.
     * @return Real 
     */
    Real distance(const Line21 &line, const Edge21 &edge) {

        // Naming.
        Line21 r{line}, s{edge};

        // r's parameters and reference point.
        Vector<Real> rv{{r(0, 0), r(1, 0), r(2, 0)}};
        Vector<Real> p0{{r(0, 1), r(1, 1), r(2, 1)}};

        // s's parameters and reference point.
        Vector<Real> sv{{s(0, 0), s(1, 0), s(2, 0)}};
        Vector<Real> q0{{s(0, 1), s(1, 1), s(2, 1)}};

        // Reference points difference.
        Vector<Real> pq = p0 - q0;

        // Some products.
        Real rv_2 = dot(rv, rv);
        Real sv_2 = dot(sv, sv);
        Real rv_pq = dot(rv, pq);
        Real sv_pq = dot(sv, pq);
        Real rv_sv = dot(rv, sv);

        // Parallelism.
        if(std::abs(rv_2 * sv_2 - rv_sv * rv_sv) <= constants::geometrical_zero)
            return distance(r, edge(0));

        // Parameters.
        Real t = (rv_sv * sv_pq - sv_2 * rv_pq) / (rv_2 * sv_2 - rv_sv * rv_sv);
        Real u = (rv_2 * sv_pq - rv_sv * rv_pq) / (rv_2 * sv_2 - rv_sv * rv_sv);

        if(contains(edge, s(u)))
            return distance(r(t), s(u));

        // Distances.
        std::vector<Real> distances{distance(r, edge(0)), distance(r, edge(1))};
        return *std::min_element(distances.begin(), distances.end());
    }

    /**
     * @brief distance(line, line).
     * 
     * @param r Line.
     * @param s Line.
     * @return Real 
     */
    Real distance(const Line21 &r, const Line21 &s) {

        // r's parameters and reference point.
        Vector<Real> rv{{r(0, 0), r(1, 0), r(2, 0)}};
        Vector<Real> p0{{r(0, 1), r(1, 1), r(2, 1)}};

        // s's parameters and reference point.
        Vector<Real> sv{{s(0, 0), s(1, 0), s(2, 0)}};
        Vector<Real> q0{{s(0, 1), s(1, 1), s(2, 1)}};

        // Reference points difference.
        Vector<Real> pq = p0 - q0;

        // Some products.
        Real rv_2 = dot(rv, rv);
        Real sv_2 = dot(sv, sv);
        Real rv_pq = dot(rv, pq);
        Real sv_pq = dot(sv, pq);
        Real rv_sv = dot(rv, sv);

        // Parallelism.
        if(std::abs(rv_2 * sv_2 - rv_sv * rv_sv) <= constants::geometrical_zero)
            return distance(r, s(0.0L));

        // Parameters.
        Real t = (rv_sv * sv_pq - sv_2 * rv_pq) / (rv_2 * sv_2 - rv_sv * rv_sv);
        Real u = (rv_2 * sv_pq - rv_sv * rv_pq) / (rv_2 * sv_2 - rv_sv * rv_sv);

        return distance(r(t), s(u));
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
            if(std::abs(line(j, 0)) > constants::geometrical_zero) {
                Real s = (point(j) - line(j, 1)) / line(j, 0);

                for(Natural k = 0; k < 3; ++k) {
                    if(j == k)
                        continue;

                    if(std::abs(point(k) - s * line(k, 0) - line(k, 1)) > constants::geometrical_zero)
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

    // Checks.

    /**
     * @brief Spatial(line).
     * 
     * @param line Line.
     * @return true 
     * @return false 
     */
    bool spatial(const Line21 &line) { return std::abs(line(2, 0)) <= constants::geometrical_zero; }

}