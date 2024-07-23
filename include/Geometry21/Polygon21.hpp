/**
 * @file Polygon21.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 2 + 1 polygons.
 * @date 2024-07-22
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef GEOMETRY_POLYGON21
#define GEOMETRY_POLYGON21

#include "./Includes.hpp"
#include "./Edge21.hpp"
#include "./Point21.hpp"

namespace ivo {

    /**
     * @brief Polygon21. {p0, p1, ..., pN}.
     * 
     */
    class Polygon21 {

        private:

            // Attributes.

            /**
             * @brief Polygon's points.
             * Counterclockwise with respect to the polygon's plane (if any).
             * 
             */
            std::vector<Point21> _points;

        public:

            // Attributes access.

            std::vector<Point21> points() const;
            std::vector<Edge21> edges() const;

            // Constructors and copy operators.

            Polygon21(const std::vector<Point21> &);
            Polygon21(const Polygon21 &);
            Polygon21 &operator =(const Polygon21 &);

            // Subscript operator, legacy scalar access (C++23).

            #if __cplusplus > 202002L
            Point21 &operator [](const Natural &);
            #endif

            // Call operator, subscript behaviour.

            Point21 operator ()(const Natural &) const;
            void operator ()(const Natural &, const Point21 &);

            // Output.
            
            friend std::ostream &operator <<(std::ostream &, const Polygon21 &);
    };

}

#endif