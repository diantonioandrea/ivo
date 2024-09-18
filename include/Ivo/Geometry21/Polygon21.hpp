/**
 * @file Polygon21.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 2 + 1 polygons.
 * @date 2024-07-22
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef GEOMETRY21_POLYGON21
#define GEOMETRY21_POLYGON21

#include "./Line21.hpp"

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

            /**
             * @brief Polygon's points.
             * 
             * @return std::vector<Point21> 
             */
            inline std::vector<Point21> points() const { return this->_points; }

            std::vector<Edge21> edges() const;

            // Constructors and copy operators.

            Polygon21(const std::vector<Point21> &);
            Polygon21(const std::initializer_list<Point21> &);
            Polygon21(const Polygon21 &);
            Polygon21 &operator =(const std::initializer_list<Point21> &);
            Polygon21 &operator =(const Polygon21 &);

            // Access.

            Point21 operator ()(const Natural &) const;
            Point21 &operator [](const Natural &);

            // Insert.

            void operator ()(const Natural &, const Point21 &);

            // Output.
            
            friend std::ostream &operator <<(std::ostream &, const Polygon21 &);
    };

}

#endif