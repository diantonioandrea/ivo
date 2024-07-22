/**
 * @file Point.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 2 + 1 Points.
 * @date 2024-07-22
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef GEOMETRY_POINT21
#define GEOMETRY_POINT21

#include "./Includes.hpp"

namespace ivo {

    /**
     * @brief Point21. (x, y, t).
     * 
     */
    class Point21 {

        private:

            // Attributes.

            /**
             * @brief Point's x.
             * 
             */
            Real _x = 0.0L;

            /**
             * @brief Point's y.
             * 
             */
            Real _y = 0.0L;

            /**
             * @brief Point's t.
             * 
             */
            Real _t = 0.0L;

        public:

            // Constructors and copy operators.

            Point21(const Real &, const Real &);
            Point21(const Real &, const Real &, const Real &);
            Point21(const Point21 &);
            Point21 &operator =(const Point21 &);

            // Comparison.

            bool operator ==(const Point21 &) const;
            bool operator !=(const Point21 &) const;

            // Subscript operator, legacy scalar access (C++23).

            #if __cplusplus > 202002L
            Real &operator [](const Natural &);
            #endif

            // Call operator, subscript behaviour.

            Real operator ()(const Natural &) const;
            void operator ()(const Natural &, const Real &);

            // Operations.

            Point21 operator +() const;
            Point21 operator -() const;

            Point21 operator +(const Real &) const;
            friend Point21 operator +(const Real &, const Point21 &);
            Point21 &operator +=(const Real &);

            Point21 operator -(const Real &) const;
            friend Point21 operator -(const Real &, const Point21 &);
            Point21 &operator -=(const Real &);

            Point21 operator *(const Real &) const;
            friend Point21 operator *(const Real &, const Point21 &);
            Point21 &operator *=(const Real &);

            Point21 operator /(const Real &) const;
            friend Point21 operator /(const Real &, const Point21 &);
            Point21 &operator /=(const Real &);

            Point21 operator +(const Point21 &) const;
            Point21 &operator +=(const Point21 &);

            Point21 operator -(const Point21 &) const;
            Point21 &operator -=(const Point21 &);

            // Output.
            
            friend std::ostream &operator <<(std::ostream &, const Point21 &);
    };

}

#endif