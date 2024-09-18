/**
 * @file Point.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 2 + 1 points.
 * @date 2024-07-22
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef GEOMETRY21_POINT21
#define GEOMETRY21_POINT21

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
            Real _x;

            /**
             * @brief Point's y.
             * 
             */
            Real _y;

            /**
             * @brief Point's t.
             * 
             */
            Real _t;

        public:

            // Constructors and copy operators.

            Point21();
            Point21(const Real &, const Real &);
            Point21(const Real &, const Real &, const Real &);
            Point21(const Point21 &);
            Point21 &operator =(const Point21 &);

            // Conversions.

            operator Vector<Real>() const;

            // Comparison.

            bool operator ==(const Point21 &) const;
            bool operator !=(const Point21 &) const;

            // Access.

            Real operator ()(const Natural &) const;
            Real &operator [](const Natural &);

            // Insert.

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

    // Literals.

    Point21 operator ""_x(const Real);
    Point21 operator ""_y(const Real);
    Point21 operator ""_t(const Real);

}

#endif