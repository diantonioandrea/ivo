/**
 * @file Line21.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 2 + 1 lines.
 * @date 2024-07-22
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef GEOMETRY21_LINE21
#define GEOMETRY21_LINE21

#include "./Edge21.hpp"

namespace ivo {

    /**
     * @brief Line21. (x, y, t) = (a, b, c) u + (x0, y0, z0).
     * 
     */
    class Line21 {

        private:

            /**
             * @brief Line's a.
             * 
             */
            Real _a;

            /**
             * @brief Line's b.
             * 
             */
            Real _b;

            /**
             * @brief Line's c.
             * 
             */
            Real _c;

            /**
             * @brief Line's x0.
             * 
             */
            Real _x0;

            /**
             * @brief Line's y0.
             * 
             */
            Real _y0;

            /**
             * @brief Line's t0.
             * 
             */
            Real _t0;

        public:

            // Constructors and copy operators.

            Line21(const Point21 &, const Point21 &);
            Line21(const Edge21 &);
            Line21(const Line21 &);
            Line21 &operator =(const Line21 &);

            // Access.

            Real operator ()(const Natural &, const Natural &) const;
            Real &operator [](const Natural &, const Natural &);

            // Insert.

            void operator ()(const Natural &, const Natural &, const Real &);

            // Call operator.

            Point21 operator ()(const Real &) const;

            // Output.

            friend std::ostream &operator <<(std::ostream &, const Line21 &);

    };

}

#endif