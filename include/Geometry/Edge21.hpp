/**
 * @file Edge21.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 2 + 1 Edges.
 * @date 2024-07-22
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef GEOMETRY_EDGE21
#define GEOMETRY_EDGE21

#include "./Includes.hpp"
#include "./Point21.hpp"

namespace ivo {

    /**
     * @brief Edge21. [a, b].
     * 
     */
    class Edge21 {

        private:

            // Attributes.

            /**
             * @brief Edge's a.
             * 
             */
            Point21 _a;

            /**
             * @brief Edge's b.
             * 
             */
            Point21 _b;

        public:

            // Constructors and copy operators.

            Edge21(const Point21 &, const Point21 &);
            Edge21(const Edge21 &);
            Edge21 &operator =(const Edge21 &);

            // Comparison.

            bool operator ==(const Edge21 &) const;
            bool operator !=(const Edge21 &) const;

            // Subscript operator, legacy scalar access (C++23).

            #if __cplusplus > 202002L
            Point21 &operator [](const Natural &);
            #endif

            // Call operator, subscript behaviour.

            Point21 operator ()(const Natural &) const;
            void operator ()(const Natural &, const Point21 &);

            // Output.
            
            friend std::ostream &operator <<(std::ostream &, const Edge21 &);
    };

}

#endif