/**
 * @file Edge21.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 2 + 1 edges.
 * @date 2024-07-22
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef GEOMETRY21_EDGE21
#define GEOMETRY21_EDGE21

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

            // Access.

            Point21 operator ()(const Natural &) const;
            Point21 &operator [](const Natural &);

            // Insert.

            void operator ()(const Natural &, const Point21 &);

            // Methods.

            inline Real size() const { return distance(this->_a, this->_b); }

            // Output.
            
            friend std::ostream &operator <<(std::ostream &, const Edge21 &);
    };

}

#endif