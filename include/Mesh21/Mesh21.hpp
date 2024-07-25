/**
 * @file Mesh21.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 2 + 1 meshes.
 * @date 2024-07-25
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef MESH21_MESH21
#define MESH21_MESH21

#include "./Includes.hpp"
#include "./Element21.hpp"

namespace ivo {

    /**
     * @brief Mesh21. 2 + 1 mesh.
     * 
     */
    class Mesh21 {

        private:

            // Attributes.

            /**
             * @brief Mesh' space cells.
             * 
             */
            Natural _cells;

            /**
             * @brief Mesh' time levels.
             * 
             */
            Natural _levels;

            /**
             * @brief Mesh' elements.
             * 
             */
            std::vector<Element21> _elements;

            /**
             * @brief Mesh' time intervals.
             * 
             */
            std::vector<Real> _times;

        public:

            // [!]

            // Call operator, subscript behaviour.

            Element21 operator ()(const Natural &, const Natural &) const;

    };

}

#endif