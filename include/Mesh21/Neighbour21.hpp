/**
 * @file Neighbour21.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 2 + 1 neighbours.
 * @date 2024-07-25
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef MESH21_NEIGHBOUR21
#define MESH21_NEIGHBOUR21

#include "./Includes.hpp"

namespace ivo {

    /**
     * @brief Neighbour21. Neighbouring structure.
     * 
     */
    class Neighbour21 {

        private:

            // Attributes.

            /**
             * @brief Top element.
             * 
             */
            Integer _top;

            /**
             * @brief Bottom element.
             * 
             */
            Integer _bottom;

            /**
             * @brief Facing elements.
             * 
             */
            std::vector<std::array<Integer, 2>> _facing;

        public:

            // Attributes access.

            /**
             * @brief Top element.
             * 
             * @return Integer 
             */
            Integer top() const { return this->_top; }

            /**
             * @brief Bottom element.
             * 
             * @return Integer 
             */
            Integer bottom() const { return this->_bottom; }

            /**
             * @brief Facing elements.
             * 
             * @return std::vector<std::array<Integer, 2>> 
             */
            std::vector<std::array<Integer, 2>> facing() const { return this->_facing; }

            // Constructors.

            Neighbour21(const Integer &, const Integer &, const std::vector<std::array<Integer, 2>> &);
    };

}

#endif