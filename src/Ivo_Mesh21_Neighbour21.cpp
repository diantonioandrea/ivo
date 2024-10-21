/**
 * @file Mesh21_Neighbour21.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-07-25
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Ivo.hpp>

namespace ivo {

    // Constructors.

    /**
     * @brief Default constructor.
     * 
     * @param top Top element index.
     * @param bottom Bottom element index.
     * @param facing Facing elemens indices.
     */
    Neighbour21::Neighbour21(const Integer &top, const Integer &bottom, const std::vector<std::array<Integer, 2>> &facing): _top{top}, _bottom{bottom}, _facing{facing} {}

}