/**
 * @file Mesher1.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 2D mesher methods.
 * @date 2024-10-20
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef MESH21_MESH21_METHODS_MESHER2
#define MESH21_MESH21_METHODS_MESHER2

#include "../Includes.hpp"

namespace ivo {

    std::vector<Polygon21> mesher2(const std::string &);
    void mesher2(const std::string &, const std::vector<Polygon21> &);

}

#endif