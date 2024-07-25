/**
 * @file Voronoi2.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Voronoi. Space only.
 * @date 2024-07-24
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef MESH_MESH2_VORONOI2
#define MESH_MESH2_VORONOI2

#include "./Includes.hpp"
#include "./Random2.hpp"

namespace ivo {

    // Random Voronoi.

    std::vector<Polygon21> __voronoi2(const Polygon21 &, const Natural &);

}

#endif