/**
 * @file Mesh_Mesh2_Voronoi2.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-07-25
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Ivo.hpp>

namespace ivo {

    // Random Voronoi.

    /**
     * @brief Returns the Voronoi diagram of a random set of points inside a polygon.
     * 
     * @param polygon Polygon.
     * @param number Number.
     * @return std::vector<Polygon21> 
     */
    std::vector<Polygon21> __voronoi2(const Polygon21 &polygon, const Natural &number) {

        // Points and cells.
        std::vector<Point21> points = __random2(polygon, number);
        std::vector<Polygon21> cells; 

        // [!]
    }

}