/**
 * @file Test_Mesh.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Simple mesh testing.
 * @date 2024-08-06
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <fstream>
#include <Ivo.hpp>

int main(int argc, char **argv) {

    // Points.
    ivo::Point21 a{0.0, 0.0};
    ivo::Point21 b{1.0, 0.0};
    ivo::Point21 c{1.0, 1.0};
    ivo::Point21 d{0.0, 1.0};

    // Polygon.
    ivo::Polygon21 abcd = {a, b, c, d};

    // Space diagram.
    std::vector<ivo::Polygon21> space = ivo::mesher2(abcd, 10);

    // Time "diagram" (intervals).
    std::vector<ivo::Real> time{0.0L, 0.25L, 0.5L, 0.75L, 1.0L};

    // Mesh.
    ivo::Mesh21 mesh{space, time};

    // Output.
    std::ofstream output("output/Test_Mesh21.p21");
    output << mesh;

    return 0;
}