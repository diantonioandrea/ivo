/**
 * @file Test_Fem.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Simple FEM testing.
 * @date 2024-08-05
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Ivo.hpp>

int main(int argc, char **argv) {

    // Gauss-Legendre nodes.
    auto [nodes, weights] = ivo::__gauss1(5, -1.0L, 1.0L);

    std::cout << nodes << std::endl;
    std::cout << weights << std::endl;

    return 0;
}