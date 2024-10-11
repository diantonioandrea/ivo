/**
 * @file Test_Quadrature.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Simple quadrature testing.
 * @date 2024-09-04
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Ivo.hpp>

int main(int argc, char **argv) {

    auto [nodes2t_x, nodes2t_y, weights2t] = ivo::quadrature2xy(3);

    std::cout << "Nodes (x): " << nodes2t_x << std::endl;
    std::cout << "Nodes (y): " << nodes2t_y << std::endl;
    std::cout << "Weights: " << weights2t << std::endl;

    return 0;
}