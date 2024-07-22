/**
 * @file Test_Geometry.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Simple geometry testing.
 * @date 2024-07-22
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Ivo.hpp>

int main(int argc, char **argv) {

    // Points.
    ivo::Point21 a{0, 0, 0};
    ivo::Point21 b{0, 1, 0};
    ivo::Point21 c{1, 1, 0};

    // Edges.
    ivo::Edge21 ab{a, b};
    ivo::Edge21 bc{b, c};
    ivo::Edge21 ca{c, a};

    // Lines.
    ivo::Line21 abl{a, b};
    ivo::Line21 bcl{b, c};
    ivo::Line21 cal{c, a};

    // Polygon.
    ivo::Polygon21 abc{{a, b, c}};

    // Output.
    std::cout << a << std::endl << std::endl;
    std::cout << b << std::endl << std::endl;
    std::cout << c << std::endl << std::endl;
    std::cout << ab << std::endl << std::endl;
    std::cout << bc << std::endl << std::endl;
    std::cout << ca << std::endl << std::endl;
    std::cout << abl << std::endl << std::endl;
    std::cout << bcl << std::endl << std::endl;
    std::cout << cal << std::endl << std::endl;
    std::cout << abc << std::endl;

    return 0;
}