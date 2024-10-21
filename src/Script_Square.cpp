/**
 * @file Script_Square.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Square domain mesher.
 * @date 2024-10-20
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include "./include/Square.hpp"

int main(int argc, char **argv) {

    if(argc != 2) {
        std::cout << "Usage: " << argv[0] << " SPACE_ELEMENTS_NUMBER [Ns]." << std::endl;
        return -1;
    }

    // Elements number.
    const ivo::Natural Ns = static_cast<ivo::Natural>(std::atoi(argv[1]));

    // Diagram.
    const std::vector<ivo::Polygon21> diagram = ivo::mesher2(ivo::square::abcd, Ns);

    // Output.
    ivo::mesher2("output/Square_" + std::to_string(Ns) + ".p2", diagram);

    return 0;
}