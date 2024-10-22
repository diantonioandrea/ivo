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

    #ifndef NVERBOSE
    std::cout << "[Ivo] SCRIPT, Building a square diagram of " << Ns << " cells\n" << std::endl;
    #else
    std::cout << "[Ivo] SCRIPT, Building a square diagram of " << Ns << " cells" << std::endl;
    #endif

    // Diagram.
    const std::vector<ivo::Polygon21> diagram = ivo::mesher2(ivo::square::abcd, Ns);

    // File.
    const std::string filename = "output/Square_" + std::to_string(Ns) + ".p2";

    // Output.
    ivo::mesher2(filename, diagram);

    #ifndef NVERBOSE
    std::cout << "\n\t[SCRIPT] Diagram saved to " << filename << "\n" << std::endl;
    #else
    std::cout << "\t[SCRIPT] Diagram saved to " << filename << std::endl;
    #endif

    std::cout << "\t[SCRIPT] Exited" << std::endl;

    return 0;
}