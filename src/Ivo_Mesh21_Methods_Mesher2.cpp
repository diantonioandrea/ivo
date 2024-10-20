/**
 * @file Mesh21_Methods_Mesher2.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Mesh21/Methods/Mesher2.hpp implementation.
 * @date 2024-10-20
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Ivo.hpp>

namespace ivo {

    /**
     * @brief Reads a space diagram from a file.
     * 
     * @param filename Filename.
     * @return std::vector<Polygon21> 
     */
    std::vector<Polygon21> mesher2(const std::string &filename) {

        // File.
        std::ifstream file{filename};
        
        return {};
    }

    /**
     * @brief Writes a space diagram to a file.
     * 
     * @param filename Filename.
     * @param diagram Space diagram.
     */
    void mesher2(const std::string &filename, const std::vector<Polygon21> &diagram) {

        // File.
        std::ofstream file{filename};

        // Header.
        file << "@ Readable space diagram." << std::endl;
        file << "@ " << diagram.size() << " cells." << std::endl;

        // Writing.
        for(const auto &cell: diagram)
            file << cell << std::endl;
    }

}