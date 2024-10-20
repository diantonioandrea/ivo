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

        #ifndef NVERBOSE
        std::cout << "[Ivo] Mesher2" << std::endl;
        std::cout << "\t[Mesher2] Reading a diagram from: " << filename << std::endl;
        #endif

        // File.
        std::ifstream file{filename};

        // Diagram.
        std::vector<Polygon21> diagram;

        // Reading.
        std::string line;

        // Reading loop.
        while(std::getline(file, line)) {

            // Skip lines starting with '@'
            if (!line.empty() && line[0] == '@') {
                continue;
            }

            // Reading inside the line.
            std::istringstream lineStream{line};

            // Points.
            std::vector<Point21> points;
            Real x, y, t;

            while(lineStream >> x >> y >> t) {
                points.emplace_back(x, y, t);
            }

            // Diagram update.
            diagram.emplace_back(points);
        }

        #ifndef NVERBOSE
        std::cout << "\t[Mesher2] Exited" << std::endl;
        #endif
        
        return diagram;
    }

    /**
     * @brief Writes a space diagram to a file.
     * 
     * @param filename Filename.
     * @param diagram Space diagram.
     */
    void mesher2(const std::string &filename, const std::vector<Polygon21> &diagram) {

        #ifndef NVERBOSE
        std::cout << "[Ivo] Mesher2" << std::endl;
        std::cout << "\t[Mesher2] Writing a diagram to: " << filename << std::endl;
        #endif

        // File.
        std::ofstream file{filename};

        // Header.
        file << "@ Readable space diagram." << std::endl;
        file << "@ " << diagram.size() << " cells." << std::endl;

        // Writing.
        for(const auto &polygon: diagram) {
            for(const auto &point: polygon.points()) {
                const Real x = point(0);
                const Real y = point(1);
                const Real t = point(2);

                file << std::setprecision(14) << x << " " << y << " " << t << " " << std::flush;
            }

            file << std::endl;
        }

        #ifndef NVERBOSE
        std::cout << "\t[Mesher2] Exited" << std::endl;
        #endif
    }

}