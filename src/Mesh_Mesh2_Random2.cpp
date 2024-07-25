/**
 * @file Mesh_Mesh2_Random2.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Mesh/Mesh2/Random2.hpp implementation.
 * @date 2024-07-25
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Ivo.hpp>

namespace ivo {

    // Random points in a polygon.

    /**
     * @brief Generates random points in a polygon
     * 
     * @param polygon Polygon.
     * @param number Number.
     * @return std::vector<Point21> 
     */
    std::vector<Point21> __random2(const Polygon21 &polygon, const Natural &number) {

        // Seeding.
        std::srand(std::time(nullptr));

        #ifndef NDEBUG // integrity check.
        assert(number > 0);
        assert(spatial(polygon));
        #endif

        // Random points.
        std::vector<Point21> points(number, Point21{});

        // Boundaries.
        auto [min_xy, max_xy] = box2(polygon);
        Real t = min_xy(2);

        // Generation.
        for(Natural j = 0; j < number; ++j) {
            bool generation = true;

            while(generation) {
                Real x = min_xy(0) + (max_xy(0) - min_xy(0)) * (static_cast<Real>(std::rand()) / static_cast<Real>(RAND_MAX));
                Real y = min_xy(1) + (max_xy(1) - min_xy(1)) * (static_cast<Real>(std::rand()) / static_cast<Real>(RAND_MAX));
                Point21 point{x, y, t};

                if(!contains2(polygon, point))
                    continue;
                
                generation = false;
                for(Natural k = 0; k < j; ++k)
                    if(point == points[k]) {
                        generation = true;
                        break;
                    }

                if(!generation)
                    points[j] = point;
            }
        }

        return points;
    }

}