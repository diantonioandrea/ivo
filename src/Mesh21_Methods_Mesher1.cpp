/**
 * @file Mesh21_Methods_Mesher1.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Mesh21/Methods/Mesher1.hpp implementation.
 * @date 2024-10-19
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Ivo.hpp>

namespace ivo {

    /**
     * @brief 1D mesher.
     * 
     * @param a Start.
     * @param b End.
     * @param n Steps.
     * @return std::vector<Real> 
     */
    std::vector<Real> mesher1(const Real &a, const Real &b, const Natural &n) {
        #ifndef NDEBUG // Integrity check.
        assert(a < b);
        assert(n > 0);
        #endif

        std::vector<Real> intervals;
        intervals.resize(n + 1, a);

        Real step = (b - a) / static_cast<Real>(n);

        for(Natural j = 1; j <= n; ++j)
            intervals[j] += j * step;

        return intervals;
    }

}