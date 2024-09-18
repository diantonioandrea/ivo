/**
 * @file Problem_Forcing.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Problem/Forcing.hpp implementation.
 * @date 2024-09-14
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Ivo.hpp>

namespace ivo {

    /**
     * @brief Builds the forcing vector for a 2+1D equation.
     * 
     * @param mesh Mesh.
     * @param equation Equation.
     * @param data Initial data.
     * @return Vector<Real> 
     */
    Vector<Real> forcing(const Mesh21 &mesh, const Equation &equation, const Data &data, const Initial &initial) {

        // Quadrature.
        auto [nodes1t, weights1] = quadrature1(constants::quadrature);
        auto [nodes2x, nodes2y, weights2] = quadrature2(constants::quadrature);

        // Forcing vector.
        Vector<Real> V{mesh.dofs()}; // Volume integrals.
        Vector<Real> I{mesh.dofs()}; // Face integrals.

        // Loop over elements.
        for(Natural j = 0; j < mesh.space() * mesh.time(); ++j) {

            // [!]
            
        }

        // Building and return
        return V + I;
    }

}