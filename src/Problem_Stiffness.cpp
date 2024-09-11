/**
 * @file Problem_Stiffness.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Problem/Stiffness.hpp implementation.
 * @date 2024-08-07
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Ivo.hpp>

namespace ivo {

    /**
     * @brief Builds the stiffness matrix for a 2+1D equation.
     * 
     * @param mesh Mesh.
     * @param equation Equation.
     * @return Matrix<Real> 
     */
    Sparse<Real> stiffness(const Mesh21 &mesh, const Equation &equation) {

        // Quadrature.
        auto [nodes1t, weights1] = quadrature1(constants::quadrature);
        auto [nodes2x, nodes2y, weights2] = quadrature2(constants::quadrature);

        // Stiffness matrix.
        Sparse<Real> A{mesh.dofs(), mesh.dofs()};

        // Loop over elements.
        for(Natural j = 0; j < mesh.space() * mesh.time(); ++j) {

            // Element.
            Element21 element = mesh.element(j);

            // Dofs.
            std::vector<Natural> dofs = mesh.dofs(j);

            // VOLUME INTEGRALS.

            // Nodes and basis.
            auto [nodes1t_j, dt_j] = internal::reference_to_element(mesh, j, nodes1t);
            auto [nodes2xy_j, dxy_j] = internal::reference_to_element(mesh, j, {nodes2x, nodes2y});

            auto [phi_s, gradx_phi_s, grady_phi_s] = basis_s(mesh, j, nodes2xy_j);
            auto [phi_t, gradt_phi_t] = basis_t(mesh, j, nodes1t);

            // [!]
            
            // FACE INTEGRALS.

            // [!]
        }
    
        return A;
    }

}