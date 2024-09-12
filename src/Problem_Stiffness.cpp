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

        // Stiffness submatrices.
        Sparse<Real> V{mesh.dofs(), mesh.dofs()}; // Volume integrals.

        // Loop over elements.
        for(Natural j = 0; j < mesh.space() * mesh.time(); ++j) {

            // Element.
            Element21 element = mesh.element(j);

            // Dofs.
            std::vector<Natural> dofs_j = mesh.dofs(j);
            Natural dofs_xy = (element.p() + 1) * (element.p() + 2) / 2;
            Natural dofs_t = element.q() + 1;

            // Equation coefficients.
            auto [convection_x, convection_y] = equation.convection(nodes1t_j);
            Vector<Real> diffusion = equation.diffusion(nodes1t_j);
            Vector<Real> reaction = equation.reaction(nodes1t_j);

            // VOLUME INTEGRALS - COMPUTING.

            // Submatrices.
            Matrix<Real> V_a_xy{dofs_xy, dofs_xy};
            Matrix<Real> V_a_t{dofs_t, dofs_t};

            Matrix<Real> V_bx_xy{dofs_xy, dofs_xy};
            Matrix<Real> V_by_xy{dofs_xy, dofs_xy};
            Matrix<Real> V_bx_t{dofs_t, dofs_t};
            Matrix<Real> V_by_t{dofs_t, dofs_t};

            Matrix<Real> V_c_xy{dofs_xy, dofs_xy};
            Matrix<Real> V_c_t{dofs_t, dofs_t};

            // Nodes and basis.
            auto [nodes1t_j, dt_j] = internal::reference_to_element(mesh, j, nodes1t);
            auto [nodes2xy_j, dxy_j] = internal::reference_to_element(mesh, j, {nodes2x, nodes2y});

            auto [phi_s, gradx_phi_s, grady_phi_s] = basis_s(mesh, j, nodes2xy_j);
            auto [phi_t, gradt_phi_t] = basis_t(mesh, j, nodes1t);

            Vector<Real> weights1_j = weights1 * dt_j;
            Vector<Real> weights2_j = weights2 * dxy_j;

            // Integrals.

            // a(*, *), diffusion.

            V_a_xy += internal::c_scale(weights2_j, gradx_phi_s).transpose() * gradx_phi_s + internal::c_scale(weights2_j, grady_phi_s).transpose() * grady_phi_s;
            V_a_t += internal::c_scale(weights1_j * diffusion, phi_t).transpose() * phi_t;

            // b(*, *), convection.

            V_bx_xy += internal::c_scale(weights2_j, gradx_phi_s).transpose() * phi_s;
            V_by_xy += internal::c_scale(weights2_j, grady_phi_s).transpose() * phi_s;
            V_bx_t += internal::c_scale(weights1_j * convection_x, phi_t).transpose() * phi_t;
            V_by_t += internal::c_scale(weights1_j * convection_y, phi_t).transpose() * phi_t;

            // c(*, *), reaction.

            V_c_xy += internal::c_scale(weights2_j,phi_s).transpose() * phi_s;
            V_c_t += internal::c_scale(weights1_j * reaction, phi_t).transpose() * phi_t;

            // VOLUME INTEGRALS - BUILDING.

            V(dofs_j, dofs_j, V(dofs_j, dofs_j) + kronecker(V_a_t, V_a_xy));
            V(dofs_j, dofs_j, V(dofs_j, dofs_j) + kronecker(V_bx_t, V_bx_xy) + kronecker(V_by_t, V_by_xy));
            V(dofs_j, dofs_j, V(dofs_j, dofs_j) + kronecker(V_c_t, V_c_xy));
            
            // FACE INTEGRALS - COMPUTING.

            // [!]

            // FACE INTEGRALS - BUILDING.

            // [!]
        }

        // Building.
        A = V;
    
        return A;
    }

}