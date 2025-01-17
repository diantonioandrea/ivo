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
    Vector<Real> forcing(const Mesh21 &mesh, const Equation &equation, const Data &data) {

        // Quadrature.
        auto [nodes1t, weights1t] = quadrature1t(constants::quadrature);
        auto [nodes1x, weights1x] = quadrature1x(constants::quadrature);
        auto [nodes2x, nodes2y, weights2] = quadrature2xy(constants::quadrature);

        // Forcing vector.
        Vector<Real> V{mesh.dofs()}; // Volume integrals.
        Vector<Real> I{mesh.dofs()}; // Face integrals.

        #ifndef NVERBOSE
        std::cout << "[Ivo] Forcing" << std::endl;
        std::cout << "\t[Forcing] Building the forcing vector" << std::endl;
        #endif

        // Loop over elements.
        for(Natural j = 0; j < mesh.space() * mesh.time(); ++j) {

            // Element.
            Element21 element = mesh.element(j);

            // Neighbours.
            Neighbour21 neighbourhood = mesh.neighbour(j);

            std::vector<std::array<Integer, 2>> facing = neighbourhood.facing();
            Natural neighbours = facing.size();

            // Dofs.
            std::vector<Natural> dofs_j = mesh.dofs(j);
            Natural dofs_xy = (element.p() + 1) * (element.p() + 2) / 2;
            Natural dofs_t = element.q() + 1;
            Natural dofs_xyt = dofs_t * dofs_xy;

            // VOLUME INTEGRALS - PRECOMPUTING.

            // Subvector.
            Vector<Real> V_xyt{dofs_xy * dofs_t};

            // Nodes and basis, time.
            auto [nodes1t_j, dt_j] = internal::reference_to_element(mesh, j, nodes1t);
            auto [phi_t, gradt_phi_t] = basis_t(mesh, j, nodes1t_j);

            Vector<Real> weights1t_j = weights1t * dt_j;

            // VOLUME INTEGRALS - COMPUTING.

            for(Natural k = 0; k < neighbours; ++k) { // Sub-triangulation.

                // Nodes and basis, space.
                auto [nodes2xy_j, dxy_j] = internal::reference_to_element(mesh, j, k, {nodes2x, nodes2y});
                auto [phi_xy, gradx_phi_xy, grady_phi_xy] = basis_xy(mesh, j, nodes2xy_j);
                auto [nodes2x_j, nodes2y_j] = nodes2xy_j;

                Vector<Real> weights2_j = weights2 * dxy_j;

                // CURRENT.

                for(Natural ht = 0; ht < dofs_t; ++ht)
                    for(Natural hxy = 0; hxy < dofs_xy; ++hxy) {
                        Real c_V_xyt = 0.0L;

                        for(Natural kt = 0; kt < phi_t.rows(); ++kt)
                            for(Natural kxy = 0; kxy < phi_xy.rows(); ++kxy) { // Brute-force integral.
                                Real x = nodes2x_j(kxy);
                                Real y = nodes2y_j(kxy);
                                Real t = nodes1t_j(kt);

                                // Data.
                                Real source = data.source(x, y, t);

                                // Source.

                                c_V_xyt += weights2_j(kxy) * weights1t_j(kt) * phi_t(kt, ht) * phi_xy(kxy, hxy) * source;
                            }

                        V_xyt(ht * dofs_xy + hxy, V_xyt(ht * dofs_xy + hxy) + c_V_xyt);
                    }
            }

            // VOLUME INTEGRALS - BUILDING.

            V(dofs_j, V(dofs_j) + V_xyt);

            // FACE INTEGRALS - PRECOMPUTING.
            
            for(Natural k = 0; k < neighbours; ++k) {

                // Boundary edges.
                if(facing[k][0] != -1)
                    continue;

                // Nodes and basis, space.
                auto [e_nodes2xy_j, normal, e_dxy_j] = internal::reference_to_element(mesh, j, k, nodes1x);
                auto [e_phi_xy, e_gradx_phi_xy, e_grady_phi_xy] = basis_xy(mesh, j, e_nodes2xy_j);
                auto [e_nodes2x_j, e_nodes2y_j] = e_nodes2xy_j;

                // Normal gradient.
                Matrix<Real> e_gradn_phi_xy = normal(0) * e_gradx_phi_xy + normal(1) * e_grady_phi_xy;

                // Weights, space.
                Vector<Real> e_weights2_j = weights1x * e_dxy_j;

                // Subvectors.
                Vector<Real> I_de_xyt{dofs_xyt}; // [!]
                Vector<Real> I_d_xyt{dofs_xyt}; // [!]
                Vector<Real> I_n_xyt{dofs_xyt}; // [!]

                // FACE INTEGRALS - COMPUTING.

                // CURRENT.

                for(Natural ht = 0; ht < dofs_t; ++ht)
                    for(Natural hxy = 0; hxy < dofs_xy; ++hxy) {
                        Real c_de_xyt = 0.0L;
                        Real c_d_xyt = 0.0L;
                        Real c_n_xyt = 0.0L;

                        for(Natural kt = 0; kt < phi_t.rows(); ++kt)
                            for(Natural kxy = 0; kxy < e_phi_xy.rows(); ++kxy) { // Brute-force integral.
                                Real x = e_nodes2x_j(kxy);
                                Real y = e_nodes2y_j(kxy);
                                Real t = nodes1t_j(kt);

                                // Equation coefficients.
                                auto [convection_x, convection_y] = equation.convection(x, y, t);
                                Real convection_n = normal(0) * convection_x + normal(1) * convection_y;
                                Real diffusion = equation.diffusion();

                                // Data.
                                Real dirichlet = data.dirichlet(x, y, t);
                                Real neumann = data.neumann(x, y, t);

                                // Boundary check.
                                Real negative = (convection_n < 0.0L) ? 1.0L : 0.0L;
                                Real positive = (convection_n >= 0.0L) ? 1.0L : 0.0L;

                                // Dirichlet.

                                c_de_xyt += negative * e_weights2_j(kxy) / e_dxy_j * weights1t_j(kt) * phi_t(kt, ht) * e_phi_xy(kxy, hxy) * dirichlet * diffusion;
                                c_de_xyt += negative * e_weights2_j(kxy) * weights1t_j(kt) * phi_t(kt, ht) * e_gradn_phi_xy(kxy, hxy) * dirichlet * diffusion;

                                c_d_xyt -= negative * e_weights2_j(kxy) * weights1t_j(kt) * phi_t(kt, ht) * e_phi_xy(kxy, hxy) * dirichlet * convection_n;

                                // Neumann.

                                c_n_xyt += positive * e_weights2_j(kxy) * weights1t_j(kt) * phi_t(kt, ht) * e_phi_xy(kxy, hxy) * neumann;
                            }
                        
                        I_de_xyt(ht * dofs_xy + hxy, I_de_xyt(ht * dofs_xy + hxy) + c_de_xyt);
                        I_d_xyt(ht * dofs_xy + hxy, I_d_xyt(ht * dofs_xy + hxy) + c_d_xyt);
                        I_n_xyt(ht * dofs_xy + hxy, I_n_xyt(ht * dofs_xy + hxy) + c_n_xyt);
                    }

                // FACE INTEGRALS - BUILDING.

                I(dofs_j, I(dofs_j) + I_de_xyt + I_d_xyt + I_n_xyt);
            }

            #ifndef NVERBOSE
            if((j + 1) % mesh.space() == 0)
                std::cout << "\t[Forcing] Progress: " << j / mesh.space() + 1 << "/" << mesh.time() << std::endl;
            #endif
        }

        #ifndef NVERBOSE
        std::cout << "\t[Forcing] Exited" << std::endl;
        #endif

        // Building and return
        return V + I;
    }

}