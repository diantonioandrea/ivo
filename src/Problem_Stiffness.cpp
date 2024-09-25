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
     * @param initial Initial condition.
     * @return Sparse<Real> 
     */
    Sparse<Real> stiffness(const Mesh21 &mesh, const Equation &equation) {

        // Quadrature.
        auto [nodes1t, weights1] = quadrature1(constants::quadrature);
        auto [nodes2x, nodes2y, weights2] = quadrature2(constants::quadrature);

        // Stiffness submatrices.
        Sparse<Real> T{mesh.dofs(), mesh.dofs()}; // Volume integrals, time.
        Sparse<Real> E{mesh.dofs(), mesh.dofs()}; // Face integrals, time.
        Sparse<Real> V{mesh.dofs(), mesh.dofs()}; // Volume integrals.
        Sparse<Real> I{mesh.dofs(), mesh.dofs()}; // Face integrals.

        #ifndef NVERBOSE
        std::cout << "[Ivo] Stiffness" << std::endl;
        std::cout << "\t[Stiffness] Building the stiffness matrix" << std::endl;
        #endif

        // Loop over elements.
        for(Natural j = 0; j < mesh.space() * mesh.time(); ++j) {

            // ELEMENT DATA.

            // Element.
            Element21 element = mesh.element(j);

            // Dofs.
            std::vector<Natural> dofs_j = mesh.dofs(j);
            Natural dofs_xy = (element.p() + 1) * (element.p() + 2) / 2;
            Natural dofs_t = element.q() + 1;
            Natural dofs_xyt = dofs_t * dofs_xy;

            // Neighbours.
            Neighbour21 neighbourhood = mesh.neighbour(j);

            std::vector<std::array<Integer, 2>> facing = neighbourhood.facing();
            Natural neighbours = facing.size();

            Integer bottom = neighbourhood.bottom();

            // VOLUME INTEGRALS - PRECOMPUTING.

            // Submatrices.
            Matrix<Real> T_xyt{dofs_xyt, dofs_xyt};
            Matrix<Real> V_a_xyt{dofs_xyt, dofs_xyt};
            Matrix<Real> V_b_xyt{dofs_xyt, dofs_xyt};
            Matrix<Real> V_c_xyt{dofs_xyt, dofs_xyt};

            // Nodes and basis, time.
            auto [nodes1t_j, dt_j] = internal::reference_to_element(mesh, j, nodes1t);
            auto [phi_t, gradt_phi_t] = basis_t(mesh, j, nodes1t);

            // Weights, time.
            Vector<Real> weights1_j = weights1 * dt_j;

            // VOLUME INTEGRALS - COMPUTING.

            for(Natural k = 0; k < neighbours; ++k) { // Sub-triangulation.

                // Nodes and basis, space.
                auto [nodes2xy_j, dxy_j] = internal::reference_to_element(mesh, j, k, {nodes2x, nodes2y});
                auto [phi_s, gradx_phi_s, grady_phi_s] = basis_s(mesh, j, nodes2xy_j);

                // Weights, space.
                Vector<Real> weights2_j = weights2 * dxy_j;

                for(Natural jt = 0; jt < dofs_t; ++jt)
                    for(Natural ht = 0; ht < dofs_t; ++ht)
                        for(Natural jxy = 0; jxy < dofs_xy; ++jxy)
                            for(Natural hxy = 0; hxy < dofs_xy; ++hxy)
                                for(Natural kt = 0; kt < phi_t.rows(); ++kt)
                                    for(Natural kxy = 0; kxy < phi_s.rows(); ++kxy) { // Brute-force integral.
                                        Real t = nodes1t_j(kt);

                                        // Equation coefficients.
                                        auto [convection_x, convection_y] = equation.convection(t);
                                        Real diffusion = equation.diffusion(t);
                                        Real reaction = equation.reaction(t);

                                        // (*', *).

                                        T_xyt(jt * dofs_xy + jxy, ht * dofs_xy + hxy, T_xyt(jt * dofs_xy + jxy, ht * dofs_xy + hxy) + weights2_j(kxy) * weights1_j(kt) * gradt_phi_t(kt, jt) * phi_s(kxy, jxy) * phi_t(kt, ht) * phi_s(kxy, hxy));

                                        // a(*, *), diffusion.

                                        V_a_xyt(jt * dofs_xy + jxy, ht * dofs_xy + hxy, V_a_xyt(jt * dofs_xy + jxy, ht * dofs_xy + hxy) + weights2_j(kxy) * weights1_j(kt) * (phi_t(kt, jt) * gradx_phi_s(kxy, jxy) * phi_t(kt, ht) * gradx_phi_s(kxy, hxy) + phi_t(kt, jt) * grady_phi_s(kxy, jxy) * phi_t(kt, ht) * grady_phi_s(kxy, hxy)) * diffusion);

                                        // b(*, *), convection.

                                        V_b_xyt(jt * dofs_xy + jxy, ht * dofs_xy + hxy, V_b_xyt(jt * dofs_xy + jxy, ht * dofs_xy + hxy) + weights2_j(kxy) * weights1_j(kt) * (phi_t(kt, jt) * gradx_phi_s(kxy, jxy) * convection_x + phi_t(kt, jt) * grady_phi_s(kxy, jxy) * convection_y) * phi_t(kt, ht) * phi_s(kxy, hxy));

                                        // c(*, *), reaction.

                                        V_c_xyt(jt * dofs_xy + jxy, ht * dofs_xy + hxy, V_c_xyt(jt * dofs_xy + jxy, ht * dofs_xy + hxy) + weights2_j(kxy) * weights1_j(kt) * phi_t(kt, jt) * phi_s(kxy, jxy) * phi_t(kt, ht) * phi_s(kxy, hxy) * reaction);
                                    }
            }

            // VOLUME INTEGRALS - BUILDING.

            T(dofs_j, dofs_j, T(dofs_j, dofs_j) + T_xyt);
            V(dofs_j, dofs_j, V(dofs_j, dofs_j) + V_a_xyt);
            V(dofs_j, dofs_j, V(dofs_j, dofs_j) + V_b_xyt);
            V(dofs_j, dofs_j, V(dofs_j, dofs_j) + V_c_xyt);

            // FACE INTEGRALS - PRECOMPUTING.

            // Submatrices.
            std::vector<Matrix<Real>> I_a;
            std::vector<Matrix<Real>> I_b;
            std::vector<Matrix<Real>> I_J;
            
            for(Natural k = 0; k < neighbours; ++k) {

                // Nodes and basis, space.
                auto [e_nodes2xy_j, normal, e_dxy_j] = internal::reference_to_element(mesh, j, k, nodes1t);
                auto [e_phi_s, e_gradx_phi_s, e_grady_phi_s] = basis_s(mesh, j, e_nodes2xy_j);

                // Weights, space.
                Vector<Real> e_weights2_j = weights1 * e_dxy_j;

                // Submatrices.
                Matrix<Real> I_a_xyt{dofs_xyt, dofs_xyt};
                Matrix<Real> I_b_xyt{dofs_xyt, dofs_xyt};
                Matrix<Real> I_J_xyt{dofs_xyt, dofs_xyt};

                if(facing[k][0] != -1) {

                    // Neighbour element.
                    Element21 n_element = mesh.element(facing[k][0]);

                    // Neighbour basis.
                    auto [n_e_phi_s, n_e_gradx_phi_s, n_e_grady_phi_s] = basis_s(mesh, facing[k][0], e_nodes2xy_j);
                    auto [n_phi_t, n_gradt_phi_t] = basis_t(mesh, facing[k][0], nodes1t);

                    // FACE INTEGRALS - COMPUTING.

                    for(Natural jt = 0; jt < dofs_t; ++jt)
                        for(Natural ht = 0; ht < dofs_t; ++ht)
                            for(Natural jxy = 0; jxy < dofs_xy; ++jxy)
                                for(Natural hxy = 0; hxy < dofs_xy; ++hxy)
                                    for(Natural kt = 0; kt < phi_t.rows(); ++kt)
                                        for(Natural kxy = 0; kxy < e_phi_s.rows(); ++kxy) { // Brute-force integral.
                                            Real t = nodes1t_j(kt);

                                            #ifndef NDEBUG // Integrity check.
                                            assert(phi_t.rows() == n_phi_t.rows());
                                            assert(e_phi_s.rows() == n_e_phi_s.rows());
                                            #endif

                                            // Equation coefficients.
                                            auto [convection_x, convection_y] = equation.convection(t);
                                            Real convection_n = normal(0) * convection_x + normal(1) * convection_y;
                                            Real diffusion = equation.diffusion(t);

                                            // Sign check.
                                            Real negative = (convection_n < 0.0L) ? 1.0L : 0.0L;

                                            // a(*, *), diffusion.

                                            I_a_xyt(jt * dofs_xy + jxy, ht * dofs_xy + hxy, I_a_xyt(jt * dofs_xy + jxy, ht * dofs_xy + hxy) + weights1_j(kt) * e_weights2_j(kxy) * (0.5L * normal(0) * (phi_t(kt, jt) * e_gradx_phi_s(kxy, jxy) + n_phi_t(kt, jt) * n_e_gradx_phi_s(kxy, jxy)) + 0.5L * normal(1) * (phi_t(kt, jt) * e_grady_phi_s(kxy, jxy) + n_phi_t(kt, jt) * n_e_grady_phi_s(kxy, jxy))) * (phi_t(kt, ht) * e_phi_s(kxy, hxy) - n_phi_t(kt, ht) * n_e_phi_s(kxy, hxy)) * diffusion);
                                            I_a_xyt(jt * dofs_xy + jxy, ht * dofs_xy + hxy, I_a_xyt(jt * dofs_xy + jxy, ht * dofs_xy + hxy) - weights1_j(kt) * e_weights2_j(kxy) * (0.5L * normal(0) * (phi_t(kt, ht) * e_gradx_phi_s(kxy, hxy) + n_phi_t(kt, ht) * n_e_gradx_phi_s(kxy, hxy)) + 0.5L * normal(1) * (phi_t(kt, ht) * e_grady_phi_s(kxy, hxy) + n_phi_t(kt, ht) * n_e_grady_phi_s(kxy, hxy))) * (phi_t(kt, jt) * e_phi_s(kxy, jxy) - n_phi_t(kt, jt) * n_e_phi_s(kxy, jxy)) * diffusion);

                                            // b(*, *), convection.

                                            I_b_xyt(jt * dofs_xy + jxy, ht * dofs_xy + hxy, I_b_xyt(jt * dofs_xy + jxy, ht * dofs_xy + hxy) + negative * weights1_j(kt) * e_weights2_j(kxy) * (phi_t(kt, jt) * e_phi_s(kxy, jxy) - n_phi_t(kt, jt) * n_e_phi_s(kxy, jxy)) * phi_t(kt, ht) * e_phi_s(kxy, hxy) * convection_n);

                                            // J(*, *). Mind the signs.

                                            I_J_xyt(jt * dofs_xy + jxy, ht * dofs_xy + hxy, I_J_xyt(jt * dofs_xy + jxy, ht * dofs_xy + hxy) - weights1_j(kt) * e_weights2_j(kxy) / e_dxy_j * (phi_t(kt, jt) * e_phi_s(kxy, jxy) - n_phi_t(kt, jt) * n_e_phi_s(kxy, jxy)) * (phi_t(kt, ht) * e_phi_s(kxy, hxy) - n_phi_t(kt, ht) * n_e_phi_s(kxy, hxy)) * diffusion);
                                        }
                } else {

                    // FACE INTEGRALS - COMPUTING.

                    for(Natural jt = 0; jt < dofs_t; ++jt)
                        for(Natural ht = 0; ht < dofs_t; ++ht)
                            for(Natural jxy = 0; jxy < dofs_xy; ++jxy)
                                for(Natural hxy = 0; hxy < dofs_xy; ++hxy)
                                    for(Natural kt = 0; kt < phi_t.rows(); ++kt)
                                        for(Natural kxy = 0; kxy < e_phi_s.rows(); ++kxy) { // Brute-force integral.
                                            Real t = nodes1t_j(kt);

                                            // Equation coefficients.
                                            auto [convection_x, convection_y] = equation.convection(t);
                                            Real convection_n = normal(0) * convection_x + normal(1) * convection_y;
                                            Real diffusion = equation.diffusion(t);

                                            // Sign check.
                                            Real negative = (convection_n < 0.0L) ? 1.0L : 0.0L;

                                            // a(*, *), diffusion.

                                            I_a_xyt(jt * dofs_xy + jxy, ht * dofs_xy + hxy, I_a_xyt(jt * dofs_xy + jxy, ht * dofs_xy + hxy) + negative * weights1_j(kt) * e_weights2_j(kxy) * (normal(0) * phi_t(kt, jt) * e_gradx_phi_s(kxy, jxy) + normal(1) * phi_t(kt, jt) * e_grady_phi_s(kxy, jxy)) * phi_t(kt, ht) * e_phi_s(kxy, hxy) * diffusion);
                                            I_a_xyt(jt * dofs_xy + jxy, ht * dofs_xy + hxy, I_a_xyt(jt * dofs_xy + jxy, ht * dofs_xy + hxy) - negative * weights1_j(kt) * e_weights2_j(kxy) * (normal(0) * phi_t(kt, ht) * e_gradx_phi_s(kxy, hxy) + normal(1) * phi_t(kt, ht) * e_grady_phi_s(kxy, hxy)) * phi_t(kt, jt) * e_phi_s(kxy, jxy) * diffusion);

                                            // b(*, *), convection.

                                            I_b_xyt(jt * dofs_xy + jxy, ht * dofs_xy + hxy, I_b_xyt(jt * dofs_xy + jxy, ht * dofs_xy + hxy) + negative * weights1_j(kt) * e_weights2_j(kxy) * phi_t(kt, jt) * e_phi_s(kxy, jxy) * phi_t(kt, ht) * e_phi_s(kxy, hxy) * convection_n);

                                            // J(*, *).

                                            I_J_xyt(jt * dofs_xy + jxy, ht * dofs_xy + hxy, I_J_xyt(jt * dofs_xy + jxy, ht * dofs_xy + hxy) - negative * weights1_j(kt) * e_weights2_j(kxy) / e_dxy_j * phi_t(kt, jt) * e_phi_s(kxy, jxy) * phi_t(kt, ht) * e_phi_s(kxy, hxy) * diffusion);
                                        }
                }

                // FACE INTEGRALS - PREBUILDING.

                I_a.emplace_back(I_a_xyt);
                I_b.emplace_back(I_b_xyt);
                I_J.emplace_back(I_J_xyt);
            }

            // FACE INTEGRALS - BUILDING.

            for(Natural k = 0; k < neighbours; ++k) {
                if(facing[k][0] != -1) {

                    // Neighbour dofs.
                    std::vector<Natural> n_dofs_j = mesh.dofs(facing[k][0]);

                    // Building.
                    I(dofs_j, n_dofs_j, I(dofs_j, n_dofs_j) + I_a[k] + I_b[k] + I_J[k]);
                    
                } else {
                    
                    // Building.
                    I(dofs_j, dofs_j, I(dofs_j, dofs_j) + I_a[k] + I_b[k] + I_J[k]);
                }
            }

            // TIME FACE INTEGRALS - PRECOMPUTING.

            // Face time basis.
            auto [f_phi_t, f_gradt_phi_t] = basis_t(mesh, j, Vector<Real>{1, -1.0L}); // [?]

            if(bottom != -1) {

                // Neighbour element.
                Element21 n_element = mesh.element(bottom);

                // Neighbour face time basis.
                auto [n_f_phi_t, n_f_gradt_phi_t] = basis_t(mesh, bottom, Vector<Real>{1, 1.0L}); // [?]

                // Neighbour dofs.
                std::vector<Natural> n_dofs_j = mesh.dofs(bottom);

                // Submatrices.
                Matrix<Real> E_xy{dofs_xy, dofs_xy};
                Matrix<Real> E_t{dofs_t, dofs_t};

                // TIME FACE INTEGRALS - COMPUTING.

                for(Natural k = 0; k < neighbours; ++k) { // Sub-triangulation.

                    // Nodes and basis.
                    auto [nodes2xy_j, dxy_j] = internal::reference_to_element(mesh, j, k, {nodes2x, nodes2y});
                    auto [phi_s, gradx_phi_s, grady_phi_s] = basis_s(mesh, j, nodes2xy_j);

                    Vector<Real> weights2_j = weights2 * dxy_j;

                    E_xy += internal::c_scale(weights2_j, phi_s).transpose() * phi_s;
                }

                E_t = (f_phi_t - n_f_phi_t).transpose() * f_phi_t;

                // TIME FACE INTEGRALS - BUILDING.

                E(n_dofs_j, dofs_j, E(n_dofs_j, dofs_j) + kronecker(E_t, E_xy)); // [??]

            } else {

                // Submatrices.
                Matrix<Real> E_xy{dofs_xy, dofs_xy};
                Matrix<Real> E_t{dofs_t, dofs_t};

                // TIME FACE INTEGRALS - COMPUTING.

                for(Natural k = 0; k < neighbours; ++k) { // Sub-triangulation.

                    // Nodes and basis.
                    auto [nodes2xy_j, dxy_j] = internal::reference_to_element(mesh, j, k, {nodes2x, nodes2y});
                    auto [phi_s, gradx_phi_s, grady_phi_s] = basis_s(mesh, j, nodes2xy_j);

                    Vector<Real> weights2_j = weights2 * dxy_j;

                    E_xy += internal::c_scale(weights2_j, phi_s).transpose() * phi_s;
                }

                E_t = f_phi_t.transpose() * f_phi_t;

                // TIME FACE INTEGRALS - BUILDING.

                E(dofs_j, dofs_j, E(dofs_j, dofs_j) + kronecker(E_t, E_xy)); // [?]
            }
        }

        #ifndef NVERBOSE
        std::cout << "\t[Stiffness] Exited" << std::endl;
        #endif

        // Building and return.
        return T + E + V - I;
    }

}