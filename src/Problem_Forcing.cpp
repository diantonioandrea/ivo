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
        Vector<Real> E{mesh.dofs()}; // Face integrals, time.
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

            Integer bottom = neighbourhood.bottom();

            // Dofs.
            std::vector<Natural> dofs_j = mesh.dofs(j);
            Natural dofs_xy = (element.p() + 1) * (element.p() + 2) / 2;
            Natural dofs_t = element.q() + 1;

            // VOLUME INTEGRALS - PRECOMPUTING.

            // Subvector.
            Vector<Real> V_xyt{dofs_xy * dofs_t};

            // Nodes and basis.
            auto [nodes1t_j, dt_j] = internal::reference_to_element(mesh, j, nodes1t);
            auto [phi_t, gradt_phi_t] = basis_t(mesh, j, nodes1t);

            Vector<Real> weights1_j = weights1 * dt_j;

            // Equation coefficients.
            auto [convection_x, convection_y] = equation.convection(nodes1t_j);
            Vector<Real> diffusion = equation.diffusion(nodes1t_j);
            Vector<Real> reaction = equation.reaction(nodes1t_j);

            // VOLUME INTEGRALS - COMPUTING.

            for(Natural k = 0; k < neighbours; ++k) { // Sub-triangulation.

                // Nodes and basis.
                auto [nodes2xy_j, dxy_j] = internal::reference_to_element(mesh, j, k, {nodes2x, nodes2y});
                auto [phi_s, gradx_phi_s, grady_phi_s] = basis_s(mesh, j, nodes2xy_j);
                auto [nodes2x_j, nodes2y_j] = nodes2xy_j;

                Vector<Real> weights2_j = weights2 * dxy_j;

                for(Natural jt = 0; jt < dofs_t; ++jt)
                    for(Natural jxy = 0; jxy < dofs_xy; ++jxy)
                        for(Natural kt = 0; kt < phi_t.rows(); ++kt)
                            for(Natural kxy = 0; kxy < phi_s.rows(); ++kxy) { // Brute-force integral.
                                Real x = nodes2x_j(kxy);
                                Real y = nodes2y_j(kxy);
                                Real t = nodes1t_j(kt);

                                V_xyt(jt * dofs_xy + jxy, V_xyt(jt * dofs_xy + jxy) + weights2_j(kxy) * weights1_j(kt) * phi_t(kt, jt) * phi_s(kxy, jxy) * data.source(x, y, t));
                            }
            }

            // VOLUME INTEGRALS - BUILDING.

            V(dofs_j, V(dofs_j) + V_xyt);

            // FACE INTEGRALS - PRECOMPUTING.

            // Subvectors.
            std::vector<Vector<Real>> I_d;
            std::vector<Vector<Real>> I_de;
            std::vector<Vector<Real>> I_n;
            
            for(Natural k = 0; k < neighbours; ++k) {

                // Nodes and basis. Using time nodes as 1D nodes.
                auto [e_nodes2xy_j, normal, e_dxy_j] = internal::reference_to_element(mesh, j, k, nodes1t);
                auto [e_phi_s, e_gradx_phi_s, e_grady_phi_s] = basis_s(mesh, j, e_nodes2xy_j);
                auto [e_nodes2x_j, e_nodes2y_j] = e_nodes2xy_j;

                Matrix<Real> e_gradn_phi_s = normal(0) * e_gradx_phi_s + normal(1) * e_grady_phi_s;
                Vector<Real> convection_n = normal(0) * convection_x + normal(1) * convection_y;

                Vector<Real> e_weights2_j = weights1 * e_dxy_j;

                // Sign check.
                Vector<Real> negative{nodes1t.size()};
                Vector<Real> positive{nodes1t.size()};

                for(Natural j = 0; j < nodes1t.size(); ++j) {
                    negative(j, (normal(0) * convection_x(j) + normal(1) * convection_y(j) < 0.0L) ? 1.0L : 0.0L);
                    positive(j, 1.0L - negative(j));
                }

                if(facing[k][0] == -1) {

                    // subvectors.
                    Vector<Real> I_de_xyt{dofs_t * dofs_xy};
                    Vector<Real> I_d_xyt{dofs_t * dofs_xy};
                    Vector<Real> I_n_xyt{dofs_t * dofs_xy};

                    // FACE INTEGRALS - COMPUTING.

                    // Dirichlet, diffusion.

                    for(Natural jt = 0; jt < dofs_t; ++jt)
                        for(Natural jxy = 0; jxy < dofs_xy; ++jxy)
                            for(Natural kt = 0; kt < phi_t.rows(); ++kt)
                                for(Natural kxy = 0; kxy < e_phi_s.rows(); ++kxy) { // Brute-force integral.
                                    Real x = e_nodes2x_j(kxy);
                                    Real y = e_nodes2y_j(kxy);
                                    Real t = nodes1t_j(kt);

                                    // Dirichlet, diffusion.

                                    I_de_xyt(jt * dofs_xy + jxy, V_xyt(jt * dofs_xy + jxy) + negative(kt) * e_weights2_j(kxy) / e_dxy_j * weights1_j(kt) * phi_t(kt, jt) * e_phi_s(kxy, jxy) * data.dirichlet(x, y, t));
                                    I_de_xyt(jt * dofs_xy + jxy, V_xyt(jt * dofs_xy + jxy) + negative(kt) * e_weights2_j(kxy) * weights1_j(kt) * phi_t(kt, jt) * e_gradn_phi_s(kxy, jxy) * data.dirichlet(x, y, t));

                                    // Dirichlet.

                                    I_d_xyt(jt * dofs_xy + jxy, V_xyt(jt * dofs_xy + jxy) + negative(kt) * e_weights2_j(kxy) * weights1_j(kt) * phi_t(kt, jt) * convection_n(kt) * e_phi_s(kxy, jxy) * data.dirichlet(x, y, t));

                                    // Neumann.

                                    I_n_xyt(jt * dofs_xy + jxy, V_xyt(jt * dofs_xy + jxy) + positive(kt) * e_weights2_j(kxy) * weights1_j(kt) * phi_t(kt, jt) * e_phi_s(kxy, jxy) * data.neumann(x, y, t));
                                }

                    // FACE INTEGRALS - PREBUILDING.

                    I_de.emplace_back(I_de_xyt);
                    I_d.emplace_back(I_d_xyt);
                    I_n.emplace_back(I_n_xyt);

                } else {
                    
                    // Empty small subvectors.
                    I_d.emplace_back(Vector<Real>{1});
                    I_de.emplace_back(Vector<Real>{1});
                    I_n.emplace_back(Vector<Real>{1});
                }
            }

            // FACE INTEGRALS - BUILDING.

            for(Natural k = 0; k < neighbours; ++k) {
                if(facing[k][0] != -1)
                    continue;

                // Building.
                I(dofs_j, I(dofs_j) + I_de[k] + I_d[k] + I_n[k]);
            }

            // TIME FACE INTEGRALS.

            if(bottom == -1) {

                // TIME FACE INTEGRALS - PRECOMPUTING.

                // Subvectors.
                Vector<Real> E_xy{dofs_xy};
                Vector<Real> E_t{dofs_t};

                // Face time basis.
                auto [f_phi_t, f_gradt_phi_t] = basis_t(mesh, j, Vector<Real>{1, -1.0L}); // [?]

                // TIME FACE INTEGRALS - COMPUTING.

                for(Natural k = 0; k < neighbours; ++k) { // Sub-triangulation.

                    // Nodes and basis.
                    auto [nodes2xy_j, dxy_j] = internal::reference_to_element(mesh, j, k, {nodes2x, nodes2y});
                    auto [phi_s, gradx_phi_s, grady_phi_s] = basis_s(mesh, j, nodes2xy_j);

                    Vector<Real> weights2_j = weights2 * dxy_j;

                    // Condition.
                    Vector<Real> condition = initial(nodes2xy_j[0], nodes2xy_j[1]);

                    E_xy += internal::c_scale(weights2_j, phi_s).transpose() * condition;
                }

                E_t = f_phi_t.row(0);

                // TIME FACE INTEGRALS - BUILDING.

                E(dofs_j, E(dofs_j) + kronecker(E_t, E_xy));
            }
        }

        #ifndef NVERBOSE
        std::cout << "\t[Forcing] Exited" << std::endl;
        #endif

        // Building and return
        return E + V + I;
    }

}