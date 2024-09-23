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

            // Subvectors.
            Vector<Real> V_xy{dofs_xy};
            Vector<Real> V_t{dofs_t};

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

                Vector<Real> weights2_j = weights2 * dxy_j;

                // Data.
                Vector<Real> source = data.source(nodes2xy_j[0], nodes2xy_j[1]);

                V_xy += internal::c_scale(weights2_j, phi_s).transpose() * source;
            }

            V_t = phi_t.transpose() * weights1_j;

            // VOLUME INTEGRALS - BUILDING.

            V(dofs_j, V(dofs_j) + kronecker(V_t, V_xy));

            // FACE INTEGRALS - PRECOMPUTING.

            // Subvectors.
            std::vector<Vector<Real>> I_d;
            std::vector<Vector<Real>> I_de;
            std::vector<Vector<Real>> I_n;
            
            for(Natural k = 0; k < neighbours; ++k) {

                // Nodes and basis. Using time nodes as 1D nodes.
                auto [e_nodes2xy_j, normal, e_dxy_j] = internal::reference_to_element(mesh, j, k, nodes1t);
                auto [e_phi_s, e_gradx_phi_s, e_grady_phi_s] = basis_s(mesh, j, e_nodes2xy_j);

                Vector<Real> e_weights2_j = weights1 * e_dxy_j;

                // Sign check.
                Vector<Real> negative{nodes1t.size()};
                Vector<Real> positive{nodes1t.size()};

                for(Natural j = 0; j < nodes1t.size(); ++j) {
                    negative(j, (normal(0) * convection_x(j) + normal(1) * convection_y(j) < 0.0L) ? 1.0L : 0.0L);
                    positive(j, (normal(0) * convection_x(j) + normal(1) * convection_y(j) >= 0.0L) ? 1.0L : 0.0L);
                }

                // Data.
                Vector<Real> dirichlet = data.dirichlet(e_nodes2xy_j[0], e_nodes2xy_j[1]);
                Vector<Real> neumann = data.neumann(e_nodes2xy_j[0], e_nodes2xy_j[1]);

                if(facing[k][0] == -1) {

                    // subvectors.
                    Vector<Real> I_d_xy{dofs_xy};
                    Vector<Real> I_d_t{dofs_t};

                    Vector<Real> I_de_xy{dofs_xy};
                    Vector<Real> I_de_t{dofs_t};

                    Vector<Real> I_n_xy{dofs_xy};
                    Vector<Real> I_n_t{dofs_t};

                    // FACE INTEGRALS - COMPUTING.

                    // Dirichlet, diffusion.

                    I_de_xy = internal::c_scale(e_weights2_j / e_dxy_j, e_phi_s).transpose() * dirichlet;
                    I_de_xy += internal::c_scale(e_weights2_j, (normal(0) * e_gradx_phi_s + normal(1) * e_grady_phi_s)).transpose() * dirichlet;
                    I_de_t = internal::c_scale(negative * diffusion, phi_t).transpose() * weights1_j;

                    // Dirichlet.

                    I_d_xy = -internal::c_scale(e_weights2_j, e_phi_s).transpose() * dirichlet;
                    I_d_t = internal::c_scale(negative * (normal(0) * convection_x + normal(1) * convection_y), phi_t).transpose() * weights1_j;

                    // Neumann.

                    I_n_xy = internal::c_scale(e_weights2_j, e_phi_s).transpose() * neumann;
                    I_n_t = internal::c_scale(positive, phi_t).transpose() * weights1_j;

                    // FACE INTEGRALS - PREBUILDING.

                    I_d.emplace_back(kronecker(I_d_t, I_d_xy));
                    I_de.emplace_back(kronecker(I_de_t, I_de_xy));
                    I_n.emplace_back(kronecker(I_n_t, I_n_xy));

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

                I(dofs_j, I(dofs_j) + I_d[k]);
                I(dofs_j, I(dofs_j) + I_de[k]);
                I(dofs_j, I(dofs_j) + I_n[k]);
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