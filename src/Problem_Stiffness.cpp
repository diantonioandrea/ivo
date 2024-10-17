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
        auto [nodes1t, weights1t] = quadrature1t(constants::quadrature);
        auto [nodes1x, weights1x] = quadrature1x(constants::quadrature);
        auto [nodes2x, nodes2y, weights2] = quadrature2xy(constants::quadrature);

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

            // Time interval.
            std::array<Real, 2> interval = element.interval();

            // Dofs.
            std::vector<Natural> dofs_j = mesh.dofs(j);
            Natural dofs_xy = (element.p() + 1) * (element.p() + 2) / 2;
            Natural dofs_t = element.q() + 1;
            Natural dofs_xyt = dofs_t * dofs_xy;

            // Neighbours.
            Neighbour21 neighbourhood = mesh.neighbour(j);

            std::vector<std::array<Integer, 2>> facing = neighbourhood.facing();
            Natural neighbours = facing.size();

            // VOLUME INTEGRALS - PRECOMPUTING.

            // Submatrices.
            Matrix<Real> V_T_xyt{dofs_xyt, dofs_xyt};
            Matrix<Real> V_a_xyt{dofs_xyt, dofs_xyt};
            Matrix<Real> V_b_xyt{dofs_xyt, dofs_xyt};
            Matrix<Real> V_c_xyt{dofs_xyt, dofs_xyt};

            // Nodes and basis, time.
            auto [nodes1t_j, dt_j] = internal::reference_to_element(mesh, j, nodes1t);
            auto [phi_t, gradt_phi_t] = basis_t(mesh, j, nodes1t_j);

            // Weights, time.
            Vector<Real> weights1t_j = weights1t * dt_j;

            // VOLUME INTEGRALS - COMPUTING.

            for(Natural k = 0; k < neighbours; ++k) { // Sub-triangulation.

                // Nodes and basis, space.
                auto [nodes2xy_j, dxy_j] = internal::reference_to_element(mesh, j, k, {nodes2x, nodes2y});
                auto [phi_xy, gradx_phi_xy, grady_phi_xy] = basis_xy(mesh, j, nodes2xy_j);

                // Weights, space.
                Vector<Real> weights2_j = weights2 * dxy_j;

                // CURRENT vs. CURRENT.

                for(Natural jt = 0; jt < dofs_t; ++jt)
                    for(Natural ht = 0; ht < dofs_t; ++ht)
                        for(Natural jxy = 0; jxy < dofs_xy; ++jxy)
                            for(Natural hxy = 0; hxy < dofs_xy; ++hxy) {
                                Real V_T_cc = 0.0L;
                                Real V_a_cc = 0.0L;
                                Real V_b_cc = 0.0L;
                                Real V_c_cc = 0.0L;

                                for(Natural kt = 0; kt < phi_t.rows(); ++kt)
                                    for(Natural kxy = 0; kxy < phi_xy.rows(); ++kxy) { // Brute-force integral.
                                        Real t = nodes1t_j(kt);

                                        // Equation coefficients.
                                        auto [convection_x, convection_y] = equation.convection(t);
                                        Real diffusion = equation.diffusion(t);
                                        Real reaction = equation.reaction(t);

                                        // (*', *).

                                        V_T_cc += weights2_j(kxy) * weights1t_j(kt) * phi_t(kt, jt) * phi_xy(kxy, jxy) * gradt_phi_t(kt, ht) * phi_xy(kxy, hxy); // [!?]

                                        // a(*, *), diffusion.

                                        V_a_cc += weights2_j(kxy) * weights1t_j(kt) * (phi_t(kt, jt) * gradx_phi_xy(kxy, jxy) * phi_t(kt, ht) * gradx_phi_xy(kxy, hxy) + phi_t(kt, jt) * grady_phi_xy(kxy, jxy) * phi_t(kt, ht) * grady_phi_xy(kxy, hxy)) * diffusion;

                                        // b(*, *), convection.

                                        V_b_cc += weights2_j(kxy) * weights1t_j(kt) * phi_t(kt, ht) * (gradx_phi_xy(kxy, hxy) * convection_x + grady_phi_xy(kxy, hxy) * convection_y) * phi_t(kt, jt) * phi_xy(kxy, jxy);

                                        // c(*, *), reaction.

                                        V_c_cc += weights2_j(kxy) * weights1t_j(kt) * phi_t(kt, jt) * phi_xy(kxy, jxy) * phi_t(kt, ht) * phi_xy(kxy, hxy) * reaction;
                                    }

                                V_T_xyt(jt * dofs_xy + jxy, ht * dofs_xy + hxy, V_T_xyt(jt * dofs_xy + jxy, ht * dofs_xy + hxy) + V_T_cc);
                                V_a_xyt(jt * dofs_xy + jxy, ht * dofs_xy + hxy, V_a_xyt(jt * dofs_xy + jxy, ht * dofs_xy + hxy) + V_a_cc);
                                V_b_xyt(jt * dofs_xy + jxy, ht * dofs_xy + hxy, V_b_xyt(jt * dofs_xy + jxy, ht * dofs_xy + hxy) + V_b_cc);
                                V_c_xyt(jt * dofs_xy + jxy, ht * dofs_xy + hxy, V_c_xyt(jt * dofs_xy + jxy, ht * dofs_xy + hxy) + V_c_cc);
                            }
            }

            // VOLUME INTEGRALS - BUILDING.

            T(dofs_j, dofs_j, T(dofs_j, dofs_j) + V_T_xyt);
            V(dofs_j, dofs_j, V(dofs_j, dofs_j) + V_a_xyt);
            V(dofs_j, dofs_j, V(dofs_j, dofs_j) + V_b_xyt);
            V(dofs_j, dofs_j, V(dofs_j, dofs_j) + V_c_xyt);

            // FACE INTEGRALS - PRECOMPUTING.

            // Submatrices, c: current, n: neighbour.
            std::vector<Matrix<Real>> I_cc;
            std::vector<Matrix<Real>> I_cn;
            std::vector<Matrix<Real>> I_nc;
            std::vector<Matrix<Real>> I_nn;
            
            for(Natural k = 0; k < neighbours; ++k) {

                // Nodes and basis, space.
                auto [e_nodes2xy_j, normal, e_dxy_j] = internal::reference_to_element(mesh, j, k, nodes1x);
                auto [e_phi_xy, e_gradx_phi_xy, e_grady_phi_xy] = basis_xy(mesh, j, e_nodes2xy_j);

                // Normal gradient.
                Matrix<Real> e_gradn_phi_xy = normal(0) * e_gradx_phi_xy + normal(1) * e_grady_phi_xy;

                // Weights, space.
                Vector<Real> e_weights2_j = weights1x * e_dxy_j;

                // Submatrices.
                Matrix<Real> I_a_cc_xyt{dofs_xyt, dofs_xyt};
                Matrix<Real> I_b_cc_xyt{dofs_xyt, dofs_xyt};
                Matrix<Real> I_J_cc_xyt{dofs_xyt, dofs_xyt};

                if(facing[k][0] != -1) {

                    // Neighbour index.
                    Natural i = facing[k][0];

                    // Neighbour element.
                    Element21 n_element = mesh.element(i);

                    // Neighbour basis.
                    auto [n_e_phi_xy, n_e_gradx_phi_xy, n_e_grady_phi_xy] = basis_xy(mesh, i, e_nodes2xy_j);
                    auto [n_phi_t, n_gradt_phi_t] = basis_t(mesh, i, nodes1t_j);

                    // Normal gradient.
                    Matrix<Real> n_e_gradn_phi_xy = normal(0) * n_e_gradx_phi_xy + normal(1) * n_e_grady_phi_xy;

                    // Dofs.
                    Natural n_dofs_xy = (n_element.p() + 1) * (n_element.p() + 2) / 2;
                    Natural n_dofs_t = n_element.q() + 1;
                    Natural n_dofs_xyt = n_dofs_t * n_dofs_xy;

                    // Submatrices.
                    Matrix<Real> I_a_cn_xyt{dofs_xyt, n_dofs_xyt};
                    Matrix<Real> I_b_cn_xyt{dofs_xyt, n_dofs_xyt};
                    Matrix<Real> I_J_cn_xyt{dofs_xyt, n_dofs_xyt};

                    Matrix<Real> I_a_nc_xyt{n_dofs_xyt, dofs_xyt};
                    Matrix<Real> I_b_nc_xyt{n_dofs_xyt, dofs_xyt};
                    Matrix<Real> I_J_nc_xyt{n_dofs_xyt, dofs_xyt};

                    Matrix<Real> I_a_nn_xyt{n_dofs_xyt, n_dofs_xyt};
                    Matrix<Real> I_b_nn_xyt{n_dofs_xyt, n_dofs_xyt};
                    Matrix<Real> I_J_nn_xyt{n_dofs_xyt, n_dofs_xyt};

                    // FACE INTEGRALS - COMPUTING.

                    // CURRENT vs. CURRENT.

                    for(Natural jt = 0; jt < dofs_t; ++jt)
                        for(Natural ht = 0; ht < dofs_t; ++ht)
                            for(Natural jxy = 0; jxy < dofs_xy; ++jxy)
                                for(Natural hxy = 0; hxy < dofs_xy; ++hxy) {
                                    Real a_cc_xyt = 0.0L;
                                    Real b_cc_xyt = 0.0L;
                                    Real J_cc_xyt = 0.0L;

                                    for(Natural kt = 0; kt < phi_t.rows(); ++kt)
                                        for(Natural kxy = 0; kxy < e_phi_xy.rows(); ++kxy) { // Brute-force integral.
                                            Real t = nodes1t_j(kt);

                                            // Equation coefficients.
                                            auto [convection_x, convection_y] = equation.convection(t);
                                            Real convection_n = normal(0) * convection_x + normal(1) * convection_y;
                                            Real diffusion = equation.diffusion(t);

                                            // Boundary check.
                                            Real negative = (convection_n < 0.0L) ? 1.0L : 0.0L;

                                            // a(*, *), diffusion.

                                            if(i < j) // [?]
                                                a_cc_xyt -= weights1t_j(kt) * e_weights2_j(kxy) * (0.5L * phi_t(kt, jt) * e_gradn_phi_xy(kxy, jxy) * phi_t(kt, ht) * e_phi_xy(kxy, hxy) - 0.5L * phi_t(kt, ht) * e_gradn_phi_xy(kxy, hxy) * phi_t(kt, jt) * e_phi_xy(kxy, jxy)) * diffusion;

                                            // b(*, *), convection.

                                            b_cc_xyt -= negative * weights1t_j(kt) * e_weights2_j(kxy) * phi_t(kt, jt) * e_phi_xy(kxy, jxy) * phi_t(kt, ht) * e_phi_xy(kxy, hxy) * convection_n;

                                            // J(*, *).

                                            J_cc_xyt += weights1t_j(kt) * e_weights2_j(kxy) / e_dxy_j * phi_t(kt, jt) * e_phi_xy(kxy, jxy) * phi_t(kt, ht) * e_phi_xy(kxy, hxy) * diffusion;
                                        }

                                    I_a_cc_xyt(jt * dofs_xy + jxy, ht * dofs_xy + hxy, I_a_cc_xyt(jt * dofs_xy + jxy, ht * dofs_xy + hxy) + a_cc_xyt);
                                    I_b_cc_xyt(jt * dofs_xy + jxy, ht * dofs_xy + hxy, I_b_cc_xyt(jt * dofs_xy + jxy, ht * dofs_xy + hxy) + b_cc_xyt);
                                    I_J_cc_xyt(jt * dofs_xy + jxy, ht * dofs_xy + hxy, I_J_cc_xyt(jt * dofs_xy + jxy, ht * dofs_xy + hxy) + J_cc_xyt);
                                }

                    // CURRENT vs. NEIGHBOUR.

                    for(Natural jt = 0; jt < dofs_t; ++jt)
                        for(Natural ht = 0; ht < n_dofs_t; ++ht)
                            for(Natural jxy = 0; jxy < dofs_xy; ++jxy)
                                for(Natural hxy = 0; hxy < n_dofs_xy; ++hxy) {
                                    Real a_cn_xyt = 0.0L;
                                    Real b_cn_xyt = 0.0L;
                                    Real J_cn_xyt = 0.0L;

                                    for(Natural kt = 0; kt < phi_t.rows(); ++kt)
                                        for(Natural kxy = 0; kxy < e_phi_xy.rows(); ++kxy) { // Brute-force integral.
                                            Real t = nodes1t_j(kt);

                                            // Equation coefficients.
                                            auto [convection_x, convection_y] = equation.convection(t);
                                            Real convection_n = normal(0) * convection_x + normal(1) * convection_y;
                                            Real diffusion = equation.diffusion(t);

                                            // Boundary check.
                                            Real negative = (convection_n < 0.0L) ? 1.0L : 0.0L;

                                            // a(*, *), diffusion.

                                            if(i < j) // [?]
                                                a_cn_xyt -= weights1t_j(kt) * e_weights2_j(kxy) * (0.5L * phi_t(kt, jt) * e_gradn_phi_xy(kxy, jxy) * (-n_phi_t(kt, ht) * n_e_phi_xy(kxy, hxy)) - 0.5L * n_phi_t(kt, ht) * n_e_gradn_phi_xy(kxy, hxy) * phi_t(kt, jt) * e_phi_xy(kxy, jxy)) * diffusion;

                                            // b(*, *), convection.

                                            b_cn_xyt += negative * weights1t_j(kt) * e_weights2_j(kxy) * phi_t(kt, jt) * e_phi_xy(kxy, jxy) * n_phi_t(kt, ht) * n_e_phi_xy(kxy, hxy) * convection_n;

                                            // J(*, *).

                                            J_cn_xyt -= weights1t_j(kt) * e_weights2_j(kxy) / e_dxy_j * phi_t(kt, jt) * e_phi_xy(kxy, jxy) * n_phi_t(kt, ht) * n_e_phi_xy(kxy, hxy) * diffusion;
                                        }

                                    I_a_cn_xyt(jt * dofs_xy + jxy, ht * n_dofs_xy + hxy, I_a_cn_xyt(jt * dofs_xy + jxy, ht * n_dofs_xy + hxy) + a_cn_xyt);
                                    I_b_cn_xyt(jt * dofs_xy + jxy, ht * n_dofs_xy + hxy, I_b_cn_xyt(jt * dofs_xy + jxy, ht * n_dofs_xy + hxy) + b_cn_xyt);
                                    I_J_cn_xyt(jt * dofs_xy + jxy, ht * n_dofs_xy + hxy, I_J_cn_xyt(jt * dofs_xy + jxy, ht * n_dofs_xy + hxy) + J_cn_xyt);
                                }

                    // NEIGHBOUR vs. CURRENT.

                    for(Natural jt = 0; jt < n_dofs_t; ++jt)
                        for(Natural ht = 0; ht < dofs_t; ++ht)
                            for(Natural jxy = 0; jxy < n_dofs_xy; ++jxy)
                                for(Natural hxy = 0; hxy < dofs_xy; ++hxy) {
                                    Real a_nc_xyt = 0.0L;
                                    Real b_nc_xyt = 0.0L;
                                    Real J_nc_xyt = 0.0L;

                                    for(Natural kt = 0; kt < phi_t.rows(); ++kt)
                                        for(Natural kxy = 0; kxy < e_phi_xy.rows(); ++kxy) { // Brute-force integral.
                                            Real t = nodes1t_j(kt);

                                            // Equation coefficients.
                                            Real diffusion = equation.diffusion(t);

                                            // a(*, *), diffusion.

                                            if(i < j) // [?]
                                                a_nc_xyt -= weights1t_j(kt) * e_weights2_j(kxy) * (0.5L * n_phi_t(kt, jt) * n_e_gradn_phi_xy(kxy, jxy) * phi_t(kt, ht) * e_phi_xy(kxy, hxy) - 0.5L * phi_t(kt, ht) * e_gradn_phi_xy(kxy, hxy) * (-n_phi_t(kt, jt) * n_e_phi_xy(kxy, jxy))) * diffusion;

                                            // J(*, *).

                                            J_nc_xyt -= weights1t_j(kt) * e_weights2_j(kxy) / e_dxy_j * n_phi_t(kt, jt) * n_e_phi_xy(kxy, jxy) * phi_t(kt, ht) * e_phi_xy(kxy, hxy) * diffusion;
                                        }

                                    I_a_nc_xyt(jt * n_dofs_xy + jxy, ht * dofs_xy + hxy, I_a_nc_xyt(jt * n_dofs_xy + jxy, ht * dofs_xy + hxy) + a_nc_xyt);
                                    I_b_nc_xyt(jt * n_dofs_xy + jxy, ht * dofs_xy + hxy, I_b_nc_xyt(jt * n_dofs_xy + jxy, ht * dofs_xy + hxy) + b_nc_xyt);
                                    I_J_nc_xyt(jt * n_dofs_xy + jxy, ht * dofs_xy + hxy, I_J_nc_xyt(jt * n_dofs_xy + jxy, ht * dofs_xy + hxy) + J_nc_xyt);
                                }

                    // NEIGHBOUR vs. NEIGHBOUR.

                    for(Natural jt = 0; jt < n_dofs_t; ++jt)
                        for(Natural ht = 0; ht < n_dofs_t; ++ht)
                            for(Natural jxy = 0; jxy < n_dofs_xy; ++jxy)
                                for(Natural hxy = 0; hxy < n_dofs_xy; ++hxy) {
                                    Real a_nn_xyt = 0.0L;
                                    Real b_nn_xyt = 0.0L;
                                    Real J_nn_xyt = 0.0L;

                                    for(Natural kt = 0; kt < phi_t.rows(); ++kt)
                                        for(Natural kxy = 0; kxy < e_phi_xy.rows(); ++kxy) { // Brute-force integral.
                                            Real t = nodes1t_j(kt);

                                            // Equation coefficients.
                                            Real diffusion = equation.diffusion(t);

                                            // a(*, *), diffusion.

                                            if(i < j) // [?]
                                                a_nn_xyt -= weights1t_j(kt) * e_weights2_j(kxy) * (0.5L * n_phi_t(kt, jt) * n_e_gradn_phi_xy(kxy, jxy) * (-n_phi_t(kt, ht) * n_e_phi_xy(kxy, hxy)) - 0.5L * n_phi_t(kt, ht) * n_e_gradn_phi_xy(kxy, hxy) * (-n_phi_t(kt, jt) * n_e_phi_xy(kxy, jxy))) * diffusion;

                                            // J(*, *).

                                            J_nn_xyt += weights1t_j(kt) * e_weights2_j(kxy) / e_dxy_j * n_phi_t(kt, jt) * n_e_phi_xy(kxy, jxy) * n_phi_t(kt, ht) * n_e_phi_xy(kxy, hxy) * diffusion;
                                        }

                                    I_a_nn_xyt(jt * n_dofs_xy + jxy, ht * n_dofs_xy + hxy, I_a_nn_xyt(jt * n_dofs_xy + jxy, ht * n_dofs_xy + hxy) + a_nn_xyt);
                                    I_b_nn_xyt(jt * n_dofs_xy + jxy, ht * n_dofs_xy + hxy, I_b_nn_xyt(jt * n_dofs_xy + jxy, ht * n_dofs_xy + hxy) + b_nn_xyt);
                                    I_J_nn_xyt(jt * n_dofs_xy + jxy, ht * n_dofs_xy + hxy, I_J_nn_xyt(jt * n_dofs_xy + jxy, ht * n_dofs_xy + hxy) + J_nn_xyt);
                                }

                    // FACE INTEGRALS - PREBUILDING.

                    I_cc.emplace_back(I_a_cc_xyt + I_b_cc_xyt + I_J_cc_xyt);
                    I_cn.emplace_back(I_a_cn_xyt + I_b_cn_xyt + I_J_cn_xyt);
                    I_nc.emplace_back(I_a_nc_xyt + I_b_nc_xyt + I_J_nc_xyt);
                    I_nn.emplace_back(I_a_nn_xyt + I_b_nn_xyt + I_J_nn_xyt);

                } else {

                    // FACE INTEGRALS - COMPUTING.

                    for(Natural jt = 0; jt < dofs_t; ++jt)
                        for(Natural ht = 0; ht < dofs_t; ++ht)
                            for(Natural jxy = 0; jxy < dofs_xy; ++jxy)
                                for(Natural hxy = 0; hxy < dofs_xy; ++hxy) {
                                    Real a_cc_xyt = 0.0L;
                                    Real b_cc_xyt = 0.0L;
                                    Real J_cc_xyt = 0.0L;

                                    for(Natural kt = 0; kt < phi_t.rows(); ++kt)
                                        for(Natural kxy = 0; kxy < e_phi_xy.rows(); ++kxy) { // Brute-force integral.
                                            Real t = nodes1t_j(kt);

                                            // Equation coefficients.
                                            auto [convection_x, convection_y] = equation.convection(t);
                                            Real convection_n = normal(0) * convection_x + normal(1) * convection_y;
                                            Real diffusion = equation.diffusion(t);

                                            // Boundary check.
                                            Real negative = (convection_n < 0.0L) ? 1.0L : 0.0L;

                                            // a(*, *), diffusion.

                                            a_cc_xyt -= negative * weights1t_j(kt) * e_weights2_j(kxy) * (phi_t(kt, jt) * e_gradn_phi_xy(kxy, jxy) * phi_t(kt, ht) * e_phi_xy(kxy, hxy) - phi_t(kt, ht) * e_gradn_phi_xy(kxy, hxy) * phi_t(kt, jt) * e_phi_xy(kxy, jxy)) * diffusion;

                                            // b(*, *), convection.

                                            b_cc_xyt -= negative * weights1t_j(kt) * e_weights2_j(kxy) * phi_t(kt, jt) * e_phi_xy(kxy, jxy) * phi_t(kt, ht) * e_phi_xy(kxy, hxy) * convection_n;

                                            // J(*, *).

                                            J_cc_xyt += negative * weights1t_j(kt) * e_weights2_j(kxy) / e_dxy_j * phi_t(kt, jt) * e_phi_xy(kxy, jxy) * phi_t(kt, ht) * e_phi_xy(kxy, hxy) * diffusion;
                                        }

                                    I_a_cc_xyt(jt * dofs_xy + jxy, ht * dofs_xy + hxy, I_a_cc_xyt(jt * dofs_xy + jxy, ht * dofs_xy + hxy) + a_cc_xyt);
                                    I_b_cc_xyt(jt * dofs_xy + jxy, ht * dofs_xy + hxy, I_b_cc_xyt(jt * dofs_xy + jxy, ht * dofs_xy + hxy) + b_cc_xyt);
                                    I_J_cc_xyt(jt * dofs_xy + jxy, ht * dofs_xy + hxy, I_J_cc_xyt(jt * dofs_xy + jxy, ht * dofs_xy + hxy) + J_cc_xyt);
                                }

                    // FACE INTEGRALS - PREBUILDING.

                    I_cc.emplace_back(I_a_cc_xyt + I_b_cc_xyt + I_J_cc_xyt);
                    I_cn.emplace_back(Matrix<Real>{1, 1});
                    I_nc.emplace_back(Matrix<Real>{1, 1});
                    I_nn.emplace_back(Matrix<Real>{1, 1});
                }
            }

            // FACE INTEGRALS - BUILDING.

            for(Natural k = 0; k < neighbours; ++k) {
                if(facing[k][0] != -1) {

                    // Neighbour dofs.
                    std::vector<Natural> n_dofs_j = mesh.dofs(facing[k][0]);

                    // Building.
                    I(dofs_j, dofs_j, I(dofs_j, dofs_j) + I_cc[k]);
                    I(dofs_j, n_dofs_j, I(dofs_j, n_dofs_j) + I_cn[k]);
                    I(n_dofs_j, dofs_j, I(n_dofs_j, dofs_j) + I_nc[k]);
                    I(n_dofs_j, n_dofs_j, I(n_dofs_j, n_dofs_j) + I_nn[k]); // [?]
                    
                } else {
                    
                    // Building.
                    I(dofs_j, dofs_j, I(dofs_j, dofs_j) + I_cc[k]);
                }
            }

            // TIME FACE INTEGRALS - PRECOMPUTING.

            // Face time basis.
            auto [f_phi_t, f_gradt_phi_t] = basis_t(mesh, j, Vector<Real>{1, interval[0]}); // [?]

            // Submatrix.
            Matrix<Real> E_cc_xyt{dofs_xyt, dofs_xyt};

            // TIME FACE INTEGRALS - COMPUTING.

            for(Natural k = 0; k < neighbours; ++k) { // Sub-triangulation.

                // Nodes and basis.
                auto [nodes2xy_j, dxy_j] = internal::reference_to_element(mesh, j, k, {nodes2x, nodes2y});
                auto [phi_xy, gradx_phi_xy, grady_phi_xy] = basis_xy(mesh, j, nodes2xy_j);

                // Weights, space.
                Vector<Real> weights2_j = weights2 * dxy_j;

                // CURRENT vs. CURRENT.

                for(Natural jt = 0; jt < dofs_t; ++jt)
                    for(Natural ht = 0; ht < dofs_t; ++ht)
                        for(Natural jxy = 0; jxy < dofs_xy; ++jxy)
                            for(Natural hxy = 0; hxy < dofs_xy; ++hxy) {
                                Real cc_xyt = 0.0L;

                                for(Natural kxy = 0; kxy < phi_xy.rows(); ++kxy) // Brute-force integral, (*, *).
                                    cc_xyt += weights2_j(kxy) * f_phi_t(0, jt) * phi_xy(kxy, jxy) * f_phi_t(0, ht) * phi_xy(kxy, hxy);
                            
                                E_cc_xyt(jt * dofs_xy + jxy, ht * dofs_xy + hxy, E_cc_xyt(jt * dofs_xy + jxy, ht * dofs_xy + hxy) + cc_xyt);
                            }
            }

            // TIME FACE INTEGRALS - BUILDING.

            E(dofs_j, dofs_j, E(dofs_j, dofs_j) + E_cc_xyt);

            #ifndef NVERBOSE
            if((j + 1) % mesh.space() == 0)
                std::cout << "\t[Stiffness] Progress: " << j / mesh.space() + 1 << "/" << mesh.time() << std::endl;
            #endif
        }

        #ifndef NVERBOSE
        std::cout << "\t[Stiffness] Exited" << std::endl;
        #endif

        // Building and return.

        return T + E + V + I;
    }

}