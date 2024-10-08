/**
 * @file Problem_Solver.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Problem/Solver.hpp implementation.
 * @date 2024-10-07
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Ivo.hpp>

namespace ivo {

    /**
     * @brief Solves Ax = b for a 2+1 problem.
     * 
     * @param mesh Mesh.
     * @param A Stiffness matrix.
     * @param b Forcing vector.
     * @param initial Initial condition.
     * @return Vector<Real> 
     */
    Vector<Real> solve(const Mesh21 &mesh, const Sparse<Real> &A, const Vector<Real> &b, const Initial &initial) {
        #ifndef NDEBUG // Integrity check.
        assert(A.rows() == b.size());
        #endif

        // Quadrature.
        auto [nodes1x, weights1x] = quadrature1x(constants::quadrature);
        auto [nodes2x, nodes2y, weights2] = quadrature2xy(constants::quadrature);

        // Solution.
        Vector<Real> x{A.columns()};

        // Time face integrals.
        Vector<Real> E{b.size()}; // Face integrals, time.

        #ifndef NVERBOSE
        std::cout << "[Ivo] Solver" << std::endl;
        std::cout << "\t[Solver] Solving the problem's linear system" << std::endl;
        #endif

        for(Natural j = 0; j < mesh.time(); ++j) {

            // Dofs.
            std::vector<Natural> dofs_j = mesh.dofs_t(j);

            // Initial condition.
            for(Natural k = 0; k < mesh.space(); ++k) {

                // Element.
                Element21 element = mesh.element(j * mesh.space() + k);

                // Neighbours.
                Neighbour21 neighbourhood = mesh.neighbour(j * mesh.space() + k);

                std::vector<std::array<Integer, 2>> facing = neighbourhood.facing();
                Natural neighbours = facing.size();

                // Dofs.
                std::vector<Natural> dofs_k = mesh.dofs(j * mesh.space() + k);
                Natural dofs_xy = (element.p() + 1) * (element.p() + 2) / 2;
                Natural dofs_t = element.q() + 1;

                // TIME FACE INTEGRALS - PRECOMPUTING.

                // Subvectors.
                Vector<Real> E_xyt{dofs_t * dofs_xy};

                // Face time basis.
                auto [f_phi_t, f_gradt_phi_t] = basis_t(mesh, j * mesh.space() + k, Vector<Real>{1, -1.0L}); // [?]

                // TIME FACE INTEGRALS - COMPUTING.

                for(Natural h = 0; h < neighbours; ++h) { // Sub-triangulation.

                    // Nodes and basis.
                    auto [nodes2xy_k, dxy_k] = internal::reference_to_element(mesh, j * mesh.space() + k, h, {nodes2x, nodes2y});
                    auto [phi_xy, gradx_phi_xy, grady_phi_xy] = basis_xy(mesh, j * mesh.space() + k, nodes2xy_k);
                    auto [nodes2x_k, nodes2y_k] = nodes2xy_k;

                    // Weights, space.
                    Vector<Real> weights2_k = weights2 * dxy_k;

                    // Condition.
                    Vector<Real> condition{phi_xy.rows()};

                    if(j == 0) { // Initial condition.

                        for(Natural kxy = 0; kxy < phi_xy.rows(); ++kxy)
                            condition(kxy, initial(nodes2x_k(kxy), nodes2y_k(kxy)));

                    } else { // Past level.

                        // Neighbour element.
                        Element21 n_element = mesh.element((j - 1) * mesh.space() + k);

                        // Neighbour basis.
                        auto [n_phi_xy, n_gradx_phi_xy, n_grady_phi_xy] = basis_xy(mesh, (j - 1) * mesh.space() + k, nodes2xy_k);
                        auto [n_f_phi_t, n_f_gradt_phi_t] = basis_t(mesh, (j - 1) * mesh.space() + k, Vector<Real>{1, 1.0L}); // [?]

                        // Dofs.
                        std::vector<Natural> n_dofs_k = mesh.dofs((j - 1) * mesh.space() + k);
                        Natural n_dofs_xy = (n_element.p() + 1) * (n_element.p() + 2) / 2;
                        Natural n_dofs_t = n_element.q() + 1;

                        // Solution.
                        Vector<Real> uh = x(n_dofs_k);

                        for(Natural kxy = 0; kxy < phi_xy.rows(); ++kxy) {
                            Real uh_xyt = 0.0L;

                            for(Natural jt = 0; jt < n_dofs_t; ++jt)
                                for(Natural jxy = 0; jxy < n_dofs_xy; ++jxy)
                                    uh_xyt += n_f_phi_t(0, jt) * n_phi_xy(kxy, jxy) * uh(jt * n_dofs_xy + jxy);

                            condition(kxy, uh_xyt);
                        }
                    }

                    // CURRENT vs. CONDITION.

                    for(Natural jt = 0; jt < dofs_t; ++jt)
                        for(Natural jxy = 0; jxy < dofs_xy; ++jxy) {
                            Real cc_xyt = 0.0L;

                            for(Natural kxy = 0; kxy < phi_xy.rows(); ++kxy) // Brute-force integral, (*, *).
                                cc_xyt += weights2_k(kxy) * f_phi_t(0, jt) * phi_xy(kxy, jxy) * condition(kxy);

                            E_xyt(jt * dofs_xy + jxy, E_xyt(jt * dofs_xy + jxy) + cc_xyt);
                        }
                }

                // TIME FACE INTEGRALS - BUILDING.

                E(dofs_k, E(dofs_k) + E_xyt);
            }

            // Sub-matrix and sub-vector.
            Sparse<Real> A_j = Sparse<Real>{A, dofs_j, dofs_j};
            Vector<Real> b_j = b(dofs_j) + E(dofs_j);

            // Solution update.
            x(dofs_j, internal::gmres(A_j, b_j));

            #ifndef NVERBOSE
            std::cout << "\t[Solver] Solved level " << j + 1 << std::endl;
            #endif
        }

        #ifndef NVERBOSE
        std::cout << "\t[Solver] Exited" << std::endl;
        #endif

        return x;
    }

}