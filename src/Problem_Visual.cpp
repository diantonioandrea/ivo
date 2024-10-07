/**
 * @file Problem_Visual.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Problem/Visual.hpp implementation.
 * @date 2024-09-18
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Ivo.hpp>

namespace ivo {

    /**
     * @brief Solution visualization.
     * 
     * @param mesh Mesh.
     * @param solution Numerical solution.
     * @param filename Output file.
     */
    void visual(const Mesh21 &mesh, const Vector<Real> &solution, const std::string &filename) {

        // Quadrature, visualization only.
        auto [nodes1t, weights1t] = quadrature1t(3);
        auto [nodes2x, nodes2y, weights2] = quadrature2xy(3);

        // Output.
        std::ofstream output(filename);

        #ifndef NVERBOSE
        std::cout << "[Ivo] Visual" << std::endl;
        std::cout << "\t[Visual] Creating the solution visualization" << std::endl;
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

            // Nodes and basis.
            auto [nodes1t_j, dt_j] = internal::reference_to_element(mesh, j, nodes1t);
            auto [phi_t, gradt_phi_t] = basis_t(mesh, j, nodes1t);

            for(Natural k = 0; k < neighbours; ++k) { // Sub-triangulation.

                // Nodes and basis.
                auto [nodes2xy_j, dxy_j] = internal::reference_to_element(mesh, j, k, {nodes2x, nodes2y});
                auto [phi_s, gradx_phi_s, grady_phi_s] = basis_xy(mesh, j, nodes2xy_j);
                auto [nodes2x_j, nodes2y_j] = nodes2xy_j;

                // Local solution.
                Vector<Real> uh_j = solution(dofs_j);

                for(Natural kt = 0; kt < phi_t.rows(); ++kt)
                    for(Natural kxy = 0; kxy < phi_s.rows(); ++kxy) { // Brute-force integral.
                        Real x = nodes2x_j(kxy);
                        Real y = nodes2y_j(kxy);
                        Real t = nodes1t_j(kt);

                        Real uh = 0.0L;

                        for(Natural jt = 0; jt < dofs_t; ++jt)
                            for(Natural jxy = 0; jxy < dofs_xy; ++jxy)
                                uh += phi_t(kt, jt) * phi_s(kxy, jxy) * uh_j(jt * dofs_xy + jxy);


                        output << x << "," << y << "," << t << "," << uh << std::endl;
                    }
            }
        }

        #ifndef NVERBOSE
        std::cout << "\t[Visual] Exited" << std::endl;
        #endif
    }

}