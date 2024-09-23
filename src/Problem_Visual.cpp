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

    // /**
    //  * @brief Solution visualization.
    //  * 
    //  * @param mesh Mesh.
    //  * @param solution Numerical solution.
    //  * @param filename Output file.
    //  */
    // void visual(const Mesh21 &mesh, const Vector<Real> &solution, const std::string &filename) {

    //     // Quadrature.
    //     auto [nodes1t, weights1] = quadrature1(constants::quadrature);
    //     auto [nodes2x, nodes2y, weights2] = quadrature2(constants::quadrature);

    //     // Output.
    //     std::ofstream output(filename);

    //     // Loop over elements.
    //     for(Natural j = 0; j < mesh.space() * mesh.time(); ++j) {

    //         // Element.
    //         Element21 element = mesh.element(j);

    //         // Dofs.
    //         std::vector<Natural> dofs_j = mesh.dofs(j);

    //         // Nodes and basis.
    //         auto [nodes1t_j, dt_j] = internal::reference_to_element(mesh, j, nodes1t);
    //         auto [nodes2xy_j, dxy_j] = internal::reference_to_element(mesh, j, {nodes2x, nodes2y});
    //         auto [nodes2x_j, nodes2y_j] = nodes2xy_j;

    //         auto [phi_s, gradx_phi_s, grady_phi_s] = basis_s(mesh, j, nodes2xy_j);
    //         auto [phi_t, gradt_phi_t] = basis_t(mesh, j, nodes1t);

    //         // Full basis.
    //         Matrix<Real> phi = kronecker(phi_t, phi_s);

    //         // Full nodes.
    //         std::vector<std::array<Real, 3>> nodes;

    //         for(Natural k = 0; k < nodes1t_j.size(); ++k)
    //             for(Natural h = 0; h < nodes2x_j.size(); ++h)
    //                 nodes.emplace_back(std::array<Real, 3>{nodes2x_j[h], nodes2y_j[h], nodes1t_j[k]});

    //         // Evaluation.
    //         Vector<Real> evaluation = phi * solution(dofs_j);

    //         // Output.
    //         for(Natural j = 0; j < nodes.size(); ++j) {
    //             output << nodes[j][0] << "," << nodes[j][1] << "," << nodes[j][2] << "," << evaluation(j) << std::endl;
    //         }
    //     }
    // }

}