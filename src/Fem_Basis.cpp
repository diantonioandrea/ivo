/**
 * @file Fem_Basis.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-09-04
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Ivo.hpp>

namespace ivo {

    namespace internal {

        // Maps.

        /**
         * @brief Reference segment to element mapping.
         * 
         * @param mesh Mesh.
         * @param j Element's index.
         * @param nodes Nodes.
         * @return std::tuple<std::array<Vector<Real>, 2>, Real> 
         */
        std::tuple<Vector<Real>, Real> reference_to_element(const Mesh21 &mesh, const Natural &j, const Vector<Real> &nodes) {

            // Element.
            Element21 element = mesh.element(j);

            // Bases.
            Polygon21 bottom = element.b_base();
            Polygon21 top = element.t_base();

            // Time.
            std::array<Real, 2> t{bottom(0)(2), top(0)(2)};

            // dt.
            Real dt = (t[1] - t[0]) / 2.0L;

            // Nodes and dt.
            return {dt * nodes + (t[0] + t[1]) / 2.0L, dt};
        }

        /**
         * @brief Reference triangle to element mapping.
         * 
         * @param mesh Mesh.
         * @param j Element's index.
         * @param nodes Nodes.
         * @return std::tuple<std::array<Vector<Real>, 2>, Real> 
         */
        std::tuple<std::array<Vector<Real>, 2>, Real> reference_to_element(const Mesh21 &mesh, const Natural &j, const std::array<Vector<Real>, 2> &nodes) {

            // Nodes.
            auto [nodesx, nodesy] = nodes;

            #ifndef NDEBUG // Integrity check.
            assert(nodesx.size() == nodesy.size());
            #endif

            // Element.
            Element21 element = mesh.element(j);

            // Base.
            Polygon21 base = element.b_base();

            // Jacobian.
            Matrix<Real> J{2, 2};

            J(0, 0, base(1)(0) - base(0)(0));
            J(0, 1, base(2)(0) - base(0)(0));
            J(1, 0, base(1)(1) - base(0)(1));
            J(1, 1, base(2)(1) - base(0)(1));

            // Translation.
            Vector<Real> T{2};

            T(0, base(0)(0));
            T(1, base(0)(1));

            // Jacobian's determinant.
            Real dxy = J(0, 0) * J(1, 1) - J(0, 1) * J(1, 0);

            // Space.
            Vector<Real> x{nodesx.size()};
            Vector<Real> y{nodesy.size()};

            for(Natural k = 0; k < nodesx.size(); ++k) {
                Vector<Real> xy = J * Vector<Real>{{nodesx(k), nodesy(k)}} + T;

                x(k, xy(0));
                y(k, xy(1));
            }

            // Nodes and dxy.
            return {{x, y}, dxy};
        }

        /**
         * @brief Reference segment to element's edge mapping.
         * 
         * @param mesh Mesh.
         * @param j Element's index.
         * @param k Edge's index.
         * @param nodes Nodes.
         * @return std::tuple<std::array<Vector<Real>, 2>, std::array<Real, 2>, Real> 
         */
        std::tuple<std::array<Vector<Real>, 2>, Vector<Real>, Real> reference_to_element(const Mesh21 &mesh, const Natural &j, const Natural &k, const Vector<Real> &nodes) {

            // Element.
            Element21 element = mesh.element(j);

            // Base.
            Polygon21 base = element.b_base();

            // Edges.
            std::vector<Edge21> edges = base.edges();

            #ifndef NDEBUG // Integrity check.
            assert(k < edges.size());
            #endif

            // Edge.
            Edge21 edge = edges[k];

            // Jacobian.
            Matrix<Real> J{2, 2};

            J(0, 0, edge(1)(0) - edge(0)(0));
            J(0, 1, 0.5L * J(0, 0));
            J(1, 0, edge(1)(1) - edge(0)(1));
            J(1, 1, 0.5L * J(1, 0));

            // Translation.
            Vector<Real> T{2};

            T(0, edge(0)(0));
            T(1, edge(0)(1));

            // Space.
            Vector<Real> x{nodes.size()};
            Vector<Real> y{nodes.size()};

            for(Natural k = 0; k < nodes.size(); ++k) {
                Vector<Real> xy = J * Vector<Real>{{nodes(k), 0.0L}} + T;

                x(k, xy(0));
                y(k, xy(1));
            }

            // Edge size.
            Real de = distance(edge(0), edge(1));

            // Edge normal.
            Vector<Real> normal{2};

            normal(0, edge(1)(1) - edge(0)(1));
            normal(1, edge(0)(0) - edge(1)(0));

            normal /= norm(normal);

            return {{x, y}, normal, de};
        }

    }

    /**
     * @brief Time basis evaluation.
     * 
     * @param mesh Mesh.
     * @param j Element's index.
     * @param nodes Nodes.
     * @return std::tuple<Matrix<Real>, Matrix<Real>, Real> 
     */
    std::array<Matrix<Real>, 2> basis_t(const Mesh21 &mesh, const Natural &j, const Vector<Real> &nodes) {

        // Element.
        Element21 element = mesh.element(j);

        // Time degree.
        const Natural q = element.q();

        // Evaluations.
        const Natural rows = nodes.size();
        const Natural columns = q + 1;

        Matrix<Real> phi{rows, columns};
        Matrix<Real> gradt_phi{rows, columns};

        // Polynomial evaluations.
        for(Natural k = 0; k < columns; ++k) {
            phi.column(k, legendre1(nodes, k));
            gradt_phi.column(k, legendre_grad1(nodes, k));
        }

        return {phi, gradt_phi};
    }

    /**
     * @brief Space basis evaluation.
     * 
     * @param mesh Mesh.
     * @param j Element's index.
     * @param nodes Nodes.
     * @return std::tuple<Matrix<Real>, Matrix<Real>, Matrix<Real>, Real> 
     */
    std::array<Matrix<Real>, 3> basis_s(const Mesh21 &mesh, const Natural &j, const std::array<Vector<Real>, 2> &nodes) {
        
        // Nodes.
        auto [nodesx, nodesy] = nodes;

        #ifndef NDEBUG // Integrity check.
        assert(nodesx.size() == nodesy.size());
        #endif

        // Element.
        Element21 element = mesh.element(j);

        // Space degree.
        const Natural p = element.p();

        // Evaluations.
        const Natural rows = nodesx.size();
        const Natural columns = (p + 1) * (p + 2) / 2;

        Matrix<Real> phi{rows, columns};
        Matrix<Real> gradx_phi{rows, columns};
        Matrix<Real> grady_phi{rows, columns};

        // Base.
        Polygon21 base = element.b_base();

        #ifndef NDEBUG // Integrity check.
        assert(base.edges().size() == 3);
        #endif

        // Box.
        auto [xy_min, xy_max] = box2(base);

        Real x_min = xy_min(0);
        Real y_min = xy_min(1);
        Real x_max = xy_max(0);
        Real y_max = xy_max(1);

        // Box map.
        Matrix<Real> M{2, 2};

        M(0, 0, (x_max - x_min) / 2.0L);
        M(1, 1, (x_max - x_min) / 2.0L);

        Real M_det = M(0, 0) * M(1, 1);

        Vector<Real> T{2};

        T(0, (x_max + x_min) / 2.0L);
        T(0, (y_max + y_min) / 2.0L);

        // Inverse map.
        Matrix<Real> M_inv{2, 2};

        M_inv(0, 0, M(1, 1) / M_det);
        M_inv(1, 1, M(0, 0) / M_det);

        Vector<Real> T_inv = -(M_inv * T);

        // Space.
        Vector<Real> x{nodesx.size()};
        Vector<Real> y{nodesy.size()};

        for(Natural k = 0; k < nodesx.size(); ++k) {
            Vector<Real> xy = M_inv * Vector<Real>{{nodesx(k), nodesy(k)}} + T_inv;

            x(k, xy(0));
            y(k, xy(1));
        }

        // Degrees.
        std::vector<Natural> px, py;

        for(Natural kx = 0; kx < p + 1; ++kx)
            for(Natural ky = 0; ky < p + 1 - kx; ++ky) {
                px.emplace_back(kx);
                py.emplace_back(ky);
            }

        // Polynomial evaluations.
        for(Natural k = 0; k < columns; ++k) {
            Vector<Real> legendre_x = legendre1(x, px[k]);
            Vector<Real> legendre_y = legendre1(y, py[k]);

            Real coefficient = std::sqrt((2.0L * px[k] + 1.0L) * (2.0L * py[k] + 1.0L)) / 2.0L;

            Vector<Real> grad_legendre_x = legendre1(x, px[k]);
            Vector<Real> grad_legendre_y = legendre1(y, py[k]);

            phi.column(k, coefficient * legendre_x * legendre_y);
            gradx_phi.column(k, coefficient * grad_legendre_x * legendre_y);
            grady_phi.column(k, coefficient * legendre_x * grad_legendre_y);
        }

        // Gradients.
        for(Natural k = 0; k < rows; ++k)
            for(Natural h = 0; h < columns; ++h) {
                Vector<Real> gradient = M_inv * Vector<Real>{{gradx_phi(k, h), grady_phi(k, h)}};

                gradx_phi(k, h, gradient(0));
                grady_phi(k, h, gradient(1));
            }

        return {phi, gradx_phi, grady_phi};
    }

}