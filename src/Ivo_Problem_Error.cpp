/**
 * @file Problem_Error.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Problem/Error.hpp implementation.
 * @date 2024-10-19
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Ivo.hpp>

namespace ivo {

    namespace internal {
        
        /**
         * @brief L2 error at given time.
         * 
         * @param mesh Mesh.
         * @param uh Numerical solution.
         * @param u Exact solution.
         * @param t Time.
         * @return Real 
         */
        Real error_at_time(const Mesh21 &mesh, const Vector<Real> &uh, const std::function<Real (Real, Real, Real)> &u, const Real &t) {
            
            // Error.
            Real l2 = 0.0L;

            // Quadrature.
            auto [nodes1x, weights1x] = quadrature1x(constants::quadrature);
            auto [nodes2x, nodes2y, weights2] = quadrature2xy(constants::quadrature);

            // Time index and ratio.
            Natural i = 0;

            for(Natural j = 0; j < mesh.space() * mesh.time(); j += mesh.space()) {

                // ELEMENT DATA.

                // Element.
                Element21 element = mesh.element(j);

                // Time interval.
                auto [t0, t1] = element.interval();

                if((t >= t0) && (t < t1)) {
                    i = j / mesh.space();
                    break;
                }
                
            }

            // Loop over elements.
            for(Natural j = mesh.space() * i; j < mesh.space() * (i + 1); ++j) {

                // ELEMENT DATA.

                // Element.
                Element21 element = mesh.element(j);

                // Dofs.
                std::vector<Natural> dofs_j = mesh.dofs(j);
                Natural dofs_xy = (element.p() + 1) * (element.p() + 2) / 2;
                Natural dofs_t = element.q() + 1;

                // Neighbours.
                Neighbour21 neighbourhood = mesh.neighbour(j);

                std::vector<std::array<Integer, 2>> facing = neighbourhood.facing();
                Natural neighbours = facing.size();

                // Nodes and basis, time.
                auto [phi_t, gradt_phi_t] = basis_t(mesh, j, Vector<Real>(1, t));

                // INTEGRALS - COMPUTING.

                for(Natural k = 0; k < neighbours; ++k) { // Sub-triangulation.

                    // Nodes and basis, space.
                    auto [nodes2xy_j, dxy_j] = internal::reference_to_element(mesh, j, k, {nodes2x, nodes2y});
                    auto [phi_xy, gradx_phi_xy, grady_phi_xy] = basis_xy(mesh, j, nodes2xy_j);
                    auto [nodes2x, nodes2y] = nodes2xy_j;

                    // Weights, space.
                    Vector<Real> weights2_j = weights2 * dxy_j;

                    // Local coefficients and solution.
                    Vector<Real> u_j = uh(dofs_j);

                    Matrix<Real> uh_j{phi_t.rows(), phi_xy.rows()};

                    for(Natural kxy = 0; kxy < phi_xy.rows(); ++kxy) {
                        Real uh_xyt = 0.0;

                        for(Natural jt = 0; jt < dofs_t; ++jt)
                            for(Natural jxy = 0; jxy < dofs_xy; ++jxy)
                                uh_xyt += phi_t(0, jt) * phi_xy(kxy, jxy) * u_j(jt * dofs_xy + jxy);

                        uh_j(0, kxy, uh_xyt);
                    }

                    // CURRENT ERROR.

                    for(Natural kxy = 0; kxy < phi_xy.rows(); ++kxy) { // Brute-force integral.
                        Real x = nodes2x[kxy];
                        Real y = nodes2y[kxy];

                        // L2 error.

                        l2 += weights2_j(kxy) * (u(x, y, t) - uh_j(0, kxy)) * (u(x, y, t) - uh_j(0, kxy));
                    }
                }
            }

            // Error.
            return std::sqrt(l2);
        }

    }

    // Constructor.
    
    /**
     * @brief Default constructor.
     * 
     * @param mesh Mesh.
     * @param equation Equation.
     * @param uh Numerical solution.
     * @param u Exact solution.
     * @param u_xy Exact solution's gradient.
     */
    Error::Error(const Mesh21 &mesh, const Equation &equation, const Vector<Real> &uh, const std::function<Real (Real, Real, Real)> &u, const std::function<std::array<Real, 2> (Real, Real, Real)> &u_xy): dofs{mesh.dofs()}, p{mesh.p()}, q{mesh.q()}, h{mesh.h()}, t{mesh.t()} {

        // Errors.
        this->l2l2s.resize(mesh.space() * mesh.time(), 0.0L);
        this->l2h1s.resize(mesh.space() * mesh.time(), 0.0L);
        this->l2Ts.resize(mesh.space(), 0.0L);
        this->l2l2 = 0.0L;
        this->l2h1 = 0.0L;
        this->l2T = 0.0L;
        this->linfl2 = 0.0L;

        // Quadrature.
        auto [nodes1t, weights1t] = quadrature1t(constants::quadrature);
        auto [nodes1x, weights1x] = quadrature1x(constants::quadrature);
        auto [nodes2x, nodes2y, weights2] = quadrature2xy(constants::quadrature);

        #ifndef NVERBOSE
        std::cout << "[Ivo] Error" << std::endl;
        std::cout << "\t[Error] Evaluating L2(L2) and L2(H1) errors." << std::endl;
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

            // Neighbours.
            Neighbour21 neighbourhood = mesh.neighbour(j);

            std::vector<std::array<Integer, 2>> facing = neighbourhood.facing();
            Natural neighbours = facing.size();

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
                auto [nodes2x, nodes2y] = nodes2xy_j;

                // Weights, space.
                Vector<Real> weights2_j = weights2 * dxy_j;

                // Local coefficients and solution.
                Vector<Real> u_j = uh(dofs_j);

                Matrix<Real> uh_j{phi_t.rows(), phi_xy.rows()};
                Matrix<Real> uh_x_j{phi_t.rows(), phi_xy.rows()};
                Matrix<Real> uh_y_j{phi_t.rows(), phi_xy.rows()};

                for(Natural kt = 0; kt < phi_t.rows(); ++kt)
                    for(Natural kxy = 0; kxy < phi_xy.rows(); ++kxy) {
                        Real uh_xyt = 0.0;
                        Real uh_x_xyt = 0.0;
                        Real uh_y_xyt = 0.0;

                        for(Natural jt = 0; jt < dofs_t; ++jt)
                            for(Natural jxy = 0; jxy < dofs_xy; ++jxy) {
                                uh_xyt += phi_t(kt, jt) * phi_xy(kxy, jxy) * u_j(jt * dofs_xy + jxy);
                                uh_x_xyt += phi_t(kt, jt) * gradx_phi_xy(kxy, jxy) * u_j(jt * dofs_xy + jxy);
                                uh_y_xyt += phi_t(kt, jt) * grady_phi_xy(kxy, jxy) * u_j(jt * dofs_xy + jxy);
                            }

                        uh_j(kt, kxy, uh_xyt);
                        uh_x_j(kt, kxy, uh_x_xyt);
                        uh_y_j(kt, kxy, uh_y_xyt);
                    }

                // CURRENT ERROR.

                for(Natural kt = 0; kt < phi_t.rows(); ++kt)
                    for(Natural kxy = 0; kxy < phi_xy.rows(); ++kxy) { // Brute-force integral.
                        Real x = nodes2x[kxy];
                        Real y = nodes2y[kxy];
                        Real t = nodes1t_j[kt];

                        // Gradient.
                        auto [u_x, u_y] = u_xy(x, y, t);

                        // L2(L2) error.

                        this->l2l2s[j] += weights2_j(kxy) * weights1t_j(kt) * (u(x, y, t) - uh_j(kt, kxy)) * (u(x, y, t) - uh_j(kt, kxy));

                        // L2(H1) error.

                        this->l2h1s[j] += weights2_j(kxy) * weights1t_j(kt) * (u_x - uh_x_j(kt, kxy)) * (u_x - uh_x_j(kt, kxy));
                        this->l2h1s[j] += weights2_j(kxy) * weights1t_j(kt) * (u_y - uh_y_j(kt, kxy)) * (u_y - uh_y_j(kt, kxy));
                    }
            }

            // Error update.
            this->l2h1s[j] *= equation.diffusion();

            this->l2l2 += this->l2l2s[j];
            this->l2h1 += this->l2h1s[j];

            #ifndef NVERBOSE
            if((j + 1) % mesh.space() == 0)
                std::cout << "\t\t[Error] Progress: " << j / mesh.space() + 1 << "/" << mesh.time() << std::endl;
            #endif
        }

        #ifndef NVERBOSE
        std::cout << "\t[Error] Evaluating L2(T) error." << std::endl;
        #endif

        // Loop over last elements.
        for(Natural j = mesh.space() * (mesh.time() - 1); j < mesh.space() * mesh.time(); ++j) {

            // Index.
            Natural i = j - mesh.space() * (mesh.time() - 1);

            // ELEMENT DATA.

            // Element.
            Element21 element = mesh.element(j);

            // Dofs.
            std::vector<Natural> dofs_j = mesh.dofs(j);
            Natural dofs_xy = (element.p() + 1) * (element.p() + 2) / 2;
            Natural dofs_t = element.q() + 1;

            // Neighbours.
            Neighbour21 neighbourhood = mesh.neighbour(j);

            std::vector<std::array<Integer, 2>> facing = neighbourhood.facing();
            Natural neighbours = facing.size();

            // Time interval.
            std::array<Real, 2> interval = element.interval();

            // Nodes and basis, time.
            auto [phi_t, gradt_phi_t] = basis_t(mesh, j, Vector<Real>(1, interval[1]));

            // INTEGRALS - COMPUTING.

            for(Natural k = 0; k < neighbours; ++k) { // Sub-triangulation.

                // Nodes and basis, space.
                auto [nodes2xy_j, dxy_j] = internal::reference_to_element(mesh, j, k, {nodes2x, nodes2y});
                auto [phi_xy, gradx_phi_xy, grady_phi_xy] = basis_xy(mesh, j, nodes2xy_j);
                auto [nodes2x, nodes2y] = nodes2xy_j;

                // Weights, space.
                Vector<Real> weights2_j = weights2 * dxy_j;

                // Local coefficients and solution.
                Vector<Real> u_j = uh(dofs_j);

                Matrix<Real> uh_j{phi_t.rows(), phi_xy.rows()};

                for(Natural kxy = 0; kxy < phi_xy.rows(); ++kxy) {
                    Real uh_xyt = 0.0;

                    for(Natural jt = 0; jt < dofs_t; ++jt)
                        for(Natural jxy = 0; jxy < dofs_xy; ++jxy)
                            uh_xyt += phi_t(0, jt) * phi_xy(kxy, jxy) * u_j(jt * dofs_xy + jxy);

                    uh_j(0, kxy, uh_xyt);
                }

                // CURRENT ERROR.

                for(Natural kxy = 0; kxy < phi_xy.rows(); ++kxy) { // Brute-force integral.
                    Real x = nodes2x[kxy];
                    Real y = nodes2y[kxy];
                    Real t = interval[1];

                    // L2(T) error.

                    this->l2Ts[i] += weights2_j(kxy) * (u(x, y, t) - uh_j(0, kxy)) * (u(x, y, t) - uh_j(0, kxy));
                }
            }

            // Error update.
            this->l2T += this->l2Ts[i];
        }

        #ifndef NVERBOSE
        std::cout << "\t[Error] Evaluating Linf(L2) error." << std::endl;
        #endif

        // Interval.
        auto [t0, t1] = mesh.element(0).interval();
        t1 = mesh.element(mesh.space() * mesh.time() - 1).interval()[1];

        // Steps.
        Natural steps = static_cast<Natural>(100.0 * (t1 - t0));

        while(std::abs(t1 - t0) > 1.0E-4L) {

            // Max.
            Real t_max = t0;

            for(Real t = t0; t <= t1; t += (t1 - t0) / steps) {
                Real error = internal::error_at_time(mesh, uh, u, t);

                if(error >= this->linfl2) {
                    this->linfl2 = error;
                    t_max = t;
                }
            }

            // Interval update.
            t0 = (t0 + t_max) / 2.0L;
            t1 = (t1 + t_max) / 2.0L;

            #ifndef NVERBOSE
            std::cout << "\t\t[Error] Progress: " << std::abs(t1 - t0) << std::endl;
            #endif
        }

        // Error.
        this->l2l2 = std::sqrt(this->l2l2);
        this->l2h1 = std::sqrt(this->l2h1);
        this->l2T = std::sqrt(this->l2T);

        #ifndef NVERBOSE
        std::cout << "\t[Error] Exited" << std::endl;
        #endif
    }

    // Output.

    /**
     * @brief Error output.
     * 
     * @param ost 
     * @param error 
     * @return std::ostream 
     */
    std::ostream &operator <<(std::ostream &ost, const Error &error) {
        ost << "[Ivo] Error" << std::endl;
        ost << "\t[Error] DoFs: " << error.dofs << std::endl;
        ost << "\t[Error] Space diagram size, h: " << error.h << std::endl;
        ost << "\t[Error] Time diagram size, t: " << error.t << std::endl;
        ost << "\t[Error] (Highest) space degree, p: " << error.p << std::endl;
        ost << "\t[Error] (Highest) time degree, q: " << error.q << std::endl;
        ost << "\t[Error] L2(L2) error, l2l2: " << error.l2l2 << std::endl;
        ost << "\t[Error] L2(T) error, l2T: " << error.l2T << std::endl;
        ost << "\t[Error] L2(H1) error, l2h1: " << error.l2h1 << std::endl;
        ost << "\t[Error] Linf(L2) error, linfl2: " << error.linfl2 << std::flush;

        return ost;
    }

}