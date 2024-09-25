/**
 * @file Solvers.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Linear solvers.
 * @date 2024-09-18
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef ALGEBRA_METHODS_SOLVERS
#define ALGEBRA_METHODS_SOLVERS

#include "../Sparse.hpp"

namespace ivo {

    namespace internal {

        /**
         * @brief Restarted GMRES, solves Ax = b for x.
         * 
         * @tparam T Numerical type.
         * @param A Sparse matrix.
         * @param b Vector.
         * @return Vector<T> 
         */
        template<Numerical T>
        Vector<T> gmres(Sparse<T> &A, const Vector<T> &b) {
            #ifndef NDEBUG // Integrity check.
            assert(A.rows() == A.columns());
            assert(A.rows() == b.size());
            #endif

            // Solution.
            Vector<T> x{A.columns()};

            // Iterations and GMRES counter.
            Natural iterations = 0;
            Natural m = 1;

            // Residual.
            Vector<T> residual = b - A * x;

            #ifndef NVERBOSE
            std::cout << "[Ivo] Restarted GMRES" << std::endl;
            std::cout << "\t[Restarted GMRES] Solving a linear system" << std::endl;
            std::cout << "\t[Restarted GMRES] Starting, residual: " << norm(residual) << std::endl;
            #endif

            do {
                ++iterations;

                // Residual's norm.
                Real residual_norm = norm(residual);

                // Hessenberg matrix.
                Matrix<T> H{m + 1, m};

                // Krylov base.
                std::vector<Vector<T>> Vs;
                Vs.emplace_back(residual / residual_norm);

                for(Natural j = 0; j < m; ++j) {
                    Vector<T> w = A * Vs[j];

                    for(Natural k = 0; k <= j; ++k) {
                        H(k, j, dot(w, Vs[k]));
                        w -= H(k, j) * Vs[k];
                    }

                    // New element for H.
                    H(j + 1, j, norm(w));

                    // New base element.
                    Vs.emplace_back(w / H(j + 1, j));
                }

                // Base matrix.
                Matrix<T> V{A.rows(), m};

                for(Natural j = 0; j < m; ++j)
                    V.column(j, Vs[j]);

                // Solution by least squares.

                // Right-hand side.
                Vector<T> rhs{m + 1};
                rhs(0, residual_norm);

                // Identity matrix.
                Matrix<T> identity{m + 1, m + 1};

                for(Natural j = 0; j < m + 1; ++j)
                    identity(j, j, static_cast<T>(1));

                // Rotations.
                for(Natural j = 0; j < m; ++j) {

                    // Rotation matrix.
                    Matrix<T> rotation = identity;

                    // Coefficients.
                    T first = H(j, j);
                    T second = H(j + 1, j);
                    T third = std::sqrt(first * first + second * second);

                    // Rotation coefficients.
                    T c = first / third;
                    T s = second / third;

                    // Rotation building.
                    rotation(j, j, c);
                    rotation(j + 1, j + 1, c);
                    rotation(j, j + 1, s);
                    rotation(j + 1, j, -s);

                    // Rotation.
                    H = rotation * H;
                    rhs = rotation * rhs;
                }

                // Backward substitution.
                Vector<T> y{m};

                for(Natural j = m; j > 0; --j) {
                    T sum = static_cast<T>(0);

                    for(Natural k = j; k < m; ++k)
                        sum += y(k) * H(j - 1, k);
                        
                    y(j - 1, (rhs(j - 1) - sum) / H(j - 1, j - 1));
                }

                // Solution estimate.
                x += V * y;

                // Residual re-evaluation.
                residual = b - A * x;

                // Exit conditions.
                if(std::abs(rhs(m)) < constants::algebra_zero)
                    break;

                // Size update.
                m = (m > constants::gmres_restart) ? 1 : m + 1;

                #ifndef NVERBOSE // Updates.
                if(iterations % 25 == 0)
                    std::cout << "\t[Restarted GMRES] Completed iteration " << iterations << ", residual: " << norm(residual) << std::endl;
                #endif

            } while(iterations < constants::solvers_stop);

            #ifndef NVERBOSE
            std::cout << "\t[Restarted GMRES] Exited, iterations: " << iterations << ", " << " residual: " << norm(residual) << std::endl;
            #endif

            return x;
        }

    }

    // Solver wrapper.

    /**
     * @brief Solves Ax = b for x.
     * 
     * @tparam T Numerical type.
     * @param A Sparse matrix.
     * @param b Vector.
     * @return Vector<T> 
     */
    template<Numerical T>
    Vector<T> solve(Sparse<T> &A, const Vector<T> &b) {
        return internal::gmres(A, b);
    }

}

#endif