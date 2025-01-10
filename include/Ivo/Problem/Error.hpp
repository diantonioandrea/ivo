/**
 * @file Error.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Error evaluation.
 * @date 2024-10-19
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef PROBLEM_ERROR
#define PROBLEM_ERROR

#include "./Includes.hpp"

namespace ivo {

    class Error {

        private:

            // Attributes.

            /**
             * @brief Degrees of freedom.
             * 
             */
            const Natural dofs;

            /**
             * @brief (Highest) space degree.
             * 
             */
            const Natural p;

            /**
             * @brief (Highest) time degree.
             * 
             */
            const Natural q;

            /**
             * @brief Space diagram size.
             * 
             */
            const Real h;

            /**
             * @brief Time diagram size.
             * 
             */
            const Real t;

            // Errors.

            /**
             * @brief L2 error.
             * 
             */
            Real l2l2;

            /**
             * @brief L2(L2) errors.
             * 
             */
            std::vector<Real> l2l2s;

            /**
             * @brief L2(T) error.
             * 
             */
            Real l2T;

            /**
             * @brief L2(T) errors.
             * 
             */
            std::vector<Real> l2Ts;

            /**
             * @brief L2(H1) error.
             * 
             */
            Real l2h1;

            /**
             * @brief L2(H1) errors.
             * 
             */
            std::vector<Real> l2h1s;

        public: 

            // Constructor.

            Error(const Mesh21 &, const Equation &, const Vector<Real> &, const std::function<Real (Real, Real, Real)> &, const std::function<std::array<Real, 2> (Real, Real, Real)> &);

            // Output.

            friend std::ostream &operator <<(std::ostream &, const Error &);
    };

}

#endif