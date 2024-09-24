/**
 * @file Data.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Data class.
 * @date 2024-09-14
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef PROBLEM_DATA
#define PROBLEM_DATA

#include "./Includes.hpp"

namespace ivo {

    /**
     * @brief Data class.
     * 
     */
    class Data {

        private:

            // Attributes.

            /**
             * @brief Source.
             * 
             */
            const std::function<Real (Real, Real, Real)> _source;

            /**
             * @brief Dirichlet boundary condition.
             * 
             */
            const std::function<Real (Real, Real, Real)> _dirichlet;

            /**
             * @brief Neumann boundary condition.
             * 
             */
            const std::function<Real (Real, Real, Real)> _neumann;

        public:

            // Attributes access.

            /**
             * @brief Source.
             * 
             * @param x 
             * @param y 
             * @param t 
             * @return constexpr Real 
             */
            constexpr Real source(const Real &x, const Real &y, const Real &t) const {
                return this->_source(x, y);
            }

            /**
             * @brief Dirichlet boundary condition.
             * 
             * @param x 
             * @param y 
             * @param t 
             * @return constexpr Real 
             */
            constexpr Real dirichlet(const Real &x, const Real &y, const Real &t) const {
                return this->_dirichlet(x, y);
            }

            /**
             * @brief Neumann boundary condition.
             * 
             * @param x 
             * @param y 
             * @param t 
             * @return constexpr Real 
             */
            constexpr Real neumann(const Real &x, const Real &y, const Real &t) const {
                return this->_neumann(x, y);
            }

            // Constructor.

            Data(const std::function<Real (Real, Real, Real)> &, const std::function<Real (Real, Real, Real)> &, const std::function<Real (Real, Real, Real)> &);

    };

}

#endif