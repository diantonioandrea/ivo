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
            const std::function<Real (Real, Real)> _source;

            /**
             * @brief Dirichlet boundary condition.
             * 
             */
            const std::function<Real (Real, Real)> _dirichlet;

            /**
             * @brief Neumann boundary condition.
             * 
             */
            const std::function<Real (Real, Real)> _neumann;

        public:

            // Attributes access.

            /**
             * @brief Source.
             * 
             * @param x 
             * @param y 
             * @return constexpr Real 
             */
            constexpr Real source(const Real &x, const Real &y) const {
                return this->_source(x, y);
            }

            /**
             * @brief Dirichlet boundary condition.
             * 
             * @param x 
             * @param y 
             * @return constexpr Real 
             */
            constexpr Real dirichlet(const Real &x, const Real &y) const {
                return this->_dirichlet(x, y);
            }

            /**
             * @brief Neumann boundary condition.
             * 
             * @param x 
             * @param y 
             * @return constexpr Real 
             */
            constexpr Real neumann(const Real &x, const Real &y) const {
                return this->_neumann(x, y);
            }

            Vector<Real> source(const Vector<Real> &, const Vector<Real> &) const;
            Vector<Real> dirichlet(const Vector<Real> &, const Vector<Real> &) const;
            Vector<Real> neumann(const Vector<Real> &, const Vector<Real> &) const;

            // Constructor.

            Data(const std::function<Real (Real, Real)> &, const std::function<Real (Real, Real)> &, const std::function<Real (Real, Real)> &);

    };

}

#endif