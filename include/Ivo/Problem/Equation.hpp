/**
 * @file Equation.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Equation class.
 * @date 2024-08-07
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef PROBLEM_EQUATION
#define PROBLEM_EQUATION

#include "./Includes.hpp"

namespace ivo {

    /**
     * @brief Equation class.
     * 
     */
    class Equation {

        private:

            // Attributes.

            /**
             * @brief Convection coefficient.
             * 
             */
            const std::function<std::array<Real, 2> (Real)> _convection;

            /**
             * @brief Diffusion coefficient.
             * 
             */
            const std::function<Real (Real)> _diffusion;

            /**
             * @brief Reaction coefficient.
             * 
             */
            const std::function<Real (Real)> _reaction;

        public:

            // Attributes access.

            /**
             * @brief Convection coefficient.
             * 
             * @param t Time.
             * @return constexpr Real 
             */
            constexpr std::array<Real, 2> convection(const Real &t) const { return this->_convection(t); }

            /**
             * @brief Diffusion coefficient.
             * 
             * @param t Time.
             * @return constexpr Real 
             */
            constexpr Real diffusion(const Real &t) const { return this->_diffusion(t); }

            /**
             * @brief Reaction coefficient.
             * 
             * @param t Time.
             * @return constexpr Real 
             */
            constexpr Real reaction(const Real &t) const { return this->_reaction(t); }

            std::array<Vector<Real>, 2> convection(const Vector<Real> &) const;
            Vector<Real> diffusion(const Vector<Real> &) const;
            Vector<Real> reaction(const Vector<Real> &) const;

            // Constructor.

            Equation(const std::function<std::array<Real, 2> (Real)> &, const std::function<Real (Real)> &, const std::function<Real (Real)> &);

    };

}

#endif