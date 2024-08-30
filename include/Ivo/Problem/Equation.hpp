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
            const std::function<Real (Real)> _convection;

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
             * @brief Convection coefficient
             * 
             * @param t Time.
             * @return constexpr Real 
             */
            constexpr Real convection(const Real &t) const { return this->_convection(t); }

            /**
             * @brief Diffusion coefficient
             * 
             * @param t Time.
             * @return constexpr Real 
             */
            constexpr Real diffusion(const Real &t) const { return this->_diffusion(t); }

            /**
             * @brief Reaction coefficient
             * 
             * @param t Time.
             * @return constexpr Real 
             */
            constexpr Real reaction(const Real &t) const { return this->_reaction(t); }

            /**
             * @brief Coefficients.
             * 
             * @param t Time.
             * @return constexpr std::array<Real, 3> 
             */
            constexpr std::array<Real, 3> coefficients(const Real &t) const {
                return {this->_convection(t), this->_diffusion(t), this->_reaction(t)};
            }

            // Constructors.

            Equation(const std::function<Real (Real)> &, const std::function<Real (Real)> &, const std::function<Real (Real)> &);

    };

}

#endif