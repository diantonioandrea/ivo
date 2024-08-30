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
            const Real _convection = 0.0L;

            /**
             * @brief Diffusion coefficient.
             * 
             */
            const Real _diffusion = 0.0L;

            /**
             * @brief Reaction coefficient.
             * 
             */
            const Real _reaction = 0.0L;

        public:

            // Attributes access.

            /**
             * @brief Convection coefficient.
             * 
             * @return constexpr Real 
             */
            constexpr Real convection() const { return this->_convection; }

            /**
             * @brief Diffusion coefficient.
             * 
             * @return constexpr Real 
             */
            constexpr Real diffusion() const { return this->_diffusion; }

            /**
             * @brief Reaction coefficient.
             * 
             * @return constexpr Real 
             */
            constexpr Real reaction() const { return this->_reaction; }

            // Constructors.

            Equation(const Real &, const Real &, const Real &);

            // Output.

            friend std::ostream &operator <<(std::ostream &, const Equation &);

    };

}

#endif