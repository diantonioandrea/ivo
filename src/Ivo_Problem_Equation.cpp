/**
 * @file Problem_Equation.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Problem/Equation.hpp implementation.
 * @date 2024-08-07
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Ivo.hpp>

namespace ivo {

    // Constructor.
    
    /**
     * @brief Default constructor.
     * 
     * @param convection Convection coefficient.
     * @param diffusion Diffusion coefficient.
     * @param reaction Reaction coefficient.
     */
    Equation::Equation(const std::function<std::array<Real, 2> (Real, Real, Real)> &convection, const Real &diffusion, const std::function<Real (Real, Real, Real)> &reaction): _convection{convection}, _diffusion{diffusion}, _reaction{reaction} {}

}