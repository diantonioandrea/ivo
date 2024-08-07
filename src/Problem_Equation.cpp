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

    // Constructors.
    
    /**
     * @brief Default constructor.
     * 
     * @param convection Convection coefficient.
     * @param diffusion Diffusion coefficient.
     * @param reaction Reaction coefficient.
     */
    Equation::Equation(const Real &convection, const Real &diffusion, const Real &reaction): _convection{convection}, _diffusion{diffusion}, _reaction{reaction} {
        #ifndef NDEBUG // Integrity check.
        assert(std::abs(convection) + std::abs(diffusion) + std::abs(reaction) > NUMERICAL_ZERO);
        #endif
    }

    // Output.
    
    /**
     * @brief Equation output.
     * 
     * @param ost 
     * @param equation Equation.
     * @return std::ostream& 
     */
    std::ostream &operator <<(std::ostream &ost, const Equation &equation) {
        return ost << "(C: " << equation._convection << ", D: " << equation._diffusion << ", R: " << equation._reaction << ")" << std::flush;
    }

}