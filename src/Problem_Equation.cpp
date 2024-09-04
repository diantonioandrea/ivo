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
    Equation::Equation(const std::function<Real (Real)> &convection, const std::function<Real (Real)> &diffusion, const std::function<Real (Real)> &reaction): _convection{convection}, _diffusion{diffusion}, _reaction{reaction} {}

    // Attributes access.

    /**
     * @brief Convection coefficient, vectorial evaluation.
     * 
     * @param T Time.
     * @return Vector<Real> 
     */
    Vector<Real> Equation::convection(const Vector<Real> &T) const {
        std::vector<Real> convection;
        std::vector<Real> time = T.entries();

        std::transform(time.begin(), time.end(), convection.begin(), [*this](const Real &t){ return this->_convection(t); });
        return Vector<Real>{convection};
    }

    /**
     * @brief Diffusion coefficient, vectorial evaluation.
     * 
     * @param T Time.
     * @return Vector<Real> 
     */
    Vector<Real> Equation::diffusion(const Vector<Real> &T) const {
        std::vector<Real> diffusion;
        std::vector<Real> time = T.entries();

        std::transform(time.begin(), time.end(), diffusion.begin(), [*this](const Real &t){ return this->_convection(t); });
        return Vector<Real>{diffusion};
    }

    /**
     * @brief Reaction coefficient, vectorial evaluation.
     * 
     * @param T Time.
     * @return Vector<Real> 
     */
    Vector<Real> Equation::reaction(const Vector<Real> &T) const {
        std::vector<Real> reaction;
        std::vector<Real> time = T.entries();

        std::transform(time.begin(), time.end(), reaction.begin(), [*this](const Real &t){ return this->_convection(t); });
        return Vector<Real>{reaction};
    }

    /**
     * @brief Coefficients, vectorial evaluation.
     * 
     * @param T Time.
     * @return std::array<Vector<Real>, 3> 
     */
    std::array<Vector<Real>, 3> Equation::coefficients(const Vector<Real> &T) const {
        std::vector<Real> convection, diffusion, reaction;
        std::vector<Real> time = T.entries();

        std::transform(time.begin(), time.end(), convection.begin(), [*this](const Real &t){ return this->_convection(t); });
        std::transform(time.begin(), time.end(), diffusion.begin(), [*this](const Real &t){ return this->_convection(t); });
        std::transform(time.begin(), time.end(), reaction.begin(), [*this](const Real &t){ return this->_convection(t); });

        return {Vector<Real>{convection}, Vector<Real>{diffusion}, Vector<Real>{reaction}};
    }

}