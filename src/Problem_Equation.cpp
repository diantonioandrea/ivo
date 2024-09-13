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
    Equation::Equation(const std::function<std::array<Real, 2> (Real)> &convection, const std::function<Real (Real)> &diffusion, const std::function<Real (Real)> &reaction): _convection{convection}, _diffusion{diffusion}, _reaction{reaction} {}

    // Attributes access.

    /**
     * @brief Convection coefficient, vectorial evaluation.
     * 
     * @param T Time.
     * @return Vector<Real> 
     */
    std::array<Vector<Real>, 2> Equation::convection(const Vector<Real> &T) const {
        std::array<std::vector<Real>, 2> convections;

        for(Natural j = 0; j < T.size(); ++j) {
            std::array<Real, 2> convection = this->_convection(T(j));

            convections[0].emplace_back(convection[0]);
            convections[1].emplace_back(convection[1]);
        }
        
        return {Vector<Real>{convections[0]}, Vector<Real>{convections[1]}};
    }

    /**
     * @brief Diffusion coefficient, vectorial evaluation.
     * 
     * @param T Time.
     * @return Vector<Real> 
     */
    Vector<Real> Equation::diffusion(const Vector<Real> &T) const {
        std::vector<Real> diffusion;
        diffusion.resize(T.size());

        std::vector<Real> time = T.entries();

        std::transform(time.begin(), time.end(), diffusion.begin(), [*this](const Real &t){ return this->_diffusion(t); });
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
        reaction.resize(T.size());

        std::vector<Real> time = T.entries();

        std::transform(time.begin(), time.end(), reaction.begin(), [*this](const Real &t){ return this->_reaction(t); });
        return Vector<Real>{reaction};
    }

}