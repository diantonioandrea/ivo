/**
 * @file Problem_Data.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Problem/Data.hpp implementation.
 * @date 2024-09-14
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
     * @param source Source.
     * @param dirichlet Dirichlet boundary condition.
     * @param neumann Neumann boundary condition.
     */
    Data::Data(const std::function<Real (Real, Real, Real)> &source, const std::function<Real (Real, Real, Real)> &dirichlet, const std::function<Real (Real, Real, Real)> &neumann): _source{source}, _dirichlet{dirichlet}, _neumann{neumann} {}
    
}