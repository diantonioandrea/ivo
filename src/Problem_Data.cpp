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
    Data::Data(const std::function<Real (Real, Real)> &source, const std::function<Real (Real, Real)> &dirichlet, const std::function<Real (Real, Real)> &neumann): _source{source}, _dirichlet{dirichlet}, _neumann{neumann} {}

    // Attributes access.

    /**
     * @brief Source.
     * 
     * @param X 
     * @param Y 
     * @return Vector<Real> 
     */
    Vector<Real> Data::source(const Vector<Real> &X, const Vector<Real> &Y) const {
        #ifndef NDEBUG // Integrity check.
        assert(X.size() == Y.size());
        #endif

        std::vector<Real> source;
        source.resize(X.size());

        std::vector<Real> x = X.entries();
        std::vector<Real> y = Y.entries();

        std::transform(x.begin(), x.end(), y.begin(), source.begin(), [*this](const Real &x, const Real &y){ return this->_source(x, y); });
        return Vector<Real>{source};
    }
    
    /**
     * @brief Dirichlet boundary condition.
     * 
     * @param X 
     * @param Y 
     * @return Vector<Real> 
     */
    Vector<Real> Data::dirichlet(const Vector<Real> &X, const Vector<Real> &Y) const {
        #ifndef NDEBUG // Integrity check.
        assert(X.size() == Y.size());
        #endif

        std::vector<Real> dirichlet;
        dirichlet.resize(X.size());

        std::vector<Real> x = X.entries();
        std::vector<Real> y = Y.entries();

        std::transform(x.begin(), x.end(), y.begin(), dirichlet.begin(), [*this](const Real &x, const Real &y){ return this->_dirichlet(x, y); });
        return Vector<Real>{dirichlet};
    }

    /**
     * @brief Neumann boundary condition.
     * 
     * @param X 
     * @param Y 
     * @return Vector<Real> 
     */
    Vector<Real> Data::neumann(const Vector<Real> &X, const Vector<Real> &Y) const {
        #ifndef NDEBUG // Integrity check.
        assert(X.size() == Y.size());
        #endif

        std::vector<Real> neumann;
        neumann.resize(X.size());

        std::vector<Real> x = X.entries();
        std::vector<Real> y = Y.entries();

        std::transform(x.begin(), x.end(), y.begin(), neumann.begin(), [*this](const Real &x, const Real &y){ return this->_neumann(x, y); });
        return Vector<Real>{neumann};
    }
    
}