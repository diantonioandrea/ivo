/**
 * @file Problem_Initial.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-09-17
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
     * @param condition Initial condition.
     */
    Initial::Initial(const std::function<Real (Real, Real)> &condition): _condition{condition} {}

    // Access.

    /**
     * @brief Initial condition.
     * 
     * @param X 
     * @param Y 
     * @return Vector<Real> 
     */
    Vector<Real> Initial::operator ()(const Vector<Real> &X, const Vector<Real> &Y) const {
        #ifndef NDEBUG // Integrity check.
        assert(X.size() == Y.size());
        #endif

        std::vector<Real> condition;
        condition.resize(X.size());

        std::vector<Real> x = X.entries();
        std::vector<Real> y = Y.entries();

        std::transform(x.begin(), x.end(), y.begin(), condition.begin(), [*this](const Real &x, const Real &y){ return this->_condition(x, y); });
        return Vector<Real>{condition};
    }
    
}