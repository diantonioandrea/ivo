/**
 * @file Initial.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Initial condition class.
 * @date 2024-09-17
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef PROBLEM_INITIAL
#define PROBLEM_INITIAL

#include "./Includes.hpp"

namespace ivo {

    /**
     * @brief Initial condition class.
     * 
     */
    class Initial {

        private:

            // Attributes.

            /**
             * @brief Initial condition.
             * 
             */
            const std::function<Real (Real, Real)> _condition;

        public:

            // Access.

            /**
             * @brief Initial condition.
             * 
             * @param x 
             * @param y 
             * @return constexpr Real 
             */
            constexpr Real operator ()(const Real &x, const Real &y) const {
                return this->_condition(x, y);
            }

            Vector<Real> operator ()(const Vector<Real> &, const Vector<Real> &) const;

            // Constructor.

            Initial(const std::function<Real (Real, Real)> &);

    };

}

#endif