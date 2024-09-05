/**
 * @file Matrix.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Matrix related methods.
 * @date 2024-09-05
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef ALGEBRA_METHODS_MATRIX
#define ALGEBRA_METHODS_MATRIX

#include "../Matrix.hpp"

namespace ivo {

    /**
     * @brief Kronecker product.
     * 
     * @tparam T Numerical type.
     * @param X First matrix.
     * @param Y Second matrix.
     * @return Matrix<T> 
     */
    template<Numerical T>
    Matrix<T> kronecker(const Matrix<T> &X, const Matrix<T> &Y) {
        Matrix<T> XY{X.rows() * Y.rows(), X.columns() * Y.columns()};

        for(Natural jx = 0; jx < X.rows(); ++jx)
            for(Natural kx = 0; kx < X.columns(); ++kx)
                for(Natural jy = 0; jy < Y.rows(); ++jy)
                    for(Natural ky = 0; ky < Y.columns; ++ky)
                        XY(jx * X.rows() + jy, kx * X.columns() + ky, X(jx, kx) * Y(jy, ky));

        return XY;
    }

}

#endif