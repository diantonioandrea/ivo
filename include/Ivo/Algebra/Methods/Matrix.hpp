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

    namespace internal {
        
        /**
         * @brief Scales a matrix' rows by a vector.
         * 
         * @tparam T Numerical type.
         * @param vector Vector.
         * @param matrix Matrix.
         * @return Matrix<T> 
         */
        template<Numerical T>
        Matrix<T> r_scale(const Vector<T> &vector, const Matrix<T> &matrix) {
            #ifndef NDEBUG // Integrity check.
            assert(vector.size() == matrix.columns());
            #endif

            Matrix<T> scaled{matrix.rows(), matrix.columns()};

            for(Natural j = 0; j < matrix.rows(); ++j)
                scaled.row(j, vector * matrix.row(j));

            return scaled;
        }

        /**
         * @brief Scales a matrix' columns by a vector.
         * 
         * @tparam T Numerical type.
         * @param vector Vector.
         * @param matrix Matrix.
         * @return Matrix<T> 
         */
        template<Numerical T>
        Matrix<T> c_scale(const Vector<T> &vector, const Matrix<T> &matrix) {
            #ifndef NDEBUG // Integrity check.
            assert(vector.size() == matrix.rows());
            #endif

            Matrix<T> scaled{matrix.rows(), matrix.columns()};

            for(Natural j = 0; j < matrix.columns(); ++j)
                scaled.column(j, vector * matrix.column(j));

            return scaled;
        }

    }

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
                    for(Natural ky = 0; ky < Y.columns(); ++ky)
                        XY(jx * Y.rows() + jy, kx * Y.columns() + ky, X(jx, kx) * Y(jy, ky));

        return XY;
    }

}

#endif