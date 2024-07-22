/**
 * @file Matrix.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Matrices.
 * @date 2024-07-20
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef ALGEBRA_MATRIX
#define ALGEBRA_MATRIX

#include "./Includes.hpp"
#include "./Vector.hpp"

namespace ivo {
    
    /**
     * @brief Dinamically allocated matrices.
     * 
     * @tparam T Numerical type.
     */
    template<Numerical T>
    class Matrix {

        private:

            // Attributes.

            /**
             * @brief Matrix' entries, by row.
             * 
             */
            std::vector<T> _entries;

            /**
             * @brief Matrix' rows.
             * 
             */
            const Natural _rows;

            /**
             * @brief Matrix' columns.
             * 
             */
            const Natural _columns;

        public:

            // Attributes access.

            /**
             * @brief Matrix' entries, by row.
             * 
             * @return std::vector<T> 
             */
            inline std::vector<T> entries() const { return this->_entries; }

            /**
             * @brief Matrix' rows.
             * 
             * @return Natural 
             */
            constexpr Natural rows() const { return this->_rows; }

            /**
             * @brief Matrix' columns.
             * 
             * @return Natural 
             */
            constexpr Natural columns() const { return this->_columns; }

            /**
             * @brief Matrix' size.
             * 
             * @return constexpr Natural 
             */
            constexpr Natural size() const { return this->_rows * this->_columns; }

            // Constructors and copy operators.

            /**
             * @brief Zero constructor.
             * 
             * @param rows Matrix' rows.
             * @param columns Matrix' columns.
             */
            Matrix(const Natural &rows, const Natural &columns): _rows{rows}, _columns{columns} {
                #ifndef NDEBUG // Integrity check.
                assert(rows > 0);
                assert(columns > 0);
                #endif

                this->_entries.resize(rows * columns, static_cast<T>(0));
            }

            /**
             * @brief Scalar constructor.
             * 
             * @param rows Matrix' rows.
             * @param columns Matrix' columns.
             * @param scalar Scalar.
             */
            Matrix(const Natural &rows, const Natural &columns, const T &scalar): _rows{rows}, _columns{columns} {
                #ifndef NDEBUG // Integrity check.
                assert(rows > 0);
                assert(columns > 0);
                #endif

                this->_entries.resize(rows * columns, scalar);
            }

            /**
             * @brief Scalar copy operator.
             * 
             * @param scalar Scalar.
             * @return Matrix& 
             */
            Matrix &operator =(const T &scalar) {
                this->_entries.resize(this->_rows * this->_columns, scalar);
                return *this;
            }

            /**
             * @brief Standard vector constructor.
             * 
             * @param rows Matrix' rows.
             * @param columns Matrix' columns.
             * @param vector Standard vector.
             */
            Matrix(const Natural &rows, const Natural &columns, const std::vector<T> &vector): _rows{rows}, _columns{columns} {
                #ifndef NDEBUG // Integrity check.
                assert(rows > 0);
                assert(columns > 0);
                assert(rows * columns == vector.size());
                #endif

                this->_entries.resize(rows * columns, static_cast<T>(0));
                std::copy(vector.begin(), vector.end(), this->_entries.begin());
            }

            /**
             * @brief Standard vector copy operator.
             * 
             * @param vector Standard vector.
             */
            Matrix &operator =(const std::vector<T> &vector) {
                #ifndef NDEBUG // Integrity check.
                assert(this->_rows * this->_columns == vector.size());
                #endif

                std::copy(vector.begin(), vector.end(), this->_entries.begin());
            }

            /**
             * @brief Vector constructor.
             * 
             * @param rows Matrix' rows.
             * @param columns Matrix' columns.
             * @param vector Vector.
             */
            Matrix(const Natural &rows, const Natural &columns, const Vector<T> &vector): Matrix(rows, columns, vector.entries()) {}

            /**
             * @brief Vector copy operator.
             * 
             * @param vector Vector.
             */
            Matrix &operator =(const Vector<T> &vector) {
                return *this = vector.entries();
            }

            /**
             * @brief Copy constructor.
             * 
             * @param matrix Matrix.
             */
            Matrix(const Matrix &matrix): _rows{matrix._rows}, _columns{matrix._columns} {
                this->_entries.resize(this->_rows * this->_columns, static_cast<T>(0));
                std::copy(matrix._entries.begin(), matrix._entries.end(), this->_entries.begin());
            }

            /**
             * @brief Copy operator.
             * 
             * @param matrix Matrix.
             * @return Matrix& 
             */
            Matrix &operator =(const Matrix &matrix) {
                #ifndef NDEBUG // Integrity check.
                assert(this->_rows == matrix._rows);
                assert(this->_columns == matrix._columns);
                #endif

                std::copy(matrix._entries.begin(), matrix._entries.end(), this->_entries.begin());
                return *this;
            }

            // Subscript operator, legacy scalar access (C++23).

            #if __cplusplus > 202002L

            /**
             * @brief Scalar reference access, legacy.
             * 
             * @param j Row index.
             * @param k Column index.
             * @return T& 
             */
            inline T &operator [](const Natural &j, const Natural &k) {
                #ifndef NDEBUG // Integrity check.
                assert(j < this->_rows);
                assert(k < this->_columns);
                #endif

                return this->_entries[j * this->_columns + k];
            }

            #endif

            // Call operator, subscript behaviour.

            /**
             * @brief Scalar access.
             * 
             * @param j Row index.
             * @param k Column index.
             * @return T 
             */
            inline T operator ()(const Natural &j, const Natural &k) const {
                #ifndef NDEBUG // Integrity check.
                assert(j < this->_rows);
                assert(k < this->_columns);
                #endif

                return this->_entries[j * this->_columns + k];
            }

            /**
             * @brief Matricial access.
             * 
             * @param J Row indices.
             * @param K Column indices.
             * @return Matrix 
             */
            Matrix operator ()(const std::vector<Natural> &J, const std::vector<Natural> &K) const {
                #ifndef NDEBUG // Integrity check.
                assert(J.size() > 0);
                assert(K.size() > 0);
                for(const auto &j: J)
                    assert(j < this->_rows);
                for(const auto &k: K)
                    assert(k < this->_columns);
                #endif

                Matrix matrix{J.size(), K.size()};
                for(Natural j = 0; j < J.size(); ++j)
                    for(Natural k = 0; k < K.size(); ++k)
                        matrix._entries[j * matrix._columns + k] = this->_entries[J[j] * this->_columns + K[k]];

                return matrix;
            }

            /**
             * @brief Scalar insert.
             * 
             * @param j Row index.
             * @param k Column index.
             * @param scalar Scalar.
             */
            void operator ()(const Natural &j, const Natural &k, const T &scalar) {
                #ifndef NDEBUG // Integrity check.
                assert(j < this->_rows);
                assert(k < this->_columns);
                #endif

                this->_entries[j * this->_columns + k] = scalar;
            }

            /**
             * @brief Matricial insert.
             * 
             * @param J Row indices.
             * @param K Column indices.
             * @param matrix Matrix.
             */
            void operator ()(const std::vector<Natural> &J, const std::vector<Natural> &K, const Matrix &matrix) {
                #ifndef NDEBUG // Integrity check.
                assert(J.size() == matrix._rows);
                assert(K.size() == matrix._columns);
                for(const auto &j: J)
                    assert(j < this->_rows);
                for(const auto &k: K)
                    assert(k < this->_columns);
                #endif

                for(Natural j = 0; j < J.size(); ++j)
                    for(Natural k = 0; k < K.size(); ++k)
                        this->_entries[J[j] * this->_columns + K[k]] = matrix._entries[j * matrix._columns + k];
            }

            // Access.

            /**
             * @brief Row access.
             * 
             * @param j Index.
             * @return Vector<T> 
             */
            Vector<T> row(const Natural &j) const {
                #ifndef NDEBUG // Integrity check.
                assert(j < this->_rows);
                #endif

                // std::vector for built-in copy.
                std::vector<T> row(this->_columns, static_cast<T>(0));
                std::copy(this->_entries.begin() + j * this->_columns, this->_entries.begin() + (j + 1) * this->_columns, row.begin());
                return row;
            }

            /**
             * @brief Row insert.
             * 
             * @param j Index.
             * @param vector 
             */
            void row(const Natural &j, const Vector<T> &vector) {
                #ifndef NDEBUG // Integrity check.
                assert(j < this->_rows);
                assert(this->_columns == vector.size());
                #endif

                // std::vector for built-in copy.
                std::vector<T> row = vector.entries();
                std::copy(row.begin(), row.end(), this->_entries.begin() + j * this->_columns);
            }

            /**
             * @brief Column access.
             * 
             * @param k Index.
             * @return Vector<T> 
             */
            Vector<T> column(const Natural &k) const {
                #ifndef NDEBUG // Integrity check.
                assert(k < this->_columns);
                #endif

                Vector<T> column{this->_rows};
                for(Natural j = 0; j < this->_rows; ++j)
                    column(j, this->_entries[j * this->_columns + k]);

                return column;
            }

            /**
             * @brief Column insert.
             * 
             * @param k 
             * @param vector 
             */
            void column(const Natural &k, const Vector<T> &vector) {
                #ifndef NDEBUG // Integrity check.
                assert(k < this->_columns);
                assert(this->_rows == vector.size());
                #endif

                for(Natural j = 0; j < this->_rows; ++j)
                    this->_entries[j * this->_columns + k] = vector(j);
            }

            // Operations.

            /**
             * @brief +Matrix.
             * 
             * @return Matrix 
             */
            Matrix operator +() const { return *this; }

            /**
             * @brief -Matrix.
             * 
             * @return Matrix 
             */
            Matrix operator -() const {
                Matrix result{this->_rows, this->_columns};
                std::transform(this->_entries.begin(), this->_entries.end(), result._entries.begin(), [](const T &entry){ return -entry; });
                return result;
            }

            /**
             * @brief Matrix + scalar.
             * 
             * @param scalar Scalar.
             * @return Matrix 
             */
            Matrix operator +(const T &scalar) const {
                Matrix result{this->_rows, this->_columns};
                std::transform(this->_entries.begin(), this->_entries.end(), result._entries.begin(), [scalar](const T &entry){ return entry + scalar; });
                return result;
            }

            /**
             * @brief Scalar + matrix.
             * 
             * @param scalar 
             * @param matrix 
             * @return Matrix 
             */
            friend Matrix operator +(const T &scalar, const Matrix &matrix) {
                Matrix result{matrix._rows, matrix._columns};
                std::transform(matrix._entries.begin(), matrix._entries.end(), result._entries.begin(), [scalar](const T &entry){ return scalar + entry; });
                return result;
            }

            /**
             * @brief Matrix += scalar.
             * 
             * @param scalar Scalar.
             * @return Matrix 
             */
            Matrix operator +=(const T &scalar) {
                std::transform(this->_entries.begin(), this->_entries.end(), this->_entries.begin(), [scalar](const T &entry){ return entry + scalar; });
                return *this;
            }

            /**
             * @brief Matrix - scalar.
             * 
             * @param scalar Scalar.
             * @return Matrix 
             */
            Matrix operator -(const T &scalar) const {
                Matrix result{this->_rows, this->_columns};
                std::transform(this->_entries.begin(), this->_entries.end(), result._entries.begin(), [scalar](const T &entry){ return entry - scalar; });
                return result;
            }

            /**
             * @brief Scalar - matrix.
             * 
             * @param scalar 
             * @param matrix 
             * @return Matrix 
             */
            friend Matrix operator -(const T &scalar, const Matrix &matrix) {
                Matrix result{matrix._rows, matrix._columns};
                std::transform(matrix._entries.begin(), matrix._entries.end(), result._entries.begin(), [scalar](const T &entry){ return scalar - entry; });
                return result;
            }

            /**
             * @brief Matrix -= scalar.
             * 
             * @param scalar Scalar.
             * @return Matrix 
             */
            Matrix operator -=(const T &scalar) {
                std::transform(this->_entries.begin(), this->_entries.end(), this->_entries.begin(), [scalar](const T &entry){ return entry - scalar; });
                return *this;
            }

            /**
             * @brief Matrix * scalar.
             * 
             * @param scalar Scalar.
             * @return Matrix 
             */
            Matrix operator *(const T &scalar) const {
                Matrix result{this->_rows, this->_columns};
                std::transform(this->_entries.begin(), this->_entries.end(), result._entries.begin(), [scalar](const T &entry){ return entry * scalar; });
                return result;
            }

            /**
             * @brief Scalar * matrix.
             * 
             * @param scalar 
             * @param matrix 
             * @return Matrix 
             */
            friend Matrix operator *(const T &scalar, const Matrix &matrix) {
                Matrix result{matrix._rows, matrix._columns};
                std::transform(matrix._entries.begin(), matrix._entries.end(), result._entries.begin(), [scalar](const T &entry){ return scalar * entry; });
                return result;
            }

            /**
             * @brief Matrix *= scalar.
             * 
             * @param scalar Scalar.
             * @return Matrix 
             */
            Matrix operator *=(const T &scalar) {
                std::transform(this->_entries.begin(), this->_entries.end(), this->_entries.begin(), [scalar](const T &entry){ return entry * scalar; });
                return *this;
            }
            
            /**
             * @brief Matrix / scalar.
             * 
             * @param scalar Scalar.
             * @return Matrix 
             */
            Matrix operator /(const T &scalar) const {
                Matrix result{this->_rows, this->_columns};
                std::transform(this->_entries.begin(), this->_entries.end(), result._entries.begin(), [scalar](const T &entry){ return entry / scalar; });
                return result;
            }

            /**
             * @brief Scalar / matrix.
             * 
             * @param scalar 
             * @param matrix 
             * @return Matrix 
             */
            friend Matrix operator /(const T &scalar, const Matrix &matrix) {
                Matrix result{matrix._rows, matrix._columns};
                std::transform(matrix._entries.begin(), matrix._entries.end(), result._entries.begin(), [scalar](const T &entry){ return scalar / entry; });
                return result;
            }

            /**
             * @brief Matrix /= scalar.
             * 
             * @param scalar Scalar.
             * @return Matrix 
             */
            Matrix operator /=(const T &scalar) {
                std::transform(this->_entries.begin(), this->_entries.end(), this->_entries.begin(), [scalar](const T &entry){ return entry / scalar; });
                return *this;
            }

            /**
             * @brief Matrix + matrix.
             * 
             * @param matrix Matrix.
             * @return Matrix 
             */
            Matrix operator +(const Matrix &matrix) const {
                #ifndef NDEBUG // Integrity check.
                assert(this->_rows == matrix._rows);
                assert(this->_columns == matrix._columns);
                #endif

                Matrix result{this->_rows, this->_columns};
                std::transform(this->_entries.begin(), this->_entries.end(), matrix._entries.begin(), result._entries.begin(), [](const T &t_entry, const T &m_entry){ return t_entry + m_entry; });
                return result;
            }

            /**
             * @brief Matrix += matrix.
             * 
             * @param matrix 
             * @return Matrix& 
             */
            Matrix &operator +=(const Matrix &matrix) {
                #ifndef NDEBUG // Integrity check.
                assert(this->_rows == matrix._rows);
                assert(this->_columns == matrix._columns);
                #endif

                std::transform(this->_entries.begin(), this->_entries.end(), matrix._entries.begin(), this->_entries.begin(), [](const T &t_entry, const T &m_entry){ return t_entry + m_entry; });
                return *this;
            }

            /**
             * @brief Matrix - matrix.
             * 
             * @param matrix Matrix.
             * @return Matrix 
             */
            Matrix operator -(const Matrix &matrix) const {
                #ifndef NDEBUG // Integrity check.
                assert(this->_rows == matrix._rows);
                assert(this->_columns == matrix._columns);
                #endif

                Matrix result{this->_rows, this->_columns};
                std::transform(this->_entries.begin(), this->_entries.end(), matrix._entries.begin(), result._entries.begin(), [](const T &t_entry, const T &m_entry){ return t_entry - m_entry; });
                return result;
            }

            /**
             * @brief Matrix -= matrix.
             * 
             * @param matrix 
             * @return Matrix& 
             */
            Matrix &operator -=(const Matrix &matrix) {
                #ifndef NDEBUG // Integrity check.
                assert(this->_rows == matrix._rows);
                assert(this->_columns == matrix._columns);
                #endif

                std::transform(this->_entries.begin(), this->_entries.end(), matrix._entries.begin(), this->_entries.begin(), [](const T &t_entry, const T &m_entry){ return t_entry - m_entry; });
                return *this;
            }

            // Operations (produts).

            /**
             * @brief Matrix * vector.
             * Row x Column product.
             * 
             * @param vector Vector.
             * @return Vector<T> 
             */
            Vector<T> operator *(const Vector<T> &vector) const {
                #ifndef NDEBUG // Integrity check.
                assert(this->_columns == vector.size());
                #endif

                Vector<T> result{this->_rows};
                for(Natural j = 0; j < this->_rows; ++j) {
                    T product = static_cast<T>(0);

                    for(Natural k = 0; k < this->_columns; ++k)
                        product += this->_entries[j * this->_columns + k] * vector(k);
                    
                    result(j, product);
                }

                return result;
            }

            /**
             * @brief Vector * matrix.
             * Row x Column product.
             * 
             * @param vector Vector.
             * @param matrix Matrix.
             * @return Vector<T> 
             */
            friend Vector<T> operator *(const Vector<T> &vector, const Matrix &matrix) {
                #ifndef NDEBUG // Integrity check.
                assert(matrix._rows == vector.size());
                #endif

                Vector<T> result{matrix._columns};
                for(Natural k = 0; k < matrix._columns; ++k) {
                    T product = static_cast<T>(0);

                    for(Natural j = 0; j < matrix._rows; ++j)
                        product += vector(j) * matrix._entries[j * matrix._columns + k];
                    
                    result(k, product);
                }

                return result;
            }

            /**
             * @brief Matrix * matrix.
             * Row x Column product.
             * 
             * @param matrix Matrix.
             * @return Matrix 
             */
            Matrix operator *(const Matrix &matrix) {
                #ifndef NDEBUG // Integrity check.
                assert(this->_columns == matrix._rows);
                #endif

                Matrix result{this->_rows, matrix._columns};
                for(Natural j = 0; j < this->_rows; ++j)
                    for(Natural k = 0; k < matrix._columns; ++k) {
                        T product = static_cast<T>(0);

                        for(Natural i = 0; i < this->_columns; ++i)
                            product += this->_entries[j * this->_columns + i] * matrix._entries[i * matrix._columns + k];

                        result._entries[j * matrix._columns + k] = product;
                    }

                return result;
            }

            // Output.

            /**
             * @brief Matrix output.
             * 
             * @param ost 
             * @param matrix Matrix.
             * @return std::ostream& 
             */
            friend std::ostream &operator <<(std::ostream &ost, const Matrix &matrix) {
                for(Natural j = 0; j < matrix._rows; ++j) {
                    for(Natural k = 0; k < matrix._columns; ++k)
                        ost << matrix._entries[j * matrix._columns + k] << " ";

                    if(j < matrix._rows - 1)
                        ost << std::endl;
                }

                return ost << std::flush;
            }
    };

}

#endif