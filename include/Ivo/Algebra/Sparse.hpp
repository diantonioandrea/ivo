/**
 * @file Sparse.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Sparse matrices.
 * @date 2024-07-21
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef ALGEBRA_SPARSE
#define ALGEBRA_SPARSE

#include "./Matrix.hpp"

namespace ivo {

    /**
     * @brief Sparse matrices.
     * Triple storage (DOK, CSR and CSC)
     * 
     * @tparam T 
     */
    template<Numerical T>
    class Sparse {

        private:

            // Attributes.

            // DOK.

            /**
             * @brief Sparse's entries, DOK.
             * 
             */
            mutable std::map<Natural, T> _entries;

            /**
             * @brief Sparse's rows.
             * 
             */
            const Natural _rows;

            /**
             * @brief Sparse's columns.
             * 
             */
            const Natural _columns;

            // CSR.

            /**
             * @brief CSR state.
             * 
             */
            bool _csr;

            /**
             * @brief CSR Inner vector.
             * 
             */
            std::vector<Natural> _csr_inner;

            /**
             * @brief CSR Outer vector.
             * 
             */
            std::vector<Natural> _csr_outer;

            /**
             * @brief CSR Entries vector.
             * 
             */
            std::vector<T> _csr_entries;

            // CSC.

            /**
             * @brief CSC state.
             * 
             */
            bool _csc;

            /**
             * @brief CSC Inner vector.
             * 
             */
            std::vector<Natural> _csc_inner;

            /**
             * @brief CSC Outer vector.
             * 
             */
            std::vector<Natural> _csc_outer;

            /**
             * @brief CSC Entries vector.
             * 
             */
            std::vector<T> _csc_entries;

        public:

            // Attributes access.

            /**
             * @brief Sparse's CSR structure.
             * 
             * @return std::tuple<std::vector<Natural>, std::vector<Natural>, std::vector<T>> 
             */
            inline std::tuple<std::vector<Natural>, std::vector<Natural>, std::vector<T>> csr() {
                this->_csr_update();
                return {this->_csr_inner, this->_csr_outer, this->_csr_entries};
            }

            /**
             * @brief Sparse's CSC structure.
             * 
             * @return std::tuple<std::vector<Natural>, std::vector<Natural>, std::vector<T>> 
             */
            inline std::tuple<std::vector<Natural>, std::vector<Natural>, std::vector<T>> csc() {
                this->_csc_update();
                return {this->_csc_inner, this->_csc_outer, this->_csc_entries};
            }

            /**
             * @brief Sparse's rows.
             * 
             * @return Natural 
             */
            constexpr Natural rows() const { return this->_rows; }

            /**
             * @brief Sparse's columns.
             * 
             * @return Natural 
             */
            constexpr Natural columns() const { return this->_columns; }

            /**
             * @brief Sparse's size.
             * 
             * @return constexpr Natural 
             */
            constexpr Natural size() const { return this->_rows * this->_columns; }

            // Constructors and copy operators.

            /**
             * @brief Zero constructor.
             * 
             * @param rows 
             * @param columns 
             */
            Sparse(const Natural &rows, const Natural &columns): _rows{rows}, _columns{columns} {
                this->_csr = false;
                this->_csc = false;
            }

            /**
             * @brief Copy constructor.
             * Copies 
             * 
             * @param sparse 
             */
            Sparse(const Sparse &sparse): _rows{sparse._rows}, _columns{sparse._columns} {
                this->_entries = sparse._entries;

                if(sparse._csr) {
                    this->_csr = true;

                    this->_csr_inner.resize(sparse._rows + 1);
                    this->_csr_outer.resize(sparse._csr_outer.size());
                    this->_csr_entries.resize(sparse._csr_entries.size());

                    std::copy(sparse._csr_inner.begin(), sparse._csr_inner.end(), this->_csr_inner.begin());
                    std::copy(sparse._csr_outer.begin(), sparse._csr_outer.end(), this->_csr_outer.begin());
                    std::copy(sparse._csr_entries.begin(), sparse._csr_entries.end(), this->_csr_entries.begin());
                } else
                    this->_csr = false;

                if(sparse._csc) {
                    this->_csc = true;

                    this->_csc_inner.resize(sparse._columns + 1);
                    this->_csc_outer.resize(sparse._csc_outer.size());
                    this->_csc_entries.resize(sparse._csc_entries.size());

                    std::copy(sparse._csc_inner.begin(), sparse._csc_inner.end(), this->_csc_inner.begin());
                    std::copy(sparse._csc_outer.begin(), sparse._csc_outer.end(), this->_csc_outer.begin());
                    std::copy(sparse._csc_entries.begin(), sparse._csc_entries.end(), this->_csc_entries.begin());
                } else
                    this->_csc = false;

            }

            /**
             * @brief Copy operator.
             * 
             * @param sparse 
             * @return Sparse& 
             */
            Sparse &operator =(const Sparse &sparse) {
                #ifndef NDEBUG
                assert(this->_rows == sparse._rows);
                assert(this->_columns == sparse._columns);
                #endif

                this->_entries = sparse._entries;

                if(sparse._csr) {
                    this->_csr = true;

                    this->_csr_inner.resize(this->_rows + 1);
                    this->_csr_outer.resize(sparse._csr_outer.size());
                    this->_csr_entries.resize(sparse._csr_entries.size());

                    std::copy(sparse._csr_inner.begin(), sparse._csr_inner.end(), this->_csr_inner.begin());
                    std::copy(sparse._csr_outer.begin(), sparse._csr_outer.end(), this->_csr_outer.begin());
                    std::copy(sparse._csr_entries.begin(), sparse._csr_entries.end(), this->_csr_entries.begin());
                } else
                    this->_csr = false;

                if(sparse._csc) {
                    this->_csc = true;

                    this->_csc_inner.resize(this->_columns + 1);
                    this->_csc_outer.resize(sparse._csc_outer.size());
                    this->_csc_entries.resize(sparse._csc_entries.size());

                    std::copy(sparse._csc_inner.begin(), sparse._csc_inner.end(), this->_csc_inner.begin());
                    std::copy(sparse._csc_outer.begin(), sparse._csc_outer.end(), this->_csc_outer.begin());
                    std::copy(sparse._csc_entries.begin(), sparse._csc_entries.end(), this->_csc_entries.begin());
                } else
                    this->_csc = false;

                return *this;
            }

            // Subscript operator, legacy scalar access (C++23).

            #if __cplusplus > 202002L

            /**
             * @brief Scalar reference access, legacy.
             * May create an element.
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

                this->_csr = false;
                this->_csc = false;

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

                if(this->_csr)
                    for(Natural h = this->_csr_inner[j]; h < this->_csr_inner[j + 1]; ++h)
                        if(this->_csr_outer[h] == k)
                            return this->_csr_entries[h];

                if(this->_entries.contains(j * this->_columns + k))
                    return this->_entries[j * this->_columns + k];

                return static_cast<T>(0);
            }

            /**
             * @brief Matricial access.
             * 
             * @param J Row indices.
             * @param K Column indices.
             * @return Matrix<T> 
             */
            Matrix<T> operator ()(const std::vector<Natural> &J, const std::vector<Natural> &K) const {
                #ifndef NDEBUG // Integrity check.
                assert(J.size() > 0);
                assert(K.size() > 0);
                for(const auto &j: J)
                    assert(j < this->_rows);
                for(const auto &k: K)
                    assert(k < this->_columns);
                #endif

                Matrix<T> matrix{J.size(), K.size()};
                for(Natural j = 0; j < J.size(); ++j)
                    for(Natural k = 0; k < K.size(); ++k) {
                        if(this->_csr) {
                            for(Natural h = this->_csr_inner[j]; h < this->_csr_inner[j + 1]; ++h)
                                if(this->_csr_outer[h] == k) {
                                    matrix(j, k, this->_csr_entries[h]);
                                }

                            continue;
                        }

                        if(!this->_entries.contains(J[j] * this->_columns + K[k]))
                            continue;

                        matrix(j, k, this->_entries[J[j] * this->_columns + K[k]]);
                    }

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

                this->_csr = false;
                this->_csc = false;

                if(std::abs(scalar) > constants::zero)
                    this->_entries[j * this->_columns + k] = scalar;
            }

            /**
             * @brief Matricial insert.
             * 
             * @param J Row indices.
             * @param K Column indices.
             * @param matrix Matrix.
             */
            void operator ()(const std::vector<Natural> &J, const std::vector<Natural> &K, const Matrix<T> &matrix) {
                #ifndef NDEBUG // Integrity check.
                assert(J.size() == matrix.rows());
                assert(K.size() == matrix.columns());
                for(const auto &j: J)
                    assert(j < this->_rows);
                for(const auto &k: K)
                    assert(k < this->_columns);
                #endif

                this->_csr = false;
                this->_csc = false;

                for(Natural j = 0; j < J.size(); ++j)
                    for(Natural k = 0; k < K.size(); ++k)
                        if(std::abs(matrix(j, k)) > constants::zero)
                            this->_entries[J[j] * this->_columns + K[k]] = matrix(j, k);
            }

            // Access.

            /**
             * @brief Row access.
             * 
             * @param j Index.
             * @return Vector<T> 
             */
            Vector<T> row(const Natural &j) {
                #ifndef NDEBUG // Integrity check.
                assert(j < this->_rows);
                #endif

                // CSR needed.
                this->_csr_update();

                Vector<T> row{this->_columns};
                for(Natural k = this->_csr_inner[j]; k < this->_csr_inner[j + 1]; ++k)
                    row(this->_csr_outer[k], this->_csr_entries[k]);

                return row;
            }

            /**
             * @brief Column access.
             * 
             * @param k Index.
             * @return Vector<T> 
             */
            Vector<T> column(const Natural &k) {
                #ifndef NDEBUG // Integrity check.
                assert(k < this->_columns);
                #endif

                // CSC needed.
                this->_csc_update();

                Vector<T> column{this->_rows};
                for(Natural j = this->_csc_inner[k]; j < this->_csc_inner[k + 1]; ++j)
                    column(this->_csc_outer[j], this->_csc_entries[j]);

                return column;
            }

            // Operations.

            /**
             * @brief +Sparse.
             * 
             * @return Sparse 
             */
            Sparse operator +() const { return *this; }

            /**
             * @brief -Sparse.
             * 
             * @return Sparse 
             */
            Sparse operator -() const {
                Sparse result{*this};

                for(auto &[index, entry]: result._entries)
                    entry = -entry;

                if(result._csr)
                    std::transform(result._csr_entries.begin(), result._csr_entries.end(), result._csr_entries.begin(), [](const T &entry){ return -entry; });
                
                if(result._csc)
                    std::transform(result._csc_entries.begin(), result._csc_entries.end(), result._csc_entries.begin(), [](const T &entry){ return -entry; });

                return result;
            }

            /**
             * @brief Sparse + scalar.
             * 
             * @param scalar Scalar.
             * @return Sparse 
             */
            Sparse operator +(const T &scalar) const {
                Sparse result{*this};

                for(auto &[index, entry]: result._entries)
                    entry += scalar;

                if(result._csr)
                    std::transform(result._csr_entries.begin(), result._csr_entries.end(), result._csr_entries.begin(), [scalar](const T &entry){ return entry + scalar; });
                
                if(result._csc)
                    std::transform(result._csc_entries.begin(), result._csc_entries.end(), result._csc_entries.begin(), [scalar](const T &entry){ return entry + scalar; });

                return result;
            }

            /**
             * @brief Scalar + sparse.
             * 
             * @param scalar Scalar.
             * @param sparse Sparse.
             * @return Sparse 
             */
            friend Sparse operator +(const T &scalar, const Sparse &sparse) {
                Sparse result{sparse};

                for(auto &[index, entry]: result._entries)
                    entry = scalar + entry;

                if(result._csr)
                    std::transform(result._csr_entries.begin(), result._csr_entries.end(), result._csr_entries.begin(), [scalar](const T &entry){ return scalar + entry; });
                
                if(result._csc)
                    std::transform(result._csc_entries.begin(), result._csc_entries.end(), result._csc_entries.begin(), [scalar](const T &entry){ return scalar + entry; });

                return result;
            }

            /**
             * @brief Sparse += scalar.
             * 
             * @param scalar Scalar.
             * @return Sparse& 
             */
            Sparse &operator +=(const T &scalar) {
                for(auto &[index, entry]: this->_entries)
                    entry += scalar;

                if(this->_csr)
                    std::transform(this->_csr_entries.begin(), this->_csr_entries.end(), this->_csr_entries.begin(), [scalar](const T &entry){ return entry + scalar; });
                
                if(this->_csc)
                    std::transform(this->_csc_entries.begin(), this->_csc_entries.end(), this->_csc_entries.begin(), [scalar](const T &entry){ return entry + scalar; });

                return *this;
            }

            /**
             * @brief Sparse - scalar.
             * 
             * @param scalar Scalar.
             * @return Sparse 
             */
            Sparse operator -(const T &scalar) const {
                Sparse result{*this};

                for(auto &[index, entry]: result._entries)
                    entry -= scalar;

                if(result._csr)
                    std::transform(result._csr_entries.begin(), result._csr_entries.end(), result._csr_entries.begin(), [scalar](const T &entry){ return entry - scalar; });
                
                if(result._csc)
                    std::transform(result._csc_entries.begin(), result._csc_entries.end(), result._csc_entries.begin(), [scalar](const T &entry){ return entry - scalar; });

                return result;
            }

            /**
             * @brief Scalar - sparse.
             * 
             * @param scalar Scalar.
             * @param sparse Sparse.
             * @return Sparse 
             */
            friend Sparse operator -(const T &scalar, const Sparse &sparse) {
                Sparse result{sparse};

                for(auto &[index, entry]: result._entries)
                    entry = scalar - entry;

                if(result._csr)
                    std::transform(result._csr_entries.begin(), result._csr_entries.end(), result._csr_entries.begin(), [scalar](const T &entry){ return scalar - entry; });
                
                if(result._csc)
                    std::transform(result._csc_entries.begin(), result._csc_entries.end(), result._csc_entries.begin(), [scalar](const T &entry){ return scalar - entry; });

                return result;
            }

            /**
             * @brief Sparse -= scalar.
             * 
             * @param scalar Scalar.
             * @return Sparse& 
             */
            Sparse &operator -=(const T &scalar) {
                for(auto &[index, entry]: this->_entries)
                    entry -= scalar;

                if(this->_csr)
                    std::transform(this->_csr_entries.begin(), this->_csr_entries.end(), this->_csr_entries.begin(), [scalar](const T &entry){ return entry - scalar; });
                
                if(this->_csc)
                    std::transform(this->_csc_entries.begin(), this->_csc_entries.end(), this->_csc_entries.begin(), [scalar](const T &entry){ return entry - scalar; });

                return *this;
            }

            /**
             * @brief Sparse * scalar.
             * 
             * @param scalar Scalar.
             * @return Sparse 
             */
            Sparse operator *(const T &scalar) const {
                Sparse result{*this};

                for(auto &[index, entry]: result._entries)
                    entry *= scalar;

                if(result._csr)
                    std::transform(result._csr_entries.begin(), result._csr_entries.end(), result._csr_entries.begin(), [scalar](const T &entry){ return entry * scalar; });
                
                if(result._csc)
                    std::transform(result._csc_entries.begin(), result._csc_entries.end(), result._csc_entries.begin(), [scalar](const T &entry){ return entry * scalar; });

                return result;
            }

            /**
             * @brief Scalar * sparse.
             * 
             * @param scalar Scalar.
             * @param sparse Sparse.
             * @return Sparse 
             */
            friend Sparse operator *(const T &scalar, const Sparse &sparse) {
                Sparse result{sparse};

                for(auto &[index, entry]: result._entries)
                    entry = scalar * entry;

                if(result._csr)
                    std::transform(result._csr_entries.begin(), result._csr_entries.end(), result._csr_entries.begin(), [scalar](const T &entry){ return scalar * entry; });
                
                if(result._csc)
                    std::transform(result._csc_entries.begin(), result._csc_entries.end(), result._csc_entries.begin(), [scalar](const T &entry){ return scalar * entry; });

                return result;
            }

            /**
             * @brief Sparse *= scalar.
             * 
             * @param scalar Scalar.
             * @return Sparse& 
             */
            Sparse &operator *=(const T &scalar) {
                for(auto &[index, entry]: this->_entries)
                    entry *= scalar;

                if(this->_csr)
                    std::transform(this->_csr_entries.begin(), this->_csr_entries.end(), this->_csr_entries.begin(), [scalar](const T &entry){ return entry * scalar; });
                
                if(this->_csc)
                    std::transform(this->_csc_entries.begin(), this->_csc_entries.end(), this->_csc_entries.begin(), [scalar](const T &entry){ return entry * scalar; });

                return *this;
            }

            /**
             * @brief Sparse / scalar.
             * 
             * @param scalar Scalar.
             * @return Sparse 
             */
            Sparse operator /(const T &scalar) const {
                Sparse result{*this};

                for(auto &[index, entry]: result._entries)
                    entry /= scalar;

                if(result._csr)
                    std::transform(result._csr_entries.begin(), result._csr_entries.end(), result._csr_entries.begin(), [scalar](const T &entry){ return entry / scalar; });
                
                if(result._csc)
                    std::transform(result._csc_entries.begin(), result._csc_entries.end(), result._csc_entries.begin(), [scalar](const T &entry){ return entry / scalar; });

                return result;
            }

            /**
             * @brief Scalar / sparse.
             * 
             * @param scalar Scalar.
             * @param sparse Sparse.
             * @return Sparse 
             */
            friend Sparse operator /(const T &scalar, const Sparse &sparse) {
                Sparse result{sparse};

                for(auto &[index, entry]: result._entries)
                    entry = scalar / entry;

                if(result._csr)
                    std::transform(result._csr_entries.begin(), result._csr_entries.end(), result._csr_entries.begin(), [scalar](const T &entry){ return scalar / entry; });
                
                if(result._csc)
                    std::transform(result._csc_entries.begin(), result._csc_entries.end(), result._csc_entries.begin(), [scalar](const T &entry){ return scalar / entry; });

                return result;
            }

            /**
             * @brief Sparse /= scalar.
             * 
             * @param scalar Scalar.
             * @return Sparse& 
             */
            Sparse &operator /=(const T &scalar) {
                for(auto &[index, entry]: this->_entries)
                    entry /= scalar;

                if(this->_csr)
                    std::transform(this->_csr_entries.begin(), this->_csr_entries.end(), this->_csr_entries.begin(), [scalar](const T &entry){ return entry / scalar; });
                
                if(this->_csc)
                    std::transform(this->_csc_entries.begin(), this->_csc_entries.end(), this->_csc_entries.begin(), [scalar](const T &entry){ return entry / scalar; });

                return *this;
            }

            /**
             * @brief Sparse + sparse.
             * 
             * @param sparse Sparse matrix.
             * @return Sparse 
             */
            Sparse operator +(const Sparse &sparse) const {
                #ifndef NDEBUG // Integrity check.
                assert(this->_rows == sparse._rows);
                assert(this->_columns == sparse._columns);
                #endif

                Sparse result{*this};

                for(const auto &[index, value]: sparse._entries) {
                    if(!result._entries.contains(index)) {
                        result._entries[index] = value;
                        continue;
                    }

                    if(std::abs(result._entries[index] + value) <= constants::zero) {
                        result._entries.erase(index);
                        continue;
                    }

                    result._entries[index] += value;
                }

                result._csr = false;
                result._csc = false;
                return result;
            }

            /**
             * @brief Sparse += sparse.
             * 
             * @param sparse 
             * @return Sparse& 
             */
            Sparse &operator +=(const Sparse &sparse) {
                #ifndef NDEBUG // Integrity check.
                assert(this->_rows == sparse._rows);
                assert(this->_columns == sparse._columns);
                #endif

                for(const auto &[index, value]: sparse._entries) {
                    if(!this->_entries.contains(index)) {
                        this->_entries[index] = value;
                        continue;
                    }

                    if(std::abs(this->_entries[index] + value) <= constants::zero) {
                        this->_entries.erase(index);
                        continue;
                    }

                    this->_entries[index] += value;
                }

                this->_csr = false;
                this->_csc = false;
                return *this;
            }

            /**
             * @brief Sparse + sparse.
             * 
             * @param sparse Sparse matrix.
             * @return Sparse 
             */
            Sparse operator -(const Sparse &sparse) const {
                #ifndef NDEBUG // Integrity check.
                assert(this->_rows == sparse._rows);
                assert(this->_columns == sparse._columns);
                #endif

                Sparse result{*this};

                for(const auto &[index, value]: sparse._entries) {
                    if(!result._entries.contains(index)) {
                        result._entries[index] = -value;
                        continue;
                    }

                    if(std::abs(result._entries[index] - value) <= constants::zero) {
                        result._entries.erase(index);
                        continue;
                    }

                    result._entries[index] -= value;
                }

                result._csr = false;
                result._csc = false;
                return result;
            }

            /**
             * @brief Sparse -= sparse.
             * 
             * @param sparse Sparse matrix.
             * @return Sparse& 
             */
            Sparse &operator -=(const Sparse &sparse) {
                #ifndef NDEBUG // Integrity check.
                assert(this->_rows == sparse._rows);
                assert(this->_columns == sparse._columns);
                #endif

                for(const auto &[index, value]: sparse._entries) {
                    if(!this->_entries.contains(index)) {
                        this->_entries[index] = -value;
                        continue;
                    }

                    if(std::abs(this->_entries[index] - value) <= constants::zero) {
                        this->_entries.erase(index);
                        continue;
                    }

                    this->_entries[index] -= value;
                }

                this->_csr = false;
                this->_csc = false;
                return *this;
            }

            // Operations (produts).

            /**
             * @brief Sparse * vector.
             * Row x Column product.
             * 
             * @param vector Vector.
             * @return Vector<T> 
             */
            Vector<T> operator *(const Vector<T> &vector) {
                #ifndef NDEBUG // Integrity check.
                assert(this->_columns == vector.size());
                #endif

                // CSR needed.
                this->_csr_update();

                Vector<T> result{this->_rows};
                for(Natural j = 0; j < this->_rows; ++j) {
                    T product = static_cast<T>(0);

                    for(Natural k = this->_csr_inner[j]; k < this->_csr_inner[j + 1]; ++k)
                        product += this->_csr_entries[k] * vector(this->_csr_outer[k]);

                    result(j, product);
                }

                return result;
            }

            /**
             * @brief Vector * matrix.
             * Row x Column product.
             * 
             * @param vector Vector.
             * @param sparse Sparse matrix.
             * @return Vector<T> 
             */
            friend Vector<T> operator *(const Vector<T> &vector, Sparse &sparse) {
                #ifndef NDEBUG // Integrity check.
                assert(sparse._rows == vector.size());
                #endif

                // CSC needed.
                sparse._csc_update();

                Vector<T> result{sparse._columns};
                for(Natural k = 0; k < sparse._columns; ++k) {
                    T product = static_cast<T>(0);

                    for(Natural j = sparse._csc_inner[k]; j < sparse._csc_inner[k + 1]; ++j)
                        product += vector(sparse._csc_outer[j]) * sparse._csc_entries[j];

                    result(k, product);
                }

                return result;
            }

            /**
             * @brief Sparse * sparse.
             * Row x Column product.
             * 
             * @param sparse Sparse matrix.
             * @return Sparse 
             */
            Sparse operator *(Sparse &sparse) const {
                #ifndef NDEBUG // Integrity check.
                assert(this->_columns == sparse._rows);
                #endif

                // CSR and CSC needed.
                this->_csr_update();
                sparse._csc_update();
                
                Sparse result{this->_rows, sparse._columns};

                // [!]

                return result;
            }

            // Output.

            /**
             * @brief Sparse output.
             * 
             * @param ost 
             * @param sparse Sparse matrix.
             * @return std::ostream& 
             */
            friend std::ostream &operator <<(std::ostream &ost, const Sparse &sparse) {
                for(const auto &[index, entry]: sparse._entries) {
                    std::cout << "(" << index / sparse._columns << ", " << index % sparse._columns << "): " << entry;

                    if(index < (*--sparse._entries.end()).first)
                        std::cout << std::endl;
                }

                return ost << std::flush;
            }

        private:

            /**
             * @brief CSR updater.
             * 
             */
            void _csr_update() {
                if(this->_csr)
                    return;

                this->_csr_inner.resize(this->_rows + 1, 0);
                this->_csr_outer.clear();
                this->_csr_entries.clear();

                for(Natural j = 0; j < this->_rows; ++j) {
                    for(const auto &[index, entry]: this->_entries) {
                        if(index < j * this->_columns)
                            continue;

                        if(index >= (j + 1) * this->_columns)
                            break;

                        if(std::abs(entry) > constants::zero) {
                            this->_csr_outer.emplace_back(index % this->_columns);
                            this->_csr_entries.emplace_back(entry);
                        }
                    }

                    this->_csr_inner[j + 1] = this->_csr_outer.size();
                }
            }

            /**
             * @brief CSC updater.
             * 
             */
            void _csc_update() {
                if(this->_csc)
                    return;

                this->_csc_inner.resize(this->_columns + 1, 0);
                this->_csc_outer.clear();
                this->_csc_entries.clear();

                for(Natural k = 0; k < this->_columns; ++k) {
                    for(const auto &[index, entry]: this->_entries) {
                        if(index % this->_columns != k)
                            continue;

                        if(std::abs(entry) > constants::zero) {
                            this->_csc_outer.emplace_back(index / this->_columns);
                            this->_csc_entries.emplace_back(entry);
                        }
                    }

                    this->_csc_inner[k + 1] = this->_csc_outer.size();
                }
            }
    };
    
}

#endif