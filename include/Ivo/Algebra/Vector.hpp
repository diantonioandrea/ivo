/**
 * @file Vector.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Vectors.
 * @date 2024-07-19
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef ALGEBRA_VECTOR
#define ALGEBRA_VECTOR

#include "./Includes.hpp"

namespace ivo {

    /**
     * @brief Dynamically allocated vectors.
     * 
     * @tparam T Numerical type.
     */
    template<Numerical T>
    class Vector {

        private:

            // Attributes.

            /**
             * @brief Vector's entries.
             * 
             */
            std::vector<T> _entries;

            /**
             * @brief Vector's size.
             * 
             */
            const Natural _size;

        public:

            // Attributes access.

            /**
             * @brief Vector's entries.
             * 
             * @return std::vector<T> 
             */
            inline std::vector<T> entries() const { return this->_entries; }

            /**
             * @brief Vector's size.
             * 
             * @return constexpr Natural 
             */
            constexpr Natural size() const { return this->_size; }

            // Constructors and copy operators.

            /**
             * @brief Zero constructor.
             * 
             * @param size Vector's size.
             */
            Vector(const Natural &size): _size{size} {
                #ifndef NDEBUG // Integrity check.
                assert(size > 0);
                #endif

                this->_entries.resize(size, static_cast<T>(0));
            }

            /**
             * @brief Scalar constructor.
             * 
             * @param size Vector's size.
             * @param scalar Scalar.
             */
            Vector(const Natural &size, const T &scalar): _size{size} {
                #ifndef NDEBUG // Integrity check.
                assert(size > 0);
                #endif

                this->_entries.resize(size, scalar);
            }

            /**
             * @brief Scalar copy operator.
             * 
             * @param scalar Scalar.
             * @return Vector& 
             */
            Vector &operator =(const T &scalar) {
                this->_entries.resize(this->_size, scalar);
                return *this;
            }

            /**
             * @brief Standard vector constructor.
             * 
             * @param vector Standard vector.
             */
            Vector(const std::vector<T> &vector): _size(vector.size()) {
                #ifndef NDEBUG // Integrity check.
                assert(vector.size() > 0);
                #endif

                this->_entries.resize(vector.size(), static_cast<T>(0));
                std::copy(vector.begin(), vector.end(), this->_entries.begin());
            }

            /**
             * @brief Standard vector copy operator.
             * 
             * @param vector Standard vector.
             * @return Vector& 
             */
            Vector &operator =(const std::vector<T> &vector) {
                #ifndef NDEBUG // Integrity check.
                assert(this->_size == vector.size());
                #endif

                std::copy(vector.begin(), vector.end(), this->_entries.begin());
                return *this;
            }

            /**
             * @brief Copy constructor.
             * 
             * @param vector Vector.
             */
            Vector(const Vector &vector): _size{vector._size} {
                this->_entries.resize(vector._size, static_cast<T>(0));
                std::copy(vector._entries.begin(), vector._entries.end(), this->_entries.begin());
            }
            
            /**
             * @brief Copy operator.
             * 
             * @param vector Vector.
             * @return Vector& 
             */
            Vector &operator =(const Vector &vector) {
                #ifndef NDEBUG // Integrity check.
                assert(this->_size == vector._size);
                #endif

                std::copy(vector._entries.begin(), vector._entries.end(), this->_entries.begin());
                return *this;
            }

            // Subscript operator, legacy scalar access (C++23).

            #if __cplusplus > 202002L

            /**
             * @brief Scalar reference access, legacy.
             * 
             * @param j Index.
             * @return T& 
             */
            T &operator [](const Natural &j) {
                #ifndef NDEBUG // Integrity check.
                assert(j < this->_size);
                #endif

                return this->_entries[j];
            }

            #endif

            // Call operator, subscript behaviour.

            /**
             * @brief Scalar access.
             * 
             * @param j 
             * @return T 
             */
            T operator ()(const Natural &j) const {
                #ifndef NDEBUG // Integrity check.
                assert(j < this->_size);
                #endif

                return this->_entries[j];
            }

            /**
             * @brief Mask access.
             * 
             * @param M Mask.
             * @return Vector 
             */
            Vector operator ()(const Mask &M) const {
                #ifndef NDEBUG // Integrity check.
                assert(this->_size == M._size);
                #endif

                std::vector<T> vector;

                for(Natural j = 0; j < this->_size; ++j)
                    if(M(j))
                        vector.emplace_back(this->_entries[j]);

                return vector;
            }

            /**
             * @brief Vectorial access.
             * 
             * @param J Indices.
             * @return Vector 
             */
            Vector operator ()(const std::vector<Natural> &J) const {
                #ifndef NDEBUG // Integrity check.
                for(const auto &j: J)
                    assert(j < this->_size);
                #endif

                std::vector<T> vector;
                for(const auto &j: J)
                    vector.emplace_back(this->_entries[j]);

                return vector;
            }

            /**
             * @brief Scalar insert.
             * 
             * @param j Index.
             * @param scalar 
             */
            void operator ()(const Natural &j, const T &scalar) {
                #ifndef NDEBUG // Integrity check.
                assert(j < this->_size);
                #endif

                this->_entries[j] = scalar;
            }

            /**
             * @brief Vectorial insert.
             * 
             * @param J Indices.
             * @param vector 
             */
            void operator ()(const std::vector<Natural> &J, const Vector &vector) {
                #ifndef NDEBUG // Integrity check.
                assert(J.size() == vector._size);
                for(const auto &j: J)
                    assert(j < this->_size);
                #endif

                for(const auto &j: J)
                    this->_entries[j] = vector._entries[j];
            }

            // Comparisons.

            /**
             * @brief Vector == scalar.
             * 
             * @param scalar Scalar.
             * @return Mask 
             */
            Mask operator ==(const T &scalar) const {
                Mask comparison{this->_size};
                std::transform(this->_entries.begin(), this->_entries.end(), comparison._entries.begin(), [scalar](const T &entry){ return std::abs(entry - scalar) <= NUMERICAL_ZERO; });
                return comparison;
            }

            /**
             * @brief Vector != scalar.
             * 
             * @param scalar Scalar.
             * @return Mask 
             */
            Mask operator !=(const T &scalar) const {
                return -(*this == scalar);
            }

            /**
             * @brief Vector < scalar.
             * 
             * @param scalar Scalar.
             * @return Mask 
             */
            Mask operator <(const T &scalar) const {
                Mask comparison{this->_size};
                std::transform(this->_entries.begin(), this->_entries.end(), comparison._entries.begin(), [scalar](const T &entry){ return entry < scalar; });
                return comparison;
            }

            /**
             * @brief Vector <= scalar.
             * 
             * @param scalar Scalar.
             * @return Mask 
             */
            Mask operator <=(const T &scalar) const {
                return (*this < scalar) + (*this == scalar);
            }

            /**
             * @brief Vector > scalar.
             * 
             * @param scalar 
             * @return Mask 
             */
            Mask operator >(const T &scalar) const {
                return -(*this <= scalar);
            }

            /**
             * @brief Vector >= scalar.
             * 
             * @param scalar 
             * @return Mask 
             */
            Mask operator >=(const T &scalar) const {
                return -(*this < scalar);
            }

            /**
             * @brief Vector == vector.
             * 
             * @param vector Vector.
             * @return Mask 
             */
            Mask operator ==(const Vector &vector) const {
                #ifndef NDEBUG // Integrity check.
                assert(this->_size == vector._size);
                #endif

                Mask comparison{this->_size};
                std::transform(this->_entries.begin(), this->_entries.end(), vector._entries.begin(), comparison._entries.begin(), [](const T &t_entry, const T &v_entry){ return std::abs(t_entry - v_entry) <= NUMERICAL_ZERO; });
                return comparison;
            }

            /**
             * @brief Vector != vector.
             * 
             * @param vector Vector.
             * @return Mask 
             */
            Mask operator !=(const Vector &vector) const {
                return -(*this == vector);
            }

            /**
             * @brief Vector < vector.
             * 
             * @param vector Vector.
             * @return Mask 
             */
            Mask operator <(const Vector &vector) const {
                #ifndef NDEBUG // Integrity check.
                assert(this->_size == vector._size);
                #endif

                Mask comparison{this->_size};
                std::transform(this->_entries.begin(), this->_entries.end(), vector._entries.begin(), comparison._entries.begin(), [](const T &t_entry, const T &v_entry){ return t_entry < v_entry; });
                return comparison;
            }

            /**
             * @brief Vector <= vector.
             * 
             * @param vector Vector.
             * @return Mask 
             */
            Mask operator <=(const Vector &vector) const {
                return (*this < vector) + (*this == vector);
            }

            /**
             * @brief Vector > vector.
             * 
             * @param vector Vector.
             * @return Mask 
             */
            Mask operator >(const Vector &vector) const {
                return -(*this <= vector);
            }

            /**
             * @brief Vector >= vector.
             * 
             * @param vector Vector.
             * @return Mask 
             */
            Mask operator >=(const Vector &vector) const {
                return -(*this < vector);
            }

            // Operations.

            /**
             * @brief +Vector.
             * 
             * @return Vector 
             */
            Vector operator +() const { return *this; }

            /**
             * @brief -Vector.
             * 
             * @return Vector
             */
            Vector operator -() const {
                Vector result{this->_size};
                std::transform(this->_entries.begin(), this->_entries.end(), result._entries.begin(), [](const T &entry){ return -entry; });
                return result;
            }

            /**
             * @brief Vector + scalar.
             * 
             * @param scalar Scalar.
             * @return Vector 
             */
            Vector operator +(const T &scalar) const {
                Vector result{this->_size};
                std::transform(this->_entries.begin(), this->_entries.end(), result._entries.begin(), [scalar](const T &entry){ return entry + scalar; });
                return result;
            }

            /**
             * @brief Scalar + vector.
             * 
             * @param scalar Scalar.
             * @param vector Vector.
             * @return Vector 
             */
            friend Vector operator +(const T &scalar, const Vector &vector) {
                Vector result{vector._size};
                std::transform(vector._entries.begin(), vector._entries.end(), result._entries.begin(), [scalar](const T &entry){ return scalar + entry; });
                return result;
            }

            /**
             * @brief Vector += scalar.
             * 
             * @param scalar Scalar.
             * @return Vector& 
             */
            Vector &operator +=(const T &scalar) {
                std::transform(this->_entries.begin(), this->_entries.end(), this->_entries.begin(), [scalar](const T &entry){ return entry + scalar; });
                return *this;
            }

            /**
             * @brief Vector - scalar.
             * 
             * @param scalar Scalar.
             * @return Vector 
             */
            Vector operator -(const T &scalar) const {
                Vector result{this->_size};
                std::transform(this->_entries.begin(), this->_entries.end(), result._entries.begin(), [scalar](const T &entry){ return entry - scalar; });
                return result;
            }

            /**
             * @brief Scalar - vector.
             * 
             * @param scalar Scalar.
             * @param vector Vector.
             * @return Vector 
             */
            friend Vector operator -(const T &scalar, const Vector &vector) {
                Vector result{vector._size};
                std::transform(vector._entries.begin(), vector._entries.end(), result._entries.begin(), [scalar](const T &entry){ return scalar - entry; });
                return result;
            }

            /**
             * @brief Vector -= scalar.
             * 
             * @param scalar Scalar.
             * @return Vector& 
             */
            Vector &operator -=(const T &scalar) {
                std::transform(this->_entries.begin(), this->_entries.end(), this->_entries.begin(), [scalar](const T &entry){ return entry - scalar; });
                return *this;
            }

            /**
             * @brief Vector * scalar.
             * 
             * @param scalar Scalar.
             * @return Vector 
             */
            Vector operator *(const T &scalar) const {
                Vector result{this->_size};
                std::transform(this->_entries.begin(), this->_entries.end(), result._entries.begin(), [scalar](const T &entry){ return entry * scalar; });
                return result;
            }

            /**
             * @brief Scalar * vector.
             * 
             * @param scalar Scalar.
             * @param vector Vector.
             * @return Vector 
             */
            friend Vector operator *(const T &scalar, const Vector &vector) {
                Vector result{vector._size};
                std::transform(vector._entries.begin(), vector._entries.end(), result._entries.begin(), [scalar](const T &entry){ return scalar * entry; });
                return result;
            }

            /**
             * @brief Vector *= scalar.
             * 
             * @param scalar Scalar.
             * @return Vector& 
             */
            Vector &operator *=(const T &scalar) {
                std::transform(this->_entries.begin(), this->_entries.end(), this->_entries.begin(), [scalar](const T &entry){ return entry * scalar; });
                return *this;
            }

            /**
             * @brief Vector / scalar.
             * 
             * @param scalar Scalar.
             * @return Vector 
             */
            Vector operator /(const T &scalar) const {
                Vector result{this->_size};
                std::transform(this->_entries.begin(), this->_entries.end(), result._entries.begin(), [scalar](const T &entry){ return entry / scalar; });
                return result;
            }

            /**
             * @brief Scalar / vector.
             * 
             * @param scalar Scalar.
             * @param vector Vector.
             * @return Vector 
             */
            friend Vector operator /(const T &scalar, const Vector &vector) {
                Vector result{vector._size};
                std::transform(vector._entries.begin(), vector._entries.end(), result._entries.begin(), [scalar](const T &entry){ return scalar / entry; });
                return result;
            }

            /**
             * @brief Vector /= scalar.
             * 
             * @param scalar Scalar.
             * @return Vector& 
             */
            Vector &operator /=(const T &scalar) {
                std::transform(this->_entries.begin(), this->_entries.end(), this->_entries.begin(), [scalar](const T &entry){ return entry / scalar; });
                return *this;
            }

            /**
             * @brief Vector + vector.
             * 
             * @param vector Vector.
             * @return Vector 
             */
            Vector operator +(const Vector &vector) const {
                #ifndef NDEBUG // Integrity check.
                assert(this->_size == vector._size);
                #endif

                Vector result{this->_size};
                std::transform(this->_entries.begin(), this->_entries.end(), vector._entries.begin(), result._entries.begin(), [](const T &t_entry, const T &v_entry){ return t_entry + v_entry; });
                return result;
            }

            /**
             * @brief Vector += vector.
             * 
             * @param vector 
             * @return Vector& 
             */
            Vector &operator +=(const Vector &vector) {
                #ifndef NDEBUG // Integrity check.
                assert(this->_size == vector._size);
                #endif

                std::transform(this->_entries.begin(), this->_entries.end(), vector._entries.begin(), this->_entries.begin(), [](T &t_entry, const T &v_entry){ return t_entry + v_entry; });
                return *this;
            }

            /**
             * @brief Vector - vector.
             * 
             * @param vector Vector.
             * @return Vector 
             */
            Vector operator -(const Vector &vector) const {
                #ifndef NDEBUG // Integrity check.
                assert(this->_size == vector._size);
                #endif

                Vector result{this->_size};
                std::transform(this->_entries.begin(), this->_entries.end(), vector._entries.begin(), result._entries.begin(), [](const T &t_entry, const T &v_entry){ return t_entry - v_entry; });
                return result;
            }

            /**
             * @brief Vector -= vector.
             * 
             * @param vector 
             * @return Vector& 
             */
            Vector &operator -=(const Vector &vector) {
                #ifndef NDEBUG // Integrity check.
                assert(this->_size == vector._size);
                #endif

                std::transform(this->_entries.begin(), this->_entries.end(), vector._entries.begin(), this->_entries.begin(), [](T &t_entry, const T &v_entry){ return t_entry - v_entry; });
                return *this;
            }

            /**
             * @brief Vector * vector (element-wise).
             * 
             * @param vector Vector.
             * @return Vector 
             */
            Vector operator *(const Vector &vector) const {
                #ifndef NDEBUG // Integrity check.
                assert(this->_size == vector._size);
                #endif

                Vector result{this->_size};
                std::transform(this->_entries.begin(), this->_entries.end(), vector._entries.begin(), result._entries.begin(), [](const T &t_entry, const T &v_entry){ return t_entry * v_entry; });
                return result;
            }

            /**
             * @brief Vector *= vector (element-wise).
             * 
             * @param vector 
             * @return Vector& 
             */
            Vector &operator *=(const Vector &vector) {
                #ifndef NDEBUG // Integrity check.
                assert(this->_size == vector._size);
                #endif

                std::transform(this->_entries.begin(), this->_entries.end(), vector._entries.begin(), this->_entries.begin(), [](T &t_entry, const T &v_entry){ return t_entry * v_entry; });
                return *this;
            }

            /**
             * @brief Vector / vector (element-wise).
             * 
             * @param vector Vector.
             * @return Vector 
             */
            Vector operator /(const Vector &vector) const {
                #ifndef NDEBUG // Integrity check.
                assert(this->_size == vector._size);
                #endif

                Vector result{this->_size};
                std::transform(this->_entries.begin(), this->_entries.end(), vector._entries.begin(), result._entries.begin(), [](const T &t_entry, const T &v_entry){ return t_entry / v_entry; });
                return result;
            }

            /**
             * @brief Vector /= vector (element-wise).
             * 
             * @param vector 
             * @return Vector& 
             */
            Vector &operator /=(const Vector &vector) {
                #ifndef NDEBUG // Integrity check.
                assert(this->_size == vector._size);
                #endif

                std::transform(this->_entries.begin(), this->_entries.end(), vector._entries.begin(), this->_entries.begin(), [](T &t_entry, const T &v_entry){ return t_entry / v_entry; });
                return *this;
            }

            // Output.

            /**
             * @brief Vector output.
             * 
             * @param ost 
             * @param vector Vector.
             * @return std::ostream& 
             */
            friend std::ostream &operator <<(std::ostream &ost, const Vector &vector) {
                ost << "(";

                for(Natural j = 0; j < vector._size - 1; ++j)
                    ost << vector(j) << ", ";

                return ost << vector(vector._size - 1) << ")" << std::flush;
            }
    };

}

#endif