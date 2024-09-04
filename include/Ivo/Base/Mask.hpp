/**
 * @file Mask.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Masks.
 * @date 2024-07-20
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef BASE_MASK
#define BASE_MASK

#include "./Primitives.hpp"

namespace ivo {

    /**
     * @brief std::vector<bool> wrapper.
     * Needs header implementation.
     * 
     */
    struct Mask {

        // Attributes.

        std::vector<bool> _entries;
        const Natural _size;

        // Constructors and copy operators.

        /**
         * @brief Bool constructor.
         * 
         * @param size Mask's size.
         * @param value Mask's value.
         */
        Mask(const Natural &size, const bool &value = false): _size{size} {
            #ifndef NDEBUG // Integrity check.
            assert(size > 0);
            #endif

            this->_entries.resize(size, value);
        }

        /**
         * @brief Standard vector constructor.
         * 
         * @param vector Standard vector.
         */
        Mask(const std::vector<bool> &vector): _size{vector.size()} {
            #ifndef NDEBUG // Integrity check.
            assert(vector.size() > 0);
            #endif

            this->_entries.resize(vector.size(), false);
            std::copy(vector.begin(), vector.end(), this->_entries.begin());
        }

        // Call operator, subscript behaviour.

        /**
         * @brief Scalar access.
         * 
         * @param j Index.
         * @return true 
         * @return false 
         */
        bool operator ()(const Natural &j) const {
            #ifndef NDEBUG // Integrity check.
            assert(j < this->_size);
            #endif

            return this->_entries[j];
        }

        /**
         * @brief Scalar insert.
         * 
         * @param j Index.
         * @param b Boolean.
         */
        void operator ()(const Natural &j, const bool &b) {
            #ifndef NDEBUG // Integrity check.
            assert(j < this->_size);
            #endif

            this->_entries[j] = b;
        }

        // Operations.

        /**
         * @brief (-) !Mask.
         * 
         * @return Mask 
         */
        Mask operator -() const {
            Mask result{this->_size};
            std::transform(this->_entries.begin(), this->_entries.end(), result._entries.begin(), [](const bool &entry){ return !entry; });
            return result;
        }

        /**
         * @brief Mask && (*) mask.
         * 
         * @param mask 
         * @return Mask 
         */
        Mask operator *(const Mask &mask) const {
            #ifndef NDEBUG // Integrity check.
            assert(this->_size == mask._size);
            #endif

            Mask result{this->_size};
            std::transform(this->_entries.begin(), this->_entries.end(), mask._entries.begin(), result._entries.begin(), [](const bool &t_entry, const bool &m_entry){ return t_entry && m_entry; });
            return result;
        }

        /**
         * @brief Mask || (+) mask.
         * 
         * @param mask 
         * @return Mask 
         */
        Mask operator +(const Mask &mask) const {
            #ifndef NDEBUG // Integrity check.
            assert(this->_size == mask._size);
            #endif

            Mask result{this->_size};
            std::transform(this->_entries.begin(), this->_entries.end(), mask._entries.begin(), result._entries.begin(), [](const bool &t_entry, const bool &m_entry){ return t_entry || m_entry; });
            return result;
        }
    };

}

#endif