/**
 * @file Mesh21.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 2 + 1 meshes.
 * @date 2024-07-25
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef MESH21_MESH21
#define MESH21_MESH21

#include "./Includes.hpp"
#include "./Element21.hpp"
#include "./Neighbour21.hpp"

namespace ivo {

    /**
     * @brief Mesh21. 2 + 1 mesh.
     * 
     */
    class Mesh21 {

        private:

            // Attributes.

            /**
             * @brief Space cells.
             * 
             */
            const Natural _space;

            /**
             * @brief Time intervals.
             * 
             */
            const Natural _time;

            /**
             * @brief Mesh' elements.
             * More memory required, faster access.
             * 
             */
            std::vector<Element21> _elements;

            /**
             * @brief Elements' neighbours.
             * 
             */
            std::vector<Neighbour21> _neighbours;

        public:

            // Attributes access.

            /**
             * @brief Space cells.
             * 
             * @return Natural 
             */
            constexpr Natural space() const { return this->_space; }

            /**
             * @brief Time intervals.
             * 
             * @return Natural 
             */
            constexpr Natural time() const { return this->_time; }

            // Constructors.

            Mesh21(const std::vector<Polygon21> &, const std::vector<Real> &);

            // Elements and neighbours access.

            /**
             * @brief Scalar element access.
             * 
             * @param j Element index.
             * @return Element21 
             */
            inline Element21 element(const Natural &j) const {
                #ifndef NDEBUG
                assert(j < this->_space * this->_time);
                #endif

                return this->_elements[j];
            }

            /**
             * @brief Scalar element access.
             * 
             * @param j Space index.
             * @param k Time index.
             * @return Element21 
             */
            inline Element21 element(const Natural &j, const Natural &k) const {
                #ifndef NDEBUG
                assert(j < this->_space);
                assert(k < this->_time);
                #endif

                return this->_elements[j + k * _time];
            }

            /**
             * @brief Scalar neighbour access.
             * 
             * @param j Neighbour index.
             * @return Neighbour21 
             */
            inline Neighbour21 neighbour(const Natural &j) const {
                #ifndef NDEBUG
                assert(j < this->_space * this->_time);
                #endif

                return this->_neighbours[j];
            }

            /**
             * @brief Scalar neighbour access.
             * 
             * @param j Space index.
             * @param k Time index.
             * @return Neighbour21 
             */
            inline Neighbour21 neighbour(const Natural &j, const Natural &k) const {
                #ifndef NDEBUG
                assert(j < this->_space);
                assert(k < this->_time);
                #endif

                return this->_neighbours[j + k * _time];
            }

            // Output.

            friend std::ostream &operator <<(std::ostream &, const Mesh21 &);

    };

}

#endif