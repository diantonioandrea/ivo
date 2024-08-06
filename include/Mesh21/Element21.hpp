/**
 * @file Element21.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 2 + 1 elements.
 * @date 2024-07-25
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef MESH21_ELEMENT21
#define MESH21_ELEMENT21

#include "./Includes.hpp"

namespace ivo {

    /**
     * @brief Element21. 2 + 1 elements.
     * Prisms.
     * 
     */
    class Element21 {

        private:

            // Attributes.

            /**
             * @brief Element's base.
             * 
             */
            const Polygon21 _base;

            /**
             * @brief Element's height.
             * 
             */
            const Real _height;

            /**
             * @brief Element's space degree.
             * 
             */
            const Natural _p;

            /**
             * @brief Element's time degree.
             * 
             */
            const Natural _q;

        public:

            // Attributes access.

            /**
             * @brief Element's height.
             * 
             * @return Real 
             */
            constexpr Real height() const { return this->_height; }

            /**
             * @brief Element's space degree.
             * 
             * @return Natural 
             */
            constexpr Natural p() const { return this->_p; }

            /**
             * @brief Element's time degree.
             * 
             * @return Natural 
             */
            constexpr Natural q() const { return this->_q; }

            // Constructors.

            Element21(const Polygon21 &, const Real &, const Natural &, const Natural &);
            Element21(const Polygon21 &, const Real &);
            Element21(const Element21 &);

            // Methods.

            Natural dofs() const;
            
            Polygon21 b_base() const;
            Polygon21 t_base() const;

            std::vector<Edge21> b_edges() const;
            std::vector<Edge21> t_edges() const;

            std::vector<Polygon21> faces() const;

            // Output.

            friend std::ostream &operator <<(std::ostream &, const Element21 &);
    };

}

#endif