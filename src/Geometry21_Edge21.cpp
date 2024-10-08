/**
 * @file Geometry21_Edge21.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-07-22
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Ivo.hpp>

namespace ivo {

    // Constructors and copy operators.

    /**
     * @brief Points constructor.
     * 
     * @param a Edge's a.
     * @param b Edge's b.
     */
    Edge21::Edge21(const Point21 &a, const Point21 &b): _a{a}, _b{b} {
        #ifndef NDEBUG // Integrity check.
        assert(a != b);
        #endif
    }

    /**
     * @brief Copy constructor.
     * 
     * @param edge Edge.
     */
    Edge21::Edge21(const Edge21 &edge): _a{edge._a}, _b{edge._b} {}

    /**
     * @brief Copy operator.
     * 
     * @param edge Edge.
     * @return Edge21& 
     */
    Edge21 &Edge21::operator =(const Edge21 &edge) {
        this->_a = edge._a;
        this->_b = edge._b;
        return *this;
    }

    // Comparison.

    /**
     * @brief Edge == edge.
     * 
     * @param edge Edge.
     * @return true 
     * @return false 
     */
    bool Edge21::operator ==(const Edge21 &edge) const {
        return ((this->_a == edge._a) && (this->_b == edge._b)) || ((this->_a == edge._b) && (this->_b == edge._a));
    }

    /**
     * @brief Edge != edge.
     * 
     * @param edge Edge.
     * @return true 
     * @return false 
     */
    bool Edge21::operator !=(const Edge21 &edge) const {
        return ((this->_a != edge._a) || (this->_b != edge._b)) && ((this->_a != edge._b) || (this->_b != edge._a));
    }
    
    // Access.

    /**
     * @brief Const scalar access.
     * 
     * @param j Index.
     * @return Point21& 
     */
    Point21 Edge21::operator ()(const Natural &j) const {
        #ifndef NDEBUG // Integrity check.
        assert(j <= 1);
        #endif

        if(j == 0)
            return this->_a;
        
        return this->_b;
    }

    /**
     * @brief Scalar access.
     * 
     * @param j Index.
     * @return Point21& 
     */
    Point21 &Edge21::operator [](const Natural &j) {
        #ifndef NDEBUG // Integrity check.
        assert(j <= 1);
        #endif

        if(j == 0)
            return this->_a;
        
        return this->_b;
    }

    // Insert.

    /**
     * @brief Scalar insert.
     * 
     * @param j Index.
     * @param point Point.
     */
    void Edge21::operator ()(const Natural &j, const Point21 &point) {
        #ifndef NDEBUG // Integrity check.
        assert(j <= 1);
        #endif

        if(j == 0)
            this->_a = point;
        else
            this->_b = point;

        #ifndef NDEBUG // Integrity check.
        assert(this->_a != this->_b);
        #endif
    }

    // Output.
    
    /**
     * @brief Edge output.
     * 
     * @param ost 
     * @param edge Edge.
     * @return std::ostream& 
     */
    std::ostream &operator <<(std::ostream &ost, const Edge21 &edge) {
        return ost << "[" << edge._a << ", " << edge._b << "]" << std::flush;
    }

}