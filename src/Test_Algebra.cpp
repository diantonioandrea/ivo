/**
 * @file Test_Algebra.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Simple algebra testing.
 * @date 2024-07-20
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Ivo.hpp>

int main(int argc, char **argv) {

    // Vectors.
    ivo::Vector<ivo::Real> v_0{{1, 2}};
    ivo::Vector<ivo::Real> v_1{{3, 4}};

    // Matrices.
    ivo::Matrix<ivo::Real> m_0{2, 2, {1, 2, 3, 4}};
    ivo::Matrix<ivo::Real> m_1{2, 3, {1, 2, 3, 4, 5, 6}};

    // Sparse matrices.
    ivo::Sparse<ivo::Real> s_0{2, 2};

    s_0(0, 0, 1.0L);
    s_0(0, 1, 2.0L);
    s_0(1, 0, 3.0L);
    s_0(1, 1, 4.0L);

    // Output.
    std::cout << v_0 << std::endl << std::endl;
    std::cout << v_1 << std::endl << std::endl;
    std::cout << m_0 << std::endl << std::endl;
    std::cout << m_1 << std::endl << std::endl;
    std::cout << s_0 << std::endl << std::endl;

    // Operations.
    std::cout << v_0 * v_0 << std::endl << std::endl;
    std::cout << v_0 * v_1 << std::endl << std::endl;
    std::cout << m_0 * v_0 << std::endl << std::endl;
    std::cout << m_0 * v_1 << std::endl << std::endl;
    std::cout << s_0 * v_0 << std::endl << std::endl;
    std::cout << s_0 * v_1 << std::endl << std::endl;
    std::cout << v_0 * m_0 << std::endl << std::endl;
    std::cout << v_1 * m_0 << std::endl << std::endl;
    std::cout << v_0 * s_0 << std::endl << std::endl;
    std::cout << v_1 * s_0 << std::endl << std::endl;
    std::cout << m_0 * m_0 << std::endl << std::endl;
    std::cout << m_0 * m_1 << std::endl;

    return 0;
}