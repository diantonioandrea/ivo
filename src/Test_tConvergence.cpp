/**
 * @file Test_hConvergence.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief t convergence test.
 * @date 2024-10-19
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include "./include/Square.hpp"

int main(int argc, char **argv) {

    if(argc != 3) {
        std::cout << "Usage: " << argv[0] << " SPACE_DEGREE [p] TIME_DEGREE [q]." << std::endl;
        return -1;
    }

    // Degrees.
    const ivo::Natural p = static_cast<ivo::Natural>(std::atoi(argv[1]));
    const ivo::Natural q = static_cast<ivo::Natural>(std::atoi(argv[2]));

    assert(p > 0);
    assert(q > 0);

    // Output.
    std::ofstream output{"output/tConvergence_" + std::to_string(p) + "_" + std::to_string(q) + ".e21"};

    #ifndef NVERBOSE
    std::cout << "[Ivo] TEST, t convergence\n" << std::endl;
    #else
    std::cout << "[Ivo] TEST, t convergence" << std::endl;
    #endif

    // Elements.
    const ivo::Natural Ns = 512;
    ivo::Natural Nt = 16;

    // Space.
    const std::vector<ivo::Polygon21> space = ivo::mesher2(ivo::square::abcd, Ns);

    // Equation.
    const ivo::Equation equation{ivo::square::convection, ivo::square::diffusion, ivo::square::reaction};
    const ivo::Initial initial{ivo::square::u0};
    const ivo::Data data{ivo::square::g, ivo::square::gd, ivo::square::gn};

    // Tests.
    for(ivo::Natural j = 0; j < 6; ++j) {

        // Time.
        const std::vector<ivo::Real> time = ivo::mesher1(0.0L, 1.0L, Nt);

        // Mesh.
        const ivo::Mesh21 mesh{space, time, p, q};

        // Matrix and vector.
        const ivo::Sparse<ivo::Real> A = ivo::stiffness(mesh, equation);
        const ivo::Vector<ivo::Real> b = ivo::forcing(mesh, equation, data);

        // Solution.
        const ivo::Vector<ivo::Real> x = ivo::solve(mesh, A, b, initial);

        // Error.
        const ivo::Error error{mesh, x, ivo::square::u, ivo::square::u_xy};

        // Output.
        output << error << "\n" << std::endl;

        #ifndef NVERBOSE
        std::cout << "\n\t[TEST] Completed iteration " << j + 1 << "\n" << std::endl;
        #else
        std::cout << "\t[TEST] Completed iteration " << j + 1 << std::endl;
        #endif

        // Update.
        Nt += 16;
    }

    std::cout << "\t[TEST] Exited" << std::endl;

    return 0;
}