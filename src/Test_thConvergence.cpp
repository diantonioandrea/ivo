/**
 * @file Test_hConvergence.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief th convergence test.
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
    std::ofstream output{"output/thConvergence_" + std::to_string(p) + "_" + std::to_string(q) + ".e21"};

    #ifndef NVERBOSE
    std::cout << "[Ivo] TEST, Testing space-time convergence\n" << std::endl;
    #else
    std::cout << "[Ivo] TEST, Testing space-time convergence" << std::endl;
    #endif

    // Space diagrams.
    std::vector<std::string> diagrams;

    diagrams.emplace_back("data/square/Square_128.p2");
    diagrams.emplace_back("data/square/Square_192.p2");
    diagrams.emplace_back("data/square/Square_256.p2");
    diagrams.emplace_back("data/square/Square_384.p2");
    diagrams.emplace_back("data/square/Square_512.p2");
    diagrams.emplace_back("data/square/Square_768.p2");
    diagrams.emplace_back("data/square/Square_1024.p2");
    diagrams.emplace_back("data/square/Square_1536.p2");
    diagrams.emplace_back("data/square/Square_2048.p2");
    diagrams.emplace_back("data/square/Square_3072.p2");

    // Equation.
    const ivo::Equation equation{ivo::square::convection, ivo::square::diffusion, ivo::square::reaction};
    const ivo::Initial initial{ivo::square::u0};
    const ivo::Data data{ivo::square::g, ivo::square::gd, ivo::square::gn};

    // Time coefficient.
    const ivo::Real Ct = std::sqrt(8.0L / (3.0L * std::sqrt(3.0L)));

    // Tests.
    const ivo::Natural tests = diagrams.size();

    // Main loop.
    for(ivo::Natural j = 0; j < tests; ++j) {

        // Timer.
        auto start = std::chrono::high_resolution_clock::now();

        // Space.
        const std::vector<ivo::Polygon21> space = ivo::mesher2(diagrams[j]);

        // Time elements, empirical.
        const ivo::Natural Nt = std::sqrt(static_cast<ivo::Real>(space.size())) / Ct;

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
        const ivo::Error error{mesh, equation, x, ivo::square::u, ivo::square::u_xy};

        // Elapsed time.
        const auto timer = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start);
        const ivo::Real elapsed = timer.count() / 1.0E3L;

        // Output.
        output << error << "\n" << std::endl;

        #ifndef NVERBOSE
        std::cout << "\n\t[TEST] Progress: " << j + 1 << "/" << tests << ", " << elapsed << "s\n" << std::endl;
        #else
        std::cout << "\t[TEST] Progress: " << j + 1 << "/" << tests << ", " << elapsed << "s" << std::endl;
        #endif
    }

    std::cout << "\t[TEST] Exited" << std::endl;

    return 0;
}