/**
 * @file Test_Forcing.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Simple forcing testing.
 * @date 2024-09-18
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <fstream>
#include <Ivo.hpp>

std::array<ivo::Real, 2> convection(const ivo::Real &);
ivo::Real diffusion(const ivo::Real &);
ivo::Real reaction(const ivo::Real &);

ivo::Real condition(const ivo::Real &, const ivo::Real &);
ivo::Real source(const ivo::Real &, const ivo::Real &);
ivo::Real dirichlet(const ivo::Real &, const ivo::Real &);
ivo::Real neumann(const ivo::Real &, const ivo::Real &);

int main(int argc, char **argv) {

    // Points.
    ivo::Point21 a{0.0, 0.0};
    ivo::Point21 b{1.0, 0.0};
    ivo::Point21 c{1.0, 1.0};
    ivo::Point21 d{0.0, 1.0};

    // Polygon.
    ivo::Polygon21 abcd = {a, b, c, d};

    // Space diagram.
    std::vector<ivo::Polygon21> space = ivo::mesher2(abcd, 5);
    space = ivo::triangulate(space); // Triangulation.

    // Time "diagram" (intervals).
    std::vector<ivo::Real> time{0.0L, 0.5L, 1.0L};

    // Mesh.
    ivo::Mesh21 mesh{space, time};

    // Mesh output.
    std::ofstream output("output/Test_Forcing.p21");
    output << mesh;

    // Equation.
    ivo::Equation equation{convection, diffusion, reaction};

    // Data.
    ivo::Data data{source, dirichlet, neumann};

    // Initial condition.
    ivo::Initial initial{condition};

    // Forcing vector.
    ivo::Vector<ivo::Real> F = ivo::forcing(mesh, equation, data, initial);

    // Stiffness matrix output.
    std::cout << F << std::endl;

    return 0;
}

std::array<ivo::Real, 2> convection(const ivo::Real &t) {
    return {0.0L, 0.0L};
}

ivo::Real diffusion(const ivo::Real &t) {
    return 1.0L;
}

ivo::Real reaction(const ivo::Real &t) {
    return 0.0L;
}

ivo::Real condition(const ivo::Real &x, const ivo::Real &y) {
    return x + y;
}

ivo::Real source(const ivo::Real &x, const ivo::Real &y) {
    return 1.0L;
}

ivo::Real dirichlet(const ivo::Real &x, const ivo::Real &y) {
    return 0.0L;
}

ivo::Real neumann(const ivo::Real &x, const ivo::Real &y) {
    return 0.0L;
}
