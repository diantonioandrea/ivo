/**
 * @file Test_Problem.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Simple problem testing.
 * @date 2024-09-24
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
ivo::Real source(const ivo::Real &, const ivo::Real &, const ivo::Real &);
ivo::Real dirichlet(const ivo::Real &, const ivo::Real &, const ivo::Real &);
ivo::Real neumann(const ivo::Real &, const ivo::Real &, const ivo::Real &);

int main(int argc, char **argv) {

    // Points.
    ivo::Point21 a{0.0, 0.0};
    ivo::Point21 b{1.0, 0.0};
    ivo::Point21 c{1.0, 1.0};
    ivo::Point21 d{0.0, 1.0};

    // Polygon.
    ivo::Polygon21 abcd = {a, b, c, d};

    // Space diagram.
    std::vector<ivo::Polygon21> space = ivo::mesher2(abcd, 15);

    // Time "diagram" (intervals).
    std::vector<ivo::Real> time{0.0L, 0.2L, 0.4L, 0.6L, 0.8L, 1.0L};

    // Mesh.
    ivo::Mesh21 mesh{space, time, 2, 2};

    // Problem.
    ivo::Equation equation{convection, diffusion, reaction};
    ivo::Data data{source, dirichlet, neumann};
    ivo::Initial initial{condition};

    // Stiffness matrix.
    ivo::Sparse<ivo::Real> A = ivo::stiffness(mesh, equation);

    // Forcing vector.
    ivo::Vector<ivo::Real> F = ivo::forcing(mesh, equation, data);

    // Solution.
    ivo::Vector<ivo::Real> X = ivo::solve(mesh, A, F, initial);

    // Solution output.
    ivo::visual(mesh, X, "output/Test_Problem.s21");

    return 0;
}

std::array<ivo::Real, 2> convection(const ivo::Real &t) {
    return {0.0L, 0.0L};
}

ivo::Real diffusion(const ivo::Real &t) {
    return 0.0L;
}

ivo::Real reaction(const ivo::Real &t) {
    return 1.0L;
}

ivo::Real condition(const ivo::Real &x, const ivo::Real &y) {
    return x * x + y * y;
}

ivo::Real source(const ivo::Real &x, const ivo::Real &y, const ivo::Real &t) {
    return condition(x, y) + t;
}

ivo::Real dirichlet(const ivo::Real &x, const ivo::Real &y, const ivo::Real &t) {
    return condition(x, y) + t;
}

ivo::Real neumann(const ivo::Real &x, const ivo::Real &y, const ivo::Real &t) {
    if(x + y <= 1.0L - ivo::constants::algebra_zero)
        return 0.0L;

    return 2.0L * diffusion(t);
}
