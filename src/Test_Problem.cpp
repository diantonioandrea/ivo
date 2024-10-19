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

ivo::Real exact(const ivo::Real &, const ivo::Real &, const ivo::Real &);

ivo::Real exact_t(const ivo::Real &, const ivo::Real &, const ivo::Real &);
std::array<ivo::Real, 2> exact_xy(const ivo::Real &, const ivo::Real &, const ivo::Real &);
ivo::Real exact_xxyy(const ivo::Real &, const ivo::Real &, const ivo::Real &);

ivo::Real condition(const ivo::Real &, const ivo::Real &);
ivo::Real source(const ivo::Real &, const ivo::Real &, const ivo::Real &);
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
    std::vector<ivo::Real> time = ivo::mesher1(0.0L, 1.0L, 5);

    // Mesh.
    ivo::Mesh21 mesh{space, time, 3, 3};

    // Problem.
    ivo::Equation equation{convection, diffusion, reaction};
    ivo::Data data{source, exact, neumann};
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
    return {0.5L, 0.5L};
}

ivo::Real diffusion(const ivo::Real &t) {
    return 1.0L;
}

ivo::Real reaction(const ivo::Real &t) {
    return 1.0L;
}

ivo::Real condition(const ivo::Real &x, const ivo::Real &y) {
    return exact(x, y, 0.0L);
}

ivo::Real exact(const ivo::Real &x, const ivo::Real &y, const ivo::Real &t) { // Dirichlet.
    return std::sin(x) * std::sin(y) * std::sin(t);
}

ivo::Real exact_t(const ivo::Real &x, const ivo::Real &y, const ivo::Real &t) {
    return std::sin(x) * std::sin(y) * std::cos(t);
}

std::array<ivo::Real, 2> exact_xy(const ivo::Real &x, const ivo::Real &y, const ivo::Real &t) {
    return {std::cos(x) * std::sin(y) * std::sin(t), std::sin(x) * std::cos(y) * std::sin(t)};
}

ivo::Real exact_xxyy(const ivo::Real &x, const ivo::Real &y, const ivo::Real &t) {
    return -2.0L * exact(x, y, t);
}

ivo::Real source(const ivo::Real &x, const ivo::Real &y, const ivo::Real &t) {
    auto [convection_x, convection_y] = convection(t);
    auto [exact_x, exact_y] = exact_xy(x, y, t);

    return exact_t(x, y, t) - diffusion(t) * exact_xxyy(x, y, t) + reaction(t) * exact(x, y, t) + convection_x * exact_x + convection_y * exact_y;
}

ivo::Real neumann(const ivo::Real &x, const ivo::Real &y, const ivo::Real &t) {
    auto [exact_x, exact_y] = exact_xy(x, y, t);

    if(x <= ivo::constants::algebra_zero)
        return -1.0L * diffusion(t) * exact_x;

    if(x >= 1.0L - ivo::constants::algebra_zero)
        return diffusion(t) * exact_x;
    
    if(y <= ivo::constants::algebra_zero)
        return -1.0L * diffusion(t) * exact_y;
    
    return diffusion(t) * exact_y;
}