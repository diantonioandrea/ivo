/**
 * @file Test_Geometry.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Simple geometry testing.
 * @date 2024-07-22
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Ivo.hpp>

int main(int argc, char **argv) {

    // Points.
    ivo::Point21 a{0, 0, 0};
    ivo::Point21 b{1, 0, 0};
    ivo::Point21 c{1, 1, 0};

    ivo::Point21 d{0, 1, 0};
    ivo::Point21 e{0.75, 0.25, 0};

    // Edges.
    ivo::Edge21 ab{a, b};
    ivo::Edge21 bc{b, c};
    ivo::Edge21 ca{c, a};

    // Lines.
    ivo::Line21 abl{a, b};
    ivo::Line21 bcl{b, c};
    ivo::Line21 cal{c, a};

    ivo::Line21 cdl{c, d};

    // Polygon.
    ivo::Polygon21 abc{{a, b, c}};
    ivo::Polygon21 abcd{{a, b, c, d}};

    // Output.
    std::cout << a << std::endl << std::endl;
    std::cout << b << std::endl << std::endl;
    std::cout << c << std::endl << std::endl;

    std::cout << ab << std::endl << std::endl;
    std::cout << bc << std::endl << std::endl;
    std::cout << ca << std::endl << std::endl;

    std::cout << abl << std::endl << std::endl;
    std::cout << bcl << std::endl << std::endl;
    std::cout << cal << std::endl << std::endl;

    std::cout << abc << std::endl << std::endl;

    std::cout << ivo::distance(a, b) << std::endl << std::endl;
    std::cout << ivo::distance(b, c) << std::endl << std::endl;
    std::cout << ivo::distance(c, a) << std::endl << std::endl;

    std::cout << ivo::distance(abl, c) << std::endl << std::endl;
    std::cout << ivo::distance(bcl, a) << std::endl << std::endl;
    std::cout << ivo::distance(cal, b) << std::endl << std::endl;

    std::cout << ivo::distance(abl, bcl) << std::endl << std::endl;
    std::cout << ivo::distance(abl, cal) << std::endl << std::endl;
    std::cout << ivo::distance(abl, cdl) << std::endl << std::endl;

    std::cout << ivo::contains2(abc, d) << std::endl << std::endl;
    std::cout << ivo::contains2(abc, e) << std::endl << std::endl;

    std::cout << ivo::bisector2(a, b) << std::endl << std::endl;

    std::cout << ivo::reduce2(abc, ivo::bisector2(b, c), e) << std::endl << std::endl;

    std::cout << area(abcd) << std::endl << std::endl;

    std::cout << centroid(abcd) << std::endl << std::endl;

    for(const auto &polygon: ivo::voronoi2(abc, 3))
        std::cout << polygon << std::endl;

    std::cout << std::endl;

    for(const auto &point: ivo::intersections(ivo::Line21{b, d}, abc))
        std::cout << point << std::endl;

    return 0;
}