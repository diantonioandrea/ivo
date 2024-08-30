/**
 * @file Line21_2.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Line21 related methods. Space only.
 * @date 2024-07-23
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef GEOMETRY21_METHODS_LINE21_2
#define GEOMETRY21_METHODS_LINE21_2

#include "../Includes.hpp"
#include "../Line21.hpp"

namespace ivo {

    // Line methods.

    Line21 bisector2(const Edge21 &);
    Line21 bisector2(const Point21 &, const Point21 &);

}

#endif