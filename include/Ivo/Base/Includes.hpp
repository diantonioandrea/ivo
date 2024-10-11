/**
 * @file Includes.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-07-22
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef BASE_INCLUDES
#define BASE_INCLUDES

// Types.
#include <cstddef>

// Containers.
#include <vector>
#include <array>

// Assertions.
#include <cassert>

// Concepts.
#include <concepts>

// Algorithms (transform, ...).
#include <algorithm>

// Math.
#include <cmath>
#include <complex>

// Parallelisation.
#ifdef _OPENMP
#include "omp.h"
#endif

#endif