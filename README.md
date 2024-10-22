# Ivo

_A `C++23` [DGFE](https://en.wikipedia.org/wiki/Discontinuous_Galerkin_method) library for `2+1` space-time problems_

Named after [**Ivo Babuška**](https://en.wikipedia.org/wiki/Ivo_Babuška).

## Table of Contents

- [Table of Contents](#table-of-contents)
- [Overview](#overview)
    - [Key Components](#key-components)
    - [Key Features](#key-features)
- [Setup](#setup)
    - [Cloning the Repository](#cloning-the-repository)

## Overview

### Key Components

- `include/`: Holds declarations and definitions for the classes and methods utilized in the library.
    - `include/Ivo/Base`: Holds declarations and definitions for the base types, concepts, and classes.
    - `include/Ivo/Algebra`: Holds definitions for the algebra classes and methods.
    - `include/Ivo/Geometry21`: Holds declarations for the `2+1` geometry classes and methods.
    - `include/Ivo/Mesh21`: Holds declarations for the `2+1` mesh.
    - `include/Ivo/Fem`: Holds declarations for the multi-purpose **FEM** tools.
    - `include/Ivo/Problem`: Holds declarations for the `2+1` problem.
- `src/`: Holds definitions for the structures and methods utilized in the library, along with custom tests and scripts.

### Key Features

- **Algebra**
    - _Support for **dense** vectors and matrices_
    - _Support for **sparse** matrices and linear systems_
- **Geometry**
    - _Support for `2+1` points, edges, lines, and polygons_
    - _Generation of `1` and `2` Mesh diagrams_
    - _Support for `2+1` meshes_
- **DGFEM**
    - _Support for Legendre polynomials_
    - _Support for Gauss-Legendre quadrature_
    - _Assembly of the `2+1` **DGFE** problem_
    - _Solution of the `2+1` **DGFE** problem_

## Setup

### Cloning the Repository

Clone the repository from [here](https://github.com/diantonioandrea/ivo):

```bash
git clone git@github.com:diantonioandrea/ivo.git
```
