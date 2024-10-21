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

- `include/`: Holds declarations for the classes and methods utilized in the library.
- `src/`: Holds definitions for the structures and methods utilized in the library, along with custom tests and scripts.

### Key Features

- **Algebra**
    - _Support for **dense** vectors and matrices_
    - _Support for **sparse** matrix and linear systems_
- **Geometry**
    - _Methods for points, edges, lines, and polygons_
    - _Mesh generation_
- **Finite Element Methods**
    - _Legendre polynomials_
    - _Gauss-Legendre quadrature methods_
    - _Stiffness matrix and force vector assembly_
    - _Custom solver for the **FEM** problem_

## Setup

### Cloning the Repository

Clone the repository from [here](https://github.com/diantonioandrea/ivo):

```bash
git clone git@github.com:diantonioandrea/ivo.git
```
