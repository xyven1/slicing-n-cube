# Abstract

We re-implement in modern C++ an old algorithm that was used to computationally assert that the minimum number of hyperplanes needed to slice all the edges of the 5-dimensional Boolean hypercube is 5. An edge is sliced by a hyperplane if its two endpoints lie on different sides of the hyperplane. Through algorithmic and programming improvements we reduce the previous running time by a factor of over one thousand. We extend this framework to computationally study other aspects of slicing the n-dimensional Boolean hypercube such as what happens if the normal vector of the hyperplanes contains only small integer weights or generalizing hyperplanes to n-variate polynomials.

# Project structure

- `data` contains an Excel spreadsheet with the charts and raw data of the MSSs.
- `extern` contains [prettyprint](https://github.com/Anmol-Singh-Jaggi/Pretty-print).
- `include` contains the framework.
- `src` contains the programs.

# Build instructions

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DN_CUBE_OUT_DIR=.
cmake --build build
```

# Build environment

- CentOS Linux 7 (Core)
- cmake 3.20.3
- gcc 8.2.0
- cgal 4.11
- gmp 6.1.2
- mpfr 4.1.0
- boost 1.74.0
