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

# Project structure

- `data` contains an Excel spreadsheet with the charts and raw data of the MSSs.
- `extern` contains [prettyprint](https://github.com/Anmol-Singh-Jaggi/Pretty-print).
- `include/` contains the framework.
- `src` contains the programs.
