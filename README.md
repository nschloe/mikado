# Mikado

[![Build Status](https://travis-ci.org/nschloe/mikado.svg?branch=master)](https://travis-ci.org/nschloe/mikado)
[![codecov](https://codecov.io/gh/nschloe/mikado/branch/master/graph/badge.svg)](https://codecov.io/gh/nschloe/mikado)

Friendly solver interfaces for Trilinos.

Trilinos is powerful, but notoriously hard to use. Mikado tries to make things
a little bit easier by providing a simple user interface for various linear and
nonlinear Trilinos solvers.

For example, given a Tpetra::CrsMatrix `A` and two vector `b` and `x, solving a
linear system is as easy as
```
mikado::linear_solve(A, b, x);
```
This uses the default solver (Amesos2 with KLU2).

If you would like to shake things up a little, just use
```
using dict = std::map<std::string, boost::any>;
mikado::linear_solve(
    A, b, x, dict{
      {"package", "Belos"},
      {"method", "Pseudo Block GMRES"},
      {"parameters", dict{
        {"Convergence Tolerance", 1.0e-10},
        {"Output Frequency", 1},
        {"Output Style", 1},
        {"Verbosity", 33}
      }},
      {"preconditioner", "MueLu"}
    }
    );
```
or
```
mikado::linear_solve(
    A, b, x, dict{{"package", "MueLu"}}
    );
```

### Installation

#### Ubuntu PPA

If you're using Ubuntu, you can get Mikado from a dedicated PPA at
https://launchpad.net/~nschloe/+archive/ubuntu/mikado-nightly/. Simply
```
sudo apt-add-repository ppa:nschloe/fenics-nightly
sudo apt update
sudo apt install libmikado-dev
```

#### Manual installation

Mikado uses CMake for configuration. Make sure to have Trilinos (optionally
with MPI) and Boost and Boost installed on your system. Get Mikado (e.g.,
[from GitHub](https://github.com/nschloe/mikado)) and configure and make as
usual with
```
cmake path/to/mikado-source  # possibly with more CMake options
make
make install
```

### License

Mikado is published under the [BSD 3-Clause License](https://opensource.org/licenses/BSD-3-Clause).
