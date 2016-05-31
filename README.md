# Mikado

[![Build Status](https://travis-ci.org/nschloe/mikado.svg?branch=master)](https://travis-ci.org/nschloe/mikado)
[![codecov](https://codecov.io/gh/nschloe/mikado/branch/master/graph/badge.svg)](https://codecov.io/gh/nschloe/mikado)

Friendly solver interfaces for Trilinos.

Trilinos is powerful, but notoriously hard to use. Mikado tries to make things
a little more friendly by providing simple user interfaces for various Trilinos
solvers.

Mikado works with the Tpetra stack.

### Usage

#### Linear solvers

Given a Tpetra::CrsMatrix `A` and two vectors `b` and `x`, solving a
linear system is as easy as
```c++
mikado::linear_solve(A, b, x);
```
This uses the default solver (Amesos2 with KLU2).


If you would like to shake things up a little, you could just use
```c++
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
```c++
mikado::linear_solve(
    A, b, x, dict{{"package", "MueLu"}}
    );
```
The packages [Amesos2](https://trilinos.org/packages/amesos2/),
[Belos](https://trilinos.org/packages/belos/), and
[MueLu](https://trilinos.org/packages/muelu/) are supported. Solver options
(such as `"method"` in the above example) can be retrieved from the respective
user manual.

#### Nonlinear solvers

Given a `Thyra::ModelEvaluatorDefaultBase<double>`, solving a nonlinear
equation system is as easy as
```c++
  const auto model = std::make_shared<your_model>(
    // ...
  );

  const auto sol = mikado::nonlinear_solve(
      model,
      {
        {"method", std::string("Newton")}
      }
      );
```
This uses
[NOX](https://trilinos.org/packages/nox-and-loca/) and
again, solver options can be taken from [the NOX
manual](https://trilinos.org/docs/dev/packages/nox/doc/html/parameters.html).

Numerical parameter continuation via LOCA can be done in a similarly
straighforward way. An example with various parameters:
```c++
  const auto model = std::make_shared<your_model>(
    // ...
  );

  // optionally specify your 
  const auto saver = std::make_shared<your_data_saver>(
    // ...
  );

  mikado::parameter_continuation(
      model, saver,
      {
        {"NOX", dict{
          {"Status Tests", dict{
            {"Test Type", "NormF"},
            {"Norm Type", "Two Norm"},
            {"Tolerance", 1.0e-8}
          }},
          {"Printing", dict{
           {"Output Information", dict{
             {"Details", true},
             {"Outer Iteration", true},
             {"Outer Iteration Status Test", true},
             {"Inner Iteration", true},
             {"Linear Solver Details", true},
             {"Parameters", true},
             {"Warning", true},
             {"Debug", true},
             {"Test Details", true},
             {"Error", true},
             {"Stepper Iteration", true},
             {"Stepper Details", true},
             {"Stepper Parameters", true}
           }}
          }}
        }},
        {"LOCA", dict{
          {"Predictor", dict{
            {"Method", "Tangent"}
          }},
          {"Stepper", dict{
            {"Continuation Method", "Arc Length"},
            {"Continuation Parameter", "alpha"},
            {"Initial Value", 1.0},
            {"Min Value", 0.0},
            {"Max Value", 3.0},
            {"Max Nonlinear Iterations", 5},
          }},
          {"Step Size", dict{
            {"Initial Step Size", 1.0e-1},
            {"Min Step Size", 1.0e-5},
            {"Max Step Size", 5.0e-1},
            {"Aggressiveness", 0.1}
          }}
        }}
      }
      );
```

### Installation

#### Ubuntu PPA

If you're using Ubuntu, you can get Mikado from a dedicated PPA at
https://launchpad.net/~nschloe/+archive/ubuntu/mikado-nightly/. Simply
```sh
sudo apt-add-repository ppa:nschloe/fenics-nightly
sudo apt update
sudo apt install libmikado-dev
```

#### Manual installation

Mikado uses CMake for configuration. Make sure to have Trilinos (optionally
with MPI) and Boost and Boost installed on your system. Get Mikado (e.g.,
[from GitHub](https://github.com/nschloe/mikado)) and configure and make as
usual with
```sh
cmake path/to/mikado-source  # possibly with more CMake options
make
make install
```

### License

Mikado is published under the [BSD 3-Clause License](https://opensource.org/licenses/BSD-3-Clause).
