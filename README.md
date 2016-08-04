# Mikado

Friendly solver interfaces for Trilinos.

[![Build Status](https://travis-ci.org/nschloe/mikado.svg?branch=master)](https://travis-ci.org/nschloe/mikado)
[![codecov](https://codecov.io/gh/nschloe/mikado/branch/master/graph/badge.svg)](https://codecov.io/gh/nschloe/mikado)
[![Coverity Scan](https://img.shields.io/coverity/scan/9037.svg?maxAge=2592000)](https://scan.coverity.com/projects/nschloe-mikado)


Trilinos is powerful, but can be hard to use. Mikado tries to make things
a little easier by providing a simple user interface for various Trilinos
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

Available solvers with examples:

(Note that for older versions of Boost, you might need to wrap the strings in
`std::string()`.)

* [*Amesos2*](https://trilinos.org/packages/amesos2/)
    ```c++
using dict = std::map<std::string, boost::any>;
mikado::linear_solve(
    A, b, x, dict{
      {"package", "Amesos2"},
      {"method", "SuperLu"},
      {"parameters", dict{
        {"IterRefine", "SLU_DOUBLE"},
        {"SymmetricMode", true}
      }}
    }
    );
    ```
    Check out
    [the Amesos2 parameter documentation](https://trilinos.org/docs/dev/packages/amesos2/doc/html/group__amesos2__solver__parameters.html)
    for more details.

* [*Belos*](https://trilinos.org/packages/belos/)
    ```c++
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
      {"preconditioner", "MueLu"},
      {"preconditioner parameters", dict{
        {"cycle type", "W"}
      }}
    }
    );
    ```
    Solver parameters are different from method to method; see, e.g.,
    [here](https://trilinos.org/docs/dev/packages/belos/doc/html/classBelos_1_1PseudoBlockGmresSolMgr.html).

* [*MueLu*](https://trilinos.org/packages/muelu/)
    ```c++
mikado::linear_solve(
    A, b, x, dict{
      {"package", "MueLu"},
      {"preconditioner parameters", dict{
        {"cycle type", "W"}
      }}
      }
    );
    ```
    [The MueLu User Guide](https://trilinos.org/wordpress/wp-content/uploads/2015/05/MueLu_Users_Guide_Trilinos12_0.pdf)
    provides full options documentation.


#### Nonlinear solvers

Given a model of type `Thyra::ModelEvaluatorDefaultBase<double>`, solving a
nonlinear equation system is as easy as
```c++
  const auto model = std::make_shared<your_model>(
    // ...
  );

  const auto sol = mikado::nonlinear_solve(
      model,
      {
        {"method", "Newton"}
      }
      );
```
This uses [NOX](https://trilinos.org/packages/nox-and-loca/); solver options
can be taken from
[the NOX manual](https://trilinos.org/docs/dev/packages/nox/doc/html/parameters.html).

Numerical parameter continuation via LOCA can be done in a similarly
straightforward way. An example with various parameters:
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

_ProTip_

Providing the solver configuration in a YAML file is as easy as
```c++
#include <yaml-cpp/yaml.h>
// [...]
const auto config = YAML::LoadFile("solver-config.yml");

mikado::parameter_continuation(
    model, saver, mikado::yaml_to_dict(config)
    );
```
For the above example, the YAML file would be
```yaml
---
NOX:
  Status Tests:
    Test Type: NormF
    Norm Type: Two Norm
    Tolerance: 1.0e-8
  Printing:
    Output Information:
      Details: true
      Outer Iteration: true
      Outer Iteration Status Test: true
      Inner Iteration: true
      Linear Solver Details: true
      Parameters: true
      Warning: true
      Debug: true
      Test Details: true
      Error: true
      Stepper Iteration: true
      Stepper Details: true
      Stepper Parameters: true
LOCA:
  Predictor:
    Method: Tangent
  Stepper:
    Continuation Method: Arc Length
    Continuation Parameter: alpha
    Initial Value: 1.0
    Min Value: 0.0
    Max Value: 3.0
    Max Nonlinear Iterations: 5
  Step Size:
    Initial Step Size: 1.0e-1
    Min Step Size: 1.0e-5
    Max Step Size: 5.0e-1
    Aggressiveness: 0.1
```

### Installation

#### Ubuntu PPA

If you're using Ubuntu, you can get Mikado from a dedicated PPA at
https://launchpad.net/~nschloe/+archive/ubuntu/mikado-nightly/. Simply
```sh
sudo apt-add-repository ppa:nschloe/mikado-nightly
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

### Development

#### Testing
Compile mikado as usual and run `ctest` in the build directory.

#### Static analyis
clang++ offers excellent very good static code analysis. To run it, configure
mikado to be built with clang++,
```
OMPI_CXX=clang++ \
cmake \
  -DCMAKE_CXX_COMPILER:PATH=/usr/bin/mpicxx \
  ../source/
```
Then run
```
scan-build make
```
in the build directory.

### License

Mikado is published under the [BSD 3-Clause License](https://opensource.org/licenses/BSD-3-Clause).
