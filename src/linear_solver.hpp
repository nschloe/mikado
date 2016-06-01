#ifndef NOSH_LINEARSOLVER_HPP
#define NOSH_LINEARSOLVER_HPP

#include <MueLu_Hierarchy.hpp>
#include <boost/any.hpp>

#include <map>
#include <memory>

using dict = std::map<std::string, boost::any>;
namespace mikado {
  // https://trilinos.org/docs/dev/packages/stratimikos/doc/html/index.html
  // http://stackoverflow.com/a/14425299/353337
  /*
  static
  const std::map<std::string, boost::any> default_linear_solver_params = {
    {"package", std::string("Belos")},
    {"method", std::string("Pseudo Block GMRES")},
    {"parameters", dict{
      {"Convergence Tolerance", 1.0e-10},
      {"Output Frequency", 1},
      {"Output Style", 1},
      {"Verbosity", 33}
    }},
    {"preconditioner", std::string("MueLu")},
    {"preconditioner parameters", dict{
    }}
  };
  */

  static
  const std::map<std::string, boost::any> default_linear_solver_params = {
    {"package", std::string("Amesos2")},
    {"method", std::string("KLU2")},
    {"parameters", dict{}}
  };

  static
  const std::map<std::string, boost::any> default_belos_params = {
    {"method", std::string("Pseudo Block GMRES")},
    {"parameters", dict{}}
  };

  //static
  //const std::map<std::string, boost::any> default_linear_solver_params = {
  //  {"package", std::string("MueLu")},
  //  {"parameters", dict{}}
  //};

  //std::map<std::string, boost::any> default_linear_solver_params =
  //{
  //  {"Linear Solver Type", std::string("Belos")},
  //  {"Linear Solver Types", dict{
  //    {"Belos", dict{
  //      {"Solver Type", std::string("Pseudo Block GMRES")},
  //      {"Solver Types", dict{
  //        {"Pseudo Block GMRES", dict{
  //          {"Convergence Tolerance", 1.0e-10},
  //          {"Output Frequency", 1},
  //          {"Output Style", 1},
  //          {"Verbosity", 33}
  //        }}
  //      }}
  //    }}
  //  }},
  //  {"Preconditioner Type", std::string("None")}
  //};

  void
  linear_solve(
      const Tpetra::CrsMatrix<double,int,int> & A,
      const Tpetra::Vector<double,int,int> & b,
      Tpetra::Vector<double,int,int> & x,
      std::map<std::string, boost::any> solver_params = mikado::default_linear_solver_params
      );

  void
  linear_solve_amesos2(
      const Tpetra::CrsMatrix<double,int,int> & A,
      const Tpetra::Vector<double,int,int> & b,
      Tpetra::Vector<double,int,int> & x,
      std::map<std::string, boost::any> solver_params = mikado::default_linear_solver_params
      );

  void
  linear_solve_belos(
      const Tpetra::Operator<double,int,int> & A,
      const Tpetra::Vector<double,int,int> & b,
      Tpetra::Vector<double,int,int> & x,
      std::map<std::string, boost::any> solver_params = mikado::default_belos_params
      );

  std::shared_ptr<MueLu::Hierarchy<double,int,int>>
  get_muelu_hierarchy(
      const Tpetra::CrsMatrix<double,int,int> & A,
      const std::map<std::string, boost::any> & muelu_params
      );

  void
  linear_solve_muelu(
      const Tpetra::CrsMatrix<double,int,int> & A,
      const Tpetra::Vector<double,int,int> & b,
      Tpetra::Vector<double,int,int> & x,
      std::map<std::string, boost::any> solver_params = mikado::default_linear_solver_params
      );

  std::map<std::string, boost::any>
  convert_to_belos_parameters(
      const std::map<std::string, boost::any> & map
      );
}
#endif // NOSH_LINEARSOLVER_HPP
