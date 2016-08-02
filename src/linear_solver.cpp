#include "linear_solver.hpp"

#include "helpers.hpp"

#include <Amesos2.hpp>
#include <MueLu_ParameterListInterpreter.hpp>
#include <Stratimikos_DefaultLinearSolverBuilder.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_RCPStdSharedPtrConversions.hpp>
#include <Thyra_Ifpack2PreconditionerFactory.hpp>
#include <Thyra_MueLuPreconditionerFactory.hpp>
#include <Thyra_TpetraThyraWrappers.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Xpetra_CrsMatrixWrap_fwd.hpp>
#include <Xpetra_TpetraCrsMatrix.hpp>

// #include <Amesos2_Details_LinearSolverFactory.hpp>
// #include <Belos_Details_LinearSolverFactory.hpp>
// #include <Ifpack2_Details_LinearSolverFactory.hpp>
// #include <Trilinos_Details_LinearSolver.hpp>
// #include <Trilinos_Details_LinearSolverFactory.hpp>

#include <map>

using SC = double;
using LO = int;
using GO = int;
using MV = Tpetra::MultiVector<SC,LO,GO>;
using OP = Tpetra::CrsMatrix<SC,LO,GO>;
using NormType = MV::mag_type;


// =============================================================================
void
mikado::
linear_solve(
    const Tpetra::CrsMatrix<double,int,int> & A,
    const Tpetra::Vector<double,int,int> & b,
    Tpetra::Vector<double,int,int> & x,
    std::map<std::string, boost::any> solver_params
    )
{
  TEUCHOS_TEST_FOR_EXCEPT_MSG(
      solver_params.find("package") == solver_params.end(),
      "Missing key \"package\" in solver parameters."
      )

  const std::string package = any_to_string(solver_params.at("package"));

  if (package == "Amesos2") {
      linear_solve_amesos2(A, b, x, solver_params);
  } else if (package == "Belos") {
      linear_solve_belos(A, b, x, solver_params);
  } else if (package == "MueLu") {
      linear_solve_muelu(A, b, x, solver_params);
  } else {
      TEUCHOS_TEST_FOR_EXCEPT_MSG(
          true,
          "Unknown linear solver package \"" << package << "\"."
          );
  }
  return;
}
// =============================================================================
void
mikado::
linear_solve_amesos2(
    const Tpetra::CrsMatrix<double,int,int> & A,
    const Tpetra::Vector<double,int,int> & b,
    Tpetra::Vector<double,int,int> & x,
    std::map<std::string, boost::any> solver_params
    )
{
  if (A.getComm()->getRank() == 0) {
    mikado::show_any(solver_params);
    std::cout << std::endl;
  }

  const std::string method = any_to_string(solver_params.at("method"));

  auto solver = Amesos2::create<OP,MV>(
        method,
        Teuchos::rcpFromRef(A),
        Teuchos::rcpFromRef(x),
        Teuchos::rcpFromRef(b)
        );

  // Create appropriate parameter list. Check out
  // <https://trilinos.org/docs/dev/packages/amesos2/doc/html/group__amesos2__solvers.html>.
  std::map<std::string, boost::any> method_params = {
    {method, solver_params.at("parameters")}
  };

  // For valid parameters, see
  // <https://trilinos.org/docs/dev/packages/amesos2/doc/html/group__amesos2__solver__parameters.html>.
  auto p = Teuchos::rcp(new Teuchos::ParameterList());
  std_map_to_teuchos_list(method_params, *p);
  solver->setParameters(p);

  solver->symbolicFactorization().numericFactorization().solve();

  auto out = Teuchos::VerboseObjectBase::getDefaultOStream();
  solver->describe(*out, Teuchos::VERB_EXTREME);

  return;
}
// =============================================================================
void
mikado::
linear_solve_belos(
    const Tpetra::Operator<double,int,int> & A,
    const Tpetra::Vector<double,int,int> & b,
    Tpetra::Vector<double,int,int> & x,
    std::map<std::string, boost::any> solver_params
    )
{
  // set x to 0
  x.putScalar(0.0);

  Stratimikos::DefaultLinearSolverBuilder builder;
  auto p = Teuchos::rcp(new Teuchos::ParameterList());
  std_map_to_teuchos_list(convert_to_belos_parameters(solver_params), *p);
  builder.setParameterList(p);

  auto lowsFactory = builder.createLinearSolveStrategy("");
#ifndef NDEBUG
  TEUCHOS_ASSERT(!lowsFactory.is_null());
#endif
  lowsFactory->setVerbLevel(Teuchos::VERB_LOW);

  const Tpetra::Operator<double,int,int> & opA = A;
  auto thyraA = Thyra::createConstLinearOp(Teuchos::rcpFromRef(opA)); // throws

  Teuchos::RCP<Thyra::LinearOpWithSolveBase<double>> lows;
  if (solver_params.find("preconditioner") == solver_params.end()) {
    // no preconditioner
    lows = Thyra::linearOpWithSolve(
        *lowsFactory,
        thyraA
        );
  } else {
    // handle preconditioner
    const std::string prec_type =
      any_to_string(solver_params.at("preconditioner"));
    Teuchos::RCP<Thyra::PreconditionerFactoryBase<double>> factory;
    // https://github.com/trilinos/Trilinos/issues/535
    // if (prec_type == "Ifpack2") {
    //   factory = Teuchos::rcp(
    //       new Thyra::Ifpack2PreconditionerFactory<Tpetra::CrsMatrix<double,int,int>>()
    //       );
    // } else
    if (prec_type == "MueLu") {
      factory = Teuchos::rcp(
          new Thyra::MueLuPreconditionerFactory<double,int,int>()
          );
    } else {
      TEUCHOS_TEST_FOR_EXCEPT_MSG(
          true,
          "Unknown preconditioner \"" << prec_type << "\". Valid values: \"MueLu\"."
          );
    }

    // TODO isCompatible
    // TEUCHOS_ASSERT(factory->isCompatible(*thyraA));
    if (solver_params.find("preconditioner parameters") != solver_params.end()) {
      const auto & prec_params = boost::any_cast<std::map<std::string,boost::any>>(
          solver_params.at("preconditioner parameters")
          );
      auto prec_p = Teuchos::rcp(new Teuchos::ParameterList());
      std_map_to_teuchos_list(prec_params, *prec_p);
      factory->setParameterList(prec_p);
    }

    const auto prec = factory->createPrec();
    Thyra::initializePrec(*factory, thyraA, prec.ptr());

    lows = lowsFactory->createOp();
    Thyra::initializePreconditionedOp<double>(*lowsFactory, thyraA, prec, lows.ptr());
  }

  const Tpetra::Vector<double,int,int> & vecF = b;
  Tpetra::Vector<double,int,int> & vecX = x;

  auto status = Thyra::solve<double>(
      *lows,
      Thyra::NOTRANS,
      *Thyra::createConstVector(Teuchos::rcpFromRef(vecF)),
      Thyra::createVector(Teuchos::rcpFromRef(vecX)).ptr()
      );

  if (A.getDomainMap()->getComm()->getRank() == 0) {
    std::cout << status << std::endl;
  }
  return;
}
// =============================================================================
std::shared_ptr<MueLu::Hierarchy<double,int,int>>
mikado::
get_muelu_hierarchy(
    const Tpetra::CrsMatrix<double,int,int> & A,
    const std::map<std::string, boost::any> & muelu_params
    )
{
  // Tpetra -> Xpetra
  Teuchos::RCP<const Tpetra::CrsMatrix<double,int,int>> ATpetra =
    Teuchos::rcpFromRef(A);
  // cast away the const from A :(
  auto nonconst_ATpetra =
    Teuchos::rcp_const_cast<Tpetra::CrsMatrix<double,int,int>>(ATpetra);
  auto AXpetra = MueLu::TpetraCrs_To_XpetraMatrix(nonconst_ATpetra);

  auto map = AXpetra->getRowMap();

  auto p = Teuchos::rcp(new Teuchos::ParameterList());
  const auto & params = boost::any_cast<std::map<std::string,boost::any>>(
      muelu_params
      );
  std_map_to_teuchos_list(params, *p);

  auto mueLuFactory =
    MueLu::ParameterListInterpreter<double,int,int>(*p, map->getComm());

  auto H = Teuchos::get_shared_ptr(mueLuFactory.CreateHierarchy());
  H->GetLevel(0)->Set("A", AXpetra);

  //// build null space vector
  //auto nullspace = Xpetra::MultiVectorFactory<double,int,int>::Build(map, 1);
  //nullspace->putScalar(1.0);
  //H->GetLevel(0)->Set("Nullspace", nullspace);

  //// TODO
  //// get the coordinates as multivector
  //RCP<MultiVector> coords = Teuchos::rcp(new Xpetra::EpetraMultiVector(epCoord));
  //H->GetLevel(0)->Set("Coordinates", coords);

  mueLuFactory.SetupHierarchy(*H);

  return H;
}
// =============================================================================
void
mikado::
linear_solve_muelu(
    const Tpetra::CrsMatrix<double,int,int> & A,
    const Tpetra::Vector<double,int,int> & b,
    Tpetra::Vector<double,int,int> & x,
    std::map<std::string, boost::any> solver_params
    )
{
  x.putScalar(0.0);
  // Tpetra -> Xpetra
  auto bXpetra = Xpetra::toXpetra(Teuchos::rcpFromRef(b));
  Tpetra::Vector<double,int,int> & xTpetra = x;
  auto xXpetra = Xpetra::toXpetra(Teuchos::rcpFromRef(xTpetra));

  std::map<std::string, boost::any> params;
  try {
    params = boost::any_cast<std::map<std::string, boost::any>>(
        solver_params.at("parameters")
        );
  } catch (std::out_of_range) {
    params = {};
  }

  // store the two custom keys
  //   "max cycles"
  // and
  //   "convergence tolerance"
  // separately.
  const auto mc_it = params.find("max cycles");
  const int max_cycles =
    mc_it != params.end() ?
    boost::any_cast<int>(params.at("max cycles")) :
    10; // default
  if (mc_it != params.end()) {
    params.erase(mc_it);
  }

  const auto ct_it = params.find("convergence tolerance");
  const double convergence_tolerance =
    ct_it != params.end() ?
    boost::any_cast<double>(params.at("convergence tolerance")) :
    1.0e-10; // default
  if (ct_it != params.end()) {
    params.erase(ct_it);
  }

  auto H = get_muelu_hierarchy(A, params);
  H->IsPreconditioner(false);

  H->Iterate(
    *bXpetra,
    *xXpetra,
    std::make_pair(max_cycles, convergence_tolerance)
    );

  return;
}
// =============================================================================
std::map<std::string, boost::any>
mikado::
convert_to_belos_parameters(
    const std::map<std::string, boost::any> & in_map
    )
{
  std::map<std::string, boost::any> out_map = {};

  if (in_map.find("method") == in_map.end()) {
    return out_map;
  }

  const std::string method = any_to_string(in_map.at("method"));

  out_map.insert({"Linear Solver Type", std::string("Belos")});
  out_map.insert({"Linear Solver Types", dict{
      {"Belos", dict{{"Solver Type", method}}}
      }});

  auto & lst = boost::any_cast<dict&>(out_map.at("Linear Solver Types"));
  auto & belos = boost::any_cast<dict&>(lst.at("Belos"));
  if (in_map.find("parameters") != in_map.end()) {
    belos.insert({
        "Solver Types",
        dict{{method, in_map.at("parameters")}}
        });
  } else {
    // insert default parameters
    belos.insert({
        "Solver Types",
        dict{{method,
          dict{
            {"Convergence Tolerance", 1.0e-10},
            {"Output Frequency", 1},
            {"Output Style", 1},
            {"Verbosity", 33}
          }
          }}
        });
  }

  // Valid values include:
  //    {
  //      "None"
  //      "ML"
  //      "Ifpack"
  //    }
  // Since we'd like to use MueLu, Ifpack2, and so forth, set it to "None" here
  // and handle the preconditioner separately.
  out_map.insert({"Preconditioner Type", std::string("None")});

  return out_map;
}
// =============================================================================
