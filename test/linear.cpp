#include <catch.hpp>

#include <Teuchos_DefaultComm.hpp>
#include <Tpetra_CrsMatrix.hpp>

#include <mikado.hpp>
#include <linear_solver.hpp>

using dict = std::map<std::string, boost::any>;
// ===========================================================================
// Create and return a simple example CrsMatrix.
Teuchos::RCP<const Tpetra::CrsMatrix<double,int,int>>
create_double_matrix(
    const Teuchos::RCP<const Teuchos::Comm<int>> & comm
    )
{
  const Tpetra::global_size_t numGlobalElements = 100;

  // Construct a Map that puts approximately the same number of
  // equations on each processor.
  const int indexBase = 0;
  auto map = Teuchos::rcp(new Tpetra::Map<int,int>(
        numGlobalElements,
        indexBase,
        comm
        ));

  const size_t numMyElements = map->getNodeNumElements();
  auto myGlobalElements = map->getNodeElementList();

  // Create a Tpetra::Matrix using the Map, with a static allocation
  // dictated by NumNz.  (We know exactly how many elements there will
  // be in each row, so we use static profile for efficiency.)
  auto A = Teuchos::rcp(new Tpetra::CrsMatrix<double,int,int>(map, 1));

  for (size_t i = 0; i < numMyElements; i++) {
    A->insertGlobalValues(
        myGlobalElements[i],
        Teuchos::tuple(myGlobalElements[i]),
        Teuchos::tuple(1.0 / (myGlobalElements[i] + 1))
        );
  }

  // Finish up the matrix.
  A->fillComplete();
  return A;
}
// ===========================================================================
// Create and return a simple example CrsMatrix.
Teuchos::RCP<const Tpetra::CrsMatrix<double,int,long long>>
create_longlong_matrix(
    const Teuchos::RCP<const Teuchos::Comm<int>> & comm
    )
{
  const Tpetra::global_size_t numGlobalElements = 100;

  // Construct a Map that puts approximately the same number of
  // equations on each processor.
  const int indexBase = 0;
  auto map = Teuchos::rcp(new Tpetra::Map<int,long long>(
        numGlobalElements,
        indexBase,
        comm
        ));

  const size_t numMyElements = map->getNodeNumElements();
  auto myGlobalElements = map->getNodeElementList();

  // Create a Tpetra::Matrix using the Map, with a static allocation
  // dictated by NumNz.  (We know exactly how many elements there will
  // be in each row, so we use static profile for efficiency.)
  auto A = Teuchos::rcp(new Tpetra::CrsMatrix<double,int,long long>(map, 1));

  for (size_t i = 0; i < numMyElements; i++) {
    A->insertGlobalValues(
        myGlobalElements[i],
        Teuchos::tuple(myGlobalElements[i]),
        Teuchos::tuple(1.0 / (myGlobalElements[i] + 1))
        );
  }

  // Finish up the matrix.
  A->fillComplete();
  return A;
}
// ===========================================================================
TEST_CASE("default solver", "[default]")
{
  const auto comm = Teuchos::DefaultComm<int>::getComm();

  const auto A = create_double_matrix(comm);

  auto b = Tpetra::Vector<double,int,int>(A->getRangeMap());
  b.putScalar(1.0);

  auto x = Tpetra::Vector<double,int,int>(A->getDomainMap());
  x.putScalar(0.0);

  mikado::linear_solve(*A, b, x);

  const auto x_data = x.getData();
  const auto myGlobalElements = x.getMap()->getNodeElementList();
  REQUIRE(myGlobalElements.size() == x_data.size());
  for (size_t i = 0; i < myGlobalElements.size(); i++) {
    REQUIRE(x_data[i] == Approx(myGlobalElements[i] + 1));
  }
}
// // ===========================================================================
TEST_CASE("default solver (long long)", "[default long long]")
{
  const auto comm = Teuchos::DefaultComm<int>::getComm();

  const auto A = create_longlong_matrix(comm);

  auto b = Tpetra::Vector<double,int,long long>(A->getRangeMap());
  b.putScalar(1.0);

  auto x = Tpetra::Vector<double,int,long long>(A->getDomainMap());
  x.putScalar(0.0);

  mikado::linear_solve(*A, b, x);

  const auto x_data = x.getData();
  const auto myGlobalElements = x.getMap()->getNodeElementList();
  REQUIRE(myGlobalElements.size() == x_data.size());
  for (size_t i = 0; i < myGlobalElements.size(); i++) {
    REQUIRE(x_data[i] == Approx(myGlobalElements[i] + 1));
  }
}
// ===========================================================================
TEST_CASE("Belos solver (no preconditioner)", "[belos no prec]")
{
  const auto comm = Teuchos::DefaultComm<int>::getComm();

  const auto A = create_double_matrix(comm);

  auto b = Tpetra::Vector<double,int,int>(A->getRangeMap());
  b.putScalar(1.0);

  auto x = Tpetra::Vector<double,int,int>(A->getDomainMap());
  x.putScalar(0.0);

  const auto x_data = x.getData();
  const auto myGlobalElements = x.getMap()->getNodeElementList();
  REQUIRE(myGlobalElements.size() == x_data.size());

  // with method
  mikado::linear_solve(
      *A, b, x, dict{
        {"package", std::string("Belos")},
        {"method", std::string("Pseudo Block GMRES")}
      }
      );
  for (size_t i = 0; i < myGlobalElements.size(); i++) {
    REQUIRE(x_data[i] == Approx(myGlobalElements[i] + 1));
  }

  // without method
  REQUIRE_THROWS_AS(
    mikado::linear_solve(
        *A, b, x, dict{
          {"package", std::string("Belos")}
        }
        ),
    std::invalid_argument
  );
}
// ===========================================================================
TEST_CASE("Belos solver with MueLu preconditioner", "[belos muelu]")
{
  const auto comm = Teuchos::DefaultComm<int>::getComm();

  const auto A = create_double_matrix(comm);

  auto b = Tpetra::Vector<double,int,int>(A->getRangeMap());
  b.putScalar(1.0);

  auto x = Tpetra::Vector<double,int,int>(A->getDomainMap());
  x.putScalar(0.0);

  mikado::linear_solve(
      *A, b, x, dict{
        {"package", std::string("Belos")},
        {"method", std::string("Pseudo Block GMRES")},
        {"preconditioner", std::string("MueLu")},
        {"preconditioner parameters", dict{
          {"cycle type", std::string("W")}
        }}
      }
      );

  const auto x_data = x.getData();
  const auto myGlobalElements = x.getMap()->getNodeElementList();
  REQUIRE(myGlobalElements.size() == x_data.size());
  for (size_t i = 0; i < myGlobalElements.size(); i++) {
    REQUIRE(x_data[i] == Approx(myGlobalElements[i] + 1));
  }
}
// ===========================================================================
TEST_CASE("Belos solver with Ifpack2 preconditioner", "[belos ifpack2]")
{
  const auto comm = Teuchos::DefaultComm<int>::getComm();

  const auto A = create_double_matrix(comm);

  auto b = Tpetra::Vector<double,int,int>(A->getRangeMap());
  b.putScalar(1.0);

  auto x = Tpetra::Vector<double,int,int>(A->getDomainMap());
  x.putScalar(0.0);

  mikado::linear_solve(
      *A, b, x, dict{
        {"package", std::string("Belos")},
        {"method", std::string("Pseudo Block GMRES")},
        {"preconditioner", std::string("Ifpack2")}
      }
      );

  const auto x_data = x.getData();
  const auto myGlobalElements = x.getMap()->getNodeElementList();
  REQUIRE(myGlobalElements.size() == x_data.size());
  for (size_t i = 0; i < myGlobalElements.size(); i++) {
    REQUIRE(x_data[i] == Approx(myGlobalElements[i] + 1));
  }
}
// ===========================================================================
TEST_CASE("MueLu solver", "[muelu]")
{
  const auto comm = Teuchos::DefaultComm<int>::getComm();

  const auto A = create_double_matrix(comm);

  auto b = Tpetra::Vector<double,int,int>(A->getRangeMap());
  b.putScalar(1.0);

  auto x = Tpetra::Vector<double,int,int>(A->getDomainMap());
  x.putScalar(0.0);

  const auto x_data = x.getData();
  const auto myGlobalElements = x.getMap()->getNodeElementList();
  REQUIRE(myGlobalElements.size() == x_data.size());

  // Without params
  mikado::linear_solve(
      *A, b, x, dict{
        {"package", std::string("MueLu")}
      }
      );
  for (size_t i = 0; i < myGlobalElements.size(); i++) {
    REQUIRE(x_data[i] == Approx(myGlobalElements[i] + 1));
  }

  // with params
  mikado::linear_solve(
      *A, b, x, dict{
        {"package", std::string("MueLu")},
        {"parameters", dict{
          {"cycle type", std::string("W")},
          {"max cycles", 1},
          {"convergence tolerance", 1.0e-10}
        }}
      }
      );
  for (size_t i = 0; i < myGlobalElements.size(); i++) {
    REQUIRE(x_data[i] == Approx(myGlobalElements[i] + 1));
  }
}
// ===========================================================================
TEST_CASE("MueLu solver (long long)", "[muelu long long]")
{
  const auto comm = Teuchos::DefaultComm<int>::getComm();

  const auto A = create_longlong_matrix(comm);

  auto b = Tpetra::Vector<double,int,long long>(A->getRangeMap());
  b.putScalar(1.0);

  auto x = Tpetra::Vector<double,int,long long>(A->getDomainMap());
  x.putScalar(0.0);

  const auto x_data = x.getData();
  const auto myGlobalElements = x.getMap()->getNodeElementList();
  REQUIRE(myGlobalElements.size() == x_data.size());

  // Without params
  mikado::linear_solve(
      *A, b, x, dict{
        {"package", std::string("MueLu")}
      }
      );
  for (size_t i = 0; i < myGlobalElements.size(); i++) {
    REQUIRE(x_data[i] == Approx(myGlobalElements[i] + 1));
  }
}
// ===========================================================================
TEST_CASE("no solver", "[no solver]")
{
  const auto comm = Teuchos::DefaultComm<int>::getComm();

  const auto A = create_double_matrix(comm);

  auto b = Tpetra::Vector<double,int,int>(A->getRangeMap());
  b.putScalar(1.0);

  auto x = Tpetra::Vector<double,int,int>(A->getDomainMap());
  x.putScalar(0.0);

  REQUIRE_THROWS(
    mikado::linear_solve(
      *A, b, x, dict{}
      );
    );
}
// ============================================================================
TEST_CASE("invalid solver", "[invalid]")
{
  const auto comm = Teuchos::DefaultComm<int>::getComm();

  const auto A = create_double_matrix(comm);

  auto b = Tpetra::Vector<double,int,int>(A->getRangeMap());
  b.putScalar(1.0);

  auto x = Tpetra::Vector<double,int,int>(A->getDomainMap());
  x.putScalar(0.0);

  REQUIRE_THROWS(
    mikado::linear_solve(
      *A, b, x, dict{
        {"package", std::string("invalid")}
      }
      );
    );
}
// ===========================================================================
TEST_CASE("Belos solver with invalid preconditioner", "[belos invalid]")
{
  const auto comm = Teuchos::DefaultComm<int>::getComm();

  const auto A = create_double_matrix(comm);

  auto b = Tpetra::Vector<double,int,int>(A->getRangeMap());
  b.putScalar(1.0);

  auto x = Tpetra::Vector<double,int,int>(A->getDomainMap());
  x.putScalar(0.0);

  REQUIRE_THROWS(
    mikado::linear_solve(
        *A, b, x, dict{
          {"package", std::string("Belos")},
          {"method", std::string("Pseudo Block GMRES")},
          {"preconditioner", std::string("invalid")}
        }
        );
  );
}
// ============================================================================
