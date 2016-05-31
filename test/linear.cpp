#include <catch.hpp>

#include <Teuchos_DefaultComm.hpp>
#include <Tpetra_CrsMatrix.hpp>

#include <mikado.hpp>

using dict = std::map<std::string, boost::any>;
// ===========================================================================
// Create and return a simple example CrsMatrix.
Teuchos::RCP<const Tpetra::CrsMatrix<double,int,int>>
create_matrix(
    const Teuchos::RCP<const Teuchos::Comm<int>> & comm
    )
{
  const Tpetra::global_size_t numGlobalElements = 1000;

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
TEST_CASE("default solver", "[default]")
{
  const auto comm = Teuchos::DefaultComm<int>::getComm();

  const auto A = create_matrix(comm);

  mikado::show_tpetra_crs_matrix(*A);

  auto b = Tpetra::Vector<double,int,int>(A->getRangeMap());
  b.putScalar(1.0);

  auto x = Tpetra::Vector<double,int,int>(A->getDomainMap());
  x.putScalar(1.0);

  mikado::linear_solve(*A, b, x);

  const auto x_data = x.getData();
  REQUIRE(x_data[0] == Approx(1.0));
  REQUIRE(x_data[1] == Approx(2.0));
  REQUIRE(x_data[2] == Approx(3.0));
}
// ===========================================================================
TEST_CASE("Belos solver", "[belos]")
{
  const auto comm = Teuchos::DefaultComm<int>::getComm();

  const auto A = create_matrix(comm);

  mikado::show_tpetra_crs_matrix(*A);

  auto b = Tpetra::Vector<double,int,int>(A->getRangeMap());
  b.putScalar(1.0);

  auto x = Tpetra::Vector<double,int,int>(A->getDomainMap());
  x.putScalar(1.0);

  mikado::linear_solve(
      *A, b, x, dict{
        {"package", std::string("Belos")},
        {"method", std::string("Pseudo Block GMRES")},
        {"preconditioner", std::string("MueLu")}
      }
      );

  const auto x_data = x.getData();
  REQUIRE(x_data[0] == Approx(1.0));
  REQUIRE(x_data[1] == Approx(2.0));
  REQUIRE(x_data[2] == Approx(3.0));
}
// ===========================================================================
TEST_CASE("MueLu solver", "[muelu]")
{
  const auto comm = Teuchos::DefaultComm<int>::getComm();

  const auto A = create_matrix(comm);

  mikado::show_tpetra_crs_matrix(*A);

  auto b = Tpetra::Vector<double,int,int>(A->getRangeMap());
  b.putScalar(1.0);

  auto x = Tpetra::Vector<double,int,int>(A->getDomainMap());
  x.putScalar(1.0);

  mikado::linear_solve(
      *A, b, x, dict{
        {"package", std::string("MueLu")}
      }
      );

  const auto x_data = x.getData();
  REQUIRE(x_data[0] == Approx(1.0));
  REQUIRE(x_data[1] == Approx(2.0));
  REQUIRE(x_data[2] == Approx(3.0));
}
// ============================================================================
