#include <catch.hpp>

#include <mikado.hpp>

using dict = std::map<std::string, boost::any>;
// ===========================================================================
TEST_CASE("maps", "[maps]")
{
  struct product {
    int weight;
    double price;
  };
  product banana;
  banana.weight = 1;
  banana.price = 2.0;

  dict d = dict{
    {"a", 1},
    {"b", 3.14},
    {"c", true},
    {"d", std::string("string")},
    {"e", banana}
  };

  mikado::show_map(d);

  auto p = Teuchos::rcp(new Teuchos::ParameterList());
  REQUIRE_THROWS(
      mikado::std_map_to_teuchos_list(d, *p)
      );
}
// ===========================================================================
// Create and return a simple example CrsMatrix.
Teuchos::RCP<const Tpetra::CrsMatrix<double,int,int>>
create_matrix(
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
TEST_CASE("show matrix", "[show matrix]")
{
  const auto comm = Teuchos::DefaultComm<int>::getComm();
  const auto A = create_matrix(comm);

  mikado::show_tpetra_crs_matrix(*A);
}
// ============================================================================
