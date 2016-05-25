#define CATCH_CONFIG_RUNNER
#include <catch.hpp>

#include <Teuchos_GlobalMPISession.hpp>

int main(int argc, char* argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);
  const int result = Catch::Session().run(argc, argv);
  return result;
}
