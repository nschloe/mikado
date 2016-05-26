#include <catch.hpp>

#include <mikado.hpp>

using dict = std::map<std::string, boost::any>;

// ===========================================================================
TEST_CASE("dummy test", "[dummy]")
{
  const dict & map =
  {
    {"method", "Pseudo Block CG"},
    {"parameters", dict{
      {"Convergence Tolerance", 1.0e-10},
      {"Output Frequency", 1},
      {"Output Style", 1},
      {"Verbosity", 33}
    }},
    {"preconditioner", "MueLu"}
  };

  const auto out = mikado::convert_to_belos_parameters(map);

  mikado::show_map(out);
}
// ============================================================================
