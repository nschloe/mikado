#include <catch.hpp>

#include <mikado.hpp>

using dict = std::map<std::string, boost::any>;

// ===========================================================================
TEST_CASE("dummy test", "[dummy]")
{
  const dict & map = {};

  const auto out = mikado::convert_to_belos_parameters(map);

  mikado::show_map(out);
}
// ============================================================================
