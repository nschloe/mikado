#include <catch.hpp>

#include <mikado.hpp>

using dict = std::map<std::string, boost::any>;
// ===========================================================================
TEST_CASE("show map", "[show map]")
{
  dict d = dict{
    {"a", 1},
    {"b", 3.14},
    {"c", true},
    {"d", std::string("string")}
  };

  mikado::show_map(d);
}
// ============================================================================
