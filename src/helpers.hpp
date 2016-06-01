#ifndef HELPER_HPP
#define HELPER_HPP

#include <boost/any.hpp>
#include <map>

#include <Tpetra_CrsMatrix.hpp>

namespace mikado
{
  void
  show_any(const boost::any & any);

  void
  show_map(
    const std::map<std::string, boost::any> & map,
    const int indent = 0
    );

  void
  show_tpetra_crs_matrix(
    const Tpetra::CrsMatrix<double,int,int> & A
    );

  void
  std_map_to_teuchos_list(
      const std::map<std::string, boost::any> & map,
      Teuchos::ParameterList & p
      );

  std::string
  any_to_string(const boost::any & in);
}
#endif
