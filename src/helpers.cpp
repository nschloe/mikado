#include "helpers.hpp"

#include <iostream>
#include <map>
#include <string>
#include <yaml-cpp/yaml.h>

#include <Teuchos_VerboseObject.hpp>

using dict = std::map<std::string, boost::any>;

namespace mikado
{
  void
  show_any(const boost::any & any)
  {
    if (any.type() == typeid(int)) {
      std::cout << boost::any_cast<int>(any);
    } else if (any.type() == typeid(double)) {
      std::cout << boost::any_cast<double>(any);
    } else if (any.type() == typeid(bool)) {
      std::cout << boost::any_cast<bool>(any);
    } else if (any.type() == typeid(std::string)) {
      std::cout << boost::any_cast<std::string>(any);
    } else if (any.type() == typeid(const char*)) {
      std::cout << boost::any_cast<const char*>(any);
    } else if (any.type() == typeid(std::map<std::string, boost::any>)) {
      show_map(
          boost::any_cast<std::map<std::string, boost::any>>(any),
          2
          );
    } else {
      std::cout << "[unhandled type]";
    }
  }

  void
  show_map(
    const std::map<std::string, boost::any> & map,
    const int indent
    )
  {
    std::cout << "{\n";
    for (const auto &p: map) {
      std::cout << std::string(indent + 2, ' ') << p.first << ": ";
      show_any(p.second);
      std::cout << std::endl;
    }
    std::cout << std::string(indent, ' ') << "}";
  }

  void
  show_tpetra_crs_matrix(
    const Tpetra::CrsMatrix<double,int,int> & A
    )
  {
    auto out = Teuchos::VerboseObjectBase::getDefaultOStream();
    A.describe(*out);
  }

  void
  std_map_to_teuchos_list(
      const std::map<std::string, boost::any> & map,
      Teuchos::ParameterList & p
      )
  {
    for (const auto & entry: map) {
        if(entry.second.type() == typeid(int)) {
          p.set(entry.first, boost::any_cast<int>(entry.second));
        } else if(entry.second.type() == typeid(double)) {
          p.set(entry.first, boost::any_cast<double>(entry.second));
        } else if(entry.second.type() == typeid(bool)) {
          p.set(entry.first, boost::any_cast<bool>(entry.second));
        } else if(entry.second.type() == typeid(const char*)) {
          p.set(entry.first, boost::any_cast<const char*>(entry.second));
        } else if(entry.second.type() == typeid(std::string)) {
          p.set(entry.first, boost::any_cast<std::string>(entry.second));
        } else if(entry.second.type() == typeid(std::map<std::string, boost::any>)) {
          std_map_to_teuchos_list(
              boost::any_cast<std::map<std::string, boost::any>>(entry.second),
              p.sublist(entry.first)
              );
        } else {
          TEUCHOS_TEST_FOR_EXCEPT_MSG(
              true,
              "Unknown value type of key \"" << entry.first << "\"."
              );
        }
    }
    return;
  }

  std::string
  any_to_string(const boost::any & in) {
    try {
      return boost::any_cast<std::string>(in);
    }
    catch (boost::bad_any_cast) {
      return boost::any_cast<const char *>(in);
    }
  }

  dict
  yaml_to_dict(const YAML::Node & node) {
    auto a = yaml_to_any(node);
    dict map = boost::any_cast<dict&>(a);
    return map;
  }

  boost::any
  yaml_to_any(const YAML::Node & node) {
    switch (node.Type()) {
      case YAML::NodeType::Scalar:
        {
          int anInt;
          double aDouble;
          bool aBool;
          if (YAML::convert<int>::decode(node, anInt))
            return anInt;
          else if (YAML::convert<double>::decode(node, aDouble))
            return aDouble;
          else if (YAML::convert<bool>::decode(node, aBool))
            return aBool;
          else {
            return node.as<std::string>();
          }
          break;
        }
      case YAML::NodeType::Map:
        {
          dict map;
          for (const auto it: node) {
            const std::string key = it.first.as<std::string>();
            map[key] = yaml_to_any(it.second);
          }
          return map;
        }
      default:
        // Sequence, Undefined, Null, whatevs
        return nullptr;
        break;
    }
  }
} // namespace mikado
