#include "mimir/program/Config.hpp"

#include <sstream>
#include <yavl-cpp/yavl.h>


namespace mimir
{
  namespace program
  {
    casadi::Dict parse_config(const YAML::Node& config, const YAML::Node& schema)
    {
      /// Validate config against schema
      YAVL::Validator checkNode(schema, config);

      if(!checkNode.validate())
      {
        std::stringstream oss;
        oss << checkNode.get_errors();
        throw YAML::Exception(YAML::Mark(), oss.str());
      }

      casadi::Dict dict;
      // Read yaml tree into casadi::Dict
      read_doc("root", dict, config, schema);
      return dict["root"];
    }


    void read_doc(const std::string& key, casadi::Dict& dict, const YAML::Node& config, const YAML::Node& schema)
    {
      YAML::Node mapNode = schema["map"];
      YAML::Node listNode = schema["list"];

      if(mapNode.IsDefined())
      {
        casadi::Dict nested_dict;
        read_map(nested_dict, config, mapNode);
        dict[key] = nested_dict;
      }
      else if (listNode.IsDefined())
        read_list(key, dict, config, listNode);
      else
        read_leaf(key, dict, config, schema);
    }

    void read_map(casadi::Dict& dict, const YAML::Node& config, const YAML::Node& schema)
    {
      for (auto i = schema.begin(); i != schema.end(); ++i)
      {
        auto key = i->first.as<std::string>();
        YAML::Node valueNode = i->second;
        read_doc(key, dict, config[key], valueNode);
      }
    }

    void read_list(
        const std::string& key, casadi::Dict& dict,
        const YAML::Node& config, const YAML::Node& schema)
    {

      YAML::Node mapNode = schema["map"];
      YAML::Node listNode = schema["list"];

      if(mapNode.IsDefined())
        throw YAML::Exception(YAML::Mark(),
         std::string("casadi::Dict does not support a vector of Dict"));

      if(listNode.IsDefined())
      {
        YAML::Node deepList = listNode["list"];
        if(deepList.IsDefined())
          throw YAML::Exception(YAML::Mark(),
           std::string("casadi::Dict only support list nesting of two, vector of vectors"));

        YAML::Node typespec_map = listNode[0];
        std::string type = typespec_map.begin()->first.as<std::string>();

        if (type == "double" || type == "float")
          dict[key] = config.as<std::vector<std::vector<double>>>();
        else if (
            type == "int" || type =="uint" || type == "int64"
            || type == "int32" || type == "uint32"
            || type == "uint16" || type == "int16"
            || type == "uint8" || type == "int8")
          dict[key] = config.as<std::vector<std::vector<casadi_int>>>();
        else
          throw YAML::Exception(YAML::Mark(),
           std::string("casadi::Dict entry of vector of vectors does not support type: ") +
           type + std::string(". check key: ") + key);
      }
      else
      {
        YAML::Node typespec_map = schema[0];
        std::string type = typespec_map.begin()->first.as<std::string>();

        if (type == "string")
          dict[key] = config.as<std::vector<std::string>>();
        else if (type == "double" || type == "float")
          dict[key] = config.as<std::vector<double>>();
        else if (type == "bool")
          dict[key] = config.as<std::vector<bool>>();
        else if (
            type == "int" || type =="uint" || type == "int64"
            || type == "int32" || type == "uint32"
            || type == "uint16" || type == "int16"
            || type == "uint8" || type == "int8")
          dict[key] = config.as<std::vector<casadi_int>>();
        else if (type == "uint64")
        {
          std::string unsup("casadi::Dict does not support: uint64.\n  ");
          unsup += key + " is of type " + type;
          throw YAML::Exception(YAML::Mark(), unsup);
        }
        else
        {
          throw YAML::Exception(YAML::Mark(),
           key + std::string(" is an unsupported type: ") + type);
        }
      }
    }

    void read_leaf(
        const std::string& key, casadi::Dict& dict,
        const YAML::Node& config, const YAML::Node& schema)
    {
      YAML::Node typespec_map = schema[0];
      std::string type = typespec_map.begin()->first.as<std::string>();

      if (type == "string")
        dict[key] = config.as<std::string>();
      else if (type == "double" || type == "float")
        dict[key] = config.as<double>();
      else if (type == "bool")
        dict[key] = config.as<bool>();
      else if (
          type == "int" || type =="uint" || type == "int64"
          || type == "int32" || type == "uint32"
          || type == "uint16" || type == "int16"
          || type == "uint8" || type == "int8")
        dict[key] = config.as<casadi_int>();
      else if (type == "uint64")
      {
        std::string unsup("casadi::Dict does not support: uint64.\n  ");
        unsup += key + " is of type " + type;
        throw YAML::Exception(YAML::Mark(), unsup);
      }
      else
      {
        throw YAML::Exception(YAML::Mark(),
         key + std::string(" is an unsupported type: ") + type);
      }
    }
  }
}
