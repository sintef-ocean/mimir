#pragma once
/**
   @file Config.hpp
   @brief Config file parser functions.
*/

#include <yaml-cpp/yaml.h>
#include <casadi/core/generic_type.hpp>

// Given two yaml nodes, one is config, the other is schema

// First, the config is checked using yavl-cpp validation
// Then, the node is traversed and read into a casadi::Dict based on types in the schema
// maps are considered to be Dicts, lists are vectors

namespace mimir
{

  namespace program
  {
    /**
       @brief Parse a YAML::Node into a casadi::Dict based on a schema description

       - The config node is validated against the provided schema.
       - The node is then parsed and converted into a casadi::Dict.
       - The schema's data types are used when instancing casadi::Dict's GenericType.

       @param [in] config Yaml config to be parsed
       @param [in] schema Yaml schema for the config

       @return The config as a casadi::Dict

    */
    casadi::Dict parse_config(const YAML::Node& config, const YAML::Node& schema);

    /**
       @brief Traverse config node using knowledge in schema

       This function uses the schema to call read_map(), read_list(), read_leaf(),
       depending on appropriate traversal of the config YAML::Node.

       @param [in] key Identifier in dict map to be populated
       @param [in, out] dict Map to be populated
       @param [in] config Node to be parsed
       @param [in] schema Schema with tree structure and data types

    */
    void read_doc(
        const std::string& key, casadi::Dict& dict,
        const YAML::Node& config, const YAML::Node& schema);

    /**
       @brief For each key in map, call read_doc()

       Traverse into a map and iterates over keys and repeatedly calls read_doc()

       @param [in,out] dict Map to be populated
       @param [in] config Node to be parsed
       @param [in] schema Schema with tree structure and data types

    */
    void read_map(
        casadi::Dict& dict, const YAML::Node& config, const YAML::Node& schema);

    /**
       @brief Read node as a list or list of lists

       casadi::Dict only supports vectors and vectors of vectors of a few types.  This
       function checks that the config YAML::Node adheres to these restrictions and read
       the YAML sequences as std::vector. It throws a YAML::Exception otherwise.

       @param [in] key Identifier in dict map to be populated as a vector of vector of vectors
       @param [in, out] dict Map to be populated
       @param [in] config Node to be parsed
       @param [in] schema Schema with tree structure and data types

    */
    void read_list(
        const std::string& key, casadi::Dict& dict,
        const YAML::Node& config, const YAML::Node& schema);

    /**
       @brief Read leaf node into dict

       Reads leaf YAML::Node as a basic type. Throws YAML::Exception if the requested type
       in schema is incompatible with casadi::Dict.

       @param [in] key Identifier in dict map to be populated as a basic type
       @param [in, out] dict Map to be populated
       @param [in] config Node to be parsed
       @param [in] schema Schema with tree structure and data types

    */
    void read_leaf(
        const std::string& key, casadi::Dict& dict,
        const YAML::Node& config, const YAML::Node& schema);
  }
}
