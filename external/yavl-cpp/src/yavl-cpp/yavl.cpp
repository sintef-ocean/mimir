#include <cstdio>
#include <cassert>
#include <cstdint>
#include <yaml-cpp/yaml.h>
#include "yavl-cpp/yavl.h"

using namespace YAVL;

namespace YAVL {
  template <>
  std::string ctype2str<std::string>()
  {
    return "string";
  }
  template <>
  std::string ctype2str<double>()
  {
    return "double";
  }
  template <>
  std::string ctype2str<float>()
  {
    return "float";
  }

  template <>
  std::string ctype2str<bool>()
  {
    return "bool";
  }

  template <>
  std::string ctype2str<std::uint64_t>()
  {
    return "uint64";
  }

  template <>
  std::string ctype2str<std::int64_t>()
  {
    return "int64";
  }

  template <>
  std::string ctype2str<std::uint32_t>()
  {
    return "uint32";
  }

  template <>
  std::string ctype2str<std::int32_t>()
  {
    return "int32";
  }

  template <>
  std::string ctype2str<std::uint16_t>()
  {
    return "uint16";
  }

  template <>
  std::string ctype2str<std::int16_t>()
  {
    return "int16";
  }

  template <>
  std::string ctype2str<std::uint8_t>()
  {
    return "uint8";
  }

  template <>
  std::string ctype2str<std::int8_t>()
  {
    return "int8";
  }

}

std::ostream& operator << (std::ostream& os, const Path& path)
{
    for (auto i = path.begin(); i != path.end(); ++i) {
        // no dot before list indexes and before first element
        if ((i != path.begin()) && ((*i)[0] != '[')) {
            os << '.';
        }
        os << *i;
    }
    return os;
}

std::ostream& operator << (std::ostream& os, const Exception& v)
{
    os << "REASON: " << v.why << std::endl;
    os << "  doc path: " << v.doc_path << std::endl;
    os << "  treespec path: " << v.gr_path << std::endl;
    os << std::endl;
    return os;
}

std::ostream& operator << (std::ostream& os, const Errors& v)
{
    for (auto i = v.begin(); i != v.end(); ++i) {
        os << *i;
    }
    return os;
}

const std::string& Validator::type2str(YAML::NodeType::value t)
{
    static std::string nonestr = "undefined";
    static std::string nullstr = "null";
    static std::string scalarstr = "scalar";
    static std::string liststr = "list";
    static std::string mapstr = "map";

    assert((t >= YAML::NodeType::value::Undefined) && (t <= YAML::NodeType::value::Map));

    switch (t) {
    case YAML::NodeType::value::Undefined:
        return nonestr;
    case YAML::NodeType::value::Null:
        return nullstr;
    case YAML::NodeType::value::Scalar:
        return scalarstr;
    case YAML::NodeType::value::Sequence:
        return liststr;
    case YAML::NodeType::value::Map:
        return mapstr;
    }
    assert(0);
    return nonestr;
}

int Validator::num_keys(const YAML::Node& doc)
{
    if (doc.Type() != YAML::NodeType::value::Map) {
        return 0;
    }
    int num = 0;
    for (auto i = doc.begin(); i != doc.end(); ++i) {
        num++;
    }
    return num;
}

bool Validator::validate_map(const YAML::Node& mapNode, const YAML::Node& doc)
{
    if (doc.Type() != YAML::NodeType::value::Map) {
        std::string reason = "expected map, but found " + type2str(doc.Type());
        gen_error(Exception(reason, gr_path, doc_path));
        return false;
    }

    bool ok = true;
    for (auto i = mapNode.begin(); i != mapNode.end(); ++i) {
        std::string key = i->first.as<std::string>();
        YAML::Node valueNode = i->second;

        // need to dive in to check if key is optional.
        YAML::Node docMapNode = doc[key];
        if (!docMapNode.IsDefined()) {
            std::string reason = "key: " + key + " not found.";
            gen_error(Exception(reason, gr_path, doc_path));
            ok = false;
        }
        else {
            doc_path.push_back(key);
            gr_path.push_back(key);

            ok = validate_doc(valueNode, docMapNode) && ok;

            gr_path.pop_back();
            doc_path.pop_back();
        }
    }
    return ok;
}

bool Validator::validate_leaf(const YAML::Node& gr, const YAML::Node& doc)
{
    if(gr.Type() != YAML::NodeType::value::Sequence){
      std::string reason = "expected type Sequence, but found " + type2str(gr.Type());
      gen_error(Exception(reason, gr_path, doc_path));
      return false;
    }
    assert(gr.Type() == YAML::NodeType::value::Sequence);

    YAML::Node typespec_map = gr[0];
    assert(num_keys(typespec_map) == 1);

    std::string type = typespec_map.begin()->first.as<std::string>();
    YAML::Node type_specifics = typespec_map.begin()->second;

    bool ok = true;
    if (type == "string") {
        attempt_to_convert<std::string>(doc, ok);
    }
    else if (type == "double") {
      attempt_to_convert<double>(doc, ok);
    }
    else if (type == "float") {
      attempt_to_convert<float>(doc, ok);
    }
    else if (type == "bool") {
      attempt_to_convert<bool>(doc, ok);
    }
    else if (type == "uint64") {
      attempt_to_convert<std::uint64_t>(doc, ok);
    }
    else if (type == "int64") {
      attempt_to_convert<std::int64_t>(doc, ok);
    }
    else if (type == "int" || type == "int32") {
      attempt_to_convert<std::int32_t>(doc, ok);
    }
    else if (type == "uint" || type == "uint32") {
      attempt_to_convert<std::uint32_t>(doc, ok);
    }
    else if (type == "int16") {
      attempt_to_convert<std::int16_t>(doc, ok);
    }
    else if (type == "uint16") {
      attempt_to_convert<std::uint16_t>(doc, ok);
    }
    else if (type == "int8") {
      attempt_to_convert<std::int8_t>(doc, ok);
    }
    else if (type == "uint8") {
      attempt_to_convert<std::uint8_t>(doc, ok);
    }
    else if (type == "enum") {
        ok = false;
        std::string docValue = doc.as<std::string>();
        for (auto i = type_specifics.begin(); i != type_specifics.end(); ++i) {
            if (i->as<std::string>() == docValue) {
                ok = true;
                break;
            }
        }
        if (!ok) {
            std::string reason = "enum string '" + docValue + "' is not allowed.";
            gen_error(Exception(reason, gr_path, doc_path));
        }
    }
    return ok;
}

bool Validator::validate_list(const YAML::Node& gr, const YAML::Node& doc)
{
    if (doc.Type() != YAML::NodeType::value::Sequence) {
        std::string reason = "expected list, but found " + type2str(doc.Type());
        gen_error(Exception(reason, gr_path, doc_path));
        return false;
    }

    bool ok = true;
    int n = 0;
    char buf[128];

    for (auto i = doc.begin(); i != doc.end(); ++i, ++n) {
        snprintf(buf, sizeof(buf), "[%d]", n);
        doc_path.push_back(buf);
        ok = validate_doc(gr, *i) && ok;
        doc_path.pop_back();
    }
    return ok;
}

bool Validator::validate_doc(const YAML::Node& gr, const YAML::Node& doc)
{
    bool ok = true;
    YAML::Node mapNode = gr["map"];
    YAML::Node listNode = gr["list"];
    if (mapNode.IsDefined()) {
        gr_path.push_back("map");
        ok = validate_map(mapNode, doc) && ok;
        gr_path.pop_back();
    }
    else if (listNode.IsDefined()) {
        gr_path.push_back("list");
        ok = validate_list(listNode, doc) && ok;
        gr_path.pop_back();
    }
    else {
        ok = validate_leaf(gr, doc) && ok;
    }
    return ok;
}
