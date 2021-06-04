#include "mimir/program/Options.hpp"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <boost/program_options.hpp>

namespace mimir {

  Options::Options(const std::string& program_descriptor)
      : m_progDesc(program_descriptor)
    { }

  bool Options::parse(int argc, char *argv[])
  {

    namespace po = boost::program_options;
    po::options_description desc("Options");
    desc.add_options()
     ("help", "Display this help message")
     ("severity,v", po::value<std::uint16_t>(&severity)->default_value(3)
      ->notifier(Options::check_range<std::uint16_t>("severity", 0, 5)),
      "Print severity: Integer in [0-5]\n"
      "5:trace, 4:debug, 3:info,\n2:warning, 1:error, 0:fatal")
     ("log,l", po::value<std::uint16_t>(&log_severity)->default_value(3)
      ->notifier(Options::check_range<std::uint16_t>("log", 0, 5)),
      "Logfile severity: Integer in [0-5]\n"
      "5:trace, 4:debug, 3:info,\n2:warning, 1:error, 0:fatal")
     ;

    po::options_description hidden;
    hidden.add_options()
     ("CONFIG", po::value<std::string>(&filename)->required(), "YAML config filename")
     ;

    po::options_description all_options;
    all_options.add(desc);
    all_options.add(hidden);

    po::positional_options_description p;
    p.add("CONFIG", -1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv)
     .options(all_options)
     .positional(p)
     .run(),
     vm);

    if(vm.count("help"))
    {
      std::cerr << make_usage_string(argv[0], desc, p) << std::endl;
      return false;
    }

    po::notify(vm);

    if(vm.count("CONFIG") > 0)
    {
      if (!std::filesystem::exists(filename))
      {
        std::cerr << "File: " << filename << " does not exist" << std::endl;
        return false;
      }
    }
    return true;
  }

  std::string Options::make_usage_string(
      const std::string& program_name,
      const boost::program_options::options_description& desc,
      boost::program_options::positional_options_description& p)
  {
    // Add program one-liner
    std::ostringstream oss;
    oss
     << std::endl
     << "Usage:  "
     << std::filesystem::path(program_name).filename().string()
     << " [OPTIONS] ";

    const std::string &last = p.name_for_position(100); // max nr of options..

    for(unsigned int i = 0; i < 100; ++i)
    {
      const std::string &n = p.name_for_position(i);
      oss << n << " ";

      if(n == last) break;
    }

    oss << std::endl << std::endl
        << m_progDesc
        << std::endl << std::endl
        << desc;

    return oss.str();
  }

}
