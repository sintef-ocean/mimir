#pragma once
/**
   @file Options.hpp
   @brief Program options class for parsing user arguments to application.
*/

#include <cinttypes>
#include <functional>
#include <string>

#include <boost/program_options.hpp>

namespace mimir {

  /**
     @brief Helper class for program options

     This class holds the implementation of program options and arguments.

  */
  class Options
  {
  protected:
    /// Helper function for checking a range.
    template<typename T>
    static std::function<void(T)> check_range(
        char const * const opt_name,
        const T& min,
        const T& max)
    {
      return [opt_name, min, max](T value)
             {
               if ( value < min || value > max)
               {
                 throw boost::program_options::validation_error(
                     boost::program_options::validation_error::invalid_option_value,
                     opt_name,
                     std::to_string(value));
               }
             };
    }

  public:
    /**
       @brief Constructor for mimir options.

       @param [in] program_descriptor String to print to the user on the command line.
    */
    Options(const std::string& program_descriptor);
    /**
       @brief Parse input arguments

       Parses program arguments into variables.

       @param [in] argc Number of arguments passed from main()
       @param [in] argv Aarguments passed from main()
       @return bool Whether parsing was successful.
    */
    bool parse(int argc , char *argv[]);

  private:
    /// Prints a usage string in case of help or wrong use.
    std::string make_usage_string(
        const std::string& program_name,
        const boost::program_options::options_description& desc,
        boost::program_options::positional_options_description& p);

  public:
    std::string filename; ///< YAML config filename
    std::uint16_t severity, ///< Print severity.
      log_severity; ///< Log file severity. (NOTE: unused.)
  private:
    const std::string m_progDesc; ///< Program description string.
  };


}
