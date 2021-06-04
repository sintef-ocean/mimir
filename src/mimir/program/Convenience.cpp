#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

#include "mimir/program/Convenience.hpp"

#include <iostream>

namespace mimir {

boost::log::trivial::severity_level transform_severity(uint16_t severity)
{
  // invert severity to follow syslog RFC5424 (lower number is higher severity..)

  switch(severity)
  {
  case 0:
    return boost::log::trivial::fatal;
  case 1:
    return boost::log::trivial::error;
  case 2:
    return boost::log::trivial::warning;
  case 3:
    return boost::log::trivial::info;
  case 4:
    return boost::log::trivial::debug;
  case 5:
    return boost::log::trivial::trace;
  default:
    return boost::log::trivial::trace;
  }
  /*
  boost::log::core::get()->set_filter(
       boost::log::trivial::severity >= severity_to_set );
  */
}

sig_atomic_t signaled = 0;

  void SIGTERMHandler(int ) { signaled = 1; }
}
