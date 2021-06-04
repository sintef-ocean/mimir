#pragma once
/**
   @file Convenience.hpp
   @brief Convenience variables and functions for mimir.
*/
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include "boost/log/utility/setup.hpp"
#include <boost/log/expressions.hpp>
#include <csignal>

namespace mimir {
  /// Global signal for user event Ctrl-C.
  extern sig_atomic_t signaled;
  /// Invert severity to follow syslog RFC5424 (lower number is higher severity..)
  boost::log::trivial::severity_level transform_severity(uint16_t severity);
  /// Handler function that sets signaled.
  void SIGTERMHandler(int signum);
}
