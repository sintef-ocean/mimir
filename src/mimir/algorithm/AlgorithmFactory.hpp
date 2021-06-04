#pragma once
/**
   @file AlgorithmFactory.hpp
   @brief Functionality for spawning instances of IAlgorithm implementations.

   When adding a new algorithm you need to:
   + Inherit from IAlgorithm and implement its virtual functions.
   + Include the header file in this factory. Link the implementation.
   + Make sure that AlgorithmCreator returns an instance of the new algorithm.

*/


#include <memory>
#include <string>

#include <boost/statechart/detail/memory.hpp>
#include <boost/statechart/fifo_scheduler.hpp>

#ifdef _MSC_VER
#pragma warning(push, 0)
#endif
#include <dds/pub/Publisher.hpp>
#include <dds/sub/Subscriber.hpp>
#ifdef _MSC_VER
#pragma warning(pop)
#endif

#include "yaml-cpp/yaml.h"
#include "mimir/IAlgorithm.hpp"
#include "mimir/algorithm/TestAlgorithm.hpp"
#include "mimir/algorithm/KinematicVessel.hpp"
#include "mimir/algorithm/FishSchool.hpp"
#include "mimir/algorithm/Leadline.hpp"
#include "mimir/algorithm/PursePlanner.hpp"

namespace mimir
{
  /**
     @brief Function for creating an instance of the mimir::IAlgorithm.

     \rstinline
     This function uses the Factory method pattern, see e.g. :cite:`GoF1994`.
     \endrstinline
     When adding new algorithm, it needs to be returned by this function.

     The factory function expects that input YAML `config` has a map
     key equal to the name identifier of the algorithm. This map is
     extracted and passed to the algorithm to be created.

     @param [in] id Identifier string for algorithm
     @param [in] config User-provided YAML configuration
     @param [in] scheduler State machine scheduler, passed from outer context (usually main)
     @param [in] machine State machine handle, passed from outer context (usually main)
     @param [in] publisher DDS publisher for the program
     @param [in] subscriber DDS subscriber for the program
     @return mimir::IAlgorithm pointer.
  */
  std::unique_ptr<mimir::IAlgorithm> AlgorithmCreator(
      const std::string& id,
      const YAML::Node& config,
      boost::statechart::fifo_scheduler<> & scheduler,
      boost::statechart::fifo_scheduler<>::processor_handle machine,
      dds::pub::Publisher publisher,
      dds::sub::Subscriber subscriber);
}
