#pragma once
/**
   @file StateMachineFwd.hpp
   @brief Forward declaration of the state machine class and its events.
*/

#include <boost/statechart/event.hpp>

namespace mimir
{
  namespace bsc = boost::statechart;

  // Events

  /// Terminate state machine.
  struct EvKill : bsc::event< EvKill > {};
  /// Start the algorithm loop.
  struct EvStart : bsc::event< EvStart > {};
  /// Stop the algorithm loop.
  struct EvStop : bsc::event< EvStop > {};
  /// A workload has finished.
  struct EvReady : bsc::event< EvReady > {};
  /// A workload has error.
  struct EvError : bsc::event< EvError > {};
  /// UNUSED currently.
  struct EvRunningOK : bsc::event< EvRunningOK > {};
  /// UNUSED currently.
  struct EvRestart : bsc::event< EvRestart > {};
  /// A timer has finished
  struct EvTimeout : bsc::event< EvTimeout > {};
  /// A workload has been canceled.
  struct EvInterrupt : bsc::event< EvInterrupt > {};

  // States
  struct Standby;      // Entry state of state machine
  struct Running;      // Parent state of:
  struct Initializing; //   Initializing
  struct Evaluating;   //   Evaluating
  struct Waiting;      //   Waiting
  struct ShutDown;     // Exit state

  
  struct StateMachine;
}

#include <boost/version.hpp>
#include <boost/intrusive_ptr.hpp>
#include <boost/statechart/asynchronous_state_machine.hpp>

#include <boost/statechart/exception_translator.hpp>

#ifndef DOXYGEN_SHOULD_SKIP_THIS

namespace boost
{
  namespace statechart
  {
#if BOOST_VERSION < 106900
    typedef void none;
#endif

    // The following class member specialization ensures that
    // state_machine<>::initiate is not instantiated at a point
    // where Standby is not yet defined.
    /*
    template<>
    inline void asynchronous_state_machine<
      mimir::StateMachine,
      mimir::Standby,
      fifo_scheduler<>,
      std::allocator< none >,
      exception_translator<> > ::initiate_impl() {}
  */ // No longer needed, or what? TODO: clean up if not needed
  }
}

template< class T>
boost::intrusive_ptr< T > make_intrusive(T * obj)
{
  return boost::intrusive_ptr<T>(obj);
}
#endif
