#pragma once
/**
   @file StateMachine.hpp
   @brief State machine declarations.
*/

#include <future>
#include <memory>

#include <boost/log/trivial.hpp>
#include <boost/mpl/list.hpp>
#include <boost/statechart/state.hpp>
#include <boost/statechart/transition.hpp>
#include <boost/statechart/custom_reaction.hpp>
#include <boost/statechart/deferral.hpp>

#include "yaml-cpp/yaml.h"
#include "mimir/StateMachineFwd.hpp"

namespace mimir { namespace control {
    class CommandResponder;
    class StateNotifier;
  } }
namespace mimir { class IAlgorithm; }
namespace mpl = boost::mpl;

namespace mimir
{

  /**
     @brief State machine for the algorithm.

     The core program logic is implemented an asynchronous state machine. The state
     diagram below describes it. Its main purpose is to call the functions of an
     IAlgorithm implementation according to a set pattern.

     \rst
     .. uml::
       :align: center

       !include ../static/style.puml
       skinparam backgroundColor white
       [*] -> Standby : DdsTriad(..)

       state Standby
       state ShutDown
       ShutDown: entry / post(EvKill)
       ShutDown: exit / terminate()
       state Running {
       state Initializing
       state Evaluating
       state Waiting
       [*] -do-> Initializing : AlgorithmCreator(..)


       }
       Standby -do-> Running : EvStart
       Standby -ri-> ShutDown : EvKill
       Standby -up-> Standby : EvTimeout
       ShutDown --> [*] : EvKill
       Running -up-> Standby : EvStop || EvInterrupt
       Running -do-> ShutDown : EvKill || EvError || exception_thrown
       Initializing : entry / algorithm->initialize(..)
       Initializing -do-> Evaluating : EvReady
       Evaluating --> Waiting : EvReady
       Evaluating : entry / algorithm->solve(..)
       Waiting -do-> Evaluating : EvTimeout
       Waiting : entry / algorithm->timer(..)
     \endrst

  */
  struct StateMachine : bsc::asynchronous_state_machine
  <
    StateMachine,
    Standby,
    bsc::fifo_scheduler<>,
    std::allocator<void>,
    bsc::exception_translator<>
    >
  {
    /**
       @brief Constructor for state machine

       The constructor creates instances of DDS-related stuff, such as
       + DdsTriad : Initialized domain participant with publisher and subscriber instances.
       + CommandResponder : Class for interacting with other programs via DDS.
       + StateNotifier : Class that notifies subscribers of which state the state machine is in.

       @param [in] ctx Context construct needed for the asynchronous state machine
       @param [in] config YAML config, user provided configuration file
       @param [in] notify_others Reference to external boolean for inter-thread notification of living state machine.
    */
    StateMachine(my_context ctx, const YAML::Node& config, std::atomic<bool>& notify_others);
    /// Destructor.
    ~StateMachine();
    class DdsTriad;
  public:
    /// Pointer to the DdsTriad.
    DdsTriad * Dds();
    /// Pointer to the StateNotifier.
    control::StateNotifier * Notifier(){ return m_notifier.get(); }
    /// Reference to YAML configuration.
    const YAML::Node & Config(){ return m_config; }
  private:
    //virtual void initiate_impl();
    std::unique_ptr<DdsTriad> m_dds; ///< Instance of participant, publisher and subscriber
    std::unique_ptr<control::CommandResponder> m_responder;
    std::unique_ptr<control::StateNotifier> m_notifier;
    std::atomic<bool>& m_notify_others; ///< External boolean to notify others of terminated state machine.
    const YAML::Node m_config; ///< YAML configuration for state machine and algorithm.
  };

  /// The state machine is in standby, not running algorithm.
  struct Standby : bsc::state< Standby, StateMachine >
  {
    typedef mpl::list
    <
      bsc::transition< EvStart, Running >,
      bsc::transition< EvKill, ShutDown >,
      bsc::custom_reaction< EvTimeout >
      > reactions; ///< Admissible reactions for the state.

    /// Constructor.
    Standby(my_context ctx);
    /// Destructor.
    ~Standby();
    /// Custom reaction (extra function calls) in the event of EvTimeout.
    bsc::result react(const EvTimeout &);
  };

  /// The state machine is shutdown, cleaning up before quitting application.
  struct ShutDown : bsc::state< ShutDown, StateMachine >
  {
    typedef mpl::list
    <
      bsc::custom_reaction< EvKill >
      > reactions; ///< Admissible reactions for the state.

    /// Constructor.
    ShutDown(my_context ctx);
    /// Destructor.
    ~ShutDown();
    /// Custom reaction (extra function calls) in the event of EvKill.
    bsc::result react(const EvKill &);
  };

  /// Super state for running algorithm.
  struct Running : bsc::state< Running, StateMachine, Initializing >
  {
    typedef mpl::list
    <
      bsc::transition< EvStop, Standby >,
      bsc::transition< EvKill, ShutDown >,
      bsc::transition< EvError, ShutDown >,
      bsc::transition< EvInterrupt, Standby >,
      bsc::custom_reaction< bsc::exception_thrown >
      > reactions; ///< Admissible reactions for the state.

    /**
       @brief Constructor for Running state

       This state calls mimir::AlgorithmCreator() by passing configuration and other
       arguments. The instance of the IAlgorithm is held by the Running state and is in
       scope for nested states too.
    */
    Running(my_context ctx);
    /// Destructor.
    ~Running();
    /// Handling of exception, which currently is posting to log and shut down.
    bsc::result react(const bsc::exception_thrown &);
    /// Pointer to algorithm instance.
    mimir::IAlgorithm* Algorithm();

  private:
    std::unique_ptr<mimir::IAlgorithm> m_algorithm; ///< Pointer to algorithm interface implementation.
  };

  /// The state machine is initializing the algorithm.
  struct Initializing : bsc::state< Initializing, Running >
  {
    typedef mpl::list
    <
      bsc::custom_reaction< EvReady >
      > reactions; ///< Admissible reactions for the state.

    /// Constructor.
    Initializing(my_context ctx);
    /// Destructor.
    ~Initializing();
    /// Needed to "retrieve void result of future".
    bsc::result react(const EvReady &);

  private:
    std::future<void> m_future; ///< Mechanism to get result of async call. Needs to exist until call is done.
    std::atomic<bool> m_canceled; ///< Canceled variable, which is read in the async call.
  };

  /// The state machine is evaluating/solving algorithm.
  struct Evaluating : bsc::state< Evaluating, Running >
  {
    typedef mpl::list
    <
      bsc::custom_reaction< EvReady >
      > reactions; ///< Admissible reactions for the state.

    /// Constructor.
    Evaluating(my_context ctx);
    /// Destructor.
    ~Evaluating();
    /// Needed to "retrieve void result of future".
    bsc::result react(const EvReady &);

  private:
    std::future<void> m_future; ///< Mechanism to get result of async call. Needs to exist until call is done.
    std::atomic<bool> m_canceled; ///< Canceled variable, which is read in the async call.
  };

  /// The state machine is waiting for a timer to timeout.
  struct Waiting : bsc::state< Waiting, Running >
  {
    typedef mpl::list
    <
      bsc::transition< EvTimeout, Evaluating >
      > reactions; ///< Admissible reactions for the state.

    /// Constructor.
    Waiting(my_context ctx);
    /// Destructor.
    ~Waiting();

  private:
    std::future<void> m_future; ///< Mechanism to get result of async call. Needs to exist until call is done.
    std::atomic<bool> m_canceled; ///< Canceled variable, which is read in the async call.
  };

}
