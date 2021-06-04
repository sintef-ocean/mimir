#pragma once
/**

   @file TestAlgorithm.hpp
   @brief Minimal implementation of mimir:IAlgorithm interface.

*/

#include <atomic>
#include <chrono>
#include <cinttypes>
#include <memory>

#include <boost/statechart/detail/memory.hpp>
#include <boost/statechart/fifo_scheduler.hpp>
#include <boost/statechart/event_base.hpp>

#include <yaml-cpp/yaml.h>

#ifdef _MSC_VER
#pragma warning(push, 0)
#endif
#include <dds/pub/Publisher.hpp>
#include <dds/sub/Subscriber.hpp>
#ifdef _MSC_VER
#pragma warning(pop)
#endif

#include <mimir/IAlgorithm.hpp>

namespace mimir
{
  namespace algorithm
  {
    /**
       @brief Test algorithm with first order response.

       \rst

       Implements a first order response :math:`\dat x=\frac{x-x_d}{T}`, where :math:`x_d`
       the is desired output, given by ``inputs:signal:topic`` DDS signal in the YAML
       file, and :math:`T=5` is the time constant in seconds. The differential equation is
       integrated using :cite:`Hindmarsh2005sundials` and its output is published to
       ``outputs:response:topic`` every ``time_step_ms`` milliseconds.

       +----------------+------------+------------------+-------------+----------------+
       | Name           | Symbol     | Description      | Causality   | Variability    |
       +================+============+==================+=============+================+
       | ``signal``     | :math:`u`  | Desired output.  | ``input``   | ``discrete``   |
       +----------------+------------+------------------+-------------+----------------+
       | ``response``   | :math:`y`  | Output response. | ``output``  | ``continuous`` |
       +----------------+------------+------------------+-------------+----------------+

       \endrst

    */
    class TestAlgorithm : public IAlgorithm
    {
    public:
      /**
         @brief TestAlgorithm constructor.

         \rst

         The constructor parses the specification from the given YAML
         node. It establishes data structures and sets up DDS readers
         and writers according to the input/output scheme of the
         algorithm. The following YAML code block shows the expected
         layout of the ``TestAlgorithm`` map of a input config file.  The
         inputs, outputs, and initial conditions are communicated with
         DDS communication. Common is their DDS topic and DDS
         identifier. The specification is deemed self-explanatory.

         .. code-block:: yaml

             time_step_ms: 200
             inputs:
               signal:                # DDS type: fkin::IdVec1d
                 topic: test_signal
                 id: Test
                 max_age_ms: -1
             outputs:                 # DDS type: fkin::IdVec1d
               response:
                 topic: test_output
                 id: Test
         \endrst

         @param [in] config YAML configuration from input file.
         @param [in] scheduler State machine scheduler, needed to post events to state machine.
         @param [in] machine State machine processor handle, needed to post events to state machine.
         @param publisher DDS data writer to send data.
         @param subscriber DDS data reader to receive data.
      */
      explicit TestAlgorithm(
          const YAML::Node& config,
          boost::statechart::fifo_scheduler<>& scheduler,
          boost::statechart::fifo_scheduler<>::processor_handle machine,
          dds::pub::Publisher publisher,
          dds::sub::Subscriber subscriber);
      /// Destructor.
      virtual ~TestAlgorithm();
      /// See base class.
      virtual void solve(const std::atomic<bool>& cancel_token);
      /// See base class.
      virtual void initialize(const std::atomic<bool>& cancel_token);
      /// See base class.
      virtual void timer(const std::atomic<bool>& cancel_token);
      /// Name identifier of algorithm.
      virtual inline const char* name(){ return "TestAlgorithm"; }
      /// Helper function that queues event to state machine.
      void event(boost::statechart::event_base * const event);

    private:
      /// March simulation time one time step ahead.
      inline void step_time() { m_next_step += m_time_step; }
      TestAlgorithm() = delete;
      class Impl;
      std::unique_ptr<Impl> m_impl; ///< Holds the implementation of the algorithm.
      boost::statechart::fifo_scheduler<> & m_scheduler; ///< Members for state machine.
      boost::statechart::fifo_scheduler<>::processor_handle m_stateMachine; ///< Member for state machine.
      const std::chrono::milliseconds m_time_step; ///< Discrete time step.
      std::chrono::steady_clock::time_point m_next_step; ///< Simulation time point.
      const YAML::Node m_config; ///< YAML configuration for algorithm.
    };

  }
}
