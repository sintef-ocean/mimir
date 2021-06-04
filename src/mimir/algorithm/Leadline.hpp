#pragma once
/**

   @file Leadline.hpp
   @brief Basic leadline response.

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
       @brief Leadline predicted response.

       \rst

       This algorithm calculates an expected depth response of the leadline based on the
       provided input parameters. The output is the solution to the initial value problem
       over a prediction horizon of ``prediction_horizon_sec`` seconds, with
       discretization step equal to ``time_step_ms`` milliseconds.

       Let :math:`x(t) \in \mathbb{R}` be the depth in meters at a given time
       :math:`t`. Let :math:`x_d \in \mathbb{R}` and :math:`\tau \in \mathbb{R}_{>0}` be a
       setpoint depth [m] and time constant [s], respectively. Suppose :math:`t_f>0` is
       the prediction horizon and :math:`\delta t` is a descretization the step
       size. Define the set of discretized time points as :math:`\mathcal{T} := \{ t : t =
       k\delta t\, \forall k \in \mathbb{Z}_{\geq 0}, t \in [0, t_f] \}`. The solution to the
       initial value problem

       .. math::
          :nowrap:

          \begin{align}
          \dat x(t) &= \frac{x(t)-x_d}{\tau} \\
          x(0) &= 0
          \end{align}

       is :math:`x(t)` and the solution set :math:`\mathcal{X} := \{ z : z = x(t)\, \forall t \in \mathcal{T} \}`.

       The algorithm solution provides :math:`(\mathcal{T},\mathcal{X})` given :math:`(x_d,\tau)`.

       \endrst

    */
    class Leadline : public IAlgorithm
    {
    public:
      /**
         @brief Leadline constructor.

         The constructor parses the specification from the given YAML
         node. It establishes data structures and sets up DDS readers
         and writers according to the input/output scheme of the
         algorithm. The following YAML code block shows the expected
         layout of the ``Leadline`` map of a input config file.  The
         inputs, outputs, and initial conditions are communicated with
         DDS communication. Common is their DDS topic and DDS
         identifier. The specification is deemed self-explanatory.

         \rst

         .. code-block:: yaml

             time_step_ms: 200
             prediction_horizon_sec: 1000
             inputs:
               parameters:                   # DDS type: fkin::IdVec2d
                 topic: leadline_parameters
                 id: Leadline
                 default: [350, 160]
             outputs:                        # DDS type: fkin::BatchIdVec1d
               depth:
                 topic: leadline_response
                 id: Leadline

         \endrst

         @param [in] config YAML configuration from input file.
         @param [in] scheduler State machine scheduler, needed to post events to state machine.
         @param [in] machine State machine processor handle, needed to post events to state machine.
         @param publisher DDS data writer to send data.
         @param subscriber DDS data reader to receive data.
      */
      explicit Leadline(
          const YAML::Node& config,
          boost::statechart::fifo_scheduler<>& scheduler,
          boost::statechart::fifo_scheduler<>::processor_handle machine,
          dds::pub::Publisher publisher,
          dds::sub::Subscriber subscriber);
      /// Destructor.
      virtual ~Leadline();
      /// See base class.
      virtual void solve(const std::atomic<bool>& cancel_token);
      /// See base class.
      virtual void initialize(const std::atomic<bool>& cancel_token);
      /// See base class.
      virtual void timer(const std::atomic<bool>& cancel_token);
      /// Name identifier of algorithm.
      virtual inline const char* name() { return "Leadline"; }
      void event(boost::statechart::event_base * const event);

    private:
      /// March simulation time one time step ahead.
      inline void step_time() { m_next_step += m_time_step; m_now += m_time_step; }
      Leadline() = delete;
      class Impl;
      std::unique_ptr<Impl> m_impl; ///< Holds the implementation of the algorithm.
      boost::statechart::fifo_scheduler<> & m_scheduler; ///< Members for state machine.
      boost::statechart::fifo_scheduler<>::processor_handle m_stateMachine; ///< Member for state machine.
      const std::chrono::milliseconds m_time_step; ///< Discrete time step.
      std::chrono::steady_clock::time_point m_next_step; ///< Simulation time point.
      std::chrono::system_clock::time_point m_now;
      const YAML::Node m_config; ///< YAML configuration for algorithm.
    };

  }
}
