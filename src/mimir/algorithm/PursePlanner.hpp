#pragma once
/**

   @file PursePlanner.hpp
   @brief Path planner algorithm wrapper.

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
       @brief Purse seine deployment planner algorithm wrapper.

       \rst

       This class wraps the optimization planner algorithm formulated in
       :cpp:class:`mimir::algorithm::PursePlannerFormulation` into the :cpp:class:`mimir::IAlgorithm`
       interface.

       .. todo::
          Add table with input, output parameters.

       \endrst

    */
    class PursePlanner : public IAlgorithm
    {
    public:

      /**
         @brief PursePlanner constructor.

         \rst

         The constructor parses the specification from the given YAML node. It establishes
         data structures and sets up DDS readers and writers according to the input/output
         scheme of the algorithm. The following YAML code block shows the expected layout
         of the ``PursePlanner`` map of a input config file.  The inputs, outputs,
         parameters, and initial conditions are communicated with DDS communication
         (settings are fixed). Common is their DDS topic and DDS identifier. The
         PursePlanner has a potentially large input file, because of the variety of
         formulation approach and solvers involved. Configuration may be propagated to
         casadi to be parsed according to a solver's supported settings. Currently, the
         YAML configuration file comes with its schema together with the specification.
         This is because some settings, for instance a solver's options must be detailed
         with data type, so that the parser is able to cast the input to appropriate type
         expected by the solver interfaces.
         See :ref:`rst/yamlconfig:Purse planner configuration schema`.

         .. code-block:: yaml

            config:
              initial_condition:     { x0: [0, 0, 0, 0, 100, 200] }
              parameters:
                setting_speed:       { topic: balder_setting_speed,        id: PursePlanner, default: 6 }
                setting_radius:      { topic: balder_setting_radius,       id: PursePlanner, default: 150 }
                aim_point_arc:       { topic: balder_aim_point_arc_length, id: PursePlanner, default: 471 }
                leadline:            { topic: balder_leadline_parameters,  id: Leadline,     default: [330, 160] }
                current_surface:     { topic: balder_current_surface,      default: [0, 0] }
                current_fish:        { topic: balder_current_fish,         default: [0, 0] }
                fish_water_velocity: { topic: balder_fish_water_velocity,  default: [2, 0] }
                fish_depth:          { topic: balder_fish_depth,           id: Fish, default: 50 }

              inputs:
                GPS_origin:          { topic: mimir_gps_origin,            default: [63.4581027, 10.3683367] }
                vessel_pos_info:     { topic: vessel_global_pos_vel }      # should max_age_ms
                vessel_gyro_info:    { topic: vessel_heading_rot }
                fish_pos_info:       { topic: fish_global_pos_vel }
                fish_relative_pos:   { topic: fish_relative_pos_3d }
              outputs:
                id: PursePlanner
                trajectory_vessel:   { topic: mimir_vessel_trajectory }
                trajectory_fish:     { topic: mimir_fish_trajectory }
                nlp_config:          { topic: mimir_nlp_config }
                nlp_stats:           { topic: mimir_nlp_stats }
                vessel_speed:        { topic: vessel_speed_command }
                vessel_heading_rate: { topic: vessel_heading_rate_command }
                deploy_position:     { topic: vessel_deploy_position }
                collide_position:    { topic: fish_collide_position }

              settings:
                plot: true
                time_step_ms: 3000
                nlp:
                  horizon_length: 240
                  formulation:
                    heading_rot_max: 0.10
                    heading_acc_max: 6.28
                    deploy_clockwise: true
                    objective:
                      vessel_fish_distance: 0.001
                      heading_rate: 0.1
                      heading_acceleration: 0.05
                      heading_post_deploy: 200
                      vessel_fish_distance_violation: 10
                      leadline_reward: 0.00
                      aim_reward: 1
                  with_callback: true
                  solver:
                    name: ipopt
                    options:
                      *ipopt-options
                  discretization:
                    *colloc-settings

                integrator:
                  name: collocation
                  options:
                    print_stats: false
                    #abstol: 1e-6
                    collocation_scheme: legendre
                    interpolation_order: 3
                    number_of_finite_elements: 60

            schema:
              *documented-elsewhere

         \endrst

         @param [in] config YAML configuration from input file.
         @param [in] scheduler State machine scheduler, needed to post events to state machine.
         @param [in] machine State machine processor handle, needed to post events to state machine.
         @param publisher DDS data writer to send data.
         @param subscriber DDS data reader to receive data.
      */
      explicit PursePlanner(
          const YAML::Node& config,
          boost::statechart::fifo_scheduler<>& scheduler,
          boost::statechart::fifo_scheduler<>::processor_handle machine,
          dds::pub::Publisher publisher,
          dds::sub::Subscriber subscriber);
      /// Destructor.
      virtual ~PursePlanner();
      /// See base class.
      virtual void solve(const std::atomic<bool>& cancel_token);
      /// See base class.
      virtual void initialize(const std::atomic<bool>& cancel_token);
      /// See base class.
      virtual void timer(const std::atomic<bool>& cancel_token);
      /// Name identifier of algorithm.
      virtual inline const char* name(){ return "PursePlanner"; }
      /// Helper function that queues event to state machine.
      void event(boost::statechart::event_base * const event);
      PursePlanner(const PursePlanner&) = delete; // noncopyable
      PursePlanner& operator=(const PursePlanner&) = delete; // noncopyable
      PursePlanner(PursePlanner&&) = delete; // noncopyable
      PursePlanner& operator=(PursePlanner&&) = delete; // noncopyable
    private:
      /// March simulation time one time step ahead.
      inline void step_time() { m_next_step += m_time_step; }
      /// Blot with Gnuplot (if enabled)
      void plot(bool do_plot);
      /// Read parameters.
      void read_parameters();
      /// Read inputs.
      void read_inputs();
      PursePlanner() = delete;
      class Impl;
      /// Holds the implementation of the algorithm.
      std::unique_ptr<Impl> m_impl; // special attention if to allow move construct/assign
      /// Pointer to implementation (model) container.
      inline Impl* model() { return m_impl.get(); }
      boost::statechart::fifo_scheduler<> & m_scheduler; ///< Members for state machine.
      boost::statechart::fifo_scheduler<>::processor_handle m_stateMachine; ///< Members for state machine.
      const std::chrono::milliseconds m_time_step; ///< Discrete time step.
      std::chrono::steady_clock::time_point m_t0,  ///< Internal clock \f$t_0\f$.
        m_next_step; ///< Simulation time point.
      std::chrono::system_clock::time_point m_t0_wall; ///< Wallclock \f$t_0\f$.
      const YAML::Node m_config; ///< YAML configuration for algorithm.
      std::uint8_t m_retries;
      bool m_keep_solution;
    };

  }
}
