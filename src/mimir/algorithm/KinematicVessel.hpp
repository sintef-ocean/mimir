#pragma once
/**

   @file KinematicVessel.hpp
   @brief Kinematic planar vessel.

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
       @brief Kinematic planar vessel numerical simulation model.

       \rst

       This class implements a kinematic planar vessel that are controlled by inputs,
       which are velocities and rate of turn.

       .. math::
          :nowrap:

          \begin{align}
          \boldsymbol{p} = \begin{bmatrix}N\\E\end{bmatrix} &\in \mathbb{R}^2 &\text{position}\\
          \boldsymbol{\psi} \in \mathbb{S} &&\text{heading}\\
          \boldsymbol{v} = \begin{bmatrix}u\\v\end{bmatrix} &\in \mathbb{R}^2 &\text{linear velocity}\\
          \boldsymbol{r} \in \mathbb{R} &&\text{angular velocity}\\
          x(t) = \begin{bmatrix}\boldsymbol{p}\\\psi\end{bmatrix}, &\quad\mathbb{R} \to \mathbb{R}^2 \times \mathbb{S}&\text{state vector}
          \end{align}

       .. math::
          :nowrap:

          \begin{align}
          \dat x &= \begin{bmatrix}\cos(\psi)& -\sin(\psi) &0\\
                                  \sin(\psi)& \cos(\psi) & 0\\
                                  0&0&1
                                  \end{bmatrix}
                   \begin{bmatrix}u\\v\\r\end{bmatrix}\\
          \boldsymbol{u} &= \begin{bmatrix}u\\v\\r\end{bmatrix}\\
          \boldsymbol{y} &= x
          \end{align}


       +----------------------+----------------------------+-------------------------------------+---------------+----------------+----------------------+-------------------+
       | Name                 | Symbol                     | Description                         | Causality     | Variability    | Default              | Unit              |
       +======================+============================+=====================================+===============+================+======================+===================+
       | ``vessel_ctrl``      | :math:`u`                  | Desired: Surge, Sway rate of turn   | ``input``     | ``discrete``   | \                    | [m/s, m/s, rad/s] |
       +----------------------+----------------------------+-------------------------------------+---------------+----------------+----------------------+-------------------+
       | ``kinematics``       | :math:`y`                  | ``fkin::Kinematics2D``              | ``output``    | ``continuous`` | \                    | \                 |
       +----------------------+----------------------------+-------------------------------------+---------------+----------------+----------------------+-------------------+
       | ``position_course``  | :math:`[p_0,\psi_0]`       | Initial; North, East, Yaw           | ``parameter`` | ``fixed``      | :math:`[0, 50, 0]`   | m                 |
       +----------------------+----------------------------+-------------------------------------+---------------+----------------+----------------------+-------------------+

       .. warning::

          This simulation model is very basic and has been abandoned in favor of a higher
          fidelity model implemented using FMI :cite:`fmi2`.

       \endrst

    */
    class KinematicVessel : public IAlgorithm
    {
    public:
      /**
         @brief KinematicVessel constructor.

         The constructor parses the specification from the given YAML
         node. It establishes data structures and sets up DDS readers
         and writers according to the input/output scheme of the
         algorithm. The following YAML code block shows the expected
         layout of the ``KinematicVessel`` map of a input config file. The
         inputs, outputs, and initial conditions are communicated with
         DDS communication. Common is their DDS topic and DDS
         identifier. The specification is deemed self-explanatory.

         \rst

         .. code-block:: yaml

             time_step_ms: 200
             inputs:
               vessel_ctrl:                 # DDS type: fkin::IdVec3d
                 topic: fkinVesselCtrl
                 id: Vessel
                 max_age_ms: -1
             outputs:                       # DDS type: fkin::Kinematics2D
               kinematics:
                 topic: fkinKinematics2D
                 id: Vessel
             initial_conditions:
               position_course:             # DDS type: fkin::IdVec3d
                 topic: fkinPositionCourse
                 id: Vessel
                 max_wait_ms: 50
                 fallback: [0, 50, 0]
         \endrst

         @param [in] config YAML configuration from input file.
         @param [in] scheduler State machine scheduler, needed to post events to state machine.
         @param [in] machine State machine processor handle, needed to post events to state machine.
         @param publisher DDS data writer to send data.
         @param subscriber DDS data reader to receive data.
      */
      explicit KinematicVessel(
          const YAML::Node& config,
          boost::statechart::fifo_scheduler<>& scheduler,
          boost::statechart::fifo_scheduler<>::processor_handle machine,
          dds::pub::Publisher publisher,
          dds::sub::Subscriber subscriber);
      /// Destructor.
      virtual ~KinematicVessel();
      /// See base class.
      virtual void solve(const std::atomic<bool>& cancel_token);
      /// See base class.
      virtual void initialize(const std::atomic<bool>& cancel_token);
      /// See base class.
      virtual void timer(const std::atomic<bool>& cancel_token);
      /// Name identifier of algorithm.
      virtual inline const char* name() { return "KinematicVessel"; }
      /// Helper function that queues event to state machine.
      void event(boost::statechart::event_base * const event);

    private:
      /// March simulation time one time step ahead.
      inline void step_time() { m_next_step += m_time_step; }
      KinematicVessel() = delete;
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
