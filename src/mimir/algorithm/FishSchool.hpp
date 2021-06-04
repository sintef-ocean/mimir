#pragma once
/**

   @file FishSchool.hpp
   @brief Fish school simulator with deterministic motion equations.

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
       @brief Kinematic fish school simulation model.

       \rst

       This class implements a simple kinematic fish school, which can be controlled by
       its inputs. The fish school is modelled as a 5 degrees of freedom kinematic
       body. It has a three dimensional position represented in an inertial reference
       frame, \{NED\}, with axes pointing North, East, and Down, respectively. The body
       can rotate about its y-axis -- pitch: :math:`\theta`, and z-axis -- yaw:
       :math:`\psi`. The fish school is modeled as an under-actuated entity, meaning that
       it can only change forward speed (surge) and its angular velocities. The input to
       the simulation model are surge :math:`u`, rate of turn :math:`r`, and desired depth
       :math:`D_d`. The desired depth is achieved with a first order response using pitch
       angular velocity as manipulated variable. The commanded surge and rate of turn is
       immidiate, i.e. there is no dynamic response. Let us define some quantities to
       describe the fish school with an ordinary differential equation.

       .. math::
          :nowrap:

          \begin{align}
            \boldsymbol{p} = \begin{bmatrix}N\\E\\D\\\end{bmatrix} &\in \mathbb{R}^3 &\text{position}\\
            \boldsymbol{\Theta} = \begin{bmatrix}\theta\\\psi\end{bmatrix} &\in \mathbb{S}^2 &\text{attitude}\\
            \boldsymbol{v} = \begin{bmatrix}u\\v\\w\end{bmatrix} &\in \mathbb{R}^3 &\text{linear velocity}\\
            \boldsymbol{w} = \begin{bmatrix}p\\q\\r\\\end{bmatrix} &\in \mathbb{R}^3 &\text{angular velocity}\\
            D_d &\in \mathbb{R} &\text{desired depth} \\
            \boldsymbol{\eta}(t) = \begin{bmatrix}\boldsymbol{p}\\\Theta\end{bmatrix}, &\quad\mathbb{R} \to \mathbb{R}^3 \times \mathbb{S}^2&\text{state vector}
          \end{align}

       Let :math:`R_y(\theta), R_z(\psi) \in SO(3)` be rotation matrixes, so that
       :math:`R(\boldsymbol{\eta})` is the chained rotation defined as follows:

       .. math::
          :nowrap:

          \begin{align*}
            R_y(\theta) &= \left[\begin{array}{ccc}\cos(\theta)&0 &\sin(\theta)\\0&1&0\\-\sin(\theta)&0&\cos(\theta) \end{array} \right] \\
            R_z(\psi) &= \left[\begin{array}{ccc}\cos(\psi)&-\sin(\psi)&0\\\sin(\psi)&\cos(\psi)&0\\0&0&1 \end{array} \right]\\
            R(\boldsymbol{\eta}) &= R_z(\psi)R_y(\theta) \\
          \end{align*}

       We then use the angular velocity transformation :math:`T(\theta)` as defined in
       :cite:`Fossen2011` and the combined state vector transformation
       :math:`J(\boldsymbol{\eta})` for our 5-dimensional system becomes:

       .. math::
          :nowrap:

          \begin{align*}
            T(\theta) &= \left[\begin{array}{cc}1&0\\0& \frac{1}{\cos(\theta)}\end{array} \right] \\
            J(\boldsymbol{\eta}) &= \left[\begin{array}{cc}R(\boldsymbol{\eta})&0\\0& T(\theta)\end{array} \right] \\
          \end{align*}

       The depth response has a proportional feedback on depth error and a stabilizing
       term to level out the pitch: :math:`q_d = 0.1(D - D_d)\cos(2\theta) - \theta`. The
       resulting ordinary differential equation and input vector are, respectively,

       .. math::
          :nowrap:

          \begin{align*}
            \dat \eta(t) &= J(\boldsymbol{\eta}) \begin{bmatrix}v\\q_d\\q\end{bmatrix} \\
            u(t_d) &= \begin{bmatrix}v_d\\r_d\\Z_d\end{bmatrix},
          \end{align*}

       where :math:`t_d` indicates discrete time points indicated by ``time_step_ms`` in
       the YAML config.

       The table below describes the inputs, outputs and parameters of the
       algorithm. These variables are specified using the YAML configuration file defined
       further below.

       +----------------+----------------------------+-------------------------------------+---------------+----------------+----------------------+-----------------+
       | Name           | Symbol                     | Description                         | Causality     | Variability    | Default              | Unit            |
       +================+============================+=====================================+===============+================+======================+=================+
       | ``fish_ctrl``  | :math:`u`                  | Desired: Surge, rate of turn, depth | ``input``     | ``discrete``   | \                    | [m/s, rad/s, m] |
       +----------------+----------------------------+-------------------------------------+---------------+----------------+----------------------+-----------------+
       | ``kinematics`` | :math:`y`                  | ``fkin::Kinematics6D``              | ``output``    | ``continuous`` | \                    | \               |
       +----------------+----------------------------+-------------------------------------+---------------+----------------+----------------------+-----------------+
       | ``position``   | :math:`p_{p,0}`            | Initial; North, East, Down          | ``parameter`` | ``fixed``      | :math:`[0, 100, 50]` | m               |
       +----------------+----------------------------+-------------------------------------+---------------+----------------+----------------------+-----------------+
       | ``euler``      | :math:`[\theta_0, \psi_0]` | Initial; Pitch, Yaw                 | ``parameter`` | ``fixed``      | :math:`[0,0]`        | rad             |
       +----------------+----------------------------+-------------------------------------+---------------+----------------+----------------------+-----------------+

       .. warning::

          This numerical model is DEPRECATED in favor of a more detailed model with
          stochastic behavior defined using FMI :cite:`fmi2`. Nevertheless, it may serve
          as an example on how to implement a simple :cpp:class:`mimir::IAlgorithm`.

       \endrst

    */
    class FishSchool : public IAlgorithm
    {
    public:
      /**
         @brief FishSchool constructor.

         The constructor parses the specification from the given YAML
         node. It establishes data structures and sets up DDS readers
         and writers according to the input/output scheme of the
         algorithm. The following YAML code block shows the expected
         layout of the ``FishSchool`` map of a input config file.  The
         inputs, outputs, and initial conditions are communicated with
         DDS communication. Common is their DDS topic and DDS
         identifier. The specification is deemed self-explanatory.

         \rst

         .. code-block:: yaml

             time_step_ms: 200
             inputs:
               fish_ctrl:                # DDS type: fkin::IdVec3d
                 topic: fkinFishCtrl
                 id: Fish
                 max_age_ms: -1
             outputs:                    # DDS type: fkin::Kinematics6D
               kinematics:
                 topic: fkinKinematics6D
                 id: Fish
             initial_conditions:
               position:                 # DDS type: fkin::IdVec3d
                 topic: fkinPosition
                 id: FishInit
                 max_wait_ms: 80
                 fallback: [0, 100, 50]
               euler:                    # DDS type: fkin::IdVec3d
                 topic: fkinEuler
                 id: FishInit
                 max_wait_ms: 80
                 fallback: [0, 0, 0]

         \endrst

         @param [in] config YAML configuration from input file.
         @param [in] scheduler State machine scheduler, needed to post events to state machine.
         @param [in] machine State machine processor handle, needed to post events to state machine.
         @param publisher DDS data writer to send data.
         @param subscriber DDS data reader to receive data.
      */
      explicit FishSchool(
          const YAML::Node& config,
          boost::statechart::fifo_scheduler<>& scheduler,
          boost::statechart::fifo_scheduler<>::processor_handle machine,
          dds::pub::Publisher publisher,
          dds::sub::Subscriber subscriber);
      /// Destructor.
      virtual ~FishSchool();
      /**
         See base class for details.

         \rst
         The algorithm uses CVODES from SUNDIALS :cite:`Hindmarsh2005sundials`.
         \endrst
      */
      virtual void solve(const std::atomic<bool>& cancel_token);
      /// See base class.
      virtual void initialize(const std::atomic<bool>& cancel_token);
      /// See base class.
      virtual void timer(const std::atomic<bool>& cancel_token);
      /// Name identifier of algorithm.
      virtual inline const char* name() { return "FishSchool"; }
      /// Helper function that queues event to state machine.
      void event(boost::statechart::event_base * const event);

    private:
      /// March simulation time one time step ahead.
      inline void step_time() { m_next_step += m_time_step; }
      FishSchool() = delete;
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
