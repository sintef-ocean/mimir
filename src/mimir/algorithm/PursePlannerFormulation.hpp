#pragma once
/**

   @file PursePlannerFormulation.hpp
   @brief Path planner formulation.

*/

#include <atomic>
#include <casadi/casadi.hpp>
#include <yaml-cpp/yaml.h>

#include <mimir/program/Config.hpp>
#include "mimir/FkinDds.hpp"

namespace mimir
{
  namespace algorithm
  {

    /**
       @brief NLP problem structure.

       \rst

       This struct holds data types and problem formulation elements for the NLP formulation.

       \endrst

    */
    struct NlpProblemBuilder
    {
      /// Holds decision parameters, point constraint functions and point objectives.
      casadi::DaeBuilder daeP;
      std::vector<std::string> constraints, ///< Function names to be used as point constraints.
                               objectives;  ///< Function names to be used as terminal constraints.
      /// Key is the name of a sub system initial condition, e.g. "x_1_0". The vector of bool is true for elements with a floating initial condition.
      std::map<std::string, std::vector<bool>> x_inits;
      /// Sub systems for which the NLP consists of. Key is e.g. "x_1"
      std::map<std::string, casadi::DaeBuilder> sub_system;
      std::vector<std::string> sub_system_names; ///< Names of the sub systems.
      std::map<std::string, casadi::Slice> decision_parameter_slice, ///< Placement of a named decision parameter in \f$v\f$.
                                           parameter_slice;          ///< Placement of a named parameter in \f$p\f$.
      casadi_int parameter_dimension, ///< Number of parameters.
        decision_parameter_dimension; ///< Number of decision parameters.
      std::map<std::string, std::map<std::string, casadi::Slice>> state_slice,   ///< Placement of a named state in \f$x_i\f$, subsystem i.
                                                                  input_slice;   ///< Placement of a named input in \f$u_i\f$, subsystem i.
      std::map<std::string, casadi_int> state_dimension, ///< Number of states in named subsystem
                                        input_dimension; ///< Number of inputs in named subsystem
    };

    /**
       @brief Collection of helper data structures for the discretized NLP problem.

       This struct holds data types, which simplifies the interaction with the discretized
       NLP problem.  In particular, it holds the discretized problem dimension, slices
       into the decision variable vector, as well as variable expressions for initial and
       terminal time points of state vectors.

    */
    struct NlpStructure {
      NlpStructure() : problem_dimension(0) {}
      /// Placement of a named variable in the decision variable vector of the NLP.
      std::map<std::string, casadi::Slice> variable_slice;
      /// A named discrete variable, not necessarily in the decision vector.
      std::map<std::string, casadi::MX> discrete_variable;
      /// Dimension of decision variable vector \f$w\f$ of NLP.
      casadi_int problem_dimension;
    };

    /**
       @brief Purse seine deployment planner algorithm formulation.

       \rst

       This class implements and holds instances of functions outlined in :ref:`rst/formulation:NLP formulation overview`
       and described in :ref:`rst/planner:Path planner for deployment`.
       Please refer to these documents for mathematical details.

       \endrst

    */
    class PursePlannerFormulation
    {
    public:

      /**
         @brief PursePlannerFormulation constructor.

         \rst

         The constructor reads the configuration and the bundled schema.
         The config is the same as the one passed to PursePlanner.
         See :ref:`rst/yamlconfig:Purse planner configuration schema`.

         To choose a specific discretization technique, use the following setting.

         .. code-block:: yaml

           config:
             settings:
               nlp:
                 discretization:
                   x_1:
                     technique: collocation

         \endrst

         @param [in] config_schema YAML configuration propagated from PursePlanner.

      */
      PursePlannerFormulation(const YAML::Node& config_schema);
      PursePlannerFormulation() = default;
      ~PursePlannerFormulation();
      PursePlannerFormulation(const PursePlannerFormulation& other) = delete;
      PursePlannerFormulation& operator=(const PursePlannerFormulation& other) = delete;
      PursePlannerFormulation(PursePlannerFormulation&& other) = default;
      PursePlannerFormulation& operator=(PursePlannerFormulation&& other) = default;

      /// Returns the slice of a variable in the concatenated vector of variables.
      casadi::Slice get_slice(const std::string& name, const std::vector<casadi::MX>& vars);
      /// Sets a named decision variable in a slice of the discretized NLP.
      void set_variable(const std::string& name, const casadi::DM& value);
      /// Get statistics and status for the optimization run.
      fkin::OptiStats stats();
      /// Sets the cancel token in the NlpCallback.
      void set_abort(const std::atomic<bool>& cancel_token);
      /// Get the NLP settings for the NLP problem.
      inline const fkin::NlpConfig& nlp_config() const { return m_nlp_config; }
      /// The NlpProblemBuilder helps construct the DAE/ODE formulation of the OCP, see @rstinline :numref:`prob:cdaeivp` and :ref:`rst/planner:Path planner for deployment`. @endrstinline
      NlpProblemBuilder nlp_builder;
      /// Holds necessary data structures for the NLP problem as stated in @rstinline :numref:`prob:nlp_revisit`. @endrstinline
      casadi::NlpBuilder nlp_problem;
      casadi::Function nlp_solver;   ///< @rstinline :ref:`rst/formulation:NLP formulation overview`. @endrstinline
      /**
         @brief Model predictive helper functions

         \rst

         The map consists of helper functions for each sub system of the NLP formulation.
         The outer key is the name of the sub system, e.g. "x_2". The inner keys refers to
         various helper functions as listed below. In some cases a helper function operates on only the portion of the solution vector relevant to the sub system. In that case, we first use :ref:`rst/formulation:NLP subsystem extractor \`\`::create_system_extractor\`\``.

         - "shifter": :ref:`rst/formulation:NLP shifter \`\`::create_horizon_shifter\`\``.
         - "unpacker": :ref:`rst/formulation:NLP tuple unpacker \`\`::create_solution_unpacker\`\``.
         - "timegrid": :ref:`rst/formulation:Time grid function \`\`::create_solution_timegrid\`\``.
         - "trajectory": :ref:`rst/formulation:NLP trajectory \`\`::create_trajectory_function\`\``.

         \endrst
      */
      std::map<std::string, std::map<std::string, casadi::Function>> mpc;
      class NlpCallback;
      /// Callback to allow cancellation of NLP algorithm solving. Inherits from casadi::Callback.
      std::unique_ptr<NlpCallback> nlp_callback;
      double omega_max;
      std::uint32_t N2;
    private:
      fkin::NlpConfig m_nlp_config; ///< Configuration settings for NLP problem.
      NlpStructure m_nlp_structure; ///< Helper variables for discretized NLP problem
    };
  }
}
