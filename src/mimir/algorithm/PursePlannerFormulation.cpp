#include <cstdint>
#include <sstream>
#include <vector>

#include <boost/log/trivial.hpp>
#include <casadi/casadi.hpp>

#include "mimir/algorithm/PursePlannerFormulation.hpp"
#include <mimir/program/Config.hpp>

class mimir::algorithm::PursePlannerFormulation::NlpCallback : public casadi::Callback {
public:
  NlpCallback(const std::string& name, const casadi::Dict& opts=casadi::Dict()){

    auto options = opts;

    m_nx = options.at("nx").as_int();
    options.erase("nx");
    m_ng = options.at("ng").as_int();
    options.erase("ng");
    m_np = options.at("np").as_int();
    options.erase("np");

    m_cancel_token = nullptr;

    construct(name, options);
  }

  ~NlpCallback() {}
  casadi_int get_n_in() override { return casadi::nlpsol_n_out(); }
  casadi_int get_n_out() override { return 1; }
  std::string get_name_in(casadi_int i) override { return casadi::nlpsol_out(i); }
  std::string get_name_out(casadi_int) override { return "ret"; }
  casadi::Sparsity get_sparsity_in(casadi_int i) override {

    using casadi::Sparsity;

    auto n = casadi::nlpsol_out(i);
    if (n == "f")
      return Sparsity::scalar();
    else if (n == "x" || n == "lam_x")
      return Sparsity::dense(m_nx);
    else if (n == "g" || n == "lam_g")
      return Sparsity::dense(m_ng);
    else if (n == "p")
      return Sparsity::dense(m_np);
    else
      return Sparsity(0,0);
  }

  std::vector<casadi::DM> eval(const std::vector<casadi::DM>& ) const override {
    casadi::DM x = 0;
    if (m_cancel_token == nullptr)
    throw std::logic_error("You must call NlpCallback::set_abort() before eval");
      else if ((*m_cancel_token))
    {
      x = 1;
      BOOST_LOG_TRIVIAL(info) << name() << " got cancel token";
    }
    return {x};
  }

  void set_abort(const std::atomic<bool>& cancel_token)
  {
    m_cancel_token = &cancel_token;
  }

private:
  casadi_int m_nx, m_ng, m_np;
  const std::atomic<bool>* m_cancel_token;
};

namespace
{

  struct NlpSettings{
    NlpSettings(const YAML::Node& discr_config, const std::vector<std::string> &sub_systems)
    {
      for ( auto && system : sub_systems)
      {
        sub_system.emplace(system, SubSystemSettings(discr_config[system]));
        BOOST_LOG_TRIVIAL(debug) << "Added system " << system;
      }
    }
    NlpSettings() = delete;

    struct SubSystemSettings{
      SubSystemSettings(){}
      SubSystemSettings(const YAML::Node& sub_system_config)
      {
        technique = sub_system_config["technique"].as<std::string>();

        if (technique == "collocation"){
          degree = sub_system_config["options"]["degree"].as<std::int32_t>();
          quadrature = sub_system_config["options"]["quadrature"].as<std::string>();
        }
        else if (technique == "multi-shooting")
          degree = 0;
        else if (technique == "single-shooting")
          degree = -1;

        checkpoints = sub_system_config["checkpoints"].as<std::uint32_t>();

        flexible_horizon = sub_system_config["horizon"]["flexible"].as<bool>();
        horizon = sub_system_config["horizon"]["length"].as<double>();
        horizon_min = sub_system_config["horizon"]["min"].as<double>();
        horizon_max = sub_system_config["horizon"]["max"].as<double>();

        elements = sub_system_config["elements"]["count"].as<std::uint32_t>();
        element_length = horizon/elements;
        element_regular_intervals = sub_system_config["elements"]["regular_intervals"].as<bool>();
        element_min = sub_system_config["elements"]["min"].as<double>();
        element_max = sub_system_config["elements"]["max"].as<double>();
      }
      std::string technique;          ///< Discretization (collocation|single-shooting|multi-shooting)
      std::string quadrature;         ///< Quadrature (Radau|Legendre) in case of collocation.
      std::int32_t degree;            ///< Degree of discretization (multi: 0, single: -1)
      std::uint32_t elements;         ///< Number of discretization elements
      std::uint32_t checkpoints;      ///< Checkpoints within one element when integrating
      bool flexible_horizon;          ///< Flexible optimization time horizon
      double horizon;                 ///< Length of horizon (typically seconds)
      double horizon_min,             ///< Minimal horizon > 0 in case of flexible horizon
        horizon_max;                  ///< Maximal horizon > min in case of flexible horizon
      bool element_regular_intervals; ///< Shall each element have identical interval
      double element_length;          ///< Nominal length of an element, e_n
      double element_min,             ///< Min fractional of e_n, >0, <1.0, when non-regular
        element_max;                  ///< Max fractional of e_n, >1.0, when non-regular
    };
    std::map<std::string, SubSystemSettings> sub_system;
  };

  bool is_custom_integrator(const std::string& candidate){
    std::vector<std::string> names{"rk4", "heun"};
    if (std::find(names.begin(), names.end(), candidate) != names.end())
      return true;
    return false;
  }

  /// Implements custom explicit Runge-Kutta methods with symbolic step size
  casadi::Function custom_integrator(
      const casadi::MX& DTin,
      const std::string &name,
      const std::string &solver,
      const casadi::MXDict &dae,
      const casadi::Dict &opts=casadi::Dict())
  {
    typedef std::vector<std::string> vecStr;
    casadi::Dict grator_options = opts;
    grator_options.erase("t0");
    grator_options.erase("tf");
    casadi::MXDict cp_dae = dae;

    try {
      cp_dae["quad"] = dae.at("quad");
    }
    catch (const std::out_of_range&) {
      cp_dae["quad"] = casadi::MX::zeros(1,1);
    }

    // make ode,quad a function
    auto f_ode = casadi::Function(
        name,
        cp_dae,
        vecStr{"x", "p", "t"},
        vecStr{"ode", "quad"},
        grator_options);

    double t0val;
    try {
      t0val = opts.at("t0");
    }
    catch (const std::out_of_range&) {
      t0val = 0.;
    }
    vecStr grator_inputs{"x0", "p"};
    vecStr grator_outputs{"xf", "qf"};
    auto t0 = casadi::MX(t0val);
    auto x0 = casadi::MX::sym("x0", cp_dae.at("x").size());
    auto p = cp_dae.at("p");
    casadi::MXDict grator;
    grator["x0"] = x0;
    grator["p"] = p;

    casadi::MX DT;

    if(DTin.is_symbolic()){
      grator_inputs.push_back("DT");
      DT = casadi::MX::sym("DT", 1);
      grator["DT"] = DT;
    } else
      DT = DTin;

    if (solver.compare("rk4") == 0)
    {
      // Runge-Kutta 4 integration, single shooting
      auto K1 = f_ode(casadi::MXDict{{"t", t0}, {"x", x0}, {"p", p}});
      auto k1 = K1.at("ode");
      auto q1 = K1.at("quad");
      auto K2 = f_ode(casadi::MXDict{{"t", t0 + DT/2}, {"x", x0 + DT/2*k1}, {"p", p}});
      auto k2 = K2.at("ode");
      auto q2 = K2.at("quad");
      auto K3 = f_ode(casadi::MXDict{{"t", t0 + DT/2}, {"x", x0 + DT/2*k2}, {"p", p}});
      auto k3 = K3.at("ode");
      auto q3 = K3.at("quad");
      auto K4 = f_ode(casadi::MXDict{{"t", t0 + DT}, {"x", x0 + DT*k3}, {"p", p}});
      auto k4 = K4.at("ode");
      auto q4 = K4.at("quad");
      grator["xf"] = x0 + DT/6*(k1 + 2*k2 + 2*k3 + k4);
      grator["qf"] = DT/6*(q1 + 2*q2 + 2*q3 + q4);
    }
    else if (solver.compare("heun") == 0)
    {
      // Heun's method
      auto K1 = f_ode(casadi::MXDict{{"t", t0}, {"x", x0}, {"p", p}});
      auto k1 = K1.at("ode");
      auto q1 = K1.at("quad");
      auto K2 = f_ode(casadi::MXDict{{"t", t0 + DT}, {"x", x0 + DT*k1}, {"p", p}});
      auto k2 = K2.at("ode");
      auto q2 = K2.at("quad");
      grator["xf"] = x0 + DT/2*(k1 + k2);
      grator["qf"] = DT/2*(q1 + q2);
    }
    else
    {
      throw std::runtime_error("Unsupported custom integrator: " + solver);
    }

    return casadi::Function(name, grator, grator_inputs, grator_outputs, grator_options);
  }

  casadi::Slice get_slice(const std::string& name, const std::vector<casadi::MX>& vars)
  {
    casadi_int start_index = 0;
    for (auto && var : vars)
    {
      if (var.name() == name)
        return casadi::Slice(start_index, start_index + var.nnz());
      else
        start_index += var.nnz();
    }
    throw std::out_of_range("Variable not found in vector: " + name);
  }

  void formulate_test_system(
      mimir::algorithm::NlpProblemBuilder &nlp,
      const ::NlpSettings&,
      const YAML::Node& )
  {
    // Formulate van der Pol oscillator, since this is given in example pack
    // The optimized solution can be compared with example pack to verify implementation.

    //typedef std::vector<double> vdob;
    using casadi::MX;
    using std::vector;

    casadi::DaeBuilder dae1, dae2;

    auto x1 = dae2.add_x("x1");
    auto x2 = dae2.add_x("x2");
    auto x = vertcat(x1, x2);
    auto u = dae2.add_u("u");
    auto p = nlp.daeP.add_p("p");
    auto z = dae2.add_q("objective"); // Dummy
    auto t01 = dae1.add_aux("t0"); // helper parameter for time shifts
    auto t02 = dae2.add_aux("t0"); // helper parameter for time shifts

    MX xdot1 = (1 - pow(x2,2))*x1 - x2 + p*u;
    MX xdot2 = x1;
    MX quad = pow(x1,2) + pow(x2,2) + pow(u,2);

    dae2.add_ode("x1dot", xdot1);
    dae2.add_ode("x2dot", xdot2);
    dae2.add_quad("L", quad);


    dae2.set_min("u", -1.);
    dae2.set_max("u", 1.);
    dae2.set_min("x1", -0.25);

    dae2.set_start("x1", 0);
    dae2.set_start("x2", 0);

    dae1.sanity_check();
    dae2.sanity_check();
    nlp.daeP.sanity_check();

    nlp.x_inits["x_1_0"] = std::vector<bool>(vertcat(dae1.x).size1(), false);
    nlp.x_inits["x_2_0"] = std::vector<bool>(vertcat(dae2.x).size1(), false);
    nlp.sub_system["x_1"] = dae1;
    nlp.sub_system["x_2"] = dae2;

    BOOST_LOG_TRIVIAL(debug) << nlp.sub_system.at("x_1").type_name() << " for van der Pol. "
                             << nlp.sub_system.at("x_1").get_str(true);
    BOOST_LOG_TRIVIAL(debug) << nlp.sub_system.at("x_2").type_name() << " for van der Pol. "
                             << nlp.sub_system.at("x_2").get_str(true);
    BOOST_LOG_TRIVIAL(debug) << nlp.daeP.type_name() << " for van der Pol. "
                             << nlp.daeP.get_str(true);
  }

  void formulate_dynamics(
      mimir::algorithm::NlpProblemBuilder &nlp,
      const ::NlpSettings& nlp_set,
      const YAML::Node& config_schema)
  {
    /* === Formulation summary ===
       ===========================
    This function is long. It consists of the following sections:
      - Config parameters: Loaded from configuration file
      - Decision parameters: Decided by the optimization problem
      - Tunable parameters: User-tunable or determined from initial conditions
      - State variables: Subsystems' differential states
      - Helper functions: some casadi functions and expressions to help formulate
      - Inputs: Subsystems' control inputs
      - Differential equations: for the subsystems
      - Path constraints: valid at every time point, intra subsystem
      - Point constraints: valid at selected time points (t0, tf), inter subsystem
      - Path objectives: Lagrange term: integral over time horizon
      - Point objectives: Mayer term: At time point, typically terminal objectives
      - Collect formulation in preparation for discretization
    */

    BOOST_LOG_TRIVIAL(trace) << "formulate_dynamics";
    using casadi::MX;
    using casadi::Function;
    using std::vector;
    casadi::DaeBuilder dae1, dae2;
    nlp.sub_system_names = std::vector<std::string>{"x_1", "x_2"};

    typedef std::vector<std::string> vstr;
    typedef std::vector<casadi::MX> vecMX;

    // === Config Parameters ===
    BOOST_LOG_TRIVIAL(trace) << "  config parameters";
    // Practitioner parameters (fixed) from config file
    auto formulation = config_schema["config"]["settings"]["nlp"]["formulation"];

    //auto psi_acc_max = formulation["heading_acc_max"].as<double>();
    //auto psi_rot_max = formulation["heading_rot_max"].as<double>();

    // Fixed parameters (tuning values read from config file)
    // Setting orientation: Clockwise (1), Anticlockwise (-1)
    double d_o = formulation["deploy_clockwise"].as<bool>() ? 1 : 0;

    // === config: decision parameters

    auto ellipse = formulation["decision_parameters"]["ellipse"];
    auto Lx_num = ellipse["along"].as<std::vector<double>>();
    auto Ly_num = ellipse["across"].as<std::vector<double>>();

    // === config: subsystem 1

    auto sys1_param = formulation["subsystem"]["x_1"];
    double gamma = sys1_param["gamma"].as<double>();                     // Gain for path constrained particles
    double k_d = sys1_param["k_d"].as<double>();                         // Proportional feedback deployment arc length
    double k_pd = sys1_param["k_pd"].as<double>();                       // Proportional feedback post deployment arc length
    double Delta_chi = sys1_param["Delta_chi"].as<double>();             // rendezvous for tilde chi (experience: should be < 1)
    double Delta = sys1_param["Delta"].as<double>();                     // Lookahead distance, larger is more conservative convergence
    double deploy_vicinity = sys1_param["deploy_vicinity"].as<double>(); // Maximal distance between vehicle convergence and deployment point
    double omega_max = sys1_param["omega_max"].as<double>();             // Maximal-ish rate of turn

    // === config: subsystem 2


    // === objective

    auto terminal_wgt = formulation["objective"]["terminal"];
    double time_penalty = terminal_wgt["time_penalty"].as<double>();
    double tdiff_penalty = terminal_wgt["tdiff_penalty"].as<double>();
    double fish_trap_slack = terminal_wgt["fish_trap_slack"].as<double>();


    // === Auxiliary and decision parameters ===
    // =========================================
    BOOST_LOG_TRIVIAL(trace) << "  config decision parameters";
    auto t01 = dae1.add_aux("t0");  dae1.set_unit("t0", "s");  // helper for time shifts
    auto t02 = dae2.add_aux("t0");  dae2.set_unit("t0", "s");  // helper for time shifts

    // === Decision parameters ===
    // Parameters that are determined by the optimization problem
    // v_min <= v <= v_max for all t in [t_0, t_f]

    // Translation along x,y-axis of rotated ellipse
    auto Lx = nlp.daeP.add_aux("w_Lx");
    nlp.daeP.set_unit("w_Lx", "m");
    nlp.daeP.set_min("w_Lx", Lx_num[0]);
    nlp.daeP.set_max("w_Lx", Lx_num[1]);
    nlp.daeP.set_start("w_Lx", (Lx_num[0] + Lx_num[1])/2.);
    auto Ly = nlp.daeP.add_aux("w_Ly");
    nlp.daeP.set_unit("w_Ly", "m");
    nlp.daeP.set_min("w_Ly", Ly_num[0]);
    nlp.daeP.set_max("w_Ly", Ly_num[1]);
    nlp.daeP.set_start("w_Ly", (Ly_num[0] + Ly_num[1])/2.);
    auto Lxy = vertcat(Lx, Ly);

    // Helper variables for t_z = max(t_diff,0);
    std::string t_z_sa("w_t_z_sa");
    auto tz_sa = nlp.daeP.add_aux(t_z_sa);
    nlp.daeP.set_unit(t_z_sa, "s");
    nlp.daeP.set_min(t_z_sa, 0);
    nlp.daeP.set_start(t_z_sa, 0);
    std::string t_z_sb("w_t_z_sb");
    auto tz_sb = nlp.daeP.add_aux(t_z_sb);
    nlp.daeP.set_unit(t_z_sb, "s");
    nlp.daeP.set_min(t_z_sb, 0);
    nlp.daeP.set_start(t_z_sb, 0);
    std::string slack_trap_fish("w_slack_trap");
    auto slack_trap = nlp.daeP.add_aux(slack_trap_fish);
    nlp.daeP.set_unit(slack_trap_fish, "m");
    nlp.daeP.set_min(slack_trap_fish, 0);
    nlp.daeP.start(slack_trap_fish, 0);


    // In case of flexible horizon or the like
    std::map<std::string, MX> Tfs;
    for (auto & sub_sys: nlp.sub_system_names)
    {
      auto sub_settings = nlp_set.sub_system.at(sub_sys);
      if (sub_settings.flexible_horizon || !sub_settings.element_regular_intervals)
      {
        casadi_int DTs = sub_settings.element_regular_intervals ?
         1 : casadi_int(sub_settings.elements);
        double elem_length =  sub_settings.element_length;
        double dt_min = sub_settings.horizon_min/sub_settings.elements;
        double dt_max = sub_settings.horizon_max/sub_settings.elements;
        auto name = sub_sys + "_DT";
        BOOST_LOG_TRIVIAL(debug)
         << "Add auxiliary for flexible horizon|non-regular intervals: " << name;
        MX DT = nlp.daeP.add_aux(name, DTs);
        nlp.daeP.guess(name, elem_length);
        nlp.daeP.set_min(name, dt_min);
        nlp.daeP.set_max(name, dt_max);

        // Create expression for Tf
        if (DTs == 1)
          Tfs[sub_sys] = sub_settings.elements*DT;
        else
          Tfs[sub_sys] = dot(MX::ones(DTs, 1), DT);

        if (!sub_settings.element_regular_intervals){
          nlp.daeP.set_min(name, elem_length*sub_settings.element_min);
          nlp.daeP.set_max(name, elem_length*sub_settings.element_max);
        }
      }
    }
    for (auto& var : Tfs)
    {
      BOOST_LOG_TRIVIAL(trace) << var;
    }
    // Collect decision parameter information
    nlp.decision_parameter_dimension = vertcat(nlp.daeP.aux).nnz();
    for (auto& var : nlp.daeP.aux)
      nlp.decision_parameter_slice[var.name()] = ::get_slice(var.name(), nlp.daeP.aux);


    // === User-defined parameters (tunable) ===
    BOOST_LOG_TRIVIAL(trace) << "  config user parameters";

    auto U_v =      nlp.daeP.add_p("U_v");      nlp.daeP.set_unit("U_v", "m/s");    // Particle speed in water
    auto q0N =      nlp.daeP.add_p("q0N");      nlp.daeP.set_unit("q0N", "m");
    auto q0E =      nlp.daeP.add_p("q0E");      nlp.daeP.set_unit("q0E", "m");
    auto q0 =       vertcat(q0N,q0E);                                               // Particle initial position (t0)
    auto V_s =      nlp.daeP.add_p("V_s");      nlp.daeP.set_unit("V_s", "m/s");    // School velocity magnitude in NED
    auto p_s0N =    nlp.daeP.add_p("p_s0N");    nlp.daeP.set_unit("p_s0N", "m");
    auto p_s0E =    nlp.daeP.add_p("p_s0E");    nlp.daeP.set_unit("p_s0E", "m");
    auto p_s0 =     vertcat(p_s0N, p_s0E);                                          // Origin of non-translated ellipse (school origin)
    auto chi_s0 =   nlp.daeP.add_p("chi_s0");   nlp.daeP.set_unit("chi_s0", "rad"); // Rotation of ellipse in NED frame. (school direction)
    // psi_s0 could be an expression or aux
    auto psi_s0 =   nlp.daeP.add_p("psi_s0");   nlp.daeP.set_unit("psi_s0", "rad"); // Rotation of ellipse in water frame. (school direction)
    // Sea current
    auto W =        nlp.daeP.add_p("W");        nlp.daeP.set_unit("W", "m/s");      // Sea current magnitude
    auto alpha =    nlp.daeP.add_p("alpha");    nlp.daeP.set_unit("alpha", "rad");  // Sea current direction
    // Sinking
    auto z_s =      nlp.daeP.add_p("z_s");      nlp.daeP.set_unit("z_s", "m");      // Fish school depth
    auto tau_ll =   nlp.daeP.add_p("tau_ll");   nlp.daeP.set_unit("tau_ll", "s");   // Time constant leadline
    auto z_d =      nlp.daeP.add_p("z_d");      nlp.daeP.set_unit("z_d", "m");      // Expected terminal sink depth
    // Deployment
    auto D_set =    nlp.daeP.add_p("D_s");      nlp.daeP.set_unit("D_s", "m");      // Setting distance before collision point
    auto Rx =       nlp.daeP.add_p("Rx");       nlp.daeP.set_unit("Rx", "m");
    auto Ry =       nlp.daeP.add_p("Ry");       nlp.daeP.set_unit("Ry", "m");
    auto Rxy =      vertcat(Rx, Ry);                                                // Ellipse Radius along x,y-axis of ellipse
    auto z_min =    nlp.daeP.add_p("z_min");    nlp.daeP.set_unit("z_min", "m");    // Minimum desired sink margin
    auto trap_min = nlp.daeP.add_p("trap_min"); nlp.daeP.set_unit("trap_min", "m"); // Minimum vessel fish distance on the backside

    // Collect parameter information
    nlp.parameter_dimension = vertcat(nlp.daeP.p).nnz();
    for (auto& var : nlp.daeP.p)
      nlp.parameter_slice[var.name()] = ::get_slice(var.name(), nlp.daeP.p);

    casadi_int idx;
    nlp.parameter_slice.emplace("setting_speed_U_v",         ::get_slice("U_v", nlp.daeP.p));
    idx = ::get_slice("q0N", nlp.daeP.p).start;
    nlp.parameter_slice.emplace("vessel_NE_q0",              casadi::Slice(idx, idx+2));
    nlp.parameter_slice.emplace("fish_speed_NE_V_s",         ::get_slice("V_s", nlp.daeP.p));
    idx = ::get_slice("p_s0N", nlp.daeP.p).start;
    nlp.parameter_slice.emplace("fish_NE_p_s0",             casadi::Slice(idx, idx+2));
    nlp.parameter_slice.emplace("fish_course_NED_chi_s",     ::get_slice("chi_s0", nlp.daeP.p));
    nlp.parameter_slice.emplace("fish_course_WF_psi_s",      ::get_slice("psi_s0", nlp.daeP.p));
    nlp.parameter_slice.emplace("current_speed_surface_W",   ::get_slice("W", nlp.daeP.p));
    nlp.parameter_slice.emplace("current_dir_surface_alpha", ::get_slice("alpha", nlp.daeP.p));
    nlp.parameter_slice.emplace("fish_depth_z_s",            ::get_slice("z_s", nlp.daeP.p));
    idx = ::get_slice("tau_ll", nlp.daeP.p).start;
    nlp.parameter_slice.emplace("leadline_tau_ll_z_d",       casadi::Slice(idx, idx+2));
    nlp.parameter_slice.emplace("aim_distance_D_s",          ::get_slice("D_s", nlp.daeP.p));
    nlp.parameter_slice.emplace("setting_Rx",                ::get_slice("Rx", nlp.daeP.p));
    nlp.parameter_slice.emplace("setting_Ry",                ::get_slice("Ry", nlp.daeP.p));
    idx = ::get_slice("Rx", nlp.daeP.p).start;
    nlp.parameter_slice.emplace("setting_radii",             casadi::Slice(idx, idx+2));
    nlp.parameter_slice.emplace("sink_margin_z_min",         ::get_slice("z_min", nlp.daeP.p));
    nlp.parameter_slice.emplace("fish_margin_d_f",           ::get_slice("trap_min", nlp.daeP.p));


    // === State variables ===
    // =======================
    BOOST_LOG_TRIVIAL(trace) << "  state variable";

    // subsystem 1
    // Particle placements
    auto wpi_0 =  dae1.add_x("wpi0");       dae1.set_unit("wpi0", "-");   // Symbolic variable for initial cond of particle projection
    auto wpi_c =  dae1.add_x("wpi_c");      dae1.set_unit("wpi_c", "-");  // Deployment collision point
    auto tau_c =  dae1.add_x("tau_c");      dae1.set_unit("tau_c", "s");  // School collision estimate
    auto wpi_d =  dae1.add_x("wpi_d");      dae1.set_unit("wpi_d", "-");  // Deployment initiation point
    auto wpi_pd = dae1.add_x("wpi_pd");     dae1.set_unit("wpi_pd", "-"); // Post deployment fish track crossing (wpi_c+pi)
    auto s_c =    dae1.add_x("s_c");        dae1.set_unit("s_c", "m");    // Deployment collision arc length
    auto s_d =    dae1.add_x("s_d");        dae1.set_unit("s_d", "m");    // Deployment initiation arc length
    auto s_pd =   dae1.add_x("s_pd");       dae1.set_unit("s_pd", "m");   // Post deployment arc length

    // Collect subsystem state information
    nlp.state_dimension["x_1"] = vertcat(dae1.x).nnz();
    nlp.state_slice["x_1"] = std::map<std::string, casadi::Slice>();
    for (auto& var : dae1.x)
      nlp.state_slice.at("x_1").emplace(var.name(), ::get_slice(var.name(), dae1.x));

    // subsystem2
    // Closed-loop particle guidance
    auto wpi_v =  dae2.add_x("wpi_v");      dae2.set_unit("wpi_v", "-");  // Particle projection parametrization
    auto qN =     dae2.add_x("qN");         dae2.set_unit("qN", "m");
    auto qE =     dae2.add_x("qE");         dae2.set_unit("qE", "m");
    auto q =      vertcat(qN, qE);          // Particle position
    auto chi =    dae2.add_x("chi");        dae2.set_unit("chi", "rad");  // Actual course vessel
    auto p_wN =   dae2.add_x("p_wN");       dae2.set_unit("p_wN", "m");
    auto p_wE =   dae2.add_x("p_wE");       dae2.set_unit("p_wE", "m");
    auto p_w =    vertcat(p_wN, p_wE);      // Origin water frame
    auto p_schN = dae2.add_x("p_schN");     dae2.set_unit("p_schN", "m");
    auto p_schE = dae2.add_x("p_schE");     dae2.set_unit("p_schE", "m");
    auto p_sch =  vertcat(p_schN, p_schE);  // School position in NED

    //auto psi_v_rot = dae2.add_x("vessel_psi_rot"); dae2.set_unit("vessel_psi_rot", "rad/s"); // rate of turn

    // Collect subsystem state information
    nlp.state_dimension["x_2"] = vertcat(dae2.x).nnz();
    nlp.state_slice["x_2"] = std::map<std::string, casadi::Slice>();
    for (auto& var : dae2.x)
      nlp.state_slice.at("x_2").emplace(var.name(), ::get_slice(var.name(), dae2.x));

    idx = ::get_slice("qN", dae2.x).start;
    nlp.state_slice.at("x_2").emplace("vessel_NE", casadi::Slice(idx, idx+2));
    idx = ::get_slice("chi", dae2.x).start;
    nlp.state_slice.at("x_2").emplace("vessel_course", casadi::Slice(idx, idx+1));
    idx = ::get_slice("p_wN", dae2.x).start;
    nlp.state_slice.at("x_2").emplace("water_NE", casadi::Slice(idx, idx+2));
    idx = ::get_slice("p_schN", dae2.x).start;
    nlp.state_slice.at("x_2").emplace("fish_NE", casadi::Slice(idx, idx+2));

    // === Helper functions and expressions ===
    // ========================================

    MX theta = MX::sym("theta",1);
    Function Rot("R", {theta},
     {
      horzcat(vertcat(cos(theta),sin(theta)),
       vertcat(-sin(theta),cos(theta)))});

    MX wpi = MX::sym("wpi", 1);              // Parametrization variable (varpi)

    // Ellipse parametrization
    MX p_p_expr = p_s0 + mtimes(
        Rot({psi_s0})[0], vertcat(
            Lxy(0) + Rxy(0)*cos(d_o*wpi),
            Lxy(1) + Rxy(1)*sin(d_o*wpi)));
    Function p_p("p_p",vecMX{wpi, p_s0, psi_s0, Lxy, Rxy},{p_p_expr});

    // Fish school parametrization (assuming constant velocity)
    MX U_ws = sqrt(V_s*V_s + W*W - 2*V_s*W*cos(chi_s0-alpha)) + 1e-3;
    MX psi_ws = atan2(V_s*sin(chi_s0)-W*sin(alpha),V_s*cos(chi_s0)-W*cos(alpha));
    /*MX p_s_expr = p_s0 + (V_s*vertcat(cos(chi_s0),sin(chi_s0))  \
     - W*vertcat(cos(alpha),sin(alpha)))*tau_c;*/
    MX p_s_expr = p_s0 + U_ws*vertcat(cos(psi_ws),sin(psi_ws))*tau_c;
    Function p_s("p_s",vecMX{tau_c, p_s0, V_s, chi_s0, W, alpha}, {p_s_expr});

    // Sink dynamics are solved explicitly
    // This expression is problematic for negative t_diff, perhaps add sink_time = max(0,t_diff)
    // after all...
    MX t_diff = MX::sym("tdiff", 1);
    Function z_sink_margin(
        "leadline_depth_margin",
        {t_diff,z_d,tau_ll,z_s},
        {z_d*(1-exp(-t_diff/tau_ll))-z_s});
    // This could be solved with a single shoot too, with predetermined steps and dt = t_diff/N3.
    // Useful in case of non-explicit solution, but increases computational complexity


    // === Control inputs ===
    // ======================
    BOOST_LOG_TRIVIAL(trace) << "  control inputs";

    //auto u_acc = dae2.add_u("vessel_psi_accel");   dae2.set_unit("vessel_psi_accel", "rad/s^2");
    //dae2.set_min("vessel_psi_accel", -psi_acc_max);  // rad/s^2 heading acceleration
    //dae2.set_max("vessel_psi_accel", psi_acc_max);   //

    nlp.input_dimension["x_1"] = vertcat(dae1.u).nnz();
    nlp.input_slice["x_1"] = std::map<std::string, casadi::Slice>();
    nlp.input_dimension["x_2"] = vertcat(dae2.u).nnz();
    nlp.input_slice["x_2"] = std::map<std::string, casadi::Slice>();


    // === Differential equations ===
    // ==============================
    BOOST_LOG_TRIVIAL(trace) << "  differential equations";

    // Expressions for wpi0 convergence
    MX p_p_start = p_p(vecMX{wpi_0, p_s0, psi_s0, Lxy, Rxy})[0];
    auto dp_p_dw0 = jacobian(p_p_start, wpi_0);  // Derivative of path wrt wpi_0
    MX chi_0 = atan2(dp_p_dw0(1), dp_p_dw0(0));
    MX epsi0 = mtimes((Rot({chi_0})[0]).T(), q0 - p_p_start);
    MX sigma0 = epsi0(0); // Along-track error for wpi0

    // Expressions for tau_c and wpi_c convergence
    MX p_p_collide = p_p(vecMX{wpi_c, p_s0, psi_s0, Lxy, Rxy})[0]; // Deployment collision
    MX p_s_collide = p_s(vecMX{tau_c, p_s0, V_s, chi_s0, W, alpha})[0]; // Fish school collision

    auto dp_p_dwc = jacobian(p_p_collide, wpi_c);  // Derivative of path wrt wpi_c
    MX chi_c = atan2(dp_p_dwc(1), dp_p_dwc(0));
    MX epsi_c = mtimes((Rot({chi_c})[0]).T(), p_s_collide - p_p_collide);
    MX sigma_c = epsi_c(0); // Along-track error for wpi_c

    auto dp_p_dtau = jacobian(p_s_collide, tau_c);
    MX chi_tau = atan2(dp_p_dtau(1), dp_p_dtau(0)); // How to avoid case
    MX epsi_tau = mtimes((Rot({chi_tau})[0]).T(), p_p_collide - p_s_collide);
    MX sigma_tau = epsi_tau(0); // Along-track error for wpi_c

    // Add arc length dynamics to find deployment position: wpi_d. Need helper state s_d
    MX deploy_expr = p_p(vecMX{wpi_d, p_s0, psi_s0, Lxy, Rxy})[0];
    MX post_deploy_expr = p_p(vecMX{wpi_pd, p_s0, psi_s0, Lxy, Rxy})[0];

    double Pi = 2*acos(0.0); // Trick to get pi

    auto dot_s_c = gamma*sigma_c;
    auto dot_s_d = -k_d*(s_d - (s_c - D_set));
    auto dot_s_pd = k_pd*(wpi_c + Pi - wpi_pd); // move ahead half an ellipse and evaluate arc length
    dae1.add_ode("dot_wpi0", gamma*sigma0/norm_2(dp_p_dw0));
    dae1.add_ode("dot_wpi_c", dot_s_c/norm_2(dp_p_dwc));
    dae1.add_ode("dot_tau_c", 0.5*gamma*sigma_tau/(norm_2(dp_p_dtau) + 1e-3)); // add 1e-3 avoid numerical issues if net fish water speed is 0
    dae1.add_ode("dot_wpi_d", dot_s_d/norm_2(jacobian(deploy_expr, wpi_d)));
    dae1.add_ode("dot_wpi_pd", dot_s_pd/norm_2(jacobian(post_deploy_expr, wpi_pd)));
    dae1.add_ode("dot_s_c", dot_s_c);
    dae1.add_ode("dot_s_d", dot_s_d);
    dae1.add_ode("dot_s_pd", dot_s_pd);

    MX p_p_proj = p_p(vecMX{wpi_v, p_s0, psi_s0, Lxy, Rxy})[0];
    auto dp_p_dw = jacobian(p_p_proj, wpi_v);  // Derivative of path wrt wpi_v
    auto d2p_p_dw2 = jacobian(dp_p_dw, wpi_v); // Derivative of path derivative wrt wpi_v
    MX chi_p = atan2(dp_p_dw(1), dp_p_dw(0));
    MX epsi = mtimes((Rot({chi_p})[0]).T(), q - p_w - p_p_proj);
    MX sigma = epsi(0); // Along-track error for wpi_v
    MX err = epsi(1);   // Cross-track error for wpi_v

    //MX not_close_to_deploy = (0.5 + atan(50*(norm_2(deploy_expr + p_w - q) - deploy_vicinity))/Pi);

    MX beta = asin((W/U_v)*sin(alpha - chi)); // Side slip.
    MX Vv = sqrt(U_v*U_v + W*W + 2*U_v*W*cos(chi - beta - alpha));//*not_close_to_deploy;
    MX kappa = (dp_p_dw(0)*d2p_p_dw2(1) - dp_p_dw(1)*d2p_p_dw2(0))/pow(norm_2(dp_p_dw),3); // Signed curvature
    MX chi_r = atan(-err/Delta);  // Relative course
    MX chi_d = chi_r + chi_p;  // Desired course
    MX sintchi = sin(chi_d)*cos(chi) - cos(chi_d)*sin(chi); // sin(chi_d - chi);// sin(tilde chi)
    MX costchi = cos(chi_d)*cos(chi) + sin(chi_d)*sin(chi); // cos(chi_d - chi);// cos(tilde chi)
    MX tildechi = atan2(sintchi,costchi); // Always in [-pi,pi]

    // It is the dynamics of the desired course for the lookahead algorithm.
    // Feedback linearizable, but appears to be ineffective, wrong, or somewhat destabilizing
    MX dot_chi_d = (-sin(2*chi_r)/(2*Delta) + kappa)*Vv*cos(chi_r) + kappa*gamma*sigma; // NOTE: unused
    MX dot_chi_d_path = kappa*(Vv*cos(chi_r) + gamma*sigma); // Feedforward path-relative rate of turn

    dae2.add_ode("dot_wpiv", (Vv*cos(chi_r) + gamma*sigma)/norm_2(dp_p_dw));
    dae2.add_ode("dot_qN", Vv*cos(chi));// + conv_N);
    dae2.add_ode("dot_qE", Vv*sin(chi));// + conv_E);
    dae2.add_ode("dot_chi", dot_chi_d_path + omega_max*(tildechi)/(sqrt(tildechi*tildechi + Delta_chi*Delta_chi)));
    dae2.add_ode("dot_p_wN", W*cos(alpha)); // Water frame origin in NED frame
    dae2.add_ode("dot_p_wE", W*sin(alpha));
    dae2.add_ode("dot_p_sN", V_s*cos(chi_s0));
    dae2.add_ode("dot_p_sE", V_s*sin(chi_s0)); // School velocity in NED

    // === Constraints ===
    BOOST_LOG_TRIVIAL(trace) << "  constraints";
    // === Path constraints ===
    // ========================
    // Defined using DaeBuilder::add_y
    // Constraints are enforced at every checkpoint of the discretization
    // f_y(p, w, x1, u1) \in Y_1
    // f_y(p, w, x2, u2) \in Y_2

    // None yet

    // === Point constraints ===
    // =========================
    // Defined using DaeBuilder::{add_fun, add_d} (slight abuse of definitions)
    // The function name and a dummy variable d have same names to get d_min <= fun <= d_max

    auto func_inputs = std::vector<std::string>{"p", "v", "x_1_f", "x_2_f", "x_1_0", "x_2_0"};
    casadi::MXDict func_map;
    auto Xp1_0 = MX::sym("x1_0", vertcat(dae1.x).nnz());
    auto Xp2_0 = MX::sym("x2_0", vertcat(dae2.x).nnz());
    func_map["p"]  = vertcat(nlp.daeP.p);
    func_map["v"]  = vertcat(nlp.daeP.aux);
    //func_map["u1"] = vertcat(NLP.dae1.u);
    //func_map["u2"] = vertcat(NLP.dae2.u);
    func_map["x_1_f"] = vertcat(dae1.x);
    func_map["x_2_f"] = vertcat(dae2.x);
    func_map["x_1_0"] = Xp1_0;
    func_map["x_2_0"] = Xp2_0; // here the get_slice would be useful

    // === Terminal/Cross constraints ===

    // Fish collision chall be in the future.
    std::string tauc_name("tc_tauc_end");
    casadi::MXDict tauc_map;
    tauc_map = func_map;
    tauc_map[tauc_name] = tau_c;
    Function tauc_end(tauc_name, tauc_map, func_inputs, vstr{tauc_name});
    nlp.daeP.add_fun(tauc_end);
    nlp.daeP.add_d(tauc_name, MX(0));
    nlp.daeP.set_min(tauc_name, 0);
    nlp.constraints.push_back(tauc_name);
    //constraint_idx.emplace(tauc_name, Slice(c_idx, c_idx+1)); c_idx +=1;

    MX deploy_to_collide_t = D_set/U_v; // Time from deploy to collide point

    // The final time of x_2 subsystem since it has non-fixed horizon, x_2_DT exists by convention
    MX Tf;
    try {
      Tf = Tfs.at("x_2");
    }
    catch (std::out_of_range&)
    {
      BOOST_LOG_TRIVIAL(error) << "The x_2 subsystem must be configured with flexible horizon";
      throw;
    }
    MX predeploy_t = Tf; // Time to arrive to deployment location
    //MX predeploy_t = tau;
    MX t_diff_expr = (tau_c - (deploy_to_collide_t + predeploy_t)) + tz_sb; // sinking duration

    std::string sinking_name("tc_sink_duration"); // constraint qualification dual formulation
    casadi::MXDict sinking_map;
    sinking_map = func_map;
    sinking_map[sinking_name] = t_diff_expr - tz_sa;
    Function sink_duration(sinking_name, sinking_map, func_inputs, vstr{sinking_name});
    nlp.daeP.add_fun(sink_duration);
    nlp.daeP.add_d(sinking_name, MX(0));
    nlp.daeP.set_min(sinking_name, 0);
    nlp.daeP.set_max(sinking_name, 0);
    nlp.constraints.push_back(sinking_name);
    //constraint_idx.emplace(sinking_name, Slice(c_idx, c_idx+1)); c_idx +=1;

    MX z_m = z_sink_margin(vecMX{t_diff_expr, z_d, tau_ll, z_s})[0] - z_min;

    // Gear sink margin shall be >= a minimal specified margin
    std::string z_name("tc_z_margin");
    casadi::MXDict z_map;
    z_map = func_map;
    z_map[z_name] = z_m;
    Function z_marginf(z_name, z_map, func_inputs, vstr{z_name});
    nlp.daeP.add_fun(z_marginf);
    nlp.daeP.add_d(z_name, MX(0));
    nlp.daeP.set_min(z_name, 0);
    nlp.constraints.push_back(z_name);
    //constraint_idx.emplace(z_name, Slice(c_idx, c_idx+1)); c_idx +=1;

    // School shall enter within ellipse before surrounding it (with margin).
    MX p_p_trap = p_p(vecMX{wpi_pd, p_s0, psi_s0, Lxy, Rxy})[0];
    auto dp_p_dwpd = jacobian(p_p_trap, wpi_pd);  // Derivative of path wrt wpi_pd
    MX chi_pd = atan2(dp_p_dwpd(1), dp_p_dwpd(0));

    MX post_deploy_t = s_pd/U_v; // Time from collide to second fish track crossing
    MX tau_pd = predeploy_t + deploy_to_collide_t + post_deploy_t;
    MX p_s_trap = p_s(vecMX{tau_pd, p_s0, V_s, chi_s0, W, alpha})[0];

    std::string fish_name("tc_fish_trap");
    casadi::MXDict fish_trap;
    fish_trap = func_map;
    // Cross-track error (depends on setting orientation d_o)
    fish_trap[fish_name] = d_o*mtimes((Rot({chi_pd})[0]).T(), p_s_trap - p_p_trap)(1) + slack_trap - trap_min;
    Function fish_trap_margin(fish_name, fish_trap, func_inputs, vstr{fish_name});
    nlp.daeP.add_fun(fish_trap_margin);
    nlp.daeP.add_d(fish_name, MX(0));
    nlp.daeP.set_min(fish_name, 0);
    nlp.constraints.push_back(fish_name);
    //constraint_idx.emplace(fish_name, Slice(c_idx, c_idx+1)); c_idx +=1;

    // Vessel shall enter within a maximal distance from the calculated deploy location
    std::string dep_name("tc_deploy_error");
    casadi::MXDict dep_map;
    dep_map = func_map;
    dep_map[dep_name] = norm_2(q - p_w - p_p(vecMX{wpi_d, p_s0, psi_s0, Lxy, Rxy})[0]);
    Function deploy_error(dep_name, dep_map, func_inputs, vstr{dep_name});
    nlp.daeP.add_fun(deploy_error);
    nlp.daeP.add_d(dep_name, MX(0));
    nlp.daeP.set_min(dep_name, 0);
    nlp.daeP.set_max(dep_name, deploy_vicinity);
    nlp.constraints.push_back(dep_name);
    //constraint_idx.emplace(dep_name, Slice(c_idx, c_idx+1)); c_idx +=1;

    // === Cross constraints ===

    // Initial value wpi_v (x2), pre-solved initial condition from (x1)
    std::string wpi_name("cc_wpi_init");
    casadi::MXDict wpi_map;
    wpi_map = func_map;
    wpi_map[wpi_name] = wpi_0 - Xp2_0(0)/*wpi_v*/;
    Function wpi_init(wpi_name, wpi_map, func_inputs, vstr{wpi_name});
    nlp.daeP.add_fun(wpi_init);
    nlp.daeP.add_d(wpi_name, MX(0));
    nlp.daeP.set_min(wpi_name, 0);
    nlp.daeP.set_max(wpi_name, 0);
    nlp.constraints.push_back(wpi_name);
    //constraint_idx.emplace(wpi_name, Slice(c_idx, c_idx+1)); c_idx +=1;


    // With non-regular intervals and fixed horizon, ensure that the sum of elements equals horizon duration
    for (auto & sub_sys: nlp.sub_system_names)
    {
      auto sub_settings = nlp_set.sub_system.at(sub_sys);
      if (!sub_settings.flexible_horizon && !sub_settings.element_regular_intervals)
      {
        auto v_name = sub_sys + "_DT";
        MX DT = nlp.daeP.var(v_name);
        double horizon = sub_settings.horizon;
        casadi_int elements = casadi_int(sub_settings.elements);

        auto final_struct = func_map;
        auto f_name = sub_sys + "_final_time";
        final_struct[f_name] = horizon - dot(MX::ones(elements, 1), DT);
        casadi::Function final_time(f_name, final_struct, func_inputs, std::vector<std::string>{f_name});
        nlp.daeP.add_fun(final_time);
        nlp.daeP.add_d(f_name, casadi::MX(0));
        nlp.daeP.set_min(f_name, 0);
        nlp.daeP.set_max(f_name, 0);
        nlp.constraints.push_back(f_name);
        BOOST_LOG_TRIVIAL(debug)
         << "Ensure fixed horizon with non-regular element time intervals: " << f_name;
      }
    }


    // === Objectives ===
    // ==================
    BOOST_LOG_TRIVIAL(trace) << "objectives";

    // "trivial" objective: minimize distance between vessel and fish.
    //MX p_diff = x_v - y_s;
    //MX dist_sq = mtimes(p_diff.T(), p_diff);

    // === Path objective ===
    // Formulation assumes the existence of a quadrature state name "quad"
    auto J_1 = dae1.add_q("objective");            dae1.set_unit("objective", "-");
    auto J_2 = dae2.add_q("objective");            dae2.set_unit("objective", "-");

    MX quad1 = MX::zeros(1,1);
    dae1.add_quad("quad", quad1);

    MX quad2 = MX::zeros(1,1); //dist_sq;
    dae2.add_quad("quad", quad2);

    // === Point (terminal) objectives ===

    std::string tf_name("obj_tf");
    casadi::MXDict tf_map;
    tf_map = func_map;
    tf_map[tf_name] = time_penalty*Tf;
    Function tf_objective(tf_name, tf_map, func_inputs, vstr{tf_name});
    nlp.daeP.add_fun(tf_objective);
    nlp.objectives.push_back(tf_name);

    std::string tdiff_name("obj_tdiff");
    casadi::MXDict tdiff_mpcc_map;
    tdiff_mpcc_map = func_map;
    tdiff_mpcc_map[tdiff_name] = tdiff_penalty*tz_sa*tz_sb;
    Function tdiff_objective(tdiff_name, tdiff_mpcc_map, func_inputs, vstr{tdiff_name});
    nlp.daeP.add_fun(tdiff_objective);
    nlp.objectives.push_back(tdiff_name);

    std::string trap_slack_name("obj_trap_slack");
    casadi::MXDict trap_slack_map;
    trap_slack_map = func_map;
    trap_slack_map[trap_slack_name] = fish_trap_slack*slack_trap;
    Function trap_slack_objective(trap_slack_name, trap_slack_map, func_inputs, vstr{trap_slack_name});
    nlp.daeP.add_fun(trap_slack_objective);
    nlp.objectives.push_back(trap_slack_name);


    // === Collect formulation ===
    // ===========================
    BOOST_LOG_TRIVIAL(trace) << "  collect formulation";

    dae1.sanity_check();
    dae2.sanity_check();
    nlp.daeP.sanity_check();

    // This allows for a non-fixed initial condition for states. the initial condition of wpi_v is free
    auto x_2_0 = std::vector<bool>(vertcat(dae2.x).size1(), false);
    auto idxes = ::get_slice("wpi_v", dae2.x); // Floating init
    x_2_0[idxes.start] = true;
    nlp.x_inits["x_1_0"] = std::vector<bool>(vertcat(dae2.x).size1(), false);;
    nlp.x_inits["x_2_0"] = x_2_0;

    nlp.sub_system["x_1"] = dae1;
    nlp.sub_system["x_2"] = dae2;


    bool more_ok = true;          // Set to false for vector states; get_str bug.

    BOOST_LOG_TRIVIAL(debug) << "Parameters: "
                             << nlp.daeP.type_name() << " for Purse Planner: "
                             << nlp.daeP.get_str(more_ok);
    for (auto & sys : nlp.sub_system_names)
      BOOST_LOG_TRIVIAL(debug)  << "Sub system " << sys << " "
                                << nlp.sub_system.at(sys).type_name() << " for Purse Planner: "
                                << nlp.sub_system.at(sys).get_str(more_ok);
  }

  void formulate_dae(
      mimir::algorithm::NlpProblemBuilder &nlp,
      const ::NlpSettings&,
      const YAML::Node& config_schema)
  {
    BOOST_LOG_TRIVIAL(trace) << "formulate_dae";
    throw std::runtime_error("formulate_dae is deprecated");
    using casadi::MX;
    using casadi::Slice;
    using std::vector;

    casadi::DaeBuilder dae1, dae2;

    auto t01 = dae1.add_aux("t0");                  // helper parameter for time shifts
    auto t02 = dae2.add_aux("t0");                  // helper parameter for time shifts

    // Formulation settings from input file
    auto formulation = config_schema["config"]["settings"]["nlp"]["formulation"];
    auto psi_rot_max = formulation["heading_rot_max"].as<double>();
    auto psi_acc_max = formulation["heading_acc_max"].as<double>();
    auto deploy_clockwise = formulation["deploy_clockwise"].as<bool>();
    auto objective = formulation["objective"];
    auto vf_dist_pen = objective["vessel_fish_distance"].as<double>();
    auto psi_rot_pen = objective["heading_rate"].as<double>();
    auto psi_acc_pen = objective["heading_acceleration"].as<double>();
    auto psi_deploy_pen = objective["heading_post_deploy"].as<double>();
    auto vf_soft_pen = objective["vessel_fish_distance_violation"].as<double>();
    //auto leadline_reward = objective["leadline_reward"].as<double>();
    auto aim_reward = objective["aim_reward"].as<double>();

    // === Parameters ===

    auto test_tau = nlp.daeP.add_aux("test_tau");  nlp.daeP.set_unit("test_tau", "s");
    nlp.daeP.set_min("test_tau", 75);
    nlp.daeP.set_max("test_tau", 150);

    // Vessel
    auto V =   nlp.daeP.add_p("vessel_speed");         nlp.daeP.set_unit("vessel_speed", "m/s"); // in water
    auto w_N = nlp.daeP.add_p("current_N");            nlp.daeP.set_unit("current_N", "m/s");
    auto w_E = nlp.daeP.add_p("current_E");            nlp.daeP.set_unit("current_E", "m/s");

    // Purse geometry parameters
    // minimal distance vessel-fish
    auto vf_min = nlp.daeP.add_p("vessel_fish_dist");  nlp.daeP.set_unit("vessel_fish_dist", "m");
    // distance from setting point to est. fish collision point along vessel track (at 0 current).
    auto D = nlp.daeP.add_p("purse_fish_collide");     nlp.daeP.set_unit("purse_fish_collide", "m");

    // Fish
    //auto Vf = nlp.daeP.add_p("fish_speed");            nlp.daeP.set_unit("fish_speed", "m/s");   // in water
    //auto psi_f = nlp.daeP.add_p("fish_heading");       nlp.daeP.set_unit("fish_heading", "rad"); // heading
    auto v_Nf = nlp.daeP.add_p("fish_velocity_N");     nlp.daeP.set_unit("fish_velocity_N", "m/s");   // in water
    auto v_Ef = nlp.daeP.add_p("fish_velocity_E");     nlp.daeP.set_unit("fish_velocity_E", "m/s");   // in water
    auto w_Nf = nlp.daeP.add_p("current_N_f");         nlp.daeP.set_unit("current_N_f", "m/s");  // at fish depth
    auto w_Ef = nlp.daeP.add_p("current_E_f");         nlp.daeP.set_unit("current_E_f", "m/s");  // at fish depth
    auto z_f =  nlp.daeP.add_p("fish_depth");          nlp.daeP.set_unit("fish_depth", "m");

    // Leadline
    auto tau = nlp.daeP.add_p("leadline_tau");         nlp.daeP.set_unit("leadline_tau", "s");     // Time constant
    auto z_d = nlp.daeP.add_p("leadline_depth_d");     nlp.daeP.set_unit("leadline_depth_d", "m"); // Expected depth

    // Add dummy parameters
    nlp.daeP.add_p("setting_Ry");
    nlp.daeP.add_p("z_min");    nlp.daeP.set_unit("z_min", "m");    // Minimum desired sink margin
    nlp.daeP.add_p("trap_min"); nlp.daeP.set_unit("trap_min", "m"); // Minimum vessel fish distance on the backside
    nlp.daeP.add_p("q0N");      nlp.daeP.set_unit("q0N", "m");
    nlp.daeP.add_p("q0E");      nlp.daeP.set_unit("q0E", "m");
    nlp.daeP.add_p("p_s0N");    nlp.daeP.set_unit("p_s0N", "m");
    nlp.daeP.add_p("p_s0E");    nlp.daeP.set_unit("p_s0E", "m");
    nlp.daeP.add_p("V_s");      nlp.daeP.set_unit("V_s", "m/s");    // School velocity magnitude in NED
    nlp.daeP.add_p("chi_s0");   nlp.daeP.set_unit("chi_s0", "rad"); // Rotation of ellipse in NED frame.
    nlp.daeP.add_p("psi_s0");   nlp.daeP.set_unit("psi_s0", "rad"); // Rotation of ellipse in water frame.
    nlp.daeP.add_p("W");        nlp.daeP.set_unit("W", "m/s");      // Sea current magnitude
    nlp.daeP.add_p("alpha");    nlp.daeP.set_unit("alpha", "rad");  // Sea current direction

    // === States ===
    auto w1 = dae1.add_x("dummy");
    dae1.add_ode("dummy", (1/test_tau)*(100-w1));

    auto x_n = dae2.add_x("vessel_N");             dae2.set_unit("vessel_N", "m");
    auto x_e = dae2.add_x("vessel_E");             dae2.set_unit("vessel_E", "m");
    auto psi_v = dae2.add_x("vessel_heading");     dae2.set_unit("vessel_heading", "rad");
    auto psi_v_rot = dae2.add_x("vessel_psi_rot"); dae2.set_unit("vessel_psi_rot", "rad/s"); // rate of turn
    auto y_n = dae2.add_x("fish_N");               dae2.set_unit("fish_N", "m");
    auto y_e = dae2.add_x("fish_E");               dae2.set_unit("fish_E", "m");

    //auto z = dae.add_x("leadline");               dae.set_unit("leadline", "m");


    auto u_acc = dae2.add_u("vessel_psi_accel");   dae2.set_unit("vessel_psi_accel", "rad/s^2");

    auto J_1 = dae1.add_q("objective");            dae1.set_unit("objective", "-");
    auto J = dae2.add_q("objective");              dae2.set_unit("objective", "-");
    //auto t_d = dae.add_aux("t_deploy");           // deployment time point

    MX t_d = 34;

    bool more_ok = true;          // Set to false for vector states; get_str bug.


    // === Define equations ===
    // Assume the following input schemes
    // f_ode(x,u,p,t)
    // f_y(x,u,p,t)   , no z or other
    // f_quad(x,u,p,t), no q or other

    // === ODE ===

    // Vessel
    dae2.add_ode("x_n_dot", V*cos(psi_v) + w_N);
    dae2.add_ode("x_e_dot", V*sin(psi_v) + w_E);
    dae2.add_ode("psi_v_dot", psi_v_rot);
    dae2.add_ode("psi_v_roc_dot", u_acc);

    // Fish
    MX y_n_dot = v_Nf + w_Nf;
    MX y_e_dot = v_Ef + w_Ef;

    dae2.add_ode("y_n_dot", y_n_dot); //Vf*cos(psi_f) + w_Nf);
    dae2.add_ode("y_e_dot", y_e_dot); //Vf*sin(psi_f) + w_Ef);

    MX fish_speed = sqrt(pow((y_n_dot),2) + pow((y_e_dot),2));
    MX fish_line = vertcat(y_n_dot/(fish_speed+0.001), y_e_dot/(fish_speed+0.001));

    // Leadline
    MX t_s = t_d + D/V; // time start of sink response
    //dae.add_ode("z_dot", -(1/tau)*(z - z_d)*(dae.t > t_s)); // maybe add/sub t0 also?


    // === Variable constraints ===
    // note: by expression with dim > 1 does not work

    // Vessel
    // x_n, x_e, psi_v are unconstrained
    dae2.set_min("vessel_psi_rot", -psi_rot_max);    // rad/s heading rate of change
    dae2.set_max("vessel_psi_rot", psi_rot_max);     //
    dae2.set_min("vessel_psi_accel", -psi_acc_max);  // rad/s^2 heading acceleration
    dae2.set_max("vessel_psi_accel", psi_acc_max);   //

    // Fish: unconstrained
    // Leadline: unconstrained

    // === Constraint expressions ===

    // Soft constr. ensuring feasible solution when vessel and fish are closer than vf_min
    auto vf_soft = dae2.add_u("vf_soft");  dae2.set_unit("vf_soft", "m");
    dae2.set_min("vf_soft", 0);

    // Distance between fish and vessel (squared) minus minimal dist vf_min
    MX vf_diff = vertcat(x_n - y_n, x_e - y_e);

    MX vf_sq_dist = mtimes(mtimes(vf_diff.T(), MX::eye(2)), vf_diff) - vf_min*vf_min + vf_soft*vf_soft;
    dae2.add_y("vf_sq_dist", vf_sq_dist);
    dae2.set_min("vf_sq_dist", 0);

    // Ensure setting orientation (assuming clockwise here)
    // Ideal rate of change?
    MX orientation = deploy_clockwise ? 1 : -1; // 1: clockwise, -1: counter-clockwise
    dae2.add_y("setting_orientation", orientation*psi_v_rot*0.5*(1 + tanh(10*(dae2.t - t_d)))); // (dae2.t > t_d)
    dae2.set_min("setting_orientation", 0);

    // === Objective function terms ===
    // gains
    double q = vf_dist_pen;               // vessel-fish distance
    double r_rot = psi_rot_pen;           // penalize aggressive maneuvers
    double r_acc = psi_acc_pen;           // heading acceleration
    double r_deployed = psi_deploy_pen;   // penalize change of heading after deploy
    double vf_soft_gain = vf_soft_pen;    // vessel-fish min distance violation penalty
    //double z_reward = leadline_reward;    // reward sinking depth
    double fish_aim = aim_reward;

    // Most gain parameters above are tuning, others may be exposed to user, esp. z_reward

    // Minimizing distance between vessel and fish vf_diff, defined above
    MX Q = q*MX::eye(2);
    MX quad_vf_diff = mtimes(mtimes(vf_diff.T(), Q), vf_diff); //+ q*vf_min*vf_min;

    // Minimizing use of heading accelerations, esp. after deploy
    MX R = r_acc*MX::eye(1);
    R = R + r_deployed*MX::eye(1)*0.5*(1 + tanh(10*(dae2.t - t_d)));    // maybe add t0 also
    MX quad_u_acc = mtimes(mtimes(u_acc.T(), R), u_acc);
    MX quad_psi_rot = r_rot*mtimes(psi_v_rot, psi_v_rot);

    // Soft constraint of vessel-fish distance violation
    MX quad_u_vf_soft = vf_soft_gain*vf_soft;


    // Distance from collision point to fish line
    // ((a-p).n)*n projected length onto the line
    // a is a point on the line, n is direction of the line
    // p is the point in question
    //(a-p) -(a-p)dot n)n, perp vector

    MX fish_pos = vertcat(y_n, y_e);
    MX vessel_pos = vertcat(x_n, x_e);
    MX perp_vec_aim = (vf_diff - mtimes(dot(vf_diff, fish_line), fish_line));

    MX quad_aim = mtimes(perp_vec_aim.T(), perp_vec_aim)*fish_aim*18*
     (1 + tanh((dae2.t - t_s) - 1))*(1 - tanh((dae2.t - t_s) + 1));

    // Need contraint to ensure that the point is in front of the fish, not behind.

    // Maximize leadline sink depth (aka. -max(.))
    //MX quad_z = -z_reward*z*(dae.t > t_s);
    //MX quad_z = 0;

    //MX quad_aim = -fish_aim*mtimes(mtimes(vf_diff.T(), MX::eye(2)), vf_diff)*(dae.t > t_d)*(dae.t <t_d + 2);

    // How to take into account sink depth as a function of setting orientation, relevant?



    // Add a terminal constraint
    auto func_inputs = std::vector<std::string>{"p", "v", "x_1_f", "x_2_f", "x_1_0", "x_2_0"};
    casadi::MXDict func_map;
    auto Xp1_0 = MX::sym("x1_0", vertcat(dae1.x).nnz());
    auto Xp2_0 = MX::sym("x2_0", vertcat(dae2.x).nnz());
    func_map["p"]  = vertcat(nlp.daeP.p);
    func_map["v"]  = vertcat(nlp.daeP.aux);
    //func_map["u1"] = vertcat(NLP.dae1.u);
    //func_map["u2"] = vertcat(NLP.dae2.u);
    func_map["x_1_f"] = vertcat(dae1.x);
    func_map["x_2_f"] = vertcat(dae2.x);
    func_map["x_1_0"] = Xp1_0;
    func_map["x_2_0"] = Xp2_0; // here the get_slice would be useful

    auto term_constr = func_map;
    term_constr["dummy_end"] = 75 - w1;
    casadi::Function dummy_end("dummy_end", term_constr, func_inputs, std::vector<std::string>{"dummy_end"});
    nlp.daeP.add_fun(dummy_end);
    nlp.daeP.add_d("dummy_end", casadi::MX(0));
    nlp.daeP.set_min("dummy_end", 0);
    nlp.daeP.set_max("dummy_end", 0);
    nlp.constraints.push_back("dummy_end");

    term_constr.erase("dummy_end");
    term_constr["dummy_cost"] = test_tau;
    casadi::Function dummy_cost("dummy_cost", term_constr, func_inputs, std::vector<std::string>{"dummy_cost"});
    nlp.daeP.add_fun(dummy_cost);
    nlp.objectives.push_back("dummy_cost");

    // === Objective quadrature ===
    MX quad = quad_vf_diff + quad_u_acc + quad_psi_rot + quad_u_vf_soft + quad_aim; // + quad_z;
    dae2.add_quad("quad", quad);

    MX quad1 = MX::zeros(1,1);
    dae1.add_quad("quad", quad1);

    dae1.sanity_check();
    dae2.sanity_check();
    nlp.daeP.sanity_check();

    auto x_1_0 = std::vector<bool>(vertcat(dae1.x).size1(), false);
    auto idx = ::get_slice("dummy", dae1.x);
    x_1_0[idx.start] = true; // we know this is a scalar, otherwise a loop repeated idx.stop times
    nlp.x_inits["x_1_0"] = x_1_0;
    nlp.x_inits["x_2_0"] = std::vector<bool>(vertcat(dae2.x).size1(), false);

    nlp.sub_system["x_1"] = dae1;
    nlp.sub_system["x_2"] = dae2;
    nlp.sub_system_names = std::vector<std::string>{"x_1", "x_2"};

    for (auto& var : nlp.daeP.p)
      nlp.parameter_slice[var.name()] = ::get_slice(var.name(), nlp.daeP.p);

    for (auto& var : nlp.daeP.aux)
      nlp.decision_parameter_slice[var.name()] = ::get_slice(var.name(), nlp.daeP.aux);

    nlp.parameter_dimension = vertcat(nlp.daeP.p).nnz();
    nlp.decision_parameter_dimension = vertcat(nlp.daeP.aux).nnz();

    casadi_int idex;
    nlp.parameter_slice.emplace("setting_speed_U_v",   Slice(0, 1));
    nlp.parameter_slice.emplace("current_surface",     Slice(1, 3));
    nlp.parameter_slice.emplace("setting_Rx",          Slice(3, 4));
    nlp.parameter_slice.emplace("aim_distance_D_s",    Slice(4, 5));
    nlp.parameter_slice.emplace("fish_water_velocity", Slice(5, 7));
    nlp.parameter_slice.emplace("current_fish",        Slice(7, 9));
    nlp.parameter_slice.emplace("fish_depth_z_s",      Slice(9, 10));
    nlp.parameter_slice.emplace("leadline_tau_ll_z_d", Slice(10, 12));
    nlp.parameter_slice.emplace("setting_Ry",                ::get_slice("setting_Ry", nlp.daeP.p));
    nlp.parameter_slice.emplace("sink_margin_z_min",         ::get_slice("z_min", nlp.daeP.p));
    nlp.parameter_slice.emplace("fish_margin_d_f",           ::get_slice("trap_min", nlp.daeP.p));
    idex = ::get_slice("q0N", nlp.daeP.p).start;
    nlp.parameter_slice.emplace("vessel_NE_q0",              Slice(idex, idex + 2));
    idex = ::get_slice("p_s0N", nlp.daeP.p).start;
    nlp.parameter_slice.emplace("fish_NE_p_s0",              Slice(idex, idex + 2));
    nlp.parameter_slice.emplace("fish_speed_NE_V_s",         ::get_slice("V_s", nlp.daeP.p));
    nlp.parameter_slice.emplace("fish_course_NED_chi_s",     ::get_slice("chi_s0", nlp.daeP.p));
    nlp.parameter_slice.emplace("fish_course_WF_psi_s",      ::get_slice("psi_s0", nlp.daeP.p));
    nlp.parameter_slice.emplace("current_speed_surface_W",   ::get_slice("W", nlp.daeP.p));
    nlp.parameter_slice.emplace("current_dir_surface_alpha", ::get_slice("alpha", nlp.daeP.p));

    BOOST_LOG_TRIVIAL(debug) << nlp.sub_system.at("x_1").type_name() << " for Purse Planner: "
                             << nlp.sub_system.at("x_1").get_str(more_ok);
    BOOST_LOG_TRIVIAL(debug) << nlp.sub_system.at("x_2").type_name() << " for Purse Planner: "
                             << nlp.sub_system.at("x_2").get_str(more_ok);
    BOOST_LOG_TRIVIAL(debug) << nlp.daeP.type_name() << " for Purse Planner: "
                             << nlp.daeP.get_str(more_ok);

    for (auto& el : nlp.parameter_slice)
      BOOST_LOG_TRIVIAL(info) << el.first << " has a " << el.second;

  }

  void make_bounds(
      const casadi::DaeBuilder& dae,
      const std::vector<casadi::MX>& vars,
      std::vector<double>& lower_bound,
      std::vector<double>& upper_bound,
      std::vector<double>& inits)
  {
    BOOST_LOG_TRIVIAL(trace) << "make_bounds";
    // set lower and upper bounds for specified var type of dae vars, u, x, z, q, y
    // vector variables are assumed to have the same bounds

    for (auto& e: vars)
    {
      auto lb = std::vector<double>(e.nnz(), dae.min(e.name()));
      auto ub = std::vector<double>(e.nnz(), dae.max(e.name()));
      auto guess = std::vector<double>(e.nnz(), dae.start(e.name()));
      lower_bound.insert(lower_bound.end(), lb.begin(), lb.end());
      upper_bound.insert(upper_bound.end(), ub.begin(), ub.end());
      inits.insert(inits.end(), guess.begin(), guess.end());
    }
  }

  void add_nlp_variable_bounds(
      casadi::NlpBuilder& prob,
      const casadi::DaeBuilder& dae,
      uint32_t elements,
      int32_t collocations,
      bool decision_parameters=false)
  {
    BOOST_LOG_TRIVIAL(trace) << "add_nlp_variable_bounds";
    // Current structure of decision variables for problem formulation are
    // w = vec(v, U, X0, Xc_1, ... Xc_N), where Xc_i is the structure of one discrete element
    // Note: this builds the decision variables for sub problems daeP or indicated dae.
    // It is assumed when decision_parameters=true, dae.u, dae.x is empty.
    // Question: Should variables such as u, z, q be ordered with each element or not?
    typedef std::vector<double> vecDob;
    using casadi::DM;
    using casadi::MX;

    int32_t d = collocations;       // For multi shooting: 0, single shooting: -1
    uint32_t N = elements;          // Number of collocation elements
    auto lbw = std::vector<DM>();   // Lower bounds on nlp variables
    auto ubw = std::vector<DM>();   // Upper bounds
    auto w = std::vector<DM>();     // Decision variable cold start

    auto low_bound = vecDob();                     // Lower bounds temporary holders
    auto up_bound = vecDob();                      // Upper bounds
    auto start_guess = vecDob();                   // Cold guess for variable

    // optionally add bounds for decision parameters (aux variables)
    if (decision_parameters)
    {
      BOOST_LOG_TRIVIAL(debug) << "Decision parameters";
      make_bounds(dae, dae.aux, low_bound, up_bound, start_guess);
      lbw.push_back(DM(low_bound));
      ubw.push_back(DM(up_bound));
      w.push_back(DM(start_guess));

      low_bound.clear();
      up_bound.clear();
      start_guess.clear();
    }

    // add bounds for u and repeat on all collocation elements
    make_bounds(dae, dae.u, low_bound, up_bound, start_guess);
    lbw.push_back(kron(DM::ones(N, 1), DM(low_bound)));
    ubw.push_back(kron(DM::ones(N, 1), DM(up_bound)));
    w.push_back(kron(DM::ones(N, 1), DM(start_guess)));

    low_bound.clear();
    up_bound.clear();
    start_guess.clear();

    // add bound for x and repeat for all collocation elements Xk_end (and X0)
    // Shooting: d=0 (seams and X0), d=-1 (only X0), for multi and single, respectively
    make_bounds(dae, dae.x, low_bound, up_bound, start_guess);
    lbw.push_back(kron(DM::ones(1 + N*(d+1), 1), DM(low_bound)));
    ubw.push_back(kron(DM::ones(1 + N*(d+1), 1), DM(up_bound)));
    w.push_back(kron(DM::ones(1 + N*(d+1), 1), DM(start_guess)));

    low_bound.clear();
    up_bound.clear();
    start_guess.clear();

    // Initial condition X0 starts at relative index N*nu for all techniques.
    prob.x_lb = vertcat(DM(prob.x_lb), vertcat(lbw)).get_elements();
    prob.x_ub = vertcat(DM(prob.x_ub), vertcat(ubw)).get_elements();
    prob.x_init = vertcat(DM(prob.x_init), vertcat(w)).get_elements();
    // Is 0 a sensible cold start for lagrange multipliers?
    prob.lambda_init = vecDob(prob.x_init.size(), 0.);


  }

  void add_point_constraints_and_objectives(
      casadi::NlpBuilder& prob,
      mimir::algorithm::NlpStructure& structure,
      const mimir::algorithm::NlpProblemBuilder& nlp)
  {
    // Adds point constraints g(x_i_0, x_i_f, ..) >= for all i in sub systems
    // Adds point objectives J += f(x_i_0, x_i_f, ..) for all i in sub systems

    // Print point expressions for x(t0) and x(tf)
    for (auto && var : structure.discrete_variable)
      BOOST_LOG_TRIVIAL(trace)
       << "Point variable: " << var.first << " is " << var.second;

    // Create input tuple: (p, v, x_1_0, x_1_f, x_2_0, x_2_f, .., x_m_0, x_m_f), i in 1,.. m.
    auto inputs = casadi::MXDict{
                                 {"p", vertcat(nlp.daeP.p)},
                                 {"v", vertcat(nlp.daeP.aux)}};
    // conventionally, terminal/cross constraints use x_1_0, x_1_f for t_1_0 and t_1_f in mapping.
    for (auto & var : structure.discrete_variable)
      inputs[var.first] = var.second;

    auto g_lb = std::vector<casadi::DM>();
    auto g_ub = std::vector<casadi::DM>();

    for (auto & func_name : nlp.constraints)
    {
      prob.g.push_back(nlp.daeP.fun(func_name) (inputs).at(func_name));
      g_lb.push_back(nlp.daeP.min(func_name));
      g_ub.push_back(nlp.daeP.max(func_name));
    }

    prob.g_lb = vertcat(casadi::DM(prob.g_lb), vertcat(g_lb)).get_elements();
    prob.g_ub = vertcat(casadi::DM(prob.g_ub), vertcat(g_ub)).get_elements();

    for (auto & func_name : nlp.objectives)
      prob.f += nlp.daeP.fun(func_name) (inputs).at(func_name);

  }

  casadi_int formulate_single_shoot(
      casadi::NlpBuilder& prob,
      mimir::algorithm::NlpStructure& structure,
      const mimir::algorithm::NlpProblemBuilder& nlp,
      const ::NlpSettings& nlp_settings,
      const std::string& sub_system,
      const casadi::Dict& opts)
  {
    BOOST_LOG_TRIVIAL(trace) << "formulate_single_shoot";
    typedef std::vector<std::string> vstr;
    typedef std::vector<double> vdob;
    auto inf = std::numeric_limits<double>::infinity();
    using namespace casadi;
    casadi_int sub_problem_size = 0;
    const casadi::DaeBuilder& dae = nlp.sub_system.at(sub_system);
    uint32_t N = nlp_settings.sub_system.at(sub_system).elements;     // Number of collocation elements
    double h_elem = nlp_settings.sub_system.at(sub_system).element_length; // Time duration of one element

    bool regular_intervals = nlp_settings.sub_system.at(sub_system).element_regular_intervals;
    bool flexible_horizon = nlp_settings.sub_system.at(sub_system).flexible_horizon;
    MX DT;
    if (flexible_horizon || !regular_intervals)
    {
      try{
        DT = nlp.daeP.var(sub_system + "_DT");
      }
      catch(casadi::CasadiException &) {
        BOOST_LOG_TRIVIAL(error) << "With flexible horizon or non-regular intervals, '"
                                 << sub_system
                                 << "_DT' must be defined in the formulation";
        throw;
      }
    }
    else
      DT = MX(h_elem);


    uint32_t nx = vertcat(dae.ode).nnz();   // Dimension of all ODE states.
    uint32_t nu = vertcat(dae.u).nnz();     // Dimension of all inputs.

    auto U = MX::sym(sub_system + "_U", nu, N);    // Inputs for all elements
    prob.x.push_back(vec(U));                      // Add vectorized U to decision variables of NLP
    sub_problem_size += nu*N;

    casadi_int idx0 = structure.problem_dimension + sub_problem_size;
    structure.variable_slice[sub_system + "_0"] = casadi::Slice(idx0, idx0 + casadi_int(nx));
    auto X = MX::sym(sub_system + "_X0", nx);      // Initial condition
    prob.x.push_back(X);
    sub_problem_size += nx;
    structure.discrete_variable[sub_system + "_0"] = X;

    auto lbg = std::vector<DM>();           // Lower bound constraint expression
    auto ubg = std::vector<DM>();           // Upper bound constraint expression

    auto V = prob.x[0];
    auto x = vertcat(dae.x);
    auto u = vertcat(dae.u);
    auto v = vertcat(nlp.daeP.aux);
    auto p = vertcat(nlp.daeP.p);
    auto t0 = dae.var("t0");

    MXDict ode_prob
     {
      {"x", vertcat(dae.x)},
      {"p", vertcat(u, v, p, t0)},
      {"t", dae.t},
      {"ode", vertcat(dae.ode)},
      {"quad", vertcat(dae.quad)}
     };

    auto grator_opts = opts.at("integrator_options").as_dict();
    auto solver = opts.at("integrator_name").as_string();

    if ((flexible_horizon || !regular_intervals) && !::is_custom_integrator(solver)) {
      throw std::runtime_error("With flexible horizon or non-regular intervals you must choose a custom integrator: rk4|heun");
    }

    casadi::Function F_int;
    if (::is_custom_integrator(solver)) {
      F_int = ::custom_integrator(DT, "F_integrator", solver, ode_prob, grator_opts);
    }
    else {
      grator_opts["tf"] = h_elem; // t0 = 0, tf = h.
      F_int = casadi::integrator("F_integrator", solver, ode_prob, grator_opts);
    }

    // Combine all output expressions into one function.
    auto F_y = Function("F_y",
     {x, u, v, p, dae.t, t0},
     {vertcat(dae.ydef)},
     vstr{"x","u","v","p","t","t0"}, vstr{"out"});
    bool skip_y = vertcat(dae.ydef).nnz() == 0;

    auto Xk = X;

    MX DTk = DT;
    for (uint32_t k = 0; k < N; ++k)
    {
      if (!regular_intervals)
        DTk = DT(k);

      MXDict int_in;
      int_in["x0"] = Xk;
      int_in["p"] = vertcat(V, U(Slice(), k), p, t0 + DTk*k);
      if (DT.is_symbolic())
        int_in["DT"] = DTk;

      auto F_int_k =  F_int(int_in);
      Xk = F_int_k["xf"];
      if (F_int_k["qf"].nnz() > 0)
        prob.f += F_int_k["qf"];

      // Add state constraints at "seams" if lower and/or upper bounds are not infinity
      int32_t indx = 0;
      for (auto && e: dae.x)
      {
        if (dae.min(e.name()) > -inf || dae.max(e.name()) < inf)
        {
          prob.g.push_back(Xk(Slice(indx, indx + e.nnz())));
          lbg.push_back(DM(vdob(e.nnz(), dae.min(e.name())))); // add lower bound for state
          ubg.push_back(DM(vdob(e.nnz(), dae.max(e.name())))); // add upper bound for state
        }
        indx += e.nnz();
      }

      // Output expressions to be satisfied, path constraints
      if(!skip_y)
      {
        auto F_y_k =
         F_y(MXDict{{"x", Xk},
                    {"u", U(Slice(), k)},
                    {"v", V},
                    {"p", p},
                    {"t", DTk*k},
                    {"t0", t0}});

        std::vector<double> low_bound_y, up_bound_y, guess_y;

        prob.g.push_back(F_y_k["out"]);
        make_bounds(dae, dae.y, low_bound_y, up_bound_y, guess_y);
        lbg.push_back(DM(low_bound_y));
        ubg.push_back(DM(up_bound_y));
      }
    }

    // X(tf)
    structure.discrete_variable[sub_system + "_f"] = Xk;

    prob.g_lb = vertcat(DM(prob.g_lb), vertcat(lbg)).get_elements();
    prob.g_ub = vertcat(DM(prob.g_ub), vertcat(ubg)).get_elements();
    ::add_nlp_variable_bounds(prob, dae, N, -1);

    BOOST_LOG_TRIVIAL(trace) << "Single shoot NLP: "
                             << "Variables: " << sub_problem_size
                             << "\n" << prob.get_str(true);

    return sub_problem_size;
  }

  casadi_int formulate_multi_shoot(
      casadi::NlpBuilder& prob,
      mimir::algorithm::NlpStructure& structure,
      const mimir::algorithm::NlpProblemBuilder& nlp,
      const ::NlpSettings& nlp_settings,
      const std::string& sub_system,
      const casadi::Dict& opts)
  {
    BOOST_LOG_TRIVIAL(debug) << "formulate_multi_shoot";
    typedef std::vector<std::string> vstr;
    typedef std::vector<double> vdob;
    using namespace casadi;
    casadi_int sub_problem_size = 0;
    const casadi::DaeBuilder& dae = nlp.sub_system.at(sub_system);
    uint32_t N = nlp_settings.sub_system.at(sub_system).elements;     // Number of collocation elements
    double h_elem = nlp_settings.sub_system.at(sub_system).element_length; // Time duration of one element

    bool regular_intervals = nlp_settings.sub_system.at(sub_system).element_regular_intervals;
    bool flexible_horizon = nlp_settings.sub_system.at(sub_system).flexible_horizon;
    MX DT;
    if (flexible_horizon || !regular_intervals)
    {
      try{
        DT = nlp.daeP.var(sub_system + "_DT");
      }
      catch(casadi::CasadiException &) {
        BOOST_LOG_TRIVIAL(error) << "With flexible horizon or non-regular intervals, '"
                                 << sub_system
                                 << "_DT' must be defined in the formulation";
        throw;
      }
    }
    else
      DT = MX(h_elem);

    uint32_t nx = vertcat(dae.ode).nnz();   // Dimension of all ODE states.
    uint32_t nu = vertcat(dae.u).nnz();     // Dimension of all inputs.

    auto U = MX::sym(sub_system + "_U", nu, N);           // Inputs for all elements
    prob.x.push_back(vec(U));               // Add vectorized U to decision variables of NLP
    sub_problem_size += nu*N;

    casadi_int idx0 = structure.problem_dimension + sub_problem_size;
    structure.variable_slice[sub_system + "_0"] =
     casadi::Slice(idx0, idx0 + casadi_int(nx));
    auto X = MX::sym(sub_system + "_X0", nx, 1 + N);       // State variables
    prob.x.push_back(vec(X));
    sub_problem_size += nx*(1 + N);
    structure.discrete_variable[sub_system + "_0"] = X(Slice(), 0);

    auto lbg = std::vector<DM>();           // Lower bound constraint expression
    auto ubg = std::vector<DM>();           // Upper bound constraint expression

    auto t0 = dae.var("t0");
    auto V = prob.x[0];
    auto x = vertcat(dae.x);
    auto u = vertcat(dae.u);
    auto v = vertcat(nlp.daeP.aux);
    auto p = vertcat(nlp.daeP.p);

    MXDict ode_prob
     {
      {"x", vertcat(dae.x)},
      {"p", vertcat(u, v, p, t0)},
      {"t", dae.t},
      {"ode", vertcat(dae.ode)},
      {"quad", vertcat(dae.quad)}
     };

    auto grator_opts = opts.at("integrator_options").as_dict();
    auto solver = opts.at("integrator_name").as_string();

    if ((flexible_horizon || !regular_intervals) && !::is_custom_integrator(solver)) {
      throw std::runtime_error("With flexible horizon or non-regular intervals you must choose a custom integrator: rk4|heun");
    }

    casadi::Function F_int;
    if (::is_custom_integrator(solver)) {
      F_int = ::custom_integrator(DT, "F", solver, ode_prob, grator_opts);
    }
    else {
      grator_opts["tf"] = h_elem; // t0 = 0, tf = h.
      F_int = casadi::integrator("F", solver, ode_prob, grator_opts);
    }
    // Combine all output expressions into one function.
    auto F_y = Function("F_y",
     {x, u, v, p, dae.t, t0},
     {vertcat(dae.ydef)},
     vstr{"x","u","v","p","t","t0"}, vstr{"out"});
    bool skip_y = vertcat(dae.ydef).nnz() == 0;

    MX p_aug;
    MX Xk, Xkp1;
    MX DTk = DT;
    // Add element constraints
    for (uint32_t k = 0; k < N; ++k)
    {
      if (!regular_intervals)
        DTk = DT(k);

      Xk = X(Slice(), k);
      MXDict int_in;
      int_in["x0"] =  X(Slice(), k);;
      int_in["p"] = vertcat(V, U(Slice(), k), p, t0 + DTk*k);
      if (DT.is_symbolic())
        int_in["DT"] = DTk;

      auto F_int_k = F_int(int_in);
      auto Xk_end = F_int_k["xf"];
      if (F_int_k["qf"].nnz() > 0)
        prob.f += F_int_k["qf"];

      // Add equality constraints at element seams
      Xkp1 = X(Slice(), k + 1);
      prob.g.push_back(Xk_end - Xkp1);
      lbg.push_back(DM::zeros(nx, 1));
      ubg.push_back(DM::zeros(nx, 1));

      // Output expressions to be satisfied
      if(!skip_y)
      {
        auto F_y_k =
         F_y(MXDict{{"x", Xk},
                    {"u", U(Slice(), k)},
                    {"v", V},
                    {"p", p},
                    {"t", DTk*k},
                    {"t0", t0}});

        vdob low_bound_y, up_bound_y, guess_y;

        prob.g.push_back(F_y_k["out"]);
        make_bounds(dae, dae.y, low_bound_y, up_bound_y, guess_y);
        lbg.push_back(DM(low_bound_y));
        ubg.push_back(DM(up_bound_y));
      }

    }

    // X(tf)
    structure.discrete_variable[sub_system + "_f"] = Xkp1; // X(:,N)

    prob.g_lb = vertcat(DM(prob.g_lb), vertcat(lbg)).get_elements();
    prob.g_ub = vertcat(DM(prob.g_ub), vertcat(ubg)).get_elements();
    ::add_nlp_variable_bounds(prob, dae, N, 0);

    BOOST_LOG_TRIVIAL(trace) << "Multi shoot NLP: "
                             << "Variables: " << sub_problem_size
                             << "\n" << prob.get_str(true);
    return sub_problem_size;
  }

  casadi_int formulate_collocation(
      casadi::NlpBuilder& prob,
      mimir::algorithm::NlpStructure& structure,
      const mimir::algorithm::NlpProblemBuilder& nlp,
      const ::NlpSettings& nlp_settings,
      const std::string& sub_system,
      const casadi::Dict& opts)
  {
    BOOST_LOG_TRIVIAL(debug) << "formulate_collocation";
    typedef std::vector<std::string> vstr;
    using namespace casadi;
    casadi_int sub_problem_size = 0;
    const casadi::DaeBuilder& dae = nlp.sub_system.at(sub_system);
    int32_t d = opts.at("degree");          // Degree of collocation polynomial
    uint32_t N = nlp_settings.sub_system.at(sub_system).elements;     // Number of collocation elements
    double h_elem = nlp_settings.sub_system.at(sub_system).element_length; // Time duration of one element

    bool regular_intervals = nlp_settings.sub_system.at(sub_system).element_regular_intervals;
    bool flexible_horizon = nlp_settings.sub_system.at(sub_system).flexible_horizon;
    MX DT;
    if (flexible_horizon || !regular_intervals)
    {
      try{
        DT = nlp.daeP.var(sub_system + "_DT");
      }
      catch(casadi::CasadiException &) {
        BOOST_LOG_TRIVIAL(error) << "With flexible horizon or non-regular intervals, '"
                                 << sub_system
                                 << "_DT' must be defined in the formulation";
        throw;
      }
    }
    else
      DT = MX(h_elem);


    uint32_t nx = vertcat(dae.ode).nnz();   // Dimension of all ODE states.
    uint32_t nu = vertcat(dae.u).nnz();     // Dimension of all inputs.

    auto U = MX::sym(sub_system + "_U", nu, N);           // Inputs for all elements
    prob.x.push_back(vec(U));               // Add vectorized U to decision variables of NLP
    sub_problem_size += nu*N;

    auto lbg = std::vector<DM>();           // Lower bound constraint expression
    auto ubg = std::vector<DM>();           // Upper bound constraint expression

    DM B,C,D;
    // Collocation points, with d=degree
    // C[d+1 x d], for slope (derivative) at collocation points
    // D[d+1 x 1], coefficients for integral at end of element
    // B[d, 1], coefficients for quadrature
    auto tau = collocation_points(d, opts.at("quadrature"));
    collocation_coeff(tau, C, D, B); // See casadi/core/integration_tools.hpp

    if(opts.at("quadrature") == "radau")
      BOOST_LOG_TRIVIAL(warning)
       << "The collocation formulation for 'radau' is flawed."
       << " A collocation point is at tau=1, duplicating the declared element seam variables."
       << " Advice to use 'legendre' for now";

    // Initial condition
    casadi_int idx0 = structure.problem_dimension + sub_problem_size;
    structure.variable_slice[sub_system + "_0"] =
     casadi::Slice(idx0, idx0 + casadi_int(nx));
    auto Xk = MX::sym(sub_system + "_X0", nx);
    prob.x.push_back(Xk);
    sub_problem_size += nx;
    structure.discrete_variable[sub_system + "_0"] = Xk;

    auto t0 = dae.var("t0");
    auto V = prob.x[0];
    auto x = vertcat(dae.x);
    auto u = vertcat(dae.u);
    auto v = vertcat(nlp.daeP.aux);
    auto p = vertcat(nlp.daeP.p);
    // ODE and quadrature
    // f(x, u, v, p, z, q, t, t0) -> [dot_x, dot_q]

    MXDict mappings;
    auto inputs = vstr{"x", "u", "v", "p", "t", "t0"}; // no z or q yet
    mappings["t"] = dae.t;
    mappings["t0"] = t0;
    mappings["x"] = x;
    mappings["u"] = u;
    mappings["v"] = v;
    mappings["p"] = p;


    auto outputs = vstr{"dot_x", "dot_q"};
    mappings["dot_x"] = vertcat(dae.ode);
    mappings["dot_q"] = vertcat(dae.quad); // currently assume quad(x,u,p,t), i.e no q

    auto ode_quad = Function("ode_quad", mappings, inputs, outputs);

    // Combine all output expressions into one function.
    auto F_y = Function("F_y",
     {x, u, v, p, dae.t, t0},
     {vertcat(dae.ydef)},
     vstr{"x","u","v","p","t","t0"}, vstr{"out"});
    bool skip_y = vertcat(dae.ydef).nnz() == 0;

    MX Z, pidot_h, Xk_end;
    MXDict res;

    auto t = mappings["t"];

    uint32_t colloc_constr = nx*d; // number of collocation constraints

    MX DTk = DT;
    // For each collocation element add dynamics and continuity constraints
    for (int32_t k = 0; k < static_cast<int64_t>(N); ++k)
    {
      if (!regular_intervals)
        DTk = DT(k);
      // Declare collocation points for element
      auto Xc = MX::sym(sub_system + "_Xc", nx, d);
      prob.x.push_back(vec(Xc));
      sub_problem_size += nx*d;

      // Perform map evaluation of ODE
      MXDict element_inputs
       {
        {"x", Xc},
        {"u", U(Slice(), Slice(k, k+1))},
        {"v", V},
        {"p", p},
        {"t0", t0},
        {"t", DTk*(MX(tau).T() + k)}};

      res = ode_quad(element_inputs);
      if (res["dot_q"].nnz() > 0)
        prob.f += mtimes(res["dot_q"], B)*DTk; // assumes quad has dim 1 and is objective

      Z = horzcat(Xk, Xc);  ///< Interpolating points of collocation polynomial
      pidot_h = mtimes(Z, C); ///< Get slope of polynomial, ZC/h = pidot

      // Polynomial slopes shall match ode
      prob.g.push_back(vec(DTk*res["dot_x"] - pidot_h)); // note the DTk
      lbg.push_back(DM::zeros(colloc_constr, 1));
      ubg.push_back(DM::zeros(colloc_constr, 1));

      // Add element end state based on collocation points numeric integral
      Xk_end = mtimes(Z, D);

      // NLP variable for state at end of element, start of the next
      Xk = MX::sym(sub_system + "_Xf", nx);
      prob.x.push_back(Xk);
      sub_problem_size += nx;

      // Equality constraint for element seams.
      prob.g.push_back(Xk_end - Xk);
      lbg.push_back(DM::zeros(nx, 1));
      ubg.push_back(DM::zeros(nx, 1));

      // Output expressions to be satisfied
      if(!skip_y)
      {
        auto F_y_k = F_y(element_inputs);
        std::vector<double> low_bound_y, up_bound_y, guess_y;
        prob.g.push_back(vec(F_y_k["out"]));                      // out is [ny x d]
        make_bounds(dae, dae.y, low_bound_y, up_bound_y, guess_y);
        lbg.push_back(kron(DM::ones(d,1), DM(low_bound_y)));      // Repeat for collocations
        ubg.push_back(kron(DM::ones(d,1), DM(up_bound_y)));
      }
    }

    structure.discrete_variable[sub_system + "_f"] = Xk;

    prob.g_lb = vertcat(DM(prob.g_lb), vertcat(lbg)).get_elements();
    prob.g_ub = vertcat(DM(prob.g_ub), vertcat(ubg)).get_elements();
    ::add_nlp_variable_bounds(prob, dae, N, d);

    BOOST_LOG_TRIVIAL(trace) << "Colloc NLP: "
                             << "Variables: " << sub_problem_size
                             << "\n" << prob.get_str(true);

    return sub_problem_size;
  }

  casadi::Function create_system_extractor(
      const mimir::algorithm::NlpProblemBuilder& nlp,
      const ::NlpSettings& nlp_settings,
      const std::string& sub_system,
      casadi_int index_start,
      casadi_int nlp_dimension)
  {
    BOOST_LOG_TRIVIAL(debug) << "create_system_extractor";
    // Defines a function that extracts w_i sub system vector from the whole NLP problem w.
    // w = col(v, w_1, ... , w_m), with w_i being the vectorized discretization of sub_system i.
    // f(w) = [0 I 0]w = w_i

    using namespace casadi;
    typedef std::vector<std::string> vstr;

    const DaeBuilder& dae = nlp.sub_system.at(sub_system);
    casadi_int N = nlp_settings.sub_system.at(sub_system).elements;
    // d>0: collocation, d=0: multi-shoot, d=-1: single-shoot
    casadi_int d = nlp_settings.sub_system.at(sub_system).degree;
    casadi_int nu = vertcat(dae.u).nnz();
    casadi_int nx = vertcat(dae.x).nnz();

    casadi_int nw_i = nu*N + nx + N*nx*(d+1);

    BOOST_LOG_TRIVIAL(trace) << "index_start: " << index_start;
    BOOST_LOG_TRIVIAL(trace) << "nw_i: " << nw_i;
    BOOST_LOG_TRIVIAL(trace) << "the_rest: " << nlp_dimension - nw_i - index_start;

    auto Av = SX(nw_i, index_start);
    auto Aw_i = SX::eye(nw_i);
    auto Aw_rest = SX(nw_i, nlp_dimension - nw_i - index_start);

    auto A = horzcat(Av, Aw_i, Aw_rest); // [0 I 0]
    auto w = SX::sym("w", nlp_dimension);

    return Function("sub_system_extractor", {w}, {mtimes(A, w)}, vstr{"w"}, vstr{"sub_out"});

  }

  casadi::Function create_horizon_shifter(
      const mimir::algorithm::NlpProblemBuilder& nlp,
      const ::NlpSettings& nlp_settings,
      const std::string& sub_system)
  {
    BOOST_LOG_TRIVIAL(debug) << "create_horizon_shifter";
    // Defines a function that takes w_i (NLP variables of sub system i) and shifts it one "element" forward.
    // The output of a function evaluation is suitable as warm start initial guess for w and lam_w.
    // w = vec(W), W = [U, X0, Xc1, ... XcN], U: nu x N, X0: nx x 1, Xci: nx x (d+1), d \in Z : d>= -1

    using namespace casadi;
    typedef std::vector<casadi_int> ivec;
    typedef std::vector<std::string> vstr;

    const casadi::DaeBuilder& dae = nlp.sub_system.at(sub_system);
    casadi_int N = nlp_settings.sub_system.at(sub_system).elements;
    // d>0: collocation, d=0: multi-shoot, d=-1: single-shoot
    casadi_int d = nlp_settings.sub_system.at(sub_system).degree;
    casadi_int nu = vertcat(dae.u).nnz();
    casadi_int nx = vertcat(dae.x).nnz();

    auto finalRight = Sparsity::triplet(N,N, ivec{N-1}, ivec{N-1}); ///< Keep last column when shifting
    /// Shifts matrix column one to the left (band(N,1) is sub diagonal, we want superdiagonal, hence transpose
    auto Ashift = SX(Sparsity::band(N,1).T() + finalRight);
    auto ushift = kron(Ashift, SX::eye(nu));                        ///< Shift operator for vec(U)
    auto x0keep = SX::eye(nx);                                      ///< No-op shift operator for vec(X0)
    auto xshift = kron(Ashift, SX::eye(nx*(d+1)));                  ///< Shift operator for vec(Xc)

    auto A = diagcat(ushift, x0keep, xshift);
    auto w = SX::sym("w", nu*N + nx + N*nx*(d+1));

    // Suggest restructure to [X0, U_1 Xc1, ... U_N XcN], shift will be easier.
    // Reason: better sparsity structure in nlp and easier to extend with more variables.
    // Note that you will then need to change add_nlp_variable_bounds.

    return casadi::Function("nlp_shifter", {w}, {mtimes(A, w)}, vstr{"w"}, vstr{"w_out"});
  }

  casadi::Function create_trajectory_function(
      const mimir::algorithm::NlpProblemBuilder& nlp,
      const ::NlpSettings& nlp_settings,
      const std::string& sub_system,
      const YAML::Node& config_schema)
  {
    BOOST_LOG_TRIVIAL(debug) << "create_trajectory_function";
    // Trajectory function (similar to single shooting). It creates expressions for a
    // solution at regular intervals within each element over the whole prediction
    // horizon, D = d+1 checkpoints for an element, N elements, total 1 + N*D elements.
    // [U, X0, p, v, t0] -> [T, [X0, X_11, .. X_1D, ... X_N1, .. X_ND]], T is time grid of X

    using namespace casadi;
    typedef std::vector<std::string> vstr;

    YAML::Node config = config_schema["config"];
    YAML::Node schema = config_schema["schema"];

    const casadi::DaeBuilder& dae = nlp.sub_system.at(sub_system);
    uint32_t N = nlp_settings.sub_system.at(sub_system).elements;          // Number of elements
    double h_elem = nlp_settings.sub_system.at(sub_system).element_length; // Time duration of one element
    uint32_t d = nlp_settings.sub_system.at(sub_system).checkpoints;       // Checkpoints within an element
    uint32_t Nd = N*(d+1);                                                 // Number of integrations

    if (!nlp_settings.sub_system.at(sub_system).element_regular_intervals)
      BOOST_LOG_TRIVIAL(warning) << "Trajectory function used on sub system with variable length!";

    auto t0 = dae.var("t0");
    MXDict ode =
     {
      {"x",   vertcat(dae.x)},
      {"p",   vertcat(vertcat(dae.u), vertcat(nlp.daeP.p), vertcat(nlp.daeP.aux), t0)},
      {"t",   dae.t},
      {"ode", vertcat(dae.ode)}
     };

    Dict config_defaults, config_opts, config_user;

    uint32_t nx = vertcat(dae.ode).nnz();   // Dimension of all ODE states.
    uint32_t nu = vertcat(dae.u).nnz();     // Dimension of all inputs.

    bool regular_intervals = nlp_settings.sub_system.at(sub_system).element_regular_intervals;
    bool flexible_horizon = nlp_settings.sub_system.at(sub_system).flexible_horizon;
    MX DT;
    if (flexible_horizon || !regular_intervals)
    {
      try{
        DT = nlp.daeP.var(sub_system + "_DT");
      }
      catch(casadi::CasadiException &) {
        BOOST_LOG_TRIVIAL(error)
         << "With flexible horizon or non-regular intervals, '"
         << sub_system << "_DT' must be defined in the formulation";
        throw;
      }
    }
    else
      DT = MX(h_elem);

    // Set integrator options
    config_user = mimir::program::parse_config(
        config["settings"]["integrator"]["options"],
        schema["map"]["settings"]["map"]["integrator"]["map"]["options"]);
    auto grator_solver = config["settings"]["integrator"]["name"].as<std::string>();

    casadi::Function elem_grator;
    if ((flexible_horizon || !regular_intervals) && !::is_custom_integrator(grator_solver)) {

        throw std::runtime_error("settings.integrator: You must choose a custom integrator: rk4|heun");
        elem_grator = ::custom_integrator(DT, "element_trajectory", grator_solver, ode, config_user);
    }
    else if (::is_custom_integrator(grator_solver)) {
      elem_grator = ::custom_integrator(DT, "element_trajectory", grator_solver, ode, config_user);
    }
    else {
      config_defaults["output_t0"] = false;
      config_opts = casadi::combine(config_user, config_defaults, true);
      elem_grator = casadi::integrator("element_trajectory", grator_solver, ode, config_opts);
    }

    auto Xki = MX::sym("X0", nx);
    auto U = MX::sym("U", nu, N);

    auto X = std::vector<MX>();
    X.push_back(Xki);
    auto T = std::vector<MX>();
    T.push_back(t0);

    MX DTk = DT;
    MX tk = t0;
    uint32_t k = 0;   // Element counter (N times)
    uint32_t ki = 0;  // Intra-element counter (d times)

    for (uint32_t kd = 0; kd < Nd; ++kd)
    {
      ki = kd % (d + 1);
      if (ki == 0 && kd > 0)
        k++;
      if (!regular_intervals)
        DTk = DT(k);
      if (ki == 0 && kd > 0)
        tk += DTk;

      MX deltaT = DTk/(d + 1);
      MX tki = tk + deltaT*(ki + 1);
      auto int_in = MXDict{{"x0", Xki},
                           {"p",
                            vertcat(
                                U(Slice(), k),
                                vertcat(nlp.daeP.p),
                                vertcat(nlp.daeP.aux),
                                tki )}};
      if (DT.is_symbolic())
        int_in["DT"] = deltaT;
      auto Fki = elem_grator(int_in);

      Xki = Fki["xf"](Slice(), Fki["xk"].size2() - 1);
      X.push_back(Xki);
      T.push_back(tki);
    }

    return casadi::Function(
        "trajectory_" + sub_system,
        {U, X[0], vertcat(nlp.daeP.p), vertcat(nlp.daeP.aux), t0},
        {vertcat(T), horzcat(X)},
        vstr{"U", "X0", "p", "v", "t0"},
        vstr{"T", "X"});
  }

  casadi::Function create_solution_unpacker(
      const mimir::algorithm::NlpProblemBuilder& nlp,
      const ::NlpSettings& nlp_settings,
      const std::string& sub_system)
  {
    BOOST_LOG_TRIVIAL(debug) << "create_solution_unpacker";
    // Unpacks the nlp decision variable vector w.
    // Let N be elements, degree: d >= -1, u: nu x 1, x: nx x 1
    // w: N*nu + nx + N*nx*(d+1), U: nu x N, X0: nx x 1, Xc: nx x N*(d+1)
    // [w] -> [U, X0, Xc]

    using namespace casadi;
    typedef std::vector<std::string> vstr;

    const casadi::DaeBuilder& dae = nlp.sub_system.at(sub_system);
    casadi_int N = nlp_settings.sub_system.at(sub_system).elements;
    casadi_int nu = vertcat(dae.u).nnz();
    casadi_int nx = vertcat(dae.x).nnz();
    casadi_int d = nlp_settings.sub_system.at(sub_system).degree;
    // d>0: collocation, d=0: multi-shoot, d=-1: single-shoot
    casadi_int Q = N*nu + nx*(1 + N*(d+1)); // Dimension of w
    casadi_int dim_var;

    auto w = SX::sym("w", Q);

    dim_var = N*nu;
    auto A_U = Sparsity::diag(dim_var);
    A_U.appendColumns(Sparsity(dim_var, Q - dim_var));

    auto U = reshape(mtimes(SX(A_U),w), nu, N);

    dim_var = nx;
    auto A_X0 = Sparsity(dim_var, N*nu);
    A_X0.appendColumns(Sparsity::diag(dim_var));
    A_X0.appendColumns(Sparsity(dim_var, Q - N*nu - dim_var));

    auto X0 = reshape(mtimes(SX(A_X0),w), nx, 1);

    dim_var = nx*N*(d+1);
    auto A_Xc = Sparsity(dim_var, Q - dim_var);
    A_Xc.appendColumns(Sparsity::diag(dim_var));

    auto Xc = reshape(mtimes(SX(A_Xc),w), nx, N*(d+1));

    return casadi::Function("nlp_unpacker", {w}, {U, X0, Xc}, vstr{"w"}, vstr{"U", "X0", "Xc"});

  }

  casadi::Function create_decision_parameter_extractor(
      const mimir::algorithm::NlpProblemBuilder& nlp,
      const mimir::algorithm::NlpStructure& nlp_structure)
  {
    BOOST_LOG_TRIVIAL(debug) << "create_decision_parameter_extractor";
    // Fetches the decision parameter vector v from w.
    // Assumes that v is placed first in w = [v, s_1, s_2, ..]
    // [w] -> [v]

    using namespace casadi;
    typedef std::vector<std::string> vstr;

    casadi_int dim_w = nlp_structure.problem_dimension;
    casadi_int dim_v = vertcat(nlp.daeP.aux).nnz();

    auto w = SX::sym("w", dim_w);
    auto A_v = Sparsity::diag(dim_v); // I: n_v x n_v
    A_v.appendColumns(Sparsity(dim_v, dim_w - dim_v)); // [I 0] n_v x n_w
    auto v = mtimes(SX(A_v),w);
    return casadi::Function("decision_parameters", {w}, {v}, vstr{"w"}, vstr{"v"});
  }

  casadi::Function create_solution_timegrid(
      const mimir::algorithm::NlpProblemBuilder& nlp,
      const ::NlpSettings& nlp_settings,
      const std::string& sub_system)
  {
    BOOST_LOG_TRIVIAL(debug) << "create_solution_timegrid";
    // Creates time grid vectors for the state vectors in the nlp
    // Let d: degree >= -1. T is time at all X, Tu is time points for U
    // [v, t0] -> [T, Tu], T: (1 + N*(d+1)) x 1, Tu: N x 1
    // It depends on the decision vector due to possible flexible horizon/non-regular time grid.
    using casadi::MX;
    using casadi::Slice;
    typedef std::vector<std::string> vstr;

    uint32_t N = nlp_settings.sub_system.at(sub_system).elements;          // Number of elements
    double h_elem = nlp_settings.sub_system.at(sub_system).element_length; // Time duration of one element
    int32_t d = nlp_settings.sub_system.at(sub_system).degree;             // -1, 0, d, single, multi, collocation
    std::string technique = nlp_settings.sub_system.at(sub_system).technique;
    std::string colloc_quadrature = nlp_settings.sub_system.at(sub_system).quadrature;

    bool regular_intervals = nlp_settings.sub_system.at(sub_system).element_regular_intervals;
    bool flexible_horizon = nlp_settings.sub_system.at(sub_system).flexible_horizon;
    MX DT;
    if (flexible_horizon || !regular_intervals)
    {
      try{
        DT = nlp.daeP.var(sub_system + "_DT");
      }
      catch(casadi::CasadiException &) {
        BOOST_LOG_TRIVIAL(error)
         << "With flexible horizon or non-regular intervals, '"
         << sub_system << "_DT' must be defined in the formulation";
        throw;
      }
    }
    else
      DT = MX(h_elem);

    auto t0 = MX::sym("t0");
    auto T = MX::zeros(1 + N*(d + 1), 1);
    auto Tu = MX::zeros(N, 1);
    auto v = vertcat(nlp.daeP.aux);
    T(0) = t0;
    Tu(0) = t0;
    MX DTk = DT;
    for(casadi_int i = 1; i < N; ++i)
    {
      if (!regular_intervals)
        DTk = DT(i-1);
      Tu(i) = Tu(i-1) + DTk;
    }
    MX tk = t0;
    if (technique == "collocation")
    {
      auto tau = casadi::collocation_points(d, colloc_quadrature);
      tau.push_back(1);
      for(casadi_int i = 0; i < N; ++i){
        if (!regular_intervals)
          DTk = DT(i);
        T(Slice(1+i*(d+1),1+(i+1)*(d+1))) = tk + DTk*MX(tau);
        tk += DTk;
      }
    }
    else if (technique == "single-shooting")
      ; // No-op
    else if (technique == "multi-shooting")
      for(casadi_int i = 1; i <= N; ++i){
        if (!regular_intervals)
          DTk = DT(i-1);
        T(i) = T(i-1) + DTk;
      }
    else
      throw std::runtime_error(
          std::string("Unable to create time grid for unknown technique: " + technique));

    return casadi::Function(
        "timegrid_" + sub_system, {t0, v}, {T, Tu}, vstr{"t0", "v"}, vstr{"T", "Tu"});
  }

  void create_helper_functions(
      const ::NlpSettings& nlp_set,
      const mimir::algorithm::NlpProblemBuilder& nlp_builder,
      const mimir::algorithm::NlpStructure& nlp_structure,
      const YAML::Node& config_schema,
      std::map<std::string, std::map<std::string, casadi::Function>>& mpc)
  {
    BOOST_LOG_TRIVIAL(debug) << "create_helper_functions";
    // Instantiate functions for interacting with NLP formulation
    typedef std::vector<std::string> vstr;

    casadi::SX w = casadi::SX::sym("w", nlp_structure.problem_dimension);

    // Function to extract decision vector
    mpc["v"]["extractor"] = ::create_decision_parameter_extractor(nlp_builder, nlp_structure);
    BOOST_LOG_TRIVIAL(debug) << mpc.at("v").at("extractor").get_str();

    // For each sub system add helper functions to MPC struct
    for (auto && sub_sys : nlp_builder.sub_system_names)
    {
      auto extractor = ::create_system_extractor(
          nlp_builder,
          nlp_set,
          sub_sys,
          nlp_structure.variable_slice.at(sub_sys).start,
          nlp_structure.problem_dimension);

      auto shifter = ::create_horizon_shifter(nlp_builder, nlp_set, sub_sys);
      casadi::SX shifter_expr = shifter({extractor({w})[0]})[0];
      mpc[sub_sys]["shifter"] = casadi::Function(
          sub_sys + "_shifter", {w}, {shifter_expr}, shifter.name_in(), shifter.name_out());

      auto unpacker = ::create_solution_unpacker(nlp_builder, nlp_set, sub_sys);
      auto unpacker_expr = unpacker({extractor({w})});
      mpc[sub_sys]["unpacker"] = casadi::Function(
          sub_sys + "_unpacker", {w}, unpacker_expr, unpacker.name_in(), unpacker.name_out());

      mpc[sub_sys]["timegrid"] = ::create_solution_timegrid(nlp_builder, nlp_set, sub_sys);

      mpc[sub_sys]["trajectory"] = ::create_trajectory_function(
          nlp_builder, nlp_set, sub_sys, config_schema);

      auto sub_settings = nlp_set.sub_system.at(sub_sys);
      if (sub_settings.flexible_horizon && sub_settings.element_regular_intervals)
      {
        auto Tf = mpc["v"]["extractor"]({w})[0](
            nlp_builder.decision_parameter_slice.at(sub_sys+"_DT"))*sub_settings.elements;
        mpc[sub_sys]["Tf"] = casadi::Function(sub_sys + "_Tf", {w}, {Tf}, vstr{"w"}, vstr{"Tf"});
      }

      BOOST_LOG_TRIVIAL(info) << "NLP sub problem " << sub_sys
                              << " is formulated using "
                              << nlp_set.sub_system.at(sub_sys).technique;
      BOOST_LOG_TRIVIAL(debug) << mpc.at(sub_sys).at("shifter").get_str();
      BOOST_LOG_TRIVIAL(debug) << mpc.at(sub_sys).at("unpacker").get_str();
      BOOST_LOG_TRIVIAL(debug) << mpc.at(sub_sys).at("timegrid").get_str();
      BOOST_LOG_TRIVIAL(debug) << mpc.at(sub_sys).at("trajectory").get_str();

    }
  }

  void formulate_nlp(
      casadi::NlpBuilder& nlp_problem,
      casadi::Function& nlp_solver,
      std::unique_ptr<mimir::algorithm::PursePlannerFormulation::NlpCallback>& nlp_callback,
      std::map<std::string, std::map<std::string, casadi::Function>>& mpc,
      fkin::NlpConfig& nlp_configuration,
      mimir::algorithm::NlpProblemBuilder& nlp_builder,
      mimir::algorithm::NlpStructure& nlp_structure,
      const YAML::Node& config_schema)
  {
    BOOST_LOG_TRIVIAL(debug) << "formulate_nlp";
    auto nlp_config = config_schema["config"]["settings"]["nlp"];
    auto nlp_schema = config_schema["schema"]["map"]["settings"]["map"]["nlp"];
    casadi::Dict nlp_dict = mimir::program::parse_config(nlp_config, nlp_schema);

    auto discr_config = nlp_config["discretization"];
    auto discr_schema = nlp_schema["map"]["discretization"];

    // This is also set independently in formulate_dynamics..
    auto sub_names = std::vector<std::string>{"x_1", "x_2"};

    ::NlpSettings nlp_set(discr_config, sub_names);

    if (true){ // Improved formulation
      ::formulate_dynamics(nlp_builder, nlp_set, config_schema);
    }
    if (false){ // Was used for fundamental testing
      ::formulate_test_system(nlp_builder, nlp_set, config_schema);
    }
    if (false){  // Old formulation
      ::formulate_dae(nlp_builder, nlp_set, config_schema);
    }

    // Confirm that sub system names set by formulate coincides with sub_names
    auto bld_names = nlp_builder.sub_system_names;
    std::sort(sub_names.begin(), sub_names.end());
    std::sort(bld_names.begin(), bld_names.end());
    bool is_equal = (sub_names.size() == bld_names.size());
    if (is_equal)
      is_equal = std::equal (bld_names.begin(), bld_names.end(), sub_names.begin());
    if(!is_equal){
      BOOST_LOG_TRIVIAL(error) << "Sub system names are not equal: "
       << bld_names << " vs. " << sub_names;
      throw std::runtime_error("Sub-system names inconsistent in formulation");
    }

    // Objective function is 0 by default
    nlp_problem.f = casadi::MX::zeros(1,1);

    // === Decision parameters ===
    nlp_problem.x.push_back(vertcat(nlp_builder.daeP.aux));
    ::add_nlp_variable_bounds(nlp_problem, nlp_builder.daeP, 0, 0, true);
    nlp_structure.variable_slice["v"] = casadi::Slice(0, nlp_problem.x[0].nnz());
    nlp_structure.problem_dimension += nlp_problem.x[0].nnz();

    bool is_hybrid = false;
    std::string prev_tech;

    // Discretize all subsystems with indicated scheme
    for (auto && sub_system : nlp_builder.sub_system_names)
    {
      casadi::Dict system_settings = mimir::program::parse_config(
          discr_config[sub_system], discr_schema["map"][sub_system]);
      casadi::Dict technique_options = system_settings.at("options").as_dict();
      casadi_int discret_size;
      std::string technique = nlp_set.sub_system.at(sub_system).technique;
      if (technique == "collocation") {
        discret_size = ::formulate_collocation(
            nlp_problem, nlp_structure, nlp_builder, nlp_set, sub_system, technique_options);
      }
      else if (technique == "single-shooting") {
        discret_size = ::formulate_single_shoot(
            nlp_problem, nlp_structure, nlp_builder, nlp_set, sub_system, technique_options);
      }
      else if (technique == "multi-shooting") {
        discret_size = ::formulate_multi_shoot(
            nlp_problem, nlp_structure, nlp_builder, nlp_set, sub_system, technique_options);
      }
      else {
        auto supported =
         std::string(".\n  Supported techniques are: collocation, single-shooting, multi-shooting");
        throw std::runtime_error(
            std::string("Could not formulate NLP: Unknown technique: ")
            + nlp_set.sub_system.at(sub_system).technique + supported);
      }

      nlp_structure.variable_slice[sub_system] =
       casadi::Slice(
           nlp_structure.problem_dimension,
           nlp_structure.problem_dimension + discret_size);
      nlp_structure.problem_dimension += discret_size;

      if(!prev_tech.empty() && prev_tech != technique)
        is_hybrid = true;
      prev_tech = technique;
    }

    ::add_point_constraints_and_objectives(nlp_problem, nlp_structure, nlp_builder);

    std::vector<casadi::MX> t0s;
    for (auto && sub_system : nlp_builder.sub_system_names)
      t0s.push_back(nlp_builder.sub_system.at(sub_system).var("t0"));

    // Formulate NLP function with MXDict
    casadi::MXDict nlp_problem_expr
     {
      {"f", nlp_problem.f},
      {"x", vertcat(nlp_problem.x)},
      {"p", vertcat(vertcat(nlp_builder.daeP.p), vertcat(t0s))},
      {"g", vertcat(nlp_problem.g)}
     };

    casadi::Dict callback_opts;
    callback_opts["nx"] = nlp_problem_expr.at("x").nnz();
    callback_opts["ng"] = nlp_problem_expr.at("g").nnz();
    callback_opts["np"] = nlp_problem_expr.at("p").nnz();
    nlp_callback.reset(new mimir::algorithm::PursePlannerFormulation::NlpCallback(
         "PursePlanner_callback", callback_opts));

    // Fetch user options for NLP solver
    auto nlp_solver_opts = nlp_dict.at("solver").as_dict();
    auto nlp_options = nlp_solver_opts.at("options").as_dict();
    std::string solver_name = nlp_solver_opts.at("name").as_string();

    if(nlp_config["with_callback"].as<bool>())
      nlp_options["iteration_callback"] = *nlp_callback.get();

    if(!casadi::has_nlpsol(solver_name))
      throw std::runtime_error(
          std::string("casadi could not find solver: ") + solver_name);

    nlp_solver = nlpsol("nlp_purse_planner", solver_name, nlp_problem_expr, nlp_options);
    BOOST_LOG_TRIVIAL(debug) << "NLP options: " << nlp_options;
    BOOST_LOG_TRIVIAL(debug) << nlp_solver.get_str();

    // Create helper functions and add them to mpc
    ::create_helper_functions(nlp_set, nlp_builder, nlp_structure, config_schema, mpc);

    // We only publish information about a sub system of the formulation (x_2)..
    std::string sub_sys = "x_2";
    nlp_configuration.id() = config_schema["config"]["outputs"]["id"].as<std::string>();
    nlp_configuration.technique() = nlp_set.sub_system.at(sub_sys).technique;
    nlp_configuration.degree() = nlp_set.sub_system.at(sub_sys).degree;
    nlp_configuration.solver() = solver_name;
    nlp_configuration.horizon() = nlp_set.sub_system.at(sub_sys).horizon;
    nlp_configuration.elements() = nlp_set.sub_system.at(sub_sys).elements;
    nlp_configuration.nx() = nlp_problem_expr.at("x").nnz();
    nlp_configuration.ng() = nlp_problem_expr.at("g").nnz();
    nlp_configuration.np() = nlp_problem_expr.at("p").nnz();

    if (is_hybrid)
      nlp_configuration.technique() = std::string("hybrid");
  }
}


namespace mimir
{

  algorithm::PursePlannerFormulation::PursePlannerFormulation(const YAML::Node& config_schema)
  {
    ::formulate_nlp(
        nlp_problem,
        nlp_solver,
        nlp_callback,
        mpc,
        m_nlp_config,
        nlp_builder,
        m_nlp_structure,
        config_schema);

    // Temporary hack to expose it for use in finite difference cap
    omega_max = config_schema["config"]["settings"]["nlp"]["formulation"]["subsystem"]["x_1"]["omega_max"].as<double>();

    // Temporary hack to expose it for use in solution trajectory integrator with custom deltaT
    N2 = config_schema["config"]["settings"]["nlp"]["discretization"]["x_2"]["elements"]["count"].as<std::uint32_t>();


  }

  algorithm::PursePlannerFormulation::~PursePlannerFormulation() { }

  casadi::Slice algorithm::PursePlannerFormulation::get_slice(
      const std::string& name, const std::vector<casadi::MX>& vars)
  {
    casadi_int start_index = 0;
    for (auto && var : vars)
    {
      if (var.name() == name)
        return casadi::Slice(start_index, start_index + var.nnz());
      else
        start_index += var.nnz();
    }
    throw std::out_of_range("Variable not found in vector: " + name);
  }

  void algorithm::PursePlannerFormulation::set_variable(
      const std::string& name, const casadi::DM& value_in)
  {

    try {
      auto slice = m_nlp_structure.variable_slice.at(name); // Slice of variable in question
      BOOST_LOG_TRIVIAL(trace) << "Slice for " << name << " is: " << slice;

      if (value_in.nnz() != slice.stop - slice.start){
        BOOST_LOG_TRIVIAL(error)
         << "Variable: "
         << name << " with inconsistent dimensions: expected: "
         << slice.stop - slice.start << ", but got: " << value_in.nnz();
        throw std::runtime_error("set_variable with inconsistent dimensions");
      }

      bool is_init_condition = false;
      std::vector<bool> floating_init;
      if (name != "v")
      {
        is_init_condition = true;
        floating_init = nlp_builder.x_inits.at(name);
      }

      size_t index = slice.start;
      for (casadi_int i = 0; i < value_in.nnz(); ++i)
      {
        double value = double(value_in.nz(i));
        nlp_problem.x_init.at(index + i) = value;

        if (is_init_condition && !floating_init[i]){

          // This is the only place lower and upper bounds are set equal to x_init
          nlp_problem.x_lb.at(index + i) = value;
          nlp_problem.x_ub.at(index + i) = value;
        }
      }
    }
    catch (const std::out_of_range& oor) {
      BOOST_LOG_TRIVIAL(error)
       << "Tried to set unknown variable: "
       << name << ": " << oor.what() << '\n';
      throw;
    }
    catch (const std::runtime_error&) {
      throw;
    }

  }


  fkin::OptiStats algorithm::PursePlannerFormulation::stats()
  {

    // Only tested for ipopt and bonmin

    casadi::Dict nlp_stats(nlp_solver.stats());

    fkin::OptiStats stats;
    fkin::NlpFuncStat func_stat;
    stats.iterations() = static_cast<std::uint32_t>(nlp_stats.at("iter_count").as_int());
    stats.status() = nlp_stats.at("success").as_bool();
    stats.status_text() = nlp_stats.at("return_status").as_string();

    BOOST_LOG_TRIVIAL(trace) << "Callback call stats";
    func_stat.callback_fcn() = nlp_stats.at("n_call_callback_fun").as_int();
    func_stat.nlp_f()        = nlp_stats.at("n_call_nlp_f").as_int();
    func_stat.nlp_g()        = nlp_stats.at("n_call_nlp_g").as_int();
    func_stat.nlp_grad()     = nlp_stats.at("n_call_nlp_grad").as_int();
    func_stat.nlp_grad_f()   = nlp_stats.at("n_call_nlp_grad_f").as_int();
    try {
      func_stat.nlp_hess_l()   = nlp_stats.at("n_call_nlp_hess_l").as_int();
    } catch(std::out_of_range&) {
      // Hessian info no available in case of hessian_approximation
      func_stat.nlp_hess_l() = 0;
    }
    func_stat.nlp_jac_g()    = nlp_stats.at("n_call_nlp_jac_g").as_int();
    func_stat.total()        = nlp_stats.at("n_call_total").as_int();
    stats.n_call() = func_stat;

    BOOST_LOG_TRIVIAL(trace) << "Callback timing process time stats";
    func_stat.callback_fcn() = nlp_stats.at("t_proc_callback_fun").as_double();
    func_stat.nlp_f()        = nlp_stats.at("t_proc_nlp_f").as_double();
    func_stat.nlp_g()        = nlp_stats.at("t_proc_nlp_g").as_double();
    func_stat.nlp_grad()     = nlp_stats.at("t_proc_nlp_grad").as_double();
    func_stat.nlp_grad_f()   = nlp_stats.at("t_proc_nlp_grad_f").as_double();
    try {
      func_stat.nlp_hess_l()   = nlp_stats.at("t_proc_nlp_hess_l").as_int();
    } catch(std::out_of_range&) {
      // Hessian info no available in case of hessian_approximation
      func_stat.nlp_hess_l() = 0;
    }
    func_stat.nlp_jac_g()    = nlp_stats.at("t_proc_nlp_jac_g").as_double();
    func_stat.total()        = nlp_stats.at("t_proc_total").as_double();
    stats.t_proc() = func_stat;

    BOOST_LOG_TRIVIAL(trace) << "Callback timing wall time stats";
    func_stat.callback_fcn() = nlp_stats.at("t_wall_callback_fun").as_double();
    func_stat.nlp_f()        = nlp_stats.at("t_wall_nlp_f").as_double();
    func_stat.nlp_g()        = nlp_stats.at("t_wall_nlp_g").as_double();
    func_stat.nlp_grad()     = nlp_stats.at("t_wall_nlp_grad").as_double();
    func_stat.nlp_grad_f()   = nlp_stats.at("t_wall_nlp_grad_f").as_double();
    try {
      func_stat.nlp_hess_l()   = nlp_stats.at("t_wall_nlp_hess_l").as_int();
    } catch(std::out_of_range&) {
      // Hessian info no available in case of hessian_approximation
      func_stat.nlp_hess_l() = 0;
    }
    func_stat.nlp_jac_g()    = nlp_stats.at("t_wall_nlp_jac_g").as_double();
    func_stat.total()        = nlp_stats.at("t_wall_total").as_double();
    stats.t_wall() = func_stat;

    BOOST_LOG_TRIVIAL(trace) << nlp_stats;

    return stats;
  }

  void algorithm::PursePlannerFormulation::set_abort(const std::atomic<bool>& cancel_token)
  {
    nlp_callback->set_abort(cancel_token);
  }


}
