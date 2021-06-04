#include <algorithm>
#include <chrono>
#include <cmath>
#include <exception>
#include <map>
#include <string>
#include <sstream>
#include <thread>
#include <utility>
#include <vector>

#include <boost/log/trivial.hpp>
#include <casadi/casadi.hpp>
#include <yavl-cpp/yavl.h>

#include "mimir/algorithm/PursePlannerFormulation.hpp"
#include "mimir/algorithm/PursePlanner.hpp"
#include "mimir/StateMachineFwd.hpp"
#include "mimir/FkinDds.hpp"

// ReaderLifespan
#include <dds/core/detail/ProprietaryApi.hpp>
#include <dds/sub/AnyDataReader.hpp>
#include <dds/pub/AnyDataWriter.hpp>


#ifdef MIMIR_WITH_GNUPLOT
#include "gnuplot-iostream.h"
#include <tuple>
#endif

namespace mimir
{
  namespace algorithm
  {

    class PursePlanner::Impl
    {
    public:
      Impl(
          const YAML::Node& conf,
          dds::pub::Publisher publisher,
          dds::sub::Subscriber subscriber)
      {

        const YAML::Node config = conf["config"];
        const YAML::Node schema = conf["schema"];
        YAVL::Validator checkNode(schema, config);

        if(!checkNode.validate())
        {
          std::stringstream oss;
          oss << checkNode.get_errors();

          BOOST_LOG_TRIVIAL(debug) << "PursePlanner config:\n" << config;
          BOOST_LOG_TRIVIAL(debug) << "PursePlanner schema:\n" << schema;

          throw YAML::Exception(YAML::Mark(), oss.str());
        }

        // set initial condition for planner state vector
        x_0 = casadi::DM(config["initial_condition"]["x0"].as<std::vector<double>>());

        // Function for mapping gps to local coordinate system
        toNED = define_gps_to_ned();

        // === Create instance of the purse planner nlp formulation ===
        BOOST_LOG_TRIVIAL(trace) << "Instantiate formulation";
        formulation.reset(new PursePlannerFormulation(conf));

        decide_idx = formulation->nlp_builder.decision_parameter_slice;
        state1_idx = formulation->nlp_builder.state_slice.at("x_1");
        state2_idx = formulation->nlp_builder.state_slice.at("x_2");
        input1_idx = formulation->nlp_builder.input_slice.at("x_1");
        input2_idx = formulation->nlp_builder.input_slice.at("x_2");
        x1_dim = formulation->nlp_builder.state_dimension.at("x_1");
        x2_dim = formulation->nlp_builder.state_dimension.at("x_2");
        v_dim = formulation->nlp_builder.decision_parameter_dimension;

        if (x_0.nnz() != x2_dim)
        {
          BOOST_LOG_TRIVIAL(error)
           << "Mismatching dimensions: x_0 in config vs. expected from formulation: x_0:"
           << x_0.nnz() << " vs. " << x2_dim;
        }

        // === DDS subscribers for parameters defined in YAML config ===
        setup_parameters(config, subscriber);

        // === DDS subscribers for inputs defined in YAML config ===
        setup_inputs(config, subscriber);

        // === DDS publishers for outputs defined in YAML config ===
        setup_outputs(config, publisher);

        // === Enable gnuplot debug plotting ===
        plot_on = config["settings"]["plot"].as<bool>();

#ifdef MIMIR_WITH_GNUPLOT
        if(plot_on)
        {
          gnuplot //<< "set yrange [-200:200]\n"
                  << "set title 'PursePlanner'\n"
                  << "set ylabel 'Various [-]'\n"
                  << "set xlabel 'Time [s]'\n";

          ne_plot << "set title 'Map plot'\n"
                  << "set size square\n"
                  << "set size ratio -1\n"
                  << "set ylabel 'North [m]'\n"
                  << "set xlabel 'East [m]'\n";
        }
#endif

      }
      ~Impl() {}
      std::map<std::string, dds::sub::AnyDataReader> parameter_inputs;
      std::map<std::string, dds::sub::AnyDataReader> input_readers;
      std::map<std::string, dds::pub::AnyDataWriter> output_writers;
      std::map<std::string, casadi::Slice>
        decide_idx, param_idx, input1_idx, input2_idx, state1_idx, state2_idx;
      std::unique_ptr<PursePlannerFormulation> formulation;
      std::string outputId;              ///< DDS outputs identifier
      casadi::DM parameters;             ///< Parameters p
      casadi::DM x_0, x_k;               ///< State at t_0 and t_k
      double t_k;                        ///< Simulation time of x_k [seconds]
      casadi::DM T, X, Tu, U;            ///< Solution trajectories
      std::vector<double> lamb_g;        ///< Lagrange multipliers: inequality constraints
      casadi::DM gps_origin;             ///< GPS coordinate of NED origin reference frame
      casadi::DM T1, X1;                 ///< To be removed
      bool plot_on;                      ///< To plot with gnuplot or not
      casadi::Function toNED;
      casadi_int x1_dim, x2_dim, v_dim, p_dim;
#ifdef MIMIR_WITH_GNUPLOT
      Gnuplot gnuplot, ne_plot;
#endif
    private:

      casadi::Function define_gps_to_ned()
      {
        using casadi::SX;
        typedef std::vector<std::string> vstr;
        typedef std::vector<SX> SXvec;

        double PI180 = std::atan(1)*4./180.;

        SX lat = SX::sym("lat");
        SX lon = SX::sym("lon");
        SX h = SX::sym("height");

        SX mu = lat*PI180;
        SX l = lon*PI180;

        // Fossen 2002 p.42
        SX re = 6378137; // Equatorial radius of ellipsoid
        SX rp = 6356752; // Polar axis radius of ellipsoid

        SX re2 = re*re;
        SX rp2 = rp*rp;

        SX N = (re2)/(sqrt(re2*pow(cos(mu),2) + rp2*pow(sin(mu),2)));

        // ECEF
        SX xe = (N + h)*cos(mu)*cos(l);
        SX ye = (N + h)*cos(mu)*sin(l);
        SX ze = ((rp2/re2)*N + h)*sin(mu);

        SX ecef = vertcat(xe, ye, ze);

        casadi::Function toECEF("toECEF", SXvec{lat,lon,h}, SXvec{ecef});

        SX lat0 = SX::sym("lat0");
        SX lon0 = SX::sym("lon0");
        SX h0 = SX::sym("h0");

        SX cosl = cos(lon*PI180);
        SX sinl = sin(lon*PI180);
        SX cosmu = cos(lat*PI180);
        SX sinmu = sin(lat*PI180);

        SX Tne = SX::zeros(3,3);
        Tne(0,0) = -cosl*sinmu;
        Tne(0,1) = -sinl;
        Tne(0,2) = -cosl*cosmu;
        Tne(1,0) = -sinl*sinmu;
        Tne(1,1) = cosl;
        Tne(1,2) = -sinl*cosmu;
        Tne(2,0) = cosmu;
        Tne(2,2) = -sinmu;

        SX diff = toECEF(SXvec{lat,lon,0})[0] - toECEF(SXvec{lat0,lon0, 0})[0];

        auto ned = mtimes(Tne.T(), diff);

        auto ret = casadi::Function("toNED",
         {vertcat(lat, lon), vertcat(lat0, lon0)},
         {ned},
         vstr{"GPS", "GPS0"}, vstr{"NED"});

        /*
        typedef std::vector<double> vdob;
        auto res = ret(casadi::DMDict{{"GPS", casadi::DM(vdob{63.45917,    10.3683367})},
                                      {"GPS0", casadi::DM(vdob{63.4581027, 10.3683367})}});
        BOOST_LOG_TRIVIAL(info) << "Zero transform is: " << res["NED"];
        */

        return ret;

      }

      void setup_parameters(const YAML::Node& config, dds::sub::Subscriber subscriber)
      {
        BOOST_LOG_TRIVIAL(trace) << "setup parameters";
        typedef std::vector<double> vdob;
        using casadi::Slice;
        auto readerQos = subscriber.default_datareader_qos();
        const auto config_param = config["parameters"];

        p_dim = formulation->nlp_builder.parameter_dimension;
        parameters = casadi::DM::zeros(p_dim, 1);
        param_idx = formulation->nlp_builder.parameter_slice;

        parameters(param_idx.at("setting_speed_U_v"))   = config_param["setting_speed_U_v"]["default"].as<double>();
        parameters(param_idx.at("setting_radii"))       = config_param["setting_radii"]["default"].as<vdob>();
        parameters(param_idx.at("aim_distance_D_s"))    = config_param["aim_distance_D_s"]["default"].as<double>();
        parameters(param_idx.at("fish_margin_d_f"))     = config_param["fish_margin_d_f"]["default"].as<double>();
        parameters(param_idx.at("fish_depth_z_s"))      = config_param["fish_depth_z_s"]["default"].as<double>();
        parameters(param_idx.at("leadline_tau_ll_z_d")) = config_param["leadline_tau_ll_z_d"]["default"].as<vdob>();
        parameters(param_idx.at("sink_margin_z_min"))   = config_param["sink_margin_z_min"]["default"].as<double>();

        auto add_parameter_IdVec1d =
         [=](const std::string& variable){

           auto topicName = config_param[variable]["topic"].as<std::string>();
           auto id = config_param[variable]["id"].as<std::string>();
           parameter_inputs.emplace(
               variable,
               dds::sub::DataReader<fkin::IdVec1d>(
                   subscriber,
                   dds::topic::ContentFilteredTopic<fkin::IdVec1d>(
                       dds::topic::Topic<fkin::IdVec1d>(
                           subscriber.participant(),
                           topicName),
                       topicName + id,
                       dds::topic::Filter("id = %0", {id})),
                   readerQos));
         };

        auto add_parameter_IdVec2d =
         [=](const std::string& variable){

           auto topicName = config_param[variable]["topic"].as<std::string>();
           auto id = config_param[variable]["id"].as<std::string>();
           parameter_inputs.emplace(
               variable,
               dds::sub::DataReader<fkin::IdVec2d>(
                   subscriber,
                   dds::topic::ContentFilteredTopic<fkin::IdVec2d>(
                       dds::topic::Topic<fkin::IdVec2d>(
                           subscriber.participant(),
                           topicName),
                       topicName + id,
                       dds::topic::Filter("id = %0", {id})),
                   readerQos));
         };

        auto add_parameter_Double1 =
         [=](const std::string& variable){

           auto topicName = config_param[variable]["topic"].as<std::string>();
           parameter_inputs.emplace(
               variable,
               dds::sub::DataReader<ratatosk::types::DoubleVal>(
                   subscriber,
                   dds::topic::Topic<ratatosk::types::DoubleVal>(
                       subscriber.participant(),
                       topicName),
                   readerQos));
         };

        auto add_parameter_Double2 =
         [=](const std::string& variable){

           auto topicName = config_param[variable]["topic"].as<std::string>();
           parameter_inputs.emplace(
               variable,
               dds::sub::DataReader<ratatosk::types::Double2>(
                   subscriber,
                   dds::topic::Topic<ratatosk::types::Double2>(
                       subscriber.participant(),
                       topicName),
                   readerQos));
         };

        auto add_parameter_Double3 =
         [=](const std::string& variable){

           auto topicName = config_param[variable]["topic"].as<std::string>();
           parameter_inputs.emplace(
               variable,
               dds::sub::DataReader<ratatosk::types::Double3>(
                   subscriber,
                   dds::topic::Topic<ratatosk::types::Double3>(
                       subscriber.participant(),
                       topicName),
                   readerQos));
         };

        // These are subscribed parameters for which we read from
        add_parameter_IdVec1d("setting_speed_U_v");
        add_parameter_Double2("current_surface");
        add_parameter_IdVec2d("setting_radii");
        add_parameter_IdVec1d("aim_distance_D_s");
        add_parameter_IdVec1d("fish_margin_d_f");
        add_parameter_Double2("fish_velocity_over_ground");
        add_parameter_Double2("current_fish");
        add_parameter_IdVec1d("fish_depth_z_s");
        add_parameter_IdVec2d("leadline_tau_ll_z_d");
        add_parameter_IdVec1d("sink_margin_z_min");

        // some parameters are not read directly from DDS, they are calculated based on other

        auto current_surface = casadi::DM(config_param["current_surface"]["default"].as<vdob>());
        double W = double(norm_2(current_surface));
        double alpha = std::atan2(
            double(current_surface(1)),
            double(current_surface(0))); // atan2(v_E, v_N)

        parameters(param_idx.at("current_speed_surface_W")) = W;
        parameters(param_idx.at("current_dir_surface_alpha")) = alpha;

        auto fish_NED_velocity =
         casadi::DM(config_param["fish_velocity_over_ground"]["default"].as<vdob>());
        auto current_fish =
         casadi::DM(config_param["current_fish"]["default"].as<vdob>());

        double V_s = double(norm_2(fish_NED_velocity));
        double chi_s = std::atan2(
            double(fish_NED_velocity(1)),
            double(fish_NED_velocity(0))); // atan2(v_sE, v_sN)

        // fish course in water frame (surface current)
        double psi_s = atan2(V_s*sin(chi_s)-W*sin(alpha),V_s*cos(chi_s)-W*cos(alpha));

        parameters(param_idx.at("fish_speed_NE_V_s")) = V_s;
        parameters(param_idx.at("fish_course_NED_chi_s")) = chi_s;
        parameters(param_idx.at("fish_course_WF_psi_s")) = psi_s;
      }

      void setup_inputs(const YAML::Node& config, dds::sub::Subscriber subscriber)
      {
        BOOST_LOG_TRIVIAL(trace) << "setup inputs";
        // inputs:
        // - GPS origin as Double2
        // - Vessel pos vel: PosInfo
        // - Vessel heading rot: GyroInfo
        // - Vessel water speed? LogInfo

        // - Fish pos vel: PosInfo
        // - Fish relative Double3
        // - Keep as Bit

        typedef std::vector<double> vdob;
        using casadi::Slice;
        auto readerQos = subscriber.default_datareader_qos();
        const auto config_input = config["inputs"];

        auto add_input_Bit =
         [=](const std::string& variable){

           auto topicName = config_input[variable]["topic"].as<std::string>();
           input_readers.emplace(
               variable,
               dds::sub::DataReader<fkin::Bit>(
                   subscriber,
                   dds::topic::Topic<fkin::Bit>(
                       subscriber.participant(),
                       topicName),
                   readerQos));
         };

        auto add_input_Double2 =
         [=](const std::string& variable){

           auto topicName = config_input[variable]["topic"].as<std::string>();
           input_readers.emplace(
               variable,
               dds::sub::DataReader<ratatosk::types::Double2>(
                   subscriber,
                   dds::topic::Topic<ratatosk::types::Double2>(
                       subscriber.participant(),
                       topicName),
                   readerQos));
         };

        auto add_input_Double3 =
         [=](const std::string& variable){

           auto topicName = config_input[variable]["topic"].as<std::string>();
           input_readers.emplace(
               variable,
               dds::sub::DataReader<ratatosk::types::Double3>(
                   subscriber,
                   dds::topic::Topic<ratatosk::types::Double3>(
                       subscriber.participant(),
                       topicName),
                   readerQos));
         };

        auto add_input_PosInfo =
         [=](const std::string& variable){

           auto topicName = config_input[variable]["topic"].as<std::string>();
           input_readers.emplace(
               variable,
               dds::sub::DataReader<ratatosk::types::PosInfo>(
                   subscriber,
                   dds::topic::Topic<ratatosk::types::PosInfo>(
                       subscriber.participant(),
                       topicName),
                   readerQos));
         };

        auto add_input_GyroInfo =
         [=](const std::string& variable){

           auto topicName = config_input[variable]["topic"].as<std::string>();
           input_readers.emplace(
               variable,
               dds::sub::DataReader<ratatosk::types::GyroInfo>(
                   subscriber,
                   dds::topic::Topic<ratatosk::types::GyroInfo>(
                       subscriber.participant(),
                       topicName),
                   readerQos));
         };

        add_input_Double2("GPS_origin");
        add_input_PosInfo("vessel_pos_info");
        add_input_GyroInfo("vessel_gyro_info");
        add_input_PosInfo("fish_pos_info");
        add_input_Double3("fish_relative_pos");
        add_input_Bit("keep_solution");

        gps_origin = casadi::DM(config_input["GPS_origin"]["default"].as<vdob>());

        // If sample with limited duration is to be used:
        //auto input_age = input_config["max_age_ms"].as<std::int64_t>();
        //readerQos << org::opensplice::core::policy::ReaderLifespan(
        //true,
        //    dds::core::Duration::from_millisecs(input_age));

      }

      void setup_outputs(const YAML::Node& config, dds::pub::Publisher publisher)
      {
        BOOST_LOG_TRIVIAL(trace) << "setup outputs";
        auto writerQos = publisher.default_datawriter_qos();
        outputId = config["outputs"]["id"].as<std::string>();

        output_writers.emplace(
            "vessel_speed",
            dds::pub::DataWriter<ratatosk::types::DoubleVal>(
                publisher,
                dds::topic::Topic<ratatosk::types::DoubleVal>(
                    publisher.participant(),
                    config["outputs"]["vessel_speed"]["topic"].as<std::string>()),
                writerQos));

        output_writers.emplace(
            "vessel_course_rate",
            dds::pub::DataWriter<ratatosk::types::DoubleVal>(
                publisher,
                dds::topic::Topic<ratatosk::types::DoubleVal>(
                    publisher.participant(),
                    config["outputs"]["vessel_course_rate"]["topic"].as<std::string>()),
                writerQos));

        output_writers.emplace(
            "deploy_position",
            dds::pub::DataWriter<ratatosk::types::Double2>(
                publisher,
                dds::topic::Topic<ratatosk::types::Double2>(
                    publisher.participant(),
                    config["outputs"]["deploy_position"]["topic"].as<std::string>()),
                writerQos));

        output_writers.emplace(
            "collide_position",
            dds::pub::DataWriter<ratatosk::types::Double2>(
                publisher,
                dds::topic::Topic<ratatosk::types::Double2>(
                    publisher.participant(),
                    config["outputs"]["collide_position"]["topic"].as<std::string>()),
                writerQos));

        output_writers.emplace(
            "deploy_time",
            dds::pub::DataWriter<ratatosk::types::DoubleVal>(
                publisher,
                dds::topic::Topic<ratatosk::types::DoubleVal>(
                    publisher.participant(),
                    config["outputs"]["deploy_time"]["topic"].as<std::string>()),
                writerQos));

        // x(t) for t: [t0, tf] for vessel
        output_writers.emplace(
            "trajectory_vessel",
            dds::pub::DataWriter<fkin::BatchKinematics2D>(
                publisher,
                dds::topic::Topic<fkin::BatchKinematics2D>(
                    publisher.participant(),
                    config["outputs"]["trajectory_vessel"]["topic"].as<std::string>()),
                writerQos));

        // rate of turn for t: [t0, tf] for vessel
        output_writers.emplace(
            "trajectory_vessel_rot",
            dds::pub::DataWriter<fkin::BatchIdVec1d>(
                publisher,
                dds::topic::Topic<fkin::BatchIdVec1d>(
                    publisher.participant(),
                    config["outputs"]["trajectory_vessel_rot"]["topic"].as<std::string>()),
                writerQos));


        // x(t) for t: [t0, tf] for fish
        output_writers.emplace(
            "trajectory_fish",
            dds::pub::DataWriter<fkin::BatchKinematics2D>(
                publisher,
                dds::topic::Topic<fkin::BatchKinematics2D>(
                    publisher.participant(),
                    config["outputs"]["trajectory_fish"]["topic"].as<std::string>()),
                writerQos));

        // nlp configuration info
        output_writers.emplace(
            "nlp_config",
            dds::pub::DataWriter<fkin::NlpConfig>(
                publisher,
                dds::topic::Topic<fkin::NlpConfig>(
                    publisher.participant(),
                    config["outputs"]["nlp_config"]["topic"].as<std::string>()),
                writerQos));

        // nlp statistics
        output_writers.emplace(
            "nlp_stats",
            dds::pub::DataWriter<fkin::OptiStats>(
                publisher,
                dds::topic::Topic<fkin::OptiStats>(
                    publisher.participant(),
                    config["outputs"]["nlp_stats"]["topic"].as<std::string>()),
                writerQos));

      }
    };

    PursePlanner::PursePlanner(
        const YAML::Node& config,
        boost::statechart::fifo_scheduler<> & scheduler,
        boost::statechart::fifo_scheduler<>::processor_handle machine,
        dds::pub::Publisher publisher,
        dds::sub::Subscriber subscriber) :
      m_impl( new PursePlanner::Impl(config, publisher, subscriber)),
      m_scheduler(scheduler),
      m_stateMachine(machine),
      m_time_step(std::chrono::milliseconds(
           config["config"]["settings"]["time_step_ms"].as<std::int32_t>())),
      m_next_step(std::chrono::steady_clock::now()),
      m_config(config),
      m_retries(10),
      m_keep_solution(false)
    {}

    void PursePlanner::solve(const std::atomic<bool>& cancel_token)
    {
      // solve() performs the following :
      // I.    Fetches DDS input and parameters, x(t0), p into local variables
      // II.   Shifts solution from previous step for warm starting optimization problem
      // III.  Sets initial condition for x(t0) in optimization problem
      // IV.   Solves the optimization problem
      // V.    Verify that the solve succeeded (publish stats) TODO: failure action
      // VI.   Updates solution in nlp problem data structure
      // VII.  Evaluates solution trajectory  (prepare for publish)
      // VIII. Publish solution trajectory and desired now


      using casadi::DM, casadi::DMDict;
      namespace sc = std::chrono;
      using namespace fkin;

      try
      {
        BOOST_LOG_TRIVIAL(trace) << "Solving " << name();

        // I. == Read DDS data =====================================================
        BOOST_LOG_TRIVIAL(debug) << "  I. read DDS data";
        // Data are put into data structures that are used further below

        read_inputs();      // Fetch inputs, sample and hold
        read_parameters();  // Fetch user-configurable parameters from DDS

        // Get constant parameters that are needed in output of BatchKinematics2D (fish course, speeds)
        double U_v = double(model()->parameters(model()->param_idx.at("setting_speed_U_v")));
        double fish_sog = double(model()->parameters(model()->param_idx.at("fish_speed_NE_V_s")));
        double W = double(model()->parameters(model()->param_idx.at("current_speed_surface_W")));
        double alpha = double(model()->parameters(model()->param_idx.at("current_dir_surface_alpha")));
        double D_s = double(model()->parameters(model()->param_idx.at("aim_distance_D_s")));
        double R_x = double(model()->parameters(model()->param_idx.at("setting_Rx")));
        double R_y = double(model()->parameters(model()->param_idx.at("setting_Ry")));
        double t_d = 0.; // set properly later

        // II. == Warm start ==================================================
        BOOST_LOG_TRIVIAL(debug) << "  II. warm start";
        auto nlp_input = DMDict{{"w", DM(model()->formulation->nlp_problem.x_init)}};

        auto d_param = model()->formulation->mpc.at("v").at("extractor")(nlp_input).at("v");
        auto sys_1 = model()->formulation->mpc.at("x_1").at("shifter")(nlp_input).at("w_out");
        auto sys_2 = model()->formulation->mpc.at("x_2").at("shifter")(nlp_input).at("w_out");

        model()->formulation->nlp_problem.x_init = vertcat(d_param, sys_1, sys_2).get_elements();

        // III. == Set variable guesses in NLP ==================================
        BOOST_LOG_TRIVIAL(debug) << "  III. set variable guesses in NLP";
        model()->formulation->set_variable("v", d_param);
        model()->formulation->set_variable("x_1_0", casadi::DM::zeros(model()->x1_dim, 1));
        model()->formulation->set_variable("x_2_0", model()->x_k);
        model()->t_k = sc::duration_cast<sc::seconds>(m_next_step - m_t0).count();

        DM t0_opt = DM(model()->t_k); // if formulation has explicit time else, 0 is ok
        //t0_opt(0) = 0; // Need to be tested more.

        // IV. == Solve NLP =====================================================
        BOOST_LOG_TRIVIAL(debug) << "  IV. solve NLP";
        model()->formulation->set_abort(cancel_token);

        BOOST_LOG_TRIVIAL(debug) << "Skip new solution: " << m_keep_solution;

        if (!m_keep_solution){

          const auto& prob = model()->formulation->nlp_problem;
          auto result = model()->formulation->nlp_solver(
              DMDict({{"x0", prob.x_init},
                      {"lbx", prob.x_lb},
                      {"ubx", prob.x_ub},
                      {"lbg", prob.g_lb},
                      {"ubg", prob.g_ub},
                      {"p", vertcat(model()->parameters, t0_opt, t0_opt)},
                      {"lam_x0", prob.lambda_init},
                      {"lam_g0", model()->lamb_g}}));

          // TODO: may need clever multi-evaluation
          // Maybe parallel evaluation with futures?

          // V. == Check status =================================================
          BOOST_LOG_TRIVIAL(debug) << "  V. check status";
          fkin::OptiStats opti_stats(model()->formulation->stats());
          opti_stats.obj() = double(result["f"]);
          opti_stats.id() = model()->outputId;
          opti_stats.p() = vertcat(model()->parameters, model()->t_k).get_elements();
          opti_stats.x0() = model()->x_k.get_elements();

          static_cast<dds::pub::DataWriter<fkin::OptiStats>>(
              model()->output_writers.at("nlp_stats"))->write(opti_stats);

          if (cancel_token)
          {
            BOOST_LOG_TRIVIAL(info) << "Solve was interrupted by user";
            event(new mimir::EvInterrupt());
            return;
          }


          if(!model()->formulation->nlp_solver.stats().at("success").as_bool())
          {
            BOOST_LOG_TRIVIAL(debug) << model()->formulation->nlp_solver.stats();

            if(--m_retries > 0) {
              // reset initial guess of Lagrange multipliers and decision parameters
              auto v_guess = casadi::DM::zeros(model()->v_dim, 1);
              model()->formulation->set_variable("v", v_guess);
              model()->lamb_g = std::vector<double>(model()->formulation->nlp_problem.g_lb.size(), 0);
              step_time();
              event(new mimir::EvReady());
              return;
            }

            // Max number of failed attempts in a row, giving up.
            event(new mimir::EvInterrupt());
            return;
          }
          else{
            m_retries = 10; // Reset retries to magic number
          }

          // VI. == Update nlp problem  =========================================
          BOOST_LOG_TRIVIAL(debug) << "  VI. update NLP problem";
          model()->formulation->nlp_problem.x_init = result["x"].get_elements();
          model()->formulation->nlp_problem.lambda_init = result["lam_x"].get_elements();
          model()->lamb_g = result["lam_g"].get_elements();
          const auto& nlp_output = model()->formulation->nlp_problem.x_init;


          // VII. == Set solution data  =========================================
          // (maybe unpacker + timegrid as optional for collocation)
          BOOST_LOG_TRIVIAL(debug) << "  VII. set solution data";

          auto decide_param =
           model()->formulation->mpc.at("v").at("extractor")(DMDict({{"w", nlp_output}})).at("v");

          auto unpacked_x1 =
           model()->formulation->mpc.at("x_1").at("unpacker")(DMDict({{"w", nlp_output}}));

          auto solT_x1 =
           model()->formulation->mpc.at("x_1").at("timegrid")(DMDict({{"t0", DM(model()->t_k)},
                                                                      {"v", decide_param}}));
          auto unpacked_x2 =
           model()->formulation->mpc.at("x_2").at("unpacker")(DMDict({{"w", nlp_output}}));

          auto solT_x2 =
           model()->formulation->mpc.at("x_2").at("timegrid")(DMDict({{"t0", DM(model()->t_k)},
                                                                      {"v", decide_param}}));

          BOOST_LOG_TRIVIAL(info) << "Decided parameters: " << decide_param;

          t_d = double(model()->formulation->mpc.at("x_2").at("Tf")(DMDict{{"w", nlp_output}}).at("Tf"));

          double tmp_h = std::pow(R_x - R_y, 2) / std::pow(R_x + R_y, 2);
          double Pi = 2*std::acos(0.0);
          // Ramanujan, second approximation:
          double ellipse_circum = Pi*(R_x + R_y)*(1 + 3*tmp_h/(10 + std::pow(4 - 3*tmp_h, 0.5)));

          double s_pd = 0.5*ellipse_circum; // We use half ellipse circumference approximation..

          // Time interval of interest to provide trajectories
          double t_f = ceil(t_d + (D_s + s_pd)/U_v); // t0 = 0, tf = t_pd = t_d + (D_s + s_pd)/Uv
          // We can use the trajectory function, since x_2 has variable step and no control input
          double deltaT = t_f/model()->formulation->N2;

          if (unpacked_x2["U"].nnz() > 0)
            throw std::runtime_error("PursePlanner assumes that x_2 subsystem has no input u");

          auto decide_param_traj = decide_param;
          decide_param_traj(model()->decide_idx.at("x_2_DT")) = deltaT;

          // Unused:
          /*auto trajx1 = model()->formulation->mpc.at("x_1").at("trajectory")(
            DMDict({{"U", unpacked_x1["U"]},
            {"X0", unpacked_x1["X0"]},
            {"p", model()->parameters},
            {"v", decide_param},
            {"t0", t0_opt}}));*/

          if(!(model()->formulation->nlp_config().technique() == "single-shooting"))
          {
            model()->x_k = unpacked_x2["Xc"](casadi::Slice(), 0); // TO BE REMOVED
            model()->T = solT_x2["T"];
            model()->X = horzcat(unpacked_x2["X0"], unpacked_x2["Xc"]);
            throw std::runtime_error("PursePlanner assumes the technique is single-shooting");
          }
          else
          {
            auto trajectory =
             model()->formulation->mpc.at("x_2").at("trajectory")(
                 DMDict({{"U", unpacked_x2["U"]},
                         {"X0", unpacked_x2["X0"]},
                         {"p", model()->parameters},
                         {"v", decide_param_traj},
                         {"t0", t0_opt}}));

            model()->x_k = trajectory["X"](casadi::Slice(), 0); // TO BE REMOVED
            model()->T = trajectory["T"];
            model()->X = trajectory["X"];

            BOOST_LOG_TRIVIAL(trace) << "T_traj: " << trajectory["T"];
            BOOST_LOG_TRIVIAL(trace) << "X_traj: " << trajectory["X"];
            BOOST_LOG_TRIVIAL(trace) << "Trajectory elements: " << trajectory["T"].nnz();

            // Update output
            model()->T1 = trajectory["T"];
            model()->X1 = trajectory["X"];

          }

          // Control input solution
          model()->U = unpacked_x2["U"];
          model()->Tu = solT_x2["Tu"];

          BOOST_LOG_TRIVIAL(trace) << "T: " << model()->T;
          BOOST_LOG_TRIVIAL(trace) << "Xc: " << unpacked_x2["Xc"];
          BOOST_LOG_TRIVIAL(trace) << "U: " << model()->U;
          BOOST_LOG_TRIVIAL(trace) << "Tu: " << model()->Tu;

        }


        // VIII. == Publish solution trajectory and desired now ================
        BOOST_LOG_TRIVIAL(debug) << "  VIII. publish solution trajectory and desired now";

        if(model()->T.nnz() > 0)
        {
          auto Tx = model()->T.get_elements();
          auto Tu = model()->Tu.get_elements();

          //auto vessel_rot_idx_u = casadi::Slice(0,1);
          fkin::BatchKinematics2D vessel_traj, fish_traj;
          fkin::BatchIdVec1d vessel_rot;//, vessel_acc;
          std::string id = model()->outputId;
          vessel_traj.id() = id;
          fish_traj.id() = id;
          vessel_rot.id() = id;
          vessel_traj.timestamps().reserve(Tx.size());
          vessel_rot.timestamps().reserve(Tx.size());
          fish_traj.timestamps().reserve(Tx.size());

          // Find the time point in solution closest to now
          // Use {Tx, Tu} resp. if a {state, input} variable is to be found
          auto t_expected = m_next_step - m_t0; // Expected time now in optimization
          auto t_now_sim_epoch = std::chrono::steady_clock::now() - m_t0; // Actual time now

          // First time in solution vector that are not less than t_now
          auto time_it = std::lower_bound(
              Tx.begin(), Tx.end(),
              sc::duration_cast<sc::seconds>(t_now_sim_epoch).count());

          // If it is found, publish vessel commands
          if(time_it != Tx.end()){

            // Index of the solution
            auto idx = std::distance(Tx.begin(), time_it);

            // Forward diff (f(idx+1) - f(idx))/deltaT
            auto idx_1 = idx+1;
            auto idx_2 = idx;
            if (static_cast<size_t>(idx + 1) == Tx.size()) {
              // Backward diff (f(idx) - f(idx-1))/deltaT
              idx_1 = idx;
              idx_2 = idx-1;
            }

            double rad_2_deg = 180./(2*std::acos(0.0)); // (180/pi)
            auto course_1 = double(model()->X(model()->state2_idx.at("vessel_course"), idx_1));
            auto course_2 = double(model()->X(model()->state2_idx.at("vessel_course"), idx_2));
            auto time_1 = Tx[idx_1];
            auto time_2 = Tx[idx_2];
            double vessel_course_rate = (course_1 - course_2)/(time_1 - time_2);
            double capped_turn_rate = std::fmin(
                std::fmax(double(vessel_course_rate), -model()->formulation->omega_max),
                model()->formulation->omega_max);
            capped_turn_rate *= rad_2_deg;

            double chi_v = double(model()->X(model()->state2_idx.at("vessel_course"), idx));
            double beta = asin((W/U_v)*sin(alpha - chi_v)); // Side slip.
            // speed over ground
            double vessel_sog = sqrt(U_v*U_v + W*W + 2*U_v*W*cos(chi_v - beta - alpha));

            BOOST_LOG_TRIVIAL(trace) << "Simulation time for commanded output: " << (*time_it);

            static_cast<dds::pub::DataWriter<ratatosk::types::DoubleVal>>(
                model()->output_writers.at("vessel_speed"))->write(
                    ratatosk::types::DoubleVal(vessel_sog));

            static_cast<dds::pub::DataWriter<ratatosk::types::DoubleVal>>(
                model()->output_writers.at("vessel_course_rate"))->write(
                    ratatosk::types::DoubleVal(capped_turn_rate));

            BOOST_LOG_TRIVIAL(debug) << "Commands are: " << vessel_sog << ", " << capped_turn_rate;

          }

          if(m_keep_solution)
          {
            step_time();
            event(new mimir::EvReady());
            return;
          }

          auto t_overtime = t_now_sim_epoch - t_expected; // Is the optimization late?
          // Find deployment position
          auto t_deploy = t_now_sim_epoch - t_overtime + sc::milliseconds(static_cast<int64_t>(t_d*1000));

          double t_late =  sc::duration_cast<sc::milliseconds>(t_overtime).count();
          BOOST_LOG_TRIVIAL(warning) << " Overtime: " << t_late;

          static_cast<dds::pub::DataWriter<ratatosk::types::DoubleVal>>(
              model()->output_writers.at("deploy_time"))->write(
                  ratatosk::types::DoubleVal(double(t_d - t_late/1000)));

          BOOST_LOG_TRIVIAL(debug) << "Deploy in: " << t_d - t_late/1000 << " sec";

          time_it = std::lower_bound(
              Tx.begin(), Tx.end(),
              sc::duration_cast<sc::seconds>(t_deploy).count());

          if(time_it != Tx.end()) {

            auto idx = std::distance(Tx.begin(), time_it);
            auto deploy_pos = model()->X(model()->state2_idx.at("vessel_NE"), idx);

            BOOST_LOG_TRIVIAL(debug) << "Deploy pos is: " << deploy_pos;

            static_cast<dds::pub::DataWriter<ratatosk::types::Double2>>(
                model()->output_writers.at("deploy_position"))->write(
                    ratatosk::types::Double2(double(deploy_pos(0)), double(deploy_pos(1))));

          }

          // Find collision point
          auto t_collide = t_deploy + sc::seconds(static_cast<int64_t>(D_s/U_v));


          time_it = std::lower_bound(
              Tx.begin(), Tx.end(),
              sc::duration_cast<sc::seconds>(t_collide).count());

          if(time_it != Tx.end()) {

            auto idx = std::distance(Tx.begin(), time_it);
            auto collide_pos = model()->X(model()->state2_idx.at("vessel_NE"), idx);

            static_cast<dds::pub::DataWriter<ratatosk::types::Double2>>(
                model()->output_writers.at("collide_position"))->write(
                    ratatosk::types::Double2(double(collide_pos(0)), double(collide_pos(1))));

          }

          auto end_idx = model()->X.size2();

          // Convert simulation time to wall clock time
          for (casadi_int i = 0; i < end_idx; ++ i) //auto& t : Tx)
          {
            auto steady_time = m_t0 + sc::milliseconds(static_cast<int64_t>(double(Tx[i])*1000));
            auto wall_time = (steady_time - m_t0 + m_t0_wall); // drifts if system_clock changes
            auto epoch_time = sc::duration_cast<sc::milliseconds>(wall_time.time_since_epoch()).count();
            vessel_traj.timestamps().push_back(fkin::Timestamp(epoch_time));
            fish_traj.timestamps().push_back(fkin::Timestamp(epoch_time));
            vessel_rot.timestamps().push_back(fkin::Timestamp(epoch_time));
          }

          // Calculate rate of turn with finite difference
          casadi_int n_elem = model()->X.size2();
          std::vector<casadi_int> last_row{n_elem - 1};
          std::vector<casadi_int> secnd_last_col{n_elem - 2};
          // Shifts one and keeps last element (band(n,1) is sub diagonal, so we transpose)
          casadi::DM shift_operator = casadi::DM(casadi::Sparsity::band(n_elem,1).T()
           + casadi::Sparsity::triplet(n_elem, n_elem, last_row, secnd_last_col));

          casadi::DM vessel_courses = model()->X(model()->state2_idx.at("vessel_course"),
           casadi::Slice()).T();
          casadi::DM courses_shifted = mtimes(shift_operator, vessel_courses);
          casadi::DM delta_courses = courses_shifted - vessel_courses;
          casadi::DM delta_times = mtimes(shift_operator, Tx) - Tx;
          casadi::DM vessel_turn_rate = delta_courses/delta_times;

          BOOST_LOG_TRIVIAL(trace) << "calc vessel_turn_rate are: " << vessel_turn_rate;

          // Loop over columns of X to fetch solution trajectory
          for(casadi_int i = 0; i < end_idx /* model()->X.size2() */; ++i)
          {
            casadi::DM stateX = model()->X(casadi::Slice(), i);

            casadi::DM vessel_NE(stateX(model()->state2_idx.at("vessel_NE")));
            casadi::DM fish_NE(stateX(model()->state2_idx.at("fish_NE")));

            casadi::DM vessel_course(stateX(model()->state2_idx.at("vessel_course")));

            double capped_turn_rate = std::fmin(
                std::fmax(double(vessel_turn_rate(i)), -model()->formulation->omega_max),
                model()->formulation->omega_max);

            double chi_v = double(vessel_course);
            double beta = asin((W/U_v)*sin(alpha - chi_v)); // Side slip.
            // speed over ground
            double vessel_sog = sqrt(U_v*U_v + W*W + 2*U_v*W*cos(chi_v - beta - alpha));

            vessel_traj.batch().push_back( Kinematics2D(
                 id,
                 Vector2d(double(vessel_NE(0)), double(vessel_NE(1))),
                 Vector1d(double(stateX(model()->state2_idx.at("vessel_course")))),
                 Vector1d(vessel_sog)));

            fish_traj.batch().push_back( Kinematics2D(
                 id,
                 Vector2d(double(fish_NE(0)), double(fish_NE(1))),
                 Vector1d(0.),
                 Vector1d(double(fish_sog))));

            vessel_rot.batch().push_back( IdVec1d(id, Vector1d(double(capped_turn_rate))) );

          }

          static_cast<dds::pub::DataWriter<fkin::BatchKinematics2D>>(
              model()->output_writers.at("trajectory_vessel"))->write(vessel_traj);

          static_cast<dds::pub::DataWriter<fkin::BatchKinematics2D>>(
              model()->output_writers.at("trajectory_fish"))->write(fish_traj);

          static_cast<dds::pub::DataWriter<fkin::BatchIdVec1d>>(
              model()->output_writers.at("trajectory_vessel_rot"))->write(vessel_rot);

        }

        plot(model()->plot_on);
        step_time();
        event(new mimir::EvReady());
      }
      catch (casadi::CasadiException &e)
      {
        BOOST_LOG_TRIVIAL(fatal) << "Casadi threw exception in " << name() << ": " << e.what();
        event(new mimir::EvError());
        throw;
      }
      catch (std::out_of_range &e)
      {
        BOOST_LOG_TRIVIAL(fatal) << "Out of bounds: " << name() << ": " << e.what();
        event(new mimir::EvError());
        throw;
      }
      catch (std::exception &e)
      {
        BOOST_LOG_TRIVIAL(fatal) << name() << ": " << e.what();
        event(new mimir::EvError());
        throw;
      }
      catch (...)
      {
        BOOST_LOG_TRIVIAL(fatal) << name() << " exception thrown";
        event(new mimir::EvError());
        throw;
      }
    }

    void PursePlanner::read_inputs()
    {
      using casadi::DM;
      using casadi::DMDict;
      typedef std::vector<double> vdob;

      auto readBit =
       [=](const std::string& key) -> std::pair<bool, bool>
       {
        auto samples = static_cast<dds::sub::DataReader<fkin::Bit>>(
            model()->input_readers.at(key))->read();
        if (samples.length() > 0){
          auto sample = (--samples.end());
          if(sample->info().valid()){
            return std::make_pair(
                true,
                sample->data().value());
          }
        }
        return std::make_pair(false, false);
       };

      auto readDouble2 =
       [=](const std::string& key) -> std::pair<bool, casadi::DM>
       {
         auto samples = static_cast<dds::sub::DataReader<ratatosk::types::Double2>>(
             model()->input_readers.at(key))->read();
         if (samples.length() > 0){
           auto sample = (--samples.end());
           if(sample->info().valid()){
             return std::make_pair(
                 true,
                 casadi::DM(std::vector<double>{sample->data().x(), sample->data().y()}));
           }
         }
         return std::make_pair(false, casadi::DM());
       };

      /*
      auto readDouble3 =
       [=](const std::string& key) -> std::pair<bool, casadi::DM>
       {
        auto samples = static_cast<dds::sub::DataReader<ratatosk::types::Double3>>(
            model()->input_readers.at(key))->read();
        if (samples.length() > 0){
          auto sample = (--samples.end());
          if(sample->info().valid()){
            return std::make_pair(
                true,
                casadi::DM(std::vector<double>{
                   sample->data().x(), sample->data().y(), sample->data().z()}));
          }
        }
        return std::make_pair(false, casadi::DM());
       }; */

      auto readPosInfo =
       [=](const std::string& key) -> std::pair<bool, ratatosk::types::PosInfo>
       {
        auto samples = static_cast<dds::sub::DataReader<ratatosk::types::PosInfo>>(
            model()->input_readers.at(key))->read();
        if (samples.length() > 0){
          auto sample = (--samples.end());
          if(sample->info().valid()){

            return std::make_pair(true, sample->data());
          }
        }
        return std::make_pair(false, ratatosk::types::PosInfo());
       };

      auto readGyroInfo =
       [=](const std::string& key) -> std::pair<bool, ratatosk::types::GyroInfo>
       {
        auto samples = static_cast<dds::sub::DataReader<ratatosk::types::GyroInfo>>(
            model()->input_readers.at(key))->read();
        if (samples.length() > 0){
          auto sample = (--samples.end());
          if(sample->info().valid()){

            return std::make_pair(true, sample->data());
          }
        }
        return std::make_pair(false, ratatosk::types::GyroInfo());
       };

      auto resD2 = readDouble2("GPS_origin");
      if(resD2.first)
        model()->gps_origin = resD2.second;

      auto vesselPosInfo = readPosInfo("vessel_pos_info");
      auto fishPosInfo = readPosInfo("fish_pos_info");
      auto vesselGyroInfo = readGyroInfo("vessel_gyro_info");

      // States of path planner are defined in ::formulate_dynamics:

      casadi::DM xk = model()->x_k;
      BOOST_LOG_TRIVIAL(debug) << "STATE before: " << xk;

      double PI180 = 2.*std::acos(0.0)/180.; //std::atan(1)*4./180.;

      if(vesselPosInfo.first)
      {
        auto res = model()->toNED(
            DMDict{{"GPS", DM(vdob{vesselPosInfo.second.lat(), vesselPosInfo.second.lon()})},
                   {"GPS0", model()->gps_origin}}).at("NED");

        xk(model()->state2_idx.at("vessel_NE")) = res(casadi::Slice(0,2));
      }

      if(vesselGyroInfo.first)
      {
        // This may be wrong. Is true heading the same as course?
        xk(model()->state2_idx.at("vessel_course")) = vesselGyroInfo.second.hdt()*PI180;
      }

      if(fishPosInfo.first)
      {
        auto resPI = model()->toNED(
            DMDict{{"GPS", DM(vdob{fishPosInfo.second.lat(), fishPosInfo.second.lon()})},
                   {"GPS0", model()->gps_origin}}).at("NED");

        BOOST_LOG_TRIVIAL(trace) << "NED fish: " << resPI;
        xk(model()->state2_idx.at("fish_NE")) = resPI(casadi::Slice(0,2));
      }

      // Set water frame to origin
      xk(model()->state2_idx.at("water_NE")) = casadi::DM::zeros(2,1);


      BOOST_LOG_TRIVIAL(debug) << "STATE after: " << xk;

      model()->x_k = xk;

      // Parameters, which in reality are initial conditions..
      model()->parameters(model()->param_idx.at("vessel_NE_q0")) =
       xk(model()->state2_idx.at("vessel_NE"));

      model()->parameters(model()->param_idx.at("fish_NE_p_s0")) =
       xk(model()->state2_idx.at("fish_NE"));

      // Whether the user requests to keep the suggested trajectory or not.
      auto keep_solution = readBit("keep_solution");
      if (keep_solution.first)
        m_keep_solution = keep_solution.second;

    }

    void PursePlanner::read_parameters()
    {

      auto readIdVec1 =
       [=](const std::string& key){
         auto samples = static_cast<dds::sub::DataReader<fkin::IdVec1d>>(
             model()->parameter_inputs.at(key))->read();
         if (samples.length() > 0){
           auto sample = (--samples.end());
           if(sample->info().valid()){
             model()->parameters(model()->param_idx.at(key)) = sample->data().vec().x();
           }
         }
       };

      auto readIdVec2 =
       [=](const std::string& key){
         auto samples = static_cast<dds::sub::DataReader<fkin::IdVec2d>>(
             model()->parameter_inputs.at(key))->read();
         if (samples.length() > 0){
           auto sample = (--samples.end());
           if(sample->info().valid()){
             model()->parameters(model()->param_idx.at(key)) =
              casadi::DM(std::vector<double>{sample->data().vec().x(), sample->data().vec().y()});
           }
         }
       };

      /*
      auto readDouble1 =
       [=](const std::string& key){
         auto samples = static_cast<dds::sub::DataReader<ratatosk::types::DoubleVal>>(
             model()->parameter_inputs.at(key))->read();
         if (samples.length() > 0){
           auto sample = (--samples.end());
           if(sample->info().valid()){
             model()->parameters(model()->param_idx.at(key)) = sample->data().val();
           }
         }
       };

      auto readDouble3 =
       [=](const std::string& key){
         auto samples = static_cast<dds::sub::DataReader<ratatosk::types::Double3>>(
             model()->parameter_inputs.at(key))->read();
         if (samples.length() > 0){
           auto sample = (--samples.end());
           if(sample->info().valid()){
             model()->parameters(model()->param_idx.at(key)) =
              casadi::DM(std::vector<double>{
                   sample->data().x(), sample->data().y(), sample->data().z()});
           }
         }
       };
      */

      auto readDouble2 =
       [=](const std::string& key) -> std::pair<bool, casadi::DM>
       {
         auto samples = static_cast<dds::sub::DataReader<ratatosk::types::Double2>>(
             model()->parameter_inputs.at(key))->read();
         if (samples.length() > 0){
           auto sample = (--samples.end());
           if(sample->info().valid()){
             //model()->parameters(model()->param_idx.at(key)) =
              return std::make_pair(
                  true,
                  casadi::DM(std::vector<double>{sample->data().x(), sample->data().y()}));
           }
         }
         return std::make_pair(false, casadi::DM());
       };

      readIdVec1("setting_speed_U_v");
      readIdVec2("setting_radii");
      readIdVec1("aim_distance_D_s");
      readIdVec1("fish_margin_d_f");
      readIdVec1("fish_depth_z_s");
      readIdVec2("leadline_tau_ll_z_d");
      readIdVec1("sink_margin_z_min");

      auto surf_current = readDouble2("current_surface");
      auto fish_vog = readDouble2("fish_velocity_over_ground");
      auto fish_current = readDouble2("current_fish");

      if(surf_current.first)
      {
        double W = double(norm_2(surf_current.second));
        double alpha = std::atan2(
            double(surf_current.second(1)),
            double(surf_current.second(0))); // atan2(v_E, v_N)

        model()->parameters(model()->param_idx.at("current_speed_surface_W")) = W;
        model()->parameters(model()->param_idx.at("current_dir_surface_alpha")) = alpha;

      }

      if(fish_vog.first && fish_current.first)
      {
        auto current_fish = fish_current.second;
        auto fish_NED_velocity = fish_vog.second;

        double V_s = double(norm_2(fish_NED_velocity));
        double chi_s = std::atan2(
            double(fish_NED_velocity(1)),
            double(fish_NED_velocity(0))); // atan2(v_sE, v_sN)

        double W = double(model()->parameters(model()->param_idx.at("current_speed_surface_W")));
        double alpha = double(model()->parameters(model()->param_idx.at("current_dir_surface_alpha")));
        // fish course in water frame (surface current)
        double psi_s = atan2(V_s*sin(chi_s)-W*sin(alpha),V_s*cos(chi_s)-W*cos(alpha));

        model()->parameters(model()->param_idx.at("fish_speed_NE_V_s")) = V_s;
        model()->parameters(model()->param_idx.at("fish_course_NED_chi_s")) = chi_s;
        model()->parameters(model()->param_idx.at("fish_course_WF_psi_s")) = psi_s;

      }

      // Note: some parameters are set in set_input, because they depend on xk (initial conditions).

      BOOST_LOG_TRIVIAL(debug) << "Parameter vector is: " << model()->parameters;

    }

    void PursePlanner::initialize(const std::atomic<bool>&)
    {
      try
      {
        using namespace std::chrono_literals;
        BOOST_LOG_TRIVIAL(debug) << name() << " is initializing";
        std::this_thread::sleep_for(2s); // use waitset instead for inputs

        m_retries = 10; // Magic number of retries before giving up;
        model()->x_k = model()->x_0;
        model()->t_k = 0.;

        // TODO: Need check to confirm that required signals are available
        read_inputs(); // sets state vector x_k for planner formulation
        read_parameters(); // sets parameter vector for planner formulation

        m_t0 = std::chrono::steady_clock::now(); // https://stackoverflow.com/a/35282833
        m_t0_wall = std::chrono::system_clock::now();

        // Cold start parameter guesses
        auto v_guess = casadi::DM::zeros(model()->v_dim, 1);
        model()->formulation->set_variable("v", v_guess);
        model()->lamb_g = std::vector<double>(model()->formulation->nlp_problem.g_lb.size(), 0);
        model()->formulation->set_variable("x_1_0", casadi::DM::zeros(model()->x1_dim, 1));
        model()->formulation->set_variable("x_2_0", model()->x_k);

        // DDS Publish configured nlp settings
        static_cast<dds::pub::DataWriter<fkin::NlpConfig>>(
            model()->output_writers.at("nlp_config"))->write(model()->formulation->nlp_config());

        m_next_step = std::chrono::steady_clock::now();
        event(new mimir::EvReady());
      }
      catch (...)
      {
        event(new mimir::EvError());
        throw;
      }
    }

    void PursePlanner::timer(const std::atomic<bool>& cancel_token)
    {

      // There is no sanity checks in terms of real-timeliness of execution.
      // probably want external time point to wait until.
      typedef std::chrono::milliseconds scm;
      auto now = std::chrono::steady_clock::now();
      auto time_p = m_next_step;
      auto timeleft = std::chrono::duration_cast<scm>(time_p - now);
      BOOST_LOG_TRIVIAL(debug) << "Time left until next step: " << timeleft.count() << " ms";

      auto fraction = (m_time_step/100 > scm(0) ? m_time_step/100 :
       (m_time_step/50 > scm(0) ? m_time_step/50 :
        (m_time_step/25 > scm(0) ? m_time_step/25 :
         (m_time_step/10 > scm(0) ? m_time_step/10 :
          (m_time_step/2 > scm(0) ?  m_time_step/2  : scm(1))))));

      while(now < time_p)
      {
        auto diff = time_p - now;
        std::this_thread::sleep_for(diff < fraction ? diff : fraction);
        if(cancel_token)
          break;
        now = std::chrono::steady_clock::now();
      }

      if(cancel_token)
      {
        BOOST_LOG_TRIVIAL(trace) << "Timer Canceled";
        event(new mimir::EvInterrupt());
      }
      else
      {
        BOOST_LOG_TRIVIAL(trace) << "Timeout";
        event(new mimir::EvTimeout());
      }
    }

    void PursePlanner::event(boost::statechart::event_base * const event)
    {
      // who is responsible for the pointer?
      m_scheduler.queue_event(
          m_stateMachine,
          make_intrusive(event));
    }

    PursePlanner::~PursePlanner() = default;


    void PursePlanner::plot(bool do_plot)
    {
      if(!do_plot)
        return;

#ifdef MIMIR_WITH_GNUPLOT
      using casadi::Slice;
      model()->gnuplot << "plot '-' with lines ls 1 title 'chi'\n"//, "
                       //<< " '-' with steps ls 2 title 'psi_{dot}'\n";//, "
                       //<< " '-' with lines ls 3 title 'slack'\n";
      model()->gnuplot.send1d(std::make_tuple(
           model()->T.get_elements(),
           model()->X(model()->state2_idx.at("vessel_course"), Slice()).get_elements()));
      /*model()->gnuplot.send1d(std::make_tuple(
           model()->Tu.get_elements(),
           model()->U(model()->input2_idx.at("vessel_acc"), Slice()).get_elements()));*/
      /*model()->gnuplot.send1d(std::make_tuple(
           model()->Tu.get_elements(),
           model()->U(model()->input2_idx.at("slack_dist"), Slice()).get_elements()));*/


      model()->gnuplot.flush();

      model()->ne_plot << "plot '-' with lines ls 1 title 'vessel', "
                       << " '-' with lines ls 2 title 'fish'\n";
      model()->ne_plot.send1d(std::make_tuple(
           model()->X(model()->state2_idx.at("vessel_NE").start + 1, Slice()).get_elements(),
           model()->X(model()->state2_idx.at("vessel_NE").start, Slice()).get_elements()));
      model()->ne_plot.send1d(std::make_tuple(
           model()->X(model()->state2_idx.at("fish_NE").start + 1, Slice()).get_elements(),
           model()->X(model()->state2_idx.at("fish_NE").start, Slice()).get_elements()));

      model()->ne_plot.flush();

#else
      BOOST_LOG_TRIVIAL(trace) << "Gnuplot not available";
#endif
    }
  }
}
