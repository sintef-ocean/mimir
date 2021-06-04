#include <chrono>
#include <cmath>
#include <exception>
#include <thread>
#include <vector>

#include <boost/log/trivial.hpp>

#include "mimir/algorithm/FishSchool.hpp"
#include "mimir/StateMachineFwd.hpp"
#include "mimir/FkinDds.hpp"

#include <casadi/casadi.hpp>

namespace mimir
{
  namespace algorithm
  {
    class FishSchool::Impl
    {
    public:
      Impl(
          const YAML::Node& config,
          dds::pub::Publisher publisher,
          dds::sub::Subscriber subscriber) :
        inputReader(dds::core::null),
        outputWriter(dds::core::null),
        initialPosition(dds::core::null),
        initialEuler(dds::core::null)
      {

        // Setup input reader
        const auto input_config = config["inputs"]["fish_ctrl"];
        const auto input_age = input_config["max_age_ms"].as<std::int64_t>();

        auto readerQos = subscriber.default_datareader_qos();
        //readerQos << dds::core::policy::Durability::TransientLocal()
        //          << dds::core::policy::Reliability::Reliable();

        if(input_age > 0)
          readerQos << org::opensplice::core::policy::ReaderLifespan(
              true,
              dds::core::Duration::from_millisecs(input_age));

        const auto topicName = input_config["topic"].as<std::string>();
        const auto id = input_config["id"].as<std::string>();
        auto topic = dds::topic::Topic<fkin::IdVec3d>(
            subscriber.participant(),
            topicName);
        auto filter = dds::topic::Filter("id = %0", {id});

        inputReader = dds::sub::DataReader<fkin::IdVec3d>(
            subscriber,
            dds::topic::ContentFilteredTopic<fkin::IdVec3d>(
                topic,
                topicName + id,
                filter),
            readerQos);

        // Setup output writer
        auto writerQos = publisher.default_datawriter_qos();
        //writerQos << dds::core::policy::Durability::TransientLocal();


        outputWriter = dds::pub::DataWriter<fkin::Kinematics6D>(
            publisher,
            dds::topic::Topic<fkin::Kinematics6D>(
                publisher.participant(),
                config["outputs"]["kinematics"]["topic"].as<std::string>()),
            writerQos);

        outputId = config["outputs"]["kinematics"]["id"].as<std::string>();

        // Setup initial condition reader
        readerQos << org::opensplice::core::policy::ReaderLifespan(
            false,
            dds::core::Duration::infinite());

        const auto initCondConfig = config["initial_conditions"]["position"];
        const auto initTopicName = initCondConfig["topic"].as<std::string>();
        const auto idInit = initCondConfig["id"].as<std::string>();
        waitPosMs = initCondConfig["max_wait_ms"].as<std::int32_t>();
        topic = dds::topic::Topic<fkin::IdVec3d>(
            subscriber.participant(),
            initTopicName);
        filter = dds::topic::Filter("id = %0", {idInit});

        initialPosition = dds::sub::DataReader<fkin::IdVec3d>(
            subscriber,
            dds::topic::ContentFilteredTopic<fkin::IdVec3d>(
                topic,
                initTopicName + idInit,
                filter),
            readerQos);

        fallbackInitPosition =
             config["initial_conditions"]["position"]["fallback"]
             .as<std::vector<double>>();

        const auto initCondConfig2 = config["initial_conditions"]["euler"];
        const auto initTopicName2 = initCondConfig2["topic"].as<std::string>();
        const auto idInit2 = initCondConfig2["id"].as<std::string>();
        waitEulerMs = initCondConfig["max_wait_ms"].as<std::int32_t>();
        topic = dds::topic::Topic<fkin::IdVec3d>(
            subscriber.participant(),
            initTopicName2);
        filter = dds::topic::Filter("id = %0", {idInit2});

        initialEuler = dds::sub::DataReader<fkin::IdVec3d>(
            subscriber,
            dds::topic::ContentFilteredTopic<fkin::IdVec3d>(
                topic,
                initTopicName2 + idInit2,
                filter),
            readerQos);

        fallbackInitEuler =
         config["initial_conditions"]["euler"]["fallback"]
         .as<std::vector<double>>();

        using namespace casadi;

        auto p = SX::sym("p", 3); // N E D
        auto eu = SX::sym("Euler", 2); // pitch and yaw
        auto x_sym = SX::vertcat({p,eu});
        auto RT = SX::zeros(5,5);
        auto vb = SX::sym("vb", 3); // Body linear velocity
        auto wb = SX::sym("wb", 2); // Body angular velocity (pitch yaw)
        auto zd = SX::sym("z_d", 1); // Desired depth

        auto rotYTheta = SX::vertcat(
            {
             SX::horzcat({cos(eu(0)),0,sin(eu(0))}),
             SX::horzcat({0,1,0}),
             SX::horzcat({-sin(eu(0)),0,cos(eu(0))})
            });

        auto rotZPsi = SX::vertcat(
            {
             SX::horzcat({cos(eu(1)),-sin(eu(1)), 0}),
             SX::horzcat({sin(eu(1)),cos(eu(1)), 0}),
             SX::horzcat({0,0,1})
            });
        auto rot = mtimes(rotZPsi, rotYTheta);

        RT(Slice(0,3),Slice(0,3)) = rot;
        RT(Slice(3,5),Slice(3,5)) =
         SX::vertcat({
                      SX::horzcat({1,0}),
                      SX::horzcat({0,1/cos(eu(0))})});

        auto pitch_rate_ctrl = 0.1*(p(2) - zd)*cos(2*eu(0)) - eu(0);
        auto u_sym = SX::vertcat({vb, wb(1), zd}); // 3d velocity, rate of turn, desired depth
        auto u_vec = SX::vertcat({vb,pitch_rate_ctrl, wb(1)});
        auto rhs = mtimes(RT,u_vec);

        SXDict ode = {{"x", x_sym}, {"p", u_sym}, {"ode", rhs}};

        model = casadi::integrator(
            "fish_school", "cvodes", ode,
            {{"tf", config["time_step_ms"].as<double>()/1000.}});
        u = casadi::DM({1.,0.,0.,0.,0.});

      }
      ~Impl() {}
      dds::sub::DataReader<fkin::IdVec3d> inputReader;
      dds::pub::DataWriter<fkin::Kinematics6D> outputWriter;
      dds::sub::DataReader<fkin::IdVec3d> initialPosition;
      dds::sub::DataReader<fkin::IdVec3d> initialEuler;
      std::chrono::steady_clock::time_point t0;
      int32_t waitPosMs, waitEulerMs;
      std::vector<double> fallbackInitPosition, fallbackInitEuler;
      std::string outputId;
      casadi::Function model;
      casadi::DM u;
      casadi::DM x_k;
    };


    FishSchool::FishSchool(
        const YAML::Node& config,
        boost::statechart::fifo_scheduler<>& scheduler,
        boost::statechart::fifo_scheduler<>::processor_handle machine,
        dds::pub::Publisher publisher,
        dds::sub::Subscriber subscriber) :
      m_impl( new FishSchool::Impl(config, publisher, subscriber)),
      m_scheduler(scheduler),
      m_stateMachine(machine),
      m_time_step(std::chrono::milliseconds(config["time_step_ms"].as<std::int32_t>())),
      m_next_step(std::chrono::steady_clock::now()),
      m_config(config)
      // Allow for simulation speed factor
    {}

    FishSchool::~FishSchool() = default;
    void FishSchool::solve(const std::atomic<bool>& cancel_token)
    {
      try
      {
        BOOST_LOG_TRIVIAL(trace) << "FishSchool solve..";

        // Publish x_k(t_now) (we do not publish future samples)
        const auto x_k = std::vector<double>(m_impl->x_k);
        const auto u_k = std::vector<double>(m_impl->u);

        m_impl->outputWriter << fkin::Kinematics6D(
            m_impl->outputId,
            fkin::Vector3d(x_k[0], x_k[1], x_k[2]),
            fkin::Vector3d(u_k[0], u_k[1], u_k[2]),
            fkin::Vector3d(0., x_k[3], x_k[4]));

        // Fetch input, sample and hold
        auto samples = m_impl->inputReader.read();

        // Set new desired value
        if(samples.length() > 0)
        {
          auto sample = (--samples.end());
          if(sample->info().valid())
          {
            auto sD = sample->data();
            m_impl->u = casadi::DM({sD.vec().x(), 0., 0., sD.vec().z(), sD.vec().y()});
            BOOST_LOG_TRIVIAL(trace) << "Got a sample: " << m_impl->u;
          }
        }
        else
        {
          BOOST_LOG_TRIVIAL(trace) << "No valid input sample for " << name() << ", using old";
        }

        if(cancel_token)
        {
          BOOST_LOG_TRIVIAL(info) << name() << " interrupted";
          event(new mimir::EvInterrupt());
          return;
        }

        // Integrate
        auto result = m_impl->model(
            casadi::DMDict( {{"x0", m_impl->x_k}, {"p", m_impl->u}} ));

        if(cancel_token)
        {
          BOOST_LOG_TRIVIAL(info) << name() << " interrupted";
          event(new mimir::EvInterrupt());
          return;
        }

        // Update output (not published yet)
        m_impl->x_k = result["xf"];
        step_time();
        event(new mimir::EvReady());

      }
      catch(casadi::CasadiException &e)
      {
        BOOST_LOG_TRIVIAL(fatal) << name() << " Casadi Exception: " << e.what();
        event(new mimir::EvError());
        throw;
      }
      catch (...)
      {
        // Need to post event in case of exception so that the state machine
        // knows that the job has finished.
        BOOST_LOG_TRIVIAL(fatal) << name() << " exception thrown";
        event(new mimir::EvError());
        throw;
      }
    }
    void FishSchool::initialize(const std::atomic<bool>&)
    {
      try
      {
        BOOST_LOG_TRIVIAL(debug) << name() << " is initializing";

        auto waitPos = dds::core::cond::WaitSet();
        auto readCondPos = dds::sub::cond::ReadCondition(
            m_impl->initialPosition,
            dds::sub::status::DataState::new_data());

        waitPos.attach_condition(readCondPos);

        auto waitEuler = dds::core::cond::WaitSet();
        auto readCondEuler = dds::sub::cond::ReadCondition(
            m_impl->initialEuler,
            dds::sub::status::DataState::new_data());

        waitEuler.attach_condition(readCondEuler);

        casadi::DM pos, euler;
        try
        {
          waitPos.wait(dds::core::Duration::from_millisecs(m_impl->waitPosMs));
          auto samples = m_impl->initialPosition.read();

          // Set new desired value
          if(samples.length() > 0)
          {
            auto sample = (--samples.end());
            if(sample->info().valid())
            {
              auto sD = sample->data();
              pos = casadi::DM({sD.vec().x(), sD.vec().y(), sD.vec().z()});
              BOOST_LOG_TRIVIAL(trace) << "Initial position received: " << pos;
            }
          }
        }
        catch (dds::core::TimeoutError &)
        {
          BOOST_LOG_TRIVIAL(trace) << "Falling back to pre-defined initial position "
                                   << m_impl->fallbackInitPosition;
          pos = casadi::DM(std::vector<double>(m_impl->fallbackInitPosition));

        }
        try
        {
          waitEuler.wait(dds::core::Duration::from_millisecs(m_impl->waitEulerMs));
          auto samples = m_impl->initialEuler.read();

          // Set new desired value
          if(samples.length() > 0)
          {
            auto sample = (--samples.end());
            if(sample->info().valid())
            {
              auto sD = sample->data();
              euler = casadi::DM(std::vector<double>{sD.vec().y(), sD.vec().z()});
              BOOST_LOG_TRIVIAL(trace) << "Initial euler received: " << euler;
            }
          }
        }
        catch (dds::core::TimeoutError &)
        {
          BOOST_LOG_TRIVIAL(trace) << "Falling back to pre-defined initial euler";
          euler = casadi::DM(std::vector<double>{m_impl->fallbackInitEuler[1], m_impl->fallbackInitEuler[2]});
        }

        BOOST_LOG_TRIVIAL(trace) << "Pos: " << pos;
        BOOST_LOG_TRIVIAL(trace) << "Euler: " << euler;

        m_impl->x_k = casadi::DM::vertcat({pos, euler});

        m_next_step = std::chrono::steady_clock::now();
        event(new mimir::EvReady());
      }
      catch (...)
      {
        BOOST_LOG_TRIVIAL(fatal) << name() << " initializing threw exception";
        event(new mimir::EvError());
        throw;
      }
    }
    void FishSchool::timer(const std::atomic<bool>& cancel_token)
    {
      // probably want external time point to wait until.
      typedef std::chrono::milliseconds scm;
      auto now = std::chrono::steady_clock::now();
      auto time_p = m_next_step;
      auto timeleft = std::chrono::duration_cast<scm>(time_p - now).count();
      BOOST_LOG_TRIVIAL(trace) << "Time left: " << timeleft;

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

    void FishSchool::event(boost::statechart::event_base * const event)
    {
      // who is responsible for the pointer?
      m_scheduler.queue_event(
          m_stateMachine,
          make_intrusive(event));
    }


  }

}
