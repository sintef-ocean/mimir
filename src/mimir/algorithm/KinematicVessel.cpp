#include <chrono>
#include <cmath>
#include <exception>
#include <thread>
#include <vector>

#include <boost/log/trivial.hpp>

#include "mimir/algorithm/KinematicVessel.hpp"
#include "mimir/StateMachineFwd.hpp"
#include "mimir/FkinDds.hpp"

#include <casadi/casadi.hpp>

namespace mimir
{
  namespace algorithm
  {
    class KinematicVessel::Impl
    {
    public:
      Impl(
          const YAML::Node& config,
          dds::pub::Publisher publisher,
          dds::sub::Subscriber subscriber) :
        inputReader(dds::core::null),
        outputWriter(dds::core::null),
        initialConditionReader(dds::core::null),
        fallbackInitCond(
            config["initial_conditions"]["position_course"]["fallback"]
            .as<std::vector<double>>())
      {

        // Setup input reader
        const auto input_config = config["inputs"]["vessel_ctrl"];
        auto input_age = input_config["max_age_ms"].as<std::int64_t>();

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


        outputWriter = dds::pub::DataWriter<fkin::Kinematics2D>(
            publisher,
            dds::topic::Topic<fkin::Kinematics2D>(
                publisher.participant(),
                config["outputs"]["kinematics"]["topic"].as<std::string>()),
            writerQos);

        outputId = config["outputs"]["kinematics"]["id"].as<std::string>();

        // Setup initial condition reader
        readerQos << org::opensplice::core::policy::ReaderLifespan(
            false,
            dds::core::Duration::infinite());

        const auto initCondConfig = config["initial_conditions"]["position_course"];
        const auto initTopicName = initCondConfig["topic"].as<std::string>();
        const auto idInit = initCondConfig["id"].as<std::string>();
        waitMs = initCondConfig["max_wait_ms"].as<std::int32_t>();
        topic = dds::topic::Topic<fkin::IdVec3d>(
            subscriber.participant(),
            initTopicName);
        filter = dds::topic::Filter("id = %0", {idInit});

        initialConditionReader = dds::sub::DataReader<fkin::IdVec3d>(
            subscriber,
            dds::topic::ContentFilteredTopic<fkin::IdVec3d>(
                topic,
                initTopicName + idInit,
                filter),
            readerQos);

        using namespace casadi;
        // Setup model description
        SX x = SX::sym("x", 3);
        SX u_s = SX::sym("u", 3);
        SX rhs = SX::zeros(3,3);

        rhs(0,0) = cos(x(2));
        rhs(1,0) = sin(x(2));
        rhs(0,1) = -sin(x(2));
        rhs(1,1) = cos(x(2));
        rhs(2,2) = SX(1);

        rhs = mtimes(rhs,u_s);

        BOOST_LOG_TRIVIAL(trace) << "RHS: " << rhs;

        SXDict ode = {{"x", x}, {"p", u_s}, {"ode", rhs}};

        model = casadi::integrator(
            "vessel", "cvodes", ode,
            {{"tf", config["time_step_ms"].as<double>()/1000.}});

        u = casadi::DM(std::vector<double>{1.,0.,0.});
      }
      ~Impl() {}
      dds::sub::DataReader<fkin::IdVec3d> inputReader;
      dds::pub::DataWriter<fkin::Kinematics2D> outputWriter;
      dds::sub::DataReader<fkin::IdVec3d> initialConditionReader;
      std::chrono::steady_clock::time_point t0;
      int32_t waitMs;
      std::vector<double> fallbackInitCond;
      std::string outputId;
      casadi::Function model;
      casadi::DM u;
      casadi::DM x_k;
    };


    KinematicVessel::KinematicVessel(
        const YAML::Node& config,
        boost::statechart::fifo_scheduler<>& scheduler,
        boost::statechart::fifo_scheduler<>::processor_handle machine,
        dds::pub::Publisher publisher,
        dds::sub::Subscriber subscriber) :
      m_impl( new KinematicVessel::Impl(config, publisher, subscriber)),
      m_scheduler(scheduler),
      m_stateMachine(machine),
      m_time_step(std::chrono::milliseconds(config["time_step_ms"].as<std::int32_t>())),
      m_next_step(std::chrono::steady_clock::now()),
      m_config(config)
      // Allow for simulation speed factor
    {}

    KinematicVessel::~KinematicVessel() = default;
    void KinematicVessel::solve(const std::atomic<bool>& cancel_token)
    {
      try
      {
        BOOST_LOG_TRIVIAL(trace) << "KinematicVessel solve..";

        // Publish x_k(t_now) (we do not publish future samples)
        const auto x_k = std::vector<double>(m_impl->x_k);
        m_impl->outputWriter << fkin::Kinematics2D(
            m_impl->outputId,
            fkin::Vector2d(x_k[0], x_k[1]),
            fkin::Vector1d(x_k[2]),
            fkin::Vector1d(std::vector<double>(m_impl->u)[0]));

        // Fetch input, sample and hold
        auto samples = m_impl->inputReader.read();

        // Set new desired value
        if(samples.length() > 0)
        {
          auto sample = (--samples.end());
          if(sample->info().valid())
          {
            auto sD = sample->data();
            m_impl->u = casadi::DM({sD.vec().x(), sD.vec().y(), sD.vec().z()});
            BOOST_LOG_TRIVIAL(trace) << "Got a sample: " << m_impl->u;
          }
        }
        else
        {
          BOOST_LOG_TRIVIAL(debug) << "No valid input sample for " << name() << ", using old";
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
    void KinematicVessel::initialize(const std::atomic<bool>&)
    {
      try
      {
        BOOST_LOG_TRIVIAL(debug) << name() << " is initializing";

        auto waitSet = dds::core::cond::WaitSet();
        auto readCondition = dds::sub::cond::ReadCondition(
            m_impl->initialConditionReader,
            dds::sub::status::DataState::new_data());

        waitSet.attach_condition(readCondition);

        try
        {
          waitSet.wait(dds::core::Duration::from_millisecs(m_impl->waitMs));

          auto samples = m_impl->initialConditionReader.read();

          // Set new desired value
          if(samples.length() > 0)
          {
            auto sample = (--samples.end());
            if(sample->info().valid())
            {
              auto sD = sample->data();
              m_impl->x_k = casadi::DM({sD.vec().x(), sD.vec().y(), sD.vec().z()});
              BOOST_LOG_TRIVIAL(trace) << "Initial condition received: " << m_impl->x_k;
            }
          }
        }
        catch (dds::core::TimeoutError &)
        {
          BOOST_LOG_TRIVIAL(trace) << "Falling back to pre-defined initial condition";
          m_impl->x_k = casadi::DM(m_impl->fallbackInitCond);
        }

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
    void KinematicVessel::timer(const std::atomic<bool>& cancel_token)
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

    void KinematicVessel::event(boost::statechart::event_base * const event)
    {
      // who is responsible for the pointer?
      m_scheduler.queue_event(
          m_stateMachine,
          make_intrusive(event));
    }


  }

}
