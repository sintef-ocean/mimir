#include <chrono>
#include <cmath>
#include <exception>
#include <thread>

#include <boost/log/trivial.hpp>
#include <casadi/casadi.hpp>

#include "mimir/algorithm/TestAlgorithm.hpp"
#include "mimir/StateMachineFwd.hpp"
#include "mimir/FkinDds.hpp"

// ReaderLifespan
#include "dds/core/detail/ProprietaryApi.hpp"

namespace mimir
{
  namespace algorithm
  {
    class TestAlgorithm::Impl
    {
    public:
      Impl(
          const YAML::Node& config,
          dds::pub::Publisher publisher,
          dds::sub::Subscriber subscriber) :
        inputReader(dds::core::null),
        outputWriter(dds::core::null)
      {

        // Setup input reader
        const auto test_config = config["inputs"]["signal"];
        auto test_age = test_config["max_age_ms"].as<std::int64_t>();

        auto readerQos = subscriber.default_datareader_qos();
        //readerQos << dds::core::policy::Durability::TransientLocal()
        //          << dds::core::policy::Reliability::Reliable()
        readerQos << org::opensplice::core::policy::ReaderLifespan(
                      true,
                      dds::core::Duration::from_millisecs(test_age));

        const auto topicName = test_config["topic"].as<std::string>();
        const auto id = test_config["id"].as<std::string>(); // For filter topic
        auto topic = dds::topic::Topic<fkin::IdVec1d>(
            subscriber.participant(),
            topicName);
        auto filter = dds::topic::Filter("id = %0", {id});

        inputReader = dds::sub::DataReader<fkin::IdVec1d>(
            subscriber,
            dds::topic::ContentFilteredTopic<fkin::IdVec1d>(
                topic,
                topicName + id,
                filter),
            readerQos);

        // Setup output writer
        auto writerQos = publisher.default_datawriter_qos();
        //writerQos << dds::core::policy::Durability::TransientLocal();

        outputWriter = dds::pub::DataWriter<fkin::IdVec1d>(
            publisher,
            dds::topic::Topic<fkin::IdVec1d>(
                publisher.participant(),
                config["outputs"]["response"]["topic"].as<std::string>()),
            writerQos);

        outputId = config["outputs"]["response"]["id"].as<std::string>();

        // Setup model description
        casadi::SX x = casadi::SX::sym("x", 1);     // state
        casadi::SX x_d = casadi::SX::sym("x_d", 1); // set-point
        casadi::SX T = casadi::SX::sym("T", 1);     // time constant
        casadi::SX rhs = -(x-x_d)/T;                // prop-ctrl with time constant T
        casadi::SXDict ode = {{"x", x}, {"p", vertcat(x_d,T)}, {"ode", rhs}};

        integrator = casadi::integrator(
            "F", "cvodes", ode,
            {{"tf", config["time_step_ms"].as<double>()/1000.}});

        p_vec = casadi::DM(std::vector<double>{5.,5.});

      }
      ~Impl() {}
      dds::sub::DataReader<fkin::IdVec1d> inputReader;
      dds::pub::DataWriter<fkin::IdVec1d> outputWriter;
      std::string outputId;
      casadi::Function integrator;
      casadi::DM p_vec;
      casadi::DM x_k;
    };

    TestAlgorithm::TestAlgorithm(
        const YAML::Node& config,
        boost::statechart::fifo_scheduler<> & scheduler,
        boost::statechart::fifo_scheduler<>::processor_handle machine,
        dds::pub::Publisher publisher,
        dds::sub::Subscriber subscriber) :
      m_impl( new TestAlgorithm::Impl(config, publisher, subscriber)),
      m_scheduler(scheduler),
      m_stateMachine(machine),
      m_time_step(std::chrono::milliseconds(config["time_step_ms"].as<std::int32_t>())),
      m_next_step(std::chrono::steady_clock::now()),
      m_config(config)
    {}
    void TestAlgorithm::solve(const std::atomic<bool>& cancel_token)
    {
      try
      {
        BOOST_LOG_TRIVIAL(trace) << "Solving " << name();

        // Publish x_k(t_now) (we do not publish future samples)
        m_impl->outputWriter << fkin::IdVec1d(
            m_impl->outputId,
            fkin::Vector1d(m_impl->x_k.get_nonzeros()[0]));

        // Fetch input, sample and hold
        auto samples = m_impl->inputReader.read();
        // waitset with short timeout, allowing run in "open loop"

        // Set new desired value
        if(samples.length() > 0)
        {
          auto sample = (--samples.end());
          if(sample->info().valid())
          {
            m_impl->p_vec(0,0) = sample->data().vec().x();
            BOOST_LOG_TRIVIAL(trace) << "Got a sample: " << sample->data().vec().x();
          }
        }
        else
        {
          // TODO: In case of measurements: Disable feedback gain
          BOOST_LOG_TRIVIAL(trace) << "No valid sample for " << name() << ", using old";
        }

        if(cancel_token)
        {
          event(new mimir::EvInterrupt());
          return;
        }

        // Integrate
        auto result = m_impl->integrator(
            casadi::DMDict( {{"x0",m_impl->x_k}, {"p", m_impl->p_vec}} ));

        if(cancel_token)
        {
          event(new mimir::EvInterrupt());
          return;
        }

        // Update output (not published yet)
        m_impl->x_k = result["xf"];
        step_time();
        event(new mimir::EvReady());

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

    void TestAlgorithm::initialize(const std::atomic<bool>&)
    {
      try
      {
        BOOST_LOG_TRIVIAL(debug) << name() << " is initializing";

        // waitset to get initial measurements that meet requirements.
        // set initial state vector..
        m_impl->x_k = casadi::DM(std::vector<double>{0.});
        // set parameters as well?
        m_next_step = std::chrono::steady_clock::now();

        event(new mimir::EvReady());
      }
      catch (...)
      {
        event(new mimir::EvError());
        throw;
      }
    }

    void TestAlgorithm::timer(const std::atomic<bool>& cancel_token)
    {

      // There is no sanity checks in terms of real-timeliness of execution.
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

    void TestAlgorithm::event(boost::statechart::event_base * const event)
    {
      // who is responsible for the pointer?
      m_scheduler.queue_event(
          m_stateMachine,
          make_intrusive(event));
    }

    TestAlgorithm::~TestAlgorithm() = default;
  }
}
