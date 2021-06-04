#include <chrono>
#include <cmath>
#include <exception>
#include <thread>
#include <vector>

#include <boost/log/trivial.hpp>

#include "mimir/algorithm/Leadline.hpp"
#include "mimir/StateMachineFwd.hpp"
#include "mimir/FkinDds.hpp"

#include <casadi/casadi.hpp>

namespace mimir
{
  namespace algorithm
  {
    class Leadline::Impl
    {
    public:
      Impl(
          const YAML::Node& config,
          dds::pub::Publisher publisher,
          dds::sub::Subscriber subscriber) :
        inputParameters(dds::core::null),
        outputWriter(dds::core::null)
      {

        // Setup input reader for parameters
        auto input_config = config["inputs"]["parameters"];
        auto readerQos = subscriber.default_datareader_qos();
        auto topicName = input_config["topic"].as<std::string>();
        auto id = input_config["id"].as<std::string>();
        auto topic = dds::topic::Topic<fkin::IdVec2d>(
            subscriber.participant(),
            topicName);
        auto filter = dds::topic::Filter("id = %0", {id});

        inputParameters = dds::sub::DataReader<fkin::IdVec2d>(
            subscriber,
            dds::topic::ContentFilteredTopic<fkin::IdVec2d>(
                topic,
                topicName + id,
                filter),
            readerQos);

        parameters = casadi::DM(config["inputs"]["parameters"]["default"].as<std::vector<double>>());

        // Setup output writer
        auto writerQos = publisher.default_datawriter_qos();

        outputWriter = dds::pub::DataWriter<fkin::BatchIdVec1d>(
            publisher,
            dds::topic::Topic<fkin::BatchIdVec1d>(
                publisher.participant(),
                config["outputs"]["depth"]["topic"].as<std::string>()),
            writerQos);

        outputId = config["outputs"]["depth"]["id"].as<std::string>();

        // Setup model description
        using namespace casadi;
        SX x = SX::sym("x", 1);
        SX tau = SX::sym("tau", 1);
        SX x_d = SX::sym("set_point", 1);
        SX param = SX::vertcat({tau,x_d});
        SX rhs = SX::zeros(1,1);
        rhs = -(x - x_d)/tau;

        SXDict ode_int = {{"x", x}, {"p", param}, {"ode", rhs}};

        double predHorizonSec = config["prediction_horizon_sec"].as<double>();
        double timeStepSec = config["time_step_ms"].as<double>()/1000.;

        int32_t steps = static_cast<int32_t>(std::ceil(predHorizonSec/timeStepSec));
        timeGrid = std::vector<double>(steps+1);
        prediction = std::vector<fkin::IdVec1d>(timeGrid.size());
        predictionTime = std::vector<fkin::Timestamp>(timeGrid.size());

        for (size_t i = 0; i < timeGrid.size(); ++i)
          timeGrid[i] = 0. + i*timeStepSec;

        ode = casadi::integrator(
            "leadline", "cvodes", ode_int,
            {{"tf", predHorizonSec}, {"grid", timeGrid}});

        xGrid = casadi::DM::zeros(1,timeGrid.size()-1);

      }
      ~Impl() {}
      dds::sub::DataReader<fkin::IdVec2d> inputParameters;
      dds::pub::DataWriter<fkin::BatchIdVec1d> outputWriter;
      std::chrono::steady_clock::time_point t0;
      std::vector<double> timeGrid;
      std::vector<fkin::IdVec1d> prediction;
      std::vector<fkin::Timestamp> predictionTime;
      int32_t waitMs;
      std::string outputId;
      casadi::Function ode;
      casadi::DM xGrid;
      casadi::DM parameters;
    };


    Leadline::Leadline(
        const YAML::Node& config,
        boost::statechart::fifo_scheduler<>& scheduler,
        boost::statechart::fifo_scheduler<>::processor_handle machine,
        dds::pub::Publisher publisher,
        dds::sub::Subscriber subscriber) :
      m_impl( new Leadline::Impl(config, publisher, subscriber)),
      m_scheduler(scheduler),
      m_stateMachine(machine),
      m_time_step(std::chrono::milliseconds(config["time_step_ms"].as<std::int32_t>())),
      m_next_step(std::chrono::steady_clock::now()),
      m_config(config)
    {}

    Leadline::~Leadline() = default;
    void Leadline::solve(const std::atomic<bool>& cancel_token)
    {
      try
      {
        BOOST_LOG_TRIVIAL(trace) << "Leadline solve..";

        // Fetch current input, sample and hold
        auto params = m_impl->inputParameters.read();


        // Set new desired value
        if(params.length() > 0)
        {
          auto sample = (--params.end());
          if(sample->info().valid())
          {
            auto sD = sample->data();
            BOOST_LOG_TRIVIAL(debug) << "Got new params: "
                                    << sD.vec().x() << ", " << sD.vec().y();
            m_impl->parameters(0) = sD.vec().x();
            m_impl->parameters(1) = sD.vec().y();
          }
        }

        if(cancel_token)
        {
          BOOST_LOG_TRIVIAL(info) << name() << " interrupted";
          event(new mimir::EvInterrupt());
          return;
        }

        BOOST_LOG_TRIVIAL(debug) << "Parameters are : " << m_impl->parameters;

        // Integrate
        auto result = m_impl->ode(
            casadi::DMDict( {{"x0", casadi::DM::zeros(1,1)}, {"p", m_impl->parameters}} ));


        if(cancel_token)
        {
          BOOST_LOG_TRIVIAL(info) << name() << " interrupted";
          event(new mimir::EvInterrupt());
          return;
        }

        using namespace casadi;
        m_impl->xGrid = result["xf"];

        const auto x_k = std::vector<double>{0};
        namespace sc = std::chrono;

        auto t0 = m_now; // system clock of time at simulation start

        m_impl->prediction[0] = fkin::IdVec1d(m_impl->outputId, fkin::Vector1d(x_k[0]));
        m_impl->predictionTime[0] = fkin::Timestamp(
            sc::duration_cast<sc::milliseconds>(t0.time_since_epoch()).count());

        for(int i = 0; i < static_cast<int>(m_impl->timeGrid.size()-1); ++i)
        {
          m_impl->prediction[i+1] = fkin::IdVec1d(
              m_impl->outputId,
              fkin::Vector1d(std::vector<double>(m_impl->xGrid(Slice(),Slice(i,i+1)))[0]));
          m_impl->predictionTime[i+1] = fkin::Timestamp(
              sc::duration_cast<sc::milliseconds>(
                  (t0 + (i+1)*m_time_step).time_since_epoch()).count());
        }

        BOOST_LOG_TRIVIAL(trace) << "Response: " << m_impl->xGrid;

        m_impl->outputWriter << fkin::BatchIdVec1d(
            m_impl->outputId,
            m_impl->prediction,
            m_impl->predictionTime);

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
    void Leadline::initialize(const std::atomic<bool>&)
    {
      try
      {
        BOOST_LOG_TRIVIAL(debug) << name() << " is initializing";
        m_next_step = std::chrono::steady_clock::now();
        m_now = std::chrono::system_clock::now() - m_time_step;
        event(new mimir::EvReady());
      }
      catch (...)
      {
        BOOST_LOG_TRIVIAL(fatal) << name() << " initializing threw exception";
        event(new mimir::EvError());
        throw;
      }
    }

    void Leadline::timer(const std::atomic<bool>& cancel_token)
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

    void Leadline::event(boost::statechart::event_base * const event)
    {
      // who is responsible for the pointer?
      m_scheduler.queue_event(
          m_stateMachine,
          make_intrusive(event));
    }
  }
}
