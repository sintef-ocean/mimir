#include "mimir/control/CommandResponder.hpp"
#include "mimir/control/StateNotifier.hpp"
#include "mimir/algorithm/AlgorithmFactory.hpp"
#include "mimir/StateMachine.hpp"

#include "mimir/FkinDds.hpp"

namespace mimir
{

  class StateMachine::DdsTriad
  {
  public:
    DdsTriad(std::uint32_t domain=0) :
      participant(domain),
      subscriber(participant),
      publisher(participant)
    {}
    ~DdsTriad() = default;
    dds::domain::DomainParticipant participant;
    dds::sub::Subscriber subscriber;
    dds::pub::Publisher publisher;
  };

  /*
  void StateMachine::initiate_impl()
  {
    bsc::state_machine<
      StateMachine,
      Standby,
      std::allocator<void>,
      bsc::exception_translator<>>::initiate();
  }*/

  StateMachine::StateMachine(
      my_context ctx,
      const YAML::Node& config,
      std::atomic<bool>& notify_others) :
    my_base(ctx),
    m_dds(new DdsTriad(config["dds"]["domain"].as<std::uint32_t>())),
    m_responder(new control::CommandResponder(
         config["algorithm"]["command"]["request_topic"].as<std::string>(),
         config["algorithm"]["command"]["reply_topic"].as<std::string>(),
         config["algorithm"]["command"]["recipient"].as<std::string>(),
         m_dds->publisher,
         m_dds->subscriber,
         my_scheduler(),
         my_handle())),
    m_notifier(new control::StateNotifier(
         config["algorithm"]["notifier"]["notify_topic"].as<std::string>(),
         config["algorithm"]["notifier"]["identifier"].as<std::string>(),
         m_dds->publisher)),
    m_notify_others(notify_others),
    m_config(config)
  {
    BOOST_LOG_TRIVIAL(trace) << "Entered : " << __FUNCTION__;
  }

  StateMachine::~StateMachine()
  {
    m_notify_others = true;
    std::this_thread::sleep_for(std::chrono::milliseconds(50));
    context< StateMachine >().Notifier()->NotifyState(fkin::ProcessStateKind::DEAD);
    BOOST_LOG_TRIVIAL(debug) << "Destructing: " << __FUNCTION__;
  }

  StateMachine::DdsTriad * StateMachine::Dds(){ return m_dds.get(); }


  Standby::Standby(my_context ctx) : my_base(ctx)
  {
    context< StateMachine >().Notifier()->NotifyState(fkin::ProcessStateKind::IDLE);
    BOOST_LOG_TRIVIAL(info) << "mimir::StateMachine in " << __FUNCTION__;
  }

  Standby::~Standby() = default;
  bsc::result Standby::react(const EvTimeout &)
  {
    return transit< Standby >();
  }

  ShutDown::ShutDown(my_context ctx) : my_base(ctx)
  {
    try
    {
      BOOST_LOG_TRIVIAL(trace) << "Entered state: " << __FUNCTION__;
      context< StateMachine >().Notifier()->NotifyState(fkin::ProcessStateKind::DEAD);
      post_event(EvKill());
    }
    catch (...)
    {
      post_event(EvKill());
    }
  }

  ShutDown::~ShutDown() = default;

  bsc::result ShutDown::react(const EvKill &)
  {
    BOOST_LOG_TRIVIAL(info) << "Killing StateMachine";
    outermost_context_type & machine = outermost_context();
    machine.my_scheduler().destroy_processor(machine.my_handle());
    machine.my_scheduler().terminate();
    return terminate();
  }

  Running::Running(my_context ctx) :
    my_base(ctx),
    m_algorithm( mimir::AlgorithmCreator(
     context< StateMachine >().Config()["algorithm"]["name"].as<std::string>(),
     context< StateMachine >().Config(),
     outermost_context().my_scheduler(),
     outermost_context().my_handle(),
     context< StateMachine >().Dds()->publisher,
     context< StateMachine >().Dds()->subscriber) )
  {
    context< StateMachine >().Notifier()->NotifyState(fkin::ProcessStateKind::RUNNING);
    BOOST_LOG_TRIVIAL(info)
     << "mimir::StateMachine in "
     << __FUNCTION__ << " with "
     << Algorithm()->name() << " algorithm";
  }

  Running::~Running() = default;

  bsc::result Running::react(const bsc::exception_thrown &)
  {
    // #TODO: add more exception handling here, or just propagate to main thread..
    try
    {
      throw;
    }
    catch (const dds::core::Exception& e)
    {
      BOOST_LOG_TRIVIAL(fatal) << "DDS exception: " + std::string(e.what());
    }
    catch (const YAML::Exception& e)
    {
      BOOST_LOG_TRIVIAL(fatal) << "YAML general exception: " + std::string(e.what());
    }
    catch (const std::runtime_error& e)
    {
      BOOST_LOG_TRIVIAL(fatal) << "std::runtime_exception: " + std::string(e.what());
    }
    catch (const std::exception& e)
    {
      BOOST_LOG_TRIVIAL(fatal) << "std::exception: " + std::string(e.what());
    }
    catch (...)
    {
      BOOST_LOG_TRIVIAL(fatal) << "Unknown error!";
    }
    return transit< ShutDown >();
  }

  mimir::IAlgorithm* Running::Algorithm(){ return m_algorithm.get(); }


  Initializing::Initializing(my_context ctx) : my_base(ctx), m_canceled(false)
  {
    BOOST_LOG_TRIVIAL(trace) << "Entered state: " << __FUNCTION__;
    context< StateMachine >().Notifier()->NotifyState(fkin::ProcessStateKind::INITIALIZING);
    m_future = std::async(
        [algorithm=context< Running >().Algorithm(),
         &cancel_token=m_canceled]
        { algorithm->initialize(cancel_token); });
  }

  Initializing::~Initializing()
  {
    BOOST_LOG_TRIVIAL(trace) << "Leaving " << __FUNCTION__;
    context< StateMachine >().Notifier()->NotifyState(fkin::ProcessStateKind::RUNNING);
    m_canceled = true;
  }

  bsc::result Initializing::react(const EvReady &)
  {
    if(m_future.valid())
      m_future.get();
    return transit< Evaluating >();
  }


  Evaluating::Evaluating(my_context ctx) : my_base(ctx), m_canceled(false)
  {
    BOOST_LOG_TRIVIAL(trace) << "Entered state: " << __FUNCTION__;
    m_future = std::async(
        [algorithm=context< Running >().Algorithm(),
         &cancel_token=m_canceled]
        { algorithm->solve(cancel_token); });
  }

  Evaluating::~Evaluating() {
    BOOST_LOG_TRIVIAL(trace) << "Leaving " << __FUNCTION__;
    m_canceled = true;
  }

  bsc::result Evaluating::react(const EvReady &)
  {
    BOOST_LOG_TRIVIAL(trace) << "Evaluation ready event";
    if(m_future.valid())
      m_future.get();
    return transit< Waiting >();
  }


  Waiting::Waiting(my_context ctx) : my_base(ctx), m_canceled(false)
  {
    BOOST_LOG_TRIVIAL(trace) << "Entered state: " << __FUNCTION__;
    m_future = std::async(
        [algorithm=context< Running >().Algorithm(),
         &cancel_token=m_canceled]
        { algorithm->timer(cancel_token); });
  }

  Waiting::~Waiting()
  {
    BOOST_LOG_TRIVIAL(trace) << "Leaving state: " << __FUNCTION__;
    m_canceled = true;
  }


}
