#pragma once
/**
   @file CommandResponder.hpp
   @brief Command request and reply over DDS to control remote applications.
*/

#include <boost/log/trivial.hpp>
#include "mimir/StateMachineFwd.hpp"

#include <memory>

#include "mimir/FkinDds.hpp"

// Note: This functionality is heavily inspired/copied from odss::reqrep
namespace mimir
{
  /// Provides implementation details in the form of helper functions.
  namespace detail
  {
    /// Convenience function for creating the default QoS for the CommandResponder's request reader.
    inline dds::sub::qos::DataReaderQos ReqReaderQos(const dds::sub::Subscriber& s)
    {
      auto qos = s.default_datareader_qos();
      qos << dds::core::policy::Durability::TransientLocal()
          << dds::core::policy::Reliability::Reliable();
      return qos;
    }
    /// Convenience function for creating the default QoS for the CommandResponder's reply writer.
    inline dds::pub::qos::DataWriterQos RepWriterQos(const dds::pub::Publisher& p)
    {
      auto qos = p.default_datawriter_qos();
      qos << dds::core::policy::Durability::TransientLocal();
      // should TransientLocal be removed?
      return qos;
    }
  }

  /// Classes that allow user control of processes, injecting events into state machine.
  namespace control
  {

    /// Listener class that calls a callback function when subscribed DDS data becomes available.
    template <typename T>
    class CommandListener :
      public dds::sub::NoOpDataReaderListener<T>
    {
    public:
      /**
         @brief Constructor.

         Data reader listener for DDS type T. When data becomes available, the
         user-provided function is called.

         @param [in] doHandle Function pointer to call, when on_data_available().
      */
      CommandListener(
          std::function<void(dds::sub::DataReader<T>&)> doHandle)
      {
        m_doHandleFcn = std::move(doHandle);
      }

      /**
         @brief Function overridden from base class to be called on data available.

         This function is called when subscribed data becomes available on the DDS bus.
         It dispatches a reference of the data reader to the handle function that the user
         provided at construction.

      */
      virtual void on_data_available(dds::sub::DataReader<T>& dataReader)
      {
        if(m_doHandleFcn)
          m_doHandleFcn(dataReader);
      }

    private:
      std::function<void(dds::sub::DataReader<T>&)> m_doHandleFcn; ///< Function object to be called.
    };


    /**
       @brief Request and reply class for user commands.

       This class ensures that commands received over DDS from other applications are
       posted as event in the specified state machine.

    */
    class CommandResponder
    {
    public:

      /**
         @brief Constructor that sets up DDS readers, writers, listeners.

         This function constructs the following:
         + A topic requestTopicName, which is filtered on the recipient key.
         + A reader for command requests of this filtered topic.
         + A writer for command responses to replyTopicName.
         + A command listener callback to post event into state machine.

         @param [in] requestTopicName Topic for requests
         @param [in] replyTopicName Topic for responses
         @param [in] recipient Recipient key
         @param [in] publisher Reference to the publisher
         @param [in] subscriber Reference to the subscriber
         @param [in] scheduler Reference to the state machine scheduler
         @param [in] machine Reference to the state machine handle
      */
      CommandResponder(
          const std::string& requestTopicName,
          const std::string& replyTopicName,
          const std::string& recipient,
          dds::pub::Publisher publisher,
          dds::sub::Subscriber subscriber,
          boost::statechart::fifo_scheduler<> & scheduler,
          boost::statechart::fifo_scheduler<>::processor_handle machine) :
        m_scheduler(scheduler),
        m_stateMachine(machine),
        m_requestReader(dds::core::null),
        m_replyWriter(dds::core::null)
      {

        try
        {
          dds::topic::Topic<fkin::Command> requestTopic(
              subscriber.participant(),
              requestTopicName);
          dds::topic::Filter filter("header.recipient = %0", {recipient});
          dds::topic::ContentFilteredTopic<fkin::Command> filteredTopic(
              requestTopic, recipient + requestTopicName, filter);

          m_requestReader = dds::sub::DataReader<fkin::Command>(
              subscriber,
              filteredTopic,
              detail::ReqReaderQos(subscriber));

          m_replyWriter = dds::pub::DataWriter<fkin::CommandResponse>(
              publisher,
              dds::topic::Topic<fkin::CommandResponse>(
                  publisher.participant(),
                  replyTopicName),
              detail::RepWriterQos(publisher));

          // clear reader queue of any old samples
          m_requestReader.wait_for_historical_data(dds::core::Duration::infinite());
          m_requestReader.take();

          m_listener = std::make_unique<CommandListener<fkin::Command>>(
              [=](dds::sub::DataReader<fkin::Command>&){ handleResponse(); });
          m_requestReader.listener(
              m_listener.get(),
              dds::core::status::StatusMask::data_available());
        }
        catch (const dds::core::Exception& e) {
          BOOST_LOG_TRIVIAL(fatal) << "DDS exception: " + std::string(e.what());
          throw std::runtime_error("CommandResponder failed to contruct");
        }

      }

      CommandResponder(const CommandResponder&) = delete;
      CommandResponder& operator=(const CommandResponder&) = delete;
      CommandResponder(CommandResponder&&) = default;
      CommandResponder& operator=(CommandResponder&&) = default;

      /// Destructor. Resets listener
      ~CommandResponder()
      {
        m_requestReader.listener(nullptr, dds::core::status::StatusMask::none());
      }

    private:
      /**
         @brief Callback function for command responses

         This function is used as the callback function in the CommandListener member variable.
         It parses the incoming command and posts appropriate events in the state machine.

      */
      void handleResponse()
      {

        const auto samples = m_requestReader.select()
         .state(dds::sub::status::DataState::new_data())
         .take();

        std::for_each(
            samples.begin(),
            samples.end(),
            [=](const dds::sub::Sample<fkin::Command>& sample)
            {
              std::string cmd("Noop");
              if (sample.info().valid())
              {
                switch(sample.data().command())
                {
                case fkin::CommandType::START_PROCESS:
                  cmd = "Started " + sample.data().header().recipient();
                  m_scheduler.queue_event(
                      m_stateMachine,
                      make_intrusive(new mimir::EvStart()));
                  break;
                case fkin::CommandType::STOP_PROCESS:
                  cmd = "Stopped " + sample.data().header().recipient();
                  m_scheduler.queue_event(
                      m_stateMachine,
                      make_intrusive(new mimir::EvStop()));
                  break;
                case fkin::CommandType::TERMINATE_PROCESS:
                  cmd = "Terminated " + sample.data().header().recipient();
                  m_scheduler.queue_event(
                      m_stateMachine,
                      make_intrusive(new mimir::EvKill()));
                  break;
                default:
                  ;
                }
              }

              // Send response
              auto response = fkin::CommandResponse(
                  fkin::ReplyHeader(sample.data().header().requestID()),
                  true,
                  cmd);

              m_replyWriter << response;
            });
      }

      CommandResponder() = delete;
      boost::statechart::fifo_scheduler<> & m_scheduler; ///< Scheduler for state machine.
      boost::statechart::fifo_scheduler<>::processor_handle m_stateMachine; ///< StateMachine reference.
      dds::sub::DataReader<fkin::Command> m_requestReader; ///< DDS reader for requests.
      dds::pub::DataWriter<fkin::CommandResponse> m_replyWriter; ///< DDS writer for command responses.
      std::unique_ptr<CommandListener<fkin::Command>> m_listener; ///< Command listener callback.

    };


  }
}
