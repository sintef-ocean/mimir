#pragma once
/**
   @file StateNotifier.hpp
   @brief State notifier class for broadcasting state machine state with DDS.
*/

#include <memory>
#include <boost/log/trivial.hpp>
#include "mimir/FkinDds.hpp"

namespace mimir
{
  namespace detail
  {
    /// Convenience function for creating the StateNootifier's default Quality of Service.
    inline dds::pub::qos::DataWriterQos NotifierQos(const dds::pub::Publisher& p)
    {
      auto qos = p.default_datawriter_qos();
      qos << dds::core::policy::Durability::TransientLocal();
      // should TransientLocal be removed?
      return qos;
    }
  }

  namespace control
  {

    /**
       @brief State machine state notifier

       With this class a state machine's state can be broadcast to other applications.
       This class is a member of a mimir::StateMachine and notifies over DDS which state
       it enters as an fkin::ProcessStateAutomaton.

    */
    class StateNotifier
    {
    public:
      /**
         @brief Constructor.

         The constructor sets up a DDS data writer for notifications.

         @param [in] notifyTopicName DDS topic to publish notifications
         @param [in] notifyIdentifier DDS key for which remote subscribers uses to identify process
         @param [in] publisher DDS publisher reference
      */
      StateNotifier(
          const std::string& notifyTopicName,
          const std::string& notifyIdentifier,
          dds::pub::Publisher publisher) :
        m_notifyWriter(dds::core::null),
        m_identifier(notifyIdentifier)
      {
        try
        {
          m_notifyWriter = dds::pub::DataWriter<fkin::ProcessStateAutomaton>(
              publisher,
              dds::topic::Topic<fkin::ProcessStateAutomaton>(
                  publisher.participant(),
                  notifyTopicName),
              detail::NotifierQos(publisher));
        }
        catch (const dds::core::Exception& e) {
          BOOST_LOG_TRIVIAL(fatal) << "DDS exception: " + std::string(e.what());
          throw std::runtime_error("StateNotifier failed to construct");
        }
      }

      StateNotifier(const StateNotifier&) = delete;
      StateNotifier& operator=(const StateNotifier&) = delete;
      StateNotifier(StateNotifier&&) = default;
      StateNotifier& operator=(StateNotifier&&) = default;
      ~StateNotifier() = default;

      /**
         @brief State notifier function that writes data on DDS publisher.
         @param [in] state Enum of state to be published.
      */
      void NotifyState(fkin::ProcessStateKind state)
      {
        // Notify state
        auto notification = fkin::ProcessStateAutomaton(m_identifier, state);
        m_notifyWriter << notification;
      }
    private:

      StateNotifier() = delete;
      dds::pub::DataWriter<fkin::ProcessStateAutomaton> m_notifyWriter; ///< DDS data writer.
      std::string m_identifier; ///< String describer for which a remote subscriber uses to identify the process.
    };


  }
}
