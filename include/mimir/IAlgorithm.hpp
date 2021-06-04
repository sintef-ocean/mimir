#pragma once
/**
   @file IAlgorithm.hpp
   @brief Algorithm interface file. The algorithm is created by an algorithm factory as an instance held by the state machine.
*/
#include <atomic>

namespace mimir
{
  /**
     @brief Algorithm interface class.

     A mimir::StateMachine holds an instance of a class that is inherited from IAlgorithm.
     The interface functions are the only directly exposed interactions with the
     algorithm. These functions drives the state machine under normal circumstances by
     post appropriate events to the state machine. Each function, except name(), is an
     asynchronous call using the std::future mechanism. The std::atomic<bool> parameter
     is an external boolean whether the algorithm should cancel its workload.

  */
  class IAlgorithm
  {
  public:
    /// Destructor.
    virtual ~IAlgorithm() {};
    /**
       @brief Function to run the actual algorithm. Called at the entry of mimir::Evaluating.

       Solves the algorithm problem. If the solution time is perceived as long,
       intermittent checks for cancel_token should be performed to provide graceful
       cancellation of the workload.

       + If solving fails, post an mimir::EvError event.
       + When evaluation is successful, post mimir::EvReady.
       + If cancel_token is true, post mimir::EvInterrupt.

       @param [in] cancel_token atomic boolean whether an external cancel has been issued.
    */
    virtual void solve(const std::atomic<bool>& cancel_token) = 0;
    /**
       @brief Function called when entering the mimir::Initializing state of the mimir::StateMachine.

       Prepares the algorithm, e.g. by retrieving input data like initial conditions from
       DDS.

       + If cancel_token is true, the function shall post a mimir::EvInterrupt event. The
       boolean should be checked regularly.
       + If initialization fails, post an mimir::EvError event.
       + If initialization is done, post mimir::EvReady.

       @param [in] cancel_token atomic boolean whether an external cancel has been issued.
    */
    virtual void initialize(const std::atomic<bool>& cancel_token) = 0;
    /**
       @brief Timer function that posts mimir::EvTimeout when the time duration is spent.

       This function is called when entering mimir::Waiting of mimir::StateMachine. How
       the timer is implemented is up to the implementor. Possible patterns can be regular
       intervals, state/logic-dependent durations. A common use case is that the timer
       waits the remainder of a fixed interval, for which solve() has spent the initial
       part.

       The function shall contain the sleep logic, using
       e.g. std::this_thread::sleep_for(), and post appropriate events:

       + If cancel_token is true, the function shall post a mimir::EvInterrupt event. The
       boolean should be checked regularly.
       + If the algorithm-specified time duration has passed, post a mimir::EvTimeout.

       @param [in] cancel_token atomic boolean whether an external cancel has been issued.
    */
    virtual void timer(const std::atomic<bool>& cancel_token) = 0;
    /// Name identifier string for the implemented algorithm.
    virtual const char* name() = 0;
  };
}
