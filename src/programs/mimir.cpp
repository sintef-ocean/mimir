#include <atomic>
#include <cstdlib>
#include <cinttypes>
#include <chrono>
#include <exception>
#include <filesystem>
#include <functional>
#include <future>
#include <iostream>
#include <string>
#include <memory>
#include <thread>

#include <yaml-cpp/yaml.h>

#include "mimir/mimir.hpp"
#include "mimir/program/Convenience.hpp"
#include "mimir/program/Options.hpp"
#include "mimir/StateMachine.hpp"

#if defined(_WIN32) || defined(_WIN64)
int setenv(const char *name, const char *value, int overwrite)
{
  int errcode = 0;
  if(!overwrite) {
    size_t envsize = 0;
    errcode = getenv_s(&envsize, NULL, 0, name);
    if(errcode || envsize) return errcode;
  }
  return _putenv_s(name, value);
}

#include <windows.h>
#include <stdio.h>
#include <signal.h>

  BOOL CtrlHandler(DWORD fdwCtrlType)
  {
    switch (fdwCtrlType)
    {
      // Handle the CTRL-C signal.
    case CTRL_C_EVENT:
      mimir::signaled = 1;
      return TRUE;
      // Pass other signals to the next handler.
    case CTRL_BREAK_EVENT:
      mimir::signaled = 1;
      return TRUE;
    default:
      return FALSE;
    }
  }
#endif

int main(int argc, char *argv[])
{

  std::signal(SIGTERM, mimir::SIGTERMHandler);
  std::signal(SIGINT, mimir::SIGTERMHandler);
#if defined(_WIN32) || defined(_WIN64)
  SetConsoleCtrlHandler((PHANDLER_ROUTINE) CtrlHandler, TRUE);
#endif

  mimir::Options options("Purse seine state estimator");
  YAML::Node config;

  // Load configuration
  try
  {
    if(!options.parse(argc, argv))
      return EXIT_FAILURE;

    boost::log::add_console_log(std::cout, boost::log::keywords::auto_flush = true);
    boost::log::core::get()->set_filter(
        boost::log::trivial::severity >=
        mimir::transform_severity(options.severity) );

    config = YAML::LoadFile(options.filename); // load the file with YAML, (future: check schema with YAVL)

#ifdef MIMIR_INSTALL_PREFIX
    // This is not relocatable, but works if installed in default location provided by CPack installer
    int do_overwrite = 0;
    setenv("OSPL_HOME", MIMIR_INSTALL_PREFIX, do_overwrite);
    BOOST_LOG_TRIVIAL(debug) << "OSPL_HOME: " << getenv("OSPL_HOME") << std::endl;
    char prefix[] = MIMIR_INSTALL_PREFIX;
    namespace fs = std::filesystem;
    fs::path ospl_path = "file://";
    ospl_path += prefix;
    ospl_path /= "etc/config/ospl.xml";
    setenv("OSPL_URI", ospl_path.string().c_str(), do_overwrite);
    BOOST_LOG_TRIVIAL(debug) << "OSPL_URI: " << getenv("OSPL_URI") << std::endl;
#endif

    std::string algorithm = config["algorithm"]["name"].as<std::string>();

    if(algorithm.compare("PursePlanner") == 0)
      BOOST_LOG_TRIVIAL(info) << mimir::PlannerBanner();
    else if(algorithm.compare("Leadline") == 0)
      BOOST_LOG_TRIVIAL(info) << mimir::LeadlineBanner();
    else if(algorithm.compare("FishSchool") == 0)
      BOOST_LOG_TRIVIAL(info) << mimir::FishSchoolBanner();
    else if(algorithm.compare("KinematicVessel") == 0)
      BOOST_LOG_TRIVIAL(info) << mimir::KinematicVesselBanner();
  }
  catch (boost::program_options::error& e)
  {
    BOOST_LOG_TRIVIAL(fatal) << e.what();
    return EXIT_FAILURE;
  }
  catch (const YAML::ParserException& e)
  {
    BOOST_LOG_TRIVIAL(fatal) << "YAML parser exception: " + std::string(e.what());
    return EXIT_FAILURE;
  }

  std::thread stateMachineThread;
  int status = EXIT_SUCCESS;
  // Start main execution
  try
  {

    std::atomic<bool> machine_destructed(false);
    mimir::bsc::fifo_scheduler<> scheduler(true);
    mimir::bsc::fifo_scheduler<>::processor_handle stateMachine;
    std::promise<void> the_promise;

    stateMachine = scheduler.create_processor<mimir::StateMachine>(config, std::ref(machine_destructed));
    scheduler.initiate_processor(stateMachine);

    // if autostart:
    //    scheduler.queue_event(stateMachine, make_intrusive(new mimir::EvStart()));

    stateMachineThread = std::thread(
        [&the_promise, &scheduler, &machine_destructed]
        {
          try
          {
            scheduler(0);
            the_promise.set_value();
          }
          catch (...)
          {
            machine_destructed = true;
            the_promise.set_exception(std::current_exception());
          }
        });

    std::future<void> the_future = the_promise.get_future();

    while(!mimir::signaled && !machine_destructed)
    {
      std::this_thread::sleep_for(std::chrono::milliseconds(1000));
      if(mimir::signaled){ BOOST_LOG_TRIVIAL(info) << "Ctrl-C"; }
    }

    if (mimir::signaled)
    {
      BOOST_LOG_TRIVIAL(info) << "Posting Kill command to State Machine";
      scheduler.queue_event(stateMachine, make_intrusive(new mimir::EvKill()));
    }

    the_future.get(); // handles exceptions

  }
  catch (const YAML::Exception& e)
  {
    BOOST_LOG_TRIVIAL(fatal) << "YAML general exception: " + std::string(e.what());
    status = EXIT_FAILURE;
  }
  catch (const std::runtime_error& e)
  {
    BOOST_LOG_TRIVIAL(fatal) << "std::runtime_error: " + std::string(e.what());
    status = EXIT_FAILURE;
  }
  catch (const std::exception& e) {
    BOOST_LOG_TRIVIAL(fatal) << "std::exception: " + std::string(e.what());
    status = EXIT_FAILURE;
  }
  catch (...) {
    BOOST_LOG_TRIVIAL(fatal) << "Unknown error!";
    status = EXIT_FAILURE;
  }

  BOOST_LOG_TRIVIAL(debug) << "Waiting for StateMachine thread destruction";
  stateMachineThread.join();
  BOOST_LOG_TRIVIAL(debug) << "Program is shutting down";
  return status;

}
