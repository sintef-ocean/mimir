#include "mimir/algorithm/AlgorithmFactory.hpp"
#include "boost/log/trivial.hpp"

namespace mimir
{

  std::unique_ptr<mimir::IAlgorithm> AlgorithmCreator(
      const std::string& id,
      const YAML::Node& config,
      boost::statechart::fifo_scheduler<> & scheduler,
      boost::statechart::fifo_scheduler<>::processor_handle machine,
      dds::pub::Publisher publisher,
      dds::sub::Subscriber subscriber)
  {
    if(id.compare("TestAlgorithm") == 0)
    {
      BOOST_LOG_TRIVIAL(trace) << "Creating TestAlgorithm instance";
      return std::make_unique<mimir::algorithm::TestAlgorithm>(
          config["TestAlgorithm"],
          scheduler,
          machine,
          publisher,
          subscriber);
    }
    else if(id.compare("PursePlanner") == 0)
    {
      BOOST_LOG_TRIVIAL(trace) << "Creating PursePlanner algorithm instance";

      return std::make_unique<mimir::algorithm::PursePlanner>(
          config["PursePlanner"],
          scheduler,
          machine,
          publisher,
          subscriber);
    }
    else if(id.compare("KinematicVessel") == 0)
    {
      BOOST_LOG_TRIVIAL(trace) << "Creating KinematicVessel algorithm instance";
      return std::make_unique<mimir::algorithm::KinematicVessel>(
          config["KinematicVessel"],
          scheduler,
          machine,
          publisher,
          subscriber);
    }
    else if(id.compare("FishSchool") == 0)
    {
      BOOST_LOG_TRIVIAL(trace) << "Creating FishSchool algorithm instance";
      return std::make_unique<mimir::algorithm::FishSchool>(
          config["FishSchool"],
          scheduler,
          machine,
          publisher,
          subscriber);
    }
    else if(id.compare("Leadline") == 0)
    {
      BOOST_LOG_TRIVIAL(trace) << "Creating Leadline algorithm instance";
      return std::make_unique<mimir::algorithm::Leadline>(
          config["Leadline"],
          scheduler,
          machine,
          publisher,
          subscriber);
    }
    else if(id.compare("Schybert") == 0)
    {
      BOOST_LOG_TRIVIAL(trace) << "Creating Schybert algorithm instance";
      throw std::logic_error("Shybert not yet implemented");
      //return std::unique_ptr<mimir::IAlgorithm>(nullptr);
    }
    else
      throw std::runtime_error(
          std::string("AlgorithmCreator got unknown algorithm identifier: ") + id);
  }
}
