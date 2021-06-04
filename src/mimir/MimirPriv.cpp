#include <iostream>
#include <cinttypes>
#include "mimir/MimirPriv.hpp"
#include "mimir/mimir.hpp"


namespace mimir
{
#if defined(_WIN32) || defined(_WIN64)
  std::string banner_start("\n");
  std::string banner_stop("\n");
#else
  std::string banner_start("\n\033[34m");
  std::string banner_stop("\n\033[0m");
#endif

  MimirPriv::MimirPriv() {}
  void MimirPriv::solve()
  {
    auto x = m_opti.variable();
    auto y = m_opti.variable();
    auto z = m_opti.variable();

    casadi::Dict ipopt_dict, nlp_dict;
    ipopt_dict["print_level"] = 0;
    ipopt_dict["print_frequency_time"] = 1;
    // ipopt_dict["linear_solver"] = "ma57";

    nlp_dict["expand"] = true;

    m_opti.minimize( x*x + 100*z*z );
    m_opti.subject_to(z+(1-x)*(1-x)-y==0);
    m_opti.solver("ipopt", nlp_dict, ipopt_dict);
    auto sol = m_opti.solve();

    std::cout << sol.value(x) << ":" << sol.value(y) << std::endl;
    std::cout << sol.stats()["return_status"] << std::endl;

  }

  MimirClass::MimirClass(): m_priv{new MimirPriv()} {}
  MimirClass::~MimirClass(){ delete m_priv; }
  void MimirClass::solve(){
    m_priv->solve();
  }

  std::string KinematicVesselBanner()
  {

    return banner_start + std::string("\
   __ ___     ___ _  _        ___      ___  _\n\
  (_   |  |\\ | | |_ |_   |\\/|  |  |\\/|  |  |_)\n\
  __) _|_ | \\| | |_ |    |  | _|_ |  | _|_ | \\\n\
  ------------> KinematicVessel <-----------")
    + banner_stop;
  }

  std::string FishSchoolBanner()
  {
    return banner_start + std::string("\
   __ ___     ___ _  _        ___      ___  _\n\
  (_   |  |\\ | | |_ |_   |\\/|  |  |\\/|  |  |_)\n\
  __) _|_ | \\| | |_ |    |  | _|_ |  | _|_ | \\\n\
  ---------------> Fish School <--------------\n\033[0m")
    + banner_stop;
  }

  std::string LeadlineBanner()
  {
    return banner_start + std::string("\
   __ ___     ___ _  _        ___      ___  _\n\
  (_   |  |\\ | | |_ |_   |\\/|  |  |\\/|  |  |_)\n\
  __) _|_ | \\| | |_ |    |  | _|_ |  | _|_ | \\\n\
  ---------------> Leadline <---------------")
    + banner_stop;
  }

  std::string PlannerBanner()
  {
    return banner_start + std::string("\
   __ ___     ___ _  _        ___      ___  _\n\
  (_   |  |\\ | | |_ |_   |\\/|  |  |\\/|  |  |_)\n\
  __) _|_ | \\| | |_ |    |  | _|_ |  | _|_ | \\\n\
  --------------> Purse Planner <-------------")
    + banner_stop;
  }


}
