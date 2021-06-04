#pragma once

#include <cinttypes>
#include <string>


/// This namespace holds the public API for mimir.
namespace mimir
{

// Test class, to be removed
#ifndef DOXYGEN_SHOULD_SKIP_THIS
  class MimirPriv;

  class MimirClass
  {
  public:
    MimirClass();
    void solve();
    ~MimirClass();
  private:
    int32_t m_count;
    mimir::MimirPriv* m_priv;
  };
#endif

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  /// Header printed on startup
  std::string PlannerBanner();
  /// Header printed on startup
  std::string KinematicVesselBanner();
  /// Header printed on startup
  std::string FishSchoolBanner();
  /// Header printed on startup
  std::string LeadlineBanner();
#endif

}
