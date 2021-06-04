#pragma once

#include <cinttypes>

#include "casadi/casadi.hpp"
//#include <casadi/core/optistack.hpp>
#include "mimir/FkinDds.hpp"

/// Private namespace for mimir application.
namespace mimir
{
#ifndef DOXYGEN_SHOULD_SKIP_THIS
  class MimirPriv
  {
  public:
    MimirPriv();
    void solve();
  private:
    casadi::Opti m_opti;
  };
#endif
}
