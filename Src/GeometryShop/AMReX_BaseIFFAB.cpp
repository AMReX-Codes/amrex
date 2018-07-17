#include <cmath>
#include <cstdlib>
#include "AMReX_BaseIFFAB.H"
#include "AMReX_REAL.H"

namespace amrex
{
  template< >
  bool BaseIFFAB<Real>::isValid(const Real& a_input) const
  {
    bool notvalid = std::isnan(a_input);
    return (!notvalid);
  }
}

