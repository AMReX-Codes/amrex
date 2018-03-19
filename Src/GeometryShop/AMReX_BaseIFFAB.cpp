#include <cmath>
#include <cstdlib>
#include "AMReX_BaseIFFAB.H"
#include "AMReX_REAL.H"

namespace amrex
{
  template< >
  bool BaseIFFAB<Real>::isValid(const int& a_iface,const int& varlocin) const
  {
    bool notvalid = std::isnan((*this)(a_iface, varlocin));
    return (!notvalid);
  }
}

