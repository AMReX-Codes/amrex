#ifndef _BASEINDEX_H_
#define _BASEINDEX_H_

#include "AMReX_REAL.H"
#include "AMReX_IntVect.H"

namespace amrex
{
  /// index for other indicies to inherit
  /**
     Aggstencil needed the concept of a unified
     index.   This is it.  Not much to it.
  */
  class BaseIndex
  {
  public:

    BaseIndex()
    {
      ;
    }

    virtual ~BaseIndex()
    {
      ;
    }
  };
}
#endif
