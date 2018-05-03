#ifndef _BASEBCFUNCEVAL_H___
#define _BASEBCFUNCEVAL_H___

#include "AMReX_RealVect.H"

namespace amrex
{
  ///
  class BaseBCFuncEval
  {
  public:
    ///
    BaseBCFuncEval()
      {
      }

    ///
    virtual ~BaseBCFuncEval()
      {
      }

    ///
    virtual Real value(const RealVect& a_point,
                       const int&      a_comp) const = 0;

    ///
    virtual Real derivative(const RealVect& a_point,
                            const int&      a_comp,
                            const int&      a_derivDir) const = 0;
  };
}

#endif
