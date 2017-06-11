#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _BASEBCFUNCEVAL_H_
#define _BASEBCFUNCEVAL_H_

#include "RealVect.H"
#include "EBArith.H"
#include "NamespaceHeader.H"

class BaseBCFuncEval
{
public:
  BaseBCFuncEval()
  {
  }

  virtual ~BaseBCFuncEval()
  {
  }

  virtual Real value(const RealVect& a_point,
                     const int&      a_comp) const = 0;

  virtual Real derivative(const RealVect& a_point,
                          const int&      a_comp,
                          const int&      a_derivDir) const = 0;
};

#include "NamespaceFooter.H"
#endif
