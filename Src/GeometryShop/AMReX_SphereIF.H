#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _SPHEREIF_H_
#define _SPHEREIF_H_

#include "MayDay.H"
#include "RealVect.H"

#include "HyperSphereIF.H"

#include "NamespaceHeader.H"

///
/**
    This implicit function specifies a sphere.
 */
class SphereIF: public HyperSphereIF
{
public:
  ///
  /**
      Constructor specifying sphere radius (a_radius), center (a_center), and
      whether the domain is on the inside (a_inside).
   */
  SphereIF(const Real&     a_radius,
           const RealVect& a_center,
           const bool&     a_inside);

  /// Destructor
  virtual ~SphereIF()
  {;}



private:
  SphereIF();

  void operator=(const SphereIF& a_inputIF)
  {
    MayDay::Abort("SphereIF doesn't allow assignment");
  }
};

#include "NamespaceFooter.H"
#endif
