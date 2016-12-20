#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _PLANEIF_H_
#define _PLANEIF_H_

#include "MayDay.H"
#include "RealVect.H"

#include "HyperPlaneIF.H"

#include "NamespaceHeader.H"

///
/**
    This implicit function specifies a plane.
 */
class PlaneIF: public HyperPlaneIF
{
public:
  ///
  /**
      Constructor specifying plane normal (a_normal), a point on the plance
      (a_point), and whether the domain is the same side as the normal
      (a_inside).
   */
  PlaneIF(const RealVect& a_normal,
          const RealVect& a_point,
          const bool&     a_inside);

  /// Destructor
  virtual ~PlaneIF()
  {
  }

protected:

private:
  PlaneIF();

};

#include "NamespaceFooter.H"
#endif
