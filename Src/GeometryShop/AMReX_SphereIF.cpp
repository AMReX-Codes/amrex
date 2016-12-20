#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "SphereIF.H"

#include "NamespaceHeader.H"

SphereIF::SphereIF(const Real&     a_radius,
                   const RealVect& a_center,
                   const bool&     a_inside)
  : HyperSphereIF(a_radius, IndexTM<Real, SpaceDim>(D_DECL(a_center[0], a_center[1], a_center[2])), a_inside)
{
}

#include "NamespaceFooter.H"
