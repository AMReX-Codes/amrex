#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "PlaneIF.H"

#include "NamespaceHeader.H"

PlaneIF::PlaneIF(const RealVect& a_normal,
                 const RealVect& a_point,
                 const bool&     a_inside)
  :HyperPlaneIF(IndexTM<Real, SpaceDim>(D_DECL(a_normal[0], a_normal[1], a_normal[2])),
                IndexTM<Real, SpaceDim>(D_DECL(a_point[0],  a_point[1],  a_point[2])),
                a_inside)
{
}

#include "NamespaceFooter.H"
