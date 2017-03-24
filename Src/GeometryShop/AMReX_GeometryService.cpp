#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

//  ANAG, LBNL, DTG

#include "GeometryService.H"
#include "NamespaceHeader.H"

GeometryService::GeometryService()
{
}
GeometryService::~GeometryService()
{
}

bool GeometryService::isIrregular(const Box&           a_region,
                                  const ProblemDomain& a_domain,
                                  const RealVect&      a_origin,
                                  const Real&          a_dx) const
{
  return !(isRegular(a_region, a_domain, a_origin, a_dx) || isCovered(a_region, a_domain, a_origin, a_dx));
}

bool GeometryService::canGenerateMultiCells() const
{
  return true;
}

GeometryService::InOut GeometryService::InsideOutside(const Box&           a_region,
                                                      const ProblemDomain& a_domain,
                                                      const RealVect&      a_origin,
                                                      const Real&          a_dx) const
{
  if (isRegular(a_region, a_domain, a_origin, a_dx)) return GeometryService::Regular;
  if (isCovered(a_region, a_domain, a_origin, a_dx)) return GeometryService::Covered;
  return GeometryService::Irregular;
}

bool GeometryService::intersection(const RealVect& a_lo1, const RealVect& a_hi1,
                                   const RealVect& a_lo2, const RealVect& a_hi2)
{
  if (a_lo1[0] >= a_hi2[0] || a_hi1[0] <= a_lo2[0]) return false;
#if CH_SPACEDIM > 1
  if (a_lo1[1] >= a_hi2[1] || a_hi1[1] <= a_lo2[1]) return false;
#endif
#if CH_SPACEDIM > 2
  if (a_lo1[2] >= a_hi2[2] || a_hi1[2] <= a_lo2[2]) return false;
#endif
  return true;
}

bool GeometryService::intersection(const Box&           a_region,
                                   const RealVect&      a_origin,
                                   const Real&          a_dx,
                                   const RealVect&      a_lower,
                                   const RealVect&      a_upper)
{
  RealVect lo(a_origin), hi(a_origin);
  lo += RealVect(a_region.smallEnd())*a_dx;
  hi += RealVect((a_region.bigEnd()+IntVect::Unit))*a_dx;
  return intersection(a_lower, a_upper, lo, hi);

}

void GeometryService::postMakeBoxLayout(const DisjointBoxLayout& a_dbl,
                                        const RealVect& a_dx)
{
  return;
}



#include "NamespaceFooter.H"
