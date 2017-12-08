/*
 *       {_       {__       {__{_______              {__      {__
 *      {_ __     {_ {__   {___{__    {__             {__   {__  
 *     {_  {__    {__ {__ { {__{__    {__     {__      {__ {__   
 *    {__   {__   {__  {__  {__{_ {__       {_   {__     {__     
 *   {______ {__  {__   {_  {__{__  {__    {_____ {__  {__ {__   
 *  {__       {__ {__       {__{__    {__  {_         {__   {__  
 * {__         {__{__       {__{__      {__  {____   {__      {__
 *
 */


#include "AMReX_IntVect.H"
#include "AMReX_GeometryService.H"

namespace amrex
{

  GeometryService::GeometryService()
  {
  }
  GeometryService::~GeometryService()
  {
  }

  bool GeometryService::isIrregular(const Box&           a_region,
                                    const Box&           a_domain,
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
                                                        const Box&           a_domain,
                                                        const RealVect&      a_origin,
                                                        const Real&          a_dx) const
  {
    if (isRegular(a_region, a_domain, a_origin, a_dx)) return GeometryService::Regular;
    if (isCovered(a_region, a_domain, a_origin, a_dx)) return GeometryService::Covered;
    return GeometryService::Irregular;
  }

  bool GeometryService::pointOutside(const RealVect& a_pt) const
  {
    amrex::Abort("GeometryService::pointOutside not implemented in the general case");
    return true;
  }

  bool GeometryService::intersection(const RealVect& a_lo1, const RealVect& a_hi1,
                                     const RealVect& a_lo2, const RealVect& a_hi2)
  {
    if (a_lo1[0] >= a_hi2[0] || a_hi1[0] <= a_lo2[0]) return false;
#if BL_SPACEDIM > 1
    if (a_lo1[1] >= a_hi2[1] || a_hi1[1] <= a_lo2[1]) return false;
#endif
#if BL_SPACEDIM > 2
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
    hi += RealVect((a_region.bigEnd()+IntVect::TheUnitVector()))*a_dx;
    return intersection(a_lower, a_upper, lo, hi);

  }

}

