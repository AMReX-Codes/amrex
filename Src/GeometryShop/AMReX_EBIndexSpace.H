
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


#ifndef _EBINDEXSPACE_H_
#define _EBINDEXSPACE_H_

#include "AMReX_REAL.H"
#include "AMReX_IntVect.H"
#include "AMReX_IntVectSet.H"
#include "AMReX_VolIndex.H"
#include "AMReX_FaceIndex.H"
#include "AMReX_EBISBox.H"
#include "AMReX_EBISLayout.H"
#include "AMReX_GeometryService.H"
#include "AMReX_EBISLevel.H"

namespace amrex
{


///
/**
   EBIndexSpace represents the geometric information
   of the domain.    It should be generated through
   the Chombo_EBIS class.  This follows the singleton
   pattern.
*/
  class EBIndexSpace
  {
  public:

    EBIndexSpace();

    ///
    /**
       If a_ncellMax is set, that is the max width of
       an internal grid.  Otherwise use defaults
       of (16 in 3D, 64 in 2d)
    */
    void
    define(const ProblemDomain   & a_domain,
           const RealVect        & a_origin,
           const Real            & a_dx,
           const GeometryService & a_geoserver,
           int                     a_nCellMax = -1,
           int                     a_maxCoarsenings = -1);


    void
    define(EBISLevel * a_level0,
           int         a_nCellMax = -1,
           int         a_maxCoarsenings = -1);

    void
    define(const ProblemDomain                        & a_entireDomain,
           const RealVect                             & a_origin,
           const Real                                 & a_dx,
           const Vector<RefCountedPtr<EBIndexSpace> > & a_patches,
           const Vector<IntVect>                      & a_offsets,
           int                                          a_maxCoarsenings = -1);

    EBISLevel* buildFirstLevel(const ProblemDomain& a_domain,
                               const RealVect& a_origin,
                               const Real& a_dx,
                               const GeometryService& a_geoserver,
                               int a_nCellMax,
                               int a_maxCoarsenings,
                               bool a_fixRegularNextToMultiValued = true);

    EBISLevel* buildNextLevel(const GeometryService & a_geoserver,
                              bool                    a_fixRegularNextToMultiValued = true);

    ///
    ~EBIndexSpace();


    ///
    /**
       return level index of domain.
       return -1 if a_domain does not
       correspond to any refinement of EBIS.  public for sanity checks
    */
    int getLevel(const ProblemDomain& a_domain) const;

    ///
    void fillEBISLayout(EBISLayout&              a_ebis,
                        const DisjointBoxLayout& a_grids,
                        const ProblemDomain&     a_domain,
                        const int&               a_nghost) const;

    ///
    /**
       Return true if the define function has been called.
    */
    bool isDefined() const;

    ///
    void clear();


  private:

    int  m_nCellMax;
    bool m_isDefined;

    std::vector<EBISLevel*> m_ebisLevel;
    std::vector<ProblemDomain> m_domainLevel;

    int m_nlevels;

    void operator=(const EBIndexSpace& ebiin)
      {
        MayDay::Error("no assignment of EBIndexSpace allowed");
      }

    EBIndexSpace(const EBIndexSpace& ebiin)
      {
        MayDay::Error("no copy construction of EBIndexSpace allowed");
      }
  };

}

#endif
