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

#ifndef _EBLEVELGRID_H_
#define _EBLEVELGRID_H_

#include "AMReX_REAL.H"
#include "AMReX_FArrayBox.H"
#include "AMReX_LevelData.H"
#include "AMReX_EBISLayout.H"
#include "AMReX_EBCellFAB.H"
#include "AMReX_Interval.H"
#include "AMReX_EBIndexSpace.H"


namespace amrex
{
///
/**
   A structure to hold the various components of an EB calculation grid.
*/
  class EBLevelGrid
  {
  public:
    ///
    /**
       Default constructor.  User must subsequently call define().
    */
    EBLevelGrid();

    ///
    ~EBLevelGrid();

    ///
    /**
       Defining constructor.  Constructs a valid object.
       Equivalent to default construction followed by define().
    */
    EBLevelGrid(const BoxArray     & a_dbl,
                const Box          & a_domain,
                const int          & a_numGhostEBISL,
                const EBIndexSpace * a_ebisPtr);

    ///
    EBLevelGrid(const BoxArray    & a_dbl,
                const EBISLayout  & a_ebisl,
                const Box         & a_domain);

    ///
    /**
       Define function that generates EBISL and CFIVS
    */
    void define(const BoxArray     & a_dbl,
                const Box          & a_domain,
                const int          & a_numGhostEBISL,
                const EBIndexSpace * a_ebisPtr);


    ///
    /**
       only generates cfivs.  set EBISPtr = chombo_ebis::instance
    */
    void define(const BoxArray   & a_dbl,
                const EBISLayout & a_ebisl,
                const Box        & a_domain);

    ///
    EBLevelGrid(const EBLevelGrid& ebcin);

    ///
    int getGhost() const
      {
        CH_assert(m_isDefined);
        return m_nghost;
      }

    ///
    BoxArray getDBL() const
      {
        CH_assert(m_isDefined);
        return m_grids;
      }

    ///
    EBISLayout getEBISL() const
      {
        CH_assert(m_isDefined);
        return m_ebisl;
      }

    ///
    RefCountedPtr<LayoutData<IntVectSet> > getCFIVS() const
      {
        CH_assert(m_isDefined);
        return m_cfivs;
      }

    ///
    bool isDefined() const
      {
        return m_isDefined;
      }

    ///
    const EBIndexSpace* getEBIS() const
      {
        CH_assert(m_isDefined);
        return m_ebisPtr;
      }

    ///
    const IntVectSet& getCoveringIVS() const
      {
        CH_assert(m_isDefined);
        if (!m_isCoveringIVSDefined)
        {
          EBLevelGrid& getAroundConst = (EBLevelGrid&)(*this);
          getAroundConst.defineCoveringIVS();
        }
        return m_coveringIVS;
      }

    ///
    const Box& getDomain() const
      {
        CH_assert(m_isDefined);
        return m_domain;
      }

    ///
    /**
       Sets the maximum level of refinement that the EBISLayout
       will have to perform.  Creates and holds new EBISLayouts
       at intermediate levels of refinement.
       Default is one (no refinement done).
    */
    void setMaxRefinementRatio(const int& a_maxRefine)
      {
        CH_assert(m_isDefined);
        m_ebisl.setMaxRefinementRatio(a_maxRefine,m_ebisPtr);
      }

    ///
    /**
       Sets the maximum level of coarsening that the
       will have to perform.  Creates and holds new EBISLayouts
       at intermediate levels of refinement.
       Default is one (no coarsening done).
    */
    void setMaxCoarseningRatio(const int&                a_maxCoarsen,
                               const EBIndexSpace* const a_ebisPtr)
      {
        CH_assert(m_isDefined);
        m_ebisl.setMaxCoarseningRatio(a_maxCoarsen,a_ebisPtr);
      }

    int getMaxCoarseningRatio() const
      {
        CH_assert(m_isDefined);
        return m_ebisl.getMaxCoarseningRatio();
      }

    int
    getMaxRefinementRatio() const
      {
        CH_assert(m_isDefined);
        return m_ebisl.getMaxRefinementRatio();
      }
    ///
    /*
      returns m_grids.coarsenable(a_ref);
    */
    bool
    coarsenable(const int& a_ref) const;

    ///
    /**
       Coarsen eblgFine by a_ref
    */
    friend void coarsen(EBLevelGrid&       a_eblgCoar,
                        const EBLevelGrid& a_eblgFine,
                        const int&         a_ref);


    friend void refine(EBLevelGrid&       a_eblgCoar,
                       const EBLevelGrid& a_eblgFine,
                       const int&         a_ref);

  protected:
    void defineCoveringIVS();
    void setDefaultValues();
    bool m_isDefined, m_isCoveringIVSDefined;

    BoxArray                      m_grids;
    EBISLayout                             m_ebisl;
    Box                          m_domain;
    RefCountedPtr<LayoutData<IntVectSet> > m_cfivs;
    IntVectSet                             m_coveringIVS;
    const EBIndexSpace*                    m_ebisPtr;
    int                                    m_nghost;

  private:
  };

}

#endif
