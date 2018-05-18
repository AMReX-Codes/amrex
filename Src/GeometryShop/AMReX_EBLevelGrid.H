#ifndef AMREX_EBLEVELGRID_H_
#define AMREX_EBLEVELGRID_H_

#include "AMReX_REAL.H"
#include "AMReX_FArrayBox.H"
#include "AMReX_LayoutData.H"
#include "AMReX_EBISLayout.H"
#include "AMReX_EBCellFAB.H"
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
    static void
    defineCFIVS(LayoutData<IntVectSet>    &  a_cfivs,
                const BoxArray            &  a_grids,
                const DistributionMapping &  a_dm,
                const Box                 &  a_probDom);
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
    EBLevelGrid(const BoxArray           & a_dbl,
                const DistributionMapping& a_dm,
                const Box                & a_domain,
                const int                & a_numGhostEBISL);

    ///
    EBLevelGrid(const BoxArray           & a_dbl,
                const DistributionMapping& a_dm,
                const EBISLayout         & a_ebisl,
                const Box                & a_domain);

    ///
    /**
       Define function that generates EBISL and CFIVS
    */
    void define(const BoxArray           & a_dbl,
                const DistributionMapping& a_dm,
                const Box                & a_domain,
                const int                & a_numGhostEBISL);


    ///
    /**
       only generates cfivs.  set EBISPtr = chombo_ebis::instance
    */
    void define(const BoxArray           & a_dbl,
                const DistributionMapping& a_dm,
                const EBISLayout         & a_ebisl,
                const Box                & a_domain);

    ///
    EBLevelGrid(const EBLevelGrid& ebcin);

    ///
    int getGhost() const
      {
        BL_ASSERT(m_isDefined);
        return m_nghost;
      }

    ///to keep porting of chombo codes less painful
    BoxArray getDBL() const
      {
        BL_ASSERT(m_isDefined);
        return getBoxArray();
      }

    ///
    BoxArray getBoxArray() const
      {
        BL_ASSERT(m_isDefined);
        return m_grids;
      }

    ///
    EBISLayout getEBISL() const
      {
        BL_ASSERT(m_isDefined);
        return m_ebisl;
      }

    ///
    std::shared_ptr<LayoutData<IntVectSet> > getCFIVS() const
      {
        BL_ASSERT(m_isDefined);
        return m_cfivs;
      }

    ///
    bool isDefined() const
      {
        return m_isDefined;
      }


    ///
    const Box& getDomain() const
      {
        BL_ASSERT(m_isDefined);
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
        BL_ASSERT(m_isDefined);
        m_ebisl.setMaxRefinementRatio(a_maxRefine);
      }

    ///
    /**
       Sets the maximum level of coarsening that the
       will have to perform.  Creates and holds new EBISLayouts
       at intermediate levels of refinement.
       Default is one (no coarsening done).
    */
    void setMaxCoarseningRatio(const int&                a_maxCoarsen)
      {
        BL_ASSERT(m_isDefined);
        m_ebisl.setMaxCoarseningRatio(a_maxCoarsen);
      }

    int getMaxCoarseningRatio() const
      {
        BL_ASSERT(m_isDefined);
        return m_ebisl.getMaxCoarseningRatio();
      }

    int
    getMaxRefinementRatio() const
      {
        BL_ASSERT(m_isDefined);
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


    friend void refine(EBLevelGrid&       a_eblgFine,
                       const EBLevelGrid& a_eblgCoar,
                       const int&         a_ref);

    DistributionMapping getDM() const
      {
        return m_dm;
      }
  protected:
    void setDefaultValues();
    bool m_isDefined;

    BoxArray                                 m_grids;
    DistributionMapping                      m_dm;
    EBISLayout                               m_ebisl;
    Box                                      m_domain;
    std::shared_ptr<LayoutData<IntVectSet> > m_cfivs;
    int                                      m_nghost;

  private:
  };

}

#endif
