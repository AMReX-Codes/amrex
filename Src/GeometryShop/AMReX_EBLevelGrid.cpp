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

#include "AMReX_EBLevelGrid.H"
#include "AMReX_parstream.H" 
namespace amrex
{
  /****/
  void
  EBLevelGrid::
  defineCFIVS(LayoutData<IntVectSet>    &  a_cfivs,
              const BoxArray            &  a_grids,
              const DistributionMapping &  a_dm,
              const Box                 &  a_probDom)
  {
    BL_PROFILE("EBLevelGrid::defineCFIVS()");

    a_cfivs.define(a_grids, a_dm);
    for (MFIter mfi(a_grids, a_dm); mfi.isValid(); ++mfi)
    {
      Box grownBox = mfi.validbox();
      grownBox.grow(1);
      grownBox &= a_probDom;
      const BoxList& bl = a_grids.complementIn(grownBox);
      for (const auto& b : bl) {
	  a_cfivs[mfi] |= b;
      }
    }
  }
  /****/
  void
  EBLevelGrid::
  setDefaultValues()
  {
    m_isDefined = false;
    m_nghost = -99;
  }
  /****/
  EBLevelGrid::
  EBLevelGrid()
  {
    setDefaultValues();
  }
  /****/
  EBLevelGrid::
  ~EBLevelGrid()
  {
  }
  /****/
  EBLevelGrid::
  EBLevelGrid(const BoxArray            & a_dbl,
              const DistributionMapping & a_dm,
              const Box                 & a_domain,
              const int                 & a_numGhostEBISL)
  {
    setDefaultValues();
    define(a_dbl, a_dm, a_domain, a_numGhostEBISL);
  }
  /****/
  EBLevelGrid::
  EBLevelGrid(const BoxArray            & a_dbl,
              const DistributionMapping & a_dm,
              const EBISLayout          & a_ebisl,
              const Box                 & a_domain)
  {
    setDefaultValues();
    define(a_dbl, a_dm, a_ebisl, a_domain);
  }
  /****/
  void
  EBLevelGrid::
  define(const BoxArray            & a_dbl,
         const DistributionMapping & a_dm,
         const EBISLayout          & a_ebisl,
         const Box                 & a_domain)
  {
    BL_PROFILE("EBLevelGrid::define()");
    m_isDefined = true;
    m_grids = a_dbl;
    m_dm    = a_dm;
    m_domain = a_domain;
    m_ebisl = a_ebisl;
    m_nghost  = m_ebisl.getGhost();
          
    m_cfivs = std::shared_ptr<LayoutData<IntVectSet> >(new LayoutData<IntVectSet>());
    defineCFIVS(*m_cfivs, m_grids, m_dm, m_domain);
  }
          
  /****/
  void
  EBLevelGrid::
  define(const BoxArray            & a_dbl,
         const DistributionMapping & a_dm,
         const Box                 & a_domain,
         const int                 & a_numGhostEBISL)
  {
    //pout() << "in eblevelgrid::define with nghost = " << a_numGhostEBISL << endl;
    m_isDefined = true;
    m_grids  = a_dbl;
    m_dm     = a_dm;
    m_domain = a_domain;
    m_nghost = a_numGhostEBISL;
    const EBIndexSpace* ebisPtr = AMReX_EBIS::instance();

    //pout() << "about to fillebislayout" << endl;

    ebisPtr->fillEBISLayout(m_ebisl, a_dbl, m_dm, m_domain, m_nghost);
    m_cfivs = std::shared_ptr<LayoutData<IntVectSet> >(new LayoutData<IntVectSet>());

    ////pout() << "about to define cfivs" << endl;

    defineCFIVS(*m_cfivs, m_grids, m_dm, m_domain);
    //pout() << "leaving eblevelgrid::define  " << endl;
  }
          
  EBLevelGrid::
  EBLevelGrid(const EBLevelGrid& a_eblg)
  {
    setDefaultValues();
    *this = a_eblg;
  }
          
  bool
  EBLevelGrid::
  coarsenable(const int& a_ref) const
  {
    return m_grids.coarsenable(a_ref);
  }
          
  void
  coarsen(EBLevelGrid&       a_eblgCoar,
          const EBLevelGrid& a_eblgFine,
          const int&         a_ref)
  {
    BL_ASSERT(a_eblgFine.coarsenable(a_ref));
          
    BoxArray coarsenedFineGrids = a_eblgFine.getDBL();
    Box              domainCoar = a_eblgFine.getDomain();
    coarsenedFineGrids.coarsen(a_ref);
    domainCoar.coarsen(a_ref);
    a_eblgCoar.define(coarsenedFineGrids, a_eblgFine.getDM(),domainCoar,a_eblgFine.getGhost());
  }
          
  void
  refine(EBLevelGrid&       a_eblgFine,
         const EBLevelGrid& a_eblgCoar,
         const int&         a_ref)
  {
    BoxArray refinedCoarseGrids = a_eblgCoar.getDBL();
    Box              domainFine = a_eblgCoar.getDomain();
    refinedCoarseGrids.refine(a_ref);
    domainFine.refine(a_ref);

    a_eblgFine.define(refinedCoarseGrids,a_eblgCoar.getDM(), domainFine,a_eblgCoar.getGhost());
  }
          
}
           
