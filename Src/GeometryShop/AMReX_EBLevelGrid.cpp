#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// dtgraves Nov 14, 2006

#include "REAL.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "DisjointBoxLayout.H"
#include "EBISLayout.H"
#include "EBCellFAB.H"
#include "Interval.H"
#include "EBIndexSpace.H"
#include "EBLevelGrid.H"
#include "EBArith.H"
#include "NamespaceHeader.H"

void
EBLevelGrid::
setDefaultValues()
{
  m_isDefined = false;
  m_isCoveringIVSDefined = false;
  m_nghost = -99;
  m_ebisPtr = NULL;
  m_cfivs = RefCountedPtr<LayoutData<IntVectSet> >(0);
}
/****/
void
EBLevelGrid::
defineCoveringIVS()
{
  CH_assert(m_isDefined);
  if (!m_isCoveringIVSDefined)
    {
      m_isCoveringIVSDefined = true;
      //define the covering ivs by first assembling its complement
      const Box& domBox = m_domain.domainBox();
      IntVectSet complementIVS(domBox);
      for (LayoutIterator lit = m_grids.layoutIterator(); lit.ok(); ++lit)
        {
          complementIVS -= m_grids.get(lit());
        }
      m_coveringIVS = IntVectSet(domBox);
      m_coveringIVS -= complementIVS;
    }
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
EBLevelGrid(const DisjointBoxLayout& a_dbl,
            const ProblemDomain&     a_domain,
            const int&               a_numGhostEBISL,
            const EBIndexSpace*      a_ebisPtr)
{
  setDefaultValues();
  define(a_dbl, a_domain, a_numGhostEBISL, a_ebisPtr);
}
/****/
EBLevelGrid::
EBLevelGrid(const DisjointBoxLayout& a_dbl,
            const EBISLayout&        a_ebisl,
            const ProblemDomain&     a_domain)
{
  setDefaultValues();
  define(a_dbl, a_ebisl, a_domain);
}
/****/
void
EBLevelGrid::
define(const DisjointBoxLayout& a_dbl,
       const EBISLayout&        a_ebisl,
       const ProblemDomain&     a_domain)
{
  m_isDefined = true;
  m_grids = a_dbl;
  m_domain = a_domain;
  m_ebisl = a_ebisl;
  m_ebisPtr = a_ebisl.getEBIS();
  m_nghost = m_ebisl.getGhost();

  m_cfivs = RefCountedPtr<LayoutData<IntVectSet> >(new LayoutData<IntVectSet>());
  EBArith::defineCFIVS(*m_cfivs, m_grids, m_domain);
}

/****/
void
EBLevelGrid::
define(const DisjointBoxLayout& a_dbl,
       const ProblemDomain&     a_domain,
       const int&               a_numGhostEBISL,
       const EBIndexSpace*      a_ebisPtr)
{
  m_isDefined = true;
  m_grids = a_dbl;
  m_domain = a_domain;
  m_ebisPtr = a_ebisPtr;
  m_nghost = a_numGhostEBISL;

  m_ebisPtr->fillEBISLayout(m_ebisl, a_dbl, m_domain, m_nghost);
  m_cfivs = RefCountedPtr<LayoutData<IntVectSet> >(new LayoutData<IntVectSet>());
  EBArith::defineCFIVS(*m_cfivs, m_grids, m_domain);
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
  CH_assert(a_eblgFine.coarsenable(a_ref));

  DisjointBoxLayout coarsenedFineGrids;
  coarsen(coarsenedFineGrids, a_eblgFine.getDBL(), a_ref);
  ProblemDomain domainCoar(a_eblgFine.getDomain());
  domainCoar.coarsen(a_ref);
  a_eblgCoar.define(coarsenedFineGrids,domainCoar,a_eblgFine.getGhost(),a_eblgFine.getEBIS());
}

void
refine(EBLevelGrid&       a_eblgFine,
       const EBLevelGrid& a_eblgCoar,
       const int&         a_ref)
{
  DisjointBoxLayout refinedCoarseGrids;
  refine(refinedCoarseGrids, a_eblgCoar.getDBL(), a_ref);
  ProblemDomain domainFine(a_eblgCoar.getDomain());
  domainFine.refine(a_ref);
  a_eblgFine.define(refinedCoarseGrids,domainFine,a_eblgCoar.getGhost(),a_eblgCoar.getEBIS());
}


#include "NamespaceFooter.H"
