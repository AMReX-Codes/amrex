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


#include "AMReX_BoxIterator.H"
#include "AMReX_EBCFCopy.H"
#include "AMReX_EBIndexSpace.H"
#include "AMReX_EBGraphFactory.H"
#include "AMReX_EBDataFactory.H"
#include "AMReX_EBISLayout.H"
#include "AMReX_VoFIterator.H"
#include "AMReX_IrregNode.H"
#include "AMReX_PolyGeom.H"


namespace amrex
{


  bool EBIndexSpace::isDefined() const
  {
    return m_isDefined;
  }


  void EBIndexSpace::define(EBISLevel * a_level0,
                            int         a_nCellMax,
                            int         a_maxCoarsenings)
  {
    BL_PROFILE("EBIndexSpace::define_ebislevel0");

    std::cout << "EBIndexSpace::define - Given level 0" << endl;

    m_nCellMax = a_nCellMax;
    m_isDefined = true;

    const ProblemDomain& fineDomain = a_level0->getDomain();

    //figure out how deep we can go
    m_nlevels = 1;
    bool canref = (fineDomain == refine(coarsen(fineDomain,2), 2));
    assert(!fineDomain.isEmpty());
    ProblemDomain refbox = fineDomain;
    while (canref)
    {
      ProblemDomain refcoarbox = coarsen(refbox,2);
      refcoarbox.refine(2);
      if (refcoarbox != refbox)
      {
        canref = false;
      }
      else
      {
        m_nlevels++;
        refbox.coarsen(2);
      }
    }
    if (a_maxCoarsenings != -1)
    {
      assert(a_maxCoarsenings >= 0);
      if (m_nlevels > a_maxCoarsenings+1) m_nlevels = a_maxCoarsenings + 1;
    }

    m_ebisLevel.resize(m_nlevels, NULL);
    m_domainLevel.resize(m_nlevels);

    ProblemDomain  domLevel = fineDomain;
    a_level0->clearMultiBoundaries();
    m_ebisLevel[0] = a_level0;
    m_domainLevel[0] = domLevel;
    AllRegularService dummy;
    for (int ilev = 1; ilev < m_nlevels; ilev++)
    {
      std::cout << "  Generating level " << ilev << endl;
      domLevel.coarsen(2);
      m_domainLevel[ilev] = domLevel;
      m_ebisLevel[ilev] = new EBISLevel(*m_ebisLevel[ilev-1], dummy, this);
      m_ebisLevel[ilev]->clearMultiBoundaries();
    }
  }

  void EBIndexSpace::define(const ProblemDomain    & a_domain,
                            const RealVect         & a_origin,
                            const Real             & a_dx,
                            const GeometryService  & a_geoserver,
                            int                      a_nCellMax,
                            int                      a_maxCoarsenings)
  {
    BL_PROFILE("EBIndexSpace::define_geoserver_domain0");

    std::cout << "EBIndexSpace::define - From domain" << endl;

    std::cout << "  Building finest level..." << endl;

    buildFirstLevel(a_domain, a_origin, a_dx, a_geoserver, a_nCellMax, a_maxCoarsenings);
    m_ebisLevel[0]->clearMultiBoundaries();

    print_memory_line("ebis_after_first_level_build");
    EBISLevel* n =  m_ebisLevel[0];
    n->printGraphSummary("    ");
    std::cout << endl;
    int level = 1;
    while (n)
    {
      std::cout << "  Building level " << level << "..." << endl;
      n->clearMultiBoundaries();
      n=buildNextLevel(a_geoserver);

      if (n)
      {
        n->printGraphSummary("    ");
      }
      else
      {
        std::cout << "    Empty" << endl;
      }
      std::cout << endl;

      level++;
    }

  }

  EBISLevel* EBIndexSpace::buildFirstLevel(const ProblemDomain&   a_domain,
                                           const RealVect&        a_origin,
                                           const Real&            a_dx,
                                           const GeometryService& a_geoserver,
                                           int a_nCellMax,
                                           int a_maxCoarsenings,
                                           bool a_fixRegularNextToMultiValued)
  {
    BL_PROFILE("EBIndexSpace::buildFirstLevel");
    clear();
    m_isDefined = true;

    if (a_nCellMax > 0)
    {
      m_nCellMax = a_nCellMax;
    }
    else
    {
      if (SpaceDim == 2)
      {
        m_nCellMax = 64;
      }
      else
      {
        m_nCellMax = 16;
      }
    }
    m_nlevels = 1;
    bool canref = (a_domain == refine(coarsen(a_domain,2), 2));

    assert(!a_domain.isEmpty());
    ProblemDomain refbox = a_domain;

    while (canref)
    {
      ProblemDomain refcoarbox = coarsen(refbox,2);
      refcoarbox.refine(2);
      if (refcoarbox != refbox)
      {
        canref = false;
      }
      else
      {
        m_nlevels++;
        refbox.coarsen(2);
      }
    }
    if (a_maxCoarsenings != -1)
    {
      assert(a_maxCoarsenings >= 0);
      if (m_nlevels > a_maxCoarsenings+1) m_nlevels = a_maxCoarsenings + 1;
    }

    m_ebisLevel.resize(m_nlevels, NULL);
    m_domainLevel.resize(m_nlevels);

    ProblemDomain  domLevel = a_domain;
    m_ebisLevel[0] = new EBISLevel(domLevel,
                                   a_origin,
                                   a_dx,
                                   a_geoserver,
                                   this,
                                   a_fixRegularNextToMultiValued);
    m_domainLevel[0] = domLevel;
    return m_ebisLevel[0];
  }

  EBISLevel* EBIndexSpace::buildNextLevel(const GeometryService & a_geoserver,
                                          bool                    a_fixRegularNextToMultiValued)
  {
    int ilev=0;
    for ( ; ilev <m_ebisLevel.size(); ++ilev)
    {
      if (m_ebisLevel[ilev] == NULL) break;
    }
    if (ilev == m_ebisLevel.size()) return NULL;

    {
      string* leak = new string("EBIndexSpace::buildNextLevel_EBISLevel_");
      char levelString[100];
      sprintf(levelString,"%d",ilev);

      *leak += levelString;

      BL_PROFILE(leak->c_str());

      m_domainLevel[ilev] = m_domainLevel[ilev-1];
      m_domainLevel[ilev].coarsen(2);
      m_ebisLevel[ilev] = new EBISLevel(*m_ebisLevel[ilev-1],
                                        a_geoserver,
                                        this,
                                        a_fixRegularNextToMultiValued);

      delete leak;
    }

    return m_ebisLevel[ilev];
  }

  void EBIndexSpace::clear()
  {
    for (int ilev = 0; ilev < m_ebisLevel.size(); ilev++)
    {
      delete m_ebisLevel[ilev];
      m_ebisLevel[ilev] = NULL;
    }
    m_ebisLevel.resize(0);
    m_domainLevel.resize(0);
    m_nlevels = 0;
    m_isDefined = false;
  }

  EBIndexSpace::EBIndexSpace()
  {
  }

  EBIndexSpace::~EBIndexSpace()
  {
    clear();
  }

  int EBIndexSpace::getLevel(const ProblemDomain& a_domain) const
  {
    bool found = false;
    int whichlev = -1;
    for (int ilev = 0; ilev < m_domainLevel.size() && !found; ilev++)
    {
      if (m_domainLevel[ilev].domainBox() == a_domain.domainBox())
      {
        found = true;
        whichlev = ilev;
      }
    }
    return whichlev;
  }

  void EBIndexSpace::fillEBISLayout(EBISLayout&              a_ebisLayout,
                                    const DisjointBoxLayout& a_grids,
                                    const ProblemDomain&     a_domain,
                                    const int&               a_nghost) const
  {
    assert(isDefined());
    BL_PROFILE("EBIndexSpace::fillEBISLayout");

    //figure out which level we are on
    int whichlev = getLevel(a_domain);
    if (whichlev < 0)
    {
      std::cout << "a_domain = " << a_domain
             << " does not correspond to any refinement of EBIS" << endl;
      amrex::Error("Bad argument to EBIndexSpace::fillEBISLayout");
    }
    m_ebisLevel[whichlev]->fillEBISLayout(a_ebisLayout, a_grids, a_nghost);
  }
}

