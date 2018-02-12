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
#include "AMReX_Print.H"
#include "AMReX_EBIndexSpace.H"
#include "AMReX_EBISLayout.H"
#include "AMReX_VoFIterator.H"
#include "AMReX_IrregNode.H"
#include "AMReX_PolyGeom.H"
#include "AMReX_parstream.H"
#include "AMReX_Utility.H"
#include "AMReX_Utility.H"
#include "AMReX_FabArrayIO.H"
#include <string>
#include <iostream>
#include <sstream>
#include <algorithm>


namespace amrex
{
  void 
  EBIndexSpace::
  getFinestLevelWithMultivaluedCells(Box& a_domain, int& a_levelNumber) const
  {
    bool found = false;
    a_levelNumber = -1;
    a_domain = Box();
//begin debug
//    for(int ilev = 0; ((ilev < m_nlevels)); ilev++)
//    {
//      bool fnerg = m_ebisLevel[ilev]->hasMultiValuedCells();      
//      if(fnerg)
//      {
//        amrex::Print() << "level " << ilev << "has multi cells" << endl;
//      }
//      else
//      {
//        amrex::Print() << "level " << ilev << "has no multi cells" << endl;
//      }
//    }
//end debug

    for(int ilev = 0; ((ilev < m_nlevels) && !found); ilev++)
    {

      if((m_ebisLevel[ilev]->hasMultiValuedCells()) && !found)
      {
        a_levelNumber = ilev;
        a_domain = m_domainLevel[ilev];
        found = true;
      }
    }
  }

  ///
  void 
  EBIndexSpace::
  write(const string& a_dirname) const
  {
    //this creates the directory of all the stuff
    UtilCreateCleanDirectory(a_dirname, true);
    writeHeader(a_dirname);
    for(int ilev = 0; ilev < m_nlevels; ilev++)
    {
      string levdirname = a_dirname + "/_lev_" + EBArith::convertInt(ilev);
//      UtilCreateCleanDirectory(levdirname, true); done inside the function
      m_ebisLevel[ilev]->write(levdirname);
    }
  }


  ///
  void 
  EBIndexSpace::
  read(const string& a_dirname)
  {
    amrex::Print() << "EBIndexSpace::define - From file input" << endl;

    readHeader(a_dirname);
    m_ebisLevel.clear();
    m_ebisLevel.resize(m_nlevels);
    for(int ilev = 0; ilev < m_nlevels; ilev++)
    {
      m_ebisLevel[ilev].reset(new EBISLevel());
      amrex::Print() << "EBIndexSpace::reading level " << ilev << endl;
      string levdirname = a_dirname + "/_lev_" + EBArith::convertInt(ilev);
      m_ebisLevel[ilev]->read(levdirname);
    }
    m_isDefined = true;
    amrex::Print() << "leaving EBIS::read" << endl;
  }

  ///
  void 
  EBIndexSpace::
  writeHeader(const string& a_dirname) const
  {
    std::ofstream headerfile;
    string filename = a_dirname + string("/headerfile");
    headerfile.open(filename.c_str(), std::ios::out | std::ios::trunc);
    headerfile << m_nlevels  << endl;
    headerfile << m_nCellMax << endl;
    for(int ilev = 0; ilev < m_domainLevel.size(); ilev++)
    {
      headerfile << m_domainLevel[ilev] << endl;
    }

    headerfile.flush();
    headerfile.close();
  }

  ///
  void 
  EBIndexSpace::
  readHeader(const string& a_dirname)
  {
    std::ifstream headerfile;
    string filename = a_dirname + string("/headerfile");
    headerfile.open(filename.c_str(), std::ios::in);
    headerfile >> m_nlevels ;
    headerfile >> m_nCellMax;
    m_domainLevel.resize(m_nlevels);
    for(int ilev = 0; ilev < m_domainLevel.size(); ilev++)
    {
      headerfile >> m_domainLevel[ilev];
    }
    headerfile.close();
    
  }
  void 
  EBIndexSpace::
  fillNodeFarrayBoxFromImplicitFunction(FArrayBox& a_fab, const RealVect& a_dx, RealVect a_origin ) const
  {
    if(!m_implicitFunction)
    {
      amrex::Error("this EBIS was not defined with an implicit function");
    }
    BL_ASSERT(a_fab.nComp() >= 1);
    for(BoxIterator bit(a_fab.box()); bit.ok(); ++bit)
    {
      RealVect loc;
      for(int idir = 0; idir < SpaceDim; idir++)
      {
        loc[idir] = a_origin[idir] + a_dx[idir]*bit()[idir];
      }
      for(int idir = 0; idir < SpaceDim; idir++)
      {
        a_fab(bit(), 0) = m_implicitFunction->value(loc);
      }
    }
  }
  ///
  void 
  EBIndexSpace::
  define(const Box              & a_domain,
         const RealVect         & a_origin,
         const Real             & a_dx,
         const GeometryService  & a_geoserver,
         int                      a_nCellMax,
         int                      a_maxCoarsenings)
  {
    BL_PROFILE("EBIndexSpace::define_geoserver_domain0");

    amrex::Print() << "EBIndexSpace::define - From domain" << a_domain << endl;

    amrex::Print() << "  Building finest level..." << endl;

    if(a_geoserver.hasImplicitFunction())
    {
      m_implicitFunction = a_geoserver.getImplicitFunction();
    }
    //this computes how many levels
    buildFirstLevel(a_domain, a_origin, a_dx, a_geoserver, a_nCellMax, a_maxCoarsenings);

    //starting at one because buildFirstLevel built 0
    for(int ilev = 1; ilev < m_nlevels; ilev++)
    {
      amrex::Print() << "  Building level " << ilev << "..." << endl;
      buildNextLevel(a_geoserver, ilev);
    }
    amrex::Print() << "leaving EBIS::define" << endl;
  }
  ///
  void 
  EBIndexSpace::
  buildFirstLevel(const Box&             a_domain,
                  const RealVect&        a_origin,
                  const Real&            a_dx,
                  const GeometryService& a_geoserver,
                  int                    a_nCellMax,
                  int                    a_maxCoarsenings)
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
    bool canref = (a_domain.coarsenable(2));
      
    AMREX_ASSERT(!a_domain.isEmpty());
    Box refbox = a_domain;
      
    while (canref)
    {
      if (!refbox.coarsenable(2))
      {
        canref = false;
      }
      else
      {
        m_nlevels++;
        refbox.coarsen(2);
      }
    }
    if (a_maxCoarsenings >= 0)
    {
      m_nlevels =  std::min(m_nlevels, a_maxCoarsenings+1);
    }
      
    m_ebisLevel.resize(m_nlevels);
    m_domainLevel.resize(m_nlevels);
      
    Box  domLevel = a_domain;
    m_ebisLevel[0].reset(new EBISLevel(domLevel,
                                       a_origin,
                                       a_dx,
                                       m_nCellMax,
                                       a_geoserver));
      
    m_domainLevel[0] = domLevel;
  }
  ///
  void
  EBIndexSpace::
  buildNextLevel(const GeometryService & a_geoserver,
                 const int             & a_whichlev)
  {
    BL_PROFILE("building_coarser_ebislevel");
    int ilev= a_whichlev;
      
    m_domainLevel[ilev] = m_domainLevel[ilev-1];
    m_domainLevel[ilev].coarsen(2);
    m_ebisLevel[ilev].reset(new EBISLevel(*m_ebisLevel[ilev-1],
                                          a_geoserver));
      
  }
  ///
  void 
  EBIndexSpace::
  clear()
  {
    m_ebisLevel.resize(0);
    m_domainLevel.resize(0);
    m_nlevels = 0;
    m_isDefined = false;
  }
  ///
  int EBIndexSpace::getLevel(const Box& a_domain) const
  {
      int whichlev = std::distance(std::begin(m_domainLevel),
                                   std::find(std::begin(m_domainLevel),
                                             std::end(m_domainLevel),
                                             a_domain));
      if (whichlev >= m_domainLevel.size()) whichlev = -1;
      return whichlev;
  }
      
  void EBIndexSpace::fillEBISLayout(EBISLayout     & a_ebisLayout,
                                    const BoxArray & a_grids,
                                    const DistributionMapping & a_dm,
                                    const Box      & a_domain,
                                    const int      & a_nghost) const
  {
    BL_PROFILE("EBIndexSpace::fillEBISLayout");
    AMREX_ASSERT(isDefined());
    AMREX_ASSERT(m_domainLevel.size() == m_ebisLevel.size());
      
    //figure out which level we are on
    int whichlev = getLevel(a_domain);
    if (whichlev < 0)
    {
      amrex::Print() << "a_domain = " << a_domain
                     << " does not correspond to any refinement of EBIS" << endl;
      amrex::Error("Bad argument to EBIndexSpace::fillEBISLayout");
    }
    m_ebisLevel[whichlev]->fillEBISLayout(a_ebisLayout, a_grids, a_dm, a_nghost);
  }
}
      
