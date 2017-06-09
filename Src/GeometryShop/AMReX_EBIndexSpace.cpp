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
#include "AMReX_Utility.H"
#include <string>
#include "AMReX_Utility.H"
#include "AMReX_FabArrayIO.H"
#include <iostream>
#include <sstream>


namespace amrex
{
  ///
  void 
  EBIndexSpace::
  write(const string& a_dirname) const
  {
    //this creates the directory of all the stuff
    UtilCreateDirectoryDestructive(a_dirname, true);
    writeHeader(a_dirname);
    for(int ilev = 0; ilev < m_nlevels; ilev++)
    {
      string levdirname = a_dirname + "/_lev_" + EBArith::convertInt(ilev);
      UtilCreateDirectoryDestructive(levdirname, true);
      m_ebisLevel[ilev]->write(levdirname);
    }
  }


  ///
  void 
  EBIndexSpace::
  read(const string& a_dirname)
  {
    readHeader(a_dirname);
    m_ebisLevel.resize(m_nlevels);
    for(int ilev = 0; ilev < m_ebisLevel.size(); ilev++)
    {
      delete m_ebisLevel[ilev];
    }
    m_ebisLevel.resize(m_nlevels, NULL);

    for(int ilev = 0; ilev < m_nlevels; ilev++)
    {
      m_ebisLevel[ilev] = new EBISLevel();
      string levdirname = a_dirname + "/_lev_" + EBArith::convertInt(ilev);
      m_ebisLevel[ilev]->read(levdirname);
    }
    m_isDefined = true;
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

    amrex::Print() << "EBIndexSpace::define - From domain" << "\n";

    amrex::Print() << "  Building finest level..." << "\n";

    //this computes how many levels
    buildFirstLevel(a_domain, a_origin, a_dx, a_geoserver, a_nCellMax, a_maxCoarsenings);

    //starting at one because buildFirstLevel built 0
    for(int ilev = 1; ilev < m_nlevels; ilev++)
    {
      amrex::Print() << "  Building level " << ilev << "..." << "\n";
      buildNextLevel(a_geoserver, ilev);
    }
  }
  ///
  void 
  EBIndexSpace::
  buildFirstLevel(const Box&   a_domain,
                  const RealVect&        a_origin,
                  const Real&            a_dx,
                  const GeometryService& a_geoserver,
                  int a_nCellMax,
                  int a_maxCoarsenings)
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
      
    assert(!a_domain.isEmpty());
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
      m_nlevels =  std::max(m_nlevels, a_maxCoarsenings+1);
    }
      
    m_ebisLevel.resize(m_nlevels, NULL);
    m_domainLevel.resize(m_nlevels);
      
    Box  domLevel = a_domain;
    m_ebisLevel[0] = new EBISLevel(domLevel,
                                   a_origin,
                                   a_dx,
                                   m_nCellMax,
                                   a_geoserver);
      
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
    m_ebisLevel[ilev] = new EBISLevel(*m_ebisLevel[ilev-1],
                                      a_geoserver);
      
  }
  ///
  void 
  EBIndexSpace::
  clear()
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
  ///
  int EBIndexSpace::getLevel(const Box& a_domain) const
  {
    bool found = false;
    int whichlev = -1;
    for (int ilev = 0; ilev < m_domainLevel.size() && !found; ilev++)
    {
      if (m_domainLevel[ilev] == a_domain)
      {
        found = true;
        whichlev = ilev;
      }
    }
    return whichlev;
  }
      
  void EBIndexSpace::fillEBISLayout(EBISLayout     & a_ebisLayout,
                                    const BoxArray & a_grids,
                                    const DistributionMapping & a_dm,
                                    const Box      & a_domain,
                                    const int      & a_nghost) const
  {
    assert(isDefined());
    BL_PROFILE("EBIndexSpace::fillEBISLayout");
      
    //figure out which level we are on
    int whichlev = getLevel(a_domain);
    if (whichlev < 0)
    {
      amrex::Print() << "a_domain = " << a_domain
                     << " does not correspond to any refinement of EBIS" << "\n";
      amrex::Error("Bad argument to EBIndexSpace::fillEBISLayout");
    }
    m_ebisLevel[whichlev]->fillEBISLayout(a_ebisLayout, a_grids, a_dm, a_nghost);
  }
}
      
