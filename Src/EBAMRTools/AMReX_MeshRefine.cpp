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


#include <AMReX_MeshRefine.H>
#include <AMReX.H>
#include <AMReX_RealBox.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_AmrMesh.H>
#include <AMReX_BoxIterator.H>




/// Factory class to produce EBCellFABs.
/**
   This makes a set of grids from a set of tags.   
   This is not as general as Chombo's meshrefine--it is useful for 
   writing tests, however.
*/
namespace amrex
{
  class MRAmr: public AmrMesh
  {
  public:
    using AmrMesh::AmrMesh;  // inherit AmrMesh constructors

    virtual ~MRAmr() {}

    //! Manually tag.  Note that tags is built on level lev grids coarsened by bf_lev.
    virtual void ManualTagsPlacement (int lev, amrex::TagBoxArray& tags,
                                      const amrex::Array<amrex::IntVect>& bf_lev) override
      {
        IntVectSet coarsenedTagsLev = (*m_vec_ivs_ptr)[lev];
        coarsenedTagsLev.coarsen(bf_lev[lev][0]);

        for (MFIter mfi(tags); mfi.isValid(); ++mfi)
        {
          TagBox& tag = tags[mfi];
          const Box& bx = tag.box();
                                      
          for(BoxIterator boxit(bx); boxit.ok(); ++boxit)
          {
            if(coarsenedTagsLev.contains(boxit()))
            {
              tag(boxit()) = TagBox::SET;
            }
          }
        }
      }

    //not safe but this is a private class hidden away.  shhh.
    void setTags(const std::vector<IntVectSet>* a_vec_ivs_ptr)
      {
        m_vec_ivs_ptr = a_vec_ivs_ptr; 
      }
      
  private:
    
    const std::vector<IntVectSet>* m_vec_ivs_ptr;

  };

  void MeshRefine(std::vector<BoxArray>           &   a_grids,
                  const std::vector<IntVectSet>   &   a_tags,
                  const std::vector<int>          &   a_refRat,
                  const int                       &   a_maxLev,
                  const int                       &   a_blockingFactor,
                  const int                       &   a_properNestingRadius,
                  const int                       &   a_maxGridSize,
                  const Real                      &   a_gridEfficiency,
                  const Box                       &   a_coarsestDomain)
  {
    RealBox prob_domain(D_DECL(0.,0.,0.), D_DECL(1.,1.,1.));

    Array<int> n_cell(SpaceDim);
    IntVect domsize =  a_coarsestDomain.size();
    for(int idir = 0; idir < SpaceDim; idir++)
    {
      n_cell[idir] = domsize[idir];
    }

    MRAmr amr(&prob_domain, a_maxLev, n_cell, -1, a_refRat);

    amr.setTags(&a_tags);
    amr.SetMaxGridSize(a_maxGridSize);
    amr.SetBlockingFactor(a_blockingFactor);
    amr.SetGridEff(a_gridEfficiency);
    amr.SetNProper(a_properNestingRadius);
    amr.MakeNewGrids();

    a_grids.resize(amr.finestLevel() + 1);
    for (int lev = 0; lev <= amr.finestLevel(); ++lev)
    {
      a_grids[lev] = amr.boxArray(lev);
    }
    
  }
  //-----------------------------------------------------------------------
  void 
  getAllIrregEBLG(Vector<EBLevelGrid>   & a_eblg,
                  const GridParameters  & a_params)
  {
    BL_PROFILE("EBLevelDataOps::getAllIrregRefinedLayouts");
    a_eblg.resize(a_params.numLevels);


    //make the tags by refining everywhere and getting the irregular cells
    std::vector<IntVectSet> tags(a_params.numLevels);
    Box domlev = a_params.coarsestDomain;
    for(int ilev = 0; ilev < a_params.numLevels; ilev++)
    {
      BoxArray bac(a_params.coarsestDomain);
      bac.maxSize(a_params.maxGridSize);
      DistributionMapping dmc(ba);
      EBLevelGrid eblglev(bac, dmc, domlev, a_params.ghostEBISBox);
      if(ilev < (a_params.numLevels-1))
      {
        domlev.refine(a_params.refRatio[ilev]);
      }
      for(MFIter mfi(bac, dmc); mfi.isValid(); ++mfi)
      {
        tags[ilev] |= eblglev.getEBISL()[mfi].getIrregIVS(grid);
      }
      //empty tags do bad things
      if(tags[ilev].isEmpty())
      {
        tags[ilev] |= IntVect::TheZeroVector();
      }
    }

    std::vector<BoxArray> grids;

    MeshRefine(grids, tags, a_params.refRatio, a_params.maxLevel, 
               a_params.blockFactor, a_params.bufferSize, a_params.maxGridSize, 
               a_params.fillRatio, a_params.coarsestDomain);

    domlev = a_params.coarsestDomain;
    for(int ilev = 0; ilev < a_params.numLevels; ilev++)
    {
      DistributionMapping dm_lev(grids[ilev]);
      a_eblg[ilev].define(grids[ilev], dm_lev, domlev, a_params.ghostEBISBox);
    }
  

    Real totalPoints = 0;
    long long totalBoxes  = 0;
    for(int ilev = 0; ilev < a_params.numLevels; ilev++)
    {
      Real pointsThisLevel = 0;
      for(LayoutIterator lit = a_grids[ilev].layoutIterator(); lit.ok(); ++lit)
      {
        pointsThisLevel += a_grids[ilev][lit()].numPts();
      }
      totalPoints += pointsThisLevel;
      totalBoxes += a_grids[ilev].size();
      pout() << "getAllIrregRefineLayouts:level[" << ilev
             << "], number of boxes = " << a_grids[ilev].size()
             << ", number of points = " << pointsThisLevel << endl;
    }
    pout() << "getAllIrregRefineLayouts:"
           <<  "   total boxes = " << totalBoxes
           <<  ", total points = " << totalPoints <<  endl;
  }

  /********/
  void 
  GridParameters::coarsen(int a_factor)
  {
    coarsestDx *= a_factor;
    coarsestDomain.coarsen(a_factor);
  }

  /********/
  void 
  GridParameters::refine(int a_factor)
  {
    coarsestDx /= a_factor;
    coarsestDomain.refine(a_factor);
  }

  /********/
  void 
  GridParameters::pout() const
  {
    ::amrex::Print() << "Input Parameters:                    "   << endl;
    ::amrex::Print() <<"whichGeom           = " << whichGeom      << endl;
    ::amrex::Print() <<"nCells              = " << nCells         << endl;
    ::amrex::Print() <<"maxGridSize         = " << maxGridSize    << endl;
    ::amrex::Print() <<"blockFactor         = " << blockFactor    << endl;
    ::amrex::Print() <<"bufferSize          = " << bufferSize     << endl;
    ::amrex::Print() <<"fillRatio           = " << fillRatio      << endl;
    ::amrex::Print() <<"maxLevel            = " << maxLevel       << endl;
    ::amrex::Print() <<"numLevels           = " << numLevels      << endl;
    ::amrex::Print() <<"refRatio            = " << refRatio       << endl;
    ::amrex::Print() <<"coarsestDomain      = " << coarsestDomain << endl;
    ::amrex::Print() <<"coarsestDx          = " << coarsestDx     << endl;
    ::amrex::Print() <<"domainLength        = " << domainLength   << endl;
    ::amrex::Print() <<"ghostPhi            = " << ghostPhi       << endl;
    ::amrex::Print() <<"ghostRHS            = " << ghostRHS       << endl;
    ::amrex::Print() <<"ghostEBISBox        = " << ghostEBISBox   << endl;
  }
  /********/
  void 
  getGridParameters(GridParameters&  a_params, 
                    int  a_forceMaxLevel,
                    bool a_verbose)
  {

    ParmParse pp;

    pp.get("which_geom"    , a_params.whichGeom            );

    std::vector<int> nCellsArray(SpaceDim);
    pp.getarr("n_cells",nCellsArray,0,SpaceDim);

    for (int idir = 0; idir < SpaceDim; idir++)
    {
      a_params.nCells[idir] = nCellsArray[idir];
      a_params.ghostPhi[idir] = 6;
      a_params.ghostRHS[idir] = 6;
    }
    a_params.ghostEBISBox = 4;
    amrex::Print() << "ghost cells phi = " << a_params.ghostPhi << ", ghost cells rhs = "  << a_params.ghostRHS << endl;
    amrex::Print() << "number of ghost cells for EBISBox = " << a_params.ghostEBISBox << endl;
  

    if(a_forceMaxLevel == 0)
    {
      a_params.maxLevel  = 0;
      a_params.numLevels = 1;
      a_params.refRatio.resize(1,2);
      a_params.blockFactor = 8;
      a_params.fillRatio = 1;
      a_params.bufferSize = 123;
    }
    else
    {
      if(a_forceMaxLevel >= 0)
      {
        a_params.maxLevel = a_forceMaxLevel;
      }
      else
      {
        pp.get("max_level", a_params.maxLevel);
      }
      a_params.numLevels = a_params.maxLevel + 1;
      pp.getarr("ref_ratio",a_params.refRatio,0,a_params.numLevels);
      pp.get("block_factor",a_params.blockFactor);
      pp.get("fill_ratio",a_params.fillRatio);
      pp.get("buffer_size",a_params.bufferSize);
    }
    IntVect lo = IntVect::Zero;
    IntVect hi = a_params.nCells;
    hi -= IntVect::Unit;

    a_params.coarsestDomain = Box(lo, hi);

    std::vector<Real> dLArray(SpaceDim);
    pp.getarr("domain_length",dLArray,0,SpaceDim);
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      a_params.domainLength[idir] = dLArray[idir];
    }

    pp.get("which_geom",   a_params.whichGeom);
    pp.get("max_grid_size",a_params.maxGridSize);

    //derived stuff
    a_params.coarsestDx = a_params.domainLength[0]/a_params.nCells[0];

    a_params.probLo = RealVect::Zero;
    a_params.probHi = RealVect::Zero;
    a_params.probHi += a_params.domainLength;
    if(a_verbose)
    {
      a_params.pout();
    }
  }
}

