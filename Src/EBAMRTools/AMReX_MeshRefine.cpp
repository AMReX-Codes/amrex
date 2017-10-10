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
#include <AMReX_ParmParse.H>
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
                                      const amrex::Vector<amrex::IntVect>& bf_lev) override
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
    void setTags(const Vector<IntVectSet>* a_vec_ivs_ptr)
      {
        m_vec_ivs_ptr = a_vec_ivs_ptr; 
      }
      
  private:
    
    const Vector<IntVectSet>* m_vec_ivs_ptr;

  };

  void MeshRefine(Vector<BoxArray>           &   a_grids,
                  const Vector<IntVectSet>   &   a_tags,
                  const Vector<int>          &   a_refRat,
                  const int                       &   a_maxLev,
                  const int                       &   a_blockingFactor,
                  const int                       &   a_properNestingRadius,
                  const int                       &   a_maxGridSize,
                  const Real                      &   a_gridEfficiency,
                  const Box                       &   a_coarsestDomain)
  {
    RealBox prob_domain(AMREX_D_DECL(0.,0.,0.), AMREX_D_DECL(1.,1.,1.));

    Vector<int> n_cell(SpaceDim);
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
                  const GridParameters       & a_params)
  {
    BL_PROFILE("EBLevelDataOps::getAllIrregRefinedLayouts");
    a_eblg.resize(a_params.numLevels);


    //make the tags by refining everywhere and getting the irregular cells
    Vector<IntVectSet> tags(a_params.numLevels);
    Box domlev = a_params.coarsestDomain;
    for(int ilev = 0; ilev < a_params.numLevels; ilev++)
    {
      BoxArray bac(a_params.coarsestDomain);
      bac.maxSize(a_params.maxGridSize);
      DistributionMapping dmc(bac);
      EBLevelGrid eblglev(bac, dmc, domlev, a_params.ghostEBISBox);
      for(MFIter mfi(bac, dmc); mfi.isValid(); ++mfi)
      {
        tags[ilev] |= eblglev.getEBISL()[mfi].getIrregIVS(bac[mfi]);
      }
      //empty tags do bad things
      if(tags[ilev].isEmpty())
      {
        tags[ilev] |= IntVect::TheZeroVector();
      }

      if(ilev < (a_params.numLevels-1))
      {
        domlev.refine(a_params.refRatio[ilev]);
      }
    }

    Vector<BoxArray> grids;

    MeshRefine(grids, tags, a_params.refRatio, a_params.maxLevel, 
               a_params.blockFactor, a_params.bufferSize, a_params.maxGridSize, 
               a_params.fillRatio, a_params.coarsestDomain);

    domlev = a_params.coarsestDomain;
    for(int ilev = 0; ilev < a_params.numLevels; ilev++)
    {
      DistributionMapping dm_lev(grids[ilev]);
      a_eblg[ilev].define(grids[ilev], dm_lev, domlev, a_params.ghostEBISBox);
      if(ilev < (a_params.numLevels-1))
      {
        domlev.refine(a_params.refRatio[ilev]);
      }
    }
  

    Real totalPoints = 0;
    long long totalBoxes  = 0;
    for(int ilev = 0; ilev < a_params.numLevels; ilev++)
    {
      Real pointsThisLevel = 0;
      for(int ibox = 0; ibox < grids[ilev].size(); ibox++)
      {
        pointsThisLevel += grids[ilev][ibox].numPts();
      }
      totalPoints += pointsThisLevel;
      totalBoxes +=  grids[ilev].size();
      long long ipoints = pointsThisLevel;
      amrex::Print() << "getAllIrregRefineLayouts:level[" << ilev
                     << "], number of boxes = " << grids[ilev].size()
                     << ", number of points = " << ipoints << "\n";
    }
    long long ipoints = totalPoints;
    amrex::Print() << "getAllIrregRefineLayouts:"
                   <<  "   total boxes = " << totalBoxes
                   <<  ", total points = " << ipoints <<  "\n";
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
    amrex::Print() << "Input Parameters:                    "   << "\n";
    amrex::Print() <<"whichGeom           = " << whichGeom      << "\n";
    amrex::Print() <<"nCells              = " << nCells         << "\n";
    amrex::Print() <<"maxGridSize         = " << maxGridSize    << "\n";
    amrex::Print() <<"blockFactor         = " << blockFactor    << "\n";
    amrex::Print() <<"bufferSize          = " << bufferSize     << "\n";
    amrex::Print() <<"fillRatio           = " << fillRatio      << "\n";
    amrex::Print() <<"maxLevel            = " << maxLevel       << "\n";
    amrex::Print() <<"numLevels           = " << numLevels      << "\n";
    amrex::Print() <<"coarsestDomain      = " << coarsestDomain << "\n";
    amrex::Print() <<"coarsestDx          = " << coarsestDx     << "\n";
    amrex::Print() <<"domainLength        = " << domainLength   << "\n";
    amrex::Print() <<"ghostEBISBox        = " << ghostEBISBox   << "\n";
    amrex::Print() <<"refRatio            = ";
    for(int iref = 0; iref << refRatio.size(); iref++)
    {
      amrex::Print() <<  refRatio[iref]       << "  ";
    }
    amrex::Print() <<   "\n";

  }
  /********/
  void 
  getGridParameters(GridParameters&  a_params, 
                    int  a_forceMaxLevel,
                    bool a_verbose)
  {

    ParmParse pp;

    pp.get("which_geom"    , a_params.whichGeom            );

    Vector<int> nCellsArray(SpaceDim);
    pp.getarr("n_cell",nCellsArray,0,SpaceDim);

    for (int idir = 0; idir < SpaceDim; idir++)
    {
      a_params.nCells[idir] = nCellsArray[idir];
    }
    a_params.ghostEBISBox = 4;
    amrex::Print() << "number of ghost cells for EBISBox = " << a_params.ghostEBISBox <<"\n";
  

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

    pp.get("domain_length",a_params.domainLength);


    pp.get("which_geom",   a_params.whichGeom);
    pp.get("max_grid_size",a_params.maxGridSize);

    //derived stuff
    a_params.coarsestDx = a_params.domainLength/a_params.nCells[0];

    a_params.probLo = RealVect::Zero;
    a_params.probHi = RealVect::Zero;
    a_params.probHi += a_params.domainLength;
    if(a_verbose)
    {
      a_params.pout();
    }
  }
}

