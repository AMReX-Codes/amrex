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


#include <cmath>
#include <cstdio>
#include <iostream>

#include "AMReX.H"
#include "AMReX_Print.H"
#include "AMReX_ParmParse.H"

#include "AMReX_EBIndexSpace.H"
#include "AMReX_EBISLayout.H"
#include "AMReX_ParmParse.H"
#include "AMReX_GeometryShop.H"
#include "AMReX_PolyGeom.H"
#include "AMReX_VoFIterator.H"
#include "AMReX_BaseIVFAB.H"
#include "AMReX_MeshRefine.H"
#include "AMReX_SphereIF.H"
#include "AMReX_RealVect.H"
#include "AMReX_EBDebugDump.H"

namespace amrex
{
  /***************/
  void
  getBaseDomain(Box  &       a_domain,
                Real &       a_dx)
  {
    ParmParse pp;
    Array<int> n_cell(SpaceDim);
    pp.getarr("n_cell",n_cell,0,SpaceDim);

    BL_ASSERT(n_cell.size() == SpaceDim);
    IntVect lo = IntVect::TheZeroVector();
    IntVect hi;
    for (int ivec = 0; ivec < SpaceDim; ivec++)
    {
      if (n_cell[ivec] <= 0)
      {
        amrex::Print() << " bogus number of cells input = " << n_cell[ivec];
        amrex::Error();
      }
      hi[ivec] = n_cell[ivec] - 1;
    }

    a_domain = Box(lo, hi);

    Real prob_hi;
    pp.get("domain_length",prob_hi);
    a_dx = prob_hi/n_cell[0];

  }
  void
  makeGeometry(Array<int> &  a_refRat,
               int              &  a_maxGridSize,
               int              &  a_maxLevel,
               const Box        &  a_domainCoar,
               const Real       &  a_dxCoar)
  {
    //input gives the *coarsest* domain
    ParmParse pp;
    pp.get("max_level",a_maxLevel);
    a_refRat.resize(a_maxLevel+1);
    pp.getarr("ref_ratio", a_refRat, 0, a_maxLevel);
    Box domainFine = a_domainCoar;
    Real dxFine = a_dxCoar;
    for(int ilev = 1; ilev <= a_maxLevel; ilev++)
    {
      dxFine /= a_refRat[ilev-1];
      domainFine.refine(a_refRat[ilev-1]);
    }
    RealVect        sphereCenter;
    Real            sphereRadius;
    Array<Real>sphereCenterVect(SpaceDim);
    pp.get("sphere_radius", sphereRadius);
    pp.get("max_grid_size", a_maxGridSize);
    pp.getarr("sphere_center",sphereCenterVect, 0, SpaceDim);
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      sphereCenter[idir] = sphereCenterVect[idir];
    }

    amrex::Print() << "using a sphere implicit function" << "\n";
    bool negativeInside = true;
    pp.get("inside_regular", negativeInside);
    SphereIF lalaBall(sphereRadius, sphereCenter, negativeInside);
    GeometryShop workshop(lalaBall);
    int ebmaxcoarsen = -1;
    RealVect origin = RealVect::Zero;
    EBIndexSpace* ebisPtr = AMReX_EBIS::instance();
    ebisPtr->define(domainFine, origin, dxFine, workshop, a_maxGridSize, ebmaxcoarsen);
  }
  /************/
  void  getTags(Array<IntVectSet>& a_tags, 
                const int              & a_maxGridSize,
                const Array<int> & a_refRat, 
                const int              & a_maxLevel, 
                const Box              & a_domainCoar)
  {
    a_tags.resize(a_maxLevel+1); //might be able to get away with maxLevel
    const EBIndexSpace* const ebisPtr = AMReX_EBIS::instance();
    
    Box domLev = a_domainCoar;
    for(int ilev = 0; ilev <= a_maxLevel; ilev++)
    {
      BoxArray ba(domLev);
      ba.maxSize(a_maxGridSize);
      DistributionMapping dm(ba);
      EBISLayout ebisl;
      ebisPtr->fillEBISLayout(ebisl, ba, dm,  domLev, 0);
      shared_ptr<FabArray<EBGraph> > allgraphs = ebisl.getAllGraphs(); //just so I can get an iterator going
      for(MFIter mfi(*allgraphs); mfi.isValid(); ++mfi)
      {
        Box validbox = ba[mfi];
        EBISBox ebisBox = ebisl[mfi];
        IntVectSet ivsIrreg = ebisBox.getIrregIVS(validbox);
        a_tags[ilev] |= ivsIrreg;
      }
      domLev.refine(a_refRat[ilev]);
    }
  }
  /************/
  void
  getSomeGrids()
  {
    //make layouts == domain
    Box domainCoar;
    Real    dxCoar;

    getBaseDomain(domainCoar,  dxCoar);
    Array<int> refRat;
    int maxGridSize, maxLevel;
    makeGeometry(refRat, maxGridSize, maxLevel, domainCoar,  dxCoar);

    Array<IntVectSet> tags;
    
    getTags(tags, maxGridSize, refRat, maxLevel, domainCoar);

    Array<BoxArray> grids;
    int block, proper;
    Real gridEff;
    ParmParse pp;
    pp.get("blocking_factor", block);
    pp.get("n_proper", proper);
    pp.get("grid_eff", gridEff);
    
    MeshRefine(grids, tags, refRat, maxLevel, block, proper, maxGridSize, gridEff, domainCoar);

    amrex::Print() << "Inputs: \n" ;
    amrex::Print() <<  "coarsest domain        = "  << domainCoar  << "\n";
    amrex::Print() <<  "max grid size          = "  << maxGridSize << "\n";
    amrex::Print() <<  "max level number       = "  << maxLevel    << "\n";
    amrex::Print() <<  "grid efficiency        = "  << gridEff     << "\n";
    amrex::Print() <<  "proper nesting radius  = "  << proper      << "\n";
    amrex::Print() <<  "refinement ratio       = "  ;
    for(int ilev = 0; ilev < maxLevel; ilev++)
    {
      amrex::Print() << refRat[ilev] <<  "   "  ;
    }
    amrex::Print() << "\n"  ;
    amrex::Print() << "Outputs: \n" ;
    for(int ilev = 0; ilev <= maxLevel; ilev++)
    {
      //cout << grids[ilev].size() << endl;
      amrex::Print() << "grids for level " << ilev << " = "  << grids[ilev];
      amrex::Print() << "\n";
    }
  }
}
/************/
/************/
int
main(int argc, char* argv[])
{
  int retval = 0;
  amrex::Initialize(argc,argv);

  amrex::getSomeGrids();

  amrex::Print() << "mesh refine test passed" << endl;
  amrex::Finalize();
  return retval;
}
