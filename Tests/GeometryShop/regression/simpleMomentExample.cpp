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
#include "AMReX_BaseIVFactory.H"
#include "AMReX_SphereIF.H"
#include "AMReX_RealVect.H"
#include "AMReX_WrappedGShop.H"

namespace amrex
{
  /***************/
  void
  getFinestDomain(Box&       a_domain,
                  RealVect&      a_dx)
  {
    ParmParse pp;
    std::vector<int> n_cell(SpaceDim);
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
    for (int idir=0;idir<SpaceDim;idir++)
    {
      a_dx[idir] = prob_hi/n_cell[idir];
    }
  }
  void
  makeGeometry(const Box&       a_domain,
               const RealVect&  a_dx,
               RealVect&        a_sphereCenter,
               Real&            a_sphereRadius)
  {
    //parse input file.  single level
    ParmParse pp;
    pp.get("sphere_radius", a_sphereRadius);
    EBIndexSpace* ebisPtr = AMReX_EBIS::instance();
    int biggridsize;
    pp.get("max_grid_size", biggridsize);
    vector<Real>  sphereCenterVect(SpaceDim);
    pp.getarr("sphere_center",sphereCenterVect, 0, SpaceDim);
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      a_sphereCenter[idir] = sphereCenterVect[idir];
    }

    amrex::Print() << "using a sphere implicit function" << "\n";
    bool negativeInside = true;
    amrex::Print() << "using WrappedGShop " << "\n";
      
    SphereIF lalaBall(a_sphereRadius, a_sphereCenter, negativeInside);
    WrappedGShop workshop(lalaBall);
    int ebmaxcoarsen = 0;
    RealVect origin = RealVect::Zero;
    ebisPtr->define(a_domain, origin, a_dx[0], workshop, biggridsize, ebmaxcoarsen);

  }
/************/
  void
  simpleMomentExample()
  {
    //make layouts == domain
    Box domainBox;
    RealVect dx;
    getFinestDomain(domainBox, dx);
    BoxArray dbl(domainBox);
    DistributionMapping dm(dbl);


    RealVect  sphereCenter;
    Real      sphereRadius;

    makeGeometry(domainBox,  dx, sphereCenter, sphereRadius);
    EBISLayout ebisl;
    const EBIndexSpace* const ebisPtr = AMReX_EBIS::instance();
    ebisPtr->fillEBISLayout(ebisl, dbl, dm, domainBox, 0);

    //just prints out the moments of the first irregular cell it finds.
    static bool printedStuff = false;
    for (MFIter mfi(dbl, dm); mfi.isValid(); ++mfi)
    {
      const EBISBox& ebisBox = ebisl[mfi];
      Box               grid =   dbl[mfi];
      IntVectSet ivs = ebisBox.getIrregIVS(grid);
      VoFIterator vofit(ivs, ebisBox.getEBGraph());
      Vector<VolIndex> vofs = vofit.getVector();
      if(vofs.size() > 0)
      {
        int ivof = 0;
        IndMomSpaceDim volmom = ebisBox.getEBData().getVolumeMoments(vofs[ivof]);
        amrex::Print() <<  "The volume moments for the vof "  << vofs[ivof] <<  " are as follows " << endl;
        for(MomItSpaceDim momit; momit.ok(); ++momit)
        {
          IvSpaceDim iv = momit();
          amrex::Print() << "integral_cell ";
          for(int idir = 0; idir < SpaceDim; idir++)
          {
            amrex::Print() << "x_" <<idir << "^" << iv[idir] << " ";
          }
          amrex::Print() << "dV = " << volmom[iv] << endl;
            
        }
        printedStuff = true;
      }
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

  amrex::simpleMomentExample();
  

  amrex::Finalize();
  return retval;
}
