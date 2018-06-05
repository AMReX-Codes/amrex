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
#include "AMReX_EBCellFAB.H"
#include "AMReX_EBArith.H"
#include "AMReX_EBCellFactory.H"
#include "AMReX_EBLevelDataOps.H"
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
  sumFineValues(Vector<Real>           & a_fineSum, 
                const EBCellFAB        & a_solutFine,
                const Vector<VolIndex> & a_fineVoFs)
  {
    BL_ASSERT(a_fineSum.size() == a_solutFine.nComp());
    for(int ivar = 0; ivar < a_solutFine.nComp(); ivar++)
    {
      Real value = 0;
      for(int ivof = 0; ivof < a_fineVoFs.size(); ivof++)
      {
        value += a_solutFine(a_fineVoFs[ivof], ivar);
      }
      a_fineSum[ivar] = value;
    }
  }
/************/
  void
  fillMomentNames(Vector<string> & a_momentNames)
  {
    int ncomp = IndMomSpaceDim::size();
    a_momentNames.resize(ncomp);
    for(int ivec = 0; ivec < ncomp; ivec++)
    {
      IvSpaceDim iv = IndMomSpaceDim::getIndex(ivec);
      string xint = string("x^") + EBArith::convertInt(iv[0]);
      string yint = string("y^") + EBArith::convertInt(iv[1]);
#if BL_SPACEDIM==3
      string zint = string("z^") + EBArith::convertInt(iv[2]);
#endif
      string integrand = D_TERM(xint, +yint, +zint);
      string name = string("int_cell") + integrand + string(" dV") ;
      a_momentNames[ivec] = name;
    }
  }
  void
  fillVolumeMoments(FabArray<EBCellFAB>& a_moments, 
                    const EBLevelGrid& a_eblg)
  {
    const BoxArray            & dbl   = a_eblg.getDBL();
    const EBISLayout          & ebisl = a_eblg.getEBISL();
    const DistributionMapping & dm    = a_eblg.getDM();
    for (MFIter mfi(dbl, dm); mfi.isValid(); ++mfi)
    {
      a_moments[mfi].setVal(0.);
      const EBISBox& ebisBox = ebisl[mfi];
      Box               grid =   dbl[mfi];
      IntVectSet ivs(grid);
      VoFIterator vofit(ivs, ebisBox.getEBGraph());
      Vector<VolIndex> vofs = vofit.getVector();
      for(int ivof = 0 ; ivof < vofs.size(); ivof++)
      {
        IndMomSpaceDim volmom = ebisBox.getEBData().getVolumeMoments(vofs[ivof]);
        for(MomItSpaceDim momit; momit.ok(); ++momit)
        {
          int ivar = IndMomSpaceDim::indexOf(momit());
          a_moments[mfi](vofs[ivof], ivar) = volmom[momit()];
        }
      }
    }
  }
/************/
  void
  sumFineMinusCoarse(FabArray<EBCellFAB>       & a_errorCoar, 
                     const FabArray<EBCellFAB> & a_momentsMedi, 
                     const FabArray<EBCellFAB> & a_momentsCoar, 
                     const EBLevelGrid         & a_eblgMedi,
                     const EBLevelGrid         & a_eblgCoar)
  {
    const BoxArray            & dbl   = a_eblgCoar.getDBL();
    const DistributionMapping & dm    = a_eblgCoar.getDM();
    const EBISLayout          & ebisl = a_eblgCoar.getEBISL();
    int ncomp = IndMomSpaceDim::size();
    for (MFIter mfi(dbl, dm); mfi.isValid(); ++mfi)
    {
      a_errorCoar[mfi].setVal(0.);
      const EBISBox& ebisBoxCoar = ebisl[mfi];
      Box               gridCoar =   dbl[mfi];
      IntVectSet ivs(gridCoar);
      VoFIterator vofit(ivs, ebisBoxCoar.getEBGraph());
      Vector<VolIndex> vofs = vofit.getVector();
      for(int ivof = 0 ; ivof < vofs.size(); ivof++)
      {
        VolIndex vofCoar = vofs[ivof];
        Vector<VolIndex> vofsFine = ebisBoxCoar.refine(vofCoar);
        Vector<Real> fineSum(ncomp, 0.0);
        sumFineValues(fineSum, a_momentsMedi[mfi], vofsFine);
        for(int ivar = 0; ivar < ncomp; ivar++)
        {
          a_errorCoar[mfi](vofCoar, ivar) = fineSum[ivar] - a_momentsCoar[mfi](vofCoar, ivar);
        }
      }
    }
  }
/************/
  void
  sphereConvTest()
  {
    //make layouts == domain
    Box domainBoxFine, domainBoxMedi, domainBoxCoar;
    RealVect dxFine, dxMedi, dxCoar;
    getFinestDomain(domainBoxFine, dxFine);
    dxMedi = 2.0*dxFine;
    dxCoar = 2.0*dxMedi;
    domainBoxMedi = coarsen(domainBoxFine, 2);
    domainBoxCoar = coarsen(domainBoxMedi, 2);

    BoxArray dblFine(domainBoxFine);
    BoxArray dblMedi(domainBoxMedi);
    BoxArray dblCoar(domainBoxCoar);

    DistributionMapping dmFine(dblFine);
    DistributionMapping dmMedi(dblMedi);
    DistributionMapping dmCoar(dblCoar);


    amrex::Print() << "=====================================================" << "\n";
    amrex::Print() << "=== Richardson Convergence Test of Volume Moments ===" << "\n";
    amrex::Print() << "=====================================================" << "\n";
    RealVect  sphereCenter;
    Real      sphereRadius;
    //doing a geometry convergence test so we have to make the geometry three times
    //this sort of test only works (because of multivalued issues) if the shape is concave.
    //since we know we are on the inside of a sphere, all is good
    makeGeometry(domainBoxFine,  dxFine, sphereCenter, sphereRadius);
    EBLevelGrid eblgFine(dblFine, dmFine, domainBoxFine, 1);
    makeGeometry(domainBoxMedi,  dxMedi, sphereCenter, sphereRadius);
    EBLevelGrid eblgMedi(dblMedi, dmMedi, domainBoxMedi, 1);
    makeGeometry(domainBoxCoar,  dxCoar, sphereCenter, sphereRadius);
    EBLevelGrid eblgCoar(dblCoar, dmCoar, domainBoxCoar, 1);

    EBCellFactory factFine(eblgFine.getEBISL());
    EBCellFactory factMedi(eblgMedi.getEBISL());
    EBCellFactory factCoar(eblgCoar.getEBISL());
    int ncomp = IndMomSpaceDim::size();
    int nghost = 0;
    //fill data holders with all the volume moments
    FabArray<EBCellFAB> momentsFine(dblFine, dmFine, ncomp, nghost, MFInfo(), factFine);
    FabArray<EBCellFAB> momentsMedi(dblMedi, dmMedi, ncomp, nghost, MFInfo(), factMedi);
    FabArray<EBCellFAB> momentsCoar(dblCoar, dmCoar, ncomp, nghost, MFInfo(), factCoar);
    fillVolumeMoments(momentsFine, eblgFine);
    fillVolumeMoments(momentsMedi, eblgMedi);
    fillVolumeMoments(momentsCoar, eblgCoar);

    //error = sum(fine) - coar (exact in this context)
    FabArray<EBCellFAB> errorMedi(dblMedi, dmMedi, ncomp, nghost, MFInfo(), factMedi);
    FabArray<EBCellFAB> errorCoar(dblCoar, dmCoar, ncomp, nghost, MFInfo(), factCoar);

    sumFineMinusCoarse(errorMedi, momentsFine, momentsMedi, eblgFine, eblgMedi);
    sumFineMinusCoarse(errorCoar, momentsMedi, momentsCoar, eblgMedi, eblgCoar);
    Vector<string> momentNames;
    fillMomentNames(momentNames);

    //compare the two errors.
    EBLevelDataOps::compareError(errorMedi, errorCoar, eblgMedi, eblgCoar, momentNames);

    amrex::Print() << "==============================================" << "\n" ;
  }
}
/************/
/************/
int
main(int argc, char* argv[])
{
  int retval = 0;
  amrex::Initialize(argc,argv);

  amrex::sphereConvTest();

  amrex::Finalize();
  return retval;
}
