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

namespace amrex
{
/***************/
  void 
  irregNorm(Real& a_ebIrregNorm,
            const BaseIVFAB<Real>& a_ebiError,
            const IntVectSet& a_ivsIrreg,
            const EBISBox& a_ebisBox,
            const int& a_comp,
            const int& a_normtype)
  {
    BL_PROFILE("EBArith::irregNorm");
    BL_ASSERT(a_normtype >= 0);
    BL_ASSERT(a_normtype <= 2);
    a_ebIrregNorm = 0.0;
    if (a_normtype == 0)
    {
      for (VoFIterator vofit(a_ivsIrreg, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
      {
        const VolIndex& vof = vofit();
        Real bdarea = a_ebisBox.bndryArea(vof);
        if (bdarea > 0)
        {
          Real valVoF = a_ebiError(vof, a_comp);
          a_ebIrregNorm = std::max(std::abs(valVoF), a_ebIrregNorm);
        }
      }
    }
    else
    {
      //integral norm
      Real areaTot = 0.0;
      Real normTot = 0.0;
      for (VoFIterator vofit(a_ivsIrreg, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
      {
        const VolIndex& vof = vofit();
        Real valVoF = a_ebiError(vof, a_comp); 
        Real bdarea = a_ebisBox.bndryArea(vof);
        areaTot += bdarea;
        if (a_normtype == 1)
          normTot += std::abs(valVoF)*bdarea;
        else
          normTot += valVoF*valVoF*bdarea;
      }
      if (areaTot > 1.0e-8)
        normTot /= areaTot;
      if (a_normtype == 2)
        normTot = sqrt(normTot);
      a_ebIrregNorm = normTot;
    }

  }
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
    SphereIF lalaBall(a_sphereRadius, a_sphereCenter, negativeInside);
    GeometryShop workshop(lalaBall);
    int ebmaxcoarsen = 0;
    RealVect origin = RealVect::Zero;
    ebisPtr->define(a_domain, origin, a_dx[0], workshop, biggridsize, ebmaxcoarsen);
  }
  /************/
  bool
  compareError(const BaseIVFAB<Real>& a_errorFine,
               const BaseIVFAB<Real>& a_errorCoar,
               const IntVectSet&      a_ivsIrregFine,
               const IntVectSet&      a_ivsIrregCoar,
               const EBISBox&         a_ebisBoxFine,
               const EBISBox&         a_ebisBoxCoar,
               const char* const a_varname)
  {
    bool failedTest = false;
    string varstring(a_varname);
    amrex::Print() << "===============================================" << "\n";
    amrex::Print() <<  varstring << " at embedded boundary test " << "\n";
    for (int inorm = 0; inorm <= 2; inorm++)
    {
      if (inorm == 0)
      {
        amrex::Print()  << "Using max norm." << "\n";
      }
      else
      {
        amrex::Print()  << "Using L-" << inorm << " norm." << "\n";
      }

      Real ebIrregNormCoar, ebIrregNormFine;

      int comp = 0;
      irregNorm(ebIrregNormCoar,
                a_errorCoar,
                a_ivsIrregCoar,
                a_ebisBoxCoar,
                comp,  inorm);

      irregNorm(ebIrregNormFine,
                a_errorFine,
                a_ivsIrregFine,
                a_ebisBoxFine,
                comp,  inorm);

      if (a_ivsIrregCoar.isEmpty())
      {
        amrex::Print() << "no irregular fine vofs" << "\n";
      }
      else
      {
        amrex::Print() << varstring << " Error Norm Coar = " << ebIrregNormCoar << "\n";
      }
      if (a_ivsIrregFine.isEmpty())
      {
        amrex::Print() << "no irregular fine vofs" << "\n";
      }
      else
      {
        amrex::Print() <<  varstring << " Error Norm Fine = " << ebIrregNormFine << "\n";
      }

      Real eps   = 1.0e-9;

      if ((std::abs(ebIrregNormCoar) > eps) && (std::abs(ebIrregNormFine) > eps))
      {
        Real order = log(ebIrregNormCoar/ebIrregNormFine)/log(2.0);
        amrex::Print() << "Order of " << varstring  <<" = " << order << "\n" << "\n";
        if (order < 0.5)
        {
          failedTest = true;
        }
      }
    }
    return failedTest;
  }
/************/
  void
  getNormalDotTrueNormM1(BaseIVFAB<Real>&  a_error,
                         const IntVectSet& a_ivsIrreg,
                         const EBISBox&    a_ebisBox,
                         const RealVect&   a_sphereCenter,
                         const RealVect&   a_dx)
  {
    //dot normal with truenormal.  right answer == 1 so we subtract 1
    for (VoFIterator vofit(a_ivsIrreg, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      RealVect centroid = a_ebisBox.bndryCentroid(vof);
      const IntVect&   iv = vof.gridIndex();
      RealVect centroidPt;
      for (int idir = 0; idir < SpaceDim; idir++)
      {
        centroidPt[idir] = a_dx[idir]*(Real(iv[idir]) + 0.5 + centroid[idir]);
      }
      //true normal points at the center
      RealVect trueNorm = a_sphereCenter;
      trueNorm -= centroidPt;
      Real sum;
      PolyGeom::unifyVector(trueNorm, sum);

      RealVect normal = a_ebisBox.normal(vof);
      Real dotProd = PolyGeom::dot(trueNorm, normal);
      Real error = std::abs(dotProd) - 1;

      a_error(vof, 0) = error ;
    }
  }

  /************/
  void
  getNormalMinuTrueNorm(BaseIVFAB<Real>&  a_error,
                        const IntVectSet& a_ivsIrreg,
                        const EBISBox&    a_ebisBox,
                        const RealVect&   a_sphereCenter,
                        const RealVect&   a_dx)
  {
    //dot normal with truenormal.  right answer == 1 so we subtract 1
    for (VoFIterator vofit(a_ivsIrreg, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      RealVect centroid = a_ebisBox.bndryCentroid(vof);
      const IntVect&   iv = vof.gridIndex();
      RealVect centroidPt;
      for (int idir = 0; idir < SpaceDim; idir++)
      {
        centroidPt[idir] = a_dx[idir]*(Real(iv[idir]) + 0.5 + centroid[idir]);
      }
      //true normal points at the center
      RealVect trueNorm = a_sphereCenter;
      trueNorm -= centroidPt;
      Real sum;
      PolyGeom::unifyVector(trueNorm, sum);

      RealVect normal = a_ebisBox.normal(vof);
      RealVect errorVect = normal;
      errorVect -= trueNorm;
      Real dotProd = PolyGeom::dot(errorVect, errorVect);
      Real error = sqrt(dotProd);

      a_error(vof, 0) = error ;
    }
  }
  /************/
  void
  getCentroidDistError(BaseIVFAB<Real>&  a_error,
                       const IntVectSet& a_ivsIrreg,
                       const EBISBox&    a_ebisBox,
                       const RealVect&   a_sphereCenter,
                       const RealVect&   a_dx,
                       const Real&       a_sphereRadius)
  {
    //dot normal with axis cross truenormal.  right answer == 0
    for (VoFIterator vofit(a_ivsIrreg, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      RealVect centroid = a_ebisBox.bndryCentroid(vof);
      const IntVect&   iv = vof.gridIndex();
      RealVect centroidPt;
      for (int idir = 0; idir < SpaceDim; idir++)
      {
        centroidPt[idir] = a_dx[idir]*(Real(iv[idir]) + 0.5 + centroid[idir]);
      }

      RealVect ebRadVec = a_sphereCenter;
      ebRadVec -= centroidPt;

      Real ebRadius = PolyGeom::dot(ebRadVec, ebRadVec);
      ebRadius = sqrt(ebRadius);

      Real error = ebRadius - a_sphereRadius;
      a_error(vof, 0) = error;
    }
  }
/************/
  void
  sphereConvTest()
  {
    //make layouts == domain
    Box domainBoxFine, domainBoxCoar;
    RealVect dxFine, dxCoar;
    getFinestDomain(domainBoxFine, dxFine);
    dxCoar = 2.0*dxFine;
    domainBoxCoar = coarsen(domainBoxFine, 2);

    BoxArray dblFine(domainBoxFine);
    BoxArray dblCoar = dblFine;
    dblCoar.coarsen(2);
    DistributionMapping dmfine(dblFine);
    DistributionMapping dmcoar(dblCoar);


    amrex::Print() << "==============================================" << "\n";
    RealVect  sphereCenter;
    Real      sphereRadius;

    makeGeometry(domainBoxFine,  dxFine, sphereCenter, sphereRadius);
    EBISLayout ebislFine, ebislCoar;
    const EBIndexSpace* const ebisPtr = AMReX_EBIS::instance();
    ebisPtr->fillEBISLayout(ebislFine, dblFine, dmfine, domainBoxFine, 0);

    makeGeometry(domainBoxCoar,  dxCoar, sphereCenter, sphereRadius);
    ebisPtr->fillEBISLayout(ebislCoar, dblCoar, dmcoar,  domainBoxCoar, 0);

    //do the whole convergence test thing.
    bool failedTest = false;

    BaseIVFactory<Real> factFine(ebislFine);
    BaseIVFactory<Real> factCoar(ebislCoar);
    DistributionMapping dmFine(dblFine);
    DistributionMapping dmCoar(dblFine);
    FabArray<BaseIVFAB<Real> >  errFine(dblFine, dmFine, 1, 0, MFInfo(), factFine);
    FabArray<BaseIVFAB<Real> >  errCoar(dblCoar, dmCoar, 1, 0, MFInfo(), factCoar);

    for (MFIter mfi(errFine); mfi.isValid(); ++mfi)
    {
      const EBISBox& ebisBoxFine = ebislFine[mfi];
      const EBISBox& ebisBoxCoar = ebislCoar[mfi];
      IntVectSet ivsIrregFine = ebisBoxFine.getIrregIVS(domainBoxFine);
      IntVectSet ivsIrregCoar = ebisBoxCoar.getIrregIVS(domainBoxCoar);

      BaseIVFAB<Real> errorFine = errFine[mfi];
      BaseIVFAB<Real> errorCoar = errCoar[mfi];

      getNormalDotTrueNormM1(errorFine, ivsIrregFine, ebisBoxFine, sphereCenter, dxFine);
      getNormalDotTrueNormM1(errorCoar, ivsIrregCoar, ebisBoxCoar, sphereCenter, dxCoar);

      failedTest |= compareError(errorFine, errorCoar,
                                 ivsIrregFine, ivsIrregCoar,
                                 ebisBoxFine, ebisBoxCoar, "Normal dotted with  true normal");


      getNormalMinuTrueNorm(errorFine, ivsIrregFine, ebisBoxFine, sphereCenter, dxFine);
      getNormalMinuTrueNorm(errorCoar, ivsIrregCoar, ebisBoxCoar, sphereCenter, dxCoar);

      failedTest |= compareError(errorFine, errorCoar,
                                 ivsIrregFine, ivsIrregCoar,
                                 ebisBoxFine, ebisBoxCoar, "mag(Normal minus with  true normal)");

      getCentroidDistError(errorFine, ivsIrregFine, ebisBoxFine, sphereCenter, dxFine, sphereRadius);
      getCentroidDistError(errorCoar, ivsIrregCoar, ebisBoxCoar, sphereCenter, dxCoar, sphereRadius);

      failedTest |= compareError(errorFine, errorCoar,
                                 ivsIrregFine, ivsIrregCoar,
                                 ebisBoxFine, ebisBoxCoar, "Centroid Dist. from Axis - rad");

    }
    amrex::Print() << "sphereConvTest passed \n";

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
