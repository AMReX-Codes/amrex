#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cmath>
#include <cstdio>
#include <iostream>

#include "ParmParse.H"
#include "EBArith.H"
#include "EBAMRIO.H"
#include "CH_HDF5.H"
#include "parstream.H"
#include "EBIndexSpace.H"
#include "EBISLayout.H"
#include "ParmParse.H"
#include "GeometryShop.H"
#include "PolyGeom.H"
#include "VoFIterator.H"
#include "SphereIF.H"
#include  "RealVect.H"
#include "UsingNamespace.H"

#ifdef CH_MPI
#include "CH_Attach.H"
#endif

/***************/
/***************/
void
getFinestDomain(Box&       a_domain,
                RealVect&      a_dx)
{
  ParmParse pp;
  Vector<int> n_cell(SpaceDim);
  pp.getarr("n_cell",n_cell,0,SpaceDim);

 CH_assert(n_cell.size() == SpaceDim);
  IntVect lo = IntVect::Zero;
  IntVect hi;
  for (int ivec = 0; ivec < SpaceDim; ivec++)
    {
      if (n_cell[ivec] <= 0)
        {
          pout() << " bogus number of cells input = " << n_cell[ivec];
          MayDay::Error();
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
  CH_XD::EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  int biggridsize;
  pp.get("max_grid_size", biggridsize);
  vector<Real>  sphereCenterVect(SpaceDim);
  pp.getarr("sphere_center",sphereCenterVect, 0, SpaceDim);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      a_sphereCenter[idir] = sphereCenterVect[idir];
    }

  pout() << "using a sphere implicit function" << endl;
  bool negativeInside = true;
  SphereIF lalaBall(a_sphereRadius, a_sphereCenter, negativeInside);
  GeometryShop workshop(lalaBall,0,a_dx);
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
  pout() << "===============================================" << endl;
  pout() <<  varstring << " at embedded boundary test " << endl;
  for (int inorm = 0; inorm <= 2; inorm++)
    {
      if (inorm == 0)
        {
          pout()  << "Using max norm." << endl;
        }
      else
        {
          pout()  << "Using L-" << inorm << " norm." << endl;
        }

      Real ebIrregNormCoar, ebIrregNormFine;

      int comp = 0;
      EBArith::irregNorm(ebIrregNormCoar,
                         a_errorCoar,
                         a_ivsIrregCoar,
                         a_ebisBoxCoar,
                         comp,  inorm);

      EBArith::irregNorm(ebIrregNormFine,
                         a_errorFine,
                         a_ivsIrregFine,
                         a_ebisBoxFine,
                         comp,  inorm);

      if (a_ivsIrregCoar.isEmpty())
        {
          pout() << "no irregular fine vofs" << endl;
        }
      else
        {
          pout() << varstring << " Error Norm Coar = " << ebIrregNormCoar << endl;
        }
      if (a_ivsIrregFine.isEmpty())
        {
          pout() << "no irregular fine vofs" << endl;
        }
      else
        {
          pout() <<  varstring << " Error Norm Fine = " << ebIrregNormFine << endl;
        }
#if defined(CH_USE_DOUBLE)
      Real eps   = 1.0e-12;
#elif defined(CH_USE_FLOAT)
      Real eps   = 1.0e-5;
#else
#error Unknown Chombo precision
#endif
      if ((Abs(ebIrregNormCoar) > eps) && (Abs(ebIrregNormFine) > eps))
        {
          Real order = log(ebIrregNormCoar/ebIrregNormFine)/log(2.0);
          pout() << "Order of " << varstring  <<" = " << order << endl << endl;
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
outputError(const BaseIVFAB<Real>& a_error,
            const IntVectSet&      a_ivsIrreg,
            const EBISBox&         a_ebisBox,
            bool                   a_isCoarse,
            int                    a_itype,
            const char* const      a_varname)
{
  char filename[100];
  if (a_isCoarse)
    {
     sprintf(filename,"%s%dCoar.hdf5", a_varname, a_itype);
    }
  else
    {
      sprintf(filename,"%s%dFine.hdf5",a_varname, a_itype);
    }

  Box domain = a_ebisBox.getDomain().domainBox();
  EBCellFAB fabData(a_ebisBox, domain, 3);
  fabData.setVal(0.);
  for (VoFIterator vofit(a_ivsIrreg, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      const Real error = a_error(vof, 0);
      const Real bndryArea = a_ebisBox.bndryArea(vof);
      fabData(vof, 0) = error;
      fabData(vof, 1) = error*bndryArea;
      fabData(vof, 2) = bndryArea;
    }

  Vector<string> names(3);
  char namewithba[80];
  sprintf(namewithba, "%sXBndryArea", a_varname);

  names[0] = string(a_varname);
  names[1] = string(namewithba);
  names[2] = string("BoundaryArea");
#ifdef CH_USE_HDF5
#ifndef CH_MPI
  Real dx=1.0;
  Real dt=1.0;
  Real time=1.0;
  Vector<Real> covval(3, 0.0);
  bool repcov = false;
  writeEBHDF5(filename, domain, fabData, names, domain, dx, dt, time, repcov, covval);
#endif
#endif
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
      Real error = Abs(dotProd) - 1;

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
  Vector<Box> boxFine(1, domainBoxFine);
  Vector<Box> boxCoar(1, domainBoxCoar);
  Vector<int> proc(1, 0);
  DisjointBoxLayout dblFine(boxFine, proc);
  DisjointBoxLayout dblCoar;
  coarsen(dblCoar, dblFine, 2);

  pout() << "==============================================" << endl;
  RealVect  sphereCenter;
  Real      sphereRadius;
  int itype = 0;

  makeGeometry(domainBoxFine,  dxFine, sphereCenter, sphereRadius);
  EBISLayout ebislFine, ebislCoar;
  const CH_XD::EBIndexSpace* const ebisPtrFine = Chombo_EBIS::instance();
  ebisPtrFine->fillEBISLayout(ebislFine, dblFine, domainBoxFine, 0);

  const CH_XD::EBIndexSpace* const ebisPtrCoar = Chombo_EBIS::instance();
  makeGeometry(domainBoxCoar,  dxCoar, sphereCenter, sphereRadius);
  ebisPtrCoar->fillEBISLayout(ebislCoar, dblCoar, domainBoxCoar, 0);

  int ifileout;
  ParmParse pp;
  pp.get("file_output", ifileout);
  //do the whole convergence test thing.
  bool failedTest = false;
  for (DataIterator dit = dblFine.dataIterator(); dit.ok(); ++dit)
    {
      const EBISBox& ebisBoxFine = ebislFine[dit()];
      const EBISBox& ebisBoxCoar = ebislCoar[dit()];
      IntVectSet ivsIrregFine = ebisBoxFine.getIrregIVS(domainBoxFine);
      IntVectSet ivsIrregCoar = ebisBoxCoar.getIrregIVS(domainBoxCoar);

      BaseIVFAB<Real> errorFine(ivsIrregFine, ebisBoxFine.getEBGraph(), 1);
      BaseIVFAB<Real> errorCoar(ivsIrregCoar, ebisBoxCoar.getEBGraph(), 1);

      getNormalDotTrueNormM1(errorFine, ivsIrregFine, ebisBoxFine, sphereCenter, dxFine);
      getNormalDotTrueNormM1(errorCoar, ivsIrregCoar, ebisBoxCoar, sphereCenter, dxCoar);

      failedTest |= compareError(errorFine, errorCoar,
                                 ivsIrregFine, ivsIrregCoar,
                                 ebisBoxFine, ebisBoxCoar, "Normal dotted with  true normal");

      if (ifileout == 1)
        {
          outputError(errorFine, ivsIrregFine, ebisBoxFine, false,itype, "sphNormDotTrueNormM1");
          outputError(errorCoar, ivsIrregCoar, ebisBoxCoar, true, itype, "sphNormDotTrueNormM1");
        }

      getNormalMinuTrueNorm(errorFine, ivsIrregFine, ebisBoxFine, sphereCenter, dxFine);
      getNormalMinuTrueNorm(errorCoar, ivsIrregCoar, ebisBoxCoar, sphereCenter, dxCoar);

      failedTest |= compareError(errorFine, errorCoar,
                                 ivsIrregFine, ivsIrregCoar,
                                 ebisBoxFine, ebisBoxCoar, "mag(Normal minus with  true normal)");

      if (ifileout == 1)
        {
          outputError(errorFine, ivsIrregFine, ebisBoxFine, false,itype, "sphNormMinuTrueNormM1");
          outputError(errorCoar, ivsIrregCoar, ebisBoxCoar, true, itype, "sphNormMinuTrueNormM1");
        }

      getCentroidDistError(errorFine, ivsIrregFine, ebisBoxFine, sphereCenter, dxFine, sphereRadius);
      getCentroidDistError(errorCoar, ivsIrregCoar, ebisBoxCoar, sphereCenter, dxCoar, sphereRadius);

      failedTest |= compareError(errorFine, errorCoar,
                                 ivsIrregFine, ivsIrregCoar,
                                 ebisBoxFine, ebisBoxCoar, "Centroid Dist. from Axis - rad");

      if (ifileout == 1)
        {
          outputError(errorFine, ivsIrregFine, ebisBoxFine, false,itype, "sphCentroidDistError");
          outputError(errorCoar, ivsIrregCoar, ebisBoxCoar, true, itype, "sphCentroidDistError");
        }
    }
  if (failedTest)
    {
      MayDay::Error("sphereConvTest failed because convergence rate not good enough");
    }
  pout() << "==============================================" << endl;
}
/************/
/************/
int
main(int a_argc, char* a_argv[])
{
  int retval = 0;
#ifdef CH_MPI
  MPI_Init(&a_argc,&a_argv);
  {
    // setChomboMPIErrorHandler();
#endif
    // Check for an input file

    ParmParse pp(0, NULL, NULL, "sphereconv.inputs");

    sphereConvTest();

    CH_XD::EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
#ifdef CH_USE_HDF5
    HDF5Handle handle("ebfile.hdf5", HDF5Handle::CREATE);
    ebisPtr->writeInfo(handle);
    handle.close();
#endif

    ebisPtr->clear();
#ifdef CH_MPI
  }
  MPI_Finalize();
#endif

  return retval;
}
