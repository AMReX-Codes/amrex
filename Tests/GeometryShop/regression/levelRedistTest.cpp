
#include "AMReX_BaseIVFactory.H"
#include "AMReX_EBIndexSpace.H"
#include "AMReX_EBISLayout.H"
#include "AMReX_BoxIterator.H"
#include "AMReX_ParmParse.H"
#include "AMReX_GeometryShop.H"
#include "AMReX_WrappedGShop.H"
#include "AMReX_EBCellFAB.H"
#include "AMReX_EBLevelGrid.H"
#include "AMReX_LayoutData.H"
#include "AMReX_EBCellFactory.H"
#include "AMReX_VoFIterator.H"
#include "AMReX_EBArith.H"
#include "AMReX_AllRegularService.H"
#include "AMReX_PlaneIF.H"
#include "AMReX_SphereIF.H"
#include "AMReX_SPMD.H"
#include "AMReX_Print.H"
#include "AMReX_EBLevelRedist.H"
#include "AMReX_RedistStencil.H"

namespace amrex
{
/***************/
  Real densityFunc(const IntVect& a_iv, const Box& a_box)
  {
    Real retval;
    int iside;
    Real problen = a_box.longside(iside);
    Real x = a_iv[0];
    Real y = a_iv[1];
    x /= problen;
    y /= problen;

    retval = 1.0 + x + y + x*x + y*y;
    return retval;
  }
/***************/
  Real massFunc(const IntVect& a_iv, const Box& a_box)
  {
    Real retval;
    int iside;
    Real problen = a_box.longside(iside);
    Real x = a_iv[0];
    Real y = a_iv[1];
    x /= problen;
    y /= problen;
    retval = x*x*x + y*y*y;

    return retval;
  }
/***************/
  int makeGeometry(Box& a_domain,
                   Real& a_dx,
                   int igeom)
  {
    int eekflag =  0;
    //parse input file
    ParmParse pp;
    RealVect origin = RealVect::Zero;
#if (CH_SPACEDIM==2)
    std::vector<int> n_cell(SpaceDim, 64);
#else
    std::vector<int> n_cell(SpaceDim, 16);
#endif

    BL_ASSERT(n_cell.size() == SpaceDim);
    IntVect lo = IntVect::TheZeroVector();
    IntVect hi;
    for (int ivec = 0; ivec < SpaceDim; ivec++)
    {
      if (n_cell[ivec] <= 0)
      {
        amrex::Print() << " bogus number of cells input = " << n_cell[ivec];
        return(-1);
      }
      hi[ivec] = n_cell[ivec] - 1;
    }

    a_domain.setSmall(lo);
    a_domain.setBig(hi);

    Real prob_hi;
    pp.get("prob_hi",prob_hi);
    a_dx = prob_hi/n_cell[0];

    int whichgeom;
    pp.get("which_geom",whichgeom);
    if (whichgeom == 0)
    {
      //allregular
      amrex::Print() << "all regular geometry" << "\n";
      AllRegularService regserv;
      EBIndexSpace* ebisPtr = AMReX_EBIS::instance();
      ebisPtr->define(a_domain, origin, a_dx, regserv);
    }
    else if (whichgeom == 1)
    {
      amrex::Print() << "ramp geometry" << "\n";
      int upDir;
      int indepVar;
      Real startPt;
      Real slope;
      pp.get("up_dir",upDir);
      pp.get("indep_var",indepVar);
      pp.get("start_pt", startPt);
      pp.get("ramp_slope", slope);

      RealVect normal = RealVect::Zero;
      normal[upDir] = 1.0;
      normal[indepVar] = -slope;

      RealVect point = RealVect::Zero;
      point[upDir] = -slope*startPt;

      bool normalInside = true;

      PlaneIF ramp(normal,point,normalInside);


      if(igeom == 0)
      {
        amrex::Print() << "using GeometryShop" << endl;
        GeometryShop workshop(ramp);
        //this generates the new EBIS
        EBIndexSpace* ebisPtr = AMReX_EBIS::instance();
        ebisPtr->define(a_domain, origin, a_dx, workshop);
      }
      else
      {
        amrex::Print() << "using WrappedGShop" << endl;
        WrappedGShop workshop(ramp);
        //this generates the new EBIS
        EBIndexSpace* ebisPtr = AMReX_EBIS::instance();
        ebisPtr->define(a_domain, origin, a_dx, workshop);
      }
    }
    else if (whichgeom == 5)
    {
      amrex::Print() << "sphere geometry\n";
      std::vector<Real> centervec(SpaceDim);
      std::vector<int>  ncellsvec(SpaceDim);
      Real radius;
      pp.getarr(  "n_cell"       , ncellsvec, 0, SpaceDim);
      pp.get(   "sphere_radius", radius);
      pp.getarr("sphere_center", centervec, 0, SpaceDim);
      RealVect center;
      for(int idir = 0; idir < SpaceDim; idir++)
      {
        center[idir] = centervec[idir];
      }
      int ireginside = 0;
      pp.get("inside_regular", ireginside);
      bool insideRegular = (ireginside == 1);
      SphereIF sphere(radius, center, insideRegular);
      GeometryShop workshop(sphere, 0);
      EBIndexSpace* ebisPtr = AMReX_EBIS::instance();
      ebisPtr->define(a_domain, origin, a_dx, workshop);
    }
    else
    {
      //bogus which_geom
      amrex::Print() << " bogus which_geom input = "
                     << whichgeom << "\n";
      eekflag = 33;
    }

    return eekflag;
  }
/***************/
  int testConservation(int igeom)
  {
    Box domain;
    Real dx;
    makeGeometry(domain, dx, igeom);
    int redistRad = 2;
    int maxboxsize;
    ParmParse pp;
    pp.get("redist_radius", redistRad);
    pp.get("maxboxsize", maxboxsize);
    BoxArray ba(domain);
    ba.maxSize(maxboxsize);
    DistributionMapping dm(ba);
    EBLevelGrid eblg(ba, dm, domain, 2);
    int ncomp = 1;

    EBCellFactory ebcellfact(eblg.getEBISL());

    BaseIVFactory<Real> baseivfact(eblg.getEBISL());
    FabArray<EBCellFAB>        solution(ba, dm,  ncomp, 0, MFInfo(), ebcellfact);
    FabArray<BaseIVFAB<Real> > massdiff(ba, dm,  ncomp, 0, MFInfo(), baseivfact);

    //initialize solution and mass difference
    //and add up the total mass in the system
    Real summassdiff = 0;
    Real sumsolmass = 0;
    for(MFIter  mfi(solution); mfi.isValid(); ++mfi)
    {
      Box grid = ba[mfi];
      EBISBox ebisBox = eblg.getEBISL()[mfi];
      IntVectSet ivsIrreg = ebisBox.getIrregIVS(grid);
      BaseIVFAB<Real>& massFAB = massdiff[mfi];
      for (VoFIterator vofit(ivsIrreg, ebisBox.getEBGraph());
           vofit.ok(); ++vofit)
      {
        const VolIndex& vof = vofit();
        const IntVect&  iv = vof.gridIndex();
        Real mass = massFunc(iv, domain);
        massFAB(vof, 0) = mass;
        summassdiff += mass;
      }

      EBCellFAB& solFAB = solution[mfi];
      IntVectSet ivsBox(grid);
      for (VoFIterator vofit(ivsBox, ebisBox.getEBGraph());
           vofit.ok(); ++vofit)
      {
        const VolIndex& vof = vofit();
        const IntVect&  iv = vof.gridIndex();
        Real density = densityFunc(iv, domain);
        Real volFrac = ebisBox.volFrac(vof);
        solFAB(vof, 0) = density;
        sumsolmass += density*volFrac;
      }
    }

    Real localSumOld = sumsolmass + summassdiff ;

    //now redistribute massdiff into solution
    EBLevelRedist distributor(eblg, ncomp, redistRad);
    distributor.setToZero();
    for(MFIter  mfi(solution); mfi.isValid(); ++mfi)
    {
      distributor.increment(massdiff[mfi], mfi, 0, 1);
    }
    distributor.redistribute(solution, 0, 0, 1);

    //now check that the solution has all the mass
    sumsolmass = 0;
    for(MFIter  mfi(solution); mfi.isValid(); ++mfi)
    {
      Box grid = ba[mfi];
      EBISBox ebisBox = eblg.getEBISL()[mfi];
      EBCellFAB& solFAB = solution[mfi];
      IntVectSet ivsBox(grid);
      for (VoFIterator vofit(ivsBox, ebisBox.getEBGraph());
           vofit.ok(); ++vofit)
      {
        const VolIndex& vof = vofit();
        Real density = solFAB(vof, 0);
        if (density > 1.0e-12)
        {
          Real volFrac = ebisBox.volFrac(vof);
          sumsolmass += density*volFrac;
        }
      }
    }
    Real localSumNew = sumsolmass;

    Real massTotOld = 0;
    Real massTotNew = 0;

    //gather what each processor thinks is the sum old and new
    std::vector<Real> sumVecOld, sumVecNew;
    int baseproc = 0;
    gather(sumVecOld, localSumOld, baseproc);
    gather(sumVecNew, localSumNew, baseproc);
    if (ParallelDescriptor::MyProc() == baseproc)
    {
      BL_ASSERT(sumVecOld.size() == ParallelDescriptor::NProcs());
      BL_ASSERT(sumVecNew.size() == ParallelDescriptor::NProcs());
      for (int ivec = 0; ivec < ParallelDescriptor::NProcs(); ivec++)
      {
        massTotOld += sumVecOld[ivec];
        massTotNew += sumVecNew[ivec];
      }
    }
    //broadcast the sum to all processors.
    broadcast(massTotOld, baseproc);
    broadcast(massTotNew, baseproc);

    amrex::Print() << "mass tot old = "  << massTotOld << endl;
    amrex::Print() << "mass tot new = "  << massTotNew << endl;
    int eekflag = 0;
    if (massTotOld > 1.0e-9)
    {
      Real relDiff = std::abs(massTotOld - massTotNew)/massTotOld;
      if (relDiff > 1.5e-6)
      {
        amrex::Print() << "doh! " << "\n";;
        amrex::Print() << "initial solution mass + diff      = "  << massTotOld << "\n";
        amrex::Print() << "total mass in solution after dist = "  << massTotNew << "\n";
        amrex::Print() << "relative difference               = "  << relDiff    << "\n"  ;
        eekflag  = 3;
      }
    }
    return eekflag;
  }
}
/***************/
int
main(int argc, char* argv[])
{
  int retval = 0;
  amrex::Initialize(argc,argv);

  for(int igeom = 0; igeom <= 1; igeom++)
  {
    retval = amrex::testConservation(igeom);
    if(retval != 0)
    {
      amrex::Print() << argv[0] << " test failed with code " << retval << "\n";
    }
  }

  amrex::Print() << argv[0] << " test passed \n";

  amrex::Finalize();
  return retval;
}
