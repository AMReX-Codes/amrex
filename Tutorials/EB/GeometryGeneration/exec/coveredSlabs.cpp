
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
#include "AMReX_VoFIterator.H"
#include "AMReX_PlaneIF.H"
#include "AMReX_SphereIF.H"
#include "AMReX_IntersectionIF.H"
#include "AMReX_UnionIF.H"
#include "CommonCode.H"
#include "DebugDump.H"
#include "AMReX_LatheIF.H"
#include "AMReX_ComplementIF.H"
#include "AMReX_TransformIF.H"
#include "WriteEBPlotFile.H"
#include "AMReX_RealVect.H"
#include <AMReX_PolynomialIF.H>
#include <AMReX_EBArith.H>
#include <AMReX_CoveredSlabs.H>

namespace amrex
{

  void
  makeEBIS(const Box&       a_domain,
           const RealVect&  a_dx)
  {
    ParmParse pp;

    int ncovered;
    pp.get("num_covered_boxes", ncovered);
    Vector<Box> coveredBoxes(ncovered);
    for(int icovered = 0; icovered < ncovered; icovered++)
    {
      string boxstr = EBArith::convertInt(icovered);
      string loString = string("cov_box_lo_") + boxstr;
      string hiString = string("cov_box_hi_") + boxstr;
      Vector<int> lo(SpaceDim), hi(SpaceDim);
      pp.getarr(loString.c_str(), lo, 0, SpaceDim);
      pp.getarr(hiString.c_str(), hi, 0, SpaceDim);
      IntVect ivlo, ivhi;
      for(int idir = 0; idir < SpaceDim; idir++)
      {
        ivlo[idir] = lo[idir];
        ivhi[idir] = hi[idir];
      }
      coveredBoxes[icovered] = Box(ivlo,ivhi);
    }
    EBIndexSpace* ebisPtr = AMReX_EBIS::instance();
    int biggridsize;
    pp.get("max_grid_size", biggridsize);
    CoveredSlabs slabs(coveredBoxes);

    int ebmaxcoarsen = 0;
    RealVect origin = RealVect::Zero;
    ebisPtr->define(a_domain, origin, a_dx[0], slabs, biggridsize, ebmaxcoarsen);
  }
/************/
  void
  makeAndOutputEBIS()
  {
    //make layouts == domain
    Box domainBox;
    RealVect dx;
    getFinestDomain(domainBox, dx);

    BoxArray ba(domainBox);
    DistributionMapping dm(ba);

    makeEBIS(domainBox,  dx);
    EBLevelGrid eblg(ba, dm, domainBox, 2);
    
    MultiFab data(ba, dm, 1, 0);
    setDataToSomething(data, eblg);
    
    std::string filename = string("coveredslabs.") + convertInt(SpaceDim) + "d.plt";
    WriteSingleLevelEBPlotFile(filename, data, eblg, Vector<string>());

  }
}
/************/
/************/
int
main(int argc, char* argv[])
{
  int retval = 0;
  amrex::Initialize(argc,argv);

  amrex::makeAndOutputEBIS();

  amrex::Finalize();
  return retval;
}
