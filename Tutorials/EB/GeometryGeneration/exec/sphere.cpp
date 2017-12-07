
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
#include "AMReX_SphereIF.H"
#include "CommonCode.H"
#include "DebugDump.H"
#include "WriteEBPlotFile.H"
#include "AMReX_RealVect.H"

namespace amrex
{
  void
  makeSphere(const Box&       a_domain,
             const RealVect&  a_dx)
  {
    RealVect  sphereCenter;
    Real      sphereRadius;
    //parse input file.  single level
    ParmParse pp;
    pp.get("sphere_radius", sphereRadius);
    EBIndexSpace* ebisPtr = AMReX_EBIS::instance();
    int biggridsize;
    pp.get("max_grid_size", biggridsize);
    vector<Real>  sphereCenterVect(SpaceDim);
    pp.getarr("sphere_center",sphereCenterVect, 0, SpaceDim);
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      sphereCenter[idir] = sphereCenterVect[idir];
    }

    amrex::Print() << "using a sphere implicit function" << "\n";

    bool negativeInside = true;
    SphereIF lalaBall(sphereRadius, sphereCenter, negativeInside);
    GeometryShop workshop(lalaBall);
    int ebmaxcoarsen = 0;
    RealVect origin = RealVect::Zero;
    ebisPtr->define(a_domain, origin, a_dx[0], workshop, biggridsize, ebmaxcoarsen);
  }
/************/
  void
  makeAndOutputSphere()
  {
    //make layouts == domain
    Box domainBox;
    RealVect dx;
    getFinestDomain(domainBox, dx);

    BoxArray ba(domainBox);
    DistributionMapping dm(ba);

    makeSphere(domainBox,  dx);
    EBLevelGrid eblg(ba, dm, domainBox, 2);
    
    MultiFab data(ba, dm, 1, 0);
    setDataToSomething(data, eblg);
    
    std::string filename = string("spheredata.") + convertInt(SpaceDim) + "d.plt";
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

  amrex::makeAndOutputSphere();

  amrex::Finalize();
  return retval;
}
