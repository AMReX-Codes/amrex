
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

namespace amrex
{

  void
  makeEBIS(const Box&       a_domain,
           const RealVect&  a_dx)
  {
    Box finest_domain = a_domain;
    Real fine_dx = a_dx[0];
    std::unique_ptr<BaseIF> impfunc;
    ParmParse pp;

    amrex::Print() << "parabola + sphere geometry\n";
    Vector<PolyTerm> poly;

    PolyTerm mono;
    Real coef;
    IntVect powers;
    Real amplitude = 1;

    // y^2 term
    coef = amplitude;
    powers = IntVect::Zero;
    powers[1] = 2;

    mono.coef   = coef;
    mono.powers = powers;

    poly.push_back(mono);

#if BL_SPACEDIM==3
    // z^2 term
    coef = amplitude;
    powers = IntVect::Zero;
    powers[2] = 2;
    mono.coef   = coef;
    mono.powers = powers;
    poly.push_back(mono);
#endif
    // x term
    coef = -1.0;
    powers = IntVect::Zero;
    powers[0] = 1;
    mono.coef   = coef;
    mono.powers = powers;

    poly.push_back(mono);

    PolynomialIF mirror(poly,true);
    RealVect translation;
      
    for(int idir = 0; idir < SpaceDim; idir++)
    {
      int finesize = finest_domain.size()[idir];
      translation[idir] = 0.5*finesize*fine_dx;
    }
    RealVect center = translation;
    translation[0] = 0;

    TransformIF transform(mirror);
    transform.translate(translation);

    Real radius = 0.2*center[0];
    SphereIF sphere(radius, center, false);
    Vector<BaseIF*> funcs(2);
    funcs[0] = &transform;
    funcs[1] = &sphere;
    IntersectionIF implicit(funcs);
    impfunc.reset(implicit.newImplicitFunction());
 
    EBIndexSpace* ebisPtr = AMReX_EBIS::instance();
    int biggridsize;

    pp.get("max_grid_size", biggridsize);
    GeometryShop workshop(*impfunc);
    int ebmaxcoarsen = 0;
    RealVect origin = RealVect::Zero;
    ebisPtr->define(a_domain, origin, a_dx[0], workshop, biggridsize, ebmaxcoarsen);
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
    
    std::string filename = string("parabwithsphere.") + convertInt(SpaceDim) + "d.plt";
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
