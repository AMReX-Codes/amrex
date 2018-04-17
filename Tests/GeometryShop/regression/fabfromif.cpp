
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
#include "AMReX_SPMD.H"
#include "AMReX_Print.H"
#include "AMReX_VisMF.H"
#include "AMReX_MultiFab.H"
#include "AMReX_SphereIF.H"
#include "AMReX_EBFluxFactory.H"
#include "AMReX_EBFluxFAB.H"
#include "AMReX_EBCellFactory.H"
#include "AMReX_EBCellFAB.H"
#include "AMReX_IrregFABFactory.H"
#include "AMReX_BaseEBCellFactory.H"
#include "AMReX_IrregFAB.H"
#include "AMReX_EBDataVarMacros.H"
#include "AMReX_FabArrayIO.H"
#include "AMReX_EllipsoidIF.H"
#include "AMReX_parstream.H"

namespace amrex
{
/***************/
  int makeGeometry(Box& a_domain,
                   Real& a_dx,
                   int igeom)
  {
    int eekflag =  0;
    int maxbox;
    //parse input file

    ParmParse pp;
    pp.get("maxboxsize", maxbox);
    RealVect origin = RealVect::Zero;
    std::vector<int> n_cell;
    pp.getarr("n_cell", n_cell, 0, SpaceDim);

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
      ebisPtr->define(a_domain, origin, a_dx, regserv, maxbox);
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
        amrex::Print() << "using GeometryShop for geom  gen" << endl;
        GeometryShop workshop(ramp,0, a_dx);
        //this generates the new EBIS
        EBIndexSpace* ebisPtr = AMReX_EBIS::instance();
        ebisPtr->define(a_domain, origin, a_dx, workshop, maxbox);
      }
      else
      {
        amrex::Print() << "using WrappedGShop for geom  gen" << endl;
        WrappedGShop workshop(ramp,0, a_dx);
        //this generates the new EBIS
        EBIndexSpace* ebisPtr = AMReX_EBIS::instance();
        ebisPtr->define(a_domain, origin, a_dx, workshop, maxbox);
      }
    }
    else if(whichgeom == 5)
    {
      Real sphereRadius;
      RealVect sphereCenter;
      pp.get("sphere_radius", sphereRadius);
      EBIndexSpace* ebisPtr = AMReX_EBIS::instance();
      
      vector<Real> sphereCenterVect;
      pp.getarr("sphere_center",sphereCenterVect, 0, SpaceDim);
      for (int idir = 0; idir < SpaceDim; idir++)
      {
        sphereCenter[idir] = sphereCenterVect[idir];
      }

      amrex::Print() << "using a sphere implicit function" << "\n";

      bool negativeInside = false;
      SphereIF lalaBall(sphereRadius, sphereCenter, negativeInside);
      if(igeom == 0)
      {
        amrex::Print() << "using GeometryShop for geom  gen" << endl;
        GeometryShop workshop(lalaBall);
      
        ebisPtr->define(a_domain, origin, a_dx, workshop, maxbox);
      }
      else
      {
        amrex::Print() << "using WrappedGShop for geom  gen" << endl;
        WrappedGShop workshop(lalaBall);
      
        ebisPtr->define(a_domain, origin, a_dx, workshop, maxbox);
      }
    }
    else if(whichgeom == 55)
    {
      RealVect ellipseRadius;
      RealVect ellipseCenter;

      EBIndexSpace* ebisPtr = AMReX_EBIS::instance();
      vector<Real> ellipseCenterVect;
      vector<Real> ellipseRadiusVect;
      pp.getarr("ellipse_center",ellipseCenterVect, 0, SpaceDim);
      pp.getarr("ellipse_radius",ellipseRadiusVect, 0, SpaceDim);
      for (int idir = 0; idir < SpaceDim; idir++)
      {
        ellipseCenter[idir] = ellipseCenterVect[idir];
        ellipseRadius[idir] = ellipseRadiusVect[idir];
      }

      amrex::Print() << "using a ellipsoid implicit function" << "\n";

      bool negativeInside = false;
      EllipsoidIF lalaBall(ellipseRadius, ellipseCenter, negativeInside);
      if(igeom == 0)
      {
        amrex::Print() << "using GeometryShop for geom  gen" << endl;
        GeometryShop workshop(lalaBall);
      
        ebisPtr->define(a_domain, origin, a_dx, workshop, maxbox);
      }
      else
      {
        amrex::Print() << "using WrappedGShop for geom  gen" << endl;
        WrappedGShop workshop(lalaBall);
      
        ebisPtr->define(a_domain, origin, a_dx, workshop, maxbox);
      }

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

  /****/
  /***************/
  int testEBIO(int igeom)
  {
    Box domain;
    Real dx;
    //make the initial geometry
    amrex::Print() << "making EBIS" << endl;
    makeGeometry(domain, dx, igeom);
    Box nodeBox = domain.surroundingNodes();
    BoxArray            ba(nodeBox);
    DistributionMapping dm(ba);


    amrex::Print() << "filling multifab with implicit function values" << endl;
    MultiFab impfuncMF(ba, dm, 1, 0);
    const EBIndexSpace* const ebisPtr = AMReX_EBIS::instance();
    
    for(MFIter mfi(impfuncMF); mfi.isValid(); ++mfi)
    {
      ebisPtr->fillNodeFarrayBoxFromImplicitFunction(impfuncMF[mfi], 1, dx*RealVect::Unit);
    }
    
    VisMF::Write(impfuncMF, string("plt.mf.impfunc"));
    return 0;
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
    retval = amrex::testEBIO(igeom);
    if(retval != 0)
    {
      amrex::Print() <<  argv[0] << " test failed with code " << retval << "\n";
    }
  }
  amrex::Print() << argv[0] << " test passed \n";

  amrex::Finalize();
  return retval;
}
