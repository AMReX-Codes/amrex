
#include <iostream>
using std::cerr;

#include "AMReX_ParmParse.H"
#include "AMReX_EBCellFAB.H"
#include "AMReX_EBCellFactory.H"
#include "AMReX_EBNormalizeByVolumeFraction.H"
#include "AMReX_GeometryShop.H"
#include "AMReX_WrappedGShop.H"
#include "AMReX_AllRegularService.H"
#include "AMReX_PlaneIF.H"
#include "AMReX_EBLevelDataOps.H"
#include "AMReX_EBCellFactory.H"
#include "AMReX_SphereIF.H"

namespace amrex
{
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
        amrex::Print() << "using GeometryShop for geom gen" << endl;
        GeometryShop workshop(ramp,0, a_dx);
        //this generates the new EBIS
        EBIndexSpace* ebisPtr = AMReX_EBIS::instance();
        ebisPtr->define(a_domain, origin, a_dx, workshop);
      }
      else
      {
        amrex::Print() << "using GeometryShop for geom gen" << endl;
        WrappedGShop workshop(ramp,0, a_dx);
        //this generates the new EBIS
        EBIndexSpace* ebisPtr = AMReX_EBIS::instance();
        ebisPtr->define(a_domain, origin, a_dx, workshop);
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

  int testNormalization(int igeom)
  {
    Box domain;
    Real dx;
  
    int eekflag = makeGeometry(domain, dx, igeom);
    if(eekflag != 0) return eekflag;
    ParmParse pp;
    int maxboxsize;
    pp.get("maxboxsize", maxboxsize);
    BoxArray ba(domain);
    ba.maxSize(maxboxsize);
    DistributionMapping dm(ba);
    int nghost = 2;
    EBLevelGrid eblg(ba, dm, domain, 3);
    // Set up a field that is 1 on the regular cells and kappa on the
    // irregular cells.
    // NOTE: The following line gets an "invalid operator" error.
    int ncomp =1;
    EBCellFactory fact(eblg.getEBISL());
    FabArray<EBCellFAB>      phi(ba, dm,  ncomp, nghost, MFInfo(), fact);
    FabArray<EBCellFAB> kappaPhi(ba, dm,  ncomp, nghost, MFInfo(), fact);
    for(MFIter  mfi(phi); mfi.isValid(); ++mfi)
    {

      EBCellFAB& phiFAB =      phi[mfi];
      EBCellFAB& kapFAB = kappaPhi[mfi];
      phiFAB.setVal(1.0); // Set all cell values to 1.
      kapFAB.setVal(1.0); // Set all cell values to 1.
    }
    // Now go over the irregular cells and set them to the
    // volume fraction.
    EBLevelDataOps::kappaWeight(phi);
    EBLevelDataOps::kappaWeight(kappaPhi);

    // Now normalize phi by the volume fractions.
    EBNormalizeByVolumeFraction normalizor(eblg, nghost, 1);
    normalizor.normalize(phi, kappaPhi,0, ncomp);

    for(MFIter  mfi(phi); mfi.isValid(); ++mfi)
    {
      // Now verify that phi == 1 on the irregular cells.
      const Box&     grid = eblg.getDBL()  [mfi];
      const EBISBox& ebis = eblg.getEBISL()[mfi];
      EBCellFAB& phiFAB = phi[mfi];
      //EBCellFAB& kapFAB = kappaPhi[mfi];
      
      const IntVectSet& irregCells = ebis.getIrregIVS(grid);
      for (VoFIterator vit(irregCells, ebis.getEBGraph()); vit.ok(); ++vit)
      {
        VolIndex vi = vit();
        //Real volfrac = ebis.volFrac(vi);
        //Real kapVal = kapFAB(vi,0);
        
        Real val = phiFAB(vi,0);
        Real tol= 1.0e-10;

        if (std::abs(val - 1.0) > tol)
        {
          eekflag = -1;
          amrex::Print() << "FAIL: phi != 1 on irregular cell!" << endl;
          return eekflag;
        }
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
    retval = amrex::testNormalization(igeom);
    if(retval != 0)
    {
      amrex::Print() << "normalization test failed with code " << retval << "\n";
    }
  }

  amrex::Print() << "normalization test passed \n";

  amrex::Finalize();
  return retval;
}

