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

#include "AMReX_ParmParse.H"
#include "AMReX_EBCellFAB.H"
#include "AMReX_EBCellFactory.H"
#include "AMReX_EBNormalizeByVolumeFraction.H"
#include "AMReX_GeometryShop.H"
#include "AMReX_AllRegularService.H"
#include "AMReX_PlaneIF.H"
#include "AMReX_EBLevelDataOps.H"
#include "AMReX_EBCellFactory.H"
#include "AMReX_SphereIF.H"
#include "AMReX_EBDebugDump.H"
#include "AMReX_MeshRefine.H"
#include "AMReX_EBFastFR.H"

namespace amrex
{
/**********/
  int makeGeometry(const GridParameters& a_params)
  {
    int eekflag =  0;
    //parse input file
    ParmParse pp;
    RealVect origin = a_params.probLo;
    Real dxFinest = a_params.coarsestDx;
    Box domainFinest = a_params.coarsestDomain;
    for(int ilev = 1; ilev < a_params.numLevels; ilev++)
    {
      dxFinest         /= a_params.refRatio[ilev-1];
      domainFinest.refine(a_params.refRatio[ilev-1]);
    }

    EBIndexSpace* ebisPtr = AMReX_EBIS::instance();
    int whichgeom;
    pp.get("which_geom",whichgeom);
    if (whichgeom == 0)
    {
      //allregular
      amrex::Print() << "all regular geometry\n";
      AllRegularService regserv;
      ebisPtr->define(domainFinest, origin, dxFinest, regserv);
    }
    else if (whichgeom == 1)
    {
      amrex::Print() << "ramp geometry\n";
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

      GeometryShop workshop(ramp,0);
      //this generates the new EBIS
      ebisPtr->define(domainFinest, origin, dxFinest, workshop);
    }
    else if (whichgeom == 5)
    {
      amrex::Print() << "sphere geometry\n";
      Vector<Real> centervec(SpaceDim);
      Vector<int>  ncellsvec(SpaceDim);
      int maxgrid;
      ParmParse pp;
      Real radius;
      pp.getarr(  "n_cell"       , ncellsvec, 0, SpaceDim);
      pp.get(   "sphere_radius", radius);
      pp.get(   "max_grid_size", maxgrid);
      pp.getarr("sphere_center", centervec, 0, SpaceDim);
      RealVect center;
      for(int idir = 0; idir < SpaceDim; idir++)
      {
        center[idir] = centervec[idir];
      }
      bool insideRegular = false;
      SphereIF sphere(radius, center, insideRegular);
      GeometryShop gshop(sphere, 0);
      ebisPtr->define(domainFinest, origin, dxFinest, gshop);
    }
    else
    {
      //bogus which_geom
      amrex::Print() << " bogus which_geom input = " << whichgeom;
      eekflag = 33;
    }

    return eekflag;
  }
  int fluxRegTest()
  {
    GridParameters params;
    amrex::Print() << "forcing maxLevel = 1 \n";

    getGridParameters(params, 1, true);
    //and defines it using a geometryservice
    int eekflag = makeGeometry(params);
    if(eekflag != 0) return eekflag;

    Vector<EBLevelGrid> eblg;
    getAllIrregEBLG(eblg, params);
    for(int ilev = 0; ilev < 2; ilev ++)
    {
      amrex::Print() << "grids[" << ilev << "] = " << eblg[ilev].getDBL() << endl;
    }
    int nvar = 1;
    int nghost = 2;
    EBCellFactory factCoar(eblg[0].getEBISL());
    FabArray<EBCellFAB> dataCoar(eblg[0].getDBL(), eblg[0].getDM(), nvar, nghost, MFInfo(), factCoar);

    //set data  to zero
    for(MFIter mfi(eblg[0].getDBL(), eblg[0].getDM()); mfi.isValid(); ++mfi)
    {
      dataCoar[mfi].setVal(0.);
    }

    int nref = params.refRatio[0];
    EBFastFR fluxReg(eblg[1],
                     eblg[0],
                     nref, nvar);

    fluxReg.setToZero();
    Real scale = 1.0;
    Real fluxVal = 4.77;
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      for(MFIter mfi(eblg[0].getDBL(), eblg[0].getDM()); mfi.isValid(); ++mfi)
      {
        Box     grid = eblg[0].getDBL()  [mfi];
        EBISBox ebis = eblg[0].getEBISL()[mfi];
        EBFaceFAB flux(ebis, grid, idir, 1);
        flux.setVal(fluxVal);
        fluxReg.incrementCoarse(flux, scale, mfi, 0, 0, 1);
      }
      for(MFIter mfi(eblg[1].getDBL(), eblg[1].getDM()); mfi.isValid(); ++mfi)
      {
        Box     grid = eblg[1].getDBL()  [mfi];
        EBISBox ebis = eblg[1].getEBISL()[mfi];
        EBFaceFAB flux(ebis, grid, idir, 1);
        flux.setVal(fluxVal);
        fluxReg.incrementFine(flux, scale, mfi, 0, 0, 1);
      }
    }

    //now reflux and see if the data changes (it should not change because the fluxes are equal)
    fluxReg.reflux(dataCoar, scale, 0, 0, 1);
    Real tol = 1.0e-10;
    for(MFIter mfi(eblg[0].getDBL(), eblg[0].getDM()); mfi.isValid(); ++mfi)
    {
      Box     grid = eblg[0].getDBL()  [mfi];
      EBISBox ebis = eblg[0].getEBISL()[mfi];
      IntVectSet ivsbox(grid);
      for(VoFIterator vofit(ivsbox, ebis.getEBGraph()); vofit.ok(); ++vofit)
      {
        const VolIndex& vof = vofit();
        Real dataval = dataCoar[mfi](vof, 0);
        if(std::abs(dataval) > tol)
        {
          amrex::Print() << "reflux failed at vof " << vof.gridIndex() << endl;
          return -1;
        }
      }
    }
    return 0;
  }
}
/***************/
/***************/
  int
  main(int argc, char* argv[])
  {
    int retval = 0;
    amrex::Initialize(argc,argv);

    retval = amrex::fluxRegTest();
    if(retval != 0)
    {
      amrex::Print() << "flux register test failed with code " << retval << "\n";
    }
    else
    {
      amrex::Print() << "flux register test passed \n";
    }
    amrex::Finalize();
    return retval;
  }

