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

#include "AMReX_EBIndexSpace.H"
#include "AMReX_EBISLayout.H"
#include "AMReX_BoxIterator.H"
#include "AMReX_ParmParse.H"
#include "AMReX_GeometryShop.H"
#include "AMReX_PolyGeom.H"
#include "AMReX_EBCellFAB.H"
#include "AMReX_EBCellFactory.H"
#include "AMReX_VoFIterator.H"
#include "AMReX_AllRegularService.H"
#include "AMReX_PlaneIF.H"
#include "AMReX_SphereIF.H"
#include "AMReX_MeshRefine.H"
#include "AMReX_EBDebugDump.H"
#include "AMReX_EBDebugOut.H"
#include "AMReX_EBLevelDataOps.H"
#include "AMReX_EBCFInterp.H"

namespace amrex
{
/***************/
  Real exactFunc(const IntVect& a_iv,
                 const Real& a_dx)
  {
    Real retval;
    Real probHi;
    ParmParse pp;
    pp.get("domain_length",probHi);
    RealVect xloc;
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      xloc[idir] = a_dx*(Real(a_iv[idir]) + 0.5);
    }
    Real x = xloc[0]/probHi;
    Real y = xloc[1]/probHi;

    retval = 1.0 + x*x + y*y + 0.1*x*x*x + 0.1*y*y*y + 0.01*x*x*x*x + 0.01*y*y*y*y;
    //retval = a_time;
    //  retval = 2.*xloc[0];
    //retval  = x;
    //  retval = xloc[0]*xloc[0];
    return retval;
  }

/***************/
  int getError(FabArray<EBCellFAB> & a_error,
               Array<EBLevelGrid> & a_eblg,
               const Real& a_dxLev1)
  {
    int nghostData = 3;
    int nghost = 2;
    int nvar   = 1;
    EBLevelGrid eblgFine = a_eblg[1];
    EBLevelGrid eblgCoar = a_eblg[0];
    Real dxCoar = a_dxLev1*2.;
    EBCellFactory factFine(a_eblg[1].getEBISL());
    EBCellFactory factCoar(a_eblg[0].getEBISL());
    FabArray<EBCellFAB> phiFineCalc(eblgFine.getDBL(), eblgFine.getDM(), nvar, nghostData, MFInfo(), factFine);
    FabArray<EBCellFAB> phiFineExac(eblgFine.getDBL(), eblgFine.getDM(), nvar, nghostData, MFInfo(), factFine);
    FabArray<EBCellFAB> phiCoarExac(eblgCoar.getDBL(), eblgCoar.getDM(), nvar, nghostData, MFInfo(), factCoar);

    for(MFIter mfi(eblgFine.getDBL(), eblgFine.getDM()); mfi.isValid(); ++mfi)
    {
      Box valid = eblgFine.getDBL()[mfi];
      Box grownBox = grow(valid, nghostData);
      grownBox &= eblgFine.getDomain();
      IntVectSet ivsBox(grownBox);
      phiFineCalc[mfi].setVal(0.0);
      phiFineExac[mfi].setVal(0.0);
      a_error[mfi].setVal(0.0);

      for (VoFIterator vofit(ivsBox, eblgFine.getEBISL()[mfi].getEBGraph());
           vofit.ok(); ++vofit)
      {
        const VolIndex& vof = vofit();
        Real rightAns = exactFunc(vof.gridIndex(), a_dxLev1);
        if(valid.contains(vof.gridIndex()))
        {
          phiFineCalc[mfi](vof, 0) = rightAns;
        }
        phiFineExac[mfi](vof, 0) = rightAns;
      }
    }

    int ibox = 0;
    for(MFIter mfi(eblgCoar.getDBL(), eblgCoar.getDM()); mfi.isValid(); ++mfi)
    {
      Box valid = eblgCoar.getDBL()[mfi];
      Box grownBox = grow(valid, nghostData);
      grownBox &= eblgCoar.getDomain();

      IntVectSet ivsBox(grownBox);
      EBCellFAB& phiCoarExacFAB = phiCoarExac[mfi];
      phiCoarExacFAB.setCoveredCellVal(0.0,0);
      for (VoFIterator vofit(ivsBox, eblgCoar.getEBISL()[mfi].getEBGraph());
           vofit.ok(); ++vofit)
      {
        const VolIndex& vof = vofit();
        Real rightAns = exactFunc(vof.gridIndex(), dxCoar);
        phiCoarExacFAB(vof, 0) = rightAns;
      }
//      dumpEBFAB(&phiCoarExacFAB);
      ibox++;
    }
    int nref = 2;
    //interpolate phiC onto phiF
    EBCFInterp interpOp(eblgFine, eblgCoar, nref, nghostData, nghost);


    interpOp.coarseFineInterp(phiFineCalc,  phiCoarExac, 0, 0, 1);

    //error = phiF - phiFExact
    for(MFIter mfi(eblgFine.getDBL(), eblgFine.getDM()); mfi.isValid(); ++mfi)
    {
      Box valid = eblgFine.getDBL()[mfi];
      Box grownBox = grow(valid, nghost);
      grownBox &= eblgFine.getDomain();
      IntVectSet ivsBox(grownBox);
      EBCellFAB& calcFAB = phiFineCalc[mfi];
      EBCellFAB& exacFAB = phiFineExac[mfi];
      for (VoFIterator vofit(ivsBox, eblgFine.getEBISL()[mfi].getEBGraph()); vofit.ok(); ++vofit)
      {
        const VolIndex& vof = vofit();
        Real diff = std::abs(calcFAB(vof,0)-exacFAB(vof,0));
        a_error[mfi](vof, 0) = diff;
      }
    }

    return 0;
  }
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
      Array<Real> centervec(SpaceDim);
      Array<int>  ncellsvec(SpaceDim);
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
/***************/
  int
  ebpwlfpTest()
  {
    GridParameters paramsFine, paramsCoar;
    amrex::Print() << "forcing maxLevel = 1 \n";

    getGridParameters(paramsCoar, 1, true);
    paramsFine = paramsCoar;
    paramsFine.refine(2);
    //and defines it using a geometryservice
    int eekflag = makeGeometry(paramsFine);

    Array<EBLevelGrid> eblgFine, eblgCoar;
    getAllIrregEBLG(eblgFine, paramsFine);
    getAllIrregEBLG(eblgCoar, paramsCoar);
    
    int nvar = 1;
    EBCellFactory factFine(eblgFine[1].getEBISL());
    EBCellFactory factCoar(eblgCoar[1].getEBISL());
    int nghost = 2;

    FabArray<EBCellFAB> errorFine(eblgFine[1].getDBL(), eblgFine[1].getDM(), nvar, nghost, MFInfo(), factFine);
    FabArray<EBCellFAB> errorCoar(eblgCoar[1].getDBL(), eblgCoar[1].getDM(), nvar, nghost, MFInfo(), factCoar);

    Real dxLev1Fine = paramsFine.coarsestDx/paramsFine.refRatio[0];
    Real dxLev1Coar = paramsCoar.coarsestDx/paramsCoar.refRatio[0];
    
    eekflag = getError(errorFine,eblgFine, dxLev1Fine);
    eekflag = getError(errorCoar,eblgCoar, dxLev1Coar);

    EBLevelDataOps::compareError(errorFine, errorCoar, eblgFine[1], eblgCoar[1], Array<string>(), true);
    return eekflag;
  }
}
/***************/
int
main(int argc, char* argv[])
{
  int retval = 0;
  amrex::Initialize(argc,argv);

  int eek = amrex::ebpwlfpTest();
  if(eek == 0)
  {
    amrex::Print() << "aggpwlfpTest passed\n";
  }
  else
  {
    amrex::Print() << "aggpwlfpTest failed with code " << eek << "\n";
  }

  amrex::Finalize();
  return retval;
}
/************/
