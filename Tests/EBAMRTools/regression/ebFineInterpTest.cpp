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
#include "AMReX_EBFineInterp.H"

namespace amrex
{
/***************/
  Real exactFunc(const IntVect& a_iv,
                 const Real& a_dx)
  {
    Real retval;
    ParmParse pp;
    RealVect xloc;
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      xloc[idir] = a_dx*(Real(a_iv[idir]) + 0.5);
    }
    Real x = xloc[0];
    Real y = xloc[1];

    retval = (1.0 + x*x + y*y + x*x*x + y*y*y +x*y +x*x*y + y*y*x);
    //retval  = x;
    //retval = x*x;
    //retval = x*y;
    //retval = x*x*x;
    return retval;
  }
/***************/
  int getError(FabArray<EBCellFAB> & a_error,
               const EBLevelGrid& a_eblgFine,
               const EBLevelGrid& a_eblgCoar,
               Real dxFine, Real dxCoar,
               int a_interpOrder)
  {
    int nghost = 0;
    int nvar = 1;


    EBCellFactory factFine(a_eblgFine.getEBISL());
    EBCellFactory factCoar(a_eblgCoar.getEBISL());
    FabArray<EBCellFAB> phiFineExac(a_eblgFine.getDBL(), a_eblgFine.getDM(), nvar, nghost, MFInfo(), factFine);
    FabArray<EBCellFAB> phiFineCalc(a_eblgFine.getDBL(), a_eblgFine.getDM(), nvar, nghost, MFInfo(), factFine);
    FabArray<EBCellFAB> phiCoarExac(a_eblgCoar.getDBL(), a_eblgCoar.getDM(), nvar, nghost, MFInfo(), factCoar);
    FabArray<EBCellFAB> phiFineQuadratic(a_eblgFine.getDBL(), a_eblgFine.getDM(), nvar, nghost, MFInfo(), factFine);
    FabArray<EBCellFAB> phiFineLinear(a_eblgFine.getDBL(), a_eblgFine.getDM(), nvar, nghost, MFInfo(), factFine);

    for(MFIter mfi(a_eblgFine.getDBL(), a_eblgFine.getDM()); mfi.isValid(); ++mfi)
    {
      a_error[mfi].setVal(0.0);
      Box valid = a_eblgFine.getDBL()[mfi];
      IntVectSet ivsBox(valid);

      for (VoFIterator vofit(ivsBox, a_eblgFine.getEBISL()[mfi].getEBGraph());
           vofit.ok(); ++vofit)
      {
        const VolIndex& vof = vofit();
        Real rightAns = exactFunc(vof.gridIndex(), dxFine);
        phiFineExac[mfi](vof, 0) = rightAns;
      }
    }

    for(MFIter mfi(a_eblgCoar.getDBL(), a_eblgCoar.getDM()); mfi.isValid(); ++mfi)
    {
      Box valid = a_eblgCoar.getDBL()[mfi];
      IntVectSet ivsBox(valid);

      for (VoFIterator vofit(ivsBox, a_eblgCoar.getEBISL()[mfi].getEBGraph());
           vofit.ok(); ++vofit)
      {
        const VolIndex& vof = vofit();
        Real rightAns = exactFunc(vof.gridIndex(), dxCoar);
        phiCoarExac[mfi](vof, 0) = rightAns;
      }

    }
    int nref = 2;
    //interpolate phiC onto phiF
    bool slowMode = true;
    EBFineInterp aveOp(a_eblgFine, a_eblgCoar, nref, nghost, a_interpOrder, slowMode);
    aveOp.interpolate(phiFineCalc, phiCoarExac, 0, 0, nvar);

    //error = phiF - phiFExact
    for(MFIter mfi(a_eblgFine.getDBL(), a_eblgFine.getDM()); mfi.isValid(); ++mfi)
    {
      Box valid = a_eblgFine.getDBL()[mfi];
      IntVectSet ivsBox(valid);
      const EBCellFAB& calcFAB = phiFineCalc[mfi];
      const EBCellFAB& exacFAB = phiFineExac[mfi];
      for (VoFIterator vofit(ivsBox, a_eblgFine.getEBISL()[mfi].getEBGraph()); vofit.ok(); ++vofit)
      {
        const VolIndex& vof = vofit();
        Real calcval = calcFAB(vof, 0);
        Real exacval = exacFAB(vof, 0);
        Real diff = std::abs(calcval - exacval);

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
      ebisPtr->define(domainFinest, origin, dxFinest, gshop, maxgrid, -1);
    }
    else
    {
      //bogus which_geom
      amrex::Print() << " bogus which_geom input = " << whichgeom;
      eekflag = 33;
    }

    return eekflag;
  }
  int fineInterpTest()
  {
    GridParameters paramsCoar, paramsMedi, paramsFine;
    amrex::Print() << "forcing maxLevel = 0 \n";

    getGridParameters(paramsCoar, 0, true);
    paramsMedi = paramsCoar;
    paramsMedi.refine(2);
    paramsFine = paramsMedi;
    paramsFine.refine(2);

    Real dxFine = paramsFine.coarsestDx;
    Real dxMedi = paramsMedi.coarsestDx;
    Real dxCoar = paramsCoar.coarsestDx;

    //and defines it using a geometryservice
    int eekflag = makeGeometry(paramsFine);
    if(eekflag != 0) return eekflag;

    Array<EBLevelGrid> veblgCoar, veblgMedi, veblgFine;
    getAllIrregEBLG(veblgFine, paramsFine);
    getAllIrregEBLG(veblgMedi, paramsMedi);
    getAllIrregEBLG(veblgCoar, paramsCoar);

    EBLevelGrid eblgCoar = veblgCoar[0];
    EBLevelGrid eblgMedi = veblgMedi[0];
    EBLevelGrid eblgFine = veblgFine[0];
    amrex::Print() << "coarse  Grids  = " << eblgCoar.getDBL() << endl;
    amrex::Print() << "medium  Grids  = " << eblgMedi.getDBL() << endl;
    amrex::Print() << "finest  Grids  = " << eblgFine.getDBL() << endl;

    int nvar = 1; int nghost = 0;
    EBCellFactory factMedi(eblgMedi.getEBISL());
    EBCellFactory factFine(eblgFine.getEBISL());

    FabArray<EBCellFAB> errorFine(eblgFine.getDBL(), eblgFine.getDM(), nvar, nghost, MFInfo(), factFine);
    FabArray<EBCellFAB> errorMedi(eblgMedi.getDBL(), eblgMedi.getDM(), nvar, nghost, MFInfo(), factMedi);

    int maxOrder = 2; //got to get this to 2
    int minOrder = 0;
    ParmParse pp; 
    pp.query("max_order", maxOrder);
    pp.query("min_order", minOrder);
    for( int iorder = minOrder; iorder <= maxOrder; iorder++)
    {
      amrex::Print() << "testing interpolation with order = " << iorder << endl;
      getError(errorFine, eblgFine, eblgMedi, dxFine, dxMedi, iorder);
      getError(errorMedi, eblgMedi, eblgCoar, dxMedi, dxCoar, iorder);

      EBLevelDataOps::compareError(errorFine, errorMedi, eblgFine, eblgMedi);
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

  retval = amrex::fineInterpTest();
  if(retval != 0)
  {
    amrex::Print() << "interpolation test failed with code " << retval << "\n";
  }
  else
  {
    amrex::Print() << "interpolation test passed \n";
  }
  amrex::Finalize();
  return retval;
}

