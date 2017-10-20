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
#include "AMReX_GradientOp.H"

namespace amrex
{
/***************/
  Real getPhiExact(const RealVect& a_loc)
  {
    Real x = a_loc[0];
    Real y = a_loc[1];

    Real retval = x*x + y*y;

    return retval;
  }
  RealVect getGradExact(const RealVect& a_loc)
  {
    Real x = a_loc[0];
    Real y = a_loc[1];

    RealVect retval = RealVect::Zero;
    retval[0] = 2.*x;
    retval[1] = 2.*y;
    return retval;
  }
/***************/
  int getError(FabArray<EBCellFAB> & a_error,
               const EBLevelGrid   & a_eblg,
               bool a_useLimiting)
  {
    int nghost = 1;
    Real dx = 1.0/a_eblg.getDomain().size()[0];
    EBCellFactory cellfact(a_eblg.getEBISL());


    FabArray<EBCellFAB> phiExac(a_eblg.getDBL(), a_eblg.getDM(),        1, nghost, MFInfo(), cellfact);
    FabArray<EBCellFAB> gphCalc(a_eblg.getDBL(), a_eblg.getDM(), SpaceDim, nghost, MFInfo(), cellfact);
    FabArray<EBCellFAB> gphExac(a_eblg.getDBL(), a_eblg.getDM(), SpaceDim, nghost, MFInfo(), cellfact);

    for(MFIter mfi(a_eblg.getDBL(), a_eblg.getDM()); mfi.isValid(); ++mfi)
    {
      Box     grid  = a_eblg.getDBL()[mfi];
      EBISBox ebis  = a_eblg.getEBISL()[mfi];
      Box grownBox = grow(grid, 1);
      Box grownBoxDom = grownBox & ebis.getDomain();
      IntVectSet ivsBox(grownBoxDom);

      for (VoFIterator vofit(ivsBox, a_eblg.getEBISL()[mfi].getEBGraph());
           vofit.ok(); ++vofit)
      {
        const VolIndex& vof = vofit();
        RealVect loc = EBArith::getVoFLocation(vof, dx, RealVect::Zero);
        Real phiexac = getPhiExact(loc);
        phiExac[mfi](vof, 0) = phiexac;
        RealVect gphexac = getGradExact(loc);
        for(int vecDir = 0; vecDir < SpaceDim; vecDir++)
        {
          gphExac[mfi](vof, vecDir) = gphexac[vecDir];
        }
      }

      //need to set ghost cells outside the domain because we use them
      IntVectSet ivsDom(grownBox);
      ivsDom -= grownBoxDom;
      BaseFab<Real>& regPhi = phiExac[mfi].getSingleValuedFAB();
      BaseFab<Real>& regGph = gphExac[mfi].getSingleValuedFAB();
      for(IVSIterator ivsit(ivsDom); ivsit.ok(); ++ivsit)
      {
        const IntVect& iv = ivsit();
        RealVect loc = EBArith::getIVLocation(iv, dx*RealVect::Unit, RealVect::Zero);
        Real phiexac = getPhiExact(loc);
        regPhi(iv, 0) = phiexac;
        RealVect gphexac = getGradExact(loc);
        for(int vecDir = 0; vecDir < SpaceDim; vecDir++)
        {
          regGph(iv, vecDir) = gphexac[vecDir];
        }
      }
      
      bool slowMode = false;
      GradientOp gradOp(a_eblg, dx, 1, nghost, a_useLimiting, slowMode);
      gradOp.gradient(gphCalc, phiExac);

      EBLevelDataOps::setVal(a_error, 0.);
      //error = calc - exact
      EBLevelDataOps::incr(a_error, gphCalc, 1.0);
      EBLevelDataOps::incr(a_error, gphExac,-1.0);

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
    Box domainFinest = a_params.coarsestDomain;
    Real dxFinest = 1.0/domainFinest.size()[0];
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
  int gradOpTest()
  {
    GridParameters paramsCoar, paramsMedi, paramsFine;
    amrex::Print() << "forcing maxLevel = 0 \n";

    getGridParameters(paramsCoar, 0, true);
    paramsFine = paramsCoar;
    paramsFine.refine(2);

    //and defines it using a geometryservice
    int eekflag = makeGeometry(paramsFine);
    if(eekflag != 0) return eekflag;

    Vector<EBLevelGrid> veblgCoar, veblgFine;
    getAllIrregEBLG(veblgCoar, paramsCoar);
    getAllIrregEBLG(veblgFine, paramsFine);

    EBLevelGrid eblgCoar = veblgCoar[0];
    EBLevelGrid eblgFine = veblgFine[0];
    amrex::Print() << "coarse  Grids  = " << eblgCoar.getDBL() << endl;
    amrex::Print() << "finest  Grids  = " << eblgFine.getDBL() << endl;

    int nvar = SpaceDim; int nghost = 1;
    EBCellFactory factFine(eblgFine.getEBISL());
    EBCellFactory factCoar(eblgCoar.getEBISL());
    FabArray<EBCellFAB> errorFine(eblgFine.getDBL(), eblgFine.getDM(), nvar, nghost, MFInfo(), factFine);
    FabArray<EBCellFAB> errorCoar(eblgCoar.getDBL(), eblgCoar.getDM(), nvar, nghost, MFInfo(), factCoar);

    for(int dolim = 0; dolim <= 1; dolim++)
    {
      bool doLimiting = (dolim==1);
      if(doLimiting)
      {
        amrex::Print() << "checking gradient operator with limiting turned on." << endl;
      }
      else
      {
        amrex::Print() << "checking gradient operator with limiting turned off." << endl;
      }
      getError(errorFine, eblgFine, doLimiting);
      getError(errorCoar, eblgCoar, doLimiting);

      EBLevelDataOps::compareError(errorFine, errorCoar, eblgFine, eblgCoar);
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

  retval = amrex::gradOpTest();
  if(retval != 0)
  {
   amrex::Print() << "gradient operator test failed with code " << retval << "\n";
  }
  else
  {
    amrex::Print() << "gradient operator  test passed \n";
  }
  amrex::Finalize();
  return retval;
}

