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
#include "AMReX_DivergenceOp.H"

namespace amrex
{
/***************/
  Real divFExact(const RealVect& a_loc)
  {
    Real x = a_loc[0];
    Real y = a_loc[1];

    Real retval = 3.*x*x + 3.*y*y;

    return retval;
  }
  RealVect fluxExact(const RealVect& a_loc)
  {
    Real x = a_loc[0];
    Real y = a_loc[1];

    RealVect retval = RealVect::Zero;
    retval[0] = x*x*x;
    retval[1] = y*y*y;
    return retval;
  }
/***************/
  int getError(FabArray<EBCellFAB> & a_error,
               const EBLevelGrid   & a_eblg)
  {
    int nghost = 1;
    int nvar = 1;
    Real dx = 1.0/a_eblg.getDomain().size()[0];
    EBCellFactory cellfact(a_eblg.getEBISL());
    EBFluxFactory fluxfact(a_eblg.getEBISL());

    FabArray<EBCellFAB> divFExac(a_eblg.getDBL(), a_eblg.getDM(), nvar, nghost, MFInfo(), cellfact);
    FabArray<EBCellFAB> divFCalc(a_eblg.getDBL(), a_eblg.getDM(), nvar, nghost, MFInfo(), cellfact);
    FabArray<EBFluxFAB> fluxExac(a_eblg.getDBL(), a_eblg.getDM(), nvar, nghost, MFInfo(), fluxfact);

    for(MFIter mfi(a_eblg.getDBL(), a_eblg.getDM()); mfi.isValid(); ++mfi)
    {
      Box     grid  = a_eblg.getDBL()[mfi];
      EBISBox ebis  = a_eblg.getEBISL()[mfi];
      IntVectSet ivsBox(grid);

      for (VoFIterator vofit(ivsBox, a_eblg.getEBISL()[mfi].getEBGraph());
           vofit.ok(); ++vofit)
      {
        const VolIndex& vof = vofit();
        RealVect centroid = ebis.centroid(vof);
        centroid *= dx;
        RealVect loc = EBArith::getVoFLocation(vof, dx, centroid);
        Real rightAns = divFExact(loc);
        divFExac[mfi](vof, 0) = rightAns;
      }
      
      IntVectSet ivsIrreg = ebis.getIrregIVS(grid);
      for (VoFIterator vofit(ivsIrreg, a_eblg.getEBISL()[mfi].getEBGraph());
           vofit.ok(); ++vofit)
      {
        const VolIndex& vof = vofit();
        RealVect centroid = ebis.bndryCentroid(vof);
        centroid *= dx;
        RealVect loc = EBArith::getVoFLocation(vof, dx, centroid);
        RealVect rightFlux = fluxExact(loc);
        RealVect normal = ebis.normal(vof);
        Real rightAns = 0;
        for(int idir = 0; idir < SpaceDim; idir++)
        {
          //- because we have the insane inward facing normal
          rightAns -= rightFlux[idir]*normal[idir];
        }
        fluxExac[mfi].getEBFlux()(vof, 0) = rightAns;
      }
      for(int faceDir = 0; faceDir < SpaceDim; faceDir++)
      {
        for(FaceIterator faceit(ivsBox, ebis.getEBGraph(),faceDir, FaceStop::SurroundingWithBoundary); faceit.ok(); ++faceit)
        {
          const FaceIndex& face = faceit();
          RealVect centroid = ebis.centroid(face);
          centroid[faceDir] = 0.;
          centroid *= dx;
          RealVect loc = EBArith::getFaceLocation(face, dx*RealVect::Unit, centroid);
          RealVect rightFlux = fluxExact(loc);
          Real rightAns = rightFlux[faceDir];
          fluxExac[mfi][faceDir](face, 0) = rightAns;
        }
      }
    }
    DivergenceOp divOp(a_eblg, dx, nvar, nghost);
    divOp.hybridDivergence(divFCalc, fluxExac, 0, 0, 1);

    //error = divFcalc - divFExact
    for(MFIter mfi(a_eblg.getDBL(), a_eblg.getDM()); mfi.isValid(); ++mfi)
    {
      Box valid = a_eblg.getDBL()[mfi];
      IntVectSet ivsBox(valid);
      for (VoFIterator vofit(ivsBox, a_eblg.getEBISL()[mfi].getEBGraph()); vofit.ok(); ++vofit)
      {
        const VolIndex& vof = vofit();
        Real calc = divFCalc[mfi](vof,0);
        Real exac = divFExac[mfi](vof,0);
        Real diff = std::abs(calc-exac);
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
  int divOpTest()
  {
    GridParameters paramsCoar, paramsMedi, paramsFine;
    amrex::Print() << "forcing maxLevel = 0 \n";

    getGridParameters(paramsCoar, 0, true);
    paramsFine = paramsCoar;
    paramsFine.refine(2);

    //and defines it using a geometryservice
    int eekflag = makeGeometry(paramsFine);
    if(eekflag != 0) return eekflag;

    Array<EBLevelGrid> veblgCoar, veblgFine;
    getAllIrregEBLG(veblgCoar, paramsCoar);
    getAllIrregEBLG(veblgFine, paramsFine);

    EBLevelGrid eblgCoar = veblgCoar[0];
    EBLevelGrid eblgFine = veblgFine[0];
    amrex::Print() << "coarse  Grids  = " << eblgCoar.getDBL() << endl;
    amrex::Print() << "finest  Grids  = " << eblgFine.getDBL() << endl;

    int nvar = 1; int nghost = 0;
    EBCellFactory factFine(eblgFine.getEBISL());
    EBCellFactory factCoar(eblgCoar.getEBISL());
    FabArray<EBCellFAB> errorFine(eblgFine.getDBL(), eblgFine.getDM(), nvar, nghost, MFInfo(), factFine);
    FabArray<EBCellFAB> errorCoar(eblgCoar.getDBL(), eblgCoar.getDM(), nvar, nghost, MFInfo(), factCoar);

    getError(errorFine, eblgFine);
    getError(errorCoar, eblgCoar);

    EBLevelDataOps::compareError(errorFine, errorCoar, eblgFine, eblgCoar);
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

  retval = amrex::divOpTest();
  if(retval != 0)
  {
   amrex::Print() << "divergence operator test failed with code " << retval << "\n";
  }
  else
  {
    amrex::Print() << "divergence operator  test passed \n";
  }
  amrex::Finalize();
  return retval;
}

