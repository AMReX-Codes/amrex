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
#include "AMReX_EBConducivityOp.H"
#include "AMReX_EBConducivityOpFactory.H"

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
      std::vector<Real> centervec(SpaceDim);
      std::vector<int>  ncellsvec(SpaceDim);
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
  //----------------
  void defineConductivitySolver(AMRMultiGrid<FabArray<EBCellFAB> >& a_solver,
                                const vector<EBLevelGrid>         & a_veblg,
                                const GridParameters              & a_params)
  {
  }
  //----------------
  int conjuctionJunction()
  {
    GridParameters params;

    getGridParameters(params);
    std::vector<EBLevelGrid> eblg;
    getAllIrregEBLG(eblg, params);
    for(int ilev = 0; ilev < 2; ilev ++)
    {
      amrex::Print() << "grids[" << ilev << "] = " << eblg[ilev].getDBL() << endl;
    }

    int eekflag = makeGeometry(params);
    if(eekflag != 0) return eekflag;

    std::vector<EBLevelGrid> veblg;
    getAllIrregEBLG(veblg, params);

    int nvar = 1; int nghost = 4;
    vector<FabArray<EBCellFAB>* > phi(veblg.size(), NULL);
    vector<FabArray<EBCellFAB>* > rhs(veblg.size(), NULL);
    for(int ilev = 0; ilev < veblg.size(); ilev++)
    {
      EBCellFactory fact(veblg[ilev].getEBISL());
      phi[ilev] = new FabArray<EBCellFAB>(veblg[ilev].getDBL(), veblg[ilev].getDM(), nvar, nghost, MFInfo(), fact);
      rhs[ilev] = new FabArray<EBCellFAB>(veblg[ilev].getDBL(), veblg[ilev].getDM(), nvar, nghost, MFInfo(), fact);
      EBLevelDataOps::setVal(*phi[ilev], 0.0);
      EBLevelDataOps::setVal(*rhs[ilev], 1.0);
    }

    AMRMultiGrid<FabArray<EBCellFAB> > solver;
    defineConductivitySolver(solver, veblg, params);
    int lbase = 0; int lmax = params.maxLevel;

    solver.solve(phi, rhs, lbase, lmax);
      
    writeEBPlotFile(string("phi.ebplt"), phi);

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

  retval = amrex::conjuctionJunction();
  if(retval != 0)
  {
    amrex::Print() << "simple conductivity exec test failed with code " << retval << "\n";
  }
  else
  {
    amrex::Print() << "simple conductivity test passed \n";
  }
  amrex::Finalize();
  return retval;
}

