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
#include "AMReX_EBConductivityOp.H"
#include "AMReX_EBConductivityOpFactory.H"
#include "AMReX_DirichletConductivityEBBC.H"
#include "AMReX_DirichletConductivityDomainBC.H"
#include "AMReX_NeumannConductivityEBBC.H"
#include "AMReX_NeumannConductivityDomainBC.H"
#include "AMReX_EBSimpleSolver.H"
#include "AMReX_LinearSolver.H"

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
  void setPhi(FabArray<EBCellFAB>& a_phi,
              const EBLevelGrid  & a_eblg,
              const Real         & a_dx)

  {
    static const Real pi = 4.0*atan(1.0);
    for(MFIter mfi(a_eblg.getDBL(), a_eblg.getDM()); mfi.isValid(); ++mfi)
    {
      const Box    & grid = a_eblg.getDBL()[mfi];
      const EBISBox& ebis = a_eblg.getEBISL()[mfi];
      a_phi[mfi].setVal(0.);
      IntVectSet ivs(grid);
      for(VoFIterator vofit(ivs, ebis.getEBGraph()); vofit.ok(); ++vofit)
      {
        const IntVect& iv = vofit().gridIndex();
        Real x = a_dx*(iv[0]+0.5);
        Real y = a_dx*(iv[1]+0.5);
        Real val = sin(pi*x)*sin(pi*y);
        a_phi[mfi](vofit(), 0) = val;
      }
    }
  }
//----------------
  void defineConductivitySolver(AMREBMultiGrid<FabArray<EBCellFAB> >& a_solver,
                                const vector<EBLevelGrid>         & a_veblg,
                                const GridParameters              & a_params,
                                int a_nghost)
  {

    bool enableLevelSolves = true; 
    Real alpha, beta, acoefval, bcoefval;
    ParmParse pp;
    pp.get("alpha", alpha);
    pp.get("beta" , beta );
    pp.get("acoef_val", acoefval);
    pp.get("bcoef_val" ,bcoefval);
    int domain_bc_type, eb_bc_type;
    pp.get("domain_bc_type" ,domain_bc_type);
    pp.get(    "eb_bc_type" ,    eb_bc_type);
    shared_ptr<ConductivityBaseDomainBCFactory> domainBCFact;
    shared_ptr<ConductivityBaseEBBCFactory>     ebBCFact;
    if(domain_bc_type == 0)
    {
      pout() << "neumann domain bcs" << endl;
      NeumannConductivityDomainBCFactory* neumptr = (new NeumannConductivityDomainBCFactory());
      neumptr->setValue(0.);
      domainBCFact = shared_ptr<ConductivityBaseDomainBCFactory>(static_cast<ConductivityBaseDomainBCFactory*>(neumptr));
    }
    else if(domain_bc_type == 1)
    {
      pout() << "dirichelet domain bcs" << endl;
      DirichletConductivityDomainBCFactory* diriptr = (new DirichletConductivityDomainBCFactory());
      diriptr->setValue(0.);
      domainBCFact = shared_ptr<ConductivityBaseDomainBCFactory>(static_cast<ConductivityBaseDomainBCFactory*>(diriptr));
    }
    else
    {
      amrex::Error("unknown domain_bc_type");
    }

    if(eb_bc_type == 0)
    {
      pout() << "neumann eb bcs" << endl;
      NeumannConductivityEBBCFactory* neumptr = (new NeumannConductivityEBBCFactory());
      neumptr->setValue(0.);
      ebBCFact = shared_ptr<ConductivityBaseEBBCFactory>(static_cast<ConductivityBaseEBBCFactory*>(neumptr));
    }
    else if(eb_bc_type == 1)
    {
      int order_ebbc;
      pp.get("order_ebbc", order_ebbc);
      
      pout() << "dirichelet eb bcs with order " << order_ebbc << endl;
      DirichletConductivityEBBCFactory* diriptr = (new DirichletConductivityEBBCFactory());
      diriptr->setValue(0.);
      diriptr->setOrder(order_ebbc);
      ebBCFact = shared_ptr<ConductivityBaseEBBCFactory>(static_cast<ConductivityBaseEBBCFactory*>(diriptr));
    }
    else
    {
      amrex::Error("unknown eb_bc_type");
    }


    int numLevels = a_veblg.size();
    vector<shared_ptr<FabArray<EBCellFAB> > > acoef(numLevels);
    vector<shared_ptr<FabArray<EBFluxFAB> > > bcoef(numLevels);
    for(int ilev = 0; ilev < numLevels; ilev++)
    {
      const BoxArray           & ba = a_veblg[ilev].getDBL();
      const DistributionMapping& dm = a_veblg[ilev].getDM();
      EBCellFactory cellfact(a_veblg[ilev].getEBISL());
      EBFluxFactory fluxfact(a_veblg[ilev].getEBISL());

      acoef[ilev] = shared_ptr<FabArray<EBCellFAB> >(new FabArray<EBCellFAB>(ba, dm, 1, a_nghost, MFInfo(), cellfact));
      bcoef[ilev] = shared_ptr<FabArray<EBFluxFAB> >(new FabArray<EBFluxFAB>(ba, dm, 1, a_nghost, MFInfo(), fluxfact));
      
      EBLevelDataOps::setVal(*acoef[ilev], acoefval);
      EBLevelDataOps::setVal(*bcoef[ilev], bcoefval);
    }
    
    EBConductivityOpFactory factory(a_veblg, alpha, beta, acoef, bcoef, a_params.coarsestDx, a_params.refRatio,
                                    domainBCFact, ebBCFact, a_nghost);

    EBSimpleSolver* simple = new EBSimpleSolver();
    int num_smooth_bottom = 10;
    pp.query("num_smooth_bottom", num_smooth_bottom);
    simple->setNumSmooths(num_smooth_bottom);
    shared_ptr<LinearSolver<FabArray<EBCellFAB> > > bottomSolverPtr(static_cast<LinearSolver<FabArray<EBCellFAB> > *>(simple));

    AMRLevelOpFactory<FabArray<EBCellFAB> >* factptr = static_cast<AMRLevelOpFactory<FabArray<EBCellFAB> >* >(&factory);
    a_solver.define(a_params.coarsestDomain,*factptr, bottomSolverPtr, a_params.refRatio, enableLevelSolves, numLevels);

    int numSmooth, numMG, maxIter;
    Real eps, hang;
    pp.get("num_smooth", numSmooth);
    pp.get("num_mg",     numMG);
    pp.get("max_iterations", maxIter);
    pp.get("tolerance", eps);
    pp.get("hang",      hang);
    Real normThresh = 1.0e-30;
    a_solver.setSolverParameters(numSmooth, numSmooth, numSmooth,
                                 numMG, maxIter, eps, hang, normThresh);
    a_solver.m_verbosity = 5;

  }
  //----------------
  int conjuctionJunction()
  {
    GridParameters params;

    getGridParameters(params);
    int eekflag = makeGeometry(params);
    if(eekflag != 0) return eekflag;

    std::vector<EBLevelGrid> veblg;
    getAllIrregEBLG(veblg, params);
    for(int ilev = 0; ilev < veblg.size(); ilev ++)
    {
      amrex::Print() << "grids[" << ilev << "] = " << veblg[ilev].getDBL() << endl;
    }

    int nvar = 1; int nghost = 4;
    vector<FabArray<EBCellFAB>* > phi(veblg.size(), NULL);
    vector<FabArray<EBCellFAB>* > rhs(veblg.size(), NULL);
    Real dxLev = params.coarsestDx;
    for(int ilev = 0; ilev < veblg.size(); ilev++)
    {
      EBCellFactory fact(veblg[ilev].getEBISL());
      phi[ilev] = new FabArray<EBCellFAB>(veblg[ilev].getDBL(), veblg[ilev].getDM(), nvar, nghost, MFInfo(), fact);
      rhs[ilev] = new FabArray<EBCellFAB>(veblg[ilev].getDBL(), veblg[ilev].getDM(), nvar, nghost, MFInfo(), fact);
      setPhi(*phi[ilev], veblg[ilev], dxLev);
      EBLevelDataOps::setVal(*rhs[ilev], 1.0);
      dxLev /= params.refRatio[ilev];
    }

    AMREBMultiGrid<FabArray<EBCellFAB> > solver;
    defineConductivitySolver(solver, veblg, params, nghost);
    int lbase = 0; int lmax = params.maxLevel;

    bool zeroPhi = false;
    solver.solve(phi, rhs, lmax, lbase,  zeroPhi);
      
    EBLevelDataOps::writeEBAMRPlotFile(string("phi.ebplt"), phi, veblg, params.refRatio, vector<string>(1, "phi"));
    EBLevelDataOps::writeEBAMRPlotFile(string("rhs.ebplt"), rhs, veblg, params.refRatio, vector<string>(1, "rhs"));

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

