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
#include "AMReX_AggEBPWLFillPatch.H"

namespace amrex
{
  Real g_coarTimeOld = 0.0;
  Real g_coarTimeNew = 1.0;
  Real g_fineTime = 0.25;
/***************/
  Real exactFunc(const IntVect& a_iv,
                 const Real& a_dx,
                 const Real& a_time)
  {
    Real retval;
    Real probHi;
    ParmParse pp;
    pp.get("prob_hi",probHi);
    RealVect xloc;
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      xloc[idir] = a_dx*(Real(a_iv[idir]) + 0.5);
    }
    Real x = xloc[0]/probHi;
    Real y = xloc[1]/probHi;

    retval = (1.0+ a_time)*(1.0 + x*x + y*y + x*x*x + y*y*y);
    //retval = a_time;
    //  retval = 2.*xloc[0];
    //  retval  = x;
    //  retval = xloc[0]*xloc[0];
    return retval;
  }
/***************/
  int
  makeLayout(BoxArray& a_dbl,
             const Box& a_domainFine)
  {
    //set up mesh refine object
    ParmParse pp;
    int eekflag= 0;
    int maxsize;
    pp.get("maxboxsize",maxsize);
    int bufferSize = 1;
    int blockFactor = 2;
    Real fillrat = 0.75;
    Box domainCoar = coarsen(a_domainFine, 2);
    std::vector<int> refRat(2,2);
    BRMeshRefine meshRefObj(domainCoar, refRat, fillrat,
                            blockFactor, bufferSize, maxsize);

    std::vector<std::vector<Box> > oldMeshes(2);
    oldMeshes[0] = std::vector<Box>(1,   domainCoar);
    oldMeshes[1] = std::vector<Box>(1, a_domainFine);

    Vector<Box>
      //set up coarse tags
      int nc = domainCoar.size(0);
    int nmi = nc/2;//16
    int nqu = nc/4;//8
    int ntf = (nc*3)/4;  //24
    int nte = (nc*3)/8; //12
    int nfe = (nc*5)/8; //20
#if (CH_SPACEDIM ==2)
    Box boxf1(IntVect(0, nqu), IntVect(nmi-1,ntf-1));
    Box boxf2(IntVect(nmi,nte), IntVect(ntf-1,nfe-1));
    Box boxf3(IntVect(nqu,0  ), IntVect(nfe-1,nqu-1));
    Box boxf4(IntVect(nfe,nqu), IntVect(nc -1,nte-1));
#else
    Box boxf1(IntVect(0, nqu,nqu), IntVect(nmi-1,ntf-1,ntf-1));
    Box boxf2(IntVect(nmi,nte,nte), IntVect(ntf-1,nfe-1,nfe-1));
    Box boxf3(IntVect(nqu,0,0  ), IntVect(nfe-1,nqu-1,nqu-1));
    Box boxf4(IntVect(nfe,nqu,nqu), IntVect(nc -1,nte-1,nte-1));
#endif
    IntVectSet tags;
    tags |= boxf1;
    tags |= boxf2;
    tags |= boxf3;
    tags |= boxf4;

    int baseLevel = 0;
    int topLevel = 0;
    std::vector<std::vector<Box> > newMeshes;
    meshRefObj.regrid(newMeshes, tags, baseLevel,
                      topLevel, oldMeshes);

    const std::vector<Box>& vbox = newMeshes[1];
    std::vector<int>  procAssign;
    eekflag = LoadBalance(procAssign,vbox);
    if (eekflag != 0) return eekflag;
    a_dbl.define(vbox, procAssign);
    return eekflag;
  }
/***************/
  int getError(FabArray<EBCellFAB> & a_error,
               Vector<EBLevelGrid> & a_eblg,
               const Real& a_dxLev1)
  {
    int nghost = 2;
    EBLevelGrid eblgFine = a_eblg[1];
    EBLevelGrid eblgCoar = a_eblg[0];
    FabArray<EBCellFAB> phi
    FabArray<EBCellFAB> errorCoar(eblgCoar.getDBL(), eblgCoar.getDM(), nvar, nghost, MFInfo(), factCoar);
    IntVect ghost = IntVect::Unit;
    EBCellFactory factFine(a_eblg[1].getEBISL());
    EBCellFactory factCoar(a_eblg[0].getEBISL());
    FabArray<EBCellFAB> phiFineCalc(eblgFine.getDBL(), eblgFine.getDM(), nvar, nghost, MFInfo(), factFine);
    FabArray<EBCellFAB> phiFineExac(eblgFine.getDBL(), eblgFine.getDM(), nvar, nghost, MFInfo(), factFine);
    FabArray<EBCellFAB> phiCoarOld (eblgCoar.getDBL(), eblgCoar.getDM(), nvar, nghost, MFInfo(), factCoar);
    FabArray<EBCellFAB> phiCoarNew (eblgCoar.getDBL(), eblgCoar.getDM(), nvar, nghost, MFInfo(), factCoar);

    for(MFIter mfi(eblgFine.getDBL(), eblgFine.getDM()); mfi.isValid(); ++mfi)
    {
      Box valid = eblgFine.getDBL()[mfi];
      Box grownBox = grow(valid, nghost);
      grownBox &= a_eblgFine.getDomain();
      IntVectSet ivsBox(grownBox);
      phiFineCalc[mfi].setVal(0.0);
      phiFineExac[mfi].setVal(0.0);
      a_errorFine[mfi].setVal(0.0);

      for (VoFIterator vofit(ivsBox, a_ebislFine[dit()].getEBGraph());
           vofit.ok(); ++vofit)
      {
        const VolIndex& vof = vofit();
        Real rightAns = exactFunc(vof.gridIndex(), a_dxFine, g_fineTime);
        if(valid.contains(vof.gridIndex()))
        {
          phiFineCalc[mfi](vof, 0) = rightAns;
        }
        phiFineExac[mfi](vof, 0) = rightAns;
      }
    }

    for(MFIter mfi(eblgCoar.getDBL(), eblgCoar.getDM()); mfi.isValid(); ++mfi)
    {
      Box valid = eblgCoar.getDBL()[mfi];
      Box grownBox = grow(valid, nghost);
      grownBox &= a_eblgCoar.getDomain();

      IntVectSet ivsBox(grownBox);
      EBCellFAB& phiCoarOldFAB = phiCoarOld[mfi];
      EBCellFAB& phiCoarNewFAB = phiCoarNew[mfi];
      phiCoarOldFAB.setCoveredCellVal(0.0,0);
      phiCoarNewFAB.setCoveredCellVal(0.0,0);
      for (VoFIterator vofit(ivsBox, ebislCoar[dit()].getEBGraph());
           vofit.ok(); ++vofit)
      {
        const VolIndex& vof = vofit();
        Real rightAnsOld = exactFunc(vof.gridIndex(), dxCoar,g_coarTimeOld);
        Real rightAnsNew = exactFunc(vof.gridIndex(), dxCoar,g_coarTimeNew);
        phiCoarOldFAB(vof, 0) = rightAnsOld;
        phiCoarNewFAB(vof, 0) = rightAnsNew;
      }
    }

    //interpolate phiC onto phiF
    AggEBPWLFillPatch interpOp(eblgFine, eblgCoar, nref, nvar, nghost, nghost);


    interpOp.interpolate(phiFineCalc, phiCoarOld, phiCoarNew,
                         g_coarTimeOld, g_coarTimeNew,
                         g_fineTime, 0, 1);

    //error = phiF - phiFExact
    for(MFIter mfi(eblgFine.getDBL(), eblgFine.getDM()); mfi.isValid(); ++mfi)
    {
      Box valid = eblgFine.getDBL()[mfi];
      Box grownBox = grow(valid, nghost);
      grownBox &= a_eblgFine.getDomain();
      IntVectSet ivsBox(grownBox);
      Real     maxDiff = 0;
      VolIndex vofDiff;
      bool found = false;
      for (VoFIterator vofit(ivsBox, a_ebislFine[dit()].getEBGraph());
           vofit.ok(); ++vofit)
      {
        const VolIndex& vof = vofit();
        Real diff = Abs(phiFineCalc[mfi](vof,0)-phFineExac[mfi](vof,0));
        a_error[mfi](vof, 0) = diff;
      }
    }

    return eekflag;
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

    int whichgeom;
    pp.get("which_geom",whichgeom);
    if (whichgeom == 0)
    {
      //allregular
      amrex::Print() << "all regular geometry" << endl;
      AllRegularService regserv;
      EBIndexSpace* ebisPtr = AMReX_EBIS::instance();
      ebisPtr->define(domainFinest, origin, dxFinest, regserv);
    }
    else if (whichgeom == 1)
    {
      amrex::Print() << "ramp geometry" << endl;
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

      RealVect vectDx = RealVect::Unit;
      vectDx *= a_dx;

      GeometryShop workshop(ramp,0);
      //this generates the new EBIS
      EBIndexSpace* ebisPtr = AMReX_EBIS::instance();
      ebisPtr->define(domainFinest, origin, dxFinest, workshop);
    }
    else if (whichgeom == 5)
    {
      std::vector<Real> centervec(SpaceDim);
      std::vector<int>  ncellsvec(SpaceDim);
      int maxgrid;
      ParmParse pp;
      pp.getarr(  "n_cell"       , ncellsvec, 0, SpaceDim);
      pp.get(   "sphere_radius", radius);
      pp.get(   "max_grid_size", maxgrid);
      pp.getarr("sphere_center", centervec, 0, SpaceDim);
      pp.get("domain_length", domlen);                     
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
  void
  ebpwlfpTest()
  {
    int eekflag = 0;
    GridParameters paramsCoar;
    amrex::Print() << "forcing maxLevel = 1 \n";

    getGridParameters(paramsCoar, 1, true);
    paramsFine = paramsCoar;
    paramsFine.refine(2);
    //and defines it using a geometryservice
    eekflag =  makeGeometry(params);

    std::vector<EBLevelGrid> eblgFine, eblgCoar;
    getAllIrregEBLG(eblgFine, paramFine);
    getAllIrregEBLG(eblgCoar, paramCoar);
    
    int nvar = 1;
    EBCellFactory factFine(eblgFine[1].getEBISL());
    EBCellFactory factCoar(eblgCoar[1].getEBISL());
    int nghost = 2;

    FabArray<EBCellFAB> errorFine(eblgFine[1].getDBL(), eblgFine[1].getDM(), nvar, nghost, MFInfo(), factFine);
    FabArray<EBCellFAB> errorCoar(eblgCaor[1].getDBL(), eblgCoar[1].getDM(), nvar, nghost, MFInfo(), factCoar);

    Real dxLev1Fine = a_paramsFine.coarsestDx/a_paramsFine.refRatio[0];
    Real dxLev1Coar = a_paramsCoar.coarsestDx/a_paramsCoar.refRatio[0];
    
    eekflag = getError(errorFine,eblgFine, dxLev1Fine);
    eekflag = getError(errorCoar,eblgCoar, dxLev1Coar);

    EBLevelGrid::compareError(errorFine, errorCoar, eblgFine[1], eblgCoar[1]);
  }
}
/***************/
int
main(int argc, char* argv[])
{
  int retval = 0;
  amrex::Initialize(argc,argv);

  amrex::ebpwlfpConvTest();

  amrex::Finalize();
  return retval;
}
/************/
