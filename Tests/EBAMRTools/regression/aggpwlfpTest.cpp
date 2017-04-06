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
#include "AMReX_LevelData.H"
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
  int compareError(const FabArray<EBCellFAB>& a_errorFine,
                   const EBISLayout& a_ebislFine,
                   const BoxArray& a_gridsFine,
                   const Box& a_domainFine,
                   const FabArray<EBCellFAB>& a_errorCoar,
                   const EBISLayout& a_ebislCoar,
                   const BoxArray& a_gridsCoar,
                   const Box& a_domainCoar)
  {
    //the norms here must go over ghost cells
    //too or the test is really silly
    BoxLayout blNormFine, blNormCoar;
    blNormFine.deepCopy(a_gridsFine);
    blNormCoar.deepCopy(a_gridsCoar);
    blNormFine.grow(1);
    blNormCoar.grow(1);
    blNormFine&= a_domainFine;
    blNormCoar&= a_domainCoar;
    blNormFine.close();
    blNormCoar.close();

    EBCellFactory factFine(a_ebislFine);
    EBCellFactory factCoar(a_ebislCoar);
    int eekflag = 0;
    for (int itype = 0; itype < 3; itype++)
    {
      EBNormType::NormMode normtype;
      if (itype == 0)
      {
        normtype = EBNormType::OverBoth;
        amrex::Print() << endl << " Using all uncovered cells." << endl  << endl;
      }
      else if (itype == 1)
      {
        normtype = EBNormType::OverOnlyRegular;
        amrex::Print() << endl << " Using only regular cells." << endl << endl;
      }
      else
      {
        normtype = EBNormType::OverOnlyIrregular;
        amrex::Print() << endl << " Using only irregular cells." << endl << endl;
      }
      int comp = 0;
      for (int inorm = 0; inorm < 3; inorm++)
      {

        if (inorm == 0)
        {
          amrex::Print() << endl << " Using max norm." << endl;
        }
        else
        {
          amrex::Print() << endl << " Using L-" << inorm << "norm." << endl;
        }
        Real coarnorm = EBArith::norm(a_errorCoar,
                                      blNormCoar, a_ebislCoar,
                                      comp, inorm, normtype);
        Real finenorm = EBArith::norm(a_errorFine,
                                      blNormFine, a_ebislFine,
                                      comp, inorm, normtype);
        amrex::Print() << "Coarse Error Norm = " << coarnorm << endl;
        amrex::Print() << "Fine   Error Norm = " << finenorm << endl;
        if (Abs(finenorm) > 1.0e-8)
        {
          Real order = log(coarnorm/finenorm)/log(2.0);
          amrex::Print() << "Order of scheme = " << order << endl;
        }
      }
    }

    return eekflag;
  }
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
  Vector<int> refRat(2,2);
  BRMeshRefine meshRefObj(domainCoar, refRat, fillrat,
                          blockFactor, bufferSize, maxsize);

  Vector<Vector<Box> > oldMeshes(2);
  oldMeshes[0] = Vector<Box>(1,   domainCoar);
  oldMeshes[1] = Vector<Box>(1, a_domainFine);

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
  Vector<Vector<Box> > newMeshes;
  meshRefObj.regrid(newMeshes, tags, baseLevel,
                    topLevel, oldMeshes);

  const Vector<Box>& vbox = newMeshes[1];
  Vector<int>  procAssign;
  eekflag = LoadBalance(procAssign,vbox);
  if (eekflag != 0) return eekflag;
  a_dbl.define(vbox, procAssign);
  return eekflag;
}
/***************/
int getError(FabArray<EBCellFAB>& a_errorFine,
             const EBISLayout& a_ebislFine,
             const BoxArray& a_gridsFine,
             const Box& a_domainFine,
             const Real& a_dxFine)
{
  ParmParse pp;
  int eekflag = 0;
  int nvar = 1;
  int nghost = 2;
  int nref = 2;
  int maxsize;
  pp.get("maxboxsize",maxsize);

  //generate coarse dx, domain
  Box domainCoar = coarsen(a_domainFine, nref);

  Real dxCoar = nref*a_dxFine;
  //make the coarse grids completely cover the domain
  Vector<Box> vbox;
  domainSplit(domainCoar, vbox,  maxsize);
  Vector<int>  procAssign;
  eekflag = LoadBalance(procAssign,vbox);
  BoxArray gridsCoar(vbox, procAssign);
  //make coarse ebisl
  EBISLayout ebislCoar;
  eekflag = makeEBISL(ebislCoar, gridsCoar, domainCoar, nghost);
  if (eekflag != 0) return eekflag;
  // amrex::Print() << "getError dom coar = " << domainCoar << endl;;
  // amrex::Print() << "getError grids coar = " << gridsCoar << endl;;
  // amrex::Print() << "getError dom fine = " << a_domainFine << endl;;
  // amrex::Print() << "getError grids fine = " << a_gridsFine << endl;;

  //create data at both refinemenets
  IntVect ghost = IntVect::Unit;
  EBCellFactory ebcellfactFine(a_ebislFine);
  EBCellFactory ebcellfactCoar(ebislCoar);
  FabArray<EBCellFAB> phiFine(a_gridsFine,  nvar, ghost, ebcellfactFine);
  FabArray<EBCellFAB> oldFine(a_gridsFine,  nvar, ghost, ebcellfactFine);
  FabArray<EBCellFAB> phiCoarOld(gridsCoar, nvar, ghost, ebcellfactCoar);
  FabArray<EBCellFAB> phiCoarNew(gridsCoar, nvar, ghost, ebcellfactCoar);

  //fill phi fine and phiCExact
  //put phiFexact into a_errorFine and phiFine
  //this should make the error at the interior = 0
  for (DataIterator dit = a_gridsFine.dataIterator();
      dit.ok(); ++dit)
    {
      Box grownBox = grow(a_gridsFine.get(dit()), 1);
      grownBox &= a_domainFine;
      IntVectSet ivsBox(grownBox);
      phiFine[dit()].setVal(0.0);
      oldFine[dit()].setVal(0.0);
      a_errorFine[dit()].setVal(0.0);
      EBCellFAB& phiFineFAB = phiFine[dit()];
      EBCellFAB& oldFineFAB = oldFine[dit()];
      EBCellFAB& errFineFAB = a_errorFine[dit()];
      phiFineFAB.setCoveredCellVal(0.0,0);
      oldFineFAB.setCoveredCellVal(0.0,0);
      for (VoFIterator vofit(ivsBox, a_ebislFine[dit()].getEBGraph());
          vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          Real rightAns = exactFunc(vof.gridIndex(), a_dxFine, g_fineTime);
          phiFineFAB(vof, 0) = rightAns;
          oldFineFAB(vof, 0) = rightAns;
          errFineFAB(vof, 0) = rightAns;
        }
    }

  for (DataIterator dit = gridsCoar.dataIterator();
      dit.ok(); ++dit)
    {
      Box grownBox = grow(gridsCoar.get(dit()), 1);
      grownBox &= domainCoar;
      IntVectSet ivsBox(grownBox);
      EBCellFAB& phiCoarOldFAB = phiCoarOld[dit()];
      EBCellFAB& phiCoarNewFAB = phiCoarNew[dit()];
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
  AggEBPWLFillPatch interpOp(a_gridsFine, gridsCoar,
                             a_ebislFine, ebislCoar,
                             domainCoar, nref, nvar,
                             1, ghost);

  //interpolate phiC onto phiF
  EBPWLFillPatch oldinterpOp(a_gridsFine, gridsCoar,
                             a_ebislFine, ebislCoar,
                             domainCoar, nref, nvar, 1);


  Interval zeroiv(0,0);
  interpOp.interpolate(phiFine, phiCoarOld, phiCoarNew,
                       g_coarTimeOld, g_coarTimeNew,
                       g_fineTime, zeroiv);

  oldinterpOp.interpolate(oldFine, phiCoarOld, phiCoarNew,
                          g_coarTimeOld, g_coarTimeNew,
                          g_fineTime, zeroiv);
  //error = phiF - phiFExact
  int ibox = 0;
  for (DataIterator dit = a_gridsFine.dataIterator();
      dit.ok(); ++dit, ++ibox)
    {
      Box grownBox = grow(a_gridsFine.get(dit()), 1);
      grownBox &= a_domainFine;
      IntVectSet ivsBox(grownBox);

      EBCellFAB& phiFineFAB = phiFine[dit()];
      EBCellFAB& oldFineFAB = oldFine[dit()];
      Real     maxDiff = 0;
      VolIndex vofDiff;
      bool found = false;
      for (VoFIterator vofit(ivsBox, a_ebislFine[dit()].getEBGraph());
          vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          Real diff = Abs(phiFineFAB(vof,0)-oldFineFAB(vof,0));
          if (diff > maxDiff)
            {
              found = true;
              vofDiff  = vof;
              maxDiff = diff;
            }
        }
      if (found)
        {
          amrex::Print() << "max diff   = " << maxDiff << endl;
          amrex::Print() << "at intvect = " << vofDiff.gridIndex() << endl;
          amrex::Print() << "at box num = " << ibox << endl;
        }

      EBCellFAB& errorFAB = a_errorFine[dit()];
      errorFAB -= phiFineFAB;
    }

  return eekflag;
}
/**********/
int makeGeometry(Box& a_domain,
                 Real& a_dx)
{
  int eekflag =  0;
  //parse input file
  ParmParse pp;
  RealVect origin = RealVect::Zero;
  Vector<int> n_cell(SpaceDim);
  pp.getarr("n_cell",n_cell,0,SpaceDim);

  CH_assert(n_cell.size() == SpaceDim);
  IntVect lo = IntVect::Zero;
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

  Vector<Real> prob_lo(SpaceDim, 1.0);
  Real prob_hi;
  pp.getarr("prob_lo",prob_lo,0,SpaceDim);
  pp.get("prob_hi",prob_hi);
  a_dx = (prob_hi-prob_lo[0])/n_cell[0];
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      origin[idir] = prob_lo[idir];
    }
  int whichgeom;
  pp.get("which_geom",whichgeom);
  if (whichgeom == 0)
    {
      //allregular
      amrex::Print() << "all regular geometry" << endl;
      AllRegularService regserv;
      EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
      ebisPtr->define(a_domain, origin, a_dx, regserv);
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

      GeometryShop workshop(ramp,0,vectDx);
      //this generates the new EBIS
      EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
      ebisPtr->define(a_domain, origin, a_dx, workshop);
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
  Box domainFine, domainCoar;
  Real dxFine, dxCoar;

    //and defines it using a geometryservice
    eekflag =  makeGeometry(domainFine,  dxFine);
    CH_assert(eekflag == 0);

    domainCoar = coarsen(domainFine, 2);
    dxCoar = 2.0*dxFine;
    //make grids
    BoxArray gridsFine, gridsCoar;
    eekflag = makeLayout(gridsFine, domainFine);
    CH_assert(eekflag == 0);
    coarsen(gridsCoar, gridsFine, 2);

    ///create ebislayout
    int nghost = 2;
    EBISLayout ebislFine, ebislCoar;
    eekflag = makeEBISL(ebislFine, gridsFine, domainFine, nghost);
    if (eekflag != 0) return eekflag;
    eekflag = makeEBISL(ebislCoar, gridsCoar, domainCoar, nghost);
    if (eekflag != 0) return eekflag;

    int nvar = 1;
    EBCellFactory ebcellfactFine(ebislFine);
    EBCellFactory ebcellfactCoar(ebislCoar);
    IntVect ivghost = nghost*IntVect::Unit;
    FabArray<EBCellFAB> errorFine(gridsFine, nvar,
                                   ivghost,ebcellfactFine);
    FabArray<EBCellFAB> errorCoar(gridsCoar, nvar,
                                   ivghost,ebcellfactCoar);

    eekflag = getError(errorFine,ebislFine, gridsFine, domainFine, dxFine);
    if (eekflag != 0) return eekflag;

    eekflag = getError(errorCoar,ebislCoar, gridsCoar, domainCoar, dxCoar);
    if (eekflag != 0) return eekflag;

    eekflag =  compareError(errorFine, ebislFine, gridsFine, domainFine,
                            errorCoar, ebislCoar, gridsCoar, domainCoar);
    if (eekflag != 0) return eekflag;

  }//end omnipresent scoping trick
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  ebisPtr->clear();
  if (eekflag == 0)
    {
      amrex::Print() << "aggpwlfpTest passed" << endl;
    }

#ifdef CH_MPI
  MPI_Finalize();
#endif
  return eekflag;
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
