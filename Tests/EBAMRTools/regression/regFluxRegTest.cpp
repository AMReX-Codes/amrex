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
#include "AMReX_MultiFab.H"

namespace amrex
{
  int fluxRegTest()
  {
    Box bCoar(IntVect::Zero, 31*IntVect::Unit);
    Box bFine(IntVect::Zero,  7*IntVect::Unit);
    Box  domCoar = bCoar;
    //Box  domFine = refine(bCoar, 2);

    int refRat = 2;
    BoxArray baCoar(bCoar);
    BoxArray baFine(bFine);
    DistributionMapping dmCoar(baCoar);
    DistributionMapping dmFine(baFine);
    MultiFab datCoar(baCoar, dmCoar, 1, 1, MFInfo());
    for(MFIter mfi(datCoar); mfi.isValid(); ++mfi)
    {
      datCoar[mfi].setVal(0.);
    }
                     
    FluxRegister fluxReg(baFine, dmFine, refRat*IntVect::Unit, 1, 1);
    Real dxCoar = 1.0/bCoar.size()[0];
    Real scale = 1.0/dxCoar;
    Real fluxVal = 4.77;
    //all that matters is that the ratio is correct
    Real coarArea = AMREX_D_TERM(1.0, *refRat, *refRat);
    Real fineArea = 1.0;
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      IntVect indexType = IntVect::Zero;
      indexType[idir] = 1;
      BoxArray baFaceCoar = convert(baCoar, indexType);
      BoxArray baFaceFine = convert(baFine, indexType);
      MultiFab fluxCoar(baFaceCoar, dmCoar, 1, 1, MFInfo());
      MultiFab fluxFine(baFaceFine, dmCoar, 1, 1, MFInfo());
      for(MFIter mfi(fluxCoar); mfi.isValid(); ++mfi)
      {
        fluxCoar[mfi].setVal(fluxVal);
      }
      for(MFIter mfi(fluxFine); mfi.isValid(); ++mfi)
      {
        fluxFine[mfi].setVal(fluxVal);
      }
      //boxlib's flux register does not scale finefaces/coarseface for you.
      fluxReg.CrseInit(fluxCoar, idir, 0, 0, 1, -coarArea);  //the -1 is default.
      fluxReg.FineAdd (fluxFine, idir, 0, 0, 1,  fineArea);
    }
    Real ratio = coarArea/fineArea;
    amrex::Print() << "ratio ="  << ratio << endl;
    RealVect rvlo = RealVect::Zero;
    RealVect rvhi = RealVect::Unit;
    RealBox rb(rvlo.dataPtr(), rvhi.dataPtr());
    //Geometry geomFine(domFine, &rb);
    Geometry geomCoar(domCoar, &rb);
    fluxReg.ClearInternalBorders(geomCoar);
    fluxReg.Reflux(datCoar, scale, 0, 0, 1, geomCoar);
    Real tol = 1.0e-10;
    for(MFIter mfi(datCoar); mfi.isValid(); ++mfi)
    {
      Box     grid = baCoar[mfi];
      for(BoxIterator boxit(grid); boxit.ok(); ++boxit)
      {
        Real dataval = datCoar[mfi](boxit(), 0);
        if(std::abs(dataval) > tol)
        {
          amrex::Print() << "reflux failed at cell " << boxit() << endl;
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
      amrex::Print() << "reflux test failed with code " << retval << "\n";
    }
    else
    {
      amrex::Print() << "reflux test passed \n";
    }
    amrex::Finalize();
    return retval;
  }

