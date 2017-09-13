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
#include "AMReX_EBCoarseAverage.H"

namespace amrex
{

/***************/
  int checkConstantAveraging(const EBLevelGrid& a_eblgFine,
                             const EBLevelGrid& a_eblgCoar)
  {
    const Real exactVal = 7.654321;
    int nghost = 0;
    int nvar = 1;
    EBFluxFactory factFine(a_eblgFine.getEBISL());
    EBFluxFactory factCoar(a_eblgCoar.getEBISL());
    FabArray<EBFluxFAB> phiFine(a_eblgFine.getDBL(), a_eblgFine.getDM(), nvar, nghost, MFInfo(), factFine);
    FabArray<EBFluxFAB> phiCoar(a_eblgCoar.getDBL(), a_eblgCoar.getDM(), nvar, nghost, MFInfo(), factCoar);


    for(MFIter mfi(a_eblgFine.getDBL(), a_eblgFine.getDM()); mfi.isValid(); ++mfi)
    {
      phiFine[mfi].setVal(exactVal);
    }

    int nref = 2;
    //interpolate phiC onto phiF
    EBCoarseAverage aveOp(a_eblgFine, a_eblgCoar, nref, nghost, true, true);


    aveOp.average(phiCoar, phiFine, 0, 0, nvar);


    //error = phiF - phiFExact
    for(MFIter mfi(a_eblgCoar.getDBL(), a_eblgCoar.getDM()); mfi.isValid(); ++mfi)
    {
      Box valid = a_eblgCoar.getDBL()[mfi];
      Box grownBox = grow(valid, nghost);
      grownBox &= a_eblgCoar.getDomain();
      IntVectSet ivsBox(grownBox);
      IntVectSet ivsIrreg = a_eblgCoar.getEBISL()[mfi].getIrregIVS(valid);
      for (VoFIterator vofit(ivsBox, a_eblgCoar.getEBISL()[mfi].getEBGraph()); vofit.ok(); ++vofit)
      {
        const VolIndex& vof = vofit();
        //check the irregular value
        if(ivsIrreg.contains(vof.gridIndex()))
        {
          Real irregVal = phiCoar[mfi].getEBFlux()(vof, 0);
          if(std::abs(irregVal - exactVal) > 1.0e-6)
          {
            amrex::Print() << "irregular value off = " << irregVal<< endl;
            return -1;
          }
        }
        //now check the coordinate faces
        for(int idir = 0; idir < SpaceDim; idir++)
        {
          for(SideIterator sit; sit.ok(); ++sit)
          {
            Array<FaceIndex> faces = a_eblgCoar.getEBISL()[mfi].getFaces(vof, idir, sit());
            for(int iface = 0; iface < faces.size(); iface++)
            {
              const EBFaceFAB& facefab = phiCoar[mfi][idir];
              Real faceVal = facefab(faces[iface], 0);
              if(std::abs(faceVal - exactVal) > 1.0e-6)
              {
                amrex::Print() << "coordinate face  value off = " << faceVal<< endl;
                return -2;
              }
            }
          }
        }
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
  int coarseAveTest()
  {
    GridParameters paramsCoar, paramsFine;
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

    eekflag = checkConstantAveraging(eblgFine, eblgCoar);

    return eekflag;
  }
}
/***************/
/***************/
int
main(int argc, char* argv[])
{
  int retval = 0;
  amrex::Initialize(argc,argv);

  retval = amrex::coarseAveTest();
  if(retval != 0)
  {
    amrex::Print() << "averaging test failed with code " << retval << "\n";
  }
  else
  {
    amrex::Print() << "averaging test passed \n";
  }
  amrex::Finalize();
  return retval;
}

