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

#include "AMReX_BaseIVFactory.H"
#include "AMReX_EBIndexSpace.H"
#include "AMReX_EBISLayout.H"
#include "AMReX_BoxIterator.H"
#include "AMReX_ParmParse.H"
#include "AMReX_GeometryShop.H"
#include "AMReX_EBCellFAB.H"
#include "AMReX_EBLevelGrid.H"
#include "AMReX_LayoutData.H"
#include "AMReX_EBCellFactory.H"
#include "AMReX_VoFIterator.H"
#include "AMReX_EBArith.H"
#include "AMReX_AllRegularService.H"
#include "AMReX_PlaneIF.H"
#include "AMReX_SPMD.H"
#include "AMReX_Print.H"
#include "AMReX_EBFluxFactory.H"
#include "AMReX_EBFluxFAB.H"
#include "AMReX_EBCellFactory.H"
#include "AMReX_EBCellFAB.H"
#include "AMReX_IrregFABFactory.H"
#include "AMReX_BaseEBCellFactory.H"
#include "AMReX_IrregFAB.H"
#include "AMReX_EBDataVarMacros.H"
#include "AMReX_FabArrayIO.H"

namespace amrex
{
/***************/
  int makeGeometry(Box& a_domain,
                   Real& a_dx)
  {
    int eekflag =  0;
    //parse input file
    ParmParse pp;
    RealVect origin = RealVect::Zero;
    std::vector<int> n_cell;
    pp.getarr("n_cell", n_cell, 0, SpaceDim);

    IntVect lo = IntVect::TheZeroVector();
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

    Real prob_hi;
    pp.get("prob_hi",prob_hi);
    a_dx = prob_hi/n_cell[0];

    int whichgeom;
    pp.get("which_geom",whichgeom);
    if (whichgeom == 0)
    {
      //allregular
      amrex::Print() << "all regular geometry" << "\n";
      AllRegularService regserv;
      EBIndexSpace* ebisPtr = AMReX_EBIS::instance();
      ebisPtr->define(a_domain, origin, a_dx, regserv);
    }
    else if (whichgeom == 1)
    {
      amrex::Print() << "ramp geometry" << "\n";
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


      GeometryShop workshop(ramp,0, a_dx);
      //this generates the new EBIS
      EBIndexSpace* ebisPtr = AMReX_EBIS::instance();
      ebisPtr->define(a_domain, origin, a_dx, workshop);
    }
    else
    {
      //bogus which_geom
      amrex::Print() << " bogus which_geom input = "
                     << whichgeom << "\n";
      eekflag = 33;
    }

    return eekflag;
  }
  /***************/
  void fillEBCFWithSomething(FabArray<BaseEBCellFAB<int> >& a_ebcf,
                             const EBLevelGrid            & a_eblg)
  {
    for(MFIter mfi(a_eblg.getDBL(), a_eblg.getDM()); mfi.isValid(); ++mfi)
    {
      Box     grid = a_eblg.getDBL()  [mfi];
      EBISBox ebis = a_eblg.getEBISL()[mfi];
      IntVectSet ivs(grid);
      a_ebcf[mfi].setVal(-1);
      int ival = 0;
      for(VoFIterator vofit(ivs, ebis.getEBGraph()); vofit.ok(); ++vofit)
      {
        for(int icomp = 0; icomp < a_ebcf.nComp(); icomp++)
        {
          a_ebcf[mfi](vofit(), icomp) = ival;
          ival++;
        }
      }
    }
  }
  /****/
  int checkEquality(const FabArray<BaseEBCellFAB<int> > & a_ebcf1,
                    const FabArray<BaseEBCellFAB<int> > & a_ebcf2,
                    const EBLevelGrid                   & a_eblg)
  {
    for(MFIter mfi(a_eblg.getDBL(), a_eblg.getDM()); mfi.isValid(); ++mfi)
    {
      Box     grid = a_eblg.getDBL()  [mfi];
      EBISBox ebis = a_eblg.getEBISL()[mfi];
      IntVectSet ivs(grid);
      for(VoFIterator vofit(ivs, ebis.getEBGraph()); vofit.ok(); ++vofit)
      {
        for(int icomp = 0; icomp < a_ebcf1.nComp(); icomp++)
        {
          int val1 = a_ebcf1[mfi](vofit(), icomp);
          int val2 = a_ebcf2[mfi](vofit(), icomp);
          if(val1 != val2)
          {
            amrex::Print() << "ebcf values do not match at " << vofit().gridIndex();
            return -1;
          }
        }
      }
    }
    return 0;
  }

  /****/
  int checkGraphs(  const FabArray<EBGraph > & a_ebg1,
                    const FabArray<EBGraph > & a_ebg2,
                    const EBLevelGrid        & a_eblg)
  {
    
    for(MFIter mfi(a_eblg.getDBL(), a_eblg.getDM()); mfi.isValid(); ++mfi)
    {
      EBGraph ebg1 = a_ebg1[mfi];
      EBGraph ebg2 = a_ebg2[mfi];
      if(ebg1.getDomain() != ebg2.getDomain())
      {
        amrex::Print() << "checkgraph: domain mismatch" << endl;
        return -1;
      }
      if(ebg1.getRegion() != ebg2.getRegion())
      {
        amrex::Print() << "checkgraph: region mismatch" << endl;
        return -2;
      }
      const Box& grid = a_eblg.getDBL()[mfi];
      IntVectSet ivs(grid);
      VoFIterator vofit1(ivs, ebg1);
      VoFIterator vofit2(ivs, ebg2);
      vector<VolIndex> vvof1 = vofit1.getVector();
      vector<VolIndex> vvof2 = vofit2.getVector();
      if(vvof1.size() != vvof2.size())
      {
        amrex::Print() << "checkgraph: vector vof size mismatch" << endl;
        return -3;
      }
      for(int ivec = 0; ivec < vvof1.size(); ivec++)
      {
        if(vvof1[ivec] != vvof2[ivec])
        {
          amrex::Print() << "checkgraph: vof mismatch at ivec = "<< ivec << endl;
          return -4;
        }
      }
      for(int facedir = 0; facedir < SpaceDim; facedir++)
      {
        FaceIterator faceit1(ivs, ebg1, facedir, FaceStop::SurroundingWithBoundary);
        FaceIterator faceit2(ivs, ebg2, facedir, FaceStop::SurroundingWithBoundary);
        vector<FaceIndex> vfac1 = faceit1.getVector();
        vector<FaceIndex> vfac2 = faceit2.getVector();
        if(vfac1.size() != vfac2.size())
        {
          amrex::Print() << "checkgraph: vector face size mismatch" << endl;
          return -5;
        }
        for(int ivec = 0; ivec < vfac1.size(); ivec++)
        {
          if(vfac1[ivec] != vfac2[ivec])
          {
            amrex::Print() << "checkgraph: face mismatch at ivec = "<< ivec << endl;
            return -6;
          }
        }
      }
    }
    return 0;
  }
  /****/
  int checkData(  const FabArray<EBData> & a_ebd1,
                  const FabArray<EBData> & a_ebd2,
                  const EBLevelGrid      & a_eblg)
  {
    
    for(MFIter mfi(a_eblg.getDBL(), a_eblg.getDM()); mfi.isValid(); ++mfi)
    {
      EBData ebd1 = a_ebd1[mfi];
      EBData ebd2 = a_ebd2[mfi];
      if(ebd1.getRegion() != ebd2.getRegion())
      {
        amrex::Print() << "checkdata: region mismatch" << endl;
        return -2;
      }
      //check the volume data
      BaseIVFAB<Real>& vdata1 = ebd1.getVolData();
      BaseIVFAB<Real>& vdata2 = ebd2.getVolData();
      if((vdata1.nComp() != V_VOLNUMBER) || (vdata2.nComp() != V_VOLNUMBER))
      {
        amrex::Print() << "checkdata: vdata comps wrong" << endl;
        return -5;
      }
      const vector<VolIndex>& vvof1 = vdata1.getVoFs();
      const vector<VolIndex>& vvof2 = vdata2.getVoFs();
      if(vvof1.size() != vvof2.size())
      {
        amrex::Print() << "checkdata: vector vof size mismatch" << endl;
        return -3;
      }
      for(int ivec = 0; ivec < vvof1.size(); ivec++)
      {
        if(vvof1[ivec] != vvof2[ivec])
        {
          amrex::Print() << "checkvof: vof mismatch at ivec = "<< ivec << endl;
          return -4;
        }
        for(int icomp = 0; icomp <  V_VOLNUMBER; icomp++)
        {
          Real val1 = vdata1(vvof1[ivec], icomp);
          Real val2 = vdata2(vvof2[ivec], icomp);
          Real tol = 1.0e-9;
          if(std::abs(val1-val2) > tol)
          {
            amrex::Print() << "checkvof: value mismatch at ivec, ivar = "<< ivec  << ","<< icomp<< endl;
            return -6;
          }
        }
      }
      //check the face data
      for(int idir = 0; idir < SpaceDim; idir++)
      {
        BaseIFFAB<Real>& fdata1 = ebd1.getFaceData(idir);
        BaseIFFAB<Real>& fdata2 = ebd2.getFaceData(idir);
        if((fdata1.nComp() != F_FACENUMBER) || (fdata2.nComp() != F_FACENUMBER))
        {
          amrex::Print() << "checkdata: fdata comps wrong" << endl;
          return -7;
        }
        const vector<FaceIndex>& vfac1 = fdata1.getFaces();
        const vector<FaceIndex>& vfac2 = fdata2.getFaces();
        if(vfac1.size() != vfac2.size())
        {
          amrex::Print() << "checkdata: vector face size mismatch" << endl;
          return -8;
        }
        for(int ivec = 0; ivec < vfac1.size(); ivec++)
        {
          if(vfac1[ivec] != vfac2[ivec])
          {
            amrex::Print() << "checkvof: vof mismatch at ivec = "<< ivec << endl;
            return -9;
          }
          for(int icomp = 0; icomp <  F_FACENUMBER; icomp++)
          {
            Real val1 = fdata1(vfac1[ivec], icomp);
            Real val2 = fdata2(vfac2[ivec], icomp);
            Real tol = 1.0e-9;
            if(std::abs(val1-val2) > tol)
            {
              amrex::Print() << "checkvof: value mismatch at ivec, ivar = "<< ivec  << ","<< icomp<< endl;
              return -10;
            }
          }
        }
      }
    }        

    return 0;
  }
  /***************/
  int testIO()
  {
    Box domain;
    Real dx;
    makeGeometry(domain, dx);
    int maxboxsize;
    ParmParse pp;
    pp.get("maxboxsize", maxboxsize);
    BoxArray ba(domain);
    ba.maxSize(maxboxsize);
    DistributionMapping dm(ba);
    EBLevelGrid eblg(ba, dm, domain, 2);
    int retval = 0;
    //ebgraph
    shared_ptr<FabArray<EBGraph> > allgraphsptr = eblg.getEBISL().getAllGraphs();
    FabArray<EBGraph>& graphsout = *allgraphsptr;
    FabArrayIO<EBGraph>::write(graphsout, string("ebgraph.plt"));

    FabArray<EBGraph> graphsin;
    FabArrayIO<EBGraph>::read(graphsin, string("ebgraph.plt"));

    retval = checkGraphs(graphsin, graphsout, eblg);
    if(retval != 0) 
    {
      amrex::Print() << "ebgraph does not match" << endl;
      return retval;
    }

    ///ebdata
    shared_ptr<FabArray<EBData> > alldataptr = eblg.getEBISL().getAllData();
    FabArray<EBData>&  dataout = *alldataptr;
    FabArrayIO<EBData>::write(dataout, string("ebdata.plt"));

    FabArray<EBData> datain;
    FabArrayIO<EBData>::read(datain, string("ebdata.plt"));

    retval = checkData(datain, dataout, eblg);
    if(retval != 0) 
    {
      amrex::Print() << "ebdata does not match" << endl;
      return retval;
    }
    //baseebcellfab
    int ncomp = 1;
    BaseEBCellFactory<int>  ebcellfact(eblg.getEBISL());
    FabArray<BaseEBCellFAB<int> >  cellout(ba, dm,  ncomp, 0, MFInfo(), ebcellfact);
    fillEBCFWithSomething(cellout, eblg);
    FabArrayIO<BaseEBCellFAB<int> >::write(cellout, string("baseebcfint_data.plt"));

    FabArray<BaseEBCellFAB<int> > cellin;
    FabArrayIO<BaseEBCellFAB<int> >::read(cellin, string("baseebcfint_data.plt"));
    retval = checkEquality(cellin, cellout, eblg);
    if(retval != 0) 
    {
      amrex::Print() << "ebcf<int> does not match" << endl;
      return retval;
    }
      

    return 0;
  }
}
/***************/
int
main(int argc, char* argv[])
{
  int retval = 0;
  amrex::Initialize(argc,argv);

  retval = amrex::testIO();
  if(retval != 0)
  {
    amrex::Print() << "simple io test failed with code " << retval << "\n";
  }
  else
  {
    amrex::Print() << "simple io test passed \n";
  }
  amrex::Finalize();
  return retval;
}
