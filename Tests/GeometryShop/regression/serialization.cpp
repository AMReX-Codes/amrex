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
#include "AMReX_EBISBox.H"

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
  void BIVF_fillWithSomething(BaseIVFAB<int>      & a_fab)
  {
    const vector<VolIndex>& vvofs = a_fab.getVoFs();
    int ival = 0;
    for(int ivof = 0; ivof < vvofs.size(); ivof++)
    {
      const VolIndex& vof = vvofs[ivof];
      for(int icomp = 0; icomp < a_fab.nComp(); icomp++)
      {
        a_fab(vof, icomp) = ival;
        ival++;
      }
    }
  }
  /****/
  int BIVF_checkEquality(BaseIVFAB<int>     & a_fab1,
                         BaseIVFAB<int>     & a_fab2,
                         const Box          & a_checkRegion,
                         const EBISBox      & a_ebis)
  {
    
    const vector<VolIndex>& vvofs1 = a_fab1.getVoFs();
    const vector<VolIndex>& vvofs2 = a_fab1.getVoFs();
    if(vvofs1.size() != vvofs2.size())
    {
      amrex::Print() << "bivfcheckequality: vector vof size mismatch" << endl;
      return -1;
    }
    if(a_fab1.nComp() != a_fab2.nComp())
    {
      amrex::Print() << "bivfcheckequality: component mismatch" << endl;
      return -2;
    }
    for(int ivof = 0 ; ivof < vvofs1.size(); ivof++)
    {
      const VolIndex& vof1 = vvofs1[ivof];
      const VolIndex& vof2 = vvofs2[ivof];
      if(vof1 != vof2)
      {
        amrex::Print() << "bivfcheckequality: vof mismatch" << endl;
        return -3;
      }
      if(a_checkRegion.contains(vof1.gridIndex()))
      {
        for(int icomp = 0; icomp < a_fab1.nComp(); icomp++)
        {
          int val1 = a_fab1(vof1, icomp);
          int val2 = a_fab2(vof2, icomp);
          if(val1 != val2)
          {
            amrex::Print() <<  "bivfcheckequality: values do not  not match at " << vof1.gridIndex();
            return -4;
          }
        }
      }
    }
    return 0;
  }

  /***************/
  int cerealTest()
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
    int nvar = SpaceDim;//just to check to see if i am assuming scalars anywhere
    for(MFIter mfi(ba, dm); mfi.isValid(); ++mfi)
    {
      const EBISBox& ebis = eblg.getEBISL()[mfi];
      const Box&     valid= eblg.getDBL  ()[mfi];
      Box grownBox = valid;
      grownBox.grow(1);
      grownBox &= ebis.getDomain();

      IntVectSet ivsGrown = ebis.getIrregIVS(grownBox);

      //BaseIVFAB      
      BaseIVFAB<int> srcBIV(ivsGrown, ebis.getEBGraph(), nvar);
      BIVF_fillWithSomething(srcBIV);
      {
        //full serialization (all meta data included)
        std::size_t nbytesbiv1 = srcBIV.nBytesFull();
        unsigned char* charbiv =  new unsigned char[nbytesbiv1];
        size_t nbytesbiv2 = srcBIV.copyToMemFull(charbiv);
        if(nbytesbiv1 != nbytesbiv2)
        {
          amrex::Print() << "byte size mismatch" << endl;
          return -10;
        }
        BaseIVFAB<int> dstBIV;
        size_t nbytesbiv3 = dstBIV.copyFromMemFull(charbiv);
        if(nbytesbiv1 != nbytesbiv3)
        {
          amrex::Print() << "byte size mismatch" << endl;
          return -11;
        }
        delete[] charbiv;

        retval = BIVF_checkEquality(srcBIV, dstBIV, grownBox, ebis);
        if(retval != 0)
        {
          amrex::Print() << "biv equality test (full) returned error" << endl;
          return retval;
        }
      }
      {
        //now test the more limited serialization stuff
        BaseIVFAB<int> dstBIV(ivsGrown, ebis.getEBGraph(), nvar);
        int startcomp = 0;
        size_t nbytesbiv1 = srcBIV.nBytes(valid, startcomp, nvar);
        unsigned char* charbiv = new unsigned char[nbytesbiv1];
        size_t nbytesbiv2 = srcBIV.copyToMem(valid, startcomp, nvar, charbiv);
        if(nbytesbiv1 != nbytesbiv2)
        {
          amrex::Print() << "byte size mismatch" << endl;
          return -12;
        }
        delete[] charbiv;

        dstBIV.copyFromMem(valid, startcomp, nvar, charbiv);
        retval = BIVF_checkEquality(srcBIV, dstBIV, valid, ebis);
        if(retval != 0)
        {
          amrex::Print() << "biv equality test (part) returned error" << endl;
          return retval;
        }

      }
      
      
    }

    return retval;
  }
}
/***************/
int
main(int argc, char* argv[])
{
  int retval = 0;
  amrex::Initialize(argc,argv);

  retval = amrex::cerealTest();
  if(retval != 0)
  {
    amrex::Print() << "serialization test failed with code " << retval << "\n";
  }
  else
  {
    amrex::Print() << "serialization test passed \n";
  }
  amrex::Finalize();
  return retval;
}
