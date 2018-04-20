
#include "AMReX_BaseIVFactory.H"
#include "AMReX_EBIndexSpace.H"
#include "AMReX_EBISLayout.H"
#include "AMReX_BoxIterator.H"
#include "AMReX_ParmParse.H"
#include "AMReX_GeometryShop.H"
#include "AMReX_WrappedGShop.H"
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
#include "AMReX_SphereIF.H"
#include "AMReX_parstream.H"
namespace amrex
{
/***************/
  int makeGeometry(Box& a_domain,
                   Real& a_dx, 
                   bool& a_hasMoments,
                   int igeom)
  {
    int eekflag =  0;

    a_hasMoments = false;
    //parse input file
    ParmParse pp;
    RealVect origin = RealVect::Zero;
    std::vector<int> n_cell;
    pp.getarr("n_cell", n_cell, 0, SpaceDim);
    int maxboxsize;
    pp.get("maxboxsize", maxboxsize);

    IntVect lo = IntVect::TheZeroVector();
    IntVect hi;
    for (int ivec = 0; ivec < SpaceDim; ivec++)
    {
      if (n_cell[ivec] <= 0)
      {
        pout() << " bogus number of cells input = " << n_cell[ivec];
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
      pout() << "all regular geometry" << "\n";
      AllRegularService regserv;
      EBIndexSpace* ebisPtr = AMReX_EBIS::instance();
      ebisPtr->define(a_domain, origin, a_dx, regserv, maxboxsize);
    }
    else if (whichgeom == 1)
    {
      pout() << "ramp geometry" << "\n";
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


      if(igeom == 0)
      {
        amrex::Print() << "using GeometryShop to generate geometric info" << endl;
        GeometryShop workshop(ramp);
        a_hasMoments = workshop.generatesHigherOrderMoments();
        //this generates the new EBIS
        EBIndexSpace* ebisPtr = AMReX_EBIS::instance();
        ebisPtr->define(a_domain, origin, a_dx, workshop, maxboxsize);
      }
      else
      {
        amrex::Print() << "using WrappedGShop to generate geometric info" << endl;
        WrappedGShop workshop(ramp);
        a_hasMoments = workshop.generatesHigherOrderMoments();
        //this generates the new EBIS
        EBIndexSpace* ebisPtr = AMReX_EBIS::instance();
        ebisPtr->define(a_domain, origin, a_dx, workshop, maxboxsize);
      }

    }
    else if(whichgeom == 5)
    {
      Real sphereRadius;
      RealVect sphereCenter;
      pp.get("sphere_radius", sphereRadius);
      EBIndexSpace* ebisPtr = AMReX_EBIS::instance();
      
      vector<Real> sphereCenterVect;
      pp.getarr("sphere_center",sphereCenterVect, 0, SpaceDim);
      for (int idir = 0; idir < SpaceDim; idir++)
      {
        sphereCenter[idir] = sphereCenterVect[idir];
      }

      pout() << "using a sphere implicit function" << "\n";

      bool negativeInside = false;
      SphereIF lalaBall(sphereRadius, sphereCenter, negativeInside);
      if(igeom == 0)
      {
        amrex::Print() << "using GeometryShop to generate geometric info" << endl;
        GeometryShop workshop(lalaBall);
        a_hasMoments = workshop.generatesHigherOrderMoments();
        ebisPtr->define(a_domain, origin, a_dx, workshop, maxboxsize);
      }
      else
      {
        amrex::Print() << "using WrappedGShop to generate geometric info" << endl;
        WrappedGShop workshop(lalaBall);
        a_hasMoments = workshop.generatesHigherOrderMoments();
        ebisPtr->define(a_domain, origin, a_dx, workshop, maxboxsize);
      }
    }
    else
    {
      //bogus which_geom
      pout() << " bogus which_geom input = "
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

  /***************/
  void BEBCF_fillWithSomething(BaseEBCellFAB<int>      & a_fab)
  {
    Box region = a_fab.box();
    EBGraph ebis = a_fab.getEBISBox().getEBGraph();
    region &= ebis.getDomain();
    IntVectSet ivs(region);
    a_fab.setVal(-1);
    
    int ival = 0;
    for(VoFIterator vofit(ivs, ebis); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      for(int icomp = 0; icomp < a_fab.nComp(); icomp++)
      {
        a_fab(vof, icomp) = ival;
        ival++;
      }
    }
  }

  /***************/
  int BEBCF_checkEquality(BaseEBCellFAB<int>      & a_fab1,
                          BaseEBCellFAB<int>      & a_fab2,
                          const Box               & a_checkRegion)
  {
    if(a_fab1.box() != a_fab2.box())
    {
      pout() << "bebcf box mismatch" << endl;
      return -1;
    }
    if(a_fab1.nComp() != a_fab2.nComp())
    {
      pout() << "bebcf box mismatch" << endl;
      return -3;
    }
    Box region = a_fab1.box();
    EBGraph ebis = a_fab1.getEBISBox().getEBGraph();
    region &= ebis.getDomain();
    region &= a_checkRegion;
    IntVectSet ivs(region);
    
    for(VoFIterator vofit(ivs, ebis); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      for(int icomp = 0; icomp < a_fab1.nComp(); icomp++)
      {
        int val1 = a_fab1(vof, icomp);
        int val2 = a_fab1(vof, icomp);
        if(val1 != val2)
        {
          pout() << "bebcf data mismatch" << endl;
          return -2;
        }
      }
    }
    return 0;
  }
  /***************/
  void BIFF_fillWithSomething(BaseIFFAB<int>      & a_fab)
  {
    const vector<FaceIndex>& vfaces = a_fab.getFaces();
    int ival = 0;
    for(int iface = 0; iface < vfaces.size(); iface++)
    {
      const FaceIndex& face = vfaces[iface];
      for(int icomp = 0; icomp < a_fab.nComp(); icomp++)
      {
        a_fab(face, icomp) = ival;
        ival++;
      }
    }
  }
  /****/
  //some slight of hand to keep from having to write code twice.
  template< class T> 
  void
  getTolerance(T& tolval)
  {
  }

  /***/
  template < > 
  void
  getTolerance(int& tolval)
  {
    tolval = 0;
  }
  /***/
  template < > 
  void
  getTolerance(Real& tolval)
  {
    tolval = 1.0e-8;
  }
  /****/
  template <class T>
  int BIVF_checkEquality(BaseIVFAB<T>     & a_fab1,
                         BaseIVFAB<T>     & a_fab2,
                         const Box        & a_checkRegion)
  {
    
    const vector<VolIndex>& vvofs1 = a_fab1.getVoFs();
    const vector<VolIndex>& vvofs2 = a_fab1.getVoFs();
    if(vvofs1.size() != vvofs2.size())
    {
      pout() << "bivfcheckequality: vector vof size mismatch" << endl;
      return -1;
    }
    if(a_fab1.nComp() != a_fab2.nComp())
    {
      pout() << "bivfcheckequality: component mismatch" << endl;
      return -2;
    }
    for(int ivof = 0 ; ivof < vvofs1.size(); ivof++)
    {
      const VolIndex& vof1 = vvofs1[ivof];
      const VolIndex& vof2 = vvofs2[ivof];
      if(vof1 != vof2)
      {
        pout() << "bivfcheckequality: vof mismatch" << endl;
        return -3;
      }
      if(a_checkRegion.contains(vof1.gridIndex()))
      {
        for(int icomp = 0; icomp < a_fab1.nComp(); icomp++)
        {
          T val1 = a_fab1(vof1, icomp);
          T val2 = a_fab2(vof2, icomp);
          T diff = std::abs(val1-val2);
          T tol;
          getTolerance<T>(tol);
          if(diff > tol)
          {
            pout() <<  "bivfcheckequality: values do not  not match at " << vof1.gridIndex();
            return -4;
          }
        }
      }
    }
    return 0;
  }

  /****/
  template <class T>
  int BIFF_checkEquality(BaseIFFAB<T>     & a_fab1,
                         BaseIFFAB<T>     & a_fab2,
                         const Box        & a_checkRegion)
  {
    
    const vector<FaceIndex>& vfaces1 = a_fab1.getFaces();
    const vector<FaceIndex>& vfaces2 = a_fab1.getFaces();
    if(vfaces1.size() != vfaces2.size())
    {
      pout() << "biff_checkequality: vector vof size mismatch" << endl;
      return -1;
    }
    if(a_fab1.nComp() != a_fab2.nComp())
    {
      pout() << "biff_checkequality: component mismatch" << endl;
      return -2;
    }
    for(int iface = 0 ; iface < vfaces1.size(); iface++)
    {
      const FaceIndex& face1 = vfaces1[iface];
      const FaceIndex& face2 = vfaces2[iface];
      if(face1 != face2)
      {
        pout() << "biff_checkequality: vof mismatch" << endl;
        return -3;
      }
      if(a_checkRegion.contains(face1.gridIndex(Side::Lo)))
      {
        for(int icomp = 0; icomp < a_fab1.nComp(); icomp++)
        {
          T val1 = a_fab1(face1, icomp);
          T val2 = a_fab2(face2, icomp);
          T diff = std::abs(val1-val2);
          T tol;
          getTolerance<T>(tol);
          if(diff > tol)
          {
            pout() <<  "biff_checkequality: values do not  not match at " << face1.gridIndex(Side::Lo);
            return -4;
          }
        }
      }
    }
    return 0;
  }

  /***************/
  int EBG_checkEquality(EBGraph            & a_ebg1,
                        EBGraph            & a_ebg2,
                        const Box          & a_checkRegion)
  {
    if(a_ebg1.getDomain() != a_ebg2.getDomain())
    {
      pout() << "ebg_checkequality: domain mismatch" << endl;
      return -1;
    }
    if(a_ebg1.getRegion() != a_ebg2.getRegion())
    {
      pout() << "ebg_checkequality: region mismatch" << endl;
      return -2;
    }
    for(BoxIterator bit(a_ebg1.getRegion()); bit.ok(); ++bit)
    {
      if(a_checkRegion.contains(bit()))
        {
          vector<VolIndex> vvofs1 = a_ebg1.getVoFs(bit());
          vector<VolIndex> vvofs2 = a_ebg1.getVoFs(bit());
          if(vvofs1.size() != vvofs2.size())
          {
            pout() << "ebg_checkequality: vector vof size mismatch" << endl;
            return -3;
          }
          for(int ivof = 0; ivof < vvofs1.size(); ivof++)
          {
            const VolIndex& vof1 = vvofs1[ivof];
            const VolIndex& vof2 = vvofs2[ivof];
            if(vof1 != vof2)
            {
              pout() << "ebg_checkequality: vof mismatch" << endl;
              return -4;
            }
            for(int idir = 0; idir < SpaceDim; idir++)
            {
              for(SideIterator sit; sit.ok(); ++sit)
              {
                vector<FaceIndex> vfaces1 = a_ebg1.getFaces(vof1, idir, sit());
                vector<FaceIndex> vfaces2 = a_ebg2.getFaces(vof2, idir, sit());
                if(vfaces1.size() != vfaces2.size())
                {
                  pout() << "ebg_checkequality: vector vof size mismatch 2" << endl;
                  return -5;
                }
                for(int iface = 0; iface < vfaces1.size(); iface++)
                {
                  const FaceIndex& face1 = vfaces1[iface];
                  const FaceIndex& face2 = vfaces2[iface];
                  if(face1 != face2)
                  {
                    pout() << "ebg_checkequality: face mismatch" << endl;
                    return -6;
                  }
                }

              }
            }
          }
        }
        
    }
    return 0;
  }
  /***************/
  int EBD_checkEquality(EBData             & a_ebd1,
                        EBData             & a_ebd2,
                        const Box          & a_checkRegion)
  {

    if(a_ebd1.getRegion() != a_ebd2.getRegion())
    {
      pout() << "ebd_checkequality: region mismatch" << endl;
      return -6;
    }

    BaseIVFAB<Real>& voldata1 = a_ebd1.getVolData();
    BaseIVFAB<Real>& voldata2 = a_ebd2.getVolData();
    int retval = BIVF_checkEquality<Real>(voldata1, voldata2, a_checkRegion);
    if(retval != 0)
    {
      pout() << "ebd_checkequality: vol data mismatch" << endl;
      return -7;
    }
    
    for(int idir = 0; idir < SpaceDim; idir++)
    {
      BaseIFFAB<Real>& facedata1 = a_ebd1.getFaceData(idir);
      BaseIFFAB<Real>& facedata2 = a_ebd2.getFaceData(idir);
      retval = BIFF_checkEquality<Real>(facedata1, facedata2, a_checkRegion);
      if(retval != 0)
      {
        pout() << "ebd_checkequality: face data mismatch" << endl;
        return -7 - 10*idir;
      }
    }
    return 0;
  }
  /***************/
  int cerealTest(int igeom)
  {
    Box domain;
    Real dx;
    pout() << "making the geometry" << endl;
    bool hasMoments;
    makeGeometry(domain, dx, hasMoments, igeom);
    int maxboxsize;
    ParmParse pp;
    pp.get("maxboxsize", maxboxsize);
    BoxArray ba(domain);
    ba.maxSize(maxboxsize);
    DistributionMapping dm(ba);
    pout() << "making the eblevelgrid" << endl;
    EBLevelGrid eblg(ba, dm, domain, 2);
    int retval = 0;
    int nvar = SpaceDim;//just to check to see if i am assuming scalars anywhere
    pout() << "about to loop over stuff" << endl;
    int ibox = 0;
    for(MFIter mfi(ba, dm); mfi.isValid(); ++mfi)
    {
      const EBISBox& ebis = eblg.getEBISL()[mfi];
      const Box&     valid= eblg.getDBL  ()[mfi];
      Box grownBox = valid;
      grownBox.grow(1);
      grownBox &= ebis.getDomain();

      IntVectSet ivsGrown = ebis.getIrregIVS(grownBox);

      //BaseIVFAB      
      {
        BaseIVFAB<int> srcBIV(ivsGrown, ebis.getEBGraph(), nvar);
        BIVF_fillWithSomething(srcBIV);
        {
          //full serialization (all meta data included)
          std::size_t nbytesbiv1 = srcBIV.nBytesFull();
          unsigned char* charbiv =  new unsigned char[nbytesbiv1];
          size_t nbytesbiv2 = srcBIV.copyToMemFull(charbiv);
          if(nbytesbiv1 != nbytesbiv2)
          {
            pout() << "byte size mismatch" << endl;
            return -10;
          }
          BaseIVFAB<int> dstBIV;
          size_t nbytesbiv3 = dstBIV.copyFromMemFull(charbiv);
          if(nbytesbiv1 != nbytesbiv3)
          {
            pout() << "byte size mismatch" << endl;
            return -11;
          }

          retval = BIVF_checkEquality<int>(srcBIV, dstBIV, grownBox);
          if(retval != 0)
          {
            pout() << "biv equality test (full) returned error" << endl;
            return retval;
          }
          delete[] charbiv;
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
            pout() << "byte size mismatch" << endl;
            return -12;
          }

          size_t nbytesbiv3 = dstBIV.copyFromMem(valid, startcomp, nvar, charbiv);
          if(nbytesbiv3 != nbytesbiv2)
          {
            pout() << "byte size mismatch" << endl;
            return -112;
          }

          retval = BIVF_checkEquality<int>(srcBIV, dstBIV, valid);
          if(retval != 0)
          {
            pout() << "biv equality test (part) returned error" << endl;
            return retval;
          }
          delete[] charbiv;

        }
      }



      //BaseEBCellFAB      
      {
        BaseEBCellFAB<int> src(ebis, grownBox, nvar);
        BEBCF_fillWithSomething(src);
        {
          //full serialization (all meta data included)
          std::size_t nbytes1 = src.nBytesFull();
          unsigned char* charbuf =  new unsigned char[nbytes1];
          size_t nbytes2 = src.copyToMemFull(charbuf);
          if(nbytes1 != nbytes2)
          {
            pout() << "byte size mismatch" << endl;
            return -210;
          }
          BaseEBCellFAB<int> dst;
          size_t nbytes3 = dst.copyFromMemFull(charbuf);
          if(nbytes1 != nbytes3)
          {
            pout() << "byte size mismatch" << endl;
            return -211;
          }

          retval = BEBCF_checkEquality(src, dst, grownBox);
          if(retval != 0)
          {
            pout() << "ebcf equality test (full) returned error" << endl;
            return retval;
          }
          delete[] charbuf;
        }
        {
          //now test the more limited serialization stuff
          BaseEBCellFAB<int> dst(ebis, grownBox, nvar);
          int startcomp = 0;
          size_t nbytes1 = src.nBytes(valid, startcomp, nvar);
          unsigned char* charbuf = new unsigned char[nbytes1];
          size_t nbytes2 = src.copyToMem(valid, startcomp, nvar, charbuf);
          if(nbytes1 != nbytes2)
          {
            pout() << "byte size mismatch" << endl;
            return -212;
          }

          size_t nbytes3 = dst.copyFromMem(valid, startcomp, nvar, charbuf);
          if(nbytes3 != nbytes2)
          {
            pout() << "byte size mismatch" << endl;
            return -213;
          }

          retval = BEBCF_checkEquality(src, dst, valid);
          if(retval != 0)
          {
            pout() << "ebcf equality test (part) returned error" << endl;
            return retval;
          }
          delete[] charbuf;

        }
      }


      //BaseIFFAB      
      for(int idir = 0; idir < SpaceDim; idir++)
      {
        //full
        {
          BaseIFFAB<int> src(ivsGrown, ebis.getEBGraph(), idir, nvar);
          BIFF_fillWithSomething(src);
          {
            //full serialization (all meta data included)
            std::size_t nbytes1 = src.nBytesFull();
            unsigned char* charbuf =  new unsigned char[nbytes1];
            size_t nbytes2 = src.copyToMemFull(charbuf);
            if(nbytes1 != nbytes2)
            {
              pout() << "byte size mismatch" << endl;
              return -110;
            }
            BaseIFFAB<int> dst;
            size_t nbytes3 = dst.copyFromMemFull(charbuf);
            if(nbytes1 != nbytes3)
            {
              pout() << "byte size mismatch" << endl;
              return -111;
            }

            retval = BIFF_checkEquality<int>(src, dst, grownBox);
            if(retval != 0)
            {
              pout() << " equality test (full) returned error" << endl;
              return retval;
            }
            delete[] charbuf;
          }
          //now test the more limited serialization stuff
          {
            BaseIFFAB<int> dst(ivsGrown, ebis.getEBGraph(), idir, nvar);
            int startcomp = 0;
            size_t nbytes1 = src.nBytes(valid, startcomp, nvar);
            unsigned char* charbuf = new unsigned char[nbytes1];
            size_t nbytes2 = src.copyToMem(valid, startcomp, nvar, charbuf);
            if(nbytes1 != nbytes2)
            {
              pout() << "byte size mismatch" << endl;
              return -112;
            }

            size_t nbytes3 = dst.copyFromMem(valid, startcomp, nvar, charbuf);
            if(nbytes3 != nbytes2)
            {
              pout() << "byte size mismatch" << endl;
              return -113;
            }

            retval = BIFF_checkEquality<int>(src, dst, valid);
            if(retval != 0)
            {
              pout() << " equality test (part) returned error" << endl;
              return retval;
            }
            delete[] charbuf;

          }
        }

      }


      //EBGraph
      {
        EBGraph src = ebis.getEBGraph();
        {
          //full serialization (all meta data included)
          std::size_t nbytes1 = src.nBytesFull();
          unsigned char* buff =  new unsigned char[nbytes1];
          size_t nbytes2 = src.copyToMemFull(buff);
          if(nbytes1 != nbytes2)
          {
            pout() << "ebg byte size mismatch" << endl;
            return -13;
          }
          EBGraph dst;
          size_t nbytes3 = dst.copyFromMemFull(buff);
          if(nbytes1 != nbytes3)
          {
            pout() << "ebg byte size mismatch 2" << endl;
            return -14;
          }

          retval = EBG_checkEquality(src, dst, grownBox);
          if(retval != 0)
          {
            pout() << " ebg equality test (full) returned error" << endl;
            return retval;
          }
          delete[] buff;
        }
        {
          //now test the more limited serialization stuff
          EBGraph dst(src.getRegion());
          dst.setDomain(src.getDomain());
          int startcomp = 0;
          size_t nbytes1 = src.nBytes(valid, startcomp, nvar);
          unsigned char* buff = new unsigned char[nbytes1];
          size_t nbytes2 = src.copyToMem(valid, startcomp, nvar, buff);
          if(nbytes1 != nbytes2)
          {
            pout() << "ebg byte size mismatch 3" << endl;
            return -15;
          }

          size_t nbytes3 = dst.copyFromMem(valid, startcomp, nvar, buff);
          if(nbytes1 != nbytes3)
          {
            pout() << "ebg byte size mismatch 4" << endl;
            return -16;
          }

          retval = EBG_checkEquality(src, dst, valid);
          if(retval != 0)
          {
            pout() << "equality test (part) returned error" << endl;
            return retval;
          }
          delete[] buff;

        }
      }


      //EBData
      {
        EBData src = ebis.getEBData();
        {
          //full serialization (all meta data included)
          std::size_t nbytes1 = src.nBytesFull();
          unsigned char* buff =  new unsigned char[nbytes1];
          size_t nbytes2 = src.copyToMemFull(buff);
          if(nbytes1 != nbytes2)
          {
            pout() << "ebg byte size mismatch" << endl;
            return -13;
          }
          EBData dst;
          size_t nbytes3 = dst.copyFromMemFull(buff);
          if(nbytes1 != nbytes3)
          {
            pout() << "ebg byte size mismatch 2" << endl;
            return -14;
          }

          retval = EBD_checkEquality(src, dst, grownBox);
          if(retval != 0)
          {
            pout() << " ebg equality test (full) returned error" << endl;
            return retval;
          }
          delete[] buff;
        }
        {
          //now test the more limited serialization stuff
          EBData dst;
          dst.define(ebis.getEBGraph(), src.getRegion(), dx, hasMoments);

          int startcomp = 0;
          size_t nbytes1 = src.nBytes(valid, startcomp, nvar);
          unsigned char* buff = new unsigned char[nbytes1];
          size_t nbytes2 = src.copyToMem(valid, startcomp, nvar, buff);
          if(nbytes1 != nbytes2)
          {
            pout() << "ebg byte size mismatch 3" << endl;
            return -15;
          }

          size_t nbytes3 = dst.copyFromMem(valid, startcomp, nvar, buff);
          if(nbytes1 != nbytes3)
          {
            pout() << "ebg byte size mismatch 4" << endl;
            return -16;
          }

          retval = EBD_checkEquality(src, dst, valid);
          if(retval != 0)
          {
            pout() << "equality test (part) returned error" << endl;
            return retval;
          }
          delete[] buff;

        }
      }
      ibox++; 
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

  for(int igeom = 0; igeom <= 1; igeom++)
  {
    retval = amrex::cerealTest(igeom);
    if(retval != 0)
    {
      amrex::Print() << "serialization test failed with code " << retval << "\n";
      return retval;
    }

  }
  amrex::Print() << "serialization test passed \n";

  amrex::Finalize();
  return retval;
}
