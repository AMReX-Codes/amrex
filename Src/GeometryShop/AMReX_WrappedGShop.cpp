
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include "AMReX_WrappedGShop.H"

#include "AMReX_NoRefinement.H"
#include "AMReX_FixedRefinement.H"
#include "AMReX_DivNormalRefinement.H"
#include "AMReX_NormalDerivativeNew.H"
#include "AMReX_MinimalCCCM.H"
#include "AMReX_MomentIterator.H"
#include "AMReX_EB_TYPEDEFS.H"
#include "AMReX_PolyGeom.H"
#include "AMReX_RealVect.H"
#include "AMReX_EBGeomDebugOut.H"

Real amrex::WrappedGShop::s_relativeTol = 0.1;

namespace amrex
{
/*********************************************/
/*********************************************/
  WrappedGShop::
  WrappedGShop(const RefCountedPtr<BaseIF>  &      a_baseIF,
               const RealVect               &      a_origin,
               const Real                   &      a_dx,
               const Box          &      a_domain,
               int                                 a_minNumberRefines ,
               int                                 a_maxNumberRefines )
  {
    CH_TIME("WrappedGShop::WrappedGShop");
    m_baseIF  = a_baseIF;

    m_minNumberRefines = a_minNumberRefines;
    m_maxNumberRefines = a_maxNumberRefines;

    m_domain  = a_domain;
    m_origin  = a_origin;
            
    m_order   = 0;
    m_degreeP = CH_EBIS_ORDER + 1;
    m_refCrit = RefCountedPtr<WGSRefinementCriterion> (new WGSRefinementCriterion());
  }

/*********************************************/
  bool 
  WrappedGShop::
  isRegular(const Box           & a_region,
            const Box & a_domain,
            const RealVect      & a_origin,
            const Real          & a_dx) const
  {
    CH_TIME("WrappedGShop::isRegular");

    // first check any of the Box corners are outside, and return false
    // right away. (bvs)
    RealVect physCorner;
    IntVect lo = a_region.smallEnd();
    IntVect len = a_region.size();
    Box unitBox(IntVect::Zero, IntVect::Unit);
    for (BoxIterator bit(unitBox); bit.ok(); ++bit)
    {
      IntVect current = lo + len*bit();
      for (int idir = 0; idir < CH_SPACEDIM; ++idir)
      {
        physCorner[idir] = a_dx*(current[idir]) + a_origin[idir];
      }
      Real functionValue = m_baseIF->value(physCorner);
      if (functionValue > 0.0 )
      {
        return false;
      }
    }

    return isRegularEveryPoint(a_region, a_domain, a_origin, a_dx);
  }
////
  bool 
  WrappedGShop::
  isRegularEveryPoint(const Box&           a_region,
                      const Box& a_domain,
                      const RealVect&      a_origin,
                      const Real&          a_dx) const
  {
    CH_TIME("WrappedGShop::isRegularEveryPoint");

    // All corner indices for the current box
    Box allCorners(a_region);
    allCorners.surroundingNodes();

    RealVect physCorner;
    BoxIterator bit(allCorners);
    // If every corner is inside, the box is regular
    for (int i=0; i<2; i++)
    {
      for (; bit.ok(); ++bit, ++bit)
      {
        // Current corner
        const IntVect& corner = bit();

        // Compute physical coordinate of corner

        for (int idir = 0; idir < CH_SPACEDIM; ++idir)
        {
          physCorner[idir] = a_dx*(corner[idir]) + a_origin[idir];
        }

        // If the implicit function value is positive then the current corner is
        // covered
        Real functionValue = m_baseIF->value(physCorner);

        if (functionValue > 0.0 )
        {
          return false;
        }
      }
      bit.reset();
      ++bit;
    }

    return true;
  }
/*********************************************/
  bool
  WrappedGShop::
  isCovered(const Box           & a_region,
            const Box & a_domain,
            const RealVect      & a_origin,
            const Real          & a_dx) const
  {
    CH_TIME("WrappedGShop::isCovered");

    // first check any of the Box corners are inside, and return false
    // right away. (bvs)
    RealVect physCorner;
    IntVect lo = a_region.smallEnd();
    IntVect len = a_region.size();
    Box unitBox(IntVect::Zero, IntVect::Unit);
    for (BoxIterator bit(unitBox); bit.ok(); ++bit)
    {
      IntVect current = lo + len*bit();
      for (int idir = 0; idir < CH_SPACEDIM; ++idir)
      {
        physCorner[idir] = a_dx*(current[idir]) + a_origin[idir];
      }
      Real functionValue = m_baseIF->value(physCorner);
      if (functionValue < 0.0 )
      {
        return false;
      }
    }

    return isCoveredEveryPoint(a_region, a_domain, a_origin, a_dx);
  }

/**********************************************/
  bool 
  WrappedGShop::
  isCoveredEveryPoint(const Box&           a_region,
                      const Box& a_domain,
                      const RealVect&      a_origin,
                      const Real&          a_dx) const
  {
    CH_TIME("WrappedGShop::isCoveredEveryPoint");

    // All corner indices for the current box
    Box allCorners(a_region);
    allCorners.surroundingNodes();

    RealVect physCorner;
    BoxIterator bit(allCorners);
    // If every corner is inside, the box is regular
    for (int i=0; i<2; i++)
    {
      for (; bit.ok(); ++bit, ++bit)
      {
        // Current corner
        IntVect corner = bit();

        // Compute physical coordinate of corner

        for (int idir = 0; idir < CH_SPACEDIM; ++idir)
        {
          physCorner[idir] = a_dx*(corner[idir]) + a_origin[idir];
        }

        // If the implicit function value is negative then the current corner is
        // not covered
        Real functionValue = m_baseIF->value(physCorner);

        if (functionValue < 0.0 )
        {
          return false;
        }
      }
      bit.reset();
      ++bit;
    }

    return true;
  }
/******/
  bool 
  WrappedGShop::
  isIrregular(const Box           & a_region,
              const Box & a_domain,
              const RealVect      & a_origin,
              const Real          & a_dx) const
  {
    CH_TIME("WrappedGShop::isIrregular");

    // first check if some Box corners are inside and some are outside, and return
    // true right away. (bvs)
    RealVect physCorner;
    IntVect lo = a_region.smallEnd();
    IntVect len = a_region.size();
    for (int idir = 0; idir < CH_SPACEDIM; ++idir)
    {
      physCorner[idir] = a_dx*(lo[idir]) + a_origin[idir];
    }
    Real originVal = m_baseIF->value(physCorner);

    Box unitBox(IntVect::Zero, IntVect::Unit);
    for (BoxIterator bit(unitBox); bit.ok(); ++bit)
    {
      IntVect current = lo + len*bit();
      for (int idir = 0; idir < CH_SPACEDIM; ++idir)
      {
        physCorner[idir] = a_dx*(current[idir]) + a_origin[idir];
      }
      Real functionValue = m_baseIF->value(physCorner);
      if (functionValue * originVal < 0.0 )
      {
        return true;
      }
    }

    return !(isRegularEveryPoint(a_region, a_domain, a_origin, a_dx) ||
             isCoveredEveryPoint(a_region, a_domain, a_origin, a_dx));
  }

/**********************************************/
  bool
  WrappedGShop::
  onBoxBoundary(const IntVect        & a_iv, 
                const Box            & a_box,
                const int            & a_dir,
                const Side::LoHiSide & a_sd) const
  {
    bool retval = false;
    if(a_sd == Side::Lo)
    {
      const IntVect& ivlo = a_box.smallEnd();
      retval = (a_iv[a_dir] == ivlo[a_dir]);
    }
    else if (a_sd == Side::Hi)
    {
      const IntVect& ivhi = a_box.bigEnd();
      retval = (a_iv[a_dir] == ivhi[a_dir]);
    }
    else
    {
      MayDay::Error("bogus side");
    }
    return retval;
  }
/**********************************************/
  void 
  WrappedGShop::
  agglomerateMoments(IrregNode              & a_node, 
                     const Vector<IrregNode>& a_refNodes,
                     const Box              & a_refBox,
                     const Real             & a_fineDx,
                     const Real             & a_coarDx) const
  {
    bool verbose = false;
    //this bit is imporant since we are adding the shifted
    //fine moments incrementally (so we have to initialize to zero)
    a_node.setMomentsToZero();
    //the origin does not matter here --- local just used for differences
    RealVect origin = RealVect::Zero;
    IntVect ivcoar = a_node.m_cell;
    VolIndex vcoar(a_node.m_cell, a_node.m_cellIndex);
    RealVect coarloc = EBArith::getVoFLocation(vcoar, a_coarDx*RealVect::Unit, origin);
    for(int inode = 0; inode < a_refNodes.size(); inode++)
    {
      const IrregNode& fineNd = a_refNodes[inode];
      IntVect ivfine = fineNd.m_cell;
      VolIndex vfine(ivfine, fineNd.m_cellIndex);
      RealVect fineloc = EBArith::getVoFLocation(vfine, a_fineDx*RealVect::Unit, origin);
      RealVect shiftRV = fineloc - coarloc;

      RvSpaceDim shiftSpaceDim;
      EBArith::convertToITM(shiftSpaceDim, shiftRV);
      {
        //first shift and add up the volume moments
        //not a reference so I can change the fine one by shifting it
        IndMomSpaceDim  volmomFine = fineNd.m_volumeMoments;
        volmomFine.shift(shiftSpaceDim);
        IndMomSpaceDim& volmomCoar = a_node.m_volumeMoments;
        volmomCoar += volmomFine;
      }

      {
        //now shift and add up the  EB moments
        //not a reference so I can change the fine one by shifting it
        IndMomSpaceDim  ebfmomFine = fineNd.m_EBMoments;
        ebfmomFine.shift(shiftSpaceDim);
        IndMomSpaceDim& ebfmomCoar = a_node.m_EBMoments;
        ebfmomCoar += ebfmomFine;
      }
        
      //finally shift and add the eb moments
      for(int facedir = 0; facedir < SpaceDim; facedir++)
      {
        RvSDMinOne shiftSDMinOne;
        EBArith::convertToITM(shiftSDMinOne, shiftRV, facedir);
        for(SideIterator sit; sit.ok(); ++sit)
        {
          //only add face moments if they are on the boundary of the  coarse cell
          if(onBoxBoundary(ivfine,a_refBox, facedir, sit()))
          {
            int iindex = IrregNode::index(facedir, sit());
            //not a reference so I can change the fine one by shifting it
            IndMomSDMinOne  facmomFine = fineNd.m_faceMoments[iindex];
            facmomFine.shift(shiftSDMinOne);
            IndMomSDMinOne& facmomCoar = a_node.m_faceMoments[iindex];
            facmomCoar += facmomFine;
          }
        }
      }
    }

    bool fixMoments = true;
    bool foundBogus = checkNodeMoments(a_node, a_coarDx, fixMoments, s_relativeTol);
    if(foundBogus && verbose)
    {
      MayDay::Warning("still found bogus after agglomeration");
    }
  }
/**********************************************/
  bool 
  WrappedGShop::
  checkNodeMoments(IrregNode & a_node, 
                   const Real& a_dx,
                   const bool& a_bindMoments,
                   const Real& a_tolerance) const
  {
    bool verbose = false;
    bool retval =  false;
    /**
       Moments are centered at the center of the cell.
       For each of these moments I shift them to the lowest 
       corner of the cell, where I know what the bounds of the integrals 
       is (lower bound always zero, upper bound = dx^d dx^px dx^py dx^pz
       If the shifted moment is out of bounds, I bound it throw a message.
    **/
    //IntVect iv = a_node.m_cell;
    //IntVect ivdebug(D_DECL(7, 4, 0));
    //if(iv == ivdebug)
    //  {
    //    pout() <<"here we are" << endl;
    //  }

    Real kappa = a_node.m_volumeMoments[IvSpaceDim::Zero];
    for(int idir = 0; idir < SpaceDim; idir++)
    {
      kappa /= a_dx;
    }

    if(verbose)
    {
      pout() << "cell = " << a_node.m_cell << ", kappa = " << kappa << endl;
    }
    Vector<IvSpaceDim> ivboguscell;

    checkMoments<SpaceDim, CH_EBIS_ORDER>(a_node.m_volumeMoments, ivboguscell, 
                                          a_dx, a_tolerance, false, a_bindMoments);
    if(ivboguscell.size() > 0) 
    {
      retval = true;
    }
    //if we are verbose, print out what powers were off
    if(verbose && (ivboguscell.size() > 0))
    {
      pout() << "for node = " << a_node.m_cell << ", kappa = " << kappa << endl; 
      pout() << "bogus volume moments found and bound were: " << endl;
      for(int ivec = 0; ivec < ivboguscell.size(); ivec++)
      {
        pout() << ivboguscell[ivec]  << " ";
      }
      pout() << endl;
    }

    /**/
    checkMoments<SpaceDim, CH_EBIS_ORDER>(a_node.m_EBMoments, ivboguscell, 
                                          a_dx, a_tolerance, true, a_bindMoments);
    if(verbose && (ivboguscell.size() > 0))
    {
      Real alpha = a_node.m_EBMoments[IvSpaceDim::Zero];
      for(int idir = 0; idir < SpaceDim-1; idir++)
      {
        alpha /= a_dx;
      }
      pout() << "for node = " << a_node.m_cell << endl; 
      pout() << "bogus EB moments found and bound were: " << endl;
      for(int ivec = 0; ivec < ivboguscell.size(); ivec++)
      {
        pout() << ivboguscell[ivec]  << " ";
      }
      pout() << endl;
    }
    if(ivboguscell.size() > 0) 
    {
      retval = true;
    }
    /**/
    for(SideIterator sit; sit.ok(); ++sit)
    {
      for (int idir = 0; idir < SpaceDim; ++idir)
      {
        IndMomSDMinOne regmom;
        regmom.setRegular(a_dx);
        Vector<IvSDMinOne> ivbogusface;
        int iindex = IrregNode::index(idir, sit());
        checkMoments<SpaceDim-1, CH_EBIS_ORDER>(a_node.m_faceMoments[iindex], ivbogusface, 
                                                a_dx, a_tolerance, false, a_bindMoments);
        Real alpha = a_node.m_faceMoments[iindex][IvSDMinOne::Zero];
        for(int jdir = 0; jdir < SpaceDim-1; jdir++)
        {
          alpha /= a_dx;
        }
        if(ivbogusface.size() > 0) 
        {
          retval = true;
        }
        if(verbose && (ivbogusface.size() > 0))
        {
          pout() << "for node = " << a_node.m_cell << endl; 
          pout() << "for direction = " << idir << " and side = " << sign(sit()) << endl;
          pout() << "alpha = " << alpha << endl;
          pout() << "bogus Face moments found and bound were: " << endl;
          for(int ivec = 0; ivec < ivbogusface.size(); ivec++)
          {
            pout() << ivbogusface[ivec]  << " ";
          }
          pout() << endl;
        }
      }
    }
    return retval;
  }
/*********************************************/
  void 
  WrappedGShop::
  fillGraph(BaseFab<int>&       a_regIrregCovered,
            Vector<IrregNode>&  a_nodes,
            const Box&          a_validRegion,
            const Box&          a_ghostRegion,
            const Box & a_domain,
            const RealVect      & a_origin,
            const Real          & a_dx) const
  {
    CH_TIME("WrappedGShop::fillGraph");

    CH_assert(a_domain.contains(a_ghostRegion));

    IntVectSet ivsirreg = IntVectSet(DenseIntVectSet(a_ghostRegion, false));

    {
      CH_TIME("boxiterator loop");
      for (BoxIterator bit(a_ghostRegion); bit.ok(); ++bit)
      {
        const IntVect iv =bit();
        RvgDim cellCenter;
        for (int idir = 0;idir < SpaceDim; ++idir)
        {
          cellCenter[idir] = a_dx*(iv[idir] +0.5) + a_origin[idir];
        }

        //member data: sign(chosen from -1,0,1) of each vertex,
        //location of each edge intersection, cellCenter,normal and gradNormal

        IndexTM<Real, SpaceDim> vectdx;
        vectdx.setAll(a_dx);

        int degreeP = m_order + m_degreeP;
        IFData<SpaceDim> edgeData(*m_baseIF, vectdx, cellCenter,  degreeP);


        //create a CutCellMoment object, in order to detect whether any face coincides with the interface
        CutCellMoments <SpaceDim> cutCell(edgeData);
        if (cutCell.isCovered())
        {
          //set covered cells to -1
          a_regIrregCovered(iv, 0) = -1;
        }
        else if (cutCell.isRegular())
        {
          //set regular cells to 1
          a_regIrregCovered(iv, 0) =  1;
        }
        else
        {
          //set irregular cells to 0
          //irregular if any face coincides with interface and edgeData.m_allVerticesIn = true
          a_regIrregCovered(iv, 0) =  0;
          if (a_validRegion.contains(iv))
          {
            ivsirreg |= iv;
          }
        }
      }
    }
    //now loop through irregular cells and make nodes for each  one.

    for (IVSIterator ivsit(ivsirreg); ivsit.ok(); ++ivsit)
    {
      IntVect iv = ivsit();
      CH_TIME("fillGraph::endOfirregularCellLoop");
      IrregNode newNode;

      fillNewNode(newNode,
                  ivsirreg,
                  a_domain,
                  a_origin,
                  a_dx,
                  ivsit());

      bool fixMoments = true;
      checkNodeMoments(newNode, a_dx, fixMoments, s_relativeTol);      
      a_nodes.push_back(newNode);
    } //end loop over cells in the box
  }
/******/
  bool
  WrappedGShop::
  needToRefine(IrregNode       & a_node,
               const Real      & a_dx,
               const int       & a_numRefSoFar) const
  {
    bool fixMoments = false;
    bool hadBadMoments = checkNodeMoments(a_node, a_dx, fixMoments, s_relativeTol);
    bool retval = false;
    if(a_numRefSoFar < m_minNumberRefines)
    {
      retval = true;
      //this one really should end it.
      return retval;
    }
    else if(a_numRefSoFar >= m_maxNumberRefines)
    {
      retval = false;
      //this one really should end it.
      return retval;
    }
    else if(hadBadMoments)
    {
      retval = true;
    }
    else if(m_refCrit->refineHere(a_node, a_node.m_cell, a_dx) == true)
    {
      retval = true;
    }
    return retval;
  }
/**********************************************/
  void 
  WrappedGShop::
  fillNewNode(IrregNode           &     a_node,
              const IntVectSet    &     a_ivsIrreg,
              const Box &     a_domain,
              const RealVect      &     a_origin,
              const Real          &     a_dx,
              const IntVect       &     a_iv) const
  {
    computeVoFInternals(a_node,
                        a_ivsIrreg,
                        a_domain,
                        a_origin,
                        a_dx,
                        a_iv);

    /**/
    Box refBox(a_iv, a_iv);
    refBox.refine(2);
    Real refDx = a_dx/2.;
    Box refDom = a_domain;
    refDom.refine(2);
    //if we need to refine the cell, refine all the cells and
    //add up their moments
    int numRefSoFar = 0;
    //  bool refinedThisOne = false;
    while(needToRefine(a_node, a_dx, numRefSoFar))
    {
      //refinedThisOne = true;
      Vector<IrregNode> refNodes;
      for(BoxIterator boxit(refBox); boxit.ok(); ++boxit)
      {
        IrregNode refNode;
        const IntVect iv =boxit();
        RvgDim cellCenter;
        for (int idir = 0;idir < SpaceDim; ++idir)
        {
          cellCenter[idir] = refDx*(iv[idir] +0.5) + a_origin[idir];
        }

        //member data: sign(chosen from -1,0,1) of each vertex,
        //location of each edge intersection, cellCenter,normal and gradNormal

        IndexTM<Real, SpaceDim> vectdx;
        vectdx.setAll(refDx);

        int degreeP = m_order + m_degreeP;
        IFData<SpaceDim> edgeData(*m_baseIF, vectdx, cellCenter,  degreeP);


        //create a CutCellMoment object, in order to detect whether any face coincides with the interface
        CutCellMoments<SpaceDim> cutCell(edgeData);
        //this gets used even in covered and regular cases
        refNode.m_cell = boxit();
        if (cutCell.isCovered())
        {
          refNode.setMomentsToZero();
        }
        else if (cutCell.isRegular())
        {
          refNode.setMomentsToRegular(refDx);
        }
        else
        {
          //these refined nodes will have the wrong
          //integers in a lot of the arcs because 
          //ivsirreg is nonsense in this context
          IntVectSet fakeIVSIrreg;
          computeVoFInternals(refNode,
                              fakeIVSIrreg,
                              refDom,
                              a_origin,
                              refDx,
                              boxit());
        }
        refNodes.push_back(refNode);
      }
      //this leaves the arcs of a_node alone because
      //a lot of the arcs in the refined case are wrong
      //because ivsirreg was wrong.
      agglomerateMoments(a_node, refNodes, refBox, refDx, a_dx);
                         
      refDom.refine(2);
      refBox.refine(2);
      refDx /= 2.;
      numRefSoFar++;
    }
//  if(refinedThisOne)
//    {
//      pout() << "refined " << a_iv << endl;
//    }
//  else
//    {
//      pout() << "did not refine " << a_iv << endl;
//    }
    /**/
  }
/*********************************************/
  void 
  WrappedGShop::
  computeVoFInternals(IrregNode                      &     a_node,
                      const IntVectSet               &     a_ivsIrreg,
                      const Box            &     a_domain,
                      const RealVect                 &     a_origin,
                      const Real                     &     a_dx,
                      const IntVect                  &     a_iv) const
  {
    CH_TIME("WrappedGShop::ComputeVofInternals");

    //for each CutCellMoments<dim>, we record the cell Center
    //(in physical coordinates at the global dimension)
    RvgDim cellCenter;
    for (int idir = 0;idir < SpaceDim; ++idir)
    {
      cellCenter[idir] = a_dx*(a_iv[idir] +0.5) + a_origin[idir];
    }

    // member data: sign (chosen from -1,0,1) of each vertex,
    // location of each edge intersection, cellCenter, normal and gradNormal
    //int degreeP = m_degreeP;
    int degreeP = m_order + m_degreeP;
    int orderP  = 0;
    IndexTM<Real, SpaceDim> vectdx;
    vectdx.setAll(a_dx);
    IFData<SpaceDim> edgeData(*m_baseIF,vectdx,cellCenter,  degreeP);

    //construct data holders for all moments
    MinimalCCCM<SpaceDim> computeThisVof(edgeData);

    //compute the moments and save answers in thisVof
    computeThisVof.computeMoments(orderP,degreeP);

    CutCellMoments<SpaceDim> thisVof = computeThisVof.m_cutCellMoments;

    Vector<int> loArc[SpaceDim];
    Vector<int> hiArc[SpaceDim];

    a_node.m_cell        = a_iv;
    a_node.m_cellIndex = 0;

    //now the arcs (0 is lo, 1 is high)
    fillArc(loArc, thisVof, 0, a_ivsIrreg, a_iv);
    fillArc(hiArc, thisVof, 1, a_ivsIrreg, a_iv);
    for(int idir = 0; idir < SpaceDim; idir++)
    {
      int indexlo = IrregNode::index(idir, Side::Lo);
      int indexhi = IrregNode::index(idir, Side::Hi);
      a_node.m_arc[indexlo] = loArc[idir];
      a_node.m_arc[indexhi] = hiArc[idir];
    }

    a_node.m_volumeMoments = thisVof.m_moments;          
    a_node.m_EBMoments     = thisVof.m_EBmoments;

    IndexTM<Real,CH_SPACEDIM> point;
    for (int idir = 0; idir < SpaceDim; ++idir)
    {
      Real cellCent = (a_iv[idir]+ 0.5)*a_dx;
      point[idir] = cellCent;
    }  

    NormalDerivativeNew<SpaceDim> normalDerivative;

    int maxOrder = CH_EBIS_ORDER;

    IndexTM<Real, SpaceDim> itmpoint;
    EBArith::convertToITM(itmpoint, point);
    IFSlicer<SpaceDim> ifSlicer(*m_baseIF);
    NormalDerivativeNew<SpaceDim>::NormalDerivativeMap ndMap = 
      normalDerivative.calculateAll(maxOrder,
                                    itmpoint,
                                    &ifSlicer);

    for (int idir = 0; idir < SpaceDim; ++idir)
    {
      MomItSpaceDim momit;
      for(momit.reset(); momit.ok(); ++momit)
      {
        Real derivVal = ndMap[momit()][idir];
        a_node.m_normalPartialDeriv[idir][momit()] = derivVal;
      }
    }

  
    Iv2 bdId;
    for(SideIterator sit; sit.ok(); ++sit)
    {
      int hilo = 0;
      if(sit() == Side::Hi) hilo = 1;
      
      for (int idir = 0; idir < SpaceDim; ++idir)
      {
        int iindex = IrregNode::index(idir, sit());
        CH_assert((a_node.m_arc[iindex].size() == 1) || (a_node.m_arc[iindex].size() == 0));
        if(a_node.m_arc[iindex].size() == 1)
        {
          bdId[BDID_HILO] = hilo;
          bdId[BDID_DIR] = idir;
          const CutCellMoments<SpaceDim-1>& bdccm = thisVof.getBdCutCellMoments(bdId);
          a_node.m_faceMoments[iindex] = bdccm.m_moments;
        }
        else
        {
          a_node.m_faceMoments[iindex].setToZero();
        }
      }
    }
    //this fills volfrac, areafrac, centroids...
    a_node.setNormalizedStuff(a_dx);
  }


//records connectivity between vofs
  void 
  WrappedGShop::
  fillArc(Vector<int>                          a_arc[SpaceDim],
          CutCellMoments<SpaceDim>       &     a_cutCellMoments,
          const int                      &     a_hilo,
          const IntVectSet               &     a_ivsIrreg,
          const IntVect                  &     a_curriv) const
  {
    Iv2 bdId;
    //a_hilo is 0 or 1
    bdId[BDID_HILO] = a_hilo;
    for (int idir = 0; idir < SpaceDim; ++idir)
    {
      bdId[BDID_DIR] = idir;
      bool covered =  a_cutCellMoments.getBdCutCellMoments(bdId).isCovered();

      if (covered)
      {
        a_arc[idir].resize(0);
      }
      else
      {
        a_arc[idir].resize(1);

        //otherIV is the iv in the idir direction on the a_hilo side
        IntVect otherIV = a_curriv;
        otherIV[idir] += (a_hilo*2) - 1;

        if (m_domain.contains(otherIV))
        {
          int otherCellIndex;
          if (a_ivsIrreg.contains(otherIV))
          {
            otherCellIndex = 0;
          }
          else
          {
            //arc to regular cell
            otherCellIndex = -2;
          }
          a_arc[idir][0]=otherCellIndex;
        }
        else if (!m_domain.contains(otherIV))
        {
          //boundary arcs always -1
          a_arc[idir][0] = -1;
        }
      }
    }
  }

}
