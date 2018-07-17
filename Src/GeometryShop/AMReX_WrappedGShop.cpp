
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include "AMReX_REAL.H"
#include "AMReX_CH_EBIS_ORDER.H"
#include "AMReX_EBArith.H"
#include "AMReX_MomentIterator.H"
#include "AMReX_WrappedGShop.H"
#include "AMReX_EB_TYPEDEFS.H"


#include "AMReX_NoRefinement.H"
#include "AMReX_FixedRefinement.H"
#include "AMReX_DivNormalRefinement.H"
#include "AMReX_NormalDerivativeNew.H"
#include "AMReX_MinimalCCCM.H"
#include "AMReX_MomentIterator.H"

#include "AMReX_PolyGeom.H"
#include "AMReX_RealVect.H"


namespace amrex
{
Real WrappedGShop::s_relativeTol = 0.1;

/*********************************************/
/*********************************************/
  WrappedGShop::
  WrappedGShop(const BaseIF  &       a_baseIF,
               int                   a_verbosity,
               Real                  a_thrshdVoF,
               int                   a_minNumberRefines ,
               int                   a_maxNumberRefines )
  {
    BL_PROFILE("WrappedGShop::WrappedGShop");
    m_baseIF.reset(a_baseIF.newImplicitFunction());;
    m_verbosity = a_verbosity;
    m_thrshdVoF = a_thrshdVoF;
    m_minNumberRefines = a_minNumberRefines;
    m_maxNumberRefines = a_maxNumberRefines;

    m_thrshdVoF = a_thrshdVoF;
            
    m_order   = 0;
    m_degreeP = CH_EBIS_ORDER + 1;
    m_refCrit = shared_ptr<WGSRefinementCriterion> (new WGSRefinementCriterion());
  }

/*********************************************/
  bool 
  WrappedGShop::
  isRegular(const Box           & a_region,
            const Box & a_domain,
            const RealVect      & a_origin,
            const Real          & a_dx) const
  {
//    BL_PROFILE("WrappedGShop::isRegular");

    // first check any of the Box corners are outside, and return false
    // right away. (bvs)
    RealVect physCorner;
    IntVect lo = a_region.smallEnd();
    IntVect len = a_region.size();
    Box unitBox(IntVect::Zero, IntVect::Unit);
    for (BoxIterator bit(unitBox); bit.ok(); ++bit)
    {
      IntVect current = lo + len*bit();
      for (int idir = 0; idir < BL_SPACEDIM; ++idir)
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
//    BL_PROFILE("WrappedGShop::isRegularEveryPoint");

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

        for (int idir = 0; idir < BL_SPACEDIM; ++idir)
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
//    BL_PROFILE("WrappedGShop::isCovered");

    // first check any of the Box corners are inside, and return false
    // right away. (bvs)
    RealVect physCorner;
    IntVect lo = a_region.smallEnd();
    IntVect len = a_region.size();
    Box unitBox(IntVect::Zero, IntVect::Unit);
    for (BoxIterator bit(unitBox); bit.ok(); ++bit)
    {
      IntVect current = lo + len*bit();
      for (int idir = 0; idir < BL_SPACEDIM; ++idir)
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
//    BL_PROFILE("WrappedGShop::isCoveredEveryPoint");

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

        for (int idir = 0; idir < BL_SPACEDIM; ++idir)
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
//    BL_PROFILE("WrappedGShop::isIrregular");

    // first check if some Box corners are inside and some are outside, and return
    // true right away. (bvs)
    RealVect physCorner;
    IntVect lo = a_region.smallEnd();
    IntVect len = a_region.size();
    for (int idir = 0; idir < BL_SPACEDIM; ++idir)
    {
      physCorner[idir] = a_dx*(lo[idir]) + a_origin[idir];
    }
    Real originVal = m_baseIF->value(physCorner);

    Box unitBox(IntVect::Zero, IntVect::Unit);
    for (BoxIterator bit(unitBox); bit.ok(); ++bit)
    {
      IntVect current = lo + len*bit();
      for (int idir = 0; idir < BL_SPACEDIM; ++idir)
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
      Error("bogus side");
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
      amrex::Warning("still found bogus after agglomeration");
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
            NodeMap&            a_intersections, //not filled here as yet.
            const Box&          a_validRegion,
            const Box&          a_ghostRegion,
            const Box & a_domain,
            const RealVect      & a_origin,
            const Real          & a_dx) const
  {
//    BL_PROFILE("WrappedGShop::fillGraph");

    BL_ASSERT(a_domain.contains(a_ghostRegion));

    AMREX_ASSERT(a_domain.contains(a_ghostRegion));
    a_nodes.resize(0);
    a_regIrregCovered.resize(a_ghostRegion, 1);

    Real thrshd = m_thrshdVoF;
    IntVectSet ivsirreg;

    {
//      BL_PROFILE("boxiterator_loop");
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

    //if a regular is next to a  covered, change to irregular with correct arcs and so on.
    for (BoxIterator bit(a_ghostRegion); bit.ok(); ++bit)
      {
        const IntVect iv =bit();

      if(a_regIrregCovered(iv, 0) == -1)
        {
          fixRegularCellsNextToCovered(a_nodes, a_regIrregCovered, a_validRegion, a_domain, iv, a_dx);
        }
      }

    //now loop through irregular cells and make nodes for each  one.
    IntVectSet ivsdrop; //volumes too small to keep
    int numIrreg = 0;
    for (IVSIterator ivsit(ivsirreg); ivsit.ok(); ++ivsit)
    {
      IntVect iv = ivsit();
      numIrreg++;
//      BL_PROFILE("fillGraph::endOfirregularCellLoop");
      IrregNode newNode;

      fillNewNode(newNode,
                  ivsirreg,
                  a_domain,
                  a_origin,
                  a_dx,
                  ivsit());
      
      Real volFrac = newNode.m_volFrac;
      if (thrshd > 0. && volFrac < thrshd)
      {
        ivsdrop |= iv;
        a_regIrregCovered(iv, 0) = -1;
      }//CP record these nodes to be removed
      else
      {
        bool fixMoments = true;
        checkNodeMoments(newNode, a_dx, fixMoments, s_relativeTol);      
        a_nodes.push_back(newNode);
      }
    } //end loop over cells in the box
    // CP: fix sweep that removes cells with volFrac less than a certain threshold
    for(IVSIterator ivsit(ivsdrop); ivsit.ok(); ++ivsit)
      {
        const IntVect& iv = ivsit();
  
        for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
          {
            for (SideIterator sit; sit.ok(); ++sit)
              {
                int isign = sign(sit());
                IntVect otherIV = iv + isign*BASISV(faceDir);
                if (a_validRegion.contains(otherIV))
                  {
                    if (a_regIrregCovered(otherIV,0) == 0)
                      {
                        // i am in the case where the other cell
                        // is also irregular.   I just made a previously
                        // irregular cell covered so I have to check to
                        // see if it had any faces pointed this way.
                        int inode = -1;
                        bool found = false;
                        for (int ivec = 0; ivec < a_nodes.size() && ! found; ivec++)
                          {
                            if (a_nodes[ivec].m_cell == otherIV)
                              {
                                inode = ivec;
                                found = true;
                              }
                          }
                        if (!found && a_validRegion.contains(otherIV))
                          {
                            amrex::Abort("something wrong in our logic");
                          }
                        if (found)
                          {
                            int arcindex = a_nodes[inode].index(faceDir, flip(sit()));
                            a_nodes[inode].m_arc[         arcindex].resize(0);
                            a_nodes[inode].m_areaFrac[    arcindex].resize(0);
                            a_nodes[inode].m_faceCentroid[arcindex].resize(0);
                            a_nodes[inode].m_faceMoments [arcindex].setToZero();
                          }
                      }
                  }//valid region
              }//sit
          }//facedir
        
        //also need to fix regular cells next to new covered cell
        fixRegularCellsNextToCovered(a_nodes, a_regIrregCovered, a_validRegion, a_domain, iv, a_dx);

      }//ivsdrop
    if(m_verbosity > 2)
    {
      amrex::AllPrint() << "WrappedGShop:num irreg vofs   = " << numIrreg << "\n";
      amrex::AllPrint() << "WrappedGShop:number of nodes  = " << a_nodes.size() << "\n";
    }
  }
  /*************/
  void
  WrappedGShop::
  fixRegularCellsNextToCovered(Vector<IrregNode>   & a_nodes, 
                               BaseFab<int>        & a_regIrregCovered,
                               const Box           & a_validRegion,
                               const Box           & a_domain,
                               const IntVect       & a_iv,
                               const Real          & a_dx) const

  {
    Box grownBox(a_iv, a_iv);
    grownBox.grow(1);
    grownBox  &= a_domain;
    IntVectSet ivstocheck(grownBox);
    ivstocheck -= a_iv;
    Box ghostRegion = a_regIrregCovered.box();
    //first check neighbors in each direction.  
    //If any of these are regular, they are replaced 
    //by irregular cells with a boundary face facing the covered cell.
    for(int idir = 0; idir < SpaceDim; idir++)
      {
        for(SideIterator sit; sit.ok(); ++sit)
          {
            int ishift = sign(sit());
            IntVect ivshift = a_iv + ishift*BASISV(idir);
            ivstocheck -= ivshift;
            int bfvalshift = -1;
            if(ghostRegion.contains(ivshift))
              {
                bfvalshift = a_regIrregCovered(ivshift, 0);
              }
            if(bfvalshift  == 1)
              {
                a_regIrregCovered(ivshift, 0) =  0;

                if(a_validRegion.contains(ivshift))
                  {
                    IrregNode newNode;
                    getFullNodeWithCoveredFace(newNode, 
                                               a_regIrregCovered,
                                               ivshift, 
                                               a_dx,
                                               a_domain);
                    a_nodes.push_back(newNode);
                  }

              }
          }
      }
    //next we loop through the remaining cells (corner cells in 2d, corner and edge cells in 3D)
    //if any of these are regular, we change them to irregular 
    for(IVSIterator ivsit(ivstocheck); ivsit.ok(); ++ivsit)
      {
        const IntVect& iv = ivsit();
        if(ghostRegion.contains(iv))
          {
            if(a_regIrregCovered(iv, 0) == 1)
              {
                a_regIrregCovered(iv, 0) = 0;
                IrregNode newNode;
                newNode.makeRegular(iv, a_dx, a_domain);
                a_nodes.push_back(newNode);
              }
          }
      }
  }
  /**********************************************/
  void
  WrappedGShop::
  getFullNodeWithCoveredFace(IrregNode            & a_newNode, 
                             const BaseFab<int>   & a_regIrregCovered,
                             const IntVect        & a_iv,
                             const Real           & a_dx,
                             const Box            & a_domain) const
  {

    a_newNode.m_cell          = a_iv;
    a_newNode.m_volFrac       = 1.0;
    a_newNode.m_cellIndex     = 0;
    a_newNode.m_volCentroid   = RealVect::Zero;
    //set regular cell values then fix up
    a_newNode.m_bndryCentroid = RealVect::Zero;
    //set all moments to regular and then zero out the appropriate face mometns
    a_newNode.setMomentsToRegular(a_dx);
    int coveredDir;
    Side::LoHiSide coveredSide;
    bool found = false;

    for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
      {
        for(SideIterator sit; sit.ok(); ++sit)
          {
            int ishift = sign(sit());
            IntVect ivshift = a_iv + ishift*BASISV(faceDir);
            Vector<int> arc;
            Vector<Real> areaFrac;
            Vector<RealVect> faceCentroid;
            if(!a_domain.contains(ivshift))
              {
                // boundary arcs always -1
                arc.resize(1,-1);
                areaFrac.resize(1, 1.0);
                faceCentroid.resize(1, RealVect::Zero);
              }
            else if (a_regIrregCovered(ivshift, 0) >= 0)
              {
                //irregular cell or regular cell
                //compute vof internals returns something special if 
                //connected to a regular cell but EBGraph treats both the  same.
                //it just  knows that the cell index of a regular cell is 0
                arc.resize(1,0);
                areaFrac.resize(1, 1.0);
                faceCentroid.resize(1, RealVect::Zero);
              }
            else if (a_regIrregCovered(ivshift, 0) < 0)
              {
                found = true;
                coveredDir= faceDir;
                coveredSide = sit();
                // covered face!
                arc.resize(0);
                areaFrac.resize(0);
                faceCentroid.resize(0);
              }
            else
              {
                amrex::Error("logic error");
              }
          
            int nodeInd = a_newNode.index(faceDir, sit());
            a_newNode.m_arc[nodeInd]          = arc;
            a_newNode.m_areaFrac[nodeInd]     = areaFrac;
            a_newNode.m_faceCentroid[nodeInd] = faceCentroid;
            a_newNode.m_faceMoments[nodeInd].setToZero();
          }
      }
    //fix boundary centroid
    if(found)
      {
        int centsign = sign(coveredSide);
        a_newNode.m_bndryCentroid[coveredDir] =  centsign*0.5;
      }
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
//    BL_PROFILE("WrappedGShop::ComputeVofInternals");

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
    fillArc(loArc, thisVof, 0, a_ivsIrreg, a_iv, a_domain);
    fillArc(hiArc, thisVof, 1, a_ivsIrreg, a_iv, a_domain);
    for(int idir = 0; idir < SpaceDim; idir++)
    {
      int indexlo = IrregNode::index(idir, Side::Lo);
      int indexhi = IrregNode::index(idir, Side::Hi);
      a_node.m_arc[indexlo] = loArc[idir];
      a_node.m_arc[indexhi] = hiArc[idir];
    }

    a_node.m_volumeMoments = thisVof.m_moments;          
    a_node.m_EBMoments     = thisVof.m_EBmoments;

    RealVect point = EBArith::getIVLocation(a_iv, a_dx*RealVect::Unit, RealVect::Zero);
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
        BL_ASSERT((a_node.m_arc[iindex].size() == 1) || (a_node.m_arc[iindex].size() == 0));
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
          const IntVect                  &     a_curriv,
          const Box                      &     a_domain) const
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

        if (a_domain.contains(otherIV))
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
        else if (!a_domain.contains(otherIV))
        {
          //boundary arcs always -1
          a_arc[idir][0] = -1;
        }
      }
    }
  }

}
