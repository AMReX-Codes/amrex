#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include "AMReX_CoveredSlabs.H"
#include "AMReX_RealVect.H"
#include "AMReX.H"
#include "AMReX_Print.H"
#include "AMReX_IntVectSet.H"
#include "AMReX_BoxIterator.H"


namespace amrex
{
    /**
       Return true if every cell in region is regular at the
       refinement described by dx.
    */
  bool 
  CoveredSlabs::
  isRegular(const Box&           a_region,
            const Box&           a_domain,
            const RealVect&      a_origin,
            const Real&          a_dx) const
  {
    bool foundCoveredBits = false;
    for(int ibox = 0; ibox < m_coveredBoxes.size(); ibox++)
    {
      //need to grow because at the border it is irregular
      Box interBox = grow(m_coveredBoxes[ibox], 1);
      interBox &= a_region;
      if(!interBox.isEmpty())
      {
        foundCoveredBits = true;
      }
    }
    return !foundCoveredBits;
  }

  ///
  /**
     Return true if every cell in region is covered at the
     refinement described by dx.
  */
  bool 
  CoveredSlabs::
  isCovered(const Box&           a_region,
            const Box     &      a_domain,
            const RealVect&      a_origin,
            const Real&          a_dx) const
  {
    IntVectSet ivsUncovered(a_region);
    for(int ibox = 0; ibox < m_coveredBoxes.size(); ibox++)
    {
      ivsUncovered -= m_coveredBoxes[ibox];
    }
    return ivsUncovered.isEmpty();
  }
  ////////////
  void
  CoveredSlabs::
  addIrregularNode(Vector<IrregNode>        & a_nodes,
                   const IntVect            & a_iv,
                   const Box                & a_domain) const
  {
    IrregNode newNode;
    newNode.m_hasMoments = false;
    newNode.m_cell = a_iv;
    newNode.m_cellIndex = 0;
    newNode.m_volFrac = 1.0;
    newNode.m_bndryCentroid = RealVect::Zero;
    newNode.m_volCentroid   = RealVect::Zero;
    for(int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      for(SideIterator sit; sit.ok(); ++sit)
      {
        int arcInd = newNode.index(faceDir, sit());
        IntVect otherIV= a_iv + sign(sit())*BASISV(faceDir);
        bool faceThere = !(cellCovered(otherIV));
        if(faceThere)
        {
          if(a_domain.contains(otherIV))
          {
            newNode.m_arc         [arcInd] = Vector<int>(1, 0);
          }
          else
          {
            newNode.m_arc         [arcInd] = Vector<int>(1, -1);
          }
          newNode.m_areaFrac    [arcInd] = Vector<Real>(1,1.0);
          newNode.m_faceCentroid[arcInd] = Vector<RealVect>(1,RealVect::Zero);
          Real cutLocDir = sign(sit()) * 0.5;
          //done as increment to cover (pun intended) the rare case where both sides of a cell are covered.
          newNode.m_bndryCentroid[faceDir] += cutLocDir;
        }
      }
    }
    a_nodes.push_back(newNode);
  }

  /////
  void
  CoveredSlabs::
  fillGraph(BaseFab<int>             & a_regIrregCovered,
            Vector<IrregNode>        & a_nodes,
            NodeMap                  & a_intersections,
            const Box                & a_validRegion,
            const Box                & a_ghostRegion,
            const Box                & a_domain,
            const RealVect           & a_origin,
            const Real               & a_dx) const
  {
    assert(a_domain.contains(a_ghostRegion));
    a_nodes.resize(0);

    Box bigGhost = grow(a_ghostRegion, 1);
    bigGhost &= a_domain;
    a_regIrregCovered.resize(bigGhost, 1);
    long int numReg=0, numIrreg=0, numCov=0;

    for (BoxIterator bit(bigGhost); bit.ok(); ++bit)
    {
      const IntVect iv =bit();
      bool cellCov = cellCovered(iv);
      bool cellReg = cellRegular(iv);
      if(cellReg)
      {
        // set regular cells to 1
        a_regIrregCovered(iv, 0) = 1;
        numReg++;
      }
      else if(cellCov)
      {
        // set regular cells to 1
        a_regIrregCovered(iv, 0) = -1;
        numCov++;
      }
      else
      {
        // set irregular cells to 0
        a_regIrregCovered(iv, 0) = 0;
        numIrreg++;
      }
    }

    for (BoxIterator bit(a_ghostRegion); bit.ok(); ++bit)
    {
      const IntVect iv =bit();
      if(a_regIrregCovered(iv, 0) == 0)
      {
        addIrregularNode(a_nodes, iv, a_domain);
      }
    }

    amrex::Print() << "CoveredSlabs: num Regular   cells  = " << numReg     << "\n";
    amrex::Print() << "CoveredSlabs: num Irregular cells  = " << numIrreg   << "\n";
    amrex::Print() << "CoveredSlabs: num Covered   cells  = " << numCov     << "\n";
    amrex::Print() << "CoveredSlabs: number of nodes  = " << a_nodes.size() << "\n";
  }
}
