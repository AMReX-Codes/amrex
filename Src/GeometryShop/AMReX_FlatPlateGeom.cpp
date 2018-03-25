#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include "AMReX_FlatPlateGeom.H"
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
  FlatPlateGeom::
  isRegular(const Box&           a_region,
            const Box& a_domain,
            const RealVect&      a_origin,
            const Real&          a_dx) const
  {
    Real hiVal = a_dx*(a_region.bigEnd()  [m_normalDir]+1);
    Real loVal = a_dx*(a_region.smallEnd()[m_normalDir]  );
    bool retval = false;
    if(hiVal < m_plateLocation)
    {
      retval = true;
    }
    else if(loVal > m_plateLocation)
    {
      retval = true;
    }
    return retval;
  }

    ///
    /**
       Return true if every cell in region is covered at the
       refinement described by dx.
    */
  bool 
  FlatPlateGeom::
  isCovered(const Box&           a_region,
            const Box& a_domain,
            const RealVect&      a_origin,
            const Real&          a_dx) const
  {
    return false;
  }
  bool
  FlatPlateGeom::
  isCellCut(const IntVect            & a_iv,
            const Box                & a_domain,
            const RealVect           & a_origin,
            const Real               & a_dx) const
  {
    bool retval = true;
    IntVect ivdiff = a_iv - a_domain.smallEnd();
    for(int idir = 0; idir < SpaceDim; idir++)
    {
      Real loclo  =  (ivdiff[idir]  )*a_dx + a_origin[idir];
      Real lochi  =  (ivdiff[idir]+1)*a_dx + a_origin[idir];
      assert(lochi >= loclo);
      assert(m_plateHi[idir] >= m_plateLo[idir]);
      if(loclo > m_plateHi[idir])
      {
        retval = false;
      }
      if(lochi < m_plateLo[idir])
      {
        retval = false;
      }
    }

    return retval;
  }

  bool
  FlatPlateGeom::
  isFaceCut(Real                     & a_areaFracLo,
            Real                     & a_areaFracHi,
            const IntVect            & a_iv, 
            const int                & a_faceDir,
            const Side::LoHiSide     & a_sd, 
            const Box                & a_domain,
            const RealVect           & a_origin,
            const Real               & a_dx) const
  {
    if(a_faceDir == m_normalDir)
      return false;

    //
    bool retval = false;
    IntVect ivdiff = a_iv - a_domain.smallEnd();
    
    Real loclo  =  (ivdiff[m_normalDir]  )*a_dx + a_origin[m_normalDir];
    Real lochi  =  (ivdiff[m_normalDir]+1)*a_dx + a_origin[m_normalDir];
    if((loclo < m_plateLocation) && (lochi >= m_plateLocation))
    {
      Real planeloc = (ivdiff[a_faceDir])*a_dx + a_origin[m_normalDir];
      if(a_sd == Side::Hi) 
      {
        planeloc += a_dx;
      }
      if((planeloc > m_plateLo[a_faceDir]) && (planeloc < m_plateHi[a_faceDir]))
      {
        retval = true;
        a_areaFracLo = (m_plateLocation - loclo)/a_dx;
        a_areaFracLo = std::max(a_areaFracLo, 0.0);
        a_areaFracLo = std::min(a_areaFracLo, 1.0);
        a_areaFracHi = 1.0 - a_areaFracLo;
      }
    }
    return retval;
  }
  /*********************************************/
  int
  FlatPlateGeom::
  getNumVolumes(const IntVect  & a_iv, 
                const Box      & a_domain, 
                const RealVect & a_origin, 
                const Real     & a_dx) const
  {
    //if the cell gets cut (SD-1)*2 times, you get two volumes.
    //otherwise, you get one.
    //these are the area fractions of the faces on the low and high side of the plane. 
    //they will be the same if both faces in particular direction are cut
      
    int numCuts = 0;
    for(int idir = 0; idir < SpaceDim; idir++)
    {
      if(idir != m_normalDir)
      {
        for(SideIterator sit; sit.ok(); ++sit)
        {
          //not used yet.
          Real areaFracLo, areaFracHi;
          bool isCut = isFaceCut(areaFracLo, areaFracHi, a_iv, idir, sit(), a_domain, a_origin, a_dx);
          if(isCut)
          {
            numCuts++;
          }
        }
      }
    }
    int retval;
    int numCutsFull = 2*(SpaceDim-1);
    if(numCuts == 0)
    {
      retval = -1; //magic number says not really cut
    }
    else if(numCuts == numCutsFull)
    {
      //cut all the way around--two volumes
      retval = 2;
    }
    else
    {
      retval = 1;
    }

    return retval;
  }
  ////////////
  void
  FlatPlateGeom::
  addIrregularNodes(Vector<IrregNode>        & a_nodes,
                    const BaseFab<int>       & a_numVolumes,
                    const IntVect            & a_iv,
                    const Box                & a_domain,
                    const RealVect           & a_origin,
                    const Real               & a_dx) const
  {
    if(a_numVolumes(a_iv, 0) == 2)
    {
      IrregNode loNode, hiNode;
      Real areaFracLo, areaFracHi;
      loNode.m_hasMoments = false;
      hiNode.m_hasMoments = false;
      loNode.m_cell = a_iv;
      hiNode.m_cell = a_iv;
      loNode.m_cellIndex = 0;
      hiNode.m_cellIndex = 1;
      for(int faceDir = 0; faceDir < SpaceDim; faceDir++)
      {
        if(faceDir != m_normalDir)
        {
          for(SideIterator sit; sit.ok(); ++sit)
          {
            int arcInd = loNode.index(faceDir, sit());
            isFaceCut(areaFracLo, areaFracHi, a_iv, faceDir, sit(), a_domain, a_origin, a_dx);
            IntVect otherIV= a_iv + sign(sit())*BASISV(faceDir);
            int numVolumesOtherCell = -1;
            //low points at low, high points at high
            //unless outside boundary so sigh
            int otherVolLo = -1;
            int otherVolHi = -1;
            if(a_domain.contains(otherIV))
            {
              numVolumesOtherCell = a_numVolumes(otherIV, 0);
              if(numVolumesOtherCell == 2)
              {
                otherVolLo = 0;
                otherVolHi = 1;
              }
              else
              {
                otherVolLo = 0;
                otherVolHi = 0;
              }
            }
            RealVect bndryCentroid  = RealVect::Zero;
            RealVect faceCentroidLo = RealVect::Zero;
            RealVect faceCentroidHi = RealVect::Zero;
            RealVect volCentroidLo  = RealVect::Zero;
            RealVect volCentroidHi  = RealVect::Zero;
            bndryCentroid [m_normalDir] = (areaFracLo- 1.0)/2.;
            faceCentroidLo[m_normalDir] = (areaFracLo- 1.0)/2.;
            faceCentroidHi[m_normalDir] = (1.0 - areaFracHi)/2.;
            volCentroidLo [m_normalDir] = (bndryCentroid[m_normalDir]-1.0)/2.;
            volCentroidHi [m_normalDir] = (1.0-bndryCentroid[m_normalDir])/2.;
            loNode.m_volFrac = areaFracLo;
            hiNode.m_volFrac = areaFracHi; 
            loNode.m_bndryCentroid = bndryCentroid;
            hiNode.m_bndryCentroid = bndryCentroid; 
            loNode.m_volCentroid   = volCentroidLo;
            hiNode.m_volCentroid   = volCentroidHi;

            loNode.m_arc         [arcInd] = Vector<int>(1,otherVolLo);
            hiNode.m_arc         [arcInd] = Vector<int>(1,otherVolHi);
            loNode.m_areaFrac    [arcInd] = Vector<Real>(1,areaFracLo);
            hiNode.m_areaFrac    [arcInd] = Vector<Real>(1,areaFracHi);
            loNode.m_faceCentroid[arcInd] = Vector<RealVect>(1,faceCentroidLo);
            hiNode.m_faceCentroid[arcInd] = Vector<RealVect>(1,faceCentroidHi);
          } // end loop over sides
        } //end if faceDir != normaldir
        else
        {
          int arcIndLo = loNode.index(faceDir, Side::Lo);
          int arcIndHi = loNode.index(faceDir, Side::Hi);
          //here there is only one face -- the low one has a high face and so on.
          //the faces point to regular cells
          loNode.m_arc[     arcIndHi].resize(0);
          loNode.m_areaFrac[arcIndHi].resize(0);
          hiNode.m_arc[     arcIndLo].resize(0);
          hiNode.m_areaFrac[arcIndLo].resize(0);

          loNode.m_arc[         arcIndLo] = Vector<int>(1,0);
          hiNode.m_arc[         arcIndHi] = Vector<int>(1,0);
          loNode.m_areaFrac[    arcIndLo] = Vector<Real>(1,1.0);
          hiNode.m_areaFrac[    arcIndHi] = Vector<Real>(1,1.0);
          loNode.m_faceCentroid[arcIndLo] = Vector<RealVect>(1,RealVect::Zero);
          hiNode.m_faceCentroid[arcIndHi] = Vector<RealVect>(1,RealVect::Zero);
        }//end faceDir == normaldir
      } //end loop over directions
      a_nodes.push_back(loNode);
      a_nodes.push_back(hiNode);
    } //end num volumes == 2
    else //i am one volume but I point into two on at least one side
    {
      IrregNode edgeNode;
      edgeNode.m_hasMoments = false;
      edgeNode.m_cell = a_iv;
      edgeNode.m_volFrac = 1.;
      edgeNode.m_cellIndex = 0;
      edgeNode.m_volCentroid   = RealVect::Zero;
      edgeNode.m_bndryCentroid = RealVect::Zero;
        
      for(int faceDir = 0; faceDir < SpaceDim; faceDir++)
      {
        for(SideIterator sit; sit.ok(); ++sit)
        {
          int arcInd = edgeNode.index(faceDir, sit());
          Real areaFracLo, areaFracHi;
          bool isCut = isFaceCut(areaFracLo, areaFracHi, a_iv, faceDir, sit(), a_domain, a_origin, a_dx);
          if(!isCut)
          {
            edgeNode.m_arc         [arcInd] = Vector<int>(1,0);
            edgeNode.m_areaFrac    [arcInd] = Vector<Real>(1,0);
            edgeNode.m_faceCentroid[arcInd] = Vector<RealVect>(1,RealVect::Zero);
            edgeNode.m_faceCentroid[arcInd] = Vector<RealVect>(1,RealVect::Zero);
          }
          else
          {
            IntVect otherIV= a_iv + sign(sit())*BASISV(faceDir);
            int numVolumesOtherCell = -1;
            //low points at low, high points at high
            //unless outside boundary so sigh
            int otherVolLo = -1;
            int otherVolHi = -1;
            if(a_domain.contains(otherIV))
            {
              numVolumesOtherCell = a_numVolumes(otherIV, 0);
              if(numVolumesOtherCell == 2)
              {
                otherVolLo = 0;
                otherVolHi = 1;
              }
              else
              {
                otherVolLo = 0;
                otherVolHi = 0;
              }
            }
            RealVect bndryCentroid = RealVect::Zero;
            RealVect faceCentroidLo = RealVect::Zero;
            RealVect faceCentroidHi = RealVect::Zero;
            RealVect volCentroidLo = RealVect::Zero;
            RealVect volCentroidHi = RealVect::Zero;
            bndryCentroid [m_normalDir] = (areaFracLo- 1.0)/2.;
            faceCentroidLo[m_normalDir] = bndryCentroid[m_normalDir];
            faceCentroidHi[m_normalDir] = (1.0 - areaFracHi)/2.;
            volCentroidLo [m_normalDir] = (bndryCentroid[m_normalDir]-1.0)/2.;
            volCentroidHi [m_normalDir] = (1.0-bndryCentroid[m_normalDir])/2.;

            edgeNode.m_arc         [arcInd].resize(2);
            edgeNode.m_areaFrac    [arcInd].resize(2);
            edgeNode.m_faceCentroid[arcInd].resize(2);;
            edgeNode.m_arc         [arcInd][0] =      otherVolLo;
            edgeNode.m_areaFrac    [arcInd][0] =      areaFracLo;
            edgeNode.m_faceCentroid[arcInd][0] =  faceCentroidLo;
            edgeNode.m_arc         [arcInd][1] =      otherVolHi;
            edgeNode.m_areaFrac    [arcInd][1] =      areaFracHi;
            edgeNode.m_faceCentroid[arcInd][1] =  faceCentroidHi;

          }//end if face is cut
        }//loop over sides
      } //end loop over face directions
      a_nodes.push_back(edgeNode);
    }
  }
  /////
  void
  FlatPlateGeom::
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
    BaseFab<int> numVolumes( bigGhost, 1);
    a_regIrregCovered.resize(bigGhost, 1);
    long int numReg=0, numIrreg=0;

    for (BoxIterator bit(bigGhost); bit.ok(); ++bit)
    {
      const IntVect iv =bit();
      bool isCut = isCellCut(iv, a_domain, a_origin, a_dx);
      if(!isCut)
      {
        // set regular cells to 1
        a_regIrregCovered(iv, 0) = 1;
        numVolumes(iv, 0) = 1;
        numReg++;
      }
      else
      {
        // set irregular cells to 0
        int numVol = getNumVolumes(iv, a_domain, a_origin, a_dx);
        //if more careful cutcell calc gives no cuts, this is set to -1 to clue me in to set to regular
        if(numVol < 1)
        {
          a_regIrregCovered(iv, 0) = 1;
          numVolumes(iv, 0) = 1;
        }
        else
        {
          a_regIrregCovered(iv, 0) = 0;
          numVolumes(iv, 0) =  numVol;
          // TODO: Add intersections here
          // ...this will break the node map data structure, since there is currently only one intersection
          // allowed per edge.
          numIrreg++;
        }

      }
    }

    for (BoxIterator bit(a_ghostRegion); bit.ok(); ++bit)
    {
      const IntVect iv =bit();
      if(a_regIrregCovered(iv, 0) == 0)
      {
        addIrregularNodes(a_nodes, numVolumes, iv, a_domain, a_origin, a_dx);
      }
    }

    amrex::Print() << "num Regular   cells  = " << numReg     << "\n";
    amrex::Print() << "num Irregular cells  = " << numIrreg   << "\n";
    amrex::Print() << "number of nodes  = " << a_nodes.size() << "\n";
  }
}
