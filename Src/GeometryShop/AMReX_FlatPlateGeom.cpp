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
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include "AMReX_FlatPlateGeom.H"
#include "AMReX_RealVect.H"
#include "AMReX.H"
#include "AMReX_IntVectSet.H"
#include "AMReX_BoxIterator.H"


namespace amrex
{
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
  addIrregularNodes(std::vector<IrregNode>   & a_nodes,
                    const BaseFab<int>       & a_numVolumes,
                    const IntVect            & a_iv,
                    const Box                & a_domain,
                    const RealVect           & a_origin,
                    const Real               & a_dx) const
  {
    if(a_numVolumes(a_iv, 0) == 2)
    {
/*
      IrregNode loNode, hiNode;
      Real areaFracLo, areaFracHi;
      for(int idir = 0; idir < SpaceDim; idir++)
      bool isCut = isFaceCut(areaFracLo, areaFracHi, a_iv, idir, sit(), a_domain, a_origin, a_dx);
*/
      //STUFF STILL TO CODE HERE.
      
    }
  }

  void
  FlatPlateGeom::
  fillGraph(BaseFab<int>             & a_regIrregCovered,
            std::vector<IrregNode>   & a_nodes,
            const Box                & a_validRegion,
            const Box                & a_ghostRegion,
            const Box                & a_domain,
            const RealVect           & a_origin,
            const Real               & a_dx) const
  {
    assert(a_domain.contains(a_ghostRegion));
    a_nodes.resize(0);
    a_regIrregCovered.resize(a_ghostRegion, 1);

    BaseFab<int> numVolumes(a_ghostRegion, 1);
    long int numReg=0, numIrreg=0;

    for (BoxIterator bit(a_ghostRegion); bit.ok(); ++bit)
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
          numIrreg++;
        }

      }
    }

    for (BoxIterator bit(a_ghostRegion); bit.ok(); ++bit)
    {
      const IntVect iv =bit();
      if(a_regIrregCovered(iv, 0) == 0)
      {
        const IntVect iv =bit();
        addIrregularNodes(a_nodes, numVolumes, iv, a_domain, a_origin, a_dx);
      }
    }

    std::cout << "num Regular   cells  = " << numIrreg << std::endl;
    std::cout << "num Irregular cells  = " << numReg   << std::endl;
    std::cout << "number of nodes  = " << a_nodes.size() << std::endl;
  }
}
