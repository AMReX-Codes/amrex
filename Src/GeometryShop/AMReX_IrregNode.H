

/*
 *      .o.       ooo        ooooo ooooooooo.             ooooooo  ooooo 
 *     .888.      `88.       .888' `888   `Y88.            `8888    d8'  
 *    .8"888.      888b     d'888   888   .d88'  .ooooo.     Y888..8P    
 *   .8' `888.     8 Y88. .P  888   888ooo88P'  d88' `88b     `8888'     
 *  .88ooo8888.    8  `888'   888   888`88b.    888ooo888    .8PY888.    
 * .8'     `888.   8    Y     888   888  `88b.  888    .o   d8'  `888b   
 *o88o     o8888o o8o        o888o o888o  o888o `Y8bod8P' o888o  o88888o 
 *
 */


#ifndef _IRREGNODE_H_
#define _IRREGNODE_H_

#include "AMReX_REAL.H"
#include "AMReX_IntVect.H"
#include "AMReX_RealVect.H"
#include "AMReX_BaseFab.H"
#include "AMReX_LoHiSide.H"

namespace amrex
{
  ///
  /**
     Node for construction of geometric information.
  */
  class IrregNode
  {
  public:
    ///
    IrregNode(){;}

    ///
    ~IrregNode(){;}

    ///return the index into the arc vectors
    static int index(int a_idir, Side::LoHiSide a_side) ;

    ///
    IntVect     m_cell;
    ///data for irregular nodes

    //volume centered data
    ///
    Real         m_volFrac;

    ///
    /**
       each vof has a unique index in the cell
    */
    int          m_cellIndex;

    ///
    RealVect     m_volCentroid;

    ///
    RealVect     m_bndryCentroid;

    //face centered irregular data

    ///
    /**
       Indicies into a_nodes to show connectivity.
       If the arc is to  an irregular
       cell, the index is the unique index of the vof
       in the cell.
       For arcs to regular cells, the arc = 0
       If the arc is to the domain boundary, arc  = -1.
    */
    Array<int>      m_arc[2*SpaceDim];

    ///
    Array<Real>     m_areaFrac[2*SpaceDim];

    ///
    Array<RealVect> m_faceCentroid[2*SpaceDim];

    Real bndryArea() const
    {
      Real retval;

      Real irregArea=0.0;
      for (int idir=0; idir < SpaceDim; idir++)
        {
          int loindex = index(idir, Side::Lo);
          int hiindex = index(idir, Side::Hi);
          Real loArea = 0;   
          Real hiArea = 0;
          Array<Real> loAreas = m_areaFrac[loindex];
          Array<Real> hiAreas = m_areaFrac[hiindex];
          for(int ilo = 0; ilo < loAreas.size(); ilo++)
            {
              loArea += loAreas[ilo];
            }
          for(int ihi = 0; ihi < hiAreas.size(); ihi++)
            {
              hiArea += hiAreas[ihi];
            }
          irregArea += (hiArea-loArea)*(hiArea-loArea);
        }
      retval = sqrt(irregArea);
      return retval;;
    }
    ///
    /**
       helper function for construction.  makes irregular cell
       that has connectivitity and values like a regular cell, this
       a person can modify as the irregular cell requires.  saves
       some coding in some cases
    */
    void makeRegular(const IntVect& iv, const Real& a_dx);

  };
}
#endif
