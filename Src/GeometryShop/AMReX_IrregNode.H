#ifndef _IRREGNODE_H_
#define _IRREGNODE_H_

#include "AMReX_REAL.H"
#include "AMReX_IntVect.H"
#include "AMReX_RealVect.H"
#include "AMReX_BaseFab.H"
#include "AMReX_LoHiSide.H"
#include "AMReX_IndexedMoments.H"

#include "AMReX_EB_TYPEDEFS.H"
#include <array>

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
    IrregNode()
      {
        m_hasMoments = false;
      }

    ///
    ~IrregNode(){;}

    ///return the index into the arc vectors
    static int index(int a_idir, Side::LoHiSide a_side) 
      {
        AMREX_ASSERT(a_idir >= 0 && a_idir < SpaceDim);
        return (a_side == Side::Lo) ? a_idir : a_idir+SpaceDim;
      }

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
    std::array<Vector<int>,2*SpaceDim>      m_arc;

    ///
    std::array<Vector<Real>,2*SpaceDim>     m_areaFrac;

    ///
    std::array<Vector<RealVect>,2*SpaceDim> m_faceCentroid;

    Real bndryArea() const;


    /// the moment at  the irregular face associated with the  monomial with the input exponents 
    /**
       Given VoF variables x, y, z,   p = mono(0), q = mono(1), r = mono(2),
       returns integral_over_irregular_area((x^p y^q z^r) dA)
    **/
    IndMomSpaceDim m_EBMoments;
  
    IndMomSpaceDim m_volumeMoments;
    IndMomSpaceDim m_normalPartialDeriv[SpaceDim];
    bool           m_hasMoments;
    ///face centered moments
    /// the moment at  the face associated with the  monomial with the input exponents 
    /**
       Given face variables x, y,    p = mono(0), q = mono(1)
       returns integral_over_face_area((x^p y^q ) dA)
    **/
    IndMomSDMinOne m_faceMoments[2*SpaceDim];

    ///
    void setMomentsToRegular(const Real& a_dx)
      {
        //no eb for regular
        m_EBMoments.setToZero();
        for(int idir  = 0; idir < SpaceDim; idir++)
        {
          m_normalPartialDeriv[idir].setToZero();
        }
        m_volumeMoments.setRegular(a_dx);
        for(int iface  = 0; iface < 2*SpaceDim; iface++)
        {
          m_faceMoments[iface].setRegular(a_dx);
        }
      }
        
    ///
    void setMomentsToZero()
      {
        m_volumeMoments.setToZero();
        m_EBMoments.setToZero();
        for(int idir  = 0; idir < SpaceDim; idir++)
        {
          m_normalPartialDeriv[idir].setToZero();
        }
        for(int iface  = 0; iface < 2*SpaceDim; iface++)
        {
          m_faceMoments[iface].setToZero();
        }
      }

    ///
    /**
       helper function for construction.  makes irregular cell
       that has connectivitity and values like a regular cell, this
       a person can modify as the irregular cell requires.  saves
       some coding in some cases
    */
    void makeRegular(const IntVect& iv, const Real& a_dx, const Box& a_domain);

/*******************************/
    void setNormalizedStuff(const Real& a_dx);
  };
}
#endif
