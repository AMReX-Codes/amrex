#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

//  ANAG, LBNL

#ifndef _IRREGNODE_H_
#define _IRREGNODE_H_

#include "REAL.H"
#include "IntVect.H"
#include "IntVectSet.H"
#include "VolIndex.H"
#include "Interval.H"
#include "FaceIndex.H"
#include "IntVectSet.H"
#include "RealVect.H"
#include "BaseFab.H"
#include "IndexedMoments.H"
#include "EB_TYPEDEFS.H"
#include "SPMD.H"
#include "NamespaceHeader.H"
///
/**
   Node for construction of geometric information.
*/
class IrregNode
{
public:
  ///
  IrregNode();

  ///
  ~IrregNode();

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
     For arcs to regular cells, the arc = -2
     If the arc is to the domain boundary, arc  = -1.
  */
  Vector<int>      m_arc[2*SpaceDim];

  ///
  Vector<Real>     m_areaFrac[2*SpaceDim];

  ///
  Vector<RealVect> m_faceCentroid[2*SpaceDim];

  ///
  /**
      helper function for construction.  makes irregular cell
      that has connectivitity and values like a regular cell, this
      a person can modify as the irregular cell requires.  saves
      some coding in some cases
  */
  void makeRegular(const IntVect& iv, const Real& a_dx);
  void faceReserve(int location, int size);



  ///for covered cells
  void setMomentsToZero();

  //for regular cells
  void setMomentsToRegular(const Real& a_dx);

  ///data for irregular nodes

  //volume centered data

  /// the moment at  the VoF associated with the  monomial with the input exponents
  /**
     Given VoF variables x, y, z,   p = mono(0), q = mono(1), r = mono(2),
     returns integral_over_VoF(x^p y^q z^r dV)
  **/
  IndMomSpaceDim m_volumeMoments;

  /// the normal*moment at  the irregular face associated with the  monomial with the input exponents 
  /**
     Given VoF variables x, y, z,   p = mono(0), q = mono(1), r = mono(2),
     returns integral_over_irregular_area((x^p y^q z^r)*n_i dA)
  **/
  //NOT FILLED YET
#if 0
  IndMomSpaceDim m_EBNormalMoments[SpaceDim];
#endif

  /// the moment at  the irregular face associated with the  monomial with the input exponents 
  /**
     Given VoF variables x, y, z,   p = mono(0), q = mono(1), r = mono(2),
     returns integral_over_irregular_area((x^p y^q z^r) dA)
  **/
  IndMomSpaceDim m_EBMoments;
  
  IndMomSpaceDim m_normalPartialDeriv[SpaceDim];
  ///face centered moments
  /// the moment at  the face associated with the  monomial with the input exponents 
  /**
     Given face variables x, y,    p = mono(0), q = mono(1)
     returns integral_over_face_area((x^p y^q ) dA)
  **/
  IndMomSDMinOne m_faceMoments[2*SpaceDim];

  void spout()
  {
    pout() << "irreg node" << endl;
    pout() << "volume moments " << endl;
    m_volumeMoments.spout();
  
    pout() << "EB moments " << endl;
    m_EBMoments.spout();
  
    for(int iindex = 0; iindex <2*SpaceDim; iindex++)
      {
        pout() << "face moments index " << iindex << endl;
        m_faceMoments[iindex].spout();
      }
  }

  ///once the moments are set, this sets volfrac, areafrac, centroid...
  void setNormalizedStuff(const Real& a_dx);
};


std::ostream& operator<< (std::ostream&  a_os,
                          const IrregNode& a_iv);

#include "NamespaceFooter.H"
#endif
