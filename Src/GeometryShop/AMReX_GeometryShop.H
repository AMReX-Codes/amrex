#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

//  ANAG, LBNL, DTG

#ifndef _GEOMETRYSHOP_H_
#define _GEOMETRYSHOP_H_

#include "REAL.H"
#include "RealVect.H"
#include "Box.H"
#include "IntVect.H"

#include "EBISBox.H"
#include "BaseIF.H"
#include "STLIF.H"
#include "GeometryService.H"
#include "Moments.H"

#include "NamespaceHeader.H"

///
/**
   This is the base class for the workshop algorithm.
   It forms the interface between the workshop classes
   and the geometryservice class.
 */
class GeometryShop: public GeometryService
{
public:
  ///
  /**
     Define the workshop using the local geometry description
  */
  GeometryShop(const BaseIF& a_localGeom,
               int           a_verbosity,
               RealVect      a_vectDx,
               Real          a_thrshdVoF = 1.0e-16);

  ///
  ~GeometryShop();

  ///
  bool twoEdgeIntersections(edgeMo a_edges[4]) const;

  /**
      Return true if every cell in region is regular at the
      refinement described by dx.
  */
  bool isRegular(const Box&           a_region,
                 const ProblemDomain& a_domain,
                 const RealVect&      a_origin,
                 const Real&          a_dx) const;

  ///
  /**
      Return true if every cell in region is covered at the
      refinement described by dx.
  */
  bool isCovered(const Box&           a_region,
                 const ProblemDomain& a_domain,
                 const RealVect&      a_origin,
                 const Real&          a_dx) const;


  virtual bool isIrregular(const Box&           a_region,
                           const ProblemDomain& a_domain,
                           const RealVect&      a_origin,
                           const Real&          a_dx) const ;

  virtual InOut InsideOutside(const Box&           a_region,
                              const ProblemDomain& a_domain,
                              const RealVect&      a_origin,
                              const Real&          a_dx) const ;


  virtual bool canGenerateMultiCells() const
  {
    return false;
  }

  ///
  /**
     Define the internals of the input ebisRegion.
  */
  virtual void fillGraph(BaseFab<int>&        a_regIrregCovered,
                         Vector<IrregNode>&   a_nodes,
                         const Box&           a_validRegion,
                         const Box&           a_ghostRegion,
                         const ProblemDomain& a_domain,
                         const RealVect&      a_origin,
                         const Real&          a_dx) const;

  /**
  */
  void computeVoFInternals(Real&                a_volFrac,
                           Vector<int>          a_loArc[SpaceDim],
                           Vector<int>          a_hiArc[SpaceDim],
                           Vector<Real>         a_loAreaFrac[SpaceDim],
                           Vector<Real>         a_hiAreaFrac[SpaceDim],
                           Real&                a_bndryArea,
                           RealVect&            a_normal,
                           RealVect&            a_volCentroid,
                           RealVect&            a_bndryCentroid,
                           Vector<RealVect>     a_loFaceCentroid[SpaceDim],
                           Vector<RealVect>     a_hiFaceCentroid[SpaceDim],
                           const BaseFab<int>&  a_regIrregCovered,
                           const IntVectSet&    a_ivsIrreg,
                           const VolIndex&      a_vof,
                           const ProblemDomain& a_domain,
                           const RealVect&      a_origin,
                           const Real&          a_dx,
                           const RealVect&      a_vectDx,
                           const IntVect&       a_iv) const;

  void
  getFullNodeNextToCovered(IrregNode            & a_newNode, 
                           const BaseFab<int>   & a_regIrregCovered,
                           const IntVect        & a_iv,
                           const ProblemDomain  & a_domain) const;

  int m_phase;

private:
  int  m_numCellsClipped;
  int  m_verbosity;
  Real m_threshold;
  Real m_thrshdVoF; //CP, threshold to remove very small VoFs.

  RealVect m_vectDx;

  const BaseIF* m_implicitFunction;
  const STLIF* m_stlIF;
  mutable bool m_STLBoxSet;

  static bool s_verbose;

  // local geometry description
  // const BaseLocalGeometry* m_localGeomPtr;

  /**
      Return true if every cell in region is regular at the
      refinement described by dx.
  */
  bool isRegularEveryPoint(const Box&           a_region,
                           const ProblemDomain& a_domain,
                           const RealVect&      a_origin,
                           const Real&          a_dx) const;

  ///
  /**
      Return true if every cell in region is covered at the
      refinement described by dx.
  */
  bool isCoveredEveryPoint(const Box&           a_region,
                           const ProblemDomain& a_domain,
                           const RealVect&      a_origin,
                           const Real&          a_dx) const;


  virtual bool isIrregularEveryPoint(const Box&           a_region,
                                     const ProblemDomain& a_domain,
                                     const RealVect&      a_origin,
                                     const Real&          a_dx,
                                     const Real&          a_originVal) const ;

  void edgeData3D(edgeMo               a_edges[4],
                  bool&                a_faceCovered,
                  bool&                a_faceRegular,
                  bool&                a_faceDontKnow,
                  const int            a_hiLoFace,
                  const int            a_faceNormal,
                  const Real&          a_dx,
                  const RealVect&      a_vectDx,
                  const IntVect&       a_coord,
                  const ProblemDomain& a_domain,
                  const RealVect&      a_origin) const;

  void edgeData2D(edgeMo               a_edges[4],
                  bool&                a_faceCovered,
                  bool&                a_faceRegular,
                  bool&                a_faceDontKnow,
                  const Real&          a_dx,
                  const RealVect&      a_vectDx,
                  const IntVect&       a_coord,
                  const ProblemDomain& a_domain,
                  const RealVect&      a_origin) const;

  void edgeType(bool& a_regular,
                bool& a_covered,
                bool& a_dontKnow,
                Real& a_signHi,
                Real& a_signLo) const;

  Real BrentRootFinder(const RealVect&      a_x1,
                       const RealVect&      a_x2,
                       const int&           a_range) const;

#if 1
  Real PrismoidalAreaCalc(RealVect& a_xVec,
                          RealVect& a_yVec) const;
#endif

  int getNumCellsClipped();

  Real Min(const Real x, const Real y) const;

  //stuff disallowed for all the usual reasons.
  GeometryShop()
  {
    MayDay::Abort("GeometryShop uses strong construction only");
  }
  GeometryShop(const GeometryShop& a_workshopin)
  {
    MayDay::Abort("GeometryShop disallows copy contruction");
  }
  void operator=(const GeometryShop& a_workshopin)
  {
    MayDay::Abort("GeometryShop disallows the assignment operator");
  }

};
#include "NamespaceFooter.H"
#endif
