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

#ifndef _EBISLEVEL_H_
#define _EBISLEVEL_H_

#include <map>

#include "REAL.H"
#include "IntVect.H"
#include "IntVectSet.H"
#include "IntVectSet.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "CH_HDF5.H"

#include "VolIndex.H"
#include "FaceIndex.H"
#include "EBGraph.H"
#include "EBData.H"
#include "EBISBox.H"
#include "EBISLayout.H"
#include "EBIndexSpace.H"
#include "GeometryService.H"

#include "NamespaceHeader.H"

typedef std::map<DisjointBoxLayout, EBISLayout> dmap;

class EBIndexSpace;

///singleton version of EBIndexSpace.
class Chombo_EBIS
{
private:
  ///
  /**
   This is NULL until the first time
     instance() is called.
     On the first time it is called it
     is set to  new EBIndexSpace();
     After that it points to the same
     thing.
   */
  static EBIndexSpace* s_instance;

  Chombo_EBIS()
  {;}

  ~Chombo_EBIS()
  {;}

  ///if you want to sneak an existing EBIS into the singleton.
  static bool s_aliased;
public:
  ///
  /**
     Return a pointer to the singleton
     EBIndexSpace.  Until this is called
     once, the singleton is null.  See
     Design patterns for why this is done
     this way.
   */
  static EBIndexSpace* instance();

  ///
  static void alias(const EBIndexSpace* a_input);
};

///
/**
   EBISLevel represents the geometric information
   of the domain at a particular level of refinement.
*/
class EBISLevel
{
public:
  static int s_ebislGhost;
  void dumpDebug(const string& a_string);

  void coarsenVoFs(EBISLevel& a_fineEBIS);

  void fixFineToCoarse(EBISLevel& a_fineEBIS);

  void coarsenFaces(EBISLevel& a_fineEBIS);

  long long numVoFsOnProc() const;

  void levelStitch(EBISLevel&       a_otherPhase,
                   const EBISLevel* a_fineThisPtr,
                   const EBISLevel* a_fineOtherPtr); // MF Addition  bvs

  void reconcileIrreg(EBISLevel & a_otherPhase);

  void cellStitch(EBData&        a_ebdataCoarA,
                  EBData&        a_ebdataCoarB,
                  const EBGraph& a_ebgrapCoarA,
                  const EBGraph& a_ebgrapCoarB,
                  const EBISBox& a_ebisbxFineA,
                  const EBISBox& a_ebisbxFineB,
                  const IntVect& a_iv,
                  const int&     a_aphase,
                  const int&     a_bphase);
  EBISLevel();

  ///
  /**
     Defines the level.  Calls the geoserver
     to fill the layout with geometric information.
   */
  EBISLevel(const ProblemDomain   & a_domain,
            const RealVect        & a_origin,
            const Real            & a_dx,
            const GeometryService & a_geoserver,
            const EBIndexSpace    * const a_ebisPtr               = Chombo_EBIS::instance(),
            const bool            & a_distributedData             = false,
            const bool            & a_fixRegularNextToMultiValued = true,
            const int             & a_level = 0);

  ///
  /**
     Defines the level from a finer level
     (factor of two refinement ratio).
     Fills in all the coarse-fine information
     in both levels.
   */
  EBISLevel(EBISLevel             & a_finerLevel,
            const GeometryService & a_geoserver,
            const EBIndexSpace    * const a_ebisPtr               = Chombo_EBIS::instance(),
            const bool            & a_distributedData             = false,
            const bool            & a_fixRegularNextToMultiValued = true);

  void checkGraph() const;
#ifdef CH_USE_HDF5
  /// read from file from top level only file
  EBISLevel(HDF5Handle& handle);

  /////  this is the one used by the writeAlllevels crowd
  EBISLevel(HDF5Handle& a_handle, 
            const int& a_levelNumber);

  //writes out top level only
  void write(HDF5Handle& a_handle) const;

  ///writes the given level for writeAllLevels
  void write(HDF5Handle& a_handle,
             const int& a_levelNumber) const;
#endif

  ///
  /**
     Checks to see the vofs are in the correct cells.
     Checks to see that the faces are over the correct cells.
     Checks that volume fractions, area fractions are positive.
     Bail out with MayDay::Error if anything fails.
  */
  void sanityCheck(const EBIndexSpace* const a_ebisPtr = Chombo_EBIS::instance());

  ///
  void fillEBISLayout(EBISLayout&              a_ebis,
                      const DisjointBoxLayout& a_grids,
                      const int&               a_nghost) const;

  ~EBISLevel();

  const ProblemDomain& getDomain() const;

  const DisjointBoxLayout& getGrids() const;

  DisjointBoxLayout getIrregGrids(const ProblemDomain& a_domain) const;

  DisjointBoxLayout getFlowGrids(const ProblemDomain& a_domain) const;

  DisjointBoxLayout getCoveredGrids(const ProblemDomain& a_domain) const;

  IntVectSet irregCells() const;

  const Real& getDX() const;

  const RealVect& getOrigin() const;

  ///
  void fixRegularNextToMultiValued();

  void clearMultiBoundaries();
  void setBoundaryPhase(int phase);

  void getGraphSummary(long long & a_irrVoFs,
                       long long & a_arcs,
                       long long & a_multiVoFs,
                       long long & a_zeroVoFs,
                       long long & a_zeroVoFsArcs);

  void printGraphSummary(char const * a_prefix);

  int m_phase;

  static bool s_recursive;

private:

  ///
  int m_level;

  ///
  DisjointBoxLayout m_grids;

  ///
  ProblemDomain m_domain;

  ///
  RealVect m_origin;

  ///
  Real m_dx;

  ///
  LevelData<EBGraph> m_graph;

  ///
  LevelData<EBData>  m_data;

  Real m_tolerance; // used in Multifluid face matching algorithm.

  //mutable std::map<int, dmap> m_cache;
  mutable dmap m_cache;
  mutable int m_cacheMisses, m_cacheHits, m_cacheStale;

  void dumpCache() const;

  static Real s_tolerance;
  static bool s_verbose;


  static void
  defineGraphFromGeo(LevelData<EBGraph>             &     a_graph,
                     LayoutData<Vector<IrregNode> > &     a_allNodes,
                     const GeometryService          &     a_geoserver,
                     const DisjointBoxLayout        &     a_grids,
                     const ProblemDomain            &     a_domain,
                     const RealVect                 &     a_origin,
                     const Real                     &     a_dx,
                     const bool                     &     a_distributedData);

  //make grids for this level.
  static void makeBoxes(Vector<Box>&                a_boxes,
                        Vector<unsigned long long>& a_loads,
                        const Box&                  a_region,
                        const ProblemDomain&        a_domain,
                        const GeometryService&      a_geoserver,
                        const RealVect&             a_origin,
                        const Real&                 a_dx,
                        const int                   a_ncellmax,
                        const EBIndexSpace* const   a_ebisPtr);

  static 
  void makeLoads(Vector<unsigned long long>& a_loads,
                 Vector<Box>&                a_boxes,
                 const Box&                  a_region,
                 const ProblemDomain&        a_domain,
                 const GeometryService&      a_geoserver,
                 const RealVect&             a_origin,
                 const Real&                 a_dx,
                 const int                   a_ncellmax,
                 const EBIndexSpace* const   a_ebisPtr);

  static void makeBoxes(std::list<Box>&        a_boxes,
                        const Box&             a_region,
                        const ProblemDomain&   a_domain,
                        const GeometryService& a_geoserver,
                        const RealVect&        a_origin,
                        const Real&            a_dx,
                        const int              a_ncellmax);

  //these are disallowed
  void operator=(const EBISLevel& ebiin)
  {
    MayDay::Error("EBISLevel::operator= not allowed");
  }

  EBISLevel(const EBISLevel& ebiin)
  {
    MayDay::Error("EBISLevel::copy constructor not allowed");
  }

  bool match(const Real& a_a,
             const Real& a_b) const;

  bool match(const RealVect& a_a,
             const RealVect& a_b) const;

  void refreshCache() const;

  friend class EBIndexSpace;
};

inline bool EBISLevel::match(const Real& a_a,
                             const Real& a_b) const
{
  return (a_a < (a_b + m_tolerance)) && (a_a > (a_b - m_tolerance));
}

inline bool EBISLevel::match(const RealVect& a_a,
                             const RealVect& a_b) const
{
  return D_TERM6( (a_a[0] < (a_b[0] + m_tolerance)) &&
                  (a_a[0] > (a_b[0] - m_tolerance)),
                  &&
                  (a_a[1] < (a_b[1] + m_tolerance)) &&
                  (a_a[1] > (a_b[1] - m_tolerance)),
                  &&
                  (a_a[2] < (a_b[2] + m_tolerance)) &&
                  (a_a[2] > (a_b[2] - m_tolerance)),
                  &&
                  (a_a[3] < (a_b[3] + m_tolerance)) &&
                  (a_a[3] > (a_b[3] - m_tolerance)),
                  &&
                  (a_a[4] < (a_b[4] + m_tolerance)) &&
                  (a_a[4] > (a_b[4] - m_tolerance)),
                  &&
                  (a_a[5] < (a_b[5] + m_tolerance)) &&
                  (a_a[5] > (a_b[5] - m_tolerance)));
}

#include "NamespaceFooter.H"

#endif
