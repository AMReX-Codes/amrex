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


#ifndef _GEOMETRYSHOP_H_
#define _GEOMETRYSHOP_H_

#include "AMReX_REAL.H"
#include "AMReX_RealVect.H"
#include "AMReX_Box.H"
#include "AMReX_IntVect.H"
#include "AMReX_BaseIF.H"
#include "AMReX_BaseFab.H"
#include "AMReX_Moments.H"
#include "AMReX_IrregNode.H"
#include "AMReX_STLIF.H"
#include "AMReX_GeometryService.H"
#include <array>

namespace amrex
{
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
                 int           a_verbosity = 0,
                 Real          a_thrshdVoF = 1.0e-16);

    ///
    ~GeometryShop() {}



    ///
    /**
       Define the internals of the input ebisRegion.
    */
    virtual void fillGraph(BaseFab<int>&        a_regIrregCovered,
                           Vector<IrregNode>&   a_nodes,
                           NodeMap&             a_intersects,
                           const Box&           a_validRegion,
                           const Box&           a_ghostRegion,
                           const Box&           a_domain,
                           const RealVect&      a_origin,
                           const Real&          a_dx) const override;

    virtual GeometryService::InOut InsideOutside(const Box&           a_region,
                                                 const Box&           a_domain,
                                                 const RealVect&      a_origin,
                                                 const Real&          a_dx) const override;

    virtual bool pointOutside(const RealVect& a_pt) const override;

    /**
     */
    void computeVoFInternals(Real&                a_volFrac,
                             std::array<Vector<int>,SpaceDim>& a_loArc,
                             std::array<Vector<int>,SpaceDim>& a_hiArc,
                             std::array<Vector<Real>,SpaceDim>& a_loAreaFrac,
                             std::array<Vector<Real>,SpaceDim>& a_hiAreaFrac,
                             Real&                a_bndryArea,
                             RealVect&            a_normal,
                             RealVect&            a_volCentroid,
                             RealVect&            a_bndryCentroid,
                             std::array<Vector<RealVect>,SpaceDim>& a_loFaceCentroid,
                             std::array<Vector<RealVect>,SpaceDim>& a_hiFaceCentroid,
                             NodeMap&             a_intersects,
                             const BaseFab<int>&  a_regIrregCovered,
                             const Box&           a_domain,
                             const RealVect&      a_origin,
                             const Real&          a_dx,
                             const IntVect&       a_iv) const;

    void
    getFullNodeWithCoveredFace(IrregNode            & a_newNode, 
                               const BaseFab<int>   & a_regIrregCovered,
                               const IntVect        & a_iv,
                               const Box            & a_domain) const;

    void
    fixRegularCellsNextToCovered(Vector<IrregNode>     & a_nodes, 
                                 BaseFab<int>         & a_regIrregCovered,
                                 const Box            & a_validRegion,
                                 const Box            & a_domain,
                                 const IntVect        & a_iv,
                                 const Real           & a_dx) const;

    ///
    /**
       Return true if every cell in region is regular at the
       refinement described by dx.
    */
    virtual bool isRegular(const Box&           a_region,
                           const Box&           a_domain,
                           const RealVect&      a_origin,
                           const Real&          a_dx) const override
      {
        return isRegularEveryPoint(a_region, a_domain, a_origin, a_dx);
      }
      
    ///
    /**
       Return true if every cell in region is covered at the
       refinement described by dx.
    */
    virtual bool isCovered(const Box&           a_region,
                           const Box&           a_domain,
                           const RealVect&      a_origin,
                           const Real&          a_dx) const override
      {
        return isCoveredEveryPoint(a_region, a_domain, a_origin, a_dx);
      }

    ///get the implicit function.
    /**
       If you have a geometry sevice that does not have an implicit function, it is a runtime error.
     */
    virtual std::shared_ptr<const BaseIF> getImplicitFunction() const
      {
        return m_implicitFunction;
      }
      
    ///True if there is an underlying implicit function
    virtual bool hasImplicitFunction() const
      {
        return true;
      }
  private:
    int  m_verbosity;
    Real m_threshold;
    Real m_thrshdVoF; //CP, threshold to remove very small VoFs.

    //for STL-based geometry descriptions
    //this just an address.  not managed by this object (do not delete)
    const STLIF* m_stlIF;
    mutable bool m_STLBoxSet;

    std::shared_ptr<const BaseIF> m_implicitFunction;


    void edgeData3D(edgeMo               a_edges[4],
                    NodeMap&             a_intersections,
                    bool&                a_faceCovered,
                    bool&                a_faceRegular,
                    bool&                a_faceDontKnow,
                    const int            a_hiLoFace,
                    const int            a_faceNormal,
                    const Real&          a_dx,
                    const IntVect&       a_coord,
                    const Box&           a_domain,
                    const RealVect&      a_origin) const;

    void edgeData2D(edgeMo          a_edges[4],
                    NodeMap&        a_intersections,
                    bool&           a_faceCovered,
                    bool&           a_faceRegular,
                    bool&           a_faceDontKnow,
                    const Real&     a_dx,
                    const IntVect&  a_iv,
                    const Box&      a_domain,
                    const RealVect& a_origin) const;

    void edgeType(bool& a_regular,
                  bool& a_covered,
                  bool& a_dontKnow,
                  Real& a_signHi,
                  Real& a_signLo) const;

    Real BrentRootFinder(const RealVect&      a_x1,
                         const RealVect&      a_x2,
                         const int&           a_range) const;

  /**
      Return true if every cell in region is regular at the
      refinement described by dx.
  */
  bool isRegularEveryPoint(const Box&           a_region,
                           const Box&           a_domain,
                           const RealVect&      a_origin,
                           const Real&          a_dx) const;

  ///
  /**
      Return true if every cell in region is covered at the
      refinement described by dx.
  */
  bool isCoveredEveryPoint(const Box&           a_region,
                           const Box&           a_domain,
                           const RealVect&      a_origin,
                           const Real&          a_dx) const;


#if 1
    Real PrismoidalAreaCalc(RealVect& a_xVec,
                            RealVect& a_yVec) const;
#endif

    GeometryShop(const GeometryShop& a_workshopin) = delete;
    void operator=(const GeometryShop& a_workshopin) = delete;

  };
}
#endif
