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


namespace amrex
{
  ///
  /**
     This is the base class for the workshop algorithm.
     It forms the interface between the workshop classes
     and the geometryservice class.
  */
  class GeometryShop
  {
  public:

    enum InOut
      {
        Regular,
        Covered,
        Irregular
      };

    ///
    /**
       Define the workshop using the local geometry description
    */
    GeometryShop(const BaseIF& a_localGeom,
                 int           a_verbosity,
                 Real          a_thrshdVoF = 1.0e-16);

    ///
    ~GeometryShop();


    ///
    /**
       Define the internals of the input ebisRegion.
    */
    virtual void fillGraph(BaseFab<int>&        a_regIrregCovered,
                           std::vector<IrregNode>&   a_nodes,
                           const Box&           a_validRegion,
                           const Box&           a_ghostRegion,
                           const Box&           a_domain,
                           const RealVect&      a_origin,
                           const Real&          a_dx) const;

    GeometryShop::InOut InsideOutside(const Box&           a_region,
                                      const Box&           a_domain,
                                      const RealVect&      a_origin,
                                      const Real&          a_dx) const;
    /**
     */
    void computeVoFInternals(Real&                a_volFrac,
                             std::vector<int>          a_loArc[SpaceDim],
                             std::vector<int>          a_hiArc[SpaceDim],
                             std::vector<Real>         a_loAreaFrac[SpaceDim],
                             std::vector<Real>         a_hiAreaFrac[SpaceDim],
                             Real&                a_bndryArea,
                             RealVect&            a_normal,
                             RealVect&            a_volCentroid,
                             RealVect&            a_bndryCentroid,
                             std::vector<RealVect>     a_loFaceCentroid[SpaceDim],
                             std::vector<RealVect>     a_hiFaceCentroid[SpaceDim],
                             const BaseFab<int>&  a_regIrregCovered,
                             const Box&           a_domain,
                             const RealVect&      a_origin,
                             const Real&          a_dx,
                             const IntVect&       a_iv) const;

    void
    getFullNodeNextToCovered(IrregNode            & a_newNode, 
                             const BaseFab<int>   & a_regIrregCovered,
                             const IntVect        & a_iv,
                             const Box            & a_domain) const;

  private:
    int  m_verbosity;
    Real m_threshold;
    Real m_thrshdVoF; //CP, threshold to remove very small VoFs.


    const BaseIF* m_implicitFunction;


    void edgeData3D(edgeMo               a_edges[4],
                    bool&                a_faceCovered,
                    bool&                a_faceRegular,
                    bool&                a_faceDontKnow,
                    const int            a_hiLoFace,
                    const int            a_faceNormal,
                    const Real&          a_dx,
                    const IntVect&       a_coord,
                    const Box& a_domain,
                    const RealVect&      a_origin) const;

    void edgeData2D(edgeMo               a_edges[4],
                    bool&                a_faceCovered,
                    bool&                a_faceRegular,
                    bool&                a_faceDontKnow,
                    const Real&          a_dx,
                    const IntVect&       a_coord,
                    const Box& a_domain,
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

    Real Min(const Real x, const Real y) const;

    //stuff disallowed for all the usual reasons.
    GeometryShop()
    {
      amrex::Abort("GeometryShop uses strong construction only");
    }
    GeometryShop(const GeometryShop& a_workshopin)
    {
      amrex::Abort("GeometryShop disallows copy contruction");
    }
    void operator=(const GeometryShop& a_workshopin)
    {
      amrex::Abort("GeometryShop disallows the assignment operator");
    }

  };
}
#endif
