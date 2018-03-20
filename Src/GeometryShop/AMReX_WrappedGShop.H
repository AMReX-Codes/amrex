
#ifndef _WRAPPEDGSHOP_H_
#define _WRAPPEDGSHOP_H_

#include "AMReX_REAL.H"
#include "AMReX_RealVect.H"
#include "AMReX_EBISBox.H"
#include "AMReX_GeometryService.H"
#include "AMReX_WrappedGShop.H"
#include "AMReX_LSquares.H"
#include "AMReX_BaseIF.H"
#include "AMReX_IrregNode.H"
#include "AMReX_CutCellMoments.H"
#include "AMReX_GeometryService.H"
#include "AMReX_IntVect.H"
#include "AMReX_IntVectSet.H"
#include "AMReX_Box.H"

using std::shared_ptr;
namespace amrex
{
///base class for special refinement
  class WGSRefinementCriterion
  {
  public:
    //default says no special refinement at any point
    virtual bool refineHere(const IrregNode& a_node, const IntVect& a_iv, const Real& a_dx) const
      {
        return false;
      }

    virtual ~WGSRefinementCriterion()
      {
      }

    WGSRefinementCriterion()
      {
      }
  };
///
/**
   This is a much simplified version of GeometryShop.
   Refinement has been taken out of the internals of computecutcellmoments and
   the constrained least squares stuff is gone alltogether.   
   This is also done for isotropic dx and there is no more stuff
   about globaldim.   
*/
  class WrappedGShop: public GeometryService
  {

  public:

    int m_phase;

    ///this does return higher order moments.
    virtual bool generatesHigherOrderMoments() const
      {
        return true; 
      }

    ///
    /**
       This class will refine a cell if it is between min and max refinements.
       If it violates bounds, that triggers refinement until max refinement.
    */
    WrappedGShop(const BaseIF  & a_baseIF,
                 int             a_verbosity = 0,
                 Real            a_thrshdVoF = 1.0e-16,
                 int             a_minNumberRefines = 0,
                 int             a_maxNumberRefines = 0);
  


    ///
    /**  
         This checks to see if moments are within sane bounds.
         If a_bindMoments = true, it enforces those bounds.  
         Only returns true if the bounds are violated by the relative tolerance.
    */
    bool checkNodeMoments(IrregNode & a_node, 
                          const Real& a_dx,
                          const bool& a_bindMoments,
                          const Real& a_tolerance) const;

    virtual const BaseIF* getBaseIFPtr() const
      {
        //used for exact geometry generation. 
        return &(*m_baseIF);
      }


    ///
    ~WrappedGShop()
      {
      }

    /**
       Return true if every cell in region is regular at the
       refinement described by dx.
    */
    virtual bool isRegular(const Box           & a_region,
                           const Box & a_domain,
                           const RealVect      & a_origin,
                           const Real          & a_dx) const;

    ///
    /**
       Return true if every cell in region is covered at the
       refinement described by dx.
    */
    virtual bool isCovered(const Box           & a_region,
                           const Box & a_domain,
                           const RealVect      & a_origin,
                           const Real          & a_dx) const;


    virtual bool isIrregular(const Box           & a_region,
                             const Box & a_domain,
                             const RealVect      & a_origin,
                             const Real          & a_dx) const;

    virtual bool canGenerateMultiCells() const
      {
        return false;
      }

    bool onBoxBoundary(const IntVect        & a_iv, 
                       const Box            & a_box,
                       const int            & a_dir,
                       const Side::LoHiSide & a_sd) const;

    ///
    /**
       Define the internals of the input ebisRegion.
    */
    virtual void fillGraph(BaseFab<int>        & a_regIrregCovered,
                           Vector<IrregNode>   & a_nodes,
                           NodeMap&             a_intersections, //not filled here as yet.
                           const Box           & a_validRegion,
                           const Box           & a_ghostRegion,
                           const Box & a_domain,
                           const RealVect      & a_origin,
                           const Real          & a_dx) const;


    /**
     */
    void
    computeVoFInternals(IrregNode                &     a_node,
                        const IntVectSet         &     a_ivsIrreg,
                        const Box      &     a_domain,
                        const RealVect           &     a_origin,
                        const Real               &     a_dx,
                        const IntVect            &     a_iv) const;



    //will bound moments and check the number of refinements
    bool needToRefine(IrregNode       & a_node, 
                      const Real      & a_dx,  
                      const int       & a_numRefSoFar)const ;

    void agglomerateMoments(IrregNode              & a_node, 
                            const Vector<IrregNode>& a_refNodes,
                            const Box              & a_refBox,
                            const Real             & a_fineDx,
                            const Real             & a_coarDx) const;

    //converts a RealVect in physical coordinates to a RealVect in coordinates relative to a cell center
    RealVect convert2RelativeCoord(const RealVect& a_rVect)const;

    //converts a IndexTM<Real,SpaceDim> in physical coordinates to a
    //RealVect in coordinates relative to a cell center
    RealVect convert2RelativeCoord(const IndexTM<Real,SpaceDim>& a_rVect)const;

    void
    fillNewNode(IrregNode                &     a_node,
                const IntVectSet         &     a_ivsIrreg,
                const Box      &     a_domain,
                const RealVect           &     a_origin,
                const Real               &     a_dx,
                const IntVect            &     a_iv) const;

    void setRefinementCriterion(const shared_ptr<WGSRefinementCriterion>& a_refCrit)
      {
        m_refCrit = a_refCrit;
      }

  private:

    shared_ptr<WGSRefinementCriterion> m_refCrit;
    int  m_order;
    int  m_verbosity;
    Real m_threshold;
    Real m_thrshdVoF; //CP, threshold to remove very small VoFs.

    int                 m_degreeP;
    int                 m_minNumberRefines;
    int                 m_maxNumberRefines;


    shared_ptr<BaseIF>       m_baseIF;
  
    static Real s_relativeTol;
    /**
       Return true if every cell in region is regular at the
       refinement described by dx.
    */
    bool isRegularEveryPoint(const Box&           a_region,
                             const Box& a_domain,
                             const RealVect&      a_origin,
                             const Real&          a_dx) const;

    ///
    /**
       Return true if every cell in region is covered at the
       refinement described by dx.
    */
    bool isCoveredEveryPoint(const Box&           a_region,
                             const Box& a_domain,
                             const RealVect&      a_origin,
                             const Real&          a_dx) const;


    void fillArc(Vector<int>                          a_arc[SpaceDim],
                 CutCellMoments<SpaceDim>       &     a_cutCellMoments,
                 const int                      &     a_hilo,
                 const IntVectSet               &     a_ivsIrreg,
                 const IntVect                  &     a_curriv,
                 const Box                      &     a_domain) const;


    void
    fixRegularCellsNextToCovered(Vector<IrregNode>   & a_nodes, 
                                 BaseFab<int>        & a_regIrregCovered,
                                 const Box           & a_validRegion,
                                 const Box           & a_domain,
                                 const IntVect       & a_iv,
                                 const Real          & a_dx) const;



    void
    getFullNodeWithCoveredFace(IrregNode            & a_newNode, 
                               const BaseFab<int>   & a_regIrregCovered,
                               const IntVect        & a_ivshift,
                               const Real           & a_dx,
                               const Box            & a_domain) const;

    //stuff disallowed for all the usual reasons.
    WrappedGShop()
      {
      }
    WrappedGShop(const WrappedGShop& a_workshopin)
      {
      }
    void operator=(const WrappedGShop& a_workshopin)
      {
      }

  };
}
#endif
