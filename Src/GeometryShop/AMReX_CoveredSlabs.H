#ifndef AMREX_COVEREDSLABS_H_
#define AMREX_COVEREDSLABS_H_

#include "AMReX_REAL.H"
#include "AMReX_RealVect.H"
#include "AMReX_Box.H"
#include "AMReX_IntVect.H"
#include "AMReX_BaseIF.H"
#include "AMReX_BaseFab.H"
#include "AMReX_Moments.H"
#include "AMReX_IrregNode.H"
#include "AMReX_GeometryService.H"


namespace amrex
{
  ///
  /**
     This defines a geometry that is an infinitely thin, finite size plate
     aligned with a coordinate plane.
  */
  class CoveredSlabs: public GeometryService
  {
  public:

    ///
    /**
    */
    CoveredSlabs(const Vector<Box> & a_coveredBoxes)
      {
        m_coveredBoxes = a_coveredBoxes;
      }

    ///
    ~CoveredSlabs()
    {
    }


    /**
       Return true if every cell in region is regular at the
       refinement described by dx.
    */
    virtual bool isRegular(const Box&           a_region,
                           const Box& a_domain,
                           const RealVect&      a_origin,
                           const Real&          a_dx) const;
    ///
    /**
       Return true if every cell in region is covered at the
       refinement described by dx.
    */
    virtual bool isCovered(const Box&           a_region,
                           const Box& a_domain,
                           const RealVect&      a_origin,
                           const Real&          a_dx) const;


    bool cellCovered(const IntVect& a_iv) const
      {
        Box ptbox(a_iv, a_iv);
        return isCovered(ptbox, Box(), RealVect(), 1.0);
      }

    bool cellRegular(const IntVect& a_iv) const
      {
        Box ptbox(a_iv, a_iv);
        return isRegular(ptbox, Box(), RealVect(), 1.0);
      }
    void
    addIrregularNode(Vector<IrregNode>        & a_nodes,
                     const IntVect            & a_iv,
                     const Box                & a_domain) const;

    /**
    */
    virtual void fillGraph(BaseFab<int>&              a_regIrregCovered,
                           Vector<IrregNode>&         a_nodes,
                           NodeMap&                   a_intersections,
                           const Box&                 a_validRegion,
                           const Box&                 a_ghostRegion,
                           const Box&                 a_domain,
                           const RealVect&            a_origin,
                           const Real&                a_dx) const override;

  private:

    Vector<Box> m_coveredBoxes;

  };
}
#endif
