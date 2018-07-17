#ifndef AMREX_REDISTSTENCIL_H_
#define AMREX_REDISTSTENCIL_H_

#include "AMReX_REAL.H"
#include "AMReX_FArrayBox.H"
#include "AMReX_FabArray.H"
#include "AMReX_BoxArray.H"
#include "AMReX_EBLevelGrid.H"
#include "AMReX_LayoutData.H"
#include "AMReX_EBCellFAB.H"
#include "AMReX_Stencils.H"
#include "AMReX_BaseIVFAB.H"

namespace amrex
{

  /// stencil for mass redistribution
  /**
     The  RedistStencil class holds the stencil for mass redistribution
     at every irregular VoF in a layout.  The default weights
     that the stencil holds are volume weights.  The class
     does allow the flexibility to redefine these weights.
  */
  class RedistStencil
  {
  public:
    ///
    /**
       Default constructor.  User must subsequently call define().
    */
    RedistStencil();

    ///
    ~RedistStencil();

    ///
    /**
       Defining constructor.  Constructs a valid object.
       Equivalent to default construction followed by define().
       ebisl needs ghost cells out to 3*redistradius
    */
    RedistStencil(const EBLevelGrid & a_eblg,
                  const int         & a_redistRadius);

    ///
    /**
       Defines this object.  Existing information is overriden.
       ebisl needs ghost cells out to 3*redistradius
    */
    void define(const EBLevelGrid & a_eblg,
                const int         & a_redistRadius);

    ///
    /**
       Returns true if this object was created with the defining
       constructor or if define() has been called.
    */
    bool isDefined() const;

    ///
    /**
       Modify the weights in the stencil by multiplying by
       the inputs in a normalized way.
       If you want mass
       weighting, send in the density.
    */
    void resetWeights(const FabArray<EBCellFAB>& a_modifier,
                      const int&                 a_ivar);

    ///
    /**
       Returns the redistribution stencil at every irregular
       point in input  Box associated with  this DataIndex.
    */
    const BaseIVFAB<VoFStencil >&
    operator[](const MFIter& datInd) const;

    ///
    /**
     */
    int getRedistRadius() const;

  protected:

    LayoutData<BaseIVFAB<VoFStencil > > m_stencil;
    LayoutData<BaseIVFAB<VoFStencil > > m_volsten;

    EBLevelGrid m_eblg;

    bool m_isDefined;
    bool m_hasDefaultWeights;
    int  m_redistRadius;

  private:
    //internal use
    void computePointStencil(VoFStencil&      a_stencil,
                             const VolIndex&  a_srcVoF,
                             const MFIter& a_datInd);

    //disallowed for all the usual reasons
    RedistStencil(const RedistStencil& ebcin);

    void operator=(const RedistStencil& fabin);

  };
}

#endif
