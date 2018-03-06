#ifndef _EBLEVELREDIST_H___
#define _EBLEVELREDIST_H___

#include "AMReX_REAL.H"
#include "AMReX_FArrayBox.H"
#include "AMReX_FabArray.H"
#include "AMReX_LayoutData.H"
#include "AMReX_EBLevelGrid.H"
#include "AMReX_EBCellFAB.H"
#include "AMReX_Stencils.H"
#include "AMReX_BaseIVFAB.H"
#include "AMReX_RedistStencil.H"

namespace amrex
{
  ///
  /**
   */
  class EBLevelRedist
  {
  public:
    ///
    /**
       Default constructor.  User must subsequently call define().
    */
    EBLevelRedist();
          
    ///
    ~EBLevelRedist();
          
    ///
    /**
       Defining constructor.  Constructs a valid object.
       Equivalent to default construction followed by define().
       ebisl needs ghost cells out to 3*redistradius of the stencil.
       This sets the stencil to volume-weighted.
       Use resetWeights to set to mass weighted or whatever.
    */
    EBLevelRedist(const EBLevelGrid& a_eblg,
                  const int& a_ncomp,
                  int a_redistRad = 1);
          
    ///
    /**
       Defines this object.  Existing information is overriden.
       ebisl needs ghost cells out to 3*redistradius
    */
    void define(const EBLevelGrid& a_eblg,
                const int& a_ncomp,
                int a_redistRad = 1);
          
    ///
    /**
       Modify the weights in the stencil by multiplying by
       the inputs in a normalized way.
       If you want mass
       weighting, send in the density.
    */
    void resetWeights(const FabArray<EBCellFAB>& modifier,
                      const int& ivar);
          
    ///
    /**
       Store the input mass difference in the buffers of the class
       by INCREMENTING the buffer. \\
       dm += massDiff
    */
    void increment(const BaseIVFAB<Real>& massDiff,
                   const MFIter         & datInd,
                   int idst, int inco);
          
    ///
    /**
       Redistribute the data contained in the internal buffers. \\
       solution += weight*dm
    */
    void redistribute(FabArray<EBCellFAB>& solution,
                      int idst, int inco);
          
    ///
    /**

       version where solution is not necessarily ncomp long.
       intervals have to be the same length.
       srcVar can be up to m_ncomp long.
       dstVar can be up to solution.nComp().
    */
    void redistribute(FabArray<EBCellFAB>& a_solution,
                      int isrc, int idst, int inco);
          
    ///
    /**
       Returns true if this object was created with the defining
       constructor or if define() has been called.
    */
    bool isDefined() const;
          
    ///
    /**
       Reset the mass buffers of the class.
    */
    void setToZero();
          
  protected:
    int redistRad;
    RedistStencil m_stencil;
    int m_ncomp;
    int m_redistRad;
    EBLevelGrid m_eblg;

    bool m_isDefined;
    FabArray<BaseIVFAB<Real> > m_buffer;
    shared_ptr<LayoutData<IntVectSet > > m_sets;
    shared_ptr<LayoutData<VoFIterator> > m_vofit;
  private:
          
    //forbidden for all the usual reasons
    void operator=(const EBLevelRedist& a_stenin);
    EBLevelRedist(const EBLevelRedist& a_stenin);

  };
}

#endif
