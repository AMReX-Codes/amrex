#ifndef _EBGRADIENTOP_H___
#define _EBGRADIENTOP_H___

#include "AMReX_FabArray.H"
#include "AMReX_EBLevelGrid.H"
#include "AMReX_EBCellFAB.H"
#include "AMReX_EBFluxFAB.H"
#include "AMReX_EBArith.H"
#include "AMReX_LayoutData.H"
#include "AMReX_AggStencil.H"



namespace amrex
{

///
/**
   compute cell centered gradient vector from cell centered scalar
   Includes van Leer limiting if you want it.
*/
  class GradientOp
  {
  public:
    ///
    /**
       Default constructor.  User must subsequently call define().
    */
    GradientOp()
      {
        m_isDefined = false;
      }

    ///
    ~GradientOp()
      {
      }

    ///
    /**
       Defining constructor.  Constructs a valid object.
       number of components is a maximum (sets buffer sizes)
       slowmode is for debugging
    */
    GradientOp(const EBLevelGrid   & a_eblg,
               const Real          & a_dx,
               const int           & a_numComp,
               const int           & a_ghostCellsInData,
               bool a_useLimiting = true,
               bool a_slowmode = false)
      {

        define(a_eblg, a_dx, a_numComp, a_ghostCellsInData, a_useLimiting, a_slowmode);
      }

    ///
    /**
       Defines this object.  Existing information is overriden.
    */
    void define(const EBLevelGrid   & a_eblg,
                const Real          & a_dx,
                const int           & a_numComp,
                const int           & a_ghostCellsInData,
                bool a_useLimiting = true,
                bool a_slowmode = false);

    ///
    /**
       Returns true if this object was created with the defining
       constructor or if define() has been called.
    */
    bool isDefined() const
      {
        return m_isDefined;
      }

    ///
    /**
       compute the cell centered gradient of a cell-centered phi.
       This DOES call a_phi.FillBoundary.
       This DOES NOT fill ghost cells outside the domain boundary.
       You need to do that.
    */
    void
    gradient(FabArray<EBCellFAB>      & a_gradPhi,
             const FabArray<EBCellFAB>& a_phi);

    ///get the index of grad phi given variable number and vector direction
    int gradIndex(int a_icomp, int a_ivecDir)
      {
        int retval = SpaceDim*a_icomp + a_ivecDir;
        return retval;
      }

  protected:

    void defineInternals();

    bool m_isDefined;
    bool m_useLimiting;

    void
    getStencils(VoFStencil         &   a_loSten,
                VoFStencil         &   a_hiSten,
                const VolIndex     &   a_vof,
                const int          &   a_idir,
                const EBISBox      &   a_ebis);

    //since this stuff involves walking the graph, put it into a  stencil
    LayoutData< shared_ptr< AggStencil<EBCellFAB, IrregFAB> > > m_loStencil[BL_SPACEDIM];
    LayoutData< shared_ptr< AggStencil<EBCellFAB, IrregFAB> > > m_hiStencil[BL_SPACEDIM];


    LayoutData< Vector<VoFStencil> > m_loStencilSlow[BL_SPACEDIM];
    LayoutData< Vector<VoFStencil> > m_hiStencilSlow[BL_SPACEDIM];
    LayoutData< Vector<VolIndex  > > m_irregVoFsSlow;
    //for now I will do the fiddly stuff that is pointwise using vofiterator
    LayoutData<VoFIterator> m_vofit;
    

    FabArray<IrregFAB>  m_loSlope[BL_SPACEDIM];
    FabArray<IrregFAB>  m_hiSlope[BL_SPACEDIM];

    int m_nComp;
    Real m_dx;
    EBLevelGrid m_eblg;
    int m_dataGhost;

    bool m_slowMode;
  private:
    //disallowed for all the usual reasons
    GradientOp(const GradientOp& ebcin);
    void operator=(const GradientOp& fabin);

  };
}

#endif
