#ifndef _EBDIVERGENCEOP_H___
#define _EBDIVERGENCEOP_H___

#include "AMReX_FabArray.H"
#include "AMReX_EBLevelGrid.H"
#include "AMReX_EBCellFAB.H"
#include "AMReX_EBFluxFAB.H"
#include "AMReX_EBArith.H"
#include "AMReX_LayoutData.H"
#include "AMReX_AggStencil.H"
#include "AMReX_EBLevelRedist.H"
#include "AMReX_EBNormalizeByVolumeFraction.H"


namespace amrex
{

///
/**
   computes hybrid divergence of face-centered fluxes 
   (1) computes over the regular grid
   (2) interpolates fluxes to centroids
   (3) computes hybrid divergence at cut cells
   using the treb average of conservative divergences as the non-conservative divergence
   (4) redistributes mass into the operator.


   THIS IS NOT MEANT TO WORK IN CASES WHERE EMBEDDED BOUNDARIES CROSS COARSE FINE BOUNDARIES;
*/
  class DivergenceOp
  {
  public:
    ///
    /**
       Default constructor.  User must subsequently call define().
    */
    DivergenceOp()
      {
        m_isDefined = false;
      }

    ///
    ~DivergenceOp()
      {
      }

    ///
    /**
       Defining constructor.  Constructs a valid object.
       number of components is a maximum (sets buffer sizes)
    */
    DivergenceOp(const EBLevelGrid   & a_eblg,
                 const Real          & a_dx,
                 const int           & a_numComp,
                 const int           & a_ghostCellsInData,
                 bool a_multiplyFluxByArea = true,
                 int a_redistRad = 1)
      {

        define(a_eblg, a_dx, a_numComp, a_ghostCellsInData, a_multiplyFluxByArea, a_redistRad);
      }

    ///
    /**
       Defines this object.  Existing information is overriden.
    */
    void define(const EBLevelGrid   & a_eblg,
                const Real          & a_dx,
                const int           & a_numComp,
                const int           & a_ghostCellsInData,
                bool a_multiplyFluxByArea = true,
                int a_redistRad = 1);


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
       Input flux at embedded boundary is assumed to be dotted 
       with OUTWARD FACING normal.  (please remember that ebisbox.normal
       returns an INWARD-facing normal for some reason).
       trustRegDivF=true  means you do not want the fortran called
    */
    void
    hybridDivergence(FabArray<EBCellFAB>      & a_divF,
                     const FabArray<EBFluxFAB>& a_Flux,
                     int isrc, int idst, int inco,
                     bool a_trustRegDivF = false);

    ///
    /**
       Modify the weights in the stencil by multiplying by
       the inputs in a normalized way.
       If you want mass
       weighting, send in the density.
    */
    void resetWeights(const FabArray<EBCellFAB>& modifier,
                      const int& ivar)
      {
        m_eblevelRedist.resetWeights(modifier, ivar);
      }

  protected:

    void defineInternals();

    bool m_isDefined;
    //since this stuff involves walking the graph, put it into a  stencil
    LayoutData< shared_ptr< AggStencil<EBFaceFAB, EBCellFAB> > > m_openStencil[BL_SPACEDIM];
    LayoutData< shared_ptr< AggStencil<IrregFAB , EBCellFAB> > > m_bdryStencil;
    //for now I will do the fiddly stuff that is pointwise using vofiterator
    LayoutData<VoFIterator> m_vofit;

    

    FabArray<EBCellFAB> m_kappaDivergence;
    FabArray<IrregFAB>  m_massDiff;

    EBLevelRedist  m_eblevelRedist;
    bool m_multiplyFluxByArea;
    EBNormalizeByVolumeFraction  m_normalizor;
    int m_nComp;
    Real m_dx;
    int m_redistRad;
    EBLevelGrid m_eblg;
    int m_dataGhost;

  private:
    //disallowed for all the usual reasons
    DivergenceOp(const DivergenceOp& ebcin);
    void operator=(const DivergenceOp& fabin);

  };
}

#endif
