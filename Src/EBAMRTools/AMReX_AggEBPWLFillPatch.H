#ifndef _AGGEBPWLFILLPATCH_H___
#define _AGGEBPWLFILLPATCH_H___

#include "AMReX_REAL.H"
#include "AMReX_FArrayBox.H"
#include "AMReX_FabArray.H"
#include "AMReX_EBISLayout.H"
#include "AMReX_EBCellFAB.H"
#include "AMReX_Stencils.H"
#include "AMReX_AggStencil.H"
#include "AMReX_EBIndexSpace.H"
#include "AMReX_EBLevelGrid.H"

namespace amrex
{
///
/**
   This is a  ghost cell filler that fills by
   piecewise linear interpolation.  The slopes used in filling
   the ghost cells are limited by delta min mod limiting.
*/
  class AggEBPWLFillPatch
  {
  public:
    ///
    /**
       Default constructor.  User must subsequently call define().
    */
    AggEBPWLFillPatch();

    ///
    virtual  ~AggEBPWLFillPatch();

    ///
    /**
       Defining constructor.  Constructs a valid object.
       Equivalent to default construction followed by define().

       {\bf Arguments:}\\
       eblgFine, eblgCoar: The fine and coarse layouts
       of the data with the geometric description.\\
       nref: The refinement ratio between the two levels. \\
       nvar: The number of variables contained in the data
       at each VoF.\\
       radius: the total number of
       ghost cells filled.
       forceNoEBCF = true promises that the embedded boundary 
       will not cross the coarse-fine interface.   This will save you  
       a  lot of compute time if you are certain you are tagging all irregulkar cells.
    */
    AggEBPWLFillPatch(const EBLevelGrid & a_eblgFine,
                      const EBLevelGrid & a_eblgCoar,
                      const int& a_nref,
                      const int& a_nvar,
                      const int& a_radius,
                      const int& a_ghostCellsInData,
                      const bool& a_forceNoEBCF = false)
      {
        setDefaultValues();
        define(a_eblgFine, a_eblgCoar,
               a_nref, a_nvar, a_radius, 
               a_ghostCellsInData, a_forceNoEBCF);
      }
                    
                    
    ///
    /**
       Defines this object.  Existing information is overriden.
       {\bf Arguments:}\\
       eblgFine, eblgCoar: The fine and coarse layouts
       of the data with the geometric description.\\
       nref: The refinement ratio between the two levels. \\
       nvar: The number of variables contained in the data
       at each VoF.\\
       radius: the total number of
       ghost cells filled.
       forceNoEBCF = true promises that the embedded boundary 
       will not cross the coarse-fine interface.   This will save you  
       a  lot of compute time if you are certain you are tagging all irregulkar cells.
    */
    void define(const EBLevelGrid & a_eblgFine,
                const EBLevelGrid & a_eblgCoar,
                const int& a_nref,
                const int& a_nvar,
                const int& a_radius,
                const int& a_ghostCellsInData,
                const bool& a_forceNoEBCF = false);


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
       Interpolate the fine data from the coarse data
       over the intersection
       of the fine layout with the refined coarse layout. \\
       {\bf Arguments:}\\
       coarDataOld: The data over the coarse layout at coarse time old.\\
       coarDatanew: The data over the coarse layout at coarse time new.\\
       fineData: The data over the fine layout.  \\
       Error occurs unless  coarTimeOld <= fineTime <= coarTimeNew \\
       Fine and coarse data must
       have the same number of variables.\\
       variables:  The variables to interpolate.  Those not
       in this range will be left alone.   This range of variables
       must be in both the coarse and fine data.
    */
    virtual void
    interpolate(FabArray<EBCellFAB>& a_fineData,
                const FabArray<EBCellFAB>& a_coarDataOld,
                const FabArray<EBCellFAB>& a_coarDataNew,
                const Real& a_coarTimeOld,
                const Real& a_coarTimeNew,
                const Real& a_fineTime,
                int idst, int inco) const;

  protected:

    //internal use only functions


    void getOffsets(const LayoutData< Vector<VolIndex> >&  a_srcVoFsCoar,
                    const LayoutData<IntVectSet>&          a_irregRegionsFine,
                    const LayoutData<Vector<VoFStencil> >  a_loStencils[SpaceDim],
                    const LayoutData<Vector<VoFStencil> >  a_hiStencils[SpaceDim],
                    const LayoutData<IntVectSet>           a_coarLoInterp[SpaceDim],
                    const LayoutData<IntVectSet>           a_coarHiInterp[SpaceDim],
                    const LayoutData<IntVectSet>           a_coarCeInterp[SpaceDim]);


    void
    interpolateFAB(EBCellFAB& a_fine,
                   const EBCellFAB& a_coarOld,
                   const EBCellFAB& a_coarNew,
                   const Real& a_coarTimeOld,
                   const Real& a_coarTimeNew,
                   const Real& a_fineTime,
                   const MFIter& a_datInd,
                   int idst, int inco) const;

    void
    makeStencils();

    void defineSlopeHolders(const LayoutData<IntVectSet>& a_irregRegionsCoFi);

    void
    getIVS(LayoutData<IntVectSet>& a_irregRegionsFine,
           LayoutData<IntVectSet>& a_irregRegionsCoFi,
           LayoutData< Vector<VolIndex> >&   a_srcVoFs);

    void
    defineAggStencils(LayoutData<Vector<VoFStencil> >  a_loStencils[SpaceDim],
                      LayoutData<Vector<VoFStencil> >  a_hiStencils[SpaceDim],
                      const LayoutData< Vector<VolIndex> >&   a_srcVoFs);
    void
    getLoHiCenIVS(LayoutData<IntVectSet>  a_coarLoInterp[SpaceDim],
                  LayoutData<IntVectSet>  a_coarHiInterp[SpaceDim],
                  LayoutData<IntVectSet>  a_coarCeInterp[SpaceDim]);

    void
    getSten(LayoutData<Vector<VoFStencil> >  a_loStencils[SpaceDim],
            LayoutData<Vector<VoFStencil> >  a_hiStencils[SpaceDim],
            LayoutData<IntVectSet>           a_coarLoInterp[SpaceDim],
            LayoutData<IntVectSet>           a_coarHiInterp[SpaceDim],
            LayoutData<IntVectSet>           a_coarCeInterp[SpaceDim],
            const LayoutData< Vector<VolIndex> >&   a_srcVoFs);

    void getSlopes(const EBCellFAB& a_coarDataOld,
                   const EBCellFAB& a_coarDataNew,
                   const MFIter& a_dit,
                   int idst, int inco) const;


    void deltaMinMod(Real& a_slope, Real& a_slopeLo, Real& a_slopeHi) const;

    void setDefaultValues();

    struct
    {
      size_t offset;
      int    dataID;
    } typedef access_t;

    struct
    {
      access_t dest_access; //offsets for fine data
      int      slop_index;  //index into vector of slope offsetsf
      IntVect  coariv;
      IntVect  fineiv;
    } typedef fine_logic_t;

    struct
    {
      access_t slop_access;      //offsets for slope data
      access_t data_access;      //offsets for coarsened fine data
      bool     has_lo[SpaceDim]; //whether there is a lo slope stencil
      bool     has_hi[SpaceDim]; //whether there is a hi slope stencil
      bool     ivs_lo[SpaceDim]; //whether inside ivs centered
      bool     ivs_hi[SpaceDim]; //whether inside ivs centered
      bool     ivs_ce[SpaceDim]; //whether inside ivs centered
    } typedef coar_logic_t;

    //data to provide fast access to fine points and their associated slopes
    LayoutData< Vector<fine_logic_t> > m_fineOffsets;
    //data so we can do the wacky delta minmod thing fast.
    LayoutData< Vector<coar_logic_t> > m_coarOffsets;

    bool m_forceNoEBCF;
    bool m_isDefined;

    //ghost cells in solution
    int  m_ghost;
    //domain at the coarse level
    Box m_coarDomain;
    //layouts
    EBLevelGrid m_eblgFine;
    EBLevelGrid m_eblgCoar;
    EBLevelGrid m_eblgCoFi;


    mutable FabArray<EBCellFAB> m_coarOnFDataOld;
    mutable FabArray<EBCellFAB> m_coarOnFDataNew;

    int m_refRat;
    int m_nComp;
    int m_radius;
    int m_coarGhostRad;

    //building this pig for speed only.   this is scratch data for slopes
    mutable LayoutData< BaseIVFAB<Real> > m_slopeLoOld[SpaceDim];
    mutable LayoutData< BaseIVFAB<Real> > m_slopeHiOld[SpaceDim];
    mutable LayoutData< BaseIVFAB<Real> > m_slopeCeOld[SpaceDim];
    mutable LayoutData< BaseIVFAB<Real> > m_slopeLoNew[SpaceDim];
    mutable LayoutData< BaseIVFAB<Real> > m_slopeHiNew[SpaceDim];
    mutable LayoutData< BaseIVFAB<Real> > m_slopeCeNew[SpaceDim];

    LayoutData< std::shared_ptr<AggStencil <EBCellFAB, BaseIVFAB<Real> > > > m_stenLo[SpaceDim];
    LayoutData< std::shared_ptr<AggStencil <EBCellFAB, BaseIVFAB<Real> > > > m_stenHi[SpaceDim];

  private:
    //disallowed for all the usual reasons.  No code specified because Jeff says that is better.
    AggEBPWLFillPatch(const AggEBPWLFillPatch& ebcin);

    void operator=(const AggEBPWLFillPatch& fabin);

  };

}
#endif
