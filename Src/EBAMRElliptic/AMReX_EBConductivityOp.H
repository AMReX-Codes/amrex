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

#ifndef _EBCONDUCTIVITYOP_H___
#define _EBCONDUCTIVITYOP_H___

#include <vector>

#include "AMReX_Box.H"
#include "AMReX_AMREBMultiGrid.H"
#include "AMReX_EBIndexSpace.H"
#include "AMReX_EBCellFAB.H"
#include "AMReX_EBCellFactory.H"
#include "AMReX_VCAggStencil.H"
#include "AMReX_EBLevelDataOps.H"
#include "AMReX_EBFastFR.H"
#include "AMReX_EBCoarseAverage.H"
#include "AMReX_EBFineInterp.H"
#include "AMReX_EBCFInterp.H"
#include "AMReX_EBLevelGrid.H"
#include "AMReX_FabArray.H"
#include "AMReX_ConductivityBaseDomainBC.H"
#include "AMReX_ConductivityBaseEBBC.H"


namespace amrex
{
  using std::shared_ptr;
  //! \class EBConductivityOp
  //! This class implements an operator that solves the equation
  //! (alpha a + beta div (b grad) )phi = rhs
  //! using the AMRLevelOp interface.
  class EBConductivityOp: public AMRLevelOp< FabArray<EBCellFAB> >
  {
  public:
    
    //! Constructs a conductivity operator using the given data. This
    //! constructor is for time-independent a and b coefficients.
    //! If you are approaching this operator from this interface, consider backing away and using
    //! EBConductivityOpFactory to generate these objects. Really.
    //! Ghost cell arguments are there for caching reasons. Once you set them,
    //! an error is thrown if you send in data that does not match.
    //! \param a_eblgFine grid at finer  level
    //! \param a_eblg grid at this  level
    //! \param a_eblgCoar grid at coarser level
    //! \param a_eblgCoarMG grid at intermediate EBMultiGrid level
    //! \param a_domainBC domain boundary conditions at this level
    //! \param a_ebBC eb boundary conditions at this level
    //! \param a_dx grid spacing at this level
    //! \param a_origin offset to lowest corner of the domain
    //! \param a_refToFine refinement ratio to finer level
    //! \param a_refToCoar refinement ratio to coarser level
    //! \param a_hasFiner true if there is a finer AMR level, false otherwise.
    //! \param a_hasCoarser true if there is a coarser AMR level.
    //! \param a_hasCoarserMG true if there is a coarser EBMultiGrid level.
    //! \param a_preCondIters number of iterations to do for pre-conditioning
    //! \param a_acoef coefficent of identity
    //! \param a_bcoef coefficient of gradient.
    //! \param a_ghostCellsPhi Number of ghost cells in phi, correction
    //! \param a_ghostCellsRhs Number of ghost cells in RHS, residual, lphi
    EBConductivityOp(const EBLevelGrid &                               a_eblgFine,
                     const EBLevelGrid &                               a_eblg,
                     const EBLevelGrid &                               a_eblgCoar,
                     const EBLevelGrid &                               a_eblgCoarMG,
                     const shared_ptr<ConductivityBaseDomainBC>&       a_domainBC,
                     const shared_ptr<ConductivityBaseEBBC>&           a_ebBC,
                     const Real    &                                   a_dx,
                     const Real    &                                   a_dxCoar,
                     const int&                                        a_refToFine,
                     const int&                                        a_refToCoar,
                     const bool&                                       a_hasFine,
                     const bool&                                       a_hasCoar,
                     const bool&                                       a_hasCoarMG,                  
                     const Real&                                       a_alpha,
                     const Real&                                       a_beta,
                     const shared_ptr<FabArray<EBCellFAB> >&           a_acoef,
                     const shared_ptr<FabArray<EBFluxFAB> >&           a_bcoef,
                     const int&                                        a_ghost);
    
    
    ~EBConductivityOp()
      {
      }
    
    
    ///
    /**
       This sets the data storage for the a coefficient to a different
       object and recalculates the stuff it depends on. Use this only if you know what you're doing.
    */
    virtual void resetACoefficient(shared_ptr<FabArray<EBCellFAB> >& a_acoef)
      {
    
        m_acoef = a_acoef;
        calculateAlphaWeight();
        calculateRelaxationCoefficient();
      }
    
    //only weights by kappa.  time and tga have their demands.
    virtual void kappaScale(FabArray<EBCellFAB> & a_rhs);
    
    ///
    /** a_residual = a_rhs - L(a_phiFine, a_phi)   no coaser AMR level*/
    void AMRResidualNC(FabArray<EBCellFAB>&       a_residual,
                       const FabArray<EBCellFAB>& a_phiFine,
                       const FabArray<EBCellFAB>& a_phi,
                       const FabArray<EBCellFAB>& a_rhs,
                       bool a_homogeneousBC,
                       AMRLevelOp<FabArray<EBCellFAB> >* a_finerOp);
    
    
    ///
    /** apply AMR operator   no coaser AMR level*/
    void AMROperatorNC(FabArray<EBCellFAB>&       a_LofPhi,
                       const FabArray<EBCellFAB>& a_phiFine,
                       const FabArray<EBCellFAB>& a_phi,
                       bool a_homogeneousBC,
                       AMRLevelOp<FabArray<EBCellFAB> >* a_finerOp);
    
    //------------------------------------
    // Overridden methods for base classes
    //------------------------------------
    
    void setAlphaAndBeta(const Real& a_alpha,
                         const Real& a_beta);
    
    void diagonalScale(FabArray<EBCellFAB> & a_rhs,
                       bool a_kappaWeighted);
    
    void divideByIdentityCoef(FabArray<EBCellFAB> & a_rhs);
    
    void fillPhiGhost(const EBCellFAB& a_phi, const MFIter& a_datInd, bool a_homog) const;
    void fillPhiGhost(const FabArray<EBCellFAB>& a_phi, bool a_homog) const;
    
    
    ///
    /**
     */
    virtual void residual(FabArray<EBCellFAB>&       a_residual,
                          const FabArray<EBCellFAB>& a_phi,
                          const FabArray<EBCellFAB>& a_rhs,
                          bool                        a_homogeneousPhysBC=false);
    
    ///
    /**
     */
    virtual void preCond(FabArray<EBCellFAB>&       a_opPhi,
                         const FabArray<EBCellFAB>& a_phi);
    
    ///
    /**
       This function assumes that coarse-fine boundary condtions have
       been dealt with.
    */
    virtual void applyOp(FabArray<EBCellFAB>&             a_opPhi,
                         const FabArray<EBCellFAB>&       a_phi,
                         const FabArray<EBCellFAB>* const a_phiCoarse,
                         const bool&                       a_homogeneousPhysBC,
                         const bool&                       a_homogeneousCFBC);
    
    /// virtual function called by LevelTGA
    virtual void applyOpNoBoundary(FabArray<EBCellFAB>&        a_opPhi,
                                   const FabArray<EBCellFAB>&  a_phi)
      {
        m_turnOffBCs = true;
        applyOp(a_opPhi, a_phi, true);
        m_turnOffBCs = false;
      }
    
    void
    applyOpIrregular(EBCellFAB&             a_lhs,
                     const EBCellFAB&       a_phi,
                     const bool&            a_homogeneous,
                     const MFIter&       a_datInd);
    
    void
    applyOpRegular(EBCellFAB&             a_lhs,
                   const EBCellFAB&       a_phi,
                   const bool&            a_homogeneous,
                   const MFIter&       a_datInd);
    
    
    ///
    /**
       this is the linearop function.  CFBC is set to homogeneous.  phic is null
    */
    virtual void applyOp(FabArray<EBCellFAB>       &        a_opPhi,
                         const FabArray<EBCellFAB> &        a_phi,
                         bool    a_homogeneousPhysBC = false);
    
    ///
    /**
     */
    virtual void create(FabArray<EBCellFAB>&       a_lhs,
                        const FabArray<EBCellFAB>& a_rhs);
    
    ///
    virtual void createCoarsened(FabArray<EBCellFAB>&       a_lhs,
                                 const FabArray<EBCellFAB>& a_rhs,
                                 const int&                  a_refRat);
    
    Real
    AMRNorm(const FabArray<EBCellFAB>& a_coarResid,
            const FabArray<EBCellFAB>& a_fineResid,
            const int& a_refRat,
            const int& a_ord);
    
    ///
    /**
     */
    virtual void assign(FabArray<EBCellFAB>&       a_lhs,
                        const FabArray<EBCellFAB>& a_rhs);
    
    ///
    /**
     */
    virtual Real dotProduct(const FabArray<EBCellFAB>& a_1,
                            const FabArray<EBCellFAB>& a_2);
    
    ///
    /**
     */
    virtual void incr(FabArray<EBCellFAB>&       a_lhs,
                      const FabArray<EBCellFAB>& a_x,
                      Real                        a_scale);
    
    ///
    /**
     */
    virtual void axby(FabArray<EBCellFAB>&       a_lhs,
                      const FabArray<EBCellFAB>& a_x,
                      const FabArray<EBCellFAB>& a_y,
                      Real                        a_a,
                      Real                        a_b);
    
    ///
    /**
     */
    virtual void scale(FabArray<EBCellFAB>& a_lhs,
                       const Real&           a_scale);
    
    ///
    /**
     */
    virtual Real norm(const FabArray<EBCellFAB>& a_rhs,
                      int                         a_ord);
    
    virtual Real localMaxNorm(const FabArray<EBCellFAB>& a_rhs);
    ///
    /**
     */
    virtual void setToZero(FabArray<EBCellFAB>& a_lhs);
    
    ///
    /**
     */
    virtual void setVal(FabArray<EBCellFAB>& a_lhs, const Real& a_value);
    
    ///
    /**
     */
    virtual void createCoarser(FabArray<EBCellFAB>&       a_coarse,
                               const FabArray<EBCellFAB>& a_fine,
                               bool                        a_ghosted);
    
    ///
    /**
     */
    virtual void relax(FabArray<EBCellFAB>&       a_e,
                       const FabArray<EBCellFAB>& a_residual,
                       int                         a_iterations);
    
    ///
    /**
       Calculate restricted residual:
       a_resCoarse[2h] = I[h->2h] (a_rhsFine[h] - L[h](a_phiFine[h]))
    */
    virtual void restrictResidual(FabArray<EBCellFAB>&       a_resCoarse,
                                  FabArray<EBCellFAB>&       a_phiFine,
                                  const FabArray<EBCellFAB>& a_rhsFine);
    
    ///
    /**
       Correct the fine solution based on coarse correction:
       a_phiThisLevel += I[2h->h] (a_correctCoarse)
    */
    virtual void prolongIncrement(FabArray<EBCellFAB>&       a_phiThisLevel,
                                  const FabArray<EBCellFAB>& a_correctCoarse);
    
    ///
    /** Refinement ratio between this level and coarser level.
        Returns 1 when there are no coarser AMRLevelOp objects */
    virtual int refToCoarser();
    
    ///
    /** Refinement ratio between this level and coarser level.
        Returns 1 when there are no coarser AMRLevelOp objects */
    virtual int refToFiner();
    
    ///
    /** a_residual = a_rhs - L(a_phi, a_phiFine, a_phiCoarse) */
    virtual void AMRResidual(FabArray<EBCellFAB>& a_residual,
                             const FabArray<EBCellFAB>& a_phiFine,
                             const FabArray<EBCellFAB>& a_phi,
                             const FabArray<EBCellFAB>& a_phiCoarse,
                             const FabArray<EBCellFAB>& a_rhs,
                             bool a_homogeneousBC,
                             AMRLevelOp<FabArray<EBCellFAB> >* a_finerOp);
    
    ///
    /** a_residual = a_rhs - L(a_phi, a_phiCoarse)  */
    virtual void AMRResidualNF(FabArray<EBCellFAB>& a_residual,
                               const FabArray<EBCellFAB>& a_phi,
                               const FabArray<EBCellFAB>& a_phiCoarse,
                               const FabArray<EBCellFAB>& a_rhs,
                               bool a_homogeneousBC);
    
    
    ///
    /** a_residual = a_rhs - L(a_phi, a_phiFine, a_phiCoarse) */
    virtual void AMROperator(FabArray<EBCellFAB>& a_LofPhi,
                             const FabArray<EBCellFAB>& a_phiFine,
                             const FabArray<EBCellFAB>& a_phi,
                             const FabArray<EBCellFAB>& a_phiCoarse,
                             bool a_homogeneousBC,
                             AMRLevelOp<FabArray<EBCellFAB> >* a_finerOp);
    
    ///
    /** a_residual = a_rhs - L(a_phi, a_phiCoarse)  */
    virtual void AMROperatorNF(FabArray<EBCellFAB>& a_LofPhi,
                               const FabArray<EBCellFAB>& a_phi,
                               const FabArray<EBCellFAB>& a_phiCoarse,
                               bool a_homogeneousBC);
    
    ///
    /** a_resCoarse = I[h-2h] (a_residual - L(a_correction, a_coarseCorrection)) */
    virtual void AMRRestrict(FabArray<EBCellFAB>& a_resCoarse,
                             const FabArray<EBCellFAB>& a_residual,
                             const FabArray<EBCellFAB>& a_correction,
                             const FabArray<EBCellFAB>& a_coarseCorrection);
    
    ///
    /** a_correction += I[2h->h](a_coarseCorrection) */
    virtual void AMRProlong(FabArray<EBCellFAB>&       a_correction,
                            const FabArray<EBCellFAB>& a_coarseCorrection);
    
    ///
    /** a_residual = a_residual - L(a_correction, a_coarseCorrection) */
    virtual void AMRUpdateResidual(FabArray<EBCellFAB>&       a_residual,
                                   const FabArray<EBCellFAB>& a_correction,
                                   const FabArray<EBCellFAB>& a_coarseCorrection);
    
    void
    gsrbColor(FabArray<EBCellFAB>&       a_phi,
              const FabArray<EBCellFAB>& a_rhs,
              const IntVect&              a_color);
    
    void reflux(FabArray<EBCellFAB>& a_residual,
                const FabArray<EBCellFAB>& a_phiFine,
                const FabArray<EBCellFAB>& a_phi,
                AMRLevelOp<FabArray<EBCellFAB> >* a_finerOp);
    
    void gsrbColor(FabArray<EBCellFAB>&       a_phi,
                   const FabArray<EBCellFAB>& a_lph,
                   const FabArray<EBCellFAB>& a_rhs,
                   const IntVect&              a_color);
    
    void getDivFStencil(VoFStencil&      a_vofStencil,
                        const VolIndex&  a_vof,
                        const MFIter  & a_dit);
    
    void getFluxStencil(VoFStencil&      a_fluxStencil,
                        const FaceIndex& a_face,
                        const MFIter   & a_dit);
    
    void getFaceCenteredFluxStencil(VoFStencil&      a_fluxStencil,
                                    const FaceIndex& a_face,
                                    const MFIter   & a_dit);
    
    void homogeneousCFInterp(FabArray<EBCellFAB>&   a_phif);
    
    //this is needed to do wacky things like reset the a coefficient.
    //not for the faint of heart
    //but none but the brave deserve the fair
    const EBLevelGrid& getEBLG() const
      {
        return m_eblg;
      }
    ///do not call this one unless you really know what you are doing
    void defineStencils();
  protected:
    
    //relaxation coeffcient safety factor
    Real getSafety();
    
    void incrOpRegularAllDirs(Box * a_loBox,
                              Box * a_hiBox,
                              int * a_hasLo,
                              int * a_hasHi,
                              Box & a_curDblBox,
                              Box & a_curPhiBox,
                              int a_nComps,
                              BaseFab<Real> & a_curOpPhiFAB,
                              const BaseFab<Real> & a_curPhiFAB,
                              bool a_homogeneousPhysBC,
                              const MFIter& a_dit);
    
    void applyDomainFlux(Box * a_loBox,
                         Box * a_hiBox,
                         int * a_hasLo,
                         int * a_hasHi,
                         Box & a_dblBox,
                         int a_nComps,
                         BaseFab<Real> & a_phiFAB,
                         bool a_homogeneousPhysBC,
                         const MFIter& a_dit);
    
    void GSColorAllIrregular(EBCellFAB&                        a_phi,
                             const EBCellFAB&                  a_rhs,
                             const int&                        a_icolor,
                             const MFIter&                  a_dit);
    
    bool   m_turnOffBCs;
    
    virtual void calculateAlphaWeight();
    virtual void calculateRelaxationCoefficient();
    
    
    void getFluxRegOnly(EBFaceFAB&                    a_fluxCentroid,
                        const EBCellFAB&              a_phi,
                        const Box&                    a_ghostedBox,
                        const MFIter&              a_datInd,
                        const int&                    a_idir);
    
    const int                   m_ghost;
    
    
    EBLevelGrid                     m_eblg;
    EBLevelGrid                     m_eblgFine;
    EBLevelGrid                     m_eblgCoar;
    EBLevelGrid                     m_eblgCoarMG;
    
    shared_ptr<ConductivityBaseDomainBC>     m_domainBC;
    shared_ptr<ConductivityBaseEBBC>         m_ebBC;
    
    Real                            m_dxFine;
    Real                            m_dx;
    Real                            m_dxCoar;
    
    shared_ptr<FabArray<EBCellFAB> >          m_acoef;
    shared_ptr<FabArray<EBFluxFAB> >          m_bcoef;
    
    Real                            m_alpha;
    Real                            m_beta;
    //weights that get multiplied by alpha
    FabArray<IrregFAB>       m_alphaDiagWeight;
    //weights that get multiplied by beta
    FabArray<IrregFAB>       m_betaDiagWeight;
    int                             m_refToFine;
    int                             m_refToCoar;

    bool                            m_hasFine;
    bool                            m_hasCoar;

    bool                            m_hasCoarMG;
    //restriction object
    EBCoarseAverage                    m_ebAverage;
    //prolongation object
    EBFineInterp                       m_ebInterp;
    //coarse-fine interpolation (with COARSER level)
    EBCFInterp                         m_cfInterp;
    //stencils for operator evaluation
    LayoutData<shared_ptr<VCAggStencil> >  m_opEBStencil;
    //stencils for operator evaluation on gauss-seidel colors
    
    //! EBMultiGrid relaxation coefficient
    FabArray<EBCellFAB>       m_relCoef;
    
    //cache the vofiterators
    //for irregular cell iteration (includes buffer around multivalued cells)
    LayoutData<VoFIterator >                     m_vofIterIrreg;
    LayoutData<VoFIterator >                     m_vofIterMulti;
    //for domain boundary conditions at ir regular cells
    LayoutData<VoFIterator >                     m_vofIterDomLo[BL_SPACEDIM];
    LayoutData<VoFIterator >                     m_vofIterDomHi[BL_SPACEDIM];
    
    
    //flux register with finer level
    EBFastFR       m_fastFR;
    
    EBCoarseAverage             m_ebAverageMG;
    EBFineInterp                m_ebInterpMG;

    
    vector<IntVect> m_colors;
    
  private:
    
    void incrementFRCoar(EBFastFR&                   a_fluxReg,
                         const FabArray<EBCellFAB>& a_phiFine,
                         const FabArray<EBCellFAB>& a_phi);
    
    void incrementFRFine(EBFastFR&                   a_fluxReg,
                         const FabArray<EBCellFAB>& a_phiFine,
                         const FabArray<EBCellFAB>& a_phi,
                         AMRLevelOp<FabArray<EBCellFAB> >* a_finerOp);
    
    void getOpVoFStencil(VoFStencil&     a_stencil,
                         const EBISBox&  a_ebisbox,
                         const VolIndex& a_vof);
    
    void getOpVoFStencil(VoFStencil&             a_stencil,
                         const int&              a_dir,
                         const vector<VolIndex>& a_allMonotoneVoFs,
                         const EBISBox&          a_ebisbox,
                         const VolIndex&         a_vof,
                         const bool&             a_lowOrder);
    
    
    void getOpFaceStencil(VoFStencil&             a_stencil,
                          const vector<VolIndex>& a_allMonotoneVofs,
                          const EBISBox&          a_ebisbox,
                          const VolIndex&         a_vof,
                          int                     a_dir,
                          const Side::LoHiSide&   a_side,
                          const FaceIndex&        a_face,
                          const bool&             a_lowOrder);
    
  private:
    
    //! Default constructor. Creates an undefined conductivity operator.
    EBConductivityOp();
    
    //copy constructor and operator= disallowed for all the usual reasons
    EBConductivityOp(const EBConductivityOp& a_opin);

    
    void operator=(const EBConductivityOp& a_opin);

  };
    
}
    
#endif
