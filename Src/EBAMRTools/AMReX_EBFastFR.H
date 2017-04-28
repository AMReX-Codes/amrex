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


#ifndef _EBFASTFR_H___
#define _EBFASTFR_H___

#include "AMReX_REAL.H"
#include "AMReX_Vector.H"
#include "AMReX_EBCellFAB.H"
#include "AMReX_EBFaceFAB.H"
#include "AMReX_EBCellFactory.H"
#include "AMReX_EBLevelDataOps.H"
#include "AMReX_EBISLayout.H"
#include "AMReX_EBISBox.H"
#include "AMReX_IntVectSet.H"
#include "AMReX_LoHiSide.H"
#include "AMReX_EBLevelGrid.H"
#include "AMReX_LayoutData.H"
class EBIndexSpace;

namespace amrex
{
///EBFastFR-A class to encapsulate a levels worth of flux registers.
/**
   A EBFastFR  handles all the data choreography
   necessary to create a levels worth of flux registers.
   This only does the elliptic flux register stuff for now (no interacting with
   redistribution, no RZ, no regular and regular separate).   
   Uses the non-EB flux register class everwhere except right near the EB.
*/
  class EBFastFR
  {
  public:
    static bool s_verbose;
    ///
    /**
       Default constructor.  Leaves object undefined.
    */
    EBFastFR();

    ///
    /**
       Full constructor.  Calls the define function which creates
       a levels worth of flux registers.
       noebcf forces no coarse-fine intersection with EB
    */
    EBFastFR(const EBLevelGrid&       a_eblgFine,
             const EBLevelGrid&       a_eblgCoar,
             const int&               a_nref,
             const int&               a_nvar,
             bool a_forceNoEBCF = false);

    ///
    virtual ~EBFastFR();

    ///
    /**
       Full define function.  Creates a levels worth of flux registers.
    */
    virtual void
    define(const EBLevelGrid&       a_eblgFine,
           const EBLevelGrid&       a_eblgCoar,
           const int&               a_nref,
           const int&               a_nvar,
           bool  a_forceNoEBCF = false);

    ///
    /**
       Initialize values of registers  to zero.
    */
    virtual void
    setToZero();


    ///
    /**
       The irregular part of this is just
       register = input flux  (instead of increment)
       To avoid double-counting.
    */
    virtual void
    incrementCoarse(const EBFaceFAB&      a_coarseFlux,
                    const Real&           a_scale,
                    const DataIndex&      a_coarseDataIndex,
                    const Interval&       a_variables,
                    const int&            a_dir,
                    const Side::LoHiSide& a_sd);


    ///
    virtual void
    incrementFine(const EBFaceFAB&      a_fineFlux,
                  const Real&           a_scale,
                  const DataIndex&      a_fineDataIndex,
                  const Interval&       a_variables,
                  const int&            a_dir,
                  const Side::LoHiSide& a_sd);

    ///
    virtual void reflux(LevelData<EBCellFAB>& a_uCoarse,
                        const Interval&       a_variables,
                        const Real&           a_scale,
                        bool a_multByKappaOneMinusKappa = false);

    ///
    /**
       Increments mass array  with left-over mass
       from reflux divergence.   this is to test this pig.
       Ordinarily this mass would go into redistribution.   this
       is used in the test suite to  check whether a constant
       flux refluxed and then unrefluxed (this function) ends up
       with a constant solution.
       Correction at each cell = (1-kappa)refluxCorrection.
    */
    void incrementDensityArray(LevelData<EBCellFAB>& a_coarMass,
                               const Interval&       a_variables,
                               const Real&           a_scale);

    ///
    bool isDefined() const;

    static int index(int a_dir, Side::LoHiSide a_side);


    bool hasEBCF() const
      {
        return m_hasEBCF;
      }

    ///undefines object
    void clear();

    ///
    /**
       Gets the set of vofs to fill a flux for on a particular side
       and direction that correspond to an EB/CF interface.   You still
       need to fill all the other fluxes at the interface but you do not
       have to do all the fancy interpolation to face centroids
       that is needed at these key points since all the other C/F points
       are treated as regular anyway.
       Do not use this unless you really know what you are doing.   Really.
       Yes this breaks encapsulation but we do some pretty ugly things
       for performance reasons, do we not?
    */
    Vector<VoFIterator>& getVoFItCoar(const DataIndex&       a_dit,
                                      const int&             a_idir,
                                      const Side::LoHiSide&  a_sd)
      {
        int iloc = index(a_idir, a_sd);
        return (m_vofiCoar[iloc])[a_dit];
      }

    ///
    /**
       Gets the set of vofs to fill a flux for on a particular side
       and direction that correspond to an EB/CF interface.   You still
       need to fill all the other fluxes at the interface but you do not
       have to do all the fancy interpolation to face centroids
       that is needed at these key points since all the other C/F points
       are treated as regular anyway.
       Do not use this unless you really know what you are doing.   Really.
       Yes this breaks encapsulation but we do some pretty ugly things
       for performance reasons, do we not?
    */
    VoFIterator& getVoFItCoFi(const DataIndex&       a_dit,
                              const int&             a_idir,
                              const Side::LoHiSide&  a_sd)
      {
        int iloc = index(a_idir, a_sd);
        return (m_vofiCoFi[iloc])[a_dit];
      }

    ///
    const EBLevelGrid& getEBLGCoFi() const
      {
        return m_eblgCoFi;
      }

  protected:

    //for the regular bits
    FluxRegister m_levelFluxReg;
    bool m_isDefined;

    void setDefaultValues();
    EBLevelGrid       m_eblgFine;
    EBLevelGrid       m_eblgCoar;
    EBLevelGrid       m_eblgCoFi;

    int               m_refRat;
    int               m_nComp;

    //nrefdmo is the number of fine faces per coar face
    //this is intrinsically dimension-dependent
    Real m_nrefdmo;
    bool m_hasEBCF;


    bool computeHasEBCF();
  private:
    //for all the usual reasons,
    //there is no copy constructor for this class
    //and there is no operator= for this class
    //(there is no rule six).
    void operator= (const EBFastFR& out);
    EBFastFR(const EBFastFR& out);
  };

}
#endif
