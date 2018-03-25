#ifndef _EBFASTFR_H___
#define _EBFASTFR_H___

#include "AMReX_REAL.H"
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
#include "AMReX_FluxRegister.H"
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
    ///
    /**
       Default constructor.  Leaves object undefined.
    */
    EBFastFR()
      {
        m_isDefined = false;
      }

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
             bool a_forceNoEBCF = true);

    ///
    virtual ~EBFastFR()
      {
      }

    ///
    /**
       Full define function.  Creates a levels worth of flux registers.
    */
    virtual void
    define(const EBLevelGrid&       a_eblgFine,
           const EBLevelGrid&       a_eblgCoar,
           const int&               a_nref,
           const int&               a_nvar,
           bool  a_forceNoEBCF = true);

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
       The face direction is known by EBFaceFAB
    */
    virtual void
    incrementCoarse(const EBFaceFAB&      a_coarseFlux,
                    const Real&           a_scale,
                    const MFIter   &      a_coarseMFIter,
                    int isrc, int idst, int inco);



    ///
    virtual void
    incrementFine(const EBFaceFAB&      a_fineFlux,
                  const Real&           a_scale,
                  const MFIter   &      a_fineMFIter   ,
                  int isrc, int idst, int inco);


    ///
    virtual void reflux(FabArray<EBCellFAB> & a_uCoarse,
                        const Real&           a_scale,
                        int isrc, int idst, int inco);

    ///
    bool isDefined() const
      {
        return m_isDefined;
      }



    bool hasEBCF() const
      {
        return m_hasEBCF;
      }

    
  protected:

    //for the regular bits
    FluxRegister m_nonEBFluxReg;

    //because amrex::FluxRegister is weird, the fluxes at the coarse level have to be stored
    FabArray<EBFluxFAB> m_coarRegister;

    bool m_isDefined;

    EBLevelGrid       m_eblgFine;
    EBLevelGrid       m_eblgCoar;
    EBLevelGrid       m_eblgCoFi;

    int               m_refRat;
    int               m_nComp;

    bool m_hasEBCF;


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
