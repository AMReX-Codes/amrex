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

#ifndef _NWOEBCONDUCTIVITYOPFACTORY_H___
#define _NWOEBCONDUCTIVITYOPFACTORY_H___

#include "REAL.H"
#include "Box.H"
#include "FArrayBox.H"
#include <vector>
#include "AMReX_AMRMultiGrid.H"
#include "AMReX_EBIndexSpace.H"
#include "AMReX_EBCellFAB.H"
#include "AMReX_EBCellFactory.H"
#include "AMReX_NWOEBConductivityOp.H"
#include "AMReX_EBLevelDataOps.H"
#include "AMReX_BaseEBBC.H"
#include "AMReX_BaseDomainBC.H"
#include "AMReX_EBFluxRegister.H"
#include "AMReX_EBMGAverage.H"
#include "AMReX_EBMGInterp.H"

namespace amrex
{

//! \class NWOEBConductivityOpFactory
//! Factory class to generate NWOEBConductivityOps.  This follows the
//! AMRLevelOpFactory interface.
class NWOEBConductivityOpFactory: public AMRLevelOpFactory<LevelData<EBCellFAB> >
{
public:

  //! Constructs a factory that builds NWOEBConductivityOps with time-independent
  //! A and B coefficients.
  NWOEBConductivityOpFactory(const Vector<EBLevelGrid>&                                  a_eblgs,
                             const Vector<RefCountedPtr<NWOEBQuadCFInterp> >&             a_quadCFI,
                             const Real&                                                 a_alpha,
                             const Real                                         &        a_beta,
                             const Vector<RefCountedPtr<LevelData<EBCellFAB> > >&        a_acoef,
                             const Vector<RefCountedPtr<LevelData<EBFluxFAB> > >&        a_bcoef,
                             const Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >& a_bcoefIrreg,
                             const Real&                                                 a_dxCoarse,
                             const Vector<int>&                                          a_refRatio,
                             const RefCountedPtr<BaseDomainBCFactory>&                   a_domainBCFactory,
                             const RefCountedPtr<BaseEBBCFactory>    &                   a_ebBcFactory,
                             const IntVect&                                              a_ghostCellsPhi,
                             const IntVect&                                              a_ghostCellsRhs,
                             const int&                                                  a_relaxType,
                             int a_numLevels = -1);


  //! Destructor.
  virtual ~NWOEBConductivityOpFactory();

  ///
  virtual void setData( Vector< RefCountedPtr<LevelData<BaseIVFAB<Real> > > >& a_data)
  {
    m_data = a_data;
    m_dataBased = true;
  }

  ///
  virtual NWOEBConductivityOp*
  MGnewOp(const ProblemDomain& a_FineindexSpace,
          int                  a_depth,
          bool                 a_homoOnly = true);

  NWOEBConductivityOp* createOperator(const EBLevelGrid&             a_eblgMGLevel,
                                   const EBLevelGrid&             a_eblgCoarMG,
                                   const bool&                    a_hasMGObjects,
                                   const RealVect&                a_dxMGLevel,
                                   const RealVect&                a_dxCoar,
                                   const int&                     a_whichLevel);
  ///
  virtual void reclaim(MGLevelOp<LevelData<EBCellFAB> >* a_reclaim);

  ///
  virtual NWOEBConductivityOp*
  AMRnewOp(const ProblemDomain& a_FineindexSpace);

  ///
  virtual void AMRreclaim(NWOEBConductivityOp* a_reclaim);

  ///
  /** Refinement ratio between this level and coarser level.
      Returns 1 when there are no coarser AMRLevelOp objects */
  virtual int refToFiner(const ProblemDomain& a_domain) const;

  ///
  /**
     testRef is the size of the coarsest domain allowed in multigrid. If testRef=2,
     then the coarsest domain in multigrid will be 2x2(x2)
   **/
  static void setTestRef(int a_testRef)
  {
    s_testRef = a_testRef;
  }

  static void setMaxBoxSize(int a_maxBoxSize)
  {
    s_maxBoxSize = a_maxBoxSize;
  }

  static int s_testRef;
  static int s_maxBoxSize;

  //! Reset the coefficients for the conductivity operator (time-independent).
  void resetCoefficients(const Real&                                                 a_alpha,
                         const Real                                         &        a_beta,
                         const Vector<RefCountedPtr<LevelData<EBCellFAB> > >&        a_acoef,
                         const Vector<RefCountedPtr<LevelData<EBFluxFAB> > >&        a_bcoef,
                         const Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >& a_bcoefIrreg)
  {
    m_alpha       = a_alpha      ;
    m_beta        = a_beta       ;
    m_acoef       = a_acoef      ;
    m_bcoef       = a_bcoef      ;
    m_bcoefIrreg  = a_bcoefIrreg ;
  }


protected:

  Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > > m_data;
  bool m_dataBased;
  Vector< Vector<EBLevelGrid> >                                   m_eblgsMG;

  //! Time-independent A multigrid coefficients
  Vector< Vector< RefCountedPtr<LevelData<EBCellFAB> > > >        m_acoefMG;

  //! Beginning-of-step (time-dependent) A multigrid coefficients
  Vector< Vector< RefCountedPtr<LevelData<EBCellFAB> > > >        m_acoefMG0;

  //! End-of-step (time-dependent) A multigrid coefficients
  Vector< Vector< RefCountedPtr<LevelData<EBCellFAB> > > >        m_acoefMG1;

  Vector< Vector< RefCountedPtr<LevelData<EBFluxFAB> > > >        m_bcoefMG;
  Vector< Vector< RefCountedPtr<LevelData<BaseIVFAB<Real> > > > > m_bcoefIrregMG;

  int                                                  m_relaxType;
  Vector<EBLevelGrid>                                  m_eblgs;
  Vector<RefCountedPtr<NWOEBQuadCFInterp> >            m_quadCFI;
  Real                                                 m_alpha;
  Real                                                 m_beta;

  //! Time-independent A coefficients.
  Vector<RefCountedPtr<LevelData<EBCellFAB> > >        m_acoef;


  Vector<RefCountedPtr<LevelData<EBFluxFAB> > >        m_bcoef;
  Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > > m_bcoefIrreg;
  Real                                                 m_dxCoarse;
  Vector<int>                                          m_refRatio;
  RefCountedPtr<BaseDomainBCFactory>                   m_domainBCFactory;
  RefCountedPtr<BaseEBBCFactory>                       m_ebBCFactory;
  IntVect                                              m_ghostCellsPhi;
  IntVect                                              m_ghostCellsRhs;
  int                                                  m_numLevels;
  std::vector< bool  >                                 m_hasMGObjects;
  Vector<Real>                                         m_dx;
private:
  ///weak construction bad
  NWOEBConductivityOpFactory()
  {
    MayDay::Error("invalid operator");
  }

  //copy constructor and operator= disallowed for all the usual reasons
  NWOEBConductivityOpFactory(const NWOEBConductivityOpFactory& a_opin)
  {
    MayDay::Error("invalid operator");
  }

  //copy constructor and operator= disallowed for all the usual reasons
  void operator=(const NWOEBConductivityOpFactory& a_opin)
  {
    MayDay::Error("invalid operator");
  }
};

}
#endif
