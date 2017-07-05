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

#ifndef _EBCONDUCTIVITYOPFACTORY_H___
#define _EBCONDUCTIVITYOPFACTORY_H___

#include "REAL.H"
#include "Box.H"
#include "FArrayBox.H"
#include <vector>
#include "AMReX_AMRMultiGrid.H"
#include "AMReX_EBIndexSpace.H"
#include "AMReX_EBCellFAB.H"
#include "AMReX_EBCellFactory.H"
#include "AMReX_EBConductivityOp.H"
#include "AMReX_EBFabArrayOps.H"
#include "AMReX_BaseEBBC.H"
#include "AMReX_BaseDomainBC.H"
#include "AMReX_EBFluxRegister.H"
#include "AMReX_EBMGAverage.H"
#include "AMReX_EBMGInterp.H"

namespace amrex
{

//! \class EBConductivityOpFactory
//! Factory class to generate EBConductivityOps.  This follows the
//! AMRLevelOpFactory interface.
class EBConductivityOpFactory: public AMRLevelOpFactory<FabArray<EBCellFAB> >
{
public:

  //! Constructs a factory that builds EBConductivityOps with time-independent
  //! A and B coefficients.
  EBConductivityOpFactory(const vector<EBLevelGrid>&                                  a_eblgs,
                             const vector<shared_ptr<EBQuadCFInterp> >&             a_quadCFI,
                             const Real&                                                 a_alpha,
                             const Real                                         &        a_beta,
                             const vector<shared_ptr<FabArray<EBCellFAB> > >&        a_acoef,
                             const vector<shared_ptr<FabArray<EBFluxFAB> > >&        a_bcoef,
                             const vector<shared_ptr<FabArray<BaseIVFAB<Real> > > >& a_bcoefIrreg,
                             const Real&                                                 a_dxCoarse,
                             const vector<int>&                                          a_refRatio,
                             const shared_ptr<BaseDomainBCFactory>&                   a_domainBCFactory,
                             const shared_ptr<BaseEBBCFactory>    &                   a_ebBcFactory,
                             const IntVect&                                              a_ghostCellsPhi,
                             const IntVect&                                              a_ghostCellsRhs,
                             const int&                                                  a_relaxType,
                             int a_numLevels = -1);


  //! Destructor.
  virtual ~EBConductivityOpFactory();

  ///
  virtual void setData( vector< shared_ptr<FabArray<BaseIVFAB<Real> > > >& a_data)
  {
    m_data = a_data;
    m_dataBased = true;
  }

  ///
  virtual EBConductivityOp*
  MGnewOp(const ProblemDomain& a_FineindexSpace,
          int                  a_depth,
          bool                 a_homoOnly = true);

  EBConductivityOp* createOperator(const EBLevelGrid&             a_eblgMGLevel,
                                   const EBLevelGrid&             a_eblgCoarMG,
                                   const bool&                    a_hasMGObjects,
                                   const RealVect&                a_dxMGLevel,
                                   const RealVect&                a_dxCoar,
                                   const int&                     a_whichLevel);
  ///
  virtual void reclaim(MGLevelOp<FabArray<EBCellFAB> >* a_reclaim);

  ///
  virtual EBConductivityOp*
  AMRnewOp(const ProblemDomain& a_FineindexSpace);

  ///
  virtual void AMRreclaim(EBConductivityOp* a_reclaim);

  ///
  /** Refinement ratio between this level and coarser level.
      Returns 1 when there are no coarser AMRLevelOp objects */
  virtual int refToFiner(const ProblemDomain& a_domain) const;



protected:

  vector<shared_ptr<FabArray<BaseIVFAB<Real> > > > m_data;
  bool m_dataBased;
  vector< vector<EBLevelGrid> >                                   m_eblgsMG;

  //! Time-independent A multigrid coefficients
  vector< vector< shared_ptr<FabArray<EBCellFAB> > > >        m_acoefMG;

  //! Beginning-of-step (time-dependent) A multigrid coefficients
  vector< vector< shared_ptr<FabArray<EBCellFAB> > > >        m_acoefMG0;

  //! End-of-step (time-dependent) A multigrid coefficients
  vector< vector< shared_ptr<FabArray<EBCellFAB> > > >        m_acoefMG1;

  vector< vector< shared_ptr<FabArray<EBFluxFAB> > > >        m_bcoefMG;
  vector< vector< shared_ptr<FabArray<BaseIVFAB<Real> > > > > m_bcoefIrregMG;

  int                                                  m_relaxType;
  vector<EBLevelGrid>                                  m_eblgs;
  vector<shared_ptr<EBQuadCFInterp> >                  m_quadCFI;
  Real                                                 m_alpha;
  Real                                                 m_beta;

  //! Time-independent A coefficients.
  vector<shared_ptr<FabArray<EBCellFAB> > >        m_acoef;


  vector<shared_ptr<FabArray<EBFluxFAB> > >        m_bcoef;
  vector<shared_ptr<FabArray<BaseIVFAB<Real> > > > m_bcoefIrreg;
  Real                                                 m_dxCoarse;
  vector<int>                                          m_refRatio;
  shared_ptr<BaseDomainBCFactory>                   m_domainBCFactory;
  shared_ptr<BaseEBBCFactory>                       m_ebBCFactory;
  IntVect                                              m_ghostCellsPhi;
  IntVect                                              m_ghostCellsRhs;
  int                                                  m_numLevels;
  std::vector< bool  >                                 m_hasMGObjects;
  vector<Real>                                         m_dx;
private:
  ///weak construction bad
  EBConductivityOpFactory()
  {
    MayDay::Error("invalid operator");
  }

  //copy constructor and operator= disallowed for all the usual reasons
  EBConductivityOpFactory(const EBConductivityOpFactory& a_opin)
  {
    MayDay::Error("invalid operator");
  }

  //copy constructor and operator= disallowed for all the usual reasons
  void operator=(const EBConductivityOpFactory& a_opin)
  {
    MayDay::Error("invalid operator");
  }
};

}
#endif
