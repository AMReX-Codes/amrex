#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _DIRICHLETCONDUCTIVITYEBBC_H_
#define _DIRICHLETCONDUCTIVITYEBBC_H_

#include "RefCountedPtr.H"

#include "DirichletPoissonEBBC.H"
#include "BaseBCValue.H"
#include "NamespaceHeader.H"

///
/**
 */
class DirichletConductivityEBBC: public ConductivityBaseEBBC
{
public:
  ///
  virtual void define(const LayoutData<IntVectSet>& a_cfivs,
                      const Real&                   a_factor);

  ///
  virtual void applyEBFlux(EBCellFAB&                    a_lphi,
                           const EBCellFAB&              a_phi,
                           VoFIterator&                  a_vofit,
                           const LayoutData<IntVectSet>& a_cfivs,
                           const DataIndex&              a_dit,
                           const RealVect&               a_probLo,
                           const RealVect&               a_dx,
                           const Real&                   a_factor,
                           const bool&                   a_useHomogeneous,
                           const Real&                   a_time);

  ///
  DirichletConductivityEBBC(const ProblemDomain& a_domain,
                            const EBISLayout&    a_layout,
                            const RealVect&      a_dx,
                            const IntVect*       a_ghostCellsPhi,
                            const IntVect*       a_ghostCellsRhs):ConductivityBaseEBBC()

  {
    m_bc.construct(a_domain, a_layout, a_dx, a_ghostCellsPhi, a_ghostCellsRhs);
  }

  virtual LayoutData<BaseIVFAB<VoFStencil> >* getFluxStencil(int ivar)
  {
    return &m_fluxStencil;
  }

  ///
  virtual ~DirichletConductivityEBBC();

  ///
  virtual void setValue(Real a_value);

  ///
  virtual int whichBC()
  {
    return 0;
  };

  virtual void setOrder(int a_order)
  {
    m_bc.setOrder(a_order);
  }

  ///
  virtual void setFunction(RefCountedPtr<BaseBCValue> a_flux);

  ///
  virtual void getEBFlux(Real&                         a_flux,
                         const VolIndex&               a_vof,
                         const LevelData<EBCellFAB>&   a_phi,
                         const LayoutData<IntVectSet>& a_cfivs,
                         const DataIndex&              a_dit,
                         const RealVect&               a_probLo,
                         const RealVect&               a_dx,
                         const bool&                   a_useHomogeneous,
                         const Real&                   a_time,
                         const pair<int,Real>*         a_cacheHint=0 );


private:
  DirichletPoissonEBBC m_bc;
  LayoutData<BaseIVFAB<VoFStencil> > m_fluxStencil;
};

///
class DirichletConductivityEBBCFactory: public BaseEBBCFactory
{
public:
  ///
  DirichletConductivityEBBCFactory();

  ///
  virtual ~DirichletConductivityEBBCFactory();

  virtual void setOrder(int a_order)
  {
    m_order = a_order;
  }

  ///
  virtual void setValue(Real a_value);

  ///
  virtual void setFunction(RefCountedPtr<BaseBCValue> a_flux);

  ///
  virtual DirichletConductivityEBBC* create(const ProblemDomain& a_domain,
                                            const EBISLayout&    a_layout,
                                            const RealVect&      a_dx,
                                            const IntVect*       a_ghostCellsPhi=0,
                                            const IntVect*       a_ghostCellsRhs=0);

private:
  bool m_onlyHomogeneous;
  bool m_isFunction;
  int m_order;
  Real m_value;
  RefCountedPtr<BaseBCValue> m_flux;
  RefCountedPtr<LevelData<BaseIVFAB<Real> > >m_data;
  bool m_dataBased;
};

#include "NamespaceFooter.H"
#endif
