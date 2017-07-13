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
  virtual void setCoef(EBLevelGrid                                &  a_eblg,
                       Real                                       &  a_beta,
                       Real                                       &  a_dx,
                       const RealVect                             &  a_probLo,
                       shared_ptr<FabArray<EBFluxFAB> >           &  a_bcoe)
    {
      ConductivityBaseEBBC::setCoef(a_eblg, a_beta, a_dx, a_probLo, a_bcoe);
      defineStencils();
    }

  ///
  virtual void applyEBFlux(EBCellFAB&                    a_lphi,
                           const EBCellFAB&              a_phi,
                           const vector<VolIndex>&       a_vofsToChange,
                           const MFIter  &               a_mfi,
                           const Real&                   a_factor,
                           const bool&                   a_useHomogeneous);

  ///
  DirichletConductivityEBBC():ConductivityBaseEBBC()

  {
  }

  virtual LayoutData<BaseIVFAB<VoFStencil> >* getFluxStencil(int ivar)
  {
    BL_ASSERT(m_coefSet);
    return &m_fluxStencil;
  }

  ///
  virtual ~DirichletConductivityEBBC()
    {
    }

  virtual void setOrder(int a_order)
  {
    m_order = a_order;
  }

private:
  void defineStencils();
  LayoutData<BaseIVFAB<VoFStencil> > m_fluxStencil;
};

///
class DirichletConductivityEBBCFactory: public ConductivityBaseEBBCFactory
{
public:
  ///
  DirichletConductivityEBBCFactory()
    {
      m_order = 2;
    }

  ///
  virtual ~DirichletConductivityEBBCFactory();

  virtual void setOrder(int a_order)
  {
    m_order = a_order;
  }

  virtual void setValue(Real a_value);

  ///
  virtual void setFunction(RefCountedPtr<BaseBCValue> a_flux);

  ///
  virtual ConductivityBaseEBBC* new_object_ptr()
    {
     DirichletConductivityEBBC* fresh = new NeumannConductivityEBBC(a_domain,a_layout,a_dx);
     fresh->setOrder(m_order);
     return static_cast<ConductivityBaseEBBC*>(fresh);
    }

private:
  int m_order;
};

#include "NamespaceFooter.H"
#endif
