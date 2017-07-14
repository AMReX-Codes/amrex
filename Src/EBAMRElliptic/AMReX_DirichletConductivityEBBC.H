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

  //second order normal gradient stencil (unless it drops order)
  void 
  getSecondOrderStencil(VoFStencil&       a_stencil,
                        Real&             a_weight,
                        const VolIndex&   a_vof,
                        const EBISBox&    a_ebisBox,
                        const RealVect&   a_dx,
                        const IntVectSet& a_cfivs);

  //first order normal gradient stencil
  void
  getFirstOrderStencil(VoFStencil&     a_stencil,
                       Real&           a_weight,
                       const VolIndex& a_vof,
                       const EBISBox&  a_ebisBox,
                       const RealVect& a_dx);

  //second order normal gradient stencil
  bool
  getSecondOrderStencil(VoFStencil&          a_stencil,
                        Real&                a_weight,
                        vector<VoFStencil>&  a_pointStencils,
                        vector<Real>&        a_distanceAlongLine,
                        const VolIndex&      a_vof,
                        const EBISBox&       a_ebisBox,
                        const RealVect&      a_dx,
                        const IntVectSet&    a_cfivs);



  LayoutData<BaseIVFAB<VoFStencil> > m_fluxStencil;
  LayoutData<BaseIVFAB<Real> > m_fluxWeight;
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
