#ifndef _DIRICHLETCONDUCTIVITYEBBC_H___
#define _DIRICHLETCONDUCTIVITYEBBC_H___

#include "AMReX_ConductivityBaseEBBC.H"

namespace amrex
{

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
                          const IntVectSet& a_cfivs);

    //first order normal gradient stencil
    void
    getFirstOrderStencil(VoFStencil&     a_stencil,
                         Real&           a_weight,
                         const VolIndex& a_vof,
                         const EBISBox&  a_ebisBox);


    //second order normal gradient stencil
    bool
    getSecondOrderStencil(VoFStencil&          a_stencil,
                          Real&                a_weight,
                          vector<VoFStencil>&  a_pointStencils,
                          vector<Real>&        a_distanceAlongLine,
                          const VolIndex&      a_vof,
                          const EBISBox&       a_ebisBox,
                          const IntVectSet&    a_cfivs);



    int m_order;
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
    virtual ~DirichletConductivityEBBCFactory()
      {
      }

    void setOrder(int a_order)
      {
        m_order = a_order;
      }

    ///
    virtual ConductivityBaseEBBC* new_object_ptr()
      {
        DirichletConductivityEBBC* fresh = new DirichletConductivityEBBC();
        fresh->setOrder(m_order);
        return static_cast<ConductivityBaseEBBC*>(fresh);
      }

  private:
    int m_order;
  };
}
#endif
