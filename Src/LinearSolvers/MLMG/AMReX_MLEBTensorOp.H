#ifndef AMREX_ML_EB_TENSOROP_H_
#define AMREX_ML_EB_TENSOROP_H_
#include <AMReX_Config.H>

#include <AMReX_MLEBABecLap.H>

namespace amrex {

// Tensor solver for high Reynolds flows with small gradient in viscosity.
// The system it solves is
//
//   alpha a v - beta div dot tau = rhs
//
// where tau = eta [grad v + (grad v)^T] +  (kappa-(2/3)eta) (div v) I.
// Here eta and kappa are shear and bulk viscosity, and I is identity tensor.
//
// The user needs to provide `a` by `setACoeffs`, eta by `setShearViscosity`,
// and kappa by `setBulkViscosity`.  If `setBulkViscosity` is not called,
// kappa is set to zero.  The user must also call `setEBShearViscosity` to set
// viscosity on EB.  Optionally, `setEBBulkViscosity` can be used to set
// bulk viscosity on EB.
//
// The scalars alpha and beta can be set with `setScalar(Real, Real)`.  If
// they are not set, their default value is 1.

class MLEBTensorOp
    : public MLEBABecLap
{
public:

    MLEBTensorOp ();
    MLEBTensorOp (const Vector<Geometry>& a_geom,
                  const Vector<BoxArray>& a_grids,
                  const Vector<DistributionMapping>& a_dmap,
                  const LPInfo& a_info,
                  const Vector<EBFArrayBoxFactory const*>& a_factory);

    ~MLEBTensorOp () override;

    MLEBTensorOp (const MLEBTensorOp&) = delete;
    MLEBTensorOp (MLEBTensorOp&&) = delete;
    MLEBTensorOp& operator= (const MLEBTensorOp&) = delete;
    MLEBTensorOp& operator= (MLEBTensorOp&&) = delete;

    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info,
                 const Vector<EBFArrayBoxFactory const*>& a_factory);

    void setShearViscosity (int amrlev, const Array<MultiFab const*,AMREX_SPACEDIM>& eta,
                            Location a_beta_loc);
    void setShearViscosity (int amrlev, Real eta);

    void setBulkViscosity (int amrlev, const Array<MultiFab const*,AMREX_SPACEDIM>& kappa);
    void setBulkViscosity (int amrlev, Real kappa);

    void setEBShearViscosity (int amrlev, MultiFab const& eta);
    void setEBShearViscosity (int amrlev, Real eta);
    void setEBShearViscosityWithInflow (int amrlev, MultiFab const& eta, MultiFab const& eb_vel);

    void setEBBulkViscosity (int amrlev, MultiFab const& kappa);
    void setEBBulkViscosity (int amrlev, Real kappa);

    [[nodiscard]] int getNComp () const final { return AMREX_SPACEDIM; }

    [[nodiscard]] bool isCrossStencil () const final { return false; }
    [[nodiscard]] bool isTensorOp () const final { return true; }

    [[nodiscard]] bool needsUpdate () const final {
        return (m_needs_update || MLEBABecLap::needsUpdate());
    }
    void update () final {
        amrex::Abort("MLEBTensorOp: update TODO");
    }

    void prepareForSolve () final;
    [[nodiscard]] bool isSingular (int /*armlev*/) const final { return false; }
    [[nodiscard]] bool isBottomSingular () const final { return false; }

    void apply (int amrlev, int mglev, MultiFab& out, MultiFab& in, BCMode bc_mode,
                StateMode s_mode, const MLMGBndry* bndry=nullptr) const final;
    void compFlux (int amrlev, const Array<MultiFab*,AMREX_SPACEDIM>& fluxes,
                   MultiFab& sol, Location loc) const override;

    void compVelGrad (int amrlev, const Array<MultiFab*,AMREX_SPACEDIM>& grads,
                      MultiFab& sol, Location loc) const;

    void setBCoeffs (int amrlev, const Array<MultiFab const*,AMREX_SPACEDIM>& beta,
                     Location a_beta_loc) = delete;
    void setEBDirichlet (int amrlev, const MultiFab& phi, const MultiFab& beta) = delete;
    void setEBHomogDirichlet (int amrlev, const MultiFab& beta) = delete;

protected:

    bool m_needs_update = true;

    bool m_has_kappa = false;
    bool m_has_eb_kappa = false;
    Vector<Vector<Array<MultiFab,AMREX_SPACEDIM> > > m_kappa;
    Vector<Vector<MultiFab> > m_eb_kappa;
    mutable Vector<Vector<Array<MultiFab,AMREX_SPACEDIM> > > m_tauflux;

public: // for cuda

    void applyBCTensor (int amrlev, int mglev, MultiFab& vel,
                        BCMode bc_mode, StateMode s_mode, const MLMGBndry* bndry) const;
    void compCrossTerms(int amrlev, int mglev, MultiFab const& mf,
                        const MLMGBndry* bndry) const;
};

}

#endif

// add ghost cells to kappa
