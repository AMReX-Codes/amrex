#ifndef AMREX_ML_TENSOR_OP_H_
#define AMREX_ML_TENSOR_OP_H_
#include <AMReX_Config.H>

#include <AMReX_MLABecLaplacian.H>
#include <AMReX_Array.H>

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
// kappa is set to zero.
//
// The scalars alpha and beta can be set with `setScalar(Real, Real)`.  If
// they are not set, their default value is 1.

class MLTensorOp
    : public MLABecLaplacian
{
public:

    MLTensorOp ();
    MLTensorOp (const Vector<Geometry>& a_geom,
                const Vector<BoxArray>& a_grids,
                const Vector<DistributionMapping>& a_dmap,
                const LPInfo& a_info = LPInfo(),
                const Vector<FabFactory<FArrayBox> const*>& a_factory = {});
    MLTensorOp (const Vector<Geometry>& a_geom,
                const Vector<BoxArray>& a_grids,
                const Vector<DistributionMapping>& a_dmap,
                const Vector<iMultiFab const*>& a_overset_mask, // 1: unknown, 0: known
                const LPInfo& a_info = LPInfo(),
                const Vector<FabFactory<FArrayBox> const*>& a_factory = {});
    ~MLTensorOp () override = default;

    MLTensorOp (const MLTensorOp&) = delete;
    MLTensorOp (MLTensorOp&&) = delete;
    MLTensorOp& operator= (const MLTensorOp&) = delete;
    MLTensorOp& operator= (MLTensorOp&&) = delete;

    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info = LPInfo(),
                 const Vector<FabFactory<FArrayBox> const*>& a_factory = {});

    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const Vector<iMultiFab const*>& a_overset_mask,
                 const LPInfo& a_info = LPInfo(),
                 const Vector<FabFactory<FArrayBox> const*>& a_factory = {});

    void setShearViscosity (int amrlev, const Array<MultiFab const*,AMREX_SPACEDIM>& eta);
    void setShearViscosity (int amrlev, Real eta);
    void setBulkViscosity (int amrlev, const Array<MultiFab const*,AMREX_SPACEDIM>& kappa);
    void setBulkViscosity (int amrlev, Real kappa);

    [[nodiscard]] int getNComp () const final { return AMREX_SPACEDIM; }

    [[nodiscard]] bool isCrossStencil () const final { return false; }
    [[nodiscard]] bool isTensorOp () const final { return true; }

    [[nodiscard]] bool needsUpdate () const final {
        return (m_needs_update || MLABecLaplacian::needsUpdate());
    }
    void update () final {
        amrex::Abort("MLTensorOp: update TODO");
    }

    void prepareForSolve () final;
    [[nodiscard]] bool isSingular (int /*armlev*/) const final { return false; }
    [[nodiscard]] bool isBottomSingular () const final { return false; }

    void apply (int amrlev, int mglev, MultiFab& out, MultiFab& in, BCMode bc_mode,
                StateMode s_mode, const MLMGBndry* bndry=nullptr) const final;

    void compFlux (int amrlev, const Array<MultiFab*,AMREX_SPACEDIM>& fluxes,
                   MultiFab& sol, Location loc) const override;

    void compVelGrad (int amrlev, const Array<MultiFab*,AMREX_SPACEDIM>& fluxes,
                      MultiFab& sol, Location loc) const;

    void setBCoeffs (int amrlev, const Array<MultiFab const*,AMREX_SPACEDIM>& beta) = delete;

protected:

    bool m_needs_update = true;

    bool m_has_kappa = false;
    Vector<Vector<Array<MultiFab,AMREX_SPACEDIM> > > m_kappa;

public: // for cuda

    void applyBCTensor (int amrlev, int mglev, MultiFab& vel,
                        BCMode bc_mode, StateMode s_mode, const MLMGBndry* bndry
                        ) const;
};

}

#endif
