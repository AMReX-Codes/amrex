#ifndef AMREX_ML_TENSOR_OP_H_
#define AMREX_ML_TENSOR_OP_H_

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
    virtual ~MLTensorOp ();

    MLTensorOp (const MLTensorOp&) = delete;
    MLTensorOp (MLTensorOp&&) = delete;
    MLTensorOp& operator= (const MLTensorOp&) = delete;
    MLTensorOp& operator= (MLTensorOp&&) = delete;

    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info = LPInfo(),
                 const Vector<FabFactory<FArrayBox> const*>& a_factory = {});

    void setShearViscosity (int amrlev, const Array<MultiFab const*,AMREX_SPACEDIM>& eta);
    void setBulkViscosity (int amrlev, const Array<MultiFab const*,AMREX_SPACEDIM>& kappa);

    virtual int getNComp () const final override { return AMREX_SPACEDIM; }

    virtual bool isCrossStencil () const final override { return false; }
    virtual bool isTensorOp () const final override { return true; }

    virtual bool needsUpdate () const final override {
        return (m_needs_update || MLABecLaplacian::needsUpdate());
    }
    virtual void update () final override {
        amrex::Abort("MLTensorOp: update TODO");
    }

    virtual void prepareForSolve () final override;
    virtual bool isSingular (int armlev) const final override { return false; }
    virtual bool isBottomSingular () const final override { return false; }

    virtual void apply (int amrlev, int mglev, MultiFab& out, MultiFab& in, BCMode bc_mode,
                        StateMode s_mode, const MLMGBndry* bndry=nullptr) const final override;

    virtual void compFlux (int amrlev, const Array<MultiFab*,AMREX_SPACEDIM>& fluxes,
                           MultiFab& sol, Location loc) const override;
  
protected:

    bool m_needs_update = true;

    bool m_has_kappa = false;
    Vector<Vector<Array<MultiFab,AMREX_SPACEDIM> > > m_kappa;

    void setBCoeffs (int amrlev, const Array<MultiFab const*,AMREX_SPACEDIM>& beta) = delete;

public: // for cuda

    void applyBCTensor (int amrlev, int mglev, MultiFab& vel,
                        BCMode bc_mode, const MLMGBndry* bndry) const;
};

}

#endif
