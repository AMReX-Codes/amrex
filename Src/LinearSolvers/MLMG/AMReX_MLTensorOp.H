#ifndef AMREX_ML_TENSOR_OP_H_
#define AMREX_ML_TENSOR_OP_H_

#include <AMReX_MLABecLaplacian.H>
#include <AMReX_Array.H>

namespace amrex {

// Tensor solver for high Reynolds flows with small gradient in shear viscosity.
//
// The system it solves is
//
//   a v + [-div dot (eta grad v) + (grad eta) (div dot v) - (grad eta) dot (grad v)^T] = rhs
//
// Note that what's inside [] is div dot tau + div dot ((kappa+eta/3) I (div dot v)),
// where tau is viscous tensor, and eta and kappa are shear and bulk viscosity, respectively.
// For incompressible flow, div dot v = 0, the system becomes
//
//   a v + div dot tau = rhs
//
// Known limitations: Using hypre or petsc as bottom solver is not supported.

class MLTensorOp
    : public MLABecLaplacian
{
public:

    MLTensorOp () {}
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

    // They are always 1 and 1.
    void setScalars (Real a, Real b) noexcept = delete;

    virtual int getNComp () const final override { return AMREX_SPACEDIM; }

    virtual bool needsUpdate () const final override {
        return (m_needs_update || MLCellABecLap::needsUpdate());
    }
    virtual void update () final override {
        amrex::Abort("MLTensorOp: update TODO");
    }

    virtual void prepareForSolve () final override;
    virtual bool isSingular (int armlev) const final override { return false; }
    virtual bool isBottomSingular () const final override { return false; }

    virtual void apply (int amrlev, int mglev, MultiFab& out, MultiFab& in, BCMode bc_mode,
                        StateMode s_mode, const MLMGBndry* bndry=nullptr) const final override;

protected:

    bool m_needs_update = true;

    Vector<MultiFab> m_gradeta;

};

}

#endif
