#ifndef AMREX_ML_TENSOR_OP_H_
#define AMREX_ML_TENSOR_OP_H_

#include <AMReX_MLABecLaplacian.H>
#include <AMReX_Array.H>

namespace amrex {

// Tensor solver for high Reynolds flows with near constant shear viscosity.
// 
// a v + [-div dot (eta grad v) + (grad eta) (div dot v) - (grad eta) dot (grad v)^T] = rhs
//
// Note that what's inside [] is div dot tau + div dot ((kappa+eta/3) I (div dot v)),
// where tau is viscous tensor, and eta and kappa are shear and bulk viscosity, respectively.
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

    virtual void solutionResidual (int amrlev, MultiFab& resid, MultiFab& x, const MultiFab& b,
                                   const MultiFab* crse_bcdata=nullptr) final override;

protected:

    bool m_needs_update = true;

    Vector<MultiFab> m_gradeta;

};

}

#endif
