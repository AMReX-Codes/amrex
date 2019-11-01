#ifndef AMREX_MLNODETENSORLAPLACIAN_H_
#define AMREX_MLNODETENSORLAPLACIAN_H_

#include <AMReX_MLNodeLinOp.H>

namespace amrex {

// del dot (sigma grad phi) = rhs
// where phi and rhs are nodal multifab, and sigma is a tensor constant.
// It is assumed that tensor is symmetric tensor of rank AMREX_SPACEDIM.
// The user can set the tensor by calling setSigma or setBeta.
// For the case of setBeta, sigma = I - beta beta^T.

class MLNodeTensorLaplacian
    : public MLNodeLinOp
{
public:

#if (AMREX_SPACEDIM == 2)
    static constexpr int nelems = 3;  // number of unique elements in sigma
#else
    static constexpr int nelems = 6;
#endif

    MLNodeTensorLaplacian () noexcept {}
    MLNodeTensorLaplacian (const Vector<Geometry>& a_geom,
                           const Vector<BoxArray>& a_grids,
                           const Vector<DistributionMapping>& a_dmap,
                           const LPInfo& a_info = LPInfo());
    virtual ~MLNodeTensorLaplacian ();

    MLNodeTensorLaplacian (const MLNodeTensorLaplacian&) = delete;
    MLNodeTensorLaplacian (MLNodeTensorLaplacian&&) = delete;
    MLNodeTensorLaplacian& operator= (const MLNodeTensorLaplacian&) = delete;
    MLNodeTensorLaplacian& operator= (MLNodeTensorLaplacian&&) = delete;

    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info = LPInfo());

    virtual std::string name () const override { return std::string("MLNodeTensorLaplacian"); }

    virtual void restriction (int amrlev, int cmglev, MultiFab& crse, MultiFab& fine) const final override;
    virtual void interpolation (int amrlev, int fmglev, MultiFab& fine, const MultiFab& crse) const final override;
    virtual void averageDownSolutionRHS (int camrlev, MultiFab& crse_sol, MultiFab& crse_rhs,
                                         const MultiFab& fine_sol, const MultiFab& fine_rhs) final override;

    virtual void reflux (int crse_amrlev,
                         MultiFab& res, const MultiFab& crse_sol, const MultiFab& crse_rhs,
                         MultiFab& fine_res, MultiFab& fine_sol, const MultiFab& fine_rhs) const final override;

    virtual void prepareForSolve () final override;
    virtual void applyBC (int amrlev, int mglev, MultiFab& phi, BCMode bc_mode, StateMode s_mode,
                          bool skip_fillboundary=false) const final override;
    virtual void Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const final override;
    virtual void Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs) const final override;
    virtual void normalize (int amrlev, int mglev, MultiFab& mf) const final override;

    virtual void fixUpResidualMask (int amrlev, iMultiFab& resmsk) final override;
};

}

#endif
