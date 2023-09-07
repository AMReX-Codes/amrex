#ifndef AMREX_MLNODETENSORLAPLACIAN_H_
#define AMREX_MLNODETENSORLAPLACIAN_H_
#include <AMReX_Config.H>

#include <AMReX_MLNodeLinOp.H>

namespace amrex {

// del dot (sigma grad phi) = rhs
// where phi and rhs are nodal multifab, and sigma is a tensor constant.
// It is assumed that tensor is symmetric tensor of rank AMREX_SPACEDIM.
// Only periodic and Dirichlet are supported.

class MLNodeTensorLaplacian
    : public MLNodeLinOp
{
public:

#if (AMREX_SPACEDIM == 2)
    static constexpr int nelems = 3;  // number of unique elements in sigma
#else
    static constexpr int nelems = 6;
#endif

    MLNodeTensorLaplacian () = default;
    MLNodeTensorLaplacian (const Vector<Geometry>& a_geom,
                           const Vector<BoxArray>& a_grids,
                           const Vector<DistributionMapping>& a_dmap,
                           const LPInfo& a_info = LPInfo());
    ~MLNodeTensorLaplacian () override = default;

    MLNodeTensorLaplacian (const MLNodeTensorLaplacian&) = delete;
    MLNodeTensorLaplacian (MLNodeTensorLaplacian&&) = delete;
    MLNodeTensorLaplacian& operator= (const MLNodeTensorLaplacian&) = delete;
    MLNodeTensorLaplacian& operator= (MLNodeTensorLaplacian&&) = delete;

    // The user can set the tensor by calling either setSigma or setBeta.
    // For 2d, a_sigma contains components xx, xy, and yy.
    // For 3d, a_sigma contains components xx, xy, xz, yy, yz, and zz.
    void setSigma (Array<Real,nelems> const& a_sigma) noexcept;
    // sigma = I - beta beta^T.
    void setBeta (Array<Real,AMREX_SPACEDIM> const& a_beta) noexcept;

    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info = LPInfo());

    [[nodiscard]] std::string name () const override { return std::string("MLNodeTensorLaplacian"); }

    void restriction (int amrlev, int cmglev, MultiFab& crse, MultiFab& fine) const final;
    void interpolation (int amrlev, int fmglev, MultiFab& fine, const MultiFab& crse) const final;
    void averageDownSolutionRHS (int camrlev, MultiFab& crse_sol, MultiFab& crse_rhs,
                                 const MultiFab& fine_sol, const MultiFab& fine_rhs) final;

    void reflux (int crse_amrlev,
                 MultiFab& res, const MultiFab& crse_sol, const MultiFab& crse_rhs,
                 MultiFab& fine_res, MultiFab& fine_sol, const MultiFab& fine_rhs) const final;

    void smooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs,
                 bool skip_fillboundary=false) const final;

    void prepareForSolve () final;
    void Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const final;
    void Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs) const final;
    void normalize (int amrlev, int mglev, MultiFab& mf) const final;

    void fixUpResidualMask (int amrlev, iMultiFab& resmsk) final;

#if defined(AMREX_USE_HYPRE) && (AMREX_SPACEDIM > 1)
    void fillIJMatrix (MFIter const& mfi,
                       Array4<HypreNodeLap::AtomicInt const> const& gid,
                       Array4<int const> const& lid,
                       HypreNodeLap::Int* ncols,
                       HypreNodeLap::Int* cols,
                       Real* mat) const override;

    void fillRHS (MFIter const& mfi,
                  Array4<int const> const& lid,
                  Real* rhs,
                  Array4<Real const> const& bfab) const override;
#endif

private:

    GpuArray<Real,nelems> scaledSigma (int amrlev, int mglev) const noexcept;

    GpuArray<Real,nelems> m_sigma;
    mutable int m_redblack = 0;
};

}

#endif
