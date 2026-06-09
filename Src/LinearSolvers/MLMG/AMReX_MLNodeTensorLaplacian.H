#ifndef AMREX_MLNODETENSORLAPLACIAN_H_
#define AMREX_MLNODETENSORLAPLACIAN_H_
#include <AMReX_Config.H>

#include <AMReX_MLNodeLinOp.H>

namespace amrex {

// del dot (sigma grad phi) = rhs
// where phi and rhs are nodal multifab, and sigma is a tensor constant.
// It is assumed that tensor is symmetric tensor of rank AMREX_SPACEDIM.
// Only periodic and Dirichlet are supported.

template <typename MF = amrex::MultiFab>
class MLNodeTensorLaplacianT
    : public MLNodeLinOpT<MF>
{
public:

    using FAB        = typename MF::fab_type;
    using RT         = typename MF::value_type;

    using BCMode     = typename MLNodeLinOpT<MF>::BCMode;
    using StateMode  = typename MLNodeLinOpT<MF>::StateMode;
    using CoarseningStrategy = typename MLNodeLinOpT<MF>::CoarseningStrategy;

#if (AMREX_SPACEDIM == 2)
    static constexpr int nelems = 3;  // number of unique elements in sigma
#else
    static constexpr int nelems = 6;
#endif

    MLNodeTensorLaplacianT () = default;
    MLNodeTensorLaplacianT (const Vector<Geometry>& a_geom,
                           const Vector<BoxArray>& a_grids,
                           const Vector<DistributionMapping>& a_dmap,
                           const LPInfo& a_info = LPInfo());
    ~MLNodeTensorLaplacianT () override = default;

    MLNodeTensorLaplacianT (const MLNodeTensorLaplacianT<MF>&) = delete;
    MLNodeTensorLaplacianT (MLNodeTensorLaplacianT<MF>&&) = delete;
    MLNodeTensorLaplacianT<MF>& operator= (const MLNodeTensorLaplacianT<MF>&) = delete;
    MLNodeTensorLaplacianT<MF>& operator= (MLNodeTensorLaplacianT<MF>&&) = delete;

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

    void restriction (int amrlev, int cmglev, MF& crse, MF& fine) const final;
    void interpolation (int amrlev, int fmglev, MF& fine, const MF& crse) const final;
    void averageDownSolutionRHS (int camrlev, MF& crse_sol, MF& crse_rhs,
                                 const MF& fine_sol, const MF& fine_rhs) final;

    void reflux (int crse_amrlev,
                 MF& res, const MF& crse_sol, const MF& crse_rhs,
                 MF& fine_res, MF& fine_sol, const MF& fine_rhs) const final;

    void smooth (int amrlev, int mglev, MF& sol, const MF& rhs,
                 bool skip_fillboundary, int niter) const final;

    void prepareForSolve () final;
    void Fapply (int amrlev, int mglev, MF& out, const MF& in) const final;
    void Fsmooth (int amrlev, int mglev, MF& sol, const MF& rhs) const final;
    void normalize (int amrlev, int mglev, MF& mf) const final;

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

using MLNodeTensorLaplacian = MLNodeTensorLaplacianT<>;

}

#endif
