#ifndef AMREX_MLNODETENSORLAPLACIAN_H_
#define AMREX_MLNODETENSORLAPLACIAN_H_
#include <AMReX_Config.H>

#include <AMReX_MLNodeLinOp.H>

namespace amrex {

/**
 * \file AMReX_MLNodeTensorLaplacian.H
 *
 * Nodal tensor Laplacian where `sigma` is a symmetric tensor
 * configured via `setSigma` or `setBeta`.
 */

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

    /**
     * \brief Configure the symmetric tensor `sigma`.
     *
     * For 2D the array stores {xx, xy, yy}; for 3D it stores
     * {xx, xy, xz, yy, yz, zz}.
     *
     * \param a_sigma Components of the symmetric tensor.
     */
    void setSigma (Array<Real,nelems> const& a_sigma) noexcept;
    /**
     * \brief Convenience helper that sets `sigma = I - beta beta^T`.
     *
     * \param a_beta Vector `beta` used to build `sigma`.
     */
    void setBeta (Array<Real,AMREX_SPACEDIM> const& a_beta) noexcept;

    /**
     * \brief Bind the tensor operator to the AMR hierarchy (no EB support).
     *
     * \param a_geom    Per-level geometries.
     * \param a_grids   Per-level grids.
     * \param a_dmap    Distribution mappings.
     * \param a_info    Optional LPInfo overrides.
     */
    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info = LPInfo());

    [[nodiscard]] std::string name () const override { return std::string("MLNodeTensorLaplacian"); }

    /**
     * \brief Restrict nodal tensor data from fine to coarse on AMR level
     * \p amrlev / MG level \p cmglev.
     */
    void restriction (int amrlev, int cmglev, MultiFab& crse, MultiFab& fine) const final;
    /**
     * \brief Add the prolongation of coarse data onto the fine grid
     * (\p fine += prolong(\p crse)).
     */
    void interpolation (int amrlev, int fmglev, MultiFab& fine, const MultiFab& crse) const final;
    /**
     * \brief Average fine AMR data (\p fine_sol, \p fine_rhs) into coarse
     * arrays (\p crse_sol, \p crse_rhs) on level \p camrlev.
     */
    void averageDownSolutionRHS (int camrlev, MultiFab& crse_sol, MultiFab& crse_rhs,
                                 const MultiFab& fine_sol, const MultiFab& fine_rhs) final;

    /**
     * \brief Flux-register reflux between coarse/fine tensor solves; updates
     * \p res and \p fine_res using solutions and RHS on AMR index \p crse_amrlev.
     */
    void reflux (int crse_amrlev,
                 MultiFab& res, const MultiFab& crse_sol, const MultiFab& crse_rhs,
                 MultiFab& fine_res, MultiFab& fine_sol, const MultiFab& fine_rhs) const final;

    /**
     * \brief Tensor-aware smoother override for nodal data (updates \p sol
     * against \p rhs for \p niter passes; set \p skip_fillboundary when
     * ghosts are valid).
     */
    void smooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs,
                 bool skip_fillboundary, int niter) const final;

    //! Finalize tensor coefficients and BC caches before solving.
    void prepareForSolve () final;
    //! Apply the tensor Laplacian to \p in at (\p amrlev,\p mglev).
    void Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const final;
    /**
     * \brief Execute the nodal tensor smoother on (\p amrlev,\p mglev).
     *
     * \param amrlev AMR level index.
     * \param mglev  Multigrid level index.
     * \param sol    Correction updated in place.
     * \param rhs    Right-hand side for the sweep.
     */
    void Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs) const final;
    //! Divide \p mf by the diagonal of the operator (used by CG-family bottom solvers).
    void normalize (int amrlev, int mglev, MultiFab& mf) const final;

    //! Adjust residual masks \p resmsk so tensor Dirichlet EB nodes remain fixed.
    void fixUpResidualMask (int amrlev, iMultiFab& resmsk) final;

#if defined(AMREX_USE_HYPRE) && (AMREX_SPACEDIM > 1)
    /**
     * \brief Fill Hypre IJ matrix entries for a nodal tensor operator.
     *
     * \param mfi  MultiFab iterator.
     * \param gid  Global ID array.
     * \param lid  Local ID array.
     * \param ncols Number of nonzeros per row (output).
     * \param cols  Column indices (output).
     * \param mat   Matrix values (output).
     */
    void fillIJMatrix (MFIter const& mfi,
                       Array4<HypreNodeLap::AtomicInt const> const& gid,
                       Array4<int const> const& lid,
                       HypreNodeLap::Int* ncols,
                       HypreNodeLap::Int* cols,
                       Real* mat) const override;

    /**
     * \brief Fill Hypre RHS vector for nodal tensor solves.
     *
     * \param mfi MultiFab iterator.
     * \param lid Local ID array.
     * \param rhs RHS values (output).
     * \param bfab Source MultiFab data.
     */
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
