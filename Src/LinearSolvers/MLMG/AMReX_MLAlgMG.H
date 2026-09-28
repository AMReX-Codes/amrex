#ifndef AMREX_ML_ALGMG_H_
#define AMREX_ML_ALGMG_H_
#include <AMReX_Config.H>

#include <AMReX_Geometry.H>
#include <AMReX_LayoutData.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>

#include <memory>

namespace amrex {

template <typename T> class AlgMG;
class MLNodeLinOp;

/**
 * \brief Algebraic system of one MG level of AMR level 0 of an MLLinOp,
 * solved with AlgMG.
 *
 * Built by MLLinOp::makeAlgMG under the ParallelContext communicator of that
 * level. The matrix is assembled once; solve() can be called repeatedly.
 */
class MLAlgMG
{
public:

    //! Nodal operator: rows are the nodes owned by this rank that are not
    //! Dirichlet nodes.
    MLAlgMG (int mglev, BoxArray const& grids, DistributionMapping const& dmap,
             Geometry const& geom, iMultiFab const& owner_mask,
             iMultiFab const& dirichlet_mask, MLNodeLinOp const& linop);

    /**
     * \brief Cell-centered operator alpha*a - beta*div(b grad): rows are all
     * cells except those of fully covered boxes.
     *
     * The a coefficient may be nullptr (zero); the b coefficients are
     * required. \p eb_bcoef is the EB Dirichlet b coefficient or nullptr.
     * \p bctype and \p bcl give the domain BC type and location per box and
     * face, and \p maxorder the order of the Dirichlet boundary stencil
     * (2 or 3).
     */
    MLAlgMG (int mglev, BoxArray const& grids, DistributionMapping const& dmap,
             Geometry const& geom, FabFactory<FArrayBox> const& factory,
             iMultiFab const* overset_mask, Real ascalar, Real bscalar,
             MultiFab const* acoef, Array<MultiFab const*,AMREX_SPACEDIM> const& bcoef,
             MultiFab const* eb_bcoef,
             LayoutData<GpuArray<int,2*AMREX_SPACEDIM>> const& bctype,
             LayoutData<GpuArray<Real,2*AMREX_SPACEDIM>> const& bcl, int maxorder);

    ~MLAlgMG ();
    MLAlgMG (MLAlgMG const&) = delete;
    MLAlgMG (MLAlgMG&&) = delete;
    MLAlgMG& operator= (MLAlgMG const&) = delete;
    MLAlgMG& operator= (MLAlgMG&&) = delete;

    /**
     * \brief Solve A x = b with zero initial guess.
     *
     * \p soln is overwritten. Converged when the residual 2-norm is below
     * max(reltol times the rhs 2-norm, abstol).
     */
    void solve (MultiFab& soln, MultiFab const& rhs, Real reltol, Real abstol, int maxiter);

    //! One AlgMG V-cycle with zero initial guess: the preconditioner form.
    void applyVcycle (MultiFab& soln, MultiFab const& rhs);

    //! The solver, for settings such as verbosity and Krylov acceleration.
    [[nodiscard]] AlgMG<Real>& solver () noexcept;

    struct Impl;

private:
    std::unique_ptr<Impl> m_impl;
};

}

#endif
