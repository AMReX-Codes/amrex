#ifndef AMREX_ML_NODE_LAPLACIAN_H_
#define AMREX_ML_NODE_LAPLACIAN_H_
#include <AMReX_Config.H>

#include <AMReX_iMultiFab.H>
#include <AMReX_MLNodeLinOp.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>

// Feature macro: when defined, MLNodeLaplacian::updateVelocity and
// MLNodeLaplacian::getFluxes correctly route through the anisotropic
// (harmonic-average) mknewu kernel when setMapped(true) was used.  Prior
// to this, both functions only used sigma component 0 (sigma_x) for all
// velocity/flux components, producing the wrong correction for
// non-identity mesh mappings.  EB is not covered by this feature (see
// comments in AMReX_MLNodeLaplacian_misc.cpp).
#define AMREX_MLNODELAP_HAS_MKNEWU_HA 1

namespace amrex {

/**
 * \file AMReX_MLNodeLaplacian.H
 *
 * Nodal Laplacian operator ``div(sigma grad phi)``.
 */

// del dot (sigma grad phi) = rhs
// where phi and rhs are nodal, and sigma is cell-centered.

/// \ingroup amrex_eb
class MLNodeLaplacian
    : public MLNodeLinOp
{
public :

    //! Construct an empty operator; call define() before solving.
    MLNodeLaplacian () = default;
    /**
     * \brief Convenience constructor without EB factories.
     *
     * \param a_geom        Per-level geometries.
     * \param a_grids       Per-level grids.
     * \param a_dmap        Distribution mappings.
     * \param a_info        LPInfo overrides.
     * \param a_factory     Optional nodal factories per level.
     * \param a_const_sigma Uniform conductivity.  Zero means variable
     *                      conductivity and setSigma() must be called.
     */
    MLNodeLaplacian (const Vector<Geometry>& a_geom,
                     const Vector<BoxArray>& a_grids,
                     const Vector<DistributionMapping>& a_dmap,
                     const LPInfo& a_info = LPInfo(),
                     const Vector<FabFactory<FArrayBox> const*>& a_factory = {},
                     Real  a_const_sigma = Real(0.0));
#ifdef AMREX_USE_EB
    //! Convenience constructor that accepts embedded-boundary factories.
    MLNodeLaplacian (const Vector<Geometry>& a_geom,
                     const Vector<BoxArray>& a_grids,
                     const Vector<DistributionMapping>& a_dmap,
                     const LPInfo& a_info,
                     const Vector<EBFArrayBoxFactory const*>& a_factory,
                     Real a_const_sigma = Real(0.0));
#endif
    ~MLNodeLaplacian () override = default;

    MLNodeLaplacian (const MLNodeLaplacian&) = delete;
    MLNodeLaplacian (MLNodeLaplacian&&) = delete;
    MLNodeLaplacian& operator= (const MLNodeLaplacian&) = delete;
    MLNodeLaplacian& operator= (MLNodeLaplacian&&) = delete;

    /**
     * \brief Bind the nodal operator (no EB data).
     *
     * \param a_geom        Per-level geometries.
     * \param a_grids       Per-level grids.
     * \param a_dmap        Distribution mappings.
     * \param a_info        LPInfo overrides.
     * \param a_factory     Optional nodal factories per level.
     * \param a_const_sigma Uniform conductivity.  Zero means variable
     *                      conductivity and setSigma() must be called.
     */
    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info = LPInfo(),
                 const Vector<FabFactory<FArrayBox> const*>& a_factory = {},
                 Real  a_const_sigma = Real(0.0));

#ifdef AMREX_USE_EB
    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info,
                 const Vector<EBFArrayBoxFactory const*>& a_factory,
                 Real a_const_sigma = Real(0.0));
#endif

    std::string name () const override { return std::string("MLNodeLaplacian"); }

    /**
     * \brief Enable cylindrical (RZ) metric terms.
     *
     * \param rz True activates RZ corrections.
     */
    void setRZCorrection (bool rz) noexcept { m_is_rz = rz; }

    /**
     * \brief Threshold for normalization (prevents dividing by tiny averages).
     *
     * \param t Minimum average magnitude allowed before rescaling.
     */
    void setNormalizationThreshold (Real t) noexcept { m_normalization_threshold = t; }

    /**
     * \brief Provide per-level conductivity (cell-centered MultiFab).
     * A multi-component input also triggers mapped-coordinate mode.
     *
     * \param amrlev AMR level receiving the coefficients.
     * \param a_sigma Cell-centered conductivity MultiFab.
     */
    void setSigma (int amrlev, const MultiFab& a_sigma);

    //! Compute `rhs = div(vel)` on nodal grids (used for projections).
    void compDivergence (const Vector<MultiFab*>& rhs, const Vector<MultiFab*>& vel);

    //! Assemble RHS including nodal/cell-centered forcing terms.
    void compRHS (const Vector<MultiFab*>& rhs, const Vector<MultiFab*>& vel,
                  const Vector<const MultiFab*>& rhnd,
                  const Vector<MultiFab*>& rhcc);

    //! Update velocities using the Poisson solution (projection correction).
    void updateVelocity (const Vector<MultiFab*>& vel, const Vector<MultiFab const*>& sol) const;

    //! Build coarse-level sync residuals for multilevel projections.
    void compSyncResidualCoarse (MultiFab& sync_resid, const MultiFab& phi,
                                 const MultiFab& vold, const MultiFab* rhcc,
                                 const BoxArray& fine_grids, const IntVect& ref_ratio);

    //! Build fine-level sync residuals for multilevel projections.
    void compSyncResidualFine (MultiFab& sync_resid, const MultiFab& phi, const MultiFab& vold,
                               const MultiFab* rhcc);

    /**
     * \brief Choose red-black Gauss-Seidel instead of Jacobi smoothing.
     *
     * \param flag True selects Gauss-Seidel.
     */
    void setGaussSeidel (bool flag) noexcept { m_use_gauss_seidel = flag; }
    /**
     * \brief Toggle harmonic averaging when coarsening sigma.
     *
     * \param flag True enables harmonic averaging.
     */
    void setHarmonicAverage (bool flag) noexcept { m_use_harmonic_average = flag; }

    /**
     * \brief Enable mapped-coordinate metrics (non-orthogonal grids).
     *
     * \param flag True enables mapped-coordinate support.
     */
    void setMapped (bool flag) noexcept { m_use_mapped = flag; }

    /**
     * \brief Pick between sigma-weighted or RAP coarsening (ignored if `sigma` constant).
     *
     * \param cs Desired coarsening strategy enum.
     */
    void setCoarseningStrategy (CoarseningStrategy cs) noexcept {
        if (m_const_sigma == Real(0.0)) { m_coarsening_strategy = cs; }
    }

    BottomSolver getDefaultBottomSolver () const final {
        return (m_coarsening_strategy == CoarseningStrategy::RAP) ?
            BottomSolver::bicgcg : BottomSolver::bicgstab;
    }

    //! Restrict nodal data from fine grids to coarse MG level \p cmglev.
    void restriction (int amrlev, int cmglev, MultiFab& crse, MultiFab& fine) const final;
    //! Add the prolongation of coarse data onto the fine grid (\p fine += prolong(\p crse)).
    void interpolation (int amrlev, int fmglev, MultiFab& fine, const MultiFab& crse) const final;
    //! Average solution/RHS pairs from fine AMR levels onto (\p camrlev).
    void averageDownSolutionRHS (int camrlev, MultiFab& crse_sol, MultiFab& crse_rhs,
                                         const MultiFab& fine_sol, const MultiFab& fine_rhs) final;

    /**
     * \brief Flux register reflux coupling coarse/fine nodal solves
     * (updates \p res using the solution and RHS arrays on both levels).
     */
    void reflux (int crse_amrlev,
                         MultiFab& res, const MultiFab& crse_sol, const MultiFab& crse_rhs,
                         MultiFab& fine_res, MultiFab& fine_sol, const MultiFab& fine_rhs) const final;

    //! Finalize metric terms, sigma averages, and BC caches before solving.
    void prepareForSolve () final;
    //! Apply the nodal Laplacian to \p in at (\p amrlev,\p mglev).
    void Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const final;
    /**
     * \brief Execute the smoother on (\p amrlev,\p mglev).
     *
     * \param amrlev AMR level index.
     * \param mglev  Multigrid level index.
     * \param sol    Correction updated in place.
     * \param rhs    RHS supplied to the smoother.
     */
    void Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs) const final;
    //! Divide \p mf by the diagonal of the operator (used by CG-family bottom solvers).
    void normalize (int amrlev, int mglev, MultiFab& mf) const final;

    //! Adjust residual masks so EB Dirichlet nodes remain constrained.
    void fixUpResidualMask (int amrlev, iMultiFab& resmsk) final;

    /**
     * \brief This overload is unused for nodal Laplacians (only cell-centered flux extraction is supported).
     *
     * \param a_flux Destination arrays (ignored).
     * \param a_sol  Solution hierarchy (ignored).
     * \param a_loc  Requested location (ignored).
     */
    void getFluxes (const Vector<Array<MultiFab*,AMREX_SPACEDIM> >& a_flux,
                            const Vector<MultiFab*>& a_sol,
                            Location a_loc) const final {
        amrex::ignore_unused(a_flux, a_sol, a_loc);
        amrex::Abort("MLNodeLaplacian::getFluxes: How did we get here?");
    }
    //! Compute face-centered fluxes given nodal solutions.
    void getFluxes (const Vector<MultiFab*>& a_flux,
                            const Vector<MultiFab*>& a_sol) const final;
    //! Remove Neumann contributions previously imposed on \p rhs.
    void unimposeNeumannBC (int amrlev, MultiFab& rhs) const final;
    //! Compute component-wise offsets needed for solvability enforcement.
    Vector<Real> getSolvabilityOffset (int amrlev, int mglev,
                                       MultiFab const& rhs) const override;
    //! Subtract solvability offsets from \p rhs.
    void fixSolvabilityByOffset (int amrlev, int mglev, MultiFab& rhs,
                                 Vector<Real> const& offset) const override;

    /**
     * \brief Unsupported flux-style gradient extractor (kept for interface completeness).
     *
     * \param amrlev AMR level index (ignored).
     * \param grad   Destination MultiFabs (ignored).
     * \param sol    Solution hierarchy (ignored).
     * \param loc    Requested location (ignored).
     */
    void compGrad (int amrlev, const Array<MultiFab*,AMREX_SPACEDIM>& grad,
                           MultiFab& sol, Location loc) const final {
        amrex::ignore_unused(amrlev, grad, sol, loc);
        amrex::Abort("MLNodeLaplacian::compGrad: How did we get here?");
    }
    /**
     * \brief Compute gradients of \p sol into nodal MultiFab \p grad on AMR level \p amrlev.
     *
     * \param amrlev AMR level index.
     * \param grad   Destination gradient MultiFab.
     * \param sol    Solution supplying nodal values.
     */
    void compGrad (int amrlev, MultiFab& grad, MultiFab& sol) const;

    //! Average conductivity down across all AMR and MG levels.
    void averageDownCoeffs ();
    //! Average coefficients from fine AMR level \p flev to coarse AMR level \p flev-1.
    void averageDownCoeffsToCoarseAmrLevel (int flev);
    //! Average coefficients across MG levels within AMR level \p amrlev.
    void averageDownCoeffsSameAmrLevel (int amrlev);

    //! Restrict interior-node RHS from fine to coarse when syncing.
    void restrictInteriorNodes (int camrlev, MultiFab& crhs, MultiFab& frhs) const;

    //! Fill ghost cells of \p sigma using \p geom for BC inference.
    void FillBoundaryCoeff (MultiFab& sigma, const Geometry& geom);

    //! Build nodal operator stencils reflecting current sigma and metrics.
    void buildStencil ();

#ifdef AMREX_USE_EB
    void buildIntegral ();
    void buildSurfaceIntegral ();

    // Tells the solver that EB boundaries have inflow specified by "eb_vel"
    void setEBInflowVelocity (int amrlev, const MultiFab& eb_vel);

#endif

#if defined(AMREX_USE_HYPRE) && (AMREX_SPACEDIM > 1)
    void fillIJMatrix (MFIter const& mfi,
                               Array4<HypreNodeLap::AtomicInt const> const& gid,
                               Array4<int const> const& lid,
                               HypreNodeLap::Int* ncols,
                               HypreNodeLap::Int* cols,
                               Real* mat) const override;

#ifdef AMREX_USE_GPU
    void fillIJMatrix_gpu (MFIter const& mfi,
                           Array4<HypreNodeLap::AtomicInt const> const& gid,
                           Array4<int const> const& lid,
                           HypreNodeLap::Int* ncols,
                           HypreNodeLap::Int* cols,
                           Real* mat) const;
#endif

    void fillIJMatrix_cpu (MFIter const& mfi,
                           Array4<HypreNodeLap::AtomicInt const> const& gid,
                           Array4<int const> const& lid,
                           HypreNodeLap::Int* ncols,
                           HypreNodeLap::Int* cols,
                           Real* mat) const;

    void fillRHS (MFIter const& mfi, Array4<int const> const& lid,
                  Real* rhs, Array4<Real const> const& bfab) const override;
#endif

protected:

    void resizeMultiGrid (int new_size) final;

private:

    int m_is_rz = 0;

    Real m_const_sigma = Real(0.0);
    Vector<Vector<Array<std::unique_ptr<MultiFab>,AMREX_SPACEDIM> > > m_sigma;
    Vector<Vector<std::unique_ptr<MultiFab> > > m_stencil;
    Vector<std::unique_ptr<MultiFab> > m_nosigma_stencil;
    Vector<Vector<Real> > m_s0_norm0;

    Real m_normalization_threshold = Real(1.e-8);

#ifdef AMREX_USE_EB
    // they could be MultiCutFab
    Vector<std::unique_ptr<MultiFab> > m_integral;
    bool m_integral_built = false;

    Vector<std::unique_ptr<MultiFab> > m_surface_integral;
    bool m_surface_integral_built = false;
    bool m_build_surface_integral = false;
    Vector<std::unique_ptr<MultiFab> > m_eb_vel_dot_n;
#endif

    bool m_use_gauss_seidel     = true;
    bool m_use_harmonic_average = false;
    bool m_use_mapped           = false;

    void checkPoint (std::string const& file_name) const final;
};

}

#endif
