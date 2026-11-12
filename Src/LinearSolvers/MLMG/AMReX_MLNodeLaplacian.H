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
// non-identity mesh mappings.  EB and RAP coarsening are not covered by
// this feature (see comments in AMReX_MLNodeLaplacian_misc.cpp).
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
template <typename MF>
class MLNodeLaplacianT
    : public MLNodeLinOpT<MF>
{
public :

    friend class HypreMLNodeLap;

    using FAB = typename MLNodeLinOpT<MF>::FAB;
    using RT  = typename MLNodeLinOpT<MF>::RT;
    using BCMode    = typename MLNodeLinOpT<MF>::BCMode;
    using StateMode = typename MLNodeLinOpT<MF>::StateMode;
    using Location  = typename MLNodeLinOpT<MF>::Location;
    using CoarseningStrategy = typename MLNodeLinOpT<MF>::CoarseningStrategy;

    //! Construct an empty operator; call define() before solving.
    MLNodeLaplacianT () = default;
    /**
     * \brief Convenience constructor without EB factories.
     *
     * \param a_geom        Per-level geometries.
     * \param a_grids       Per-level grids.
     * \param a_dmap        Distribution mappings.
     * \param a_info        LPInfo overrides.
     * \param a_factory     Optional nodal factories per level.
     * \param a_const_sigma Uniform conductivity.  Zero selects variable
     *                      conductivity; setSigma() must then supply nonzero coefficients.
     */
    MLNodeLaplacianT (const Vector<Geometry>& a_geom,
                      const Vector<BoxArray>& a_grids,
                      const Vector<DistributionMapping>& a_dmap,
                      const LPInfo& a_info = LPInfo(),
                      const Vector<FabFactory<FAB> const*>& a_factory = {},
                      RT  a_const_sigma = RT(0.0));
#ifdef AMREX_USE_EB
    //! Convenience constructor that accepts embedded-boundary factories.
    MLNodeLaplacianT (const Vector<Geometry>& a_geom,
                      const Vector<BoxArray>& a_grids,
                      const Vector<DistributionMapping>& a_dmap,
                      const LPInfo& a_info,
                      const Vector<EBFArrayBoxFactory const*>& a_factory,
                      RT a_const_sigma = RT(0.0));
#endif
    ~MLNodeLaplacianT () override = default;

    MLNodeLaplacianT (const MLNodeLaplacianT&) = delete;
    MLNodeLaplacianT (MLNodeLaplacianT&&) = delete;
    MLNodeLaplacianT& operator= (const MLNodeLaplacianT&) = delete;
    MLNodeLaplacianT& operator= (MLNodeLaplacianT&&) = delete;

    /**
     * \brief Bind the nodal operator (no EB data).
     *
     * \param a_geom        Per-level geometries.
     * \param a_grids       Per-level grids.
     * \param a_dmap        Distribution mappings.
     * \param a_info        LPInfo overrides.
     * \param a_factory     Optional nodal factories per level.
     * \param a_const_sigma Uniform conductivity.  Zero selects variable
     *                      conductivity; setSigma() must then supply nonzero coefficients.
     */
    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info = LPInfo(),
                 const Vector<FabFactory<FAB> const*>& a_factory = {},
                 RT  a_const_sigma = RT(0.0));

#ifdef AMREX_USE_EB
    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info,
                 const Vector<EBFArrayBoxFactory const*>& a_factory,
                 RT a_const_sigma = RT(0.0));
#endif

    [[nodiscard]] std::string name () const override { return std::string("MLNodeLaplacian"); }

    [[nodiscard]] bool supportsMixedPrecision () const override {
        return typeid(*this) == typeid(MLNodeLaplacianT<MF>);
    }

    /**
     * \brief Enable cylindrical (RZ) metric terms.
     *
     * \param rz True activates RZ corrections.
     */
    void setRZCorrection (bool rz) noexcept {
        m_is_rz = rz;
        this->template spForward<MLNodeLaplacianT>([&] (auto* sp) { sp->setRZCorrection(rz); });
    }

    /**
     * \brief Threshold for normalization (prevents dividing by tiny averages).
     *
     * \param t Minimum average magnitude allowed before rescaling.
     */
    void setNormalizationThreshold (RT t) noexcept {
        if (m_normalization_threshold != t) {
            m_normalization_threshold = t;
            m_needs_update = true;
        }
        this->template spForward<MLNodeLaplacianT>([&] (auto* sp) {
            sp->setNormalizationThreshold(static_cast<float>(t)); });
    }

    /**
     * \brief Provide per-level conductivity (cell-centered MultiFab).
     * A multi-component input also triggers mapped-coordinate mode.
     *
     * \param amrlev AMR level receiving the coefficients.
     * \param a_sigma Cell-centered MultiFab containing nonzero conductivity coefficients.
     */
    void setSigma (int amrlev, const MF& a_sigma);

    //! Mixed precision: send sigma, scaled, to the single precision copy.
    //! Public for CUDA.
    void spForwardSigma ();

    //! Compute `rhs = div(vel)` on nodal grids (used for projections).
    void compDivergence (const Vector<MF*>& rhs, const Vector<MF*>& vel);

    //! Assemble RHS including nodal/cell-centered forcing terms.
    void compRHS (const Vector<MF*>& rhs, const Vector<MF*>& vel,
                  const Vector<const MF*>& rhnd,
                  const Vector<MF*>& rhcc);

    //! Update velocities using the Poisson solution (projection correction).
    void updateVelocity (const Vector<MF*>& vel, const Vector<MF const*>& sol) const;

    //! Build coarse-level sync residuals for multilevel projections.
    void compSyncResidualCoarse (MF& sync_resid, const MF& phi,
                                 const MF& vold, const MF* rhcc,
                                 const BoxArray& fine_grids, const IntVect& ref_ratio);

    //! Build fine-level sync residuals for multilevel projections.
    void compSyncResidualFine (MF& sync_resid, const MF& phi, const MF& vold,
                               const MF* rhcc);

    /**
     * \brief Choose red-black Gauss-Seidel instead of Jacobi smoothing.
     *
     * \param flag True selects Gauss-Seidel.
     */
    void setGaussSeidel (bool flag) noexcept {
        m_use_gauss_seidel = flag;
        this->template spForward<MLNodeLaplacianT>([&] (auto* sp) { sp->setGaussSeidel(flag); });
    }
    /**
     * \brief Toggle harmonic averaging when coarsening sigma.
     *
     * \param flag True enables harmonic averaging.
     */
    void setHarmonicAverage (bool flag) noexcept {
        if (m_use_harmonic_average != flag) {
            m_use_harmonic_average = flag;
            m_needs_update = true;
        }
        this->template spForward<MLNodeLaplacianT>([&] (auto* sp) { sp->setHarmonicAverage(flag); });
    }

    /**
     * \brief Enable mapped-coordinate metrics (non-orthogonal grids).
     *
     * \param flag True enables mapped-coordinate support.
     */
    void setMapped (bool flag) noexcept {
        if (m_use_mapped != flag) {
            m_use_mapped = flag;
            m_needs_update = true;
        }
        this->template spForward<MLNodeLaplacianT>([&] (auto* sp) { sp->setMapped(flag); });
    }

    /**
     * \brief Pick between sigma-weighted or RAP coarsening (ignored if `sigma` constant).
     *
     * \param cs Desired coarsening strategy enum.
     */
    void setCoarseningStrategy (CoarseningStrategy cs) noexcept {
        if (m_const_sigma == RT(0.0) && this->m_coarsening_strategy != cs) {
            this->m_coarsening_strategy = cs;
            this->m_masks_built = false;
            this->m_precond_weight_mask.clear();
            m_needs_update = true;
        }
        this->template spForward<MLNodeLaplacianT>([&] (auto* sp) { sp->setCoarseningStrategy(cs); });
    }

    [[nodiscard]] BottomSolver getDefaultBottomSolver () const final {
        return (this->m_coarsening_strategy == CoarseningStrategy::RAP) ?
            BottomSolver::bicgcg : BottomSolver::bicgstab;
    }

    //! Restrict nodal data from fine grids to coarse MG level \p cmglev.
    void restriction (int amrlev, int cmglev, MF& crse, MF& fine) const final;
    //! Add the prolongation of coarse data onto the fine grid (\p fine += prolong(\p crse)).
    void interpolation (int amrlev, int fmglev, MF& fine, const MF& crse) const final;
    //! Average solution/RHS pairs from fine AMR levels onto (\p camrlev).
    void averageDownSolutionRHS (int camrlev, MF& crse_sol, MF& crse_rhs,
                                         const MF& fine_sol, const MF& fine_rhs) final;

    /**
     * \brief Set the coarse/fine nodes of fine AMR level \p famrlev (Dirichlet
     * mask set, not on a Dirichlet domain face) to the interpolated coarse
     * solution, so that the composite solution does not depend on the
     * initial guess there.  Known nodes of an overset mask are kept.
     *
     * \param famrlev  Fine AMR level index.
     * \param fine_sol Fine solution, modified at the coarse/fine nodes.
     * \param crse_sol Coarse solution.
     */
    void interpCoarseFineSolution (int famrlev, MF& fine_sol,
                                   const MF& crse_sol) const final;

    /**
     * \brief Flux register reflux coupling coarse/fine nodal solves
     * (updates \p res using the solution and RHS arrays on both levels).
     */
    void reflux (int crse_amrlev,
                         MF& res, const MF& crse_sol, const MF& crse_rhs,
                         MF& fine_res, MF& fine_sol, const MF& fine_rhs) const final;

    //! Finalize metric terms, sigma averages, and BC caches before solving.
    void prepareForSolve () final;
    //! True if sigma-derived data need to be refreshed before reuse.
    [[nodiscard]] bool needsUpdate () const final { return m_needs_update; }
    //! Refresh sigma averages/stencils and clear the dirty flag.
    void update () final;
    //! Apply the nodal Laplacian to \p in at (\p amrlev,\p mglev).
    void Fapply (int amrlev, int mglev, MF& out, const MF& in) const final;
    /**
     * \brief Execute the smoother on (\p amrlev,\p mglev).
     *
     * \param amrlev AMR level index.
     * \param mglev  Multigrid level index.
     * \param sol    Correction updated in place.
     * \param rhs    RHS supplied to the smoother.
     */
    void Fsmooth (int amrlev, int mglev, MF& sol, const MF& rhs) const final;
    //! Divide \p mf by the diagonal of the operator (used by CG-family bottom solvers).
    void normalize (int amrlev, int mglev, MF& mf) const final;

    //! Adjust residual masks so EB Dirichlet nodes remain constrained.
    void fixUpResidualMask (int amrlev, iMultiFab& resmsk) final;

    /**
     * \brief This overload is unused for nodal Laplacians (only cell-centered flux extraction is supported).
     *
     * \param a_flux Destination arrays (ignored).
     * \param a_sol  Solution hierarchy (ignored).
     * \param a_loc  Requested location (ignored).
     */
    void getFluxes (const Vector<Array<MF*,AMREX_SPACEDIM> >& a_flux,
                            const Vector<MF*>& a_sol,
                            Location a_loc) const final {
        amrex::ignore_unused(a_flux, a_sol, a_loc);
        amrex::Abort("MLNodeLaplacian::getFluxes: How did we get here?");
    }
    //! Compute face-centered fluxes given nodal solutions.
    void getFluxes (const Vector<MF*>& a_flux,
                            const Vector<MF*>& a_sol) const final;
    //! Remove Neumann contributions previously imposed on \p rhs.
    void unimposeNeumannBC (int amrlev, MF& rhs) const final;
    //! Compute component-wise offsets needed for solvability enforcement.
    Vector<RT> getSolvabilityOffset (int amrlev, int mglev,
                                     MF const& rhs) const override;
    //! Subtract solvability offsets from \p rhs.
    void fixSolvabilityByOffset (int amrlev, int mglev, MF& rhs,
                                 Vector<RT> const& offset) const override;

    /**
     * \brief Unsupported flux-style gradient extractor (kept for interface completeness).
     *
     * \param amrlev AMR level index (ignored).
     * \param grad   Destination MultiFabs (ignored).
     * \param sol    Solution hierarchy (ignored).
     * \param loc    Requested location (ignored).
     */
    void compGrad (int amrlev, const Array<MF*,AMREX_SPACEDIM>& grad,
                           MF& sol, Location loc) const final {
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
    void compGrad (int amrlev, MF& grad, MF& sol) const;

    //! Average conductivity down across all AMR and MG levels.
    void averageDownCoeffs ();
    //! Average coefficients from fine AMR level \p flev to coarse AMR level \p flev-1.
    void averageDownCoeffsToCoarseAmrLevel (int flev);
    //! Average coefficients across MG levels within AMR level \p amrlev.
    void averageDownCoeffsSameAmrLevel (int amrlev);

    //! Restrict interior-node RHS from fine to coarse when syncing.
    void restrictInteriorNodes (int camrlev, MF& crhs, MF& frhs) const;

    //! Fill ghost cells of \p sigma using \p geom for BC inference.
    void FillBoundaryCoeff (MF& sigma, const Geometry& geom);

    //! Build nodal operator stencils reflecting current sigma and metrics.
    void buildStencil ();

#ifdef AMREX_USE_EB
    void buildIntegral ();
    void buildSurfaceIntegral ();

    // Tells the solver that EB boundaries have inflow specified by "eb_vel"
    void setEBInflowVelocity (int amrlev, const MF& eb_vel);

#endif

#if (AMREX_SPACEDIM > 1)
#if defined(AMREX_USE_HYPRE)
    void fillIJMatrix (MFIter const& mfi,
                       Array4<HypreNodeLap::AtomicInt const> const& gid,
                       Array4<int const> const& lid,
                       HypreNodeLap::Int* ncols,
                       HypreNodeLap::Int* cols,
                       Real* mat) const override;
#endif

    [[nodiscard]] bool supportsAlgMG () const override { return true; }
    void fillAlgMatrix (int mglev, MFIter const& mfi,
                        Array4<Long const> const& gid,
                        Array4<int const> const& lid,
                        Long* ncols, Long* cols, Real* mat) const override;

    void fillRHS (int mglev, MFIter const& mfi, Array4<int const> const& lid,
                  Real* rhs, Array4<Real const> const& bfab) const override;

    // Public for nvcc. AlgGid is the global row id type: Long for AlgMG,
    // HypreNodeLap::AtomicInt for hypre.
    template <typename AlgInt, typename AlgGid>
    void fillMatrix_doit (int mglev, MFIter const& mfi,
                          Array4<AlgGid const> const& gid,
                          Array4<int const> const& lid,
                          AlgInt* ncols, AlgInt* cols, Real* mat) const;
#ifdef AMREX_USE_GPU
    template <typename AlgInt, typename AlgGid>
    void fillMatrix_gpu (int mglev, MFIter const& mfi,
                         Array4<AlgGid const> const& gid,
                         Array4<int const> const& lid,
                         AlgInt* ncols, AlgInt* cols, Real* mat) const;
#endif
    template <typename AlgInt, typename AlgGid>
    void fillMatrix_cpu (int mglev, MFIter const& mfi,
                         Array4<AlgGid const> const& gid,
                         Array4<int const> const& lid,
                         AlgInt* ncols, AlgInt* cols, Real* mat) const;
#endif

protected:

    void spCopySettings () override
    {
        MLNodeLinOpT<MF>::spCopySettings();
        if (m_is_rz) { setRZCorrection(true); }
        setNormalizationThreshold(m_normalization_threshold);
        setGaussSeidel(m_use_gauss_seidel);
        setHarmonicAverage(m_use_harmonic_average);
        setMapped(m_use_mapped);
        setCoarseningStrategy(this->m_coarsening_strategy);
    }

    void resizeMultiGrid (int new_size) final;

    //! Mixed precision: the copy's sigma is m_sp_sigma_scale times ours and
    //! its lengths m_sp_length_scale times ours.  The RZ terms do not depend
    //! on the length.
    void spUpdateScale () noexcept {
        if (this->m_sp_linop) {
            this->m_sp_scale = m_sp_sigma_scale
                * static_cast<RT>(std::pow(this->m_sp_length_scale, -2));
        }
    }

#ifdef AMREX_USE_EB
    //! EB factory of AMR level amrlev, or nullptr.  EB data are double
    //! precision, so the float operator never has one.
    [[nodiscard]] EBFArrayBoxFactory const* ebFactory (int amrlev) const
    {
        if constexpr (std::is_same_v<MF,MultiFab>) {
            return dynamic_cast<EBFArrayBoxFactory const*>(this->m_factory[amrlev][0].get());
        } else {
            amrex::ignore_unused(amrlev);
            return nullptr;
        }
    }

    //! Volume fraction of EB factory f, or nullptr.
    [[nodiscard]] static MF const* ebVolFrac (EBFArrayBoxFactory const* f)
    {
        if constexpr (std::is_same_v<MF,MultiFab>) {
            return (f) ? &(f->getVolFrac()) : nullptr;
        } else {
            amrex::ignore_unused(f);
            return nullptr;
        }
    }
#endif

private:

    int m_is_rz = 0;

    RT m_const_sigma = RT(0.0);
    RT m_sp_sigma_scale = RT(1.0);
    Vector<Vector<Array<std::unique_ptr<MF>,AMREX_SPACEDIM> > > m_sigma;
    Vector<Vector<std::unique_ptr<MF> > > m_stencil;
    Vector<std::unique_ptr<MF> > m_nosigma_stencil;
    Vector<Vector<RT> > m_s0_norm0;

    RT m_normalization_threshold = RT(1.e-8);

#ifdef AMREX_USE_EB
    // they could be MultiCutFab
    Vector<std::unique_ptr<MF> > m_integral;
    bool m_integral_built = false;

    Vector<std::unique_ptr<MF> > m_surface_integral;
    bool m_surface_integral_built = false;
    bool m_build_surface_integral = false;
    Vector<std::unique_ptr<MF> > m_eb_vel_dot_n;
#endif

    bool m_use_gauss_seidel     = true;
    bool m_use_harmonic_average = false;
    bool m_use_mapped           = false;
    bool m_needs_update         = true;

    void checkPoint (std::string const& file_name) const final;
};

extern template class MLNodeLaplacianT<MultiFab>;
#ifdef AMREX_MLMG_HAS_FLOAT
extern template class MLNodeLaplacianT<fMultiFab>;
#endif

/// \cond DOXYGEN_IGNORE
namespace detail {
    // One empty function per AMReX_MLNodeLaplacian_*.cpp, for mlndlap_link_parts.
    void mlndlap_anchor_sync ();
    void mlndlap_anchor_sten ();
    void mlndlap_anchor_misc ();
    void mlndlap_anchor_matrix ();
#ifdef AMREX_USE_EB
    void mlndlap_anchor_eb ();
#endif
    void mlndlap_link_parts ();
}
/// \endcond

//! Double precision nodal Laplacian.
using MLNodeLaplacian = MLNodeLaplacianT<MultiFab>;

}

#endif
