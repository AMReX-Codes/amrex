#ifndef AMREX_ML_CURL_CURL_H_
#define AMREX_ML_CURL_CURL_H_
#include <AMReX_Config.H>

#include <AMReX_MLLinOp.H>
#include <AMReX_MLCurlCurl_K.H>

namespace amrex {

/**
 * \file AMReX_MLCurlCurl.H
 *
 * Vector-valued operator for `curl(alpha curl E) + beta E = rhs`, implemented
 * on a staggered grid with optional variable coefficients and PCG smoothing.
 */

/**
 * \brief curl (alpha curl E) + beta E = rhs
 *
 * Here E is an Array of 3 MultiFabs on the staggered grid.  The coefficient
 * \f$\alpha\f$ can be supplied either as a positive scalar (via setScalars)
 * or as a nodal MultiFab (via setAlpha), and \f$\beta\f$ can be either a
 * non-negative scalar or an edge-centered MultiFab (via setBeta).
 *
 * It's the caller's responsibility to make sure rhs has consistent nodal
 * data. If needed, one could call prepareRHS for this.
 *
 * The smoother is based on the 4-color Gauss-Seidel smoother of Li
 * et. al. 2020.  "An Efficient Preconditioner for 3-D Finite Difference
 * Modeling of the Electromagnetic Diffusion Process in the Frequency
 * Domain", IEEE Transactions on Geoscience and Remote Sensing, 58, 500-509.
 */
class MLCurlCurl
    : public MLLinOpT<Array<MultiFab,3> >
{
public:
    using MF        = Array<MultiFab,3>;
    using RT        = typename MLLinOpT<MF>::RT;
    using BCType    = typename MLLinOpT<MF>::BCType;
    using BCMode    = typename MLLinOpT<MF>::BCMode;
    using StateMode = typename MLLinOpT<MF>::StateMode;
    using Location  = typename MLLinOpT<MF>::Location;

    //! Construct an empty operator; call define() before solving.
    MLCurlCurl () = default;
    /**
     * \brief Convenience constructor forwarding to define().
     *
     * \param a_geom  Per-level geometries.
     * \param a_grids Per-level grids.
     * \param a_dmap  Distribution mappings.
     * \param a_info  LPInfo overrides (e.g., agglomeration).
     * \param a_coord Coordinate-system index (0 = Cartesian).
     */
    MLCurlCurl (const Vector<Geometry>& a_geom,
                const Vector<BoxArray>& a_grids,
                const Vector<DistributionMapping>& a_dmap,
                const LPInfo& a_info = LPInfo(),
                int a_coord = 0);

    /**
     * \brief Bind the operator to the AMR hierarchy.
     *
     * \param a_geom  Per-level geometries.
     * \param a_grids Per-level grids.
     * \param a_dmap  Distribution mappings.
     * \param a_info  LPInfo overrides.
     * \param a_coord Coordinate-system index (0 Cartesian, 1 RZ, etc.).
     */
    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info = LPInfo(),
                 int a_coord = 0);

    //! Set scalar coefficients.  This clears any previously supplied nodal alpha
    //! or edge-centered beta MultiFabs, so call setAlpha/setBeta afterwards if
    //! you still need spatially varying coefficients.
    void setScalars (RT a_alpha, RT a_beta) noexcept;

    //! This is needed only if there is variable beta coefficient.  Call this
    //! after setScalars if both APIs are used, because setScalars clears the
    //! cached beta MultiFabs.
    void setBeta (const Vector<Array<MultiFab const*,3>>& a_bcoefs);

    //! Set variable alpha defined on nodal MultiFab.  Call this after
    //! setScalars if both APIs are used, because setScalars clears the cached
    //! nodal alpha fields.
    void setAlpha (const Vector<MultiFab const*>& a_acoeffs);

    //! Synchronize RHS on nodal points. If the user can guarantee it, this
    //! function does not need to be called.
    void prepareRHS (Vector<MF*> const& rhs) const;

    //! Zero out Dirichlet nodes on (\p amrlev,\p mglev) (useful for GMRES).
    void setDirichletNodesToZero (int amrlev, int mglev, MF& a_mf) const override;

    [[nodiscard]] std::string name () const override {
        return std::string("curl of curl");
    }

    /**
     * \brief Toggle the PCG smoother; returns the previous value.
     *
     * \param flag True to enable the PCG smoother.
     * \return Previous state (true = PCG was enabled).
     */
    bool setUsePCG (bool flag) { return std::exchange(m_use_pcg, flag); }

    /**
     * \brief Provide level-specific BC data (compatible with MLLinOp interface).
     *
     * \param amrlev      AMR level index.
     * \param levelbcdata Optional BC data for this level.
     * \param robinbc_a   Optional Robin `a` coefficients defined on nodes.
     * \param robinbc_b   Optional Robin `b` coefficients defined on nodes.
     * \param robinbc_f   Optional Robin `f` coefficients defined on nodes.
     */
    void setLevelBC (int amrlev, const MF* levelbcdata,
                     const MF* robinbc_a = nullptr,
                     const MF* robinbc_b = nullptr,
                     const MF* robinbc_f = nullptr) override;

    /**
     * \brief Restrict edge-centered fields from fine to coarse MG level.
     *
     * \param amrlev AMR level index.
     * \param cmglev Coarse multigrid level receiving the data.
     * \param crse   Destination staggered field (3-component \p MF) on (\p amrlev,\p cmglev).
     * \param fine   Source staggered field on the next finer MG level.
     */
    void restriction (int amrlev, int cmglev, MF& crse, MF& fine) const override;

    /**
     * \brief Add the prolongation of coarse data onto the fine grid
     * (\p fine += prolong(\p crse)).
     *
     * \param amrlev AMR level index.
     * \param fmglev Fine multigrid level to populate.
     * \param fine   Destination MultiFab on (\p amrlev,\p fmglev).
     * \param crse   Source MultiFab on the next coarser MG level.
     */
    void interpolation (int amrlev, int fmglev, MF& fine, const MF& crse) const override;

    /**
     * \brief Apply the curl-curl operator with boundary conditions.
     *
     * \param amrlev AMR level index.
     * \param mglev  Multigrid level index.
     * \param out    Destination field updated with L(\p in).
     * \param in     Input field for the apply.
     * \param bc_mode Boundary-condition interpretation (solution vs correction).
     * \param s_mode  State interpretation (solution or correction).
     * \param bndry   Optional cached boundary metadata.
     */
    void apply (int amrlev, int mglev, MF& out, MF& in, BCMode bc_mode,
                StateMode s_mode, const MLMGBndryT<MF>* bndry=nullptr) const override;

    /**
     * \brief Run the Gauss-Seidel / PCG smoother for \p niter sweeps.
     *
     * \param amrlev AMR level index.
     * \param mglev  Multigrid level index.
     * \param sol    correction vector updated in place.
     * \param rhs    Right-hand side.
     * \param skip_fillboundary True if ghost nodes are already valid.
     * \param niter  Number of smoothing passes to execute.
     */
    void smooth (int amrlev, int mglev, MF& sol, const MF& rhs,
                 bool skip_fillboundary, int niter) const override;

    /**
     * \brief Compute residuals using solution BCs.
     *
     * \param amrlev      AMR level index.
     * \param resid       Destination residual.
     * \param x           Current solution.
     * \param b           Right-hand side.
     * \param crse_bcdata Optional coarse data to seed BCs.
     */
    void solutionResidual (int amrlev, MF& resid, MF& x, const MF& b,
                           const MF* crse_bcdata=nullptr) override;

    /**
     * \brief Compute residuals for correction solves.
     *
     * \param amrlev      AMR level index.
     * \param mglev       Multigrid level index.
     * \param resid       Destination residual.
     * \param x           Current correction.
     * \param b           Right-hand side.
     * \param bc_mode     Boundary-condition mode (usually homogeneous).
     * \param crse_bcdata Optional coarse data to inform BCs.
     */
    void correctionResidual (int amrlev, int mglev, MF& resid, MF& x,
                             const MF& b, BCMode bc_mode,
                             const MF* crse_bcdata=nullptr) override;

    [[nodiscard]] bool needsUpdate () const override {
        return (m_needs_update || MLLinOpT<Array<MultiFab,3>>::needsUpdate());
    }
    void update () override;

    //! Finish BC caches, mask updates, and metric setup prior to solves.
    void prepareForSolve () override;

    //! Prepare auxiliary data used by the multigrid preconditioner.
    void preparePrecond () override;

    [[nodiscard]] bool isSingular (int /*amrlev*/) const override { return false; }
    [[nodiscard]] bool isBottomSingular () const override { return false; }

    /**
     * \brief Dot-product helper.
     *
     * \param amrlev AMR level index.
     * \param mglev  Multigrid level index.
     * \param x      First vector.
     * \param y      Second vector.
     * \param local  True to limit the calculation to local data (skip MPI reduction).
     */
    RT xdoty (int amrlev, int mglev, const MF& x, const MF& y, bool local) const override;

    /**
     * \brief Infinity norm helper for vector-valued MultiFabs.
     *
     * \param amrlev AMR level index.
     * \param mf     Field whose infinity norm is computed.
     * \param local  True to avoid an MPI reduction.
     */
    [[nodiscard]] RT normInf (int amrlev, MF const& mf, bool local) const override;

    //! Average vector fields down the AMR hierarchy and synchronize interfaces using \p sol.
    void averageDownAndSync (Vector<MF>& sol) const override;

    [[nodiscard]] IntVect getNGrowVectRestriction () const override {
        return IntVect(1);
    }

    //! Allocate a hierarchy of \p MF s stored in \p mf with grow cells \p ng.
    void make (Vector<Vector<MF> >& mf, IntVect const& ng) const override;

    //! Allocate an \p MF on (\p amrlev,\p mglev) with grow cells \p ng.
    [[nodiscard]] MF make (int amrlev, int mglev, IntVect const& ng) const override;

    //! Create a shallow alias of \p mf (no allocation).
    [[nodiscard]] MF makeAlias (MF const& mf) const override;

    //! Allocate an \p MF on the next coarser MG level (\p mglev+1) of
    //! AMR level \p amrlev with grow cells \p ng.
    [[nodiscard]] MF makeCoarseMG (int amrlev, int mglev, IntVect const& ng) const override;

    //! Allocate an \p MF on the next coarser AMR level (\p famrlev-1)
    //! with grow cells \p ng.
    [[nodiscard]] MF makeCoarseAmr (int famrlev, IntVect const& ng) const override;

// public for cuda

#if (AMREX_SPACEDIM > 1)
    /**
     * \brief One color sweep of the 4-color smoother (2D/3D).
     *
     * \param amrlev AMR level index.
     * \param mglev  Multigrid level index.
     * \param sol    Correction updated in place.
     * \param rhs    Right-hand side for the sweep.
     * \param color  Color index (0-3) selecting the sub-lattice.
     */
    void smooth4 (int amrlev, int mglev, MF& sol, MF const& rhs, int color) const;
#else
    /**
     * \brief Two-color smoother for 1D cases.
     *
     * \param amrlev AMR level index.
     * \param mglev  Multigrid level index.
     * \param sol    Correction updated in place.
     * \param rhs    Right-hand side for the sweep.
     * \param color  Color index (0 or 1) selecting the sub-lattice.
     */
    void smooth1D (int amrlev, int mglev, MF& sol, MF const& rhs, int color) const;
#endif

    /**
     * \brief Residual computation helper on (\p amrlev,\p mglev).
     *
     * \param amrlev AMR level index.
     * \param mglev  Multigrid level index.
     * \param resid  Destination residual.
     * \param b      Right-hand side.
     */
    void compresid (int amrlev, int mglev, MF& resid, MF const& b) const;

    /**
     * \brief Apply physical BCs to \p mf for the chosen state type.
     *
     * \param amrlev AMR level index.
     * \param mglev  Multigrid level index.
     * \param mf     Field updated in place.
     * \param type   State interpretation (solution vs correction).
     */
    void applyPhysBC (int amrlev, int mglev, MultiFab& mf, CurlCurlStateType type) const;

private:

    void applyBC (int amrlev, int mglev, MF& in, CurlCurlStateType type) const;

    [[nodiscard]] iMultiFab const& getDotMask (int amrlev, int mglev, int idim) const;

    [[nodiscard]] CurlCurlDirichletInfo getDirichletInfo (int amrlev, int mglev) const;
    [[nodiscard]] CurlCurlSymmetryInfo getSymmetryInfo (int amrlev, int mglev) const;

    void update_lusolver ();

    void set_curvilinear_domain_bc ();
    void clearAlphaMultiFab ();
    void clearBetaMultiFab ();

    int m_coord = 0;

    RT m_alpha = std::numeric_limits<RT>::lowest();
    RT m_beta  = std::numeric_limits<RT>::lowest();

    Array<IntVect,3> m_etype
#if (AMREX_SPACEDIM == 3)
         {IntVect(0,1,1), IntVect(1,0,1), IntVect(1,1,0)};
#elif (AMREX_SPACEDIM == 2)
         {IntVect(0,1), IntVect(1,0), IntVect(1,1)};
#else
         {IntVect(0), IntVect(1), IntVect(1)};
#endif

    mutable Vector<Vector<Array<std::unique_ptr<iMultiFab>,3>>> m_dotmask;
    static constexpr int m_ncomp = 1;
    Vector<Vector<std::unique_ptr<Gpu::DeviceScalar
                                  <LUSolver<AMREX_SPACEDIM*2,RT>>>>> m_lusolver;
    Vector<Vector<Array<std::unique_ptr<MultiFab>,3>>> m_bcoefs;
    Vector<Vector<Array<std::unique_ptr<MultiFab>,3>>> m_acoefs;
    bool m_use_pcg = false;
    bool m_needs_update = true;
};

}

#endif
