#ifndef AMREX_MLNODEABECLAPLACIAN_H_
#define AMREX_MLNODEABECLAPLACIAN_H_
#include <AMReX_Config.H>

#include <AMReX_MLNodeLinOp.H>

namespace amrex {

/**
 * \file AMReX_MLNodeABecLaplacian.H
 *
 * Nodal ABec Laplacian ``alpha a phi - beta div(b grad phi)``.
 */

// (alpha * a - beta * (del dot b grad)) phi = rhs
// a, phi and rhs are nodal. b is cell-centered.

class MLNodeABecLaplacian
    : public MLNodeLinOp
{
public:

    MLNodeABecLaplacian () = default;
    MLNodeABecLaplacian (const Vector<Geometry>& a_geom,
                         const Vector<BoxArray>& a_grids,
                         const Vector<DistributionMapping>& a_dmap,
                         const LPInfo& a_info = LPInfo(),
                         const Vector<FabFactory<FArrayBox> const*>& a_factory = {});
    ~MLNodeABecLaplacian () override = default;

    MLNodeABecLaplacian (const MLNodeABecLaplacian&) = delete;
    MLNodeABecLaplacian (MLNodeABecLaplacian&&) = delete;
    MLNodeABecLaplacian& operator= (const MLNodeABecLaplacian&) = delete;
    MLNodeABecLaplacian& operator= (MLNodeABecLaplacian&&) = delete;

    /**
     * \brief Bind the nodal ABec Laplacian to the AMR hierarchy.
     *
     * \param a_geom    Per-level geometries.
     * \param a_grids   Per-level grids.
     * \param a_dmap    Distribution mappings.
     * \param a_info    Optional LPInfo overrides.
     * \param a_factory Optional nodal factories per level.
     */
    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info = LPInfo(),
                 const Vector<FabFactory<FArrayBox> const*>& a_factory = {});

    std::string name () const override { return std::string("MLNodeABecLaplacian"); }

    //! Set constant scalars `a` and `b` for the ABec operator.
    void setScalars (Real a, Real b) {
        m_a_scalar = a;
        m_b_scalar = b;
    }

    /**
     * \brief Provide scalar nodal `a` coefficients for AMR level \p amrlev.
     *
     * \param amrlev AMR level receiving the scalar value.
     * \param a_acoef Scalar nodal coefficient.
     */
    void setACoeffs (int amrlev, Real a_acoef);
    /**
     * \brief Provide nodal `a` coefficients from a MultiFab.
     *
     * \param amrlev AMR level receiving the data.
     * \param a_acoef Nodal MultiFab describing `a`.
     */
    void setACoeffs (int amrlev, const MultiFab& a_acoef);

    /**
     * \brief Provide scalar cell-centered `b` coefficients.
     *
     * \param amrlev AMR level receiving the scalar value.
     * \param a_bcoef Scalar `b` coefficient.
     */
    void setBCoeffs (int amrlev, Real a_bcoef);
    /**
     * \brief Provide cell-centered `b` coefficients from a MultiFab.
     *
     * \param amrlev AMR level receiving the data.
     * \param a_bcoef Cell-centered MultiFab describing `b`.
     */
    void setBCoeffs (int amrlev, const MultiFab& a_bcoef);

    //! Apply the nodal ABec operator on (\p amrlev,\p mglev) to \p in and store the result in \p out.
    void Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const final;
    /**
     * \brief Execute a smoothing sweep on (\p amrlev,\p mglev).
     *
     * \param amrlev AMR level index.
     * \param mglev  Multigrid level index.
     * \param sol    Correction updated in place.
     * \param rhs    Right-hand side for the sweep.
     */
    void Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs) const final;

    //! Enforce EB Dirichlet nodes by modifying the residual mask \p resmsk.
    void fixUpResidualMask (int amrlev, iMultiFab& resmsk) final;

    bool isSingular (int /*amrlev*/) const final { return false; }
    bool isBottomSingular () const final { return false; }

    //! Restrict nodal data from fine to coarse MG levels (AMR level \p amrlev, coarse MG index \p cmglev).
    void restriction (int amrlev, int cmglev, MultiFab& crse, MultiFab& fine) const final;
    //! Add the prolongation of coarse data onto the fine grid (\p fine += prolong(\p crse)).
    void interpolation (int amrlev, int fmglev, MultiFab& fine, const MultiFab& crse) const final;
    //! Average solution/RHS pairs from fine AMR level to coarse storage.
    void averageDownSolutionRHS (int camrlev, MultiFab& crse_sol, MultiFab& crse_rhs,
                                 const MultiFab& fine_sol, const MultiFab& fine_rhs) final;

    //! Flux register reflux step coupling coarse and fine nodal solves.
    void reflux (int crse_amrlev,
                 MultiFab& res, const MultiFab& crse_sol, const MultiFab& crse_rhs,
                 MultiFab& fine_res, MultiFab& fine_sol, const MultiFab& fine_rhs) const final;

    //! Finalize coefficient averages and BC caches prior to solving.
    void prepareForSolve () final;

    [[nodiscard]] bool needsUpdate () const final { return m_needs_update; }

    //! Refresh coefficients (averages them down) and clear the dirty flag.
    void update () final;

    //! Average `a` and `b` coefficients down across all AMR and MG levels.
    void averageDownCoeffs ();
    //! Average coefficients from fine AMR level \p flev to coarse AMR level \p flev-1.
    void averageDownCoeffsToCoarseAmrLevel (int flev);
    //! Average coefficients across MG levels within the same AMR level.
    void averageDownCoeffsSameAmrLevel (int amrlev);

private:

    bool m_needs_update = true;

    Real m_a_scalar = std::numeric_limits<Real>::quiet_NaN();
    Real m_b_scalar = std::numeric_limits<Real>::quiet_NaN();
    Vector<Vector<MultiFab> > m_a_coeffs;
    Vector<Vector<MultiFab> > m_b_coeffs;
};

}

#endif
