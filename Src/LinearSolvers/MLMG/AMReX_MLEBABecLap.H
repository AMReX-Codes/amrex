#ifndef AMREX_MLEBABECLAP_H_
#define AMREX_MLEBABECLAP_H_
#include <AMReX_Config.H>

#include <AMReX_EBFabFactory.H>
#include <AMReX_MLCellABecLap.H>
#include <AMReX_Array.H>
#include <limits>

namespace amrex {

/**
 * \file AMReX_MLEBABecLap.H
 *
 * Embedded-boundary-aware ABec Laplacian operator that extends `MLCellABecLap`
 * with EB-specific coefficient and boundary handling.
 */

template <typename T>
struct MLMGABCEBTag;

// (alpha * a - beta * (del dot b grad)) phi

/// \ingroup amrex_eb
/// \brief Embedded-boundary ABec Laplacian (`alpha a - beta div(b grad phi)`).
class MLEBABecLap
    : public MLCellABecLap
{
public:

    //! Construct an empty operator; call define() before use.
    MLEBABecLap () = default;
    /**
     * \brief Convenience constructor that forwards to define().
     *
     * \param a_geom    Per-level geometries.
     * \param a_grids   Per-level grids.
     * \param a_dmap    Distribution mappings.
     * \param a_info    LPInfo overrides (agglomeration, coarsening, etc.).
     * \param a_factory EB factories per level.
     * \param a_ncomp   Number of components handled.
     */
    MLEBABecLap (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info,
                 const Vector<EBFArrayBoxFactory const*>& a_factory,
                 int a_ncomp = 1);

    MLEBABecLap (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const Vector<iMultiFab const*>& a_overset_mask, // 1: unknown, 0: known
                 const LPInfo& a_info,
                 const Vector<EBFArrayBoxFactory const*>& a_factory,
                 int a_ncomp = 1);

    ~MLEBABecLap () override;

    MLEBABecLap (const MLEBABecLap&) = delete;
    MLEBABecLap (MLEBABecLap&&) = delete;
    MLEBABecLap& operator= (const MLEBABecLap&) = delete;
    MLEBABecLap& operator= (MLEBABecLap&&) = delete;

    /**
     * \brief Bind the EB operator to the AMR hierarchy.
     *
     * \param a_geom    Per-level geometries.
     * \param a_grids   Per-level grids.
     * \param a_dmap    Distribution mappings.
     * \param a_info    LPInfo overrides.
     * \param a_factory EB factories per level.
     * \param a_ncomp   Number of components handled.
     */
    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info,
                 const Vector<EBFArrayBoxFactory const*>& a_factory,
                 int a_ncomp = 1);

    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const Vector<iMultiFab const*>& a_overset_mask,
                 const LPInfo& a_info,
                 const Vector<EBFArrayBoxFactory const*>& a_factory,
                 int a_ncomp = 1);

    //! Interprets `phi` on EB centroids instead of cell centers when set.
    void setPhiOnCentroid ();

    /**
     * \brief Set constant scalar coefficients (see MLCellABecLap semantics).
     *
     * \param a Scalar multiplying the zero-th order term.
     * \param b Scalar multiplying the diffusion term.
     */
    void setScalars (Real a, Real b);
    /**
     * \brief Provide per-cell `a` coefficients for AMR level \p amrlev.
     *
     * \param amrlev AMR level index.
     * \param alpha  Cell-centered MultiFab containing `a`.
     */
    void setACoeffs (int amrlev, const MultiFab& alpha);
    /**
     * \brief Provide a constant `a` coefficient.
     *
     * \param amrlev AMR level index.
     * \param alpha  Scalar value assigned everywhere.
     */
    void setACoeffs (int amrlev, Real alpha);

    /**
     * \brief Provide spatially varying `b` coefficients and their location.
     *
     * \param amrlev    AMR level index.
     * \param beta      Array of face-centered MultiFabs.
     * \param a_beta_loc Location of the coefficients (face center/centroid).
     */
    void setBCoeffs (int amrlev, const Array<MultiFab const*,AMREX_SPACEDIM>& beta,
                     Location a_beta_loc);
    void setBCoeffs (int amrlev, const Array<MultiFab const*,AMREX_SPACEDIM>& beta)
        {setBCoeffs (amrlev, beta, Location::FaceCenter);}

    /**
     * \brief Provide constant scalar `b` coefficients.
     *
     * \param amrlev AMR level index.
     * \param beta   Scalar applied to every face.
     */
    void setBCoeffs (int amrlev, Real beta);
    /**
     * \brief Provide component-specific constant `b` coefficients.
     *
     * \param amrlev AMR level index.
     * \param beta   Vector of scalars (one per component).
     */
    void setBCoeffs (int amrlev, Vector<Real> const& beta);

    /**
     * \brief Specify Dirichlet values on EB faces using `phi`.
     *
     * \param amrlev AMR level index.
     * \param phi    Prescribed values at EB face centroids.
     * \param beta   Face coefficients.
     */
    void setEBDirichlet      (int amrlev, const MultiFab& phi, const MultiFab& beta);
    /**
     * \brief Specify Dirichlet values on EB faces.
     *
     * \param amrlev AMR level index.
     * \param phi    Prescribed values at EB face centroids.
     * \param beta   Single scalar coefficient applied to every EB face.
     */
    void setEBDirichlet      (int amrlev, const MultiFab& phi, Real beta);
    /**
     * \brief Specify Dirichlet values on EB faces.
     *
     * \param amrlev AMR level index.
     * \param phi    Prescribed values at EB face centroids.
     * \param beta   Per-component scalar coefficients.
     */
    void setEBDirichlet      (int amrlev, const MultiFab& phi, Vector<Real> const& beta);

    /**
     * \brief Mark EB faces as homogeneous Dirichlet (phi=0) with optional coefficients.
     *
     * \param amrlev AMR level index.
     * \param beta   Coefficients applied to EB faces.
     */
    void setEBHomogDirichlet (int amrlev,                      const MultiFab& beta);
    /**
     * \brief Mark EB faces as homogeneous Dirichlet (phi=0).
     *
     * \param amrlev AMR level index.
     * \param beta   Single scalar coefficient applied to every EB face.
     */
    void setEBHomogDirichlet (int amrlev,                      Real beta);
    /**
     * \brief Mark EB faces as homogeneous Dirichlet (phi=0).
     *
     * \param amrlev AMR level index.
     * \param beta   Per-component scalar coefficients.
     */
    void setEBHomogDirichlet (int amrlev,                      Vector<Real> const& beta);

    int getNComp () const override { return m_ncomp; }

    bool needsUpdate () const override {
        return (m_needs_update || MLCellABecLap::needsUpdate());
    }
    void update () override;

    //! EB-aware factory used when the solver needs auxiliary storage on (\p amrlev,\p mglev).
    std::unique_ptr<FabFactory<FArrayBox> > makeFactory (int amrlev, int mglev) const final;

    //! EB discretization uses a non-cross stencil.
    bool isCrossStencil () const override { return false; }

    //! Apply physical BCs to \p in on (\p amrlev,\p mglev) using the supplied BC metadata.
    void applyBC (int amrlev, int mglev, MultiFab& in, BCMode bc_mode, StateMode s_mode,
                          const MLMGBndry* bndry=nullptr, bool skip_fillboundary=false) const final;
    //! Apply the discrete operator to \p in on (\p amrlev,\p mglev) and store the result in \p out.
    void apply (int amrlev, int mglev, MultiFab& out, MultiFab& in, BCMode bc_mode,
                        StateMode s_mode, const MLMGBndry* bndry=nullptr) const override;
    //! Compute gradients from the solution \p sol on level \p amrlev at location \p loc.
    void compGrad (int amrlev, const Array<MultiFab*,AMREX_SPACEDIM>& grad,
                           MultiFab& sol, Location loc) const final;

    //! Finalize EB masks, averages, and singularity flags prior to solving.
    void prepareForSolve () override;
    //! True if AMR level \p amrlev is singular.
    bool isSingular (int amrlev) const override { return m_is_singular[amrlev]; }
    //! True if the bottom level is singular.
    bool isBottomSingular () const override { return m_is_singular[0]; }
    //! Apply the operator on (\p amrlev,\p mglev) and write to \p out.
    void Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const final;
    //! Execute a smoothing sweep on (\p amrlev,\p mglev).  \p redblack selects the red (0) or black (1) half of the grid.
    void Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs, int redblack) const final;
    //! Compute EB-aware fluxes for \p sol on the tilebox described by \p mfi, storing to \p flux at \p loc (respecting \p face_only).
    void FFlux (int amrlev, const MFIter& mfi,
                        const Array<FArrayBox*,AMREX_SPACEDIM>& flux,
                        const FArrayBox& sol, Location loc,
                        int face_only=0) const final;

    //! Divide \p mf by the diagonal of the operator (used by CG-family bottom solvers).
    void normalize (int amrlev, int mglev, MultiFab& mf) const final;

    //! Scalar applied to the zeroth-order term.
    Real getAScalar () const final { return m_a_scalar; }
    //! Scalar applied to the diffusion term.
    Real getBScalar () const final { return m_b_scalar; }
    //! Access the `a` coefficient MultiFab.
    MultiFab const* getACoeffs (int amrlev, int mglev) const final
        { return &(m_a_coeffs[amrlev][mglev]); }
    //! Access the `b` coefficients for (\p amrlev,\p mglev).
    Array<MultiFab const*,AMREX_SPACEDIM> getBCoeffs (int amrlev, int mglev) const final
        { return amrex::GetArrOfConstPtrs(m_b_coeffs[amrlev][mglev]); }

    std::unique_ptr<MLLinOp> makeNLinOp (int /*grid_size*/) const final {
        amrex::Abort("MLABecLaplacian::makeNLinOp: Not implemented");
        return std::unique_ptr<MLLinOp>{};
    }

    void averageDownCoeffsSameAmrLevel (int amrlev, Vector<MultiFab>& a,
                                        Vector<Array<MultiFab,AMREX_SPACEDIM> >& b,
                                        const Vector<MultiFab*>& b_eb);

    /**
     * \brief Restrict fine-grid data onto coarse storage.
     *
     * \param amrlev AMR level index requested by the base interface.
     * \param cmglev Multigrid level index on that AMR level.
     * \param crse   Destination coarse MultiFab.
     * \param fine   Source fine MultiFab (modified if ghost fills are needed).
     */
    void restriction (int amrlev, int cmglev, MultiFab& crse, MultiFab& fine) const final;

    /**
     * \brief Add the prolongation of coarse data onto the fine grid
     * (\p fine += prolong(\p crse)).
     *
     * \param amrlev AMR level index.
     * \param fmglev Fine MG level index.
     * \param fine   Destination fine MultiFab (updated in place).
     * \param crse   Source coarse MultiFab.
     */
    void interpolation (int amrlev, int fmglev, MultiFab& fine, const MultiFab& crse) const final;

    //! Average solution/RHS from fine AMR levels into (\p camrlev) storage.
    void averageDownSolutionRHS (int camrlev, MultiFab& crse_sol, MultiFab& crse_rhs,
                                         const MultiFab& fine_sol, const MultiFab& fine_rhs) final;

    /**
     * \brief Extract EB face fluxes.
     *
     * \param a_flux Destination MultiFabs (one per level) that receive EB fluxes.
     * \param a_sol  Solution hierarchy sampled at EB faces.
     */
    void getEBFluxes (const Vector<MultiFab*>& a_flux,
                              const Vector<MultiFab*>& a_sol) const override;

    //! Apply Robin BC terms to the stored coefficients when requested.
    void applyRobinBCTermsCoeffs ();

#if defined(AMREX_USE_HYPRE) && (AMREX_SPACEDIM > 1)
    //! Instantiate a Hypre backend using interface selector \p hypre_interface.
    [[nodiscard]] std::unique_ptr<Hypre> makeHypre (Hypre::Interface hypre_interface) const override;
#endif

#if defined(AMREX_USE_PETSC) && (AMREX_SPACEDIM > 1)
    [[nodiscard]] std::unique_ptr<PETScABecLap> makePETSc () const override;
#endif

    //! Scalar multiplier for the `a` term (set via setScalars()).
    Real m_a_scalar = std::numeric_limits<Real>::quiet_NaN();
    Real m_b_scalar = std::numeric_limits<Real>::quiet_NaN();
    Vector<Vector<MultiFab> > m_a_coeffs;
    Vector<Vector<Array<MultiFab,AMREX_SPACEDIM> > > m_b_coeffs;

    bool m_scalars_set = false;
    bool m_acoef_set = false;

protected:

    int m_ncomp = 1;

    Location m_beta_loc; // Location of coefficients: face centers or face centroids
    Location m_phi_loc;  // Location of solution variable: cell centers or cell centroids

    Vector<Vector<iMultiFab> > m_cc_mask;

    Vector<std::unique_ptr<MultiFab> > m_eb_phi;
    Vector<Vector<std::unique_ptr<MultiFab> > > m_eb_b_coeffs;

    Vector<int> m_is_singular;

    mutable int m_is_eb_inhomog;

    //
    // functions
    //
    bool isEBDirichlet   () const noexcept { return m_eb_phi[0] != nullptr; }

    //! Average coefficients down across all AMR and MG levels.
    void averageDownCoeffs ();
    //! Average inhomogeneous EB Dirichlet values down across AMR levels.
    void averageDownEBPhi ();
    //! Average coefficients from fine AMR level \p flev to coarse AMR level \p flev-1.
    void averageDownCoeffsToCoarseAmrLevel (int flev);

    [[nodiscard]] bool supportRobinBC () const noexcept override { return true; }

    mutable Vector<Vector<TagVector<MLMGABCEBTag<RT>>>> m_eb_bc_tags;

private:

    bool m_needs_update = true;

    void define_ab_coeffs ();
};

/// \cond DOXYGEN_IGNORE

// Tag for applyBC EB version
template <typename T>
struct MLMGABCEBTag {
    amrex::Array4<int const> mask;
    amrex::Array4<T const> area;
    T bcloc;
    amrex::Box bx;
    amrex::BoundCond bctype;
    int blen;
    int comp;
    int dir;
    int is_eb;
    int side;
    int local_index;

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Box const& box() const noexcept { return bx; }
};

/// \endcond

}

#endif
