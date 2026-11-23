#ifndef AMREX_ML_ABECLAP_DRIFT_H_
#define AMREX_ML_ABECLAP_DRIFT_H_
#include <AMReX_Config.H>

#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLABecLap_Drift_K.H>

namespace amrex {

/**
 * \file AMReX_MLABecLapDrift.H
 *
 * Cell-centered multi-level operator for
 * ``alpha a - beta div(b grad phi) - C div(c phi)``, with an upwind drift
 * term on top of MLABecLaplacian.
 */

/**
 * \ingroup amrex_solver_mlmg
 * \brief `(A \alpha - B \nabla \cdot \beta \nabla - C \nabla \cdot c) \phi = f`
 *
 * MLABecLaplacian with an upwind drift term: the flux through a face is c
 * times the cell value upwind of it (the ghost cell for the outside). Stable
 * when C has the sign opposite to B's; not symmetric, so the bottom solver
 * must be bicgstab or the smoother. c must not enter through a Neumann face.
 * With C = 0 or without setDriftCoeffs, this is MLABecLaplacian. See the
 * Linear Solvers chapter of the documentation for the restrictions.
 */
template <typename MF>
class MLABecLapDriftT
    : public MLABecLaplacianT<MF>
{
public:

    using FAB = typename MF::fab_type;
    using RT  = typename MF::value_type;

    using Location = typename MLLinOpT<MF>::Location;

    MLABecLapDriftT () = default;
    MLABecLapDriftT (const Vector<Geometry>& a_geom,
                     const Vector<BoxArray>& a_grids,
                     const Vector<DistributionMapping>& a_dmap,
                     const LPInfo& a_info = LPInfo(),
                     const Vector<FabFactory<FAB> const*>& a_factory = {},
                     int a_ncomp = 1);

    ~MLABecLapDriftT () override = default;

    MLABecLapDriftT (const MLABecLapDriftT<MF>&) = delete;
    MLABecLapDriftT (MLABecLapDriftT<MF>&&) = delete;
    MLABecLapDriftT<MF>& operator= (const MLABecLapDriftT<MF>&) = delete;
    MLABecLapDriftT<MF>& operator= (MLABecLapDriftT<MF>&&) = delete;

    //! See MLABecLaplacianT::define. Removes the drift coefficients.
    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info = LPInfo(),
                 const Vector<FabFactory<FAB> const*>& a_factory = {},
                 int a_ncomp = 1);

    //! Set the scalars A and B, with C = 0 (no drift). Not virtual: call it
    //! on this class, not through a reference to MLABecLaplacian.
    template <typename T1, typename T2>
    requires (std::is_convertible_v<T1, typename MF::value_type> &&
              std::is_convertible_v<T2, typename MF::value_type>)
    void setScalars (T1 a, T2 b) noexcept;

    //! Set the scalars A, B and C. C != 0 needs setDriftCoeffs.
    template <typename T1, typename T2, typename T3>
    requires (std::is_convertible_v<T1, typename MF::value_type> &&
              std::is_convertible_v<T2, typename MF::value_type> &&
              std::is_convertible_v<T3, typename MF::value_type>)
    void setScalars (T1 a, T2 b, T3 c) noexcept;

    /**
     * \brief Set the face velocities c of the drift term on AMR level
     * \p amrlev (call it on every level; see setScalars for C).
     *
     * \param amrlev AMR level.
     * \param c      Face velocities per direction, with one or getNComp()
     *               components.
     */
    template <FabArrayType AMF>
    requires (std::is_convertible_v<typename AMF::value_type, typename MF::value_type>)
    void setDriftCoeffs (int amrlev, const Array<AMF const*,AMREX_SPACEDIM>& c);

    //! Relaxation factor of the smoother of the operator with drift (default 1)
    void setDriftRelaxation (RT omega) noexcept { m_drift_omega = omega; }

    //! Whether the operator has the drift term of setDriftCoeffs
    [[nodiscard]] bool hasDrift () const noexcept {
        return m_c_scalar != RT(0.0) && !m_c_coeffs.empty();
    }

    [[nodiscard]] bool needsUpdate () const override {
        return m_drift_needs_update || MLABecLaplacianT<MF>::needsUpdate();
    }
    void update () override;
    void prepareForSolve () override;

    //! The matrix-based solvers assemble the ABec part only: not with drift
    [[nodiscard]] bool supportsAlgMG () const override {
        return !hasDrift() && MLABecLaplacianT<MF>::supportsAlgMG();
    }
    //! The operator with drift is not symmetric
    [[nodiscard]] bool isSymmetric () const override { return !hasDrift(); }
#if defined(AMREX_USE_HYPRE) && (AMREX_SPACEDIM > 1)
    [[nodiscard]] std::unique_ptr<Hypre> makeHypre (Hypre::Interface hypre_interface) const override {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!hasDrift(),
            "MLABecLapDrift: no hypre bottom solver with drift");
        return MLABecLaplacianT<MF>::makeHypre(hypre_interface);
    }
#endif
#if defined(AMREX_USE_PETSC) && (AMREX_SPACEDIM > 1)
    [[nodiscard]] std::unique_ptr<PETScABecLap> makePETSc () const override {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!hasDrift(),
            "MLABecLapDrift: no PETSc bottom solver with drift");
        return MLABecLaplacianT<MF>::makePETSc();
    }
#endif

    void Fapply (int amrlev, int mglev, MF& out, const MF& in) const override;
    void Fsmooth (int amrlev, int mglev, MF& sol, const MF& rhs, int redblack) const override;
    void normalize (int amrlev, int mglev, MF& mf) const override;

    using MLABecLaplacianT<MF>::FFlux;
    //! The fluxes of MLABecLaplacian plus the drift flux `-C c \phi_{upwind}`
    void FFlux (int amrlev, const MFIter& mfi,
                const Array<FAB*,AMREX_SPACEDIM>& flux,
                const FAB& sol, Location loc,
                int face_only=0) const override;

    // public for cuda

    /**
     * \brief The drift through the boundary of AMR level \p amrlev: the faces
     * with a ghost cell from the boundary conditions (mask not covered), so
     * not the periodic ones. Aborts if it enters through a Neumann face.
     *
     * \return The number of components whose drift crosses the boundary.
     */
    [[nodiscard]] int boundaryDrift (int amrlev) const;

private:

    RT m_c_scalar = RT(0.0);
    //! Face velocities on every AMR and MG level, empty without setDriftCoeffs
    Vector<Vector<Array<MF,AMREX_SPACEDIM> > > m_c_coeffs;
    RT m_drift_omega = RT(1.0);
    bool m_drift_needs_update = true;

    //! Averaging down and singular flags of the drift; after the update of
    //! MLABecLaplacian
    void updateDrift ();

    //! The restrictions of the drift term
    void checkDrift () const;

    //! Whether the coarse/fine boundary of AMR level \p amrlev is Neumann:
    //! only on a level 0 that does not cover the domain; on the finer levels
    //! the ghost cells come from the coarse data
    [[nodiscard]] bool cfNeumann (int amrlev) const {
        return amrlev == 0 && this->m_needs_coarse_data_for_bc
            && this->m_coarse_fine_bc_type == LinOpBCType::Neumann;
    }

    //! The coefficients of the operator with drift on (\p amrlev, \p mglev)
    [[nodiscard]] MLABecLapDriftCoeffsMA<RT> driftArrays (int amrlev, int mglev) const;
};

template <typename MF>
MLABecLapDriftT<MF>::MLABecLapDriftT (const Vector<Geometry>& a_geom,
                                      const Vector<BoxArray>& a_grids,
                                      const Vector<DistributionMapping>& a_dmap,
                                      const LPInfo& a_info,
                                      const Vector<FabFactory<FAB> const*>& a_factory,
                                      int a_ncomp)
{
    define(a_geom, a_grids, a_dmap, a_info, a_factory, a_ncomp);
}

template <typename MF>
void
MLABecLapDriftT<MF>::define (const Vector<Geometry>& a_geom,
                             const Vector<BoxArray>& a_grids,
                             const Vector<DistributionMapping>& a_dmap,
                             const LPInfo& a_info,
                             const Vector<FabFactory<FAB> const*>& a_factory,
                             int a_ncomp)
{
    MLABecLaplacianT<MF>::define(a_geom, a_grids, a_dmap, a_info, a_factory, a_ncomp);
    m_c_scalar = RT(0.0);
    m_c_coeffs.clear();
    m_drift_needs_update = true;
}

template <typename MF>
template <typename T1, typename T2>
requires (std::is_convertible_v<T1, typename MF::value_type> &&
          std::is_convertible_v<T2, typename MF::value_type>)
void
MLABecLapDriftT<MF>::setScalars (T1 a, T2 b) noexcept
{
    MLABecLaplacianT<MF>::setScalars(a, b);
    m_c_scalar = RT(0.0);
}

template <typename MF>
template <typename T1, typename T2, typename T3>
requires (std::is_convertible_v<T1, typename MF::value_type> &&
          std::is_convertible_v<T2, typename MF::value_type> &&
          std::is_convertible_v<T3, typename MF::value_type>)
void
MLABecLapDriftT<MF>::setScalars (T1 a, T2 b, T3 c) noexcept
{
    MLABecLaplacianT<MF>::setScalars(a, b);
    m_c_scalar = RT(c);
}

template <typename MF>
template <FabArrayType AMF>
requires (std::is_convertible_v<typename AMF::value_type, typename MF::value_type>)
void
MLABecLapDriftT<MF>::setDriftCoeffs (int amrlev, const Array<AMF const*,AMREX_SPACEDIM>& c)
{
    const int ncomp = this->getNComp();
    if (m_c_coeffs.empty()) {
        m_c_coeffs.resize(this->m_num_amr_levels);
        for (int alev = 0; alev < this->m_num_amr_levels; ++alev) {
            m_c_coeffs[alev].resize(this->m_num_mg_levels[alev]);
            for (int mglev = 0; mglev < this->m_num_mg_levels[alev]; ++mglev) {
                for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                    const BoxArray& ba = amrex::convert(this->m_grids[alev][mglev],
                                                        IntVect::TheDimensionVector(idim));
                    m_c_coeffs[alev][mglev][idim].define
                        (ba, this->m_dmap[alev][mglev], ncomp, 0, MFInfo(),
                         *(this->m_factory[alev][mglev]));
                    m_c_coeffs[alev][mglev][idim].setVal(RT(0.0));
                }
            }
        }
    }
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        AMREX_ALWAYS_ASSERT(c[idim]->nComp() == 1 || c[idim]->nComp() == ncomp);
        for (int icomp = 0; icomp < ncomp; ++icomp) {
            int const scomp = (c[idim]->nComp() == 1) ? 0 : icomp;
            m_c_coeffs[amrlev][0][idim].LocalCopy(*c[idim], scomp, icomp, 1, IntVect(0));
        }
#if (AMREX_SPACEDIM != 3)
        this->applyMetricTerm(amrlev, 0, m_c_coeffs[amrlev][0][idim]);
#endif
    }
    m_drift_needs_update = true;
}

template <typename MF>
void
MLABecLapDriftT<MF>::update ()
{
    MLABecLaplacianT<MF>::update();
    updateDrift();
}

template <typename MF>
void
MLABecLapDriftT<MF>::prepareForSolve ()
{
    BL_PROFILE("MLABecLapDrift::prepareForSolve()");
    MLABecLaplacianT<MF>::prepareForSolve();
    updateDrift();
}

template <typename MF>
void
MLABecLapDriftT<MF>::updateDrift ()
{
    m_drift_needs_update = false;
    if (m_c_scalar == RT(0.0)) { return; }
    checkDrift();

    // Average down as MLABecLaplacian does b: the MG levels of each AMR level,
    // then from the coarsest MG level of an AMR level to the next AMR level
    for (int alev = this->m_num_amr_levels-1; alev >= 0; --alev) {
        auto& c = m_c_coeffs[alev];
        for (int mglev = 1; mglev < this->m_num_mg_levels[alev]; ++mglev) {
            IntVect const ratio = (alev > 0) ? IntVect(this->mg_coarsen_ratio)
                                             : this->mg_coarsen_ratio_vec[mglev-1];
            amrex::average_down_faces(amrex::GetArrOfConstPtrs(c[mglev-1]),
                                      amrex::GetArrOfPtrs(c[mglev]), ratio, 0);
        }
        if (alev > 0) {
            amrex::average_down_faces(amrex::GetArrOfConstPtrs(c.back()),
                                      amrex::GetArrOfPtrs(m_c_coeffs[alev-1].front()),
                                      IntVect(this->mg_coarsen_ratio),
                                      this->m_geom[alev-1][0]);
        }
    }

    // The drift fluxes cancel in the sum over the domain, so the constants
    // stay in the left null space (the one the right-hand side must be
    // orthogonal to), unless the drift crosses the boundary of the level.
    // boundaryDrift also checks the Neumann faces.
    bool domain_neumann = false;
    for (int n = 0; n < this->getNComp(); ++n) {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            domain_neumann = domain_neumann
                || this->m_lobc[n][idim] == LinOpBCType::Neumann
                || this->m_hibc[n][idim] == LinOpBCType::Neumann;
        }
    }
    for (int alev = 0; alev < this->m_num_amr_levels; ++alev) {
        if (!this->m_is_singular[alev] && !domain_neumann && !cfNeumann(alev)) { continue; }
        int const ncross = boundaryDrift(alev);
        if (this->m_is_singular[alev]) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ncross == 0 || ncross == this->getNComp(),
                "MLABecLapDrift: the drift of some components crosses the boundary and "
                "of others not, so the operator is singular for some components only");
            this->m_is_singular[alev] = (ncross == 0);
        }
    }
}

template <typename MF>
void
MLABecLapDriftT<MF>::checkDrift () const
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!m_c_coeffs.empty(),
        "MLABecLapDrift: setScalars(A, B, C) with C != 0 needs setDriftCoeffs");
    bool ok = !this->hasRobinBC() && !this->hasInhomogNeumannBC() && !this->hasHiddenDimension();
    for (int mglev = 1; mglev < this->m_num_mg_levels[0]; ++mglev) {
        ok = ok && (this->mg_coarsen_ratio_vec[mglev-1] == this->mg_coarsen_ratio);
    }
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ok,
        "MLABecLapDrift with drift: no semi-coarsening, hidden dimension, "
        "Robin or inhomogeneous Neumann boundary conditions");
}

template <typename MF>
MLABecLapDriftCoeffsMA<typename MF::value_type>
MLABecLapDriftT<MF>::driftArrays (int amrlev, int mglev) const
{
    MLABecLapDriftCoeffsMA<RT> r;
    r.a = this->m_a_coeffs[amrlev][mglev].const_arrays();
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        r.b[idim] = this->m_b_coeffs[amrlev][mglev][idim].const_arrays();
        r.c[idim] = m_c_coeffs[amrlev][mglev][idim].const_arrays();
        r.dxinv[idim] = static_cast<RT>(this->m_geom[amrlev][mglev].InvCellSize(idim));
    }
    r.alpha = this->m_a_scalar;
    r.beta = this->m_b_scalar;
    r.cscalar = m_c_scalar;
    return r;
}

template <typename MF>
void
MLABecLapDriftT<MF>::FFlux (int amrlev, const MFIter& mfi,
                            const Array<FAB*,AMREX_SPACEDIM>& flux,
                            const FAB& sol, Location loc, int face_only) const
{
    MLABecLaplacianT<MF>::FFlux(amrlev, mfi, flux, sol, loc, face_only);
    if (!hasDrift()) { return; }

    // Add the drift flux -C c phi_upwind on the same faces
    const int mglev = 0;
    const Box& box = mfi.tilebox();
    const int ncomp = this->getNComp();
    const RT cscalar = m_c_scalar;
    const auto& solarr = sol.const_array();
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        const auto& c = m_c_coeffs[amrlev][mglev][idim].const_array(mfi);
        const auto& f = flux[idim]->array();
        // face_only: the low faces, each with the high face blen cells up
        Box const fbx = face_only ? amrex::bdryLo(box, idim)
                                  : amrex::surroundingNodes(box, idim);
        IntVect const hi_shift = IntVect::TheDimensionVector(idim) * box.length(idim);
        ParallelFor(fbx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            amrex::ignore_unused(k);
            IntVect const iv(AMREX_D_DECL(i,j,k));
            mlabeclap_flux_drift(iv, idim, n, f, solarr, c, cscalar);
            if (face_only) {
                mlabeclap_flux_drift(iv+hi_shift, idim, n, f, solarr, c, cscalar);
            }
        });
    }
}

template <typename MF>
void
MLABecLapDriftT<MF>::Fapply (int amrlev, int mglev, MF& out, const MF& in) const
{
    if (!hasDrift()) {
        MLABecLaplacianT<MF>::Fapply(amrlev, mglev, out, in);
        return;
    }
    BL_PROFILE("MLABecLapDrift::Fapply()");
    const auto op = driftArrays(amrlev, mglev);
    const auto& yma = out.arrays();
    const auto& xma = in.const_arrays();
    ParallelFor(out, IntVect(0), this->getNComp(),
    [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
    {
        mlabeclap_adotx_drift(i, j, k, n, yma[box_no], xma[box_no], op[box_no]);
    });
    if (!Gpu::inNoSyncRegion()) { Gpu::streamSynchronize(); }
}

template <typename MF>
void
MLABecLapDriftT<MF>::Fsmooth (int amrlev, int mglev, MF& sol, const MF& rhs,
                              int redblack) const
{
    if (!hasDrift()) {
        MLABecLaplacianT<MF>::Fsmooth(amrlev, mglev, sol, rhs, redblack);
        return;
    }
    BL_PROFILE("MLABecLapDrift::Fsmooth()");
    MF Ax;
    MultiArray4<RT const> axma; // Jacobi only
    if (!this->m_use_gauss_seidel) {
        Ax.define(sol.boxArray(), sol.DistributionMap(), sol.nComp(), 0,
                  MFInfo().SetArena(The_Async_Arena()));
        Fapply(amrlev, mglev, Ax, sol);
        axma = Ax.const_arrays();
    }

    const auto& undrrelxr = this->m_undrrelxr[amrlev][mglev];
    const auto& maskvals  = this->m_maskvals [amrlev][mglev];
    GpuArray<MultiArray4<int const>,2*AMREX_SPACEDIM> mma;
    GpuArray<MultiArray4<RT const>,2*AMREX_SPACEDIM> fma;
    for (OrientationIter oit; oit; ++oit) {
        Orientation const ori = oit();
        int const cdir(ori);
        mma[cdir] = maskvals[cdir].const_arrays();
        fma[cdir] = undrrelxr[ori].const_arrays();
    }

    const auto op = driftArrays(amrlev, mglev);
    const auto& solma = sol.arrays();
    const auto& rhsma = rhs.const_arrays();
    // Jacobi takes 2/3 of the update, as without drift
    const RT omega = this->m_use_gauss_seidel ? m_drift_omega : RT(2.0/3.0)*m_drift_omega;
    auto const relax = [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
    {
        GpuArray<Array4<int const>,2*AMREX_SPACEDIM> m;
        GpuArray<Array4<RT const>,2*AMREX_SPACEDIM> f;
        for (int cdir = 0; cdir < 2*AMREX_SPACEDIM; ++cdir) {
            m[cdir] = mma[cdir][box_no];
            f[cdir] = fma[cdir][box_no];
        }
        mlabeclap_relax_drift(i, j, k, n, solma[box_no], rhsma[box_no],
                              axma ? axma[box_no] : Array4<RT const>{}, op[box_no],
                              m, f, omega, Box(op.a[box_no]));
    };
    const int ncomp = this->getNComp();
    if (this->m_use_gauss_seidel) {
        ParallelForRedBlack(sol, ncomp, redblack, relax);
    } else {
        ParallelFor(sol, IntVect(0), ncomp, relax);
    }
    if (!Gpu::inNoSyncRegion()) { Gpu::streamSynchronize(); }
}

template <typename MF>
void
MLABecLapDriftT<MF>::normalize (int amrlev, int mglev, MF& mf) const
{
    if (!hasDrift()) {
        MLABecLaplacianT<MF>::normalize(amrlev, mglev, mf);
        return;
    }
    BL_PROFILE("MLABecLapDrift::normalize()");
    const auto op = driftArrays(amrlev, mglev);
    const auto& xma = mf.arrays();
    ParallelFor(mf, IntVect(0), this->getNComp(),
    [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
    {
        mlabeclap_normalize_drift(i, j, k, n, xma[box_no], op[box_no]);
    });
    if (!Gpu::inNoSyncRegion()) { Gpu::streamSynchronize(); }
}

template <typename MF>
int
MLABecLapDriftT<MF>::boundaryDrift (int amrlev) const
{
    const int ncomp = this->getNComp();
    const auto& maskvals = this->m_maskvals[amrlev][0];
    bool const cf_neumann = cfNeumann(amrlev);
    // Per component: the max over the boundary faces of |c| and of the inflow
    // through the Neumann ones, and the max |c| over all the faces
    Vector<RT> bmax(3*ncomp, RT(0.0));
    for (int n = 0; n < ncomp; ++n) {
        ReduceOps<ReduceOpMax, ReduceOpMax, ReduceOpMax> reduce_op;
        ReduceData<RT, RT, RT> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;
        for (MFIter mfi(this->m_a_coeffs[amrlev][0]); mfi.isValid(); ++mfi) {
            Box const& vbox = mfi.validbox();
            for (OrientationIter oit; oit; ++oit) {
                Orientation const ori = oit();
                int const idim = ori.coordDir();
                auto const& c = m_c_coeffs[amrlev][0][idim].const_array(mfi);
                auto const& mask = maskvals[ori].const_array(mfi);
                // The outside cell of a face, and the sign of the inflow
                IntVect const out = ori.isLow() ? -IntVect::TheDimensionVector(idim) : IntVect(0);
                RT const sign = ori.isLow() ? RT(1.0) : RT(-1.0);
                bool const domain_neumann = (ori.isLow() ? this->m_lobc : this->m_hibc)[n][idim]
                    == LinOpBCType::Neumann;
                reduce_op.eval(ori.isLow() ? amrex::bdryLo(vbox, idim) : amrex::bdryHi(vbox, idim),
                               reduce_data,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
                {
                    amrex::ignore_unused(k);
                    IntVect const iv(AMREX_D_DECL(i,j,k));
                    int const m = mask(iv+out);
                    if (m == BndryData::covered) { return {RT(0.0), RT(0.0), RT(0.0)}; }
                    bool const neumann = (m == BndryData::outside_domain) ? domain_neumann
                                                                         : cf_neumann;
                    return {std::abs(c(iv,n)), neumann ? sign*c(iv,n) : RT(0.0), RT(0.0)};
                });
                if (ori.isLow()) {
                    reduce_op.eval(amrex::surroundingNodes(vbox, idim), reduce_data,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
                    {
                        return {RT(0.0), RT(0.0), std::abs(c(i,j,k,n))};
                    });
                }
            }
        }
        auto const& v = reduce_data.value(reduce_op);
        bmax[3*n  ] = amrex::get<0>(v);
        bmax[3*n+1] = amrex::get<1>(v);
        bmax[3*n+2] = amrex::get<2>(v);
    }
    ParallelAllReduce::Max(bmax.data(), 3*ncomp, ParallelContext::CommunicatorSub());

    // Round-off velocities on the boundary do not count
    int ncross = 0;
    for (int n = 0; n < ncomp; ++n) {
        RT const tol = bmax[3*n+2] * RT(1.e-12);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(bmax[3*n+1] <= tol,
            "MLABecLapDrift: the drift enters through a Neumann face, where the "
            "ghost cell is the cell itself; use a Dirichlet face or c = 0 there");
        if (bmax[3*n] > tol) { ++ncross; }
    }
    return ncross;
}

using MLABecLapDrift = MLABecLapDriftT<MultiFab>;

}

#endif
