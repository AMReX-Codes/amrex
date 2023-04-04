#ifndef AMREX_ML_ABECLAPLACIAN_H_
#define AMREX_ML_ABECLAPLACIAN_H_
#include <AMReX_Config.H>

#include <AMReX_MLCellABecLap.H>
#include <AMReX_MLABecLap_K.H>

namespace amrex {

// (alpha * a - beta * (del dot b grad)) phi

template <typename MF>
class MLABecLaplacianT
    : public MLCellABecLapT<MF>
{
public:

    using FAB = typename MF::fab_type;
    using RT  = typename MF::value_type;

    using BCType = LinOpBCType;
    using Location  = typename MLLinOpT<MF>::Location;

    MLABecLaplacianT () {}
    MLABecLaplacianT (const Vector<Geometry>& a_geom,
                      const Vector<BoxArray>& a_grids,
                      const Vector<DistributionMapping>& a_dmap,
                      const LPInfo& a_info = LPInfo(),
                      const Vector<FabFactory<FAB> const*>& a_factory = {},
                      const int a_ncomp = 1);

    MLABecLaplacianT (const Vector<Geometry>& a_geom,
                      const Vector<BoxArray>& a_grids,
                      const Vector<DistributionMapping>& a_dmap,
                      const Vector<iMultiFab const*>& a_overset_mask, // 1: unknown, 0: known
                      const LPInfo& a_info = LPInfo(),
                      const Vector<FabFactory<FAB> const*>& a_factory = {},
                      const int a_ncomp = 1);

    virtual ~MLABecLaplacianT ();

    MLABecLaplacianT (const MLABecLaplacianT<MF>&) = delete;
    MLABecLaplacianT (MLABecLaplacianT<MF>&&) = delete;
    MLABecLaplacianT<MF>& operator= (const MLABecLaplacianT<MF>&) = delete;
    MLABecLaplacianT<MF>& operator= (MLABecLaplacianT<MF>&&) = delete;

    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info = LPInfo(),
                 const Vector<FabFactory<FAB> const*>& a_factory = {},
                 const int a_ncomp = 1);

    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const Vector<iMultiFab const*>& a_overset_mask,
                 const LPInfo& a_info = LPInfo(),
                 const Vector<FabFactory<FAB> const*>& a_factory = {},
                 const int a_ncomp = 1);

    /**
     * Set scalar constants A and B in the equation:
     * (A \alpha - B \nabla \cdot \beta \nabla ) \phi = f
     * for the Multi-Level AB Laplacian Solver.
     */
    template <typename T1, typename T2,
              std::enable_if_t<std::is_convertible_v<T1,typename MF::value_type> &&
                               std::is_convertible_v<T2,typename MF::value_type>,
                               int> = 0>
    void setScalars (T1 a, T2 b) noexcept;

    /**
     * Sets alpha as a scalar field to values from a single component
     * mutlifab.
     *
     * \param [in] amrlev The level of the multifab for the solver, with
     *                    \p amrlev = 0 always being the lowest level in the
     *                    AMR hierarchy represented in the solve.
     * \param [in] alpha  Multifab of alpha values.
     */
    template <typename AMF,
              std::enable_if_t<IsFabArray<AMF>::value &&
                               std::is_convertible_v<typename AMF::value_type,
                                                     typename  MF::value_type>,
                               int> = 0>
    void setACoeffs (int amrlev, const AMF& alpha);

    /**
     * Sets alpha as a single scalar constant value across
     * the mutlifab.
     *
     * \param [in] amrlev The level of the multifab for the solver, with
     *                    \p amrlev = 0 always being the lowest level in the
     *                    AMR hierarchy represented in the solve.
     * \param [in] alpha  Single scalar value to populate across mutlifab.
     */
    template <typename T,
              std::enable_if_t<std::is_convertible_v<T,typename MF::value_type>,
                               int> = 0>
    void setACoeffs (int amrlev, T alpha);

    /**
     * Sets beta as a scalar field to be the values defined
     * in the supplied multifabs (one for each space dimension).
     *
     * \param [in] amrlev The level of the multifab for the solver, with
     *                    \p amrlev = 0 always being the lowest level in the
     *                    AMR hierarchy represented in the solve.
     * \param [in] beta   Array of Multifabs of beta values.
     */
    template <typename AMF,
              std::enable_if_t<IsFabArray<AMF>::value &&
                               std::is_convertible_v<typename AMF::value_type,
                                                     typename  MF::value_type>,
                               int> = 0>
    void setBCoeffs (int amrlev, const Array<AMF const*,AMREX_SPACEDIM>& beta);

    /**
     * Sets beta as a single scalar constant value across
     * the mutlifabs (one for each dimension).
     *
     * \param [in] amrlev The level of the multifab for the solver, with
     *                    \p amrlev = 0 always being the lowest level in the
     *                    AMR hierarchy represented in the solve.
     * \param [in] beta   Single scalar value to populate across mutlifabs.
     */
    template <typename T,
              std::enable_if_t<std::is_convertible_v<T,typename MF::value_type>,
                               int> = 0>
    void setBCoeffs (int amrlev, T beta);

    /**
     * Set each beta component to a single scalar constant value corresponding to the
     * respective component of the supplied vector.
     *
     * \param [in] amrlev The level of the multifab for the solver, with
     *                    \p amrlev = 0 always being the lowest level in the
     *                    AMR hierarchy represented in the solve.
     * \param [in] beta   Vector of scalar constant values.
     */
    template <typename T,
              std::enable_if_t<std::is_convertible_v<T,typename MF::value_type>,
                               int> = 0>
    void setBCoeffs (int amrlev, Vector<T> const& beta);

    virtual int getNComp () const override { return m_ncomp; }

    virtual bool needsUpdate () const override {
        return (m_needs_update || MLCellABecLapT<MF>::needsUpdate());
    }
    virtual void update () override;

    virtual void prepareForSolve () override;
    virtual bool isSingular (int amrlev) const override { return m_is_singular[amrlev]; }
    virtual bool isBottomSingular () const override { return m_is_singular[0]; }
    virtual void Fapply (int amrlev, int mglev, MF& out, const MF& in) const final override;
    virtual void Fsmooth (int amrlev, int mglev, MF& sol, const MF& rhs, int redblack) const final override;
    virtual void FFlux (int amrlev, const MFIter& mfi,
                        const Array<FAB*,AMREX_SPACEDIM>& flux,
                        const FAB& sol, Location /* loc */,
                        const int face_only=0) const final override;

    virtual void normalize (int amrlev, int mglev, MF& mf) const final override;

    virtual RT getAScalar () const final override { return m_a_scalar; }
    virtual RT getBScalar () const final override { return m_b_scalar; }
    virtual MF const* getACoeffs (int amrlev, int mglev) const final override
        { return &(m_a_coeffs[amrlev][mglev]); }
    virtual Array<MF const*,AMREX_SPACEDIM> getBCoeffs (int amrlev, int mglev) const final override
        { return amrex::GetArrOfConstPtrs(m_b_coeffs[amrlev][mglev]); }

    virtual std::unique_ptr<MLLinOpT<MF>> makeNLinOp (int /*grid_size*/) const final override;

    virtual bool supportNSolve () const final override;

    virtual void copyNSolveSolution (MF& dst, MF const& src) const final override;

    void averageDownCoeffsSameAmrLevel (int amrlev, Vector<MF>& a,
                                        Vector<Array<MF,AMREX_SPACEDIM> >& b);
    void averageDownCoeffs ();
    void averageDownCoeffsToCoarseAmrLevel (int flev);

    void applyMetricTermsCoeffs ();

    void applyRobinBCTermsCoeffs ();

    static void FFlux (Box const& box, Real const* dxinv, RT bscalar,
                       Array<FAB const*, AMREX_SPACEDIM> const& bcoef,
                       Array<FAB*,AMREX_SPACEDIM> const& flux,
                       FAB const& sol, int face_only, int ncomp);

protected:

    bool m_needs_update = true;

    RT m_a_scalar = std::numeric_limits<RT>::quiet_NaN();
    RT m_b_scalar = std::numeric_limits<RT>::quiet_NaN();
    Vector<Vector<MF> > m_a_coeffs;
    Vector<Vector<Array<MF,AMREX_SPACEDIM> > > m_b_coeffs;

    Vector<int> m_is_singular;

    virtual bool supportRobinBC () const noexcept override { return true; }

private:

    int m_ncomp = 1;

    void define_ab_coeffs ();

    void update_singular_flags ();
};

template <typename MF>
MLABecLaplacianT<MF>::MLABecLaplacianT (const Vector<Geometry>& a_geom,
                                        const Vector<BoxArray>& a_grids,
                                        const Vector<DistributionMapping>& a_dmap,
                                        const LPInfo& a_info,
                                        const Vector<FabFactory<FAB> const*>& a_factory,
                                        const int a_ncomp)
{
    define(a_geom, a_grids, a_dmap, a_info, a_factory, a_ncomp);
}

template <typename MF>
MLABecLaplacianT<MF>::MLABecLaplacianT (const Vector<Geometry>& a_geom,
                                        const Vector<BoxArray>& a_grids,
                                        const Vector<DistributionMapping>& a_dmap,
                                        const Vector<iMultiFab const*>& a_overset_mask,
                                        const LPInfo& a_info,
                                        const Vector<FabFactory<FAB> const*>& a_factory,
                                        const int a_ncomp)
{
    define(a_geom, a_grids, a_dmap, a_overset_mask, a_info, a_factory, a_ncomp);
}

template <typename MF> MLABecLaplacianT<MF>::~MLABecLaplacianT () {}

template <typename MF>
void
MLABecLaplacianT<MF>::define (const Vector<Geometry>& a_geom,
                              const Vector<BoxArray>& a_grids,
                              const Vector<DistributionMapping>& a_dmap,
                              const LPInfo& a_info,
                              const Vector<FabFactory<FAB> const*>& a_factory,
                              const int a_ncomp)
{
    BL_PROFILE("MLABecLaplacian::define()");
    this->m_ncomp = a_ncomp;
    MLCellABecLapT<MF>::define(a_geom, a_grids, a_dmap, a_info, a_factory);
    define_ab_coeffs();
}

template <typename MF>
void
MLABecLaplacianT<MF>::define (const Vector<Geometry>& a_geom,
                              const Vector<BoxArray>& a_grids,
                              const Vector<DistributionMapping>& a_dmap,
                              const Vector<iMultiFab const*>& a_overset_mask,
                              const LPInfo& a_info,
                              const Vector<FabFactory<FAB> const*>& a_factory,
                              const int a_ncomp)
{
    BL_PROFILE("MLABecLaplacian::define(overset)");
    this->m_ncomp = a_ncomp;
    MLCellABecLapT<MF>::define(a_geom, a_grids, a_dmap, a_overset_mask, a_info, a_factory);
    define_ab_coeffs();
}

template <typename MF>
void
MLABecLaplacianT<MF>::define_ab_coeffs ()
{
    m_a_coeffs.resize(this->m_num_amr_levels);
    m_b_coeffs.resize(this->m_num_amr_levels);
    for (int amrlev = 0; amrlev < this->m_num_amr_levels; ++amrlev)
    {
        m_a_coeffs[amrlev].resize(this->m_num_mg_levels[amrlev]);
        m_b_coeffs[amrlev].resize(this->m_num_mg_levels[amrlev]);
        for (int mglev = 0; mglev < this->m_num_mg_levels[amrlev]; ++mglev)
        {
            m_a_coeffs[amrlev][mglev].define
                (this->m_grids[amrlev][mglev], this->m_dmap[amrlev][mglev],
                 1, 0, MFInfo(), *(this->m_factory[amrlev][mglev]));
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
            {
                const BoxArray& ba = amrex::convert(this->m_grids[amrlev][mglev],
                                                    IntVect::TheDimensionVector(idim));
                m_b_coeffs[amrlev][mglev][idim].define
                    (ba, this->m_dmap[amrlev][mglev], m_ncomp, 0, MFInfo(),
                     *(this->m_factory[amrlev][mglev]));
            }
        }
    }
}

template <typename MF>
template <typename T1, typename T2,
          std::enable_if_t<std::is_convertible_v<T1,typename MF::value_type> &&
                           std::is_convertible_v<T2,typename MF::value_type>, int>>
void
MLABecLaplacianT<MF>::setScalars (T1 a, T2 b) noexcept
{
    m_a_scalar = RT(a);
    m_b_scalar = RT(b);
    if (m_a_scalar == RT(0.0)) {
        for (int amrlev = 0; amrlev < this->m_num_amr_levels; ++amrlev) {
            m_a_coeffs[amrlev][0].setVal(RT(0.0));
        }
    }
}

template <typename MF>
template <typename AMF,
          std::enable_if_t<IsFabArray<AMF>::value &&
                           std::is_convertible_v<typename AMF::value_type,
                                                 typename  MF::value_type>, int>>
void
MLABecLaplacianT<MF>::setACoeffs (int amrlev, const AMF& alpha)
{
    AMREX_ASSERT_WITH_MESSAGE(alpha.nComp() == 1,
                              "MLABecLaplacian::setACoeffs: alpha is supposed to be single component.");
    m_a_coeffs[amrlev][0].LocalCopy(alpha, 0, 0, 1, IntVect(0));
    m_needs_update = true;
}

template <typename MF>
template <typename T,
          std::enable_if_t<std::is_convertible_v<T,typename MF::value_type>,int>>
void
MLABecLaplacianT<MF>::setACoeffs (int amrlev, T alpha)
{
    m_a_coeffs[amrlev][0].setVal(RT(alpha));
    m_needs_update = true;
}


template <typename MF>
template <typename AMF,
          std::enable_if_t<IsFabArray<AMF>::value &&
                           std::is_convertible_v<typename AMF::value_type,
                                                 typename  MF::value_type>, int>>
void
MLABecLaplacianT<MF>::setBCoeffs (int amrlev,
                                  const Array<AMF const*,AMREX_SPACEDIM>& beta)
{
    const int ncomp = this->getNComp();
    AMREX_ASSERT(beta[0]->nComp() == 1 || beta[0]->nComp() == ncomp);
    if (beta[0]->nComp() == ncomp) {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            for (int icomp = 0; icomp < ncomp; ++icomp) {
                m_b_coeffs[amrlev][0][idim].LocalCopy(*beta[idim], icomp, icomp, 1, IntVect(0));
            }
        }
    } else {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            for (int icomp = 0; icomp < ncomp; ++icomp) {
                m_b_coeffs[amrlev][0][idim].LocalCopy(*beta[idim], 0, icomp, 1, IntVect(0));
            }
        }
    }
    m_needs_update = true;
}

template <typename MF>
template <typename T,
          std::enable_if_t<std::is_convertible_v<T,typename MF::value_type>,int>>
void
MLABecLaplacianT<MF>::setBCoeffs (int amrlev, T beta)
{
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        m_b_coeffs[amrlev][0][idim].setVal(RT(beta));
    }
    m_needs_update = true;
}

template <typename MF>
template <typename T,
          std::enable_if_t<std::is_convertible_v<T,typename MF::value_type>,int>>
void
MLABecLaplacianT<MF>::setBCoeffs (int amrlev, Vector<T> const& beta)
{
    const int ncomp = this->getNComp();
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        for (int icomp = 0; icomp < ncomp; ++icomp) {
            m_b_coeffs[amrlev][0][idim].setVal(RT(beta[icomp]));
        }
    }
    m_needs_update = true;
}

template <typename MF>
void
MLABecLaplacianT<MF>::update ()
{
    if (MLCellABecLapT<MF>::needsUpdate()) {
        MLCellABecLapT<MF>::update();
    }

#if (AMREX_SPACEDIM != 3)
    applyMetricTermsCoeffs();
#endif

    averageDownCoeffs();

    update_singular_flags();

    m_needs_update = false;
}

template <typename MF>
void
MLABecLaplacianT<MF>::prepareForSolve ()
{
    BL_PROFILE("MLABecLaplacian::prepareForSolve()");

    MLCellABecLapT<MF>::prepareForSolve();

#if (AMREX_SPACEDIM != 3)
    applyMetricTermsCoeffs();
#endif

    applyRobinBCTermsCoeffs();

    averageDownCoeffs();

    update_singular_flags();

    m_needs_update = false;
}

template <typename MF>
void
MLABecLaplacianT<MF>::applyMetricTermsCoeffs ()
{
#if (AMREX_SPACEDIM != 3)
    for (int alev = 0; alev < this->m_num_amr_levels; ++alev)
    {
        const int mglev = 0;
        this->applyMetricTerm(alev, mglev, m_a_coeffs[alev][mglev]);
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            this->applyMetricTerm(alev, mglev, m_b_coeffs[alev][mglev][idim]);
        }
    }
#endif
}

//
// Suppose we are solving `alpha u - del (beta grad u) = rhs` (Scalar
// coefficients can be easily added back in the end) and there is Robin BC
// `a u + b du/dn = f` at the upper end of the x-direction.  The 1D
// discretization at the last cell i is
//
//    alpha u_i + (beta_{i-1/2} (du/dx)_{i-1/2} - beta_{i+1/2} (du/dx)_{i+1/2}) / h = rhs_i
//
// where h is the cell size.  At `i+1/2` (i.e., the boundary), we have
//
//    a (u_i + u_{i+1})/2 + b (u_{i+1}-u_i)/h = f,
//
// according to the Robin BC.  This gives
//
//    u_{i+1} = A + B u_i,
//
// where `A = f/(b/h + a/2)` and `B = (b/h - a/2) / (b/h + a/2).  We then
// use `u_i` and `u_{i+1}` to compute `(du/dx)_{i+1/2}`.  The discretization
// at cell i then becomes
//
//    \tilde{alpha}_i u_i + (beta_{i-1/2} (du/dx)_{i-1/2} - 0) / h = \tilde{rhs}_i
//
// This is equivalent to having homogeneous Neumann BC with modified alpha and rhs.
//
//    \tilde{alpha}_i = alpha_i + (1-B) beta_{i+1/2} / h^2
//    \tilde{rhs}_i = rhs_i + A beta_{i+1/2} / h^2
//
template <typename MF>
void
MLABecLaplacianT<MF>::applyRobinBCTermsCoeffs ()
{
    if (!(this->hasRobinBC())) return;

    const int ncomp = this->getNComp();
    bool reset_alpha = false;
    if (m_a_scalar == RT(0.0)) {
        m_a_scalar = RT(1.0);
        reset_alpha = true;
    }
    const RT bovera = m_b_scalar/m_a_scalar;

    for (int amrlev = 0; amrlev < this->m_num_amr_levels; ++amrlev) {
        const int mglev = 0;
        const Box& domain = this->m_geom[amrlev][mglev].Domain();
        const RT dxi = static_cast<RT>(this->m_geom[amrlev][mglev].InvCellSize(0));
        const RT dyi = static_cast<RT>((AMREX_SPACEDIM >= 2) ? this->m_geom[amrlev][mglev].InvCellSize(1) : Real(1.0));
        const RT dzi = static_cast<RT>((AMREX_SPACEDIM == 3) ? this->m_geom[amrlev][mglev].InvCellSize(2) : Real(1.0));

        if (reset_alpha) {
            m_a_coeffs[amrlev][mglev].setVal(RT(0.0));
        }

        MFItInfo mfi_info;
        if (Gpu::notInLaunchRegion()) mfi_info.SetDynamic(true);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(m_a_coeffs[amrlev][mglev], mfi_info); mfi.isValid(); ++mfi)
        {
            const Box& vbx = mfi.validbox();
            auto const& afab = m_a_coeffs[amrlev][mglev].array(mfi);
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                auto const& bfab = m_b_coeffs[amrlev][mglev][idim].const_array(mfi);
                const Box& blo = amrex::adjCellLo(vbx,idim);
                const Box& bhi = amrex::adjCellHi(vbx,idim);
                bool outside_domain_lo = !(domain.contains(blo));
                bool outside_domain_hi = !(domain.contains(bhi));
                if ((!outside_domain_lo) && (!outside_domain_hi)) continue;
                for (int icomp = 0; icomp < ncomp; ++icomp) {
                    auto const& rbc = (*(this->m_robin_bcval[amrlev]))[mfi].const_array(icomp*3);
                    if (this->m_lobc_orig[icomp][idim] == LinOpBCType::Robin && outside_domain_lo)
                    {
                        if (idim == 0) {
                            RT fac = bovera*dxi*dxi;
                            AMREX_HOST_DEVICE_FOR_3D(blo, i, j, k,
                            {
                                RT B = (rbc(i,j,k,1)*dxi - rbc(i,j,k,0)*RT(0.5))
                                  /    (rbc(i,j,k,1)*dxi + rbc(i,j,k,0)*RT(0.5));
                                afab(i+1,j,k,icomp) += fac*bfab(i+1,j,k,icomp)*(RT(1.0)-B);
                            });
                        } else if (idim == 1) {
                            RT fac = bovera*dyi*dyi;
                            AMREX_HOST_DEVICE_FOR_3D(blo, i, j, k,
                            {
                                RT B = (rbc(i,j,k,1)*dyi - rbc(i,j,k,0)*RT(0.5))
                                  /    (rbc(i,j,k,1)*dyi + rbc(i,j,k,0)*RT(0.5));
                                afab(i,j+1,k,icomp) += fac*bfab(i,j+1,k,icomp)*(RT(1.0)-B);
                            });
                        } else {
                            RT fac = bovera*dzi*dzi;
                            AMREX_HOST_DEVICE_FOR_3D(blo, i, j, k,
                            {
                                RT B = (rbc(i,j,k,1)*dzi - rbc(i,j,k,0)*RT(0.5))
                                  /    (rbc(i,j,k,1)*dzi + rbc(i,j,k,0)*RT(0.5));
                                afab(i,j,k+1,icomp) += fac*bfab(i,j,k+1,icomp)*(RT(1.0)-B);
                            });
                        }
                    }
                    if (this->m_hibc_orig[icomp][idim] == LinOpBCType::Robin && outside_domain_hi)
                    {
                        if (idim == 0) {
                            RT fac = bovera*dxi*dxi;
                            AMREX_HOST_DEVICE_FOR_3D(bhi, i, j, k,
                            {
                                RT B = (rbc(i,j,k,1)*dxi - rbc(i,j,k,0)*RT(0.5))
                                  /    (rbc(i,j,k,1)*dxi + rbc(i,j,k,0)*RT(0.5));
                                afab(i-1,j,k,icomp) += fac*bfab(i,j,k,icomp)*(RT(1.0)-B);
                            });
                        } else if (idim == 1) {
                            RT fac = bovera*dyi*dyi;
                            AMREX_HOST_DEVICE_FOR_3D(bhi, i, j, k,
                            {
                                RT B = (rbc(i,j,k,1)*dyi - rbc(i,j,k,0)*RT(0.5))
                                  /    (rbc(i,j,k,1)*dyi + rbc(i,j,k,0)*RT(0.5));
                                afab(i,j-1,k,icomp) += fac*bfab(i,j,k,icomp)*(RT(1.0)-B);
                            });
                        } else {
                            RT fac = bovera*dzi*dzi;
                            AMREX_HOST_DEVICE_FOR_3D(bhi, i, j, k,
                            {
                                RT B = (rbc(i,j,k,1)*dzi - rbc(i,j,k,0)*RT(0.5))
                                  /    (rbc(i,j,k,1)*dzi + rbc(i,j,k,0)*RT(0.5));
                                afab(i,j,k-1,icomp) += fac*bfab(i,j,k,icomp)*(RT(1.0)-B);
                            });
                        }
                    }
                }
            }
        }
    }
}

template <typename MF>
void
MLABecLaplacianT<MF>::averageDownCoeffs ()
{
    BL_PROFILE("MLABecLaplacian::averageDownCoeffs()");

    for (int amrlev = this->m_num_amr_levels-1; amrlev > 0; --amrlev)
    {
        auto& fine_a_coeffs = m_a_coeffs[amrlev];
        auto& fine_b_coeffs = m_b_coeffs[amrlev];

        averageDownCoeffsSameAmrLevel(amrlev, fine_a_coeffs, fine_b_coeffs);
        averageDownCoeffsToCoarseAmrLevel(amrlev);
    }

    averageDownCoeffsSameAmrLevel(0, m_a_coeffs[0], m_b_coeffs[0]);
}

template <typename MF>
void
MLABecLaplacianT<MF>::averageDownCoeffsSameAmrLevel (int amrlev, Vector<MF>& a,
                                                     Vector<Array<MF,AMREX_SPACEDIM> >& b)
{
    int nmglevs = a.size();
    for (int mglev = 1; mglev < nmglevs; ++mglev)
    {
        IntVect ratio = (amrlev > 0) ? IntVect(this->mg_coarsen_ratio) : this->mg_coarsen_ratio_vec[mglev-1];

        if (m_a_scalar == 0.0)
        {
            a[mglev].setVal(RT(0.0));
        }
        else
        {
            amrex::average_down(a[mglev-1], a[mglev], 0, 1, ratio);
        }

        Vector<const MF*> fine {AMREX_D_DECL(&(b[mglev-1][0]),
                                             &(b[mglev-1][1]),
                                             &(b[mglev-1][2]))};
        Vector<MF*> crse {AMREX_D_DECL(&(b[mglev][0]),
                                       &(b[mglev][1]),
                                       &(b[mglev][2]))};

        amrex::average_down_faces(fine, crse, ratio, 0);
    }

    for (int mglev = 1; mglev < nmglevs; ++mglev)
    {
        if (this->m_overset_mask[amrlev][mglev]) {
            const RT fac = static_cast<RT>(1 << mglev); // 2**mglev
            const RT osfac = RT(2.0)*fac/(fac+RT(1.0));
            const int ncomp = this->getNComp();
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(a[mglev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                AMREX_D_TERM(Box const& xbx = mfi.nodaltilebox(0);,
                             Box const& ybx = mfi.nodaltilebox(1);,
                             Box const& zbx = mfi.nodaltilebox(2));
                AMREX_D_TERM(auto const& bx = b[mglev][0].array(mfi);,
                             auto const& by = b[mglev][1].array(mfi);,
                             auto const& bz = b[mglev][2].array(mfi));
                auto const& osm = this->m_overset_mask[amrlev][mglev]->const_array(mfi);
                AMREX_LAUNCH_HOST_DEVICE_LAMBDA_DIM
                    (xbx, t_xbx,
                     {
                         overset_rescale_bcoef_x(t_xbx, bx, osm, ncomp, osfac);
                     },
                     ybx, t_ybx,
                     {
                         overset_rescale_bcoef_y(t_ybx, by, osm, ncomp, osfac);
                     },
                     zbx, t_zbx,
                     {
                         overset_rescale_bcoef_z(t_zbx, bz, osm, ncomp, osfac);
                     });
            }
        }
    }
}

template <typename MF>
void
MLABecLaplacianT<MF>::averageDownCoeffsToCoarseAmrLevel (int flev)
{
    auto& fine_a_coeffs = m_a_coeffs[flev  ].back();
    auto& fine_b_coeffs = m_b_coeffs[flev  ].back();
    auto& crse_a_coeffs = m_a_coeffs[flev-1].front();
    auto& crse_b_coeffs = m_b_coeffs[flev-1].front();

    if (m_a_scalar != 0.0) {
        // We coarsen from the back of flev to the front of flev-1.
        // So we use mg_coarsen_ratio.
        amrex::average_down(fine_a_coeffs, crse_a_coeffs, 0, 1, this->mg_coarsen_ratio);
    }

    amrex::average_down_faces(amrex::GetArrOfConstPtrs(fine_b_coeffs),
                              amrex::GetArrOfPtrs(crse_b_coeffs),
                              IntVect(this->mg_coarsen_ratio),
                              this->m_geom[flev-1][0]);
}

template <typename MF>
void
MLABecLaplacianT<MF>::update_singular_flags ()
{
    m_is_singular.clear();
    m_is_singular.resize(this->m_num_amr_levels, false);
    auto itlo = std::find(this->m_lobc[0].begin(), this->m_lobc[0].end(), BCType::Dirichlet);
    auto ithi = std::find(this->m_hibc[0].begin(), this->m_hibc[0].end(), BCType::Dirichlet);
    if (itlo == this->m_lobc[0].end() && ithi == this->m_hibc[0].end())
    {  // No Dirichlet
        for (int alev = 0; alev < this->m_num_amr_levels; ++alev)
        {
            // For now this assumes that overset regions are treated as Dirichlet bc's
            if (this->m_domain_covered[alev] && !this->m_overset_mask[alev][0])
            {
                if (m_a_scalar == 0.0)
                {
                    m_is_singular[alev] = true;
                }
                else
                {
                    RT asum = m_a_coeffs[alev].back().sum(0,IntVect(0));
                    RT amax = m_a_coeffs[alev].back().norminf(0,1,IntVect(0));
                    m_is_singular[alev] = (std::abs(asum) <= amax * RT(1.e-12));
                }
            }
        }
    }
}

template <typename MF>
void
MLABecLaplacianT<MF>::Fapply (int amrlev, int mglev, MF& out, const MF& in) const
{
    BL_PROFILE("MLABecLaplacian::Fapply()");

    const MF& acoef = m_a_coeffs[amrlev][mglev];
    AMREX_D_TERM(const MF& bxcoef = m_b_coeffs[amrlev][mglev][0];,
                 const MF& bycoef = m_b_coeffs[amrlev][mglev][1];,
                 const MF& bzcoef = m_b_coeffs[amrlev][mglev][2];);

    const GpuArray<RT,AMREX_SPACEDIM> dxinv
        {AMREX_D_DECL(static_cast<RT>(this->m_geom[amrlev][mglev].InvCellSize(0)),
                      static_cast<RT>(this->m_geom[amrlev][mglev].InvCellSize(1)),
                      static_cast<RT>(this->m_geom[amrlev][mglev].InvCellSize(2)))};

    const RT ascalar = m_a_scalar;
    const RT bscalar = m_b_scalar;

    const int ncomp = this->getNComp();

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion() && out.isFusingCandidate()) {
        const auto& xma = in.const_arrays();
        const auto& yma = out.arrays();
        const auto& ama = acoef.arrays();
        AMREX_D_TERM(const auto& bxma = bxcoef.const_arrays();,
                     const auto& byma = bycoef.const_arrays();,
                     const auto& bzma = bzcoef.const_arrays(););
        if (this->m_overset_mask[amrlev][mglev]) {
            const auto& osmma = this->m_overset_mask[amrlev][mglev]->const_arrays();
            ParallelFor(out, IntVect(0), ncomp,
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
            {
                mlabeclap_adotx_os(i,j,k,n, yma[box_no], xma[box_no], ama[box_no],
                                   AMREX_D_DECL(bxma[box_no],byma[box_no],bzma[box_no]),
                                   osmma[box_no], dxinv, ascalar, bscalar);
            });
        } else {
            ParallelFor(out, IntVect(0), ncomp,
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
            {
                mlabeclap_adotx(i,j,k,n, yma[box_no], xma[box_no], ama[box_no],
                                AMREX_D_DECL(bxma[box_no],byma[box_no],bzma[box_no]),
                                dxinv, ascalar, bscalar);
            });
        }
        Gpu::streamSynchronize();
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(out, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            const auto& xfab = in.array(mfi);
            const auto& yfab = out.array(mfi);
            const auto& afab = acoef.array(mfi);
            AMREX_D_TERM(const auto& bxfab = bxcoef.array(mfi);,
                         const auto& byfab = bycoef.array(mfi);,
                         const auto& bzfab = bzcoef.array(mfi););
            if (this->m_overset_mask[amrlev][mglev]) {
                const auto& osm = this->m_overset_mask[amrlev][mglev]->const_array(mfi);
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D(bx, ncomp, i, j, k, n,
                {
                    mlabeclap_adotx_os(i,j,k,n, yfab, xfab, afab, AMREX_D_DECL(bxfab,byfab,bzfab),
                                       osm, dxinv, ascalar, bscalar);
                });
            } else {
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D(bx, ncomp, i, j, k, n,
                {
                    mlabeclap_adotx(i,j,k,n, yfab, xfab, afab, AMREX_D_DECL(bxfab,byfab,bzfab),
                                    dxinv, ascalar, bscalar);
                });
            }
        }
    }
}

template <typename MF>
void
MLABecLaplacianT<MF>::Fsmooth (int amrlev, int mglev, MF& sol, const MF& rhs, int redblack) const
{
    BL_PROFILE("MLABecLaplacian::Fsmooth()");

    bool regular_coarsening = true;
    if (amrlev == 0 && mglev > 0) {
        regular_coarsening = this->mg_coarsen_ratio_vec[mglev-1] == this->mg_coarsen_ratio;
    }

    const MF& acoef = m_a_coeffs[amrlev][mglev];
    AMREX_ALWAYS_ASSERT(acoef.nGrowVect() == 0);
    AMREX_D_TERM(const MF& bxcoef = m_b_coeffs[amrlev][mglev][0];,
                 const MF& bycoef = m_b_coeffs[amrlev][mglev][1];,
                 const MF& bzcoef = m_b_coeffs[amrlev][mglev][2];);
    const auto& undrrelxr = this->m_undrrelxr[amrlev][mglev];
    const auto& maskvals  = this->m_maskvals [amrlev][mglev];

    OrientationIter oitr;

    const auto& f0 = undrrelxr[oitr()]; ++oitr;
    const auto& f1 = undrrelxr[oitr()]; ++oitr;
#if (AMREX_SPACEDIM > 1)
    const auto& f2 = undrrelxr[oitr()]; ++oitr;
    const auto& f3 = undrrelxr[oitr()]; ++oitr;
#if (AMREX_SPACEDIM > 2)
    const auto& f4 = undrrelxr[oitr()]; ++oitr;
    const auto& f5 = undrrelxr[oitr()]; ++oitr;
#endif
#endif

    const MultiMask& mm0 = maskvals[0];
    const MultiMask& mm1 = maskvals[1];
#if (AMREX_SPACEDIM > 1)
    const MultiMask& mm2 = maskvals[2];
    const MultiMask& mm3 = maskvals[3];
#if (AMREX_SPACEDIM > 2)
    const MultiMask& mm4 = maskvals[4];
    const MultiMask& mm5 = maskvals[5];
#endif
#endif

    const int nc = this->getNComp();
    const Real* h = this->m_geom[amrlev][mglev].CellSize();
    AMREX_D_TERM(const RT dhx = m_b_scalar/static_cast<RT>(h[0]*h[0]);,
                 const RT dhy = m_b_scalar/static_cast<RT>(h[1]*h[1]);,
                 const RT dhz = m_b_scalar/static_cast<RT>(h[2]*h[2]));
    const RT alpha = m_a_scalar;

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion() && sol.isFusingCandidate()
        && (this->m_overset_mask[amrlev][mglev] || regular_coarsening))
    {
        const auto& m0ma = mm0.const_arrays();
        const auto& m1ma = mm1.const_arrays();
#if (AMREX_SPACEDIM > 1)
        const auto& m2ma = mm2.const_arrays();
        const auto& m3ma = mm3.const_arrays();
#if (AMREX_SPACEDIM > 2)
        const auto& m4ma = mm4.const_arrays();
        const auto& m5ma = mm5.const_arrays();
#endif
#endif

        const auto& solnma = sol.arrays();
        const auto& rhsma = rhs.const_arrays();
        const auto& ama = acoef.const_arrays();

        AMREX_D_TERM(const auto& bxma = bxcoef.const_arrays();,
                     const auto& byma = bycoef.const_arrays();,
                     const auto& bzma = bzcoef.const_arrays(););

        const auto& f0ma = f0.const_arrays();
        const auto& f1ma = f1.const_arrays();
#if (AMREX_SPACEDIM > 1)
        const auto& f2ma = f2.const_arrays();
        const auto& f3ma = f3.const_arrays();
#if (AMREX_SPACEDIM > 2)
        const auto& f4ma = f4.const_arrays();
        const auto& f5ma = f5.const_arrays();
#endif
#endif

        if (this->m_overset_mask[amrlev][mglev]) {
            const auto& osmma = this->m_overset_mask[amrlev][mglev]->const_arrays();
            ParallelFor(sol, IntVect(0), nc,
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
            {
                Box vbx(ama[box_no]);
                abec_gsrb_os(i,j,k,n, solnma[box_no], rhsma[box_no], alpha, ama[box_no],
                             AMREX_D_DECL(dhx, dhy, dhz),
                             AMREX_D_DECL(bxma[box_no],byma[box_no],bzma[box_no]),
                             AMREX_D_DECL(m0ma[box_no],m2ma[box_no],m4ma[box_no]),
                             AMREX_D_DECL(m1ma[box_no],m3ma[box_no],m5ma[box_no]),
                             AMREX_D_DECL(f0ma[box_no],f2ma[box_no],f4ma[box_no]),
                             AMREX_D_DECL(f1ma[box_no],f3ma[box_no],f5ma[box_no]),
                             osmma[box_no], vbx, redblack);
            });
        } else if (regular_coarsening) {
            ParallelFor(sol, IntVect(0), nc,
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
            {
                Box vbx(ama[box_no]);
                abec_gsrb(i,j,k,n, solnma[box_no], rhsma[box_no], alpha, ama[box_no],
                          AMREX_D_DECL(dhx, dhy, dhz),
                          AMREX_D_DECL(bxma[box_no],byma[box_no],bzma[box_no]),
                          AMREX_D_DECL(m0ma[box_no],m2ma[box_no],m4ma[box_no]),
                          AMREX_D_DECL(m1ma[box_no],m3ma[box_no],m5ma[box_no]),
                          AMREX_D_DECL(f0ma[box_no],f2ma[box_no],f4ma[box_no]),
                          AMREX_D_DECL(f1ma[box_no],f3ma[box_no],f5ma[box_no]),
                          vbx, redblack);
            });
        }
        Gpu::streamSynchronize();
    } else
#endif
    {
        MFItInfo mfi_info;
        if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling().SetDynamic(true);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(sol,mfi_info); mfi.isValid(); ++mfi)
        {
            const auto& m0 = mm0.array(mfi);
            const auto& m1 = mm1.array(mfi);
#if (AMREX_SPACEDIM > 1)
            const auto& m2 = mm2.array(mfi);
            const auto& m3 = mm3.array(mfi);
#if (AMREX_SPACEDIM > 2)
            const auto& m4 = mm4.array(mfi);
            const auto& m5 = mm5.array(mfi);
#endif
#endif

            const Box& tbx = mfi.tilebox();
            const Box& vbx = mfi.validbox();
            const auto& solnfab = sol.array(mfi);
            const auto& rhsfab  = rhs.const_array(mfi);
            const auto& afab    = acoef.const_array(mfi);

            AMREX_D_TERM(const auto& bxfab = bxcoef.const_array(mfi);,
                         const auto& byfab = bycoef.const_array(mfi);,
                         const auto& bzfab = bzcoef.const_array(mfi););

            const auto& f0fab = f0.const_array(mfi);
            const auto& f1fab = f1.const_array(mfi);
#if (AMREX_SPACEDIM > 1)
            const auto& f2fab = f2.const_array(mfi);
            const auto& f3fab = f3.const_array(mfi);
#if (AMREX_SPACEDIM > 2)
            const auto& f4fab = f4.const_array(mfi);
            const auto& f5fab = f5.const_array(mfi);
#endif
#endif

            if (this->m_overset_mask[amrlev][mglev]) {
                const auto& osm = this->m_overset_mask[amrlev][mglev]->const_array(mfi);
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D(tbx, nc, i, j, k, n,
                {
                    abec_gsrb_os(i,j,k,n, solnfab, rhsfab, alpha, afab,
                                 AMREX_D_DECL(dhx, dhy, dhz),
                                 AMREX_D_DECL(bxfab, byfab, bzfab),
                                 AMREX_D_DECL(m0,m2,m4),
                                 AMREX_D_DECL(m1,m3,m5),
                                 AMREX_D_DECL(f0fab,f2fab,f4fab),
                                 AMREX_D_DECL(f1fab,f3fab,f5fab),
                                 osm, vbx, redblack);
                });
            } else if (regular_coarsening) {
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D(tbx, nc, i, j, k, n,
                {
                    abec_gsrb(i,j,k,n, solnfab, rhsfab, alpha, afab,
                              AMREX_D_DECL(dhx, dhy, dhz),
                              AMREX_D_DECL(bxfab, byfab, bzfab),
                              AMREX_D_DECL(m0,m2,m4),
                              AMREX_D_DECL(m1,m3,m5),
                              AMREX_D_DECL(f0fab,f2fab,f4fab),
                              AMREX_D_DECL(f1fab,f3fab,f5fab),
                              vbx, redblack);
                });
            } else {
                Gpu::LaunchSafeGuard lsg(false); // xxxxx gpu todo
                // line solve does not with with GPU
                AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( tbx, thread_box,
                {
                    abec_gsrb_with_line_solve(thread_box, solnfab, rhsfab, alpha, afab,
                                              AMREX_D_DECL(dhx, dhy, dhz),
                                              AMREX_D_DECL(bxfab, byfab, bzfab),
                                              AMREX_D_DECL(m0,m2,m4),
                                              AMREX_D_DECL(m1,m3,m5),
                                              AMREX_D_DECL(f0fab,f2fab,f4fab),
                                              AMREX_D_DECL(f1fab,f3fab,f5fab),
                                              vbx, redblack, nc);
                });
            }
        }
    }
}

template <typename MF>
void
MLABecLaplacianT<MF>::FFlux (int amrlev, const MFIter& mfi,
                             const Array<FAB*,AMREX_SPACEDIM>& flux,
                             const FAB& sol, Location, const int face_only) const
{
    BL_PROFILE("MLABecLaplacian::FFlux()");

    const int mglev = 0;
    const Box& box = mfi.tilebox();
    const Real* dxinv = this->m_geom[amrlev][mglev].InvCellSize();
    const int ncomp = this->getNComp();
    FFlux(box, dxinv, m_b_scalar,
          Array<FAB const*,AMREX_SPACEDIM>{{AMREX_D_DECL(&(m_b_coeffs[amrlev][mglev][0][mfi]),
                                                         &(m_b_coeffs[amrlev][mglev][1][mfi]),
                                                         &(m_b_coeffs[amrlev][mglev][2][mfi]))}},
          flux, sol, face_only, ncomp);
}

template <typename MF>
void
MLABecLaplacianT<MF>::FFlux (Box const& box, Real const* dxinv, RT bscalar,
                             Array<FAB const*, AMREX_SPACEDIM> const& bcoef,
                             Array<FAB*,AMREX_SPACEDIM> const& flux,
                             FAB const& sol, int face_only, int ncomp)
{
    AMREX_D_TERM(const auto bx = bcoef[0]->const_array();,
                 const auto by = bcoef[1]->const_array();,
                 const auto bz = bcoef[2]->const_array(););
    AMREX_D_TERM(const auto& fxarr = flux[0]->array();,
                 const auto& fyarr = flux[1]->array();,
                 const auto& fzarr = flux[2]->array(););
    const auto& solarr = sol.array();

    if (face_only)
    {
        RT fac = bscalar*static_cast<RT>(dxinv[0]);
        Box blo = amrex::bdryLo(box, 0);
        int blen = box.length(0);
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( blo, tbox,
        {
            mlabeclap_flux_xface(tbox, fxarr, solarr, bx, fac, blen, ncomp);
        });
#if (AMREX_SPACEDIM >= 2)
        fac = bscalar*static_cast<RT>(dxinv[1]);
        blo = amrex::bdryLo(box, 1);
        blen = box.length(1);
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( blo, tbox,
        {
            mlabeclap_flux_yface(tbox, fyarr, solarr, by, fac, blen, ncomp);
        });
#endif
#if (AMREX_SPACEDIM == 3)
        fac = bscalar*static_cast<RT>(dxinv[2]);
        blo = amrex::bdryLo(box, 2);
        blen = box.length(2);
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( blo, tbox,
        {
            mlabeclap_flux_zface(tbox, fzarr, solarr, bz, fac, blen, ncomp);
        });
#endif
    }
    else
    {
        RT fac = bscalar*static_cast<RT>(dxinv[0]);
        Box bflux = amrex::surroundingNodes(box, 0);
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bflux, tbox,
        {
            mlabeclap_flux_x(tbox, fxarr, solarr, bx, fac, ncomp);
        });
#if (AMREX_SPACEDIM >= 2)
        fac = bscalar*static_cast<RT>(dxinv[1]);
        bflux = amrex::surroundingNodes(box, 1);
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bflux, tbox,
        {
            mlabeclap_flux_y(tbox, fyarr, solarr, by, fac, ncomp);
        });
#endif
#if (AMREX_SPACEDIM == 3)
        fac = bscalar*static_cast<RT>(dxinv[2]);
        bflux = amrex::surroundingNodes(box, 2);
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bflux, tbox,
        {
            mlabeclap_flux_z(tbox, fzarr, solarr, bz, fac, ncomp);
        });
#endif
    }
}

template <typename MF>
void
MLABecLaplacianT<MF>::normalize (int amrlev, int mglev, MF& mf) const
{
    BL_PROFILE("MLABecLaplacian::normalize()");

    const auto& acoef = m_a_coeffs[amrlev][mglev];
    AMREX_D_TERM(const auto& bxcoef = m_b_coeffs[amrlev][mglev][0];,
                 const auto& bycoef = m_b_coeffs[amrlev][mglev][1];,
                 const auto& bzcoef = m_b_coeffs[amrlev][mglev][2];);

    const GpuArray<RT,AMREX_SPACEDIM> dxinv
        {AMREX_D_DECL(static_cast<RT>(this->m_geom[amrlev][mglev].InvCellSize(0)),
                      static_cast<RT>(this->m_geom[amrlev][mglev].InvCellSize(1)),
                      static_cast<RT>(this->m_geom[amrlev][mglev].InvCellSize(2)))};

    const RT ascalar = m_a_scalar;
    const RT bscalar = m_b_scalar;

    const int ncomp = getNComp();

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion() && mf.isFusingCandidate()) {
        const auto& ma = mf.arrays();
        const auto& ama = acoef.const_arrays();
        AMREX_D_TERM(const auto& bxma = bxcoef.const_arrays();,
                     const auto& byma = bycoef.const_arrays();,
                     const auto& bzma = bzcoef.const_arrays(););
        ParallelFor(mf, IntVect(0), ncomp,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
        {
            mlabeclap_normalize(i,j,k,n, ma[box_no], ama[box_no],
                                AMREX_D_DECL(bxma[box_no],byma[box_no],bzma[box_no]),
                                dxinv, ascalar, bscalar);
        });
        Gpu::streamSynchronize();
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(mf, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            const auto& fab = mf.array(mfi);
            const auto& afab = acoef.array(mfi);
            AMREX_D_TERM(const auto& bxfab = bxcoef.array(mfi);,
                         const auto& byfab = bycoef.array(mfi);,
                         const auto& bzfab = bzcoef.array(mfi););

            AMREX_HOST_DEVICE_PARALLEL_FOR_4D(bx, ncomp, i, j, k, n,
            {
                mlabeclap_normalize(i,j,k,n, fab, afab, AMREX_D_DECL(bxfab,byfab,bzfab),
                                    dxinv, ascalar, bscalar);
            });
        }
    }
}

template <typename MF>
bool
MLABecLaplacianT<MF>::supportNSolve () const
{
    bool support = false;
    if (this->m_overset_mask[0][0]) {
        if (this->m_geom[0].back().Domain().coarsenable(MLLinOp::mg_coarsen_ratio,
                                                        MLLinOp::mg_domain_min_width)
            && this->m_grids[0].back().coarsenable(MLLinOp::mg_coarsen_ratio, MLLinOp::mg_box_min_width))
        {
            support = true;
        }
    }
    return support;
}

template <typename MF>
std::unique_ptr<MLLinOpT<MF>>
MLABecLaplacianT<MF>::makeNLinOp (int /*grid_size*/) const
{
    if (this->m_overset_mask[0][0] == nullptr) return nullptr;

    const Geometry& geom = this->m_geom[0].back();
    const BoxArray& ba = this->m_grids[0].back();
    const DistributionMapping& dm = this->m_dmap[0].back();

    std::unique_ptr<MLLinOpT<MF>> r
        {new MLABecLaplacianT<MF>({geom}, {ba}, {dm}, this->m_lpinfo_arg)};

    auto nop = dynamic_cast<MLABecLaplacianT<MF>*>(r.get());
    if (!nop) {
        return nullptr;
    }

    nop->m_parent = this;

    nop->setMaxOrder(this->maxorder);
    nop->setVerbose(this->verbose);

    nop->setDomainBC(this->m_lobc, this->m_hibc);

    if (this->needsCoarseDataForBC())
    {
        const Real* dx0 = this->m_geom[0][0].CellSize();
        const Real fac = Real(0.5)*this->m_coarse_data_crse_ratio;
        RealVect cbloc {AMREX_D_DECL(dx0[0]*fac, dx0[1]*fac, dx0[2]*fac)};
        nop->setCoarseFineBCLocation(cbloc);
    }

    nop->setScalars(m_a_scalar, m_b_scalar);

    MF const& alpha_bottom = m_a_coeffs[0].back();
    iMultiFab const& osm_bottom = *(this->m_overset_mask[0].back());
    const int ncomp = alpha_bottom.nComp();
    MF alpha(ba, dm, ncomp, 0);

    RT a_max = alpha_bottom.norminf(0, ncomp, IntVect(0), true, true);
    const int ncomp_b = m_b_coeffs[0].back()[0].nComp();
    AMREX_D_TERM(RT bx_max = m_b_coeffs[0].back()[0].norminf(0,ncomp_b,IntVect(0),true,true);,
                 RT by_max = m_b_coeffs[0].back()[1].norminf(0,ncomp_b,IntVect(0),true,true);,
                 RT bz_max = m_b_coeffs[0].back()[2].norminf(0,ncomp_b,IntVect(0),true,true));
    const GpuArray<RT,AMREX_SPACEDIM> dxinv
        {AMREX_D_DECL(static_cast<RT>(geom.InvCellSize(0)),
                      static_cast<RT>(geom.InvCellSize(1)),
                      static_cast<RT>(geom.InvCellSize(2)))};
    RT huge_alpha = RT(1.e30) *
        amrex::max(a_max*std::abs(m_a_scalar),
                   AMREX_D_DECL(std::abs(m_b_scalar)*bx_max*dxinv[0]*dxinv[0],
                                std::abs(m_b_scalar)*by_max*dxinv[1]*dxinv[1],
                                std::abs(m_b_scalar)*bz_max*dxinv[2]*dxinv[2]));
    ParallelAllReduce::Max(huge_alpha, ParallelContext::CommunicatorSub());

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion() && alpha.isFusingCandidate()) {
        auto const& ama = alpha.arrays();
        auto const& abotma = alpha_bottom.const_arrays();
        auto const& mma = osm_bottom.const_arrays();
        ParallelFor(alpha, IntVect(0), ncomp,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
        {
            if (mma[box_no](i,j,k)) {
                ama[box_no](i,j,k,n) = abotma[box_no](i,j,k,n);
            } else {
                ama[box_no](i,j,k,n) = huge_alpha;
            }
        });
        Gpu::streamSynchronize();
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(alpha,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            Box const& bx = mfi.tilebox();
            auto const& a = alpha.array(mfi);
            auto const& abot = alpha_bottom.const_array(mfi);
            auto const& m = osm_bottom.const_array(mfi);
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D(bx, ncomp, i, j, k, n,
            {
                if (m(i,j,k)) {
                    a(i,j,k,n) = abot(i,j,k,n);
                } else {
                    a(i,j,k,n) = huge_alpha;
                }
            });
        }
    }

    nop->setACoeffs(0, alpha);
    nop->setBCoeffs(0, GetArrOfConstPtrs(m_b_coeffs[0].back()));

    return r;
}

template <typename MF>
void
MLABecLaplacianT<MF>::copyNSolveSolution (MF& dst, MF const& src) const
{
    if (this->m_overset_mask[0].back() == nullptr) return;

    const int ncomp = dst.nComp();

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion() && dst.isFusingCandidate()) {
        auto const& dstma = dst.arrays();
        auto const& srcma = src.const_arrays();
        auto const& mma = this->m_overset_mask[0].back()->const_arrays();
        ParallelFor(dst, IntVect(0), ncomp,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
        {
            if (mma[box_no](i,j,k)) {
                dstma[box_no](i,j,k,n) = srcma[box_no](i,j,k,n);
            } else {
                dstma[box_no](i,j,k,n) = RT(0.0);
            }
        });
        Gpu::streamSynchronize();
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(dst,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            Box const& bx = mfi.tilebox();
            auto const& dfab = dst.array(mfi);
            auto const& sfab = src.const_array(mfi);
            auto const& m = this->m_overset_mask[0].back()->const_array(mfi);
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D(bx, ncomp, i, j, k, n,
            {
                if (m(i,j,k)) {
                    dfab(i,j,k,n) = sfab(i,j,k,n);
                } else {
                    dfab(i,j,k,n) = RT(0.0);
                }
            });
        }
    }
}

extern template class MLABecLaplacianT<MultiFab>;

using MLABecLaplacian = MLABecLaplacianT<MultiFab>;

}

#endif
