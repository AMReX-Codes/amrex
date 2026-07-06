#ifndef AMREX_FFT_STOKES_H_
#define AMREX_FFT_STOKES_H_

#include <AMReX_FFT.H>
#include <AMReX_Geometry.H>

namespace amrex::FFT
{

/**
 * \ingroup amrex_fft
 * \brief Stokes solver for periodic domains using FFT.
 */
template <typename MF = MultiFab>
class Stokes
{
public:

    static_assert(AMREX_SPACEDIM >= 2, "FFT::Stokes requires 2D or 3D");

    using cMF = FabArray<BaseFab<GpuComplex<typename MF::value_type>>>;

    Stokes () = default;

    /**
     * \brief Construct a Stokes solver with explicit boundary types.
     *
     * Only periodic BCs are supported.
     *
     * \param geom Geometry describing the domain and metric information.
     * \param bc   Pair of boundary descriptors (low/high) per coordinate direction.
     */
    Stokes (Geometry const& geom,
             Array<std::pair<Boundary,Boundary>,AMREX_SPACEDIM> const& bc)
        requires (IsFabArray_v<MF>)
    {
        define(geom, bc);
    }

    /**
     * \brief Construct a purely periodic Stokes solver.
     *
     * \param geom Periodic geometry (all directions must be periodic).
     */
    explicit Stokes (Geometry const& geom)
        requires (IsFabArray_v<MF>)
    {
        define(geom);
    }

    /**
     * \brief Define a Stokes solver with explicit boundary types.
     *
     * Only periodic BCs are supported.
     *
     * \param geom Geometry describing the domain and metric information.
     * \param bc   Pair of boundary descriptors (low/high) per coordinate direction.
     */
    void define (Geometry const& geom,
                 Array<std::pair<Boundary,Boundary>,AMREX_SPACEDIM> const& bc)
        requires (IsFabArray_v<MF>)
    {
        m_domain_lo = geom.Domain().smallEnd();
        m_geom = detail::shift_geom(geom);
        m_bc = bc;

        bool all_periodic = true;
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            all_periodic = all_periodic
                && (bc[idim].first == Boundary::periodic)
                && (bc[idim].second == Boundary::periodic);
        }
        if (!all_periodic) {
            amrex::Abort("FFT::Stokes: only supports periodic BC");
        }

        m_r2c = std::make_unique<R2C<typename MF::value_type>>(m_geom.Domain());
    }

    /**
     * \brief Define a purely periodic Stokes solver.
     *
     * \param geom Periodic geometry (all directions must be periodic).
     */
    void define (Geometry const& geom)
        requires (IsFabArray_v<MF>)
    {
        define(geom, {AMREX_D_DECL(std::make_pair(Boundary::periodic,Boundary::periodic),
                                    std::make_pair(Boundary::periodic,Boundary::periodic),
                                    std::make_pair(Boundary::periodic,Boundary::periodic))});
    }

    /**
     * \brief Solve the generalized Stokes problem in spectral space.
     *
     * Solves \f$(\alpha - \eta \nabla^2) \mathbf{u} + \nabla p = \mathbf{f}\f$
     * with the divergence-free constraint enforced in Fourier space.
     * This solver is for 2D and 3D only.
     *
     * This solver uses a MAC staggered grid. \p U and \p rhsx must be
     * at x-face centers, \p V and \p rhsy at y-face centers,
     * \p W and \p rhsz at z-face centers, and \p p cell-centered.
     *
     * \param U     x-component of velocity (output, at x-face centers).
     * \param V     y-component of velocity (output, at y-face centers).
     * \param W     z-component of velocity (output, at z-face centers, 3D only).
     * \param p     Pressure (output, cell-centered).
     * \param rhsx  x-component of the right-hand side (at x-face centers).
     * \param rhsy  y-component of the right-hand side (at y-face centers).
     * \param rhsz  z-component of the right-hand side (at z-face centers, 3D only).
     * \param alpha Coefficient multiplying the velocity.
     * \param eta   Coefficient multiplying the velocity Laplacian.
     */
    void solve (AMREX_D_DECL(MF& U, MF& V, MF& W), MF& p,
                AMREX_D_DECL(MF const& rhsx, MF const& rhsy, MF const& rhsz),
                typename MF::value_type alpha, typename MF::value_type eta);

private:
    IntVect m_domain_lo;
    Geometry m_geom;
    Array<std::pair<Boundary,Boundary>,AMREX_SPACEDIM> m_bc;
    std::unique_ptr<R2C<typename MF::value_type>> m_r2c;
};


template <typename MF>
void Stokes<MF>::solve (MF& U,
                        MF& V,
#if (BL_SPACEDIM == 3)
                        MF& W,
#endif
                        MF& p,
                        MF const& rhsx,
                        MF const& rhsy,
#if (BL_SPACEDIM == 3)
                        MF const& rhsz,
#endif
                        typename MF::value_type alpha,
                        typename MF::value_type eta)
{
    BL_PROFILE("FFT::Stokes::solve");

    using T = typename MF::value_type;

    if (!m_r2c) {
        amrex::Abort("FFT::Stokes::solve called before define()");
    }

    AMREX_ASSERT(p.ixType() == IndexType::TheCellType() &&
                 AMREX_D_TERM(U.ixType() == IndexType(IntVect::TheDimensionVector(0)) &&
                              rhsx.ixType() == U.ixType(),
                           && V.ixType() == IndexType(IntVect::TheDimensionVector(1)) &&
                              rhsy.ixType() == V.ixType(),
                           && W.ixType() == IndexType(IntVect::TheDimensionVector(2)) &&
                              rhsz.ixType() == W.ixType()));

    MF* Umf = &U;
    MF* Vmf = &V;
#if (BL_SPACEDIM == 3)
    MF* Wmf = &W;
#endif
    MF* pmf = &p;
    MF const* rhsxmf = &rhsx;
    MF const* rhsymf = &rhsy;
#if (BL_SPACEDIM == 3)
    MF const* rhszmf = &rhsz;
#endif
    MF Utmp, rhsxtmp;
    MF Vtmp, rhsytmp;
#if (BL_SPACEDIM == 3)
    MF Wtmp, rhsztmp;
#endif
    MF ptmp;
    if (m_domain_lo != 0) {
        detail::shift_mfs(m_domain_lo, U, rhsx, Utmp, rhsxtmp);
        detail::shift_mfs(m_domain_lo, V, rhsy, Vtmp, rhsytmp);
#if (BL_SPACEDIM == 3)
        detail::shift_mfs(m_domain_lo, W, rhsz, Wtmp, rhsztmp);
#endif
        detail::shift_mf(m_domain_lo, p, ptmp);
        Umf = &Utmp;
        Vmf = &Vtmp;
#if (BL_SPACEDIM == 3)
        Wmf = &Wtmp;
#endif
        pmf = &ptmp;
        rhsxmf = &rhsxtmp;
        rhsymf = &rhsytmp;
#if (BL_SPACEDIM == 3)
        rhszmf = &rhsztmp;
#endif
    }

    auto& r2c = *m_r2c;
    auto const& dxinv = m_geom.InvCellSizeArray();
    auto const scaling = r2c.scalingFactor();
    auto const& [cba, cdm] = r2c.getSpectralDataLayout();

    cMF phat(cba, cdm, 1, 0);

    cMF rxhat(cba,cdm,1,0);
    r2c.forward(*rhsxmf, rxhat);

    cMF ryhat(cba,cdm,1,0);
    r2c.forward(*rhsymf, ryhat);

#if (BL_SPACEDIM == 3)
    cMF rzhat(cba,cdm,1,0);
    r2c.forward(*rhszmf, rzhat);
#endif

    using Complex = GpuComplex<T>;
    T constexpr tol = std::numeric_limits<T>::epsilon() * T(10);
    int const nx = m_geom.Domain().length(0);
#if (AMREX_SPACEDIM >= 2)
    int const ny = m_geom.Domain().length(1);
#endif
#if (AMREX_SPACEDIM == 3)
    int const nz = m_geom.Domain().length(2);
#endif
    for (MFIter mfi(phat); mfi.isValid(); ++mfi) {
        auto const& pb = phat[mfi].box();
        auto const& parr = phat[mfi].array();
        auto const& rx = rxhat[mfi].array();
        auto const& ry = ryhat[mfi].array();
#if (BL_SPACEDIM == 3)
        auto const& rz = rzhat[mfi].array();
#endif
        AMREX_D_TERM(T kwx = T(2)*Math::pi<T>()/T(nx);,
                     T kwy = T(2)*Math::pi<T>()/T(ny);,
                     T kwz = T(2)*Math::pi<T>()/T(nz);)
        ParallelFor(pb, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            AMREX_D_TERM(int ik = (i <= nx/2) ? i : i - nx;,
                         int jk = (j <= ny/2) ? j : j - ny;,
                         int kk = (k <= nz/2) ? k : k - nz);

            GpuArray<T,AMREX_SPACEDIM> kwave{AMREX_D_DECL(ik*kwx,jk*kwy,kk*kwz)};

            T delsqk = AMREX_D_TERM(T(2)*(std::cos(kwave[0])-T(1))*(dxinv[0]*dxinv[0]),
                                  + T(2)*(std::cos(kwave[1])-T(1))*(dxinv[1]*dxinv[1]),
                                  + T(2)*(std::cos(kwave[2])-T(1))*(dxinv[2]*dxinv[2]));

            GpuArray<Complex, AMREX_SPACEDIM> delkp
                {AMREX_D_DECL(Complex((std::cos(kwave[0])-T(1))*dxinv[0],
                                    std::sin(kwave[0])         *dxinv[0]),
                              Complex((std::cos(kwave[1])-T(1))*dxinv[1],
                                    std::sin(kwave[1])         *dxinv[1]),
                              Complex((std::cos(kwave[2])-T(1))*dxinv[2],
                                    std::sin(kwave[2])         *dxinv[2]))};

            GpuArray<Complex, AMREX_SPACEDIM> delkm
                {AMREX_D_DECL(Complex((T(1)-std::cos(kwave[0]))*dxinv[0],
                                    std::sin(kwave[0])         *dxinv[0]),
                              Complex((T(1)-std::cos(kwave[1]))*dxinv[1],
                                    std::sin(kwave[1])         *dxinv[1]),
                              Complex((T(1)-std::cos(kwave[2]))*dxinv[2],
                                    std::sin(kwave[2])         *dxinv[2]))};

            AMREX_D_TERM(Complex const rxk = rx(i,j,k);,
                         Complex const ryk = ry(i,j,k);,
                         Complex const rzk = rz(i,j,k);)

            Complex rhsdotdp = scaling * (AMREX_D_TERM(rxk*delkp[0],
                                                      +ryk*delkp[1],
                                                      +rzk*delkp[2]));

            if (std::abs(delsqk) > tol) {
                parr(i,j,k)= rhsdotdp / delsqk;
            } else {
                parr(i,j,k)= T(0.);
            }

            T diffop = alpha - eta*delsqk;
            if (diffop > tol) {
                rx(i,j,k) = (scaling*rxk - parr(i,j,k)*delkm[0])/diffop;
                ry(i,j,k) = (scaling*ryk - parr(i,j,k)*delkm[1])/diffop;
#if (BL_SPACEDIM == 3)
                rz(i,j,k) = (scaling*rzk - parr(i,j,k)*delkm[2])/diffop;
#endif
            } else {
                rx(i,j,k) = T(0.);
                ry(i,j,k) = T(0.);
#if (BL_SPACEDIM == 3)
                rz(i,j,k) = T(0.);
#endif
            }

        });
    }

    r2c.backward(phat, *pmf);
    r2c.backward(rxhat, *Umf);
    r2c.backward(ryhat, *Vmf);
#if (BL_SPACEDIM == 3)
    r2c.backward(rzhat, *Wmf);
#endif
}

} // namespace amrex::FFT

#endif
