#ifndef AMREX_FFT_CROSSPROJ_H_
#define AMREX_FFT_CROSSPROJ_H_

#include <AMReX_FFT.H>
#include <AMReX_Geometry.H>

namespace amrex::FFT
{

/**
 * \file AMReX_FFT_CrossProj.H
 *
 * FFT-based CrossProj solvers (periodic, mixed BC, hybrid, and open BC variants).
 */

/// \cond DOXYGEN_IGNORE
namespace detail {
template <typename MF>
void fill_physbc (MF& mf, Geometry const& geom,
                  Array<std::pair<Boundary,Boundary>,AMREX_SPACEDIM> const& bc);

}
/// \endcond

/**
 * \ingroup amrex_fft
 * \brief FFT-based projection of a cell-centered velocity field onto the
 * space of fields that are divergence free with respect to the node-centered
 * cross-stencil divergence.
 *
 * The divergence at each node is computed from the 4 (2D) or 8 (3D) cells
 * surrounding it. Currently only all-periodic boundaries are supported, in
 * 2D and 3D. The velocity MF must have AMREX_SPACEDIM components.
 */
template <typename MF = MultiFab>
class CrossProj
{
    static_assert(AMREX_SPACEDIM >= 2, "FFT::CrossProj requires AMREX_SPACEDIM >= 2");
public:

    /**
     * \brief Construct a CrossProj solver with explicit boundary types.
     *
     * \param geom Geometry describing the domain and metric information.
     * \param bc   Pair of boundary descriptors (low/high) per coordinate direction.
     */
    CrossProj (Geometry const& geom,
             Array<std::pair<Boundary,Boundary>,AMREX_SPACEDIM> const& bc)
        requires (IsFabArray_v<MF>)
        : m_domain_lo(geom.Domain().smallEnd()),
          m_geom(detail::shift_geom(geom)),
          m_bc(bc)
    {
        bool all_periodic = true;
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            all_periodic = all_periodic
                && (bc[idim].first == Boundary::periodic)
                && (bc[idim].second == Boundary::periodic);
        }
        if (all_periodic) {
            Info info = Info{};
            info.setBatchSize(AMREX_SPACEDIM);
            m_r2c = std::make_unique<R2C<typename MF::value_type>>(m_geom.Domain(), info);
        } else {
            m_r2x = std::make_unique<R2X<typename MF::value_type>> (m_geom.Domain(), m_bc);
        }
    }

    /**
     * \brief Construct a purely periodic CrossProj solver.
     *
     * \param geom Periodic geometry (all directions must be periodic).
     */
    explicit CrossProj (Geometry const& geom)
        requires (IsFabArray_v<MF>)
        : m_domain_lo(geom.Domain().smallEnd()),
          m_geom(detail::shift_geom(geom)),
          m_bc{AMREX_D_DECL(std::make_pair(Boundary::periodic,Boundary::periodic),
                            std::make_pair(Boundary::periodic,Boundary::periodic),
                            std::make_pair(Boundary::periodic,Boundary::periodic))}
    {
        if (m_geom.isAllPeriodic()) {
            Info info = Info{};
            info.setBatchSize(AMREX_SPACEDIM);
            m_r2c = std::make_unique<R2C<typename MF::value_type>>(m_geom.Domain(), info);
        } else {
            amrex::Abort("FFT::CrossProj: wrong BC");
        }
    }

    /**
     * \brief Project the velocity in \p a_rhs and store the result in \p a_soln.
     *
     * If \p soln has ghost cells, one layer of ghost cells will be filled.
     * The two MFs could be the same MF.
     *
     * \param a_soln Projected velocity.
     * \param a_rhs  Velocity to be projected.
     */
    void solve (MF& a_soln, MF const& a_rhs);

private:
    IntVect m_domain_lo;
    Geometry m_geom;
    Array<std::pair<Boundary,Boundary>,AMREX_SPACEDIM> m_bc;
    std::unique_ptr<R2X<typename MF::value_type>> m_r2x;
    std::unique_ptr<R2C<typename MF::value_type>> m_r2c;
};

template <typename MF>
void CrossProj<MF>::solve (MF& a_soln, MF const& a_rhs)
{
    BL_PROFILE("FFT::CrossProj::solve");

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_r2c != nullptr,
        "FFT::CrossProj: only all-periodic boundaries are currently supported");

    MF* soln = &a_soln;
    MF const* rhs = &a_rhs;
    MF solntmp, rhstmp;
    if (m_domain_lo != 0) {
        detail::shift_mfs(m_domain_lo, a_soln, a_rhs, solntmp, rhstmp);
        soln = &solntmp;
        rhs = &rhstmp;
    }

    AMREX_ASSERT(soln->is_cell_centered() && rhs->is_cell_centered());

    using T = typename MF::value_type;

    GpuArray<T,AMREX_SPACEDIM> fac
        {AMREX_D_DECL(Math::pi<T>()/(T(m_geom.Domain().length(0))),
                      Math::pi<T>()/(T(m_geom.Domain().length(1))),
                      Math::pi<T>()/(T(m_geom.Domain().length(2))))};
    GpuArray<T,AMREX_SPACEDIM> dxfac
        {AMREX_D_DECL(T(1)/T(m_geom.CellSize(0)),
                      T(1)/T(m_geom.CellSize(1)),
                      T(1)/T(m_geom.CellSize(2)))};
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (m_geom.Domain().length(idim) == 1) {
            dxfac[idim] = 0;
        }
    }
    auto scale = (m_r2x) ? m_r2x->scalingFactor() : m_r2c->scalingFactor();

    GpuArray<T,AMREX_SPACEDIM> offset{AMREX_D_DECL(T(0),T(0),T(0))};
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (m_bc[idim].first == Boundary::odd &&
            m_bc[idim].second == Boundary::odd)
        {
            offset[idim] = T(1);
        }
        else if ((m_bc[idim].first == Boundary::odd &&
                  m_bc[idim].second == Boundary::even) ||
                 (m_bc[idim].first == Boundary::even &&
                  m_bc[idim].second == Boundary::odd))
        {
            offset[idim] = T(0.5);
        }
    }

    GpuArray<int,AMREX_SPACEDIM> len{AMREX_D_DECL(m_geom.Domain().length(0),
                                                  m_geom.Domain().length(1),
                                                  m_geom.Domain().length(2))};

    auto f = [=] AMREX_GPU_DEVICE (int i, int j, int k, CellData < GpuComplex < T >> spectral_data)
    {
        amrex::ignore_unused(j,k);
        AMREX_D_TERM(T a = fac[0]*(T(i)+offset[0]);,
                     T b = fac[1]*(T(j)+offset[1]);,
                     T c = fac[2]*(T(k)+offset[2]));

        // Fourier symbol of the node-centered cross-stencil divergence, up
        // to a common factor that cancels in the projection.  The half
        // angles a, b, c make the symbol only change sign (as a whole
        // vector) for the unshifted high wavenumbers stored by R2C.
        // The symbol vanishes at the zero mode and where at least two half
        // angles are pi/2. Detect these by index, since roundoff in cos(pi/2)
        // makes the computed symbol nonzero.
#if (AMREX_SPACEDIM == 2)
        bool const null_mode = (i == 0 && j == 0)
            || (2*i == len[0] && 2*j == len[1]);
        T sx = std::sin(a)*std::cos(b)*dxfac[0];
        T sy = std::cos(a)*std::sin(b)*dxfac[1];
        if (null_mode) { sx = T(0); sy = T(0); }

        T const s2 = sx*sx + sy*sy;
        // Select the denominator instead of guarding the division (see Poisson).
        T const s2d = (s2 != T(0)) ? s2 : T(1);
        auto cdiv = (sx * spectral_data[0] + sy * spectral_data[1]) / s2d;
        spectral_data[0] -= sx * cdiv;
        spectral_data[1] -= sy * cdiv;
#else
        int const nhalf = int(2*i == len[0]) + int(2*j == len[1]) + int(2*k == len[2]);
        bool const null_mode = (i == 0 && j == 0 && k == 0) || (nhalf >= 2);
        T sx = std::sin(a)*std::cos(b)*std::cos(c)*dxfac[0];
        T sy = std::cos(a)*std::sin(b)*std::cos(c)*dxfac[1];
        T sz = std::cos(a)*std::cos(b)*std::sin(c)*dxfac[2];
        if (null_mode) { sx = T(0); sy = T(0); sz = T(0); }

        T const s2 = sx*sx + sy*sy + sz*sz;
        // Select the denominator instead of guarding the division (see Poisson).
        T const s2d = (s2 != T(0)) ? s2 : T(1);
        auto cdiv = (sx * spectral_data[0] + sy * spectral_data[1]
                     + sz * spectral_data[2]) / s2d;
        spectral_data[0] -= sx * cdiv;
        spectral_data[1] -= sy * cdiv;
        spectral_data[2] -= sz * cdiv;
#endif
        for (int n = 0; n < AMREX_SPACEDIM; ++n) {
            spectral_data[n] *= scale;
        }
    };

    IntVect const& ng = amrex::elemwiseMin(soln->nGrowVect(), IntVect(1));

    m_r2c->forward(*rhs);
    m_r2c->post_forward_doit_0(f);
    m_r2c->backward_doit(*soln, ng, m_geom.periodicity());
}


/// \cond DOXYGEN_IGNORE
namespace detail {

template <class T>
struct FFTPhysBCTag {
    Array4<T>   dfab;
    Box         dbox;
    Boundary    bc;
    Orientation face;

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Box const& box () const noexcept { return dbox; }
};

template <typename MF>
void fill_physbc (MF& mf, Geometry const& geom,
                  Array<std::pair<Boundary,Boundary>,AMREX_SPACEDIM> const& bc)
{
    using T = typename MF::value_type;
    using Tag = FFTPhysBCTag<T>;
    Vector<Tag> tags;

    for (MFIter mfi(mf, MFItInfo{}.DisableDeviceSync()); mfi.isValid(); ++mfi)
    {
        auto const& box = mfi.fabbox();
        auto const& arr = mf.array(mfi);
        for (OrientationIter oit; oit; ++oit) {
            Orientation face = oit();
            int idim = face.coordDir();
            Box b = geom.Domain();
            Boundary fbc;
            if (face.isLow()) {
                b.setRange(idim,geom.Domain().smallEnd(idim)-1);
                fbc = bc[idim].first;
            } else {
                b.setRange(idim,geom.Domain().bigEnd(idim)+1);
                fbc = bc[idim].second;
            }
            b &= box;
            if (b.ok() && fbc != Boundary::periodic) {
                tags.push_back(Tag{.dfab = arr, .dbox = b, .bc = fbc, .face = face});
            }
        }
    }

#if defined(AMREX_USE_GPU)
    amrex::ParallelFor(tags, [=] AMREX_GPU_DEVICE (int i, int j, int k,
                                                   Tag const& tag) noexcept
#else
    auto ntags = int(tags.size());
#ifdef AMREX_USE_OMP
#pragma omp parallel for
#endif
    for (int itag = 0; itag < ntags; ++itag) {
        Tag const& tag = tags[itag];
        amrex::LoopOnCpu(tag.dbox, [&] (int i, int j, int k)
#endif
        {
            int sgn = tag.face.isLow() ? 1 : -1;
            IntVect siv = IntVect(AMREX_D_DECL(i,j,k))
                +  sgn * IntVect::TheDimensionVector(tag.face.coordDir());
            if (tag.bc == Boundary::odd) {
                tag.dfab(i,j,k) = -tag.dfab(siv);
            } else { // even
                tag.dfab(i,j,k) = tag.dfab(siv);
            }
        });
#if !defined(AMREX_USE_GPU)
    }
#endif
}
}
/// \endcond

}

#endif
