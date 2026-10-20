#ifndef AMREX_FFT_CROSSPROJ_H_
#define AMREX_FFT_CROSSPROJ_H_

#include <AMReX_FFT.H>
#include <AMReX_Geometry.H>

namespace amrex::FFT
{

/**
 * \file AMReX_FFT_CrossProj.H
 *
 * FFT-based cross-stencil projection for periodic domains.
 */

/**
 * \ingroup amrex_fft
 * \brief FFT-based projection of a cell-centered velocity field onto the
 * space of fields that are divergence free with respect to the node-centered
 * cross-stencil divergence.
 *
 * The divergence at each node is computed from the 4 (2D) or 8 (3D) cells
 * surrounding it. Currently only all-periodic boundaries are supported, in
 * 2D and 3D, with at least 2 cells in each direction. The first
 * AMREX_SPACEDIM components of the MFs are the velocity.
 */
template <typename MF = MultiFab>
class CrossProj
{
    static_assert(AMREX_SPACEDIM >= 2, "FFT::CrossProj requires AMREX_SPACEDIM >= 2");
public:

    /**
     * \brief Construct a periodic CrossProj solver.
     *
     * \param geom Periodic geometry (all directions must be periodic).
     */
    explicit CrossProj (Geometry const& geom)
        requires (IsFabArray_v<MF>)
        : m_domain_lo(geom.Domain().smallEnd()),
          m_geom(detail::shift_geom(geom))
    {
        if (!m_geom.isAllPeriodic()) {
            amrex::Abort("FFT::CrossProj: only supports periodic BC");
        }
        if (m_geom.Domain().shortside() < 2) {
            amrex::Abort("FFT::CrossProj: domain must have at least 2 cells in each direction");
        }
        Info info{};
        info.setBatchSize(AMREX_SPACEDIM);
        m_r2c = std::make_unique<R2C<typename MF::value_type>>(m_geom.Domain(), info);
    }

    /**
     * \brief Construct a CrossProj solver with explicit boundary types.
     *
     * Only periodic BCs are supported.
     *
     * \param geom Periodic geometry.
     * \param bc   Pair of boundary descriptors (low/high) per coordinate direction.
     */
    CrossProj (Geometry const& geom,
               Array<std::pair<Boundary,Boundary>,AMREX_SPACEDIM> const& bc)
        requires (IsFabArray_v<MF>)
        : CrossProj(geom)
    {
        for (auto const& b : bc) {
            if (b.first != Boundary::periodic || b.second != Boundary::periodic) {
                amrex::Abort("FFT::CrossProj: only supports periodic BC");
            }
        }
    }

    /**
     * \brief Project the velocity in \p a_rhs and store the result in \p a_soln.
     *
     * If \p a_soln has ghost cells, one layer of ghost cells will be filled.
     * The two MFs could be the same MF.
     *
     * \param a_soln Projected velocity.
     * \param a_rhs  Velocity to be projected.
     */
    void solve (MF& a_soln, MF const& a_rhs);

private:
    IntVect m_domain_lo;
    Geometry m_geom;
    std::unique_ptr<R2C<typename MF::value_type>> m_r2c;
};

template <typename MF>
void CrossProj<MF>::solve (MF& a_soln, MF const& a_rhs)
{
    BL_PROFILE("FFT::CrossProj::solve");

    AMREX_ALWAYS_ASSERT(a_soln.is_cell_centered() && a_rhs.is_cell_centered() &&
                        a_soln.nComp() >= AMREX_SPACEDIM &&
                        a_rhs.nComp() >= AMREX_SPACEDIM);

    MF* soln = &a_soln;
    MF const* rhs = &a_rhs;
    MF solntmp, rhstmp;
    if (m_domain_lo != 0) {
        detail::shift_mfs(m_domain_lo, a_soln, a_rhs, solntmp, rhstmp);
        soln = &solntmp;
        rhs = &rhstmp;
    }

    using T = typename MF::value_type;

    GpuArray<T,AMREX_SPACEDIM> fac
        {AMREX_D_DECL(Math::pi<T>()/(T(m_geom.Domain().length(0))),
                      Math::pi<T>()/(T(m_geom.Domain().length(1))),
                      Math::pi<T>()/(T(m_geom.Domain().length(2))))};
    GpuArray<T,AMREX_SPACEDIM> const dxfac
        {AMREX_D_DECL(T(1)/T(m_geom.CellSize(0)),
                      T(1)/T(m_geom.CellSize(1)),
                      T(1)/T(m_geom.CellSize(2)))};
    auto scale = m_r2c->scalingFactor();

    GpuArray<int,AMREX_SPACEDIM> len{AMREX_D_DECL(m_geom.Domain().length(0),
                                                  m_geom.Domain().length(1),
                                                  m_geom.Domain().length(2))};

    auto f = [=] AMREX_GPU_DEVICE (int i, int j, int k, CellData < GpuComplex < T >> spectral_data)
    {
        amrex::ignore_unused(j,k);
        AMREX_D_TERM(T a = fac[0]*T(i);,
                     T b = fac[1]*T(j);,
                     T c = fac[2]*T(k));
        // Use exact cos(pi/2) = 0 so that null modes get a zero symbol.
        AMREX_D_TERM(T ca = (2*i == len[0]) ? T(0) : std::cos(a);,
                     T cb = (2*j == len[1]) ? T(0) : std::cos(b);,
                     T cc = (2*k == len[2]) ? T(0) : std::cos(c));

        // Fourier symbol of the node-centered cross-stencil divergence, up
        // to a common factor that cancels in the projection.
#if (AMREX_SPACEDIM == 2)
        T const sx = std::sin(a)*cb*dxfac[0];
        T const sy = ca*std::sin(b)*dxfac[1];

        T const s2 = sx*sx + sy*sy;
        // Select the denominator instead of guarding the division (see Poisson).
        T const s2d = (s2 != T(0)) ? s2 : T(1);
        auto cdiv = (sx * spectral_data[0] + sy * spectral_data[1]) / s2d;
        spectral_data[0] -= sx * cdiv;
        spectral_data[1] -= sy * cdiv;
#else
        T const sx = std::sin(a)*cb*cc*dxfac[0];
        T const sy = ca*std::sin(b)*cc*dxfac[1];
        T const sz = ca*cb*std::sin(c)*dxfac[2];

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

}

#endif
