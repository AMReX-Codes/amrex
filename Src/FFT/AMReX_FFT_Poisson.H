#ifndef AMREX_FFT_POISSON_H_
#define AMREX_FFT_POISSON_H_

#include <AMReX_FFT.H>
#include <AMReX_Geometry.H>

namespace amrex::FFT
{

namespace detail {
template <typename MF>
void fill_physbc (MF& mf, Geometry const& geom,
                  Array<std::pair<Boundary,Boundary>,AMREX_SPACEDIM> const& bc);
}

/**
 * \brief Poisson solver for periodic, Dirichlet & Neumann boundaries using
 * FFT.
 */
template <typename MF = MultiFab>
class Poisson
{
public:

    template <typename FA=MF, std::enable_if_t<IsFabArray_v<FA>,int> = 0>
    Poisson (Geometry const& geom,
             Array<std::pair<Boundary,Boundary>,AMREX_SPACEDIM> const& bc)
        : m_geom(geom), m_bc(bc)
    {
        bool all_periodic = true;
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            all_periodic = all_periodic
                && (bc[idim].first == Boundary::periodic)
                && (bc[idim].second == Boundary::periodic);
        }
        if (all_periodic) {
            m_r2c = std::make_unique<R2C<typename MF::value_type>>(m_geom.Domain());
        } else {
            m_r2x = std::make_unique<R2X<typename MF::value_type>> (m_geom.Domain(), m_bc);
        }
    }

    template <typename FA=MF, std::enable_if_t<IsFabArray_v<FA>,int> = 0>
    explicit Poisson (Geometry const& geom)
        : m_geom(geom),
          m_bc{AMREX_D_DECL(std::make_pair(Boundary::periodic,Boundary::periodic),
                            std::make_pair(Boundary::periodic,Boundary::periodic),
                            std::make_pair(Boundary::periodic,Boundary::periodic))}
    {
        if (m_geom.isAllPeriodic()) {
            m_r2c = std::make_unique<R2C<typename MF::value_type>>(m_geom.Domain());
        } else {
            amrex::Abort("FFT::Poisson: wrong BC");
        }
    }

    /*
     * \brief Solve del dot grad soln = rhs
     *
     * If soln has ghost cells, one layer of ghost cells will be filled
     * except for the corners of the physical domain where they are not used
     * in the cross stencil of the operator. The two MFs could be the same MF.
     */
    void solve (MF& soln, MF const& rhs);

private:
    Geometry m_geom;
    Array<std::pair<Boundary,Boundary>,AMREX_SPACEDIM> m_bc;
    std::unique_ptr<R2X<typename MF::value_type>> m_r2x;
    std::unique_ptr<R2C<typename MF::value_type>> m_r2c;
};

#if (AMREX_SPACEDIM == 3)
/**
 * \brief Poisson solve for Open BC using FFT.
 */
template <typename MF = MultiFab>
class PoissonOpenBC
{
public:

    template <typename FA=MF, std::enable_if_t<IsFabArray_v<FA>,int> = 0>
    explicit PoissonOpenBC (Geometry const& geom,
                            IndexType ixtype = IndexType::TheCellType(),
                            IntVect const& ngrow = IntVect(0));

    void solve (MF& soln, MF const& rhs);

    void define_doit (); // has to be public for cuda

private:
    Geometry m_geom;
    Box m_grown_domain;
    IntVect m_ngrow;
    OpenBCSolver<typename MF::value_type> m_solver;
};
#endif

/**
 * \brief 3D Poisson solver for periodic, Dirichlet & Neumann boundaries in
 * the first two dimensions, and Neumann in the last dimension. The last
 * dimension could have non-uniform mesh.
 */
template <typename MF = MultiFab>
class PoissonHybrid
{
public:
    using T = typename MF::value_type;

    template <typename FA=MF, std::enable_if_t<IsFabArray_v<FA>,int> = 0>
    PoissonHybrid (Geometry const& geom,
                   Array<std::pair<Boundary,Boundary>,AMREX_SPACEDIM> const& bc)
        : m_geom(geom), m_bc(bc)
    {
#if (AMREX_SPACEDIM < 3)
        amrex::Abort("FFT::PoissonHybrid: 1D & 2D todo");
        return;
#endif
        bool periodic_xy = true;
        for (int idim = 0; idim < 2; ++idim) {
            if (m_geom.Domain().length(idim) > 1) {
                periodic_xy = periodic_xy && (bc[idim].first == Boundary::periodic);
                AMREX_ALWAYS_ASSERT((bc[idim].first == Boundary::periodic &&
                                     bc[idim].second == Boundary::periodic) ||
                                    (bc[idim].first != Boundary::periodic &&
                                     bc[idim].second != Boundary::periodic));
            }
        }
        Info info{};
        info.setBatchMode(true);
        if (periodic_xy) {
            m_r2c = std::make_unique<R2C<typename MF::value_type>>(m_geom.Domain(),
                                                                   info);
        } else {
            m_r2x = std::make_unique<R2X<typename MF::value_type>> (m_geom.Domain(),
                                                                    m_bc, info);
        }
        build_spmf();
    }

    template <typename FA=MF, std::enable_if_t<IsFabArray_v<FA>,int> = 0>
    explicit PoissonHybrid (Geometry const& geom)
        : m_geom(geom),
          m_bc{AMREX_D_DECL(std::make_pair(Boundary::periodic,Boundary::periodic),
                            std::make_pair(Boundary::periodic,Boundary::periodic),
                            std::make_pair(Boundary::even,Boundary::even))},
          m_r2c(std::make_unique<R2C<typename MF::value_type>>
                (geom.Domain(), Info().setBatchMode(true)))
    {
#if (AMREX_SPACEDIM == 3)
        AMREX_ALWAYS_ASSERT(geom.isPeriodic(0) && geom.isPeriodic(1));
#else
        amrex::Abort("FFT::PoissonHybrid: 1D & 2D todo");
        return;
#endif
        build_spmf();
    }

    /*
     * \brief Solve del dot grad soln = rhs
     *
     * If soln has ghost cells, one layer of ghost cells will be filled
     * except for the corners of the physical domain where they are not used
     * in the cross stencil of the operator. The two MFs could be the same MF.
     */
    void solve (MF& soln, MF const& rhs);
    void solve (MF& soln, MF const& rhs, Vector<T> const& dz);
    void solve (MF& soln, MF const& rhs, Gpu::DeviceVector<T> const& dz);

    template <typename TRIA, typename TRIC>
    void solve (MF& soln, MF const& rhs, TRIA const& tria, TRIC const& tric);

    // This is public for cuda
    template <typename FA, typename TRIA, typename TRIC>
    void solve_z (FA& spmf, TRIA const& tria, TRIC const& tric);

    [[nodiscard]] std::pair<BoxArray,DistributionMapping> getSpectralDataLayout () const;

private:

    void build_spmf ();

    Geometry m_geom;
    Array<std::pair<Boundary,Boundary>,AMREX_SPACEDIM> m_bc;
    std::unique_ptr<R2X<typename MF::value_type>> m_r2x;
    std::unique_ptr<R2C<typename MF::value_type>> m_r2c;
    MF m_spmf_r;
    using cMF = FabArray<BaseFab<GpuComplex<T>>>;
    cMF m_spmf_c;
};

template <typename MF>
void Poisson<MF>::solve (MF& soln, MF const& rhs)
{
    BL_PROFILE("FFT::Poisson::solve");

    AMREX_ASSERT(soln.is_cell_centered() && rhs.is_cell_centered());

    using T = typename MF::value_type;

    GpuArray<T,AMREX_SPACEDIM> fac
        {AMREX_D_DECL(Math::pi<T>()/T(m_geom.Domain().length(0)),
                      Math::pi<T>()/T(m_geom.Domain().length(1)),
                      Math::pi<T>()/T(m_geom.Domain().length(2)))};
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (m_bc[idim].first == Boundary::periodic) {
            fac[idim] *= T(2);
        }
    }
    GpuArray<T,AMREX_SPACEDIM> dxfac
        {AMREX_D_DECL(T(2)/T(m_geom.CellSize(0)*m_geom.CellSize(0)),
                      T(2)/T(m_geom.CellSize(1)*m_geom.CellSize(1)),
                      T(2)/T(m_geom.CellSize(2)*m_geom.CellSize(2)))};
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

    auto f = [=] AMREX_GPU_DEVICE (int i, int j, int k, auto& spectral_data)
    {
        amrex::ignore_unused(j,k);
        AMREX_D_TERM(T a = fac[0]*(i+offset[0]);,
                     T b = fac[1]*(j+offset[1]);,
                     T c = fac[2]*(k+offset[2]));
        T k2 = AMREX_D_TERM(dxfac[0]*(std::cos(a)-T(1)),
                           +dxfac[1]*(std::cos(b)-T(1)),
                           +dxfac[2]*(std::cos(c)-T(1)));
        if (k2 != T(0)) {
            spectral_data /= k2;
        }
        spectral_data *= scale;
    };

    IntVect const& ng = amrex::elemwiseMin(soln.nGrowVect(), IntVect(1));

    if (m_r2x) {
        m_r2x->template forwardThenBackward_doit<0>(rhs, soln, f, ng, m_geom.periodicity());
        detail::fill_physbc(soln, m_geom, m_bc);
    } else {
        m_r2c->forward(rhs);
        m_r2c->template post_forward_doit<0>(f);
        m_r2c->backward_doit(soln, ng, m_geom.periodicity());
    }
}

#if (AMREX_SPACEDIM == 3)

template <typename MF>
template <typename FA, std::enable_if_t<IsFabArray_v<FA>,int> FOO>
PoissonOpenBC<MF>::PoissonOpenBC (Geometry const& geom, IndexType ixtype,
                                  IntVect const& ngrow)
    : m_geom(geom),
      m_grown_domain(amrex::grow(amrex::convert(geom.Domain(),ixtype),ngrow)),
      m_ngrow(ngrow),
      m_solver(m_grown_domain)
{
    define_doit();
}

template <typename MF>
void PoissonOpenBC<MF>::define_doit ()
{
    using T = typename MF::value_type;
    auto const& lo = m_grown_domain.smallEnd();
    auto const dx = T(m_geom.CellSize(0));
    auto const dy = T(m_geom.CellSize(1));
    auto const dz = T(m_geom.CellSize(2));
    auto const gfac = T(1)/T(std::sqrt(T(12)));
    // 0.125 comes from that there are 8 Gauss quadrature points
    auto const fac = T(-0.125) * (dx*dy*dz) / (T(4)*Math::pi<T>());
    m_solver.setGreensFunction([=] AMREX_GPU_DEVICE (int i, int j, int k) -> T
    {
        auto x = (T(i-lo[0]) - gfac) * dx; // first Gauss quadrature point
        auto y = (T(j-lo[1]) - gfac) * dy;
        auto z = (T(k-lo[2]) - gfac) * dz;
        T r = 0;
        for (int gx = 0; gx < 2; ++gx) {
        for (int gy = 0; gy < 2; ++gy) {
        for (int gz = 0; gz < 2; ++gz) {
            auto xg = x + 2*gx*gfac*dx;
            auto yg = y + 2*gy*gfac*dy;
            auto zg = z + 2*gz*gfac*dz;
            r += T(1)/std::sqrt(xg*xg+yg*yg+zg*zg);
        }}}
        return fac * r;
    });
}

template <typename MF>
void PoissonOpenBC<MF>::solve (MF& soln, MF const& rhs)
{
    AMREX_ASSERT(m_grown_domain.ixType() == soln.ixType() && m_grown_domain.ixType() == rhs.ixType());
    m_solver.solve(soln, rhs);
}

#endif /* AMREX_SPACEDIM == 3 */

namespace fft_poisson_detail {
    template <typename T>
    struct Tri_Uniform {
        [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        T operator() (int, int, int) const
        {
            return m_dz2inv;
        }
        T m_dz2inv;
    };

    template <typename T>
    struct TriA {
        [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        T operator() (int, int, int k) const
        {
            return T(2.0) /(m_dz[k]*(m_dz[k]+m_dz[k-1]));
        }
        T const* m_dz;
    };

    template <typename T>
    struct TriC {
        [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        T operator() (int, int, int k) const
        {
            return T(2.0) /(m_dz[k]*(m_dz[k]+m_dz[k+1]));
        }
        T const* m_dz;
    };
}

template <typename MF>
std::pair<BoxArray,DistributionMapping>
PoissonHybrid<MF>::getSpectralDataLayout () const
{
    if (!m_spmf_r.empty()) {
        return std::make_pair(m_spmf_r.boxArray(), m_spmf_r.DistributionMap());
    } else {
        return std::make_pair(m_spmf_c.boxArray(), m_spmf_c.DistributionMap());
    }
}

template <typename MF>
void PoissonHybrid<MF>::build_spmf ()
{
#if (AMREX_SPACEDIM == 3)
    AMREX_ALWAYS_ASSERT(m_geom.Domain().length(2) > 1 &&
                        (m_geom.Domain().length(0) > 1 ||
                         m_geom.Domain().length(1) > 1));

    if (m_r2c) {
        Box cdomain = m_geom.Domain();
        if (cdomain.length(0) > 1) {
            cdomain.setBig(0,cdomain.length(0)/2);
        } else {
            cdomain.setBig(1,cdomain.length(1)/2);
        }
        auto cba = amrex::decompose(cdomain, ParallelContext::NProcsSub(),
                                    {AMREX_D_DECL(true,true,false)});
        DistributionMapping dm = detail::make_iota_distromap(cba.size());
        m_spmf_c.define(cba, dm, 1, 0);
    } else if (m_geom.Domain().length(0) > 1 &&
               m_geom.Domain().length(1) > 1) {
        if (m_r2x->m_cy.empty()) { // spectral data is real
            auto sba = amrex::decompose(m_geom.Domain(),ParallelContext::NProcsSub(),
                                        {AMREX_D_DECL(true,true,false)});
            DistributionMapping dm = detail::make_iota_distromap(sba.size());
            m_spmf_r.define(sba, dm, 1, 0);
        } else { // spectral data is complex. one of the first two dimensions is periodic.
            Box cdomain = m_geom.Domain();
            if (m_bc[0].first == Boundary::periodic) {
                cdomain.setBig(0,cdomain.length(0)/2);
            } else {
                cdomain.setBig(1,cdomain.length(1)/2);
            }
            auto cba = amrex::decompose(cdomain, ParallelContext::NProcsSub(),
                                        {AMREX_D_DECL(true,true,false)});
            DistributionMapping dm = detail::make_iota_distromap(cba.size());
            m_spmf_c.define(cba, dm, 1, 0);
        }
    } else {
        // spectral data is real
        auto sba = amrex::decompose(m_geom.Domain(),ParallelContext::NProcsSub(),
                                    {AMREX_D_DECL(true,true,false)});
        DistributionMapping dm = detail::make_iota_distromap(sba.size());
        m_spmf_r.define(sba, dm, 1, 0);
    }
#else
    amrex::ignore_unused(this);
#endif
}

template <typename MF>
void PoissonHybrid<MF>::solve (MF& soln, MF const& rhs)
{
    auto delz = T(m_geom.CellSize(AMREX_SPACEDIM-1));
    solve(soln, rhs,
          fft_poisson_detail::Tri_Uniform<T>{T(1)/(delz*delz)},
          fft_poisson_detail::Tri_Uniform<T>{T(1)/(delz*delz)});
}

template <typename MF>
void PoissonHybrid<MF>::solve (MF& soln, MF const& rhs, Gpu::DeviceVector<T> const& dz)
{
    auto const* pdz = dz.dataPtr();
    solve(soln, rhs,
          fft_poisson_detail::TriA<T>{pdz},
          fft_poisson_detail::TriC<T>{pdz});
}

template <typename MF>
void PoissonHybrid<MF>::solve (MF& soln, MF const& rhs, Vector<T> const& dz)
{
    AMREX_ASSERT(soln.is_cell_centered() && rhs.is_cell_centered());

#ifdef AMREX_USE_GPU
    Gpu::DeviceVector<T> d_dz(dz.size());
    Gpu::htod_memcpy_async(d_dz.data(), dz.data(), dz.size()*sizeof(T));
    auto const* pdz = d_dz.data();
#else
    auto const* pdz = dz.data();
#endif
    solve(soln, rhs,
          fft_poisson_detail::TriA<T>{pdz},
          fft_poisson_detail::TriC<T>{pdz});
}

template <typename MF>
template <typename TRIA, typename TRIC>
void PoissonHybrid<MF>::solve (MF& soln, MF const& rhs, TRIA const& tria,
                               TRIC const& tric)
{
    BL_PROFILE("FFT::PoissonHybrid::solve");

    AMREX_ASSERT(soln.is_cell_centered() && rhs.is_cell_centered());

#if (AMREX_SPACEDIM < 3)
    amrex::ignore_unused(soln, rhs, tria, tric);
#else

    IntVect const& ng = amrex::elemwiseMin(soln.nGrowVect(), IntVect(1));

    if (m_r2c)
    {
        m_r2c->forward(rhs, m_spmf_c);
        solve_z(m_spmf_c, tria, tric);
        m_r2c->backward_doit(m_spmf_c, soln, ng, m_geom.periodicity());
    }
    else
    {
        if (m_r2x->m_cy.empty()) { // spectral data is real
            m_r2x->forward(rhs, m_spmf_r);
            solve_z(m_spmf_r, tria, tric);
            m_r2x->backward(m_spmf_r, soln, ng, m_geom.periodicity());
        } else { // spectral data is complex.
            m_r2x->forward(rhs, m_spmf_c);
            solve_z(m_spmf_c, tria, tric);
            m_r2x->backward(m_spmf_c, soln, ng, m_geom.periodicity());
        }
    }

    detail::fill_physbc(soln, m_geom, m_bc);
#endif
}

template <typename MF>
template <typename FA, typename TRIA, typename TRIC>
void PoissonHybrid<MF>::solve_z (FA& spmf, TRIA const& tria, TRIC const& tric)
{
    BL_PROFILE("PoissonHybrid::solve_z");

#if (AMREX_SPACEDIM < 3)
    amrex::ignore_unused(spmf, tria, tric);
#else
    auto facx = Math::pi<T>()/T(m_geom.Domain().length(0));
    auto facy = Math::pi<T>()/T(m_geom.Domain().length(1));
    if (m_bc[0].first == Boundary::periodic) { facx *= T(2); }
    if (m_bc[1].first == Boundary::periodic) { facy *= T(2); }
    auto dxfac = T(2)/T(m_geom.CellSize(0)*m_geom.CellSize(0));
    auto dyfac = T(2)/T(m_geom.CellSize(1)*m_geom.CellSize(1));
    auto scale = (m_r2x) ? m_r2x->scalingFactor() : m_r2c->scalingFactor();

    if (m_geom.Domain().length(0) == 1) { dxfac = 0; }
    if (m_geom.Domain().length(1) == 1) { dyfac = 0; }

    GpuArray<T,AMREX_SPACEDIM-1> offset{T(0),T(0)};
    for (int idim = 0; idim < AMREX_SPACEDIM-1; ++idim) {
        if (m_geom.Domain().length(idim) > 1) {
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
    }

    bool has_dirichlet = (offset[0] != T(0)) || (offset[1] != T(0));

    auto nz = m_geom.Domain().length(2);

    for (MFIter mfi(spmf); mfi.isValid(); ++mfi)
    {
        auto const& spectral = spmf.array(mfi);
        auto const& box = mfi.validbox();
        auto const& xybox = amrex::makeSlab(box, 2, 0);

#ifdef AMREX_USE_GPU
        // xxxxx TODO: We need to explore how to optimize this
        // function. Maybe we can use cusparse. Maybe we should make
        // z-direction to be the unit stride direction.

        FArrayBox tridiag_workspace(box,4);
        auto const& ald = tridiag_workspace.array(0);
        auto const& bd = tridiag_workspace.array(1);
        auto const& cud = tridiag_workspace.array(2);
        auto const& scratch = tridiag_workspace.array(3);

        amrex::ParallelFor(xybox, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            T a = facx*(i+offset[0]);
            T b = facy*(j+offset[1]);
            T k2 = dxfac * (std::cos(a)-T(1))
                +  dyfac * (std::cos(b)-T(1));

            // Tridiagonal solve with homogeneous Neumann
            for(int k=0; k < nz; k++) {
                if(k==0) {
                    ald(i,j,k) = T(0.);
                    cud(i,j,k) = tric(i,j,k);
                    bd(i,j,k) = k2 -ald(i,j,k)-cud(i,j,k);
                } else if (k == nz-1) {
                    ald(i,j,k) = tria(i,j,k);
                    cud(i,j,k) = T(0.);
                    bd(i,j,k) = k2 -ald(i,j,k)-cud(i,j,k);
                    if (i == 0 && j == 0 && !has_dirichlet) {
                        bd(i,j,k) *= T(2.0);
                    }
                } else {
                    ald(i,j,k) = tria(i,j,k);
                    cud(i,j,k) = tric(i,j,k);
                    bd(i,j,k) = k2 -ald(i,j,k)-cud(i,j,k);
                }
            }

            scratch(i,j,0) = cud(i,j,0)/bd(i,j,0);
            spectral(i,j,0) = spectral(i,j,0)/bd(i,j,0);

            for (int k = 1; k < nz; k++) {
                if (k < nz-1) {
                    scratch(i,j,k) = cud(i,j,k) / (bd(i,j,k) - ald(i,j,k) * scratch(i,j,k-1));
                }
                spectral(i,j,k) = (spectral(i,j,k) - ald(i,j,k) * spectral(i,j,k - 1))
                    / (bd(i,j,k) - ald(i,j,k) * scratch(i,j,k-1));
            }

            for (int k = nz - 2; k >= 0; k--) {
                spectral(i,j,k) -= scratch(i,j,k) * spectral(i,j,k + 1);
            }

            for (int k = 0; k < nz; ++k) {
                spectral(i,j,k) *= scale;
            }
        });
        Gpu::streamSynchronize();

#else

        Gpu::DeviceVector<T> ald(nz);
        Gpu::DeviceVector<T> bd(nz);
        Gpu::DeviceVector<T> cud(nz);
        Gpu::DeviceVector<T> scratch(nz);

        amrex::LoopOnCpu(xybox, [&] (int i, int j, int)
        {
            T a = facx*(i+offset[0]);
            T b = facy*(j+offset[1]);
            T k2 = dxfac * (std::cos(a)-T(1))
                +  dyfac * (std::cos(b)-T(1));

            // Tridiagonal solve with homogeneous Neumann
            for(int k=0; k < nz; k++) {
                if(k==0) {
                    ald[k] = T(0.);
                    cud[k] = tric(i,j,k);
                    bd[k] = k2 -ald[k]-cud[k];
                } else if (k == nz-1) {
                    ald[k] = tria(i,j,k);
                    cud[k] = T(0.);
                    bd[k] = k2 -ald[k]-cud[k];
                    if (i == 0 && j == 0 && !has_dirichlet) {
                        bd[k] *= T(2.0);
                    }
                } else {
                    ald[k] = tria(i,j,k);
                    cud[k] = tric(i,j,k);
                    bd[k] = k2 -ald[k]-cud[k];
                }
            }

            scratch[0] = cud[0]/bd[0];
            spectral(i,j,0) = spectral(i,j,0)/bd[0];

            for (int k = 1; k < nz; k++) {
                if (k < nz-1) {
                    scratch[k] = cud[k] / (bd[k] - ald[k] * scratch[k-1]);
                }
                spectral(i,j,k) = (spectral(i,j,k) - ald[k] * spectral(i,j,k - 1))
                    / (bd[k] - ald[k] * scratch[k-1]);
            }

            for (int k = nz - 2; k >= 0; k--) {
                spectral(i,j,k) -= scratch[k] * spectral(i,j,k + 1);
            }

            for (int k = 0; k < nz; ++k) {
                spectral(i,j,k) *= scale;
            }
        });
#endif
    }
#endif
}

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
                tags.push_back({arr, b, fbc, face});
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

}

#endif
