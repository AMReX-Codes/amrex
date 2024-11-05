#ifndef AMREX_FFT_POISSON_H_
#define AMREX_FFT_POISSON_H_

#include <AMReX_FFT.H>
#include <AMReX_Geometry.H>

namespace amrex::FFT
{

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
        : m_geom(geom), m_bc(bc), m_r2x(geom.Domain(),bc)
    {}

    template <typename FA=MF, std::enable_if_t<IsFabArray_v<FA>,int> = 0>
    explicit Poisson (Geometry const& geom)
        : m_geom(geom),
          m_bc{AMREX_D_DECL(std::make_pair(Boundary::periodic,Boundary::periodic),
                            std::make_pair(Boundary::periodic,Boundary::periodic),
                            std::make_pair(Boundary::periodic,Boundary::periodic))},
          m_r2x(geom.Domain(),m_bc)
    {
        AMREX_ALWAYS_ASSERT(m_geom.isAllPeriodic());
    }

    void solve (MF& soln, MF const& rhs);

private:
    Geometry m_geom;
    Array<std::pair<Boundary,Boundary>,AMREX_SPACEDIM> m_bc;
    R2X<typename MF::value_type> m_r2x;
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
 * \brief 3D Poisson solver for periodic boundaries in the first two
 * dimensions and Neumann in the last dimension.
 */
template <typename MF = MultiFab>
class PoissonHybrid
{
public:

    template <typename FA=MF, std::enable_if_t<IsFabArray_v<FA>,int> = 0>
    explicit PoissonHybrid (Geometry const& geom)
        : m_geom(geom), m_r2c(geom.Domain(), Info().setBatchMode(true))
    {
#if (AMREX_SPACEDIM == 3)
        AMREX_ALWAYS_ASSERT(geom.isPeriodic(0) && geom.isPeriodic(1));
#else
        amrex::Abort("FFT::PoissonHybrid: 1D & 2D todo");
#endif
    }

    void solve (MF& soln, MF const& rhs);

private:
    Geometry m_geom;
    R2C<typename MF::value_type, Direction::both> m_r2c;
};

template <typename MF>
void Poisson<MF>::solve (MF& soln, MF const& rhs)
{
    BL_PROFILE("FFT::Poisson::solve");

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
    auto scale = m_r2x.scalingFactor();

    GpuArray<T,AMREX_SPACEDIM> offset{AMREX_D_DECL(T(0),T(0),T(0))};
    // Not sure about odd-even and even-odd yet
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

    m_r2x.forwardThenBackward(rhs, soln,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, auto& spectral_data)
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
    });
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

template <typename MF>
void PoissonHybrid<MF>::solve (MF& soln, MF const& rhs)
{
    BL_PROFILE("FFT::PoissonHybrid::solve");

#if (AMREX_SPACEDIM < 3)
    amrex::ignore_unused(soln, rhs);
#else
    using T = typename MF::value_type;

    auto facx = T(2)*Math::pi<T>()/T(m_geom.ProbLength(0));
    auto facy = T(2)*Math::pi<T>()/T(m_geom.ProbLength(1));
    auto dx = T(m_geom.CellSize(0));
    auto dy = T(m_geom.CellSize(1));
    auto scale = T(1.0)/(T(m_geom.Domain().length(0)) *
                         T(m_geom.Domain().length(1)));
    auto ny = m_geom.Domain().length(1);
    auto nz = m_geom.Domain().length(2);

    Gpu::DeviceVector<T> delzv(nz, T(m_geom.CellSize(2)));
    auto const* delz = delzv.data();

    Box cdomain = m_geom.Domain();
    cdomain.setBig(0,cdomain.length(0)/2);
    auto cba = amrex::decompose(cdomain, ParallelContext::NProcsSub(),
                                {AMREX_D_DECL(true,true,false)});
    DistributionMapping dm = detail::make_iota_distromap(cba.size());
    FabArray<BaseFab<GpuComplex<T> > > spmf(cba, dm, 1, 0);

    m_r2c.forward(rhs, spmf);

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
            T a = facx*i;
            T b = (j < ny/2) ? facy*j : facy*(ny-j);

            T k2 = T(2)*(std::cos(a*dx)-T(1))/(dx*dx)
                +  T(2)*(std::cos(b*dy)-T(1))/(dy*dy);

            // Tridiagonal solve with homogeneous Neumann
            for(int k=0; k < nz; k++) {
                if(k==0) {
                    ald(i,j,k) = 0.;
                    cud(i,j,k) = 2.0 /(delz[k]*(delz[k]+delz[k+1]));
                    bd(i,j,k) = k2 -ald(i,j,k)-cud(i,j,k);
                } else if (k == nz-1) {
                    ald(i,j,k) = 2.0 /(delz[k]*(delz[k]+delz[k-1]));
                    cud(i,j,k) = 0.;
                    bd(i,j,k) = k2 -ald(i,j,k)-cud(i,j,k);
                    if (i == 0 && j == 0) {
                        bd(i,j,k) *= 2.0;
                    }
                } else {
                    ald(i,j,k) = 2.0 /(delz[k]*(delz[k]+delz[k-1]));
                    cud(i,j,k) = 2.0 /(delz[k]*(delz[k]+delz[k+1]));
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

        Gpu::DeviceVector<GpuComplex<Real>> ald(nz);
        Gpu::DeviceVector<GpuComplex<Real>> bd(nz);
        Gpu::DeviceVector<GpuComplex<Real>> cud(nz);
        Gpu::DeviceVector<GpuComplex<Real>> scratch(nz);

        amrex::LoopOnCpu(xybox, [&] (int i, int j, int)
        {
            T a = facx*i;
            T b = (j < ny/2) ? facy*j : facy*(ny-j);

            T k2 = T(2)*(std::cos(a*dx)-T(1))/(dx*dx)
                +  T(2)*(std::cos(b*dy)-T(1))/(dy*dy);

            // Tridiagonal solve with homogeneous Neumann
            for(int k=0; k < nz; k++) {
                if(k==0) {
                    ald[k] = 0.;
                    cud[k] = 2.0 /(delz[k]*(delz[k]+delz[k+1]));
                    bd[k] = k2 -ald[k]-cud[k];
                } else if (k == nz-1) {
                    ald[k] = 2.0 /(delz[k]*(delz[k]+delz[k-1]));
                    cud[k] = 0.;
                    bd[k] = k2 -ald[k]-cud[k];
                    if (i == 0 && j == 0) {
                        bd[k] *= 2.0;
                    }
                } else {
                    ald[k] = 2.0 /(delz[k]*(delz[k]+delz[k-1]));
                    cud[k] = 2.0 /(delz[k]*(delz[k]+delz[k+1]));
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

    m_r2c.backward(spmf, soln);
#endif
}

}

#endif
