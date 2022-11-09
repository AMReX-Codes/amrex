#ifndef INIT_PROB_K_H_
#define INIT_PROB_K_H_

#include <AMReX_FArrayBox.H>

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void actual_init_poisson (int i, int j, int k,
                          amrex::Array4<amrex::Real> const& rhs,
                          amrex::Array4<amrex::Real> const& exact,
                          amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& prob_lo,
                          amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx)
{
    constexpr amrex::Real tpi = 2.*3.1415926535897932;
    constexpr amrex::Real fpi = 4.*3.1415926535897932;
    constexpr amrex::Real fac = tpi*tpi*static_cast<amrex::Real>(AMREX_SPACEDIM);
    amrex::Real x = prob_lo[0] + dx[0] * (i + 0.5);
    amrex::Real y = prob_lo[1] + dx[1] * (j + 0.5);

#if (AMREX_SPACEDIM == 2)

    exact(i,j,k) = (std::sin(tpi*x) * std::sin(tpi*y))
           + .25 * (std::sin(fpi*x) * std::sin(fpi*y));

    rhs(i,j,k) = -fac * (std::sin(tpi*x) * std::sin(tpi*y))
                 -fac * (std::sin(fpi*x) * std::sin(fpi*y));

#else

    amrex::Real z = prob_lo[2] + dx[2] * (k + 0.5);

    exact(i,j,k) = (std::sin(tpi*x) * std::sin(tpi*y) * std::sin(tpi*z))
           + .25 * (std::sin(fpi*x) * std::sin(fpi*y) * std::sin(fpi*z));

    rhs(i,j,k) = -fac * (std::sin(tpi*x) * std::sin(tpi*y) * std::sin(tpi*z))
                 -fac * (std::sin(fpi*x) * std::sin(fpi*y) * std::sin(fpi*z));

#endif
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void actual_init_bcoef (int i, int j, int k,
                        amrex::Array4<amrex::Real> const& beta,
                        amrex::Array4<amrex::Real> const& sol,
                        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& prob_lo,
                        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& prob_hi,
                        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
                        amrex::Box const& domain)
{
    constexpr amrex::Real w = 0.05;
    constexpr amrex::Real sigma = 10.;
    const amrex::Real theta = 0.5*std::log(3.) / (w + 1.e-50);

    constexpr amrex::Real pi = 3.1415926535897932;
    constexpr amrex::Real tpi =  2.*pi;
    constexpr amrex::Real fpi =  4.*pi;

    amrex::Real xc = (prob_hi[0] + prob_lo[0])*0.5;
    amrex::Real yc = (prob_hi[1] + prob_lo[1])*0.5;
#if (AMREX_SPACEDIM == 2)
    amrex::Real zc = 0.0;
#else
    amrex::Real zc = (prob_hi[2] + prob_lo[2])*0.5;
#endif

    amrex::Real x = prob_lo[0] + dx[0] * (i + 0.5);
    amrex::Real y = prob_lo[1] + dx[1] * (j + 0.5);
#if (AMREX_SPACEDIM == 2)
    amrex::Real z = 0.0;
#else
    amrex::Real z = prob_lo[2] + dx[2] * (k + 0.5);
#endif

    amrex::Real r = std::sqrt((x-xc)*(x-xc) + (y-yc)*(y-yc) + (z-zc)*(z-zc));
    beta(i,j,k) = (sigma-1.)/2.*std::tanh(theta*(r-0.25)) + (sigma+1.)/2.;

    if (domain.contains(i,j,k)) {
        sol(i,j,k) = 0.0;
    } else {
        AMREX_D_TERM(x = std::clamp(x, prob_lo[0], prob_hi[0]);,
                     y = std::clamp(y, prob_lo[1], prob_hi[1]);,
                     z = std::clamp(z, prob_lo[2], prob_hi[2]););
        sol(i,j,k) = std::cos(tpi*x) * std::cos(tpi*y) * std::cos(tpi*z)
            +  .25 * std::cos(fpi*x) * std::cos(fpi*y) * std::cos(fpi*z);
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void actual_init_abeclap (int i, int j, int k,
                          amrex::Array4<amrex::Real      > const& rhs,
                          amrex::Array4<amrex::Real      > const& exact,
                          amrex::Array4<amrex::Real      > const& alpha,
                          amrex::Array4<amrex::Real const> const& beta,
                          amrex::Real a, amrex::Real b,
                          amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& prob_lo,
                          amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& prob_hi,
                          amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx)
{
    constexpr amrex::Real w = 0.05;
    constexpr amrex::Real sigma = 10.;
    const amrex::Real theta = 0.5*std::log(3.) / (w + 1.e-50);

    constexpr amrex::Real pi = 3.1415926535897932;
    constexpr amrex::Real tpi =  2.*pi;
    constexpr amrex::Real fpi =  4.*pi;
    constexpr amrex::Real fac = static_cast<amrex::Real>(AMREX_SPACEDIM*4)*pi*pi;

    amrex::Real xc = (prob_hi[0] + prob_lo[0])*0.5;
    amrex::Real yc = (prob_hi[1] + prob_lo[1])*0.5;
#if (AMREX_SPACEDIM == 2)
    amrex::Real zc = 0.0;
#else
    amrex::Real zc = (prob_hi[2] + prob_lo[2])*0.5;
#endif

    amrex::Real x = prob_lo[0] + dx[0] * (i + 0.5);
    amrex::Real y = prob_lo[1] + dx[1] * (j + 0.5);
#if (AMREX_SPACEDIM == 2)
    amrex::Real z = 0.0;
#else
    amrex::Real z = prob_lo[2] + dx[2] * (k + 0.5);
#endif

    amrex::Real r = std::sqrt((x-xc)*(x-xc) + (y-yc)*(y-yc) + (z-zc)*(z-zc));
    amrex::Real tmp = std::cosh(theta*(r-0.25));
    amrex::Real dbdrfac = (sigma-1.)/2./(tmp*tmp) * theta/r;
    dbdrfac *= b;

    alpha(i,j,k) = 1.;

    exact(i,j,k) = std::cos(tpi*x) * std::cos(tpi*y) * std::cos(tpi*z)
           + .25 * std::cos(fpi*x) * std::cos(fpi*y) * std::cos(fpi*z);

    rhs(i,j,k) = beta(i,j,k)*b*fac*(std::cos(tpi*x) * std::cos(tpi*y) * std::cos(tpi*z)
                                  + std::cos(fpi*x) * std::cos(fpi*y) * std::cos(fpi*z))
             + dbdrfac*((x-xc)*(tpi*std::sin(tpi*x) * std::cos(tpi*y) * std::cos(tpi*z)
                               + pi*std::sin(fpi*x) * std::cos(fpi*y) * std::cos(fpi*z))
                      + (y-yc)*(tpi*std::cos(tpi*x) * std::sin(tpi*y) * std::cos(tpi*z)
                               + pi*std::cos(fpi*x) * std::sin(fpi*y) * std::cos(fpi*z))
                      + (z-zc)*(tpi*std::cos(tpi*x) * std::cos(tpi*y) * std::sin(tpi*z)
                               + pi*std::cos(fpi*x) * std::cos(fpi*y) * std::sin(fpi*z)))
                             + a * (std::cos(tpi*x) * std::cos(tpi*y) * std::cos(tpi*z)
                           + 0.25 * std::cos(fpi*x) * std::cos(fpi*y) * std::cos(fpi*z));
}

#endif
