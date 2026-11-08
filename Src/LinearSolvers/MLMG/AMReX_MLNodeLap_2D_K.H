#ifndef AMREX_MLNODELAP_2D_K_H_
#define AMREX_MLNODELAP_2D_K_H_
#include <AMReX_Config.H>

namespace amrex {

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_zero_fine (int i, int j, int, Array4<T> const& phi,
                        Array4<int const> const& msk, int fine_flag) noexcept
{
    // Testing if the node is covered by a fine level in computing
    // coarse sync residual
    if (msk(i-1,j-1,0) == fine_flag &&
        msk(i  ,j-1,0) == fine_flag &&
        msk(i-1,j  ,0) == fine_flag &&
        msk(i  ,j  ,0) == fine_flag)
    {
        phi(i,j,0) = T(0.0);
    }
}

//
// coeffs
//

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_avgdown_coeff_x (int i, int j, int k, Array4<T> const& crse,
                              Array4<std::type_identity_t<T> const> const& fine) noexcept
{
    T a = fine(2*i  ,2*j,k) + fine(2*i  ,2*j+1,k);
    T b = fine(2*i+1,2*j,k) + fine(2*i+1,2*j+1,k);
    crse(i,j,k) = a*b/(a+b);
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_avgdown_coeff_y (int i, int j, int k, Array4<T> const& crse,
                              Array4<std::type_identity_t<T> const> const& fine) noexcept
{
    T a = fine(2*i,2*j  ,k) + fine(2*i+1,2*j  ,k);
    T b = fine(2*i,2*j+1,k) + fine(2*i+1,2*j+1,k);
    crse(i,j,k) = a*b/(a+b);
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_semi_avgdown_coeff (int i, int j, int k, Array4<T> const& crse,
                              Array4<std::type_identity_t<T> const> const& fine, int idir) noexcept
{
    if (idir == 1) {
        T a = fine(2*i  ,j,k);
        T b = fine(2*i+1,j,k);
        crse(i,j,k) = T(2.0)*a*b/(a+b);
    } else {
        T a = fine(i,2*j  ,k);
        T b = fine(i,2*j+1,k);
        crse(i,j,k) = T(2.0)*a*b/(a+b);
    }
}

//
// operator
//

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T mlndlap_adotx_ha (int i, int j, int k, Array4<T const> const& x,
                       Array4<std::type_identity_t<T> const> const& sx, Array4<std::type_identity_t<T> const> const& sy,
                       Array4<int const> const& msk, bool is_rz,
                       GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    if (msk(i,j,k)) {
        return T(0.0);
    } else {
        T facx = T(1./6.)*dxinv[0]*dxinv[0];
        T facy = T(1./6.)*dxinv[1]*dxinv[1];
        T y   = x(i-1,j-1,k)*(facx*sx(i-1,j-1,k)+facy*sy(i-1,j-1,k))
               +   x(i+1,j-1,k)*(facx*sx(i  ,j-1,k)+facy*sy(i  ,j-1,k))
               +   x(i-1,j+1,k)*(facx*sx(i-1,j  ,k)+facy*sy(i-1,j  ,k))
               +   x(i+1,j+1,k)*(facx*sx(i  ,j  ,k)+facy*sy(i  ,j  ,k))
               +   x(i-1,j,k)*(T(2.0)*facx*(sx(i-1,j-1,k)+sx(i-1,j,k))
                                   -     facy*(sy(i-1,j-1,k)+sy(i-1,j,k)))
               +   x(i+1,j,k)*(T(2.0)*facx*(sx(i  ,j-1,k)+sx(i  ,j,k))
                                   -     facy*(sy(i  ,j-1,k)+sy(i  ,j,k)))
               +   x(i,j-1,k)*(   -facx*(sx(i-1,j-1,k)+sx(i,j-1,k))
                        +T(2.0)*facy*(sy(i-1,j-1,k)+sy(i,j-1,k)))
               +   x(i,j+1,k)*(   -facx*(sx(i-1,j  ,k)+sx(i,j  ,k))
                        +T(2.0)*facy*(sy(i-1,j  ,k)+sy(i,j  ,k)))
               +   x(i,j,k)*(-T(2.0))*(facx*(sx(i-1,j-1,k)+sx(i,j-1,k)+sx(i-1,j,k)+sx(i,j,k))
                                         +facy*(sy(i-1,j-1,k)+sy(i,j-1,k)+sy(i-1,j,k)+sy(i,j,k)));
        if (is_rz) {
            T fp = facy / static_cast<T>(2*i+1);
            T fm = facy / static_cast<T>(2*i-1);
            y += (fm*sy(i-1,j  ,k)-fp*sy(i,j  ,k))*(x(i,j+1,k)-x(i,j,k))
               + (fm*sy(i-1,j-1,k)-fp*sy(i,j-1,k))*(x(i,j-1,k)-x(i,j,k));
        }
        return y;
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T mlndlap_adotx_aa (int i, int j, int k, Array4<T const> const& x,
                       Array4<std::type_identity_t<T> const> const& sig, Array4<int const> const& msk,
                       bool is_rz, GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    if (msk(i,j,k)) {
        return T(0.0);
    } else {
        T facx = T(1.0/6.0)*dxinv[0]*dxinv[0];
        T facy = T(1.0/6.0)*dxinv[1]*dxinv[1];
        T fxy = facx + facy;
        T f2xmy = T(2.0)*facx - facy;
        T fmx2y = T(2.0)*facy - facx;
        T y   = x(i-1,j-1,k)*fxy*sig(i-1,j-1,k)
               +   x(i+1,j-1,k)*fxy*sig(i  ,j-1,k)
               +   x(i-1,j+1,k)*fxy*sig(i-1,j  ,k)
               +   x(i+1,j+1,k)*fxy*sig(i  ,j  ,k)
               +   x(i-1,j,k)*f2xmy*(sig(i-1,j-1,k)+sig(i-1,j,k))
               +   x(i+1,j,k)*f2xmy*(sig(i  ,j-1,k)+sig(i  ,j,k))
               +   x(i,j-1,k)*fmx2y*(sig(i-1,j-1,k)+sig(i,j-1,k))
               +   x(i,j+1,k)*fmx2y*(sig(i-1,j  ,k)+sig(i,j  ,k))
               +   x(i,j,k)*(-T(2.0))*fxy*(sig(i-1,j-1,k)+sig(i,j-1,k)
                                             +sig(i-1,j,k)+sig(i,j,k));
        if (is_rz) {
            T fp = facy / static_cast<T>(2*i+1);
            T fm = facy / static_cast<T>(2*i-1);
            y += (fm*sig(i-1,j  ,k)-fp*sig(i,j  ,k))*(x(i,j+1,k)-x(i,j,k))
               + (fm*sig(i-1,j-1,k)-fp*sig(i,j-1,k))*(x(i,j-1,k)-x(i,j,k));
        }
        return y;
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T mlndlap_adotx_c (int i, int j, int k, Array4<T const> const& x,
                      std::type_identity_t<T> sigma, Array4<int const> const& msk,
                      bool is_rz, GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    if (msk(i,j,k)) {
        return T(0.0);
    } else {
        T facx = T(1.0/6.0)*dxinv[0]*dxinv[0];
        T facy = T(1.0/6.0)*dxinv[1]*dxinv[1];
        T fxy = facx + facy;
        T f2xmy = T(2.0)*facx - facy;
        T fmx2y = T(2.0)*facy - facx;
        T y  = (x(i-1,j-1,k)*fxy
               +   x(i+1,j-1,k)*fxy
               +   x(i-1,j+1,k)*fxy
               +   x(i+1,j+1,k)*fxy
               +   x(i-1,j,k)*f2xmy*T(2.)
               +   x(i+1,j,k)*f2xmy*T(2.)
               +   x(i,j-1,k)*fmx2y*T(2.)
               +   x(i,j+1,k)*fmx2y*T(2.)
               +   x(i,j,k)*(-T(2.0))*fxy*T(4.));
        if (is_rz) {
            T fp = facy / static_cast<T>(2*i+1);
            T fm = facy / static_cast<T>(2*i-1);
            y += ((fm-fp)*(x(i,j+1,k)-x(i,j,k))
               +  (fm-fp)*(x(i,j-1,k)-x(i,j,k)));
        }
        return y * sigma;
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_normalize_ha (int i, int j, int k, Array4<T> const& x, Array4<std::type_identity_t<T> const> const& sx,
                           Array4<std::type_identity_t<T> const> const& sy, Array4<int const> const& msk,
                           GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    T facx = T(1.0/6.0)*dxinv[0]*dxinv[0];
    T facy = T(1.0/6.0)*dxinv[1]*dxinv[1];

    if (!msk(i,j,k)) {
        x(i,j,k) /= (-T(2.0))*(facx*(sx(i-1,j-1,k)+sx(i,j-1,k)+sx(i-1,j,k)+sx(i,j,k))
                                 +facy*(sy(i-1,j-1,k)+sy(i,j-1,k)+sy(i-1,j,k)+sy(i,j,k)));
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_normalize_aa (int i, int j, int k, Array4<T> const& x, Array4<std::type_identity_t<T> const> const& sig,
                           Array4<int const> const& msk, GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    T facx = T(1.0/6.0)*dxinv[0]*dxinv[0];
    T facy = T(1.0/6.0)*dxinv[1]*dxinv[1];
    T fxy = facx + facy;

    if (!msk(i,j,k)) {
        x(i,j,k) /= (-T(2.0))*fxy*(sig(i-1,j-1,k)+sig(i,j-1,k)+sig(i-1,j,k)+sig(i,j,k));
    }
}

template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void mlndlap_jacobi_ha (int i, int j, int k, Array4<T> const& sol, std::type_identity_t<T> Ax,
                        Array4<std::type_identity_t<T> const> const& rhs, Array4<std::type_identity_t<T> const> const& sx,
                        Array4<std::type_identity_t<T> const> const& sy, Array4<int const> const& msk,
                        GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    T facx = -T(2.0/6.0)*dxinv[0]*dxinv[0];
    T facy = -T(2.0/6.0)*dxinv[1]*dxinv[1];

    if (msk(i,j,k)) {
        sol(i,j,k) = T(0.0);
    } else {
        sol(i,j,k) += T(2.0/3.0) * (rhs(i,j,k) - Ax)
            / (facx*(sx(i-1,j-1,k)+sx(i,j-1,k)+sx(i-1,j,k)+sx(i,j,k))
            +  facy*(sy(i-1,j-1,k)+sy(i,j-1,k)+sy(i-1,j,k)+sy(i,j,k)));
    }
}

template <typename T>
inline
void mlndlap_jacobi_ha (Box const& bx, Array4<T> const& sol, Array4<std::type_identity_t<T> const> const& Ax,
                        Array4<std::type_identity_t<T> const> const& rhs, Array4<std::type_identity_t<T> const> const& sx,
                        Array4<std::type_identity_t<T> const> const& sy, Array4<int const> const& msk,
                        GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    T facx = -T(2.0/6.0)*dxinv[0]*dxinv[0];
    T facy = -T(2.0/6.0)*dxinv[1]*dxinv[1];

    amrex::LoopConcurrentOnCpu(bx, [&] (int i, int j, int k) noexcept
    {
        if (msk(i,j,k)) {
            sol(i,j,k) = T(0.0);
        } else {
            sol(i,j,k) += T(2.0/3.0) * (rhs(i,j,k) - Ax(i,j,k))
                / (facx*(sx(i-1,j-1,k)+sx(i,j-1,k)+sx(i-1,j,k)+sx(i,j,k))
                +  facy*(sy(i-1,j-1,k)+sy(i,j-1,k)+sy(i-1,j,k)+sy(i,j,k)));
        }
    });
}

template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void mlndlap_jacobi_aa (int i, int j, int k, Array4<T> const& sol, std::type_identity_t<T> Ax,
                        Array4<std::type_identity_t<T> const> const& rhs, Array4<std::type_identity_t<T> const> const& sig,
                        Array4<int const> const& msk, GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    T fac = -T(2.0/6.0)*(dxinv[0]*dxinv[0] + dxinv[1]*dxinv[1]);

    if (msk(i,j,k)) {
        sol(i,j,k) = T(0.0);
    } else {
        sol(i,j,k) += T(2.0/3.0) * (rhs(i,j,k) - Ax)
            / (fac*(sig(i-1,j-1,k)+sig(i,j-1,k)+sig(i-1,j,k)+sig(i,j,k)));
    }
}

template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void mlndlap_jacobi_c (int i, int j, int k, Array4<T> const& sol, std::type_identity_t<T> Ax,
                       Array4<std::type_identity_t<T> const> const& rhs, std::type_identity_t<T> sig,
                       Array4<int const> const& msk, GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    T fac = -T(2.0/6.0)*(dxinv[0]*dxinv[0] + dxinv[1]*dxinv[1]);

    if (msk(i,j,k)) {
        sol(i,j,k) = T(0.0);
    } else {
        sol(i,j,k) += T(2.0/3.0) * (rhs(i,j,k) - Ax)
            / (fac*T(4.)*sig);
    }
}

template <typename T>
inline
void mlndlap_jacobi_aa (Box const& bx, Array4<T> const& sol, Array4<std::type_identity_t<T> const> const& Ax,
                        Array4<std::type_identity_t<T> const> const& rhs, Array4<std::type_identity_t<T> const> const& sig,
                        Array4<int const> const& msk, GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    T fac = -T(2.0/6.0)*(dxinv[0]*dxinv[0] + dxinv[1]*dxinv[1]);

    amrex::LoopConcurrentOnCpu(bx, [&] (int i, int j, int k) noexcept
    {
        if (msk(i,j,k)) {
            sol(i,j,k) = T(0.0);
        } else {
            sol(i,j,k) += T(2.0/3.0) * (rhs(i,j,k) - Ax(i,j,k))
                / (fac*(sig(i-1,j-1,k)+sig(i,j-1,k)+sig(i-1,j,k)+sig(i,j,k)));
        }
    });
}

template <typename T>
inline
void mlndlap_jacobi_c (Box const& bx, Array4<T> const& sol, Array4<std::type_identity_t<T> const> const& Ax,
                       Array4<std::type_identity_t<T> const> const& rhs, std::type_identity_t<T> sig,
                       Array4<int const> const& msk, GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    T fac = -T(2.0/6.0)*(dxinv[0]*dxinv[0] + dxinv[1]*dxinv[1]);

    amrex::LoopConcurrentOnCpu(bx, [&] (int i, int j, int k) noexcept
    {
        if (msk(i,j,k)) {
            sol(i,j,k) = T(0.0);
        } else {
            sol(i,j,k) += T(2.0/3.0) * (rhs(i,j,k) - Ax(i,j,k))
                / (fac*T(4.)*sig);
        }
    });
}

template <typename T>
inline
void mlndlap_gauss_seidel_ha (Box const& bx, Array4<T> const& sol,
                              Array4<std::type_identity_t<T> const> const& rhs, Array4<std::type_identity_t<T> const> const& sx,
                              Array4<std::type_identity_t<T> const> const& sy, Array4<int const> const& msk,
                              GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv,
                              bool is_rz) noexcept
{
    T facx = T(1.0/6.0)*dxinv[0]*dxinv[0];
    T facy = T(1.0/6.0)*dxinv[1]*dxinv[1];

    amrex::LoopOnCpu(bx, [&] (int i, int j, int k) noexcept
    {
        if (msk(i,j,k)) {
            sol(i,j,k) = T(0.0);
        } else {
            T s0 = T(-2.0)*(facx*(sx(i-1,j-1,k)+sx(i,j-1,k)+sx(i-1,j,k)+sx(i,j,k))
                                 +facy*(sy(i-1,j-1,k)+sy(i,j-1,k)+sy(i-1,j,k)+sy(i,j,k)));

            T Ax = sol(i-1,j-1,k)*(facx*sx(i-1,j-1,k)+facy*sy(i-1,j-1,k))
                    + sol(i+1,j-1,k)*(facx*sx(i  ,j-1,k)+facy*sy(i  ,j-1,k))
                    + sol(i-1,j+1,k)*(facx*sx(i-1,j  ,k)+facy*sy(i-1,j  ,k))
                    + sol(i+1,j+1,k)*(facx*sx(i  ,j  ,k)+facy*sy(i  ,j  ,k))
                    + sol(i-1,j,k)*(T(2.0)*facx*(sx(i-1,j-1,k)+sx(i-1,j,k))
                                        -     facy*(sy(i-1,j-1,k)+sy(i-1,j,k)))
                    + sol(i+1,j,k)*(T(2.0)*facx*(sx(i  ,j-1,k)+sx(i  ,j,k))
                                        -     facy*(sy(i  ,j-1,k)+sy(i  ,j,k)))
                    + sol(i,j-1,k)*(   -facx*(sx(i-1,j-1,k)+sx(i,j-1,k))
                             +T(2.0)*facy*(sy(i-1,j-1,k)+sy(i,j-1,k)))
                    + sol(i,j+1,k)*(   -facx*(sx(i-1,j  ,k)+sx(i,j  ,k))
                             +T(2.0)*facy*(sy(i-1,j  ,k)+sy(i,j  ,k)))
                    + sol(i,j,k)*s0;

            if (is_rz) {
                T fp = facy / static_cast<T>(2*i+1);
                T fm = facy / static_cast<T>(2*i-1);
                T frzlo = fm*sy(i-1,j-1,k)-fp*sy(i,j-1,k);
                T frzhi = fm*sy(i-1,j  ,k)-fp*sy(i,j  ,k);
                s0 += - frzhi - frzlo;
                Ax += frzhi*(sol(i,j+1,k)-sol(i,j,k))
                    + frzlo*(sol(i,j-1,k)-sol(i,j,k));
            }

            sol(i,j,k) += (rhs(i,j,k) - Ax) / s0;
        }
    });
}

template <typename T>
inline
void mlndlap_gauss_seidel_aa (Box const& bx, Array4<T> const& sol,
                              Array4<std::type_identity_t<T> const> const& rhs, Array4<std::type_identity_t<T> const> const& sig,
                              Array4<int const> const& msk,
                              GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv,
                              bool is_rz) noexcept
{
    T facx = T(1.0/6.0)*dxinv[0]*dxinv[0];
    T facy = T(1.0/6.0)*dxinv[1]*dxinv[1];
    T fxy = facx + facy;
    T f2xmy = T(2.0)*facx - facy;
    T fmx2y = T(2.0)*facy - facx;

    amrex::LoopOnCpu(bx, [&] (int i, int j, int k) noexcept
    {
        if (msk(i,j,k)) {
            sol(i,j,k) = T(0.0);
        } else {
            T s0 = (-T(2.0))*fxy*(sig(i-1,j-1,k)+sig(i,j-1,k)+sig(i-1,j,k)+sig(i,j,k));
            T Ax =   sol(i-1,j-1,k)*fxy*sig(i-1,j-1,k)
                      + sol(i+1,j-1,k)*fxy*sig(i  ,j-1,k)
                      + sol(i-1,j+1,k)*fxy*sig(i-1,j  ,k)
                      + sol(i+1,j+1,k)*fxy*sig(i  ,j  ,k)
                      + sol(i-1,j,k)*f2xmy*(sig(i-1,j-1,k)+sig(i-1,j,k))
                      + sol(i+1,j,k)*f2xmy*(sig(i  ,j-1,k)+sig(i  ,j,k))
                      + sol(i,j-1,k)*fmx2y*(sig(i-1,j-1,k)+sig(i,j-1,k))
                      + sol(i,j+1,k)*fmx2y*(sig(i-1,j  ,k)+sig(i,j  ,k))
                      + sol(i,j,k)*s0;

            if (is_rz) {
                T fp = facy / static_cast<T>(2*i+1);
                T fm = facy / static_cast<T>(2*i-1);
                T frzlo = fm*sig(i-1,j-1,k)-fp*sig(i,j-1,k);
                T frzhi = fm*sig(i-1,j  ,k)-fp*sig(i,j  ,k);
                s0 += - frzhi - frzlo;
                Ax += frzhi*(sol(i,j+1,k)-sol(i,j,k))
                    + frzlo*(sol(i,j-1,k)-sol(i,j,k));
            }

            sol(i,j,k) += (rhs(i,j,k) - Ax) / s0;
        }
    });
}

template <typename T>
inline
void mlndlap_gauss_seidel_c (Box const& bx, Array4<T> const& sol,
                             Array4<std::type_identity_t<T> const> const& rhs, std::type_identity_t<T> sig,
                             Array4<int const> const& msk,
                             GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv,
                             bool is_rz) noexcept
{
    T facx = T(1.0/6.0)*dxinv[0]*dxinv[0];
    T facy = T(1.0/6.0)*dxinv[1]*dxinv[1];
    T fxy = facx + facy;
    T f2xmy = T(2.0)*facx - facy;
    T fmx2y = T(2.0)*facy - facx;

    amrex::LoopOnCpu(bx, [&] (int i, int j, int k) noexcept
    {
        if (msk(i,j,k)) {
            sol(i,j,k) = T(0.0);
        } else {
            T s0 = (-T(2.0))*fxy*T(4.);
            T Ax =   sol(i-1,j-1,k)*fxy
                      + sol(i+1,j-1,k)*fxy
                      + sol(i-1,j+1,k)*fxy
                      + sol(i+1,j+1,k)*fxy
                      + sol(i-1,j,k)*f2xmy*T(2.)
                      + sol(i+1,j,k)*f2xmy*T(2.)
                      + sol(i,j-1,k)*fmx2y*T(2.)
                      + sol(i,j+1,k)*fmx2y*T(2.)
                      + sol(i,j,k)*s0;

            if (is_rz) {
                T fp = facy / static_cast<T>(2*i+1);
                T fm = facy / static_cast<T>(2*i-1);
                T frzlo = fm-fp;
                T frzhi = fm-fp;
                s0 += - frzhi - frzlo;
                Ax += frzhi*(sol(i,j+1,k)-sol(i,j,k))
                    + frzlo*(sol(i,j-1,k)-sol(i,j,k));
            }

            sol(i,j,k) += (rhs(i,j,k) - Ax*sig) / (s0*sig);
        }
    });
}

template <typename T>
AMREX_FORCE_INLINE
void tridiagonal_solve (Array1D<T,0,31>& a_ls, Array1D<T,0,31>& b_ls, Array1D<T,0,31>& c_ls,
                        Array1D<T,0,31>& r_ls, Array1D<T,0,31>& u_ls, Array1D<T,0,31>& gam,
                        int ilen )
{
    T bet = b_ls(0);
    u_ls(0) = r_ls(0) / bet;

    for (int i = 1; i <= ilen - 1; i++) {
        gam(i) = c_ls(i-1) / bet;
        bet = b_ls(i) - a_ls(i)*gam(i);
        if (bet == 0) { amrex::Abort(">>>TRIDIAG FAILED"); }
        u_ls(i) = (r_ls(i)-a_ls(i)*u_ls(i-1)) / bet;
    }
    for (int i = ilen-2; i >= 0; i--) {
        u_ls(i) = u_ls(i) - gam(i+1)*u_ls(i+1);
    }
}

template <typename T>
inline
void mlndlap_gauss_seidel_with_line_solve_aa (Box const& bx, Array4<T> const& sol,
                              Array4<std::type_identity_t<T> const> const& rhs, Array4<std::type_identity_t<T> const> const& sig,
                              Array4<int const> const& msk,
                              GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv,
                              bool is_rz)
{
    T facx = T(1.0/6.0)*dxinv[0]*dxinv[0];
    T facy = T(1.0/6.0)*dxinv[1]*dxinv[1];
    T fxy = facx + facy;
    T f2xmy = T(2.0)*facx - facy;
    T fmx2y = T(2.0)*facy - facx;

    if (is_rz) {
        amrex::Abort("mlndlap_gauss_seidel_with_line_solve_aa is not implemented in r-z 2D ");
    }

    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    int idir = -1;
    int ilen = 33;
    int k = 0;
    if (dxinv[0] <= dxinv[1])  {
        idir = 1;
        ilen = hi.y - lo.y + 1;
    }
    if (dxinv[1] <= dxinv[0]) {
        idir = 0;
        ilen = hi.x - lo.x + 1;
    }

    if (ilen > 32) {
        amrex::Abort("mlndlap_gauss_seidel_with_line_solve_aa is hard-wired to be no longer than 32");
    }

    Array1D<T,0,31> a_ls,b_ls,c_ls,u_ls,r_ls,gam;

    if (idir == 1) {
        for (int i = lo.x; i <= hi.x; ++i)
        {
            for (int j = lo.y; j <= hi.y; ++j)
            {
                if (msk(i,j,k)) {
                    a_ls(j-lo.y) = 0.;
                    b_ls(j-lo.y) = 1.;
                    c_ls(j-lo.y) = 0.;
                    u_ls(j-lo.y) = 0.;
                    r_ls(j-lo.y) = 0.;
                }
                else
                {
                    T s0 = (-T(2.0))*fxy*(sig(i-1,j-1,k)+sig(i,j-1,k)+sig(i-1,j,k)+sig(i,j,k));

                    T Ax = sol(i-1,j-1,k)*fxy*sig(i-1,j-1,k)
                            + sol(i+1,j-1,k)*fxy*sig(i  ,j-1,k)
                            + sol(i-1,j+1,k)*fxy*sig(i-1,j  ,k)
                            + sol(i+1,j+1,k)*fxy*sig(i  ,j  ,k)
                            + sol(i-1,j,k)*f2xmy*(sig(i-1,j-1,k)+sig(i-1,j,k))
                            + sol(i+1,j,k)*f2xmy*(sig(i  ,j-1,k)+sig(i  ,j,k));

                    a_ls(j-lo.y) = fmx2y*(sig(i-1,j-1,k)+sig(i,j-1,k));
                    b_ls(j-lo.y) = s0;
                    c_ls(j-lo.y) = fmx2y*(sig(i-1,j  ,k)+sig(i,j  ,k));
                    u_ls(j-lo.y) = 0.;
                    r_ls(j-lo.y) = rhs(i,j,k) - Ax;
                }
            }
            tridiagonal_solve(a_ls, b_ls, c_ls, r_ls, u_ls, gam, ilen);

            for (int j = lo.y; j <= hi.y; ++j)
            {
                sol(i,j,k) = u_ls(j-lo.y);
            }
        }
    } else if (idir == 0) {
        for (int j = lo.y ;j <= hi.y; ++j)
        {
            for (int i = lo.x; i <= hi.x; ++i)
            {
                if (msk(i,j,k)) {
                    a_ls(i-lo.x) = 0.;
                    b_ls(i-lo.x) = 1.;
                    c_ls(i-lo.x) = 0.;
                    u_ls(i-lo.x) = 0.;
                    r_ls(i-lo.x) = 0.;
                }
                else
                {
                    T s0 = (-T(2.0))*fxy*(sig(i-1,j-1,k)+sig(i,j-1,k)+sig(i-1,j,k)+sig(i,j,k));

                    T Ax = sol(i-1,j-1,k)*fxy*sig(i-1,j-1,k)
                            + sol(i+1,j-1,k)*fxy*sig(i  ,j-1,k)
                            + sol(i-1,j+1,k)*fxy*sig(i-1,j  ,k)
                            + sol(i+1,j+1,k)*fxy*sig(i  ,j  ,k)
                            + sol(i,j-1,k)*fmx2y*(sig(i-1,j-1,k)+sig(i,j-1,k))
                            + sol(i,j+1,k)*fmx2y*(sig(i-1,j  ,k)+sig(i,j  ,k));

                    a_ls(i-lo.x) = f2xmy*(sig(i-1,j-1,k)+sig(i-1,j,k));
                    b_ls(i-lo.x) = s0;
                    c_ls(i-lo.x) = f2xmy*(sig(i  ,j-1,k)+sig(i  ,j,k));
                    u_ls(i-lo.x) = 0.;
                    r_ls(i-lo.x) = rhs(i,j,k) - Ax;

                }
            }
            tridiagonal_solve(a_ls, b_ls, c_ls, r_ls, u_ls, gam, ilen);

            for (int i = lo.x; i <= hi.x; ++i)
            {
                sol(i,j,k) = u_ls(i-lo.x);
            }
        }
    } else {
        amrex::Abort("mlndlap_gauss_seidel_with_line_solve_aa is wrong direction.");
    }

}

//
// interpolation
//

    template <typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T aa_interp_line_x (Array4<T const> const& crse, Array4<std::type_identity_t<T> const> const& sig,
                           int i, int j, int ic, int jc) noexcept
    {
        T w1 = sig(i-1,j-1,0) + sig(i-1,j,0);
        T w2 = sig(i  ,j-1,0) + sig(i  ,j,0);
        return (w1*crse(ic,jc,0)+w2*crse(ic+1,jc,0))/(w1+w2);
    }

    template <typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T aa_interp_line_y (Array4<T const> const& crse, Array4<std::type_identity_t<T> const> const& sig,
                           int i, int j, int ic, int jc) noexcept
    {
        T w1 = sig(i-1,j-1,0) + sig(i,j-1,0);
        T w2 = sig(i-1,j  ,0) + sig(i,j  ,0);
        return (w1*crse(ic,jc,0)+w2*crse(ic,jc+1,0))/(w1+w2);
    }

    template <typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T aa_interp_face_xy (Array4<T const> const& crse, Array4<std::type_identity_t<T> const> const& sig,
                            int i, int j, int ic, int jc) noexcept
    {
        T w1 = sig(i-1,j-1,0) + sig(i-1,j,0);
        T w2 = sig(i  ,j-1,0) + sig(i  ,j,0);
        T w3 = sig(i-1,j-1,0) + sig(i,j-1,0);
        T w4 = sig(i-1,j  ,0) + sig(i,j  ,0);
        return (w1 * aa_interp_line_y(crse,sig,i-1,j  ,ic  ,jc  ) +
                w2 * aa_interp_line_y(crse,sig,i+1,j  ,ic+1,jc  ) +
                w3 * aa_interp_line_x(crse,sig,i  ,j-1,ic  ,jc  ) +
                w4 * aa_interp_line_x(crse,sig,i  ,j+1,ic  ,jc+1)) / (w1+w2+w3+w4);
    }

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_interpadd_c (int i, int j, int, Array4<T> const& fine,
                          Array4<std::type_identity_t<T> const> const& crse,
                           Array4<int const> const& msk) noexcept
{
    if (!msk(i,j,0)) {
        int ic = amrex::coarsen(i,2);
        int jc = amrex::coarsen(j,2);
        bool i_is_odd = (ic*2 != i);
        bool j_is_odd = (jc*2 != j);
        if (i_is_odd && j_is_odd) {
            // Node on a X-Y face
            fine(i,j,0) += T(0.25) * (crse(ic  ,jc  ,0) +
                                         crse(ic+1,jc  ,0) +
                                         crse(ic  ,jc+1,0) +
                                         crse(ic+1,jc+1,0));
        } else if (i_is_odd) {
            // Node on X line
            fine(i,j,0) += T(0.5) * (crse(ic,jc,0)+crse(ic+1,jc,0));
        } else if (j_is_odd) {
            // Node on Y line
            fine(i,j,0) += T(0.5) * (crse(ic,jc,0)+crse(ic,jc+1,0));
        } else {
            // Node coincident with coarse node
            fine(i,j,0) += crse(ic,jc,0);
        }
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_interpadd_aa (int i, int j, int, Array4<T> const& fine,
                           Array4<std::type_identity_t<T> const> const& crse, Array4<std::type_identity_t<T> const> const& sig,
                           Array4<int const> const& msk) noexcept
{
    if (!msk(i,j,0)) {
        int ic = amrex::coarsen(i,2);
        int jc = amrex::coarsen(j,2);
        bool i_is_odd = (ic*2 != i);
        bool j_is_odd = (jc*2 != j);
        if (i_is_odd && j_is_odd) {
            // Node on a X-Y face
            fine(i,j,0) += aa_interp_face_xy(crse,sig,i,j,ic,jc);
        } else if (i_is_odd) {
            // Node on X line
            fine(i,j,0) += aa_interp_line_x(crse,sig,i,j,ic,jc);
        } else if (j_is_odd) {
            // Node on Y line
            fine(i,j,0) += aa_interp_line_y(crse,sig,i,j,ic,jc);
        } else {
            // Node coincident with coarse node
            fine(i,j,0) += crse(ic,jc,0);
        }
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_semi_interpadd_aa (int i, int j, int, Array4<T> const& fine,
                           Array4<std::type_identity_t<T> const> const& crse, Array4<std::type_identity_t<T> const> const& sig,
                           Array4<int const> const& msk, int idir) noexcept
{
    if (idir == 1) {
        if (!msk(i,j,0)) {
            int ic = amrex::coarsen(i,2);
            int jc = j;
            bool i_is_odd = (ic*2 != i);
            if (i_is_odd) {
                // Node on X line
                fine(i,j,0) += aa_interp_line_x(crse,sig,i,j,ic,jc);
            } else {
                //Node coincident with coarse node
                fine(i,j,0) += crse(ic,jc,0);
            }
        }
    } else if (idir == 0 ) {
        if (!msk(i,j,0)) {
            int ic = i;
            int jc = amrex::coarsen(j,2);
            bool j_is_odd = (jc*2 != j);
            if (j_is_odd) {
                // Node on Y line
                fine(i,j,0) += aa_interp_line_y(crse,sig,i,j,ic,jc);
            } else {
                //Node coincident with coarse node
                fine(i,j,0) += crse(ic,jc,0);
            }
        }
    }
}

    template <typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T ha_interp_face_xy (Array4<T const> const& crse,
                            Array4<std::type_identity_t<T> const> const& sigx, Array4<std::type_identity_t<T> const> const& sigy,
                            int i, int j, int ic, int jc) noexcept
    {
        T w1 = sigx(i-1,j-1,0) + sigx(i-1,j,0);
        T w2 = sigx(i  ,j-1,0) + sigx(i  ,j,0);
        T w3 = sigy(i-1,j-1,0) + sigy(i,j-1,0);
        T w4 = sigy(i-1,j  ,0) + sigy(i,j  ,0);
        return (w1 * aa_interp_line_y(crse,sigy,i-1,j  ,ic  ,jc  ) +
                w2 * aa_interp_line_y(crse,sigy,i+1,j  ,ic+1,jc  ) +
                w3 * aa_interp_line_x(crse,sigx,i  ,j-1,ic  ,jc  ) +
                w4 * aa_interp_line_x(crse,sigx,i  ,j+1,ic  ,jc+1)) / (w1+w2+w3+w4);
    }

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_interpadd_ha (int i, int j, int,
                           Array4<T> const& fine, Array4<std::type_identity_t<T> const> const& crse,
                           Array4<std::type_identity_t<T> const> const& sigx, Array4<std::type_identity_t<T> const> const& sigy,
                           Array4<int const> const& msk) noexcept
{
    if (!msk(i,j,0)) {
        int ic = amrex::coarsen(i,2);
        int jc = amrex::coarsen(j,2);
        bool i_is_odd = (ic*2 != i);
        bool j_is_odd = (jc*2 != j);
        if (i_is_odd && j_is_odd) {
            // Node on a X-Y face
            fine(i,j,0) += ha_interp_face_xy(crse,sigx,sigy,i,j,ic,jc);
        } else if (i_is_odd) {
            // Node on X line
            fine(i,j,0) += aa_interp_line_x(crse,sigx,i,j,ic,jc);
        } else if (j_is_odd) {
            // Node on Y line
            fine(i,j,0) += aa_interp_line_y(crse,sigy,i,j,ic,jc);
        } else {
            // Node coincident with coarse node
            fine(i,j,0) += crse(ic,jc,0);
        }
    }
}

//
// rhs & u
//

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_divu (int i, int j, int k, Array4<T> const& rhs, Array4<std::type_identity_t<T> const> const& vel,
                   Array4<int const> const& msk, GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv,
                   Box const& nodal_domain,
                   GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bclo,
                   GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bchi,
                   bool is_rz) noexcept
{
    T facx = T(0.5)*dxinv[0];
    T facy = T(0.5)*dxinv[1];

    const auto domlo = amrex::lbound(nodal_domain);
    const auto domhi = amrex::ubound(nodal_domain);

    if (msk(i,j,k)) {
        rhs(i,j,k) = T(0.0);
    } else {

        T zero_ilo = T(1.0);
        T zero_ihi = T(1.0);
        T zero_jlo = T(1.0);
        T zero_jhi = T(1.0);

        // The nodal divergence operator should not see the tangential velocity
        //     at an inflow face
        if ((bclo[0] == LinOpBCType::Neumann || bclo[0] == LinOpBCType::inflow)
            && i == domlo.x)
        {
            zero_ilo = T(0.0);
        }
        if ((bchi[0] == LinOpBCType::Neumann || bchi[0] == LinOpBCType::inflow)
            && i == domhi.x)
        {
            zero_ihi = T(0.0);
        }
        if ((bclo[1] == LinOpBCType::Neumann || bclo[1] == LinOpBCType::inflow)
            && j == domlo.y)
        {
            zero_jlo = T(0.0);
        }
        if ((bchi[1] == LinOpBCType::Neumann || bchi[1] == LinOpBCType::inflow)
            && j == domhi.y)
        {
            zero_jhi = T(0.0);
        }

        rhs(i,j,k) = facx*(-vel(i-1,j-1,k,0)*zero_jlo + vel(i,j-1,k,0)*zero_jlo
                           -vel(i-1,j  ,k,0)*zero_jhi + vel(i,j  ,k,0)*zero_jhi)
                   + facy*(-vel(i-1,j-1,k,1)*zero_ilo - vel(i,j-1,k,1)*zero_ihi
                           +vel(i-1,j  ,k,1)*zero_ilo + vel(i,j  ,k,1)*zero_ihi);
        if (is_rz) {
            T fm = facy / static_cast<T>(6*i-3);
            T fp = facy / static_cast<T>(6*i+3);
            rhs(i,j,k) += fm*(vel(i-1,j,k,1)-vel(i-1,j-1,k,1))*zero_ilo
                        - fp*(vel(i  ,j,k,1)-vel(i  ,j-1,k,1))*zero_ihi;
        }
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T mlndlap_rhcc (int i, int j, int k, Array4<T const> const& rhcc,
                   Array4<int const> const& msk) noexcept
{
    T r;
    if (msk(i,j,k)) {
        r = T(0.0);
    } else {
        r = T(0.25) * (rhcc(i-1,j-1,k)+rhcc(i,j-1,k)+rhcc(i-1,j,k)+rhcc(i,j,k));
    }
    return r;
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_mknewu (int i, int j, int k, Array4<T> const& u, Array4<std::type_identity_t<T> const> const& p,
                     Array4<std::type_identity_t<T> const> const& sig,  GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv,
                     bool is_rz) noexcept
{
    T facx = T(0.5)*dxinv[0];
    T facy = T(0.5)*dxinv[1];
    u(i,j,k,0) -= sig(i,j,k)*facx*(-p(i,j,k)+p(i+1,j,k)-p(i,j+1,k)+p(i+1,j+1,k));
    u(i,j,k,1) -= sig(i,j,k)*facy*(-p(i,j,k)-p(i+1,j,k)+p(i,j+1,k)+p(i+1,j+1,k));
    if (is_rz) {
        T frz = sig(i,j,k)*facy / static_cast<T>(6*i+3);
        u(i,j,k,1) -= frz*(p(i,j,k)-p(i+1,j,k)-p(i,j+1,k)+p(i+1,j+1,k));
    }
}

// Harmonic-average (anisotropic sigma) variant.  Used when the node-
// centered Poisson operator was set up with a multi-component sigma
// (see MLNodeLaplacian::setSigma and setMapped) so each velocity
// component is corrected with its own sigma direction.  Under RZ the
// axial correction is applied to the y-component (which is the axial
// direction in AMReX's 2D RZ convention).
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_mknewu_ha (int i, int j, int k, Array4<T> const& u, Array4<std::type_identity_t<T> const> const& p,
                        Array4<std::type_identity_t<T> const> const& sx,
                        Array4<std::type_identity_t<T> const> const& sy,
                        GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv,
                        bool is_rz) noexcept
{
    T facx = T(0.5)*dxinv[0];
    T facy = T(0.5)*dxinv[1];
    u(i,j,k,0) -= sx(i,j,k)*facx*(-p(i,j,k)+p(i+1,j,k)-p(i,j+1,k)+p(i+1,j+1,k));
    u(i,j,k,1) -= sy(i,j,k)*facy*(-p(i,j,k)-p(i+1,j,k)+p(i,j+1,k)+p(i+1,j+1,k));
    if (is_rz) {
        T frz = sy(i,j,k)*facy / static_cast<T>(6*i+3);
        u(i,j,k,1) -= frz*(p(i,j,k)-p(i+1,j,k)-p(i,j+1,k)+p(i+1,j+1,k));
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_mknewu_c (int i, int j, int k, Array4<T> const& u, Array4<std::type_identity_t<T> const> const& p,
                       std::type_identity_t<T> sig,  GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv,
                       bool is_rz) noexcept
{
    T facx = T(0.5)*dxinv[0];
    T facy = T(0.5)*dxinv[1];
    u(i,j,k,0) -= sig*facx*(-p(i,j,k)+p(i+1,j,k)-p(i,j+1,k)+p(i+1,j+1,k));
    u(i,j,k,1) -= sig*facy*(-p(i,j,k)-p(i+1,j,k)+p(i,j+1,k)+p(i+1,j+1,k));
    if (is_rz) {
        T frz = sig*facy / static_cast<T>(6*i+3);
        u(i,j,k,1) -= frz*(p(i,j,k)-p(i+1,j,k)-p(i,j+1,k)+p(i+1,j+1,k));
    }
}

    template <typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T mlndlap_sum_Df (int ii, int jj, std::type_identity_t<T> facx, std::type_identity_t<T> facy,
                         Array4<T const> const& vel, Box const& velbx, bool is_rz) noexcept
    {
        T fm = is_rz ? facy / static_cast<T>(6*ii-3) : T(0.0);
        T fp = is_rz ? facy / static_cast<T>(6*ii+3) : T(0.0);

        T Df = T(0.0);
        if (velbx.contains(ii-1,jj-1,0)) {
            Df += -facx*vel(ii-1,jj-1,0,0) - (facy+fm)*vel(ii-1,jj-1,0,1);
        }
        if (velbx.contains(ii,jj-1,0)) {
            Df += facx*vel(ii,jj-1,0,0) - (facy-fp)*vel(ii,jj-1,0,1);
        }
        if (velbx.contains(ii-1,jj,0)) {
            Df += -facx*vel(ii-1,jj,0,0) + (facy+fm)*vel(ii-1,jj,0,1);
        }
        if (velbx.contains(ii,jj,0)) {
            Df += facx*vel(ii,jj,0,0) + (facy-fp)*vel(ii,jj,0,1);
        }
        return Df;
    }

template <int rr, typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_divu_fine_contrib (int i, int j, int /*k*/, Box const& fvbx, Box const& velbx,
                                Array4<T> const& rhs, Array4<std::type_identity_t<T> const> const& vel,
                                Array4<std::type_identity_t<T> const> const& frhs, Array4<int const> const& msk,
                                GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv, bool is_rz) noexcept
{
    const int ii = rr*i;
    const int jj = rr*j;
    if (msk(ii,jj,0)) {
        const T facx = T(0.5)*dxinv[0];
        const T facy = T(0.5)*dxinv[1];

        T Df = T(0.0);

        const int ilo = amrex::max(ii-rr+1, fvbx.smallEnd(0));
        const int ihi = amrex::min(ii+rr-1, fvbx.bigEnd  (0));
        const int jlo = amrex::max(jj-rr+1, fvbx.smallEnd(1));
        const int jhi = amrex::min(jj+rr-1, fvbx.bigEnd  (1));

        for (int joff = jlo; joff <= jhi; ++joff) {
        for (int ioff = ilo; ioff <= ihi; ++ioff) {
            T scale = static_cast<T>((rr-std::abs(ii-ioff)) *
                                           (rr-std::abs(jj-joff)));
            if (fvbx.strictly_contains(ioff,joff,0)) {
                Df += scale * frhs(ioff,joff,0);
            } else {
                Df += scale * mlndlap_sum_Df(ioff, joff, facx, facy, vel, velbx, is_rz);
            }
        }}

        rhs(i,j,0) = Df * (T(1.0)/static_cast<T>(rr*rr*rr*rr));
    } else {
        rhs(i,j,0) = T(0.0);
    }
}

template <int rr, typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_rhcc_fine_contrib (int i, int j, int, Box const& ccbx,
                                Array4<T> const& rhs, Array4<std::type_identity_t<T> const> const& cc,
                                Array4<int const> const& msk) noexcept
{
    const int ii = rr*i;
    const int jj = rr*j;
    if (msk(ii,jj,0)) {
        T tmp = T(0.0);

        const int ilo = amrex::max(ii-rr  , ccbx.smallEnd(0));
        const int ihi = amrex::min(ii+rr-1, ccbx.bigEnd  (0));
        const int jlo = amrex::max(jj-rr  , ccbx.smallEnd(1));
        const int jhi = amrex::min(jj+rr-1, ccbx.bigEnd  (1));

        for (int joff = jlo; joff <= jhi; ++joff) {
        for (int ioff = ilo; ioff <= ihi; ++ioff) {
            T scale = (static_cast<T>(rr)-std::abs(static_cast<T>(ioff-ii)+T(0.5)))
                *        (static_cast<T>(rr)-std::abs(static_cast<T>(joff-jj)+T(0.5)));
            tmp += cc(ioff,joff,0) * scale;
        }}

        rhs(i,j,0) += tmp * (T(1.0)/T(rr*rr*rr*rr));
    }
}

    template <typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T neumann_scale (int i, int j, Box const& nddom,
                        GpuArray<LinOpBCType,AMREX_SPACEDIM> const& bclo,
                        GpuArray<LinOpBCType,AMREX_SPACEDIM> const& bchi) noexcept
    {
        T val = T(1.0);

        const auto ndlo = amrex::lbound(nddom);
        const auto ndhi = amrex::ubound(nddom);

        if ((i == ndlo.x && ( bclo[0] == LinOpBCType::Neumann ||
                              bclo[0] == LinOpBCType::inflow)) ||
            (i == ndhi.x && ( bchi[0] == LinOpBCType::Neumann ||
                              bchi[0] == LinOpBCType::inflow))) {
            val *= T(2.0);
        }

        if ((j == ndlo.y && ( bclo[1] == LinOpBCType::Neumann ||
                              bclo[1] == LinOpBCType::inflow)) ||
            (j == ndhi.y && ( bchi[1] == LinOpBCType::Neumann ||
                              bchi[1] == LinOpBCType::inflow))) {
            val *= T(2.0);
        }

        return val;
    }

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_divu_cf_contrib (int i, int j, int, Array4<T> const& rhs,
                              Array4<std::type_identity_t<T> const> const& vel, Array4<std::type_identity_t<T> const> const& fc,
                              Array4<std::type_identity_t<T> const> const& rhcc, Array4<int const> const& dmsk,
                              Array4<int const> const& ndmsk, Array4<int const> const& ccmsk,
                              bool is_rz, GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv,
                              Box const& ccdom_p, Box const& veldom, Box const& nddom,
                              GpuArray<LinOpBCType,AMREX_SPACEDIM> const& bclo,
                              GpuArray<LinOpBCType,AMREX_SPACEDIM> const& bchi) noexcept
{
    using namespace nodelap_detail;

    if (!dmsk(i,j,0) && ndmsk(i,j,0) == crse_fine_node) {
        T facx = T(0.5) * dxinv[0];
        T facy = T(0.5) * dxinv[1];
        T tmp = fc(i,j,0);

        T fm = is_rz ? facy / static_cast<T>(6*i-3) : T(0.0);
        T fp = is_rz ? facy / static_cast<T>(6*i+3) : T(0.0);

        // Where there is inflow, veldom there is bigger than ccdom_p by one cell.
        // ccdom_p is cc domain grown at periodic boundaries.

        if (ccmsk(i-1,j-1,0) == crse_cell && veldom.contains(i-1,j-1,0)) {
            tmp += -facx*vel(i-1,j-1,0,0) - (facy+fm)*vel(i-1,j-1,0,1);
            if (rhcc && ccdom_p.contains(i-1,j-1,0)) {
                tmp += T(0.25) * rhcc(i-1,j-1,0);
            }
        }

        if (ccmsk(i,j-1,0) == crse_cell && veldom.contains(i,j-1,0)) {
            tmp += facx*vel(i,j-1,0,0) - (facy-fp)*vel(i,j-1,0,1);
            if (rhcc && ccdom_p.contains(i,j-1,0)) {
                tmp += T(0.25) * rhcc(i,j-1,0);
            }
        }

        if (ccmsk(i-1,j,0) == crse_cell && veldom.contains(i-1,j,0)) {
            tmp += -facx*vel(i-1,j,0,0) + (facy+fm)*vel(i-1,j,0,1);
            if (rhcc && ccdom_p.contains(i-1,j,0)) {
                tmp += T(0.25) * rhcc(i-1,j,0);
            }
        }

        if (ccmsk(i,j,0) == crse_cell && veldom.contains(i,j,0)) {
            tmp += facx*vel(i,j,0,0) + (facy-fp)*vel(i,j,0,1);
            if (rhcc && ccdom_p.contains(i,j,0)) {
                tmp += T(0.25) * rhcc(i,j,0);
            }
        }

        rhs(i,j,0) = tmp * neumann_scale<T>(i, j, nddom, bclo, bchi);
    }
}

//
// residual
//

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_crse_resid (int i, int j, int k, Array4<T> const& resid,
                         Array4<std::type_identity_t<T> const> const& rhs, Array4<int const> const& msk,
                         Box const& nddom, GpuArray<LinOpBCType,AMREX_SPACEDIM> const& bclo,
                         GpuArray<LinOpBCType,AMREX_SPACEDIM> const& bchi,
                         bool neumann_doubling) noexcept
{
    if ( msk(i-1,j-1,k  ) == 0 ||
         msk(i  ,j-1,k  ) == 0 ||
         msk(i-1,j  ,k  ) == 0 ||
         msk(i  ,j  ,k  ) == 0 )
    {
        T fac = T(1.0);
        if (neumann_doubling) {
            const auto ndlo = amrex::lbound(nddom);
            const auto ndhi = amrex::ubound(nddom);
            if ((i == ndlo.x && ( bclo[0] == LinOpBCType::Neumann ||
                                  bclo[0] == LinOpBCType::inflow)) ||
                (i == ndhi.x && ( bchi[0] == LinOpBCType::Neumann ||
                                  bchi[0] == LinOpBCType::inflow))) {
                fac *= T(2.);
            }
            if ((j == ndlo.y && ( bclo[1] == LinOpBCType::Neumann ||
                                  bclo[1] == LinOpBCType::inflow)) ||
                (j == ndhi.y && ( bchi[1] == LinOpBCType::Neumann ||
                                  bchi[1] == LinOpBCType::inflow))) {
                fac *= T(2.);
            }
        }

        resid(i,j,k) = (rhs(i,j,k) - resid(i,j,k)) * fac;
    } else {
        resid(i,j,k) = T(0.);
    }
}

//
// sync residual
//

    template <typename P, typename S, typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T mlndlap_sum_Ax (P const& pred, S const& sig, int i, int j, std::type_identity_t<T> facx, std::type_identity_t<T> facy,
                         Array4<T const> const& phi, bool is_rz) noexcept
    {
        T Ax = T(0.0);
        T fp = T(0.0), fm = T(0.0);
        if (is_rz) {
            fp = facy / static_cast<T>(2*i+1);
            fm = facy / static_cast<T>(2*i-1);
        }
        if (pred(i-1,j-1)) {
            Ax += sig(i-1,j-1,0)*(facx*(T(2.)*(phi(i-1,j  ,0)-phi(i  ,j  ,0))
                                           +     (phi(i-1,j-1,0)-phi(i  ,j-1,0)))
                                + facy*(T(2.)*(phi(i  ,j-1,0)-phi(i  ,j  ,0))
                                           +     (phi(i-1,j-1,0)-phi(i-1,j  ,0)))
                                + fm  *          (phi(i  ,j-1,0)-phi(i  ,j  ,0)));
        }
        if (pred(i,j-1)) {
            Ax += sig(i,j-1,0)*(facx*(T(2.)*(phi(i+1,j  ,0)-phi(i  ,j  ,0))
                                         +     (phi(i+1,j-1,0)-phi(i  ,j-1,0)))
                              + facy*(T(2.)*(phi(i  ,j-1,0)-phi(i  ,j  ,0))
                                         +     (phi(i+1,j-1,0)-phi(i+1,j  ,0)))
                              - fp  *          (phi(i  ,j-1,0)-phi(i  ,j  ,0)));
        }
        if (pred(i-1,j)) {
            Ax += sig(i-1,j,0)*(facx*(T(2.)*(phi(i-1,j  ,0)-phi(i  ,j  ,0))
                                         +     (phi(i-1,j+1,0)-phi(i  ,j+1,0)))
                              + facy*(T(2.)*(phi(i  ,j+1,0)-phi(i  ,j  ,0))
                                         +     (phi(i-1,j+1,0)-phi(i-1,j  ,0)))
                              + fm  *          (phi(i  ,j+1,0)-phi(i  ,j  ,0)));
        }
        if (pred(i,j)) {
            Ax += sig(i,j,0)*(facx*(T(2.)*(phi(i+1,j  ,0)-phi(i  ,j  ,0))
                                      +      (phi(i+1,j+1,0)-phi(i  ,j+1,0)))
                            + facy*(T(2.)*(phi(i  ,j+1,0)-phi(i  ,j  ,0))
                                      +      (phi(i+1,j+1,0)-phi(i+1,j  ,0)))
                            - fp  *          (phi(i  ,j+1,0)-phi(i  ,j  ,0)));
        }
        return Ax;
    }

    template <int rr, typename S, typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void mlndlap_Ax_fine_contrib_doit (S const& sig, int i, int j, Box const& ndbx, Box const& ccbx,
                                       Array4<T> const& f, Array4<std::type_identity_t<T> const> const& res,
                                       Array4<std::type_identity_t<T> const> const& rhs,
                                       Array4<std::type_identity_t<T> const> const& phi,
                                       Array4<int const> const& msk, bool is_rz,
                                       GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
    {
        const int ii = rr*i;
        const int jj = rr*j;
        if (msk(ii,jj,0)) {
            T facx = T(1./6.)*dxinv[0]*dxinv[0];
            T facy = T(1./6.)*dxinv[1]*dxinv[1];

            auto is_fine = [&ccbx] (int ix, int iy) -> bool {
                return ccbx.contains(ix,iy,0);
            };

            T Df = T(0.0);

            const int ilo = amrex::max(ii-rr+1, ndbx.smallEnd(0));
            const int ihi = amrex::min(ii+rr-1, ndbx.bigEnd  (0));
            const int jlo = amrex::max(jj-rr+1, ndbx.smallEnd(1));
            const int jhi = amrex::min(jj+rr-1, ndbx.bigEnd  (1));

            for (int joff = jlo; joff <= jhi; ++joff) {
            for (int ioff = ilo; ioff <= ihi; ++ioff) {
                T scale = static_cast<T>((rr-std::abs(ii-ioff)) *
                                               (rr-std::abs(jj-joff)));
                if (ndbx.strictly_contains(ioff,joff,0)) {
                    Df += scale * (rhs(ioff,joff,0)-res(ioff,joff,0));
                } else {
                    Df += scale * mlndlap_sum_Ax
                        (is_fine, sig, ioff, joff, facx, facy, phi, is_rz);
                }
            }}

            f(i,j,0) = Df * (T(1.0)/static_cast<T>(rr*rr*rr*rr));
        } else {
            f(i,j,0) = T(0.0);
        }
    }

template <int rr, typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_Ax_fine_contrib (int i, int j, int /*k*/, Box const& ndbx, Box const& ccbx,
                              Array4<T> const& f, Array4<std::type_identity_t<T> const> const& res,
                              Array4<std::type_identity_t<T> const> const& rhs, Array4<std::type_identity_t<T> const> const& phi,
                              Array4<std::type_identity_t<T> const> const& sig, Array4<int const> const& msk,
                              bool is_rz,
                              GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    mlndlap_Ax_fine_contrib_doit<rr>
        ([&sig] (int ix, int iy, int) -> T const& { return sig(ix,iy,0); },
         i,j,ndbx,ccbx,f,res,rhs,phi,msk,is_rz,dxinv);
}

template <int rr, typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_Ax_fine_contrib_cs (int i, int j, int /*k*/, Box const& ndbx, Box const& ccbx,
                                 Array4<T> const& f, Array4<std::type_identity_t<T> const> const& res,
                                 Array4<std::type_identity_t<T> const> const& rhs, Array4<std::type_identity_t<T> const> const& phi,
                                 std::type_identity_t<T> sig, Array4<int const> const& msk,
                                 bool is_rz,
                                 GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    mlndlap_Ax_fine_contrib_doit<rr>
        ([=] (int, int, int) -> T { return sig; },
         i,j,ndbx,ccbx,f,res,rhs,phi,msk,is_rz,dxinv);
}

// Anisotropic-sigma variants of mlndlap_sum_Ax, mlndlap_Ax_fine_contrib,
// and mlndlap_res_cf_contrib used by the AMR-sync paths when m_use_mapped.
// Per-axis sigma: facx-terms use sx, facy-terms (including RZ fp/fm) use sy.
// Recover the scalar kernels when sx == sy.
    template <typename P, typename SX, typename SY, typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T mlndlap_sum_Ax_ha (P const& pred, SX const& sx, SY const& sy,
                            int i, int j, std::type_identity_t<T> facx, std::type_identity_t<T> facy,
                            Array4<T const> const& phi, bool is_rz) noexcept
    {
        T Ax = T(0.0);
        T fp = T(0.0), fm = T(0.0);
        if (is_rz) {
            fp = facy / static_cast<T>(2*i+1);
            fm = facy / static_cast<T>(2*i-1);
        }
        if (pred(i-1,j-1)) {
            Ax += facx*sx(i-1,j-1,0)*(T(2.)*(phi(i-1,j  ,0)-phi(i  ,j  ,0))
                                          +    (phi(i-1,j-1,0)-phi(i  ,j-1,0)))
                + facy*sy(i-1,j-1,0)*(T(2.)*(phi(i  ,j-1,0)-phi(i  ,j  ,0))
                                          +    (phi(i-1,j-1,0)-phi(i-1,j  ,0)))
                + fm  *sy(i-1,j-1,0)*          (phi(i  ,j-1,0)-phi(i  ,j  ,0));
        }
        if (pred(i,j-1)) {
            Ax += facx*sx(i,j-1,0)*(T(2.)*(phi(i+1,j  ,0)-phi(i  ,j  ,0))
                                        +    (phi(i+1,j-1,0)-phi(i  ,j-1,0)))
                + facy*sy(i,j-1,0)*(T(2.)*(phi(i  ,j-1,0)-phi(i  ,j  ,0))
                                        +    (phi(i+1,j-1,0)-phi(i+1,j  ,0)))
                - fp  *sy(i,j-1,0)*          (phi(i  ,j-1,0)-phi(i  ,j  ,0));
        }
        if (pred(i-1,j)) {
            Ax += facx*sx(i-1,j,0)*(T(2.)*(phi(i-1,j  ,0)-phi(i  ,j  ,0))
                                        +    (phi(i-1,j+1,0)-phi(i  ,j+1,0)))
                + facy*sy(i-1,j,0)*(T(2.)*(phi(i  ,j+1,0)-phi(i  ,j  ,0))
                                        +    (phi(i-1,j+1,0)-phi(i-1,j  ,0)))
                + fm  *sy(i-1,j,0)*          (phi(i  ,j+1,0)-phi(i  ,j  ,0));
        }
        if (pred(i,j)) {
            Ax += facx*sx(i,j,0)*(T(2.)*(phi(i+1,j  ,0)-phi(i  ,j  ,0))
                                     +     (phi(i+1,j+1,0)-phi(i  ,j+1,0)))
                + facy*sy(i,j,0)*(T(2.)*(phi(i  ,j+1,0)-phi(i  ,j  ,0))
                                     +     (phi(i+1,j+1,0)-phi(i+1,j  ,0)))
                - fp  *sy(i,j,0)*          (phi(i  ,j+1,0)-phi(i  ,j  ,0));
        }
        return Ax;
    }

    template <int rr, typename SX, typename SY, typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void mlndlap_Ax_fine_contrib_doit_ha (SX const& sx, SY const& sy,
                                          int i, int j, Box const& ndbx, Box const& ccbx,
                                          Array4<T> const& f, Array4<std::type_identity_t<T> const> const& res,
                                          Array4<std::type_identity_t<T> const> const& rhs,
                                          Array4<std::type_identity_t<T> const> const& phi,
                                          Array4<int const> const& msk, bool is_rz,
                                          GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
    {
        const int ii = rr*i;
        const int jj = rr*j;
        if (msk(ii,jj,0)) {
            T facx = T(1./6.)*dxinv[0]*dxinv[0];
            T facy = T(1./6.)*dxinv[1]*dxinv[1];

            auto is_fine = [&ccbx] (int ix, int iy) -> bool {
                return ccbx.contains(ix,iy,0);
            };

            T Df = T(0.0);

            const int ilo = amrex::max(ii-rr+1, ndbx.smallEnd(0));
            const int ihi = amrex::min(ii+rr-1, ndbx.bigEnd  (0));
            const int jlo = amrex::max(jj-rr+1, ndbx.smallEnd(1));
            const int jhi = amrex::min(jj+rr-1, ndbx.bigEnd  (1));

            for (int joff = jlo; joff <= jhi; ++joff) {
            for (int ioff = ilo; ioff <= ihi; ++ioff) {
                T scale = static_cast<T>((rr-std::abs(ii-ioff)) *
                                               (rr-std::abs(jj-joff)));
                if (ndbx.strictly_contains(ioff,joff,0)) {
                    Df += scale * (rhs(ioff,joff,0)-res(ioff,joff,0));
                } else {
                    Df += scale * mlndlap_sum_Ax_ha
                        (is_fine, sx, sy, ioff, joff, facx, facy, phi, is_rz);
                }
            }}

            f(i,j,0) = Df * (T(1.0)/static_cast<T>(rr*rr*rr*rr));
        } else {
            f(i,j,0) = T(0.0);
        }
    }

template <int rr, typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_Ax_fine_contrib_ha (int i, int j, int /*k*/, Box const& ndbx, Box const& ccbx,
                                 Array4<T> const& f, Array4<std::type_identity_t<T> const> const& res,
                                 Array4<std::type_identity_t<T> const> const& rhs, Array4<std::type_identity_t<T> const> const& phi,
                                 Array4<std::type_identity_t<T> const> const& sx, Array4<std::type_identity_t<T> const> const& sy,
                                 Array4<int const> const& msk,
                                 bool is_rz,
                                 GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    mlndlap_Ax_fine_contrib_doit_ha<rr>
        ([&sx] (int ix, int iy, int) -> T const& { return sx(ix,iy,0); },
         [&sy] (int ix, int iy, int) -> T const& { return sy(ix,iy,0); },
         i,j,ndbx,ccbx,f,res,rhs,phi,msk,is_rz,dxinv);
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_res_cf_contrib (int i, int j, int, Array4<T> const& res,
                             Array4<std::type_identity_t<T> const> const& phi, Array4<std::type_identity_t<T> const> const& rhs,
                             Array4<std::type_identity_t<T> const> const& sig, Array4<int const> const& dmsk,
                             Array4<int const> const& ndmsk, Array4<int const> const& ccmsk,
                             Array4<std::type_identity_t<T> const> const& fc,
                             GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv,
                             Box const& ccdom_p, Box const& nddom,
                             bool is_rz,
                             GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bclo,
                             GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bchi,
                             bool neumann_doubling) noexcept
{
    using namespace nodelap_detail;

    if (!dmsk(i,j,0) && ndmsk(i,j,0) == crse_fine_node) {
        T facx = T(1./6.)*dxinv[0]*dxinv[0];
        T facy = T(1./6.)*dxinv[1]*dxinv[1];

        T Ax = mlndlap_sum_Ax([&ccmsk, &ccdom_p] (int ix, int iy) -> bool
                                 {
                                     return ccdom_p.contains(ix,iy,0)
                                         && (ccmsk(ix,iy,0) == crse_cell);
                                 },
                                 [&sig] (int ix, int iy, int) -> T const&
                                 {
                                     return sig(ix,iy,0);
                                 },
                                 i, j, facx, facy, phi, is_rz);
        Ax += fc(i,j,0);
        T const ns = (neumann_doubling) ? neumann_scale<T>(i,j,nddom,bclo,bchi) : T(1.0);
        res(i,j,0) = rhs(i,j,0) - Ax*ns;
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_res_cf_contrib_ha (int i, int j, int, Array4<T> const& res,
                                Array4<std::type_identity_t<T> const> const& phi, Array4<std::type_identity_t<T> const> const& rhs,
                                Array4<std::type_identity_t<T> const> const& sx, Array4<std::type_identity_t<T> const> const& sy,
                                Array4<int const> const& dmsk,
                                Array4<int const> const& ndmsk, Array4<int const> const& ccmsk,
                                Array4<std::type_identity_t<T> const> const& fc,
                                GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv,
                                Box const& ccdom_p, Box const& nddom,
                                bool is_rz,
                                GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bclo,
                                GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bchi,
                                bool neumann_doubling) noexcept
{
    using namespace nodelap_detail;

    if (!dmsk(i,j,0) && ndmsk(i,j,0) == crse_fine_node) {
        T facx = T(1./6.)*dxinv[0]*dxinv[0];
        T facy = T(1./6.)*dxinv[1]*dxinv[1];

        T Ax = mlndlap_sum_Ax_ha([&ccmsk, &ccdom_p] (int ix, int iy) -> bool
                                    {
                                        return ccdom_p.contains(ix,iy,0)
                                            && (ccmsk(ix,iy,0) == crse_cell);
                                    },
                                    [&sx] (int ix, int iy, int) -> T const&
                                    {
                                        return sx(ix,iy,0);
                                    },
                                    [&sy] (int ix, int iy, int) -> T const&
                                    {
                                        return sy(ix,iy,0);
                                    },
                                    i, j, facx, facy, phi, is_rz);
        Ax += fc(i,j,0);
        T const ns = (neumann_doubling) ? neumann_scale<T>(i,j,nddom,bclo,bchi) : T(1.0);
        res(i,j,0) = rhs(i,j,0) - Ax*ns;
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_res_cf_contrib_cs (int i, int j, int, Array4<T> const& res,
                                Array4<std::type_identity_t<T> const> const& phi, Array4<std::type_identity_t<T> const> const& rhs,
                                std::type_identity_t<T> sig, Array4<int const> const& dmsk,
                                Array4<int const> const& ndmsk, Array4<int const> const& ccmsk,
                                Array4<std::type_identity_t<T> const> const& fc,
                                GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv,
                                Box const& ccdom_p, Box const& nddom,
                                bool is_rz,
                                GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bclo,
                                GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bchi,
                                bool neumann_doubling) noexcept
{
    using namespace nodelap_detail;

    if (!dmsk(i,j,0) && ndmsk(i,j,0) == crse_fine_node) {
        T facx = T(1./6.)*dxinv[0]*dxinv[0];
        T facy = T(1./6.)*dxinv[1]*dxinv[1];

        T Ax = mlndlap_sum_Ax([&ccmsk, &ccdom_p] (int ix, int iy) -> bool
                                 {
                                     return ccdom_p.contains(ix,iy,0)
                                         && (ccmsk(ix,iy,0) == crse_cell);
                                 },
                                 [=] (int, int, int) -> T
                                 {
                                     return sig;
                                 },
                                 i, j, facx, facy, phi, is_rz);
        Ax += fc(i,j,0);
        T const ns = (neumann_doubling) ? neumann_scale<T>(i,j,nddom,bclo,bchi) : T(1.0);
        res(i,j,0) = rhs(i,j,0) - Ax*ns;
    }
}

//
// RAP
//

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_set_stencil (Box const& bx, Array4<T> const& sten,
                          Array4<std::type_identity_t<T> const> const& sigma,
                          GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    T facx = T(1.0/6.0)*dxinv[0]*dxinv[0];
    T facy = T(1.0/6.0)*dxinv[1]*dxinv[1];
    T fxy = facx + facy;
    T f2xmy = T(2.0)*facx - facy;
    T fmx2y = T(2.0)*facy - facx;

    amrex::LoopConcurrent(bx, [=] (int i, int j, int k) noexcept
    {
        sten(i,j,k,1) = f2xmy*(sigma(i,j-1,k)+sigma(i,j,k));
        sten(i,j,k,2) = fmx2y*(sigma(i-1,j,k)+sigma(i,j,k));
        sten(i,j,k,3) = fxy*sigma(i,j,k);
    });
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_set_stencil_s0 (int i, int j, int k, Array4<T> const& sten) noexcept
{
    using namespace nodelap_detail;

    sten(i,j,k,0) = -(sten(i-1,j  ,k,1) + sten(i  ,j  ,k,1)
                    + sten(i  ,j-1,k,2) + sten(i  ,j  ,k,2)
                    + sten(i-1,j-1,k,3) + sten(i  ,j-1,k,3)
                    + sten(i-1,j  ,k,3) + sten(i  ,j  ,k,3));
    sten(i,j,k,4) = T(1.0) / (std::abs(sten(i-1,j  ,k,1)) + std::abs(sten(i,j  ,k,1))
                               + std::abs(sten(i  ,j-1,k,2)) + std::abs(sten(i,j  ,k,2))
                               + std::abs(sten(i-1,j-1,k,3)) + std::abs(sten(i,j-1,k,3))
                               + std::abs(sten(i-1,j  ,k,3)) + std::abs(sten(i,j  ,k,3)) + eps_v<T>);
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_stencil_rap (int i, int j, int, Array4<T> const& csten,
                          Array4<std::type_identity_t<T> const> const& fsten) noexcept
{
    using namespace nodelap_detail;

    constexpr int k = 0;

#if 0
    auto interp_from_mm_to = [&fsten] (int i_, int j_) -> T {
        T wxm = std::abs(fsten(i_-1,j_  ,0,1))/(std::abs(fsten(i_-1,j_-1,0,3))+std::abs(fsten(i_-1,j_  ,0,3))+eps_v<T>);
        T wym = std::abs(fsten(i_  ,j_-1,0,2))/(std::abs(fsten(i_-1,j_-1,0,3))+std::abs(fsten(i_  ,j_-1,0,3))+eps_v<T>);
        T wmm = std::abs(fsten(i_-1,j_-1,0,3)) * (T(1.) + wxm + wym);
        return wmm * fsten(i_,j_,0,4);
    };
#endif

    auto interp_from_mp_to = [&fsten] (int i_, int j_) -> T {
        T wxm = std::abs(fsten(i_-1,j_  ,0,1))/(std::abs(fsten(i_-1,j_-1,0,3))+std::abs(fsten(i_-1,j_  ,0,3))+eps_v<T>);
        T wyp = std::abs(fsten(i_  ,j_  ,0,2))/(std::abs(fsten(i_-1,j_  ,0,3))+std::abs(fsten(i_  ,j_  ,0,3))+eps_v<T>);
        T wmp = std::abs(fsten(i_-1,j_  ,0,3)) *(T(1.) + wxm + wyp);
        return wmp * fsten(i_,j_,0,4);
    };

    auto interp_from_pm_to = [&fsten] (int i_, int j_) -> T {
        T wxp = std::abs(fsten(i_  ,j_  ,0,1))/(std::abs(fsten(i_  ,j_-1,0,3))+std::abs(fsten(i_  ,j_  ,0,3))+eps_v<T>);
        T wym = std::abs(fsten(i_  ,j_-1,0,2))/(std::abs(fsten(i_-1,j_-1,0,3))+std::abs(fsten(i_  ,j_-1,0,3))+eps_v<T>);
        T wpm = std::abs(fsten(i_  ,j_-1,0,3)) * (T(1.) + wxp + wym);
        return wpm * fsten(i_,j_,0,4);
    };

    auto interp_from_pp_to = [&fsten] (int i_, int j_) -> T {
        T wxp = std::abs(fsten(i_  ,j_  ,0,1))/(std::abs(fsten(i_  ,j_-1,0,3))+std::abs(fsten(i_  ,j_  ,0,3))+eps_v<T>);
        T wyp = std::abs(fsten(i_  ,j_  ,0,2))/(std::abs(fsten(i_-1,j_  ,0,3))+std::abs(fsten(i_  ,j_  ,0,3))+eps_v<T>);
        T wpp = std::abs(fsten(i_  ,j_  ,0,3)) * (T(1.) + wxp + wyp);
        return wpp * fsten(i_,j_,0,4);
    };

    auto interp_from_m0_to = [&fsten] (int i_, int j_) -> T {
        return std::abs(fsten(i_-1,j_,0,1))/(std::abs(fsten(i_-1,j_,0,1))+std::abs(fsten(i_,j_,0,1))+eps_v<T>);
    };

    auto interp_from_p0_to = [&fsten] (int i_, int j_) -> T {
        return std::abs(fsten(i_,j_,0,1))/(std::abs(fsten(i_-1,j_,0,1))+std::abs(fsten(i_,j_,0,1))+eps_v<T>);
    };

    auto interp_from_0m_to = [&fsten] (int i_, int j_) -> T {
        return std::abs(fsten(i_,j_-1,0,2))/(std::abs(fsten(i_,j_-1,0,2))+std::abs(fsten(i_,j_,0,2))+eps_v<T>);
    };

    auto interp_from_0p_to = [&fsten] (int i_, int j_) -> T {
        return std::abs(fsten(i_,j_,0,2))/(std::abs(fsten(i_,j_-1,0,2))+std::abs(fsten(i_,j_,0,2))+eps_v<T>);
    };

#if 0
    auto Amm = [&fsten] (int i_, int j_) -> T {
        return fsten(i_-1,j_-1,0,3);
    };
#endif

    auto A0m = [&fsten] (int i_, int j_) -> T {
        return fsten(i_,j_-1,0,2);
    };

    auto Apm = [&fsten] (int i_, int j_) -> T {
        return fsten(i_,j_-1,0,3);
    };

    auto Am0 = [&fsten] (int i_, int j_) -> T {
        return fsten(i_-1,j_,0,1);
    };

    auto A00 = [&fsten] (int i_, int j_) -> T {
        return fsten(i_,j_,0,0);
    };

    auto Ap0 = [&fsten] (int i_, int j_) -> T {
        return fsten(i_,j_,0,1);
    };

    auto Amp = [&fsten] (int i_, int j_) -> T {
        return fsten(i_-1,j_,0,3);
    };

    auto A0p = [&fsten] (int i_, int j_) -> T {
        return fsten(i_,j_,0,2);
    };

    auto App = [&fsten] (int i_, int j_) -> T {
        return fsten(i_,j_,0,3);
    };

#if 0
    auto restrict_from_mm_to = [&fsten] (int ii_, int jj_) -> T {
        T wxp = std::abs(fsten(ii_-1,jj_-1,0,1))/(std::abs(fsten(ii_-1,jj_-2,0,3))+std::abs(fsten(ii_-1,jj_-1,0,3))+eps_v<T>);
        T wyp = std::abs(fsten(ii_-1,jj_-1,0,2))/(std::abs(fsten(ii_-2,jj_-1,0,3))+std::abs(fsten(ii_-1,jj_-1,0,3))+eps_v<T>);
        T wpp = std::abs(fsten(ii_-1,jj_-1,0,3))*(T(1.)+wxp+wyp);
        return wpp * fsten(ii_-1,jj_-1,0,4);
    };
#endif

    auto restrict_from_0m_to = [&fsten] (int ii_, int jj_) -> T {
        return std::abs(fsten(ii_,jj_-1,0,2))/(std::abs(fsten(ii_,jj_-2,0,2))+std::abs(fsten(ii_,jj_-1,0,2))+eps_v<T>);
    };

    auto restrict_from_pm_to = [&fsten] (int ii_, int jj_) -> T {
        T wxm = std::abs(fsten(ii_  ,jj_-1,0,1))/(std::abs(fsten(ii_,jj_-2,0,3))+std::abs(fsten(ii_  ,jj_-1,0,3))+eps_v<T>);
        T wyp = std::abs(fsten(ii_+1,jj_-1,0,2))/(std::abs(fsten(ii_,jj_-1,0,3))+std::abs(fsten(ii_+1,jj_-1,0,3))+eps_v<T>);
        T wmp = std::abs(fsten(ii_  ,jj_-1,0,3)) *(T(1.) + wxm + wyp);
        return wmp * fsten(ii_+1,jj_-1,0,4);
    };

    auto restrict_from_m0_to = [&fsten] (int ii_, int jj_) -> T {
        return std::abs(fsten(ii_-1,jj_,0,1))/(std::abs(fsten(ii_-2,jj_,0,1))+std::abs(fsten(ii_-1,jj_,0,1))+eps_v<T>);
    };

    auto restrict_from_p0_to = [&fsten] (int ii_, int jj_) -> T {
        return std::abs(fsten(ii_,jj_,0,1))/(std::abs(fsten(ii_,jj_,0,1))+std::abs(fsten(ii_+1,jj_,0,1))+eps_v<T>);
    };

    auto restrict_from_mp_to = [&fsten] (int ii_, int jj_) -> T {
        T wxp = std::abs(fsten(ii_-1,jj_+1,0,1))/(std::abs(fsten(ii_-1,jj_,0,3))+std::abs(fsten(ii_-1,jj_+1,0,3))+eps_v<T>);
        T wym = std::abs(fsten(ii_-1,jj_  ,0,2))/(std::abs(fsten(ii_-2,jj_,0,3))+std::abs(fsten(ii_-1,jj_  ,0,3))+eps_v<T>);
        T wpm = std::abs(fsten(ii_-1,jj_  ,0,3)) * (T(1.) + wxp + wym);
        return wpm * fsten(ii_-1,jj_+1,0,4);
    };

    auto restrict_from_0p_to = [&fsten] (int ii_, int jj_) -> T {
        return std::abs(fsten(ii_,jj_,0,2))/(std::abs(fsten(ii_,jj_,0,2))+std::abs(fsten(ii_,jj_+1,0,2))+eps_v<T>);
    };

    auto restrict_from_pp_to = [&fsten] (int ii_, int jj_) -> T {
        T wxm = std::abs(fsten(ii_  ,jj_+1,0,1))/(std::abs(fsten(ii_  ,jj_  ,0,3))+std::abs(fsten(ii_  ,jj_+1,0,3))+eps_v<T>);
        T wym = std::abs(fsten(ii_+1,jj_  ,0,2))/(std::abs(fsten(ii_  ,jj_  ,0,3))+std::abs(fsten(ii_+1,jj_  ,0,3))+eps_v<T>);
        T wmm = std::abs(fsten(ii_  ,jj_  ,0,3)) * (T(1.) + wxm + wym);
        return wmm * fsten(ii_+1,jj_+1,0,4);
    };

    int ii = 2*i;
    int jj = 2*j;
    Array2D<T,-1,1,-1,1> ap, p;

    // csten(i,j,k,1)
    p(-1,-1) = interp_from_pp_to(ii+1,jj-1);
    p( 0,-1) = interp_from_0p_to(ii+2,jj-1);
    p(-1, 0) = interp_from_p0_to(ii+1,jj  );
    p( 0, 0) = T(1.);
    p(-1, 1) = interp_from_pm_to(ii+1,jj+1);
    p( 0, 1) = interp_from_0m_to(ii+2,jj+1);

    ap(0,-1) = Ap0(ii,jj-1)*p(-1,-1) + App(ii,jj-1)*p(-1,0);
    ap(1,-1) = A00(ii+1,jj-1)*p(-1,-1) + Ap0(ii+1,jj-1)*p(0,-1)
        +      A0p(ii+1,jj-1)*p(-1,0) + App(ii+1,jj-1)*p(0,0);
    ap(0,0) = Apm(ii,jj)*p(-1,-1) + Ap0(ii,jj)*p(-1,0) + App(ii,jj)*p(-1,1);
    ap(1,0) = A0m(ii+1,jj)*p(-1,-1) + Apm(ii+1,jj)*p(0,-1)
        +     A00(ii+1,jj)*p(-1,0) + Ap0(ii+1,jj)*p(0,0)
        +     A0p(ii+1,jj)*p(-1,1) + App(ii+1,jj)*p(0,1);
    ap(0,1) = Apm(ii,jj+1)*p(-1,0) + Ap0(ii,jj+1)*p(-1,1);
    ap(1,1) = A0m(ii+1,jj+1)*p(-1,0) + Apm(ii+1,jj+1)*p(0,0)
        +     A00(ii+1,jj+1)*p(-1,1) + Ap0(ii+1,jj+1)*p(0,1);

    csten(i,j,k,1) = T(0.25)*(restrict_from_0m_to(ii,jj)*ap(0,-1)
                               + restrict_from_pm_to(ii,jj)*ap(1,-1)
                               + ap(0,0)
                               + restrict_from_p0_to(ii,jj)*ap(1,0)
                               + restrict_from_0p_to(ii,jj)*ap(0,1)
                               + restrict_from_pp_to(ii,jj)*ap(1,1));

    // csten(i,j,k,2)
    p(-1,-1) = interp_from_pp_to(ii-1,jj+1);
    p( 0,-1) = interp_from_0p_to(ii  ,jj+1);
    p( 1,-1) = interp_from_mp_to(ii+1,jj+1);
    p(-1, 0) = interp_from_p0_to(ii-1,jj+2);
    p( 0, 0) = T(1.);
    p( 1, 0) = interp_from_m0_to(ii+1,jj+2);

    ap(-1,0) = A0p(ii-1,jj)*p(-1,-1) + App(ii-1,jj)*p(0,-1);
    ap(0,0) = Amp(ii,jj)*p(-1,-1) + A0p(ii,jj)*p(0,-1) + App(ii,jj)*p(1,-1);
    ap(1,0) = Amp(ii+1,jj)*p(0,-1) + A0p(ii+1,jj)*p(1,-1);
    ap(-1,1) = A00(ii-1,jj+1)*p(-1,-1) + Ap0(ii-1,jj+1)*p(0,-1)
        +      A0p(ii-1,jj+1)*p(-1,0) + App(ii-1,jj+1)*p(0,0);
    ap(0,1) = Am0(ii,jj+1)*p(-1,-1) + A00(ii,jj+1)*p(0,-1) + Ap0(ii,jj+1)*p(1,-1)
        +     Amp(ii,jj+1)*p(-1,0) + A0p(ii,jj+1)*p(0,0) + App(ii,jj+1)*p(1,0);
    ap(1,1) = Am0(ii+1,jj+1)*p(0,-1) + A00(ii+1,jj+1)*p(1,-1)
        +     Amp(ii+1,jj+1)*p(0,0) + A0p(ii+1,jj+1)*p(1,0);

    csten(i,j,k,2) = T(0.25)*(restrict_from_m0_to(ii,jj)*ap(-1,0)
                            + ap(0,0)
                            + restrict_from_p0_to(ii,jj)*ap(1,0)
                            + restrict_from_mp_to(ii,jj)*ap(-1,1)
                            + restrict_from_0p_to(ii,jj)*ap(0,1)
                            + restrict_from_pp_to(ii,jj)*ap(1,1));

    // csten(i,j,k,3)
    p(-1,-1) = interp_from_pp_to(ii+1,jj+1);
    p( 0,-1) = interp_from_0p_to(ii+2,jj+1);
    p(-1, 0) = interp_from_p0_to(ii+1,jj+2);
    p( 0, 0) = T(1.);

    ap(0,0) = App(ii,jj)*p(-1,-1);
    ap(1,0) = A0p(ii+1,jj)*p(-1,-1) + App(ii+1,jj)*p(0,-1);
    ap(0,1) = Ap0(ii,jj+1)*p(-1,-1) + App(ii,jj+1)*p(-1,0);
    ap(1,1) = A00(ii+1,jj+1)*p(-1,-1) + Ap0(ii+1,jj+1)*p(0,-1)
        +     A0p(ii+1,jj+1)*p(-1,0) + App(ii+1,jj+1)*p(0,0);

    T cross1 = T(0.25)*(ap(0,0)
                           + restrict_from_p0_to(ii,jj)*ap(1,0)
                           + restrict_from_0p_to(ii,jj)*ap(0,1)
                           + restrict_from_pp_to(ii,jj)*ap(1,1));

    p(0,-1) = interp_from_0p_to(ii,jj+1);
    p(1,-1) = interp_from_mp_to(ii+1,jj+1);
    p(0, 0) = T(1.);
    p(1, 0) = interp_from_m0_to(ii+1,jj+2);

    ap(-1,0) = Amp(ii+1,jj)*p(0,-1) + A0p(ii+1,jj)*p(1,-1);
    ap( 0,0) = Amp(ii+2,jj)*p(1,-1);
    ap(-1,1) = Am0(ii+1,jj+1)*p(0,-1) + A00(ii+1,jj+1)*p(1,-1) + Amp(ii+1,jj+1)*p(0,0)
        + A0p(ii+1,jj+1)*p(1,0);
    ap( 0,1) = Am0(ii+2,jj+1)*p(1,-1) + Amp(ii+2,jj+1)*p(1,0);

    T cross2 = T(0.25)*(ap(0,0)
                           + restrict_from_m0_to(ii+2,jj)*ap(-1,0)
                           + restrict_from_mp_to(ii+2,jj)*ap(-1,1)
                           + restrict_from_0p_to(ii+2,jj)*ap( 0,1));

    csten(i,j,k,3) = T(0.5)*(cross1+cross2);
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T mlndlap_adotx_sten_doit (int i, int j, int k, Array4<T const> const& x,
                              Array4<std::type_identity_t<T> const> const& sten) noexcept
{
    return     x(i-1,j-1,k)*sten(i-1,j-1,k,3)
        +      x(i  ,j-1,k)*sten(i  ,j-1,k,2)
        +      x(i+1,j-1,k)*sten(i  ,j-1,k,3)
        +      x(i-1,j  ,k)*sten(i-1,j  ,k,1)
        +      x(i  ,j  ,k)*sten(i  ,j  ,k,0)
        +      x(i+1,j  ,k)*sten(i  ,j  ,k,1)
        +      x(i-1,j+1,k)*sten(i-1,j  ,k,3)
        +      x(i  ,j+1,k)*sten(i  ,j  ,k,2)
        +      x(i+1,j+1,k)*sten(i  ,j  ,k,3);
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T mlndlap_adotx_sten (int i, int j, int k, Array4<T const> const& x,
                         Array4<std::type_identity_t<T> const> const& sten, Array4<int const> const& msk) noexcept
{
    if (msk(i,j,k)) {
        return T(0.0);
    } else {
        return mlndlap_adotx_sten_doit(i,j,k,x,sten);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_gauss_seidel_sten (int i, int j, int k, Array4<T> const& sol,
                                Array4<std::type_identity_t<T> const> const& rhs,
                                Array4<std::type_identity_t<T> const> const& sten,
                                Array4<int const> const& msk) noexcept
{
    if (msk(i,j,k)) {
        sol(i,j,k) = T(0.0);
    } else if (sten(i,j,k,0) != T(0.0)) {
        T Ax = mlndlap_adotx_sten_doit<T>(i,j,k,sol,sten);
        sol(i,j,k) += (rhs(i,j,k) - Ax) / sten(i,j,k,0);
    }
}

template <typename T>
inline
void mlndlap_gauss_seidel_sten (Box const& bx, Array4<T> const& sol,
                                Array4<std::type_identity_t<T> const> const& rhs,
                                Array4<std::type_identity_t<T> const> const& sten,
                                Array4<int const> const& msk) noexcept
{
    AMREX_LOOP_3D(bx, i, j, k,
    {
        mlndlap_gauss_seidel_sten(i,j,k,sol,rhs,sten,msk);
    });
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_interpadd_rap (int i, int j, int, Array4<T> const& fine,
                            Array4<std::type_identity_t<T> const> const& crse, Array4<std::type_identity_t<T> const> const& sten,
                            Array4<int const> const& msk) noexcept
{
    using namespace nodelap_detail;

    if (!msk(i,j,0) && sten(i,j,0,0) != T(0.0)) {
        int ic = amrex::coarsen(i,2);
        int jc = amrex::coarsen(j,2);
        bool ieven = ic*2 == i;
        bool jeven = jc*2 == j;
        T fv;
        if (ieven && jeven) {
            fv = crse(ic,jc,0);
        } else if (ieven) {
            T wym = std::abs(sten(i,j-1,0,2));
            T wyp = std::abs(sten(i,j  ,0,2));
            fv = (wym*crse(ic,jc,0) + wyp*crse(ic,jc+1,0)) / (wym+wyp+eps_v<T>);
        } else if (jeven) {
            T wxm = std::abs(sten(i-1,j,0,1));
            T wxp = std::abs(sten(i  ,j,0,1));
            fv = (wxm*crse(ic,jc,0) + wxp*crse(ic+1,jc,0)) / (wxm+wxp+eps_v<T>);
        } else {
            T wxm = std::abs(sten(i-1,j  ,0,1)) /
                (std::abs(sten(i-1,j-1,0,3))+std::abs(sten(i-1,j  ,0,3))+eps_v<T>);
            T wxp = std::abs(sten(i  ,j  ,0,1)) /
                (std::abs(sten(i  ,j-1,0,3))+std::abs(sten(i  ,j  ,0,3))+eps_v<T>);
            T wym = std::abs(sten(i  ,j-1,0,2)) /
                (std::abs(sten(i-1,j-1,0,3))+std::abs(sten(i  ,j-1,0,3))+eps_v<T>);
            T wyp = std::abs(sten(i  ,j  ,0,2)) /
                (std::abs(sten(i-1,j  ,0,3))+std::abs(sten(i  ,j  ,0,3))+eps_v<T>);
            T wmm = std::abs(sten(i-1,j-1,0,3)) * (T(1.0) + wxm + wym);
            T wpm = std::abs(sten(i,j-1,0,3)) * (T(1.0) + wxp + wym);
            T wmp = std::abs(sten(i-1,j,0,3)) *(T(1.0) + wxm + wyp);
            T wpp = std::abs(sten(i,j,0,3)) * (T(1.0) + wxp + wyp);
            fv = (wmm*crse(ic,jc,0) + wpm*crse(ic+1,jc,0)
                  + wmp*crse(ic,jc+1,0) + wpp*crse(ic+1,jc+1,0))
                / (wmm+wpm+wmp+wpp+eps_v<T>);
        }

        fine(i,j,0) += fv;
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_restriction_rap (int i, int j, int /*k*/, Array4<T> const& crse,
                              Array4<std::type_identity_t<T> const> const& fine, Array4<std::type_identity_t<T> const> const& sten,
                              Array4<int const> const& msk) noexcept
{
    using namespace nodelap_detail;

    int ii = i*2;
    int jj = j*2;
    if (msk(ii,jj,0)) {
        crse(i,j,0) = T(0.0);
    } else {

        T cv = fine(ii,jj,0)
            + fine(ii-1,jj  ,0)*std::abs(sten(ii-1,jj  ,0,1))
            /                  (std::abs(sten(ii-2,jj  ,0,1))
                               +std::abs(sten(ii-1,jj  ,0,1))+eps_v<T>)
            + fine(ii+1,jj  ,0)*std::abs(sten(ii  ,jj  ,0,1))
            /                  (std::abs(sten(ii  ,jj  ,0,1))
                               +std::abs(sten(ii+1,jj  ,0,1))+eps_v<T>)
            + fine(ii  ,jj-1,0)*std::abs(sten(ii  ,jj-1,0,2))
            /                  (std::abs(sten(ii  ,jj-2,0,2))
                               +std::abs(sten(ii  ,jj-1,0,2))+eps_v<T>)
            + fine(ii  ,jj+1,0)*std::abs(sten(ii  ,jj  ,0,2))
            /                  (std::abs(sten(ii  ,jj  ,0,2))
                               +std::abs(sten(ii  ,jj+1,0,2))+eps_v<T>);

        T wxp = std::abs(sten(ii-1,jj-1,0,1))
            /     (std::abs(sten(ii-1,jj-2,0,3))
                  +std::abs(sten(ii-1,jj-1,0,3))+eps_v<T>);
        T wyp = std::abs(sten(ii-1,jj-1,0,2))
            /     (std::abs(sten(ii-2,jj-1,0,3))
                  +std::abs(sten(ii-1,jj-1,0,3))+eps_v<T>);
        T wpp = std::abs(sten(ii-1,jj-1,0,3))*(T(1.0) + wxp + wyp);
        cv +=           wpp*sten(ii-1,jj-1,0,4)*fine(ii-1,jj-1,0);

        T wxm = std::abs(sten(ii  ,jj-1,0,1))
            /     (std::abs(sten(ii  ,jj-2,0,3))
                  +std::abs(sten(ii  ,jj-1,0,3))+eps_v<T>);
        wyp      = std::abs(sten(ii+1,jj-1,0,2))
            /     (std::abs(sten(ii  ,jj-1,0,3))
                  +std::abs(sten(ii+1,jj-1,0,3))+eps_v<T>);
        T wmp = std::abs(sten(ii  ,jj-1,0,3))*(T(1.0) + wxm + wyp);
        cv +=           wmp*sten(ii+1,jj-1,0,4)*fine(ii+1,jj-1,0);

        wxp      = std::abs(sten(ii-1,jj+1,0,1))
            /     (std::abs(sten(ii-1,jj  ,0,3))
                  +std::abs(sten(ii-1,jj+1,0,3))+eps_v<T>);
        T wym = std::abs(sten(ii-1,jj  ,0,2))
            /     (std::abs(sten(ii-2,jj  ,0,3))
                  +std::abs(sten(ii-1,jj  ,0,3))+eps_v<T>);
        T wpm = std::abs(sten(ii-1,jj  ,0,3)) * (T(1.0) + wxp + wym);
        cv +=           wpm*sten(ii-1,jj+1,0,4)*fine(ii-1,jj+1,0);

        wxm      = std::abs(sten(ii  ,jj+1,0,1))
            /     (std::abs(sten(ii  ,jj  ,0,3))
                  +std::abs(sten(ii  ,jj+1,0,3))+eps_v<T>);
        wym      = std::abs(sten(ii+1,jj  ,0,2))
            /     (std::abs(sten(ii  ,jj  ,0,3))
                  +std::abs(sten(ii+1,jj  ,0,3))+eps_v<T>);
        T wmm = std::abs(sten(ii  ,jj  ,0,3)) * (T(1.0) + wxm + wym);
        cv +=           wmm*sten(ii+1,jj+1,0,4)*fine(ii+1,jj+1,0);

        crse(i,j,0) = cv * T(0.25);
    }
}

#ifdef AMREX_USE_EB

/// \cond DOXYGEN_IGNORE
namespace nodelap_detail {
    constexpr int i_S_x     = 0;
    constexpr int i_S_y     = 1;
    constexpr int i_S_x2    = 2;
    constexpr int i_S_y2    = 3;
    constexpr int i_S_xy    = 4;
    constexpr int n_Sintg   = 5;

    constexpr int i_B_x     = 0;
    constexpr int i_B_y     = 1;
    constexpr int i_B_xy    = 2;
    constexpr int n_Bintg   = 3;
}
/// \endcond

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_set_connection (int i, int j, int, Array4<T> const& conn,
                             Array4<std::type_identity_t<T> const> const& intg, Array4<std::type_identity_t<T> const> const& vol,
                             Array4<EBCellFlag const> const& flag) noexcept
{
    using namespace nodelap_detail;

    if (flag(i,j,0).isCovered()) {
        for (int n = 0; n < 6; ++n) { conn(i,j,0,n) = T(0.); }
    } else if (flag(i,j,0).isRegular() || vol(i,j,0) >= almostone_v<T>) {
        for (int n = 0; n < 6; ++n) { conn(i,j,0,n) = T(1.); }
    } else {
        // Note that these are normalized so that they equal 1 in the case of a regular cell

        conn(i,j,0,0) = T(3.)*(T(0.25)*vol(i,j,0) + intg(i,j,0,i_S_y2) - intg(i,j,0,i_S_y));
        conn(i,j,0,1) = T(6.)*(T(0.25)*vol(i,j,0) - intg(i,j,0,i_S_y2));
        conn(i,j,0,2) = T(3.)*(T(0.25)*vol(i,j,0) + intg(i,j,0,i_S_y2) + intg(i,j,0,i_S_y));

        conn(i,j,0,3) = T(3.)*(T(0.25)*vol(i,j,0) + intg(i,j,0,i_S_x2) - intg(i,j,0,i_S_x));
        conn(i,j,0,4) = T(6.)*(T(0.25)*vol(i,j,0) - intg(i,j,0,i_S_x2));
        conn(i,j,0,5) = T(3.)*(T(0.25)*vol(i,j,0) + intg(i,j,0,i_S_x2) + intg(i,j,0,i_S_x));
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_set_stencil_eb (int i, int j, int, Array4<T> const& sten,
                             Array4<std::type_identity_t<T> const> const& sig, Array4<std::type_identity_t<T> const> const& conn,
                             GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    T facx = T(1./6.)*dxinv[0]*dxinv[0];
    T facy = T(1./6.)*dxinv[1]*dxinv[1];

    sten(i,j,0,1) = T(2.)*facx*(sig(i,j-1,0)*conn(i,j-1,0,2)+sig(i,j,0)*conn(i,j,0,0))
                            -facy*(sig(i,j-1,0)*conn(i,j-1,0,4)+sig(i,j,0)*conn(i,j,0,4));
    sten(i,j,0,2) = T(2.)*facy*(sig(i-1,j,0)*conn(i-1,j,0,5)+sig(i,j,0)*conn(i,j,0,3))
                            -facx*(sig(i-1,j,0)*conn(i-1,j,0,1)+sig(i,j,0)*conn(i,j,0,1));
    sten(i,j,0,3) = (facx*conn(i,j,0,1)+facy*conn(i,j,0,4))*sig(i,j,0);
}


template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_divu_eb (int i, int j, int, Array4<T> const& rhs, Array4<std::type_identity_t<T> const> const& vel,
                      Array4<std::type_identity_t<T> const> const& vfrac, Array4<std::type_identity_t<T> const> const& intg,
                      Array4<int const> const& msk, GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv,
                      Box const& nodal_domain,
                      GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bclo,
                      GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bchi) noexcept
{
    T facx = T(0.5)*dxinv[0];
    T facy = T(0.5)*dxinv[1];

    const auto domlo = amrex::lbound(nodal_domain);
    const auto domhi = amrex::ubound(nodal_domain);

    if (!msk(i,j,0)) {

        T zero_ilo = T(1.0);
        T zero_ihi = T(1.0);
        T zero_jlo = T(1.0);
        T zero_jhi = T(1.0);

        // The nodal divergence operator should not see the tangential velocity
        //     at an inflow face
        if ((bclo[0] == LinOpBCType::Neumann || bclo[0] == LinOpBCType::inflow)
            && i == domlo.x)
        {
            zero_ilo = T(0.0);
        }
        if ((bchi[0] == LinOpBCType::Neumann || bchi[0] == LinOpBCType::inflow)
            && i == domhi.x)
        {
            zero_ihi = T(0.0);
        }
        if ((bclo[1] == LinOpBCType::Neumann || bclo[1] == LinOpBCType::inflow)
            && j == domlo.y)
        {
            zero_jlo = T(0.0);
        }
        if ((bchi[1] == LinOpBCType::Neumann || bchi[1] == LinOpBCType::inflow)
            && j == domhi.y)
        {
            zero_jhi = T(0.0);
        }

        rhs(i,j,0) = facx*(-vel(i-1,j-1,0,0)*(vfrac(i-1,j-1,0)+T(2.)*intg(i-1,j-1,0,1))*zero_jlo
                           +vel(i  ,j-1,0,0)*(vfrac(i  ,j-1,0)+T(2.)*intg(i  ,j-1,0,1))*zero_jlo
                           -vel(i-1,j  ,0,0)*(vfrac(i-1,j  ,0)-T(2.)*intg(i-1,j  ,0,1))*zero_jhi
                           +vel(i  ,j  ,0,0)*(vfrac(i  ,j  ,0)-T(2.)*intg(i  ,j  ,0,1))*zero_jhi)
                   + facy*(-vel(i-1,j-1,0,1)*(vfrac(i-1,j-1,0)+T(2.)*intg(i-1,j-1,0,0))*zero_ilo
                           -vel(i  ,j-1,0,1)*(vfrac(i  ,j-1,0)-T(2.)*intg(i  ,j-1,0,0))*zero_ihi
                           +vel(i-1,j  ,0,1)*(vfrac(i-1,j  ,0)+T(2.)*intg(i-1,j  ,0,0))*zero_ilo
                           +vel(i  ,j  ,0,1)*(vfrac(i  ,j  ,0)-T(2.)*intg(i  ,j  ,0,0))*zero_ihi);
    } else {
        rhs(i,j,0) = T(0.);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void add_eb_flow_contribution (int i, int j, int, Array4<T> const& rhs,
                      Array4<int const> const& msk, GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv,
                      Array4<std::type_identity_t<T> const> const& bareaarr,
                      Array4<std::type_identity_t<T> const> const& sintg,
                      Array4<std::type_identity_t<T> const> const& eb_vel_dot_n) noexcept
{
    using namespace nodelap_detail;

    T fac_eb = T(0.25) * dxinv[0];

    if (!msk(i,j,0)) {
        rhs(i,j,0) += fac_eb*(
           eb_vel_dot_n(i-1,j-1,0)*(         bareaarr(i-1,j-1,0)
                                      +T(2.)*sintg(i-1,j-1,0,i_B_x )
                                      +T(2.)*sintg(i-1,j-1,0,i_B_y )
                                      +T(4.)*sintg(i-1,j-1,0,i_B_xy))
          +eb_vel_dot_n(i  ,j-1,0)*(         bareaarr(i  ,j-1,0)
                                      -T(2.)*sintg(i  ,j-1,0,i_B_x )
                                      +T(2.)*sintg(i  ,j-1,0,i_B_y )
                                      -T(4.)*sintg(i  ,j-1,0,i_B_xy))
          +eb_vel_dot_n(i-1,j  ,0)*(         bareaarr(i-1,j  ,0)
                                      +T(2.)*sintg(i-1,j  ,0,i_B_x )
                                      -T(2.)*sintg(i-1,j  ,0,i_B_y )
                                      -T(4.)*sintg(i-1,j  ,0,i_B_xy))
          +eb_vel_dot_n(i  ,j  ,0)*(         bareaarr(i  ,j  ,0)
                                      -T(2.)*sintg(i  ,j  ,0,i_B_x )
                                      -T(2.)*sintg(i  ,j  ,0,i_B_y )
                                      +T(4.)*sintg(i  ,j  ,0,i_B_xy)));
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_mknewu_eb (int i, int j, int, Array4<T> const& u, Array4<std::type_identity_t<T> const> const& p,
                        Array4<std::type_identity_t<T> const> const& sig, Array4<std::type_identity_t<T> const> const& vfrac,
                        Array4<std::type_identity_t<T> const> const& intg, GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    T facx = T(0.5)*dxinv[0];
    T facy = T(0.5)*dxinv[1];
    if (vfrac(i,j,0) == T(0.)) {
        u(i,j,0,0) = u(i,j,0,1) = T(0.);
    } else {
        T dpdx = facx*(-p(i,j,0)+p(i+1,j,0)-p(i,j+1,0)+p(i+1,j+1,0));
        T dpdy = facy*(-p(i,j,0)-p(i+1,j,0)+p(i,j+1,0)+p(i+1,j+1,0));
        T dpp = (p(i,j,0)+p(i+1,j+1,0)-p(i+1,j,0)-p(i,j+1,0))/vfrac(i,j,0);
        u(i,j,0,0) -= sig(i,j,0)*(dpdx + dxinv[0]*intg(i,j,0,1)*dpp);
        u(i,j,0,1) -= sig(i,j,0)*(dpdy + dxinv[1]*intg(i,j,0,0)*dpp);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_mknewu_eb_c (int i, int j, int, Array4<T> const& u, Array4<std::type_identity_t<T> const> const& p,
                          std::type_identity_t<T> sig, Array4<std::type_identity_t<T> const> const& vfrac,
                          Array4<std::type_identity_t<T> const> const& intg, GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    T facx = T(0.5)*dxinv[0];
    T facy = T(0.5)*dxinv[1];
    if (vfrac(i,j,0) == T(0.)) {
        u(i,j,0,0) = u(i,j,0,1) = T(0.);
    } else {
        T dpdx = facx*(-p(i,j,0)+p(i+1,j,0)-p(i,j+1,0)+p(i+1,j+1,0));
        T dpdy = facy*(-p(i,j,0)-p(i+1,j,0)+p(i,j+1,0)+p(i+1,j+1,0));
        T dpp = (p(i,j,0)+p(i+1,j+1,0)-p(i+1,j,0)-p(i,j+1,0))/vfrac(i,j,0);
        u(i,j,0,0) -= sig*(dpdx + dxinv[0]*intg(i,j,0,1)*dpp);
        u(i,j,0,1) -= sig*(dpdy + dxinv[1]*intg(i,j,0,0)*dpp);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T mlndlap_rhcc_eb (int i, int j, int, Array4<T const> const& rhcc,
                      Array4<std::type_identity_t<T> const> const& vfrac, Array4<std::type_identity_t<T> const> const& intg,
                      Array4<int const> const& msk) noexcept
{
    using namespace nodelap_detail;

    if (!msk(i,j,0)) {
        return
            rhcc(i  ,j  ,0)*(T(0.25)*vfrac(i  ,j  ,0)
                             + T(0.5)*(-intg(i  ,j  ,0,i_S_x)-intg(i  ,j  ,0,i_S_y))
                             +             intg(i  ,j  ,0,i_S_xy)) +
            rhcc(i-1,j  ,0)*(T(0.25)*vfrac(i-1,j  ,0)
                             + T(0.5)*( intg(i-1,j  ,0,i_S_x)-intg(i-1,j  ,0,i_S_y))
                             -             intg(i-1,j  ,0,i_S_xy)) +
            rhcc(i-1,j-1,0)*(T(0.25)*vfrac(i-1,j-1,0)
                             + T(0.5)*( intg(i-1,j-1,0,i_S_x)+intg(i-1,j-1,0,i_S_y))
                             +             intg(i-1,j-1,0,i_S_xy)) +
            rhcc(i  ,j-1,0)*(T(0.25)*vfrac(i  ,j-1,0)
                             + T(0.5)*(-intg(i  ,j-1,0,i_S_x)+intg(i  ,j-1,0,i_S_y))
                             -             intg(i  ,j-1,0,i_S_xy));
    } else {
        return T(0.);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_set_integral (int i, int j, int, Array4<T> const& intg) noexcept
{
    using namespace nodelap_detail;

    intg(i,j,0,i_S_x ) = T(0.);
    intg(i,j,0,i_S_y ) = T(0.);
    intg(i,j,0,i_S_x2) = T(1./12.);
    intg(i,j,0,i_S_y2) = T(1./12.);
    intg(i,j,0,i_S_xy) = T(0.);
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_set_surface_integral (int i, int j, int, Array4<T> const& sintg) noexcept
{
    using namespace nodelap_detail;

    sintg(i,j,0,i_B_x ) = T(0.);
    sintg(i,j,0,i_B_y ) = T(0.);
    sintg(i,j,0,i_B_xy) = T(0.);
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_set_integral_eb (int i, int j, int, Array4<T> const& intg,
                              Array4<EBCellFlag const> const& flag, Array4<std::type_identity_t<T> const> const& vol,
                              Array4<std::type_identity_t<T> const> const& ax, Array4<std::type_identity_t<T> const> const& ay,
                              Array4<std::type_identity_t<T> const> const& bcen)
{
    using namespace nodelap_detail;

    if (flag(i,j,0).isCovered()) {
        intg(i,j,0,i_S_x ) = T(0.);
        intg(i,j,0,i_S_y ) = T(0.);
        intg(i,j,0,i_S_x2) = T(0.);
        intg(i,j,0,i_S_y2) = T(0.);
        intg(i,j,0,i_S_xy) = T(0.);
    } else if (flag(i,j,0).isRegular() || vol(i,j,0) >= almostone_v<T>) {
        intg(i,j,0,i_S_x ) = T(0.);
        intg(i,j,0,i_S_y ) = T(0.);
        intg(i,j,0,i_S_x2) = T(1./12.);
        intg(i,j,0,i_S_y2) = T(1./12.);
        intg(i,j,0,i_S_xy) = T(0.);
    } else {
        T axm = ax(i,j,0);
        T axp = ax(i+1,j,0);
        T aym = ay(i,j,0);
        T ayp = ay(i,j+1,0);

        T apnorm = std::sqrt((axm-axp)*(axm-axp) + (aym-ayp)*(aym-ayp));
        if (apnorm == T(0.)) {
            amrex::Abort("amrex_mlndlap_set_integral: we are in trouble");
        }

        T apnorminv = T(1.)/apnorm;
        T anrmx = (axm-axp) * apnorminv;  // pointing to the wall
        T anrmy = (aym-ayp) * apnorminv;

        T bcx = bcen(i,j,0,0);
        T bcy = bcen(i,j,0,1);

        T Sx, Sy, Sx2, Sy2, Sxy;
        if (std::abs(anrmx) <= almostzero_v<T>) {
            Sx = T(0.);
            Sx2 = T(1./24.)*(axm+axp);
            Sxy = T(0.);
        } else if (std::abs(anrmy) <= almostzero_v<T>) {
            Sx  = T(1./8.) *(axp-axm) + anrmx*T(0.5)*(bcx*bcx);
            Sx2 = T(1./24.)*(axp+axm) + anrmx*T(1./3.)*(bcx*bcx*bcx);
            Sxy = T(0.);
        } else {
            T xmin, xmax;
            if (anrmx > T(0.)) {
                xmin = T(-0.5) + amrex::min(aym,ayp);
                xmax = T(-0.5) + amrex::max(aym,ayp);
            } else {
                xmin = T(0.5) - amrex::max(aym,ayp);
                xmax = T(0.5) - amrex::min(aym,ayp);
            }
            T xmin3 = xmin*xmin*xmin;
            T xmin4 = xmin3*xmin;
            T xmax3 = xmax*xmax*xmax;
            T xmax4 = xmax3*xmax;
            Sx  = T(1./8.) *(axp-axm) + (anrmx/std::abs(anrmy))*T(1./6.) *(xmax3-xmin3);
            Sx2 = T(1./24.)*(axp+axm) + (anrmx/std::abs(anrmy))*T(1./12.)*(xmax4-xmin4);

            T kk = -anrmx/anrmy;
            T bb = bcy-kk*bcx;
            Sxy = T(1./8.)*kk*kk*(xmax4-xmin4) + T(1./3.)*kk*bb*(xmax3-xmin3)
                + (T(0.25)*bb*bb-T(1./16.))*(xmax*xmax-xmin*xmin);
            if (anrmy < T(0.)) { Sxy = -Sxy; }
        }

        if (std::abs(anrmy) <= almostzero_v<T>) {
            Sy = T(0.);
            Sy2 = T(1./24.)*(aym+ayp);
        } else if (std::abs(anrmx) <= almostzero_v<T>) {
            Sy  = T(1./8.) *(ayp-aym) + anrmy*T(0.5)*(bcy*bcy);
            Sy2 = T(1./24.)*(ayp+aym) + anrmy*T(1./3.)*(bcy*bcy*bcy);
        } else {
            T ymin, ymax;
            if (anrmy > T(0.)) {
                ymin = T(-0.5) + amrex::min(axm,axp);
                ymax = T(-0.5) + amrex::max(axm,axp);
            } else {
                ymin = T(0.5) - amrex::max(axm,axp);
                ymax = T(0.5) - amrex::min(axm,axp);
            }
            T ymin3 = ymin*ymin*ymin;
            T ymin4 = ymin3*ymin;
            T ymax3 = ymax*ymax*ymax;
            T ymax4 = ymax3*ymax;
            Sy  = T(1./8.) *(ayp-aym) + (anrmy/std::abs(anrmx))*T(1./6.) *(ymax3-ymin3);
            Sy2 = T(1./24.)*(ayp+aym) + (anrmy/std::abs(anrmx))*T(1./12.)*(ymax4-ymin4);
        }

        intg(i,j,0,i_S_x ) = Sx;
        intg(i,j,0,i_S_y ) = Sy;
        intg(i,j,0,i_S_x2) = Sx2;
        intg(i,j,0,i_S_y2) = Sy2;
        intg(i,j,0,i_S_xy) = Sxy;
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_set_surface_integral_eb (int i, int j, int, Array4<T> const& sintg,
                              Array4<EBCellFlag const> const& flag,
                              Array4<std::type_identity_t<T> const> const& bcen,
                              Array4<std::type_identity_t<T> const> const& barea,
                              Array4<std::type_identity_t<T> const> const& bnorm) noexcept
{
    using namespace nodelap_detail;

    if (flag(i,j,0).isCovered() || flag(i,j,0).isRegular()) {
        sintg(i,j,0,i_B_x ) = T(0.);
        sintg(i,j,0,i_B_y ) = T(0.);
        sintg(i,j,0,i_B_xy) = T(0.);
    } else {
        T bcx = bcen(i,j,0,0);
        T bcy = bcen(i,j,0,1);

        T btanx =  bnorm(i,j,0,1);
        T btany = -bnorm(i,j,0,0);

        T x0 =  bcx - T(0.5)*barea(i,j,0)*btanx;
        T x1 =  bcx + T(0.5)*barea(i,j,0)*btanx;

        T y0 =  bcy - T(0.5)*barea(i,j,0)*btany;
        T y1 =  bcy + T(0.5)*barea(i,j,0)*btany;

        T Bx  = barea(i,j,0)*T(0.5)*(x1 + x0);
        T By  = barea(i,j,0)*T(0.5)*(y1 + y0);
        T Bxy = barea(i,j,0)*(x0*y0 + (x0*(y1 - y0) + y0*(x1 - x0))/T(2.) + (x1 - x0)*(y1 - y0)/T(3.));

        sintg(i,j,0,i_B_x ) = Bx;
        sintg(i,j,0,i_B_y ) = By;
        sintg(i,j,0,i_B_xy) = Bxy;
    }
}

#endif


template <typename AlgInt, typename AlgGid>
void mlndlap_fillijmat_sten_cpu (Box const& ndbx,
                                 Array4<AlgGid const> const& gid,
                                 Array4<int const> const& lid,
                                 AlgInt* ncols, AlgInt* cols,
                                 Real* mat, // NOLINT(readability-non-const-parameter)
                                 Array4<Real const> const& sten) noexcept
{
    constexpr auto gidmax = std::numeric_limits<AlgGid>::max();
    AlgInt nelems = 0;
    amrex::LoopOnCpu(ndbx, [&] (int i, int j, int k) noexcept
    {
        if (lid(i,j,k) >= 0)
        {
            cols[nelems] = gid(i,j,k);
            mat[nelems] = sten(i,j,k,0);
            AlgInt nelems_old = nelems;
            ++nelems;

            if                (gid(i-1,j-1,k) < gidmax) {
                cols[nelems] = gid(i-1,j-1,k);
                mat[nelems] = sten(i-1,j-1,k,3);
                ++nelems;
            }

            if                (gid(i,j-1,k) < gidmax) {
                cols[nelems] = gid(i,j-1,k);
                mat[nelems] = sten(i,j-1,k,2);
                ++nelems;
            }

            if                (gid(i+1,j-1,k) < gidmax) {
                cols[nelems] = gid(i+1,j-1,k);
                mat[nelems] = sten(i  ,j-1,k,3);
                ++nelems;
            }

            if                (gid(i-1,j,k) < gidmax) {
                cols[nelems] = gid(i-1,j,k);
                mat[nelems] = sten(i-1,j,k,1);
                ++nelems;
            }

            if                (gid(i+1,j,k) < gidmax) {
                cols[nelems] = gid(i+1,j,k);
                mat[nelems] = sten(i  ,j,k,1);
                ++nelems;
            }

            if                (gid(i-1,j+1,k) < gidmax) {
                cols[nelems] = gid(i-1,j+1,k);
                mat[nelems] = sten(i-1,j  ,k,3);
                ++nelems;
            }

            if                (gid(i,j+1,k) < gidmax) {
                cols[nelems] = gid(i,j+1,k);
                mat[nelems] = sten(i,j  ,k,2);
                ++nelems;
            }

            if                (gid(i+1,j+1,k) < gidmax) {
                cols[nelems] = gid(i+1,j+1,k);
                mat[nelems] = sten(i  ,j  ,k,3);
                ++nelems;
            }

            ncols[lid(i,j,k)] = nelems - nelems_old;
        }
    });
}

template <typename AlgInt, typename AlgGid>
void mlndlap_fillijmat_aa_cpu (Box const& ndbx,
                               Array4<AlgGid const> const& gid,
                               Array4<int const> const& lid,
                               AlgInt* ncols, AlgInt* cols,
                               Real* mat, // NOLINT(readability-non-const-parameter)
                               Array4<Real const> const& sig,
                               GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                               Box const& ccdom, bool is_rz) noexcept
{
    Real facx = Real(1.0/6.0)*dxinv[0]*dxinv[0];
    Real facy = Real(1.0/6.0)*dxinv[1]*dxinv[1];
    Real fxy = facx + facy;
    Real f2xmy = Real(2.0)*facx - facy;
    Real fmx2y = Real(2.0)*facy - facx;

    // Note that ccdom has been grown at periodic boundaries.
    const Box& nddom = amrex::surroundingNodes(ccdom);

    constexpr auto gidmax = std::numeric_limits<AlgGid>::max();
    AlgInt nelems = 0;
    amrex::LoopOnCpu(ndbx, [&] (int i, int j, int k) noexcept
    {
        if (lid(i,j,k) >= 0)
        {
            Real fp, fm;
            if (is_rz) {
                fp = facy / static_cast<Real>(2*i+1);
                fm = facy / static_cast<Real>(2*i-1);
            } else {
                fp = fm = Real(0.0);
            }

            AlgInt nelems_old = nelems;
            cols[nelems_old] = gid(i,j,k);
            Real m0 = Real(0.);
            ++nelems;

            if (nddom.contains(i-1,j-1,k)) {
                Real tmp = fxy*sig(i-1,j-1,k);
                m0 -= tmp;
                if (               gid(i-1,j-1,k) < gidmax) {
                    cols[nelems] = gid(i-1,j-1,k);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i,j-1,k)) {
                Real tmp = Real(0.0);
                if (  ccdom.contains(i-1,j-1,k)) {
                    tmp +=       sig(i-1,j-1,k) * (fmx2y + fm);
                }
                if (  ccdom.contains(i,j-1,k)) {
                    tmp +=       sig(i,j-1,k) * (fmx2y - fp);
                }
                m0 -= tmp;
                if                (gid(i,j-1,k) < gidmax) {
                    cols[nelems] = gid(i,j-1,k);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i+1,j-1,k)) {
                Real tmp = fxy*sig(i  ,j-1,k);
                m0 -= tmp;
                if                (gid(i+1,j-1,k) < gidmax) {
                    cols[nelems] = gid(i+1,j-1,k);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i-1,j,k)) {
                Real tmp = Real(0.0);
                if (  ccdom.contains(i-1,j-1,k)) {
                    tmp += f2xmy*sig(i-1,j-1,k);
                }
                if (  ccdom.contains(i-1,j,k)) {
                    tmp += f2xmy*sig(i-1,j,k);
                }
                m0 -= tmp;
                if                (gid(i-1,j,k) < gidmax) {
                    cols[nelems] = gid(i-1,j,k);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i+1,j,k)) {
                Real tmp = Real(0.0);
                if (  ccdom.contains(i  ,j-1,k)) {
                    tmp += f2xmy*sig(i  ,j-1,k);
                }
                if (  ccdom.contains(i  ,j,k)) {
                    tmp += f2xmy*sig(i  ,j,k);
                }
                m0 -= tmp;
                if                (gid(i+1,j,k) < gidmax) {
                    cols[nelems] = gid(i+1,j,k);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i-1,j+1,k)) {
                Real tmp = fxy*sig(i-1,j  ,k);
                m0 -= tmp;
                if                (gid(i-1,j+1,k) < gidmax) {
                    cols[nelems] = gid(i-1,j+1,k);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i,j+1,k)) {
                Real tmp = Real(0.0);
                if (  ccdom.contains(i-1,j  ,k)) {
                    tmp +=       sig(i-1,j  ,k) * (fmx2y + fm);
                }
                if (  ccdom.contains(i,j  ,k)) {
                    tmp +=       sig(i,j  ,k) * (fmx2y - fp);
                }
                m0 -= tmp;
                if                (gid(i,j+1,k) < gidmax) {
                    cols[nelems] = gid(i,j+1,k);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i+1,j+1,k)) {
                Real tmp = fxy*sig(i  ,j  ,k);
                m0 -= tmp;
                if                (gid(i+1,j+1,k) < gidmax) {
                    cols[nelems] = gid(i+1,j+1,k);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            mat[nelems_old] = m0;
            ncols[lid(i,j,k)] = nelems - nelems_old;
        }
    });
}

template <typename AlgInt, typename AlgGid>
void mlndlap_fillijmat_ha_cpu (Box const& ndbx,
                               Array4<AlgGid const> const& gid,
                               Array4<int const> const& lid,
                               AlgInt* ncols, AlgInt* cols,
                               Real* mat, // NOLINT(readability-non-const-parameter)
                               Array4<Real const> const& sx,
                               Array4<Real const> const& sy,
                               GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                               Box const& ccdom, bool is_rz) noexcept
{
    Real facx = Real(1.0/6.0)*dxinv[0]*dxinv[0];
    Real facy = Real(1.0/6.0)*dxinv[1]*dxinv[1];

    // Note that ccdom has been grown at periodic boundaries.
    const Box& nddom = amrex::surroundingNodes(ccdom);

    constexpr auto gidmax = std::numeric_limits<AlgGid>::max();
    AlgInt nelems = 0;
    amrex::LoopOnCpu(ndbx, [&] (int i, int j, int k) noexcept
    {
        if (lid(i,j,k) >= 0)
        {
            Real fp, fm;
            if (is_rz) {
                fp = facy / static_cast<Real>(2*i+1);
                fm = facy / static_cast<Real>(2*i-1);
            } else {
                fp = fm = Real(0.0);
            }

            AlgInt nelems_old = nelems;
            cols[nelems_old] = gid(i,j,k);
            Real m0 = Real(0.);
            ++nelems;

            if (nddom.contains(i-1,j-1,k)) {
                Real tmp = facx*sx(i-1,j-1,k) + facy*sy(i-1,j-1,k);
                m0 -= tmp;
                if (               gid(i-1,j-1,k) < gidmax) {
                    cols[nelems] = gid(i-1,j-1,k);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i,j-1,k)) {
                Real tmp = Real(0.0);
                if (  ccdom.contains(i-1,j-1,k)) {
                    tmp +=        sy(i-1,j-1,k) * (facy * Real(2.0) + fm)
                                - sx(i-1,j-1,k) *  facx;
                }
                if (  ccdom.contains(i,j-1,k)) {
                    tmp +=        sy(i,j-1,k) * (facy * Real(2.0) - fp)
                                - sx(i,j-1,k) *  facx;
                }
                m0 -= tmp;
                if                (gid(i,j-1,k) < gidmax) {
                    cols[nelems] = gid(i,j-1,k);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i+1,j-1,k)) {
                Real tmp = facx*sx(i  ,j-1,k) + facy*sy(i  ,j-1,k);
                m0 -= tmp;
                if                (gid(i+1,j-1,k) < gidmax) {
                    cols[nelems] = gid(i+1,j-1,k);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i-1,j,k)) {
                Real tmp = Real(0.0);
                if (  ccdom.contains(i-1,j-1,k)) {
                    tmp +=        sx(i-1,j-1,k) * facx*Real(2.0)
                                - sy(i-1,j-1,k) * facy;
                }
                if (  ccdom.contains(i-1,j,k)) {
                    tmp +=        sx(i-1,j,k) * facx*Real(2.0)
                                - sy(i-1,j,k) * facy;
                }
                m0 -= tmp;
                if                (gid(i-1,j,k) < gidmax) {
                    cols[nelems] = gid(i-1,j,k);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i+1,j,k)) {
                Real tmp = Real(0.0);
                if (  ccdom.contains(i  ,j-1,k)) {
                    tmp +=        sx(i  ,j-1,k) * facx*Real(2.0)
                                - sy(i  ,j-1,k) * facy;
                }
                if (  ccdom.contains(i  ,j,k)) {
                    tmp +=        sx(i  ,j,k) * facx*Real(2.0)
                                - sy(i  ,j,k) * facy;
                }
                m0 -= tmp;
                if                (gid(i+1,j,k) < gidmax) {
                    cols[nelems] = gid(i+1,j,k);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i-1,j+1,k)) {
                Real tmp = facx*sx(i-1,j  ,k) + facy*sy(i-1,j  ,k);
                m0 -= tmp;
                if                (gid(i-1,j+1,k) < gidmax) {
                    cols[nelems] = gid(i-1,j+1,k);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i,j+1,k)) {
                Real tmp = Real(0.0);
                if (  ccdom.contains(i-1,j  ,k)) {
                    tmp +=        sy(i-1,j  ,k) * (facy*Real(2.0) + fm)
                                - sx(i-1,j  ,k) *  facx;
                }
                if (  ccdom.contains(i,j  ,k)) {
                    tmp +=       sy(i,j  ,k) * (facy*Real(2.0) - fp)
                               - sx(i,j  ,k) *  facx;
                }
                m0 -= tmp;
                if                (gid(i,j+1,k) < gidmax) {
                    cols[nelems] = gid(i,j+1,k);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i+1,j+1,k)) {
                Real tmp = facx*sx(i  ,j  ,k) + facy*sy(i  ,j  ,k);
                m0 -= tmp;
                if                (gid(i+1,j+1,k) < gidmax) {
                    cols[nelems] = gid(i+1,j+1,k);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            mat[nelems_old] = m0;
            ncols[lid(i,j,k)] = nelems - nelems_old;
        }
    });
}

template <typename AlgInt, typename AlgGid>
void mlndlap_fillijmat_cs_cpu (Box const& ndbx,
                               Array4<AlgGid const> const& gid,
                               Array4<int const> const& lid,
                               AlgInt* ncols, AlgInt* cols,
                               Real* mat, // NOLINT(readability-non-const-parameter)
                               Real sig,
                               GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                               Box const& ccdom, bool is_rz) noexcept
{
    Real facx = Real(1.0/6.0)*dxinv[0]*dxinv[0] * sig;
    Real facy = Real(1.0/6.0)*dxinv[1]*dxinv[1] * sig;
    Real fxy = facx + facy;
    Real f2xmy = Real(2.0)*facx - facy;
    Real fmx2y = Real(2.0)*facy - facx;

    // Note that ccdom has been grown at periodic boundaries.
    const Box& nddom = amrex::surroundingNodes(ccdom);

    constexpr auto gidmax = std::numeric_limits<AlgGid>::max();
    AlgInt nelems = 0;
    amrex::LoopOnCpu(ndbx, [&] (int i, int j, int k) noexcept
    {
        if (lid(i,j,k) >= 0)
        {
            Real fp, fm;
            if (is_rz) {
                fp = facy / static_cast<Real>(2*i+1);
                fm = facy / static_cast<Real>(2*i-1);
            } else {
                fp = fm = Real(0.0);
            }

            AlgInt nelems_old = nelems;
            cols[nelems_old] = gid(i,j,k);
            Real m0 = Real(0.);
            ++nelems;

            if (nddom.contains(i-1,j-1,k)) {
                Real tmp = fxy;
                m0 -= tmp;
                if (               gid(i-1,j-1,k) < gidmax) {
                    cols[nelems] = gid(i-1,j-1,k);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i,j-1,k)) {
                Real tmp = Real(0.0);
                if (  ccdom.contains(i-1,j-1,k)) {
                    tmp += fmx2y + fm;
                }
                if (  ccdom.contains(i,j-1,k)) {
                    tmp += fmx2y - fp;
                }
                m0 -= tmp;
                if                (gid(i,j-1,k) < gidmax) {
                    cols[nelems] = gid(i,j-1,k);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i+1,j-1,k)) {
                Real tmp = fxy;
                m0 -= tmp;
                if                (gid(i+1,j-1,k) < gidmax) {
                    cols[nelems] = gid(i+1,j-1,k);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i-1,j,k)) {
                Real tmp = Real(0.0);
                if (  ccdom.contains(i-1,j-1,k)) {
                    tmp += f2xmy;
                }
                if (  ccdom.contains(i-1,j,k)) {
                    tmp += f2xmy;
                }
                m0 -= tmp;
                if                (gid(i-1,j,k) < gidmax) {
                    cols[nelems] = gid(i-1,j,k);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i+1,j,k)) {
                Real tmp = Real(0.0);
                if (  ccdom.contains(i  ,j-1,k)) {
                    tmp += f2xmy;
                }
                if (  ccdom.contains(i  ,j,k)) {
                    tmp += f2xmy;
                }
                m0 -= tmp;
                if                (gid(i+1,j,k) < gidmax) {
                    cols[nelems] = gid(i+1,j,k);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i-1,j+1,k)) {
                Real tmp = fxy;
                m0 -= tmp;
                if                (gid(i-1,j+1,k) < gidmax) {
                    cols[nelems] = gid(i-1,j+1,k);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i,j+1,k)) {
                Real tmp = Real(0.0);
                if (  ccdom.contains(i-1,j  ,k)) {
                    tmp += fmx2y + fm;
                }
                if (  ccdom.contains(i,j  ,k)) {
                    tmp += fmx2y - fp;
                }
                m0 -= tmp;
                if                (gid(i,j+1,k) < gidmax) {
                    cols[nelems] = gid(i,j+1,k);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i+1,j+1,k)) {
                Real tmp = fxy;
                m0 -= tmp;
                if                (gid(i+1,j+1,k) < gidmax) {
                    cols[nelems] = gid(i+1,j+1,k);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            mat[nelems_old] = m0;
            ncols[lid(i,j,k)] = nelems - nelems_old;
        }
    });
}

#ifdef AMREX_USE_GPU

template <typename AlgInt, typename AlgGid>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void mlndlap_fillijmat_sten_gpu (const int ps, const int i, const int j, const int k,
                                 const int offset,
                                 Array4<AlgGid const> const& gid,
                                 Array4<int const> const& lid,
                                 AlgInt* ncols, AlgInt* cols,
                                 Real* mat, // NOLINT(readability-non-const-parameter)
                                 Array4<Real const> const& sten) noexcept
{
    if (lid(i,j,k) >= 0)
    {
        constexpr auto gidmax = std::numeric_limits<AlgGid>::max();
        int nelems = 0;

        if (offset == 1 || offset == 0) {
            if                (gid(i-1,j-1,k) < gidmax) {
                if (offset != 0) {
                    cols[ps] = gid(i-1,j-1,k);
                    mat[ps] = sten(i-1,j-1,k,3);
                }
                ++nelems;
            }
            if (offset != 0) { return; }
        }

        if (offset == 2 || offset == 0) {
            if                (gid(i,j-1,k) < gidmax) {
                if (offset != 0) {
                    cols[ps] = gid(i,j-1,k);
                    mat[ps] = sten(i,j-1,k,2);
                }
                ++nelems;
            }
            if (offset != 0) { return; }
        }

        if (offset == 3 || offset == 0) {
            if                (gid(i+1,j-1,k) < gidmax) {
                if (offset != 0) {
                    cols[ps] = gid(i+1,j-1,k);
                    mat[ps] = sten(i  ,j-1,k,3);
                }
                ++nelems;
            }
            if (offset != 0) { return; }
        }

        if (offset == 4 || offset == 0) {
            if                (gid(i-1,j,k) < gidmax) {
                if (offset != 0) {
                    cols[ps] = gid(i-1,j,k);
                    mat[ps] = sten(i-1,j,k,1);
                }
                ++nelems;
            }
            if (offset != 0) { return; }
        }

        if (offset == 5 || offset == 0) {
            if                (gid(i+1,j,k) < gidmax) {
                if (offset != 0) {
                    cols[ps] = gid(i+1,j,k);
                    mat[ps] = sten(i  ,j,k,1);
                }
                ++nelems;
            }
            if (offset != 0) { return; }
        }

        if (offset == 6 || offset == 0) {
            if                (gid(i-1,j+1,k) < gidmax) {
                if (offset != 0) {
                    cols[ps] = gid(i-1,j+1,k);
                    mat[ps] = sten(i-1,j  ,k,3);
                }
                ++nelems;
            }
            if (offset != 0) { return; }
        }

        if (offset == 7 || offset == 0) {
            if                (gid(i,j+1,k) < gidmax) {
                if (offset != 0) {
                    cols[ps] = gid(i,j+1,k);
                    mat[ps] = sten(i,j  ,k,2);
                }
                ++nelems;
            }
            if (offset != 0) { return; }
        }

        if (offset == 8 || offset == 0) {
            if                (gid(i+1,j+1,k) < gidmax) {
                if (offset != 0) {
                    cols[ps] = gid(i+1,j+1,k);
                    mat[ps] = sten(i  ,j  ,k,3);
                }
                ++nelems;
            }
            if (offset != 0) { return; }
        }

        // Only offset == 0 could get this far.
        cols[ps] = gid(i,j,k);
        mat[ps] = sten(i,j,k,0);
        ncols[lid(i,j,k)] = nelems+1;
    }
}

template <typename AlgInt, typename AlgGid>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void mlndlap_fillijmat_aa_gpu (const int ps, const int i, const int j, const int k,
                               const int offset,
                               Box const& /*ndbx*/, Array4<AlgGid const> const& gid,
                               Array4<int const> const& lid,
                               AlgInt* ncols, AlgInt* cols,
                               Real* mat, // NOLINT(readability-non-const-parameter)
                               Array4<Real const> const& sig,
                               GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                               Box const& ccdom, bool is_rz) noexcept
{
    if (lid(i,j,k) >= 0)
    {
        Real facx = Real(1.0/6.0)*dxinv[0]*dxinv[0];
        Real facy = Real(1.0/6.0)*dxinv[1]*dxinv[1];
        Real fxy = facx + facy;
        Real f2xmy = Real(2.0)*facx - facy;
        Real fmx2y = Real(2.0)*facy - facx;

        Real fp, fm;
        if (is_rz) {
            fp = facy / static_cast<Real>(2*i+1);
            fm = facy / static_cast<Real>(2*i-1);
        } else {
            fp = fm = Real(0.0);
        }

        const Box& nddom = amrex::surroundingNodes(ccdom);

        constexpr auto gidmax = std::numeric_limits<AlgGid>::max();
        int nelems = 0;
        Real m0 = Real(0.);

        if (offset == 1 || offset == 0) {
            if (nddom.contains(i-1,j-1,k)) {
                Real tmp = fxy*sig(i-1,j-1,k);
                m0 -= tmp;
                if                    (gid(i-1,j-1,k) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i-1,j-1,k);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 2 || offset == 0) {
            if (nddom.contains(i,j-1,k)) {
                Real tmp = Real(0.0);
                if (  ccdom.contains(i-1,j-1,k)) {
                    tmp +=       sig(i-1,j-1,k) * (fmx2y + fm);
                }
                if (  ccdom.contains(i,j-1,k)) {
                    tmp +=       sig(i,j-1,k) * (fmx2y - fp);
                }
                m0 -= tmp;
                if                (gid(i,j-1,k) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i,j-1,k);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 3 || offset == 0) {
            if (nddom.contains(i+1,j-1,k)) {
                Real tmp = fxy*sig(i  ,j-1,k);
                m0 -= tmp;
                if                (gid(i+1,j-1,k) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i+1,j-1,k);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 4 || offset == 0) {
            if (nddom.contains(i-1,j,k)) {
                Real tmp = Real(0.0);
                if (  ccdom.contains(i-1,j-1,k)) {
                    tmp += f2xmy*sig(i-1,j-1,k);
                }
                if (  ccdom.contains(i-1,j,k)) {
                    tmp += f2xmy*sig(i-1,j,k);
                }
                m0 -= tmp;
                if                (gid(i-1,j,k) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i-1,j,k);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 5 || offset == 0) {
            if (nddom.contains(i+1,j,k)) {
                Real tmp = Real(0.0);
                if (  ccdom.contains(i  ,j-1,k)) {
                    tmp += f2xmy*sig(i  ,j-1,k);
                }
                if (  ccdom.contains(i  ,j,k)) {
                    tmp += f2xmy*sig(i  ,j,k);
                }
                m0 -= tmp;
                if                (gid(i+1,j,k) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i+1,j,k);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 6 || offset == 0) {
            if (nddom.contains(i-1,j+1,k)) {
                Real tmp = fxy*sig(i-1,j  ,k);
                m0 -= tmp;
                if                (gid(i-1,j+1,k) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i-1,j+1,k);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 7 || offset == 0) {
            if (nddom.contains(i,j+1,k)) {
                Real tmp = Real(0.0);
                if (  ccdom.contains(i-1,j  ,k)) {
                    tmp +=       sig(i-1,j  ,k) * (fmx2y + fm);
                }
                if (  ccdom.contains(i,j  ,k)) {
                    tmp +=       sig(i,j  ,k) * (fmx2y - fp);
                }
                m0 -= tmp;
                if                (gid(i,j+1,k) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i,j+1,k);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 8 || offset == 0) {
            if (nddom.contains(i+1,j+1,k)) {
                Real tmp = fxy*sig(i  ,j  ,k);
                m0 -= tmp;
                if                (gid(i+1,j+1,k) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i+1,j+1,k);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        // Only offset == 0 could get this far.
        cols[ps] = gid(i,j,k);
        mat[ps] = m0;
        ncols[lid(i,j,k)] = nelems+1;
    }
}

template <typename AlgInt, typename AlgGid>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void mlndlap_fillijmat_ha_gpu (const int ps, const int i, const int j, const int k,
                               const int offset,
                               Box const& /*ndbx*/, Array4<AlgGid const> const& gid,
                               Array4<int const> const& lid,
                               AlgInt* ncols, AlgInt* cols,
                               Real* mat, // NOLINT(readability-non-const-parameter)
                               Array4<Real const> const& sx,
                               Array4<Real const> const& sy,
                               GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                               Box const& ccdom, bool is_rz) noexcept
{
    if (lid(i,j,k) >= 0)
    {
        Real facx = Real(1.0/6.0)*dxinv[0]*dxinv[0];
        Real facy = Real(1.0/6.0)*dxinv[1]*dxinv[1];

        Real fp, fm;
        if (is_rz) {
            fp = facy / static_cast<Real>(2*i+1);
            fm = facy / static_cast<Real>(2*i-1);
        } else {
            fp = fm = Real(0.0);
        }

        const Box& nddom = amrex::surroundingNodes(ccdom);

        constexpr auto gidmax = std::numeric_limits<AlgGid>::max();
        int nelems = 0;
        Real m0 = Real(0.);

        if (offset == 1 || offset == 0) {
            if (nddom.contains(i-1,j-1,k)) {
                Real tmp = facx*sx(i-1,j-1,k) + facy*sy(i-1,j-1,k);
                m0 -= tmp;
                if (               gid(i-1,j-1,k) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i-1,j-1,k);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 2 || offset == 0) {
            if (nddom.contains(i,j-1,k)) {
                Real tmp = Real(0.0);
                if (  ccdom.contains(i-1,j-1,k)) {
                    tmp +=        sy(i-1,j-1,k) * (facy * Real(2.0) + fm)
                                - sx(i-1,j-1,k) *  facx;
                }
                if (  ccdom.contains(i,j-1,k)) {
                    tmp +=        sy(i,j-1,k) * (facy * Real(2.0) - fp)
                                - sx(i,j-1,k) *  facx;
                }
                m0 -= tmp;
                if                (gid(i,j-1,k) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i,j-1,k);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 3 || offset == 0) {
            if (nddom.contains(i+1,j-1,k)) {
                Real tmp = facx*sx(i  ,j-1,k) + facy*sy(i  ,j-1,k);
                m0 -= tmp;
                if                (gid(i+1,j-1,k) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i+1,j-1,k);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 4 || offset == 0) {
            if (nddom.contains(i-1,j,k)) {
                Real tmp = Real(0.0);
                if (  ccdom.contains(i-1,j-1,k)) {
                    tmp +=        sx(i-1,j-1,k) * facx*Real(2.0)
                                - sy(i-1,j-1,k) * facy;
                }
                if (  ccdom.contains(i-1,j,k)) {
                    tmp +=        sx(i-1,j,k) * facx*Real(2.0)
                                - sy(i-1,j,k) * facy;
                }
                m0 -= tmp;
                if                (gid(i-1,j,k) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i-1,j,k);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 5 || offset == 0) {
            if (nddom.contains(i+1,j,k)) {
                Real tmp = Real(0.0);
                if (  ccdom.contains(i  ,j-1,k)) {
                    tmp +=        sx(i  ,j-1,k) * facx*Real(2.0)
                                - sy(i  ,j-1,k) * facy;
                }
                if (  ccdom.contains(i  ,j,k)) {
                    tmp +=        sx(i  ,j,k) * facx*Real(2.0)
                                - sy(i  ,j,k) * facy;
                }
                m0 -= tmp;
                if                (gid(i+1,j,k) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i+1,j,k);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 6 || offset == 0) {
            if (nddom.contains(i-1,j+1,k)) {
                Real tmp = facx*sx(i-1,j  ,k) + facy*sy(i-1,j  ,k);
                m0 -= tmp;
                if                (gid(i-1,j+1,k) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i-1,j+1,k);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 7 || offset == 0) {
            if (nddom.contains(i,j+1,k)) {
                Real tmp = Real(0.0);
                if (  ccdom.contains(i-1,j  ,k)) {
                    tmp +=        sy(i-1,j  ,k) * (facy*Real(2.0) + fm)
                                - sx(i-1,j  ,k) *  facx;
                }
                if (  ccdom.contains(i,j  ,k)) {
                    tmp +=       sy(i,j  ,k) * (facy*Real(2.0) - fp)
                               - sx(i,j  ,k) *  facx;
                }
                m0 -= tmp;
                if                (gid(i,j+1,k) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i,j+1,k);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 8 || offset == 0) {
            if (nddom.contains(i+1,j+1,k)) {
                Real tmp = facx*sx(i  ,j  ,k) + facy*sy(i  ,j  ,k);
                m0 -= tmp;
                if                (gid(i+1,j+1,k) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i+1,j+1,k);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        // Only offset == 0 could get this far.
        cols[ps] = gid(i,j,k);
        mat[ps] = m0;
        ncols[lid(i,j,k)] = nelems+1;
    }
}

template <typename AlgInt, typename AlgGid>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void mlndlap_fillijmat_cs_gpu (const int ps, const int i, const int j, const int k,
                               const int offset,
                               Box const& /*ndbx*/, Array4<AlgGid const> const& gid,
                               Array4<int const> const& lid,
                               AlgInt* ncols, AlgInt* cols,
                               Real* mat, // NOLINT(readability-non-const-parameter)
                               Real sig, GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                               Box const& ccdom, bool is_rz) noexcept
{
    if (lid(i,j,k) >= 0)
    {
        Real facx = Real(1.0/6.0)*dxinv[0]*dxinv[0] * sig;
        Real facy = Real(1.0/6.0)*dxinv[1]*dxinv[1] * sig;
        Real fxy = facx + facy;
        Real f2xmy = Real(2.0)*facx - facy;
        Real fmx2y = Real(2.0)*facy - facx;

        Real fp, fm;
        if (is_rz) {
            fp = facy / static_cast<Real>(2*i+1);
            fm = facy / static_cast<Real>(2*i-1);
        } else {
            fp = fm = Real(0.0);
        }

        // Note that nddom has been grown at periodic boundaries.
        const Box& nddom = amrex::surroundingNodes(ccdom);

        constexpr auto gidmax = std::numeric_limits<AlgGid>::max();
        int nelems = 0;
        Real m0 = Real(0.);

        if (offset == 1 || offset == 0) {
            if (nddom.contains(i-1,j-1,k)) {
                Real tmp = fxy;
                m0 -= tmp;
                if (               gid(i-1,j-1,k) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i-1,j-1,k);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 2 || offset == 0) {
            if (nddom.contains(i,j-1,k)) {
                Real tmp = Real(0.0);
                if (  ccdom.contains(i-1,j-1,k)) {
                    tmp += fmx2y + fm;
                }
                if (  ccdom.contains(i,j-1,k)) {
                    tmp += fmx2y - fp;
                }
                m0 -= tmp;
                if                (gid(i,j-1,k) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i,j-1,k);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 3 || offset == 0) {
            if (nddom.contains(i+1,j-1,k)) {
                Real tmp = fxy;
                m0 -= tmp;
                if                (gid(i+1,j-1,k) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i+1,j-1,k);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 4 || offset == 0) {
            if (nddom.contains(i-1,j,k)) {
                Real tmp = Real(0.0);
                if (  ccdom.contains(i-1,j-1,k)) {
                    tmp += f2xmy;
                }
                if (  ccdom.contains(i-1,j,k)) {
                    tmp += f2xmy;
                }
                m0 -= tmp;
                if                (gid(i-1,j,k) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i-1,j,k);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 5 || offset == 0) {
            if (nddom.contains(i+1,j,k)) {
                Real tmp = Real(0.0);
                if (  ccdom.contains(i  ,j-1,k)) {
                    tmp += f2xmy;
                }
                if (  ccdom.contains(i  ,j,k)) {
                    tmp += f2xmy;
                }
                m0 -= tmp;
                if                (gid(i+1,j,k) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i+1,j,k);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 6 || offset == 0) {
            if (nddom.contains(i-1,j+1,k)) {
                Real tmp = fxy;
                m0 -= tmp;
                if                (gid(i-1,j+1,k) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i-1,j+1,k);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 7 || offset == 0) {
            if (nddom.contains(i,j+1,k)) {
                Real tmp = Real(0.0);
                if (  ccdom.contains(i-1,j  ,k)) {
                    tmp += fmx2y + fm;
                }
                if (  ccdom.contains(i,j  ,k)) {
                    tmp += fmx2y - fp;
                }
                m0 -= tmp;
                if                (gid(i,j+1,k) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i,j+1,k);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 8 || offset == 0) {
            if (nddom.contains(i+1,j+1,k)) {
                Real tmp = fxy;
                m0 -= tmp;
                if                (gid(i+1,j+1,k) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i+1,j+1,k);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        // Only offset == 0 could get this far.
        cols[ps] = gid(i,j,k);
        mat[ps] = m0;
        ncols[lid(i,j,k)] = nelems+1;
    }
}

#endif


template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void mlndlap_gscolor_ha (int i, int j, int k, Array4<T> const& sol,
                         Array4<std::type_identity_t<T> const> const& rhs, Array4<std::type_identity_t<T> const> const& sx,
                         Array4<std::type_identity_t<T> const> const& sy, Array4<int const> const& msk,
                         GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv,
                         bool is_rz) noexcept
{
    if (msk(i,j,k)) {
        sol(i,j,k) = T(0.0);
    } else {
        T facx = T(1.0/6.0)*dxinv[0]*dxinv[0];
        T facy = T(1.0/6.0)*dxinv[1]*dxinv[1];

        T s0 = T(-2.0)*(facx*(sx(i-1,j-1,k)+sx(i,j-1,k)+sx(i-1,j,k)+sx(i,j,k))
                             +facy*(sy(i-1,j-1,k)+sy(i,j-1,k)+sy(i-1,j,k)+sy(i,j,k)));

        T Ax = sol(i-1,j-1,k)*(facx*sx(i-1,j-1,k)+facy*sy(i-1,j-1,k))
                + sol(i+1,j-1,k)*(facx*sx(i  ,j-1,k)+facy*sy(i  ,j-1,k))
                + sol(i-1,j+1,k)*(facx*sx(i-1,j  ,k)+facy*sy(i-1,j  ,k))
                + sol(i+1,j+1,k)*(facx*sx(i  ,j  ,k)+facy*sy(i  ,j  ,k))
                + sol(i-1,j,k)*(T(2.0)*facx*(sx(i-1,j-1,k)+sx(i-1,j,k))
                                    -     facy*(sy(i-1,j-1,k)+sy(i-1,j,k)))
                + sol(i+1,j,k)*(T(2.0)*facx*(sx(i  ,j-1,k)+sx(i  ,j,k))
                                    -     facy*(sy(i  ,j-1,k)+sy(i  ,j,k)))
                + sol(i,j-1,k)*(   -facx*(sx(i-1,j-1,k)+sx(i,j-1,k))
                         +T(2.0)*facy*(sy(i-1,j-1,k)+sy(i,j-1,k)))
                + sol(i,j+1,k)*(   -facx*(sx(i-1,j  ,k)+sx(i,j  ,k))
                         +T(2.0)*facy*(sy(i-1,j  ,k)+sy(i,j  ,k)))
                + sol(i,j,k)*s0;

        if (is_rz) {
            T fp = facy / static_cast<T>(2*i+1);
            T fm = facy / static_cast<T>(2*i-1);
            T frzlo = fm*sy(i-1,j-1,k)-fp*sy(i,j-1,k);
            T frzhi = fm*sy(i-1,j  ,k)-fp*sy(i,j  ,k);
            s0 += - frzhi - frzlo;
            Ax += frzhi*(sol(i,j+1,k)-sol(i,j,k))
                + frzlo*(sol(i,j-1,k)-sol(i,j,k));
        }

        sol(i,j,k) += (rhs(i,j,k) - Ax) / s0;
    }
}

template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void mlndlap_gscolor_aa (int i, int j, int k, Array4<T> const& sol,
                         Array4<std::type_identity_t<T> const> const& rhs, Array4<std::type_identity_t<T> const> const& sig,
                         Array4<int const> const& msk,
                         GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv,
                         bool is_rz) noexcept
{
    if (msk(i,j,k)) {
        sol(i,j,k) = T(0.0);
    } else {
        T facx = T(1.0/6.0)*dxinv[0]*dxinv[0];
        T facy = T(1.0/6.0)*dxinv[1]*dxinv[1];
        T fxy = facx + facy;
        T f2xmy = T(2.0)*facx - facy;
        T fmx2y = T(2.0)*facy - facx;

        T s0 = (-T(2.0))*fxy*(sig(i-1,j-1,k)+sig(i,j-1,k)+sig(i-1,j,k)+sig(i,j,k));
        T Ax =   sol(i-1,j-1,k)*fxy*sig(i-1,j-1,k)
                  + sol(i+1,j-1,k)*fxy*sig(i  ,j-1,k)
                  + sol(i-1,j+1,k)*fxy*sig(i-1,j  ,k)
                  + sol(i+1,j+1,k)*fxy*sig(i  ,j  ,k)
                  + sol(i-1,j,k)*f2xmy*(sig(i-1,j-1,k)+sig(i-1,j,k))
                  + sol(i+1,j,k)*f2xmy*(sig(i  ,j-1,k)+sig(i  ,j,k))
                  + sol(i,j-1,k)*fmx2y*(sig(i-1,j-1,k)+sig(i,j-1,k))
                  + sol(i,j+1,k)*fmx2y*(sig(i-1,j  ,k)+sig(i,j  ,k))
                  + sol(i,j,k)*s0;

        if (is_rz) {
            T fp = facy / static_cast<T>(2*i+1);
            T fm = facy / static_cast<T>(2*i-1);
            T frzlo = fm*sig(i-1,j-1,k)-fp*sig(i,j-1,k);
            T frzhi = fm*sig(i-1,j  ,k)-fp*sig(i,j  ,k);
            s0 += - frzhi - frzlo;
            Ax += frzhi*(sol(i,j+1,k)-sol(i,j,k))
                + frzlo*(sol(i,j-1,k)-sol(i,j,k));
        }

        sol(i,j,k) += (rhs(i,j,k) - Ax) / s0;
    }
}

template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void mlndlap_gscolor_c (int i, int j, int k, Array4<T> const& sol,
                        Array4<std::type_identity_t<T> const> const& rhs, std::type_identity_t<T> sig,
                        Array4<int const> const& msk,
                        GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv,
                        bool is_rz) noexcept
{
    if (msk(i,j,k)) {
        sol(i,j,k) = T(0.0);
    } else {
        T facx = T(1.0/6.0)*dxinv[0]*dxinv[0];
        T facy = T(1.0/6.0)*dxinv[1]*dxinv[1];
        T fxy = facx + facy;
        T f2xmy = T(2.0)*facx - facy;
        T fmx2y = T(2.0)*facy - facx;

        T s0 = (-T(2.0))*fxy*T(4.);
        T Ax =   sol(i-1,j-1,k)*fxy
                  + sol(i+1,j-1,k)*fxy
                  + sol(i-1,j+1,k)*fxy
                  + sol(i+1,j+1,k)*fxy
                  + sol(i-1,j,k)*f2xmy*T(2.)
                  + sol(i+1,j,k)*f2xmy*T(2.)
                  + sol(i,j-1,k)*fmx2y*T(2.)
                  + sol(i,j+1,k)*fmx2y*T(2.)
                  + sol(i,j,k)*s0;

        if (is_rz) {
            T fp = facy / static_cast<T>(2*i+1);
            T fm = facy / static_cast<T>(2*i-1);
            T frzlo = fm-fp;
            T frzhi = fm-fp;
            s0 += - frzhi - frzlo;
            Ax += frzhi*(sol(i,j+1,k)-sol(i,j,k))
                + frzlo*(sol(i,j-1,k)-sol(i,j,k));
        }

        sol(i,j,k) += (rhs(i,j,k) - Ax*sig) / (s0*sig);
    }
}

}
#endif
