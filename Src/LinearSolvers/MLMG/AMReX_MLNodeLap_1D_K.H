#ifndef AMREX_MLNODELAP_1D_K_H_
#define AMREX_MLNODELAP_1D_K_H_
#include <AMReX_Config.H>

namespace amrex {

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_zero_fine (int i, int, int, Array4<T> const& phi,
                        Array4<int const> const& msk, int fine_flag) noexcept
{
    // Testing if the node is covered by a fine level in computing
    // coarse sync residual
    if (msk(i-1,0,0) == fine_flag &&
        msk(i  ,0,0) == fine_flag)
    {
        phi(i,0,0) = T(0.0);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_avgdown_coeff_x (int i, int, int, Array4<T> const& crse,
                              Array4<std::type_identity_t<T> const> const& fine) noexcept
{
    T a = fine(2*i  ,0,0);
    T b = fine(2*i+1,0,0);
    crse(i,0,0) = T(2.0)*a*b/(a+b);
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_semi_avgdown_coeff (int i, int j, int k, Array4<T> const& crse,
                                 Array4<std::type_identity_t<T> const> const& fine, int) noexcept
{
    mlndlap_avgdown_coeff_x(i,j,k,crse,fine);
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T mlndlap_adotx_c (int i, int, int, Array4<T const> const& x,
                      std::type_identity_t<T> sigma, Array4<int const> const& msk,
                      GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    if (msk(i,0,0)) {
        return T(0.0);
    } else {
        T y  = x(i-1,0,0)
               +  x(i+1,0,0)
               +  x(i  ,0,0)*T(-2.0);
        return dxinv[0]*dxinv[0] * y * sigma;
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T mlndlap_adotx_ha (int i, int, int, Array4<T const> const& x,
                       Array4<std::type_identity_t<T> const> const& sx, Array4<int const> const& msk,
                       GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    if (msk(i,0,0)) {
        return T(0.0);
    } else {
        T y = (x(i+1,0,0) - x(i  ,0,0)) * sx(i  ,0,0)
            -    (x(i  ,0,0) - x(i-1,0,0)) * sx(i-1,0,0);
        return dxinv[0]*dxinv[0] * y;
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T mlndlap_adotx_aa (int i, int j, int k, Array4<T const> const& x,
                       Array4<std::type_identity_t<T> const> const& sx, Array4<int const> const& msk,
                       GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    return mlndlap_adotx_ha(i,j,k,x,sx,msk,dxinv);
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_normalize_ha (int i, int, int, Array4<T> const& x,
                           Array4<std::type_identity_t<T> const> const& sx, Array4<int const> const& msk,
                           GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    if (!msk(i,0,0)) {
        x(i,0,0) /= -dxinv[0]*dxinv[0] * (sx(i-1,0,0)+sx(i,0,0));
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_normalize_aa (int i, int j, int k, Array4<T> const& x,
                           Array4<std::type_identity_t<T> const> const& sx, Array4<int const> const& msk,
                           GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    mlndlap_normalize_ha(i,j,k,x,sx,msk,dxinv);
}

template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void mlndlap_jacobi_ha (int i, int, int, Array4<T> const& sol, std::type_identity_t<T> Ax,
                        Array4<std::type_identity_t<T> const> const& rhs, Array4<std::type_identity_t<T> const> const& sx,
                        Array4<int const> const& msk,
                        GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    if (msk(i,0,0)) {
        sol(i,0,0) = T(0.0);
    } else {
        sol(i,0,0) += T(2.0/3.0) * (rhs(i,0,0) - Ax)
            / (-dxinv[0]*dxinv[0] * (sx(i-1,0,0)+sx(i,0,0)));
    }
}

template <typename T>
inline
void mlndlap_jacobi_ha (Box const& bx, Array4<T> const& sol, Array4<std::type_identity_t<T> const> const& Ax,
                        Array4<std::type_identity_t<T> const> const& rhs, Array4<std::type_identity_t<T> const> const& sx,
                        Array4<int const> const& msk,
                        GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    T fac = -dxinv[0]*dxinv[0];

    amrex::LoopConcurrentOnCpu(bx, [&] (int i, int, int) noexcept
    {
        if (msk(i,0,0)) {
            sol(i,0,0) = T(0.0);
        } else {
            sol(i,0,0) += T(2.0/3.0) * (rhs(i,0,0) - Ax(i,0,0))
                / (fac * (sx(i-1,0,0)+sx(i,0,0)));
        }
    });
}

template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void mlndlap_jacobi_aa (int i, int j, int k, Array4<T> const& sol, std::type_identity_t<T> Ax,
                        Array4<std::type_identity_t<T> const> const& rhs, Array4<std::type_identity_t<T> const> const& sig,
                        Array4<int const> const& msk, GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    mlndlap_jacobi_ha(i,j,k,sol,Ax,rhs,sig,msk,dxinv);
}

template <typename T>
inline
void mlndlap_jacobi_aa (Box const& bx, Array4<T> const& sol, Array4<std::type_identity_t<T> const> const& Ax,
                        Array4<std::type_identity_t<T> const> const& rhs, Array4<std::type_identity_t<T> const> const& sig,
                        Array4<int const> const& msk, GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    mlndlap_jacobi_ha(bx,sol,Ax,rhs,sig,msk,dxinv);
}

template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void mlndlap_jacobi_c (int i, int, int, Array4<T> const& sol, std::type_identity_t<T> Ax,
                       Array4<std::type_identity_t<T> const> const& rhs, std::type_identity_t<T> sig,
                       Array4<int const> const& msk, GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    if (msk(i,0,0)) {
        sol(i,0,0) = T(0.0);
    } else {
        sol(i,0,0) += T(2.0/3.0) * (rhs(i,0,0) - Ax)
            / (-dxinv[0]*dxinv[0]*T(2.0)*sig);
    }
}

template <typename T>
inline
void mlndlap_jacobi_c (Box const& bx, Array4<T> const& sol, Array4<std::type_identity_t<T> const> const& Ax,
                       Array4<std::type_identity_t<T> const> const& rhs, std::type_identity_t<T> sig,
                       Array4<int const> const& msk, GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    amrex::LoopConcurrentOnCpu(bx, [&] (int i, int, int) noexcept
    {
        if (msk(i,0,0)) {
            sol(i,0,0) = T(0.0);
        } else {
            sol(i,0,0) += T(2.0/3.0) * (rhs(i,0,0) - Ax(i,0,0))
                / (-dxinv[0]*dxinv[0]*T(2.0)*sig);
        }
    });
}

template <typename T>
inline
void mlndlap_gauss_seidel_ha (Box const& bx, Array4<T> const& sol,
                              Array4<std::type_identity_t<T> const> const& rhs,
                              Array4<std::type_identity_t<T> const> const& sx,
                              Array4<int const> const& msk,
                              GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    T fac = dxinv[0]*dxinv[0];

    amrex::LoopOnCpu(bx, [&] (int i, int, int) noexcept
    {
        if (msk(i,0,0)) {
            sol(i,0,0) = T(0.0);
        } else {
            T s0 = T(-1.0) * fac * (sx(i-1,0,0)+sx(i,0,0));
            T Ax = sol(i-1,0,0)*fac*sx(i-1,0,0)
                +     sol(i+1,0,0)*fac*sx(i  ,0,0)
                +     sol(i  ,0,0)*s0;
            sol(i,0,0) += (rhs(i,0,0) - Ax) / s0;
        }
    });
}

template <typename T>
inline
void mlndlap_gauss_seidel_aa (Box const& bx, Array4<T> const& sol,
                              Array4<std::type_identity_t<T> const> const& rhs,
                              Array4<std::type_identity_t<T> const> const& sx,
                              Array4<int const> const& msk,
                              GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    mlndlap_gauss_seidel_ha(bx,sol,rhs,sx,msk,dxinv);
}

template <typename T>
inline
void mlndlap_gauss_seidel_c (Box const& bx, Array4<T> const& sol,
                             Array4<std::type_identity_t<T> const> const& rhs, std::type_identity_t<T> sig,
                             Array4<int const> const& msk,
                             GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    T fac = dxinv[0]*dxinv[0];

    amrex::LoopOnCpu(bx, [&] (int i, int, int) noexcept
    {
        if (msk(i,0,0)) {
            sol(i,0,0) = T(0.0);
        } else {
            T s0 = T(-2.0) * fac * sig;
            T Ax = sol(i-1,0,0)*fac*sig
                +     sol(i+1,0,0)*fac*sig
                +     sol(i  ,0,0)*s0;
            sol(i,0,0) += (rhs(i,0,0) - Ax) / s0;
        }
    });
}

template <typename T>
inline
void mlndlap_gauss_seidel_with_line_solve_aa(Box const&, Array4<T> const&,
                                             Array4<std::type_identity_t<T> const> const&, Array4<std::type_identity_t<T> const> const&,
                                             Array4<int const> const&, GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const&)
{
    amrex::Abort("mlndlap_gauss_seidel_with_line_solve_aa: not implemented in 1D");
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_interpadd_c (int i, int , int, Array4<T> const& fine,
                          Array4<std::type_identity_t<T> const> const& crse,
                          Array4<int const> const& msk) noexcept
{
    if (!msk(i,0,0)) {
        int ic = amrex::coarsen(i,2);
        bool i_is_odd = (ic*2 != i);
        if (i_is_odd) {
            // Node on X line
            fine(i,0,0) += T(0.5) * (crse(ic,0,0)+crse(ic+1,0,0));
        } else {
            // Node coincident with coarse node
            fine(i,0,0) += crse(ic,0,0);
        }
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_interpadd_aa (int i, int , int, Array4<T> const& fine,
                           Array4<std::type_identity_t<T> const> const& crse, Array4<std::type_identity_t<T> const> const& sig,
                           Array4<int const> const& msk) noexcept
{
    if (!msk(i,0,0)) {
        int ic = amrex::coarsen(i,2);
        bool i_is_odd = (ic*2 != i);
        if (i_is_odd) {
            // Node on X line
            fine(i,0,0) += (sig(i-1,0,0)*crse(ic,0,0) + sig(i,0,0)*crse(ic+1,0,0))
                / (sig(i-1,0,0) + sig(i,0,0));
        } else {
            // Node coincident with coarse node
            fine(i,0,0) += crse(ic,0,0);
        }
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_semi_interpadd_aa (int /*i*/, int /*j*/, int /*k*/, Array4<T> const&,
                           Array4<std::type_identity_t<T> const> const&, Array4<std::type_identity_t<T> const> const&,
                           Array4<int const> const&, int) noexcept
{}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_interpadd_ha (int i, int j, int k, Array4<T> const& fine,
                           Array4<std::type_identity_t<T> const> const& crse, Array4<std::type_identity_t<T> const> const& sig,
                           Array4<int const> const& msk) noexcept
{
    mlndlap_interpadd_aa(i,j,k,fine,crse,sig,msk);
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_divu (int i, int, int, Array4<T> const& rhs, Array4<std::type_identity_t<T> const> const& vel,
                   Array4<int const> const& msk, GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv,
                   Box const&,
                   GpuArray<LinOpBCType, AMREX_SPACEDIM> const&,
                   GpuArray<LinOpBCType, AMREX_SPACEDIM> const&) noexcept
{
    T fac = dxinv[0];
    if (msk(i,0,0)) {
        rhs(i,0,0) = T(0.0);
    } else {
        rhs(i,0,0) = fac*(vel(i,0,0)-vel(i-1,0,0));
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T mlndlap_rhcc (int i, int, int, Array4<T const> const& rhcc,
                   Array4<int const> const& msk) noexcept
{
    T r;
    if (msk(i,0,0)) {
        r = T(0.0);
    } else {
        r = T(0.5) * (rhcc(i-1,0,0)+rhcc(i,0,0));
    }
    return r;
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_mknewu (int i, int, int, Array4<T> const& u,
                     Array4<std::type_identity_t<T> const> const& p,
                     Array4<std::type_identity_t<T> const> const& sig,
                     GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    T fac = dxinv[0];
    u(i,0,0) -= sig(i,0,0)*fac*(p(i+1,0,0)-p(i,0,0));
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_mknewu_c (int i, int, int, Array4<T> const& u,
                       Array4<std::type_identity_t<T> const> const& p,
                       std::type_identity_t<T> sig, GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    T fac = dxinv[0];
    u(i,0,0) -= sig*fac*(p(i+1,0,0)-p(i,0,0));
}

template <int rr, typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_rhcc_fine_contrib (int /*i*/, int /*j*/, int /*k*/, Box const&,
                                Array4<T> const&, Array4<std::type_identity_t<T> const> const&,
                                Array4<int const> const&) noexcept
{}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_divu_cf_contrib (int, int, int, Array4<T> const&,
                              Array4<std::type_identity_t<T> const> const&, Array4<std::type_identity_t<T> const> const&,
                              Array4<std::type_identity_t<T> const> const&, Array4<int const> const&,
                              Array4<int const> const&, Array4<int const> const&,
                              GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const&,
                              Box const&, GpuArray<LinOpBCType,AMREX_SPACEDIM> const&,
                              GpuArray<LinOpBCType,AMREX_SPACEDIM> const&,
                              bool) noexcept
{}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_crse_resid (int /*i*/, int /*j*/, int /*k*/, Array4<T> const&,
                         Array4<std::type_identity_t<T> const> const&, Array4<int const> const&,
                         Box const&, GpuArray<LinOpBCType,AMREX_SPACEDIM> const&,
                         GpuArray<LinOpBCType,AMREX_SPACEDIM> const&,
                         bool) noexcept
{}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_res_cf_contrib (int /*i*/, int /*j*/, int /*k*/, Array4<T> const&,
                             Array4<std::type_identity_t<T> const> const&, Array4<std::type_identity_t<T> const> const&,
                             Array4<std::type_identity_t<T> const> const&, Array4<int const> const&,
                             Array4<int const> const&, Array4<int const> const&,
                             Array4<std::type_identity_t<T> const> const&,
                             GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const&, Box const&,
                             GpuArray<LinOpBCType, AMREX_SPACEDIM> const&,
                             GpuArray<LinOpBCType, AMREX_SPACEDIM> const&,
                             bool) noexcept
{}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_res_cf_contrib_cs (int /*i*/, int /*j*/, int /*k*/, Array4<T> const&,
                                Array4<std::type_identity_t<T> const> const&, Array4<std::type_identity_t<T> const> const&,
                                std::type_identity_t<T>, Array4<int const> const&,
                                Array4<int const> const&, Array4<int const> const&,
                                Array4<std::type_identity_t<T> const> const&,
                                GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const&, Box const&,
                                GpuArray<LinOpBCType, AMREX_SPACEDIM> const&,
                                GpuArray<LinOpBCType, AMREX_SPACEDIM> const&,
                                bool) noexcept
{}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_set_stencil (Box const&, Array4<T> const&,
                          Array4<std::type_identity_t<T> const> const&,
                          GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const&) noexcept
{}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_set_stencil_s0 (int /*i*/, int /*j*/, int /*k*/, Array4<T> const&) noexcept
{}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_stencil_rap (int /*i*/, int /*j*/, int /*k*/, Array4<T> const&,
                          Array4<std::type_identity_t<T> const> const&) noexcept
{}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T mlndlap_adotx_sten (int /*i*/, int /*j*/, int /*k*/, Array4<T const> const&,
                         Array4<std::type_identity_t<T> const> const&, Array4<int const> const&) noexcept
{ return T(0.0); }

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_gauss_seidel_sten (int /*i*/, int /*j*/, int /*k*/, Array4<T> const&,
                                Array4<std::type_identity_t<T> const> const&,
                                Array4<std::type_identity_t<T> const> const&,
                                Array4<int const> const&) noexcept
{}

template <typename T>
inline
void mlndlap_gauss_seidel_sten (Box const&, Array4<T> const&,
                                Array4<std::type_identity_t<T> const> const&,
                                Array4<std::type_identity_t<T> const> const&,
                                Array4<int const> const&) noexcept
{}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_interpadd_rap (int /*i*/, int /*j*/, int /*k*/, Array4<T> const&,
                            Array4<std::type_identity_t<T> const> const&, Array4<std::type_identity_t<T> const> const&,
                            Array4<int const> const&) noexcept
{}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_restriction_rap (int /*i*/, int /*j*/, int /*k*/, Array4<T> const&,
                              Array4<std::type_identity_t<T> const> const&, Array4<std::type_identity_t<T> const> const&,
                              Array4<int const> const&) noexcept
{}

template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void mlndlap_gscolor_ha (int i, int, int, Array4<T> const& sol,
                         Array4<std::type_identity_t<T> const> const& rhs,
                         Array4<std::type_identity_t<T> const> const& sx,
                         Array4<int const> const& msk,
                         GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    if (msk(i,0,0)) {
        sol(i,0,0) = T(0.0);
    } else {
        T fac = dxinv[0]*dxinv[0];

        T s0 = T(-1.0) * fac * (sx(i-1,0,0)+sx(i,0,0));
        T Ax = sol(i-1,0,0)*fac*sx(i-1,0,0)
            +     sol(i+1,0,0)*fac*sx(i  ,0,0)
            +     sol(i  ,0,0)*s0;
        sol(i,0,0) += (rhs(i,0,0) - Ax) / s0;
    }
}

template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void mlndlap_gscolor_aa (int i, int j, int k, Array4<T> const& sol,
                         Array4<std::type_identity_t<T> const> const& rhs,
                         Array4<std::type_identity_t<T> const> const& sx,
                         Array4<int const> const& msk,
                         GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    mlndlap_gscolor_ha(i,j,k,sol,rhs,sx,msk,dxinv);
}

template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void mlndlap_gscolor_c (int i, int, int, Array4<T> const& sol,
                        Array4<std::type_identity_t<T> const> const& rhs, std::type_identity_t<T> sig,
                        Array4<int const> const& msk,
                        GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    if (msk(i,0,0)) {
        sol(i,0,0) = T(0.0);
    } else {
        T fac = dxinv[0]*dxinv[0];

        T s0 = T(-2.0) * fac * sig;
        T Ax = sol(i-1,0,0)*fac*sig
            +     sol(i+1,0,0)*fac*sig
            +     sol(i  ,0,0)*s0;
        sol(i,0,0) += (rhs(i,0,0) - Ax) / s0;
    }
}

}

#endif
