#ifndef AMREX_MLABECLAP_3D_K_H_
#define AMREX_MLABECLAP_3D_K_H_
#include <AMReX_Config.H>

namespace amrex {

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlabeclap_adotx (int i, int j, int k, int n, Array4<T> const& y,
                      Array4<T const> const& x,
                      Array4<T const> const& a,
                      Array4<T const> const& bX,
                      Array4<T const> const& bY,
                      Array4<T const> const& bZ,
                      GpuArray<T,AMREX_SPACEDIM> const& dxinv,
                      T alpha, T beta) noexcept
{
    const T dhx = beta*dxinv[0]*dxinv[0];
    const T dhy = beta*dxinv[1]*dxinv[1];
    const T dhz = beta*dxinv[2]*dxinv[2];
    y(i,j,k,n) = alpha*a(i,j,k)*x(i,j,k,n)
        - dhx * (bX(i+1,j,k,n)*(x(i+1,j,k,n) - x(i  ,j,k,n))
               - bX(i  ,j,k,n)*(x(i  ,j,k,n) - x(i-1,j,k,n)))
        - dhy * (bY(i,j+1,k,n)*(x(i,j+1,k,n) - x(i,j  ,k,n))
               - bY(i,j  ,k,n)*(x(i,j  ,k,n) - x(i,j-1,k,n)))
        - dhz * (bZ(i,j,k+1,n)*(x(i,j,k+1,n) - x(i,j,k  ,n))
               - bZ(i,j,k  ,n)*(x(i,j,k  ,n) - x(i,j,k-1,n)));
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlabeclap_adotx_os (int i, int j, int k, int n, Array4<T> const& y,
                         Array4<T const> const& x,
                         Array4<T const> const& a,
                         Array4<T const> const& bX,
                         Array4<T const> const& bY,
                         Array4<T const> const& bZ,
                         Array4<int const> const& osm,
                         GpuArray<T,AMREX_SPACEDIM> const& dxinv,
                         T alpha, T beta) noexcept
{
    const T dhx = beta*dxinv[0]*dxinv[0];
    const T dhy = beta*dxinv[1]*dxinv[1];
    const T dhz = beta*dxinv[2]*dxinv[2];
    if (osm(i,j,k) == 0) {
        y(i,j,k,n) = T(0.0);
    } else {
        y(i,j,k,n) = alpha*a(i,j,k)*x(i,j,k,n)
            - dhx * (bX(i+1,j,k,n)*(x(i+1,j,k,n) - x(i  ,j,k,n))
                   - bX(i  ,j,k,n)*(x(i  ,j,k,n) - x(i-1,j,k,n)))
            - dhy * (bY(i,j+1,k,n)*(x(i,j+1,k,n) - x(i,j  ,k,n))
                   - bY(i,j  ,k,n)*(x(i,j  ,k,n) - x(i,j-1,k,n)))
            - dhz * (bZ(i,j,k+1,n)*(x(i,j,k+1,n) - x(i,j,k  ,n))
                   - bZ(i,j,k  ,n)*(x(i,j,k  ,n) - x(i,j,k-1,n)));
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlabeclap_normalize (int i, int j, int k, int n, Array4<T> const& x,
                          Array4<T const> const& a,
                          Array4<T const> const& bX,
                          Array4<T const> const& bY,
                          Array4<T const> const& bZ,
                          GpuArray<T,AMREX_SPACEDIM> const& dxinv,
                          T alpha, T beta) noexcept
{
    const T dhx = beta*dxinv[0]*dxinv[0];
    const T dhy = beta*dxinv[1]*dxinv[1];
    const T dhz = beta*dxinv[2]*dxinv[2];
    x(i,j,k,n) /= alpha*a(i,j,k)
        + dhx*(bX(i,j,k,n)+bX(i+1,j,k,n))
        + dhy*(bY(i,j,k,n)+bY(i,j+1,k,n))
        + dhz*(bZ(i,j,k,n)+bZ(i,j,k+1,n));
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlabeclap_flux_x (Box const& box, Array4<T> const& fx, Array4<T const> const& sol,
                       Array4<T const> const& bx, T fac, int ncomp) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int n = 0; n < ncomp; ++n) {
    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                fx(i,j,k,n) = -fac*bx(i,j,k,n)*(sol(i,j,k,n)-sol(i-1,j,k,n));
            }
        }
    }
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlabeclap_flux_xface (Box const& box, Array4<T> const& fx, Array4<T const> const& sol,
                           Array4<T const> const& bx, T fac, int xlen, int ncomp) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int n = 0; n < ncomp; ++n) {
    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            int i = lo.x;
            fx(i,j,k,n) = -fac*bx(i,j,k,n)*(sol(i,j,k,n)-sol(i-1,j,k,n));
            i += xlen;
            fx(i,j,k,n) = -fac*bx(i,j,k,n)*(sol(i,j,k,n)-sol(i-1,j,k,n));
        }
    }
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlabeclap_flux_y (Box const& box, Array4<T> const& fy, Array4<T const> const& sol,
                       Array4<T const> const& by, T fac, int ncomp) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int n = 0; n < ncomp; ++n) {
    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                fy(i,j,k,n) = -fac*by(i,j,k,n)*(sol(i,j,k,n)-sol(i,j-1,k,n));
            }
        }
    }
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlabeclap_flux_yface (Box const& box, Array4<T> const& fy, Array4<T const> const& sol,
                           Array4<T const> const& by, T fac, int ylen, int ncomp) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int n = 0; n < ncomp; ++n) {
    for     (int k = lo.z; k <= hi.z; ++k) {
        int j = lo.y;
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            fy(i,j,k,n) = -fac*by(i,j,k,n)*(sol(i,j,k,n)-sol(i,j-1,k,n));
        }
        j += ylen;
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            fy(i,j,k,n) = -fac*by(i,j,k,n)*(sol(i,j,k,n)-sol(i,j-1,k,n));
        }
    }
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlabeclap_flux_z (Box const& box, Array4<T> const& fz, Array4<T const> const& sol,
                       Array4<T const> const& bz, T fac, int ncomp) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int n = 0; n < ncomp; ++n) {
    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                fz(i,j,k,n) = -fac*bz(i,j,k,n)*(sol(i,j,k,n)-sol(i,j,k-1,n));
            }
        }
    }
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlabeclap_flux_zface (Box const& box, Array4<T> const& fz, Array4<T const> const& sol,
                           Array4<T const> const& bz, T fac, int zlen, int ncomp) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int n = 0; n < ncomp; ++n) {
    int k = lo.z;
    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            fz(i,j,k,n) = -fac*bz(i,j,k,n)*(sol(i,j,k,n)-sol(i,j,k-1,n));
        }
    }

    k += zlen;
    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            fz(i,j,k,n) = -fac*bz(i,j,k,n)*(sol(i,j,k,n)-sol(i,j,k-1,n));
        }
    }
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void abec_gsrb (int i, int j, int k, int n, Array4<T> const& phi, Array4<T const> const& rhs,
                T alpha, Array4<T const> const& a,
                T dhx, T dhy, T dhz,
                Array4<T const> const& bX, Array4<T const> const& bY,
                Array4<T const> const& bZ,
                Array4<int const> const& m0, Array4<int const> const& m2,
                Array4<int const> const& m4,
                Array4<int const> const& m1, Array4<int const> const& m3,
                Array4<int const> const& m5,
                Array4<T const> const& f0, Array4<T const> const& f2,
                Array4<T const> const& f4,
                Array4<T const> const& f1, Array4<T const> const& f3,
                Array4<T const> const& f5,
                Box const& vbox, int redblack) noexcept
{
    constexpr T omega = T(1.15);

    if ((i+j+k+redblack)%2 == 0) {
        const auto vlo = amrex::lbound(vbox);
        const auto vhi = amrex::ubound(vbox);

        T cf0 = (i == vlo.x && m0(vlo.x-1,j,k) > 0)
            ? f0(vlo.x,j,k,n) : T(0.0);
        T cf1 = (j == vlo.y && m1(i,vlo.y-1,k) > 0)
            ? f1(i,vlo.y,k,n) : T(0.0);
        T cf2 = (k == vlo.z && m2(i,j,vlo.z-1) > 0)
            ? f2(i,j,vlo.z,n) : T(0.0);
        T cf3 = (i == vhi.x && m3(vhi.x+1,j,k) > 0)
            ? f3(vhi.x,j,k,n) : T(0.0);
        T cf4 = (j == vhi.y && m4(i,vhi.y+1,k) > 0)
            ? f4(i,vhi.y,k,n) : T(0.0);
        T cf5 = (k == vhi.z && m5(i,j,vhi.z+1) > 0)
            ? f5(i,j,vhi.z,n) : T(0.0);

        T gamma = alpha*a(i,j,k)
            +   dhx*(bX(i,j,k,n)+bX(i+1,j,k,n))
            +   dhy*(bY(i,j,k,n)+bY(i,j+1,k,n))
            +   dhz*(bZ(i,j,k,n)+bZ(i,j,k+1,n));

        T g_m_d = gamma
            - (dhx*(bX(i,j,k,n)*cf0 + bX(i+1,j,k,n)*cf3)
            +  dhy*(bY(i,j,k,n)*cf1 + bY(i,j+1,k,n)*cf4)
            +  dhz*(bZ(i,j,k,n)*cf2 + bZ(i,j,k+1,n)*cf5));

        T rho =  dhx*( bX(i  ,j,k,n)*phi(i-1,j,k,n)
               +       bX(i+1,j,k,n)*phi(i+1,j,k,n) )
               + dhy*( bY(i,j  ,k,n)*phi(i,j-1,k,n)
               +       bY(i,j+1,k,n)*phi(i,j+1,k,n) )
               + dhz*( bZ(i,j,k  ,n)*phi(i,j,k-1,n)
               +       bZ(i,j,k+1,n)*phi(i,j,k+1,n) );

        T res =  rhs(i,j,k,n) - (gamma*phi(i,j,k,n) - rho);
        phi(i,j,k,n) = phi(i,j,k,n) + omega/g_m_d * res;
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void abec_gsrb_os (int i, int j, int k, int n,
                   Array4<T> const& phi, Array4<T const> const& rhs,
                   T alpha, Array4<T const> const& a,
                   T dhx, T dhy, T dhz,
                   Array4<T const> const& bX, Array4<T const> const& bY,
                   Array4<T const> const& bZ,
                   Array4<int const> const& m0, Array4<int const> const& m2,
                   Array4<int const> const& m4,
                   Array4<int const> const& m1, Array4<int const> const& m3,
                   Array4<int const> const& m5,
                   Array4<T const> const& f0, Array4<T const> const& f2,
                   Array4<T const> const& f4,
                   Array4<T const> const& f1, Array4<T const> const& f3,
                   Array4<T const> const& f5,
                   Array4<int const> const& osm,
                   Box const& vbox, int redblack) noexcept
{
    constexpr T omega = T(1.15);

    if ((i+j+k+redblack)%2 == 0) {
        if (osm(i,j,k) == 0) {
            phi(i,j,k,n) = T(0.0);
        } else {
            const auto vlo = amrex::lbound(vbox);
            const auto vhi = amrex::ubound(vbox);

            T cf0 = (i == vlo.x && m0(vlo.x-1,j,k) > 0)
                ? f0(vlo.x,j,k,n) : T(0.0);
            T cf1 = (j == vlo.y && m1(i,vlo.y-1,k) > 0)
                ? f1(i,vlo.y,k,n) : T(0.0);
            T cf2 = (k == vlo.z && m2(i,j,vlo.z-1) > 0)
                ? f2(i,j,vlo.z,n) : T(0.0);
            T cf3 = (i == vhi.x && m3(vhi.x+1,j,k) > 0)
                ? f3(vhi.x,j,k,n) : T(0.0);
            T cf4 = (j == vhi.y && m4(i,vhi.y+1,k) > 0)
                ? f4(i,vhi.y,k,n) : T(0.0);
            T cf5 = (k == vhi.z && m5(i,j,vhi.z+1) > 0)
                ? f5(i,j,vhi.z,n) : T(0.0);

            T gamma = alpha*a(i,j,k)
                +   dhx*(bX(i,j,k,n)+bX(i+1,j,k,n))
                +   dhy*(bY(i,j,k,n)+bY(i,j+1,k,n))
                +   dhz*(bZ(i,j,k,n)+bZ(i,j,k+1,n));

            T g_m_d = gamma
                - (dhx*(bX(i,j,k,n)*cf0 + bX(i+1,j,k,n)*cf3)
                +  dhy*(bY(i,j,k,n)*cf1 + bY(i,j+1,k,n)*cf4)
                +  dhz*(bZ(i,j,k,n)*cf2 + bZ(i,j,k+1,n)*cf5));

            T rho =  dhx*( bX(i  ,j,k,n)*phi(i-1,j,k,n)
                   +       bX(i+1,j,k,n)*phi(i+1,j,k,n) )
                   + dhy*( bY(i,j  ,k,n)*phi(i,j-1,k,n)
                   +       bY(i,j+1,k,n)*phi(i,j+1,k,n) )
                   + dhz*( bZ(i,j,k  ,n)*phi(i,j,k-1,n)
                   +       bZ(i,j,k+1,n)*phi(i,j,k+1,n) );

            T res =  rhs(i,j,k,n) - (gamma*phi(i,j,k,n) - rho);
            phi(i,j,k,n) = phi(i,j,k,n) + omega/g_m_d * res;
        }
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void tridiagonal_solve (Array1D<T,0,31>& a_ls, Array1D<T,0,31>& b_ls, Array1D<T,0,31>& c_ls,
                        Array1D<T,0,31>& r_ls, Array1D<T,0,31>& u_ls, Array1D<T,0,31>& gam,
                        int ilen ) noexcept
{
    T bet = b_ls(0);
    u_ls(0) = r_ls(0) / bet;

    for (int i = 1; i <= ilen - 1; i++) {
        gam(i) = c_ls(i-1) / bet;
        bet = b_ls(i) - a_ls(i)*gam(i);
        if (bet == 0) amrex::Abort(">>>TRIDIAG FAILED");
        u_ls(i) = (r_ls(i)-a_ls(i)*u_ls(i-1)) / bet;
    }
    for (int i = ilen-2; i >= 0; i--) {
        u_ls(i) = u_ls(i) - gam(i+1)*u_ls(i+1);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void abec_gsrb_with_line_solve (
                Box const& box, Array4<T> const& phi, Array4<T const> const& rhs,
                T alpha, Array4<T const> const& a,
                T dhx, T dhy, T dhz,
                Array4<T const> const& bX, Array4<T const> const& bY,
                Array4<T const> const& bZ,
                Array4<int const> const& m0, Array4<int const> const& m2,
                Array4<int const> const& m4,
                Array4<int const> const& m1, Array4<int const> const& m3,
                Array4<int const> const& m5,
                Array4<T const> const& f0, Array4<T const> const& f2,
                Array4<T const> const& f4,
                Array4<T const> const& f1, Array4<T const> const& f3,
                Array4<T const> const& f5,
                Box const& vbox, int redblack, int nc) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    const auto vlo = amrex::lbound(vbox);
    const auto vhi = amrex::ubound(vbox);

    // amrex::Print() << "GSRB LS " << box << " " << dhx << " " << dhy << " " << dhz << std::endl;

    // idir is the direction in which we will do the tridiagonal solve --
    //    it should be the direction in which the mesh spacing is much larger
    //    than in the other directions
    int idir = 2;
    int ilen = hi.z - lo.z + 1;

    if ( (dhx <= dhy) && (dhz <= dhy) ) {
        idir = 1;
        ilen = hi.y - lo.y + 1;
    }
    if ( (dhy <= dhx) && (dhz <= dhx) ) {
        idir = 0;
        ilen = hi.x - lo.x + 1;
    }

    // This assertion should be moved outside the kernel for performance!
    if (ilen > 32) amrex::Abort("abec_gsrb_with_line_solve is hard-wired to be no longer than 32");

    Array1D<T,0,31> a_ls;
    Array1D<T,0,31> b_ls;
    Array1D<T,0,31> c_ls;
    Array1D<T,0,31> r_ls;
    Array1D<T,0,31> u_ls;
    Array1D<T,0,31> gam;

    if (idir == 2) {
        for (int n = 0; n < nc; ++n) {
            for (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    if ((i+j+redblack)%2 == 0) {

                        for (int k = lo.z; k <= hi.z; ++k)
                        {
                            T gamma = alpha*a(i,j,k)
                                +   dhx*(bX(i,j,k,n)+bX(i+1,j,k,n))
                                +   dhy*(bY(i,j,k,n)+bY(i,j+1,k,n))
                                +   dhz*(bZ(i,j,k,n)+bZ(i,j,k+1,n));

                            T cf0 = (i == vlo.x && m0(vlo.x-1,j,k) > 0)
                                ? f0(vlo.x,j,k,n) : T(0.0);
                            T cf1 = (j == vlo.y && m1(i,vlo.y-1,k) > 0)
                                ? f1(i,vlo.y,k,n) : T(0.0);
                            T cf2 = (k == vlo.z && m2(i,j,vlo.z-1) > 0)
                                ? f2(i,j,vlo.z,n) : T(0.0);
                            T cf3 = (i == vhi.x && m3(vhi.x+1,j,k) > 0)
                                ? f3(vhi.x,j,k,n) : T(0.0);
                            T cf4 = (j == vhi.y && m4(i,vhi.y+1,k) > 0)
                                ? f4(i,vhi.y,k,n) : T(0.0);
                            T cf5 = (k == vhi.z && m5(i,j,vhi.z+1) > 0)
                                ? f5(i,j,vhi.z,n) : T(0.0);

                            T g_m_d = gamma
                                - (dhx*(bX(i,j,k,n)*cf0 + bX(i+1,j,k,n)*cf3)
                                +  dhy*(bY(i,j,k,n)*cf1 + bY(i,j+1,k,n)*cf4)
                                +  dhz*(bZ(i,j,k,n)*cf2 + bZ(i,j,k+1,n)*cf5));

                            T rho =  dhx*( bX(i  ,j,k,n)*phi(i-1,j,k,n)
                                   +       bX(i+1,j,k,n)*phi(i+1,j,k,n) )
                                   + dhy*( bY(i,j  ,k,n)*phi(i,j-1,k,n)
                                   +       bY(i,j+1,k,n)*phi(i,j+1,k,n) );

                            // We have already accounted for this external boundary in the coefficient of phi(i,j,k,n)
                            if (i == vlo.x && m0(vlo.x-1,j,k) > 0)
                                rho -= dhx*bX(i  ,j,k,n)*phi(i-1,j,k,n);
                            if (i == vhi.x && m3(vhi.x+1,j,k) > 0)
                                rho -= dhx*bX(i+1,j,k,n)*phi(i+1,j,k,n);
                            if (j == vlo.y && m1(i,vlo.y-1,k) > 0)
                                rho -= dhy*bY(i,j  ,k,n)*phi(i,j-1,k,n);
                            if (j == vhi.y && m4(i,vhi.y+1,k) > 0)
                                rho -= dhy*bY(i,j+1,k,n)*phi(i,j+1,k,n);

                            a_ls(k-lo.z) = -dhz*bZ(i,j,k,n);
                            b_ls(k-lo.z) =  g_m_d;
                            c_ls(k-lo.z) = -dhz*bZ(i,j,k+1,n);
                            u_ls(k-lo.z) = T(0.);
                            r_ls(k-lo.z) = rhs(i,j,k,n) + rho;
                            // r_ls(k-lo.z) = g_m_d*phi(i,j,k,n) -gamma*phi(i,j,k,n) + rhs(i,j,k,n) + rho;

                            if (k == lo.z)
                            {
                                a_ls(k-lo.z) = T(0.);
                                if (!(m2(i,j,vlo.z-1) > 0)) r_ls(k-lo.z) += dhz*bZ(i,j,k,n)*phi(i,j,k-1,n);
                            }
                            if (k == hi.z)
                            {
                                c_ls(k-lo.z) = T(0.);
                                if (!(m5(i,j,vhi.z+1) > 0)) r_ls(k-lo.z) += dhz*bZ(i,j,k+1,n)*phi(i,j,k+1,n);
                            }
                        }

                        tridiagonal_solve(a_ls, b_ls, c_ls, r_ls, u_ls, gam, ilen);

                        for (int k = lo.z; k <= hi.z; ++k)
                        {
                            phi(i,j,k,n) = u_ls(k-lo.z);
                        }
                    }
                }
            }
        }
    } else if (idir == 1) {
        for (int n = 0; n < nc; ++n) {
            for (int i = lo.x; i <= hi.x; ++i) {
                AMREX_PRAGMA_SIMD
                for (int k = lo.z; k <= hi.z; ++k) {
                    if ((i+k+redblack)%2 == 0) {

                        for (int j = lo.y; j <= hi.y; ++j)
                        {
                            T gamma = alpha*a(i,j,k)
                                +   dhx*(bX(i,j,k,n)+bX(i+1,j,k,n))
                                +   dhy*(bY(i,j,k,n)+bY(i,j+1,k,n))
                                +   dhz*(bZ(i,j,k,n)+bZ(i,j,k+1,n));

                            T cf0 = (i == vlo.x && m0(vlo.x-1,j,k) > 0)
                                ? f0(vlo.x,j,k,n) : T(0.0);
                            T cf1 = (j == vlo.y && m1(i,vlo.y-1,k) > 0)
                                ? f1(i,vlo.y,k,n) : T(0.0);
                            T cf2 = (k == vlo.z && m2(i,j,vlo.z-1) > 0)
                                ? f2(i,j,vlo.z,n) : T(0.0);
                            T cf3 = (i == vhi.x && m3(vhi.x+1,j,k) > 0)
                                ? f3(vhi.x,j,k,n) : T(0.0);
                            T cf4 = (j == vhi.y && m4(i,vhi.y+1,k) > 0)
                                ? f4(i,vhi.y,k,n) : T(0.0);
                            T cf5 = (k == vhi.z && m5(i,j,vhi.z+1) > 0)
                                ? f5(i,j,vhi.z,n) : T(0.0);

                            T g_m_d = gamma
                                - (dhx*(bX(i,j,k,n)*cf0 + bX(i+1,j,k,n)*cf3)
                                +  dhy*(bY(i,j,k,n)*cf1 + bY(i,j+1,k,n)*cf4)
                                +  dhz*(bZ(i,j,k,n)*cf2 + bZ(i,j,k+1,n)*cf5));

                            T rho =  dhx*( bX(i  ,j,k,n)*phi(i-1,j,k,n)
                                   +       bX(i+1,j,k,n)*phi(i+1,j,k,n) )
                                   + dhz*( bZ(i,j  ,k,n)*phi(i,j,k-1,n)
                                   +       bZ(i,j,k+1,n)*phi(i,j,k+1,n) );

                            // We have already accounted for this external boundary in the coefficient of phi(i,j,k,n)
                            if (i == vlo.x && m0(vlo.x-1,j,k) > 0)
                                rho -= dhx*bX(i  ,j,k,n)*phi(i-1,j,k,n);
                            if (i == vhi.x && m3(vhi.x+1,j,k) > 0)
                                rho -= dhx*bX(i+1,j,k,n)*phi(i+1,j,k,n);
                            if (k == vlo.z && m2(i,j,vlo.z-1) > 0)
                                rho -= dhz*bZ(i,j  ,k,n)*phi(i,j,k-1,n);
                            if (k == vhi.z && m5(i,j,vhi.z+1) > 0)
                                rho -= dhz*bZ(i,j,k+1,n)*phi(i,j,k+1,n);

                            a_ls(j-lo.y) = -dhy*bY(i,j,k,n);
                            b_ls(j-lo.y) =  g_m_d;
                            c_ls(j-lo.y) = -dhy*bY(i,j+1,k,n);
                            u_ls(j-lo.y) = T(0.);
                            r_ls(j-lo.y) = rhs(i,j,k,n) + rho;

                            if (j == lo.y)
                            {
                                a_ls(j-lo.y) = T(0.);
                                if (!(m1(i,vlo.y-1,k) > 0)) r_ls(j-lo.y) += dhy*bY(i,j,k,n)*phi(i,j-1,k,n);
                            }
                            if (j == hi.y)
                            {
                                c_ls(j-lo.y) = T(0.);
                                if (!(m4(i,vhi.y+1,k) > 0)) r_ls(j-lo.y) += dhy*bY(i,j+1,k,n)*phi(i,j+1,k,n);
                            }
                        }

                        tridiagonal_solve(a_ls, b_ls, c_ls, r_ls, u_ls, gam, ilen);

                        for (int j = lo.y; j <= hi.y; ++j)
                        {
                            phi(i,j,k,n) = u_ls(j-lo.y);
                        }
                    }
                }
            }
        }
    } else if (idir == 0) {
        for (int n = 0; n < nc; ++n) {
            for (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int k = lo.z; k <= hi.z; ++k) {
                    if ((j+k+redblack)%2 == 0) {

                        for (int i = lo.x; i <= hi.x; ++i)
                        {
                            T gamma = alpha*a(i,j,k)
                                +   dhx*(bX(i,j,k,n)+bX(i+1,j,k,n))
                                +   dhy*(bY(i,j,k,n)+bY(i,j+1,k,n))
                                +   dhz*(bZ(i,j,k,n)+bZ(i,j,k+1,n));

                            T cf0 = (i == vlo.x && m0(vlo.x-1,j,k) > 0)
                                ? f0(vlo.x,j,k,n) : T(0.0);
                            T cf1 = (j == vlo.y && m1(i,vlo.y-1,k) > 0)
                                ? f1(i,vlo.y,k,n) : T(0.0);
                            T cf2 = (k == vlo.z && m2(i,j,vlo.z-1) > 0)
                                ? f2(i,j,vlo.z,n) : T(0.0);
                            T cf3 = (i == vhi.x && m3(vhi.x+1,j,k) > 0)
                                ? f3(vhi.x,j,k,n) : T(0.0);
                            T cf4 = (j == vhi.y && m4(i,vhi.y+1,k) > 0)
                                ? f4(i,vhi.y,k,n) : T(0.0);
                            T cf5 = (k == vhi.z && m5(i,j,vhi.z+1) > 0)
                                ? f5(i,j,vhi.z,n) : T(0.0);

                            T g_m_d = gamma
                                - (dhx*(bX(i,j,k,n)*cf0 + bX(i+1,j,k,n)*cf3)
                                +  dhy*(bY(i,j,k,n)*cf1 + bY(i,j+1,k,n)*cf4)
                                +  dhz*(bZ(i,j,k,n)*cf2 + bZ(i,j,k+1,n)*cf5));

                            T rho =  dhy*( bY(i,j  ,k,n)*phi(i,j-1,k,n)
                                   +       bY(i,j+1,k,n)*phi(i,j+1,k,n) )
                                   + dhz*( bZ(i,j  ,k,n)*phi(i,j,k-1,n)
                                   +       bZ(i,j,k+1,n)*phi(i,j,k+1,n) );

                            // We have already accounted for this external boundary in the coefficient of phi(i,j,k,n)
                            if (j == vlo.y && m1(i,vlo.y-1,k) > 0)
                                rho -= dhy*bY(i,j  ,k,n)*phi(i,j-1,k,n);
                            if (j == vhi.y && m4(i,vhi.y+1,k) > 0)
                                rho -= dhy*bY(i,j+1,k,n)*phi(i,j+1,k,n);
                            if (k == vlo.z && m2(i,j,vlo.z-1) > 0)
                                rho -= dhz*bZ(i,j  ,k,n)*phi(i,j,k-1,n);
                            if (k == vhi.z && m5(i,j,vhi.z+1) > 0)
                                rho -= dhz*bZ(i,j,k+1,n)*phi(i,j,k+1,n);

                            a_ls(i-lo.x) = -dhx*bX(i,j,k,n);
                            b_ls(i-lo.x) =  g_m_d;
                            c_ls(i-lo.x) = -dhx*bX(i+1,j,k,n);
                            u_ls(i-lo.x) = T(0.);
                            r_ls(i-lo.x) = rhs(i,j,k,n) + rho;

                            if (i == lo.x)
                            {
                                a_ls(i-lo.x) = T(0.);
                                if (!(m0(vlo.x-1,j,k) > 0)) r_ls(i-lo.x) += dhx*bX(i,j,k,n)*phi(i-1,j,k,n);
                            }
                            if (i == hi.x)
                            {
                                c_ls(i-lo.x) = T(0.);
                                if (!(m3(vhi.x+1,j,k) > 0)) r_ls(i-lo.x) += dhx*bX(i+1,j,k,n)*phi(i+1,j,k,n);
                            }
                        }

                        tridiagonal_solve(a_ls, b_ls, c_ls, r_ls, u_ls, gam, ilen);

                        for (int i = lo.x; i <= hi.x; ++i)
                        {
                            phi(i,j,k,n) = u_ls(i-lo.x);
                        }
                    }
                }
            }
        }
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void overset_rescale_bcoef_x (Box const& box, Array4<T> const& bX, Array4<int const> const& osm,
                              int ncomp, T osfac) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    for (int n = 0; n < ncomp; ++n) {
        for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
            if ((osm(i-1,j,k)+osm(i,j,k)) == 1) {
                bX(i,j,k,n) *= osfac;
            }
        }}}
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void overset_rescale_bcoef_y (Box const& box, Array4<T> const& bY, Array4<int const> const& osm,
                              int ncomp, T osfac) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    for (int n = 0; n < ncomp; ++n) {
        for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
            if ((osm(i,j-1,k)+osm(i,j,k)) == 1) {
                bY(i,j,k,n) *= osfac;
            }
        }}}
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void overset_rescale_bcoef_z (Box const& box, Array4<T> const& bZ, Array4<int const> const& osm,
                              int ncomp, T osfac) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    for (int n = 0; n < ncomp; ++n) {
        for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
            if ((osm(i,j,k-1)+osm(i,j,k)) == 1) {
                bZ(i,j,k,n) *= osfac;
            }
        }}}
    }
}

}
#endif
