#ifndef AMREX_ML_EB_TENSOR_3D_K_H_
#define AMREX_ML_EB_TENSOR_3D_K_H_
#include <AMReX_Config.H>

#include <AMReX_MLEBABecLap_K.H>

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mlebtensor_dz_on_xface (int i, int j, int, int n,
                             Array4<Real const> const& vel, Real dzi,
                             Real whi, Real wlo,
                             int khip, int khim, int klop, int klom) noexcept
{
    return Real(0.5)*dzi * ((vel(i  ,j,khip,n)-vel(i  ,j,khim,n))*whi +
                            (vel(i-1,j,klop,n)-vel(i-1,j,klom,n))*wlo);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mlebtensor_dz_on_yface (int i, int j, int, int n,
                             Array4<Real const> const& vel, Real dzi,
                             Real whi, Real wlo,
                             int khip, int khim, int klop, int klom) noexcept
{
    return Real(0.5)*dzi * ((vel(i,j  ,khip,n)-vel(i,j  ,khim,n))*whi +
                            (vel(i,j-1,klop,n)-vel(i,j-1,klom,n))*wlo);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mlebtensor_dx_on_zface (int, int j, int k, int n,
                             Array4<Real const> const& vel, Real dxi,
                             Real whi, Real wlo,
                             int ihip, int ihim, int ilop, int ilom) noexcept
{
    return Real(0.5)*dxi * ((vel(ihip,j,k  ,n)-vel(ihim,j,k  ,n))*whi +
                            (vel(ilop,j,k-1,n)-vel(ilom,j,k-1,n))*wlo);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mlebtensor_dy_on_zface (int i, int, int k, int n,
                             Array4<Real const> const& vel, Real dyi,
                             Real whi, Real wlo,
                             int jhip, int jhim, int jlop, int jlom) noexcept
{
    return Real(0.5)*dyi * ((vel(i,jhip,k  ,n)-vel(i,jhim,k  ,n))*whi +
                            (vel(i,jlop,k-1,n)-vel(i,jlom,k-1,n))*wlo);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebtensor_cross_terms_fx (Box const& box, Array4<Real> const& fx,
                                Array4<Real const> const& vel,
                                Array4<Real const> const& etax,
                                Array4<Real const> const& kapx,
                                Array4<Real const> const& apx,
                                Array4<EBCellFlag const> const& flag,
                                GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    const Real dyi = dxinv[1];
    const Real dzi = dxinv[2];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    constexpr Real twoThirds = 2./3.;

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                if (apx(i,j,k) == Real(0.0))
                {
                    fx(i,j,k,0) = Real(0.0);
                    fx(i,j,k,1) = Real(0.0);
                    fx(i,j,k,2) = Real(0.0);
                }
                else
                {
                    int jhip = j + flag(i  ,j,k).isConnected(0, 1,0);
                    int jhim = j - flag(i  ,j,k).isConnected(0,-1,0);
                    int jlop = j + flag(i-1,j,k).isConnected(0, 1,0);
                    int jlom = j - flag(i-1,j,k).isConnected(0,-1,0);
                    Real whi = mlebtensor_weight(jhip-jhim);
                    Real wlo = mlebtensor_weight(jlop-jlom);
                    Real dudy = mlebtensor_dy_on_xface(i,j,k,0,vel,dyi,
                                                       whi,wlo,jhip,jhim,jlop,jlom);
                    Real dvdy = mlebtensor_dy_on_xface(i,j,k,1,vel,dyi,
                                                       whi,wlo,jhip,jhim,jlop,jlom);
                    int khip = k + flag(i  ,j,k).isConnected(0,0, 1);
                    int khim = k - flag(i  ,j,k).isConnected(0,0,-1);
                    int klop = k + flag(i-1,j,k).isConnected(0,0, 1);
                    int klom = k - flag(i-1,j,k).isConnected(0,0,-1);
                    whi = mlebtensor_weight(khip-khim);
                    wlo = mlebtensor_weight(klop-klom);
                    Real dudz = mlebtensor_dz_on_xface(i,j,k,0,vel,dzi,
                                                       whi,wlo,khip,khim,klop,klom);
                    Real dwdz = mlebtensor_dz_on_xface(i,j,k,2,vel,dzi,
                                                       whi,wlo,khip,khim,klop,klom);
                    Real divu = dvdy + dwdz;
                    Real xif = kapx(i,j,k);
                    Real mun = Real(0.75)*(etax(i,j,k,0)-xif);// restore the original eta
                    Real mut =             etax(i,j,k,1);
                    fx(i,j,k,0) = -mun*(-twoThirds*divu) - xif*divu;
                    fx(i,j,k,1) = -mut*dudy;
                    fx(i,j,k,2) = -mut*dudz;
                }
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebtensor_cross_terms_fy (Box const& box, Array4<Real> const& fy,
                                Array4<Real const> const& vel,
                                Array4<Real const> const& etay,
                                Array4<Real const> const& kapy,
                                Array4<Real const> const& apy,
                                Array4<EBCellFlag const> const& flag,
                                GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    const Real dxi = dxinv[0];
    const Real dzi = dxinv[2];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    constexpr Real twoThirds = 2./3.;

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                if (apy(i,j,k) == Real(0.0))
                {
                    fy(i,j,k,0) = Real(0.0);
                    fy(i,j,k,1) = Real(0.0);
                    fy(i,j,k,2) = Real(0.0);
                }
                else
                {
                    int ihip = i + flag(i,j  ,k).isConnected( 1,0,0);
                    int ihim = i - flag(i,j  ,k).isConnected(-1,0,0);
                    int ilop = i + flag(i,j-1,k).isConnected( 1,0,0);
                    int ilom = i - flag(i,j-1,k).isConnected(-1,0,0);
                    Real whi = mlebtensor_weight(ihip-ihim);
                    Real wlo = mlebtensor_weight(ilop-ilom);
                    Real dudx = mlebtensor_dx_on_yface(i,j,k,0,vel,dxi,
                                                       whi,wlo,ihip,ihim,ilop,ilom);
                    Real dvdx = mlebtensor_dx_on_yface(i,j,k,1,vel,dxi,
                                                       whi,wlo,ihip,ihim,ilop,ilom);
                    int khip = k + flag(i,j  ,k).isConnected(0,0, 1);
                    int khim = k - flag(i,j  ,k).isConnected(0,0,-1);
                    int klop = k + flag(i,j-1,k).isConnected(0,0, 1);
                    int klom = k - flag(i,j-1,k).isConnected(0,0,-1);
                    whi = mlebtensor_weight(khip-khim);
                    wlo = mlebtensor_weight(klop-klom);
                    Real dvdz = mlebtensor_dz_on_yface(i,j,k,1,vel,dzi,
                                                       whi,wlo,khip,khim,klop,klom);
                    Real dwdz = mlebtensor_dz_on_yface(i,j,k,2,vel,dzi,
                                                       whi,wlo,khip,khim,klop,klom);
                    Real divu = dudx + dwdz;
                    Real xif = kapy(i,j,k);
                    Real mun = Real(0.75)*(etay(i,j,k,1)-xif);// restore the original eta
                    Real mut =             etay(i,j,k,0);
                    fy(i,j,k,0) = -mut*dvdx;
                    fy(i,j,k,1) = -mun*(-twoThirds*divu) - xif*divu;
                    fy(i,j,k,2) = -mut*dvdz;
                }
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebtensor_cross_terms_fz (Box const& box, Array4<Real> const& fz,
                                Array4<Real const> const& vel,
                                Array4<Real const> const& etaz,
                                Array4<Real const> const& kapz,
                                Array4<Real const> const& apz,
                                Array4<EBCellFlag const> const& flag,
                                GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    const Real dxi = dxinv[0];
    const Real dyi = dxinv[1];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    constexpr Real twoThirds = 2./3.;

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                if (apz(i,j,k) == Real(0.0))
                {
                    fz(i,j,k,0) = Real(0.0);
                    fz(i,j,k,1) = Real(0.0);
                    fz(i,j,k,2) = Real(0.0);
                }
                else
                {
                    int ihip = i + flag(i,j,k  ).isConnected( 1,0,0);
                    int ihim = i - flag(i,j,k  ).isConnected(-1,0,0);
                    int ilop = i + flag(i,j,k-1).isConnected( 1,0,0);
                    int ilom = i - flag(i,j,k-1).isConnected(-1,0,0);
                    Real whi = mlebtensor_weight(ihip-ihim);
                    Real wlo = mlebtensor_weight(ilop-ilom);
                    Real dudx = mlebtensor_dx_on_zface(i,j,k,0,vel,dxi,
                                                       whi,wlo,ihip,ihim,ilop,ilom);
                    Real dwdx = mlebtensor_dx_on_zface(i,j,k,2,vel,dxi,
                                                       whi,wlo,ihip,ihim,ilop,ilom);
                    int jhip = j + flag(i,j,k  ).isConnected(0, 1,0);
                    int jhim = j - flag(i,j,k  ).isConnected(0,-1,0);
                    int jlop = j + flag(i,j,k-1).isConnected(0, 1,0);
                    int jlom = j - flag(i,j,k-1).isConnected(0,-1,0);
                    whi = mlebtensor_weight(jhip-jhim);
                    wlo = mlebtensor_weight(jlop-jlom);
                    Real dvdy = mlebtensor_dy_on_zface(i,j,k,1,vel,dyi,
                                                       whi,wlo,jhip,jhim,jlop,jlom);
                    Real dwdy = mlebtensor_dy_on_zface(i,j,k,2,vel,dyi,
                                                       whi,wlo,jhip,jhim,jlop,jlom);
                    Real divu = dudx + dvdy;
                    Real xif = kapz(i,j,k);
                    Real mun = Real(0.75)*(etaz(i,j,k,2)-xif);// restore the original eta
                    Real mut =             etaz(i,j,k,0);

                    fz(i,j,k,0) = -mut*dwdx;
                    fz(i,j,k,1) = -mut*dwdy;
                    fz(i,j,k,2) = -mun*(-twoThirds*divu) - xif*divu;
                }
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mlebtensor_dz_on_xface (int i, int j, int k, int n,
                             Array4<Real const> const& vel, Real dzi,
                             Array4<Real const> const& bvxlo,
                             Array4<Real const> const& bvxhi,
                             Array2D<BoundCond,
                                     0,2*AMREX_SPACEDIM,
                                     0,AMREX_SPACEDIM> const& bct,
                             Dim3 const& dlo, Dim3 const& dhi,
                             Real whi, Real wlo,
                             int khip, int khim, int klop, int klom) noexcept
{
    Real ddz;
    if (i == dlo.x) {
        if (bct(Orientation::xlo(),n) == AMREX_LO_DIRICHLET && bvxlo) {
            if (k == dlo.z) {
                ddz = (bvxlo(i-1,j,k  ,n) * Real(-1.5) +
                       bvxlo(i-1,j,k+1,n) * Real(2.) +
                       bvxlo(i-1,j,k+2,n) * Real(-0.5)) * dzi;
            } else if (k == dhi.z) {
                ddz = -(bvxlo(i-1,j,k  ,n) * Real(-1.5) +
                        bvxlo(i-1,j,k-1,n) * Real(2.) +
                        bvxlo(i-1,j,k-2,n) * Real(-0.5)) * dzi;
            } else {
                ddz = wlo*dzi*(bvxlo(i-1,j,klop,n)-bvxlo(i-1,j,klom,n));
            }
        } else if (bct(Orientation::xlo(),n) == AMREX_LO_NEUMANN) {
            ddz = whi*dzi*(vel(i,j,khip,n)-vel(i,j,khim,n));
        } else { // AMREX_LO_REFLECT_ODD or homogeneous Dirichlet
            ddz = Real(0.);
        }
    } else if (i == dhi.x+1) {
        if (bct(Orientation::xhi(),n) == AMREX_LO_DIRICHLET && bvxhi) {
            if (k == dlo.z) {
                ddz = (bvxhi(i,j,k  ,n) * Real(-1.5) +
                       bvxhi(i,j,k+1,n) * Real(2.) +
                       bvxhi(i,j,k+2,n) * Real(-0.5)) * dzi;
            } else if (k == dhi.z) {
                ddz = -(bvxhi(i,j,k  ,n) * Real(-1.5) +
                        bvxhi(i,j,k-1,n) * Real(2.) +
                        bvxhi(i,j,k-2,n) * Real(-0.5)) * dzi;
            } else {
                ddz = whi*dzi*(bvxhi(i,j,khip,n)-bvxhi(i,j,khim,n));
            }
        } else if (bct(Orientation::xhi(),n) == AMREX_LO_NEUMANN) {
            ddz = wlo*dzi*(vel(i-1,j,klop,n)-vel(i-1,j,klom,n));
        } else { // AMREX_LO_REFLECT_ODD or homogeneous Dirichlet
            ddz = Real(0.);
        }
    } else {
        ddz = mlebtensor_dz_on_xface(i,j,k,n,vel,dzi,whi,wlo,khip,khim,klop,klom);
    }
    return ddz;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mlebtensor_dz_on_yface (int i, int j, int k, int n,
                             Array4<Real const> const& vel, Real dzi,
                             Array4<Real const> const& bvylo,
                             Array4<Real const> const& bvyhi,
                             Array2D<BoundCond,
                                     0,2*AMREX_SPACEDIM,
                                     0,AMREX_SPACEDIM> const& bct,
                             Dim3 const& dlo, Dim3 const& dhi,
                             Real whi, Real wlo,
                             int khip, int khim, int klop, int klom) noexcept
{
    Real ddz;
    if (j == dlo.y) {
        if (bct(Orientation::ylo(),n) == AMREX_LO_DIRICHLET && bvylo) {
            if (k == dlo.z) {
                ddz = (bvylo(i,j-1,k  ,n) * Real(-1.5) +
                       bvylo(i,j-1,k+1,n) * Real(2.) +
                       bvylo(i,j-1,k+2,n) * Real(-0.5)) * dzi;
            } else if (k == dhi.z) {
                ddz = -(bvylo(i,j-1,k  ,n) * Real(-1.5) +
                        bvylo(i,j-1,k-1,n) * Real(2.) +
                        bvylo(i,j-1,k-2,n) * Real(-0.5)) * dzi;
            } else {
                ddz = wlo*dzi*(bvylo(i,j-1,klop,n)-bvylo(i,j-1,klom,n));
            }
        } else if (bct(Orientation::ylo(),n) == AMREX_LO_NEUMANN) {
            ddz = whi*dzi*(vel(i,j,khip,n)-vel(i,j,khim,n));
        } else { // AMREX_LO_REFLECT_ODD or homogeneous Dirichlet
            ddz = Real(0.);
        }
    } else if (j == dhi.y+1) {
        if (bct(Orientation::yhi(),n) == AMREX_LO_DIRICHLET && bvyhi) {
            if (k == dlo.z) {
                ddz = (bvyhi(i,j,k  ,n) * Real(-1.5) +
                       bvyhi(i,j,k+1,n) * Real(2.) +
                       bvyhi(i,j,k+2,n) * Real(-0.5)) * dzi;
            } else if (k == dhi.z) {
                ddz = -(bvyhi(i,j,k  ,n) * Real(-1.5) +
                        bvyhi(i,j,k-1,n) * Real(2.) +
                        bvyhi(i,j,k-2,n) * Real(-0.5)) * dzi;
            } else {
                ddz = whi*dzi*(bvyhi(i,j,khip,n)-bvyhi(i,j,khim,n));
            }
        } else if (bct(Orientation::yhi(),n) == AMREX_LO_NEUMANN) {
            ddz = wlo*dzi*(vel(i,j-1,klop,n)-vel(i,j-1,klom,n));
        } else { // AMREX_LO_REFLECT_ODD or homogeneous Dirichlet
            ddz = Real(0.);
        }
    } else {
        ddz = mlebtensor_dz_on_yface(i,j,k,n,vel,dzi,whi,wlo,khip,khim,klop,klom);
    }
    return ddz;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mlebtensor_dx_on_zface (int i, int j, int k, int n,
                             Array4<Real const> const& vel, Real dxi,
                             Array4<Real const> const& bvzlo,
                             Array4<Real const> const& bvzhi,
                             Array2D<BoundCond,
                                     0,2*AMREX_SPACEDIM,
                                     0,AMREX_SPACEDIM> const& bct,
                             Dim3 const& dlo, Dim3 const& dhi,
                             Real whi, Real wlo,
                             int ihip, int ihim, int ilop, int ilom) noexcept
{
    Real ddx;
    if (k == dlo.z) {
        if (bct(Orientation::zlo(),n) == AMREX_LO_DIRICHLET && bvzlo) {
            if (i == dlo.x) {
                ddx = (bvzlo(i  ,j,k-1,n) * Real(-1.5) +
                       bvzlo(i+1,j,k-1,n) * Real(2.) +
                       bvzlo(i+2,j,k-1,n) * Real(-0.5)) * dxi;
            } else if (i == dhi.x) {
                ddx = -(bvzlo(i  ,j,k-1,n) * Real(-1.5) +
                        bvzlo(i-1,j,k-1,n) * Real(2.) +
                        bvzlo(i-2,j,k-1,n) * Real(-0.5)) * dxi;
            } else {
                ddx = wlo*dxi*(bvzlo(ilop,j,k-1,n)-bvzlo(ilom,j,k-1,n));
            }
        } else if (bct(Orientation::zlo(),n) == AMREX_LO_NEUMANN) {
            ddx = whi*dxi*(vel(ihip,j,k,n)-vel(ihim,j,k,n));
        } else { // AMREX_LO_REFLECT_ODD or homogeneous Dirichlet
            ddx = Real(0.);
        }
    } else if (k == dhi.z+1) {
        if (bct(Orientation::zhi(),n) == AMREX_LO_DIRICHLET && bvzhi) {
            if (i == dlo.x) {
                ddx = (bvzhi(i  ,j,k,n) * Real(-1.5) +
                       bvzhi(i+1,j,k,n) * Real(2.) +
                       bvzhi(i+2,j,k,n) * Real(-0.5)) * dxi;
            } else if (i == dhi.x) {
                ddx = -(bvzhi(i  ,j,k,n) * Real(-1.5) +
                        bvzhi(i-1,j,k,n) * Real(2.) +
                        bvzhi(i-2,j,k,n) * Real(-0.5)) * dxi;
            } else {
                ddx = whi*dxi*(bvzhi(ihip,j,k,n)-bvzhi(ihim,j,k,n));
            }
        } else if (bct(Orientation::zhi(),n) == AMREX_LO_NEUMANN) {
            ddx = wlo*dxi*(vel(ilop,j,k-1,n)-vel(ilom,j,k-1,n));
        } else { // AMREX_LO_REFLECT_ODD or homogeneous Dirichlet
            ddx = Real(0.);
        }
    } else {
        ddx = mlebtensor_dx_on_zface(i,j,k,n,vel,dxi,whi,wlo,ihip,ihim,ilop,ilom);

    }
    return ddx;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mlebtensor_dy_on_zface (int i, int j, int k, int n,
                             Array4<Real const> const& vel, Real dyi,
                             Array4<Real const> const& bvzlo,
                             Array4<Real const> const& bvzhi,
                             Array2D<BoundCond,
                                     0,2*AMREX_SPACEDIM,
                                     0,AMREX_SPACEDIM> const& bct,
                             Dim3 const& dlo, Dim3 const& dhi,
                             Real whi, Real wlo,
                             int jhip, int jhim, int jlop, int jlom) noexcept
{
    Real ddy;
    if (k == dlo.z) {
        if (bct(Orientation::zlo(),n) == AMREX_LO_DIRICHLET && bvzlo) {
            if (j == dlo.y) {
                ddy = (bvzlo(i,j  ,k-1,n) * Real(-1.5) +
                       bvzlo(i,j+1,k-1,n) * Real(2.) +
                       bvzlo(i,j+2,k-1,n) * Real(-0.5)) * dyi;
            } else if (j == dhi.y) {
                ddy = -(bvzlo(i,j  ,k-1,n) * Real(-1.5) +
                        bvzlo(i,j-1,k-1,n) * Real(2.) +
                        bvzlo(i,j-2,k-1,n) * Real(-0.5)) * dyi;
            } else {
                ddy = wlo*dyi*(bvzlo(i,jlop,k-1,n)-bvzlo(i,jlom,k-1,n));
            }
        } else if (bct(Orientation::zlo(),n) == AMREX_LO_NEUMANN) {
            ddy = whi*dyi*(vel(i,jhip,k,n)-vel(i,jhim,k,n));
        } else { // AMREX_LO_REFLECT_ODD or homogeneous Dirichlet
            ddy = Real(0.);
        }
    } else if (k == dhi.z+1) {
        if (bct(Orientation::zhi(),n) == AMREX_LO_DIRICHLET && bvzhi) {
            if (j == dlo.y) {
                ddy = (bvzhi(i,j  ,k,n) * Real(-1.5) +
                       bvzhi(i,j+1,k,n) * Real(2.) +
                       bvzhi(i,j+2,k,n) * Real(-0.5)) * dyi;
            } else if (j == dhi.y) {
                ddy = -(bvzhi(i,j  ,k,n) * Real(-1.5) +
                        bvzhi(i,j-1,k,n) * Real(2.) +
                        bvzhi(i,j-2,k,n) * Real(-0.5)) * dyi;
            } else {
                ddy = whi*dyi*(bvzhi(i,jhip,k,n)-bvzhi(i,jhim,k,n));
            }
        } else if (bct(Orientation::zhi(),n) == AMREX_LO_NEUMANN) {
            ddy = wlo*dyi*(vel(i,jlop,k-1,n)-vel(i,jlom,k-1,n));
        } else { // AMREX_LO_REFLECT_ODD or homogeneous Dirichlet
            ddy = Real(0.);
        }
    } else {
        ddy = mlebtensor_dy_on_zface(i,j,k,n,vel,dyi,whi,wlo,jhip,jhim,jlop,jlom);
    }
    return ddy;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebtensor_cross_terms_fx (Box const& box, Array4<Real> const& fx,
                                Array4<Real const> const& vel,
                                Array4<Real const> const& etax,
                                Array4<Real const> const& kapx,
                                Array4<Real const> const& apx,
                                Array4<EBCellFlag const> const& flag,
                                GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                                Array4<Real const> const& bvxlo,
                                Array4<Real const> const& bvxhi,
                                Array2D<BoundCond,
                                        0,2*AMREX_SPACEDIM,
                                        0,AMREX_SPACEDIM> const& bct,
                                Dim3 const& dlo, Dim3 const& dhi) noexcept

{
    const Real dyi = dxinv[1];
    const Real dzi = dxinv[2];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    constexpr Real twoThirds = 2./3.;

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                if (apx(i,j,k) == Real(0.0))
                {
                    fx(i,j,k,0) = Real(0.0);
                    fx(i,j,k,1) = Real(0.0);
                    fx(i,j,k,2) = Real(0.0);
                }
                else
                {
                    int jhip = j + flag(i  ,j,k).isConnected(0, 1,0);
                    int jhim = j - flag(i  ,j,k).isConnected(0,-1,0);
                    int jlop = j + flag(i-1,j,k).isConnected(0, 1,0);
                    int jlom = j - flag(i-1,j,k).isConnected(0,-1,0);
                    Real whi = mlebtensor_weight(jhip-jhim);
                    Real wlo = mlebtensor_weight(jlop-jlom);
                    Real dudy = mlebtensor_dy_on_xface(i,j,k,0,vel,dyi,
                                                       bvxlo,bvxhi,bct,dlo,dhi,
                                                       whi,wlo,jhip,jhim,jlop,jlom);
                    Real dvdy = mlebtensor_dy_on_xface(i,j,k,1,vel,dyi,
                                                       bvxlo,bvxhi,bct,dlo,dhi,
                                                       whi,wlo,jhip,jhim,jlop,jlom);
                    int khip = k + flag(i  ,j,k).isConnected(0,0, 1);
                    int khim = k - flag(i  ,j,k).isConnected(0,0,-1);
                    int klop = k + flag(i-1,j,k).isConnected(0,0, 1);
                    int klom = k - flag(i-1,j,k).isConnected(0,0,-1);
                    whi = mlebtensor_weight(khip-khim);
                    wlo = mlebtensor_weight(klop-klom);
                    Real dudz = mlebtensor_dz_on_xface(i,j,k,0,vel,dzi,
                                                       bvxlo,bvxhi,bct,dlo,dhi,
                                                       whi,wlo,khip,khim,klop,klom);
                    Real dwdz = mlebtensor_dz_on_xface(i,j,k,2,vel,dzi,
                                                       bvxlo,bvxhi,bct,dlo,dhi,
                                                       whi,wlo,khip,khim,klop,klom);
                    Real divu = dvdy + dwdz;
                    Real xif = kapx(i,j,k);
                    Real mun = Real(0.75)*(etax(i,j,k,0)-xif);// restore the original eta
                    Real mut =             etax(i,j,k,1);
                    fx(i,j,k,0) = -mun*(-twoThirds*divu) - xif*divu;
                    fx(i,j,k,1) = -mut*dudy;
                    fx(i,j,k,2) = -mut*dudz;
                }
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebtensor_cross_terms_fy (Box const& box, Array4<Real> const& fy,
                                Array4<Real const> const& vel,
                                Array4<Real const> const& etay,
                                Array4<Real const> const& kapy,
                                Array4<Real const> const& apy,
                                Array4<EBCellFlag const> const& flag,
                                GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                                Array4<Real const> const& bvylo,
                                Array4<Real const> const& bvyhi,
                                Array2D<BoundCond,
                                        0,2*AMREX_SPACEDIM,
                                        0,AMREX_SPACEDIM> const& bct,
                                Dim3 const& dlo, Dim3 const& dhi) noexcept
{
    const Real dxi = dxinv[0];
    const Real dzi = dxinv[2];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    constexpr Real twoThirds = 2./3.;

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                if (apy(i,j,k) == Real(0.0))
                {
                    fy(i,j,k,0) = Real(0.0);
                    fy(i,j,k,1) = Real(0.0);
                    fy(i,j,k,2) = Real(0.0);
                }
                else
                {
                    int ihip = i + flag(i,j  ,k).isConnected( 1,0,0);
                    int ihim = i - flag(i,j  ,k).isConnected(-1,0,0);
                    int ilop = i + flag(i,j-1,k).isConnected( 1,0,0);
                    int ilom = i - flag(i,j-1,k).isConnected(-1,0,0);
                    Real whi = mlebtensor_weight(ihip-ihim);
                    Real wlo = mlebtensor_weight(ilop-ilom);
                    Real dudx = mlebtensor_dx_on_yface(i,j,k,0,vel,dxi,
                                                       bvylo,bvyhi,bct,dlo,dhi,
                                                       whi,wlo,ihip,ihim,ilop,ilom);
                    Real dvdx = mlebtensor_dx_on_yface(i,j,k,1,vel,dxi,
                                                       bvylo,bvyhi,bct,dlo,dhi,
                                                       whi,wlo,ihip,ihim,ilop,ilom);
                    int khip = k + flag(i,j  ,k).isConnected(0,0, 1);
                    int khim = k - flag(i,j  ,k).isConnected(0,0,-1);
                    int klop = k + flag(i,j-1,k).isConnected(0,0, 1);
                    int klom = k - flag(i,j-1,k).isConnected(0,0,-1);
                    whi = mlebtensor_weight(khip-khim);
                    wlo = mlebtensor_weight(klop-klom);
                    Real dvdz = mlebtensor_dz_on_yface(i,j,k,1,vel,dzi,
                                                       bvylo,bvyhi,bct,dlo,dhi,
                                                       whi,wlo,khip,khim,klop,klom);
                    Real dwdz = mlebtensor_dz_on_yface(i,j,k,2,vel,dzi,
                                                       bvylo,bvyhi,bct,dlo,dhi,
                                                       whi,wlo,khip,khim,klop,klom);
                    Real divu = dudx + dwdz;
                    Real xif = kapy(i,j,k);
                    Real mun = Real(0.75)*(etay(i,j,k,1)-xif);// restore the original eta
                    Real mut =             etay(i,j,k,0);
                    fy(i,j,k,0) = -mut*dvdx;
                    fy(i,j,k,1) = -mun*(-twoThirds*divu) - xif*divu;
                    fy(i,j,k,2) = -mut*dvdz;
                }
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebtensor_cross_terms_fz (Box const& box, Array4<Real> const& fz,
                                Array4<Real const> const& vel,
                                Array4<Real const> const& etaz,
                                Array4<Real const> const& kapz,
                                Array4<Real const> const& apz,
                                Array4<EBCellFlag const> const& flag,
                                GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                                Array4<Real const> const& bvzlo,
                                Array4<Real const> const& bvzhi,
                                Array2D<BoundCond,
                                        0,2*AMREX_SPACEDIM,
                                        0,AMREX_SPACEDIM> const& bct,
                                Dim3 const& dlo, Dim3 const& dhi) noexcept
{
    const Real dxi = dxinv[0];
    const Real dyi = dxinv[1];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    constexpr Real twoThirds = 2./3.;

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                if (apz(i,j,k) == Real(0.0))
                {
                    fz(i,j,k,0) = Real(0.0);
                    fz(i,j,k,1) = Real(0.0);
                    fz(i,j,k,2) = Real(0.0);
                }
                else
                {
                    int ihip = i + flag(i,j,k  ).isConnected( 1,0,0);
                    int ihim = i - flag(i,j,k  ).isConnected(-1,0,0);
                    int ilop = i + flag(i,j,k-1).isConnected( 1,0,0);
                    int ilom = i - flag(i,j,k-1).isConnected(-1,0,0);
                    Real whi = mlebtensor_weight(ihip-ihim);
                    Real wlo = mlebtensor_weight(ilop-ilom);
                    Real dudx = mlebtensor_dx_on_zface(i,j,k,0,vel,dxi,
                                                       bvzlo,bvzhi,bct,dlo,dhi,
                                                       whi,wlo,ihip,ihim,ilop,ilom);
                    Real dwdx = mlebtensor_dx_on_zface(i,j,k,2,vel,dxi,
                                                       bvzlo,bvzhi,bct,dlo,dhi,
                                                       whi,wlo,ihip,ihim,ilop,ilom);
                    int jhip = j + flag(i,j,k  ).isConnected(0, 1,0);
                    int jhim = j - flag(i,j,k  ).isConnected(0,-1,0);
                    int jlop = j + flag(i,j,k-1).isConnected(0, 1,0);
                    int jlom = j - flag(i,j,k-1).isConnected(0,-1,0);
                    whi = mlebtensor_weight(jhip-jhim);
                    wlo = mlebtensor_weight(jlop-jlom);
                    Real dvdy = mlebtensor_dy_on_zface(i,j,k,1,vel,dyi,
                                                       bvzlo,bvzhi,bct,dlo,dhi,
                                                       whi,wlo,jhip,jhim,jlop,jlom);
                    Real dwdy = mlebtensor_dy_on_zface(i,j,k,2,vel,dyi,
                                                       bvzlo,bvzhi,bct,dlo,dhi,
                                                       whi,wlo,jhip,jhim,jlop,jlom);
                    Real divu = dudx + dvdy;
                    Real xif = kapz(i,j,k);
                    Real mun = Real(0.75)*(etaz(i,j,k,2)-xif);// restore the original eta
                    Real mut =             etaz(i,j,k,0);

                    fz(i,j,k,0) = -mut*dwdx;
                    fz(i,j,k,1) = -mut*dwdy;
                    fz(i,j,k,2) = -mun*(-twoThirds*divu) - xif*divu;
                }
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebtensor_cross_terms (Box const& box, Array4<Real> const& Ax,
                             Array4<Real const> const& fx,
                             Array4<Real const> const& fy,
                             Array4<Real const> const& fz,
                             Array4<Real const> const& vel,
                             Array4<Real const> const& velb,
                             Array4<Real const> const& etab,
                             Array4<Real const> const& kapb,
                             Array4<int const> const& ccm,
                             Array4<EBCellFlag const> const& flag,
                             Array4<Real const> const& vol,
                             Array4<Real const> const& apx,
                             Array4<Real const> const& apy,
                             Array4<Real const> const& apz,
                             Array4<Real const> const& fcx,
                             Array4<Real const> const& fcy,
                             Array4<Real const> const& fcz,
                             Array4<Real const> const& bc,
                             bool is_dirichlet, bool is_inhomog,
                             GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                             Real bscalar) noexcept
{
    const Real dxi = dxinv[0];
    const Real dyi = dxinv[1];
    const Real dzi = dxinv[2];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                if (flag(i,j,k).isRegular())
                {
                    Ax(i,j,k,0) += bscalar*(dxi*(fx(i+1,j  ,k  ,0) - fx(i,j,k,0))
                                          + dyi*(fy(i  ,j+1,k  ,0) - fy(i,j,k,0))
                                          + dzi*(fz(i  ,j  ,k+1,0) - fz(i,j,k,0)));
                    Ax(i,j,k,1) += bscalar*(dxi*(fx(i+1,j  ,k  ,1) - fx(i,j,k,1))
                                          + dyi*(fy(i  ,j+1,k  ,1) - fy(i,j,k,1))
                                          + dzi*(fz(i  ,j  ,k+1,1) - fz(i,j,k,1)));
                    Ax(i,j,k,2) += bscalar*(dxi*(fx(i+1,j  ,k  ,2) - fx(i,j,k,2))
                                          + dyi*(fy(i  ,j+1,k  ,2) - fy(i,j,k,2))
                                          + dzi*(fz(i  ,j  ,k+1,2) - fz(i,j,k,2)));
                }
                else if (flag(i,j,k).isSingleValued())
                {
                    Real fxm_0 = fx(i,j,k,0);
                    Real fxm_1 = fx(i,j,k,1);
                    Real fxm_2 = fx(i,j,k,2);
                    if (apx(i,j,k) > Real(0.0) && apx(i,j,k) < Real(1.0)) {
                        int jj = j + (fcx(i,j,k,0) >= Real(0.0) ? 1 : -1);
                        int kk = k + (fcx(i,j,k,1) >= Real(0.0) ? 1 : -1);
                        Real fracy = (ccm(i-1,jj,k) || ccm(i,jj,k)) ? std::abs(fcx(i,j,k,0)) : Real(0.0);
                        Real fracz = (ccm(i-1,j,kk) || ccm(i,j,kk)) ? std::abs(fcx(i,j,k,1)) : Real(0.0);
                        fxm_0 = (Real(1.0)-fracy)*(Real(1.0)-fracz) *fxm_0
                            + fracy*(Real(1.0)-fracz) * fx(i,jj,k ,0)
                            + fracz*(Real(1.0)-fracy) * fx(i,j ,kk,0)
                            + fracy*     fracz  * fx(i,jj,kk,0);
                        fxm_1 = (Real(1.0)-fracy)*(Real(1.0)-fracz) *fxm_1
                            + fracy*(Real(1.0)-fracz) * fx(i,jj,k ,1)
                            + fracz*(Real(1.0)-fracy) * fx(i,j ,kk,1)
                            + fracy*     fracz  * fx(i,jj,kk,1);
                        fxm_2 = (Real(1.0)-fracy)*(Real(1.0)-fracz) *fxm_2
                            + fracy*(Real(1.0)-fracz) * fx(i,jj,k ,2)
                            + fracz*(Real(1.0)-fracy) * fx(i,j ,kk,2)
                            + fracy*     fracz  * fx(i,jj,kk,2);
                    }

                    Real fxp_0 = fx(i+1,j,k,0);
                    Real fxp_1 = fx(i+1,j,k,1);
                    Real fxp_2 = fx(i+1,j,k,2);
                    if (apx(i+1,j,k) > Real(0.0) && apx(i+1,j,k) < Real(1.0)) {
                        int jj = j + (fcx(i+1,j,k,0) >= Real(0.0) ? 1 : -1);
                        int kk = k + (fcx(i+1,j,k,1) >= Real(0.0) ? 1 : -1);
                        Real fracy = (ccm(i,jj,k) || ccm(i+1,jj,k)) ? std::abs(fcx(i+1,j,k,0)) : Real(0.0);
                        Real fracz = (ccm(i,j,kk) || ccm(i+1,j,kk)) ? std::abs(fcx(i+1,j,k,1)) : Real(0.0);
                        fxp_0 = (Real(1.0)-fracy)*(Real(1.0)-fracz) *  fxp_0
                            + fracy*(Real(1.0)-fracz) * fx(i+1,jj,k ,0)
                            + fracz*(Real(1.0)-fracy) * fx(i+1,j ,kk,0)
                            + fracy*     fracz  * fx(i+1,jj,kk,0);
                        fxp_1 = (Real(1.0)-fracy)*(Real(1.0)-fracz) *fxp_1
                            + fracy*(Real(1.0)-fracz) * fx(i+1,jj,k ,1)
                            + fracz*(Real(1.0)-fracy) * fx(i+1,j ,kk,1)
                            + fracy*     fracz  * fx(i+1,jj,kk,1);
                        fxp_2 = (Real(1.0)-fracy)*(Real(1.0)-fracz) *fxp_2
                            + fracy*(Real(1.0)-fracz) * fx(i+1,jj,k ,2)
                            + fracz*(Real(1.0)-fracy) * fx(i+1,j ,kk,2)
                            + fracy*     fracz  * fx(i+1,jj,kk,2);
                    }

                    Real fym_0 = fy(i,j,k,0);
                    Real fym_1 = fy(i,j,k,1);
                    Real fym_2 = fy(i,j,k,2);
                    if (apy(i,j,k) > Real(0.0) && apy(i,j,k) < Real(1.0)) {
                        int ii = i + (fcy(i,j,k,0) >= Real(0.0) ? 1 : -1);
                        int kk = k + (fcy(i,j,k,1) >= Real(0.0) ? 1 : -1);
                        Real fracx = (ccm(ii,j-1,k) || ccm(ii,j,k)) ? std::abs(fcy(i,j,k,0)) : Real(0.0);
                        Real fracz = (ccm(i,j-1,kk) || ccm(i,j,kk)) ? std::abs(fcy(i,j,k,1)) : Real(0.0);
                        fym_0 = (Real(1.0)-fracx)*(Real(1.0)-fracz) *fym_0
                            + fracx*(Real(1.0)-fracz) * fy(ii,j,k ,0)
                            + fracz*(Real(1.0)-fracx) * fy(i ,j,kk,0)
                            + fracx*     fracz  * fy(ii,j,kk,0);
                        fym_1 = (Real(1.0)-fracx)*(Real(1.0)-fracz) *fym_1
                            + fracx*(Real(1.0)-fracz) * fy(ii,j,k ,1)
                            + fracz*(Real(1.0)-fracx) * fy(i ,j,kk,1)
                            + fracx*     fracz  * fy(ii,j,kk,1);
                        fym_2 = (Real(1.0)-fracx)*(Real(1.0)-fracz) *fym_2
                            + fracx*(Real(1.0)-fracz) * fy(ii,j,k ,2)
                            + fracz*(Real(1.0)-fracx) * fy(i ,j,kk,2)
                            + fracx*     fracz  * fy(ii,j,kk,2);
                    }

                    Real fyp_0 = fy(i,j+1,k,0);
                    Real fyp_1 = fy(i,j+1,k,1);
                    Real fyp_2 = fy(i,j+1,k,2);
                    if (apy(i,j+1,k) > Real(0.0) && apy(i,j+1,k) < Real(1.0)) {
                        int ii = i + (fcy(i,j+1,k,0) >= Real(0.0) ? 1 : -1);
                        int kk = k + (fcy(i,j+1,k,1) >= Real(0.0) ? 1 : -1);
                        Real fracx = (ccm(ii,j,k) || ccm(ii,j+1,k)) ? std::abs(fcy(i,j+1,k,0)) : Real(0.0);
                        Real fracz = (ccm(i,j,kk) || ccm(i,j+1,kk)) ? std::abs(fcy(i,j+1,k,1)) : Real(0.0);
                        fyp_0 = (Real(1.0)-fracx)*(Real(1.0)-fracz) *  fyp_0
                            + fracx*(Real(1.0)-fracz) * fy(ii,j+1,k ,0)
                            + fracz*(Real(1.0)-fracx) * fy(i ,j+1,kk,0)
                            + fracx*     fracz  * fy(ii,j+1,kk,0);
                        fyp_1 = (Real(1.0)-fracx)*(Real(1.0)-fracz) *fyp_1
                            + fracx*(Real(1.0)-fracz) * fy(ii,j+1,k ,1)
                            + fracz*(Real(1.0)-fracx) * fy(i ,j+1,kk,1)
                            + fracx*     fracz  * fy(ii,j+1,kk,1);
                        fyp_2 = (Real(1.0)-fracx)*(Real(1.0)-fracz) *fyp_2
                            + fracx*(Real(1.0)-fracz) * fy(ii,j+1,k ,2)
                            + fracz*(Real(1.0)-fracx) * fy(i ,j+1,kk,2)
                            + fracx*     fracz  * fy(ii,j+1,kk,2);
                    }

                    Real fzm_0 = fz(i,j,k,0);
                    Real fzm_1 = fz(i,j,k,1);
                    Real fzm_2 = fz(i,j,k,2);
                    if (apz(i,j,k) > Real(0.0) && apz(i,j,k) < Real(1.0)) {
                        int ii = i + (fcz(i,j,k,0) >= Real(0.0) ? 1 : -1);
                        int jj = j + (fcz(i,j,k,1) >= Real(0.0) ? 1 : -1);
                        Real fracx = (ccm(ii,j,k-1) || ccm(ii,j,k)) ? std::abs(fcz(i,j,k,0)) : Real(0.0);
                        Real fracy = (ccm(i,jj,k-1) || ccm(i,jj,k)) ? std::abs(fcz(i,j,k,1)) : Real(0.0);
                        fzm_0 = (Real(1.0)-fracx)*(Real(1.0)-fracy) *fzm_0
                            + fracx*(Real(1.0)-fracy) * fz(ii,j ,k,0)
                            + fracy*(Real(1.0)-fracx) * fz(i ,jj,k,0)
                            + fracx*     fracy  * fz(ii,jj,k,0);
                        fzm_1 = (Real(1.0)-fracx)*(Real(1.0)-fracy) *fzm_1
                            + fracx*(Real(1.0)-fracy) * fz(ii,j ,k,1)
                            + fracy*(Real(1.0)-fracx) * fz(i ,jj,k,1)
                            + fracx*     fracy  * fz(ii,jj,k,1);
                        fzm_2 = (Real(1.0)-fracx)*(Real(1.0)-fracy) *fzm_2
                            + fracx*(Real(1.0)-fracy) * fz(ii,j ,k,2)
                            + fracy*(Real(1.0)-fracx) * fz(i ,jj,k,2)
                            + fracx*     fracy  * fz(ii,jj,k,2);
                    }

                    Real fzp_0 = fz(i,j,k+1,0);
                    Real fzp_1 = fz(i,j,k+1,1);
                    Real fzp_2 = fz(i,j,k+1,2);
                    if (apz(i,j,k+1) > Real(0.0) && apz(i,j,k+1) < Real(1.0)) {
                        int ii = i + (fcz(i,j,k+1,0) >= Real(0.0) ? 1 : -1);
                        int jj = j + (fcz(i,j,k+1,1) >= Real(0.0) ? 1 : -1);
                        Real fracx = (ccm(ii,j,k) || ccm(ii,j,k)) ? std::abs(fcz(i,j,k+1,0)) : Real(0.0);
                        Real fracy = (ccm(i,jj,k) || ccm(i,jj,k)) ? std::abs(fcz(i,j,k+1,1)) : Real(0.0);
                        fzp_0 = (Real(1.0)-fracx)*(Real(1.0)-fracy) *  fzp_0
                            + fracx*(Real(1.0)-fracy) * fz(ii,j ,k+1,0)
                            + fracy*(Real(1.0)-fracx) * fz(i ,jj,k+1,0)
                            + fracx*     fracy  * fz(ii,jj,k+1,0);
                        fzp_1 = (Real(1.0)-fracx)*(Real(1.0)-fracy) *fzp_1
                            + fracx*(Real(1.0)-fracy) * fz(ii,j ,k+1,1)
                            + fracy*(Real(1.0)-fracx) * fz(i ,jj,k+1,1)
                            + fracx*     fracy  * fz(ii,jj,k+1,1);
                        fzp_2 = (Real(1.0)-fracx)*(Real(1.0)-fracy) *fzp_2
                            + fracx*(Real(1.0)-fracy) * fz(ii,j ,k+1,2)
                            + fracy*(Real(1.0)-fracx) * fz(i ,jj,k+1,2)
                            + fracx*     fracy  * fz(ii,jj,k+1,2);
                    }

                    Real kappa = vol(i,j,k);
                    Real feb_0 = Real(0.0);
                    Real feb_1 = Real(0.0);
                    Real feb_2 = Real(0.0);

                    if (is_dirichlet)
                    {
                        Real dapx = apx(i+1,j,k)-apx(i,j,k);
                        Real dapy = apy(i,j+1,k)-apy(i,j,k);
                        Real dapz = apz(i,j,k+1)-apz(i,j,k);
                        Real anorminv = Real(1.0)/std::sqrt(dapx*dapx+dapy*dapy+dapz*dapz);
                        Real anrmx = -dapx * anorminv;
                        Real anrmy = -dapy * anorminv;
                        Real anrmz = -dapz * anorminv;

                        Real velb_0 = 0, velb_1 = 0, velb_2 = 0;

                        if (is_inhomog) {
                            velb_0 = velb(i,j,k,0);
                            velb_1 = velb(i,j,k,1);
                            velb_2 = velb(i,j,k,2);
                        }

                        Real dx_eb = amrex::max(Real(0.3), (kappa*kappa-Real(0.25))/(Real(2.)*kappa));
                        Real dg = dx_eb / amrex::max(std::abs(anrmx),
                                                     std::abs(anrmy),
                                                     std::abs(anrmz));
                        Real dginv = Real(1.0)/dg;
                        Real gx = bc(i,j,k,0) - dg*anrmx;
                        Real gy = bc(i,j,k,1) - dg*anrmy;
                        Real gz = bc(i,j,k,2) - dg*anrmz;
                        int isx = (anrmx > Real(0.0)) ? 1 : -1;
                        int isy = (anrmy > Real(0.0)) ? 1 : -1;
                        int isz = (anrmz > Real(0.0)) ? 1 : -1;
                        int ii = i - isx;
                        int jj = j - isy;
                        int kk = k - isz;

                        gx *= static_cast<Real>(isx);
                        gy *= static_cast<Real>(isy);
                        gz *= static_cast<Real>(isz);
                        Real gxy = gx*gy;
                        Real gxz = gx*gz;
                        Real gyz = gy*gz;
                        Real gxyz = gx*gy*gz;
                        Real oneggg = Real(1.0)+gx+gy+gz+gxy+gxz+gyz+gxyz;

                        Real velg = oneggg             * vel(i ,j ,k ,0)
                            + (-gz - gxz - gyz - gxyz) * vel(i ,j ,kk,0)
                            + (-gy - gxy - gyz - gxyz) * vel(i ,jj,k ,0)
                            + (gyz + gxyz)             * vel(i ,jj,kk,0)
                            + (-gx - gxy - gxz - gxyz) * vel(ii,j ,k ,0)
                            + (gxz + gxyz)             * vel(ii,j ,kk,0)
                            + (gxy + gxyz)             * vel(ii,jj,k ,0)
                            + (-gxyz)                  * vel(ii,jj,kk,0);
                        Real dudn = (velb_0-velg) * dginv;

                        velg      = oneggg             * vel(i ,j ,k ,1)
                            + (-gz - gxz - gyz - gxyz) * vel(i ,j ,kk,1)
                            + (-gy - gxy - gyz - gxyz) * vel(i ,jj,k ,1)
                            + (gyz + gxyz)             * vel(i ,jj,kk,1)
                            + (-gx - gxy - gxz - gxyz) * vel(ii,j ,k ,1)
                            + (gxz + gxyz)             * vel(ii,j ,kk,1)
                            + (gxy + gxyz)             * vel(ii,jj,k ,1)
                            + (-gxyz)                  * vel(ii,jj,kk,1);
                        Real dvdn = (velb_1-velg) * dginv;

                        velg      = oneggg             * vel(i ,j ,k ,2)
                            + (-gz - gxz - gyz - gxyz) * vel(i ,j ,kk,2)
                            + (-gy - gxy - gyz - gxyz) * vel(i ,jj,k ,2)
                            + (gyz + gxyz)             * vel(i ,jj,kk,2)
                            + (-gx - gxy - gxz - gxyz) * vel(ii,j ,k ,2)
                            + (gxz + gxyz)             * vel(ii,j ,kk,2)
                            + (gxy + gxyz)             * vel(ii,jj,k ,2)
                            + (-gxyz)                  * vel(ii,jj,kk,2);
                        Real dwdn = (velb_2-velg) * dginv;

                        // transverse derivatives are zero on the wall
                        Real dudx = dudn * anrmx;
                        Real dudy = dudn * anrmy;
                        Real dudz = dudn * anrmz;
                        Real dvdx = dvdn * anrmx;
                        Real dvdy = dvdn * anrmy;
                        Real dvdz = dvdn * anrmz;
                        Real dwdx = dwdn * anrmx;
                        Real dwdy = dwdn * anrmy;
                        Real dwdz = dwdn * anrmz;
                        Real divu = dudx + dvdy + dwdz;
                        Real xi = kapb(i,j,k);
                        Real mu = etab(i,j,k);
                        Real tautmp = (xi-Real(2./3.)*mu)*divu;
                        // Note that mu*(grad v) has been included already in MLEBABecLap
                        Real tauxx = mu*dudx + tautmp;
                        Real tauyx = mu*dvdx;
                        Real tauzx = mu*dwdx;
                        Real tauyy = mu*dvdy + tautmp;
                        Real tauxy = mu*dudy;
                        Real tauzy = mu*dwdy;
                        Real tauzz = mu*dwdz + tautmp;
                        Real tauxz = mu*dudz;
                        Real tauyz = mu*dvdz;
                        // assuming dxi == dyi == dzi
                        feb_0 = dxi*(dapx*tauxx + dapy*tauyx + dapz*tauzx);
                        feb_1 = dxi*(dapx*tauxy + dapy*tauyy + dapz*tauzy);
                        feb_2 = dxi*(dapx*tauxz + dapy*tauyz + dapz*tauzz);
                    }

                    Real volinv = bscalar / kappa;
                    Ax(i,j,k,0) += volinv * (dxi*(apx(i+1,j,k)*fxp_0-apx(i,j,k)*fxm_0)
                                            +dyi*(apy(i,j+1,k)*fyp_0-apy(i,j,k)*fym_0)
                                            +dzi*(apz(i,j,k+1)*fzp_0-apz(i,j,k)*fzm_0)
                                            +dxi*                               feb_0);
                    Ax(i,j,k,1) += volinv * (dxi*(apx(i+1,j,k)*fxp_1-apx(i,j,k)*fxm_1)
                                            +dyi*(apy(i,j+1,k)*fyp_1-apy(i,j,k)*fym_1)
                                            +dzi*(apz(i,j,k+1)*fzp_1-apz(i,j,k)*fzm_1)
                                            +dxi*                               feb_1);
                    Ax(i,j,k,2) += volinv * (dxi*(apx(i+1,j,k)*fxp_2-apx(i,j,k)*fxm_2)
                                            +dyi*(apy(i,j+1,k)*fyp_2-apy(i,j,k)*fym_2)
                                            +dzi*(apz(i,j,k+1)*fzp_2-apz(i,j,k)*fzm_2)
                                            +dxi*                               feb_2);
                }
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebtensor_flux_0 (Box const& box,
                        Array4<Real> const& Ax,
                        Array4<Real const> const& fx,
                        Array4<Real const> const& ap,
                        Real bscalar) noexcept
  {
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                if (ap(i,j,k) != Real(0.0)) {
                    for (int n=0; n<AMREX_SPACEDIM; n++)
                      Ax(i,j,k,n) += bscalar*fx(i,j,k,n);
                }
                //else MLEBABec::compFlux already set Ax = 0
            }
        }
    }
}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebtensor_flux_x (Box const& box, Array4<Real> const& Ax,
                        Array4<Real const> const& fx, Array4<Real const> const& apx,
                        Array4<Real const> const& fcx,
                        Real const bscalar, Array4<int const> const& ccm,
                        int face_only, Box const& xbox) noexcept
{
    int lof = xbox.smallEnd(0);
    int hif = xbox.bigEnd(0);
    amrex::LoopConcurrent(box, [=] (int i, int j, int k) noexcept
    {
        if (!face_only || lof == i || hif == i) {
          if (apx(i,j,k) == Real(1.0)) {
            for (int n=0; n<AMREX_SPACEDIM; n++)
                Ax(i,j,k,n) += bscalar*fx(i,j,k,n);
          }
          else if (apx(i,j,k) != Real(0.0)) {
                Real fxm_0 = fx(i,j,k,0);
                Real fxm_1 = fx(i,j,k,1);
                Real fxm_2 = fx(i,j,k,2);

                int jj = j + (fcx(i,j,k,0) >= Real(0.0) ? 1 : -1);
                int kk = k + (fcx(i,j,k,1) >= Real(0.0) ? 1 : -1);
                Real fracy = (ccm(i-1,jj,k) || ccm(i,jj,k)) ? std::abs(fcx(i,j,k,0)) : Real(0.0);
                Real fracz = (ccm(i-1,j,kk) || ccm(i,j,kk)) ? std::abs(fcx(i,j,k,1)) : Real(0.0);
                fxm_0 = (Real(1.0)-fracy)*(Real(1.0)-fracz) *fxm_0
                        + fracy*(Real(1.0)-fracz) * fx(i,jj,k ,0)
                        + fracz*(Real(1.0)-fracy) * fx(i,j ,kk,0)
                        + fracy*     fracz  * fx(i,jj,kk,0);
                fxm_1 = (Real(1.0)-fracy)*(Real(1.0)-fracz) *fxm_1
                        + fracy*(Real(1.0)-fracz) * fx(i,jj,k ,1)
                        + fracz*(Real(1.0)-fracy) * fx(i,j ,kk,1)
                        + fracy*     fracz  * fx(i,jj,kk,1);
                fxm_2 = (Real(1.0)-fracy)*(Real(1.0)-fracz) *fxm_2
                        + fracy*(Real(1.0)-fracz) * fx(i,jj,k ,2)
                        + fracz*(Real(1.0)-fracy) * fx(i,j ,kk,2)
                        + fracy*     fracz  * fx(i,jj,kk,2);

                Ax(i,j,k,0) += bscalar*fxm_0;
                Ax(i,j,k,1) += bscalar*fxm_1;
                Ax(i,j,k,2) += bscalar*fxm_2;
          }
        }
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebtensor_flux_y (Box const& box, Array4<Real> const& Ay,
                        Array4<Real const> const& fy, Array4<Real const> const& apy,
                        Array4<Real const> const& fcy,
                        Real const bscalar, Array4<int const> const& ccm,
                        int face_only, Box const& ybox) noexcept
{
    int lof = ybox.smallEnd(1);
    int hif = ybox.bigEnd(1);
    amrex::LoopConcurrent(box, [=] (int i, int j, int k) noexcept
    {
        if (!face_only || lof == j || hif == j) {
          if (apy(i,j,k) == Real(1.0)) {
            for (int n=0; n<AMREX_SPACEDIM; n++)
                  Ay(i,j,k,n) += bscalar*fy(i,j,k,n);
          } else if (apy(i,j,k) != 0) {
                  Real fym_0 = fy(i,j,k,0);
                  Real fym_1 = fy(i,j,k,1);
                  Real fym_2 = fy(i,j,k,2);

                  int ii = i + (fcy(i,j,k,0) >= Real(0.0) ? 1 : -1);
                  int kk = k + (fcy(i,j,k,1) >= Real(0.0) ? 1 : -1);
                  Real fracx = (ccm(ii,j-1,k) || ccm(ii,j,k)) ? std::abs(fcy(i,j,k,0)) : Real(0.0);
                  Real fracz = (ccm(i,j-1,kk) || ccm(i,j,kk)) ? std::abs(fcy(i,j,k,1)) : Real(0.0);
                  fym_0 = (Real(1.0)-fracx)*(Real(1.0)-fracz) *fym_0
                          + fracx*(Real(1.0)-fracz) * fy(ii,j,k ,0)
                          + fracz*(Real(1.0)-fracx) * fy(i ,j,kk,0)
                          + fracx*     fracz  * fy(ii,j,kk,0);
                  fym_1 = (Real(1.0)-fracx)*(Real(1.0)-fracz) *fym_1
                          + fracx*(Real(1.0)-fracz) * fy(ii,j,k ,1)
                          + fracz*(Real(1.0)-fracx) * fy(i ,j,kk,1)
                          + fracx*     fracz  * fy(ii,j,kk,1);
                  fym_2 = (Real(1.0)-fracx)*(Real(1.0)-fracz) *fym_2
                          + fracx*(Real(1.0)-fracz) * fy(ii,j,k ,2)
                          + fracz*(Real(1.0)-fracx) * fy(i ,j,kk,2)
                          + fracx*     fracz  * fy(ii,j,kk,2);

                  Ay(i,j,k,0) += bscalar*fym_0;
                  Ay(i,j,k,1) += bscalar*fym_1;
                  Ay(i,j,k,2) += bscalar*fym_2;
          }
          //else MLEBABec::compFlux already set Ay = 0
        }
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebtensor_flux_z (Box const& box, Array4<Real> const& Az,
                        Array4<Real const> const& fz, Array4<Real const> const& apz,
                        Array4<Real const> const& fcz,
                        Real const bscalar, Array4<int const> const& ccm,
                        int face_only, Box const& zbox) noexcept
{
    int lof = zbox.smallEnd(2);
    int hif = zbox.bigEnd(2);
    amrex::LoopConcurrent(box, [=] (int i, int j, int k) noexcept
    {
        if (!face_only || lof == k || hif == k) {
          if (apz(i,j,k) == Real(1.0)) {
            for (int n=0; n<AMREX_SPACEDIM; n++)
                  Az(i,j,k,n) += bscalar*fz(i,j,k,n);
          }
          else if (apz(i,j,k) != 0.) {
                  Real fzm_0 = fz(i,j,k,0);
                  Real fzm_1 = fz(i,j,k,1);
                  Real fzm_2 = fz(i,j,k,2);

                  int ii = i + (fcz(i,j,k,0) >= Real(0.0) ? 1 : -1);
                  int jj = j + (fcz(i,j,k,1) >= Real(0.0) ? 1 : -1);
                  Real fracx = (ccm(ii,j,k-1) || ccm(ii,j,k)) ? std::abs(fcz(i,j,k,0)) : Real(0.0);
                  Real fracy = (ccm(i,jj,k-1) || ccm(i,jj,k)) ? std::abs(fcz(i,j,k,1)) : Real(0.0);
                  fzm_0 = (Real(1.0)-fracx)*(Real(1.0)-fracy) *fzm_0
                          + fracx*(Real(1.0)-fracy) * fz(ii,j ,k,0)
                          + fracy*(Real(1.0)-fracx) * fz(i ,jj,k,0)
                          + fracx*     fracy  * fz(ii,jj,k,0);
                  fzm_1 = (Real(1.0)-fracx)*(Real(1.0)-fracy) *fzm_1
                          + fracx*(Real(1.0)-fracy) * fz(ii,j ,k,1)
                          + fracy*(Real(1.0)-fracx) * fz(i ,jj,k,1)
                          + fracx*     fracy  * fz(ii,jj,k,1);
                  fzm_2 = (Real(1.0)-fracx)*(Real(1.0)-fracy) *fzm_2
                          + fracx*(Real(1.0)-fracy) * fz(ii,j ,k,2)
                          + fracy*(Real(1.0)-fracx) * fz(i ,jj,k,2)
                          + fracx*     fracy  * fz(ii,jj,k,2);

                  Az(i,j,k,0) += bscalar*fzm_0;
                  Az(i,j,k,1) += bscalar*fzm_1;
                  Az(i,j,k,2) += bscalar*fzm_2;
          }
          //else MLEBABec::compFlux already set Az = 0
        }
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebtensor_vel_grads_fx (Box const& box, Array4<Real> const& fx,
                              Array4<Real const> const& vel, Array4<Real const> const& apx,
                              Array4<EBCellFlag const> const& flag,
                              GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    const Real dxi = dxinv[0];
    const Real dyi = dxinv[1];
    const Real dzi = dxinv[2];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                if (apx(i,j,k) == Real(0.0))
                {
                    fx(i,j,k,0) = Real(0.0);
                    fx(i,j,k,1) = Real(0.0);
                    fx(i,j,k,2) = Real(0.0);
                    fx(i,j,k,3) = Real(0.0);
                    fx(i,j,k,4) = Real(0.0);
                    fx(i,j,k,5) = Real(0.0);
                    fx(i,j,k,6) = Real(0.0);
                    fx(i,j,k,7) = Real(0.0);
                    fx(i,j,k,8) = Real(0.0);
                }
                else
                {
                    Real dudx = (vel(i,j,k,0) - vel(i-1,j,k,0))*dxi;
                    Real dvdx = (vel(i,j,k,1) - vel(i-1,j,k,1))*dxi;
                    Real dwdx = (vel(i,j,k,2) - vel(i-1,j,k,2))*dxi;

                    int jhip = j + flag(i  ,j,k).isConnected(0, 1,0);
                    int jhim = j - flag(i  ,j,k).isConnected(0,-1,0);
                    int jlop = j + flag(i-1,j,k).isConnected(0, 1,0);
                    int jlom = j - flag(i-1,j,k).isConnected(0,-1,0);
                    Real whi = mlebtensor_weight(jhip-jhim);
                    Real wlo = mlebtensor_weight(jlop-jlom);
                    Real dudy = mlebtensor_dy_on_xface(i,j,k,0,vel,dyi,
                                                       whi,wlo,jhip,jhim,jlop,jlom);
                    Real dvdy = mlebtensor_dy_on_xface(i,j,k,1,vel,dyi,
                                                       whi,wlo,jhip,jhim,jlop,jlom);
                    Real dwdy = mlebtensor_dy_on_xface(i,j,k,2,vel,dyi,
                                                       whi,wlo,jhip,jhim,jlop,jlom);

                    int khip = k + flag(i  ,j,k).isConnected(0,0, 1);
                    int khim = k - flag(i  ,j,k).isConnected(0,0,-1);
                    int klop = k + flag(i-1,j,k).isConnected(0,0, 1);
                    int klom = k - flag(i-1,j,k).isConnected(0,0,-1);
                    whi = mlebtensor_weight(khip-khim);
                    wlo = mlebtensor_weight(klop-klom);
                    Real dudz = mlebtensor_dz_on_xface(i,j,k,0,vel,dzi,
                                                       whi,wlo,khip,khim,klop,klom);
                    Real dvdz = mlebtensor_dz_on_xface(i,j,k,1,vel,dzi,
                                                       whi,wlo,khip,khim,klop,klom);
                    Real dwdz = mlebtensor_dz_on_xface(i,j,k,2,vel,dzi,
                                                       whi,wlo,khip,khim,klop,klom);

                    fx(i,j,k,0) = dudx;
                    fx(i,j,k,1) = dvdx;
                    fx(i,j,k,2) = dwdx;
                    fx(i,j,k,3) = dudy;
                    fx(i,j,k,4) = dvdy;
                    fx(i,j,k,5) = dwdy;
                    fx(i,j,k,6) = dudz;
                    fx(i,j,k,7) = dvdz;
                    fx(i,j,k,8) = dwdz;
                }
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebtensor_vel_grads_fy (Box const& box, Array4<Real> const& fy,
                              Array4<Real const> const& vel, Array4<Real const> const& apy,
                              Array4<EBCellFlag const> const& flag,
                              GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    const Real dxi = dxinv[0];
    const Real dyi = dxinv[1];
    const Real dzi = dxinv[2];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                if (apy(i,j,k) == Real(0.0))
                {
                    fy(i,j,k,0) = Real(0.0);
                    fy(i,j,k,1) = Real(0.0);
                    fy(i,j,k,2) = Real(0.0);
                    fy(i,j,k,3) = Real(0.0);
                    fy(i,j,k,4) = Real(0.0);
                    fy(i,j,k,5) = Real(0.0);
                    fy(i,j,k,6) = Real(0.0);
                    fy(i,j,k,7) = Real(0.0);
                    fy(i,j,k,8) = Real(0.0);
                }
                else
                {
                    int ihip = i + flag(i,j  ,k).isConnected( 1,0,0);
                    int ihim = i - flag(i,j  ,k).isConnected(-1,0,0);
                    int ilop = i + flag(i,j-1,k).isConnected( 1,0,0);
                    int ilom = i - flag(i,j-1,k).isConnected(-1,0,0);
                    Real whi = mlebtensor_weight(ihip-ihim);
                    Real wlo = mlebtensor_weight(ilop-ilom);
                    Real dudx = mlebtensor_dx_on_yface(i,j,k,0,vel,dxi,
                                                       whi,wlo,ihip,ihim,ilop,ilom);
                    Real dvdx = mlebtensor_dx_on_yface(i,j,k,1,vel,dxi,
                                                       whi,wlo,ihip,ihim,ilop,ilom);
                    Real dwdx = mlebtensor_dx_on_yface(i,j,k,2,vel,dxi,
                                                       whi,wlo,ihip,ihim,ilop,ilom);

                    Real dudy = (vel(i,j,k,0) - vel(i,j-1,k,0))*dyi;
                    Real dvdy = (vel(i,j,k,1) - vel(i,j-1,k,1))*dyi;
                    Real dwdy = (vel(i,j,k,2) - vel(i,j-1,k,2))*dyi;

                    int khip = k + flag(i,j  ,k).isConnected(0,0, 1);
                    int khim = k - flag(i,j  ,k).isConnected(0,0,-1);
                    int klop = k + flag(i,j-1,k).isConnected(0,0, 1);
                    int klom = k - flag(i,j-1,k).isConnected(0,0,-1);
                    whi = mlebtensor_weight(khip-khim);
                    wlo = mlebtensor_weight(klop-klom);
                    Real dudz = mlebtensor_dz_on_yface(i,j,k,0,vel,dzi,
                                                       whi,wlo,khip,khim,klop,klom);
                    Real dvdz = mlebtensor_dz_on_yface(i,j,k,1,vel,dzi,
                                                       whi,wlo,khip,khim,klop,klom);
                    Real dwdz = mlebtensor_dz_on_yface(i,j,k,2,vel,dzi,
                                                       whi,wlo,khip,khim,klop,klom);

                    fy(i,j,k,0) = dudx;
                    fy(i,j,k,1) = dvdx;
                    fy(i,j,k,2) = dwdx;
                    fy(i,j,k,3) = dudy;
                    fy(i,j,k,4) = dvdy;
                    fy(i,j,k,5) = dwdy;
                    fy(i,j,k,6) = dudz;
                    fy(i,j,k,7) = dvdz;
                    fy(i,j,k,8) = dwdz;
                }
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebtensor_vel_grads_fz (Box const& box, Array4<Real> const& fz,
                              Array4<Real const> const& vel, Array4<Real const> const& apz,
                              Array4<EBCellFlag const> const& flag,
                              GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    const Real dxi = dxinv[0];
    const Real dyi = dxinv[1];
    const Real dzi = dxinv[2];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                if (apz(i,j,k) == Real(0.0))
                {
                    fz(i,j,k,0) = Real(0.0);
                    fz(i,j,k,1) = Real(0.0);
                    fz(i,j,k,2) = Real(0.0);
                    fz(i,j,k,3) = Real(0.0);
                    fz(i,j,k,4) = Real(0.0);
                    fz(i,j,k,5) = Real(0.0);
                    fz(i,j,k,6) = Real(0.0);
                    fz(i,j,k,7) = Real(0.0);
                    fz(i,j,k,8) = Real(0.0);
                }
                else
                {
                    int ihip = i + flag(i,j,k  ).isConnected( 1,0,0);
                    int ihim = i - flag(i,j,k  ).isConnected(-1,0,0);
                    int ilop = i + flag(i,j,k-1).isConnected( 1,0,0);
                    int ilom = i - flag(i,j,k-1).isConnected(-1,0,0);
                    Real whi = mlebtensor_weight(ihip-ihim);
                    Real wlo = mlebtensor_weight(ilop-ilom);
                    Real dudx = mlebtensor_dx_on_zface(i,j,k,0,vel,dxi,
                                                       whi,wlo,ihip,ihim,ilop,ilom);
                    Real dvdx = mlebtensor_dx_on_zface(i,j,k,1,vel,dxi,
                                                       whi,wlo,ihip,ihim,ilop,ilom);
                    Real dwdx = mlebtensor_dx_on_zface(i,j,k,2,vel,dxi,
                                                       whi,wlo,ihip,ihim,ilop,ilom);

                    int jhip = j + flag(i,j,k  ).isConnected(0, 1,0);
                    int jhim = j - flag(i,j,k  ).isConnected(0,-1,0);
                    int jlop = j + flag(i,j,k-1).isConnected(0, 1,0);
                    int jlom = j - flag(i,j,k-1).isConnected(0,-1,0);
                    whi = mlebtensor_weight(jhip-jhim);
                    wlo = mlebtensor_weight(jlop-jlom);
                    Real dudy = mlebtensor_dy_on_zface(i,j,k,0,vel,dyi,
                                                       whi,wlo,jhip,jhim,jlop,jlom);
                    Real dvdy = mlebtensor_dy_on_zface(i,j,k,1,vel,dyi,
                                                       whi,wlo,jhip,jhim,jlop,jlom);
                    Real dwdy = mlebtensor_dy_on_zface(i,j,k,2,vel,dyi,
                                                       whi,wlo,jhip,jhim,jlop,jlom);

                    Real dudz = (vel(i,j,k,0) - vel(i,j,k-1,0))*dzi;
                    Real dvdz = (vel(i,j,k,1) - vel(i,j,k-1,1))*dzi;
                    Real dwdz = (vel(i,j,k,2) - vel(i,j,k-1,2))*dzi;

                    fz(i,j,k,0) = dudx;
                    fz(i,j,k,1) = dvdx;
                    fz(i,j,k,2) = dwdx;
                    fz(i,j,k,3) = dudy;
                    fz(i,j,k,4) = dvdy;
                    fz(i,j,k,5) = dwdy;
                    fz(i,j,k,6) = dudz;
                    fz(i,j,k,7) = dvdz;
                    fz(i,j,k,8) = dwdz;
                }
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebtensor_vel_grads_fx (Box const& box, Array4<Real> const& fx,
                              Array4<Real const> const& vel,
                              Array4<Real const> const& apx,
                              Array4<EBCellFlag const> const& flag,
                              GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                              Array4<Real const> const& bvxlo,
                              Array4<Real const> const& bvxhi,
                              Array2D<BoundCond,
                                      0,2*AMREX_SPACEDIM,
                                      0,AMREX_SPACEDIM> const& bct,
                              Dim3 const& dlo, Dim3 const& dhi) noexcept
{
    const Real dxi = dxinv[0];
    const Real dyi = dxinv[1];
    const Real dzi = dxinv[2];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                if (apx(i,j,k) == Real(0.0))
                {
                    fx(i,j,k,0) = Real(0.0);
                    fx(i,j,k,1) = Real(0.0);
                    fx(i,j,k,2) = Real(0.0);
                    fx(i,j,k,3) = Real(0.0);
                    fx(i,j,k,4) = Real(0.0);
                    fx(i,j,k,5) = Real(0.0);
                    fx(i,j,k,6) = Real(0.0);
                    fx(i,j,k,7) = Real(0.0);
                    fx(i,j,k,8) = Real(0.0);
                }
                else
                {
                    Real dudx = (vel(i,j,k,0) - vel(i-1,j,k,0))*dxi;
                    Real dvdx = (vel(i,j,k,1) - vel(i-1,j,k,1))*dxi;
                    Real dwdx = (vel(i,j,k,2) - vel(i-1,j,k,2))*dxi;

                    int jhip = j + flag(i  ,j,k).isConnected(0, 1,0);
                    int jhim = j - flag(i  ,j,k).isConnected(0,-1,0);
                    int jlop = j + flag(i-1,j,k).isConnected(0, 1,0);
                    int jlom = j - flag(i-1,j,k).isConnected(0,-1,0);
                    Real whi = mlebtensor_weight(jhip-jhim);
                    Real wlo = mlebtensor_weight(jlop-jlom);
                    Real dudy = mlebtensor_dy_on_xface(i,j,k,0,vel,dyi,
                                                       bvxlo,bvxhi,bct,dlo,dhi,
                                                       whi,wlo,jhip,jhim,jlop,jlom);
                    Real dvdy = mlebtensor_dy_on_xface(i,j,k,1,vel,dyi,
                                                       bvxlo,bvxhi,bct,dlo,dhi,
                                                       whi,wlo,jhip,jhim,jlop,jlom);
                    Real dwdy = mlebtensor_dy_on_xface(i,j,k,2,vel,dyi,
                                                       bvxlo,bvxhi,bct,dlo,dhi,
                                                       whi,wlo,jhip,jhim,jlop,jlom);

                    int khip = k + flag(i  ,j,k).isConnected(0,0, 1);
                    int khim = k - flag(i  ,j,k).isConnected(0,0,-1);
                    int klop = k + flag(i-1,j,k).isConnected(0,0, 1);
                    int klom = k - flag(i-1,j,k).isConnected(0,0,-1);
                    whi = mlebtensor_weight(khip-khim);
                    wlo = mlebtensor_weight(klop-klom);
                    Real dudz = mlebtensor_dz_on_xface(i,j,k,0,vel,dzi,
                                                       bvxlo,bvxhi,bct,dlo,dhi,
                                                       whi,wlo,khip,khim,klop,klom);
                    Real dvdz = mlebtensor_dz_on_xface(i,j,k,1,vel,dzi,
                                                       bvxlo,bvxhi,bct,dlo,dhi,
                                                       whi,wlo,khip,khim,klop,klom);
                    Real dwdz = mlebtensor_dz_on_xface(i,j,k,2,vel,dzi,
                                                       bvxlo,bvxhi,bct,dlo,dhi,
                                                       whi,wlo,khip,khim,klop,klom);

                    fx(i,j,k,0) = dudx;
                    fx(i,j,k,1) = dvdx;
                    fx(i,j,k,2) = dwdx;
                    fx(i,j,k,3) = dudy;
                    fx(i,j,k,4) = dvdy;
                    fx(i,j,k,5) = dwdy;
                    fx(i,j,k,6) = dudz;
                    fx(i,j,k,7) = dvdz;
                    fx(i,j,k,8) = dwdz;
                }
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebtensor_vel_grads_fy (Box const& box, Array4<Real> const& fy,
                              Array4<Real const> const& vel,
                              Array4<Real const> const& apy,
                              Array4<EBCellFlag const> const& flag,
                              GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                              Array4<Real const> const& bvylo,
                              Array4<Real const> const& bvyhi,
                              Array2D<BoundCond,
                                      0,2*AMREX_SPACEDIM,
                                      0,AMREX_SPACEDIM> const& bct,
                              Dim3 const& dlo, Dim3 const& dhi) noexcept
{
    const Real dxi = dxinv[0];
    const Real dyi = dxinv[1];
    const Real dzi = dxinv[2];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                if (apy(i,j,k) == Real(0.0))
                {
                    fy(i,j,k,0) = Real(0.0);
                    fy(i,j,k,1) = Real(0.0);
                    fy(i,j,k,2) = Real(0.0);
                    fy(i,j,k,3) = Real(0.0);
                    fy(i,j,k,4) = Real(0.0);
                    fy(i,j,k,5) = Real(0.0);
                    fy(i,j,k,6) = Real(0.0);
                    fy(i,j,k,7) = Real(0.0);
                    fy(i,j,k,8) = Real(0.0);
                }
                else
                {
                    int ihip = i + flag(i,j  ,k).isConnected( 1,0,0);
                    int ihim = i - flag(i,j  ,k).isConnected(-1,0,0);
                    int ilop = i + flag(i,j-1,k).isConnected( 1,0,0);
                    int ilom = i - flag(i,j-1,k).isConnected(-1,0,0);
                    Real whi = mlebtensor_weight(ihip-ihim);
                    Real wlo = mlebtensor_weight(ilop-ilom);
                    Real dudx = mlebtensor_dx_on_yface(i,j,k,0,vel,dxi,
                                                       bvylo,bvyhi,bct,dlo,dhi,
                                                       whi,wlo,ihip,ihim,ilop,ilom);
                    Real dvdx = mlebtensor_dx_on_yface(i,j,k,1,vel,dxi,
                                                       bvylo,bvyhi,bct,dlo,dhi,
                                                       whi,wlo,ihip,ihim,ilop,ilom);
                    Real dwdx = mlebtensor_dx_on_yface(i,j,k,2,vel,dxi,
                                                       bvylo,bvyhi,bct,dlo,dhi,
                                                       whi,wlo,ihip,ihim,ilop,ilom);

                    Real dudy = (vel(i,j,k,0) - vel(i,j-1,k,0))*dyi;
                    Real dvdy = (vel(i,j,k,1) - vel(i,j-1,k,1))*dyi;
                    Real dwdy = (vel(i,j,k,2) - vel(i,j-1,k,2))*dyi;

                    int khip = k + flag(i,j  ,k).isConnected(0,0, 1);
                    int khim = k - flag(i,j  ,k).isConnected(0,0,-1);
                    int klop = k + flag(i,j-1,k).isConnected(0,0, 1);
                    int klom = k - flag(i,j-1,k).isConnected(0,0,-1);
                    whi = mlebtensor_weight(khip-khim);
                    wlo = mlebtensor_weight(klop-klom);
                    Real dudz = mlebtensor_dz_on_yface(i,j,k,0,vel,dzi,
                                                       bvylo,bvyhi,bct,dlo,dhi,
                                                       whi,wlo,khip,khim,klop,klom);
                    Real dvdz = mlebtensor_dz_on_yface(i,j,k,1,vel,dzi,
                                                       bvylo,bvyhi,bct,dlo,dhi,
                                                       whi,wlo,khip,khim,klop,klom);
                    Real dwdz = mlebtensor_dz_on_yface(i,j,k,2,vel,dzi,
                                                       bvylo,bvyhi,bct,dlo,dhi,
                                                       whi,wlo,khip,khim,klop,klom);

                    fy(i,j,k,0) = dudx;
                    fy(i,j,k,1) = dvdx;
                    fy(i,j,k,2) = dwdx;
                    fy(i,j,k,3) = dudy;
                    fy(i,j,k,4) = dvdy;
                    fy(i,j,k,5) = dwdy;
                    fy(i,j,k,6) = dudz;
                    fy(i,j,k,7) = dvdz;
                    fy(i,j,k,8) = dwdz;
                }
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebtensor_vel_grads_fz (Box const& box, Array4<Real> const& fz,
                              Array4<Real const> const& vel,
                              Array4<Real const> const& apz,
                              Array4<EBCellFlag const> const& flag,
                              GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                              Array4<Real const> const& bvzlo,
                              Array4<Real const> const& bvzhi,
                              Array2D<BoundCond,
                                      0,2*AMREX_SPACEDIM,
                                      0,AMREX_SPACEDIM> const& bct,
                              Dim3 const& dlo, Dim3 const& dhi) noexcept
{
    const Real dxi = dxinv[0];
    const Real dyi = dxinv[1];
    const Real dzi = dxinv[2];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                if (apz(i,j,k) == Real(0.0))
                {
                    fz(i,j,k,0) = Real(0.0);
                    fz(i,j,k,1) = Real(0.0);
                    fz(i,j,k,2) = Real(0.0);
                    fz(i,j,k,3) = Real(0.0);
                    fz(i,j,k,4) = Real(0.0);
                    fz(i,j,k,5) = Real(0.0);
                    fz(i,j,k,6) = Real(0.0);
                    fz(i,j,k,7) = Real(0.0);
                    fz(i,j,k,8) = Real(0.0);
                }
                else
                {
                    int ihip = i + flag(i,j,k  ).isConnected( 1,0,0);
                    int ihim = i - flag(i,j,k  ).isConnected(-1,0,0);
                    int ilop = i + flag(i,j,k-1).isConnected( 1,0,0);
                    int ilom = i - flag(i,j,k-1).isConnected(-1,0,0);
                    Real whi = mlebtensor_weight(ihip-ihim);
                    Real wlo = mlebtensor_weight(ilop-ilom);
                    Real dudx = mlebtensor_dx_on_zface(i,j,k,0,vel,dxi,
                                                       bvzlo,bvzhi,bct,dlo,dhi,
                                                       whi,wlo,ihip,ihim,ilop,ilom);
                    Real dvdx = mlebtensor_dx_on_zface(i,j,k,1,vel,dxi,
                                                       bvzlo,bvzhi,bct,dlo,dhi,
                                                       whi,wlo,ihip,ihim,ilop,ilom);
                    Real dwdx = mlebtensor_dx_on_zface(i,j,k,2,vel,dxi,
                                                       bvzlo,bvzhi,bct,dlo,dhi,
                                                       whi,wlo,ihip,ihim,ilop,ilom);

                    int jhip = j + flag(i,j,k  ).isConnected(0, 1,0);
                    int jhim = j - flag(i,j,k  ).isConnected(0,-1,0);
                    int jlop = j + flag(i,j,k-1).isConnected(0, 1,0);
                    int jlom = j - flag(i,j,k-1).isConnected(0,-1,0);
                    whi = mlebtensor_weight(jhip-jhim);
                    wlo = mlebtensor_weight(jlop-jlom);
                    Real dudy = mlebtensor_dy_on_zface(i,j,k,0,vel,dyi,
                                                       bvzlo,bvzhi,bct,dlo,dhi,
                                                       whi,wlo,jhip,jhim,jlop,jlom);
                    Real dvdy = mlebtensor_dy_on_zface(i,j,k,1,vel,dyi,
                                                       bvzlo,bvzhi,bct,dlo,dhi,
                                                       whi,wlo,jhip,jhim,jlop,jlom);
                    Real dwdy = mlebtensor_dy_on_zface(i,j,k,2,vel,dyi,
                                                       bvzlo,bvzhi,bct,dlo,dhi,
                                                       whi,wlo,jhip,jhim,jlop,jlom);

                    Real dudz = (vel(i,j,k,0) - vel(i,j,k-1,0))*dzi;
                    Real dvdz = (vel(i,j,k,1) - vel(i,j,k-1,1))*dzi;
                    Real dwdz = (vel(i,j,k,2) - vel(i,j,k-1,2))*dzi;

                    fz(i,j,k,0) = dudx;
                    fz(i,j,k,1) = dvdx;
                    fz(i,j,k,2) = dwdx;
                    fz(i,j,k,3) = dudy;
                    fz(i,j,k,4) = dvdy;
                    fz(i,j,k,5) = dwdy;
                    fz(i,j,k,6) = dudz;
                    fz(i,j,k,7) = dvdz;
                    fz(i,j,k,8) = dwdz;
                }
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebtensor_vel_grads_fx (Box const& box, Array4<Real> const& fx,
                              Array4<Real const> const& vel, Array4<Real const> const& apx,
                              Array4<EBCellFlag const> const& flag,
                              Array4<int const> const& ccm,
                              Array4<Real const> const& fcx,
                              GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    const Real dxi = dxinv[0];
    const Real dyi = dxinv[1];
    const Real dzi = dxinv[2];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                if (apx(i,j,k) == Real(0.0))
                {
                    fx(i,j,k,0) = Real(0.0);
                    fx(i,j,k,1) = Real(0.0);
                    fx(i,j,k,2) = Real(0.0);
                    fx(i,j,k,3) = Real(0.0);
                    fx(i,j,k,4) = Real(0.0);
                    fx(i,j,k,5) = Real(0.0);
                    fx(i,j,k,6) = Real(0.0);
                    fx(i,j,k,7) = Real(0.0);
                    fx(i,j,k,8) = Real(0.0);
                }
                else
                {
                    Real dudx = (vel(i,j,k,0) - vel(i-1,j,k,0))*dxi;
                    Real dvdx = (vel(i,j,k,1) - vel(i-1,j,k,1))*dxi;
                    Real dwdx = (vel(i,j,k,2) - vel(i-1,j,k,2))*dxi;
                    if (apx(i,j,k) < Real(1.0)) {
                        int jj = j + (fcx(i,j,k,0) >= Real(0.0) ? 1 : -1);
                        int kk = k + (fcx(i,j,k,1) >= Real(0.0) ? 1 : -1);
                        Real fracy = (ccm(i-1,jj,k) || ccm(i,jj,k)) ? std::abs(fcx(i,j,k,0)) : Real(0.0);
                        Real fracz = (ccm(i-1,j,kk) || ccm(i,j,kk)) ? std::abs(fcx(i,j,k,1)) : Real(0.0);
                        dudx = (Real(1.0)-fracy)*(Real(1.0)-fracz) * dudx
                            + fracy*(Real(1.0)-fracz) * (vel(i,jj,k ,0) - vel(i-1,jj,k ,0))*dxi
                            + fracz*(Real(1.0)-fracy) * (vel(i,j ,kk,0) - vel(i-1,j ,kk,0))*dxi
                            + fracy*     fracz  * (vel(i,jj,kk,0) - vel(i-1,jj,kk,0))*dxi;
                        dvdx = (Real(1.0)-fracy)*(Real(1.0)-fracz) * dvdx
                            + fracy*(Real(1.0)-fracz) * (vel(i,jj,k ,1) - vel(i-1,jj,k ,1))*dxi
                            + fracz*(Real(1.0)-fracy) * (vel(i,j ,kk,1) - vel(i-1,j ,kk,1))*dxi
                            + fracy*     fracz  * (vel(i,jj,kk,1) - vel(i-1,jj,kk,1))*dxi;
                        dwdx = (Real(1.0)-fracy)*(Real(1.0)-fracz) * dwdx
                            + fracy*(Real(1.0)-fracz) * (vel(i,jj,k ,2) - vel(i-1,jj,k ,2))*dxi
                            + fracz*(Real(1.0)-fracy) * (vel(i,j ,kk,2) - vel(i-1,j ,kk,2))*dxi
                            + fracy*     fracz  * (vel(i,jj,kk,2) - vel(i-1,jj,kk,2))*dxi;
                    }

                    int jhip = j + flag(i  ,j,k).isConnected(0, 1,0);
                    int jhim = j - flag(i  ,j,k).isConnected(0,-1,0);
                    int jlop = j + flag(i-1,j,k).isConnected(0, 1,0);
                    int jlom = j - flag(i-1,j,k).isConnected(0,-1,0);
                    Real whi = mlebtensor_weight(jhip-jhim);
                    Real wlo = mlebtensor_weight(jlop-jlom);
                    Real dudy = mlebtensor_dy_on_xface(i,j,k,0,vel,dyi,
                                                       whi,wlo,jhip,jhim,jlop,jlom);
                    Real dvdy = mlebtensor_dy_on_xface(i,j,k,1,vel,dyi,
                                                       whi,wlo,jhip,jhim,jlop,jlom);
                    Real dwdy = mlebtensor_dy_on_xface(i,j,k,2,vel,dyi,
                                                       whi,wlo,jhip,jhim,jlop,jlom);

                    int khip = k + flag(i  ,j,k).isConnected(0,0, 1);
                    int khim = k - flag(i  ,j,k).isConnected(0,0,-1);
                    int klop = k + flag(i-1,j,k).isConnected(0,0, 1);
                    int klom = k - flag(i-1,j,k).isConnected(0,0,-1);
                    whi = mlebtensor_weight(khip-khim);
                    wlo = mlebtensor_weight(klop-klom);
                    Real dudz = mlebtensor_dz_on_xface(i,j,k,0,vel,dzi,
                                                       whi,wlo,khip,khim,klop,klom);
                    Real dvdz = mlebtensor_dz_on_xface(i,j,k,1,vel,dzi,
                                                       whi,wlo,khip,khim,klop,klom);
                    Real dwdz = mlebtensor_dz_on_xface(i,j,k,2,vel,dzi,
                                                       whi,wlo,khip,khim,klop,klom);

                    fx(i,j,k,0) = dudx;
                    fx(i,j,k,1) = dvdx;
                    fx(i,j,k,2) = dwdx;
                    fx(i,j,k,3) = dudy;
                    fx(i,j,k,4) = dvdy;
                    fx(i,j,k,5) = dwdy;
                    fx(i,j,k,6) = dudz;
                    fx(i,j,k,7) = dvdz;
                    fx(i,j,k,8) = dwdz;
                }
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebtensor_vel_grads_fy (Box const& box, Array4<Real> const& fy,
                              Array4<Real const> const& vel, Array4<Real const> const& apy,
                              Array4<EBCellFlag const> const& flag,
                              Array4<int const> const& ccm,
                              Array4<Real const> const& fcy,
                              GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    const Real dxi = dxinv[0];
    const Real dyi = dxinv[1];
    const Real dzi = dxinv[2];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                if (apy(i,j,k) == Real(0.0))
                {
                    fy(i,j,k,0) = Real(0.0);
                    fy(i,j,k,1) = Real(0.0);
                    fy(i,j,k,2) = Real(0.0);
                    fy(i,j,k,3) = Real(0.0);
                    fy(i,j,k,4) = Real(0.0);
                    fy(i,j,k,5) = Real(0.0);
                    fy(i,j,k,6) = Real(0.0);
                    fy(i,j,k,7) = Real(0.0);
                    fy(i,j,k,8) = Real(0.0);
                }
                else
                {
                    int ihip = i + flag(i,j  ,k).isConnected( 1,0,0);
                    int ihim = i - flag(i,j  ,k).isConnected(-1,0,0);
                    int ilop = i + flag(i,j-1,k).isConnected( 1,0,0);
                    int ilom = i - flag(i,j-1,k).isConnected(-1,0,0);
                    Real whi = mlebtensor_weight(ihip-ihim);
                    Real wlo = mlebtensor_weight(ilop-ilom);
                    Real dudx = mlebtensor_dx_on_yface(i,j,k,0,vel,dxi,
                                                       whi,wlo,ihip,ihim,ilop,ilom);
                    Real dvdx = mlebtensor_dx_on_yface(i,j,k,1,vel,dxi,
                                                       whi,wlo,ihip,ihim,ilop,ilom);
                    Real dwdx = mlebtensor_dx_on_yface(i,j,k,2,vel,dxi,
                                                       whi,wlo,ihip,ihim,ilop,ilom);

                    Real dudy = (vel(i,j,k,0) - vel(i,j-1,k,0))*dyi;
                    Real dvdy = (vel(i,j,k,1) - vel(i,j-1,k,1))*dyi;
                    Real dwdy = (vel(i,j,k,2) - vel(i,j-1,k,2))*dyi;
                    if (apy(i,j,k) < Real(1.0)) {
                        int ii = i + (fcy(i,j,k,0) >= Real(0.0) ? 1 : -1);
                        int kk = k + (fcy(i,j,k,1) >= Real(0.0) ? 1 : -1);
                        Real fracx = (ccm(ii,j-1,k) || ccm(ii,j,k)) ? std::abs(fcy(i,j,k,0)) : Real(0.0);
                        Real fracz = (ccm(i,j-1,kk) || ccm(i,j,kk)) ? std::abs(fcy(i,j,k,1)) : Real(0.0);
                        dudy = (Real(1.0)-fracx)*(Real(1.0)-fracz) * dudy
                            + fracx*(Real(1.0)-fracz) * (vel(ii,j,k ,0) - vel(ii,j-1,k ,0))*dyi
                            + fracz*(Real(1.0)-fracx) * (vel(i ,j,kk,0) - vel(i ,j-1,kk,0))*dyi
                            + fracx*     fracz  * (vel(ii,j,kk,0) - vel(ii,j-1,kk,0))*dyi;
                        dvdy = (Real(1.0)-fracx)*(Real(1.0)-fracz) * dvdy
                            + fracx*(Real(1.0)-fracz) * (vel(ii,j,k ,1) - vel(ii,j-1,k ,1))*dyi
                            + fracz*(Real(1.0)-fracx) * (vel(i ,j,kk,1) - vel(i ,j-1,kk,1))*dyi
                            + fracx*     fracz  * (vel(ii,j,kk,1) - vel(ii,j-1,kk,1))*dyi;
                        dwdy = (Real(1.0)-fracx)*(Real(1.0)-fracz) * dwdy
                            + fracx*(Real(1.0)-fracz) * (vel(ii,j,k ,2) - vel(ii,j-1,k ,2))*dyi
                            + fracz*(Real(1.0)-fracx) * (vel(i ,j,kk,2) - vel(i ,j-1,kk,2))*dyi
                            + fracx*     fracz  * (vel(ii,j,kk,2) - vel(ii,j-1,kk,2))*dyi;
                    }

                    int khip = k + flag(i,j  ,k).isConnected(0,0, 1);
                    int khim = k - flag(i,j  ,k).isConnected(0,0,-1);
                    int klop = k + flag(i,j-1,k).isConnected(0,0, 1);
                    int klom = k - flag(i,j-1,k).isConnected(0,0,-1);
                    whi = mlebtensor_weight(khip-khim);
                    wlo = mlebtensor_weight(klop-klom);
                    Real dudz = mlebtensor_dz_on_yface(i,j,k,0,vel,dzi,
                                                       whi,wlo,khip,khim,klop,klom);
                    Real dvdz = mlebtensor_dz_on_yface(i,j,k,1,vel,dzi,
                                                       whi,wlo,khip,khim,klop,klom);
                    Real dwdz = mlebtensor_dz_on_yface(i,j,k,2,vel,dzi,
                                                       whi,wlo,khip,khim,klop,klom);

                    fy(i,j,k,0) = dudx;
                    fy(i,j,k,1) = dvdx;
                    fy(i,j,k,2) = dwdx;
                    fy(i,j,k,3) = dudy;
                    fy(i,j,k,4) = dvdy;
                    fy(i,j,k,5) = dwdy;
                    fy(i,j,k,6) = dudz;
                    fy(i,j,k,7) = dvdz;
                    fy(i,j,k,8) = dwdz;
                }
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebtensor_vel_grads_fz (Box const& box, Array4<Real> const& fz,
                              Array4<Real const> const& vel, Array4<Real const> const& apz,
                              Array4<EBCellFlag const> const& flag,
                              Array4<int const> const& ccm,
                              Array4<Real const> const& fcz,
                              GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    const Real dxi = dxinv[0];
    const Real dyi = dxinv[1];
    const Real dzi = dxinv[2];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                if (apz(i,j,k) == Real(0.0))
                {
                    fz(i,j,k,0) = Real(0.0);
                    fz(i,j,k,1) = Real(0.0);
                    fz(i,j,k,2) = Real(0.0);
                    fz(i,j,k,3) = Real(0.0);
                    fz(i,j,k,4) = Real(0.0);
                    fz(i,j,k,5) = Real(0.0);
                    fz(i,j,k,6) = Real(0.0);
                    fz(i,j,k,7) = Real(0.0);
                    fz(i,j,k,8) = Real(0.0);
                }
                else
                {
                    int ihip = i + flag(i,j,k  ).isConnected( 1,0,0);
                    int ihim = i - flag(i,j,k  ).isConnected(-1,0,0);
                    int ilop = i + flag(i,j,k-1).isConnected( 1,0,0);
                    int ilom = i - flag(i,j,k-1).isConnected(-1,0,0);
                    Real whi = mlebtensor_weight(ihip-ihim);
                    Real wlo = mlebtensor_weight(ilop-ilom);
                    Real dudx = mlebtensor_dx_on_zface(i,j,k,0,vel,dxi,
                                                       whi,wlo,ihip,ihim,ilop,ilom);
                    Real dvdx = mlebtensor_dx_on_zface(i,j,k,1,vel,dxi,
                                                       whi,wlo,ihip,ihim,ilop,ilom);
                    Real dwdx = mlebtensor_dx_on_zface(i,j,k,2,vel,dxi,
                                                       whi,wlo,ihip,ihim,ilop,ilom);

                    int jhip = j + flag(i,j,k  ).isConnected(0, 1,0);
                    int jhim = j - flag(i,j,k  ).isConnected(0,-1,0);
                    int jlop = j + flag(i,j,k-1).isConnected(0, 1,0);
                    int jlom = j - flag(i,j,k-1).isConnected(0,-1,0);
                    whi = mlebtensor_weight(jhip-jhim);
                    wlo = mlebtensor_weight(jlop-jlom);
                    Real dudy = mlebtensor_dy_on_zface(i,j,k,0,vel,dyi,
                                                       whi,wlo,jhip,jhim,jlop,jlom);
                    Real dvdy = mlebtensor_dy_on_zface(i,j,k,1,vel,dyi,
                                                       whi,wlo,jhip,jhim,jlop,jlom);
                    Real dwdy = mlebtensor_dy_on_zface(i,j,k,2,vel,dyi,
                                                       whi,wlo,jhip,jhim,jlop,jlom);

                    Real dudz = (vel(i,j,k,0) - vel(i,j,k-1,0))*dzi;
                    Real dvdz = (vel(i,j,k,1) - vel(i,j,k-1,1))*dzi;
                    Real dwdz = (vel(i,j,k,2) - vel(i,j,k-1,2))*dzi;
                    if (apz(i,j,k) < Real(1.0)) {
                        int ii = i + (fcz(i,j,k,0) >= Real(0.0) ? 1 : -1);
                        int jj = j + (fcz(i,j,k,1) >= Real(0.0) ? 1 : -1);
                        Real fracx = (ccm(ii,j,k-1) || ccm(ii,j,k)) ? std::abs(fcz(i,j,k,0)) : Real(0.0);
                        Real fracy = (ccm(i,jj,k-1) || ccm(i,jj,k)) ? std::abs(fcz(i,j,k,1)) : Real(0.0);
                        dudz = (Real(1.0)-fracx)*(Real(1.0)-fracy) * dudz
                            + fracx*(Real(1.0)-fracy) * (vel(ii,j ,k,0) - vel(ii,j ,k-1,0))*dzi
                            + fracy*(Real(1.0)-fracx) * (vel(i ,jj,k,0) - vel(i ,jj,k-1,0))*dzi
                            + fracx*     fracy  * (vel(ii,jj,k,0) - vel(ii,jj,k-1,0))*dzi;
                        dvdz = (Real(1.0)-fracx)*(Real(1.0)-fracy) * dvdz
                            + fracx*(Real(1.0)-fracy) * (vel(ii,j ,k,1) - vel(ii,j ,k-1,1))*dzi
                            + fracy*(Real(1.0)-fracx) * (vel(i ,jj,k,1) - vel(i ,jj,k-1,1))*dzi
                            + fracx*     fracy  * (vel(ii,jj,k,1) - vel(ii,jj,k-1,1))*dzi;
                        dwdz = (Real(1.0)-fracx)*(Real(1.0)-fracy) * dwdz
                            + fracx*(Real(1.0)-fracy) * (vel(ii,j ,k,2) - vel(ii,j ,k-1,2))*dzi
                            + fracy*(Real(1.0)-fracx) * (vel(i ,jj,k,2) - vel(i ,jj,k-1,2))*dzi
                            + fracx*     fracy  * (vel(ii,jj,k,2) - vel(ii,jj,k-1,2))*dzi;
                    }

                    fz(i,j,k,0) = dudx;
                    fz(i,j,k,1) = dvdx;
                    fz(i,j,k,2) = dwdx;
                    fz(i,j,k,3) = dudy;
                    fz(i,j,k,4) = dvdy;
                    fz(i,j,k,5) = dwdy;
                    fz(i,j,k,6) = dudz;
                    fz(i,j,k,7) = dvdz;
                    fz(i,j,k,8) = dwdz;
                }
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebtensor_vel_grads_fx (Box const& box, Array4<Real> const& fx,
                              Array4<Real const> const& vel,
                              Array4<Real const> const& apx,
                              Array4<EBCellFlag const> const& flag,
                              Array4<int const> const& ccm,
                              Array4<Real const> const& fcx,
                              GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                              Array4<Real const> const& bvxlo,
                              Array4<Real const> const& bvxhi,
                              Array2D<BoundCond,
                                      0,2*AMREX_SPACEDIM,
                                      0,AMREX_SPACEDIM> const& bct,
                              Dim3 const& dlo, Dim3 const& dhi) noexcept
{
    const Real dxi = dxinv[0];
    const Real dyi = dxinv[1];
    const Real dzi = dxinv[2];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                if (apx(i,j,k) == Real(0.0))
                {
                    fx(i,j,k,0) = Real(0.0);
                    fx(i,j,k,1) = Real(0.0);
                    fx(i,j,k,2) = Real(0.0);
                    fx(i,j,k,3) = Real(0.0);
                    fx(i,j,k,4) = Real(0.0);
                    fx(i,j,k,5) = Real(0.0);
                    fx(i,j,k,6) = Real(0.0);
                    fx(i,j,k,7) = Real(0.0);
                    fx(i,j,k,8) = Real(0.0);
                }
                else
                {
                    Real dudx = (vel(i,j,k,0) - vel(i-1,j,k,0))*dxi;
                    Real dvdx = (vel(i,j,k,1) - vel(i-1,j,k,1))*dxi;
                    Real dwdx = (vel(i,j,k,2) - vel(i-1,j,k,2))*dxi;
                    if (apx(i,j,k) < Real(1.0)) {
                        int jj = j + (fcx(i,j,k,0) >= Real(0.0) ? 1 : -1);
                        int kk = k + (fcx(i,j,k,1) >= Real(0.0) ? 1 : -1);
                        Real fracy = (ccm(i-1,jj,k) || ccm(i,jj,k)) ? std::abs(fcx(i,j,k,0)) : Real(0.0);
                        Real fracz = (ccm(i-1,j,kk) || ccm(i,j,kk)) ? std::abs(fcx(i,j,k,1)) : Real(0.0);
                        dudx = (Real(1.0)-fracy)*(Real(1.0)-fracz) * dudx
                            + fracy*(Real(1.0)-fracz) * (vel(i,jj,k ,0) - vel(i-1,jj,k ,0))*dxi
                            + fracz*(Real(1.0)-fracy) * (vel(i,j ,kk,0) - vel(i-1,j ,kk,0))*dxi
                            + fracy*     fracz  * (vel(i,jj,kk,0) - vel(i-1,jj,kk,0))*dxi;
                        dvdx = (Real(1.0)-fracy)*(Real(1.0)-fracz) * dvdx
                            + fracy*(Real(1.0)-fracz) * (vel(i,jj,k ,1) - vel(i-1,jj,k ,1))*dxi
                            + fracz*(Real(1.0)-fracy) * (vel(i,j ,kk,1) - vel(i-1,j ,kk,1))*dxi
                            + fracy*     fracz  * (vel(i,jj,kk,1) - vel(i-1,jj,kk,1))*dxi;
                        dwdx = (Real(1.0)-fracy)*(Real(1.0)-fracz) * dwdx
                            + fracy*(Real(1.0)-fracz) * (vel(i,jj,k ,2) - vel(i-1,jj,k ,2))*dxi
                            + fracz*(Real(1.0)-fracy) * (vel(i,j ,kk,2) - vel(i-1,j ,kk,2))*dxi
                            + fracy*     fracz  * (vel(i,jj,kk,2) - vel(i-1,jj,kk,2))*dxi;
                    }

                    int jhip = j + flag(i  ,j,k).isConnected(0, 1,0);
                    int jhim = j - flag(i  ,j,k).isConnected(0,-1,0);
                    int jlop = j + flag(i-1,j,k).isConnected(0, 1,0);
                    int jlom = j - flag(i-1,j,k).isConnected(0,-1,0);
                    Real whi = mlebtensor_weight(jhip-jhim);
                    Real wlo = mlebtensor_weight(jlop-jlom);
                    Real dudy = mlebtensor_dy_on_xface(i,j,k,0,vel,dyi,
                                                       bvxlo,bvxhi,bct,dlo,dhi,
                                                       whi,wlo,jhip,jhim,jlop,jlom);
                    Real dvdy = mlebtensor_dy_on_xface(i,j,k,1,vel,dyi,
                                                       bvxlo,bvxhi,bct,dlo,dhi,
                                                       whi,wlo,jhip,jhim,jlop,jlom);
                    Real dwdy = mlebtensor_dy_on_xface(i,j,k,2,vel,dyi,
                                                       bvxlo,bvxhi,bct,dlo,dhi,
                                                       whi,wlo,jhip,jhim,jlop,jlom);

                    int khip = k + flag(i  ,j,k).isConnected(0,0, 1);
                    int khim = k - flag(i  ,j,k).isConnected(0,0,-1);
                    int klop = k + flag(i-1,j,k).isConnected(0,0, 1);
                    int klom = k - flag(i-1,j,k).isConnected(0,0,-1);
                    whi = mlebtensor_weight(khip-khim);
                    wlo = mlebtensor_weight(klop-klom);
                    Real dudz = mlebtensor_dz_on_xface(i,j,k,0,vel,dzi,
                                                       bvxlo,bvxhi,bct,dlo,dhi,
                                                       whi,wlo,khip,khim,klop,klom);
                    Real dvdz = mlebtensor_dz_on_xface(i,j,k,1,vel,dzi,
                                                       bvxlo,bvxhi,bct,dlo,dhi,
                                                       whi,wlo,khip,khim,klop,klom);
                    Real dwdz = mlebtensor_dz_on_xface(i,j,k,2,vel,dzi,
                                                       bvxlo,bvxhi,bct,dlo,dhi,
                                                       whi,wlo,khip,khim,klop,klom);

                    fx(i,j,k,0) = dudx;
                    fx(i,j,k,1) = dvdx;
                    fx(i,j,k,2) = dwdx;
                    fx(i,j,k,3) = dudy;
                    fx(i,j,k,4) = dvdy;
                    fx(i,j,k,5) = dwdy;
                    fx(i,j,k,6) = dudz;
                    fx(i,j,k,7) = dvdz;
                    fx(i,j,k,8) = dwdz;
                }
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebtensor_vel_grads_fy (Box const& box, Array4<Real> const& fy,
                              Array4<Real const> const& vel,
                              Array4<Real const> const& apy,
                              Array4<EBCellFlag const> const& flag,
                              Array4<int const> const& ccm,
                              Array4<Real const> const& fcy,
                              GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                              Array4<Real const> const& bvylo,
                              Array4<Real const> const& bvyhi,
                              Array2D<BoundCond,
                                      0,2*AMREX_SPACEDIM,
                                      0,AMREX_SPACEDIM> const& bct,
                              Dim3 const& dlo, Dim3 const& dhi) noexcept
{
    const Real dxi = dxinv[0];
    const Real dyi = dxinv[1];
    const Real dzi = dxinv[2];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                if (apy(i,j,k) == Real(0.0))
                {
                    fy(i,j,k,0) = Real(0.0);
                    fy(i,j,k,1) = Real(0.0);
                    fy(i,j,k,2) = Real(0.0);
                    fy(i,j,k,3) = Real(0.0);
                    fy(i,j,k,4) = Real(0.0);
                    fy(i,j,k,5) = Real(0.0);
                    fy(i,j,k,6) = Real(0.0);
                    fy(i,j,k,7) = Real(0.0);
                    fy(i,j,k,8) = Real(0.0);
                }
                else
                {
                    int ihip = i + flag(i,j  ,k).isConnected( 1,0,0);
                    int ihim = i - flag(i,j  ,k).isConnected(-1,0,0);
                    int ilop = i + flag(i,j-1,k).isConnected( 1,0,0);
                    int ilom = i - flag(i,j-1,k).isConnected(-1,0,0);
                    Real whi = mlebtensor_weight(ihip-ihim);
                    Real wlo = mlebtensor_weight(ilop-ilom);
                    Real dudx = mlebtensor_dx_on_yface(i,j,k,0,vel,dxi,
                                                       bvylo,bvyhi,bct,dlo,dhi,
                                                       whi,wlo,ihip,ihim,ilop,ilom);
                    Real dvdx = mlebtensor_dx_on_yface(i,j,k,1,vel,dxi,
                                                       bvylo,bvyhi,bct,dlo,dhi,
                                                       whi,wlo,ihip,ihim,ilop,ilom);
                    Real dwdx = mlebtensor_dx_on_yface(i,j,k,2,vel,dxi,
                                                       bvylo,bvyhi,bct,dlo,dhi,
                                                       whi,wlo,ihip,ihim,ilop,ilom);

                    Real dudy = (vel(i,j,k,0) - vel(i,j-1,k,0))*dyi;
                    Real dvdy = (vel(i,j,k,1) - vel(i,j-1,k,1))*dyi;
                    Real dwdy = (vel(i,j,k,2) - vel(i,j-1,k,2))*dyi;
                    if (apy(i,j,k) < Real(1.0)) {
                        int ii = i + (fcy(i,j,k,0) >= Real(0.0) ? 1 : -1);
                        int kk = k + (fcy(i,j,k,1) >= Real(0.0) ? 1 : -1);
                        Real fracx = (ccm(ii,j-1,k) || ccm(ii,j,k)) ? std::abs(fcy(i,j,k,0)) : Real(0.0);
                        Real fracz = (ccm(i,j-1,kk) || ccm(i,j,kk)) ? std::abs(fcy(i,j,k,1)) : Real(0.0);
                        dudy = (Real(1.0)-fracx)*(Real(1.0)-fracz) * dudy
                            + fracx*(Real(1.0)-fracz) * (vel(ii,j,k ,0) - vel(ii,j-1,k ,0))*dyi
                            + fracz*(Real(1.0)-fracx) * (vel(i ,j,kk,0) - vel(i ,j-1,kk,0))*dyi
                            + fracx*     fracz  * (vel(ii,j,kk,0) - vel(ii,j-1,kk,0))*dyi;
                        dvdy = (Real(1.0)-fracx)*(Real(1.0)-fracz) * dvdy
                            + fracx*(Real(1.0)-fracz) * (vel(ii,j,k ,1) - vel(ii,j-1,k ,1))*dyi
                            + fracz*(Real(1.0)-fracx) * (vel(i ,j,kk,1) - vel(i ,j-1,kk,1))*dyi
                            + fracx*     fracz  * (vel(ii,j,kk,1) - vel(ii,j-1,kk,1))*dyi;
                        dwdy = (Real(1.0)-fracx)*(Real(1.0)-fracz) * dwdy
                            + fracx*(Real(1.0)-fracz) * (vel(ii,j,k ,2) - vel(ii,j-1,k ,2))*dyi
                            + fracz*(Real(1.0)-fracx) * (vel(i ,j,kk,2) - vel(i ,j-1,kk,2))*dyi
                            + fracx*     fracz  * (vel(ii,j,kk,2) - vel(ii,j-1,kk,2))*dyi;
                    }

                    int khip = k + flag(i,j  ,k).isConnected(0,0, 1);
                    int khim = k - flag(i,j  ,k).isConnected(0,0,-1);
                    int klop = k + flag(i,j-1,k).isConnected(0,0, 1);
                    int klom = k - flag(i,j-1,k).isConnected(0,0,-1);
                    whi = mlebtensor_weight(khip-khim);
                    wlo = mlebtensor_weight(klop-klom);
                    Real dudz = mlebtensor_dz_on_yface(i,j,k,0,vel,dzi,
                                                       bvylo,bvyhi,bct,dlo,dhi,
                                                       whi,wlo,khip,khim,klop,klom);
                    Real dvdz = mlebtensor_dz_on_yface(i,j,k,1,vel,dzi,
                                                       bvylo,bvyhi,bct,dlo,dhi,
                                                       whi,wlo,khip,khim,klop,klom);
                    Real dwdz = mlebtensor_dz_on_yface(i,j,k,2,vel,dzi,
                                                       bvylo,bvyhi,bct,dlo,dhi,
                                                       whi,wlo,khip,khim,klop,klom);

                    fy(i,j,k,0) = dudx;
                    fy(i,j,k,1) = dvdx;
                    fy(i,j,k,2) = dwdx;
                    fy(i,j,k,3) = dudy;
                    fy(i,j,k,4) = dvdy;
                    fy(i,j,k,5) = dwdy;
                    fy(i,j,k,6) = dudz;
                    fy(i,j,k,7) = dvdz;
                    fy(i,j,k,8) = dwdz;
                }
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebtensor_vel_grads_fz (Box const& box, Array4<Real> const& fz,
                              Array4<Real const> const& vel,
                              Array4<Real const> const& apz,
                              Array4<EBCellFlag const> const& flag,
                              Array4<int const> const& ccm,
                              Array4<Real const> const& fcz,
                              GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                              Array4<Real const> const& bvzlo,
                              Array4<Real const> const& bvzhi,
                              Array2D<BoundCond,
                                      0,2*AMREX_SPACEDIM,
                                      0,AMREX_SPACEDIM> const& bct,
                              Dim3 const& dlo, Dim3 const& dhi) noexcept
{
    const Real dxi = dxinv[0];
    const Real dyi = dxinv[1];
    const Real dzi = dxinv[2];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                if (apz(i,j,k) == Real(0.0))
                {
                    fz(i,j,k,0) = Real(0.0);
                    fz(i,j,k,1) = Real(0.0);
                    fz(i,j,k,2) = Real(0.0);
                    fz(i,j,k,3) = Real(0.0);
                    fz(i,j,k,4) = Real(0.0);
                    fz(i,j,k,5) = Real(0.0);
                    fz(i,j,k,6) = Real(0.0);
                    fz(i,j,k,7) = Real(0.0);
                    fz(i,j,k,8) = Real(0.0);
                }
                else
                {
                    int ihip = i + flag(i,j,k  ).isConnected( 1,0,0);
                    int ihim = i - flag(i,j,k  ).isConnected(-1,0,0);
                    int ilop = i + flag(i,j,k-1).isConnected( 1,0,0);
                    int ilom = i - flag(i,j,k-1).isConnected(-1,0,0);
                    Real whi = mlebtensor_weight(ihip-ihim);
                    Real wlo = mlebtensor_weight(ilop-ilom);
                    Real dudx = mlebtensor_dx_on_zface(i,j,k,0,vel,dxi,
                                                       bvzlo,bvzhi,bct,dlo,dhi,
                                                       whi,wlo,ihip,ihim,ilop,ilom);
                    Real dvdx = mlebtensor_dx_on_zface(i,j,k,1,vel,dxi,
                                                       bvzlo,bvzhi,bct,dlo,dhi,
                                                       whi,wlo,ihip,ihim,ilop,ilom);
                    Real dwdx = mlebtensor_dx_on_zface(i,j,k,2,vel,dxi,
                                                       bvzlo,bvzhi,bct,dlo,dhi,
                                                       whi,wlo,ihip,ihim,ilop,ilom);

                    int jhip = j + flag(i,j,k  ).isConnected(0, 1,0);
                    int jhim = j - flag(i,j,k  ).isConnected(0,-1,0);
                    int jlop = j + flag(i,j,k-1).isConnected(0, 1,0);
                    int jlom = j - flag(i,j,k-1).isConnected(0,-1,0);
                    whi = mlebtensor_weight(jhip-jhim);
                    wlo = mlebtensor_weight(jlop-jlom);
                    Real dudy = mlebtensor_dy_on_zface(i,j,k,0,vel,dyi,
                                                       bvzlo,bvzhi,bct,dlo,dhi,
                                                       whi,wlo,jhip,jhim,jlop,jlom);
                    Real dvdy = mlebtensor_dy_on_zface(i,j,k,1,vel,dyi,
                                                       bvzlo,bvzhi,bct,dlo,dhi,
                                                       whi,wlo,jhip,jhim,jlop,jlom);
                    Real dwdy = mlebtensor_dy_on_zface(i,j,k,2,vel,dyi,
                                                       bvzlo,bvzhi,bct,dlo,dhi,
                                                       whi,wlo,jhip,jhim,jlop,jlom);

                    Real dudz = (vel(i,j,k,0) - vel(i,j,k-1,0))*dzi;
                    Real dvdz = (vel(i,j,k,1) - vel(i,j,k-1,1))*dzi;
                    Real dwdz = (vel(i,j,k,2) - vel(i,j,k-1,2))*dzi;
                    if (apz(i,j,k) < Real(1.0)) {
                        int ii = i + (fcz(i,j,k,0) >= Real(0.0) ? 1 : -1);
                        int jj = j + (fcz(i,j,k,1) >= Real(0.0) ? 1 : -1);
                        Real fracx = (ccm(ii,j,k-1) || ccm(ii,j,k)) ? std::abs(fcz(i,j,k,0)) : Real(0.0);
                        Real fracy = (ccm(i,jj,k-1) || ccm(i,jj,k)) ? std::abs(fcz(i,j,k,1)) : Real(0.0);
                        dudz = (Real(1.0)-fracx)*(Real(1.0)-fracy) * dudz
                            + fracx*(Real(1.0)-fracy) * (vel(ii,j ,k,0) - vel(ii,j ,k-1,0))*dzi
                            + fracy*(Real(1.0)-fracx) * (vel(i ,jj,k,0) - vel(i ,jj,k-1,0))*dzi
                            + fracx*     fracy  * (vel(ii,jj,k,0) - vel(ii,jj,k-1,0))*dzi;
                        dvdz = (Real(1.0)-fracx)*(Real(1.0)-fracy) * dvdz
                            + fracx*(Real(1.0)-fracy) * (vel(ii,j ,k,1) - vel(ii,j ,k-1,1))*dzi
                            + fracy*(Real(1.0)-fracx) * (vel(i ,jj,k,1) - vel(i ,jj,k-1,1))*dzi
                            + fracx*     fracy  * (vel(ii,jj,k,1) - vel(ii,jj,k-1,1))*dzi;
                        dwdz = (Real(1.0)-fracx)*(Real(1.0)-fracy) * dwdz
                            + fracx*(Real(1.0)-fracy) * (vel(ii,j ,k,2) - vel(ii,j ,k-1,2))*dzi
                            + fracy*(Real(1.0)-fracx) * (vel(i ,jj,k,2) - vel(i ,jj,k-1,2))*dzi
                            + fracx*     fracy  * (vel(ii,jj,k,2) - vel(ii,jj,k-1,2))*dzi;
                    }

                    fz(i,j,k,0) = dudx;
                    fz(i,j,k,1) = dvdx;
                    fz(i,j,k,2) = dwdx;
                    fz(i,j,k,3) = dudy;
                    fz(i,j,k,4) = dvdy;
                    fz(i,j,k,5) = dwdy;
                    fz(i,j,k,6) = dudz;
                    fz(i,j,k,7) = dvdz;
                    fz(i,j,k,8) = dwdz;
                }
            }
        }
    }
}

}

#endif
