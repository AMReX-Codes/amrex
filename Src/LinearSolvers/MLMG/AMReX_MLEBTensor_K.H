#ifndef AMREX_MLEBTENSOR_K_H_
#define AMREX_MLEBTENSOR_K_H_
#include <AMReX_Config.H>

#include <AMReX_FArrayBox.H>

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mlebtensor_weight (int d) {
    return (d==2) ? 0.5 : ((d==1) ? 1.0 : 0.0);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mlebtensor_dy_on_xface (int i, int, int k, int n,
                             Array4<Real const> const& vel, Real dyi,
                             Real whi, Real wlo,
                             int jhip, int jhim, int jlop, int jlom) noexcept
{
    return Real(0.5)*dyi * ((vel(i  ,jhip,k,n)-vel(i  ,jhim,k,n))*whi +
                            (vel(i-1,jlop,k,n)-vel(i-1,jlom,k,n))*wlo);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mlebtensor_dx_on_yface (int, int j, int k, int n,
                             Array4<Real const> const& vel, Real dxi,
                             Real whi, Real wlo,
                             int ihip, int ihim, int ilop, int ilom) noexcept
{
    return Real(0.5)*dxi * ((vel(ihip,j  ,k,n)-vel(ihim,j  ,k,n))*whi +
                            (vel(ilop,j-1,k,n)-vel(ilom,j-1,k,n))*wlo);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mlebtensor_dy_on_xface (int i, int j, int k, int n,
                             Array4<Real const> const& vel, Real dyi,
                             Array4<Real const> const& bvxlo,
                             Array4<Real const> const& bvxhi,
                             Array2D<BoundCond,
                                     0,2*AMREX_SPACEDIM,
                                     0,AMREX_SPACEDIM> const& bct,
                             Dim3 const& dlo, Dim3 const& dhi,
                             Real whi, Real wlo,
                             int jhip, int jhim, int jlop, int jlom) noexcept
{
    Real ddy;
    if (i == dlo.x) {
        if (bct(Orientation::xlo(),n) == AMREX_LO_DIRICHLET && bvxlo) {
            if (j == dlo.y) {
                ddy = (bvxlo(i-1,j  ,k,n) * Real(-1.5) +
                       bvxlo(i-1,j+1,k,n) * Real(2.) +
                       bvxlo(i-1,j+2,k,n) * Real(-0.5)) * dyi;
            } else if (j == dhi.y) {
                ddy = -(bvxlo(i-1,j  ,k,n) * Real(-1.5) +
                        bvxlo(i-1,j-1,k,n) * Real(2.) +
                        bvxlo(i-1,j-2,k,n) * Real(-0.5)) * dyi;
            } else {
                ddy = wlo*dyi*(bvxlo(i-1,jlop,k,n)-bvxlo(i-1,jlom,k,n));
            }
        } else if (bct(Orientation::xlo(),n) == AMREX_LO_NEUMANN) {
            ddy = whi*dyi*(vel(i,jhip,k,n)-vel(i,jhim,k,n));
        } else { // AMREX_LO_REFLECT_ODD or homogeneous Dirichlet
            ddy = Real(0.);
        }
    } else if (i == dhi.x+1) {
        if (bct(Orientation::xhi(),n) == AMREX_LO_DIRICHLET && bvxhi) {
            if (j == dlo.y) {
                ddy = (bvxhi(i,j  ,k,n) * Real(-1.5) +
                       bvxhi(i,j+1,k,n) * Real(2.) +
                       bvxhi(i,j+2,k,n) * Real(-0.5)) * dyi;
            } else if (j == dhi.y) {
                ddy = -(bvxhi(i,j  ,k,n) * Real(-1.5) +
                        bvxhi(i,j-1,k,n) * Real(2.) +
                        bvxhi(i,j-2,k,n) * Real(-0.5)) * dyi;
            } else {
                ddy = whi*dyi*(bvxhi(i,jhip,k,n)-bvxhi(i,jhim,k,n));
            }
        } else if (bct(Orientation::xhi(),n) == AMREX_LO_NEUMANN) {
            ddy = wlo*dyi*(vel(i-1,jlop,k,n)-vel(i-1,jlom,k,n));
        } else { // AMREX_LO_REFLECT_ODD or homogeneous Dirichlet
            ddy = Real(0.);
        }
    } else {
        ddy = mlebtensor_dy_on_xface(i,j,k,n,vel,dyi,whi,wlo,jhip,jhim,jlop,jlom);
    }
    return ddy;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mlebtensor_dx_on_yface (int i, int j, int k, int n,
                             Array4<Real const> const& vel, Real dxi,
                             Array4<Real const> const& bvylo,
                             Array4<Real const> const& bvyhi,
                             Array2D<BoundCond,
                                     0,2*AMREX_SPACEDIM,
                                     0,AMREX_SPACEDIM> const& bct,
                             Dim3 const& dlo, Dim3 const& dhi,
                             Real whi, Real wlo,
                             int ihip, int ihim, int ilop, int ilom) noexcept
{
    Real ddx;
    if (j == dlo.y) {
        if (bct(Orientation::ylo(),n) == AMREX_LO_DIRICHLET && bvylo) {
            if (i == dlo.x) {
                ddx = (bvylo(i  ,j-1,k,n) * Real(-1.5) +
                       bvylo(i+1,j-1,k,n) * Real(2.) +
                       bvylo(i+2,j-1,k,n) * Real(-0.5)) * dxi;
            } else if (i == dhi.x) {
                ddx = -(bvylo(i  ,j-1,k,n) * Real(-1.5) +
                        bvylo(i-1,j-1,k,n) * Real(2.) +
                        bvylo(i-2,j-1,k,n) * Real(-0.5)) * dxi;
            } else {
                ddx = wlo*dxi*(bvylo(ilop,j-1,k,n)-bvylo(ilom,j-1,k,n));
            }
        } else if (bct(Orientation::ylo(),n) == AMREX_LO_NEUMANN) {
            ddx = whi*dxi*(vel(ihip,j,k,n)-vel(ihim,j,k,n));
        } else { // AMREX_LO_REFLECT_ODD or homogeneous Dirichlet
            ddx = Real(0.);
        }
    } else if (j == dhi.y+1) {
        if (bct(Orientation::yhi(),n) == AMREX_LO_DIRICHLET && bvyhi) {
            if (i == dlo.x) {
                ddx = (bvyhi(i  ,j,k,n) * Real(-1.5) +
                       bvyhi(i+1,j,k,n) * Real(2.) +
                       bvyhi(i+2,j,k,n) * Real(-0.5)) * dxi;
            } else if (i == dhi.x) {
                ddx = -(bvyhi(i  ,j,k,n) * Real(-1.5) +
                        bvyhi(i-1,j,k,n) * Real(2.) +
                        bvyhi(i-2,j,k,n) * Real(-0.5)) * dxi;
            } else {
                ddx = whi*dxi*(bvyhi(ihip,j,k,n)-bvyhi(ihim,j,k,n));
            }
        } else if (bct(Orientation::yhi(),n) == AMREX_LO_NEUMANN) {
            ddx = wlo*dxi*(vel(ilop,j-1,k,n)-vel(ilom,j-1,k,n));
        } else { // AMREX_LO_REFLECT_ODD or homogeneous Dirichlet
            ddx = Real(0.);
        }
    } else {
        ddx = mlebtensor_dx_on_yface(i,j,k,n,vel,dxi,whi,wlo,ihip,ihim,ilop,ilom);
    }
    return ddx;
}

}

#if (AMREX_SPACEDIM == 1)
#elif (AMREX_SPACEDIM == 2)
#include <AMReX_MLEBTensor_2D_K.H>
#else
#include <AMReX_MLEBTensor_3D_K.H>
#endif

#endif
