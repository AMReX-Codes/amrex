#ifndef EM_PIC_K_H_
#define EM_PIC_K_H_

#include "Constants.H"
#include "EMParticleContainer.H"
#include <AMReX_FArrayBox.H>
#include <cmath>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void gather_fields (EMParticleContainer::ParticleType const& p,
                    amrex::Real& Exp, amrex::Real& Eyp, amrex::Real& Ezp,
                    amrex::Real& Bxp, amrex::Real& Byp, amrex::Real& Bzp,
                    amrex::Array4<amrex::Real const> const& Exarr,
                    amrex::Array4<amrex::Real const> const& Eyarr,
                    amrex::Array4<amrex::Real const> const& Ezarr,
                    amrex::Array4<amrex::Real const> const& Bxarr,
                    amrex::Array4<amrex::Real const> const& Byarr,
                    amrex::Array4<amrex::Real const> const& Bzarr,
                    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi)
{
    amrex::Real x = (p.pos(0) - plo[0]) * dxi[0];
    amrex::Real y = (p.pos(1) - plo[1]) * dxi[1];
    amrex::Real z = (p.pos(2) - plo[2]) * dxi[2];

    int j = std::floor(x);
    int k = std::floor(y);
    int l = std::floor(z);
    int j0 = j;
    int k0 = k;
    int l0 = l;

    amrex::Real xint = x - j;
    amrex::Real yint = y - k;
    amrex::Real zint = z - l;

    amrex::Real sx[] = {1.-xint, xint};
    amrex::Real sy[] = {1.-yint, yint};
    amrex::Real sz[] = {1.-zint, zint};

    xint = x - 0.5 - j0;
    yint = y - 0.5 - k0;
    zint = z - 0.5 - l0;

    constexpr amrex::Real sx0[] = {1.,0.};
    constexpr amrex::Real sy0[] = {1.,0.};
    constexpr amrex::Real sz0[] = {1.,0.};
    
    constexpr int ixmin = 0;
    constexpr int ixmax = 0;
    constexpr int iymin = 0;
    constexpr int iymax = 0;
    constexpr int izmin = 0;
    constexpr int izmax = 0;

    constexpr int ixmin0 = 0;
    constexpr int ixmax0 = 0;
    constexpr int iymin0 = 0;
    constexpr int iymax0 = 0;
    constexpr int izmin0 = 0;
    constexpr int izmax0 = 0;

    Exp = 0.;
    for         (int ll = izmin ; ll <= izmax +1; ++ll) {
        for     (int kk = iymin ; kk <= iymax +1; ++kk) {
            for (int jj = ixmin0; jj <= ixmax0  ; ++jj) {
                Exp += sx0[jj]*sy[kk]*sz[ll]*Exarr(j0+jj,k+kk,l+ll);
            }
        }
    }            

    Eyp = 0.;
    for         (int ll = izmin ; ll <= izmax +1; ++ll) {
        for     (int kk = iymin0; kk <= iymax0  ; ++kk) {
            for (int jj = ixmin ; jj <= ixmax +1; ++jj) {
                Eyp += sx[jj]*sy0[kk]*sz[ll]*Eyarr(j+jj,k0+kk,l+ll);
            }
        }
    }            

    Ezp = 0.;
    for         (int ll = izmin0; ll <= izmax0  ; ++ll) {
        for     (int kk = iymin ; kk <= iymax +1; ++kk) {
            for (int jj = ixmin ; jj <= ixmax +1; ++jj) {
                Ezp += sx[jj]*sy[kk]*sz0[ll]*Ezarr(j+jj,k+kk,l0+ll);
            }
        }
    }            

    Bxp = 0.;
    for         (int ll = izmin0 ; ll <= izmax0 ; ++ll) {
        for     (int kk = iymin0 ; kk <= iymax0 ; ++kk) {
            for (int jj = ixmin  ; jj <= ixmax+1; ++jj) {
                Bxp += sx[jj]*sy0[kk]*sz0[ll]*Bxarr(j+jj,k0+kk,l0+ll);
            }
        }
    }            

    Byp = 0.;
    for         (int ll = izmin0; ll <= izmax0  ; ++ll) {
        for     (int kk = iymin ; kk <= iymax +1; ++kk) {
            for (int jj = ixmin0; jj <= ixmax0  ; ++jj) {
                Byp += sx0[jj]*sy[kk]*sz0[ll]*Byarr(j0+jj,k+kk,l0+ll);
            }
        }
    }            

    Bzp = 0.;
    for         (int ll = izmin ; ll <= izmax +1; ++ll) {
        for     (int kk = iymin0; kk <= iymax0  ; ++kk) {
            for (int jj = ixmin0; jj <= ixmax0  ; ++jj) {
                Bzp += sx0[jj]*sy0[kk]*sz[ll]*Bzarr(j0+jj,k0+kk,l+ll);
            }
        }
    }                
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void push_momentum_boris(amrex::Real& uxp, amrex::Real& uyp, amrex::Real& uzp, amrex::Real& gaminv,
                         amrex::Real Exp, amrex::Real Eyp, amrex::Real Ezp,
                         amrex::Real Bxp, amrex::Real Byp, amrex::Real Bzp,
                         amrex::Real q, amrex::Real m, amrex::Real dt)
{
    amrex::Real cst = q*dt*0.5/m;
    constexpr amrex::Real clight = 2.99792458e8;
    constexpr amrex::Real clightisq = 1./(clight*clight);
    
    uxp += Exp * cst;
    uyp += Eyp * cst;
    uzp += Ezp * cst;

    amrex::Real usq = (uxp*uxp + uyp*uyp + uzp*uzp) * clightisq;
    amrex::Real gaminvtmp = 1.0/std::sqrt(1.0+usq);

    amrex::Real tx = gaminvtmp * Bxp * cst;
    amrex::Real ty = gaminvtmp * Byp * cst;
    amrex::Real tz = gaminvtmp * Bzp * cst;
    amrex::Real tsqi = 2.0/(1.0 + tx*tx + ty*ty + tz*tz);
    amrex::Real sx = tx*tsqi;
    amrex::Real sy = ty*tsqi;
    amrex::Real sz = tz*tsqi;
    amrex::Real uxppr = uxp + uyp*tz - uzp*ty;
    amrex::Real uyppr = uyp + uzp*tx - uxp*tz;
    amrex::Real uzppr = uzp + uxp*ty - uyp*tx;
    uxp += uyppr*sz - uzppr*sy;
    uyp += uzppr*sx - uxppr*sz;
    uzp += uxppr*sy - uyppr*sx;

    uxp += Exp*cst;
    uyp += Eyp*cst;
    uzp += Ezp*cst;

    usq = (uxp*uxp + uyp*uyp + uzp*uzp) * clightisq;
    gaminv = 1.0/std::sqrt(1.0+usq);    
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void push_position_boris(EMParticleContainer::ParticleType& p,
                         amrex::Real uxp, amrex::Real uyp, amrex::Real uzp,
                         amrex::Real gaminv, amrex::Real dt)
{
    p.pos(0) += uxp*gaminv*dt;
    p.pos(1) += uyp*gaminv*dt;
    p.pos(2) += uzp*gaminv*dt;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void deposit_current (amrex::Array4<amrex::Real> const& jx,
                      amrex::Array4<amrex::Real> const& jy,
                      amrex::Array4<amrex::Real> const& jz,
                      EMParticleContainer::ParticleType const& p,
                      amrex::Real uxp, amrex::Real uyp, amrex::Real uzp,
                      amrex::Real gaminv, amrex::Real w, amrex::Real q, amrex::Real dt,
                      amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                      amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi)
{
    amrex::Real invvol = dxi[0]*dxi[1]*dxi[2];
    amrex::Real dts2dx = 0.5*dt*dxi[0];
    amrex::Real dts2dy = 0.5*dt*dxi[1];
    amrex::Real dts2dz = 0.5*dt*dxi[2];

    amrex::Real x = (p.pos(0) - plo[0]) * dxi[0];
    amrex::Real y = (p.pos(1) - plo[1]) * dxi[1];
    amrex::Real z = (p.pos(2) - plo[2]) * dxi[2];

    amrex::Real vx = uxp*gaminv;
    amrex::Real vy = uyp*gaminv;
    amrex::Real vz = uzp*gaminv;

    amrex::Real wq=q*w;
    amrex::Real wqx=wq*invvol*vx;
    amrex::Real wqy=wq*invvol*vy;
    amrex::Real wqz=wq*invvol*vz;

    amrex::Real xmid=x-dts2dx*vx;
    amrex::Real ymid=y-dts2dy*vy;
    amrex::Real zmid=z-dts2dz*vz;

    int j  = std::floor(xmid);
    int k  = std::floor(ymid);
    int l  = std::floor(zmid);
    int j0 = std::floor(xmid-0.5);
    int k0 = std::floor(ymid-0.5);
    int l0 = std::floor(zmid-0.5);

    amrex::Real xint = xmid-j;
    amrex::Real yint = ymid-k;
    amrex::Real zint = zmid-l;
    amrex::Real sx[] = {1.-xint, xint};
    amrex::Real sy[] = {1.-yint, yint};
    amrex::Real sz[] = {1.-zint, zint};

    xint = xmid-j0-0.5;
    yint = ymid-k0-0.5;
    zint = zmid-l0-0.5;
    amrex::Real sx0[] = {1.-xint, xint};
    amrex::Real sy0[] = {1.-yint, yint};
    amrex::Real sz0[] = {1.-zint, zint};

    for         (int loff = 0; loff < 2; ++loff) {
        for     (int koff = 0; koff < 2; ++koff) {
            for (int joff = 0; joff < 2; ++joff) {
                amrex::Gpu::Atomic::Add(&jx(j0+joff,k +koff,l +loff), sx0[joff]*sy [koff]*sz [loff]*wqx);
                amrex::Gpu::Atomic::Add(&jy(j +joff,k0+koff,l +loff), sx [joff]*sy0[koff]*sz [loff]*wqy);
                amrex::Gpu::Atomic::Add(&jz(j +joff,k +koff,l0+loff), sx [joff]*sy [koff]*sz0[loff]*wqz);
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void push_electric_field_x (int j, int k, int l,
                            amrex::Array4<amrex::Real      > const& Ex,
                            amrex::Array4<amrex::Real const> const& By,
                            amrex::Array4<amrex::Real const> const& Bz,
                            amrex::Array4<amrex::Real const> const& jx,
                            amrex::Real mudt, amrex::Real dtsdy, amrex::Real dtsdz)
{
    Ex(j,k,l) = Ex(j,k,l) + dtsdy * (Bz(j,k,l) - Bz(j,k-1,l  ))
                          - dtsdz * (By(j,k,l) - By(j,k  ,l-1))
        - mudt  * jx(j,k,l);
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void push_electric_field_y (int j, int k, int l,
                            amrex::Array4<amrex::Real      > const& Ey,
                            amrex::Array4<amrex::Real const> const& Bx,
                            amrex::Array4<amrex::Real const> const& Bz,
                            amrex::Array4<amrex::Real const> const& jy,
                            amrex::Real mudt, amrex::Real dtsdx, amrex::Real dtsdz)
{
    Ey(j,k,l) = Ey(j,k,l) - dtsdx * (Bz(j,k,l) - Bz(j-1,k,l))
                          + dtsdz * (Bx(j,k,l) - Bx(j,k,l-1))
        - mudt  * jy(j,k,l);
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void push_electric_field_z (int j, int k, int l,
                            amrex::Array4<amrex::Real      > const& Ez,
                            amrex::Array4<amrex::Real const> const& Bx,
                            amrex::Array4<amrex::Real const> const& By,
                            amrex::Array4<amrex::Real const> const& jz,
                            amrex::Real mudt, amrex::Real dtsdx, amrex::Real dtsdy)
{
    Ez(j,k,l) = Ez(j,k,l) + dtsdx * (By(j,k,l) - By(j-1,k  ,l))
                          - dtsdy * (Bx(j,k,l) - Bx(j  ,k-1,l))
        - mudt  * jz(j,k,l);
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void push_magnetic_field_x (int j, int k, int l,
                            amrex::Array4<amrex::Real      > const& Bx,
                            amrex::Array4<amrex::Real const> const& Ey,
                            amrex::Array4<amrex::Real const> const& Ez,
                            amrex::Real dtsdy, amrex::Real dtsdz)
{
    Bx(j,k,l) = Bx(j,k,l) - dtsdy * (Ez(j  ,k+1,l  ) - Ez(j,k,l))
                          + dtsdz * (Ey(j  ,k  ,l+1) - Ey(j,k,l));
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void push_magnetic_field_y (int j, int k, int l,
                            amrex::Array4<amrex::Real      > const& By,
                            amrex::Array4<amrex::Real const> const& Ex,
                            amrex::Array4<amrex::Real const> const& Ez,
                            amrex::Real dtsdx, amrex::Real dtsdz)
{
    By(j,k,l) = By(j,k,l) + dtsdx * (Ez(j+1,k  ,l  ) - Ez(j,k,l))
                          - dtsdz * (Ex(j  ,k  ,l+1) - Ex(j,k,l));
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void push_magnetic_field_z (int j, int k, int l,
                            amrex::Array4<amrex::Real      > const& Bz,
                            amrex::Array4<amrex::Real const> const& Ex,
                            amrex::Array4<amrex::Real const> const& Ey,
                            amrex::Real dtsdx, amrex::Real dtsdy)
{
    Bz(j,k,l) = Bz(j,k,l) - dtsdx * (Ey(j+1,k  ,l  ) - Ey(j,k,l))
                          + dtsdy * (Ex(j  ,k+1,l  ) - Ex(j,k,l));
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real
check_langmuir_solution (amrex::Box const& bx, amrex::Box const& test_box,
                         amrex::FArrayBox const& jxfab, amrex::Real j_exact)
{
    amrex::Real error = 0.0;

    const amrex::Box b = bx & test_box;
    const auto lo  = amrex::lbound(b);
    const auto len = amrex::length(b);
    const auto jx = jxfab.view(lo);

    for         (int l = 0; l < len.z; ++l) {
        for     (int k = 0; k < len.y; ++k) {
            for (int j = 0; j < len.x; ++j) {
                error = amrex::max(error, std::abs((jx(j,k,l) - j_exact) / j_exact));
            }
        }
    }

    return error;
}

#endif
