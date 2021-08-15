#ifndef AMREX_SLOPES_K_H_
#define AMREX_SLOPES_K_H_
#include <AMReX_Config.H>

#include <AMReX_Array4.H>
#include <AMReX_Extension.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Gpu.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>


namespace amrex {

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
Real amrex_calc_xslope (int i, int j, int k, int n, int order,
                               amrex::Array4<Real const> const& q) noexcept
{
    if (order == 2)
    {
        Real dl = 2.0_rt*(q(i  ,j,k,n) - q(i-1,j,k,n));
        Real dr = 2.0_rt*(q(i+1,j,k,n) - q(i  ,j,k,n));
        Real dc = 0.5_rt*(q(i+1,j,k,n) - q(i-1,j,k,n));
        Real slope = amrex::min(Math::abs(dl),Math::abs(dc),Math::abs(dr));
        slope = (dr*dl > 0.0_rt) ? slope : 0.0_rt;
        return (dc > 0.0_rt) ? slope : -slope;

    } else if (order == 4) {

        Real dlft, drgt, dcen, dfm, dfp, dlim, dsgn, dtemp;
        Real qm, qp, qi;
        qi = q(i,j,k,n);
        qm = q(i-1,j,k,n);
        qp = q(i+1,j,k,n);

        dlft = qm - q(i-2,j,k,n);
        drgt = qi - qm;
        dcen = 0.5_rt*(dlft+drgt);
        dsgn = Math::copysign(1.e0_rt, dcen);
        dlim = (dlft*drgt >= 0.0_rt) ? 2.0_rt*amrex::min(Math::abs(dlft), Math::abs(drgt)) : 0.0_rt;
        dfm = dsgn*amrex::min(dlim, Math::abs(dcen));

        dlft = qp - qi;
        drgt = q(i+2,j,k,n) - qp;
        dcen = 0.5_rt*(dlft+drgt);
        dsgn = Math::copysign(1.e0_rt, dcen);
        dlim = (dlft*drgt >= 0.0_rt) ? 2.0_rt*amrex::min(Math::abs(dlft), Math::abs(drgt)) : 0.0_rt;
        dfp = dsgn*amrex::min(dlim, Math::abs(dcen));

        dlft = qi - qm;
        drgt = qp - qi;
        dcen = 0.5_rt*(dlft+drgt);
        dsgn = Math::copysign(1.e0_rt, dcen);
        dlim = (dlft*drgt >= 0.0_rt) ? 2.0_rt*amrex::min(Math::abs(dlft), Math::abs(drgt)) : 0.0_rt;

        dtemp  = 4.0_rt/3.0_rt*dcen - 1.0_rt/6.0_rt*(dfp + dfm);

        return dsgn*amrex::min(dlim, Math::abs(dtemp));

    } else {
        return 0.;
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
Real amrex_calc_xslope_extdir (int i, int j, int k, int n, int order,
                                      amrex::Array4<Real const> const& q,
                                      bool edlo, bool edhi, int domlo, int domhi) noexcept
{
    if (order == 2)
    {
        Real dl = 2.0_rt*(q(i  ,j,k,n) - q(i-1,j,k,n));
        Real dr = 2.0_rt*(q(i+1,j,k,n) - q(i  ,j,k,n));
        Real dc = 0.5_rt*(q(i+1,j,k,n) - q(i-1,j,k,n));

        if (edlo && i == domlo) {
            dc = (q(i+1,j,k,n)+3.0_rt*q(i,j,k,n)-4.0_rt*q(i-1,j,k,n))/3.0_rt;
        } else if (edhi && i == domhi) {
            dc = (4.0_rt*q(i+1,j,k,n)-3.0_rt*q(i,j,k,n)-q(i-1,j,k,n))/3.0_rt;
        }

        Real slope = amrex::min(Math::abs(dl),Math::abs(dc),Math::abs(dr));
        slope = (dr*dl > 0.0_rt) ? slope : 0.0_rt;
        return (dc > 0.0_rt) ? slope : -slope;

    } else if (order == 4) {

        Real dlft, drgt, dcen, dfm, dfp, dlim, dsgn, dtemp, dlimsh, dsgnsh;
        Real qm, qp, qi;
        qi = q(i,j,k,n);
        qm = q(i-1,j,k,n);
        qp = q(i+1,j,k,n);

        dlft = qm - q(i-2,j,k,n);
        drgt = qi - qm;
        dcen = 0.5_rt*(dlft+drgt);
        dsgn = Math::copysign(1.e0_rt, dcen);
        dlim = (dlft*drgt >= 0.0_rt) ? 2.0_rt*amrex::min(Math::abs(dlft), Math::abs(drgt)) : 0.0_rt;
        dfm = dsgn*amrex::min(dlim, Math::abs(dcen));

        dlft = qp - qi;
        drgt = q(i+2,j,k,n) - qp;
        dcen = 0.5_rt*(dlft+drgt);
        dsgn = Math::copysign(1.e0_rt, dcen);
        dlim = (dlft*drgt >= 0.0_rt) ? 2.0_rt*amrex::min(Math::abs(dlft), Math::abs(drgt)) : 0.0_rt;
        dfp = dsgn*amrex::min(dlim, Math::abs(dcen));

        dlft = qi - qm;
        drgt = qp - qi;
        dcen = 0.5_rt*(dlft+drgt);
        dsgn = Math::copysign(1.e0_rt, dcen);
        dlim = (dlft*drgt >= 0.0_rt) ? 2.0_rt*amrex::min(Math::abs(dlft), Math::abs(drgt)) : 0.0_rt;

        dtemp  = 4.0_rt/3.0_rt*dcen - 1.0_rt/6.0_rt*(dfp + dfm);

        if (edlo && i == domlo) {
           dtemp  = -16._rt/15._rt*q(i-1,j,k,n) + .5_rt*q(i,j,k,n) + 2._rt/3._rt*q(i+1,j,k,n) -  0.1_rt*q(i+2,j,k,n);
           dlft = 2._rt*(q(i  ,j,k,n)-q(i-1,j,k,n));
           drgt = 2._rt*(q(i+1,j,k,n)-q(i  ,j,k,n));
           dlim = (dlft*drgt >= 0.0_rt) ? amrex::min(Math::abs(dlft), Math::abs(drgt)) : 0.0_rt;
           dsgn = Math::copysign(1.e0_rt, dtemp);
        } else if (edlo && i == domlo+1) {
           dfm  = -16._rt/15._rt*q(domlo-1,j,k,n) + .5_rt*q(domlo,j,k,n) + 2._rt/3._rt*q(domlo+1,j,k,n) -  0.1_rt*q(domlo+2,j,k,n);
           dlft = 2._rt*(q(domlo  ,j,k,n)-q(domlo-1,j,k,n));
           drgt = 2._rt*(q(domlo+1,j,k,n)-q(domlo  ,j,k,n));
           dlimsh = (dlft*drgt >= 0.0_rt) ? amrex::min(Math::abs(dlft), Math::abs(drgt)) : 0.0_rt;
           dsgnsh = Math::copysign(1.e0_rt, dfm);
           dfm = dsgnsh*amrex::min(dlimsh, Math::abs(dfm));
           dtemp  = 4.0_rt/3.0_rt*dcen - 1.0_rt/6.0_rt*(dfp + dfm);
        }

        if (edhi && i == domhi) {
           dtemp  = 16._rt/15._rt*q(i+1,j,k,n) - .5_rt*q(i,j,k,n) - 2._rt/3._rt*q(i-1,j,k,n) +  0.1_rt*q(i-2,j,k,n);
           dlft = 2._rt*(q(i  ,j,k,n)-q(i-1,j,k,n));
           drgt = 2._rt*(q(i+1,j,k,n)-q(i  ,j,k,n));
           dlim = (dlft*drgt >= 0.0_rt) ? amrex::min(Math::abs(dlft), Math::abs(drgt)) : 0.0_rt;
           dsgn = Math::copysign(1.e0_rt, dtemp);
        } else if (edhi && i == domhi-1) {
           dfp  = 16._rt/15._rt*q(domhi+1,j,k,n) - .5_rt*q(domhi,j,k,n) - 2._rt/3._rt*q(domhi-1,j,k,n) +  0.1_rt*q(domhi-2,j,k,n);
           dlft = 2._rt*(q(domhi  ,j,k,n)-q(domhi-1,j,k,n));
           drgt = 2._rt*(q(domhi+1,j,k,n)-q(domhi  ,j,k,n));
           dlimsh = (dlft*drgt >= 0.0_rt) ? amrex::min(Math::abs(dlft), Math::abs(drgt)) : 0.0_rt;
           dsgnsh = Math::copysign(1.e0_rt, dfp);
           dfp = dsgnsh*amrex::min(dlimsh, Math::abs(dfp));
           dtemp  = 4.0_rt/3.0_rt*dcen - 1.0_rt/6.0_rt*(dfp + dfm);
        }

        return dsgn*amrex::min(dlim, Math::abs(dtemp));

    } else {
        return 0.;
    }

}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
Real amrex_calc_yslope (int i, int j, int k, int n, int order,
                               amrex::Array4<Real const> const& q) noexcept
{
    if (order == 2)
    {
        Real dl = 2.0_rt*(q(i,j  ,k,n) - q(i,j-1,k,n));
        Real dr = 2.0_rt*(q(i,j+1,k,n) - q(i,j  ,k,n));
        Real dc = 0.5_rt*(q(i,j+1,k,n) - q(i,j-1,k,n));
        Real slope = amrex::min(Math::abs(dl),Math::abs(dc),Math::abs(dr));
        slope = (dr*dl > 0.0_rt) ? slope : 0.0_rt;
        return (dc > 0.0_rt) ? slope : -slope;

    } else if (order == 4) {

        Real dlft, drgt, dcen, dfm, dfp, dlim, dsgn, dtemp;
        Real qm, qp, qj;
        qj = q(i,j,k,n);
        qm = q(i,j-1,k,n);
        qp = q(i,j+1,k,n);

        dlft = qm - q(i,j-2,k,n);
        drgt = qj - qm;
        dcen = 0.5_rt*(dlft+drgt);
        dsgn = Math::copysign(1.e0_rt, dcen);
        dlim = (dlft*drgt >= 0.0_rt) ? 2.0_rt*amrex::min(Math::abs(dlft), Math::abs(drgt)) : 0.0_rt;
        dfm = dsgn*amrex::min(dlim, Math::abs(dcen));

        dlft = qp - qj;
        drgt = q(i,j+2,k,n) - qp;
        dcen = 0.5_rt*(dlft+drgt);
        dsgn = Math::copysign(1.e0_rt, dcen);
        dlim = (dlft*drgt >= 0.0_rt) ? 2.0_rt*amrex::min(Math::abs(dlft), Math::abs(drgt)) : 0.0_rt;
        dfp = dsgn*amrex::min(dlim, Math::abs(dcen));

        dlft = qj - qm;
        drgt = qp - qj;
        dcen = 0.5_rt*(dlft+drgt);
        dsgn = Math::copysign(1.e0_rt, dcen);
        dlim = (dlft*drgt >= 0.0_rt) ? 2.0_rt*amrex::min(Math::abs(dlft), Math::abs(drgt)) : 0.0_rt;

        dtemp  = 4.0_rt/3.0_rt*dcen - 1.0_rt/6.0_rt*(dfp + dfm);
        return dsgn*amrex::min(dlim, Math::abs(dtemp));

    } else {
        return 0.;
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
Real amrex_calc_yslope_extdir (int i, int j, int k, int n, int order,
                                      amrex::Array4<Real const> const& q,
                                      bool edlo, bool edhi, int domlo, int domhi) noexcept
{
    if (order == 2)
    {
        Real dl = 2.0_rt*(q(i,j  ,k,n) - q(i,j-1,k,n));
        Real dr = 2.0_rt*(q(i,j+1,k,n) - q(i,j  ,k,n));
        Real dc = 0.5_rt*(q(i,j+1,k,n) - q(i,j-1,k,n));
        if (edlo && j == domlo) {
            dc = (q(i,j+1,k,n)+3.0_rt*q(i,j,k,n)-4.0_rt*q(i,j-1,k,n))/3.0_rt;
        } else if (edhi && j == domhi) {
            dc = (4.0_rt*q(i,j+1,k,n)-3.0_rt*q(i,j,k,n)-q(i,j-1,k,n))/3.0_rt;
        }
        Real slope = amrex::min(Math::abs(dl),Math::abs(dc),Math::abs(dr));
        slope = (dr*dl > 0.0_rt) ? slope : 0.0_rt;
        return (dc > 0.0_rt) ? slope : -slope;

    } else if (order == 4) {

        Real dlft, drgt, dcen, dfm, dfp, dlim, dsgn, dtemp, dlimsh,dsgnsh;
        Real qm, qp, qj;
        qj = q(i,j,k,n);
        qm = q(i,j-1,k,n);
        qp = q(i,j+1,k,n);

        dlft = qm - q(i,j-2,k,n);
        drgt = qj - qm;
        dcen = 0.5_rt*(dlft+drgt);
        dsgn = Math::copysign(1.e0_rt, dcen);
        dlim = (dlft*drgt >= 0.0_rt) ? 2.0_rt*amrex::min(Math::abs(dlft), Math::abs(drgt)) : 0.0_rt;
        dfm = dsgn*amrex::min(dlim, Math::abs(dcen));

        dlft = qp - qj;
        drgt = q(i,j+2,k,n) - qp;
        dcen = 0.5_rt*(dlft+drgt);
        dsgn = Math::copysign(1.e0_rt, dcen);
        dlim = (dlft*drgt >= 0.0_rt) ? 2.0_rt*amrex::min(Math::abs(dlft), Math::abs(drgt)) : 0.0_rt;
        dfp = dsgn*amrex::min(dlim, Math::abs(dcen));

        dlft = qj - qm;
        drgt = qp - qj;
        dcen = 0.5_rt*(dlft+drgt);
        dsgn = Math::copysign(1.e0_rt, dcen);
        dlim = (dlft*drgt >= 0.0_rt) ? 2.0_rt*amrex::min(Math::abs(dlft), Math::abs(drgt)) : 0.0_rt;

        dtemp  = 4.0_rt/3.0_rt*dcen - 1.0_rt/6.0_rt*(dfp + dfm);

        if (edlo && j == domlo) {
           dtemp  = -16._rt/15._rt*q(i,j-1,k,n) + .5_rt*q(i,j,k,n) + 2._rt/3._rt*q(i,j+1,k,n) -  0.1_rt*q(i,j+2,k,n);
           dlft = 2._rt*(q(i  ,j,k,n)-q(i,j-1,k,n));
           drgt = 2._rt*(q(i,j+1,k,n)-q(i  ,j,k,n));
           dlim = (dlft*drgt >= 0.0_rt) ? amrex::min(Math::abs(dlft), Math::abs(drgt)) : 0.0_rt;
           dsgn = Math::copysign(1.e0_rt, dtemp);
        } else if (edlo && j == domlo+1) {
           dfm  = -16._rt/15._rt*q(i,domlo-1,k,n) + .5_rt*q(i,domlo,k,n) + 2._rt/3._rt*q(i,domlo+1,k,n) -  0.1_rt*q(i,domlo+2,k,n);
           dlft = 2._rt*(q(i  ,domlo,k,n)-q(i,domlo-1,k,n));
           drgt = 2._rt*(q(i,domlo+1,k,n)-q(i  ,domlo,k,n));
           dlimsh = (dlft*drgt >= 0.0_rt) ? amrex::min(Math::abs(dlft), Math::abs(drgt)) : 0.0_rt;
           dsgnsh = Math::copysign(1.e0_rt, dfm);
           dfm = dsgnsh*amrex::min(dlimsh, Math::abs(dfm));
           dtemp  = 4.0_rt/3.0_rt*dcen - 1.0_rt/6.0_rt*(dfp + dfm);
        }

        if (edhi && j == domhi) {
           dtemp  = 16._rt/15._rt*q(i,j+1,k,n) - .5_rt*q(i,j,k,n) - 2._rt/3._rt*q(i,j-1,k,n) +  0.1_rt*q(i,j-2,k,n);
           dlft = 2._rt*(q(i  ,j,k,n)-q(i,j-1,k,n));
           drgt = 2._rt*(q(i,j+1,k,n)-q(i  ,j,k,n));
           dlim = (dlft*drgt >= 0.0_rt) ? amrex::min(Math::abs(dlft), Math::abs(drgt)) : 0.0_rt;
           dsgn = Math::copysign(1.e0_rt, dtemp);
        } else if (edhi && j == domhi-1) {
           dfp  = 16._rt/15._rt*q(i,domhi+1,k,n) - .5_rt*q(i,domhi,k,n) - 2._rt/3._rt*q(i,domhi-1,k,n) +  0.1_rt*q(i,domhi-2,k,n);
           dlft = 2._rt*(q(i  ,domhi,k,n)-q(i,domhi-1,k,n));
           drgt = 2._rt*(q(i,domhi+1,k,n)-q(i  ,domhi,k,n));
           dlimsh = (dlft*drgt >= 0.0_rt) ? amrex::min(Math::abs(dlft), Math::abs(drgt)) : 0.0_rt;
           dsgnsh = Math::copysign(1.e0_rt, dfp);
           dfp = dsgnsh*amrex::min(dlimsh, Math::abs(dfp));
           dtemp  = 4.0_rt/3.0_rt*dcen - 1.0_rt/6.0_rt*(dfp + dfm);
        }

        return dsgn*amrex::min(dlim, Math::abs(dtemp));

    } else {
        return 0.;
    }
}

#if (AMREX_SPACEDIM == 3)
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
Real amrex_calc_zslope (int i, int j, int k, int n, int order,
                               amrex::Array4<Real const> const& q) noexcept
{
    if (order == 2)
    {
        Real dl = 2.0_rt*(q(i,j,k  ,n) - q(i,j,k-1,n));
        Real dr = 2.0_rt*(q(i,j,k+1,n) - q(i,j,k  ,n));
        Real dc = 0.5_rt*(q(i,j,k+1,n) - q(i,j,k-1,n));
        Real slope = amrex::min(Math::abs(dl),Math::abs(dc),Math::abs(dr));
        slope = (dr*dl > 0.0_rt) ? slope : 0.0_rt;
        return (dc > 0.0_rt) ? slope : -slope;

    } else if (order == 4) {

        Real dlft, drgt, dcen, dfm, dfp, dlim, dsgn, dtemp;
        Real qm, qp, qk;
        qk = q(i,j,k,n);
        qm = q(i,j,k-1,n);
        qp = q(i,j,k+1,n);

        dlft = qm - q(i,j,k-2,n);
        drgt = qk - qm;
        dcen = 0.5_rt*(dlft+drgt);
        dsgn = Math::copysign(1.e0_rt, dcen);
        dlim = (dlft*drgt >= 0.0_rt) ? 2.0_rt*amrex::min(Math::abs(dlft), Math::abs(drgt)) : 0.0_rt;
        dfm = dsgn*amrex::min(dlim, Math::abs(dcen));

        dlft = qp - qk;
        drgt = q(i,j,k+2,n) - qp;
        dcen = 0.5_rt*(dlft+drgt);
        dsgn = Math::copysign(1.e0_rt, dcen);
        dlim = (dlft*drgt >= 0.0_rt) ? 2.0_rt*amrex::min(Math::abs(dlft), Math::abs(drgt)) : 0.0_rt;
        dfp = dsgn*amrex::min(dlim, Math::abs(dcen));

        dlft = qk - qm;
        drgt = qp - qk;
        dcen = 0.5_rt*(dlft+drgt);
        dsgn = Math::copysign(1.e0_rt, dcen);
        dlim = (dlft*drgt >= 0.0_rt) ? 2.0_rt*amrex::min(Math::abs(dlft), Math::abs(drgt)) : 0.0_rt;

        dtemp  = 4.0_rt/3.0_rt*dcen - 1.0_rt/6.0_rt*(dfp + dfm);
        return dsgn*amrex::min(dlim, Math::abs(dtemp));

    } else {
        return 0.;
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
Real amrex_calc_zslope_extdir (int i, int j, int k, int n, int order,
                                      amrex::Array4<Real const> const& q,
                                      bool edlo, bool edhi, int domlo, int domhi) noexcept
{
    if (order == 2)
    {

        Real dl = 2.0_rt*(q(i,j,k  ,n) - q(i,j,k-1,n));
        Real dr = 2.0_rt*(q(i,j,k+1,n) - q(i,j,k  ,n));
        Real dc = 0.5_rt*(q(i,j,k+1,n) - q(i,j,k-1,n));
        if (edlo && k == domlo) {
            dc = (q(i,j,k+1,n)+3.0_rt*q(i,j,k,n)-4.0_rt*q(i,j,k-1,n))/3.0_rt;
        } else if (edhi && k == domhi) {
            dc = (4.0_rt*q(i,j,k+1,n)-3.0_rt*q(i,j,k,n)-q(i,j,k-1,n))/3.0_rt;
        }
        Real slope = amrex::min(Math::abs(dl),Math::abs(dc),Math::abs(dr));
        slope = (dr*dl > 0.0_rt) ? slope : 0.0_rt;
        return (dc > 0.0_rt) ? slope : -slope;

    } else if (order == 4) {

        Real dlft, drgt, dcen, dfm, dfp, dlim, dsgn, dtemp, dlimsh, dsgnsh;
        Real qm, qp, qk;
        qk = q(i,j,k,n);
        qm = q(i,j,k-1,n);
        qp = q(i,j,k+1,n);

        dlft = qm - q(i,j,k-2,n);
        drgt = qk - qm;
        dcen = 0.5_rt*(dlft+drgt);
        dsgn = Math::copysign(1.e0_rt, dcen);
        dlim = (dlft*drgt >= 0.0_rt) ? 2.0_rt*amrex::min(Math::abs(dlft), Math::abs(drgt)) : 0.0_rt;
        dfm = dsgn*amrex::min(dlim, Math::abs(dcen));

        dlft = qp - qk;
        drgt = q(i,j,k+2,n) - qp;
        dcen = 0.5_rt*(dlft+drgt);
        dsgn = Math::copysign(1.e0_rt, dcen);
        dlim = (dlft*drgt >= 0.0_rt) ? 2.0_rt*amrex::min(Math::abs(dlft), Math::abs(drgt)) : 0.0_rt;
        dfp = dsgn*amrex::min(dlim, Math::abs(dcen));

        dlft = qk - qm;
        drgt = qp - qk;
        dcen = 0.5_rt*(dlft+drgt);
        dsgn = Math::copysign(1.e0_rt, dcen);
        dlim = (dlft*drgt >= 0.0_rt) ? 2.0_rt*amrex::min(Math::abs(dlft), Math::abs(drgt)) : 0.0_rt;

        dtemp  = 4.0_rt/3.0_rt*dcen - 1.0_rt/6.0_rt*(dfp + dfm);

        if (edlo && k == domlo) {
           dtemp  = -16._rt/15._rt*q(i,j,k-1,n) + .5_rt*q(i,j,k,n) + 2._rt/3._rt*q(i,j,k+1,n) -  0.1_rt*q(i,j,k+2,n);
           dlft = 2._rt*(q(i  ,j,k,n)-q(i,j,k-1,n));
           drgt = 2._rt*(q(i,j,k+1,n)-q(i  ,j,k,n));
           dlim = (dlft*drgt >= 0.0_rt) ? amrex::min(Math::abs(dlft), Math::abs(drgt)) : 0.0_rt;
           dsgn = Math::copysign(1.e0_rt, dtemp);
        } else if (edlo && k == domlo+1) {
           dfm  = -16._rt/15._rt*q(i,j,domlo-1,n) + .5_rt*q(i,j,domlo,n) + 2._rt/3._rt*q(i,j,domlo+1,n) -  0.1_rt*q(i,j,domlo+2,n);
           dlft = 2._rt*(q(i  ,j,domlo,n)-q(i,j,domlo-1,n));
           drgt = 2._rt*(q(i,j,domlo+1,n)-q(i  ,j,domlo,n));
           dlimsh = (dlft*drgt >= 0.0_rt) ? amrex::min(Math::abs(dlft), Math::abs(drgt)) : 0.0_rt;
           dsgnsh = Math::copysign(1.e0_rt, dfm);
           dfm = dsgnsh*amrex::min(dlimsh, Math::abs(dfm));
           dtemp  = 4.0_rt/3.0_rt*dcen - 1.0_rt/6.0_rt*(dfp + dfm);
        }

        if (edhi && k == domhi) {
           dtemp  = 16._rt/15._rt*q(i,j,k+1,n) - .5_rt*q(i,j,k,n) - 2._rt/3._rt*q(i,j,k-1,n) +  0.1_rt*q(i,j,k-2,n);
           dlft = 2._rt*(q(i  ,j,k,n)-q(i,j,k-1,n));
           drgt = 2._rt*(q(i,j,k+1,n)-q(i  ,j,k,n));
           dlim = (dlft*drgt >= 0.0_rt) ? amrex::min(Math::abs(dlft), Math::abs(drgt)) : 0.0_rt;
           dsgn = Math::copysign(1.e0_rt, dtemp);
        } else if (edhi && k == domhi-1) {
           dfp  = 16._rt/15._rt*q(i,j,domhi+1,n) - .5_rt*q(i,j,domhi,n) - 2._rt/3._rt*q(i,j,domhi-1,n) +  0.1_rt*q(i,j,domhi-2,n);
           dlft = 2._rt*(q(i  ,j,domhi,n)-q(i,j,domhi-1,n));
           drgt = 2._rt*(q(i,j,domhi+1,n)-q(i  ,j,domhi,n));
           dlimsh = (dlft*drgt >= 0.0_rt) ? amrex::min(Math::abs(dlft), Math::abs(drgt)) : 0.0_rt;
           dsgnsh = Math::copysign(1.e0_rt, dfp);
           dfp = dsgnsh*amrex::min(dlimsh, Math::abs(dfp));
           dtemp  = 4.0_rt/3.0_rt*dcen - 1.0_rt/6.0_rt*(dfp + dfm);
        }
        return dsgn*amrex::min(dlim, Math::abs(dtemp));

    } else {
        return 0.;
    }
}
#endif

}
#endif
