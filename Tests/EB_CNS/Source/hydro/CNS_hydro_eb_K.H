#ifndef CNS_HYDRO_EB_K_H_
#define CNS_HYDRO_EB_K_H_

#include "CNS_index_macros.H"
#include "CNS.H"
#include "CNS_parm.H"
#include "CNS_hydro_K.H"
#include <AMReX_FArrayBox.H>
#include <cmath>

namespace {

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real limiter_eb (amrex::Real dlft, amrex::Real drgt, amrex::Real plm_theta) noexcept
{
    using amrex::Real;

    Real dcen = Real(0.5)*(dlft+drgt);
    Real dsgn = std::copysign(Real(1.0), dcen);
    Real slop = plm_theta * amrex::min(std::abs(dlft),std::abs(drgt));
    Real dlim = (dlft*drgt >= Real(0.0)) ? slop : Real(0.0);
    return dsgn * amrex::min(dlim,std::abs(dcen));
}

}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
cns_slope_eb_x (int i, int j, int k,
                amrex::Array4<amrex::Real> const& dq,
                amrex::Array4<amrex::Real const> const& q,
                amrex::Array4<amrex::EBCellFlag const> const& flag,
                int plm_iorder, amrex::Real plm_theta) noexcept
{
    using amrex::Real;

    Real dlft0 = Real(0.0);
    Real dlft1 = Real(0.0);
    Real dlft2 = Real(0.0);
    Real dlft3 = Real(0.0);
#if (AMREX_SPACEDIM == 3)
    Real dlft4 = Real(0.0);
#endif

    Real drgt0 = Real(0.0);
    Real drgt1 = Real(0.0);
    Real drgt2 = Real(0.0);
    Real drgt3 = Real(0.0);
#if (AMREX_SPACEDIM == 3)
    Real drgt4 = Real(0.0);
#endif

    if (plm_iorder == 1)
    {
        dq(i,j,k,0) = Real(0.);
        dq(i,j,k,1) = Real(0.);
        dq(i,j,k,2) = Real(0.);
        dq(i,j,k,3) = Real(0.);
#if (AMREX_SPACEDIM == 3)
        dq(i,j,k,4) = Real(0.);
#endif

    } else {

        if (flag(i,j,k).isConnected(-1,0,0))
        {
            dlft0 = Real(0.5)*(q(i,j,k,QPRES)-q(i-1,j,k,QPRES))/q(i,j,k,QCS) - Real(0.5)*q(i,j,k,QRHO)*(q(i,j,k,QU) - q(i-1,j,k,QU));
            dlft1 = (q(i,j,k,QRHO)-q(i-1,j,k,QRHO)) - (q(i,j,k,QPRES) - q(i-1,j,k,QPRES))/(q(i,j,k,QCS)*q(i,j,k,QCS));
            dlft2 = Real(0.5)*(q(i,j,k,QPRES)-q(i-1,j,k,QPRES))/q(i,j,k,QCS) + Real(0.5)*q(i,j,k,QRHO)*(q(i,j,k,QU) - q(i-1,j,k,QU));
            dlft3 = q(i,j,k,QV) - q(i-1,j,k,QV);
#if (AMREX_SPACEDIM == 3)
            dlft4 = q(i,j,k,QW) - q(i-1,j,k,QW);
#endif
        }

        if (flag(i,j,k).isConnected(1,0,0))
        {
            drgt0 = Real(0.5)*(q(i+1,j,k,QPRES)-q(i,j,k,QPRES))/q(i,j,k,QCS) - Real(0.5)*q(i,j,k,QRHO)*(q(i+1,j,k,QU) - q(i,j,k,QU));
            drgt1 = (q(i+1,j,k,QRHO)-q(i,j,k,QRHO)) - (q(i+1,j,k,QPRES) - q(i,j,k,QPRES))/(q(i,j,k,QCS)*q(i,j,k,QCS));
            drgt2 = Real(0.5)*(q(i+1,j,k,QPRES)-q(i,j,k,QPRES))/q(i,j,k,QCS) + Real(0.5)*q(i,j,k,QRHO)*(q(i+1,j,k,QU) - q(i,j,k,QU));
            drgt3 = q(i+1,j,k,QV) - q(i,j,k,QV);
#if (AMREX_SPACEDIM == 3)
            drgt4 = q(i+1,j,k,QW) - q(i,j,k,QW);
#endif
        }

        dq(i,j,k,0) = limiter_eb(dlft0, drgt0, plm_theta);
        dq(i,j,k,1) = limiter_eb(dlft1, drgt1, plm_theta);
        dq(i,j,k,2) = limiter_eb(dlft2, drgt2, plm_theta);
        dq(i,j,k,3) = limiter_eb(dlft3, drgt3, plm_theta);
#if (AMREX_SPACEDIM == 3)
        dq(i,j,k,4) = limiter_eb(dlft4, drgt4, plm_theta);
#endif
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
cns_slope_eb_y (int i, int j, int k,
                amrex::Array4<amrex::Real> const& dq,
                amrex::Array4<amrex::Real const> const& q,
                amrex::Array4<amrex::EBCellFlag const> const& flag,
                int plm_iorder, amrex::Real plm_theta) noexcept
{
    using amrex::Real;

    Real dlft0 = Real(0.0);
    Real dlft1 = Real(0.0);
    Real dlft2 = Real(0.0);
    Real dlft3 = Real(0.0);
#if (AMREX_SPACEDIM == 3)
    Real dlft4 = Real(0.0);
#endif

    Real drgt0 = Real(0.0);
    Real drgt1 = Real(0.0);
    Real drgt2 = Real(0.0);
    Real drgt3 = Real(0.0);
#if (AMREX_SPACEDIM == 3)
    Real drgt4 = Real(0.0);
#endif

    if (plm_iorder == 1)
    {
        dq(i,j,k,0) = Real(0.);
        dq(i,j,k,1) = Real(0.);
        dq(i,j,k,2) = Real(0.);
        dq(i,j,k,3) = Real(0.);
#if (AMREX_SPACEDIM == 3)
        dq(i,j,k,4) = Real(0.);
#endif

    } else {

        if (flag(i,j,k).isConnected(0,-1,0))
        {
            dlft0 = Real(0.5)*(q(i,j,k,QPRES)-q(i,j-1,k,QPRES))/q(i,j,k,QCS) - Real(0.5)*q(i,j,k,QRHO)*(q(i,j,k,QV) - q(i,j-1,k,QV));
            dlft1 = (q(i,j,k,QRHO)-q(i,j-1,k,QRHO)) - (q(i,j,k,QPRES) - q(i,j-1,k,QPRES))/(q(i,j,k,QCS)*q(i,j,k,QCS));
            dlft2 = Real(0.5)*(q(i,j,k,QPRES)-q(i,j-1,k,QPRES))/q(i,j,k,QCS) + Real(0.5)*q(i,j,k,QRHO)*(q(i,j,k,QV) - q(i,j-1,k,QV));
            dlft3 =q(i,j,k,QU) - q(i,j-1,k,QU);
#if (AMREX_SPACEDIM == 3)
            dlft4 = q(i,j,k,QW) - q(i,j-1,k,QW);
#endif
        }

        if (flag(i,j,k).isConnected(0,1,0))
        {
            drgt0 = Real(0.5)*(q(i,j+1,k,QPRES)-q(i,j,k,QPRES))/q(i,j,k,QCS) - Real(0.5)*q(i,j,k,QRHO)*(q(i,j+1,k,QV) - q(i,j,k,QV));
            drgt1 = (q(i,j+1,k,QRHO)-q(i,j,k,QRHO)) - (q(i,j+1,k,QPRES) - q(i,j,k,QPRES))/(q(i,j,k,QCS)*q(i,j,k,QCS));
            drgt2 = Real(0.5)*(q(i,j+1,k,QPRES)-q(i,j,k,QPRES))/q(i,j,k,QCS) + Real(0.5)*q(i,j,k,QRHO)*(q(i,j+1,k,QV) - q(i,j,k,QV));
            drgt3 = q(i,j+1,k,QU) - q(i,j,k,QU);
#if (AMREX_SPACEDIM == 3)
            drgt4 = q(i,j+1,k,QW) - q(i,j,k,QW);
#endif
        }

        dq(i,j,k,0) = limiter_eb(dlft0, drgt0, plm_theta);
        dq(i,j,k,1) = limiter_eb(dlft1, drgt1, plm_theta);
        dq(i,j,k,2) = limiter_eb(dlft2, drgt2, plm_theta);
        dq(i,j,k,3) = limiter_eb(dlft3, drgt3, plm_theta);
#if (AMREX_SPACEDIM == 3)
        dq(i,j,k,4) = limiter_eb(dlft4, drgt4, plm_theta);
#endif
   }
}

#if (AMREX_SPACEDIM == 3)
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
cns_slope_eb_z (int i, int j, int k,
                amrex::Array4<amrex::Real> const& dq,
                amrex::Array4<amrex::Real const> const& q,
                amrex::Array4<amrex::EBCellFlag const> const& flag,
                int plm_iorder, amrex::Real plm_theta) noexcept
{
    using amrex::Real;

    Real dlft0, dlft1, dlft2, dlft3, dlft4;
    Real drgt0, drgt1, drgt2, drgt3, drgt4;

    if (plm_iorder == 1)
    {
        dq(i,j,k,0) = Real(0.);
        dq(i,j,k,1) = Real(0.);
        dq(i,j,k,2) = Real(0.);
        dq(i,j,k,3) = Real(0.);
        dq(i,j,k,4) = Real(0.);

    } else {

        if (flag(i,j,k).isConnected(0,0,-1))
        {
            dlft0 = Real(0.5)*(q(i,j,k,QPRES)-q(i,j,k-1,QPRES))/q(i,j,k,QCS) - Real(0.5)*q(i,j,k,QRHO)*(q(i,j,k,QW) - q(i,j,k-1,QW));
            dlft1 = (q(i,j,k,QRHO)-q(i,j,k-1,QRHO)) - (q(i,j,k,QPRES) - q(i,j,k-1,QPRES))/(q(i,j,k,QCS)*q(i,j,k,QCS));
            dlft2 = Real(0.5)*(q(i,j,k,QPRES)-q(i,j,k-1,QPRES))/q(i,j,k,QCS) + Real(0.5)*q(i,j,k,QRHO)*(q(i,j,k,QW) - q(i,j,k-1,QW));
            dlft3 = q(i,j,k,QU) - q(i,j,k-1,QU);
            dlft4 = q(i,j,k,QV) - q(i,j,k-1,QV);
        } else {
            dlft0 = Real(0.0);
            dlft1 = Real(0.0);
            dlft2 = Real(0.0);
            dlft3 = Real(0.0);
            dlft4 = Real(0.0);
        }

        if (flag(i,j,k).isConnected(0,0,1))
        {
            drgt0 = Real(0.5)*(q(i,j,k+1,QPRES)-q(i,j,k,QPRES))/q(i,j,k,QCS) - Real(0.5)*q(i,j,k,QRHO)*(q(i,j,k+1,QW) - q(i,j,k,QW));
            drgt1 = (q(i,j,k+1,QRHO)-q(i,j,k,QRHO)) - (q(i,j,k+1,QPRES) - q(i,j,k,QPRES))/(q(i,j,k,QCS)*q(i,j,k,QCS));
            drgt2 = Real(0.5)*(q(i,j,k+1,QPRES)-q(i,j,k,QPRES))/q(i,j,k,QCS) + Real(0.5)*q(i,j,k,QRHO)*(q(i,j,k+1,QW) - q(i,j,k,QW));
            drgt3 = q(i,j,k+1,QU) - q(i,j,k,QU);
            drgt4 = q(i,j,k+1,QV) - q(i,j,k,QV);
        } else {
            drgt0 = Real(0.0);
            drgt1 = Real(0.0);
            drgt2 = Real(0.0);
            drgt3 = Real(0.0);
            drgt4 = Real(0.0);
        }


        dq(i,j,k,0) = limiter_eb(dlft0, drgt0, plm_theta);
        dq(i,j,k,1) = limiter_eb(dlft1, drgt1, plm_theta);
        dq(i,j,k,2) = limiter_eb(dlft2, drgt2, plm_theta);
        dq(i,j,k,3) = limiter_eb(dlft3, drgt3, plm_theta);
        dq(i,j,k,4) = limiter_eb(dlft4, drgt4, plm_theta);
   }
}
#endif
#endif
