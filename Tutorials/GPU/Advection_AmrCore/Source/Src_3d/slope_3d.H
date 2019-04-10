#ifndef _slope_3d_H_
#define _slope_3d_H_


#include <AMReX_Box.H>
#include <AMReX_Geometry.H>
#include <AMReX_Gpu.H>
#include <AmrCoreAdv_F.H>

using namespace amrex;

AMREX_GPU_DEVICE
AMREX_INLINE
void slopex2(Box const& bx,
             Array4<Real> const& q,
             Array4<Real> const& dq)
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                Real dlft = q(i,j,k) - q(i-1,j,k);
                Real drgt = q(i+1,j,k) - q(i,j,k);
                Real dcen = 0.5*(dlft+drgt);
                Real dsgn = copysign(1.0, dcen);
                Real dslop = 2.0 * ((std::abs(dlft) < std::abs(drgt)) ?
                                     std::abs(dlft) : std::abs(drgt));
                Real dlim = (dlft*drgt >= 0.0) ? dslop : 0.0;
                dq(i,j,k) = dsgn*amrex::min(dlim, std::abs(dcen));
            }
        }
    }
}

AMREX_GPU_DEVICE
AMREX_INLINE
void slopex4(Box const& bx,
             Array4<Real> const& q,
             Array4<Real> const& dq,
             Array4<Real> const& dq4)
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                Real dlft = q(i,j,k) - q(i-1,j,k);
                Real drgt = q(i+1,j,k) - q(i,j,k);
                Real dcen = 0.5*(dlft+drgt);
                Real dsgn = copysign(1.0, dcen);
                Real dslop = 2.0 * ((std::abs(dlft) < std::abs(drgt)) ?
                                     std::abs(dlft) : std::abs(drgt));
                Real dlim = (dlft*drgt >= 0.0) ? dslop : 0.0;
                Real dq1 = 4.0/3.0*dcen - (1.0/6.0)*(dq(i+1,j,k) + dq(i-1,j,k));
                dq4(i,j,k) = dsgn*amrex::min(dlim, std::abs(dq1));
            }
        }
    }
}

// ***********************************************************

AMREX_GPU_DEVICE
AMREX_INLINE
void slopey2(Box const& bx,
             Array4<Real> const& q,
             Array4<Real> const& dq)
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                Real dlft = q(i,j,k) - q(i,j-1,k);
                Real drgt = q(i,j+1,k) - q(i,j,k);
                Real dcen = 0.5*(dlft+drgt);
                Real dsgn = copysign(1.0, dcen);
                Real dslop = 2.0 * ((std::abs(dlft) < std::abs(drgt)) ?
                                     std::abs(dlft) : std::abs(drgt));
                Real dlim = (dlft*drgt >= 0.0) ? dslop : 0.0;
                dq(i,j,k) = dsgn*amrex::min(dlim, std::abs(dcen));
            }
        }
    }
}

AMREX_GPU_DEVICE
AMREX_INLINE
void slopey4(Box const& bx,
             Array4<Real> const& q,
             Array4<Real> const& dq,
             Array4<Real> const& dq4)
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                Real dlft = q(i,j,k) - q(i,j-1,k);
                Real drgt = q(i,j+1,k) - q(i,j,k);
                Real dcen = 0.5*(dlft+drgt);
                Real dsgn = copysign(1.0, dcen);
                Real dslop = 2.0 * ((std::abs(dlft) < std::abs(drgt)) ?
                                     std::abs(dlft) : std::abs(drgt));
                Real dlim = (dlft*drgt >= 0.0) ? dslop : 0.0;
                Real dq1 = 4.0/3.0*dcen - (1.0/6.0)*(dq(i,j+1,k) + dq(i,j-1,k));
                dq4(i,j,k) = dsgn*amrex::min(dlim, std::abs(dq1));
            }
        }
    }
}

// ***********************************************************

AMREX_GPU_DEVICE
AMREX_INLINE
void slopez2(Box const& bx,
             Array4<Real> const&  q,
             Array4<Real> const& dq)
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                Real dlft = q(i,j,k) - q(i,j,k-1);
                Real drgt = q(i,j,k+1) - q(i,j,k);
                Real dcen = 0.5*(dlft+drgt);
                Real dsgn = copysign(1.0, dcen);
                Real dslop = 2.0 * ((std::abs(dlft) < std::abs(drgt)) ?
                                     std::abs(dlft) : std::abs(drgt));
                Real dlim = (dlft*drgt >= 0.0) ? dslop : 0.0;
                dq(i,j,k) = dsgn*amrex::min(dlim, std::abs(dcen));
            }
        }
    }
}

AMREX_GPU_DEVICE
AMREX_INLINE
void slopez4(Box const& bx,
             Array4<Real> const& q,
             Array4<Real> const& dq,
             Array4<Real> const& dq4)
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                Real dlft = q(i,j,k) - q(i,j,k-1);
                Real drgt = q(i,j,k+1) - q(i,j,k);
                Real dcen = 0.5*(dlft+drgt);
                Real dsgn = copysign(1.0, dcen);
                Real dslop = 2.0 * ((std::abs(dlft) < std::abs(drgt)) ?
                                     std::abs(dlft) : std::abs(drgt));
                Real dlim = (dlft*drgt >= 0.0) ? dslop : 0.0;
                Real dq1 = 4.0/3.0*dcen - (1.0/6.0)*(dq(i,j,k+1) + dq(i,j,k-1));
                dq4(i,j,k) = dsgn*amrex::min(dlim, std::abs(dq1));
            }
        }
    }
}

#endif
