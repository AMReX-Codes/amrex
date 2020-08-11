#ifndef slope_K_H_
#define slope_K_H_

#include <AMReX_Box.H>
#include <AMReX_Geometry.H>
#include <AMReX_Gpu.H>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void slopex2(amrex::Box const& bx,
             amrex::Array4<amrex::Real> const& dq,
             amrex::Array4<amrex::Real const> const& q)
{
    using namespace amrex;

    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                Real dlft = q(i,j,k) - q(i-1,j,k);
                Real drgt = q(i+1,j,k) - q(i,j,k);
                Real dcen = 0.5*(dlft+drgt);
                Real dsgn = amrex::Math::copysign(1.0, dcen);
                Real dslop = 2.0 * ((amrex::Math::abs(dlft) < amrex::Math::abs(drgt)) ?
                                     amrex::Math::abs(dlft) : amrex::Math::abs(drgt));
                Real dlim = (dlft*drgt >= 0.0) ? dslop : 0.0;
                dq(i,j,k) = dsgn*amrex::min(dlim, amrex::Math::abs(dcen));
            }
        }
    }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void slopex4(amrex::Box const& bx,
             amrex::Array4<amrex::Real> const& dq4,
             amrex::Array4<amrex::Real const> const& q,
             amrex::Array4<amrex::Real const> const& dq)
{
    using namespace amrex;

    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                Real dlft = q(i,j,k) - q(i-1,j,k);
                Real drgt = q(i+1,j,k) - q(i,j,k);
                Real dcen = 0.5*(dlft+drgt);
                Real dsgn = amrex::Math::copysign(1.0, dcen);
                Real dslop = 2.0 * ((amrex::Math::abs(dlft) < amrex::Math::abs(drgt)) ?
                                     amrex::Math::abs(dlft) : amrex::Math::abs(drgt));
                Real dlim = (dlft*drgt >= 0.0) ? dslop : 0.0;
                Real dq1 = 4.0/3.0*dcen - (1.0/6.0)*(dq(i+1,j,k) + dq(i-1,j,k));
                dq4(i,j,k) = dsgn*amrex::min(dlim, amrex::Math::abs(dq1));
            }
        }
    }
}

// ***********************************************************

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void slopey2(amrex::Box const& bx,
             amrex::Array4<amrex::Real> const& dq,
             amrex::Array4<amrex::Real const> const& q)
{
    using namespace amrex;

    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                Real dlft = q(i,j,k) - q(i,j-1,k);
                Real drgt = q(i,j+1,k) - q(i,j,k);
                Real dcen = 0.5*(dlft+drgt);
                Real dsgn = amrex::Math::copysign(1.0, dcen);
                Real dslop = 2.0 * ((amrex::Math::abs(dlft) < amrex::Math::abs(drgt)) ?
                                     amrex::Math::abs(dlft) : amrex::Math::abs(drgt));
                Real dlim = (dlft*drgt >= 0.0) ? dslop : 0.0;
                dq(i,j,k) = dsgn*amrex::min(dlim, amrex::Math::abs(dcen));
            }
        }
    }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void slopey4(amrex::Box const& bx,
             amrex::Array4<amrex::Real> const& dq4,
             amrex::Array4<amrex::Real const> const& q,
             amrex::Array4<amrex::Real const> const& dq)
{
    using namespace amrex;

    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                Real dlft = q(i,j,k) - q(i,j-1,k);
                Real drgt = q(i,j+1,k) - q(i,j,k);
                Real dcen = 0.5*(dlft+drgt);
                Real dsgn = amrex::Math::copysign(1.0, dcen);
                Real dslop = 2.0 * ((amrex::Math::abs(dlft) < amrex::Math::abs(drgt)) ?
                                     amrex::Math::abs(dlft) : amrex::Math::abs(drgt));
                Real dlim = (dlft*drgt >= 0.0) ? dslop : 0.0;
                Real dq1 = 4.0/3.0*dcen - (1.0/6.0)*(dq(i,j+1,k) + dq(i,j-1,k));
                dq4(i,j,k) = dsgn*amrex::min(dlim, amrex::Math::abs(dq1));
            }
        }
    }
}

// ***********************************************************

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void slopez2(amrex::Box const& bx,
             amrex::Array4<amrex::Real> const& dq,
             amrex::Array4<amrex::Real const> const& q)
{
    using namespace amrex;

    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                Real dlft = q(i,j,k) - q(i,j,k-1);
                Real drgt = q(i,j,k+1) - q(i,j,k);
                Real dcen = 0.5*(dlft+drgt);
                Real dsgn = amrex::Math::copysign(1.0, dcen);
                Real dslop = 2.0 * ((amrex::Math::abs(dlft) < amrex::Math::abs(drgt)) ?
                                     amrex::Math::abs(dlft) : amrex::Math::abs(drgt));
                Real dlim = (dlft*drgt >= 0.0) ? dslop : 0.0;
                dq(i,j,k) = dsgn*amrex::min(dlim, amrex::Math::abs(dcen));
            }
        }
    }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void slopez4(amrex::Box const& bx,
             amrex::Array4<amrex::Real> const& dq4,
             amrex::Array4<amrex::Real const> const& q,
             amrex::Array4<amrex::Real const> const& dq)
{
    using namespace amrex;

    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                Real dlft = q(i,j,k) - q(i,j,k-1);
                Real drgt = q(i,j,k+1) - q(i,j,k);
                Real dcen = 0.5*(dlft+drgt);
                Real dsgn = amrex::Math::copysign(1.0, dcen);
                Real dslop = 2.0 * ((amrex::Math::abs(dlft) < amrex::Math::abs(drgt)) ?
                                     amrex::Math::abs(dlft) : amrex::Math::abs(drgt));
                Real dlim = (dlft*drgt >= 0.0) ? dslop : 0.0;
                Real dq1 = 4.0/3.0*dcen - (1.0/6.0)*(dq(i,j,k+1) + dq(i,j,k-1));
                dq4(i,j,k) = dsgn*amrex::min(dlim, amrex::Math::abs(dq1));
            }
        }
    }
}

#endif
