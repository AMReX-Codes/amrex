

#include <slope_3d.H>
#include <AMReX_Gpu.H>
#include <AmrCoreAdv_F.H>

using namespace amrex;

AMREX_GPU_DEVICE
void slopex2(Box const& bx,
            const FArrayBox &qfab,
            FArrayBox &dqfab)
{
    const auto len = length(bx);
    const auto lo  = lbound(bx);
    const auto q   =  qfab.view(lo);
    const auto dq  = dqfab.view(lo);

    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            for (int i = 0; i < len.x; ++i) {
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
void slopex4(Box const& bx,
             const FArrayBox &qfab,
             const FArrayBox &dqfab,
             FArrayBox &dq4fab)
{
    const auto len = length(bx);
    const auto lo  = lbound(bx);
    const auto q   =   qfab.view(lo);
    const auto dq  =  dqfab.view(lo);
    const auto dq4 = dq4fab.view(lo);

    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            for (int i = 0; i < len.x; ++i) {
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
void slopey2(Box const& bx,
            const FArrayBox &qfab,
            FArrayBox &dqfab)
{
    const auto len = length(bx);
    const auto lo  = lbound(bx);
    const auto q   =  qfab.view(lo);
    const auto dq  = dqfab.view(lo);

    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            for (int i = 0; i < len.x; ++i) {
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
void slopey4(Box const& bx,
             const FArrayBox &qfab,
             const FArrayBox &dqfab,
             FArrayBox &dq4fab)
{
    const auto len = length(bx);
    const auto lo  = lbound(bx);
    const auto q   =   qfab.view(lo);
    const auto dq  =  dqfab.view(lo);
    const auto dq4 = dq4fab.view(lo);

    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            for (int i = 0; i < len.x; ++i) {
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
void slopez2(Box const& bx,
            const FArrayBox &qfab,
            FArrayBox &dqfab)
{
    const auto len = length(bx);
    const auto lo  = lbound(bx);
    const auto q   =  qfab.view(lo);
    const auto dq  = dqfab.view(lo);

    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            for (int i = 0; i < len.x; ++i) {
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
void slopez4(Box const& bx,
             const FArrayBox &qfab,
             const FArrayBox &dqfab,
             FArrayBox &dq4fab)
{
    const auto len = length(bx);
    const auto lo  = lbound(bx);
    const auto q   =   qfab.view(lo);
    const auto dq  =  dqfab.view(lo);
    const auto dq4 = dq4fab.view(lo);

    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            for (int i = 0; i < len.x; ++i) {
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

