
#include "CNS_K.H"
#include "CNS_index_macros.H"

using namespace amrex;

AMREX_GPU_DEVICE
void cns_ctoprim (Box const& bx, FArrayBox const& ufab, FArrayBox & qfab)
{
    const auto len = length(bx);
    const auto lo  = lbound(bx);
    const auto u = ufab.view(lo);
    const auto q = qfab.view(lo);
    const Real smallr = 1.e-19;
    const Real smallp = 1.e-10;
    const Real gamma = 1.4;

    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = 0; i < len.x; ++i) {
                Real rho = u(i,j,k,URHO);
                rho = (rho > smallr) ? rho : smallr;
                Real rhoinv = 1.0/rho;
                Real ux = u(i,j,k,UMX)*rhoinv;
                Real uy = u(i,j,k,UMY)*rhoinv;
                Real uz = u(i,j,k,UMZ)*rhoinv;
                Real kineng = 0.5*rho*(ux*ux+uy*uy+uz*uz);
                Real ei = u(i,j,k,UEDEN) - kineng;
                if (ei <= 0.0) ei = u(i,j,k,UEINT);
                Real p = (gamma-1.0)*ei;
                p = (p > smallp) ? p : smallp;
                ei *= rhoinv;

                q(i,j,k,QRHO) = rho;
                q(i,j,k,QU) = ux;
                q(i,j,k,QV) = uy;
                q(i,j,k,QW) = uz;
                q(i,j,k,QEINT) = ei;
                q(i,j,k,QPRES) = p;
                q(i,j,k,QCS) = std::sqrt(gamma*p*rhoinv);
                q(i,j,k,QTEMP) = 0.0;
            }
        }
    }
}

#if (AMREX_SPACEDIM == 3)
AMREX_GPU_DEVICE
void cns_flux_to_dudt (Box const& bx, FArrayBox& dudtfab,
                       FArrayBox const& fxfab, FArrayBox const& fyfab, FArrayBox const& fzfab,
                       GpuArray<Real,AMREX_SPACEDIM> const& a_dxinv)
{
    const auto len = length(bx);
    const auto lo  = lbound(bx);
    const int ncomp = dudtfab.nComp();
    const auto dudt = dudtfab.view(lo);
    const auto fx   =   fxfab.view(lo);
    const auto fy   =   fyfab.view(lo);
    const auto fz   =   fzfab.view(lo);
    const Real dxinv = a_dxinv[0];
    const Real dyinv = a_dxinv[1];
    const Real dzinv = a_dxinv[2];

    for (int n = 0; n < NCONS; ++n) {
        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = 0; i < len.x; ++i) {
                    dudt(i,j,k,n) = dxinv * (fx(i,j,k,n) - fx(i+1,j,k,n))
                        +           dyinv * (fy(i,j,k,n) - fy(i,j+1,k,n))
                        +           dzinv * (fz(i,j,k,n) - fz(i,j,k+1,n));
                }
            }
        }
    }
}
#endif

AMREX_GPU_DEVICE
void cns_slope_x (Box const& bx, FArrayBox& dqfab, FArrayBox const& qfab)
{
    const auto len = length(bx);
    const auto lo  = lbound(bx);
    const auto dq = dqfab.view(lo);
    const auto  q =  qfab.view(lo);

    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = 0; i < len.x; ++i) {
                auto f = [] (Real dlft, Real drgt) -> Real {
                    Real dcen = 0.5*(dlft+drgt);
                    Real dsgn = std::copysign(1.0, dcen);
                    Real slop = 2.0 * ((std::abs(dlft) < std::abs(drgt)) ?
                                        std::abs(dlft) : std::abs(drgt));
                    Real dlim = (dlft*drgt >= 0.0) ? slop : 0.0;
                    Real d = dsgn*((dlim < std::abs(dcen)) ?
                                    dlim : std::abs(dcen));
                    return d;
                };

                Real dlft = 0.5*(q(i,j,k,QPRES)-q(i-1,j,k,QPRES))/q(i,j,k,QCS) - 0.5*q(i,j,k,QRHO)*(q(i,j,k,QU) - q(i-1,j,k,QU));
                Real drgt = 0.5*(q(i+1,j,k,QPRES)-q(i,j,k,QPRES))/q(i,j,k,QCS) - 0.5*q(i,j,k,QRHO)*(q(i+1,j,k,QU) - q(i,j,k,QU));
                Real d0 = f(dlft, drgt);

                Real cs2 = q(i,j,k,QCS)*q(i,j,k,QCS);
                dlft = (q(i,j,k,QRHO)-q(i-1,j,k,QRHO)) - (q(i,j,k,QPRES) - q(i-1,j,k,QPRES))/cs2;
                drgt = (q(i+1,j,k,QRHO)-q(i,j,k,QRHO)) - (q(i+1,j,k,QPRES) - q(i,j,k,QPRES))/cs2;
                Real d1 = f(dlft, drgt);

                dlft = 0.5*(q(i,j,k,QPRES)-q(i-1,j,k,QPRES))/q(i,j,k,QCS) + 0.5*q(i,j,k,QRHO)*(q(i,j,k,QU) - q(i-1,j,k,QU));
                drgt = 0.5*(q(i+1,j,k,QPRES)-q(i,j,k,QPRES))/q(i,j,k,QCS) + 0.5*q(i,j,k,QRHO)*(q(i+1,j,k,QU) - q(i,j,k,QU));
                Real d2 = f(dlft, drgt);

                dlft = q(i,j,k,QV) - q(i-1,j,k,QV);
                drgt = q(i+1,j,k,QV) - q(i,j,k,QV);
                Real d3 = f(dlft, drgt);

                dlft = q(i,j,k,QW) - q(i-1,j,k,QW);
                drgt = q(i+1,j,k,QW) - q(i,j,k,QW);
                Real d4 = f(dlft, drgt);

                dq(i,j,k,0) = d0;
                dq(i,j,k,1) = d1;
                dq(i,j,k,2) = d2;
                dq(i,j,k,3) = d3;
                dq(i,j,k,4) = d4;
            }
        }
    }
}

AMREX_GPU_DEVICE
void cns_slope_y (Box const& bx, FArrayBox& dqfab, FArrayBox const& qfab)
{
    const auto len = length(bx);
    const auto lo  = lbound(bx);
    const auto dq = dqfab.view(lo);
    const auto  q =  qfab.view(lo);

    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = 0; i < len.x; ++i) {
                auto f = [] (Real dlft, Real drgt) -> Real {
                    Real dcen = 0.5*(dlft+drgt);
                    Real dsgn = std::copysign(1.0, dcen);
                    Real slop = 2.0 * ((std::abs(dlft) < std::abs(drgt)) ?
                                        std::abs(dlft) : std::abs(drgt));
                    Real dlim = (dlft*drgt >= 0.0) ? slop : 0.0;
                    Real d = dsgn*((dlim < std::abs(dcen)) ?
                                    dlim : std::abs(dcen));
                    return d;
                };

                Real dlft = 0.5*(q(i,j,k,QPRES)-q(i,j-1,k,QPRES))/q(i,j,k,QCS) - 0.5*q(i,j,k,QRHO)*(q(i,j,k,QV) - q(i,j-1,k,QV));
                Real drgt = 0.5*(q(i,j+1,k,QPRES)-q(i,j,k,QPRES))/q(i,j,k,QCS) - 0.5*q(i,j,k,QRHO)*(q(i,j+1,k,QV) - q(i,j,k,QV));
                Real d0 = f(dlft, drgt);

                Real cs2 = q(i,j,k,QCS)*q(i,j,k,QCS);
                dlft = (q(i,j,k,QRHO)-q(i,j-1,k,QRHO)) - (q(i,j,k,QPRES) - q(i,j-1,k,QPRES))/cs2;
                drgt = (q(i,j+1,k,QRHO)-q(i,j,k,QRHO)) - (q(i,j+1,k,QPRES) - q(i,j,k,QPRES))/cs2;
                Real d1 = f(dlft, drgt);

                dlft = 0.5*(q(i,j,k,QPRES)-q(i,j-1,k,QPRES))/q(i,j,k,QCS) + 0.5*q(i,j,k,QRHO)*(q(i,j,k,QV) - q(i,j-1,k,QV));
                drgt = 0.5*(q(i,j+1,k,QPRES)-q(i,j,k,QPRES))/q(i,j,k,QCS) + 0.5*q(i,j,k,QRHO)*(q(i,j+1,k,QV) - q(i,j,k,QV));
                Real d2 = f(dlft, drgt);

                dlft = q(i,j,k,QU) - q(i,j-1,k,QU);
                drgt = q(i,j+1,k,QU) - q(i,j,k,QU);
                Real d3 = f(dlft, drgt);

                dlft = q(i,j,k,QW) - q(i,j-1,k,QW);
                drgt = q(i,j+1,k,QW) - q(i,j,k,QW);
                Real d4 = f(dlft, drgt);

                dq(i,j,k,0) = d0;
                dq(i,j,k,1) = d1;
                dq(i,j,k,2) = d2;
                dq(i,j,k,3) = d3;
                dq(i,j,k,4) = d4;
            }
        }
    }
}

AMREX_GPU_DEVICE
void cns_slope_z (Box const& bx, FArrayBox& dqfab, FArrayBox const& qfab)
{
    const auto len = length(bx);
    const auto lo  = lbound(bx);
    const auto dq = dqfab.view(lo);
    const auto  q =  qfab.view(lo);

    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = 0; i < len.x; ++i) {
                auto f = [] (Real dlft, Real drgt) -> Real {
                    Real dcen = 0.5*(dlft+drgt);
                    Real dsgn = std::copysign(1.0, dcen);
                    Real slop = 2.0 * ((std::abs(dlft) < std::abs(drgt)) ?
                                        std::abs(dlft) : std::abs(drgt));
                    Real dlim = (dlft*drgt >= 0.0) ? slop : 0.0;
                    Real d = dsgn*((dlim < std::abs(dcen)) ?
                                    dlim : std::abs(dcen));
                    return d;
                };

                Real dlft = 0.5*(q(i,j,k,QPRES)-q(i,j,k-1,QPRES))/q(i,j,k,QCS) - 0.5*q(i,j,k,QRHO)*(q(i,j,k,QW) - q(i,j,k-1,QW));
                Real drgt = 0.5*(q(i,j,k+1,QPRES)-q(i,j,k,QPRES))/q(i,j,k,QCS) - 0.5*q(i,j,k,QRHO)*(q(i,j,k+1,QW) - q(i,j,k,QW));
                Real d0 = f(dlft, drgt);

                Real cs2 = q(i,j,k,QCS)*q(i,j,k,QCS);
                dlft = (q(i,j,k,QRHO)-q(i,j,k-1,QRHO)) - (q(i,j,k,QPRES) - q(i,j,k-1,QPRES))/cs2;
                drgt = (q(i,j,k+1,QRHO)-q(i,j,k,QRHO)) - (q(i,j,k+1,QPRES) - q(i,j,k,QPRES))/cs2;
                Real d1 = f(dlft, drgt);

                dlft = 0.5*(q(i,j,k,QPRES)-q(i,j,k-1,QPRES))/q(i,j,k,QCS) + 0.5*q(i,j,k,QRHO)*(q(i,j,k,QW) - q(i,j,k-1,QW));
                drgt = 0.5*(q(i,j,k+1,QPRES)-q(i,j,k,QPRES))/q(i,j,k,QCS) + 0.5*q(i,j,k,QRHO)*(q(i,j,k+1,QW) - q(i,j,k,QW));
                Real d2 = f(dlft, drgt);

                dlft = q(i,j,k,QU) - q(i,j,k-1,QU);
                drgt = q(i,j,k+1,QU) - q(i,j,k,QU);
                Real d3 = f(dlft, drgt);

                dlft = q(i,j,k,QV) - q(i,j,k-1,QV);
                drgt = q(i,j,k+1,QV) - q(i,j,k,QV);
                Real d4 = f(dlft, drgt);

                dq(i,j,k,0) = d0;
                dq(i,j,k,1) = d1;
                dq(i,j,k,2) = d2;
                dq(i,j,k,3) = d3;
                dq(i,j,k,4) = d4;
            }
        }
    }
}

