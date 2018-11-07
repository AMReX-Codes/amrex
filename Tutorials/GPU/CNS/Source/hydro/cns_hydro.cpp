
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
                Real rho = amrex::max(u(i,j,k,URHO),smallr);
                Real rhoinv = 1.0/rho;
                Real ux = u(i,j,k,UMX)*rhoinv;
                Real uy = u(i,j,k,UMY)*rhoinv;
                Real uz = u(i,j,k,UMZ)*rhoinv;
                Real kineng = 0.5*rho*(ux*ux+uy*uy+uz*uz);
                Real ei = u(i,j,k,UEDEN) - kineng;
                if (ei <= 0.0) ei = u(i,j,k,UEINT);
                Real p = amrex::max((gamma-1.0)*ei,smallp);
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
                    Real slop = 2.0 * amrex::min(std::abs(dlft),std::abs(drgt));
                    Real dlim = (dlft*drgt >= 0.0) ? slop : 0.0;
                    Real d = dsgn * amrex::min(dlim,std::abs(dcen));
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
                    Real slop = 2.0 * amrex::min(std::abs(dlft), std::abs(drgt));
                    Real dlim = (dlft*drgt >= 0.0) ? slop : 0.0;
                    Real d = dsgn * amrex::min(dlim, std::abs(dcen));
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
                    Real slop = 2.0 * amrex::min(std::abs(dlft), std::abs(drgt));
                    Real dlim = (dlft*drgt >= 0.0) ? slop : 0.0;
                    Real d = dsgn * amrex::min(dlim, std::abs(dcen));
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

namespace {

AMREX_GPU_DEVICE
void riemann (const Real gamma, const Real smallp, const Real smallr,
              const Real rl, const Real ul, const Real pl, const Real ut1l, const Real ut2l,
              const Real rr, const Real ur, const Real pr, const Real ut1r, const Real ut2r,
              Real& flxrho, Real& flxu, Real& flxut, Real& flxutt, Real& flxe)
{
    constexpr Real weakwv = 1.e-3;
    constexpr Real small = 1.e-6;

    Real clsql = gamma*pl*rl;
    Real clsqr = gamma*pr*rr;
    Real wl = std::sqrt(clsql);
    Real wr = std::sqrt(clsqr);
    Real cleft = wl/rl;
    Real cright = wr/rr;
    Real ccsmall = small*(cleft+cright);

    Real pstar = (wl*pr + wr*pl - wr*wl*(ur-ul))/(wl+wr);
    pstar = amrex::max(pstar,smallp);
    Real pstnm1 = pstar;

    Real wlsq = (.5*(gamma-1.)*(pstar+pl)+pstar)*rl;
    Real wrsq = (.5*(gamma-1.)*(pstar+pr)+pstar)*rr;

    wl = std::sqrt(wlsq);
    wr = std::sqrt(wrsq);
    Real ustarp = ul - (pstar-pl)/wl;
    Real ustarm = ur + (pstar-pr)/wr;

    pstar = (wl*pr + wr*pl - wr*wl*(ur-ul))/(wl+wr);
    pstar = amrex::max(pstar,smallp);

    Real ustar;
    for (int iter = 0; iter < 3; ++iter)
    {
        wlsq = (.5*(gamma-1.)*(pstar+pl)+pstar)*rl;
        wrsq = (.5*(gamma-1.)*(pstar+pr)+pstar)*rr;

        wl = 1./sqrt(wlsq);
        wr = 1./sqrt(wrsq);

        Real ustnm1 = ustarm;
        Real ustnp1 = ustarp;

        ustarm = ur - (pr - pstar)*wr;
        ustarp = ul + (pl - pstar)*wl;

        Real dpditer = std::abs(pstnm1-pstar);
        Real zp = std::abs(ustarp-ustnp1);
        if (zp-weakwv*cleft < 0.0 ) {
            zp = dpditer*wl;
        }
        Real zm = std::abs(ustarm-ustnm1);
        if (zm-weakwv*cright < 0.0 ) {
            zm = dpditer*wr;
        }

        Real zz = zp+zm;
        Real denom = dpditer/ amrex::max(zz,ccsmall);
        pstnm1 = pstar;
        pstar = pstar - denom*(ustarm-ustarp);
        pstar = amrex::max(pstar,smallp);
        ustar = 0.5*(ustarm+ustarp);
    }

    Real ro, uo, po, sgnm, utrans1, utrans2;
    if (ustar > 0.) {
        ro = rl;
        uo = ul;
        po = pl;
        sgnm = 1.;
        utrans1 = ut1l;
        utrans2 = ut2l;
    } else if (ustar < 0.) {
        ro = rr;
        uo = ur;
        po = pr;
        sgnm = -1.;
        utrans1 = ut1r;
        utrans2 = ut2r;
    } else {
        uo = 0.5*(ur+ul);
        po = 0.5*(pr+pl);
        ro = 2.*(rl*rr)/(rl+rr);
        sgnm = 1.;
        utrans1 = 0.5*(ut1l+ut1r);
        utrans2 = 0.5*(ut2l+ut2r);
    }
    Real wosq = (.5*(gamma-1.)*(pstar+po)+pstar)*ro;
    Real co = std::sqrt(gamma * po / ro);
    Real wo = std::sqrt(wosq);
    Real dpjmp = pstar-po;
    Real rstar = ro/(1.-ro*dpjmp/wosq);
    Real cstar = sqrt(gamma * pstar / rstar);
    Real spout = co-sgnm*uo;
    Real spin = cstar - sgnm*uo;
    if(pstar >= po) {
        spin = wo/ro-sgnm*uo;
        spout = spin;
    }
    Real ss = amrex::max(spout-spin, spout+spin);
    Real frac = 0.5*(1.+(spin+spout)/amrex::max(ss,ccsmall));

    Real rgdnv, ugdnv, pgdnv;
    if (spout < 0.) {
        rgdnv = ro;
        ugdnv = uo;
        pgdnv = po;
    } else if(spin >= 0.) {
        rgdnv = rstar;
        ugdnv = ustar;
        pgdnv = pstar;
    } else {
        rgdnv = frac*rstar + (1. - frac)* ro;
        ugdnv = frac*ustar + (1. - frac)* uo;
        pgdnv = frac*pstar + (1. - frac)* po;
    }
    
    flxrho = rgdnv*ugdnv;
    flxu = rgdnv*ugdnv*ugdnv+pgdnv;
    flxut = rgdnv*ugdnv*utrans1;
    flxutt = rgdnv*ugdnv*utrans2;
    flxe = ugdnv*(0.5*rgdnv*(ugdnv*ugdnv+utrans1*utrans1+utrans2*utrans2) + pgdnv/(gamma -1.) + pgdnv);
}

}


AMREX_GPU_DEVICE
void cns_riemann_x (Box const& bx, FArrayBox& fluxfab, FArrayBox const& dqfab, FArrayBox const& qfab)
{
    const auto len = length(bx);
    const auto lo  = lbound(bx);
    const auto fx = fluxfab.view(lo);
    const auto dq = dqfab.view(lo);
    const auto  q =  qfab.view(lo);

    const Real smallr = 1.e-19;
    const Real smallp = 1.e-10;
    const Real gamma = 1.4;

    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            for (int i = 0; i < len.x; ++i) {
                Real cspeed = q(i-1,j,k,QCS);
                Real rl = q(i-1,j,k,QRHO) + 0.5 * ( (dq(i-1,j,k,0)+dq(i-1,j,k,2))/cspeed + dq(i-1,j,k,1));
                rl = amrex::max(rl, smallr);
                Real ul = q(i-1,j,k,QU) + 0.5 * ( (dq(i-1,j,k,2)-dq(i-1,j,k,0))/q(i-1,j,k,QRHO));
                Real pl = q(i-1,j,k,QPRES) + 0.5 *  (dq(i-1,j,k,0)+dq(i-1,j,k,2))*cspeed;
                pl = amrex::max(pl, smallp);
                Real ut1l = q(i-1,j,k,QV) + 0.5 * dq(i-1,j,k,3);
                Real ut2l = q(i-1,j,k,QW) + 0.5 * dq(i-1,j,k,4);
             
                cspeed = q(i,j,k,QCS);
                Real rr = q(i,j,k,QRHO) - 0.5 * ( (dq(i,j,k,0)+dq(i,j,k,2))/cspeed + dq(i,j,k,1));
                rr = amrex::max(rr, smallr);
                Real ur = q(i,j,k,QU) - 0.5 * ( (dq(i,j,k,2)-dq(i,j,k,0))/q(i,j,k,QRHO));
                Real pr = q(i,j,k,QPRES) - 0.5 * (dq(i,j,k,0)+dq(i,j,k,2))*cspeed;
                pr = amrex::max(pr, smallp);
                Real ut1r = q(i,j,k,QV) - 0.5 * dq(i,j,k,3);
                Real ut2r = q(i,j,k,QW) - 0.5 * dq(i,j,k,4);

                riemann(gamma, smallp, smallr, rl, ul, pl, ut1l, ut2l, rr, ur, pr, ut1r, ut2r,
                        fx(i,j,k,URHO), fx(i,j,k,UMX), fx(i,j,k,UMY), fx(i,j,k,UMZ), fx(i,j,k,UEDEN));
            }
        }
    }
}

AMREX_GPU_DEVICE
void cns_riemann_y (Box const& bx, FArrayBox& fluxfab, FArrayBox const& dqfab, FArrayBox const& qfab)
{
    const auto len = length(bx);
    const auto lo  = lbound(bx);
    const auto fy = fluxfab.view(lo);
    const auto dq = dqfab.view(lo);
    const auto  q =  qfab.view(lo);

    const Real smallr = 1.e-19;
    const Real smallp = 1.e-10;
    const Real gamma = 1.4;

    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            for (int i = 0; i < len.x; ++i) {
                Real cspeed = q(i,j-1,k,QCS);
                Real rl = q(i,j-1,k,QRHO) + 0.5 * ( (dq(i,j-1,k,0)+dq(i,j-1,k,2))/cspeed + dq(i,j-1,k,1));
                rl = amrex::max(rl, smallr);
                Real ul = q(i,j-1,k,QV) + 0.5 * ( (dq(i,j-1,k,2)-dq(i,j-1,k,0))/q(i,j-1,k,QRHO));
                Real pl = q(i,j-1,k,QPRES) + 0.5 *  (dq(i,j-1,k,0)+dq(i,j-1,k,2))*cspeed;
                pl = amrex::max(pl, smallp);
                Real ut1l = q(i,j-1,k,QU) + 0.5 * dq(i,j-1,k,3);
                Real ut2l = q(i,j-1,k,QW) + 0.5 * dq(i,j-1,k,4);

                cspeed = q(i,j,k,QCS);
                Real rr = q(i,j,k,QRHO) - 0.5 * ( (dq(i,j,k,0)+dq(i,j,k,2))/cspeed + dq(i,j,k,1));
                rr = amrex::max(rr, smallr);
                Real ur = q(i,j,k,QV) - 0.5 * ( (dq(i,j,k,2)-dq(i,j,k,0))/q(i,j,k,QRHO));
                Real pr = q(i,j,k,QPRES) - 0.5 * (dq(i,j,k,0)+dq(i,j,k,2))*cspeed;
                pr = amrex::max(pr, smallp);
                Real ut1r = q(i,j,k,QU) - 0.5 * dq(i,j,k,3);
                Real ut2r = q(i,j,k,QW) - 0.5 * dq(i,j,k,4);

                riemann(gamma, smallp, smallr, rl, ul, pl, ut1l, ut2l, rr, ur, pr, ut1r, ut2r,
                        fy(i,j,k,URHO), fy(i,j,k,UMY), fy(i,j,k,UMX), fy(i,j,k,UMZ), fy(i,j,k,UEDEN));
            }
        }
    }
}

AMREX_GPU_DEVICE
void cns_riemann_z (Box const& bx, FArrayBox& fluxfab, FArrayBox const& dqfab, FArrayBox const& qfab)
{
    const auto len = length(bx);
    const auto lo  = lbound(bx);
    const auto fz = fluxfab.view(lo);
    const auto dq = dqfab.view(lo);
    const auto  q =  qfab.view(lo);

    const Real smallr = 1.e-19;
    const Real smallp = 1.e-10;
    const Real gamma = 1.4;

    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            for (int i = 0; i < len.x; ++i) {
                Real cspeed = q(i,j,k-1,QCS);
                Real rl = q(i,j,k-1,QRHO) + 0.5 * ( (dq(i,j,k-1,0)+dq(i,j,k-1,2))/cspeed + dq(i,j,k-1,1));
                rl = amrex::max(rl, smallr);
                Real ul = q(i,j,k-1,QW) + 0.5 * ( (dq(i,j,k-1,2)-dq(i,j,k-1,0))/q(i,j,k-1,QRHO));
                Real pl = q(i,j,k-1,QPRES) + 0.5 *  (dq(i,j,k-1,0)+dq(i,j,k-1,2))*cspeed;
                pl = amrex::max(pl, smallp);
                Real ut1l = q(i,j,k-1,QU) + 0.5 * dq(i,j,k-1,3);
                Real ut2l = q(i,j,k-1,QV) + 0.5 * dq(i,j,k-1,4);
                
                cspeed = q(i,j,k,QCS);
                Real rr = q(i,j,k,QRHO) - 0.5 * ( (dq(i,j,k,0)+dq(i,j,k,2))/cspeed + dq(i,j,k,1));
                rr = amrex::max(rr, smallr);
                Real ur = q(i,j,k,QW) - 0.5 * ( (dq(i,j,k,2)-dq(i,j,k,0))/q(i,j,k,QRHO));
                Real pr = q(i,j,k,QPRES) - 0.5 *  (dq(i,j,k,0)+dq(i,j,k,2))*cspeed;
                pr = amrex::max(pr, smallp);
                Real ut1r = q(i,j,k,QU) - 0.5 * dq(i,j,k,3);
                Real ut2r = q(i,j,k,QV) - 0.5 * dq(i,j,k,4);
                
                riemann(gamma, smallp, smallr, rl, ul, pl, ut1l, ut2l, rr, ur, pr, ut1r, ut2r,
                        fz(i,j,k,URHO), fz(i,j,k,UMZ), fz(i,j,k,UMX), fz(i,j,k,UMY), fz(i,j,k,UEDEN));
            }
        }
    }
}
