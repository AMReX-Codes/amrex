#ifndef CNS_HYDRO_K_H_
#define CNS_HYDRO_K_H_

#include "CNS_index_macros.H"
#include "CNS_parm.H"
#include <AMReX_FArrayBox.H>
#include <cmath>

AMREX_GPU_DEVICE
inline
void
cns_ctoprim (int i, int j, int k,
             amrex::Array4<amrex::Real const> const& u,
             amrex::Array4<amrex::Real> const& q) noexcept
{
    amrex::Real rho = amrex::max(u(i,j,k,URHO),Parm::smallr);
    amrex::Real rhoinv = 1.0/rho;
    amrex::Real ux = u(i,j,k,UMX)*rhoinv;
    amrex::Real uy = u(i,j,k,UMY)*rhoinv;
    amrex::Real uz = u(i,j,k,UMZ)*rhoinv;
    amrex::Real kineng = 0.5*rho*(ux*ux+uy*uy+uz*uz);
    amrex::Real ei = u(i,j,k,UEDEN) - kineng;
    if (ei <= 0.0) ei = u(i,j,k,UEINT);
    amrex::Real p = amrex::max((Parm::eos_gamma-1.0)*ei,Parm::smallp);
    ei *= rhoinv;

    q(i,j,k,QRHO) = rho;
    q(i,j,k,QU) = ux;
    q(i,j,k,QV) = uy;
    q(i,j,k,QW) = uz;
    q(i,j,k,QEINT) = ei;
    q(i,j,k,QPRES) = p;
    q(i,j,k,QCS) = std::sqrt(Parm::eos_gamma*p*rhoinv);
    q(i,j,k,QTEMP) = 0.0;
}

AMREX_GPU_DEVICE
inline
void
cns_flux_to_dudt (int i, int j, int k, int n,
                  amrex::Array4<amrex::Real> const& dudt,
                  amrex::Array4<amrex::Real const> const& fx,
                  amrex::Array4<amrex::Real const> const& fy,
                  amrex::Array4<amrex::Real const> const& fz,
                  amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    dudt(i,j,k,n) = dxinv[0] * (fx(i,j,k,n) - fx(i+1,j,k,n))
        +           dxinv[1] * (fy(i,j,k,n) - fy(i,j+1,k,n))
        +           dxinv[2] * (fz(i,j,k,n) - fz(i,j,k+1,n));
}

namespace {

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real limiter (amrex::Real dlft, amrex::Real drgt) noexcept
{
    amrex::Real dcen = 0.5*(dlft+drgt);
    amrex::Real dsgn = std::copysign(1.0, dcen);
    amrex::Real slop = 2.0 * amrex::min(std::abs(dlft),std::abs(drgt));
    amrex::Real dlim = (dlft*drgt >= 0.0) ? slop : 0.0;
    return dsgn * amrex::min(dlim,std::abs(dcen));
}

}

AMREX_GPU_DEVICE
inline
void
cns_slope_x (int i, int j, int k,
             amrex::Array4<amrex::Real> const& dq,
             amrex::Array4<amrex::Real const> const& q) noexcept
{
    amrex::Real dlft = 0.5*(q(i,j,k,QPRES)-q(i-1,j,k,QPRES))/q(i,j,k,QCS) - 0.5*q(i,j,k,QRHO)*(q(i,j,k,QU) - q(i-1,j,k,QU));
    amrex::Real drgt = 0.5*(q(i+1,j,k,QPRES)-q(i,j,k,QPRES))/q(i,j,k,QCS) - 0.5*q(i,j,k,QRHO)*(q(i+1,j,k,QU) - q(i,j,k,QU));
    amrex::Real d0 = limiter(dlft, drgt);

    amrex::Real cs2 = q(i,j,k,QCS)*q(i,j,k,QCS);
    dlft = (q(i,j,k,QRHO)-q(i-1,j,k,QRHO)) - (q(i,j,k,QPRES) - q(i-1,j,k,QPRES))/cs2;
    drgt = (q(i+1,j,k,QRHO)-q(i,j,k,QRHO)) - (q(i+1,j,k,QPRES) - q(i,j,k,QPRES))/cs2;
    amrex::Real d1 = limiter(dlft, drgt);

    dlft = 0.5*(q(i,j,k,QPRES)-q(i-1,j,k,QPRES))/q(i,j,k,QCS) + 0.5*q(i,j,k,QRHO)*(q(i,j,k,QU) - q(i-1,j,k,QU));
    drgt = 0.5*(q(i+1,j,k,QPRES)-q(i,j,k,QPRES))/q(i,j,k,QCS) + 0.5*q(i,j,k,QRHO)*(q(i+1,j,k,QU) - q(i,j,k,QU));
    amrex::Real d2 = limiter(dlft, drgt);

    dlft = q(i,j,k,QV) - q(i-1,j,k,QV);
    drgt = q(i+1,j,k,QV) - q(i,j,k,QV);
    amrex::Real d3 = limiter(dlft, drgt);

    dlft = q(i,j,k,QW) - q(i-1,j,k,QW);
    drgt = q(i+1,j,k,QW) - q(i,j,k,QW);
    amrex::Real d4 = limiter(dlft, drgt);

    dq(i,j,k,0) = d0;
    dq(i,j,k,1) = d1;
    dq(i,j,k,2) = d2;
    dq(i,j,k,3) = d3;
    dq(i,j,k,4) = d4;
}

AMREX_GPU_DEVICE
inline
void
cns_slope_y (int i, int j, int k,
             amrex::Array4<amrex::Real> const& dq,
             amrex::Array4<amrex::Real const> const& q) noexcept
{
    amrex::Real dlft = 0.5*(q(i,j,k,QPRES)-q(i,j-1,k,QPRES))/q(i,j,k,QCS) - 0.5*q(i,j,k,QRHO)*(q(i,j,k,QV) - q(i,j-1,k,QV));
    amrex::Real drgt = 0.5*(q(i,j+1,k,QPRES)-q(i,j,k,QPRES))/q(i,j,k,QCS) - 0.5*q(i,j,k,QRHO)*(q(i,j+1,k,QV) - q(i,j,k,QV));
    amrex::Real d0 = limiter(dlft, drgt);

    amrex::Real cs2 = q(i,j,k,QCS)*q(i,j,k,QCS);
    dlft = (q(i,j,k,QRHO)-q(i,j-1,k,QRHO)) - (q(i,j,k,QPRES) - q(i,j-1,k,QPRES))/cs2;
    drgt = (q(i,j+1,k,QRHO)-q(i,j,k,QRHO)) - (q(i,j+1,k,QPRES) - q(i,j,k,QPRES))/cs2;
    amrex::Real d1 = limiter(dlft, drgt);

    dlft = 0.5*(q(i,j,k,QPRES)-q(i,j-1,k,QPRES))/q(i,j,k,QCS) + 0.5*q(i,j,k,QRHO)*(q(i,j,k,QV) - q(i,j-1,k,QV));
    drgt = 0.5*(q(i,j+1,k,QPRES)-q(i,j,k,QPRES))/q(i,j,k,QCS) + 0.5*q(i,j,k,QRHO)*(q(i,j+1,k,QV) - q(i,j,k,QV));
    amrex::Real d2 = limiter(dlft, drgt);

    dlft = q(i,j,k,QU) - q(i,j-1,k,QU);
    drgt = q(i,j+1,k,QU) - q(i,j,k,QU);
    amrex::Real d3 = limiter(dlft, drgt);

    dlft = q(i,j,k,QW) - q(i,j-1,k,QW);
    drgt = q(i,j+1,k,QW) - q(i,j,k,QW);
    amrex::Real d4 = limiter(dlft, drgt);

    dq(i,j,k,0) = d0;
    dq(i,j,k,1) = d1;
    dq(i,j,k,2) = d2;
    dq(i,j,k,3) = d3;
    dq(i,j,k,4) = d4;
}

AMREX_GPU_DEVICE
inline
void
cns_slope_z (int i, int j, int k,
             amrex::Array4<amrex::Real> const& dq,
             amrex::Array4<amrex::Real const> const& q) noexcept
{
    amrex::Real dlft = 0.5*(q(i,j,k,QPRES)-q(i,j,k-1,QPRES))/q(i,j,k,QCS) - 0.5*q(i,j,k,QRHO)*(q(i,j,k,QW) - q(i,j,k-1,QW));
    amrex::Real drgt = 0.5*(q(i,j,k+1,QPRES)-q(i,j,k,QPRES))/q(i,j,k,QCS) - 0.5*q(i,j,k,QRHO)*(q(i,j,k+1,QW) - q(i,j,k,QW));
    amrex::Real d0 = limiter(dlft, drgt);

    amrex::Real cs2 = q(i,j,k,QCS)*q(i,j,k,QCS);
    dlft = (q(i,j,k,QRHO)-q(i,j,k-1,QRHO)) - (q(i,j,k,QPRES) - q(i,j,k-1,QPRES))/cs2;
    drgt = (q(i,j,k+1,QRHO)-q(i,j,k,QRHO)) - (q(i,j,k+1,QPRES) - q(i,j,k,QPRES))/cs2;
    amrex::Real d1 = limiter(dlft, drgt);

    dlft = 0.5*(q(i,j,k,QPRES)-q(i,j,k-1,QPRES))/q(i,j,k,QCS) + 0.5*q(i,j,k,QRHO)*(q(i,j,k,QW) - q(i,j,k-1,QW));
    drgt = 0.5*(q(i,j,k+1,QPRES)-q(i,j,k,QPRES))/q(i,j,k,QCS) + 0.5*q(i,j,k,QRHO)*(q(i,j,k+1,QW) - q(i,j,k,QW));
    amrex::Real d2 = limiter(dlft, drgt);

    dlft = q(i,j,k,QU) - q(i,j,k-1,QU);
    drgt = q(i,j,k+1,QU) - q(i,j,k,QU);
    amrex::Real d3 = limiter(dlft, drgt);

    dlft = q(i,j,k,QV) - q(i,j,k-1,QV);
    drgt = q(i,j,k+1,QV) - q(i,j,k,QV);
    amrex::Real d4 = limiter(dlft, drgt);

    dq(i,j,k,0) = d0;
    dq(i,j,k,1) = d1;
    dq(i,j,k,2) = d2;
    dq(i,j,k,3) = d3;
    dq(i,j,k,4) = d4;
}

namespace {

AMREX_GPU_DEVICE
inline
void
riemann (const amrex::Real gamma, const amrex::Real smallp, const amrex::Real smallr,
         const amrex::Real rl, const amrex::Real ul, const amrex::Real pl,
         const amrex::Real ut1l, const amrex::Real ut2l,
         const amrex::Real rr, const amrex::Real ur, const amrex::Real pr,
         const amrex::Real ut1r, const amrex::Real ut2r,
         amrex::Real& flxrho, amrex::Real& flxu, amrex::Real& flxut,
         amrex::Real& flxutt, amrex::Real& flxe) noexcept
{
    constexpr amrex::Real weakwv = 1.e-3;
    constexpr amrex::Real small = 1.e-6;

    amrex::Real clsql = gamma*pl*rl;
    amrex::Real clsqr = gamma*pr*rr;
    amrex::Real wl = std::sqrt(clsql);
    amrex::Real wr = std::sqrt(clsqr);
    amrex::Real cleft = wl/rl;
    amrex::Real cright = wr/rr;
    amrex::Real ccsmall = small*(cleft+cright);

    amrex::Real pstar = (wl*pr + wr*pl - wr*wl*(ur-ul))/(wl+wr);
    pstar = amrex::max(pstar,smallp);
    amrex::Real pstnm1 = pstar;

    amrex::Real wlsq = (.5*(gamma-1.)*(pstar+pl)+pstar)*rl;
    amrex::Real wrsq = (.5*(gamma-1.)*(pstar+pr)+pstar)*rr;

    wl = std::sqrt(wlsq);
    wr = std::sqrt(wrsq);
    amrex::Real ustarp = ul - (pstar-pl)/wl;
    amrex::Real ustarm = ur + (pstar-pr)/wr;

    pstar = (wl*pr + wr*pl - wr*wl*(ur-ul))/(wl+wr);
    pstar = amrex::max(pstar,smallp);

    amrex::Real ustar;
    for (int iter = 0; iter < 3; ++iter)
    {
        wlsq = (.5*(gamma-1.)*(pstar+pl)+pstar)*rl;
        wrsq = (.5*(gamma-1.)*(pstar+pr)+pstar)*rr;

        wl = 1./sqrt(wlsq);
        wr = 1./sqrt(wrsq);

        amrex::Real ustnm1 = ustarm;
        amrex::Real ustnp1 = ustarp;

        ustarm = ur - (pr - pstar)*wr;
        ustarp = ul + (pl - pstar)*wl;

        amrex::Real dpditer = std::abs(pstnm1-pstar);
        amrex::Real zp = std::abs(ustarp-ustnp1);
        if (zp-weakwv*cleft < 0.0 ) {
            zp = dpditer*wl;
        }
        amrex::Real zm = std::abs(ustarm-ustnm1);
        if (zm-weakwv*cright < 0.0 ) {
            zm = dpditer*wr;
        }

        amrex::Real zz = zp+zm;
        amrex::Real denom = dpditer/ amrex::max(zz,ccsmall);
        pstnm1 = pstar;
        pstar = pstar - denom*(ustarm-ustarp);
        pstar = amrex::max(pstar,smallp);
        ustar = 0.5*(ustarm+ustarp);
    }

    amrex::Real ro, uo, po, sgnm, utrans1, utrans2;
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
    amrex::Real wosq = (.5*(gamma-1.)*(pstar+po)+pstar)*ro;
    amrex::Real co = std::sqrt(gamma * po / ro);
    amrex::Real wo = std::sqrt(wosq);
    amrex::Real dpjmp = pstar-po;
    amrex::Real rstar = ro/(1.-ro*dpjmp/wosq);
    amrex::Real cstar = sqrt(gamma * pstar / rstar);
    amrex::Real spout = co-sgnm*uo;
    amrex::Real spin = cstar - sgnm*uo;
    if(pstar >= po) {
        spin = wo/ro-sgnm*uo;
        spout = spin;
    }
    amrex::Real ss = amrex::max(spout-spin, spout+spin);
    amrex::Real frac = 0.5*(1.+(spin+spout)/amrex::max(ss,ccsmall));

    amrex::Real rgdnv, ugdnv, pgdnv;
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
inline
void
cns_riemann_x (int i, int j, int k,
               amrex::Array4<amrex::Real> const& fx,
               amrex::Array4<amrex::Real const> const& dq,
               amrex::Array4<amrex::Real const> const& q) noexcept
{
    amrex::Real cspeed = q(i-1,j,k,QCS);
    amrex::Real rl = q(i-1,j,k,QRHO) + 0.5 * ( (dq(i-1,j,k,0)+dq(i-1,j,k,2))/cspeed + dq(i-1,j,k,1));
    rl = amrex::max(rl, Parm::smallr);
    amrex::Real ul = q(i-1,j,k,QU) + 0.5 * ( (dq(i-1,j,k,2)-dq(i-1,j,k,0))/q(i-1,j,k,QRHO));
    amrex::Real pl = q(i-1,j,k,QPRES) + 0.5 *  (dq(i-1,j,k,0)+dq(i-1,j,k,2))*cspeed;
    pl = amrex::max(pl, Parm::smallp);
    amrex::Real ut1l = q(i-1,j,k,QV) + 0.5 * dq(i-1,j,k,3);
    amrex::Real ut2l = q(i-1,j,k,QW) + 0.5 * dq(i-1,j,k,4);

    cspeed = q(i,j,k,QCS);
    amrex::Real rr = q(i,j,k,QRHO) - 0.5 * ( (dq(i,j,k,0)+dq(i,j,k,2))/cspeed + dq(i,j,k,1));
    rr = amrex::max(rr, Parm::smallr);
    amrex::Real ur = q(i,j,k,QU) - 0.5 * ( (dq(i,j,k,2)-dq(i,j,k,0))/q(i,j,k,QRHO));
    amrex::Real pr = q(i,j,k,QPRES) - 0.5 * (dq(i,j,k,0)+dq(i,j,k,2))*cspeed;
    pr = amrex::max(pr, Parm::smallp);
    amrex::Real ut1r = q(i,j,k,QV) - 0.5 * dq(i,j,k,3);
    amrex::Real ut2r = q(i,j,k,QW) - 0.5 * dq(i,j,k,4);

    riemann(Parm::eos_gamma, Parm::smallp, Parm::smallr,
            rl, ul, pl, ut1l, ut2l, rr, ur, pr, ut1r, ut2r,
            fx(i,j,k,URHO), fx(i,j,k,UMX), fx(i,j,k,UMY), fx(i,j,k,UMZ), fx(i,j,k,UEDEN));
}

AMREX_GPU_DEVICE
inline
void
cns_riemann_y (int i, int j, int k,
               amrex::Array4<amrex::Real> const& fy,
               amrex::Array4<amrex::Real const> const& dq,
               amrex::Array4<amrex::Real const> const& q) noexcept
{
    amrex::Real cspeed = q(i,j-1,k,QCS);
    amrex::Real rl = q(i,j-1,k,QRHO) + 0.5 * ( (dq(i,j-1,k,0)+dq(i,j-1,k,2))/cspeed + dq(i,j-1,k,1));
    rl = amrex::max(rl, Parm::smallr);
    amrex::Real ul = q(i,j-1,k,QV) + 0.5 * ( (dq(i,j-1,k,2)-dq(i,j-1,k,0))/q(i,j-1,k,QRHO));
    amrex::Real pl = q(i,j-1,k,QPRES) + 0.5 *  (dq(i,j-1,k,0)+dq(i,j-1,k,2))*cspeed;
    pl = amrex::max(pl, Parm::smallp);
    amrex::Real ut1l = q(i,j-1,k,QU) + 0.5 * dq(i,j-1,k,3);
    amrex::Real ut2l = q(i,j-1,k,QW) + 0.5 * dq(i,j-1,k,4);

    cspeed = q(i,j,k,QCS);
    amrex::Real rr = q(i,j,k,QRHO) - 0.5 * ( (dq(i,j,k,0)+dq(i,j,k,2))/cspeed + dq(i,j,k,1));
    rr = amrex::max(rr, Parm::smallr);
    amrex::Real ur = q(i,j,k,QV) - 0.5 * ( (dq(i,j,k,2)-dq(i,j,k,0))/q(i,j,k,QRHO));
    amrex::Real pr = q(i,j,k,QPRES) - 0.5 * (dq(i,j,k,0)+dq(i,j,k,2))*cspeed;
    pr = amrex::max(pr, Parm::smallp);
    amrex::Real ut1r = q(i,j,k,QU) - 0.5 * dq(i,j,k,3);
    amrex::Real ut2r = q(i,j,k,QW) - 0.5 * dq(i,j,k,4);

    riemann(Parm::eos_gamma, Parm::smallp, Parm::smallr,
            rl, ul, pl, ut1l, ut2l, rr, ur, pr, ut1r, ut2r,
            fy(i,j,k,URHO), fy(i,j,k,UMY), fy(i,j,k,UMX), fy(i,j,k,UMZ), fy(i,j,k,UEDEN));
}

AMREX_GPU_DEVICE
inline
void
cns_riemann_z (int i, int j, int k,
               amrex::Array4<amrex::Real> const& fz,
               amrex::Array4<amrex::Real const> const& dq,
               amrex::Array4<amrex::Real const> const& q) noexcept
{
    amrex::Real cspeed = q(i,j,k-1,QCS);
    amrex::Real rl = q(i,j,k-1,QRHO) + 0.5 * ( (dq(i,j,k-1,0)+dq(i,j,k-1,2))/cspeed + dq(i,j,k-1,1));
    rl = amrex::max(rl, Parm::smallr);
    amrex::Real ul = q(i,j,k-1,QW) + 0.5 * ( (dq(i,j,k-1,2)-dq(i,j,k-1,0))/q(i,j,k-1,QRHO));
    amrex::Real pl = q(i,j,k-1,QPRES) + 0.5 *  (dq(i,j,k-1,0)+dq(i,j,k-1,2))*cspeed;
    pl = amrex::max(pl, Parm::smallp);
    amrex::Real ut1l = q(i,j,k-1,QU) + 0.5 * dq(i,j,k-1,3);
    amrex::Real ut2l = q(i,j,k-1,QV) + 0.5 * dq(i,j,k-1,4);

    cspeed = q(i,j,k,QCS);
    amrex::Real rr = q(i,j,k,QRHO) - 0.5 * ( (dq(i,j,k,0)+dq(i,j,k,2))/cspeed + dq(i,j,k,1));
    rr = amrex::max(rr, Parm::smallr);
    amrex::Real ur = q(i,j,k,QW) - 0.5 * ( (dq(i,j,k,2)-dq(i,j,k,0))/q(i,j,k,QRHO));
    amrex::Real pr = q(i,j,k,QPRES) - 0.5 *  (dq(i,j,k,0)+dq(i,j,k,2))*cspeed;
    pr = amrex::max(pr, Parm::smallp);
    amrex::Real ut1r = q(i,j,k,QU) - 0.5 * dq(i,j,k,3);
    amrex::Real ut2r = q(i,j,k,QV) - 0.5 * dq(i,j,k,4);

    riemann(Parm::eos_gamma, Parm::smallp, Parm::smallr,
            rl, ul, pl, ut1l, ut2l, rr, ur, pr, ut1r, ut2r,
            fz(i,j,k,URHO), fz(i,j,k,UMZ), fz(i,j,k,UMX), fz(i,j,k,UMY), fz(i,j,k,UEDEN));
}

#endif
