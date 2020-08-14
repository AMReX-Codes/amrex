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
             amrex::Array4<amrex::Real> const& q,
             Parm const& parm) noexcept
{
    using amrex::Real;

    Real rho = amrex::max(u(i,j,k,URHO),parm.smallr);
    Real rhoinv = Real(1.0)/rho;
    Real ux = u(i,j,k,UMX)*rhoinv;
    Real uy = u(i,j,k,UMY)*rhoinv;
    Real uz = u(i,j,k,UMZ)*rhoinv;
    Real kineng = Real(0.5)*rho*(ux*ux+uy*uy+uz*uz);
    Real ei = u(i,j,k,UEDEN) - kineng;
    if (ei <= Real(0.0)) ei = u(i,j,k,UEINT);
    Real p = amrex::max((parm.eos_gamma-Real(1.0))*ei,parm.smallp);
    ei *= rhoinv;

    q(i,j,k,QRHO) = rho;
    q(i,j,k,QU) = ux;
    q(i,j,k,QV) = uy;
    q(i,j,k,QW) = uz;
    q(i,j,k,QEINT) = ei;
    q(i,j,k,QPRES) = p;
    q(i,j,k,QCS) = std::sqrt(parm.eos_gamma*p*rhoinv);
    q(i,j,k,QTEMP) = Real(0.0);
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
    using amrex::Real;

    Real dcen = Real(0.5)*(dlft+drgt);
    Real dsgn = amrex::Math::copysign(Real(1.0), dcen);
    Real slop = Real(2.0) * amrex::min(amrex::Math::abs(dlft),amrex::Math::abs(drgt));
    Real dlim = (dlft*drgt >= Real(0.0)) ? slop : Real(0.0);
    return dsgn * amrex::min(dlim,amrex::Math::abs(dcen));
}

}

AMREX_GPU_DEVICE
inline
void
cns_slope_x (int i, int j, int k,
             amrex::Array4<amrex::Real> const& dq,
             amrex::Array4<amrex::Real const> const& q) noexcept
{
    using amrex::Real;

    Real dlft = Real(0.5)*(q(i,j,k,QPRES)-q(i-1,j,k,QPRES))/q(i,j,k,QCS) - Real(0.5)*q(i,j,k,QRHO)*(q(i,j,k,QU) - q(i-1,j,k,QU));
    Real drgt = Real(0.5)*(q(i+1,j,k,QPRES)-q(i,j,k,QPRES))/q(i,j,k,QCS) - Real(0.5)*q(i,j,k,QRHO)*(q(i+1,j,k,QU) - q(i,j,k,QU));
    Real d0 = limiter(dlft, drgt);

    Real cs2 = q(i,j,k,QCS)*q(i,j,k,QCS);
    dlft = (q(i,j,k,QRHO)-q(i-1,j,k,QRHO)) - (q(i,j,k,QPRES) - q(i-1,j,k,QPRES))/cs2;
    drgt = (q(i+1,j,k,QRHO)-q(i,j,k,QRHO)) - (q(i+1,j,k,QPRES) - q(i,j,k,QPRES))/cs2;
    Real d1 = limiter(dlft, drgt);

    dlft = Real(0.5)*(q(i,j,k,QPRES)-q(i-1,j,k,QPRES))/q(i,j,k,QCS) + Real(0.5)*q(i,j,k,QRHO)*(q(i,j,k,QU) - q(i-1,j,k,QU));
    drgt = Real(0.5)*(q(i+1,j,k,QPRES)-q(i,j,k,QPRES))/q(i,j,k,QCS) + Real(0.5)*q(i,j,k,QRHO)*(q(i+1,j,k,QU) - q(i,j,k,QU));
    Real d2 = limiter(dlft, drgt);

    dlft = q(i,j,k,QV) - q(i-1,j,k,QV);
    drgt = q(i+1,j,k,QV) - q(i,j,k,QV);
    Real d3 = limiter(dlft, drgt);

    dlft = q(i,j,k,QW) - q(i-1,j,k,QW);
    drgt = q(i+1,j,k,QW) - q(i,j,k,QW);
    Real d4 = limiter(dlft, drgt);

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
    using amrex::Real;

    Real dlft = Real(0.5)*(q(i,j,k,QPRES)-q(i,j-1,k,QPRES))/q(i,j,k,QCS) - Real(0.5)*q(i,j,k,QRHO)*(q(i,j,k,QV) - q(i,j-1,k,QV));
    Real drgt = Real(0.5)*(q(i,j+1,k,QPRES)-q(i,j,k,QPRES))/q(i,j,k,QCS) - Real(0.5)*q(i,j,k,QRHO)*(q(i,j+1,k,QV) - q(i,j,k,QV));
    Real d0 = limiter(dlft, drgt);

    Real cs2 = q(i,j,k,QCS)*q(i,j,k,QCS);
    dlft = (q(i,j,k,QRHO)-q(i,j-1,k,QRHO)) - (q(i,j,k,QPRES) - q(i,j-1,k,QPRES))/cs2;
    drgt = (q(i,j+1,k,QRHO)-q(i,j,k,QRHO)) - (q(i,j+1,k,QPRES) - q(i,j,k,QPRES))/cs2;
    Real d1 = limiter(dlft, drgt);

    dlft = Real(0.5)*(q(i,j,k,QPRES)-q(i,j-1,k,QPRES))/q(i,j,k,QCS) + Real(0.5)*q(i,j,k,QRHO)*(q(i,j,k,QV) - q(i,j-1,k,QV));
    drgt = Real(0.5)*(q(i,j+1,k,QPRES)-q(i,j,k,QPRES))/q(i,j,k,QCS) + Real(0.5)*q(i,j,k,QRHO)*(q(i,j+1,k,QV) - q(i,j,k,QV));
    Real d2 = limiter(dlft, drgt);

    dlft = q(i,j,k,QU) - q(i,j-1,k,QU);
    drgt = q(i,j+1,k,QU) - q(i,j,k,QU);
    Real d3 = limiter(dlft, drgt);

    dlft = q(i,j,k,QW) - q(i,j-1,k,QW);
    drgt = q(i,j+1,k,QW) - q(i,j,k,QW);
    Real d4 = limiter(dlft, drgt);

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
    using amrex::Real;

    Real dlft = Real(0.5)*(q(i,j,k,QPRES)-q(i,j,k-1,QPRES))/q(i,j,k,QCS) - Real(0.5)*q(i,j,k,QRHO)*(q(i,j,k,QW) - q(i,j,k-1,QW));
    Real drgt = Real(0.5)*(q(i,j,k+1,QPRES)-q(i,j,k,QPRES))/q(i,j,k,QCS) - Real(0.5)*q(i,j,k,QRHO)*(q(i,j,k+1,QW) - q(i,j,k,QW));
    Real d0 = limiter(dlft, drgt);

    Real cs2 = q(i,j,k,QCS)*q(i,j,k,QCS);
    dlft = (q(i,j,k,QRHO)-q(i,j,k-1,QRHO)) - (q(i,j,k,QPRES) - q(i,j,k-1,QPRES))/cs2;
    drgt = (q(i,j,k+1,QRHO)-q(i,j,k,QRHO)) - (q(i,j,k+1,QPRES) - q(i,j,k,QPRES))/cs2;
    Real d1 = limiter(dlft, drgt);

    dlft = Real(0.5)*(q(i,j,k,QPRES)-q(i,j,k-1,QPRES))/q(i,j,k,QCS) + Real(0.5)*q(i,j,k,QRHO)*(q(i,j,k,QW) - q(i,j,k-1,QW));
    drgt = Real(0.5)*(q(i,j,k+1,QPRES)-q(i,j,k,QPRES))/q(i,j,k,QCS) + Real(0.5)*q(i,j,k,QRHO)*(q(i,j,k+1,QW) - q(i,j,k,QW));
    Real d2 = limiter(dlft, drgt);

    dlft = q(i,j,k,QU) - q(i,j,k-1,QU);
    drgt = q(i,j,k+1,QU) - q(i,j,k,QU);
    Real d3 = limiter(dlft, drgt);

    dlft = q(i,j,k,QV) - q(i,j,k-1,QV);
    drgt = q(i,j,k+1,QV) - q(i,j,k,QV);
    Real d4 = limiter(dlft, drgt);

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
    using amrex::Real;

    constexpr Real weakwv = Real(1.e-3);
    constexpr Real small = Real(1.e-6);

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

    Real wlsq = (Real(0.5)*(gamma-Real(1.))*(pstar+pl)+pstar)*rl;
    Real wrsq = (Real(0.5)*(gamma-Real(1.))*(pstar+pr)+pstar)*rr;

    wl = std::sqrt(wlsq);
    wr = std::sqrt(wrsq);
    Real ustarp = ul - (pstar-pl)/wl;
    Real ustarm = ur + (pstar-pr)/wr;

    pstar = (wl*pr + wr*pl - wr*wl*(ur-ul))/(wl+wr);
    pstar = amrex::max(pstar,smallp);

    Real ustar;
    for (int iter = 0; iter < 3; ++iter)
    {
        wlsq = (Real(0.5)*(gamma-Real(1.))*(pstar+pl)+pstar)*rl;
        wrsq = (Real(0.5)*(gamma-Real(1.))*(pstar+pr)+pstar)*rr;

        wl = Real(1.)/std::sqrt(wlsq);
        wr = Real(1.)/std::sqrt(wrsq);

        Real ustnm1 = ustarm;
        Real ustnp1 = ustarp;

        ustarm = ur - (pr - pstar)*wr;
        ustarp = ul + (pl - pstar)*wl;

        Real dpditer = amrex::Math::abs(pstnm1-pstar);
        Real zp = amrex::Math::abs(ustarp-ustnp1);
        if (zp-weakwv*cleft < Real(0.0) ) {
            zp = dpditer*wl;
        }
        Real zm = amrex::Math::abs(ustarm-ustnm1);
        if (zm-weakwv*cright < Real(0.0) ) {
            zm = dpditer*wr;
        }

        Real zz = zp+zm;
        Real denom = dpditer/ amrex::max(zz,ccsmall);
        pstnm1 = pstar;
        pstar = pstar - denom*(ustarm-ustarp);
        pstar = amrex::max(pstar,smallp);
        ustar = Real(0.5)*(ustarm+ustarp);
    }

    Real ro, uo, po, sgnm, utrans1, utrans2;
    if (ustar > Real(0.)) {
        ro = rl;
        uo = ul;
        po = pl;
        sgnm = Real(1.);
        utrans1 = ut1l;
        utrans2 = ut2l;
    } else if (ustar < Real(0.)) {
        ro = rr;
        uo = ur;
        po = pr;
        sgnm = Real(-1.);
        utrans1 = ut1r;
        utrans2 = ut2r;
    } else {
        uo = Real(0.5)*(ur+ul);
        po = Real(0.5)*(pr+pl);
        ro = Real(2.)*(rl*rr)/(rl+rr);
        sgnm = Real(1.);
        utrans1 = Real(0.5)*(ut1l+ut1r);
        utrans2 = Real(0.5)*(ut2l+ut2r);
    }
    Real wosq = (Real(0.5)*(gamma-Real(1.))*(pstar+po)+pstar)*ro;
    Real co = std::sqrt(gamma * po / ro);
    Real wo = std::sqrt(wosq);
    Real dpjmp = pstar-po;
    Real rstar = ro/(Real(1.)-ro*dpjmp/wosq);
    Real cstar = std::sqrt(gamma * pstar / rstar);
    Real spout = co-sgnm*uo;
    Real spin = cstar - sgnm*uo;
    if(pstar >= po) {
        spin = wo/ro-sgnm*uo;
        spout = spin;
    }
    Real ss = amrex::max(spout-spin, spout+spin);
    Real frac = Real(0.5)*(Real(1.)+(spin+spout)/amrex::max(ss,ccsmall));

    Real rgdnv, ugdnv, pgdnv;
    if (spout < Real(0.)) {
        rgdnv = ro;
        ugdnv = uo;
        pgdnv = po;
    } else if(spin >= Real(0.)) {
        rgdnv = rstar;
        ugdnv = ustar;
        pgdnv = pstar;
    } else {
        rgdnv = frac*rstar + (Real(1.) - frac)* ro;
        ugdnv = frac*ustar + (Real(1.) - frac)* uo;
        pgdnv = frac*pstar + (Real(1.) - frac)* po;
    }
    
    flxrho = rgdnv*ugdnv;
    flxu = rgdnv*ugdnv*ugdnv+pgdnv;
    flxut = rgdnv*ugdnv*utrans1;
    flxutt = rgdnv*ugdnv*utrans2;
    flxe = ugdnv*(Real(0.5)*rgdnv*(ugdnv*ugdnv+utrans1*utrans1+utrans2*utrans2) + pgdnv/(gamma -Real(1.)) + pgdnv);
}
}

AMREX_GPU_DEVICE
inline
void
cns_riemann_x (int i, int j, int k,
               amrex::Array4<amrex::Real> const& fx,
               amrex::Array4<amrex::Real const> const& dq,
               amrex::Array4<amrex::Real const> const& q,
               Parm const& parm) noexcept
{
    using amrex::Real;

    Real cspeed = q(i-1,j,k,QCS);
    Real rl = q(i-1,j,k,QRHO) + Real(0.5) * ( (dq(i-1,j,k,0)+dq(i-1,j,k,2))/cspeed + dq(i-1,j,k,1));
    rl = amrex::max(rl, parm.smallr);
    Real ul = q(i-1,j,k,QU) + Real(0.5) * ( (dq(i-1,j,k,2)-dq(i-1,j,k,0))/q(i-1,j,k,QRHO));
    Real pl = q(i-1,j,k,QPRES) + Real(0.5) *  (dq(i-1,j,k,0)+dq(i-1,j,k,2))*cspeed;
    pl = amrex::max(pl, parm.smallp);
    Real ut1l = q(i-1,j,k,QV) + Real(0.5) * dq(i-1,j,k,3);
    Real ut2l = q(i-1,j,k,QW) + Real(0.5) * dq(i-1,j,k,4);

    cspeed = q(i,j,k,QCS);
    Real rr = q(i,j,k,QRHO) - Real(0.5) * ( (dq(i,j,k,0)+dq(i,j,k,2))/cspeed + dq(i,j,k,1));
    rr = amrex::max(rr, parm.smallr);
    Real ur = q(i,j,k,QU) - Real(0.5) * ( (dq(i,j,k,2)-dq(i,j,k,0))/q(i,j,k,QRHO));
    Real pr = q(i,j,k,QPRES) - Real(0.5) * (dq(i,j,k,0)+dq(i,j,k,2))*cspeed;
    pr = amrex::max(pr, parm.smallp);
    Real ut1r = q(i,j,k,QV) - Real(0.5) * dq(i,j,k,3);
    Real ut2r = q(i,j,k,QW) - Real(0.5) * dq(i,j,k,4);

    riemann(parm.eos_gamma, parm.smallp, parm.smallr,
            rl, ul, pl, ut1l, ut2l, rr, ur, pr, ut1r, ut2r,
            fx(i,j,k,URHO), fx(i,j,k,UMX), fx(i,j,k,UMY), fx(i,j,k,UMZ), fx(i,j,k,UEDEN));
}

AMREX_GPU_DEVICE
inline
void
cns_riemann_y (int i, int j, int k,
               amrex::Array4<amrex::Real> const& fy,
               amrex::Array4<amrex::Real const> const& dq,
               amrex::Array4<amrex::Real const> const& q,
               Parm const& parm) noexcept
{
    using amrex::Real;

    Real cspeed = q(i,j-1,k,QCS);
    Real rl = q(i,j-1,k,QRHO) + Real(0.5) * ( (dq(i,j-1,k,0)+dq(i,j-1,k,2))/cspeed + dq(i,j-1,k,1));
    rl = amrex::max(rl, parm.smallr);
    Real ul = q(i,j-1,k,QV) + Real(0.5) * ( (dq(i,j-1,k,2)-dq(i,j-1,k,0))/q(i,j-1,k,QRHO));
    Real pl = q(i,j-1,k,QPRES) + Real(0.5) *  (dq(i,j-1,k,0)+dq(i,j-1,k,2))*cspeed;
    pl = amrex::max(pl, parm.smallp);
    Real ut1l = q(i,j-1,k,QU) + Real(0.5) * dq(i,j-1,k,3);
    Real ut2l = q(i,j-1,k,QW) + Real(0.5) * dq(i,j-1,k,4);

    cspeed = q(i,j,k,QCS);
    Real rr = q(i,j,k,QRHO) - Real(0.5) * ( (dq(i,j,k,0)+dq(i,j,k,2))/cspeed + dq(i,j,k,1));
    rr = amrex::max(rr, parm.smallr);
    Real ur = q(i,j,k,QV) - Real(0.5) * ( (dq(i,j,k,2)-dq(i,j,k,0))/q(i,j,k,QRHO));
    Real pr = q(i,j,k,QPRES) - Real(0.5) * (dq(i,j,k,0)+dq(i,j,k,2))*cspeed;
    pr = amrex::max(pr, parm.smallp);
    Real ut1r = q(i,j,k,QU) - Real(0.5) * dq(i,j,k,3);
    Real ut2r = q(i,j,k,QW) - Real(0.5) * dq(i,j,k,4);

    riemann(parm.eos_gamma, parm.smallp, parm.smallr,
            rl, ul, pl, ut1l, ut2l, rr, ur, pr, ut1r, ut2r,
            fy(i,j,k,URHO), fy(i,j,k,UMY), fy(i,j,k,UMX), fy(i,j,k,UMZ), fy(i,j,k,UEDEN));
}

AMREX_GPU_DEVICE
inline
void
cns_riemann_z (int i, int j, int k,
               amrex::Array4<amrex::Real> const& fz,
               amrex::Array4<amrex::Real const> const& dq,
               amrex::Array4<amrex::Real const> const& q,
               Parm const& parm) noexcept
{
    using amrex::Real;

    Real cspeed = q(i,j,k-1,QCS);
    Real rl = q(i,j,k-1,QRHO) + Real(0.5) * ( (dq(i,j,k-1,0)+dq(i,j,k-1,2))/cspeed + dq(i,j,k-1,1));
    rl = amrex::max(rl, parm.smallr);
    Real ul = q(i,j,k-1,QW) + Real(0.5) * ( (dq(i,j,k-1,2)-dq(i,j,k-1,0))/q(i,j,k-1,QRHO));
    Real pl = q(i,j,k-1,QPRES) + Real(0.5) *  (dq(i,j,k-1,0)+dq(i,j,k-1,2))*cspeed;
    pl = amrex::max(pl, parm.smallp);
    Real ut1l = q(i,j,k-1,QU) + Real(0.5) * dq(i,j,k-1,3);
    Real ut2l = q(i,j,k-1,QV) + Real(0.5) * dq(i,j,k-1,4);

    cspeed = q(i,j,k,QCS);
    Real rr = q(i,j,k,QRHO) - Real(0.5) * ( (dq(i,j,k,0)+dq(i,j,k,2))/cspeed + dq(i,j,k,1));
    rr = amrex::max(rr, parm.smallr);
    Real ur = q(i,j,k,QW) - Real(0.5) * ( (dq(i,j,k,2)-dq(i,j,k,0))/q(i,j,k,QRHO));
    Real pr = q(i,j,k,QPRES) - Real(0.5) *  (dq(i,j,k,0)+dq(i,j,k,2))*cspeed;
    pr = amrex::max(pr, parm.smallp);
    Real ut1r = q(i,j,k,QU) - Real(0.5) * dq(i,j,k,3);
    Real ut2r = q(i,j,k,QV) - Real(0.5) * dq(i,j,k,4);

    riemann(parm.eos_gamma, parm.smallp, parm.smallr,
            rl, ul, pl, ut1l, ut2l, rr, ur, pr, ut1r, ut2r,
            fz(i,j,k,URHO), fz(i,j,k,UMZ), fz(i,j,k,UMX), fz(i,j,k,UMY), fz(i,j,k,UEDEN));
}

#endif
