#ifndef AMREX_MF_INTERP_1D_C_H_
#define AMREX_MF_INTERP_1D_C_H_

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mf_cell_cons_lin_interp_limit_minmax_llslope (int i, int, int, Array4<Real> const& slope,
                                      Array4<Real const> const& u, int scomp, int ncomp,
                                      Box const& domain, IntVect const& ratio, BCRec const* bc) noexcept
{
    Real sfx = Real(1.0);

    for (int ns = 0; ns < ncomp; ++ns) {
        int nu = ns + scomp;

        // x-direction
        Real dc = mf_compute_slopes_x(i, 0, 0, u, nu, domain, bc[ns]);
        Real df = Real(2.0) * (u(i+1,0,0,nu) - u(i  ,0,0,nu));
        Real db = Real(2.0) * (u(i  ,0,0,nu) - u(i-1,0,0,nu));
        Real sx = (df*db >= Real(0.0)) ? amrex::min(std::abs(df),std::abs(db)) : Real(0.);
        sx = std::copysign(Real(1.),dc)*amrex::min(sx,std::abs(dc));
        slope(i,0,0,ns) = dc;

        // additional limiting is unnecessary in 1D

        if (dc != Real(0.0)) {
            sfx = amrex::min(sfx, sx / dc);
        }
    }

    for (int ns = 0; ns < ncomp; ++ns) {
        slope(i,0,0,ns) *= sfx;
    }

    amrex::ignore_unused(ratio);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mf_cell_cons_lin_interp_llslope (int i, int, int, Array4<Real> const& slope,
                                      Array4<Real const> const& u, int scomp, int ncomp,
                                      Box const& domain, BCRec const* bc) noexcept
{
    Real sfx = Real(1.0);

    for (int ns = 0; ns < ncomp; ++ns) {
        int nu = ns + scomp;

        // x-direction
        Real dc = mf_compute_slopes_x(i, 0, 0, u, nu, domain, bc[ns]);
        Real df = Real(2.0) * (u(i+1,0,0,nu) - u(i  ,0,0,nu));
        Real db = Real(2.0) * (u(i  ,0,0,nu) - u(i-1,0,0,nu));
        Real sx = (df*db >= Real(0.0)) ? amrex::min(std::abs(df),std::abs(db)) : Real(0.);
        sx = std::copysign(Real(1.),dc)*amrex::min(sx,std::abs(dc));
        if (dc != Real(0.0)) {
            sfx = amrex::min(sfx, sx / dc);
        }
        slope(i,0,0,ns) = dc;
    }

    for (int ns = 0; ns < ncomp; ++ns) {
        slope(i,0,0,ns) *= sfx;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mf_cell_cons_lin_interp_mcslope (int i, int /*j*/, int /*k*/, int ns,
                                      Array4<Real> const& slope,
                                      Array4<Real const> const& u, int scomp, int /*ncomp*/,
                                      Box const& domain, IntVect const& ratio,
                                      BCRec const* bc) noexcept
{
    int nu = ns + scomp;

    // x-direction
    Real dc = mf_compute_slopes_x(i, 0, 0, u, nu, domain, bc[ns]);
    Real df = Real(2.0) * (u(i+1,0,0,nu) - u(i  ,0,0,nu));
    Real db = Real(2.0) * (u(i  ,0,0,nu) - u(i-1,0,0,nu));
    Real sx = (df*db >= Real(0.0)) ? amrex::min(std::abs(df),std::abs(db)) : Real(0.);
    sx = std::copysign(Real(1.),dc)*amrex::min(sx,std::abs(dc));

    Real alpha = Real(1.0);
    if (sx != Real(0.0)) {
        Real dumax = std::abs(sx) * Real(ratio[0]-1)/Real(2*ratio[0]);
        Real umax = u(i,0,0,nu);
        Real umin = u(i,0,0,nu);
        for (int ioff = -1; ioff <= 1; ++ioff) {
            umin = amrex::min(umin, u(i+ioff,0,0,nu));
            umax = amrex::max(umax, u(i+ioff,0,0,nu));
        }
        if (dumax * alpha > (umax - u(i,0,0,nu))) {
            alpha = (umax - u(i,0,0,nu)) / dumax;
        }
        if (dumax * alpha > (u(i,0,0,nu) - umin)) {
            alpha = (u(i,0,0,nu) - umin) / dumax;
        }
    }

    slope(i,0,0,ns) = sx * alpha;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mf_cell_cons_lin_interp (int i, int /*j*/, int /*k*/, int ns,
                              Array4<Real> const& fine, int fcomp,
                              Array4<Real const> const& slope, Array4<Real const> const& crse,
                              int ccomp, int /*ncomp*/, IntVect const& ratio) noexcept
{
    const int ic = amrex::coarsen(i, ratio[0]);
    const Real xoff = (static_cast<Real>(i - ic*ratio[0]) + Real(0.5)) / Real(ratio[0]) - Real(0.5);
    fine(i,0,0,fcomp+ns) = crse(ic,0,0,ccomp+ns)
        + xoff * slope(ic,0,0,ns);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mf_cell_cons_lin_interp_mcslope_sph (int i, int ns, Array4<Real> const& slope,
                                          Array4<Real const> const& u, int scomp, int /*ncomp*/,
                                          Box const& domain, IntVect const& ratio,
                                          BCRec const* bc, Real drf, Real rlo) noexcept
{
    int nu = ns + scomp;

    // x-direction
    Real dc = mf_compute_slopes_x(i, 0, 0, u, nu, domain, bc[ns]);
    Real df = Real(2.0) * (u(i+1,0,0,nu) - u(i  ,0,0,nu));
    Real db = Real(2.0) * (u(i  ,0,0,nu) - u(i-1,0,0,nu));
    Real sx = (df*db >= Real(0.0)) ? amrex::min(std::abs(df),std::abs(db)) : Real(0.);
    sx = std::copysign(Real(1.),dc)*amrex::min(sx,std::abs(dc));

    Real alpha = Real(1.0);
    if (sx != Real(0.0)) {
        const Real drc = drf * ratio[0];
        const Real rcm =  i    * drc + rlo;
        const Real rcp = (i+1) * drc + rlo;
        const Real vcm = rcm*rcm*rcm;
        const Real vcp = rcp*rcp*rcp;
        Real rfm =  i*ratio[0]      * drf + rlo;
        Real rfp = (i*ratio[0] + 1) * drf + rlo;
        Real vfm = rfm*rfm*rfm;
        Real vfp = rfp*rfp*rfp;
        Real xlo = Real(0.5) * ((vfm+vfp) - (vcm+vcp)) / (vcp - vcm);
        rfm = ((i+1)*ratio[0] - 1) * drf + rlo;
        rfp =  (i+1)*ratio[0]      * drf + rlo;
        vfm = rfm*rfm*rfm;
        vfp = rfp*rfp*rfp;
        Real xhi = Real(0.5) * ((vfm+vfp) - (vcm+vcp)) / (vcp - vcm);
        Real dumax =  amrex::max(sx*xlo, sx*xhi);
        Real dumin = -amrex::min(sx*xlo, sx*xhi);
        Real umax = u(i,0,0,nu);
        Real umin = u(i,0,0,nu);
        for (int ioff = -1; ioff <= 1; ++ioff) {
            umin = amrex::min(umin, u(i+ioff,0,0,nu));
            umax = amrex::max(umax, u(i+ioff,0,0,nu));
        }
        if (dumax * alpha > (umax - u(i,0,0,nu))) {
            alpha = (umax - u(i,0,0,nu)) / dumax;
        }
        if (dumin * alpha > (u(i,0,0,nu) - umin)) {
            alpha = (u(i,0,0,nu) - umin) / dumin;
        }
    }

    slope(i,0,0,ns) = sx * alpha;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mf_cell_cons_lin_interp_sph (int i, int ns, Array4<Real> const& fine, int fcomp,
                                  Array4<Real const> const& slope, Array4<Real const> const& crse,
                                  int ccomp, int /*ncomp*/, IntVect const& ratio, Real drf, Real rlo) noexcept
{
    const int ic = amrex::coarsen(i, ratio[0]);
    const Real drc =  drf * ratio[0];
    const Real rcm =  ic    * drc + rlo;
    const Real rcp = (ic+1) * drc + rlo;
    const Real rfm =  i     * drf + rlo;
    const Real rfp = (i +1) * drf + rlo;
    const Real vcm = rcm*rcm*rcm;
    const Real vcp = rcp*rcp*rcp;
    const Real vfm = rfm*rfm*rfm;
    const Real vfp = rfp*rfp*rfp;
    const Real xoff = Real(0.5) * ((vfm+vfp) - (vcm+vcp)) / (vcp - vcm);
    fine(i,0,0,fcomp+ns) = crse(ic,0,0,ccomp+ns)
        + xoff * slope(ic,0,0,ns);
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mf_cell_bilin_interp (int i, int, int, int n, Array4<T> const& fine, int fcomp,
                           Array4<T const> const& crse, int ccomp, IntVect const& ratio) noexcept
{
    int ic = amrex::coarsen(i,ratio[0]);
    int ioff = i - ic*ratio[0];
    int sx;
    Real wx;
    if (ioff*2 < ratio[0]) {
        sx = -1;
        wx = Real(ratio[0]+1+2*ioff) / Real(2*ratio[0]);
    } else {
        sx = 1;
        wx = Real(3*ratio[0]-1-2*ioff) / Real(2*ratio[0]);
    }
    fine(i,0,0,n+fcomp) =
        crse(ic   ,0,0,n+ccomp)*           wx  +
        crse(ic+sx,0,0,n+ccomp)*(Real(1.0)-wx);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mf_nodebilin_interp (int i, int, int, int n, Array4<Real> const& fine, int fcomp,
                          Array4<Real const> const& crse, int ccomp, IntVect const& ratio) noexcept
{
    int ic = amrex::coarsen(i,ratio[0]);
    int ioff = i - ic*ratio[0];
    if (ioff != 0) {
        Real rxinv = Real(1.0) / Real(ratio[0]);
        fine(i,0,0,n+fcomp) = rxinv*(static_cast<Real>(ratio[0]-ioff)*crse(ic  ,0,0,n+ccomp) +
                                     static_cast<Real>(         ioff)*crse(ic+1,0,0,n+ccomp));
    } else {
        fine(i,0,0,n+fcomp) = crse(ic,0,0,n+ccomp);
    }
}

}

#endif
