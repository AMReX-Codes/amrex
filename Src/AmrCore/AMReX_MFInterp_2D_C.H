#ifndef AMREX_MF_INTERP_2D_C_H_
#define AMREX_MF_INTERP_2D_C_H_

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mf_cell_cons_lin_interp_llslope (int i, int j, int, Array4<Real> const& slope,
                                      Array4<Real const> const& u, int scomp, int ncomp,
                                      Box const& domain, BCRec const* bc) noexcept
{
    Real sfx = Real(1.0);
    Real sfy = Real(1.0);

    for (int ns = 0; ns < ncomp; ++ns) {
        int nu = ns + scomp;

        // x-direction
        Real dc = mf_compute_slopes_x(i, j, 0, u, nu, domain, bc[ns]);
        Real df = Real(2.0) * (u(i+1,j,0,nu) - u(i  ,j,0,nu));
        Real db = Real(2.0) * (u(i  ,j,0,nu) - u(i-1,j,0,nu));
        Real sx = (df*db >= Real(0.0)) ? amrex::min(amrex::Math::abs(df),amrex::Math::abs(db)) : Real(0.);
        sx = amrex::Math::copysign(Real(1.),dc)*amrex::min(sx,amrex::Math::abs(dc));
        if (dc != Real(0.0)) {
            sfx = amrex::min(sfx, sx / dc);
        }
        slope(i,j,0,ns        ) = dc;

        // y-direction
        dc = mf_compute_slopes_y(i, j, 0, u, nu, domain, bc[ns]);
        df = Real(2.0) * (u(i,j+1,0,nu) - u(i,j  ,0,nu));
        db = Real(2.0) * (u(i,j  ,0,nu) - u(i,j-1,0,nu));
        Real sy = (df*db >= Real(0.0)) ? amrex::min(amrex::Math::abs(df),amrex::Math::abs(db)) : Real(0.);
        sy = amrex::Math::copysign(Real(1.),dc)*amrex::min(sy,amrex::Math::abs(dc));
        if (dc != Real(0.0)) {
            sfy = amrex::min(sfy, sy / dc);
        }
        slope(i,j,0,ns+  ncomp) = dc;
    }

    for (int ns = 0; ns < ncomp; ++ns) {
        slope(i,j,0,ns        ) *= sfx;
        slope(i,j,0,ns+  ncomp) *= sfy;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mf_cell_cons_lin_interp_mcslope (int i, int j, int /*k*/, int ns, Array4<Real> const& slope,
                                      Array4<Real const> const& u, int scomp, int ncomp,
                                      Box const& domain, IntVect const& ratio,
                                      BCRec const* bc) noexcept
{
    int nu = ns + scomp;

    // x-direction
    Real dc = mf_compute_slopes_x(i, j, 0, u, nu, domain, bc[ns]);
    Real df = Real(2.0) * (u(i+1,j,0,nu) - u(i  ,j,0,nu));
    Real db = Real(2.0) * (u(i  ,j,0,nu) - u(i-1,j,0,nu));
    Real sx = (df*db >= Real(0.0)) ? amrex::min(amrex::Math::abs(df),amrex::Math::abs(db)) : Real(0.);
    sx = amrex::Math::copysign(Real(1.),dc)*amrex::min(sx,amrex::Math::abs(dc));

    // y-direction
    dc = mf_compute_slopes_y(i, j, 0, u, nu, domain, bc[ns]);
    df = Real(2.0) * (u(i,j+1,0,nu) - u(i,j  ,0,nu));
    db = Real(2.0) * (u(i,j  ,0,nu) - u(i,j-1,0,nu));
    Real sy = (df*db >= Real(0.0)) ? amrex::min(amrex::Math::abs(df),amrex::Math::abs(db)) : Real(0.);
    sy = amrex::Math::copysign(Real(1.),dc)*amrex::min(sy,amrex::Math::abs(dc));

    Real alpha = Real(1.0);
    if (sx != Real(0.0) || sy != Real(0.0)) {
        Real dumax = amrex::Math::abs(sx) * Real(ratio[0]-1)/Real(2*ratio[0])
            +        amrex::Math::abs(sy) * Real(ratio[1]-1)/Real(2*ratio[1]);
        Real umax = u(i,j,0,nu);
        Real umin = u(i,j,0,nu);
        for (int joff = -1; joff <= 1; ++joff) {
        for (int ioff = -1; ioff <= 1; ++ioff) {
            umin = amrex::min(umin, u(i+ioff,j+joff,0,nu));
            umax = amrex::max(umax, u(i+ioff,j+joff,0,nu));
        }}
        if (dumax * alpha > (umax - u(i,j,0,nu))) {
            alpha = (umax - u(i,j,0,nu)) / dumax;
        }
        if (dumax * alpha > (u(i,j,0,nu) - umin)) {
            alpha = (u(i,j,0,nu) - umin) / dumax;
        }
    }

    slope(i,j,0,ns        ) = sx * alpha;
    slope(i,j,0,ns+  ncomp) = sy * alpha;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mf_cell_cons_lin_interp (int i, int j, int /*k*/, int ns, Array4<Real> const& fine, int fcomp,
                              Array4<Real const> const& slope, Array4<Real const> const& crse,
                              int ccomp, int ncomp, IntVect const& ratio) noexcept
{
    const int ic = amrex::coarsen(i, ratio[0]);
    const int jc = amrex::coarsen(j, ratio[1]);
    const Real xoff = (i - ic*ratio[0] + Real(0.5)) / Real(ratio[0]) - Real(0.5);
    const Real yoff = (j - jc*ratio[1] + Real(0.5)) / Real(ratio[1]) - Real(0.5);
    fine(i,j,0,fcomp+ns) = crse(ic,jc,0,ccomp+ns)
        + xoff * slope(ic,jc,0,ns)
        + yoff * slope(ic,jc,0,ns+ncomp);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mf_cell_cons_lin_interp_mcslope_rz (int i, int j, int ns, Array4<Real> const& slope,
                                         Array4<Real const> const& u, int scomp, int ncomp,
                                         Box const& domain, IntVect const& ratio,
                                         BCRec const* bc, Real drf, Real rlo) noexcept
{
    int nu = ns + scomp;

    // x-direction
    Real dc = mf_compute_slopes_x(i, j, 0, u, nu, domain, bc[ns]);
    Real df = Real(2.0) * (u(i+1,j,0,nu) - u(i  ,j,0,nu));
    Real db = Real(2.0) * (u(i  ,j,0,nu) - u(i-1,j,0,nu));
    Real sx = (df*db >= Real(0.0)) ? amrex::min(amrex::Math::abs(df),amrex::Math::abs(db)) : Real(0.);
    sx = amrex::Math::copysign(Real(1.),dc)*amrex::min(sx,amrex::Math::abs(dc));

    // y-direction
    dc = mf_compute_slopes_y(i, j, 0, u, nu, domain, bc[ns]);
    df = Real(2.0) * (u(i,j+1,0,nu) - u(i,j  ,0,nu));
    db = Real(2.0) * (u(i,j  ,0,nu) - u(i,j-1,0,nu));
    Real sy = (df*db >= Real(0.0)) ? amrex::min(amrex::Math::abs(df),amrex::Math::abs(db)) : Real(0.);
    sy = amrex::Math::copysign(Real(1.),dc)*amrex::min(sy,amrex::Math::abs(dc));

    Real alpha = Real(1.0);
    if (sx != Real(0.0) || sy != Real(0.0)) {
        const Real drc = drf * ratio[0];
        const Real rcm =  i    * drc + rlo;
        const Real rcp = (i+1) * drc + rlo;
        const Real vcm = rcm*rcm;
        const Real vcp = rcp*rcp;
        Real rfm =  i*ratio[0]      * drf + rlo;
        Real rfp = (i*ratio[0] + 1) * drf + rlo;
        Real vfm = rfm*rfm;
        Real vfp = rfp*rfp;
        Real xlo = Real(0.5) * ((vfm+vfp) - (vcm+vcp)) / (vcp - vcm);
        rfm = ((i+1)*ratio[0] - 1) * drf + rlo;
        rfp =  (i+1)*ratio[0]      * drf + rlo;
        vfm = rfm*rfm;
        vfp = rfp*rfp;
        Real xhi = Real(0.5) * ((vfm+vfp) - (vcm+vcp)) / (vcp - vcm);
        Real dumax =  amrex::max(sx*xlo, sx*xhi)
            + amrex::Math::abs(sy) * Real(ratio[1]-1)/Real(2*ratio[1]);
        Real dumin = -amrex::min(sx*xlo, sx*xhi)
            + amrex::Math::abs(sy) * Real(ratio[1]-1)/Real(2*ratio[1]);
        Real umax = u(i,j,0,nu);
        Real umin = u(i,j,0,nu);
        for (int joff = -1; joff <= 1; ++joff) {
        for (int ioff = -1; ioff <= 1; ++ioff) {
            umin = amrex::min(umin, u(i+ioff,j+joff,0,nu));
            umax = amrex::max(umax, u(i+ioff,j+joff,0,nu));
        }}
        if (dumax * alpha > (umax - u(i,j,0,nu))) {
            alpha = (umax - u(i,j,0,nu)) / dumax;
        }
        if (dumin * alpha > (u(i,j,0,nu) - umin)) {
            alpha = (u(i,j,0,nu) - umin) / dumin;
        }
    }

    slope(i,j,0,ns        ) = sx * alpha;
    slope(i,j,0,ns+  ncomp) = sy * alpha;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mf_cell_cons_lin_interp_rz (int i, int j, int ns, Array4<Real> const& fine, int fcomp,
                                 Array4<Real const> const& slope, Array4<Real const> const& crse,
                                 int ccomp, int ncomp, IntVect const& ratio, Real drf, Real rlo) noexcept
{
    const int ic = amrex::coarsen(i, ratio[0]);
    const int jc = amrex::coarsen(j, ratio[1]);
    const Real drc =  drf * ratio[0];
    const Real rcm =  ic    * drc + rlo;
    const Real rcp = (ic+1) * drc + rlo;
    const Real rfm =  i     * drf + rlo;
    const Real rfp = (i +1) * drf + rlo;
    const Real vcm = rcm*rcm;
    const Real vcp = rcp*rcp;
    const Real vfm = rfm*rfm;
    const Real vfp = rfp*rfp;
    const Real xoff = Real(0.5) * ((vfm+vfp) - (vcm+vcp)) / (vcp - vcm);
    const Real yoff = (j - jc*ratio[1] + Real(0.5)) / Real(ratio[1]) - Real(0.5);
    fine(i,j,0,fcomp+ns) = crse(ic,jc,0,ccomp+ns)
        + xoff * slope(ic,jc,0,ns)
        + yoff * slope(ic,jc,0,ns+ncomp);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mf_cell_bilin_interp (int i, int j, int, int n, Array4<Real> const& fine, int fcomp,
                           Array4<Real const> const& crse, int ccomp, IntVect const& ratio) noexcept
{
    int ic = amrex::coarsen(i,ratio[0]);
    int jc = amrex::coarsen(j,ratio[1]);
    int ioff = i - ic*ratio[0];
    int joff = j - jc*ratio[1];
    int sx, sy;
    Real wx, wy;
    if (ioff*2 < ratio[0]) {
        sx = -1;
        wx = Real(ratio[0]+1+2*ioff) / Real(2*ratio[0]);
    } else {
        sx = 1;
        wx = Real(3*ratio[0]-1-2*ioff) / Real(2*ratio[0]);
    }
    if (joff*2 < ratio[1]) {
        sy = -1;
        wy = Real(ratio[1]+1+2*joff) / Real(2*ratio[1]);
    } else {
        sy = 1;
        wy = Real(3*ratio[1]-1-2*joff) / Real(2*ratio[1]);
    }
    fine(i,j,0,n+fcomp) =
        crse(ic   ,jc   ,0,n+ccomp)*           wx *           wy  +
        crse(ic+sx,jc   ,0,n+ccomp)*(Real(1.0)-wx)*           wy  +
        crse(ic   ,jc+sy,0,n+ccomp)*           wx *(Real(1.0)-wy) +
        crse(ic+sx,jc+sy,0,n+ccomp)*(Real(1.0)-wx)*(Real(1.0)-wy);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mf_nodebilin_interp (int i, int j, int, int n, Array4<Real> const& fine, int fcomp,
                          Array4<Real const> const& crse, int ccomp, IntVect const& ratio) noexcept
{
    int ic = amrex::coarsen(i,ratio[0]);
    int jc = amrex::coarsen(j,ratio[1]);
    int ioff = i - ic*ratio[0];
    int joff = j - jc*ratio[1];
    Real rxinv = Real(1.0) / Real(ratio[0]);
    Real ryinv = Real(1.0) / Real(ratio[1]);
    if (ioff != 0 && joff != 0) {
        // Node on a X-Y face
        fine(i,j,0,n+fcomp) = rxinv * ryinv *
            (crse(ic  ,jc  ,0,n+ccomp) * (ratio[0]-ioff) * (ratio[1]-joff) +
             crse(ic+1,jc  ,0,n+ccomp) * (         ioff) * (ratio[1]-joff) +
             crse(ic  ,jc+1,0,n+ccomp) * (ratio[0]-ioff) * (         joff) +
             crse(ic+1,jc+1,0,n+ccomp) * (         ioff) * (         joff));
    } else if (ioff != 0) {
        // Node on X line
        fine(i,j,0,n+fcomp) = rxinv*((ratio[0]-ioff)*crse(ic  ,jc,0,n+ccomp) +
                                     (         ioff)*crse(ic+1,jc,0,n+ccomp));
    } else if (joff != 0) {
        // Node on Y line
        fine(i,j,0,n+fcomp) = ryinv*((ratio[1]-joff)*crse(ic,jc  ,0,n+ccomp) +
                                     (         joff)*crse(ic,jc+1,0,n+ccomp));
    } else {
        // Node coincident with coarse node
        fine(i,j,0,n+fcomp) = crse(ic,jc,0,n+ccomp);
    }
}

}

#endif
