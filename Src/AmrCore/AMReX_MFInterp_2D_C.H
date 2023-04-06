#ifndef AMREX_MF_INTERP_2D_C_H_
#define AMREX_MF_INTERP_2D_C_H_

#include <AMReX_Geometry.H>

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mf_cell_cons_lin_interp_limit_minmax_llslope (int i, int j, int, Array4<Real> const& slope,
                                      Array4<Real const> const& u, int scomp, int ncomp,
                                      Box const& domain, IntVect const& ratio, BCRec const* bc) noexcept
{
    Real sfx = Real(1.0);
    Real sfy = Real(1.0);

    for (int ns = 0; ns < ncomp; ++ns) {
        int nu = ns + scomp;

        // x-direction
        Real dcx = mf_compute_slopes_x(i, j, 0, u, nu, domain, bc[ns]);
        Real df = Real(2.0) * (u(i+1,j,0,nu) - u(i  ,j,0,nu));
        Real db = Real(2.0) * (u(i  ,j,0,nu) - u(i-1,j,0,nu));
        Real sx = (df*db >= Real(0.0)) ? amrex::min(std::abs(df),std::abs(db)) : Real(0.);
        sx = std::copysign(Real(1.),dcx)*amrex::min(sx,std::abs(dcx));
        slope(i,j,0,ns        ) = dcx;

        // y-direction
        Real dcy = mf_compute_slopes_y(i, j, 0, u, nu, domain, bc[ns]);
        df = Real(2.0) * (u(i,j+1,0,nu) - u(i,j  ,0,nu));
        db = Real(2.0) * (u(i,j  ,0,nu) - u(i,j-1,0,nu));
        Real sy = (df*db >= Real(0.0)) ? amrex::min(std::abs(df),std::abs(db)) : Real(0.);
        sy = std::copysign(Real(1.),dcy)*amrex::min(sy,std::abs(dcy));
        slope(i,j,0,ns+  ncomp) = dcy;

        // adjust limited slopes to prevent new min/max for this component
        Real alpha = Real(1.0);
        if (sx != Real(0.0) || sy != Real(0.0)) {
            Real dumax = std::abs(sx) * Real(ratio[0]-1)/Real(2*ratio[0])
                +        std::abs(sy) * Real(ratio[1]-1)/Real(2*ratio[1]);
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
        sx *= alpha;
        sy *= alpha;

        // for each direction, compute minimum of the ratio of limited to unlimited slopes
        if (dcx != Real(0.0)) {
            sfx = amrex::min(sfx, sx / dcx);
        }
        if (dcy != Real(0.0)) {
            sfy = amrex::min(sfy, sy / dcy);
        }
    }

    // multiply unlimited slopes by the minimum of the ratio of limited to unlimited slopes
    for (int ns = 0; ns < ncomp; ++ns) {
        slope(i,j,0,ns        ) *= sfx;
        slope(i,j,0,ns+  ncomp) *= sfy;
    }
}

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
        Real sx = (df*db >= Real(0.0)) ? amrex::min(std::abs(df),std::abs(db)) : Real(0.);
        sx = std::copysign(Real(1.),dc)*amrex::min(sx,std::abs(dc));
        if (dc != Real(0.0)) {
            sfx = amrex::min(sfx, sx / dc);
        }
        slope(i,j,0,ns        ) = dc;

        // y-direction
        dc = mf_compute_slopes_y(i, j, 0, u, nu, domain, bc[ns]);
        df = Real(2.0) * (u(i,j+1,0,nu) - u(i,j  ,0,nu));
        db = Real(2.0) * (u(i,j  ,0,nu) - u(i,j-1,0,nu));
        Real sy = (df*db >= Real(0.0)) ? amrex::min(std::abs(df),std::abs(db)) : Real(0.);
        sy = std::copysign(Real(1.),dc)*amrex::min(sy,std::abs(dc));
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
    Real sx = (df*db >= Real(0.0)) ? amrex::min(std::abs(df),std::abs(db)) : Real(0.);
    sx = std::copysign(Real(1.),dc)*amrex::min(sx,std::abs(dc));

    // y-direction
    dc = mf_compute_slopes_y(i, j, 0, u, nu, domain, bc[ns]);
    df = Real(2.0) * (u(i,j+1,0,nu) - u(i,j  ,0,nu));
    db = Real(2.0) * (u(i,j  ,0,nu) - u(i,j-1,0,nu));
    Real sy = (df*db >= Real(0.0)) ? amrex::min(std::abs(df),std::abs(db)) : Real(0.);
    sy = std::copysign(Real(1.),dc)*amrex::min(sy,std::abs(dc));

    Real alpha = Real(1.0);
    if (sx != Real(0.0) || sy != Real(0.0)) {
        Real dumax = std::abs(sx) * Real(ratio[0]-1)/Real(2*ratio[0])
            +        std::abs(sy) * Real(ratio[1]-1)/Real(2*ratio[1]);
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
    const Real xoff = (static_cast<Real>(i - ic*ratio[0]) + Real(0.5)) / Real(ratio[0]) - Real(0.5);
    const Real yoff = (static_cast<Real>(j - jc*ratio[1]) + Real(0.5)) / Real(ratio[1]) - Real(0.5);
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
    Real sx = (df*db >= Real(0.0)) ? amrex::min(std::abs(df),std::abs(db)) : Real(0.);
    sx = std::copysign(Real(1.),dc)*amrex::min(sx,std::abs(dc));

    // y-direction
    dc = mf_compute_slopes_y(i, j, 0, u, nu, domain, bc[ns]);
    df = Real(2.0) * (u(i,j+1,0,nu) - u(i,j  ,0,nu));
    db = Real(2.0) * (u(i,j  ,0,nu) - u(i,j-1,0,nu));
    Real sy = (df*db >= Real(0.0)) ? amrex::min(std::abs(df),std::abs(db)) : Real(0.);
    sy = std::copysign(Real(1.),dc)*amrex::min(sy,std::abs(dc));

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
            + std::abs(sy) * Real(ratio[1]-1)/Real(2*ratio[1]);
        Real dumin = -amrex::min(sx*xlo, sx*xhi)
            + std::abs(sy) * Real(ratio[1]-1)/Real(2*ratio[1]);
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

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mf_cell_bilin_interp (int i, int j, int, int n, Array4<T> const& fine, int fcomp,
                           Array4<T const> const& crse, int ccomp, IntVect const& ratio) noexcept
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
            (crse(ic  ,jc  ,0,n+ccomp) * static_cast<Real>((ratio[0]-ioff) * (ratio[1]-joff)) +
             crse(ic+1,jc  ,0,n+ccomp) * static_cast<Real>((         ioff) * (ratio[1]-joff)) +
             crse(ic  ,jc+1,0,n+ccomp) * static_cast<Real>((ratio[0]-ioff) * (         joff)) +
             crse(ic+1,jc+1,0,n+ccomp) * static_cast<Real>((         ioff) * (         joff)));
    } else if (ioff != 0) {
        // Node on X line
        fine(i,j,0,n+fcomp) = rxinv*(static_cast<Real>(ratio[0]-ioff)*crse(ic  ,jc,0,n+ccomp) +
                                     static_cast<Real>(         ioff)*crse(ic+1,jc,0,n+ccomp));
    } else if (joff != 0) {
        // Node on Y line
        fine(i,j,0,n+fcomp) = ryinv*(static_cast<Real>(ratio[1]-joff)*crse(ic,jc  ,0,n+ccomp) +
                                     static_cast<Real>(         joff)*crse(ic,jc+1,0,n+ccomp));
    } else {
        // Node coincident with coarse node
        fine(i,j,0,n+fcomp) = crse(ic,jc,0,n+ccomp);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mf_cell_quadratic_calcslope (int i, int j, int /*k*/, int n,
                                  Array4<Real const> const& crse,  int ccomp,
                                  Array4<Real>       const& slope,
                                  Box const& domain,
                                  BCRec const* bc) noexcept
{
    int nu = ccomp + n;
    Real sx  = mf_compute_slopes_x(i, j, 0, crse, nu, domain, bc[n]);
    Real sy  = mf_compute_slopes_y(i, j, 0, crse, nu, domain, bc[n]);
    Real sxx = mf_cell_quadratic_compute_slopes_xx (i, j, 0, crse, nu, domain, bc[n]);
    Real syy = mf_cell_quadratic_compute_slopes_yy (i, j, 0, crse, nu, domain, bc[n]);
    Real sxy = mf_cell_quadratic_compute_slopes_xy (i, j, 0, crse, nu, domain, bc[n]);

    slope(i,j,0,5*n  ) = sx;  // x
    slope(i,j,0,5*n+1) = sy;  // y
    slope(i,j,0,5*n+2) = sxx; // x^2
    slope(i,j,0,5*n+3) = syy; // y^2
    slope(i,j,0,5*n+4) = sxy; // xy
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mf_cell_quadratic_interp (int i, int j, int /*k*/, int n,
                               Array4<Real>       const& fine,  int fcomp,
                               Array4<Real const> const& crse,  int ccomp,
                               Array4<Real const> const& slope,
                               IntVect const& ratio) noexcept
{
    int ic = amrex::coarsen(i, ratio[0]);
    int jc = amrex::coarsen(j, ratio[1]);
    int irx = i - ic*ratio[0]; // = abs(i % ratio[0])
    int jry = j - jc*ratio[1]; // = abs(j % ratio[1])

    // Compute offsets.
    Real xoff = ( Real(irx) + 0.5_rt ) / Real(ratio[0]) - 0.5_rt;
    Real yoff = ( Real(jry) + 0.5_rt ) / Real(ratio[1]) - 0.5_rt;

    fine(i,j,0,fcomp+n) = crse(ic,jc,0,ccomp+n)
                          +          xoff        * slope(ic,jc,0,5*n  )  // x
                          +                 yoff * slope(ic,jc,0,5*n+1)  // y
                          + 0.5_rt * xoff * xoff * slope(ic,jc,0,5*n+2)  // x^2
                          + 0.5_rt * yoff * yoff * slope(ic,jc,0,5*n+3)  // y^2
                          +          xoff * yoff * slope(ic,jc,0,5*n+4); // xy
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mf_cell_quadratic_interp_rz (int i, int j, int /*k*/, int n,
                                  Array4<Real>       const& fine,  int fcomp,
                                  Array4<Real const> const& crse,  int ccomp,
                                  Array4<Real const> const& slope,
                                  IntVect const& ratio,
                                  GeometryData const& cs_geomdata,
                                  GeometryData const& fn_geomdata) noexcept
{
    int ic = amrex::coarsen(i, ratio[0]);
    int jc = amrex::coarsen(j, ratio[1]);
    int jry = j - jc*ratio[1]; // = abs(j % ratio[1])

    // Extract geometry data for radial direction.
    Real fn_dr  = fn_geomdata.CellSize(0);
    Real fn_rlo = fn_geomdata.ProbLo(0);
    Real cs_dr  = cs_geomdata.CellSize(0);
    Real cs_rlo = cs_geomdata.ProbLo(0);

    // Compute radial mesh.
    Real fn_rm = fn_rlo + Real(i)    * fn_dr;
    Real fn_rp = fn_rlo + Real(i+1)  * fn_dr;
    Real cs_rm = cs_rlo + Real(ic)   * cs_dr;
    Real cs_rp = cs_rlo + Real(ic+1) * cs_dr;

    // Compute radial cell edge coords.
    // Note: omit the 0.5 prefactor here -- we'll add it once at the last step
    Real fn_vm = fn_rm*fn_rm;
    Real fn_vp = fn_rp*fn_rp;
    Real fcen  = fn_vm + fn_vp;
    Real cs_vm = cs_rm*cs_rm;
    Real cs_vp = cs_rp*cs_rp;
    Real ccen  = cs_vm + cs_vp;

    // Compute r offset
    Real xoff = 0.5_rt * ( fcen - ccen ) / ( cs_vp - cs_vm ) ;

    // Compute z offset (same formula as Cartesian)
    Real yoff = ( Real(jry) + 0.5_rt ) / Real(ratio[1]) - 0.5_rt;

    fine(i,j,0,fcomp+n) = crse(ic,jc,0,ccomp+n)
                          +          xoff        * slope(ic,jc,0,5*n  )  // x
                          +                 yoff * slope(ic,jc,0,5*n+1)  // y
                          + 0.5_rt * xoff * xoff * slope(ic,jc,0,5*n+2)  // x^2
                          + 0.5_rt * yoff * yoff * slope(ic,jc,0,5*n+3)  // y^2
                          +          xoff * yoff * slope(ic,jc,0,5*n+4); // xy
}

} // namespace amrex

#endif
