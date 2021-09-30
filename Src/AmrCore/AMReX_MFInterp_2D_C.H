#ifndef AMREX_MF_INTERP_2D_C_H_
#define AMREX_MF_INTERP_2D_C_H_

#include <AMReX_Geometry.H>

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

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mf_cell_quadratic_compute_slopes_xx (int i, int j, int /*k*/,
                                          Array4<Real const> const& u, int nu,
                                          Box const& domain, BCRec const& bc)
{
    Real xx = u(i-1,j,0,nu) - 2.0_rt * u(i,j,0,nu) + u(i+1,j,0,nu);
    if (i == domain.smallEnd(0) && (bc.lo(0) == BCType::ext_dir ||
                                    bc.lo(0) == BCType::hoextrap)) {
        if (i+2 < u.end.x) {
            xx = 0._rt;
        }
    }
    if (i == domain.bigEnd(0) && (bc.hi(0) == BCType::ext_dir ||
                                  bc.hi(0) == BCType::hoextrap)) {
        if (i-2 >= u.begin.x) {
            xx = 0._rt;
        }
    }
    return xx;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mf_cell_quadratic_compute_slopes_yy (int i, int j, int /*k*/,
                                          Array4<Real const> const& u, int nu,
                                          Box const& domain, BCRec const& bc)
{
    Real yy = u(i,j-1,0,nu) - 2.0_rt * u(i,j,0,nu) + u(i,j+1,0,nu);
    if (j == domain.smallEnd(1) && (bc.lo(1) == BCType::ext_dir ||
                                    bc.lo(1) == BCType::hoextrap)) {
        if (j+2 < u.end.y) {
            yy = 0._rt;
        }
    }
    if (j == domain.bigEnd(1) && (bc.hi(1) == BCType::ext_dir ||
                                  bc.hi(1) == BCType::hoextrap)) {
        if (j-2 >= u.begin.y) {
            yy = 0._rt;
        }
    }
    return yy;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mf_cell_quadratic_compute_slopes_xy (int i, int j, int /*k*/,
                                          Array4<Real const> const& u, int nu,
                                          Box const& domain, BCRec const& bc)
{
    Real xy = (1._rt/4._rt) * ( u(i-1,j-1,0,nu) - u(i+1,j-1,0,nu)
                              - u(i-1,j+1,0,nu) + u(i+1,j+1,0,nu) );
    if (i == domain.smallEnd(0) && (bc.lo(0) == BCType::ext_dir ||
                                    bc.lo(0) == BCType::hoextrap)) {
        if (i+2 < u.end.x) {
            xy = 0._rt;
        }
    }
    if (i == domain.bigEnd(0) && (bc.hi(0) == BCType::ext_dir ||
                                  bc.hi(0) == BCType::hoextrap)) {
        if (i-2 >= u.begin.x) {
            xy = 0._rt;
        }
    }
    if (j == domain.smallEnd(1) && (bc.lo(1) == BCType::ext_dir ||
                                    bc.lo(1) == BCType::hoextrap)) {
        if (j+2 < u.end.y) {
            xy = 0._rt;
        }
    }
    if (j == domain.bigEnd(1) && (bc.hi(1) == BCType::ext_dir ||
                                  bc.hi(1) == BCType::hoextrap)) {
        if (j-2 >= u.begin.y) {
            xy = 0._rt;
        }
    }
    return xy;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mf_cell_quadratic_calcslope (int i, int j, int /*k*/, int n,
                                  Array4<Real const> const& crse,  int ccomp,
                                  Array4<Real>       const& slope, int ncomp,
                                  Box const& domain,
                                  BCRec const* bc) noexcept
{
    int nu = ccomp + n;
    Real sx  = mf_compute_slopes_x(i, j, 0, crse, nu, domain, bc[n]);
    Real sy  = mf_compute_slopes_y(i, j, 0, crse, nu, domain, bc[n]);
    Real sxx = mf_cell_quadratic_compute_slopes_xx (i, j, 0, crse, nu, domain, bc[n]);
    Real syy = mf_cell_quadratic_compute_slopes_yy (i, j, 0, crse, nu, domain, bc[n]);
    Real sxy = mf_cell_quadratic_compute_slopes_xy (i, j, 0, crse, nu, domain, bc[n]);

    slope(i,j,0,n*ncomp  ) = sx;  // x
    slope(i,j,0,n*ncomp+1) = sy;  // y
    slope(i,j,0,n*ncomp+2) = sxx; // x^2
    slope(i,j,0,n*ncomp+3) = syy; // y^2
    slope(i,j,0,n*ncomp+4) = sxy; // xy
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mf_cell_quadratic_interp (int i, int j, int /*k*/, int n,
                               Array4<Real>       const& fine,  int fcomp,
                               Array4<Real const> const& crse,  int ccomp,
                               Array4<Real const> const& slope, int ncomp,
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
                          +          xoff        * slope(ic,jc,0,n*ncomp  )  // x
                          +                 yoff * slope(ic,jc,0,n*ncomp+1)  // y
                          + 0.5_rt * xoff * xoff * slope(ic,jc,0,n*ncomp+2)  // x^2
                          + 0.5_rt * yoff * yoff * slope(ic,jc,0,n*ncomp+3)  // y^2
                          +          xoff * yoff * slope(ic,jc,0,n*ncomp+4); // xy
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mf_cell_quadratic_interp_rz (int i, int j, int /*k*/, int n,
                                  Array4<Real>       const& fine,  int fcomp,
                                  Array4<Real const> const& crse,  int ccomp,
                                  Array4<Real const> const& slope, int ncomp,
                                  IntVect const& ratio,
                                  GeometryData cs_geomdata,
                                  GeometryData fn_geomdata) noexcept
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
    Real fn_vm = 0.5_rt * fn_rm * fn_rm;
    Real fn_vp = 0.5_rt * fn_rp * fn_rp;
    Real cs_vm = 0.5_rt * fn_rm * fn_rm;
    Real cs_vp = 0.5_rt * fn_rm * fn_rm;

    // Compute r offset
    Real fcen = 0.5_rt * ( fn_vm + fn_vp );
    Real ccen = 0.5_rt * ( cs_vm + cs_vp );
    Real xoff = ( fcen - ccen ) / ( cs_vp - cs_vm ) ;

    // Compute z offset (same as Cartesian)
    Real yoff = ( Real(jry) + 0.5_rt ) / Real(ratio[1]) - 0.5_rt;

    fine(i,j,0,fcomp+n) = crse(ic,jc,0,ccomp+n)
                          +          xoff        * slope(ic,jc,0,n*ncomp  )  // x
                          +                 yoff * slope(ic,jc,0,n*ncomp+1)  // y
                          + 0.5_rt * xoff * xoff * slope(ic,jc,0,n*ncomp+2)  // x^2
                          + 0.5_rt * yoff * yoff * slope(ic,jc,0,n*ncomp+3)  // y^2
                          +          xoff * yoff * slope(ic,jc,0,n*ncomp+4); // xy
}

} // namespace amrex

#endif
