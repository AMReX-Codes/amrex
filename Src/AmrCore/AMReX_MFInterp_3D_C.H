#ifndef AMREX_MF_INTERP_3D_C_H_
#define AMREX_MF_INTERP_3D_C_H_

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mf_cell_cons_lin_interp_limit_minmax_llslope (int i, int j, int k, Array4<Real> const& slope,
                                      Array4<Real const> const& u, int scomp, int ncomp,
                                      Box const& domain, IntVect const& ratio, BCRec const* bc) noexcept
{
    Real sfx = Real(1.0);
    Real sfy = Real(1.0);
    Real sfz = Real(1.0);

    for (int ns = 0; ns < ncomp; ++ns) {
        int nu = ns + scomp;

        // x-direction
        Real dcx = mf_compute_slopes_x(i, j, k, u, nu, domain, bc[ns]);
        Real df = Real(2.0) * (u(i+1,j,k,nu) - u(i  ,j,k,nu));
        Real db = Real(2.0) * (u(i  ,j,k,nu) - u(i-1,j,k,nu));
        Real sx = (df*db >= Real(0.0)) ? amrex::min(std::abs(df),std::abs(db)) : Real(0.);
        sx = std::copysign(Real(1.),dcx)*amrex::min(sx,std::abs(dcx));
        slope(i,j,k,ns        ) = dcx;  // unlimited slope

        // y-direction
        Real dcy = mf_compute_slopes_y(i, j, k, u, nu, domain, bc[ns]);
        df = Real(2.0) * (u(i,j+1,k,nu) - u(i,j  ,k,nu));
        db = Real(2.0) * (u(i,j  ,k,nu) - u(i,j-1,k,nu));
        Real sy = (df*db >= Real(0.0)) ? amrex::min(std::abs(df),std::abs(db)) : Real(0.);
        sy = std::copysign(Real(1.),dcy)*amrex::min(sy,std::abs(dcy));
        slope(i,j,k,ns+  ncomp) = dcy;  // unlimited slope

        // z-direction
        Real dcz = mf_compute_slopes_z(i, j, k, u, nu, domain, bc[ns]);
        df = Real(2.0) * (u(i,j,k+1,nu) - u(i,j,k  ,nu));
        db = Real(2.0) * (u(i,j,k  ,nu) - u(i,j,k-1,nu));
        Real sz = (df*db >= Real(0.0)) ? amrex::min(std::abs(df),std::abs(db)) : Real(0.);
        sz = std::copysign(Real(1.),dcz)*amrex::min(sz,std::abs(dcz));
        slope(i,j,k,ns+2*ncomp) = dcz;  // unlimited slope

        // adjust limited slopes to prevent new min/max for this component
        Real alpha = 1.0;
        if (sx != Real(0.0) || sy != Real(0.0) || sz != Real(0.0)) {
            Real dumax = std::abs(sx) * Real(ratio[0]-1)/Real(2*ratio[0])
                +        std::abs(sy) * Real(ratio[1]-1)/Real(2*ratio[1])
                +        std::abs(sz) * Real(ratio[2]-1)/Real(2*ratio[2]);
            Real umax = u(i,j,k,nu);
            Real umin = u(i,j,k,nu);
            for (int koff = -1; koff <= 1; ++koff) {
            for (int joff = -1; joff <= 1; ++joff) {
            for (int ioff = -1; ioff <= 1; ++ioff) {
                umin = amrex::min(umin, u(i+ioff,j+joff,k+koff,nu));
                umax = amrex::max(umax, u(i+ioff,j+joff,k+koff,nu));
            }}}
            if (dumax * alpha > (umax - u(i,j,k,nu))) {
                alpha = (umax - u(i,j,k,nu)) / dumax;
            }
            if (dumax * alpha > (u(i,j,k,nu) - umin)) {
                alpha = (u(i,j,k,nu) - umin) / dumax;
            }
        }
        sx *= alpha;
        sy *= alpha;
        sz *= alpha;

        // for each direction, compute minimum of the ratio of limited to unlimited slopes
        if (dcx != Real(0.0)) {
            sfx = amrex::min(sfx, sx / dcx);
        }
        if (dcy != Real(0.0)) {
            sfy = amrex::min(sfy, sy / dcy);
        }
        if (dcz != Real(0.0)) {
            sfz = amrex::min(sfz, sz / dcz);
        }
    }

    // multiply unlimited slopes by the minimum of the ratio of limited to unlimited slopes
    for (int ns = 0; ns < ncomp; ++ns) {
        slope(i,j,k,ns        ) *= sfx;
        slope(i,j,k,ns+  ncomp) *= sfy;
        slope(i,j,k,ns+2*ncomp) *= sfz;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mf_cell_cons_lin_interp_llslope (int i, int j, int k, Array4<Real> const& slope,
                                      Array4<Real const> const& u, int scomp, int ncomp,
                                      Box const& domain, BCRec const* bc) noexcept
{
    Real sfx = Real(1.0);
    Real sfy = Real(1.0);
    Real sfz = Real(1.0);

    for (int ns = 0; ns < ncomp; ++ns) {
        int nu = ns + scomp;

        // x-direction
        Real dc = mf_compute_slopes_x(i, j, k, u, nu, domain, bc[ns]);
        Real df = Real(2.0) * (u(i+1,j,k,nu) - u(i  ,j,k,nu));
        Real db = Real(2.0) * (u(i  ,j,k,nu) - u(i-1,j,k,nu));
        Real sx = (df*db >= Real(0.0)) ? amrex::min(std::abs(df),std::abs(db)) : Real(0.);
        sx = std::copysign(Real(1.),dc)*amrex::min(sx,std::abs(dc));
        if (dc != Real(0.0)) {
            sfx = amrex::min(sfx, sx / dc);
        }
        slope(i,j,k,ns        ) = dc;

        // y-direction
        dc = mf_compute_slopes_y(i, j, k, u, nu, domain, bc[ns]);
        df = Real(2.0) * (u(i,j+1,k,nu) - u(i,j  ,k,nu));
        db = Real(2.0) * (u(i,j  ,k,nu) - u(i,j-1,k,nu));
        Real sy = (df*db >= Real(0.0)) ? amrex::min(std::abs(df),std::abs(db)) : Real(0.);
        sy = std::copysign(Real(1.),dc)*amrex::min(sy,std::abs(dc));
        if (dc != Real(0.0)) {
            sfy = amrex::min(sfy, sy / dc);
        }
        slope(i,j,k,ns+  ncomp) = dc;

        // z-direction
        dc = mf_compute_slopes_z(i, j, k, u, nu, domain, bc[ns]);
        df = Real(2.0) * (u(i,j,k+1,nu) - u(i,j,k  ,nu));
        db = Real(2.0) * (u(i,j,k  ,nu) - u(i,j,k-1,nu));
        Real sz = (df*db >= Real(0.0)) ? amrex::min(std::abs(df),std::abs(db)) : Real(0.);
        sz = std::copysign(Real(1.),dc)*amrex::min(sz,std::abs(dc));
        if (dc != Real(0.0)) {
            sfz = amrex::min(sfz, sz / dc);
        }
        slope(i,j,k,ns+2*ncomp) = dc;
    }

    for (int ns = 0; ns < ncomp; ++ns) {
        slope(i,j,k,ns        ) *= sfx;
        slope(i,j,k,ns+  ncomp) *= sfy;
        slope(i,j,k,ns+2*ncomp) *= sfz;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mf_cell_cons_lin_interp_mcslope (int i, int j, int k, int ns, Array4<Real> const& slope,
                                      Array4<Real const> const& u, int scomp, int ncomp,
                                      Box const& domain, IntVect const& ratio,
                                      BCRec const* bc) noexcept
{
    int nu = ns + scomp;

    // x-direction
    Real dc = mf_compute_slopes_x(i, j, k, u, nu, domain, bc[ns]);
    Real df = Real(2.0) * (u(i+1,j,k,nu) - u(i  ,j,k,nu));
    Real db = Real(2.0) * (u(i  ,j,k,nu) - u(i-1,j,k,nu));
    Real sx = (df*db >= Real(0.0)) ? amrex::min(std::abs(df),std::abs(db)) : Real(0.);
    sx = std::copysign(Real(1.),dc)*amrex::min(sx,std::abs(dc));

    // y-direction
    dc = mf_compute_slopes_y(i, j, k, u, nu, domain, bc[ns]);
    df = Real(2.0) * (u(i,j+1,k,nu) - u(i,j  ,k,nu));
    db = Real(2.0) * (u(i,j  ,k,nu) - u(i,j-1,k,nu));
    Real sy = (df*db >= Real(0.0)) ? amrex::min(std::abs(df),std::abs(db)) : Real(0.);
    sy = std::copysign(Real(1.),dc)*amrex::min(sy,std::abs(dc));

    // z-direction
    dc = mf_compute_slopes_z(i, j, k, u, nu, domain, bc[ns]);
    df = Real(2.0) * (u(i,j,k+1,nu) - u(i,j,k  ,nu));
    db = Real(2.0) * (u(i,j,k  ,nu) - u(i,j,k-1,nu));
    Real sz = (df*db >= Real(0.0)) ? amrex::min(std::abs(df),std::abs(db)) : Real(0.);
    sz = std::copysign(Real(1.),dc)*amrex::min(sz,std::abs(dc));

    Real alpha = 1.0;
    if (sx != Real(0.0) || sy != Real(0.0) || sz != Real(0.0)) {
        Real dumax = std::abs(sx) * Real(ratio[0]-1)/Real(2*ratio[0])
            +        std::abs(sy) * Real(ratio[1]-1)/Real(2*ratio[1])
            +        std::abs(sz) * Real(ratio[2]-1)/Real(2*ratio[2]);
        Real umax = u(i,j,k,nu);
        Real umin = u(i,j,k,nu);
        for (int koff = -1; koff <= 1; ++koff) {
        for (int joff = -1; joff <= 1; ++joff) {
        for (int ioff = -1; ioff <= 1; ++ioff) {
            umin = amrex::min(umin, u(i+ioff,j+joff,k+koff,nu));
            umax = amrex::max(umax, u(i+ioff,j+joff,k+koff,nu));
        }}}
        if (dumax * alpha > (umax - u(i,j,k,nu))) {
            alpha = (umax - u(i,j,k,nu)) / dumax;
        }
        if (dumax * alpha > (u(i,j,k,nu) - umin)) {
            alpha = (u(i,j,k,nu) - umin) / dumax;
        }
    }

    slope(i,j,k,ns        ) = sx * alpha;
    slope(i,j,k,ns+  ncomp) = sy * alpha;
    slope(i,j,k,ns+2*ncomp) = sz * alpha;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mf_cell_cons_lin_interp (int i, int j, int k, int ns, Array4<Real> const& fine, int fcomp,
                              Array4<Real const> const& slope, Array4<Real const> const& crse,
                              int ccomp, int ncomp, IntVect const& ratio) noexcept
{
    const int ic = amrex::coarsen(i, ratio[0]);
    const int jc = amrex::coarsen(j, ratio[1]);
    const int kc = amrex::coarsen(k, ratio[2]);
    const Real xoff = (static_cast<Real>(i - ic*ratio[0]) + Real(0.5)) / Real(ratio[0]) - Real(0.5);
    const Real yoff = (static_cast<Real>(j - jc*ratio[1]) + Real(0.5)) / Real(ratio[1]) - Real(0.5);
    const Real zoff = (static_cast<Real>(k - kc*ratio[2]) + Real(0.5)) / Real(ratio[2]) - Real(0.5);
    fine(i,j,k,fcomp+ns) = crse(ic,jc,kc,ccomp+ns)
        + xoff * slope(ic,jc,kc,ns)
        + yoff * slope(ic,jc,kc,ns+ncomp)
        + zoff * slope(ic,jc,kc,ns+ncomp*2);
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mf_cell_bilin_interp (int i, int j, int k, int n, Array4<T> const& fine, int fcomp,
                           Array4<T const> const& crse, int ccomp, IntVect const& ratio) noexcept
{
    int ic = amrex::coarsen(i,ratio[0]);
    int jc = amrex::coarsen(j,ratio[1]);
    int kc = amrex::coarsen(k,ratio[2]);
    int ioff = i - ic*ratio[0];
    int joff = j - jc*ratio[1];
    int koff = k - kc*ratio[2];
    int sx, sy, sz;
    Real wx, wy, wz;
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
    if (koff*2 < ratio[2]) {
        sz = -1;
        wz = Real(ratio[2]+1+2*koff) / Real(2*ratio[2]);
    } else {
        sz = 1;
        wz = Real(3*ratio[2]-1-2*koff) / Real(2*ratio[2]);
    }
    fine(i,j,k,n+fcomp) =
        crse(ic   ,jc   ,kc   ,n+ccomp)*           wx *           wy *           wz  +
        crse(ic+sx,jc   ,kc   ,n+ccomp)*(Real(1.0)-wx)*           wy *           wz  +
        crse(ic   ,jc+sy,kc   ,n+ccomp)*           wx *(Real(1.0)-wy)*           wz  +
        crse(ic+sx,jc+sy,kc   ,n+ccomp)*(Real(1.0)-wx)*(Real(1.0)-wy)*           wz  +
        crse(ic   ,jc   ,kc+sz,n+ccomp)*           wx *           wy *(Real(1.0)-wz) +
        crse(ic+sx,jc   ,kc+sz,n+ccomp)*(Real(1.0)-wx)*           wy *(Real(1.0)-wz) +
        crse(ic   ,jc+sy,kc+sz,n+ccomp)*           wx *(Real(1.0)-wy)*(Real(1.0)-wz) +
        crse(ic+sx,jc+sy,kc+sz,n+ccomp)*(Real(1.0)-wx)*(Real(1.0)-wy)*(Real(1.0)-wz);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mf_nodebilin_interp (int i, int j, int k, int n, Array4<Real> const& fine, int fcomp,
                          Array4<Real const> const& crse, int ccomp, IntVect const& ratio) noexcept
{
    int ic = amrex::coarsen(i,ratio[0]);
    int jc = amrex::coarsen(j,ratio[1]);
    int kc = amrex::coarsen(k,ratio[2]);
    int ioff = i - ic*ratio[0];
    int joff = j - jc*ratio[1];
    int koff = k - kc*ratio[2];
    Real rxinv = Real(1.0) / Real(ratio[0]);
    Real ryinv = Real(1.0) / Real(ratio[1]);
    Real rzinv = Real(1.0) / Real(ratio[2]);
    if (ioff != 0 && joff != 0 && koff != 0) {
        // Fine node at center of cell
        fine(i,j,k,n+fcomp) = rxinv * ryinv * rzinv *
            (crse(ic  ,jc  ,kc  ,n+ccomp) * static_cast<Real>((ratio[0]-ioff) * (ratio[1]-joff) * (ratio[2]-koff)) +
             crse(ic+1,jc  ,kc  ,n+ccomp) * static_cast<Real>((         ioff) * (ratio[1]-joff) * (ratio[2]-koff)) +
             crse(ic  ,jc+1,kc  ,n+ccomp) * static_cast<Real>((ratio[0]-ioff) * (         joff) * (ratio[2]-koff)) +
             crse(ic+1,jc+1,kc  ,n+ccomp) * static_cast<Real>((         ioff) * (         joff) * (ratio[2]-koff)) +
             crse(ic  ,jc  ,kc+1,n+ccomp) * static_cast<Real>((ratio[0]-ioff) * (ratio[1]-joff) * (         koff)) +
             crse(ic+1,jc  ,kc+1,n+ccomp) * static_cast<Real>((         ioff) * (ratio[1]-joff) * (         koff)) +
             crse(ic  ,jc+1,kc+1,n+ccomp) * static_cast<Real>((ratio[0]-ioff) * (         joff) * (         koff)) +
             crse(ic+1,jc+1,kc+1,n+ccomp) * static_cast<Real>((         ioff) * (         joff) * (         koff)));
    } else if (joff != 0 && koff != 0) {
        // Node on a Y-Z face
        fine(i,j,k,n+fcomp) = ryinv * rzinv *
            (crse(ic,jc  ,kc  ,n+ccomp) * static_cast<Real>((ratio[1]-joff) * (ratio[2]-koff)) +
             crse(ic,jc+1,kc  ,n+ccomp) * static_cast<Real>((         joff) * (ratio[2]-koff)) +
             crse(ic,jc  ,kc+1,n+ccomp) * static_cast<Real>((ratio[1]-joff) * (         koff)) +
             crse(ic,jc+1,kc+1,n+ccomp) * static_cast<Real>((         joff) * (         koff)));
    } else if (ioff != 0 && koff != 0) {
        // Node on a Z-X face
        fine(i,j,k,n+fcomp) = rxinv * rzinv *
            (crse(ic  ,jc,kc  ,n+ccomp) * static_cast<Real>((ratio[0]-ioff) * (ratio[2]-koff)) +
             crse(ic+1,jc,kc  ,n+ccomp) * static_cast<Real>((         ioff) * (ratio[2]-koff)) +
             crse(ic  ,jc,kc+1,n+ccomp) * static_cast<Real>((ratio[0]-ioff) * (         koff)) +
             crse(ic+1,jc,kc+1,n+ccomp) * static_cast<Real>((         ioff) * (         koff)));
    } else if (ioff != 0 && joff != 0) {
        // Node on a X-Y face
        fine(i,j,k,n+fcomp) = rxinv * ryinv *
            (crse(ic  ,jc  ,kc,n+ccomp) * static_cast<Real>((ratio[0]-ioff) * (ratio[1]-joff)) +
             crse(ic+1,jc  ,kc,n+ccomp) * static_cast<Real>((         ioff) * (ratio[1]-joff)) +
             crse(ic  ,jc+1,kc,n+ccomp) * static_cast<Real>((ratio[0]-ioff) * (         joff)) +
             crse(ic+1,jc+1,kc,n+ccomp) * static_cast<Real>((         ioff) * (         joff)));
    } else if (ioff != 0) {
        // Node on X line
        fine(i,j,k,n+fcomp) = rxinv*(static_cast<Real>(ratio[0]-ioff)*crse(ic  ,jc,kc,n+ccomp) +
                                     static_cast<Real>(         ioff)*crse(ic+1,jc,kc,n+ccomp));
    } else if (joff != 0) {
        // Node on Y line
        fine(i,j,k,n+fcomp) = ryinv*(static_cast<Real>(ratio[1]-joff)*crse(ic,jc  ,kc,n+ccomp) +
                                     static_cast<Real>(         joff)*crse(ic,jc+1,kc,n+ccomp));
    } else if (koff != 0) {
        // Node on Z line
        fine(i,j,k,n+fcomp) = rzinv*(static_cast<Real>(ratio[2]-koff)*crse(ic,jc,kc  ,n+ccomp) +
                                     static_cast<Real>(         koff)*crse(ic,jc,kc+1,n+ccomp));
    } else {
        // Node coincident with coarse node
        fine(i,j,k,n+fcomp) = crse(ic,jc,kc,n+ccomp);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mf_cell_quadratic_calcslope (int i, int j, int k, int n,
                                  Array4<Real const> const& crse,  int ccomp,
                                  Array4<Real>       const& slope,
                                  Box const& domain,
                                  BCRec const* bc) noexcept
{
    int nu = ccomp + n;
    Real sx  = mf_compute_slopes_x(i, j, k, crse, nu, domain, bc[n]);
    Real sy  = mf_compute_slopes_y(i, j, k, crse, nu, domain, bc[n]);
    Real sz  = mf_compute_slopes_z(i, j, k, crse, nu, domain, bc[n]);
    Real sxx = mf_cell_quadratic_compute_slopes_xx (i, j, k, crse, nu, domain, bc[n]);
    Real syy = mf_cell_quadratic_compute_slopes_yy (i, j, k, crse, nu, domain, bc[n]);
    Real szz = mf_cell_quadratic_compute_slopes_zz (i, j, k, crse, nu, domain, bc[n]);
    Real sxy = mf_cell_quadratic_compute_slopes_xy (i, j, k, crse, nu, domain, bc[n]);
    Real sxz = mf_cell_quadratic_compute_slopes_xz (i, j, k, crse, nu, domain, bc[n]);
    Real syz = mf_cell_quadratic_compute_slopes_yz (i, j, k, crse, nu, domain, bc[n]);

    slope(i,j,k,9*n  ) = sx;  // x
    slope(i,j,k,9*n+1) = sy;  // y
    slope(i,j,k,9*n+2) = sz;  // z
    slope(i,j,k,9*n+3) = sxx; // x^2
    slope(i,j,k,9*n+4) = syy; // y^2
    slope(i,j,k,9*n+5) = szz; // z^2
    slope(i,j,k,9*n+6) = sxy; // xy
    slope(i,j,k,9*n+7) = sxz; // xz
    slope(i,j,k,9*n+8) = syz; // yz
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mf_cell_quadratic_interp (int i, int j, int k, int n,
                               Array4<Real>       const& fine,  int fcomp,
                               Array4<Real const> const& crse,  int ccomp,
                               Array4<Real const> const& slope,
                               IntVect const& ratio) noexcept
{
    int ic = amrex::coarsen(i, ratio[0]);
    int jc = amrex::coarsen(j, ratio[1]);
    int kc = amrex::coarsen(k, ratio[2]);
    int irx = i - ic*ratio[0]; // = abs(i % ratio[0])
    int jry = j - jc*ratio[1]; // = abs(j % ratio[1])
    int krz = k - kc*ratio[2]; // = abs(k % ratio[2])

    // Compute offsets.
    Real xoff = ( Real(irx) + 0.5_rt ) / Real(ratio[0]) - 0.5_rt;
    Real yoff = ( Real(jry) + 0.5_rt ) / Real(ratio[1]) - 0.5_rt;
    Real zoff = ( Real(krz) + 0.5_rt ) / Real(ratio[2]) - 0.5_rt;

    fine(i,j,k,fcomp+n) = crse(ic,jc,kc,ccomp+n)
                          +   xoff               * slope(ic,jc,kc,9*n  )  // x
                          +          yoff        * slope(ic,jc,kc,9*n+1)  // y
                          +                 zoff * slope(ic,jc,kc,9*n+2)  // z
                          + 0.5_rt * xoff * xoff * slope(ic,jc,kc,9*n+3)  // x^2
                          + 0.5_rt * yoff * yoff * slope(ic,jc,kc,9*n+4)  // y^2
                          + 0.5_rt * zoff * zoff * slope(ic,jc,kc,9*n+5)  // z^2
                          +   xoff * yoff        * slope(ic,jc,kc,9*n+6)  // xy
                          +   xoff        * zoff * slope(ic,jc,kc,9*n+7)  // xz
                          +          yoff * zoff * slope(ic,jc,kc,9*n+8); // yz
}

} // namespace amrex

#endif
