#ifndef AMREX_MF_INTERP_3D_C_H_
#define AMREX_MF_INTERP_3D_C_H_

namespace amrex {

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
        Real sx = (df*db >= Real(0.0)) ? amrex::min(amrex::Math::abs(df),amrex::Math::abs(db)) : Real(0.);
        sx = amrex::Math::copysign(Real(1.),dc)*amrex::min(sx,amrex::Math::abs(dc));
        if (dc != Real(0.0)) {
            sfx = amrex::min(sfx, sx / dc);
        } else {
            sfx = Real(0.0);
        }

        // y-direction
        dc = mf_compute_slopes_y(i, j, k, u, nu, domain, bc[ns]);
        df = Real(2.0) * (u(i,j+1,k,nu) - u(i,j  ,k,nu));
        db = Real(2.0) * (u(i,j  ,k,nu) - u(i,j-1,k,nu));
        Real sy = (df*db >= Real(0.0)) ? amrex::min(amrex::Math::abs(df),amrex::Math::abs(db)) : Real(0.);
        sy = amrex::Math::copysign(Real(1.),dc)*amrex::min(sy,amrex::Math::abs(dc));
        if (dc != Real(0.0)) {
            sfy = amrex::min(sfy, sy / dc);
        } else {
            sfy = Real(0.0);
        }

        // z-direction
        dc = mf_compute_slopes_z(i, j, k, u, nu, domain, bc[ns]);
        df = Real(2.0) * (u(i,j,k+1,nu) - u(i,j,k  ,nu));
        db = Real(2.0) * (u(i,j,k  ,nu) - u(i,j,k-1,nu));
        Real sz = (df*db >= Real(0.0)) ? amrex::min(amrex::Math::abs(df),amrex::Math::abs(db)) : Real(0.);
        sz = amrex::Math::copysign(Real(1.),dc)*amrex::min(sz,amrex::Math::abs(dc));
        if (dc != Real(0.0)) {
            sfz = amrex::min(sfz, sz / dc);
        } else {
            sfz = Real(0.0);
        }

        slope(i,j,k,ns        ) = sx;
        slope(i,j,k,ns+  ncomp) = sy;
        slope(i,j,k,ns+2*ncomp) = sz;
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
    Real sx = (df*db >= Real(0.0)) ? amrex::min(amrex::Math::abs(df),amrex::Math::abs(db)) : Real(0.);
    sx = amrex::Math::copysign(Real(1.),dc)*amrex::min(sx,amrex::Math::abs(dc));

    // y-direction
    dc = mf_compute_slopes_y(i, j, k, u, nu, domain, bc[ns]);
    df = Real(2.0) * (u(i,j+1,k,nu) - u(i,j  ,k,nu));
    db = Real(2.0) * (u(i,j  ,k,nu) - u(i,j-1,k,nu));
    Real sy = (df*db >= Real(0.0)) ? amrex::min(amrex::Math::abs(df),amrex::Math::abs(db)) : Real(0.);
    sy = amrex::Math::copysign(Real(1.),dc)*amrex::min(sy,amrex::Math::abs(dc));

    // z-direction
    dc = mf_compute_slopes_z(i, j, k, u, nu, domain, bc[ns]);
    df = Real(2.0) * (u(i,j,k+1,nu) - u(i,j,k  ,nu));
    db = Real(2.0) * (u(i,j,k  ,nu) - u(i,j,k-1,nu));
    Real sz = (df*db >= Real(0.0)) ? amrex::min(amrex::Math::abs(df),amrex::Math::abs(db)) : Real(0.);
    sz = amrex::Math::copysign(Real(1.),dc)*amrex::min(sz,amrex::Math::abs(dc));

    Real alpha = 1.0;
    if (sx != Real(0.0) || sy != Real(0.0) || sz != Real(0.0)) {
        Real dumax = amrex::Math::abs(sx) * Real(ratio[0]-1)/Real(2*ratio[0])
            +        amrex::Math::abs(sy) * Real(ratio[1]-1)/Real(2*ratio[1])
            +        amrex::Math::abs(sz) * Real(ratio[2]-1)/Real(2*ratio[2]);
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
    const Real xoff = (i - ic*ratio[0] + Real(0.5)) / Real(ratio[0]) - Real(0.5);
    const Real yoff = (j - jc*ratio[1] + Real(0.5)) / Real(ratio[1]) - Real(0.5);
    const Real zoff = (k - kc*ratio[2] + Real(0.5)) / Real(ratio[2]) - Real(0.5);
    fine(i,j,k,fcomp+ns) = crse(ic,jc,kc,ccomp+ns)
        + xoff * slope(ic,jc,kc,ns)
        + yoff * slope(ic,jc,kc,ns+ncomp)
        + zoff * slope(ic,jc,kc,ns+ncomp*2);
}

}

#endif
