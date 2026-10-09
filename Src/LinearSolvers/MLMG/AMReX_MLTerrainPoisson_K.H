#ifndef AMREX_ML_TERRAIN_POISSON_K_H_
#define AMREX_ML_TERRAIN_POISSON_K_H_
#include <AMReX_Config.H>

#include <AMReX_Array4.H>
#include <AMReX_GpuContainers.H>

namespace amrex::mlterrain {

// Domain boundary information.  bc is 0 for periodic, 1 for even
// reflection (Neumann) and -1 for odd reflection (Dirichlet).  Neumann
// faces carry zero flux.  Derivatives in the hidden direction are zero.
struct BCInfo
{
    GpuArray<int,3> lo;
    GpuArray<int,3> hi;
    GpuArray<int,3> bclo;
    GpuArray<int,3> bchi;
    int hidden = -1;

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool zero_flux (int dir, int face) const noexcept {
        return (face == lo[dir]   && bclo[dir] == 1)
            || (face == hi[dir]+1 && bchi[dir] == 1);
    }
};

// Metric terms, with averaging factors folded in: rx = h_xi/(4 h_zeta) on
// y-edges (x-face, z-node), ry = h_eta/(4 h_zeta) on x-edges (y-face,
// z-node), and on z-faces 1/h_zeta, h_xi/2 and h_eta/2.
struct Metrics
{
    Array4<Real const> rx;
    Array4<Real const> ry;
    Array4<Real const> zf;
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real metric_rx (int i, int j, int k, Array4<Real const> const& zp, Real dxinv) noexcept
{
    Real h_xi   = Real(0.25) * ( zp(i+1,j  ,k) - zp(i-1,j  ,k)
                                +zp(i+1,j+1,k) - zp(i-1,j+1,k) ) * dxinv;
    Real h_zeta = ( zp(i,j  ,k+1) - zp(i,j  ,k-1)
                   +zp(i,j+1,k+1) - zp(i,j+1,k-1) );
    return h_xi / h_zeta;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real metric_ry (int i, int j, int k, Array4<Real const> const& zp, Real dyinv) noexcept
{
    Real h_eta  = Real(0.25) * ( zp(i  ,j+1,k) - zp(i  ,j-1,k)
                                +zp(i+1,j+1,k) - zp(i+1,j-1,k) ) * dyinv;
    Real h_zeta = ( zp(i  ,j,k+1) - zp(i  ,j,k-1)
                   +zp(i+1,j,k+1) - zp(i+1,j,k-1) );
    return h_eta / h_zeta;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void metric_zf (int i, int j, int k, Array4<Real> const& zf,
                Array4<Real const> const& zp, Real dxinv, Real dyinv) noexcept
{
    zf(i,j,k,0) = Real(8.0) / ( (zp(i,j,k+1) + zp(i+1,j,k+1) + zp(i,j+1,k+1) + zp(i+1,j+1,k+1))
                               -(zp(i,j,k-1) + zp(i+1,j,k-1) + zp(i,j+1,k-1) + zp(i+1,j+1,k-1)) );
    zf(i,j,k,1) = Real(0.25) * ( zp(i+1,j+1,k) + zp(i+1,j,k) - zp(i,j+1,k) - zp(i,j,k) ) * dxinv;
    zf(i,j,k,2) = Real(0.25) * ( zp(i+1,j+1,k) + zp(i,j+1,k) - zp(i+1,j,k) - zp(i,j,k) ) * dyinv;
}

// The x-terms on x-face i: the x-gradient at constant height, and px[l],
// the x-term of the z-gradient on the edge at z-node k+l.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real face_x (int i, int j, int k, Array4<Real const> const& x, Metrics const& m,
             Real dxinv, BCInfo const& bci, Real* px) noexcept
{
    if (bci.hidden == 0) {
        px[0] = px[1] = Real(0.0);
        return Real(0.0);
    }
    Real const d1 = x(i,j,k) - x(i-1,j,k);
    Real const d0 = x(i,j,k-1) - x(i-1,j,k-1);
    Real const d2 = x(i,j,k+1) - x(i-1,j,k+1);
    Real const s0 = x(i,j,k-1) + x(i-1,j,k-1);
    Real const s1 = x(i,j,k  ) + x(i-1,j,k  );
    Real const s2 = x(i,j,k+1) + x(i-1,j,k+1);
    Real const e0 = (s1-s0) * m.rx(i,j,k);
    Real const e1 = (s2-s1) * m.rx(i,j,k+1);
    Real const hdx = Real(0.5) * dxinv;
    px[0] = (d1+d0)*hdx - Real(2.0)*e0;
    px[1] = (d2+d1)*hdx - Real(2.0)*e1;
    return d1*dxinv - (e0+e1);
}

// The y-terms on y-face j, as face_x.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real face_y (int i, int j, int k, Array4<Real const> const& x, Metrics const& m,
             Real dyinv, BCInfo const& bci, Real* py) noexcept
{
    if (bci.hidden == 1) {
        py[0] = py[1] = Real(0.0);
        return Real(0.0);
    }
    Real const d1 = x(i,j,k) - x(i,j-1,k);
    Real const d0 = x(i,j,k-1) - x(i,j-1,k-1);
    Real const d2 = x(i,j,k+1) - x(i,j-1,k+1);
    Real const s0 = x(i,j,k-1) + x(i,j-1,k-1);
    Real const s1 = x(i,j,k  ) + x(i,j-1,k  );
    Real const s2 = x(i,j,k+1) + x(i,j-1,k+1);
    Real const e0 = (s1-s0) * m.ry(i,j,k);
    Real const e1 = (s2-s1) * m.ry(i,j,k+1);
    Real const hdy = Real(0.5) * dyinv;
    py[0] = (d1+d0)*hdy - Real(2.0)*e0;
    py[1] = (d2+d1)*hdy - Real(2.0)*e1;
    return d1*dyinv - (e0+e1);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real flux_x (int i, int j, int k, Array4<Real const> const& x,
             Metrics const& m, Real dxinv, BCInfo const& bci) noexcept
{
    if (bci.zero_flux(0,i)) { return Real(0.0); }
    Real px[2];
    return -face_x(i,j,k,x,m,dxinv,bci,px);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real flux_y (int i, int j, int k, Array4<Real const> const& x,
             Metrics const& m, Real dyinv, BCInfo const& bci) noexcept
{
    if (bci.zero_flux(1,j)) { return Real(0.0); }
    Real py[2];
    return -face_y(i,j,k,x,m,dyinv,bci,py);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real flux_z (int i, int j, int k, Array4<Real const> const& x,
             Metrics const& m, Real dxinv, Real dyinv,
             BCInfo const& bci) noexcept
{
    if (bci.zero_flux(2,k)) { return Real(0.0); }
    Real g = (x(i,j,k) - x(i,j,k-1)) * m.zf(i,j,k,0);
    if (bci.hidden != 0) {
        Real p = Real(0.0);
        for (int ii = i; ii <= i+1; ++ii) {
            Real const d = (x(ii,j,k) - x(ii-1,j,k)) + (x(ii,j,k-1) - x(ii-1,j,k-1));
            Real const s = (x(ii,j,k) + x(ii-1,j,k)) - (x(ii,j,k-1) + x(ii-1,j,k-1));
            p += d*(Real(0.5)*dxinv) - Real(2.0)*(s*m.rx(ii,j,k));
        }
        g -= m.zf(i,j,k,1) * p;
    }
    if (bci.hidden != 1) {
        Real p = Real(0.0);
        for (int jj = j; jj <= j+1; ++jj) {
            Real const d = (x(i,jj,k) - x(i,jj-1,k)) + (x(i,jj,k-1) - x(i,jj-1,k-1));
            Real const s = (x(i,jj,k) + x(i,jj-1,k)) - (x(i,jj,k-1) + x(i,jj-1,k-1));
            p += d*(Real(0.5)*dyinv) - Real(2.0)*(s*m.ry(i,jj,k));
        }
        g -= m.zf(i,j,k,2) * p;
    }
    return -g;
}

// Divergence of the area-weighted fluxes in cell (i,j,k).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real adotx (int i, int j, int k, Array4<Real const> const& x,
            Metrics const& m, Array4<Real const> const& ax,
            Array4<Real const> const& ay, Array4<Real const> const& az,
            Real dxinv, Real dyinv, Real dzinv, BCInfo const& bci) noexcept
{
    Real pxl[2], pxh[2], pyl[2], pyh[2];
    Real gxlo = face_x(i  ,j,k,x,m,dxinv,bci,pxl);
    Real gxhi = face_x(i+1,j,k,x,m,dxinv,bci,pxh);
    Real gylo = face_y(i,j  ,k,x,m,dyinv,bci,pyl);
    Real gyhi = face_y(i,j+1,k,x,m,dyinv,bci,pyh);
    if (bci.zero_flux(0,i  )) { gxlo = Real(0.0); }
    if (bci.zero_flux(0,i+1)) { gxhi = Real(0.0); }
    if (bci.zero_flux(1,j  )) { gylo = Real(0.0); }
    if (bci.zero_flux(1,j+1)) { gyhi = Real(0.0); }
    Real const gzlo = bci.zero_flux(2,k  ) ? Real(0.0)
        : (x(i,j,k  ) - x(i,j,k-1)) * m.zf(i,j,k  ,0)
        - m.zf(i,j,k  ,1) * (pxl[0] + pxh[0]) - m.zf(i,j,k  ,2) * (pyl[0] + pyh[0]);
    Real const gzhi = bci.zero_flux(2,k+1) ? Real(0.0)
        : (x(i,j,k+1) - x(i,j,k  )) * m.zf(i,j,k+1,0)
        - m.zf(i,j,k+1,1) * (pxl[1] + pxh[1]) - m.zf(i,j,k+1,2) * (pyl[1] + pyh[1]);

    return (ax(i+1,j,k)*gxhi - ax(i,j,k)*gxlo) * dxinv
        +  (ay(i,j+1,k)*gyhi - ay(i,j,k)*gylo) * dyinv
        +  (az(i,j,k+1)*gzhi - az(i,j,k)*gzlo) * dzinv;
}

// Map index n in direction dir into the domain the way the ghost cell fill
// does, and update the sign accordingly.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int map_index (int n, int dir, BCInfo const& bci, int& sign) noexcept
{
    if (n < bci.lo[dir]) {
        if (bci.bclo[dir] == 0) {
            n += bci.hi[dir] - bci.lo[dir] + 1;
        } else {
            n = 2*bci.lo[dir] - 1 - n;
            sign *= bci.bclo[dir];
        }
    } else if (n > bci.hi[dir]) {
        if (bci.bchi[dir] == 0) {
            n -= bci.hi[dir] - bci.lo[dir] + 1;
        } else {
            n = 2*bci.hi[dir] + 1 - n;
            sign *= bci.bchi[dir];
        }
    }
    return n;
}

// The tridiagonal coefficients of row (i,j,k) of the column block
// klo..khi, with the boundary conditions folded in, stored in lu.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void column_coefs (int i, int j, int k, int klo, int khi, Array4<Real> const& lu,
                   Metrics const& m, Array4<Real const> const& ax,
                   Array4<Real const> const& ay, Array4<Real const> const& az,
                   Real dxinv, Real dyinv, Real dzinv, BCInfo const& bci) noexcept
{
    int im[3], jm[3], km[3], si[3], sj[3], sk[3];
    for (int n = 0; n < 3; ++n) {
        si[n] = sj[n] = sk[n] = 1;
        im[n] = map_index(i+n-1, 0, bci, si[n]);
        jm[n] = map_index(j+n-1, 1, bci, sj[n]);
        km[n] = map_index(k+n-1, 2, bci, sk[n]);
    }
    for (int t = 0; t < 3; ++t) {
        int const kt = k+t-1;
        Real coef = Real(0.0);
        if (kt >= klo && kt <= khi) {
            Real xs[27];
            for (int kk = 0; kk < 3; ++kk) {
            for (int jj = 0; jj < 3; ++jj) {
            for (int ii = 0; ii < 3; ++ii) {
                xs[ii+3*(jj+3*kk)] = (im[ii] == i && jm[jj] == j && km[kk] == kt)
                    ? Real(si[ii]*sj[jj]*sk[kk]) : Real(0.0);
            }}}
            Array4<Real const> const xa(xs, Dim3{.x=i-1, .y=j-1, .z=k-1},
                                        Dim3{.x=i+2, .y=j+2, .z=k+2}, 1);
            coef = adotx(i,j,k,xa,m,ax,ay,az,dxinv,dyinv,dzinv,bci);
        }
        lu(i,j,k,t) = coef;
    }
}

// Factor in place the column block at (i,j), given its coefficients from
// column_coefs.  lu then holds the lower coefficient, the inverse pivot and
// the upper factor.  If pin_top, the block is singular and the top
// correction is set to zero.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void factor_column (int i, int j, int klo, int khi, Array4<Real> const& lu,
                    bool pin_top) noexcept
{
    Real cprev = Real(0.0);
    Real ibprev = Real(0.0);
    for (int k = klo; k <= khi; ++k) {
        Real const a = lu(i,j,k,0);
        Real const b = lu(i,j,k,1);
        Real const c = lu(i,j,k,2);
        Real const gam = cprev * ibprev;
        Real const ib = (pin_top && k == khi) ? Real(0.0) : Real(1.0) / (b - a*gam);
        lu(i,j,k,1) = ib;
        lu(i,j,k,2) = gam;
        cprev = c;
        ibprev = ib;
    }
}

// Solve the column block for the correction, given the residual in r, and
// add the correction to sol.  r is overwritten.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void solve_column (int i, int j, int klo, int khi, Array4<Real> const& sol,
                   Array4<Real> const& r, Array4<Real const> const& lu) noexcept
{
    Real uprev = Real(0.0);
    for (int k = klo; k <= khi; ++k) {
        uprev = (r(i,j,k) - lu(i,j,k,0)*uprev) * lu(i,j,k,1);
        r(i,j,k) = uprev;
    }
    sol(i,j,khi) += uprev;
    for (int k = khi-1; k >= klo; --k) {
        uprev = r(i,j,k) - lu(i,j,k+1,2) * uprev;
        sol(i,j,k) += uprev;
    }
}

}

#endif
