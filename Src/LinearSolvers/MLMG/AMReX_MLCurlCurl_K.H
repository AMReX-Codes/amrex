#ifndef AMREX_ML_CURL_CURL_K_H_
#define AMREX_ML_CURL_CURL_K_H_
#include <AMReX_Config.H>

#include <AMReX_Array4.H>

namespace amrex {

/* Index types
 * E_x        : (0,1,1)
 * E_y        : (1,0,1)
 * E_z        : (1,1,0)
 * (curl E)_x : (1,0,0)
 * (curl E)_y : (0,1,0)
 * (curl E)_z : (0,0,1)
 */

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void mlcurlcurl_adotx_x (int i, int j, int k, Array4<Real> const& Ax,
                         Array4<Real const> const& ex,
                         Array4<Real const> const& ey,
                         Array4<Real const> const& ez,
                         Real alpha, Real beta,
                         GpuArray<Real,AMREX_SPACEDIM> const& dxinv)
{
#if (AMREX_SPACEDIM == 2)
    amrex::ignore_unused(ez);
    Real dyy = dxinv[1] * dxinv[1];
    Real dxy = dxinv[0] * dxinv[1];
    Real ccex =  ex(i  ,j  ,k  ) * dyy * Real(2.0)
        - dyy * (ex(i  ,j-1,k  ) +
                 ex(i  ,j+1,k  ))
        + dxy * (ey(i  ,j-1,k  )
               - ey(i  ,j  ,k  )
               - ey(i+1,j-1,k  )
               + ey(i+1,j  ,k  ));
#else
    Real dyy = dxinv[1] * dxinv[1];
    Real dzz = dxinv[2] * dxinv[2];
    Real dxy = dxinv[0] * dxinv[1];
    Real dxz = dxinv[0] * dxinv[2];
    Real ccex =  ex(i  ,j  ,k  ) * (dyy+dzz)*Real(2.0)
        - dyy * (ex(i  ,j-1,k  ) +
                 ex(i  ,j+1,k  ))
        - dzz * (ex(i  ,j  ,k+1) +
                 ex(i  ,j  ,k-1))
        + dxy * (ey(i  ,j-1,k  )
               - ey(i  ,j  ,k  )
               - ey(i+1,j-1,k  )
               + ey(i+1,j  ,k  ))
        + dxz * (ez(i  ,j  ,k-1)
               - ez(i  ,j  ,k  )
               - ez(i+1,j  ,k-1)
               + ez(i+1,j  ,k  ));
#endif
    Ax(i,j,k) = alpha * ccex + beta * ex(i,j,k);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void mlcurlcurl_adotx_y (int i, int j, int k, Array4<Real> const& Ay,
                         Array4<Real const> const& ex,
                         Array4<Real const> const& ey,
                         Array4<Real const> const& ez,
                         Real alpha, Real beta,
                         GpuArray<Real,AMREX_SPACEDIM> const& dxinv)
{
#if (AMREX_SPACEDIM == 2)
    amrex::ignore_unused(ez);
    Real dxx = dxinv[0] * dxinv[0];
    Real dxy = dxinv[0] * dxinv[1];
    Real ccey =  ey(i  ,j  ,k  ) * dxx * Real(2.0)
        - dxx * (ey(i-1,j  ,k  ) +
                 ey(i+1,j  ,k  ))
        + dxy * (ex(i-1,j  ,k  )
               - ex(i  ,j  ,k  )
               - ex(i-1,j+1,k  )
               + ex(i  ,j+1,k  ));
#else
    Real dxx = dxinv[0] * dxinv[0];
    Real dzz = dxinv[2] * dxinv[2];
    Real dxy = dxinv[0] * dxinv[1];
    Real dyz = dxinv[1] * dxinv[2];
    Real ccey =  ey(i  ,j  ,k  ) * (dxx+dzz)*Real(2.0)
        - dxx * (ey(i-1,j  ,k  ) +
                 ey(i+1,j  ,k  ))
        - dzz * (ey(i  ,j  ,k-1) +
                 ey(i  ,j  ,k+1))
        + dxy * (ex(i-1,j  ,k  )
               - ex(i  ,j  ,k  )
               - ex(i-1,j+1,k  )
               + ex(i  ,j+1,k  ))
        + dyz * (ez(i  ,j  ,k-1)
               - ez(i  ,j  ,k  )
               - ez(i  ,j+1,k-1)
               + ez(i  ,j+1,k  ));
#endif
    Ay(i,j,k) = alpha * ccey + beta * ey(i,j,k);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void mlcurlcurl_adotx_z (int i, int j, int k, Array4<Real> const& Az,
                         Array4<Real const> const& ex,
                         Array4<Real const> const& ey,
                         Array4<Real const> const& ez,
                         Real alpha, Real beta,
                         GpuArray<Real,AMREX_SPACEDIM> const& dxinv)
{
#if (AMREX_SPACEDIM == 2)
    amrex::ignore_unused(ex,ey);
    Real dxx = dxinv[0] * dxinv[0];
    Real dyy = dxinv[1] * dxinv[1];
    Real ccez =  ez(i  ,j  ,k  ) * (dxx+dyy)*Real(2.0)
        - dxx * (ez(i-1,j  ,k  ) +
                 ez(i+1,j  ,k  ))
        - dyy * (ez(i  ,j-1,k  ) +
                 ez(i  ,j+1,k  ));
#else
    Real dxx = dxinv[0] * dxinv[0];
    Real dyy = dxinv[1] * dxinv[1];
    Real dxz = dxinv[0] * dxinv[2];
    Real dyz = dxinv[1] * dxinv[2];
    Real ccez =  ez(i  ,j  ,k  ) * (dxx+dyy)*Real(2.0)
        - dxx * (ez(i-1,j  ,k  ) +
                 ez(i+1,j  ,k  ))
        - dyy * (ez(i  ,j-1,k  ) +
                 ez(i  ,j+1,k  ))
        + dxz * (ex(i-1,j  ,k  )
               - ex(i  ,j  ,k  )
               - ex(i-1,j  ,k+1)
               + ex(i  ,j  ,k+1))
        + dyz * (ey(i  ,j-1,k  )
               - ey(i  ,j  ,k  )
               - ey(i  ,j-1,k+1)
               + ey(i  ,j  ,k+1));
#endif
    Az(i,j,k) = alpha * ccez + beta * ez(i,j,k);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void mlcurlcurl_gsrb_x (int i, int j, int k,
                        Array4<Real      > const& ex,
                        Array4<Real const> const& ey,
                        Array4<Real const> const& ez,
                        Array4<Real const> const& rhs,
                        Real alpha, Real beta,
                        GpuArray<Real,AMREX_SPACEDIM> const& dxinv, int redblack)
{
    constexpr Real omega = Real(1.15);

    if ((i+j+k+redblack) % 2 != 0) { return; }

#if (AMREX_SPACEDIM == 2)
    amrex::ignore_unused(ez);
    Real dyy = dxinv[1] * dxinv[1];
    Real dxy = dxinv[0] * dxinv[1];
    Real gamma = alpha * (dyy)*Real(2.0) + beta;
    Real ccex =
        - dyy * (ex(i  ,j-1,k  ) +
                 ex(i  ,j+1,k  ))
        + dxy * (ey(i  ,j-1,k  )
               - ey(i  ,j  ,k  )
               - ey(i+1,j-1,k  )
               + ey(i+1,j  ,k  ));
#else
    Real dyy = dxinv[1] * dxinv[1];
    Real dzz = dxinv[2] * dxinv[2];
    Real dxy = dxinv[0] * dxinv[1];
    Real dxz = dxinv[0] * dxinv[2];
    Real gamma = alpha * (dyy+dzz)*Real(2.0) + beta;
    Real ccex =
        - dyy * (ex(i  ,j-1,k  ) +
                 ex(i  ,j+1,k  ))
        - dzz * (ex(i  ,j  ,k+1) +
                 ex(i  ,j  ,k-1))
        + dxy * (ey(i  ,j-1,k  )
               - ey(i  ,j  ,k  )
               - ey(i+1,j-1,k  )
               + ey(i+1,j  ,k  ))
        + dxz * (ez(i  ,j  ,k-1)
               - ez(i  ,j  ,k  )
               - ez(i+1,j  ,k-1)
               + ez(i+1,j  ,k  ));
#endif
    Real res = rhs(i,j,k) - (gamma*ex(i,j,k) + alpha*ccex);
    ex(i,j,k) += omega/gamma * res;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void mlcurlcurl_gsrb_y (int i, int j, int k,
                        Array4<Real const> const& ex,
                        Array4<Real      > const& ey,
                        Array4<Real const> const& ez,
                        Array4<Real const> const& rhs,
                        Real alpha, Real beta,
                        GpuArray<Real,AMREX_SPACEDIM> const& dxinv, int redblack)
{
    constexpr Real omega = Real(1.15);

    if ((i+j+k+redblack) % 2 != 0) { return; }

#if (AMREX_SPACEDIM == 2)
    amrex::ignore_unused(ez);
    Real dxx = dxinv[0] * dxinv[0];
    Real dxy = dxinv[0] * dxinv[1];
    Real gamma = alpha * (dxx)*Real(2.0) + beta;
    Real ccey =
        - dxx * (ey(i-1,j  ,k  ) +
                 ey(i+1,j  ,k  ))
        + dxy * (ex(i-1,j  ,k  )
               - ex(i  ,j  ,k  )
               - ex(i-1,j+1,k  )
               + ex(i  ,j+1,k  ));
#else
    Real dxx = dxinv[0] * dxinv[0];
    Real dzz = dxinv[2] * dxinv[2];
    Real dxy = dxinv[0] * dxinv[1];
    Real dyz = dxinv[1] * dxinv[2];
    Real gamma = alpha * (dxx+dzz)*Real(2.0) + beta;
    Real ccey =
        - dxx * (ey(i-1,j  ,k  ) +
                 ey(i+1,j  ,k  ))
        - dzz * (ey(i  ,j  ,k-1) +
                 ey(i  ,j  ,k+1))
        + dxy * (ex(i-1,j  ,k  )
               - ex(i  ,j  ,k  )
               - ex(i-1,j+1,k  )
               + ex(i  ,j+1,k  ))
        + dyz * (ez(i  ,j  ,k-1)
               - ez(i  ,j  ,k  )
               - ez(i  ,j+1,k-1)
               + ez(i  ,j+1,k  ));
#endif
    Real res = rhs(i,j,k) - (gamma*ey(i,j,k) + alpha*ccey);
    ey(i,j,k) += omega/gamma * res;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void mlcurlcurl_gsrb_z (int i, int j, int k,
                        Array4<Real const> const& ex,
                        Array4<Real const> const& ey,
                        Array4<Real      > const& ez,
                        Array4<Real const> const& rhs,
                        Real alpha, Real beta,
                        GpuArray<Real,AMREX_SPACEDIM> const& dxinv, int redblack)
{
    constexpr Real omega = Real(1.15);

    if ((i+j+k+redblack) % 2 != 0) { return; }

#if (AMREX_SPACEDIM == 2)
    amrex::ignore_unused(ex,ey);
    Real dxx = dxinv[0] * dxinv[0];
    Real dyy = dxinv[1] * dxinv[1];
    Real gamma = alpha * (dxx+dyy)*Real(2.0) + beta;
    Real ccez =
        - dxx * (ez(i-1,j  ,k  ) +
                 ez(i+1,j  ,k  ))
        - dyy * (ez(i  ,j-1,k  ) +
                 ez(i  ,j+1,k  ));
#else
    Real dxx = dxinv[0] * dxinv[0];
    Real dyy = dxinv[1] * dxinv[1];
    Real dxz = dxinv[0] * dxinv[2];
    Real dyz = dxinv[1] * dxinv[2];
    Real gamma = alpha * (dxx+dyy)*Real(2.0) + beta;
    Real ccez =
        - dxx * (ez(i-1,j  ,k  ) +
                 ez(i+1,j  ,k  ))
        - dyy * (ez(i  ,j-1,k  ) +
                 ez(i  ,j+1,k  ))
        + dxz * (ex(i-1,j  ,k  )
               - ex(i  ,j  ,k  )
               - ex(i-1,j  ,k+1)
               + ex(i  ,j  ,k+1))
        + dyz * (ey(i  ,j-1,k  )
               - ey(i  ,j  ,k  )
               - ey(i  ,j-1,k+1)
               + ey(i  ,j  ,k+1));
#endif
    Real res = rhs(i,j,k) - (gamma*ez(i,j,k) + alpha*ccez);
    ez(i,j,k) += omega/gamma * res;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void mlcurlcurl_interpadd (int dir, int i, int j, int k,
                           Array4<Real> const& fine,
                           Array4<Real const> const& crse)
{
    int ic = amrex::coarsen(i,2);
    int jc = amrex::coarsen(j,2);
    int kc = amrex::coarsen(k,2);
    if (dir == 0) {
        bool j_is_odd = (jc*2 != j);
        bool k_is_odd = (kc*2 != k);
        if (j_is_odd && k_is_odd) {
            fine(i,j,k) += Real(0.25) *
                (crse(ic,jc,kc  ) + crse(ic,jc+1,kc   ) +
                 crse(ic,jc,kc+1) + crse(ic,jc+1,kc+1));
        } else if (j_is_odd) {
            fine(i,j,k) += Real(0.5) * (crse(ic,jc,kc) + crse(ic,jc+1,kc));
        } else if (k_is_odd) {
            fine(i,j,k) += Real(0.5) * (crse(ic,jc,kc) + crse(ic,jc,kc+1));
        } else {
            fine(i,j,k) += crse(ic,jc,kc);
        }
    } else if (dir == 1) {
        bool i_is_odd = (ic*2 != i);
        bool k_is_odd = (kc*2 != k);
        if (i_is_odd && k_is_odd) {
            fine(i,j,k) += Real(0.25) *
                (crse(ic  ,jc,kc  ) + crse(ic+1,jc,kc  ) +
                 crse(ic  ,jc,kc+1) + crse(ic+1,jc,kc+1));
        } else if (i_is_odd) {
            fine(i,j,k) += Real(0.5) * (crse(ic,jc,kc) + crse(ic+1,jc,kc));
        } else if (k_is_odd) {
            fine(i,j,k) += Real(0.5) * (crse(ic,jc,kc) + crse(ic,jc,kc+1));
        } else {
            fine(i,j,k) += crse(ic,jc,kc);
        }
    } else {
        bool i_is_odd = (ic*2 != i);
        bool j_is_odd = (jc*2 != j);
        if (i_is_odd && j_is_odd) {
            fine(i,j,k) += Real(0.25) *
                (crse(ic  ,jc  ,kc) + crse(ic+1,jc  ,kc) +
                 crse(ic  ,jc+1,kc) + crse(ic+1,jc+1,kc));
        } else if (i_is_odd) {
            fine(i,j,k) += Real(0.5) * (crse(ic,jc,kc) + crse(ic+1,jc,kc));
        } else if (j_is_odd) {
            fine(i,j,k) += Real(0.5) * (crse(ic,jc,kc) + crse(ic,jc+1,kc));
        } else {
            fine(i,j,k) += crse(ic,jc,kc);
        }
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void mlcurlcurl_restriction (int dir, int i, int j, int k,
                             Array4<Real> const& crse,
                             Array4<Real const> const& fine,
                             Array4<int const> const& dmsk)
{
    int ii = i*2;
    int jj = j*2;
    int kk = k*2;
    if (dmsk(ii,jj,kk)) {
        crse(i,j,k) = Real(0.0);
    }
    else
    {
#if (AMREX_SPACEDIM == 3)
        if (dir == 0) {
            crse(i,j,k) = Real(1./32.) * (fine(ii  ,jj-1,kk-1)             +
                                          fine(ii  ,jj  ,kk-1) * Real(2.0) +
                                          fine(ii  ,jj+1,kk-1)             +
                                          fine(ii  ,jj-1,kk  ) * Real(2.0) +
                                          fine(ii  ,jj  ,kk  ) * Real(4.0) +
                                          fine(ii  ,jj+1,kk  ) * Real(2.0) +
                                          fine(ii  ,jj-1,kk+1)             +
                                          fine(ii  ,jj  ,kk+1) * Real(2.0) +
                                          fine(ii  ,jj+1,kk+1)             +
                                          fine(ii+1,jj-1,kk-1)             +
                                          fine(ii+1,jj  ,kk-1) * Real(2.0) +
                                          fine(ii+1,jj+1,kk-1)             +
                                          fine(ii+1,jj-1,kk  ) * Real(2.0) +
                                          fine(ii+1,jj  ,kk  ) * Real(4.0) +
                                          fine(ii+1,jj+1,kk  ) * Real(2.0) +
                                          fine(ii+1,jj-1,kk+1)             +
                                          fine(ii+1,jj  ,kk+1) * Real(2.0) +
                                          fine(ii+1,jj+1,kk+1)             );
        } else if (dir == 1) {
            crse(i,j,k) = Real(1./32.) * (fine(ii-1,jj  ,kk-1)             +
                                          fine(ii  ,jj  ,kk-1) * Real(2.0) +
                                          fine(ii+1,jj  ,kk-1)             +
                                          fine(ii-1,jj  ,kk  ) * Real(2.0) +
                                          fine(ii  ,jj  ,kk  ) * Real(4.0) +
                                          fine(ii+1,jj  ,kk  ) * Real(2.0) +
                                          fine(ii-1,jj  ,kk+1)             +
                                          fine(ii  ,jj  ,kk+1) * Real(2.0) +
                                          fine(ii+1,jj  ,kk+1)             +
                                          fine(ii-1,jj+1,kk-1)             +
                                          fine(ii  ,jj+1,kk-1) * Real(2.0) +
                                          fine(ii+1,jj+1,kk-1)             +
                                          fine(ii-1,jj+1,kk  ) * Real(2.0) +
                                          fine(ii  ,jj+1,kk  ) * Real(4.0) +
                                          fine(ii+1,jj+1,kk  ) * Real(2.0) +
                                          fine(ii-1,jj+1,kk+1)             +
                                          fine(ii  ,jj+1,kk+1) * Real(2.0) +
                                          fine(ii+1,jj+1,kk+1)             );
        } else {
            crse(i,j,k) = Real(1./32.) * (fine(ii-1,jj-1,kk  )             +
                                          fine(ii  ,jj-1,kk  ) * Real(2.0) +
                                          fine(ii+1,jj-1,kk  )             +
                                          fine(ii-1,jj  ,kk  ) * Real(2.0) +
                                          fine(ii  ,jj  ,kk  ) * Real(4.0) +
                                          fine(ii+1,jj  ,kk  ) * Real(2.0) +
                                          fine(ii-1,jj+1,kk  )             +
                                          fine(ii  ,jj+1,kk  ) * Real(2.0) +
                                          fine(ii+1,jj+1,kk  )             +
                                          fine(ii-1,jj-1,kk+1)             +
                                          fine(ii  ,jj-1,kk+1) * Real(2.0) +
                                          fine(ii+1,jj-1,kk+1)             +
                                          fine(ii-1,jj  ,kk+1) * Real(2.0) +
                                          fine(ii  ,jj  ,kk+1) * Real(4.0) +
                                          fine(ii+1,jj  ,kk+1) * Real(2.0) +
                                          fine(ii-1,jj+1,kk+1)             +
                                          fine(ii  ,jj+1,kk+1) * Real(2.0) +
                                          fine(ii+1,jj+1,kk+1)             );
        }
#else
        if (dir == 0) {
            crse(i,j,0) = Real(0.125) * (fine(ii  ,jj-1,0)             +
                                         fine(ii  ,jj  ,0) * Real(2.0) +
                                         fine(ii  ,jj+1,0)             +
                                         fine(ii+1,jj-1,0)             +
                                         fine(ii+1,jj  ,0) * Real(2.0) +
                                         fine(ii+1,jj+1,0)             );
        } else if (dir == 1) {
            crse(i,j,0) = Real(0.125) * (fine(ii-1,jj  ,0)             +
                                         fine(ii  ,jj  ,0) * Real(2.0) +
                                         fine(ii+1,jj  ,0)             +
                                         fine(ii-1,jj+1,0)             +
                                         fine(ii  ,jj+1,0) * Real(2.0) +
                                         fine(ii+1,jj+1,0)             );
        } else {
            crse(i,j,0) = Real(1./16.)*(fine(ii-1,jj-1,0) + Real(2.)*fine(ii  ,jj-1,0) +          fine(ii+1,jj-1,0)
                             + Real(2.)*fine(ii-1,jj  ,0) + Real(4.)*fine(ii  ,jj  ,0) + Real(2.)*fine(ii+1,jj  ,0)
                                      + fine(ii-1,jj+1,0) + Real(2.)*fine(ii  ,jj+1,0) +          fine(ii+1,jj+1,0));

        }
#endif
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void mlcurlcurl_bc_symmetry (int i, int j, int k, Orientation face, IndexType it,
                             Array4<Real> const& a)
{
    int const idir = face.coordDir();
    int offset = face.isLow() ? 1 : -1;
    Real sign;
    if (it.cellCentered(idir)) {
        sign = Real(-1.0);
    } else {
        sign = Real(1.0);
        offset *= 2;
    }

    if (idir == 0) {
        a(i,j,k) = sign * a(i+offset,j,k);
    } else if (idir == 1) {
        a(i,j,k) = sign * a(i,j+offset,k);
    } else {
        a(i,j,k) = sign * a(i,j,k+offset);
    }
}

}

#endif
