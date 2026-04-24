#ifndef AMREX_ML_CURL_CURL_K_H_
#define AMREX_ML_CURL_CURL_K_H_
#include <AMReX_Config.H>

#include <AMReX_Array4.H>
#include <AMReX_LUSolver.H>
#include <AMReX_PCGSolver.H>

namespace amrex {

/* Index types
 * E_x        : (0,1,1)
 * E_y        : (1,0,1)
 * E_z        : (1,1,0)
 * (curl E)_x : (1,0,0)
 * (curl E)_y : (0,1,0)
 * (curl E)_z : (0,0,1)
 */

/*
  Notes for gs4:

  Interior nodes:

    For v = [Ex(i-1,j,k), Ex(i,j,k), Ey(i,j-1,k), Ey(i,j,k), Ez(i,j,k-1), Ez(i,j,k)]^T,
    we have A*v = b, where

    A00 = alpha*(dyy+dzz)*2 + beta
    A01 = 0
    A02 = -alpha*dxy
    A03 =  alpha*dxy
    A04 = -alpha*dxz
    A05 =  alpha*dxz

    A10 = 0
    A11 = alpha*(dyy+dzz)*2 + beta
    A12 =  alpha*dxy
    A13 = -alpha*dxy
    A14 =  alpha*dxz
    A15 = -alpha*dxz

    A20 = -alpha*dxy
    A21 =  alpha*dxy
    A22 = alpha*(dxx+dzz)*2 + beta
    A23 = 0
    A24 = -alpha*dyz
    A25 =  alpha*dyz

    A30 =  alpha*dxy
    A31 = -alpha*dxy
    A32 = 0
    A33 = alpha*(dxx+dzz)*2 + beta
    A34 =  alpha*dyz
    A35 = -alpha*dyz

    A40 = -alpha*dxz
    A41 =  alpha*dxz
    A42 = -alpha*dyz
    A43 =  alpha*dyz
    A44 = alpha*(dxx+dyy)*2 + beta
    A45 = 0

    A50 =  alpha*dxz
    A51 = -alpha*dxz
    A52 =  alpha*dyz
    A53 = -alpha*dyz
    A54 = 0
    A55 = alpha*(dxx+dyy)*2 + beta

    b0 = rhsx(i-1,j,k) - (alpha*ccex), where
      ccex = - dyy * (ex(i-1,j-1,k  ) +
                      ex(i-1,j+1,k  ))
             - dzz * (ex(i-1,j  ,k+1) +
                      ex(i-1,j  ,k-1))
             + dxy * (ey(i-1,j-1,k  )
                    - ey(i-1,j  ,k  ))
             + dxz * (ez(i-1,j  ,k-1)
                    - ez(i-1,j  ,k  ))
    b1 = rhsx(i,j,k) - (alpha*ccex), where
      ccex = - dyy * ( ex(i  ,j-1,k  ) +
                       ex(i  ,j+1,k  ))
             - dzz * ( ex(i  ,j  ,k+1) +
                       ex(i  ,j  ,k-1))
             + dxy * (-ey(i+1,j-1,k  )
                     + ey(i+1,j  ,k  ))
             + dxz * (-ez(i+1,j  ,k-1)
                     + ez(i+1,j  ,k  ));
    b2 = rhsy(i,j-1,k) - alpha*ccey, where
      ccey = - dxx * (ey(i-1,j-1,k  ) +
                      ey(i+1,j-1,k  ))
             - dzz * (ey(i  ,j-1,k-1) +
                      ey(i  ,j-1,k+1))
             + dxy * (ex(i-1,j-1,k  )
                    - ex(i  ,j-1,k  ))
             + dyz * (ez(i  ,j-1,k-1)
                    - ez(i  ,j-1,k  ))
    b3 = rhsy(i,j,k) - alpha*ccey, where
      ccey = - dxx * ( ey(i-1,j  ,k  ) +
                       ey(i+1,j  ,k  ))
             - dzz * ( ey(i  ,j  ,k-1) +
                       ey(i  ,j  ,k+1))
             + dxy * (-ex(i-1,j+1,k  )
                     + ex(i  ,j+1,k  ))
             + dyz * (-ez(i  ,j+1,k-1)
                     + ez(i  ,j+1,k  ));
    b4 = rhsz(i,j,k-1) - alpha*ccez, where
      ccez = - dxx * (ez(i-1,j  ,k-1) +
                      ez(i+1,j  ,k-1))
             - dyy * (ez(i  ,j-1,k-1) +
                      ez(i  ,j+1,k-1))
             + dxz * (ex(i-1,j  ,k-1)
                    - ex(i  ,j  ,k-1))
             + dyz * (ey(i  ,j-1,k-1)
                    - ey(i  ,j  ,k-1))
    b5 = rhsz(i,j,k) - alpha*ccez, where
      ccez = - dxx * ( ez(i-1,j  ,k  ) +
                       ez(i+1,j  ,k  ))
             - dyy * ( ez(i  ,j-1,k  ) +
                       ez(i  ,j+1,k  ))
             + dxz * (-ex(i-1,j  ,k+1)
                     + ex(i  ,j  ,k+1))
             + dyz * (-ey(i  ,j-1,k+1)
                     + ey(i  ,j  ,k+1));

    dxx = 1/(dx*dx)
    dyy = 1/(dy*dy)
    dzz = 1/(dz*dz)
    dxy = 1/(dx*dy)
    dxz = 1/(dx*dz)
    dyz = 1/(dy*dz)

  For Dirichlet boundary nodes, we don't do anything.

  For symmetric boundary nodes, we treat it as interior nodes because the
  rhs outside the domain has been filled properly.

  In 2D,

    For v = [Ex(i-1,j,k), Ex(i,j,k), Ey(i,j-1,k), Ey(i,j,k)]^T,
    we have A*v = b, where

    A00 = alpha*dyy*2 + beta
    A01 = 0
    A02 = -alpha*dxy
    A03 =  alpha*dxy

    A10 = 0
    A11 = alpha*dyy*2 + beta
    A12 =  alpha*dxy
    A13 = -alpha*dxy

    A20 = -alpha*dxy
    A21 =  alpha*dxy
    A22 = alpha*dxx*2 + beta
    A23 = 0

    A30 =  alpha*dxy
    A31 = -alpha*dxy
    A32 = 0
    A33 = alpha*dxx*2 + beta

    b0 = rhsx(i-1,j,k) - (alpha*ccex), where
      ccex = - dyy * (ex(i-1,j-1,k  ) +
                      ex(i-1,j+1,k  ))
             + dxy * (ey(i-1,j-1,k  )
                    - ey(i-1,j  ,k  ))
    b1 = rhsx(i,j,k) - (alpha*ccex), where
      ccex = - dyy * ( ex(i  ,j-1,k  ) +
                       ex(i  ,j+1,k  ))
             + dxy * (-ey(i+1,j-1,k  )
                     + ey(i+1,j  ,k  ))
    b2 = rhsy(i,j-1,k) - alpha*ccey, where
      ccey = - dxx * (ey(i-1,j-1,k  ) +
                      ey(i+1,j-1,k  ))
             + dxy * (ex(i-1,j-1,k  )
                    - ex(i  ,j-1,k  ))
    b3 = rhsy(i,j,k) - alpha*ccey, where
      ccey = - dxx * ( ey(i-1,j  ,k  ) +
                       ey(i+1,j  ,k  ))
             + dxy * (-ex(i-1,j+1,k  )
                     + ex(i  ,j+1,k  ))
*/

/// \cond DOXYGEN_IGNORE
struct CurlCurlDirichletInfo
{
    IntVect dirichlet_lo;
    IntVect dirichlet_hi;
#if (AMREX_SPACEDIM < 3)
    int coord = 0;
#endif

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool is_dirichlet_node (int i, int j, int k) const
    {
#if (AMREX_SPACEDIM == 1)
        amrex::ignore_unused(j,k);
        return (i == dirichlet_lo[0]) || (i == dirichlet_hi[0]);
#elif (AMREX_SPACEDIM == 2)
        amrex::ignore_unused(k);
        return (i == dirichlet_lo[0]) || (i == dirichlet_hi[0])
            || (j == dirichlet_lo[1]) || (j == dirichlet_hi[1]);
#else
        return (i == dirichlet_lo[0]) || (i == dirichlet_hi[0])
            || (j == dirichlet_lo[1]) || (j == dirichlet_hi[1])
            || (k == dirichlet_lo[2]) || (k == dirichlet_hi[2]);
#endif
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool is_dirichlet_x_edge (int, int j, int k) const // NOLINT(readability-convert-member-functions-to-static)
    {
#if (AMREX_SPACEDIM == 1)
        amrex::ignore_unused(j, k);
        return false;
#elif (AMREX_SPACEDIM == 2)
        amrex::ignore_unused(k);
        return (j == dirichlet_lo[1]) || (j == dirichlet_hi[1]);
#else
        return (j == dirichlet_lo[1]) || (j == dirichlet_hi[1])
            || (k == dirichlet_lo[2]) || (k == dirichlet_hi[2]);
#endif
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool is_dirichlet_y_edge (int i, int, int k) const
    {
#if (AMREX_SPACEDIM < 3)
        amrex::ignore_unused(k);
        return (i == dirichlet_lo[0]) || (i == dirichlet_hi[0]);
#else
        return (i == dirichlet_lo[0]) || (i == dirichlet_hi[0])
            || (k == dirichlet_lo[2]) || (k == dirichlet_hi[2]);
#endif
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool is_dirichlet_z_edge (int i, int j, int) const
    {
#if (AMREX_SPACEDIM == 1)
        amrex::ignore_unused(j);
        if ((i == 0) && (coord == 1)) {
            return false; // Ez in 1d cyl is not dirichlet
        } else {
            return (i == dirichlet_lo[0]) || (i == dirichlet_hi[0]);
        }
#else
        return (i == dirichlet_lo[0]) || (i == dirichlet_hi[0])
            || (j == dirichlet_lo[1]) || (j == dirichlet_hi[1]);
#endif
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool is_dirichlet_edge (int dim, int i, int j, int k) const
    {
        if (dim == 0) {
            return is_dirichlet_x_edge(i,j,k);
        } else if (dim == 1) {
            return is_dirichlet_y_edge(i,j,k);
        } else {
            return is_dirichlet_z_edge(i,j,k);
        }
    }
};

struct CurlCurlSymmetryInfo
{
    IntVect symmetry_lo;
    IntVect symmetry_hi;

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool xlo_is_symmetric (int i) const
    {
        return i == symmetry_lo[0];
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool xhi_is_symmetric (int i) const
    {
        return i == symmetry_hi[0];
    }

#if (AMREX_SPACEDIM > 1)
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool ylo_is_symmetric (int j) const
    {
        return j == symmetry_lo[1];
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool yhi_is_symmetric (int j) const
    {
        return j == symmetry_hi[1];
    }
#endif

#if (AMREX_SPACEDIM == 3)
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool zlo_is_symmetric (int k) const
    {
        return k == symmetry_lo[2];
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool zhi_is_symmetric (int k) const
    {
        return k == symmetry_hi[2];
    }
#endif

    [[nodiscard]] bool is_symmetric (int dir, int side, int idx) const
    {
#if (AMREX_SPACEDIM == 1)
        amrex::ignore_unused(dir);
        return (side == 0) ? xlo_is_symmetric(idx) : xhi_is_symmetric(idx);
#elif (AMREX_SPACEDIM == 2)
        if (dir == 0) {
            return (side == 0) ? xlo_is_symmetric(idx) : xhi_is_symmetric(idx);
        } else {
            return (side == 0) ? ylo_is_symmetric(idx) : yhi_is_symmetric(idx);
        }
#else
        if (dir == 0) {
            return (side == 0) ? xlo_is_symmetric(idx) : xhi_is_symmetric(idx);
        } else if (dir == 1) {
            return (side == 0) ? ylo_is_symmetric(idx) : yhi_is_symmetric(idx);
        } else {
            return (side == 0) ? zlo_is_symmetric(idx) : zhi_is_symmetric(idx);
        }
#endif
    }
};

enum struct CurlCurlStateType { x, b, r }; // x, b & r=b-Ax

/// \endcond

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void mlcurlcurl_adotx_x (int i, int j, int k, Array4<Real> const& Ax,
                         Array4<Real const> const& ex,
                         Array4<Real const> const& ey,
                         Array4<Real const> const& ez,
                         Real beta, GpuArray<Real,AMREX_SPACEDIM> const& adxinv)
{
#if (AMREX_SPACEDIM == 1)
    amrex::ignore_unused(ey,ez,adxinv);
    Real ccex = 0;
#elif (AMREX_SPACEDIM == 2)
    amrex::ignore_unused(ez);
    Real dyy = adxinv[1] * adxinv[1];
    Real dxy = adxinv[0] * adxinv[1];
    Real ccex =  ex(i  ,j  ,k  ) * dyy * Real(2.0)
        - dyy * (ex(i  ,j-1,k  ) +
                 ex(i  ,j+1,k  ))
        + dxy * (ey(i  ,j-1,k  )
               - ey(i  ,j  ,k  )
               - ey(i+1,j-1,k  )
               + ey(i+1,j  ,k  ));
#elif (AMREX_SPACEDIM == 3)
    Real dyy = adxinv[1] * adxinv[1];
    Real dzz = adxinv[2] * adxinv[2];
    Real dxy = adxinv[0] * adxinv[1];
    Real dxz = adxinv[0] * adxinv[2];
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
    Ax(i,j,k) = ccex + beta * ex(i,j,k);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void mlcurlcurl_adotx_y (int i, int j, int k, Array4<Real> const& Ay,
                         Array4<Real const> const& ex,
                         Array4<Real const> const& ey,
                         Array4<Real const> const& ez,
                         Real beta, GpuArray<Real,AMREX_SPACEDIM> const& adxinv
#if (AMREX_SPACEDIM < 3)
                         , int coord
#endif
                        )
{
#if (AMREX_SPACEDIM == 1)
    amrex::ignore_unused(ex,ez);
    Real ccey;
    if (coord == 0) {
        ccey =  ey(i  ,j,k) * Real(2.0)
            -   ey(i-1,j,k)
            -   ey(i+1,j,k);
    } else if (coord == 1) {
        // i=0 won't get here
        ccey =  ey(i  ,j,k) * (Real( 1.0)/(Real(i)*Real(i)) + Real(2.0))
            +   ey(i-1,j,k) * (Real( 0.5)/Real(i) - Real(1.0))
            +   ey(i+1,j,k) * (Real(-0.5)/Real(i) - Real(1.0));
    } else {
        // i=0 won't get here
        ccey =  ey(i  ,j,k) * Real(2.0)
            +   ey(i-1,j,k) * (Real( 1.0)/Real(i) - Real(1.0))
            +   ey(i+1,j,k) * (Real(-1.0)/Real(i) - Real(1.0));
    }
    ccey *= adxinv[0] * adxinv[0];
#elif (AMREX_SPACEDIM == 2)
    amrex::ignore_unused(ez,coord);
    Real dxx = adxinv[0] * adxinv[0];
    Real dxy = adxinv[0] * adxinv[1];
    Real ccey =  ey(i  ,j  ,k  ) * dxx * Real(2.0)
        - dxx * (ey(i-1,j  ,k  ) +
                 ey(i+1,j  ,k  ))
        + dxy * (ex(i-1,j  ,k  )
               - ex(i  ,j  ,k  )
               - ex(i-1,j+1,k  )
               + ex(i  ,j+1,k  ));
#elif (AMREX_SPACEDIM == 3)
    Real dxx = adxinv[0] * adxinv[0];
    Real dzz = adxinv[2] * adxinv[2];
    Real dxy = adxinv[0] * adxinv[1];
    Real dyz = adxinv[1] * adxinv[2];
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
    Ay(i,j,k) = ccey + beta * ey(i,j,k);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void mlcurlcurl_adotx_z (int i, int j, int k, Array4<Real> const& Az,
                         Array4<Real const> const& ex,
                         Array4<Real const> const& ey,
                         Array4<Real const> const& ez,
                         Real beta, GpuArray<Real,AMREX_SPACEDIM> const& adxinv
#if (AMREX_SPACEDIM < 3)
                         , int coord
#endif
                        )
{
#if (AMREX_SPACEDIM == 1)
    amrex::ignore_unused(ex,ey);
    Real ccez;
    if (coord == 0) {
        ccez =  ez(i  ,j,k) * Real(2.0)
            -   ez(i-1,j,k)
            -   ez(i+1,j,k);
    } else if (coord == 1) {
        if (i == 0) {
            ccez = Real(4.0) * (ez(0,0,0) - ez(1,0,0));
        } else {
            ccez =  ez(i  ,j,k) * Real(2.0)
                +   ez(i-1,j,k) * (Real( 0.5)/Real(i) - Real(1.0))
                +   ez(i+1,j,k) * (Real(-0.5)/Real(i) - Real(1.0));
        }
    } else {
        ccez =  ez(i  ,j,k) * Real(2.0)
            +   ez(i-1,j,k) * (Real( 1.0)/Real(i) - Real(1.0))
            +   ez(i+1,j,k) * (Real(-1.0)/Real(i) - Real(1.0));
    }
    ccez *= adxinv[0] * adxinv[0];
#elif (AMREX_SPACEDIM == 2)
    amrex::ignore_unused(ex,ey,coord);
    Real dxx = adxinv[0] * adxinv[0];
    Real dyy = adxinv[1] * adxinv[1];
    Real ccez =  ez(i  ,j  ,k  ) * (dxx+dyy)*Real(2.0)
        - dxx * (ez(i-1,j  ,k  ) +
                 ez(i+1,j  ,k  ))
        - dyy * (ez(i  ,j-1,k  ) +
                 ez(i  ,j+1,k  ));
#elif (AMREX_SPACEDIM == 3)
    Real dxx = adxinv[0] * adxinv[0];
    Real dyy = adxinv[1] * adxinv[1];
    Real dxz = adxinv[0] * adxinv[2];
    Real dyz = adxinv[1] * adxinv[2];
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
    Az(i,j,k) = ccez + beta * ez(i,j,k);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void mlcurlcurl_adotx_x_alpha (int i, int j, int k, Array4<Real> const& Ax,
                               Array4<Real const> const& ex,
                               Array4<Real const> const& ey,
                               Array4<Real const> const& ez,
                               Array4<Real const> const& afy,
                               Array4<Real const> const& afz,
                               Real beta, GpuArray<Real,AMREX_SPACEDIM> const& dxinv)
{
#if (AMREX_SPACEDIM == 1)
    amrex::ignore_unused(ey,ez,afy,afz,dxinv);
    Ax(i,j,k) = beta * ex(i,j,k);
#elif (AMREX_SPACEDIM == 2)
    amrex::ignore_unused(ez,afy);

    Real const alpha_hi = afz(i,j,k);
    Real const alpha_lo = afz(i,j-1,k);

    Real const curlz_hi = dxinv[0] * (ey(i+1,j  ,k) - ey(i  ,j  ,k))
                        - dxinv[1] * (ex(i  ,j+1,k) - ex(i  ,j  ,k));
    Real const curlz_lo = dxinv[0] * (ey(i+1,j-1,k) - ey(i  ,j-1,k))
                        - dxinv[1] * (ex(i  ,j  ,k) - ex(i  ,j-1,k));

    Real const ccex = dxinv[1] * (alpha_hi * curlz_hi - alpha_lo * curlz_lo);

    Ax(i,j,k) = ccex + beta * ex(i,j,k);
#else
    Real const ay_lo = afy(i,j  ,k-1);
    Real const ay_hi = afy(i,j  ,k  );
    Real const az_lo = afz(i,j-1,k  );
    Real const az_hi = afz(i,j  ,k  );

    Real const curlz_hi = dxinv[0] * (ey(i+1,j  ,k  ) - ey(i  ,j  ,k  ))
                        - dxinv[1] * (ex(i  ,j+1,k  ) - ex(i  ,j  ,k  ));
    Real const curlz_lo = dxinv[0] * (ey(i+1,j-1,k  ) - ey(i  ,j-1,k  ))
                        - dxinv[1] * (ex(i  ,j  ,k  ) - ex(i  ,j-1,k  ));

    Real const curly_hi = dxinv[2] * (ex(i  ,j  ,k+1) - ex(i  ,j  ,k  ))
                        - dxinv[0] * (ez(i+1,j  ,k  ) - ez(i  ,j  ,k  ));
    Real const curly_lo = dxinv[2] * (ex(i  ,j  ,k  ) - ex(i  ,j  ,k-1))
                        - dxinv[0] * (ez(i+1,j  ,k-1) - ez(i  ,j  ,k-1));

    Real const ccex = dxinv[1] * (az_hi * curlz_hi - az_lo * curlz_lo)
        -             dxinv[2] * (ay_hi * curly_hi - ay_lo * curly_lo);

    Ax(i,j,k) = ccex + beta * ex(i,j,k);
#endif
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void mlcurlcurl_adotx_y_alpha (int i, int j, int k, Array4<Real> const& Ax,
                               Array4<Real const> const& ex,
                               Array4<Real const> const& ey,
                               Array4<Real const> const& ez,
                               Array4<Real const> const& afx,
                               Array4<Real const> const& afz,
                               Real beta, GpuArray<Real,AMREX_SPACEDIM> const& dxinv
#if (AMREX_SPACEDIM < 3)
                               , int coord
#endif
                               )
{
#if (AMREX_SPACEDIM == 1)
    amrex::ignore_unused(ex,ez,afx);
    Real const dxinv0 = dxinv[0];
    Real const alpha_hi = afz(i  ,j,k);
    Real const alpha_lo = afz(i-1,j,k);

    // The discrezation for curvilinear coordinates does not reduce to the
    // same form used by the scalar alpha case. But I think either form is fine.

    if (coord == 0) {
        Real const curl_lo = dxinv0 * (ey(i  ,j,k) - ey(i-1,j,k));
        Real const curl_hi = dxinv0 * (ey(i+1,j,k) - ey(i  ,j,k));
        Real const ccey = dxinv0 * (alpha_lo * curl_lo - alpha_hi * curl_hi);
        Ax(i,j,k) = ccey + beta * ey(i,j,k);
    } else if (coord == 1) {
        Real const curl_lo = dxinv0 * (Real(i)*ey(i  ,j,k) - Real(i-1)*ey(i-1,j,k))
            / (Real(i)-Real(0.5));
        Real const curl_hi = dxinv0 * (Real(i+1)*ey(i+1,j,k) - Real(i)*ey(i  ,j,k))
            / (Real(i)+Real(0.5));
        Real const ccey = dxinv0 * (alpha_lo * curl_lo - alpha_hi * curl_hi);
        Ax(i,j,k) = ccey + beta * ey(i,j,k);
    } else {
        Real const curl_lo = dxinv0 * (Real(i)*ey(i  ,j,k) - Real(i-1)*ey(i-1,j,k));
        Real const curl_hi = dxinv0 * (Real(i+1)*ey(i+1,j,k) - Real(i)*ey(i  ,j,k));
        Real const ccey = dxinv0 * (alpha_lo*curl_lo - alpha_hi*curl_hi) / Real(i);
        Ax(i,j,k) = ccey + beta * ey(i,j,k);
    }
#elif (AMREX_SPACEDIM == 2)
    amrex::ignore_unused(ez,afx,coord);

    Real const alpha_rt = afz(i  ,j,k);
    Real const alpha_lf = afz(i-1,j,k);

    Real const curlz_rt = dxinv[0] * (ey(i+1,j  ,k) - ey(i  ,j  ,k))
                        - dxinv[1] * (ex(i  ,j+1,k) - ex(i  ,j  ,k));
    Real const curlz_lf = dxinv[0] * (ey(i  ,j  ,k) - ey(i-1,j  ,k))
                        - dxinv[1] * (ex(i-1,j+1,k) - ex(i-1,j  ,k));

    Real const ccey = dxinv[0] * (alpha_lf * curlz_lf - alpha_rt * curlz_rt);

    Ax(i,j,k) = ccey + beta * ey(i,j,k);
#else
    Real const ax_lo = afx(i  ,j,k-1);
    Real const ax_hi = afx(i  ,j,k  );
    Real const az_lo = afz(i-1,j,k  );
    Real const az_hi = afz(i  ,j,k  );

    Real const curlz_hi = dxinv[0] * (ey(i+1,j  ,k  ) - ey(i  ,j  ,k  ))
                        - dxinv[1] * (ex(i  ,j+1,k  ) - ex(i  ,j  ,k  ));
    Real const curlz_lo = dxinv[0] * (ey(i  ,j  ,k  ) - ey(i-1,j  ,k  ))
                        - dxinv[1] * (ex(i-1,j+1,k  ) - ex(i-1,j  ,k  ));

    Real const curlx_hi = dxinv[1] * (ez(i  ,j+1,k  ) - ez(i  ,j  ,k  ))
                        - dxinv[2] * (ey(i  ,j  ,k+1) - ey(i  ,j  ,k  ));
    Real const curlx_lo = dxinv[1] * (ez(i  ,j+1,k-1) - ez(i  ,j  ,k-1))
                        - dxinv[2] * (ey(i  ,j  ,k  ) - ey(i  ,j  ,k-1));

    Real const ccey = dxinv[2] * (ax_hi * curlx_hi - ax_lo * curlx_lo)
        -             dxinv[0] * (az_hi * curlz_hi - az_lo * curlz_lo);

    Ax(i,j,k) = ccey + beta * ey(i,j,k);
#endif
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void mlcurlcurl_adotx_z_alpha (int i, int j, int k, Array4<Real> const& Ax,
                               Array4<Real const> const& ex,
                               Array4<Real const> const& ey,
                               Array4<Real const> const& ez,
                               Array4<Real const> const& afx,
                               Array4<Real const> const& afy,
                               Real beta, GpuArray<Real,AMREX_SPACEDIM> const& dxinv
#if (AMREX_SPACEDIM < 3)
                               , int coord
#endif
                               )
{
#if (AMREX_SPACEDIM == 1)
    amrex::ignore_unused(ex,ey,afx);
    Real const dxinv0 = dxinv[0];
    Real const alpha_hi = afy(i  ,j,k);
    Real const alpha_lo = afy(i-1,j,k);

    // The discrezation for curvilinear coordinates does not reduce to the
    // same form used by the scalar alpha case. But I think either form is fine.

    if (coord == 0) {
        Real const curl_lo = dxinv0 * (ez(i  ,j,k) - ez(i-1,j,k));
        Real const curl_hi = dxinv0 * (ez(i+1,j,k) - ez(i  ,j,k));
        Real const ccez = dxinv0 * (alpha_lo * curl_lo - alpha_hi * curl_hi);
        Ax(i,j,k) = ccez + beta * ez(i,j,k);
    } else if (coord == 1) {
        Real ccez;
        if (i == 0) {
            ccez = Real(4.0) * (ez(0,0,0) - ez(1,0,0)) * alpha_hi * dxinv0*dxinv0;
        } else {
            Real const curl_lo = dxinv0 * (ez(i  ,j,k) - ez(i-1,j,k));
            Real const curl_hi = dxinv0 * (ez(i+1,j,k) - ez(i  ,j,k));
            ccez = dxinv0 * ((Real(i)-Real(0.5))*alpha_lo*curl_lo -
                             (Real(i)+Real(0.5))*alpha_hi*curl_hi)/Real(i);
        }
        Ax(i,j,k) = ccez + beta * ez(i,j,k);
    } else {
        Real const curl_lo = dxinv0 * (Real(i)*ez(i  ,j,k) - Real(i-1)*ez(i-1,j,k));
        Real const curl_hi = dxinv0 * (Real(i+1)*ez(i+1,j,k) - Real(i)*ez(i  ,j,k));
        Real const ccez = dxinv0 * (alpha_lo*curl_lo - alpha_hi*curl_hi) / Real(i);
        Ax(i,j,k) = ccez + beta * ez(i,j,k);
    }
#elif (AMREX_SPACEDIM == 2)
    amrex::ignore_unused(ex,ey,coord);

    Real const cx_up =  dxinv[1] * (ez(i  ,j+1,k) - ez(i  ,j  ,k));
    Real const cx_dn =  dxinv[1] * (ez(i  ,j  ,k) - ez(i  ,j-1,k));
    Real const cy_rt = -dxinv[0] * (ez(i+1,j  ,k) - ez(i  ,j  ,k));
    Real const cy_lf = -dxinv[0] * (ez(i  ,j  ,k) - ez(i-1,j  ,k));

    Real const ccez = dxinv[0] * (afy(i  ,j  ,k) * cy_rt - afy(i-1,j  ,k) * cy_lf)
                    - dxinv[1] * (afx(i  ,j  ,k) * cx_up - afx(i  ,j-1,k) * cx_dn);

    Ax(i,j,k) = ccez + beta * ez(i,j,k);
#else
    Real const ay_lo = afy(i-1,j  ,k);
    Real const ay_hi = afy(i  ,j  ,k);
    Real const ax_lo = afx(i  ,j-1,k);
    Real const ax_hi = afx(i  ,j  ,k);

    Real const curlx_hi = dxinv[1] * (ez(i  ,j+1,k  ) - ez(i  ,j  ,k  ))
                        - dxinv[2] * (ey(i  ,j  ,k+1) - ey(i  ,j  ,k  ));
    Real const curlx_lo = dxinv[1] * (ez(i  ,j  ,k  ) - ez(i  ,j-1,k  ))
                        - dxinv[2] * (ey(i  ,j-1,k+1) - ey(i  ,j-1,k  ));

    Real const curly_hi = dxinv[2] * (ex(i  ,j  ,k+1) - ex(i  ,j  ,k  ))
                        - dxinv[0] * (ez(i+1,j  ,k  ) - ez(i  ,j  ,k  ));
    Real const curly_lo = dxinv[2] * (ex(i-1,j  ,k+1) - ex(i-1,j  ,k  ))
                        - dxinv[0] * (ez(i  ,j  ,k  ) - ez(i-1,j  ,k  ));

    Real const ccez = dxinv[0] * (ay_hi * curly_hi - ay_lo * curly_lo)
        -             dxinv[1] * (ax_hi * curlx_hi - ax_lo * curlx_lo);

    Ax(i,j,k) = ccez + beta * ez(i,j,k);
#endif
}

#if (AMREX_SPACEDIM == 1)
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void mlcurlcurl_smooth_1d (int i, int j, int k,
                           Array4<Real> const& ex,
                           Array4<Real> const& ey,
                           Array4<Real> const& ez,
                           Array4<Real const> const& rhsx,
                           Array4<Real const> const& rhsy,
                           Array4<Real const> const& rhsz,
                           Real beta,
                           GpuArray<Real,AMREX_SPACEDIM> const& adxinv,
                           int color,
                           CurlCurlDirichletInfo const& dinfo, bool valid_x,
                           int coord)
{
    Real dxx = adxinv[0] * adxinv[0];

    int my_color = i%2;

    if (my_color == color)
    {
        if (valid_x) {
            ex(i,j,k) = rhsx(i,j,k) / beta;
        }

        if (coord == 0)
        {
            if (dinfo.is_dirichlet_node(i,j,k)) { return; }

            Real cm = -dxx;
            Real cp = -dxx;

            Real gamma = dxx * Real(2.0) + beta;

            Real ccey = cm * ey(i-1,j,k) + cp * ey(i+1,j,k);
            Real res_y = rhsy(i,j,k) - ( gamma * ey(i,j,k) + ccey );
            ey(i,j,k) += res_y/gamma;

            Real ccez = cm * ez(i-1,j,k) + cp * ez(i+1,j,k);
            Real res_z = rhsz(i,j,k) - ( gamma * ez(i,j,k) + ccez );
            ez(i,j,k) += res_z/gamma;
        }
        else if (coord == 1)
        {
            if ((i > 0) && ! dinfo.is_dirichlet_node(i,j,k)) {
                Real cm = dxx * (Real( 0.5)/Real(i) - Real(1.0));
                Real cp = dxx * (Real(-0.5)/Real(i) - Real(1.0));

                Real gamma = dxx * (Real(1.0)/(Real(i)*Real(i)) + Real(2.0)) + beta;

                Real ccey = cm * ey(i-1,j,k) + cp * ey(i+1,j,k);
                Real res_y = rhsy(i,j,k) - ( gamma * ey(i,j,k) + ccey );
                ey(i,j,k) += res_y/gamma;

                gamma = dxx * Real(2.0) + beta;

                Real ccez = cm * ez(i-1,j,k) + cp * ez(i+1,j,k);
                Real res_z = rhsz(i,j,k) - ( gamma * ez(i,j,k) + ccez );
                ez(i,j,k) += res_z/gamma;
            }
            if ((i == 0) || (i == 1)) {
                ez(0,0,0) = (rhsz(0,0,0) + Real(4.0)*dxx*ez(1,0,0))
                    / (beta + Real(4.0)*dxx);
            }
        }
        else
        {
            if (i == 0 || dinfo.is_dirichlet_node(i,j,k)) { return; }

            Real cm = dxx * (Real( 1.0)/Real(i) - Real(1.0));
            Real cp = dxx * (Real(-1.0)/Real(i) - Real(1.0));

            Real gamma = dxx * Real(2.0) + beta;

            Real ccey = cm * ey(i-1,j,k) + cp * ey(i+1,j,k);
            Real res_y = rhsy(i,j,k) - ( gamma * ey(i,j,k) + ccey );
            ey(i,j,k) += res_y/gamma;

            Real ccez = cm * ez(i-1,j,k) + cp * ez(i+1,j,k);
            Real res_z = rhsz(i,j,k) - ( gamma * ez(i,j,k) + ccez );
            ez(i,j,k) += res_z/gamma;
        }
    }
}


AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void mlcurlcurl_smooth_1d (int i, int j, int k,
                           Array4<Real> const& ex,
                           Array4<Real> const& ey,
                           Array4<Real> const& ez,
                           Array4<Real const> const& rhsx,
                           Array4<Real const> const& rhsy,
                           Array4<Real const> const& rhsz,
                           Array4<Real const> const& betax,
                           Array4<Real const> const& betay,
                           Array4<Real const> const& betaz,
                           GpuArray<Real,AMREX_SPACEDIM> const& adxinv,
                           int color,
                           CurlCurlDirichletInfo const& dinfo, bool valid_x,
                           int coord)
{
    Real dxx = adxinv[0] * adxinv[0];

    int my_color = i%2;

    if (my_color == color)
    {
        if (valid_x) {
            ex(i,j,k) = rhsx(i,j,k) / betax(i,j,k);
        }

        if (coord == 0)
        {
            if (dinfo.is_dirichlet_node(i,j,k)) { return; }

            Real cm = -dxx;
            Real cp = -dxx;

            Real gamma_y = dxx * Real(2.0) + betay(i,j,k);
            Real ccey = cm * ey(i-1,j,k) + cp * ey(i+1,j,k);
            Real res_y = rhsy(i,j,k) - ( gamma_y * ey(i,j,k) + ccey );
            ey(i,j,k) += res_y/gamma_y;

            Real gamma_z = dxx * Real(2.0) + betaz(i,j,k);
            Real ccez = cm * ez(i-1,j,k) + cp * ez(i+1,j,k);
            Real res_z = rhsz(i,j,k) - ( gamma_z * ez(i,j,k) + ccez );
            ez(i,j,k) += res_z/gamma_z;
        }
        else if (coord == 1)
        {
            if ((i > 0) && ! dinfo.is_dirichlet_node(i,j,k)) {
                Real cm = dxx * (Real( 0.5)/Real(i) - Real(1.0));
                Real cp = dxx * (Real(-0.5)/Real(i) - Real(1.0));

                Real gamma = dxx * (Real(1.0)/(Real(i)*Real(i)) + Real(2.0)) + betay(i,j,k);

                Real ccey = cm * ey(i-1,j,k) + cp * ey(i+1,j,k);
                Real res_y = rhsy(i,j,k) - ( gamma * ey(i,j,k) + ccey );
                ey(i,j,k) += res_y/gamma;

                gamma = dxx * Real(2.0) + betaz(i,j,k);

                Real ccez = cm * ez(i-1,j,k) + cp * ez(i+1,j,k);
                Real res_z = rhsz(i,j,k) - ( gamma * ez(i,j,k) + ccez );
                ez(i,j,k) += res_z/gamma;
            }
            if ((i == 0) || (i == 1)) {
                ez(0,0,0) = (rhsz(0,0,0) + Real(4.0)*dxx*ez(1,0,0)) / (betaz(0,0,0) + Real(4.0)*dxx);
            }
        }
        else
        {
            if (dinfo.is_dirichlet_node(i,j,k)) { return; }

            Real cm = dxx * (Real( 1.0)/Real(i) - Real(1.0));
            Real cp = dxx * (Real(-1.0)/Real(i) - Real(1.0));

            Real gamma_y = dxx * Real(2.0) + betay(i,j,k);
            Real ccey = cm * ey(i-1,j,k) + cp * ey(i+1,j,k);
            Real res_y = rhsy(i,j,k) - ( gamma_y * ey(i,j,k) + ccey );
            ey(i,j,k) += res_y/gamma_y;

            Real gamma_z = dxx * Real(2.0) + betaz(i,j,k);
            Real ccez = cm * ez(i-1,j,k) + cp * ez(i+1,j,k);
            Real res_z = rhsz(i,j,k) - ( gamma_z * ez(i,j,k) + ccez );
            ez(i,j,k) += res_z/gamma_z;
        }
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void mlcurlcurl_smooth_1d_alpha_beta (int i, int j, int k,
                                      Array4<Real> const& ex,
                                      Array4<Real> const& ey,
                                      Array4<Real> const& ez,
                                      Array4<Real const> const& rhsx,
                                      Array4<Real const> const& rhsy,
                                      Array4<Real const> const& rhsz,
                                      Array4<Real const> const& betax,
                                      Array4<Real const> const& betay,
                                      Array4<Real const> const& betaz,
                                      GpuArray<Real,AMREX_SPACEDIM> const& adxinv,
                                      int color,
                                      CurlCurlDirichletInfo const& dinfo, bool valid_x,
                                      int coord,
                                      Array4<Real const> const& alphay,
                                      Array4<Real const> const& alphaz)
{
    Real dxx = adxinv[0] * adxinv[0];

    int my_color = i%2;
    if (my_color != color) { return; }

    if (valid_x) {
        ex(i,j,k) = rhsx(i,j,k) / betax(i,j,k);
    }

    auto cartesian_update = [&] (Array4<Real> const& comp,
                                 Array4<Real const> const& rhs,
                                 Array4<Real const> const& beta,
                                 Array4<Real const> const& alpha)
    {
        if (dinfo.is_dirichlet_node(i,j,k)) { return; }
        Real alpha_lo = alpha(i-1,j,k);
        Real alpha_hi = alpha(i  ,j,k);
        Real cm = -dxx * alpha_lo;
        Real cp = -dxx * alpha_hi;
        Real gamma = dxx * (alpha_lo + alpha_hi) + beta(i,j,k);
        Real cce = cm * comp(i-1,j,k) + cp * comp(i+1,j,k);
        Real res = rhs(i,j,k) - (gamma * comp(i,j,k) + cce);
        comp(i,j,k) += res / gamma;
    };

    auto cylindrical_y = [&] ()
    {
        if ((i > 0) && !dinfo.is_dirichlet_node(i,j,k)) {
            Real ri = Real(i);
            Real rimh = ri - Real(0.5);
            Real riph = ri + Real(0.5);
            Real alpha_lo = alphaz(i-1,j,k);
            Real alpha_hi = alphaz(i  ,j,k);
            Real cm = -dxx * alpha_lo * (ri-Real(1.0)) / rimh;
            Real cp = -dxx * alpha_hi * (ri+Real(1.0)) / riph;
            Real gamma = dxx * (alpha_lo * ri/rimh + alpha_hi * ri/riph)
                + betay(i,j,k);
            Real ccey = cm * ey(i-1,j,k) + cp * ey(i+1,j,k);
            Real res_y = rhsy(i,j,k) - (gamma * ey(i,j,k) + ccey);
            ey(i,j,k) += res_y / gamma;
        }
    };

    auto cylindrical_z = [&] ()
    {
        if ((i > 0) && !dinfo.is_dirichlet_node(i,j,k)) {
            Real ri = Real(i);
            Real alpha_lo = alphay(i-1,j,k);
            Real alpha_hi = alphay(i  ,j,k);
            Real cm = -(dxx/ri) * (ri-Real(0.5)) * alpha_lo;
            Real cp = -(dxx/ri) * (ri+Real(0.5)) * alpha_hi;
            Real gamma = (dxx/ri) * ((ri-Real(0.5))*alpha_lo +
                                     (ri+Real(0.5))*alpha_hi)
                + betaz(i,j,k);
            Real ccez = cm * ez(i-1,j,k) + cp * ez(i+1,j,k);
            Real res_z = rhsz(i,j,k) - (gamma * ez(i,j,k) + ccez);
            ez(i,j,k) += res_z / gamma;
        }
        if ((i == 0) || (i == 1)) {
            Real alpha_axis = alphay(0,j,k);
            Real gamma = betaz(0,j,k) + Real(4.0)*dxx*alpha_axis;
            Real num = rhsz(0,j,k) + Real(4.0)*dxx*alpha_axis * ez(1,j,k);
            ez(0,j,k) = num / gamma;
        }
    };

    auto spherical_update = [&] (Array4<Real> const& comp,
                                 Array4<Real const> const& rhs,
                                 Array4<Real const> const& beta,
                                 Array4<Real const> const& alpha)
    {
        if (i == 0 || dinfo.is_dirichlet_node(i,j,k)) { return; }
        Real ri = Real(i);
        Real alpha_lo = alpha(i-1,j,k);
        Real alpha_hi = alpha(i  ,j,k);
        Real cm = -dxx * alpha_lo * Real(i-1) / ri;
        Real cp = -dxx * alpha_hi * Real(i+1) / ri;
        Real gamma = dxx * (alpha_lo + alpha_hi) + beta(i,j,k);
        Real cce = cm * comp(i-1,j,k) + cp * comp(i+1,j,k);
        Real res = rhs(i,j,k) - (gamma * comp(i,j,k) + cce);
        comp(i,j,k) += res / gamma;
    };

    if (coord == 0) {
        cartesian_update(ey, rhsy, betay, alphaz);
        cartesian_update(ez, rhsz, betaz, alphay);
    } else if (coord == 1) {
        cylindrical_y();
        cylindrical_z();
    } else {
        spherical_update(ey, rhsy, betay, alphaz);
        spherical_update(ez, rhsz, betaz, alphay);
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void mlcurlcurl_smooth_1d_alpha (int i, int j, int k,
                                 Array4<Real> const& ex,
                                 Array4<Real> const& ey,
                                 Array4<Real> const& ez,
                                 Array4<Real const> const& rhsx,
                                 Array4<Real const> const& rhsy,
                                 Array4<Real const> const& rhsz,
                                 Real beta,
                                 GpuArray<Real,AMREX_SPACEDIM> const& adxinv,
                                 int color,
                                 CurlCurlDirichletInfo const& dinfo, bool valid_x,
                                 int coord,
                                 Array4<Real const> const& alphay,
                                 Array4<Real const> const& alphaz)
{
    Real dxx = adxinv[0] * adxinv[0];

    int my_color = i%2;
    if (my_color != color) { return; }

    if (valid_x) {
        ex(i,j,k) = rhsx(i,j,k) / beta;
    }

    auto cartesian_update = [&] (Array4<Real> const& comp,
                                 Array4<Real const> const& rhs,
                                 Array4<Real const> const& alpha)
    {
        if (dinfo.is_dirichlet_node(i,j,k)) { return; }
        Real alpha_lo = alpha(i-1,j,k);
        Real alpha_hi = alpha(i  ,j,k);
        Real cm = -dxx * alpha_lo;
        Real cp = -dxx * alpha_hi;
        Real gamma = dxx * (alpha_lo + alpha_hi) + beta;
        Real cce = cm * comp(i-1,j,k) + cp * comp(i+1,j,k);
        Real res = rhs(i,j,k) - (gamma * comp(i,j,k) + cce);
        comp(i,j,k) += res / gamma;
    };

    auto cylindrical_y = [&] ()
    {
        if ((i > 0) && !dinfo.is_dirichlet_node(i,j,k)) {
            Real ri = Real(i);
            Real rimh = ri - Real(0.5);
            Real riph = ri + Real(0.5);
            Real alpha_lo = alphaz(i-1,j,k);
            Real alpha_hi = alphaz(i  ,j,k);
            Real cm = -dxx * alpha_lo * (ri-Real(1.0)) / rimh;
            Real cp = -dxx * alpha_hi * (ri+Real(1.0)) / riph;
            Real gamma = dxx * (alpha_lo * ri/rimh + alpha_hi * ri/riph) + beta;
            Real ccey = cm * ey(i-1,j,k) + cp * ey(i+1,j,k);
            Real res_y = rhsy(i,j,k) - (gamma * ey(i,j,k) + ccey);
            ey(i,j,k) += res_y / gamma;
        }
    };

    auto cylindrical_z = [&] ()
    {
        if ((i > 0) && !dinfo.is_dirichlet_node(i,j,k)) {
            Real ri = Real(i);
            Real alpha_lo = alphay(i-1,j,k);
            Real alpha_hi = alphay(i  ,j,k);
            Real cm = -(dxx/ri) * (ri-Real(0.5)) * alpha_lo;
            Real cp = -(dxx/ri) * (ri+Real(0.5)) * alpha_hi;
            Real gamma = (dxx/ri) * ((ri-Real(0.5))*alpha_lo +
                                     (ri+Real(0.5))*alpha_hi) + beta;
            Real ccez = cm * ez(i-1,j,k) + cp * ez(i+1,j,k);
            Real res_z = rhsz(i,j,k) - (gamma * ez(i,j,k) + ccez);
            ez(i,j,k) += res_z / gamma;
        }
        if ((i == 0) || (i == 1)) {
            Real alpha_axis = alphay(0,j,k);
            Real gamma = beta + Real(4.0)*dxx*alpha_axis;
            Real num = rhsz(0,j,k) + Real(4.0)*dxx*alpha_axis * ez(1,j,k);
            ez(0,j,k) = num / gamma;
        }
    };

    auto spherical_update = [&] (Array4<Real> const& comp,
                                 Array4<Real const> const& rhs,
                                 Array4<Real const> const& alpha)
    {
        if (i == 0 || dinfo.is_dirichlet_node(i,j,k)) { return; }
        Real ri = Real(i);
        Real alpha_lo = alpha(i-1,j,k);
        Real alpha_hi = alpha(i  ,j,k);
        Real cm = -dxx * alpha_lo * Real(i-1) / ri;
        Real cp = -dxx * alpha_hi * Real(i+1) / ri;
        Real gamma = dxx * (alpha_lo + alpha_hi) + beta;
        Real cce = cm * comp(i-1,j,k) + cp * comp(i+1,j,k);
        Real res = rhs(i,j,k) - (gamma * comp(i,j,k) + cce);
        comp(i,j,k) += res / gamma;
    };

    if (coord == 0) {
        cartesian_update(ey, rhsy, alphaz);
        cartesian_update(ez, rhsz, alphay);
    } else if (coord == 1) {
        cylindrical_y();
        cylindrical_z();
    } else {
        spherical_update(ey, rhsy, alphaz);
        spherical_update(ez, rhsz, alphay);
    }
}


#endif

#if (AMREX_SPACEDIM > 1)

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void mlcurlcurl_gs4_lu (int i, int j, int k,
                        Array4<Real> const& ex,
                        Array4<Real> const& ey,
                        Array4<Real> const& ez,
                        Array4<Real const> const& rhsx,
                        Array4<Real const> const& rhsy,
                        Array4<Real const> const& rhsz,
#if (AMREX_SPACEDIM == 2)
                        Real beta,
#endif
                        GpuArray<Real,AMREX_SPACEDIM> const& adxinv,
                        int color, LUSolver<AMREX_SPACEDIM*2,Real> const& lusolver,
                        CurlCurlDirichletInfo const& dinfo,
                        CurlCurlSymmetryInfo const& sinfo)
{
    if (dinfo.is_dirichlet_node(i,j,k)) { return; }

    int my_color = i%2 + 2*(j%2);
    if (k%2 != 0) {
        my_color = 3 - my_color;
    }

#if (AMREX_SPACEDIM == 2)

    Real dxx = adxinv[0] * adxinv[0];
    Real dyy = adxinv[1] * adxinv[1];

    if (((my_color == 0 || my_color == 3) && (color == 0 || color == 3)) ||
        ((my_color == 1 || my_color == 2) && (color == 1 || color == 2)))
    {
        Real gamma = (dxx+dyy)*Real(2.0) + beta;
        Real ccez = - dxx * (ez(i-1,j  ,k  ) +
                             ez(i+1,j  ,k  ))
                    - dyy * (ez(i  ,j-1,k  ) +
                             ez(i  ,j+1,k  ));
        Real res = rhsz(i,j,k) - (gamma*ez(i,j,k) + ccez);
        constexpr Real omega = Real(1.15);
        ez(i,j,k) += omega/gamma * res;
    }

    if (my_color != color) { return; }

    Real dxy = adxinv[0]*adxinv[1];

    GpuArray<Real,4> b
        {rhsx(i-1,j,k) - (-dyy * ( ex(i-1,j-1,k  ) +
                                   ex(i-1,j+1,k  ))
                         + dxy * ( ey(i-1,j-1,k  )
                                  -ey(i-1,j  ,k  ))),
         rhsx(i  ,j,k) - (-dyy * ( ex(i  ,j-1,k  ) +
                                   ex(i  ,j+1,k  ))
                         + dxy * (-ey(i+1,j-1,k  )
                                  +ey(i+1,j  ,k  ))),
         rhsy(i,j-1,k) - (-dxx * ( ey(i-1,j-1,k  ) +
                                   ey(i+1,j-1,k  ))
                         + dxy * ( ex(i-1,j-1,k  )
                                  -ex(i  ,j-1,k  ))),
         rhsy(i,j  ,k) - (-dxx * ( ey(i-1,j  ,k  ) +
                                   ey(i+1,j  ,k  ))
                         + dxy * (-ex(i-1,j+1,k  )
                                  +ex(i  ,j+1,k  )))};

    if (sinfo.xlo_is_symmetric(i)) {
        b[0] = -b[1];
    } else if (sinfo.xhi_is_symmetric(i)) {
        b[1] = -b[0];
    }

    if (sinfo.ylo_is_symmetric(j)) {
        b[2] = -b[3];
    } else if (sinfo.yhi_is_symmetric(j)) {
        b[3] = -b[2];
    }

    GpuArray<Real,4> x;
    lusolver(x.data(), b.data());
    ex(i-1,j  ,k  ) = x[0];
    ex(i  ,j  ,k  ) = x[1];
    ey(i  ,j-1,k  ) = x[2];
    ey(i  ,j  ,k  ) = x[3];

#elif (AMREX_SPACEDIM == 3)

    if (my_color != color) { return; }

    Real dxx = adxinv[0]*adxinv[0];
    Real dyy = adxinv[1]*adxinv[1];
    Real dzz = adxinv[2]*adxinv[2];
    Real dxy = adxinv[0]*adxinv[1];
    Real dxz = adxinv[0]*adxinv[2];
    Real dyz = adxinv[1]*adxinv[2];

    GpuArray<Real,6> b
        {rhsx(i-1,j,k) - (-dyy * ( ex(i-1,j-1,k  ) +
                                   ex(i-1,j+1,k  ))
                         - dzz * ( ex(i-1,j  ,k+1) +
                                   ex(i-1,j  ,k-1))
                         + dxy * ( ey(i-1,j-1,k  )
                                  -ey(i-1,j  ,k  ))
                         + dxz * ( ez(i-1,j  ,k-1)
                                  -ez(i-1,j  ,k  ))),
         rhsx(i  ,j,k) - (-dyy * ( ex(i  ,j-1,k  ) +
                                   ex(i  ,j+1,k  ))
                         - dzz * ( ex(i  ,j  ,k+1) +
                                   ex(i  ,j  ,k-1))
                         + dxy * (-ey(i+1,j-1,k  )
                                  +ey(i+1,j  ,k  ))
                         + dxz * (-ez(i+1,j  ,k-1)
                                  +ez(i+1,j  ,k  ))),
         rhsy(i,j-1,k) - (-dxx * ( ey(i-1,j-1,k  ) +
                                   ey(i+1,j-1,k  ))
                         - dzz * ( ey(i  ,j-1,k-1) +
                                   ey(i  ,j-1,k+1))
                         + dxy * ( ex(i-1,j-1,k  )
                                  -ex(i  ,j-1,k  ))
                         + dyz * ( ez(i  ,j-1,k-1)
                                  -ez(i  ,j-1,k  ))),
         rhsy(i,j  ,k) - (-dxx * ( ey(i-1,j  ,k  ) +
                                   ey(i+1,j  ,k  ))
                         - dzz * ( ey(i  ,j  ,k-1) +
                                   ey(i  ,j  ,k+1))
                         + dxy * (-ex(i-1,j+1,k  )
                                  +ex(i  ,j+1,k  ))
                         + dyz * (-ez(i  ,j+1,k-1)
                                  +ez(i  ,j+1,k  ))),
         rhsz(i,j,k-1) - (-dxx * ( ez(i-1,j  ,k-1) +
                                   ez(i+1,j  ,k-1))
                         - dyy * ( ez(i  ,j-1,k-1) +
                                   ez(i  ,j+1,k-1))
                         + dxz * ( ex(i-1,j  ,k-1)
                                  -ex(i  ,j  ,k-1))
                         + dyz * ( ey(i  ,j-1,k-1)
                                  -ey(i  ,j  ,k-1))),
         rhsz(i,j,k  ) - (-dxx * ( ez(i-1,j  ,k  ) +
                                   ez(i+1,j  ,k  ))
                         - dyy * ( ez(i  ,j-1,k  ) +
                                   ez(i  ,j+1,k  ))
                         + dxz * (-ex(i-1,j  ,k+1)
                                  +ex(i  ,j  ,k+1))
                         + dyz * (-ey(i  ,j-1,k+1)
                                  +ey(i  ,j  ,k+1)))};

    if (sinfo.xlo_is_symmetric(i)) {
        b[0] = -b[1];
    } else if (sinfo.xhi_is_symmetric(i)) {
        b[1] = -b[0];
    }

    if (sinfo.ylo_is_symmetric(j)) {
        b[2] = -b[3];
    } else if (sinfo.yhi_is_symmetric(j)) {
        b[3] = -b[2];
    }

    if (sinfo.zlo_is_symmetric(k)) {
        b[4] = -b[5];
    } else if (sinfo.zhi_is_symmetric(k)) {
        b[5] = -b[4];
    }

    GpuArray<Real,6> x;
    lusolver(x.data(), b.data());
    ex(i-1,j  ,k  ) = x[0];
    ex(i  ,j  ,k  ) = x[1];
    ey(i  ,j-1,k  ) = x[2];
    ey(i  ,j  ,k  ) = x[3];
    ez(i  ,j  ,k-1) = x[4];
    ez(i  ,j  ,k  ) = x[5];
#endif
}

template <bool PCG>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void mlcurlcurl_gs4 (int i, int j, int k,
                     Array4<Real> const& ex,
                     Array4<Real> const& ey,
                     Array4<Real> const& ez,
                     Array4<Real const> const& rhsx,
                     Array4<Real const> const& rhsy,
                     Array4<Real const> const& rhsz,
                     GpuArray<Real,AMREX_SPACEDIM> const& adxinv,
                     int color,
                     Array4<Real const> const& betax,
                     Array4<Real const> const& betay,
                     Array4<Real const> const& betaz,
                     CurlCurlDirichletInfo const& dinfo,
                     CurlCurlSymmetryInfo const& sinfo)
{
    if (dinfo.is_dirichlet_node(i,j,k)) { return; }

    int my_color = i%2 + 2*(j%2);
    if (k%2 != 0) {
        my_color = 3 - my_color;
    }

#if (AMREX_SPACEDIM == 2)

    Real dxx = adxinv[0] * adxinv[0];
    Real dyy = adxinv[1] * adxinv[1];

    if (((my_color == 0 || my_color == 3) && (color == 0 || color == 3)) ||
        ((my_color == 1 || my_color == 2) && (color == 1 || color == 2)))
    {
        Real gamma = (dxx+dyy)*Real(2.0) + betaz(i,j,k);
        Real ccez = - dxx * (ez(i-1,j  ,k  ) +
                             ez(i+1,j  ,k  ))
                    - dyy * (ez(i  ,j-1,k  ) +
                             ez(i  ,j+1,k  ));
        Real res = rhsz(i,j,k) - (gamma*ez(i,j,k) + ccez);
        constexpr Real omega = Real(1.15);
        ez(i,j,k) += omega/gamma * res;
    }

    if (my_color != color) { return; }

    Real dxy = adxinv[0]*adxinv[1];

    GpuArray<Real,4> b
        {rhsx(i-1,j,k) - (-dyy * ( ex(i-1,j-1,k  ) +
                                   ex(i-1,j+1,k  ))
                         + dxy * ( ey(i-1,j-1,k  )
                                  -ey(i-1,j  ,k  ))),
         rhsx(i  ,j,k) - (-dyy * ( ex(i  ,j-1,k  ) +
                                   ex(i  ,j+1,k  ))
                         + dxy * (-ey(i+1,j-1,k  )
                                  +ey(i+1,j  ,k  ))),
         rhsy(i,j-1,k) - (-dxx * ( ey(i-1,j-1,k  ) +
                                   ey(i+1,j-1,k  ))
                         + dxy * ( ex(i-1,j-1,k  )
                                  -ex(i  ,j-1,k  ))),
         rhsy(i,j  ,k) - (-dxx * ( ey(i-1,j  ,k  ) +
                                   ey(i+1,j  ,k  ))
                         + dxy * (-ex(i-1,j+1,k  )
                                  +ex(i  ,j+1,k  )))};

    GpuArray<Real,4> beta;

    if (sinfo.xlo_is_symmetric(i)) {
        b[0] = -b[1];
        beta[0] = beta[1] = betax(i,j,k);
    } else if (sinfo.xhi_is_symmetric(i)) {
        b[1] = -b[0];
        beta[0] = beta[1] = betax(i-1,j,k);
    } else {
        beta[0] = betax(i-1,j,k);
        beta[1] = betax(i  ,j,k);
    }

    if (sinfo.ylo_is_symmetric(j)) {
        b[2] = -b[3];
        beta[2] = beta[3] = betay(i,j,k);
    } else if (sinfo.yhi_is_symmetric(j)) {
        b[3] = -b[2];
        beta[2] = beta[3] = betay(i,j-1,k);
    } else {
        beta[2] = betay(i,j-1,k);
        beta[3] = betay(i,j  ,k);
    }

    if constexpr (PCG) {
        Real diagInv[4] = {Real(1.0) / (dyy*Real(2.0) + beta[0]),
                           Real(1.0) / (dyy*Real(2.0) + beta[1]),
                           Real(1.0) / (dxx*Real(2.0) + beta[2]),
                           Real(1.0) / (dxx*Real(2.0) + beta[3])};
        auto precond = [&] (Real * AMREX_RESTRICT z,
                            Real const* AMREX_RESTRICT r)
        {
            for (int m = 0; m < 4; ++m) { z[m] = r[m] * diagInv[m]; }
        };
        auto mat = [&] (Real * AMREX_RESTRICT Av,
                        Real const* AMREX_RESTRICT v)
        {
            Av[0] = (dyy*Real(2.0) + beta[0]) * v[0] - dxy * v[2] + dxy * v[3];
            Av[1] = (dyy*Real(2.0) + beta[1]) * v[1] + dxy * v[2] - dxy * v[3];
            Av[2] = -dxy * v[0] + dxy * v[1] + (dxx*Real(2.0) + beta[2]) * v[2];
            Av[3] =  dxy * v[0] - dxy * v[1] + (dxx*Real(2.0) + beta[3]) * v[3];
        };
        Real sol[4] = {0, 0, 0, 0};
        pcg_solve<4>(sol, b.data(), mat, precond, 8, Real(1.e-8));
        ex(i-1,j  ,k  ) = sol[0];
        ex(i  ,j  ,k  ) = sol[1];
        ey(i  ,j-1,k  ) = sol[2];
        ey(i  ,j  ,k  ) = sol[3];
    } else {
        LUSolver<4,Real> lusolver
            ({dyy*Real(2.0) + beta[0],
              Real(0.0),
             -dxy,
              dxy,
              //
              Real(0.0),
              dyy*Real(2.0) + beta[1],
              dxy,
             -dxy,
              //
             -dxy,
              dxy,
              dxx*Real(2.0) + beta[2],
              Real(0.0),
              //
              dxy,
             -dxy,
              Real(0.0),
              dxx*Real(2.0) + beta[3]});
        lusolver(beta.data(), b.data());
        ex(i-1,j  ,k  ) = beta[0];
        ex(i  ,j  ,k  ) = beta[1];
        ey(i  ,j-1,k  ) = beta[2];
        ey(i  ,j  ,k  ) = beta[3];
    }

#elif (AMREX_SPACEDIM == 3)

    if (my_color != color) { return; }

    Real dxx = adxinv[0]*adxinv[0];
    Real dyy = adxinv[1]*adxinv[1];
    Real dzz = adxinv[2]*adxinv[2];
    Real dxy = adxinv[0]*adxinv[1];
    Real dxz = adxinv[0]*adxinv[2];
    Real dyz = adxinv[1]*adxinv[2];

    GpuArray<Real,6> b
        {rhsx(i-1,j,k) - (-dyy * ( ex(i-1,j-1,k  ) +
                                   ex(i-1,j+1,k  ))
                         - dzz * ( ex(i-1,j  ,k+1) +
                                   ex(i-1,j  ,k-1))
                         + dxy * ( ey(i-1,j-1,k  )
                                  -ey(i-1,j  ,k  ))
                         + dxz * ( ez(i-1,j  ,k-1)
                                  -ez(i-1,j  ,k  ))),
         rhsx(i  ,j,k) - (-dyy * ( ex(i  ,j-1,k  ) +
                                   ex(i  ,j+1,k  ))
                         - dzz * ( ex(i  ,j  ,k+1) +
                                   ex(i  ,j  ,k-1))
                         + dxy * (-ey(i+1,j-1,k  )
                                  +ey(i+1,j  ,k  ))
                         + dxz * (-ez(i+1,j  ,k-1)
                                  +ez(i+1,j  ,k  ))),
         rhsy(i,j-1,k) - (-dxx * ( ey(i-1,j-1,k  ) +
                                   ey(i+1,j-1,k  ))
                         - dzz * ( ey(i  ,j-1,k-1) +
                                   ey(i  ,j-1,k+1))
                         + dxy * ( ex(i-1,j-1,k  )
                                  -ex(i  ,j-1,k  ))
                         + dyz * ( ez(i  ,j-1,k-1)
                                  -ez(i  ,j-1,k  ))),
         rhsy(i,j  ,k) - (-dxx * ( ey(i-1,j  ,k  ) +
                                   ey(i+1,j  ,k  ))
                         - dzz * ( ey(i  ,j  ,k-1) +
                                   ey(i  ,j  ,k+1))
                         + dxy * (-ex(i-1,j+1,k  )
                                  +ex(i  ,j+1,k  ))
                         + dyz * (-ez(i  ,j+1,k-1)
                                  +ez(i  ,j+1,k  ))),
         rhsz(i,j,k-1) - (-dxx * ( ez(i-1,j  ,k-1) +
                                   ez(i+1,j  ,k-1))
                         - dyy * ( ez(i  ,j-1,k-1) +
                                   ez(i  ,j+1,k-1))
                         + dxz * ( ex(i-1,j  ,k-1)
                                  -ex(i  ,j  ,k-1))
                         + dyz * ( ey(i  ,j-1,k-1)
                                  -ey(i  ,j  ,k-1))),
         rhsz(i,j,k  ) - (-dxx * ( ez(i-1,j  ,k  ) +
                                   ez(i+1,j  ,k  ))
                         - dyy * ( ez(i  ,j-1,k  ) +
                                   ez(i  ,j+1,k  ))
                         + dxz * (-ex(i-1,j  ,k+1)
                                  +ex(i  ,j  ,k+1))
                         + dyz * (-ey(i  ,j-1,k+1)
                                  +ey(i  ,j  ,k+1)))};

    GpuArray<Real,6> beta;

    if (sinfo.xlo_is_symmetric(i)) {
        b[0] = -b[1];
        beta[0] = beta[1] = betax(i,j,k);
    } else if (sinfo.xhi_is_symmetric(i)) {
        b[1] = -b[0];
        beta[0] = beta[1] = betax(i-1,j,k);
    } else {
        beta[0] = betax(i-1,j,k);
        beta[1] = betax(i  ,j,k);
    }

    if (sinfo.ylo_is_symmetric(j)) {
        b[2] = -b[3];
        beta[2] = beta[3] = betay(i,j,k);
    } else if (sinfo.yhi_is_symmetric(j)) {
        b[3] = -b[2];
        beta[2] = beta[3] = betay(i,j-1,k);
    } else {
        beta[2] = betay(i,j-1,k);
        beta[3] = betay(i,j  ,k);
    }

    if (sinfo.zlo_is_symmetric(k)) {
        b[4] = -b[5];
        beta[4] = beta[5] = betaz(i,j,k);
    } else if (sinfo.zhi_is_symmetric(k)) {
        b[5] = -b[4];
        beta[4] = beta[5] = betaz(i,j,k-1);
    } else {
        beta[4] = betaz(i,j,k-1);
        beta[5] = betaz(i,j,k  );
    }

    if constexpr (PCG) {
        Real diagInv[6] = {Real(1.0) / ((dyy+dzz)*Real(2.0) + beta[0]),
                           Real(1.0) / ((dyy+dzz)*Real(2.0) + beta[1]),
                           Real(1.0) / ((dxx+dzz)*Real(2.0) + beta[2]),
                           Real(1.0) / ((dxx+dzz)*Real(2.0) + beta[3]),
                           Real(1.0) / ((dxx+dyy)*Real(2.0) + beta[4]),
                           Real(1.0) / ((dxx+dyy)*Real(2.0) + beta[5])};
        auto precond = [&] (Real * AMREX_RESTRICT z,
                            Real const* AMREX_RESTRICT r)
        {
            for (int m = 0; m < 6; ++m) { z[m] = r[m] * diagInv[m]; }
        };
        auto mat = [&] (Real * AMREX_RESTRICT Av,
                        Real const* AMREX_RESTRICT v)
        {
            Av[0] = ((dyy+dzz)*Real(2.0) + beta[0]) * v[0] - dxy * v[2]
                + dxy * v[3] - dxz * v[4] + dxz * v[5];
            Av[1] = ((dyy+dzz)*Real(2.0) + beta[1]) * v[1] + dxy * v[2]
                - dxy * v[3] + dxz * v[4] - dxz * v[5];
            Av[2] = -dxy * v[0] + dxy * v[1] + ((dxx+dzz)*Real(2.0) + beta[2]) * v[2]
                - dyz * v[4] + dyz * v[5];
            Av[3] = dxy * v[0] - dxy * v[1] + ((dxx+dzz)*Real(2.0) + beta[3]) * v[3]
                + dyz * v[4] - dyz * v[5];
            Av[4] = -dxz * v[0] + dxz * v[1] - dyz * v[2] + dyz * v[3]
                + ((dxx+dyy)*Real(2.0) + beta[4]) * v[4];
            Av[5] = dxz * v[0] - dxz * v[1] + dyz * v[2] - dyz * v[3]
                + ((dxx+dyy)*Real(2.0) + beta[5]) * v[5];
        };
        Real sol[6] = {0, 0, 0, 0, 0, 0};
        pcg_solve<6>(sol, b.data(), mat, precond, 8, Real(1.e-8));
        ex(i-1,j  ,k  ) = sol[0];
        ex(i  ,j  ,k  ) = sol[1];
        ey(i  ,j-1,k  ) = sol[2];
        ey(i  ,j  ,k  ) = sol[3];
        ez(i  ,j  ,k-1) = sol[4];
        ez(i  ,j  ,k  ) = sol[5];
    } else {
        LUSolver<6,Real> lusolver
            ({(dyy+dzz)*Real(2.0) + beta[0],
              Real(0.0),
             -dxy,
              dxy,
             -dxz,
              dxz,
              //
              Real(0.0),
              (dyy+dzz)*Real(2.0) + beta[1],
              dxy,
             -dxy,
              dxz,
             -dxz,
              //
             -dxy,
              dxy,
              (dxx+dzz)*Real(2.0) + beta[2],
              Real(0.0),
             -dyz,
              dyz,
              //
              dxy,
             -dxy,
              Real(0.0),
              (dxx+dzz)*Real(2.0) + beta[3],
              dyz,
             -dyz,
              //
             -dxz,
              dxz,
             -dyz,
              dyz,
              (dxx+dyy)*Real(2.0) + beta[4],
              Real(0.0),
              //
              dxz,
             -dxz,
              dyz,
             -dyz,
              Real(0.0),
              (dxx+dyy)*Real(2.0) + beta[5]});
        lusolver(beta.data(), b.data());
        ex(i-1,j  ,k  ) = beta[0];
        ex(i  ,j  ,k  ) = beta[1];
        ey(i  ,j-1,k  ) = beta[2];
        ey(i  ,j  ,k  ) = beta[3];
        ez(i  ,j  ,k-1) = beta[4];
        ez(i  ,j  ,k  ) = beta[5];
    }
#endif
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void mlcurlcurl_gs4_alpha (int i, int j, int k,
                           Array4<Real> const& ex,
                           Array4<Real> const& ey,
                           Array4<Real> const& ez,
                           Array4<Real const> const& rhsx,
                           Array4<Real const> const& rhsy,
                           Array4<Real const> const& rhsz,
                           GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                           int color,
                           Array4<Real const> const& alphax,
                           Array4<Real const> const& alphay,
                           Array4<Real const> const& alphaz,
                           Array4<Real const> const& betax,
                           Array4<Real const> const& betay,
                           Array4<Real const> const& betaz, Real bscalar,
                           CurlCurlDirichletInfo const& dinfo,
                           CurlCurlSymmetryInfo const& sinfo)
{
    if (dinfo.is_dirichlet_node(i,j,k)) { return; }

    int my_color = i%2 + 2*(j%2);
    if (k%2 != 0) {
        my_color = 3 - my_color;
    }

#if (AMREX_SPACEDIM == 2)

    Real dxx = dxinv[0] * dxinv[0];
    Real dyy = dxinv[1] * dxinv[1];

    if (((my_color == 0 || my_color == 3) && (color == 0 || color == 3)) ||
        ((my_color == 1 || my_color == 2) && (color == 1 || color == 2)))
    {
        Real beta = betaz ? betaz(i,j,k) : bscalar;
        Real gamma = dxx*(alphay(i-1,j  ,k)+alphay(i,j,k))
            +        dyy*(alphax(i  ,j-1,k)+alphax(i,j,k)) + beta;
        Real ccez = - dxx * (ez(i-1,j  ,k  ) * alphay(i-1,j  ,k) +
                             ez(i+1,j  ,k  ) * alphay(i  ,j  ,k))
                    - dyy * (ez(i  ,j-1,k  ) * alphax(i  ,j-1,k) +
                             ez(i  ,j+1,k  ) * alphax(i  ,j  ,k));
        Real res = rhsz(i,j,k) - (gamma*ez(i,j,k) + ccez);
        constexpr Real omega = Real(1.15);
        ez(i,j,k) += omega/gamma * res;
    }

    if (my_color != color) { return; }

    Real dxy = dxinv[0]*dxinv[1];

    Real amm = alphaz(i-1,j-1,k);
    Real am0 = alphaz(i-1,j  ,k);
    Real a0m = alphaz(i  ,j-1,k);
    Real a00 = alphaz(i  ,j  ,k);

    GpuArray<Real,4> b
        {rhsx(i-1,j,k) - (-dyy * ( ex(i-1,j-1,k  ) * amm +
                                   ex(i-1,j+1,k  ) * am0 )
                         + dxy * ( ey(i-1,j-1,k  ) * amm
                                  -ey(i-1,j  ,k  ) * am0 )),
         rhsx(i  ,j,k) - (-dyy * ( ex(i  ,j-1,k  ) * a0m +
                                   ex(i  ,j+1,k  ) * a00 )
                         + dxy * (-ey(i+1,j-1,k  ) * a0m
                                  +ey(i+1,j  ,k  ) * a00 )),
         rhsy(i,j-1,k) - (-dxx * ( ey(i-1,j-1,k  ) * amm +
                                   ey(i+1,j-1,k  ) * a0m )
                         + dxy * ( ex(i-1,j-1,k  ) * amm
                                  -ex(i  ,j-1,k  ) * a0m )),
         rhsy(i,j  ,k) - (-dxx * ( ey(i-1,j  ,k  ) * am0 +
                                   ey(i+1,j  ,k  ) * a00 )
                         + dxy * (-ex(i-1,j+1,k  ) * am0
                                  +ex(i  ,j+1,k  ) * a00))};

    GpuArray<Real,4> beta;

    Real bxm = bscalar;
    Real bx0 = bscalar;
    if (betax) {
        bxm = betax(i-1,j,k);
        bx0 = betax(i  ,j,k);
    }
    Real bym = bscalar;
    Real by0 = bscalar;
    if (betay) {
        bym = betay(i,j-1,k);
        by0 = betay(i,j  ,k);
    }

    if (sinfo.xlo_is_symmetric(i)) {
        b[0] = -b[1];
        beta[0] = beta[1] = bx0;
    } else if (sinfo.xhi_is_symmetric(i)) {
        b[1] = -b[0];
        beta[0] = beta[1] = bxm;
    } else {
        beta[0] = bxm;
        beta[1] = bx0;
    }

    if (sinfo.ylo_is_symmetric(j)) {
        b[2] = -b[3];
        beta[2] = beta[3] = by0;
    } else if (sinfo.yhi_is_symmetric(j)) {
        b[3] = -b[2];
        beta[2] = beta[3] = bym;
    } else {
        beta[2] = bym;
        beta[3] = by0;
    }

    Real diagInv[4] = {Real(1.0) / (dyy*(amm+am0) + beta[0]),
                       Real(1.0) / (dyy*(a0m+a00) + beta[1]),
                       Real(1.0) / (dxx*(amm+a0m) + beta[2]),
                       Real(1.0) / (dxx*(am0+a00) + beta[3])};
    auto precond = [&] (Real * AMREX_RESTRICT z,
                        Real const* AMREX_RESTRICT r)
    {
        for (int m = 0; m < 4; ++m) { z[m] = r[m] * diagInv[m]; }
    };
    auto mat = [&] (Real * AMREX_RESTRICT Av,
                    Real const* AMREX_RESTRICT v)
    {
        Av[0] = (dyy*(am0+amm) + beta[0]) * v[0] - dxy*amm * v[2] + dxy*am0 * v[3];
        Av[1] = (dyy*(a00+a0m) + beta[1]) * v[1] + dxy*a0m * v[2] - dxy*a00 * v[3];
        Av[2] = -dxy*amm * v[0] + dxy*a0m * v[1] + (dxx*(amm+a0m) + beta[2]) * v[2];
        Av[3] =  dxy*am0 * v[0] - dxy*a00 * v[1] + (dxx*(am0+a00) + beta[3]) * v[3];
    };
    Real sol[4] = {0, 0, 0, 0};
    pcg_solve<4>(sol, b.data(), mat, precond, 8, Real(1.e-8));
    ex(i-1,j  ,k  ) = sol[0];
    ex(i  ,j  ,k  ) = sol[1];
    ey(i  ,j-1,k  ) = sol[2];
    ey(i  ,j  ,k  ) = sol[3];

#elif (AMREX_SPACEDIM == 3)

    if (my_color != color) { return; }

    Real dxx = dxinv[0]*dxinv[0];
    Real dyy = dxinv[1]*dxinv[1];
    Real dzz = dxinv[2]*dxinv[2];
    Real dxy = dxinv[0]*dxinv[1];
    Real dxz = dxinv[0]*dxinv[2];
    Real dyz = dxinv[1]*dxinv[2];

    Real const az_mm = alphaz(i-1,j-1,k);
    Real const az_m0 = alphaz(i-1,j  ,k);
    Real const az_0m = alphaz(i  ,j-1,k);
    Real const az_00 = alphaz(i  ,j  ,k);

    Real const ay_mm = alphay(i-1,j  ,k-1);
    Real const ay_m0 = alphay(i-1,j  ,k  );
    Real const ay_0m = alphay(i  ,j  ,k-1);
    Real const ay_00 = alphay(i  ,j  ,k  );

    Real const ax_mm = alphax(i  ,j-1,k-1);
    Real const ax_m0 = alphax(i  ,j-1,k  );
    Real const ax_0m = alphax(i  ,j  ,k-1);
    Real const ax_00 = alphax(i  ,j  ,k  );

    GpuArray<Real,6> b
        {rhsx(i-1,j,k) - (-dyy * ( ex(i-1,j-1,k  ) * az_mm +
                                   ex(i-1,j+1,k  ) * az_m0)
                         - dzz * ( ex(i-1,j  ,k-1) * ay_mm +
                                   ex(i-1,j  ,k+1) * ay_m0)
                         + dxy * ( ey(i-1,j-1,k  ) * az_mm
                                  -ey(i-1,j  ,k  ) * az_m0)
                         + dxz * ( ez(i-1,j  ,k-1) * ay_mm
                                  -ez(i-1,j  ,k  ) * ay_m0)),
         rhsx(i  ,j,k) - (-dyy * ( ex(i  ,j-1,k  ) * az_0m +
                                   ex(i  ,j+1,k  ) * az_00)
                         - dzz * ( ex(i  ,j  ,k-1) * ay_0m +
                                   ex(i  ,j  ,k+1) * ay_00)
                         + dxy * (-ey(i+1,j-1,k  ) * az_0m
                                  +ey(i+1,j  ,k  ) * az_00)
                         + dxz * (-ez(i+1,j  ,k-1) * ay_0m
                                  +ez(i+1,j  ,k  ) * ay_00)),
         rhsy(i,j-1,k) - (-dxx * ( ey(i-1,j-1,k  ) * az_mm +
                                   ey(i+1,j-1,k  ) * az_0m)
                         - dzz * ( ey(i  ,j-1,k-1) * ax_mm +
                                   ey(i  ,j-1,k+1) * ax_m0)
                         + dxy * ( ex(i-1,j-1,k  ) * az_mm
                                  -ex(i  ,j-1,k  ) * az_0m)
                         + dyz * ( ez(i  ,j-1,k-1) * ax_mm
                                  -ez(i  ,j-1,k  ) * ax_m0)),
         rhsy(i,j  ,k) - (-dxx * ( ey(i-1,j  ,k  ) * az_m0 +
                                   ey(i+1,j  ,k  ) * az_00)
                         - dzz * ( ey(i  ,j  ,k-1) * ax_0m +
                                   ey(i  ,j  ,k+1) * ax_00)
                         + dxy * (-ex(i-1,j+1,k  ) * az_m0
                                  +ex(i  ,j+1,k  ) * az_00)
                         + dyz * (-ez(i  ,j+1,k-1) * ax_0m
                                  +ez(i  ,j+1,k  ) * ax_00)),
         rhsz(i,j,k-1) - (-dxx * ( ez(i-1,j  ,k-1) * ay_mm +
                                   ez(i+1,j  ,k-1) * ay_0m)
                         - dyy * ( ez(i  ,j-1,k-1) * ax_mm +
                                   ez(i  ,j+1,k-1) * ax_0m)
                         + dxz * ( ex(i-1,j  ,k-1) * ay_mm
                                  -ex(i  ,j  ,k-1) * ay_0m)
                         + dyz * ( ey(i  ,j-1,k-1) * ax_mm
                                  -ey(i  ,j  ,k-1) * ax_0m)),
         rhsz(i,j,k  ) - (-dxx * ( ez(i-1,j  ,k  ) * ay_m0 +
                                   ez(i+1,j  ,k  ) * ay_00)
                         - dyy * ( ez(i  ,j-1,k  ) * ax_m0 +
                                   ez(i  ,j+1,k  ) * ax_00)
                         + dxz * (-ex(i-1,j  ,k+1) * ay_m0
                                  +ex(i  ,j  ,k+1) * ay_00)
                         + dyz * (-ey(i  ,j-1,k+1) * ax_m0
                                  +ey(i  ,j  ,k+1) * ax_00))};

    auto get_bx = [&] (int ii, int jj, int kk) -> Real
    {
        return betax ? betax(ii,jj,kk) : bscalar;
    };
    auto get_by = [&] (int ii, int jj, int kk) -> Real
    {
        return betay ? betay(ii,jj,kk) : bscalar;
    };
    auto get_bz = [&] (int ii, int jj, int kk) -> Real
    {
        return betaz ? betaz(ii,jj,kk) : bscalar;
    };

    GpuArray<Real,6> beta;

    if (sinfo.xlo_is_symmetric(i)) {
        b[0] = -b[1];
        Real bx = get_bx(i,j,k);
        beta[0] = beta[1] = bx;
    } else if (sinfo.xhi_is_symmetric(i)) {
        b[1] = -b[0];
        Real bx = get_bx(i-1,j,k);
        beta[0] = beta[1] = bx;
    } else {
        beta[0] = get_bx(i-1,j,k);
        beta[1] = get_bx(i  ,j,k);
    }

    if (sinfo.ylo_is_symmetric(j)) {
        b[2] = -b[3];
        Real by = get_by(i,j,k);
        beta[2] = beta[3] = by;
    } else if (sinfo.yhi_is_symmetric(j)) {
        b[3] = -b[2];
        Real by = get_by(i,j-1,k);
        beta[2] = beta[3] = by;
    } else {
        beta[2] = get_by(i,j-1,k);
        beta[3] = get_by(i,j  ,k);
    }

    if (sinfo.zlo_is_symmetric(k)) {
        b[4] = -b[5];
        Real bz = get_bz(i,j,k);
        beta[4] = beta[5] = bz;
    } else if (sinfo.zhi_is_symmetric(k)) {
        b[5] = -b[4];
        Real bz = get_bz(i,j,k-1);
        beta[4] = beta[5] = bz;
    } else {
        beta[4] = get_bz(i,j,k-1);
        beta[5] = get_bz(i,j,k  );
    }

    Real diag0 = dyy*(az_mm+az_m0) + dzz*(ay_mm+ay_m0) + beta[0];
    Real diag1 = dyy*(az_0m+az_00) + dzz*(ay_0m+ay_00) + beta[1];
    Real diag2 = dxx*(az_mm+az_0m) + dzz*(ax_mm+ax_m0) + beta[2];
    Real diag3 = dxx*(az_m0+az_00) + dzz*(ax_0m+ax_00) + beta[3];
    Real diag4 = dxx*(ay_mm+ay_0m) + dyy*(ax_mm+ax_0m) + beta[4];
    Real diag5 = dxx*(ay_m0+ay_00) + dyy*(ax_m0+ax_00) + beta[5];

    Real diagInv[6] = {Real(1.0) / diag0,
                       Real(1.0) / diag1,
                       Real(1.0) / diag2,
                       Real(1.0) / diag3,
                       Real(1.0) / diag4,
                       Real(1.0) / diag5};
    auto precond = [&] (Real * AMREX_RESTRICT z,
                        Real const* AMREX_RESTRICT r)
    {
        for (int m = 0; m < 6; ++m) { z[m] = r[m] * diagInv[m]; }
    };
    auto mat = [&] (Real * AMREX_RESTRICT Av,
                    Real const* AMREX_RESTRICT v)
    {
        Av[0] = diag0 * v[0]
            - dxy * az_mm * v[2] + dxy * az_m0 * v[3]
            - dxz * ay_mm * v[4] + dxz * ay_m0 * v[5];
        Av[1] = diag1 * v[1]
            + dxy * az_0m * v[2] - dxy * az_00 * v[3]
            + dxz * ay_0m * v[4] - dxz * ay_00 * v[5];
        Av[2] = diag2 * v[2]
            - dxy * az_mm * v[0] + dxy * az_0m * v[1]
            - dyz * ax_mm * v[4] + dyz * ax_m0 * v[5];
        Av[3] = diag3 * v[3]
            + dxy * az_m0 * v[0] - dxy * az_00 * v[1]
            + dyz * ax_0m * v[4] - dyz * ax_00 * v[5];
        Av[4] = diag4 * v[4]
            - dxz * ay_mm * v[0] + dxz * ay_0m * v[1]
            - dyz * ax_mm * v[2] + dyz * ax_0m * v[3];
        Av[5] = diag5 * v[5]
            + dxz * ay_m0 * v[0] - dxz * ay_00 * v[1]
            + dyz * ax_m0 * v[2] - dyz * ax_00 * v[3];
    };
    Real sol[6] = {0, 0, 0, 0, 0, 0};
    pcg_solve<6>(sol, b.data(), mat, precond, 8, Real(1.e-8));
    ex(i-1,j  ,k  ) = sol[0];
    ex(i  ,j  ,k  ) = sol[1];
    ey(i  ,j-1,k  ) = sol[2];
    ey(i  ,j  ,k  ) = sol[3];
    ez(i  ,j  ,k-1) = sol[4];
    ez(i  ,j  ,k  ) = sol[5];
#endif
}

#endif

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
                             CurlCurlDirichletInfo const& dinfo)
{
    int ii = i*2;
    int jj = j*2;
    int kk = k*2;
    if (dinfo.is_dirichlet_edge(dir,ii,jj,kk)) {
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
#elif (AMREX_SPACEDIM == 2)
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

#elif (AMREX_SPACEDIM == 1)
        if (dir == 0) {
            crse(i,0,0) = Real(0.5) * (fine(ii,0,0) + fine(ii+1,0,0));
        } else if (i == 0 && dir == 2 && dinfo.coord == 1) {
            crse(0,0,0) = Real(0.5) * (fine(0,0,0) + fine(1,0,0));
        } else {
            crse(i,0,0) = Real(0.25) * (fine(ii-1,0,0)             +
                                        fine(ii  ,0,0) * Real(2.0) +
                                        fine(ii+1,0,0));
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
