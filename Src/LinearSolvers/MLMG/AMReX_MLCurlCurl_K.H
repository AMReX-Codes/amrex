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

struct CurlCurlDirichletInfo
{
    IntVect dirichlet_lo;
    IntVect dirichlet_hi;

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool is_dirichlet_node (int i, int j, int k) const
    {
#if (AMREX_SPACEDIM == 2)
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
    bool is_dirichlet_x_edge (int, int j, int k) const
    {
#if (AMREX_SPACEDIM == 2)
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
#if (AMREX_SPACEDIM == 2)
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
        return (i == dirichlet_lo[0]) || (i == dirichlet_hi[0])
            || (j == dirichlet_lo[1]) || (j == dirichlet_hi[1]);
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
#if (AMREX_SPACEDIM == 2)
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

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void mlcurlcurl_adotx_x (int i, int j, int k, Array4<Real> const& Ax,
                         Array4<Real const> const& ex,
                         Array4<Real const> const& ey,
                         Array4<Real const> const& ez,
                         Real beta, GpuArray<Real,AMREX_SPACEDIM> const& adxinv)
{
#if (AMREX_SPACEDIM == 2)
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
#else
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
                         Real beta, GpuArray<Real,AMREX_SPACEDIM> const& adxinv)
{
#if (AMREX_SPACEDIM == 2)
    amrex::ignore_unused(ez);
    Real dxx = adxinv[0] * adxinv[0];
    Real dxy = adxinv[0] * adxinv[1];
    Real ccey =  ey(i  ,j  ,k  ) * dxx * Real(2.0)
        - dxx * (ey(i-1,j  ,k  ) +
                 ey(i+1,j  ,k  ))
        + dxy * (ex(i-1,j  ,k  )
               - ex(i  ,j  ,k  )
               - ex(i-1,j+1,k  )
               + ex(i  ,j+1,k  ));
#else
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
                         Real beta, GpuArray<Real,AMREX_SPACEDIM> const& adxinv)
{
#if (AMREX_SPACEDIM == 2)
    amrex::ignore_unused(ex,ey);
    Real dxx = adxinv[0] * adxinv[0];
    Real dyy = adxinv[1] * adxinv[1];
    Real ccez =  ez(i  ,j  ,k  ) * (dxx+dyy)*Real(2.0)
        - dxx * (ez(i-1,j  ,k  ) +
                 ez(i+1,j  ,k  ))
        - dyy * (ez(i  ,j-1,k  ) +
                 ez(i  ,j+1,k  ));
#else
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

    GpuArray<Real,6> b
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

#else

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

    GpuArray<Real,6> b
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

#else

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
