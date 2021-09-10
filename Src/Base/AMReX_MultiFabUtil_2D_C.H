#ifndef AMREX_MULTIFAB_UTIL_2D_C_H_
#define AMREX_MULTIFAB_UTIL_2D_C_H_
#include <AMReX_Config.H>

#include <AMReX_Gpu.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_IArrayBox.H>
#include <cmath>

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void amrex_avg_nd_to_cc (int i, int j, int, int n,
                         Array4<Real      > const& cc,
                         Array4<Real const> const& nd,
                         int cccomp, int ndcomp) noexcept
{
    cc(i,j,0,n+cccomp) = Real(0.25)*( nd(i,j  ,0,n+ndcomp) + nd(i+1,j  ,0,n+ndcomp)
                                    + nd(i,j+1,0,n+ndcomp) + nd(i+1,j+1,0,n+ndcomp));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void amrex_avg_eg_to_cc (int i, int j, int,
                         Array4<Real      > const& cc,
                         Array4<Real const> const& Ex,
                         Array4<Real const> const& Ey,
                         int cccomp) noexcept
{
    cc(i,j,0,0+cccomp) = Real(0.5) * ( Ex(i,j,0) + Ex(i,j+1,0) );
    cc(i,j,0,1+cccomp) = Real(0.5) * ( Ey(i,j,0) + Ey(i+1,j,0) );
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void amrex_avg_fc_to_cc (int i, int j, int,
                         Array4<Real      > const& cc,
                         Array4<Real const> const& fx,
                         Array4<Real const> const& fy,
                         int cccomp) noexcept
{
    cc(i,j,0,0+cccomp) = Real(0.5) * ( fx(i,j,0) + fx(i+1,j,0) );
    cc(i,j,0,1+cccomp) = Real(0.5) * ( fy(i,j,0) + fy(i,j+1,0) );
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void amrex_avg_cc_to_fc (int i, int j, int, int n, Box const& xbx, Box const& ybx,
                         Array4<Real> const& fx, Array4<Real> const& fy,
                         Array4<Real const> const& cc, bool use_harmonic_averaging)  noexcept
{
    if (use_harmonic_averaging)
    {
        if (xbx.contains(i,j,0)) {
            if (cc(i-1,j,0,n) == Real(0.0) || cc(i,j,0,n) == Real(0.0)) {
                fx(i,j,0,n) = Real(0.0);
            } else {
                fx(i,j,0,n) = Real(2.0) / (Real(1.0) / cc(i-1,j,0,n) + Real(1.0) / cc(i,j,0,n));
            }
        }
        if (ybx.contains(i,j,0)) {
            if (cc(i,j-1,0,n) == Real(0.0) || cc(i,j,0,n) == Real(0.0)) {
                fy(i,j,0,n) = Real(0.0);
            } else {
                fy(i,j,0,n) = Real(2.0) / (Real(1.0) / cc(i,j-1,0,n) + Real(1.0) / cc(i,j,0,n));
            }
        }
    } else {
        if (xbx.contains(i,j,0)) {
            fx(i,j,0,n) = Real(0.5)*(cc(i-1,j,0,n) + cc(i,j,0,n));
        }
        if (ybx.contains(i,j,0)) {
            fy(i,j,0,n) = Real(0.5)*(cc(i,j-1,0,n) + cc(i,j,0,n));
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void amrex_avgdown_faces (Box const& bx, Array4<Real> const& crse,
                          Array4<Real const> const& fine,
                          int ccomp, int fcomp, int ncomp,
                          IntVect const& ratio, int idir) noexcept
{
    const auto clo = lbound(bx);
    const auto chi = ubound(bx);
    const int facx = ratio[0];
    const int facy = ratio[1];

    switch (idir) {
    case 0:
    {
        Real facInv =  Real(1.0) / facy;
        for (int n = 0; n < ncomp; ++n) {
            for (int j = clo.y; j <= chi.y; ++j) {
            for (int i = clo.x; i <= chi.x; ++i) {
                int ii = i*facx;
                int jj = j*facy;
                Real c = 0.;
                for (int jref = 0; jref < facy; ++jref) {
                    c += fine(ii,jj+jref,0,n+fcomp);
                }
                crse(i,j,0,n+ccomp) = c * facInv;
            }}
        }
        break;
    }
    case 1:
    {
        Real facInv = Real(1.0) / facx;
        for (int n = 0; n < ncomp; ++n) {
            for (int j = clo.y; j <= chi.y; ++j) {
            for (int i = clo.x; i <= chi.x; ++i) {
                int ii = i*facx;
                int jj = j*facy;
                Real c = 0.;
                for (int iref = 0; iref < facx; ++iref) {
                    c += fine(ii+iref,jj,0,n+fcomp);
                }
                crse(i,j,0,n+ccomp) = c * facInv;
            }}
        }
        break;
    }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void amrex_avgdown_faces (int i, int j, int, int n, Array4<Real> const& crse,
                          Array4<Real const> const& fine,
                          int ccomp, int fcomp, IntVect const& ratio, int idir) noexcept
{
    const int facx = ratio[0];
    const int facy = ratio[1];

    switch (idir) {
    case 0:
    {
        const Real facInv =  Real(1.0) / facy;
        const int ii = i*facx;
        const int jj = j*facy;
        Real c = 0.;
        for (int jref = 0; jref < facy; ++jref) {
            c += fine(ii,jj+jref,0,n+fcomp);
        }
        crse(i,j,0,n+ccomp) = c * facInv;
        break;
    }
    case 1:
    {
        const Real facInv = Real(1.0) / facx;
        const int ii = i*facx;
        const int jj = j*facy;
        Real c = 0.;
        for (int iref = 0; iref < facx; ++iref) {
            c += fine(ii+iref,jj,0,n+fcomp);
        }
        crse(i,j,0,n+ccomp) = c * facInv;
        break;
    }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void amrex_avgdown_edges (Box const& bx, Array4<Real> const& crse,
                          Array4<Real const> const& fine,
                          int ccomp, int fcomp, int ncomp,
                          IntVect const& ratio, int idir) noexcept
{
    const auto clo = lbound(bx);
    const auto chi = ubound(bx);
    const int facx = ratio[0];
    const int facy = ratio[1];

    switch (idir) {
    case 0:
    {
        Real facInv = Real(1.0) / facx;
        for (int n = 0; n < ncomp; ++n) {
            for (int j = clo.y; j <= chi.y; ++j) {
            for (int i = clo.x; i <= chi.x; ++i) {
                int ii = i*facx;
                int jj = j*facy;
                Real c = 0.;
                for (int iref = 0; iref < facx; ++iref) {
                    c += fine(ii+iref,jj,0,n+fcomp);
                }
                crse(i,j,0,n+ccomp) = c * facInv;
            }}
        }
        break;
    }
    case 1:
    {
        Real facInv = Real(1.0) / facy;
        for (int n = 0; n < ncomp; ++n) {
            for (int j = clo.y; j <= chi.y; ++j) {
            for (int i = clo.x; i <= chi.x; ++i) {
                int ii = i*facx;
                int jj = j*facy;
                Real c = 0.;
                for (int jref = 0; jref < facx; ++jref) {
                    c += fine(ii,jj+jref,0,n+fcomp);
                }
                crse(i,j,0,n+ccomp) = c * facInv;
            }}
        }
        break;
    }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void amrex_avgdown_edges (int i, int j, int, int n, Array4<Real> const& crse,
                          Array4<Real const> const& fine,
                          int ccomp, int fcomp, IntVect const& ratio, int idir) noexcept
{
    const int facx = ratio[0];
    const int facy = ratio[1];

    switch (idir) {
    case 0:
    {
        const Real facInv = Real(1.0) / facx;
        const int ii = i*facx;
        const int jj = j*facy;
        Real c = 0.;
        for (int iref = 0; iref < facx; ++iref) {
            c += fine(ii+iref,jj,0,n+fcomp);
        }
        crse(i,j,0,n+ccomp) = c * facInv;
        break;
    }
    case 1:
    {
        const Real facInv = Real(1.0) / facy;
        const int ii = i*facx;
        const int jj = j*facy;
        Real c = 0.;
        for (int jref = 0; jref < facx; ++jref) {
            c += fine(ii,jj+jref,0,n+fcomp);
        }
        crse(i,j,0,n+ccomp) = c * facInv;
        break;
    }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void amrex_avgdown (Box const& bx, Array4<Real> const& crse,
                    Array4<Real const> const& fine,
                    int ccomp, int fcomp, int ncomp,
                    IntVect const& ratio) noexcept
{
    const auto clo = lbound(bx);
    const auto chi = ubound(bx);
    const int facx = ratio[0];
    const int facy = ratio[1];
    const Real volfrac = static_cast<Real>(1.0)/static_cast<Real>(facx*facy);

    for (int n = 0; n < ncomp; ++n) {
        for (int j = clo.y; j <= chi.y; ++j) {
        for (int i = clo.x; i <= chi.x; ++i) {
            int ii = i*facx;
            int jj = j*facy;
            Real c = 0.;
            for (int jref = 0; jref < facy; ++jref) {
            for (int iref = 0; iref < facx; ++iref) {
                c += fine(ii+iref,jj+jref,0,n+fcomp);
            }}
            crse(i,j,0,n+ccomp) = volfrac * c;
        }}
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void amrex_avgdown (int i, int j, int, int n, Array4<Real> const& crse,
                    Array4<Real const> const& fine,
                    int ccomp, int fcomp, IntVect const& ratio) noexcept
{
    const int facx = ratio[0];
    const int facy = ratio[1];
    const Real volfrac = static_cast<Real>(1.0)/static_cast<Real>(facx*facy);
    const int ii = i*facx;
    const int jj = j*facy;
    Real c = Real(0.);
    for (int jref = 0; jref < facy; ++jref) {
    for (int iref = 0; iref < facx; ++iref) {
        c += fine(ii+iref,jj+jref,0,n+fcomp);
    }}
    crse(i,j,0,n+ccomp) = volfrac * c;
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void amrex_avgdown_nodes (Box const& bx, Array4<T> const& crse,
                          Array4<T const> const& fine,
                          int ccomp, int fcomp, int ncomp,
                          IntVect const& ratio) noexcept
{
    const auto clo = lbound(bx);
    const auto chi = ubound(bx);
    const int facx = ratio[0];
    const int facy = ratio[1];

    for (int n = 0; n < ncomp; ++n) {
        for     (int j = clo.y; j <= chi.y; ++j) {
            int jj = j*facy;
            AMREX_PRAGMA_SIMD
            for (int i = clo.x; i <= chi.x; ++i) {
                crse(i,j,0,n+ccomp) = fine(i*facx,jj,0,n+fcomp);
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void amrex_avgdown_nodes (int i, int j, int, int n, Array4<Real> const& crse,
                          Array4<Real const> const& fine,
                          int ccomp, int fcomp, IntVect const& ratio) noexcept
{
    crse(i,j,0,n+ccomp) = fine(i*ratio[0],j*ratio[1],0,n+fcomp);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void amrex_avgdown_with_vol (int i, int j, int, int n, Array4<Real> const& crse,
                             Array4<Real const> const& fine,
                             Array4<Real const> const& fv,
                             int ccomp, int fcomp, IntVect const& ratio) noexcept
{
    const int facx = ratio[0];
    const int facy = ratio[1];
    const int ii = i*facx;
    const int jj = j*facy;
    Real cd = Real(0.), cv = Real(0.);
    for (int jref = 0; jref < facy; ++jref) {
    for (int iref = 0; iref < facx; ++iref) {
        cv +=                                 fv(ii+iref,jj+jref,0);
        cd += fine(ii+iref,jj+jref,0,n+fcomp)*fv(ii+iref,jj+jref,0);
    }}
    crse(i,j,0,n+ccomp) = cd/cv;
}

/* All DG methods assume refinement ratio of 2 */
const Real VolumeRatio_2D = 1.0 / 4.0;
const Real Zero_2D = 0.0;

/* First order DG methods */
const int nDOFX_Order1_2D = 1;
const Real GaussLegendreWeights_Order1_2D[nDOFX_Order1_2D]
{
  +1.0000000000000000E+000
};
const Real ProjectionMatrix_T_Order1_2D[4][nDOFX_Order1_2D][nDOFX_Order1_2D]
{
  {
    {
      +1.0000000000000000E+000
    }
  },
  {
    {
      +1.0000000000000000E+000
    }
  },
  {
    {
      +1.0000000000000000E+000
    }
  },
  {
    {
      +1.0000000000000000E+000
    }
  }
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void amrex_avgdown_dg_order1
       ( int iCrse, int jCrse, int, int nComp,
         Array4<Real> const& crse,
         Array4<Real const> const& fine,
         Array4<Real const> const& fv,
         int ccomp, int fcomp, IntVect const& ratio ) noexcept
{

  amrex::ignore_unused( fv, fcomp, ccomp );

  const int nDOFX = nDOFX_Order1_2D;

  int iRefRatio = ratio[0];
  int jRefRatio = ratio[1];

  // Number of fields in MultiFab
  int nFields = nComp / nDOFX;

  // Index for projection matrix
  int iProj;

  // Offset between fine and coarse elements
  const int iOS = iCrse * iRefRatio;
  const int jOS = jCrse * jRefRatio;

  // Loop over fields
  for( int iField = 0; iField < nFields; iField++ )
  {

    // Loop over DOF of coarse element
    for( int iNX = 0; iNX < nDOFX; iNX++ )
    {

      // Initialize coarse data to zero
      crse(iCrse,jCrse,0,nDOFX*iField+iNX) = Zero_2D;

      iProj = -1;

      // Loop over fine cells in y-direction
      for( int jFine = 0; jFine < jRefRatio; jFine++ )
      {

        // Loop over fine cells in x-direction
        for( int iFine = 0; iFine < iRefRatio; iFine++ )
        {

          iProj += 1;

          // Project fine data onto coarse data
          for( int jNX = 0; jNX < nDOFX; jNX++ )
          {

            crse(iCrse,jCrse,0,nDOFX*iField+iNX)
              = crse(iCrse,jCrse,0,nDOFX*iField+iNX)
                  + ProjectionMatrix_T_Order1_2D[iProj][iNX][jNX]
                      * fine(iOS+iFine,jOS+jFine,0,nDOFX*iField+jNX);

          } // jNX

        } // iFine

      } // jFine

      crse(iCrse,jCrse,0,nDOFX*iField+iNX)
        = crse(iCrse,jCrse,0,nDOFX*iField+iNX)
            * VolumeRatio_2D / GaussLegendreWeights_Order1_2D[iNX];

    } // iNX
  } // iField

} // End of amrex_avgdown_dg_order1

/* Second order DG methods */
const int nDOFX_Order2_2D = 4;
const Real GaussLegendreWeights_Order2_2D[nDOFX_Order2_2D]
{
  +2.5000000000000000E-001,
  +2.5000000000000000E-001,
  +2.5000000000000000E-001,
  +2.5000000000000000E-001
};
const Real ProjectionMatrix_T_Order2_2D[4][nDOFX_Order2_2D][nDOFX_Order2_2D]
{
  {
    {
      +3.4987976320958231E-001,
      +2.0200317547305482E-001,
      +2.0200317547305482E-001,
      +1.1662658773652741E-001
    },
    {
      -5.4126587736527426E-002,
      +9.3750000000000000E-002,
      -3.1250000000000007E-002,
      +5.4126587736527405E-002
    },
    {
      -5.4126587736527426E-002,
      -3.1250000000000007E-002,
      +9.3750000000000000E-002,
      +5.4126587736527405E-002
    },
    {
      +8.3734122634725877E-003,
      -1.4503175473054832E-002,
      -1.4503175473054832E-002,
      +2.5120236790417749E-002
    }
  },
  {
    {
      +9.3750000000000000E-002,
      -5.4126587736527426E-002,
      +5.4126587736527405E-002,
      -3.1250000000000007E-002
    },
    {
      +2.0200317547305482E-001,
      +3.4987976320958231E-001,
      +1.1662658773652741E-001,
      +2.0200317547305482E-001
    },
    {
      -1.4503175473054832E-002,
      +8.3734122634725877E-003,
      +2.5120236790417749E-002,
      -1.4503175473054832E-002
    },
    {
      -3.1250000000000007E-002,
      -5.4126587736527426E-002,
      +5.4126587736527405E-002,
      +9.3750000000000000E-002
    }
  },
  {
    {
      +9.3750000000000000E-002,
      +5.4126587736527405E-002,
      -5.4126587736527426E-002,
      -3.1250000000000007E-002
    },
    {
      -1.4503175473054832E-002,
      +2.5120236790417749E-002,
      +8.3734122634725877E-003,
      -1.4503175473054832E-002
    },
    {
      +2.0200317547305482E-001,
      +1.1662658773652741E-001,
      +3.4987976320958231E-001,
      +2.0200317547305482E-001
    },
    {
      -3.1250000000000007E-002,
      +5.4126587736527405E-002,
      -5.4126587736527426E-002,
      +9.3750000000000000E-002
    }
  },
  {
    {
      +2.5120236790417749E-002,
      -1.4503175473054832E-002,
      -1.4503175473054832E-002,
      +8.3734122634725877E-003
    },
    {
      +5.4126587736527405E-002,
      +9.3750000000000000E-002,
      -3.1250000000000007E-002,
      -5.4126587736527426E-002
    },
    {
      +5.4126587736527405E-002,
      -3.1250000000000007E-002,
      +9.3750000000000000E-002,
      -5.4126587736527426E-002
    },
    {
      +1.1662658773652741E-001,
      +2.0200317547305482E-001,
      +2.0200317547305482E-001,
      +3.4987976320958231E-001
    }
  }
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void amrex_avgdown_dg_order2
       ( int iCrse, int jCrse, int, int nComp,
         Array4<Real> const& crse,
         Array4<Real const> const& fine,
         Array4<Real const> const& fv,
         int ccomp, int fcomp, IntVect const& ratio ) noexcept
{

  amrex::ignore_unused( fv, fcomp, ccomp );

  const int nDOFX = nDOFX_Order2_2D;

  int iRefRatio = ratio[0];
  int jRefRatio = ratio[1];

  // Number of fields in MultiFab
  int nFields = nComp / nDOFX;

  // Index for projection matrix
  int iProj;

  // Offset between fine and coarse elements
  const int iOS = iCrse * iRefRatio;
  const int jOS = jCrse * jRefRatio;

  // Loop over fields
  for( int iField = 0; iField < nFields; iField++ )
  {

    // Loop over DOF of coarse element
    for( int iNX = 0; iNX < nDOFX; iNX++ )
    {

      // Initialize coarse data to zero
      crse(iCrse,jCrse,0,nDOFX*iField+iNX) = Zero_2D;

      iProj = -1;

      // Loop over fine cells in y-direction
      for( int jFine = 0; jFine < jRefRatio; jFine++ )
      {

        // Loop over fine cells in x-direction
        for( int iFine = 0; iFine < iRefRatio; iFine++ )
        {

          iProj += 1;

          // Project fine data onto coarse data
          for( int jNX = 0; jNX < nDOFX; jNX++ )
          {

            crse(iCrse,jCrse,0,nDOFX*iField+iNX)
              = crse(iCrse,jCrse,0,nDOFX*iField+iNX)
                  + ProjectionMatrix_T_Order2_2D[iProj][iNX][jNX]
                      * fine(iOS+iFine,jOS+jFine,0,nDOFX*iField+jNX);

          } // jNX

        } // iFine

      } // jFine

      crse(iCrse,jCrse,0,nDOFX*iField+iNX)
        = crse(iCrse,jCrse,0,nDOFX*iField+iNX)
            * VolumeRatio_2D / GaussLegendreWeights_Order2_2D[iNX];

    } // iNX
  } // iField

} // End of amrex_avgdown_dg_order2

/* Third order DG methods */
const int nDOFX_Order3_2D = 9;
const Real GaussLegendreWeights_Order3_2D[nDOFX_Order3_2D]
{
  +7.7160493827160503E-002,
  +1.2345679012345678E-001,
  +7.7160493827160503E-002,
  +1.2345679012345678E-001,
  +1.9753086419753085E-001,
  +1.2345679012345678E-001,
  +7.7160493827160503E-002,
  +1.2345679012345678E-001,
  +7.7160493827160503E-002
};
const Real ProjectionMatrix_T_Order3_2D[4][nDOFX_Order3_2D][nDOFX_Order3_2D]
{
  {
    {
      +1.1651423916206255E-001,
      +8.0569101641365715E-002,
      +7.9014310551776996E-003,
      +8.0569101641365715E-002,
      +5.5713191674947982E-002,
      +5.4638060238403757E-003,
      +7.9014310551776996E-003,
      +5.4638060238403757E-003,
      +5.3583676268861434E-004
    },
    {
      -2.9598497555107860E-002,
      +8.8496027817990250E-002,
      +9.2809945996529461E-002,
      -2.0467235379121582E-002,
      +6.1194627467012229E-002,
      +6.4177683569844601E-002,
      -2.0072266656029412E-003,
      +6.0013717421124818E-003,
      +6.2939207671118564E-003
    },
    {
      +7.9014310551776996E-003,
      -1.7357653199944128E-002,
      -5.8942043895747588E-003,
      +5.4638060238403757E-003,
      -1.2002743484224964E-002,
      -4.0758173076004501E-003,
      +5.3583676268861434E-004,
      -1.1771119223314612E-003,
      -3.9971637753709788E-004
    },
    {
      -2.9598497555107857E-002,
      -2.0467235379121578E-002,
      -2.0072266656029412E-003,
      +8.8496027817990236E-002,
      +6.1194627467012229E-002,
      +6.0013717421124810E-003,
      +9.2809945996529461E-002,
      +6.4177683569844587E-002,
      +6.2939207671118564E-003
    },
    {
      +7.5190042334755056E-003,
      -2.2480938654752945E-002,
      -2.3576817558299035E-002,
      -2.2480938654752942E-002,
      +6.7215363511659812E-002,
      +7.0491912591652800E-002,
      -2.3576817558299039E-002,
      +7.0491912591652814E-002,
      +7.3928183695193897E-002
    },
    {
      -2.0072266656029412E-003,
      +4.4094220542980509E-003,
      +1.4973242366666817E-003,
      +6.0013717421124810E-003,
      -1.3183653530112371E-002,
      -4.4768234284154971E-003,
      +6.2939207671118564E-003,
      -1.3826317432949743E-002,
      -4.6950552569634653E-003
    },
    {
      +7.9014310551776996E-003,
      +5.4638060238403757E-003,
      +5.3583676268861434E-004,
      -1.7357653199944124E-002,
      -1.2002743484224964E-002,
      -1.1771119223314612E-003,
      -5.8942043895747588E-003,
      -4.0758173076004501E-003,
      -3.9971637753709788E-004
    },
    {
      -2.0072266656029412E-003,
      +6.0013717421124818E-003,
      +6.2939207671118564E-003,
      +4.4094220542980509E-003,
      -1.3183653530112371E-002,
      -1.3826317432949743E-002,
      +1.4973242366666817E-003,
      -4.4768234284154971E-003,
      -4.6950552569634653E-003
    },
    {
      +5.3583676268861434E-004,
      -1.1771119223314612E-003,
      -3.9971637753709788E-004,
      -1.1771119223314612E-003,
      +2.5858481055732723E-003,
      +8.7808628730366700E-004,
      -3.9971637753709788E-004,
      +8.7808628730366700E-004,
      +2.9817510405538792E-004
    }
  },
  {
    {
      -5.8942043895747597E-003,
      -1.7357653199944128E-002,
      +7.9014310551776979E-003,
      -4.0758173076004510E-003,
      -1.2002743484224964E-002,
      +5.4638060238403749E-003,
      -3.9971637753709799E-004,
      -1.1771119223314612E-003,
      +5.3583676268861434E-004
    },
    {
      +9.2809945996529461E-002,
      +8.8496027817990250E-002,
      -2.9598497555107853E-002,
      +6.4177683569844601E-002,
      +6.1194627467012229E-002,
      -2.0467235379121578E-002,
      +6.2939207671118564E-003,
      +6.0013717421124818E-003,
      -2.0072266656029408E-003
    },
    {
      +7.9014310551776979E-003,
      +8.0569101641365715E-002,
      +1.1651423916206259E-001,
      +5.4638060238403749E-003,
      +5.5713191674947982E-002,
      +8.0569101641365728E-002,
      +5.3583676268861434E-004,
      +5.4638060238403757E-003,
      +7.9014310551777014E-003
    },
    {
      +1.4973242366666819E-003,
      +4.4094220542980509E-003,
      -2.0072266656029408E-003,
      -4.4768234284154971E-003,
      -1.3183653530112371E-002,
      +6.0013717421124801E-003,
      -4.6950552569634662E-003,
      -1.3826317432949743E-002,
      +6.2939207671118555E-003
    },
    {
      -2.3576817558299035E-002,
      -2.2480938654752945E-002,
      +7.5190042334755038E-003,
      +7.0491912591652800E-002,
      +6.7215363511659812E-002,
      -2.2480938654752938E-002,
      +7.3928183695193897E-002,
      +7.0491912591652814E-002,
      -2.3576817558299032E-002
    },
    {
      -2.0072266656029408E-003,
      -2.0467235379121578E-002,
      -2.9598497555107864E-002,
      +6.0013717421124801E-003,
      +6.1194627467012229E-002,
      +8.8496027817990264E-002,
      +6.2939207671118555E-003,
      +6.4177683569844587E-002,
      +9.2809945996529489E-002
    },
    {
      -3.9971637753709799E-004,
      -1.1771119223314612E-003,
      +5.3583676268861434E-004,
      +8.7808628730366711E-004,
      +2.5858481055732723E-003,
      -1.1771119223314610E-003,
      +2.9817510405538797E-004,
      +8.7808628730366700E-004,
      -3.9971637753709782E-004
    },
    {
      +6.2939207671118564E-003,
      +6.0013717421124818E-003,
      -2.0072266656029408E-003,
      -1.3826317432949743E-002,
      -1.3183653530112371E-002,
      +4.4094220542980501E-003,
      -4.6950552569634653E-003,
      -4.4768234284154971E-003,
      +1.4973242366666812E-003
    },
    {
      +5.3583676268861434E-004,
      +5.4638060238403757E-003,
      +7.9014310551777014E-003,
      -1.1771119223314610E-003,
      -1.2002743484224964E-002,
      -1.7357653199944128E-002,
      -3.9971637753709782E-004,
      -4.0758173076004501E-003,
      -5.8942043895747606E-003
    }
  },
  {
    {
      -5.8942043895747588E-003,
      -4.0758173076004501E-003,
      -3.9971637753709793E-004,
      -1.7357653199944124E-002,
      -1.2002743484224964E-002,
      -1.1771119223314612E-003,
      +7.9014310551776979E-003,
      +5.4638060238403749E-003,
      +5.3583676268861423E-004
    },
    {
      +1.4973242366666819E-003,
      -4.4768234284154980E-003,
      -4.6950552569634662E-003,
      +4.4094220542980509E-003,
      -1.3183653530112371E-002,
      -1.3826317432949743E-002,
      -2.0072266656029408E-003,
      +6.0013717421124810E-003,
      +6.2939207671118555E-003
    },
    {
      -3.9971637753709793E-004,
      +8.7808628730366700E-004,
      +2.9817510405538792E-004,
      -1.1771119223314612E-003,
      +2.5858481055732723E-003,
      +8.7808628730366700E-004,
      +5.3583676268861423E-004,
      -1.1771119223314612E-003,
      -3.9971637753709782E-004
    },
    {
      +9.2809945996529461E-002,
      +6.4177683569844587E-002,
      +6.2939207671118564E-003,
      +8.8496027817990236E-002,
      +6.1194627467012229E-002,
      +6.0013717421124810E-003,
      -2.9598497555107850E-002,
      -2.0467235379121575E-002,
      -2.0072266656029408E-003
    },
    {
      -2.3576817558299039E-002,
      +7.0491912591652814E-002,
      +7.3928183695193897E-002,
      -2.2480938654752942E-002,
      +6.7215363511659812E-002,
      +7.0491912591652800E-002,
      +7.5190042334755047E-003,
      -2.2480938654752942E-002,
      -2.3576817558299032E-002
    },
    {
      +6.2939207671118564E-003,
      -1.3826317432949743E-002,
      -4.6950552569634653E-003,
      +6.0013717421124810E-003,
      -1.3183653530112371E-002,
      -4.4768234284154971E-003,
      -2.0072266656029408E-003,
      +4.4094220542980501E-003,
      +1.4973242366666814E-003
    },
    {
      +7.9014310551776979E-003,
      +5.4638060238403749E-003,
      +5.3583676268861423E-004,
      +8.0569101641365715E-002,
      +5.5713191674947982E-002,
      +5.4638060238403757E-003,
      +1.1651423916206258E-001,
      +8.0569101641365728E-002,
      +7.9014310551777014E-003
    },
    {
      -2.0072266656029408E-003,
      +6.0013717421124810E-003,
      +6.2939207671118555E-003,
      -2.0467235379121582E-002,
      +6.1194627467012229E-002,
      +6.4177683569844601E-002,
      -2.9598497555107867E-002,
      +8.8496027817990264E-002,
      +9.2809945996529475E-002
    },
    {
      +5.3583676268861423E-004,
      -1.1771119223314612E-003,
      -3.9971637753709782E-004,
      +5.4638060238403757E-003,
      -1.2002743484224964E-002,
      -4.0758173076004501E-003,
      +7.9014310551777014E-003,
      -1.7357653199944128E-002,
      -5.8942043895747597E-003
    }
  },
  {
    {
      +2.9817510405538797E-004,
      +8.7808628730366700E-004,
      -3.9971637753709788E-004,
      +8.7808628730366711E-004,
      +2.5858481055732723E-003,
      -1.1771119223314610E-003,
      -3.9971637753709788E-004,
      -1.1771119223314612E-003,
      +5.3583676268861423E-004
    },
    {
      -4.6950552569634662E-003,
      -4.4768234284154980E-003,
      +1.4973242366666814E-003,
      -1.3826317432949743E-002,
      -1.3183653530112371E-002,
      +4.4094220542980501E-003,
      +6.2939207671118555E-003,
      +6.0013717421124810E-003,
      -2.0072266656029404E-003
    },
    {
      -3.9971637753709788E-004,
      -4.0758173076004501E-003,
      -5.8942043895747606E-003,
      -1.1771119223314610E-003,
      -1.2002743484224964E-002,
      -1.7357653199944128E-002,
      +5.3583676268861423E-004,
      +5.4638060238403749E-003,
      +7.9014310551776996E-003
    },
    {
      -4.6950552569634662E-003,
      -1.3826317432949743E-002,
      +6.2939207671118555E-003,
      -4.4768234284154971E-003,
      -1.3183653530112371E-002,
      +6.0013717421124801E-003,
      +1.4973242366666817E-003,
      +4.4094220542980501E-003,
      -2.0072266656029404E-003
    },
    {
      +7.3928183695193897E-002,
      +7.0491912591652814E-002,
      -2.3576817558299032E-002,
      +7.0491912591652800E-002,
      +6.7215363511659812E-002,
      -2.2480938654752938E-002,
      -2.3576817558299032E-002,
      -2.2480938654752942E-002,
      +7.5190042334755021E-003
    },
    {
      +6.2939207671118555E-003,
      +6.4177683569844587E-002,
      +9.2809945996529489E-002,
      +6.0013717421124801E-003,
      +6.1194627467012229E-002,
      +8.8496027817990264E-002,
      -2.0072266656029404E-003,
      -2.0467235379121575E-002,
      -2.9598497555107860E-002
    },
    {
      -3.9971637753709788E-004,
      -1.1771119223314612E-003,
      +5.3583676268861423E-004,
      -4.0758173076004510E-003,
      -1.2002743484224964E-002,
      +5.4638060238403749E-003,
      -5.8942043895747606E-003,
      -1.7357653199944128E-002,
      +7.9014310551776996E-003
    },
    {
      +6.2939207671118555E-003,
      +6.0013717421124810E-003,
      -2.0072266656029404E-003,
      +6.4177683569844601E-002,
      +6.1194627467012229E-002,
      -2.0467235379121578E-002,
      +9.2809945996529475E-002,
      +8.8496027817990264E-002,
      -2.9598497555107857E-002
    },
    {
      +5.3583676268861423E-004,
      +5.4638060238403749E-003,
      +7.9014310551776996E-003,
      +5.4638060238403749E-003,
      +5.5713191674947982E-002,
      +8.0569101641365728E-002,
      +7.9014310551776996E-003,
      +8.0569101641365728E-002,
      +1.1651423916206261E-001
    }
  }
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void amrex_avgdown_dg_order3
       ( int iCrse, int jCrse, int, int nComp,
         Array4<Real> const& crse,
         Array4<Real const> const& fine,
         Array4<Real const> const& fv,
         int ccomp, int fcomp, IntVect const& ratio ) noexcept
{

  amrex::ignore_unused( fv, fcomp, ccomp );

  const int nDOFX = nDOFX_Order3_2D;

  int iRefRatio = ratio[0];
  int jRefRatio = ratio[1];

  // Number of fields in MultiFab
  int nFields = nComp / nDOFX;

  // Index for projection matrix
  int iProj;

  // Offset between fine and coarse elements
  const int iOS = iCrse * iRefRatio;
  const int jOS = jCrse * jRefRatio;

  // Loop over fields
  for( int iField = 0; iField < nFields; iField++ )
  {

    // Loop over DOF of coarse element
    for( int iNX = 0; iNX < nDOFX; iNX++ )
    {

      // Initialize coarse data to zero
      crse(iCrse,jCrse,0,nDOFX*iField+iNX) = Zero_2D;

      iProj = -1;

      // Loop over fine cells in y-direction
      for( int jFine = 0; jFine < jRefRatio; jFine++ )
      {

        // Loop over fine cells in x-direction
        for( int iFine = 0; iFine < iRefRatio; iFine++ )
        {

          iProj += 1;

          // Project fine data onto coarse data
          for( int jNX = 0; jNX < nDOFX; jNX++ )
          {

            crse(iCrse,jCrse,0,nDOFX*iField+iNX)
              = crse(iCrse,jCrse,0,nDOFX*iField+iNX)
                  + ProjectionMatrix_T_Order3_2D[iProj][iNX][jNX]
                      * fine(iOS+iFine,jOS+jFine,0,nDOFX*iField+jNX);

          } // jNX

        } // iFine

      } // jFine

      crse(iCrse,jCrse,0,nDOFX*iField+iNX)
        = crse(iCrse,jCrse,0,nDOFX*iField+iNX)
            * VolumeRatio_2D / GaussLegendreWeights_Order3_2D[iNX];

    } // iNX
  } // iField

} // End of amrex_avgdown_dg_order3

AMREX_GPU_HOST_DEVICE
inline
void amrex_compute_divergence (Box const& bx, Array4<Real> const& divu,
                               Array4<Real const> const& u,
                               Array4<Real const> const& v,
                               GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    const auto lo = lbound(bx);
    const auto hi = ubound(bx);
    const Real dxi = dxinv[0];
    const Real dyi = dxinv[1];

    for     (int n = 0; n < divu.ncomp; ++n) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    divu(i,j,0,n) = dxi * (u(i+1,j,0,n)-u(i,j,0,n))
                        +           dyi * (v(i,j+1,0,n)-v(i,j,0,n));
                }
        }
    }
}

AMREX_GPU_HOST_DEVICE
inline
void amrex_compute_gradient (Box const& bx, Array4<Real> const& grad,
                             Array4<Real const> const& u,
                             Array4<Real const> const& v,
                             GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    const auto lo = lbound(bx);
    const auto hi = ubound(bx);
    const Real dxi = dxinv[0];
    const Real dyi = dxinv[1];

    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            grad(i,j,0,0) = dxi * (u(i+1,j,0)-u(i,j,0));
            grad(i,j,0,1) = dyi * (v(i,j+1,0)-v(i,j,0));
        }
    }
}

AMREX_GPU_HOST_DEVICE
inline
void amrex_compute_convective_difference (Box const& bx, Array4<amrex::Real> const& diff,
                                          Array4<Real const> const& u_face,
                                          Array4<Real const> const& v_face,
                                          Array4<Real const> const& s_on_x_face,
                                          Array4<Real const> const& s_on_y_face,
                                          GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    const auto lo = lbound(bx);
    const auto hi = ubound(bx);
    const Real dxi = dxinv[0];
    const Real dyi = dxinv[1];

    for         (int n = 0; n <  diff.ncomp; ++n) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    diff(i,j,n) = Real(0.5)*dxi * (u_face(i+1,j,0,0)+u_face(i,j,0,0)) *
                                            (s_on_x_face(i+1,j,0,n)-s_on_x_face(i,j,0,n))
                        +         Real(0.5)*dyi * (v_face(i,j+1,0,0)+v_face(i,j,0,0)) *
                                            (s_on_y_face(i,j+1,0,n)-s_on_y_face(i,j,0,n));
                }
        }
    }
}

AMREX_GPU_HOST_DEVICE
inline
void amrex_compute_divergence_rz (Box const& bx, Array4<Real> const& divu,
                                  Array4<Real const> const& u,
                                  Array4<Real const> const& v,
                                  Array4<Real const> const& ax,
                                  Array4<Real const> const& ay,
                                  Array4<Real const> const& vol) noexcept
{
    const auto lo = lbound(bx);
    const auto hi = ubound(bx);

    for     (int n = 0; n < divu.ncomp; ++n) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    divu(i,j,0,n) = (ax(i+1,j,0,0)*u(i+1,j,0,n)-ax(i,j,0,0)*u(i,j,0,n))
                        +           (ay(i,j+1,0,0)*v(i,j+1,0,n)-ay(i,j,0,0)*v(i,j,0,n));
                    divu(i,j,0,n) /= vol(i,j,0,0);
                }
        }
    }
}

AMREX_GPU_HOST_DEVICE
inline
void amrex_compute_gradient_rz (Box const& bx, Array4<Real> const& grad,
                                Array4<Real const> const& u,
                                Array4<Real const> const& v,
                                Array4<Real const> const& ax,
                                Array4<Real const> const& ay,
                                Array4<Real const> const& vol) noexcept
{
    const auto lo = lbound(bx);
    const auto hi = ubound(bx);

    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            grad(i,j,0,0) = (ax(i+1,j,0,0)*u(i+1,j,0)-ax(i,j,0,0)*u(i,j,0))/vol(i,j,0,0);
            grad(i,j,0,1) = (ay(i,j+1,0,0)*v(i,j+1,0)-ay(i,j,0,0)*v(i,j,0))/vol(i,j,0,0);
        }
    }
}

}

#endif
