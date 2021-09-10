#ifndef AMREX_INTERP_2D_C_H_
#define AMREX_INTERP_2D_C_H_
#include <AMReX_Config.H>

#include <AMReX_FArrayBox.H>
#include <AMReX_BCRec.H>
#include <AMReX_Vector.H>
#include <cmath>

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
pcinterp_interp (Box const& bx,
                 Array4<Real> const& fine, const int fcomp, const int ncomp,
                 Array4<Real const> const& crse, const int ccomp,
                 IntVect const& ratio) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    for (int n = 0; n < ncomp; ++n) {
        for (int j = lo.y; j <= hi.y; ++j) {
            const int jc = amrex::coarsen(j,ratio[1]);
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                const int ic = amrex::coarsen(i,ratio[0]);
                fine(i,j,0,n+fcomp) = crse(ic,jc,0,n+ccomp);
            }
        }
    }
}

const Real Zero_2D = 0.0;

/* First order DG methods */
const int nDOFX_Order1_2D = 1;
const Real GaussLegendreWeights_Order1_2D[nDOFX_Order1_2D] =
{
  +1.0000000000000000E+000
};
const Real ProjectionMatrix_Order1_2D[4][nDOFX_Order1_2D][nDOFX_Order1_2D] =
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

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
dginterp_interp_order1
  ( Box const& bx,
    Array4<Real> const& fine, const int fcomp, const int ncomp,
    Array4<Real const> const& crse, const int ccomp,
    IntVect const& ratio) noexcept
{

  amrex::ignore_unused( fcomp, ccomp );

  const int nDOFX = nDOFX_Order1_2D;

  // Get bounds on fine grid
  const auto loFine = amrex::lbound(bx);
  const auto hiFine = amrex::ubound(bx);

  // Number of fields in MultiFab
  int nFields = ncomp / nDOFX;

  // Index for projection matrix
  int iProj;

  // Loop over fields
  for( int iField = 0; iField < nFields; iField++ )
  {

    // Loop over fine cells in y-direction
    for( int jFine = loFine.y; jFine <= hiFine.y; jFine++ )
    {

      // Get coarse element corresponding to fine element jFine
      const int jCrse = amrex::coarsen(jFine,ratio[1]);
      AMREX_PRAGMA_SIMD

      // Loop over fine cells in x-direction
      for( int iFine = loFine.x; iFine <= hiFine.x; iFine++ )
      {

        if     ( ( iFine % 2 == 0 ) & ( jFine % 2 == 0 ) ) { iProj = 0; }
        else if( ( jFine % 2 == 0 )                      ) { iProj = 1; }
        else if( ( iFine % 2 == 0 )                      ) { iProj = 2; }
        else                                               { iProj = 3; }

        // Get coarse element corresponding to fine element iFine
        const int iCrse = amrex::coarsen(iFine,ratio[0]);

        // Loop over DOF of fine element
        for( int iNX = 0; iNX < nDOFX; iNX++ )
        {

          // Initialize fine data to zero
          fine(iFine,jFine,0,nDOFX*iField+iNX) = Zero_2D;

          // Project coarse data onto fine data
          for( int jNX = 0; jNX < nDOFX; jNX++ )
          {

            fine(iFine,jFine,0,nDOFX*iField+iNX)
              = fine(iFine,jFine,0,nDOFX*iField+iNX)
                  + ProjectionMatrix_Order1_2D[iProj][iNX][jNX]
                      * crse(iCrse,jCrse,0,nDOFX*iField+jNX);

          } // jNX

          fine(iFine,jFine,0,nDOFX*iField+iNX)
            = fine(iFine,jFine,0,nDOFX*iField+iNX)
                / GaussLegendreWeights_Order1_2D[iNX];

        } // iNX
      } // iFine
    } // jFine
  } // iField

} // End of function dginterp_interp_order1

/* Second order DG methods */
const int nDOFX_Order2_2D = 4;
const Real GaussLegendreWeights_Order2_2D[nDOFX_Order2_2D] =
{
  +2.5000000000000000E-001,
  +2.5000000000000000E-001,
  +2.5000000000000000E-001,
  +2.5000000000000000E-001
};
const Real ProjectionMatrix_Order2_2D[4][nDOFX_Order2_2D][nDOFX_Order2_2D] =
{
  {
    {
      +3.4987976320958231E-001,
      -5.4126587736527426E-002,
      -5.4126587736527426E-002,
      +8.3734122634725877E-003
    },
    {
      +2.0200317547305482E-001,
      +9.3750000000000000E-002,
      -3.1250000000000007E-002,
      -1.4503175473054832E-002
    },
    {
      +2.0200317547305482E-001,
      -3.1250000000000007E-002,
      +9.3750000000000000E-002,
      -1.4503175473054832E-002
    },
    {
      +1.1662658773652741E-001,
      +5.4126587736527405E-002,
      +5.4126587736527405E-002,
      +2.5120236790417749E-002
    }
  },
  {
    {
      +9.3750000000000000E-002,
      +2.0200317547305482E-001,
      -1.4503175473054832E-002,
      -3.1250000000000007E-002
    },
    {
      -5.4126587736527426E-002,
      +3.4987976320958231E-001,
      +8.3734122634725877E-003,
      -5.4126587736527426E-002
    },
    {
      +5.4126587736527405E-002,
      +1.1662658773652741E-001,
      +2.5120236790417749E-002,
      +5.4126587736527405E-002
    },
    {
      -3.1250000000000007E-002,
      +2.0200317547305482E-001,
      -1.4503175473054832E-002,
      +9.3750000000000000E-002
    }
  },
  {
    {
      +9.3750000000000000E-002,
      -1.4503175473054832E-002,
      +2.0200317547305482E-001,
      -3.1250000000000007E-002
    },
    {
      +5.4126587736527405E-002,
      +2.5120236790417749E-002,
      +1.1662658773652741E-001,
      +5.4126587736527405E-002
    },
    {
      -5.4126587736527426E-002,
      +8.3734122634725877E-003,
      +3.4987976320958231E-001,
      -5.4126587736527426E-002
    },
    {
      -3.1250000000000007E-002,
      -1.4503175473054832E-002,
      +2.0200317547305482E-001,
      +9.3750000000000000E-002
    }
  },
  {
    {
      +2.5120236790417749E-002,
      +5.4126587736527405E-002,
      +5.4126587736527405E-002,
      +1.1662658773652741E-001
    },
    {
      -1.4503175473054832E-002,
      +9.3750000000000000E-002,
      -3.1250000000000007E-002,
      +2.0200317547305482E-001
    },
    {
      -1.4503175473054832E-002,
      -3.1250000000000007E-002,
      +9.3750000000000000E-002,
      +2.0200317547305482E-001
    },
    {
      +8.3734122634725877E-003,
      -5.4126587736527426E-002,
      -5.4126587736527426E-002,
      +3.4987976320958231E-001
    }
  }
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
dginterp_interp_order2
  ( Box const& bx,
    Array4<Real> const& fine, const int fcomp, const int ncomp,
    Array4<Real const> const& crse, const int ccomp,
    IntVect const& ratio) noexcept
{

  amrex::ignore_unused( fcomp, ccomp );

  const int nDOFX = nDOFX_Order2_2D;

  // Get bounds on fine grid
  const auto loFine = amrex::lbound(bx);
  const auto hiFine = amrex::ubound(bx);

  // Number of fields in MultiFab
  int nFields = ncomp / nDOFX;

  // Index for projection matrix
  int iProj;

  // Loop over fields
  for( int iField = 0; iField < nFields; iField++ )
  {

    // Loop over fine cells in y-direction
    for( int jFine = loFine.y; jFine <= hiFine.y; jFine++ )
    {

      // Get coarse element corresponding to fine element jFine
      const int jCrse = amrex::coarsen(jFine,ratio[1]);
      AMREX_PRAGMA_SIMD

      // Loop over fine cells in x-direction
      for( int iFine = loFine.x; iFine <= hiFine.x; iFine++ )
      {

        if     ( ( iFine % 2 == 0 ) & ( jFine % 2 == 0 ) ) { iProj = 0; }
        else if( ( jFine % 2 == 0 )                      ) { iProj = 1; }
        else if( ( iFine % 2 == 0 )                      ) { iProj = 2; }
        else                                               { iProj = 3; }

        // Get coarse element corresponding to fine element iFine
        const int iCrse = amrex::coarsen(iFine,ratio[0]);

        // Loop over DOF of fine element
        for( int iNX = 0; iNX < nDOFX; iNX++ )
        {

          // Initialize fine data to zero
          fine(iFine,jFine,0,nDOFX*iField+iNX) = Zero_2D;

          // Project coarse data onto fine data
          for( int jNX = 0; jNX < nDOFX; jNX++ )
          {

            fine(iFine,jFine,0,nDOFX*iField+iNX)
              = fine(iFine,jFine,0,nDOFX*iField+iNX)
                  + ProjectionMatrix_Order2_2D[iProj][iNX][jNX]
                      * crse(iCrse,jCrse,0,nDOFX*iField+jNX);

          } // jNX

          fine(iFine,jFine,0,nDOFX*iField+iNX)
            = fine(iFine,jFine,0,nDOFX*iField+iNX)
                / GaussLegendreWeights_Order2_2D[iNX];

        } // iNX
      } // iFine
    } // jFine
  } // iField

} // End of function dginterp_interp_order2

/* Third order DG methods */
const int nDOFX_Order3_2D = 9;
const Real GaussLegendreWeights_Order3_2D[nDOFX_Order3_2D] =
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
const Real ProjectionMatrix_Order3_2D[4][nDOFX_Order3_2D][nDOFX_Order3_2D] =
{
  {
    {
      +1.1651423916206255E-001,
      -2.9598497555107860E-002,
      +7.9014310551776996E-003,
      -2.9598497555107857E-002,
      +7.5190042334755056E-003,
      -2.0072266656029412E-003,
      +7.9014310551776996E-003,
      -2.0072266656029412E-003,
      +5.3583676268861434E-004
    },
    {
      +8.0569101641365715E-002,
      +8.8496027817990250E-002,
      -1.7357653199944128E-002,
      -2.0467235379121578E-002,
      -2.2480938654752945E-002,
      +4.4094220542980509E-003,
      +5.4638060238403757E-003,
      +6.0013717421124818E-003,
      -1.1771119223314612E-003
    },
    {
      +7.9014310551776996E-003,
      +9.2809945996529461E-002,
      -5.8942043895747588E-003,
      -2.0072266656029412E-003,
      -2.3576817558299035E-002,
      +1.4973242366666817E-003,
      +5.3583676268861434E-004,
      +6.2939207671118564E-003,
      -3.9971637753709788E-004
    },
    {
      +8.0569101641365715E-002,
      -2.0467235379121582E-002,
      +5.4638060238403757E-003,
      +8.8496027817990236E-002,
      -2.2480938654752942E-002,
      +6.0013717421124810E-003,
      -1.7357653199944124E-002,
      +4.4094220542980509E-003,
      -1.1771119223314612E-003
    },
    {
      +5.5713191674947982E-002,
      +6.1194627467012229E-002,
      -1.2002743484224964E-002,
      +6.1194627467012229E-002,
      +6.7215363511659812E-002,
      -1.3183653530112371E-002,
      -1.2002743484224964E-002,
      -1.3183653530112371E-002,
      +2.5858481055732723E-003
    },
    {
      +5.4638060238403757E-003,
      +6.4177683569844601E-002,
      -4.0758173076004501E-003,
      +6.0013717421124810E-003,
      +7.0491912591652800E-002,
      -4.4768234284154971E-003,
      -1.1771119223314612E-003,
      -1.3826317432949743E-002,
      +8.7808628730366700E-004
    },
    {
      +7.9014310551776996E-003,
      -2.0072266656029412E-003,
      +5.3583676268861434E-004,
      +9.2809945996529461E-002,
      -2.3576817558299039E-002,
      +6.2939207671118564E-003,
      -5.8942043895747588E-003,
      +1.4973242366666817E-003,
      -3.9971637753709788E-004
    },
    {
      +5.4638060238403757E-003,
      +6.0013717421124818E-003,
      -1.1771119223314612E-003,
      +6.4177683569844587E-002,
      +7.0491912591652814E-002,
      -1.3826317432949743E-002,
      -4.0758173076004501E-003,
      -4.4768234284154971E-003,
      +8.7808628730366700E-004
    },
    {
      +5.3583676268861434E-004,
      +6.2939207671118564E-003,
      -3.9971637753709788E-004,
      +6.2939207671118564E-003,
      +7.3928183695193897E-002,
      -4.6950552569634653E-003,
      -3.9971637753709788E-004,
      -4.6950552569634653E-003,
      +2.9817510405538792E-004
    }
  },
  {
    {
      -5.8942043895747597E-003,
      +9.2809945996529461E-002,
      +7.9014310551776979E-003,
      +1.4973242366666819E-003,
      -2.3576817558299035E-002,
      -2.0072266656029408E-003,
      -3.9971637753709799E-004,
      +6.2939207671118564E-003,
      +5.3583676268861434E-004
    },
    {
      -1.7357653199944128E-002,
      +8.8496027817990250E-002,
      +8.0569101641365715E-002,
      +4.4094220542980509E-003,
      -2.2480938654752945E-002,
      -2.0467235379121578E-002,
      -1.1771119223314612E-003,
      +6.0013717421124818E-003,
      +5.4638060238403757E-003
    },
    {
      +7.9014310551776979E-003,
      -2.9598497555107853E-002,
      +1.1651423916206259E-001,
      -2.0072266656029408E-003,
      +7.5190042334755038E-003,
      -2.9598497555107864E-002,
      +5.3583676268861434E-004,
      -2.0072266656029408E-003,
      +7.9014310551777014E-003
    },
    {
      -4.0758173076004510E-003,
      +6.4177683569844601E-002,
      +5.4638060238403749E-003,
      -4.4768234284154971E-003,
      +7.0491912591652800E-002,
      +6.0013717421124801E-003,
      +8.7808628730366711E-004,
      -1.3826317432949743E-002,
      -1.1771119223314610E-003
    },
    {
      -1.2002743484224964E-002,
      +6.1194627467012229E-002,
      +5.5713191674947982E-002,
      -1.3183653530112371E-002,
      +6.7215363511659812E-002,
      +6.1194627467012229E-002,
      +2.5858481055732723E-003,
      -1.3183653530112371E-002,
      -1.2002743484224964E-002
    },
    {
      +5.4638060238403749E-003,
      -2.0467235379121578E-002,
      +8.0569101641365728E-002,
      +6.0013717421124801E-003,
      -2.2480938654752938E-002,
      +8.8496027817990264E-002,
      -1.1771119223314610E-003,
      +4.4094220542980501E-003,
      -1.7357653199944128E-002
    },
    {
      -3.9971637753709799E-004,
      +6.2939207671118564E-003,
      +5.3583676268861434E-004,
      -4.6950552569634662E-003,
      +7.3928183695193897E-002,
      +6.2939207671118555E-003,
      +2.9817510405538797E-004,
      -4.6950552569634653E-003,
      -3.9971637753709782E-004
    },
    {
      -1.1771119223314612E-003,
      +6.0013717421124818E-003,
      +5.4638060238403757E-003,
      -1.3826317432949743E-002,
      +7.0491912591652814E-002,
      +6.4177683569844587E-002,
      +8.7808628730366700E-004,
      -4.4768234284154971E-003,
      -4.0758173076004501E-003
    },
    {
      +5.3583676268861434E-004,
      -2.0072266656029408E-003,
      +7.9014310551777014E-003,
      +6.2939207671118555E-003,
      -2.3576817558299032E-002,
      +9.2809945996529489E-002,
      -3.9971637753709782E-004,
      +1.4973242366666812E-003,
      -5.8942043895747606E-003
    }
  },
  {
    {
      -5.8942043895747588E-003,
      +1.4973242366666819E-003,
      -3.9971637753709793E-004,
      +9.2809945996529461E-002,
      -2.3576817558299039E-002,
      +6.2939207671118564E-003,
      +7.9014310551776979E-003,
      -2.0072266656029408E-003,
      +5.3583676268861423E-004
    },
    {
      -4.0758173076004501E-003,
      -4.4768234284154980E-003,
      +8.7808628730366700E-004,
      +6.4177683569844587E-002,
      +7.0491912591652814E-002,
      -1.3826317432949743E-002,
      +5.4638060238403749E-003,
      +6.0013717421124810E-003,
      -1.1771119223314612E-003
    },
    {
      -3.9971637753709793E-004,
      -4.6950552569634662E-003,
      +2.9817510405538792E-004,
      +6.2939207671118564E-003,
      +7.3928183695193897E-002,
      -4.6950552569634653E-003,
      +5.3583676268861423E-004,
      +6.2939207671118555E-003,
      -3.9971637753709782E-004
    },
    {
      -1.7357653199944124E-002,
      +4.4094220542980509E-003,
      -1.1771119223314612E-003,
      +8.8496027817990236E-002,
      -2.2480938654752942E-002,
      +6.0013717421124810E-003,
      +8.0569101641365715E-002,
      -2.0467235379121582E-002,
      +5.4638060238403757E-003
    },
    {
      -1.2002743484224964E-002,
      -1.3183653530112371E-002,
      +2.5858481055732723E-003,
      +6.1194627467012229E-002,
      +6.7215363511659812E-002,
      -1.3183653530112371E-002,
      +5.5713191674947982E-002,
      +6.1194627467012229E-002,
      -1.2002743484224964E-002
    },
    {
      -1.1771119223314612E-003,
      -1.3826317432949743E-002,
      +8.7808628730366700E-004,
      +6.0013717421124810E-003,
      +7.0491912591652800E-002,
      -4.4768234284154971E-003,
      +5.4638060238403757E-003,
      +6.4177683569844601E-002,
      -4.0758173076004501E-003
    },
    {
      +7.9014310551776979E-003,
      -2.0072266656029408E-003,
      +5.3583676268861423E-004,
      -2.9598497555107850E-002,
      +7.5190042334755047E-003,
      -2.0072266656029408E-003,
      +1.1651423916206258E-001,
      -2.9598497555107867E-002,
      +7.9014310551777014E-003
    },
    {
      +5.4638060238403749E-003,
      +6.0013717421124810E-003,
      -1.1771119223314612E-003,
      -2.0467235379121575E-002,
      -2.2480938654752942E-002,
      +4.4094220542980501E-003,
      +8.0569101641365728E-002,
      +8.8496027817990264E-002,
      -1.7357653199944128E-002
    },
    {
      +5.3583676268861423E-004,
      +6.2939207671118555E-003,
      -3.9971637753709782E-004,
      -2.0072266656029408E-003,
      -2.3576817558299032E-002,
      +1.4973242366666814E-003,
      +7.9014310551777014E-003,
      +9.2809945996529475E-002,
      -5.8942043895747597E-003
    }
  },
  {
    {
      +2.9817510405538797E-004,
      -4.6950552569634662E-003,
      -3.9971637753709788E-004,
      -4.6950552569634662E-003,
      +7.3928183695193897E-002,
      +6.2939207671118555E-003,
      -3.9971637753709788E-004,
      +6.2939207671118555E-003,
      +5.3583676268861423E-004
    },
    {
      +8.7808628730366700E-004,
      -4.4768234284154980E-003,
      -4.0758173076004501E-003,
      -1.3826317432949743E-002,
      +7.0491912591652814E-002,
      +6.4177683569844587E-002,
      -1.1771119223314612E-003,
      +6.0013717421124810E-003,
      +5.4638060238403749E-003
    },
    {
      -3.9971637753709788E-004,
      +1.4973242366666814E-003,
      -5.8942043895747606E-003,
      +6.2939207671118555E-003,
      -2.3576817558299032E-002,
      +9.2809945996529489E-002,
      +5.3583676268861423E-004,
      -2.0072266656029404E-003,
      +7.9014310551776996E-003
    },
    {
      +8.7808628730366711E-004,
      -1.3826317432949743E-002,
      -1.1771119223314610E-003,
      -4.4768234284154971E-003,
      +7.0491912591652800E-002,
      +6.0013717421124801E-003,
      -4.0758173076004510E-003,
      +6.4177683569844601E-002,
      +5.4638060238403749E-003
    },
    {
      +2.5858481055732723E-003,
      -1.3183653530112371E-002,
      -1.2002743484224964E-002,
      -1.3183653530112371E-002,
      +6.7215363511659812E-002,
      +6.1194627467012229E-002,
      -1.2002743484224964E-002,
      +6.1194627467012229E-002,
      +5.5713191674947982E-002
    },
    {
      -1.1771119223314610E-003,
      +4.4094220542980501E-003,
      -1.7357653199944128E-002,
      +6.0013717421124801E-003,
      -2.2480938654752938E-002,
      +8.8496027817990264E-002,
      +5.4638060238403749E-003,
      -2.0467235379121578E-002,
      +8.0569101641365728E-002
    },
    {
      -3.9971637753709788E-004,
      +6.2939207671118555E-003,
      +5.3583676268861423E-004,
      +1.4973242366666817E-003,
      -2.3576817558299032E-002,
      -2.0072266656029404E-003,
      -5.8942043895747606E-003,
      +9.2809945996529475E-002,
      +7.9014310551776996E-003
    },
    {
      -1.1771119223314612E-003,
      +6.0013717421124810E-003,
      +5.4638060238403749E-003,
      +4.4094220542980501E-003,
      -2.2480938654752942E-002,
      -2.0467235379121575E-002,
      -1.7357653199944128E-002,
      +8.8496027817990264E-002,
      +8.0569101641365728E-002
    },
    {
      +5.3583676268861423E-004,
      -2.0072266656029404E-003,
      +7.9014310551776996E-003,
      -2.0072266656029404E-003,
      +7.5190042334755021E-003,
      -2.9598497555107860E-002,
      +7.9014310551776996E-003,
      -2.9598497555107857E-002,
      +1.1651423916206261E-001
    }
  }
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
dginterp_interp_order3
  ( Box const& bx,
    Array4<Real> const& fine, const int fcomp, const int ncomp,
    Array4<Real const> const& crse, const int ccomp,
    IntVect const& ratio) noexcept
{

  amrex::ignore_unused( fcomp, ccomp );

  const int nDOFX = nDOFX_Order3_2D;

  // Get bounds on fine grid
  const auto loFine = amrex::lbound(bx);
  const auto hiFine = amrex::ubound(bx);

  // Number of fields in MultiFab
  int nFields = ncomp / nDOFX;

  // Index for projection matrix
  int iProj;

  // Loop over fields
  for( int iField = 0; iField < nFields; iField++ )
  {

    // Loop over fine cells in y-direction
    for( int jFine = loFine.y; jFine <= hiFine.y; jFine++ )
    {

      // Get coarse element corresponding to fine element jFine
      const int jCrse = amrex::coarsen(jFine,ratio[1]);
      AMREX_PRAGMA_SIMD

      // Loop over fine cells in x-direction
      for( int iFine = loFine.x; iFine <= hiFine.x; iFine++ )
      {

        if     ( ( iFine % 2 == 0 ) & ( jFine % 2 == 0 ) ) { iProj = 0; }
        else if( ( jFine % 2 == 0 )                      ) { iProj = 1; }
        else if( ( iFine % 2 == 0 )                      ) { iProj = 2; }
        else                                               { iProj = 3; }

        // Get coarse element corresponding to fine element iFine
        const int iCrse = amrex::coarsen(iFine,ratio[0]);

        // Loop over DOF of fine element
        for( int iNX = 0; iNX < nDOFX; iNX++ )
        {

          // Initialize fine data to zero
          fine(iFine,jFine,0,nDOFX*iField+iNX) = Zero_2D;

          // Project coarse data onto fine data
          for( int jNX = 0; jNX < nDOFX; jNX++ )
          {

            fine(iFine,jFine,0,nDOFX*iField+iNX)
              = fine(iFine,jFine,0,nDOFX*iField+iNX)
                  + ProjectionMatrix_Order3_2D[iProj][iNX][jNX]
                      * crse(iCrse,jCrse,0,nDOFX*iField+jNX);

          } // jNX

          fine(iFine,jFine,0,nDOFX*iField+iNX)
            = fine(iFine,jFine,0,nDOFX*iField+iNX)
                / GaussLegendreWeights_Order3_2D[iNX];

        } // iNX
      } // iFine
    } // jFine
  } // iField

} // End of function dginterp_interp_order3

namespace {
    static constexpr int ix   = 0;
    static constexpr int iy   = 1;
    static constexpr int ixy  = 2;
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
nodebilin_slopes (Box const& bx, Array4<T> const& slope, Array4<T const> const& u,
                  const int icomp, const int ncomp, IntVect const& ratio) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    const Real rx = Real(1.)/ratio[0];
    const Real ry = Real(1.)/ratio[1];

    for (int n = 0; n < ncomp; ++n) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                T dx0 = u(i+1,j,0,n+icomp) - u(i,j,0,n+icomp);
                T d0x = u(i,j+1,0,n+icomp) - u(i,j,0,n+icomp);
                T dx1 = u(i+1,j+1,0,n+icomp) - u(i,j+1,0,n+icomp);

                slope(i,j,0,n+ncomp*ix ) = rx*dx0;
                slope(i,j,0,n+ncomp*iy ) = ry*d0x;
                slope(i,j,0,n+ncomp*ixy) = rx*ry*(dx1 - dx0);
            }
        }
    }
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
nodebilin_interp (Box const& bx, Array4<T> const& fine, const int fcomp, const int ncomp,
                  Array4<T const> const& slope, Array4<T const> const& crse,
                  const int ccomp, IntVect const& ratio) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);
    const auto chi = amrex::ubound(slope);

    for (int n = 0; n < ncomp; ++n) {
        for (int j = lo.y; j <= hi.y; ++j) {
            const int jc = amrex::min(amrex::coarsen(j,ratio[1]),chi.y);
            const Real fy = j - jc*ratio[1];
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                const int ic = amrex::min(amrex::coarsen(i,ratio[0]),chi.x);
                const Real fx = i - ic*ratio[0];
                fine(i,j,0,n+fcomp) = crse(ic,jc,0,n+ccomp)
                    + fx*slope(ic,jc,0,n+ncomp*ix)
                    + fy*slope(ic,jc,0,n+ncomp*iy)
                    + fx*fy*slope(ic,jc,0,n+ncomp*ixy);
            }
        }
    }
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
facediv_face_interp (int ci, int cj, int /*ck*/,
                     int nc, int nf, int idir,
                     Array4<T const> const& crse,
                     Array4<T> const& fine,
                     Array4<int const> const& mask,
                     IntVect const& ratio) noexcept
{
    if (mask) {
        if (!mask(ci, cj, 0, nc))
            { return; }
    }

    const int fi = ci*ratio[0];
    const int fj = cj*ratio[1];

    switch (idir) {
        case 0:
        {
            const Real neg = crse(ci, cj-1, 0, nc);
            const Real cen = crse(ci, cj  , 0, nc);
            const Real pos = crse(ci, cj+1, 0, nc);

            fine(fi, fj  , 0, nf) = Real(0.125)*(8*cen + neg - pos);
            fine(fi, fj+1, 0, nf) = Real(0.125)*(8*cen + pos - neg);

            break;
        }
        case 1:
        {
            const Real neg = crse(ci-1, cj, 0, nc);
            const Real cen = crse(ci  , cj, 0, nc);
            const Real pos = crse(ci+1, cj, 0, nc);

            fine(fi  , fj, 0, nf) = Real(0.125)*(8*cen + neg - pos);
            fine(fi+1, fj, 0, nf) = Real(0.125)*(8*cen + pos - neg);

            break;
        }
    }
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
facediv_int (int ci, int cj, int /*ck*/, int nf,
             GpuArray<Array4<T>, AMREX_SPACEDIM> const& fine,
             IntVect const& ratio,
             GpuArray<Real, AMREX_SPACEDIM> const& cellSize) noexcept
{
    const int fi = ci*ratio[0];
    const int fj = cj*ratio[1];

    // References to fine exterior values.
    const Real umm = fine[0](fi,   fj,   0, nf);
    const Real ump = fine[0](fi,   fj+1, 0, nf);
    const Real upm = fine[0](fi+2, fj,   0, nf);
    const Real upp = fine[0](fi+2, fj+1, 0, nf);

    const Real vmm = fine[1](fi,   fj,   0, nf);
    const Real vmp = fine[1](fi+1, fj,   0, nf);
    const Real vpm = fine[1](fi,   fj+2, 0, nf);
    const Real vpp = fine[1](fi+1, fj+2, 0, nf);

    const Real dxdy = cellSize[0]/cellSize[1];
    const Real x_corr = Real(0.25)*dxdy * (vpp+vmm-vmp-vpm);
    const Real y_corr = Real(0.25)/dxdy * (upp+umm-ump-upm);

    // Calc fine faces on interior of coarse cells.
    fine[0](fi+1,fj  ,0,nf) = Real(0.5)*(umm+upm) + x_corr;
    fine[0](fi+1,fj+1,0,nf) = Real(0.5)*(ump+upp) + x_corr;
    fine[1](fi,  fj+1,0,nf) = Real(0.5)*(vmm+vpm) + y_corr;
    fine[1](fi+1,fj+1,0,nf) = Real(0.5)*(vmp+vpp) + y_corr;
}



template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
face_linear_interp_x (int i, int j, int /*k*/, int n, Array4<T> const& fine,
                      Array4<T const> const& crse, IntVect const& ratio) noexcept
{
    const int ii = amrex::coarsen(i,ratio[0]);
    const int jj = amrex::coarsen(j,ratio[1]);
    if (i-ii*ratio[0] == 0) {
        fine(i,j,0,n) = crse(ii,jj,0,n);
    } else {
        Real const w = static_cast<Real>(i-ii*ratio[0]) * (Real(1.)/ratio[0]);
        fine(i,j,0,n) = (Real(1.)-w) * crse(ii,jj,0,n) + w * crse(ii+1,jj,0,n);
    }
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
face_linear_interp_y (int i, int j, int /*k*/, int n, Array4<T> const& fine,
                      Array4<T const> const& crse, IntVect const& ratio) noexcept
{
    const int ii = amrex::coarsen(i,ratio[0]);
    const int jj = amrex::coarsen(j,ratio[1]);
    if (j-jj*ratio[1] == 0) {
        fine(i,j,0,n) = crse(ii,jj,0,n);
    } else {
        Real const w = static_cast<Real>(j-jj*ratio[1]) * (Real(1.)/ratio[1]);
        fine(i,j,0,n) = (Real(1.)-w) * crse(ii,jj,0,n) + w * crse(ii,jj+1,0,n);
    }
}

}

#endif
