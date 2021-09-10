#ifndef AMREX_INTERP_1D_C_H_
#define AMREX_INTERP_1D_C_H_
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
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            const int ic = amrex::coarsen(i,ratio[0]);
            fine(i,0,0,n+fcomp) = crse(ic,0,0,n+ccomp);
        }
    }
}

const Real Zero_1D = 0.0;

/* First order DG methods */
const int nDOFX_Order1_1D = 1;
const Real GaussLegendreWeights_Order1_1D[nDOFX_Order1_1D] =
{
  +1.0000000000000000E+000
};
const Real ProjectionMatrix_Order1_1D[2][nDOFX_Order1_1D][nDOFX_Order1_1D] =
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
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
dginterp_interp_order1
  ( Box const& bx,
    Array4<Real> const& fine, const int fcomp, const int ncomp,
    Array4<Real const> const& crse, const int ccomp,
    IntVect const& ratio) noexcept
{

  amrex::ignore_unused( fcomp, ccomp );

  const int nDOFX = nDOFX_Order1_1D;

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

    // Loop over fine cells in x-direction
    for( int iFine = loFine.x; iFine <= hiFine.x; iFine++ )
    {

      if( iFine % 2 == 0 ) { iProj = 0; }
      else                 { iProj = 1; }

      // Get coarse element corresponding to fine element iFine
      const int iCrse = amrex::coarsen(iFine,ratio[0]);

      // Loop over DOF of fine element
      for( int iNX = 0; iNX < nDOFX; iNX++ )
      {

        // Initialize fine data to zero
        fine(iFine,0,0,nDOFX*iField+iNX) = Zero_1D;

        // Project coarse data onto fine data
        for( int jNX = 0; jNX < nDOFX; jNX++ )
        {

          fine(iFine,0,0,nDOFX*iField+iNX)
            = fine(iFine,0,0,nDOFX*iField+iNX)
                + ProjectionMatrix_Order1_1D[iProj][iNX][jNX]
                    * crse(iCrse,0,0,nDOFX*iField+jNX);

        } // jNX

        fine(iFine,0,0,nDOFX*iField+iNX)
          = fine(iFine,0,0,nDOFX*iField+iNX)
              / GaussLegendreWeights_Order1_1D[iNX];

      } // iNX
    } // iFine
  } // iField

} // End of function dginterp_interp_order1

/* Second order DG methods */
const int nDOFX_Order2_1D = 2;
const Real GaussLegendreWeights_Order2_1D[nDOFX_Order2_1D] =
{
  +5.0000000000000000E-001,
  +5.0000000000000000E-001
};
const Real ProjectionMatrix_Order2_1D[2][nDOFX_Order2_1D][nDOFX_Order2_1D] =
{
  {
    {
      +5.9150635094610970E-001,
      -9.1506350946109677E-002
    },
    {
      +3.4150635094610965E-001,
      +1.5849364905389032E-001
    }
  },
  {
    {
      +1.5849364905389032E-001,
      +3.4150635094610965E-001
    },
    {
      -9.1506350946109677E-002,
      +5.9150635094610970E-001
    }
  },
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
dginterp_interp_order2
  ( Box const& bx,
    Array4<Real> const& fine, const int fcomp, const int ncomp,
    Array4<Real const> const& crse, const int ccomp,
    IntVect const& ratio) noexcept
{

  amrex::ignore_unused( fcomp, ccomp );

  const int nDOFX = nDOFX_Order2_1D;

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

    // Loop over fine cells in x-direction
    for( int iFine = loFine.x; iFine <= hiFine.x; iFine++ )
    {

      if( iFine % 2 == 0 ) { iProj = 0; }
      else                 { iProj = 1; }

      // Get coarse element corresponding to fine element iFine
      const int iCrse = amrex::coarsen(iFine,ratio[0]);

      // Loop over DOF of fine element
      for( int iNX = 0; iNX < nDOFX; iNX++ )
      {

        // Initialize fine data to zero
        fine(iFine,0,0,nDOFX*iField+iNX) = Zero_1D;

        // Project coarse data onto fine data
        for( int jNX = 0; jNX < nDOFX; jNX++ )
        {

          fine(iFine,0,0,nDOFX*iField+iNX)
            = fine(iFine,0,0,nDOFX*iField+iNX)
                + ProjectionMatrix_Order2_1D[iProj][iNX][jNX]
                    * crse(iCrse,0,0,nDOFX*iField+jNX);

        } // jNX

        fine(iFine,0,0,nDOFX*iField+iNX)
          = fine(iFine,0,0,nDOFX*iField+iNX)
              / GaussLegendreWeights_Order2_1D[iNX];

      } // iNX
    } // iFine
  } // iField

} // End of function dginterp_interp_order2

/* Third order DG methods */
const int nDOFX_Order3_1D = 3;
const Real GaussLegendreWeights_Order3_1D[nDOFX_Order3_1D] =
{
  +2.7777777777777779E-001,
  +4.4444444444444442E-001,
  +2.7777777777777779E-001
};
const Real ProjectionMatrix_Order3_1D[2][nDOFX_Order3_1D][nDOFX_Order3_1D] =
{
  {
    {
      +3.4134182158367665E-001,
      -8.6712191954047071E-002,
      +2.3148148148148143E-002
    },
    {
      +2.3603642022990431E-001,
      +2.5925925925925924E-001,
      -5.0851235044719140E-002
    },
    {
      +2.3148148148148143E-002,
      +2.7189737713923223E-001,
      -1.7267747509602631E-002
    }
  },
  {
    {
      -1.7267747509602634E-002,
      +2.7189737713923223E-001,
      +2.3148148148148140E-002
    },
    {
      -5.0851235044719140E-002,
      +2.5925925925925924E-001,
      +2.3603642022990431E-001
    },
    {
      +2.3148148148148140E-002,
      -8.6712191954047058E-002,
      +3.4134182158367671E-001
    }
  },
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
dginterp_interp_order3
  ( Box const& bx,
    Array4<Real> const& fine, const int fcomp, const int ncomp,
    Array4<Real const> const& crse, const int ccomp,
    IntVect const& ratio) noexcept
{

  amrex::ignore_unused( fcomp, ccomp );

  const int nDOFX = nDOFX_Order3_1D;

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

    // Loop over fine cells in x-direction
    for( int iFine = loFine.x; iFine <= hiFine.x; iFine++ )
    {

      if( iFine % 2 == 0 ) { iProj = 0; }
      else                 { iProj = 1; }

      // Get coarse element corresponding to fine element iFine
      const int iCrse = amrex::coarsen(iFine,ratio[0]);

      // Loop over DOF of fine element
      for( int iNX = 0; iNX < nDOFX; iNX++ )
      {

        // Initialize fine data to zero
        fine(iFine,0,0,nDOFX*iField+iNX) = Zero_1D;

        // Project coarse data onto fine data
        for( int jNX = 0; jNX < nDOFX; jNX++ )
        {

          fine(iFine,0,0,nDOFX*iField+iNX)
            = fine(iFine,0,0,nDOFX*iField+iNX)
                + ProjectionMatrix_Order3_1D[iProj][iNX][jNX]
                    * crse(iCrse,0,0,nDOFX*iField+jNX);

        } // jNX

        fine(iFine,0,0,nDOFX*iField+iNX)
          = fine(iFine,0,0,nDOFX*iField+iNX)
              / GaussLegendreWeights_Order3_1D[iNX];

      } // iNX
    } // iFine
  } // iField

} // End of function dginterp_interp_order3

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
nodebilin_slopes (Box const& bx, Array4<T> const& slope, Array4<T const> const& u,
                  const int icomp, const int ncomp, IntVect const& ratio) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    const Real rx = Real(1.)/ratio[0];

    for (int n = 0; n < ncomp; ++n) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            T dx0 = u(i+1,0,0,n+icomp) - u(i,0,0,n+icomp);
            slope(i,0,0,n) = rx*dx0;
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
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            const int ic = amrex::min(amrex::coarsen(i,ratio[0]),chi.x);
            const Real fx = i - ic*ratio[0];
            fine(i,0,0,n+fcomp) = crse(ic,0,0,n+ccomp) + fx*slope(ic,0,0,0);
        }
    }
}

// Or, compile time error this?
// Remove these functions in a way that gives a good error message?

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
facediv_face_interp (int /*ci*/, int /*cj*/, int /*ck*/,
                     int /*nc*/, int /*nf*/, int /*idir*/,
                     Array4<T const> const& /*crse*/,
                     Array4<T> const& /*fine*/,
                     Array4<const int> const& /*mask*/,
                     IntVect const& /*ratio*/) noexcept
{
    amrex::Abort("No 1D version of FaceDiv exists.\n");
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
facediv_int (int /*ci*/, int /*cj*/, int /*ck*/, int /*nf*/,
             GpuArray<Array4<T>, AMREX_SPACEDIM> const& /*fine*/,
             IntVect const& /*ratio*/,
             GpuArray<Real, AMREX_SPACEDIM> const& /*cellSize*/) noexcept
{
    amrex::Abort("No 1D version of FaceDiv exists.\n");
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
face_linear_interp_x (int i, int /*j*/, int /*k*/, int n, Array4<T> const& fine,
                      Array4<T const> const& crse, IntVect const& ratio) noexcept
{
    const int ii = amrex::coarsen(i,ratio[0]);
    if (i-ii*ratio[0] == 0) {
        fine(i,0,0,n) = crse(ii,0,0,n);
    } else {
        Real const w = static_cast<Real>(i-ii*ratio[0]) * (Real(1.)/ratio[0]);
        fine(i,0,0,n) = (Real(1.)-w) * crse(ii,0,0,n) + w * crse(ii+1,0,0,n);
    }
}

}

#endif
