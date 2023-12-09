#ifndef AMREX_INTERP_1D_C_H_
#define AMREX_INTERP_1D_C_H_
#include <AMReX_Config.H>

#include <AMReX_FArrayBox.H>
#include <AMReX_BCRec.H>
#include <AMReX_Vector.H>
#include <AMReX_Array.H>
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

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
dginterpConservative_interp
  ( Box                const & FineBox,
    Array4<Real>       const & FineArr,
    Array4<Real const> const & FineArrG,
    int                const   nComp,
    Array4<Real const> const & CrseArr,
    Array4<Real const> const & CrseArrG,
    IntVect            const & RefRatio,
    int                        nDOFX,
    Array4<Real const>         CoarseToFineProjectionMatrix ) noexcept
{
  const Real Zero = 0.0;

  // Get bounds on fine grid
  const auto loFine = amrex::lbound( FineBox );
  const auto hiFine = amrex::ubound( FineBox );

  // Number of fields in MultiFab
  int nFields = nComp / nDOFX;

  // Index for projection matrix
  int iProj;

  // Loop over fields
  for( int iField = 0; iField < nFields; iField++ )
  {
    // Loop over fine cells in x-direction
    for( int iFine = loFine.x; iFine <= hiFine.x; iFine++ )
    {
      // Get coarse element corresponding to fine element iFine
      const int iCrse = amrex::coarsen( iFine, RefRatio[0] );

      if( iFine % 2 == 0 ) { iProj = 0; }
      else                 { iProj = 1; }

      // Loop over DOF of fine element
      for( int iNX = 0; iNX < nDOFX; iNX++ )
      {
        // Initialize fine data to zero
        FineArr(iFine,0,0,nDOFX*iField+iNX) = Zero;

        // Project coarse data onto fine data
        for( int jNX = 0; jNX < nDOFX; jNX++ )
        {
          FineArr(iFine,0,0,nDOFX*iField+iNX)
            += CoarseToFineProjectionMatrix(0,iProj,iNX,jNX)
                 * CrseArr (iCrse,0,0,nDOFX*iField+jNX)
                 * CrseArrG(iCrse,0,0,jNX);
        } // jNX

        FineArr(iFine,0,0,nDOFX*iField+iNX)
          /= FineArrG(iFine,0,0,iNX);

      } // iNX
    } // iFine
  } // iField
} // end void dginterpConservative_interp

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
dginterpPointWise_interp
  ( Box                const & FineBox,
    Array4<Real>       const & FineArr,
    int                const   nComp,
    Array4<Real const> const & CrseArr,
    IntVect            const & RefRatio,
    int                        nDOFX,
    Array4<Real const>         CoarseToFineProjectionMatrix ) noexcept
{
  const Real Zero = 0.0;

  // Get bounds on fine grid
  const auto loFine = amrex::lbound( FineBox );
  const auto hiFine = amrex::ubound( FineBox );

  // Number of fields in MultiFab
  int nFields = nComp / nDOFX;

  // Index for projection matrix
  int iProj;

  // Loop over fields
  for( int iField = 0; iField < nFields; iField++ )
  {
    // Loop over fine cells in x-direction
    for( int iFine = loFine.x; iFine <= hiFine.x; iFine++ )
    {
      // Get coarse element corresponding to fine element iFine
      const int iCrse = amrex::coarsen( iFine, RefRatio[0] );

      if( iFine % 2 == 0 ) { iProj = 0; }
      else                 { iProj = 1; }

      // Loop over DOF of fine element
      for( int iNX = 0; iNX < nDOFX; iNX++ )
      {
        // Initialize fine data to zero
        FineArr(iFine,0,0,nDOFX*iField+iNX) = Zero;

        // Project coarse data onto fine data
        for( int jNX = 0; jNX < nDOFX; jNX++ )
        {
          FineArr(iFine,0,0,nDOFX*iField+iNX)
            += CoarseToFineProjectionMatrix(0,iProj,iNX,jNX)
                 * CrseArr(iCrse,0,0,nDOFX*iField+jNX);
        } // jNX
      } // iNX
    } // iFine
  } // iField
} // end void dginterpPointWise_interp

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
cginterp_interp
  ( Box const& bx,
    Array4<Real> const& fine, const int fcomp, const int ncomp,
    Array4<Real const> const& crse, const int ccomp,
    IntVect const& ratio,
    int nDOFX, Array4<Real> ProjectionMatrix,
    Array4<Real> WeightsX_q ) noexcept
{
  amrex::ignore_unused
  ( bx,
    fine, fcomp, ncomp,
    crse, ccomp,
    ratio,
    nDOFX, ProjectionMatrix,
    WeightsX_q );
  amrex::Abort( "Src/AmrCore/AMReX_Interp_1D_C.H: cginterp_interp not implemented. Calling amrex::Abort()" );
} // end void cginterp_interp

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
nodebilin_slopes (Box const& bx, Array4<T> const& slope, Array4<T const> const& u,
                  const int icomp, const int ncomp, IntVect const& ratio) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    const Real rx = Real(1.)/Real(ratio[0]);

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
        Real const w = static_cast<Real>(i-ii*ratio[0]) * (Real(1.)/Real(ratio[0]));
        fine(i,0,0,n) = (Real(1.)-w) * crse(ii,0,0,n) + w * crse(ii+1,0,0,n);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void ccquartic_interp (int i, int /*j*/, int /*k*/, int n,
                       Array4<Real const> const& crse,
                       Array4<Real>       const& fine ) noexcept

{
    // Note: there are asserts in CellConservativeQuartic::interp()
    //       to check whether ratio is all equal to 2.

    constexpr Array1D<Real, -2, 2> cL = { -0.01171875_rt,  0.0859375_rt, 0.5_rt, -0.0859375_rt, 0.01171875_rt };

    int ic = amrex::coarsen(i,2);
    int irx = i - 2*ic; // = abs(i % 2)

    Real ftmp = 2.0_rt * ( cL(-2)*crse(ic-2,0,0,n)
                         + cL(-1)*crse(ic-1,0,0,n)
                         + cL( 0)*crse(ic  ,0,0,n)
                         + cL( 1)*crse(ic+1,0,0,n)
                         + cL( 2)*crse(ic+2,0,0,n) );
    if (irx) {
        ftmp = 2.0_rt * crse(ic,0,0,n) - ftmp;
    }

    fine(i,0,0,n) = ftmp;

}

} // namespace amrex

#endif
