#ifndef AMReX_extrapolater_1D_K_H_
#define AMReX_extrapolater_1D_K_H_

namespace amrex {

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
amrex_first_order_extrap_cpu(amrex::Box const& bx,
                             int               nComp,
                             amrex::Array4<const int>   const& mask,
                             amrex::Array4<amrex::Real> const& data) noexcept
{
   constexpr int crsecell = 0;

   const auto lo = amrex::lbound(bx);
   const auto hi = amrex::ubound(bx);

   if (mask(lo.x-1,lo.y,lo.z) == crsecell) {
      for (int n = 0; n < nComp; n++) {
         data(lo.x-1,lo.y,lo.z,n) = data(lo.x,lo.y,lo.z,n);
      }
   }
   if (mask(hi.x+1,hi.y,hi.z) == crsecell) {
      for (int n = 0; n < nComp; n++) {
         data(hi.x+1,hi.y,hi.z,n) = data(hi.x,hi.y,hi.z,n);
      }
   }

}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
amrex_first_order_extrap_gpu(int i, int j, int k, int n,
                             amrex::Box const& bx,
                             amrex::Array4<const int>   const& mask,
                             amrex::Array4<amrex::Real> const& data) noexcept
{
   constexpr int crsecell = 0;

   const auto lo = amrex::lbound(bx);
   const auto hi = amrex::ubound(bx);

   if ( (i == lo.x-1) && ( mask(i,j,k) == crsecell ) ) {
      data(i,j,k,n) = data(i+1,j,k,n);
   }

   if ( (i == hi.x+1) && ( mask(i,j,k) == crsecell ) ) {
      data(i,j,k,n) = data(i-1,j,k,n);
   }
}

}
#endif
