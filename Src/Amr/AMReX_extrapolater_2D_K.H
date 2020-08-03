#ifndef AMReX_extrapolater_2D_K_H_
#define AMReX_extrapolater_2D_K_H_

namespace amrex {

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
amrex_first_order_extrap(amrex::Box const& bx,
                         int               nComp,
                         amrex::Array4<const int>   const& mask,
                         amrex::Array4<amrex::Real> const& data) noexcept
{
   using namespace amrex::literals;

   const int finecell = 1;
   const int crsecell = 0;

   const auto lo = amrex::lbound(bx);
   const auto hi = amrex::ubound(bx);
   int idlo[3] = {AMREX_D_DECL(lo.x,lo.y,lo.z)};            // base box bounds
   int idhi[3] = {AMREX_D_DECL(hi.x,hi.y,hi.z)};
   int idlo_e[3] = {AMREX_D_DECL(lo.x,lo.y,lo.z)};          // extended box bounds
   int idhi_e[3] = {AMREX_D_DECL(hi.x,hi.y,hi.z)};
   for (int dim = 0; dim < AMREX_SPACEDIM; dim++) {   // Extend only AMREX_SPACEDIMs
      idlo_e[dim] -= 1;
      idhi_e[dim] += 1;
   }

   int k = idlo[2];

   for (int n = 0; n < nComp; n++) {

      // set all crse cells to zero first
      for (int j = idlo_e[1]; j <= idhi_e[1]; ++j) {
         for (int i = idlo_e[0]; i <= idhi_e[0]; ++i) {
            if (mask(i,j,k) == crsecell) data(i,j,k,n) = 0.0_rt;
         }
      }

      // Corners
      // xlo, ylo
      {
         int i = idlo_e[0];
         int j = idlo_e[1];
         if ( mask(i,j,k) == crsecell ) {
            if ( ( mask(i+1,j,k) == finecell ) ||
                 ( mask(i,j+1,k) == finecell ) ) {
               data(i,j,k,n) = ( mask(i+1,j,k) * data(i+1,j,k,n) + mask(i,j+1,k) * data(i,j+1,k,n) )
                               / ( mask(i+1,j,k) + mask(i,j+1,k) );
            } else {
               data(i,j,k,n) = data(i+1,j+1,k,n);
            }
         }
      }
      // xlo, yhi
      {
         int i = idlo_e[0];
         int j = idhi_e[1];
         if ( mask(i,j,k) == crsecell ) {
            if ( ( mask(i+1,j,k) == finecell ) ||
                 ( mask(i,j-1,k) == finecell ) ) {
               data(i,j,k,n) = ( mask(i+1,j,k) * data(i+1,j,k,n) + mask(i,j-1,k) * data(i,j-1,k,n) )
                               / ( mask(i+1,j,k) + mask(i,j-1,k) );
            } else {
               data(i,j,k,n) = data(i+1,j-1,k,n);
            }
         }
      }
      // xhi, ylo
      {
         int i = idhi_e[0];
         int j = idlo_e[1];
         if ( mask(i,j,k) == crsecell ) {
            if ( ( mask(i-1,j,k) == finecell ) ||
                 ( mask(i,j+1,k) == finecell ) ) {
               data(i,j,k,n) = ( mask(i-1,j,k) * data(i-1,j,k,n) + mask(i,j+1,k) * data(i,j+1,k,n) )
                               / ( mask(i-1,j,k) + mask(i,j+1,k) );
            } else {
               data(i,j,k,n) = data(i-1,j+1,k,n);
            }
         }
      }
      // xhi, yhi
      {
         int i = idhi_e[0];
         int j = idhi_e[1];
         if ( mask(i,j,k) == crsecell ) {
            if ( ( mask(i-1,j,k) == finecell ) ||
                 ( mask(i,j-1,k) == finecell ) ) {
               data(i,j,k,n) = ( mask(i-1,j,k) * data(i-1,j,k,n) + mask(i,j-1,k) * data(i,j-1,k,n) )
                               / ( mask(i-1,j,k) + mask(i,j-1,k) );
            } else {
               data(i,j,k,n) = data(i-1,j-1,k,n);
            }
         }
      }
      // Edges
      // xlo, y-valid
      {
         int i = idlo_e[0];
         for (int j = idlo[1]; j <= idhi[1]; ++j) {
            if ( mask(i,j,k) == crsecell ) {
               data(i,j,k,n) = ( mask(i,j-1,k) * data(i,j-1,k,n) + mask(i+1,j,k) + mask(i,j+1,k) * data(i,j+1,k,n) )
                               / ( mask(i,j-1,k) + 1 + mask(i,j+1,k) );
            }
         }
      }
      // xhi, y-valid
      {
         int i = idhi_e[0];
         for (int j = idlo[1]; j <= idhi[1]; ++j) {
            if ( mask(i,j,k) == crsecell ) {
               data(i,j,k,n) = ( mask(i,j-1,k) * data(i,j-1,k,n) + mask(i-1,j,k) + mask(i,j+1,k) * data(i,j+1,k,n) )
                               / ( mask(i,j-1,k) + 1 + mask(i,j+1,k) );
            }
         }
      }
      // x-valid, ylo
      {
         int j = idlo_e[1];
         for (int i = idlo[0]; i <= idhi[0]; ++i) {
            if ( mask(i,j,k) == crsecell ) {
               data(i,j,k,n) = ( mask(i-1,j,k) * data(i-1,j,k,n) + mask(i+1,j,k) * data(i+1,j,k,n) + data(i,j+1,k,n) )
                               / ( mask(i-1,j,k) + mask(i+1,j,k) + 1 );
            }
         }
      }
      // x-valid, yhi
      {
         int j = idhi_e[1];
         for (int i = idlo[0]; i <= idhi[0]; ++i) {
            if ( mask(i,j,k) == crsecell ) {
               data(i,j,k,n) = ( mask(i-1,j,k) * data(i-1,j,k,n) + mask(i+1,j,k) * data(i+1,j,k,n) + data(i,j-1,k,n) )
                               / ( mask(i-1,j,k) + mask(i+1,j,k) + 1 );
            }
         }
      }
   }
}

}
#endif
