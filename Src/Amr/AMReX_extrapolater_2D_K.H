#ifndef AMReX_extrapolater_2D_K_H_
#define AMReX_extrapolater_2D_K_H_

namespace amrex {

AMREX_GPU_HOST
AMREX_FORCE_INLINE
void
amrex_first_order_extrap_cpu(amrex::Box const& bx,
                             int               nComp,
                             amrex::Array4<const int>   const& mask,
                             amrex::Array4<amrex::Real> const& data) noexcept
{
   using namespace amrex::literals;

   constexpr int finecell = 1;
   constexpr int crsecell = 0;

   const auto lo = amrex::lbound(bx);
   const auto hi = amrex::ubound(bx);

   int k = lo.z;

   for (int n = 0; n < nComp; n++) {

      // set all crse cells to zero first
      for (int j = lo.y-1; j <= hi.y+1; ++j) {
         for (int i = lo.x-1; i <= hi.x+1; ++i) {
             if (mask(i,j,k) == crsecell) data(i,j,k,n) = Real(0.0);
         }
      }

      // Corners
      // xlo, ylo
      {
         int i = lo.x-1;
         int j = lo.y-1;
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
         int i = lo.x-1;
         int j = hi.y+1;
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
         int i = hi.x+1;
         int j = lo.y-1;
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
         int i = hi.x+1;
         int j = hi.y+1;
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
         int i = lo.x-1;
         for (int j = lo.y; j <= hi.y; ++j) {
            if ( mask(i,j,k) == crsecell ) {
               data(i,j,k,n) = ( mask(i,j-1,k) * data(i,j-1,k,n) + data(i+1,j,k,n) + mask(i,j+1,k) * data(i,j+1,k,n) )
                               / ( mask(i,j-1,k) + 1 + mask(i,j+1,k) );
            }
         }
      }
      // xhi, y-valid
      {
         int i = hi.x+1;
         for (int j = lo.y; j <= hi.y; ++j) {
            if ( mask(i,j,k) == crsecell ) {
               data(i,j,k,n) = ( mask(i,j-1,k) * data(i,j-1,k,n) + data(i-1,j,k,n) + mask(i,j+1,k) * data(i,j+1,k,n) )
                               / ( mask(i,j-1,k) + 1 + mask(i,j+1,k) );
            }
         }
      }
      // x-valid, ylo
      {
         int j = lo.y-1;
         for (int i = lo.x; i <= hi.x; ++i) {
            if ( mask(i,j,k) == crsecell ) {
               data(i,j,k,n) = ( mask(i-1,j,k) * data(i-1,j,k,n) + mask(i+1,j,k) * data(i+1,j,k,n) + data(i,j+1,k,n) )
                               / ( mask(i-1,j,k) + mask(i+1,j,k) + 1 );
            }
         }
      }
      // x-valid, yhi
      {
         int j = hi.y+1;
         for (int i = lo.x; i <= hi.x; ++i) {
            if ( mask(i,j,k) == crsecell ) {
               data(i,j,k,n) = ( mask(i-1,j,k) * data(i-1,j,k,n) + mask(i+1,j,k) * data(i+1,j,k,n) + data(i,j-1,k,n) )
                               / ( mask(i-1,j,k) + mask(i+1,j,k) + 1 );
            }
         }
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
   using namespace amrex::literals;

   constexpr int finecell = 1;
   constexpr int crsecell = 0;

   const auto lo = amrex::lbound(bx);
   const auto hi = amrex::ubound(bx);

   if ( mask(i,j,k) == crsecell ) {
      // Corners
      // xlo, ylo
      if ( (i == lo.x-1) && (j == lo.y-1) ) {
         if ( ( mask(i+1,j,k) == finecell ) ||
              ( mask(i,j+1,k) == finecell ) ) {
            data(i,j,k,n) = ( mask(i+1,j,k) * data(i+1,j,k,n) + mask(i,j+1,k) * data(i,j+1,k,n) )
                          / ( mask(i+1,j,k) + mask(i,j+1,k) );
         } else {
            data(i,j,k,n) = data(i+1,j+1,k,n);
         }
      // xlo, yhi
      } else if ( (i == lo.x-1) && (j == hi.y+1) ) {
         if ( ( mask(i+1,j,k) == finecell ) ||
              ( mask(i,j-1,k) == finecell ) ) {
            data(i,j,k,n) = ( mask(i+1,j,k) * data(i+1,j,k,n) + mask(i,j-1,k) * data(i,j-1,k,n) )
                          / ( mask(i+1,j,k) + mask(i,j-1,k) );
         } else {
            data(i,j,k,n) = data(i+1,j-1,k,n);
         }
      // xhi, ylo
      } else if ( (i == hi.x+1) && (j == lo.y-1) ) {
         if ( ( mask(i-1,j,k) == finecell ) ||
              ( mask(i,j+1,k) == finecell ) ) {
            data(i,j,k,n) = ( mask(i-1,j,k) * data(i-1,j,k,n) + mask(i,j+1,k) * data(i,j+1,k,n) )
                          / ( mask(i-1,j,k) + mask(i,j+1,k) );
         } else {
            data(i,j,k,n) = data(i-1,j+1,k,n);
         }
      // xhi, yhi
      } else if ( (i == hi.x+1) && (j == hi.y+1) ) {
         if ( ( mask(i-1,j,k) == finecell ) ||
              ( mask(i,j-1,k) == finecell ) ) {
            data(i,j,k,n) = ( mask(i-1,j,k) * data(i-1,j,k,n) + mask(i,j-1,k) * data(i,j-1,k,n) )
                          / ( mask(i-1,j,k) + mask(i,j-1,k) );
         } else {
            data(i,j,k,n) = data(i-1,j-1,k,n);
         }
      // Edges
      // xlo, y-valid
      } else if ( (i == lo.x-1) && (j >= lo.y) && (j <= hi.y) ) {
         data(i,j,k,n) = ( mask(i,j-1,k) * data(i,j-1,k,n) + data(i+1,j,k,n) + mask(i,j+1,k) * data(i,j+1,k,n) )
                       / ( mask(i,j-1,k) + 1 + mask(i,j+1,k) );
      // xhi, y-valid
      } else if ( (i == hi.x+1) && (j >= lo.y) && (j <= hi.y) ) {
         data(i,j,k,n) = ( mask(i,j-1,k) * data(i,j-1,k,n) + data(i-1,j,k,n) + mask(i,j+1,k) * data(i,j+1,k,n) )
                       / ( mask(i,j-1,k) + 1 + mask(i,j+1,k) );
      // x-valid, ylo
      } else if ( (i >= lo.x) && (i <= hi.x) && (j == lo.y-1) ) {
         data(i,j,k,n) = ( mask(i-1,j,k) * data(i-1,j,k,n) + mask(i+1,j,k) * data(i+1,j,k,n) + data(i,j+1,k,n) )
                       / ( mask(i-1,j,k) + mask(i+1,j,k) + 1 );
      // x-valid, yhi
      } else if ( (i >= lo.x) && (i <= hi.x) && (j == hi.y+1) ) {
         data(i,j,k,n) = ( mask(i-1,j,k) * data(i-1,j,k,n) + mask(i+1,j,k) * data(i+1,j,k,n) + data(i,j-1,k,n) )
                         / ( mask(i-1,j,k) + mask(i+1,j,k) + 1 );
      }
   }
}

}
#endif
