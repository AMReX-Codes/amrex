#ifndef AMReX_extrapolater_3D_K_H_
#define AMReX_extrapolater_3D_K_H_

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

   for (int n = 0; n < nComp; n++) {

      // set all crse cells to zero first
      for (int k = lo.z-1; k <= hi.z+1; ++k) {
         for (int j = lo.y-1; j <= hi.y+1; ++j) {
            for (int i = lo.x-1; i <= hi.x+1; ++i) {
               if (mask(i,j,k) == crsecell) data(i,j,k,n) = Real(0.0);
            }
         }
      }

      // Corners
      // xlo, ylo, zlo
      {
         int i = lo.x-1;
         int j = lo.y-1;
         int k = lo.z-1;
         if ( mask(i,j,k) == crsecell ) {
            if ( ( mask(i+1,j,k) == finecell ) ||
                 ( mask(i,j+1,k) == finecell ) ||
                 ( mask(i,j,k+1) == finecell ) ) {
               data(i,j,k,n) = (  mask(i+1,j,k) * data(i+1,j,k,n) 
                                + mask(i,j+1,k) * data(i,j+1,k,n) 
                                + mask(i,j,k+1) * data(i,j,k+1,n) )
                               / ( mask(i+1,j,k) + mask(i,j+1,k) + mask(i,j,k+1) );
            } else if ( ( mask(i+1,j+1,k) == finecell ) ||
                        ( mask(i+1,j,k+1) == finecell ) ||
                        ( mask(i,j+1,k+1) == finecell ) ) {
               data(i,j,k,n) = (  mask(i+1,j+1,k) * data(i+1,j+1,k,n) 
                                + mask(i+1,j,k+1) * data(i+1,j,k+1,n) 
                                + mask(i,j+1,k+1) * data(i,j+1,k+1,n) )
                               / ( mask(i+1,j+1,k) + mask(i+1,j,k+1) + mask(i,j+1,k+1) );
            } else {
               data(i,j,k,n) = data(i+1,j+1,k+1,n);
            }
         }
      }
      // xlo, ylo, zhi
      {
         int i = lo.x-1;
         int j = lo.y-1;
         int k = hi.z+1;
         if ( mask(i,j,k) == crsecell ) {
            if ( ( mask(i+1,j,k) == finecell ) ||
                 ( mask(i,j+1,k) == finecell ) ||
                 ( mask(i,j,k-1) == finecell ) ) {
               data(i,j,k,n) = (  mask(i+1,j,k) * data(i+1,j,k,n) 
                                + mask(i,j+1,k) * data(i,j+1,k,n) 
                                + mask(i,j,k-1) * data(i,j,k-1,n) )
                               / ( mask(i+1,j,k) + mask(i,j+1,k) + mask(i,j,k-1) );
            } else if ( ( mask(i+1,j+1,k) == finecell ) ||
                        ( mask(i+1,j,k-1) == finecell ) ||
                        ( mask(i,j+1,k-1) == finecell ) ) {
               data(i,j,k,n) = (  mask(i+1,j+1,k) * data(i+1,j+1,k,n) 
                                + mask(i+1,j,k-1) * data(i+1,j,k-1,n) 
                                + mask(i,j+1,k-1) * data(i,j+1,k-1,n) )
                               / ( mask(i+1,j+1,k) + mask(i+1,j,k-1) + mask(i,j+1,k-1) );
            } else {
               data(i,j,k,n) = data(i+1,j+1,k-1,n);
            }
         }
      }
      // xlo, yhi, zlo
      {
         int i = lo.x-1;
         int j = hi.y+1;
         int k = lo.z-1;
         if ( mask(i,j,k) == crsecell ) {
            if ( ( mask(i+1,j,k) == finecell ) ||
                 ( mask(i,j-1,k) == finecell ) ||
                 ( mask(i,j,k+1) == finecell ) ) {
               data(i,j,k,n) = (  mask(i+1,j,k) * data(i+1,j,k,n) 
                                + mask(i,j-1,k) * data(i,j-1,k,n) 
                                + mask(i,j,k+1) * data(i,j,k+1,n) )
                               / ( mask(i+1,j,k) + mask(i,j-1,k) + mask(i,j,k+1) );
            } else if ( ( mask(i+1,j-1,k) == finecell ) ||
                        ( mask(i+1,j,k+1) == finecell ) ||
                        ( mask(i,j-1,k+1) == finecell ) ) {
               data(i,j,k,n) = (  mask(i+1,j-1,k) * data(i+1,j-1,k,n) 
                                + mask(i+1,j,k+1) * data(i+1,j,k+1,n) 
                                + mask(i,j-1,k+1) * data(i,j-1,k+1,n) )
                               / ( mask(i+1,j-1,k) + mask(i+1,j,k+1) + mask(i,j-1,k+1) );
            } else {
               data(i,j,k,n) = data(i+1,j-1,k+1,n);
            }
         }
      }
      // xlo, yhi, zhi
      {
         int i = lo.x-1;
         int j = hi.y+1;
         int k = hi.z+1;
         if ( mask(i,j,k) == crsecell ) {
            if ( ( mask(i+1,j,k) == finecell ) ||
                 ( mask(i,j-1,k) == finecell ) ||
                 ( mask(i,j,k-1) == finecell ) ) {
               data(i,j,k,n) = (  mask(i+1,j,k) * data(i+1,j,k,n) 
                                + mask(i,j-1,k) * data(i,j-1,k,n) 
                                + mask(i,j,k-1) * data(i,j,k-1,n) )
                               / ( mask(i+1,j,k) + mask(i,j-1,k) + mask(i,j,k-1) );
            } else if ( ( mask(i+1,j-1,k) == finecell ) ||
                        ( mask(i+1,j,k-1) == finecell ) ||
                        ( mask(i,j-1,k-1) == finecell ) ) {
               data(i,j,k,n) = (  mask(i+1,j-1,k) * data(i+1,j-1,k,n) 
                                + mask(i+1,j,k-1) * data(i+1,j,k-1,n) 
                                + mask(i,j-1,k-1) * data(i,j-1,k-1,n) )
                               / ( mask(i+1,j-1,k) + mask(i+1,j,k-1) + mask(i,j-1,k-1) );
            } else {
               data(i,j,k,n) = data(i+1,j-1,k-1,n);
            }
         }
      }
      // xhi, ylo, zlo
      {
         int i = hi.x+1;
         int j = lo.y-1;
         int k = lo.z-1;
         if ( mask(i,j,k) == crsecell ) {
            if ( ( mask(i-1,j,k) == finecell ) ||
                 ( mask(i,j+1,k) == finecell ) ||
                 ( mask(i,j,k+1) == finecell ) ) {
               data(i,j,k,n) = (  mask(i-1,j,k) * data(i-1,j,k,n) 
                                + mask(i,j+1,k) * data(i,j+1,k,n) 
                                + mask(i,j,k+1) * data(i,j,k+1,n) )
                               / ( mask(i-1,j,k) + mask(i,j+1,k) + mask(i,j,k+1) );
            } else if ( ( mask(i-1,j+1,k) == finecell ) ||
                        ( mask(i-1,j,k+1) == finecell ) ||
                        ( mask(i,j+1,k+1) == finecell ) ) {
               data(i,j,k,n) = (  mask(i-1,j+1,k) * data(i-1,j+1,k,n) 
                                + mask(i-1,j,k+1) * data(i-1,j,k+1,n) 
                                + mask(i,j+1,k+1) * data(i,j+1,k+1,n) )
                               / ( mask(i-1,j+1,k) + mask(i-1,j,k+1) + mask(i,j+1,k+1) );
            } else {
               data(i,j,k,n) = data(i-1,j+1,k+1,n);
            }
         }
      }
      // xhi, ylo, zhi
      {
         int i = hi.x+1;
         int j = lo.y-1;
         int k = hi.z+1;
         if ( mask(i,j,k) == crsecell ) {
            if ( ( mask(i-1,j,k) == finecell ) ||
                 ( mask(i,j+1,k) == finecell ) ||
                 ( mask(i,j,k-1) == finecell ) ) {
               data(i,j,k,n) = (  mask(i-1,j,k) * data(i-1,j,k,n) 
                                + mask(i,j+1,k) * data(i,j+1,k,n) 
                                + mask(i,j,k-1) * data(i,j,k-1,n) )
                               / ( mask(i-1,j,k) + mask(i,j+1,k) + mask(i,j,k-1) );
            } else if ( ( mask(i-1,j+1,k) == finecell ) ||
                        ( mask(i-1,j,k-1) == finecell ) ||
                        ( mask(i,j+1,k-1) == finecell ) ) {
               data(i,j,k,n) = (  mask(i-1,j+1,k) * data(i-1,j+1,k,n) 
                                + mask(i-1,j,k-1) * data(i-1,j,k-1,n) 
                                + mask(i,j+1,k-1) * data(i,j+1,k-1,n) )
                               / ( mask(i-1,j+1,k) + mask(i-1,j,k-1) + mask(i,j+1,k-1) );
            } else {
               data(i,j,k,n) = data(i-1,j+1,k-1,n);
            }
         }
      }
      // xhi, yhi, zlo
      {
         int i = hi.x+1;
         int j = hi.y+1;
         int k = lo.z-1;
         if ( mask(i,j,k) == crsecell ) {
            if ( ( mask(i-1,j,k) == finecell ) ||
                 ( mask(i,j-1,k) == finecell ) ||
                 ( mask(i,j,k+1) == finecell ) ) {
               data(i,j,k,n) = (  mask(i-1,j,k) * data(i-1,j,k,n) 
                                + mask(i,j-1,k) * data(i,j-1,k,n) 
                                + mask(i,j,k+1) * data(i,j,k+1,n) )
                               / ( mask(i-1,j,k) + mask(i,j-1,k) + mask(i,j,k+1) );
            } else if ( ( mask(i-1,j-1,k) == finecell ) ||
                        ( mask(i-1,j,k+1) == finecell ) ||
                        ( mask(i,j-1,k+1) == finecell ) ) {
               data(i,j,k,n) = (  mask(i-1,j-1,k) * data(i-1,j-1,k,n) 
                                + mask(i-1,j,k+1) * data(i-1,j,k+1,n) 
                                + mask(i,j-1,k+1) * data(i,j-1,k+1,n) )
                               / ( mask(i-1,j-1,k) + mask(i-1,j,k+1) + mask(i,j-1,k+1) );
            } else {
               data(i,j,k,n) = data(i-1,j-1,k+1,n);
            }
         }
      }
      // xhi, yhi, zhi
      {
         int i = hi.x+1;
         int j = hi.y+1;
         int k = hi.z+1;
         if ( mask(i,j,k) == crsecell ) {
            if ( ( mask(i-1,j,k) == finecell ) ||
                 ( mask(i,j-1,k) == finecell ) ||
                 ( mask(i,j,k-1) == finecell ) ) {
               data(i,j,k,n) = (  mask(i-1,j,k) * data(i-1,j,k,n) 
                                + mask(i,j-1,k) * data(i,j-1,k,n) 
                                + mask(i,j,k-1) * data(i,j,k-1,n) )
                               / ( mask(i-1,j,k) + mask(i,j-1,k) + mask(i,j,k-1) );
            } else if ( ( mask(i-1,j-1,k) == finecell ) ||
                        ( mask(i-1,j,k-1) == finecell ) ||
                        ( mask(i,j-1,k-1) == finecell ) ) {
               data(i,j,k,n) = (  mask(i-1,j-1,k) * data(i-1,j-1,k,n) 
                                + mask(i-1,j,k-1) * data(i-1,j,k-1,n) 
                                + mask(i,j-1,k-1) * data(i,j-1,k-1,n) )
                               / ( mask(i-1,j-1,k) + mask(i-1,j,k-1) + mask(i,j-1,k-1) );
            } else {
               data(i,j,k,n) = data(i-1,j-1,k-1,n);
            }
         }
      }
      // Edges
      // xlo, ylo, z-valid
      {
         int i = lo.x-1;
         int j = lo.y-1;
         for (int k = lo.z; k <= hi.z; ++k) {
            if ( mask(i,j,k) == crsecell ) {
               if ( ( mask(i+1,j,k) == finecell ) ||
                    ( mask(i,j+1,k) == finecell ) ||
                    ( mask(i,j,k-1) == finecell ) ||
                    ( mask(i,j,k+1) == finecell ) ) {
                  data(i,j,k,n) = (  mask(i+1,j,k) * data(i+1,j,k,n)
                                   + mask(i,j+1,k) * data(i,j+1,k,n) 
                                   + mask(i,j,k-1) * data(i,j,k-1,n) 
                                   + mask(i,j,k+1) * data(i,j,k+1,n) ) 
                                 / ( mask(i+1,j,k) + mask(i,j+1,k) + mask(i,j,k-1) + mask(i,j,k+1) );
               } else {
                  data(i,j,k,n) = (                    data(i+1,j+1,k,n)
                                   + mask(i+1,j,k-1) * data(i+1,j,k-1,n)
                                   + mask(i+1,j,k+1) * data(i+1,j,k+1,n)
                                   + mask(i,j+1,k-1) * data(i,j+1,k-1,n)
                                   + mask(i,j+1,k+1) * data(i,j+1,k+1,n) )
                                 / ( 1.0 + mask(i+1,j,k-1) + mask(i+1,j,k+1) + mask(i,j+1,k-1) + mask(i,j+1,k+1) );
               }
            }
         }
      }
      // xlo, yhi, z-valid
      {
         int i = lo.x-1;
         int j = hi.y+1;
         for (int k = lo.z; k <= hi.z; ++k) {
            if ( mask(i,j,k) == crsecell ) {
               if ( ( mask(i+1,j,k) == finecell ) ||
                    ( mask(i,j-1,k) == finecell ) ||
                    ( mask(i,j,k-1) == finecell ) ||
                    ( mask(i,j,k+1) == finecell ) ) {
                  data(i,j,k,n) = (  mask(i+1,j,k) * data(i+1,j,k,n)
                                   + mask(i,j-1,k) * data(i,j-1,k,n) 
                                   + mask(i,j,k-1) * data(i,j,k-1,n) 
                                   + mask(i,j,k+1) * data(i,j,k+1,n) ) 
                                 / ( mask(i+1,j,k) + mask(i,j-1,k) + mask(i,j,k-1) + mask(i,j,k+1) );
               } else {
                  data(i,j,k,n) = (                    data(i+1,j-1,k,n)
                                   + mask(i+1,j,k-1) * data(i+1,j,k-1,n)
                                   + mask(i+1,j,k+1) * data(i+1,j,k+1,n)
                                   + mask(i,j-1,k-1) * data(i,j-1,k-1,n)
                                   + mask(i,j-1,k+1) * data(i,j-1,k+1,n) )
                                 / ( 1.0 + mask(i+1,j,k-1) + mask(i+1,j,k+1) + mask(i,j-1,k-1) + mask(i,j-1,k+1) );
               }
            }
         }
      }
      // xhi, ylo, z-valid
      {
         int i = hi.x+1;
         int j = lo.y-1;
         for (int k = lo.z; k <= hi.z; ++k) {
            if ( mask(i,j,k) == crsecell ) {
               if ( ( mask(i-1,j,k) == finecell ) ||
                    ( mask(i,j+1,k) == finecell ) ||
                    ( mask(i,j,k-1) == finecell ) ||
                    ( mask(i,j,k+1) == finecell ) ) {
                  data(i,j,k,n) = (  mask(i-1,j,k) * data(i-1,j,k,n)
                                   + mask(i,j+1,k) * data(i,j+1,k,n) 
                                   + mask(i,j,k-1) * data(i,j,k-1,n) 
                                   + mask(i,j,k+1) * data(i,j,k+1,n) ) 
                                 / ( mask(i-1,j,k) + mask(i,j+1,k) + mask(i,j,k-1) + mask(i,j,k+1) );
               } else {
                  data(i,j,k,n) = (                    data(i-1,j+1,k,n)
                                   + mask(i-1,j,k-1) * data(i-1,j,k-1,n)
                                   + mask(i-1,j,k+1) * data(i-1,j,k+1,n)
                                   + mask(i,j+1,k-1) * data(i,j+1,k-1,n)
                                   + mask(i,j+1,k+1) * data(i,j+1,k+1,n) )
                                 / ( 1.0 + mask(i-1,j,k-1) + mask(i-1,j,k+1) + mask(i,j+1,k-1) + mask(i,j+1,k+1) );
               }
            }
         }
      }
      // xhi, yhi, z-valid
      {
         int i = hi.x+1;
         int j = hi.y+1;
         for (int k = lo.z; k <= hi.z; ++k) {
            if ( mask(i,j,k) == crsecell ) {
               if ( ( mask(i-1,j,k) == finecell ) ||
                    ( mask(i,j-1,k) == finecell ) ||
                    ( mask(i,j,k-1) == finecell ) ||
                    ( mask(i,j,k+1) == finecell ) ) {
                  data(i,j,k,n) = (  mask(i-1,j,k) * data(i-1,j,k,n)
                                   + mask(i,j-1,k) * data(i,j-1,k,n) 
                                   + mask(i,j,k-1) * data(i,j,k-1,n) 
                                   + mask(i,j,k+1) * data(i,j,k+1,n) ) 
                                 / ( mask(i-1,j,k) + mask(i,j-1,k) + mask(i,j,k-1) + mask(i,j,k+1) );
               } else {
                  data(i,j,k,n) = (                    data(i-1,j-1,k,n)
                                   + mask(i-1,j,k-1) * data(i-1,j,k-1,n)
                                   + mask(i-1,j,k+1) * data(i-1,j,k+1,n)
                                   + mask(i,j-1,k-1) * data(i,j-1,k-1,n)
                                   + mask(i,j-1,k+1) * data(i,j-1,k+1,n) )
                                 / ( 1.0 + mask(i-1,j,k-1) + mask(i-1,j,k+1) + mask(i,j-1,k-1) + mask(i,j-1,k+1) );
               }
            }
         }
      }
      // xlo, y-valid, zlo
      {
         int i = lo.x-1;
         int k = lo.z-1;
         for (int j = lo.y; j <= hi.y; ++j) {
            if ( mask(i,j,k) == crsecell ) {
               if ( ( mask(i+1,j,k) == finecell ) ||
                    ( mask(i,j-1,k) == finecell ) ||
                    ( mask(i,j+1,k) == finecell ) ||
                    ( mask(i,j,k+1) == finecell ) ) {
                  data(i,j,k,n) = (  mask(i+1,j,k) * data(i+1,j,k,n)
                                   + mask(i,j-1,k) * data(i,j-1,k,n) 
                                   + mask(i,j+1,k) * data(i,j+1,k,n) 
                                   + mask(i,j,k+1) * data(i,j,k+1,n) ) 
                                 / ( mask(i+1,j,k) + mask(i,j-1,k) + mask(i,j+1,k) + mask(i,j,k+1) );
               } else {
                  data(i,j,k,n) = (  mask(i+1,j-1,k) * data(i+1,j-1,k,n)
                                   + mask(i+1,j+1,k) * data(i+1,j+1,k,n)
                                   +                   data(i+1,j,k+1,n)
                                   + mask(i,j-1,k+1) * data(i,j-1,k+1,n)
                                   + mask(i,j+1,k+1) * data(i,j+1,k+1,n) )
                                 / ( mask(i+1,j-1,k) + mask(i+1,j+1,k) + 1.0 + mask(i,j-1,k+1) + mask(i,j+1,k+1) );
               }
            }
         }
      }
      // xlo, y-valid, zhi
      {
         int i = lo.x-1;
         int k = hi.z+1;
         for (int j = lo.y; j <= hi.y; ++j) {
            if ( mask(i,j,k) == crsecell ) {
               if ( ( mask(i+1,j,k) == finecell ) ||
                    ( mask(i,j-1,k) == finecell ) ||
                    ( mask(i,j+1,k) == finecell ) ||
                    ( mask(i,j,k-1) == finecell ) ) {
                  data(i,j,k,n) = (  mask(i+1,j,k) * data(i+1,j,k,n)
                                   + mask(i,j-1,k) * data(i,j-1,k,n) 
                                   + mask(i,j+1,k) * data(i,j+1,k,n) 
                                   + mask(i,j,k-1) * data(i,j,k-1,n) ) 
                                 / ( mask(i+1,j,k) + mask(i,j-1,k) + mask(i,j+1,k) + mask(i,j,k-1) );
               } else {
                  data(i,j,k,n) = (  mask(i+1,j-1,k) * data(i+1,j-1,k,n)
                                   + mask(i+1,j+1,k) * data(i+1,j+1,k,n)
                                   +                   data(i+1,j,k-1,n)
                                   + mask(i,j-1,k-1) * data(i,j-1,k-1,n)
                                   + mask(i,j+1,k-1) * data(i,j+1,k-1,n) )
                                 / ( mask(i+1,j-1,k) + mask(i+1,j+1,k) + 1.0 + mask(i,j-1,k-1) + mask(i,j+1,k-1) );
               }
            }
         }
      }
      // xhi, y-valid, zlo
      {
         int i = hi.x+1;
         int k = lo.z-1;
         for (int j = lo.y; j <= hi.y; ++j) {
            if ( mask(i,j,k) == crsecell ) {
               if ( ( mask(i-1,j,k) == finecell ) ||
                    ( mask(i,j-1,k) == finecell ) ||
                    ( mask(i,j+1,k) == finecell ) ||
                    ( mask(i,j,k+1) == finecell ) ) {
                  data(i,j,k,n) = (  mask(i-1,j,k) * data(i-1,j,k,n)
                                   + mask(i,j-1,k) * data(i,j-1,k,n) 
                                   + mask(i,j+1,k) * data(i,j+1,k,n) 
                                   + mask(i,j,k+1) * data(i,j,k+1,n) ) 
                                 / ( mask(i-1,j,k) + mask(i,j-1,k) + mask(i,j+1,k) + mask(i,j,k+1) );
               } else {
                  data(i,j,k,n) = (  mask(i-1,j-1,k) * data(i-1,j-1,k,n)
                                   + mask(i-1,j+1,k) * data(i-1,j+1,k,n)
                                   +                   data(i-1,j,k+1,n)
                                   + mask(i,j-1,k+1) * data(i,j-1,k+1,n)
                                   + mask(i,j+1,k+1) * data(i,j+1,k+1,n) )
                                 / ( mask(i-1,j-1,k) + mask(i-1,j+1,k) + 1.0 + mask(i,j-1,k+1) + mask(i,j+1,k+1) );
               }
            }
         }
      }
      // xhi, y-valid, zhi
      {
         int i = hi.x+1;
         int k = hi.z+1;
         for (int j = lo.y; j <= hi.y; ++j) {
            if ( mask(i,j,k) == crsecell ) {
               if ( ( mask(i-1,j,k) == finecell ) ||
                    ( mask(i,j-1,k) == finecell ) ||
                    ( mask(i,j+1,k) == finecell ) ||
                    ( mask(i,j,k-1) == finecell ) ) {
                  data(i,j,k,n) = (  mask(i-1,j,k) * data(i-1,j,k,n)
                                   + mask(i,j-1,k) * data(i,j-1,k,n) 
                                   + mask(i,j+1,k) * data(i,j+1,k,n) 
                                   + mask(i,j,k-1) * data(i,j,k-1,n) ) 
                                 / ( mask(i-1,j,k) + mask(i,j-1,k) + mask(i,j+1,k) + mask(i,j,k-1) );
               } else {
                  data(i,j,k,n) = (  mask(i-1,j-1,k) * data(i-1,j-1,k,n)
                                   + mask(i-1,j+1,k) * data(i-1,j+1,k,n)
                                   +                   data(i-1,j,k-1,n)
                                   + mask(i,j-1,k-1) * data(i,j-1,k-1,n)
                                   + mask(i,j+1,k-1) * data(i,j+1,k-1,n) )
                                 / ( mask(i-1,j-1,k) + mask(i-1,j+1,k) + 1.0 + mask(i,j-1,k-1) + mask(i,j+1,k-1) );
               }
            }
         }
      }
      // x-valid, ylo, zlo
      {
         int j = lo.y-1;
         int k = lo.z-1;
         for (int i = lo.x; i <= hi.x; ++i) {
            if ( mask(i,j,k) == crsecell ) {
               if ( ( mask(i-1,j,k) == finecell ) ||
                    ( mask(i+1,j,k) == finecell ) ||
                    ( mask(i,j+1,k) == finecell ) ||
                    ( mask(i,j,k+1) == finecell ) ) {
                  data(i,j,k,n) = (  mask(i-1,j,k) * data(i-1,j,k,n)
                                   + mask(i+1,j,k) * data(i+1,j,k,n) 
                                   + mask(i,j+1,k) * data(i,j+1,k,n) 
                                   + mask(i,j,k+1) * data(i,j,k+1,n) ) 
                                 / ( mask(i-1,j,k) + mask(i+1,j,k) + mask(i,j+1,k) + mask(i,j,k+1) );
               } else {
                  data(i,j,k,n) = (  mask(i-1,j+1,k) * data(i-1,j+1,k,n)
                                   + mask(i+1,j+1,k) * data(i+1,j+1,k,n)
                                   + mask(i-1,j,k+1) * data(i-1,j,k+1,n)
                                   + mask(i+1,j,k+1) * data(i+1,j,k+1,n)
                                   +                   data(i,j+1,k+1,n) )
                                 / ( mask(i-1,j+1,k) + mask(i+1,j+1,k) + mask(i-1,j,k+1) + mask(i+1,j,k+1) + 1.0 );
               }
            }
         }
      }
      // x-valid, ylo, zhi
      {
         int j = lo.y-1;
         int k = hi.z+1;
         for (int i = lo.x; i <= hi.x; ++i) {
            if ( mask(i,j,k) == crsecell ) {
               if ( ( mask(i-1,j,k) == finecell ) ||
                    ( mask(i+1,j,k) == finecell ) ||
                    ( mask(i,j+1,k) == finecell ) ||
                    ( mask(i,j,k-1) == finecell ) ) {
                  data(i,j,k,n) = (  mask(i-1,j,k) * data(i-1,j,k,n)
                                   + mask(i+1,j,k) * data(i+1,j,k,n) 
                                   + mask(i,j+1,k) * data(i,j+1,k,n) 
                                   + mask(i,j,k-1) * data(i,j,k-1,n) ) 
                                 / ( mask(i-1,j,k) + mask(i+1,j,k) + mask(i,j+1,k) + mask(i,j,k-1) );
               } else {
                  data(i,j,k,n) = (  mask(i-1,j+1,k) * data(i-1,j+1,k,n)
                                   + mask(i+1,j+1,k) * data(i+1,j+1,k,n)
                                   + mask(i-1,j,k-1) * data(i-1,j,k-1,n)
                                   + mask(i+1,j,k-1) * data(i+1,j,k-1,n)
                                   +                   data(i,j+1,k-1,n) )
                                 / ( mask(i-1,j+1,k) + mask(i+1,j+1,k) + mask(i-1,j,k-1) + mask(i+1,j,k-1) + 1.0 );
               }
            }
         }
      }
      // x-valid, yhi, zlo
      {
         int j = hi.y+1;
         int k = lo.z-1;
         for (int i = lo.x; i <= hi.x; ++i) {
            if ( mask(i,j,k) == crsecell ) {
               if ( ( mask(i-1,j,k) == finecell ) ||
                    ( mask(i+1,j,k) == finecell ) ||
                    ( mask(i,j-1,k) == finecell ) ||
                    ( mask(i,j,k+1) == finecell ) ) {
                  data(i,j,k,n) = (  mask(i-1,j,k) * data(i-1,j,k,n)
                                   + mask(i+1,j,k) * data(i+1,j,k,n) 
                                   + mask(i,j-1,k) * data(i,j-1,k,n) 
                                   + mask(i,j,k+1) * data(i,j,k+1,n) ) 
                                 / ( mask(i-1,j,k) + mask(i+1,j,k) + mask(i,j-1,k) + mask(i,j,k+1) );
               } else {
                  data(i,j,k,n) = (  mask(i-1,j-1,k) * data(i-1,j-1,k,n)
                                   + mask(i+1,j-1,k) * data(i+1,j-1,k,n)
                                   + mask(i-1,j,k+1) * data(i-1,j,k+1,n)
                                   + mask(i+1,j,k+1) * data(i+1,j,k+1,n)
                                   +                   data(i,j-1,k+1,n) )
                                 / ( mask(i-1,j-1,k) + mask(i+1,j-1,k) + mask(i-1,j,k+1) + mask(i+1,j,k+1) + 1.0 );
               }
            }
         }
      }
      // x-valid, yhi, zhi
      {
         int j = hi.y+1;
         int k = hi.z+1;
         for (int i = lo.x; i <= hi.x; ++i) {
            if ( mask(i,j,k) == crsecell ) {
               if ( ( mask(i-1,j,k) == finecell ) ||
                    ( mask(i+1,j,k) == finecell ) ||
                    ( mask(i,j-1,k) == finecell ) ||
                    ( mask(i,j,k-1) == finecell ) ) {
                  data(i,j,k,n) = (  mask(i-1,j,k) * data(i-1,j,k,n)
                                   + mask(i+1,j,k) * data(i+1,j,k,n) 
                                   + mask(i,j-1,k) * data(i,j-1,k,n) 
                                   + mask(i,j,k-1) * data(i,j,k-1,n) ) 
                                 / ( mask(i-1,j,k) + mask(i+1,j,k) + mask(i,j-1,k) + mask(i,j,k-1) );
               } else {
                  data(i,j,k,n) = (  mask(i-1,j-1,k) * data(i-1,j-1,k,n)
                                   + mask(i+1,j-1,k) * data(i+1,j-1,k,n)
                                   + mask(i-1,j,k-1) * data(i-1,j,k-1,n)
                                   + mask(i+1,j,k-1) * data(i+1,j,k-1,n)
                                   +                   data(i,j-1,k-1,n) )
                                 / ( mask(i-1,j-1,k) + mask(i+1,j-1,k) + mask(i-1,j,k-1) + mask(i+1,j,k-1) + 1.0 );
               }
            }
         }
      }
      // Faces
      // xlo, y-valid, z-valid
      {
         int i = lo.x-1;
         for (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
               if ( mask(i,j,k) == crsecell ) {
                  data(i,j,k,n) = (                  data(i+1,j,k,n)
                                   + mask(i,j-1,k) * data(i,j-1,k,n)
                                   + mask(i,j+1,k) * data(i,j+1,k,n)
                                   + mask(i,j,k-1) * data(i,j,k-1,n)
                                   + mask(i,j,k+1) * data(i,j,k+1,n) )
                                 / ( 1.0 + mask(i,j-1,k) + mask(i,j+1,k) + mask(i,j,k-1) + mask(i,j,k+1) );
               }
            }
         }
      }
      // xhi, y-valid, z-valid
      {
         int i = hi.x+1;
         for (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
               if ( mask(i,j,k) == crsecell ) {
                  data(i,j,k,n) = (                  data(i-1,j,k,n)
                                   + mask(i,j-1,k) * data(i,j-1,k,n)
                                   + mask(i,j+1,k) * data(i,j+1,k,n)
                                   + mask(i,j,k-1) * data(i,j,k-1,n)
                                   + mask(i,j,k+1) * data(i,j,k+1,n) )
                                 / ( 1.0 + mask(i,j-1,k) + mask(i,j+1,k) + mask(i,j,k-1) + mask(i,j,k+1) );
               }
            }
         }
      }
      // x-valid, ylo, z-valid
      {
         int j = lo.y-1;
         for (int k = lo.z; k <= hi.z; ++k) {
            for (int i = lo.x; i <= hi.x; ++i) {
               if ( mask(i,j,k) == crsecell ) {
                  data(i,j,k,n) = (  mask(i-1,j,k) * data(i-1,j,k,n)
                                   + mask(i+1,j,k) * data(i+1,j,k,n)
                                   +                 data(i,j+1,k,n)
                                   + mask(i,j,k-1) * data(i,j,k-1,n)
                                   + mask(i,j,k+1) * data(i,j,k+1,n) )
                                 / ( mask(i-1,j,k) + mask(i+1,j,k) + 1.0 + mask(i,j,k-1) + mask(i,j,k+1) );
               }
            }
         }
      }
      // x-valid, yhi, z-valid
      {
         int j = hi.y+1;
         for (int k = lo.z; k <= hi.z; ++k) {
            for (int i = lo.x; i <= hi.x; ++i) {
               if ( mask(i,j,k) == crsecell ) {
                  data(i,j,k,n) = (  mask(i-1,j,k) * data(i-1,j,k,n)
                                   + mask(i+1,j,k) * data(i+1,j,k,n)
                                   +                 data(i,j-1,k,n)
                                   + mask(i,j,k-1) * data(i,j,k-1,n)
                                   + mask(i,j,k+1) * data(i,j,k+1,n) )
                                 / ( mask(i-1,j,k) + mask(i+1,j,k) + 1.0 + mask(i,j,k-1) + mask(i,j,k+1) );
               }
            }
         }
      }
      // x-valid, y-valid, zlo
      {
         int k = lo.z-1;
         for (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
               if ( mask(i,j,k) == crsecell ) {
                  data(i,j,k,n) = (  mask(i-1,j,k) * data(i-1,j,k,n)
                                   + mask(i+1,j,k) * data(i+1,j,k,n)
                                   + mask(i,j-1,k) * data(i,j-1,k,n)
                                   + mask(i,j+1,k) * data(i,j+1,k,n)
                                   +                 data(i,j,k+1,n) )
                                 / ( mask(i-1,j,k) + mask(i+1,j,k) + mask(i,j-1,k) + mask(i,j+1,k) + 1.0 );
               }
            }
         }
      }
      // x-valid, y-valid, zhi
      {
         int k = hi.z+1;
         for (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
               if ( mask(i,j,k) == crsecell ) {
                  data(i,j,k,n) = (  mask(i-1,j,k) * data(i-1,j,k,n)
                                   + mask(i+1,j,k) * data(i+1,j,k,n)
                                   + mask(i,j-1,k) * data(i,j-1,k,n)
                                   + mask(i,j+1,k) * data(i,j+1,k,n)
                                   +                 data(i,j,k-1,n) )
                                 / ( mask(i-1,j,k) + mask(i+1,j,k) + mask(i,j-1,k) + mask(i,j+1,k) + 1.0 );
               }
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
      // xlo, ylo, zlo
      if ( ( i == lo.x-1) && ( j == lo.y-1 ) && ( k == lo.z-1 ) ) {  
         if ( ( mask(i+1,j,k) == finecell ) ||
              ( mask(i,j+1,k) == finecell ) ||
              ( mask(i,j,k+1) == finecell ) ) {
            data(i,j,k,n) = (  mask(i+1,j,k) * data(i+1,j,k,n) 
                             + mask(i,j+1,k) * data(i,j+1,k,n) 
                             + mask(i,j,k+1) * data(i,j,k+1,n) )
                            / ( mask(i+1,j,k) + mask(i,j+1,k) + mask(i,j,k+1) );
         } else if ( ( mask(i+1,j+1,k) == finecell ) ||
                     ( mask(i+1,j,k+1) == finecell ) ||
                     ( mask(i,j+1,k+1) == finecell ) ) {
            data(i,j,k,n) = (  mask(i+1,j+1,k) * data(i+1,j+1,k,n) 
                             + mask(i+1,j,k+1) * data(i+1,j,k+1,n) 
                             + mask(i,j+1,k+1) * data(i,j+1,k+1,n) )
                            / ( mask(i+1,j+1,k) + mask(i+1,j,k+1) + mask(i,j+1,k+1) );
         } else {
            data(i,j,k,n) = data(i+1,j+1,k+1,n);
         }
      // xlo, ylo, zhi
      } else if ( ( i == lo.x-1) && ( j == lo.y-1 ) && ( k == hi.z+1 ) ) {  
         if ( ( mask(i+1,j,k) == finecell ) ||
              ( mask(i,j+1,k) == finecell ) ||
              ( mask(i,j,k-1) == finecell ) ) {
            data(i,j,k,n) = (  mask(i+1,j,k) * data(i+1,j,k,n) 
                             + mask(i,j+1,k) * data(i,j+1,k,n) 
                             + mask(i,j,k-1) * data(i,j,k-1,n) )
                            / ( mask(i+1,j,k) + mask(i,j+1,k) + mask(i,j,k-1) );
         } else if ( ( mask(i+1,j+1,k) == finecell ) ||
                     ( mask(i+1,j,k-1) == finecell ) ||
                     ( mask(i,j+1,k-1) == finecell ) ) {
            data(i,j,k,n) = (  mask(i+1,j+1,k) * data(i+1,j+1,k,n) 
                             + mask(i+1,j,k-1) * data(i+1,j,k-1,n) 
                             + mask(i,j+1,k-1) * data(i,j+1,k-1,n) )
                            / ( mask(i+1,j+1,k) + mask(i+1,j,k-1) + mask(i,j+1,k-1) );
         } else {
            data(i,j,k,n) = data(i+1,j+1,k-1,n);
         }
      // xlo, yhi, zlo
      } else if ( ( i == lo.x-1) && ( j == hi.y+1 ) && ( k == lo.z-1 ) ) {  
         if ( ( mask(i+1,j,k) == finecell ) ||
              ( mask(i,j-1,k) == finecell ) ||
              ( mask(i,j,k+1) == finecell ) ) {
            data(i,j,k,n) = (  mask(i+1,j,k) * data(i+1,j,k,n) 
                             + mask(i,j-1,k) * data(i,j-1,k,n) 
                             + mask(i,j,k+1) * data(i,j,k+1,n) )
                            / ( mask(i+1,j,k) + mask(i,j-1,k) + mask(i,j,k+1) );
         } else if ( ( mask(i+1,j-1,k) == finecell ) ||
                     ( mask(i+1,j,k+1) == finecell ) ||
                     ( mask(i,j-1,k+1) == finecell ) ) {
            data(i,j,k,n) = (  mask(i+1,j-1,k) * data(i+1,j-1,k,n) 
                             + mask(i+1,j,k+1) * data(i+1,j,k+1,n) 
                             + mask(i,j-1,k+1) * data(i,j-1,k+1,n) )
                            / ( mask(i+1,j-1,k) + mask(i+1,j,k+1) + mask(i,j-1,k+1) );
         } else {
            data(i,j,k,n) = data(i+1,j-1,k+1,n);
         }
      // xlo, yhi, zhi
      } else if ( ( i == lo.x-1) && ( j == hi.y+1 ) && ( k == hi.z+1 ) ) {  
         if ( ( mask(i+1,j,k) == finecell ) ||
              ( mask(i,j-1,k) == finecell ) ||
              ( mask(i,j,k-1) == finecell ) ) {
            data(i,j,k,n) = (  mask(i+1,j,k) * data(i+1,j,k,n) 
                             + mask(i,j-1,k) * data(i,j-1,k,n) 
                             + mask(i,j,k-1) * data(i,j,k-1,n) )
                            / ( mask(i+1,j,k) + mask(i,j-1,k) + mask(i,j,k-1) );
         } else if ( ( mask(i+1,j-1,k) == finecell ) ||
                     ( mask(i+1,j,k-1) == finecell ) ||
                     ( mask(i,j-1,k-1) == finecell ) ) {
            data(i,j,k,n) = (  mask(i+1,j-1,k) * data(i+1,j-1,k,n) 
                             + mask(i+1,j,k-1) * data(i+1,j,k-1,n) 
                             + mask(i,j-1,k-1) * data(i,j-1,k-1,n) )
                            / ( mask(i+1,j-1,k) + mask(i+1,j,k-1) + mask(i,j-1,k-1) );
         } else {
            data(i,j,k,n) = data(i+1,j-1,k-1,n);
         }
      // xhi, ylo, zlo
      } else if ( ( i == hi.x+1) && ( j == lo.y-1 ) && ( k == lo.z-1 ) ) {  
         if ( ( mask(i-1,j,k) == finecell ) ||
              ( mask(i,j+1,k) == finecell ) ||
              ( mask(i,j,k+1) == finecell ) ) {
            data(i,j,k,n) = (  mask(i-1,j,k) * data(i-1,j,k,n) 
                             + mask(i,j+1,k) * data(i,j+1,k,n) 
                             + mask(i,j,k+1) * data(i,j,k+1,n) )
                            / ( mask(i-1,j,k) + mask(i,j+1,k) + mask(i,j,k+1) );
         } else if ( ( mask(i-1,j+1,k) == finecell ) ||
                     ( mask(i-1,j,k+1) == finecell ) ||
                     ( mask(i,j+1,k+1) == finecell ) ) {
            data(i,j,k,n) = (  mask(i-1,j+1,k) * data(i-1,j+1,k,n) 
                             + mask(i-1,j,k+1) * data(i-1,j,k+1,n) 
                             + mask(i,j+1,k+1) * data(i,j+1,k+1,n) )
                            / ( mask(i-1,j+1,k) + mask(i-1,j,k+1) + mask(i,j+1,k+1) );
         } else {
            data(i,j,k,n) = data(i-1,j+1,k+1,n);
         }
      // xhi, ylo, zhi
      } else if ( ( i == hi.x+1) && ( j == lo.y-1 ) && ( k == hi.z+1 ) ) {  
         if ( ( mask(i-1,j,k) == finecell ) ||
              ( mask(i,j+1,k) == finecell ) ||
              ( mask(i,j,k-1) == finecell ) ) {
            data(i,j,k,n) = (  mask(i-1,j,k) * data(i-1,j,k,n) 
                             + mask(i,j+1,k) * data(i,j+1,k,n) 
                             + mask(i,j,k-1) * data(i,j,k-1,n) )
                            / ( mask(i-1,j,k) + mask(i,j+1,k) + mask(i,j,k-1) );
         } else if ( ( mask(i-1,j+1,k) == finecell ) ||
                     ( mask(i-1,j,k-1) == finecell ) ||
                     ( mask(i,j+1,k-1) == finecell ) ) {
            data(i,j,k,n) = (  mask(i-1,j+1,k) * data(i-1,j+1,k,n) 
                             + mask(i-1,j,k-1) * data(i-1,j,k-1,n) 
                             + mask(i,j+1,k-1) * data(i,j+1,k-1,n) )
                            / ( mask(i-1,j+1,k) + mask(i-1,j,k-1) + mask(i,j+1,k-1) );
         } else {
            data(i,j,k,n) = data(i-1,j+1,k-1,n);
         }
      // xhi, yhi, zlo
      } else if ( ( i == hi.x+1) && ( j == hi.y+1 ) && ( k == lo.z-1 ) ) {  
         if ( ( mask(i-1,j,k) == finecell ) ||
              ( mask(i,j-1,k) == finecell ) ||
              ( mask(i,j,k+1) == finecell ) ) {
            data(i,j,k,n) = (  mask(i-1,j,k) * data(i-1,j,k,n) 
                             + mask(i,j-1,k) * data(i,j-1,k,n) 
                             + mask(i,j,k+1) * data(i,j,k+1,n) )
                            / ( mask(i-1,j,k) + mask(i,j-1,k) + mask(i,j,k+1) );
         } else if ( ( mask(i-1,j-1,k) == finecell ) ||
                     ( mask(i-1,j,k+1) == finecell ) ||
                     ( mask(i,j-1,k+1) == finecell ) ) {
            data(i,j,k,n) = (  mask(i-1,j-1,k) * data(i-1,j-1,k,n) 
                             + mask(i-1,j,k+1) * data(i-1,j,k+1,n) 
                             + mask(i,j-1,k+1) * data(i,j-1,k+1,n) )
                            / ( mask(i-1,j-1,k) + mask(i-1,j,k+1) + mask(i,j-1,k+1) );
         } else {
            data(i,j,k,n) = data(i-1,j-1,k+1,n);
         }
      // xhi, yhi, zhi
      } else if ( ( i == hi.x+1) && ( j == hi.y+1 ) && ( k == hi.z+1 ) ) {  
         if ( ( mask(i-1,j,k) == finecell ) ||
              ( mask(i,j-1,k) == finecell ) ||
              ( mask(i,j,k-1) == finecell ) ) {
            data(i,j,k,n) = (  mask(i-1,j,k) * data(i-1,j,k,n) 
                             + mask(i,j-1,k) * data(i,j-1,k,n) 
                             + mask(i,j,k-1) * data(i,j,k-1,n) )
                            / ( mask(i-1,j,k) + mask(i,j-1,k) + mask(i,j,k-1) );
         } else if ( ( mask(i-1,j-1,k) == finecell ) ||
                     ( mask(i-1,j,k-1) == finecell ) ||
                     ( mask(i,j-1,k-1) == finecell ) ) {
            data(i,j,k,n) = (  mask(i-1,j-1,k) * data(i-1,j-1,k,n) 
                             + mask(i-1,j,k-1) * data(i-1,j,k-1,n) 
                             + mask(i,j-1,k-1) * data(i,j-1,k-1,n) )
                            / ( mask(i-1,j-1,k) + mask(i-1,j,k-1) + mask(i,j-1,k-1) );
         } else {
            data(i,j,k,n) = data(i-1,j-1,k-1,n);
         }
      // Edges
      // xlo, ylo, z-valid
      } else if ( ( i == lo.x-1) && ( j == lo.y-1 ) && 
                  ( k >= lo.z )  && ( k <= hi.z ) ) {
         if ( ( mask(i+1,j,k) == finecell ) ||
              ( mask(i,j+1,k) == finecell ) ||
              ( mask(i,j,k-1) == finecell ) ||
              ( mask(i,j,k+1) == finecell ) ) {
            data(i,j,k,n) = (  mask(i+1,j,k) * data(i+1,j,k,n)
                             + mask(i,j+1,k) * data(i,j+1,k,n) 
                             + mask(i,j,k-1) * data(i,j,k-1,n) 
                             + mask(i,j,k+1) * data(i,j,k+1,n) ) 
                           / ( mask(i+1,j,k) + mask(i,j+1,k) + mask(i,j,k-1) + mask(i,j,k+1) );
         } else {
            data(i,j,k,n) = (                    data(i+1,j+1,k,n)
                             + mask(i+1,j,k-1) * data(i+1,j,k-1,n)
                             + mask(i+1,j,k+1) * data(i+1,j,k+1,n)
                             + mask(i,j+1,k-1) * data(i,j+1,k-1,n)
                             + mask(i,j+1,k+1) * data(i,j+1,k+1,n) )
                           / ( 1.0 + mask(i+1,j,k-1) + mask(i+1,j,k+1) + mask(i,j+1,k-1) + mask(i,j+1,k+1) );
         }
      // xlo, yhi, z-valid
      } else if ( ( i == lo.x-1) && ( j == hi.y+1 ) && 
                  ( k >= lo.z )  && ( k <= hi.z ) ) {
         if ( ( mask(i+1,j,k) == finecell ) ||
              ( mask(i,j-1,k) == finecell ) ||
              ( mask(i,j,k-1) == finecell ) ||
              ( mask(i,j,k+1) == finecell ) ) {
            data(i,j,k,n) = (  mask(i+1,j,k) * data(i+1,j,k,n)
                             + mask(i,j-1,k) * data(i,j-1,k,n) 
                             + mask(i,j,k-1) * data(i,j,k-1,n) 
                             + mask(i,j,k+1) * data(i,j,k+1,n) ) 
                           / ( mask(i+1,j,k) + mask(i,j-1,k) + mask(i,j,k-1) + mask(i,j,k+1) );
         } else {
            data(i,j,k,n) = (                    data(i+1,j-1,k,n)
                             + mask(i+1,j,k-1) * data(i+1,j,k-1,n)
                             + mask(i+1,j,k+1) * data(i+1,j,k+1,n)
                             + mask(i,j-1,k-1) * data(i,j-1,k-1,n)
                             + mask(i,j-1,k+1) * data(i,j-1,k+1,n) )
                           / ( 1.0 + mask(i+1,j,k-1) + mask(i+1,j,k+1) + mask(i,j-1,k-1) + mask(i,j-1,k+1) );
         }
      // xhi, ylo, z-valid
      } else if ( ( i == hi.x+1) && ( j == lo.y-1 ) && 
                  ( k >= lo.z )  && ( k <= hi.z ) ) {
         if ( ( mask(i-1,j,k) == finecell ) ||
              ( mask(i,j+1,k) == finecell ) ||
              ( mask(i,j,k-1) == finecell ) ||
              ( mask(i,j,k+1) == finecell ) ) {
            data(i,j,k,n) = (  mask(i-1,j,k) * data(i-1,j,k,n)
                             + mask(i,j+1,k) * data(i,j+1,k,n) 
                             + mask(i,j,k-1) * data(i,j,k-1,n) 
                             + mask(i,j,k+1) * data(i,j,k+1,n) ) 
                           / ( mask(i-1,j,k) + mask(i,j+1,k) + mask(i,j,k-1) + mask(i,j,k+1) );
         } else {
            data(i,j,k,n) = (                    data(i-1,j+1,k,n)
                             + mask(i-1,j,k-1) * data(i-1,j,k-1,n)
                             + mask(i-1,j,k+1) * data(i-1,j,k+1,n)
                             + mask(i,j+1,k-1) * data(i,j+1,k-1,n)
                             + mask(i,j+1,k+1) * data(i,j+1,k+1,n) )
                           / ( 1.0 + mask(i-1,j,k-1) + mask(i-1,j,k+1) + mask(i,j+1,k-1) + mask(i,j+1,k+1) );
         }
      // xhi, yhi, z-valid
      } else if ( ( i == hi.x+1) && ( j == hi.y+1 ) && 
                  ( k >= lo.z )  && ( k <= hi.z ) ) {
         if ( ( mask(i-1,j,k) == finecell ) ||
              ( mask(i,j-1,k) == finecell ) ||
              ( mask(i,j,k-1) == finecell ) ||
              ( mask(i,j,k+1) == finecell ) ) {
            data(i,j,k,n) = (  mask(i-1,j,k) * data(i-1,j,k,n)
                             + mask(i,j-1,k) * data(i,j-1,k,n) 
                             + mask(i,j,k-1) * data(i,j,k-1,n) 
                             + mask(i,j,k+1) * data(i,j,k+1,n) ) 
                           / ( mask(i-1,j,k) + mask(i,j-1,k) + mask(i,j,k-1) + mask(i,j,k+1) );
         } else {
            data(i,j,k,n) = (                    data(i-1,j-1,k,n)
                             + mask(i-1,j,k-1) * data(i-1,j,k-1,n)
                             + mask(i-1,j,k+1) * data(i-1,j,k+1,n)
                             + mask(i,j-1,k-1) * data(i,j-1,k-1,n)
                             + mask(i,j-1,k+1) * data(i,j-1,k+1,n) )
                           / ( 1.0 + mask(i-1,j,k-1) + mask(i-1,j,k+1) + mask(i,j-1,k-1) + mask(i,j-1,k+1) );
         }
      // xlo, y-valid, zlo
      } else if ( ( i == lo.x-1) && ( j >= lo.y ) && 
                  ( j <= hi.y )  && ( k == lo.z-1 ) ) {
         if ( ( mask(i+1,j,k) == finecell ) ||
              ( mask(i,j-1,k) == finecell ) ||
              ( mask(i,j+1,k) == finecell ) ||
              ( mask(i,j,k+1) == finecell ) ) {
            data(i,j,k,n) = (  mask(i+1,j,k) * data(i+1,j,k,n)
                             + mask(i,j-1,k) * data(i,j-1,k,n) 
                             + mask(i,j+1,k) * data(i,j+1,k,n) 
                             + mask(i,j,k+1) * data(i,j,k+1,n) ) 
                           / ( mask(i+1,j,k) + mask(i,j-1,k) + mask(i,j+1,k) + mask(i,j,k+1) );
         } else {
            data(i,j,k,n) = (  mask(i+1,j-1,k) * data(i+1,j-1,k,n)
                             + mask(i+1,j+1,k) * data(i+1,j+1,k,n)
                             +                   data(i+1,j,k+1,n)
                             + mask(i,j-1,k+1) * data(i,j-1,k+1,n)
                             + mask(i,j+1,k+1) * data(i,j+1,k+1,n) )
                           / ( mask(i+1,j-1,k) + mask(i+1,j+1,k) + 1.0 + mask(i,j-1,k+1) + mask(i,j+1,k+1) );
         }
      // xlo, y-valid, zhi
      } else if ( ( i == lo.x-1) && ( j >= lo.y ) && 
                  ( j <= hi.y )  && ( k == hi.z+1 ) ) {
         if ( ( mask(i+1,j,k) == finecell ) ||
              ( mask(i,j-1,k) == finecell ) ||
              ( mask(i,j+1,k) == finecell ) ||
              ( mask(i,j,k-1) == finecell ) ) {
            data(i,j,k,n) = (  mask(i+1,j,k) * data(i+1,j,k,n)
                             + mask(i,j-1,k) * data(i,j-1,k,n) 
                             + mask(i,j+1,k) * data(i,j+1,k,n) 
                             + mask(i,j,k-1) * data(i,j,k-1,n) ) 
                           / ( mask(i+1,j,k) + mask(i,j-1,k) + mask(i,j+1,k) + mask(i,j,k-1) );
         } else {
            data(i,j,k,n) = (  mask(i+1,j-1,k) * data(i+1,j-1,k,n)
                             + mask(i+1,j+1,k) * data(i+1,j+1,k,n)
                             +                   data(i+1,j,k-1,n)
                             + mask(i,j-1,k-1) * data(i,j-1,k-1,n)
                             + mask(i,j+1,k-1) * data(i,j+1,k-1,n) )
                           / ( mask(i+1,j-1,k) + mask(i+1,j+1,k) + 1.0 + mask(i,j-1,k-1) + mask(i,j+1,k-1) );
         }
      // xhi, y-valid, zlo
      } else if ( ( i == hi.x+1) && ( j >= lo.y ) && 
                  ( j <= hi.y )  && ( k == lo.z-1 ) ) {
         if ( ( mask(i-1,j,k) == finecell ) ||
              ( mask(i,j-1,k) == finecell ) ||
              ( mask(i,j+1,k) == finecell ) ||
              ( mask(i,j,k+1) == finecell ) ) {
            data(i,j,k,n) = (  mask(i-1,j,k) * data(i-1,j,k,n)
                             + mask(i,j-1,k) * data(i,j-1,k,n) 
                             + mask(i,j+1,k) * data(i,j+1,k,n) 
                             + mask(i,j,k+1) * data(i,j,k+1,n) ) 
                           / ( mask(i-1,j,k) + mask(i,j-1,k) + mask(i,j+1,k) + mask(i,j,k+1) );
         } else {
            data(i,j,k,n) = (  mask(i-1,j-1,k) * data(i-1,j-1,k,n)
                             + mask(i-1,j+1,k) * data(i-1,j+1,k,n)
                             +                   data(i-1,j,k+1,n)
                             + mask(i,j-1,k+1) * data(i,j-1,k+1,n)
                             + mask(i,j+1,k+1) * data(i,j+1,k+1,n) )
                           / ( mask(i-1,j-1,k) + mask(i-1,j+1,k) + 1.0 + mask(i,j-1,k+1) + mask(i,j+1,k+1) );
         }
      // xhi, y-valid, zhi
      } else if ( ( i == hi.x+1) && ( j >= lo.y ) && 
                  ( j <= hi.y )  && ( k == hi.z+1 ) ) {
         if ( ( mask(i-1,j,k) == finecell ) ||
              ( mask(i,j-1,k) == finecell ) ||
              ( mask(i,j+1,k) == finecell ) ||
              ( mask(i,j,k-1) == finecell ) ) {
            data(i,j,k,n) = (  mask(i-1,j,k) * data(i-1,j,k,n)
                             + mask(i,j-1,k) * data(i,j-1,k,n) 
                             + mask(i,j+1,k) * data(i,j+1,k,n) 
                             + mask(i,j,k-1) * data(i,j,k-1,n) ) 
                           / ( mask(i-1,j,k) + mask(i,j-1,k) + mask(i,j+1,k) + mask(i,j,k-1) );
         } else {
            data(i,j,k,n) = (  mask(i-1,j-1,k) * data(i-1,j-1,k,n)
                             + mask(i-1,j+1,k) * data(i-1,j+1,k,n)
                             +                   data(i-1,j,k-1,n)
                             + mask(i,j-1,k-1) * data(i,j-1,k-1,n)
                             + mask(i,j+1,k-1) * data(i,j+1,k-1,n) )
                           / ( mask(i-1,j-1,k) + mask(i-1,j+1,k) + 1.0 + mask(i,j-1,k-1) + mask(i,j+1,k-1) );
         }
      // x-valid, ylo, zlo
      } else if ( ( i >= lo.x) && ( i <= hi.x ) && 
                  ( j == lo.y-1 )  && ( k == lo.z-1 ) ) {
         if ( ( mask(i-1,j,k) == finecell ) ||
              ( mask(i+1,j,k) == finecell ) ||
              ( mask(i,j+1,k) == finecell ) ||
              ( mask(i,j,k+1) == finecell ) ) {
            data(i,j,k,n) = (  mask(i-1,j,k) * data(i-1,j,k,n)
                             + mask(i+1,j,k) * data(i+1,j,k,n) 
                             + mask(i,j+1,k) * data(i,j+1,k,n) 
                             + mask(i,j,k+1) * data(i,j,k+1,n) ) 
                           / ( mask(i-1,j,k) + mask(i+1,j,k) + mask(i,j+1,k) + mask(i,j,k+1) );
         } else {
            data(i,j,k,n) = (  mask(i-1,j+1,k) * data(i-1,j+1,k,n)
                             + mask(i+1,j+1,k) * data(i+1,j+1,k,n)
                             + mask(i-1,j,k+1) * data(i-1,j,k+1,n)
                             + mask(i+1,j,k+1) * data(i+1,j,k+1,n)
                             +                   data(i,j+1,k+1,n) )
                           / ( mask(i-1,j+1,k) + mask(i+1,j+1,k) + mask(i-1,j,k+1) + mask(i+1,j,k+1) + 1.0 );
         }
      // x-valid, ylo, zhi
      } else if ( ( i >= lo.x) && ( i <= hi.x ) && 
                  ( j == lo.y-1 )  && ( k == hi.z+1 ) ) {
         if ( ( mask(i-1,j,k) == finecell ) ||
              ( mask(i+1,j,k) == finecell ) ||
              ( mask(i,j+1,k) == finecell ) ||
              ( mask(i,j,k-1) == finecell ) ) {
            data(i,j,k,n) = (  mask(i-1,j,k) * data(i-1,j,k,n)
                             + mask(i+1,j,k) * data(i+1,j,k,n) 
                             + mask(i,j+1,k) * data(i,j+1,k,n) 
                             + mask(i,j,k-1) * data(i,j,k-1,n) ) 
                           / ( mask(i-1,j,k) + mask(i+1,j,k) + mask(i,j+1,k) + mask(i,j,k-1) );
         } else {
            data(i,j,k,n) = (  mask(i-1,j+1,k) * data(i-1,j+1,k,n)
                             + mask(i+1,j+1,k) * data(i+1,j+1,k,n)
                             + mask(i-1,j,k-1) * data(i-1,j,k-1,n)
                             + mask(i+1,j,k-1) * data(i+1,j,k-1,n)
                             +                   data(i,j+1,k-1,n) )
                           / ( mask(i-1,j+1,k) + mask(i+1,j+1,k) + mask(i-1,j,k-1) + mask(i+1,j,k-1) + 1.0 );
         }
      // x-valid, yhi, zlo
      } else if ( ( i >= lo.x) && ( i <= hi.x ) && 
                  ( j == hi.y+1 )  && ( k == lo.z-1 ) ) {
         if ( ( mask(i-1,j,k) == finecell ) ||
              ( mask(i+1,j,k) == finecell ) ||
              ( mask(i,j-1,k) == finecell ) ||
              ( mask(i,j,k+1) == finecell ) ) {
            data(i,j,k,n) = (  mask(i-1,j,k) * data(i-1,j,k,n)
                             + mask(i+1,j,k) * data(i+1,j,k,n) 
                             + mask(i,j-1,k) * data(i,j-1,k,n) 
                             + mask(i,j,k+1) * data(i,j,k+1,n) ) 
                           / ( mask(i-1,j,k) + mask(i+1,j,k) + mask(i,j-1,k) + mask(i,j,k+1) );
         } else {
            data(i,j,k,n) = (  mask(i-1,j-1,k) * data(i-1,j-1,k,n)
                             + mask(i+1,j-1,k) * data(i+1,j-1,k,n)
                             + mask(i-1,j,k+1) * data(i-1,j,k+1,n)
                             + mask(i+1,j,k+1) * data(i+1,j,k+1,n)
                             +                   data(i,j-1,k+1,n) )
                           / ( mask(i-1,j-1,k) + mask(i+1,j-1,k) + mask(i-1,j,k+1) + mask(i+1,j,k+1) + 1.0 );
         }
      // x-valid, yhi, zhi
      } else if ( ( i >= lo.x) && ( i <= hi.x ) && 
                  ( j == hi.y+1 )  && ( k == hi.z+1 ) ) {
         if ( ( mask(i-1,j,k) == finecell ) ||
              ( mask(i+1,j,k) == finecell ) ||
              ( mask(i,j-1,k) == finecell ) ||
              ( mask(i,j,k-1) == finecell ) ) {
            data(i,j,k,n) = (  mask(i-1,j,k) * data(i-1,j,k,n)
                             + mask(i+1,j,k) * data(i+1,j,k,n) 
                             + mask(i,j-1,k) * data(i,j-1,k,n) 
                             + mask(i,j,k-1) * data(i,j,k-1,n) ) 
                           / ( mask(i-1,j,k) + mask(i+1,j,k) + mask(i,j-1,k) + mask(i,j,k-1) );
         } else {
            data(i,j,k,n) = (  mask(i-1,j-1,k) * data(i-1,j-1,k,n)
                             + mask(i+1,j-1,k) * data(i+1,j-1,k,n)
                             + mask(i-1,j,k-1) * data(i-1,j,k-1,n)
                             + mask(i+1,j,k-1) * data(i+1,j,k-1,n)
                             +                   data(i,j-1,k-1,n) )
                           / ( mask(i-1,j-1,k) + mask(i+1,j-1,k) + mask(i-1,j,k-1) + mask(i+1,j,k-1) + 1.0 );
         }
      // Faces
      // xlo, y-valid, z-valid
      } else if ( ( i == lo.x-1) &&
                  ( j >= lo.y )  && ( j <= hi.y ) &&
                  ( k >= lo.z )  && ( k <= hi.z ) ) {
         data(i,j,k,n) = (                  data(i+1,j,k,n)
                          + mask(i,j-1,k) * data(i,j-1,k,n)
                          + mask(i,j+1,k) * data(i,j+1,k,n)
                          + mask(i,j,k-1) * data(i,j,k-1,n)
                          + mask(i,j,k+1) * data(i,j,k+1,n) )
                        / ( 1.0 + mask(i,j-1,k) + mask(i,j+1,k) + mask(i,j,k-1) + mask(i,j,k+1) );
      // xhi, y-valid, z-valid
      } else if ( ( i == hi.x+1) &&
                  ( j >= lo.y )  && ( j <= hi.y ) &&
                  ( k >= lo.z )  && ( k <= hi.z ) ) {
         data(i,j,k,n) = (                  data(i-1,j,k,n)
                          + mask(i,j-1,k) * data(i,j-1,k,n)
                          + mask(i,j+1,k) * data(i,j+1,k,n)
                          + mask(i,j,k-1) * data(i,j,k-1,n)
                          + mask(i,j,k+1) * data(i,j,k+1,n) )
                        / ( 1.0 + mask(i,j-1,k) + mask(i,j+1,k) + mask(i,j,k-1) + mask(i,j,k+1) );
      // x-valid, ylo, z-valid
      } else if ( ( i >= lo.x )  && ( i <= hi.x ) &&
                  ( j == lo.y-1) &&
                  ( k >= lo.z )  && ( k <= hi.z ) ) {
         data(i,j,k,n) = (  mask(i-1,j,k) * data(i-1,j,k,n)
                          + mask(i+1,j,k) * data(i+1,j,k,n)
                          +                 data(i,j+1,k,n)
                          + mask(i,j,k-1) * data(i,j,k-1,n)
                          + mask(i,j,k+1) * data(i,j,k+1,n) )
                        / ( mask(i-1,j,k) + mask(i+1,j,k) + 1.0 + mask(i,j,k-1) + mask(i,j,k+1) );
      // x-valid, yhi, z-valid
      } else if ( ( i >= lo.x )  && ( i <= hi.x ) &&
                  ( j == hi.y+1) &&
                  ( k >= lo.z )  && ( k <= hi.z ) ) {
         data(i,j,k,n) = (  mask(i-1,j,k) * data(i-1,j,k,n)
                          + mask(i+1,j,k) * data(i+1,j,k,n)
                          +                 data(i,j-1,k,n)
                          + mask(i,j,k-1) * data(i,j,k-1,n)
                          + mask(i,j,k+1) * data(i,j,k+1,n) )
                        / ( mask(i-1,j,k) + mask(i+1,j,k) + 1.0 + mask(i,j,k-1) + mask(i,j,k+1) );
      // x-valid, y-valid, zlo
      } else if ( ( i >= lo.x )  && ( i <= hi.x ) &&
                  ( j >= lo.y )  && ( j <= hi.y ) &&
                  ( k == lo.z-1) ) {
         data(i,j,k,n) = (  mask(i-1,j,k) * data(i-1,j,k,n)
                          + mask(i+1,j,k) * data(i+1,j,k,n)
                          + mask(i,j-1,k) * data(i,j-1,k,n)
                          + mask(i,j+1,k) * data(i,j+1,k,n)
                          +                 data(i,j,k+1,n) )
                        / ( mask(i-1,j,k) + mask(i+1,j,k) + mask(i,j-1,k) + mask(i,j+1,k) + 1.0 );
      // x-valid, y-valid, zhi
      } else if ( ( i >= lo.x )  && ( i <= hi.x ) &&
                  ( j >= lo.y )  && ( j <= hi.y ) &&
                  ( k == hi.z+1) ) {
         data(i,j,k,n) = (  mask(i-1,j,k) * data(i-1,j,k,n)
                          + mask(i+1,j,k) * data(i+1,j,k,n)
                          + mask(i,j-1,k) * data(i,j-1,k,n)
                          + mask(i,j+1,k) * data(i,j+1,k,n)
                          +                 data(i,j,k-1,n) )
                        / ( mask(i-1,j,k) + mask(i+1,j,k) + mask(i,j-1,k) + mask(i,j+1,k) + 1.0 );
      }
   }
}

}
#endif
