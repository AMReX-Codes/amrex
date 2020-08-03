#ifndef AMReX_extrapolater_3D_K_H_
#define AMReX_extrapolater_3D_K_H_

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
   int idhi[3] = {ARMEX_D_DECL(hi.x,hi.y,hi.z)};
   int idlo_e[3] = {AMREX_D_DECL(lo.x,lo.y,lo.z)};          // extended box bounds
   int idhi_e[3] = {AMREX_D_DECL(hi.x,hi.y,hi.z)};
   for (int dim = 0; dim < AMREX_SPACEDIM; dim++) {   // Extend only AMREX_SPACEDIMs
      idlo_e[dim] -= 1;
      idhi_e[dim] += 1;
   }

   for (int n = 0; n < nComp; n++) {

      // set all crse cells to zero first
      for (int k = idlo_e[2]; k <= idhi_e[2]; ++k) {
         for (int j = idlo_e[1]; j <= idhi_e[1]; ++j) {
            for (int i = idlo_e[0]; i <= idhi_e[0]; ++i) {
               if (mask(i,j,k) == crsecell) data(i,j,k,n) = 0.0_rt;
            }
         }
      }

      // Corners
      // xlo, ylo, zlo
      {
         int i = idlo_e[0];
         int j = idlo_e[1];
         int k = idlo_e[2];
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
         int i = idlo_e[0];
         int j = idlo_e[1];
         int k = idhi_e[2];
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
         int i = idlo_e[0];
         int j = idhi_e[1];
         int k = idlo_e[2];
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
         int i = idlo_e[0];
         int j = idhi_e[1];
         int k = idhi_e[2];
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
         int i = idhi_e[0];
         int j = idlo_e[1];
         int k = idlo_e[2];
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
                               / ( mask(i-1,j+1,k) + mask(i+1,j,k+1) + mask(i,j+1,k+1) );
            } else {
               data(i,j,k,n) = data(i-1,j+1,k+1,n);
            }
         }
      }
      // xhi, ylo, zhi
      {
         int i = idhi_e[0];
         int j = idlo_e[1];
         int k = idhi_e[2];
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
         int i = idhi_e[0];
         int j = idhi_e[1];
         int k = idlo_e[2];
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
         int i = idhi_e[0];
         int j = idhi_e[1];
         int k = idhi_e[2];
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
         int i = idlo_e[0];
         int j = idlo_e[1];
         for (int k = idlo[2]; k <= idhi[2]; ++k) {
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
         int i = idlo_e[0];
         int j = idhi_e[1];
         for (int k = idlo[2]; k <= idhi[2]; ++k) {
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
         int i = idhi_e[0];
         int j = idlo_e[1];
         for (int k = idlo[2]; k <= idhi[2]; ++k) {
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
         int i = idhi_e[0];
         int j = idhi_e[1];
         for (int k = idlo[2]; k <= idhi[2]; ++k) {
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
         int i = idlo_e[0];
         int k = idlo_e[2];
         for (int j = idlo[1]; j <= idhi[1]; ++j) {
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
         int i = idlo_e[0];
         int k = idhi_e[2];
         for (int j = idlo[1]; j <= idhi[1]; ++j) {
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
         int i = idhi_e[0];
         int k = idlo_e[2];
         for (int j = idlo[1]; j <= idhi[1]; ++j) {
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
         int i = idhi_e[0];
         int k = idhi_e[2];
         for (int j = idlo[1]; j <= idhi[1]; ++j) {
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
         int j = idlo_e[1];
         int k = idlo_e[2];
         for (int i = idlo[0]; i <= idhi[0]; ++i) {
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
         int j = idlo_e[1];
         int k = idhi_e[2];
         for (int i = idlo[0]; i <= idhi[0]; ++i) {
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
         int j = idhi_e[1];
         int k = idlo_e[2];
         for (int i = idlo[0]; i <= idhi[0]; ++i) {
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
         int j = idhi_e[1];
         int k = idhi_e[2];
         for (int i = idlo[0]; i <= idhi[0]; ++i) {
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
         int i = idlo_e[0];
         for (int k = idlo[2]; k <= idhi[2]; ++k) {
            for (int j = idlo[1]; j <= idhi[1]; ++j) {
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
         int i = idhi_e[0];
         for (int k = idlo[2]; k <= idhi[2]; ++k) {
            for (int j = idlo[1]; j <= idhi[1]; ++j) {
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
         int j = idlo_e[1];
         for (int k = idlo[2]; k <= idhi[2]; ++k) {
            for (int i = idlo[0]; i <= idhi[0]; ++i) {
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
         int j = idhi_e[1];
         for (int k = idlo[2]; k <= idhi[2]; ++k) {
            for (int i = idlo[0]; i <= idhi[0]; ++i) {
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
         int k = idlo_e[2];
         for (int j = idlo[1]; j <= idhi[1]; ++j) {
            for (int i = idlo[0]; i <= idhi[0]; ++i) {
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
         int k = idhi_e[2];
         for (int j = idlo[1]; j <= idhi[1]; ++j) {
            for (int i = idlo[0]; i <= idhi[0]; ++i) {
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

}
#endif
