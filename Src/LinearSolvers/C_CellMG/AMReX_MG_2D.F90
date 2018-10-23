module amrex_mg_module

  use amrex_fort_module
  use amrex_constants_module

  implicit none

contains

    subroutine amrex_mg_average ( &
           c, c_l1,c_l2,c_h1,c_h2, &
           f, f_l1,f_l2,f_h1,f_h2, &
           lo, hi, nc) bind(c,name='amrex_mg_average')

      implicit none

      integer nc
      integer f_l1,f_l2,f_h1,f_h2
      integer c_l1,c_l2,c_h1,c_h2
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      real(amrex_real) f(f_l1:f_h1,f_l2:f_h2,nc)
      real(amrex_real) c(c_l1:c_h1,c_l2:c_h2,nc)

      integer i
      integer j
      integer n
      real(amrex_real) denom
      parameter(denom=fourth)

      do n = 1, nc
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               c(i,j,n) =  ( &
                    f(2*i+1,2*j+1,n) + f(2*i  ,2*j+1,n) &
                    + f(2*i+1,2*j,n ) + f(2*i  ,2*j ,n))*denom
            end do
         end do
      end do

    end subroutine amrex_mg_average

    subroutine amrex_mg_interp ( &
           f, f_l1,f_l2,f_h1,f_h2, &
           c, c_l1,c_l2,c_h1,c_h2, &
           lo, hi, nc) bind(c,name='amrex_mg_interp')

      implicit none

      integer nc
      integer f_l1,f_l2,f_h1,f_h2
      integer c_l1,c_l2,c_h1,c_h2
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      real(amrex_real) f(f_l1:f_h1,f_l2:f_h2,nc)
      real(amrex_real) c(c_l1:c_h1,c_l2:c_h2,nc)

      integer i, j, n, twoi, twoj, twoip1, twojp1

!     MultiGrid::relax(...) does only V-cycles (not F-cycles), and for V-cycles, 
!     piecewise-constant interpolation performs better than linear interpolation,
!     as measured both by run-time and number of V-cycles for convergence.

      do n = 1, nc
         do j = lo(2),hi(2)
            twoj   = 2*j
            twojp1 = twoj+1

            do i = lo(1),hi(1)

               twoi   = 2*i
               twoip1 = twoi+1

               f(twoi,   twoj  ,n) = f(twoi,   twoj  ,n) + c(i,j,n)
               f(twoip1, twoj  ,n) = f(twoip1, twoj  ,n) + c(i,j,n)
               f(twoi,   twojp1,n) = f(twoi,   twojp1,n) + c(i,j,n)
               f(twoip1, twojp1,n) = f(twoip1, twojp1,n) + c(i,j,n)

            end do
         end do
      end do

    end subroutine amrex_mg_interp

end module amrex_mg_module
