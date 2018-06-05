module amrex_mg_module

  use amrex_fort_module
  use amrex_constants_module

  implicit none

contains

    subroutine amrex_mg_average ( &
           c, c_l1,c_h1, &
           f, f_l1,f_h1, &
           lo, hi, nc) bind(c,name='amrex_mg_average')

      integer nc
      integer f_l1,f_h1
      integer c_l1,c_h1
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      real(amrex_real) f(f_l1:f_h1,nc)
      real(amrex_real) c(c_l1:c_h1,nc)

      integer i,n

      do n = 1, nc
         do i = lo(1), hi(1)
            c(i,n) =  half * ( f(2*i+1,n) + f(2*i,n) )
         end do
      end do

    end subroutine amrex_mg_average

    subroutine amrex_mg_interp ( &
           f, f_l1,f_h1, &
           c, c_l1,c_h1, &
           lo, hi, nc) bind(c,name='amrex_mg_interp')

      integer nc
      integer f_l1,f_h1
      integer c_l1,c_h1
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      real(amrex_real) f(f_l1:f_h1,nc)
      real(amrex_real) c(c_l1:c_h1,nc)

      integer i,n

      do n = 1, nc
         do i = lo(1), hi(1)
            f(2*i+1,n) = c(i,n) + f(2*i+1,n)
            f(2*i  ,n) = c(i,n) + f(2*i  ,n)
         end do
      end do

    end subroutine amrex_mg_interp

end module amrex_mg_module
