! ***************************************************************************************
! subroutine bl_avgdown_with_vol
! ***************************************************************************************

subroutine bl_avgdown_with_vol (lo,hi,&
     fine,f_l1,f_l2,f_h1,f_h2, &
     crse,c_l1,c_l2,c_h1,c_h2, &
     fv,fv_l1,fv_l2,fv_h1,fv_h2, &
     lrat,ncomp)

  use amrex_fort_module, only : amrex_real
  implicit none
  
  integer f_l1,f_l2,f_h1,f_h2
  integer c_l1,c_l2,c_h1,c_h2
  integer fv_l1,fv_l2,fv_h1,fv_h2
  integer lo(2), hi(2)
  integer lrat(2), ncomp
  real(amrex_real) fine(f_l1:f_h1,f_l2:f_h2,ncomp)
  real(amrex_real) crse(c_l1:c_h1,c_l2:c_h2,ncomp)
  real(amrex_real) fv(fv_l1:fv_h1,fv_l2:fv_h2)

  integer :: i, j, ii, jj, n, iref, jref
  real(amrex_real) :: cv

  do n = 1, ncomp
     do j     = lo(2), hi(2)
        jj    = j * lrat(2)
        do i  = lo(1), hi(1)
           ii = i * lrat(1)
           crse(i,j,n) = 0.d0
           cv          = 0.d0
           do    jref = 0, lrat(2)-1
              do iref = 0, lrat(1)-1
                 cv          = cv          +                         fv(ii+iref,jj+jref)
                 crse(i,j,n) = crse(i,j,n) + fine(ii+iref,jj+jref,n)*fv(ii+iref,jj+jref)
              end do
           end do
           crse(i,j,n) = crse(i,j,n) / cv
        end do
     end do
  end do

end subroutine bl_avgdown_with_vol

