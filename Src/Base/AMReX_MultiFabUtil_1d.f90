
! ***************************************************************************************
! subroutine bl_avgdown_faces
! ***************************************************************************************

subroutine bl_avgdown_faces (lo, hi, &
     f, f_l1, f_h1, &
     c, c_l1, c_h1, &
     ratio,idir,nc)

  use amrex_fort_module, only : amrex_real
  implicit none
  integer          :: lo(1),hi(1)
  integer          :: f_l1, f_h1
  integer          :: c_l1, c_h1
  integer          :: ratio(1), idir, nc
  real(amrex_real) :: f(f_l1:f_h1,nc)
  real(amrex_real) :: c(c_l1:c_h1,nc)

  ! Local variables
  integer i,n,facx

  facx = ratio(1)

   ! lo(1)..hi(1) are edge base indices
  do n = 1, nc
     do i = lo(1), hi(1)
        c(i,n) = f(facx*i,n)
     end do
  end do

end subroutine bl_avgdown_faces

! ***************************************************************************************
! subroutine bl_avgdown_edges
! ***************************************************************************************

subroutine bl_avgdown_edges (lo, hi, &
     f, f_l1, f_h1, &
     c, c_l1, c_h1, &
     ratio,idir,nc)

  use amrex_fort_module, only : amrex_real
  implicit none
  integer          :: lo(1),hi(1)
  integer          :: f_l1, f_h1
  integer          :: c_l1, c_h1
  integer          :: ratio(1), idir, nc
  real(amrex_real) :: f(f_l1:f_h1,nc)
  real(amrex_real) :: c(c_l1:c_h1,nc)

  ! Local variables
  integer i,n,facx,iref
  real(amrex_real) :: facInv

  facx = ratio(1)

  facInv = 1.d0 / facx

  do n = 1, nc
     do i = lo(1), hi(1)
        c(i,n) = 0._amrex_real
        do iref = 0, facx-1
           c(i,n) = c(i,n) + f(facx*i+iref,n)
        end do
        c(i,n) = c(i,n) * facInv
     end do
  end do

end subroutine bl_avgdown_edges

! ***************************************************************************************
! subroutine bl_avgdown - THIS VERSION DOES NOT USE VOLUME WEIGHTING
! ***************************************************************************************

subroutine bl_avgdown (lo,hi,&
     fine,f_l1,f_h1, &
     crse,c_l1,c_h1, &
     lrat,ncomp)
  
  use amrex_fort_module, only : amrex_real
  implicit none
  
  integer f_l1,f_h1
  integer c_l1,c_h1
  integer lo(1), hi(1)
  integer lrat(1), ncomp
  real(amrex_real) fine(f_l1:f_h1,ncomp)
  real(amrex_real) crse(c_l1:c_h1,ncomp)
  
  integer :: i, ii, n, iref
  real(amrex_real) :: volfrac
  
  volfrac = 1.d0 / dble(lrat(1))
  
  do n = 1, ncomp
     do i = lo(1), hi(1)
        ii = i * lrat(1)
        crse(i,n) = 0.d0
        do iref = 0, lrat(1)-1
           crse(i,n) = crse(i,n) + fine(ii+iref,n)
        end do
        crse(i,n) = volfrac*crse(i,n)
     end do
  end do
  
end subroutine bl_avgdown

! ***************************************************************************************
! subroutine bl_avgdown_nodes
! ***************************************************************************************

subroutine bl_avgdown_nodes (lo,hi,&
     fine,f_l1,f_h1, &
     crse,c_l1,c_h1, &
     lrat,ncomp)
  
  use amrex_fort_module, only : amrex_real
  implicit none
  
  integer f_l1,f_h1
  integer c_l1,c_h1
  integer lo(1), hi(1)
  integer lrat(1), ncomp
  real(amrex_real) fine(f_l1:f_h1,ncomp)
  real(amrex_real) crse(c_l1:c_h1,ncomp)
  
  integer :: i, ii, n
  
  do n = 1, ncomp
     do i = lo(1), hi(1)
        ii = i * lrat(1)
        crse(i,n) = fine(ii, n)
     end do
  end do
  
end subroutine bl_avgdown_nodes

! ***************************************************************************************
! subroutine bl_avgdown_with_vol
! ***************************************************************************************

subroutine bl_avgdown_with_vol (lo,hi,&
     fine,f_l1,f_h1, &
     crse,c_l1,c_h1, &
     fv,fv_l1,fv_h1, &
     lrat,ncomp)

  use amrex_fort_module, only : amrex_real
  implicit none

  integer f_l1,f_h1
  integer c_l1,c_h1
  integer fv_l1,fv_h1
  integer lo(1), hi(1)
  integer lrat(1), ncomp
  real(amrex_real) fine(f_l1:f_h1,ncomp)
  real(amrex_real) crse(c_l1:c_h1,ncomp)
  real(amrex_real) fv(fv_l1:fv_h1)

  integer :: i, ii, n, iref
  real(amrex_real) :: cv

  do n = 1, ncomp
     do i = lo(1), hi(1)
        ii = i * lrat(1)
        crse(i,n) = 0.d0
        cv        = 0.d0
        do iref = 0, lrat(1)-1
           cv        = cv        +                 fv(ii+iref)
           crse(i,n) = crse(i,n) + fine(ii+iref,n)*fv(ii+iref)
        end do
        crse(i,n) = crse(i,n) / cv
     end do
  end do

end subroutine bl_avgdown_with_vol


subroutine amrex_compute_divergence (lo, hi, divu, dlo, dhi, u, ulo, uhi, dxinv) bind(c)
  use amrex_fort_module, only : amrex_real
  implicit none
  integer, dimension(1), intent(in) :: lo, hi, dlo, dhi, ulo, uhi
  real(amrex_real), intent(inout) :: divu(dlo(1):dhi(1))
  real(amrex_real), intent(in   ) ::    u(ulo(1):uhi(1))
  real(amrex_real), intent(in) :: dxinv(1)
  integer :: i
  do i = lo(1), hi(1)
     divu(i) = dxinv(1) * (u(i+1)-u(i))
  end do
end subroutine amrex_compute_divergence
