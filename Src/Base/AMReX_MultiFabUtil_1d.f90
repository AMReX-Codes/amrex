
subroutine amrex_fort_avg_nd_to_cc (lo, hi, ncomp, &
     cc, ccl1, cch1, &
     nd, ndl1, ndh1) bind(c)
  use amrex_fort_module, only : amrex_real
  implicit none
  integer          :: lo(1),hi(1), ncomp
  integer          :: ccl1, cch1
  integer          :: ndl1, ndh1
  real(amrex_real) :: cc(ccl1:cch1,ncomp)
  real(amrex_real) :: nd(ndl1:ndh1,ncomp)
  integer :: i, n
  do n = 1, ncomp
     do i=lo(1),hi(1)
        cc(i,n) = 0.5_amrex_real * ( nd(i,n) + nd(i+1,n) )
     end do
  end do
end subroutine amrex_fort_avg_nd_to_cc

! ***************************************************************************************
! subroutine bl_avg_eg_to_cc 
! ***************************************************************************************

subroutine bl_avg_eg_to_cc (lo, hi, &
     cc, ccl1, cch1, &
     Ex, Exl1, Exh1 )
  use amrex_fort_module, only : amrex_real
  implicit none
  integer          :: lo(1),hi(1)
  integer          :: ccl1, cch1
  integer          :: Exl1, Exh1
  real(amrex_real) :: cc(ccl1:cch1)
  real(amrex_real) :: Ex(Exl1:Exh1)

  cc(lo(1):hi(1)) = Ex(lo(1):hi(1))
end subroutine bl_avg_eg_to_cc


! ***************************************************************************************
! subroutine bl_avg_fc_to_cc
! ***************************************************************************************

subroutine bl_avg_fc_to_cc (lo, hi, &
     cc, ccl1, cch1, &
     fx, fxl1, fxh1, &
     dx, problo, coord_type)

  use amrex_fort_module, only : amrex_real
  implicit none
  integer          :: lo(1),hi(1), coord_type
  integer          :: ccl1, cch1
  integer          :: fxl1, fxh1
  real(amrex_real) :: cc(ccl1:cch1)
  real(amrex_real) :: fx(fxl1:fxh1)
  real(amrex_real) :: dx(1), problo(1)

  ! Local variables
  integer          :: i
  real(amrex_real) :: rlo,rhi,rcen

  if (coord_type .eq. 0) then ! Cartesian

     do i=lo(1),hi(1)
        cc(i) = 0.5d0 * ( fx(i) + fx(i+1) )
     enddo
  
  else if (coord_type .eq. 1) then ! Cylindrical

     do i=lo(1),hi(1)
        rlo = problo(1) + (dble(i)  )*dx(1)
        rhi = problo(1) + (dble(i+1))*dx(1)
        rcen = 0.5d0 * (rlo + rhi)
        cc(i) = 0.5d0 * ( rlo*fx(i) + rhi*fx(i+1) ) / rcen
     enddo
  
  else  ! Spherical

     do i=lo(1),hi(1)
        rlo = problo(1) + (dble(i)  )*dx(1)
        rhi = problo(1) + (dble(i+1))*dx(1)
        rcen = 0.5d0 * (rlo + rhi)
        cc(i) = 0.5d0 * ( rlo**2 * fx(i) + rhi**2 * fx(i+1) ) / rcen**2
     enddo

  end if

end subroutine bl_avg_fc_to_cc

! ***************************************************************************************
! subroutine bl_avg_cc_to_fc
! ***************************************************************************************

subroutine bl_avg_cc_to_fc (lo, hi, &
     fx, fxl1, fxh1, &
     cc, ccl1, cch1, &
     dx, problo, coord_type)

  use amrex_fort_module, only : amrex_real
  implicit none
  integer          :: lo(1),hi(1), coord_type
  integer          :: ccl1, cch1
  integer          :: fxl1, fxh1
  real(amrex_real) :: cc(ccl1:cch1)
  real(amrex_real) :: fx(fxl1:fxh1)
  real(amrex_real) :: dx(1), problo(1)

  ! Local variables
  integer          :: i
  real(amrex_real) :: rlo,rhi,rcen

  if (coord_type .eq. 0) then ! Cartesian

     do i = lo(1), hi(1)
        fx(i) = 0.5d0 * (cc(i-1) + cc(i))
     end do

  else if (coord_type .eq. 1) then ! Cylindrical

     do i = lo(1), hi(1)
        rlo = problo(1) + (dble(i-1)+0.5d0)*dx(1)
        rhi = problo(1) + (dble(i  )+0.5d0)*dx(1)
        rcen = 0.5d0 * (rlo + rhi)
        fx(i) = 0.5d0 * (rlo*cc(i-1) + rhi*cc(i)) / rcen
     end do

  else  ! Spherical

     do i = lo(1), hi(1)
        rlo = problo(1) + (dble(i-1)+0.5d0)*dx(1)
        rhi = problo(1) + (dble(i  )+0.5d0)*dx(1)
        rcen = 0.5d0 * (rlo + rhi)
        fx(i) = 0.5d0 * (rlo**2 * cc(i-1) + rhi**2 * cc(i)) / rcen**2
     end do

  end if

end subroutine bl_avg_cc_to_fc

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
