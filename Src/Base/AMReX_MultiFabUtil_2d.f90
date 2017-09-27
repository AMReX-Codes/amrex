subroutine amrex_fort_avg_nd_to_cc (lo, hi, ncomp, &
     cc, ccl1, ccl2, cch1, cch2, &
     nd, ndl1, ndl2, ndh1, ndh2 ) bind(c)

  use amrex_fort_module, only : amrex_real
  implicit none
  integer          :: lo(2),hi(2), ncomp
  integer          :: ccl1, ccl2, cch1, cch2
  integer          :: ndl1, ndl2, ndh1, ndh2
  real(amrex_real) :: cc(ccl1:cch1, ccl2:cch2, ncomp)
  real(amrex_real) :: nd(ndl1:ndh1, ndl2:ndh2, ncomp)

  ! Local variables
  integer          :: i,j,n
  
  do n = 1, ncomp
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)
           cc(i,j,n) = 0.25_amrex_real * ( nd(i,j,n) + nd(i+1,j,n) + nd(i,j+1,n) + nd(i+1,j+1,n) )
        end do
     end do
  end do

end subroutine amrex_fort_avg_nd_to_cc

! ***************************************************************************************
! subroutine bl_avg_eg_to_cc 
! ***************************************************************************************

subroutine bl_avg_eg_to_cc (lo, hi, &
     cc, ccl1, ccl2, cch1, cch2, &
     Ex, Exl1, Exl2, Exh1, Exh2, &
     Ey, Eyl1, Eyl2, Eyh1, Eyh2 )

  use amrex_fort_module, only : amrex_real
  implicit none
  integer          :: lo(2),hi(2)
  integer          :: ccl1, ccl2, cch1, cch2
  integer          :: Exl1, Exl2, Exh1, Exh2
  integer          :: Eyl1, Eyl2, Eyh1, Eyh2
  real(amrex_real) :: cc(ccl1:cch1, ccl2:cch2, 2)
  real(amrex_real) :: Ex(Exl1:Exh1, Exl2:Exh2)
  real(amrex_real) :: Ey(Eyl1:Eyh1, Eyl2:Eyh2)

  ! Local variables
  integer          :: i,j
  
  do j=lo(2),hi(2)
     do i=lo(1),hi(1)
        cc(i,j,1) = 0.5d0 * ( Ex(i,j) + Ex(i,j+1))
        cc(i,j,2) = 0.5d0 * ( Ey(i,j) + Ey(i+1,j))
     enddo
  enddo

end subroutine bl_avg_eg_to_cc

! ***************************************************************************************
! subroutine bl_avg_fc_to_cc
! ***************************************************************************************

subroutine bl_avg_fc_to_cc (lo, hi, &
     cc, ccl1, ccl2, cch1, cch2, &
     fx, fxl1, fxl2, fxh1, fxh2, &
     fy, fyl1, fyl2, fyh1, fyh2, &
     dx, problo, coord_type)

  use amrex_fort_module, only : amrex_real
  implicit none
  integer          :: lo(2),hi(2), coord_type
  integer          :: ccl1, ccl2, cch1, cch2
  integer          :: fxl1, fxl2, fxh1, fxh2
  integer          :: fyl1, fyl2, fyh1, fyh2
  real(amrex_real) :: cc(ccl1:cch1, ccl2:cch2, 2)
  real(amrex_real) :: fx(fxl1:fxh1, fxl2:fxh2)
  real(amrex_real) :: fy(fyl1:fyh1, fyl2:fyh2)
  real(amrex_real) :: dx(2), problo(2)

  ! Local variables
  integer          :: i,j
  
  do j=lo(2),hi(2)
     do i=lo(1),hi(1)
        cc(i,j,1) = 0.5d0 * ( fx(i,j) + fx(i+1,j) )
        cc(i,j,2) = 0.5d0 * ( fy(i,j) + fy(i,j+1) )
     enddo
  enddo

end subroutine bl_avg_fc_to_cc

! ***************************************************************************************
! subroutine bl_avg_cc_to_fc
! ***************************************************************************************

subroutine bl_avg_cc_to_fc (xlo, xhi, ylo, yhi, &
     fx, fxl1, fxl2, fxh1, fxh2, &
     fy, fyl1, fyl2, fyh1, fyh2, &
     cc, ccl1, ccl2, cch1, cch2, &
     dx, problo, coord_type)

  use amrex_fort_module, only : amrex_real
  implicit none
  integer          :: xlo(2),xhi(2), ylo(2),yhi(2), coord_type
  integer          :: fxl1, fxl2, fxh1, fxh2
  integer          :: fyl1, fyl2, fyh1, fyh2
  integer          :: ccl1, ccl2, cch1, cch2
  real(amrex_real) :: cc(ccl1:cch1, ccl2:cch2)
  real(amrex_real) :: fx(fxl1:fxh1, fxl2:fxh2)
  real(amrex_real) :: fy(fyl1:fyh1, fyl2:fyh2)
  real(amrex_real) :: dx(2), problo(2)

  ! Local variables
  integer          :: i,j
  
  do j=xlo(2),xhi(2)
     do i=xlo(1),xhi(1)
        fx(i,j) = 0.5d0 * (cc(i-1,j) + cc(i,j))
     enddo
  enddo

  do j=ylo(2),yhi(2)
     do i=ylo(1),yhi(1)
        fy(i,j) = 0.5d0 * (cc(i,j-1) + cc(i,j))
     enddo
  enddo

end subroutine bl_avg_cc_to_fc

! ***************************************************************************************
! subroutine bl_avgdown_faces
! ***************************************************************************************

subroutine bl_avgdown_faces (lo, hi, &
     f, f_l1, f_l2, f_h1, f_h2, &
     c, c_l1, c_l2, c_h1, c_h2, &
     ratio, idir,nc)

  use amrex_fort_module, only : amrex_real
  implicit none
  integer          :: lo(2),hi(2)
  integer          :: f_l1, f_l2, f_h1, f_h2
  integer          :: c_l1, c_l2, c_h1, c_h2
  integer          :: ratio(3), idir,nc
  real(amrex_real) :: f(f_l1:f_h1, f_l2:f_h2, nc)
  real(amrex_real) :: c(c_l1:c_h1, c_l2:c_h2, nc)

  ! Local variables
  integer i,j,n,facx,facy,iref,jref, ii, jj
  real(amrex_real) :: facInv

  facx = ratio(1)
  facy = ratio(2)

  if (idir .eq. 0) then

     facInv = 1.d0 / facy

     do n = 1, nc
        do j     = lo(2), hi(2)
           jj    = j * facy
           do i  = lo(1), hi(1)
              ii = i * facx
              c(i,j,n) = 0.d0
              do jref = 0, facy-1
                 c(i,j,n) = c(i,j,n) + f(ii,jj+jref,n)
              end do
              c(i,j,n) = c(i,j,n) * facInv
           end do
        end do
     end do

  else 

     facInv = 1.d0 / facx

     do n = 1, nc
        do j     = lo(2), hi(2)
           jj    = j * facy
           do i  = lo(1), hi(1)
              ii = i * facx
              c(i,j,n) = 0.d0
              do iref = 0,facx-1
                 c(i,j,n) = c(i,j,n) + f(ii+iref,jj,n)
              end do
              c(i,j,n) = c(i,j,n) *facInv
           end do
        end do
     end do

  end if

end subroutine bl_avgdown_faces

! ***************************************************************************************
! subroutine bl_avgdown_edges
! ***************************************************************************************

subroutine bl_avgdown_edges (lo, hi, &
     f, f_l1, f_l2, f_h1, f_h2, &
     c, c_l1, c_l2, c_h1, c_h2, &
     ratio, idir,nc)

  use amrex_fort_module, only : amrex_real
  implicit none
  integer          :: lo(2),hi(2)
  integer          :: f_l1, f_l2, f_h1, f_h2
  integer          :: c_l1, c_l2, c_h1, c_h2
  integer          :: ratio(3), idir,nc
  real(amrex_real) :: f(f_l1:f_h1, f_l2:f_h2, nc)
  real(amrex_real) :: c(c_l1:c_h1, c_l2:c_h2, nc)

  ! Local variables
  integer i,j,n,facx,facy,iref,jref, ii, jj
  real(amrex_real) :: facInv

  facx = ratio(1)
  facy = ratio(2)

  if (idir .eq. 0) then

     facInv = 1.d0 / facx

     do n = 1, nc
        do j     = lo(2), hi(2)
           jj    = j * facy
           do i  = lo(1), hi(1)
              ii = i * facx
              c(i,j,n) = 0.d0
              do iref = 0, facx-1
                 c(i,j,n) = c(i,j,n) + f(ii+iref,jj,n)
              end do
              c(i,j,n) = c(i,j,n) * facInv
           end do
        end do
     end do

  else 

     facInv = 1.d0 / facy

     do n = 1, nc
        do j     = lo(2), hi(2)
           jj    = j * facy
           do i  = lo(1), hi(1)
              ii = i * facx
              c(i,j,n) = 0.d0
              do jref = 0,facy-1
                 c(i,j,n) = c(i,j,n) + f(ii,jj+jref,n)
              end do
              c(i,j,n) = c(i,j,n) *facInv
           end do
        end do
     end do

  end if

end subroutine bl_avgdown_edges

! ***************************************************************************************
! subroutine bl_avgdown - THIS VERISON DOES NOT DO VOLUME WEIGHTING
! ***************************************************************************************

subroutine bl_avgdown (lo,hi,&
     fine,f_l1,f_l2,f_h1,f_h2, &
     crse,c_l1,c_l2,c_h1,c_h2, &
     lrat,ncomp)
  
  use amrex_fort_module, only : amrex_real
  implicit none
  
  integer f_l1,f_l2,f_h1,f_h2
  integer c_l1,c_l2,c_h1,c_h2
  integer lo(2), hi(2)
  integer lrat(2), ncomp
  real(amrex_real) fine(f_l1:f_h1,f_l2:f_h2,ncomp)
  real(amrex_real) crse(c_l1:c_h1,c_l2:c_h2,ncomp)

  integer :: i, j, ii, jj, n, iref, jref
  real(amrex_real) :: volfrac

  volfrac = 1.d0 / dble(lrat(1)*lrat(2))

  do n = 1, ncomp
     do j     = lo(2), hi(2)
        jj    = j * lrat(2)
        do i  = lo(1), hi(1)
           ii = i * lrat(1)
           crse(i,j,n) = 0.d0
           do    jref = 0, lrat(2)-1
              do iref = 0, lrat(1)-1
                 crse(i,j,n) = crse(i,j,n) + fine(ii+iref,jj+jref,n)
              end do
           end do
           crse(i,j,n) = volfrac * crse(i,j,n)
        end do
     end do
  end do

end subroutine bl_avgdown

! ***************************************************************************************
! subroutine bl_avgdown_nodes
! ***************************************************************************************

subroutine bl_avgdown_nodes (lo,hi,&
     fine,f_l1,f_l2,f_h1,f_h2, &
     crse,c_l1,c_l2,c_h1,c_h2, &
     lrat,ncomp)
  
  use amrex_fort_module, only : amrex_real
  implicit none
  
  integer f_l1,f_l2,f_h1,f_h2
  integer c_l1,c_l2,c_h1,c_h2
  integer lo(2), hi(2)
  integer lrat(2), ncomp
  real(amrex_real) fine(f_l1:f_h1,f_l2:f_h2,ncomp)
  real(amrex_real) crse(c_l1:c_h1,c_l2:c_h2,ncomp)

  integer :: i, j, ii, jj, n

  do n = 1, ncomp
     do j     = lo(2), hi(2)
        jj    = j * lrat(2)
        do i  = lo(1), hi(1)
           ii = i * lrat(1)
           crse(i,j,n) = fine(ii, jj, n)
        end do
     end do
  end do

end subroutine bl_avgdown_nodes

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


