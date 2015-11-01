! ***************************************************************************************
! subroutine bl_avg_fc_to_cc
! ***************************************************************************************

subroutine bl_avg_fc_to_cc (lo, hi, &
     cc, ccl1, ccl2, cch1, cch2, &
     fx, fxl1, fxl2, fxh1, fxh2, &
     fy, fyl1, fyl2, fyh1, fyh2, &
     dx, problo, coord_type)

  implicit none
  integer          :: lo(2),hi(2), coord_type
  integer          :: ccl1, ccl2, cch1, cch2
  integer          :: fxl1, fxl2, fxh1, fxh2
  integer          :: fyl1, fyl2, fyh1, fyh2
  double precision :: cc(ccl1:cch1, ccl2:cch2, 2)
  double precision :: fx(fxl1:fxh1, fxl2:fxh2)
  double precision :: fy(fyl1:fyh1, fyl2:fyh2)
  double precision :: dx(2), problo(2)

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

  implicit none
  integer          :: xlo(2),xhi(2), ylo(2),yhi(2), coord_type
  integer          :: fxl1, fxl2, fxh1, fxh2
  integer          :: fyl1, fyl2, fyh1, fyh2
  integer          :: ccl1, ccl2, cch1, cch2
  double precision :: cc(ccl1:cch1, ccl2:cch2)
  double precision :: fx(fxl1:fxh1, fxl2:fxh2)
  double precision :: fy(fyl1:fyh1, fyl2:fyh2)
  double precision :: dx(2), problo(2)

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
     ratio, idir)

  implicit none
  integer          :: lo(2),hi(2)
  integer          :: f_l1, f_l2, f_h1, f_h2
  integer          :: c_l1, c_l2, c_h1, c_h2
  integer          :: ratio(3), idir
  double precision :: f(f_l1:f_h1, f_l2:f_h2)
  double precision :: c(c_l1:c_h1, c_l2:c_h2)

  ! Local variables
  integer i,j,n,facx,facy

  facx = ratio(1)
  facy = ratio(2)

  if (idir .eq. 0) then
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               c(i,j) = 0.d0
               do n = 0,facy-1
                  c(i,j) = c(i,j) + f(facx*i,facy*j+n)
               end do
               c(i,j) = c(i,j) / facy
            end do
         end do
  else 
         do i = lo(1), hi(1)
            do j = lo(2), hi(2)
               c(i,j) = 0.d0
               do n = 0,facx-1
                  c(i,j) = c(i,j) + f(facx*i+n,facy*j)
               end do
               c(i,j) = c(i,j) / facx
            end do
         end do
  end if

end subroutine bl_avgdown_faces

! ***************************************************************************************
! subroutine bl_avgdown - THIS VERISON DOES NOT DO VOLUME WEIGHTING
! ***************************************************************************************

      subroutine bl_avgdown (lo,hi,&
                             fine,f_l1,f_l2,f_h1,f_h2, &
                             crse,c_l1,c_l2,c_h1,c_h2, &
                             lrat,ncomp)

      use bl_constants_module

      implicit none

      integer f_l1,f_l2,f_h1,f_h2
      integer c_l1,c_l2,c_h1,c_h2
      integer lo(2), hi(2)
      integer lrat(2), ncomp
      double precision fine(f_l1:f_h1,f_l2:f_h2,ncomp)
      double precision crse(c_l1:c_h1,c_l2:c_h2,ncomp)

      integer clo(2), chi(2)
      integer i, j, ic, jc, ioff, joff
      double precision volfrac

      clo(1:2) = lo(1:2) / lrat(1:2)
      chi(1:2) = hi(1:2) / lrat(1:2)

      !
      ! ::::: set coarse grid to zero on overlap
      !
      do jc = clo(2), chi(2)
         do ic = clo(1), chi(1)
            crse(ic,jc,:) = ZERO
         enddo
      enddo
      !
      ! ::::: sum fine data
      !
      do joff = 0, lrat(2)-1
        do jc = clo(2), chi(2)
          j = jc*lrat(2) + joff
          do ioff = 0, lrat(1)-1
            do ic = clo(1), chi(1)
              i = ic*lrat(1) + ioff
              crse(ic,jc,1:ncomp) = crse(ic,jc,1:ncomp) + fine(i,j,1:ncomp)
            enddo
          enddo
        enddo
      enddo

      volfrac = ONE/dble(lrat(1)*lrat(2))
      do jc = clo(2), chi(2)
         do ic = clo(1), chi(1)
            crse(ic,jc,1:ncomp) = volfrac*crse(ic,jc,1:ncomp)
         enddo
      enddo

      end subroutine bl_avgdown

! ***************************************************************************************
! subroutine bl_avgdown_with_vol
! ***************************************************************************************

      subroutine bl_avgdown_with_vol (lo,hi,&
                             fine,f_l1,f_l2,f_h1,f_h2, &
                             crse,c_l1,c_l2,c_h1,c_h2, &
                             fv,fv_l1,fv_l2,fv_h1,fv_h2, &
                             cv,cv_l1,cv_l2,cv_h1,cv_h2, &
                             lrat,ncomp)

      use bl_constants_module

      implicit none

      integer f_l1,f_l2,f_h1,f_h2
      integer c_l1,c_l2,c_h1,c_h2
      integer fv_l1,fv_l2,fv_h1,fv_h2
      integer cv_l1,cv_l2,cv_h1,cv_h2
      integer lo(2), hi(2)
      integer lrat(2), ncomp
      double precision fine(f_l1:f_h1,f_l2:f_h2,ncomp)
      double precision crse(c_l1:c_h1,c_l2:c_h2,ncomp)
      double precision fv(fv_l1:fv_h1,fv_l2:fv_h2)
      double precision cv(cv_l1:cv_h1,cv_l2:cv_h2)

      integer clo(2), chi(2)
      integer i, j, ic, jc, ioff, joff

      clo(1:2) = lo(1:2) / lrat(1:2)
      chi(1:2) = hi(1:2) / lrat(1:2)

      !
      ! ::::: set coarse grid to zero on overlap
      !
      do jc = clo(2), chi(2)
         do ic = clo(1), chi(1)
            crse(ic,jc,:) = ZERO
         enddo
      enddo
      !
      ! ::::: sum fine data
      !
      do joff = 0, lrat(2)-1
        do jc = clo(2), chi(2)
          j = jc*lrat(2) + joff
          do ioff = 0, lrat(1)-1
            do ic = clo(1), chi(1)
              i = ic*lrat(1) + ioff
              crse(ic,jc,1:ncomp) = crse(ic,jc,1:ncomp) + fv(i,j)*fine(i,j,1:ncomp)
            enddo
          enddo
        enddo
      enddo

      do jc = clo(2), chi(2)
         do ic = clo(1), chi(1)
            crse(ic,jc,1:ncomp) = crse(ic,jc,1:ncomp) / cv(ic,jc)
         enddo
      enddo

      end subroutine bl_avgdown_with_vol

