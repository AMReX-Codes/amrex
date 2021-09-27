
module amrex_interp_module

  use amrex_fort_module
  use amrex_constants_module
  use amrex_bc_types_module
  use amrex_error_module

  implicit none

contains

#define IX_PROJ(A,B) (A+B*iabs(A))/B-iabs(A)

! :::
! ::: --------------------------------------------------------------
! ::: cqinterp:  cell-centered quadratic interpolation
! :::
! ::: NOTE: it is assumed that the coarse grid array is
! ::: large enough to define interpolated values
! ::: in the region fblo:fbhi on the fine grid
! :::
! :::

     subroutine AMREX_CQINTERP (fine, fine_l1,fine_l2,fine_h1,fine_h2, &
                               fb_l1, fb_l2, fb_h1, fb_h2, &
                               nvar, lratiox, lratioy, crse, clo, chi, &
                               cb_l1, cb_l2, cb_h1, cb_h2, &
                               fslo, fshi, cslope, clen, fslope, fdat, &
                               flen, voff, bc, limslope, &
                               fvcx, fvcy, cvcx, cvcy, &
                               actual_comp,actual_state) bind(c,name='amrex_cqinterp')

      implicit none

      integer fine_l1,fine_l2,fine_h1,fine_h2
      integer fslo(2), fshi(2)
      integer fb_l1, fb_l2, fb_h1, fb_h2
      integer cb_l1, cb_l2, cb_h1, cb_h2
      integer clo, chi
      integer lratiox, lratioy, nvar, clen, flen, limslope
      integer bc(2,2,nvar)
      integer actual_comp,actual_state
      real(amrex_real) fine(fine_l1:fine_h1,fine_l2:fine_h2,nvar)
      real(amrex_real) crse(clo:chi, nvar)
      real(amrex_real) cslope(clo:chi, 5)
      real(amrex_real) fslope(flen, 5)
      real(amrex_real) fdat(flen)
      real(amrex_real) voff(flen)
      real(amrex_real) fvcx(fb_l1:fb_h1+1)
      real(amrex_real) fvcy(fb_l2:fb_h2+1)
      real(amrex_real) cvcx(cb_l1:cb_h1+1)
      real(amrex_real) cvcy(cb_l2:cb_h2+1)

#define bclo(i,n) bc(i,1,n)
#define bchi(i,n) bc(i,2,n)

      integer n, fn
      integer i, ic, ioff
      integer j, jc, joff
      integer ist, jst
      real(amrex_real) cen
      real(amrex_real) fcen, ccen
      real(amrex_real) diffxy,diffxx,diffyy
      real(amrex_real) yoff
      integer ncbx, ncby
      integer ncsx
      integer jslo
      integer icc, istart, iend
      logical xok, yok

      ncbx = cb_h1-cb_l1+1
      ncby = cb_h2-cb_l2+1
      xok = (ncbx .ge. 2)
      yok = (ncby .ge. 2)
      ncsx = ncbx+2
      ist = 1
      jst = ncsx
      jslo = cb_l2-1

      do i = fb_l1, fb_h1
         fn = i-fslo(1)+1
         ic = IX_PROJ(i,lratiox)
         fcen = half*(fvcx(i)+fvcx(i+1))
         ccen = half*(cvcx(ic)+cvcx(ic+1))
         voff(fn) = (fcen-ccen)/(cvcx(ic+1)-cvcx(ic))
      end do

      do n = 1, nvar
        do i = clo, chi
          crse(i,n) = merge(crse(i,n),zero,abs(crse(i,n)).gt.1.0d-50)
        end do
      end do

      do 290 n = 1, nvar

            do i = 1, clen
               cen = half*(crse(i+ist,n)-crse(i-ist,n))
               diffxy = fourth*(crse(i+ist+jst,n)+crse(i-ist-jst,n) &
                               -crse(i-ist+jst,n)-crse(i+ist-jst,n))
               diffxx = crse(i+ist,n)-two*crse(i,n)+crse(i-ist,n)
               cslope(i,1)=cen
               cslope(i,3)=diffxx
               cslope(i,5)=diffxy
            end do
            if (xok) then
               if (bclo(1,n) .eq. amrex_bc_ext_dir .or. bclo(1,n).eq.amrex_bc_hoextrap) then
                  do i = 1, clen, jst
                     cen  = -sixteen/fifteen*crse(i-ist,n) + half*crse(i,n) &
                          + two3rd*crse(i+ist,n) - tenth*crse(i+2*ist,n)
                     cslope(i,1)=cen
                     cslope(i,3)=zero
                     cslope(i,5)=zero
                  end do
               end if
               if (bchi(1,n) .eq. amrex_bc_ext_dir .or. bchi(1,n).eq.amrex_bc_hoextrap) then
                  do i = ncbx, clen, jst
                     cen = sixteen/fifteen*crse(i+ist,n) - half*crse(i,n) &
                          - two3rd*crse(i-ist,n) + tenth*crse(i-2*ist,n)
                     cslope(i,1)=cen
                     cslope(i,3)=zero
                     cslope(i,5)=zero
                  end do
               end if
            end if

            do i = 1, clen
               cen  = half*(crse(i+jst,n)-crse(i-jst,n))
               diffyy = crse(i+jst,n)-two*crse(i,n)+crse(i-jst,n)
               cslope(i,2)=cen
               cslope(i,4)=diffyy
            end do
            if (yok) then
               if (bclo(2,n) .eq. amrex_bc_ext_dir .or. bclo(2,n).eq.amrex_bc_hoextrap) then
                  do i = 1, ncbx
                     cen  = -sixteen/fifteen*crse(i-jst,n) + half*crse(i,n) &
                          + two3rd*crse(i+jst,n) - tenth*crse(i+2*jst,n)
                     cslope(i,2)=cen
                     cslope(i,4)=zero
                     cslope(i,5)=zero
                  end do
               end if
               if (bchi(2,n) .eq. amrex_bc_ext_dir .or. bchi(2,n).eq.amrex_bc_hoextrap) then
                  do i = clen-ncbx,clen
                     cen = sixteen/fifteen*crse(i+jst,n) - half*crse(i,n) &
                          - two3rd*crse(i-jst,n) + tenth*crse(i-2*jst,n)
                     cslope(i,2)=cen
                     cslope(i,4)=zero
                     cslope(i,5)=zero
                  end do
               end if
            end if

            do 360 jc = cb_l2, cb_h2
               do 370 ioff = 1, lratiox
                  icc = clo + ist + jst*(jc-jslo)
                  istart = ioff
                  iend = ioff + (ncbx-1)*lratiox
                  do 380 fn = istart, iend, lratiox
                     fslope(fn,1) = cslope(icc,1)
                     fslope(fn,2) = cslope(icc,2)
                     fslope(fn,3) = cslope(icc,3)
                     fslope(fn,4) = cslope(icc,4)
                     fslope(fn,5) = cslope(icc,5)
                     fdat(fn) = crse(icc,n)
                     icc = icc + ist
380               continue
370            continue

               do 390 joff = 0, lratioy-1
                  j = lratioy*jc + joff
                  if ((j.lt.fb_l2).or.(j.gt.fb_h2)) goto 390
                  fcen = half*(fvcy(j)+fvcy(j+1))
                  ccen = half*(cvcy(jc)+cvcy(jc+1))
                  yoff = (fcen-ccen)/(cvcy(jc+1)-cvcy(jc))

                  do 400 i = fb_l1, fb_h1
                     fn = i-fslo(1)+1
                     fine(i,j,n) = fdat(fn) &
                          + voff(fn)              *fslope(fn,1) &
                          + yoff                  *fslope(fn,2) &
                          + half*voff(fn)*voff(fn)*fslope(fn,3) &
                          + half*yoff    *yoff    *fslope(fn,4) &
                          + voff(fn)*yoff         *fslope(fn,5)
400               continue
390            continue
360         continue

290   continue

    end subroutine AMREX_CQINTERP

! :::
! ::: --------------------------------------------------------------
! ::: quartinterp: quartic conservative interpolation from coarse grid to
! ::: subregion of fine grid defined by (fblo,fbhi)
! :::
! ::: Inputs/Outputs
! ::: fine        <=>  (modify) fine grid array
! ::: flo,fhi      =>  (const)  index limits of fine grid
! ::: fblo,fbhi    =>  (const)  subregion of fine grid to get values
! ::: nvar         =>  (const)  number of variables in state vector
! ::: lratio[xy]   =>  (const)  refinement ratio between levels
! :::
! ::: crse         =>  (const)  coarse grid data
! ::: clo,chi      =>  (const)  index limits of crse grid
! ::: cblo,cbhi    =>  (const)  coarse grid region containing fblo,fbhi and widen by 2 or 4 cells
! :::
! ::: cb2lo,cb2hi  =>  (const)  coarse grid region containing fblo,fbhi
! ::: fb2lo,fb2hi  =>  (const)  fine version of cb2. It could be wider than fb
! :::
! ::: TEMPORARY ARRAYS
! ::: ftmp         =>  1-D temp array
! ::: ctmp         =>  2-D temp array
! ::: --------------------------------------------------------------
! :::
     subroutine AMREX_QUARTINTERP (fine, fine_l1,fine_l2,fine_h1,fine_h2, &
                                  fblo, fbhi, fb2lo, fb2hi, &
                                  crse, crse_l1,crse_l2,crse_h1,crse_h2, &
                                  cblo, cbhi, cb2lo, cb2hi, &
                                  nvar, &
                                  lratiox, lratioy, &
                                  ftmp, ctmp, &
                                  bc,actual_comp,actual_state) bind(c,name='amrex_quartinterp')

       implicit none

       integer fine_l1,fine_l2,fine_h1,fine_h2
       integer crse_l1,crse_l2,crse_h1,crse_h2
       integer fblo(2), fbhi(2), fb2lo(2), fb2hi(2)
       integer cblo(2), cbhi(2), cb2lo(2), cb2hi(2)
       integer nvar,lratiox,lratioy
       integer bc(2,2,nvar)
       integer actual_comp,actual_state
       real(amrex_real) fine(fine_l1:fine_h1,fine_l2:fine_h2,nvar)
       real(amrex_real) crse(crse_l1:crse_h1,crse_l2:crse_h2,nvar)
       real(amrex_real) ftmp(fb2lo(1):fb2hi(1))
       real(amrex_real) ctmp(cblo(1):cbhi(1),0:lratioy-1)

!      Local variables
       integer i,j,ii,jj,n,iry
       real(amrex_real) cL(-2:2)
!       real(amrex_real) cR(-2:2)
       data cL/ -0.01171875D0,  0.0859375D0, 0.5d0, -0.0859375D0, &
                 0.01171875D0 /
!       data cR/  0.01171875D0, -0.0859375D0, 0.5d0,  0.0859375D0, &
!                -0.01171875D0 /

       if (lratiox.eq.2 .and. lratioy.eq.2) then
          do n = 1, nvar
             do j = cb2lo(2), cb2hi(2)
                do i = cblo(1), cbhi(1)
                   ctmp(i,0) = 2.d0*(cL(-2)*crse(i,j-2,n) &
                        +            cL(-1)*crse(i,j-1,n) &
                        +            cL( 0)*crse(i,j  ,n) &
                        +            cL( 1)*crse(i,j+1,n) &
                        +            cL( 2)*crse(i,j+2,n))
                   ctmp(i,1) = 2.d0*crse(i,j,n) - ctmp(i,0)
!$$$                   ctmp(i,1) = 2.d0*(cR(-2)*crse(i,j-2,n) &
!$$$                        +            cR(-1)*crse(i,j-1,n) &
!$$$                        +            cR( 0)*crse(i,j  ,n) &
!$$$                        +            cR( 1)*crse(i,j+1,n) &
!$$$                        +            cR( 2)*crse(i,j+2,n))
                enddo
                do iry = 0, 1
                   jj = j*2+iry
                   if (jj.ge.fblo(2).and.jj.le.fbhi(2)) then
                      do i = cb2lo(1), cb2hi(1)
                         ii = 2*i
                         ftmp(ii  ) = 2.d0*(cL(-2)*ctmp(i-2,iry) &
                              +             cL(-1)*ctmp(i-1,iry) &
                              +             cL( 0)*ctmp(i  ,iry) &
                              +             cL( 1)*ctmp(i+1,iry) &
                              +             cL( 2)*ctmp(i+2,iry))
                         ftmp(ii+1) = 2.d0*ctmp(i,iry) - ftmp(ii)
!$$$                         ftmp(ii+1) = 2.d0*(cR(-2)*ctmp(i-2,iry) &
!$$$                              +             cR(-1)*ctmp(i-1,iry) &
!$$$                              +             cR( 0)*ctmp(i  ,iry) &
!$$$                              +             cR( 1)*ctmp(i+1,iry) &
!$$$                              +             cR( 2)*ctmp(i+2,iry))
                      enddo
                      do ii = fblo(1), fbhi(1)
                         fine(ii,jj,n) = ftmp(ii)
                      enddo
                   endif
                enddo
             enddo
          enddo
       else if (lratiox.eq.4 .and. lratioy.eq.4) then
          call amrex_error('AMREX_QUARTINTERP: refinement ratio = 4 TODO')
       else
          call amrex_error('AMREX_QUARTINTERP: unsupported refinement ratio')
       endif

     end subroutine AMREX_QUARTINTERP

end module amrex_interp_module

