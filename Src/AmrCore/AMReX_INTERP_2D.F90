
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

end module amrex_interp_module

