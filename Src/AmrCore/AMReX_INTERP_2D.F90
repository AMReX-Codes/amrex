
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
! ::: cbinterp:  cell centered bilinear interpolation
! ::: 
! ::: NOTE: it is assumed that the coarse grid array is
! ::: large enough to define interpolated values
! ::: in the region fblo:fbhi on the fine grid
! ::: 
! ::: Inputs/Outputs
! ::: fine        <=>  (modify) fine grid array
! ::: fine_l1,fine_l2,fine_h1,fine_h2   =>  (const)  index limits of fine grid
! ::: fb_l1,fb_l2,fb_h1,fb_h2     =>  (const)  subregion of fine grid to get values
! ::: 
! ::: crse         =>  (const)  coarse grid data 
! ::: crse_l1,crse_l2,crse_h1,crse_h2   =>  (const)  index limits of coarse grid
! ::: 
! ::: lratio(2)    =>  (const)  refinement ratio between levels
! ::: nvar         =>  (const)  number of components in array
! ::: 
! ::: TEMPORARY ARRAYS
! ::: slx,sly,slxy =>  1-D slope arrays
! ::: strip        =>  1-D temp array
! ::: --------------------------------------------------------------
! ::: 
    subroutine AMREX_CBINTERP (crse, crse_l1,crse_l2,crse_h1,crse_h2, cb_l1,cb_l2,cb_h1,cb_h2, &
                              fine, fine_l1,fine_l2,fine_h1,fine_h2, fb_l1,fb_l2,fb_h1,fb_h2, &
                              lratiox, lratioy, nvar, &
                              sl, num_slp, strip, strip_lo, strip_hi, &
                              actual_comp,actual_state) bind(c,name='amrex_cbinterp')

      implicit none

      integer crse_l1,crse_l2,crse_h1,crse_h2
      integer cb_l1,cb_l2,cb_h1,cb_h2
      integer fine_l1,fine_l2,fine_h1,fine_h2
      integer fb_l1,fb_l2,fb_h1,fb_h2
      integer lratiox, lratioy, nvar
      integer num_slp
      integer actual_comp,actual_state
      integer strip_lo, strip_hi
      real(amrex_real)  fine(fine_l1:fine_h1,fine_l2:fine_h2, nvar)
      real(amrex_real)  crse(crse_l1:crse_h1,crse_l2:crse_h2, nvar)
      real(amrex_real)  sl(cb_l1:cb_h1,num_slp)
      real(amrex_real)  strip(strip_lo:strip_hi)

#define SLX 1
#define SLY 2
#define SLXY 3

      integer lx, ly, hratx, hraty, ic, jc, jfn, jfc, i, n
      real(amrex_real) x, y, denomx, denomy

      denomx = one/dble(2*lratiox)
      denomy = one/dble(2*lratioy)

      hratx = lratiox/2
      hraty = lratioy/2

      do n = 1, nvar 
         do jc = cb_l2, cb_h2-1 
            do ic = cb_l1, cb_h1-1
               sl(ic,SLX) = crse(ic+1,jc,n)-crse(ic,jc,n)
               sl(ic,SLY) = crse(ic,jc+1,n)-crse(ic,jc,n)
               sl(ic,SLXY) = crse(ic+1,jc+1,n)-crse(ic+1,jc,n) &
                    - crse(ic  ,jc+1,n)+crse(ic  ,jc,n)
            end do
            do ly = 0, lratioy-1 
               jfn = jc*lratioy + ly
               jfc = jfn + hraty
               if (jfc .ge. fb_l2  .and.  jfc .le. fb_h2) then
                  y = denomy*(two*ly + one)
                  do lx = 0, lratiox-1
                     do ic = cb_l1, cb_h1-1
                        i = ic*lratiox + lx
                        x = denomx*(two*lx + one)
                        strip(i) = crse(ic,jc,n) + x*sl(ic,SLX) + &
                                   y*sl(ic,SLY) + x*y*sl(ic,SLXY)
                     end do
                  end do
                  do i = fb_l1, fb_h1 
                     fine(i,jfc,n) = strip(i-hratx)
                  end do
               end if
            end do
         end do
      end do

    end subroutine AMREX_CBINTERP

#undef  SLX
#undef  SLY
#undef  SLXY

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
! ::: protect_interp:   redo interpolation if the result of linccinterp
! ::: generates under- or overshoots.
! ::: 
! ::: 
! ::: Inputs/Outputs
! ::: fine        <=>  (modify) fine grid array
! ::: flo,fhi      =>  (const)  index limits of fine grid
! ::: fblo,fbhi    =>  (const)  subregion of fine grid to get values
! ::: cblo,cbhi    =>  (const)  coarse equivalent of fblo,fbhi
! ::: nvar         =>  (const)  number of variables in state vector
! ::: lratio(3)    =>  (const)  refinement ratio between levels
! ::: 
! ::: crse         =>  (const)  coarse grid data widended by 1 zone
! ::: clo,chi      =>  (const)  index limits of crse grid
! :::
! ::: --------------------------------------------------------------
! ::: 

    subroutine AMREX_PROTECT_INTERP (fine, fine_l1,fine_l2,fine_h1,fine_h2, fblo, fbhi, &
                                    crse, crse_l1,crse_l2,crse_h1,crse_h2, cblo, cbhi, &
                                    fvcx, fvcy, &
                                    fb_l1, fb_l2, fb_h1, fb_h2, &
                                    cvcx, cvcy, &
                                    cb_l1, cb_l2, cb_h1, cb_h2, &
                                    fine_state, state_l1,state_l2,state_h1,state_h2, &
                                    nvar, lratiox, lratioy, bc) bind(c,name='amrex_protect_interp')

      implicit none

      integer fine_l1,fine_l2,fine_h1,fine_h2
      integer crse_l1,crse_l2,crse_h1,crse_h2
      integer state_l1,state_l2,state_h1,state_h2
      integer fblo(2), fbhi(2)
      integer cblo(2), cbhi(2)
      integer fb_l1, fb_l2, fb_h1, fb_h2
      integer cb_l1, cb_l2, cb_h1, cb_h2
      integer lratiox, lratioy, nvar
      integer bc(2,2,nvar)
      real(amrex_real) fine(fine_l1:fine_h1,fine_l2:fine_h2,nvar)
      real(amrex_real) crse(crse_l1:crse_h1,crse_l2:crse_h2, nvar)
      real(amrex_real) fine_state(state_l1:state_h1,state_l2:state_h2, nvar)
      real(amrex_real) fvcx(fb_l1:fb_h1)
      real(amrex_real) fvcy(fb_l2:fb_h2)
      real(amrex_real) cvcx(cb_l1:cb_h1)
      real(amrex_real) cvcy(cb_l2:cb_h2)

      integer rMAX
      parameter (rMAX = 32)
      real(amrex_real) alpha, sumN, sumP, negVal, posVal
      real(amrex_real) crseTot, crseTotnew
      real(amrex_real) orig_fine(0:rMAX-1,0:rMAX-1)
      real(amrex_real) fvol,cvol
      integer redo_me
      integer ilo,ihi,jlo,jhi
      integer i,j,ic,jc,n
      integer icase

      if (MAX(lratiox,lratioy).gt.rMAX) then
#ifdef AMREX_DEBUG
         print *,'rMAX in INTERP_2D::AMREX_PROTECT_INTERP must be >= ',MAX(lratiox,lratioy)
#endif
         call bl_abort("rMAX in INTERP_2D")
      endif

      do jc = cblo(2), cbhi(2)
      do ic = cblo(1), cbhi(1)

         ilo = max(lratiox*ic            ,fine_l1)
         ihi = min(lratiox*ic+(lratiox-1),fine_h1)
         jlo = max(lratioy*jc            ,fine_l2)
         jhi = min(lratioy*jc+(lratioy-1),fine_h2)

         do n = 2, nvar-1

            redo_me = 0
            do j = jlo,jhi
            do i = ilo,ihi
               if ((fine_state(i,j,n)+fine(i,j,n)) .lt. 0.d0) redo_me = 1
            enddo
            enddo

! ****************************************************************************************
!
!           If all the fine values are non-negative after the original interpolated 
!            correction, then we do nothing here.
!
!           If any of the fine values are negative after the original interpolated
!            correction, then we do our best.
!
!           Special cases:
!
!             1) Coarse correction > 0, and fine_state has some cells with 
!                negative values which will be filled before adding to the other cells.
!                Use the correction to bring negative cells to zero, then
!                distribute the remaining positive proportionally.
!
!             2) Coarse correction > 0, and correction can not make them all
!                positive.  Add correction only to the negative cells, in proportion
!                to their magnitude.
!
!             3) Coarse correction < 0, and fine_state DOES NOT have enough
!                  have enough positive state to absorb it.  Here we bring
!                  all the positive fine cells to zero then distribute the remaining
!                  negative amount in such a way as to make them all as close to the
!                  same negative value as possible.
!
!             4) Coarse correction < 0, fine_state has enough
!                  positive state to absorb it without making any fine 
!                  cells negative, BUT fine_state+fine is currently negative
!                  in at least one fine cell.  Here just take a constant percentage
!                  away from each positive and don't touch the negatives.
!
!             crseTot = volume-weighted sum of all interpolated values of the correction,
!                       which is equivalent to the total volume-weighted coarse correction
!             SumN = volume-weighted sum of all negative values of fine_state
!             SumP = volume-weighted sum of all positive values of fine_state
!
! ****************************************************************************************

            if (redo_me .eq. 1) then

               icase = 0

               do j = jlo,jhi
               do i = ilo,ihi
                  orig_fine(i-ilo,j-jlo) = fine(i,j,n)
               enddo
               enddo

               crseTot = 0.d0
               do j = jlo,jhi
               do i = ilo,ihi
                  fvol = (fvcx(i+1)-fvcx(i)) * (fvcy(j+1)-fvcy(j))
                  crseTot = crseTot + fvol * fine(i,j,n)
               enddo
               enddo

               cvol = (cvcx(ic+1)-cvcx(ic)) * (cvcy(jc+1)-cvcy(jc))

               sumN = zero
               sumP = zero
               do j = jlo,jhi
               do i = ilo,ihi
                  fvol = (fvcx(i+1)-fvcx(i)) * (fvcy(j+1)-fvcy(j))
                  if (fine_state(i,j,n) .le. 0.d0) then
                    sumN = SumN + fvol * fine_state(i,j,n)
                  else
                    sumP = sumP + fvol * fine_state(i,j,n)
                  endif
               enddo
               enddo

               if (crseTot .gt. 0.d0 .and. crseTot .ge. abs(sumN)) then
!              Here we want to fill in the negative values first, then add
!                the remaining positive proportionally.

                   icase = 1
                   do j = jlo,jhi
                   do i = ilo,ihi
                      if (fine_state(i,j,n) .le. 0.d0) then
                        fine(i,j,n) = -fine_state(i,j,n)
                      endif
                   enddo
                   enddo

                   if (sumP > 0.d0) then

                     alpha = (crseTot - abs(sumN)) / sumP

                     do j = jlo,jhi
                     do i = ilo,ihi
                       if (fine_state(i,j,n) .ge. 0.d0) then
                         fine(i,j,n) = alpha * fine_state(i,j,n)
                       endif
                     enddo
                     enddo

                   else

                     posVal = (crseTot - abs(sumN)) / cvol

                     do j = jlo,jhi
                     do i = ilo,ihi
                       fine(i,j,n) = fine(i,j,n) + posVal
                     enddo
                     enddo

                   endif
            
                 endif

               if (crseTot .gt. 0.d0 .and. crseTot .lt. abs(sumN)) then
!              Here we don't have enough positive correction to fill all the
!                negative values of state, so we just try to fill them proportionally
!                and don't add any correction to the states already positive.

                   icase = 2
                   alpha = crseTot / abs(sumN)

                   do j = jlo,jhi
                   do i = ilo,ihi
                     if (fine_state(i,j,n) .lt. 0.d0) then
                       fine(i,j,n) = alpha * abs(fine_state(i,j,n))
                     else 
                       fine(i,j,n) = 0.d0
                     endif
                   enddo
                   enddo

               endif

               if (crseTot .lt. 0.d0 .and. abs(crseTot) .gt. sumP) then
!              Here we don't have enough positive states to absorb all the
!                negative correction, so we want to end up with all the fine
!                cells having the same negative value.

                   icase = 3
                   negVal = (sumP + sumN + crseTot)/cvol

                   do j = jlo,jhi
                   do i = ilo,ihi
                      fine(i,j,n) = negVal - fine_state(i,j,n)
                   enddo
                   enddo

               endif

               if (crseTot .lt. 0.d0 .and. abs(crseTot) .lt. sumP &
                                     .and. (sumP+sumN+crseTot) .gt. 0.d0) then
!              Here we have enough positive states to absorb all the
!                negative correction *and* redistribute to make negative cells
!                positive. 

                   icase = 4
                   alpha = (crseTot + sumN) / sumP

                   do j = jlo,jhi
                   do i = ilo,ihi
                      if (fine_state(i,j,n) .lt. 0.d0) then
                        fine(i,j,n) = -fine_state(i,j,n)
                      else
                        fine(i,j,n) = alpha * fine_state(i,j,n)
                      endif  
                   enddo
                   enddo

               endif

               if (crseTot .lt. 0.d0 .and. abs(crseTot) .lt. sumP &
                                     .and. (sumP+sumN+crseTot) .le. 0.d0) then
!              Here we have enough positive states to absorb all the
!                negative correction, but not to fix the states already negative. 
!                We bring all the positive states to zero, and use whatever 
!                remaining positiveness from the states to help the negative states.

                   icase = 5
                   alpha = (crseTot + sumP) / sumN

                   do j = jlo,jhi
                   do i = ilo,ihi
                      if (fine_state(i,j,n) .gt. 0.d0) then
                        fine(i,j,n) = -fine_state(i,j,n)
                      else 
                        fine(i,j,n) = alpha * fine_state(i,j,n)
                      endif
                   enddo
                   enddo

               endif

               crseTotnew   = 0.d0
               do j = jlo,jhi
               do i = ilo,ihi
                  fvol = (fvcx(i+1)-fvcx(i)) * (fvcy(j+1)-fvcy(j))
                  crseTotnew   = crseTotnew   + fvol * fine(i,j,n)
               enddo
               enddo

#ifdef AMREX_DEBUG
               if (abs(crseTotnew - crseTot)/cvol .gt. 1.e-8) then
                  print *,' '
                  print *,'BLEW CONSERVATION with ICASE = ',icase
                  print *,'AT COARSE CELL ',ic,jc,' AND COMPONENT ',n
                  print *,'CRSETOT NEW OLD ',crseTotnew, crseTot
                  print *,'CVOL ',cvol
                  print *,'SUMP SUMN ',sumP,sumN
                  do j = jlo,jhi
                  do i = ilo,ihi
                     fvol = (fvcx(i+1)-fvcx(i)) * (fvcy(j+1)-fvcy(j))
                     print *,'FINE OLD NEW ',i,j,orig_fine(i-ilo,j-jlo), &
                                             fine(i,j,n), fine_state(i,j,n), &
                                             fvol
                     if (abs(fvol) .lt. 1.d-50) then
                       print *,'MAKING FVOL ',fvcx(i+1),fvcx(i),fvcy(j+1),fvcy(j)
                     endif
                  enddo
                  enddo
               endif
#endif

!              do j = jlo,jhi
!              do i = ilo,ihi
!                 if ((fine_state(i,j,n) + fine(i,j,n)) .lt. 0.d0) then
!                    print *,'STILL NEGATIVE AT ',i,j,n
!                    print *,'AT COARSE CELL ',ic,jc
!                    print *,'FINE STATE ',fine_state(i,j,n)
!                    print *,'FINE CORRECTION ',fine(i,j,n)
!                    print *,'CRSETOT ',crseTot
!                    print *,'SUMN / SUMP ',sumN, sumP
!                    print *,' '
!                 endif
!              enddo
!              enddo
!              enddo

!           End (if redo .eq. 1)
            endif

         enddo

!     Set sync for density (n=1) to sum of spec sync (2:nvar-1)
         do j = jlo,jhi
         do i = ilo,ihi
            fine(i,j,1) = 0.d0
            do n = 2,nvar-1
               fine(i,j,1) = fine(i,j,1) + fine(i,j,n)
            enddo
         enddo
         enddo

!     End of coarse index loops
      enddo
      enddo

    end subroutine AMREX_PROTECT_INTERP

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

