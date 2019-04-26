
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
! ::: fine_l1,fine_l2,fine_l3,fine_h1,fine_h2,fine_h3   =>  (const)  index limits of fine grid
! ::: fb_l1,fb_l2,fb_l3,fb_h1,fb_h2,fb_h3     =>  (const)  subregion of fine grid to get values
! ::: 
! ::: crse         =>  (const)  coarse grid data 
! ::: crse_l1,crse_l2,crse_l3,crse_h1,crse_h2,crse_h3   =>  (const)  index limits of coarse grid
! ::: 
! ::: lratio(3)    =>  (const)  refinement ratio between levels
! ::: nvar         =>  (const)  number of components in array
! ::: 
! ::: TEMPORARY ARRAYS
! ::: slx,sly,slxy =>  1-D slope arrays
! ::: strip        =>  1-D temp array
! ::: --------------------------------------------------------------
! ::: 
    subroutine AMREX_CBINTERP (crse, crse_l1,crse_l2,crse_l3,crse_h1,crse_h2,crse_h3, cb_l1,cb_l2,cb_l3,cb_h1,cb_h2,cb_h3, &
                              fine, fine_l1,fine_l2,fine_l3,fine_h1,fine_h2,fine_h3, fb_l1,fb_l2,fb_l3,fb_h1,fb_h2,fb_h3, &
                              lratiox, lratioy, lratioz, nvar, &
                              sl, num_slp, strip, strip_lo, strip_hi, &
                              actual_comp, actual_state) bind(c,name='amrex_cbinterp')

      implicit none

      integer crse_l1,crse_l2,crse_l3,crse_h1,crse_h2,crse_h3
      integer cb_l1,cb_l2,cb_l3,cb_h1,cb_h2,cb_h3
      integer fine_l1,fine_l2,fine_l3,fine_h1,fine_h2,fine_h3
      integer fb_l1,fb_l2,fb_l3,fb_h1,fb_h2,fb_h3
      integer lratiox, lratioy, lratioz, nvar
      integer num_slp
      integer strip_lo, strip_hi
      integer actual_comp,actual_state
      real(amrex_real) fine(fine_l1:fine_h1,fine_l2:fine_h2,fine_l3:fine_h3,nvar)
      real(amrex_real) crse(crse_l1:crse_h1,crse_l2:crse_h2,crse_l3:crse_h3,nvar)
      real(amrex_real) strip(strip_lo:strip_hi)
      real(amrex_real) sl(cb_l1:cb_h1, num_slp)

#define SLX 1
#define SLY 2
#define SLZ 3
#define SLXY 4
#define SLXZ 5
#define SLYZ 6
#define SLXYZ 7

      integer lx, ly, lz, hratx, hraty, hratz, ic, jc, kc, jfn, jfc, kfn, kfc, i, n
      real(amrex_real) x, y, z, denomx, denomy, denomz

      denomx = one/dble(2*lratiox)
      denomy = one/dble(2*lratioy)
      denomz = one/dble(2*lratioz)

      hratx = lratiox/2
      hraty = lratioy/2
      hratz = lratioz/2

      do n = 1, nvar
         do kc = cb_l3, cb_h3-1
         do jc = cb_l2, cb_h2-1
            do ic = cb_l1, cb_h1-1
               sl(ic,SLX  ) = crse(ic+1,jc,  kc,  n)-crse(ic,  jc,  kc,  n)
               sl(ic,SLY  ) = crse(ic,  jc+1,kc,  n)-crse(ic,  jc,  kc,  n)
               sl(ic,SLZ  ) = crse(ic,  jc,  kc+1,n)-crse(ic,  jc,  kc,  n)
               sl(ic,SLXY ) = crse(ic+1,jc+1,kc,  n)-crse(ic+1,jc,  kc,  n) &
                            - crse(ic  ,jc+1,kc,  n)+crse(ic  ,jc,  kc,  n)
               sl(ic,SLXZ ) = crse(ic+1,jc,  kc+1,n)-crse(ic+1,jc,  kc,  n) &
                            - crse(ic,  jc,  kc+1,n)+crse(ic,  jc,  kc,  n)
               sl(ic,SLYZ ) = crse(ic,  jc+1,kc+1,n)-crse(ic,  jc,  kc+1,n) &
                            - crse(ic  ,jc+1,kc,  n)+crse(ic  ,jc,  kc,  n)
               sl(ic,SLXYZ) = crse(ic+1,jc,  kc,  n)+crse(ic,  jc+1,kc,  n) &
                            + crse(ic,  jc,  kc+1,n)-crse(ic,  jc+1,kc+1,n) &
                            - crse(ic+1,jc,  kc+1,n)-crse(ic+1,jc+1,kc,  n) &
                            + crse(ic+1,jc+1,kc+1,n)-crse(ic,  jc,  kc,  n)

            end do

            do lz = 0, lratioz-1
               kfn = kc*lratioz + lz
               kfc = kfn + hratz
               if (kfc .ge. fb_l3 .and. kfc .le. fb_h3) then
                  z = denomz*(two*lz + one)
                  do ly = 0, lratioy-1
                     jfn = jc*lratioy + ly
                     jfc = jfn + hraty
                     if (jfc .ge. fb_l2  .and.  jfc .le. fb_h2) then
                        y = denomy*(two*ly + one)
                        do lx = 0, lratiox-1
                           do ic = cb_l1, cb_h1-1
                              i = ic*lratiox + lx
                              x = denomx*(two*lx + one)
                              strip(i) = crse(ic,jc,kc,n) +     x*sl(ic,SLX  ) + &
                                           y*sl(ic,SLY )  +     z*sl(ic,SLZ  ) + &
                                         x*y*sl(ic,SLXY)  +   x*z*sl(ic,SLXZ ) + &
                                         y*z*sl(ic,SLYZ)  + x*y*z*sl(ic,SLXYZ)
                           end do
                        end do
                        do i = fb_l1, fb_h1
                           fine(i,jfc,kfc,n) = strip(i-hratx)
                        end do
                     end if
                  end do
               end if
            end do
         end do
         end do
      end do

    end subroutine AMREX_CBINTERP

    subroutine AMREX_CQINTERP (fine, fine_l1,fine_l2,fine_l3,fine_h1,fine_h2,fine_h3, &
                              fb_l1,fb_l2,fb_l3,fb_h1,fb_h2,fb_h3, &
                              nvar, lratiox, lratioy, lratioz, crse, &
                              clo, chi, cb_l1,cb_l2,cb_l3,cb_h1,cb_h2,cb_h3, &
                              fslo, fshi, cslope, clen, fslope, fdat, &
                              flen, voff, bc, limslope, &
                              fvcx, fvcy, fvcz, cvcx, cvcy, cvcz, &
                              actual_comp, actual_state) bind(c,name='amrex_cqinterp')

      implicit none

      integer fine_l1,fine_l2,fine_l3,fine_h1,fine_h2,fine_h3
      integer fb_l1,fb_l2,fb_l3,fb_h1,fb_h2,fb_h3
      integer cb_l1,cb_l2,cb_l3,cb_h1,cb_h2,cb_h3
      integer fslo(3), fshi(3)
      integer nvar, lratiox, lratioy, lratioz
      integer bc(3,2,nvar)
      integer clen, flen, clo, chi, limslope
      integer actual_comp,actual_state
      real(amrex_real) fine(fine_l1:fine_h1,fine_l2:fine_h2,fine_l3:fine_h3,nvar)
      real(amrex_real) crse(clo:chi, nvar)
      real(amrex_real) cslope(clo:chi, 3)
      real(amrex_real) fslope(flen, 3)
      real(amrex_real) fdat(flen)
      real(amrex_real) voff(flen)
      real(amrex_real) fvcx(fb_l1:fb_h1+1)
      real(amrex_real) fvcy(fb_l2:fb_h2+1)
      real(amrex_real) fvcz(fb_l3:fb_h3+1)
      real(amrex_real) cvcx(cb_l1:cb_h1+1)
      real(amrex_real) cvcy(cb_l2:cb_h2+1)
      real(amrex_real) cvcz(cb_l3:cb_h3+1)

      call bl_abort('QUADRATIC INTERP NOT IMPLEMEMNTED IN 3-D')

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
    subroutine AMREX_PROTECT_INTERP (fine, fine_l1,fine_l2,fine_l3,fine_h1,fine_h2,fine_h3, fblo, fbhi, &
                                    crse, crse_l1,crse_l2,crse_l3,crse_h1,crse_h2,crse_h3, cblo, cbhi, &
                                    fine_state, state_l1,state_l2,state_l3,state_h1,state_h2,state_h3, &
                                    nvar, lratiox,lratioy,lratioz, bc) bind(c,name='amrex_protect_interp')

      implicit none

      integer fine_l1,fine_l2,fine_l3,fine_h1,fine_h2,fine_h3
      integer crse_l1,crse_l2,crse_l3,crse_h1,crse_h2,crse_h3
      integer state_l1,state_l2,state_l3,state_h1,state_h2,state_h3
      integer fblo(3), fbhi(3)
      integer cblo(3), cbhi(3)
      integer lratiox, lratioy, lratioz, nvar
      integer bc(3,2,nvar)
      real(amrex_real) fine(fine_l1:fine_h1,fine_l2:fine_h2,fine_l3:fine_h3,nvar)
      real(amrex_real) crse(crse_l1:crse_h1,crse_l2:crse_h2,crse_l3:crse_h3, nvar)
      real(amrex_real) fine_state(state_l1:state_h1,state_l2:state_h2,state_l3:state_h3, nvar)

      integer rMAX
      parameter (rMAX = 16)
      real(amrex_real) alpha, sumN, sumP, crseTot, negVal, posVal
      real(amrex_real) sum_fine_new,sum_fine_old
      real(amrex_real) orig_fine(0:rMAX-1,0:rMAX-1,0:rMAX-1)
      integer redo_me
      integer ilo,ihi,jlo,jhi,klo,khi
      integer i,j,k,ic,jc,kc,n
      integer numFineCells
      integer icase

      if (MAX(lratiox,lratioy,lratioz).gt.rMAX) then
!         print *,'rMAX in INTERP_3D::AMREX_PROTECT_INTERP must be >= ', &
!              MAX(lratiox,lratioy,lratioz)
         call bl_abort("rMAX in INTERP_3D")
      endif

      do kc = cblo(3), cbhi(3)
      do jc = cblo(2), cbhi(2)
      do ic = cblo(1), cbhi(1)

         ilo = max(lratiox*ic            ,fine_l1)
         ihi = min(lratiox*ic+(lratiox-1),fine_h1)
         jlo = max(lratioy*jc            ,fine_l2)
         jhi = min(lratioy*jc+(lratioy-1),fine_h2)
         klo = max(lratioz*kc            ,fine_l3)
         khi = min(lratioz*kc+(lratioz-1),fine_h3)

         do n = 2, nvar-1

            redo_me = 0
            do k = klo,khi
            do j = jlo,jhi
            do i = ilo,ihi
               if ((fine_state(i,j,k,n)+fine(i,j,k,n)) .lt. 0.d0) redo_me = 1
            enddo
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
!             crseTot = sum of all interpolated values of the correction,
!                       which is equivalent to the coarse correction * ratio**3
!             SumN = sum of all negative values of fine_state
!             SumP = sum of all positive values of fine_state
!
! ****************************************************************************************

            if (redo_me .eq. 1) then

               icase = 0
               sum_fine_old = 0.d0
               do k = klo,khi
               do j = jlo,jhi
               do i = ilo,ihi
                  sum_fine_old = sum_fine_old + fine(i,j,k,n)
                  orig_fine(i-ilo,j-jlo,k-klo) = fine(i,j,k,n)
               enddo
               enddo
               enddo

               crseTot = sum_fine_old
               numFineCells = (ihi-ilo+1) * (jhi-jlo+1) * (khi-klo+1)

               sumN = zero
               sumP = zero
               do k = klo,khi
               do j = jlo,jhi
               do i = ilo,ihi
                  if (fine_state(i,j,k,n) .le. 0.d0) then
                    sumN = SumN + fine_state(i,j,k,n)
                  else
                    sumP = sumP + fine_state(i,j,k,n)
                  endif
               enddo
               enddo
               enddo

               if (crseTot .gt. 0.d0 .and. crseTot .ge. abs(sumN)) then
!              Here we want to fill in the negative values first, then add
!                the remaining positive proportionally.

                   icase = 1
                   do k = klo,khi
                   do j = jlo,jhi
                   do i = ilo,ihi
                      if (fine_state(i,j,k,n) .le. 0.d0) then
                        fine(i,j,k,n) = -fine_state(i,j,k,n)
                      endif
                   enddo
                   enddo
                   enddo

                   if (sumP .gt. 0.d0) then

                    alpha = (crseTot - abs(sumN)) / sumP

                    do k = klo,khi
                    do j = jlo,jhi
                    do i = ilo,ihi
                       if (fine_state(i,j,k,n) .ge. 0.d0) then
                         fine(i,j,k,n) = alpha * fine_state(i,j,k,n)
                       endif
                    enddo
                    enddo
                    enddo

                  else

                    posVal = (crseTot - abs(sumN)) / float(numFineCells)

                    do k = klo,khi
                    do j = jlo,jhi
                    do i = ilo,ihi
                       fine(i,j,k,n) = fine(i,j,k,n) + posVal
                    enddo
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

                   do k = klo,khi
                   do j = jlo,jhi
                   do i = ilo,ihi
                      if (fine_state(i,j,k,n) .lt. 0.d0) then
                        fine(i,j,k,n) = alpha * abs(fine_state(i,j,k,n))
                      else 
                        fine(i,j,k,n) = 0.d0
                      endif
                   enddo
                   enddo
                   enddo

               endif

               if (crseTot .lt. 0.d0 .and. abs(crseTot) .gt. sumP) then
!              Here we don't have enough positive states to absorb all the
!                negative correction, so we want to end up with all the fine
!                cells having the same negative value.

                   icase = 3
                   negVal = (sumP + sumN + crseTot)/float(numFineCells)

                   do k = klo,khi
                   do j = jlo,jhi
                   do i = ilo,ihi
                      fine(i,j,k,n) = negVal - fine_state(i,j,k,n)
                   enddo
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

                   do k = klo,khi
                   do j = jlo,jhi
                   do i = ilo,ihi
                      if (fine_state(i,j,k,n) .lt. 0.d0) then
                        fine(i,j,k,n) = -fine_state(i,j,k,n)
                      else
                        fine(i,j,k,n) = alpha * fine_state(i,j,k,n)
                      endif  
                   enddo
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

                   do k = klo,khi
                   do j = jlo,jhi
                   do i = ilo,ihi
                      if (fine_state(i,j,k,n) .gt. 0.d0) then
                        fine(i,j,k,n) = -fine_state(i,j,k,n)
                      else 
                        fine(i,j,k,n) = alpha * fine_state(i,j,k,n)
                      endif
                   enddo
                   enddo
                   enddo

               endif

               sum_fine_new = 0.d0
               do k = klo,khi
               do j = jlo,jhi
               do i = ilo,ihi
                  sum_fine_new = sum_fine_new + fine(i,j,k,n)
               enddo
               enddo
               enddo

#ifdef AMREX_DEBUG
               if (abs(sum_fine_new - sum_fine_old) .gt. 1.e-8) then
                  print *,' '
                  print *, &
                       'PROTECT_INTERP: BLEW CONSERVATION with ICASE = ' &
                       ,icase
                  print *,'AT COARSE CELL ',ic,jc,kc,' AND COMPONENT ',n
                  print *,'NEW SUM / OLD SUM ', &
                       sum_fine_new, sum_fine_old
                  print *,'CRSETOT ',crseTot
                  print *,'SUMP SUMN ',sumP,sumN
                  do k = klo,khi
                  do j = jlo,jhi
                  do i = ilo,ihi
                     print *,'FINE OLD NEW ', &
                          i,j,k,orig_fine(i-ilo,j-jlo,k-klo), &
                          fine(i,j,k,n), fine_state(i,j,k,n)
                  enddo
                  enddo
                  enddo
               endif
#endif

!              do k = klo,khi
!              do j = jlo,jhi
!              do i = ilo,ihi
!                 if ((fine_state(i,j,k,n) + fine(i,j,k,n)) .lt. 0.d0) then
!                    print *,'STILL NEGATIVE AT ',i,j,k,n
!                    print *,'AT COARSE CELL ',ic,jc,kc
!                    print *,'FINE STATE ',fine_state(i,j,k,n)
!                    print *,'FINE CORRECTION ',fine(i,j,k,n)
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
         do k = klo,khi
         do j = jlo,jhi
         do i = ilo,ihi
            fine(i,j,k,1) = 0.d0
            do n = 2,nvar-1
               fine(i,j,k,1) = fine(i,j,k,1) + fine(i,j,k,n)
            enddo
         enddo
         enddo
         enddo

!     End of coarse index loops
      enddo
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
! ::: lratio[xyz]  =>  (const)  refinement ratio between levels
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
! ::: ctmp2        =>  2-D temp array
! ::: --------------------------------------------------------------
! ::: 
     subroutine AMREX_QUARTINTERP (fine, fine_l1,fine_l2,fine_l3,fine_h1,fine_h2,fine_h3, &
                                  fblo, fbhi, fb2lo, fb2hi, &
                                  crse, crse_l1,crse_l2,crse_l3,crse_h1,crse_h2,crse_h3, &
                                  cblo, cbhi, cb2lo, cb2hi, &
                                  nvar, &
                                  lratiox, lratioy, lratioz, &
                                  ftmp, ctmp, ctmp2, &
                                  bc,actual_comp,actual_state) bind(c,name='amrex_quartinterp')

       implicit none

       integer fine_l1,fine_l2,fine_l3,fine_h1,fine_h2,fine_h3
       integer crse_l1,crse_l2,crse_l3,crse_h1,crse_h2,crse_h3
       integer fblo(3), fbhi(3), fb2lo(3), fb2hi(3)
       integer cblo(3), cbhi(3), cb2lo(3), cb2hi(3)
       integer nvar,lratiox,lratioy,lratioz
       integer bc(3,2,nvar)
       integer actual_comp,actual_state
       real(amrex_real) fine(fine_l1:fine_h1,fine_l2:fine_h2,fine_l3:fine_h3,nvar)
       real(amrex_real) crse(crse_l1:crse_h1,crse_l2:crse_h2,crse_l3:crse_h3,nvar)
       real(amrex_real) ftmp(fb2lo(1):fb2hi(1))
       real(amrex_real) ctmp(cblo(1):cbhi(1),0:lratioy-1)
       real(amrex_real) ctmp2(cblo(1):cbhi(1),cblo(2):cbhi(2),0:lratioz-1)

!      Local variables
       integer i,j,k,ii,jj,kk,n,iry,irz
       real(amrex_real) cL(-2:2)
!       real(amrex_real) cR(-2:2)
       data cL/ -0.01171875D0,  0.0859375D0, 0.5d0, -0.0859375D0, &
                 0.01171875D0 /
!       data cR/  0.01171875D0, -0.0859375D0, 0.5d0,  0.0859375D0, &
!                -0.01171875D0 /
       
       if (lratiox.eq.2 .and. lratioy.eq.2 .and. lratioz.eq.2) then

          do n = 1, nvar
          do k = cb2lo(3), cb2hi(3)

             do j = cblo(2), cbhi(2)
             do i = cblo(1), cbhi(1)
                ctmp2(i,j,0) = 2.d0*(cL(-2)*crse(i,j,k-2,n) &
                     +               cL(-1)*crse(i,j,k-1,n) &
                     +               cL( 0)*crse(i,j,k  ,n) &
                     +               cL( 1)*crse(i,j,k+1,n) &
                     +               cL( 2)*crse(i,j,k+2,n))
                ctmp2(i,j,1) = 2.d0*crse(i,j,k,n) - ctmp2(i,j,0)
!$$$                ctmp2(i,j,1) = 2.d0*(cR(-2)*crse(i,j,k-2,n) & 
!$$$                     +               cR(-1)*crse(i,j,k-1,n) &
!$$$                     +               cR( 0)*crse(i,j,k  ,n) &
!$$$                     +               cR( 1)*crse(i,j,k+1,n) &
!$$$                     +               cR( 2)*crse(i,j,k+2,n))
             enddo
             enddo

             do irz = 0, 1
                kk = k*2+irz
                if (kk.ge.fblo(3) .and. kk.le.fbhi(3)) then

                   do j = cb2lo(2), cb2hi(2)

                      do i = cblo(1), cbhi(1)
                         ctmp(i,0) = 2.d0*(cL(-2)*ctmp2(i,j-2,irz) & 
                              +            cL(-1)*ctmp2(i,j-1,irz) &
                              +            cL( 0)*ctmp2(i,j  ,irz) &
                              +            cL( 1)*ctmp2(i,j+1,irz) &
                              +            cL( 2)*ctmp2(i,j+2,irz))
                         ctmp(i,1) = 2.d0*ctmp2(i,j,irz) - ctmp(i,0)
!$$$                         ctmp(i,1) = 2.d0*(cR(-2)*ctmp2(i,j-2,irz) & 
!$$$                              +            cR(-1)*ctmp2(i,j-1,irz) &
!$$$                              +            cR( 0)*ctmp2(i,j  ,irz) &
!$$$                              +            cR( 1)*ctmp2(i,j+1,irz) &
!$$$                              +            cR( 2)*ctmp2(i,j+2,irz))
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
!$$$                               ftmp(ii+1) = 2.d0*(cR(-2)*ctmp(i-2,iry) &
!$$$                                    +             cR(-1)*ctmp(i-1,iry) &
!$$$                                    +             cR( 0)*ctmp(i  ,iry) &
!$$$                                    +             cR( 1)*ctmp(i+1,iry) &
!$$$                                    +             cR( 2)*ctmp(i+2,iry))
                            enddo
                            do ii = fblo(1), fbhi(1)
                               fine(ii,jj,kk,n) = ftmp(ii)
                            enddo
                         endif  ! if (jj.ge.......
                      enddo  ! do iry

                   enddo  ! do j

                endif  ! if (kk.ge.......
             enddo  ! do irz

          enddo  ! do k
          enddo  ! do n

       else if (lratiox.eq.4 .and. lratioy.eq.4 .and. lratioz.eq.4) then
          call amrex_error('AMREX_QUARTINTERP: refinement ratio = 4 TODO')
       else
          call amrex_error('AMREX_QUARTINTERP: unsupported refinement ratio')
       endif

     end subroutine AMREX_QUARTINTERP

end module amrex_interp_module
