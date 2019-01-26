! ::
! :: ----------------------------------------------------------
! ::
    subroutine enforce_minimum_density(uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
                                       uout,uout_l1,uout_l2,uout_l3, &
                                       uout_h1,uout_h2,uout_h3, &
                                       lo,hi,print_fortran_warnings)

      use amrex_fort_module, only : rt => amrex_real
      use network, only : nspec, naux
      use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, &
                                     UFS, UFA, small_dens, nadv

      implicit none

      integer          :: lo(3), hi(3), print_fortran_warnings
      integer          ::  uin_l1,  uin_l2,  uin_l3,  uin_h1,  uin_h2,  uin_h3
      integer          :: uout_l1, uout_l2, uout_l3, uout_h1, uout_h2, uout_h3
      real(rt) ::  uin( uin_l1: uin_h1, uin_l2: uin_h2, uin_l3: uin_h3,NVAR)
      real(rt) :: uout(uout_l1:uout_h1,uout_l2:uout_h2,uout_l3:uout_h3,NVAR)

      ! Local variables
      integer          :: i,ii,ilo,ihi
      integer          :: j,jj,jlo,jhi
      integer          :: k,kk,klo,khi
      integer          :: n,nmax
      logical          :: do_diag
      real(rt) :: min_dens
      real(rt) :: sum_state (NVAR)
      real(rt) :: sum_before(NVAR)
      real(rt) :: sum_after (NVAR)
      real(rt) :: min_vel(3), max_vel(3)
      real(rt) :: min_e     , max_e
      real(rt) :: delta_mass,frac,omfrac
      real(rt) :: delta_xmom, delta_ymom, delta_zmom
      real(rt) :: delta_rhoe
      real(rt) :: new_xvel, new_yvel, new_zvel, new_e
      real(rt) :: temp_sum, temp_num

      if (UFS .gt. 0) then
          nmax = UFS+nspec+naux-1
      else if (UFS .le. 0 .and. UFA .gt. 0) then
          nmax = UFA+nadv-1
      else 
          nmax = UEINT
      end if

      ! do_diag = .true.
      do_diag = .false.

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               if (uout(i,j,k,URHO) < small_dens) then

                  ilo = max(i-1,lo(1))
                  jlo = max(j-1,lo(2))
                  klo = max(k-1,lo(3))
                  ihi = min(i+1,hi(1))
                  jhi = min(j+1,hi(2))
                  khi = min(k+1,hi(3))

                  if (do_diag) then
                      print *,'old xvel      ',uout(i,j,k,UMX  )/uout(i,j,k,URHO)
                      print *,'old yvel      ',uout(i,j,k,UMY  )/uout(i,j,k,URHO)
                      print *,'old zvel      ',uout(i,j,k,UMZ  )/uout(i,j,k,URHO)
                      print *,'old   e       ',uout(i,j,k,UEINT)/uout(i,j,k,URHO)
                      if (UFS .gt. 0) then
                          print *,'old ufs       ',uout(i,j,k,UFS  )/uout(i,j,k,URHO)
                          print *,'old ufsp1     ',uout(i,j,k,UFS+1)/uout(i,j,k,URHO)
                      end if
                  end if

                  ! Find the minimum density and the sum of all the conserved
                  !      quantities of the non-corner neighbors
                  min_dens      =  1.d100
                  sum_state(:)  =  0.d0
                  sum_before(:) =  0.d0
                  min_vel(:)    =  1.d100
                  max_vel(:)    = -1.d100
                  min_e         =  1.d100
                  max_e         = -1.d100
                  do kk = klo,khi
                  do jj = jlo,jhi
                  do ii = ilo,ihi
                      if ( (i-ii)*(j-jj)*(k-kk) .eq. 0 .and. &
                          uout(ii,jj,kk,URHO).gt.small_dens) then
                          min_dens = min(min_dens,uout(ii,jj,kk,URHO))
                          min_vel(1:3) = min(min_vel(1:3),uout(ii,jj,kk,UMX:UMZ)/uout(ii,jj,kk,URHO))
                          max_vel(1:3) = max(max_vel(1:3),uout(ii,jj,kk,UMX:UMZ)/uout(ii,jj,kk,URHO))
                          min_e        = min(min_e       ,uout(ii,jj,kk,UEINT)  /uout(ii,jj,kk,URHO))
                          max_e        = max(max_e       ,uout(ii,jj,kk,UEINT)  /uout(ii,jj,kk,URHO))
                          sum_state(URHO:nmax) = sum_state(URHO:nmax) + uout(ii,jj,kk,URHO:nmax)
                      end if
                      sum_before(URHO:nmax) = sum_before(URHO:nmax) + uout(ii,jj,kk,URHO:nmax)
                  end do
                  end do
                  end do

                  ! This is how much mass we are adding to cell (i,j,k), therefore
                  !         how much total mass we must subtract from the neighbors
                  delta_mass = min_dens - uout(i,j,k,URHO)

                  ! Subtract from the neighbors in proportion to their own mass -- 
                  !     this conserves the conserved quantities and keeps velocities,
                  !     e, and E unchanged.
                    frac = (delta_mass / sum_state(URHO))
                  omfrac = 1.d0 - frac
                  do kk = klo,khi
                  do jj = jlo,jhi
                  do ii = ilo,ihi
                      if ( (i-ii)*(j-jj)*(k-kk) .eq. 0 .and. &
                           uout(ii,jj,kk,URHO).gt.small_dens ) then
                         uout(ii,jj,kk,URHO:nmax) = uout(ii,jj,kk,URHO:nmax) * omfrac
                      end if
                  end do
                  end do
                  end do

                  if (print_fortran_warnings .gt. 0) then
                     !
                     ! A critical region since we usually can't write from threads.
                     !
                     if (uout(i,j,k,URHO) < 0.d0) then
                        print *,'   '
                        print *,'>>> RESETTING STATE with NEG.  DENSITY AT ',i,j,k
                        print *,'>>> FROM ',uout(i,j,k,URHO),' TO ',min_dens
                        print *,'   '
                     else
                        print *,'   '
                        print *,'>>> RESETTING STATE with SMALL DENSITY AT ',i,j,k
                        print *,'>>> FROM ',uout(i,j,k,URHO),' TO ',min_dens
                        print *,'   '
                     end if
                  end if

                  ! Now define the new state at (i,j,k)
                  uout(i,j,k,URHO    ) = min_dens
                  uout(i,j,k,UMX:nmax) = uout(i,j,k,UMX:nmax) + frac * sum_state(UMX:nmax)

                  ! ***** Done fixing the density -- now worry about the other quantities ***

                  ! Re-set the velocities to be the average of the neighbors in the same direction
                  ! For now don't worry about conservation of momentum
                  if (i-1.ge.lo(1) .and. i+1.le.hi(1)) then
                     new_xvel = 0.5d0*(uout(i+1,j,k,UMX)/uout(i+1,j,k,URHO) + &
                                       uout(i-1,j,k,UMX)/uout(i-1,j,k,URHO) )
                  else
                      temp_sum = 0.d0
                      temp_num = 0.d0
                      do kk = klo,khi
                      do jj = jlo,jhi
                      do ii = ilo,ihi
                         temp_sum = temp_sum +  uout(ii,jj,kk,UMX)/uout(ii,jj,kk,URHO)
                         temp_num = temp_num +  1.d0
                      end do
                      end do
                      end do
                      new_xvel = temp_sum / temp_num
                  end if
                  if (j-1.ge.lo(2) .and. j+1.le.hi(2)) then
                      new_yvel = 0.5d0*(uout(i,j+1,k,UMY)/uout(i,j+1,k,URHO) + &
                                        uout(i,j-1,k,UMY)/uout(i,j-1,k,URHO) )
                  else
                      temp_sum = 0.d0
                      temp_num = 0.d0
                      do kk = klo,khi
                      do jj = jlo,jhi
                      do ii = ilo,ihi
                         temp_sum = temp_sum +  uout(ii,jj,kk,UMY)/uout(ii,jj,kk,URHO)
                         temp_num = temp_num +  1.d0
                      end do
                      end do
                      end do
                      new_yvel = temp_sum / temp_num
                  end if
                  if (k-1.ge.lo(3) .and. k+1.le.hi(3)) then
                      new_zvel = 0.5d0*(uout(i,j,k+1,UMZ)/uout(i,j,k+1,URHO) + &
                                        uout(i,j,k-1,UMZ)/uout(i,j,k-1,URHO) )
                  else
                      temp_sum = 0.d0
                      temp_num = 0.d0
                      do kk = klo,khi
                      do jj = jlo,jhi
                      do ii = ilo,ihi
                         temp_sum = temp_sum +  uout(ii,jj,kk,UMZ)/uout(ii,jj,kk,URHO)
                         temp_num = temp_num +  1.d0
                      end do
                      end do
                      end do
                      new_zvel = temp_sum / temp_num
                  end if

                  ! Sum over all 27 cells except corners and center.
                  temp_sum = 0.d0
                  temp_num = 0.d0
                  do kk = klo,khi
                  do jj = jlo,jhi
                  do ii = ilo,ihi
                      if ( (i-ii)*(j-jj)*(k-kk) .eq. 0 .and. &
                         .not. (i.eq.ii .and. j.eq.jj .and. k.eq.kk) ) then
                         temp_sum = temp_sum +  uout(ii,jj,kk,UEINT)/uout(ii,jj,kk,URHO)
                         temp_num = temp_num +  1.d0
                      end if
                  end do
                  end do
                  end do
                  new_e = temp_sum / temp_num

                  delta_xmom = uout(i,j,k,URHO) * new_xvel  - uout(i,j,k,UMX)
                  delta_ymom = uout(i,j,k,URHO) * new_yvel  - uout(i,j,k,UMY)
                  delta_zmom = uout(i,j,k,URHO) * new_zvel  - uout(i,j,k,UMZ)
                  delta_rhoe = uout(i,j,k,URHO) * new_e     - uout(i,j,k,UEINT)

                  uout(i,j,k,UMX) =  uout(i,j,k,URHO) * new_xvel
                  uout(i,j,k,UMY) =  uout(i,j,k,URHO) * new_yvel
                  uout(i,j,k,UMZ) =  uout(i,j,k,URHO) * new_zvel

                  uout(i,j,k,UEINT) = uout(i,j,k,URHO) * new_e
                  uout(i,j,k,UEDEN) = uout(i,j,k,UEINT) + 0.5d0 / uout(i,j,k,URHO) * &
                    ( uout(i,j,k,UMX)**2 + uout(i,j,k,UMY)**2 + uout(i,j,k,UMZ)**2 )
                  
                  ! Make sure the velocities didn't go nutty
                  if (do_diag) then
                      print *,'min / new / max xvel ',min_vel(1), uout(i,j,k,UMX  )/uout(i,j,k,URHO), max_vel(1)
                      print *,'min / new / max yvel ',min_vel(2), uout(i,j,k,UMY  )/uout(i,j,k,URHO), max_vel(2)
                      print *,'min / new / max zvel ',min_vel(3), uout(i,j,k,UMZ  )/uout(i,j,k,URHO), max_vel(3)
                      print *,'min / new / max  e   ',min_e     , uout(i,j,k,UEINT)/uout(i,j,k,URHO), max_e

                      print *,'Adding to xmom ',delta_xmom
                      print *,'Adding to ymom ',delta_ymom
                      print *,'Adding to zmom ',delta_zmom
                      print *,'Adding to re   ',delta_rhoe

                      if (UFS .gt. 0) then
                      !     print *,'new ufs       ',uout(i,j,k,UFS  )/uout(i,j,k,URHO)
                      !     print *,'new ufsp1     ',uout(i,j,k,UFS+1)/uout(i,j,k,URHO)
                      end if

                      sum_after(:) = 0.d0
                      do kk = klo,khi
                      do jj = jlo,jhi
                      do ii = ilo,ihi
                          sum_after(URHO:nmax) = sum_after(URHO:nmax) + uout(ii,jj,kk,URHO:nmax)
                      end do
                      end do
                      end do

                      do n = URHO, nmax
                          print *,"SUMS: BEFORE AFTER ",n,sum_before(n), sum_after(n)
                      end do
                  end if  ! end do_diag

               end if  ! end if rho < rho_min
            enddo
         enddo
      enddo

    end subroutine enforce_minimum_density
