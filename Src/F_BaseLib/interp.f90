module interp_module

  use bl_types
  use bl_error_module
  use bl_constants_module
!  use bc_module

  implicit none

  private :: ix_proj

contains

  elemental function ix_proj(a,b) result(r)
    integer :: r
    integer, intent(in) :: a, b
    r = (a+b*abs(a))/b-abs(a)
  end function ix_proj

  !! bl_nd_interp:  node based bilinear interpolation

  subroutine bl_nd_interp_2d (crse, fine, lratio)

    integer, intent(in) :: lratio(:)
    real(kind=dp_t), intent(out) ::  fine(0:,0:,:)
    real(kind=dp_t), intent(in ) ::  crse(0:,0:,:)
    integer, parameter ::  NUM_SLP = 3
    real(kind=dp_t)  :: sl(0:size(crse,1)-1, NUM_SLP)

    integer :: lx, ly
    integer :: i, j, ifn, jfn, n
    integer :: ilo, ihi, jlo, jhi

    real(kind=dp_t) :: fx, fy
    real(kind=dp_t) :: RX, RY, RXY
    real(kind=dp_t) :: dx0, d0x, dx1

    RX = ONE/real(lratio(1),kind=dp_t)
    RY = ONE/real(lratio(2),kind=dp_t)
    RXY = RX*RY

    do  n = 1, size(crse,3)
       do j = 0, size(crse,2)-1
          jlo = j*lratio(2)
          jhi = min(jlo + lratio(2) - 1, size(fine,2)-1)
          do i = 0, size(crse,1)-1
             dx0 = ZERO
             if (i  /=  size(crse,1)-1 ) dx0 = crse(i+1,j,n) - crse(i,j,n)

             d0x = ZERO
             if (j  /=  size(crse,2)-1 ) d0x = crse(i,j+1,n) - crse(i,j,n)

             dx1 = ZERO
             if (i  /=  size(crse,1)-1 .and. j  /=  size(crse,2)-1 ) dx1 = crse(i+1,j+1,n) - crse(i,j+1,n)

             sl(i,1) = RX*dx0
             sl(i,2) = RY*d0x
             sl(i,3) = RXY*(dx1 - dx0)
          end do

          !  compute fine strip of interpolated data

          do ly = jlo, jhi
             jfn = lratio(2)*j + ly 
             fy = real(ly,kind=dp_t)

             do i = 0, size(crse,1)-1
                ilo = i*lratio(1)
                ihi = min(ilo + lratio(1) - 1, size(fine,1)-1)
                do lx = ilo, ihi
                   ifn = lratio(1)*i + lx
                   fx = real(lx,kind=dp_t)

                   fine(ifn,jfn,n) = crse(i,j,n) +fx*sl(i,1) + fy*sl(i,2) + fx*fy*sl(i,3)

                end do
             end do
          end do
       end do
    end do

  end subroutine bl_nd_interp_2d

  !! lin_cc_interp:   linear conservative interpolation from coarse grid to fine

  subroutine lin_cc_interp_2d (fine, crse, lratio, bc, &
       fvcx, fvcy, cvcx, cvcy, &
       lim_slope, lin_limit)

    integer, intent(in) :: lratio(:)
    integer, intent(in) :: bc(:,:,:)
    logical, intent(in) :: lim_slope, lin_limit
    real(kind=dp_t), intent(out) :: fine(0:,0:,:)
    real(kind=dp_t), intent(inout) :: crse(-1:,-1:,:)
    real(kind=dp_t), intent(in) :: fvcx(0:)
    real(kind=dp_t), intent(in) :: fvcy(0:)
    real(kind=dp_t), intent(in) :: cvcx(0:)
    real(kind=dp_t), intent(in) :: cvcy(0:)


    real(kind=dp_t) ::     uc_xslope(0:ubound(crse,1)-1, 0:ubound(crse,2)-1,size(fine,3))
    real(kind=dp_t) ::     lc_xslope(0:ubound(crse,1)-1, 0:ubound(crse,2)-1,size(fine,3))
    real(kind=dp_t) :: xslope_factor(0:ubound(crse,1)-1, 0:ubound(crse,2)-1)
    real(kind=dp_t) ::     uc_yslope(0:ubound(crse,1)-1, 0:ubound(crse,2)-1,size(fine,3))
    real(kind=dp_t) ::     lc_yslope(0:ubound(crse,1)-1, 0:ubound(crse,2)-1,size(fine,3))
    real(kind=dp_t) :: yslope_factor(0:ubound(crse,1)-1, 0:ubound(crse,2)-1)
    real(kind=dp_t) ::         alpha(0:ubound(crse,1)-1, 0:ubound(crse,2)-1,size(fine,3))
    real(kind=dp_t) ::          cmax(0:ubound(crse,1)-1, 0:ubound(crse,2)-1,size(fine,3))
    real(kind=dp_t) ::          cmin(0:ubound(crse,1)-1, 0:ubound(crse,2)-1,size(fine,3))
    real(kind=dp_t) ::         voffx(0:size(fine,1)-1)
    real(kind=dp_t) ::         voffy(0:size(fine,2)-1)

    integer :: n
    integer :: i, ic
    integer :: j, jc
    real(kind=dp_t) :: fxcen, cxcen, fycen, cycen
    real(kind=dp_t) :: orig_corr_fact, corr_fact
    logical :: xok(2)
    integer :: nxc(2)
    integer :: ioff, joff

    forall (i =1:2) nxc(i) = size(crse,i)-2

    xok = (nxc >=  2)

    do j = 0, size(fine,2)-1
       jc = IX_PROJ(j,lratio(2))
       fycen = HALF*(fvcy(j)+fvcy(j+1))
       cycen = HALF*(cvcy(jc)+cvcy(jc+1))
       voffy(j) = (fycen-cycen)/(cvcy(jc+1)-cvcy(jc))
    end do
    do i = 0, size(fine,1)-1
       ic = IX_PROJ(i,lratio(1))
       fxcen = HALF*(fvcx(i)+fvcx(i+1))
       cxcen = HALF*(cvcx(ic)+cvcx(ic+1))
       voffx(i) = (fxcen-cxcen)/(cvcx(ic+1)-cvcx(ic))
    end do

    ! Prevent underflow for small crse values.
    where ( abs(crse) <= 1.0e-20_dp_t ) crse = ZERO
    alpha = ONE
    cmax = crse(0:nxc(1)-1,0:nxc(2)-1,:)
    cmin = crse(0:nxc(1)-1,0:nxc(2)-1,:)

    do n = 1, size(fine,3)
       ! Initialize alpha = 1 and define cmax and cmin as neighborhood max/mins.
       do j = 0, nxc(2)-1
          do i = 0, nxc(1)-1
             do joff = -1, 1
                do ioff = -1, 1
                   cmax(i,j,n) = max(cmax(i,j,n),crse(i+ioff,j+joff,n))
                   cmin(i,j,n) = min(cmin(i,j,n),crse(i+ioff,j+joff,n))
                end do
             end do
          end do
       end do

    end do
    !
    ! Compute unlimited and limited slopes
    !
    do n = 1, size(fine,3)
       do j = 0, nxc(2)-1
          do i = 0, nxc(1)-1
             uc_xslope(i,j,n) = HALF*(crse(i+1,j,n)-crse(i-1,j,n))
             lc_xslope(i,j,n) = uclc_slope(uc_xslope(i,j,n),crse,i,j,n,dim=1)
          end do
       end do

       if (bc(1,1,n)  ==  EXT_DIR .or. bc(1,1,n) == HOEXTRAP) then
          i = 0
          if ( xok(1) ) then
             do j = 0, nxc(2)-1
                uc_xslope(i,j,n) = -SIXTEEN/FIFTEEN*crse(i-1,j,n)+ HALF*crse(i,j,n) &
                     + TWO3RD*crse(i+1,j,n) - TENTH*crse(i+2,j,n)
             end do
          else
             do j = 0, nxc(2)-1
                uc_xslope(i,j,n) = FOURTH * (crse(i+1,j,n) + FIVE*crse(i,j,n) - SIX*crse(i-1,j,n) )
             end do
          end if
          do j = 0, nxc(2)-1
             lc_xslope(i,j,n) = uclc_slope(uc_xslope(i,j,n),crse,i,j,n,dim=1)
          end do
       end if

       if (bc(1,2,n)  ==  EXT_DIR .or. bc(1,2,n) == HOEXTRAP) then
          i = nxc(1)-1
          if ( xok(1) ) then
             do j = 0, nxc(2)-1
                uc_xslope(i,j,n) = SIXTEEN/FIFTEEN*crse(i+1,j,n)- HALF*crse(i,j,n) &
                     - TWO3RD*crse(i-1,j,n) + TENTH*crse(i-2,j,n)
             end do
          else
             do j = 0, nxc(2)-1
                uc_xslope(i,j,n) = -FOURTH * (crse(i-1,j,n) + FIVE*crse(i,j,n) - SIX*crse(i+1,j,n) )
             end do
          end if
          do j = 0, nxc(2)-1
             lc_xslope(i,j,n) = uclc_slope(uc_xslope(i,j,n),crse,i,j,n, dim=1)
          end do
       end if

       do j = 0, nxc(2)-1
          do i = 0, nxc(1)-1
             uc_yslope(i,j,n) = HALF*(crse(i,j+1,n)-crse(i,j-1,n))
             lc_yslope(i,j,n) = uclc_slope(uc_yslope(i,j,n),crse,i,j,n,dim=2)
          end do
       end do

       if (bc(2,1,n)  == EXT_DIR .or. bc(2,1,n) == HOEXTRAP) then
          j = 0
          if ( xok(1) ) then
             do i = 0, nxc(1)-1
                uc_yslope(i,j,n) = -SIXTEEN/FIFTEEN*crse(i,j-1,n)+ HALF*crse(i,j,n) &
                     + TWO3RD*crse(i,j+1,n) - TENTH*crse(i,j+2,n)
             end do
          else
             do i = 0, nxc(1)-1
                uc_yslope(i,j,n) = FOURTH * (crse(i,j+1,n) + FIVE*crse(i,j,n) - SIX*crse(i,j-1,n) )
             end do
          end if
          do i = 0, nxc(1)-1
             lc_yslope(i,j,n) = uclc_slope(uc_yslope(i,j,n),crse,i,j,n,dim=2)
          end do
       end if

       if (bc(2,2,n)  ==  EXT_DIR .or. bc(2,2,n) == HOEXTRAP) then
          j = nxc(2)-1
          if ( xok(2) ) then
             do i = 0, nxc(1)-1
                uc_yslope(i,j,n) = SIXTEEN/FIFTEEN*crse(i,j+1,n)- HALF*crse(i,j,n) &
                     - TWO3RD*crse(i,j-1,n) + TENTH*crse(i,j-2,n)
             end do
          else
             do i = 0, nxc(1)-1
                uc_yslope(i,j,n) = -FOURTH * (crse(i,j-1,n) + FIVE*crse(i,j,n) - SIX*crse(i,j+1,n) )
             end do
          end if
          do i = 0, nxc(1)-1
             lc_yslope(i,j,n) = uclc_slope(uc_yslope(i,j,n),crse,i,j,n,dim=2)
          end do
       end if

    end do

    if ( lim_slope ) then

       ! Do the interpolation using unlimited slopes.

       do n = 1, size(fine,3)
          do j = 0, size(fine,2)-1
             jc = IX_PROJ(j,lratio(2))
             do i = 0, size(fine,1)-1
                ic = IX_PROJ(i,lratio(1))
                fine(i,j,n) = crse(ic,jc,n) + voffx(i)*uc_xslope(ic,jc,n)+ voffy(j)*uc_yslope(ic,jc,n)
             end do
          end do
       end do

    else 

       if ( lin_limit ) then

          ! compute linear limited slopes
          ! Note that the limited and the unlimited slopes
          ! have the same sign, and it is assumed that they do.

          ! compute slope factors

          xslope_factor = ONE
          yslope_factor = ONE

          do n = 1, size(fine,3)
             where( uc_xslope(:,:,n) /= 0 )
                xslope_factor = min(xslope_factor,lc_xslope(:,:,n)/uc_xslope(:,:,n),ONE)
             end where
             where( uc_yslope(:,:,n) /= 0 )
                yslope_factor = min(yslope_factor,lc_yslope(:,:,n)/uc_yslope(:,:,n),ONE)
             end where
          end do

          ! compute linear limited slopes

          do n = 1, size(fine,3)
             lc_xslope(:,:,n) = xslope_factor*uc_xslope(:,:,n)
             lc_yslope(:,:,n) = yslope_factor*uc_yslope(:,:,n)
          end do

       else

          ! Limit slopes so as to not introduce new maxs or mins.

          do n = 1, size(fine,3)
             do j = 0, size(fine,2)-1

                jc = IX_PROJ(j,lratio(2))

                do i = 0, size(fine,1)-1
                   ic = IX_PROJ(i,lratio(1))

                   orig_corr_fact = voffx(i)*lc_xslope(ic,jc,n)+ voffy(j)*lc_yslope(ic,jc,n) 
                   fine(i,j,n) = crse(ic,jc,n) + orig_corr_fact
                   if ((fine(i,j,n)  >  cmax(ic,jc,n)) &
                        .and.(abs(orig_corr_fact)  >  1.e-10*abs(crse(ic,jc,n)))) then
                      corr_fact = (cmax(ic,jc,n) - crse(ic,jc,n)) / orig_corr_fact
                      alpha(ic,jc,n) = min(alpha(ic,jc,n),corr_fact)
                   end if
                   if ((fine(i,j,n)  <  cmin(ic,jc,n)) &
                        .and.(abs(orig_corr_fact)  >  1.e-10*abs(crse(ic,jc,n)))) then
                      corr_fact = (cmin(ic,jc,n) - crse(ic,jc,n)) / orig_corr_fact
                      alpha(ic,jc,n) = min(alpha(ic,jc,n),corr_fact)
                   end if

                end do
             end do
          end do

       end if

       ! Do the interpolation with limited slopes.

       do n = 1, size(fine,3)
          do j = 0, size(fine,2)-1
             jc = IX_PROJ(j,lratio(2))

             do i = 0, size(fine,1)-1
                ic = IX_PROJ(i,lratio(1))

                fine(i,j,n) = crse(ic,jc,n) + alpha(ic,jc,n)*( &
                     + voffx(i)*lc_xslope(ic,jc,n) &
                     + voffy(j)*lc_yslope(ic,jc,n) &
                     )
             end do
          end do
       end do

    end if

  contains

    function uclc_slope(uslope, crse, i, j, n, dim) result(lc)
      real(kind = dp_t) :: lc
      real(kind=dp_t), intent(in) :: uslope
      real(kind=dp_t), intent(in)  :: crse(-1:,-1:,:)
      real(kind=dp_t) :: cen, forw, back, slp
      integer, intent(in) :: i, j, n, dim
      if ( dim == 1 ) then
         cen  = uslope
         forw = TWO*(crse(i+1,j,n)-crse(i,j,n))
         back = TWO*(crse(i,j,n)-crse(i-1,j,n))
      else if ( dim == 2 ) then
         cen  = uslope
         forw = TWO*(crse(i,j+1,n)-crse(i,j,n))
         back = TWO*(crse(i,j,n)-crse(i,j-1,n))
      end if
      slp  = min(abs(forw),abs(back))
      if (forw*back  <  ZERO) then
         slp = ZERO
      end if
      lc = sign(ONE,cen)*min(slp,abs(cen))
    end function uclc_slope

  end subroutine lin_cc_interp_2d

  !! pcinterp:  cell centered piecewise constant interpolation

  subroutine pc_cc_interp_2d (crse, fine, lratio)
    integer, intent(in) :: lratio(2)
    real(kind=dp_t), intent(in) :: crse(0:,0:,:)
    real(kind=dp_t), intent(out) :: fine(0:,0:,:)

    integer :: i, j, ic, jc, ioff, joff, n

    do n = 1, size(crse,3)
       do jc = 0, ubound(crse,2)
          do ic = 0, ubound(crse,1)
             do joff = 0, lratio(2)-1
                j = lratio(2)*jc + joff
                do ioff = 0, lratio(1)-1
                   i = lratio(1)*ic + ioff
                   fine(i,j,n) = crse(ic,jc,n)
                end do
             end do
          end do
       end do
    end do

  end subroutine pc_cc_interp_2d

  !! protect_interp:   redo interpolation if the result of linccinterp
  !! generates under- or overshoots.

  subroutine pr_cc_interp_2d (fine, crse, lratio, fine_state, &
       fvcx, fvcy, cvcx, cvcy)

    integer, intent(in) :: lratio(:)
    real(kind=dp_t), intent(out) :: fine(0:,0:,0:)
    real(kind=dp_t), intent(inout) :: crse(-1:,-1:,0:) ! FIXME
    real(kind=dp_t), intent(in) :: fine_state(0:,0:,0:)

    real(kind=dp_t), intent(in) :: fvcx(0:)
    real(kind=dp_t), intent(in) :: fvcy(0:)
    real(kind=dp_t), intent(in) :: cvcx(0:)
    real(kind=dp_t), intent(in) :: cvcy(0:)

    real(kind=dp_t) :: alpha, sumN, sumP, negVal, posVal
    real(kind=dp_t) :: crseTot, crseTotnew
    real(kind=dp_t) :: orig_fine(0:lratio(1)-1,0:lratio(2)-1)
    real(kind=dp_t) :: fvol,cvol
    logical :: redo_me
    integer :: ilo, ihi, jlo, jhi
    integer :: i, j, ic, jc, n
    integer :: icase, nxc(2)

    nxc(1) = size(crse,1)-2
    nxc(2) = size(crse,2)-2

    do jc = 0, size(crse,2)-1
       do ic = 0, size(crse,1)-1

          ilo = lratio(1)*ic
          ihi = ilo + lratio(1) - 1
          jlo = lratio(2)*jc
          jhi = jlo + lratio(2) - 1

          do n = 1, ubound(fine,3)

             redo_me = .false.
             do j = jlo,jhi
                do i = ilo,ihi
                   if ((fine_state(i,j,n)+fine(i,j,n))  <  ZERO) redo_me = .true.
                end do
             end do

!!! ****************************************************************************************
!!!
!!! If all the fine values are non-negative after the original interpolated 
!!! correction, then we do nothing here.
!!!
!!! If any of the fine values are negative after the original interpolated
!!! correction, then we do our best.
!!!
!!! Special cases:
!!!
!!! 1) Coarse correction > 0, and fine_state has some cells with 
!!! negative values which will be filled before adding to the other cells.
!!! Use the correction to bring negative cells to ZERO, then
!!! distribute the remaining positive proportionally.
!!!
!!! 2) Coarse correction > 0, and correction can not make them all
!!! positive.  Add correction only to the negative cells, in proportion
!!! to their magnitude.
!!!
!!! 3) Coarse correction < 0, and fine_state DOES NOT have enough
!!! have enough positive state to absorb it.  Here we bring
!!! all the positive fine cells to ZERO then distribute the remaining
!!! negative amount in such a way as to make them all as close to the
!!! same negative value as possible.
!!!
!!! 4) Coarse correction < 0, fine_state has enough
!!! positive state to absorb it without making any fine 
!!! cells negative, BUT fine_state+fine is currently negative
!!! in at least ONE fine cell.  Here just take a constant percentage
!!! away from each positive and don't touch the negatives.
!!!
!!! crseTot = volume-weighted sum of all interpolated values of the correction,
!!! which is equivalent to the total volume-weighted coarse correction
!!! SumN = volume-weighted sum of all negative values of fine_state
!!! SumP = volume-weighted sum of all positive values of fine_state
!!!
!!! ****************************************************************************************


             if ( redo_me ) then

                icase = 0

                orig_fine = fine(ilo:ihi,jlo:jhi,n)

                crseTot = ZERO
                do j = jlo,jhi
                   do i = ilo,ihi
                      fvol = (fvcx(i+1)-fvcx(i)) * (fvcy(j+1)-fvcy(j))
                      crseTot = crseTot + fvol * fine(i,j,n)
                   end do
                end do

                cvol = (cvcx(ic+1)-cvcx(ic)) * (cvcy(jc+1)-cvcy(jc))

                sumN = ZERO
                sumP = ZERO
                do j = jlo,jhi
                   do i = ilo,ihi
                      fvol = (fvcx(i+1)-fvcx(i)) * (fvcy(j+1)-fvcy(j))
                      if (fine_state(i,j,n)  <=  ZERO) then
                         sumN = SumN + fvol * fine_state(i,j,n)
                      else
                         sumP = sumP + fvol * fine_state(i,j,n)
                      end if
                   end do
                end do

                if (crseTot  >  ZERO .and. crseTot  >=  abs(sumN)) then

                   ! Here we want to fill in the negative values first, then add
                   ! the remaining positive proportionally.

                   icase = 1
                   where(fine_state(ilo:ihi,jlo:jhi,n) <= ZERO ) &
                        fine(ilo:ihi,jlo:jhi,n) = - fine_state(ilo:ihi,jlo:jhi,n)

                   if (sumP > ZERO) then

                      alpha = (crseTot - abs(sumN)) / sumP

                      where(fine_state(ilo:ihi,jlo:jhi,n) >= ZERO ) &
                           fine(ilo:ihi,jlo:jhi,n) = alpha*fine_state(ilo:ihi,jlo:jhi,n)

                   else

                      posVal = (crseTot - abs(sumN)) / cvol

                      fine(ilo:ihi,jlo:jhi,n) = fine(ilo:ihi,jlo:jhi,n) + posVal

                   end if

                end if

                if (crseTot  >  ZERO .and. crseTot  <  abs(sumN)) then
                   ! Here we don't have enough positive correction to fill all the
                   ! negative values of state, so we just try to fill them proportionally
                   ! and don't add any correction to the states already positive.

                   icase = 2
                   alpha = crseTot / abs(sumN)


                   where(fine_state(ilo:ihi,jlo:jhi,n) < ZERO )
                      fine(ilo:ihi,jlo:jhi,n) = alpha*abs(fine_state(ilo:ihi,jlo:jhi,n))
                   elsewhere
                      fine(ilo:ihi,jlo:jhi,n)  = ZERO
                   end where

                end if

                if (crseTot  <  ZERO .and. abs(crseTot)  >  sumP) then
                   ! Here we don't have enough positive states to absorb all the
                   ! negative correction, so we want to end up with all the fine
                   ! cells having the same negative value.

                   icase = 3
                   negVal = (sumP + sumN + crseTot)/cvol

                   fine(ilo:ihi,jlo:jhi,n) = negVal - fine_state(ilo:ihi,jlo:jhi,n)

                end if

                if (crseTot  <  ZERO .and. abs(crseTot)  <  sumP.and. (sumP+sumN+crseTot)  >  ZERO) then
                   ! Here we have enough positive states to absorb all the
                   ! negative correction *and* redistribute to make negative cells
                   ! positive. 

                   icase = 4
                   alpha = (crseTot + sumN) / sumP


                   where(fine_state(ilo:ihi,jlo:jhi,n) < ZERO )
                      fine(ilo:ihi,jlo:jhi,n) = -fine_state(ilo:ihi,jlo:jhi,n)
                   elsewhere
                      fine(ilo:ihi,jlo:jhi,n)  = alpha*fine_state(ilo:ihi,jlo:jhi,n)
                   end where

                end if

                if (crseTot  <  ZERO .and. abs(crseTot)  <  sumP.and. (sumP+sumN+crseTot)  <=  ZERO) then
                   ! Here we have enough positive states to absorb all the
                   ! negative correction, but not to fix the states already negative. 
                   ! We bring all the positive states to ZERO, and use whatever 
                   ! remaining positiveness from the states to help the negative states.

                   icase = 5
                   alpha = (crseTot + sumP) / sumN

                   where(fine_state(ilo:ihi,jlo:jhi,n) > ZERO )
                      fine(ilo:ihi,jlo:jhi,n) = -fine_state(ilo:ihi,jlo:jhi,n)
                   elsewhere
                      fine(ilo:ihi,jlo:jhi,n)  = alpha*fine_state(ilo:ihi,jlo:jhi,n)
                   end where

                end if

                crseTotnew   = ZERO
                do j = jlo,jhi
                   do i = ilo,ihi
                      fvol = (fvcx(i+1)-fvcx(i)) * (fvcy(j+1)-fvcy(j))
                      crseTotnew   = crseTotnew   + fvol * fine(i,j,n)
                   end do
                end do

                if (abs(crseTotnew - crseTot)/cvol  >  1.e-8) then
                   print *,' '
                   print *,'BLEW CONSERVATION with ICASE = ',icase
                   print *,'AT COARSE CELL ',ic,jc,' AND COMPONENT ',n
                   print *,'CRSETOT NEW OLD ',crseTotnew, crseTot
                   print *,'CVOL ',cvol
                   print *,'SUMP SUMN ',sumP,sumN
                   do j = jlo,jhi
                      do i = ilo,ihi
                         fvol = (fvcx(i+1)-fvcx(i)) * (fvcy(j+1)-fvcy(j))
                         print *,'FINE OLD NEW ', &
                              i,j,orig_fine(i-ilo,j-jlo),fine(i,j,n), fine_state(i,j,n),fvol
                         if (abs(fvol)  <  1.e-20) then
                            print *,'MAKING FVOL ',fvcx(i+1),fvcx(i),fvcy(j+1),fvcy(j)
                         end if
                      end do
                   end do
                end if
             end if

          end do

          ! Set sync for density (n=1) to sum of spec sync (2:nvar-1)
          fine(:,:,0) = sum(fine(:,:,1:),dim=3)

          ! End of coarse index loops
       end do
    end do
  end subroutine pr_cc_interp_2d


  !! bl_nd_interp:  node based bilinear interpolation

  subroutine bl_nd_interp_3d (crse, fine, lratio)

    integer, intent(in) :: lratio(:)
    integer, parameter :: NUM_SLP = 7
    REAL(kind=dp_t) :: fine(0:,0:,0:,:)
    REAL(kind=dp_t) :: crse(0:,0:,0:,:)
    REAL(kind=dp_t) :: sl(0:ubound(crse,1), NUM_SLP)

    integer :: lx, ly, lz
    integer :: i, j, k, ifn, jfn, kfn, n
    integer :: ilo, ihi, jlo, jhi, klo, khi

    REAL(kind=dp_t) :: fx, fy,fz
    REAL(kind=dp_t) :: RX, RY, RZ, RXY, RXZ, RYZ, RXYZ
    REAL(kind=dp_t) :: dx00, d0x0, d00x, dx10, dx01, d0x1, dx11

    RX = ONE/real(lratio(1),kind=dp_t)
    RY = ONE/real(lratio(2),kind=dp_t)
    RZ = ONE/real(lratio(3),kind=dp_t)
    RXY = RX*RY
    RXZ = RX*RZ
    RYZ = RY*RZ
    RXYZ = RX*RY*RZ
    !
    ! NOTES:
    ! 1) (i, j, k) loop over the coarse cells
    ! 2) ?strtFine and ?stopFine are the beginning and ending fine cell
    ! indices corresponding to the current coarse cell.  ?stopFine
    ! is restricted for the last coarse cell in each direction since
    ! for this cell we only need to do the face and not the fine nodes
    ! inside this cell.
    ! 3) (lx, ly, lz) as well as ?lo and ?hi refer to the fine node indices
    ! as an offset from ?strtFine.
    !
    do n = 1, size(crse,4)
       do  k = 0, ubound(crse,3)
          klo = k * lratio(3)
          khi = min(klo + lratio(3) - 1, ubound(fine,3))
          do  j = 0, ubound(crse,2)
             jlo = j * lratio(2)
             jhi = min(jlo + lratio(2) - 1, ubound(fine,2))
             !
             !  compute slopes 
             !
             ! NOTE: The IF logic in the calculation of the slopes is to
             ! prevent stepping out of bounds on the coarse data when
             ! computing the slopes on the ARG_H?(cb) cells.  These
             ! slopes actually are not used since they are multiplied by
             ! ZERO.
             !
             do i = 0, ubound(crse,1)
                dx00 = ZERO
                if (i  /=  ubound(crse,1) ) dx00 = crse(i+1,j,k,n) - crse(i,j,k,n)

                d0x0 = ZERO
                if (j  /=  ubound(crse,2) ) d0x0 = crse(i,j+1,k,n) - crse(i,j,k,n)

                d00x = ZERO
                if (k  /=  ubound(crse,3) ) d00x = crse(i,j,k+1,n) - crse(i,j,k,n)

                dx10 = ZERO
                if (i  /=  ubound(crse,1) .and. j  /=  ubound(crse,2) ) &
                     dx10 = crse(i+1,j+1,k,n) - crse(i,j+1,k,n)

                dx01 = ZERO
                if (i  /=  ubound(crse,1) .and. k  /=  ubound(crse,3) ) &
                     dx01 = crse(i+1,j,k+1,n) - crse(i,j,k+1,n)

                d0x1 = ZERO
                if (j  /=  ubound(crse,2) .and. k  /=  ubound(crse,3) ) &
                     d0x1 = crse(i,j+1,k+1,n) - crse(i,j,k+1,n)

                dx11 = ZERO
                if (i  /=  ubound(crse,1) .and. j  /=  ubound(crse,2) .and. k  /=  ubound(crse,3) ) &
                    dx11 = crse(i+1,j+1,k+1,n) - crse(i,j+1,k+1,n)

                sl(i,1) = RX*dx00
                sl(i,2) = RY*d0x0
                sl(i,3) = RZ*d00x
                sl(i,4) = RXY*(dx10 - dx00)
                sl(i,5) = RXZ*(dx01 - dx00)
                sl(i,6) = RYZ*(d0x1 - d0x0)
                sl(i,7) = RXYZ*(dx11 - dx01 - dx10 + dx00)
             end do
             !
             ! compute fine strip of interpolated data
             !
             do lz = klo, khi
                kfn = lratio(3) * k + lz
                fz = real(lz, kind=dp_t)
                do ly = jlo, jhi
                   jfn = lratio(2) * j + ly
                   fy = real(ly,kind=dp_t)
                   do i = 0, ubound(crse,1)
                      ilo = i * lratio(1)
                      ihi = min(ilo + lratio(1) - 1, ubound(fine,1))
                      do lx = ilo, ihi
                         ifn = lratio(1) * i + lx
                         fx = real(lx,kind=dp_t)
                         fine(ifn,jfn,kfn,n) = &
                              crse(i,j,k,n) &
                              +fx*sl(i,1) &
                              + fy*sl(i,2) &
                              + fz*sl(i,3) &
                              +fx*fy*sl(i,4) &
                              + fx*fz*sl(i,5) &
                              + fy*fz*sl(i,6) &
                              + fx*fy*fz*sl(i,7)
                      end do
                   end do
                end do
             end do

          end do
       end do
    end do

  end subroutine bl_nd_interp_3d


  !! linccinterp:   linear conservative interpolation from coarse grid to

  subroutine lin_cc_interp_3d (fine,  &
       crse, &
       lratio, &
       bc,fvcx, fvcy, fvcz, cvcx, cvcy, cvcz, lim_slope, lin_limit)

    integer, intent(in) :: bc(:,:,:)
    integer, intent(in) :: lratio(:)
    REAL(kind=dp_t), intent(out) :: fine(0:,0:,0:,:)
    REAL(kind=dp_t), intent(inout) :: crse(-1:,-1:,-1:,:)
    REAL(kind=dp_t), intent(in) :: fvcx(0:)
    REAL(kind=dp_t), intent(in) :: fvcy(0:)
    REAL(kind=dp_t), intent(in) :: fvcz(0:)
    REAL(kind=dp_t), intent(in) :: cvcx(0:)
    REAL(kind=dp_t), intent(in) :: cvcy(0:)
    REAL(kind=dp_t), intent(in) :: cvcz(0:)
    logical, intent(in) ::  lim_slope, lin_limit

    REAL(kind=dp_t) ::     uc_xslope(0:ubound(crse,1)-1, 0:ubound(crse,2)-1, 0:ubound(crse,3)-1,size(crse,4))
    REAL(kind=dp_t) ::     lc_xslope(0:ubound(crse,1)-1, 0:ubound(crse,2)-1, 0:ubound(crse,3)-1,size(crse,4))
    REAL(kind=dp_t) :: xslope_factor(0:ubound(crse,1)-1, 0:ubound(crse,2)-1, 0:ubound(crse,3)-1)
    REAL(kind=dp_t) ::     uc_yslope(0:ubound(crse,1)-1, 0:ubound(crse,2)-1, 0:ubound(crse,3)-1,size(crse,4))
    REAL(kind=dp_t) ::     lc_yslope(0:ubound(crse,1)-1, 0:ubound(crse,2)-1, 0:ubound(crse,3)-1,size(crse,4))
    REAL(kind=dp_t) :: yslope_factor(0:ubound(crse,1)-1, 0:ubound(crse,2)-1, 0:ubound(crse,3)-1)
    REAL(kind=dp_t) ::     uc_zslope(0:ubound(crse,1)-1, 0:ubound(crse,2)-1, 0:ubound(crse,3)-1,size(crse,4))
    REAL(kind=dp_t) ::     lc_zslope(0:ubound(crse,1)-1, 0:ubound(crse,2)-1, 0:ubound(crse,3)-1,size(crse,4))
    REAL(kind=dp_t) :: zslope_factor(0:ubound(crse,1)-1, 0:ubound(crse,2)-1, 0:ubound(crse,3)-1)
    REAL(kind=dp_t) ::         alpha(0:ubound(crse,1)-1, 0:ubound(crse,2)-1, 0:ubound(crse,3)-1,size(crse,4))
    REAL(kind=dp_t) ::          cmax(0:ubound(crse,1)-1, 0:ubound(crse,2)-1, 0:ubound(crse,3)-1,size(crse,4))
    REAL(kind=dp_t) ::          cmin(0:ubound(crse,1)-1, 0:ubound(crse,2)-1, 0:ubound(crse,3)-1,size(crse,4))
    REAL(kind=dp_t) :: voffx(0:size(fine,1)-1)
    REAL(kind=dp_t) :: voffy(0:size(fine,2)-1)
    REAL(kind=dp_t) :: voffz(0:size(fine,3)-1)

    integer :: n 
    integer :: i, ic
    integer :: j, jc
    integer :: k, kc
    REAL(kind=dp_t) :: fxcen, cxcen, fycen, cycen, fzcen, czcen
    REAL(kind=dp_t) :: corr_fact,orig_corr_fact
    logical :: xok(3)
    integer :: nxc(3)
    integer :: ioff, joff, koff

    forall(i = 1:3) nxc(i) = size(crse,dim=i) - 2
    xok = (nxc >=  2)

    do k = 0, size(fine,3)-1
       kc = IX_PROJ(k,lratio(3))
       fzcen = HALF*(fvcz(k)+fvcz(k+1))
       czcen = HALF*(cvcz(kc)+cvcz(kc+1))
       voffz(k) = (fzcen-czcen)/(cvcz(kc+1)-cvcz(kc))
    end do
    do j = 0, size(fine,2)-1
       jc = IX_PROJ(j,lratio(2))
       fycen = HALF*(fvcy(j)+fvcy(j+1))
       cycen = HALF*(cvcy(jc)+cvcy(jc+1))
       voffy(j) = (fycen-cycen)/(cvcy(jc+1)-cvcy(jc))
    end do
    do i = 0, size(fine,1)-1
       ic = IX_PROJ(i,lratio(1))
       fxcen = HALF*(fvcx(i)+fvcx(i+1))
       cxcen = HALF*(cvcx(ic)+cvcx(ic+1))
       voffx(i) = (fxcen-cxcen)/(cvcx(ic+1)-cvcx(ic))
    end do

    ! Prevent underflow for small crse values.
    where(abs(crse) < 1.0d-20) crse = ZERO
    alpha = 1.d0
    cmax = crse(0:nxc(1)-1,0:nxc(2)-1,0:nxc(3)-1,:)
    cmin = crse(0:nxc(1)-1,0:nxc(2)-1,0:nxc(3)-1,:)
    
    do n = 1, size(crse,4) 
       !
       !
       !
       ! Initialize alpha = 1 and define cmax and cmin as neighborhood max/mins.
       !
       do k = 0, nxc(3)-1
          do j = 0, nxc(2) - 1
             do i = 0, nxc(1) - 1
                do koff = -1,1
                   do joff = -1,1
                      do ioff = -1,1
                         cmax(i,j,k,n) = max(cmax(i,j,k,n),crse(i+ioff,j+joff,k+koff,n))
                         cmin(i,j,k,n) = min(cmin(i,j,k,n),crse(i+ioff,j+joff,k+koff,n))
                      end do
                   end do
                end do
             end do
          end do
       end do

    end do
    !
    ! Computed unlimited and limited slopes
    !
    do n = 1, size(crse,4) 

       do k=0, nxc(3)-1
          do j=0, nxc(2)-1
             do i=0, nxc(1)-1
                uc_xslope(i,j,k,n) = HALF*(crse(i+1,j,k,n)-crse(i-1,j,k,n))
                lc_xslope(i,j,k,n) = uclc_slope(uc_xslope(i,j,k,n),crse,i,j,k,n,dim=1)
             end do
          end do
       end do

       if (bc(1,1,n)  ==  EXT_DIR .or. bc(1,1,n) == HOEXTRAP) then
          i = 0
          if ( xok(1) ) then
             do k=0, nxc(3)-1
                do j=0, nxc(2) - 1
                   uc_xslope(i,j,k,n) = -SIXTEEN/FIFTEEN*crse(i-1,j,k,n) &
                        + HALF*crse(i,j,k,n)+ TWO3RD*crse(i+1,j,k,n) - TENTH*crse(i+2,j,k,n)
                end do
             end do
          else
             do k=0, nxc(3) - 1
                do j=0, nxc(2) - 1
                   uc_xslope(i,j,k,n) = FOURTH * (crse(i+1,j,k,n) + FIVE*crse(i,j,k,n) - SIX*crse(i-1,j,k,n) )
                end do
             end do
          endif
          do k=0, nxc(3) - 1
             do j=0, nxc(2) - 1
                lc_xslope(i,j,k,n) = uclc_slope(uc_xslope(i,j,k,n),crse,i,j,k,n,dim=1)
             end do
          end do
       end if

       if (bc(1,2,n)  ==  EXT_DIR .or. bc(1,2,n) == HOEXTRAP) then
          i = nxc(1) - 1
          if ( xok(1) ) then
             do k=0, nxc(3)-1
                do j=0, nxc(2)-1
                   uc_xslope(i,j,k,n) = SIXTEEN/FIFTEEN*crse(i+1,j,k,n) &
                        - HALF*crse(i,j,k,n)- TWO3RD*crse(i-1,j,k,n) + TENTH*crse(i-2,j,k,n)
                end do
             end do
          else 
             do k=0, nxc(3)-1
                do j=0, nxc(2)-1
                   uc_xslope(i,j,k,n) = -FOURTH * (crse(i-1,j,k,n) &
                        + FIVE*crse(i,j,k,n) - SIX*crse(i+1,j,k,n) )
                end do
             end do
          endif
          do k=0, nxc(3)-1
             do j=0, nxc(2)-1
                lc_xslope(i,j,k,n) = uclc_slope(uc_xslope(i,j,k,n),crse,i,j,k,n,dim=1)
             end do
          end do
       end if

       do k=0, nxc(3)-1
          do j=0, nxc(2)-1
             do i=0, nxc(1)-1
                uc_yslope(i,j,k,n) = HALF*(crse(i,j+1,k,n)-crse(i,j-1,k,n))
                lc_yslope(i,j,k,n) = uclc_slope(uc_yslope(i,j,k,n),crse,i,j,k,n,dim=2)
             end do
          end do
       end do

       if (bc(2,1,n)  ==  EXT_DIR .or. bc(2,1,n) == HOEXTRAP) then
          j = 0
          if ( xok(2) ) then
             do k=0, nxc(3)-1
                do i=0, nxc(1)-1
                   uc_yslope(i,j,k,n) = -SIXTEEN/FIFTEEN*crse(i,j-1,k,n) &
                        + HALF*crse(i,j,k,n)+ TWO3RD*crse(i,j+1,k,n) - TENTH*crse(i,j+2,k,n)
                end do
             end do
          else
             do k=0, nxc(3)-1
                do i=0, nxc(1)-1
                   uc_yslope(i,j,k,n) = FOURTH * (crse(i,j+1,k,n) + FIVE*crse(i,j,k,n) - SIX*crse(i,j-1,k,n) )
                end do
             end do
          endif
          do k=0, nxc(3)-1
             do i=0, nxc(1)-1
                lc_yslope(i,j,k,n) = uclc_slope(uc_yslope(i,j,k,n),crse,i,j,k,n,dim=2)
             end do
          end do
       end if

       if (bc(2,2,n)  ==  EXT_DIR .or. bc(2,2,n) == HOEXTRAP) then
          j = nxc(2)-1
          if ( xok(2) ) then
             do k=0, nxc(3)-1
                do i=0, nxc(1)-1
                   uc_yslope(i,j,k,n) = SIXTEEN/FIFTEEN*crse(i,j+1,k,n) &
                        - HALF*crse(i,j,k,n)- TWO3RD*crse(i,j-1,k,n) + TENTH*crse(i,j-2,k,n)
                end do
             end do
          else
             do k=0, nxc(3)-1
                do i=0, nxc(1)-1
                   uc_yslope(i,j,k,n) = -FOURTH * (crse(i,j-1,k,n) &
                        + FIVE*crse(i,j,k,n) - SIX*crse(i,j+1,k,n) )
                end do
             end do
          endif
          do k=0, nxc(3)-1
             do i=0, nxc(1)-1
                lc_yslope(i,j,k,n) = uclc_slope(uc_yslope(i,j,k,n),crse,i,j,k,n,dim=2)
             end do
          end do
       end if

       do k=0, nxc(3)-1
          do j=0, nxc(2)-1
             do i=0, nxc(1)-1
                uc_zslope(i,j,k,n) = HALF*(crse(i,j,k+1,n)-crse(i,j,k-1,n))
                lc_zslope(i,j,k,n) = uclc_slope(uc_zslope(i,j,k,n),crse,i,j,k,n,dim=3)
             end do
          end do
       end do

       if (bc(3,1,n)  == EXT_DIR .or. bc(3,1,n) == HOEXTRAP ) then
          k = 0
          if ( xok(3) ) then
             do j=0, nxc(2)-1
                do i=0, nxc(1)-1
                   uc_zslope(i,j,k,n) = -SIXTEEN/FIFTEEN*crse(i,j,k-1,n) &
                        + HALF*crse(i,j,k,n)+ TWO3RD*crse(i,j,k+1,n) - TENTH*crse(i,j,k+2,n)
                end do
             end do
          else
             do j=0, nxc(2)-1
                do i=0, nxc(1)-1
                   uc_zslope(i,j,k,n) = FOURTH * (crse(i,j,k+1,n) + FIVE*crse(i,j,k,n) - SIX*crse(i,j,k-1,n) )
                end do
             end do
          endif
          do j=0, nxc(2)-1
             do i=0, nxc(1)-1
                lc_zslope(i,j,k,n) = uclc_slope(uc_zslope(i,j,k,n),crse,i,j,k,n,dim=3)
             end do
          end do
       end if

       if (bc(3,2,n)  ==  EXT_DIR .or. bc(3,2,n)  ==  HOEXTRAP) then
          k = nxc(3)-1
          if ( xok(3) ) then
             do j=0, nxc(2)-1
                do i=0, nxc(1)-1
                   uc_zslope(i,j,k,n) = SIXTEEN/FIFTEEN*crse(i,j,k+1,n) &
                        - HALF*crse(i,j,k,n)- TWO3RD*crse(i,j,k-1,n) + TENTH*crse(i,j,k-2,n)
                end do
             end do
          else
             do j=0, nxc(2)-1
                do i=0, nxc(1)-1
                   uc_zslope(i,j,k,n) = -FOURTH * (crse(i,j,k-1,n) &
                        + FIVE*crse(i,j,k,n) - SIX*crse(i,j,k+1,n) )
                end do
             end do
          endif
          do j=0, nxc(2)-1
             do i=0, nxc(1)-1
                lc_zslope(i,j,k,n) = uclc_slope(uc_zslope(i,j,k,n),crse,i,j,k,n,dim=3)
             end do
          end do
       end if

    end do

    if ( lim_slope ) then

       ! Do the interpolation using unlimited slopes.

       do n = 1, size(crse,4)
          do k = 0, size(fine,3) - 1
             kc = IX_PROJ(k,lratio(3))
             do j = 0, size(fine, 2) - 1
                jc = IX_PROJ(j,lratio(2))
                do i = 0, size(fine,1) - 1
                   ic = IX_PROJ(i,lratio(1))
                   fine(i,j,k,n) = crse(ic,jc,kc,n) &
                        + voffx(i)*uc_xslope(ic,jc,kc,n) &
                        + voffy(j)*uc_yslope(ic,jc,kc,n) &
                        + voffz(k)*uc_zslope(ic,jc,kc,n)
                end do
             end do
          end do
       end do

    else

       if ( lin_limit ) then

          ! compute linear limited slopes
          ! Note that the limited and the unlimited slopes
          ! have the same sign, and it is assumed that they do.
          !
          ! compute slope factors

          xslope_factor = ONE
          yslope_factor = ONE
          zslope_factor = ONE

          do n = 1, size(crse,4) 
             where( uc_xslope(:,:,:,n) /= 0 )
                xslope_factor = min(xslope_factor,lc_xslope(:,:,:,n)/uc_xslope(:,:,:,n),ONE)
             end where
             where( uc_yslope(:,:,:,n) /= 0 )
                yslope_factor = min(yslope_factor,lc_yslope(:,:,:,n)/uc_yslope(:,:,:,n),ONE)
             end where
             where( uc_zslope(:,:,:,n) /= 0 )
                zslope_factor = min(zslope_factor,lc_zslope(:,:,:,n)/uc_zslope(:,:,:,n),ONE)
             end where
          end do

          ! Compute linear limited slopes

          do n = 1, size(crse,4) 
             lc_xslope(:,:,:,n) = xslope_factor*uc_xslope(:,:,:,n)
             lc_yslope(:,:,:,n) = yslope_factor*uc_yslope(:,:,:,n)
             lc_zslope(:,:,:,n) = zslope_factor*uc_zslope(:,:,:,n)
          end do

       else

          ! Limit slopes so as to not introduce new maxs or mins.

          do n = 1,size(crse,4)

             do k = 0, size(fine, 3) - 1
                kc = IX_PROJ(k,lratio(3))
                do j = 0, size(fine, 2) - 1
                   jc = IX_PROJ(j,lratio(2))
                   do i = 0, size(fine, 1) - 1
                      ic = IX_PROJ(i,lratio(1))
                      orig_corr_fact = &
                           + voffx(i)*lc_xslope(ic,jc,kc,n) &
                           + voffy(j)*lc_yslope(ic,jc,kc,n) &
                           + voffz(k)*lc_zslope(ic,jc,kc,n)
                      fine(i,j,k,n) = crse(ic,jc,kc,n) + orig_corr_fact
                      if ((fine(i,j,k,n)  >  cmax(ic,jc,kc,n)) &
                           .and. (abs(orig_corr_fact) > 1.e-10*abs(crse(ic,jc,kc,n)))) then
                         corr_fact = (cmax(ic,jc,kc,n) - crse(ic,jc,kc,n)) / orig_corr_fact
                         alpha(ic,jc,kc,n) = min(alpha(ic,jc,kc,n),corr_fact)
                      endif
                      if ((fine(i,j,k,n)  <  cmin(ic,jc,kc,n)) &
                           .and.(abs(orig_corr_fact)  >  1.e-10*abs(crse(ic,jc,kc,n)))) then
                         corr_fact = (cmin(ic,jc,kc,n) - crse(ic,jc,kc,n)) / orig_corr_fact
                         alpha(ic,jc,kc,n) = min(alpha(ic,jc,kc,n),corr_fact)
                      endif

                   end do
                end do
             end do
          end do
       end if

       ! Do the interpolation with limited slopes.

       do n = 1, size(crse,4)
          do k = 0, size(fine, 3) - 1
             kc = IX_PROJ(k,lratio(3))
             do j = 0, size(fine, 2) - 1
                jc = IX_PROJ(j,lratio(2))
                do i = 0, size(fine, 1) - 1
                   ic = IX_PROJ(i,lratio(1))
                   fine(i,j,k,n) = crse(ic,jc,kc,n) &
                        + alpha(ic,jc,kc,n) *( &
                        + voffx(i)*lc_xslope(ic,jc,kc,n)&
                        + voffy(j)*lc_yslope(ic,jc,kc,n)&
                        + voffz(k)*lc_zslope(ic,jc,kc,n) &
                        )
                end do
             end do
          end do
       end do

    end if
  contains

    function uclc_slope(uslope, crse, i, j, k, n, dim) result(lc)
      real(kind = dp_t) :: lc
      real(kind=dp_t), intent(in) :: uslope
      real(kind=dp_t), intent(in)  :: crse(-1:,-1:,-1:,:)
      real(kind=dp_t) :: cen, forw, back, slp
      integer, intent(in) :: i, j, k, n, dim
      if ( dim == 1 ) then
         cen  = uslope
         forw = TWO*(crse(i+1,j,k,n)-crse(i,j,k,n))
         back = TWO*(crse(i,j,k,n)-crse(i-1,j,k,n))
      else if ( dim == 2 ) then
         cen  = uslope
         forw = TWO*(crse(i,j+1,k,n)-crse(i,j,k,n))
         back = TWO*(crse(i,j,k,n)-crse(i,j-1,k,n))
      else if ( dim == 3 ) then
         cen  = uslope
         forw = TWO*(crse(i,j,k+1,n)-crse(i,j,k,n))
         back = TWO*(crse(i,j,k,n)-crse(i,j,k-1,n))
      end if
      slp  = min(abs(forw),abs(back))
      if ( forw*back  <  ZERO ) then
         slp = ZERO
      end if
      lc = sign(ONE,cen)*min(slp,abs(cen))
    end function uclc_slope

  end subroutine lin_cc_interp_3d

  !! pc_cc_interp_3d

  subroutine pc_cc_interp_3d (crse, fine, lratio)
    integer lratio(:)
    REAL(kind=dp_t), intent(in) :: crse(0:,0:,0:,:)
    REAL(kind=dp_t), intent(out) ::  fine(0:,0:,0:,:)

    integer i, j, k, ic, jc, kc, ioff, joff, koff, n

    do n = 1, size(crse,4)
       do kc = 0, ubound(crse,4)
          do jc = 0, ubound(crse,2)
             do ic = 0, ubound(crse,1)
                do koff = 0, lratio(3) - 1
                   k = lratio(3)*kc + koff
                   do joff = 0, lratio(2) - 1
                      j = lratio(2)*jc + joff
                      do ioff = 0, lratio(1)-1
                         i = lratio(1)*ic + ioff
                         fine(i,j,k,n) = crse(ic,jc,kc,n)
                      end do
                   end do
                end do
             end do
          end do
       end do
    end do

  end subroutine pc_cc_interp_3d

  !!  protect_interp:   redo interpolation if the result of linccinterp

  subroutine pr_cc_interp_3d (fine, crse, lratio, fine_state)

    integer, intent(in) :: lratio(:)
    REAL(kind=dp_t) :: fine(0:,0:,0:,0:)
    REAL(kind=dp_t) :: crse(-1:,-1:,-1:,0:)
    REAL(kind=dp_t) :: fine_state(0:,0:,0:,0:)

    REAL(kind=dp_t) :: alpha, sumN, sumP, crseTot, negVal, posVal
    REAL(kind=dp_t) :: sum_fine_new, sum_fine_old
    REAL(kind=dp_t) :: orig_fine(0:lratio(1)-1,0:lratio(2)-1,0:lratio(3)-1)
    logical :: redo_me
    integer :: ilo, ihi, jlo, jhi, klo, khi
    integer :: i, j, k, ic, jc, kc, n
    integer :: numFineCells
    integer :: icase, nxc(3)

    forall(i=1:3) nxc(i) = size(crse,dim=i)-2

    do kc = 0, nxc(3) - 1
       do jc = 0, nxc(2) - 1
          do ic = 0, nxc(1) - 1

             ilo = lratio(1)*ic
             ihi = ilo + lratio(1) - 1 
             jlo = lratio(2)*jc 
             jhi = jlo + lratio(2) - 1
             klo = lratio(3)*kc
             khi = klo + lratio(3) - 1

             do n = 1, ubound(fine,4)

                redo_me = any(fine_state(ilo:ihi,jlo:jhi,klo:khi,n) < -fine(ilo:ihi,jlo:jhi,klo:khi,n))

                ! **************************************************************************
                !
                ! If all the fine values are non-negative after the original interpolated 
                ! correction, then we do nothing here.
                !
                ! If any of the fine values are negative after the original interpolated
                ! correction, then we do our best.
                !
                ! Special cases:
                !
                ! 1) Coarse correction > 0, and fine_state has some cells with 
                ! negative values which will be filled before adding to the other cells.
                ! Use the correction to bring negative cells to ZERO, then
                ! distribute the remaining positive proportionally.
                !
                ! 2) Coarse correction > 0, and correction can not make them all
                ! positive.  Add correction only to the negative cells, in proportion
                ! to their magnitude.
                !
                ! 3) Coarse correction < 0, and fine_state DOES NOT have enough
                ! have enough positive state to absorb it.  Here we bring
                ! all the positive fine cells to ZERO then distribute the remaining
                ! negative amount in such a way as to make them all as close to the
                ! same negative value as possible.
                !
                ! 4) Coarse correction < 0, fine_state has enough
                ! positive state to absorb it without making any fine 
                ! cells negative, BUT fine_state+fine is currently negative
                ! in at least ONE fine cell.  Here just take a constant percentage
                ! away from each positive and don't touch the negatives.
                !
                ! crseTot = sum of all interpolated values of the correction,
                ! which is equivalent to the coarse correction * ratio**3
                ! SumN = sum of all negative values of fine_state
                ! SumP = sum of all positive values of fine_state
                !
                ! ****************************************************************************************
                !

                if ( redo_me ) then

                   icase = 0
                   sum_fine_old = sum(fine(ilo:ihi,jlo:jhi,klo:khi,n))
                   orig_fine = fine(ilo:ihi,jlo:jhi,klo:khi,n)

                   crseTot = sum_fine_old
                   numFineCells = (ihi-ilo+1) * (jhi-jlo+1) * (khi-klo+1)

                   sumN = sum(fine_state(ilo:ihi,jlo:jhi,klo:khi,n), &
                        mask = fine_state(ilo:ihi,jlo:jhi,klo:khi,n)<=ZERO)
                   sumP = sum(fine_state(ilo:ihi,jlo:jhi,klo:khi,n), &
                        mask = fine_state(ilo:ihi,jlo:jhi,klo:khi,n) >ZERO)

                   if (crseTot  >  ZERO .and. crseTot  >=  abs(sumN)) then
                      ! Here we want to fill in the negative values first, then add
                      ! the remaining positive proportionally.

                      icase = 1
                      where(fine_state(ilo:ihi,jlo:jhi,klo:khi,n)<=ZERO) &
                           fine(ilo:ihi,jlo:jhi,klo:khi,n) = -fine_state(ilo:ihi,jlo:jhi,klo:khi,n)

                      if (sumP > ZERO) then

                         alpha = (crseTot - abs(sumN)) / sumP
                         
                         where(fine_state(ilo:ihi,jlo:jhi,klo:khi,n)>=ZERO) &
                           fine(ilo:ihi,jlo:jhi,klo:khi,n) = alpha*fine_state(ilo:ihi,jlo:jhi,klo:khi,n)

                      else

                         posVal = (crseTot - abs(sumN)) / float(numFineCells)

                         fine(ilo:ihi,jlo:jhi,klo:khi,n) = fine(ilo:ihi,jlo:jhi,klo:khi,n) + posVal

                      endif

                   endif

                   if (crseTot  >  ZERO .and. crseTot  <  abs(sumN)) then

                      ! Here we don't have enough positive correction to fill all the
                      ! negative values of state, so we just try to fill them proportionally
                      ! and don't add any correction to the states already positive.

                      icase = 2
                      alpha = crseTot / abs(sumN)

                      where(fine_state(ilo:ihi,jlo:jhi,klo:khi,n) < ZERO)
                         fine(ilo:ihi,jlo:jhi,klo:khi,n) = alpha*abs(fine_state(ilo:ihi,jlo:jhi,klo:khi,n))
                      elsewhere
                         fine(ilo:ihi,jlo:jhi,klo:khi,n) = ZERO
                      end where

                   endif

                   if (crseTot  <  ZERO .and. abs(crseTot)  >  sumP) then
                      ! Here we don't have enough positive states to absorb all the
                      ! negative correction, so we want to end up with all the fine
                      ! cells having the same negative value.

                      icase = 3
                      negVal = (sumP + sumN + crseTot)/float(numFineCells)

                      fine(ilo:ihi,jlo:jhi,klo:khi,n) = negVal - fine_state(ilo:ihi,jlo:jhi,klo:khi,n)

                   endif

                   if (crseTot  <  ZERO .and. abs(crseTot)  <  sumP.and. (sumP+sumN+crseTot)  >  ZERO) then
                      ! Here we have enough positive states to absorb all the
                      ! negative correction *and* redistribute to make negative cells
                      ! positive. 

                      icase = 4
                      alpha = (crseTot + sumN) / sumP

                      where(fine_state(ilo:ihi,jlo:jhi,klo:khi,n) < ZERO)
                         fine(ilo:ihi,jlo:jhi,klo:khi,n) = -fine_state(ilo:ihi,jlo:jhi,klo:khi,n)
                      elsewhere
                         fine(ilo:ihi,jlo:jhi,klo:khi,n) = alpha*fine_state(ilo:ihi,jlo:jhi,klo:khi,n)
                      end where

                   endif

                   if (crseTot  <  ZERO .and. abs(crseTot)  <  sumP.and. (sumP+sumN+crseTot)  <= ZERO) then
                      ! Here we have enough positive states to absorb all the
                      ! negative correction, but not to fix the states already negative. 
                      ! We bring all the positive states to ZERO, and use whatever 
                      ! remaining positiveness from the states to help the negative states.

                      icase = 5
                      alpha = (crseTot + sumP) / sumN

                      where(fine_state(ilo:ihi,jlo:jhi,klo:khi,n) > ZERO)
                         fine(ilo:ihi,jlo:jhi,klo:khi,n) = -fine_state(ilo:ihi,jlo:jhi,klo:khi,n)
                      elsewhere
                         fine(ilo:ihi,jlo:jhi,klo:khi,n) = alpha*fine_state(ilo:ihi,jlo:jhi,klo:khi,n)
                      end where

                   endif

                   sum_fine_new = sum(fine(ilo:ihi,jlo:jhi,klo:khi,n))

                   if (abs(sum_fine_new - sum_fine_old)  >  1.e-8) then
                      print *,' '
                      print *,'PROTECT_INTERP: BLEW CONSERVATION with ICASE = ',icase
                      print *,'AT COARSE CELL ',ic,jc,kc,' AND COMPONENT ',n
                      print *,'NEW SUM / OLD SUM ',sum_fine_new, sum_fine_old
                      print *,'CRSETOT ',crseTot
                      print *,'SUMP SUMN ',sumP,sumN
                      do k = klo,khi
                         do j = jlo,jhi
                            do i = ilo,ihi
                               print *,'FINE OLD NEW ', &
                                    i, j, k, orig_fine(i-ilo,j-jlo,k-klo), fine(i,j,k,n), fine_state(i,j,k,n)
                            end do
                         end do
                      end do
                   endif

                endif

             end do

             ! Set sync for density (n=1) to sum of spec sync (2:nvar-1)
             fine(:,:,:,0) = sum(fine(:,:,:,1:),4)

          end do
       end do
    end do
  end subroutine pr_cc_interp_3d

end module interp_module
