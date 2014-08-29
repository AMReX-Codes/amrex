module interp_module

  use bl_types
  use bl_constants_module
  use bc_module

  implicit none

  private :: ix_proj

contains

  elemental function ix_proj(a,b) result(r)
    integer :: r
    integer, intent(in) :: a, b
    r = (a+b*abs(a))/b-abs(a)
  end function ix_proj

  !! lin_cc_interp:   linear conservative interpolation from coarse grid to fine

  subroutine lin_cc_interp_1d (fine, fine_lo, crse, crse_lo, &
                               lratio, bc, &
                               fvcx, fvcx_lo, cvcx, cvcx_lo, &
                               cslope_lo, cslope_hi, lim_slope, lin_limit)

    integer, intent(in) ::   fine_lo(:),  crse_lo(:)
    integer, intent(in) :: cslope_lo(:),cslope_hi(:)
    integer, intent(in) :: lratio(:)
    integer, intent(in) :: bc(:,:,:)
    integer, intent(in) :: fvcx_lo,cvcx_lo
    logical, intent(in) :: lim_slope, lin_limit
    real(kind=dp_t), intent(  out) :: fine(fine_lo(1):,:)
    real(kind=dp_t), intent(inout) :: crse(crse_lo(1):,:)
    real(kind=dp_t), intent(in   ) :: fvcx(fvcx_lo:)
    real(kind=dp_t), intent(in   ) :: cvcx(cvcx_lo:)


    real(kind=dp_t) ::     uc_xslope(cslope_lo(1):cslope_hi(1),size(fine,2))
    real(kind=dp_t) ::     lc_xslope(cslope_lo(1):cslope_hi(1),size(fine,2))
    real(kind=dp_t) :: xslope_factor(cslope_lo(1):cslope_hi(1))
    real(kind=dp_t) ::         alpha(cslope_lo(1):cslope_hi(1),size(fine,2))
    real(kind=dp_t) ::          cmax(cslope_lo(1):cslope_hi(1),size(fine,2))
    real(kind=dp_t) ::          cmin(cslope_lo(1):cslope_hi(1),size(fine,2))
    real(kind=dp_t) ::         voffx(fvcx_lo:fvcx_lo+size(fine,1)-1)

    integer :: n
    integer :: i, ic
    real(kind=dp_t) :: fxcen, cxcen
    real(kind=dp_t) :: orig_corr_fact, corr_fact
    logical :: xok(1)
    integer :: nxc(1)
    integer :: ioff

    nxc(1) = size(crse,1)-2

    xok = (nxc >=  2)

    do i = fvcx_lo, fvcx_lo+size(fine,1)-1
       ic = IX_PROJ(i,lratio(1))
       fxcen = HALF*(fvcx(i)+fvcx(i+1))
       cxcen = HALF*(cvcx(ic)+cvcx(ic+1))
       voffx(i) = (fxcen-cxcen)/(cvcx(ic+1)-cvcx(ic))
    end do

    ! Prevent underflow for small crse values.
    where ( abs(crse) <= 1.0e-20_dp_t ) crse = ZERO
    alpha = ONE
    cmax = crse(cslope_lo(1):cslope_hi(1),:)
    cmin = crse(cslope_lo(1):cslope_hi(1),:)

    do n = 1, size(fine,2)
       !
       ! Initialize alpha = 1 and define cmax and cmin as neighborhood max/mins.
       !
       do i = cslope_lo(1), cslope_hi(1)
          do ioff = -1, 1
             cmax(i,n) = max(cmax(i,n),crse(i+ioff,n))
             cmin(i,n) = min(cmin(i,n),crse(i+ioff,n))
          end do
       end do

    end do
    !
    ! Compute unlimited and limited slopes
    !
    do n = 1, size(fine,2)
       do i = cslope_lo(1), cslope_hi(1)
          uc_xslope(i,n) = HALF*(crse(i+1,n)-crse(i-1,n))
          lc_xslope(i,n) = uclc_slope_1d(uc_xslope(i,n),crse,crse_lo,i,n)
       end do

       if ( bc(1,1,n)  ==  EXT_DIR .or. bc(1,1,n) == HOEXTRAP ) then
          i = cslope_lo(1)
          if ( xok(1) ) then
             uc_xslope(i,n) = -SIXTEEN/FIFTEEN*crse(i-1,n)+ HALF*crse(i,n) &
                  + TWO3RD*crse(i+1,n) - TENTH*crse(i+2,n)
          else
             uc_xslope(i,n) = &
                  FOURTH * (crse(i+1,n) + FIVE*crse(i,n) - SIX*crse(i-1,n) )
          end if
          lc_xslope(i,n) = uclc_slope_1d(uc_xslope(i,n),crse,crse_lo,i,n)
       end if

       if ( bc(1,2,n)  ==  EXT_DIR .or. bc(1,2,n) == HOEXTRAP ) then
          i = cslope_hi(1)
          if ( xok(1) ) then
             uc_xslope(i,n) = SIXTEEN/FIFTEEN*crse(i+1,n)- HALF*crse(i,n) &
                  - TWO3RD*crse(i-1,n) + TENTH*crse(i-2,n)
          else
             uc_xslope(i,n) = &
                  -FOURTH * (crse(i-1,n) + FIVE*crse(i,n) - SIX*crse(i+1,n) )
          end if
          lc_xslope(i,n) = uclc_slope_1d(uc_xslope(i,n),crse,crse_lo,i,n)
       end if
    end do

    if ( lim_slope ) then

       if ( lin_limit ) then

          ! compute linear limited slopes
          ! Note that the limited and the unlimited slopes
          ! have the same sign, and it is assumed that they do.

          ! compute slope factors

          xslope_factor = ONE

          do n = 1, size(fine,2)
             where( uc_xslope(:,n) /= 0 )
                xslope_factor = min(xslope_factor,lc_xslope(:,n)/uc_xslope(:,n),ONE)
             end where
          end do

          ! compute linear limited slopes

          do n = 1, size(fine,2)
             lc_xslope(:,n) = xslope_factor*uc_xslope(:,n)
          end do

       else

          ! Limit slopes so as to not introduce new maxs or mins.

          do n = 1, size(fine,2)
             do i = fine_lo(1), fine_lo(1)+size(fine,1)-1
                ic = IX_PROJ(i,lratio(1))

                orig_corr_fact = voffx(i)*lc_xslope(ic,n)
                fine(i,n) = crse(ic,n) + orig_corr_fact
                if ( fine(i,n)  >  cmax(ic,n) ) then
                   if ( abs(orig_corr_fact)  >  1.e-10_dp_t*abs(crse(ic,n)) ) then
                      corr_fact = (cmax(ic,n) - crse(ic,n)) / orig_corr_fact
                      alpha(ic,n) = min(alpha(ic,n),corr_fact)
                   end if
                end if
                if ( fine(i,n)  <  cmin(ic,n) ) then
                   if ( abs(orig_corr_fact)  >  1.e-10_dp_t*abs(crse(ic,n)) ) then
                      corr_fact = (cmin(ic,n) - crse(ic,n)) / orig_corr_fact
                      alpha(ic,n) = min(alpha(ic,n),corr_fact)
                      end if
                end if

             end do
          end do

       end if

       ! Do the interpolation with limited slopes.

       do n = 1, size(fine,2)
          do i = fine_lo(1), fine_lo(1)+size(fine,1)-1
             ic = IX_PROJ(i,lratio(1))

             fine(i,n) = crse(ic,n) + alpha(ic,n)*( &
                  + voffx(i)*lc_xslope(ic,n) )
          end do
       end do

    else

       ! Do the interpolation using unlimited slopes.

       do n = 1, size(fine,2)
          do i = fine_lo(1), fine_lo(1)+size(fine,1)-1
             ic = IX_PROJ(i,lratio(1))
             fine(i,n) = &
                  crse(ic,n) + voffx(i)*uc_xslope(ic,n)
          end do
       end do
    end if

  contains

    function uclc_slope_1d(uslope, crse, crse_lo, i, n) result(lc)

      integer, intent(in) :: crse_lo(:)
      real(kind = dp_t) :: lc
      real(kind=dp_t), intent(in) :: uslope
      real(kind=dp_t), intent(in)  :: crse(crse_lo(1):,:)
      real(kind=dp_t) :: cen, forw, back, slp
      integer, intent(in) :: i,  n
      cen  = uslope
      forw = TWO*(crse(i+1,n)-crse(i,n))
      back = TWO*(crse(i,n)-crse(i-1,n))
      slp  = min(abs(forw),abs(back))
      if ( forw*back  <  ZERO ) then
         slp = ZERO
      end if
      lc = sign(ONE,cen)*min(slp,abs(cen))
    end function uclc_slope_1d

  end subroutine lin_cc_interp_1d

  !! lin_cc_interp:   linear conservative interpolation from coarse grid to fine

  subroutine lin_cc_interp_2d (fine, fine_lo, crse, crse_lo, &
                               lratio, bc, &
                               fvcx, fvcx_lo, fvcy, fvcy_lo, cvcx, cvcx_lo, cvcy, cvcy_lo, &
                               cslope_lo, cslope_hi, lim_slope, lin_limit)

    integer, intent(in) :: fine_lo(:),crse_lo(:)
    integer, intent(in) :: cslope_lo(:),cslope_hi(:)
    integer, intent(in) :: lratio(:)
    integer, intent(in) :: bc(:,:,:)
    integer, intent(in) :: fvcx_lo,fvcy_lo,cvcx_lo,cvcy_lo
    logical, intent(in) :: lim_slope, lin_limit
    real(kind=dp_t), intent(  out) :: fine(fine_lo(1):,fine_lo(2):,:)
    real(kind=dp_t), intent(inout) :: crse(crse_lo(1):,crse_lo(2):,:)
    real(kind=dp_t), intent(in   ) :: fvcx(fvcx_lo:)
    real(kind=dp_t), intent(in   ) :: fvcy(fvcy_lo:)
    real(kind=dp_t), intent(in   ) :: cvcx(cvcx_lo:)
    real(kind=dp_t), intent(in   ) :: cvcy(cvcy_lo:)


    real(kind=dp_t) ::     uc_xslope(cslope_lo(1):cslope_hi(1), &
         cslope_lo(2):cslope_hi(2),size(fine,3))
    real(kind=dp_t) ::     lc_xslope(cslope_lo(1):cslope_hi(1), &
         cslope_lo(2):cslope_hi(2),size(fine,3))
    real(kind=dp_t) :: xslope_factor(cslope_lo(1):cslope_hi(1), &
         cslope_lo(2):cslope_hi(2))
    real(kind=dp_t) ::     uc_yslope(cslope_lo(1):cslope_hi(1), &
         cslope_lo(2):cslope_hi(2),size(fine,3))
    real(kind=dp_t) ::     lc_yslope(cslope_lo(1):cslope_hi(1), &
         cslope_lo(2):cslope_hi(2),size(fine,3))
    real(kind=dp_t) :: yslope_factor(cslope_lo(1):cslope_hi(1), &
         cslope_lo(2):cslope_hi(2))
    real(kind=dp_t) ::         alpha(cslope_lo(1):cslope_hi(1), &
         cslope_lo(2):cslope_hi(2),size(fine,3))
    real(kind=dp_t) ::          cmax(cslope_lo(1):cslope_hi(1), &
         cslope_lo(2):cslope_hi(2),size(fine,3))
    real(kind=dp_t) ::          cmin(cslope_lo(1):cslope_hi(1), &
         cslope_lo(2):cslope_hi(2),size(fine,3))
    real(kind=dp_t) ::         voffx(fvcx_lo:fvcx_lo+size(fine,1)-1)
    real(kind=dp_t) ::         voffy(fvcy_lo:fvcy_lo+size(fine,2)-1)

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

    do j = fvcy_lo, fvcy_lo+size(fine,2)-1
       jc = IX_PROJ(j,lratio(2))
       fycen = HALF*(fvcy(j)+fvcy(j+1))
       cycen = HALF*(cvcy(jc)+cvcy(jc+1))
       voffy(j) = (fycen-cycen)/(cvcy(jc+1)-cvcy(jc))
    end do
    do i = fvcx_lo, fvcx_lo+size(fine,1)-1
       ic = IX_PROJ(i,lratio(1))
       fxcen = HALF*(fvcx(i)+fvcx(i+1))
       cxcen = HALF*(cvcx(ic)+cvcx(ic+1))
       voffx(i) = (fxcen-cxcen)/(cvcx(ic+1)-cvcx(ic))
    end do

    ! Prevent underflow for small crse values.
    where ( abs(crse) <= 1.0e-20_dp_t ) crse = ZERO
    alpha = ONE
    cmax = crse(cslope_lo(1):cslope_hi(1),cslope_lo(2):cslope_hi(2),:)
    cmin = crse(cslope_lo(1):cslope_hi(1),cslope_lo(2):cslope_hi(2),:)

    do n = 1, size(fine,3)
       !
       ! Initialize alpha = 1 and define cmax and cmin as neighborhood max/mins.
       !
       do j = cslope_lo(2), cslope_hi(2)
          do i = cslope_lo(1), cslope_hi(1)
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
       do j = cslope_lo(2), cslope_hi(2)
          do i = cslope_lo(1), cslope_hi(1)
             uc_xslope(i,j,n) = HALF*(crse(i+1,j,n)-crse(i-1,j,n))
             lc_xslope(i,j,n) = uclc_slope(uc_xslope(i,j,n),crse,crse_lo,i,j,n,dim=1)
          end do
       end do

       if ( bc(1,1,n)  ==  EXT_DIR .or. bc(1,1,n) == HOEXTRAP ) then
          i = cslope_lo(1)
          if ( xok(1) ) then
             do j = cslope_lo(2), cslope_hi(2)
                uc_xslope(i,j,n) = -SIXTEEN/FIFTEEN*crse(i-1,j,n)+ HALF*crse(i,j,n) &
                     + TWO3RD*crse(i+1,j,n) - TENTH*crse(i+2,j,n)
             end do
          else
             do j = cslope_lo(2), cslope_hi(2)
                uc_xslope(i,j,n) = &
                     FOURTH * (crse(i+1,j,n) + FIVE*crse(i,j,n) - SIX*crse(i-1,j,n) )
             end do
          end if
          do j = cslope_lo(2), cslope_hi(2)
             lc_xslope(i,j,n) = uclc_slope(uc_xslope(i,j,n),crse,crse_lo,i,j,n,dim=1)
          end do
       end if

       if ( bc(1,2,n)  ==  EXT_DIR .or. bc(1,2,n) == HOEXTRAP ) then
          i = cslope_hi(1)
          if ( xok(1) ) then
             do j = cslope_lo(2), cslope_hi(2)
                uc_xslope(i,j,n) = SIXTEEN/FIFTEEN*crse(i+1,j,n)- HALF*crse(i,j,n) &
                     - TWO3RD*crse(i-1,j,n) + TENTH*crse(i-2,j,n)
             end do
          else
             do j = cslope_lo(2), cslope_hi(2)
                uc_xslope(i,j,n) = &
                     -FOURTH * (crse(i-1,j,n) + FIVE*crse(i,j,n) - SIX*crse(i+1,j,n) )
             end do
          end if
          do j = cslope_lo(2), cslope_hi(2)
             lc_xslope(i,j,n) = uclc_slope(uc_xslope(i,j,n),crse,crse_lo,i,j,n, dim=1)
          end do
       end if

       do j = cslope_lo(2), cslope_hi(2)
          do i = cslope_lo(1), cslope_hi(1)
             uc_yslope(i,j,n) = HALF*(crse(i,j+1,n)-crse(i,j-1,n))
             lc_yslope(i,j,n) = uclc_slope(uc_yslope(i,j,n),crse,crse_lo,i,j,n,dim=2)
          end do
       end do

       if ( bc(2,1,n)  == EXT_DIR .or. bc(2,1,n) == HOEXTRAP ) then
          j = cslope_lo(2)
          if ( xok(2) ) then
             do i = cslope_lo(1), cslope_hi(1)
                uc_yslope(i,j,n) = -SIXTEEN/FIFTEEN*crse(i,j-1,n)+ HALF*crse(i,j,n) &
                     + TWO3RD*crse(i,j+1,n) - TENTH*crse(i,j+2,n)
             end do
          else
             do i = cslope_lo(1), cslope_hi(1)
                uc_yslope(i,j,n) = &
                     FOURTH * (crse(i,j+1,n) + FIVE*crse(i,j,n) - SIX*crse(i,j-1,n) )
             end do
          end if
          do i = cslope_lo(1), cslope_hi(1)
             lc_yslope(i,j,n) = uclc_slope(uc_yslope(i,j,n),crse,crse_lo,i,j,n,dim=2)
          end do
       end if

       if ( bc(2,2,n)  ==  EXT_DIR .or. bc(2,2,n) == HOEXTRAP ) then
          j = cslope_hi(2)
          if ( xok(2) ) then
             do i = cslope_lo(1), cslope_hi(1)
                uc_yslope(i,j,n) = SIXTEEN/FIFTEEN*crse(i,j+1,n)- HALF*crse(i,j,n) &
                     - TWO3RD*crse(i,j-1,n) + TENTH*crse(i,j-2,n)
             end do
          else
             do i = cslope_lo(1), cslope_hi(1)
                uc_yslope(i,j,n) = &
                     -FOURTH * (crse(i,j-1,n) + FIVE*crse(i,j,n) - SIX*crse(i,j+1,n) )
             end do
          end if
          do i = cslope_lo(1), cslope_hi(1)
             lc_yslope(i,j,n) = uclc_slope(uc_yslope(i,j,n),crse,crse_lo,i,j,n,dim=2)
          end do
       end if

    end do

    if ( lim_slope ) then

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
             do j = fine_lo(2), fine_lo(2)+size(fine,2)-1

                jc = IX_PROJ(j,lratio(2))

                do i = fine_lo(1), fine_lo(1)+size(fine,1)-1
                   ic = IX_PROJ(i,lratio(1))

                   orig_corr_fact = voffx(i)*lc_xslope(ic,jc,n)+ voffy(j)*lc_yslope(ic,jc,n) 
                   fine(i,j,n) = crse(ic,jc,n) + orig_corr_fact
                   if ( fine(i,j,n)  >  cmax(ic,jc,n) ) then
                      if ( abs(orig_corr_fact)  >  1.e-10_dp_t*abs(crse(ic,jc,n)) ) then
                         corr_fact = (cmax(ic,jc,n) - crse(ic,jc,n)) / orig_corr_fact
                         alpha(ic,jc,n) = min(alpha(ic,jc,n),corr_fact)
                      end if
                   end if
                   if ( fine(i,j,n)  <  cmin(ic,jc,n) ) then
                      if ( abs(orig_corr_fact)  >  1.e-10_dp_t*abs(crse(ic,jc,n)) ) then
                         corr_fact = (cmin(ic,jc,n) - crse(ic,jc,n)) / orig_corr_fact
                         alpha(ic,jc,n) = min(alpha(ic,jc,n),corr_fact)
                      end if
                   end if

                end do
             end do
          end do

       end if

       ! Do the interpolation with limited slopes.

       do n = 1, size(fine,3)
          do j = fine_lo(2), fine_lo(2)+size(fine,2)-1
             jc = IX_PROJ(j,lratio(2))

             do i = fine_lo(1), fine_lo(1)+size(fine,1)-1
                ic = IX_PROJ(i,lratio(1))

                fine(i,j,n) = crse(ic,jc,n) + alpha(ic,jc,n)*( &
                     + voffx(i)*lc_xslope(ic,jc,n) &
                     + voffy(j)*lc_yslope(ic,jc,n) &
                     )
             end do
          end do
       end do

    else

       ! Do the interpolation using unlimited slopes.

       do n = 1, size(fine,3)
          do j = fine_lo(2), fine_lo(2)+size(fine,2)-1
             jc = IX_PROJ(j,lratio(2))
             do i = fine_lo(1), fine_lo(1)+size(fine,1)-1
                ic = IX_PROJ(i,lratio(1))
                fine(i,j,n) = &
                     crse(ic,jc,n) + voffx(i)*uc_xslope(ic,jc,n)+ voffy(j)*uc_yslope(ic,jc,n)
             end do
          end do
       end do
    end if

  contains

    function uclc_slope(uslope, crse, crse_lo, i, j, n, dim) result(lc)

      integer, intent(in) :: crse_lo(:)
      real(kind = dp_t) :: lc
      real(kind=dp_t), intent(in) :: uslope
      real(kind=dp_t), intent(in)  :: crse(crse_lo(1):,crse_lo(2):,:)
      real(kind=dp_t) :: cen, forw, back, slp
      integer, intent(in) :: i, j, n, dim
      cen  = uslope
      forw = TWO
      back = TWO
      if ( dim == 1 ) then
         forw = forw*(crse(i+1,j,n)-crse(i,j,n))
         back = back*(crse(i,j,n)-crse(i-1,j,n))
      else if ( dim == 2 ) then
         forw = forw*(crse(i,j+1,n)-crse(i,j,n))
         back = back*(crse(i,j,n)-crse(i,j-1,n))
      end if
      slp  = min(abs(forw),abs(back))
      if ( forw*back  <  ZERO ) then
         slp = ZERO
      end if
      lc = sign(ONE,cen)*min(slp,abs(cen))
    end function uclc_slope

  end subroutine lin_cc_interp_2d

  !! linccinterp:   linear conservative interpolation from coarse grid to

  subroutine lin_cc_interp_3d (fine, fine_lo, crse, crse_lo, &
                               lratio, bc, &
                               fvcx, fvcx_lo, fvcy, fvcy_lo, fvcz, fvcz_lo, &
                               cvcx, cvcx_lo, cvcy, cvcy_lo, cvcz, cvcz_lo, &
                               cslope_lo, cslope_hi, lim_slope, lin_limit)

    integer, intent(in) :: fine_lo(:),crse_lo(:)
    integer, intent(in) :: cslope_lo(:),cslope_hi(:)
    integer, intent(in) :: lratio(:)
    integer, intent(in) :: bc(:,:,:)
    integer, intent(in) :: fvcx_lo,fvcy_lo,fvcz_lo,cvcx_lo,cvcy_lo,cvcz_lo
    REAL(kind=dp_t), intent(out) :: fine(fine_lo(1):,fine_lo(2):,fine_lo(3):,:)
    REAL(kind=dp_t), intent(inout) :: crse(crse_lo(1):,crse_lo(2):,crse_lo(3):,:)
    REAL(kind=dp_t), intent(in) :: fvcx(fvcx_lo:)
    REAL(kind=dp_t), intent(in) :: fvcy(fvcy_lo:)
    REAL(kind=dp_t), intent(in) :: fvcz(fvcz_lo:)
    REAL(kind=dp_t), intent(in) :: cvcx(cvcx_lo:)
    REAL(kind=dp_t), intent(in) :: cvcy(cvcy_lo:)
    REAL(kind=dp_t), intent(in) :: cvcz(cvcz_lo:)
    logical        , intent(in) :: lim_slope, lin_limit

    REAL(kind=dp_t) :: voffx(fvcx_lo:fvcx_lo+size(fine,1)-1)
    REAL(kind=dp_t) :: voffy(fvcy_lo:fvcy_lo+size(fine,2)-1)
    REAL(kind=dp_t) :: voffz(fvcz_lo:fvcz_lo+size(fine,3)-1)

    integer :: n, i, ic, j, jc, k, kc, ioff, joff, koff
    REAL(kind=dp_t) :: fxcen, cxcen, fycen, cycen, fzcen, czcen
    REAL(kind=dp_t) :: corr_fact,orig_corr_fact
    logical :: xok(3)
    integer :: nxc(3)

    REAL(kind=dp_t), allocatable ::     uc_xslope(:,:,:,:)
    REAL(kind=dp_t), allocatable ::     lc_xslope(:,:,:,:)
    REAL(kind=dp_t), allocatable ::  slope_factor(:,:,:)
    REAL(kind=dp_t), allocatable ::     uc_yslope(:,:,:,:)
    REAL(kind=dp_t), allocatable ::     lc_yslope(:,:,:,:)
    REAL(kind=dp_t), allocatable ::     uc_zslope(:,:,:,:)
    REAL(kind=dp_t), allocatable ::     lc_zslope(:,:,:,:)
    REAL(kind=dp_t), allocatable ::         alpha(:,:,:,:)
    REAL(kind=dp_t), allocatable ::          cmax(:,:,:)
    REAL(kind=dp_t), allocatable ::          cmin(:,:,:)

    allocate(     uc_xslope(cslope_lo(1):cslope_hi(1),cslope_lo(2):cslope_hi(2), &
                            cslope_lo(3):cslope_hi(3),size(fine,4)))
    allocate(     lc_xslope(cslope_lo(1):cslope_hi(1),cslope_lo(2):cslope_hi(2), &
                            cslope_lo(3):cslope_hi(3),size(fine,4)))
    allocate(     uc_yslope(cslope_lo(1):cslope_hi(1),cslope_lo(2):cslope_hi(2), &
                            cslope_lo(3):cslope_hi(3),size(fine,4)))
    allocate(     lc_yslope(cslope_lo(1):cslope_hi(1),cslope_lo(2):cslope_hi(2), &
                            cslope_lo(3):cslope_hi(3),size(fine,4)))
    allocate(     uc_zslope(cslope_lo(1):cslope_hi(1),cslope_lo(2):cslope_hi(2), &
                            cslope_lo(3):cslope_hi(3),size(fine,4)))
    allocate(     lc_zslope(cslope_lo(1):cslope_hi(1),cslope_lo(2):cslope_hi(2), &
                            cslope_lo(3):cslope_hi(3),size(fine,4)))

    forall(i = 1:3) nxc(i) = size(crse,dim=i) - 2
    xok = (nxc >=  2)

    do k = fvcz_lo, fvcz_lo+size(fine,3)-1
       kc = IX_PROJ(k,lratio(3))
       fzcen = HALF*(fvcz(k)+fvcz(k+1))
       czcen = HALF*(cvcz(kc)+cvcz(kc+1))
       voffz(k) = (fzcen-czcen)/(cvcz(kc+1)-cvcz(kc))
    end do
    do j = fvcy_lo, fvcy_lo+size(fine,2)-1
       jc = IX_PROJ(j,lratio(2))
       fycen = HALF*(fvcy(j)+fvcy(j+1))
       cycen = HALF*(cvcy(jc)+cvcy(jc+1))
       voffy(j) = (fycen-cycen)/(cvcy(jc+1)-cvcy(jc))
    end do
    do i = fvcx_lo, fvcx_lo+size(fine,1)-1
       ic = IX_PROJ(i,lratio(1))
       fxcen = HALF*(fvcx(i)+fvcx(i+1))
       cxcen = HALF*(cvcx(ic)+cvcx(ic+1))
       voffx(i) = (fxcen-cxcen)/(cvcx(ic+1)-cvcx(ic))
    end do
    !
    ! Prevent underflow for small crse values.
    !
    where(abs(crse) < 1.e-20_dp_t) crse = ZERO
    !
    ! Computed unlimited and limited slopes
    !
    do n = 1, size(crse,4) 

       do k = cslope_lo(3),cslope_hi(3)
          do j = cslope_lo(2),cslope_hi(2)
             do i = cslope_lo(1),cslope_hi(1)
                uc_xslope(i,j,k,n) = HALF*(crse(i+1,j,k,n)-crse(i-1,j,k,n))
                lc_xslope(i,j,k,n) = &
                     uclc_slope(uc_xslope(i,j,k,n),crse,crse_lo,i,j,k,n,dim=1)
             end do
          end do
       end do

       if ( bc(1,1,n)  ==  EXT_DIR .or. bc(1,1,n) == HOEXTRAP ) then
          i = cslope_lo(1)
          if ( xok(1) ) then
             do k = cslope_lo(3),cslope_hi(3)
                do j = cslope_lo(2),cslope_hi(2)
                   uc_xslope(i,j,k,n) = -SIXTEEN/FIFTEEN*crse(i-1,j,k,n) &
                        + HALF*crse(i,j,k,n)+ TWO3RD*crse(i+1,j,k,n) - TENTH*crse(i+2,j,k,n)
                end do
             end do
          else
             do k = cslope_lo(3),cslope_hi(3)
                do j = cslope_lo(2),cslope_hi(2)
                   uc_xslope(i,j,k,n) = &
                        FOURTH * (crse(i+1,j,k,n) + FIVE*crse(i,j,k,n) - SIX*crse(i-1,j,k,n))
                end do
             end do
          end if
          do k = cslope_lo(3),cslope_hi(3)
             do j = cslope_lo(2),cslope_hi(2)
                lc_xslope(i,j,k,n) = &
                     uclc_slope(uc_xslope(i,j,k,n),crse,crse_lo,i,j,k,n,dim=1)
             end do
          end do
       end if

       if ( bc(1,2,n)  ==  EXT_DIR .or. bc(1,2,n) == HOEXTRAP ) then
          i = cslope_hi(1)
          if ( xok(1) ) then
             do k = cslope_lo(3),cslope_hi(3)
                do j = cslope_lo(2),cslope_hi(2)
                   uc_xslope(i,j,k,n) = SIXTEEN/FIFTEEN*crse(i+1,j,k,n) &
                        - HALF*crse(i,j,k,n)- TWO3RD*crse(i-1,j,k,n) + TENTH*crse(i-2,j,k,n)
                end do
             end do
          else 
             do k = cslope_lo(3),cslope_hi(3)
                do j = cslope_lo(2),cslope_hi(2)
                   uc_xslope(i,j,k,n) = -FOURTH * (crse(i-1,j,k,n) &
                        + FIVE*crse(i,j,k,n) - SIX*crse(i+1,j,k,n) )
                end do
             end do
          end if
          do k = cslope_lo(3),cslope_hi(3)
             do j = cslope_lo(2),cslope_hi(2)
                lc_xslope(i,j,k,n) = &
                     uclc_slope(uc_xslope(i,j,k,n),crse,crse_lo,i,j,k,n,dim=1)
             end do
          end do
       end if

       do k = cslope_lo(3),cslope_hi(3)
          do j = cslope_lo(2),cslope_hi(2)
             do i = cslope_lo(1),cslope_hi(1)
                uc_yslope(i,j,k,n) = HALF*(crse(i,j+1,k,n)-crse(i,j-1,k,n))
                lc_yslope(i,j,k,n) = &
                     uclc_slope(uc_yslope(i,j,k,n),crse,crse_lo,i,j,k,n,dim=2)
             end do
          end do
       end do

       if ( bc(2,1,n)  ==  EXT_DIR .or. bc(2,1,n) == HOEXTRAP ) then
          j = cslope_lo(2)
          if ( xok(2) ) then
             do k = cslope_lo(3),cslope_hi(3)
                do i = cslope_lo(1),cslope_hi(1)
                   uc_yslope(i,j,k,n) = -SIXTEEN/FIFTEEN*crse(i,j-1,k,n) &
                        + HALF*crse(i,j,k,n)+ TWO3RD*crse(i,j+1,k,n) - TENTH*crse(i,j+2,k,n)
                end do
             end do
          else
             do k = cslope_lo(3),cslope_hi(3)
                do i = cslope_lo(1),cslope_hi(1)
                   uc_yslope(i,j,k,n) = &
                        FOURTH * (crse(i,j+1,k,n) + FIVE*crse(i,j,k,n) - SIX*crse(i,j-1,k,n))
                end do
             end do
          end if
          do k = cslope_lo(3),cslope_hi(3)
             do i = cslope_lo(1),cslope_hi(1)
                lc_yslope(i,j,k,n) = &
                     uclc_slope(uc_yslope(i,j,k,n),crse,crse_lo,i,j,k,n,dim=2)
             end do
          end do
       end if

       if ( bc(2,2,n)  ==  EXT_DIR .or. bc(2,2,n) == HOEXTRAP ) then
          j = cslope_hi(2)
          if ( xok(2) ) then
             do k = cslope_lo(3),cslope_hi(3)
                do i = cslope_lo(1),cslope_hi(1)
                   uc_yslope(i,j,k,n) = SIXTEEN/FIFTEEN*crse(i,j+1,k,n) &
                        - HALF*crse(i,j,k,n)- TWO3RD*crse(i,j-1,k,n) + TENTH*crse(i,j-2,k,n)
                end do
             end do
          else
             do k = cslope_lo(3),cslope_hi(3)
                do i = cslope_lo(1),cslope_hi(1)
                   uc_yslope(i,j,k,n) = -FOURTH * (crse(i,j-1,k,n) &
                        + FIVE*crse(i,j,k,n) - SIX*crse(i,j+1,k,n) )
                end do
             end do
          end if
          do k = cslope_lo(3),cslope_hi(3)
             do i = cslope_lo(1),cslope_hi(1)
                lc_yslope(i,j,k,n) = &
                     uclc_slope(uc_yslope(i,j,k,n),crse,crse_lo,i,j,k,n,dim=2)
             end do
          end do
       end if

       do k = cslope_lo(3),cslope_hi(3)
          do j = cslope_lo(2),cslope_hi(2)
             do i = cslope_lo(1),cslope_hi(1)
                uc_zslope(i,j,k,n) = HALF*(crse(i,j,k+1,n)-crse(i,j,k-1,n))
                lc_zslope(i,j,k,n) = &
                     uclc_slope(uc_zslope(i,j,k,n),crse,crse_lo,i,j,k,n,dim=3)
             end do
          end do
       end do

       if ( bc(3,1,n)  == EXT_DIR .or. bc(3,1,n) == HOEXTRAP ) then
          k = cslope_lo(3)
          if ( xok(3) ) then
             do j = cslope_lo(2),cslope_hi(2)
                do i = cslope_lo(1),cslope_hi(1)
                   uc_zslope(i,j,k,n) = -SIXTEEN/FIFTEEN*crse(i,j,k-1,n) &
                        + HALF*crse(i,j,k,n)+ TWO3RD*crse(i,j,k+1,n) - TENTH*crse(i,j,k+2,n)
                end do
             end do
          else
             do j = cslope_lo(2),cslope_hi(2)
                do i = cslope_lo(1),cslope_hi(1)
                   uc_zslope(i,j,k,n) = &
                        FOURTH * (crse(i,j,k+1,n) + FIVE*crse(i,j,k,n) - SIX*crse(i,j,k-1,n))
                end do
             end do
          end if
          do j = cslope_lo(2),cslope_hi(2)
             do i = cslope_lo(1),cslope_hi(1)
                lc_zslope(i,j,k,n) = &
                     uclc_slope(uc_zslope(i,j,k,n),crse,crse_lo,i,j,k,n,dim=3)
             end do
          end do
       end if

       if ( bc(3,2,n)  ==  EXT_DIR .or. bc(3,2,n)  ==  HOEXTRAP ) then
          k = cslope_hi(3)
          if ( xok(3) ) then
             do j = cslope_lo(2),cslope_hi(2)
                do i = cslope_lo(1),cslope_hi(1)
                   uc_zslope(i,j,k,n) = SIXTEEN/FIFTEEN*crse(i,j,k+1,n) &
                        - HALF*crse(i,j,k,n)- TWO3RD*crse(i,j,k-1,n) + TENTH*crse(i,j,k-2,n)
                end do
             end do
          else
             do j = cslope_lo(2),cslope_hi(2)
                do i = cslope_lo(1),cslope_hi(1)
                   uc_zslope(i,j,k,n) = -FOURTH * (crse(i,j,k-1,n) &
                        + FIVE*crse(i,j,k,n) - SIX*crse(i,j,k+1,n) )
                end do
             end do
          end if
          do j = cslope_lo(2),cslope_hi(2)
             do i = cslope_lo(1),cslope_hi(1)
                lc_zslope(i,j,k,n) = &
                     uclc_slope(uc_zslope(i,j,k,n),crse,crse_lo,i,j,k,n,dim=3)
             end do
          end do
       end if

    end do

    if ( lim_slope ) then

       allocate(alpha(cslope_lo(1):cslope_hi(1),cslope_lo(2):cslope_hi(2), &
            cslope_lo(3):cslope_hi(3),size(fine,4)))

       alpha = ONE

       if ( lin_limit ) then

          ! compute linear limited slopes
          ! Note that the limited and the unlimited slopes
          ! have the same sign, and it is assumed that they do.
          !
          allocate( slope_factor(cslope_lo(1):cslope_hi(1),cslope_lo(2):cslope_hi(2), &
               cslope_lo(3):cslope_hi(3)) )
          !
          ! x slopes
          !
          slope_factor = ONE
          do n = 1, size(crse,4) 
             where( uc_xslope(:,:,:,n) /= 0 )
                slope_factor = min(slope_factor,lc_xslope(:,:,:,n)/uc_xslope(:,:,:,n),ONE)
             end where
          end do
          do n = 1, size(crse,4) 
             lc_xslope(:,:,:,n) = slope_factor*uc_xslope(:,:,:,n)
          end do
          !
          ! y slopes
          !
          slope_factor = ONE
          do n = 1, size(crse,4) 
             where( uc_yslope(:,:,:,n) /= 0 )
                slope_factor = min(slope_factor,lc_yslope(:,:,:,n)/uc_yslope(:,:,:,n),ONE)
             end where
          end do
          do n = 1, size(crse,4) 
             lc_yslope(:,:,:,n) = slope_factor*uc_yslope(:,:,:,n)
          end do
          !
          ! z slopes
          !
          slope_factor = ONE
          do n = 1, size(crse,4) 
             where( uc_zslope(:,:,:,n) /= 0 )
                slope_factor = min(slope_factor,lc_zslope(:,:,:,n)/uc_zslope(:,:,:,n),ONE)
             end where
          end do
          do n = 1, size(crse,4) 
             lc_zslope(:,:,:,n) = slope_factor*uc_zslope(:,:,:,n)
          end do

          deallocate(slope_factor)

       else

          allocate( cmax(cslope_lo(1):cslope_hi(1),cslope_lo(2):cslope_hi(2), &
               cslope_lo(3):cslope_hi(3)) )
          allocate( cmin(cslope_lo(1):cslope_hi(1),cslope_lo(2):cslope_hi(2), &
               cslope_lo(3):cslope_hi(3)) )

          do n = 1, size(crse,4) 
             !
             ! Define cmax and cmin as neighborhood max/mins.
             !
             do k = cslope_lo(3),cslope_hi(3)
                do j = cslope_lo(2),cslope_hi(2)
                   do i = cslope_lo(1),cslope_hi(1)
                      cmax(i,j,k) = crse(i,j,k,n)
                      cmin(i,j,k) = crse(i,j,k,n)
                      do koff = -1,1
                         do joff = -1,1
                            do ioff = -1,1
                               cmax(i,j,k) = max(cmax(i,j,k),crse(i+ioff,j+joff,k+koff,n))
                               cmin(i,j,k) = min(cmin(i,j,k),crse(i+ioff,j+joff,k+koff,n))
                            end do
                         end do
                      end do
                   end do
                end do
             end do
             !
             ! Limit slopes so as to not introduce new maxs or mins.
             !
             do k = fine_lo(3), fine_lo(3)+size(fine, 3) - 1
                kc = IX_PROJ(k,lratio(3))
                do j = fine_lo(2), fine_lo(2)+size(fine, 2) - 1
                   jc = IX_PROJ(j,lratio(2))
                   do i = fine_lo(1), fine_lo(1)+size(fine, 1) - 1
                      ic = IX_PROJ(i,lratio(1))
                      orig_corr_fact = &
                           + voffx(i)*lc_xslope(ic,jc,kc,n) &
                           + voffy(j)*lc_yslope(ic,jc,kc,n) &
                           + voffz(k)*lc_zslope(ic,jc,kc,n)
                      fine(i,j,k,n) = crse(ic,jc,kc,n) + orig_corr_fact
                      if ( fine(i,j,k,n) > cmax(ic,jc,kc) ) then
                         if ( abs(orig_corr_fact) > 1.e-10_dp_t*abs(crse(ic,jc,kc,n)) ) then
                            corr_fact = (cmax(ic,jc,kc) - crse(ic,jc,kc,n)) / orig_corr_fact
                            alpha(ic,jc,kc,n) = min(alpha(ic,jc,kc,n),corr_fact)
                         end if
                      end if
                      if ( fine(i,j,k,n)  <  cmin(ic,jc,kc) ) then
                         if ( abs(orig_corr_fact)  >  1.e-10_dp_t*abs(crse(ic,jc,kc,n)) ) then
                            corr_fact = (cmin(ic,jc,kc) - crse(ic,jc,kc,n)) / orig_corr_fact
                            alpha(ic,jc,kc,n) = min(alpha(ic,jc,kc,n),corr_fact)
                         end if
                      end if
                   end do
                end do
             end do
          end do

          deallocate(cmax,cmin)
       end if
       !
       ! Do the interpolation with limited slopes.
       !
       do n = 1, size(crse,4)
          do k = fine_lo(3), fine_lo(3)+size(fine, 3) - 1
             kc = IX_PROJ(k,lratio(3))
             do j = fine_lo(2), fine_lo(2)+size(fine, 2) - 1
                jc = IX_PROJ(j,lratio(2))
                do i = fine_lo(1), fine_lo(1)+size(fine, 1) - 1
                   ic = IX_PROJ(i,lratio(1))
                   fine(i,j,k,n) = crse(ic,jc,kc,n) &
                        + alpha(ic,jc,kc,n) * &
                        ( &
                        + voffx(i)*lc_xslope(ic,jc,kc,n)&
                        + voffy(j)*lc_yslope(ic,jc,kc,n)&
                        + voffz(k)*lc_zslope(ic,jc,kc,n)&
                        )
                end do
             end do
          end do
       end do

       deallocate(alpha)

    else
       !
       ! Do the interpolation using unlimited slopes.
       !
       do n = 1, size(crse,4)
          do k = fine_lo(3), fine_lo(3)+size(fine,3) - 1
             kc = IX_PROJ(k,lratio(3))
             do j = fine_lo(2), fine_lo(2)+size(fine, 2) - 1
                jc = IX_PROJ(j,lratio(2))
                do i = fine_lo(1), fine_lo(1)+size(fine,1) - 1
                   ic = IX_PROJ(i,lratio(1))
                   fine(i,j,k,n) = crse(ic,jc,kc,n) &
                        + voffx(i)*uc_xslope(ic,jc,kc,n) &
                        + voffy(j)*uc_yslope(ic,jc,kc,n) &
                        + voffz(k)*uc_zslope(ic,jc,kc,n)
                end do
             end do
          end do
       end do

    end if

    deallocate(uc_xslope,lc_xslope)
    deallocate(uc_yslope,lc_yslope)
    deallocate(uc_zslope,lc_zslope)


  contains

    function uclc_slope(uslope, crse, crse_lo, i, j, k, n, dim) result(lc)
      integer, intent(in) :: crse_lo(:)
      real(kind = dp_t) :: lc
      real(kind=dp_t), intent(in) :: uslope
      real(kind=dp_t), intent(in)  :: crse(crse_lo(1):,crse_lo(2):,crse_lo(3):,:)
      real(kind=dp_t) :: cen, forw, back, slp
      integer, intent(in) :: i, j, k, n, dim
      cen  = uslope
      forw = TWO
      back = TWO
      if ( dim == 1 ) then
         forw = forw*(crse(i+1,j,k,n)-crse(i,j,k,n))
         back = back*(crse(i,j,k,n)-crse(i-1,j,k,n))
      else if ( dim == 2 ) then
         forw = forw*(crse(i,j+1,k,n)-crse(i,j,k,n))
         back = back*(crse(i,j,k,n)-crse(i,j-1,k,n))
      else if ( dim == 3 ) then
         forw = forw*(crse(i,j,k+1,n)-crse(i,j,k,n))
         back = back*(crse(i,j,k,n)-crse(i,j,k-1,n))
      end if
      slp  = min(abs(forw),abs(back))
      if ( forw*back  <  ZERO ) then
         slp = ZERO
      end if
      lc = sign(ONE,cen)*min(slp,abs(cen))
    end function uclc_slope

  end subroutine lin_cc_interp_3d

  !! fourth_order_interp:   fourth order unlimited interpolation from coarse grid to fine
  !!                    :   NOTE: this is not necessarily conservative!

  subroutine fourth_order_interp_2d (fine, fine_lo, crse, crse_lo, lratio, &
       cg_lo, cg_hi)

    use fourth_order_interp_coef_module

    implicit none

    integer, intent(in) :: fine_lo(:), crse_lo(:), cg_lo(:), cg_hi(:)
    integer, intent(in) :: lratio(:)
    real(kind=dp_t), intent(  out) :: fine(fine_lo(1):,fine_lo(2):,:)
    real(kind=dp_t), intent(inout) :: crse(crse_lo(1):,crse_lo(2):,:)

    real(kind=dp_t) :: b(21), A2T(size(A2,2),0:size(A2,1)-1), c(15)

    integer :: n, i, j, ic, jc, k
    !
    ! Prevent underflow for small crse values.
    !
    where ( abs(crse) <= 1.0e-20_dp_t ) crse = ZERO
    !
    ! Use A2T instead of A2 for more efficient memory access.
    !
    A2T = Transpose(A2)
    !
    ! Do interpolation.  For this method of interpolation, it is more
    ! efficient to loop over the coarse cells than the fine cells
    ! since the b and c vectors depend on the coarse cell, not on the
    ! fine cell.
    !
    do n = 1, size(fine,3)
       do jc = cg_lo(2), cg_hi(2)
          do ic = cg_lo(1), cg_hi(1)
             b(  1) = crse(ic-2, jc-1, n)
             b(  2) = crse(ic-2, jc+0, n)
             b(  3) = crse(ic-2, jc+1, n)
             b(  4) = crse(ic-1, jc-2, n)
             b(  5) = crse(ic-1, jc-1, n)
             b(  6) = crse(ic-1, jc+0, n)
             b(  7) = crse(ic-1, jc+1, n)
             b(  8) = crse(ic-1, jc+2, n)
             b(  9) = crse(ic+0, jc-2, n)
             b( 10) = crse(ic+0, jc-1, n)
             b( 11) = crse(ic+0, jc+1, n)
             b( 12) = crse(ic+0, jc+2, n)
             b( 13) = crse(ic+1, jc-2, n)
             b( 14) = crse(ic+1, jc-1, n)
             b( 15) = crse(ic+1, jc+0, n)
             b( 16) = crse(ic+1, jc+1, n)
             b( 17) = crse(ic+1, jc+2, n)
             b( 18) = crse(ic+2, jc-1, n)
             b( 19) = crse(ic+2, jc+0, n)
             b( 20) = crse(ic+2, jc+1, n)
             b( 21) = 1000.0_dp_t*crse(ic+0, jc+0, n)

             c = matmul(P2, b)

             do j = jc*lratio(2),(jc+1)*lratio(2)-1
                do i = ic*lratio(1),(ic+1)*lratio(1)-1

                   k = 2*(i-ic*lratio(1)) + (j-jc*lratio(2))

                   fine(i,j,n) = dot_product(c, A2T(:,k))*4
                   
                end do
             end do

          end do
       end do
    end do

    ! print *, 'COARSE'
    ! print *, crse

    ! print *, 'FINE'
    ! print *, fine

  end subroutine fourth_order_interp_2d

end module interp_module
