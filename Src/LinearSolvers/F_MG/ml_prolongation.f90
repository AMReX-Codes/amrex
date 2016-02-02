module ml_prolongation_module

  use bl_constants_module
  use bl_types
  use multifab_module

  implicit none

  real(dp_t), private, parameter :: THREE_EIGHTHS= 0.375_dp_t

  private

  public :: ml_cc_prolongation, ml_nodal_prolongation, ml_interp_bcs

contains

  subroutine ml_cc_prolongation(fine, crse, ir, lininterp, ptype)
    use bl_prof_module
    use mg_prolongation_module
    use multifab_module, only: get_dim

    type(multifab), intent(inout) :: fine
    type(multifab), intent(in   ) :: crse
    integer,        intent(in   ) :: ir(:), ptype
    logical,        intent(in   ) :: lininterp

    integer             :: lo (get_dim(fine)), hi (get_dim(fine))
    integer             :: loc(get_dim(fine)), lof(get_dim(fine))
    integer             :: i, n, dm
    real(dp_t), pointer :: fp(:,:,:,:), cp(:,:,:,:)
    type(boxarray)      :: ba
    type(layout)        :: lacfine, laf, latmp
    type(multifab)      :: cfine, cfgrow
    type(bl_prof_timer), save :: bpt

    if ( ncomp(crse) .ne. ncomp(fine) ) then
       call bl_error('ml_cc_prolongation: crse & fine must have same # of components')
    end if

    call build(bpt, "ml_cc_prolongation")

    laf = get_layout(fine)

    call layout_build_coarse(lacfine, laf, ir)

    if ( .false. ) then
       !
       ! I'm checking this in disabled.  I still want to do more testing.
       ! The point here is to get better meshing together to grids when
       ! we prolongate by using one grow cell, which "should" always be
       ! available due to proper nesting of grids.
       !
       call boxarray_build_copy(ba, get_boxarray(lacfine))

       call boxarray_grow(ba,1)

       call layout_build_ba(latmp, ba, get_pd(lacfine), explicit_mapping = lacfine%lap%prc)

       call boxarray_destroy(ba)
       !
       ! "crse" should always cover "cfine" (including ghost cells) on
       ! the valid region due to proper nesting.
       !
       call multifab_build(cfine, lacfine, nc = ncomp(crse), ng = 1)
       !
       ! Use this to fill cfine including ghost cells.
       !
       call multifab_build(cfgrow, latmp, nc = ncomp(crse), ng = 0)

       call copy(cfgrow, 1, crse, 1, ncomp(crse))
       !
       ! Finally copy FAB-wise from cfgrow -> cfine to get good ghost cells.
       !
       do i = 1, nfabs(cfgrow)
          fp => dataptr(cfgrow, i)
          cp => dataptr(cfine,  i)
          call cpy_d(cp,fp)
       end do

       call multifab_destroy(cfgrow)
       call layout_destroy(latmp)
    else
       call multifab_build(cfine, lacfine, nc = ncomp(crse), ng = 0)
       call copy(cfine, 1, crse, 1, ncomp(crse))
    end if

    dm = get_dim(crse)

    !$OMP PARALLEL DO PRIVATE(i,loc,lof,lo,hi,n,fp,cp)
    do i = 1, nfabs(fine)
       loc =  lwb(get_pbox(cfine,i))
       lof =  lwb(get_pbox(fine, i))
       lo  =  lwb(get_ibox(fine, i))
       hi  =  upb(get_ibox(fine, i))
       fp  => dataptr(fine, i)
       cp  => dataptr(cfine,i)
       do n = 1, ncomp(crse)
          select case ( dm )
          case (1)
             call cc_prolongation(fp(:,1,1,n), lof, cp(:,1,1,n), loc, lo, hi, ir, lininterp)
          case (2)
             call cc_prolongation(fp(:,:,1,n), lof, cp(:,:,1,n), loc, lo, hi, ir, lininterp, ptype)
          case (3)
             call cc_prolongation(fp(:,:,:,n), lof, cp(:,:,:,n), loc, lo, hi, ir, lininterp, ptype)
          end select
       end do
    end do
    !$OMP END PARALLEL DO

    call multifab_destroy(cfine)
    call destroy(bpt)

  end subroutine ml_cc_prolongation

  subroutine ml_nodal_prolongation(fine, crse, ir)

    use fabio_module
    use bl_prof_module
    use mg_prolongation_module

    type(multifab), intent(inout) :: fine
    type(multifab), intent(in   ) :: crse
    integer,        intent(in   ) :: ir(:)

    integer             :: lo (get_dim(fine)), hi (get_dim(fine))
    integer             :: loc(get_dim(fine)), lof(get_dim(fine))
    integer             :: i, n, dm
    real(dp_t), pointer :: fp(:,:,:,:), cp(:,:,:,:)
    type(boxarray)      :: ba
    type(layout)        :: lacfine, laf, latmp
    type(multifab)      :: cfine, cfgrow
    character(len=3)    :: number
    character(len=20)   :: filename

!    logical oldval

    integer, save :: cnt = 0

    type(bl_prof_timer), save :: bpt

    if ( ncomp(crse) .ne. ncomp(fine) ) then
       call bl_error('ml_nodal_prolongation: crse & fine must have same # of components')
    end if

    call build(bpt, "ml_nodal_prolongation")

    laf = get_layout(fine)

    call layout_build_coarse(lacfine, laf, ir)

!    oldval = set_using_nodal_cubic(.true.)

    if ( using_nodal_cubic() .and. get_dim(crse) > 1 ) then
       !
       ! I'm checking this in disabled.  I still want to do more testing.
       ! The point here is to get better meshing together to grids when
       ! we prolongate by using one grow cell, which "should" always be
       ! available due to proper nesting of grids.
       !
       call boxarray_build_copy(ba, get_boxarray(lacfine))

       call boxarray_grow(ba,1)

       call layout_build_ba(latmp, ba, get_pd(lacfine), explicit_mapping = lacfine%lap%prc)

       call boxarray_destroy(ba)
       !
       ! "crse" should always cover "cfine" (including ghost cells) on
       ! the valid region due to proper nesting.
       !
       call multifab_build(cfine, lacfine, nc = ncomp(crse), ng = 1, nodal = nodal_flags(fine))
       !
       ! Use this to fill cfine including ghost cells.
       !
       call multifab_build(cfgrow, latmp, nc = ncomp(crse), ng = 0, nodal = nodal_flags(fine))

       call copy(cfgrow, 1, crse, 1, ncomp(crse))
       !
       ! Finally copy FAB-wise from cfgrow -> cfine to get good ghost cells.
       !
       do i = 1, nfabs(cfgrow)
          fp => dataptr(cfgrow, i)
          cp => dataptr(cfine,  i)
          call cpy_d(cp,fp)
       end do

       call multifab_destroy(cfgrow)
       call layout_destroy(latmp)
    else
       call multifab_build(cfine, lacfine, nc = ncomp(crse), ng = 0, nodal = nodal_flags(fine))
       call copy(cfine, 1, crse, 1, ncomp(crse))
    end if

    dm = get_dim(crse)

    !$OMP PARALLEL DO PRIVATE(i,loc,lof,lo,hi,n,fp,cp)
    do i = 1, nfabs(fine)
       loc =  lwb(get_pbox(cfine,i))
       lof =  lwb(get_pbox(fine, i))
       lo  =  lwb(get_ibox(fine, i))
       hi  =  upb(get_ibox(fine, i))
       fp  => dataptr(fine,  i)
       cp  => dataptr(cfine, i)
       do n = 1, ncomp(crse)
          select case (dm)
          case (1)
             call nodal_prolongation_1d(fp(:,1,1,n), lof, cp(:,1,1,n), loc, lo, hi, ir)
          case (2)
             call nodal_prolongation_2d(fp(:,:,1,n), lof, cp(:,:,1,n), loc, lo, hi, ir)
          case (3)
             call nodal_prolongation_3d(fp(:,:,:,n), lof, cp(:,:,:,n), loc, lo, hi, ir)
          end select
       end do
    end do
    !$OMP END PARALLEL DO

!    oldval = set_using_nodal_cubic(oldval)

    if ( .false. ) then
       !
       ! Some debugging code I want to keep around for a while.
       ! I don't want to have to recreate this all the time :-)
       !
       cnt = cnt + 1
       write(number,fmt='(i3.3)') cnt
       filename = 'prolong_' // number
       call fabio_write(fine, 'debug', trim(filename))
    end if

    call multifab_destroy(cfine)
    call destroy(bpt)

  end subroutine ml_nodal_prolongation

  subroutine ml_interp_bcs_c(fine, cf, crse, cc, fine_domain, ir, ng_fill, nc, facemap, indxmap, uncovered)
    use bl_prof_module
    type(multifab), intent(inout) :: fine
    type(multifab), intent(in   ) :: crse
    type(box),      intent(in)    :: fine_domain
    integer,        intent(in)    :: ir(:)
    integer,        intent(in)    :: cf, cc, ng_fill, nc, facemap(:), indxmap(:)
    type(lmultifab),intent(in)    :: uncovered

    type(box) :: fbox, fbox_grown, cbox_refined, isect
    integer   :: lo (get_dim(fine)), hi (get_dim(fine))
    integer   :: loc(get_dim(fine)), hic(get_dim(fine))
    integer   :: lof(get_dim(fine))
    integer   :: dm, i, j, n, dir, face, side
    logical   :: pmask(get_dim(fine))

    real(dp_t), pointer :: fp(:,:,:,:), cp(:,:,:,:)
    logical, pointer :: mp(:,:,:,:)
    
    type(bl_prof_timer), save :: bpt

    call build(bpt, "ml_interp_bcs_c")

    dm   = get_dim(fine)
    pmask = get_pmask(get_layout(fine))

    ! this loop cannot be OMP'd because there are overlaps in fbox
    ! If we really need to OMP, we could pass in bndry_reg instead of crse mf
    do i = 1, nfabs(crse)
       side = facemap(i)
       dir  = abs(side)
       face = sign(1, side)

       j = indxmap(i)  ! index for fine

       loc          = lwb(get_pbox(crse, i))
       hic          = upb(get_pbox(crse, i))
       cbox_refined = refine(get_ibox(crse,i), ir)

       fbox         = get_ibox(fine, j)
       lof          = lwb(get_pbox(fine, j))

       if ( .not. nodal_q(crse) ) then
          fbox_grown = grow(fbox, ng_fill, dir, face)
          if (.not. pmask(dir)) fbox_grown = intersection(fbox_grown, fine_domain)
       else
          fbox_grown = fbox
       end if

       isect = intersection(fbox_grown, cbox_refined)

       if ( empty(isect) ) cycle

       lo =  lwb(isect)
       hi =  upb(isect)
       fp => dataptr(fine, j)
       cp => dataptr(crse, i)

       mp => dataptr(uncovered, i)

       do n = 0, nc - 1
          select case (dm)
          case (1)
             if ( cell_centered_q(crse) ) then
                call ml_interp_bcs_1d(fp(:,1,1,cf+n), lof, cp(:,1,1,cc+n), loc, lo, ir)
             else if ( nodal_q(crse) ) then
                call bl_error("ML_INTERP_BCS: nodal 1d not provided")
             else
                call bl_error("ML_INTERP_BCS: nodal or cell centered only")
             end if
          case (2)
             if ( cell_centered_q(crse) ) then
                call ml_interp_bcs_2d(fp(:,:,1,cf+n), lof, cp(:,:,1,cc+n), mp(:,:,1,1), loc, hic, lo, hi, ir, side)
             else if ( nodal_q(crse) ) then
!                call ml_interp_bcs_2d_nodal(fp(:,:,1,cf+n), lof, cp(:,:,1,cc+n), loc, lo, hi, ir, side)
                call bl_error("ML_INTERP_BCS: nodal 2d not provided")
             else
                call bl_error("ML_INTERP_BCS: nodal or cell centered only")
             end if
          case (3)
             if ( cell_centered_q(crse) ) then
                call ml_interp_bcs_3d(fp(:,:,:,cf+n), lof, cp(:,:,:,cc+n), mp(:,:,:,1), loc, hic, lo, hi, ir, side)
             else if ( nodal_q(crse) )  then
!                call ml_interp_bcs_3d_nodal(fp(:,:,:,cf+n), lof, cp(:,:,:,cc+n), loc, lo, hi, ir, side)
                call bl_error("ML_INTERP_BCS: nodal 3d not provided")
             else
                call bl_error("ML_INTERP_BCS: nodal or cell centered only")
             end if
          end select
       end do
    end do

!   We have moved this call to the routine which calls ml_interp_bcs.
!   call multifab_fill_boundary(fine)

    call destroy(bpt)

  end subroutine ml_interp_bcs_c

  subroutine ml_interp_bcs(fine, crse, fine_domain, ir, ng_fill, facemap, indxmap, uncovered)
    type(multifab), intent(inout) :: fine
    type(multifab), intent(in   ) :: crse
    type(box),      intent(in)    :: fine_domain
    integer,        intent(in)    :: ir(:), ng_fill, facemap(:), indxmap(:)
    type(lmultifab),intent(in)   :: uncovered  ! could make this optional if needed
    if ( ncomp(crse) .ne. ncomp(fine) ) then
       call bl_error('ml_interp_bcs: crse & fine must have same # of components')
    end if
!    if (.not. present(uncovered) .and. cell_centered_q(crse) ) then
!       call bl_error("ml_interp_bcs_c: must have uncovered in cell-centered case")
!    end if
    call ml_interp_bcs_c(fine, 1, crse, 1, fine_domain, ir, ng_fill, ncomp(fine), facemap, indxmap, uncovered)
  end subroutine ml_interp_bcs

  subroutine ml_interp_bcs_1d(ff, lof, cc, loc, lo, ir)
    integer, intent(in) :: loc(:)
    integer, intent(in) :: lof(:)
    integer, intent(in) :: lo(:)
    real (dp_t), intent(inout) :: ff(lof(1):)
    real (dp_t), intent(in) :: cc(loc(1):)
    integer, intent(in) :: ir(:)
    integer :: i, ic
    i     = lo(1)
    ic    = i / ir(1)
    ff(i) = cc(ic)
  end subroutine ml_interp_bcs_1d

  subroutine ml_interp_bcs_2d(ff, lof, cc, uc, loc, hic, lo, hi, ir, side)
    integer, intent(in) :: loc(:), hic(:)
    integer, intent(in) :: lof(:)
    integer, intent(in) :: lo(:), hi(:)
    real (dp_t), intent(inout) :: ff(lof(1):,lof(2):)
    real (dp_t), intent(inout) :: cc(loc(1):,loc(2):)
    logical    , intent(inout) :: uc(loc(1):,loc(2):)
    integer, intent(in) :: side
    integer, intent(in) :: ir(:)
    integer :: i, j, ic, jc, m
    real (dp_t) :: xloc(0:3)
    real (dp_t) :: first_der, second_der

    if (ir(1) == 2) then
       xloc(0) = -FOURTH
       xloc(1) =  FOURTH
    else
       xloc(0) = -THREE_EIGHTHS
       xloc(1) = -EIGHTH
       xloc(2) =  EIGHTH
       xloc(3) =  THREE_EIGHTHS
    end if

!   This is the same interpolation as done in IAMR.
    if ( abs(side) == 1 ) then
      call bl_assert(size(cc,dim=1) .eq. 1, "cc must be one cell wide in ml_interp_bcs_2d")
      ic = loc(1)
      do jc = lo(2)/ir(2), hi(2)/ir(2)
        first_der  = ZERO
        second_der = ZERO
        if ( uc(ic,jc-1) .and. uc(ic,jc+1) ) then
          first_der  = HALF * (cc(ic,jc+1)-cc(ic,jc-1))
          second_der = (cc(ic,jc+1)+cc(ic,jc-1)-TWO*cc(ic,jc))
        else if ( uc(ic,jc-1) ) then
          if (uc(ic,jc-2)) then
            first_der  = HALF*( THREE*cc(ic,jc)-FOUR*cc(ic,jc-1)+cc(ic,jc-2))
            second_der = cc(ic,jc-2) - TWO*cc(ic,jc-1) + cc(ic,jc)
          else
            first_der  =         (cc(ic,jc  )-cc(ic,jc-1))
            second_der = ZERO
          end if
        else if ( uc(ic,jc+1) ) then
          if ( uc(ic,jc+2) ) then
            first_der  = HALF*(-THREE*cc(ic,jc)+FOUR*cc(ic,jc+1)-cc(ic,jc+2))
            second_der = cc(ic,jc+2) - TWO*cc(ic,jc+1) + cc(ic,jc)
          else
            first_der  =        (cc(ic,jc+1)-cc(ic,jc  ))
            second_der = ZERO
          end if
        end if

        j = ir(2)*jc
        do i = lo(1), hi(1)
           do m = 0, ir(2)-1
             ff(i,j+m) = cc(ic,jc) + xloc(m)*first_der + HALF*xloc(m)**2*second_der
           end do
        end do
      end do

    else if ( abs(side) == 2 ) then

      call bl_assert(size(cc,dim=2) .eq. 1, "cc must be one cell wide in ml_interp_bcs_2d")
      jc = loc(2)
      do ic = lo(1)/ir(1), hi(1)/ir(1)
        first_der  = ZERO
        second_der = ZERO
        if ( uc(ic-1,jc) .and. uc(ic+1,jc) ) then
          first_der  = HALF * (cc(ic+1,jc)-cc(ic-1,jc))
          second_der = (cc(ic+1,jc)+cc(ic-1,jc)-TWO*cc(ic,jc))
        else if ( uc(ic-1,jc) ) then
          if ( uc(ic-2,jc) ) then
            first_der  = HALF*( THREE*cc(ic,jc)-FOUR*cc(ic-1,jc)+cc(ic-2,jc))
            second_der = cc(ic-2,jc) - TWO*cc(ic-1,jc) + cc(ic,jc)
          else
            first_der  = (cc(ic,jc)-cc(ic-1,jc))
            second_der = ZERO
          end if
        else if ( uc(ic+1,jc) ) then
          if ( uc(ic+2,jc) ) then
            first_der  = HALF*(-THREE*cc(ic,jc)+FOUR*cc(ic+1,jc)-cc(ic+2,jc))
            second_der = cc(ic+2,jc) - TWO*cc(ic+1,jc) + cc(ic,jc)
          else
            first_der  = (cc(ic+1,jc)-cc(ic,jc  ))
            second_der = ZERO
          end if
        end if

        i = ir(1)*ic
        do j = lo(2),hi(2)
           do m = 0,ir(1)-1
             ff(i+m,j) = cc(ic,jc) + xloc(m)*first_der + HALF*xloc(m)**2*second_der
           end do
        end do
      end do

    end if

  end subroutine ml_interp_bcs_2d

  subroutine ml_interp_bcs_2d_nodal(ff, lof, cc, loc, lo, hi, ir, side)
    integer, intent(in) :: loc(:)
    integer, intent(in) :: lof(:)
    integer, intent(in) :: lo(:), hi(:)
    real (dp_t), intent(inout) :: ff(lof(1):,lof(2):)
    real (dp_t), intent(inout) :: cc(loc(1):,loc(2):)
    integer, intent(in) :: side
    integer, intent(in) :: ir(:)
    integer :: i, j, ic, jc, m
    real (dp_t) :: xloc(0:3)
    real (dp_t) :: first_der

    if (ir(1) == 2) then
       xloc(0) = ZERO
       xloc(1) = HALF
    else
       xloc(0) = ZERO
       xloc(1) = FOURTH
       xloc(2) = HALF
       xloc(3) = THREE*FOURTH
    end if

!   This is the same interpolation as done in IAMR.
    if (side == 1 .or. side == -1) then
       i  = lo(1)
       ic = i / ir(1)
       do jc = lo(2)/ir(2), hi(2)/ir(2)-1
          j = ir(2)*jc
          first_der = cc(ic,jc+1)-cc(ic,jc)
          do m = 0,ir(2)-1
             ff(i,j+m) = cc(ic,jc) + xloc(m)*first_der
          end do
          ff(i,hi(2)) = cc(ic,hi(2)/ir(2))
       end do

    else if (side == 2 .or. side == -2) then

       j  = lo(2)
       jc = j / ir(2)
       do ic = lo(1)/ir(1), hi(1)/ir(1)-1
          i = ir(1)*ic
          first_der = cc(ic+1,jc)-cc(ic,jc)
          do m = 0,ir(1)-1
             ff(i+m,j) = cc(ic,jc) + xloc(m)*first_der
          end do
          ff(hi(1),j) = cc(hi(1)/ir(1),jc)
       end do

    end if

  end subroutine ml_interp_bcs_2d_nodal

  subroutine ml_interp_bcs_3d(ff, lof, cc, uc, loc, hic, lo, hi, ir, side)
    integer, intent(in) :: loc(:), hic(:)
    integer, intent(in) :: lof(:)
    integer, intent(in) :: lo(:), hi(:)
    real (dp_t), intent(inout) :: ff(lof(1):,lof(2):,lof(3):)
    real (dp_t), intent(in) :: cc(loc(1):,loc(2):,loc(3):)
    logical    , intent(in) :: uc(loc(1):,loc(2):,loc(3):)  ! uncovered
    integer, intent(in) :: side
    integer, intent(in) :: ir(:)
    integer :: i, j, k, ic, jc, kc, m, n
    real (dp_t) :: xloc(0:3)
    real (dp_t) :: xder,x2der,yder,y2der,zder,z2der
    real (dp_t) :: xyder,yzder,xzder

    if (ir(1) == 2) then
       xloc(0) = -FOURTH
       xloc(1) =  FOURTH
    else
       xloc(0) = -THREE_EIGHTHS
       xloc(1) = -EIGHTH
       xloc(2) =  EIGHTH
       xloc(3) =  THREE_EIGHTHS
    end if

    xder = 0; x2der = 0; yder = 0; y2der = 0; zder = 0; z2der = 0

!   This is the same interpolation as done in IAMR.
    if (side == 1 .or. side == -1) then
      i  = lo(1)
      ic = loc(1)

      do kc = lo(3)/ir(3), hi(3)/ir(3)
         do jc = lo(2)/ir(2), hi(2)/ir(2)
            if (uc(ic,jc-1,kc) .and. uc(ic,jc+1,kc)) then
               yder  = HALF  * (cc(ic,jc+1,kc) - cc(ic,jc-1,kc))
               y2der = HALF  * (cc(ic,jc+1,kc) + cc(ic,jc-1,kc) - TWO*cc(ic,jc,kc))
            else if (uc(ic,jc+1,kc)) then
               yder  = cc(ic,jc+1,kc) - cc(ic,jc,kc)
               y2der = ZERO
            else if (uc(ic,jc-1,kc)) then
               yder  = cc(ic,jc,kc) - cc(ic,jc-1,kc)
               y2der = ZERO
            end if

            if (uc(ic,jc,kc-1) .and. uc(ic,jc,kc+1)) then
               zder  = HALF  * (cc(ic,jc,kc+1) - cc(ic,jc,kc-1))
               z2der = HALF  * (cc(ic,jc,kc+1) + cc(ic,jc,kc-1) - TWO*cc(ic,jc,kc))
            else if (uc(ic,jc,kc+1)) then
               zder  = cc(ic,jc,kc+1) - cc(ic,jc,kc)
               z2der = ZERO
            else if (uc(ic,jc,kc-1)) then
               zder  = cc(ic,jc,kc) - cc(ic,jc,kc-1)
               z2der = ZERO
            end if

            if (uc(ic,jc+1,kc+1) .and. uc(ic,jc-1,kc-1) .and. uc(ic,jc-1,kc+1) .and. uc(ic,jc+1,kc-1)) then
               yzder = FOURTH * (cc(ic,jc+1,kc+1) + cc(ic,jc-1,kc-1) &
                    -cc(ic,jc-1,kc+1) - cc(ic,jc+1,kc-1) )
            else
               yzder = ZERO
            end if

            j = ir(2)*jc
            k = ir(3)*kc
            do n = 0,ir(3)-1
               do m = 0,ir(2)-1
                  ff(i,j+m,k+n) = cc(ic,jc,kc) + xloc(m)*yder + xloc(n)*zder &
                       + xloc(m)**2*y2der + xloc(n)**2*z2der &
                       + xloc(m)*xloc(n)*yzder
               end do
            end do
         end do
      end do

    else if (side == 2 .or. side == -2) then

      j  = lo(2)
      jc = loc(2)
      do kc = lo(3)/ir(3), hi(3)/ir(3)
         do ic = lo(1)/ir(1), hi(1)/ir(1)
            if (uc(ic-1,jc,kc) .and. uc(ic+1,jc,kc)) then
               xder  = HALF  * (cc(ic+1,jc,kc) - cc(ic-1,jc,kc))
               x2der = HALF  * (cc(ic+1,jc,kc) + cc(ic-1,jc,kc) - TWO*cc(ic,jc,kc))
            else if (uc(ic+1,jc,kc)) then
               xder  = cc(ic+1,jc,kc) - cc(ic,jc,kc)
               x2der = ZERO
            else if (uc(ic-1,jc,kc)) then
               xder  = cc(ic,jc,kc) - cc(ic-1,jc,kc)
               x2der = ZERO
            end if

            if (uc(ic,jc,kc-1) .and. uc(ic,jc,kc+1)) then
               zder  = HALF  * (cc(ic,jc,kc+1) - cc(ic,jc,kc-1))
               z2der = HALF  * (cc(ic,jc,kc+1) + cc(ic,jc,kc-1) - TWO*cc(ic,jc,kc))
            else if (uc(ic,jc,kc+1)) then
               zder  = cc(ic,jc,kc+1) - cc(ic,jc,kc)
               z2der = ZERO
            else if (uc(ic,jc,kc-1)) then
               zder  = cc(ic,jc,kc) - cc(ic,jc,kc-1)
               z2der = ZERO
            end if

            if (uc(ic+1,jc,kc+1) .and. uc(ic-1,jc,kc-1) .and. uc(ic-1,jc,kc+1) .and. uc(ic+1,jc,kc-1)) then
               xzder = FOURTH * (cc(ic+1,jc,kc+1) + cc(ic-1,jc,kc-1) &
                    -cc(ic-1,jc,kc+1) - cc(ic+1,jc,kc-1) )
            else
               xzder = ZERO
            end if

            i = ir(1)*ic
            k = ir(3)*kc
            do n = 0,ir(3)-1
               do m = 0,ir(1)-1
                  ff(i+m,j,k+n) = cc(ic,jc,kc) + xloc(m)*xder + xloc(n)*zder &
                       + xloc(m)**2*x2der + xloc(n)**2*z2der &
                       + xloc(m)*xloc(n)*xzder
               end do
            end do
         end do
      end do

    else if (side == 3 .or. side == -3) then

      k  = lo(3)
      kc = loc(3)
      do jc = lo(2)/ir(2), hi(2)/ir(2)
         do ic = lo(1)/ir(1), hi(1)/ir(1)
            if (uc(ic,jc-1,kc) .and. uc(ic,jc+1,kc)) then
               yder  = HALF  * (cc(ic,jc+1,kc) - cc(ic,jc-1,kc))
               y2der = HALF  * (cc(ic,jc+1,kc) + cc(ic,jc-1,kc) - TWO*cc(ic,jc,kc))
            else if (uc(ic,jc+1,kc)) then
               yder  = cc(ic,jc+1,kc) - cc(ic,jc,kc)
               y2der = ZERO
            else if (uc(ic,jc-1,kc)) then
               yder  = cc(ic,jc,kc) - cc(ic,jc-1,kc)
               y2der = ZERO
            end if

            if (uc(ic-1,jc,kc) .and. uc(ic+1,jc,kc)) then
               xder  = HALF  * (cc(ic+1,jc,kc) - cc(ic-1,jc,kc))
               x2der = HALF  * (cc(ic+1,jc,kc) + cc(ic-1,jc,kc) - TWO*cc(ic,jc,kc))
            else if (uc(ic+1,jc,kc)) then
               xder  = cc(ic+1,jc,kc) - cc(ic,jc,kc)
               x2der = ZERO
            else if (uc(ic-1,jc,kc)) then
               xder  = cc(ic,jc,kc) - cc(ic-1,jc,kc)
               x2der = ZERO
            end if

            if (uc(ic+1,jc+1,kc) .and. uc(ic-1,jc-1,kc) .and. uc(ic+1,jc-1,kc) .and. uc(ic-1,jc+1,kc) ) then
               xyder = FOURTH * (cc(ic+1,jc+1,kc) + cc(ic-1,jc-1,kc) &
                    -cc(ic+1,jc-1,kc) - cc(ic-1,jc+1,kc) )
            else
               xyder = ZERO
            end if

            j = ir(2)*jc
            i = ir(1)*ic
            do n = 0,ir(1)-1
               do m = 0,ir(2)-1
                  ff(i+n,j+m,k) = cc(ic,jc,kc) + xloc(m)*yder + xloc(n)*xder &
                       + xloc(m)**2*y2der + xloc(n)**2*x2der &
                       + xloc(m)*xloc(n)*xyder
               end do
            end do
         end do
      end do

    end if

  end subroutine ml_interp_bcs_3d

  subroutine ml_interp_bcs_3d_nodal(ff, lof, cc, loc, lo, hi, ir, side)
    integer, intent(in) :: loc(:)
    integer, intent(in) :: lof(:)
    integer, intent(in) :: lo(:), hi(:)
    real (dp_t), intent(inout) :: ff(lof(1):,lof(2):,lof(3):)
    real (dp_t), intent(in) :: cc(loc(1):,loc(2):,loc(3):)
    integer, intent(in) :: side
    integer, intent(in) :: ir(:)
    integer :: i, j, k, ic, jc, kc, m, n
    real (dp_t) :: xloc(0:3)
    real (dp_t) :: xder, yder, zder
    real (dp_t) :: xyder,yzder,xzder

    if (ir(1) == 2) then
       xloc(0) = ZERO
       xloc(1) = HALF
    else
       xloc(0) = ZERO
       xloc(1) = FOURTH
       xloc(2) = HALF
       xloc(3) = THREE*FOURTH
    end if

    i = 0; j = 0; k = 0

!   This is the same interpolation as done in IAMR.
    if (side == 1 .or. side == -1) then
      i  = lo(1)
      ic = i / ir(1)

      do kc = lo(3)/ir(3), hi(3)/ir(3)-1
         do jc = lo(2)/ir(2), hi(2)/ir(2)-1

            yder  = cc(ic,jc+1,kc) - cc(ic,jc,kc)
            zder  = cc(ic,jc,kc+1) - cc(ic,jc,kc)
            yzder = cc(ic,jc,kc) + cc(ic,jc+1,kc+1) - cc(ic,jc+1,kc) - cc(ic,jc,kc+1)

            j = ir(2)*jc
            k = ir(3)*kc
            do n = 0,ir(3)-1
               do m = 0,ir(2)-1
                  ff(i,j+m,k+n) = cc(ic,jc,kc) + xloc(m)*yder + xloc(n)*zder &
                       + xloc(m)*xloc(n)*yzder
               end do
            end do
         end do

         jc = hi(2)/ir(2)
         j  = hi(2)
         zder  = cc(ic,jc,kc+1) - cc(ic,jc,kc)
         do m = 0,ir(3)-1
            ff(i,j,k+m) = cc(ic,jc,kc) + xloc(m)*zder
         end do

      end do

      kc = hi(3)/ir(3)
      k  = hi(3)
      do jc = lo(2)/ir(2), hi(2)/ir(2)-1
         yder  = cc(ic,jc+1,kc) - cc(ic,jc,kc)
         do m = 0,ir(2)-1
            ff(i,j+m,k) = cc(ic,jc,kc) + xloc(m)*yder
         end do
      end do

      ff(i,hi(2),hi(3)) = cc(ic,hi(2)/ir(2),hi(3)/ir(3))

    else if (side == 2 .or. side == -2) then

      j  = lo(2)
      jc = j / ir(2)

      do kc = lo(3)/ir(3), hi(3)/ir(3)-1
         do ic = lo(1)/ir(1), hi(1)/ir(1)-1

            xder  = cc(ic+1,jc,kc) - cc(ic,jc,kc)
            zder  = cc(ic,jc,kc+1) - cc(ic,jc,kc)
            xzder = cc(ic,jc,kc) + cc(ic+1,jc,kc+1) - cc(ic+1,jc,kc) - cc(ic,jc,kc+1)

            i = ir(1)*ic
            k = ir(3)*kc
            do n = 0,ir(3)-1
               do m = 0,ir(1)-1
                  ff(i+m,j,k+n) = cc(ic,jc,kc) + xloc(m)*xder + xloc(n)*zder &
                       + xloc(m)*xloc(n)*xzder
               end do
            end do
         end do

         ic = hi(1)/ir(1)
         i  = hi(1)
         zder  = cc(ic,jc,kc+1) - cc(ic,jc,kc)
         do m = 0,ir(3)-1
            ff(i,j,k+m) = cc(ic,jc,kc) + xloc(m)*zder
         end do
      end do

      kc = hi(3)/ir(3)
      k  = hi(3)
      do ic = lo(1)/ir(1), hi(1)/ir(1)-1
         xder  = cc(ic+1,jc,kc) - cc(ic,jc,kc)
         do m = 0,ir(1)-1
            ff(i+m,j,k) = cc(ic,jc,kc) + xloc(m)*xder
         end do
      end do

    else if (side == 3 .or. side == -3) then

      k  = lo(3)
      kc = k / ir(3)

      do jc = lo(2)/ir(2), hi(2)/ir(2)-1
         do ic = lo(1)/ir(1), hi(1)/ir(1)-1

            xder  = cc(ic+1,jc,kc) - cc(ic,jc,kc)
            yder  = cc(ic,jc+1,kc) - cc(ic,jc,kc)
            xyder = cc(ic,jc,kc) + cc(ic+1,jc+1,kc) - cc(ic+1,jc,kc) - cc(ic,jc+1,kc)

            i = ir(1)*ic
            j = ir(2)*jc
            do n = 0,ir(2)-1
               do m = 0,ir(1)-1
                  ff(i+m,j+n,k) = cc(ic,jc,kc) + xloc(m)*xder + xloc(n)*yder &
                       + xloc(m)*xloc(n)*xyder
               end do
            end do
         end do

         ic = hi(1)/ir(1)
         i  = hi(1)
         yder  = cc(ic,jc+1,kc) - cc(ic,jc,kc)
         do m = 0,ir(2)-1
            ff(i,j+m,k) = cc(ic,jc,kc) + xloc(m)*yder
         end do
      end do

      jc = hi(2)/ir(2)
      j  = hi(2)
      do ic = lo(1)/ir(1), hi(1)/ir(1)-1
         xder  = cc(ic+1,jc,kc) - cc(ic,jc,kc)
         do m = 0,ir(1)-1
            ff(i+m,j,k) = cc(ic,jc,kc) + xloc(m)*xder
         end do
      end do

    end if

  end subroutine ml_interp_bcs_3d_nodal

end module ml_prolongation_module
