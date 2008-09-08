module fillpatch_module

  use multifab_module
  use define_bc_module

  implicit none

contains

  subroutine fillpatch(fine, crse, ng, ir, bc_crse, bc_fine, icomp_fine, icomp_crse, &
                       bcomp, nc, no_final_physbc_input, lim_slope_input, lin_limit_input)

    use bc_module
    use setbc_module
    use layout_module
    use interp_module
    use bl_constants_module
    use multifab_physbc_module

    type(multifab), intent(inout)           :: fine
    type(multifab), intent(inout),target    :: crse
    integer       , intent(in   )           :: ng
    integer       , intent(in   )           :: ir(:)
    type(bc_level), intent(in   )           :: bc_crse, bc_fine
    integer       , intent(in   )           :: icomp_fine, icomp_crse, bcomp, nc
    logical       , intent(in   ), optional :: no_final_physbc_input
    logical       , intent(in   ), optional :: lim_slope_input
    logical       , intent(in   ), optional :: lin_limit_input


    integer         :: i, j, dm, local_bc(fine%dim,2,nc), shft(3**fine%dim,fine%dim), cnt
    integer         :: lo_f(3), lo_c(3), hi_f(3), hi_c(3), cslope_lo(3), cslope_hi(3)
    integer         :: n_extra_valid_regions, np
    type(layout)    :: la, fla, tmpla
    type(multifab)  :: cfine, tmpcrse, tmpfine
    type(box)       :: bx, fbx, cbx, fine_box, fdomain, cdomain, bxs(3**fine%dim)
    type(list_box)  :: bl, pbl, pieces, leftover, extra
    type(boxarray)  :: ba, tmpba
    real(kind=dp_t) :: dx(3)
    logical         :: lim_slope, lin_limit, pmask(fine%dim), have_periodic_gcells
    logical         :: no_final_physbc

    type(list_box_node),   pointer     :: bln
    type(box_intersector), pointer     :: bi(:)
    real(kind=dp_t),       allocatable :: fvcx(:), fvcy(:), fvcz(:)
    real(kind=dp_t),       allocatable :: cvcx(:), cvcy(:), cvcz(:)
    integer,               allocatable :: procmap(:)
    real(kind=dp_t),       pointer     :: src(:,:,:,:), dst(:,:,:,:), fp(:,:,:,:)

    type(multifab), target  :: gcrse
    type(multifab), pointer :: pcrse

    type(bl_prof_timer), save :: bpt

    call build(bpt, "fillpatch")

    if (nghost(fine) < ng) call bl_error('fillpatch: fine does NOT have enough ghost cells')
    if (nghost(crse) < ng) call bl_error('fillpatch: crse does NOT have enough ghost cells')

    if ( .not. cell_centered_q(fine) ) call bl_error('fillpatch: fine is NOT cell centered')
    if ( .not. cell_centered_q(crse) ) call bl_error('fillpatch: crse is NOT cell centered')

    dx                   = ONE
    dm                   = crse%dim
    lim_slope            = .true.
    lin_limit            = .false.
    have_periodic_gcells = .false.
    no_final_physbc      = .false.

    if ( present(lim_slope_input)       ) lim_slope       = lim_slope_input
    if ( present(lin_limit_input)       ) lin_limit       = lin_limit_input
    if ( present(no_final_physbc_input) ) no_final_physbc = no_final_physbc_input
    !
    ! Force crse to have good data in ghost cells (only the ng that are needed 
    ! in case has more than ng).
    !
    call fill_boundary(crse, icomp_crse, nc, ng)

    call multifab_physbc(crse,icomp_crse,bcomp,nc,bc_crse)
    !
    ! Build coarsened version of fine such that the fabs @ i are owned by the same CPUs.
    ! We don't try to directly fill anything at fine level outside of the domain.
    !
    fdomain = get_pd(fine%la)

    do i = 1, nboxes(fine)
       !
       ! We don't use get_pbox here as we only want to fill ng ghost cells of 
       ! fine & it may have more ghost cells than that.
       !
       ! Note: we let bl contain empty boxes so that we keep the same number of boxes
       !       in bl as in fine, but in the interpolation later we cycle if it's empty
       bx = intersection(grow(get_ibox(fine,i),ng),fdomain)
       call push_back(bl, bx)
    end do

    pmask(1:dm) = layout_get_pmask(fine%la)

    if ( any(pmask) ) then
       !
       ! Collect additional boxes that contribute to periodically filling fine ghost cells.
       !
       do i = 1, nboxes(fine)
          bx = get_ibox(fine,i)
          call box_periodic_shift(fdomain, bx, fine%nodal, pmask, ng, shft, cnt, bxs)
          if ( cnt > 0 ) have_periodic_gcells = .true.
          do j = 1, cnt
             call push_back(pbl, bxs(j))
          end do
          bln => begin(pbl)
          do while (associated(bln))
             bx =  value(bln)
             bi => layout_get_box_intersector(fine%la, bx)
             do j = 1, size(bi)
                call push_back(pieces, bi(j)%bx)
             end do
             deallocate(bi)
             leftover = boxlist_boxlist_diff(bx, pieces)
             call splice(extra, leftover)
             call destroy(pieces)
             bln => next(bln)
          end do
          call destroy(pbl)
       end do
    end if
    !
    ! n_extra_valid_regions > 0 implies:
    !
    !     Must add additional valid regions to 'fine' to enable the setting
    !     of periodic ghost cells via fill_boundary().  We do this by using
    !     an intermediate multifab 'tmpfine'.
    !
    n_extra_valid_regions = size(extra)
    !
    ! Force first nboxes(fine) in 'tmpfine' to have the same distribution as 'fine'.
    !
    allocate(procmap(1:nboxes(fine)+n_extra_valid_regions))

    procmap(1:nboxes(fine)) = get_proc(fine%la)

    if ( n_extra_valid_regions > 0 ) then
       np = parallel_nprocs()
       do i = 1, n_extra_valid_regions
          !
          ! Distribute extra boxes round-robin.
          !
          procmap(nboxes(fine)+i) = mod(i,np)
       end do
       !
       ! tmpfine looks like fine with extra boxes added.
       !
       do i = 1, nboxes(fine)
          call push_back(pbl, get_ibox(fine,i))
       end do
       bln => begin(extra)
       do while (associated(bln))
          call push_back(pbl, value(bln))
          bln => next(bln)
       end do
       call build(ba, pbl, sort = .false.)
       call destroy(pbl)
       call build(fla, ba, pd = fdomain, pmask = pmask, explicit_mapping = procmap)
       call destroy(ba)
       call build(tmpfine, fla, nc = nc, ng = ng)
       !
       ! Now grow the boxes in extra in preparation for building cfine.
       !
       bln => begin(extra)
       do while (associated(bln))
          call set(bln, intersection(grow(value(bln),ng),fdomain))
          bln => next(bln)
       end do
    end if

    call splice(bl, extra)
    call build(ba, bl, sort = .false.)
    call destroy(bl)
    call boxarray_coarsen(ba, ir)
    call boxarray_grow(ba, 1) ! Grow by one for stencil in lin_cc_interp.
    call build(la, ba, pd = cdomain, pmask = pmask, explicit_mapping = procmap)
    call destroy(ba)
    call build(cfine, la, nc = nc, ng = 0)
    !
    ! Fill cfine from crse.
    ! Got to do it in stages as parallel copy only goes from valid -> valid.
    ! If ng==0, got to build and use a version of crse that has a grow cell.
    !
    if (ng .eq. 0) then
       call build(gcrse, crse%la, nc = nc, ng = 1)
       call copy(gcrse, crse)
       call fill_boundary(gcrse, icomp_crse, nc)
       call multifab_physbc(gcrse,icomp_crse,bcomp,nc,bc_crse)
       pcrse => gcrse
    else
       pcrse => crse
    end if

    do i = 1, nboxes(pcrse)
       call push_back(bl, get_pbox(pcrse,i))
    end do

    call build(tmpba, bl, sort = .false.)
    call destroy(bl)
    call build(tmpla, tmpba, explicit_mapping = get_proc(pcrse%la))
    call destroy(tmpba)
    call build(tmpcrse, tmpla, nc = nc, ng = 0)

    do i = 1, nboxes(pcrse)
       if ( remote(pcrse, i) ) cycle
       src => dataptr(pcrse,   i, icomp_crse, nc)
       dst => dataptr(tmpcrse, i, 1         , nc)
       dst = src
    end do

    if (ng .eq. 0) call destroy(gcrse)

    call copy(cfine, 1, tmpcrse, 1, nc)

    call destroy(tmpcrse)
    call destroy(tmpla)

    cdomain = get_pd(crse%la)

    do i = 1, nboxes(cfine)
       if ( remote(cfine, i) ) cycle

       cbx = get_ibox(cfine,i)

       if ( n_extra_valid_regions > 0 ) then
          fine_box = get_ibox(tmpfine,i)
       else
          fine_box = get_ibox(fine,   i)
       end if

       fbx = intersection(grow(fine_box,ng),fdomain)
       if (empty(fbx)) cycle

       cslope_lo(1:dm) = lwb(grow(cbx, -1))
       cslope_hi(1:dm) = upb(grow(cbx, -1))

       local_bc(:,:,1:nc) = INTERIOR

       if ( cslope_lo(1) == cdomain%lo(1) ) then
          local_bc(1,1,1:nc) = bc_crse%adv_bc_level_array(0,1,1,bcomp:bcomp+nc-1)
       end if
       if ( cslope_hi(1) == cdomain%hi(1) ) then
          local_bc(1,2,1:nc) = bc_crse%adv_bc_level_array(0,1,2,bcomp:bcomp+nc-1)
       end if
       if ( dm > 1 ) then
          if ( cslope_lo(2) == cdomain%lo(2) ) then
             local_bc(2,1,1:nc) = bc_crse%adv_bc_level_array(0,2,1,bcomp:bcomp+nc-1)
          end if
          if ( cslope_hi(2) == cdomain%hi(2) ) then
             local_bc(2,2,1:nc) = bc_crse%adv_bc_level_array(0,2,2,bcomp:bcomp+nc-1)
          end if
       end if
       if ( dm > 2 ) then
          if ( cslope_lo(dm) == cdomain%lo(dm) ) then
             local_bc(dm,1,1:nc) = bc_crse%adv_bc_level_array(0,dm,1,bcomp:bcomp+nc-1)
          end if
          if ( cslope_hi(dm) == cdomain%hi(dm) ) then
             local_bc(dm,2,1:nc) = bc_crse%adv_bc_level_array(0,dm,2,bcomp:bcomp+nc-1)
          end if
       end if

       lo_c(1:dm) = lwb(cbx)
       hi_c(1:dm) = upb(cbx)
       lo_f(1:dm) = lwb(fbx)
       hi_f(1:dm) = upb(fbx)

       allocate(fvcx(lo_f(1):hi_f(1)+1))
       forall (j = lo_f(1):hi_f(1)+1) fvcx(j) = dble(j)
       if ( dm > 1 ) then
          allocate(fvcy(lo_f(2):hi_f(2)+1))
          forall (j = lo_f(2):hi_f(2)+1) fvcy(j) = dble(j) 
          if ( dm > 2 ) then
             allocate(fvcz(lo_f(3):hi_f(3)+1))
             forall (j = lo_f(3):hi_f(3)+1) fvcz(j) = dble(j)
          end if
       end if

       allocate(cvcx(lo_c(1):hi_c(1)+1))
       forall (j = lo_c(1):hi_c(1)+1) cvcx(j) = dble(j) * TWO
       if ( dm > 1 ) then
          allocate(cvcy(lo_c(2):hi_c(2)+1))
          forall (j = lo_c(2):hi_c(2)+1) cvcy(j) = dble(j) * TWO
          if ( dm > 2 ) then
             allocate(cvcz(lo_c(3):hi_c(3)+1))
             forall (j = lo_c(3):hi_c(3)+1) cvcz(j) = dble(j) * TWO
          end if
       end if

       src => dataptr(cfine, i)

       select case (dm)
       case (2)
          allocate(fp(lo_f(1):hi_f(1),lo_f(2):hi_f(2),1:1,1:nc))
          call lin_cc_interp_2d(fp(:,:,1,:), lo_f, src(:,:,1,:), lo_c, ir, local_bc, &
             fvcx, lo_f(1), fvcy, lo_f(2), &
             cvcx, lo_c(1), cvcy, lo_c(2), &
             cslope_lo, cslope_hi, lim_slope, lin_limit)
       case (3)
          allocate(fp(lo_f(1):hi_f(1),lo_f(2):hi_f(2),lo_f(3):hi_f(3),1:nc))
          call lin_cc_interp_3d(fp(:,:,:,:), lo_f, src(:,:,:,:), lo_c, ir, local_bc, &
               fvcx, lo_f(1), fvcy, lo_f(2), fvcz, lo_f(3), &
               cvcx, lo_c(1), cvcy, lo_c(2), cvcz, lo_c(3), &
               cslope_lo, cslope_hi, lim_slope, lin_limit)
       end select

       if ( n_extra_valid_regions > 0 ) then
          dst => dataptr(tmpfine,  i, fbx, 1         , nc)
       else
          dst => dataptr(fine,     i, fbx, icomp_fine, nc)
       end if

       dst = fp

       deallocate(cvcx, fvcx, fp)
       if ( dm > 1 ) deallocate(cvcy, fvcy)
       if ( dm > 2 ) deallocate(cvcz, fvcz)

    end do

    if ( have_periodic_gcells ) then
       if ( n_extra_valid_regions > 0 ) then
          call fill_boundary(tmpfine, 1, nc, ng)
          do i = 1, nboxes(fine)
             if ( remote(fine, i) ) cycle
             bx  =  grow(get_ibox(fine,i), ng)
             dst => dataptr(fine,    i, bx, icomp_fine, nc)
             src => dataptr(tmpfine, i, bx, 1         , nc)
             dst =  src
          end do
       else
          !
          ! We can fill periodic ghost cells simply by calling fill_boundary().
          !
          call fill_boundary(fine, icomp_fine, nc, ng)
       end if
    end if

    if(.not. no_final_physbc) call multifab_physbc(fine, icomp_fine, bcomp, nc, bc_fine)

    call destroy(cfine)
    call destroy(la)

    if ( n_extra_valid_regions > 0 ) then
       call destroy(tmpfine)
       call destroy(fla) 
    end if

    call destroy(bpt)

  end subroutine

end module fillpatch_module
