module fillpatch_module

  use multifab_module
  use define_bc_module

  implicit none

contains

  subroutine fillpatch(fine, crse, ng, ir, bc_crse, bc_fine, icomp_fine, icomp_crse, &
                       bcomp, nc, no_final_physbc_input, lim_slope_input, lin_limit_input, &
                       fill_crse_input, stencil_width_input, fourth_order_input)

    use bc_module
    use layout_module
    use interp_module
    use bl_constants_module
    use multifab_physbc_module

    type(multifab), intent(inout)           :: fine
    type(multifab), intent(inout)           :: crse
    integer       , intent(in   )           :: ng
    integer       , intent(in   )           :: ir(:)
    type(bc_level), intent(in   )           :: bc_crse, bc_fine
    integer       , intent(in   )           :: icomp_fine, icomp_crse, bcomp, nc
    logical       , intent(in   ), optional :: no_final_physbc_input
    logical       , intent(in   ), optional :: lim_slope_input
    logical       , intent(in   ), optional :: lin_limit_input
    logical       , intent(in   ), optional :: fill_crse_input
    integer       , intent(in   ), optional :: stencil_width_input
    logical       , intent(in   ), optional :: fourth_order_input


    integer         :: i, j, dm, local_bc(multifab_get_dim(fine),2,nc), shft(3**multifab_get_dim(fine),multifab_get_dim(fine)), cnt
    integer         :: lo_f(3), lo_c(3), hi_f(3), hi_c(3), cslope_lo(3), cslope_hi(3)
    integer         :: n_extra_valid_regions, np
    type(layout)    :: la, fla, tmpla, fine_la
    type(multifab)  :: cfine, tmpcrse, tmpfine
    type(box)       :: bx, fbx, cbx, fine_box, fdomain, cdomain, bxs(3**multifab_get_dim(fine))
    type(list_box)  :: bl, pbl, pieces, leftover, extra
    type(boxarray)  :: ba, tmpba
    real(kind=dp_t) :: dx(3)
    logical         :: lim_slope, lin_limit, pmask(multifab_get_dim(fine)), have_periodic_gcells
    logical         :: no_final_physbc, fill_crse, nodalflags(multifab_get_dim(fine)), fourth_order
    logical         :: done
    integer         :: stencil_width

    type(list_box_node),   pointer     :: bln
    type(box_intersector), pointer     :: bi(:)
    real(kind=dp_t),       allocatable :: fvcx(:), fvcy(:), fvcz(:)
    real(kind=dp_t),       allocatable :: cvcx(:), cvcy(:), cvcz(:)
    integer,               allocatable :: procmap(:)
    real(kind=dp_t),       pointer     :: src(:,:,:,:), dst(:,:,:,:), fp(:,:,:,:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "fillpatch")

    if ( nghost(fine) < ng ) call bl_error('fillpatch: fine does NOT have enough ghost cells')
    if ( nghost(crse) < ng ) call bl_error('fillpatch: crse does NOT have enough ghost cells')

    if ( .not. cell_centered_q(fine) ) call bl_error('fillpatch: fine is NOT cell centered')
    if ( .not. cell_centered_q(crse) ) call bl_error('fillpatch: crse is NOT cell centered')

    dx                   = ONE
    dm                   = get_dim(crse)
    lim_slope            = .true.
    lin_limit            = .false.
    have_periodic_gcells = .false.
    no_final_physbc      = .false.
    fill_crse            = .true.
    stencil_width        = 1
    fourth_order         = .false.

    if ( present(stencil_width_input)   ) stencil_width   = stencil_width_input
    if ( present(fourth_order_input)    ) fourth_order    = fourth_order_input
    if ( present(lim_slope_input)       ) lim_slope       = lim_slope_input
    if ( present(lin_limit_input)       ) lin_limit       = lin_limit_input
    if ( present(no_final_physbc_input) ) no_final_physbc = no_final_physbc_input
    if ( present(fill_crse_input)       ) fill_crse       = fill_crse_input

    if (fourth_order .and. (stencil_width < 2)) &
       call bl_error('fillpatch: need at least stencil_width = 2 for fourth order interp')

    !
    ! Force crse to have good data in ghost cells (only the ng that are needed 
    ! in case has more than ng).
    !
    if ( fill_crse ) call fill_boundary(crse, icomp_crse, nc, ng=nghost(crse))

    call multifab_physbc(crse,icomp_crse,bcomp,nc,bc_crse)
    !
    ! Build coarsened version of fine such that the fabs @ i are owned by the same CPUs.
    ! We don't try to directly fill anything at fine level outside of the domain.
    !
    fine_la = get_layout(fine)
    fdomain = get_pd(fine_la)

    do i = 1, nboxes(fine%la)
       !
       ! We don't use get_pbox here as we only want to fill ng ghost cells of 
       ! fine & it may have more ghost cells than that.
       !
       ! Note: we let bl contain empty boxes so that we keep the same number of boxes
       !       in bl as in fine, but in the interpolation later we cycle if it's empty.
       !       We do this to keep the mapping of boxes between the various
       !       multifabs consistent with one another and because it means we
       !       don't have to manage changes to the processor map.
       !
       bx = intersection(grow(box_nodalize(get_box(fine%la,i),fine%nodal),ng),fdomain)
       call push_back(bl, bx)
    end do

    pmask(1:dm) = get_pmask(fine_la)

    if ( any(pmask) ) then
       !
       ! Collect additional boxes that contribute to periodically filling fine ghost cells.
       !
       nodalflags = nodal_flags(fine)

       do i = 1, nboxes(fine%la)
          bx = box_nodalize(get_box(fine%la,i),fine%nodal)
          call box_periodic_shift(fdomain, bx, nodalflags, pmask, ng, shft, cnt, bxs)
          if ( cnt > 0 ) have_periodic_gcells = .true.
          do j = 1, cnt
             call push_back(pbl, bxs(j))
          end do
          bln => begin(pbl)
          do while (associated(bln))
             bx =  value(bln)
             bi => layout_get_box_intersector(fine_la, bx)
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
    allocate(procmap(1:nboxes(fine%la)+n_extra_valid_regions))

    procmap(1:nboxes(fine%la)) = get_proc(fine_la)

    if ( n_extra_valid_regions > 0 ) then
       np = parallel_nprocs()
       do i = 1, n_extra_valid_regions
          !
          ! Distribute extra boxes round-robin.
          !
          procmap(nboxes(fine%la)+i) = mod(i,np)
       end do
       !
       ! tmpfine looks like fine with extra boxes added.
       !
       do i = 1, nboxes(fine%la)
          call push_back(pbl, box_nodalize(get_box(fine%la,i),fine%nodal))
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
    ! Grow by stencil_width for stencil in interpolation routine. 
    ! Don't grow empty boxes!
    call boxarray_grow(ba, stencil_width, allow_empty=.true.) 
    call build(la, ba, pd = get_pd(get_layout(crse)), pmask = pmask, explicit_mapping = procmap)
    call destroy(ba)
    call build(cfine, la, nc = nc, ng = 0)
 
    ! Set all of cfine to huge so that we can make sure it gets completely filled below.
    ! Empty boxes aren't setval()'d
    call setval(cfine, Huge(ONE), allow_empty=.true.)
    !
    ! Fill cfine from crse.
    ! Got to do it in stages as parallel copy only goes from valid -> valid.

    do i = 1, nboxes(crse%la)
       call push_back(bl, grow(box_nodalize(get_box(crse%la,i),crse%nodal),crse%ng))
    end do

    call build(tmpba, bl, sort = .false.)
    call destroy(bl)
    call build(tmpla, tmpba, boxarray_bbox(tmpba), explicit_mapping = get_proc(get_layout(crse)))
    call destroy(tmpba)
    call build(tmpcrse, tmpla, nc = nc, ng = 0)

    !$OMP PARALLEL DO PRIVATE(i,src,dst)
    do i = 1, nfabs(crse)
       src => dataptr(crse,    i, icomp_crse, nc)
       dst => dataptr(tmpcrse, i, 1         , nc)
       ! dst = src failed using Intel compiler 9.1.043
       call cpy_d(dst,src)
    end do
    !$OMP END PARALLEL DO

    call copy(cfine, 1, tmpcrse, 1, nc)

    if (multifab_max(cfine, allow_empty=.true., local=.true.) .gt. Huge(ONE)-ONE) then
       if (multifab_max(tmpcrse, allow_empty=.true., local=.true.) .gt. Huge(ONE)-ONE) then
          call bl_error('fillpatch: tmpcrse greater than Huge-1 before trying to fill cfine')
       else
          call bl_error('fillpatch: cfine was not completely filled by tmpcrse' // &
               ' (likely because grids are not properly nested)')
       end if
    end if

    call destroy(tmpcrse)
    call destroy(tmpla)

    cdomain = get_pd(get_layout(crse))

    !$OMP PARALLEL DO PRIVATE(i,j,cbx,fine_box,fbx,cslope_lo,cslope_hi,local_bc,lo_c,hi_c,lo_f,hi_f) &
    !$OMP PRIVATE(fvcx,fvcy,fvcz,cvcx,cvcy,cvcz,src,fp,dst)
    do i = 1, nfabs(cfine)

       cbx = get_ibox(cfine,i)

       if ( n_extra_valid_regions > 0 ) then
          fine_box = get_ibox(tmpfine,i)
       else
          fine_box = get_ibox(fine,   i)
       end if

       fbx = intersection(grow(fine_box,ng),fdomain)
       if (empty(fbx)) cycle

       cslope_lo(1:dm) = lwb(grow(cbx,-stencil_width))
       cslope_hi(1:dm) = upb(grow(cbx,-stencil_width))

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
       forall (j = lo_f(1):hi_f(1)+1) fvcx(j) = j
       if ( dm > 1 ) then
          allocate(fvcy(lo_f(2):hi_f(2)+1))
          forall (j = lo_f(2):hi_f(2)+1) fvcy(j) = j 
          if ( dm > 2 ) then
             allocate(fvcz(lo_f(3):hi_f(3)+1))
             forall (j = lo_f(3):hi_f(3)+1) fvcz(j) = j
          end if
       end if

       allocate(cvcx(lo_c(1):hi_c(1)+1))
       forall (j = lo_c(1):hi_c(1)+1) cvcx(j) = j * TWO
       if ( dm > 1 ) then
          allocate(cvcy(lo_c(2):hi_c(2)+1))
          forall (j = lo_c(2):hi_c(2)+1) cvcy(j) = j * TWO
          if ( dm > 2 ) then
             allocate(cvcz(lo_c(3):hi_c(3)+1))
             forall (j = lo_c(3):hi_c(3)+1) cvcz(j) = j * TWO
          end if
       end if

       src => dataptr(cfine, i)

       select case (dm)
       case (1)
          allocate(fp(lo_f(1):hi_f(1),1:1,1:1,1:nc))
          if (fourth_order) then
             call bl_error('fillpatch: fourth_order_interp not implemented in 1d')
          else
             call lin_cc_interp_1d(fp(:,1,1,:), lo_f, src(:,1,1,:), lo_c, ir, local_bc, &
                fvcx, lo_f(1), cvcx, lo_c(1), &
                cslope_lo, cslope_hi, lim_slope, lin_limit)
          endif
       case (2)
          allocate(fp(lo_f(1):hi_f(1),lo_f(2):hi_f(2),1:1,1:nc))
          if (fourth_order) then
             call fourth_order_interp_2d(fp(:,:,1,:), lo_f, src(:,:,1,:), lo_c, ir, &
                  cslope_lo, cslope_hi)
          else
             call lin_cc_interp_2d(fp(:,:,1,:), lo_f, src(:,:,1,:), lo_c, ir, local_bc, &
                fvcx, lo_f(1), fvcy, lo_f(2), &
                cvcx, lo_c(1), cvcy, lo_c(2), &
                cslope_lo, cslope_hi, lim_slope, lin_limit)
          endif
       case (3)
          allocate(fp(lo_f(1):hi_f(1),lo_f(2):hi_f(2),lo_f(3):hi_f(3),1:nc))
          if (fourth_order) then
             call bl_error('fillpatch: fourth_order_interp not implemented in 3d')
          else
             call lin_cc_interp_3d(fp(:,:,:,:), lo_f, src(:,:,:,:), lo_c, ir, local_bc, &
                  fvcx, lo_f(1), fvcy, lo_f(2), fvcz, lo_f(3), &
                  cvcx, lo_c(1), cvcy, lo_c(2), cvcz, lo_c(3), &
                  cslope_lo, cslope_hi, lim_slope, lin_limit)
          endif
       end select

       if ( n_extra_valid_regions > 0 ) then
          dst => dataptr(tmpfine,  i, fbx, 1         , nc)
       else
          dst => dataptr(fine,     i, fbx, icomp_fine, nc)
       end if

       ! dst = fp failed using Intel 9.1.043
       call cpy_d(dst,fp)

       deallocate(cvcx, fvcx, fp)
       if ( dm > 1 ) deallocate(cvcy, fvcy)
       if ( dm > 2 ) deallocate(cvcz, fvcz)

    end do
    !$OMP END PARALLEL DO

    if ( have_periodic_gcells ) then
       if ( n_extra_valid_regions > 0 ) then
          call fill_boundary(tmpfine, 1, nc, ng)
          !$OMP PARALLEL DO PRIVATE(i,bx,dst,src)
          do i = 1, nfabs(fine)
             bx  =  grow(get_ibox(fine,i), ng)
             dst => dataptr(fine,    i, bx, icomp_fine, nc)
             src => dataptr(tmpfine, i, bx, 1         , nc)
             call cpy_d(dst,src)
          end do
          !$OMP END PARALLEL DO
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

  end subroutine fillpatch

  !
  ! This version of fillpatch takes both a crse_old and a crse_new
  !     and does interpolation in time as well as in space.  The coefficient alpha
  !     is used to define:   fine = interp (alpha*crse_old + (1-alpha)*crse_new )
  !
  subroutine fillpatch_t(fine, crse_old, crse_new, alpha, &
                         ng, ir, bc_crse, bc_fine, icomp_fine, icomp_crse, &
                         bcomp, nc, no_final_physbc_input, lim_slope_input, lin_limit_input, &
                         fill_crse_input, stencil_width_input, fourth_order_input)

    use bc_module
    use layout_module
    use interp_module
    use bl_constants_module
    use multifab_physbc_module

    type(multifab), intent(inout)           :: fine
    type(multifab), intent(inout)           :: crse_old
    type(multifab), intent(inout)           :: crse_new
    real(kind=dp_t), intent(in   )          :: alpha
    integer       , intent(in   )           :: ng
    integer       , intent(in   )           :: ir(:)
    type(bc_level), intent(in   )           :: bc_crse, bc_fine
    integer       , intent(in   )           :: icomp_fine, icomp_crse, bcomp, nc
    logical       , intent(in   ), optional :: no_final_physbc_input
    logical       , intent(in   ), optional :: lim_slope_input
    logical       , intent(in   ), optional :: lin_limit_input
    logical       , intent(in   ), optional :: fill_crse_input
    integer       , intent(in   ), optional :: stencil_width_input
    logical       , intent(in   ), optional :: fourth_order_input


    integer         :: i, j, dm, local_bc(multifab_get_dim(fine),2,nc), shft(3**multifab_get_dim(fine),multifab_get_dim(fine)), cnt
    integer         :: lo_f(3), lo_c(3), hi_f(3), hi_c(3), cslope_lo(3), cslope_hi(3)
    integer         :: n_extra_valid_regions, np
    type(layout)    :: la, fla, tmpla, fine_la
    type(multifab)  :: cfine, tmpcrse, tmpfine
    type(box)       :: bx, fbx, cbx, fine_box, fdomain, cdomain, bxs(3**multifab_get_dim(fine))
    type(list_box)  :: bl, pbl, pieces, leftover, extra
    type(boxarray)  :: ba, tmpba
    real(kind=dp_t) :: dx(3)
    logical         :: lim_slope, lin_limit, pmask(multifab_get_dim(fine)), have_periodic_gcells
    logical         :: no_final_physbc, fill_crse, nodalflags(multifab_get_dim(fine)), fourth_order
    integer         :: stencil_width
    integer         :: ng_crse

    type(list_box_node),   pointer     :: bln
    type(box_intersector), pointer     :: bi(:)
    real(kind=dp_t),       allocatable :: fvcx(:), fvcy(:), fvcz(:)
    real(kind=dp_t),       allocatable :: cvcx(:), cvcy(:), cvcz(:)
    integer,               allocatable :: procmap(:)
    real(kind=dp_t),       pointer     :: src_old(:,:,:,:), src_new(:,:,:,:)
    real(kind=dp_t),       pointer     :: src(:,:,:,:), dst(:,:,:,:), fp(:,:,:,:)

    real(kind=dp_t) :: omalpha

    type(bl_prof_timer), save :: bpt

    call build(bpt, "fillpatch_t")

    ng_crse = nghost(crse_old)
    if (nghost(crse_new) .ne. ng_crse) &
       call bl_error('fillpatch: nghost(crse_old) does not equal nghost(crse_new)')

    if ( nghost(fine) < ng ) call bl_error('fillpatch: fine does NOT have enough ghost cells')
    if ( ng_crse      < ng ) call bl_error('fillpatch: crse_old does NOT have enough ghost cells')

    if ( .not. cell_centered_q(fine)     ) call bl_error('fillpatch: fine is NOT cell centered')
    if ( .not. cell_centered_q(crse_old) ) call bl_error('fillpatch: crse_old is NOT cell centered')
    if ( .not. cell_centered_q(crse_new) ) call bl_error('fillpatch: crse_new is NOT cell centered')

    dx                   = ONE
    dm                   = get_dim(crse_old)
    lim_slope            = .true.
    lin_limit            = .false.
    have_periodic_gcells = .false.
    no_final_physbc      = .false.
    fill_crse            = .true.
    stencil_width        = 1
    fourth_order         = .false.

    if ( present(stencil_width_input)   ) stencil_width   = stencil_width_input
    if ( present(fourth_order_input)    ) fourth_order    = fourth_order_input
    if ( present(lim_slope_input)       ) lim_slope       = lim_slope_input
    if ( present(lin_limit_input)       ) lin_limit       = lin_limit_input
    if ( present(no_final_physbc_input) ) no_final_physbc = no_final_physbc_input
    if ( present(fill_crse_input)       ) fill_crse       = fill_crse_input

    if (fourth_order .and. (stencil_width < 2)) &
       call bl_error('fillpatch: need at least stencil_width = 2 for fourth order interp')

    !
    ! Force crse to have good data in ghost cells (only the ng that are needed 
    ! in case has more than ng).
    !
    if ( fill_crse ) then
        call fill_boundary(crse_old, icomp_crse, nc, ng=ng_crse)
        call fill_boundary(crse_new, icomp_crse, nc, ng=ng_crse)
    end if

    call multifab_physbc(crse_old,icomp_crse,bcomp,nc,bc_crse)
    call multifab_physbc(crse_new,icomp_crse,bcomp,nc,bc_crse)
    !
    ! Build coarsened version of fine such that the fabs @ i are owned by the same CPUs.
    ! We don't try to directly fill anything at fine level outside of the domain.
    !
    fine_la = get_layout(fine)
    fdomain = get_pd(fine_la)

    do i = 1, nboxes(fine%la)
       !
       ! We don't use get_pbox here as we only want to fill ng ghost cells of 
       ! fine & it may have more ghost cells than that.
       !
       ! Note: we let bl contain empty boxes so that we keep the same number of boxes
       !       in bl as in fine, but in the interpolation later we cycle if it's empty.
       !       We do this to keep the mapping of boxes between the various
       !       multifabs consistent with one another and because it means we
       !       don't have to manage changes to the processor map.
       !
       bx = intersection(grow(box_nodalize(get_box(fine%la,i),fine%nodal),ng),fdomain)
       call push_back(bl, bx)
    end do

    pmask(1:dm) = get_pmask(fine_la)

    if ( any(pmask) ) then
       !
       ! Collect additional boxes that contribute to periodically filling fine ghost cells.
       !
       nodalflags = nodal_flags(fine)

       do i = 1, nboxes(fine%la)
          bx = box_nodalize(get_box(fine%la,i),fine%nodal)
          call box_periodic_shift(fdomain, bx, nodalflags, pmask, ng, shft, cnt, bxs)
          if ( cnt > 0 ) have_periodic_gcells = .true.
          do j = 1, cnt
             call push_back(pbl, bxs(j))
          end do
          bln => begin(pbl)
          do while (associated(bln))
             bx =  value(bln)
             bi => layout_get_box_intersector(fine_la, bx)
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
    allocate(procmap(1:nboxes(fine%la)+n_extra_valid_regions))

    procmap(1:nboxes(fine%la)) = get_proc(fine_la)

    if ( n_extra_valid_regions > 0 ) then
       np = parallel_nprocs()
       do i = 1, n_extra_valid_regions
          !
          ! Distribute extra boxes round-robin.
          !
          procmap(nboxes(fine%la)+i) = mod(i,np)
       end do
       !
       ! tmpfine looks like fine with extra boxes added.
       !
       do i = 1, nboxes(fine%la)
          call push_back(pbl, box_nodalize(get_box(fine%la,i),fine%nodal))
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
    ! Grow by stencil_width for stencil in interpolation routine. 
    ! Don't grow empty boxes!
    call boxarray_grow(ba, stencil_width, allow_empty=.true.) 
    call build(la, ba, pd = get_pd(get_layout(crse_old)), pmask = pmask, explicit_mapping = procmap)
    call destroy(ba)
    call build(cfine, la, nc = nc, ng = 0)
 
    ! Set all of cfine to huge so that we can make sure it gets completely filled below.
    ! Empty boxes aren't setval()'d
    call setval(cfine, Huge(ONE), allow_empty=.true.)
    !
    ! Fill cfine from crse.
    ! Got to do it in stages as parallel copy only goes from valid -> valid.

    do i = 1, nboxes(crse_old%la)
       call push_back(bl, grow(box_nodalize(get_box(crse_old%la,i),crse_old%nodal),crse_old%ng))
    end do

    call build(tmpba, bl, sort = .false.)
    call destroy(bl)
    call build(tmpla, tmpba, boxarray_bbox(tmpba), explicit_mapping = get_proc(get_layout(crse_old)))
    call destroy(tmpba)
    call build(tmpcrse, tmpla, nc = nc, ng = 0)

    ! tmpcrse = alpha * crse_old + (ONE-alpha) * crse_new
    omalpha = ONE - alpha 

    !$OMP PARALLEL DO PRIVATE(i)
    do i = 1, nfabs(crse_old)
       src_old => dataptr(crse_old,i,icomp_crse, nc)
       src_new => dataptr(crse_new,i,icomp_crse, nc)
       dst     => dataptr(tmpcrse, i, 1         , nc)
       dst = alpha * src_old + omalpha * src_new 
    end do
    !$OMP END PARALLEL DO

    call copy(cfine, 1, tmpcrse, 1, nc)

    if (multifab_max(cfine, allow_empty=.true., local=.true.) .ge. Huge(ONE)-ONE) then
       if (multifab_max(tmpcrse, allow_empty=.true., local=.true.) .ge. Huge(ONE)-ONE) then
          call bl_error('fillpatch: tmpcrse greater than Huge-1 before trying to fill cfine')
       else
          call bl_error('fillpatch: cfine was not completely filled by tmpcrse' // &
               ' (likely because grids are not properly nested)')
       end if
    end if

    call destroy(tmpcrse)
    call destroy(tmpla)

    cdomain = get_pd(get_layout(crse_old))

    !$OMP PARALLEL DO PRIVATE(i,j,cbx,fine_box,fbx,cslope_lo,cslope_hi,local_bc,lo_c,hi_c,lo_f,hi_f) &
    !$OMP PRIVATE(fvcx,fvcy,fvcz,cvcx,cvcy,cvcz,src,fp,dst)
    do i = 1, nfabs(cfine)

       cbx = get_ibox(cfine,i)

       if ( n_extra_valid_regions > 0 ) then
          fine_box = get_ibox(tmpfine,i)
       else
          fine_box = get_ibox(fine,   i)
       end if

       fbx = intersection(grow(fine_box,ng),fdomain)
       if (empty(fbx)) cycle

       cslope_lo(1:dm) = lwb(grow(cbx,-stencil_width))
       cslope_hi(1:dm) = upb(grow(cbx,-stencil_width))

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
       forall (j = lo_f(1):hi_f(1)+1) fvcx(j) = j
       if ( dm > 1 ) then
          allocate(fvcy(lo_f(2):hi_f(2)+1))
          forall (j = lo_f(2):hi_f(2)+1) fvcy(j) = j 
          if ( dm > 2 ) then
             allocate(fvcz(lo_f(3):hi_f(3)+1))
             forall (j = lo_f(3):hi_f(3)+1) fvcz(j) = j
          end if
       end if

       allocate(cvcx(lo_c(1):hi_c(1)+1))
       forall (j = lo_c(1):hi_c(1)+1) cvcx(j) = j * TWO
       if ( dm > 1 ) then
          allocate(cvcy(lo_c(2):hi_c(2)+1))
          forall (j = lo_c(2):hi_c(2)+1) cvcy(j) = j * TWO
          if ( dm > 2 ) then
             allocate(cvcz(lo_c(3):hi_c(3)+1))
             forall (j = lo_c(3):hi_c(3)+1) cvcz(j) = j * TWO
          end if
       end if

       src => dataptr(cfine, i)

       select case (dm)
       case (1)
          allocate(fp(lo_f(1):hi_f(1),1:1,1:1,1:nc))
          if (fourth_order) then
             call bl_error('fillpatch: fourth_order_interp not implemented in 1d')
          else
             call lin_cc_interp_1d(fp(:,1,1,:), lo_f, src(:,1,1,:), lo_c, ir, local_bc, &
                fvcx, lo_f(1), cvcx, lo_c(1), &
                cslope_lo, cslope_hi, lim_slope, lin_limit)
          endif
       case (2)
          allocate(fp(lo_f(1):hi_f(1),lo_f(2):hi_f(2),1:1,1:nc))
          if (fourth_order) then
             call fourth_order_interp_2d(fp(:,:,1,:), lo_f, src(:,:,1,:), lo_c, ir, &
                  cslope_lo, cslope_hi)
          else
             call lin_cc_interp_2d(fp(:,:,1,:), lo_f, src(:,:,1,:), lo_c, ir, local_bc, &
                fvcx, lo_f(1), fvcy, lo_f(2), &
                cvcx, lo_c(1), cvcy, lo_c(2), &
                cslope_lo, cslope_hi, lim_slope, lin_limit)
          endif
       case (3)
          allocate(fp(lo_f(1):hi_f(1),lo_f(2):hi_f(2),lo_f(3):hi_f(3),1:nc))
          if (fourth_order) then
             call bl_error('fillpatch: fourth_order_interp not implemented in 3d')
          else
             call lin_cc_interp_3d(fp(:,:,:,:), lo_f, src(:,:,:,:), lo_c, ir, local_bc, &
                  fvcx, lo_f(1), fvcy, lo_f(2), fvcz, lo_f(3), &
                  cvcx, lo_c(1), cvcy, lo_c(2), cvcz, lo_c(3), &
                  cslope_lo, cslope_hi, lim_slope, lin_limit)
          endif
       end select

       if ( n_extra_valid_regions > 0 ) then
          dst => dataptr(tmpfine,  i, fbx, 1         , nc)
       else
          dst => dataptr(fine,     i, fbx, icomp_fine, nc)
       end if

       ! dst = fp failed using Intel 9.1.043
       call cpy_d(dst,fp)

       deallocate(cvcx, fvcx, fp)
       if ( dm > 1 ) deallocate(cvcy, fvcy)
       if ( dm > 2 ) deallocate(cvcz, fvcz)

    end do
    !$OMP END PARALLEL DO

    if ( have_periodic_gcells ) then
       if ( n_extra_valid_regions > 0 ) then
          call fill_boundary(tmpfine, 1, nc, ng)
          !$OMP PARALLEL DO PRIVATE(i,bx,dst,src)
          do i = 1, nfabs(fine)
             bx  =  grow(get_ibox(fine,i), ng)
             dst => dataptr(fine,    i, bx, icomp_fine, nc)
             src => dataptr(tmpfine, i, bx, 1         , nc)
             call cpy_d(dst,src)
          end do
          !$OMP END PARALLEL DO
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

  end subroutine fillpatch_t

end module fillpatch_module
