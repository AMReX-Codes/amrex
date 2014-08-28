module multifab_fill_ghost_module

  use fillpatch_module

  implicit none

contains

  subroutine multifab_fill_ghost_cells(fine,crse,ng,ir,bc_crse,bc_fine,icomp,bcomp,nc, &
                                       fill_crse_input,stencil_width_input,fourth_order_input)

    use layout_module
    use bl_constants_module
    use multifab_physbc_module

    type(multifab), intent(inout) :: fine
    type(multifab), intent(inout) :: crse
    integer       , intent(in   ) :: ng
    integer       , intent(in   ) :: ir(:)
    type(bc_level), intent(in   ) :: bc_crse,bc_fine
    integer       , intent(in   ) :: icomp,bcomp,nc
    logical       , intent(in   ), optional :: fill_crse_input
    integer       , intent(in   ), optional :: stencil_width_input
    logical       , intent(in   ), optional :: fourth_order_input

    integer         :: i, j
    type(multifab)  :: ghost, tmpfine
    type(box)       :: bx
    type(boxarray)  :: ba
    type(list_box)  :: bl
    type(layout)    :: la, tmpla, fine_la
    real(kind=dp_t) :: dx(3)
    type(fgassoc)   :: fgasc
    logical         :: fill_crse
    integer         :: stencil_width
    logical         :: fourth_order

    real(kind=dp_t),     pointer :: src(:,:,:,:), dst(:,:,:,:)

    type(bl_prof_timer), save :: bpt

    if (ng == 0) return

    call build(bpt, "mf_fill_ghost_cells")

    fill_crse    = .true.
    fourth_order = .false.

    if ( present(fill_crse_input    ) )    fill_crse = fill_crse_input
    if ( present(fourth_order_input ) ) fourth_order = fourth_order_input

    if ( nghost(fine) <  ng          ) &
         call bl_error('multifab_fill_ghost_cells: fine does NOT have enough ghost cells')

    if ( .not. cell_centered_q(fine) ) &
         call bl_error('fillpatch: fine is NOT cell centered')

    if (present(stencil_width_input)) then
       if ( fourth_order ) then
          if ( stencil_width_input < 2) &
            call bl_error('fillpatch: fourth_order but stencil_width < 2')
       end if
    else
       if ( fourth_order) then
          stencil_width = 2
       else
          stencil_width = 1
       end if
    end if

    fine_la = get_layout(fine)
    !
    ! Grab the cached boxarray of all ghost cells not covered by valid region.
    !
    fgasc = layout_fgassoc(fine_la, ng)
    !
    ! Now fillpatch a temporary multifab on those ghost cells.
    !
    ! We ask for a grow cell so we get wide enough strips to enable HOEXTRAP.
    !
    call build(la, fgasc%ba, get_pd(fine_la), get_pmask(fine_la))

    ! Don't need to make any ghost cells.
    call build(ghost, la, nc, ng = 0)

    ! Don't ask fillpatch to fill any ghost cells.
    call fillpatch(ghost, crse, 0, ir, bc_crse, bc_fine, 1, icomp, bcomp, nc, &
                   no_final_physbc_input=.true., fill_crse_input=fill_crse, &
                   stencil_width_input=stencil_width, fourth_order_input=fourth_order)

    !
    ! Copy fillpatch()d ghost cells to fine.
    ! We want to copy the valid region of ghost -> valid + ghost region of fine.
    ! Got to do it in two stages since copy()s only go from valid -> valid.
    !
    do i = 1, nboxes(fine%la)
       call push_back(bl, grow(box_nodalize(get_box(fine%la,i),fine%nodal),ng))
    end do

    call build(ba, bl, sort = .false.)
    call destroy(bl)
    call build(tmpla, ba, get_pd(fine_la), get_pmask(fine_la), explicit_mapping = get_proc(fine_la))
    call destroy(ba)
    call build(tmpfine, tmpla, nc = nc, ng = 0)
    call setval(tmpfine, 0.0_dp_t, all = .true. )

    call copy(tmpfine, 1, ghost, 1, nc)  ! parallel copy

    !$OMP PARALLEL DO PRIVATE(i,ba,j,bx,dst,src)
    do i = 1, nfabs(fine)
       call boxarray_box_diff(ba, get_ibox(tmpfine,i), get_ibox(fine,i))
       do j = 1, nboxes(ba)
          bx  =  get_box(ba,j)
          dst => dataptr(fine,    i, bx, icomp, nc)
          src => dataptr(tmpfine, i, bx, 1    , nc)
          call cpy_d(dst,src)
       end do
       call destroy(ba)
    end do
    !$OMP END PARALLEL DO
    !
    ! Finish up.
    !
    call fill_boundary(fine, icomp, nc, ng)

    dx = ONE

    call multifab_physbc(fine, icomp, bcomp, nc, bc_fine)

    call destroy(ghost)
    call destroy(tmpfine)

    call destroy(la)
    call destroy(tmpla)

    call destroy(bpt)

  end subroutine multifab_fill_ghost_cells

  !
  ! This version of multifab_fill_ghost_cells takes both a crse_old and a crse_new
  !     and does interpolation in time as well as in space.  The coefficient alpha
  !     is used to define:   fine = interp (alpha*crse_old + (1-alpha)*crse_new )
  !
  subroutine multifab_fill_ghost_cells_t(fine,crse_old,crse_new,alpha,ng, &
                                         ir,bc_crse,bc_fine,icomp,bcomp,nc, &
                                         fill_crse_input,stencil_width_input,fourth_order_input)

    use layout_module
    use bl_constants_module
    use multifab_physbc_module

    type(multifab) , intent(inout) :: fine
    type(multifab) , intent(inout) :: crse_old
    type(multifab) , intent(inout) :: crse_new
    integer        , intent(in   ) :: ng
    integer        , intent(in   ) :: ir(:)
    real(kind=dp_t), intent(in   ) :: alpha
    type(bc_level) , intent(in   ) :: bc_crse,bc_fine
    integer        , intent(in   ) :: icomp,bcomp,nc
    logical        , intent(in   ), optional :: fill_crse_input
    integer        , intent(in   ), optional :: stencil_width_input
    logical        , intent(in   ), optional :: fourth_order_input

    integer         :: i, j
    type(multifab)  :: ghost, tmpfine
    type(box)       :: bx
    type(boxarray)  :: ba
    type(list_box)  :: bl
    type(layout)    :: la, tmpla, fine_la
    real(kind=dp_t) :: dx(3)
    type(fgassoc)   :: fgasc
    logical         :: fill_crse
    integer         :: stencil_width
    logical         :: fourth_order

    real(kind=dp_t),     pointer :: src(:,:,:,:), dst(:,:,:,:)

    type(bl_prof_timer), save :: bpt

    if (ng == 0) return

    call build(bpt, "mf_fill_ghost_cells")

    fill_crse    = .true.
    fourth_order = .false.

    if ( present(fill_crse_input    ) )    fill_crse = fill_crse_input
    if ( present(fourth_order_input ) ) fourth_order = fourth_order_input

    if ( nghost(fine) <  ng          ) &
         call bl_error('multifab_fill_ghost_cells: fine does NOT have enough ghost cells')

    if ( .not. cell_centered_q(fine) ) &
         call bl_error('fillpatch: fine is NOT cell centered')

    if (present(stencil_width_input)) then
       if ( fourth_order ) then
          if ( stencil_width_input < 2) &
            call bl_error('fillpatch: fourth_order but stencil_width < 2')
       end if
    else
       if ( fourth_order) then
          stencil_width = 2
       else
          stencil_width = 1
       end if
    end if

    fine_la = get_layout(fine)
    !
    ! Grab the cached boxarray of all ghost cells not covered by valid region.
    !
    fgasc = layout_fgassoc(fine_la, ng)
    !
    ! Now fillpatch a temporary multifab on those ghost cells.
    !
    ! We ask for a grow cell so we get wide enough strips to enable HOEXTRAP.
    !
    call build(la, fgasc%ba, get_pd(fine_la), get_pmask(fine_la))

    ! Don't need to make any ghost cells.
    call build(ghost, la, nc, ng = 0)

    ! Don't ask fillpatch to fill any ghost cells.
    call fillpatch_t(ghost, crse_old, crse_new, alpha, &
                     0, ir, bc_crse, bc_fine, 1, icomp, bcomp, nc, &
                     no_final_physbc_input=.true., fill_crse_input=fill_crse, &
                     stencil_width_input=stencil_width, fourth_order_input=fourth_order)

    !
    ! Copy fillpatch()d ghost cells to fine.
    ! We want to copy the valid region of ghost -> valid + ghost region of fine.
    ! Got to do it in two stages since copy()s only go from valid -> valid.
    !
    do i = 1, nboxes(fine%la)
       call push_back(bl, grow(box_nodalize(get_box(fine%la,i),fine%nodal),ng))
    end do

    call build(ba, bl, sort = .false.)
    call destroy(bl)
    call build(tmpla, ba, get_pd(fine_la), get_pmask(fine_la), explicit_mapping = get_proc(fine_la))
    call destroy(ba)
    call build(tmpfine, tmpla, nc = nc, ng = 0)
    call setval(tmpfine, 0.0_dp_t, all = .true. )

    call copy(tmpfine, 1, ghost, 1, nc)  ! parallel copy

    !$OMP PARALLEL DO PRIVATE(i,ba,j,bx,dst,src)
    do i = 1, nfabs(fine)
       call boxarray_box_diff(ba, get_ibox(tmpfine,i), get_ibox(fine,i))
       do j = 1, nboxes(ba)
          bx  =  get_box(ba,j)
          dst => dataptr(fine,    i, bx, icomp, nc)
          src => dataptr(tmpfine, i, bx, 1    , nc)
          call cpy_d(dst,src)
       end do
       call destroy(ba)
    end do
    !$OMP END PARALLEL DO
    !
    ! Finish up.
    !
    call fill_boundary(fine, icomp, nc, ng)

    dx = ONE

    call multifab_physbc(fine, icomp, bcomp, nc, bc_fine)

    call destroy(ghost)
    call destroy(tmpfine)

    call destroy(la)
    call destroy(tmpla)

    call destroy(bpt)

  end subroutine multifab_fill_ghost_cells_t

end module multifab_fill_ghost_module
