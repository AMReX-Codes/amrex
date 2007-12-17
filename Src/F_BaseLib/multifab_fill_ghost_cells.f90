module multifab_fill_ghost_module

  use layout_module
  use fab_module
  use bl_mem_stat_module
  use multifab_module
  use bc_module
  use setbc_module
  use interp_module
  use fillpatch_module
  use bl_prof_module

  implicit none

contains

  subroutine multifab_fill_ghost_cells(fine,crse,ng,ir,bc_crse,bc_fine,icomp,bcomp,nc)

    type(multifab), intent(inout) :: fine
    type(multifab), intent(inout) :: crse
    integer       , intent(in   ) :: ng
    integer       , intent(in   ) :: ir(:)
    type(bc_level), intent(in   ) :: bc_crse,bc_fine
    integer       , intent(in   ) :: icomp,bcomp,nc

    integer         :: i
    type(multifab)  :: ghost, tmpfine
    type(box)       :: bx
    type(boxarray)  :: ba
    type(list_box)  :: bl
    type(layout)    :: la, tmpla
    real(kind=dp_t) :: dx(3)
    type(fgassoc)   :: fgasc

    real(kind=dp_t), pointer :: src(:,:,:,:), dst(:,:,:,:)

    type(bl_prof_timer), save :: bpt

    if (ng == 0) return

    call build(bpt, "mf_fill_ghost_cells")

    if ( nghost(fine) <  ng          ) &
         call bl_error('fillpatch: fine does NOT have enough ghost cells')

    if ( .not. cell_centered_q(fine) ) &
         call bl_error('fillpatch: fine is NOT cell centered')
    !
    ! Grab the cached boxarray of all ghost cells not covered by valid region.
    !
    fgasc = layout_fgassoc(fine%la, ng)
    !
    ! Now fillpatch a temporary multifab on those ghost cells.
    !
    ! We ask for a grow cell so we get wide enough strips to enable HOEXTRAP.
    !
    call build(la, fgasc%ba, get_pd(fine%la), get_pmask(fine%la))

    call build(ghost, la, nc, ng = 1)

    call fillpatch(ghost, crse, 1, ir, bc_crse, bc_fine, 1, icomp, bcomp, nc, &
                   no_final_physbc_input = .true.)
    !
    ! Copy fillpatch()d ghost cells to fine.
    ! We want to copy the valid region of ghost -> valid + ghost region of fine.
    ! Got to do it in two stages since copy()s only go from valid -> valid.
    !
    do i = 1, nboxes(fine)
       call push_back(bl, grow(get_ibox(fine,i),ng))
    end do

    call build(ba, bl, sort = .false.)
    call destroy(bl)
    call build(tmpla, ba, get_pd(fine%la), get_pmask(fine%la), &
               explicit_mapping = get_proc(fine%la))
    call destroy(ba)
    call build(tmpfine, tmpla, nc = nc, ng = 0)

    do i = 1, nboxes(fine)
       if ( remote(fine, i) ) cycle
       bx  =  get_ibox(tmpfine,i)
       src => dataptr(fine,    i, bx, icomp, nc)
       dst => dataptr(tmpfine, i, bx, 1    , nc)
       dst =  src
    end do

    call copy(tmpfine, 1, ghost, 1, nc)  ! parallel copy

    do i = 1, nboxes(fine)
       if ( remote(fine, i) ) cycle
       bx  =  get_ibox(tmpfine,i)
       dst => dataptr(fine,    i, bx, icomp, nc)
       src => dataptr(tmpfine, i, bx, 1    , nc)
       dst =  src
    end do
    !
    ! Finish up.
    !
    call fill_boundary(fine, icomp, nc, ng)

    dx = ONE

    call multifab_physbc(fine, icomp, bcomp, nc, dx, bc_fine)

    call destroy(ghost)
    call destroy(tmpfine)

    call destroy(la)
    call destroy(tmpla)

    call destroy(bpt)

  end subroutine

end module multifab_fill_ghost_module
