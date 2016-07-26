module multifab_fill_ghost_module

  use fillpatch_module

  implicit none

contains

  subroutine multifab_fill_ghost_cells(fine,crse,ng,ir,bc_crse,bc_fine,icomp,bcomp,nc, &
                                       stencil_width_input,fourth_order_input, &
                                       fill_crse_boundary_input,fill_crse_physbc_input)

    ! This subroutine used to have an argument called fill_crse_input.  It is renamed
    ! to fill_crse_boundary_input, and its position has changed.  This is intended to
    ! break codes that call this subroutine with "fil_crse_input=.false..
    ! The reason is as follows.  This subroutine is usually called after ml_cc_restriction
    ! is called.  In an old version of ml_cc_restriction, multifab_fill_boundary wass
    ! called on the coarse multifab.  In that case, we did not need to call fill_boundary
    ! here.  However, the new version of ml_cc_restriction no longer calls fill_boundary
    ! for the coarse multifab.  Therefore, we have to call fill_boundary here.  So we
    ! like to break codes that rely on the old behaivor of ml_cc_restriction.

    use layout_module
    use bl_constants_module
    use multifab_physbc_module

    type(multifab), intent(inout) :: fine
    type(multifab), intent(inout) :: crse
    integer       , intent(in   ) :: ng
    integer       , intent(in   ) :: ir(:)
    type(bc_level), intent(in   ) :: bc_crse,bc_fine
    integer       , intent(in   ) :: icomp,bcomp,nc
    logical       , intent(in   ), optional :: fill_crse_boundary_input
    integer       , intent(in   ), optional :: stencil_width_input
    logical       , intent(in   ), optional :: fourth_order_input
    logical       , intent(in   ), optional :: fill_crse_physbc_input

    integer         :: i, j
    type(multifab)  :: ghost
    type(box)       :: bx
    type(layout)    :: la, fine_la
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

    if ( present(fill_crse_boundary_input ) )    fill_crse = fill_crse_boundary_input
    if ( present(fourth_order_input       ) ) fourth_order = fourth_order_input

    if ( nghost(fine) <  ng          ) &
         call bl_error('multifab_fill_ghost_cells: fine does NOT have enough ghost cells')

    if ( .not. cell_centered_q(fine) ) &
         call bl_error('fillpatch: fine is NOT cell centered')

    if (present(stencil_width_input)) then
       if ( fourth_order ) then
          if ( stencil_width_input < 2) &
            call bl_error('fillpatch: fourth_order but stencil_width < 2')
       end if
       stencil_width = stencil_width_input
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
    call layout_build_ba(la, fgasc%ba, get_pd(fine_la), get_pmask(fine_la), explicit_mapping=fgasc%prc)

    ! Don't need to make any ghost cells.
    call multifab_build(ghost, la, nc, ng = 0)

    ! Don't ask fillpatch to fill any ghost cells.
    call fillpatch(ghost, crse, 0, ir, bc_crse, bc_fine, 1, icomp, bcomp, nc, &
                   no_final_physbc_input=.true., fill_crse_input=fill_crse, &
                   stencil_width_input=stencil_width, fourth_order_input=fourth_order, &
                   fill_crse_physbc_input=fill_crse_physbc_input)

    !
    ! Copy fillpatch()d ghost cells to fine.
    !
    !$OMP PARALLEL DO PRIVATE(i,j,bx,src,dst)
    do i = 1, nfabs(ghost)
       j = fgasc%idx(i)
       bx = box_intersection(get_ibox(ghost,i), get_pbox(fine,j))
       src => dataptr(ghost, i, bx, 1    , nc)
       dst => dataptr(fine , j, bx, icomp, nc)
       call cpy_d(dst,src)
    end do
    !$OMP END PARALLEL DO

    call multifab_destroy(ghost)
    call layout_destroy(la)

    call fill_boundary(fine, icomp, nc, ng)
    call multifab_physbc(fine, icomp, bcomp, nc, bc_fine)

    call destroy(bpt)

  end subroutine multifab_fill_ghost_cells

  !
  ! This version of multifab_fill_ghost_cells takes both a crse_old and a crse_new
  !     and does interpolation in time as well as in space.  The coefficient alpha
  !     is used to define:   fine = interp (alpha*crse_old + (1-alpha)*crse_new )
  !
  subroutine multifab_fill_ghost_cells_t(fine,crse_old,crse_new,alpha,ng, &
                                         ir,bc_crse,bc_fine,icomp,bcomp,nc, &
                                         stencil_width_input,fourth_order_input, &
                                         fill_crse_boundary_input)

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
    logical        , intent(in   ), optional :: fill_crse_boundary_input
    integer        , intent(in   ), optional :: stencil_width_input
    logical        , intent(in   ), optional :: fourth_order_input

    logical :: fill_crse
    type(multifab) :: crse

    if (ng == 0) return

    if (crse_new%la /= crse_old%la) then
       call bl_error("multifab_fill_ghost_cells_t: crse_new and crse_old have different layout")
    end if

    if (crse_new%ng /= crse_old%ng) then
       call bl_error("multifab_fill_ghost_cells_t: crse_new and crse_old have different number of ghost cells")
    end if

    fill_crse = .true.
    if ( present(fill_crse_boundary_input) ) fill_crse = fill_crse_boundary_input

    if (fill_crse) then
       call fill_boundary(crse_new, icomp, nc, ng=nghost(crse_new))
       call fill_boundary(crse_old, icomp, nc, ng=nghost(crse_old))
    end if

    call multifab_physbc(crse_new,icomp,bcomp,nc,bc_crse)
    call multifab_physbc(crse_old,icomp,bcomp,nc,bc_crse)

    call multifab_build(crse, crse_old%la, nc=crse_old%nc, ng=crse_old%ng)
    call saxpy(crse, alpha, crse_old, (ONE-alpha), crse_new, all=.true.)

    call multifab_fill_ghost_cells(fine,crse,ng,ir,bc_crse,bc_fine,icomp,bcomp,nc, &
         stencil_width_input,fourth_order_input,.false., .false.)

    call multifab_destroy(crse)

  end subroutine multifab_fill_ghost_cells_t

end module multifab_fill_ghost_module
