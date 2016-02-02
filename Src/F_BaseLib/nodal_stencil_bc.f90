module nodal_stencil_bc_module

  use bl_constants_module
  use bl_types
  use multifab_module
  use bc_functions_module

  implicit none

contains

  subroutine stencil_set_bc_nodal(sdim, bx, nbx, idx, mask, face_type, pd_periodic, la_periodic)
    integer,         intent(in   ) :: sdim
    type(box),       intent(in   ) :: bx, nbx
    type(imultifab), intent(inout) :: mask
    integer,         intent(in   ) :: idx
    integer,         intent(in   ) :: face_type(:,:,:)
    type(box),       intent(in   ) :: pd_periodic
    type(layout),    intent(inout) :: la_periodic

    integer, pointer :: mp(:,:,:,:)
    type(box)        :: bx1
    type(boxarray)   :: ba
    integer          :: ii, dm, ib, jb, kb, jb_lo, kb_lo
    logical          :: nodal(sdim)

    nodal = .true.
    !
    ! Set the mask to BC_DIR or BC_NEU based on face_type at a physical boundary.
    !
    do dm = 1, sdim
       !
       ! Lo side
       !
       bx1 = nbx
       bx1%hi(dm) = bx1%lo(dm)
       mp => dataptr(mask, idx, bx1)
       if (face_type(idx,dm,1) == BC_NEU) then
          mp = ibset(mp, BC_BIT(BC_NEU, dm, -1))
       else if (face_type(idx,dm,1) == BC_DIR) then
          mp = ibset(mp, BC_BIT(BC_DIR, 1, 0))
       end if
       !
       ! Hi side
       !
       bx1 = nbx
       bx1%lo(dm) = bx1%hi(dm)
       mp => dataptr(mask, idx, bx1)
       if (face_type(idx,dm,2) == BC_NEU) then
          mp = ibset(mp, BC_BIT(BC_NEU, dm, +1))
       else if (face_type(idx,dm,2) == BC_DIR) then
          mp = ibset(mp, BC_BIT(BC_DIR, 1, 0))
       end if
    end do
    !
    ! Set the mask to BC_DIR at coarse-fine boundaries.
    !
    jb_lo = -1; if (sdim .lt. 2) jb_lo = 1
    kb_lo = -1; if (sdim .lt. 3) kb_lo = 1

    do kb = kb_lo, 1
       do jb = jb_lo, 1
          do ib = -1, 1
             bx1 = shift(bx,ib,1)
             if (sdim > 1) bx1 = shift(bx1,jb,2)
             if (sdim > 2) bx1 = shift(bx1,kb,3)
             bx1 = intersection(bx1, pd_periodic)
             if ( empty(bx1) ) cycle
             call layout_boxarray_diff(ba, bx1, la_periodic)
             do ii = 1, nboxes(ba)
                bx1 = intersection(box_nodalize(get_box(ba,ii),nodal), nbx)
                if ( empty(bx1) ) cycle
                mp => dataptr(mask, idx, bx1)
                mp = ibset(mp, BC_BIT(BC_DIR,1,0))
             end do
             call boxarray_destroy(ba)
          end do
       end do
    end do

  end subroutine stencil_set_bc_nodal

end module nodal_stencil_bc_module
