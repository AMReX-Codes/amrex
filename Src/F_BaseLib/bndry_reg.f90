module bndry_reg_module

  use layout_module
  use multifab_module

  implicit none

  type, public :: bndry_reg
     integer :: dim = 0
     integer :: nc  = 1
     type(multifab), pointer :: bmf(:,:) => Null()
     type(layout) :: la
     type(layout), pointer :: laf(:,:) => Null()
  end type bndry_reg

  interface destroy
     module procedure bndry_reg_destroy
  end interface

contains

  subroutine bndry_reg_destroy(br)
    type(bndry_reg), intent(inout) :: br
    integer :: i, f
    if ( br%dim /= 0 ) then
       do i = 1, br%dim
          do f = 0, 1
             call destroy(br%bmf(i,f))
             call destroy(br%laf(i,f))
          end do
       end do
       deallocate(br%bmf)
       deallocate(br%laf)
    end if
    br%dim = 0
  end subroutine bndry_reg_destroy

  subroutine bndry_reg_build(br, la, rr, pd, nc, width, nodal)
    type(layout), intent(inout) :: la
    type(bndry_reg), intent(out) :: br
    integer, intent(in) :: rr(:)
    type(box), intent(in) :: pd
    integer, intent(in), optional :: nc
    integer, intent(in), optional :: width
    logical, intent(in), optional :: nodal(:)
    type(box) :: rbox
    integer :: i, j, id, f
    integer :: dm
    integer :: nb
    integer :: lw
    integer :: lo(size(rr)), hi(size(rr))
    type(box), allocatable :: bxs(:)
    type(box) :: lpd
    type(boxarray) :: baa
    logical :: nd_flag(size(rr))
    integer :: lnc

    lnc = 1; if ( present(nc) ) lnc = nc
    lw = 0; if ( present(width) ) lw = width

    dm = layout_dim(la)
    nb = layout_nboxes(la)

    allocate(bxs(nb))
    allocate(br%bmf(dm,0:1))
    allocate(br%laf(dm,0:1))

    br%dim = dm
    br%la  = la
    br%nc  = lnc

    lpd = box_nodalize(pd, nodal)

    nd_flag = .false.; if ( present(nodal) ) nd_flag = nodal

    if ( dm /= layout_dim(la) .or. &
         dm /= box_dim(pd)) then
       call bl_error("BNDRY_REG_BUILD: DIM inconsistent")
    end if

    do i = 1, dm
       do f = 0, 1
          do j = 1, nb

             rbox = coarsen(box_nodalize(get_box(la,j), nodal), rr)

             lo = lwb(rbox)
             hi = upb(rbox)

             select case (f)
             case ( 0 )
                !! LO SIDE OBJECTS
                if ( nd_flag(i) ) then
                   lo(i) = lo(i)
                else
                   lo(i) = lo(i)-1
                end if
                hi(i) = lo(i)
             case ( 1 )
                !! Build hi-side objects
                if ( nd_flag(i) ) then
                   hi(i) = hi(i)
                else
                   hi(i) = hi(i)+1
                end if
                lo(i) = hi(i)
             case default
                call bl_error("BUILD_BNDRY_REG: This can't be happening")
             end select


             ! Grow in the other directions for interping bc's
             ! Grow by lw if possible; if not lw then lw-1, etc; then none.
             ! Note that this makes sure not to leave the physical boundary,
             ! but doesn't see the other grids.  NEEDS TO BE FIXED.

             do id = 1, dm
                if ( id /= i ) then
                   lo(id) = max(lo(id)-lw, lpd%lo(id))
                   hi(id) = min(hi(id)+lw, lpd%hi(id))
                end if
             end do

             call build(bxs(j), lo, hi)

          end do

          call build(baa, bxs, sort = .false.)
          call build(br%laf(i,f), baa)
          call build(br%bmf(i,f), br%laf(i,f), nc = lnc, ng = 0)
          call destroy(baa)

       end do
    end do
  end subroutine bndry_reg_build

  subroutine bndry_reg_copy(br, mf, all)
    type(multifab), intent(in) :: mf
    type(bndry_reg), intent(inout) :: br
    logical, intent(in), optional :: all
    integer :: i, f
    do i = 1, br%dim
       do f = 0, 1
          call copy(br%bmf(i,f), mf, all=all)
       end do
    end do
  end subroutine bndry_reg_copy
  subroutine bndry_reg_copy_c(br, cb, mf, cm, nc, all)
    type(multifab), intent(in) :: mf
    type(bndry_reg), intent(inout) :: br
    integer, intent(in) :: cb, cm
    integer, intent(in), optional :: nc
    logical, intent(in), optional :: all
    integer :: i, f, lnc
    lnc = 1; if ( present(nc) ) lnc = nc
    do i = 1, br%dim
       do f = 0, 1
          call copy(br%bmf(i,f), cb, mf, cm, all=all)
       end do
    end do
  end subroutine bndry_reg_copy_c

  subroutine bndry_reg_setval(br, val, all)
    type(bndry_reg), intent(inout) :: br
    real(kind=dp_t), intent(in) :: val
    logical, intent(in), optional :: all
    integer :: i, f
    do i = 1, br%dim
       do f = 0, 1
          call setval(br%bmf(i,f), val, all=all)
       end do
    end do
  end subroutine bndry_reg_setval

  function bndry_reg_get_boxarray(br, i, f) result(r)
    type(boxarray) :: r
    type(bndry_reg), intent(in) :: br
    integer, intent(in) :: i, f
    r = get_boxarray(br%bmf(i,f))
  end function bndry_reg_get_boxarray

  function bndry_reg_get_layout(br, i, f) result(r)
    type(layout) :: r
    type(bndry_reg), intent(in) :: br
    integer, intent(in) :: i, f
    r = get_layout(br%bmf(i,f))
  end function bndry_reg_get_layout

  subroutine bndry_reg_print(br, str, unit, all, data, skip)
    use bl_IO_module
    type(bndry_reg), intent(in) :: br
    character (len=*), intent(in), optional :: str
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: all, data
    integer, intent(in), optional :: skip
    integer :: i, f, un
    un = unit_stdout(unit)
    call unit_skip(un, skip)
    write(unit=un, fmt='("BNDRY_REG")', advance = 'no')
    if ( present(str) ) then
       write(unit=un, fmt='(": ",A)') str
    else
       write(unit=un, fmt='()')
    end if
    do i = 1, br%dim
       do f = 0, 1
          call multifab_print(br%bmf(i,f), &
               unit = unit, &
               all = all, &
               data = data,  &
               skip = unit_get_skip(skip) + 2)
       end do
    end do
  end subroutine bndry_reg_print

end module bndry_reg_module
