module flux_reg_module

  use layout_module
  use multifab_module

  implicit none

  type, public :: flux_reg
     integer :: dim = 0
     integer :: nc  = 1
     type(multifab), pointer :: bmf(:,:) => Null()
     type(layout) :: la
     type(layout), pointer :: laf(:,:) => Null()
  end type flux_reg

  interface destroy
     module procedure flux_reg_destroy
  end interface

contains

  subroutine flux_reg_destroy(br)
    type(flux_reg), intent(inout) :: br
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
  end subroutine flux_reg_destroy

  subroutine flux_reg_build(br, la, pd, nc, nodal)
    type(layout), intent(inout) :: la
    type(flux_reg), intent(out) :: br
    type(box), intent(in) :: pd
    integer, intent(in), optional :: nc
    logical, intent(in), optional :: nodal(:)
    type(box) :: rbox
    integer :: i, j, id, f
    integer :: dm,nb
    type(box), allocatable :: bxs(:)
    type(box) :: lpd
    type(boxarray) :: baa
    integer, allocatable :: lo(:),hi(:)
    logical, allocatable :: nd_flag(:)
    integer :: lnc

    lnc = 1; if ( present(nc) ) lnc = nc

    dm = layout_dim(la)
    nb = layout_nboxes(la)

    allocate(lo(dm),hi(dm),nd_flag(dm))

    allocate(bxs(nb))
    allocate(br%bmf(dm,0:1))
    allocate(br%laf(dm,0:1))

    br%dim = dm
    br%la  = la
    br%nc  = lnc

    lpd = box_nodalize(pd, nodal)

    nd_flag = .false.; if ( present(nodal) ) nd_flag = nodal

    if (dm /= box_dim(pd)) call bl_error("BNDRY_REG_BUILD: DIM inconsistent")

    do i = 1, dm
       do f = 0, 1
          do j = 1, nb

             rbox = get_box(la,j)

             lo = lwb(rbox)
             hi = upb(rbox)

             select case (f)
             case ( 0 )
                !! Build lo-side objects
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
                call bl_error("BUILD_FLUX_REG: This can't be happening")
             end select

             call build(bxs(j), lo, hi)

          end do

          call build(baa, bxs, sort = .false.)
          call build(br%laf(i,f), baa)
          call build(br%bmf(i,f), br%laf(i,f), nc = lnc, ng = 0)
          call destroy(baa)

       end do
    end do

    deallocate(lo,hi,nd_flag)

  end subroutine flux_reg_build

  subroutine flux_reg_copy(br, mf, all)
    type(multifab), intent(in) :: mf
    type(flux_reg), intent(inout) :: br
    logical, intent(in), optional :: all
    integer :: i, f
    do i = 1, br%dim
       do f = 0, 1
          call copy(br%bmf(i,f), mf, all=all)
       end do
    end do
  end subroutine flux_reg_copy
  subroutine flux_reg_copy_c(br, cb, mf, cm, nc, all)
    type(multifab), intent(in) :: mf
    type(flux_reg), intent(inout) :: br
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
  end subroutine flux_reg_copy_c

  subroutine flux_reg_setval(br, val, all)
    type(flux_reg), intent(inout) :: br
    real(kind=dp_t), intent(in) :: val
    logical, intent(in), optional :: all
    integer :: i, f
    do i = 1,br%dim
       do f = 0, 1
          call setval(br%bmf(i,f), val, all=all)
       end do
    end do
  end subroutine flux_reg_setval

  function flux_reg_get_boxarray(br, i, f) result(r)
    type(boxarray) :: r
    type(flux_reg), intent(in) :: br
    integer, intent(in) :: i, f
    r = get_boxarray(br%bmf(i,f))
  end function flux_reg_get_boxarray

  function flux_reg_get_layout(br, i, f) result(r)
    type(layout) :: r
    type(flux_reg), intent(in) :: br
    integer, intent(in) :: i, f
    r = get_layout(br%bmf(i,f))
  end function flux_reg_get_layout

  subroutine flux_reg_print(br, str, unit, all, data, skip)
    use bl_IO_module
    type(flux_reg), intent(in) :: br
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
  end subroutine flux_reg_print

end module flux_reg_module
