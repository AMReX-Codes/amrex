module bndry_reg_module

  use layout_module
  use multifab_module

  implicit none

  type, public :: bndry_reg
     integer :: dim = 0
     integer :: nc  = 1
     type(multifab), pointer ::  bmf(:,:) => Null()
     type(multifab), pointer :: obmf(:,:) => Null()
     type(layout) :: la
     type(layout), pointer :: laf(:,:) => Null()
     type(layout), pointer :: olaf(:,:) => Null()
  end type bndry_reg

  interface build
     module procedure bndry_reg_rr_build
     module procedure bndry_reg_build
  end interface

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
             call destroy(br%obmf(i,f))
             call destroy(br%laf(i,f))
             call destroy(br%olaf(i,f))
          end do
       end do
       deallocate(br%bmf)
       deallocate(br%obmf)
       deallocate(br%laf)
       deallocate(br%olaf)
    end if
    br%dim = 0
  end subroutine bndry_reg_destroy

  subroutine bndry_reg_rr_build_1(br, la, lac, rr, pd, nc, width, nodal)
    type(layout), intent(inout) :: la, lac
    type(bndry_reg), intent(out) :: br
    integer, intent(in) :: rr(:)
    type(box), intent(in) :: pd
    integer, intent(in), optional :: nc
    integer, intent(in), optional :: width
    logical, intent(in), optional :: nodal(:)
    type(box) :: rbox
    integer :: i, j, id, f, ff
    integer :: dm
    integer :: nb
    integer :: lw
    integer :: lo(size(rr)), hi(size(rr))
    integer :: lo1(size(rr)), hi1(size(rr))
    type(box), allocatable :: bxs(:), bxs1(:), bxsc(:)
    integer, allocatable :: prcc(:)
    type(box) :: lpd
    type(box) :: bx, bx1
    type(boxarray) :: baa
    logical :: nd_flag(size(rr))
    integer :: lnc, cnt, k

    lnc = 1; if ( present(nc) ) lnc = nc
    lw = 0; if ( present(width) ) lw = width

    dm = get_dim(la)
    nb = nboxes(la)

    allocate(bxs(nb), bxs1(nb))
    allocate(br%bmf(dm,0:1))
    allocate(br%obmf(dm,0:1))
    allocate(br%laf(dm,0:1))
    allocate(br%olaf(dm,0:1))

    br%dim = dm
    br%la  = la
    br%nc  = lnc

    lpd = box_nodalize(pd, nodal)

    nd_flag = .false.; if ( present(nodal) ) nd_flag = nodal

    if ( dm /= get_dim(la) .or. &
         dm /= box_dim(pd)) then
       call bl_error("BNDRY_REG_BUILD: DIM inconsistent")
    end if

    do i = 1, dm
       do f = 0, 1
          ff = -1; if ( f == 1 ) ff = 1
          cnt = 0
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

             lo1 = lo
             hi1 = hi

             do id = 1, dm
                if ( id /= i ) then
                   lo1(id) = max(lo1(id), lpd%lo(id))
                   hi1(id) = min(hi1(id), lpd%hi(id))
                end if
             end do
             call build(bxs1(j), lo1, hi1)
             do k = 1, nboxes(lac)
! print *, i, f, j, k
! call print(bxs1(j), "bxs1")
! call print(get_box(lac,k), "gblk")
! call print(grow(get_box(lac,k),1,i,ff), "grow")
                bx = intersection(bxs1(j), grow(get_box(lac, k),1,i,ff))
! call print(bx, "bx")
                if ( .not. empty(bx) ) then
                   cnt = cnt + 1
                end if
             end do

             do id = 1, dm
                if ( id /= i ) then
                   lo(id) = max(lo(id)-lw, lpd%lo(id))
                   hi(id) = min(hi(id)+lw, lpd%hi(id))
                end if
             end do

             call build(bxs(j), lo, hi)
! call print(bxs(j), "bxs")

          end do

! print *, 'CNT = ', cnt
          call build(baa, bxs, sort = .false.)
          call build(br%laf(i,f), baa, explicit_mapping = get_proc(la))
          call build(br%bmf(i,f), br%laf(i,f), nc = lnc, ng = 0)
          call destroy(baa)

          allocate(bxsc(cnt))
          allocate(prcc(cnt))
          cnt = 1
          do j = 1, nb
             do k = 1, nboxes(lac)
                bx = intersection(bxs1(j), grow(get_box(lac,k),1,i,ff))
                if ( .not. empty(bx) ) then
                   lo = lwb(bx)
                   hi = upb(bx)
                   do id = 1, dm
                      if ( id /= i ) then
                         lo(id) = max(lo(id)-lw, lpd%lo(id))
                         hi(id) = min(hi(id)+lw, lpd%hi(id))
                      end if
                   end do
                   call build(bx1, lo, hi)
                   bxsc(cnt) = bx1
                   prcc(cnt) = get_proc(lac,k)
                   cnt = cnt + 1
                end if
             end do
          end do
! print *, 'PRCC = ', prcc
          call build(baa, bxsc, sort = .false.)
          call build(br%olaf(i,f), baa, explicit_mapping = prcc)
          deallocate(bxsc, prcc)
          call destroy(baa)
          call build(br%obmf(i,f), br%olaf(i,f), nc = lnc, ng = 0)
       end do
    end do
! if ( parallel_ioprocessor() ) then
! do i = 1, dm
!    do f = 0, 1
!       call print(br%laf(i,f), unit = 50)
!       call print(br%olaf(i,f), unit = 51)
! end do
! end do
! end if
! call parallel_barrier()
  end subroutine bndry_reg_rr_build_1

  subroutine bndry_reg_rr_build(br, la, rr, pd, nc, width, nodal)
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

    dm = get_dim(la)
    nb = nboxes(la)

    allocate(bxs(nb))
    allocate(br%bmf(dm,0:1))
    allocate(br%obmf(dm,0:1))
    allocate(br%laf(dm,0:1))
    allocate(br%olaf(dm,0:1))

    br%dim = dm
    br%la  = la
    br%nc  = lnc

    lpd = box_nodalize(pd, nodal)

    nd_flag = .false.; if ( present(nodal) ) nd_flag = nodal

    if ( dm /= get_dim(la) .or. &
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
          call build(br%olaf(i,f), baa) ! FIXME
          call build(br%bmf(i,f), br%laf(i,f), nc = lnc, ng = 0)
          call build(br%obmf(i,f), br%olaf(i,f), nc = lnc, ng = 0)
          call destroy(baa)

       end do
    end do
  end subroutine bndry_reg_rr_build

  subroutine bndry_reg_build(br, la, pd, nc, nodal)
    type(layout), intent(inout) :: la
    type(bndry_reg), intent(out) :: br
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

    dm = get_dim(la)
    nb = nboxes(la)

    allocate(lo(dm),hi(dm),nd_flag(dm))

    allocate(bxs(nb))
    allocate(br%bmf(dm,0:1))
    allocate(br%obmf(dm,0:1))
    allocate(br%laf(dm,0:1))
    allocate(br%olaf(dm,0:1))

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
          call build(br%olaf(i,f), baa) ! FIXME
          call build(br%bmf(i,f), br%laf(i,f), nc = lnc, ng = 0)
          call build(br%obmf(i,f), br%olaf(i,f), nc = lnc, ng = 0)
          call destroy(baa)

       end do
    end do
  end subroutine bndry_reg_build

  subroutine bndry_reg_copy(br, mf, all)
    type(multifab), intent(in) :: mf
    type(bndry_reg), intent(inout) :: br
    logical, intent(in), optional :: all
    integer :: i, f
    type(bl_prof_timer), save :: bpt

    call build(bpt, "br_copy")

    do i = 1, br%dim
       do f = 0, 1
          call copy(br%bmf(i,f), mf, all=all)
       end do
    end do

    call destroy(bpt)
  end subroutine bndry_reg_copy

  subroutine bndry_reg_copy_to_other(br)
    type(bndry_reg), intent(inout) :: br
    integer :: i, f
    do i = 1, br%dim
       do f = 0, 1
          call copy(br%obmf(i,f), br%bmf(i,f))
       end do
    end do
  end subroutine bndry_reg_copy_to_other

  subroutine bndry_reg_copy_from_other(br)
    type(bndry_reg), intent(inout) :: br
    integer :: i, f
    do i = 1, br%dim
       do f = 0, 1
          call copy(br%bmf(i,f), br%obmf(i,f))
       end do
    end do
  end subroutine bndry_reg_copy_from_other

  subroutine bndry_reg_copy_c(br, cb, mf, cm, nc, all)
    type(multifab), intent(in) :: mf
    type(bndry_reg), intent(inout) :: br
    integer, intent(in) :: cb, cm
    integer, intent(in), optional :: nc
    logical, intent(in), optional :: all
    integer :: i, f
    type(bl_prof_timer), save :: bpt

    call build(bpt, "br_copy_c")

    do i = 1, br%dim
       do f = 0, 1
          call copy(br%bmf(i,f), cb, mf, cm, nc = nc, all=all)
       end do
    end do

    call destroy(bpt)

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
          call print(br%bmf(i,f), &
               unit = unit, &
               all = all, &
               data = data,  &
               skip = unit_get_skip(skip) + 2)
       end do
    end do
  end subroutine bndry_reg_print

end module bndry_reg_module
