module bndry_reg_module

  use layout_module
  use multifab_module

  implicit none

  type :: bndry_reg
     integer :: dim   = 0
     integer :: nc    = 1
     logical :: other = .false.
     type(multifab), pointer ::  bmf(:,:) => Null()
     type(multifab), pointer :: obmf(:,:) => Null()
     type(layout),   pointer ::  laf(:,:) => Null()
     type(layout),   pointer :: olaf(:,:) => Null()
     integer, pointer :: indexmap(:,:,:) => Null()
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
             if ( br%other ) then
                call destroy(br%obmf(i,f))
                call destroy(br%olaf(i,f))
             end if
          end do
       end do
       deallocate(br%bmf)
       deallocate(br%laf)
       if ( br%other ) then
          deallocate(br%obmf)
          deallocate(br%olaf)
          deallocate(br%indexmap)
       end if
    end if
    br%dim = 0
  end subroutine bndry_reg_destroy

  subroutine bndry_reg_rr_build(br, la, lac, rr, pdc, nc, width, nodal, other)
    use bl_error_module
    type(layout),    intent(inout)           :: la, lac
    type(bndry_reg), intent(out  )           :: br
    integer,         intent(in   )           :: rr(:)
    type(box),       intent(in   )           :: pdc
    integer,         intent(in   ), optional :: nc
    integer,         intent(in   ), optional :: width
    logical,         intent(in   ), optional :: nodal(:)
    logical,         intent(in   ), optional :: other

    integer                        :: i, j, id, f, kk, dm, nb, lw, lnc, cnto, cnt
    integer                        :: nl, ncell, myproc, nlthin
    integer                        :: lo(size(rr)), hi(size(rr))
    integer                        :: lo1(size(rr)), hi1(size(rr))
    type(box), allocatable         :: bxs(:), bxs1(:), bxsc(:), bxc(:)
    type(box_intersector), pointer :: bi(:)
    integer, allocatable           :: prcc(:), prf(:)
    type(box)                      :: rbox, lpdc, bx, bx1
    type(boxarray)                 :: baa, bac
    type(layout)                   :: lactmp, laftmp
    logical                        :: lnodal(size(rr)), lother
    logical                        :: pmask(size(rr))
    integer, allocatable           :: shifted(:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "bndry_reg_rr_build")

    lnc    = 1 ;       if ( present(nc)    ) lnc    = nc
    lw     = 0 ;       if ( present(width) ) lw     = width
    lnodal = .false. ; if ( present(nodal) ) lnodal = nodal
    lother = .true.  ; if ( present(other) ) lother = other

    call bl_assert(.not.(lnodal.and.lother), "bndry_reg_rr_build(): nodal and other cannot be both true.")

    myproc   = parallel_myproc()
    dm       = get_dim(la)
    nb       = nboxes(la)
    nl       = nlocal(la)
    br%dim   = dm
    br%nc    = lnc
    br%other = lother
    lpdc     = box_nodalize(pdc, lnodal)

    allocate(bxs(nb))
    allocate(br%bmf(dm,0:1), br%laf(dm,0:1))

    if (br%other) then
       allocate(bxs1(nb))
       allocate(br%obmf(dm,0:1), br%olaf(dm,0:1))
       allocate(shifted(nb))
       allocate(br%indexmap(nl,dm,0:1))
    end if

    if ( dm /= get_dim(la) .or. dm /= box_dim(pdc) ) call bl_error("BNDRY_REG_BUILD: DIM inconsistent")

    if ( br%other ) then
       !
       ! Build a layout to be used in intersection tests below.
       !
       allocate(bxc(nboxes(lac)))

       !$OMP PARALLEL DO PRIVATE(i,j,bx,lo,hi)
       do i = 1, nboxes(lac)
          bx = box_nodalize(get_box(lac,i),lnodal)
          lo = lwb(bx)
          hi = upb(bx)
          do j = 1, dm
             if ( .not. lnodal(j) ) then
                if ( lo(j) == lwb(lpdc,j) ) lo(j) = lo(j) - 1
                if ( hi(j) == upb(lpdc,j) ) hi(j) = hi(j) + 1
             end if
          end do
          call build(bxc(i), lo, hi)
       end do
       !$OMP END PARALLEL DO

       call build(bac, bxc, sort = .false.)
       deallocate(bxc)
       call build(lactmp, bac, boxarray_bbox(bac), explicit_mapping = get_proc(lac))
       call destroy(bac)

       pmask = get_pmask(lac)

       ! Build a coarsen version of the fine boxarray
       allocate(bxc(nb))
       do i = 1, nb
          bxc(i) = coarsen(get_box(la,i), rr)
       end do
       call build(bac, bxc, sort = .false.)
       deallocate(bxc)
       call build(laftmp, bac, boxarray_bbox(bac), explicit_mapping = get_proc(la))
       call destroy(bac)
    end if

    allocate(prf(nb))

    do i = 1, dm
       do f = 0, 1
          nlthin = 0
          cnt = 0
          cnto = 0
          if ( br%other ) shifted = 0

          do j = 1, nb
             rbox = coarsen(box_nodalize(get_box(la,j),lnodal), rr)
             lo   = lwb(rbox)
             hi   = upb(rbox)
             if ( f == 0 ) then
                if ( .not. lnodal(i) ) lo(i) = lo(i) - 1
                hi(i) = lo(i)
             else
                if ( .not. lnodal(i) ) hi(i) = hi(i) + 1
                lo(i) = hi(i)
             end if

             do id = 1, dm
                if ( id /= i ) then
                   lo(id) = max(lo(id)-lw, lpdc%lo(id))
                   hi(id) = min(hi(id)+lw, lpdc%hi(id))
                end if
             end do

             call build(bx, lo, hi)

             if (br%other) then
                bi => layout_get_box_intersector(laftmp, bx)
                ncell = 0
                do kk=1, size(bi)
                   ncell = ncell + volume(bi(kk)%bx)
                end do
                deallocate(bi)
                if (ncell .eq. volume(bx)) then
                   cycle  ! bx is entirely covered by fine grids
                end if
             end if

             cnt = cnt+1
             prf(cnt) = get_proc(la,j)
             if (br%other) then
                if (myproc .eq. prf(cnt)) then
                   nlthin = nlthin + 1
                   br%indexmap(nlthin,i,f) = local_index(la,j)
                end if
             end if
             call build(bxs(cnt), lo, hi)

             !
             ! Grow in the other directions for interping bc's
             ! Grow by lw if possible; if not lw then lw-1, etc; then none.
             ! Note that this makes sure not to leave the physical boundary,
             ! but doesn't see the other grids.
             !
             if ( br%other ) then
                lo1 = lo
                hi1 = hi
                if (pmask(i)) then
                   if (lo(i) .lt. pdc%lo(i)) then
                      shifted(cnt) = extent(pdc,i)
                      lo1(i) = lo1(i) + shifted(cnt)
                      hi1(i) = lo1(i)
                   else if (hi(i) .gt. pdc%hi(i)) then
                      shifted(cnt) = -extent(pdc,i)
                      lo1(i) = lo1(i) + shifted(cnt)
                      hi1(i) = lo1(i)
                   end if
                end if

                call build(bxs1(cnt), lo1, hi1)
                bi => layout_get_box_intersector(lactmp, bxs1(cnt))
                cnto = cnto + size(bi)
                deallocate(bi)
             end if
          end do

          call build(baa, bxs(1:cnt), sort = .false.)
          call build(br%laf(i,f), baa, boxarray_bbox(baa), explicit_mapping = prf(1:cnt))
          call build(br%bmf(i,f), br%laf(i,f), nc = lnc, ng = 0)
          call destroy(baa)

          if ( br%other ) then
             allocate(bxsc(cnto), prcc(cnto))
             cnto = 1
             do j = 1, cnt

                bi => layout_get_box_intersector(lactmp, bxs1(j))

                do kk = 1, size(bi)
                   lo = lwb(bi(kk)%bx)
                   hi = upb(bi(kk)%bx)
                   do id = 1, dm
                      if ( id /= i ) then
                         lo(id) = max(lo(id)-lw, lpdc%lo(id))
                         hi(id) = min(hi(id)+lw, lpdc%hi(id))
                      end if
                   end do
                   if (shifted(j) .ne. 0) then
                      lo(i) = lo(i) - shifted(j)
                      hi(i) = lo(i)
                   end if
                   call build(bx1, lo, hi)
                   bxsc(cnto) = bx1
                   prcc(cnto) = get_proc(lac,bi(kk)%i)
                   cnto = cnto + 1
                end do
                deallocate(bi)
             end do
             call build(baa, bxsc, sort = .false.)
             call build(br%olaf(i,f), baa, boxarray_bbox(baa), explicit_mapping = prcc)
             deallocate(bxsc, prcc)
             call destroy(baa)
             call build(br%obmf(i,f), br%olaf(i,f), nc = lnc, ng = 0)
          end if
       end do
    end do

    if ( br%other ) then
       deallocate(bxs1,shifted)
       call destroy(lactmp)
       call destroy(laftmp)
    end if

    deallocate(prf)
    call destroy(bpt)
  end subroutine bndry_reg_rr_build

  subroutine bndry_reg_build(br, la, pd, nc, nodal)
    use bl_error_module
    type(layout),    intent(in   )           :: la
    type(bndry_reg), intent(out  )           :: br
    type(box),       intent(in   )           :: pd
    integer,         intent(in   ), optional :: nc
    logical,         intent(in   ), optional :: nodal(:)

    type(box)              :: rbox
    integer                :: i, j, f, dm, nb, lnc
    type(box), allocatable :: bxs(:)
    type(boxarray)         :: baa
    integer                :: lo(la%lap%dim), hi(la%lap%dim)
    logical                :: lnodal(la%lap%dim)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "bndry_reg_build")

    lnc    = 1;       if ( present(nc)    ) lnc    = nc
    lnodal = .false.; if ( present(nodal) ) lnodal = nodal

    dm       = get_dim(la)
    nb       = nboxes(la)
    br%dim   = dm
    br%nc    = lnc
    br%other = .false.

    allocate(bxs(nb))
    allocate(br%bmf(dm,0:1), br%laf(dm,0:1))

    if (dm /= box_dim(pd)) call bl_error("BNDRY_REG_BUILD: DIM inconsistent")

    do i = 1, dm
       do f = 0, 1

          !$OMP PARALLEL DO PRIVATE(j,rbox,lo,hi)
          do j = 1, nb
             rbox = get_box(la,j)
             lo   = lwb(rbox)
             hi   = upb(rbox)

             if ( f == 0 ) then
                !! Build lo-side objects
                hi(i) = lo(i)
             else
                !! Build hi-side objects
                if ( .not. lnodal(i) ) hi(i) = hi(i)+1
                lo(i) = hi(i)
             end if

             call build(bxs(j), lo, hi)
          end do
          !$OMP END PARALLEL DO

          call build(baa, bxs, sort = .false.)
          call build(br%laf(i,f), baa, boxarray_bbox(baa), explicit_mapping = get_proc(la))
          call build(br%bmf(i,f), br%laf(i,f), nc = lnc, ng = 0)
          call destroy(baa)

       end do
    end do
    call destroy(bpt)
  end subroutine bndry_reg_build

  subroutine bndry_reg_copy(br, mf)
    type(multifab) , intent(inout) :: mf
    type(bndry_reg), intent(inout) :: br

    integer                   :: i, j, f
    type(list_box)            :: bl
    type(multifab)            :: tmf
    type(boxarray)            :: ba
    type(layout)              :: la,mf_la
    type(box)                 :: domain
    logical                   :: have_periodic_boxes
    real(kind=dp_t), pointer  :: src(:,:,:,:), dst(:,:,:,:)
    type(bl_prof_timer), save :: bpt

    call build(bpt, "br_copy")

    mf_la = get_layout(mf)

    have_periodic_boxes = .false.

    if ( any(get_pmask(mf_la)) ) then

       domain = grow(get_pd(mf_la), nghost(mf), .not. get_pmask(mf_la))

       loop: do i = 1, br%dim
          do f = 0, 1
             do j = 1, nboxes(br%bmf(i,f)%la)
                if ( .not. contains(domain, get_box(br%bmf(i,f)%la,j)) ) then
                   have_periodic_boxes = .true.
                   exit loop
                end if
             end do
          end do
       end do loop

    end if

    if ( have_periodic_boxes ) then
       !
       ! We're periodic & have boxes that extend outside the domain in periodic
       ! direction.  In order to fill those boxes we do the usual trick of copy()ing
       ! from a multifab whose valid region has been extended to cover the ghost region.
       !
       do i = 1, nboxes(mf%la)
          call push_back(bl, grow(box_nodalize(get_box(mf%la,i),mf%nodal),nghost(mf)))
       end do
       !
       ! Need to fill the ghost cells of the crse array before copying from them.
       !
       call multifab_fill_boundary(mf)

       call build(ba, bl, sort = .false.)
       call destroy(bl)
       call build(la, ba, get_pd(mf_la), get_pmask(mf_la), explicit_mapping = get_proc(mf_la))
       call destroy(ba)
       call build(tmf, la, nc = ncomp(mf), ng = 0)

       do i = 1, nfabs(mf)
          src => dataptr(mf,  i)
          dst => dataptr(tmf, i)
          call cpy_d(dst,src)
       end do

       do i = 1, br%dim
          do f = 0, 1
             call copy(br%bmf(i,f), tmf)
          end do
       end do

       call destroy(tmf)
       call destroy(la)
    else
       do i = 1, br%dim
          do f = 0, 1
             call copy(br%bmf(i,f), mf)
          end do
       end do
    end if

    call destroy(bpt)
  end subroutine bndry_reg_copy

  subroutine bndry_reg_copy_to_other(br)
    type(bndry_reg), intent(inout) :: br
    integer :: i, f
    type(bl_prof_timer), save :: bpt
    if ( .not. br%other ) call bl_error('bndry_reg_copy_to_other: other not defined')
    call build(bpt, "br_copy_to_other")
    do i = 1, br%dim
       do f = 0, 1
          call copy(br%obmf(i,f), br%bmf(i,f))
       end do
    end do
    call destroy(bpt)
  end subroutine bndry_reg_copy_to_other

  subroutine bndry_reg_copy_from_other(br)
    type(bndry_reg), intent(inout) :: br
    integer :: i, f
    type(bl_prof_timer), save :: bpt
    if ( .not. br%other ) call bl_error('bndry_reg_copy_from_other: other not defined')
    call build(bpt, "br_copy_from_other")
    do i = 1, br%dim
       do f = 0, 1
          call copy(br%bmf(i,f), br%obmf(i,f))
       end do
    end do
    call destroy(bpt)
  end subroutine bndry_reg_copy_from_other

  subroutine bndry_reg_copy_c(br, cb, mf, cm, nc)
    type(multifab), intent(in) :: mf
    type(bndry_reg), intent(inout) :: br
    integer, intent(in) :: cb, cm
    integer, intent(in), optional :: nc
    integer :: i, f
    type(bl_prof_timer), save :: bpt
    call build(bpt, "br_copy_c")
    do i = 1, br%dim
       do f = 0, 1
          call copy(br%bmf(i,f), cb, mf, cm, nc = nc)
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
