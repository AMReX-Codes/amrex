module bndry_reg_module

  use layout_module
  use multifab_module
  use bl_error_module

  implicit none

  logical, private :: bndry_reg_thin = .false.

  type :: bndry_reg
     integer :: dim   = 0
     integer :: nc    = 1
     logical :: other = .false.
     type(multifab), pointer ::  bmf(:,:) => Null()
     type(multifab), pointer :: obmf(:,:) => Null()
     type(layout),   pointer ::  laf(:,:) => Null()
     type(layout),   pointer :: olaf(:,:) => Null()
     integer, pointer ::  indxmap(:,:) => Null()
     integer, pointer :: oindxmap(:,:) => Null()
     integer, pointer ::  facemap(:,:) => Null()
     integer, pointer :: ofacemap(:,:) => Null()
  end type bndry_reg

  interface destroy
     module procedure bndry_reg_destroy
  end interface

  private :: rr_build, rr_build_other

contains

  subroutine bndry_reg_set_thin(v)
    logical, intent(in) :: v
    bndry_reg_thin = v
  end subroutine bndry_reg_set_thin
  pure function bndry_reg_get_thin() result(r)
    logical :: r
    r = bndry_reg_thin
  end function bndry_reg_get_thin

  subroutine bndry_reg_destroy(br)
    type(bndry_reg), intent(inout) :: br
    integer :: i, f, n(2)
    if ( br%dim /= 0 ) then
       n = shape(br%bmf)
       do f = 0, n(2)-1
          do i = 1, n(1)
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
          deallocate(br%obmf, br%olaf)
          deallocate(br%indxmap, br%facemap)
          deallocate(br%oindxmap, br%ofacemap)
       end if
    end if
    br%dim = 0
  end subroutine bndry_reg_destroy

  subroutine bndry_reg_rr_build(br, la, lac, rr, pdc, nc, width, nodal, other)
    type(layout),    intent(inout)           :: la, lac
    type(bndry_reg), intent(out  )           :: br
    integer,         intent(in   )           :: rr(:)
    type(box),       intent(in   )           :: pdc
    integer,         intent(in   ), optional :: nc
    integer,         intent(in   ), optional :: width
    logical,         intent(in   ), optional :: nodal(:)
    logical,         intent(in   ), optional :: other

    integer :: lnc, lw
    logical :: lnodal(la%lap%dim), lother
    
    type(bl_prof_timer), save :: bpt

    call build(bpt, "bndry_reg_rr_build")

    lnc    = 1 ;       if ( present(nc)    ) lnc    = nc
    lw     = 0 ;       if ( present(width) ) lw     = width
    lnodal = .false. ; if ( present(nodal) ) lnodal = nodal
    lother = .true.  ; if ( present(other) ) lother = other

    if (lother) then
       call bl_assert(.not.lnodal, "bndry_reg_rr_build(): nodal and other cannot be both true.")
       call bl_assert(lw.eq.0, "bndry_reg_rr_build(): width must be zero when other is true.")
       call rr_build_other(br, la, lac, rr, pdc, lnc)
    else
       call rr_build(br, la, rr, pdc, lnc, lw, lnodal)
    end if

    call destroy(bpt)
  end subroutine bndry_reg_rr_build

  subroutine rr_build(br, la, rr, pdc, nc, width, nodal)
    type(layout),    intent(inout) :: la
    type(bndry_reg), intent(out  ) :: br
    integer,         intent(in   ) :: rr(:)
    type(box),       intent(in   ) :: pdc
    integer,         intent(in   ) :: nc
    integer,         intent(in   ) :: width
    logical,         intent(in   ) :: nodal(:)

    integer                        :: i, j, id, f, dm, nb
    integer                        :: lo(la%lap%dim), hi(la%lap%dim)
    type(box), allocatable         :: bxs(:)
    type(box)                      :: rbox, lpdc
    type(boxarray)                 :: baa

    dm       = get_dim(la)
    nb       = nboxes(la)

    br%dim   = dm
    br%nc    = nc
    br%other = .false.

    lpdc     = box_nodalize(pdc, nodal)

    allocate(bxs(nb))
    allocate(br%bmf(dm,0:1), br%laf(dm,0:1))

    if ( dm /= get_dim(la) .or. dm /= box_dim(pdc) ) call bl_error("BNDRY_REG_BUILD: DIM inconsistent")

    do i = 1, dm
       do f = 0, 1
          do j = 1, nb
             rbox = coarsen(box_nodalize(get_box(la,j),nodal), rr)
             lo   = lwb(rbox)
             hi   = upb(rbox)
             if ( f == 0 ) then
                if ( .not. nodal(i) ) lo(i) = lo(i) - 1
                hi(i) = lo(i)
             else
                if ( .not. nodal(i) ) hi(i) = hi(i) + 1
                lo(i) = hi(i)
             end if

             do id = 1, dm
                if ( id /= i ) then
                   lo(id) = max(lo(id)-width, lpdc%lo(id))
                   hi(id) = min(hi(id)+width, lpdc%hi(id))
                end if
             end do

             call build(bxs(j), lo, hi)
          end do

          call build(baa, bxs, sort = .false.)
          call build(br%laf(i,f), baa, boxarray_bbox(baa), explicit_mapping = get_proc(la))
          call build(br%bmf(i,f), br%laf(i,f), nc = nc, ng = 0)
          call destroy(baa)
       end do
    end do
  end subroutine rr_build

  subroutine rr_build_other(br, la, lac, rr, pdc, nc)
    use vector_i_module
    type(layout),    intent(inout) :: la, lac
    type(bndry_reg), intent(out  ) :: br
    integer,         intent(in   ) :: rr(:)
    type(box),       intent(in   ) :: pdc
    integer,         intent(in   ) :: nc

    logical                        :: pmask(size(rr))
    integer                        :: i, j, kk, dm, nb, cnto, cnt, f, ilocal, nlmax
    integer                        :: nl, ncell, myproc, nlthin, nlthino(size(rr))
    integer                        :: lo(size(rr)), hi(size(rr))
    integer                        :: lof(size(rr),0:1), hif(size(rr),0:1)
    integer, allocatable           :: prcc(:), prf(:), pshift(:)
    type(box)                      :: rbox, bx
    type(box), allocatable         :: bxs(:), bxso(:), bxsc(:)
    type(box_intersector), pointer :: bi(:)
    type(boxarray)                 :: baa
    type(layout)                   :: laftmp
    integer, pointer :: tindxmap(:), tfacemap(:)
    integer, allocatable, target :: oindxmap1(:), oindxmap2(:), oindxmap3(:)
    integer, allocatable, target :: ofacemap1(:), ofacemap2(:), ofacemap3(:)
    type(vector_i) :: oproc, oface

    myproc = parallel_myproc()
    dm     = get_dim(la)
    nb     = nboxes(la)
    nl     = nlocal(la)
    pmask  = get_pmask(lac)

    br%dim   = dm
    br%nc    = nc
    br%other = .true.

    if (bndry_reg_thin) then
       ! Build a coarsen version of the fine boxarray
       allocate(bxsc(nb))
       do i = 1, nb
          bxsc(i) = coarsen(get_box(la,i), rr)
       end do
       call build(baa, bxsc, sort = .false.)
       deallocate(bxsc)
       call build(laftmp, baa, boxarray_bbox(baa), explicit_mapping = get_proc(la))
       call destroy(baa)
    end if

    allocate(bxs(2*nb))
    allocate(bxso(2*nb))
    allocate(pshift(2*nb))
    allocate(prf(2*nb))

    allocate(br%bmf(dm,0:0), br%laf(dm,0:0))
    allocate(br%obmf(dm,0:0), br%olaf(dm,0:0))
    allocate(br%indxmap(2*nl,dm), br%facemap(2*nl,dm))

    do i = 1, dm
       nlthin = 0
       nlthino(i) = 0
       cnt = 0
       cnto = 0
       pshift = 0

       call build(oproc)
       call build(oface)

       do j = 1, nb
          rbox = coarsen(get_box(la,j), rr)
          lo   = lwb(rbox)
          hi   = upb(rbox)
          ! lo face
          lof(:,0) = lo  
          hif(:,0) = hi
          lof(i,0) = lof(i,0) - 1
          hif(i,0) = lof(i,0)
          ! hi face
          lof(:,1) = lo  
          hif(:,1) = hi
          lof(i,1) = hif(i,1) + 1
          hif(i,1) = lof(i,1)

          do f = 0, 1
             call build(bx, lof(:,f), hif(:,f))

             if (bndry_reg_thin) then
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
             if (myproc .eq. prf(cnt)) then
                nlthin = nlthin + 1
                br%indxmap(nlthin,i) = local_index(la,j)
                br%facemap(nlthin,i) = (2*f-1)*i  ! possible values: -1,+1,-2,+2,-3:+3
             end if

             bxs(cnt) = bx

             if (pmask(i)) then
                if (lof(i,f) .lt. pdc%lo(i)) then
                   pshift(cnt) = extent(pdc,i)
                else if (hif(i,f) .gt. pdc%hi(i)) then
                   pshift(cnt) = -extent(pdc,i)
                end if
             end if

             bxso(cnt) = shift(bx, pshift(cnt), i)

             bi => layout_get_box_intersector(lac, bxso(cnt))
             cnto = cnto + size(bi)
             do kk=1, size(bi)
                if (myproc .eq. get_proc(lac,bi(kk)%i)) then
                   call push_back(oproc, local_index(lac,bi(kk)%i))
                   call push_back(oface, (2*f-1)*i)  ! possible values: -1,+1,-2,+2,-3:+3 
                   nlthino(i) = nlthino(i)+1
                end if
             end do
             deallocate(bi)
          end do
       end do

       call build(baa, bxs(1:cnt), sort = .false.)
       call build(br%laf(i,0), baa, boxarray_bbox(baa), explicit_mapping = prf(1:cnt))
       call build(br%bmf(i,0), br%laf(i,0), nc = nc, ng = 0)
       call destroy(baa)

       allocate(bxsc(cnto), prcc(cnto))
       if (i .eq. 1) then
          allocate(oindxmap1(nlthino(i)), ofacemap1(nlthino(i)))
          tindxmap => oindxmap1
          tfacemap => ofacemap1
       else if (i .eq. 2) then
          allocate(oindxmap2(nlthino(i)), ofacemap2(nlthino(i)))
          tindxmap => oindxmap2
          tfacemap => ofacemap2
       else
          allocate(oindxmap3(nlthino(i)), ofacemap3(nlthino(i)))
          tindxmap => oindxmap3
          tfacemap => ofacemap3
       end if

       ilocal = 0
       cnto = 0
       do j = 1, cnt

          bi => layout_get_box_intersector(lac, bxso(j))

          do kk = 1, size(bi)
             cnto = cnto + 1
             bxsc(cnto) = shift(bi(kk)%bx, -pshift(j), i)
             prcc(cnto) = get_proc(lac,bi(kk)%i)

             if (myproc .eq. get_proc(lac,bi(kk)%i)) then
                ilocal = ilocal + 1
                tindxmap(ilocal) = at(oproc, ilocal)
                tfacemap(ilocal) = at(oface, ilocal)
             end if
          end do

          deallocate(bi)
       end do

       call build(baa, bxsc, sort = .false.)
       call build(br%olaf(i,0), baa, boxarray_bbox(baa), explicit_mapping = prcc)
       deallocate(bxsc, prcc)
       call destroy(baa)
       call build(br%obmf(i,0), br%olaf(i,0), nc = nc, ng = 0)

       nullify(tindxmap,tfacemap)
       call destroy(oproc)
       call destroy(oface)
    end do

    nlmax = maxval(nlthino)
    allocate(br%oindxmap(nlmax, dm))
    allocate(br%ofacemap(nlmax, dm))

    br%oindxmap(1:nlthino(1),1) = oindxmap1
    br%ofacemap(1:nlthino(1),1) = ofacemap1
    deallocate(oindxmap1,ofacemap1)
    if (br%dim .gt. 1) then
       br%oindxmap(1:nlthino(2),2) = oindxmap2
       br%ofacemap(1:nlthino(2),2) = ofacemap2
       deallocate(oindxmap2,ofacemap2)
    end if
    if (br%dim .gt. 2) then
       br%oindxmap(1:nlthino(3),3) = oindxmap3
       br%ofacemap(1:nlthino(3),3) = ofacemap3
       deallocate(oindxmap3,ofacemap3)
    end if

    deallocate(bxs,bxso,pshift,prf)
    if (bndry_reg_thin) call destroy(laftmp)
  end subroutine rr_build_other


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

  subroutine bndry_reg_copy(br, mf, filled)
    type(multifab) , intent(inout) :: mf
    type(bndry_reg), intent(inout) :: br
    logical, intent(in), optional  :: filled

    integer                   :: i, j, f
    type(list_box)            :: bl
    type(multifab)            :: tmf
    type(boxarray)            :: ba
    type(layout)              :: la,mf_la
    type(box)                 :: domain
    logical                   :: have_periodic_boxes
    real(kind=dp_t), pointer  :: src(:,:,:,:), dst(:,:,:,:)
    type(bl_prof_timer), save :: bpt
    logical                   :: lfilled

    call build(bpt, "br_copy")

    lfilled = .false.;  if (present(filled)) lfilled = filled

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
       if (.not.lfilled) call multifab_fill_boundary(mf)

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
    integer :: i, f, n(2)
    type(bl_prof_timer), save :: bpt
    if ( .not. br%other ) call bl_error('bndry_reg_copy_to_other: other not defined')
    call build(bpt, "br_copy_to_other")
    n = shape(br%bmf)
    do f = 0, n(2)-1
       do i = 1, n(1)
          call copy(br%obmf(i,f), br%bmf(i,f), bndry_reg_to_other=.true.)
       end do
    end do
    call destroy(bpt)
  end subroutine bndry_reg_copy_to_other

  ! subroutine bndry_reg_copy_from_other(br)
  !   type(bndry_reg), intent(inout) :: br
  !   integer :: i, f
  !   type(bl_prof_timer), save :: bpt
  !   if ( .not. br%other ) call bl_error('bndry_reg_copy_from_other: other not defined')
  !   call build(bpt, "br_copy_from_other")
  !   do i = 1, br%dim
  !      do f = 0, 1
  !         call copy(br%bmf(i,f), br%obmf(i,f))
  !      end do
  !   end do
  !   call destroy(bpt)
  ! end subroutine bndry_reg_copy_from_other

  subroutine bndry_reg_copy_c(br, cb, mf, cm, nc)
    type(multifab), intent(in) :: mf
    type(bndry_reg), intent(inout) :: br
    integer, intent(in) :: cb, cm
    integer, intent(in), optional :: nc
    integer :: i, f, n(2)
    type(bl_prof_timer), save :: bpt
    call build(bpt, "br_copy_c")
    n = shape(br%bmf)
    do f = 0, n(2)-1
       do i = 1, n(1)
          call copy(br%bmf(i,f), cb, mf, cm, nc = nc)
       end do
    end do
    call destroy(bpt)
  end subroutine bndry_reg_copy_c

  subroutine bndry_reg_setval(br, val, all)
    type(bndry_reg), intent(inout) :: br
    real(kind=dp_t), intent(in) :: val
    logical, intent(in), optional :: all
    integer :: i, f, n(2)
    n = shape(br%bmf)
    do f = 0, n(2)-1
       do i = 1, n(1)
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
    integer :: i, f, n(2), un
    un = unit_stdout(unit)
    call unit_skip(un, skip)
    write(unit=un, fmt='("BNDRY_REG")', advance = 'no')
    if ( present(str) ) then
       write(unit=un, fmt='(": ",A)') str
    else
       write(unit=un, fmt='()')
    end if
    n = shape(br%bmf)
    do f = 0, n(2)-1
       do i = 1, n(1)
          call print(br%bmf(i,f), &
               unit = unit, &
               all = all, &
               data = data,  &
               skip = unit_get_skip(skip) + 2)
       end do
    end do
  end subroutine bndry_reg_print

end module bndry_reg_module
