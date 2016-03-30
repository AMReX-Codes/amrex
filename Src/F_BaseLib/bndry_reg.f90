module bndry_reg_module

  use layout_module
  use multifab_module
  use bl_error_module
  use bl_constants_module

  implicit none

  type :: bndry_reg
     integer :: dim   = 0
     integer :: nc    = 1
     integer ::  nbegin(3),  nend(3)
     integer :: onbegin(3), onend(3)
     integer :: ref_ratio(3)
     type(box) :: crse_domain
     logical :: pmask(3) = .false.
     logical :: other    = .false.
     logical :: mask     = .false.
     logical :: nodal(3) = .false.
     integer :: width    = 0
     type(multifab), pointer ::  bmf(:,:) => Null()
     type(multifab), pointer :: obmf(:,:) => Null()
     type(lmultifab) :: uncovered
     type(layout),   pointer ::  laf(:,:) => Null()
     type(layout),   pointer :: olaf(:,:) => Null()
     integer, pointer ::  indxmap(:) => Null()
     integer, pointer :: oindxmap(:) => Null()
     integer, pointer ::  facemap(:) => Null()
     integer, pointer :: ofacemap(:) => Null()
  end type bndry_reg

  interface destroy
     module procedure bndry_reg_destroy
  end interface

contains

  subroutine bndry_reg_destroy(br)
    type(bndry_reg), intent(inout) :: br
    integer :: i, f, n(2)
    if ( br%dim /= 0 ) then
       if (built_q(br%uncovered)) call destroy(br%uncovered)
       n = shape(br%bmf)
       do f = 0, n(2)-1
          do i = 1, n(1)
             call multifab_destroy(br%bmf(i,f))
             call layout_destroy(br%laf(i,f))
             if ( br%other .and. associated(br%obmf) ) then
                call multifab_destroy(br%obmf(i,f))
                call layout_destroy(br%olaf(i,f))
             end if
          end do
       end do
       deallocate(br%bmf,br%laf)
       if (associated(br%indxmap)) deallocate(br%indxmap)
       if (associated(br%facemap)) deallocate(br%facemap)
       if ( br%other ) then
          if (associated(br%obmf)) deallocate(br%obmf)
          if (associated(br%olaf)) deallocate(br%olaf)
          if (associated(br%oindxmap)) deallocate(br%oindxmap)
          if (associated(br%ofacemap)) deallocate(br%ofacemap)
       end if
    end if
    br%dim = 0
  end subroutine bndry_reg_destroy


  subroutine bndry_reg_rr_build_nd(br, la, rr, pdc, nodal, nc)
    type(layout),    intent(inout) :: la
    type(bndry_reg), intent(out  ) :: br
    integer,         intent(in   ) :: rr(:)
    type(box),       intent(in   ) :: pdc
    integer, intent(in), optional  :: nc
    logical, intent(in), optional  :: nodal(:)

    integer                        :: i, j, f, dm, nb, lnc
    integer                        :: lo(la%lap%dim), hi(la%lap%dim)
    logical                        :: lnodal(la%lap%dim)
    type(box), allocatable         :: bxs(:)
    type(box)                      :: rbox, lpdc
    type(boxarray)                 :: baa

    type(bl_prof_timer), save :: bpt

    call build(bpt, "bndry_reg_rr_build_nd")

    lnc    = 1     ;  if ( present(   nc) ) lnc    =    nc
    lnodal = .true.;  if ( present(nodal) ) lnodal = nodal
    
    dm       = get_dim(la)
    nb       = nboxes(la)

    br%dim         = dm
    br%nc          = lnc
    br%other       = .false.
    br%mask        = .false.
    br%nodal(1:dm) = lnodal

    lpdc = box_nodalize(pdc, lnodal)

    allocate(bxs(nb))
    allocate(br%bmf(dm,0:1), br%laf(dm,0:1))

    if ( dm /= get_dim(la) .or. dm /= box_dim(pdc) ) call bl_error("BNDRY_REG_BUILD: DIM inconsistent")

    do i = 1, dm
       do f = 0, 1
          !$omp parallel do private(j,rbox,lo,hi)
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
             call build(bxs(j), lo, hi)
          end do
          !$omp end parallel do

          call boxarray_build_v(baa, bxs, sort = .false.)
          call layout_build_ba(br%laf(i,f), baa, boxarray_bbox(baa), explicit_mapping = get_proc(la))
          call multifab_build(br%bmf(i,f), br%laf(i,f), nc = br%nc, ng = 0)
          call boxarray_destroy(baa)
       end do
    end do

    call destroy(bpt)
  end subroutine bndry_reg_rr_build_nd


  subroutine bndry_reg_rr_build(br, la, lac, rr, pdc, nc, width, other)
    use vector_i_module
    type(layout),    intent(inout)           :: la, lac
    type(bndry_reg), intent(out  )           :: br
    integer,         intent(in   )           :: rr(:)
    type(box),       intent(in   )           :: pdc
    integer,         intent(in   ), optional :: nc
    integer,         intent(in   ), optional :: width
    logical,         intent(in   ), optional :: other

    integer :: lnc, lwidth
    logical :: lother, mask

    integer :: i, j, kk, id, dm, nb, f, nl, myproc, iproc, iface, nlbtot, onlbtot, pshift
    integer :: lo(size(rr)), hi(size(rr)), nlb(3), onlb(3)
    type(box) :: rbox, bx, obx, bxtmp
    type(box), allocatable :: bxsc(:)
    type(box_intersector), pointer :: bi(:)
    type(boxarray) :: baa
    type(layout) :: laftmp
    type(vector_i) :: indx, face, oindx, oface, procf, procc
    type(list_box) :: bxl, obxl
    integer :: nbxs, revshft(size(rr))
    integer, allocatable :: shft(:,:)
    type(box), allocatable :: bxs(:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "bndry_reg_rr_build")

    lnc     = 1 ;       if ( present(nc)    ) lnc    = nc
    lwidth  = 0 ;       if ( present(width) ) lwidth = width
    lother  = .false.;  if ( present(other) ) lother = other

    mask = .not. lother

    if (lother .and. lwidth>0) call bl_error("bndry_reg_rr_build: lother .and. lwidth>0")

    myproc = parallel_myproc()
    dm     = get_dim(la)
    nb     = nboxes(la)
    nl     = nlocal(la)

    br%dim = dm
    br%nc  = lnc
    br%ref_ratio(1:dm) = rr 
    br%crse_domain     = pdc
    br%pmask(1:dm)     = get_pmask(lac)
    br%other = lother
    br%mask  = mask
    br%nodal = .false.
    br%width = lwidth

    nlb = 0
    call reserve(indx, 2*dm*nb)
    call reserve(face, 2*dm*nb)

    if (br%other) then
       onlb = 0
       call reserve(oindx, 2*dm*nb)
       call reserve(oface, 2*dm*nb)
    end if
       
    do i = 1, dm
       do f = 0, 1
          do j = 1, nb

             rbox = coarsen(get_box(la,j), rr)

             lo = lwb(rbox)
             hi = upb(rbox)
             do id = 1, dm
                if ( id .eq. i ) then
                   if (f.eq.0) then
                      lo(i) = lo(i) - 1
                      hi(i) = lo(i)
                   else
                      lo(id) = hi(id) + 1
                      hi(id) = lo(id)
                   end if
                else
                   lo(id) = lo(id)-lwidth
                   hi(id) = hi(id)+lwidth
                end if
             end do

             call build(bx, lo, hi)

             iproc = get_proc(la,j)
             iface = (2*f-1)*i  ! possible values: -1,+1,-2,+2,-3:+3

             call push_back(procf, iproc)

             if (myproc .eq. iproc) then
                nlb(i) = nlb(i) + 1
                call push_back(indx, local_index(la,j))
                call push_back(face,iface)
             end if
             
             call push_back(bxl, bx)

             if (br%other) then

                pshift = 0
                if (br%pmask(i)) then
                   if (bx%lo(i) .lt. pdc%lo(i)) then
                      pshift = extent(pdc,i)
                   else if (bx%hi(i) .gt. pdc%hi(i)) then
                      pshift = -extent(pdc,i)
                   end if
                end if
                
                obx = shift(bx, pshift, i)
                
                bi => layout_get_box_intersector(lac, obx)
                
                do kk=1, size(bi)
                   iproc = get_proc(lac,bi(kk)%i)
                   call push_back(procc, iproc)
                   if (myproc .eq. iproc) then
                      onlb(i) = onlb(i) + 1
                      call push_back(oindx, local_index(lac,bi(kk)%i))
                      call push_back(oface, iface)
                   end if
                   call push_back(obxl, shift(bi(kk)%bx, -pshift, i))
                end do

                deallocate(bi)

             end if

          end do
       end do
    end do

    allocate(br%laf(1,0:0))
    allocate(br%bmf(1,0:0))

    call boxarray_build_l(baa, bxl, sort = .false.)
    call destroy(bxl)
    call layout_build_ba(br%laf(1,0), baa, boxarray_bbox(baa), explicit_mapping = dataptr(procf))
    call boxarray_destroy(baa)
    call multifab_build(br%bmf(1,0), br%laf(1,0), nc = br%nc, ng = 0)

    nlbtot = sum(nlb(1:dm))
    if (nlbtot .ne. size(indx)) call bl_error("bndry_reg_rr_build: how did this happen")

    allocate(br%indxmap(nlbtot))
    allocate(br%facemap(nlbtot))

    br%indxmap(1:nlbtot) = dataptr(indx)
    br%facemap(1:nlbtot) = dataptr(face)

    call destroy(procf)
    call destroy(indx)
    call destroy(face)

    if (br%other .and. .not.empty(obxl)) then

       allocate(br%olaf(1,0:0))
       allocate(br%obmf(1,0:0))

       call boxarray_build_l(baa, obxl, sort = .false.)
       call destroy(obxl)
       call layout_build_ba(br%olaf(1,0), baa, boxarray_bbox(baa), explicit_mapping = dataptr(procc))
       call boxarray_destroy(baa)
       call multifab_build(br%obmf(1,0), br%olaf(1,0), nc = br%nc, ng = 0)
       
       onlbtot = sum(onlb(1:dm))
       if (onlbtot .ne. size(oindx)) call bl_error("bndry_reg_rr_build other: how did this happen")

       allocate(br%oindxmap(onlbtot))
       allocate(br%ofacemap(onlbtot))

       br%oindxmap(1:onlbtot) = dataptr(oindx)
       br%ofacemap(1:onlbtot) = dataptr(oface)

       call destroy(procc)
       call destroy(oindx)
       call destroy(oface)

    end if
    
    br%nbegin(1) = 1
    br%nend(1) = br%nbegin(1) + nlb(1) - 1
    do i=2,dm
       br%nbegin(i) = br%nend(i-1) + 1
       br%nend(i) = br%nbegin(i) + nlb(i) - 1
    end do

    if (br%other .and. associated(br%obmf)) then
       br%onbegin(1) = 1
       br%onend(1) = br%onbegin(1) + onlb(1) - 1
       do i=2,dm
          br%onbegin(i) = br%onend(i-1) + 1
          br%onend(i) = br%onbegin(i) + onlb(i) - 1
       end do
    end if

    if (mask) then
       call lmultifab_build(br%uncovered, br%laf(1,0), nc=1, ng=0)
       call lmultifab_setval(br%uncovered, .true.)

       ! Build a coarsen version of the fine boxarray
       allocate(bxsc(nb))
       do i = 1, nb
          bxsc(i) = coarsen(get_box(la,i), rr)
       end do
       call boxarray_build_v(baa, bxsc, sort = .false.)
       deallocate(bxsc)
       call layout_build_ba(laftmp, baa, boxarray_bbox(baa), explicit_mapping = get_proc(la))
       call boxarray_destroy(baa)

       if (any(br%pmask(1:dm))) then
          nbxs = 3**dm
          allocate(shft(nbxs,dm))
          allocate(bxs(nbxs))
       end if

       do i = 1, nfabs(br%uncovered)

          bx = get_box(br%uncovered,i)
          bi => layout_get_box_intersector(laftmp, bx)
          do kk=1, size(bi)
             ! Cells under fine are marked with false.
             call lfab_setval_bx(br%uncovered%fbs(i), .false., bi(kk)%bx)
          end do
          deallocate(bi)

          ! Cells outside physical domain are marked with false.
          bxl = boxlist_box_diff(bx,pdc)
          do while (.not. empty(bxl))
             bxtmp = front(bxl)
             call pop_front(bxl)
             call lfab_setval_bx(br%uncovered%fbs(i), .false., bxtmp)
          end do
          call destroy(bxl)

          if (any(br%pmask(1:dm))) then
             call box_periodic_shift(pdc, bx, br%nodal(1:dm), br%pmask(1:dm), 0, shft, nbxs, bxs)
             do j = 1, nbxs
                revshft = -shft(j,:)
                bxtmp = shift(bxs(j), revshft) ! this is the original unshifted sub-box
                ! these cells are marked true because they are in the 'periodic' domain
                call lfab_setval_bx(br%uncovered%fbs(i), .true., bxtmp)
                !unless they are under fine
                bi => layout_get_box_intersector(laftmp, bxs(j))
                do kk = 1, size(bi)
                   call lfab_setval_bx(br%uncovered%fbs(i), .false., shift(bi(kk)%bx,revshft))
                end do
                deallocate(bi)
             end do
          end if

       end do

       call layout_destroy(laftmp)
    end if

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

    dm             = get_dim(la)
    nb             = nboxes(la)
    br%dim         = dm
    br%nc          = lnc
    br%other       = .false.
    br%mask        = .false.
    br%nodal(1:dm) = lnodal

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

          call boxarray_build_v(baa, bxs, sort = .false.)
          call layout_build_ba(br%laf(i,f), baa, boxarray_bbox(baa), explicit_mapping = get_proc(la))
          call multifab_build(br%bmf(i,f), br%laf(i,f), nc = lnc, ng = 0)
          call boxarray_destroy(baa)

       end do
    end do
    call destroy(bpt)
  end subroutine bndry_reg_build

  subroutine bndry_reg_copy(br, mf)
    type(multifab) , intent(inout) :: mf
    type(bndry_reg), intent(inout) :: br

    integer                   :: i, j, f, n(2)
    type(multifab)            :: mftmp
    type(layout)              :: mf_la
    type(box)                 :: domain
    logical                   :: have_periodic_boxes
    type(bl_prof_timer), save :: bpt

    call build(bpt, "br_copy")

    n = shape(br%bmf)

    mf_la = get_layout(mf)

    have_periodic_boxes = .false.

    if ( any(get_pmask(mf_la)) ) then

       domain = grow(get_pd(mf_la), nghost(mf), .not. get_pmask(mf_la))

       loop: do f = 0, n(2)-1
          do i = 1, n(1)
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
       ! Need to fill the ghost cells of the crse array before copying from them.
       !
       if (nghost(mf) .ge. br%width) then
          call multifab_fill_boundary(mf)
          do f = 0, n(2)-1
             do i = 1, n(1)
                call copy(br%bmf(i,f), mf, ngsrc=br%width)
             end do
          end do
       else
          call multifab_build(mftmp, mf_la, ncomp(mf), br%width, mf%nodal)
          call copy(mftmp, mf)
          call multifab_fill_boundary(mftmp)
          do f = 0, n(2)-1
             do i = 1, n(1)
                call copy(br%bmf(i,f), mftmp, ngsrc=br%width)
             end do
          end do
          call destroy(mftmp)
       end if

    else
       do f = 0, n(2)-1
          do i = 1, n(1)
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
    if (.not.associated(br%obmf)) return
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

  ! subroutine bndry_reg_copy_c(br, cb, mf, cm, nc)
  !   type(multifab), intent(in) :: mf
  !   type(bndry_reg), intent(inout) :: br
  !   integer, intent(in) :: cb, cm
  !   integer, intent(in), optional :: nc
  !   integer :: i, f, n(2)
  !   type(bl_prof_timer), save :: bpt
  !   call build(bpt, "br_copy_c")
  !   n = shape(br%bmf)
  !   do f = 0, n(2)-1
  !      do i = 1, n(1)
  !         call copy(br%bmf(i,f), cb, mf, cm, nc = nc)
  !      end do
  !   end do
  !   call destroy(bpt)
  ! end subroutine bndry_reg_copy_c

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

  subroutine flux_reg_build(br, la, lac, rr, pdc, nc)
    type(layout),    intent(inout)           :: la, lac
    type(bndry_reg), intent(out  )           :: br
    integer,         intent(in   )           :: rr(:)
    type(box),       intent(in   )           :: pdc
    integer,         intent(in   ), optional :: nc
    call bndry_reg_rr_build(br, la, lac, rr, pdc, nc, other=.true.)
  end subroutine flux_reg_build

  subroutine flux_reg_crse_init(br, flux, s)
    type(bndry_reg), intent(inout) :: br
    type(multifab) , intent(in   ) :: flux(:)
    real(kind=dp_t), intent(in   ) :: s 
    
    integer :: ibr, iflx, idim, face 
    integer, dimension(3) :: blo, bhi
    type(box) :: br_box
    real(kind=dp_t), pointer :: bp(:,:,:,:), fp(:,:,:,:)

    call bndry_reg_setval(br, ZERO)  ! zero bmf

    if (.not.associated(br%obmf)) return

    !$omp parallel private(ibr,iflx,idim,face,blo,bhi,br_box,bp,fp)
    blo = 1; bhi = 1
    !$omp do
    do ibr = 1, nfabs(br%obmf(1,0))
       iflx = br%oindxmap(ibr)
       idim = abs (   br%ofacemap(ibr))
       face = sign(1, br%ofacemap(ibr))

       bp     => dataptr(br%obmf(1,0),ibr)
       br_box = get_ibox(br%obmf(1,0),ibr)

       if (br%pmask(idim) .and. (.not. contains(br%crse_domain,br_box)) ) then
          if ( face .eq. -1 ) then
             br_box = shift(br_box,  extent(br%crse_domain,idim), idim)
          else
             br_box = shift(br_box, -extent(br%crse_domain,idim), idim)
          end if
       end if

       if (face .eq. -1) then
          br_box = shift(br_box, 1, idim)
       end if

       blo(1:br%dim) = lwb(br_box)
       bhi(1:br%dim) = upb(br_box)
       
       fp => dataptr(flux(idim),iflx)
       
       bp = s * fp(blo(1):bhi(1),blo(2):bhi(2),blo(3):bhi(3),:)
    end do
    !$omp end do
    !$omp end parallel
  end subroutine flux_reg_crse_init

  subroutine flux_reg_fine_add(br, flux, s)
    type(bndry_reg), intent(inout) :: br
    type(multifab) , intent(in   ) :: flux(:)
    real(kind=dp_t), intent(in   ) :: s 

    integer :: ibr, iflx, idim, face, i, j, k, ii, jj, kk, n, ir1, ir2, ir3
    integer, dimension(br%dim) :: blo, bhi, flo, fhi
    type(box) :: br_box, flx_box
    real(kind=dp_t), pointer :: bp(:,:,:,:), fp(:,:,:,:)

    !$omp parallel do private(ibr,iflx,idim,face,i,j,k,ii,jj,kk,n,ir1,ir2,ir3) &
    !$omp private(blo,bhi,flo,fhi,br_box,flx_box,bp,fp)
    do ibr = 1, nfabs(br%bmf(1,0))
       iflx = br%indxmap(ibr)
       idim = abs (   br%facemap(ibr))
       face = sign(1, br%facemap(ibr))

       bp     => dataptr(br%bmf(1,0),ibr)
       br_box = get_ibox(br%bmf(1,0),ibr)
       blo    = lwb(br_box)
       bhi    = upb(br_box)
          
       fp      => dataptr(flux(idim),iflx)
       flx_box = get_ibox(flux(idim),iflx)
       flo     = lwb(flx_box)
       fhi     = upb(flx_box)
       
       ! boxes are globally zero based in active dimensions
       
       select case(br%dim)
       case (1) ! 1D
          
          if (face .eq. -1) then
             bp(blo(1),1,1,:) = bp(blo(1),1,1,:) + s*fp(flo(1),1,1,:)
          else
             bp(bhi(1),1,1,:) = bp(bhi(1),1,1,:) + s*fp(fhi(1),1,1,:)
          end if
          
       case (2) ! 2D
          
          select case (idim)
          case (1) ! 2D-x

             if (face .eq. -1) then
                i  = blo(1)
                ii = flo(1)
             else
                i  = bhi(1)
                ii = fhi(1) 
             end if
             ir1 = 0

             do n=1,size(fp,4)
                do j=blo(2),bhi(2)
                   jj = j * br%ref_ratio(2)
                   do ir2 = 0, br%ref_ratio(2)-1
                      bp(i,j,1,n) = bp(i,j,1,n) + s*fp(ii+ir1,jj+ir2,1,n)
                   end do
                end do
             end do

          case (2) ! 2D-y

             if (face .eq. -1) then
                j  = blo(2)
                jj = flo(2)
             else
                j  = bhi(2)
                jj = fhi(2) 
             end if
             ir2 = 0
             
             do n=1,size(fp,4)
                do i=blo(1),bhi(1)
                   ii = i * br%ref_ratio(1)
                   do ir1 = 0, br%ref_ratio(1)-1
                      bp(i,j,1,n) = bp(i,j,1,n) + s*fp(ii+ir1,jj+ir2,1,n)
                   end do
                end do
             end do
             
          end select

       case (3) ! 3D

          select case (idim)
          case (1) ! 3D-x

             if (face .eq. -1) then
                i  = blo(1)
                ii = flo(1)
             else
                i  = bhi(1)
                ii = fhi(1) 
             end if
             ir1 = 0

             do n=1,size(fp,4)
                do k=blo(3),bhi(3)
                   kk = k * br%ref_ratio(3)
                   do j=blo(2),bhi(2)
                      jj = j * br%ref_ratio(2)
                      do ir3 = 0, br%ref_ratio(3)-1
                         do ir2 = 0, br%ref_ratio(2)-1
                            bp(i,j,k,n) = bp(i,j,k,n) + s*fp(ii+ir1,jj+ir2,kk+ir3,n)
                         end do
                      end do
                   end do
                end do
             end do

          case (2) ! 3D-y

             if (face .eq. -1) then
                j  = blo(2)
                jj = flo(2)
             else
                j  = bhi(2)
                jj = fhi(2) 
             end if
             ir2 = 0

             do n=1,size(fp,4)
                do k=blo(3),bhi(3)
                   kk = k * br%ref_ratio(3)
                   do i=blo(1),bhi(1)
                      ii = i * br%ref_ratio(1)
                      do ir3 = 0, br%ref_ratio(3)-1
                         do ir1 = 0, br%ref_ratio(1)-1
                            bp(i,j,k,n) = bp(i,j,k,n) + s*fp(ii+ir1,jj+ir2,kk+ir3,n)
                         end do
                      end do
                   end do
                end do
             end do

          case (3)

             if (face .eq. -1) then
                k  = blo(3)
                kk = flo(3)
             else
                k  = bhi(3)
                kk = fhi(3) 
             end if
             ir3 = 0

             do n=1,size(fp,4)
                do j=blo(2),bhi(2)
                   jj = j * br%ref_ratio(2)
                   do i=blo(1),bhi(1)
                      ii = i * br%ref_ratio(1)
                      do ir2 = 0, br%ref_ratio(2)-1
                         do ir1 = 0, br%ref_ratio(1)-1
                            bp(i,j,k,n) = bp(i,j,k,n) + s*fp(ii+ir1,jj+ir2,kk+ir3,n)
                         end do
                      end do
                   end do
                end do
             end do

          end select
       end select
    end do
    !$omp end parallel do
  end subroutine flux_reg_fine_add
    
  subroutine reflux(u, br, s)
    type(multifab) , intent(inout) :: u
    type(bndry_reg), intent(in   ) :: br
    real(kind=dp_t), intent(in   ) :: s 

    type(multifab) :: obmf2
    integer :: blo(3), bhi(3)
    integer :: idim, ibr, iu
    real(kind=dp_t) :: face
    type(box) :: br_box
    real(kind=dp_t), dimension(:,:,:,:), pointer :: b1p, b2p, up

    if (.not.associated(br%obmf)) return

    call multifab_build(obmf2, br%olaf(1,0), nc=ncomp(br%obmf(1,0)), ng=0)
    call copy (obmf2, br%bmf(1,0), bndry_reg_to_other=.true.)

    !$omp parallel private(blo,bhi,idim,ibr,iu,face,br_box,b1p,b2p,up)
    blo = 1; bhi = 1
    do idim = 1, br%dim
       !$omp do
       do ibr = br%onbegin(idim), br%onend(idim)
          iu = br%oindxmap(ibr)
          face = sign(1, br%ofacemap(ibr))

          b1p    => dataptr(br%obmf(1,0), ibr)
          b2p    => dataptr(   obmf2    , ibr)
          br_box = get_ibox(   obmf2    , ibr)

          if (br%pmask(idim) .and. (.not. contains(br%crse_domain,br_box)) ) then
             if ( face .eq. -1 ) then
                br_box = shift(br_box,  extent(br%crse_domain,idim), idim)
             else
                br_box = shift(br_box, -extent(br%crse_domain,idim), idim)
             end if
          end if

          blo(1:br%dim) = lwb(br_box)
          bhi(1:br%dim) = upb(br_box)

          up => dataptr(u, iu)

          up     (blo(1):bhi(1),blo(2):bhi(2),blo(3):bhi(3),:) = &
               up(blo(1):bhi(1),blo(2):bhi(2),blo(3):bhi(3),:)   &
               - face * s * (b1p + b2p)
       end do
       !$omp end do
    end do
    !$omp end parallel

    call multifab_destroy(obmf2)
  end subroutine reflux

end module bndry_reg_module
