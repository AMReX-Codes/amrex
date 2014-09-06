module bndry_reg_module

  use layout_module
  use multifab_module
  use bl_error_module
  use bl_constants_module

  implicit none

  logical, private :: bndry_reg_thin = .false.

  type :: bndry_reg
     integer :: dim   = 0
     integer :: nc    = 1
     integer ::  nbegin(3),  nend(3)
     integer :: onbegin(3), onend(3)
     integer :: ref_ratio(3)
     type(box) :: crse_domain
     logical :: pmask(3)
     logical :: other = .false.
     type(multifab), pointer ::  bmf(:,:) => Null()
     type(multifab), pointer :: obmf(:,:) => Null()
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

  private :: rr_build_nd, rr_build_cc

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
       deallocate(br%bmf,br%laf)
       if (associated(br%indxmap)) deallocate(br%indxmap)
       if (associated(br%facemap)) deallocate(br%facemap)
       if ( br%other ) then
          deallocate(br%obmf,br%olaf,br%oindxmap,br%ofacemap)
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

    if (any(lnodal)) then
       call bl_assert(.not.lother, "bndry_reg_rr_build(): nodal and other cannot be both true.")
       call rr_build_nd(br, la, rr, pdc, lnc, lw, lnodal)
    else
       call rr_build_cc(br, la, lac, rr, pdc, lnc, lw, lother)
    end if

    call destroy(bpt)
  end subroutine bndry_reg_rr_build

  subroutine rr_build_nd(br, la, rr, pdc, nc, width, nodal)
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
  end subroutine rr_build_nd


  subroutine rr_build_cc(br, la, lac, rr, pdc, nc, width, other)
    use vector_i_module
    type(layout),    intent(inout) :: la, lac
    type(bndry_reg), intent(out  ) :: br
    integer,         intent(in   ) :: rr(:)
    type(box),       intent(in   ) :: pdc
    integer,         intent(in   ) :: nc
    integer,         intent(in   ) :: width
    logical,         intent(in   ) :: other

    integer                        :: i, j, kk, id, dm, nb, cnto, cnt, f, ilocal
    integer                        :: nl, ncell, myproc, nlthin, nlthino
    integer                        :: lo(size(rr)), hi(size(rr)), nfb(3), onfb(3)
    integer, allocatable           :: prcc(:), prf(:), pshift(:,:)
    type(box)                      :: rbox, bx
    type(box), allocatable         :: bxs(:), bxso(:), bxsc(:)
    type(box_intersector), pointer :: bi(:)
    type(boxarray)                 :: baa
    type(layout)                   :: laftmp
    type(vector_i) :: oproc, oface

    myproc = parallel_myproc()
    dm     = get_dim(la)
    nb     = nboxes(la)
    nl     = nlocal(la)

    br%dim = dm
    br%nc  = nc
    br%ref_ratio(1:dm) = rr 
    br%crse_domain     = pdc
    br%pmask(1:dm)     = get_pmask(lac)
    br%other = other

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

    allocate(bxs(2*dm*nb))
    allocate(prf(2*dm*nb))
    allocate(br%bmf(1,0:0), br%laf(1,0:0))
    allocate(br%indxmap(2*nl*dm), br%facemap(2*nl*dm))

    if (other) then
       allocate(bxso(2*dm*nb))
       allocate(pshift(dm,2*dm*nb))
       allocate(br%obmf(1,0:0), br%olaf(1,0:0))
    end if

    nlthin = 0
    nfb = 0
    cnt = 0

    if (other) then
       nlthino = 0
       onfb = 0
       cnto = 0
       pshift = 0
       
       call build(oproc)
       call build(oface)
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
                      lo(i) = hi(i) + 1
                      hi(i) = lo(i)
                   end if
                else
                   lo(id) = max(lo(id)-width, pdc%lo(id))
                   hi(id) = min(hi(id)+width, pdc%hi(id))
                end if
             end do

             call build(bx, lo, hi)

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
                nfb(i) = nfb(i) + 1
                br%indxmap(nlthin) = local_index(la,j)
                br%facemap(nlthin) = (2*f-1)*i  ! possible values: -1,+1,-2,+2,-3:+3
             end if

             bxs(cnt) = bx

             if (other) then
                if (br%pmask(i)) then
                   if (lo(i) .lt. pdc%lo(i)) then
                      pshift(i,cnt) = extent(pdc,i)
                   else if (hi(i) .gt. pdc%hi(i)) then
                      pshift(i,cnt) = -extent(pdc,i)
                   end if
                end if
                
                bxso(cnt) = shift(bx, pshift(:,cnt))
                
                bi => layout_get_box_intersector(lac, bxso(cnt))
                cnto = cnto + size(bi)
                do kk=1, size(bi)
                   if (myproc .eq. get_proc(lac,bi(kk)%i)) then
                      call push_back(oproc, local_index(lac,bi(kk)%i))
                      call push_back(oface, (2*f-1)*i)  ! possible values: -1,+1,-2,+2,-3:+3 
                      nlthino = nlthino+1
                      onfb(i) = onfb(i) + 1
                   end if
                end do
                deallocate(bi)
             end if
          end do
       end do
    end do

    call build(baa, bxs(1:cnt), sort = .false.)
    call build(br%laf(1,0), baa, boxarray_bbox(baa), explicit_mapping = prf(1:cnt))
    call build(br%bmf(1,0), br%laf(1,0), nc = nc, ng = 0)
    call destroy(baa)

    if (other) then
       allocate(br%oindxmap(nlthino))
       allocate(br%ofacemap(nlthino))

       allocate(bxsc(cnto), prcc(cnto))
          
       ilocal = 0
       cnto = 0
       do j = 1, cnt
             
          bi => layout_get_box_intersector(lac, bxso(j))
             
          do kk = 1, size(bi)
             cnto = cnto + 1
             bxsc(cnto) = shift(bi(kk)%bx, -pshift(:,j))
             prcc(cnto) = get_proc(lac,bi(kk)%i)
                
             if (myproc .eq. get_proc(lac,bi(kk)%i)) then
                ilocal = ilocal + 1
                br%oindxmap(ilocal) = at(oproc, ilocal)
                br%ofacemap(ilocal) = at(oface, ilocal)
             end if
          end do
             
          deallocate(bi)
       end do
          
       call build(baa, bxsc, sort = .false.)
       call build(br%olaf(1,0), baa, boxarray_bbox(baa), explicit_mapping = prcc)
       deallocate(bxsc, prcc)
       call destroy(baa)
       call build(br%obmf(1,0), br%olaf(1,0), nc = nc, ng = 0)
          
       call destroy(oproc)
       call destroy(oface)

       deallocate(bxso,pshift)       
    end if
    
    deallocate(bxs,prf)
    if (bndry_reg_thin) call destroy(laftmp)

    br%nbegin(1) = 1
    br%nend(1) = br%nbegin(1) + nfb(1) - 1
    do i=2,dm
       br%nbegin(i) = br%nend(i-1) + 1
       br%nend(i) = br%nbegin(i) + nfb(i) - 1
    end do

    if (other) then
       br%onbegin(1) = 1
       br%onend(1) = br%onbegin(1) + onfb(1) - 1
       do i=2,dm
          br%onbegin(i) = br%onend(i-1) + 1
          br%onend(i) = br%onbegin(i) + onfb(i) - 1
       end do
    end if

  end subroutine rr_build_cc


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

    integer                   :: i, j, f, n(2)
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

       do f = 0, n(2)-1
          do i = 1, n(1)
             call copy(br%bmf(i,f), tmf)
          end do
       end do

       call destroy(tmf)
       call destroy(la)
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

  subroutine flux_reg_build(br, la, lac, rr, pdc, nc)
    type(layout),    intent(inout)           :: la, lac
    type(bndry_reg), intent(out  )           :: br
    integer,         intent(in   )           :: rr(:)
    type(box),       intent(in   )           :: pdc
    integer,         intent(in   ), optional :: nc
    integer :: lnc
    lnc = 1;  if (present(nc)) lnc = nc
    call rr_build_cc(br, la, lac, rr, pdc, lnc, 0, .true.)
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

    call build(obmf2, br%olaf(1,0), nc=ncomp(br%obmf(1,0)), ng=0)
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

    call destroy(obmf2)
  end subroutine reflux

end module bndry_reg_module
