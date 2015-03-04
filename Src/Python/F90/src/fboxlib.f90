module fboxlib

  use iso_c_binding
  use blobjects
  implicit none

contains

  subroutine pybl_cstring(cstr, slen, fstr)
    integer(c_int),    intent(in) :: slen
    character(c_char), intent(in) :: cstr(slen)
    character(len=slen), intent(out) :: fstr

    integer :: i

    do i = 1, slen
       fstr(i:i) = cstr(i)
    end do
  end subroutine pybl_cstring

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! open, close, set_comm

  subroutine pybl_hello() bind(c)
    print *, 'Hello from Fortran, with FBoxLib!'
  end subroutine pybl_hello

  subroutine pybl_open() bind(c)
    use boxlib
    call boxlib_initialize()
  end subroutine pybl_open

  subroutine pybl_close() bind(c)
    use boxlib
    call boxlib_finalize()
  end subroutine pybl_close

  subroutine pybl_mpi_rank(r) bind(c)
    use parallel
    integer(c_int), intent(out) :: r
    r = parallel_myproc()
  end subroutine pybl_mpi_rank

  subroutine pybl_mpi_size(r) bind(c)
    use parallel
    integer(c_int), intent(out) :: r
    r = parallel_nprocs()
  end subroutine pybl_mpi_size

  subroutine pybl_mpi_reduce_max(r) bind(c)
    use parallel
    real(c_double), intent(inout) :: r
    real(c_double) :: s
    s = r
    call parallel_reduce(r, s, MPI_MAX)
  end subroutine pybl_mpi_reduce_max

  ! subroutine pybl_set_comm(comm)
  !   use parallel
  !   implicit none
  !   integer, intent(in) :: comm
  !   call parallel_set_m_comm(comm)
  ! end subroutine pybl_set_comm

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! boxarray routines

  subroutine pybl_boxarray_create_from_boxes(boxes,nboxes,dim,cptr) bind(c)
    integer(c_int), intent(in   ), value :: dim, nboxes
    integer(c_int), intent(in   )        :: boxes(dim,2,nboxes)
    type(c_ptr),    intent(  out)        :: cptr

    type(boxarray), pointer :: ba
    integer   :: i
    type(box) :: bs(nboxes)

    do i=1,nboxes
       bs(i) = make_box(boxes(:,1,i), boxes(:,2,i))
    end do

    call pybl_boxarray_new(cptr,ba)
    call build(ba, bs)
  end subroutine pybl_boxarray_create_from_boxes

  subroutine pybl_boxarray_print(cptr) bind(c)
    type(c_ptr), intent(in   ), value :: cptr
    type(boxarray), pointer :: ba
    call pybl_boxarray_get(cptr, ba)
    call print(ba)
  end subroutine pybl_boxarray_print

  subroutine pybl_boxarray_nboxes(cptr,nb) bind(c)
    type(c_ptr),    intent(in   ), value :: cptr
    integer(c_int), intent(  out)        :: nb
    type(boxarray), pointer :: ba
    call pybl_boxarray_get(cptr, ba)
    nb = nboxes(ba)
  end subroutine pybl_boxarray_nboxes

  subroutine pybl_boxarray_dim(cptr,dim) bind(c)
    type(c_ptr),    intent(in   ), value :: cptr
    integer(c_int), intent(  out)        :: dim
    type(boxarray), pointer :: ba
    call pybl_boxarray_get(cptr, ba)
    dim = ba%dim
  end subroutine pybl_boxarray_dim

  subroutine pybl_boxarray_maxsize(cptr,maxsize,ndim) bind(c)
    type(c_ptr),    intent(in   ), value :: cptr
    integer(c_int), intent(in   ), value :: ndim
    integer(c_int), intent(in   )        :: maxsize(ndim)
    type(boxarray), pointer :: ba
    call pybl_boxarray_get(cptr, ba)
    call boxarray_maxsize(ba, maxsize)
  end subroutine pybl_boxarray_maxsize

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! layout routines

  subroutine pybl_layout_create_from_boxarray(bacptr,dim,pmask,cptr) bind(c)
    type(c_ptr),    intent(in   ),  value :: bacptr
    type(c_ptr),    intent(  out)         :: cptr
    integer(c_int), intent(in   ),  value :: dim
    integer(c_int), intent(in   )         :: pmask(dim)

    type(boxarray), pointer :: ba
    type(layout), pointer   :: la
    logical                 :: lpmask(dim)

    call pybl_boxarray_get(bacptr, ba)
    call pybl_layout_new(cptr, la)

    lpmask = pmask == 1
    call build(la, ba, boxarray_bbox(ba), pmask=lpmask)
  end subroutine pybl_layout_create_from_boxarray

  subroutine pybl_layout_print(cptr) bind(c)
    type(c_ptr), intent(in   ), value :: cptr
    type(layout), pointer :: la
    call pybl_layout_get(cptr, la)
    call print(la)
  end subroutine pybl_layout_print

  subroutine pybl_layout_nboxes(cptr, boxes) bind(c)
    type(c_ptr),    intent(in   ), value :: cptr
    integer(c_int), intent(  out)        :: boxes
    type(layout), pointer :: la
    call pybl_layout_get(cptr,la)
    boxes = nboxes(la)
  end subroutine pybl_layout_nboxes

  subroutine pybl_layout_get_box(cptr, intbuf, i) bind(c)
    type(c_ptr),    intent(in   ), value :: cptr, intbuf
    integer(c_int), intent(in   ), value :: i
    type(layout), pointer :: la
    integer, pointer :: buf(:)
    type(box) :: lbx
    call c_f_pointer(intbuf, buf, [ 7 ])
    call pybl_layout_get(cptr,la)
    lbx = get_box(la, i)
    buf(1) = lbx%dim
    buf(2:4) = lbx%lo
    buf(5:7) = lbx%hi
  end subroutine pybl_layout_get_box

  subroutine pybl_layout_local(cptr, nbox, islocal) bind(c)
    type(c_ptr),    intent(in   ), value :: cptr
    integer(c_int), intent(in   ), value :: nbox
    integer(c_int), intent(  out)        :: islocal
    type(layout), pointer :: la

    call pybl_layout_get(cptr,la)
    islocal = 0
    if (local(la, nbox)) then
       islocal = 1
    end if
  end subroutine pybl_layout_local

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! multifab routines

  subroutine pybl_lmultifab_create_from_layout(la_cptr,cptr) bind(c)
    type(c_ptr),    intent(in   ), value :: la_cptr
    type(c_ptr),    intent(  out)        :: cptr

    type(layout), pointer   :: la
    type(lmultifab), pointer :: mfab

    call pybl_layout_get(la_cptr,la)
    call pybl_lmultifab_new(cptr,mfab)

    call build(mfab, la)
    call setval(mfab, .false.)
  end subroutine pybl_lmultifab_create_from_layout

  subroutine pybl_multifab_create_from_layout(la_cptr,nc,ng,cptr) bind(c)
    type(c_ptr),    intent(in   ), value :: la_cptr
    integer(c_int), intent(in   ), value :: nc, ng
    type(c_ptr),    intent(  out)        :: cptr

    type(layout), pointer   :: la
    type(multifab), pointer :: mfab

    call pybl_layout_get(la_cptr,la)
    call pybl_multifab_new(cptr,mfab)

    call build(mfab, la, nc=nc, ng=ng) !, interleave=interleave)
    call setval(mfab, 0.0d0)
  end subroutine pybl_multifab_create_from_layout

  subroutine pybl_multifab_print(cptr) bind(c)
    type(c_ptr), intent(in   ), value :: cptr
    type(multifab), pointer :: mfab
    call pybl_multifab_get(cptr, mfab)
    call print(mfab)
  end subroutine pybl_multifab_print

  subroutine pybl_multifab_destroy(cptr) bind(c)
    type(c_ptr), intent(in   ), value :: cptr
    type(multifab), pointer :: mfab
    call pybl_multifab_get(cptr, mfab)
    call destroy(mfab)
  end subroutine pybl_multifab_destroy

  subroutine pybl_multifab_fill_boundary(cptr) bind(c)
    type(c_ptr), intent(in   ), value :: cptr
    type(multifab), pointer :: mfab
    call pybl_multifab_get(cptr, mfab)
    call fill_boundary(mfab)
  end subroutine pybl_multifab_fill_boundary

  subroutine pybl_multifab_write(cptr, dirname, dlen, header, hlen) bind(c)
    type(c_ptr),       intent(in   ), value :: cptr
    integer(c_int),    intent(in   ), value :: dlen, hlen
    character(c_char), intent(in   )        :: dirname(dlen)
    character(c_char), intent(in   )        :: header(hlen)

    type(multifab), pointer :: mfab
    character(len=dlen) :: dname
    character(len=hlen) :: hname

    call pybl_cstring(dirname, dlen, dname)
    call pybl_cstring(header, hlen, hname)

    call pybl_multifab_get(cptr, mfab)
    call fabio_write(mfab, dname, hname)
  end subroutine pybl_multifab_write

  subroutine pybl_multifab_read(dirname, dlen, header, hlen, cptr) bind(c)
    integer(c_int),    intent(in   ), value :: dlen, hlen
    character(c_char), intent(in   )  :: dirname(dlen)
    character(c_char), intent(in   )  :: header(hlen)
    type(c_ptr),       intent(  out) :: cptr

    type(multifab), pointer :: mfab
    character(len=dlen) :: dname
    character(len=hlen) :: hname

    call pybl_cstring(dirname, dlen, dname)
    call pybl_cstring(header, hlen, hname)
    call pybl_multifab_new(cptr, mfab)
    call fabio_multifab_read_d(mfab, dname, hname)
  end subroutine pybl_multifab_read

  subroutine pybl_multifab_copy(dcptr, scptr) bind(c)
    type(c_ptr), intent(in   ), value :: dcptr, scptr
    type(multifab), pointer :: dmfab, smfab
    call pybl_multifab_get(dcptr, dmfab)
    call pybl_multifab_get(scptr, smfab)
    call copy(dmfab, 1, smfab, 1, ncomp(smfab))
  end subroutine pybl_multifab_copy

  subroutine pybl_multifab_info(cptr, dim, nboxes, nc, ng) bind(c)
    type(c_ptr),    intent(in   ), value :: cptr
    integer(c_int), intent(  out)        :: dim, nboxes, nc, ng
    type(multifab), pointer :: mfab
    call pybl_multifab_get(cptr,mfab)
    dim = mfab%dim
    nboxes = nfabs(mfab)
    nc = mfab%nc
    ng = mfab%ng
  end subroutine pybl_multifab_info

  subroutine pybl_multifab_layout(cptr, laptr) bind(c)
    type(c_ptr), intent(in   ), value :: cptr
    type(c_ptr), intent(  out)        :: laptr
    type(multifab), pointer :: mfab
    type(layout), pointer :: la
    call pybl_multifab_get(cptr,mfab)
    la => mfab%la
    laptr = c_loc(la)
  end subroutine pybl_multifab_layout

  subroutine pybl_lmultifab_info(cptr, dim, nboxes, nc, ng) bind(c)
    type(c_ptr),    intent(in   ), value :: cptr
    integer(c_int), intent(  out)        :: dim, nboxes, nc, ng
    type(lmultifab), pointer :: mfab
    call pybl_lmultifab_get(cptr,mfab)
    dim = mfab%dim
    nboxes = nfabs(mfab)
    nc = mfab%nc
    ng = mfab%ng
  end subroutine pybl_lmultifab_info

  subroutine pybl_fab_info(cptr, nbox, dim, nc, bx_lo, bx_hi, pbx_lo, pbx_hi, ibx_lo, ibx_hi) bind(c)
    type(c_ptr),            intent(in   ), value        :: cptr
    integer(c_int),         intent(in   ), value        :: nbox
    integer(c_int),         intent(  out)               :: dim, nc
    integer(c_int),         intent(  out), dimension(3) :: bx_lo, bx_hi, pbx_lo, pbx_hi, ibx_lo, ibx_hi

    type(multifab), pointer :: mfab
    type(box) :: bx

    call pybl_multifab_get(cptr,mfab)
    dim = get_dim(mfab%fbs(nbox))
    nc = ncomp(mfab%fbs(nbox))

    bx = get_box(mfab%fbs(nbox))
    bx_lo = bx%lo
    bx_hi = bx%hi

    bx = get_pbox(mfab%fbs(nbox))
    pbx_lo = bx%lo
    pbx_hi = bx%hi

    bx = get_ibox(mfab%fbs(nbox))
    ibx_lo = bx%lo
    ibx_hi = bx%hi
  end subroutine pybl_fab_info

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! regrid

  subroutine pybl_regrid(lacptrs, mfcptrs, nlevs, max_levs, tag_boxes_cb) bind(c)
    use bc_module
    use define_bc_module
    use regrid_module

    integer(c_int), intent(in   ), value :: max_levs
    integer(c_int), intent(inout)        :: nlevs
    type(c_ptr),    intent(inout)        :: lacptrs(max_levs), mfcptrs(max_levs)
    type(c_ptr),    intent(in   ), value :: tag_boxes_cb

    type(multifab), pointer :: mf
    type(layout), pointer   :: la
    type(multifab), target  :: mfs(max_levs)

    integer           :: i, n, dim, new_nlevs, phys_bc_in(3,2)
    type(ml_boxarray) :: mba
    type(ml_layout)   :: mla
    type(bc_tower)    :: bct
    real(dp_t)        :: dx(max_levs)

    ! build multilevel layout
    call pybl_layout_get(lacptrs(1),la)
    dim = get_dim(la)

    call build(mba, max_levs, dim)
    mba%rr = 2
    dx = 1
    do i=1,max_levs-1
       dx(i+1) = dx(i)/2
    end do

    mba%pd(1) = get_box(la, 1)
    do n=2,max_levs
       mba%pd(n) = refine(mba%pd(n-1),mba%rr((n-1),:))
    enddo

    do i=1,nlevs
       call pybl_layout_get(lacptrs(i),la)
       call pybl_multifab_get(mfcptrs(i),mf)
       mba%bas(i) = get_boxarray(la)
       call boxarray_verify_dim(mba%bas(i))
       mfs(i) = mf
    end do

    call ml_layout_restricted_build(mla,mba,nlevs)
    ! call destroy(mba)

    phys_bc_in = PERIODIC
    call bc_tower_init(bct, max_levs, dim, phys_bc_in)
    do n = 1,nlevs
       call bc_tower_level_build(bct,n,mla%la(n))
    end do

    call regrid(mla, mfs, nlevs, max_levs, dx, bct, 2, 64, tag_boxes_cb)

    do i=1,nlevs
       call pybl_multifab_new(mfcptrs(i),mf)
       mf = mfs(i)
    end do

  end subroutine pybl_regrid

end module fboxlib
