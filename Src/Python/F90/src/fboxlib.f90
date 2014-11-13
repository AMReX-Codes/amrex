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

  subroutine pybl_get_multifab_info(cptr, dim, nboxes, nc, ng) bind(c)
    type(c_ptr),    intent(in   ), value :: cptr
    integer(c_int), intent(  out)        :: dim, nboxes, nc, ng
    type(multifab), pointer :: mfab
    call pybl_multifab_get(cptr,mfab)
    dim = mfab%dim
    nboxes = nfabs(mfab)
    nc = mfab%nc
    ng = mfab%ng
  end subroutine pybl_get_multifab_info

  subroutine pybl_get_multifab_fab_info(cptr, nbox, dim, nc, bx_lo, bx_hi, pbx_lo, pbx_hi, ibx_lo, ibx_hi) bind(c)
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
  end subroutine pybl_get_multifab_fab_info

  subroutine pybl_create_multifab_from_layout(la_cptr,nc,ng,cptr) bind(c)
    type(c_ptr),    intent(in   ), value :: la_cptr
    integer(c_int), intent(in   ), value :: nc, ng
    type(c_ptr),    intent(  out)        :: cptr

    type(layout), pointer   :: la
    type(multifab), pointer :: mfab

    call pybl_layout_get(la_cptr,la)
    call pybl_multifab_new(cptr,mfab)

    call build(mfab, la, nc=nc, ng=ng) !, interleave=interleave)
    call setval(mfab, 0.0d0)
  end subroutine pybl_create_multifab_from_layout

  subroutine pybl_create_multifab_from_bbox(cptr1,nc,ng,cptr) bind(c)
    type(c_ptr),    intent(in   ), value :: cptr1
    integer(c_int), intent(in   )        :: nc, ng
    type(c_ptr),    intent(  out)        :: cptr

    type(multifab), pointer :: mfab1, mfab
    type(layout)   :: la
    type(boxarray) :: ba

    call pybl_multifab_get(cptr1, mfab1)
    call build(ba, boxarray_bbox(get_boxarray(get_layout(mfab1))))
    call build(la, ba, boxarray_bbox(ba))
    call pybl_multifab_new(cptr,mfab)
    call build(mfab, la, nc=nc, ng=ng) !, interleave=interleave)
    call setval(mfab, 0.0d0)
  end subroutine pybl_create_multifab_from_bbox

  subroutine pybl_print_multifab(cptr) bind(c)
    type(c_ptr), intent(in   ), value :: cptr
    type(multifab), pointer :: mfab
    call pybl_multifab_get(cptr, mfab)
    call print(mfab)
  end subroutine pybl_print_multifab

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

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! plotfile

  subroutine pybl_create_plotfile(dname, dlen, cptr) bind(c)
    use bl_io_module
    use plotfile_module
    type(c_ptr),    intent(in   ), value :: dname
    integer(c_int), intent(in   ), value :: dlen
    type(c_ptr),    intent(  out)        :: cptr

    type(plotfile), pointer :: pf
    character(len=dlen), pointer :: root
    integer :: un

    call c_f_pointer(dname, root)
    call pybl_plotfile_new(cptr, pf)
    un = unit_new()
    call build(pf, root, un)
  end subroutine pybl_create_plotfile

  subroutine pybl_get_plotfile_info(cptr, dim, nvars, flevel) bind(c)
    use plotfile_module, only: plotfile
    type(c_ptr),    intent(in   ), value :: cptr
    integer(c_int), intent(  out)        :: dim, nvars, flevel

    type(plotfile), pointer :: pf
    call pybl_plotfile_get(cptr, pf)
    dim    = pf%dim
    nvars  = pf%nvars
    flevel = pf%flevel
  end subroutine pybl_get_plotfile_info

  subroutine pybl_get_plotfile_grid_info(cptr, level, nboxes) bind(c)
    use plotfile_module, only: plotfile
    type(c_ptr),    intent(in   ), value :: cptr
    integer(c_int), intent(in   ), value :: level
    integer(c_int), intent(  out)        :: nboxes
    type(plotfile), pointer :: pf
    call pybl_plotfile_get(cptr, pf)
    nboxes = plotfile_nboxes_n(pf, level)
  end subroutine pybl_get_plotfile_grid_info

  subroutine pybl_get_plotfile_name(cptr, nvar, nlen, nameptr) bind(c)
    use plotfile_module
    type(c_ptr),    intent(in   ), value :: cptr
    integer(c_int), intent(in   ), value :: nvar, nlen
    type(c_ptr),    intent(in   ), value :: nameptr

    type(plotfile), pointer :: pf
    character(len=nlen), pointer :: name

    call pybl_plotfile_get(cptr, pf)
    call c_f_pointer(nameptr, name)
    name = pf%names(nvar)
  end subroutine pybl_get_plotfile_name

  subroutine pybl_plotfile_bind(cptr, i, j, c) bind(c)
    use plotfile_module
    type(c_ptr),    intent(in), value :: cptr
    integer(c_int), intent(in), value :: i, j, c

    type(plotfile), pointer :: pf

    call pybl_plotfile_get(cptr, pf)
    call fab_bind_comp_vec(pf, i, j, (/ c /) )
  end subroutine pybl_plotfile_bind

  subroutine pybl_plotfile_unbind(cptr, i, j) &
       bind(c, name='pybl_plotfile_unbind')
    use plotfile_module
    type(c_ptr),    intent(in), value :: cptr
    integer(c_int), intent(in), value :: i, j

    type(plotfile), pointer :: pf

    call pybl_plotfile_get(cptr, pf)
    call fab_unbind(pf, i, j)
  end subroutine pybl_plotfile_unbind

  subroutine pybl_get_plotfile_fab_info(cptr, level, nbox, dim, nc, &
       pbx_lo, pbx_hi) &
       bind(c, name='pybl_get_plotfile_fab_info')
    use plotfile_module
    implicit none
    type(c_ptr), intent(in), value :: cptr
    integer(c_int), intent(in), value :: level, nbox
    integer(c_int), intent(out) :: dim, nc
    integer(c_int), intent(out), dimension(3) :: pbx_lo, pbx_hi

    type(plotfile), pointer :: pf
    type(box) :: bx
    call pybl_plotfile_get(cptr, pf)
    dim = pf%dim
    nc  = pf%nvars

    bx = get_box(pf, level, nbox)
    pbx_lo = bx%lo
    pbx_hi = bx%hi
  end subroutine pybl_get_plotfile_fab_info

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! cluster/regrid

  subroutine pybl_cluster(tags_cptr, buffer_width, boxes_cptr) bind(c)
    use cluster_module
    implicit none

    type(c_ptr),    intent(in), value :: tags_cptr
    integer(c_int), intent(in), value :: buffer_width
    type(c_ptr),    intent(out) :: boxes_cptr

    type(lmultifab), pointer :: tags
    type(boxarray),  pointer :: boxes

    call pybl_boxarray_new(boxes_cptr, boxes)
    call pybl_lmultifab_get(tags_cptr, tags)

    call cluster(boxes, tags, buffer_width)
  end subroutine pybl_cluster

  subroutine pybl_regrid(lacptrs, mfcptrs, nlevs, max_levs, tag_boxes_cb) bind(c)
    use bc_module
    use define_bc_module
    use regrid_module

    integer(c_int), intent(in   ), value :: max_levs
    integer(c_int), intent(inout)        :: nlevs
    type(c_ptr),    intent(in   )        :: lacptrs(max_levs), mfcptrs(max_levs)
    type(c_ptr),    intent(in   ), value :: tag_boxes_cb

    type(multifab), pointer :: mf
    type(layout), pointer   :: la

    integer           :: i, n, dim, new_nlevs, phys_bc_in(3,2)
    type(multifab)    :: mfs(max_levs)
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
  end subroutine pybl_regrid

end module fboxlib
