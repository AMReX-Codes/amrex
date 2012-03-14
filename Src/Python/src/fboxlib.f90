module fboxlib

  use iso_c_binding
  use blobjects
  implicit none

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! open, close, set_comm

  subroutine pybl_hello() bind(c, name='pybl_hello')
    print *, 'Hello from Fortran, with FBoxLib!'
  end subroutine pybl_hello

  subroutine pybl_open() bind(c, name='pybl_open')
    use parallel
    call parallel_initialize()
  end subroutine pybl_open

  subroutine pybl_close() bind(c, name='pybl_close')
    use parallel
    call parallel_finalize()
  end subroutine pybl_close

  subroutine pybl_mpi_rank(r) bind(c, name='pybl_mpi_rank')
    use parallel
    integer(c_int), intent(out) :: r
    r = parallel_myproc()
  end subroutine pybl_mpi_rank

  subroutine pybl_mpi_size(r) bind(c, name='pybl_mpi_size')
    use parallel
    integer(c_int), intent(out) :: r
    r = parallel_nprocs()
  end subroutine pybl_mpi_size

  ! subroutine pybl_set_comm(comm)
  !   use parallel
  !   implicit none
  !   integer, intent(in) :: comm
  !   call parallel_set_m_comm(comm)
  ! end subroutine pybl_set_comm

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! boxarray routines

  subroutine pybl_create_boxarray_from_boxes(boxes,nboxes,dim,cptr) &
       bind(c, name='pybl_create_boxarray_from_boxes')
    implicit none
    integer(c_int), intent(in), value  :: dim, nboxes
    integer(c_int), intent(in) :: boxes(nboxes,2,dim)
    type(c_ptr), intent(out) :: cptr

    integer :: i
    type(box) :: bs(nboxes)
    type(boxarray), pointer :: ba

    do i=1,nboxes
       bs(i) = make_box(boxes(i,1,:), boxes(i,2,:))
    end do

    call pybl_boxarray_new(cptr,ba)
    call build(ba, bs)
  end subroutine pybl_create_boxarray_from_boxes

  subroutine pybl_print_boxarray(cptr) bind(c, name='pybl_print_boxarray')
    implicit none
    type(c_ptr), intent(in), value :: cptr
    type(boxarray), pointer :: ba

    call pybl_boxarray_get(cptr, ba)
    call print(ba)
  end subroutine pybl_print_boxarray

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! layout routines

  subroutine pybl_create_layout_from_boxarray(ba_cptr,cptr) &
       bind(c, name='pybl_create_layout_from_boxarray')
    implicit none
    type(c_ptr), intent(in), value :: ba_cptr
    type(c_ptr), intent(out) :: cptr

    type(boxarray), pointer :: ba
    type(layout), pointer :: la

    call pybl_boxarray_get(ba_cptr, ba)
    call pybl_layout_new(cptr, la)

    call build(la, ba)
  end subroutine pybl_create_layout_from_boxarray

  subroutine pybl_create_layout_from_boxes(boxes,nboxes,dim,cptr) &
       bind(c, name='pybl_create_layout_from_boxes')
    implicit none
    integer(c_int), intent(in), value :: dim, nboxes
    integer(c_int), intent(in) :: boxes(nboxes,2,dim)
    type(c_ptr), intent(out) :: cptr

    integer :: i
    type(box) :: bs(nboxes)
    type(boxarray) :: ba
    type(layout), pointer :: la

    do i=1,nboxes
       bs(i) = make_box(boxes(i,1,:), boxes(i,2,:))
    end do

    call pybl_layout_new(cptr,la)

    call build(ba, bs)
    call build(la, ba)
  end subroutine pybl_create_layout_from_boxes

  subroutine pybl_create_ml_layout_from_layouts(lacptrs,nlevels,cptr) &
       bind(c, name='pybl_create_ml_layout_from_layouts')
    implicit none
    integer(c_int), intent(in), value  :: nlevels
    type(c_ptr), intent(in) :: lacptrs(nlevels)
    type(c_ptr), intent(out) :: cptr

    integer :: i, dim
    type(ml_boxarray) :: mba
    type(ml_layout), pointer :: mla
    type(layout), pointer :: la
    type(boxarray) :: bas

    call pybl_ml_layout_new(cptr,mla)
    call pybl_layout_get(lacptrs(1),la)

    dim = get_dim(la)

    call build(mba, nlevels, dim)

    do i=1,nlevels
       call pybl_layout_get(lacptrs(i),la)
       mba%bas(i) = get_boxarray(la)
    end do

    call build(mla,mba)
  end subroutine pybl_create_ml_layout_from_layouts

  subroutine pybl_print_layout(cptr) &
       bind(c, name='pybl_print_layout')
    implicit none
    type(c_ptr), intent(in), value :: cptr
    type(layout), pointer :: la

    call pybl_layout_get(cptr, la)
    call print(la)
  end subroutine pybl_print_layout

  subroutine pybl_print_mllayout(cptr) &
       bind(c, name='pybl_print_mllayout')
    implicit none
    type(c_ptr), intent(in), value :: cptr
    type(ml_layout), pointer :: mla

    call pybl_ml_layout_get(cptr, mla)
    call print(mla)
  end subroutine pybl_print_mllayout

  subroutine pybl_layout_nboxes(cptr, boxes) &
       bind(c, name='pybl_layout_nboxes')
    implicit none
    type(c_ptr),    intent(in), value :: cptr
    integer(c_int), intent(out)       :: boxes

    type(layout), pointer :: la

    call pybl_layout_get(cptr,la)

    boxes = nboxes(la)
  end subroutine pybl_layout_nboxes

  subroutine pybl_layout_local(cptr, nbox, islocal) &
       bind(c, name='pybl_layout_local')
    implicit none
    type(c_ptr),    intent(in), value :: cptr
    integer(c_int), intent(in), value :: nbox
    logical,        intent(out)       :: islocal

    type(layout), pointer :: la

    call pybl_layout_get(cptr,la)

    islocal = local(la, nbox)
  end subroutine pybl_layout_local

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! multifab routines

  subroutine pybl_get_multifab_info(cptr, dim, nboxes, nc, ng) &
       bind(c, name='pybl_get_multifab_info')
    implicit none
    type(c_ptr), intent(in), value :: cptr
    integer(c_int), intent(out) :: dim, nboxes, nc, ng
    type(multifab), pointer :: mfab

    call pybl_multifab_get(cptr,mfab)
    dim = mfab%dim
    nboxes = mfab%nboxes
    nc = mfab%nc
    ng = mfab%ng
  end subroutine pybl_get_multifab_info

  subroutine pybl_get_multifab_fab_info(cptr, nbox, dim, nc, &
       bx_lo, bx_hi, pbx_lo, pbx_hi, ibx_lo, ibx_hi) &
       bind(c, name='pybl_get_multifab_fab_info')
    implicit none
    type(c_ptr), intent(in), value :: cptr
    integer(c_int), intent(in), value :: nbox
    integer(c_int), intent(out) :: dim, nc
    integer(c_int), intent(out), dimension(3) :: bx_lo, bx_hi, &
         pbx_lo, pbx_hi, ibx_lo, ibx_hi
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

  subroutine pybl_create_multifab_from_layout(la_cptr,nc,ng,interleave,cptr) &
       bind(c, name='pybl_create_multifab_from_layout')
    implicit none
    type(c_ptr),    intent(in), value :: la_cptr
    integer(c_int), intent(in), value :: nc, ng
    logical,        intent(in), value :: interleave
    type(c_ptr),    intent(out)       :: cptr

    type(layout), pointer :: la
    type(multifab), pointer :: mfab

    call pybl_layout_get(la_cptr,la)
    call pybl_multifab_new(cptr,mfab)

    call build(mfab, la, nc=nc, ng=ng) !, interleave=interleave)
    call setval(mfab, 0.0d0)
  end subroutine pybl_create_multifab_from_layout

  subroutine pybl_create_multifab_from_bbox(cptr1, nc,ng,interleave, cptr) &
       bind(c, name='pybl_create_multifab_from_bbox')
    implicit none
    type(c_ptr), intent(in), value :: cptr1
    integer(c_int), intent(in) :: nc, ng
    logical, intent(in)  :: interleave
    type(c_ptr), intent(out) :: cptr
    type(multifab), pointer :: mfab1, mfab
    type(layout) :: la
    type(boxarray) :: ba
    type(box) :: bx

    call pybl_multifab_get(cptr1, mfab1)
    call build(ba, boxarray_bbox(get_boxarray(get_layout(mfab1))))
    call build(la, ba)

    call pybl_multifab_new(cptr,mfab)

    call build(mfab, la, nc=nc, ng=ng) !, interleave=interleave)
    call setval(mfab, 0.0d0)
  end subroutine pybl_create_multifab_from_bbox

  subroutine pybl_print_multifab(cptr) &
       bind(c, name='pybl_print_multifab')
    implicit none
    type(c_ptr), intent(in), value :: cptr
    type(multifab), pointer :: mfab

    call pybl_multifab_get(cptr, mfab)
    call print(mfab)
  end subroutine pybl_print_multifab

  subroutine pybl_multifab_fill_boundary(cptr) &
       bind(c, name='pybl_multifab_fill_boundary')
    implicit none
    type(c_ptr), intent(in), value :: cptr
    type(multifab), pointer :: mfab

    call pybl_multifab_get(cptr, mfab)

    call fill_boundary(mfab, 1, ncomp(mfab))
  end subroutine pybl_multifab_fill_boundary

  ! subroutine pybl_multifab_write(cptr, dirname, header) &
  !      bind(c, name='pybl_multifab_write')
  !   use fabio_module
  !   implicit none
  !   type(c_ptr), intent(in), value :: cptr
  !   ! character(len=*), intent(in) :: dirname, header
  !   type(multifab), pointer :: mfab

  !   call pybl_multifab_get(cptr, mfab)

  !   call fabio_write(mfab, dirname, header)
  ! end subroutine pybl_multifab_write

  ! subroutine pybl_multifab_read(dirname, header, cptr) &
  !      bind(c, name='pybl_multifab_read')
  !   use fabio_module
  !   implicit none
  !   character(len=*), intent(in) :: dirname, header
  !   type(multifab), pointer :: mfab
  !   type(c_ptr), intent(out) :: cptr

  !   call pybl_multifab_new(cptr, mfab)

  !   call fabio_multifab_read_d(mfab, dirname, header)
  ! end subroutine pybl_multifab_read

  subroutine pybl_multifab_copy(dcptr, scptr) &
       bind(c, name='pybl_multifab_copy')
    use fabio_module
    implicit none
    type(c_ptr), intent(in), value :: dcptr, scptr
    type(multifab), pointer :: dmfab, smfab

    call pybl_multifab_get(dcptr, dmfab)
    call pybl_multifab_get(scptr, smfab)

    call copy(dmfab, 1, smfab, 1, ncomp(smfab))
  end subroutine pybl_multifab_copy


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! regrid

  subroutine pybl_regrid(tags_cptr, buffer_width, boxes_cptr) &
       bind(c, name='pybl_regrid')
    use cluster_module
    implicit none

    type(c_ptr), intent(in), value :: tags_cptr
    integer, intent(in) :: buffer_width
    type(c_ptr), intent(out) :: boxes_cptr

    type(lmultifab), pointer :: tags
    type(boxarray), pointer :: boxes

    call pybl_boxarray_new(boxes_cptr, boxes)
    call pybl_lmultifab_get(tags_cptr, tags)

    call cluster(boxes, tags, buffer_width)
  end subroutine pybl_regrid

end module fboxlib
