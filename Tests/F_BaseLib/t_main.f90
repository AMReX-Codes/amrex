subroutine t_mt_random_numbers
  use mt19937_module
  implicit none
  integer n
  integer, allocatable :: gt(:)
  real a

  print *, 'MT_VALIDATE() = ', mt_validate()
  call mt_random_seed(size = n)
  print *, 'MT_RANDOM_SEED(SIZE) ', n
  allocate(gt(n))
  print *, 'MT_RANDOM_SEED(GET)'
  call mt_random_seed(get = gt)
  call mt_random_number(a)
  print *, 'MT_RANDOM_NUMBER(a) = ', a
  call mt_random_number(a)
  print *, 'MT_RANDOM_NUMBER(a) = ', a
  print *, 'MT_RANDOM_SEED(PUT)'
  call mt_random_seed(put = gt)
  call mt_random_number(a)
  print *, 'MT_RANDOM_NUMBER(a) = ', a
  call mt_random_number(a)
  print *, 'MT_RANDOM_NUMBER(a) = ', a
end subroutine t_mt_random_numbers

subroutine t_plotfile
  use f2kcli
  use plotfile_module
  use bl_IO_module
  use bl_error_module
  implicit none
  type(plotfile) pf
  integer :: unit
  integer :: i
  character(len=256) :: pfname

  if ( command_argument_count() /= 1 ) &
       call bl_error('number of args')
  call get_command_argument(1, pfname)

  unit = unit_new()
  call build(pf, pfname, unit)
  write(*,'("PLOTFILE: ", a)') pfname
  write(*,'("Time               : ", g15.10)') plotfile_time(pf)
  write(*,'("Number of Levels     ", i15)') plotfile_nlevels(pf)
  if ( plotfile_nlevels(pf) > 0 ) then
     write(*,'("Refinement ratio   : ")', advance='no')
     do i = 1, plotfile_nlevels(pf)-1
        write(*, '(i5)', advance='no') plotfile_refrat_n(pf, i)
     end do
     write(*,fmt='()')
  end if
  write(*,'("Number of Variables: ", i15)') plotfile_nvars(pf)
  do i = 1, plotfile_nvars(pf)
     write(*,'(i3,1x, a20)') i, adjustr(plotfile_var_name(pf, i))
  end do
  do i = 1, plotfile_nlevels(pf)
     write(*,*) 'level ', i, ' has ', plotfile_nboxes_n(pf, i), ' grids'
  end do
  call destroy(pf)

end subroutine t_plotfile

subroutine t_fabio
  use fabio_module
  use fab_module
  use multifab_module
  implicit none
  integer :: fd
  type(fab) :: fb
  type(box) :: bx
  integer :: offset
  integer :: sz, i, j
  real(kind=dp_t), pointer :: fbp(:,:,:,:)

  sz = 32
  bx = refine(unit_box(2), sz)
  call fab_build(fb, bx)
  fbp => dataptr(fb)
  do j = 0, sz-1
     do i = 0, sz-1
        fbp(i,j,1,1) = 1000*i + j
     end do
  end do

  call fabio_mkdir("tdir")
  call fabio_open(fd, "tdir/foo", FABIO_WRONLY)

  call fabio_write(fd, offset, fb)

  call destroy(fb)

  print *, 'offset = ', offset
  
end subroutine t_fabio

subroutine t_boxassoc

  use layout_module
  type(boxassoc) :: bxasc
  type(boxarray) :: ba
  integer :: dm
  integer :: ng
  integer :: pd_xyz(MAX_SPACEDIM)
  integer :: ba_maxsize
  type(layout) :: la
  logical :: nodal(MAX_SPACEDIM)
  type(box) :: bxs(4)
  logical pmask(MAX_SPACEDIM)
  type(box) :: pd

  dm = 2
  ng = 1
  pd_xyz     = 32
  ba_maxsize = 16

  nodal = .true.
  nodal = .false.
  pmask = .true.

  if ( .false. ) then
     bxs(1) = make_box((/0,0/),(/3,3/))
     bxs(2) = make_box((/4,4/),(/7,7/))
     bxs(3) = make_box((/0,8/),(/7,16/))
     call build(ba, bxs(1:2))
  else
     call build(ba, make_box((/(0,i=1,dm)/), pd_xyz(1:dm)-1))
     call boxarray_maxsize(ba, ba_maxsize)
  end if

  pd = bbox(ba)
  call build(la, ba, pd = pd, pmask = pmask(1:dm))

  call boxassoc_build(bxasc, la%lap, ng, nodal = nodal(1:dm))

  call boxassoc_print(bxasc)

  call boxassoc_destroy(bxasc)
  call destroy(la)
  call boxarray_destroy(ba)

end subroutine t_boxassoc

subroutine t_box_conn
  use mboxarray_module
  use box_util_module
  implicit none
  type(mboxarray) :: mba
  integer dm, i
  integer :: ml
  real(kind=dp_t) :: d

  call mboxarray_read_boxes(mba, "../../data/grids/3D_4_level_96x96x96")

  call box_conn(mba%bas(4))

  dm = mba%dim
  do i = 1, mba%nlevel
     ml =  local_max_mg_levels(mba%bas(i),1)
     print *, i, "mg_levels: ", ml
     d = boxarray_dvolume(mba%bas(i)) 
     print *, i, ": ", d, " => ", d/(2**3)**(ml-1)
  end do

contains
  subroutine box_conn(ba)
    use bl_IO_module
    type(boxarray), intent(in) :: ba
    integer :: i, j
    type(box):: b1, b2
    integer cnt, vol
    integer :: un
    un = unit_new()
    open(un,file='conn', status = 'replace', action = 'write')
    write(un,fmt='(1x,i10)') ba%nboxes
    do i = 1, ba%nboxes
       b1 = grow(ba%bxs(i),1)
       cnt = 0
       do j = 1, ba%nboxes
          if ( i == j ) cycle
          b2 = ba%bxs(j)
          if ( intersects(b1, b2) ) then
             cnt = cnt + 1
          end if
       end do
       vol = box_dvolume(ba%bxs(i))
       write(un,fmt='(1x,i10,1x,i10,1x,i10)', advance='no') i, vol, cnt
       do j = 1, ba%nboxes
          if ( i == j ) cycle
          b2 = ba%bxs(j)
          if ( intersects(b1, b2) ) then
             write(un,fmt='(1x,i5)',advance='no') j
          end if
       end do
       write(un,fmt='()')
    end do
    close(un)
  end subroutine box_conn
  function local_max_mg_levels(ba, min_size) result(r)
    type(boxarray), intent(in) :: ba
    integer, intent(in), optional :: min_size
    integer :: r
    integer, parameter :: rrr = 2
    type(box) :: bx, bx1
    integer :: i, rr, lmn, dm
    lmn = 1; if ( present(min_size) ) lmn = min_size
    r = 1
    rr = rrr
    dm = ba%dim
    do
       do i = 1, size(ba%bxs)
          bx = ba%bxs(i)
          bx1 = coarsen(bx, rr)
          if ( any(extent(bx1) < lmn) ) return
          if ( bx /= refine(bx1, rr) ) then
             return
          end if
       end do
       rr = rr*rrr
       r  = r + 1
    end do
  end function local_max_mg_levels
end subroutine t_box_conn

subroutine t_boxarray
  use boxarray_module
  type(boxarray) :: bao, ba
  type(box) :: pd
  integer :: dm
  integer :: pd_xyz(MAX_SPACEDIM)
  integer :: ba_maxsize

  dm = 2
  pd_xyz     = 32
  ba_maxsize = 32

  call box_build_2(pd, (/(0,i=1,dm)/), pd_xyz(1:dm)-1)
  call boxarray_build_bx(ba, pd)
  call boxarray_maxsize(ba, ba_maxsize)
  call boxarray_print(ba, "boxarray")

  call boxarray_boundary_n_d_f(bao, ba, 1, 1, 1)

  call boxarray_print(bao, "boundary")

end subroutine t_boxarray

subroutine t_mf
  use multifab_module
  use bl_IO_module
  use fabio_module
  implicit none
  integer, parameter :: dm = 2
  type (box) :: bx, pd
  type (box), dimension(4) :: bxs
  type (boxarray):: ba
  type (layout) :: la
  type (multifab) :: mf
  type (boxassoc) :: bxasc
  logical :: pmask(2)
  integer :: ms
  integer :: sz
  integer :: i

  sz = 8
  ms = 4

!  pmask = .true.
pmask = .false.
  if ( dm > 1 ) pmask(2) = .false.

  if ( .false. ) then
     if ( .false. ) then
        bx = refine(unit_box(dim=2),2)
        bxs(1) = bx
        bxs(2) = shift(bx,2,dim=1)
        bxs(3) = shift(bx,2,dim=2)
        bxs(4) = shift(bx,2)
        call build(ba, bxs)
     else
        call build(ba, make_box((/(0,i=1,dm)/), (/(sz-1,i=1,dm)/)))
        call boxarray_maxsize(ba, ms)
     end if
     pd = bbox(ba)
  else
     bxs(1) = make_box((/160, 48/), (/215,71/))
     bxs(2) = make_box((/152, 72/), (/231,103/))
     call build(ba, bxs(1:2))
     pd = make_box((/0,0/), (/511,255/))
  end if
  call parallel_barrier()
  call boxarray_print(ba, "BOXARRAY")
  call build(la, ba, pd = pd, pmask = pmask)
  call parallel_barrier()
  bxasc = layout_boxassoc(la, 1)
  call parallel_barrier()
  call boxassoc_print(bxasc, "BOXASSOC")

  call multifab_build(mf, la, nc = 1, ng = 1)
  call multifab_debug_fill(mf, loc = .true.)
  call parallel_barrier()
!  call mf_print(mf, "before")
  call multifab_fill_boundary(mf)
  call parallel_barrier()
!  call mf_print(mf, "after")
  call destroy(mf)
  call destroy(la)
  call destroy(ba)

contains
  
  subroutine mf_print(mf, str)
    type(multifab), intent(in) :: mf
    character(len=*), intent(in) :: str
    integer n
    real(kind=dp_t), pointer :: fp(:,:,:,:)
    type(box) :: pbx, bx
    print *, "MF: ", str
    do n = 1, mf%nboxes; if ( remote(mf, n) ) cycle
       fp => dataptr(mf, n)
       bx = get_box(mf, n)
       pbx = get_pbox(mf, n)
       call array_print_2d(fp(:,:,1,1), 6, lwb(pbx), lwb(bx), upb(bx))
    end do
  end subroutine mf_print

  subroutine array_print_2d(ar, unit, plo, lo, hi)
    use bl_types
    integer, intent(in) :: plo(:), lo(:), hi(:)
    real(kind=dp_t), intent(in) :: ar(plo(1):,plo(2):)
    integer, intent(in) :: unit
    integer i, j
    integer ext(2), wid
    character c
    ext = hi(1:2)-lo(1:2) + 1
    wid = size(ar,dim=1)*10
    do j = ubound(ar,2), lbound(ar,2), -1
       if ( j == hi(2) .or. j == lo(2)-1 ) then
          write(unit=unit,fmt='(a)') repeat('-',wid)
       end if
       do i = lbound(ar,1), ubound(ar,1)
          c = ' '
          if ( i == lo(1)-1 .or. i == hi(1) ) c = '|'
          write(unit=unit,fmt='(f9.4,a1)',advance = 'no') ar(i,j), c
       end do
       write(unit=unit,fmt='()')
    end do
    write(unit=unit,fmt='()')
  end subroutine array_print_2d
end subroutine t_mf

subroutine t_mf_fabio
  use multifab_module
  use bl_IO_module
  use fabio_module
  implicit none
  type (box) :: bx
  type (box), dimension(4) :: bxs
  type (boxarray):: ba
  type (layout) :: la
  type (multifab) :: mf
  type (boxassoc) :: bxasc
  integer :: rrs(1)

  bx = refine(unit_box(dim=2),2)
  bxs(1) = bx
  bxs(2) = shift(bx,2,dim=1)
  bxs(3) = shift(bx,2,dim=2)
  bxs(4) = shift(bx,2)
  call build(ba, bxs)
  call boxarray_print(ba, "BOXARRAY")
  call build(la, ba)
  bxasc = layout_boxassoc(la, 1)
  call boxassoc_print(bxasc, "BOXASSOC")
  call multifab_build(mf, la, nc = 1, ng=1)
  call setval(mf, -1.0_dp_t, ALL=.True.)
  call multifab_debug_fill(mf, loc = .True.)
!   do n = 1, mf%nboxes; if ( remote(mf, n) ) cycle
!      fp => dataptr(mf, n, get_ibox(mf, n))
!      nx = size(fp,1)
!      ny = size(fp,2)
!      do j = 1, ny
!         do i = 1, nx
!            fp(i,j,1,1) = (i-1) + 10*(j-1)
!         end do
!      end do
!   end do
  call print(mf, "before")
  call fabio_multifab_write_d(mf, "tdir", "flan")
  call fabio_ml_multifab_write_d((/mf/), rrs(1:0), "tdir1")
  call destroy(mf)
  call destroy(la)
  call destroy(ba)

end subroutine t_mf_fabio

subroutine t_nodal_mf_fabio
  use multifab_module
  use bl_IO_module
  use fabio_module
  implicit none
  type (box) :: bx
  type (box), dimension(1) :: bxs
  type (boxarray):: ba
  type (layout) :: la
  type (multifab) :: mf
  type (multifab) :: mfc
  type (boxassoc) :: bxasc
  integer :: i, j, n, nx, ny, dm
  logical, allocatable :: nodal(:)
  real(kind=dp_t), pointer :: fp(:,:,:,:)

  dm = 2
  allocate(nodal(dm)); nodal = .true.
  bx = refine(unit_box(dim=dm),4)
  bxs(1) = bx
!  bxs(2) = shift(bx,2,dim=1)
!  bxs(3) = shift(bx,2,dim=2)
!  bxs(4) = shift(bx,2)
  call build(ba, bxs)
  call boxarray_print(ba, "BOXARRAY")
  call build(la, ba)
  bxasc = layout_boxassoc(la, 1)
  call boxassoc_print(bxasc, "BOXASSOC")
  call multifab_build(mf, la, nc = 1, ng=1, nodal = nodal)
  call multifab_build(mfc, la, nc = 1, ng=1)
  call setval(mf, -1.0_dp_t, ALL=.True.)
  call setval(mfc, -1.0_dp_t, ALL=.True.)
  do n = 1, mf%nboxes; if ( remote(mf, n) ) cycle
     fp => dataptr(mf, n, get_ibox(mf, n))
     nx = size(fp,1)
     ny = size(fp,2)
     do j = 1, ny
        do i = 1, nx
           fp(i,j,1,1) = (i-1) !+ 10*(j-1)
        end do
     end do
     fp => dataptr(mfc, n, get_ibox(mf, n))
     nx = size(fp,1)
     ny = size(fp,2)
     do j = 1, ny
        do i = 1, nx
           fp(i,j,1,1) = (i-1) !+ 10*(j-1)
        end do
     end do
  end do
  call print(mf, "before")
  call print(mfc, "before")
  call fabio_multifab_write_d(mf, "tdir", "nd_flan")
  call fabio_multifab_write_d(mfc, "tdir", "cc_flan")
  call destroy(mf)
  call destroy(mfc)
  call destroy(la)
end subroutine t_nodal_mf_fabio

subroutine t_ba
  use boxarray_module
  use list_box_module
  use sort_box_module
  implicit none
  integer chunk(2), i
  type(boxarray) bxa

  call build(bxa, make_box((/2,8/), (/127,127/)))
  chunk = 64
  call print(bxa, "before maxsize")
  call boxarray_maxsize(bxa, chunk)
  call print(bxa, "after maxsize")
  print *, 'volume(bxa) = ', volume(bxa)
  do i = 1, nboxes(bxa)
     print *, 'vol', volume(get_box(bxa,i))
  end do
  call boxarray_simplify(bxa)
  print *, 'simplified volume(bxa) = ', volume(bxa)
  do i = 1, nboxes(bxa)
     print *, 'simplified vol', volume(get_box(bxa,i))
  end do
  call destroy(bxa)
end subroutine t_ba

subroutine t_cluster
  use cluster_module
  call cluster_set_verbose(.true.)
  call t_cls
end subroutine t_cluster

subroutine t_box_mod
  use box_module
  use boxarray_module
  implicit none
  type(box) :: bx, pd, bxi
  type(box), allocatable :: r(:,:)
  type(boxarray) :: ba
  integer, parameter :: dm = 2
  integer :: lo(dm), hi(dm), hi1(dm)
  integer :: i, n, m
  logical :: pmask(dm)

  allocate(r(3**dm,2))
  lo = 0
  hi = 3
  hi1 = 7
  pmask = .true.
  if ( dm > 1 ) pmask(2) = .false.

  call build(bxi, lo, hi)
  call build(pd, lo, hi1)

  call print(bxi, "bx")
  call print(pd, "pd")

  call boxarray_box_boundary_n(ba, bxi, 1)
  do m = 1, nboxes(ba)
     bx = get_box(ba, m)
     print *, 'm = ', m
     call print(bx, "bx")
     call box_decompose(r(:,1), n, bx, bxi)
     print *, 'BOX DECOMPOSE'
     do i = 1, n
        print *, 'i = ', i
        call print(r(i,1))
     end do
     
     call box_decompose_mod(r, n, bx, bxi, pd, pmask)
     
     print *, 'BOX DECOMPOSE'
     do i = 1, n
        print *, 'i = ', i
        call print(r(i,1), advance ='no'); write(*,fmt='("-->")', advance='no'); call print(r(i,2))
     end do
  end do
  call destroy(ba)

end subroutine t_box_mod


subroutine t_domain
  use f2kcli
  use bl_IO_module
  use mboxarray_module
  use box_util_module
  use mt19937_module
  implicit none
  type(mboxarray) :: mba
  character(len=128) :: dfile, test_set
  logical :: test_set_mglib
  logical :: test_set_hgproj
  integer :: narg, farg
  integer :: test_lev
  character(len=128) :: fname
  logical :: test_random_boxes
  integer :: random_min_size, random_max_size
  integer :: random_blocking_factor, random_num_boxes, random_iseed
  integer :: ba_maxsize, pd_xyz(MAX_SPACEDIM)
  integer :: dm
  integer :: i, un
  logical :: verbose
  logical :: lexist
  type(box) :: pd
  type(boxarray) :: ba, bado

  namelist /probin/ verbose
  namelist /probin/ dm
  namelist /probin/ test_random_boxes
  namelist /probin/ random_blocking_factor, random_min_size
  namelist /probin/ random_max_size, random_num_boxes, random_iseed
  namelist /probin/ dfile
  namelist /probin/ test_set
  namelist /probin/ test_set_mglib
  namelist /probin/ test_set_hgproj
  namelist /probin/ test_lev
  namelist /probin/ pd_xyz

  verbose = .true.
  dm = 2
  pd_xyz         = 32

  test_set_mglib = .FALSE.
  test_set_hgproj = .FALSE.
  test_lev    = 0
  test_set    = ''

  test_random_boxes = .false.
  random_blocking_factor = 8
  random_blocking_factor = 1
  random_min_size =  8
  random_max_size = 32
  random_num_boxes = 1
  random_iseed = 1

  narg = command_argument_count()
  farg = 1

  if ( narg >= 1 ) then
     call get_command_argument(farg, value = fname)
     inquire(file = fname, exist = lexist )
     if ( lexist ) then
        farg = farg + 1
        un = unit_new()
        open(unit=un, file = fname, status = 'old', action = 'read')
        read(unit=un, nml = probin)
        close(unit=un)
     end if
  end if

  do while ( farg <= narg )

     call get_command_argument(farg, value = fname)

     select case (fname)
     case ('--dim')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname, *) dm
     case ('--verbose')
        verbose = .true.

     case ('--test_set')
        farg = farg + 1
        call get_command_argument(farg, value = test_set)
     case ('--test_lev')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname, *) test_lev
     case ('--test_set_mglib')
        test_set_mglib = .True.
     case ('--test_set_hgproj')
        test_set_hgproj = .True.

     case ('--test_random_boxes')
        test_random_boxes = .True.
     case ('--random_blocking_factor')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname,*) random_blocking_factor
     case ('--random_min_size')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname,*) random_min_size
     case ('--random_max_size')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname,*) random_max_size
     case ('--random_num_boxes')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname,*) random_num_boxes
     case ('--random_iseed')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname,*) random_iseed

     case ('--pd_x')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname, *) pd_xyz(1)
     case ('--pd_y')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname, *) pd_xyz(2)
     case ('--pd_z')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname, *) pd_xyz(3)
     case ('--pd_xyz')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname, *) i
        pd_xyz = i

     case ('--')
        farg = farg + 1
        exit

     case default
        if ( .not. parallel_q() ) then
           write(*,*) 'UNKNOWN option = ', fname
           call bl_error("MAIN")
        end if
     end select

     farg = farg + 1
  end do

  if ( test_set_hgproj ) then
     call read_a_hgproj_grid(mba, test_set)
  else if ( test_set_mglib ) then
     call read_a_mglib_grid(mba, test_set)
  else if ( test_random_boxes ) then
     call box_build_2(pd, (/(0,i=1,dm)/), pd_xyz(1:dm)-1)
     call init_genrand(random_iseed)
     call build_random_boxarray(ba, pd, &
          random_num_boxes, random_min_size, &
          random_max_size, random_blocking_factor)
     call build(mba, ba, pd)
     call destroy(ba)
  else if ( test_set /= '' ) then
     call mboxarray_read_boxes(mba, test_set)
  else
     call box_build_2(pd, (/(0,i=1,dm)/), pd_xyz(1:dm)-1)
     call build(ba, pd)
     call build(mba, ba)
     call destroy(ba)
  end if

  if ( test_lev == 0 ) test_lev = mba%nlevel

  if ( .false. .and. verbose ) then
     call print(mba)
  end if

  call copy(ba, mba%bas(test_lev))
  pd = mba%pd(test_lev)


  call boxarray_simplify(ba)

  call boxarray_sort(ba)

! call print(ba,  "BA")

  if ( .not. boxarray_clean(ba%bxs) ) then
     print *, 'ba is not clean'
  end if

  call boxarray_decompose(bado, ba)

  call boxarray_sort(bado)

  if ( .not. boxarray_clean(bado%bxs) ) then
     print *, 'bado, decomposed is not clean'
  end if
! call print(bado, "BADO")

  print *, 'volume(bado) = ', volume(bado), ' nboxes = ', nboxes(bado)
  print *, 'volume(ba)   = ', volume(ba), ' nboxes = ', nboxes(ba)
  call boxarray_simplify(bado)
  print *, 'volume(bado) = ', volume(bado), ' nboxes = ', nboxes(bado)
  if ( .not. boxarray_clean(bado%bxs) ) then
     print *, 'bado is not clean'
  end if

  if ( .false. ) then
     call boxarray_pn_domain_bx(bado, ba, pd, 2)
  end if

! call print(bado, "BADO")

  call destroy(bado)
  call destroy(mba)
  call destroy(ba)

  call print(boxarray_mem_stats(),  " boxarray")
  call print(mboxarray_mem_stats(), "mboxarray")

end subroutine t_domain

subroutine t_box_chop
  use box_module
  implicit none
  type(box) :: bx, bxl, bxr

  bx = make_box((/0,0/), (/3,3/))
  call box_chop(bx, bxl, bxr, 1, 0)
  call print(bxl)
  call print(bxr)

  bx = make_box((/0,0/), (/3,3/))
  call box_chop(bx, bxl, bxr, 1, 3)
  call print(bxl)
  call print(bxr)
  
end subroutine t_box_chop

subroutine t_knap
  use knapsack_module
  implicit none
  call t_knapsack
end subroutine t_knap

subroutine t_timer
  use bl_timer_module
  implicit none
  real(kind=dp_t) d
  print *, 'MY_CPU_SECOND_TICK = ', MY_CPU_SECOND_TICK()
  print *, 'MY_WALL_SECOND_TICK = ', MY_WALL_SECOND_TICK()
  call cpu_second_tick(d)
  print *, 'CPU_SECOND_TICK = ', d
  call wall_second_tick(d)
  print *, 'WALL_SECOND_TICK = ', d
end subroutine t_timer

subroutine t_bl_types
  use bl_types
  use bl_IO_module
  implicit none
  call bl_types_info(unit_stdout())
end subroutine t_bl_types
