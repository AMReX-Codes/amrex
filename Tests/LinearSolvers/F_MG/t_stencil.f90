subroutine t_stencil
  use stencil_module
  use multifab_module
  use itsol_module
  use bl_IO_module
  use fabio_module
  use bl_constants_module
  use f2kcli
  implicit none
  type(box) :: pd
  type(boxarray):: ba
  type(layout) :: la
  type(stencil) :: st_2, st_4
  type(multifab) :: uu
  type(multifab) :: ff_2, ff_4
  type(multifab) :: tm, tms
  type(multifab) :: exact_rhs, exact_uu
  type(multifab) :: coeffs

  integer :: dm, ms, sz
  integer :: type_2, type_4
  integer :: exo_2,  exo_4
  real(dp_t), allocatable :: dh(:), xa(:), xb(:), hh(:)
  logical, allocatable :: pmask(:)
  real(dp_t) :: hhh
  logical, allocatable :: nodal(:)
  integer, allocatable :: bc_face(:,:)
  real(dp_t) :: alpha

  integer :: verbose
  integer :: max_iter
  real(dp_t) :: eps

  character(len=128) :: fname
  logical :: lexist, pd_pmask(MAX_SPACEDIM), tl
  integer :: un
  integer :: narg, farg

  namelist /psten/ eps
  namelist /psten/ max_iter
  namelist /psten/ verbose
  namelist /psten/ dm
  namelist /psten/ sz
  namelist /psten/ ms
  namelist /psten/ hhh
  namelist /psten/ pd_pmask

  verbose = 1
  dm = 2
  sz = 128
  ms = -Huge(ms)
  hhh = 1
  eps       = 1.0e-6_dp_t
  max_iter = -Huge(max_iter)
  pd_pmask = .false.

  narg = command_argument_count()
  farg = 1

  if ( narg >= 1 ) then
     call get_command_argument(1, value = fname)
     inquire(file = fname, exist = lexist )
     if ( lexist ) then
        un = unit_new()
        farg = farg + 1
        open(unit=un, file = fname, status = 'old', action = 'read')
        read(unit=un, nml = psten)
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
     case ('--size')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname, *) sz
     case ('--maxsize')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname, *) ms
     case ('--eps')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname, *) eps
     case ('--pd_pmask_x')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname, *) pd_pmask(1)
     case ('--pd_pmask_y')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname, *) pd_pmask(2)
     case ('--pd_pmask_z')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname, *) pd_pmask(3)
     case ('--pd_pmask')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname, *) tl
        pd_pmask = tl
     case ('--max_iter')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname, *) max_iter
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

  ms = max(ms, sz)
  if ( max_iter < 0 ) max_iter = sz**dm

  allocate(dh(dm), nodal(dm), xa(dm), xb(dm), hh(dm), pmask(dm))
  allocate(bc_face(dm,2))
  pmask = pd_pmask(1:dm)
  bc_face = BC_DIR
  hh = hhh
  alpha = hhh/20
  dh = hh/sz
  nodal = .false.
  xb = HALF*dh
  xa = HALF*dh

  pd = refine(unit_box(dim=dm),sz)
  call build(ba, pd)
  call boxarray_maxsize(ba, ms)
  if ( verbose > 2 ) then
     call boxarray_print(ba)
  end if
  call build(la, ba, pd = pd, pmask = pmask)

  call build(uu, la, 1, 1)
  call build(ff_2, la, 1, 0)
  call build(ff_4, la, 1, 0)
  call build(tm, la, 1, 1)
  call build(tms, la, 1, 0)
  call build(exact_rhs, la, 1, 0)
  call build(exact_uu, la, 1, 1)
  call build(coeffs, la, 1+dm, 1)
  call setval(coeffs, ZERO, 1, all=.true.)
  call setval(coeffs, ONE, 2, dm, all=.true.)


  !! -LAPLACIAN
  exo_2 = 2
  type_2 = ST_CROSS
  call stencil_build(st_2, la, dh, type_2)
  st_2%xa = xa
  st_2%xb = xb
  call stencil_set_bc_st(st_2, bc_face)
  call stencil_set_extrap_bc(st_2, exo_2)
  call stencil_fill(st_2, ba, bc_face, coeffs = coeffs, fill = ST_FILL_LAPLACE_2)

  !! -MEHERSTALLEN
  exo_4 = 5
  type_4 = ST_DENSE
  call stencil_build(st_4, la, dh, type_4)
  st_4%xa = xa
  st_4%xb = xb
  call stencil_set_bc_st(st_4, bc_face)
  call stencil_set_extrap_bc(st_4, exo_4)
  call stencil_fill(st_4, ba, bc_face, coeffs = coeffs, fill = ST_FILL_LAPLACE_4)

  if ( verbose > 2 ) then
     call stencil_print(st_2, str = "ST_2")
     call stencil_print(st_4, str = "ST_4")
  end if

  call fabio_write(st_2%ss, "tdir", "ss_2")
  call fabio_write(st_4%ss, "tdir", "ss_4")

  !! Stencil check ....
  call mf_init(exact_rhs, pd, .true.)
  call fabio_write(exact_rhs, "tdir", "exact_rhs")
  call mf_init(exact_uu,  pd, .false.)
  call fabio_write(exact_uu, "tdir", "exact_uu")

  call stencil_defect_st(st_2, tms, exact_rhs, exact_uu)
  print *, sz, 'norm_l2(res_2) = ', norm_l2(tms)*sqrt(product(dh))
  print *, sz, 'norm_inf(res_2) = ', norm_inf(tms)
  call fabio_write(tms, "tdir", "res_2")
  
  call mehr_rhs(st_2, exact_rhs)
  call fabio_write(exact_rhs, "tdir", "mehr_exact_rhs")
  call stencil_defect_st(st_4, tms, exact_rhs, exact_uu)
  print *, sz, 'norm_l2(res_4) = ', norm_l2(tms)*sqrt(product(dh))
  print *, sz, 'norm_inf(res_4) = ', norm_inf(tms)
  call fabio_write(tms, "tdir", "res_4")

  !! CG solver tests
  if ( .false. ) then
     call mf_init(ff_2, pd, .true.)
     call fabio_write(ff_2, "tdir", "ff_2")
     call setval(uu, ZERO, all = .true. )
     verbose   =  1
     call CG_Solve_st(st_2, uu, ff_2, eps, max_iter, verbose)
     call fabio_write(uu, "tdir", "uu_2")
     print *, sz, 'norm_inf(uu_2) = ', norm_inf(uu)
     print *, sz, 'norm_l2(uu_2) = ', norm_l2(uu)

     call mf_init(ff_4, pd, .true.)
     call mehr_rhs(st_2, ff_4)
     call fabio_write(ff_4, "tdir", "ff_4")
     call setval(uu, ZERO, all = .true. )
     verbose   =  1
     call CG_Solve_st(st_4, uu, ff_4, eps, max_iter, verbose)
     call fabio_write(uu, "tdir", "uu_4")
     print *, sz, 'norm_inf(uu_4) = ', norm_inf(uu)
     print *, sz, 'norm_l2(uu_4) = ', norm_l2(uu)
  end if

  call stencil_destroy(st_2)
  call stencil_destroy(st_4)

  call destroy(uu)
  call destroy(ff_2)
  call destroy(ff_4)
  call destroy(tm)
  call destroy(exact_rhs)
  call destroy(exact_uu)
  call destroy(coeffs)
  call destroy(la)
  call destroy(ba)

contains

  subroutine finish
    call print(multifab_mem_stats(),  " multifab at end")
    call print(imultifab_mem_stats(), "imultifab at end")
    call print(lmultifab_mem_stats(), "lmultifab at end")
    call print(fab_mem_stats(),       "      fab at end")
    call print(ifab_mem_stats(),      "     ifab at end")
    call print(lfab_mem_stats(),      "     lfab at end")
    call print(boxarray_mem_stats(),  " boxarray at end")
    call print(boxassoc_mem_stats(),  " boxassoc at end")
    call print(layout_mem_stats(),    "   layout at end")
    stop
  end subroutine finish

  subroutine mf_init(mf, pd, deriv)
    type(multifab), intent(inout) :: mf
    integer i, ng
    type(box) bx
    type(box), intent(in) :: pd
    real(kind=dp_t), pointer :: mp(:,:,:,:)
    logical :: deriv

    ng = mf%ng
    call setval(mf, 0.0_dp_t, all = .true.)
    do i = 1, mf%nboxes; if ( remote(mf,i) ) cycle
       mp => dataptr(mf, i, get_ibox(mf, i))
       bx = get_ibox(mf, i)
       select case ( mf%dim ) 
       case (2)
          call init_2d(mp(:,:,1,1), bx, pd, deriv)
       case (3)
          call init_3d(mp(:,:,:,1), bx, pd, deriv)
       end select
    end do
  end subroutine mf_init

  subroutine init_2d(po, bx, pd, deriv)
    type(box), intent(in) :: bx, pd
    real(dp_t), intent(inout) :: po(bx%lo(1):,bx%lo(2):)
    integer :: i, j
    logical :: deriv
    real(dp_t) :: xcen(2)
    xcen =      (upb(pd) + 1)*dh/2
    do j = bx%lo(2), bx%hi(2)
       do i = bx%lo(1), bx%hi(1)
          po(i,j) = gau_2d(i,j, xcen, alpha, deriv)
       end do
    end do
  end subroutine init_2d
  subroutine init_3d(po, bx, pd, deriv)
    type(box), intent(in) :: bx, pd
    real(dp_t), intent(inout) :: po(bx%lo(1):,bx%lo(2):,bx%lo(3):)
    integer i, j, k
    logical :: deriv
    real(dp_t) :: xcen(3)
    xcen =      (upb(pd) + 1)*dh/2
    do k = bx%lo(3), bx%hi(3)
       do j = bx%lo(2), bx%hi(2)
          do i = bx%lo(1), bx%hi(1)
             po(i,j,k) = gau_3d(i,j,k,xcen,alpha, deriv)
          end do
       end do
    end do
  end subroutine init_3d
  function gau_2d(i, j, xcen, alpha,deriv) result(r)
    real(dp_t) :: r, x
    integer, intent(in) :: i, j
    logical :: deriv
    real(dp_t) :: xcen(2)
    real(dp_t) :: alpha
    x =  ( &
         +(((i+HALF)*dh(1) - xcen(1))**2) &
         +(((j+HALF)*dh(2) - xcen(2))**2) &
           )/(alpha**2)
    if ( deriv ) then
      r = exp(-x*HALF)*(x - ( 1 + dm - 1))
    else
      r = alpha**2*exp(-x*HALF)
    endif
  end function gau_2d
  function gau_3d(i, j, k, xcen, alpha,deriv) result(r)
    real(dp_t) :: r, x
    integer, intent(in) :: i, j, k
    logical :: deriv
    real(dp_t) :: xcen(3)
    real(dp_t) :: alpha
    x =  ( &
         +(((i+HALF)*dh(1) - xcen(1))**2) &
         +(((j+HALF)*dh(2) - xcen(2))**2) &
         +(((k+HALF)*dh(3) - xcen(3))**2) &
           )/(alpha**2)
    if ( deriv) then
       r = exp(-x*HALF)*(x - ( 1 + dm - 1))
    else
       r = alpha**2 * exp(-x*HALF)
    end if
  end function gau_3d


  subroutine mehr_rhs(st_2, aa)
    type(stencil), intent(in) :: st_2
    type(multifab), intent(inout) :: aa
    call copy(tm, aa)
    call stencil_apply_st(st_2, aa, tm)
    call saxpy(tm, -st_2%dh(1)**2/12.0_dp_t, aa)
    call copy(aa, tm)
  end subroutine mehr_rhs

end subroutine t_stencil

subroutine t_stencil_map
  use stencil_module
  implicit none
  integer i, j, k
  print *, 'DENSE_MAP_3D'
  do k = -1, 1; do j = -1, 1; do i = -1, 1
     print *, i,j,k, ST_DENSE_MAP_3D(i,j,k)
  end do; end do; end do
  print *, 'DENSE_MAP_2D'
  do j = -1, 1; do i = -1, 1
     print *, i,j, ST_DENSE_MAP_2D(i,j)
  end do; end do
  print *, 'DENSE_MAP_1D'
  do i = -1, 1
     print *, i, ST_DENSE_MAP_1D(i)
  end do
end subroutine t_stencil_map
