module cpp_mg_module
  use mg_module
  use ml_layout_module
  use ml_multifab_module
  use bndry_reg_module

  implicit none

  type mg_server
     logical        :: final = .false.
     logical        :: nodal
     integer        :: dim  = 0
     integer        :: han  = -1
     integer        :: nlevel
     integer        :: stencil_order = 2
     type(ml_layout) :: mla
     type(mg_tower), pointer :: mgt(:) => Null()
     type(box), pointer :: pd(:) => Null()
     integer, pointer :: bc(:,:) => Null()
     integer, pointer :: rr(:,:)
     type(multifab), pointer :: rh(:) => Null()
     type(multifab), pointer :: uu(:) => Null()
     type(multifab), pointer :: coeffs(:) => Null()
     real(dp_t), pointer :: xa(:) => Null()
     real(dp_t), pointer :: xb(:) => Null()
     real(dp_t), pointer :: pxa(:) => Null()
     real(dp_t), pointer :: pxb(:) => Null()
  end type mg_server

  integer, parameter :: MGT_NU1 = 0
  integer, parameter :: MGT_NU2 = 1
  integer, parameter :: MGT_GAMMA = 2
  integer, parameter :: MGT_MAX_ITER = 3
  integer, parameter :: MGT_MAX_NLEVEL = 4
  integer, parameter :: MGT_MIN_WIDTH = 5
  integer, parameter :: MGT_VERBOSE = 6
  integer, parameter :: MGT_BOTTOM_SOLVER = 7
  integer, parameter :: MGT_BOTTOM_MAX_ITER = 8
  integer, parameter :: MGT_SMOOTHER = 9
  integer, parameter :: MGT_CYCLE = 10

  integer, parameter :: MGT_EPS = 0
  integer, parameter :: MGT_OMEGA = 1
  integer, parameter :: MGT_BOTTOM_SOLVER_EPS = 2

  integer, parameter :: max_mgt = 10
  type(mg_server), save :: mgts(max_mgt)

contains
  
  subroutine mgt_verify(mgt, str)
    integer, intent(in) :: mgt
    character(len=*), intent(in) :: str

    if ( mgt < 1 .or. mgt > max_mgt ) then
       call bl_error( trim(str) // ": MGT out of bounds: mgt: ", mgt)
    end if
    if ( mgts(mgt)%dim == 0 ) then
       call bl_error( trim(str) // ": MGT invalid DIM: not allocated: mgt: ", mgt)
    end if
    if ( mgts(mgt)%han /= mgt ) then
       call bl_error( trim(str) // ": MGT invalid han: mgt: ", mgt)
    end if
    
  end subroutine mgt_verify

  subroutine mgt_verify_lev(mgt, str, lev)
    integer, intent(in) :: mgt, lev
    character(len=*), intent(in) :: str
    call mgt_verify(mgt, str)
    if ( lev < 1 .or. lev > mgts(mgt)%nlevel ) then
       call bl_error( trim(str) // ": Level out of bounds", lev)
    end if
  end subroutine mgt_verify_lev

  subroutine mgt_verify_n(mgt, str, lev, n, lo, hi)
    integer, intent(in) :: mgt, lev, n, lo(:), hi(:)
    character(len=*), intent(in) :: str
    type(box) :: bx

    call mgt_verify_lev(mgt, str, lev)
    if ( n < 1 .or. n > nboxes(mgts(mgt)%mla, lev) ) then
       call bl_error( trim(str) // ": Box out of bounds", n)
    end if
    bx = make_box(lo, hi)
    if ( bx /= get_box(mgts(mgt)%mla, lev, n) ) then
       call bl_error( trim(str) // ": Box no filling")
    end if

  end subroutine mgt_verify_n

  subroutine mgt_not_final(mgt, str)
    integer, intent(in) :: mgt
    character(len=*), intent(in) :: str

    call mgt_verify(mgt, str)
    if ( mgts(mgt)%final ) then
       call bl_error( trim(str) // ": Changes made to finalized solver!")
    end if
  end subroutine mgt_not_final

end module cpp_mg_module

subroutine mgt_alloc(mgt, dm, nlevel, nodal)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: dm, nlevel
  integer, intent(out) :: mgt
  integer :: nodal
  integer i

  do i = 1, max_mgt
     if ( mgts(i)%dim == 0 ) then
        mgt = i
        mgts(i)%dim = dm
        mgts(i)%han =  i
        mgts(i)%nlevel = nlevel
     end if
  end do
  mgts(mgt)%nodal = nodal /= 0
  allocate(mgts(mgt)%rr(nlevel-1,dm))
  allocate(mgts(mgt)%rh(nlevel))
  allocate(mgts(mgt)%pd(nlevel))
  allocate(mgts(mgt)%uu(nlevel))
  allocate(mgts(mgt)%mgt(nlevel))
  call build(mgts(mgt)%mla, nlevel, dm)

end subroutine mgt_alloc

subroutine mgt_set_level(mgt, lev, nb, dm, lo, hi, pd_lo, pd_hi, bc, pm, pmap)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: mgt, lev, nb, dm
  integer, intent(in) :: lo(nb,dm), hi(nb,dm), pd_lo(dm), pd_hi(dm), bc(2,dm), pm(dm), pmap(nb+1)

  type(box) :: bxs(nb)
  integer   :: i
  integer   :: nc
  logical   :: pmask(dm)
  type(box) :: pd
  integer   :: flev
  logical, allocatable :: nodal(:)

  flev = lev + 1
  call mgt_verify_lev(mgt, "MGT_SET_LEVEL", flev)

  pmask = (pm /= 0)

  allocate(nodal(dm))

  if ( dm /= mgts(mgt)%dim ) then
     call bl_error("MGT_SET_LEVEL: Input DIM doesn't match internal DIM")
  end if
  call build(mgts(mgt)%pd(flev), pd_lo(1:dm), pd_hi(1:dm))
  do i = 1, nb
     bxs(i) = make_box(lo(i,:), hi(i,:))
  end do
  call build(mgts(mgt)%mla%mba%bas(flev), bxs)
  call build(mgts(mgt)%mla%la(flev),  &
       mgts(mgt)%mla%mba%bas(flev), &
       mgts(mgt)%pd(flev), pmask = pmask, &
       mapping = LA_EXPLICIT, explicit_mapping = pmap(1:nb))

  print *, 'flev = ', flev
  call print(mgts(mgt)%pd(flev))
  call print(mgts(mgt)%mla%la(flev))

  call print(mgts(mgt)%mla%la(flev))

  allocate(mgts(mgt)%bc(dm,2))

  mgts(mgt)%bc = transpose(bc)

end subroutine mgt_set_level

subroutine mgt_finalize(mgt, nipar, ipar, nrpar, rpar)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: mgt, nipar, nrpar, ipar(0:*)
  real(dp_t) :: rpar(0:*)
  integer :: dm, i, nlev, n
  integer :: ns
  integer :: nc
  logical, allocatable :: nodal(:)

  integer :: max_nlevel_in
  integer :: bottom_solver_in
  integer :: bottom_max_iter_in

  integer :: nu1
  integer :: nu2
  integer :: gamma
  integer :: max_iter
  integer :: max_nlevel
  integer :: min_width
  integer :: verbose
  integer :: bottom_solver
  integer :: bottom_max_iter
  integer :: smoother
  integer :: cycle

  real(dp_t) :: eps
  real(dp_t) :: omega
  real(dp_t) :: bottom_solver_eps
  type(boxarray) :: bac

  call mgt_verify(mgt, "MGT_FINALIZE")

  nu1 = ipar(MGT_NU1)
  nu2 = ipar(MGT_NU2)
  gamma = ipar(MGT_GAMMA)
  max_iter = ipar(MGT_MAX_ITER)
  max_nlevel = ipar(MGT_MAX_NLEVEL)
  min_width = ipar(MGT_MIN_WIDTH)
  verbose = ipar(MGT_VERBOSE)
  bottom_solver = ipar(MGT_BOTTOM_SOLVER)
  bottom_max_iter = ipar(MGT_BOTTOM_MAX_ITER)
  smoother = ipar(MGT_SMOOTHER)
  cycle = ipar(MGT_CYCLE)

  eps = rpar(MGT_EPS)
  omega = rpar(MGT_OMEGA)
  bottom_solver_eps = rpar(MGT_BOTTOM_SOLVER_EPS)

! print *, 'rpar ', rpar(0:2)

  dm = mgts(mgt)%dim

  nc = 1
  nlev = mgts(mgt)%nlevel

  do i = 1, nlev-1
     mgts(mgt)%rr(i,:) = mgts(mgt)%mla%mba%rr(i,:)
  end do

  do i = 1, nlev
     call build(mgts(mgt)%uu(i), mgts(mgt)%mla%la(i), nc, ng = 1)
     call build(mgts(mgt)%rh(i), mgts(mgt)%mla%la(i), nc, ng = 0)
  end do

  do i = nlev-1, 1, -1
     call build(mgts(mgt)%mla%mask(i), mgts(mgt)%mla%la(i), nc = 1, ng = 0)
     call setval(mgts(mgt)%mla%mask(n), val = .TRUE.)
     call copy(bac, mgts(mgt)%mla%mba%bas(n+1))
     call boxarray_coarsen(bac, mgts(mgt)%rr(n,:))
     call setval(mgts(mgt)%mla%mask(n), .false., bac)
     call destroy(bac)
  end do

  allocate(nodal(1:dm))
  nodal = mgts(mgt)%nodal
  ns = 1 + dm*3

  do n = nlev, 1, -1
     if ( n == 1 ) then
        max_nlevel_in = max_nlevel
        bottom_solver_in = bottom_solver
        bottom_max_iter_in = bottom_max_iter
     else
        if ( all(mgts(mgt)%rr == 2) ) then
           max_nlevel_in = 1
        else if ( all(mgts(mgt)%rr == 4) ) then
           max_nlevel_in = 2
        else
           call bl_error("MGT_FINALIZE: confused about ref_ratio")
        end if
        bottom_solver_in = 0
        bottom_max_iter_in = nu1
     end if
     call mg_tower_build(mgts(mgt)%mgt(n), mgts(mgt)%mla%la(n), mgts(mgt)%pd(n), mgts(mgt)%bc, &
          ns = ns, &
          smoother = smoother, &
          nu1 = nu1, nu2 = nu2, gamma = gamma, cycle = cycle, omega = omega, &
          bottom_solver = bottom_solver_in, &
          bottom_max_iter = bottom_max_iter_in, &
          bottom_solver_eps = bottom_solver_eps, &
          max_iter = max_iter, &
          max_nlevel = max_nlevel_in, &
          min_width = min_width, &
          eps = eps, &
          verbose = verbose, &
          nodal = nodal &
          )
  end do

end subroutine mgt_finalize

subroutine mgt_init_coeffs_lev(mgt, lev)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: mgt, lev
  integer :: nlev, dm, i
  integer :: flev
  flev = lev + 1
  call mgt_verify_lev(mgt, "MGT_INIT_STENCIL_LEV", flev)

  dm = mgts(mgt)%dim
  nlev = mgts(mgt)%mgt(flev)%nlevels
  allocate(mgts(mgt)%coeffs(nlev))
print *, 'size(coeffs) ', size(mgts(mgt)%coeffs)
call print(mgts(mgt)%mla%la(1))

  call build(mgts(mgt)%coeffs(nlev), mgts(mgt)%mgt(flev)%ss(nlev)%la, 1+dm, 1)
  do i = nlev - 1, 1, -1
     call build(mgts(mgt)%coeffs(i), mgts(mgt)%mgt(flev)%ss(i)%la, 1+dm, 1)
     call setval(mgts(mgt)%coeffs(i), 0.0_dp_t, all=.true.)
  end do

  do i = 1, nlev
     print *, 'i ', i, nboxes(mgts(mgt)%coeffs(i))
  end do

end subroutine mgt_init_coeffs_lev

subroutine mgt_set_coeffs_lev(mgt, lev, n)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: mgt, lev, n
  integer :: nlev
  integer :: flev
  flev = lev + 1
  call mgt_verify_lev(mgt, "MGT_SET_COEFS_LEV", flev)
end subroutine mgt_set_coeffs_lev

subroutine mgt_finalize_stencil_lev(mgt, lev, xa, xb, pxa, pxb)
  use cpp_mg_module
  use coeffs_module
  implicit none
  integer, intent(in) :: mgt, lev
  real(dp_t), intent(in) :: xa(*), xb(*), pxa(*), pxb(*)
  integer :: nlev, i, dm
  integer :: flev
  type(box) :: pd
  type(boxarray) :: pdv
  dm = mgts(mgt)%dim
  flev = lev + 1
  call mgt_verify_lev(mgt, "MGT_SET_COEFS_LEV", flev)
  
  pd = mgts(mgt)%mla%mba%pd(flev)
  nlev = mgts(mgt)%mgt(flev)%nlevels
  do i = nlev-1, 1, -1
     call coarsen_coeffs(mgts(mgt)%coeffs(i+1), mgts(mgt)%coeffs(i))
  end do
  do i = nlev, 1, -1
     pdv = layout_boxarray(mgts(mgt)%mgt(flev)%ss(i)%la)
print *, nboxes(pdv)
print *, nboxes(mgts(mgt)%coeffs(i))
     call stencil_fill_cc(mgts(mgt)%mgt(flev)%ss(i), mgts(mgt)%coeffs(i), mgts(mgt)%mgt(flev)%dh(:,i), &
          pdv, mgts(mgt)%mgt(flev)%mm(i), xa(1:dm), xb(1:dm), pxa(1:dm), pxb(1:dm), pd, &
          mgts(mgt)%stencil_order, mgts(mgt)%bc)
     call destroy(mgts(mgt)%coeffs(i))
  end do
  deallocate(mgts(mgt)%coeffs)

end subroutine mgt_finalize_stencil_lev

subroutine mgt_finalize_stencil(mgt)
   use cpp_mg_module
   implicit none
   integer, intent(in) :: mgt
   call mgt_verify(mgt, "MGT_FINALIZE_STENCIL")
   mgts(mgt)%final = .true.
end subroutine mgt_finalize_stencil


subroutine mgt_set_ss_2d(mgt, lev, n, a, b, aa, bb1, bb2, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: mgt, lev, n, lo(2), hi(2), plo(2), phi(2)
  real(kind=dp_t), intent(in) :: a, b
  real(kind=dp_t), intent(in) :: aa(plo(1):phi(1), plo(2):phi(2))
  real(kind=dp_t), intent(in) :: bb1(plo(1):phi(1), plo(2):phi(2))
  real(kind=dp_t), intent(in) :: bb2(plo(1):phi(1), plo(2):phi(2))
end subroutine mgt_set_ss_2d
subroutine mgt_set_ss_3d(mgt, lev, n, a, b, aaa, bb1, bb2, bb3, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: mgt, lev, n, lo(2), hi(2), plo(2), phi(2)
  real(kind=dp_t), intent(in) :: a, b
  real(kind=dp_t), intent(in) :: aaa(plo(1):phi(1), plo(2):phi(2))
  real(kind=dp_t), intent(in) :: bb1(plo(1):phi(1), plo(2):phi(2))
  real(kind=dp_t), intent(in) :: bb2(plo(1):phi(1), plo(2):phi(2))
  real(kind=dp_t), intent(in) :: bb3(plo(1):phi(1), plo(2):phi(2))
end subroutine mgt_set_ss_3d

subroutine mgt_set_rh_1d(mgt, lev, n, rh, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: mgt, lev, n, lo(1), hi(1), plo(1), phi(1)
  real(kind=dp_t), intent(in) :: rh(plo(1):phi(1))
  real(kind=dp_t), pointer :: rp(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1
  call mgt_verify_n(mgt, "MGT_SET_RH", flev, fn, lo, hi)

  rp => dataptr(mgts(mgt)%rh(flev), fn)
  rp(lo(1):hi(1), 1,1,1) = rh(lo(1):hi(1))

end subroutine mgt_set_rh_1d
subroutine mgt_set_rh_2d(mgt, lev, n, rh, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: mgt, lev, n, lo(2), hi(2), plo(2), phi(2)
  real(kind=dp_t), intent(in) :: rh(plo(1):phi(1), plo(2):phi(2))
  real(kind=dp_t), pointer :: rp(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1
  
  call mgt_verify_n(mgt, "MGT_SET_RH", flev, fn, lo, hi)

  rp => dataptr(mgts(mgt)%rh(flev), fn)
  rp(lo(1):hi(1), lo(2):hi(2),1,1) = rh(lo(1):hi(1), lo(2):hi(2))

end subroutine mgt_set_rh_2d
subroutine mgt_set_rh_3d(mgt, lev, n, rh, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: mgt, lev, n, lo(3), hi(3), plo(3), phi(3)
  real(kind=dp_t), intent(in) :: rh(plo(1):phi(1), plo(2):phi(2), plo(3):phi(3))
  real(kind=dp_t), pointer :: rp(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1
  
  call mgt_verify_n(mgt, "MGT_SET_RH", flev, fn, lo, hi)

  rp => dataptr(mgts(mgt)%rh(lev), fn)
  rp(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3),1) = rh(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))

end subroutine mgt_set_rh_3d

subroutine mgt_get_rh_1d(mgt, lev, n, rh, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: mgt, lev, n, lo(1), hi(1), plo(1), phi(1)
  real(kind=dp_t), intent(inout) :: rh(plo(1):phi(1))
  real(kind=dp_t), pointer :: rp(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1
  
  call mgt_verify_n(mgt, "MGT_GET_RH", flev, fn, lo, hi)

  rp => dataptr(mgts(mgt)%rh(flev), fn)
  rh(lo(1):hi(1)) = rp(lo(1):hi(1), 1,1,1)

end subroutine mgt_get_rh_1d
subroutine mgt_get_rh_2d(mgt, lev, n, rh, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: mgt, lev, n, lo(2), hi(2), plo(2), phi(2)
  real(kind=dp_t), intent(inout) :: rh(plo(1):phi(1), plo(2):phi(2))
  real(kind=dp_t), pointer :: rp(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1
  
  call mgt_verify_n(mgt, "MGT_GET_RH", flev, fn, lo, hi)

  rp => dataptr(mgts(mgt)%rh(flev), fn)
  rh(lo(1):hi(1), lo(2):hi(2)) = rp(lo(1):hi(1), lo(2):hi(2),1,1)

end subroutine mgt_get_rh_2d
subroutine mgt_get_rh_3d(mgt, lev, n, rh, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: mgt, lev, n, lo(3), hi(3), plo(3), phi(3)
  real(kind=dp_t), intent(inout) :: rh(plo(1):phi(1), plo(2):phi(2), plo(3):phi(3))
  real(kind=dp_t), pointer :: rp(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1
  
  call mgt_verify_n(mgt, "MGT_GET_RH", flev, fn, lo, hi)

  rp => dataptr(mgts(mgt)%rh(flev), fn)
  rh(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)) = rp(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 1)

end subroutine mgt_get_rh_3d

subroutine mgt_set_cf_1d(mgt, lev, n, cf, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: mgt, lev, n, lo(1), hi(1), plo(1), phi(1)
  real(kind=dp_t), intent(in) :: cf(plo(1):phi(1),*)
  real(kind=dp_t), pointer :: cp(:,:,:,:)
  integer :: flev, fn, nlev
  fn = n + 1
  flev = lev+1
  nlev = size(mgts(mgt)%coeffs)
  call mgt_verify_n(mgt, "MGT_SET_UU", lev, n, lo, hi)

  cp => dataptr(mgts(mgt)%coeffs(nlev), fn)
  cp(lo(1):hi(1), 1,1, 1:4) = cf(lo(1):hi(1), 1:4)

end subroutine mgt_set_cf_1d

subroutine mgt_set_uu_1d(mgt, lev, n, uu, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: mgt, lev, n, lo(1), hi(1), plo(1), phi(1)
  real(kind=dp_t), intent(in) :: uu(plo(1):phi(1))
  real(kind=dp_t), pointer :: up(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1
  
  call mgt_verify_n(mgt, "MGT_SET_UU", lev, n, lo, hi)

  up => dataptr(mgts(mgt)%uu(lev), fn)
  up(lo(1):hi(1), 1,1,1) = uu(lo(1):hi(1))

end subroutine mgt_set_uu_1d

subroutine mgt_set_cf_2d(mgt, lev, n, cf, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: mgt, lev, n, lo(2), hi(2), plo(2), phi(2)
  real(kind=dp_t), intent(in) :: cf(plo(1):phi(1), plo(2):phi(2), *)
  real(kind=dp_t), pointer :: cp(:,:,:,:)
  integer :: flev, fn, nlev
  fn = n + 1
  flev = lev+1
  nlev = size(mgts(mgt)%coeffs)
  call mgt_verify_n(mgt, "MGT_SET_UU", flev, fn, lo, hi)

  cp => dataptr(mgts(mgt)%coeffs(nlev), fn)
  cp(lo(1):hi(1), lo(2):hi(2),1, 1:4) = cf(lo(1):hi(1), lo(2):hi(2), 1:4)

end subroutine mgt_set_cf_2d

subroutine mgt_set_uu_2d(mgt, lev, n, uu, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: mgt, lev, n, lo(2), hi(2), plo(2), phi(2)
  real(kind=dp_t), intent(in) :: uu(plo(1):phi(1), plo(2):phi(2))
  real(kind=dp_t), pointer :: up(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1
  
  call mgt_verify_n(mgt, "MGT_SET_UU", flev, fn, lo, hi)

  up => dataptr(mgts(mgt)%uu(flev), fn)
  up(lo(1):hi(1), lo(2):hi(2),1,1) = uu(lo(1):hi(1), lo(2):hi(2))

end subroutine mgt_set_uu_2d

subroutine mgt_set_cf_3d(mgt, lev, n, cf, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: mgt, lev, n, lo(3), hi(3), plo(3), phi(3)
  real(kind=dp_t), intent(in) :: cf(plo(1):phi(1), plo(2):phi(2), plo(3):phi(3), 4)
  real(kind=dp_t), pointer :: cp(:,:,:,:)
  integer :: flev, fn, nlev
  fn = n + 1
  flev = lev+1
  nlev = size(mgts(mgt)%coeffs)
  call mgt_verify_n(mgt, "MGT_SET_UU", flev, fn, lo, hi)

  cp => dataptr(mgts(mgt)%coeffs(nlev), fn)
  cp(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 1:4) = cf(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 1:4)

end subroutine mgt_set_cf_3d

subroutine mgt_set_uu_3d(mgt, lev, n, uu, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: mgt, lev, n, lo(3), hi(3), plo(3), phi(3)
  real(kind=dp_t), intent(in) :: uu(plo(1):phi(1), plo(2):phi(2), plo(3):phi(3))
  real(kind=dp_t), pointer :: up(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1
  
  call mgt_verify_n(mgt, "MGT_SET_UU", flev, fn, lo, hi)

  up => dataptr(mgts(mgt)%uu(flev), fn)
  up(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 1) = uu(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))

end subroutine mgt_set_uu_3d

subroutine mgt_get_uu_1d(mgt, lev, n, uu, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: mgt, lev, n, lo(1), hi(1), plo(1), phi(1)
  real(kind=dp_t), intent(inout) :: uu(plo(1):phi(1))
  real(kind=dp_t), pointer :: up(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1
  
  call mgt_verify_n(mgt, "MGT_GET_UU", flev, fn, lo, hi)

  up => dataptr(mgts(mgt)%uu(flev), fn)
  uu(lo(1):hi(1)) = up(lo(1):hi(1), 1,1,1)

end subroutine mgt_get_uu_1d
subroutine mgt_get_uu_2d(mgt, lev, n, uu, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: mgt, lev, n, lo(2), hi(2), plo(2), phi(2)
  real(kind=dp_t), intent(inout) :: uu(plo(1):phi(1), plo(2):phi(2))
  real(kind=dp_t), pointer :: up(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1
  
  call mgt_verify_n(mgt, "MGT_GET_UU", flev, fn, lo, hi)

  up => dataptr(mgts(mgt)%uu(flev), fn)
  uu(lo(1):hi(1), lo(2):hi(2)) = up(lo(1):hi(1), lo(2):hi(2),1,1)

end subroutine mgt_get_uu_2d
subroutine mgt_get_uu_3d(mgt, lev, n, uu, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: mgt, lev, n, lo(3), hi(3), plo(3), phi(3)
  real(kind=dp_t), intent(inout) :: uu(plo(1):phi(1), plo(2):phi(2), plo(3):phi(3))
  real(kind=dp_t), pointer :: up(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1
  
  call mgt_verify_n(mgt, "MGT_GET_UU", flev, fn, lo, hi)

  up => dataptr(mgts(mgt)%uu(flev), fn)
  uu(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)) = up(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 1)

end subroutine mgt_get_uu_3d

subroutine mgt_dealloc(mgt)
  use cpp_mg_module
  implicit none
  integer, intent(inout) :: mgt
  integer :: i
  
  call mgt_verify(mgt, "MGT_DEALLOC")
  if ( .not. mgts(mgt)%final ) then
     call bl_error("MGT_DEALLOC: MGT not finalized")
  end if

  call mg_tower_destroy(mgts(mgt)%mgt(1))
  do i = mgts(mgt)%nlevel, 1, -1
     call destroy(mgts(mgt)%rh(i))
     call destroy(mgts(mgt)%uu(i))
  end do
  call destroy(mgts(mgt)%mla)
  mgts(mgt)%han = -1
  mgts(mgt)%dim = 0
  mgt = -1

end subroutine mgt_dealloc

subroutine mgt_solve_cc(mgt)
  use cpp_mg_module
  use ml_cc_module
  use ml_nd_module
  use fabio_module
  implicit none
  integer, intent(in) :: mgt
  integer :: do_diagnostics
  real(dp_t) :: eps

  call mgt_verify(mgt, "MGT_SOLVE")
  if ( .not. mgts(mgt)%final ) then
     call bl_error("MGT_SOLVE: MGT not finalized")
  end if

call fabio_ml_write(mgts(mgt)%rh, mgts(mgt)%rr(:,1), "mgt_rhs")

  if ( mgts(mgt)%nodal ) then
     call ml_nd(mgts(mgt)%mla, mgts(mgt)%mgt, &
          mgts(mgt)%rh, mgts(mgt)%uu, &
          mgts(mgt)%mla%mask, mgts(mgt)%rr, &
          do_diagnostics, eps)
  else
     call ml_cc(mgts(mgt)%mla, mgts(mgt)%mgt, &
          mgts(mgt)%rh, mgts(mgt)%uu, &
          mgts(mgt)%mla%mask, mgts(mgt)%rr, &
          do_diagnostics, eps)
  end if

call fabio_ml_write(mgts(mgt)%uu, mgts(mgt)%rr(:,1), "mgt_uu")

end subroutine mgt_solve_cc

subroutine mgt_set_nu1(mgt, nu1)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: mgt, nu1

  call mgt_not_final(mgt, "MGT_SET_NU1")
  mgts(mgt)%mgt%nu1 = nu1

end subroutine mgt_set_nu1

subroutine mgt_set_nu2(mgt, nu2)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: mgt, nu2

  call mgt_not_final(mgt, "MGT_SET_NU2")
  mgts(mgt)%mgt%nu2 = nu2

end subroutine mgt_set_nu2

subroutine mgt_set_eps(mgt, eps)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: mgt
  real(kind=dp_t), intent(in) :: eps

  call mgt_not_final(mgt, "MGT_SET_EPS")
  mgts(mgt)%mgt%eps = eps

end subroutine mgt_set_eps

subroutine mgt_set_bottom_solver(mgt, bottom_solver)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: mgt
  integer, intent(in) :: bottom_solver

  call mgt_not_final(mgt, "MGT_SET_BOTTOM_SOLVER") 
  mgts(mgt)%mgt%bottom_solver = bottom_solver

end subroutine mgt_set_bottom_solver

subroutine mgt_set_bottom_max_iter(mgt, bottom_max_iter)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: mgt
  integer, intent(in) :: bottom_max_iter

  call mgt_not_final(mgt, "MGT_SET_BOTTOM_MAX_ITER")
  mgts(mgt)%mgt%bottom_max_iter = bottom_max_iter

end subroutine mgt_set_bottom_max_iter

subroutine mgt_set_bottom_solver_eps(mgt, bottom_solver_eps)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: mgt
  real(kind=dp_t), intent(in) :: bottom_solver_eps

  call mgt_not_final(mgt, "MGT_SET_BOTTOM_SOLVER_EPS")
  mgts(mgt)%mgt%bottom_solver_eps = bottom_solver_eps

end subroutine mgt_set_bottom_solver_eps

subroutine mgt_set_max_nlevel(mgt, max_nlevel)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: mgt
  integer, intent(in) :: max_nlevel

  call mgt_not_final(mgt, "MGT_SET_MAX_NLEVEL")
  mgts(mgt)%mgt%max_nlevel = max_nlevel

end subroutine mgt_set_max_nlevel

subroutine mgt_set_min_width(mgt, min_width)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: mgt
  integer, intent(in) :: min_width

  call mgt_not_final(mgt, "MGT_SET_MIN_WIDTH")
  mgts(mgt)%mgt%min_width = min_width

end subroutine mgt_set_min_width

subroutine mgt_set_verbose(mgt, verbose)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: mgt
  integer, intent(in) :: verbose

  call mgt_not_final(mgt, "MGT_SET_MIN_WIDTH")
  mgts(mgt)%mgt%verbose = verbose

end subroutine mgt_set_verbose

subroutine mgt_get_rpar_defaults(rpar, n)
  Use cpp_mg_module
  implicit none
  integer :: n
  real(dp_t), intent(out) :: rpar(0:n-1)
  type(mg_tower) :: mgt
  rpar(MGT_EPS) = mgt%eps
  rpar(MGT_OMEGA) = mgt%eps
  rpar(MGT_BOTTOM_SOLVER_EPS) = mgt%bottom_solver_eps
end subroutine mgt_get_rpar_defaults
subroutine mgt_get_ipar_defaults(ipar, n)
  Use cpp_mg_module
  implicit none
  integer, intent(in) :: n
  integer, intent(out) :: ipar(0:n-1)
  type(mg_tower) :: mgt
  if ( n < MGT_CYCLE-1 ) then
     call bl_Error("ipar array too small: ", n);
  end if
  ipar(MGT_NU1) = mgt%nu1
  ipar(MGT_NU2) = mgt%nu2
  ipar(MGT_GAMMA) = mgt%gamma
  ipar(MGT_MAX_ITER) = mgt%max_iter
  ipar(MGT_MAX_NLEVEL) = mgt%max_nlevel
  ipar(MGT_MIN_WIDTH) = mgt%min_width
  ipar(MGT_VERBOSE) = mgt%verbose
  ipar(MGT_BOTTOM_SOLVER) = mgt%bottom_solver
  ipar(MGT_BOTTOM_MAX_ITER) = mgt%bottom_max_iter
  ipar(MGT_SMOOTHER) = mgt%smoother
  ipar(MGT_CYCLE) = mgt%cycle
end subroutine mgt_get_ipar_defaults
