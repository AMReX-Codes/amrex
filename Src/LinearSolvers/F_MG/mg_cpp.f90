module cpp_mg_module
  use mg_module
  use ml_layout_module
  use ml_multifab_module
  use bndry_reg_module

  implicit none

  type mg_server
     logical        :: final = .false.
     integer        :: dim  = 0
     integer        :: han  = -1
     type(ml_boxarray) :: mba
     type(ml_layout)   :: mla
     type(mg_tower), pointer :: mgt(:)
     type(box)      :: pd
     integer, pointer :: bc(:,:) => Null()
     type(ml_multifab) :: rh
     type(ml_multifab) :: uu
     type(multifab), pointer :: coeffs(:) => Null()
     real(dp_t), pointer :: xa(:) => Null()
     real(dp_t), pointer :: xb(:) => Null()
     real(dp_t), pointer :: pxa(:) => Null()
     real(dp_t), pointer :: pxb(:) => Null()
     type(bndry_reg), pointer :: flx(:) => Null()
     type(bndry_reg), pointer :: bcs(:) => Null()
  end type mg_server

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

  subroutine mgt_verify_n(mgt, str, lev, n, lo, hi)
    integer, intent(in) :: mgt, lev, n, lo(:), hi(:)
    character(len=*), intent(in) :: str
    type(box) :: bx

    call mgt_verify(mgt, str)
    if ( n < 1 .or. n > nboxes(mgts(mgt)%mla%la(lev)) ) then
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

subroutine mgt_alloc(mgt, dm, nlevel)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: dm, nlevel
  integer, intent(out) :: mgt
  integer i

  do i = 1, max_mgt
     if ( mgts(i)%dim == 0 ) then
        mgt = i
        mgts(i)%dim = dm
        mgts(i)%han =  i
     end if
  end do
  call build(mgts(mgt)%mba, nlevel, dm)

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

  integer   :: flev

  call mgt_verify(mgt, "MGT_SET_LEVEL")

  flev = lev + 1

  pmask = (pm /= 0)

  if ( dm /= mgts(mgt)%dim ) then
     call bl_error("MGT_SET_LEVEL: Input DIM doesn't match internal DIM")
  end if
  call build(mgts(mgt)%pd, pd_lo(1:dm), pd_hi(1:dm))
  do i = 1, nb
     bxs(i) = make_box(lo(i,:), hi(i,:))
  end do
  call build(mgts(mgt)%mba%bas(flev), bxs)
  allocate(mgts(mgt)%bc(dm,2))

  print *, shape(bc)
  print *, shape(transpose(bc))
  print *, shape(mgts(mgt)%bc)

  do i = 1, dm
     print *, 'i', i, bc(:,i)
  end do

2  mgts(mgt)%bc = transpose(bc)

end subroutine mgt_set_level

subroutine mgt_finalize(mgt)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: mgt
  type(ml_layout) :: mla
  integer :: dm, i, nlev
  integer :: ns
  integer :: nc

  call mgt_verify(mgt, "MGT_FINALIZE")
  if ( .not. built_q(mgts(mgt)%mla) ) then
     call bl_error("MGT_FINALIZE: ml_layout not built")
  end if

  dm = mgts(mgt)%dim
  mla = mgts(mgt)%mla

!  call build(mgts(mgt)%mla%la(1), ba, pd = mgts(mgt)%pd, pmask = pmask, &
!       mapping = LA_EXPLICIT, explicit_mapping = pmap(1:nb))

  nc = 1

  call build(mgts(mgt)%uu, mgts(mgt)%mla, nc, ng = 1)
  call build(mgts(mgt)%rh, mgts(mgt)%mla, nc, ng = 0)

  ns = 1 + dm*3
  call mg_tower_build(mgts(mgt)%mgt(1), mla%la(1), mgts(mgt)%pd, mgts(mgt)%bc, &
       ns = ns &
       )

  nlev = nlevels(mgts(mgt)%mla)

  allocate(mgts(mgt)%coeffs(nlev))
  call build(mgts(mgt)%coeffs(nlev), mla%la(1), 1+dm, 1)
  do i = nlev - 1, 1, -1
     call multifab_build(mgts(mgt)%coeffs(i), mgts(mgt)%mgt(1)%ss(i)%la, 1+dm, 1)
     call setval(mgts(mgt)%coeffs(i), 0.0_dp_t, all=.true.)
  end do

  mgts(mgt)%final = .true.

end subroutine mgt_finalize

subroutine mgt_finalize_stencil(mgt, xa, xb, pxa, pxb)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: mgt
  real(dp_t), intent(in) :: xa(*), xb(*), pxa(*), pxb(*)
  integer :: nlev, i
  type(boxarray) :: pdv
  integer :: dm
  type(box) :: pd
  integer :: stencil_order

  call mgt_verify(mgt, "MGT_SOLVE")
  if ( .not. mgts(mgt)%final ) then
     call bl_error("MGT_SOLVE: MGT not finalized")
  end if

  dm = mgts(mgt)%dim
  allocate(mgts(mgt)%xa(dm), mgts(mgt)%xb(dm))
  allocate(mgts(mgt)%pxa(dm), mgts(mgt)%pxb(dm))

  ! May be un-needed but salt away for now.
  mgts(mgt)%pxa = pxa(1:dm)
  mgts(mgt)%pxb = pxb(1:dm)
  mgts(mgt)%xa  = xa(1:dm)
  mgts(mgt)%xb  = xb(1:dm)

print *, 'pxa', mgts(mgt)%pxa
print *, 'pxb', mgts(mgt)%pxb
print *, 'xa ', mgts(mgt)%xa
print *, 'xb ', mgts(mgt)%xb

  nlev = nlevels(mgts(mgt)%mla)
  pd = get_pd(mgts(mgt)%mla%la(1))

  do i = nlev - 1, 1, -1
     pdv = get_boxarray(mgts(mgt)%mgt(i)%ss(i)%la)
     call stencil_fill_cc(mgts(mgt)%mgt(1)%ss(i), mgts(mgt)%coeffs(i), mgts(mgt)%mgt(1)%dh(:,i), &
          pdv, mgts(mgt)%mgt(1)%mm(i), xa(1:dm), xb(1:dm), pxa(1:dm), pxb(1:dm), pd, &
          stencil_order, mgts(mgt)%bc)
  end do

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
  
  call mgt_verify_n(mgt, "MGT_SET_RH", lev, n, lo, hi)

  rp => dataptr(mgts(mgt)%rh, lev, n)
  rp(lo(1):hi(1), 1,1,1) = rh(lo(1):hi(1))

end subroutine mgt_set_rh_1d
subroutine mgt_set_rh_2d(mgt, lev, n, rh, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: mgt, lev, n, lo(2), hi(2), plo(2), phi(2)
  real(kind=dp_t), intent(in) :: rh(plo(1):phi(1), plo(2):phi(2))
  real(kind=dp_t), pointer :: rp(:,:,:,:)
  
  call mgt_verify_n(mgt, "MGT_SET_RH", lev, n, lo, hi)

  rp => dataptr(mgts(mgt)%rh, lev, n)
  rp(lo(1):hi(1), lo(2):hi(2),1,1) = rh(lo(1):hi(1), lo(2):hi(2))

end subroutine mgt_set_rh_2d
subroutine mgt_set_rh_3d(mgt, lev, n, rh, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: mgt, lev, n, lo(3), hi(3), plo(3), phi(3)
  real(kind=dp_t), intent(in) :: rh(plo(1):phi(1), plo(2):phi(2), plo(3):phi(3))
  real(kind=dp_t), pointer :: rp(:,:,:,:)
  
  call mgt_verify_n(mgt, "MGT_SET_RH", lev, n, lo, hi)

  rp => dataptr(mgts(mgt)%rh, lev, n)
  rp(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3),1) = rh(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))

end subroutine mgt_set_rh_3d

subroutine mgt_get_rh_1d(mgt, lev, n, rh, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: mgt, lev, n, lo(1), hi(1), plo(1), phi(1)
  real(kind=dp_t), intent(inout) :: rh(plo(1):phi(1))
  real(kind=dp_t), pointer :: rp(:,:,:,:)
  
  call mgt_verify_n(mgt, "MGT_GET_RH", lev, n, lo, hi)

  rp => dataptr(mgts(mgt)%rh, lev, n)
  rh(lo(1):hi(1)) = rp(lo(1):hi(1), 1,1,1)

end subroutine mgt_get_rh_1d
subroutine mgt_get_rh_2d(mgt, lev, n, rh, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: mgt, lev, n, lo(2), hi(2), plo(2), phi(2)
  real(kind=dp_t), intent(inout) :: rh(plo(1):phi(1), plo(2):phi(2))
  real(kind=dp_t), pointer :: rp(:,:,:,:)
  
  call mgt_verify_n(mgt, "MGT_GET_RH", lev, n, lo, hi)

  rp => dataptr(mgts(mgt)%rh, lev, n)
  rh(lo(1):hi(1), lo(2):hi(2)) = rp(lo(1):hi(1), lo(2):hi(2),1,1)

end subroutine mgt_get_rh_2d
subroutine mgt_get_rh_3d(mgt, lev, n, rh, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: mgt, lev, n, lo(3), hi(3), plo(3), phi(3)
  real(kind=dp_t), intent(inout) :: rh(plo(1):phi(1), plo(2):phi(2), plo(3):phi(3))
  real(kind=dp_t), pointer :: rp(:,:,:,:)
  
  call mgt_verify_n(mgt, "MGT_GET_RH", lev, n, lo, hi)

  rp => dataptr(mgts(mgt)%rh, lev, n)
  rh(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)) = rp(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 1)

end subroutine mgt_get_rh_3d

subroutine mgt_set_uu_1d(mgt, lev, n, uu, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: mgt, lev, n, lo(1), hi(1), plo(1), phi(1)
  real(kind=dp_t), intent(in) :: uu(plo(1):phi(1))
  real(kind=dp_t), pointer :: up(:,:,:,:)
  
  call mgt_verify_n(mgt, "MGT_SET_UU", lev, n, lo, hi)

  up => dataptr(mgts(mgt)%uu, lev, n)
  up(lo(1):hi(1), 1,1,1) = uu(lo(1):hi(1))

end subroutine mgt_set_uu_1d
subroutine mgt_set_uu_2d(mgt, lev, n, uu, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: mgt, lev, n, lo(2), hi(2), plo(2), phi(2)
  real(kind=dp_t), intent(in) :: uu(plo(1):phi(1), plo(2):phi(2))
  real(kind=dp_t), pointer :: up(:,:,:,:)
  
  call mgt_verify_n(mgt, "MGT_SET_UU", lev, n, lo, hi)

  up => dataptr(mgts(mgt)%uu, lev, n)
  up(lo(1):hi(1), lo(2):hi(2),1,1) = uu(lo(1):hi(1), lo(2):hi(2))

end subroutine mgt_set_uu_2d
subroutine mgt_set_uu_3d(mgt, lev, n, uu, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: mgt, lev, n, lo(3), hi(3), plo(3), phi(3)
  real(kind=dp_t), intent(in) :: uu(plo(1):phi(1), plo(2):phi(2), plo(3):phi(3))
  real(kind=dp_t), pointer :: up(:,:,:,:)
  
  call mgt_verify_n(mgt, "MGT_SET_UU", lev, n, lo, hi)

  up => dataptr(mgts(mgt)%uu, lev, n)
  up(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 1) = uu(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))

end subroutine mgt_set_uu_3d

subroutine mgt_get_uu_1d(mgt, lev, n, uu, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: mgt, lev, n, lo(1), hi(1), plo(1), phi(1)
  real(kind=dp_t), intent(inout) :: uu(plo(1):phi(1))
  real(kind=dp_t), pointer :: up(:,:,:,:)
  
  call mgt_verify_n(mgt, "MGT_GET_UU", lev, n, lo, hi)

  up => dataptr(mgts(mgt)%uu, lev, n)
  uu(lo(1):hi(1)) = up(lo(1):hi(1), 1,1,1)

end subroutine mgt_get_uu_1d
subroutine mgt_get_uu_2d(mgt, lev, n, uu, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: mgt, lev, n, lo(2), hi(2), plo(2), phi(2)
  real(kind=dp_t), intent(inout) :: uu(plo(1):phi(1), plo(2):phi(2))
  real(kind=dp_t), pointer :: up(:,:,:,:)
  
  call mgt_verify_n(mgt, "MGT_GET_UU", lev, n, lo, hi)

  up => dataptr(mgts(mgt)%uu, lev, n)
  uu(lo(1):hi(1), lo(2):hi(2)) = up(lo(1):hi(1), lo(2):hi(2),1,1)

end subroutine mgt_get_uu_2d
subroutine mgt_get_uu_3d(mgt, lev, n, uu, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: mgt, lev, n, lo(3), hi(3), plo(3), phi(3)
  real(kind=dp_t), intent(inout) :: uu(plo(1):phi(1), plo(2):phi(2), plo(3):phi(3))
  real(kind=dp_t), pointer :: up(:,:,:,:)
  
  call mgt_verify_n(mgt, "MGT_GET_UU", lev, n, lo, hi)

  up => dataptr(mgts(mgt)%uu, lev, n)
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

  do i = nlevels(mgts(mgt)%mla), 1, -1
     call destroy(mgts(mgt)%coeffs(i))
  end do
  call mg_tower_destroy(mgts(mgt)%mgt(1))
  call destroy(mgts(mgt)%rh)
  call destroy(mgts(mgt)%uu)
  call destroy(mgts(mgt)%mla)
  mgts(mgt)%han = -1
  mgts(mgt)%dim = 0
  mgt = -1

end subroutine mgt_dealloc

subroutine mgt_solve_cc(mgt)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: mgt

  call mgt_verify(mgt, "MGT_SOLVE")
  if ( .not. mgts(mgt)%final ) then
     call bl_error("MGT_SOLVE: MGT not finalized")
  end if

  call mg_tower_solve(mgts(mgt)%mgt(1), mgts(mgt)%uu%mf(1), mgts(mgt)%rh%mf(1))

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
