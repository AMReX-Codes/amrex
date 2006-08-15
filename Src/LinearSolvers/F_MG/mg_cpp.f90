module cpp_mg_module
  use mg_module
  use ml_layout_module
  use ml_multifab_module
  use bndry_reg_module

  implicit none

  type mg_server
     logical         :: final = .false.
     logical         :: nodal
     integer         :: dim  = 0
     integer         :: nlevel
!    integer         :: stencil_order = 3
!    integer         :: stencil_order = 2
     integer         :: stencil_order = 1
     integer         :: nu1, nu2, nuf, nub
     integer         :: gamma
     real(dp_t)      :: omega
     integer         :: max_iter
     integer         :: max_nlevel
     integer         :: min_width
     integer         :: smoother
     integer         :: cycle
     integer         :: verbose
     integer         :: cg_verbose
     integer         :: bottom_max_iter
     integer         :: bottom_solver
     real(dp_t)      :: bottom_solver_eps
     real(dp_t)      :: eps
     type(ml_layout) :: mla
     type(mg_tower)  :: mg_tower_default
     type(mg_tower), pointer :: mgt(:) => Null()
     type(box), pointer :: pd(:) => Null()
     integer, pointer :: bc(:,:) => Null()
     integer, pointer :: rr(:,:)
     type(multifab), pointer :: rh(:) => Null()
     type(multifab), pointer :: uu(:) => Null()
     type(multifab), pointer :: coeffs(:) => Null()
  end type mg_server

  type(mg_server), save   :: mgts

contains
  
  subroutine mgt_verify(str)
    character(len=*), intent(in) :: str

  2 if ( mgts%dim == 0 ) then
       call bl_error( trim(str) // ": MGT invalid DIM: not allocated: ")
    end if
    
  end subroutine mgt_verify

  subroutine mgt_verify_lev(str, lev)
    integer, intent(in) :: lev
    character(len=*), intent(in) :: str
    call mgt_verify(str)
    if ( lev < 1 .or. lev > mgts%nlevel ) then
       call bl_error( trim(str) // ": Level out of bounds", lev)
    end if
  end subroutine mgt_verify_lev

  subroutine mgt_verify_n(str, lev, n, lo, hi)
    integer, intent(in) :: lev, n, lo(:), hi(:)
    character(len=*), intent(in) :: str
    type(box) :: bx

    call mgt_verify_lev(str, lev)
    if ( n < 1 .or. n > nboxes(mgts%mla, lev) ) then
       call bl_error( trim(str) // ": Box out of bounds", n)
    end if
    bx = make_box(lo, hi)
    if ( bx /= get_box(mgts%mla, lev, n) ) then
       call bl_error( trim(str) // ": Box no filling")
    end if

  end subroutine mgt_verify_n

  subroutine mgt_not_final(str)
    character(len=*), intent(in) :: str

    call mgt_verify(str)
    if ( mgts%final ) then
       call bl_error( trim(str) // ": Changes made to finalized solver!")
    end if
  end subroutine mgt_not_final

end module cpp_mg_module

subroutine mgt_init ()
  use parallel
  logical, save :: first = .true.
  if ( first ) then
     call parallel_initialize()
     first = .false.
  end if
end subroutine mgt_init

subroutine mgt_alloc(dm, nlevel, nodal)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: dm, nlevel
  integer :: nodal
  integer i

  if ( mgts%dim == 0 ) then
     mgts%dim = dm
     mgts%nlevel = nlevel
     mgts%nodal = (nodal /= 0)
  end if

  allocate(mgts%rr(nlevel-1,dm))
  allocate(mgts%rh(nlevel))
  allocate(mgts%pd(nlevel))
  allocate(mgts%uu(nlevel))
  allocate(mgts%mgt(nlevel))

  call build(mgts%mla, nlevel, dm)

end subroutine mgt_alloc

subroutine mgt_set_level(lev, nb, dm, lo, hi, pd_lo, pd_hi, bc, pm, pmap)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: lev, nb, dm
  integer, intent(in) :: lo(nb,dm), hi(nb,dm), pd_lo(dm), pd_hi(dm), bc(2,dm), pm(dm), pmap(nb+1)

  type(box) :: bxs(nb)
  integer   :: i
  integer   :: nc
  logical   :: pmask(dm)
  type(box) :: pd
  integer   :: flev
  logical, allocatable :: nodal(:)

  flev = lev + 1
  call mgt_verify_lev("MGT_SET_LEVEL", flev)

  pmask = (pm /= 0)

  allocate(nodal(dm))

  if ( dm /= mgts%dim ) then
     call bl_error("MGT_SET_LEVEL: Input DIM doesn't match internal DIM")
  end if
  call build(mgts%pd(flev), pd_lo(1:dm), pd_hi(1:dm))
  do i = 1, nb
     bxs(i) = make_box(lo(i,:), hi(i,:))
  end do
  call build(mgts%mla%mba%bas(flev), bxs)
  call build(mgts%mla%la(flev),  &
       mgts%mla%mba%bas(flev), &
       mgts%pd(flev), pmask = pmask, &
       mapping = LA_EXPLICIT, explicit_mapping = pmap(1:nb))

!   if (parallel_IOProcessor()) call print(mgts%pd(flev))
!   if (parallel_IOProcessor()) call print(mgts%mla%la(flev))

  allocate(mgts%bc(dm,2))

  mgts%bc = transpose(bc)

end subroutine mgt_set_level

subroutine mgt_finalize(dx)
  use cpp_mg_module
  implicit none
  real(dp_t), intent(in) :: dx(mgts%nlevel,mgts%dim)
  integer :: i, dm, nlev, n
  integer :: ns
  integer :: nc
  logical, allocatable :: nodal(:)

  integer :: max_nlevel_in
  integer :: bottom_solver_in

  type(boxarray) :: bac

  integer :: bottom_max_iter_in

  call mgt_verify("MGT_FINALIZE")

  dm = mgts%dim

!  if (parallel_IOProcessor()) print *,'DX IN MGT_FINALIZE ', dx(1,:)

  nc = 1
  nlev = mgts%nlevel

  do i = 1, nlev-1
     mgts%rr(i,:) = mgts%mla%mba%rr(i,:)
  end do

  do i = 1, nlev
     call build(mgts%uu(i), mgts%mla%la(i), nc, ng = 1)
     call build(mgts%rh(i), mgts%mla%la(i), nc, ng = 0)
  end do

  do i = nlev-1, 1, -1
     call build(mgts%mla%mask(i), mgts%mla%la(i), nc = 1, ng = 0)
     call setval(mgts%mla%mask(i), val = .TRUE.)
     call copy(bac, mgts%mla%mba%bas(i+1))
     call boxarray_coarsen(bac, mgts%rr(i,:))
     call setval(mgts%mla%mask(i), .false., bac)
     call destroy(bac)
  end do

  allocate(nodal(1:dm))
  nodal = mgts%nodal
  ns = 1 + dm*3

  do n = nlev, 1, -1
     if ( n == 1 ) then
        max_nlevel_in = mgts%max_nlevel
        bottom_solver_in = mgts%bottom_solver
        bottom_max_iter_in = mgts%bottom_max_iter
     else
        if ( all(mgts%rr == 2) ) then
           max_nlevel_in = 1
        else if ( all(mgts%rr == 4) ) then
           max_nlevel_in = 2
        else
           call bl_error("MGT_FINALIZE: confused about ref_ratio")
        end if
        bottom_solver_in = 0
        bottom_max_iter_in = mgts%nu1
     end if

!     if ( parallel_IOProcessor() ) print *,'DX INTO MG_TOWER ', dx(n,1), dx(n,2)

     call mg_tower_build(mgts%mgt(n), mgts%mla%la(n), mgts%pd(n), mgts%bc, &
          dh                = dx(n,:), &
          ns                = ns, &
          smoother          = mgts%smoother, &
          nu1               = mgts%nu1, &
          nu2               = mgts%nu2, &
          nuf               = mgts%nuf, &
          nub               = mgts%nub, &
          gamma             = mgts%gamma, &
          cycle             = mgts%cycle, &
          omega             = mgts%omega, &
          bottom_solver     = bottom_solver_in, &
          bottom_max_iter   = bottom_max_iter_in, &
          bottom_solver_eps = mgts%bottom_solver_eps, &
          max_iter          = mgts%max_iter, &
          max_nlevel        = max_nlevel_in, &
          min_width         = mgts%min_width, &
          verbose           = mgts%verbose, &
          cg_verbose        = mgts%cg_verbose, &
          nodal             = nodal &
          )
  end do

  if ( parallel_IOProcessor() .and. mgts%verbose > 0) then
    call mg_tower_print(mgts%mgt(nlev))
  end if

end subroutine mgt_finalize

subroutine mgt_init_coeffs_lev(lev)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: lev
  integer :: nlev, dm, i
  integer :: flev
  flev = lev + 1
  call mgt_verify_lev("MGT_INIT_STENCIL_LEV", flev)

  dm = mgts%dim
  nlev = mgts%mgt(flev)%nlevels
  allocate(mgts%coeffs(nlev))

  call  build(mgts%coeffs(nlev), mgts%mgt(flev)%ss(nlev)%la, 1+dm, 1)
  call setval(mgts%coeffs(nlev), 0.0_dp_t, all=.true.)
  do i = nlev - 1, 1, -1
     call build(mgts%coeffs(i), mgts%mgt(flev)%ss(i)%la, 1+dm, 1)
     call setval(mgts%coeffs(i), 0.0_dp_t, all=.true.)
  end do

end subroutine mgt_init_coeffs_lev

subroutine mgt_finalize_stencil_lev(lev, xa, xb, pxa, pxb)
  use cpp_mg_module
  use coeffs_module
  implicit none
  integer, intent(in) :: lev
  real(dp_t), intent(in) :: xa(*), xb(*), pxa(*), pxb(*)
  integer :: nlev, i, dm
  integer :: flev
  type(box) :: pd
  type(boxarray) :: pdv
  dm = mgts%dim
  flev = lev + 1
  call mgt_verify_lev("MGT_SET_COEFS_LEV", flev)
  
  pd = mgts%pd(flev)
  nlev = mgts%mgt(flev)%nlevels
  do i = nlev-1, 1, -1
     call coarsen_coeffs(mgts%coeffs(i+1), mgts%coeffs(i))
  end do
  do i = nlev, 1, -1
     pdv = layout_boxarray(mgts%mgt(flev)%ss(i)%la)
     call stencil_fill_cc(mgts%mgt(flev)%ss(i), mgts%coeffs(i), mgts%mgt(flev)%dh(:,i), &
          pdv, mgts%mgt(flev)%mm(i), xa(1:dm), xb(1:dm), pxa(1:dm), pxb(1:dm), &
          mgts%mgt(flev)%pd(i), &
          mgts%stencil_order, mgts%bc)

     if (i .eq. 1 .and. mgts%bottom_solver .eq. 3) then
        call copy(mgts%mgt(i)%ss1, mgts%mgt(i)%ss(1))
        call copy(mgts%mgt(i)%mm1, mgts%mgt(i)%mm(1))
        if ( parallel_IOProcessor() ) then
           call sparse_build(mgts%mgt(i)%sparse_object, mgts%mgt(i)%ss1, &
                mgts%mgt(i)%mm1, mgts%mgt(i)%ss1%la, mgts%stencil_order, mgts%mgt(i)%verbose)
        end if
     end if

     call destroy(mgts%coeffs(i))

  end do

  deallocate(mgts%coeffs)

end subroutine mgt_finalize_stencil_lev

subroutine mgt_finalize_stencil()
   use cpp_mg_module
   implicit none
   call mgt_verify("MGT_FINALIZE_STENCIL")
   mgts%final = .true.
end subroutine mgt_finalize_stencil

subroutine mgt_set_rh_1d(lev, n, rh, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(1), hi(1), plo(1), phi(1)
  real(kind=dp_t), intent(in) :: rh(plo(1):phi(1))
  real(kind=dp_t), pointer :: rp(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1
  call mgt_verify_n("MGT_SET_RH", flev, fn, lo, hi)

  rp => dataptr(mgts%rh(flev), fn)
  rp(lo(1):hi(1), 1,1,1) = rh(lo(1):hi(1))

end subroutine mgt_set_rh_1d
subroutine mgt_set_rh_2d(lev, n, rh, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(2), hi(2), plo(2), phi(2)
  real(kind=dp_t), intent(in) :: rh(plo(1):phi(1), plo(2):phi(2))
  real(kind=dp_t), pointer :: rp(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1
  
  call mgt_verify_n("MGT_SET_RH", flev, fn, lo, hi)

  rp => dataptr(mgts%rh(flev), fn)
  rp(lo(1):hi(1), lo(2):hi(2),1,1) = rh(lo(1):hi(1), lo(2):hi(2))

end subroutine mgt_set_rh_2d
subroutine mgt_set_rh_3d(lev, n, rh, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(3), hi(3), plo(3), phi(3)
  real(kind=dp_t), intent(in) :: rh(plo(1):phi(1), plo(2):phi(2), plo(3):phi(3))
  real(kind=dp_t), pointer :: rp(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1
  
  call mgt_verify_n("MGT_SET_RH", flev, fn, lo, hi)

  rp => dataptr(mgts%rh(flev), fn)
  rp(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3),1) = rh(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))

end subroutine mgt_set_rh_3d

subroutine mgt_set_cf_1d(lev, n, cf, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(1), hi(1), plo(1), phi(1)
  real(kind=dp_t), intent(in) :: cf(plo(1):phi(1),*)
  real(kind=dp_t), pointer :: cp(:,:,:,:)
  integer :: flev, fn, nlev
  fn = n + 1
  flev = lev+1
  nlev = size(mgts%coeffs)
  call mgt_verify_n("MGT_SET_UU", lev, n, lo, hi)

  cp => dataptr(mgts%coeffs(nlev), fn)
  cp(lo(1):hi(1), 1,1, 1:4) = cf(lo(1):hi(1), 1:4)

end subroutine mgt_set_cf_1d

subroutine mgt_set_uu_1d(lev, n, uu, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(1), hi(1), plo(1), phi(1)
  real(kind=dp_t), intent(in) :: uu(plo(1):phi(1))
  real(kind=dp_t), pointer :: up(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1
  
  call mgt_verify_n("MGT_SET_UU", lev, n, lo, hi)

  up => dataptr(mgts%uu(lev), fn)
  up(lo(1)-1:hi(1)+1, 1,1,1) = uu(lo(1)-1:hi(1)+1)

end subroutine mgt_set_uu_1d

subroutine mgt_set_cfa_2d(lev, n, cf, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(2), hi(2), plo(2), phi(2)
  real(kind=dp_t), intent(in) :: cf(plo(1):phi(1), plo(2):phi(2), *)
  real(kind=dp_t), pointer :: cp(:,:,:,:)
  integer :: flev, fn, nlev, i, j

  fn = n + 1
  flev = lev+1
  nlev = size(mgts%coeffs)
  call mgt_verify_n("MGT_SET_CF", flev, fn, lo, hi)

  cp => dataptr(mgts%coeffs(nlev), fn)
  cp(lo(1):hi(1), lo(2):hi(2), 1, 1) = cf(lo(1):hi(1), lo(2):hi(2), 1)

end subroutine mgt_set_cfa_2d

subroutine mgt_set_cfbx_2d(lev, n, cf, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(2), hi(2), plo(2), phi(2)
  real(kind=dp_t), intent(in) :: cf(plo(1):phi(1), plo(2):phi(2))
  real(kind=dp_t), pointer :: cp(:,:,:,:)
  integer :: flev, fn, nlev, i, j

  fn = n + 1
  flev = lev+1
  nlev = size(mgts%coeffs)
  call mgt_verify_n("MGT_SET_CF", flev, fn, lo, hi)

  cp => dataptr(mgts%coeffs(nlev), fn)
  cp(lo(1):hi(1)+1, lo(2):hi(2), 1, 2) = cf(lo(1):hi(1)+1, lo(2):hi(2))

end subroutine mgt_set_cfbx_2d

subroutine mgt_set_cfby_2d(lev, n, cf, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(2), hi(2), plo(2), phi(2)
  real(kind=dp_t), intent(in) :: cf(plo(1):phi(1), plo(2):phi(2))
  real(kind=dp_t), pointer :: cp(:,:,:,:)
  integer :: flev, fn, nlev, i, j 
  fn = n + 1
  flev = lev+1
  nlev = size(mgts%coeffs)
  call mgt_verify_n("MGT_SET_CF", flev, fn, lo, hi)

  cp => dataptr(mgts%coeffs(nlev), fn)
  cp(lo(1):hi(1), lo(2):hi(2)+1, 1, 3) = cf(lo(1):hi(1), lo(2):hi(2)+1)

end subroutine mgt_set_cfby_2d

subroutine mgt_set_uu_2d(lev, n, uu, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(2), hi(2), plo(2), phi(2)
  real(kind=dp_t), intent(in) :: uu(plo(1):phi(1), plo(2):phi(2))
  real(kind=dp_t), pointer :: up(:,:,:,:)
  integer :: flev, fn, i, j
  fn = n + 1
  flev = lev+1

  call mgt_verify_n("MGT_SET_UU", flev, fn, lo, hi)

  up => dataptr(mgts%uu(flev), fn)
  up(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1,1,1) = uu(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1)

end subroutine mgt_set_uu_2d

subroutine mgt_set_cfa_3d(lev, n, cf, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(3), hi(3), plo(3), phi(3)
  real(kind=dp_t), intent(in) :: cf(plo(1):phi(1), plo(2):phi(2), plo(3):phi(3))
  real(kind=dp_t), pointer :: cp(:,:,:,:)
  integer :: flev, fn, nlev, i, j, k
  fn = n + 1
  flev = lev+1
  nlev = size(mgts%coeffs)
  call mgt_verify_n("MGT_SET_CF", flev, fn, lo, hi)

  cp => dataptr(mgts%coeffs(nlev), fn)
  cp(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 1) = cf(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))

end subroutine mgt_set_cfa_3d

subroutine mgt_set_cfbx_3d(lev, n, cf, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(3), hi(3), plo(3), phi(3)
  real(kind=dp_t), intent(in) :: cf(plo(1):phi(1), plo(2):phi(2), plo(3):phi(3))
  real(kind=dp_t), pointer :: cp(:,:,:,:)
  integer :: flev, fn, nlev, i, j, k
  fn = n + 1
  flev = lev+1
  nlev = size(mgts%coeffs)
  call mgt_verify_n("MGT_SET_CF", flev, fn, lo, hi)

  cp => dataptr(mgts%coeffs(nlev), fn)
  cp(lo(1):hi(1)+1, lo(2):hi(2), lo(3):hi(3), 2) = cf(lo(1):hi(1)+1, lo(2):hi(2), lo(3):hi(3))

end subroutine mgt_set_cfbx_3d

subroutine mgt_set_cfby_3d(lev, n, cf, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(3), hi(3), plo(3), phi(3)
  real(kind=dp_t), intent(in) :: cf(plo(1):phi(1), plo(2):phi(2), plo(3):phi(3))
  real(kind=dp_t), pointer :: cp(:,:,:,:)
  integer :: flev, fn, nlev, i, j, k 
  fn = n + 1
  flev = lev+1
  nlev = size(mgts%coeffs)
  call mgt_verify_n("MGT_SET_CF", flev, fn, lo, hi)

  cp => dataptr(mgts%coeffs(nlev), fn)
  cp(lo(1):hi(1), lo(2):hi(2)+1, lo(3):hi(3), 3) = cf(lo(1):hi(1), lo(2):hi(2)+1, lo(3):hi(3))

end subroutine mgt_set_cfby_3d

subroutine mgt_set_cfbz_3d(lev, n, cf, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(3), hi(3), plo(3), phi(3)
  real(kind=dp_t), intent(in) :: cf(plo(1):phi(1), plo(2):phi(2), plo(3):phi(3))
  real(kind=dp_t), pointer :: cp(:,:,:,:)
  integer :: flev, fn, nlev, i, j, k 
  fn = n + 1
  flev = lev+1
  nlev = size(mgts%coeffs)
  call mgt_verify_n("MGT_SET_CF", flev, fn, lo, hi)

  cp => dataptr(mgts%coeffs(nlev), fn)
  cp(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)+1, 4) = cf(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)+1)

end subroutine mgt_set_cfbz_3d

subroutine mgt_set_uu_3d(lev, n, uu, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(3), hi(3), plo(3), phi(3)
  real(kind=dp_t), intent(in) :: uu(plo(1):phi(1), plo(2):phi(2), plo(3):phi(3))
  real(kind=dp_t), pointer :: up(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1
  
  call mgt_verify_n("MGT_SET_UU", flev, fn, lo, hi)

  up => dataptr(mgts%uu(flev), fn)
  up(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, lo(3)-1:hi(3)+1, 1) = &
     uu(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, lo(3)-1:hi(3)+1)

end subroutine mgt_set_uu_3d

subroutine mgt_get_uu_1d(lev, n, uu, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(1), hi(1), plo(1), phi(1)
  real(kind=dp_t), intent(inout) :: uu(plo(1):phi(1))
  real(kind=dp_t), pointer :: up(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1
  
  call mgt_verify_n("MGT_GET_UU", flev, fn, lo, hi)

  up => dataptr(mgts%uu(flev), fn)
  uu(lo(1):hi(1)) = up(lo(1):hi(1), 1,1,1)

end subroutine mgt_get_uu_1d
subroutine mgt_get_uu_2d(lev, n, uu, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(2), hi(2), plo(2), phi(2)
  real(kind=dp_t), intent(inout) :: uu(plo(1):phi(1), plo(2):phi(2))
  real(kind=dp_t), pointer :: up(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1
  
  call mgt_verify_n("MGT_GET_UU", flev, fn, lo, hi)

  up => dataptr(mgts%uu(flev), fn)
  uu(lo(1):hi(1), lo(2):hi(2)) = up(lo(1):hi(1), lo(2):hi(2),1,1)

end subroutine mgt_get_uu_2d
subroutine mgt_get_uu_3d(lev, n, uu, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(3), hi(3), plo(3), phi(3)
  real(kind=dp_t), intent(inout) :: uu(plo(1):phi(1), plo(2):phi(2), plo(3):phi(3))
  real(kind=dp_t), pointer :: up(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1
  
  call mgt_verify_n("MGT_GET_UU", flev, fn, lo, hi)

  up => dataptr(mgts%uu(flev), fn)
  uu(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)) = up(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 1)

end subroutine mgt_get_uu_3d

subroutine mgt_dealloc()
  use cpp_mg_module
  implicit none
  integer :: i
  
  call mgt_verify("MGT_DEALLOC")
  if ( .not. mgts%final ) then
     call bl_error("MGT_DEALLOC: MGT not finalized")
  end if

  do i = 1, mgts%nlevel
    call mg_tower_destroy(mgts%mgt(i))
  end do

  do i = mgts%nlevel, 1, -1
     call destroy(mgts%rh(i))
     call destroy(mgts%uu(i))
  end do
  call destroy(mgts%mla)
  mgts%dim = 0
  mgts%final = .false.

end subroutine mgt_dealloc

subroutine mgt_solve(tol,abs_tol)
  use cpp_mg_module
  use ml_cc_module
  use fabio_module
  implicit none
  real(kind=dp_t), intent(in) :: tol, abs_tol

  integer :: do_diagnostics

  call mgt_verify("MGT_SOLVE")
  if ( .not. mgts%final ) then
     call bl_error("MGT_SOLVE: MGT not finalized")
  end if

  do_diagnostics = 0
  call ml_cc(mgts%mla, mgts%mgt, &
       mgts%rh, mgts%uu, &
       mgts%mla%mask, mgts%rr, &
       do_diagnostics, tol)

end subroutine mgt_solve

subroutine mgt_set_defaults(nu_1,nu_2,nu_b,nu_f,gamma,omega,max_iter,bottom_max_iter, &
                            bottom_solver,bottom_solver_eps, &
                            verbose,cg_verbose,max_nlevel,min_width,cycle,smoother)
  use cpp_mg_module
  implicit none
  integer   , intent(in) :: nu_1,nu_2,nu_b,nu_f,gamma,max_iter,bottom_max_iter,bottom_solver
  integer   , intent(in) :: verbose, cg_verbose, max_nlevel, min_width, cycle, smoother
  real(dp_t), intent(in) :: omega, bottom_solver_eps

  call mgt_not_final("MGT_SET_DEFAULTS")

  mgts%nu1             = nu_1
  mgts%nu2             = nu_2
  mgts%nuf             = nu_f
  mgts%nub             = nu_b
  mgts%gamma           = gamma
  mgts%omega           = omega
  mgts%max_iter        = max_iter
  mgts%verbose         = verbose
  mgts%cg_verbose      = cg_verbose
  mgts%smoother        = smoother
  mgts%cycle           = cycle
  mgts%bottom_max_iter = bottom_max_iter
  mgts%bottom_solver   = bottom_solver
  mgts%bottom_solver_eps = bottom_solver_eps
  mgts%max_nlevel      = max_nlevel
  mgts%min_width       = min_width

end subroutine mgt_set_defaults

subroutine mgt_get_defaults(nu_1,nu_2,nu_b,nu_f,gamma,omega,max_iter,bottom_max_iter, &
                            bottom_solver, &
                            verbose,cg_verbose,max_nlevel,min_width,cycle,smoother)
  use cpp_mg_module
  implicit none
  integer   , intent(out) :: nu_1,nu_2,nu_b,nu_f,gamma,max_iter,bottom_max_iter,bottom_solver
  integer   , intent(out) :: verbose, cg_verbose, max_nlevel, min_width, cycle, smoother
  real(dp_t), intent(out) :: omega

  nu_1       = mgts%mg_tower_default%nu1
  nu_2       = mgts%mg_tower_default%nu2
  nu_f       = mgts%mg_tower_default%nuf
  nu_b       = mgts%mg_tower_default%nub
  gamma      = mgts%mg_tower_default%gamma
  omega      = mgts%mg_tower_default%omega
  max_iter   = mgts%mg_tower_default%max_iter
  verbose    = mgts%mg_tower_default%verbose
  cg_verbose = mgts%mg_tower_default%cg_verbose
  smoother   = mgts%mg_tower_default%smoother
  cycle      = mgts%mg_tower_default%cycle

  bottom_max_iter = mgts%mg_tower_default%bottom_max_iter
  bottom_solver    = mgts%mg_tower_default%bottom_solver
  max_nlevel      = mgts%mg_tower_default%max_nlevel
  min_width       = mgts%mg_tower_default%min_width

end subroutine mgt_get_defaults
