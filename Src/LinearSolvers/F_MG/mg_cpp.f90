module cpp_mg_module

  use mg_module
  use multifab_module
  use ml_layout_module

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
     integer         :: cycle_type
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
     type(multifab), pointer :: res(:) => Null()
     type(multifab), pointer :: uu(:) => Null()
     type(multifab), pointer :: gp(:,:) => Null()
     type(multifab), pointer :: coeffs(:) => Null()
  end type mg_server

  type(mg_server), save   :: mgts

contains
  
  subroutine mgt_verify(str)
    character(len=*), intent(in) :: str

    if ( mgts%dim == 0 ) then
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

  if ( mgts%dim == 0 ) then
     mgts%dim = dm
     mgts%nlevel = nlevel
     mgts%nodal = (nodal /= 0)
  end if

  allocate(mgts%rr(nlevel-1,dm))
  allocate(mgts%rh(nlevel))
  allocate(mgts%res(nlevel))
  allocate(mgts%pd(nlevel))
  allocate(mgts%uu(nlevel))
  allocate(mgts%gp(nlevel,dm))
  allocate(mgts%mgt(nlevel))
  allocate(mgts%bc(dm,2))

  call build(mgts%mla, nlevel, dm)

end subroutine mgt_alloc

subroutine mgt_set_level(lev, nb, dm, lo, hi, pd_lo, pd_hi, pm, pmap)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: lev, nb, dm
  integer, intent(in) :: lo(nb,dm), hi(nb,dm), pd_lo(dm), pd_hi(dm), pm(dm), pmap(nb+1)

  type(box) :: bxs(nb)
  integer   :: i
  logical   :: pmask(dm)
  integer   :: flev

  flev = lev + 1
  call mgt_verify_lev("MGT_SET_LEVEL", flev)

  pmask = (pm /= 0)

  if ( dm /= mgts%dim ) then
     call bl_error("MGT_SET_LEVEL: Input DIM doesn't match internal DIM")
  end if
  call build(mgts%mla%mba%pd(flev), pd_lo(1:dm), pd_hi(1:dm))
  call build(mgts%pd(flev), pd_lo(1:dm), pd_hi(1:dm))
  if (flev > 1) then
    do i = 1, dm
      mgts%mla%mba%rr(flev-1,i) = extent(mgts%mla%mba%pd(flev),i) / extent(mgts%mla%mba%pd(flev-1),i)
    end do
  end if
  do i = 1, nb 
    bxs(i) = make_box(lo(i,:), hi(i,:))
  end do
  call build(mgts%mla%mba%bas(flev), bxs)
  call build(mgts%mla%la(flev),  &
       mgts%mla%mba%bas(flev), &
       mgts%pd(flev), pmask = pmask, &
       mapping = LA_EXPLICIT, explicit_mapping = pmap(1:nb))

end subroutine mgt_set_level

subroutine mgt_finalize(dx,bc)
  use cpp_mg_module
  implicit none
  real(dp_t), intent(in) :: dx(mgts%nlevel,mgts%dim)
  integer   , intent(in) :: bc(2,mgts%dim)
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

  nc = 1
  nlev = mgts%nlevel

  mgts%bc = transpose(bc)

  do i = 1, nlev-1
     mgts%rr(i,:) = mgts%mla%mba%rr(i,:)
  end do

  do i = 1, nlev
     call build(mgts%uu(i) , mgts%mla%la(i), nc, ng = 1)
     call build(mgts%rh(i) , mgts%mla%la(i), nc, ng = 0)
     call build(mgts%res(i), mgts%mla%la(i), nc, ng = 0)
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
     call mg_tower_build(mgts%mgt(n), mgts%mla%la(n), mgts%pd(n), mgts%bc, &
          dh                = dx(n,:), &
          ns                = ns, &
          smoother          = mgts%smoother, &
          nu1               = mgts%nu1, &
          nu2               = mgts%nu2, &
          nuf               = mgts%nuf, &
          nub               = mgts%nub, &
          gamma             = mgts%gamma, &
          cycle_type        = mgts%cycle_type, &
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

! if ( parallel_IOProcessor() .and. mgts%verbose > 1) then
!   call mg_tower_print(mgts%mgt(nlev))
! end if

end subroutine mgt_finalize

subroutine mgt_init_coeffs_lev(lev)
  use cpp_mg_module
  implicit none
  integer          , intent(in) :: lev

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

subroutine mgt_init_mc_coeffs_lev(lev,nc)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: lev,nc

  integer :: nlev, dm, i
  integer :: flev
  flev = lev + 1
  call mgt_verify_lev("MGT_INIT_STENCIL_LEV", flev)

  dm = mgts%dim
  nlev = mgts%mgt(flev)%nlevels
  allocate(mgts%coeffs(nlev))

  ! If nc = 1 then the first coefficient is alpha, the next coefficient
  !   is betax, the next is betay, etc.
  ! If nc > 1 then the first coefficient is alpha, the next nc coefficients
  !   are betax with components 1:nc, the next nc coefficients
  !   are betay with components 1:nc, etc.

  call  build(mgts%coeffs(nlev), mgts%mgt(flev)%ss(nlev)%la, 1+nc*dm, 1)
  call setval(mgts%coeffs(nlev), 0.0_dp_t, all=.true.)
  do i = nlev - 1, 1, -1
     call build(mgts%coeffs(i), mgts%mgt(flev)%ss(i)%la, 1+nc*dm, 1)
     call setval(mgts%coeffs(i), 0.0_dp_t, all=.true.)
  end do

end subroutine mgt_init_mc_coeffs_lev

subroutine mgt_finalize_stencil_lev(lev, xa, xb, pxa, pxb)
  use cpp_mg_module
  use coeffs_module
  implicit none
  integer, intent(in) :: lev
  real(dp_t), intent(in) :: xa(*), xb(*), pxa(*), pxb(*)
  integer :: nlev, i, dm
  integer :: flev
  type(boxarray) :: pdv
  dm = mgts%dim
  flev = lev + 1
  call mgt_verify_lev("MGT_SET_COEFS_LEV", flev)
  
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

subroutine mgt_set_cfa_1d(lev, n, cf, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(2), hi(2), plo(2), phi(2)
  real(kind=dp_t), intent(in) :: cf(plo(1):phi(1), *)
  real(kind=dp_t), pointer :: cp(:,:,:,:)
  integer :: flev, fn, nlev

  fn = n + 1
  flev = lev+1
  nlev = size(mgts%coeffs)
  call mgt_verify_n("MGT_SET_CFA_2D", flev, fn, lo, hi)

  cp => dataptr(mgts%coeffs(nlev), fn)
  cp(lo(1):hi(1), 1, 1, 1) = cf(lo(1):hi(1), 1)

end subroutine mgt_set_cfa_1d

subroutine mgt_set_cfa_1d_const(lev, n, lo, hi, coeff_value)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(1), hi(1)
  real(kind=dp_t), intent(in) :: coeff_value
  real(kind=dp_t), pointer :: cp(:,:,:,:)
  integer :: flev, fn, nlev
  fn = n + 1
  flev = lev+1
  nlev = size(mgts%coeffs)
  call mgt_verify_n("MGT_SET_CFA_1D_CONST", flev, fn, lo, hi)

  cp => dataptr(mgts%coeffs(nlev), fn)
  cp(lo(1):hi(1), 1, 1, 1) = coeff_value

end subroutine mgt_set_cfa_1d_const

subroutine mgt_set_cfbx_1d(lev, n, cf, b, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none 
  integer, intent(in) :: lev, n, lo(1), hi(1), plo(1), phi(1)
  real(kind=dp_t), intent(in) :: cf(plo(1):phi(1))
  real(kind=dp_t), intent(in) :: b
  real(kind=dp_t), pointer :: cp(:,:,:,:)
  integer :: flev, fn, nlev

  fn = n + 1
  flev = lev+1
  nlev = size(mgts%coeffs)
  call mgt_verify_n("MGT_SET_CFBX_2D", flev, fn, lo, hi)

  cp => dataptr(mgts%coeffs(nlev), fn)
  cp(lo(1):hi(1)+1, 1, 1, 2) = b * cf(lo(1):hi(1)+1)

end subroutine mgt_set_cfbx_1d

subroutine mgt_set_cfbnx_1d(lev, n, cf, b, plo, phi, lo, hi, nc)
  use cpp_mg_module
  implicit none 
  integer, intent(in) :: lev, n, lo(1), hi(1), plo(1), phi(1), nc
  real(kind=dp_t), intent(in) :: cf(plo(1):phi(1),nc)
  real(kind=dp_t), intent(in) :: b
  real(kind=dp_t), pointer :: cp(:,:,:,:)
  integer :: flev, fn, nlev

  fn = n + 1
  flev = lev+1
  nlev = size(mgts%coeffs)
  call mgt_verify_n("MGT_SET_CFBX_2D", flev, fn, lo, hi)

  cp => dataptr(mgts%coeffs(nlev), fn)
  cp(lo(1):hi(1)+1, 1, 1, 2:nc+1) = b * cf(lo(1):hi(1)+1,1:nc)

end subroutine mgt_set_cfbnx_1d

subroutine mgt_set_cfbx_1d_const(lev, n, lo, hi, coeff_value)
  use cpp_mg_module
  implicit none
  integer        , intent(in) :: lev, n, lo(1), hi(1)
  real(kind=dp_t), intent(in) :: coeff_value
  real(kind=dp_t), pointer :: cp(:,:,:,:)
  integer :: flev, fn, nlev

  fn = n + 1
  flev = lev+1
  nlev = size(mgts%coeffs)
  call mgt_verify_n("MGT_SET_CFBX_1D_CONST", flev, fn, lo, hi)

  cp => dataptr(mgts%coeffs(nlev), fn)
  cp(lo(1):hi(1)+1, 1, 1, 2) = coeff_value

end subroutine mgt_set_cfbx_1d_const

subroutine mgt_set_cfa_2d(lev, n, cf, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(2), hi(2), plo(2), phi(2)
  real(kind=dp_t), intent(in) :: cf(plo(1):phi(1), plo(2):phi(2), *)
  real(kind=dp_t), pointer :: cp(:,:,:,:)
  integer :: flev, fn, nlev

  fn = n + 1
  flev = lev+1
  nlev = size(mgts%coeffs)
  call mgt_verify_n("MGT_SET_CFA_2D", flev, fn, lo, hi)

  cp => dataptr(mgts%coeffs(nlev), fn)
  cp(lo(1):hi(1), lo(2):hi(2), 1, 1) = cf(lo(1):hi(1), lo(2):hi(2), 1)

end subroutine mgt_set_cfa_2d

subroutine mgt_set_cfa_2d_const(lev, n, lo, hi, coeff_value)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(2), hi(2)
  real(kind=dp_t), intent(in) :: coeff_value
  real(kind=dp_t), pointer :: cp(:,:,:,:)
  integer :: flev, fn, nlev

  fn = n + 1
  flev = lev+1
  nlev = size(mgts%coeffs)
  call mgt_verify_n("MGT_SET_CFA_2D_CONST", flev, fn, lo, hi)

  cp => dataptr(mgts%coeffs(nlev), fn)
  cp(lo(1):hi(1), lo(2):hi(2), 1, 1) = coeff_value

end subroutine mgt_set_cfa_2d_const

subroutine mgt_set_cfbx_2d(lev, n, cf, b, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(2), hi(2), plo(2), phi(2)
  real(kind=dp_t), intent(in) :: cf(plo(1):phi(1), plo(2):phi(2))
  real(kind=dp_t), intent(in) :: b
  real(kind=dp_t), pointer :: cp(:,:,:,:)
  integer :: flev, fn, nlev
  integer :: i,j

  fn = n + 1
  flev = lev+1 
  nlev = size(mgts%coeffs)
  call mgt_verify_n("MGT_SET_CFBX_2D", flev, fn, lo, hi)

  cp => dataptr(mgts%coeffs(nlev), fn)
  cp(lo(1):hi(1)+1, lo(2):hi(2), 1, 2) = b * cf(lo(1):hi(1)+1, lo(2):hi(2))
  do j = lo(2),hi(2)
  do i = lo(1),hi(1)+1
     print *,'CF ',cf(i,j)
  end do
  end do

end subroutine mgt_set_cfbx_2d

subroutine mgt_set_cfbnx_2d(lev, n, cf, b, plo, phi, lo, hi, nc)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, nc, lo(2), hi(2), plo(2), phi(2)
  real(kind=dp_t), intent(in) :: cf(plo(1):phi(1), plo(2):phi(2), nc)
  real(kind=dp_t), intent(in) :: b
  real(kind=dp_t), pointer :: cp(:,:,:,:)
  integer :: flev, fn, nlev
  integer :: i,j
  fn = n + 1
  flev = lev+1 
  nlev = size(mgts%coeffs)
  call mgt_verify_n("MGT_SET_CFBX_2D", flev, fn, lo, hi)

  ! betax(i,j,1:nc) goes into coeffs(i,j,2:nc+1)
  cp => dataptr(mgts%coeffs(nlev), fn)
  cp(lo(1):hi(1)+1, lo(2):hi(2), 1, 2:nc+1) = b * cf(lo(1):hi(1)+1, lo(2):hi(2),1:nc)
  do j = lo(2),hi(2)
  do i = lo(1),hi(1)+1
     print *,'CF ',cf(i,j,1)+cf(i,j,2)
  end do
  end do

end subroutine mgt_set_cfbnx_2d

subroutine mgt_set_cfbx_2d_const(lev, n, lo, hi, coeff_value)
  use cpp_mg_module
  implicit none
  integer        , intent(in) :: lev, n, lo(2), hi(2)
  real(kind=dp_t), intent(in) :: coeff_value
  real(kind=dp_t), pointer :: cp(:,:,:,:)
  integer :: flev, fn, nlev

  fn = n + 1
  flev = lev+1
  nlev = size(mgts%coeffs)
  call mgt_verify_n("MGT_SET_CFBX_2D_CONST", flev, fn, lo, hi)

  cp => dataptr(mgts%coeffs(nlev), fn)
  cp(lo(1):hi(1)+1, lo(2):hi(2), 1, 2) = coeff_value

end subroutine mgt_set_cfbx_2d_const

subroutine mgt_set_cfby_2d(lev, n, cf, b, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(2), hi(2), plo(2), phi(2)
  real(kind=dp_t), intent(in) :: cf(plo(1):phi(1), plo(2):phi(2))
  real(kind=dp_t), intent(in) :: b
  real(kind=dp_t), pointer :: cp(:,:,:,:)
  integer :: flev, fn, nlev
  fn = n + 1
  flev = lev+1
  nlev = size(mgts%coeffs)
  call mgt_verify_n("MGT_SET_CFBY_2D", flev, fn, lo, hi)

  cp => dataptr(mgts%coeffs(nlev), fn)
  cp(lo(1):hi(1), lo(2):hi(2)+1, 1, 3) = b * cf(lo(1):hi(1), lo(2):hi(2)+1)

end subroutine mgt_set_cfby_2d

subroutine mgt_set_cfbny_2d(lev, n, cf, b, plo, phi, lo, hi, nc)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, nc, lo(2), hi(2), plo(2), phi(2)
  real(kind=dp_t), intent(in) :: cf(plo(1):phi(1), plo(2):phi(2), nc)
  real(kind=dp_t), intent(in) :: b
  real(kind=dp_t), pointer :: cp(:,:,:,:)
  integer :: flev, fn, nlev

  fn = n + 1
  flev = lev+1
  nlev = size(mgts%coeffs)
  call mgt_verify_n("MGT_SET_CFBY_2D", flev, fn, lo, hi)

  ! betay(i,j,1:nc) goes into coeffs(i,j,nc+2:2*nc+1)
  cp => dataptr(mgts%coeffs(nlev), fn)
  cp(lo(1):hi(1), lo(2):hi(2)+1, 1, nc+2:2*nc+1) = b * cf(lo(1):hi(1), lo(2):hi(2)+1, 1:nc)

end subroutine mgt_set_cfbny_2d


subroutine mgt_set_cfby_2d_const(lev, n, lo, hi, coeff_value)
  use cpp_mg_module
  implicit none
  integer        , intent(in) :: lev, n, lo(2), hi(2)
  real(kind=dp_t), intent(in) :: coeff_value
  real(kind=dp_t), pointer :: cp(:,:,:,:)
  integer :: flev, fn, nlev
  fn = n + 1
  flev = lev+1
  nlev = size(mgts%coeffs)
  call mgt_verify_n("MGT_SET_CFBY_2D_CONST", flev, fn, lo, hi)

  cp => dataptr(mgts%coeffs(nlev), fn)
  cp(lo(1):hi(1), lo(2):hi(2)+1, 1, 3) = coeff_value

end subroutine mgt_set_cfby_2d_const

subroutine mgt_set_cfa_3d(lev, n, cf, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(3), hi(3), plo(3), phi(3)
  real(kind=dp_t), intent(in) :: cf(plo(1):phi(1), plo(2):phi(2), plo(3):phi(3))
  real(kind=dp_t), pointer :: cp(:,:,:,:)
  integer :: flev, fn, nlev
  fn = n + 1
  flev = lev+1
  nlev = size(mgts%coeffs)
  call mgt_verify_n("MGT_SET_CF", flev, fn, lo, hi)

  cp => dataptr(mgts%coeffs(nlev), fn)
  cp(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 1) = cf(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))

end subroutine mgt_set_cfa_3d

subroutine mgt_set_cfa_3d_const(lev, n, lo, hi, coeff_value)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(3), hi(3)
  real(kind=dp_t), intent(in) :: coeff_value
  real(kind=dp_t), pointer :: cp(:,:,:,:)
  integer :: flev, fn, nlev
  fn = n + 1
  flev = lev+1
  nlev = size(mgts%coeffs)
  call mgt_verify_n("MGT_SET_CF", flev, fn, lo, hi)

  cp => dataptr(mgts%coeffs(nlev), fn)
  cp(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 1) = coeff_value

end subroutine mgt_set_cfa_3d_const

subroutine mgt_set_cfbx_3d(lev, n, cf, b, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(3), hi(3), plo(3), phi(3)
  real(kind=dp_t), intent(in) :: cf(plo(1):phi(1), plo(2):phi(2), plo(3):phi(3))
  real(kind=dp_t), intent(in) :: b
  real(kind=dp_t), pointer :: cp(:,:,:,:)
  integer :: flev, fn, nlev
  fn = n + 1
  flev = lev+1
  nlev = size(mgts%coeffs)
  call mgt_verify_n("MGT_SET_CF", flev, fn, lo, hi)

  cp => dataptr(mgts%coeffs(nlev), fn)
  cp(lo(1):hi(1)+1, lo(2):hi(2), lo(3):hi(3), 2) = b * cf(lo(1):hi(1)+1, lo(2):hi(2), lo(3):hi(3))

end subroutine mgt_set_cfbx_3d

subroutine mgt_set_cfbnx_3d(lev, n, cf, b, plo, phi, lo, hi, nc)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, nc, lo(3), hi(3), plo(3), phi(3)
  real(kind=dp_t), intent(in) :: cf(plo(1):phi(1), plo(2):phi(2), plo(3):phi(3), nc)
  real(kind=dp_t), intent(in) :: b
  real(kind=dp_t), pointer :: cp(:,:,:,:)
  integer :: flev, fn, nlev
  fn = n + 1
  flev = lev+1
  nlev = size(mgts%coeffs)
  call mgt_verify_n("MGT_SET_CF", flev, fn, lo, hi)

  ! betax(i,j,k,1:nc) goes into coeffs(i,j,k,2:nc+1)
  cp => dataptr(mgts%coeffs(nlev), fn)
  cp(lo(1):hi(1)+1, lo(2):hi(2), lo(3):hi(3), 2:nc+1) = &
     b * cf(lo(1):hi(1)+1, lo(2):hi(2), lo(3):hi(3), 1:nc)

end subroutine mgt_set_cfbnx_3d

subroutine mgt_set_cfby_3d(lev, n, cf, b, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(3), hi(3), plo(3), phi(3)
  real(kind=dp_t), intent(in) :: cf(plo(1):phi(1), plo(2):phi(2), plo(3):phi(3))
  real(kind=dp_t), intent(in) :: b
  real(kind=dp_t), pointer :: cp(:,:,:,:)
  integer :: flev, fn, nlev
  fn = n + 1
  flev = lev+1
  nlev = size(mgts%coeffs)
  call mgt_verify_n("MGT_SET_CF", flev, fn, lo, hi)

  cp => dataptr(mgts%coeffs(nlev), fn)
  cp(lo(1):hi(1), lo(2):hi(2)+1, lo(3):hi(3), 3) = b * cf(lo(1):hi(1), lo(2):hi(2)+1, lo(3):hi(3))

end subroutine mgt_set_cfby_3d

subroutine mgt_set_cfbny_3d(lev, n, cf, b, plo, phi, lo, hi, nc)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, nc, lo(3), hi(3), plo(3), phi(3)
  real(kind=dp_t), intent(in) :: cf(plo(1):phi(1), plo(2):phi(2), plo(3):phi(3), nc)
  real(kind=dp_t), intent(in) :: b
  real(kind=dp_t), pointer :: cp(:,:,:,:)
  integer :: flev, fn, nlev
  fn = n + 1
  flev = lev+1
  nlev = size(mgts%coeffs)
  call mgt_verify_n("MGT_SET_CF", flev, fn, lo, hi)

  ! betay(i,j,k,1:nc) goes into coeffs(i,j,k,nc+2:2*nc+1)
  cp => dataptr(mgts%coeffs(nlev), fn)
  cp(lo(1):hi(1), lo(2):hi(2)+1, lo(3):hi(3), nc+2:2*nc+1) = &
      b * cf(lo(1):hi(1), lo(2):hi(2)+1, lo(3):hi(3), 1:nc)

end subroutine mgt_set_cfbny_3d

subroutine mgt_set_cfbz_3d(lev, n, cf, b, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(3), hi(3), plo(3), phi(3)
  real(kind=dp_t), intent(in) :: cf(plo(1):phi(1), plo(2):phi(2), plo(3):phi(3))
  real(kind=dp_t), intent(in) :: b
  real(kind=dp_t), pointer :: cp(:,:,:,:)
  integer :: flev, fn, nlev
  fn = n + 1
  flev = lev+1
  nlev = size(mgts%coeffs)
  call mgt_verify_n("MGT_SET_CF", flev, fn, lo, hi)

  cp => dataptr(mgts%coeffs(nlev), fn)
  cp(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)+1, 4) = b * cf(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)+1)

end subroutine mgt_set_cfbz_3d

subroutine mgt_set_cfbnz_3d(lev, n, cf, b, plo, phi, lo, hi, nc)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, nc, lo(3), hi(3), plo(3), phi(3)
  real(kind=dp_t), intent(in) :: cf(plo(1):phi(1), plo(2):phi(2), plo(3):phi(3), nc)
  real(kind=dp_t), intent(in) :: b
  real(kind=dp_t), pointer :: cp(:,:,:,:)
  integer :: flev, fn, nlev
  fn = n + 1
  flev = lev+1
  nlev = size(mgts%coeffs)
  call mgt_verify_n("MGT_SET_CF", flev, fn, lo, hi)

  ! betaz(i,j,k,1:nc) goes into coeffs(i,j,k,2*nc+2:3*nc+1)
  cp => dataptr(mgts%coeffs(nlev), fn)
  cp(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)+1, 2*nc+2:3*nc+1) = &
      b * cf(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)+1, 1:nc)

end subroutine mgt_set_cfbnz_3d

subroutine mgt_set_cfbx_3d_const(lev, n, lo, hi, coeff_value)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(3), hi(3)
  real(kind=dp_t), intent(in) :: coeff_value
  real(kind=dp_t), pointer :: cp(:,:,:,:)
  integer :: flev, fn, nlev
  fn = n + 1
  flev = lev+1
  nlev = size(mgts%coeffs)
  call mgt_verify_n("MGT_SET_CF", flev, fn, lo, hi)

  cp => dataptr(mgts%coeffs(nlev), fn)
  cp(lo(1):hi(1)+1, lo(2):hi(2), lo(3):hi(3), 2) = coeff_value

end subroutine mgt_set_cfbx_3d_const

subroutine mgt_set_cfby_3d_const(lev, n, lo, hi, coeff_value)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(3), hi(3)
  real(kind=dp_t), intent(in) :: coeff_value
  real(kind=dp_t), pointer :: cp(:,:,:,:)
  integer :: flev, fn, nlev
  fn = n + 1
  flev = lev+1
  nlev = size(mgts%coeffs)
  call mgt_verify_n("MGT_SET_CF", flev, fn, lo, hi)

  cp => dataptr(mgts%coeffs(nlev), fn)
  cp(lo(1):hi(1), lo(2):hi(2)+1, lo(3):hi(3), 3) = coeff_value

end subroutine mgt_set_cfby_3d_const

subroutine mgt_set_cfbz_3d_const(lev, n, lo, hi, coeff_value)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(3), hi(3)
  real(kind=dp_t), intent(in) :: coeff_value
  real(kind=dp_t), pointer :: cp(:,:,:,:)
  integer :: flev, fn, nlev
  fn = n + 1
  flev = lev+1
  nlev = size(mgts%coeffs)
  call mgt_verify_n("MGT_SET_CF", flev, fn, lo, hi)

  cp => dataptr(mgts%coeffs(nlev), fn)
  cp(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)+1, 4) = coeff_value

end subroutine mgt_set_cfbz_3d_const

subroutine mgt_set_uu_1d(lev, n, uu, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(1), hi(1), plo(1), phi(1)
  real(kind=dp_t), intent(in) :: uu(plo(1):phi(1))
  real(kind=dp_t), pointer :: up(:,:,:,:)
  integer :: flev,fn
  fn = n + 1
  flev = lev+1

  call mgt_verify_n("MGT_SET_UU_1D", flev, fn, lo, hi)

  up => dataptr(mgts%uu(flev), fn)
  up(lo(1)-1:hi(1)+1,1,1,1) = uu(lo(1)-1:hi(1)+1)

end subroutine mgt_set_uu_1d

subroutine mgt_set_uu_2d(lev, n, uu, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(2), hi(2), plo(2), phi(2)
  real(kind=dp_t), intent(in) :: uu(plo(1):phi(1), plo(2):phi(2))
  real(kind=dp_t), pointer :: up(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1

  call mgt_verify_n("MGT_SET_UU_2D", flev, fn, lo, hi)

  up => dataptr(mgts%uu(flev), fn)
  up(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1,1,1) = &
  uu(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1    )

end subroutine mgt_set_uu_2d

subroutine mgt_set_uu_3d(lev, n, uu, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(3), hi(3), plo(3), phi(3)
  real(kind=dp_t), intent(in) :: uu(plo(1):phi(1), plo(2):phi(2), plo(3):phi(3))
  real(kind=dp_t), pointer :: up(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1
  
  call mgt_verify_n("MGT_SET_UU_3D", flev, fn, lo, hi)

  up => dataptr(mgts%uu(flev), fn)
  up(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, lo(3)-1:hi(3)+1, 1) = &
  uu(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, lo(3)-1:hi(3)+1   )

end subroutine mgt_set_uu_3d

subroutine mgt_get_uu_1d(lev, n, uu, plo, phi, lo, hi, ng)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(1), hi(1), plo(1), phi(1), ng
  real(kind=dp_t), intent(inout) :: uu(plo(1):phi(1))
  real(kind=dp_t), pointer :: up(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1
  
  call mgt_verify_n("MGT_GET_UU_1D", flev, fn, lo, hi)

  up => dataptr(mgts%uu(flev), fn)
  uu(lo(1)-ng:hi(1)+ng) = up(lo(1)-ng:hi(1)+ng, 1,1,1)

end subroutine mgt_get_uu_1d
subroutine mgt_get_uu_2d(lev, n, uu, plo, phi, lo, hi, ng)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(2), hi(2), plo(2), phi(2), ng
  real(kind=dp_t), intent(inout) :: uu(plo(1):phi(1), plo(2):phi(2))
  real(kind=dp_t), pointer :: up(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1
  
  call mgt_verify_n("MGT_GET_UU_2D", flev, fn, lo, hi)

  up => dataptr(mgts%uu(flev), fn)
  uu(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng) = up(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,1,1)

end subroutine mgt_get_uu_2d
subroutine mgt_get_uu_3d(lev, n, uu, plo, phi, lo, hi ,ng)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(3), hi(3), plo(3), phi(3), ng
  real(kind=dp_t), intent(inout) :: uu(plo(1):phi(1), plo(2):phi(2), plo(3):phi(3))
  real(kind=dp_t), pointer :: up(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1
  
  call mgt_verify_n("MGT_GET_UU_3D", flev, fn, lo, hi)

  up => dataptr(mgts%uu(flev), fn)
  uu(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng) = &
  up(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng,1)

end subroutine mgt_get_uu_3d

subroutine mgt_get_gp_1d(lev, dir, n, gp, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: lev, dir, n, lo(1), hi(1), plo(1), phi(1)
  real(kind=dp_t), intent(inout) :: gp(plo(1):phi(1))
  real(kind=dp_t), pointer :: gpp(:,:,:,:)
  integer :: flev, fdir, fn
  fn = n + 1
  flev = lev+1
  fdir = dir+1

  gpp => dataptr(mgts%gp(flev,fdir), fn)
  gp(lo(1):hi(1)) = gpp(lo(1):hi(1),1,1,1)

end subroutine mgt_get_gp_1d

subroutine mgt_get_gp_2d(lev, dir, n, gp, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: lev, dir, n, lo(2), hi(2), plo(2), phi(2)
  real(kind=dp_t), intent(inout) :: gp(plo(1):phi(1), plo(2):phi(2))
  real(kind=dp_t), pointer :: gpp(:,:,:,:)
  integer :: flev, fdir, fn
  fn = n + 1
  flev = lev+1
  fdir = dir+1

  gpp => dataptr(mgts%gp(flev,fdir), fn)
  gp(lo(1):hi(1),lo(2):hi(2)) = gpp(lo(1):hi(1),lo(2):hi(2),1,1)

end subroutine mgt_get_gp_2d

subroutine mgt_get_gp_3d(lev, dir, n, gp, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: lev, dir, n, lo(3), hi(3), plo(3), phi(3)
  real(kind=dp_t), intent(inout) :: gp(plo(1):phi(1), plo(2):phi(2), plo(3):phi(3))
  real(kind=dp_t), pointer :: gpp(:,:,:,:)
  integer :: flev, fdir,fn
  fn = n + 1
  flev = lev+1
  fdir = dir+1

  gpp => dataptr(mgts%gp(flev,fdir), fn)
  gp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = &
         gpp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1)

end subroutine mgt_get_gp_3d

subroutine mgt_get_res_1d(lev, n, res, rlo, rhi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(1), hi(1), rlo(1), rhi(1)
  real(kind=dp_t), intent(inout) :: res(rlo(1):rhi(1))
  real(kind=dp_t), pointer :: rp(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1

  call mgt_verify_n("MGT_GET_RES", flev, fn, lo, hi)

  rp => dataptr(mgts%res(flev), fn)
  res(lo(1):hi(1)) = rp(lo(1):hi(1),1,1,1)

end subroutine mgt_get_res_1d

subroutine mgt_get_res_2d(lev, n, res, rlo, rhi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(2), hi(2), rlo(2), rhi(2)
  real(kind=dp_t), intent(inout) :: res(rlo(1):rhi(1), rlo(2):rhi(2))
  real(kind=dp_t), pointer :: rp(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1

  call mgt_verify_n("MGT_GET_RES", flev, fn, lo, hi)

  rp => dataptr(mgts%res(flev), fn)
  res(lo(1):hi(1),lo(2):hi(2)) = rp(lo(1):hi(1),lo(2):hi(2),1,1)

end subroutine mgt_get_res_2d
subroutine mgt_get_res_3d(lev, n, res, rlo, rhi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(3), hi(3), rlo(3), rhi(3)
  real(kind=dp_t), intent(inout) :: res(rlo(1):rhi(1), rlo(2):rhi(2), rlo(3):rhi(3))
  real(kind=dp_t), pointer :: rp(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1

  call mgt_verify_n("MGT_GET_RES", flev, fn, lo, hi)

  rp => dataptr(mgts%res(flev), fn)
  res(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = rp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1)

end subroutine mgt_get_res_3d


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
     call destroy(mgts%res(i))
     call destroy(mgts%uu(i))
  end do

  deallocate(mgts%rr)
  deallocate(mgts%rh)
  deallocate(mgts%res)
  deallocate(mgts%pd)
  deallocate(mgts%uu)
  deallocate(mgts%gp)
  deallocate(mgts%mgt)
  deallocate(mgts%bc)

  call destroy(mgts%mla)
  mgts%dim = 0
  mgts%final = .false.

end subroutine mgt_dealloc

subroutine mgt_solve(tol,abs_tol,needgradphi,final_resnorm)
  use cpp_mg_module
  use ml_cc_module
  use fabio_module
  implicit none
  real(kind=dp_t), intent(in   ) :: tol, abs_tol
  integer        , intent(in   ), optional :: needgradphi
  real(kind=dp_t), intent(  out), optional :: final_resnorm

  integer :: do_diagnostics

  call mgt_verify("MGT_SOLVE")
  if ( .not. mgts%final ) then
     call bl_error("MGT_SOLVE: MGT not finalized")
  end if

  do_diagnostics = 0
  if (present(needgradphi)) then
    if (needgradphi .eq. 1) then
      call ml_cc(mgts%mla, mgts%mgt, &
           mgts%rh, mgts%uu, &
           mgts%mla%mask, mgts%rr, &
           do_diagnostics, tol, &
           abs_eps_in = abs_tol, &
           need_grad_phi_in = .true.,&
           final_resnorm = final_resnorm)
    else
      call ml_cc(mgts%mla, mgts%mgt, &
           mgts%rh, mgts%uu, &
           mgts%mla%mask, mgts%rr, &
           do_diagnostics, tol, &
           abs_eps_in = abs_tol, &
           final_resnorm = final_resnorm )
    end if
  else
    call ml_cc(mgts%mla, mgts%mgt, &
         mgts%rh, mgts%uu, &
         mgts%mla%mask, mgts%rr, &
         do_diagnostics, tol,&
         abs_eps_in = abs_tol, &
         final_resnorm = final_resnorm )
  end if

end subroutine mgt_solve

subroutine mgt_applyop()
  use cpp_mg_module
  use ml_cc_module
  use fabio_module
  implicit none

  call mgt_verify("MGT_APPLYOP")
  if ( .not. mgts%final ) then
     call bl_error("MGT_APPLYOP: MGT not finalized")
  end if

  call ml_cc_applyop(mgts%mla, mgts%mgt, &
                     mgts%res, mgts%uu, &
                     mgts%rr)

end subroutine mgt_applyop

subroutine mgt_compute_residual()
  use cpp_mg_module
  use ml_cc_module
  use fabio_module
  implicit none

  call mgt_verify("MGT_COMPUTE_RESIDUAL")
  if ( .not. mgts%final ) then
     call bl_error("MGT_COMPUTE_RESIDUAL: MGT not finalized")
  end if

  call ml_resid(mgts%mla, mgts%mgt, &
                mgts%rh, mgts%res, mgts%uu, &
                mgts%rr)

end subroutine mgt_compute_residual

subroutine mgt_compute_flux(lev)
  use cpp_mg_module
  use ml_cc_module
  use ml_util_module
  use fabio_module
  implicit none
 
  integer, intent(in) :: lev
  integer             :: mglev,flev
  integer             :: dir

  flev = lev+1

  call mgt_verify("MGT_COMPUTE_FLUX")
  if ( .not. mgts%final ) then
     call bl_error("MGT_COMPUTE_FLUX: MGT not finalized")
  end if

  do dir = 1, mgts%dim
     call build(mgts%gp(flev,dir), mgts%mla%la(flev), nc = 1, ng = 1)
  end do

  mglev = mgts%mgt(flev)%nlevels
  call ml_fill_all_fluxes(mgts%mgt(flev)%ss(mglev), mgts%gp(flev,:), &
                          mgts%uu(flev), mgts%mgt(flev)%mm(mglev))

end subroutine mgt_compute_flux

subroutine mgt_delete_flux(lev)

  use cpp_mg_module
! use ml_cc_module
! use ml_util_module
  implicit none

  integer, intent(in) :: lev
  integer             :: flev
  integer             :: dir

  flev = lev+1

  do dir = 1, mgts%dim
     call destroy(mgts%gp(flev,dir))
  end do

end subroutine mgt_delete_flux

subroutine mgt_set_defaults(nu_1,nu_2,nu_b,nu_f,gamma,omega,max_iter,bottom_max_iter, &
                            bottom_solver,bottom_solver_eps, &
                            verbose,cg_verbose,max_nlevel,min_width,cycle_type,smoother)
  use cpp_mg_module
  implicit none
  integer   , intent(in) :: nu_1,nu_2,nu_b,nu_f,gamma,max_iter,bottom_max_iter,bottom_solver
  integer   , intent(in) :: verbose, cg_verbose, max_nlevel, min_width, cycle_type, smoother
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
  mgts%cycle_type      = cycle_type
  mgts%bottom_max_iter = bottom_max_iter
  mgts%bottom_solver   = bottom_solver
  mgts%bottom_solver_eps = bottom_solver_eps
  mgts%max_nlevel      = max_nlevel
  mgts%min_width       = min_width

end subroutine mgt_set_defaults

subroutine mgt_get_defaults(nu_1,nu_2,nu_b,nu_f,gamma,omega,max_iter,bottom_max_iter, &
                            bottom_solver, &
                            verbose,cg_verbose,max_nlevel,min_width,cycle_type,smoother)
  use cpp_mg_module
  implicit none
  integer   , intent(out) :: nu_1,nu_2,nu_b,nu_f,gamma,max_iter,bottom_max_iter,bottom_solver
  integer   , intent(out) :: verbose, cg_verbose, max_nlevel, min_width, cycle_type, smoother
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
  cycle_type = mgts%mg_tower_default%cycle_type

  bottom_max_iter = mgts%mg_tower_default%bottom_max_iter
  bottom_solver    = mgts%mg_tower_default%bottom_solver
  max_nlevel      = mgts%mg_tower_default%max_nlevel
  min_width       = mgts%mg_tower_default%min_width

end subroutine mgt_get_defaults
