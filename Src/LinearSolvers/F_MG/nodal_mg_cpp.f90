module nodal_cpp_mg_module

  use mg_module
  use multifab_module
  use ml_layout_module
  use bl_constants_module

  implicit none

  type mg_server

     logical         :: final = .false.
     integer         :: dim  = 0
     integer         :: nlevel
     integer         :: nu1, nu2, nuf, nub
     integer         :: max_iter
     integer         :: max_nlevel
     integer         :: min_width
     integer         :: smoother
     integer         :: cycle_type
     integer         :: verbose
     integer         :: cg_verbose
     integer         :: bottom_max_iter
     integer         :: bottom_solver
     integer         :: stencil_type
     real(dp_t)      :: bottom_solver_eps
     real(dp_t)      :: eps
     type(ml_layout) :: mla
     type(mg_tower)  :: mg_tower_default
     type(mg_tower), pointer :: mgt(:) => Null()
     type(box), pointer :: pd(:) => Null()
     logical, pointer :: nodal(:) => Null()
     integer, pointer :: bc(:,:) => Null()
     integer, pointer :: rr(:,:)
     type(multifab) , pointer ::           rh(:) => Null()
     type(multifab) , pointer ::           uu(:) => Null()
     type(multifab) , pointer ::          vel(:) => Null()
     type(multifab) , pointer ::   amr_coeffs(:) => Null()
     type(lmultifab), pointer ::    fine_mask(:) => Null()
     type(multifab) , pointer ::     sync_res(:) => Null()
     type(multifab) , pointer ::     sync_msk(:) => Null()
     type(multifab) , pointer ::         vold(:) => Null()
     type(multifab) , pointer ::         rhcc(:) => Null()
  end type mg_server

  type(mg_server), save :: mgts

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

end module nodal_cpp_mg_module

subroutine mgt_nodal_alloc(dm, nlevel, stencil_type)
  use nodal_cpp_mg_module
  implicit none
  integer, intent(in) :: dm, nlevel, stencil_type

  if ( mgts%dim == 0 ) then
     mgts%dim = dm
     mgts%nlevel = nlevel
     allocate(mgts%nodal(dm))
     mgts%nodal = .true.
  end if

  mgts%stencil_type = stencil_type

  allocate(mgts%rr(nlevel-1,dm))
  allocate(mgts%rh(nlevel))
  allocate(mgts%vel(nlevel))
  allocate(mgts%pd(nlevel))
  allocate(mgts%uu(nlevel))
  allocate(mgts%mgt(nlevel))
  allocate(mgts%amr_coeffs(nlevel))
  allocate(mgts%fine_mask(nlevel))

  call build(mgts%mla, nlevel, dm)

end subroutine mgt_nodal_alloc

subroutine mgt_set_nodal_level(lev, nb, dm, lo, hi, pd_lo, pd_hi, pm, pmap)
  use nodal_cpp_mg_module
  implicit none
  integer, intent(in) :: lev, nb, dm
  integer, intent(in) :: lo(nb,dm), hi(nb,dm), pd_lo(dm), pd_hi(dm), pm(dm), pmap(nb+1)

  integer   :: i
  logical   :: pmask(dm)
  integer   :: flev

  type(box), allocatable :: bxs(:)

  allocate(bxs(nb))

  flev = lev + 1
  call mgt_verify_lev("MGT_SET_NODAL_LEVEL", flev)

  pmask = (pm /= 0)

  if ( dm /= mgts%dim ) then
     call bl_error("MGT_SET_NODAL_LEVEL: Input DIM doesn't match internal DIM")
  end if
  call build(mgts%pd(flev), pd_lo(1:dm), pd_hi(1:dm))
  do i = 1, nb
     bxs(i) = make_box(lo(i,:), hi(i,:))
  end do
  call boxarray_build_v(mgts%mla%mba%bas(flev), bxs)
  call layout_build_ba(mgts%mla%la(flev),  &
                       mgts%mla%mba%bas(flev), &
                       mgts%pd(flev), pmask = pmask, &
                       mapping = LA_EXPLICIT, explicit_mapping = pmap(1:nb))

end subroutine mgt_set_nodal_level

subroutine mgt_nodal_finalize(dx,bc)

  use nodal_cpp_mg_module
  use stencil_types_module

  implicit none
  real(dp_t), intent(in) :: dx(mgts%nlevel,mgts%dim)
  integer   , intent(in) :: bc(2,mgts%dim)
  integer :: dm, i, nlev, n
  logical, allocatable :: nodal(:)

  integer :: max_nlevel_in
  integer :: bottom_solver_in
  integer :: bottom_max_iter_in

  call mgt_verify("MGT_FINALIZE")

  dm = mgts%dim

  allocate(mgts%bc(dm,2))
  mgts%bc = transpose(bc)

  allocate(nodal(1:dm))
  nodal = mgts%nodal

  nlev = mgts%nlevel

  do i = 1, nlev-1
     mgts%mla%mba%rr(i,:) = (upb(mgts%pd(i+1)) -  lwb(mgts%pd(i+1)) + 1) / &
                            (upb(mgts%pd(i  )) -  lwb(mgts%pd(i  )) + 1) 
             mgts%rr(i,:) = mgts%mla%mba%rr(i,:)
  end do

  do i = 1, nlev
     call multifab_build(mgts%uu(i) , mgts%mla%la(i), nc = 1, ng = 1, nodal = nodal)
     call multifab_build(mgts%rh(i) , mgts%mla%la(i), nc = 1, ng = 1, nodal = nodal)
     call multifab_build(mgts%vel(i), mgts%mla%la(i), nc =dm, ng = 1)

     call setval(mgts%uu(i),ZERO,all=.true.)
     call setval(mgts%rh(i),ZERO,all=.true.)
     call setval(mgts%vel(i),ZERO,all=.true.)
  end do

  if (mgts%stencil_type .eq. ND_DENSE_STENCIL) then
     if ( parallel_ioprocessor() .and. mgts%verbose > 0 ) &
         print *,'Using dense stencil in nodal solver ...'
  else if (mgts%stencil_type .eq. ND_CROSS_STENCIL) then
     if ( parallel_ioprocessor() .and. mgts%verbose > 0 ) &
         print *,'Using cross stencil in nodal solver ...'
  else if (mgts%stencil_type .eq. ND_VATER_STENCIL) then
     if ( parallel_ioprocessor() .and. mgts%verbose > 0 ) &
         print *,'Using Vater dense stencil in nodal solver ...'
  else
     if ( parallel_ioprocessor()) &
         print *,'Dont know this stencil type ',mgts%stencil_type
     call bl_error("MGT_FINALIZE: stuck")
  end if

  do n = nlev, 1, -1
     if ( n == 1 ) then
        max_nlevel_in = mgts%max_nlevel
        bottom_solver_in = mgts%bottom_solver
        bottom_max_iter_in = mgts%bottom_max_iter
     else
        if ( all(mgts%rr(n-1,:) == 2) ) then
           max_nlevel_in = 1
        else if ( all(mgts%rr(n-1,:) == 4) ) then
           max_nlevel_in = 2
        else
           call bl_error("MGT_FINALIZE: confused about ref_ratio")
        end if
        bottom_solver_in = 0
        bottom_max_iter_in = mgts%nu1
     end if

     call mg_tower_build(mgts%mgt(n), mgts%mla%la(n), mgts%pd(n), mgts%bc, mgts%stencil_type, &
          dh                = dx(n,:), &
          smoother          = mgts%smoother, &
          nu1               = mgts%nu1, &
          nu2               = mgts%nu2, &
          nuf               = mgts%nuf, &
          nub               = mgts%nub, &
          cycle_type        = mgts%cycle_type, &
          bottom_solver     = bottom_solver_in, &
          bottom_max_iter   = bottom_max_iter_in, &
          bottom_solver_eps = mgts%bottom_solver_eps, &
          max_iter          = mgts%max_iter, &
          max_nlevel        = max_nlevel_in, &
          min_width         = mgts%min_width, &
          verbose           = mgts%verbose, &
          cg_verbose        = mgts%cg_verbose, &
          nodal             = nodal)

  end do

  if ( parallel_IOProcessor() .and. mgts%verbose > 0) then
    call mg_tower_print(mgts%mgt(nlev))
  end if

  deallocate(nodal)
  deallocate(mgts%bc)

end subroutine mgt_nodal_finalize

subroutine mgt_init_nodal_coeffs_lev(lev,val)
  use nodal_cpp_mg_module
  implicit none
  integer   , intent(in)           :: lev
  real(dp_t), intent(in), optional :: val

  integer :: nlev
  integer :: flev
  flev = lev + 1
  call mgt_verify_lev("MGT_INIT_NODAL_COEFFS_LEV", flev)

  nlev = mgts%mgt(flev)%nlevels

  ! These only exist at amr levels, not the lower multigrid levels
  call multifab_build( mgts%amr_coeffs(flev), mgts%mgt(flev)%ss(nlev)%la, 1, 1)
  call multifab_setval(mgts%amr_coeffs(flev), ZERO, all=.true.)

end subroutine mgt_init_nodal_coeffs_lev

subroutine mgt_init_const_nodal_coeffs_lev(lev,val)
  use nodal_cpp_mg_module
  implicit none
  integer   , intent(in) :: lev
  real(dp_t), intent(in) :: val

  integer :: nlev
  integer :: flev
  flev = lev + 1
  call mgt_verify_lev("MGT_INIT_CONST_NODAL_COEFFS_LEV", flev)

  nlev = mgts%mgt(flev)%nlevels

  ! These only exist at amr levels, not the lower multigrid levels
  call multifab_build( mgts%amr_coeffs(flev), mgts%mgt(flev)%ss(nlev)%la, 1, 1)
  call multifab_setval(mgts%amr_coeffs(flev), val, all=.true.)

end subroutine mgt_init_const_nodal_coeffs_lev

subroutine mgt_finalize_nodal_stencil_lev(lev)
  use nodal_cpp_mg_module
  use nodal_stencil_fill_module
  implicit none
  integer, intent(in) :: lev
  integer :: nlev
  integer :: flev
  type(multifab), allocatable :: cell_coeffs(:)

  flev = lev + 1

  call mgt_verify_lev("MGT_FINALIZE_NODAL_STENCIL_LEV", flev)
  
  nlev = mgts%mgt(flev)%nlevels

  allocate(cell_coeffs(nlev))
  call multifab_build_copy(cell_coeffs(nlev), mgts%amr_coeffs(flev))

  call multifab_fill_boundary(cell_coeffs(nlev))

  call stencil_fill_nodal_all_mglevels(mgts%mgt(flev), cell_coeffs)

  call multifab_destroy(cell_coeffs(nlev))
  deallocate(cell_coeffs)

end subroutine mgt_finalize_nodal_stencil_lev

subroutine mgt_divu(lo_inflow, hi_inflow)

  use nodal_divu_module
  use nodal_cpp_mg_module

  integer, intent(in) :: lo_inflow(3), hi_inflow(3)

  integer    :: n
  real(dp_t) :: r, rhmax

  call divu(mgts%nlevel,mgts%mgt,mgts%vel,mgts%rh,mgts%rr,mgts%nodal,  &
       lo_inflow, hi_inflow)

  if (mgts%verbose > 0) then
     rhmax = norm_inf(mgts%rh(mgts%nlevel),local=.true.)
     do n = mgts%nlevel-1, 1, -1
       r = norm_inf(mgts%rh(n),mgts%fine_mask(n),local=.true.)
       rhmax = max(r,rhmax) 
     end do
     call parallel_reduce(r, rhmax, MPI_MAX, proc = parallel_IOProcessorNode())
     rhmax = r
     if (parallel_IOProcessor()) then
        print *,'F90: Source norm is ',rhmax
     endif
  end if

end subroutine mgt_divu

subroutine mgt_newu()
  use nodal_cpp_mg_module
  use nodal_newu_module

  integer :: n

  do n = 1,mgts%nlevel
     call mkunew(mgts%vel(n),mgts%uu(n),mgts%amr_coeffs(n), &
                 mgts%mgt(n)%dh(:,mgts%mgt(n)%nlevels), &
                 mgts%vel(n)%ng)
  end do
       
end subroutine mgt_newu

subroutine mgt_finalize_nodal_stencil()
   use nodal_mask_module
   use nodal_cpp_mg_module
   implicit none
   integer :: i
   call mgt_verify("MGT_FINALIZE_NODAL_STENCIL")

   do i = 1, mgts%nlevel-1
     call lmultifab_build(mgts%fine_mask(i), mgts%mla%la(i), nc = 1, ng = 0, nodal = mgts%nodal)
     call lmultifab_setval(mgts%fine_mask(i), val = .TRUE., all = .TRUE.)
     call create_nodal_mask(i,mgts%fine_mask(i), &
                            mgts%mgt(i  )%mm(mgts%mgt(i  )%nlevels), &
                            mgts%mgt(i+1)%mm(mgts%mgt(i+1)%nlevels), &
                            mgts%mla)
   end do

   mgts%final = .true.
end subroutine mgt_finalize_nodal_stencil

subroutine mgt_set_vel_1d(lev, n, vel_in, plo, phi, lo, hi, nv, iv)
  use nodal_cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(1), hi(1), plo(1), phi(1), nv, iv
  real(kind=dp_t), intent(in) :: vel_in(plo(1):phi(1), nv)
  real(kind=dp_t), pointer :: vp(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1
  vp => dataptr(mgts%vel(flev), fn)
  vp(lo(1):hi(1),1,1,1) = vel_in(lo(1):hi(1),iv+1)
end subroutine mgt_set_vel_1d

subroutine mgt_get_vel_1d(lev, n, vel_out, plo, phi, lo, hi, nv, iv)
  use nodal_cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(1), hi(1), plo(1), phi(1), nv, iv
  real(kind=dp_t), intent(inout) :: vel_out(plo(1):phi(1), nv)
  real(kind=dp_t), pointer :: vp(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1
  vp => dataptr(mgts%vel(flev), fn)
  vel_out(lo(1):hi(1),iv+1) = vp(lo(1):hi(1),1,1,1)
end subroutine mgt_get_vel_1d

subroutine mgt_set_vel_2d(lev, n, vel_in, plo, phi, lo, hi, nv, iv)
  use nodal_cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(2), hi(2), plo(2), phi(2), nv, iv
  real(kind=dp_t), intent(in) :: vel_in(plo(1):phi(1), plo(2):phi(2), nv)
  real(kind=dp_t), pointer :: vp(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1
  vp => dataptr(mgts%vel(flev), fn)
  vp         (lo(1):hi(1), lo(2):hi(2), 1, 1:2) =  &
       vel_in(lo(1):hi(1), lo(2):hi(2), iv+1:iv+2)
end subroutine mgt_set_vel_2d

subroutine mgt_get_vel_2d(lev, n, vel_out, plo, phi, lo, hi, nv, iv)
  use nodal_cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(2), hi(2), plo(2), phi(2), nv, iv
  real(kind=dp_t), intent(inout) :: vel_out(plo(1):phi(1), plo(2):phi(2), nv)
  real(kind=dp_t), pointer :: vp(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1
  vp => dataptr(mgts%vel(flev), fn)
  vel_out(lo(1):hi(1), lo(2):hi(2), iv+1:iv+2) = vp(lo(1):hi(1), lo(2):hi(2), 1, 1:2)
end subroutine mgt_get_vel_2d

subroutine mgt_set_vel_3d(lev, n, vel_in, plo, phi, lo, hi, nv, iv)
  use nodal_cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(3), hi(3), plo(3), phi(3), nv, iv
  real(kind=dp_t), intent(in) :: vel_in(plo(1):phi(1), plo(2):phi(2), plo(3):phi(3), nv)
  real(kind=dp_t), pointer :: vp(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1  
  vp => dataptr(mgts%vel(flev), fn)
  vp         (lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 1:3) =  &
       vel_in(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), iv+1:iv+3)

end subroutine mgt_set_vel_3d

subroutine mgt_get_vel_3d(lev, n, vel_out, plo, phi, lo, hi, nv, iv)
  use nodal_cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(3), hi(3), plo(3), phi(3), nv, iv
  real(kind=dp_t), intent(inout) :: vel_out(plo(1):phi(1), plo(2):phi(2), plo(3):phi(3), nv)
  real(kind=dp_t), pointer :: vp(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1
  vp => dataptr(mgts%vel(flev), fn)
  vel_out(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), iv+1:iv+3) = &
       vp(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 1:3)
end subroutine mgt_get_vel_3d

subroutine mgt_set_cfs_1d(lev, n, cf, plo, phi, lo, hi)
  use nodal_cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(1), hi(1), plo(1), phi(1)
  real(kind=dp_t), intent(in) :: cf(plo(1):phi(1))
  real(kind=dp_t), pointer :: acp(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1
  acp => dataptr(mgts%amr_coeffs(flev), fn)
  acp(lo(1):hi(1), 1, 1, 1) = cf(lo(1):hi(1))
end subroutine mgt_set_cfs_1d

subroutine mgt_set_cfs_2d(lev, n, cf, plo, phi, lo, hi)
  use nodal_cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(2), hi(2), plo(2), phi(2)
  real(kind=dp_t), intent(in) :: cf(plo(1):phi(1), plo(2):phi(2))
  real(kind=dp_t), pointer :: acp(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1
  acp => dataptr(mgts%amr_coeffs(flev), fn)
  acp(lo(1):hi(1), lo(2):hi(2), 1, 1) = cf(lo(1):hi(1), lo(2):hi(2))
end subroutine mgt_set_cfs_2d

subroutine mgt_set_cfs_3d(lev, n, cf, plo, phi, lo, hi)
  use nodal_cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(3), hi(3), plo(3), phi(3)
  real(kind=dp_t), intent(in) :: cf(plo(1):phi(1), plo(2):phi(2), plo(3):phi(3))
  real(kind=dp_t), pointer :: acp(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1
  acp => dataptr(mgts%amr_coeffs(flev), fn)
  acp(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 1) = cf(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))
end subroutine mgt_set_cfs_3d

subroutine mgt_set_pr_1d(lev, n, uu, plo, phi, lo, hi, np, ip)
  use nodal_cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(1), hi(1), plo(1), phi(1), np, ip
  real(kind=dp_t), intent(in) :: uu(plo(1):phi(1), np)
  real(kind=dp_t), pointer :: up(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1
  up => dataptr(mgts%uu(flev), fn)
  up(lo(1):hi(1),1,1,1) = uu(lo(1):hi(1), ip+1)
end subroutine mgt_set_pr_1d

subroutine mgt_set_pr_2d(lev, n, uu, plo, phi, lo, hi, np, ip)
  use nodal_cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(2), hi(2), plo(2), phi(2), np, ip
  real(kind=dp_t), intent(in) :: uu(plo(1):phi(1), plo(2):phi(2), np)
  real(kind=dp_t), pointer :: up(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1
  up => dataptr(mgts%uu(flev), fn)
  up     (lo(1):hi(1), lo(2):hi(2), 1,1) = &
       uu(lo(1):hi(1), lo(2):hi(2), ip+1)
end subroutine mgt_set_pr_2d

subroutine mgt_set_pr_3d(lev, n, uu, plo, phi, lo, hi, np, ip)
  use nodal_cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(3), hi(3), plo(3), phi(3), np, ip
  real(kind=dp_t), intent(in) :: uu(plo(1):phi(1), plo(2):phi(2), plo(3):phi(3), np)
  real(kind=dp_t), pointer :: up(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1
  up => dataptr(mgts%uu(flev), fn)
  up     (lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 1) = &
       uu(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), ip+1)
end subroutine mgt_set_pr_3d

subroutine mgt_get_pr_1d(lev, n, uu, plo, phi, lo, hi, np, ip)
  use nodal_cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(1), hi(1), plo(1), phi(1), np, ip
  real(kind=dp_t), intent(inout) :: uu(plo(1):phi(1), np)
  real(kind=dp_t), pointer :: up(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1
  up => dataptr(mgts%uu(flev), fn)
  uu(lo(1):hi(1), ip+1) = up(lo(1):hi(1),1,1,1)
end subroutine mgt_get_pr_1d

subroutine mgt_get_pr_2d(lev, n, uu, plo, phi, lo, hi, np, ip)
  use nodal_cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(2), hi(2), plo(2), phi(2), np, ip
  real(kind=dp_t), intent(inout) :: uu(plo(1):phi(1), plo(2):phi(2), np)
  real(kind=dp_t), pointer :: up(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1
  up => dataptr(mgts%uu(flev), fn)
  uu(lo(1):hi(1), lo(2):hi(2), ip+1) = up(lo(1):hi(1), lo(2):hi(2),1,1)
end subroutine mgt_get_pr_2d

subroutine mgt_get_pr_3d(lev, n, uu, plo, phi, lo, hi, np, ip)
  use nodal_cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(3), hi(3), plo(3), phi(3), np, ip
  real(kind=dp_t), intent(inout) :: uu(plo(1):phi(1), plo(2):phi(2), plo(3):phi(3), np)
  real(kind=dp_t), pointer :: up(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1
  up => dataptr(mgts%uu(flev), fn)
  uu     (lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), ip+1) =  &
       up(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 1)
end subroutine mgt_get_pr_3d

subroutine mgt_add_rh_nodal_1d(lev, n, rh_in, plo, phi, lo, hi, rhmax)
  use nodal_cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(1), hi(1), plo(1), phi(1)
  real(kind=dp_t), intent(inout) :: rhmax
  real(kind=dp_t), intent(in   ) :: rh_in(plo(1):phi(1))
  real(kind=dp_t), pointer       :: rp(:,:,:,:)
  integer        , pointer       :: mp(:,:,:,:)
  integer                        :: flev, fn, i
  fn = n + 1
  flev = lev+1

  rp => dataptr(mgts%rh(flev), fn)
  mp => dataptr(mgts%mgt(flev)%mm(mgts%mgt(flev)%nlevels), fn)

  ! Only add in the nodal RHS if it is on a non-Dirichlet node
  do i = lo(1),hi(1)
      if (.not. bc_dirichlet(mp(i,1,1,1),1,0)) &
          rp(i,1,1,1) = rp(i,1,1,1) + rh_in(i)
      rhmax = max(rhmax, abs(rp(i,1,1,1)))
  end do

end subroutine mgt_add_rh_nodal_1d

subroutine mgt_add_rh_nodal_2d(lev, n, rh_in, plo, phi, lo, hi, rhmax)
  use nodal_cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(2), hi(2), plo(2), phi(2)
  real(kind=dp_t), intent(inout) :: rhmax
  real(kind=dp_t), intent(in   ) :: rh_in(plo(1):phi(1), plo(2):phi(2))
  real(kind=dp_t), pointer       :: rp(:,:,:,:)
  integer        , pointer       :: mp(:,:,:,:)
  integer                        :: flev, fn, i, j
  fn = n + 1
  flev = lev+1

  rp => dataptr(mgts%rh(flev), fn)
  mp => dataptr(mgts%mgt(flev)%mm(mgts%mgt(flev)%nlevels), fn)

  ! Only add in the nodal RHS if it is on a non-Dirichlet node
  do j = lo(2),hi(2)
  do i = lo(1),hi(1)
      if (.not. bc_dirichlet(mp(i,j,1,1),1,0)) &
          rp(i,j,1,1) = rp(i,j,1,1) + rh_in(i,j)
      rhmax = max(rhmax, abs(rp(i,j,1,1)))
  end do
  end do

end subroutine mgt_add_rh_nodal_2d

subroutine mgt_add_rh_nodal_3d(lev, n, rh_in, plo, phi, lo, hi, rhmax)
  use nodal_cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(3), hi(3), plo(3), phi(3)
  real(kind=dp_t), intent(inout) :: rhmax
  real(kind=dp_t), intent(in   ) :: rh_in(plo(1):phi(1),plo(2):phi(2),plo(3):phi(3))
  real(kind=dp_t), pointer       :: rp(:,:,:,:)
  integer        , pointer       :: mp(:,:,:,:)
  integer                        :: flev, fn, i, j, k
  fn = n + 1
  flev = lev+1

  rp => dataptr(mgts%rh(flev), fn)
  mp => dataptr(mgts%mgt(flev)%mm(mgts%mgt(flev)%nlevels), fn)

  ! Only add in the nodal RHS if it is on a non-Dirichlet node
  do k = lo(3),hi(3)
  do j = lo(2),hi(2)
  do i = lo(1),hi(1)
      if (.not. bc_dirichlet(mp(i,j,k,1),1,0)) &
          rp(i,j,k,1) = rp(i,j,k,1) + rh_in(i,j,k)
      rhmax = max(rhmax, abs(rp(i,j,k,1)))
  end do
  end do
  end do

end subroutine mgt_add_rh_nodal_3d

subroutine mgt_nodal_dealloc()
  use nodal_cpp_mg_module
  implicit none
  integer :: i
  
  call mgt_verify("MGT_NODAL_DEALLOC")
  if ( .not. mgts%final ) then
     call bl_error("MGT_NODAL_DEALLOC: MGT not finalized")
  end if

  do i = 1, mgts%nlevel
    call mg_tower_destroy(mgts%mgt(i))
  end do

  deallocate(mgts%nodal)

  do i = 1,mgts%nlevel
     call multifab_destroy(mgts%rh(i))
     call multifab_destroy(mgts%uu(i))
     call multifab_destroy(mgts%vel(i))
     call multifab_destroy(mgts%amr_coeffs(i))
  end do
  do i = 1,mgts%nlevel-1
     call lmultifab_destroy(mgts%fine_mask(i))
  end do
  call destroy(mgts%mla)
  mgts%dim = 0
  mgts%final = .false.

  deallocate(mgts%rr)
  deallocate(mgts%rh)
  deallocate(mgts%vel)
  deallocate(mgts%pd)
  deallocate(mgts%uu)
  deallocate(mgts%mgt)
  deallocate(mgts%amr_coeffs)
  deallocate(mgts%fine_mask)

end subroutine mgt_nodal_dealloc

subroutine mgt_nodal_solve(tol, abs_tol)
  use nodal_cpp_mg_module
  use ml_nd_module

  implicit none
  real(kind=dp_t), intent(in) :: tol, abs_tol

  integer :: do_diagnostics

  call mgt_verify("MGT_SOLVE")
  if ( .not. mgts%final ) then
     call bl_error("MGT_SOLVE: MGT not finalized")
  end if

  if (mgts%verbose >= 4) then 
     do_diagnostics = 1
  else
     do_diagnostics = 0
  endif

  mgts%mgt%eps     = tol
  mgts%mgt%abs_eps = abs_tol

  call ml_nd(mgts%mla, mgts%mgt, &
       mgts%rh, mgts%uu, &
       mgts%fine_mask, &
       do_diagnostics)

end subroutine mgt_nodal_solve

subroutine mgt_set_nodal_defaults(nu_1,nu_2,nu_b,nu_f,max_iter,bottom_max_iter, &
                                  bottom_solver,bottom_solver_eps, &
                                  verbose,cg_verbose,max_nlevel,min_width,cycle_type,smoother,stencil_type)
  use nodal_cpp_mg_module
  implicit none
  integer   , intent(in) :: nu_1,nu_2,nu_b,nu_f,max_iter,bottom_max_iter,bottom_solver
  integer   , intent(in) :: verbose, cg_verbose, max_nlevel, min_width, cycle_type, smoother, stencil_type
  real(dp_t), intent(in) :: bottom_solver_eps

  call mgt_not_final("MGT_SET_NODAL_DEFAULTS")

  mgts%nu1             = nu_1
  mgts%nu2             = nu_2
  mgts%nuf             = nu_f
  mgts%nub             = nu_b
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

  mgts%stencil_type    = stencil_type

end subroutine mgt_set_nodal_defaults

subroutine mgt_get_nodal_defaults(nu_1,nu_2,nu_b,nu_f,max_iter,bottom_max_iter, &
                                  bottom_solver, &
                                  verbose,cg_verbose,max_nlevel,min_width,cycle_type,smoother)
  use nodal_cpp_mg_module
  implicit none
  integer   , intent(out) :: nu_1,nu_2,nu_b,nu_f,max_iter,bottom_max_iter,bottom_solver
  integer   , intent(out) :: verbose, cg_verbose, max_nlevel, min_width, cycle_type, smoother

  nu_1       = mgts%mg_tower_default%nu1
  nu_2       = mgts%mg_tower_default%nu2
  nu_f       = mgts%mg_tower_default%nuf
  nu_b       = mgts%mg_tower_default%nub
  max_iter   = mgts%mg_tower_default%max_iter
  verbose    = mgts%mg_tower_default%verbose
  cg_verbose = mgts%mg_tower_default%cg_verbose
  smoother   = mgts%mg_tower_default%smoother
  cycle_type = mgts%mg_tower_default%cycle_type

  bottom_max_iter = mgts%mg_tower_default%bottom_max_iter
  bottom_solver    = mgts%mg_tower_default%bottom_solver
  max_nlevel      = mgts%mg_tower_default%max_nlevel
  min_width       = mgts%mg_tower_default%min_width

end subroutine mgt_get_nodal_defaults


subroutine mgt_alloc_rhcc_nodal()
  use nodal_cpp_mg_module
  implicit none
  integer :: i

  allocate(mgts%rhcc(mgts%nlevel))

  do i = 1, mgts%nlevel
     call multifab_build(mgts%rhcc(i), mgts%mla%la(i), nc = 1, ng = 1)
     call multifab_setval(mgts%rhcc(i),ZERO,all=.true.)
  end do

end subroutine mgt_alloc_rhcc_nodal

subroutine mgt_dealloc_rhcc_nodal()
  use nodal_cpp_mg_module
  implicit none
  integer :: i

  do i = 1, mgts%nlevel
     call multifab_destroy(mgts%rhcc(i))
  end do

  deallocate(mgts%rhcc)

end subroutine mgt_dealloc_rhcc_nodal

subroutine mgt_set_rhcc_nodal_1d(lev, n, rh, plo, phi, lo, hi)
  use nodal_cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(1), hi(1), plo(1), phi(1)
  real(kind=dp_t), intent(in) :: rh(plo(1):phi(1))
  real(kind=dp_t), pointer :: rp(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1
  rp => dataptr(mgts%rhcc(flev), fn)
  rp(lo(1):hi(1), 1,1,1) = rh(lo(1):hi(1))
end subroutine mgt_set_rhcc_nodal_1d

subroutine mgt_set_rhcc_nodal_2d(lev, n, rh, plo, phi, lo, hi)
  use nodal_cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(2), hi(2), plo(2), phi(2)
  real(kind=dp_t), intent(in) :: rh(plo(1):phi(1), plo(2):phi(2))
  real(kind=dp_t), pointer :: rp(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1
  rp => dataptr(mgts%rhcc(flev), fn)
  rp(lo(1):hi(1), lo(2):hi(2),1,1) = rh(lo(1):hi(1), lo(2):hi(2))
end subroutine mgt_set_rhcc_nodal_2d

subroutine mgt_set_rhcc_nodal_3d(lev, n, rh, plo, phi, lo, hi)
  use nodal_cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(3), hi(3), plo(3), phi(3)
  real(kind=dp_t), intent(in) :: rh(plo(1):phi(1), plo(2):phi(2), plo(3):phi(3))
  real(kind=dp_t), pointer :: rp(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1
  rp => dataptr(mgts%rhcc(flev), fn)
  rp(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3),1) = rh(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))
end subroutine mgt_set_rhcc_nodal_3d

subroutine mgt_add_divucc()

  use nodal_divu_module
  use nodal_cpp_mg_module

  integer    :: n
  real(dp_t) :: r, rhmax

  call divucc(mgts%nlevel,mgts%mgt,mgts%rhcc,mgts%rh,mgts%rr,mgts%nodal)

  if (mgts%verbose > 0) then
     rhmax = norm_inf(mgts%rh(mgts%nlevel),local=.true.)
     do n = mgts%nlevel-1, 1, -1
       r = norm_inf(mgts%rh(n),mgts%fine_mask(n),local=.true.)
       rhmax = max(r,rhmax) 
     end do 
     call parallel_reduce(r, rhmax, MPI_MAX, proc = parallel_IOProcessorNode())
     rhmax = r
     if (parallel_IOProcessor()) then
        print *,'F90: Source norm after adding cc    RHS is ',rhmax
     endif
  end if

end subroutine mgt_add_divucc
