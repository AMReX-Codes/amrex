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
     type(multifab) , pointer ::  cell_coeffs(:) => Null()
     type(multifab) , pointer ::   amr_coeffs(:) => Null()
     type(multifab) , pointer :: one_sided_ss(:) => Null()
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

subroutine mgt_nodal_alloc(dm, nlevel, stencil_type_in)
  use nodal_cpp_mg_module
  implicit none
  integer, intent(in) :: dm, nlevel, stencil_type_in

  if ( mgts%dim == 0 ) then
     mgts%dim = dm
     mgts%nlevel = nlevel
     allocate(mgts%nodal(dm))
     mgts%nodal = .true.
  end if

  mgts%stencil_type = stencil_type_in

  allocate(mgts%rr(nlevel-1,dm))
  allocate(mgts%rh(nlevel))
  allocate(mgts%vel(nlevel))
  allocate(mgts%pd(nlevel))
  allocate(mgts%uu(nlevel))
  allocate(mgts%mgt(nlevel))
  allocate(mgts%amr_coeffs(nlevel))
  allocate(mgts%one_sided_ss(nlevel))
  allocate(mgts%fine_mask(nlevel))

  call build(mgts%mla, nlevel, dm)

end subroutine mgt_nodal_alloc

subroutine mgt_set_nodal_level(lev, nb, dm, lo, hi, pd_lo, pd_hi, pm, pmap)
  use nodal_cpp_mg_module
  implicit none
  integer, intent(in) :: lev, nb, dm
  integer, intent(in) :: lo(nb,dm), hi(nb,dm), pd_lo(dm), pd_hi(dm), pm(dm), pmap(nb+1)

  type(box) :: bxs(nb)
  integer   :: i
  logical   :: pmask(dm)
  integer   :: flev

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
  call build(mgts%mla%mba%bas(flev), bxs)
  call build(mgts%mla%la(flev),  &
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
  integer :: ns
  logical, allocatable :: nodal(:)

  integer :: max_nlevel_in
  integer :: bottom_solver_in

  type(boxarray) :: bac

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
     call build(mgts%uu(i) , mgts%mla%la(i), nc = 1, ng = 1, nodal = nodal)
     call build(mgts%rh(i) , mgts%mla%la(i), nc = 1, ng = 1, nodal = nodal)
     call build(mgts%vel(i), mgts%mla%la(i), nc =dm, ng = 1)

     call setval(mgts%uu(i),ZERO,all=.true.)
     call setval(mgts%rh(i),ZERO,all=.true.)
     call setval(mgts%vel(i),ZERO,all=.true.)
  end do

  do i = nlev-1, 1, -1
     call build(mgts%mla%mask(i), mgts%mla%la(i), nc = 1, ng = 0)
     call setval(mgts%mla%mask(i), val = .TRUE.) 
     call copy(bac, mgts%mla%mba%bas(i+1))
     call boxarray_coarsen(bac, mgts%rr(i,:))
     call setval(mgts%mla%mask(i), .false., bac)
     call destroy(bac) 
  end do

  if (mgts%stencil_type .eq. ND_DENSE_STENCIL) then

     if ( parallel_ioprocessor() .and. mgts%verbose > 0 ) &
         print *,'Using dense stencil in nodal solver ...'

     if (dm .eq. 3) then
       if ( dx(nlev,1) .eq. dx(nlev,2) .and. dx(nlev,1) .eq. dx(nlev,3) ) then
         ns = 21
       else
         ns = 27
       end if
     else if (dm .eq. 2) then
       ns = 9
     end if

  else if (mgts%stencil_type .eq. ND_CROSS_STENCIL) then

     if ( parallel_ioprocessor() .and. mgts%verbose > 0 ) &
         print *,'Using cross stencil in nodal solver ...'

     ns = 2*dm+1
     do n = nlev, 2, -1
       call multifab_build(mgts%one_sided_ss(n), mgts%mla%la(n), ns, 0, nodal, stencil=.true.)
     end do

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
          ns                = ns, &
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

subroutine mgt_init_nodal_coeffs_lev(lev)
  use nodal_cpp_mg_module
  implicit none
  integer, intent(in) :: lev
  integer :: nlev
  integer :: flev
  flev = lev + 1
  call mgt_verify_lev("MGT_INIT_STENCIL_LEV", flev)

  nlev = mgts%mgt(flev)%nlevels
  allocate(mgts%cell_coeffs(nlev))

  call  build(mgts%cell_coeffs(nlev), mgts%mgt(flev)%ss(nlev)%la, 1, 1)
  call setval(mgts%cell_coeffs(nlev), 0.0_dp_t, all=.true.)

  ! These only exist at amr levels, not the lower multigrid levels
  call  build(mgts%amr_coeffs(flev), mgts%mgt(flev)%ss(nlev)%la, 1, 1)
  call setval(mgts%amr_coeffs(flev), 0.0_dp_t, all=.true.)

end subroutine mgt_init_nodal_coeffs_lev

subroutine mgt_finalize_nodal_stencil_lev(lev)
  use nodal_cpp_mg_module
  use nodal_stencil_fill_module
  implicit none
  integer, intent(in) :: lev
  integer :: nlev
  integer :: flev

  flev = lev + 1

  call mgt_verify_lev("MGT_FINALIZE_NODAL_STENCIL_LEV", flev)
  
  nlev = mgts%mgt(flev)%nlevels

  call multifab_fill_boundary(mgts%cell_coeffs(nlev))

  call stencil_fill_nodal_all_mglevels(mgts%mgt(flev), mgts%cell_coeffs, mgts%stencil_type)

  if (mgts%stencil_type .eq. ND_CROSS_STENCIL .and. flev .gt. 1) then
     call stencil_fill_one_sided(mgts%one_sided_ss(flev), mgts%cell_coeffs(nlev), &
                                 mgts%mgt(flev)%dh(:,nlev), &
                                 mgts%mgt(flev)%mm(nlev), mgts%mgt(flev)%face_type)
  end if

  call destroy(mgts%cell_coeffs(nlev))
  deallocate(mgts%cell_coeffs)

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
     call build(mgts%fine_mask(i), mgts%mla%la(i), nc = 1, ng = 0, nodal = mgts%nodal)
     call setval(mgts%fine_mask(i), val = .TRUE., all = .TRUE.)
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

  call mgt_verify_n("MGT_SET_VEL_1D", flev, fn, lo, hi)

  vp => dataptr(mgts%vel(flev), local_index(mgts%vel(flev),fn))
  vp(lo(1)-1:hi(1)+1,1,1,1) = vel_in(lo(1)-1:hi(1)+1,iv+1)

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

  vp => dataptr(mgts%vel(flev), local_index(mgts%vel(flev),fn))
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
  
  call mgt_verify_n("MGT_SET_VEL_2D", flev, fn, lo, hi)

  vp => dataptr(mgts%vel(flev), local_index(mgts%vel(flev),fn))
  vp(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, 1, 1:2) =  &
       vel_in(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, iv+1:iv+2)

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

  vp => dataptr(mgts%vel(flev), local_index(mgts%vel(flev),fn))
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
  
  call mgt_verify_n("MGT_SET_VEL_3D", flev, fn, lo, hi)

  vp => dataptr(mgts%vel(flev), local_index(mgts%vel(flev),fn))
  vp(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, lo(3)-1:hi(3)+1, 1:3) =  &
       vel_in(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, lo(3)-1:hi(3)+1, iv+1:iv+3)

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

  vp => dataptr(mgts%vel(flev), local_index(mgts%vel(flev),fn))
  vel_out(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), iv+1:iv+3) = &
       vp(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 1:3)

end subroutine mgt_get_vel_3d

subroutine mgt_set_cfs_1d(lev, n, cf, plo, phi, lo, hi)
  use nodal_cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(1), hi(1), plo(1), phi(1)
  real(kind=dp_t), intent(in) :: cf(plo(1):phi(1))
  real(kind=dp_t), pointer :: cp(:,:,:,:)
  real(kind=dp_t), pointer :: acp(:,:,:,:)
  integer :: flev, fn, nlev
  fn = n + 1
  flev = lev+1
  nlev = size(mgts%cell_coeffs)
  call mgt_verify_n("MGT_SET_CFS_1D", flev, fn, lo, hi)

  cp => dataptr(mgts%cell_coeffs(nlev), local_index(mgts%cell_coeffs(nlev),fn))
  cp(lo(1):hi(1), 1, 1, 1) = cf(lo(1):hi(1))

  acp => dataptr(mgts%amr_coeffs(flev), local_index(mgts%amr_coeffs(flev),fn))
  acp(lo(1):hi(1), 1, 1, 1) = cf(lo(1):hi(1))

end subroutine mgt_set_cfs_1d

subroutine mgt_set_cfs_2d(lev, n, cf, plo, phi, lo, hi)
  use nodal_cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(2), hi(2), plo(2), phi(2)
  real(kind=dp_t), intent(in) :: cf(plo(1):phi(1), plo(2):phi(2))
  real(kind=dp_t), pointer :: cp(:,:,:,:)
  real(kind=dp_t), pointer :: acp(:,:,:,:)
  integer :: flev, fn, nlev
  fn = n + 1
  flev = lev+1
  nlev = size(mgts%cell_coeffs)
  call mgt_verify_n("MGT_SET_CFS_2D", flev, fn, lo, hi)

  cp => dataptr(mgts%cell_coeffs(nlev), local_index(mgts%cell_coeffs(nlev),fn))
  cp(lo(1):hi(1), lo(2):hi(2), 1, 1) = cf(lo(1):hi(1), lo(2):hi(2))

  acp => dataptr(mgts%amr_coeffs(flev), local_index(mgts%amr_coeffs(flev),fn))
  acp(lo(1):hi(1), lo(2):hi(2), 1, 1) = cf(lo(1):hi(1), lo(2):hi(2))

end subroutine mgt_set_cfs_2d

subroutine mgt_set_cfs_3d(lev, n, cf, plo, phi, lo, hi)
  use nodal_cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(3), hi(3), plo(3), phi(3)
  real(kind=dp_t), intent(in) :: cf(plo(1):phi(1), plo(2):phi(2), plo(3):phi(3))
  real(kind=dp_t), pointer :: cp(:,:,:,:)
  real(kind=dp_t), pointer :: acp(:,:,:,:)
  integer :: flev, fn, nlev
  fn = n + 1
  flev = lev+1
  nlev = size(mgts%cell_coeffs)
  call mgt_verify_n("MGT_SET_CFS_3D", flev, fn, lo, hi)

  cp => dataptr(mgts%cell_coeffs(nlev), local_index(mgts%cell_coeffs(nlev),fn))
  cp(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 1) = cf(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))

  acp => dataptr(mgts%amr_coeffs(flev), local_index(mgts%amr_coeffs(flev),fn))
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

  up => dataptr(mgts%uu(flev), local_index(mgts%uu(flev),fn))
  up(lo(1)-1:hi(1)+1,1,1,1) = uu(lo(1)-1:hi(1)+1, ip+1)

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

  up => dataptr(mgts%uu(flev), local_index(mgts%uu(flev),fn))
  up(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1,1,1) = uu(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, ip+1)

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

  up => dataptr(mgts%uu(flev), local_index(mgts%uu(flev),fn))
  up(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, lo(3)-1:hi(3)+1, 1) = &
     uu(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, lo(3)-1:hi(3)+1, ip+1)

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

  up => dataptr(mgts%uu(flev), local_index(mgts%uu(flev),fn))
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

  up => dataptr(mgts%uu(flev), local_index(mgts%uu(flev),fn))
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

  up => dataptr(mgts%uu(flev), local_index(mgts%uu(flev),fn))
  uu(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), ip+1) =  &
       up(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 1)

end subroutine mgt_get_pr_3d

subroutine mgt_add_rh_nodal_1d(lev, n, rh_in, plo, phi, lo, hi)
  use nodal_cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(1), hi(1), plo(1), phi(1)
  real(kind=dp_t), intent(in) :: rh_in(plo(1):phi(1))
  real(kind=dp_t), pointer :: rp(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1

  rp => dataptr(mgts%rh(flev), local_index(mgts%rh(flev),fn))
  rp(lo(1):hi(1),1,1,1) = rp(lo(1):hi(1),1,1,1) + rh_in(lo(1):hi(1))

end subroutine mgt_add_rh_nodal_1d

subroutine mgt_add_rh_nodal_2d(lev, n, rh_in, plo, phi, lo, hi)
  use nodal_cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(2), hi(2), plo(2), phi(2)
  real(kind=dp_t), intent(in) :: rh_in(plo(1):phi(1), plo(2):phi(2))
  real(kind=dp_t), pointer :: rp(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1

  rp => dataptr(mgts%rh(flev), local_index(mgts%rh(flev),fn))
  rp(lo(1):hi(1),lo(2):hi(2),1,1) = &
       rp(lo(1):hi(1),lo(2):hi(2),1,1) +  &
       rh_in(lo(1):hi(1),lo(2):hi(2))

end subroutine mgt_add_rh_nodal_2d

subroutine mgt_add_rh_nodal_3d(lev, n, rh_in, plo, phi, lo, hi)
  use nodal_cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(3), hi(3), plo(3), phi(3)
  real(kind=dp_t), intent(in) :: rh_in(plo(1):phi(1),plo(2):phi(2),plo(3):phi(3))
  real(kind=dp_t), pointer :: rp(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1

  rp => dataptr(mgts%rh(flev), local_index(mgts%rh(flev),fn))
  rp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1) = &
       rp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1) + &
       rh_in(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))

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
     call destroy(mgts%rh(i))
     call destroy(mgts%uu(i))
     call destroy(mgts%vel(i))
     call destroy(mgts%amr_coeffs(i))
  end do
  do i = 1,mgts%nlevel-1
     call destroy(mgts%fine_mask(i))
  end do
  if (mgts%stencil_type .eq. ND_CROSS_STENCIL) then
     do i = 2,mgts%nlevel
        call destroy(mgts%one_sided_ss(i))
     end do
  endif
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
  deallocate(mgts%one_sided_ss)
  deallocate(mgts%fine_mask)

  call parallel_finalize(.false.) ! do not finalize MPI but free communicator

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

  call ml_nd(mgts%mla, mgts%mgt, &
       mgts%rh, mgts%uu, &
       mgts%fine_mask, &
       mgts%one_sided_ss(2:), mgts%rr, &
       do_diagnostics, tol, abs_tol)

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
     call build(mgts%rhcc(i), mgts%mla%la(i), nc = 1, ng = 1)
     call setval(mgts%rhcc(i),ZERO,all=.true.)
  end do

end subroutine mgt_alloc_rhcc_nodal

subroutine mgt_dealloc_rhcc_nodal()
  use nodal_cpp_mg_module
  implicit none
  integer :: i

  do i = 1, mgts%nlevel
     call destroy(mgts%rhcc(i))
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
  call mgt_verify_n("MGT_SET_RHCC_NODAL_1D", flev, fn, lo, hi)

  rp => dataptr(mgts%rhcc(flev), local_index(mgts%rhcc(flev),fn))
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
  
  call mgt_verify_n("MGT_SET_RHCC_NODAL_2D", flev, fn, lo, hi)

  rp => dataptr(mgts%rhcc(flev), local_index(mgts%rhcc(flev),fn))
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
  
  call mgt_verify_n("MGT_SET_RHCC_NODAL_3D", flev, fn, lo, hi)

  rp => dataptr(mgts%rhcc(flev), local_index(mgts%rhcc(flev),fn))
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
        print *,'F90: Source norm after adding rhs is ',rhmax
     endif
  end if

end subroutine mgt_add_divucc
