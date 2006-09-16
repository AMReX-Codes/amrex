module nodal_cpp_mg_module
  use mg_module
  use ml_layout_module
  use ml_multifab_module
  use bndry_reg_module
  use nodal_mask_module
  use nodal_divu_module

  implicit none

  type mg_server
     logical         :: final = .false.
     integer         :: dim  = 0
     integer         :: nlevel
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
     type(multifab) , pointer ::       coeffs(:) => Null()
     type(multifab) , pointer ::   amr_coeffs(:) => Null()
     type(multifab) , pointer :: one_sided_ss(:) => Null()
     type(lmultifab), pointer ::    fine_mask(:) => Null()
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

subroutine mgt_nodal_alloc(dm, nlevel, nodal, stencil_type_in)
  use nodal_cpp_mg_module
  implicit none
  integer, intent(in) :: dm, nlevel, stencil_type_in
  integer :: nodal
  integer i

  if ( mgts%dim == 0 ) then
     mgts%dim = dm
     mgts%nlevel = nlevel
     allocate(mgts%nodal(dm))
     mgts%nodal = (nodal /= 0)
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
  integer   :: nc
  logical   :: pmask(dm)
  type(box) :: pd
  integer   :: flev
  logical, allocatable :: nodal(:)

  flev = lev + 1
  call mgt_verify_lev("MGT_SET_NODAL_LEVEL", flev)

  pmask = (pm /= 0)

  allocate(nodal(dm))

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
  implicit none
  real(dp_t), intent(in) :: dx(mgts%nlevel,mgts%dim)
  integer   , intent(in) :: bc(2,mgts%dim)
  integer :: dm, i, nlev, n
  integer :: ns
  integer :: nc
  logical, allocatable :: nodal(:)

  integer :: max_nlevel_in
  integer :: bottom_solver_in

  type(boxarray) :: bac
  type(layout)   :: la

  integer :: bottom_max_iter_in

  call mgt_verify("MGT_FINALIZE")

  dm = mgts%dim

  allocate(mgts%bc(dm,2))
  mgts%bc = transpose(bc)

  allocate(nodal(1:dm))
  nodal = mgts%nodal

  nc = 1
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

  if (mgts%stencil_type .eq. ST_DENSE) then
     if (dm .eq. 3) then
       i = mgts%mgt(nlev)%nlevels
       if ( (mgts%mgt(nlev)%dh(1,i) .eq. mgts%mgt(nlev)%dh(2,i)) .and. &
            (mgts%mgt(nlev)%dh(1,i) .eq. mgts%mgt(nlev)%dh(3,i)) ) then
         ns = 21
       else
         ns = 27
       end if
     else if (dm .eq. 2) then
       ns = 9
     end if
     if ( parallel_ioprocessor() ) print *,'SETTING UP DENSE STENCIL WITH NS = ',ns
  else
    ns = 2*dm+1
    do n = nlev, 2, -1
      call multifab_build(mgts%one_sided_ss(n), mgts%mla%la(n), ns, 0, nodal)
      call setval(mgts%one_sided_ss(n), 0.0_dp_t, all=.true.)
    end do
  end if

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

end subroutine mgt_nodal_finalize

subroutine mgt_init_nodal_coeffs_lev(lev)
  use nodal_cpp_mg_module
  implicit none
  integer, intent(in) :: lev
  integer :: nlev, dm, i
  integer :: flev
  flev = lev + 1
  call mgt_verify_lev("MGT_INIT_STENCIL_LEV", flev)

  dm = mgts%dim
  nlev = mgts%mgt(flev)%nlevels
  allocate(mgts%coeffs(nlev))

  call  build(mgts%coeffs(nlev), mgts%mgt(flev)%ss(nlev)%la, 1, 1)
  call setval(mgts%coeffs(nlev), 0.0_dp_t, all=.true.)

  do i = nlev-1, 1, -1
     call build(mgts%coeffs(i), mgts%mgt(flev)%ss(i)%la, 1, 1)
     call setval(mgts%coeffs(i), 0.0_dp_t, all=.true.)
  end do

  ! These only exist at amr levels, not the lower multigrid levels
  call  build(mgts%amr_coeffs(flev), mgts%mgt(flev)%ss(nlev)%la, 1, 1)
  call setval(mgts%amr_coeffs(flev), 0.0_dp_t, all=.true.)

end subroutine mgt_init_nodal_coeffs_lev

subroutine mgt_finalize_nodal_stencil_lev(lev)
  use nodal_cpp_mg_module
  use coeffs_module
  implicit none
  integer, intent(in) :: lev
  integer :: nlev, i, dm
  integer :: flev
  type(box) :: pd
  type(boxarray) :: pdv
  dm = mgts%dim
  flev = lev + 1
  call mgt_verify_lev("MGT_FINALIZE_NODAL_STENCIL_LEV", flev)
  
  pd = mgts%pd(flev)
  nlev = mgts%mgt(flev)%nlevels

  call multifab_fill_boundary(mgts%coeffs(nlev))

  do i = nlev, 1, -1

     if (i < nlev) then
       call coarsen_coeffs(mgts%coeffs(i+1), mgts%coeffs(i))
       call multifab_fill_boundary(mgts%coeffs(i))
     end if

     call stencil_fill_nodal(mgts%mgt(flev)%ss(i), mgts%coeffs(i), &
                             mgts%mgt(flev    )%dh(:,i), &
                             mgts%mgt(flev)%dh(:,mgts%mgt(flev)%nlevels), &
                             mgts%mgt(flev)%mm(i), mgts%mgt(flev)%face_type, mgts%stencil_type)

  end do

  if (mgts%stencil_type .eq. ST_CROSS .and. flev .gt. 1) then
     call stencil_fill_one_sided(mgts%one_sided_ss(flev), mgts%coeffs(nlev), &
                                 mgts%mgt(flev)%dh(:,nlev), &
                                 mgts%mgt(flev)%mm(nlev), mgts%mgt(flev)%face_type)
  end if

  deallocate(mgts%coeffs)

end subroutine mgt_finalize_nodal_stencil_lev

subroutine mgt_divu()
  use nodal_cpp_mg_module

  integer    :: n
  real(dp_t) :: r, rhmax

  call divu(mgts%nlevel,mgts%mgt,mgts%vel,mgts%rh,mgts%rr,mgts%verbose,mgts%nodal)

  if (parallel_IOProcessor() .and. mgts%verbose > 0) then
     rhmax = norm_inf(mgts%rh(mgts%nlevel))
     do n = mgts%nlevel-1, 1, -1
       r = norm_inf(mgts%rh(n), mgts%fine_mask(n))
       rhmax = max(r, rhmax) 
     end do 
     print *,'F90: Source norm is ',rhmax
  end if

end subroutine mgt_divu

subroutine mgt_newu()
  use nodal_cpp_mg_module
  use nodal_newu_module

  integer :: i,n

  do n = 1,mgts%nlevel
     call mkunew(mgts%vel(n),mgts%uu(n),mgts%amr_coeffs(n), &
                 mgts%mgt(n)%dh(:,mgts%mgt(n)%nlevels), &
                 mgts%vel(n)%ng)
  end do
       
end subroutine mgt_newu

subroutine mgt_finalize_nodal_stencil()
   use nodal_cpp_mg_module
   implicit none
   integer :: i
   call mgt_verify("MGT_FINALIZE_STENCIL")

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

subroutine mgt_set_vel_2d(lev, n, vel_in, plo, phi, lo, hi)
  use nodal_cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(2), hi(2), plo(2), phi(2)
  real(kind=dp_t), intent(in) :: vel_in(plo(1):phi(1), plo(2):phi(2), 2)
  real(kind=dp_t), pointer :: vp(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1
  
  call mgt_verify_n("MGT_SET_VEL_2D", flev, fn, lo, hi)

  vp => dataptr(mgts%vel(flev), fn)
  vp(:,:,:,:) = ZERO
  vp(lo(1):hi(1), lo(2):hi(2),1,1) = vel_in(lo(1):hi(1), lo(2):hi(2),1)
  vp(lo(1):hi(1), lo(2):hi(2),1,2) = vel_in(lo(1):hi(1), lo(2):hi(2),2)

end subroutine mgt_set_vel_2d

subroutine mgt_get_vel_2d(lev, n, vel_out, plo, phi, lo, hi)
  use nodal_cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(2), hi(2), plo(2), phi(2)
  real(kind=dp_t), intent(inout) :: vel_out(plo(1):phi(1), plo(2):phi(2), 2)
  real(kind=dp_t), pointer :: vp(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1

  vp => dataptr(mgts%vel(flev), fn)
  vel_out(lo(1):hi(1), lo(2):hi(2),1) = vp(lo(1):hi(1), lo(2):hi(2),1,1)
  vel_out(lo(1):hi(1), lo(2):hi(2),2) = vp(lo(1):hi(1), lo(2):hi(2),1,2)

end subroutine mgt_get_vel_2d

subroutine mgt_set_vel_3d(lev, n, vel_in, plo, phi, lo, hi)
  use nodal_cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(3), hi(3), plo(3), phi(3)
  real(kind=dp_t), intent(in) :: vel_in(plo(1):phi(1), plo(2):phi(2), plo(3):phi(3), 3)
  real(kind=dp_t), pointer :: vp(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1
  
  call mgt_verify_n("MGT_SET_VEL_3D", flev, fn, lo, hi)

  vp => dataptr(mgts%vel(flev), fn)
  vp(:,:,:,:) = ZERO
  vp(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3),1) = vel_in(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3),1)
  vp(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3),2) = vel_in(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3),2)
  vp(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3),3) = vel_in(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3),3)

end subroutine mgt_set_vel_3d

subroutine mgt_get_vel_3d(lev, n, vel_out, plo, phi, lo, hi)
  use nodal_cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(3), hi(3), plo(3), phi(3)
  real(kind=dp_t), intent(inout) :: vel_out(plo(1):phi(1), plo(2):phi(2), plo(3):phi(3), 3)
  real(kind=dp_t), pointer :: vp(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1

  vp => dataptr(mgts%vel(flev), fn)
  vel_out(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3),1) = vp(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3),1)
  vel_out(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3),2) = vp(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3),2)
  vel_out(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3),3) = vp(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3),3)

end subroutine mgt_get_vel_3d

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
  nlev = size(mgts%coeffs)
  call mgt_verify_n("MGT_SET_CFS_2D", flev, fn, lo, hi)

  cp => dataptr(mgts%coeffs(nlev), fn)
  cp(lo(1):hi(1), lo(2):hi(2), 1, 1) = cf(lo(1):hi(1), lo(2):hi(2))

  acp => dataptr(mgts%amr_coeffs(flev), fn)
  acp(lo(1):hi(1), lo(2):hi(2), 1, 1) = cf(lo(1):hi(1), lo(2):hi(2))

end subroutine mgt_set_cfs_2d

subroutine mgt_set_pr_2d(lev, n, uu, plo, phi, lo, hi)
  use nodal_cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(2), hi(2), plo(2), phi(2)
  real(kind=dp_t), intent(in) :: uu(plo(1):phi(1), plo(2):phi(2))
  real(kind=dp_t), pointer :: up(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1

  up => dataptr(mgts%uu(flev), fn)
  up(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1,1,1) = uu(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1)

end subroutine mgt_set_pr_2d

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
  nlev = size(mgts%coeffs)
  call mgt_verify_n("MGT_SET_CFS_3D", flev, fn, lo, hi)

  cp => dataptr(mgts%coeffs(nlev), fn)
  cp(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 1) = cf(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))

  acp => dataptr(mgts%amr_coeffs(flev), fn)
  acp(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 1) = cf(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))

end subroutine mgt_set_cfs_3d

subroutine mgt_set_pr_3d(lev, n, uu, plo, phi, lo, hi)
  use nodal_cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(3), hi(3), plo(3), phi(3)
  real(kind=dp_t), intent(in) :: uu(plo(1):phi(1), plo(2):phi(2), plo(3):phi(3))
  real(kind=dp_t), pointer :: up(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1

  up => dataptr(mgts%uu(flev), fn)
  up(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, lo(3)-1:hi(3)+1, 1) = &
     uu(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, lo(3)-1:hi(3)+1)

end subroutine mgt_set_pr_3d

subroutine mgt_get_pr_2d(lev, n, uu, plo, phi, lo, hi)
  use nodal_cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(2), hi(2), plo(2), phi(2)
  real(kind=dp_t), intent(inout) :: uu(plo(1):phi(1), plo(2):phi(2))
  real(kind=dp_t), pointer :: up(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1

  up => dataptr(mgts%uu(flev), fn)
  uu(lo(1):hi(1), lo(2):hi(2)) = up(lo(1):hi(1), lo(2):hi(2),1,1)

end subroutine mgt_get_pr_2d

subroutine mgt_get_pr_3d(lev, n, uu, plo, phi, lo, hi)
  use nodal_cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(3), hi(3), plo(3), phi(3)
  real(kind=dp_t), intent(inout) :: uu(plo(1):phi(1), plo(2):phi(2), plo(3):phi(3))
  real(kind=dp_t), pointer :: up(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1

  up => dataptr(mgts%uu(flev), fn)
  uu(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)) = up(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 1)

end subroutine mgt_get_pr_3d

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
     call destroy(mgts%amr_coeffs(i))
  end do
  do i = 1,mgts%nlevel-1
     call destroy(mgts%fine_mask(i))
  end do
  do i = 2,mgts%nlevel
     call destroy(mgts%one_sided_ss(i))
  end do
  call destroy(mgts%mla)
  mgts%dim = 0
  mgts%final = .false.

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

  do_diagnostics = 0
  call ml_nd(mgts%mla, mgts%mgt, &
       mgts%rh, mgts%uu, &
       mgts%fine_mask, &
       mgts%one_sided_ss(2:), mgts%rr, &
       do_diagnostics, tol)

end subroutine mgt_nodal_solve

subroutine mgt_set_nodal_defaults(nu_1,nu_2,nu_b,nu_f,gamma,omega,max_iter,bottom_max_iter, &
                                  bottom_solver,bottom_solver_eps, &
                                  verbose,cg_verbose,max_nlevel,min_width,cycle,smoother,stencil_type)
  use nodal_cpp_mg_module
  implicit none
  integer   , intent(in) :: nu_1,nu_2,nu_b,nu_f,gamma,max_iter,bottom_max_iter,bottom_solver
  integer   , intent(in) :: verbose, cg_verbose, max_nlevel, min_width, cycle, smoother, stencil_type
  real(dp_t), intent(in) :: omega, bottom_solver_eps

  call mgt_not_final("MGT_SET_NODAL_DEFAULTS")

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

  mgts%stencil_type    = stencil_type

end subroutine mgt_set_nodal_defaults

subroutine mgt_get_nodal_defaults(nu_1,nu_2,nu_b,nu_f,gamma,omega,max_iter,bottom_max_iter, &
                                  bottom_solver, &
                                  verbose,cg_verbose,max_nlevel,min_width,cycle,smoother)
  use nodal_cpp_mg_module
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

end subroutine mgt_get_nodal_defaults
