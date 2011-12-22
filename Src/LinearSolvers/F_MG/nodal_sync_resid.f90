module nodal_sync_resid_module
  use bl_constants_module
  implicit none

contains

  subroutine get_phi_crse_2d(pc, msk, uu)
    real(kind=dp_t), intent(inout) :: pc(-1:,-1:)
    real(kind=dp_t), intent(in   ) :: msk(-1:,-1:)
    real(kind=dp_t), intent(in   ) :: uu (-1:,-1:)
    
    integer :: i, j, nx, ny

    nx = size(msk,dim=1) - 2
    ny = size(msk,dim=2) - 2

    do j = 0, ny
       do i = 0, nx
          if (all(msk(i-1:i,j-1:j) .eq. ZERO)) then
             ! covered by fine
             pc(i,j) = ZERO
          else
             pc(i,j) = uu(i,j)             
          end if
       end do
    end do

  end subroutine get_phi_crse_2d

  subroutine get_phi_fine_2d(pf, msk, uu)
    real(kind=dp_t), intent(inout) :: pf(-1:,-1:)
    real(kind=dp_t), intent(in   ) :: msk(-1:,-1:)
    real(kind=dp_t), intent(in   ) :: uu (-1:,-1:)
    
    integer :: i, j, nx, ny

    nx = size(msk,dim=1) - 2
    ny = size(msk,dim=2) - 2

    do j = 0, ny
       do i = 0, nx
          if (all(msk(i-1:i,j-1:j) .eq. ONE)) then
             ! interior node
             pf(i,j) = uu(i,j)
          else
             pf(i,j) = ZERO
          end if
       end do
    end do

  end subroutine get_phi_fine_2d

  subroutine add_Dv_crse_2d(res, msk, vel, dx)
    real(kind=dp_t), intent(inout) :: res(-1:,-1:)
    real(kind=dp_t), intent(in   ) :: msk(-1:,-1:)
    real(kind=dp_t), intent(in   ) :: vel(-1:,-1:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    
    integer :: i, j, nx, ny
    real(kind=dp_t) :: divv

    nx = size(msk,dim=1) - 2
    ny = size(msk,dim=2) - 2

    do j = 0, ny
       do i = 0, nx
          if ( any(msk(i-1:i,j-1:j) .eq. ZERO) .and. &
               any(msk(i-1:i,j-1:j) .eq. ONE)) then
             divv = (vel(i  ,j,1)*msk(i , j) + vel(i  ,j-1,1)*msk(i  ,j-1) &
                  -  vel(i-1,j,1)*msk(i-1,j) - vel(i-1,j-1,1)*msk(i-1,j-1)) / dx(1) &
                  + (vel(i,j  ,2)*msk(i,j  ) + vel(i-1,j  ,2)*msk(i-1,j  ) &
                  -  vel(i,j-1,2)*msk(i,j-1) - vel(i-1,j-1,2)*msk(i-1,j-1)) / dx(2)
             divv = HALF * divv
             res(i,j) = divv  ! res(i,j) + divv
          else
             res(i,j) = ZERO
          end if
       end do
    end do

  end subroutine add_Dv_crse_2d

  subroutine add_Dv_fine_2d(res, msk, vel, dx)
    real(kind=dp_t), intent(inout) :: res(-1:,-1:)
    real(kind=dp_t), intent(in   ) :: msk(-1:,-1:)
    real(kind=dp_t), intent(in   ) :: vel(-1:,-1:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    
    integer :: i, j, nx, ny
    real(kind=dp_t) :: divv

    nx = size(msk,dim=1) - 2
    ny = size(msk,dim=2) - 2

    do j = 0, ny
       do i = 0, nx
          if ( any(msk(i-1:i,j-1:j) .eq. ZERO) .and. &
               any(msk(i-1:i,j-1:j) .eq. ONE)) then
             divv = (vel(i  ,j,1)*msk(i , j) + vel(i  ,j-1,1)*msk(i  ,j-1) &
                  -  vel(i-1,j,1)*msk(i-1,j) - vel(i-1,j-1,1)*msk(i-1,j-1)) / dx(1) &
                  + (vel(i,j  ,2)*msk(i,j  ) + vel(i-1,j  ,2)*msk(i-1,j  ) &
                  -  vel(i,j-1,2)*msk(i,j-1) - vel(i-1,j-1,2)*msk(i-1,j-1)) / dx(2)
             divv = HALF * divv
             res(i,j) = divv ! res(i,j) + divv
          else
             res(i,j) = ZERO
          end if
       end do
    end do

  end subroutine add_Dv_fine_2d

end module nodal_sync_resid_module


subroutine mgt_alloc_nodal_sync()
  use nodal_cpp_mg_module
  implicit none
  logical,dimension(3) :: nodal

  allocate(mgts%sync_res(1))
  allocate(mgts%sync_msk(1))
  allocate(mgts%vold(1))

  nodal = .true.

  call build(mgts%sync_res(1) , mgts%mla%la(1), nc = 1, ng = 1, nodal = nodal)
  call build(mgts%sync_msk(1) , mgts%mla%la(1), nc = 1, ng = 1)
  call build(mgts%vold(1)     , mgts%mla%la(1), nc = mgts%dim, ng = 1)

  call setval(mgts%sync_res(1),ZERO,all=.true.)

end subroutine mgt_alloc_nodal_sync

subroutine mgt_dealloc_nodal_sync()
  use nodal_cpp_mg_module
  implicit none
  
  call destroy(mgts%sync_res(1))
  call destroy(mgts%sync_msk(1))
  call destroy(mgts%vold(1))

  deallocate(mgts%sync_res)
  deallocate(mgts%sync_msk)
  deallocate(mgts%vold)

end subroutine mgt_dealloc_nodal_sync

subroutine mgt_set_sync_msk_1d(lev, n, msk_in, plo, phi, lo, hi)
  use nodal_cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(1), hi(1), plo(1), phi(1)
  real(kind=dp_t), intent(in) :: msk_in(plo(1):phi(1))
  real(kind=dp_t), pointer :: mskp(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1

  mskp => dataptr(mgts%sync_msk(flev), fn)
  mskp(plo(1):phi(1),1,1,1) = msk_in(plo(1):phi(1))
end subroutine mgt_set_sync_msk_1d

subroutine mgt_set_sync_msk_2d(lev, n, msk_in, plo, phi, lo, hi)
  use nodal_cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(2), hi(2), plo(2), phi(2)
  real(kind=dp_t), intent(in) :: msk_in(plo(1):phi(1),plo(2):phi(2))
  real(kind=dp_t), pointer :: mskp(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1

  mskp => dataptr(mgts%sync_msk(flev), fn)
  mskp(plo(1):phi(1),plo(2):phi(2),1,1) = msk_in(plo(1):phi(1),plo(2):phi(2))
end subroutine mgt_set_sync_msk_2d

subroutine mgt_set_sync_msk_3d(lev, n, msk_in, plo, phi, lo, hi)
  use nodal_cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(3), hi(3), plo(3), phi(3)
  real(kind=dp_t), intent(in) :: msk_in(plo(1):phi(1),plo(2):phi(2),plo(3):phi(3))
  real(kind=dp_t), pointer :: mskp(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1

  mskp => dataptr(mgts%sync_msk(flev), fn)
  mskp(plo(1):phi(1),plo(2):phi(2),plo(3):phi(3),1) = &
       msk_in(plo(1):phi(1),plo(2):phi(2),plo(3):phi(3))
end subroutine mgt_set_sync_msk_3d

subroutine mgt_set_vold_1d(lev, n, v_in, plo, phi, lo, hi)
  use nodal_cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(1), hi(1), plo(1), phi(1)
  real(kind=dp_t), intent(in) :: v_in(plo(1):phi(1))
  real(kind=dp_t), pointer :: vp(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1

  vp => dataptr(mgts%vold(flev), fn)
  vp(plo(1):phi(1),1,1,1) = v_in(plo(1):phi(1))
end subroutine mgt_set_vold_1d

subroutine mgt_set_vold_2d(lev, n, v_in, plo, phi, lo, hi)
  use nodal_cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(2), hi(2), plo(2), phi(2)
  real(kind=dp_t), intent(in) :: v_in(plo(1):phi(1),plo(2):phi(2),1:2)
  real(kind=dp_t), pointer :: vp(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1

  vp => dataptr(mgts%vold(flev), fn)
  vp(plo(1):phi(1),plo(2):phi(2),1,1:2) = v_in(plo(1):phi(1),plo(2):phi(2),1:2)
end subroutine mgt_set_vold_2d

subroutine mgt_set_vold_3d(lev, n, v_in, plo, phi, lo, hi)
  use nodal_cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(3), hi(3), plo(3), phi(3)
  real(kind=dp_t), intent(in) :: v_in(plo(1):phi(1),plo(2):phi(2),plo(3):phi(3),1:3)
  real(kind=dp_t), pointer :: vp(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1

  vp => dataptr(mgts%vold(flev), fn)
  vp(plo(1):phi(1),plo(2):phi(2),plo(3):phi(3),1:3) = &
       v_in(plo(1):phi(1),plo(2):phi(2),plo(3):phi(3),1:3)
end subroutine mgt_set_vold_3d

subroutine mgt_get_sync_res_1d(lev, n, res, plo, phi, lo, hi)
  use nodal_cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(1), hi(1), plo(1), phi(1)
  real(kind=dp_t), intent(inout) :: res(plo(1):phi(1))
  real(kind=dp_t), pointer :: rp(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1

  rp => dataptr(mgts%sync_res(flev), fn)
  res(plo(1):phi(1)) = rp(plo(1):phi(1), 1, 1, 1)

end subroutine mgt_get_sync_res_1d

subroutine mgt_get_sync_res_2d(lev, n, res, plo, phi, lo, hi)
  use nodal_cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(2), hi(2), plo(2), phi(2)
  real(kind=dp_t), intent(inout) :: res(plo(1):phi(1), plo(2):phi(2))
  real(kind=dp_t), pointer :: rp(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1

  rp => dataptr(mgts%sync_res(flev), fn)
  res(plo(1):phi(1), plo(2):phi(2)) = rp(plo(1):phi(1), plo(2):phi(2), 1, 1)

end subroutine mgt_get_sync_res_2d

subroutine mgt_get_sync_res_3d(lev, n, res, plo, phi, lo, hi)
  use nodal_cpp_mg_module
  implicit none
  integer, intent(in) :: lev, n, lo(3), hi(3), plo(3), phi(3)
  real(kind=dp_t), intent(inout) :: res(plo(1):phi(1), plo(2):phi(2), plo(3):phi(3))
  real(kind=dp_t), pointer :: rp(:,:,:,:)
  integer :: flev, fn
  fn = n + 1
  flev = lev+1

  rp => dataptr(mgts%sync_res(flev), fn)
  res(plo(1):phi(1), plo(2):phi(2), plo(3):phi(3)) =  &
       rp(plo(1):phi(1), plo(2):phi(2), plo(3):phi(3), 1)

end subroutine mgt_get_sync_res_3d

subroutine mgt_compute_sync_resid_crse()
  use nodal_cpp_mg_module
  use nodal_sync_resid_module
  use itsol_module, only : itsol_stencil_apply
  implicit none

  integer :: i, dm, mglev
  real(kind=dp_t), pointer :: res(:,:,:,:) 
  real(kind=dp_t), pointer :: msk(:,:,:,:) 
  real(kind=dp_t), pointer :: vo(:,:,:,:) 
  real(kind=dp_t), pointer :: uu(:,:,:,:) 
  real(kind=dp_t), pointer :: pc(:,:,:,:) 
  type(multifab) :: phi_crse
  logical, dimension(3) :: nodal

  dm = get_dim(mgts%sync_res(1))
  mglev = mgts%mgt(1)%nlevels

  nodal = .true.

  call build(phi_crse, mgts%mla%la(1), nc=1, ng=1, nodal=nodal)
  call setval(phi_crse,ZERO,all=.true.)

  do i = 1, nboxes(phi_crse)
     if (remote(phi_crse, i)) cycle
     pc  => dataptr(phi_crse, i)
     msk => dataptr(mgts%sync_msk(1), i)
     uu  => dataptr(mgts%uu      (1), i)
     select case (dm)
     case (1)
        call bl_error('mgt_compute_sync_resid_crse: 1d not done')
     case (2)
        call get_phi_crse_2d(pc(:,:,1,1), msk(:,:,1,1), uu(:,:,1,1))
     case (3)
        call bl_error('mgt_compute_sync_resid_crse: 3d not done')
     end select
  end do

  call multifab_fill_boundary(phi_crse)

  call itsol_stencil_apply(mgts%mgt(1)%ss(mglev), mgts%sync_res(1), phi_crse, &
       mgts%mgt(1)%mm(mglev), mgts%mgt(1)%uniform_dh)

  do i = 1, nboxes(mgts%sync_res(1))
     if (remote(mgts%sync_res(1), i)) cycle
     res => dataptr(mgts%sync_res(1), i)
     msk => dataptr(mgts%sync_msk(1), i)
     vo  => dataptr(mgts%vold    (1), i)
     select case (dm)
     case (1)
        call bl_error('mgt_compute_sync_resid_crse: 1d not done')
     case (2)
        call add_Dv_crse_2d(res(:,:,1,1), msk(:,:,1,1), vo(:,:,1,:), mgts%mgt(1)%dh(:,mglev))
     case (3)
        call bl_error('mgt_compute_sync_resid_crse: 3d not done')
     end select
  end do

  call destroy(phi_crse)

end subroutine mgt_compute_sync_resid_crse

subroutine mgt_compute_sync_resid_fine()
  use nodal_cpp_mg_module
  use nodal_sync_resid_module
  use itsol_module, only : itsol_stencil_apply
  implicit none

  integer :: i, dm, mglev
  real(kind=dp_t), pointer :: res(:,:,:,:) 
  real(kind=dp_t), pointer :: msk(:,:,:,:) 
  real(kind=dp_t), pointer :: vo(:,:,:,:) 
  real(kind=dp_t), pointer :: uu(:,:,:,:) 
  real(kind=dp_t), pointer :: pf(:,:,:,:) 
  type(multifab) :: phi_fine
  logical, dimension(3) :: nodal

  dm = get_dim(mgts%sync_res(1))
  mglev = mgts%mgt(1)%nlevels

  nodal = .true.

  call build(phi_fine, mgts%mla%la(1), nc=1, ng=1, nodal=nodal)
  call setval(phi_fine,ZERO,all=.true.)

  do i = 1, nboxes(phi_fine)
     if (remote(phi_fine, i)) cycle
     pf  => dataptr(phi_fine, i)
     msk => dataptr(mgts%sync_msk(1), i)
     uu  => dataptr(mgts%uu      (1), i)
     select case (dm)
     case (1)
        call bl_error('mgt_compute_sync_resid_fine: 1d not done')
     case (2)
        call get_phi_fine_2d(pf(:,:,1,1), msk(:,:,1,1), uu(:,:,1,1))
     case (3)
        call bl_error('mgt_compute_sync_resid_fine: 3d not done')
     end select
  end do

  call multifab_fill_boundary(phi_fine)

  call itsol_stencil_apply(mgts%mgt(1)%ss(mglev), mgts%sync_res(1), phi_fine, &
       mgts%mgt(1)%mm(mglev), mgts%mgt(1)%uniform_dh)

  do i = 1, nboxes(mgts%sync_res(1))
     if (remote(mgts%sync_res(1), i)) cycle
     res => dataptr(mgts%sync_res(1), i)
     msk => dataptr(mgts%sync_msk(1), i)
     vo  => dataptr(mgts%vel     (1), i)
     select case (dm)
     case (1)
        call bl_error('mgt_compute_sync_resid_fine: 1d not done')
     case (2)
        call add_Dv_fine_2d(res(:,:,1,1), msk(:,:,1,1), vo(:,:,1,:), mgts%mgt(1)%dh(:,mglev))
     case (3)
        call bl_error('mgt_compute_sync_resid_fine: 3d not done')
     end select
  end do

  call destroy(phi_fine)

end subroutine mgt_compute_sync_resid_fine
