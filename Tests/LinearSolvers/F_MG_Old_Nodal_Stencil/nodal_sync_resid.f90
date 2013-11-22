module nodal_sync_resid_module
  use bl_constants_module
  use bc_functions_module
  use multifab_module
  implicit none

contains

  subroutine compute_divuo(divuo, mask, vold, dx, face_type)
    type(multifab) , intent(inout) :: divuo
    type(multifab) , intent(in   ) :: mask
    type(multifab) , intent(in   ) :: vold
    real(kind=dp_t), intent(in   ) :: dx(:)
    integer        , intent(in   ) :: face_type(:,:,:)

    integer :: i, dm
    real(kind=dp_t), pointer :: msk(:,:,:,:) 
    real(kind=dp_t), pointer :: vo(:,:,:,:) 
    real(kind=dp_t), pointer :: dvo(:,:,:,:) 
    
    dm = get_dim(vold)

    do i = 1, nfabs(divuo)
       dvo => dataptr(divuo, i)
       msk => dataptr(mask , i)
       vo  => dataptr(vold , i)
       select case (dm)
       case (1)
          call bl_error('divuo: 1d not done')
       case (2)
          call divuo_2d(dvo(:,:,1,1), msk(:,:,1,1), vo(:,:,1,:), dx, face_type(i,:,:))
       case (3)
          call divuo_3d(dvo(:,:,:,1), msk(:,:,:,1), vo(:,:,:,:), dx, face_type(i,:,:))
       end select
    end do

  end subroutine compute_divuo

  subroutine divuo_2d(dvo, msk, vel, dx, face_type)
    real(kind=dp_t), intent(inout) :: dvo(-1:,-1:)
    real(kind=dp_t), intent(in   ) :: msk(-1:,-1:)
    real(kind=dp_t), intent(in   ) :: vel(-1:,-1:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    integer        , intent(in   ) :: face_type(:,:)
    
    integer :: i, j, nx, ny
    real(kind=dp_t) :: divv

    nx = size(msk,dim=1) - 2
    ny = size(msk,dim=2) - 2

    do j = 0, ny
       do i = 0, nx
          divv = (vel(i  ,j,1)*msk(i , j) + vel(i  ,j-1,1)*msk(i  ,j-1) &
               -  vel(i-1,j,1)*msk(i-1,j) - vel(i-1,j-1,1)*msk(i-1,j-1)) / dx(1) &
               + (vel(i,j  ,2)*msk(i,j  ) + vel(i-1,j  ,2)*msk(i-1,j  ) &
               -  vel(i,j-1,2)*msk(i,j-1) - vel(i-1,j-1,2)*msk(i-1,j-1)) / dx(2)
          dvo(i,j) = HALF * divv
       end do
    end do

    if (face_type(1,1) == BC_NEU) dvo( 0,:) = TWO*dvo( 0,:)
    if (face_type(1,2) == BC_NEU) dvo(nx,:) = TWO*dvo(nx,:)
    if (face_type(2,1) == BC_NEU) dvo(:, 0) = TWO*dvo(:, 0)
    if (face_type(2,2) == BC_NEU) dvo(:,ny) = TWO*dvo(:,ny)

  end subroutine divuo_2d

  subroutine divuo_3d(dvo, msk, vel, dx, face_type)
    real(kind=dp_t), intent(inout) :: dvo(-1:,-1:,-1:)
    real(kind=dp_t), intent(in   ) :: msk(-1:,-1:,-1:)
    real(kind=dp_t), intent(in   ) :: vel(-1:,-1:,-1:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    integer        , intent(in   ) :: face_type(:,:)
    
    integer         :: i, j, k, nx, ny, nz
    real(kind=dp_t) :: ivdx,ivdy,ivdz

    nx = size(msk,dim=1) - 2
    ny = size(msk,dim=2) - 2
    nz = size(msk,dim=3) - 2

    ivdx = 1.0d0 / dx(1)
    ivdy = 1.0d0 / dx(2)
    ivdz = 1.0d0 / dx(3)

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = 0, nz
    do j = 0, ny
       do i = 0, nx
          dvo(i,j,k) = FOURTH * (  &
               & (vel(i  ,j  ,k  ,1) * msk(i  ,j  ,k  ) &
               +  vel(i  ,j-1,k  ,1) * msk(i  ,j-1,k  ) &
               +  vel(i  ,j  ,k-1,1) * msk(i  ,j  ,k-1) &
               +  vel(i  ,j-1,k-1,1) * msk(i  ,j-1,k-1) &
               -  vel(i-1,j  ,k  ,1) * msk(i-1,j  ,k  ) &
               -  vel(i-1,j-1,k  ,1) * msk(i-1,j-1,k  ) &
               -  vel(i-1,j  ,k-1,1) * msk(i-1,j  ,k-1) &
               -  vel(i-1,j-1,k-1,1) * msk(i-1,j-1,k-1)) * ivdx &
               + (vel(i  ,j  ,k  ,2) * msk(i  ,j  ,k  ) &
               +  vel(i-1,j  ,k  ,2) * msk(i-1,j  ,k  ) &
               +  vel(i  ,j  ,k-1,2) * msk(i  ,j  ,k-1) &
               +  vel(i-1,j  ,k-1,2) * msk(i-1,j  ,k-1) &
               -  vel(i  ,j-1,k  ,2) * msk(i  ,j-1,k  ) &
               -  vel(i-1,j-1,k  ,2) * msk(i-1,j-1,k  ) &
               -  vel(i  ,j-1,k-1,2) * msk(i  ,j-1,k-1) & 
               -  vel(i-1,j-1,k-1,2) * msk(i-1,j-1,k-1)) * ivdy &
               + (vel(i  ,j  ,k  ,3) * msk(i  ,j  ,k  ) &
               +  vel(i-1,j  ,k  ,3) * msk(i-1,j  ,k  ) &
               +  vel(i  ,j-1,k  ,3) * msk(i  ,j-1,k  ) & 
               +  vel(i-1,j-1,k  ,3) * msk(i-1,j-1,k  ) &
               -  vel(i  ,j  ,k-1,3) * msk(i  ,j  ,k-1) & 
               -  vel(i-1,j  ,k-1,3) * msk(i-1,j  ,k-1) &
               -  vel(i  ,j-1,k-1,3) * msk(i  ,j-1,k-1) &
               -  vel(i-1,j-1,k-1,3) * msk(i-1,j-1,k-1)) * ivdz )
       end do
    end do
    end do
    !$OMP END PARALLEL DO

    if (face_type(1,1) == BC_NEU) dvo( 0,:,:) = TWO*dvo( 0,:,:)
    if (face_type(1,2) == BC_NEU) dvo(nx,:,:) = TWO*dvo(nx,:,:)
    if (face_type(2,1) == BC_NEU) dvo(:, 0,:) = TWO*dvo(:, 0,:)
    if (face_type(2,2) == BC_NEU) dvo(:,ny,:) = TWO*dvo(:,ny,:)
    if (face_type(3,1) == BC_NEU) dvo(:,:, 0) = TWO*dvo(:,:, 0)
    if (face_type(3,2) == BC_NEU) dvo(:,:,nz) = TWO*dvo(:,:,nz)

  end subroutine divuo_3d

  subroutine comp_sync_res(sync_res, divuo, mask, sign_res)
    type(multifab) , intent(inout) :: sync_res
    type(multifab) , intent(in   ) :: divuo
    type(multifab) , intent(in   ) :: mask
    real(kind=dp_t), intent(in   ) :: sign_res

    integer :: i, dm
    real(kind=dp_t), pointer :: res(:,:,:,:) 
    real(kind=dp_t), pointer :: dvo(:,:,:,:) 
    real(kind=dp_t), pointer :: msk(:,:,:,:) 

    dm = get_dim(sync_res)

    do i = 1, nfabs(sync_res)
       res => dataptr(sync_res, i)
       dvo => dataptr(divuo   , i)
       msk => dataptr(mask    , i)
       select case (dm)
       case (1)
          call bl_error('comp_sync_res: 1d not done')
       case (2)
          call comp_sync_res_2d(res(:,:,1,1), dvo(:,:,1,1), msk(:,:,1,1), sign_res)
       case (3)
          call comp_sync_res_3d(res(:,:,:,1), dvo(:,:,:,1), msk(:,:,:,1), sign_res)
       end select
    end do

  end subroutine comp_sync_res

  subroutine comp_sync_res_2d(res, dvo, msk, sgnr)
    real(kind=dp_t), intent(inout) :: res(-1:,-1:)
    real(kind=dp_t), intent(in   ) :: dvo(-1:,-1:)
    real(kind=dp_t), intent(in   ) :: msk(-1:,-1:)
    real(kind=dp_t), intent(in   ) :: sgnr
    
    integer :: i, j, nx, ny

    nx = size(msk,dim=1) - 2
    ny = size(msk,dim=2) - 2

    do j = 0, ny
       do i = 0, nx
          if ( any(msk(i-1:i,j-1:j) .eq. ONE) .and. &
               any(msk(i-1:i,j-1:j) .eq. ZERO) ) then
             res(i,j) = sgnr*res(i,j) + dvo(i,j)
          else
             res(i,j) = ZERO
          end if
       end do
    end do

  end subroutine comp_sync_res_2d

  subroutine comp_sync_res_3d(res, dvo, msk, sgnr)
    real(kind=dp_t), intent(inout) :: res(-1:,-1:,-1:)
    real(kind=dp_t), intent(in   ) :: dvo(-1:,-1:,-1:)
    real(kind=dp_t), intent(in   ) :: msk(-1:,-1:,-1:)
    real(kind=dp_t), intent(in   ) :: sgnr
    
    integer         :: i, j, k, nx, ny, nz
    real(kind=dp_t) :: result

    nx = size(msk,dim=1) - 2
    ny = size(msk,dim=2) - 2
    nz = size(msk,dim=3) - 2

    !$OMP PARALLEL DO PRIVATE(i,j,k,result)
    do k = 0, nz
       do j = 0, ny
          do i = 0, nx
             result = ZERO
             if ( any(msk(i-1:i,j-1:j,k-1:k) .eq. ONE) ) then
                if ( any(msk(i-1:i,j-1:j,k-1:k) .eq. ZERO) ) then
                   result = sgnr*res(i,j,k) + dvo(i,j,k)
                endif
             endif
             res(i,j,k) = result
          end do
       end do
    end do
    !$OMP END PARALLEL DO

  end subroutine comp_sync_res_3d

  subroutine divuo_add_rhcc(divuo, rhcc, mask, face_type)
    type(multifab) , intent(inout) :: divuo
    type(multifab) , intent(inout) :: rhcc
    type(multifab) , intent(in   ) :: mask
    integer        , intent(in   ) :: face_type(:,:,:)
    
    integer                  :: i, dm
    real(kind=dp_t), pointer :: msk(:,:,:,:) 
    real(kind=dp_t), pointer :: dvo(:,:,:,:) 
    real(kind=dp_t), pointer :: rc(:,:,:,:) 
    
    dm = get_dim(rhcc)

    do i = 1, nfabs(divuo)
       dvo => dataptr(divuo, i)
       msk => dataptr(mask , i)
       rc  => dataptr(rhcc , i)
       select case (dm)
       case (1)
          call bl_error('divuo_rhcc_1d: 1d not done')
       case (2)
          call divuo_rhcc_2d(dvo(:,:,1,1), msk(:,:,1,1), rc(:,:,1,1), face_type(i,:,:))
       case (3)
          call divuo_rhcc_3d(dvo(:,:,:,1), msk(:,:,:,1), rc(:,:,:,1), face_type(i,:,:))
       end select
    end do

  end subroutine divuo_add_rhcc

  subroutine divuo_rhcc_2d(dvo, msk, rc, face_type)
    real(kind=dp_t), intent(inout) :: dvo(-1:,-1:)
    real(kind=dp_t), intent(in   ) :: msk(-1:,-1:)
    real(kind=dp_t), intent(inout) ::  rc(-1:,-1:)
    integer        , intent(in   ) :: face_type(:,:)
    
    integer :: i, j, nx, ny
    real(kind=dp_t), pointer :: tmp(:,:)

    nx = size(msk,dim=1) - 2
    ny = size(msk,dim=2) - 2

    rc(-1,:) = ZERO
    rc(nx,:) = ZERO
    rc(:,-1) = ZERO
    rc(:,ny) = ZERO

    allocate(tmp(0:nx,0:ny))

    do j = 0, ny
       do i = 0, nx
          tmp(i,j) = FOURTH *   &
               ( rc(i-1,j-1) * msk(i-1,j-1)  &
               + rc(i  ,j-1) * msk(i  ,j-1)  &
               + rc(i-1,j  ) * msk(i-1,j  )  &
               + rc(i  ,j  ) * msk(i  ,j  ) )
       end do
    end do

    if (face_type(1,1) == BC_NEU) tmp( 0,:) = TWO*tmp( 0,:)
    if (face_type(1,2) == BC_NEU) tmp(nx,:) = TWO*tmp(nx,:)
    if (face_type(2,1) == BC_NEU) tmp(:, 0) = TWO*tmp(:, 0)
    if (face_type(2,2) == BC_NEU) tmp(:,ny) = TWO*tmp(:,ny)

    do j = 0, ny
       do i = 0, nx
          dvo(i,j) = dvo(i,j) + tmp(i,j)
       end do
    end do

    deallocate(tmp)

  end subroutine divuo_rhcc_2d

  subroutine divuo_rhcc_3d(dvo, msk, rc, face_type)
    real(kind=dp_t), intent(inout) :: dvo(-1:,-1:,-1:)
    real(kind=dp_t), intent(in   ) :: msk(-1:,-1:,-1:)
    real(kind=dp_t), intent(inout) ::  rc(-1:,-1:,-1:)
    integer        , intent(in   ) :: face_type(:,:)
    
    integer :: i, j, k, nx, ny, nz
    real(kind=dp_t), pointer :: tmp(:,:,:)

    nx = size(msk,dim=1) - 2
    ny = size(msk,dim=2) - 2
    nz = size(msk,dim=3) - 2

    rc(-1,:,:) = ZERO
    rc(nx,:,:) = ZERO
    rc(:,-1,:) = ZERO
    rc(:,ny,:) = ZERO
    rc(:,:,-1) = ZERO
    rc(:,:,nz) = ZERO

    allocate(tmp(0:nx,0:ny,0:nz))

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = 0, nz
    do j = 0, ny
    do i = 0, nx
       tmp(i,j,k) = EIGHTH *   &
            ( rc(i-1,j-1,k-1) * msk(i-1,j-1,k-1)  &
            + rc(i  ,j-1,k-1) * msk(i  ,j-1,k-1)  &
            + rc(i-1,j  ,k-1) * msk(i-1,j  ,k-1)  &
            + rc(i  ,j  ,k-1) * msk(i  ,j  ,k-1)  &
            + rc(i-1,j-1,k  ) * msk(i-1,j-1,k  )  &
            + rc(i  ,j-1,k  ) * msk(i  ,j-1,k  )  &
            + rc(i-1,j  ,k  ) * msk(i-1,j  ,k  )  &
            + rc(i  ,j  ,k  ) * msk(i  ,j  ,k  )  )
    end do
    end do
    end do
    !$OMP END PARALLEL DO

    if (face_type(1,1) == BC_NEU) tmp( 0,:,:) = TWO*tmp( 0,:,:)
    if (face_type(1,2) == BC_NEU) tmp(nx,:,:) = TWO*tmp(nx,:,:)
    if (face_type(2,1) == BC_NEU) tmp(:, 0,:) = TWO*tmp(:, 0,:)
    if (face_type(2,2) == BC_NEU) tmp(:,ny,:) = TWO*tmp(:,ny,:)
    if (face_type(3,1) == BC_NEU) tmp(:,:, 0) = TWO*tmp(:,:, 0)
    if (face_type(3,2) == BC_NEU) tmp(:,:,nz) = TWO*tmp(:,:,nz)

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = 0, nz
    do j = 0, ny
    do i = 0, nx
       dvo(i,j,k) = dvo(i,j,k) + tmp(i,j,k)
    end do
    end do
    end do
    !$OMP END PARALLEL DO

    deallocate(tmp)

  end subroutine divuo_rhcc_3d

  subroutine sync_res_fine_bndry(res_fine, face_type)
    type(multifab), intent(inout) :: res_fine
    integer, intent(in) :: face_type(:,:,:)

    integer :: i, dm
    real(kind=dp_t), pointer :: resp(:,:,:,:) 

    dm = get_dim(res_fine)

    do i = 1, nfabs(res_fine)
       resp => dataptr(res_fine, i)
       select case (dm)
       case (1)
          call bl_error('sync_res_fine_bndry: 1d not done')
       case (2)
          call res_fine_bndry_2d(resp(:,:,1,1), face_type(i,:,:))
       case (3)
          call res_fine_bndry_3d(resp(:,:,:,1), face_type(i,:,:))
       end select
    end do

  end subroutine sync_res_fine_bndry

  subroutine res_fine_bndry_2d(res, face_type)
    real(kind=dp_t), intent(inout) :: res(-1:,-1:)
    integer        , intent(in   ) :: face_type(:,:)
    integer :: nx, ny

    nx = size(res,dim=1) - 3
    ny = size(res,dim=2) - 3

    if (face_type(1,1) == BC_NEU) res( 0,:) = HALF*res( 0,:)
    if (face_type(1,2) == BC_NEU) res(nx,:) = HALF*res(nx,:)
    if (face_type(2,1) == BC_NEU) res(:, 0) = HALF*res(:, 0)
    if (face_type(2,2) == BC_NEU) res(:,ny) = HALF*res(:,ny)

  end subroutine res_fine_bndry_2d

  subroutine res_fine_bndry_3d(res, face_type)
    real(kind=dp_t), intent(inout) :: res(-1:,-1:,-1:)
    integer        , intent(in   ) :: face_type(:,:)
    integer :: nx, ny, nz

    nx = size(res,dim=1) - 3
    ny = size(res,dim=2) - 3
    nz = size(res,dim=3) - 3

    if (face_type(1,1) == BC_NEU) res( 0,:,:) = HALF*res( 0,:,:)
    if (face_type(1,2) == BC_NEU) res(nx,:,:) = HALF*res(nx,:,:)
    if (face_type(2,1) == BC_NEU) res(:, 0,:) = HALF*res(:, 0,:)
    if (face_type(2,2) == BC_NEU) res(:,ny,:) = HALF*res(:,ny,:)
    if (face_type(3,1) == BC_NEU) res(:,:, 0) = HALF*res(:,:, 0)
    if (face_type(3,2) == BC_NEU) res(:,:,nz) = HALF*res(:,:,nz)

  end subroutine res_fine_bndry_3d

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
! wqz. unnecessary to call setval(mgts%vold(1),ZERO,all=.true.)
  
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

  mskp => dataptr(mgts%sync_msk(flev), local_index(mgts%sync_msk(flev),fn))
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

  mskp => dataptr(mgts%sync_msk(flev), local_index(mgts%sync_msk(flev),fn))
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

  mskp => dataptr(mgts%sync_msk(flev), local_index(mgts%sync_msk(flev),fn))
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

  vp => dataptr(mgts%vold(flev), local_index(mgts%vold(flev),fn))
  vp(lo(1)-1:hi(1)+1,1,1,1) = v_in(lo(1)-1:hi(1)+1)
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

  vp => dataptr(mgts%vold(flev), local_index(mgts%vold(flev),fn))
  vp(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,1,1:2) =   &
       v_in(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,1:2)
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

  vp => dataptr(mgts%vold(flev), local_index(mgts%vold(flev),fn))
  vp(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,1:3) = &
       v_in(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,1:3)
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

  rp => dataptr(mgts%sync_res(flev), local_index(mgts%sync_res(flev),fn))
  res(lo(1):hi(1)) = rp(lo(1):hi(1), 1, 1, 1)

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

  rp => dataptr(mgts%sync_res(flev), local_index(mgts%sync_res(flev),fn))
  res(lo(1):hi(1), lo(2):hi(2)) = rp(lo(1):hi(1), lo(2):hi(2), 1, 1)

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

  rp => dataptr(mgts%sync_res(flev), local_index(mgts%sync_res(flev),fn))
  res(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)) =  &
       rp(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 1)

end subroutine mgt_get_sync_res_3d

subroutine mgt_compute_sync_resid_crse()
  use nodal_cpp_mg_module
  use nodal_sync_resid_module
  use nodal_stencil_fill_module, only : stencil_fill_nodal
  use stencil_defect_module, only : stencil_apply
  implicit none

  integer :: dm, mglev
  real(kind=dp_t) :: sign_res
  type(multifab) :: divuo
  logical :: nodal(3)

  nodal = .true.
  dm = get_dim(mgts%sync_res(1))
  mglev = mgts%mgt(1)%nlevels

  call build(divuo, mgts%mla%la(1), nc=1, ng=1, nodal=nodal)

  call compute_divuo(divuo, mgts%sync_msk(1), mgts%vold(1), mgts%mgt(1)%dh(:,mglev), &
       mgts%mgt(1)%face_type)

  if (associated(mgts%rhcc)) then  ! only single level solve could get here
     call divuo_add_rhcc(divuo, mgts%rhcc(1), mgts%sync_msk(1), mgts%mgt(1)%face_type)
  end if

  call multifab_mult_mult(mgts%amr_coeffs(1), mgts%sync_msk(1))

  call stencil_fill_nodal(mgts%mgt(1)%ss(mglev), mgts%amr_coeffs(1), &
       mgts%mgt(1)%dh(:,mglev), mgts%mgt(1)%mm(mglev), &
       mgts%mgt(1)%face_type, mgts%stencil_type)

  call stencil_apply(mgts%mgt(1)%ss(mglev), mgts%sync_res(1), mgts%uu(1), &
                     mgts%mgt(1)%mm(mglev), mgts%mgt(1)%stencil_type, &
                     mgts%mgt(1)%lcross, mgts%mgt(1)%uniform_dh)

  sign_res = -ONE
  call comp_sync_res(mgts%sync_res(1), divuo, mgts%sync_msk(1), sign_res)

  call destroy(divuo)

end subroutine mgt_compute_sync_resid_crse

subroutine mgt_compute_sync_resid_fine()
  use nodal_cpp_mg_module
  use nodal_sync_resid_module
  use ml_nd_module
  use nodal_stencil_fill_module, only : stencil_fill_nodal, stencil_fill_one_sided
  implicit none

  integer :: dm, mglev
  real(kind=dp_t) :: sign_res
  type(multifab) :: ss1
  type(multifab) :: divuo
  type(multifab) :: rh0  
  logical :: nodal(3)

  nodal = .true.
  dm = get_dim(mgts%sync_res(1))
  mglev = mgts%mgt(1)%nlevels

  call build(divuo, mgts%mla%la(1), nc=1, ng=1, nodal=nodal)

  call build(rh0, mgts%mla%la(1), nc=1, ng=1, nodal=nodal)
  call setval(rh0,ZERO,all=.true.)

  call compute_divuo(divuo, mgts%sync_msk(1), mgts%vold(1), mgts%mgt(1)%dh(:,mglev), &
       mgts%mgt(1)%face_type)

  if (associated(mgts%rhcc)) then  
     call divuo_add_rhcc(divuo, mgts%rhcc(1), mgts%sync_msk(1), mgts%mgt(1)%face_type)
  end if

  if (mgts%stencil_type .eq. ND_CROSS_STENCIL) then
     call multifab_build(ss1, mgts%mla%la(1), 2*dm+1, 0, nodal, stencil=.true.)
     call stencil_fill_one_sided(ss1, mgts%amr_coeffs(1), mgts%mgt(1)%dh(:,mglev), &
          mgts%mgt(1)%mm(mglev), mgts%mgt(1)%face_type)

     call grid_res(ss1, &
          mgts%sync_res(1), rh0, mgts%uu(1), mgts%mgt(1)%mm(mglev), &
          mgts%mgt(1)%face_type, &
          mgts%mgt(1)%lcross, mgts%mgt(1)%uniform_dh)
  else
     call grid_res(mgts%mgt(1)%ss(mglev), &
          mgts%sync_res(1), rh0, mgts%uu(1), mgts%mgt(1)%mm(mglev), &
          mgts%mgt(1)%face_type, &
          mgts%mgt(1)%lcross, mgts%mgt(1)%uniform_dh)
  endif

  sign_res = ONE
  call comp_sync_res(mgts%sync_res(1), divuo, mgts%sync_msk(1), sign_res)

  call sync_res_fine_bndry(mgts%sync_res(1), mgts%mgt(1)%face_type)

  if (mgts%stencil_type .eq. ND_CROSS_STENCIL) then
     call destroy(ss1)
  endif

  call destroy(divuo)
  call destroy(rh0)

end subroutine mgt_compute_sync_resid_fine

