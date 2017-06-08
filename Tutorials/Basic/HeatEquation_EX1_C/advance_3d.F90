module advance_module

  implicit none

contains

#ifdef CUDA
  attributes(device) &
#endif
  subroutine compute_flux_doit (lo, hi, phi, philo, phihi, &
       fluxx, fxlo, fxhi, fluxy, fylo, fyhi, fluxz, fzlo, fzhi, &
       dx)

    use amrex_fort_module, only : amrex_real
    implicit none

    integer lo(3), hi(3), philo(3), phihi(3), fxlo(3), fxhi(3), fylo(3), fyhi(3), fzlo(3), fzhi(3)
    real(amrex_real), intent(in)    :: phi  (philo(1):phihi(1),philo(2):phihi(2),philo(3):phihi(3))
    real(amrex_real), intent(inout) :: fluxx( fxlo(1): fxhi(1), fxlo(2): fxhi(2), fxlo(3): fxhi(3))
    real(amrex_real), intent(inout) :: fluxy( fylo(1): fyhi(1), fylo(2): fyhi(2), fylo(3): fyhi(3))
    real(amrex_real), intent(inout) :: fluxz( fzlo(1): fzhi(1), fzlo(2): fzhi(2), fzlo(3): fzhi(3))
    real(amrex_real), intent(in)    :: dx(3)
  
    ! local variables
    integer i,j,k

    ! x-fluxes
    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)+1
             fluxx(i,j,k) = ( phi(i,j,k) - phi(i-1,j,k) ) / dx(1)
          end do
       end do
    end do

    ! y-fluxes
    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)+1
          do i = lo(1), hi(1)
             fluxy(i,j,k) = ( phi(i,j,k) - phi(i,j-1,k) ) / dx(2)
          end do
       end do
    end do

    ! z-fluxes
    do       k = lo(3), hi(3)+1
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             fluxz(i,j,k) = ( phi(i,j,k) - phi(i,j,k-1) ) / dx(3)
          end do
       end do
    end do

  end subroutine compute_flux_doit



#ifdef CUDA
  attributes(device) &
#endif
  subroutine update_phi_doit (lo, hi, phiold, polo, pohi, phinew, pnlo, pnhi, &
       fluxx, fxlo, fxhi, fluxy, fylo, fyhi, fluxz, fzlo, fzhi, &
       dx, dt)

    use amrex_fort_module, only : amrex_real
    implicit none

    integer lo(3), hi(3), polo(3), pohi(3), pnlo(3), pnhi(3), &
         fxlo(3), fxhi(3), fylo(3), fyhi(3), fzlo(3), fzhi(3)
    real(amrex_real), intent(in)    :: phiold(polo(1):pohi(1),polo(2):pohi(2),polo(3):pohi(3))
    real(amrex_real), intent(inout) :: phinew(pnlo(1):pnhi(1),pnlo(2):pnhi(2),pnlo(3):pnhi(3))
    real(amrex_real), intent(in   ) :: fluxx (fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3))
    real(amrex_real), intent(in   ) :: fluxy (fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3))
    real(amrex_real), intent(in   ) :: fluxz (fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3))
    real(amrex_real), intent(in)    :: dx(3)
    real(amrex_real), value         :: dt

    ! local variables
    integer i,j,k
    real(amrex_real) :: dtdx(3)

    dtdx = dt/dx

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)

             phinew(i,j,k) = phiold(i,j,k) &
                  + dtdx(1) * (fluxx(i+1,j  ,k  ) - fluxx(i,j,k)) &
                  + dtdx(2) * (fluxy(i  ,j+1,k  ) - fluxy(i,j,k)) &
                  + dtdx(3) * (fluxz(i  ,j  ,k+1) - fluxz(i,j,k))

          end do
       end do
    end do

  end subroutine update_phi_doit

end module advance_module



subroutine compute_flux (lo, hi, phi, philo, phihi, &
     fluxx, fxlo, fxhi, fluxy, fylo, fyhi, fluxz, fzlo, fzhi, &
     dx, idx) bind(C, name="compute_flux")

  use amrex_fort_module, only : amrex_real
  use advance_module, only: compute_flux_doit
#ifdef CUDA
  use cuda_module, only: cuda_streams, stream_from_index, threads_and_blocks
  use cudafor, only: cudaMemcpyAsync, cudaMemcpyHostToDevice, cudaStreamSynchronize, &
                     cuda_stream_kind, dim3
#endif


  implicit none

  integer lo(3), hi(3), philo(3), phihi(3), fxlo(3), fxhi(3), fylo(3), fyhi(3), fzlo(3), fzhi(3), idx
  real(amrex_real), intent(in)    :: phi  (philo(1):phihi(1),philo(2):phihi(2),philo(3):phihi(3))
  real(amrex_real), intent(inout) :: fluxx( fxlo(1): fxhi(1), fxlo(2): fxhi(2), fxlo(3): fxhi(3))
  real(amrex_real), intent(inout) :: fluxy( fylo(1): fyhi(1), fylo(2): fyhi(2), fylo(3): fyhi(3))
  real(amrex_real), intent(inout) :: fluxz( fzlo(1): fzhi(1), fzlo(2): fzhi(2), fzlo(3): fzhi(3))
  real(amrex_real), intent(in)    :: dx(3)

#ifdef CUDA
  attributes(device) :: phi, fluxx, fluxy, fluxz

  integer :: cuda_result
  integer(kind=cuda_stream_kind) :: stream
  type(dim3) :: numThreads, numBlocks

  integer, device :: lo_d(3), hi_d(3)
  integer, device :: philo_d(3), phihi_d(3) 
  integer, device :: fxlo_d(3), fxhi_d(3) 
  integer, device :: fylo_d(3), fyhi_d(3)
  integer, device :: fzlo_d(3), fzhi_d(3)
  real(amrex_real), device :: dx_d(3)

  stream = cuda_streams(stream_from_index(idx)+1) 

  cuda_result = cudaMemcpyAsync(lo_d, lo, 3, cudaMemcpyHostToDevice, stream)
  cuda_result = cudaMemcpyAsync(hi_d, hi, 3, cudaMemcpyHostToDevice, stream)

  cuda_result = cudaMemcpyAsync(philo_d, philo, 3, cudaMemcpyHostToDevice, stream)
  cuda_result = cudaMemcpyAsync(phihi_d, phihi, 3, cudaMemcpyHostToDevice, stream)

  cuda_result = cudaMemcpyAsync(fxlo_d, fxlo, 3, cudaMemcpyHostToDevice, stream)
  cuda_result = cudaMemcpyAsync(fxhi_d, fxhi, 3, cudaMemcpyHostToDevice, stream)

  cuda_result = cudaMemcpyAsync(fylo_d, fylo, 3, cudaMemcpyHostToDevice, stream)
  cuda_result = cudaMemcpyAsync(fyhi_d, fyhi, 3, cudaMemcpyHostToDevice, stream)

  cuda_result = cudaMemcpyAsync(fzlo_d, fzlo, 3, cudaMemcpyHostToDevice, stream)
  cuda_result = cudaMemcpyAsync(fzhi_d, fzhi, 3, cudaMemcpyHostToDevice, stream)

  cuda_result = cudaMemcpyAsync(dx_d, dx, 3, cudaMemcpyHostToDevice, stream)

  call threads_and_blocks(lo, hi, numBlocks, numThreads)

  call cuda_compute_flux<<<numBlocks, numThreads, 0, stream>>>(lo_d, hi_d, phi, philo_d, phihi_d, &
                                                               fluxx, fxlo_d, fxhi_d, &
                                                               fluxy, fylo_d, fyhi_d, &
                                                               fluxz, fzlo_d, fzhi_d, dx_d)

#else

  call compute_flux_doit(lo, hi, phi, philo, phihi, fluxx, fxlo, fxhi, fluxy, fylo, fyhi, &
                          fzlo, fzhi, dx)

#endif

end subroutine compute_flux


#ifdef CUDA
attributes(global) &
subroutine cuda_compute_flux (lo, hi, phi, philo, phihi, &
     fluxx, fxlo, fxhi, fluxy, fylo, fyhi, fluxz, fzlo, fzhi, dx)

  use amrex_fort_module, only : amrex_real
  use advance_module, only: compute_flux_doit

  implicit none

  integer lo(3), hi(3), philo(3), phihi(3), fxlo(3), fxhi(3), fylo(3), fyhi(3), fzlo(3), fzhi(3)
  real(amrex_real), intent(in)    :: phi  (philo(1):phihi(1),philo(2):phihi(2),philo(3):phihi(3))
  real(amrex_real), intent(inout) :: fluxx( fxlo(1): fxhi(1), fxlo(2): fxhi(2), fxlo(3): fxhi(3))
  real(amrex_real), intent(inout) :: fluxy( fylo(1): fyhi(1), fylo(2): fyhi(2), fylo(3): fyhi(3))
  real(amrex_real), intent(inout) :: fluxz( fzlo(1): fzhi(1), fzlo(2): fzhi(2), fzlo(3): fzhi(3))
  real(amrex_real), intent(in)    :: dx(3)

  integer :: idx(3)

  ! Get our spatial index based on the CUDA thread index

  idx(1) = lo(1) + (threadIdx%x - 1) + blockDim%x * (blockIdx%x - 1)
  idx(2) = lo(2) + (threadIdx%y - 1) + blockDim%y * (blockIdx%y - 1)
  idx(3) = lo(3) + (threadIdx%z - 1) + blockDim%z * (blockIdx%z - 1)

  if (idx(1) .gt. hi(1) .or. idx(2) .gt. hi(2) .or. idx(3) .gt. hi(3)) return

  call compute_flux_doit(idx, idx, phi, philo, phihi, fluxx, fxlo, fxhi, fluxy, fylo, fyhi, &
                         fluxz, fzlo, fzhi, dx)

end subroutine cuda_compute_flux
#endif


subroutine update_phi (lo, hi, phiold, polo, pohi, phinew, pnlo, pnhi, &
     fluxx, fxlo, fxhi, fluxy, fylo, fyhi, fluxz, fzlo, fzhi, &
     dx, dt, idx) bind(C, name="update_phi")

  use amrex_fort_module, only : amrex_real
  use advance_module, only: update_phi_doit
#ifdef CUDA
  use cuda_module, only: cuda_streams, stream_from_index, threads_and_blocks
  use cudafor, only: cudaMemcpyAsync, cudaMemcpyHostToDevice, cudaStreamSynchronize, &
                     cuda_stream_kind, dim3
#endif

  implicit none

    integer lo(3), hi(3), polo(3), pohi(3), pnlo(3), pnhi(3), &
       fxlo(3), fxhi(3), fylo(3), fyhi(3), fzlo(3), fzhi(3), idx
  real(amrex_real), intent(in)    :: phiold(polo(1):pohi(1),polo(2):pohi(2),polo(3):pohi(3))
  real(amrex_real), intent(inout) :: phinew(pnlo(1):pnhi(1),pnlo(2):pnhi(2),pnlo(3):pnhi(3))
  real(amrex_real), intent(in   ) :: fluxx (fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3))
  real(amrex_real), intent(in   ) :: fluxy (fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3))
  real(amrex_real), intent(in   ) :: fluxz (fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3))
  real(amrex_real), intent(in)    :: dx(3)
  real(amrex_real), value         :: dt

#ifdef CUDA
  attributes(device) :: phi, fluxx, fluxy, fluxz

  integer :: cuda_result
  integer(kind=cuda_stream_kind) :: stream
  type(dim3) :: numThreads, numBlocks

  integer, device :: lo_d(3), hi_d(3)
  integer, device :: polo_d(3), pohi_d(3)
  integer, device :: pnlo_d(3), pnhi_d(3)
  integer, device :: fxlo_d(3), fxhi_d(3)
  integer, device :: fylo_d(3), fyhi_d(3)
  integer, device :: fzlo_d(3), fzhi_d(3)
  real(amrex_real), device :: dx_d(3), dt_d

  stream = cuda_streams(stream_from_index(idx)+1) 

  cuda_result = cudaMemcpyAsync(lo_d, lo, 3, cudaMemcpyHostToDevice, stream)
  cuda_result = cudaMemcpyAsync(hi_d, hi, 3, cudaMemcpyHostToDevice, stream)

  cuda_result = cudaMemcpyAsync(polo_d, polo, 3, cudaMemcpyHostToDevice, stream)
  cuda_result = cudaMemcpyAsync(pohi_d, pohi, 3, cudaMemcpyHostToDevice, stream)

  cuda_result = cudaMemcpyAsync(pnlo_d, pnlo, 3, cudaMemcpyHostToDevice, stream)
  cuda_result = cudaMemcpyAsync(pnhi_d, pnhi, 3, cudaMemcpyHostToDevice, stream)

  cuda_result = cudaMemcpyAsync(fxlo_d, fxlo, 3, cudaMemcpyHostToDevice, stream)
  cuda_result = cudaMemcpyAsync(fxhi_d, fxhi, 3, cudaMemcpyHostToDevice, stream)

  cuda_result = cudaMemcpyAsync(fylo_d, fylo, 3, cudaMemcpyHostToDevice, stream)
  cuda_result = cudaMemcpyAsync(fyhi_d, fyhi, 3, cudaMemcpyHostToDevice, stream)

  cuda_result = cudaMemcpyAsync(fzlo_d, fzlo, 3, cudaMemcpyHostToDevice, stream)
  cuda_result = cudaMemcpyAsync(fzhi_d, fzhi, 3, cudaMemcpyHostToDevice, stream)

  cuda_result = cudaMemcpyAsync(dx_d, dx, 3, cudaMemcpyHostToDevice, stream)

  cuda_result = cudaMemcpyAsync(dt_d, dt, 1, cudaMemcpyHostToDevice, stream)

  call threads_and_blocks(lo, hi, numBlocks, numThreads)

  call cuda_update_phi<<<numBlocks, numThreads, 0, stream>>>(lo_d, hi_d, phiold, polo_d, pohi_d, &
                                                             phinew, pnlo_d, pnhi_d, fluxx, fxlo_d, fxhi_d, &
                                                             fluxy, fylo_d, fyhi_d, &
                                                             fluxz, fzlo_d, fzhi_d, &
                                                             dx_d, dt_d)
    
#else

  call update_phi_doit(idx, idx, phiold, polo, pohi, phinew, pnlo, pnhi, &
                       fluxx, fxlo, fxhi, fluxy, fylo, fyhi, fluxz, fzlo, fzhi, &
                       dx, dt)
#endif 

  end subroutine update_phi


#ifdef CUDA
attributes(global) &
subroutine cuda_update_phi (lo, hi, phiold, polo, pohi, phinew, pnlo, pnhi, &
     fluxx, fxlo, fxhi, fluxy, fylo, fyhi, fluxz, fzlo, fzhi, &
     dx, dt)

  use amrex_fort_module, only : amrex_real
  use advance_module, only: update_phi_doit

  implicit none

    integer lo(3), hi(3), polo(3), pohi(3), pnlo(3), pnhi(3), &
       fxlo(3), fxhi(3), fylo(3), fyhi(3), fzlo(3), fzhi(3)
  real(amrex_real), intent(in)    :: phiold(polo(1):pohi(1),polo(2):pohi(2),polo(3):pohi(3))
  real(amrex_real), intent(inout) :: phinew(pnlo(1):pnhi(1),pnlo(2):pnhi(2),pnlo(3):pnhi(3))
  real(amrex_real), intent(in   ) :: fluxx (fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3))
  real(amrex_real), intent(in   ) :: fluxy (fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3))
  real(amrex_real), intent(in   ) :: fluxz (fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3))
  real(amrex_real), intent(in)    :: dx(3)
  real(amrex_real), value         :: dt

  integer :: idx(3)

  ! Get our spatial index based on the CUDA thread index

  idx(1) = lo(1) + (threadIdx%x - 1) + blockDim%x * (blockIdx%x - 1)
  idx(2) = lo(2) + (threadIdx%y - 1) + blockDim%y * (blockIdx%y - 1)
  idx(3) = lo(3) + (threadIdx%z - 1) + blockDim%z * (blockIdx%z - 1)

  if (idx(1) .gt. hi(1) .or. idx(2) .gt. hi(2) .or. idx(3) .gt. hi(3)) return


  call update_phi_doit(idx, idx, phiold, polo, pohi, phinew, pnlo, pnhi, &
                       fluxx, fxlo, fxhi, fluxy, fylo, fyhi, fluxz, fzlo, fzhi, &
                       dx, dt)

  end subroutine cuda_update_phi
#endif
