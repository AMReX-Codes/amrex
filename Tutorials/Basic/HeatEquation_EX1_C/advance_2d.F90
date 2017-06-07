module advance_module

  implicit none

contains

#ifdef CUDA
  attributes(device) &
#endif
  subroutine compute_flux_doit (lo, hi, phi, philo, phihi, &
       fluxx, fxlo, fxhi, fluxy, fylo, fyhi, dx)

    use amrex_fort_module, only : amrex_real
    implicit none

    integer lo(2), hi(2), philo(2), phihi(2), fxlo(2), fxhi(2), fylo(2), fyhi(2)
    real(amrex_real), intent(in)    :: phi  (philo(1):phihi(1),philo(2):phihi(2))
    real(amrex_real), intent(inout) :: fluxx( fxlo(1): fxhi(1), fxlo(2): fxhi(2))
    real(amrex_real), intent(inout) :: fluxy( fylo(1): fyhi(1), fylo(2): fyhi(2))
    real(amrex_real), intent(in)    :: dx(2)

    ! local variables
    integer i,j

    ! x-fluxes
    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)+1
          fluxx(i,j) = ( phi(i,j) - phi(i-1,j) ) / dx(1)
       end do
    end do

    ! y-fluxes
    do    j = lo(2), hi(2)+1
       do i = lo(1), hi(1)
          fluxy(i,j) = ( phi(i,j) - phi(i,j-1) ) / dx(2)
       end do
    end do

  end subroutine compute_flux_doit


  subroutine update_phi (lo, hi, phiold, polo, pohi, phinew, pnlo, pnhi, &
       fluxx, fxlo, fxhi, fluxy, fylo, fyhi, dx, dt) bind(C, name="update_phi")

    use amrex_fort_module, only : amrex_real
    implicit none

    integer lo(2), hi(2), polo(2), pohi(2), pnlo(2), pnhi(2), fxlo(2), fxhi(2), fylo(2), fyhi(2)
    real(amrex_real), intent(in)    :: phiold(polo(1):pohi(1),polo(2):pohi(2))
    real(amrex_real), intent(inout) :: phinew(pnlo(1):pnhi(1),pnlo(2):pnhi(2))
    real(amrex_real), intent(in   ) :: fluxx (fxlo(1):fxhi(1),fxlo(2):fxhi(2))
    real(amrex_real), intent(in   ) :: fluxy (fylo(1):fyhi(1),fylo(2):fyhi(2))
    real(amrex_real), intent(in)    :: dx(2)
    real(amrex_real), value         :: dt

    ! local variables
    integer i,j
    real(amrex_real) :: dtdx(2)

    dtdx = dt/dx

    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)

          phinew(i,j) = phiold(i,j) &
               + dtdx(1) * (fluxx(i+1,j  ) - fluxx(i,j)) &
               + dtdx(2) * (fluxy(i  ,j+1) - fluxy(i,j))

       end do
    end do

  end subroutine update_phi

end module advance_module



subroutine compute_flux (lo, hi, phi, philo, phihi, &
     fluxx, fxlo, fxhi, fluxy, fylo, fyhi, dx, idx) bind(C, name="compute_flux")

  use amrex_fort_module, only : amrex_real
  use advance_module, only: compute_flux_doit
#ifdef CUDA
  use cuda_module, only: cuda_streams, stream_from_index, threads_and_blocks
  use cudafor, only: cudaMemcpyAsync, cudaMemcpyHostToDevice, cudaStreamSynchronize, &
                     cuda_stream_kind, dim3
#endif

  implicit none

  integer lo(2), hi(2), philo(2), phihi(2), fxlo(2), fxhi(2), fylo(2), fyhi(2), idx
  real(amrex_real), intent(in)    :: phi  (philo(1):phihi(1),philo(2):phihi(2))
  real(amrex_real), intent(inout) :: fluxx( fxlo(1): fxhi(1), fxlo(2): fxhi(2))
  real(amrex_real), intent(inout) :: fluxy( fylo(1): fyhi(1), fylo(2): fyhi(2))
  real(amrex_real), intent(in)    :: dx(2)

#ifdef CUDA
  attributes(device) :: phi, fluxx, fluxy

  integer :: cuda_result
  integer(kind=cuda_stream_kind) :: stream
  type(dim3) :: numThreads, numBlocks

  integer, device :: lo_d(2), hi_d(2)
  integer, device :: philo_d(2), phihi_d(2)
  integer, device :: fxlo_d(2), fxhi_d(2)
  integer, device :: fylo_d(2), fyhi_d(2)
  real(amrex_real), device :: dx_d(2)

  stream = cuda_streams(stream_from_index(idx)+1)

  cuda_result = cudaMemcpyAsync(lo_d, lo, 2, cudaMemcpyHostToDevice, stream)
  cuda_result = cudaMemcpyAsync(hi_d, hi, 2, cudaMemcpyHostToDevice, stream)

  cuda_result = cudaMemcpyAsync(philo_d, philo, 2, cudaMemcpyHostToDevice, stream)
  cuda_result = cudaMemcpyAsync(phihi_d, phihi, 2, cudaMemcpyHostToDevice, stream)

  cuda_result = cudaMemcpyAsync(fxlo_d, fxlo, 2, cudaMemcpyHostToDevice, stream)
  cuda_result = cudaMemcpyAsync(fxhi_d, fxhi, 2, cudaMemcpyHostToDevice, stream)

  cuda_result = cudaMemcpyAsync(fylo_d, fylo, 2, cudaMemcpyHostToDevice, stream)
  cuda_result = cudaMemcpyAsync(fyhi_d, fyhi, 2, cudaMemcpyHostToDevice, stream)

  cuda_result = cudaMemcpyAsync(dx_d, dx, 2, cudaMemcpyHostToDevice, stream)

  call threads_and_blocks(lo, hi, numBlocks, numThreads)

  call cuda_compute_flux<<<numBlocks, numThreads, 0, stream>>>(lo_d, hi_d, phi, philo_d, phihi_d, &
                                                               fluxx, fxlo_d, fxhi_d, &
                                                               fluxy, fylo_d, fyhi_d, dx_d)

  cuda_result = cudaStreamSynchronize(stream)

#else

  call compute_flux_doit(lo, hi, phi, philo, phihi, fluxx, fxlo, fxhi, fluxy, fylo, fyhi, dx)

#endif

end subroutine compute_flux



#ifdef CUDA
attributes(global) &
subroutine cuda_compute_flux (lo, hi, phi, philo, phihi, &
     fluxx, fxlo, fxhi, fluxy, fylo, fyhi, dx)

  use amrex_fort_module, only : amrex_real
  use advance_module, only: compute_flux_doit

  implicit none

  integer lo(2), hi(2), philo(2), phihi(2), fxlo(2), fxhi(2), fylo(2), fyhi(2)
  real(amrex_real), intent(in)    :: phi  (philo(1):phihi(1),philo(2):phihi(2))
  real(amrex_real), intent(inout) :: fluxx( fxlo(1): fxhi(1), fxlo(2): fxhi(2))
  real(amrex_real), intent(inout) :: fluxy( fylo(1): fyhi(1), fylo(2): fyhi(2))
  real(amrex_real), intent(in)    :: dx(2)

  integer :: idx(2)

  ! Get our spatial index based on the CUDA thread index

  idx(1) = lo(1) + (threadIdx%x - 1) + blockDim%x * (blockIdx%x - 1)
  idx(2) = lo(2) + (threadIdx%y - 1) + blockDim%y * (blockIdx%y - 1)

  if (idx(1) .gt. hi(1) .or. idx(2) .gt. hi(2)) return

  call compute_flux_doit(idx, idx, phi, philo, phihi, fluxx, fxlo, fxhi, fluxy, fylo, fyhi, dx)

end subroutine cuda_compute_flux
#endif
