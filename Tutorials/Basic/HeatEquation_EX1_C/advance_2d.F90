module advance_module

  implicit none

contains

  subroutine compute_flux (lo, hi, phi, philo, phihi, &
       fluxx, fxlo, fxhi, fluxy, fylo, fyhi, dx, idx) bind(C, name="compute_flux")

    use amrex_fort_module, only: rt => amrex_real
#ifdef CUDA
    use cuda_module, only: cuda_streams, stream_from_index, threads_and_blocks
    use cudafor, only: cudaMemcpyAsync, cudaMemcpyHostToDevice, cudaStreamSynchronize, &
                       cuda_stream_kind, dim3
#endif

    implicit none

    integer lo(2), hi(2), philo(2), phihi(2), fxlo(2), fxhi(2), fylo(2), fyhi(2), idx
    real(rt), intent(in   ) :: phi  (philo(1):phihi(1),philo(2):phihi(2))
    real(rt), intent(inout) :: fluxx( fxlo(1): fxhi(1), fxlo(2): fxhi(2))
    real(rt), intent(inout) :: fluxy( fylo(1): fyhi(1), fylo(2): fyhi(2))
    real(rt), intent(in   ) :: dx(2)

    integer :: blo(2), bhi(2), idir

#ifdef CUDA
    attributes(device) :: phi, fluxx, fluxy

    integer :: cuda_result
    integer(kind=cuda_stream_kind) :: stream
    type(dim3) :: numThreads, numBlocks

    integer, device :: blo_d(2), bhi_d(2)
    integer, device :: philo_d(2), phihi_d(2)
    integer, device :: fxlo_d(2), fxhi_d(2)
    integer, device :: fylo_d(2), fyhi_d(2)
    real(rt), device :: dx_d(2)
    integer, device :: idir_d

    stream = cuda_streams(stream_from_index(idx)+1)

    cuda_result = cudaMemcpyAsync(philo_d, philo, 2, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(phihi_d, phihi, 2, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(fxlo_d, fxlo, 2, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(fxhi_d, fxhi, 2, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(fylo_d, fylo, 2, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(fyhi_d, fyhi, 2, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(dx_d, dx, 2, cudaMemcpyHostToDevice, stream)

    blo = [lo(1),   lo(2)]
    bhi = [hi(1)+1, hi(2)]

    cuda_result = cudaMemcpyAsync(blo_d, blo, 2, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(bhi_d, bhi, 2, cudaMemcpyHostToDevice, stream)

    idir = 1
    cuda_result = cudaMemcpyAsync(idir_d, idir, 1, cudaMemcpyHostToDevice, stream)

    call threads_and_blocks(blo, bhi, numBlocks, numThreads)

    call compute_flux_doit<<<numBlocks, numThreads, 0, stream>>>(blo_d, bhi_d, phi, philo_d, phihi_d, &
                                                                 fluxx, fxlo_d, fxhi_d, dx_d, idir_d)

    bhi = [hi(1), hi(2)+1]

    cuda_result = cudaMemcpyAsync(bhi_d, bhi, 2, cudaMemcpyHostToDevice, stream)

    idir = 2
    cuda_result = cudaMemcpyAsync(idir_d, idir, 1, cudaMemcpyHostToDevice, stream)

    call threads_and_blocks(blo, bhi, numBlocks, numThreads)

    call compute_flux_doit<<<numBlocks, numThreads, 0, stream>>>(blo_d, bhi_d, phi, philo_d, phihi_d, &
                                                                 fluxy, fylo_d, fyhi_d, dx_d, idir_d)

#else

    blo = [lo(1),   lo(2)]
    bhi = [hi(1)+1, hi(2)]

    idir = 1

    call compute_flux_doit(blo, bhi, phi, philo, phihi, fluxx, fxlo, fxhi, dx, idir)

    bhi = [hi(1), hi(2)+1]

    idir = 2

    call compute_flux_doit(blo, bhi, phi, philo, phihi, fluxy, fylo, fyhi, dx, idir)

#endif

  end subroutine compute_flux



  subroutine update_phi (lo, hi, phiold, polo, pohi, phinew, pnlo, pnhi, &
       fluxx, fxlo, fxhi, fluxy, fylo, fyhi, dx, dt, idx) bind(C, name="update_phi")

    use amrex_fort_module, only: rt => amrex_real
#ifdef CUDA
    use cuda_module, only: cuda_streams, stream_from_index, threads_and_blocks
    use cudafor, only: cudaMemcpyAsync, cudaMemcpyHostToDevice, cudaStreamSynchronize, &
                       cuda_stream_kind, dim3
#endif

    implicit none

    integer lo(2), hi(2), polo(2), pohi(2), pnlo(2), pnhi(2), fxlo(2), fxhi(2), fylo(2), fyhi(2), idx
    real(rt), intent(in)    :: phiold(polo(1):pohi(1),polo(2):pohi(2))
    real(rt), intent(inout) :: phinew(pnlo(1):pnhi(1),pnlo(2):pnhi(2))
    real(rt), intent(in   ) :: fluxx (fxlo(1):fxhi(1),fxlo(2):fxhi(2))
    real(rt), intent(in   ) :: fluxy (fylo(1):fyhi(1),fylo(2):fyhi(2))
    real(rt), intent(in   ) :: dx(2)
    real(rt), intent(in   ) :: dt

#ifdef CUDA
    attributes(device) :: phiold, phinew, fluxx, fluxy

    integer :: cuda_result
    integer(kind=cuda_stream_kind) :: stream
    type(dim3) :: numThreads, numBlocks

    integer, device :: lo_d(2), hi_d(2)
    integer, device :: polo_d(2), pohi_d(2)
    integer, device :: pnlo_d(2), pnhi_d(2)
    integer, device :: fxlo_d(2), fxhi_d(2)
    integer, device :: fylo_d(2), fyhi_d(2)
    real(rt), device :: dx_d(2), dt_d

    stream = cuda_streams(stream_from_index(idx)+1)

    cuda_result = cudaMemcpyAsync(lo_d, lo, 2, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(hi_d, hi, 2, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(polo_d, polo, 2, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(pohi_d, pohi, 2, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(pnlo_d, pnlo, 2, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(pnhi_d, pnhi, 2, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(fxlo_d, fxlo, 2, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(fxhi_d, fxhi, 2, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(fylo_d, fylo, 2, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(fyhi_d, fyhi, 2, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(dx_d, dx, 2, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(dt_d, dt, 1, cudaMemcpyHostToDevice, stream)

    call threads_and_blocks(lo, hi, numBlocks, numThreads)

    call update_phi_doit<<<numBlocks, numThreads, 0, stream>>>(lo_d, hi_d, phiold, polo_d, pohi_d, &
                                                               phinew, pnlo_d, pnhi_d, fluxx, fxlo_d, fxhi_d, &
                                                               fluxy, fylo_d, fyhi_d, dx_d, dt_d)

#else

    call update_phi_doit(lo, hi, phiold, polo, pohi, phinew, pnlo, pnhi, &
                         fluxx, fxlo, fxhi, fluxy, fylo, fyhi, dx, dt)

#endif

  end subroutine update_phi

#ifdef CUDA
  attributes(global) &
#endif
  subroutine compute_flux_doit (lo, hi, phi, p_lo, p_hi, flx, f_lo, f_hi, dx, idir)

    use amrex_fort_module, only: rt => amrex_real, get_loop_bounds

    implicit none

    integer,  intent(in   ) :: lo(2), hi(2), p_lo(2), p_hi(2), f_lo(2), f_hi(2)
    real(rt), intent(in   ) :: phi(p_lo(1):p_hi(1),p_lo(2):p_hi(2))
    real(rt), intent(inout) :: flx(f_lo(1):f_hi(1),f_lo(2):f_hi(2))
    real(rt), intent(in   ) :: dx(2)
    integer,  intent(in   ) :: idir

    ! local variables
    integer :: i, j
    integer :: blo(3), bhi(3)

    call get_loop_bounds(blo, bhi, [lo(1), lo(2), 0], [hi(1), hi(2), 0])

    do    j = blo(2), bhi(2)
       do i = blo(1), bhi(1)
          if (idir .eq. 1) then
             ! x-flux
             flx(i,j) = ( phi(i,j) - phi(i-1,j) ) / dx(1)
          else
             ! y-flux
             flx(i,j) = ( phi(i,j) - phi(i,j-1) ) / dx(2)
          endif
       end do
    end do

  end subroutine compute_flux_doit



#ifdef CUDA
  attributes(global) &
#endif
  subroutine update_phi_doit (lo, hi, phiold, polo, pohi, phinew, pnlo, pnhi, &
                              fluxx, fxlo, fxhi, fluxy, fylo, fyhi, dx, dt)

    use amrex_fort_module, only: rt => amrex_real, get_loop_bounds

    implicit none

    integer lo(2), hi(2), polo(2), pohi(2), pnlo(2), pnhi(2), fxlo(2), fxhi(2), fylo(2), fyhi(2)
    real(rt), intent(in)    :: phiold(polo(1):pohi(1),polo(2):pohi(2))
    real(rt), intent(inout) :: phinew(pnlo(1):pnhi(1),pnlo(2):pnhi(2))
    real(rt), intent(in   ) :: fluxx (fxlo(1):fxhi(1),fxlo(2):fxhi(2))
    real(rt), intent(in   ) :: fluxy (fylo(1):fyhi(1),fylo(2):fyhi(2))
    real(rt), intent(in   ) :: dx(2)
    real(rt), intent(in   ) :: dt

    ! local variables
    integer i,j
    integer :: blo(2), bhi(2)
    real(rt) :: dtdx(2)

    call get_loop_bounds(blo, bhi, [lo(1), lo(2), 0], [hi(1), hi(2), 0])

    dtdx = dt/dx

    do    j = blo(2), bhi(2)
       do i = blo(1), bhi(1)

          phinew(i,j) = phiold(i,j) &
               + dtdx(1) * (fluxx(i+1,j  ) - fluxx(i,j)) &
               + dtdx(2) * (fluxy(i  ,j+1) - fluxy(i,j))

       end do
    end do

  end subroutine update_phi_doit

end module advance_module
