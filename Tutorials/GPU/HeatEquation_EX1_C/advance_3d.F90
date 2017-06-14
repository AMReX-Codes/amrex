module advance_module

  implicit none

contains

  subroutine compute_flux (lo, hi, phi, philo, phihi, &
                           flux, flo, fhi, dx, idir, idx) bind(C, name="compute_flux")

    use amrex_fort_module, only: rt => amrex_real
#ifdef CUDA
    use cuda_module, only: cuda_streams, stream_from_index, threads_and_blocks
    use cudafor, only: cudaMemcpyAsync, cudaMemcpyHostToDevice, cuda_stream_kind, dim3
#endif

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3), philo(3), phihi(3), flo(3), fhi(3)
    real(rt), intent(in   ) :: phi (philo(1):phihi(1),philo(2):phihi(2),philo(3):phihi(3))
    real(rt), intent(inout) :: flux(  flo(1):  fhi(1),  flo(2):  fhi(2),  flo(3):  fhi(3))
    real(rt), intent(in   ) :: dx(3)
    integer,  intent(in   ), value :: idx, idir

#ifdef CUDA
    attributes(device) :: phi, flux, philo, phihi, flo, fhi, dx
    attributes(managed) :: lo, hi

    integer :: cuda_result
    integer(kind=cuda_stream_kind) :: stream
    type(dim3) :: numThreads, numBlocks

    stream = cuda_streams(stream_from_index(idx)+1)

    call threads_and_blocks(lo, hi, numBlocks, numThreads)
#endif

    call compute_flux_doit &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, stream>>> &
#endif
         (lo, hi, phi, philo, phihi, flux, flo, fhi, dx, idir)

  end subroutine compute_flux



  subroutine update_phi (lo, hi, phiold, polo, pohi, phinew, pnlo, pnhi, &
                         fluxx, fxlo, fxhi, &
                         fluxy, fylo, fyhi, &
                         fluxz, fzlo, fzhi, &
                         dx, dt, idx) bind(C, name="update_phi")

    use amrex_fort_module, only: rt => amrex_real
#ifdef CUDA
    use cuda_module, only: cuda_streams, stream_from_index, threads_and_blocks
    use cudafor, only: cudaMemcpyAsync, cudaMemcpyHostToDevice, cuda_stream_kind, dim3
#endif

    implicit none

    integer lo(3), hi(3), polo(3), pohi(3), pnlo(3), pnhi(3), &
        fxlo(3), fxhi(3), &
        fylo(3), fyhi(3), &
        fzlo(3), fzhi(3), &
        idx
    real(rt), intent(in)    :: phiold(polo(1):pohi(1),polo(2):pohi(2),polo(3):pohi(3))
    real(rt), intent(inout) :: phinew(pnlo(1):pnhi(1),pnlo(2):pnhi(2),pnlo(3):pnhi(3))
    real(rt), intent(in   ) :: fluxx (fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3))
    real(rt), intent(in   ) :: fluxy (fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3))
    real(rt), intent(in   ) :: fluxz (fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3))
    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(in   ) :: dt

#ifdef CUDA
    attributes(device) :: phiold, phinew, fluxx, fluxy, fluxz, polo, pohi, pnlo, pnhi, &
                          fxlo, fxhi, fylo, fyhi, fzlo, fzhi, dx
    attributes(managed) :: lo, hi

    integer :: cuda_result
    integer(kind=cuda_stream_kind) :: stream
    type(dim3) :: numThreads, numBlocks

    real(rt), device :: dx_d(3), dt_d

    stream = cuda_streams(stream_from_index(idx)+1)

    call threads_and_blocks(lo, hi, numBlocks, numThreads)

    call update_phi_doit<<<numBlocks, numThreads, 0, stream>>>(lo, hi, phiold, polo, pohi, &
                                                               phinew, pnlo, pnhi, &
                                                               fluxx, fxlo, fxhi, &
                                                               fluxy, fylo, fyhi, &
                                                               fluxz, fzlo, fzhi, &
                                                               dx, dt)

#else

    call update_phi_doit(lo, hi, phiold, polo, pohi, phinew, pnlo, pnhi, &
                         fluxx, fxlo, fxhi, &
                         fluxy, fylo, fyhi, &
                         fluxz, fzlo, fzhi, &
                         dx, dt)

#endif

  end subroutine update_phi

#ifdef CUDA
  attributes(global) &
#endif
  subroutine compute_flux_doit (lo, hi, phi, p_lo, p_hi, flx, f_lo, f_hi, dx, idir)

    use amrex_fort_module, only: rt => amrex_real, get_loop_bounds

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3), p_lo(3), p_hi(3), f_lo(3), f_hi(3)
    real(rt), intent(in   ) :: phi(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3))
    real(rt), intent(inout) :: flx(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3))
    real(rt), intent(in   ) :: dx(3)
    integer,  intent(in   ), value :: idir

    ! local variables
    integer :: i, j, k
    integer :: blo(3), bhi(3)

    call get_loop_bounds(blo, bhi, lo, hi)

    do         k = blo(3), bhi(3)
        do     j = blo(2), bhi(2)
            do i = blo(1), bhi(1)
                if (idir .eq. 1) then
                    ! x-flux
                    flx(i,j,k) = ( phi(i,j,k) - phi(i-1,j,k) ) / dx(1)
                else if (idir .eq. 2) then
                    ! y-flux
                    flx(i,j,k) = ( phi(i,j,k) - phi(i,j-1,k) ) / dx(2)
                else
                    ! z-flux
                    flx(i,j,k) = ( phi(i,j,k) - phi(i,j,k-1) ) / dx(3)
                endif
            end do
        end do
    end do

  end subroutine compute_flux_doit



#ifdef CUDA
  attributes(global) &
#endif
  subroutine update_phi_doit (lo, hi, phiold, polo, pohi, phinew, pnlo, pnhi, &
                              fluxx, fxlo, fxhi, &
                              fluxy, fylo, fyhi, &
                              fluxz, fzlo, fzhi, &
                              dx, dt)

    use amrex_fort_module, only: rt => amrex_real, get_loop_bounds

    implicit none

    integer lo(3), hi(3), polo(3), pohi(3), pnlo(3), pnhi(3), &
        fxlo(3), fxhi(3), &
        fylo(3), fyhi(3), &
        fzlo(3), fzhi(3)
    real(rt), intent(in)    :: phiold(polo(1):pohi(1),polo(2):pohi(2),polo(3):pohi(3))
    real(rt), intent(inout) :: phinew(pnlo(1):pnhi(1),pnlo(2):pnhi(2),pnlo(3):pnhi(3))
    real(rt), intent(in   ) :: fluxx (fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3))
    real(rt), intent(in   ) :: fluxy (fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3))
    real(rt), intent(in   ) :: fluxz (fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3))
    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(in   ), value :: dt

    ! local variables
    integer i,j,k
    integer :: blo(3), bhi(3)
    real(rt) :: dtdx(3)

    call get_loop_bounds(blo, bhi, lo, hi)

    dtdx = dt/dx

    do        k = blo(3), bhi(3)
        do    j = blo(2), bhi(2)
           do i = blo(1), bhi(1)

              phinew(i,j,k) = phiold(i,j,k) &
                   + dtdx(1) * (fluxx(i+1,j  ,k  ) - fluxx(i,j,k)) &
                   + dtdx(2) * (fluxy(i  ,j+1,k  ) - fluxy(i,j,k)) &
                   + dtdx(3) * (fluxz(i  ,j  ,k+1) - fluxz(i,j,k))

           end do
        end do
    end do

  end subroutine update_phi_doit

end module advance_module
