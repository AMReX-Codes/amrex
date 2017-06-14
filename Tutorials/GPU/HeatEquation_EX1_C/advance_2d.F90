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

    integer,  intent(in   ) :: lo(2), hi(2), philo(2), phihi(2), flo(2), fhi(2)
    real(rt), intent(in   ) :: phi (philo(1):phihi(1),philo(2):phihi(2))
    real(rt), intent(inout) :: flux(  flo(1):  fhi(1),  flo(2):  fhi(2))
    real(rt), intent(in   ) :: dx(2)
    integer,  intent(in   ), value :: idx, idir

#ifdef CUDA
    attributes(device) :: phi, flux, dx, philo, phihi, flo, fhi
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
                         fluxx, fxlo, fxhi, fluxy, fylo, fyhi, dx, dt, idx) bind(C, name="update_phi")

    use amrex_fort_module, only: rt => amrex_real
#ifdef CUDA
    use cuda_module, only: cuda_streams, stream_from_index, threads_and_blocks
    use cudafor, only: cuda_stream_kind, dim3
#endif

    implicit none

    integer,  intent(in   ) :: lo(2), hi(2), polo(2), pohi(2), pnlo(2), pnhi(2), fxlo(2), fxhi(2), fylo(2), fyhi(2)
    real(rt), intent(in   ) :: phiold(polo(1):pohi(1),polo(2):pohi(2))
    real(rt), intent(inout) :: phinew(pnlo(1):pnhi(1),pnlo(2):pnhi(2))
    real(rt), intent(in   ) :: fluxx (fxlo(1):fxhi(1),fxlo(2):fxhi(2))
    real(rt), intent(in   ) :: fluxy (fylo(1):fyhi(1),fylo(2):fyhi(2))
    real(rt), intent(in   ) :: dx(2)
    integer,  intent(in   ), value :: idx
    real(rt), intent(in   ), value :: dt

#ifdef CUDA
    attributes(device) :: phiold, phinew, fluxx, fluxy, dx, polo, pohi, pnlo, pnhi, fxlo, fxhi, fylo, fyhi
    attributes(managed) :: lo, hi

    integer :: cuda_result
    integer(kind=cuda_stream_kind) :: stream
    type(dim3) :: numThreads, numBlocks

    stream = cuda_streams(stream_from_index(idx)+1)

    call threads_and_blocks(lo, hi, numBlocks, numThreads)
#endif

    call update_phi_doit &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, stream>>> &
#endif
         (lo, hi, phiold, polo, pohi, &
          phinew, pnlo, pnhi, fluxx, fxlo, fxhi, &
          fluxy, fylo, fyhi, dx, dt)

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
    integer,  intent(in   ), value :: idir

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

    integer,  intent(in   ) ::  lo(2), hi(2), polo(2), pohi(2), pnlo(2), pnhi(2), fxlo(2), fxhi(2), fylo(2), fyhi(2)
    real(rt), intent(in   ) :: phiold(polo(1):pohi(1),polo(2):pohi(2))
    real(rt), intent(inout) :: phinew(pnlo(1):pnhi(1),pnlo(2):pnhi(2))
    real(rt), intent(in   ) :: fluxx (fxlo(1):fxhi(1),fxlo(2):fxhi(2))
    real(rt), intent(in   ) :: fluxy (fylo(1):fyhi(1),fylo(2):fyhi(2))
    real(rt), intent(in   ) :: dx(2)
    real(rt), intent(in   ), value :: dt

    ! local variables
    integer i,j
    integer :: blo(3), bhi(3)
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
