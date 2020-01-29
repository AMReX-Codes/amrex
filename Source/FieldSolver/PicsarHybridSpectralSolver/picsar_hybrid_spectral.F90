! Copyright 2019 Andrew Myers, Maxence Thevenet, Remi Lehe
! Weiqun Zhang
!
! This file is part of WarpX.
!
! License: BSD-3-Clause-LBNL


module warpx_fft_module
  use amrex_error_module, only : amrex_error, amrex_abort
  use amrex_fort_module, only : amrex_real
  use iso_c_binding
  implicit none

  include 'fftw3-mpi.f03'

  private
  public :: warpx_fft_mpi_init, warpx_fft_domain_decomp, warpx_fft_dataplan_init, warpx_fft_nullify, &
       warpx_fft_push_eb

contains

!> @brief
!! Set the communicator of the PICSAR module to the one that is passed in argument
  subroutine warpx_fft_mpi_init (fcomm) bind(c,name='warpx_fft_mpi_init')
    use shared_data, only : comm, rank, nproc
    integer, intent(in), value :: fcomm

    integer :: ierr, lnproc, lrank

    comm = fcomm

    call mpi_comm_size(comm, lnproc, ierr)
    nproc = lnproc

    call mpi_comm_rank(comm, lrank, ierr)
    rank = lrank

#ifdef _OPENMP
    ierr = fftw_init_threads()
    if (ierr.eq.0) call amrex_error("fftw_init_threads failed")
#endif
    call fftw_mpi_init()
#ifdef _OPENMP
    call dfftw_init_threads(ierr)
    if (ierr.eq.0) call amrex_error("dfftw_init_threads failed")
#endif
  end subroutine warpx_fft_mpi_init

!> @brief
!! Ask FFTW to do domain decomposition.
!
! This is always a 1d domain decomposition along z ; it is typically
! done on the *FFT sub-groups*, not the all domain
  subroutine warpx_fft_domain_decomp (warpx_local_nz, warpx_local_z0, global_lo, global_hi) &
       bind(c,name='warpx_fft_domain_decomp')
    use picsar_precision, only : idp
    use shared_data, only : comm, &
         nx_global, ny_global, nz_global, & ! size of global FFT
         nx, ny, nz ! size of local subdomains
    use mpi_fftw3, only : local_nz, local_z0, fftw_mpi_local_size_3d, alloc_local

    integer, intent(out) :: warpx_local_nz, warpx_local_z0
    integer, dimension(3), intent(in) :: global_lo, global_hi

    nx_global = INT(global_hi(1)-global_lo(1)+1,idp)
    ny_global = INT(global_hi(2)-global_lo(2)+1,idp)
    nz_global = INT(global_hi(3)-global_lo(3)+1,idp)

    alloc_local = fftw_mpi_local_size_3d( &
         INT(nz_global,C_INTPTR_T), &
         INT(ny_global,C_INTPTR_T), &
         INT(nx_global,C_INTPTR_T)/2+1, &
         comm, local_nz, local_z0)

    nx = nx_global
    ny = ny_global
    nz = local_nz

    warpx_local_nz = local_nz
    warpx_local_z0 = local_z0
  end subroutine warpx_fft_domain_decomp


!> @brief
!! Set all the flags and metadata of the PICSAR FFT module.
!! Allocate the auxiliary arrays of `fft_data`
!!
!! Note: fft_data is a stuct containing 22 pointers to arrays
!! 1-11: padded arrays in real space ; 12-22 arrays for the fields in Fourier space
  subroutine warpx_fft_dataplan_init (nox, noy, noz, fft_data, ndata, dx_wrpx, dt_wrpx, fftw_measure, do_nodal) &
       bind(c,name='warpx_fft_dataplan_init')
    USE picsar_precision, only: idp
    use shared_data, only : c_dim,  p3dfft_flag, fftw_plan_measure, &
         fftw_with_mpi, fftw_threads_ok, fftw_hybrid, fftw_mpi_transpose, &
         nx, ny, nz, & ! size of local subdomains
         nkx, nky, nkz, & ! size of local ffts
         iz_min_r, iz_max_r, iy_min_r, iy_max_r, ix_min_r, ix_max_r, & ! loop bounds
         dx, dy, dz
    use fields, only : nxguards, nyguards, nzguards,  & ! size of guard regions
         ex_r, ey_r, ez_r, bx_r, by_r, bz_r, &
         jx_r, jy_r, jz_r, rho_r, rhoold_r, &
         exf, eyf, ezf, bxf, byf, bzf, &
         jxf, jyf, jzf, rhof, rhooldf, &
         l_spectral, l_staggered, norderx, nordery, norderz
    use mpi_fftw3, only : alloc_local
    use omp_lib, only: omp_get_max_threads
    USE gpstd_solver, only: init_gpstd
    USE fourier_psaotd, only: init_plans_fourier_mpi
    use params, only : dt

    integer, intent(in) :: nox, noy, noz, ndata
    integer, intent(in) :: fftw_measure
    integer, intent(in) :: do_nodal
    type(c_ptr), intent(inout) :: fft_data(ndata)
    real(c_double), intent(in) :: dx_wrpx(3), dt_wrpx

    integer(idp) :: nopenmp
    integer :: nx_padded
    integer, dimension(3) :: shp
    integer(kind=c_size_t) :: sz

    ! No need to distinguish physical and guard cells for the global FFT;
    ! only nx+2*nxguards counts. Thus we declare 0 guard cells for simplicity
    nxguards = 0_idp
    nyguards = 0_idp
    nzguards = 0_idp

    ! For the calculation of the modified [k] vectors
    if (do_nodal == 0) then
       l_staggered = .TRUE.
    else
       l_staggered = .FALSE.
    endif
    norderx = int(nox, idp)
    nordery = int(noy, idp)
    norderz = int(noz, idp)
    ! Define parameters of FFT plans
    c_dim = INT(AMREX_SPACEDIM,idp)   ! Dimensionality of the simulation (2d/3d)
    fftw_with_mpi = .TRUE. ! Activate MPI FFTW
    fftw_hybrid = .FALSE.   ! FFT per MPI subgroup (instead of global)
    fftw_mpi_transpose = .FALSE. ! Do not transpose the data
    fftw_plan_measure = (fftw_measure .ne. 0)
    p3dfft_flag = .FALSE.
    l_spectral  = .TRUE.   ! Activate spectral Solver, using FFT
#ifdef _OPENMP
    fftw_threads_ok = .TRUE.
    nopenmp = OMP_GET_MAX_THREADS()
#else
    fftw_threads_ok = .FALSE.
    nopenmp = 1
#endif

    ! Allocate padded arrays for MPI FFTW
    nx_padded = 2*(nx/2 + 1)
    shp = [nx_padded, int(ny), int(nz)]
    sz = 2*alloc_local
    fft_data(1) = fftw_alloc_real(sz)
    call c_f_pointer(fft_data(1), ex_r, shp)
    fft_data(2) = fftw_alloc_real(sz)
    call c_f_pointer(fft_data(2), ey_r, shp)
    fft_data(3) = fftw_alloc_real(sz)
    call c_f_pointer(fft_data(3), ez_r, shp)
    fft_data(4) = fftw_alloc_real(sz)
    call c_f_pointer(fft_data(4), bx_r, shp)
    fft_data(5) = fftw_alloc_real(sz)
    call c_f_pointer(fft_data(5), by_r, shp)
    fft_data(6) = fftw_alloc_real(sz)
    call c_f_pointer(fft_data(6), bz_r, shp)
    fft_data(7) = fftw_alloc_real(sz)
    call c_f_pointer(fft_data(7), jx_r, shp)
    fft_data(8) = fftw_alloc_real(sz)
    call c_f_pointer(fft_data(8), jy_r, shp)
    fft_data(9) = fftw_alloc_real(sz)
    call c_f_pointer(fft_data(9), jz_r, shp)
    fft_data(10) = fftw_alloc_real(sz)
    call c_f_pointer(fft_data(10), rho_r, shp)
    fft_data(11) = fftw_alloc_real(sz)
    call c_f_pointer(fft_data(11), rhoold_r, shp)

    ! Set array bounds when copying ex to ex_r in PICSAR
    ix_min_r = 1; ix_max_r = nx
    iy_min_r = 1; iy_max_r = ny
    iz_min_r = 1; iz_max_r = nz
    ! Allocate Fourier space fields of the same size
    nkx = nx/2 + 1
    nky = ny
    nkz = nz
    shp = [int(nkx), int(nky), int(nkz)]
    sz = alloc_local
    fft_data(12) = fftw_alloc_complex(sz)
    call c_f_pointer(fft_data(12), exf, shp)
    fft_data(13) = fftw_alloc_complex(sz)
    call c_f_pointer(fft_data(13), eyf, shp)
    fft_data(14) = fftw_alloc_complex(sz)
    call c_f_pointer(fft_data(14), ezf, shp)
    fft_data(15) = fftw_alloc_complex(sz)
    call c_f_pointer(fft_data(15), bxf, shp)
    fft_data(16) = fftw_alloc_complex(sz)
    call c_f_pointer(fft_data(16), byf, shp)
    fft_data(17) = fftw_alloc_complex(sz)
    call c_f_pointer(fft_data(17), bzf, shp)
    fft_data(18) = fftw_alloc_complex(sz)
    call c_f_pointer(fft_data(18), jxf, shp)
    fft_data(19) = fftw_alloc_complex(sz)
    call c_f_pointer(fft_data(19), jyf, shp)
    fft_data(20) = fftw_alloc_complex(sz)
    call c_f_pointer(fft_data(20), jzf, shp)
    fft_data(21) = fftw_alloc_complex(sz)
    call c_f_pointer(fft_data(21), rhof, shp)
    fft_data(22) = fftw_alloc_complex(sz)
    call c_f_pointer(fft_data(22), rhooldf, shp)
!$acc enter data create (exf,eyf,ezf,bxf,byf,bzf,jxf,jyf,jzf,rhof,rhooldf)
    if (ndata < 22) then
       call amrex_abort("size of fft_data is too small")
    end if

    dx = dx_wrpx(1)
    dy = dx_wrpx(2)
    dz = dx_wrpx(3)
    dt = dt_wrpx

    ! Initialize the matrix blocks for the PSATD solver
    CALL init_gpstd()
    ! Initialize the plans for fftw with MPI
    CALL init_plans_fourier_mpi(nopenmp)

  end subroutine warpx_fft_dataplan_init


  subroutine warpx_fft_nullify () bind(c,name='warpx_fft_nullify')
    use fields, only: ex_r, ey_r, ez_r, bx_r, by_r, bz_r, &
         jx_r, jy_r, jz_r, rho_r, rhoold_r, &
         exf, eyf, ezf, bxf, byf, bzf, &
         jxf, jyf, jzf, rhof, rhooldf
    use mpi_fftw3, only : plan_r2c_mpi, plan_c2r_mpi
    nullify(ex_r)
    nullify(ey_r)
    nullify(ez_r)
    nullify(bx_r)
    nullify(by_r)
    nullify(bz_r)
    nullify(jx_r)
    nullify(jy_r)
    nullify(jz_r)
    nullify(rho_r)
    nullify(rhoold_r)
    nullify(exf)
    nullify(eyf)
    nullify(ezf)
    nullify(bxf)
    nullify(byf)
    nullify(bzf)
    nullify(jxf)
    nullify(jyf)
    nullify(jzf)
    nullify(rhof)
    nullify(rhooldf)
    call fftw_destroy_plan(plan_r2c_mpi)
    call fftw_destroy_plan(plan_c2r_mpi)
    call fftw_mpi_cleanup()
  end subroutine warpx_fft_nullify


  subroutine warpx_fft_push_eb ( &
       ex_wrpx, exlo, exhi, &
       ey_wrpx, eylo, eyhi, &
       ez_wrpx, ezlo, ezhi, &
       bx_wrpx, bxlo, bxhi, &
       by_wrpx, bylo, byhi, &
       bz_wrpx, bzlo, bzhi, &
       jx_wrpx, jxlo, jxhi, &
       jy_wrpx, jylo, jyhi, &
       jz_wrpx, jzlo, jzhi, &
       rhoold_wrpx, r1lo, r1hi, &
       rho_wrpx, r2lo, r2hi) &
       bind(c,name='warpx_fft_push_eb')

    use fields, only: ex, ey, ez, bx, by, bz, jx, jy, jz
    use shared_data, only: rhoold, rho
    use constants, only: num

    integer, dimension(3), intent(in) :: exlo, exhi, eylo, eyhi, ezlo, ezhi, bxlo, bxhi, &
         bylo, byhi, bzlo, bzhi, jxlo, jxhi, jylo, jyhi, jzlo, jzhi, r1lo, r1hi, r2lo, r2hi
    REAL(num), INTENT(INOUT), TARGET :: ex_wrpx(0:exhi(1)-exlo(1),0:exhi(2)-exlo(2),0:exhi(3)-exlo(3))
    REAL(num), INTENT(INOUT), TARGET :: ey_wrpx(0:eyhi(1)-eylo(1),0:eyhi(2)-eylo(2),0:eyhi(3)-eylo(3))
    REAL(num), INTENT(INOUT), TARGET :: ez_wrpx(0:ezhi(1)-ezlo(1),0:ezhi(2)-ezlo(2),0:ezhi(3)-ezlo(3))
    REAL(num), INTENT(INOUT), TARGET :: bx_wrpx(0:bxhi(1)-bxlo(1),0:bxhi(2)-bxlo(2),0:bxhi(3)-bxlo(3))
    REAL(num), INTENT(INOUT), TARGET :: by_wrpx(0:byhi(1)-bylo(1),0:byhi(2)-bylo(2),0:byhi(3)-bylo(3))
    REAL(num), INTENT(INOUT), TARGET :: bz_wrpx(0:bzhi(1)-bzlo(1),0:bzhi(2)-bzlo(2),0:bzhi(3)-bzlo(3))
    REAL(num), INTENT(INOUT), TARGET :: jx_wrpx(0:jxhi(1)-jxlo(1),0:jxhi(2)-jxlo(2),0:jxhi(3)-jxlo(3))
    REAL(num), INTENT(INOUT), TARGET :: jy_wrpx(0:jyhi(1)-jylo(1),0:jyhi(2)-jylo(2),0:jyhi(3)-jylo(3))
    REAL(num), INTENT(INOUT), TARGET :: jz_wrpx(0:jzhi(1)-jzlo(1),0:jzhi(2)-jzlo(2),0:jzhi(3)-jzlo(3))
    REAL(num), INTENT(INOUT), TARGET :: rhoold_wrpx(0:r1hi(1)-r1lo(1),0:r1hi(2)-r1lo(2),0:r1hi(3)-r1lo(3))
    REAL(num), INTENT(INOUT), TARGET :: rho_wrpx(0:r2hi(1)-r2lo(1),0:r2hi(2)-r2lo(2),0:r2hi(3)-r2lo(3))

    ! Point the fields in the PICSAR modules to the fields provided by WarpX
    ex => ex_wrpx
    ey => ey_wrpx
    ez => ez_wrpx
    bx => bx_wrpx
    by => by_wrpx
    bz => bz_wrpx
    jx => jx_wrpx
    jy => jy_wrpx
    jz => jz_wrpx
    rho => rho_wrpx
    rhoold => rhoold_wrpx

    ! Call the corresponding PICSAR function
    CALL push_psatd_ebfield()

    ex => null()
    ey => null()
    ez => null()
    bx => null()
    by => null()
    bz => null()
    jx => null()
    jy => null()
    jz => null()
    rho => null()
    rhoold => null()
  end subroutine warpx_fft_push_eb

end module warpx_fft_module
