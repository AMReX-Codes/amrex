
subroutine bar (comm) bind(c)
  use amrex_base_module
  implicit none
  integer, intent(in), value :: comm

  integer :: rank, ierr

  call amrex_init(comm)

  if (amrex_parallel_ioprocessor()) then
     print *, "bar: AMReX Fortran has been initialized."
  end if

  call fsolve()

  call amrex_finalize()

  ! After amrex_finalize(), amrex can no longer be used.
  call mpi_comm_rank(comm, rank, ierr)
  if (rank .eq. 0) then
     print *, "bar: AMReX Fortran has been finalized."
  end if
  
contains

  subroutine fsolve
    use amrex_linear_solver_module

    type(amrex_box) :: domain, bx
    type(amrex_geometry) :: geom
    type(amrex_boxarray) :: grids
    type(amrex_distromap) :: dm
    type(amrex_multifab) :: rhs, phi
    type(amrex_mfiter) :: mfi
    real(amrex_real), pointer :: prhs(:,:,:,:)
    integer :: lo(4), hi(4), i, j, k
    real(amrex_real) :: r, error
    type(amrex_poisson) :: poisson
    type(amrex_multigrid) :: multigrid

    call amrex_geometry_set_coord_sys(0)  ! Cartesian
    call amrex_geometry_set_prob_domain([0._amrex_real,0._amrex_real,0._amrex_real], &
         &                              [1._amrex_real,1._amrex_real,1._amrex_real])
    call amrex_geometry_set_periodic([.true., .true., .true.])

    domain = amrex_box([0,0,0], [63,63,63])  ! # of cells

    call amrex_geometry_build(geom, domain)
    call amrex_boxarray_build(grids, domain)
    call grids % maxSize(32)
    call amrex_distromap_build(dm, grids)
    call amrex_multifab_build(rhs, grids, dm, nc=1, ng=0)
    call amrex_multifab_build(phi, grids, dm, nc=1, ng=1)

    ! set right hand side to random numbers
    call amrex_mfiter_build(mfi, rhs)
    do while (mfi%next())
       bx = mfi%tilebox()
       prhs => rhs%dataptr(mfi)
       lo = lbound(prhs)
       hi = ubound(prhs)
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                call random_number(r)
                prhs(i,j,k,1) = r
             end do
          end do
       end do
    end do
    call amrex_mfiter_destroy(mfi)

    call phi % setVal(0.0_amrex_real) ! intial guess

    call amrex_poisson_build(poisson, [geom], [grids], [dm])
    call poisson % set_domain_bc([amrex_lo_periodic,amrex_lo_periodic,amrex_lo_periodic], &
         &                       [amrex_lo_periodic,amrex_lo_periodic,amrex_lo_periodic]);
    call poisson % set_level_bc(0, phi)

    call amrex_multigrid_build(multigrid, poisson)
    call multigrid % set_verbose(1)

    error = multigrid % solve([phi], [rhs], 1.e-10_amrex_real, 0.0_amrex_real)

    call amrex_poisson_destroy(poisson)
    call amrex_multigrid_destroy(multigrid)

    call amrex_multifab_destroy(rhs)
    call amrex_multifab_destroy(phi)
    call amrex_geometry_destroy(geom)
    call amrex_boxarray_destroy(grids)
    call amrex_distromap_destroy(dm)
  end subroutine fsolve

end subroutine bar
