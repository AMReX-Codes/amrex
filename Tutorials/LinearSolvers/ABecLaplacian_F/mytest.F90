
module mytest_module

  use amrex_base_module
  use amrex_linear_solver_module

  implicit none

  ! parameters

  integer, save :: max_level = 1
  integer, save :: ref_ratio = 2
  integer, save :: n_cell = 128
  integer, save :: max_grid_size = 64

  logical, save :: composite_solve = .true.

  ! prob_type 1 here is Poisson with homogeneous Dirichlet boundary.
  ! prob_type 2 here is ABecLaplacian with homogeneous Neumann boundary.
  integer, save :: prob_type = 1

  integer, save :: verbose = 2
  integer, save :: cg_verbose = 0
  integer, save :: max_iter = 100
  integer, save :: max_fmg_iter = 0
  integer, save :: bottom_solver = amrex_bottom_default
  integer, save :: linop_maxorder = 2
  logical, save :: agglomeration = .true.
  logical, save :: consolidation = .true.
  integer, save :: max_coarsening_level = 30

  ! data
  type(amrex_geometry), allocatable, save :: geom(:)
  type(amrex_boxarray), allocatable, save :: ba(:)
  type(amrex_distromap), allocatable, save :: dm(:)

  type(amrex_multifab), allocatable, save :: solution(:)
  type(amrex_multifab), allocatable, save :: rhs(:)
  type(amrex_multifab), allocatable, save :: exact_solution(:)
  type(amrex_multifab), allocatable, save :: acoef(:)
  type(amrex_multifab), allocatable, save :: bcoef(:)
  real(amrex_real), save :: ascalar, bscalar
  

  private
  public :: init, finalize, solve, write_plotfile

contains

  subroutine init ()
    call init_parameters()
    call init_data()
  end subroutine init


  subroutine init_parameters()
    type(amrex_parmparse) pp

    call amrex_parmparse_build(pp)

    call pp % query("max_level", max_level)
    call pp % query("ref_ratio", ref_ratio)
    call pp % query("n_cell", n_cell);
    call pp % query("max_grid_size", max_grid_size)

    call pp % query("composite_solve", composite_solve)

    call pp % query("prob_type", prob_type)

    call pp % query("verbose", verbose)
    call pp % query("cg_verbose", cg_verbose)
    call pp % query("max_iter", max_iter)
    call pp % query("max_fmg_iter", max_fmg_iter)
    call pp % query("bottom_solver", bottom_solver)
    call pp % query("linop_maxorder", linop_maxorder)
    call pp % query("agglomeration", agglomeration)
    call pp % query("consolidation", consolidation)
    call pp % query("max_coarsening_level", max_coarsening_level)

    call amrex_parmparse_destroy(pp)
  end subroutine init_parameters


  subroutine init_data ()
    use init_prob_module
    allocate(geom(0:max_level))
    allocate(ba(0:max_level))
    allocate(dm(0:max_level))
    allocate(solution(0:max_level))
    allocate(rhs(0:max_level))
    allocate(exact_solution(0:max_level))
    if (prob_type .eq. 2) then
       allocate(acoef(0:max_level))
       allocate(bcoef(0:max_level))
    end if

    call init_geom()
    call init_grids()
    call init_mf()
    if (prob_type .eq. 1) then
       call init_prob_poisson(geom,solution, rhs, exact_solution)
    else
       call init_prob_abeclaplacian(geom,solution, rhs, exact_solution, acoef, bcoef, ascalar, bscalar)
    end if
  end subroutine init_data


  subroutine init_geom ()
    integer :: ilev
    type(amrex_box) :: domain
    call amrex_geometry_set_coord_sys(0)  ! Cartesian
    call amrex_geometry_set_prob_domain([0._amrex_real,0._amrex_real,0._amrex_real], &
         &                              [1._amrex_real,1._amrex_real,1._amrex_real])
    call amrex_geometry_set_periodic([.false., .false., .false.])
    domain = amrex_box([0,0,0], [n_cell-1,n_cell-1,n_cell-1])
    do ilev = 0, max_level
       call amrex_geometry_build(geom(ilev), domain)
       call domain % refine(ref_ratio)
    end do
  end subroutine init_geom


  subroutine init_grids ()
    integer :: ilev
    type(amrex_box) :: dom
    dom = geom(0) % domain
    do ilev = 0, max_level
       call amrex_boxarray_build(ba(ilev), dom)
       call ba(ilev) % maxSize(max_grid_size)
       call dom % grow(-n_cell/4)     ! fine level cover the middle of the coarse domain
       call dom % refine(ref_ratio)
    end do
  end subroutine init_grids


  subroutine init_mf ()
    integer :: ilev
    do ilev = 0, max_level
       call amrex_distromap_build(dm(ilev),ba(ilev))
       ! one ghost cell to store boundary conditions
       call amrex_multifab_build(solution(ilev), ba(ilev), dm(ilev), nc=1, ng=1)
       call amrex_multifab_build(rhs(ilev), ba(ilev), dm(ilev), nc=1, ng=0)
       call amrex_multifab_build(exact_solution(ilev), ba(ilev), dm(ilev), nc=1, ng=0)
       if (allocated(acoef)) then
          call amrex_multifab_build(acoef(ilev), ba(ilev), dm(ilev), nc=1, ng=0)
          ! 1 ghost cell for averaging from cell centers to faces
          call amrex_multifab_build(bcoef(ilev), ba(ilev), dm(ilev), nc=1, ng=1)
       end if
    end do
  end subroutine init_mf


  subroutine finalize ()
    integer :: ilev, idim
    do ilev = 0, max_level
       call amrex_geometry_destroy(geom(ilev))
       call amrex_boxarray_destroy(ba(ilev))
       call amrex_distromap_destroy(dm(ilev))
       call amrex_multifab_destroy(solution(ilev))
       call amrex_multifab_destroy(rhs(ilev))
       call amrex_multifab_destroy(exact_solution(ilev))
       if (allocated(acoef)) then
          call amrex_multifab_destroy(acoef(ilev))
          call amrex_multifab_destroy(bcoef(ilev))
       end if
    end do
  end subroutine finalize


  subroutine solve ()
    if (prob_type .eq. 1) then
       call solve_poisson()
    else
       call solve_abeclaplacian()
    end if
  end subroutine solve


  subroutine solve_poisson ()
    type(amrex_poisson) :: poisson
    type(amrex_multigrid) :: multigrid
    integer :: ilev
    real(amrex_real) :: err

    if (composite_solve) then

       call amrex_poisson_build(poisson, geom, ba, dm, &
            metric_term=.false., agglomeration=agglomeration, consolidation=consolidation, &
            max_coarsening_level=max_coarsening_level)
       
       call poisson % set_maxorder(linop_maxorder)

       ! This is a 3d problem with Dirichlet BC
       call poisson % set_domain_bc([amrex_lo_dirichlet, amrex_lo_dirichlet, amrex_lo_dirichlet], &
            &                       [amrex_lo_dirichlet, amrex_lo_dirichlet, amrex_lo_dirichlet])

       do ilev = 0, max_level
          ! solution multifab's ghost cells at physical boundaries have been set to bc values.
          call poisson % set_level_bc(ilev, solution(ilev))
       end do

       call amrex_multigrid_build(multigrid, poisson)
       call multigrid % set_verbose(verbose)
       call multigrid % set_cg_verbose(cg_verbose)
       call multigrid % set_max_iter(max_iter)
       call multigrid % set_max_fmg_iter(max_fmg_iter)
       call multigrid % set_bottom_solver(bottom_solver)

       err = multigrid % solve(solution, rhs, 1.e-10_amrex_real, 0.0_amrex_real)

       call amrex_poisson_destroy(poisson)
       call amrex_multigrid_destroy(multigrid)

    else
       do ilev = 0, max_level

          call amrex_poisson_build(poisson, [geom(ilev)], [ba(ilev)], [dm(ilev)], &
               metric_term=.false., agglomeration=agglomeration, consolidation=consolidation, &
               max_coarsening_level=max_coarsening_level)
       
          call poisson % set_maxorder(linop_maxorder)

          ! The order of the following set bc calls matters.

          ! This is a 3d problem with Dirichlet BC
          call poisson % set_domain_bc([amrex_lo_dirichlet, amrex_lo_dirichlet, amrex_lo_dirichlet], &
               &                       [amrex_lo_dirichlet, amrex_lo_dirichlet, amrex_lo_dirichlet])
               
          if (ilev > 0) then
             ! use coarse level data to set up bc at corase/fine boundary
             call poisson % set_coarse_fine_bc(solution(ilev-1), ref_ratio)
          end if

          ! Note that to the linear solver, the level is ZERO.  In
          ! this test problem, when lev > 0, solution(lev) is going to
          ! be ignored because fine level grids are completed
          ! surrounded by coarse level.  If fine level grids do touch
          ! phyical domain, the multifab must have bc values at
          ! physical boundaries stored in ghost cells.
          call poisson % set_level_bc(0, solution(ilev))

          call amrex_multigrid_build(multigrid, poisson);
          call multigrid % set_verbose(verbose)
          call multigrid % set_cg_verbose(cg_verbose)
          call multigrid % set_max_iter(max_iter)
          call multigrid % set_max_fmg_iter(max_fmg_iter)
          call multigrid % set_bottom_solver(bottom_solver)

          err = multigrid % solve([solution(ilev)], [rhs(ilev)], 1.e-10_amrex_real, 0.0_amrex_real)

          call amrex_poisson_destroy(poisson)
          call amrex_multigrid_destroy(multigrid)

       end do
    end if

  end subroutine solve_poisson


  subroutine solve_abeclaplacian ()
    type(amrex_abeclaplacian) :: abeclap
    type(amrex_multigrid) :: multigrid
    integer :: ilev, idim
    integer(amrex_long) :: npts
    real(amrex_real) :: err, avg1, avg2, offset
    type(amrex_multifab), allocatable :: beta(:,:)
    logical :: nodal(3)
    type(amrex_multifab) :: nullmf

    ! For ABecLaplacian, the b coefficents are on faces
    allocate(beta(amrex_spacedim,0:max_level))
    do ilev = 0, max_level
       do idim = 1, amrex_spacedim
          nodal = .false.
          nodal(idim) = .true.
          call amrex_multifab_build(beta(idim,ilev), ba(ilev), dm(ilev), 1, 0, nodal)
       end do
       call amrex_average_cellcenter_to_face(beta(:,ilev), bcoef(ilev), geom(ilev))
    end do

    if (composite_solve) then

       call amrex_abeclaplacian_build(abeclap, geom, ba, dm, &
            metric_term=.false., agglomeration=agglomeration, consolidation=consolidation, &
            max_coarsening_level=max_coarsening_level)

       call abeclap % set_maxorder(linop_maxorder)

       ! This is set up to have homogeneous Neumann BC
       call abeclap % set_domain_bc([amrex_lo_neumann, amrex_lo_neumann, amrex_lo_neumann], &
            &                       [amrex_lo_neumann, amrex_lo_neumann, amrex_lo_neumann])

       do ilev = 0, max_level
          ! for problem with pure homogeneous Neumann BC, we could pass an empty multifab
          call abeclap % set_level_bc(ilev, nullmf)
       end do

       call abeclap % set_scalars(ascalar, bscalar)
       do ilev = 0, max_level
          call abeclap % set_acoeffs(ilev, acoef(ilev))
          call abeclap % set_bcoeffs(ilev, beta(:,ilev))
       end do

       call amrex_multigrid_build(multigrid, abeclap)
       call multigrid % set_verbose(verbose)
       call multigrid % set_cg_verbose(cg_verbose)
       call multigrid % set_max_iter(max_iter)
       call multigrid % set_max_fmg_iter(max_fmg_iter)
       call multigrid % set_bottom_solver(bottom_solver)

       err = multigrid % solve(solution, rhs, 1.e-10_amrex_real, 0.0_amrex_real)

       call amrex_abeclaplacian_destroy(abeclap)
       call amrex_multigrid_destroy(multigrid)

    else
       do ilev = 0, max_level
          
          call amrex_abeclaplacian_build(abeclap, [geom(ilev)], [ba(ilev)], [dm(ilev)], &
               metric_term=.false., agglomeration=agglomeration, consolidation=consolidation, &
               max_coarsening_level=max_coarsening_level)

       call abeclap % set_maxorder(linop_maxorder)

       ! This is set up to have homogeneous Neumann BC
       call abeclap % set_domain_bc([amrex_lo_neumann, amrex_lo_neumann, amrex_lo_neumann], &
            &                       [amrex_lo_neumann, amrex_lo_neumann, amrex_lo_neumann])

       if (ilev > 0) then
          ! use coarse level data to set up bc at corase/fine boundary
          call abeclap % set_coarse_fine_bc(solution(ilev-1), ref_ratio)
       end if

       ! for problem with pure homogeneous Neumann BC, we could pass an empty multifab
       call abeclap % set_level_bc(0, nullmf)

       call abeclap % set_scalars(ascalar, bscalar)
       call abeclap % set_acoeffs(0, acoef(ilev))
       call abeclap % set_bcoeffs(0, beta(:,ilev))

       call amrex_multigrid_build(multigrid, abeclap)
       call multigrid % set_verbose(verbose)
       call multigrid % set_cg_verbose(cg_verbose)
       call multigrid % set_max_iter(max_iter)
       call multigrid % set_max_fmg_iter(max_fmg_iter)
       call multigrid % set_bottom_solver(bottom_solver)

       err = multigrid % solve([solution(ilev)], [rhs(ilev)], 1.e-10_amrex_real, 0.0_amrex_real)

       call amrex_abeclaplacian_destroy(abeclap)
       call amrex_multigrid_destroy(multigrid)

       end do
    end if

    ! Since this problem has Neumann BC, solution + constant is also a
    ! solution.  So we are going to shift the solution by a constant
    ! for comparison with the "exact solution".
    npts = ba(0)%num_pts()
    avg1 = exact_solution(0) % sum()
    avg2 = solution(0) % sum()
    offset = (avg1-avg2)/npts
    do ilev = 0, max_level
       call solution(ilev)%plus(offset, icomp=1, ncomp=1, nghost=0)
    end do

    ! cannot trust Fortran compiler to do this correctly
    do ilev = 0, max_level
       do idim = 1, amrex_spacedim
          call amrex_multifab_destroy(beta(idim,ilev))
       end do
    end do
  end subroutine solve_abeclaplacian

  subroutine write_plotfile ()
    type(amrex_multifab) :: plotmf(0:max_level)
    type(amrex_string), allocatable :: varname(:)
    integer, dimension(0:max_level) :: steps, rr
    integer :: ilev, nc

    if (allocated(acoef)) then
       nc = 6
    else
       nc = 4
    end if
    allocate(varname(nc))

    call amrex_string_build(varname(1), "solution")
    call amrex_string_build(varname(2), "rhs")
    call amrex_string_build(varname(3), "exact_solution")
    call amrex_string_build(varname(4), "error")
    if (allocated(acoef)) then
       call amrex_string_build(varname(5), "acoef")
       call amrex_string_build(varname(6), "bcoef")
    end if

    do ilev = 0, max_level
       call amrex_multifab_build(plotmf(ilev), ba(ilev), dm(ilev), nc, 0)
       call plotmf(ilev) % copy(      solution(ilev), 1, 1, 1, 0)
       call plotmf(ilev) % copy(           rhs(ilev), 1, 2, 1, 0)
       call plotmf(ilev) % copy(exact_solution(ilev), 1, 3, 1, 0)
       call plotmf(ilev) % copy(      solution(ilev), 1, 4, 1, 0)
       call plotmf(ilev) % subtract(exact_solution(ilev),1,4,1,0)
       if (allocated(acoef)) then
          call plotmf(ilev) % copy(acoef(ilev), 1, 5, 1, 0)
          call plotmf(ilev) % copy(bcoef(ilev), 1, 6, 1, 0)
       end if
    end do

    steps = 1
    rr = ref_ratio

    call amrex_write_plotfile("plot", max_level+1, plotmf, varname, geom, 0._amrex_real, steps, rr)

    ! let's not realy on finalizer, which is feature not all compilers support properly.
    do ilev = 0, max_level
       call amrex_multifab_destroy(plotmf(ilev))
    end do
  end subroutine write_plotfile

end module mytest_module
