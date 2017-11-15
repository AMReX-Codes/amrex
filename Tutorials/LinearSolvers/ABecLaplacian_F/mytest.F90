
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
  
  integer, save :: prob_type = 1

  integer, save :: verbose = 2
  integer, save :: cg_verbose = 0
  integer, save :: max_iter = 100
  integer, save :: max_fmg_iter = 0
  integer, save :: linop_maxorder = 2
  logical, save :: agglomeration = .true.
  logical, save :: consolidation = .true.

  ! data
  type(amrex_geometry), allocatable :: geom(:)
  type(amrex_boxarray), allocatable :: ba(:)
  type(amrex_distromap), allocatable :: dm(:)
  type(amrex_multifab), allocatable :: solution(:)
  type(amrex_multifab), allocatable :: rhs(:)
  type(amrex_multifab), allocatable :: exact_solution(:)
  type(amrex_multifab), allocatable :: acoef(:)
  type(amrex_multifab), allocatable :: bcoef(:,:)
  

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

    call pp%query("max_level", max_level)
    call pp%query("ref_ratio", ref_ratio)
    call pp%query("n_cell", n_cell);
    call pp%query("max_grid_size", max_grid_size)

    call pp%query("composite_solve", composite_solve)

    call pp%query("prob_type", prob_type)

    call pp%query("verbose", verbose)
    call pp%query("cg_verbose", cg_verbose)
    call pp%query("max_iter", max_iter)
    call pp%query("max_fmg_iter", max_fmg_iter)
    call pp%query("linop_maxorder", linop_maxorder)
    call pp%query("agglomeration", agglomeration)
    call pp%query("consolidation", consolidation)

    call amrex_parmparse_destroy(pp)
  end subroutine init_parameters


  subroutine init_data ()
    allocate(geom(0:max_level))
    allocate(ba(0:max_level))
    allocate(dm(0:max_level))
    allocate(solution(0:max_level))
    allocate(rhs(0:max_level))
    allocate(exact_solution(0:max_level))
    if (prob_type .eq. 2) then
       allocate(acoef(0:max_level))
       allocate(bcoef(3,0:max_level))
    end if

    ! buid ....
  end subroutine init_data

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
          do idim = 1, 3
             call amrex_multifab_destroy(bcoef(idim,ilev))
          end do
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

    if (composite_solve) then

       call amrex_poisson_build(poisson, geom, ba, dm, metric_term=.false., &
            agglomeration=agglomeration, consolidation=consolidation)
       
       call poisson%set_maxorder(linop_maxorder)

       ! This is a 3d problem with Dirichlet BC
       call poisson%set_domain_bc([amrex_lo_dirichlet, amrex_lo_dirichlet, amrex_lo_dirichlet], &
                                  [amrex_lo_dirichlet, amrex_lo_dirichlet, amrex_lo_dirichlet])

       do ilev = 0, max_level
          ! Input: solution multifab's ghost cells at physical boundaries contain bc values.
          call poisson%set_level_bc(ilev, solution(ilev))
       end do

       call amrex_multigrid_build(multigrid,poisson)

       call amrex_poisson_destroy(poisson)
       call amrex_multigrid_destroy(multigrid)
    else

    end if
  end subroutine solve_poisson


  subroutine solve_abeclaplacian ()
  end subroutine solve_abeclaplacian

  subroutine write_plotfile ()
  end subroutine write_plotfile

end module mytest_module
