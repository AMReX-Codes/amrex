module solve_helmholtz

  use helmholtz 

  private
  public :: solve

contains  

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
            metric_term=.false., agglomeration=agglomeration, consolidation=consolidation)
       
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

       err = multigrid % solve(solution, rhs, 1.e-10_amrex_real, 0.0_amrex_real)

       call amrex_poisson_destroy(poisson)
       call amrex_multigrid_destroy(multigrid)

    else
       do ilev = 0, max_level

          call amrex_poisson_build(poisson, [geom(ilev)], [ba(ilev)], [dm(ilev)], &
               metric_term=.false., agglomeration=agglomeration, consolidation=consolidation)
       
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
    integer(c_long) :: npts
    real(amrex_real) :: err, avg1, avg2, offset
    type(amrex_multifab), allocatable :: beta(:,:)
    logical :: nodal(3)
    type(amrex_multifab) :: null

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
            metric_term=.false., agglomeration=agglomeration, consolidation=consolidation)

       call abeclap % set_maxorder(linop_maxorder)

       ! This is set up to have homogeneous Neumann BC
       call abeclap % set_domain_bc([amrex_lo_neumann, amrex_lo_neumann, amrex_lo_neumann], &
            &                       [amrex_lo_neumann, amrex_lo_neumann, amrex_lo_neumann])

       do ilev = 0, max_level
          ! for problem with pure homogeneous Neumann BC, we could pass an empty multifab
          call abeclap % set_level_bc(ilev, null)
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

       err = multigrid % solve(solution, rhs, 1.e-10_amrex_real, 0.0_amrex_real)

       call amrex_abeclaplacian_destroy(abeclap)
       call amrex_multigrid_destroy(multigrid)

    else
       do ilev = 0, max_level
          
          call amrex_abeclaplacian_build(abeclap, [geom(ilev)], [ba(ilev)], [dm(ilev)], &
               metric_term=.false., agglomeration=agglomeration, consolidation=consolidation)

       call abeclap % set_maxorder(linop_maxorder)

       ! This is set up to have homogeneous Neumann BC
       call abeclap % set_domain_bc([amrex_lo_neumann, amrex_lo_neumann, amrex_lo_neumann], &
            &                       [amrex_lo_neumann, amrex_lo_neumann, amrex_lo_neumann])

       if (ilev > 0) then
          ! use coarse level data to set up bc at corase/fine boundary
          call abeclap % set_coarse_fine_bc(solution(ilev-1), ref_ratio)
       end if

       ! for problem with pure homogeneous Neumann BC, we could pass an empty multifab
       call abeclap % set_level_bc(0, null)

       call abeclap % set_scalars(ascalar, bscalar)
       call abeclap % set_acoeffs(0, acoef(ilev))
       call abeclap % set_bcoeffs(0, beta(:,ilev))

       call amrex_multigrid_build(multigrid, abeclap)
       call multigrid % set_verbose(verbose)
       call multigrid % set_cg_verbose(cg_verbose)
       call multigrid % set_max_iter(max_iter)
       call multigrid % set_max_fmg_iter(max_fmg_iter)

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

end module solve_helmholtz
