module solve_helmholtz

 use amrex_amr_module
 use amr_data_module

 private
 public :: solve

contains  

  subroutine solve ()
    type(amrex_multigrid) :: multigrid
    integer :: ilev, idim
    integer(c_long) :: npts
    real(amrex_real) :: err, avg1, avg2, offset
    logical :: nodal(3)
    type(amrex_multifab) :: null
    integer :: nlevs,do_solve, total_levels

    integer :: myrank, ierr
    type(amrex_boxarray) :: batemp
    type(amrex_distromap) :: dmtemp
    logical :: periodic(3)

    nlevs=amrex_get_finest_level()

    ! For ABecLaplacian, the b coefficents are on faces

    do ilev = 0, nlevs
       call amrex_average_cellcenter_to_face(beta(:,ilev), bcoef(ilev), amrex_geom(ilev))
    end do

			
    call amrex_abeclaplacian_build(abeclap, amrex_geom(0:nlevs), phi_new(0:nlevs)%ba, phi_new(0:nlevs)%dm, &
            metric_term=.false., agglomeration=agglomeration, consolidation=consolidation)
	 
    call abeclap % set_maxorder(linop_maxorder)

    ! This is set up to have homogeneous Neumann BC
     call abeclap % set_domain_bc([amrex_lo_neumann, amrex_lo_neumann, amrex_lo_neumann], &
          &                       [amrex_lo_neumann, amrex_lo_neumann, amrex_lo_neumann])

     do ilev = 0, nlevs
         ! for problem with pure homogeneous Neumann BC, we could pass an empty multifab
         call abeclap % set_level_bc(ilev, solution(ilev))
     end do

     call abeclap % set_scalars(ascalar, bscalar)

     do ilev = 0, nlevs
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

    ! Since this problem has Neumann BC, solution + constant is also a
    ! solution.  So we are going to shift the solution by a constant
    ! for comparison with the "exact solution".
    npts = phi_new(0)%ba%num_pts()
    avg1 = exact_solution(0) % sum()
    avg2 = solution(0) % sum()
    offset = (avg1-avg2)/npts
    do ilev = 0, nlevs
       call solution(ilev)%plus(offset, icomp=1, ncomp=1, nghost=0)
    end do

  end subroutine solve

end module solve_helmholtz
