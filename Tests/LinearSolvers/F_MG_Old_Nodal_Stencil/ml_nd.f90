module ml_nd_module

  use bl_constants_module
  use bl_prof_module
  use mg_module
  use fabio_module
  use ml_layout_module
  use bndry_reg_module

  implicit none

!  private :: grid_res, grid_laplace_1d, grid_laplace_2d, grid_laplace_3d

contains

  subroutine ml_nd(mla,mgt,rh,full_soln,fine_mask,one_sided_ss,ref_ratio, &
                   do_diagnostics,rel_eps,abs_eps_in)

    use ml_norm_module        , only : ml_norm_inf
    use ml_restriction_module , only : ml_restriction, periodic_add_copy
    use ml_prolongation_module, only : ml_nodal_prolongation

    type(ml_layout), intent(in   ) :: mla
    type(mg_tower ), intent(inout) :: mgt(:)
    type( multifab), intent(inout) :: rh(:)
    type( multifab), intent(inout) :: full_soln(:)
    type(lmultifab), intent(in   ) :: fine_mask(:)
    type( multifab), intent(in   ) :: one_sided_ss(2:)
    integer        , intent(in   ) :: ref_ratio(:,:)
    integer        , intent(in   ) :: do_diagnostics 
    real(dp_t)     , intent(in   ) :: rel_eps

    real(dp_t)     , intent(in   ), optional :: abs_eps_in

    integer    :: nlevs

    type(multifab), allocatable  ::      soln(:)
    type(multifab), allocatable  ::        uu(:)
    type(multifab), allocatable  ::   uu_hold(:)
    type(multifab), allocatable  ::       res(:)
    type(multifab), allocatable  ::  temp_res(:)

    ! This is just a holder for a zero array at the same size as rh --
    !   we only need it to pass in to crse_fine_residual
    type(multifab), allocatable  ::   zero_rh(:)

    type(bndry_reg), allocatable :: brs_flx(:)

    type(box   ) :: pd,pdc
    type(layout) :: la
    integer :: n, dm
    integer :: mglev, mglev_crse, iter
    logical :: fine_converged

    real(dp_t) :: Anorm, bnorm, fac, tres, ttres, tres0, abs_eps, t1(3), t2(3)
    real(dp_t) :: stime, bottom_solve_time

    character(len=3)          :: number
    character(len=20)         :: filename

    logical nodal(get_dim(rh(1)))

    type(bl_prof_timer), save :: bpt

    call build(bpt, "ml_nd")

    dm                = get_dim(rh(1))
    stime             = parallel_wtime()
    nodal             = .True.
    nlevs             = mla%nlevel
    bottom_solve_time = zero

    if (nghost(rh(nlevs)) .ne. 1) &
        call bl_error("ml_nd requires one ghost cell for the RHS");

    if ( present(abs_eps_in) ) then
       abs_eps = abs_eps_in
    else
       abs_eps = mgt(nlevs)%abs_eps
    end if

    allocate(soln(nlevs), uu(nlevs), uu_hold(2:nlevs-1), res(nlevs))
    allocate(temp_res(nlevs))
    allocate(brs_flx(2:nlevs))
    allocate(zero_rh(2:nlevs))

    do n = 2,nlevs-1
       la = mla%la(n)
       call multifab_build( uu_hold(n), la, 1, 1, nodal)
       call setval( uu_hold(n), ZERO,all=.true.)
    end do

    do n = nlevs, 1, -1

       la = mla%la(n)
       call multifab_build(    soln(n), la, 1, 1, nodal)
       call multifab_build(      uu(n), la, 1, 1, nodal)
       call multifab_build(     res(n), la, 1, 1, nodal)
       call multifab_build(temp_res(n), la, 1, 1, nodal)
       call setval(    soln(n), ZERO,all=.true.)
       call setval(      uu(n), ZERO,all=.true.)
       call setval(     res(n), ZERO,all=.true.)
       call setval(temp_res(n), ZERO,all=.true.)

       if ( n == 1 ) exit

       call multifab_build(zero_rh(n), la, 1, nghost(rh(n)), nodal)
       call setval(zero_rh(n), ZERO,all=.true.)

       ! Build the (coarse resolution) flux registers to be used in computing
       !  the residual at a non-finest AMR level.

       pdc = layout_get_pd(mla%la(n-1))
       call bndry_reg_rr_build_nd(brs_flx(n), la, ref_ratio(n-1,:), pdc, nodal)

    end do
    !
    ! Let's elide some reductions by doing these reductions together.
    !
    bnorm = ml_norm_inf(rh,fine_mask,local=.true.)

    Anorm = stencil_norm(mgt(nlevs)%ss(mgt(nlevs)%nlevels),local=.true.)
    do n = 1, nlevs-1
       Anorm = max(stencil_norm(mgt(n)%ss(mgt(n)%nlevels),fine_mask(n),local=.true.), Anorm)
    end do

    do n = nlevs,1,-1
       mglev = mgt(n)%nlevels
       call compute_defect(mgt(n)%ss(mglev),res(n),rh(n),full_soln(n),mgt(n)%mm(mglev), &
                           mgt(n)%stencil_type, mgt(n)%lcross, mgt(n)%uniform_dh)
    end do

    do n = 1,nlevs
       call multifab_copy(rh(n),res(n),ng=nghost(rh(n)))
    end do

    t1(1) = bnorm
    t1(2) = Anorm
    t1(3) = ml_norm_inf(rh,fine_mask,local=.true.)

    call parallel_reduce(t2, t1, MPI_MAX)

    bnorm = t2(1)
    Anorm = t2(2)
    tres0 = t2(3)

    if ( parallel_IOProcessor() .and. mgt(nlevs)%verbose > 0 ) then
       write(unit=*, &
             fmt='("F90mg: Initial rhs                  = ",g15.8)') bnorm
       write(unit=*, &
             fmt='("F90mg: Initial residual (resid0)    = ",g15.8)') tres0
    end if

    bnorm = max(bnorm, tres0)  
    ! wqz.  Otherwise, we may sometimes have trouble if tres0 is much larger than bnorm

    ! ****************************************************************************

    fine_converged = .false.

    if ( ml_converged(res, fine_mask, bnorm, &
                      rel_eps, abs_eps, mgt(nlevs)%verbose) ) then
       if ( parallel_IOProcessor() .and. mgt(nlevs)%verbose > 0 ) &
            write(unit=*, fmt='("F90mg: No iterations needed ")')

    else   ! Not already converged

     do iter = 1, mgt(nlevs)%max_iter

       if ( (iter .eq. 1) .or. fine_converged ) then
          if ( ml_converged(res, fine_mask, bnorm, rel_eps, abs_eps, mgt(nlevs)%verbose) ) exit
       end if

       ! Set: uu = 0
       do n = 1,nlevs
          call setval(uu(n), ZERO, all=.true.)
       end do

       ! Set: uu_hold = 0
       do n = 2,nlevs-1
          call setval(uu_hold(n), ZERO, all=.true.)
       end do

       !   Down the V-cycle
       do n = nlevs,1,-1

          mglev = mgt(n)%nlevels

          if ( do_diagnostics == 1 ) then
             tres = norm_inf(res(n))
             if ( parallel_ioprocessor() ) then
                print *,'DWN: RES BEFORE GSRB AT LEVEL ',n, tres
             end if
          end if

          ! Relax ...
          if (n > 1) then
             call mini_cycle(mgt(n), mglev, mgt(n)%ss(mglev), &
                  uu(n), res(n), mgt(n)%mm(mglev), mgt(n)%nu1, mgt(n)%nu2)
          else 
             call mg_tower_cycle(mgt(n), mgt(n)%cycle_type, mglev, mgt(n)%ss(mglev), &
                  uu(n), res(n), mgt(n)%mm(mglev), mgt(n)%nu1, mgt(n)%nu2, &
                  bottom_solve_time = bottom_solve_time)
          end if

          ! Add: soln += uu
          call plus_plus(soln(n),uu(n))

          if (n > 1) then
             mglev_crse = mgt(n-1)%nlevels

             ! Compute COARSE Res = Rh - Lap(Soln)
             call compute_defect(mgt(n-1)%ss(mglev_crse),res(n-1), &
                  rh(n-1),soln(n-1),mgt(n-1)%mm(mglev_crse), &
                  mgt(n-1)%stencil_type, mgt(n-1)%lcross, mgt(n-1)%uniform_dh)

             if ( dm .eq. 3 ) then
                fac = (8.0_dp_t)**(ref_ratio(n-1,1)/2)
                call multifab_mult_mult_s(res(n-1),fac,nghost(res(n-1)))
             else if ( dm .eq. 2 ) then
                fac = (4.0_dp_t)**(ref_ratio(n-1,1)/2)
                call multifab_mult_mult_s(res(n-1),fac,nghost(res(n-1)))
             end if

             ! Compute FINE Res = Res - Lap(uu)
             mglev = mgt(n)%nlevels
             call compute_defect(mgt(n)%ss(mglev), temp_res(n), &
                  res(n),uu(n),mgt(n)%mm(mglev), &
                  mgt(n)%stencil_type, mgt(n)%lcross, mgt(n)%uniform_dh)
             call multifab_copy(res(n),temp_res(n),ng=nghost(res(n)))

             if ( do_diagnostics == 1 ) then
                tres = norm_inf(res(n))
                if ( parallel_ioprocessor() ) &
                   print *,'DWN: RES AFTER  GSRB AT LEVEL ',n, tres 
             end if

             ! Restrict FINE Res to COARSE Res
             call ml_restriction(res(n-1), res(n), mgt(n)%mm(mglev),&
                                 mgt(n-1)%mm(mglev_crse), ref_ratio(n-1,:))

             ! Compute CRSE-FINE Res = Rh - Lap(Soln)
             pdc = layout_get_pd(mla%la(n-1))
             call crse_fine_residual_nodal(n,mgt,brs_flx(n),res(n-1),zero_rh(n),temp_res(n),temp_res(n-1), &
                  soln(n-1),soln(n),one_sided_ss(n),ref_ratio(n-1,:),pdc)

             ! Copy u_hold = uu
             if (n < nlevs) call multifab_copy(uu_hold(n),uu(n),ng=nghost(uu(n)))

             ! Set: uu = 0
             call setval(uu(n),ZERO,all=.true.)

             if ( dm .eq. 3 ) then
                fac = 1.0_dp_t / (8.0_dp_t)**(ref_ratio(n-1,1)/2)
                call multifab_mult_mult_s(res(n-1),fac,nghost(res(n-1)))
             else if ( dm .eq. 2 ) then
                fac = 1.0_dp_t / (4.0_dp_t)**(ref_ratio(n-1,1)/2)
                call multifab_mult_mult_s(res(n-1),fac,nghost(res(n-1)))
             end if

          else

             if (do_diagnostics == 1 ) then
                call compute_defect(mgt(n)%ss(mglev),temp_res(n), res(n),uu(n),mgt(n)%mm(mglev), &
                                    mgt(n)%stencil_type, mgt(n)%lcross, mgt(n)%uniform_dh)
                tres = norm_inf(temp_res(n))
                if ( parallel_ioprocessor() ) then
                   print *,'DWN: RES AFTER  GSRB AT LEVEL ',n, tres
                end if
             end if

          end if

        end do

        !   Back up the V-cycle
        do n = 2, nlevs

          pd = layout_get_pd(mla%la(n))
          mglev = mgt(n)%nlevels

          ! Interpolate uu from coarser level
          if (iter == 1) call plus_plus(uu(n-1),  full_soln(n-1))
          call ml_nodal_prolongation(uu(n), uu(n-1), ref_ratio(n-1,:))
          if (iter == 1) call sub_sub(uu(n-1), full_soln(n-1))

          ! Subtract: uu -= full_soln
          !     Must do this in order to remove interpolated full_soln...
          if (iter == 1) call sub_sub(uu(n),full_soln(n))

          ! Add: soln += uu
          call plus_plus(soln(n), uu(n), nghost(uu(n)))

          ! Add: uu_hold += uu 
          if (n < nlevs) call plus_plus(uu_hold(n), uu(n), nghost(uu(n)))

          ! Compute Res = Res - Lap(uu)
          call compute_defect(mgt(n)%ss(mglev),temp_res(n),res(n),uu(n),mgt(n)%mm(mglev), &
                              mgt(n)%stencil_type, mgt(n)%lcross, mgt(n)%uniform_dh)
          call multifab_copy(res(n),temp_res(n),ng=nghost(res(n)))

          if ( do_diagnostics == 1 ) then
             tres = norm_inf(res(n))
             if ( parallel_ioprocessor() ) then
                print *,'UP : RES BEFORE GSRB AT LEVEL ',n, tres
             end if
          end if

          ! Set: uu = 0
          call setval(uu(n),ZERO,all=.true.)

          ! Relax ...
          call mini_cycle(mgt(n), mglev, mgt(n)%ss(mglev), &
               uu(n), res(n), mgt(n)%mm(mglev), mgt(n)%nu1, mgt(n)%nu2)

          ! Compute Res = Res - Lap(uu)
          call compute_defect(mgt(n)%ss(mglev),temp_res(n),res(n),uu(n),mgt(n)%mm(mglev), &
                              mgt(n)%stencil_type, mgt(n)%lcross, mgt(n)%uniform_dh)
          call multifab_copy(res(n),temp_res(n),ng=nghost(res(n)))

          if ( do_diagnostics == 1 ) then
             tres = norm_inf(res(n))
             if ( parallel_ioprocessor() ) then
                print *,'UP : RES AFTER  GSRB AT LEVEL ',n, tres
                if (n == nlevs) print *,' '
             end if
          end if

          ! Add: soln += uu
          call plus_plus(soln(n), uu(n), nghost(uu(n)))

          ! Add: uu += uu_hold so that it will be interpolated too.
          if (n < nlevs) call plus_plus(uu(n), uu_hold(n), nghost(uu(n)))

       end do

       !    Inject the solution to the coarser grids.
       do n = nlevs,2,-1
          mglev      = mgt(n)%nlevels
          mglev_crse = mgt(n-1)%nlevels
          call ml_restriction(soln(n-1), soln(n), mgt(n)%mm(mglev), &
               mgt(n-1)%mm(mglev_crse), ref_ratio(n-1,:), inject = .true.)
       end do

        do n = 1,nlevs
           call multifab_fill_boundary(soln(n), cross = mgt(n)%lcross)
        end do

       !    Optimization so don't have to do multilevel convergence test each time

       !    Compute the residual on just the finest level
        n = nlevs
        mglev = mgt(n)%nlevels
        call compute_defect(mgt(n)%ss(mglev),res(n),rh(n),soln(n),mgt(n)%mm(mglev), &
                            mgt(n)%stencil_type, mgt(n)%lcross, mgt(n)%uniform_dh)

        if ( ml_fine_converged(res, bnorm, rel_eps, abs_eps) ) then

          fine_converged = .true.

          !      Compute the residual on every level
          do n = 1,nlevs-1
             mglev = mgt(n)%nlevels
             call compute_defect(mgt(n)%ss(mglev),res(n),rh(n),soln(n),mgt(n)%mm(mglev), &
                                 mgt(n)%stencil_type, mgt(n)%lcross, mgt(n)%uniform_dh)
          end do

          do n = nlevs,2,-1
             !  Restrict the finer residual onto the coarser grid
             mglev      = mgt(n  )%nlevels
             mglev_crse = mgt(n-1)%nlevels
             if ( dm .eq. 3 ) then
                fac = (8.0_dp_t)**(ref_ratio(n-1,1)/2)
                call multifab_mult_mult_s(res(n-1),fac,nghost(res(n-1)))
             else if ( dm .eq. 2 ) then
                fac = (4.0_dp_t)**(ref_ratio(n-1,1)/2)
                call multifab_mult_mult_s(res(n-1),fac,nghost(res(n-1)))
             end if
             call ml_restriction(res(n-1), res(n), mgt(n)%mm(mglev),&
                  mgt(n-1)%mm(mglev_crse), ref_ratio(n-1,:))

             !  Compute the coarse-fine residual at coarse-fine nodes
             pdc = layout_get_pd(mla%la(n-1))
             call crse_fine_residual_nodal(n,mgt,brs_flx(n),res(n-1), &
                  zero_rh(n),temp_res(n),temp_res(n-1), &
                  soln(n-1),soln(n),one_sided_ss(n),ref_ratio(n-1,:),pdc)
             if ( dm .eq. 3 ) then
                fac = 1.0_dp_t / (8.0_dp_t)**(ref_ratio(n-1,1)/2)
                call multifab_mult_mult_s(res(n-1),fac,nghost(res(n-1)))
             else if ( dm .eq. 2 ) then
                fac = 1.0_dp_t / (4.0_dp_t)**(ref_ratio(n-1,1)/2)
                call multifab_mult_mult_s(res(n-1),fac,nghost(res(n-1)))
             end if
          end do

          if ( mgt(nlevs)%verbose > 1 ) then
             do n = 1,nlevs
                ttres = norm_inf(res(n),local=.true.)
                call parallel_reduce(tres, ttres, MPI_MAX, proc = parallel_IOProcessorNode())
                if ( parallel_IOProcessor() ) then
                   write(unit=*, fmt='("F90mg: Iteration   ",i3," Lev ",i1," resid/resid0 = ",g15.8)') &
                        iter,n,tres/tres0
                end if
             end do
          end if

       else

          fine_converged = .false.
          if ( mgt(nlevs)%verbose > 1 ) then

             if ( .false. ) then
                !
                ! Some debugging code I want to keep around for a while.
                ! I don't want to have to recreate this all the time :-)
                !
                write(number,fmt='(i3.3)') iter
                filename = 'res_fine_iter=' // number
                call fabio_write(res(nlevs), 'debug', trim(filename))
             end if

             ttres = norm_inf(res(nlevs),local=.true.)
             call parallel_reduce(tres, ttres, MPI_MAX, proc = parallel_IOProcessorNode())
             if ( parallel_IOProcessor() ) then
                write(unit=*, fmt='("F90mg: Iteration   ",i3," Fine  resid/resid0 = ",g15.8)') iter,tres/tres0
             end if
          end if

       end if

     end do

     iter = iter-1
     if (iter < mgt(nlevs)%max_iter) then
        if ( mgt(nlevs)%verbose > 0 ) then
           !
           ! Consolidate these reductions.
           !
           t1(1) = ml_norm_inf(res,fine_mask,local=.true.) 
           t1(2) = (parallel_wtime() - stime)
           t1(3) = bottom_solve_time

           call parallel_reduce(t2, t1, MPI_MAX, proc = parallel_IOProcessorNode())

           if ( parallel_IOProcessor() ) then
              if ( tres0 .gt. 0.0_dp_t) then
                 write(unit=*, fmt='("F90mg: Final Iter. ",i3," resid/resid0 = ",g15.8)') iter,t2(1)/tres0
                 write(unit=*, fmt='("F90mg: Solve time: ",g13.6, " Bottom Solve time: ", g13.6)') t2(2), t2(3)
                 write(unit=*, fmt='("")')
              else
                 write(unit=*, fmt='("F90mg: Final Iter. ",i3," resid/resid0 = ",g15.8)') iter,0.0_dp_t
                 write(unit=*, fmt='("F90mg: Solve time: ",g13.6, " Bottom Solve time: ", g13.6)') t2(2), t2(3)
                 write(unit=*, fmt='("")')
              end if
           end if
        end if
     else
        if (mgt(nlevs)%abort_on_max_iter) then
           call bl_error("Multigrid Solve: failed to converge in max_iter iterations")
        end if
     end if

    ! Add: full_soln += soln
    do n = 1,nlevs
       call plus_plus(full_soln(n),soln(n))
    end do

    end if

    ! ****************************************************************************

    do n = 2,nlevs-1
       call multifab_destroy(uu_hold(n))
    end do

    do n = nlevs, 1, -1
       call multifab_destroy(    soln(n))
       call multifab_destroy(      uu(n))
       call multifab_destroy(     res(n))
       call multifab_destroy(temp_res(n))
       if ( n == 1 ) exit
       call multifab_destroy( zero_rh(n))
       call bndry_reg_destroy(brs_flx(n))
    end do

    call destroy(bpt)

  contains

    subroutine crse_fine_residual_nodal(n,mgt,brs_flx,crse_res,fine_rhs,temp_res,temp_crse_res, &
         crse_soln,fine_soln,one_sided_ss,ref_ratio,pdc)

      use nodal_interface_stencil_module , only : ml_crse_contrib, ml_fine_contrib

      integer        , intent(in   ) :: n
      type(mg_tower) , intent(inout) :: mgt(:)
      type(bndry_reg), intent(inout) :: brs_flx
      type(multifab) , intent(inout) :: crse_res
      type(multifab) , intent(in   ) :: fine_rhs
      type(multifab) , intent(inout) :: temp_res
      type(multifab) , intent(inout) :: temp_crse_res
      type(multifab) , intent(inout) :: crse_soln
      type(multifab) , intent(inout) :: fine_soln
      type(multifab) , intent(in   ) :: one_sided_ss
      integer        , intent(in   ) :: ref_ratio(:)
      type(box)      , intent(in   ) :: pdc

      integer :: i,dm,mglev_crse,mglev_fine

      mglev_crse = mgt(n-1)%nlevels
      mglev_fine = mgt(n  )%nlevels
      dm        = get_dim(temp_res)

      !    Compute the fine contributions at faces, edges and corners.

      !    First compute a residual which only takes contributions from the
      !       grid on which it is calculated.

      if (mgt(n)%lcross) then
        call grid_res(one_sided_ss,temp_res, &
             fine_rhs,fine_soln,mgt(n)%mm(mglev_fine),mgt(n)%face_type, &
             mgt(n)%lcross, mgt(n)%uniform_dh)
      else
        call grid_res(mgt(n)%ss(mglev_fine),temp_res, &
             fine_rhs,fine_soln,mgt(n)%mm(mglev_fine),mgt(n)%face_type, &
             mgt(n)%lcross, mgt(n)%uniform_dh)
      end if

      !    Zero out the flux registers which will hold the fine contributions
      call bndry_reg_setval(brs_flx, ZERO, all = .true.)

      do i = 1,dm
         call ml_fine_contrib(brs_flx%bmf(i,0), &
              temp_res,mgt(n)%mm(mglev_fine),ref_ratio,pdc,-i)
         call ml_fine_contrib(brs_flx%bmf(i,1), &
              temp_res,mgt(n)%mm(mglev_fine),ref_ratio,pdc,+i)
      end do

!     Compute the crse contributions at edges and corners and add to fine contributions
!        in temp_crse_res (need to do this in a temporary for periodic issues)
      call setval(temp_crse_res,ZERO,all=.true.)
      do i = 1,dm
         call ml_crse_contrib(temp_crse_res, brs_flx%bmf(i,0), crse_soln, &
              mgt(n-1)%ss(mgt(n-1)%nlevels), &
              mgt(n-1)%mm(mglev_crse), &
              mgt(n  )%mm(mglev_fine), &
              pdc,ref_ratio, -i)
         call ml_crse_contrib(temp_crse_res, brs_flx%bmf(i,1), crse_soln, &
              mgt(n-1)%ss(mgt(n-1)%nlevels), &
              mgt(n-1)%mm(mglev_crse), &
              mgt(n  )%mm(mglev_fine), &
              pdc,ref_ratio, +i)
      end do

!     Add to res(n-1).
      call plus_plus(crse_res,temp_crse_res)

      call periodic_add_copy(crse_res,temp_crse_res,synced=.true.)

!     Clear temp_crse_res (which is temp_res(n-1) from calling routine) just in case...
      call setval(temp_crse_res,ZERO,all=.true.)

    end subroutine crse_fine_residual_nodal

    function ml_fine_converged(res, bnorm, rel_eps, abs_eps) result(r)
      logical :: r
      type(multifab), intent(in) :: res(:)
      real(dp_t), intent(in) :: rel_eps, abs_eps, bnorm
      real(dp_t) :: ni_res
      integer    :: nlevs
      nlevs = size(res)
      ni_res = norm_inf(res(nlevs))
      r = ( ni_res <= rel_eps*(bnorm) .or. ni_res <= abs_eps )
    end function ml_fine_converged

    function ml_converged(res, mask, bnorm, rel_eps, abs_eps, verbose) result(r)
      use ml_norm_module, only : ml_norm_inf
      logical :: r
      integer :: verbose
      type(multifab), intent(in) :: res(:)
      type(lmultifab), intent(in) :: mask(:)
      real(dp_t), intent(in) :: rel_eps, abs_eps, bnorm
      real(dp_t) :: ni_res
      ni_res = ml_norm_inf(res, mask)
      r = ( ni_res <= rel_eps*(bnorm) .or. ni_res <= abs_eps )
      if ( r .and. parallel_IOProcessor() .and. verbose > 1) then
         if (ni_res <= rel_eps*bnorm) then
            print *,'Converged res < rel_eps*bnorm '
         else if (ni_res <= abs_eps) then
            print *,'Converged res < abs_eps '
         end if
      end if
    end function ml_converged

  end subroutine ml_nd

  subroutine grid_res(ss, dd, ff, uu, mm, face_type, lcross, uniform_dh)

    type(multifab), intent(in)    :: ff, ss
    type(multifab), intent(inout) :: dd, uu
    type(imultifab), intent(in)   :: mm
    integer, intent(in) :: face_type(:,:,:)
    logical, intent(in) :: lcross
    logical, intent(in) :: uniform_dh

    integer :: i, n
    real(kind=dp_t), pointer :: dp(:,:,:,:)
    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: up(:,:,:,:)
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    integer        , pointer :: mp(:,:,:,:)
    integer :: dm,nodal_ng
    type(bl_prof_timer), save :: bpt

    call build(bpt, "grid_res")

    nodal_ng = 0; if ( nodal_q(uu) ) nodal_ng = 1

    dm = get_dim(uu)

    call multifab_fill_boundary(uu, cross = lcross)

    do i = 1, nfabs(uu)
       dp  => dataptr(dd, i)
       fp  => dataptr(ff, i)
       up  => dataptr(uu, i)
       sp  => dataptr(ss, i)
       mp  => dataptr(mm, i)
       do n = 1, ncomp(uu)
          select case(dm)
          case (1)
             call grid_laplace_1d(sp(:,:,1,1), dp(:,1,1,n), fp(:,1,1,n), up(:,1,1,n), &
                                  mp(:,1,1,1),nghost(uu))
          case (2)
             call grid_laplace_2d(sp(:,:,:,1), dp(:,:,1,n), fp(:,:,1,n), up(:,:,1,n), &
                                  mp(:,:,1,1), nghost(uu), face_type(i,:,:))
          case (3)
             call grid_laplace_3d(sp(:,:,:,:), dp(:,:,:,n), fp(:,:,:,n), up(:,:,:,n), &
                                  mp(:,:,:,1), nghost(uu), face_type(i,:,:), uniform_dh)
          end select
       end do
    end do

    call destroy(bpt)

  end subroutine grid_res

  subroutine grid_laplace_1d(ss, dd, ff, uu, mm, ng)

    use impose_neumann_bcs_module

    integer           , intent(in   ) :: ng
    real (kind = dp_t), intent(in   ) :: ss(0:,:)
    real (kind = dp_t), intent(inout) :: dd(0:)
    real (kind = dp_t), intent(in   ) :: ff(0:)
    integer,            intent(in   ) :: mm(:)
    real (kind = dp_t), intent(inout) :: uu(1-ng:)

    integer :: i,nx,lo(1)
    nx = size(ss,dim=2)-1

    lo = 1
    call impose_neumann_bcs_1d(uu,mm,lo,ng)

    i = 1
    dd(i) = HALF*ff(i) - (ss(0,i)*uu(i) + ss(1,i)*(uu(i+1)-uu(i)))

    do i = 2,nx
      dd(i) = ff(i) - &
              (ss(0,i)*uu(i) + ss(1,i) * uu(i+1) + ss(2,i) * uu(i-1))
    end do

    i = nx+1
    dd(i) = HALF*ff(i) - (ss(2,i)*(uu(i-1)-uu(i)))

  end subroutine grid_laplace_1d

  subroutine grid_laplace_2d(ss, dd, ff, uu, mm, ng, face_type)

    use bc_module
    use impose_neumann_bcs_module

    integer           , intent(in   ) :: ng
    real (kind = dp_t), intent(in   ) :: ss(0:,:,:)
    real (kind = dp_t), intent(inout) :: dd(0:,0:) 
    real (kind = dp_t), intent(in   ) :: ff(0:,0:)
    real (kind = dp_t), intent(inout) :: uu(1-ng:,1-ng:)
    integer,            intent(in   ) :: mm(:,:)
    integer,            intent(in ) :: face_type(:,:)
    integer :: i,j,nx,ny,lo(2)
    integer :: istart,iend,jstart,jend

    nx = size(ss,dim=2)-1
    ny = size(ss,dim=3)-1

    lo = 1
    call impose_neumann_bcs_2d(uu,mm,lo,ng)

    if (face_type(1,1) == BC_NEU) then
      istart = 1
    else
      istart = 2
    end if
    if (face_type(1,2) == BC_NEU) then
      iend = nx+1
    else
      iend = nx
    end if
    if (face_type(2,1) == BC_NEU) then
      jstart = 1
    else
      jstart = 2
    end if
    if (face_type(2,2) == BC_NEU) then
      jend = ny+1
    else
      jend = ny
    end if

    if (size(ss,dim=1) .eq. 9) then
!     Corners
      i = 1
      j = 1
      dd(i,j) = ss(8,i,j)*(uu(i+1,j+1)+HALF*uu(i,j+1)+HALF*uu(i+1,j)-TWO*uu(i,j))
      dd(i,j) = FOURTH*ff(i,j) - dd(i,j)

      i = 1
      j = ny+1
      dd(i,j) = ss(3,i,j)*(uu(i+1,j-1)+HALF*uu(i,j-1)+HALF*uu(i+1,j)-TWO*uu(i,j))
      dd(i,j) = FOURTH*ff(i,j) - dd(i,j)
 
      i = nx+1
      j = 1
      dd(i,j) = ss(6,i,j)*(uu(i-1,j+1)+HALF*uu(i,j+1)+HALF*uu(i-1,j)-TWO*uu(i,j))
      dd(i,j) = FOURTH*ff(i,j) - dd(i,j)

      i = nx+1
      j = ny+1
      dd(i,j) = ss(1,i,j)*(uu(i-1,j-1)+HALF*uu(i,j-1)+HALF*uu(i-1,j)-TWO*uu(i,j))
      dd(i,j) = FOURTH*ff(i,j) - dd(i,j)
 
!     Lo-x edge
      i = 1
      do j = jstart,jend
         dd(i,j) = ss(3,i,j)*(uu(i+1,j-1)+HALF*uu(i,j-1)+HALF*uu(i+1,j)-TWO*uu(i,j)) &
                  +ss(8,i,j)*(uu(i+1,j+1)+HALF*uu(i,j+1)+HALF*uu(i+1,j)-TWO*uu(i,j))
         dd(i,j) = HALF*ff(i,j) - dd(i,j)
      end do

!     Hi-x edge
      i = nx+1
      do j = jstart,jend
         dd(i,j) = ss(1,i,j)*(uu(i-1,j-1)+HALF*uu(i,j-1)+HALF*uu(i-1,j)-TWO*uu(i,j)) &
                  +ss(6,i,j)*(uu(i-1,j+1)+HALF*uu(i,j+1)+HALF*uu(i-1,j)-TWO*uu(i,j))
         dd(i,j) = HALF*ff(i,j) - dd(i,j)
      end do
 
!     Lo-y edge
      j = 1
      do i = istart,iend
         dd(i,j) = ss(6,i,j)*(uu(i-1,j+1)+HALF*uu(i,j+1)+HALF*uu(i-1,j)-TWO*uu(i,j)) &
                  +ss(8,i,j)*(uu(i+1,j+1)+HALF*uu(i,j+1)+HALF*uu(i+1,j)-TWO*uu(i,j))
         dd(i,j) = HALF*ff(i,j) - dd(i,j)
      end do

!     Hi-y edge
      j = ny+1
      do i = istart,iend
         dd(i,j) = ss(1,i,j)*(uu(i-1,j-1)+HALF*uu(i,j-1)+HALF*uu(i-1,j)-TWO*uu(i,j)) &
                  +ss(3,i,j)*(uu(i+1,j-1)+HALF*uu(i,j-1)+HALF*uu(i+1,j)-TWO*uu(i,j))
         dd(i,j) = HALF*ff(i,j) - dd(i,j)
      end do
 
!     Interior
      do j = jstart,jend
      do i = istart,iend
         dd(i,j) = ss(0,i,j)*uu(i,j) + ss(1,i,j) * uu(i-1,j-1) &
                                     + ss(2,i,j) * uu(i  ,j-1) &
                                     + ss(3,i,j) * uu(i+1,j-1) &
                                     + ss(4,i,j) * uu(i-1,j  ) &
                                     + ss(5,i,j) * uu(i+1,j  ) &
                                     + ss(6,i,j) * uu(i-1,j+1) &
                                     + ss(7,i,j) * uu(i  ,j+1) &
                                     + ss(8,i,j) * uu(i+1,j+1)
         dd(i,j) = ff(i,j) - dd(i,j)
      end do
      end do

    else if (size(ss,dim=1) .eq. 5) then

!     Corners
      i = 1
      j = 1
      dd(i,j) = ss(1,i,j)*(uu(i+1,j)-uu(i,j)) + ss(3,i,j)*(uu(i,j+1)-uu(i,j))
      dd(i,j) = FOURTH*ff(i,j) - dd(i,j)

      i = 1
      j = ny+1
      dd(i,j) = ss(1,i,j)*(uu(i+1,j)-uu(i,j)) + ss(4,i,j)*(uu(i,j-1)-uu(i,j))
      dd(i,j) = FOURTH*ff(i,j) - dd(i,j)

      i = nx+1
      j = 1
      dd(i,j) = ss(2,i,j)*(uu(i-1,j)-uu(i,j)) + ss(3,i,j)*(uu(i,j+1)-uu(i,j))
      dd(i,j) = FOURTH*ff(i,j) - dd(i,j)

      i = nx+1
      j = ny+1
      dd(i,j) = ss(2,i,j)*(uu(i-1,j)-uu(i,j)) + ss(4,i,j)*(uu(i,j-1)-uu(i,j))
      dd(i,j) = FOURTH*ff(i,j) - dd(i,j)

!     Lo-x edge
      i = 1
      do j = jstart,jend
         dd(i,j) =  ss(1,i,j)*(uu(i+1,j)-uu(i,j)) &
                   +ss(3,i,j)*(uu(i,j+1)-uu(i,j)) &
                   +ss(4,i,j)*(uu(i,j-1)-uu(i,j)) 
         dd(i,j) = HALF*ff(i,j) - dd(i,j)
      end do

!     Hi-x edge
      i = nx+1
      do j = jstart,jend
         dd(i,j) = ss(2,i,j)*(uu(i-1,j)-uu(i,j)) &
                  +ss(3,i,j)*(uu(i,j+1)-uu(i,j)) &
                  +ss(4,i,j)*(uu(i,j-1)-uu(i,j)) 
         dd(i,j) = HALF*ff(i,j) - dd(i,j)
      end do
 
!     Lo-y edge
      j = 1
      do i = istart,iend
         dd(i,j) = ss(3,i,j)*(uu(i,j+1)-uu(i,j)) &
                  +ss(1,i,j)*(uu(i+1,j)-uu(i,j)) &
                  +ss(2,i,j)*(uu(i-1,j)-uu(i,j)) 
         dd(i,j) = HALF*ff(i,j) - dd(i,j)
      end do

!     Hi-y edge
      j = ny+1
      do i = istart,iend
         dd(i,j) = ss(4,i,j)*(uu(i,j-1)-uu(i,j)) &
                  +ss(1,i,j)*(uu(i+1,j)-uu(i,j)) &
                  +ss(2,i,j)*(uu(i-1,j)-uu(i,j)) 
         dd(i,j) = HALF*ff(i,j) - dd(i,j)
      end do
 
!     Interior
      do j = jstart,jend
      do i = istart,iend
         dd(i,j) = ss(0,i,j)*uu(i,j) + ss(1,i,j) * uu(i+1,j  ) &
                                     + ss(2,i,j) * uu(i-1,j  ) &
                                     + ss(3,i,j) * uu(i  ,j+1) &
                                     + ss(4,i,j) * uu(i  ,j-1)  
         dd(i,j) = ff(i,j) - dd(i,j)
      end do
      end do

    end if

  end subroutine grid_laplace_2d

  subroutine grid_laplace_3d(ss, dd, ff, uu, mm, ng, face_type, uniform_dh)
    use bc_module
    use impose_neumann_bcs_module
    integer, intent(in) :: ng
    real (kind = dp_t), intent(in   ) :: ff(0:,0:,0:)
    real (kind = dp_t), intent(inout) :: uu(1-ng:,1-ng:,1-ng:)
    real (kind = dp_t), intent(inout) :: dd(0:,0:,0:)
    real (kind = dp_t), intent(in   ) :: ss(0:,:,:,:)
    integer,            intent(in   ) :: mm(:,:,:)
    integer,            intent(in   ) :: face_type(:,:)
    logical,            intent(in)    :: uniform_dh

    integer :: i, j, k, lo(3)
    integer :: istart, iend, jstart, jend, kstart, kend
    integer :: nx, ny, nz

    nx = size(ss,dim=2)-1
    ny = size(ss,dim=3)-1
    nz = size(ss,dim=4)-1

    lo = 1
    call impose_neumann_bcs_3d(uu,mm,lo,ng)

    if (face_type(1,1) == BC_NEU) then
      istart = 1
    else
      istart = 2
    end if
    if (face_type(1,2) == BC_NEU) then
      iend = nx+1
    else
      iend = nx
    end if
    if (face_type(2,1) == BC_NEU) then
      jstart = 1
    else
      jstart = 2
    end if
    if (face_type(2,2) == BC_NEU) then
      jend = ny+1
    else
      jend = ny
    end if
    if (face_type(3,1) == BC_NEU) then
      kstart = 1
    else
      kstart = 2
    end if
    if (face_type(3,2) == BC_NEU) then
      kend = nz+1
    else
      kend = nz
    end if

    if (size(ss,dim=1) .eq. 27 .or. size(ss,dim=1) .eq. 21) then
       !
       !     Corners
       !
      i = 1
      j = 1
      k = 1
      dd(i,j,k) = ss(20,i,j,k)*(uu(i+1,j+1,k+1) + uu(i+1,j+1,k  ) &
                               +uu(i+1,j  ,k+1) + uu(i  ,j+1,k+1) &
                         - FOUR*uu(i  ,j  ,k) )
      dd(i,j,k) = EIGHTH*ff(i,j,k) - dd(i,j,k)

      i = 1
      j = ny+1
      k = 1
      dd(i,j,k) = ss(15,i,j,k)*(uu(i+1,j-1,k+1) + uu(i+1,j-1,k  ) &
                               +uu(i+1,j  ,k+1) + uu(i  ,j-1,k+1) &
                         - FOUR*uu(i  ,j  ,k) )
      dd(i,j,k) = EIGHTH*ff(i,j,k) - dd(i,j,k)
   
      i = nx+1
      j = 1
      k = 1
      dd(i,j,k) = ss(18,i,j,k)*(uu(i-1,j+1,k+1) + uu(i-1,j+1,k  ) &
                               +uu(i-1,j  ,k+1) + uu(i  ,j+1,k+1) &
                         - FOUR*uu(i  ,j  ,k) )
      dd(i,j,k) = EIGHTH*ff(i,j,k) - dd(i,j,k)
  
      i = nx+1
      j = ny+1
      k = 1
      dd(i,j,k) = ss(13,i,j,k)*(uu(i-1,j-1,k+1) + uu(i-1,j-1,k  ) &
                               +uu(i-1,j  ,k+1) + uu(i  ,j-1,k+1) &
                         - FOUR*uu(i  ,j  ,k) )
      dd(i,j,k) = EIGHTH*ff(i,j,k) - dd(i,j,k)

      i = 1
      j = 1
      k = nz+1
      dd(i,j,k) = ss( 8,i,j,k)*(uu(i+1,j+1,k-1) + uu(i+1,j+1,k  ) &
                               +uu(i+1,j  ,k-1) + uu(i  ,j+1,k-1) &
                         - FOUR*uu(i  ,j  ,k) )
      dd(i,j,k) = EIGHTH*ff(i,j,k) - dd(i,j,k)
  
      i = 1
      j = ny+1
      k = nz+1
      dd(i,j,k) = ss( 3,i,j,k)*(uu(i+1,j-1,k-1) + uu(i+1,j-1,k  ) &
                               +uu(i+1,j  ,k-1) + uu(i  ,j-1,k-1) &
                         - FOUR*uu(i  ,j  ,k) )
      dd(i,j,k) = EIGHTH*ff(i,j,k) - dd(i,j,k)
   
      i = nx+1
      j = 1
      k = nz+1
      dd(i,j,k) = ss( 6,i,j,k)*(uu(i-1,j+1,k-1) + uu(i-1,j+1,k  ) &
                               +uu(i-1,j  ,k-1) + uu(i  ,j+1,k-1) &
                         - FOUR*uu(i  ,j  ,k) )
      dd(i,j,k) = EIGHTH*ff(i,j,k) - dd(i,j,k)
  
      i = nx+1
      j = ny+1
      k = nz+1
      dd(i,j,k) = ss( 1,i,j,k)*(uu(i-1,j-1,k-1) + uu(i-1,j-1,k  ) &
                               +uu(i-1,j  ,k-1) + uu(i  ,j-1,k-1) &
                         - FOUR*uu(i  ,j  ,k) )
      dd(i,j,k) = EIGHTH*ff(i,j,k) - dd(i,j,k)
 
 
!     Lo-x / Lo-y edge
      i = 1
      j = 1
      do k = kstart,kend
         dd(i,j,k) = ss(20,i,j,k)*(uu(i+1,j+1,k+1) + uu(i+1,j+1,k  ) &
                                  +uu(i+1,j  ,k+1) + uu(i  ,j+1,k+1) &
                            - FOUR*uu(i  ,j  ,k) ) &
                    +ss( 8,i,j,k)*(uu(i+1,j+1,k-1) + uu(i+1,j+1,k  ) &
                                  +uu(i+1,j  ,k-1) + uu(i  ,j+1,k-1) &
                            - FOUR*uu(i  ,j  ,k) )
         dd(i,j,k) = FOURTH*ff(i,j,k) - dd(i,j,k)
      end do
      !
      !     Hi-x / Lo-y edge
      !
      i = nx+1
      j = 1
      do k = kstart,kend
         dd(i,j,k) = ss(18,i,j,k)*(uu(i-1,j+1,k+1) + uu(i-1,j+1,k  ) &
                                   +uu(i-1,j  ,k+1) + uu(i  ,j+1,k+1) &
                            - FOUR*uu(i  ,j  ,k) ) &
                    +ss( 6,i,j,k)*(uu(i-1,j+1,k-1) + uu(i-1,j+1,k  ) &
                                  +uu(i-1,j  ,k-1) + uu(i  ,j+1,k-1) &
                            - FOUR*uu(i  ,j  ,k) )
         dd(i,j,k) = FOURTH*ff(i,j,k) - dd(i,j,k)
      end do
      !
      !     Lo-x / Hi-y edge
      !
      i = 1
      j = ny+1
      do k = kstart,kend
         dd(i,j,k) = ss(15,i,j,k)*(uu(i+1,j-1,k+1) + uu(i+1,j-1,k  ) &
                                  +uu(i+1,j  ,k+1) + uu(i  ,j-1,k+1) &
                            - FOUR*uu(i  ,j  ,k) ) &
                    +ss( 3,i,j,k)*(uu(i+1,j-1,k-1) + uu(i+1,j-1,k  ) &
                                  +uu(i+1,j  ,k-1) + uu(i  ,j-1,k-1) &
                            - FOUR*uu(i  ,j  ,k) )
         dd(i,j,k) = FOURTH*ff(i,j,k) - dd(i,j,k)
      end do
      !
      !     Hi-x / Hi-y edge
      !
      i = nx+1
      j = ny+1
      do k = kstart,kend
         dd(i,j,k) = ss(13,i,j,k)*(uu(i-1,j-1,k+1) + uu(i-1,j-1,k  ) &
                                  +uu(i-1,j  ,k+1) + uu(i  ,j-1,k+1) &
                            - FOUR*uu(i  ,j  ,k) ) &
                    +ss( 1,i,j,k)*(uu(i-1,j-1,k-1) + uu(i-1,j-1,k  ) &
                                  +uu(i-1,j  ,k-1) + uu(i  ,j-1,k-1) &
                            - FOUR*uu(i  ,j  ,k) )
         dd(i,j,k) = FOURTH*ff(i,j,k) - dd(i,j,k)
      end do
      !
      !     Lo-x / Lo-z edge
      !
      i = 1
      k = 1
      do j = jstart,jend
         dd(i,j,k) = ss(20,i,j,k)*(uu(i+1,j+1,k+1) + uu(i+1,j+1,k  ) &
                                  +uu(i+1,j  ,k+1) + uu(i  ,j+1,k+1) &
                            - FOUR*uu(i  ,j  ,k) ) &
                    +ss(15,i,j,k)*(uu(i+1,j-1,k+1) + uu(i+1,j-1,k  ) &
                                  +uu(i+1,j  ,k+1) + uu(i  ,j-1,k+1) &
                            - FOUR*uu(i  ,j  ,k) )
         dd(i,j,k) = FOURTH*ff(i,j,k) - dd(i,j,k)
      end do
      !
      !     Hi-x / Lo-z edge
      !
      i = nx+1
      k = 1
      do j = jstart,jend
         dd(i,j,k) = ss(18,i,j,k)*(uu(i-1,j+1,k+1) + uu(i-1,j+1,k  ) &
                                  +uu(i-1,j  ,k+1) + uu(i  ,j+1,k+1) &
                            - FOUR*uu(i  ,j  ,k) ) &
                    +ss(13,i,j,k)*(uu(i-1,j-1,k+1) + uu(i-1,j  ,k+1) &
                                  +uu(i-1,j-1,k  ) + uu(i  ,j-1,k+1) &
                            - FOUR*uu(i  ,j  ,k) )
         dd(i,j,k) = FOURTH*ff(i,j,k) - dd(i,j,k)
      end do
      !
      !     Lo-x / Hi-z edge
      !
      i = 1
      k = nz+1
      do j = jstart,jend
         dd(i,j,k) = ss( 8,i,j,k)*(uu(i+1,j+1,k-1) + uu(i+1,j  ,k-1) &
                                  +uu(i+1,j+1,k  ) + uu(i  ,j+1,k-1) &
                            - FOUR*uu(i  ,j  ,k) ) &
                    +ss( 3,i,j,k)*(uu(i+1,j-1,k-1) + uu(i+1,j-1,k  ) &
                                  +uu(i+1,j  ,k-1) + uu(i  ,j-1,k-1) &
                            - FOUR*uu(i  ,j  ,k) )
         dd(i,j,k) = FOURTH*ff(i,j,k) - dd(i,j,k)
      end do
      !
      !     Hi-x / Hi-z edge
      !
      i = nx+1
      k = nz+1
      do j = jstart,jend
         dd(i,j,k) = ss( 6,i,j,k)*(uu(i-1,j+1,k-1) + uu(i-1,j  ,k-1) &
                                  +uu(i-1,j+1,k  ) + uu(i  ,j+1,k-1) &
                            - FOUR*uu(i  ,j  ,k) ) &
                    +ss( 1,i,j,k)*(uu(i-1,j-1,k-1) + uu(i-1,j-1,k  ) &
                                +uu(i-1,j  ,k-1) + uu(i  ,j-1,k-1) &
                            - FOUR*uu(i  ,j  ,k) )
         dd(i,j,k) = FOURTH*ff(i,j,k) - dd(i,j,k)
      end do
      !
      !     Lo-y / Lo-z edge
      !
      j = 1
      k = 1
      do i = istart,iend
         dd(i,j,k) = ss(20,i,j,k)*(uu(i+1,j+1,k+1) + uu(i+1,j+1,k  ) &
                                  +uu(i+1,j  ,k+1) + uu(i  ,j+1,k+1) &
                            - FOUR*uu(i  ,j  ,k) ) &
                    +ss(18,i,j,k)*(uu(i-1,j+1,k+1) + uu(i-1,j+1,k  ) &
                                  +uu(i  ,j+1,k+1) + uu(i-1,j  ,k+1) &
                            - FOUR*uu(i  ,j  ,k) )
         dd(i,j,k) = FOURTH*ff(i,j,k) - dd(i,j,k)
      end do
      !
      !     Hi-y / Lo-z edge
      !
      j = ny+1
      k = 1
      do i = istart,iend
         dd(i,j,k) = ss(15,i,j,k)*(uu(i+1,j-1,k+1) + uu(i+1,j-1,k  ) &
                                  +uu(i  ,j-1,k+1) + uu(i+1,j  ,k+1) &
                            - FOUR*uu(i  ,j  ,k) ) &
                    +ss(13,i,j,k)*(uu(i-1,j-1,k+1) + uu(i  ,j-1,k+1) &
                                  +uu(i-1,j-1,k  ) + uu(i-1,j  ,k+1) &
                            - FOUR*uu(i  ,j  ,k) )
         dd(i,j,k) = FOURTH*ff(i,j,k) - dd(i,j,k)
      end do
      !
      !     Lo-y / Hi-z edge
      !
      j = 1
      k = nz+1
      do i = istart,iend
         dd(i,j,k) = ss( 8,i,j,k)*(uu(i+1,j+1,k-1) + uu(i+1,j  ,k-1) &
                                  +uu(i+1,j+1,k  ) + uu(i  ,j+1,k-1) &
                          - FOUR*uu(i  ,j  ,k) ) &
                      +ss( 6,i,j,k)*(uu(i-1,j+1,k-1) + uu(i-1,j+1,k  ) &
                                  +uu(i  ,j+1,k-1) + uu(i-1,j  ,k-1) &
                            - FOUR*uu(i  ,j  ,k) )
         dd(i,j,k) = FOURTH*ff(i,j,k) - dd(i,j,k)
      end do
      !
      !     Hi-y / Hi-z edge
      !
      j = ny+1
      k = nz+1
      do i = istart,iend
         dd(i,j,k) = ss( 3,i,j,k)*(uu(i+1,j-1,k-1) + uu(i  ,j-1,k-1) &
                                  +uu(i+1,j-1,k  ) + uu(i+1,j  ,k-1) &
                            - FOUR*uu(i  ,j  ,k) ) &
                    +ss( 1,i,j,k)*(uu(i-1,j-1,k-1) + uu(i-1,j-1,k  ) &
                                  +uu(i-1,j  ,k-1) + uu(i  ,j-1,k-1) &
                            - FOUR*uu(i  ,j  ,k) )
         dd(i,j,k) = FOURTH*ff(i,j,k) - dd(i,j,k)
      end do
      !
      !     Lo-x face
      !
      i = 1
      !$OMP PARALLEL DO PRIVATE(j,k)
      do k = kstart,kend
      do j = jstart,jend
         dd(i,j,k) = ss(20,i,j,k)*(uu(i+1,j+1,k+1) + uu(i+1,j+1,k  ) &
                                  +uu(i+1,j  ,k+1) + uu(i  ,j+1,k+1) &
                            - FOUR*uu(i  ,j  ,k) ) &
                    +ss( 8,i,j,k)*(uu(i+1,j+1,k-1) + uu(i+1,j+1,k  ) &
                                  +uu(i+1,j  ,k-1) + uu(i  ,j+1,k-1) &
                            - FOUR*uu(i  ,j  ,k) ) &
                    +ss(15,i,j,k)*(uu(i+1,j-1,k+1) + uu(i+1,j-1,k  ) &
                                  +uu(i+1,j  ,k+1) + uu(i  ,j-1,k+1) &
                            - FOUR*uu(i  ,j  ,k) ) &
                    +ss( 3,i,j,k)*(uu(i+1,j-1,k-1) + uu(i+1,j-1,k  ) &
                                  +uu(i+1,j  ,k-1) + uu(i  ,j-1,k-1) &
                            - FOUR*uu(i  ,j  ,k) )
         dd(i,j,k) = HALF*ff(i,j,k) - dd(i,j,k)
      end do
      end do
      !$OMP END PARALLEL DO
      !
      !     Hi-x face
      !
      i = nx+1
      !$OMP PARALLEL DO PRIVATE(j,k)
      do k = kstart,kend
      do j = jstart,jend
         dd(i,j,k) = ss(18,i,j,k)*(uu(i-1,j+1,k+1) + uu(i-1,j+1,k  ) &
                                  +uu(i-1,j  ,k+1) + uu(i  ,j+1,k+1) &
                            - FOUR*uu(i  ,j  ,k) ) &
                    +ss(13,i,j,k)*(uu(i-1,j-1,k+1) + uu(i-1,j  ,k+1) &
                                  +uu(i-1,j-1,k  ) + uu(i  ,j-1,k+1) &
                            - FOUR*uu(i  ,j  ,k) ) &
                    +ss( 6,i,j,k)*(uu(i-1,j+1,k-1) + uu(i-1,j  ,k-1) &
                                  +uu(i-1,j+1,k  ) + uu(i  ,j+1,k-1) &
                            - FOUR*uu(i  ,j  ,k) ) &
                    +ss( 1,i,j,k)*(uu(i-1,j-1,k-1) + uu(i-1,j-1,k  ) &
                                  +uu(i-1,j  ,k-1) + uu(i  ,j-1,k-1) &
                          - FOUR*uu(i  ,j  ,k) )
         dd(i,j,k) = HALF*ff(i,j,k) - dd(i,j,k)
      end do
      end do
      !$OMP END PARALLEL DO
      !
      !     Lo-y face
      !
      j = 1
      !$OMP PARALLEL DO PRIVATE(i,k)
      do k = kstart,kend
      do i = istart,iend
         dd(i,j,k) = ss(20,i,j,k)*(uu(i+1,j+1,k+1) + uu(i+1,j+1,k  ) &
                                  +uu(i+1,j  ,k+1) + uu(i  ,j+1,k+1) &
                            - FOUR*uu(i  ,j  ,k) ) &
                    +ss( 8,i,j,k)*(uu(i+1,j+1,k-1) + uu(i+1,j+1,k  ) &
                                  +uu(i+1,j  ,k-1) + uu(i  ,j+1,k-1) &
                            - FOUR*uu(i  ,j  ,k) ) &
                    +ss(18,i,j,k)*(uu(i-1,j+1,k+1) + uu(i-1,j+1,k  ) &
                                  +uu(i-1,j  ,k+1) + uu(i  ,j+1,k+1) &
                            - FOUR*uu(i  ,j  ,k) ) &
                    +ss( 6,i,j,k)*(uu(i-1,j+1,k-1) + uu(i-1,j+1,k  ) &
                                  +uu(i-1,j  ,k-1) + uu(i  ,j+1,k-1) &
                            - FOUR*uu(i  ,j  ,k) )
         dd(i,j,k) = HALF*ff(i,j,k) - dd(i,j,k)
      end do
      end do
      !$OMP END PARALLEL DO
      !
      !     Hi-y face
      !
      j = ny+1
      !$OMP PARALLEL DO PRIVATE(i,k)
      do k = kstart,kend
      do i = istart,iend
         dd(i,j,k) = ss( 3,i,j,k)*(uu(i+1,j-1,k-1) + uu(i  ,j-1,k-1) &
                                  +uu(i+1,j-1,k  ) + uu(i+1,j  ,k-1) &
                            - FOUR*uu(i  ,j  ,k) ) &
                    +ss( 1,i,j,k)*(uu(i-1,j-1,k-1) + uu(i-1,j-1,k  ) &
                                  +uu(i-1,j  ,k-1) + uu(i  ,j-1,k-1) &
                            - FOUR*uu(i  ,j  ,k) ) &
                    +ss(15,i,j,k)*(uu(i+1,j-1,k+1) + uu(i+1,j-1,k  ) &
                                  +uu(i  ,j-1,k+1) + uu(i+1,j  ,k+1) &
                            - FOUR*uu(i  ,j  ,k) ) &
                    +ss(13,i,j,k)*(uu(i-1,j-1,k+1) + uu(i  ,j-1,k+1) &
                                  +uu(i-1,j-1,k  ) + uu(i-1,j  ,k+1) &
                            - FOUR*uu(i  ,j  ,k) )
         dd(i,j,k) = HALF*ff(i,j,k) - dd(i,j,k)
      end do
      end do
      !$OMP END PARALLEL DO
      !
      !     Lo-z face
      !
      k = 1
      !$OMP PARALLEL DO PRIVATE(i,j)
      do j = jstart,jend
      do i = istart,iend
         dd(i,j,k) = ss(15,i,j,k)*(uu(i+1,j-1,k+1) + uu(i+1,j-1,k  ) &
                                  +uu(i  ,j-1,k+1) + uu(i+1,j  ,k+1) &
                            - FOUR*uu(i  ,j  ,k) ) &
                    +ss(13,i,j,k)*(uu(i-1,j-1,k+1) + uu(i  ,j-1,k+1) &
                                  +uu(i-1,j-1,k  ) + uu(i-1,j  ,k+1) &
                            - FOUR*uu(i  ,j  ,k) ) &
                    +ss(20,i,j,k)*(uu(i+1,j+1,k+1) + uu(i+1,j+1,k  ) &
                                  +uu(i+1,j  ,k+1) + uu(i  ,j+1,k+1) &
                            - FOUR*uu(i  ,j  ,k) ) &
                    +ss(18,i,j,k)*(uu(i-1,j+1,k+1) + uu(i-1,j+1,k  ) &
                                  +uu(i  ,j+1,k+1) + uu(i-1,j  ,k+1) &
                            - FOUR*uu(i  ,j  ,k) )
         dd(i,j,k) = HALF*ff(i,j,k) - dd(i,j,k)
      end do
      end do
      !$OMP END PARALLEL DO
      !
      !     Hi-z face
      !
      k = nz+1
      !$OMP PARALLEL DO PRIVATE(i,j)
      do j = jstart,jend
      do i = istart,iend
         dd(i,j,k) = ss( 3,i,j,k)*(uu(i+1,j-1,k-1) + uu(i  ,j-1,k-1) &
                                  +uu(i+1,j-1,k  ) + uu(i+1,j  ,k-1) &
                            - FOUR*uu(i  ,j  ,k) ) &
                    +ss( 1,i,j,k)*(uu(i-1,j-1,k-1) + uu(i-1,j-1,k  ) &
                                  +uu(i-1,j  ,k-1) + uu(i  ,j-1,k-1) &
                            - FOUR*uu(i  ,j  ,k) ) &
                    +ss( 8,i,j,k)*(uu(i+1,j+1,k-1) + uu(i+1,j  ,k-1) &
                                  +uu(i+1,j+1,k  ) + uu(i  ,j+1,k-1) &
                            - FOUR*uu(i  ,j  ,k) ) &
                    +ss( 6,i,j,k)*(uu(i-1,j+1,k-1) + uu(i-1,j+1,k  ) &
                                  +uu(i  ,j+1,k-1) + uu(i-1,j  ,k-1) &
                            - FOUR*uu(i  ,j  ,k) )
         dd(i,j,k) = HALF*ff(i,j,k) - dd(i,j,k)
      end do
      end do
      !$OMP END PARALLEL DO
      !
      !     Interior
      !
      !$OMP PARALLEL DO PRIVATE(i,j,k)
      do k = kstart,kend
      do j = jstart,jend
      do i = istart,iend

        dd(i,j,k) = ss(0,i,j,k)*uu(i,j,k) &
            + ss( 1,i,j,k) * uu(i-1,j-1,k-1) + ss( 2,i,j,k) * uu(i  ,j-1,k-1) &
            + ss( 3,i,j,k) * uu(i+1,j-1,k-1) + ss( 4,i,j,k) * uu(i-1,j  ,k-1) &
            + ss( 5,i,j,k) * uu(i+1,j  ,k-1) + ss( 6,i,j,k) * uu(i-1,j+1,k-1) &
            + ss( 7,i,j,k) * uu(i  ,j+1,k-1) + ss( 8,i,j,k) * uu(i+1,j+1,k-1) &
            + ss( 9,i,j,k) * uu(i-1,j-1,k  ) + ss(10,i,j,k) * uu(i+1,j-1,k  ) &
            + ss(11,i,j,k) * uu(i-1,j+1,k  ) + ss(12,i,j,k) * uu(i+1,j+1,k  ) &
            + ss(13,i,j,k) * uu(i-1,j-1,k+1) + ss(14,i,j,k) * uu(i  ,j-1,k+1) &
            + ss(15,i,j,k) * uu(i+1,j-1,k+1) + ss(16,i,j,k) * uu(i-1,j  ,k+1) &
            + ss(17,i,j,k) * uu(i+1,j  ,k+1) + ss(18,i,j,k) * uu(i-1,j+1,k+1) &
            + ss(19,i,j,k) * uu(i  ,j+1,k+1) + ss(20,i,j,k) * uu(i+1,j+1,k+1)

        if ((size(ss,dim=1) .eq. 27) .and. (.not. uniform_dh) ) then
            !
            ! Add faces (only non-zero for non-uniform dx)
            !
            dd(i,j,k) = dd(i,j,k) + &
                  ss(21,i,j,k) * uu(i-1,j  ,k  ) + ss(22,i,j,k) * uu(i+1,j  ,k  ) &
                  + ss(23,i,j,k) * uu(i  ,j-1,k  ) + ss(24,i,j,k) * uu(i  ,j+1,k  ) &
                  + ss(25,i,j,k) * uu(i  ,j  ,k-1) + ss(26,i,j,k) * uu(i  ,j  ,k+1)
         end if
  
         dd(i,j,k) = ff(i,j,k) - dd(i,j,k)
      end do
      end do
      end do
      !$OMP END PARALLEL DO

    else if (size(ss,dim=1) .eq. 7) then
       !
       !     Corners
       !
      i = 1
      j = 1
      k = 1
      dd(i,j,k) =  ss(1,i,j,k)*(uu(i+1,j,k)-uu(i,j,k)) &
                 + ss(3,i,j,k)*(uu(i,j+1,k)-uu(i,j,k)) &
                 + ss(5,i,j,k)*(uu(i,j,k+1)-uu(i,j,k))
      dd(i,j,k) = EIGHTH*ff(i,j,k) - dd(i,j,k)

      i = 1
      j = ny+1
      k = 1
      dd(i,j,k) =  ss(1,i,j,k)*(uu(i+1,j,k)-uu(i,j,k)) &
                 + ss(4,i,j,k)*(uu(i,j-1,k)-uu(i,j,k)) &
                 + ss(5,i,j,k)*(uu(i,j,k+1)-uu(i,j,k))
      dd(i,j,k) = EIGHTH*ff(i,j,k) - dd(i,j,k)
   
      i = nx+1
      j = 1
      k = 1
      dd(i,j,k) =  ss(2,i,j,k)*(uu(i-1,j,k)-uu(i,j,k)) &
                 + ss(3,i,j,k)*(uu(i,j+1,k)-uu(i,j,k)) &
                 + ss(5,i,j,k)*(uu(i,j,k+1)-uu(i,j,k))
      dd(i,j,k) = EIGHTH*ff(i,j,k) - dd(i,j,k)
  
      i = nx+1
      j = ny+1
      k = 1
      dd(i,j,k) =  ss(2,i,j,k)*(uu(i-1,j,k)-uu(i,j,k)) &
                 + ss(4,i,j,k)*(uu(i,j-1,k)-uu(i,j,k)) &
                 + ss(5,i,j,k)*(uu(i,j,k+1)-uu(i,j,k))
      dd(i,j,k) = EIGHTH*ff(i,j,k) - dd(i,j,k)

      i = 1
      j = 1
      k = nz+1
      dd(i,j,k) =  ss(1,i,j,k)*(uu(i+1,j,k)-uu(i,j,k)) &
                 + ss(3,i,j,k)*(uu(i,j+1,k)-uu(i,j,k)) &
                 + ss(6,i,j,k)*(uu(i,j,k-1)-uu(i,j,k))
      dd(i,j,k) = EIGHTH*ff(i,j,k) - dd(i,j,k)
  
      i = 1
      j = ny+1
      k = nz+1
      dd(i,j,k) =  ss(1,i,j,k)*(uu(i+1,j,k)-uu(i,j,k)) &
                 + ss(4,i,j,k)*(uu(i,j-1,k)-uu(i,j,k)) &
                 + ss(6,i,j,k)*(uu(i,j,k-1)-uu(i,j,k))
      dd(i,j,k) = EIGHTH*ff(i,j,k) - dd(i,j,k)
   
      i = nx+1
      j = 1
      k = nz+1
      dd(i,j,k) =  ss(2,i,j,k)*(uu(i-1,j,k)-uu(i,j,k)) &
                 + ss(3,i,j,k)*(uu(i,j+1,k)-uu(i,j,k)) &
                 + ss(6,i,j,k)*(uu(i,j,k-1)-uu(i,j,k))
      dd(i,j,k) = EIGHTH*ff(i,j,k) - dd(i,j,k)
  
      i = nx+1
      j = ny+1
      k = nz+1
      dd(i,j,k) =  ss(2,i,j,k)*(uu(i-1,j,k)-uu(i,j,k)) &
                 + ss(4,i,j,k)*(uu(i,j-1,k)-uu(i,j,k)) &
                 + ss(6,i,j,k)*(uu(i,j,k-1)-uu(i,j,k))
      dd(i,j,k) = EIGHTH*ff(i,j,k) - dd(i,j,k)
      !
      !     Lo-x / Lo-y edge
      !
      i = 1
      j = 1
      do k = kstart,kend
         dd(i,j,k) =      ss(1,i,j,k)*(uu(i+1,j,k)-uu(i,j,k)) &
                        + ss(3,i,j,k)*(uu(i,j+1,k)-uu(i,j,k)) &
                        + ss(5,i,j,k)*(uu(i,j,k+1)-uu(i,j,k)) &
                        + ss(6,i,j,k)*(uu(i,j,k-1)-uu(i,j,k))
         dd(i,j,k) = FOURTH*ff(i,j,k) - dd(i,j,k)
      end do
      !
      !     Hi-x / Lo-y edge
      !
      i = nx+1
      j = 1
      do k = kstart,kend
         dd(i,j,k) =      ss(2,i,j,k)*(uu(i-1,j,k)-uu(i,j,k)) &
                        + ss(3,i,j,k)*(uu(i,j+1,k)-uu(i,j,k)) &
                        + ss(5,i,j,k)*(uu(i,j,k+1)-uu(i,j,k)) &
                        + ss(6,i,j,k)*(uu(i,j,k-1)-uu(i,j,k))
         dd(i,j,k) = FOURTH*ff(i,j,k) - dd(i,j,k)
      end do
      !
      !     Lo-x / Hi-y edge
      !
      i = 1
      j = ny+1
      do k = kstart,kend
         dd(i,j,k) =      ss(1,i,j,k)*(uu(i+1,j,k)-uu(i,j,k)) &
                        + ss(4,i,j,k)*(uu(i,j-1,k)-uu(i,j,k)) &
                        + ss(5,i,j,k)*(uu(i,j,k+1)-uu(i,j,k)) &
                        + ss(6,i,j,k)*(uu(i,j,k-1)-uu(i,j,k))
         dd(i,j,k) = FOURTH*ff(i,j,k) - dd(i,j,k)
      end do
      !
      !     Hi-x / Hi-y edge
      !
      i = nx+1
      j = ny+1
      do k = kstart,kend
         dd(i,j,k) =      ss(2,i,j,k)*(uu(i-1,j,k)-uu(i,j,k)) &
                        + ss(4,i,j,k)*(uu(i,j-1,k)-uu(i,j,k)) &
                        + ss(5,i,j,k)*(uu(i,j,k+1)-uu(i,j,k)) &
                        + ss(6,i,j,k)*(uu(i,j,k-1)-uu(i,j,k))
         dd(i,j,k) = FOURTH*ff(i,j,k) - dd(i,j,k)
      end do
      !
      !     Lo-x / Lo-z edge
      !
      i = 1
      k = 1
      do j = jstart,jend
         dd(i,j,k) =      ss(1,i,j,k)*(uu(i+1,j,k)-uu(i,j,k)) &
                        + ss(3,i,j,k)*(uu(i,j+1,k)-uu(i,j,k)) &
                        + ss(4,i,j,k)*(uu(i,j-1,k)-uu(i,j,k)) &
                        + ss(5,i,j,k)*(uu(i,j,k+1)-uu(i,j,k)) 
         dd(i,j,k) = FOURTH*ff(i,j,k) - dd(i,j,k)
      end do
      !
      !     Hi-x / Lo-z edge
      !
      i = nx+1
      k = 1
      do j = jstart,jend
         dd(i,j,k) =      ss(2,i,j,k)*(uu(i-1,j,k)-uu(i,j,k)) &
                        + ss(3,i,j,k)*(uu(i,j+1,k)-uu(i,j,k)) &
                        + ss(4,i,j,k)*(uu(i,j-1,k)-uu(i,j,k)) &
                        + ss(5,i,j,k)*(uu(i,j,k+1)-uu(i,j,k)) 
         dd(i,j,k) = FOURTH*ff(i,j,k) - dd(i,j,k)
      end do
      !
      !     Lo-x / Hi-z edge
      !
      i = 1
      k = nz+1
      do j = jstart,jend
         dd(i,j,k) =      ss(1,i,j,k)*(uu(i+1,j,k)-uu(i,j,k)) &
                        + ss(3,i,j,k)*(uu(i,j+1,k)-uu(i,j,k)) &
                        + ss(4,i,j,k)*(uu(i,j-1,k)-uu(i,j,k)) &
                        + ss(6,i,j,k)*(uu(i,j,k-1)-uu(i,j,k)) 
         dd(i,j,k) = FOURTH*ff(i,j,k) - dd(i,j,k)
      end do
      !
      !     Hi-x / Hi-z edge
      !
      i = nx+1
      k = nz+1
      do j = jstart,jend
         dd(i,j,k) =      ss(2,i,j,k)*(uu(i-1,j,k)-uu(i,j,k)) &
                        + ss(3,i,j,k)*(uu(i,j+1,k)-uu(i,j,k)) &
                        + ss(4,i,j,k)*(uu(i,j-1,k)-uu(i,j,k)) &
                        + ss(6,i,j,k)*(uu(i,j,k-1)-uu(i,j,k)) 
         dd(i,j,k) = FOURTH*ff(i,j,k) - dd(i,j,k)
      end do
      !
      !     Lo-y / Lo-z edge
      !
      j = 1
      k = 1
      do i = istart,iend
         dd(i,j,k) =      ss(1,i,j,k)*(uu(i+1,j,k)-uu(i,j,k)) &
                        + ss(2,i,j,k)*(uu(i-1,j,k)-uu(i,j,k)) &
                        + ss(3,i,j,k)*(uu(i,j+1,k)-uu(i,j,k)) &
                        + ss(5,i,j,k)*(uu(i,j,k+1)-uu(i,j,k)) 
         dd(i,j,k) = FOURTH*ff(i,j,k) - dd(i,j,k)
      end do
      !
      !     Hi-y / Lo-z edge
      !
      j = ny+1
      k = 1
      do i = istart,iend
         dd(i,j,k) =      ss(1,i,j,k)*(uu(i+1,j,k)-uu(i,j,k)) &
                        + ss(2,i,j,k)*(uu(i-1,j,k)-uu(i,j,k)) &
                        + ss(4,i,j,k)*(uu(i,j-1,k)-uu(i,j,k)) &
                        + ss(5,i,j,k)*(uu(i,j,k+1)-uu(i,j,k)) 
         dd(i,j,k) = FOURTH*ff(i,j,k) - dd(i,j,k)
      end do
      !
      !     Lo-y / Hi-z edge
      !
      j = 1
      k = nz+1
      do i = istart,iend
         dd(i,j,k) =      ss(1,i,j,k)*(uu(i+1,j,k)-uu(i,j,k)) &
                        + ss(2,i,j,k)*(uu(i-1,j,k)-uu(i,j,k)) &
                        + ss(3,i,j,k)*(uu(i,j+1,k)-uu(i,j,k)) &
                        + ss(6,i,j,k)*(uu(i,j,k-1)-uu(i,j,k)) 
         dd(i,j,k) = FOURTH*ff(i,j,k) - dd(i,j,k)
      end do
      !
      !     Hi-y / Hi-z edge
      !
      j = ny+1
      k = nz+1
      do i = istart,iend
         dd(i,j,k) =      ss(1,i,j,k)*(uu(i+1,j,k)-uu(i,j,k)) &
                        + ss(2,i,j,k)*(uu(i-1,j,k)-uu(i,j,k)) &
                        + ss(4,i,j,k)*(uu(i,j-1,k)-uu(i,j,k)) &
                        + ss(6,i,j,k)*(uu(i,j,k-1)-uu(i,j,k)) 
         dd(i,j,k) = FOURTH*ff(i,j,k) - dd(i,j,k)
      end do
      !
      !     Lo-x face
      !
      i = 1
      !$OMP PARALLEL DO PRIVATE(j,k)
      do k = kstart,kend
      do j = jstart,jend
         dd(i,j,k) =      ss(1,i,j,k)*(uu(i+1,j,k)-uu(i,j,k)) &
                        + ss(3,i,j,k)*(uu(i,j+1,k)-uu(i,j,k)) & 
                        + ss(4,i,j,k)*(uu(i,j-1,k)-uu(i,j,k)) & 
                        + ss(5,i,j,k)*(uu(i,j,k+1)-uu(i,j,k)) & 
                        + ss(6,i,j,k)*(uu(i,j,k-1)-uu(i,j,k))  
         dd(i,j,k) = HALF*ff(i,j,k) - dd(i,j,k)
      end do
      end do
      !$OMP END PARALLEL DO
      !
      !     Hi-x face
      !
      i = nx+1
      !$OMP PARALLEL DO PRIVATE(j,k)
      do k = kstart,kend
      do j = jstart,jend
         dd(i,j,k) =      ss(2,i,j,k)*(uu(i-1,j,k)-uu(i,j,k)) &
                        + ss(3,i,j,k)*(uu(i,j+1,k)-uu(i,j,k)) & 
                        + ss(4,i,j,k)*(uu(i,j-1,k)-uu(i,j,k)) & 
                        + ss(5,i,j,k)*(uu(i,j,k+1)-uu(i,j,k)) & 
                        + ss(6,i,j,k)*(uu(i,j,k-1)-uu(i,j,k))  
         dd(i,j,k) = HALF*ff(i,j,k) - dd(i,j,k)
      end do
      end do
      !$OMP END PARALLEL DO
      !
      !     Lo-y face
      !
      j = 1
      !$OMP PARALLEL DO PRIVATE(i,k)
      do k = kstart,kend
      do i = istart,iend
         dd(i,j,k) =      ss(3,i,j,k)*(uu(i,j+1,k)-uu(i,j,k)) &
                        + ss(1,i,j,k)*(uu(i+1,j,k)-uu(i,j,k)) & 
                        + ss(2,i,j,k)*(uu(i-1,j,k)-uu(i,j,k)) & 
                        + ss(5,i,j,k)*(uu(i,j,k+1)-uu(i,j,k)) & 
                        + ss(6,i,j,k)*(uu(i,j,k-1)-uu(i,j,k))  
         dd(i,j,k) = HALF*ff(i,j,k) - dd(i,j,k)
      end do
      end do
      !$OMP END PARALLEL DO
      !
      !     Hi-y face
      !
      j = ny+1
      !$OMP PARALLEL DO PRIVATE(i,k)
      do k = kstart,kend
      do i = istart,iend
         dd(i,j,k) =      ss(4,i,j,k)*(uu(i,j-1,k)-uu(i,j,k)) &
                        + ss(1,i,j,k)*(uu(i+1,j,k)-uu(i,j,k)) & 
                        + ss(2,i,j,k)*(uu(i-1,j,k)-uu(i,j,k)) & 
                        + ss(5,i,j,k)*(uu(i,j,k+1)-uu(i,j,k)) & 
                        + ss(6,i,j,k)*(uu(i,j,k-1)-uu(i,j,k))  
         dd(i,j,k) = HALF*ff(i,j,k) - dd(i,j,k)
      end do
      end do
      !$OMP END PARALLEL DO
      !
      !     Lo-z face
      !
      k = 1
      !$OMP PARALLEL DO PRIVATE(i,j)
      do j = jstart,jend
      do i = istart,iend
         dd(i,j,k) =      ss(5,i,j,k)*(uu(i,j,k+1)-uu(i,j,k)) &
                        + ss(1,i,j,k)*(uu(i+1,j,k)-uu(i,j,k)) & 
                        + ss(2,i,j,k)*(uu(i-1,j,k)-uu(i,j,k)) & 
                        + ss(3,i,j,k)*(uu(i,j+1,k)-uu(i,j,k)) & 
                        + ss(4,i,j,k)*(uu(i,j-1,k)-uu(i,j,k)) 
         dd(i,j,k) = HALF*ff(i,j,k) - dd(i,j,k)
      end do
      end do
      !$OMP END PARALLEL DO
      !
      !     Hi-z face
      !
      k = nz+1
      !$OMP PARALLEL DO PRIVATE(i,j)
      do j = jstart,jend
      do i = istart,iend
         dd(i,j,k) =      ss(6,i,j,k)*(uu(i,j,k-1)-uu(i,j,k)) &
                        + ss(1,i,j,k)*(uu(i+1,j,k)-uu(i,j,k)) & 
                        + ss(2,i,j,k)*(uu(i-1,j,k)-uu(i,j,k)) & 
                        + ss(3,i,j,k)*(uu(i,j+1,k)-uu(i,j,k)) & 
                        + ss(4,i,j,k)*(uu(i,j-1,k)-uu(i,j,k)) 
         dd(i,j,k) = HALF*ff(i,j,k) - dd(i,j,k)
      end do
      end do
      !$OMP END PARALLEL DO
      !
      !     Interior
      !
      !$OMP PARALLEL DO PRIVATE(i,j,k)
      do k = kstart,kend
      do j = jstart,jend
      do i = istart,iend
          dd(i,j,k) = ss(0,i,j,k)*uu(i,j,k) &
            + ss(1,i,j,k) * uu(i+1,j  ,k  ) + ss(2,i,j,k) * uu(i-1,j  ,k  ) &
            + ss(3,i,j,k) * uu(i  ,j+1,k  ) + ss(4,i,j,k) * uu(i  ,j-1,k  ) &
            + ss(5,i,j,k) * uu(i  ,j  ,k+1) + ss(6,i,j,k) * uu(i  ,j  ,k-1) 
          dd(i,j,k) = ff(i,j,k) - dd(i,j,k)
      end do
      end do
      end do
      !$OMP END PARALLEL DO

    end if

  end subroutine grid_laplace_3d

end module ml_nd_module
