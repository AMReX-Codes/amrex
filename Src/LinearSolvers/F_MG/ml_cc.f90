module ml_cc_module

  use bl_constants_module
  use mg_module
  use ml_layout_module
  use bndry_reg_module
  ! use mg_hypre_module
  
  implicit none

contains

  subroutine ml_cc(mla, mgt, rh, full_soln, fine_mask, &
                   do_diagnostics, need_grad_phi_in, final_resnorm, status)

    use bl_prof_module
    use ml_norm_module          , only : ml_norm_inf
    use ml_cc_restriction_module, only : ml_cc_restriction
    use ml_prolongation_module  , only : ml_cc_prolongation, ml_interp_bcs
    use cc_ml_resid_module      , only : crse_fine_residual_cc
    use backtrace_module        , only : fpe_trap, quiet_nan

    type(ml_layout), intent(in   ) :: mla
    type(mg_tower) , intent(inout) :: mgt(:)
    type( multifab), intent(inout) :: rh(:)
    type( multifab), intent(inout) :: full_soln(:)
    type(lmultifab), intent(in   ) :: fine_mask(:)
    integer        , intent(in   ) :: do_diagnostics

    logical        , intent(in   ), optional :: need_grad_phi_in
    real(dp_t)     , intent(  out), optional :: final_resnorm
    integer        , intent(  out), optional :: status

    integer :: nlevs
    type(multifab), allocatable  ::        uu(:)
    type(multifab), allocatable  ::   uu_hold(:)
    type(multifab), allocatable  ::       res(:)
    type(multifab), allocatable  ::  temp_res(:)

    type(bndry_reg), allocatable :: brs_flx(:)
    type(bndry_reg), allocatable :: brs_bcs(:)

    type(box) :: pd, pdc
    type(layout) :: la, lac
    integer :: n, ng_fill
    integer :: mglev, mglev_crse, iter, iter_solved
    logical :: fine_converged,need_grad_phi
    logical :: using_bnorm

    real(dp_t) :: bnorm, ni_res
    real(dp_t) :: tres, tres0, max_norm
    real(dp_t) :: sum, coeff_sum, coeff_max

    real(dp_t) :: r1,r2,t1(3),t2(3),stime,bottom_solve_time
    logical :: solved

    type(bl_prof_timer), save :: bpt
    
    call bl_proffortfuncstart("ml_cc")
    call bl_proffortfuncstart("ml_cc:0")

    call build(bpt, "ml_cc")

    nlevs             = mla%nlevel
    stime             = parallel_wtime()
    bottom_solve_time = zero

    if ( present(need_grad_phi_in) ) then
       need_grad_phi = need_grad_phi_in 
    else
       need_grad_phi = .false.
    end if

    if (fpe_trap()) then
       do n = 2, nlevs 
          if (mgt(n)%lcross) then
             call multifab_set_corner(full_soln(n), quiet_nan())
          end if
       end do
    end if

    allocate(uu(nlevs), res(nlevs), temp_res(nlevs))
    allocate(uu_hold(2:nlevs-1))
    allocate(brs_flx(2:nlevs))
    allocate(brs_bcs(2:nlevs))

    do n = 2,nlevs-1
       la = mla%la(n)
       call multifab_build( uu_hold(n),la,1,1)
       call multifab_setval(uu_hold(n), ZERO,all=.true.)
    end do

    do n = nlevs, 1, -1

       la = mla%la(n)
       call multifab_build(      uu(n), la, 1, nghost(full_soln(1)))
       call multifab_build(     res(n), la, 1, 0)
       call multifab_build(temp_res(n), la, 1, 0)
       call multifab_setval(      uu(n), ZERO,all=.true.)
       call multifab_setval(     res(n), ZERO,all=.true.)
       call multifab_setval(temp_res(n), ZERO,all=.true.)

       if ( n == 1 ) exit

       ! Build the (coarse resolution) flux registers to be used in computing
       !  the residual at a non-finest AMR level.

       pdc = layout_get_pd(mla%la(n-1))
       lac = mla%la(n-1)
       call bndry_reg_rr_build(brs_bcs(n), la, lac, mla%mba%rr(n-1,:), pdc, width = 2)
       call  flux_reg_build   (brs_flx(n), la, lac, mla%mba%rr(n-1,:), pdc)

    end do

    do n = nlevs,2,-1
       mglev      = mgt(n  )%nlevels
       mglev_crse = mgt(n-1)%nlevels
       call ml_cc_restriction(rh(n-1), rh(n), mla%mba%rr(n-1,:))
    end do

    !
    ! First we must restrict the final solution onto coarser levels, because those coarse
    ! cells may be used to construct slopes in the interpolation of the boundary conditions.
    !
    do n = nlevs,2,-1
       mglev      = mgt(n  )%nlevels
       mglev_crse = mgt(n-1)%nlevels

       pdc = layout_get_pd(mla%la(n-1))

       call ml_cc_restriction(full_soln(n-1), full_soln(n), mla%mba%rr(n-1,:))
    enddo

    !  Now make sure full_soln at fine grid has the correct coarse grid bc's in 
    !  its ghost cells before we evaluate the initial residual  
    do n = 2,nlevs
       ng_fill = nghost(full_soln(n))
       pd = layout_get_pd(mla%la(n))
       call bndry_reg_copy(brs_bcs(n), full_soln(n-1))
       call ml_interp_bcs(full_soln(n), brs_bcs(n)%bmf(1,0), pd, &
            mla%mba%rr(n-1,:), ng_fill, brs_bcs(n)%facemap, brs_bcs(n)%indxmap, &
            brs_bcs(n)%uncovered)
    end do

    !
    ! Elide some reduction.
    !
    t1(1) = max_of_stencil_sum(mgt(1)%ss(1),local=.true.) 
    t1(2) = stencil_norm(mgt(1)%ss(1),local=.true.) 

    call parallel_reduce(t2(1:2), t1(1:2), MPI_MAX)

    coeff_sum = t2(1)
    coeff_max = t2(2)

    !
    ! Test on whether coefficients sum to zero in order to know whether to enforce solvability
    ! Only test on lowest mg level of lowest AMR level -- this should be cheapest
    !
    if ( coeff_sum .lt. (1.0e-12_dp_t * coeff_max) ) then
       mgt(1)%coeffs_sum_to_zero = .true.
       if ( parallel_IOProcessor() .and. (do_diagnostics == 1) ) &
            print *,'Coefficients sum to zero '
    else
       if ( parallel_IOProcessor() .and. (do_diagnostics == 1) ) then
          print *,'Coefficients sum to ', coeff_sum
          print *,' ... coeff_max   is ', coeff_max
          print *,'Not setting singular flag '
       end if
    end if

    ! Make sure to pass this flag through to the bottom_mgt object if there is one. 
    ! Otherwise the BiCG/CG bottom solver will not see it.
    if (associated(mgt(1)%bottom_mgt)) then
       mgt(1)%bottom_mgt%coeffs_sum_to_zero = mgt(1)%coeffs_sum_to_zero
    end if

    ! Enforce solvability if appropriate -- note that we need to make "rh" solvable, not just "res" as before,
    ! because "res" is recomputed as rh - L(full_soln) after each cycle
    if (nlevs .eq. 1 .and. mgt(1)%bottom_singular .and. mgt(1)%coeffs_sum_to_zero .and. &
        mgt(1)%ok_to_fix_singular) then

       sum = multifab_sum(rh(1))  / boxarray_dvolume(get_boxarray(rh(1)))

       ! Subtract "sum" from rh(1) in order to make this solvable
       call sub_sub(rh(1), sum)
          
       if ( parallel_IOProcessor() .and. (do_diagnostics == 1) ) then
          write(unit=*, fmt='("F90mg: Subtracting from rh  ",g15.8)') sum
       end if

    end if

    do n = 1,nlevs,1
       mglev = mgt(n)%nlevels
       call compute_defect(mgt(n)%ss(mglev),res(n),rh(n),full_soln(n), &
                      mgt(n)%mm(mglev), mgt(n)%stencil_type, mgt(n)%lcross)
    end do

    do n = nlevs,2,-1
       mglev      = mgt(n  )%nlevels
       mglev_crse = mgt(n-1)%nlevels

       pdc = layout_get_pd(mla%la(n-1))
       call crse_fine_residual_cc(n,mgt,full_soln,res(n-1),brs_flx(n),pdc, &
            mla%mba%rr(n-1,:), filled=.true.)

       call ml_cc_restriction(res(n-1), res(n), mla%mba%rr(n-1,:))
    enddo

    !
    ! Elide some reduction.
    !
    t1(1) = ml_norm_inf(res,fine_mask,local=.true.)
    t1(2) = ml_norm_inf(rh,fine_mask,local=.true.)

    call parallel_reduce(t2(1:2), t1(1:2), MPI_MAX)

    call bl_proffortfuncstop("ml_cc:0")
    call bl_proffortfuncstart("ml_cc:1")

    tres0     = t2(1)
    bnorm     = t2(2)

    if ( parallel_IOProcessor() .and. mgt(nlevs)%verbose > 0 ) then
       write(unit=*, &
            fmt='("F90mg: Initial rhs                  = ",g15.8)') bnorm
       write(unit=*, &
            fmt='("F90mg: Initial residual (resid0)    = ",g15.8)') tres0
    end if

    ! ************************************************************************
    !  Define norm to be used for convergence testing that is the maximum
    !    of bnorm (norm of rhs) and tres0 (norm of resid0)
    ! ************************************************************************
    if (mgt(1)%always_use_bnorm .or. bnorm .ge. tres0) then
       max_norm = bnorm
       using_bnorm = .true.
    else
       max_norm = tres0
       using_bnorm = .false.
    end if

    fine_converged = .false.

    r1 = parallel_wtime() 

    solved = .false.

    call bl_proffortfuncstart("ml_cc:1.1")

    ! Set flag "optimistically", 0 indicates no problems (1: smoother failed, <0: too many mlmg iterations)
    if ( present(status) ) status = 0

    if ( ml_converged(res, fine_mask, max_norm, mgt(nlevs)%eps, mgt(nlevs)%abs_eps, ni_res, mgt(nlevs)%verbose) ) then

       solved = .true.

       if ( parallel_IOProcessor() .and. mgt(nlevs)%verbose > 0 ) &
            write(unit=*, fmt='("F90mg: No iterations needed ")')

       n_mg_iters = 0

       !   else if (mgt%use_hypre .eq. 1) then

       !      if (nlevs .gt. 1) then
       !         call bl_error("ml_cc: can't use HYPRE with nlevs > 1")
       !      else
       !         call cc_hypre_solve(mgt, res, full_soln, eps)
       !      end if

    else 

       ! Only print this once per solve.
       do n = 2,nlevs
           if (parallel_IOProcessor() .and. &
                mla%mba%rr(n-1,1) /= 2 .and. mgt(n-1)%use_lininterp) then
              call bl_warn('ml_cc: linear prolongation not supported since ir /= 2')
           end if
       end do

       do iter = 1, mgt(nlevs)%max_iter

          iter_solved = iter 

          if ( fine_converged ) then
             if ( ml_converged(res, fine_mask, max_norm, mgt(nlevs)%eps, &
                  mgt(nlevs)%abs_eps, ni_res, mgt(nlevs)%verbose) ) then

                ! Subtract one from "iter" here because we have already converged
                iter_solved = iter-1
                solved = .true.
                exit

             end if
          end if

          ! Set: uu = 0
          do n = 1,nlevs
             call setval(uu(n), ZERO, all=.true.)
          end do

          ! Set: uu_hold = 0
          do n = 2,nlevs-1
             call setval(uu_hold(n), ZERO, all=.true.)
          end do

          call bl_proffortfuncstart("ml_cc:DownV")
          !   Down the V-cycle
          do n = nlevs,1,-1

             if ( do_diagnostics == 1 .and.  parallel_ioprocessor() ) then
                print *,' '
                print *,'AT AMR LEVEL ',n
             end if

             mglev = mgt(n)%nlevels

             ! Enforce solvability if appropriate
             if (n .eq. 1 .and. mgt(1)%bottom_singular .and. mgt(1)%coeffs_sum_to_zero) then

                sum = multifab_sum(res(1))  / boxarray_dvolume(get_boxarray(res(1)))

                ! Subtract "sum" from res(1) in order to make this solvable
                call sub_sub(res(1), sum)

                if ( parallel_IOProcessor() .and. (do_diagnostics == 1) ) then
                   write(unit=*, fmt='("F90mg: Subtracting from res ",g15.8)') sum
                end if

             end if

             if ( do_diagnostics == 1 ) then
                tres = norm_inf(res(n))
                if ( parallel_ioprocessor() ) then
                   print *,'DWN: Norm before relax ',tres
                end if
             end if

             ! Relax ...
             call bl_proffortfuncstart("ml_cc:Relax")
             if (n > 1) then
                call mini_cycle(mgt(n), mglev, &
                     mgt(n)%ss(mglev), uu(n), res(n), &
                     mgt(n)%mm(mglev), mgt(n)%nu1, mgt(n)%nu2)
             else
                if (mgt(n)%cycle_type == MG_FVCycle) then
                    if (iter == 1) then
                        mgt(1)%use_lininterp = .true.
                        call mg_tower_cycle(mgt(n), MG_FCycle, mglev, &
                             mgt(n)%ss(mglev), uu(n), res(n), &
                             mgt(n)%mm(mglev), mgt(n)%nu1, mgt(n)%nu2, &
                             bottom_solve_time = bottom_solve_time)
                    else
                        mgt(1)%use_lininterp = .false.
                        call mg_tower_cycle(mgt(n), MG_VCycle, mglev, &
                             mgt(n)%ss(mglev), uu(n), res(n), &
                             mgt(n)%mm(mglev), mgt(n)%nu1, mgt(n)%nu2, &
                             bottom_solve_time = bottom_solve_time)
                    end if
                else if (mgt(n)%cycle_type == MG_VCycle) then
                    mgt(1)%use_lininterp = .false.
                    call mg_tower_cycle(mgt(n), mgt(n)%cycle_type, mglev, &
                         mgt(n)%ss(mglev), uu(n), res(n), &
                         mgt(n)%mm(mglev), mgt(n)%nu1, mgt(n)%nu2, &
                         bottom_solve_time = bottom_solve_time)
                else if (mgt(n)%cycle_type == MG_FCycle) then
                    mgt(1)%use_lininterp = .true.
                    call mg_tower_cycle(mgt(n), mgt(n)%cycle_type, mglev, &
                         mgt(n)%ss(mglev), uu(n), res(n), &
                         mgt(n)%mm(mglev), mgt(n)%nu1, mgt(n)%nu2, &
                         bottom_solve_time = bottom_solve_time)
                end if
             end if
             call bl_proffortfuncstop("ml_cc:Relax")

             ! Add: Soln += uu
             call plus_plus(full_soln(n), uu(n))

             if (n > 1) then

                mglev_crse = mgt(n-1)%nlevels

                ! Compute COARSE Res = Rh - Lap(Soln)
                call compute_defect(mgt(n-1)%ss(mglev_crse),res(n-1), &
                               rh(n-1),full_soln(n-1),mgt(n-1)%mm(mglev_crse), &
                               mgt(n-1)%stencil_type, mgt(n-1)%lcross)

                ! Compute FINE Res = Res - Lap(uu)
                call compute_defect(mgt(n)%ss(mglev),temp_res(n), &
                               res(n),uu(n),mgt(n)%mm(mglev), &
                               mgt(n)%stencil_type, mgt(n)%lcross)
                call multifab_copy(res(n), temp_res(n), ng = nghost(res(n)))

                if (do_diagnostics == 1 ) then
                   tres = norm_inf(res(n))
                   if ( parallel_ioprocessor()) then
                      print *,'DWN: Norm  after relax ',tres
                   end if
                end if

                ! Compute CRSE-FINE Res = Res - Crse Flux(soln) + Fine Flux(soln)
                pdc = layout_get_pd(mla%la(n-1))
                call crse_fine_residual_cc(n,mgt,full_soln,res(n-1),brs_flx(n), &
                     pdc,mla%mba%rr(n-1,:))

                ! Restrict FINE Res to COARSE Res (important to do this last
                !     so we overwrite anything extra which may have been defined
                !     above near fine-fine interfaces)
                call ml_cc_restriction(res(n-1), res(n), mla%mba%rr(n-1,:))

                ! Copy u_hold = uu
                if (n < nlevs) call multifab_copy(uu_hold(n), uu(n), ng = nghost(uu(n)))

                ! Set: uu = 0
                call setval(uu(n), ZERO, all=.true.)

             else

                if (do_diagnostics == 1 ) then
                   call compute_defect(mgt(n)%ss(mglev), temp_res(n), res(n), uu(n), &
                        mgt(n)%mm(mglev), mgt(n)%stencil_type, mgt(n)%lcross)
                   tres = norm_inf(temp_res(n))
                   if ( parallel_ioprocessor()) then
                      print *,'DWN: RES AFTER  GSRB AT LEVEL ',n, tres
                   end if
                else
                   ! Seem to need this in periodic case to get right answer...
                   call multifab_fill_boundary(uu(n), cross = mgt(n)%lcross)
                end if

             end if

          end do
          call bl_proffortfuncstop("ml_cc:DownV")

          call bl_proffortfuncstart("ml_cc:UpV")
          !   Back up the V-cycle
          do n = 2, nlevs

             if ( do_diagnostics == 1 .and.  parallel_ioprocessor() ) &
                  print *,'AT AMR LEVEL ',n

             pd = layout_get_pd(mla%la(n))
             mglev = mgt(n)%nlevels

             ! Interpolate uu from coarser level
             call ml_cc_prolongation(uu(n), uu(n-1), mla%mba%rr(n-1,:), &
                                     mgt(n-1)%use_lininterp, mgt(n-1)%ptype)

             ! Add: soln(n) += uu
             call plus_plus(full_soln(n), uu(n), nghost(uu(n)))

             ! Add: uu_hold += uu
             if (n < nlevs) call plus_plus(uu_hold(n), uu(n), nghost(uu(n)))

             ! Interpolate uu to supply boundary conditions for new 
             ! residual calculation
             call bndry_reg_copy(brs_bcs(n), uu(n-1))
             ng_fill = nghost(uu(n))
             call ml_interp_bcs(uu(n), brs_bcs(n)%bmf(1,0), pd, &
                  mla%mba%rr(n-1,:), ng_fill, brs_bcs(n)%facemap, brs_bcs(n)%indxmap, &
                  brs_bcs(n)%uncovered)
             call multifab_fill_boundary(uu(n))

             ! Compute Res = Res - Lap(uu)
             call compute_defect(mgt(n)%ss(mglev), temp_res(n), res(n), uu(n), &
                            mgt(n)%mm(mglev), mgt(n)%stencil_type, mgt(n)%lcross, filled=.true.)
             call multifab_copy(res(n), temp_res(n), ng = nghost(res(n)))

             if (do_diagnostics == 1 ) then
                tres = norm_inf(temp_res(n))
                if ( parallel_ioprocessor() ) then
                   print *,'UP : RES BEFORE GSRB AT LEVEL ',n, tres
                end if
             end if

             ! Set: uu = 0
             call setval(uu(n), ZERO, all=.true.)

             ! Relax ...
             call mini_cycle(mgt(n), mglev, mgt(n)%ss(mglev), &
                  uu(n), res(n), mgt(n)%mm(mglev), mgt(n)%nu1, mgt(n)%nu2)

             ! Compute Res = Res - Lap(uu)

             if (do_diagnostics == 1 ) then
                call compute_defect(mgt(n)%ss(mglev), temp_res(n), res(n), uu(n), &
                               mgt(n)%mm(mglev), mgt(n)%stencil_type, mgt(n)%lcross)
                call multifab_copy(res(n), temp_res(n), ng = nghost(res(n)))
                tres = norm_inf(res(n))
                if ( parallel_ioprocessor() ) then
                   print *,'UP : RES AFTER  GSRB AT LEVEL ',n, tres
                   if (n == nlevs) print *,' '
                end if
             end if

             ! Add: soln += uu
             call plus_plus(full_soln(n), uu(n), nghost(uu(n)))

             ! Add: uu += uu_hold so that it will be interpolated too
             if (n < nlevs) call plus_plus(uu(n), uu_hold(n), nghost(uu_hold(n)))

             ! Only do this as long as tangential interp looks under fine grids
             mglev_crse = mgt(n-1)%nlevels
             call ml_cc_restriction(full_soln(n-1), full_soln(n), mla%mba%rr(n-1,:))

          end do
          call bl_proffortfuncstop("ml_cc:UpV")

          !    Average the solution to the coarser grids.
          do n = nlevs,2,-1
             mglev      = mgt(n)%nlevels
             mglev_crse = mgt(n-1)%nlevels
             call ml_cc_restriction(full_soln(n-1), full_soln(n), mla%mba%rr(n-1,:))
          end do

          ! Interpolate soln to supply boundary conditions 
          do n = 2,nlevs
             ng_fill = nghost(full_soln(n))
             pd = layout_get_pd(mla%la(n))
             call bndry_reg_copy(brs_bcs(n), full_soln(n-1))
             call ml_interp_bcs(full_soln(n), brs_bcs(n)%bmf(1,0), pd, &
                  mla%mba%rr(n-1,:), ng_fill, brs_bcs(n)%facemap, brs_bcs(n)%indxmap, &
                  brs_bcs(n)%uncovered)
          end do

          !    Optimization so don't have to do multilevel convergence test 
          !    each time

          !    Compute the residual on just the finest level
          n = nlevs
          mglev = mgt(n)%nlevels
          call compute_defect(mgt(n)%ss(mglev),res(n),rh(n),full_soln(n), &
                         mgt(n)%mm(mglev), mgt(n)%stencil_type, mgt(n)%lcross)

          if ( ml_fine_converged(res, max_norm, mgt(nlevs)%eps, mgt(nlevs)%abs_eps) ) then

             fine_converged = .true.

             !      Compute the residual on every level
             do n = 1,nlevs-1
                mglev = mgt(n)%nlevels
                call compute_defect(mgt(n)%ss(mglev),res(n),rh(n),full_soln(n), &
                               mgt(n)%mm(mglev), mgt(n)%stencil_type, mgt(n)%lcross)
             end do

             !      Compute the coarse-fine residual 
             do n = nlevs,2,-1
                pdc = layout_get_pd(mla%la(n-1))
                call crse_fine_residual_cc(n,mgt,full_soln,res(n-1),brs_flx(n),pdc, &
                     mla%mba%rr(n-1,:), filled=.true.)
             end do

             !      Average the fine residual onto the coarser level
             do n = nlevs,2,-1
                mglev      = mgt(n  )%nlevels
                mglev_crse = mgt(n-1)%nlevels
                call ml_cc_restriction(res(n-1), res(n), mla%mba%rr(n-1,:))
             end do

             if ( mgt(nlevs)%verbose > 1 ) then
                do n = 1,nlevs
                   tres = norm_inf(res(n))
                   if ( parallel_ioprocessor() ) then
                      if (using_bnorm) then
                          write(unit=*, fmt='("F90mg: Iteration   ",i3," Lev ",i1," resid/bnorm = ",g15.8)') &
                               iter,n,tres/max_norm
                      else
                          write(unit=*, fmt='("F90mg: Iteration   ",i3," Lev ",i1," resid/resid0 = ",g15.8)') &
                               iter,n,tres/max_norm
                      end if
                   end if
                end do
             end if

          else 

             fine_converged = .false.

             ! Only compute this once if we're going to use it
             if (mgt(nlevs)%max_L0_growth > 1 .or. mgt(nlevs)%verbose > 1 ) &
                tres = norm_inf(res(nlevs))

             ! Test on growth only if max_L0_growth > 1
             if ( (mgt(nlevs)%max_L0_growth > 1) .and. (tres/tres0 > mgt(nlevs)%max_L0_growth) ) then
                if ( mgt(nlevs)%verbose > 1  .and. parallel_IOProcessor() ) &
                     write(unit=*, fmt='("F90mg: Iteration blowing up, bailing...")')

                solved = .false.

                ! set flag to indicate that smoother failed
                if ( present(status) ) status = 1
                exit

             end if

             if ( mgt(nlevs)%verbose > 1 .and. parallel_IOProcessor() ) then
                if (using_bnorm) then
                    write(unit=*, fmt='("F90mg: Iteration   ",i3," Fine  resid/bnorm  = ",g15.8)') iter,tres/max_norm
                else
                    write(unit=*, fmt='("F90mg: Iteration   ",i3," Fine  resid/resid0 = ",g15.8)') iter,tres/max_norm
                end if
             end if

          end if

       enddo

       ! if status==0, but not solved then we ran out of iterations
       ! Set status to (neg)num_iters
       if ( present(status) ) then 
          if (status .eq. 0  .and.  .not. solved) then 
             status = -iter
          end if
       end if

       ! ****************************************************************************
       if ( solved ) then
          if ( mgt(nlevs)%verbose > 0 ) then
             !
             ! Consolidate these reductions.
             !
             t1(1) = ml_norm_inf(res,fine_mask,local=.true.) 
             t1(2) = (parallel_wtime() - stime)
             t1(3) = bottom_solve_time

             call parallel_reduce(t2, t1, MPI_MAX, proc = parallel_IOProcessorNode())

             n_mg_iters = iter_solved

             if ( parallel_IOProcessor() ) then

                if (using_bnorm) then
                   write(unit=*, fmt='("F90mg: Final Iter. ",i3," resid/bnorm  = ",g15.8)') iter_solved,t2(1)/max_norm
                else
                   write(unit=*, fmt='("F90mg: Final Iter. ",i3," resid/resid0 = ",g15.8)') iter_solved,t2(1)/max_norm
                end if

                write(unit=*, fmt='("F90mg: Solve time: ",g13.6, " Bottom Solve time: ", g13.6)') t2(2), t2(3)
                write(unit=*, fmt='("")')
             end if
          end if
       else
          if (.not. present(status) .and. mgt(nlevs)%abort_on_max_iter) then
             call bl_error("Multigrid Solve: failed to converge in max_iter iterations")
          end if
       end if

    end if

    call bl_proffortfuncstop("ml_cc:1.1")

    do n = 2,nlevs-1
       call multifab_destroy(uu_hold(n))
    end do

    if (solved .and. need_grad_phi) then

       !   First we must restrict the final solution onto coarser levels, because those coarse
       !   cells may be used to construct slopes in the interpolation of the boundary conditions
       do n = nlevs,2,-1
          mglev      = mgt(n  )%nlevels
          mglev_crse = mgt(n-1)%nlevels

          pdc = layout_get_pd(mla%la(n-1))

          call ml_cc_restriction(full_soln(n-1), full_soln(n), mla%mba%rr(n-1,:))
       enddo

       !   Interpolate boundary conditions of soln in order to get correct grad(phi) at
       !   crse-fine boundaries
       do n = 2,nlevs
          ng_fill = nghost(full_soln(n))
          pd = layout_get_pd(mla%la(n))
          call bndry_reg_copy(brs_bcs(n), full_soln(n-1))
          call ml_interp_bcs(full_soln(n), brs_bcs(n)%bmf(1,0), pd, &
               mla%mba%rr(n-1,:), ng_fill, brs_bcs(n)%facemap, brs_bcs(n)%indxmap, &
               brs_bcs(n)%uncovered)
       end do

    end if

    if (solved) then
       do n = 1,nlevs
          call multifab_fill_boundary(full_soln(n))
       end do
    end if

    do n = nlevs, 1, -1
       call multifab_destroy(      uu(n))
       call multifab_destroy(     res(n))
       call multifab_destroy(temp_res(n))
       if ( n == 1 ) exit
       call bndry_reg_destroy(brs_flx(n))
       call bndry_reg_destroy(brs_bcs(n))
    end do

    call destroy(bpt)

    if (solved) then
       if ( present(final_resnorm) ) final_resnorm = ni_res
       r2 = (parallel_wtime() - r1)
       call parallel_reduce(r1, r2, MPI_MAX, proc = parallel_IOProcessorNode())
       if ( parallel_IOProcessor() .and. mgt(nlevs)%verbose > 0 ) &
            print*, 'Solve Time = ', r1
    end if

    call bl_proffortfuncstop("ml_cc:1")
    call bl_proffortfuncstop("ml_cc")

  contains

    function ml_fine_converged(res, max_norm, rel_eps, abs_eps) result(r)
      logical                    :: r
      type(multifab), intent(in) :: res(:)
      real(dp_t),     intent(in) :: rel_eps, abs_eps, max_norm
      real(dp_t)                 :: ni_res
      integer                    :: nlevs
      nlevs = size(res)
      ni_res = norm_inf(res(nlevs))
      r = ( ni_res <= rel_eps*(max_norm) .or. ni_res <= abs_eps)
    end function ml_fine_converged

    function ml_converged(res, mask, max_norm, rel_eps, abs_eps, ni_res, verbose) result(r)

      use ml_norm_module, only : ml_norm_inf

      logical                      :: r
      integer                      :: verbose
      type(multifab),  intent(in ) :: res(:)
      type(lmultifab), intent(in ) :: mask(:)
      real(dp_t),      intent(in ) :: rel_eps, abs_eps, max_norm
      real(dp_t),      intent(out) :: ni_res
      real(dp_t) :: l_ni_res
      l_ni_res = ml_norm_inf(res, mask, local=.true.)
      call parallel_reduce(ni_res, l_ni_res, MPI_MAX) 
      r = ( ni_res <= rel_eps*(max_norm) .or. ni_res <= abs_eps )
      if ( r .and. parallel_IOProcessor() .and. verbose > 1) then
         if ( ni_res <= rel_eps*max_norm ) then
            print *,'Converged res < rel_eps*max_norm ', ni_res, rel_eps
         else if ( ni_res <= abs_eps ) then
            print *,'Converged res < abs_eps ', ni_res, abs_eps
         end if
      end if
    end function ml_converged

  end subroutine ml_cc

end module ml_cc_module
