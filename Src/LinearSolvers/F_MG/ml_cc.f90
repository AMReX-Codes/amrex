module ml_cc_module

  use bl_constants_module
  use mg_module
  use ml_layout_module
  use bndry_reg_module

  implicit none

contains

  subroutine ml_cc(mla, mgt, rh, full_soln, fine_mask, ref_ratio, &
                   do_diagnostics, eps, abs_eps_in, need_grad_phi_in, final_resnorm, &
                   bottom_mgt)

    use bl_prof_module
    use ml_util_module, only: ml_norm_inf
    use ml_restriction_module, only: ml_restriction
    use ml_prolongation_module, only: ml_prolongation, ml_interp_bcs

    type(ml_layout), intent(in   ) :: mla
    type(mg_tower) , intent(inout) :: mgt(:)
    type( multifab), intent(inout) :: rh(:)
    type( multifab), intent(inout) :: full_soln(:)
    type(lmultifab), intent(in   ) :: fine_mask(:)
    integer        , intent(in   ) :: ref_ratio(:,:)
    integer        , intent(in   ) :: do_diagnostics
    real(dp_t)     , intent(in   ) :: eps

    type(mg_tower) , intent(inout), optional :: bottom_mgt

    real(dp_t)     , intent(in   ), optional :: abs_eps_in
    logical        , intent(in   ), optional :: need_grad_phi_in
    real(dp_t)     , intent(  out), optional :: final_resnorm

    integer :: nlevs
    type(multifab), allocatable  ::      soln(:)
    type(multifab), allocatable  ::        uu(:)
    type(multifab), allocatable  ::   uu_hold(:)
    type(multifab), allocatable  ::       res(:)
    type(multifab), allocatable  ::  temp_res(:)

    type(bndry_reg), allocatable :: brs_flx(:)
    type(bndry_reg), allocatable :: brs_bcs(:)

    type(box) :: pd, pdc
    type(layout) :: la, lac
    integer :: i, n, dm
    integer :: mglev, mglev_crse, iter
    logical :: fine_converged,need_grad_phi,lcross

    real(dp_t) :: Anorm, bnorm, abs_eps, ni_res
    real(dp_t) :: tres, tres0

    type(bl_prof_timer), save :: bpt

    call build(bpt, "ml_cc")

    if ( present(abs_eps_in) ) then
       abs_eps = abs_eps_in 
    else
       abs_eps = 0.d0
    end if

    if ( present(need_grad_phi_in) ) then
       need_grad_phi = need_grad_phi_in 
    else
       need_grad_phi = .false.
    end if

    dm = rh(1)%dim

    nlevs = mla%nlevel

    allocate(soln(nlevs), uu(nlevs), res(nlevs), temp_res(nlevs))
    allocate(uu_hold(2:nlevs-1))
    allocate(brs_flx(2:nlevs))
    allocate(brs_bcs(2:nlevs))

    do n = 2,nlevs-1
       la = mla%la(n)
       call build(uu_hold(n),la,1,1)
       call setval( uu_hold(n), ZERO,all=.true.)
    end do

    do n = nlevs, 1, -1

       la = mla%la(n)
       call build(    soln(n), la, 1, full_soln(1)%ng)
       call build(      uu(n), la, 1, full_soln(1)%ng)
       call build(     res(n), la, 1, 0)
       call build(temp_res(n), la, 1, 0)
       call setval(    soln(n), ZERO,all=.true.)
       call setval(      uu(n), ZERO,all=.true.)
       call setval(     res(n), ZERO,all=.true.)
       call setval(temp_res(n), ZERO,all=.true.)

       if ( n == 1 ) exit

       ! Build the (coarse resolution) flux registers to be used in computing
       !  the residual at a non-finest AMR level.

       pdc = layout_get_pd(mla%la(n-1))
       lac = mla%la(n-1)
       call bndry_reg_rr_build(brs_flx(n), la, lac, ref_ratio(n-1,:), pdc, width = 0)
       call bndry_reg_rr_build(brs_bcs(n), la, lac, ref_ratio(n-1,:), pdc, width = 2, other = .false.)

    end do

    do n = nlevs,2,-1
       mglev      = mgt(n  )%nlevels
       mglev_crse = mgt(n-1)%nlevels
       call ml_restriction(rh(n-1), rh(n), mgt(n)%mm(mglev),&
            mgt(n-1)%mm(mglev_crse), mgt(n)%face_type, ref_ratio(n-1,:))
    end do
    bnorm = ml_norm_inf(rh,fine_mask)

    lcross = ((ncomp(mgt(nlevs)%ss(mgt(nlevs)%nlevels)) == 5) .or. &
         (ncomp(mgt(nlevs)%ss(mgt(nlevs)%nlevels)) == 7))

    Anorm = stencil_norm(mgt(nlevs)%ss(mgt(nlevs)%nlevels))
    do n = 1, nlevs-1
       Anorm = max(stencil_norm(mgt(n)%ss(mgt(n)%nlevels),fine_mask(n)),Anorm)
    end do

    !  Make sure full_soln at fine grid has the correct coarse grid bc's in 
    !  its ghost cells before we evaluate the initial residual  
    do n = 2,nlevs
       pd = layout_get_pd(mla%la(n))
       call bndry_reg_copy(brs_bcs(n), full_soln(n-1))
       do i = 1, dm
          call ml_interp_bcs(full_soln(n), brs_bcs(n)%bmf(i,0), pd, &
                             ref_ratio(n-1,:), -i)
          call ml_interp_bcs(full_soln(n), brs_bcs(n)%bmf(i,1), pd, &
                             ref_ratio(n-1,:), +i)
       end do
       call multifab_fill_boundary(full_soln(n))
    end do

    !   Make sure all periodic and internal boundaries are filled as well
    do n = 1,nlevs   
       call multifab_fill_boundary(full_soln(n))
    end do

    do n = 1,nlevs,1
       mglev = mgt(n)%nlevels
       call mg_defect(mgt(n)%ss(mglev),res(n),rh(n),full_soln(n), &
                      mgt(n)%mm(mglev))
    end do

    do n = nlevs,2,-1
       mglev      = mgt(n  )%nlevels
       mglev_crse = mgt(n-1)%nlevels

       pdc = layout_get_pd(mla%la(n-1))
       call crse_fine_residual_cc(n,mgt,full_soln,res(n-1),brs_flx(n),pdc, &
                                  ref_ratio(n-1,:))

       call ml_restriction(res(n-1), res(n), mgt(n)%mm(mglev),&
            mgt(n-1)%mm(mglev_crse), mgt(n)%face_type, ref_ratio(n-1,:))
    enddo

    do n = 1,nlevs
       call multifab_copy(rh(n),res(n),ng=rh(n)%ng)
    end do

    tres0 = ml_norm_inf(rh,fine_mask)
    if ( parallel_IOProcessor() .and. mgt(nlevs)%verbose > 0 ) then
       write(unit=*, &
             fmt='("F90mg: Initial rhs                  = ",g15.8)') bnorm
       write(unit=*, &
             fmt='("F90mg: Initial error (error0)       = ",g15.8)') tres0
    end if

    ! ************************************************************************

    fine_converged = .false.

    if ( ml_converged(res, soln, fine_mask, bnorm, Anorm, eps, abs_eps, ni_res, mgt(nlevs)%verbose) ) then
       if ( parallel_IOProcessor() .and. mgt(nlevs)%verbose > 0 ) &
            write(unit=*, fmt='("F90mg: No iterations needed ")')

    else

     do iter = 1, mgt(nlevs)%max_iter

       if ( fine_converged ) then
          if ( ml_converged(res, soln, fine_mask, bnorm, Anorm, eps, abs_eps, ni_res, mgt(nlevs)%verbose) ) exit
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
          if (iter < mgt(nlevs)%max_iter) then
             if (n > 1) then
                call mini_cycle(mgt(n), mgt(n)%cycle, mglev, &
                                mgt(n)%ss(mglev), uu(n), res(n), &
                                mgt(n)%mm(mglev), mgt(n)%nu1, mgt(n)%nu2, &
                                mgt(n)%gamma)
             else 
                if (present(bottom_mgt)) then
                   call mg_tower_cycle(mgt(n), mgt(n)%cycle, mglev, &
                                       mgt(n)%ss(mglev), uu(n), res(n), &
                                       mgt(n)%mm(mglev), mgt(n)%nu1, mgt(n)%nu2, &
                                       mgt(n)%gamma, bottom_mgt=bottom_mgt)
                else 
                   call mg_tower_cycle(mgt(n), mgt(n)%cycle, mglev, &
                                       mgt(n)%ss(mglev), uu(n), res(n), &
                                       mgt(n)%mm(mglev), mgt(n)%nu1, mgt(n)%nu2, &
                                       mgt(n)%gamma)
                end if
             end if
          end if

          ! Add: Soln += uu
          call saxpy(soln(n), ONE, uu(n))

          if (n > 1) then

             mglev_crse = mgt(n-1)%nlevels

             ! Compute COARSE Res = Rh - Lap(Soln)
             call mg_defect(mgt(n-1)%ss(mglev_crse),res(n-1), &
                  rh(n-1),soln(n-1),mgt(n-1)%mm(mglev_crse))

             ! Compute FINE Res = Res - Lap(uu)
             call mg_defect(mgt(n)%ss(mglev),temp_res(n), &
                  res(n),uu(n),mgt(n)%mm(mglev))
             call multifab_copy(res(n), temp_res(n), ng=res(n)%ng)

             if (do_diagnostics == 1 ) then
                tres = norm_inf(res(n))
                if ( parallel_ioprocessor()) then
                   print *,'DWN: RES AFTER  GSRB AT LEVEL ',n, tres
                end if
             end if

             ! Compute CRSE-FINE Res = Res - Crse Flux(soln) + Fine Flux(soln)
             pdc = layout_get_pd(mla%la(n-1))
             call crse_fine_residual_cc(n,mgt,soln,res(n-1),brs_flx(n), &
                                        pdc,ref_ratio(n-1,:))

             ! Restrict FINE Res to COARSE Res (important to do this last
             !     so we overwrite anything extra which may have been defined
             !     above near fine-fine interfaces)
             call ml_restriction(res(n-1), res(n), mgt(n)%mm(mglev), &
                  mgt(n-1)%mm(mglev_crse), mgt(n)%face_type, ref_ratio(n-1,:))

             ! Copy u_hold = uu
             if (n < nlevs) call multifab_copy(uu_hold(n), uu(n), ng=uu(n)%ng)

             ! Set: uu = 0
             call setval(uu(n), ZERO, all=.true.)

          else

             if (do_diagnostics == 1 ) then
                call mg_defect(mgt(n)%ss(mglev), temp_res(n), res(n), uu(n), &
                               mgt(n)%mm(mglev))
                tres = norm_inf(temp_res(n))
                if ( parallel_ioprocessor()) then
                   print *,'DWN: RES AFTER  GSRB AT LEVEL ',n, tres
                end if
             else
                ! Seem to need this in periodic case to get right answer...
                call multifab_fill_boundary(uu(n), cross = lcross)
             end if

          end if

       end do

       !   Back up the V-cycle
       do n = 2, nlevs

          pd = layout_get_pd(mla%la(n))
          mglev = mgt(n)%nlevels

          ! Interpolate uu from coarser level
          call ml_prolongation(uu(n), uu(n-1), pd, ref_ratio(n-1,:))

          ! Add: soln += uu
          call saxpy(soln(n), ONE, uu(n), .true.)

          ! Add: uu_hold += uu
          if (n < nlevs) call saxpy(uu_hold(n), ONE, uu(n), .true.)

          ! Interpolate uu to supply boundary conditions for new 
          ! residual calculation
          call bndry_reg_copy(brs_bcs(n), uu(n-1))
          do i = 1, dm
             call ml_interp_bcs(uu(n), brs_bcs(n)%bmf(i,0), pd, &
                                ref_ratio(n-1,:), -i)
             call ml_interp_bcs(uu(n), brs_bcs(n)%bmf(i,1), pd, &
                                ref_ratio(n-1,:), +i)
          end do
          call multifab_fill_boundary(uu(n))

          ! Compute Res = Res - Lap(uu)
          call mg_defect(mgt(n)%ss(mglev), temp_res(n), res(n), uu(n), &
                         mgt(n)%mm(mglev))
          call multifab_copy(res(n), temp_res(n), ng=res(n)%ng)

          if (do_diagnostics == 1 ) then
             tres = norm_inf(temp_res(n))
             if ( parallel_ioprocessor() ) then
                print *,'UP : RES BEFORE GSRB AT LEVEL ',n, tres
             end if
          end if

          ! Set: uu = 0
          call setval(uu(n), ZERO, all=.true.)

          ! Relax ...
          call mini_cycle(mgt(n), mgt(n)%cycle, mglev, mgt(n)%ss(mglev), &
               uu(n), res(n), mgt(n)%mm(mglev), mgt(n)%nu1, mgt(n)%nu2, &
               mgt(n)%gamma)

          ! Compute Res = Res - Lap(uu)

          if (do_diagnostics == 1 ) then
             call mg_defect(mgt(n)%ss(mglev), temp_res(n), res(n), uu(n), &
                            mgt(n)%mm(mglev))
             call multifab_copy(res(n), temp_res(n), ng=res(n)%ng)
             tres = norm_inf(res(n))
             if ( parallel_ioprocessor() ) then
                print *,'UP : RES AFTER  GSRB AT LEVEL ',n, tres
                if (n == nlevs) print *,' '
             end if
          end if

          ! Add: soln += uu
          call saxpy(soln(n), ONE, uu(n), .true.)

          ! Add: uu += uu_hold so that it will be interpolated too
          if (n < nlevs) call saxpy(uu(n), ONE, uu_hold(n), .true.)

          ! Only do this as long as tangential interp looks under fine grids
          mglev_crse = mgt(n-1)%nlevels
          call ml_restriction(soln(n-1), soln(n), mgt(n)%mm(mglev), &
               mgt(n-1)%mm(mglev_crse), mgt(n)%face_type, ref_ratio(n-1,:))

       end do

       !    Average the solution to the coarser grids.
       do n = nlevs,2,-1
          mglev      = mgt(n)%nlevels
          mglev_crse = mgt(n-1)%nlevels
          call ml_restriction(soln(n-1), soln(n), mgt(n)%mm(mglev), &
               mgt(n-1)%mm(mglev_crse), mgt(n)%face_type, &
               ref_ratio(n-1,:))
       end do

       do n = 1,nlevs
          call multifab_fill_boundary(soln(n), cross = lcross)
       end do

       ! Interpolate soln to supply boundary conditions 
       do n = 2,nlevs
          pd = layout_get_pd(mla%la(n))
          call bndry_reg_copy(brs_bcs(n), soln(n-1))
          do i = 1, dm
             call ml_interp_bcs(soln(n), brs_bcs(n)%bmf(i,0), pd, &
                                ref_ratio(n-1,:), -i)
             call ml_interp_bcs(soln(n), brs_bcs(n)%bmf(i,1), pd, &
                                ref_ratio(n-1,:), +i)
          end do
          call multifab_fill_boundary(soln(n))
       end do

       !    Optimization so don't have to do multilevel convergence test 
       !    each time

       !    Compute the residual on just the finest level
       n = nlevs
       mglev = mgt(n)%nlevels
       call mg_defect(mgt(n)%ss(mglev),res(n),rh(n),soln(n),mgt(n)%mm(mglev))

       if ( ml_fine_converged(res, soln, bnorm, Anorm, eps, abs_eps) ) then

          fine_converged = .true.

          !      Compute the residual on every level
          do n = 1,nlevs-1
             mglev = mgt(n)%nlevels
             call mg_defect(mgt(n)%ss(mglev),res(n),rh(n),soln(n), &
                            mgt(n)%mm(mglev))
          end do

          !      Compute the coarse-fine residual 
          do n = nlevs,2,-1
             pdc = layout_get_pd(mla%la(n-1))
             call crse_fine_residual_cc(n,mgt,soln,res(n-1),brs_flx(n),pdc, &
                                        ref_ratio(n-1,:))
          end do

          !      Average the fine residual onto the coarser level
          do n = nlevs,2,-1
             mglev      = mgt(n  )%nlevels
             mglev_crse = mgt(n-1)%nlevels
             call ml_restriction(res(n-1), res(n), mgt(n)%mm(mglev),&
                  mgt(n-1)%mm(mglev_crse), mgt(n)%face_type, ref_ratio(n-1,:))
          end do

          if ( mgt(nlevs)%verbose > 1 ) then
             do n = 1,nlevs
                tres = norm_inf(res(n))
                if ( parallel_ioprocessor() ) then
!                  write(unit=*, fmt='(i3,": Level ",i2,"  : SL_Ninf(defect) = ",g15.8)') &
!                       iter,n,tres
                   write(unit=*, fmt='("F90mg: Iteration   ",i3," Lev ",i1," error/error0 = ",g15.8)') &
                        iter,n,tres/tres0
                end if
             end do
!            tres = ml_norm_inf(res,fine_mask)
!            if ( parallel_ioprocessor() ) then
!               write(unit=*, fmt='(i3,": All Levels: ML_Ninf(defect) = ",g15.8)') iter, tres
!            end if
          end if

       else 

          fine_converged = .false.
          if ( mgt(nlevs)%verbose > 1 ) then 
             tres = norm_inf(res(nlevs))
             if ( parallel_IOProcessor() ) then
!               write(unit=*, fmt='(i3,": FINE_Ninf(defect) = ",g15.8)') iter, tres
                write(unit=*, fmt='("F90mg: Iteration   ",i3," Fine  error/error0 = ",g15.8)') iter,tres/tres0
             end if
          end if

       end if


     end do

    ! ****************************************************************************

     iter = iter-1
     if (iter < mgt(nlevs)%max_iter) then
        if ( mgt(nlevs)%verbose > 0 ) then
          tres = ml_norm_inf(res,fine_mask)
          if ( parallel_IOProcessor() ) then
             if (tres0 .gt. 0.0_dp_t) then
               write(unit=*, fmt='("F90mg: Final Iter. ",i3," error/error0 = ",g15.8)') iter,tres/tres0
             else
               write(unit=*, fmt='("F90mg: Final Iter. ",i3," error/error0 = ",g15.8)') iter,0.0_dp_t
             end if
          end if
        end if
     else
        call bl_error("Multigrid Solve: failed to converge in max_iter iterations")
     end if

    end if

    ! Add: soln += full_soln
    do n = 1,nlevs
       call saxpy(full_soln(n),ONE,soln(n))
       call multifab_fill_boundary(full_soln(n))
    end do

    do n = 2,nlevs-1
       call multifab_destroy(uu_hold(n))
    end do

    if (need_grad_phi) then

       !   Interpolate boundary conditions of soln in order to get correct grad(phi) at
       !   crse-fine boundaries
       do n = 2,nlevs
          pd = layout_get_pd(mla%la(n))
          call bndry_reg_copy(brs_bcs(n), full_soln(n-1))
          do i = 1, dm
             call ml_interp_bcs(full_soln(n), brs_bcs(n)%bmf(i,0), pd, ref_ratio(n-1,:), -i)
             call ml_interp_bcs(full_soln(n), brs_bcs(n)%bmf(i,1), pd, ref_ratio(n-1,:), +i)
          end do
          call multifab_fill_boundary(full_soln(n))
       end do

    end if

    !   Make sure all periodic and internal boundaries are filled
    do n = 1,nlevs
       call multifab_fill_boundary(full_soln(n))
    end do

    do n = nlevs, 1, -1
       call multifab_destroy(    soln(n))
       call multifab_destroy(      uu(n))
       call multifab_destroy(     res(n))
       call multifab_destroy(temp_res(n))
       if ( n == 1 ) exit
       call bndry_reg_destroy(brs_flx(n))
       call bndry_reg_destroy(brs_bcs(n))
    end do

    call destroy(bpt)

    if ( present(final_resnorm) ) &
       final_resnorm = ni_res

  contains

    function ml_fine_converged(res, sol, bnorm, Anorm, eps, abs_eps) result(r)
      logical :: r
      type(multifab), intent(in) :: res(:), sol(:)
      real(dp_t), intent(in) :: Anorm, eps, abs_eps, bnorm
      real(dp_t) :: ni_res, ni_sol
      integer    :: nlevs
      nlevs = size(res)
      ni_res = norm_inf(res(nlevs))
      ni_sol = norm_inf(sol(nlevs))
      r =  ni_res <= eps*(Anorm*ni_sol + bnorm) .or. &
           ni_res <= abs_eps .or. &
           ni_res <= epsilon(Anorm)*Anorm
    end function ml_fine_converged

    function ml_converged(res, sol, mask, bnorm, Anorm, eps, abs_eps, ni_res, verbose) result(r)

      use ml_util_module

      logical :: r
      integer :: verbose
      type(multifab), intent(in) :: res(:), sol(:)
      type(lmultifab), intent(in) :: mask(:)
      real(dp_t), intent(in   ) :: Anorm, eps, abs_eps, bnorm
      real(dp_t), intent(  out) :: ni_res
      real(dp_t) :: ni_sol

      ni_res = ml_norm_inf(res, mask)
      ni_sol = ml_norm_inf(sol, mask)
      r =  ni_res <= eps*(Anorm*ni_sol + bnorm) .or. &
           ni_res <= abs_eps .or. &
           ni_res <= epsilon(Anorm)*Anorm
      if ( r .and. parallel_IOProcessor() .and. verbose > 1) then
         if (ni_res <= eps*Anorm*ni_sol) then
            print *,'Converged res < eps*Anorm*sol'
         else if (ni_res <= abs_eps) then
            print *,'Converged res < abs_eps '
         else if (ni_res <= eps*bnorm) then
            print *,'Converged res < eps*bnorm '
         else 
            print *,'Converged res < epsilon(Anorm)*Anorm'
         end if
      end if
    end function ml_converged

  end subroutine ml_cc

  subroutine crse_fine_residual_cc(n, mgt, uu, crse_res, brs_flx, pdc, ref_ratio)

      use ml_util_module
      use ml_interface_stencil_module

      integer        , intent(in   ) :: n
      type(mg_tower) , intent(inout) :: mgt(:)
      type(bndry_reg), intent(inout) :: brs_flx
      type(multifab) , intent(inout) :: uu(:)
      type(multifab) , intent(inout) :: crse_res
      type(box)      , intent(in   ) :: pdc
      integer        , intent(in   ) :: ref_ratio(:)

      integer :: i, dm, mglev

      dm = brs_flx%dim
      mglev = mgt(n)%nlevels

      call multifab_fill_boundary(uu(n))

      do i = 1, dm
         call ml_fill_fluxes(mgt(n)%ss(mglev), brs_flx%bmf(i,0), &
              uu(n), mgt(n)%mm(mglev), ref_ratio(i), -1, i)
         call ml_fill_fluxes(mgt(n)%ss(mglev), brs_flx%bmf(i,1), &
              uu(n), mgt(n)%mm(mglev), ref_ratio(i), 1, i)
      end do
      call bndry_reg_copy_to_other(brs_flx)
      do i = 1, dm
         call ml_interface(crse_res, brs_flx%obmf(i,0), uu(n-1), &
              mgt(n-1)%ss(mgt(n-1)%nlevels), pdc, -1, i, ONE)
         call ml_interface(crse_res, brs_flx%obmf(i,1), uu(n-1), &
              mgt(n-1)%ss(mgt(n-1)%nlevels), pdc, +1, i, ONE)
      end do

  end subroutine crse_fine_residual_cc

  subroutine ml_resid(mla, mgt, rh, res, full_soln, ref_ratio)

    use ml_restriction_module
    use ml_prolongation_module

    type(ml_layout), intent(in)    :: mla
    type(mg_tower) , intent(inout) :: mgt(:)
    type( multifab), intent(inout) :: rh(:)
    type( multifab), intent(inout) :: res(:)
    type( multifab), intent(inout) :: full_soln(:)
    integer        , intent(in   ) :: ref_ratio(:,:)

    type(bndry_reg), allocatable :: brs_flx(:)
    type(bndry_reg), allocatable :: brs_bcs(:)

    type(box)    :: pd, pdc
    type(layout) :: la, lac
    integer      :: i, n, dm, nlevs, mglev, mglev_crse

    dm = rh(1)%dim

    nlevs = mla%nlevel

    allocate(brs_bcs(2:nlevs))
    allocate(brs_flx(2:nlevs))

    do n = nlevs, 1, -1

       la = mla%la(n)

       if ( n == 1 ) exit

       ! Build the (coarse resolution) flux registers to be used in computing
       !  the residual at a non-finest AMR level.

       pdc = layout_get_pd(mla%la(n-1))
       lac = mla%la(n-1)
       call bndry_reg_rr_build(brs_flx(n), la, lac, ref_ratio(n-1,:), pdc, width = 0)
       call bndry_reg_rr_build(brs_bcs(n), la, lac, ref_ratio(n-1,:), pdc, width = 2, other = .false.)

    end do

    !  Make sure full_soln at fine grid has the correct coarse grid bc's in its ghost cells 
    !   before we evaluate the initial residual  
    do n = 2,nlevs
       pd = layout_get_pd(mla%la(n))
       call bndry_reg_copy(brs_bcs(n), full_soln(n-1))
       do i = 1, dm
          call ml_interp_bcs(full_soln(n), brs_bcs(n)%bmf(i,0), pd, ref_ratio(n-1,:), -i)
          call ml_interp_bcs(full_soln(n), brs_bcs(n)%bmf(i,1), pd, ref_ratio(n-1,:), +i)
       end do
       call multifab_fill_boundary(full_soln(n))
    end do

    !   Make sure all periodic and internal boundaries are filled as well
    do n = 1,nlevs   
       call multifab_fill_boundary(full_soln(n))
    end do

    do n = 1,nlevs,1
       mglev = mgt(n)%nlevels
       call mg_defect(mgt(n)%ss(mglev),res(n),rh(n),full_soln(n),mgt(n)%mm(mglev))
    end do

    do n = nlevs,2,-1
       mglev      = mgt(n  )%nlevels
       mglev_crse = mgt(n-1)%nlevels

       pdc = layout_get_pd(mla%la(n-1))
       call crse_fine_residual_cc(n,mgt,full_soln,res(n-1),brs_flx(n),pdc,ref_ratio(n-1,:))

       call ml_restriction(res(n-1), res(n), mgt(n)%mm(mglev),&
            mgt(n-1)%mm(mglev_crse), mgt(n)%face_type, ref_ratio(n-1,:))
    enddo

    do n = nlevs, 1, -1
       if ( n == 1 ) exit
       call bndry_reg_destroy(brs_flx(n))
       call bndry_reg_destroy(brs_bcs(n))
    end do

  end subroutine ml_resid

  subroutine ml_cc_applyop(mla, mgt, res, full_soln, ref_ratio)

    use bl_prof_module
    use ml_util_module
    use ml_restriction_module, only: ml_restriction
    use ml_prolongation_module, only: ml_prolongation, ml_interp_bcs

    type(ml_layout), intent(in)    :: mla
    type(mg_tower) , intent(inout) :: mgt(:)
    type( multifab), intent(inout) :: res(:)
    type( multifab), intent(inout) :: full_soln(:)
    integer        , intent(in   ) :: ref_ratio(:,:)

    integer :: nlevs
    type(multifab), allocatable  ::      soln(:)
    type(multifab), allocatable  ::        uu(:)
    type(multifab), allocatable  ::   uu_hold(:)
    type(multifab), allocatable  ::        rh(:) ! this will be set to zero
    type(multifab), allocatable  ::  temp_res(:)

    type(bndry_reg), allocatable :: brs_bcs(:)

    type(box) :: pd, pdc
    type(layout) :: la, lac
    integer :: i, n, dm
    integer :: mglev, mglev_crse

    type(bl_prof_timer), save :: bpt
    integer                   :: lo(res(1)%dim),hi(res(1)%dim),ng
    real(kind=dp_t),  pointer :: resp(:,:,:,:)

    call build(bpt, "ml_cc_applyop")

    nlevs = mla%nlevel

    allocate(soln(nlevs), uu(nlevs), rh(nlevs), temp_res(nlevs))
    allocate(uu_hold(2:nlevs-1))
    allocate(brs_bcs(2:nlevs))

    do n = 2,nlevs-1
       la = mla%la(n)
       call build(uu_hold(n),la,1,1)
       call setval( uu_hold(n), ZERO,all=.true.)
    end do

    do n = nlevs, 1, -1

       la = mla%la(n)
       call build(    soln(n), la, 1, 1)
       call build(      uu(n), la, 1, 1)
       call build(      rh(n), la, 1, 0)
       call build(temp_res(n), la, 1, 0)
       call setval(    soln(n), ZERO,all=.true.)
       call setval(      uu(n), ZERO,all=.true.)
       call setval(      rh(n), ZERO,all=.true.)
       call setval(temp_res(n), ZERO,all=.true.)

       ! zero residual just to be safe
       call setval(     res(n), ZERO,all=.true.)

       if ( n == 1 ) exit

       ! Build the (coarse resolution) flux registers to be used in computing
       !  the residual at a non-finest AMR level.

       pdc = layout_get_pd(mla%la(n-1))
       lac = mla%la(n-1)
       call bndry_reg_rr_build(brs_bcs(n), la, lac, ref_ratio(n-1,:), pdc, width = 2, other = .false.)

    end do

    dm = rh(1)%dim

    !  Make sure full_soln at fine grid has the correct coarse grid bc's in 
    !  its ghost cells before we evaluate the initial residual  
    do n = 2,nlevs
       pd = layout_get_pd(mla%la(n))
       call bndry_reg_copy(brs_bcs(n), full_soln(n-1))
       do i = 1, dm
          call ml_interp_bcs(full_soln(n), brs_bcs(n)%bmf(i,0), pd, &
                             ref_ratio(n-1,:), -i)
          call ml_interp_bcs(full_soln(n), brs_bcs(n)%bmf(i,1), pd, &
                             ref_ratio(n-1,:), +i)
       end do
       call multifab_fill_boundary(full_soln(n))
    end do

    !   Make sure all periodic and internal boundaries are filled as well
    do n = 1,nlevs   
       call multifab_fill_boundary(full_soln(n))
    end do


    do n = 1,nlevs,1
       mglev = mgt(n)%nlevels
       if (n.eq.1) call print(mgt(n)%ss(mglev),'SS')
       call mg_defect(mgt(n)%ss(mglev),res(n),rh(n),full_soln(n), &
                      mgt(n)%mm(mglev))
    end do

    call print(res(1),'RES')
    stop

    ! still need to multiply residual by -1 to get (alpha - del dot beta grad)
    do n=1,nlevs
       ng = res(n)%ng
       
       do i=1,res(n)%nboxes
          if (multifab_remote(res(n),i)) cycle
          resp  => dataptr(res(n),i)
          lo =  lwb(get_box(res(n), i))
          hi =  upb(get_box(res(n), i))
          select case (dm)
          case (1)
             call scale_residual_1d(lo,hi,ng,resp(:,1,1,1))
          case (2)
             call scale_residual_2d(lo,hi,ng,resp(:,:,1,1))
          case (3)
             call scale_residual_3d(lo,hi,ng,resp(:,:,:,1))
          end select
       end do
    enddo

    do n = 2,nlevs-1
       call destroy(uu_hold(n))
    end do

    do n = nlevs, 1, -1
       call destroy(soln(n))
       call destroy(uu(n))
       call destroy(rh(n))
       call destroy(temp_res(n))
       if ( n == 1 ) exit
       call bndry_reg_destroy(brs_bcs(n))
    end do

    call destroy(bpt)

  end subroutine ml_cc_applyop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Multiply residual by -1 in 1d
  subroutine scale_residual_1d(lo,hi,ng,res)

  integer        , intent(in   ) :: lo(:),hi(:),ng
  real(kind=dp_t), intent(inout) :: res(lo(1)-ng:)

! Local
  integer :: i

  do i=lo(1),hi(1)
     res(i) = -res(i)
  enddo
  
  end subroutine scale_residual_1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Multiply residual by -1 in 2d
  subroutine scale_residual_2d(lo,hi,ng,res)

  integer        , intent(in   ) :: lo(:),hi(:),ng
  real(kind=dp_t), intent(inout) :: res(lo(1)-ng:,lo(2)-ng:)

! Local
  integer :: i,j

  do j=lo(2),hi(2)
     do i=lo(1),hi(1)
        res(i,j) = -res(i,j)
     enddo
  enddo
  
  end subroutine scale_residual_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Multiply residual by -1 in 3d
  subroutine scale_residual_3d(lo,hi,ng,res)

  integer        , intent(in   ) :: lo(:),hi(:),ng
  real(kind=dp_t), intent(inout) :: res(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)

! Local
  integer :: i,j,k

  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)
           res(i,j,k) = -res(i,j,k)
        enddo
     enddo
  enddo

  end subroutine scale_residual_3d

end module ml_cc_module
