module ml_nd_module

  use bl_constants_module
  use bl_prof_module
  use mg_module
  use fabio_module
  use ml_layout_module
  use bndry_reg_module

  implicit none

contains

  subroutine ml_nd(mla,mgt,rh,full_soln,fine_mask,do_diagnostics)

    use ml_norm_module             , only : ml_norm_inf
    use ml_nodal_restriction_module, only : ml_nodal_restriction, periodic_add_copy
    use ml_prolongation_module     , only : ml_nodal_prolongation

    type(ml_layout), intent(in   ) :: mla
    type(mg_tower ), intent(inout) :: mgt(:)
    type( multifab), intent(inout) :: rh(:)
    type( multifab), intent(inout) :: full_soln(:)
    type(lmultifab), intent(in   ) :: fine_mask(:)
    integer        , intent(in   ) :: do_diagnostics 

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

    real(dp_t) :: Anorm, bnorm, fac, tres, ttres, tres0, t1(3), t2(3)
    real(dp_t) :: mltres(mla%nlevel), lmltres(mla%nlevel)
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
       call bndry_reg_rr_build_nd(brs_flx(n), la, mla%mba%rr(n-1,:), pdc, nodal)

    end do
    !
    ! Let's elide some reductions by doing these reductions together.
    !
    bnorm = ml_norm_inf(rh,fine_mask,local=.true.)

    Anorm = nodal_stencil_norm(mgt(nlevs)%ss(mgt(nlevs)%nlevels),mgt(n)%stencil_type,&
                               mgt(n)%uniform_dh, local=.true.)
    do n = 1, nlevs-1
       Anorm = max(nodal_stencil_norm(mgt(n)%ss(mgt(n)%nlevels), &
                   mgt(n)%stencil_type, mgt(n)%uniform_dh, &
                   fine_mask(n), local=.true.), Anorm)
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
                      mgt(nlevs)%eps, mgt(nlevs)%abs_eps, mgt(nlevs)%verbose) ) then
       if ( parallel_IOProcessor() .and. mgt(nlevs)%verbose > 0 ) &
            write(unit=*, fmt='("F90mg: No iterations needed ")')

       n_mg_iters = 0

    else   ! Not already converged

     do iter = 1, mgt(nlevs)%max_iter

       if ( (iter .eq. 1) .or. fine_converged ) then
          if ( ml_converged(res, fine_mask, bnorm, mgt(nlevs)%eps, mgt(nlevs)%abs_eps, mgt(nlevs)%verbose) ) exit
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
             if (mgt(n)%cycle_type == MG_FVCycle) then
                if (iter == 1) then
                    call mg_tower_cycle(mgt(n), MG_FCycle, mglev, mgt(n)%ss(mglev), &
                         uu(n), res(n), mgt(n)%mm(mglev), mgt(n)%nu1, mgt(n)%nu2, &
                         bottom_solve_time = bottom_solve_time)
                else
                    call mg_tower_cycle(mgt(n), MG_VCycle, mglev, mgt(n)%ss(mglev), &
                         uu(n), res(n), mgt(n)%mm(mglev), mgt(n)%nu1, mgt(n)%nu2, &
                         bottom_solve_time = bottom_solve_time)
                end if
             else 
                call mg_tower_cycle(mgt(n), mgt(n)%cycle_type, mglev, mgt(n)%ss(mglev), &
                     uu(n), res(n), mgt(n)%mm(mglev), mgt(n)%nu1, mgt(n)%nu2, &
                     bottom_solve_time = bottom_solve_time)
             end if
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
                fac = EIGHT**(mla%mba%rr(n-1,1)/2)
                call multifab_mult_mult_s(res(n-1),fac,nghost(res(n-1)))
             else if ( dm .eq. 2 ) then
                fac = FOUR**(mla%mba%rr(n-1,1)/2)
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
             call ml_nodal_restriction(res(n-1), res(n), mgt(n)%mm(mglev),&
                                       mgt(n-1)%mm(mglev_crse), mla%mba%rr(n-1,:))

             ! Compute CRSE-FINE Res = Rh - Lap(Soln)
             pdc = layout_get_pd(mla%la(n-1))
             call crse_fine_residual_nodal(n,mgt,brs_flx(n),res(n-1),zero_rh(n),temp_res(n),temp_res(n-1), &
                  soln(n-1),soln(n),mla%mba%rr(n-1,:),pdc)

             ! Copy u_hold = uu
             if (n < nlevs) call multifab_copy(uu_hold(n),uu(n),ng=nghost(uu(n)))

             ! Set: uu = 0
             call setval(uu(n),ZERO,all=.true.)

             if ( dm .eq. 3 ) then
                fac = ONE / EIGHT**(mla%mba%rr(n-1,1)/2)
                call multifab_mult_mult_s(res(n-1),fac,nghost(res(n-1)))
             else if ( dm .eq. 2 ) then
                fac = ONE / FOUR**(mla%mba%rr(n-1,1)/2)
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
          call ml_nodal_prolongation(uu(n), uu(n-1), mla%mba%rr(n-1,:))
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
          call ml_nodal_restriction(soln(n-1), soln(n), mgt(n)%mm(mglev), &
                                    mgt(n-1)%mm(mglev_crse), mla%mba%rr(n-1,:), inject = .true.)
       end do

        do n = 1,nlevs
           call multifab_fill_boundary(soln(n), cross = mgt(n)%lcross)
        end do

       !    Optimization so don't have to do multilevel convergence test each time

       !    Compute the residual on just the finest level
        n = nlevs
        mglev = mgt(n)%nlevels
        call compute_defect(mgt(n)%ss(mglev),res(n),rh(n),soln(n),mgt(n)%mm(mglev), &
                       mgt(n)%stencil_type, mgt(n)%lcross, mgt(n)%uniform_dh, filled=.true.)

        if ( ml_fine_converged(res, bnorm, mgt(nlevs)%eps, mgt(nlevs)%abs_eps) ) then

          fine_converged = .true.

          !      Compute the residual on every level
          do n = 1,nlevs-1
             mglev = mgt(n)%nlevels
             call compute_defect(mgt(n)%ss(mglev),res(n),rh(n),soln(n),mgt(n)%mm(mglev), &
                            mgt(n)%stencil_type, mgt(n)%lcross, mgt(n)%uniform_dh, filled=.true.)
          end do

          do n = nlevs,2,-1
             !  Restrict the finer residual onto the coarser grid
             mglev      = mgt(n  )%nlevels
             mglev_crse = mgt(n-1)%nlevels
             if ( dm .eq. 3 ) then
                fac = EIGHT**(mla%mba%rr(n-1,1)/2)
                call multifab_mult_mult_s(res(n-1),fac,nghost(res(n-1)))
             else if ( dm .eq. 2 ) then
                fac = FOUR**(mla%mba%rr(n-1,1)/2)
                call multifab_mult_mult_s(res(n-1),fac,nghost(res(n-1)))
             end if
             call ml_nodal_restriction(res(n-1), res(n), mgt(n)%mm(mglev),&
                  mgt(n-1)%mm(mglev_crse), mla%mba%rr(n-1,:))

             !  Compute the coarse-fine residual at coarse-fine nodes
             pdc = layout_get_pd(mla%la(n-1))
             call crse_fine_residual_nodal(n,mgt,brs_flx(n),res(n-1), &
                  zero_rh(n),temp_res(n),temp_res(n-1), &
                  soln(n-1),soln(n),mla%mba%rr(n-1,:),pdc,filled=.true.) ! soln(n) is filled
             if ( dm .eq. 3 ) then
                fac = ONE / EIGHT**(mla%mba%rr(n-1,1)/2)
                call multifab_mult_mult_s(res(n-1),fac,nghost(res(n-1)))
             else if ( dm .eq. 2 ) then
                fac = ONE / FOUR**(mla%mba%rr(n-1,1)/2)
                call multifab_mult_mult_s(res(n-1),fac,nghost(res(n-1)))
             end if
          end do

          if ( mgt(nlevs)%verbose > 1 ) then
             do n = 1,nlevs
                lmltres(n) = norm_inf(res(n),local=.true.)
             end do
             call parallel_reduce(mltres, lmltres, MPI_MAX, proc = parallel_IOProcessorNode())
             do n = 1,nlevs
                if ( parallel_IOProcessor() ) then
                   write(unit=*, fmt='("F90mg: Iteration   ",i3," Lev ",i1," resid/resid0 = ",g15.8)') &
                        iter,n,mltres(n)/tres0
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

           n_mg_iters = iter

           if ( parallel_IOProcessor() ) then
              if ( tres0 .gt. ZERO) then
                 write(unit=*, fmt='("F90mg: Final Iter. ",i3," resid/resid0 = ",g15.8)') iter,t2(1)/tres0
                 write(unit=*, fmt='("F90mg: Solve time: ",g13.6, " Bottom Solve time: ", g13.6)') t2(2), t2(3)
                 write(unit=*, fmt='("")')
              else
                 write(unit=*, fmt='("F90mg: Final Iter. ",i3," resid/resid0 = ",g15.8)') iter,ZERO
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
         crse_soln,fine_soln,ref_ratio,pdc,filled)

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
      integer        , intent(in   ) :: ref_ratio(:)
      type(box)      , intent(in   ) :: pdc
      logical, intent(in), optional  :: filled  ! is fine_soln filled?

      integer :: i,dm,mglev_crse,mglev_fine

      mglev_crse = mgt(n-1)%nlevels
      mglev_fine = mgt(n  )%nlevels
      dm        = get_dim(temp_res)

      !    Compute the fine contributions at faces, edges and corners.

      !    First compute a residual which only takes contributions from the
      !       grid on which it is calculated.

      call grid_res(mgt(n)%ss(mglev_fine),temp_res, &
                    fine_rhs,fine_soln,mgt(n)%mm(mglev_fine), &
                    mgt(n)%lcross,  mgt(n)%stencil_type, mgt(n)%uniform_dh, filled)

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
              pdc,ref_ratio, mgt(n)%stencil_type, -i)
         call ml_crse_contrib(temp_crse_res, brs_flx%bmf(i,1), crse_soln, &
              mgt(n-1)%ss(mgt(n-1)%nlevels), &
              mgt(n-1)%mm(mglev_crse), &
              mgt(n  )%mm(mglev_fine), &
              pdc,ref_ratio, mgt(n)%stencil_type, +i)
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
      real(dp_t) :: ni_res, l_ni_res
      l_ni_res = ml_norm_inf(res, mask, local=.true.)
      call parallel_reduce(ni_res, l_ni_res, MPI_MAX) 
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

  subroutine grid_res(ss, dd, ff, uu, mm, lcross, stencil_type, uniform_dh, filled)

    type(multifab), intent(in)    :: ff, ss
    type(multifab), intent(inout) :: dd, uu
    type(imultifab), intent(in)   :: mm
    logical, intent(in) :: lcross
    integer, intent(in) :: stencil_type
    logical, intent(in) :: uniform_dh
    logical, intent(in), optional :: filled

    integer :: i, n
    real(kind=dp_t), pointer :: dp(:,:,:,:)
    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: up(:,:,:,:)
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    integer        , pointer :: mp(:,:,:,:)
    integer :: dm,nodal_ng
    logical :: lfilled
    type(bl_prof_timer), save :: bpt

    call build(bpt, "grid_res")

    lfilled = .false.;  if (present(filled)) lfilled = filled

    nodal_ng = 0; if ( nodal_q(uu) ) nodal_ng = 1

    dm = get_dim(uu)

    if (.not.lfilled) call multifab_fill_boundary(uu, cross = lcross)

    do i = 1, nfabs(uu)
       dp  => dataptr(dd, i)
       fp  => dataptr(ff, i)
       up  => dataptr(uu, i)
       sp  => dataptr(ss, i)
       mp  => dataptr(mm, i)
       do n = 1, ncomp(uu)
          select case(dm)
          case (1)
             call grid_laplace_1d(sp(1,:,1,1), dp(:,1,1,n), fp(:,1,1,n), up(:,1,1,n), &
                                  mp(:,1,1,1),nghost(uu))
          case (2)
             call grid_laplace_2d(sp(1,:,:,1), dp(:,:,1,n), fp(:,:,1,n), up(:,:,1,n), &
                                  mp(:,:,1,1), nghost(uu), stencil_type)
          case (3)
             call grid_laplace_3d(sp(1,:,:,:), dp(:,:,:,n), fp(:,:,:,n), up(:,:,:,n), &
                                  mp(:,:,:,1), nghost(uu), stencil_type, uniform_dh)
          end select
       end do
    end do

    call destroy(bpt)

  end subroutine grid_res

  subroutine grid_laplace_1d(sg, dd, ff, uu, mm, ng)

    use impose_neumann_bcs_module

    integer           , intent(in   ) :: ng
    real (kind = dp_t), intent(in   ) :: sg(0:)
    real (kind = dp_t), intent(inout) :: dd(0:)
    real (kind = dp_t), intent(in   ) :: ff(0:)
    integer,            intent(in   ) :: mm(:)
    real (kind = dp_t), intent(inout) :: uu(1-ng:)

    real (kind = dp_t), allocatable :: sg_int(:)

    integer :: i,nx,lo(1)
    nx = size(sg,dim=1)-2

    lo = 1
    call impose_neumann_bcs_1d(uu,mm,lo,ng)

    allocate(sg_int(0:nx+1))
    sg_int(:) = ZERO 

    ! Copy the interior values only, we do *not* want to fill ghost cells that 
    !      intersect other fine grids.
    sg_int(1:nx) = sg(1:nx)

    do i = 1,nx+1
      dd(i) = ff(i) - &
        (sg_int(i) * (uu(i+1) - uu(i)) + sg_int(i-1) * (uu(i-1) - uu(i)))
    end do
    
    deallocate(sg_int)

  end subroutine grid_laplace_1d

  subroutine grid_laplace_2d(sg, dd, ff, uu, mm, ng, stencil_type)

    use bc_module
    use impose_neumann_bcs_module
    use nodal_stencil_module

    integer           , intent(in   ) :: ng
    real (kind = dp_t), intent(in   ) :: sg(0:,0:)
    real (kind = dp_t), intent(inout) :: dd(0:,0:) 
    real (kind = dp_t), intent(in   ) :: ff(0:,0:)
    real (kind = dp_t), intent(inout) :: uu(1-ng:,1-ng:)
    integer,            intent(in   ) :: mm(:,:)
    integer           , intent(in   ) :: stencil_type
    integer :: i,j,nx,ny,lo(2),hi(2)

    real (kind = dp_t), allocatable :: sg_int(:,:)

    nx = size(sg,dim=1)-2
    ny = size(sg,dim=2)-2

    lo(:) = 1
    hi(1) = nx
    hi(2) = ny

    call impose_neumann_bcs_2d(uu,mm,lo,ng)

    allocate(sg_int(0:nx+1,0:ny+1))
    sg_int(:,:) = ZERO 

    ! Copy the interior values only, we do *not* want to fill ghost cells that 
    !      intersect other fine grids.
    do j = 1,ny
    do i = 1,nx
        sg_int(i,j) = sg(i,j)
    end do
    end do

    ! Set values across Neumann boundaries
    call set_faces_edges_corners_2d(sg_int, mm, lo, hi)

    if (stencil_type .eq. ND_VATER_STENCIL) then
        call bl_error('ND_VATER_STENCIL not implemented in grid_laplace_2d')
    else if (stencil_type .eq. ND_DENSE_STENCIL) then
 
!     Interior
      do j = 1,ny+1
      do i = 1,nx+1
         dd(i,j) = ff(i,j) - THIRD * ( &
                   sg_int(i-1,j-1) * uu(i-1,j-1) + &
                   sg_int(i  ,j-1) * uu(i+1,j-1) + &
                   sg_int(i-1,j  ) * uu(i-1,j+1) + &
                   sg_int(i  ,j  ) * uu(i+1,j+1) + &
           HALF * ( (sg_int(i-1,j-1) + sg_int(i  ,j-1)) * uu(i,j-1) + &
                    (sg_int(i-1,j-1) + sg_int(i-1,j  )) * uu(i-1,j) + &
                    (sg_int(i  ,j-1) + sg_int(i  ,j  )) * uu(i+1,j) + &
                    (sg_int(i-1,j  ) + sg_int(i  ,j  )) * uu(i,j+1) ) - &
            TWO * (sg_int(i-1,j-1) + sg_int(i,j-1) + sg_int(i-1,j) + sg_int(i,j)) * uu(i,j) )
      end do
      end do

    else if (stencil_type .eq. ND_CROSS_STENCIL) then

!     Interior
      do j = 1,ny+1
      do i = 1,nx+1
         dd(i,j) = ff(i,j) - &
                   HALF * ( (sg_int(i  ,j-1)+sg_int(i  ,j  )) * uu(i+1,j) + &
                            (sg_int(i-1,j-1)+sg_int(i-1,j  )) * uu(i-1,j) + &
                            (sg_int(i-1,j  )+sg_int(i  ,j  )) * uu(i,j+1) + &
                            (sg_int(i-1,j-1)+sg_int(i  ,j-1)) * uu(i,j-1) ) + &
            (sg_int(i-1,j-1)+sg_int(i,j-1)+sg_int(i-1,j)+sg_int(i,j)) * uu(i,j)   
      end do
      end do

    else 
        call bl_error('unknown stencil_type in grid_laplace_2d')
    end if

    deallocate(sg_int)

  end subroutine grid_laplace_2d

  subroutine grid_laplace_3d(sg, dd, ff, uu, mm, ng, stencil_type, uniform_dh)
    use bc_module
    use impose_neumann_bcs_module
    use nodal_stencil_module
    integer           , intent(in   ) :: ng
    real (kind = dp_t), intent(in   ) :: sg(0:,0:,0:)
    real (kind = dp_t), intent(inout) :: dd(0:,0:,0:)
    real (kind = dp_t), intent(in   ) :: ff(0:,0:,0:)
    real (kind = dp_t), intent(inout) :: uu(1-ng:,1-ng:,1-ng:)
    integer,            intent(in   ) :: mm(:,:,:)
    logical,            intent(in)    :: uniform_dh
    integer           , intent(in   ) :: stencil_type

    real (kind = dp_t), allocatable   :: sg_int(:,:,:)
    real (kind = dp_t) :: f0, fx, fy, fz, fxyz, f2y2zx, f2x2zy, f2x2yz, ss0
    real (kind = dp_t) :: ff_fac

    integer :: i, j, k, lo(3), hi(3)
    integer :: nx, ny, nz

    nx = size(sg,dim=1)-2
    ny = size(sg,dim=2)-2
    nz = size(sg,dim=3)-2

    lo(:) = 1
    hi(1) = nx
    hi(2) = ny
    hi(3) = nz
    call impose_neumann_bcs_3d(uu,mm,lo,ng)

    allocate(sg_int(0:nx+1,0:ny+1,0:nz+1))
    sg_int(:,:,:) = ZERO 

    ! Copy the interior values only, we do *not* want to fill ghost cells that 
    !      intersect other fine grids.
    sg_int(1:nx,1:ny,1:nz) = sg(1:nx,1:ny,1:nz)

    ! Set values across Neumann boundaries
    call set_faces_edges_corners_3d(sg_int, mm, lo, hi)

    fx     = ONE/36._dp_t
    fy     = fx
    fz     = fx
    f0     = FOUR * (fx + fy + fz)
    fxyz   = (fx+fy+fz)
    f2y2zx = (TWO*fy+TWO*fz-fx)
    f2x2zy = (TWO*fx+TWO*fz-fy)
    f2x2yz = (TWO*fx+TWO*fy-fz)

    if (stencil_type .eq. ND_VATER_STENCIL) then
        call bl_error('ND_VATER_STENCIL not implemented in grid_laplace_3d')
    else if (stencil_type .eq. ND_DENSE_STENCIL) then
      !$OMP PARALLEL DO PRIVATE(i,j,k,ff_fac,ss0)
      do k = 1,nz+1
      do j = 1,ny+1
      do i = 1,nx+1

          ff_fac = ONE
          if (i.eq.1 .or. i.eq.nx+1) ff_fac = HALF * ff_fac
          if (j.eq.1 .or. j.eq.ny+1) ff_fac = HALF * ff_fac
          if (k.eq.1 .or. k.eq.nz+1) ff_fac = HALF * ff_fac

          ss0 =  -( sg_int(i-1,j-1,k-1) + sg_int(i,j-1,k-1) &
                   +sg_int(i-1,j  ,k-1) + sg_int(i,j  ,k-1) &
                   +sg_int(i-1,j-1,k  ) + sg_int(i,j-1,k  ) &
                   +sg_int(i-1,j  ,k  ) + sg_int(i,j  ,k  ) ) * f0

          dd(i,j,k) = ff_fac * ff(i,j,k) - ( &
               fxyz * ( &   ! Corners
               sg_int(i-1,j-1,k-1) * uu(i-1,j-1,k-1) + sg_int(i  ,j-1,k-1) * uu(i+1,j-1,k-1) + &
               sg_int(i-1,j  ,k-1) * uu(i-1,j+1,k-1) + sg_int(i  ,j  ,k-1) * uu(i+1,j+1,k-1) + &
               sg_int(i-1,j-1,k  ) * uu(i-1,j-1,k+1) + sg_int(i  ,j-1,k  ) * uu(i+1,j-1,k+1) + &
               sg_int(i-1,j  ,k  ) * uu(i-1,j+1,k+1) + sg_int(i  ,j  ,k  ) * uu(i+1,j+1,k+1)) &
               + f2y2zx * ( & ! Edges in x-direction
               (sg_int(i  ,j-1,k-1) + sg_int(i-1,j-1,k-1)) * uu(i  ,j-1,k-1) + &
               (sg_int(i  ,j  ,k-1) + sg_int(i-1,j  ,k-1)) * uu(i  ,j+1,k-1) + &
               (sg_int(i  ,j-1,k  ) + sg_int(i-1,j-1,k  )) * uu(i  ,j-1,k+1) + &
               (sg_int(i  ,j  ,k  ) + sg_int(i-1,j  ,k  )) * uu(i  ,j+1,k+1)) &
               + f2x2zy * ( & ! Edges in y-direction
               (sg_int(i-1,j-1,k-1) + sg_int(i-1,j  ,k-1)) * uu(i-1,j  ,k-1) + &
               (sg_int(i  ,j-1,k-1) + sg_int(i  ,j  ,k-1)) * uu(i+1,j  ,k-1) + &
               (sg_int(i-1,j-1,k  ) + sg_int(i-1,j  ,k  )) * uu(i-1,j  ,k+1) + &
               (sg_int(i  ,j-1,k  ) + sg_int(i  ,j  ,k  )) * uu(i+1,j  ,k+1)) &
               + f2x2yz * ( & ! Edges in z-direction
               (sg_int(i-1,j-1,k-1) + sg_int(i-1,j-1,k  )) * uu(i-1,j-1,k  ) + &
               (sg_int(i  ,j-1,k-1) + sg_int(i  ,j-1,k  )) * uu(i+1,j-1,k  ) + &
               (sg_int(i-1,j  ,k-1) + sg_int(i-1,j  ,k  )) * uu(i-1,j+1,k  ) + &
               (sg_int(i  ,j  ,k-1) + sg_int(i  ,j  ,k  )) * uu(i+1,j+1,k  )) &
               + ss0 * uu(i,j,k) )
  
          if (.not. uniform_dh) then
          !
          ! Add faces (only non-zero for non-uniform dx)
          !   (minus sign here because this is really (ff - lap(u))
          !
          dd(i,j,k) = dd(i,j,k) - ( &
               ( (FOUR*fx-TWO*fy-TWO*fz)*(sg_int(i-1,j-1,k-1) + sg_int(i-1,j-1,k  ) &
                +sg_int(i-1,j  ,k-1) + sg_int(i-1,j  ,k  )) ) * uu(i-1,j  ,k  ) + &
               ( (FOUR*fx-TWO*fy-TWO*fz)*(sg_int(i  ,j-1,k-1) + sg_int(i  ,j-1,k  ) &
                +sg_int(i  ,j  ,k-1) + sg_int(i  ,j  ,k  )) ) * uu(i+1,j  ,k  ) + &
               ( (FOUR*fy-TWO*fx-TWO*fz)*(sg_int(i-1,j-1,k-1) + sg_int(i-1,j-1,k  ) &
                +sg_int(i  ,j-1,k-1) + sg_int(i  ,j-1,k  )) ) * uu(i  ,j-1,k  ) + &
               ( (FOUR*fy-TWO*fx-TWO*fz)*(sg_int(i-1,j  ,k-1) + sg_int(i-1,j  ,k  ) &
                +sg_int(i  ,j  ,k-1) + sg_int(i  ,j  ,k  )) ) * uu(i  ,j+1,k  ) + &
               ( (FOUR*fz-TWO*fx-TWO*fy)*(sg_int(i-1,j-1,k-1) + sg_int(i-1,j  ,k-1) &
                +sg_int(i  ,j-1,k-1) + sg_int(i  ,j  ,k-1)) ) * uu(i  ,j  ,k-1) + &
               ( (FOUR*fz-TWO*fx-TWO*fy)*(sg_int(i-1,j-1,k  ) + sg_int(i-1,j  ,k  ) &
                +sg_int(i  ,j-1,k  ) + sg_int(i  ,j  ,k  )) ) * uu(i  ,j  ,k+1)     )

          end if

          ! This accounts for the fact that fac = 1/(4*dx*dx) to be compatible with
          !      the cross stencil
          dd(i,j,k) = FOUR * dd(i,j,k)
      end do
      end do
      end do
      !$OMP END PARALLEL DO

    else if (stencil_type .eq. ND_CROSS_STENCIL) then

      !$OMP PARALLEL DO PRIVATE(i,j,k,ff_fac)
      do k = 1,nz+1
      do j = 1,ny+1
      do i = 1,nx+1

          ff_fac = ONE
          if (i.eq.1 .or. i.eq.nx+1) ff_fac = HALF * ff_fac
          if (j.eq.1 .or. j.eq.ny+1) ff_fac = HALF * ff_fac
          if (k.eq.1 .or. k.eq.nz+1) ff_fac = HALF * ff_fac

          dd(i,j,k) = ff_fac*ff(i,j,k) - ( &
             (sg_int(i-1,j-1,k-1) + sg_int(i-1,j-1,k  ) &
             +sg_int(i-1,j  ,k-1) + sg_int(i-1,j  ,k  )) * uu(i-1,j  ,k  ) + &
             (sg_int(i  ,j-1,k-1) + sg_int(i  ,j-1,k  ) &
             +sg_int(i  ,j  ,k-1) + sg_int(i  ,j  ,k  )) * uu(i+1,j  ,k  ) + &
             (sg_int(i-1,j-1,k-1) + sg_int(i-1,j-1,k  ) &
             +sg_int(i  ,j-1,k-1) + sg_int(i  ,j-1,k  )) * uu(i  ,j-1,k  ) + &
             (sg_int(i-1,j  ,k-1) + sg_int(i-1,j  ,k  ) &
             +sg_int(i  ,j  ,k-1) + sg_int(i  ,j  ,k  )) * uu(i  ,j+1,k  ) + &
             (sg_int(i-1,j-1,k-1) + sg_int(i-1,j  ,k-1) &
             +sg_int(i  ,j-1,k-1) + sg_int(i  ,j  ,k-1)) * uu(i  ,j  ,k-1) + &
             (sg_int(i-1,j-1,k  ) + sg_int(i-1,j  ,k  ) &
             +sg_int(i  ,j-1,k  ) + sg_int(i  ,j  ,k  )) * uu(i  ,j  ,k+1) - &
             (sg_int(i-1,j-1,k-1) + sg_int(i-1,j  ,k-1) &
             +sg_int(i  ,j-1,k-1) + sg_int(i  ,j  ,k-1) &
             +sg_int(i-1,j-1,k  ) + sg_int(i-1,j  ,k  ) &
             +sg_int(i  ,j-1,k  ) + sg_int(i  ,j  ,k  )) * THREE * uu(i  ,j  ,k  ) )
      end do
      end do
      end do
      !$OMP END PARALLEL DO

    else 
        call bl_error('unknown stencil_type in grid_laplace_3d')
    end if

    deallocate(sg_int)

  end subroutine grid_laplace_3d

  function nodal_stencil_norm(sg, stencil_type, uniform_dh, mask, local) result(r)

    use bl_prof_module
    real(kind=dp_t) :: r
    type( multifab), intent(inout)        :: sg
    integer        , intent(in)           :: stencil_type
    type(lmultifab), intent(in), optional :: mask
    logical        , intent(in), optional :: uniform_dh, local

    ! Local variables
    integer                   :: i,j,k,b,dm
    real(kind=dp_t)           :: r1, sum_comps, ss0
    real (kind = dp_t) :: f0, fx, fy, fz, fxyz, f2y2zx, f2x2zy, f2x2yz
    real(kind=dp_t), pointer  :: sp(:,:,:,:)
    logical        , pointer  :: lp(:,:,:,:)
    logical                   :: llocal
    type(bl_prof_timer), save :: bpt

    type(box) :: bx

    call build(bpt, "st_norm_nd")

    llocal = .false.; if ( present(local) ) llocal = local

    dm = get_dim(sg)

    if (.not. cell_centered_q(sg)) &
        call bl_error("sg is supposed to be cell-centered in nodal_stencil_norm")

    r1 = -Huge(r1)

    if (stencil_type .eq. ND_CROSS_STENCIL) then

      if (present(mask)) then

       !$OMP PARALLEL DO PRIVATE(i,j,k,bx,sum_comps,ss0,sp,lp) REDUCTION(max:r1)
       do b = 1, nfabs(sg)
          bx = get_box(sg,b)

          sp => dataptr(sg, b)
          lp => dataptr(mask, b)

          if (dm.eq.2) then
             do j = bx%lo(2), bx%hi(2)+1
             do i = bx%lo(1), bx%hi(1)+1
                if ( lp(i,j,1,1) ) then
                   ss0 = -( sp(1,i-1,j-1,1)+sp(1,i,j-1,1) &
                           +sp(1,i-1,j,1)+sp(1,i,j,1) )
                   sum_comps = abs(ss0) + HALF * ( &
                      abs(sp(1,i  ,j-1,1)+sp(1,i  ,j  ,1)) + &
                      abs(sp(1,i-1,j-1,1)+sp(1,i-1,j  ,1)) + & 
                      abs(sp(1,i-1,j  ,1)+sp(1,i  ,j  ,1)) + &
                      abs(sp(1,i-1,j-1,1)+sp(1,i  ,j-1,1)) )
                   r1 = max(r1,sum_comps)
                end if ! end mask test
             end do
             end do

          else if (dm.eq.3) then

             do k = bx%lo(3), bx%hi(3)+1
             do j = bx%lo(2), bx%hi(2)+1
             do i = bx%lo(1), bx%hi(1)+1
                if ( lp(i,j,k,1) ) then
                   ss0 = -( sp(1,i-1,j-1,k-1) + sp(1,i-1,j  ,k-1) &
                           +sp(1,i  ,j-1,k-1) + sp(1,i  ,j  ,k-1) &
                           +sp(1,i-1,j-1,k  ) + sp(1,i-1,j  ,k  ) &
                           +sp(1,i  ,j-1,k  ) + sp(1,i  ,j  ,k  )) * THREE
                   sum_comps = abs(ss0) + &
                     abs(sp(1,i-1,j-1,k-1) + sp(1,i-1,j-1,k  )    &
                        +sp(1,i-1,j  ,k-1) + sp(1,i-1,j  ,k  )) + &
                     abs(sp(1,i  ,j-1,k-1) + sp(1,i  ,j-1,k  )    &
                        +sp(1,i  ,j  ,k-1) + sp(1,i  ,j  ,k  )) + &
                     abs(sp(1,i-1,j-1,k-1) + sp(1,i-1,j-1,k  )    &
                        +sp(1,i  ,j-1,k-1) + sp(1,i  ,j-1,k  )) + &
                     abs(sp(1,i-1,j  ,k-1) + sp(1,i-1,j  ,k  )    &
                        +sp(1,i  ,j  ,k-1) + sp(1,i  ,j  ,k  )) + &
                     abs(sp(1,i-1,j-1,k-1) + sp(1,i-1,j  ,k-1)    &
                        +sp(1,i  ,j-1,k-1) + sp(1,i  ,j  ,k-1)) + &
                     abs(sp(1,i-1,j-1,k  ) + sp(1,i-1,j  ,k  )    &
                        +sp(1,i  ,j-1,k  ) + sp(1,i  ,j  ,k  )) 
                   r1 = max(r1,sum_comps)
                end if ! end mask test
             end do ! end i
             end do ! end j
             end do ! end k

          end if ! (dm.eq.3)

       end do ! end nfabs
       !$OMP END PARALLEL DO

      else   ! .not. present(mask)

       !$OMP PARALLEL DO PRIVATE(i,j,k,bx,sum_comps,ss0,sp) REDUCTION(max:r1)
       do b = 1, nfabs(sg)
          bx = get_box(sg,b)

          sp => dataptr(sg, b)

          if (dm.eq.2) then
             do j = bx%lo(2), bx%hi(2)+1
             do i = bx%lo(1), bx%hi(1)+1
                ss0 = -(sp(1,i-1,j-1,1)+sp(1,i,j-1,1)+sp(1,i-1,j,1)+sp(1,i,j,1))
                sum_comps = abs(ss0) + HALF * ( &
                   abs(sp(1,i  ,j-1,1)+sp(1,i  ,j  ,1)) + &
                   abs(sp(1,i-1,j-1,1)+sp(1,i-1,j  ,1)) + & 
                   abs(sp(1,i-1,j  ,1)+sp(1,i  ,j  ,1)) + &
                   abs(sp(1,i-1,j-1,1)+sp(1,i  ,j-1,1)) )
                r1 = max(r1,sum_comps)
             end do
             end do

          else if (dm.eq.3) then

             do k = bx%lo(3), bx%hi(3)+1
             do j = bx%lo(2), bx%hi(2)+1
             do i = bx%lo(1), bx%hi(1)+1
                ss0 = -( sp(1,i-1,j-1,k-1) + sp(1,i-1,j  ,k-1) &
                        +sp(1,i  ,j-1,k-1) + sp(1,i  ,j  ,k-1) &
                        +sp(1,i-1,j-1,k  ) + sp(1,i-1,j  ,k  ) &
                        +sp(1,i  ,j-1,k  ) + sp(1,i  ,j  ,k  )) * THREE
                sum_comps = abs(ss0) + &
                  abs(sp(1,i-1,j-1,k-1) + sp(1,i-1,j-1,k  )    &
                     +sp(1,i-1,j  ,k-1) + sp(1,i-1,j  ,k  )) + &
                  abs(sp(1,i  ,j-1,k-1) + sp(1,i  ,j-1,k  )    &
                     +sp(1,i  ,j  ,k-1) + sp(1,i  ,j  ,k  )) + &
                  abs(sp(1,i-1,j-1,k-1) + sp(1,i-1,j-1,k  )    &
                     +sp(1,i  ,j-1,k-1) + sp(1,i  ,j-1,k  )) + &
                  abs(sp(1,i-1,j  ,k-1) + sp(1,i-1,j  ,k  )    &
                     +sp(1,i  ,j  ,k-1) + sp(1,i  ,j  ,k  )) + &
                  abs(sp(1,i-1,j-1,k-1) + sp(1,i-1,j  ,k-1)    &
                     +sp(1,i  ,j-1,k-1) + sp(1,i  ,j  ,k-1)) + &
                  abs(sp(1,i-1,j-1,k  ) + sp(1,i-1,j  ,k  )    &
                     +sp(1,i  ,j-1,k  ) + sp(1,i  ,j  ,k  )) 
                r1 = max(r1,sum_comps)
             end do ! end i
             end do ! end j
             end do ! end k
          end if ! (dm.eq.3)
       end do ! end nfabs
       !$OMP END PARALLEL DO

      end if ! present(mask) test

    else if (stencil_type .eq. ND_DENSE_STENCIL) then

      fx     = ONE/36._dp_t
      fy     = fx
      fz     = fx
      f0     = FOUR * (fx + fy + fz)
      fxyz   = (fx+fy+fz)
      f2y2zx = (TWO*fy+TWO*fz-fx)
      f2x2zy = (TWO*fx+TWO*fz-fy)
      f2x2yz = (TWO*fx+TWO*fy-fz)

      if (present(mask)) then
      
       !$OMP PARALLEL DO PRIVATE(i,j,k,bx,sum_comps,ss0,sp,lp) REDUCTION(max:r1)
       do b = 1, nfabs(sg)
          bx = get_box(sg,b)

          sp => dataptr(sg, b)
          lp => dataptr(mask, b)

          if (dm.eq.2) then
             do j = bx%lo(2), bx%hi(2)+1
             do i = bx%lo(1), bx%hi(1)+1
                if ( lp(i,j,1,1) ) then
                  ss0 = -TWO * THIRD * (sp(1,i-1,j-1,1) + sp(1,i,j-1,1) + &
                                        sp(1,i-1,j  ,1) + sp(1,i,j  ,1))
                  sum_comps = abs(ss0) + THIRD * ( &
                         abs(sp(1,i-1,j-1,1)) + abs(sp(1,i  ,j-1,1)) + &
                         abs(sp(1,i-1,j  ,1)) + abs(sp(1,i  ,j  ,1)) + &
                      HALF * ( abs(sp(1,i-1,j-1,1) + sp(1,i  ,j-1,1)) &
                              +abs(sp(1,i-1,j-1,1) + sp(1,i-1,j  ,1)) &
                              +abs(sp(1,i  ,j-1,1) + sp(1,i  ,j  ,1)) &
                              +abs(sp(1,i-1,j  ,1) + sp(1,i  ,j  ,1)) ) )
                  r1 = max(r1,sum_comps)
                end if ! mask test
             end do
             end do

          else if (dm.eq.3) then

             do k = bx%lo(3), bx%hi(3)+1
             do j = bx%lo(2), bx%hi(2)+1
             do i = bx%lo(1), bx%hi(1)+1
                if ( lp(i,j,k,1) ) then
                  ss0 =  -( sp(1,i-1,j-1,k-1) + sp(1,i,j-1,k-1) &
                           +sp(1,i-1,j  ,k-1) + sp(1,i,j  ,k-1) &
                           +sp(1,i-1,j-1,k  ) + sp(1,i,j-1,k  ) &
                           +sp(1,i-1,j  ,k  ) + sp(1,i,j  ,k  ) ) * f0
                  sum_comps = abs(ss0) + &
                        fxyz * ( &   ! Corners
                    abs(sp(1,i-1,j-1,k-1)) + abs(sp(1,i  ,j-1,k-1)) + &
                    abs(sp(1,i-1,j  ,k-1)) + abs(sp(1,i  ,j  ,k-1)) + &
                    abs(sp(1,i-1,j-1,k  )) + abs(sp(1,i  ,j-1,k  )) + &
                    abs(sp(1,i-1,j  ,k  )) + abs(sp(1,i  ,j  ,k  )) ) &
                      + f2y2zx * ( & ! Edges in x-direction
                    abs(sp(1,i  ,j-1,k-1) + sp(1,i-1,j-1,k-1)) + & 
                    abs(sp(1,i  ,j  ,k-1) + sp(1,i-1,j  ,k-1)) + & 
                    abs(sp(1,i  ,j-1,k  ) + sp(1,i-1,j-1,k  )) + & 
                    abs(sp(1,i  ,j  ,k  ) + sp(1,i-1,j  ,k  )) ) & 
                      + f2x2zy * ( & ! Edges in y-direction
                    abs(sp(1,i-1,j-1,k-1) + sp(1,i-1,j  ,k-1)) + &
                    abs(sp(1,i  ,j-1,k-1) + sp(1,i  ,j  ,k-1)) + &
                    abs(sp(1,i-1,j-1,k  ) + sp(1,i-1,j  ,k  )) + &
                    abs(sp(1,i  ,j-1,k  ) + sp(1,i  ,j  ,k  )) ) &
                      + f2x2yz * ( & ! Edges in z-direction
                    abs(sp(1,i-1,j-1,k-1) + sp(1,i-1,j-1,k  )) + &
                    abs(sp(1,i  ,j-1,k-1) + sp(1,i  ,j-1,k  )) + &
                    abs(sp(1,i-1,j  ,k-1) + sp(1,i-1,j  ,k  )) + &
                    abs(sp(1,i  ,j  ,k-1) + sp(1,i  ,j  ,k  )) )

                  if (.not. uniform_dh) then
                  !
                  ! Add faces (only non-zero for non-uniform dx)
                  !
                  sum_comps = abs(ss0) + &
                    abs( (FOUR*fx-TWO*fy-TWO*fz)*(sp(1,i-1,j-1,k-1) + &
                     sp(1,i-1,j-1,k  ) + sp(1,i-1,j  ,k-1) + sp(1,i-1,j,k)) ) + &
                    abs( (FOUR*fx-TWO*fy-TWO*fz)*(sp(1,i  ,j-1,k-1) + &
                     sp(1,i  ,j-1,k  ) + sp(1,i  ,j  ,k-1) + sp(1,i,j,k)) ) + &
                    abs( (FOUR*fy-TWO*fx-TWO*fz)*(sp(1,i-1,j-1,k-1) + &
                     sp(1,i-1,j-1,k  ) + sp(1,i  ,j-1,k-1) + sp(1,i,j-1,k)) ) + &
                    abs( (FOUR*fy-TWO*fx-TWO*fz)*(sp(1,i-1,j  ,k-1) + &
                     sp(1,i-1,j  ,k  ) + sp(1,i  ,j  ,k-1) + sp(1,i,j,k)) ) + &
                    abs( (FOUR*fz-TWO*fx-TWO*fy)*(sp(1,i-1,j-1,k-1) + &
                     sp(1,i-1,j  ,k-1)  + sp(1,i  ,j-1,k-1) + sp(1,i,j,k-1)) ) + &
                    abs( (FOUR*fz-TWO*fx-TWO*fy)*(sp(1,i-1,j-1,k  ) + &
                     sp(1,i-1,j  ,k  ) + sp(1,i  ,j-1,k  ) + sp(1,i,j,k)) ) 
                  end if  ! end .not. uniform_dh 
               r1 = max(r1,sum_comps)
             end if  ! end mask test
          end do ! end i
          end do ! end j
          end do ! end k
          end if ! (dm.eq.3)
       end do ! end nfabs
       !$OMP END PARALLEL DO

      else   ! present(mask) test

       !$OMP PARALLEL DO PRIVATE(i,j,k,bx,sum_comps,ss0,sp) REDUCTION(max:r1)
       do b = 1, nfabs(sg)
          ! This is the cell-centered box with one ghost cell
          bx = get_box(sg,b)
          sp => dataptr(sg, b)

          if (dm.eq.2) then

             do j = bx%lo(2), bx%hi(2)+1
             do i = bx%lo(1), bx%hi(1)+1
               ss0 = -TWO * THIRD * (sp(1,i-1,j-1,1) + sp(1,i,j-1,1) + &
                                     sp(1,i-1,j  ,1) + sp(1,i,j  ,1))
               sum_comps = abs(ss0) + THIRD * ( &
                      abs(sp(1,i-1,j-1,1)) + abs(sp(1,i  ,j-1,1)) + &
                      abs(sp(1,i-1,j  ,1)) + abs(sp(1,i  ,j  ,1)) + &
                   HALF * ( abs(sp(1,i-1,j-1,1) + sp(1,i  ,j-1,1)) &
                           +abs(sp(1,i-1,j-1,1) + sp(1,i-1,j  ,1)) &
                           +abs(sp(1,i  ,j-1,1) + sp(1,i  ,j  ,1)) &
                           +abs(sp(1,i-1,j  ,1) + sp(1,i  ,j  ,1)) ) )
               r1 = max(r1,sum_comps)
             end do ! end i
             end do ! end j

          else if (dm.eq.3) then

             do k = bx%lo(3), bx%hi(3)+1
             do j = bx%lo(2), bx%hi(2)+1
             do i = bx%lo(1), bx%hi(1)+1
               ss0 =  -( sp(1,i-1,j-1,k-1) + sp(1,i,j-1,k-1) &
                        +sp(1,i-1,j  ,k-1) + sp(1,i,j  ,k-1) &
                        +sp(1,i-1,j-1,k  ) + sp(1,i,j-1,k  ) &
                        +sp(1,i-1,j  ,k  ) + sp(1,i,j  ,k  ) ) * f0
               sum_comps = abs(ss0) + &
                     fxyz * ( &   ! Corners
                 abs(sp(1,i-1,j-1,k-1)) + abs(sp(1,i  ,j-1,k-1)) + &
                 abs(sp(1,i-1,j  ,k-1)) + abs(sp(1,i  ,j  ,k-1)) + &
                 abs(sp(1,i-1,j-1,k  )) + abs(sp(1,i  ,j-1,k  )) + &
                 abs(sp(1,i-1,j  ,k  )) + abs(sp(1,i  ,j  ,k  )) ) &
                   + f2y2zx * ( & ! Edges in x-direction
                 abs(sp(1,i  ,j-1,k-1) + sp(1,i-1,j-1,k-1)) + & 
                 abs(sp(1,i  ,j  ,k-1) + sp(1,i-1,j  ,k-1)) + & 
                 abs(sp(1,i  ,j-1,k  ) + sp(1,i-1,j-1,k  )) + & 
                 abs(sp(1,i  ,j  ,k  ) + sp(1,i-1,j  ,k  )) ) & 
                   + f2x2zy * ( & ! Edges in y-direction
                 abs(sp(1,i-1,j-1,k-1) + sp(1,i-1,j  ,k-1)) + &
                 abs(sp(1,i  ,j-1,k-1) + sp(1,i  ,j  ,k-1)) + &
                 abs(sp(1,i-1,j-1,k  ) + sp(1,i-1,j  ,k  )) + &
                 abs(sp(1,i  ,j-1,k  ) + sp(1,i  ,j  ,k  )) ) &
                   + f2x2yz * ( & ! Edges in z-direction
                 abs(sp(1,i-1,j-1,k-1) + sp(1,i-1,j-1,k  )) + &
                 abs(sp(1,i  ,j-1,k-1) + sp(1,i  ,j-1,k  )) + &
                 abs(sp(1,i-1,j  ,k-1) + sp(1,i-1,j  ,k  )) + &
                 abs(sp(1,i  ,j  ,k-1) + sp(1,i  ,j  ,k  )) )

               if (.not. uniform_dh) then
               !
               ! Add faces (only non-zero for non-uniform dx)
               !
               sum_comps = abs(ss0) + &
                    abs( (FOUR*fx-TWO*fy-TWO*fz)*(sp(1,i-1,j-1,k-1) + sp(1,i-1,j-1,k  ) &
                     +sp(1,i-1,j  ,k-1) + sp(1,i-1,j  ,k  )) ) + &
                    abs( (FOUR*fx-TWO*fy-TWO*fz)*(sp(1,i  ,j-1,k-1) + sp(1,i  ,j-1,k  ) &
                     +sp(1,i  ,j  ,k-1) + sp(1,i  ,j  ,k  )) ) + &
                    abs( (FOUR*fy-TWO*fx-TWO*fz)*(sp(1,i-1,j-1,k-1) + sp(1,i-1,j-1,k  ) &
                     +sp(1,i  ,j-1,k-1) + sp(1,i  ,j-1,k  )) ) + &
                    abs( (FOUR*fy-TWO*fx-TWO*fz)*(sp(1,i-1,j  ,k-1) + sp(1,i-1,j  ,k  ) &
                     +sp(1,i  ,j  ,k-1) + sp(1,i  ,j  ,k  )) ) + &
                    abs( (FOUR*fz-TWO*fx-TWO*fy)*(sp(1,i-1,j-1,k-1) + sp(1,i-1,j  ,k-1) &
                     +sp(1,i  ,j-1,k-1) + sp(1,i  ,j  ,k-1)) ) + &
                    abs( (FOUR*fz-TWO*fx-TWO*fy)*(sp(1,i-1,j-1,k  ) + sp(1,i-1,j  ,k  ) &
                     +sp(1,i  ,j-1,k  ) + sp(1,i  ,j  ,k  )) ) 
               end if  ! end .not. uniform_dh 
            r1 = max(r1,sum_comps)
          end do
          end do
          end do
          end if
       end do
       !$OMP END PARALLEL DO

      end if   ! present(mask) test

    else if (stencil_type .eq. ND_VATER_STENCIL) then

      if (dm .eq. 3) &
          call bl_error("nodal_stencil_norm: ND_VATER_STENCIL not implemented in 3-d")

      if (present(mask)) then
      
       !$OMP PARALLEL DO PRIVATE(i,j,k,bx,sum_comps,ss0,sp,lp) REDUCTION(max:r1)
       do b = 1, nfabs(sg)
          bx = get_box(sg,b)

          sp => dataptr(sg, b)
          lp => dataptr(mask, b)

          if (dm.eq.2) then
             do j = bx%lo(2), bx%hi(2)+1
             do i = bx%lo(1), bx%hi(1)+1
                if ( lp(i,j,1,1) ) then
                  ss0 = -THREE4TH * (sp(1,i-1,j-1,1) + sp(1,i,j-1,1) + &
                                     sp(1,i-1,j  ,1) + sp(1,i,j  ,1))
                  sum_comps = abs(ss0) + FOURTH * ( &
                         abs(sp(1,i-1,j-1,1)) + abs(sp(1,i  ,j-1,1)) + &
                         abs(sp(1,i-1,j  ,1)) + abs(sp(1,i  ,j  ,1)) + &
                         abs(sp(1,i-1,j-1,1) + sp(1,i  ,j-1,1)) + &
                         abs(sp(1,i-1,j-1,1) + sp(1,i-1,j  ,1)) + &
                         abs(sp(1,i  ,j-1,1) + sp(1,i  ,j  ,1)) + &
                         abs(sp(1,i-1,j  ,1) + sp(1,i  ,j  ,1)) )
                  r1 = max(r1,sum_comps)
                end if ! mask test
             end do
             end do
          end if

       end do ! end nfabs
       !$OMP END PARALLEL DO

      else   ! present(mask) test

       !$OMP PARALLEL DO PRIVATE(i,j,k,bx,sum_comps,ss0,sp) REDUCTION(max:r1)
       do b = 1, nfabs(sg)
          ! This is the cell-centered box with one ghost cell
          bx = get_box(sg,b)
          sp => dataptr(sg, b)

          if (dm.eq.2) then

             do j = bx%lo(2), bx%hi(2)+1
             do i = bx%lo(1), bx%hi(1)+1
               ss0 = -THREE4TH * (sp(1,i-1,j-1,1) + sp(1,i,j-1,1) + &
                                  sp(1,i-1,j  ,1) + sp(1,i,j  ,1))
               sum_comps = abs(ss0) + FOURTH * ( &
                      abs(sp(1,i-1,j-1,1)) + abs(sp(1,i  ,j-1,1)) + &
                      abs(sp(1,i-1,j  ,1)) + abs(sp(1,i  ,j  ,1)) + &
                      abs(sp(1,i-1,j-1,1) + sp(1,i  ,j-1,1)) + &
                      abs(sp(1,i-1,j-1,1) + sp(1,i-1,j  ,1)) + &
                      abs(sp(1,i  ,j-1,1) + sp(1,i  ,j  ,1)) + &
                      abs(sp(1,i-1,j  ,1) + sp(1,i  ,j  ,1)) )
               r1 = max(r1,sum_comps)
             end do ! end i
             end do ! end j

          end if
       end do
       !$OMP END PARALLEL DO

      end if   ! present(mask) test

    else 
       call bl_error("Unknown stencil_type in nodal_stencil_norm")
    end if

    r = r1

    if ( .not. llocal ) call parallel_reduce(r, r1, MPI_MAX)

    call destroy(bpt)

  end function nodal_stencil_norm

end module ml_nd_module
