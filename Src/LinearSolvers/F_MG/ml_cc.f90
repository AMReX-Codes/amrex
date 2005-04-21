module ml_cc_module

  use BoxLib
  use bl_constants_module
  use omp_module
  use f2kcli
  use stencil_module
  use mg_module
  use list_box_module
  use ml_boxarray_module
  use itsol_module
  use box_util_module
  use bl_IO_module

  use ml_restriction_module
  use ml_prolongation_module
  use ml_interface_stencil_module
  use ml_util_module
  use bndry_reg_module

contains

subroutine ml_cc(la_tower,mgt,rh,full_soln,fine_mask,ref_ratio,do_diagnostics,eps)

  implicit none

  type(layout), intent(in)       :: la_tower(:)
  type(mg_tower) , intent(inout) :: mgt(:)
  type( multifab), intent(inout) :: rh(:)
  type( multifab), intent(inout) :: full_soln(:)
  type(lmultifab), intent(in   ) :: fine_mask(:)
  integer        , intent(in   ) :: ref_ratio(:,:)
  integer        , intent(in   ) :: do_diagnostics 
  real(dp_t)     , intent(in   ) :: eps

  integer :: nlevs
  type(multifab), pointer  ::      soln(:) => Null()
  type(multifab), pointer  ::        uu(:) => Null()
  type(multifab), pointer  ::   uu_hold(:) => Null()
  type(multifab), pointer  ::       res(:) => Null()
  type(multifab), pointer  ::  temp_res(:) => Null()

  type(bndry_reg), pointer :: brs_flx(:) => Null()
  type(bndry_reg), pointer :: brs_bcs(:) => Null()

  type(box   ) :: pd,pdc
  type(layout) :: la
  integer :: i, n, dm
  integer :: mglev, mglev_crse, iter, it
  logical :: fine_converged

  real(dp_t) :: Anorm, bnorm, res_norm
  real(dp_t) :: snrm(2)

  dm = rh(1)%dim

  nlevs = size(la_tower)

  allocate(soln(nlevs), uu(nlevs), res(nlevs), temp_res(nlevs))
  allocate(uu_hold(2:nlevs-1))
  allocate(brs_flx(2:nlevs))
  allocate(brs_bcs(2:nlevs))

  do n = 2,nlevs-1
     la = la_tower(n)
     call multifab_build(uu_hold(n),la,1,1)
     call setval( uu_hold(n), ZERO,all=.true.)
  end do

  do n = nlevs, 1, -1

     la = la_tower(n)
     call multifab_build(    soln(n), la, 1, 1)
     call multifab_build(      uu(n), la, 1, 1)
     call multifab_build(     res(n), la, 1, 0)
     call multifab_build(temp_res(n), la, 1, 0)
     call setval(    soln(n), ZERO,all=.true.)
     call setval(      uu(n), ZERO,all=.true.)
     call setval(     res(n), ZERO,all=.true.)
     call setval(temp_res(n), ZERO,all=.true.)
 
     if ( n == 1 ) exit

     ! Build the (coarse resolution) flux registers to be used in computing
     !  the residual at a non-finest AMR level.

     pdc = layout_get_pd(la_tower(n-1))
     call bndry_reg_build(brs_flx(n), la, ref_ratio(n-1,:), pdc, width = 0)
     call bndry_reg_build(brs_bcs(n), la, ref_ratio(n-1,:), pdc, width = 2)

  end do

  do n = nlevs,2,-1
    mglev      = mgt(n  )%nlevels
    mglev_crse = mgt(n-1)%nlevels
    call ml_restriction(rh(n-1), rh(n), mgt(n)%mm(mglev),&
                        mgt(n-1)%mm(mglev_crse), mgt(n)%face_type, ref_ratio(n-1,:))
  end do
  bnorm = ml_norm_inf(rh,fine_mask)

  Anorm = stencil_norm(mgt(nlevs)%ss(mgt(nlevs)%nlevels))
  do n = 1, nlevs-1
     Anorm = max(stencil_norm(mgt(n)%ss(mgt(n)%nlevels), fine_mask(n)), Anorm)
  end do

  do n = 1,nlevs,1
    mglev = mgt(n)%nlevels
    call mg_defect(mgt(n)%ss(mglev),res(n),rh(n),full_soln(n),mgt(n)%mm(mglev))
  end do

  do n = nlevs,2,-1
    mglev      = mgt(n  )%nlevels
    mglev_crse = mgt(n-1)%nlevels

    pdc = layout_get_pd(la_tower(n-1))
    call crse_fine_residual_cc(n,mgt,full_soln,res(n-1),brs_flx(n),pdc,ref_ratio(n-1,:))

    call ml_restriction(res(n-1), res(n), mgt(n)%mm(mglev),&
                        mgt(n-1)%mm(mglev_crse), mgt(n)%face_type, ref_ratio(n-1,:))
  enddo

  do n = 1,nlevs
     call multifab_copy(rh(n),res(n),all=.true.)
  end do

! ****************************************************************************

  fine_converged = .false.

  do iter = 1, mgt(nlevs)%max_iter

     if ( (iter .eq. 1) .or. fine_converged ) then
       if ( ml_converged(res, soln, fine_mask, bnorm, Anorm, eps) ) exit
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

        if (do_diagnostics == 1 .and. parallel_ioprocessor() ) &
           print *,'DWN: RES BEFORE GSRB AT LEVEL ',n, norm_inf(res(n))

        ! Relax ...
        if (iter < mgt(nlevs)%max_iter) &
        call mg_tower_cycle(mgt(n), mgt(n)%cycle, mglev, mgt(n)%ss(mglev), &
                            uu(n), res(n), mgt(n)%mm(mglev), mgt(n)%nu1, mgt(n)%nu2, &
                            mgt(n)%gamma)

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
           call multifab_copy(res(n), temp_res(n), all=.true.)

           if (do_diagnostics == 1 .and. parallel_ioprocessor()) &
             print *,'DWN: RES AFTER  GSRB AT LEVEL ',n, norm_inf(res(n))

           ! Compute CRSE-FINE Res = Res - Lap(uu)
           pdc = layout_get_pd(la_tower(n-1))
           call crse_fine_residual_cc(n,mgt,uu,res(n-1),brs_flx(n),pdc,ref_ratio(n-1,:))

           ! Restrict FINE Res to COARSE Res (important to do this last
           !     so we overwrite anything extra which may have been defined
           !     above near fine-fine interfaces)
           call ml_restriction(res(n-1), res(n), mgt(n)%mm(mglev), &
                               mgt(n-1)%mm(mglev_crse), mgt(n)%face_type, ref_ratio(n-1,:))

           ! Copy u_hold = uu
           if (n < nlevs) call multifab_copy(uu_hold(n), uu(n), all=.true.)

           ! Set: uu = 0
           call setval(uu(n), ZERO, all=.true.)

        else

           if (do_diagnostics == 1 .and. parallel_ioprocessor()) then
             call mg_defect(mgt(n)%ss(mglev), temp_res(n), res(n), uu(n), mgt(n)%mm(mglev))
             print *,'DWN: RES AFTER  GSRB AT LEVEL ',n, norm_inf(temp_res(n))
           end if

        end if

     end do

     !   Back up the V-cycle
     do n = 2, nlevs

        pd = layout_get_pd(la_tower(n))
        mglev = mgt(n)%nlevels

        ! Interpolate uu from coarser level
        call ml_prolongation(uu(n), uu(n-1), pd, ref_ratio(n-1,:))

        ! Add: soln += uu
        call saxpy(soln(n), ONE, uu(n), .true.)

        ! Add: uu_hold += uu
        if (n < nlevs) call saxpy(uu_hold(n), ONE, uu(n), .true.)

        ! Interpolate uu to supply boundary conditions for new residual calculation
        call bndry_reg_copy(brs_bcs(n), uu(n-1))
        do i = 1, dm
           call ml_interp_bcs(uu(n), brs_bcs(n)%bmf(i,0), pd, ref_ratio(n-1,:), -i)
           call ml_interp_bcs(uu(n), brs_bcs(n)%bmf(i,1), pd, ref_ratio(n-1,:), +i)
        end do

        ! Compute Res = Res - Lap(uu)
        call mg_defect(mgt(n)%ss(mglev), temp_res(n), res(n), uu(n), mgt(n)%mm(mglev))
        call multifab_copy(res(n), temp_res(n), all=.true.)

        if (do_diagnostics == 1 .and. parallel_ioprocessor() ) &
           print *,'UP : RES BEFORE GSRB AT LEVEL ',n, norm_inf(temp_res(n))

        ! Set: uu = 0
        call setval(uu(n), ZERO, all=.true.)

        ! Relax ...
        call mg_tower_cycle(mgt(n), mgt(n)%cycle, mglev, mgt(n)%ss(mglev), &
                            uu(n), res(n), mgt(n)%mm(mglev), mgt(n)%nu1, mgt(n)%nu2, &
                            mgt(n)%gamma)

        ! Compute Res = Res - Lap(uu)
        call mg_defect(mgt(n)%ss(mglev), temp_res(n), res(n), uu(n), mgt(n)%mm(mglev))
        call multifab_copy(res(n), temp_res(n), all=.true.)

        if (do_diagnostics == 1 .and. parallel_ioprocessor() ) then
           print *,'UP : RES AFTER  GSRB AT LEVEL ',n, norm_inf(res(n))
           if (n == nlevs) print *,' '
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
       call multifab_fill_boundary(soln(n))
     end do 

     ! Interpolate soln to supply boundary conditions 
     do n = 2,nlevs
       pd = layout_get_pd(la_tower(n))
       call bndry_reg_copy(brs_bcs(n), soln(n-1))
       do i = 1, dm
          call ml_interp_bcs(soln(n), brs_bcs(n)%bmf(i,0), pd, ref_ratio(n-1,:), -i)
          call ml_interp_bcs(soln(n), brs_bcs(n)%bmf(i,1), pd, ref_ratio(n-1,:), +i)
       end do
     end do

!    Optimization so don't have to do multilevel convergence test each time

!    Compute the residual on just the finest level
     n = nlevs
     mglev = mgt(n)%nlevels
     call mg_defect(mgt(n)%ss(mglev),res(n),rh(n),soln(n),mgt(n)%mm(mglev))

     if ( ml_fine_converged(res, soln, bnorm, Anorm, eps) ) then

       fine_converged = .true.

!      Compute the residual on every level
       do n = 1,nlevs-1
          mglev = mgt(n)%nlevels
          call mg_defect(mgt(n)%ss(mglev),res(n),rh(n),soln(n),mgt(n)%mm(mglev))
       end do
!      Compute the coarse-fine residual 
       do n = nlevs,2,-1
          pdc = layout_get_pd(la_tower(n-1))
          call crse_fine_residual_cc(n,mgt,soln,res(n-1),brs_flx(n),pdc,ref_ratio(n-1,:))
       end do

!      Average the fine residual onto the coarser level
       do n = nlevs,2,-1
          mglev      = mgt(n  )%nlevels
          mglev_crse = mgt(n-1)%nlevels
          call ml_restriction(res(n-1), res(n), mgt(n)%mm(mglev),&
                              mgt(n-1)%mm(mglev_crse), mgt(n)%face_type, ref_ratio(n-1,:))
       end do
       if ( mgt(nlevs)%verbose > 0 .and. parallel_IOProcessor() ) then
         do n = 1,nlevs
           write(unit=*, fmt='(i3,": Level ",i2,"  : SL_Ninf(defect) = ",g15.8)') iter,n,norm_inf(res(n))
         end do
         write(unit=*, fmt='(i3,": All Levels: ML_Ninf(defect) = ",g15.8)') iter, ml_norm_inf(res,fine_mask)
       end if

     else 

        fine_converged = .false.
        if ( mgt(nlevs)%verbose > 0 .and. parallel_IOProcessor() ) &
          write(unit=*, fmt='(i3,": FINE_Ninf(defect) = ",g15.8)') iter, norm_inf(res(nlevs))

     end if


  end do

! ****************************************************************************

  if ( mgt(nlevs)%verbose > 0 .AND. parallel_IOProcessor() ) &
     write(unit=*, fmt='("MG finished at ", i3, " iterations")') iter-1

  ! Add: soln += full_soln
  do n = 1,nlevs
     call saxpy(full_soln(n),ONE,soln(n))
  end do

  deallocate(soln, uu, uu_hold, res, temp_res)

contains

  subroutine crse_fine_residual_cc(n,mgt,uu,crse_res,brs_flx,pdc,ref_ratio)

     integer        , intent(in   ) :: n
     type(mg_tower) , intent(inout) :: mgt(:)
     type(bndry_reg), intent(inout) :: brs_flx
     type(multifab) , intent(inout) :: uu(:)
     type(multifab) , intent(inout) :: crse_res
     type(box)      , intent(in   ) :: pdc
     integer        , intent(in   ) :: ref_ratio(:)

     integer :: i,dm,mglev

     dm = brs_flx%dim
     mglev = mgt(n)%nlevels

     do i = 1, dm
        call ml_fill_fluxes(mgt(n)%ss(mglev), brs_flx%bmf(i,0), &
             uu(n), mgt(n)%mm(mglev), ref_ratio(i), -1, i)
        call ml_interface(crse_res, brs_flx%bmf(i,0), uu(n-1), &
             mgt(n-1)%ss(mgt(n-1)%nlevels), pdc, -1, i, ONE)

        call ml_fill_fluxes(mgt(n)%ss(mglev), brs_flx%bmf(i,1), &
             uu(n), mgt(n)%mm(mglev), ref_ratio(i), 1, i)
        call ml_interface(crse_res, brs_flx%bmf(i,1), uu(n-1), &
             mgt(n-1)%ss(mgt(n-1)%nlevels), pdc, +1, i, ONE)
     end do

  end subroutine crse_fine_residual_cc

  function ml_fine_converged(res, sol, bnorm, Anorm, eps) result(r)
    logical :: r
    type(multifab), intent(in) :: res(:), sol(:)
    real(dp_t), intent(in) :: Anorm, eps, bnorm
    real(dp_t) :: ni_res, ni_sol
    integer    :: nlevs
    nlevs = size(res)
    ni_res = norm_inf(res(nlevs))
    ni_sol = norm_inf(sol(nlevs))
    r =  ni_res <= eps*(Anorm*ni_sol + bnorm) .or. &
         ni_res <= spacing(Anorm)
  end function ml_fine_converged

  function ml_converged(res, sol, mask, bnorm, Anorm, eps) result(r)
    logical :: r
    type(multifab), intent(in) :: res(:), sol(:)
    type(lmultifab), intent(in) :: mask(:)
    real(dp_t), intent(in) :: Anorm, eps, bnorm
    real(dp_t) :: ni_res, ni_sol
    ni_res = ml_norm_inf(res, mask)
    ni_sol = ml_norm_inf(sol, mask)
    r =  ni_res <= eps*(Anorm*ni_sol + bnorm) .or. &
         ni_res <= spacing(Anorm)
  end function ml_converged

  function ml_norm_inf(rr, mask) result(r)
    real(dp_t)  :: r
    type(multifab), intent(in) :: rr(:)
    type(lmultifab), intent(in) :: mask(:)
    integer n
    r = 0
    do n = 1, size(rr)
       r = max(norm_inf(rr(n),mask(n)), r)
    end do
  end function ml_norm_inf

end subroutine ml_cc

end module ml_cc_module
