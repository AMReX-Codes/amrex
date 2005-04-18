module ml_nd_module

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

subroutine ml_nd(la_tower,mgt,rh,full_soln,fine_mask,ref_ratio,do_diagnostics,eps)

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

  type(box   ) :: pd,pdc
  type(layout) :: la
  integer :: i, n, dm
  integer :: mglev, mglev_crse, iter, it

  real(dp_t) :: Anorm, bnorm, res_norm
  real(dp_t) :: snrm(2)

  logical :: all_done
  logical, allocatable :: nodal(:)

  dm = rh(1)%dim
  allocate(nodal(dm))
  nodal = .True.

  nlevs = size(la_tower)

  allocate(soln(nlevs), uu(nlevs), uu_hold(nlevs), res(nlevs))
  allocate(temp_res(nlevs))
  allocate(brs_flx(2:nlevs))

  do n = nlevs, 1, -1

     la = la_tower(n)
     call multifab_build(    soln(n), la, 1, 1, nodal)
     call multifab_build(      uu(n), la, 1, 1, nodal)
     call multifab_build( uu_hold(n), la, 1, 1, nodal)
     call multifab_build(     res(n), la, 1, 1, nodal)
     call multifab_build(temp_res(n), la, 1, 1, nodal)
     call setval(    soln(n), ZERO,all=.true.)
     call setval(      uu(n), ZERO,all=.true.)
     call setval( uu_hold(n), ZERO,all=.true.)
     call setval(     res(n), ZERO,all=.true.)
     call setval(temp_res(n), ZERO,all=.true.)
 
     if ( n == 1 ) exit

     ! Build the (coarse resolution) flux registers to be used in computing
     !  the residual at a non-finest AMR level.

     pdc = layout_get_pd(la_tower(n-1))
     call bndry_reg_build(brs_flx(n), la, ref_ratio(n-1,:), pdc, nodal = nodal)

  end do

! ****************************************************************************

  bnorm = ml_norm_inf(rh,fine_mask)

  Anorm = stencil_norm(mgt(nlevs)%ss(mgt(nlevs)%nlevels))
  do n = 1, nlevs-1
     Anorm = max(stencil_norm(mgt(n)%ss(mgt(n)%nlevels), fine_mask(n)), Anorm)
  end do

  do n = nlevs,1,-1
    mglev = mgt(n)%nlevels
    call mg_defect(mgt(n)%ss(mglev),res(n),rh(n),full_soln(n),mgt(n)%mm(mglev))
  end do

  do n = nlevs,2,-1
    mglev      = mgt(n  )%nlevels
    mglev_crse = mgt(n-1)%nlevels
    call ml_restriction(res(n-1), res(n), mgt(n)%mm(mglev),&
                        mgt(n-1)%mm(mglev_crse), mgt(n)%face_type, ref_ratio(n-1,:))
    pdc = layout_get_pd(la_tower(n-1))
    call crse_fine_residual_nodal(n,mgt,brs_flx(n),res(n-1),temp_res(n), &
                                  full_soln(n-1),full_soln(n),ref_ratio(n-1,:),pdc)
  enddo

  do n = 1,nlevs
     call multifab_copy(rh(n),res(n),all=.true.)
  end do

  do iter = 1, mgt(nlevs)%max_iter

     if ( ml_converged(res, soln, fine_mask, bnorm, Anorm, eps) ) exit

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

        if (do_diagnostics == 1 .and.  parallel_ioprocessor() ) &
           print *,'DWN: RES BEFORE GSRB AT LEVEL ',n, norm_inf(res(n))

        ! Relax ...
        if (n > 1) then
          call mini_cycle(mgt(n), mgt(n)%cycle, mglev, mgt(n)%ss(mglev), &
               uu(n), res(n), mgt(n)%mm(mglev), mgt(n)%nu1, mgt(n)%nu2, &
               mgt(n)%gamma)
        else 
          call mg_tower_cycle(mgt(n), mgt(n)%cycle, mglev, mgt(n)%ss(mglev), &
               uu(n), res(n), mgt(n)%mm(mglev), mgt(n)%nu1, mgt(n)%nu2, &
               mgt(n)%gamma)
        end if

        ! Add: Soln += uu
        call saxpy(soln(n),ONE,uu(n))

        if (n > 1) then
           mglev_crse = mgt(n-1)%nlevels

           ! Compute COARSE Res = Rh - Lap(Soln)
           call mg_defect(mgt(n-1)%ss(mglev_crse),res(n-1), &
                          rh(n-1),soln(n-1),mgt(n-1)%mm(mglev_crse))

           ! Compute FINE Res = Res - Lap(uu)
           mglev = mgt(n)%nlevels
           call mg_defect(mgt(n)%ss(mglev), temp_res(n), &
                          res(n),uu(n),mgt(n)%mm(mglev))
           call multifab_copy(res(n),temp_res(n),all=.true.)

           if (do_diagnostics == 1 .and.  parallel_ioprocessor() ) &
              print *,'DWN: RES AFTER  GSRB AT LEVEL ',n, norm_inf(res(n))

           ! Restrict FINE Res to COARSE Res
           call ml_restriction(res(n-1), res(n), mgt(n)%mm(mglev),& 
                               mgt(n-1)%mm(mglev_crse), mgt(n)%face_type, ref_ratio(n-1,:))
           
           ! Compute CRSE-FINE Res = Rh - Lap(Soln)
           pdc = layout_get_pd(la_tower(n-1))
           call crse_fine_residual_nodal(n,mgt,brs_flx(n),res(n-1),temp_res(n), &
                                         soln(n-1),soln(n),ref_ratio(n-1,:),pdc)

           ! Copy u_hold = uu
           if (n < nlevs) call multifab_copy(uu_hold(n),uu(n),all=.true.)

           ! Set: uu = 0
           call setval(uu(n),ZERO,all=.true.)

        else

           if (do_diagnostics == 1 .and.  parallel_ioprocessor() ) then
              call mg_defect(mgt(n)%ss(mglev),temp_res(n), &
                             res(n),uu(n),mgt(n)%mm(mglev))
              print *,'DWN: RES AFTER  GSRB AT LEVEL ',n, norm_inf(temp_res(n))
           end if

        end if


     end do

     !   Back up the V-cycle
     do n = 2, nlevs

        pd = layout_get_pd(la_tower(n))
        mglev = mgt(n)%nlevels

        ! Interpolate uu from coarser level
        if (iter == 1) call saxpy(uu(n-1),  ONE, full_soln(n-1))
        call ml_prolongation(uu(n), uu(n-1), pd, ref_ratio(n-1,:))
        if (iter == 1) call saxpy(uu(n-1), -ONE, full_soln(n-1))

        ! Subtract: uu -= full_soln
        !     Must do this in order to remove interpolated full_soln...
        if (iter == 1) call saxpy(uu(n),-ONE,full_soln(n))

        ! Add: Soln += uu
        call saxpy(soln(n), ONE, uu(n), .true.)

        ! Add: uu_hold += uu so that it interpolated uu be interpolated too.
        if (n < nlevs) call saxpy(uu_hold(n), ONE, uu(n), .true.)

        ! Compute Res = Res - Lap(uu)
        call mg_defect(mgt(n)%ss(mglev),temp_res(n),res(n),uu(n),mgt(n)%mm(mglev))
        call multifab_copy(res(n),temp_res(n),all=.true.)

        if (do_diagnostics == 1 .and.  parallel_ioprocessor() ) &
           print *,'UP : RES BEFORE GSRB AT LEVEL ',n, norm_inf(res(n))

        ! Set: uu = 0
        call setval(uu(n),ZERO,all=.true.)

        ! Relax ...
        call mini_cycle(mgt(n), mgt(n)%cycle, mglev, mgt(n)%ss(mglev), &
             uu(n), res(n), mgt(n)%mm(mglev), mgt(n)%nu1, mgt(n)%nu2, &
             mgt(n)%gamma)

        ! Compute Res = Res - Lap(uu)
        call mg_defect(mgt(n)%ss(mglev),temp_res(n),res(n),uu(n),mgt(n)%mm(mglev))
        call multifab_copy(res(n),temp_res(n),all=.true.)

        if (do_diagnostics == 1 .and.  parallel_ioprocessor() ) then
           print *,'UP : RES AFTER  GSRB AT LEVEL ',n, norm_inf(res(n))
           if (n == nlevs) print *,' '
        end if

        ! Add: soln += uu
        call saxpy(soln(n), ONE, uu(n), .true.)

        ! Add: uu += uu_hold so that it will be interpolated too.
        if (n < nlevs) call saxpy(  uu(n), ONE, uu_hold(n), .true.)

     end do

!    Inject the solution to the coarser grids.
     do n = nlevs,2,-1
       mglev      = mgt(n)%nlevels
       mglev_crse = mgt(n-1)%nlevels
       call ml_restriction(soln(n-1), soln(n), mgt(n)%mm(mglev), &
                           mgt(n-1)%mm(mglev_crse), mgt(n)%face_type, &
                           ref_ratio(n-1,:), inject = .true.)
     end do 

     do n = 1,nlevs
       call multifab_fill_boundary(soln(n))
     end do 

     do n = 1,nlevs
        mglev = mgt(n)%nlevels
        call mg_defect(mgt(n)%ss(mglev),res(n),rh(n),soln(n),mgt(n)%mm(mglev))
     end do

     do n = nlevs,2,-1
        mglev      = mgt(n  )%nlevels
        mglev_crse = mgt(n-1)%nlevels
        call ml_restriction(res(n-1), res(n), mgt(n)%mm(mglev),&
                            mgt(n-1)%mm(mglev_crse), mgt(n)%face_type, ref_ratio(n-1,:))
     end do

     do n = nlevs,2,-1
        pdc = layout_get_pd(la_tower(n-1))
        call crse_fine_residual_nodal(n,mgt,brs_flx(n),res(n-1),temp_res(n), &
                                      soln(n-1),soln(n),ref_ratio(n-1,:),pdc)
     end do

     if ( mgt(nlevs)%verbose > 0 .and. parallel_IOProcessor() ) &
       write(unit=*, fmt='(i3,": Ninf(defect) = ",g15.8)') iter, ml_norm_inf(res, fine_mask)

  end do

  if ( mgt(nlevs)%verbose > 0 .AND. parallel_IOProcessor() ) &
     write(unit=*, fmt='("MG finished at ", i3, " iterations")') iter-1

! Add: soln += full_soln
  do n = 1,nlevs
     call saxpy(full_soln(n),ONE,soln(n))
  end do


! ****************************************************************************

  deallocate(soln, uu, uu_hold, res, temp_res)

contains

  subroutine crse_fine_residual_nodal(n,mgt,brs_flx,crse_res,temp_res, &
                                      crse_soln,fine_soln,ref_ratio,pdc)

     integer        , intent(in   ) :: n
     type(mg_tower) , intent(inout) :: mgt(:)
     type(bndry_reg), intent(inout) :: brs_flx
     type(multifab) , intent(inout) :: crse_res
     type(multifab) , intent(inout) :: temp_res
     type(multifab) , intent(inout) :: crse_soln
     type(multifab) , intent(inout) :: fine_soln
     integer        , intent(in   ) :: ref_ratio(:)
     type(box)      , intent(in   ) :: pdc

     type(layout)   :: la
     type(multifab) :: dummy_rhs
     integer :: i,dm,mglev_crse,mglev_fine

     mglev_crse = mgt(n-1)%nlevels
     mglev_fine = mgt(n  )%nlevels
     dm = temp_res%dim

     la = multifab_get_layout(temp_res)
     call multifab_build(dummy_rhs, la, 1, 1)
     call setval(dummy_rhs,ZERO,all=.true.)

!    Zero out the flux registers which will hold the fine contributions
     call bndry_reg_setval(brs_flx, ZERO, all = .true.)

!    Compute the fine contributions at faces, edges and corners.

!    First compute a residual which only takes contributions from the
!       grid on which it is calculated.
           call grid_res(mgt(n),mglev_fine,mgt(n)%ss(mglev_fine),temp_res, &
                         dummy_rhs,fine_soln,mgt(n)%mm(mglev_fine),mgt(n)%face_type)

     do i = 1,dm
        call ml_fine_contrib(brs_flx%bmf(i,0), &
                             temp_res,mgt(n)%mm(mglev_fine),ref_ratio,pdc,-i)
        call ml_fine_contrib(brs_flx%bmf(i,1), &
                             temp_res,mgt(n)%mm(mglev_fine),ref_ratio,pdc,+i)
     end do

!    Compute the crse contributions at edges and corners and add to res(n-1).
     do i = 1,dm
        call ml_crse_contrib(crse_res, brs_flx%bmf(i,0), crse_soln, &
             mgt(n-1)%ss(mgt(n-1)%nlevels), &
             mgt(n-1)%mm(mglev_crse), &
             mgt(n  )%mm(mglev_fine), &
             pdc,ref_ratio, -i)
        call ml_crse_contrib(crse_res, brs_flx%bmf(i,1), crse_soln, &
             mgt(n-1)%ss(mgt(n-1)%nlevels), &
             mgt(n-1)%mm(mglev_crse), &
             mgt(n  )%mm(mglev_fine), &
             pdc,ref_ratio, +i)
     end do
  end subroutine crse_fine_residual_nodal

  function ml_converged(res, sol, mask, bnorm, Anorm, eps) result(r)
    logical :: r
    type(multifab), intent(in) :: res(:), sol(:)
    type(lmultifab), intent(in) :: mask(:)
    real(dp_t), intent(in) :: Anorm, eps, bnorm
    real(dp_t) :: ni_res, ni_sol
    ni_res = ml_norm_inf(res, mask)
    ni_sol = ml_norm_inf(sol, mask)
    r = ni_res <= eps*(Anorm*ni_sol + bnorm)
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

end subroutine ml_nd

end module ml_nd_module
