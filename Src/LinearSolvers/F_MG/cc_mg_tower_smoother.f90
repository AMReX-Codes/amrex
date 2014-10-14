module cc_mg_tower_smoother_module
 
  use multifab_module
  use cc_stencil_module
  use mg_tower_module
  use bl_timer_module
 
  implicit none

contains

  subroutine cc_mg_tower_smoother(mgt, lev, ss, uu, ff, mm)

    use omp_module
    use bl_prof_module
    use cc_smoothers_module, only: gs_line_solve_1d, gs_rb_smoother_1d, &
         jac_smoother_2d, jac_smoother_3d, &
         gs_rb_smoother_2d,  gs_rb_smoother_3d, &
         fourth_order_smoother_2d, fourth_order_smoother_3d
    use itsol_module, only: itsol_bicgstab_solve

    integer        , intent(in   ) :: lev
    type( mg_tower), intent(inout) :: mgt
    type( multifab), intent(inout) :: uu
    type( multifab), intent(in   ) :: ff
    type( multifab), intent(in   ) :: ss
    type(imultifab), intent(in   ) :: mm

    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: up(:,:,:,:)
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    integer        , pointer :: mp(:,:,:,:)
    integer :: i, n, ng, nn, stat, npts
    integer :: lo(mgt%dim)
    type(bl_prof_timer), save :: bpt
    logical :: pmask(mgt%dim), singular_test
    real(kind=dp_t) :: local_eps

    singular_test = ( mgt%bottom_singular .and. mgt%coeffs_sum_to_zero )

    pmask = get_pmask(get_layout(uu))
    !
    ! Make sure to define this here so we don't assume a certain number of ghost cells for uu.
    !
    ng = nghost(uu)

    call build(bpt, "cc_mgt_smoother")

    if ( mgt%skewed_not_set(lev) ) then
       !$OMP PARALLEL DO PRIVATE(i,mp)
       do i = 1, nfabs(mm)
          mp => dataptr(mm, i)
          mgt%skewed(lev,i) = skewed_q(mp)
       end do
       !$OMP END PARALLEL DO
       mgt%skewed_not_set(lev) = .false.
    end if

    if ( (mgt%dim .eq. 1) .and. (nboxes(uu%la) .eq. 1) ) then

       call fill_boundary(uu, cross = mgt%lcross)
       !
       ! We do these line solves as a preconditioner.
       !
       do i = 1, nfabs(ff)
          up => dataptr(uu, i)
          fp => dataptr(ff, i)
          sp => dataptr(ss, i)
          mp => dataptr(mm, i)
          lo =  lwb(get_box(ss, i))
          do n = 1, mgt%nc
             call gs_line_solve_1d(sp(:,:,1,1), up(:,1,1,n), fp(:,1,1,n), &
                  mp(:,1,1,1), lo, ng)
          end do
       end do
       
       call fill_boundary(uu, cross = mgt%lcross)
       
       local_eps = 0.001 
       npts = multifab_volume(ff)
       call itsol_bicgstab_solve(ss, uu, ff, mm, &
            local_eps, npts, mgt%cg_verbose, &
            mgt%stencil_type, mgt%lcross, &
            stat = stat, &
            singular_in = singular_test, &
            uniform_dh = mgt%uniform_dh)
       
       call fill_boundary(uu, cross = mgt%lcross)
       
       if ( stat /= 0 ) then
          if ( parallel_IOProcessor() ) call bl_error("BREAKDOWN in 1d CG solve")
       end if
       
    else
       !
       ! Cell-centered stencils.
       !
       select case ( mgt%smoother )
          
       case ( MG_SMOOTHER_GS_RB )
          
          do nn = 0, 1
             call fill_boundary(uu, cross = mgt%lcross)
             do i = 1, nfabs(ff)
                up => dataptr(uu, i)
                fp => dataptr(ff, i)
                sp => dataptr(ss, i)
                mp => dataptr(mm, i)
                lo =  lwb(get_box(ss, i))
                do n = 1, mgt%nc
                   select case ( mgt%dim)
                   case (1)
                      call gs_rb_smoother_1d(sp(:,:,1,1), up(:,1,1,n), &
                           fp(:,1,1,n), mp(:,1,1,1), lo, ng, nn, &
                           mgt%skewed(lev,i))
                   case (2)
                      call gs_rb_smoother_2d(sp(:,:,:,1), up(:,:,1,n), &
                           fp(:,:,1,n), mp(:,:,1,1), lo, ng, nn, &
                           mgt%skewed(lev,i))
                   case (3)
                      call gs_rb_smoother_3d(sp(:,:,:,:), up(:,:,:,n), &
                           fp(:,:,:,n), mp(:,:,:,1), lo, ng, nn, &
                           mgt%skewed(lev,i))
                   end select
                end do
             end do
          end do
          
       case ( MG_SMOOTHER_EFF_RB )
          
          call fill_boundary(uu, cross = mgt%lcross)
          do i = 1, nfabs(ff)
             up => dataptr(uu, i)
             fp => dataptr(ff, i)
             sp => dataptr(ss, i)
             mp => dataptr(mm, i)
             lo =  lwb(get_box(ss, i))
             do n = 1, mgt%nc
                select case ( mgt%dim)
                case (1)
                   call gs_rb_smoother_1d(sp(:,:,1,1), up(:,1,1,n), &
                        fp(:,1,1,n), mp(:,1,1,1), lo, ng, 0, &
                        mgt%skewed(lev,i))
                   call gs_rb_smoother_1d(sp(:,:,1,1), up(:,1,1,n), &
                        fp(:,1,1,n), mp(:,1,1,1), lo, ng, 1, &
                        mgt%skewed(lev,i))
                case (2)
                   call gs_rb_smoother_2d(sp(:,:,:,1), up(:,:,1,n), &
                        fp(:,:,1,n), mp(:,:,1,1), lo, ng, 0, &
                        mgt%skewed(lev,i))
                   call gs_rb_smoother_2d(sp(:,:,:,1), up(:,:,1,n), &
                        fp(:,:,1,n), mp(:,:,1,1), lo, ng, 1, &
                        mgt%skewed(lev,i))
                case (3)
                   call gs_rb_smoother_3d(sp(:,:,:,:), up(:,:,:,n), &
                        fp(:,:,:,n), mp(:,:,:,1), lo, ng, 0, &
                        mgt%skewed(lev,i))
                   call gs_rb_smoother_3d(sp(:,:,:,:), up(:,:,:,n), &
                        fp(:,:,:,n), mp(:,:,:,1), lo, ng, 1, &
                        mgt%skewed(lev,i))
                end select
             end do
          end do
          
       case ( MG_SMOOTHER_MINION_CROSS )
          
          do nn = 0, 1
             if ( (nn == 1) .and. (mgt%dim == 3) ) &
                  !
                  ! Only the 2D fourth order stencil does red-black GS sweeps.
                  !
                  exit
             call fill_boundary(uu, cross = mgt%lcross)
             do i = 1, nfabs(ff)
                up => dataptr(uu, i)
                fp => dataptr(ff, i)
                sp => dataptr(ss, i)
                mp => dataptr(mm, i)
                lo =  lwb(get_box(ss, i))
                do n = 1, mgt%nc
                   select case (mgt%dim)
                   case (2)
                      call fourth_order_smoother_2d(sp(:,:,:,1), up(:,:,1,1), &
                           fp(:,:,1,1), lo, ng, mgt%stencil_type, nn)
                   case (3)
                      call fourth_order_smoother_3d(sp(:,:,:,:), up(:,:,:,1), &
                           fp(:,:,:,1), lo, ng, mgt%stencil_type)
                   end select
                end do
             end do
          end do
          
       case ( MG_SMOOTHER_MINION_FULL )
          
          call fill_boundary(uu, cross = mgt%lcross)
          do i = 1, nfabs(ff)
             up => dataptr(uu, i)
             fp => dataptr(ff, i)
             sp => dataptr(ss, i)
             mp => dataptr(mm, i)
             lo =  lwb(get_box(ss, i))
             do n = 1, mgt%nc
                select case (mgt%dim)
                case (2)
                   call fourth_order_smoother_2d(sp(:,:,:,1), up(:,:,1,1), &
                        fp(:,:,1,1), lo, ng, mgt%stencil_type, 0)
                case (3)
                   call fourth_order_smoother_3d(sp(:,:,:,:), up(:,:,:,1), &
                        fp(:,:,:,1), lo, ng, mgt%stencil_type)
                end select
             end do
          end do
          
       case ( MG_SMOOTHER_JACOBI )
          
          call fill_boundary(uu, cross = mgt%lcross)
          do i = 1, nfabs(ff)
             up => dataptr(uu, i)
             fp => dataptr(ff, i)
             sp => dataptr(ss, i)
             mp => dataptr(mm, i)
             lo =  lwb(get_box(ss, i))
             do n = 1, mgt%nc
                select case ( mgt%dim)
                case (2)
                   call jac_smoother_2d(sp(:,:,:,1), up(:,:,1,n), fp(:,:,1,n), &
                        mp(:,:,1,1), ng)
                case (3)
                   call jac_smoother_3d(sp(:,:,:,:), up(:,:,:,n), fp(:,:,:,n), &
                        mp(:,:,:,1), ng)
                end select
             end do
          end do
          
       case default
          call bl_error("MG_TOWER_SMOOTHER: no such smoother")
       end select
       
    end if ! if (mgt%dim > 1)
    
    call destroy(bpt)

  end subroutine cc_mg_tower_smoother

  subroutine mg_jacobi_smoother(mgt, lev, ss, uu, ff, mm)

    use bl_prof_module
    use cc_smoothers_module, only: jac_smoother_1d, jac_smoother_2d, jac_smoother_3d

    type(mg_tower), intent(inout) :: mgt
    type(multifab), intent(inout) :: uu
    type(multifab), intent(in) :: ff
    type(multifab), intent(in) :: ss
    type(imultifab), intent(in) :: mm
    integer, intent(in) :: lev
    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: up(:,:,:,:)
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    integer        , pointer :: mp(:,:,:,:)
    integer :: i, n, ng, iter
    type(bl_prof_timer), save :: bpt

    call build(bpt, "mgt_jacobi_smoother")

    ng = nghost(uu)

    if ( mgt%skewed_not_set(lev) ) then 
       do i = 1, nfabs(mm)
          mp => dataptr(mm, i)
          mgt%skewed(lev,i) = skewed_q(mp)
       end do
       mgt%skewed_not_set(lev) = .false.
    end if

    if ( .not. cell_centered_q(uu) ) then
       call bl_error('MG_JACOBI_SMOOTHER ONLY DESIGNED FOR CELL-CENTERED RIGHT NOW')
    end if

    do iter = 1, mgt%nu1

       call fill_boundary(uu, cross = mgt%lcross)

       do i = 1, nfabs(ff)
          up => dataptr(uu, i)
          fp => dataptr(ff, i)
          sp => dataptr(ss, i)
          mp => dataptr(mm, i)
          do n = 1, mgt%nc
             select case ( mgt%dim)
             case (1)
                call jac_smoother_1d(sp(:,:,1,1), up(:,1,1,n), fp(:,1,1,n), &
                                     mp(:,1,1,1), ng)
             case (2)
                call jac_smoother_2d(sp(:,:,:,1), up(:,:,1,n), fp(:,:,1,n), &
                                     mp(:,:,1,1), ng)
             case (3)
                call jac_smoother_3d(sp(:,:,:,:), up(:,:,:,n), fp(:,:,:,n), &
                                     mp(:,:,:,1), ng)
             end select
          end do
       end do
    end do

    call destroy(bpt)

  end subroutine mg_jacobi_smoother

end module cc_mg_tower_smoother_module
