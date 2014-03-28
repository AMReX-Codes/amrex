module nodal_mg_tower_smoother_module
 
  use multifab_module
  use cc_stencil_module
  use mg_tower_module
  use bl_timer_module
 
  implicit none

contains

  subroutine nodal_mg_tower_smoother(mgt, lev, ss, uu, ff, mm)

    use omp_module
    use bl_prof_module
    use nodal_smoothers_module, only: nodal_smoother_1d, &
         nodal_smoother_2d, nodal_smoother_3d, nodal_line_solve_1d

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
    integer :: i, k, n, ng
    integer :: lo(mgt%dim)
    type(bl_prof_timer), save :: bpt
    logical :: pmask(mgt%dim)

    pmask = get_pmask(get_layout(uu))
    !
    ! Make sure to define this here so we don't assume a certain number of ghost cells for uu.
    !
    ng = nghost(uu)

    call build(bpt, "mgt_smoother")

    !
    ! Nodal stencils.
    !
    if (mgt%lcross) then
          !
          ! Cross stencils.
          !
          if ( get_dim(ff) == 1 ) then

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
                      call nodal_line_solve_1d(sp(:,:,1,1), up(:,1,1,n), fp(:,1,1,n), &
                           mp(:,1,1,1), lo, ng)
                   end select
                end do
             end do

          else

             do k = 0, 1
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
                         if ( k.eq.0 ) &
                              call nodal_line_solve_1d(sp(:,:,1,1), up(:,1,1,n), &
                              fp(:,1,1,n), mp(:,1,1,1), lo, ng)
                      case (2)
                         call nodal_smoother_2d(sp(:,:,:,1), up(:,:,1,n), &
                              fp(:,:,1,n), mp(:,:,1,1), lo, ng, &
                              pmask, mgt%stencil_type, k)
                      case (3)
                         call nodal_smoother_3d(sp(:,:,:,:), up(:,:,:,n), &
                              fp(:,:,:,n), mp(:,:,:,1), lo, ng, &
                              mgt%uniform_dh, pmask, mgt%stencil_type, k)
                      end select
                   end do
                end do
             end do
          end if
       else
          !
          ! Nodal dense stencils.
          !
          call fill_boundary(uu, cross = mgt%lcross)
          !
          ! Thread over FABS.
          !
          !$OMP PARALLEL DO PRIVATE(i,n,up,fp,sp,mp,lo)
          do i = 1, nfabs(ff)
             up => dataptr(uu, i)
             fp => dataptr(ff, i)
             sp => dataptr(ss, i)
             mp => dataptr(mm, i)
             lo =  lwb(get_box(ss, i))
             do n = 1, mgt%nc
                select case (mgt%dim)
                case (1)
                   call nodal_line_solve_1d(sp(:,:,1,1), up(:,1,1,n), &
                        fp(:,1,1,n), mp(:,1,1,1), lo, ng)
                case (2)
                   call nodal_smoother_2d(sp(:,:,:,1), up(:,:,1,n), &
                        fp(:,:,1,n), mp(:,:,1,1), lo, ng, &
                        pmask, mgt%stencil_type, 0)
                case (3)
                   call nodal_smoother_3d(sp(:,:,:,:), up(:,:,:,n), &
                        fp(:,:,:,n), mp(:,:,:,1), lo, ng, &
                        mgt%uniform_dh, pmask, mgt%stencil_type, 0)
                end select
             end do
          end do
          !$OMP END PARALLEL DO
       endif

       call multifab_internal_sync(uu)

    call destroy(bpt)

  end subroutine nodal_mg_tower_smoother

end module nodal_mg_tower_smoother_module
