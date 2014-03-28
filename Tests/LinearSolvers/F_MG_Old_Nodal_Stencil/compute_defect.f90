module stencil_defect_module

  use bl_types
  use multifab_module
  use cc_stencil_module
  use stencil_types_module

  implicit none

contains

  ! Computes dd = ff - ss * uu
  subroutine compute_defect(ss, dd, ff, uu, mm, stencil_type, lcross, &
                            uniform_dh, bottom_solver)

    use bl_prof_module

    type(multifab), intent(in)    :: ff, ss
    type(multifab), intent(inout) :: dd, uu
    type(imultifab), intent(in)   :: mm
    integer, intent(in)           :: stencil_type
    logical, intent(in)           :: lcross
    logical, intent(in), optional :: uniform_dh, bottom_solver
    type(bl_prof_timer), save     :: bpt

    call build(bpt, "compute_defect")
    call stencil_apply(ss, dd, uu, mm, stencil_type, lcross, &
                       uniform_dh, bottom_solver)
    call saxpy(dd, ff, -1.0_dp_t, dd)
    call destroy(bpt)

  end subroutine compute_defect

  !
  ! Computes rr = aa * uu
  !
  subroutine stencil_apply(aa, rr, uu, mm, stencil_type, lcross, &
                           uniform_dh, bottom_solver)

    use bl_prof_module

    use cc_stencil_apply_module
    use nodal_stencil_apply_module, only: stencil_apply_1d_nodal, &
                                          stencil_apply_2d_nodal, &
                                          stencil_apply_3d_nodal

    type(multifab), intent(in)    :: aa
    type(multifab), intent(inout) :: rr
    type(multifab), intent(inout) :: uu
    type(imultifab), intent(in)   :: mm
    integer, intent(in)           :: stencil_type
    logical, intent(in)           :: lcross
    logical, intent(in),optional  :: uniform_dh, bottom_solver

    real(kind=dp_t), pointer :: rp(:,:,:,:), up(:,:,:,:), ap(:,:,:,:)
    integer        , pointer :: mp(:,:,:,:)
    integer                  :: i, n, lo(get_dim(rr)), hi(get_dim(rr)), dm
    logical                  :: nodal_flag, luniform_dh, lbottom_solver

    type(bl_prof_timer), save :: bpt

    call build(bpt, "its_stencil_apply")

    luniform_dh    = .false. ; if ( present(uniform_dh)    ) luniform_dh    = uniform_dh
    lbottom_solver = .false. ; if ( present(bottom_solver) ) lbottom_solver = bottom_solver

    call bl_assert(ncomp(uu).ge.ncomp(rr), 'uu must have at least as many components as rr')

    call fill_boundary(uu, 1, ncomp(rr), cross = lcross)

    dm         = get_dim(rr)
    nodal_flag = nodal_q(uu)

    do i = 1, nfabs(rr)
       rp => dataptr(rr, i)
       up => dataptr(uu, i)
       ap => dataptr(aa, i)
       mp => dataptr(mm, i)
       lo = lwb(get_box(uu,i))
       hi = upb(get_box(uu,i))
       do n = 1, ncomp(rr)
          select case(dm)
          case (1)
             if ( .not. nodal_flag) then
                call stencil_apply_1d(ap(:,:,1,1), rp(:,1,1,n), nghost(rr), up(:,1,1,n), nghost(uu),  &
                                      mp(:,1,1,1), lo, hi)
             else
                call stencil_apply_1d_nodal(ap(:,:,1,1), rp(:,1,1,n), up(:,1,1,n),  &
                     mp(:,1,1,1), nghost(uu))
             end if
          case (2)
             if ( .not. nodal_flag) then
                call stencil_apply_2d(ap(:,:,:,1), rp(:,:,1,n), nghost(rr), up(:,:,1,n), nghost(uu),  &
                     mp(:,:,1,1), lo, hi)
             else
                call stencil_apply_2d_nodal(ap(:,:,:,1), rp(:,:,1,n), up(:,:,1,n),  &
                     mp(:,:,1,1), nghost(uu), stencil_type)
             end if
          case (3)
             if ( .not. nodal_flag) then
                call stencil_apply_3d(ap(:,:,:,:), rp(:,:,:,n), nghost(rr), up(:,:,:,n), nghost(uu),  &
                                      mp(:,:,:,1), bottom_solver=lbottom_solver)
             else
                call stencil_apply_3d_nodal(ap(:,:,:,:), rp(:,:,:,n), up(:,:,:,n),  &
                     mp(:,:,:,1), nghost(uu), stencil_type, luniform_dh, lbottom_solver)
             end if
          end select
       end do
    end do

    call destroy(bpt)

  end subroutine stencil_apply

end module stencil_defect_module
