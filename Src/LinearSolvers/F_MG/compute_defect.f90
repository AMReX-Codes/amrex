module stencil_defect_module

  use bl_constants_module
  use bl_types
  use multifab_module
  use cc_stencil_module
  use stencil_types_module

  implicit none

  private
  public :: compute_defect, stencil_apply

contains

  ! Computes dd = ff - ss * uu
  subroutine compute_defect(ss, dd, ff, uu, mm, stencil_type, lcross, &
                            uniform_dh, bottom_solver, diagonalize, filled)

    use bl_prof_module

    type(multifab), intent(in)    :: ff, ss
    type(multifab), intent(inout) :: dd, uu
    type(imultifab), intent(in)   :: mm
    integer, intent(in)           :: stencil_type
    logical, intent(in)           :: lcross
    logical, intent(in), optional :: uniform_dh, bottom_solver, diagonalize, filled
    type(bl_prof_timer), save     :: bpt

    call bl_proffortfuncstart("compute_defect")
    call build(bpt, "compute_defect")
    call stencil_apply(ss, dd, uu, mm, stencil_type, lcross, &
                       uniform_dh, bottom_solver, diagonalize, filled)
    call sub_sub(dd, ff)
    call rescale(dd, -one)
    call destroy(bpt)
    call bl_proffortfuncstop("compute_defect")

  end subroutine compute_defect

  !
  ! Computes rr = aa * uu
  !
  subroutine stencil_apply(aa, rr, uu, mm, stencil_type, lcross, &
                           uniform_dh, bottom_solver, diagonalize, filled)

    use bl_prof_module

    use cc_stencil_apply_module, only : stencil_apply_1d, stencil_apply_2d, stencil_apply_3d, &
         stencil_apply_ibc_2d, stencil_apply_ibc_3d
    use nodal_stencil_apply_module, only: stencil_apply_1d_nodal, &
                                          stencil_apply_2d_nodal, &
                                          stencil_apply_3d_nodal
    use stencil_util_module, only : is_ibc_stencil

    type(multifab), intent(in)    :: aa
    type(multifab), intent(inout) :: rr
    type(multifab), intent(inout) :: uu
    type(imultifab), intent(in)   :: mm
    integer, intent(in)           :: stencil_type
    logical, intent(in)           :: lcross
    logical, intent(in),optional  :: uniform_dh, bottom_solver, diagonalize, filled

    real(kind=dp_t), pointer :: rp(:,:,:,:), up(:,:,:,:), ap(:,:,:,:)
    integer        , pointer :: mp(:,:,:,:)
    integer                  :: i, n, lo(get_dim(rr)), hi(get_dim(rr)), dm
    logical                  :: nodal_flag, luniform_dh, lbottom_solver, ldiagonalize, lfilled

    type(bl_prof_timer), save :: bpt

    call build(bpt, "its_stencil_apply")

    luniform_dh    = .false. ; if ( present(uniform_dh)    ) luniform_dh    = uniform_dh
    lbottom_solver = .false. ; if ( present(bottom_solver) ) lbottom_solver = bottom_solver
    ldiagonalize   = .false. ; if ( present(diagonalize)   ) ldiagonalize   = diagonalize
    lfilled        = .false. ; if ( present(filled)        ) lfilled        = filled

    if (ldiagonalize .and. .not. nodal_q(rr)) then
       call bl_error("Dont set diagonalize flag  = true for cell-centered in stencil_apply")
    end if

    call bl_assert(ncomp(uu).ge.ncomp(rr), 'uu must have at least as many components as rr')

    if (.not.lfilled) call fill_boundary(uu, 1, ncomp(rr), cross = lcross)

    dm         = get_dim(rr)
    nodal_flag = nodal_q(uu)

    do i = 1, nfabs(rr)
       rp => dataptr(rr, i)
       up => dataptr(uu, i)
       ap => dataptr(aa, i)
       lo = lwb(get_box(uu,i))
       hi = upb(get_box(uu,i))

       if (is_ibc_stencil(aa,i)) then
          do n = 1, ncomp(rr)
             select case(dm)
             case (2)
                call stencil_apply_ibc_2d(ap(:,1,1,1), rp(:,:,1,n), nghost(rr), &
                     up(:,:,1,n), nghost(uu), lo, hi)
             case (3)
                call stencil_apply_ibc_3d(ap(:,1,1,1), rp(:,:,:,n), nghost(rr), &
                     up(:,:,:,n), nghost(uu), lo, hi)
             end select
          end do
       else
          mp => dataptr(mm, i)
          do n = 1, ncomp(rr)
             select case(dm)
             case (1)
                if ( .not. nodal_flag) then
                   call stencil_apply_1d(ap(:,:,1,1), rp(:,1,1,n), nghost(rr), up(:,1,1,n), nghost(uu),  &
                        mp(:,1,1,1), lo, hi)
                else
                   call stencil_apply_1d_nodal(ap(1,:,1,1), rp(:,1,1,n), up(:,1,1,n),  &
                        mp(:,1,1,1), nghost(uu), nghost(rr), ldiagonalize)
                end if
             case (2)
                if ( .not. nodal_flag) then
                   call stencil_apply_2d(ap(:,:,:,1), rp(:,:,1,n), nghost(rr), up(:,:,1,n), nghost(uu),  &
                        mp(:,:,1,1), lo, hi)
                else
                   call stencil_apply_2d_nodal(ap(1,:,:,1), rp(:,:,1,n), up(:,:,1,n),  &
                        mp(:,:,1,1), nghost(uu), nghost(rr), stencil_type, ldiagonalize)
                end if
             case (3)
                if ( .not. nodal_flag) then
                   call stencil_apply_3d(ap(:,:,:,:), rp(:,:,:,n), nghost(rr), up(:,:,:,n), nghost(uu),  &
                        mp(:,:,:,1), bottom_solver=lbottom_solver)
                else
                   call stencil_apply_3d_nodal(ap(1,:,:,:), rp(:,:,:,n), up(:,:,:,n),  &
                        mp(:,:,:,1), nghost(uu), nghost(rr), stencil_type, luniform_dh, lbottom_solver, ldiagonalize)
                end if
             end select
          end do
       end if
    end do

    call destroy(bpt)

  end subroutine stencil_apply

end module stencil_defect_module
