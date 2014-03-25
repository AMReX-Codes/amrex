module nodal_applyop_module

  use ml_layout_module
  use mg_module
  use define_bc_module
  use nodal_stencil_fill_module

  implicit none

  private

  public :: nodal_applyop

contains

  subroutine nodal_applyop(mla,res,phi,beta,dx,the_bc_tower,bc_comp)

    ! compute res = del dot beta grad phi

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: res(:)  ! nodal
    type(multifab) , intent(inout) :: phi(:)  ! nodal
    type(multifab) , intent(in   ) :: beta(:) ! cell-centered
    real(dp_t)     , intent(in   ) :: dx(:,:)
    type(bc_tower) , intent(in   ) :: the_bc_tower
    integer        , intent(in   ) :: bc_comp

    ! local
    integer :: n,nlevs,dm
    integer :: stencil_order

    type(box) :: pd
    type(mg_tower) :: mgt(mla%nlevel)

    type(multifab), allocatable :: coeffs(:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "nodal_applyop")

    dm = mla%dim
    nlevs = mla%nlevel

    stencil_order = 2

    do n = nlevs, 1, -1

       pd = layout_get_pd(mla%la(n))

       call mg_tower_build(mgt(n), mla%la(n), pd, &
                           the_bc_tower%bc_tower_array(n)%ell_bc_level_array(0,:,:,bc_comp),&
                           stencil_type = ND_DENSE_STENCIL, &
                           dh = dx(n,:), &
                           nodal = nodal_flags(res(nlevs)))

    end do

    do n = nlevs,1,-1

       allocate(coeffs(mgt(n)%nlevels))

       call multifab_build(coeffs(mgt(n)%nlevels), mla%la(n), 1, 1)
       call multifab_copy_c(coeffs(mgt(n)%nlevels),1,beta(n),1,1,beta(n)%ng)

       call stencil_fill_nodal_all_mglevels(mgt(n), coeffs)

       call destroy(coeffs(mgt(n)%nlevels))
       deallocate(coeffs)

    end do

    call ml_nodal_applyop(mla,mgt,res,phi)

    do n=1,nlevs
       call mg_tower_destroy(mgt(n))
    end do

    call destroy(bpt)

  end subroutine nodal_applyop

  subroutine ml_nodal_applyop(mla,mgt,res,full_soln)

    type(ml_layout), intent(in)    :: mla
    type(mg_tower) , intent(inout) :: mgt(:)
    type( multifab), intent(inout) :: res(:)
    type( multifab), intent(inout) :: full_soln(:)

    ! local
    integer    :: n,nlevs,mglev

    type(multifab), allocatable  ::        rh(:) ! this will be set to zero

    nlevs = mla%nlevel

    do n=1,nlevs
       call multifab_build_nodal(rh(n), mla%la(n), 1, 0)
       call setval(rh(n), 0.d0, all=.true.)
    end do

    do n = nlevs,1,-1
       mglev = mgt(n)%nlevels
       call compute_defect(mgt(n)%ss(mglev),res(n),rh(n),full_soln(n),mgt(n)%mm(mglev), &
                           mgt(n)%stencil_type, mgt(n)%lcross, mgt(n)%uniform_dh)
    end do

    do n=1,nlevs
       call multifab_destroy(rh(n))
    end do

  end subroutine ml_nodal_applyop

end module nodal_applyop_module
