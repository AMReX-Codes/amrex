module initialize_module

  use multifab_module
  use variables_module

  implicit none

  private
  public :: initialize_from_scratch

contains

  subroutine initialize_from_scratch(la,dt,courno,dx,U)

    use init_data_module, only : init_data
    use time_module, only : time

    use probin_module, only: n_cellx, n_celly, n_cellz, prob_lo, prob_hi,&
         max_grid_size, pmask
    use derivative_stencil_module, only : stencil_ng

    type(layout),intent(inout) :: la
    real(dp_t), intent(inout) :: dt, courno
    real(dp_t), pointer :: dx(:)
    type(multifab), intent(inout) :: U

    ! local
    integer :: lo(3), hi(3), ng
    type(box)          :: bx
    type(boxarray)     :: ba

    time = ZERO
    dt   = 1.d20
    courno = -1.d20

    lo = 0
    hi(1) = n_cellx-1
    hi(2) = n_celly-1
    hi(3) = n_cellz - 1

    bx = make_box(lo,hi)
    
    call boxarray_build_bx(ba,bx)
    call boxarray_maxsize(ba,max_grid_size)
    call layout_build_ba(la,ba,boxarray_bbox(ba),pmask=pmask)
    call destroy(ba)

    allocate(dx(3))
    dx(1) = (prob_hi(1)-prob_lo(1)) / n_cellx
    dx(2) = (prob_hi(2)-prob_lo(2)) / n_celly
    dx(3) = (prob_hi(3)-prob_lo(3)) / n_cellz

    ng = stencil_ng

    call multifab_build(U,la,ncons,ng)
  
    call init_data(U,dx,prob_lo,prob_hi)

  end subroutine initialize_from_scratch

end module initialize_module
