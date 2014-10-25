module ml_optimization_module

  implicit none

  ! 0: do nothing
  ! 1: sfc on each level;  ignore fine when distribute;  keep sfc order
  ! 2: do sfc on the finest level first and keep its order; work our way
  !    down; mark coarse grids with fine proc id; honor them if not over
  !    volpercpu; do rest with sfc
  ! 3: work our way up; mark fine grids with coarse proc id.  Do sfc and
  !    cut into chunks.  Let the ones that can benefit most pick first.  Then
  !    let the ones with most works pick.  Try to think how to minimize mpi
  !    gather.
  integer, save :: ml_optimization_strategy = 1

  private

  public :: set_ml_optimization_strategy

contains

  subroutine set_ml_optimization_strategy(i)
    integer, intent(in) :: i
    ml_optimization_strategy = i
  end subroutine set_ml_optimization_strategy

  

end module ml_optimization_module

