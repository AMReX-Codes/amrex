subroutine t_bx
  use box_module
  implicit none
  type(box) :: bx
  bx = allbox(2)
  print *, bx
  bx = allbox_grc(bx, grow = 2)
  print *, bx
contains
  function allbox_grc(bx, grow, refine, coarsen) result(r)
    type(box) :: r
    type(box), intent(in) :: bx
    integer, intent(in), optional :: grow, refine, coarsen
    integer :: dm
    r = bx
    dm = r%dim
    if ( present(grow) ) then
    else if ( present(refine) ) then
    else if ( present(coarsen) ) then
    end if
  end function allbox_grc
end subroutine t_bx
