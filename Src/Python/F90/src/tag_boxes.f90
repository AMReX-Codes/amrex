module tag_boxes_module
  use multifab_module
  use bl_error_module
  use iso_c_binding
  implicit none
  logical, save :: tagging_needs_ghost_cells = .true.

  interface
     subroutine tag_boxes_p(mf,tagboxes,dx,lev)
       use iso_c_binding
       type(c_ptr), intent(in), value :: mf, tagboxes
       real(c_double), intent(in), value :: dx
       integer(c_int), intent(in), value :: lev
     end subroutine tag_boxes_p
  end interface
contains

  subroutine tag_boxes(mf,tagboxes,dx,lev,tag_boxes_cb)
    type( multifab)         , intent(in   ), target :: mf
    type(lmultifab)         , intent(inout), target :: tagboxes
    real(dp_t)              , intent(in   ) :: dx
    integer                 , intent(in   ) :: lev
    type(c_funptr)          , intent(in   ) :: tag_boxes_cb

    procedure(tag_boxes_p), pointer :: tag
    call c_f_procpointer(tag_boxes_cb, tag)
    call tag(c_loc(mf), c_loc(tagboxes), dx, lev)
  end subroutine tag_boxes

end module tag_boxes_module
