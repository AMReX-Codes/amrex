module blobjects
  use multifab_module
  use layout_module
  use plotfile_module
  use ml_layout_module
contains
  
subroutine pybl_lmultifab_get(cptr,fptr)
  use iso_c_binding
  type(c_ptr), intent(in), value :: cptr
  type(lmultifab), pointer, intent(out) :: fptr

  call c_f_pointer(cptr, fptr)
end subroutine pybl_lmultifab_get

subroutine pybl_lmultifab_new(cptr,fptr)
  use iso_c_binding
  type(c_ptr), intent(out) :: cptr
  type(lmultifab), pointer, intent(out) :: fptr

  allocate(fptr)
  cptr = c_loc(fptr)
end subroutine pybl_lmultifab_new


subroutine pybl_multifab_get(cptr,fptr)
  use iso_c_binding
  type(c_ptr), intent(in), value :: cptr
  type(multifab), pointer, intent(out) :: fptr

  call c_f_pointer(cptr, fptr)
end subroutine pybl_multifab_get

subroutine pybl_multifab_new(cptr,fptr)
  use iso_c_binding
  type(c_ptr), intent(out) :: cptr
  type(multifab), pointer, intent(out) :: fptr

  allocate(fptr)
  cptr = c_loc(fptr)
end subroutine pybl_multifab_new


subroutine pybl_ml_layout_get(cptr,fptr)
  use iso_c_binding
  type(c_ptr), intent(in), value :: cptr
  type(ml_layout), pointer, intent(out) :: fptr

  call c_f_pointer(cptr, fptr)
end subroutine pybl_ml_layout_get

subroutine pybl_ml_layout_new(cptr,fptr)
  use iso_c_binding
  type(c_ptr), intent(out) :: cptr
  type(ml_layout), pointer, intent(out) :: fptr

  allocate(fptr)
  cptr = c_loc(fptr)
end subroutine pybl_ml_layout_new


subroutine pybl_layout_get(cptr,fptr)
  use iso_c_binding
  type(c_ptr), intent(in), value :: cptr
  type(layout), pointer, intent(out) :: fptr

  call c_f_pointer(cptr, fptr)
end subroutine pybl_layout_get

subroutine pybl_layout_new(cptr,fptr)
  use iso_c_binding
  type(c_ptr), intent(out) :: cptr
  type(layout), pointer, intent(out) :: fptr

  allocate(fptr)
  cptr = c_loc(fptr)
end subroutine pybl_layout_new


subroutine pybl_boxarray_get(cptr,fptr)
  use iso_c_binding
  type(c_ptr), intent(in), value :: cptr
  type(boxarray), pointer, intent(out) :: fptr

  call c_f_pointer(cptr, fptr)
end subroutine pybl_boxarray_get

subroutine pybl_boxarray_new(cptr,fptr)
  use iso_c_binding
  type(c_ptr), intent(out) :: cptr
  type(boxarray), pointer, intent(out) :: fptr

  allocate(fptr)
  cptr = c_loc(fptr)
end subroutine pybl_boxarray_new


subroutine pybl_plotfile_get(cptr,fptr)
  use iso_c_binding
  type(c_ptr), intent(in), value :: cptr
  type(plotfile), pointer, intent(out) :: fptr

  call c_f_pointer(cptr, fptr)
end subroutine pybl_plotfile_get

subroutine pybl_plotfile_new(cptr,fptr)
  use iso_c_binding
  type(c_ptr), intent(out) :: cptr
  type(plotfile), pointer, intent(out) :: fptr

  allocate(fptr)
  cptr = c_loc(fptr)
end subroutine pybl_plotfile_new

end module blobjects
