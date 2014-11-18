module boxlib_numpy
  use iso_c_binding
  use blobjects
  implicit none
contains

  subroutine multifab_as_numpy_f(cptr, nbox, dims, aptr) bind(c)
    type(c_ptr),    intent(in   ), value :: cptr
    integer(c_int), intent(in   ), value :: nbox
    type(c_ptr),    intent(  out)        :: aptr
    integer(c_int), intent(  out)        :: dims(4)

    real(8), pointer :: fptr(:,:,:,:)
    type(multifab), pointer :: mfab
    integer :: lo(4)

    call pybl_multifab_get(cptr,mfab)

    fptr => dataptr(mfab,nbox)
    lo = lbound(fptr)

    aptr = c_loc(fptr(lo(1),lo(2),lo(3),lo(4)))
    dims(1) = size(fptr,1)
    dims(2) = size(fptr,2)
    dims(3) = size(fptr,3)
    dims(4) = size(fptr,4)
  end subroutine multifab_as_numpy_f

  subroutine lmultifab_as_numpy_f(cptr, nbox, dims, aptr) bind(c)
    type(c_ptr),    intent(in   ), value :: cptr
    integer(c_int), intent(in   ), value :: nbox
    type(c_ptr),    intent(  out)        :: aptr
    integer(c_int), intent(  out)        :: dims(4)

    logical, pointer :: fptr(:,:,:,:)
    type(lmultifab), pointer :: mfab
    integer :: lo(4)

    call pybl_lmultifab_get(cptr,mfab)

    fptr => dataptr(mfab,nbox)
    lo = lbound(fptr)

    aptr = c_loc(fptr(lo(1),lo(2),lo(3),lo(4)))
    dims(1) = size(fptr,1)
    dims(2) = size(fptr,2)
    dims(3) = size(fptr,3)
    dims(4) = size(fptr,4)
  end subroutine lmultifab_as_numpy_f

  ! subroutine plotfile_as_numpy_f(cptr, level, nbox, aptr, n1, n2, n3, n4) &
  !      bind(c, name='plotfile_as_numpy_f')
  !   use iso_c_binding
  !   implicit none

  !   type(c_ptr),    intent(in)  :: cptr
  !   integer(c_int), intent(in)  :: level, nbox
  !   type(c_ptr),    intent(out) :: aptr
  !   integer(c_int), intent(out) :: n1, n2, n3, n4

  !   real(8), pointer :: fptr(:,:,:,:)
  !   type(plotfile), pointer :: pf
  !   integer :: lo(4)

  !   call pybl_plotfile_get(cptr,pf)

  !   fptr => dataptr(pf,level,nbox)
  !   lo = lbound(fptr)

  !   aptr = c_loc(fptr(lo(1),lo(2),lo(3),lo(4)))
  !   n1   = size(fptr,1)
  !   n2   = size(fptr,2)
  !   n3   = size(fptr,3)
  !   n4   = size(fptr,4)
  ! end subroutine plotfile_as_numpy_f

end module boxlib_numpy
