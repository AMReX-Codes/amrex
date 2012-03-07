module boxlib_numpy

  use blobjects

contains

  subroutine multifab_as_numpy_f(cptr, nbox, aptr, n1, n2, n3, n4) &
       bind(c, name='multifab_as_numpy_f')
    use iso_c_binding
    implicit none

    type(c_ptr),    intent(in)  :: cptr
    integer(c_int), intent(in)  :: nbox
    type(c_ptr),    intent(out) :: aptr
    integer(c_int), intent(out) :: n1, n2, n3, n4

    real(8), pointer :: fptr(:,:,:,:)
    type(multifab), pointer :: mfab
    integer :: lo(4)

    call pybl_multifab_get(cptr,mfab)

    fptr => dataptr(mfab,nbox)
    lo = lbound(fptr)

    aptr = c_loc(fptr(lo(1),lo(2),lo(3),lo(4)))
    n1   = size(fptr,1)
    n2   = size(fptr,2)
    n3   = size(fptr,3)
    n4   = size(fptr,4)
  end subroutine multifab_as_numpy_f

  ! subroutine lmultifab_as_numpy_f(cptr, nbox, aptr, n1, n2, n3, n4) bind(c)
  !   use iso_c_binding
  !   implicit none

  !   integer(c_int), intent(in)  :: cptr, nbox
  !   type(c_ptr), intent(out)    :: aptr
  !   integer(c_int), intent(out) :: n1, n2, n3, n4

  !   logical, pointer :: fptr(:,:,:,:)
  !   type(lmultifab), pointer :: mfab
  !   integer :: lo(4)

  !   call pybl_lmultifab_get(cptr,mfab)

  !   fptr => dataptr(mfab,nbox)
  !   lo = lbound(fptr)

  !   aptr = c_loc(fptr(lo(1),lo(2),lo(3),lo(4)))
  !   n1   = size(fptr,1)
  !   n2   = size(fptr,2)
  !   n3   = size(fptr,3)
  !   n4   = size(fptr,4)
  ! end subroutine lmultifab_as_numpy_f

end module boxlib_numpy

