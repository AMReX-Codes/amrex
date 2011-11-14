module boxlib_numpy

  use blobjects

contains

  subroutine multifab_as_numpy_f(oid, nbox, aptr, nx, ny, nz, nc) bind(c)
    use iso_c_binding
    implicit none

    integer(c_int), intent(in)  :: oid, nbox
    type(c_ptr), intent(out)    :: aptr
    integer(c_int), intent(out) :: nx, ny, nz, nc

    real(8), pointer :: fptr(:,:,:,:)
    type(multifab), pointer :: mfab
    integer :: lo(4)

    call pybl_multifab_get(oid,mfab)

    fptr => dataptr(mfab,nbox)
    lo = lbound(fptr)

    aptr = c_loc(fptr(lo(1),lo(2),lo(3),lo(4)))
    nx   = size(fptr,1)
    ny   = size(fptr,2)
    nz   = size(fptr,3)
    nc   = size(fptr,4)

  end subroutine multifab_as_numpy_f

  subroutine lmultifab_as_numpy_f(oid, nbox, aptr, nx, ny, nz, nc) bind(c)
    use iso_c_binding
    implicit none

    integer(c_int), intent(in)  :: oid, nbox
    type(c_ptr), intent(out)    :: aptr
    integer(c_int), intent(out) :: nx, ny, nz, nc

    logical, pointer :: fptr(:,:,:,:)
    type(lmultifab), pointer :: mfab
    integer :: lo(4)

    call pybl_lmultifab_get(oid,mfab)

    fptr => dataptr(mfab,nbox)
    lo = lbound(fptr)

    aptr = c_loc(fptr(lo(1),lo(2),lo(3),lo(4)))
    nx   = size(fptr,1)
    ny   = size(fptr,2)
    nz   = size(fptr,3)
    nc   = size(fptr,4)

  end subroutine lmultifab_as_numpy_f

end module boxlib_numpy

