module transform_module

  use amrex_fort_module, only : amrex_real
  implicit none

  public

contains

  subroutine transform(lo, hi, sIn, sInlo, sInhi, ncIn, sOut, sOutlo, sOuthi, ncOut) bind(C, name="transform")
    integer, intent(in) :: lo(3),hi(3)
    integer, intent(in) :: sInlo(3),sInhi(3)
    integer, intent(in) :: sOutlo(3),sOuthi(3)
    integer, intent(in) :: ncIn, ncOut

    real (kind=amrex_real),intent(in   ) :: sIn(sInlo(1):sInhi(1),sInlo(2):sInhi(2),sInlo(3):sInhi(3),ncIn)
    real (kind=amrex_real),intent(inout) :: sOut(sOutlo(1):sOuthi(1),sOutlo(2):sOuthi(2),sOutlo(3):sOuthi(3),ncOut)

    integer :: i, j, k

    ! This is an example pointwise transformation
    ! Here, sIn(i,j,k,1LncIn) contains data from the plotfile, stacked up in the order that the 
    ! user called the function with.  The output data has (hardwired) ncOut components
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             sOut(i,j,k,1:ncOut-1) = 0.d0
             sOut(i,j,k,ncOut) = sIn(i,j,k,ncIn)
          enddo
       enddo
    enddo

  end subroutine transform
end module transform_module
