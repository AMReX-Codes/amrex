module amrex_eb_amr_util_nd_module

  implicit none
  private

  public :: amrex_tag_cutcells

contains

  subroutine amrex_tag_cutcells (lo, hi, tag, tlo, thi, flag, flo, fhi, tagval, clearval) &
       bind(c,name='amrex_tag_cutcells')
    use iso_c_binding, only : c_char
    use amrex_ebcellflag_module, only : is_single_valued_cell
    integer, dimension(3), intent(in) :: lo, hi, tlo, thi, flo, fhi
    character(kind=c_char), intent(inout) :: tag(tlo(1):thi(1),tlo(2):thi(2),tlo(3):thi(3))
    integer,                intent(in)   :: flag(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))
    character(kind=c_char), value :: tagval, clearval

    integer :: i,j,k

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (is_single_valued_cell(flag(i,j,k))) then
                tag(i,j,k) = tagval
             end if
          end do
       end do
    end do
  end subroutine amrex_tag_cutcells

end module amrex_eb_amr_util_nd_module
