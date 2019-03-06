module amrex_eb_amr_util_nd_module

  implicit none
  private

  public :: amrex_tag_cutcells
  public :: amrex_tag_volfrac

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


  subroutine amrex_tag_volfrac(lo,    hi,                 &
       &                       tag,   tag_lo,   tag_hi,   &
       &                       vfrac, vfrac_lo, vfrac_hi, &
       &                       set,   clear,    tol      )&
       &     bind(C, name="amrex_tag_volfrac")

    use iso_c_binding, only: c_int, c_char
    use amrex_fort_module, only : amrex_real

    integer(c_int),  dimension(3), intent(in   ) :: lo, hi, vfrac_lo, vfrac_hi, tag_lo, tag_hi
    character(c_char),             value         :: set, clear
    real(amrex_real),              value         :: tol

    real(amrex_real),  intent(in   ) :: vfrac(vfrac_lo(1):vfrac_hi(1), &
         &                                    vfrac_lo(2):vfrac_hi(2), &
         &                                    vfrac_lo(3):vfrac_hi(3))
    character(c_char), intent(inout) :: tag(tag_lo(1):tag_hi(1), &
         &                                  tag_lo(2):tag_hi(2), &
         &                                  tag_lo(3):tag_hi(3))

    integer  :: i, j, k

    ! Tag on regions of high phi
    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if (( vfrac(i, j, k) .le. (1. - tol)).and.(vfrac(i, j, k) .ge. tol)) tag(i, j, k) = set

          enddo
       enddo
    enddo

  end subroutine amrex_tag_volfrac



end module amrex_eb_amr_util_nd_module
