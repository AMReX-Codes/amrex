module turb_module
  use amrex_fort_module, only : amrex_real, amrex_spacedim
  implicit none

  public

contains


  subroutine fliprowsy(u, ulo, uhi) bind(C, name="fliprowsy")
    integer, intent(in) :: ulo(2), uhi(2)
    real (kind=amrex_real),intent(inout) :: u(ulo(1):uhi(1),ulo(2):uhi(2))
    real (kind=amrex_real) :: tmp
    integer :: i,j,jmid,k

    if (amrex_spacedim.eq.2) then
       k    = 0
       jmid = (ulo(2) + uhi(2)) / 2

       do j = ulo(2), jmid
          do i = ulo(1), uhi(1)
             tmp           = u(i,j)
             u(i,j)        = u(i,uhi(2)-k)
             u(i,uhi(2)-k) = tmp
          end do
          k = k + 1
       end do
    else
       call bl_pd_abort('Should not be here for 3D')
    endif
  end subroutine fliprowsy

end module turb_module
