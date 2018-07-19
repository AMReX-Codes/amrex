module amrex_habec_module

  ! habec is Hypre abec, where abec is the form of the linear equation
  ! we are solving:
  ! 
  ! alpha*phi - div(beta*grad phi) + div(\vec{c}*phi) 

  use iso_c_binding
  use amrex_fort_module, only : rt => amrex_real
  use amrex_lo_bctypes_module, only : amrex_lo_dirichlet, amrex_lo_neumann
  use amrex_error_module, only : amrex_error
  use amrex_constants_module, only : zero, half
  implicit none

contains

  subroutine amrex_hpacoef (lo, hi, mat, a, alo, ahi, sa) bind(c,name='amrex_hpacoef')
    integer, intent(in) :: lo(2), hi(2), alo(2), ahi(2)
    real(rt), intent(inout) :: mat(0:4, lo(1): hi(1), lo(2): hi(2))
    real(rt), intent(in   ) ::   a(    alo(1):ahi(1),alo(2):ahi(2))
    real(rt), intent(in) :: sa
    integer :: i, j
    if (sa .eq. zero) then
       !$omp parallel do private(i,j)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             mat(0,i,j) = zero
          enddo
       enddo
       !$omp end parallel do
    else
       !$omp parallel do private(i,j)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             mat(0,i,j) = sa * a(i,j)
          enddo
       enddo
       !$omp end parallel do
    end if
  end subroutine amrex_hpacoef


  subroutine amrex_hpbcoef (lo, hi, mat, b, blo, bhi, sb, dx, idim) &
       bind(c,name='amrex_hpbcoef')
    integer, intent(in) :: lo(2), hi(2), blo(2), bhi(2), idim
    real(rt), intent(inout) :: mat(0:4, lo(1): hi(1), lo(2): hi(2))
    real(rt), intent(in   ) ::   b(    blo(1):bhi(1),blo(2):bhi(2))
    real(rt), intent(in) :: sb, dx(2)

    integer :: i, j
    real(rt) :: fac

    fac = sb / dx(idim+1)**2 

    if (idim .eq. 0) then
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             mat(0,i,j) = mat(0,i,j) + fac * (b(i,j) + b(i+1,j))
             mat(1,i,j) = - fac * b(i,j)
             mat(2,i,j) = - fac * b(i+1,j)
          enddo
       enddo
    else
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             mat(0,i,j) = mat(0,i,j) + fac * (b(i,j) + b(i,j+1))
             mat(3,i,j) = - fac * b(i,j)
             mat(4,i,j) = - fac * b(i,j+1)
          enddo
       enddo
    endif
    
  end subroutine amrex_hpbcoef

end module amrex_habec_module
