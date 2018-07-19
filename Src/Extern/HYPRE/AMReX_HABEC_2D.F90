module amrex_habec_module

  ! habec is Hypre abec, where abec is the form of the linear equation
  ! we are solving:
  ! 
  ! alpha*phi - div(beta*grad phi) + div(\vec{c}*phi) 

  use iso_c_binding
  use amrex_fort_module, only : rt => amrex_real
  use amrex_lo_bctypes_module, only : amrex_lo_dirichlet, amrex_lo_neumann
  use amrex_error_module, only : amrex_error
  use amrex_constants_module, only : zero, one, half, three
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
       !$omp parallel do private(i,j)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             mat(0,i,j) = mat(0,i,j) + fac * (b(i,j) + b(i+1,j))
             mat(1,i,j) = - fac * b(i,j)
             mat(2,i,j) = - fac * b(i+1,j)
          enddo
       enddo
       !$omp end parallel do
    else
       !$omp parallel do private(i,j)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             mat(0,i,j) = mat(0,i,j) + fac * (b(i,j) + b(i,j+1))
             mat(3,i,j) = - fac * b(i,j)
             mat(4,i,j) = - fac * b(i,j+1)
          enddo
       enddo
       !$omp end parallel do
    endif
    
  end subroutine amrex_hpbcoef


  subroutine amrex_hpmat (lo, hi, mat, b, blo, bhi, mask, mlo, mhi, &
       sb, dx, cdir, bct, bcl, bho) bind(c,name='amrex_hpmat')
    integer, intent(in) :: lo(2), hi(2), blo(2), bhi(2), mlo(2), mhi(2), cdir, bct, bho
    real(rt), intent(inout) ::  mat(0:4, lo(1): hi(1), lo(2): hi(2))
    real(rt), intent(in   ) ::    b(    blo(1):bhi(1),blo(2):bhi(2))
    integer , intent(in   ) :: mask(    mlo(1):mhi(1),mlo(2):mhi(2))
    real(rt), intent(in) :: sb, dx(2), bcl

    integer :: i, j
    real(rt) :: fac, h, h2, h3, bf1, bf2

    if (cdir .eq. 0 .or. cdir .eq. 2) then
       h = dx(1)
    else
       h = dx(2)
    endif
    fac = sb / (h**2)

    if (bct .eq. amrex_lo_dirichlet) then
       h2 = half * h
       if (bho.ge.1) then
          h3 = three * h2
          bf1 = fac * ((h3 - bcl) / (bcl + h2) - one)
          bf2 = fac * (bcl - h2) / (bcl + h3)
       else
          bf1 = fac * ( h / (bcl + h2) - one)          
          bf2 = zero
       end if
    else if (bct .eq. amrex_lo_neumann) then
       bf1 = -fac
       bf2 = zero
    else
       call amrex_error("hpmat: unsupported boundary type")       
    end if
    
    if (cdir .eq. 0) then
       i = lo(1)
       !$omp parallel do private(j)
       do j = lo(2), hi(2)
          if (mask(i-1,j) .gt. 0) then
             mat(0,i,j) = mat(0,i,j) + bf1 * b(i,j)
             mat(1,i,j) = zero
             mat(2,i,j) = mat(2,i,j) + bf2 * b(i,j)
          endif
       enddo
    else if (cdir .eq. 2) then
       i = hi(1)
       !$omp parallel do private(j)
       do j = lo(2), hi(2)
          if (mask(i+1,j) .gt. 0) then
             mat(0,i,j) = mat(0,i,j) + bf1 * b(i+1,j)
             mat(2,i,j) = zero
             mat(1,i,j) = mat(1,i,j) + bf2 * b(i+1,j)
          endif
       enddo
    else if (cdir .eq. 1) then
       j = lo(2)
       do i = lo(1), hi(1)
          if (mask(i,j-1) .gt. 0) then
             mat(0,i,j) = mat(0,i,j) + bf1 * b(i,j)
             mat(3,i,j) = zero
             mat(4,i,j) = mat(4,i,j) + bf2 * b(i,j)
          endif
       enddo
    else if (cdir .eq. 3) then
       j = hi(2)
       do i = lo(1), hi(1)
          if (mask(i,j+1) .gt. 0) then
             mat(0,i,j) = mat(0,i,j) + bf1 * b(i,j+1)
             mat(4,i,j) = zero
             mat(3,i,j) = mat(3,i,j) + bf2 * b(i,j+1)
          endif
       enddo
    else
       call amrex_error("hpmat: impossible face orientation")
    endif
    
  end subroutine amrex_hpmat

end module amrex_habec_module
