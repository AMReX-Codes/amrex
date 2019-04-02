
module my_kernel_module
  use iso_c_binding
  use amrex_fort_module, only : amrex_real
  implicit none

contains

#ifdef AMREX_USE_CUDA_FORTRAN   
  AMREX_CUDA_FORT_DEVICE subroutine plusone_cudafort (lo, hi, dat, dlo, dhi) &
       bind(c,name='plusone_cudafort')
    integer(c_int), intent(in) :: lo(3), hi(3), dlo(3), dhi(3)
    real(amrex_real), intent(inout) :: dat(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3))

    integer(c_int) :: i,j,k

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             dat(i,j,k) = dat(i,j,k) + 1.0_amrex_real
          end do
       end do
    end do
  end subroutine plusone_cudafort
#endif
  
#ifdef AMREX_USE_ACC
  subroutine plusone_acc (lo, hi, dat, dlo, dhi) &
       bind(c,name='plusone_acc')
    integer(c_int), intent(in) :: lo(3), hi(3), dlo(3), dhi(3)
    real(amrex_real), intent(inout) :: dat(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3))

    integer(c_int) :: i,j,k

    !$acc kernels deviceptr(dat)
    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             dat(i,j,k) = dat(i,j,k) + 1.0_amrex_real
          end do
       end do
    end do
    !$acc end kernels
 end subroutine plusone_acc
#endif
 
#ifdef AMREX_OMP_OFFLOAD
  subroutine plusone_omp (lo, hi, dat, dlo, dhi) &
       bind(c,name='plusone_omp')
    integer(c_int), intent(in) :: lo(3), hi(3), dlo(3), dhi(3)
    real(amrex_real), intent(inout) :: dat(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3))

    integer(c_int) :: i,j,k

    !$omp target teams distribute parallel do collapse(3) schedule(static,1) is_device_ptr(dat)
    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             dat(i,j,k) = dat(i,j,k) + 1.0_amrex_real
          end do
       end do
    end do
 end subroutine plusone_omp
#endif

end module my_kernel_module
