subroutine set_simd (simd_width_in) bind(C, name='set_simd')

   use misc_params, only: simd_width
   implicit none

   integer, intent(in) :: simd_width_in

   simd_width = simd_width_in

end subroutine set_simd

subroutine fort_alloc_simd_vec() bind(C, name='fort_alloc_simd_vec')
  use misc_params, only: simd_width
  use vode_aux_module, only: T_vode_vec, ne_vode_vec, rho_vode_vec
  use amrex_error_module, only: amrex_abort
  implicit none

  !$omp parallel
  if (allocated(T_vode_vec) .or. allocated(ne_vode_vec) .or. allocated(rho_vode_vec)) then
    !$omp single
    call amrex_abort("Why are VODE SIMD vectors already allocated??")
    !$omp end single
  end if

  allocate(T_vode_vec(simd_width), ne_vode_vec(simd_width), rho_vode_vec(simd_width))
  !$omp end parallel
end subroutine fort_alloc_simd_vec


subroutine fort_dealloc_simd_vec() bind(C, name='fort_dealloc_simd_vec')
  use vode_aux_module, only: T_vode_vec, ne_vode_vec, rho_vode_vec
  use amrex_error_module, only: amrex_abort
  implicit none

  !$omp parallel
  if (.not. (allocated(T_vode_vec) .and. allocated(ne_vode_vec) .and. allocated(rho_vode_vec))) then
    !$omp single
    call amrex_abort("Why are VODE SIMD vectors already deallocated??")
    !$omp end single
  end if

  deallocate(T_vode_vec, ne_vode_vec, rho_vode_vec)
  !$omp end parallel
end subroutine fort_dealloc_simd_vec
