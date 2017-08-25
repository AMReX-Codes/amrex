module warpx_filter_module

  ! 1D, 3-point filter:
  use amrex_fort_module, only : amrex_real
  real(amrex_real), dimension(0:2), parameter :: ff1_1D = &
       (/ 0.25e0_rt, 0.5e0_rt 0.25e0_rt/)

  ! 2D, 3-point filter:
  real(amrex_real), dimension(0:2, 0:2), parameter :: ff1_2D = reshape( &
       [ 0.0625e0_rt, 0.125e0_rt, 0.0625e0_rt, &
         0.125e0_rt,  0.25e0_rt,  0.125e0_rt,  &
         0.0625e0_rt, 0.125e0_rt, 0.0625e0_rt], &
       shape(ff1_2D) )

  ! 3D, 3-point filter:
  real(amrex_real), dimension(0:2, 0:2, 0:2), parameter :: ff1_3D = reshape( &
       [ 0.015625e0_rt, 0.03125e0_rt,  0.015625e0_rt,  &
         0.03125e0_rt,  0.0625e0_rt,   0.03125e0_rt,   &
         0.015625e0_rt, 0.03125e0_rt,  0.015625e0_rt,  &
         0.03125e0_rt,  0.0625e0_rt,   0.03125e0_rt,   &
         0.0625e0_rt,   0.125e0_rt,    0.0625e0_rt,    &
         0.03125e0_rt,  0.0625e0_rt,   0.03125e0_rt,   &
         0.015625e0_rt, 0.03125e0_rt,  0.015625e0_rt,  &
         0.03125e0_rt,  0.0625e0_rt,   0.03125e0_rt,   &
         0.015625e0_rt, 0.03125e0_rt,  0.015625e0_rt], &
       shape(ff1_3D) )

end module warpx_filter_module
