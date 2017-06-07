
attributes(global) subroutine double_mass(particles)
  
  use amrex_fort_module, only : amrex_real  
  implicit none
  
  integer              :: ns, np
  real(amrex_real)     :: particles(5, 20971520)
  integer i

  i = blockDim%x * (blockIdx%x - 1) + threadIdx%x
  if (i <= np) particles(4, i) = particles(4, i)*2.d0
 
end subroutine double_mass

subroutine process_particles(particles, ns, np) &
     bind(c,name='process_particles')
  
  use iso_c_binding
  use amrex_fort_module, only : amrex_real
  use cudafor
  
  implicit none
  
  integer, value   :: ns, np
  real(amrex_real) :: particles(5, 20971520)

  attributes(device) :: particles

  integer :: cuda_result
  integer, device    :: ns_d, np_d
  type(dim3) :: numThreads, numBlocks

  cuda_result = cudaMemcpy(ns_d, ns, 1, cudaMemcpyHostToDevice)
  cuda_result = cudaMemcpy(np_d, np,  1, cudaMemcpyHostToDevice)

  numBlocks  = dim3(256,1,1)
  numThreads = dim3(ceiling(real(np)/numBlocks%x),1,1)

  call double_mass<<<numBlocks, numThreads>>>(particles, ns_d, np_d)

end subroutine process_particles

attributes(global) subroutine deposit_kernel(particles, ns, np, rho, lo, hi, plo, dx)
  
  use amrex_fort_module, only : amrex_real  
  implicit none

  integer, value       :: ns, np  
  real(amrex_real)     :: particles(ns, np)
  integer, intent(in)  :: lo(3), hi(3)
  real(amrex_real)     :: rho(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))
  real(amrex_real)     :: plo(3)
  real(amrex_real)     :: dx(3)
  
  integer i, j, k, n
  real(amrex_real) wx_lo, wy_lo, wz_lo, wx_hi, wy_hi, wz_hi
  real(amrex_real) lx, ly, lz
  real(amrex_real) inv_dx(3)
  inv_dx = 1.0d0/dx

  n = blockDim%x * (blockIdx%x - 1) + threadIdx%x
  if (n <= np) then
     lx = (particles(1, n) - plo(1))*inv_dx(1) + 0.5d0
     ly = (particles(2, n) - plo(2))*inv_dx(2) + 0.5d0
     lz = (particles(3, n) - plo(3))*inv_dx(3) + 0.5d0
     
     i = floor(lx)
     j = floor(ly)
     k = floor(lz)
     
     wx_hi = lx - i
     wy_hi = ly - j
     wz_hi = lz - k
     
     wx_lo = 1.0d0 - wx_hi
     wy_lo = 1.0d0 - wy_hi
     wz_lo = 1.0d0 - wz_hi
     
     rho(i-1, j-1, k-1) = rho(i-1, j-1, k-1) + wx_lo*wy_lo*wz_lo*particles(4, n)
     rho(i-1, j-1, k)   = rho(i-1, j-1, k)   + wx_lo*wy_lo*wz_hi*particles(4, n)
     rho(i-1, j,   k-1) = rho(i-1, j,   k-1) + wx_lo*wy_hi*wz_lo*particles(4, n)
     rho(i-1, j,   k)   = rho(i-1, j,   k)   + wx_lo*wy_hi*wz_hi*particles(4, n)
     rho(i,   j-1, k-1) = rho(i,   j-1, k-1) + wx_hi*wy_lo*wz_lo*particles(4, n)
     rho(i,   j-1, k)   = rho(i,   j-1, k)   + wx_hi*wy_lo*wz_hi*particles(4, n)
     rho(i,   j,   k-1) = rho(i,   j,   k-1) + wx_hi*wy_hi*wz_lo*particles(4, n)
     rho(i,   j,   k)   = rho(i,   j,   k)   + wx_hi*wy_hi*wz_hi*particles(4, n)
  end if
 
end subroutine deposit_kernel

subroutine deposit_cic(particles, ns, np, rho, lo, hi, plo, dx) &
     bind(c,name='deposit_cic')
  
  use iso_c_binding
  use amrex_fort_module, only : amrex_real
  use cudafor
  
  implicit none
  
  integer, value       :: ns, np
  real(amrex_real)     :: particles(ns,np)
  integer              :: lo(3)
  integer              :: hi(3)
  real(amrex_real)     :: rho(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))
  real(amrex_real)     :: plo(3)
  real(amrex_real)     :: dx(3)
  
  attributes(device) :: particles
  attributes(device) :: rho
  integer,  device           :: lo_d(3), hi_d(3)
  real(amrex_real), device   :: plo_d(3), dx_d(3)
  integer, device            :: ns_d, np_d

  integer :: cuda_result
  integer(kind=cuda_stream_kind) :: stream
  type(dim3) :: numThreads, numBlocks

  cuda_result = cudaMemcpy(lo_d, lo, 3, cudaMemcpyHostToDevice)
  cuda_result = cudaMemcpy(hi_d, hi, 3, cudaMemcpyHostToDevice)

  cuda_result = cudaMemcpy(plo_d, plo, 3, cudaMemcpyHostToDevice)
  cuda_result = cudaMemcpy(dx_d,  dx,  3, cudaMemcpyHostToDevice)

  cuda_result = cudaMemcpy(ns_d, ns, 1, cudaMemcpyHostToDevice)
  cuda_result = cudaMemcpy(np_d, np,  1, cudaMemcpyHostToDevice)

  numBlocks  = dim3(256,1,1)
  numThreads = dim3(ceiling(real(np)/numBlocks%x),1,1)

  call deposit_kernel<<<numBlocks, numThreads>>>(particles, ns_d, np_d, &
       rho, lo_d, hi_d, &
       plo_d, dx_d)
  
end subroutine deposit_cic
