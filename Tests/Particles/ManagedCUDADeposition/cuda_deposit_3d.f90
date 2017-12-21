attributes(global) subroutine push_kernel(particles)
  
  use amrex_fort_module, only : amrex_real  
  implicit none
  
  integer              :: ns, np
  real(amrex_real)     :: particles(ns, np)
  integer i

  real(amrex_real) mass, charge, dt, fac
  mass = 1.d0
  charge = 1.d0
  dt = 1.d-6
  fac = charge*dt / mass

  i = blockDim%x * (blockIdx%x - 1) + threadIdx%x
  if (i <= np) then 

     particles(5, i) = particles(5, i) + fac*particles(8,  i)
     particles(6, i) = particles(6, i) + fac*particles(9,  i)
     particles(7, i) = particles(7, i) + fac*particles(10, i)

     particles(1, i) = particles(1, i) + dt * particles(5, i)
     particles(2, i) = particles(2, i) + dt * particles(6, i)
     particles(3, i) = particles(3, i) + dt * particles(7, i)

  end if
 
end subroutine push_kernel

subroutine cuda_push_particles(particles, ns, np) &
     bind(c,name='cuda_push_particles')
  
  use iso_c_binding
  use amrex_fort_module, only : amrex_real
  use cudafor, only: dim3, cudaMemcpy, cudaMemcpyHostToDevice
  
  implicit none
  
  integer, value   :: ns, np
  real(amrex_real) :: particles(ns, np)

  attributes(device) :: particles

  integer :: cuda_result
  integer, device    :: ns_d, np_d
  type(dim3) :: numThreads, numBlocks

  cuda_result = cudaMemcpy(ns_d, ns, 1, cudaMemcpyHostToDevice)
  cuda_result = cudaMemcpy(np_d, np, 1, cudaMemcpyHostToDevice)

  numThreads = dim3(256,1,1)
  numBlocks  = dim3(ceiling(real(np)/numThreads%x),1,1)

  call push_kernel<<<numBlocks, numThreads>>>(particles, ns_d, np_d)

end subroutine cuda_push_particles

attributes(global) subroutine deposit_kernel(particles, ns, np, &
                                             rho, lo, hi, plo, dx)
  
  use amrex_fort_module, only : amrex_real
  implicit none

  integer              :: ns, np
  real(amrex_real)     :: particles(ns, np)
  integer, intent(in)  :: lo(3), hi(3)
  real(amrex_real)     :: rho(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))
  real(amrex_real)     :: plo(3)
  real(amrex_real)     :: dx(3)
  
  integer i, j, k, n
  real(amrex_real) retval
  real(amrex_real) wx_lo, wy_lo, wz_lo, wx_hi, wy_hi, wz_hi
  real(amrex_real) lx, ly, lz
  real(amrex_real) inv_dx(3)

  lx = particles(10, n)

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

     retval = atomicadd(rho(i-1, j-1, k-1), wx_lo*wy_lo*wz_lo*particles(4, n))
     retval = atomicadd(rho(i-1, j-1, k),   wx_lo*wy_lo*wz_hi*particles(4, n))
     retval = atomicadd(rho(i-1, j,   k-1), wx_lo*wy_hi*wz_lo*particles(4, n))
     retval = atomicadd(rho(i-1, j,   k),   wx_lo*wy_hi*wz_hi*particles(4, n))
     retval = atomicadd(rho(i,   j-1, k-1), wx_hi*wy_lo*wz_lo*particles(4, n))
     retval = atomicadd(rho(i,   j-1, k),   wx_hi*wy_lo*wz_hi*particles(4, n))
     retval = atomicadd(rho(i,   j,   k-1), wx_hi*wy_hi*wz_lo*particles(4, n))
     retval = atomicadd(rho(i,   j,   k),   wx_hi*wy_hi*wz_hi*particles(4, n))

  end if
 
end subroutine deposit_kernel

subroutine cuda_deposit_cic(particles, ns, np, gid, &
                            rho, lo, hi, plo, dx) &
                            bind(c,name='cuda_deposit_cic')
  
  use iso_c_binding
  use amrex_fort_module, only : amrex_real
  use cudafor
  use cuda_module, only: cuda_streams, max_cuda_streams

  implicit none
  
  integer, value       :: ns, np, gid
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

  stream = cuda_streams(mod(gid, max_cuda_streams) + 1)

  cuda_result = cudaMemcpyAsync(lo_d, lo, 3, cudaMemcpyHostToDevice, stream)

  cuda_result = cudaMemcpyAsync(hi_d, hi, 3, cudaMemcpyHostToDevice, stream)

  cuda_result = cudaMemcpyAsync(plo_d, plo, 3, cudaMemcpyHostToDevice, stream)

  cuda_result = cudaMemcpyAsync(dx_d,  dx,  3, cudaMemcpyHostToDevice, stream)

  cuda_result = cudaMemcpyAsync(ns_d, ns, 1, cudaMemcpyHostToDevice, stream)

  cuda_result = cudaMemcpyAsync(np_d, np, 1, cudaMemcpyHostToDevice, stream)

  numThreads = dim3(256,1,1)
  numBlocks  = dim3(ceiling(real(np)/numThreads%x),1,1)

  call deposit_kernel<<<numBlocks, numThreads, 0, stream>>>(particles, ns_d, np_d, &
       rho, lo_d, hi_d, &
       plo_d, dx_d)

end subroutine cuda_deposit_cic

attributes(global) subroutine interpolate_kernel(particles, ns, np, &
                                                 acc, lo, hi, plo, dx)
  
  use amrex_fort_module, only : amrex_real  
  implicit none

  integer              :: ns, np
  real(amrex_real)     :: particles(ns, np)
  integer, intent(in)  :: lo(3), hi(3)
  real(amrex_real)     :: acc(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 3)
  real(amrex_real)     :: plo(3)
  real(amrex_real)     :: dx(3)
  
  integer i, j, k, n, nc
  real(amrex_real) retval
  real(amrex_real) wx_lo, wy_lo, wz_lo, wx_hi, wy_hi, wz_hi
  real(amrex_real) lx, ly, lz
  real(amrex_real) inv_dx(3)
  inv_dx = 1.0d0/dx

  n = blockDim%x * (blockIdx%x - 1) + threadIdx%x
  if (n <= np) then
     lx = (particles(1, n) - plo(1))*inv_dx(1) + 0.5d0
     ly = (particles(2, n) - plo(2))*inv_dx(2) + 0.5d0
     lz = (particles(3, n) - plo(3))*inv_dx(3) + 0.5d0

     lx = particles(10, n)

     i = floor(lx)
     j = floor(ly)
     k = floor(lz)
     
     wx_hi = lx - i
     wy_hi = ly - j
     wz_hi = lz - k
     
     wx_lo = 1.0d0 - wx_hi
     wy_lo = 1.0d0 - wy_hi
     wz_lo = 1.0d0 - wz_hi

     do nc = 1, 3
        particles(7+nc, n)  = wx_lo*wy_lo*wz_lo*acc(i-1, j-1, k-1, nc) + &
                              wx_lo*wy_lo*wz_hi*acc(i-1, j-1, k  , nc) + &
                              wx_lo*wy_hi*wz_lo*acc(i-1, j,   k-1, nc) + &
                              wx_lo*wy_hi*wz_hi*acc(i-1, j,   k  , nc) + &
                              wx_hi*wy_lo*wz_lo*acc(i,   j-1, k-1, nc) + &
                              wx_hi*wy_lo*wz_hi*acc(i,   j-1, k  , nc) + &
                              wx_hi*wy_hi*wz_lo*acc(i,   j,   k-1, nc) + &
                              wx_hi*wy_hi*wz_hi*acc(i,   j,   k  , nc)       
       end do
  end if
 
end subroutine interpolate_kernel

subroutine cuda_interpolate_cic(particles, ns, np, gid, &
                                acc, lo, hi, plo, dx)   &
                                bind(c,name='cuda_interpolate_cic')
  
  use iso_c_binding
  use amrex_fort_module, only : amrex_real
  use cudafor
  use cuda_module, only: cuda_streams, max_cuda_streams
  
  implicit none
  
  integer, value       :: ns, np, gid
  real(amrex_real)     :: particles(ns,np)
  integer              :: lo(3)
  integer              :: hi(3)
  real(amrex_real)     :: acc(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 3)
  real(amrex_real)     :: plo(3)
  real(amrex_real)     :: dx(3)
  
  attributes(device)         :: particles
  attributes(device)         :: acc
  integer,  device           :: lo_d(3), hi_d(3)
  real(amrex_real), device   :: plo_d(3), dx_d(3)
  integer, device            :: ns_d, np_d

  integer :: cuda_result
  integer(kind=cuda_stream_kind) :: stream
  type(dim3) :: numThreads, numBlocks

  stream = cuda_streams(mod(gid, max_cuda_streams) + 1)

  cuda_result = cudaMemcpyAsync(lo_d, lo, 3, cudaMemcpyHostToDevice, stream)
  cuda_result = cudaMemcpyAsync(hi_d, hi, 3, cudaMemcpyHostToDevice, stream)

  cuda_result = cudaMemcpyAsync(plo_d, plo, 3, cudaMemcpyHostToDevice, stream)
  cuda_result = cudaMemcpyAsync(dx_d,  dx,  3, cudaMemcpyHostToDevice, stream)

  cuda_result = cudaMemcpyAsync(ns_d,  ns,  1, cudaMemcpyHostToDevice, stream)
  cuda_result = cudaMemcpyAsync(np_d,  np,  1, cudaMemcpyHostToDevice, stream)

  numThreads = dim3(256,1,1)
  numBlocks  = dim3(ceiling(real(np)/numThreads%x),1,1)

  call interpolate_kernel<<<numBlocks, numThreads, 0, stream>>>(particles, ns_d, np_d, &
       acc, lo_d, hi_d, &
       plo_d, dx_d)

end subroutine cuda_interpolate_cic
