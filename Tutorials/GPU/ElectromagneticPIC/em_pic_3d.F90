module constants
  use amrex_fort_module, only : amrex_real
  
  real(amrex_real), parameter :: clight  = 2.99792458d8
  real(amrex_real), parameter :: epsilon_0 = 8.85418782d-12
  real(amrex_real), parameter :: electron_charge = 1.60217662d-19
  real(amrex_real), parameter :: electron_mass = 9.10938356d-31

end module

#ifdef AMREX_USE_ACC
subroutine push_momentum_boris(np, uxp, uyp, uzp, gaminv, &
     ex, ey, ez, bx, by, bz, q, m, dt)
#else
AMREX_CUDA_FORT_DEVICE subroutine push_momentum_boris(np, uxp, uyp, uzp, gaminv, &
     ex, ey, ez, bx, by, bz, q, m, dt) &
     bind(c,name='push_momentum_boris')
#endif

  use amrex_fort_module, only : amrex_real
  use constants, only : clight
  implicit none

  integer,          intent(in), value :: np
  real(amrex_real), intent(inout)     :: uxp(np), uyp(np), uzp(np), gaminv(np)
  real(amrex_real), intent(in)        :: ex(np), ey(np), ez(np)
  real(amrex_real), intent(in)        :: bx(np), by(np), bz(np)
  real(amrex_real), intent(in), value :: q, m, dt

  integer                         :: ip
  real(amrex_real)                :: const
  real(amrex_real)                :: clghtisq, usq, tsqi
  real(amrex_real)                :: tx, ty, tz
  real(amrex_real)                :: sx, sy, sz
  real(amrex_real)                :: uxppr, uyppr, uzppr
  real(amrex_real)                :: gaminvtmp

  const = q*dt*0.5d0/m
  clghtisq = 1.d0/clight**2

!#if defined (AMREX_USE_CUDA) && !defined(AMREX_USE_ACC)
!  ip = blockDim%x * (blockIdx%x - 1) + threadIdx%x
!  if (ip .le. np) then 
!#else
# ifdef AMREX_USE_ACC
!$acc parallel deviceptr(uxp, uyp, uzp, gaminv, ex, ey, ez, bx, by, bz)
!$acc loop gang vector
# endif
  do ip = 1, np
!#endif
    ! Push using the electric field
    uxp(ip) = uxp(ip) + ex(ip)*const
    uyp(ip) = uyp(ip) + ey(ip)*const
    uzp(ip) = uzp(ip) + ez(ip)*const

    ! Compute temporary gamma
    usq = (uxp(ip)**2 + uyp(ip)**2+ uzp(ip)**2)*clghtisq
    gaminvtmp = 1.d0/sqrt(1.d0 + usq)

    ! Magnetic rotation
    tx = gaminvtmp*bx(ip)*const
    ty = gaminvtmp*by(ip)*const
    tz = gaminvtmp*bz(ip)*const
    tsqi = 2.d0/(1.d0 + tx**2 + ty**2 + tz**2)
    sx = tx*tsqi
    sy = ty*tsqi
    sz = tz*tsqi
    uxppr = uxp(ip) + uyp(ip)*tz - uzp(ip)*ty
    uyppr = uyp(ip) + uzp(ip)*tx - uxp(ip)*tz
    uzppr = uzp(ip) + uxp(ip)*ty - uyp(ip)*tx
    uxp(ip) = uxp(ip) + uyppr*sz - uzppr*sy
    uyp(ip) = uyp(ip) + uzppr*sx - uxppr*sz
    uzp(ip) = uzp(ip) + uxppr*sy - uyppr*sx

    ! Push using the electric field
    uxp(ip) = uxp(ip) + ex(ip)*const
    uyp(ip) = uyp(ip) + ey(ip)*const
    uzp(ip) = uzp(ip) + ez(ip)*const

    ! Compute final gamma
    usq = (uxp(ip)**2 + uyp(ip)**2+ uzp(ip)**2)*clghtisq
    gaminv(ip) = 1.d0/sqrt(1.d0 + usq)

!#if defined (AMREX_USE_CUDA) && !defined(AMREX_USE_ACC)
!  end if
!#else
  end do
# ifdef AMREX_USE_ACC
!$acc end loop
!$acc end parallel
# endif
!#endif
  
end subroutine push_momentum_boris


#ifdef AMREX_USE_ACC
subroutine push_position_boris(np, xp, yp, zp, uxp, uyp, uzp, gaminv, dt)
#else
AMREX_CUDA_FORT_DEVICE subroutine push_position_boris(np, xp, yp, zp, uxp, uyp, uzp, gaminv, dt) &
     bind(c,name='push_position_boris')
#endif

  use amrex_fort_module, only : amrex_real  
  implicit none

  integer,          intent(in), value  :: np
  real(amrex_real), intent(inout)      :: xp(np), yp(np), zp(np)
  real(amrex_real), intent(in)         :: uxp(np), uyp(np), uzp(np), gaminv(np)
  real(amrex_real), intent(in), value  :: dt
  
  integer                              :: ip

!#if defined (AMREX_USE_CUDA) && !defined(AMREX_USE_ACC)
!  ip = blockDim%x * (blockIdx%x - 1) + threadIdx%x
!  if (ip .le. np) then 
!#else
# ifdef AMREX_USE_ACC
!$acc parallel deviceptr(xp, yp, zp, uxp, uyp, uzp, gaminv)
!$acc loop gang vector
# endif
   do ip = 1, np
!#endif

    xp(ip) = xp(ip) + uxp(ip)*gaminv(ip)*dt
    yp(ip) = yp(ip) + uyp(ip)*gaminv(ip)*dt
    zp(ip) = zp(ip) + uzp(ip)*gaminv(ip)*dt

!#if defined (AMREX_USE_CUDA) && !defined(AMREX_USE_ACC)
!  end if
!#else
  end do
# ifdef AMREX_USE_ACC
!$acc end loop
!$acc end parallel
# endif
!#endif

end subroutine push_position_boris

AMREX_CUDA_FORT_DEVICE subroutine set_gamma(np, uxp, uyp, uzp, gaminv) &
     bind(c,name='set_gamma')
  
  use amrex_fort_module, only : amrex_real
  use constants, only : clight 
  implicit none
  
  integer,          intent(in), value :: np
  real(amrex_real), intent(in)        :: uxp(np), uyp(np), uzp(np)
  real(amrex_real), intent(inout)     :: gaminv(np)
  
  integer                           :: ip
  real(amrex_real)                  :: clghtisq, usq
  clghtisq = 1.d0/clight**2

!#ifdef AMREX_USE_CUDA
!  ip = blockDim%x * (blockIdx%x - 1) + threadIdx%x
!  if (ip .le. np) then 
!#else
  do ip = 1, np
!#endif

    usq = (uxp(ip)**2 + uyp(ip)**2+ uzp(ip)**2)*clghtisq
    gaminv(ip) = 1.d0/sqrt(1.d0 + usq)

!#ifdef AMREX_USE_CUDA
!  end if
!#else
  end do
!#endif

end subroutine set_gamma

AMREX_CUDA_FORT_DEVICE subroutine enforce_periodic(np, xp, yp, zp, plo, phi) &
     bind(c,name='enforce_periodic')
  
  use amrex_fort_module, only : amrex_real
  implicit none
  
  integer,          intent(in), value      :: np
  real(amrex_real), intent(inout)   :: xp(np), yp(np), zp(np)
  real(amrex_real), intent(in)      :: plo(3), phi(3)
  
  integer                           :: ip
  real(amrex_real)                  :: domain_size(3)

  domain_size = phi - plo

!#ifdef AMREX_USE_CUDA
!  ip = blockDim%x * (blockIdx%x - 1) + threadIdx%x
!  if (ip .le. np) then 
!#else
  do ip = 1, np
!#endif
     
     if (xp(ip) .gt. phi(1)) then
        xp(ip) = xp(ip) - domain_size(1)
     else if (xp(ip) .lt. plo(1)) then
        xp(ip) = xp(ip) + domain_size(1)
     end if

     if (yp(ip) .gt. phi(2)) then
        yp(ip) = yp(ip) - domain_size(2)
     else if (yp(ip) .lt. plo(2)) then
        yp(ip) = yp(ip) + domain_size(2)
     end if

     if (zp(ip) .gt. phi(3)) then
        zp(ip) = zp(ip) - domain_size(3)
     else if (zp(ip) .lt. plo(3)) then
        zp(ip) = zp(ip) + domain_size(3)
     end if

!#ifdef AMREX_USE_CUDA
!  end if
!#else
  end do
!#endif

end subroutine enforce_periodic

#ifdef AMREX_USE_ACC
subroutine deposit_current(jx, jxlo, jxhi, jy, jylo, jyhi, jz, jzlo, jzhi, np, xp, yp, zp, & 
     uxp, uyp, uzp, gaminv, w, q, plo, dt, dx)
#else
AMREX_CUDA_FORT_DEVICE subroutine deposit_current(jx, jxlo, jxhi, jy, jylo, jyhi, jz, jzlo, jzhi, np, xp, yp, zp, & 
     uxp, uyp, uzp, gaminv, w, q, plo, dt, dx) & 
     bind(c,name="deposit_current")
#endif

  use amrex_fort_module, only : amrex_real
  use constants, only : clight
  implicit none
  
  integer,          intent(in), value    :: np
  integer,          intent(in)    :: jxlo(3), jxhi(3)
  integer,          intent(in)    :: jylo(3), jyhi(3)
  integer,          intent(in)    :: jzlo(3), jzhi(3)
  real(amrex_real), intent(inout) :: jx(jxlo(1):jxhi(1), jxlo(2):jxhi(2), jxlo(3):jxhi(3))
  real(amrex_real), intent(inout) :: jy(jylo(1):jyhi(1), jylo(2):jyhi(2), jylo(3):jyhi(3))
  real(amrex_real), intent(inout) :: jz(jzlo(1):jzhi(1), jzlo(2):jzhi(2), jzlo(3):jzhi(3))
  real(amrex_real), intent(in)    :: xp(np), yp(np), zp(np), uxp(np), uyp(np), uzp(np)
  real(amrex_real), intent(in)    :: w(np), gaminv(np)
  real(amrex_real), value         :: q, dt  
  real(amrex_real), intent(in)    :: dx(3), plo(3)
  
  real(amrex_real)                :: dxi, dyi, dzi, xint, yint, zint, retval
  real(amrex_real)                :: x, y, z, xmid, ymid, zmid, vx, vy, vz
  real(amrex_real)                :: invvol, dts2dx, dts2dy, dts2dz
  real(amrex_real)                :: wq, wqx, wqy, wqz, clightsq  
  real(amrex_real), dimension(2)  :: sx(0:1), sy(0:1), sz(0:1), sx0(0:1), sy0(0:1), sz0(0:1)
  real(amrex_real), parameter     :: onesixth=1.d0/6.d0, twothird=2.d0/3.d0
  integer                         :: j, k, l, j0, k0, l0, ip
  
  dxi = 1.d0/dx(1)
  dyi = 1.d0/dx(2)
  dzi = 1.d0/dx(3)
  invvol = dxi*dyi*dzi
  dts2dx = 0.5d0*dt*dxi
  dts2dy = 0.5d0*dt*dyi
  dts2dz = 0.5d0*dt*dzi
  clightsq = 1.d0/clight**2
  ! CD: Not needed - initialized in parallel loop. It must be done this way
  ! because firstprivate is not supported on a loop construct.
  !sx=0.d0; sy=0.d0; sz=0.d0;
  !sx0=0.d0;sy0=0.d0;sz0=0.d0;

!#if defined (AMREX_USE_CUDA) && !defined(AMREX_USE_ACC)
!  ip = blockDim%x * (blockIdx%x - 1) + threadIdx%x
!  if (ip .le. np) then 
!#else
# ifdef AMREX_USE_ACC
! CD: The "private" clause is placed on the loop construct (containing a
! vector schedule) to ensure that there is a private copy of the data
! per thread - see 2.9.10. If "private" was applied to the parallel
! construct then there would be a private copy per gang - see 2.5.10.
!$acc parallel deviceptr(jxlo, jxhi, jylo, jyhi, jzlo, jzhi, jx, jy, jz, xp, yp, zp, uxp, uyp, uzp, w, gaminv, dx, plo) 
!$acc loop gang vector private(sx(0:1), sy(0:1), sz(0:1), sx0(0:1), sy0(0:1), sz0(0:1))
# endif
  do ip = 1, np
!#endif

    ! --- computes position in grid units at (n+1)
    x = (xp(ip)-plo(1))*dxi
    y = (yp(ip)-plo(2))*dyi
    z = (zp(ip)-plo(3))*dzi

    ! Computes velocity
    vx = uxp(ip)*gaminv(ip)
    vy = uyp(ip)*gaminv(ip)
    vz = uzp(ip)*gaminv(ip)

    ! --- computes particles weights
    wq=q*w(ip)
    wqx=wq*invvol*vx
    wqy=wq*invvol*vy
    wqz=wq*invvol*vz

    ! Gets position in grid units at (n+1/2) for computing rho(n+1/2)
    xmid=x-dts2dx*vx
    ymid=y-dts2dy*vy
    zmid=z-dts2dz*vz

    ! --- finds node of cell containing particles for current positions
    j  = floor(xmid)
    k  = floor(ymid)
    l  = floor(zmid)
    j0 = floor(xmid-0.5d0)
    k0 = floor(ymid-0.5d0)
    l0 = floor(zmid-0.5d0)

    ! --- computes set of coefficients for node centered quantities
    xint = xmid-j
    yint = ymid-k
    zint = zmid-l
    sx( 0) = 1.d0-xint
    sx( 1) = xint
    sy( 0) = 1.d0-yint
    sy( 1) = yint
    sz( 0) = 1.d0-zint
    sz( 1) = zint

    ! --- computes set of coefficients for staggered quantities
    xint = xmid-j0-0.5d0
    yint = ymid-k0-0.5d0
    zint = zmid-l0-0.5d0
    sx0( 0) = 1.d0-xint
    sx0( 1) = xint
    sy0( 0) = 1.d0-yint
    sy0( 1) = yint
    sz0( 0) = 1.d0-zint
    sz0( 1) = zint

    ! --- add current contributions in the form rho(n+1/2)v(n+1/2)
#ifdef AMREX_USE_CUDA

# ifdef AMREX_USE_ACC
    !$acc atomic update
    jx(j0, k, l  )      = jx(j0, k, l      )  +   sx0(0)*sy(0)*sz(0)*wqx
    !$acc atomic update
    jx(j0+1, k, l  )    = jx(j0+1, k, l    )  +   sx0(1)*sy(0)*sz(0)*wqx
    !$acc atomic update
    jx(j0, k+1, l  )    = jx(j0, k+1, l    )  +   sx0(0)*sy(1)*sz(0)*wqx
    !$acc atomic update
    jx(j0+1, k+1, l  )  = jx(j0+1, k+1, l  )  +   sx0(1)*sy(1)*sz(0)*wqx
    !$acc atomic update
    jx(j0, k, l+1)      = jx(j0, k, l+1    )  +   sx0(0)*sy(0)*sz(1)*wqx
    !$acc atomic update
    jx(j0+1, k, l+1)    = jx(j0+1, k, l+1  )  +   sx0(1)*sy(0)*sz(1)*wqx
    !$acc atomic update
    jx(j0, k+1, l+1)    = jx(j0, k+1, l+1  )  +   sx0(0)*sy(1)*sz(1)*wqx
    !$acc atomic update
    jx(j0+1, k+1, l+1)  = jx(j0+1, k+1, l+1)  +   sx0(1)*sy(1)*sz(1)*wqx

    !$acc atomic update
    jy(j, k0, l  )      = jy(j, k0, l      )  +   sx(0)*sy0(0)*sz(0)*wqy
    !$acc atomic update
    jy(j+1, k0, l  )    = jy(j+1, k0, l    )  +   sx(1)*sy0(0)*sz(0)*wqy
    !$acc atomic update
    jy(j, k0+1, l  )    = jy(j, k0+1, l    )  +   sx(0)*sy0(1)*sz(0)*wqy
    !$acc atomic update
    jy(j+1, k0+1, l  )  = jy(j+1, k0+1, l  )  +   sx(1)*sy0(1)*sz(0)*wqy
    !$acc atomic update
    jy(j, k0, l+1)      = jy(j, k0, l+1    )  +   sx(0)*sy0(0)*sz(1)*wqy
    !$acc atomic update
    jy(j+1, k0, l+1)    = jy(j+1, k0, l+1  )  +   sx(1)*sy0(0)*sz(1)*wqy
    !$acc atomic update
    jy(j, k0+1, l+1)    = jy(j, k0+1, l+1  )  +   sx(0)*sy0(1)*sz(1)*wqy
    !$acc atomic update
    jy(j+1, k0+1, l+1)  = jy(j+1, k0+1, l+1)  +   sx(1)*sy0(1)*sz(1)*wqy

    !$acc atomic update
    jz(j, k, l0  )      = jz(j, k, l0      )  +   sx(0)*sy(0)*sz0(0)*wqz
    !$acc atomic update
    jz(j+1, k, l0  )    = jz(j+1, k, l0    )  +   sx(1)*sy(0)*sz0(0)*wqz
    !$acc atomic update
    jz(j, k+1, l0  )    = jz(j, k+1, l0    )  +   sx(0)*sy(1)*sz0(0)*wqz
    !$acc atomic update
    jz(j+1, k+1, l0  )  = jz(j+1, k+1, l0  )  +   sx(1)*sy(1)*sz0(0)*wqz
    !$acc atomic update
    jz(j, k, l0+1)      = jz(j, k, l0+1    )  +   sx(0)*sy(0)*sz0(1)*wqz
    !$acc atomic update
    jz(j+1, k, l0+1)    = jz(j+1, k, l0+1  )  +   sx(1)*sy(0)*sz0(1)*wqz
    !$acc atomic update
    jz(j, k+1, l0+1)    = jz(j, k+1, l0+1  )  +   sx(0)*sy(1)*sz0(1)*wqz
    !$acc atomic update
    jz(j+1, k+1, l0+1)  = jz(j+1, k+1, l0+1)  +   sx(1)*sy(1)*sz0(1)*wqz
# else
    ! https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#atomicadd
    retval = atomicadd(jx(j0  , k  , l  ), sx0(0)*sy(0)*sz(0)*wqx)
    retval = atomicadd(jx(j0+1, k  , l  ), sx0(1)*sy(0)*sz(0)*wqx)
    retval = atomicadd(jx(j0  , k+1, l  ), sx0(0)*sy(1)*sz(0)*wqx)
    retval = atomicadd(jx(j0+1, k+1, l  ), sx0(1)*sy(1)*sz(0)*wqx)
    retval = atomicadd(jx(j0  , k  , l+1), sx0(0)*sy(0)*sz(1)*wqx)
    retval = atomicadd(jx(j0+1, k  , l+1), sx0(1)*sy(0)*sz(1)*wqx)
    retval = atomicadd(jx(j0  , k+1, l+1), sx0(0)*sy(1)*sz(1)*wqx)
    retval = atomicadd(jx(j0+1, k+1, l+1), sx0(1)*sy(1)*sz(1)*wqx)

    retval = atomicadd(jy(j  , k0  , l  ), sx(0)*sy0(0)*sz(0)*wqy)
    retval = atomicadd(jy(j+1, k0  , l  ), sx(1)*sy0(0)*sz(0)*wqy)
    retval = atomicadd(jy(j,   k0+1, l  ), sx(0)*sy0(1)*sz(0)*wqy)
    retval = atomicadd(jy(j+1, k0+1, l  ), sx(1)*sy0(1)*sz(0)*wqy)
    retval = atomicadd(jy(j  , k0  , l+1), sx(0)*sy0(0)*sz(1)*wqy)
    retval = atomicadd(jy(j+1, k0  , l+1), sx(1)*sy0(0)*sz(1)*wqy)
    retval = atomicadd(jy(j  , k0+1, l+1), sx(0)*sy0(1)*sz(1)*wqy)
    retval = atomicadd(jy(j+1, k0+1, l+1), sx(1)*sy0(1)*sz(1)*wqy)

    retval = atomicadd(jz(j  , k  , l0  ), sx(0)*sy(0)*sz0(0)*wqz)
    retval = atomicadd(jz(j+1, k  , l0  ), sx(1)*sy(0)*sz0(0)*wqz)
    retval = atomicadd(jz(j  , k+1, l0  ), sx(0)*sy(1)*sz0(0)*wqz)
    retval = atomicadd(jz(j+1, k+1, l0  ), sx(1)*sy(1)*sz0(0)*wqz)
    retval = atomicadd(jz(j  , k  , l0+1), sx(0)*sy(0)*sz0(1)*wqz)
    retval = atomicadd(jz(j+1, k  , l0+1), sx(1)*sy(0)*sz0(1)*wqz)
    retval = atomicadd(jz(j  , k+1, l0+1), sx(0)*sy(1)*sz0(1)*wqz)
    retval = atomicadd(jz(j+1, k+1, l0+1), sx(1)*sy(1)*sz0(1)*wqz)
# endif

#else

    jx(j0, k, l  )      = jx(j0, k, l      )  +   sx0(0)*sy(0)*sz(0)*wqx
    jx(j0+1, k, l  )    = jx(j0+1, k, l    )  +   sx0(1)*sy(0)*sz(0)*wqx
    jx(j0, k+1, l  )    = jx(j0, k+1, l    )  +   sx0(0)*sy(1)*sz(0)*wqx
    jx(j0+1, k+1, l  )  = jx(j0+1, k+1, l  )  +   sx0(1)*sy(1)*sz(0)*wqx
    jx(j0, k, l+1)      = jx(j0, k, l+1    )  +   sx0(0)*sy(0)*sz(1)*wqx
    jx(j0+1, k, l+1)    = jx(j0+1, k, l+1  )  +   sx0(1)*sy(0)*sz(1)*wqx
    jx(j0, k+1, l+1)    = jx(j0, k+1, l+1  )  +   sx0(0)*sy(1)*sz(1)*wqx
    jx(j0+1, k+1, l+1)  = jx(j0+1, k+1, l+1)  +   sx0(1)*sy(1)*sz(1)*wqx

    jy(j, k0, l  )      = jy(j, k0, l      )  +   sx(0)*sy0(0)*sz(0)*wqy
    jy(j+1, k0, l  )    = jy(j+1, k0, l    )  +   sx(1)*sy0(0)*sz(0)*wqy
    jy(j, k0+1, l  )    = jy(j, k0+1, l    )  +   sx(0)*sy0(1)*sz(0)*wqy
    jy(j+1, k0+1, l  )  = jy(j+1, k0+1, l  )  +   sx(1)*sy0(1)*sz(0)*wqy
    jy(j, k0, l+1)      = jy(j, k0, l+1    )  +   sx(0)*sy0(0)*sz(1)*wqy
    jy(j+1, k0, l+1)    = jy(j+1, k0, l+1  )  +   sx(1)*sy0(0)*sz(1)*wqy
    jy(j, k0+1, l+1)    = jy(j, k0+1, l+1  )  +   sx(0)*sy0(1)*sz(1)*wqy
    jy(j+1, k0+1, l+1)  = jy(j+1, k0+1, l+1)  +   sx(1)*sy0(1)*sz(1)*wqy

    jz(j, k, l0  )      = jz(j, k, l0      )  +   sx(0)*sy(0)*sz0(0)*wqz
    jz(j+1, k, l0  )    = jz(j+1, k, l0    )  +   sx(1)*sy(0)*sz0(0)*wqz
    jz(j, k+1, l0  )    = jz(j, k+1, l0    )  +   sx(0)*sy(1)*sz0(0)*wqz
    jz(j+1, k+1, l0  )  = jz(j+1, k+1, l0  )  +   sx(1)*sy(1)*sz0(0)*wqz
    jz(j, k, l0+1)      = jz(j, k, l0+1    )  +   sx(0)*sy(0)*sz0(1)*wqz
    jz(j+1, k, l0+1)    = jz(j+1, k, l0+1  )  +   sx(1)*sy(0)*sz0(1)*wqz
    jz(j, k+1, l0+1)    = jz(j, k+1, l0+1  )  +   sx(0)*sy(1)*sz0(1)*wqz
    jz(j+1, k+1, l0+1)  = jz(j+1, k+1, l0+1)  +   sx(1)*sy(1)*sz0(1)*wqz

#endif

!#if defined (AMREX_USE_CUDA) && !defined(AMREX_USE_ACC)
!  end if
!#else
  end do
# ifdef AMREX_USE_ACC
!$acc end loop
!$acc end parallel
# endif
!#endif

end subroutine deposit_current

#ifdef AMREX_USE_ACC
subroutine gather_magnetic_field(np, xp, yp, zp, bx, by, bz, &
     bxg, bxglo, bxghi, byg, byglo, byghi, bzg, bzglo, bzghi, plo, dx)
#else
AMREX_CUDA_FORT_DEVICE subroutine gather_magnetic_field(np, xp, yp, zp, bx, by, bz, &
     bxg, bxglo, bxghi, byg, byglo, byghi, bzg, bzglo, bzghi, plo, dx) & 
     bind(c,name="gather_magnetic_field")
#endif

  use amrex_fort_module, only : amrex_real
  implicit none
  
  integer,          intent(in), value    :: np
  integer,          intent(in)    :: bxglo(3), bxghi(3)
  integer,          intent(in)    :: byglo(3), byghi(3)
  integer,          intent(in)    :: bzglo(3), bzghi(3)
  real(amrex_real), intent(in)    :: bxg(bxglo(1):bxghi(1), bxglo(2):bxghi(2), bxglo(3):bxghi(3))
  real(amrex_real), intent(in)    :: byg(byglo(1):byghi(1), byglo(2):byghi(2), byglo(3):byghi(3))
  real(amrex_real), intent(in)    :: bzg(bzglo(1):bzghi(1), bzglo(2):bzghi(2), bzglo(3):bzghi(3))
  real(amrex_real), intent(in)    :: xp(np), yp(np), zp(np)
  real(amrex_real), intent(inout) :: bx(np), by(np), bz(np)
  real(amrex_real), intent(in)    :: dx(3), plo(3)
  
  real(amrex_real)                :: x, y, z, dxi, dyi, dzi, xint, yint, zint
  real(amrex_real), dimension(2)  :: sx(0:1), sy(0:1), sz(0:1), sx0(0:1), sy0(0:1), sz0(0:1)
  real(amrex_real), parameter     :: onesixth=1.d0/6.d0, twothird=2.d0/3.d0
  integer                         :: j, k, l, ip, jj, kk, ll, j0, k0, l0
  integer                         :: ixmin, ixmax, iymin, iymax, izmin, izmax
  integer                         :: ixmin0, ixmax0, iymin0, iymax0, izmin0, izmax0

  dxi = 1.d0 / dx(1)
  dyi = 1.d0 / dx(2)
  dzi = 1.d0 / dx(3)

  ! CD: Not needed - initialized in parallel loop
  !sx  = 0.d0
  !sy  = 0.d0
  !sz  = 0.d0
  !sx0 = 0.d0
  !sy0 = 0.d0
  !sz0 = 0.d0

  ixmin = 0
  ixmax = 0
  iymin = 0
  iymax = 0
  izmin = 0
  izmax = 0

  ixmin0 = 0
  ixmax0 = 0
  iymin0 = 0
  iymax0 = 0
  izmin0 = 0
  izmax0 = 0

!#if defined (AMREX_USE_CUDA) && !defined(AMREX_USE_ACC)
!  ip = blockDim%x * (blockIdx%x - 1) + threadIdx%x
!  if (ip .le. np) then 
!#else
# ifdef AMREX_USE_ACC
!$acc parallel deviceptr(bxglo, bxghi, byglo, byghi, bzglo, bzghi, bxg, byg, bzg, xp, yp, zp, bx, by, bz, dx, plo)
!$acc loop gang vector private(sx(0:1), sy(0:1), sz(0:1), sx0(0:1), sy0(0:1), sz0(0:1))
# endif
  do ip = 1, np
!#endif
     
     bx(ip) = 0.d0
     by(ip) = 0.d0
     bz(ip) = 0.d0

     x = (xp(ip)-plo(1))*dxi
     y = (yp(ip)-plo(2))*dyi
     z = (zp(ip)-plo(3))*dzi
     
     ! Compute index of particle
     j  = floor(x)
     j0 = floor(x)
     k  = floor(y)
     k0 = floor(y)
     l  = floor(z)
     l0 = floor(z)
     
     xint = x - j
     yint = y - k
     zint = z - l
     
     ! Compute shape factors
     sx(0) = 1.d0-xint
     sx(1) = xint
     sy(0) = 1.d0-yint
     sy(1) = yint
     sz(0) = 1.d0-zint
     sz(1) = zint
     
     xint=x-0.5d0-j0
     yint=y-0.5d0-k0
     zint=z-0.5d0-l0
     
     sx0(0) = 1.d0
     sy0(0) = 1.d0
     sz0(0) = 1.d0
     sx0(1) = 0.d0 ! Added by CD
     sy0(1) = 0.d0 ! Added by CD
     sz0(1) = 0.d0 ! Added by CD
     
#ifdef AMREX_USE_ACC
! CD: The independent clause silences the compiler about the loop
! carried dependence of bx. There is no loop carried dependence
! because each bx index is only ever accessed by 1 thread.
!$acc loop seq independent collapse(3)
#endif
     do ll = izmin0, izmax0
        do kk = iymin0, iymax0
           do jj = ixmin, ixmax+1
              bx(ip) = bx(ip) + sx(jj)*sy0(kk)*sz0(ll)*bxg(j+jj, k0+kk, l0+ll)
           end do
        end do
     end do
#ifdef AMREX_USE_ACC
!$acc end loop    
#endif

#ifdef AMREX_USE_ACC
!$acc loop seq independent collapse(3)
#endif
     do ll = izmin0, izmax0
        do kk = iymin, iymax+1
           do jj = ixmin0, ixmax0
              by(ip) = by(ip) + sx0(jj)*sy(kk)*sz0(ll)*byg(j0+jj, k+kk, l0+ll)
           end do
        end do
     end do
#ifdef AMREX_USE_ACC
!$acc end loop
#endif

#ifdef AMREX_USE_ACC
!$acc loop seq independent collapse(3)
#endif
     do ll = izmin, izmax+1
        do kk = iymin0, iymax0
           do jj = ixmin0, ixmax0
              bz(ip) = bz(ip) + sx0(jj)*sy0(kk)*sz(ll)*bzg(j0+jj, k0+kk, l+ll)
           end do
        end do
     end do
#ifdef AMREX_USE_ACC
!$acc end loop
#endif

!#if defined (AMREX_USE_CUDA) && !defined(AMREX_USE_ACC)
!  end if
!#else
  end do
# ifdef AMREX_USE_ACC
!$acc end loop
!$acc end parallel
# endif
!#endif

end subroutine gather_magnetic_field

#ifdef AMREX_USE_ACC
subroutine gather_electric_field(np, xp, yp, zp, ex, ey, ez, &
     exg, exglo, exghi, eyg, eyglo, eyghi, ezg, ezglo, ezghi, plo, dx)
#else
AMREX_CUDA_FORT_DEVICE subroutine gather_electric_field(np, xp, yp, zp, ex, ey, ez, &
     exg, exglo, exghi, eyg, eyglo, eyghi, ezg, ezglo, ezghi, plo, dx) & 
     bind(c,name="gather_electric_field")
#endif

  use amrex_fort_module, only : amrex_real
  use constants, only : clight
  implicit none
  
  integer,          intent(in), value    :: np
  integer,          intent(in)    :: exglo(3), exghi(3)
  integer,          intent(in)    :: eyglo(3), eyghi(3)
  integer,          intent(in)    :: ezglo(3), ezghi(3)
  real(amrex_real), intent(in)    :: exg(exglo(1):exghi(1), exglo(2):exghi(2), exglo(3):exghi(3))
  real(amrex_real), intent(in)    :: eyg(eyglo(1):eyghi(1), eyglo(2):eyghi(2), eyglo(3):eyghi(3))
  real(amrex_real), intent(in)    :: ezg(ezglo(1):ezghi(1), ezglo(2):ezghi(2), ezglo(3):ezghi(3))
  real(amrex_real), intent(in)    :: xp(np), yp(np), zp(np)
  real(amrex_real), intent(inout) :: ex(np), ey(np), ez(np)
  real(amrex_real), intent(in)    :: dx(3), plo(3)
  
  real(amrex_real)                :: x, y, z, dxi, dyi, dzi, xint, yint, zint
  real(amrex_real), dimension(2)  :: sx(0:1), sy(0:1), sz(0:1), sx0(0:1), sy0(0:1), sz0(0:1)
  real(amrex_real), parameter     :: onesixth=1.d0/6.d0, twothird=2.d0/3.d0
  integer                         :: j, k, l, ip, jj, kk, ll, j0, k0, l0
  integer                         :: ixmin, ixmax, iymin, iymax, izmin, izmax
  integer                         :: ixmin0, ixmax0, iymin0, iymax0, izmin0, izmax0

  dxi = 1.d0 / dx(1)
  dyi = 1.d0 / dx(2)
  dzi = 1.d0 / dx(3)

  ixmin = 0
  ixmax = 0
  iymin = 0
  iymax = 0
  izmin = 0
  izmax = 0

  ! CD: Not needed - initialized in parallel loop
  !sx  = 0.d0
  !sy  = 0.d0
  !sz  = 0.d0
  !sx0 = 0.d0
  !sy0 = 0.d0
  !sz0 = 0.d0

  ixmin0 = 0
  ixmax0 = 0
  iymin0 = 0
  iymax0 = 0
  izmin0 = 0
  izmax0 = 0

!#if defined (AMREX_USE_CUDA) && !defined(AMREX_USE_ACC)
!  ip = blockDim%x * (blockIdx%x - 1) + threadIdx%x
!  if (ip .le. np) then 
!#else
# ifdef AMREX_USE_ACC
!$acc parallel deviceptr(exglo, exghi, eyglo, eyghi, ezglo, ezghi, exg, eyg, ezg, xp, yp, zp, ex, ey, ez, dx, plo)
!$acc loop gang vector private(sx(0:1), sy(0:1), sz(0:1), sx0(0:1), sy0(0:1), sz0(0:1))
# endif
  do ip = 1, np
!#endif
     
     ex(ip) = 0.d0
     ey(ip) = 0.d0
     ez(ip) = 0.d0

     x = (xp(ip)-plo(1))*dxi
     y = (yp(ip)-plo(2))*dyi
     z = (zp(ip)-plo(3))*dzi
     
     ! Compute index of particle
     j  = floor(x)
     j0 = floor(x)
     k  = floor(y)
     k0 = floor(y)
     l  = floor(z)
     l0 = floor(z)
     
     xint = x - j
     yint = y - k
     zint = z - l
     
     ! Compute shape factors
     sx(0) = 1.d0-xint
     sx(1) = xint
     sy(0) = 1.d0-yint
     sy(1) = yint
     sz(0) = 1.d0-zint
     sz(1) = zint
     
     xint=x-0.5d0-j0
     yint=y-0.5d0-k0
     zint=z-0.5d0-l0
     
     sx0(0) = 1.d0
     sy0(0) = 1.d0
     sz0(0) = 1.d0
     sx0(1) = 0.d0 ! Added by CD
     sy0(1) = 0.d0 ! Added by CD
     sz0(1) = 0.d0 ! Added by CD
     
#ifdef AMREX_USE_ACC
!$acc loop seq independent collapse(3)
#endif
     do ll = izmin, izmax+1
        do kk = iymin, iymax+1
           do jj = ixmin0, ixmax0
              ex(ip) = ex(ip) + sx0(jj)*sy(kk)*sz(ll)*exg(j0+jj, k+kk, l+ll)
           end do
        end do
     end do
#ifdef AMREX_USE_ACC
!$acc end loop    
#endif

#ifdef AMREX_USE_ACC
!$acc loop seq independent collapse(3)
#endif
     do ll = izmin, izmax+1
        do kk = iymin0, iymax0
           do jj = ixmin, ixmax+1
              ey(ip) = ey(ip) + sx(jj)*sy0(kk)*sz(ll)*eyg(j+jj, k0+kk, l+ll)
           end do
        end do
     end do
#ifdef AMREX_USE_ACC
!$acc end loop    
#endif
     
#ifdef AMREX_USE_ACC
!$acc loop seq independent collapse(3)
#endif
     do ll = izmin0, izmax0
        do kk = iymin, iymax+1
           do jj = ixmin, ixmax+1
              ez(ip) = ez(ip) + sx(jj)*sy(kk)*sz0(ll)*ezg(j+jj, k+kk, l0+ll)
           end do
        end do
     end do
#ifdef AMREX_USE_ACC
!$acc end loop    
#endif

!#if defined (AMREX_USE_CUDA) && !defined(AMREX_USE_ACC)
!  end if
!#else
  end do
# ifdef AMREX_USE_ACC
!$acc end loop
!$acc end parallel
# endif
!#endif
   
end subroutine gather_electric_field

AMREX_CUDA_FORT_DEVICE subroutine push_electric_field_x(xlo, xhi, ex, exlo, exhi,               &
     by, bylo, byhi, bz, bzlo, bzhi, jx, jxlo, jxhi, mudt, dtsdy, dtsdz) & 
     bind(c,name="push_electric_field_x")
  
  use amrex_fort_module, only : amrex_real, get_loop_bounds_device
  implicit none
  
  integer,          intent(in)    :: xlo(3),  xhi(3)
  integer,          intent(in)    :: exlo(3),exhi(3)
  integer,          intent(in)    :: bylo(3),byhi(3),bzlo(3),bzhi(3)
  integer,          intent(in)    :: jxlo(3),jxhi(3)
  real(amrex_real), intent(inout) :: ex(exlo(1):exhi(1),exlo(2):exhi(2),exlo(3):exhi(3))
  real(amrex_real), intent(in)    :: by(bylo(1):byhi(1),bylo(2):byhi(2),bylo(3):byhi(3))
  real(amrex_real), intent(in)    :: bz(bzlo(1):bzhi(1),bzlo(2):bzhi(2),bzlo(3):bzhi(3))
  real(amrex_real), intent(in)    :: jx(jxlo(1):jxhi(1),jxlo(2):jxhi(2),jxlo(3):jxhi(3))
  real(amrex_real), value         :: mudt,dtsdy,dtsdz
  
  integer :: j,k,l
  integer :: blo(3), bhi(3)

!  call get_loop_bounds_device(blo, bhi, xlo, xhi)

  do l       = xlo(3), xhi(3)
     do k    = xlo(2), xhi(2)
        do j = xlo(1), xhi(1)
           Ex(j,k,l) = Ex(j,k,l) + dtsdy * (Bz(j,k,l) - Bz(j,k-1,l  )) &
                - dtsdz * (By(j,k,l) - By(j,k  ,l-1)) &
                - mudt  * jx(j,k,l)
        end do
     end do
  end do
    
end subroutine push_electric_field_x

AMREX_CUDA_FORT_DEVICE subroutine push_electric_field_y(ylo, yhi, &
     ey, eylo, eyhi, bx, bxlo, bxhi, bz, bzlo, bzhi, &
     jy, jylo, jyhi, mudt, dtsdx, dtsdz)     & 
     bind(c,name="push_electric_field_y")
  
  use amrex_fort_module, only : amrex_real, get_loop_bounds_device
  implicit none
  
  integer,          intent(in)    :: ylo(3), yhi(3)
  integer,          intent(in)    :: eylo(3),eyhi(3)
  integer,          intent(in)    :: bxlo(3),bxhi(3), bzlo(3),bzhi(3)
  integer,          intent(in)    :: jylo(3),jyhi(3)
  real(amrex_real), intent(inout) :: ey(eylo(1):eyhi(1),eylo(2):eyhi(2),eylo(3):eyhi(3))
  real(amrex_real), intent(in)    :: bx(bxlo(1):bxhi(1),bxlo(2):bxhi(2),bxlo(3):bxhi(3))
  real(amrex_real), intent(in)    :: bz(bzlo(1):bzhi(1),bzlo(2):bzhi(2),bzlo(3):bzhi(3))
  real(amrex_real), intent(in)    :: jy(jylo(1):jyhi(1),jylo(2):jyhi(2),jylo(3):jyhi(3))
  real(amrex_real), value         :: mudt,dtsdx,dtsdz
  
  integer :: j,k,l
  integer :: blo(3), bhi(3)

!  call get_loop_bounds_device(blo, bhi, ylo, yhi)

  do l       = ylo(3), yhi(3)
     do k    = ylo(2), yhi(2)
        do j = ylo(1), yhi(1)
           Ey(j,k,l) = Ey(j,k,l) - dtsdx * (Bz(j,k,l) - Bz(j-1,k,l)) &
                + dtsdz * (Bx(j,k,l) - Bx(j,k,l-1)) &
                - mudt  * jy(j,k,l)
        end do
     end do
  end do

end subroutine push_electric_field_y

AMREX_CUDA_FORT_DEVICE subroutine push_electric_field_z(zlo, zhi, &
     ez,ezlo, ezhi, bx, bxlo, bxhi, by, bylo, byhi, &
     jz, jzlo, jzhi, mudt, dtsdx, dtsdy)     & 
     bind(c,name="push_electric_field_z")
  
  use amrex_fort_module, only : amrex_real, get_loop_bounds_device
  implicit none
  
  integer,          intent(in)    :: zlo(3), zhi(3)
  integer,          intent(in)    :: ezlo(3),ezhi(3)
  integer,          intent(in)    :: bxlo(3),bxhi(3),bylo(3),byhi(3)
  integer,          intent(in)    :: jzlo(3),jzhi(3)
  real(amrex_real), intent(inout) :: ez(ezlo(1):ezhi(1),ezlo(2):ezhi(2),ezlo(3):ezhi(3))
  real(amrex_real), intent(in)    :: bx(bxlo(1):bxhi(1),bxlo(2):bxhi(2),bxlo(3):bxhi(3))
  real(amrex_real), intent(in)    :: by(bylo(1):byhi(1),bylo(2):byhi(2),bylo(3):byhi(3))
  real(amrex_real), intent(in)    :: jz(jzlo(1):jzhi(1),jzlo(2):jzhi(2),jzlo(3):jzhi(3))
  real(amrex_real), value         :: mudt,dtsdx,dtsdy
  
  integer :: j,k,l
  integer :: blo(3), bhi(3)

!  call get_loop_bounds_device(blo, bhi, zlo, zhi)

  do l       = zlo(3), zhi(3)
     do k    = zlo(2), zhi(2)
        do j = zlo(1), zhi(1)
           Ez(j,k,l) = Ez(j,k,l) + dtsdx * (By(j,k,l) - By(j-1,k  ,l)) &
                - dtsdy * (Bx(j,k,l) - Bx(j  ,k-1,l)) &
                - mudt  * jz(j,k,l)
        end do
     end do
  end do
  
end subroutine push_electric_field_z

AMREX_CUDA_FORT_DEVICE subroutine push_magnetic_field_x(xlo, xhi, bx, bxlo, bxhi, ey, eylo, eyhi, & 
     ez, ezlo, ezhi, dtsdy, dtsdz)    & 
     bind(c,name='push_magnetic_field_x')
  
  use amrex_fort_module, only : amrex_real, get_loop_bounds_device
  implicit none
  
  integer,          intent(in)    :: xlo(3), xhi(3)
  integer,          intent(in)    :: bxlo(3),bxhi(3)
  integer,          intent(in)    :: eylo(3),eyhi(3),ezlo(3),ezhi(3)
  real(amrex_real), intent(inout) :: bx(bxlo(1):bxhi(1),bxlo(2):bxhi(2),bxlo(3):bxhi(3))
  real(amrex_real), intent(in)    :: ey(eylo(1):eyhi(1),eylo(2):eyhi(2),eylo(3):eyhi(3))
  real(amrex_real), intent(in)    :: ez(ezlo(1):ezhi(1),ezlo(2):ezhi(2),ezlo(3):ezhi(3))
  real(amrex_real), value         :: dtsdy,dtsdz
  
  integer :: j,k,l
  integer :: blo(3), bhi(3)

!  call get_loop_bounds_device(blo, bhi, xlo, xhi)

  do l       = xlo(3), xhi(3)
     do k    = xlo(2), xhi(2)
        do j = xlo(1), xhi(1)
           Bx(j,k,l) = Bx(j,k,l) - dtsdy *(Ez(j, k+1, l) - Ez(j,k,l)) &
                + dtsdz *(Ey(j  ,k  ,l+1) - Ey(j,k,l))
        end do
     end do
  end do

end subroutine push_magnetic_field_x
  
AMREX_CUDA_FORT_DEVICE subroutine push_magnetic_field_y(ylo, yhi, by, bylo, byhi, ex, exlo, exhi, & 
     ez, ezlo, ezhi, dtsdx, dtsdz)    & 
     bind(c,name='push_magnetic_field_y')
  
  use amrex_fort_module, only : amrex_real, get_loop_bounds_device
  implicit none
  
  integer,          intent(in)    :: ylo(3), yhi(3)
  integer,          intent(in)    :: bylo(3),byhi(3)
  integer,          intent(in)    :: exlo(3),exhi(3),ezlo(3),ezhi(3)
  real(amrex_real), intent(inout) :: by(bylo(1):byhi(1),bylo(2):byhi(2),bylo(3):byhi(3))
  real(amrex_real), intent(in)    :: ex(exlo(1):exhi(1),exlo(2):exhi(2),exlo(3):exhi(3))
  real(amrex_real), intent(in)    :: ez(ezlo(1):ezhi(1),ezlo(2):ezhi(2),ezlo(3):ezhi(3))
  real(amrex_real), value         :: dtsdx,dtsdz
  
  integer :: j,k,l
  integer :: blo(3), bhi(3)

!  call get_loop_bounds_device(blo, bhi, ylo, yhi)
  
  do l       = ylo(3), yhi(3)
     do k    = ylo(2), yhi(2)
        do j = ylo(1), yhi(1)
           By(j,k,l) = By(j,k,l) + dtsdx * (Ez(j+1,k  ,l  ) - Ez(j,k,l)) &
                - dtsdz * (Ex(j  ,k  ,l+1) - Ex(j,k,l))
        end do
     end do
  end do

end subroutine push_magnetic_field_y

AMREX_CUDA_FORT_DEVICE subroutine push_magnetic_field_z(zlo, zhi, bz, bzlo, bzhi, ex, exlo, exhi, & 
     ey, eylo, eyhi, dtsdx, dtsdy)    & 
     bind(c,name='push_magnetic_field_z')
  
  use amrex_fort_module, only : amrex_real, get_loop_bounds_device
  implicit none
  
  integer,          intent(in)    :: zlo(3), zhi(3)
  integer,          intent(in)    :: bzlo(3),bzhi(3)
  integer,          intent(in)    :: exlo(3),exhi(3),eylo(3),eyhi(3)
  real(amrex_real), intent(inout) :: bz(bzlo(1):bzhi(1),bzlo(2):bzhi(2),bzlo(3):bzhi(3))
  real(amrex_real), intent(in)    :: ex(exlo(1):exhi(1),exlo(2):exhi(2),exlo(3):exhi(3))
  real(amrex_real), intent(in)    :: ey(eylo(1):eyhi(1),eylo(2):eyhi(2),eylo(3):eyhi(3))
  real(amrex_real), value         :: dtsdx,dtsdy
  
  integer :: j,k,l
  integer :: blo(3), bhi(3)
  
!  call get_loop_bounds_device(blo, bhi, zlo, zhi)
  
  do l       = zlo(3), zhi(3)
     do k    = zlo(2), zhi(2)
        do j = zlo(1), zhi(1)
           Bz(j,k,l) = Bz(j,k,l) - dtsdx * (Ey(j+1,k  ,l  ) - Ey(j,k,l)) &
                + dtsdy * (Ex(j  ,k+1,l  ) - Ex(j,k,l))
        end do
     end do
  end do
  
end subroutine push_magnetic_field_z

subroutine check_langmuir_solution(boxlo, boxhi, testlo, testhi, jx, jxlo, jxhi, & 
     time, max_error) bind(c,name='check_langmuir_solution')
  
  use amrex_fort_module, only : amrex_real
  use constants

  implicit none
  
  integer,          intent(in)    :: boxlo(3),  boxhi(3)
  integer,          intent(in)    :: testlo(3), testhi(3)
  integer,          intent(in)    :: jxlo(3),   jxhi(3)
  real(amrex_real), intent(in)    :: jx(jxlo(1):jxhi(1),jxlo(2):jxhi(2),jxlo(3):jxhi(3))
  real(amrex_real), intent(in), value :: time
  real(amrex_real), intent(inout) :: max_error

  integer :: j,k,l
  integer :: lo(3), hi(3)

  real(amrex_real) error
  real(amrex_real) exact

  real(amrex_real), parameter :: u  = 0.01
  real(amrex_real), parameter :: n0 = 1.d25
  real(amrex_real) wp

  wp = sqrt(n0*electron_charge**2/(electron_mass*epsilon_0))

  error = 0.0
  exact = -n0*electron_charge*clight*u*cos(wp*time)

  lo = max(boxlo, testlo)
  hi = min(boxhi, testhi)

  do l       = lo(3), hi(3)
     do k    = lo(2), hi(2)
        do j = lo(1), hi(1)
           error = max(error, abs(jx(j,k,l) - exact) / abs(exact))
        end do
     end do
  end do

  max_error = error
  
end subroutine check_langmuir_solution

