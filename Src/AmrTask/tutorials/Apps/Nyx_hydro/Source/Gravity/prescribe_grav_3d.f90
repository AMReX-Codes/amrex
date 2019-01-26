module ps_grav

  use amrex_fort_module, only : rt => amrex_real

  implicit none

  real(rt), private, save :: dm_const,dm_const2
  logical, private, save :: isnt_init = .true.
  public :: ps_grav_init,ps_grav_accel,ps_grav_phi
  contains

      subroutine ps_grav_init()
        use fundamental_constants_module, only : Gconst,Mpc,M_s
        use probdata_module, only : dcenx,dceny,dcenz,dmconc,dmmass,dmscale
        implicit none
        if(isnt_init) then
          dm_const = -(dmmass/M_s)*Gconst / (log(1.0d0+dmconc) - dmconc/(1.0d0+dmconc))
          dm_const2 = dmconc/(dmscale/Mpc)
          isnt_init = .false.
        else
          write(*,*)'ps_grav is already initialized, why are we here?'
        endif
      end subroutine ps_grav_init

!     This is an example of how to specify a radial profile
!     Note that r_c and rho_c must be specified in probdata_module
      pure function ps_grav_accel(r) result(g)
        implicit none
        real(rt), intent(in) :: r!,smbh_const
        real(rt)             :: g
!         ! put function for g(r) in here
        g = dm_const*(dlog(1.0d0+dm_const2*r)/(r*r)-dm_const2 &
               /(r*(1.0d0+r*dm_const2)))
      end function ps_grav_accel

      pure function ps_grav_phi(r) result(g)
        implicit none
        real(rt), intent(in) :: r
        real(rt)             :: g
!         ! put function for g(r) in here
        g = -dm_const*dlog(1.0d0+r*dm_const2)/r
      end function ps_grav_phi
end module ps_grav
 
subroutine fort_prescribe_grav (lo,hi,dx, &
     grav,g_l1,g_l2,g_l3,g_h1,g_h2,g_h3,&
     problo,add)

  use probdata_module, only : dcenx,dceny,dcenz
  use bl_constants_module, only : HALF
  use ps_grav, only : ps_grav_accel
  implicit none
  integer          :: g_l1,g_l2,g_l3,g_h1,g_h2,g_h3,add
  integer          :: lo(3),hi(3)
  real(rt) :: grav(g_l1:g_h1,g_l2:g_h2,g_l3:g_h3,3)
  real(rt) :: dx(3)
  real(rt) :: problo(3)

  integer          :: i,j,k
  real(rt) :: fort_prescribe_grav_gravityprofile
  real(rt) :: x,y,z,y2,z2
  real(rt) :: r,maggrav
  real(rt) :: dxm,r1
  dxm = min(dx(1),dx(2),dx(3))
!
! This is an example of how to use the radial profile above.
!
  if (add.eq.0) then
     do k = lo(3), hi(3)
        z = problo(3) + (dble(k)+HALF) * dx(3) - dcenz
        z2 = z*z
        do j = lo(2), hi(2)
           y = problo(2) + (dble(j)+HALF) * dx(2) - dceny
           y2 = y*y
           do i = lo(1), hi(1)
              x = problo(1) + (dble(i)+HALF) * dx(1) - dcenx
              r1 = dsqrt(x*x+y2+z2)
              r = max(dxm,r1)

              maggrav = ps_grav_accel(r)

              !              Put in angular dependence
              grav(i,j,k,1) = maggrav*x/r
              grav(i,j,k,2) = maggrav*y/r
              grav(i,j,k,3) = maggrav*z/r
           enddo
        enddo
     enddo
  else
     do k = lo(3), hi(3)
        z = problo(3) + (dble(k)+HALF) * dx(3) - dcenz
        z2 = z*z
        do j = lo(2), hi(2)
           y = problo(2) + (dble(j)+HALF) * dx(2) - dceny
           y2 = y*y
           do i = lo(1), hi(1)
              x = problo(1) + (dble(i)+HALF) * dx(1) - dcenx
              r1 = dsqrt(x*x+y2+z2)
              r = max(dxm,r1)

              maggrav = ps_grav_accel(r)

              !              Put in angular dependence
              grav(i,j,k,1) = grav(i,j,k,1) + maggrav*x/r
              grav(i,j,k,2) = grav(i,j,k,2) + maggrav*y/r
              grav(i,j,k,3) = grav(i,j,k,3) + maggrav*z/r
           enddo
        enddo
     enddo
  endif

end subroutine fort_prescribe_grav

