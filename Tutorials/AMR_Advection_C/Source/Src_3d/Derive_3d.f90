!-----------------------------------------------------------------------

      subroutine dermagvel(magvel,vel_l1,vel_l2,vel_l3,vel_h1,vel_h2,vel_h3,nv, &
                           dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc,lo,hi,domlo, &
                           domhi,delta,xlo,time,dt,bc,level,grid_no)
      !
      ! This routine will derive magnitude of velocity.
      !
      implicit none

      integer          lo(3), hi(3)
      integer          vel_l1,vel_l2,vel_l3,vel_h1,vel_h2,vel_h3,nv
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3)
      integer          bc(3,2,nc)
      double precision delta(3), xlo(3), time, dt
      double precision magvel(vel_l1:vel_h1,vel_l2:vel_h2,vel_l3:vel_h3,nv)
      double precision    dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
      integer    level, grid_no

      integer i,j,k

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               magvel(i,j,k,1) = sqrt( (dat(i,j,k,2) / dat(i,j,k,1))**2 + &
                                       (dat(i,j,k,3) / dat(i,j,k,1))**2 + &
                                       (dat(i,j,k,4) / dat(i,j,k,1))**2 )
            end do
         end do
      end do

      end subroutine dermagvel

!-----------------------------------------------------------------------

      subroutine derlapvar(var,var_l1,var_l2,var_h1,var_h2,nv, &
                           dat,dat_l1,dat_l2,dat_h1,dat_h2,nc,lo,hi,domlo, &
                           domhi,delta,xlo,time,dt,bc,level,grid_no)
!
!     This routine will derive the weighted-Laplacian of the variable for
!       the purposes of error estimation
!
      implicit none

      integer          lo(2), hi(2)
      integer          var_l1,var_l2,var_h1,var_h2,nv
      integer          dat_l1,dat_l2,dat_h1,dat_h2,nc
      integer          domlo(2), domhi(2)
      integer          bc(2,2,nc)
      double precision delta(2), xlo(2), time, dt
      double precision var(var_l1:var_h1,var_l2:var_h2,nv)
      double precision dat(dat_l1:dat_h1,dat_l2:dat_h2,nc)
      integer    level, grid_no

      double precision ::  delu(2,var_l1:var_h1,var_l2:var_h2)
      double precision :: delua(2,var_l1:var_h1,var_l2:var_h2)
      double precision :: delu2(4), delu3(4), delu4(4)
      double precision :: num, denom
      integer          :: i,j

      ! This value is taken from FLASH
      double precision, parameter:: epsil=0.02

      ! adapted from ref_marking.f90 in FLASH2.5

      ! d/dx
      do j=lo(2)-1,hi(2)+1
      do i=lo(1)-1,hi(1)+1
          delu(1,i,j) =     dat(i+1,j,1) -      dat(i-1,j,1)
         delua(1,i,j) = abs(dat(i+1,j,1)) + abs(dat(i-1,j,1))
      end do
      end do

      ! d/dy
      do j=lo(2)-1,hi(2)+1
      do i=lo(1)-1,hi(1)+1
          delu(2,i,j) =     dat(i,j+1,1) -      dat(i,j-1,1)
         delua(2,i,j) = abs(dat(i,j+1,1)) + abs(dat(i,j-1,1))
      end do
      end do

      do j = lo(2),hi(2)
      do i = lo(1),hi(1)

         ! d/dxdx
         delu2(1) =     delu(1,i+1,j)  -     delu(1,i-1,j)
         delu3(1) = abs(delu(1,i+1,j)) + abs(delu(1,i-1,j))
         delu4(1) =    delua(1,i+1,j)  +    delua(1,i-1,j)

         ! d/dydx
         delu2(2) =     delu(1,i,j+1)  -     delu(1,i,j-1)
         delu3(2) = abs(delu(1,i,j+1)) + abs(delu(1,i,j-1))
         delu4(2) =    delua(1,i,j+1)  +    delua(1,i,j-1)

         ! d/dxdy
         delu2(3) =     delu(2,i+1,j)  -     delu(2,i-1,j)
         delu3(3) = abs(delu(2,i+1,j)) + abs(delu(2,i-1,j))
         delu4(3) =    delua(2,i+1,j)  +    delua(2,i-1,j)

         ! d/dydy
         delu2(4) =     delu(2,i,j+1)  -     delu(2,i,j-1)
         delu3(4) = abs(delu(2,i,j+1)) + abs(delu(2,i,j-1))
         delu4(4) =    delua(2,i,j+1)  +    delua(2,i,j-1)

         ! compute the error
         num   =  delu2(1)**2 + delu2(2)**2 + delu2(3)**2 + delu2(4)**2
         denom = (delu3(1) + (epsil*delu4(1)+1.d-99))**2 + &
                 (delu3(2) + (epsil*delu4(2)+1.d-99))**2 + &
                 (delu3(3) + (epsil*delu4(3)+1.d-99))**2 + &
                 (delu3(4) + (epsil*delu4(4)+1.d-99))**2

         var(i,j,1) = sqrt(num/denom)

      end do
      end do


      end subroutine derlapvar

!-----------------------------------------------------------------------

      subroutine dermagvort(vort,v_l1,v_l2,v_l3,v_h1,v_h2,v_h3,nv, & 
                            dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc,lo,hi,domlo, &
                            domhi,delta,xlo,time,dt,bc,level,grid_no)
      !
      ! This routine will calculate vorticity
      !     

      use bl_constants_module

      implicit none

      integer          lo(3), hi(3)
      integer            v_l1,  v_l2,  v_l3,  v_h1,  v_h2,  v_h3,nv
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3), level, grid_no
      integer          bc(3,2,nc)
      double precision delta(3), xlo(3), time, dt
      double precision vort(  v_l1:  v_h1,  v_l2:  v_h2,  v_l3:  v_h3,nv)
      double precision, intent(in) :: dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)

      integer          :: i,j,k
      double precision :: uy,uz,vx,vz,wx,wy,v1,v2,v3
      double precision :: ldat(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,2:4)
      !
      ! Convert momentum to velocity.
      !
      do k = lo(3)-1, hi(3)+1
         do j = lo(2)-1, hi(2)+1
            do i = lo(1)-1, hi(1)+1
               ldat(i,j,k,2) = dat(i,j,k,2) / dat(i,j,k,1)
               ldat(i,j,k,3) = dat(i,j,k,3) / dat(i,j,k,1)
               ldat(i,j,k,4) = dat(i,j,k,4) / dat(i,j,k,1)
            end do
         end do
      end do
      !
      ! Calculate vorticity.
      !
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               uy = HALF * (ldat(i,j+1,k,2) - ldat(i,j-1,k,2)) / delta(2)
               uz = HALF * (ldat(i,j,k+1,2) - ldat(i,j,k-1,2)) / delta(3)
               vx = HALF * (ldat(i+1,j,k,3) - ldat(i-1,j,k,3)) / delta(1)
               vz = HALF * (ldat(i,j,k+1,3) - ldat(i,j,k-1,3)) / delta(3)
               wx = HALF * (ldat(i+1,j,k,4) - ldat(i-1,j,k,4)) / delta(1)
               wy = HALF * (ldat(i,j+1,k,4) - ldat(i,j-1,k,4)) / delta(2)
               v1 = wy - vz
               v2 = uz - wx
               v3 = vx - uy
               vort(i,j,k,1) = sqrt(v1*v1 + v2*v2 + v3*v3)
            end do
         end do
      end do

      end subroutine dermagvort

!-----------------------------------------------------------------------

      subroutine derspec(spec,spec_l1,spec_l2,spec_l3,spec_h1,spec_h2,spec_h3,nv, &
                         dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc,lo,hi,domlo, &
                         domhi,delta,xlo,time,dt,bc,level,grid_no)
      !
      ! This routine will derive the velocity from the momentum.
      !
      implicit none

      integer          lo(3), hi(3)
      integer          spec_l1,spec_l2,spec_l3,spec_h1,spec_h2,spec_h3,nv
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3)
      integer          bc(3,2,nc)
      double precision delta(3), xlo(3), time, dt
      double precision spec(spec_l1:spec_h1,spec_l2:spec_h2,spec_l3:spec_h3,nv)
      double precision dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
      integer    level, grid_no
 
      integer i,j,k
 
      print *,'IN DERSPEC '
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               spec(i,j,k,1) = dat(i,j,k,2) / dat(i,j,k,1)
            end do
         end do
      end do
 
      end subroutine derspec


!-----------------------------------------------------------------------

      subroutine derstate(state,state_l1,state_l2,state_l3,state_h1,state_h2,state_h3,nv,&
                          dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc,lo,hi,domlo,&
                          domhi,delta,xlo,time,dt,bc,level,grid_no)
!
!     This routine will derive the X from the (rho X)
!
      use network, only : nspec
      implicit none 

      integer          lo(3), hi(3)
      integer          state_l1,state_l2,state_l3,state_h1,state_h2,state_h3,nv
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3)
      integer          bc(3,2,nc)
      double precision delta(3), xlo(3), time, dt
      double precision state(state_l1:state_h1,state_l2:state_h2,state_l3:state_h3,nv)
      double precision   dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
      integer    level, grid_no
 
      integer    i,j,k
 
      ! Density 
      do k = lo(3), hi(3)
      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            state(i,j,k,1:nspec) = dat(i,j,k,2:nspec+1) / dat(i,j,k,1)
         end do
      end do
      end do
 
      end subroutine derstate
