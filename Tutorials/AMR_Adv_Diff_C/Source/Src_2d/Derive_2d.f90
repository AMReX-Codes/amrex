!-----------------------------------------------------------------------

      subroutine dermagvel(magvel,vel_l1,vel_l2,vel_h1,vel_h2,nv, &
                           dat,dat_l1,dat_l2,dat_h1,dat_h2,nc,lo,hi,domlo, &
                           domhi,dx,xlo,time,dt,bc,level,grid_no)
!
!     This routine will derive magnitude of velocity.
!
      implicit none 

      integer          lo(2), hi(2)
      integer          vel_l1,vel_l2,vel_h1,vel_h2,nv
      integer          dat_l1,dat_l2,dat_h1,dat_h2,nc
      integer          domlo(2), domhi(2)
      integer          bc(2,2,nc)
      double precision dx(2), xlo(2), time, dt
      double precision magvel(vel_l1:vel_h1,vel_l2:vel_h2,nv)
      double precision    dat(dat_l1:dat_h1,dat_l2:dat_h2,nc)
      integer    level, grid_no

      integer    i,j

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            magvel(i,j,1) = sqrt( dat(i,j,2)**2 + dat(i,j,3)**2 )
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

      subroutine dermagvort(vort,v_l1,v_l2,v_h1,v_h2,nv, &
                            dat,dat_l1,dat_l2,dat_h1,dat_h2,nc,lo,hi,domlo, &
                            domhi,delta,xlo,time,dt,bc,level,grid_no)
!
!     This routine will calculate vorticity
!
      implicit none

      integer          lo(2), hi(2)
      integer            v_l1,  v_l2,  v_h1,  v_h2,nv
      integer          dat_l1,dat_l2,dat_h1,dat_h2,nc
      integer          domlo(2), domhi(2)
      integer          bc(2,2,nc)
      double precision delta(2), xlo(2), time, dt
      double precision vort(  v_l1:  v_h1,  v_l2:  v_h2,nv)
      double precision  dat(dat_l1:dat_h1,dat_l2:dat_h2,nc)
      integer    level, grid_no

      integer          :: i,j
      double precision :: vx,uy

      ! Calculate vorticity
      do j = lo(2), hi(2)
      do i = lo(1), hi(1)
         vx = 0.5d0 * (dat(i+1,j,2) - dat(i-1,j,2)) / delta(1) 
         uy = 0.5d0 * (dat(i,j+1,1) - dat(i,j-1,1)) / delta(2) 
         vort(i,j,1) = abs(vx - uy)
      end do
      end do

      end subroutine dermagvort

