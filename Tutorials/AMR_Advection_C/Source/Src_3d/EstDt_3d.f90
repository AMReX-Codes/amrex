     subroutine estdt(u,u_l1,u_l2,u_l3,u_h1,u_h2,u_h3,lo,hi,dx,dt)

     use meth_params_module, only : NVAR, UX, UY, UZ

     implicit none

     integer          :: u_l1,u_l2,u_l3,u_h1,u_h2,u_h3
     integer          :: lo(3), hi(3)
     double precision :: u(u_l1:u_h1,u_l2:u_h2,u_l3:u_h3,NVAR)
     double precision :: dx(3),dt

     double precision :: xvel,yvel,zvel,dt1,dt2,dt3
     integer          :: i,j,k

      do k = lo(3),hi(3)
      do j = lo(2),hi(2)
      do i = lo(1),hi(1)

         xvel = u(i,j,k,UX)
         yvel = u(i,j,k,UY)
         zvel = u(i,j,k,UZ)

         dt1 = dx(1)/abs(xvel)
         dt2 = dx(2)/abs(yvel)
         dt2 = dx(3)/abs(zvel)

         dt = min(min(min(dt,dt1),dt2),dt3)

      enddo
      enddo
      enddo

      end subroutine estdt
