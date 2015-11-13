     subroutine estdt(u,u_l1,u_l2,u_h1,u_h2,lo,hi,dx,dt)

     use meth_params_module, only : NVAR, UX, UY

     implicit none

     integer          :: u_l1,u_l2,u_h1,u_h2
     integer          :: lo(2), hi(2)
     double precision :: u(u_l1:u_h1,u_l2:u_h2,NVAR)
     double precision :: dx(2),dt

     double precision :: xvel,yvel,dt1,dt2
     integer          :: i,j

      do j = lo(2),hi(2)
         do i = lo(1),hi(1)

            xvel = u(i,j,UX)
            yvel = u(i,j,UY)

            dt1 = dx(1)/abs(xvel)
            dt2 = dx(2)/abs(yvel)

            dt = min(dt,dt1,dt2)

         enddo
      enddo

      end subroutine estdt
