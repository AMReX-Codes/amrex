     subroutine estdt(u,u_l1,u_l2,u_l3,u_h1,u_h2,u_h3,lo,hi,dx,dt)

     use meth_params_module, only : NVAR, UX, UY, UZ

     implicit none

     integer          :: u_l1,u_l2,u_l3,u_h1,u_h2,u_h3
     integer          :: lo(3), hi(3)
     double precision :: u(u_l1:u_h1,u_l2:u_h2,u_l3:u_h3,NVAR)
     double precision :: dx(3),dt

     double precision :: umax,vmax,wmax
     integer          :: i,j,k

     umax = 0.d0
     vmax = 0.d0
     wmax = 0.d0

     do k = lo(3),hi(3)
     do j = lo(2),hi(2)
     do i = lo(1),hi(1)

         umax = max(umax,abs(u(i,j,k,UX)))
         vmax = max(vmax,abs(u(i,j,k,UY)))
         wmax = max(wmax,abs(u(i,j,k,UZ)))

     enddo
     enddo
     enddo
     
     umax = max(umax/dx(1),max(vmax/dx(2),wmax/dx(3)))

     if (umax .gt. 0.d0) then
         dt = 1.d0 / umax
     else
        print *,'ZERO VELOCITY -- TIME STEP SET TO BE DX '
        dt = dx(1)
     end if

     end subroutine estdt
