! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine trace(q,qd_l1,qd_l2,qd_h1,qd_h2, &
                       dq,qxm,qxp,qym,qyp,qpd_l1,qpd_l2,qpd_h1,qpd_h2, &
                       ilo1,ilo2,ihi1,ihi2,dx,dy,dt)

      use meth_params_module, only : QVAR, QFA, nadv

      implicit none

      integer ilo1,ilo2,ihi1,ihi2
      integer qd_l1,qd_l2,qd_h1,qd_h2
      integer qpd_l1,qpd_l2,qpd_h1,qpd_h2

      double precision dx, dy, dt
      double precision     q(qd_l1:qd_h1,qd_l2:qd_h2,QVAR)
      double precision  dq(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
      double precision qxm(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
      double precision qxp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
      double precision qym(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
      double precision qyp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)

      ! Local variables
      integer i, j
      integer n, iadv

      ! Default dq=0.
      dq(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:QVAR) = 0.d0

      do iadv = 1, nadv
         n = QFA + iadv - 1
         call trace_x(n,q,qd_l1,qd_l2,qd_h1,qd_h2, &
                      dq,qxm,qxp,ilo1-1,ilo2-1,ihi1+2,ihi2+2, &
                      ilo1,ilo2,ihi1,ihi2,dx,dt)
      end do

!     ------------------------------------------------------------------

      do iadv = 1, nadv
         n = QFA + iadv - 1
         call trace_y(n,q,qd_l1,qd_l2,qd_h1,qd_h2, &
                      dq,qym,qyp,ilo1-1,ilo2-1,ihi1+2,ihi2+2, &
                      ilo1,ilo2,ihi1,ihi2,dy,dt)
      end do

      end subroutine trace

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine trace_x(n,q,qd_l1,qd_l2,qd_h1,qd_h2, &
                         dq,qxm,qxp,qpd_l1,qpd_l2,qpd_h1,qpd_h2, &
                         ilo1,ilo2,ihi1,ihi2,dx,dt)

      use meth_params_module, only : QVAR, QU, QFA, nadv

      implicit none

      integer n
      integer ilo1,ilo2,ihi1,ihi2
      integer qd_l1,qd_l2,qd_h1,qd_h2
      integer qpd_l1,qpd_l2,qpd_h1,qpd_h2

      double precision :: dx, dt
      double precision ::   q(qd_l1:qd_h1,qd_l2:qd_h2,QVAR)
      double precision ::  dq(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
      double precision :: qxm(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
      double precision :: qxp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)

      ! Local variables
      integer          :: i, j
      double precision :: dtdx
      double precision :: u, spzero, acmprght, acmpleft

      dtdx = dt/dx

      ! Compute slopes in x-direction
      call slope( q, qd_l1, qd_l2, qd_h1, qd_h2, &
                 dq, qpd_l1,qpd_l2,qpd_h1,qpd_h2, &
                 ilo1,ilo2,ihi1,ihi2,QVAR,1)

      do j = ilo2-1, ihi2+1

         ! Right state
         do i = ilo1, ihi1+1
            u = q(i,j,QU)
            if (u .gt. 0.d0) then
               spzero = -1.d0
            else
               spzero = u*dtdx
            endif
            acmprght = 0.5d0*(-1.d0 - spzero )*dq(i,j,n)
            qxp(i,j,n) = q(i,j,n) + acmprght

            if (abs(ugdx(i,j)) < eps) then
                phi_edge = 0.5d0 * (phi(i,j) + phi(i-1,j))
            else if (ugdx(i,j) > 0.d0) then
                phi_edge = phi(i-1,j) + 0.5d0 * (1.d0 - uadv*dt/dx)*dq(i-1,j)
            else
                phi_edge = phi(i  ,j) - 0.5d0 * (1.d0 + uadv*dt/dx)*dq(i  ,j)
            endif
  
         enddo

         ! Left state
         do i = ilo1-1, ihi1
            u = q(i,j,QU)
            if (u .ge. 0.d0) then
               spzero = u*dtdx
            else
               spzero = 1.d0
            endif
            acmpleft = 0.5d0*(1.d0 - spzero )*dq(i,j,n)
            qxm(i+1,j,n) = q(i,j,n) + acmpleft
         enddo

      enddo

      end subroutine trace_x

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine trace_y(n,q,qd_l1,qd_l2,qd_h1,qd_h2, &
                         dq,qym,qyp,qpd_l1,qpd_l2,qpd_h1,qpd_h2, &
                         ilo1,ilo2,ihi1,ihi2,dy,dt)

      use meth_params_module, only : QVAR, QV, QFA, nadv

      implicit none

      integer n
      integer ilo1,ilo2,ihi1,ihi2
      integer qd_l1,qd_l2,qd_h1,qd_h2
      integer qpd_l1,qpd_l2,qpd_h1,qpd_h2

      double precision ::  dy, dt
      double precision ::   q(qd_l1:qd_h1,qd_l2:qd_h2,QVAR)
      double precision ::  dq(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
      double precision :: qym(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
      double precision :: qyp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)

      ! Local variables
      integer          :: i, j
      double precision :: v, spzero, acmptop, acmpbot
      double precision :: dtdy

      dtdy = dt/dy

      ! Compute slopes in y-direction
      call slope(q,qd_l1,qd_l2,qd_h1,qd_h2, &
                 dq,qpd_l1,qpd_l2,qpd_h1,qpd_h2, &
                 ilo1,ilo2,ihi1,ihi2,QVAR,2)

      do i = ilo1-1, ihi1+1

         ! Top state
         do j = ilo2, ihi2+1
            v = q(i,j,QV)
            if (v .gt. 0.d0) then
               spzero = -1.d0
            else
               spzero = v*dtdy
            endif
            acmptop = 0.5d0*(-1.d0 - spzero )*dq(i,j,n)
            qyp(i,j,n) = q(i,j,n) + acmptop
         enddo

         ! Bottom state
         do j = ilo2-1, ihi2
            v = q(i,j,QV)
            if (v .ge. 0.d0) then
               spzero = v*dtdy
            else
               spzero = 1.d0
            endif
            acmpbot = 0.5d0*(1.d0 - spzero )*dq(i,j,n)
            qym(i,j+1,n) = q(i,j,n) + acmpbot
         enddo

      enddo

      end subroutine trace_y
