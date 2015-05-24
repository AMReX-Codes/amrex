! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine trace(q,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                       dq,qxm,qxp,qym,qyp,qzm,qzp,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                       ilo1,ilo2,ilo3,ihi1,ihi2,ihi3,dx,dy,dz,dt)

      use network           , only : nspec
      use meth_params_module, only : iorder, QVAR, QRHO, QFA, QFS, nadv
      use slope_module

      implicit none

      integer ilo1,ilo2,ilo3,ihi1,ihi2,ihi3
      integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
      integer qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3

      double precision dx, dy, dz, dt
      double precision   q( qd_l1: qd_h1, qd_l2: qd_h2, qd_l3: qd_h3,QVAR)
      double precision  dq(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3,qpd_h3,QVAR)
      double precision qxm(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3,qpd_h3,QVAR)
      double precision qxp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3,qpd_h3,QVAR)
      double precision qym(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3,qpd_h3,QVAR)
      double precision qyp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3,qpd_h3,QVAR)
      double precision qzm(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3,qpd_h3,QVAR)
      double precision qzp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3,qpd_h3,QVAR)

      ! Local variables
      integer n, iadv, ispec

      ! Default dq=0.
      dq = 0.d0

      ! Compute slopes in x-direction
      if (iorder .ne. 1) then

         call uslope( q,  qd_l1, qd_l2, qd_l3, qd_h1, qd_h2, qd_h3, &
                     dq, qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                     ilo1,ilo2,ilo3,ihi1,ihi2,ihi3,QVAR,1)

      endif

      n = QRHO
      call trace_x(n,q,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                   dq,qxm,qxp,ilo1-1,ilo2-1,ilo3-1,ihi1+2,ihi2+2,ihi3+2, &
                   ilo1,ilo2,ilo3,ihi1,ihi2,ihi3,dx,dt)

      do ispec = 1, nspec
         n = QFS + ispec - 1
         call trace_x(n,q,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                      dq,qxm,qxp,ilo1-1,ilo2-1,ilo3-1,ihi1+2,ihi2+2,ihi3+2, &
                      ilo1,ilo2,ilo3,ihi1,ihi2,ihi3,dx,dt)
      end do

      do iadv = 1, nadv
         n = QFA + iadv - 1
         call trace_x(n,q,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                      dq,qxm,qxp,ilo1-1,ilo2-1,ilo3-1,ihi1+2,ihi2+2,ihi3+2, &
                      ilo1,ilo2,ilo3,ihi1,ihi2,ihi3,dx,dt)
      end do

!     ------------------------------------------------------------------

      ! Compute slopes in y-direction
      if (iorder .ne. 1) then

         call uslope( q,  qd_l1, qd_l2, qd_l3, qd_h1, qd_h2, qd_h3, &
                     dq, qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                     ilo1,ilo2,ilo3,ihi1,ihi2,ihi3,QVAR,2)

      endif

      n = QRHO
      call trace_y(n,q,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                   dq,qym,qyp,ilo1-1,ilo2-1,ilo3-1,ihi1+2,ihi2+2,ihi3+2, &
                   ilo1,ilo2,ilo3,ihi1,ihi2,ihi3,dy,dt)

      do ispec = 1, nspec
         n = QFS + ispec - 1
         call trace_y(n,q,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                      dq,qym,qyp,ilo1-1,ilo2-1,ilo3-1,ihi1+2,ihi2+2,ihi3+2, &
                      ilo1,ilo2,ilo3,ihi1,ihi2,ihi3,dy,dt)
      end do

      do iadv = 1, nadv
         n = QFA + iadv - 1
         call trace_y(n,q,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                      dq,qym,qyp,ilo1-1,ilo2-1,ilo3-1,ihi1+2,ihi2+2,ihi3+2, &
                      ilo1,ilo2,ilo3,ihi1,ihi2,ihi3,dy,dt)
      end do

!     ------------------------------------------------------------------

      ! Compute slopes in z-direction
      if (iorder .ne. 1) then

         call uslope( q, qd_l1, qd_l2, qd_l3, qd_h1, qd_h2, qd_h3, &
                     dq, qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                     ilo1,ilo2,ilo3,ihi1,ihi2,ihi3,QVAR,3)

      endif

      n = QRHO
      call trace_z(n,q,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                   dq,qzm,qzp,ilo1-1,ilo2-1,ilo3-1,ihi1+2,ihi2+2,ihi3+2, &
                   ilo1,ilo2,ilo3,ihi1,ihi2,ihi3,dz,dt)

      do ispec = 1, nspec
         n = QFS + ispec - 1
         call trace_z(n,q,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                      dq,qzm,qzp,ilo1-1,ilo2-1,ilo3-1,ihi1+2,ihi2+2,ihi3+2, &
                      ilo1,ilo2,ilo3,ihi1,ihi2,ihi3,dz,dt)
      end do

      do iadv = 1, nadv
         n = QFA + iadv - 1
         call trace_z(n,q,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                      dq,qzm,qzp,ilo1-1,ilo2-1,ilo3-1,ihi1+2,ihi2+2,ihi3+2, &
                      ilo1,ilo2,ilo3,ihi1,ihi2,ihi3,dz,dt)
      end do

      end subroutine trace

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine trace_x(n,q,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                         dq,qxm,qxp,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                         ilo1,ilo2,ilo3,ihi1,ihi2,ihi3,dx,dt)

      use network           , only : nspec
      use meth_params_module, only : QVAR, QU

      implicit none

      integer n
      integer ilo1,ilo2,ilo3,ihi1,ihi2,ihi3
      integer  qd_l1, qd_l2, qd_l3, qd_h1, qd_h2, qd_h3
      integer qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3

      double precision :: dx, dt
      double precision ::   q( qd_l1: qd_h1, qd_l2: qd_h2, qd_l3: qd_h3,QVAR)
      double precision ::  dq(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
      double precision :: qxm(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
      double precision :: qxp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)

      ! Local variables
      integer          :: i, j, k
      double precision :: dtdx
      double precision :: u, spzero, acmprght, acmpleft

      dtdx = dt/dx

      do k = ilo3-1, ihi3+1
      do j = ilo2-1, ihi2+1

         ! Right state
         do i = ilo1, ihi1+1
            u = q(i,j,k,QU)
            if (u .gt. 0.d0) then
               spzero = -1.d0
            else
               spzero = u*dtdx
            endif
            acmprght = 0.5d0*(-1.d0 - spzero )*dq(i,j,k,n)
            qxp(i,j,k,n) = q(i,j,k,n) + acmprght
         enddo

         ! Left state
         do i = ilo1-1, ihi1
            u = q(i,j,k,QU)
            if (u .ge. 0.d0) then
               spzero = u*dtdx
            else
               spzero = 1.d0
            endif
            acmpleft = 0.5d0*(1.d0 - spzero )*dq(i,j,k,n)
            qxm(i+1,j,k,n) = q(i,j,k,n) + acmpleft
         enddo

      enddo
      enddo

      end subroutine trace_x

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine trace_y(n,q,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                         dq,qym,qyp,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                         ilo1,ilo2,ilo3,ihi1,ihi2,ihi3,dy,dt)

      use network, only : nspec
      use meth_params_module, only : QVAR, QV

      implicit none

      integer n
      integer ilo1,ilo2,ilo3,ihi1,ihi2,ihi3
      integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
      integer qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3

      double precision ::  dy, dt
      double precision ::   q(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision ::  dq(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
      double precision :: qym(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
      double precision :: qyp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)

      ! Local variables
      integer          :: i, j, k
      double precision :: v, spzero, acmptop, acmpbot
      double precision :: dtdy

      dtdy = dt/dy

      do k = ilo3-1, ihi3+1
      do i = ilo1-1, ihi1+1

         ! Top state
         do j = ilo2, ihi2+1
            v = q(i,j,k,QV)
            if (v .gt. 0.d0) then
               spzero = -1.d0
            else
               spzero = v*dtdy
            endif

            acmptop = 0.5d0*(-1.d0 - spzero )*dq(i,j,k,n)

            qyp(i,j,k,n) = q(i,j,k,n) + acmptop

         enddo

         ! Bottom state
         do j = ilo2-1, ihi2
            v = q(i,j,k,QV)
            if (v .ge. 0.d0) then
               spzero = v*dtdy
            else
               spzero = 1.d0
            endif
            acmpbot = 0.5d0*(1.d0 - spzero )*dq(i,j,k,n)
            qym(i,j+1,k,n) = q(i,j,k,n) + acmpbot
         enddo

      enddo
      enddo

      end subroutine trace_y

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine trace_z(n,q,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                         dq,qzm,qzp,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                         ilo1,ilo2,ilo3,ihi1,ihi2,ihi3,dz,dt)

      use network, only : nspec
      use meth_params_module, only : QVAR, QW

      implicit none

      integer n
      integer ilo1,ilo2,ilo3,ihi1,ihi2,ihi3
      integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
      integer qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3

      double precision ::  dz, dt
      double precision ::   q(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision ::  dq(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
      double precision :: qzm(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
      double precision :: qzp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)

      ! Local variables
      integer          :: i, j, k
      double precision :: w, spzero, acmptop, acmpbot
      double precision :: dtdz

      dtdz = dt/dz

      do j = ilo2-1, ihi2+1
      do i = ilo1-1, ihi1+1

         ! Top state
         do k = ilo3, ihi3+1
            w = q(i,j,k,QW)
            if (w .gt. 0.d0) then
               spzero = -1.d0
            else
               spzero = w*dtdz
            endif
            acmptop = 0.5d0*(-1.d0 - spzero )*dq(i,j,k,n)
            qzp(i,j,k,n) = q(i,j,k,n) + acmptop
         enddo

         ! Bottom state
         do k = ilo3-1, ihi3
            w = q(i,j,k,QW)
            if (w .ge. 0.d0) then
               spzero = w*dtdz
            else
               spzero = 1.d0
            endif
            acmpbot = 0.5d0*(1.d0 - spzero )*dq(i,j,k,n)
            qzm(i,j,k+1,n) = q(i,j,k,n) + acmpbot
         enddo

      enddo
      enddo

      end subroutine trace_z
