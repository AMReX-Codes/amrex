! :::
! ::: ------------------------------------------------------------------
! :::
    subroutine tracexy_cg(q,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                         dqx,dqy,dq_l1,dq_l2,dq_l3,dq_h1,dq_h2,dq_h3, &
                         qxm,qxp,qym,qyp,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                         ilo1,ilo2,ihi1,ihi2,dx,dy,dt,kc,k3d,a_old)

      use amrex_fort_module, only : rt => amrex_real
      use meth_params_module, only : QVAR, QRHO, QU, QV, QW, &
                                     QREINT, QPRES, &
                                     npassive, qpass_map, &
                                     small_dens, small_pres, gamma_minus_1
      implicit none

      integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
      integer dq_l1,dq_l2,dq_l3,dq_h1,dq_h2,dq_h3
      integer qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3
      integer ilo1,ilo2,ihi1,ihi2
      integer kc,k3d

      real(rt)    q(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      real(rt)    c(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
      real(rt)  dqx(dq_l1:dq_h1,dq_l2:dq_h2,dq_l3:dq_h3,QVAR)
      real(rt)  dqy(dq_l1:dq_h1,dq_l2:dq_h2,dq_l3:dq_h3,QVAR)

      real(rt) qxm(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
      real(rt) qxp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
      real(rt) qym(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
      real(rt) qyp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
      real(rt) a_old
      real(rt) dx, dy, dt
      real(rt), parameter :: eps = 1.d-3

      ! Local variables
      integer i, j, n

      real(rt) :: dtdx, dtdy
      real(rt) :: cc, csq, rho, u, v, w, p
      real(rt) :: drho, du, dv, dw, dp
      real(rt) :: spzero

      real(rt) :: t1,t2,t3
      real(rt) :: beta0l,beta0r,betaml,betapr
      real(rt) :: clag,dd,dd0,dleft,dleft0,drght,drght0
      real(rt) :: ptr0l,ptr0r,rtilde,rtr0l,rtr0r
      integer          :: ipassive

      dtdx = dt/(dx*a_old)
      dtdy = dt/(dy*a_old)

      !!!!!!!!!!!!!!!
      ! Colella-Glaz algorithm
      !!!!!!!!!!!!!!!

      ! Compute left and right traced states
      do j = ilo2-1, ihi2+1
         do i = ilo1-1, ihi1+1

            cc   = c(i,j,k3d)
            csq  = cc*cc
            rho  = q(i,j,k3d,QRHO)
            u    = q(i,j,k3d,QU)
            v    = q(i,j,k3d,QV)
            w    = q(i,j,k3d,QW)
            p    = q(i,j,k3d,QPRES)

            drho  = dqx(i,j,kc,QRHO)
            du    = dqx(i,j,kc,QU)
            dv    = dqx(i,j,kc,QV)
            dw    = dqx(i,j,kc,QW)
            dp    = dqx(i,j,kc,QPRES)

            clag = rho*cc
            t1   = 1.0d0 / clag
            t2   = 0.5d0 * dtdx/rho
            t3   = t1*t1

            ! **************************************************************************

            drght  = -0.5d0 - 0.5d0* min(0.d0,u-cc)*dtdx
            drght0 = -0.5d0 * (1.d0+u*dtdx)
            dd     = drght
            dd0    = drght0
            rtilde = rho + drght*drho
            rtr0r  = rho + drght0*drho
            ptr0r  = p + drght0*dp

            betapr = -(dp*t1+du)*t2
            beta0r =  (dp*t3 -drho/(rtilde*rtr0r)) *0.5d0*dtdx*cc

            if (u+cc .ge. 0.d0) then
                betapr = 0.d0
            endif
            if (u .ge. 0.d0) then
                beta0r = 0.d0
            endif

            if (u .gt. 0.d0) then
               spzero = -1.d0
            else
               spzero = u*dtdx
            endif

            if (i .ge. ilo1) then
               qxp(i,j,kc,QRHO)  = 1.d0/(1.d0/rtilde - beta0r - betapr)
               qxp(i,j,kc,QPRES) = p + drght*dp + clag*clag*betapr
               qxp(i,j,kc,QU)    = u + drght*du + clag*     betapr
               qxp(i,j,kc,QV)    = v + 0.5d0*(-1.d0 - spzero )*dv
               qxp(i,j,kc,QW)    = w + 0.5d0*(-1.d0 - spzero )*dw

               ! If rho or p too small, set all the slopes to zero
               if (qxp(i,j,kc,QRHO ) .lt. small_dens .or. &
                   qxp(i,j,kc,QPRES) .lt. small_pres) then
                  qxp(i,j,kc,QPRES) = p
                  qxp(i,j,kc,QRHO)  = rho
                  qxp(i,j,kc,QU)    = u
               end if

               qxp(i,j,kc,QREINT) = qxp(i,j,kc,QPRES) / gamma_minus_1

            end if

            ! **************************************************************************

            dleft  = 0.5d0 - max(0.d0,u + cc)*0.5d0*dtdx
            dleft0 = 0.5d0 * (1.d0-u*dtdx)
            dd0    = dleft0
            dd     = dleft
            rtilde = rho + dleft*drho
            rtr0l  = rho + dleft0*drho
            ptr0l  = p + dleft0*dp

            betaml =  (dp*t1-du)*t2
            beta0l = -(dp*t3 - drho/(rtilde*rtr0l))*0.5d0*dtdx*cc

            if (u-cc .lt. 0.d0) then
              betaml = 0.d0
            endif
            if (u .lt. 0.d0) then
              beta0l = 0.d0
            endif

            if (u .ge. 0.d0) then
               spzero = u*dtdx
            else
               spzero = 1.d0
            endif

            if (i .le. ihi1) then

               qxm(i+1,j,kc,QRHO)  = 1.d0/(1.d0/rtilde - beta0l  - betaml)
               qxm(i+1,j,kc,QPRES) = p + dleft*dp + clag*clag*betaml
               qxm(i+1,j,kc,QU)    = u + dleft*du - clag*     betaml
               qxm(i+1,j,kc,QV)    = v + 0.5d0*(1.d0 - spzero )*dv
               qxm(i+1,j,kc,QW)    = w + 0.5d0*(1.d0 - spzero )*dw

               ! If rho or p too small, set all the slopes to zero
               if (qxm(i+1,j,kc,QRHO ) .lt. small_dens .or. &
                   qxm(i+1,j,kc,QPRES) .lt. small_pres) then
                  qxm(i+1,j,kc,QRHO)  = rho
                  qxm(i+1,j,kc,QPRES) = p
                  qxm(i+1,j,kc,QU)    = u
               end if

               qxm(i+1,j,kc,QREINT) = qxm(i+1,j,kc,QPRES) / gamma_minus_1

            endif

         enddo
      enddo

      ! Do all of the passively advected quantities in one loop
      do ipassive = 1, npassive
         n = qpass_map(ipassive)

         do j = ilo2-1, ihi2+1
            ! Right state
            do i = ilo1, ihi1+1
               u = q(i,j,k3d,QU)
               if (u .gt. 0.d0) then
                  spzero = -1.d0
               else
                  spzero = u*dtdx
               endif
               qxp(i,j,kc,n) = q(i,j,k3d,n) + 0.5d0*(-1.d0 - spzero )*dqx(i,j,kc,n)
            enddo

            ! Left state
            do i = ilo1-1, ihi1
               u = q(i,j,k3d,QU)
               if (u .ge. 0.d0) then
                  spzero = u*dtdx
               else
                  spzero = 1.d0
               endif
               qxm(i+1,j,kc,n) = q(i,j,k3d,n) + 0.5d0*(1.d0 - spzero )*dqx(i,j,kc,n)
            enddo
         enddo
      enddo

      do j = ilo2-1, ihi2+1
         do i = ilo1-1, ihi1+1

            cc = c(i,j,k3d)
            csq = cc*cc
            rho = q(i,j,k3d,QRHO)
            u = q(i,j,k3d,QU)
            v = q(i,j,k3d,QV)
            w = q(i,j,k3d,QW)
            p = q(i,j,k3d,QPRES)

            drho  = dqy(i,j,kc,QRHO)
            du    = dqy(i,j,kc,QU)
            dv    = dqy(i,j,kc,QV)
            dw    = dqy(i,j,kc,QW)
            dp    = dqy(i,j,kc,QPRES)

            clag = rho*cc
            t1   = 1.0d0 / clag
            t2   = 0.5d0 * dtdy/rho
            t3   = t1*t1

            ! **************************************************************************

            drght  = -0.5d0 - 0.5d0 * min(0.d0,v-cc)*dtdy
            drght0 = -0.5d0 * (1.d0+v*dtdy)
            dd     = drght
            dd0    = drght0
            rtilde = rho + drght*drho
            rtr0r  = rho + drght0*drho
            ptr0r  = p + drght0*dp

            betapr = -(dp*t1+dv)*t2
            beta0r =  (dp*t3 -drho/(rtilde*rtr0r)) *0.5d0*dtdy*cc

            if (v+cc .ge. 0.d0) then
                betapr = 0.d0
            endif
            if (v .ge. 0.d0) then
                beta0r = 0.d0
            endif

            if (v .gt. 0.d0) then
               spzero = -1.d0
            else
               spzero = v*dtdy
            endif

            if (j .ge. ilo2) then
               qyp(i,j,kc,QRHO)  = 1.d0/(1.d0/rtilde - beta0r - betapr)
               qyp(i,j,kc,QPRES) = p + drght*dp +   clag*clag*betapr
               qyp(i,j,kc,QV)    = v + drght*dv + clag*betapr
               qyp(i,j,kc,QU)    = u + 0.5d0*(-1.d0 - spzero )*du
               qyp(i,j,kc,QW)    = w + 0.5d0*(-1.d0 - spzero )*dw

               ! If rho or p too small, set all the slopes to zero
               if (qyp(i,j,kc,QRHO ) .lt. small_dens .or. &
                   qyp(i,j,kc,QPRES) .lt. small_pres) then
                  qyp(i,j,kc,QRHO)  = rho
                  qyp(i,j,kc,QPRES) = p
                  qyp(i,j,kc,QV)    = v
               end if

               qyp(i,j,kc,QREINT) = qyp(i,j,kc,QPRES) / gamma_minus_1

            end if

            ! **************************************************************************

            dleft  = 0.5d0 - max(0.d0,v + cc)*0.5d0*dtdy
            dleft0 = 0.5d0 * (1.d0-v*dtdy)
            dd0    = dleft0
            dd     = dleft
            rtilde = rho + dleft*drho
            rtr0l  = rho + dleft0*drho
            ptr0l  = p + dleft0*dp

            betaml =  (dp*t1-dv)*t2
            beta0l = -(dp*t3 - drho/(rtilde*rtr0l))*0.5d0*dtdy*cc

            if (v-cc .lt. 0.d0) then
              betaml = 0.d0
            endif
            if (v .lt. 0.d0) then
              beta0l = 0.d0
            endif

            if (v .ge. 0.d0) then 
               spzero = v*dtdy
            else
               spzero = 1.d0
            endif

            if (j .le. ihi2) then

               qym(i,j+1,kc,QRHO)  = 1.d0/(1.d0/rtilde - beta0l  - betaml)
               qym(i,j+1,kc,QPRES) = p + dleft*dp + clag*clag*betaml
               qym(i,j+1,kc,QV)    = v + dleft*du - clag*betaml
               qym(i,j+1,kc,QU)    = u + 0.5d0*(1.d0 - spzero )*du
               qym(i,j+1,kc,QW)    = w + 0.5d0*(1.d0 - spzero )*dw

               ! If rho or p too small, set all the slopes to zero
               if (qym(i,j+1,kc,QRHO ) .lt. small_dens .or. &
                   qym(i,j+1,kc,QPRES) .lt. small_pres) then
                  qym(i,j+1,kc,QRHO)  = rho
                  qym(i,j+1,kc,QPRES) = p
                  qym(i,j+1,kc,QV)    = v
               end if

               qym(i,j+1,kc,QREINT) = qym(i,j+1,kc,QPRES) / gamma_minus_1

            endif

         enddo
      enddo

      ! Do all of the passively advected quantities in one loop
      do ipassive = 1, npassive
         n = qpass_map(ipassive)

         ! Top state
         do j = ilo2, ihi2+1
            do i = ilo1-1, ihi1+1
               v = q(i,j,k3d,QV)
               if (v .gt. 0.d0) then
                  spzero = -1.d0
               else
                  spzero = v*dtdy
               endif
               qyp(i,j,kc,n) = q(i,j,k3d,n) + 0.5d0*(-1.d0 - spzero )*dqy(i,j,kc,n)
            enddo
         end do

         ! Bottom state
         do j = ilo2-1, ihi2
            do i = ilo1-1, ihi1+1
               v = q(i,j,k3d,QV)
               if (v .ge. 0.d0) then
                  spzero = v*dtdy
               else
                  spzero = 1.d0
               endif
               qym(i,j+1,kc,n) = q(i,j,k3d,n) + 0.5d0*(1.d0 - spzero )*dqy(i,j,kc,n)
            enddo
         enddo
      enddo

    end subroutine tracexy_cg

! :::
! ::: ------------------------------------------------------------------
! :::

    subroutine tracez_cg(q,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
         dqz,dq_l1,dq_l2,dq_l3,dq_h1,dq_h2,dq_h3, &
         qzm,qzp,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
         ilo1,ilo2,ihi1,ihi2,dz,dt,km,kc,k3d,a_old)

      use eos_module
      use amrex_fort_module, only : rt => amrex_real
      use meth_params_module, only : QVAR, QRHO, QU, QV, QW, &
                                     QREINT, QPRES, &
                                     npassive, qpass_map, &
                                     small_dens, small_pres, gamma_minus_1

      implicit none

      integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
      integer dq_l1,dq_l2,dq_l3,dq_h1,dq_h2,dq_h3
      integer qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3
      integer ilo1,ilo2,ihi1,ihi2
      integer km,kc,k3d

      real(rt)    q(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      real(rt)    c(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)

      real(rt)  dqz(dq_l1:dq_h1,dq_l2:dq_h2,dq_l3:dq_h3,QVAR)
      real(rt) qzm(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
      real(rt) qzp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
      real(rt) a_old
      real(rt) dz, dt

      ! Local variables
      integer :: i, j, n

      real(rt) :: dtdz
      real(rt) :: cc, csq, rho, u, v, w, p
      real(rt) :: drho, du, dv, dw, dp
      real(rt) :: spzero

      real(rt) :: t1,t2,t3
      real(rt) :: beta0l,beta0r,betaml,betapr
      real(rt) :: clag,dd,dd0,dleft,dleft0,drght,drght0
      real(rt) :: ptr0l,ptr0r,rtilde,rtr0l,rtr0r
      integer          :: ipassive

      real(rt), parameter :: eps = 1.d-3

      dtdz = dt/(dz*a_old)

      !!!!!!!!!!!!!!!
      ! Colella-Glaz algorithm
      !!!!!!!!!!!!!!!

      do j = ilo2-1, ihi2+1
         do i = ilo1-1, ihi1+1

            ! **************************************************************************
            ! This is all for qzp
            ! **************************************************************************

            cc = c(i,j,k3d)
            csq = cc*cc
            rho = q(i,j,k3d,QRHO)
            u   = q(i,j,k3d,QU)
            v   = q(i,j,k3d,QV)
            w   = q(i,j,k3d,QW)
            p   = q(i,j,k3d,QPRES)

            drho = dqz(i,j,kc,QRHO)
            du   = dqz(i,j,kc,QU)
            dv   = dqz(i,j,kc,QV)
            dw   = dqz(i,j,kc,QW)
            dp   = dqz(i,j,kc,QPRES)

            clag = rho*cc
            t1   = 1.0d0 / clag
            t2   = 0.5d0 * dtdz/rho
            t3   = t1*t1

            drght  = -0.5d0 - 0.5d0* min(0.d0,w-cc)*dtdz
            drght0 = -0.5d0 * (1.d0+w*dtdz)
            dd     = drght
            dd0    = drght0
            rtilde = rho + drght*drho
            rtr0r  = rho + drght0*drho
            ptr0r  = p + drght0*dp

            if (w+cc .ge. 0.d0) then
               betapr = 0.d0
            else
               betapr = -(dp*t1+dw)*t2
            endif

            if (w .ge. 0.d0) then
               beta0r = 0.d0
            else
               beta0r =  (dp*t3 -drho/(rtilde*rtr0r)) *0.5d0*dtdz*cc
            endif

            if (w .gt. 0.d0) then
               spzero = -1.d0
            else
               spzero = w*dtdz
            endif

            qzp(i,j,kc,QRHO)  = 1.d0/(1.d0/rtilde - beta0r - betapr)
            qzp(i,j,kc,QPRES) = p + drght*dp + clag*clag*betapr
            qzp(i,j,kc,QW)    = w + drght*dw + clag     *betapr
            qzp(i,j,kc,QU)    = u + 0.5d0*(-1.d0 - spzero)*du
            qzp(i,j,kc,QV)    = v + 0.5d0*(-1.d0 - spzero)*dv

            ! If rho or p too small, set all the slopes to zero
            if (qzp(i,j,kc,QRHO ) .lt. small_dens .or. &
                qzp(i,j,kc,QPRES) .lt. small_pres) then
               qzp(i,j,kc,QRHO)  = rho
               qzp(i,j,kc,QPRES) = p
               qzp(i,j,kc,QW)    = w
            end if

            qzp(i,j,kc,QREINT) = qzp(i,j,kc,QPRES) / gamma_minus_1

            ! **************************************************************************
            ! This is all for qzm
            ! **************************************************************************

            ! repeat above with km (k3d-1) to get qzm at kc
            cc  = c(i,j,k3d-1)
            csq = cc*cc
            rho = q(i,j,k3d-1,QRHO)
            u   = q(i,j,k3d-1,QU)
            v   = q(i,j,k3d-1,QV)
            w   = q(i,j,k3d-1,QW)
            p   = q(i,j,k3d-1,QPRES)

            drho = dqz(i,j,km,QRHO)
            du   = dqz(i,j,km,QU)
            dv   = dqz(i,j,km,QV)
            dw   = dqz(i,j,km,QW)
            dp   = dqz(i,j,km,QPRES)

            clag = rho*cc
            t1   = 1.0d0 / clag
            t2   = 0.5d0 * dtdz/rho
            t3   = t1*t1

            dleft  = 0.5d0 - max(0.d0,w + cc)*0.5d0*dtdz
            dleft0 = 0.5d0 * (1.d0-w*dtdz)
            dd0    = dleft0
            dd     = dleft
            rtilde = rho + dleft*drho
            rtr0l  = rho + dleft0*drho
            ptr0l  = p + dleft0*dp

            if (w-cc .lt. 0.d0) then
               betaml = 0.d0
            else
               betaml =  (dp*t1-dw)*t2
            endif

            if (w .lt. 0.d0) then
               beta0l = 0.d0
            else
               beta0l = -(dp*t3 - drho/(rtilde*rtr0l))*0.5d0*dtdz*cc
            endif

            if (w .ge. 0.d0) then
               spzero = w*dtdz
            else
               spzero = 1.d0
            endif

            qzm(i,j,kc,QRHO)  = 1.d0/(1.d0/rtilde - beta0l  - betaml)
            qzm(i,j,kc,QPRES) = p + dleft*dp + clag*clag*betaml
            qzm(i,j,kc,QW)    = w + dleft*dw - clag     *betaml
            qzm(i,j,kc,QU)    = u + 0.5d0*(1.d0 - spzero )*du
            qzm(i,j,kc,QV)    = v + 0.5d0*(1.d0 - spzero )*dv

            ! If rho or p too small, set all the slopes to zero
            if (qzm(i,j,kc,QRHO ) .lt. small_dens .or. &
                qzm(i,j,kc,QPRES) .lt. small_pres) then
               qzm(i,j,kc,QRHO)  = rho
               qzm(i,j,kc,QPRES) = p
               qzm(i,j,kc,QW)    = w
            end if

            qzm(i,j,kc,QREINT) = qzm(i,j,kc,QPRES) / gamma_minus_1

         enddo
      enddo

      ! Do all of the passively advected quantities in one loop
      do ipassive = 1, npassive
         n = qpass_map(ipassive)

         do j = ilo2-1, ihi2+1
            do i = ilo1-1, ihi1+1

               ! Top state
               w = q(i,j,k3d,QW)
               if (w .gt. 0.d0) then
                  spzero = -1.d0
               else
                  spzero = w*dtdz
               endif
               qzp(i,j,kc,n) = q(i,j,k3d,n) + 0.5d0*(-1.d0 - spzero )*dqz(i,j,kc,n)

               ! Bottom state
               w = q(i,j,k3d-1,QW)
               if (w .ge. 0.d0) then
                  spzero = w*dtdz
               else
                  spzero = 1.d0
               endif
               qzm(i,j,kc,n) = q(i,j,k3d-1,n) + 0.5d0*(1.d0 - spzero )*dqz(i,j,km,n)

            enddo
         enddo
      enddo

    end subroutine tracez_cg
