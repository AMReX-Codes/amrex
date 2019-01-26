! :::
! ::: ------------------------------------------------------------------
! :::
      subroutine tracexy(q,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                         dqx,dqy,dq_l1,dq_l2,dq_l3,dq_h1,dq_h2,dq_h3, &
                         qxm,qxp,qym,qyp,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                         ilo1,ilo2,ihi1,ihi2,dx,dy,dt,kc,k3d,a_old)

      use amrex_error_module
      use amrex_fort_module, only : rt => amrex_real
      use meth_params_module, only : QVAR, QRHO, QU, QV, QW, &
                                     QREINT, QPRES, &
                                     ppm_type, small_dens, small_pres, &
                                     npassive, qpass_map, gamma_minus_1
      use amrex_constants_module
      implicit none

      integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
      integer dq_l1,dq_l2,dq_l3,dq_h1,dq_h2,dq_h3
      integer qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3
      integer ilo1,ilo2,ihi1,ihi2
      integer kc,k3d

      real(rt)     q(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      real(rt)     c(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)

      real(rt)  dqx(dq_l1:dq_h1,dq_l2:dq_h2,dq_l3:dq_h3,QVAR)
      real(rt)  dqy(dq_l1:dq_h1,dq_l2:dq_h2,dq_l3:dq_h3,QVAR)

      real(rt) qxm(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
      real(rt) qxp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
      real(rt) qym(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
      real(rt) qyp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
      real(rt) a_old
      real(rt) dx, dy, dt

      ! Local variables
      integer i, j, n

      real(rt) dtdx, dtdy
      real(rt) cc, csq, rho, u, v, w, p, rhoe
      real(rt) drho, du, dv, dw, dp, drhoe

      real(rt) enth, alpham, alphap, alpha0r, alpha0e
      real(rt) alpha0u, alpha0v, alpha0w
      real(rt) spminus, spplus, spzero
      real(rt) apright, amright, azrright, azeright
      real(rt) azu1rght, azv1rght, azw1rght
      real(rt) apleft, amleft, azrleft, azeleft
      real(rt) :: azu1left, azv1left, azw1left
      integer          :: ipassive

      dtdx = dt/(dx*a_old)
      dtdy = dt/(dy*a_old)

      if (ppm_type .ne. 0) then
        print *,'Oops -- shouldnt be in tracexy with ppm_type != 0'
        call amrex_error("Error:: Nyx_advection_3d.f90 :: tracexy")
      end if

      !!!!!!!!!!!!!!!
      ! NON-PPM CODE
      !!!!!!!!!!!!!!!

      ! Compute left and right traced states
      do j = ilo2-1, ihi2+1
         do i = ilo1-1, ihi1+1

            cc = c(i,j,k3d)
            csq = cc**2
            rho = q(i,j,k3d,QRHO)
            u = q(i,j,k3d,QU)
            v = q(i,j,k3d,QV)
            w = q(i,j,k3d,QW)
            p = q(i,j,k3d,QPRES)
            rhoe = q(i,j,k3d,QREINT)
            enth = ( (rhoe+p)/rho )/csq

            drho = dqx(i,j,kc,QRHO)
            du = dqx(i,j,kc,QU)
            dv = dqx(i,j,kc,QV)
            dw = dqx(i,j,kc,QW)
            dp = dqx(i,j,kc,QPRES)
            drhoe = dqx(i,j,kc,QREINT)

            alpham = HALF*(dp/(rho*cc) - du)*rho/cc
            alphap = HALF*(dp/(rho*cc) + du)*rho/cc
            alpha0r = drho - dp/csq
            alpha0e = drhoe - dp*enth
            alpha0v = dv
            alpha0w = dw

            if (u-cc .gt. ZERO) then
               spminus = -ONE
            else
               spminus = (u-cc)*dtdx
            endif
            if (u+cc .gt. ZERO) then
               spplus = -ONE
            else
               spplus = (u+cc)*dtdx
            endif
            if (u .gt. ZERO) then
               spzero = -ONE
            else
               spzero = u*dtdx
            endif

            apright = HALF*(-ONE - spplus )*alphap
            amright = HALF*(-ONE - spminus)*alpham
            azrright = HALF*(-ONE - spzero )*alpha0r
            azeright = HALF*(-ONE - spzero )*alpha0e
            azv1rght = HALF*(-ONE - spzero )*alpha0v
            azw1rght = HALF*(-ONE - spzero )*alpha0w

            if (i .ge. ilo1) then
               qxp(i,j,kc,QRHO  ) = rho  + apright + amright + azrright
               qxp(i,j,kc,QU    ) = u    + (apright - amright)*cc/rho
               qxp(i,j,kc,QV    ) = v    +  azv1rght
               qxp(i,j,kc,QW    ) = w    +  azw1rght
               qxp(i,j,kc,QPRES ) = p    + (apright + amright)*csq
!              qxp(i,j,kc,QREINT) = rhoe + (apright + amright)*enth*csq + azeright

               ! If rho or p too small, set all the slopes to zero
               if (qxp(i,j,kc,QRHO ) .lt. small_dens .or. &
                   qxp(i,j,kc,QPRES) .lt. small_pres) then
                  qxp(i,j,kc,QPRES) = p
                  qxp(i,j,kc,QRHO)  = rho
                  qxp(i,j,kc,QU)    = u
               end if

               qxp(i,j,kc,QREINT) = qxp(i,j,kc,QPRES) / gamma_minus_1

            end if

            if (u-cc .ge. ZERO) then
               spminus = (u-cc)*dtdx
            else
               spminus = ONE
            endif
            if (u+cc .ge. ZERO) then
               spplus = (u+cc)*dtdx
            else
               spplus = ONE
            endif
            if (u .ge. ZERO) then
               spzero = u*dtdx
            else
               spzero = ONE
            endif

            apleft = HALF*(ONE - spplus )*alphap
            amleft = HALF*(ONE - spminus)*alpham
            azrleft = HALF*(ONE - spzero )*alpha0r
            azeleft = HALF*(ONE - spzero )*alpha0e
            azv1left = HALF*(ONE - spzero )*alpha0v
            azw1left = HALF*(ONE - spzero )*alpha0w

            if (i .le. ihi1) then
               qxm(i+1,j,kc,QRHO  ) = rho  +  apleft + amleft + azrleft
               qxm(i+1,j,kc,QU    ) = u    + (apleft - amleft)*cc/rho
               qxm(i+1,j,kc,QV    ) = v    +  azv1left
               qxm(i+1,j,kc,QW    ) = w    +  azw1left
               qxm(i+1,j,kc,QPRES ) = p    + (apleft + amleft)*csq
            !  qxm(i+1,j,kc,QREINT) = rhoe + (apleft + amleft)*enth*csq + azeleft

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
               if (u .gt. ZERO) then
                  spzero = -ONE
               else
                  spzero = u*dtdx
               endif
               qxp(i,j,kc,n) = q(i,j,k3d,n) + HALF*(-ONE - spzero )*dqx(i,j,kc,n)
            enddo

            ! Left state
            do i = ilo1-1, ihi1
               u = q(i,j,k3d,QU)
               if (u .ge. ZERO) then
                  spzero = u*dtdx
               else
                  spzero = ONE
               endif
               qxm(i+1,j,kc,n) = q(i,j,k3d,n) + HALF*(ONE - spzero )*dqx(i,j,kc,n)
            enddo
         enddo
      enddo

      do j = ilo2-1, ihi2+1
         do i = ilo1-1, ihi1+1

            cc = c(i,j,k3d)
            csq = cc**2
            rho = q(i,j,k3d,QRHO)
            u = q(i,j,k3d,QU)
            v = q(i,j,k3d,QV)
            w = q(i,j,k3d,QW)
            p = q(i,j,k3d,QPRES)
            rhoe = q(i,j,k3d,QREINT)
            enth = ( (rhoe+p)/rho )/csq

            drho = dqy(i,j,kc,QRHO)
            du = dqy(i,j,kc,QU)
            dv = dqy(i,j,kc,QV)
            dw = dqy(i,j,kc,QW)
            dp = dqy(i,j,kc,QPRES)
            drhoe = dqy(i,j,kc,QREINT)

            alpham = HALF*(dp/(rho*cc) - dv)*rho/cc
            alphap = HALF*(dp/(rho*cc) + dv)*rho/cc
            alpha0r = drho - dp/csq
            alpha0e = drhoe - dp*enth
            alpha0u = du
            alpha0w = dw

            if (v-cc .gt. ZERO) then
               spminus = -ONE
            else
               spminus = (v-cc)*dtdy
            endif
            if (v+cc .gt. ZERO) then
               spplus = -ONE
            else
               spplus = (v+cc)*dtdy
            endif
            if (v .gt. ZERO) then
               spzero = -ONE
            else
               spzero = v*dtdy
            endif

            apright = HALF*(-ONE - spplus )*alphap
            amright = HALF*(-ONE - spminus)*alpham
            azrright = HALF*(-ONE - spzero )*alpha0r
            azeright = HALF*(-ONE - spzero )*alpha0e
            azu1rght = HALF*(-ONE - spzero )*alpha0u
            azw1rght = HALF*(-ONE - spzero )*alpha0w

            if (j .ge. ilo2) then
               qyp(i,j,kc,QRHO) = rho + apright + amright + azrright
               qyp(i,j,kc,QV) = v + (apright - amright)*cc/rho
               qyp(i,j,kc,QU) = u + azu1rght
               qyp(i,j,kc,QW) = w + azw1rght
               qyp(i,j,kc,QPRES) = p + (apright + amright)*csq
!              qyp(i,j,kc,QREINT) = rhoe + (apright + amright)*enth*csq + azeright

               ! If rho or p too small, set all the slopes to zero
               if (qyp(i,j,kc,QRHO ) .lt. small_dens .or. &
                   qyp(i,j,kc,QPRES) .lt. small_pres) then
                  qyp(i,j,kc,QRHO)  = rho
                  qyp(i,j,kc,QPRES) = p
                  qyp(i,j,kc,QV)    = v
               end if

               qyp(i,j,kc,QREINT) = qyp(i,j,kc,QPRES) / gamma_minus_1

            end if

            if (v-cc .ge. ZERO) then
               spminus = (v-cc)*dtdy
            else
               spminus = ONE
            endif
            if (v+cc .ge. ZERO) then
               spplus = (v+cc)*dtdy
            else
               spplus = ONE
            endif
            if (v .ge. ZERO) then
               spzero = v*dtdy
            else
               spzero = ONE
            endif

            apleft = HALF*(ONE - spplus )*alphap
            amleft = HALF*(ONE - spminus)*alpham
            azrleft = HALF*(ONE - spzero )*alpha0r
            azeleft = HALF*(ONE - spzero )*alpha0e
            azu1left = HALF*(ONE - spzero )*alpha0u
            azw1left = HALF*(ONE - spzero )*alpha0w

            if (j .le. ihi2) then
               qym(i,j+1,kc,QRHO) = rho + apleft + amleft + azrleft
               qym(i,j+1,kc,QV) = v + (apleft - amleft)*cc/rho
               qym(i,j+1,kc,QU) = u + azu1left
               qym(i,j+1,kc,QW) = w + azw1left
               qym(i,j+1,kc,QPRES) = p + (apleft + amleft)*csq
!              qym(i,j+1,kc,QREINT) = rhoe + (apleft + amleft)*enth*csq + azeleft

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
               if (v .gt. ZERO) then
                  spzero = -ONE
               else
                  spzero = v*dtdy
               endif
               qyp(i,j,kc,n) = q(i,j,k3d,n) + HALF*(-ONE - spzero )*dqy(i,j,kc,n)
            enddo
         end do

         ! Bottom state
         do j = ilo2-1, ihi2
            do i = ilo1-1, ihi1+1
               v = q(i,j,k3d,QV)
               if (v .ge. ZERO) then
                  spzero = v*dtdy
               else
                  spzero = ONE
               endif
               qym(i,j+1,kc,n) = q(i,j,k3d,n) + HALF*(ONE - spzero )*dqy(i,j,kc,n)
            enddo
         enddo
      enddo

    end subroutine tracexy

! :::
! ::: ------------------------------------------------------------------
! :::

      subroutine tracez(q,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
           dqz,dq_l1,dq_l2,dq_l3,dq_h1,dq_h2,dq_h3, &
           qzm,qzp,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
           ilo1,ilo2,ihi1,ihi2,dz,dt,km,kc,k3d,a_old)

      use amrex_constants_module
      use amrex_error_module
      use amrex_fort_module, only : rt => amrex_real
      use meth_params_module, only : QVAR, QRHO, QU, QV, QW, &
                                     QREINT, QPRES, &
                                     ppm_type, small_dens, small_pres, &
                                     npassive, qpass_map, gamma_minus_1

      implicit none

      integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
      integer dq_l1,dq_l2,dq_l3,dq_h1,dq_h2,dq_h3
      integer qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3
      integer ilo1,ilo2,ihi1,ihi2
      integer km,kc,k3d

      real(rt)     q(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      real(rt)     c(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)

      real(rt)  dqz(dq_l1:dq_h1,dq_l2:dq_h2,dq_l3:dq_h3,QVAR)
      real(rt) qzm(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
      real(rt) qzp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
      real(rt) a_old
      real(rt) dz, dt

      ! Local variables
      integer i, j
      integer n

      real(rt) dtdz
      real(rt) cc, csq, rho, u, v, w, p, rhoe

      real(rt) drho, du, dv, dw, dp, drhoe
      real(rt) enth, alpham, alphap, alpha0r, alpha0e
      real(rt) alpha0u, alpha0v
      real(rt) spminus, spplus, spzero
      real(rt) apright, amright, azrright, azeright
      real(rt) azu1rght, azv1rght
      real(rt) apleft, amleft, azrleft, azeleft
      real(rt) azu1left, azv1left

      integer ipassive

      if (ppm_type .ne. 0) then
        print *,'Oops -- shouldnt be in tracez with ppm_type != 0'
        call amrex_error("Error:: Nyx_advection_3d.f90 :: tracez")
      end if

      dtdz = dt/(dz*a_old)

      !!!!!!!!!!!!!!!
      ! NON-PPM CODE
      !!!!!!!!!!!!!!!

      do j = ilo2-1, ihi2+1
         do i = ilo1-1, ihi1+1

            cc = c(i,j,k3d)
            csq = cc**2
            rho = q(i,j,k3d,QRHO)
            u = q(i,j,k3d,QU)
            v = q(i,j,k3d,QV)
            w = q(i,j,k3d,QW)
            p = q(i,j,k3d,QPRES)
            rhoe = q(i,j,k3d,QREINT)
            enth = ( (rhoe+p)/rho )/csq

            drho = dqz(i,j,kc,QRHO)
            du = dqz(i,j,kc,QU)
            dv = dqz(i,j,kc,QV)
            dw = dqz(i,j,kc,QW)
            dp = dqz(i,j,kc,QPRES)
            drhoe = dqz(i,j,kc,QREINT)

            alpham = HALF*(dp/(rho*cc) - dw)*rho/cc
            alphap = HALF*(dp/(rho*cc) + dw)*rho/cc
            alpha0r = drho - dp/csq
            alpha0e = drhoe - dp*enth
            alpha0u = du
            alpha0v = dv

            if (w-cc .gt. ZERO) then
               spminus = -ONE
            else
               spminus = (w-cc)*dtdz
            endif
            if (w+cc .gt. ZERO) then
               spplus = -ONE
            else
               spplus = (w+cc)*dtdz
            endif
            if (w .gt. ZERO) then
               spzero = -ONE
            else
               spzero = w*dtdz
            endif

            apright = HALF*(-ONE - spplus )*alphap
            amright = HALF*(-ONE - spminus)*alpham
            azrright = HALF*(-ONE - spzero )*alpha0r
            azeright = HALF*(-ONE - spzero )*alpha0e
            azu1rght = HALF*(-ONE - spzero )*alpha0u
            azv1rght = HALF*(-ONE - spzero )*alpha0v

            qzp(i,j,kc,QRHO) = rho + apright + amright + azrright
            qzp(i,j,kc,QW) = w + (apright - amright)*cc/rho
            qzp(i,j,kc,QU) = u + azu1rght
            qzp(i,j,kc,QV) = v + azv1rght
            qzp(i,j,kc,QPRES) = p + (apright + amright)*csq
!           qzp(i,j,kc,QREINT) = rhoe + (apright + amright)*enth*csq + azeright

            ! If rho or p too small, set all the slopes to zero
            if (qzp(i,j,kc,QRHO ) .lt. small_dens .or. &
                qzp(i,j,kc,QPRES) .lt. small_pres) then
               qzp(i,j,kc,QRHO)  = rho
               qzp(i,j,kc,QPRES) = p
               qzp(i,j,kc,QW)    = w
            end if
 
            qzp(i,j,kc,QREINT) = qzp(i,j,kc,QPRES) / gamma_minus_1

            ! **************************************************************************

            ! repeat above with km (k3d-1) to get qzm at kc
            cc  = c(i,j,k3d-1)
            csq = cc**2
            rho = q(i,j,k3d-1,QRHO)
            u = q(i,j,k3d-1,QU)
            v = q(i,j,k3d-1,QV)
            w = q(i,j,k3d-1,QW)
            p = q(i,j,k3d-1,QPRES)
            rhoe = q(i,j,k3d-1,QREINT)
            enth = ( (rhoe+p)/rho )/csq

            drho = dqz(i,j,km,QRHO)
            du = dqz(i,j,km,QU)
            dv = dqz(i,j,km,QV)
            dw = dqz(i,j,km,QW)
            dp = dqz(i,j,km,QPRES)
            drhoe = dqz(i,j,km,QREINT)

            alpham = HALF*(dp/(rho*cc) - dw)*rho/cc
            alphap = HALF*(dp/(rho*cc) + dw)*rho/cc
            alpha0r = drho - dp/csq
            alpha0e = drhoe - dp*enth
            alpha0u = du
            alpha0v = dv

            if (w-cc .ge. ZERO) then
               spminus = (w-cc)*dtdz
            else
               spminus = ONE
            endif
            if (w+cc .ge. ZERO) then
               spplus = (w+cc)*dtdz
            else
               spplus = ONE
            endif
            if (w .ge. ZERO) then
               spzero = w*dtdz
            else
               spzero = ONE
            endif

            apleft = HALF*(ONE - spplus )*alphap
            amleft = HALF*(ONE - spminus)*alpham
            azrleft = HALF*(ONE - spzero )*alpha0r
            azeleft = HALF*(ONE - spzero )*alpha0e
            azu1left = HALF*(ONE - spzero )*alpha0u
            azv1left = HALF*(ONE - spzero )*alpha0v

            qzm(i,j,kc,QRHO) = rho + apleft + amleft + azrleft
            qzm(i,j,kc,QW) = w + (apleft - amleft)*cc/rho
            qzm(i,j,kc,QU) = u + azu1left
            qzm(i,j,kc,QV) = v + azv1left
            qzm(i,j,kc,QPRES) = p + (apleft + amleft)*csq
!           qzm(i,j,kc,QREINT) = rhoe + (apleft + amleft)*enth*csq + azeleft

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
               if (w .gt. ZERO) then
                  spzero = -ONE
               else
                  spzero = w*dtdz
               endif
               qzp(i,j,kc,n) = q(i,j,k3d,n) + HALF*(-ONE - spzero )*dqz(i,j,kc,n)

               ! Bottom state
               w = q(i,j,k3d-1,QW)
               if (w .ge. ZERO) then
                  spzero = w*dtdz
               else
                  spzero = ONE
               endif
               qzm(i,j,kc,n) = q(i,j,k3d-1,n) + HALF*(ONE - spzero )*dqz(i,j,km,n)

            enddo
         enddo
      enddo

    end subroutine tracez
