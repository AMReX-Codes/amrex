module trace_src_module

  implicit none

  private

  public tracex_src, tracey_src, tracez_src

contains

    subroutine tracex_src(q,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                          qxm,qxp,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                          srcQ,srcq_l1,srcq_l2,srcq_l3,srcq_h1,srcq_h2,srcq_h3, &
                          ilo1,ilo2,ihi1,ihi2,dt,a_old,kc,k3d)

    use amrex_fort_module, only : rt => amrex_real
    use meth_params_module, only : QVAR, QRHO, QU, QV, QW, &
                                   QREINT, QPRES, gamma_minus_1

    implicit none

    integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
    integer qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3
    integer srcq_l1,srcq_l2,srcq_l3,srcq_h1,srcq_h2,srcq_h3
    integer ilo1,ilo2,ihi1,ihi2
    integer kc,k3d

    real(rt)     q(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
    real(rt)     c(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)

    real(rt) qxm ( qpd_l1: qpd_h1, qpd_l2: qpd_h2, qpd_l3: qpd_h3,QVAR)
    real(rt) qxp ( qpd_l1: qpd_h1, qpd_l2: qpd_h2, qpd_l3: qpd_h3,QVAR)
    real(rt) srcQ(srcq_l1:srcq_h1,srcq_l2:srcq_h2,srcq_l3:srcq_h3,QVAR)

    real(rt) dt, a_old

    ! Local variables
    integer i, j

    real(rt) cc, csq, rho, u, v, w, p, rhoe

    real(rt) drho, dv, dw, dp, drhoe
    real(rt) dup, dpp
    real(rt) dum, dpm

    real(rt) enth, alpham, alphap, alpha0r, alpha0e
    real(rt) alpha0v, alpha0w
    real(rt) apright, amright, azrright, azeright
    real(rt) azv1rght, azw1rght
    real(rt) apleft, amleft, azrleft, azeleft
    real(rt) azv1left, azw1left
    real(rt) halfdt

    halfdt = 0.5d0 * dt

    ! *********************************************************************************************
    ! x-direction
    ! *********************************************************************************************

    ! Trace to left and right edges using upwind PPM
    do j = ilo2-1, ihi2+1
       do i = ilo1-1, ihi1+1

          rho  = q(i,j,k3d,QRHO)
          u    = q(i,j,k3d,QU)
          v    = q(i,j,k3d,QV)
          w    = q(i,j,k3d,QW)
          p    = q(i,j,k3d,QPRES)
          rhoe = q(i,j,k3d,QREINT)

          cc = c(i,j,k3d)
          csq = cc**2
          enth = ( (rhoe+p)/rho )/csq

          ! ******************************************************************************

          if (i .ge. ilo1) then

             ! Plus state on face i

             ! *m are the jumps carried by u-c
             ! *p are the jumps carried by u+c

             ! Note: for the transverse velocities, the jump is carried
             !       only by the u wave (the contact)
   
             dum    = -halfdt*(srcQ(i,j,k3d,QU)/a_old)
             dpm    = 0.d0
   
             drho  = 0.d0
             dv    = -halfdt*(srcQ(i,j,k3d,QV)/a_old)
             dw    = -halfdt*(srcQ(i,j,k3d,QW)/a_old)
             dp    = 0.d0
             drhoe = 0.d0
   
             dup    = -halfdt*(srcQ(i,j,k3d,QU)/a_old)
             dpp    = 0.d0
   
             alpham = 0.5d0*(dpm/(rho*cc) - dum)*rho/cc
             alphap = 0.5d0*(dpp/(rho*cc) + dup)*rho/cc
             alpha0r = drho - dp/csq
             alpha0e = drhoe - dp*enth
             alpha0v = dv
             alpha0w = dw

             if (u-cc .gt. 0.d0) then
                amright = 0.d0
             else if (u-cc .lt. 0.d0) then
                amright = -alpham
             else
                amright = -0.5d0*alpham
             endif

             if (u+cc .gt. 0.d0) then
                apright = 0.d0
             else if (u+cc .lt. 0.d0) then
                apright = -alphap
             else
                apright = -0.5d0*alphap
             endif

             if (u .gt. 0.d0) then
                azrright = 0.d0
                azeright = 0.d0
                azv1rght = 0.d0
                azw1rght = 0.d0
             else if (u .lt. 0.d0) then
                azrright = -alpha0r
                azeright = -alpha0e
                azv1rght = -alpha0v
                azw1rght = -alpha0w
             else
                azrright = -0.5d0*alpha0r
                azeright = -0.5d0*alpha0e
                azv1rght = -0.5d0*alpha0v
                azw1rght = -0.5d0*alpha0w
             endif

             qxp(i,j,kc,QRHO  ) = qxp(i,j,kc,QRHO  ) + apright + amright + azrright
             qxp(i,j,kc,QU    ) = qxp(i,j,kc,QU    ) + (apright - amright)*cc/rho
             qxp(i,j,kc,QV    ) = qxp(i,j,kc,QV    ) + azv1rght
             qxp(i,j,kc,QW    ) = qxp(i,j,kc,QW    ) + azw1rght
             qxp(i,j,kc,QPRES ) = qxp(i,j,kc,QPRES ) + (apright + amright)*csq

!            qxp(i,j,kc,QREINT) = qxp(i,j,kc,QREINT) + (apright + amright)*enth*csq + azeright
             qxp(i,j,kc,QREINT) = qxp(i,j,kc,QPRES) / gamma_minus_1
          end if

          ! ******************************************************************************

          if (i .le. ihi1) then

             ! Minus state on face i+1

             ! *m are the jumps carried by u-c
             ! *p are the jumps carried by u+c
   
             ! Note: for the transverse velocities, the jump is carried
             !       only by the u wave (the contact)

             dum   = -halfdt*(srcQ(i,j,k3d,QU)/a_old)
             dpm    = 0.d0
   
             drho  = 0.d0
             dv    = -halfdt*srcQ(i,j,k3d,QV)/a_old
             dw    = -halfdt*srcQ(i,j,k3d,QW)/a_old
             dp    = 0.d0
             drhoe = 0.d0

             dup   = -halfdt*srcQ(i,j,k3d,QU)/a_old
             dpp   = 0.d0
   
             alpham = 0.5d0*(dpm/(rho*cc) - dum)*rho/cc
             alphap = 0.5d0*(dpp/(rho*cc) + dup)*rho/cc
             alpha0r = drho - dp/csq
             alpha0e = drhoe - dp*enth
             alpha0v = dv
             alpha0w = dw

             if (u-cc .gt. 0.d0) then
                amleft = -alpham
             else if (u-cc .lt. 0.d0) then
                amleft = 0.d0
             else
                amleft = -0.5d0*alpham
             endif

             if (u+cc .gt. 0.d0) then
                apleft = -alphap
             else if (u+cc .lt. 0.d0) then
                apleft = 0.d0
             else
                apleft = -0.5d0*alphap
             endif
   
             if (u .gt. 0.d0) then
                azrleft = -alpha0r
                azeleft = -alpha0e
                azv1left = -alpha0v
                azw1left = -alpha0w
             else if (u .lt. 0.d0) then
                azrleft = 0.d0
                azeleft = 0.d0
                azv1left = 0.d0
                azw1left = 0.d0
             else
                azrleft = -0.5d0*alpha0r
                azeleft = -0.5d0*alpha0e
                azv1left = -0.5d0*alpha0v
                azw1left = -0.5d0*alpha0w
             endif

             qxm(i+1,j,kc,QRHO  ) = qxm(i+1,j,kc,QRHO  ) + apleft + amleft + azrleft
             qxm(i+1,j,kc,QU    ) = qxm(i+1,j,kc,QU    ) + (apleft - amleft)*cc/rho
             qxm(i+1,j,kc,QV    ) = qxm(i+1,j,kc,QV    ) + azv1left
             qxm(i+1,j,kc,QW    ) = qxm(i+1,j,kc,QW    ) + azw1left
             qxm(i+1,j,kc,QPRES ) = qxm(i+1,j,kc,QPRES ) + (apleft + amleft)*csq

!            qxm(i+1,j,kc,QREINT) = qxm(i+1,j,kc,QREINT) + (apleft + amleft)*enth*csq + azeleft
             qxm(i+1,j,kc,QREINT) = qxm(i+1,j,kc,QPRES) / gamma_minus_1
          end if

       end do
    end do

    end subroutine tracex_src

    subroutine tracey_src(q,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                          qym,qyp,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                          srcQ,srcq_l1,srcq_l2,srcq_l3,srcq_h1,srcq_h2,srcq_h3, &
                          ilo1,ilo2,ihi1,ihi2,dt,a_old,kc,k3d)

    use amrex_fort_module, only : rt => amrex_real
    use meth_params_module, only : QVAR, QRHO, QU, QV, QW, &
                                   QREINT, QPRES, gamma_minus_1

    implicit none

    integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
    integer qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3
    integer srcq_l1,srcq_l2,srcq_l3,srcq_h1,srcq_h2,srcq_h3
    integer ilo1,ilo2,ihi1,ihi2
    integer kc,k3d

    real(rt)     q(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
    real(rt)     c(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)

    real(rt) qym ( qpd_l1: qpd_h1, qpd_l2: qpd_h2, qpd_l3: qpd_h3,QVAR)
    real(rt) qyp ( qpd_l1: qpd_h1, qpd_l2: qpd_h2, qpd_l3: qpd_h3,QVAR)
    real(rt) srcQ(srcq_l1:srcq_h1,srcq_l2:srcq_h2,srcq_l3:srcq_h3,QVAR)

    real(rt) dt, a_old

    ! Local variables
    integer i, j

    real(rt) cc, csq, rho, u, v, w, p, rhoe

    real(rt) drho, du, dw, dp, drhoe
    real(rt) dvp, dpp
    real(rt) dvm, dpm

    real(rt) enth, alpham, alphap, alpha0r, alpha0e
    real(rt) alpha0w, alpha0u
    real(rt) apright, amright, azrright, azeright
    real(rt) azu1rght, azw1rght
    real(rt) apleft, amleft, azrleft, azeleft
    real(rt) azu1left, azw1left
    real(rt) halfdt

    halfdt = 0.5d0 * dt

    ! *********************************************************************************************
    ! y-direction
    ! *********************************************************************************************

    ! Trace to bottom and top edges using upwind PPM
    do j = ilo2-1, ihi2+1
       do i = ilo1-1, ihi1+1

          rho = q(i,j,k3d,QRHO)
          u = q(i,j,k3d,QU)
          v = q(i,j,k3d,QV)
          w = q(i,j,k3d,QW)
          p = q(i,j,k3d,QPRES)
          rhoe = q(i,j,k3d,QREINT)

          cc = c(i,j,k3d)
          csq = cc**2
          enth = ( (rhoe+p)/rho )/csq

          if (j .ge. ilo2) then

             ! Plus state on face j

             ! *m are the jumps carried by v-c
             ! *p are the jumps carried by v+c
   
             ! Note: for the transverse velocities, the jump is carried
             !       only by the v wave (the contact)

             dvm   = -halfdt*(srcQ(i,j,k3d,QV)/a_old)
             dpm   = 0.d0
   
             drho  = 0.d0
             du    = -halfdt*(srcQ(i,j,k3d,QU)/a_old)
             dw    = -halfdt*(srcQ(i,j,k3d,QW)/a_old)
             dp    = 0.d0
             drhoe = 0.d0

             dvp   = -halfdt*(srcQ(i,j,k3d,QV)/a_old)
             dpp   = 0.d0
   
             alpham = 0.5d0*(dpm/(rho*cc) - dvm)*rho/cc
             alphap = 0.5d0*(dpp/(rho*cc) + dvp)*rho/cc
             alpha0r = drho - dp/csq
             alpha0e = drhoe - dp*enth
             alpha0u = du
             alpha0w = dw

             if (v-cc .gt. 0.d0) then
                amright = 0.d0
             else if (v-cc .lt. 0.d0) then
                amright = -alpham
             else
                amright = -0.5d0*alpham
             endif

             if (v+cc .gt. 0.d0) then
                apright = 0.d0
             else if (v+cc .lt. 0.d0) then
                apright = -alphap
             else
                apright = -0.5d0*alphap
             endif

             if (v .gt. 0.d0) then
                azrright = 0.d0
                azeright = 0.d0
                azu1rght = 0.d0
                azw1rght = 0.d0
             else if (v .lt. 0.d0) then
                azrright = -alpha0r
                azeright = -alpha0e
                azu1rght = -alpha0u
                azw1rght = -alpha0w
             else
                azrright = -0.5d0*alpha0r
                azeright = -0.5d0*alpha0e
                azu1rght = -0.5d0*alpha0u
                azw1rght = -0.5d0*alpha0w
             endif

             qyp(i,j,kc,QRHO  ) = qyp(i,j,kc,QRHO  ) + apright + amright + azrright
             qyp(i,j,kc,QV    ) = qyp(i,j,kc,QV    ) + (apright - amright)*cc/rho
             qyp(i,j,kc,QU    ) = qyp(i,j,kc,QU    ) + azu1rght
             qyp(i,j,kc,QW    ) = qyp(i,j,kc,QW    ) + azw1rght
             qyp(i,j,kc,QPRES ) = qyp(i,j,kc,QPRES ) + (apright + amright)*csq

!            qyp(i,j,kc,QREINT) = qyp(i,j,kc,QREINT) + (apright + amright)*enth*csq + azeright
             qyp(i,j,kc,QREINT) = qyp(i,j,kc,QPRES) / gamma_minus_1
          end if

          ! ******************************************************************************

          if (j .le. ihi2) then

             ! Minus state on face j+1

             ! *m are the jumps carried by v-c
             ! *p are the jumps carried by v+c
   
             ! Note: for the transverse velocities, the jump is carried
             !       only by the v wave (the contact)

             dvm   = -halfdt*(srcQ(i,j,k3d,QV)/a_old)
             dpm   = 0.d0
   
             drho  = 0.d0
             du    = -halfdt*(srcQ(i,j,k3d,QU)/a_old)
             dw    = -halfdt*(srcQ(i,j,k3d,QW)/a_old)
             dp    = 0.d0
             drhoe = 0.d0

             dvp   = -halfdt*(srcQ(i,j,k3d,QV)/a_old)
             dpp   = 0.d0

             alpham = 0.5d0*(dpm/(rho*cc) - dvm)*rho/cc
             alphap = 0.5d0*(dpp/(rho*cc) + dvp)*rho/cc
             alpha0r = drho - dp/csq
             alpha0e = drhoe - dp*enth
             alpha0u = du
             alpha0w = dw
   
             if (v-cc .gt. 0.d0) then
                amleft = -alpham
             else if (v-cc .lt. 0.d0) then
                amleft = 0.d0
             else
                amleft = -0.5d0*alpham
             endif

             if (v+cc .gt. 0.d0) then
                apleft = -alphap
             else if (v+cc .lt. 0.d0) then
                apleft = 0.d0
             else
                apleft = -0.5d0*alphap
             endif

             if (v .gt. 0.d0) then
                azrleft = -alpha0r
                azeleft = -alpha0e
                azu1left = -alpha0u
                azw1left = -alpha0w
             else if (v .lt. 0.d0) then
                azrleft = 0.d0
                azeleft = 0.d0
                azu1left = 0.d0
                azw1left = 0.d0
             else
                azrleft = -0.5d0*alpha0r
                azeleft = -0.5d0*alpha0e
                azu1left = -0.5d0*alpha0u
                azw1left = -0.5d0*alpha0w
             endif

             qym(i,j+1,kc,QRHO  ) = qym(i,j+1,kc,QRHO  )+ apleft + amleft + azrleft
             qym(i,j+1,kc,QV    ) = qym(i,j+1,kc,QV    )+ (apleft - amleft)*cc/rho
             qym(i,j+1,kc,QU    ) = qym(i,j+1,kc,QU    )+ azu1left
             qym(i,j+1,kc,QW    ) = qym(i,j+1,kc,QW    )+ azw1left
             qym(i,j+1,kc,QPRES ) = qym(i,j+1,kc,QPRES )+ (apleft + amleft)*csq

!            qym(i,j+1,kc,QREINT) = qym(i,j+1,kc,QREINT)+ (apleft + amleft)*enth*csq + azeleft
             qym(i,j+1,kc,QREINT) = qym(i,j+1,kc,QPRES) / gamma_minus_1

          end if

       end do
    end do

    end subroutine tracey_src

  ! ::: 
  ! ::: ------------------------------------------------------------------
  ! ::: 

    subroutine tracez_src(q,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                          qzm,qzp,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                          srcQ,srcq_l1,srcq_l2,srcq_l3,srcq_h1,srcq_h2,srcq_h3, &
                          ilo1,ilo2,ihi1,ihi2,dt,a_old,km,kc,k3d)

    use amrex_fort_module, only : rt => amrex_real
    use meth_params_module, only : QVAR, QRHO, QU, QV, QW, &
                                   QREINT, QPRES, gamma_minus_1

    implicit none

    integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
    integer qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3
    integer srcq_l1,srcq_l2,srcq_l3,srcq_h1,srcq_h2,srcq_h3
    integer ilo1,ilo2,ihi1,ihi2
    integer km,kc,k3d

    real(rt)     q(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
    real(rt)     c(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)

    real(rt) qzm ( qpd_l1: qpd_h1, qpd_l2: qpd_h2, qpd_l3: qpd_h3,QVAR)
    real(rt) qzp ( qpd_l1: qpd_h1, qpd_l2: qpd_h2, qpd_l3: qpd_h3,QVAR)
    real(rt) srcQ(srcq_l1:srcq_h1,srcq_l2:srcq_h2,srcq_l3:srcq_h3,QVAR)

    real(rt) dt, a_old

    !     Local variables
    integer i, j

    real(rt) cc, csq
    real(rt) rho, u, v, w, p, rhoe
    real(rt) dwp, dpp
    real(rt) dwm, dpm

    real(rt) drho, du, dv, dp, drhoe
    real(rt) enth, alpham, alphap, alpha0r, alpha0e
    real(rt) alpha0u, alpha0v
    real(rt) apright, amright, azrright, azeright
    real(rt) azu1rght, azv1rght
    real(rt) apleft, amleft, azrleft, azeleft
    real(rt) azu1left, azv1left
    real(rt) halfdt

    halfdt = 0.5d0 * dt

    !!!!!!!!!!!!!!!
    ! PPM CODE
    !!!!!!!!!!!!!!!

    ! Trace to left and right edges using upwind PPM

    ! Note: in contrast to the above code for x and y, here the loop
    ! is over interfaces, not over cell-centers.

    do j = ilo2-1, ihi2+1
       do i = ilo1-1, ihi1+1

          ! **************************************************************************
          ! This is all for qzp
          ! **************************************************************************

          rho  = q(i,j,k3d,QRHO)
          u    = q(i,j,k3d,QU)
          v    = q(i,j,k3d,QV)
          w    = q(i,j,k3d,QW)
          p    = q(i,j,k3d,QPRES)
          rhoe = q(i,j,k3d,QREINT)

          cc   = c(i,j,k3d)
          csq  = cc**2
          enth = ( (rhoe+p)/rho )/csq

          ! Plus state on face kc

          ! *m are the jumps carried by w-c
          ! *p are the jumps carried by w+c

          ! Note: for the transverse velocities, the jump is carried
          !       only by the w wave (the contact)

          dwm   = -halfdt*(srcQ(i,j,k3d,QW)/a_old)
          dpm   = 0.d0 

          drho  =  0.d0
          du    = -halfdt*(srcQ(i,j,k3d,QU)/a_old)
          dv    = -halfdt*(srcQ(i,j,k3d,QV)/a_old)
          dp    =  0.d0
          drhoe =  0.d0

          dwp   = -halfdt*(srcQ(i,j,k3d,QW)/a_old)
          dpp   =  0.d0

          alpham = 0.5d0*(dpm/(rho*cc) - dwm)*rho/cc
          alphap = 0.5d0*(dpp/(rho*cc) + dwp)*rho/cc
          alpha0r = drho - dp/csq
          alpha0e = drhoe - dp*enth
          alpha0u = du
          alpha0v = dv

          if (w-cc .gt. 0.d0) then
             amright = 0.d0
          else if (w-cc .lt. 0.d0) then
             amright = -alpham
          else
             amright = -0.5d0*alpham
          endif
          if (w+cc .gt. 0.d0) then
             apright = 0.d0
          else if (w+cc .lt. 0.d0) then
             apright = -alphap
          else
             apright = -0.5d0*alphap
          endif
          if (w .gt. 0.d0) then
             azrright = 0.d0
             azeright = 0.d0
             azu1rght = 0.d0
             azv1rght = 0.d0
          else if (w .lt. 0.d0) then
             azrright = -alpha0r
             azeright = -alpha0e
             azu1rght = -alpha0u
             azv1rght = -alpha0v
          else
             azrright = -0.5d0*alpha0r
             azeright = -0.5d0*alpha0e
             azu1rght = -0.5d0*alpha0u
             azv1rght = -0.5d0*alpha0v
          endif

          qzp(i,j,kc,QRHO  ) = qzp(i,j,kc,QRHO  ) + apright + amright + azrright
          qzp(i,j,kc,QW    ) = qzp(i,j,kc,QW    ) + (apright - amright)*cc/rho
          qzp(i,j,kc,QU    ) = qzp(i,j,kc,QU    ) + azu1rght
          qzp(i,j,kc,QV    ) = qzp(i,j,kc,QV    ) + azv1rght
          qzp(i,j,kc,QPRES ) = qzp(i,j,kc,QPRES ) + (apright + amright)*csq

!         qzp(i,j,kc,QREINT) = qzp(i,j,kc,QREINT) + (apright + amright)*enth*csq + azeright
          qzp(i,j,kc,QREINT) = qzp(i,j,kc,QPRES) / gamma_minus_1

          ! **************************************************************************
          ! This is all for qzm
          ! **************************************************************************

          ! Minus state on face kc

          ! Note this is different from how we do 1D, 2D, and the
          ! x and y-faces in 3D, where the analogous thing would have
          ! been to find the minus state on face kc+1

          rho  = q(i,j,k3d-1,QRHO)
          u    = q(i,j,k3d-1,QU)
          v    = q(i,j,k3d-1,QV)
          w    = q(i,j,k3d-1,QW)
          p    = q(i,j,k3d-1,QPRES)
          rhoe = q(i,j,k3d-1,QREINT)

          cc   = c(i,j,k3d-1)
          csq  = cc**2
          enth = ( (rhoe+p)/rho )/csq

          ! *m are the jumps carried by w-c
          ! *p are the jumps carried by w+c

          ! Note: for the transverse velocities, the jump is carried
          !       only by the w wave (the contact)

          dwm   = -halfdt*(srcQ(i,j,k3d-1,QW)/a_old)
          dpm   = 0.d0

          drho  = 0.d0
          du    = -halfdt*(srcQ(i,j,k3d-1,QU)/a_old)
          dv    = -halfdt*(srcQ(i,j,k3d-1,QV)/a_old)
          dp    = 0.d0
          drhoe = 0.d0

          dwp   = -halfdt*(srcQ(i,j,k3d-1,QW)/a_old)
          dpp   = 0.d0

          alpham = 0.5d0*(dpm/(rho*cc) - dwm)*rho/cc
          alphap = 0.5d0*(dpp/(rho*cc) + dwp)*rho/cc
          alpha0r = drho - dp/csq
          alpha0e = drhoe - dp*enth
          alpha0u = du
          alpha0v = dv

          if (w-cc .gt. 0.d0) then
             amleft = -alpham
          else if (w-cc .lt. 0.d0) then
             amleft = 0.d0
          else
             amleft = -0.5d0*alpham
          endif
          if (w+cc .gt. 0.d0) then
             apleft = -alphap
          else if (w+cc .lt. 0.d0) then
             apleft = 0.d0
          else
             apleft = -0.5d0*alphap
          endif
          if (w .gt. 0.d0) then
             azrleft = -alpha0r
             azeleft = -alpha0e
             azu1left = -alpha0u
             azv1left = -alpha0v
          else if (w .lt. 0.d0) then
             azrleft = 0.d0
             azeleft = 0.d0
             azu1left = 0.d0
             azv1left = 0.d0
          else
             azrleft = -0.5d0*alpha0r
             azeleft = -0.5d0*alpha0e
             azu1left = -0.5d0*alpha0u
             azv1left = -0.5d0*alpha0v
          endif

          qzm(i,j,kc,QRHO  ) = qzm(i,j,kc,QRHO  ) + apleft + amleft + azrleft
          qzm(i,j,kc,QW    ) = qzm(i,j,kc,QW    ) + (apleft - amleft)*cc/rho
          qzm(i,j,kc,QU    ) = qzm(i,j,kc,QU    ) + azu1left
          qzm(i,j,kc,QV    ) = qzm(i,j,kc,QV    ) + azv1left
          qzm(i,j,kc,QPRES ) = qzm(i,j,kc,QPRES ) + (apleft + amleft)*csq

!         qzm(i,j,kc,QREINT) = qzm(i,j,kc,QREINT) + (apleft + amleft)*enth*csq + azeleft
          qzm(i,j,kc,QREINT) = qzm(i,j,kc,QPRES) / gamma_minus_1

       end do
    end do

    end subroutine tracez_src

end module trace_src_module
