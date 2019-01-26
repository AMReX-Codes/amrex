module trace_ppm_module

  implicit none

  private

  public tracexy_ppm, tracez_ppm

contains

    subroutine tracexy_ppm(q,c,flatn,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                           Ip,Im,Ip_g,Im_g, &
                           qxm,qxp,qym,qyp,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                           srcQ,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3, &
                           ilo1,ilo2,ihi1,ihi2,dt,a_old,kc,k3d)

    use amrex_error_module
    use amrex_fort_module, only : rt => amrex_real
    use meth_params_module, only : QVAR, QRHO, QU, QV, QW, &
                                   QREINT, QPRES, version_2, &
                                   npassive, qpass_map, ppm_type, ppm_reference, &
                                   ppm_flatten_before_integrals, &
                                   small_dens, small_pres, gamma_minus_1
    use amrex_constants_module

    implicit none

    integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
    integer qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3
    integer src_l1,src_l2,src_l3,src_h1,src_h2,src_h3
    integer ilo1,ilo2,ihi1,ihi2
    integer kc,k3d

    real(rt)     q(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
    real(rt)     c(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
    real(rt) flatn(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)

    real(rt)   Ip(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3,QVAR)
    real(rt)   Im(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3,QVAR)

    real(rt) Ip_g(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3,3)
    real(rt) Im_g(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3,3)

    real(rt) qxm (qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
    real(rt) qxp (qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
    real(rt) qym (qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
    real(rt) qyp (qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
    real(rt) srcQ(src_l1:src_h1,src_l2:src_h2,src_l3:src_h3,QVAR)

    real(rt) dt, a_old

    ! Local variables
    integer i, j
    integer n, ipassive

    real(rt) cc, csq, rho, u, v, w, p, rhoe
    real(rt) rho_ref, u_ref, v_ref, w_ref, p_ref, rhoe_ref

    real(rt) drho, du, dv, dw, dp, drhoe
    real(rt) dup, dvp, dpp
    real(rt) dum, dvm, dpm

    real(rt) enth, alpham, alphap, alpha0r, alpha0e
    real(rt) apright, amright, azrright, azeright
    real(rt) apleft, amleft, azrleft, azeleft
    real(rt) xi, xi1
    real(rt) halfdt

    integer, parameter :: igx = 1
    integer, parameter :: igy = 2
    integer, parameter :: igz = 3

    if (ppm_type .eq. 0) then
       print *,'Oops -- shouldnt be in tracexy_ppm with ppm_type = 0'
       call amrex_error("Error:: trace_ppm_3d.f90 :: tracexy_ppm")
    end if

    halfdt = HALF * dt

    !!!!!!!!!!!!!!!
    ! PPM CODE
    !!!!!!!!!!!!!!!

    ! This does the characteristic tracing to build the interface
    ! states using the normal predictor only (no transverse terms).
    
    ! We come in with the Im and Ip arrays -- these are the averages
    ! of the various primitive state variables under the parabolic
    ! interpolant over the region swept out by one of the 3 different
    ! characteristic waves.
    
    ! Im is integrating to the left interface of the current zone
    ! (which will be used to build the right ("p") state at that interface)
    ! and Ip is integrating to the right interface of the current zone
    ! (which will be used to build the left ("m") state at that interface).

    ! The choice of reference state is designed to minimize the
    ! effects of the characteristic projection.  We subtract the I's
    ! off of the reference state, project the quantity such that it is
    ! in terms of the characteristic varaibles, and then add all the
    ! jumps that are moving toward the interface to the reference
    ! state to get the full state on that interface.

    ! Version 2 includes the source terms in the jump in 
    ! the velocities that goes through the characteristic projection.

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

             ! Set the reference state
             if (ppm_reference == 0 .or. &
                (ppm_reference == 1 .and. u - cc >= 0.0d0) ) then
                 rho_ref = rho
                   u_ref = u
                   v_ref = v
                   w_ref = w
                   p_ref = p
                rhoe_ref = rhoe
             else
                 ! This will be the fastest moving state to the left
                 rho_ref = Im(i,j,kc,1,1,QRHO)
                   u_ref = Im(i,j,kc,1,1,QU)
                   v_ref = Im(i,j,kc,1,1,QV)
                   w_ref = Im(i,j,kc,1,1,QW)
                   p_ref = Im(i,j,kc,1,1,QPRES)
                rhoe_ref = Im(i,j,kc,1,1,QREINT)
             endif
   
             ! *m are the jumps carried by u-c
             ! *p are the jumps carried by u+c

             ! Note: for the transverse velocities, the jump is carried
             !       only by the u wave (the contact)
   
             dum   =   u_ref - Im(i,j,kc,1,1,QU)
             dpm   =   p_ref - Im(i,j,kc,1,1,QPRES)
   
             drho  =  rho_ref - Im(i,j,kc,1,2,QRHO)
             dv    =    v_ref - Im(i,j,kc,1,2,QV)
             dw    =    w_ref - Im(i,j,kc,1,2,QW)
             dp    =    p_ref - Im(i,j,kc,1,2,QPRES)
             drhoe = rhoe_ref - Im(i,j,kc,1,2,QREINT)
   
             dup   =    u_ref - Im(i,j,kc,1,3,QU)
             dpp   =    p_ref - Im(i,j,kc,1,3,QPRES)
   
             if (version_2 .eq. 1) then
                 dum = dum - halfdt*srcQ(i,j,k3d,QU)/a_old
                 dup = dup - halfdt*srcQ(i,j,k3d,QU)/a_old
             else if (version_2 .eq. 2) then
                 dum = dum - halfdt*Im_g(i,j,kc,1,1,igx)/a_old
                 dup = dup - halfdt*Im_g(i,j,kc,1,3,igx)/a_old
             end if

            ! These are analogous to the beta's from the original PPM
            ! paper (except we work with rho instead of tau).  This is
            ! simply (l . dq), where dq = qref - I(q)

             alpham = HALF*(dpm/(rho*cc) - dum)*rho/cc
             alphap = HALF*(dpp/(rho*cc) + dup)*rho/cc
             alpha0r = drho - dp/csq
             alpha0e = drhoe - dp*enth

             if (u-cc .gt. ZERO) then
                amright = ZERO
             else if (u-cc .lt. ZERO) then
                amright = -alpham
             else
                amright = -HALF*alpham
             endif

             if (u+cc .gt. ZERO) then
                apright = ZERO
             else if (u+cc .lt. ZERO) then
                apright = -alphap
             else
                apright = -HALF*alphap
             endif

             if (u .gt. ZERO) then
                azrright = ZERO
                azeright = ZERO
             else if (u .lt. ZERO) then
                azrright = -alpha0r
                azeright = -alpha0e
             else
                azrright = -HALF*alpha0r
                azeright = -HALF*alpha0e
             endif

             ! The final interface states are just
             ! q_s = q_ref - sum(l . dq) r
             qxp(i,j,kc,QRHO  ) =  rho_ref +  apright + amright + azrright
             qxp(i,j,kc,QU    ) =    u_ref + (apright - amright)*cc/rho
             qxp(i,j,kc,QREINT) = rhoe_ref + (apright + amright)*enth*csq + azeright
             qxp(i,j,kc,QPRES ) =    p_ref + (apright + amright)*csq

             ! Transverse velocities -- there's no projection here, so we don't
             ! need a reference state.  We only care about the state traced under
             ! the middle wave
             dv = Im(i,j,kc,1,2,QV)
             dw = Im(i,j,kc,1,2,QW)
   
             if (version_2 .eq. 1) then
                dv  = dv  - halfdt*srcQ(i,j,k3d,QV)/a_old
                dw  = dw  - halfdt*srcQ(i,j,k3d,QW)/a_old
             else if (version_2 .eq. 2) then
                dv  = dv  - halfdt*Im_g(i,j,kc,1,2,igy)/a_old
                dw  = dw  - halfdt*Im_g(i,j,kc,1,2,igz)/a_old
             end if

             ! Recall that I already takes the limit of the parabola
             ! in the event that the wave is not moving toward the
             ! interface
             if (u > ZERO) then
                qxp(i,j,kc,QV    ) = v
                qxp(i,j,kc,QW    ) = w
             else ! wave moving toward the interface
                qxp(i,j,kc,QV    ) = dv
                qxp(i,j,kc,QW    ) = dw
             endif

             ! We may have already dealt with the flattening in the construction
             ! of the parabola
             if (ppm_flatten_before_integrals == 0) then
                xi  = flatn(i,j,k3d)
                xi1 = ONE-flatn(i,j,k3d)
 
                qxp(i,j,kc,QRHO  ) = xi1*rho  + xi*qxp(i,j,kc,QRHO  )
                qxp(i,j,kc,QU    ) = xi1*u    + xi*qxp(i,j,kc,QU    )
                qxp(i,j,kc,QV    ) = xi1*v    + xi*qxp(i,j,kc,QV    )
                qxp(i,j,kc,QW    ) = xi1*w    + xi*qxp(i,j,kc,QW    )
                qxp(i,j,kc,QREINT) = xi1*rhoe + xi*qxp(i,j,kc,QREINT)
                qxp(i,j,kc,QPRES ) = xi1*p    + xi*qxp(i,j,kc,QPRES )
             endif

             ! If rho or p too small, set all the slopes to zero
             if (qxp(i,j,kc,QRHO ) .lt. small_dens .or. &
                 qxp(i,j,kc,QPRES) .lt. small_pres) then
                qxp(i,j,kc,QPRES) = p
                qxp(i,j,kc,QRHO)  = rho
                qxp(i,j,kc,QU)    = u
             end if

             qxp(i,j,kc,QREINT) = qxp(i,j,kc,QPRES) / gamma_minus_1

          end if

          ! ******************************************************************************

          if (i .le. ihi1) then

             ! Minus state on face i+1

             ! Set the reference state
             if (ppm_reference == 0 .or. &
                (ppm_reference == 1 .and. u + cc <= 0.0d0) ) then
                 rho_ref = rho
                   u_ref = u
                   v_ref = v
                   w_ref = w
                   p_ref = p
                rhoe_ref = rhoe
             else
                 ! This will be the fastest moving state to the right
                 rho_ref = Ip(i,j,kc,1,3,QRHO)
                   u_ref = Ip(i,j,kc,1,3,QU)
                   v_ref = Ip(i,j,kc,1,3,QV)
                   w_ref = Ip(i,j,kc,1,3,QW)
                   p_ref = Ip(i,j,kc,1,3,QPRES)
                rhoe_ref = Ip(i,j,kc,1,3,QREINT)
             endif
   
             ! *m are the jumps carried by u-c
             ! *p are the jumps carried by u+c
   
             ! Note: for the transverse velocities, the jump is carried
             !       only by the u wave (the contact)

             dum   =    u_ref - Ip(i,j,kc,1,1,QU)
             dpm   =    p_ref - Ip(i,j,kc,1,1,QPRES)
   
             drho  =  rho_ref - Ip(i,j,kc,1,2,QRHO)
             dp    =    p_ref - Ip(i,j,kc,1,2,QPRES)
             drhoe = rhoe_ref - Ip(i,j,kc,1,2,QREINT)

             dup   =    u_ref - Ip(i,j,kc,1,3,QU)
             dpp   =    p_ref - Ip(i,j,kc,1,3,QPRES)
   
             if (version_2 .eq. 1) then
                 dum = dum - halfdt*srcQ(i,j,k3d,QU)/a_old
                 dup = dup - halfdt*srcQ(i,j,k3d,QU)/a_old
             else if (version_2 .eq. 2) then
                 dum = dum - halfdt*Ip_g(i,j,kc,1,1,igx)/a_old
                 dup = dup - halfdt*Ip_g(i,j,kc,1,3,igx)/a_old
             end if

             ! These are analogous to the beta's from the original PPM
             ! paper (except we work with rho instead of tau).  This is
             ! simply (l . dq), where dq = qref - I(q)

             alpham = HALF*(dpm/(rho*cc) - dum)*rho/cc
             alphap = HALF*(dpp/(rho*cc) + dup)*rho/cc
             alpha0r = drho - dp/csq
             alpha0e = drhoe - dp*enth

             if (u-cc .gt. ZERO) then
                amleft = -alpham
             else if (u-cc .lt. ZERO) then
                amleft = ZERO
             else
                amleft = -HALF*alpham
             endif

             if (u+cc .gt. ZERO) then
                apleft = -alphap
             else if (u+cc .lt. ZERO) then
                apleft = ZERO 
             else
                apleft = -HALF*alphap
             endif
   
             if (u .gt. ZERO) then
                azrleft = -alpha0r
                azeleft = -alpha0e
             else if (u .lt. ZERO) then
                azrleft = ZERO
                azeleft = ZERO
             else
                azrleft = -HALF*alpha0r
                azeleft = -HALF*alpha0e
             endif

             ! The final interface states are just
             ! q_s = q_ref - sum (l . dq) r
             qxm(i+1,j,kc,QRHO  ) =  rho_ref + apleft + amleft + azrleft
             qxm(i+1,j,kc,QU    ) =    u_ref + (apleft - amleft)*cc/rho
             qxm(i+1,j,kc,QREINT) = rhoe_ref + (apleft + amleft)*enth*csq + azeleft
             qxm(i+1,j,kc,QPRES ) =    p_ref + (apleft + amleft)*csq

             ! Transverse velocities
             dv    = Ip(i,j,kc,1,2,QV)
             dw    = Ip(i,j,kc,1,2,QW)
   
             if (version_2 .eq. 1) then
                 dv  = dv  - halfdt*srcQ(i,j,k3d,QV)/a_old
                 dw  = dw  - halfdt*srcQ(i,j,k3d,QW)/a_old
             else if (version_2 .eq. 2) then
                 dv  = dv  - halfdt*Ip_g(i,j,kc,1,2,igy)/a_old
                 dw  = dw  - halfdt*Ip_g(i,j,kc,1,2,igz)/a_old
             end if

             if (u < ZERO) then
                qxm(i+1,j,kc,QV    ) = v
                qxm(i+1,j,kc,QW    ) = w
             else
                qxm(i+1,j,kc,QV    ) = dv
                qxm(i+1,j,kc,QW    ) = dw
             endif
 
             ! We may have already dealt with flattening in the parabolas
             if (ppm_flatten_before_integrals == 0) then
                xi  = flatn(i,j,k3d)
                xi1 = ONE - flatn(i,j,k3d)
 
                qxm(i+1,j,kc,QRHO  ) = xi1*rho  + xi*qxm(i+1,j,kc,QRHO  )
                qxm(i+1,j,kc,QU    ) = xi1*u    + xi*qxm(i+1,j,kc,QU    )
                qxm(i+1,j,kc,QV    ) = xi1*v    + xi*qxm(i+1,j,kc,QV    )
                qxm(i+1,j,kc,QW    ) = xi1*w    + xi*qxm(i+1,j,kc,QW    )
                qxm(i+1,j,kc,QREINT) = xi1*rhoe + xi*qxm(i+1,j,kc,QREINT)
                qxm(i+1,j,kc,QPRES ) = xi1*p    + xi*qxm(i+1,j,kc,QPRES )
             endif
 
             ! If rho or p too small, set all the slopes to zero
             if (qxm(i+1,j,kc,QRHO ) .lt. small_dens .or. &
                 qxm(i+1,j,kc,QPRES) .lt. small_pres) then
                qxm(i+1,j,kc,QRHO)  = rho
                qxm(i+1,j,kc,QPRES) = p
                qxm(i+1,j,kc,QU)    = u
             end if

             qxm(i+1,j,kc,QREINT) = qxm(i+1,j,kc,QPRES) / gamma_minus_1
          end if

       end do
    end do

    ! ******************************************************************************
    ! Passively advected quantities 
    ! ******************************************************************************

    ! Do all of the passively advected quantities in one loop
    do ipassive = 1, npassive
       n = qpass_map(ipassive)
       do j = ilo2-1, ihi2+1

          ! Plus state on face i
          do i = ilo1, ihi1+1
             u = q(i,j,k3d,QU)

             if (ppm_flatten_before_integrals == 0) then
                xi = flatn(i,j,k3d)
             else
                xi = ONE
             endif
 
             ! The flattening here is a little confusing.  If
             ! ppm_flatten_before_integrals = 0, then we are blending
             ! the cell centered state and the edge state here through
             ! the flattening procedure.  Otherwise, we've already
             ! took care of flattening.  What we want to do is:
             !
             ! q_l*  (1-xi)*q_i + xi*q_l
             !
             ! where
             !
             ! q_l = q_ref - Proj{(q_ref - I)}
             !
             ! and Proj{} represents the characteristic projection.
             ! But for these, there is only 1-wave that matters, the u
             ! wave, so no projection is needed.  Since we are not
             ! projecting, the reference state doesn't matter, so we
             ! take it to be q_i, therefore, we reduce to
             !
             ! q_l* = (1-xi)*q_i + xi*[q_i - (q_i - I)]
             !      = q_i + xi*(I - q_i)
 
             if (u .gt. ZERO) then
                qxp(i,j,kc,n) = q(i,j,k3d,n)
             else if (u .lt. ZERO) then
                qxp(i,j,kc,n) = q(i,j,k3d,n) + xi*(Im(i,j,kc,1,2,n) - q(i,j,k3d,n))
             else
                qxp(i,j,kc,n) = q(i,j,k3d,n) + HALF*xi*(Im(i,j,kc,1,2,n) - q(i,j,k3d,n))
             endif

	  enddo

          ! Minus state on face i+1
          do i = ilo1-1, ihi1
             u = q(i,j,k3d,QU)
 
             if (ppm_flatten_before_integrals == 0) then
                xi = flatn(i,j,k3d)
             else
                xi = ONE
             endif
 
             if (u .gt. ZERO) then
                qxm(i+1,j,kc,n) = q(i,j,k3d,n) + xi*(Ip(i,j,kc,1,2,n) - q(i,j,k3d,n))
             else if (u .lt. ZERO) then
                qxm(i+1,j,kc,n) = q(i,j,k3d,n)
             else
                qxm(i+1,j,kc,n) = q(i,j,k3d,n) + HALF*xi*(Ip(i,j,kc,1,2,n) - q(i,j,k3d,n))
             endif
          enddo

       enddo
    enddo

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

             ! Set the reference state
             if (ppm_reference == 0 .or. &
                (ppm_reference == 1 .and. v - cc >= 0.0d0) ) then
                 rho_ref = rho
                   u_ref = u
                   v_ref = v
                   w_ref = w
                   p_ref = p
                rhoe_ref = rhoe
             else
                 ! This will be the fastest moving state to the left
                 rho_ref = Im(i,j,kc,2,1,QRHO)
                   u_ref = Im(i,j,kc,2,1,QU)
                   v_ref = Im(i,j,kc,2,1,QV)
                   w_ref = Im(i,j,kc,2,1,QW)
                   p_ref = Im(i,j,kc,2,1,QPRES)
                rhoe_ref = Im(i,j,kc,2,1,QREINT)
             endif
   
             ! *m are the jumps carried by v-c
             ! *p are the jumps carried by v+c
   
             ! Note: for the transverse velocities, the jump is carried
             !       only by the v wave (the contact)

             dvm   =    v_ref - Im(i,j,kc,2,1,QV)
             dpm   =    p_ref - Im(i,j,kc,2,1,QPRES)
   
             drho  =  rho_ref - Im(i,j,kc,2,2,QRHO)
             dp    =    p_ref - Im(i,j,kc,2,2,QPRES)
             drhoe = rhoe_ref - Im(i,j,kc,2,2,QREINT)

             dvp   =    v_ref - Im(i,j,kc,2,3,QV)
             dpp   =    p_ref - Im(i,j,kc,2,3,QPRES)
   
             if (version_2 .eq. 1) then
                 dvm = dvm - halfdt*srcQ(i,j,k3d,QV)/a_old
                 dvp = dvp - halfdt*srcQ(i,j,k3d,QV)/a_old
             else if (version_2 .eq. 2) then
                 dvm = dvm - halfdt*Im_g(i,j,kc,2,1,igy)/a_old
                 dvp = dvp - halfdt*Im_g(i,j,kc,2,3,igy)/a_old
             end if

             ! These are analogous to the beta's from the original PPM
             ! paper (except we work with rho instead of tau).  This
             ! is simply (l . dq), where dq = qref - I(q)
 
             alpham = HALF*(dpm/(rho*cc) - dvm)*rho/cc
             alphap = HALF*(dpp/(rho*cc) + dvp)*rho/cc
             alpha0r = drho - dp/csq
             alpha0e = drhoe - dp*enth
 
             if (v-cc .gt. ZERO) then
                amright = ZERO
             else if (v-cc .lt. ZERO) then
                amright = -alpham
             else
                amright = -HALF*alpham
             endif
 
             if (v+cc .gt. ZERO) then
                apright = ZERO
             else if (v+cc .lt. ZERO) then
                apright = -alphap
             else
                apright = -HALF*alphap
             endif
 
             if (v .gt. ZERO) then
                azrright = ZERO
                azeright = ZERO
             else if (v .lt. ZERO) then
                azrright = -alpha0r
                azeright = -alpha0e
             else
                azrright = -HALF*alpha0r
                azeright = -HALF*alpha0e
             endif
 
             ! The final interface states are just
             ! q_s = q_ref - sum (l . dq) r
             qyp(i,j,kc,QRHO  ) =  rho_ref +  apright + amright + azrright
             qyp(i,j,kc,QV    ) =    v_ref + (apright - amright)*cc/rho
             qyp(i,j,kc,QREINT) = rhoe_ref + (apright + amright)*enth*csq + azeright
             qyp(i,j,kc,QPRES ) =    p_ref + (apright + amright)*csq

             ! Transverse velocities
             du    = Im(i,j,kc,2,2,QU)
             dw    = Im(i,j,kc,2,2,QW)
   
             if (version_2 .eq. 1) then
                 du  = du  - halfdt*srcQ(i,j,k3d,QU)/a_old
                 dw  = dw  - halfdt*srcQ(i,j,k3d,QW)/a_old
             else if (version_2 .eq. 2) then
                 du  = du  - halfdt*Im_g(i,j,kc,2,2,igx)/a_old
                 dw  = dw  - halfdt*Im_g(i,j,kc,2,2,igz)/a_old
             end if
 
             if (v > ZERO) then
                qyp(i,j,kc,QU    ) = u
                qyp(i,j,kc,QW    ) = w
             else ! wave moving toward the interface
                qyp(i,j,kc,QU    ) = du
                qyp(i,j,kc,QW    ) = dw
             endif
 
             ! We may have already dealt with flattening in the parabola
             if (ppm_flatten_before_integrals == 0) then
                xi  = flatn(i,j,k3d)
                xi1 = ONE - flatn(i,j,k3d)
 
                qyp(i,j,kc,QRHO  ) = xi1*rho  + xi*qyp(i,j,kc,QRHO  )
                qyp(i,j,kc,QV    ) = xi1*v    + xi*qyp(i,j,kc,QV    )
                qyp(i,j,kc,QU    ) = xi1*u    + xi*qyp(i,j,kc,QU    )
                qyp(i,j,kc,QW    ) = xi1*w    + xi*qyp(i,j,kc,QW    )
                qyp(i,j,kc,QREINT) = xi1*rhoe + xi*qyp(i,j,kc,QREINT)
                qyp(i,j,kc,QPRES ) = xi1*p    + xi*qyp(i,j,kc,QPRES )
             endif

             ! If rho or p too small, set all the slopes to zero
             if (qyp(i,j,kc,QRHO ) .lt. small_dens .or. &
                 qyp(i,j,kc,QPRES) .lt. small_pres) then
                qyp(i,j,kc,QRHO)  = rho
                qyp(i,j,kc,QPRES) = p
                qyp(i,j,kc,QV)    = v
             end if

             qyp(i,j,kc,QREINT) = qyp(i,j,kc,QPRES) / gamma_minus_1
          end if

          ! ******************************************************************************

          if (j .le. ihi2) then

             ! Minus state on face j+1

             ! Set the reference state
             if (ppm_reference == 0 .or. &
                (ppm_reference == 1 .and. v + cc <= 0.0d0) ) then
                 rho_ref = rho
                   u_ref = u
                   v_ref = v
                   w_ref = w
                   p_ref = p
                rhoe_ref = rhoe
             else
                 ! This will be the fastest moving state to the right
                 rho_ref = Ip(i,j,kc,2,3,QRHO)
                   u_ref = Ip(i,j,kc,2,3,QU)
                   v_ref = Ip(i,j,kc,2,3,QV)
                   w_ref = Ip(i,j,kc,2,3,QW)
                   p_ref = Ip(i,j,kc,2,3,QPRES)
                rhoe_ref = Ip(i,j,kc,2,3,QREINT)
             endif
   
             ! *m are the jumps carried by v-c
             ! *p are the jumps carried by v+c
   
             ! Note: for the transverse velocities, the jump is carried
             !       only by the v wave (the contact)

             dvm   =    v_ref - Ip(i,j,kc,2,1,QV)
             dpm   =    p_ref - Ip(i,j,kc,2,1,QPRES)
   
             drho  =  rho_ref - Ip(i,j,kc,2,2,QRHO)
             dp    =    p_ref - Ip(i,j,kc,2,2,QPRES)
             drhoe = rhoe_ref - Ip(i,j,kc,2,2,QREINT)

             dvp   =    v_ref - Ip(i,j,kc,2,3,QV)
             dpp   =    p_ref - Ip(i,j,kc,2,3,QPRES)

             if (version_2 .eq. 1) then
                 dvm = dvm - halfdt*srcQ(i,j,k3d,QV)/a_old
                 dvp = dvp - halfdt*srcQ(i,j,k3d,QV)/a_old
             else if (version_2 .eq. 2) then
                 dvm = dvm - halfdt*Ip_g(i,j,kc,2,1,igy)/a_old
                 dvp = dvp - halfdt*Ip_g(i,j,kc,2,3,igy)/a_old
             end if

             ! These are analogous to the beta's from the original PPM
             ! paper.  This is simply (l . dq), where dq = qref - I(q)
 
             alpham = HALF*(dpm/(rho*cc) - dvm)*rho/cc
             alphap = HALF*(dpp/(rho*cc) + dvp)*rho/cc
             alpha0r = drho - dp/csq
             alpha0e = drhoe - dp*enth
 
             if (v-cc .gt. ZERO) then
                amleft = -alpham
             else if (v-cc .lt. ZERO) then
                amleft = ZERO
             else
                amleft = -HALF*alpham
             endif

             if (v+cc .gt. ZERO) then
                apleft = -alphap
             else if (v+cc .lt. ZERO) then
                apleft = ZERO
             else
                apleft = -HALF*alphap
             endif

             if (v .gt. ZERO) then
                azrleft = -alpha0r
                azeleft = -alpha0e
             else if (v .lt. ZERO) then
                azrleft = ZERO
                azeleft = ZERO
             else
                azrleft = -HALF*alpha0r
                azeleft = -HALF*alpha0e
             endif

             ! The final interface states are just
             ! q_s = q_ref - sum (l . dq) r
             qym(i,j+1,kc,QRHO  ) =  rho_ref +  apleft + amleft + azrleft
             qym(i,j+1,kc,QV    ) =    v_ref + (apleft - amleft)*cc/rho
             qym(i,j+1,kc,QREINT) = rhoe_ref + (apleft + amleft)*enth*csq + azeleft
             qym(i,j+1,kc,QPRES ) =    p_ref + (apleft + amleft)*csq

             ! Transverse velocities
             du    = Ip(i,j,kc,2,2,QU)
             dw    = Ip(i,j,kc,2,2,QW)

             if (version_2 .eq. 1) then
                 du  = du  - halfdt*srcQ(i,j,k3d,QU)/a_old
                 dw  = dw  - halfdt*srcQ(i,j,k3d,QW)/a_old
             else if (version_2 .eq. 2) then
                 du  = du  - halfdt*Ip_g(i,j,kc,2,2,igx)/a_old
                 dw  = dw  - halfdt*Ip_g(i,j,kc,2,2,igz)/a_old
             end if
 
             if (v < ZERO) then
                qym(i,j+1,kc,QU    ) = u
                qym(i,j+1,kc,QW    ) = w
             else ! wave is moving toward the interface
                qym(i,j+1,kc,QU    ) = du
                qym(i,j+1,kc,QW    ) = dw
             endif

             ! We may have already dealt with flattening in the parabola
             if (ppm_flatten_before_integrals == 0) then
                xi  = flatn(i,j,k3d)
                xi1 = ONE - flatn(i,j,k3d)
 
                qym(i,j+1,kc,QRHO  ) = xi1*rho  + xi*qym(i,j+1,kc,QRHO  )
                qym(i,j+1,kc,QV    ) = xi1*v    + xi*qym(i,j+1,kc,QV    )
                qym(i,j+1,kc,QU    ) = xi1*u    + xi*qym(i,j+1,kc,QU    )
                qym(i,j+1,kc,QW    ) = xi1*w    + xi*qym(i,j+1,kc,QW    )
                qym(i,j+1,kc,QREINT) = xi1*rhoe + xi*qym(i,j+1,kc,QREINT)
                qym(i,j+1,kc,QPRES ) = xi1*p    + xi*qym(i,j+1,kc,QPRES )
             endif

             ! If rho or p too small, set all the slopes to zero
             if (qym(i,j+1,kc,QRHO ) .lt. small_dens .or. &
                 qym(i,j+1,kc,QPRES) .lt. small_pres) then
                qym(i,j+1,kc,QRHO)  = rho
                qym(i,j+1,kc,QPRES) = p
                qym(i,j+1,kc,QV)    = v
             end if

             qym(i,j+1,kc,QREINT) = qym(i,j+1,kc,QPRES) / gamma_minus_1
          end if

       end do
    end do

    !--------------------------------------------------------------------------
    ! Passively advected quantities
    !--------------------------------------------------------------------------

    ! Do all of the passively advected quantities in one loop
    do ipassive = 1, npassive
       n = qpass_map(ipassive)
       do i = ilo1-1, ihi1+1

          ! Plus state on face j
          do j = ilo2, ihi2+1
             v = q(i,j,k3d,QV)

             if (ppm_flatten_before_integrals == 0) then
                xi = flatn(i,j,k3d)
             else
                xi = ONE
             endif

             if (v .gt. ZERO) then
                qyp(i,j,kc,n) = q(i,j,k3d,n)
             else if (v .lt. ZERO) then
                qyp(i,j,kc,n) = q(i,j,k3d,n) + xi*(Im(i,j,kc,2,2,n) - q(i,j,k3d,n))
             else
                qyp(i,j,kc,n) = q(i,j,k3d,n) + HALF*xi*(Im(i,j,kc,2,2,n) - q(i,j,k3d,n))
             endif
          enddo
          
          ! Minus state on face j+1
          do j = ilo2-1, ihi2
             v = q(i,j,k3d,QV)

             if (ppm_flatten_before_integrals == 0) then
                xi = flatn(i,j,k3d)
             else
                xi = ONE
             endif

             if (v .gt. ZERO) then
                qym(i,j+1,kc,n) = q(i,j,k3d,n) + xi*(Ip(i,j,kc,2,2,n) - q(i,j,k3d,n))
             else if (v .lt. ZERO) then
                qym(i,j+1,kc,n) = q(i,j,k3d,n)
             else
                qym(i,j+1,kc,n) = q(i,j,k3d,n) + HALF*xi*(Ip(i,j,kc,2,2,n) - q(i,j,k3d,n))
             endif
          enddo
          
       enddo
    enddo

    end subroutine tracexy_ppm

  ! ::: 
  ! ::: ------------------------------------------------------------------
  ! ::: 

    subroutine tracez_ppm(q,c,flatn,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                          Ip,Im,Ip_g,Im_g, &
                          qzm,qzp,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                          srcQ,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3, &
                          ilo1,ilo2,ihi1,ihi2,dt,a_old,km,kc,k3d)

    use amrex_error_module
    use amrex_fort_module, only : rt => amrex_real
    use meth_params_module, only : QVAR, QRHO, QU, QV, QW, &
                                   QREINT, QPRES, version_2, &
                                   npassive, qpass_map, ppm_type, &
                                   npassive, qpass_map, ppm_type, ppm_reference, &
                                   ppm_flatten_before_integrals, &
                                   small_dens, small_pres, gamma_minus_1
    use amrex_constants_module

    implicit none

    integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
    integer qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3
    integer src_l1,src_l2,src_l3,src_h1,src_h2,src_h3
    integer ilo1,ilo2,ihi1,ihi2
    integer km,kc,k3d

    real(rt)     q(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
    real(rt)     c(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
    real(rt) flatn(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)

    real(rt)   Ip(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3,QVAR)
    real(rt)   Im(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3,QVAR)

    real(rt) Ip_g(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3,3)
    real(rt) Im_g(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3,3)

    real(rt) qzm (qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
    real(rt) qzp (qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
    real(rt) srcQ(src_l1:src_h1,src_l2:src_h2,src_l3:src_h3,QVAR)
    real(rt) dt, a_old

    !     Local variables
    integer i, j
    integer n, ipassive

    real(rt) cc, csq
    real(rt) rho, u, v, w, p, rhoe
    real(rt) rho_ref, u_ref, v_ref, w_ref, p_ref, rhoe_ref
    real(rt) dwp, dpp
    real(rt) dwm, dpm

    real(rt) drho, du, dv, dp, drhoe
    real(rt) enth, alpham, alphap, alpha0r, alpha0e
    real(rt) apright, amright, azrright, azeright
    real(rt) apleft, amleft, azrleft, azeleft
    real(rt) halfdt
    real(rt) xi, xi1

    integer, parameter :: igx = 1
    integer, parameter :: igy = 2
    integer, parameter :: igz = 3

    if (ppm_type .eq. 0) then
       print *,'Oops -- shouldnt be in tracez_ppm with ppm_type = 0'
       call amrex_error("Error:: trace_ppm_3d.f90 :: tracez_ppm")
    end if

    halfdt = HALF * dt

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

          ! Set the reference state
          if (ppm_reference == 0 .or. &
             (ppm_reference == 1 .and. w - cc >= 0.0d0) ) then
              rho_ref = rho
                u_ref = u
                v_ref = v
                w_ref = w
                p_ref = p
             rhoe_ref = rhoe
          else
                 ! This will be the fastest moving state to the left
              rho_ref = Im(i,j,kc,3,1,QRHO)
                u_ref = Im(i,j,kc,3,1,QU)
                v_ref = Im(i,j,kc,3,1,QV)
                w_ref = Im(i,j,kc,3,1,QW)
                p_ref = Im(i,j,kc,3,1,QPRES)
             rhoe_ref = Im(i,j,kc,3,1,QREINT)
          endif

          ! *m are the jumps carried by w-c
          ! *p are the jumps carried by w+c

          ! Note: for the transverse velocities, the jump is carried
          !       only by the w wave (the contact)

          dwm   =    w_ref - Im(i,j,kc,3,1,QW)
          dpm   =    p_ref - Im(i,j,kc,3,1,QPRES)

          drho  =  rho_ref - Im(i,j,kc,3,2,QRHO)
          dp    =    p_ref - Im(i,j,kc,3,2,QPRES)
          drhoe = rhoe_ref - Im(i,j,kc,3,2,QREINT)

          dwp   =    w_ref - Im(i,j,kc,3,3,QW)
          dpp   =    p_ref - Im(i,j,kc,3,3,QPRES)

          if (version_2 .eq. 1) then
              dwm = dwm - halfdt*srcQ(i,j,k3d,QW)/a_old
              dwp = dwp - halfdt*srcQ(i,j,k3d,QW)/a_old
          else if (version_2 .eq. 2) then
              dwm = dwm - halfdt*Im_g(i,j,kc,3,1,igz)/a_old
              dwp = dwp - halfdt*Im_g(i,j,kc,3,3,igz)/a_old
          end if

          ! These are analogous to the beta's from the original PPM
          ! paper.  This is simply (l . dq), where dq = qref - I(q)
          alpham = HALF*(dpm/(rho*cc) - dwm)*rho/cc
          alphap = HALF*(dpp/(rho*cc) + dwp)*rho/cc
          alpha0r = drho - dp/csq
          alpha0e = drhoe - dp*enth

          if (w-cc .gt. ZERO) then
             amright = ZERO
          else if (w-cc .lt. ZERO) then
             amright = -alpham
          else
             amright = -HALF*alpham
          endif
          if (w+cc .gt. ZERO) then
             apright = ZERO
          else if (w+cc .lt. ZERO) then
             apright = -alphap
          else
             apright = -HALF*alphap
          endif
          if (w .gt. ZERO) then
             azrright = ZERO
             azeright = ZERO
          else if (w .lt. ZERO) then
             azrright = -alpha0r
             azeright = -alpha0e
          else
             azrright = -HALF*alpha0r
             azeright = -HALF*alpha0e
          endif

          ! The final interface states are just
          ! q_s = q_ref - sum (l . dq) r
          qzp(i,j,kc,QRHO  ) =  rho_ref +  apright + amright + azrright
          qzp(i,j,kc,QW    ) =    w_ref + (apright - amright)*cc/rho
          qzp(i,j,kc,QREINT) = rhoe_ref + (apright + amright)*enth*csq + azeright
          qzp(i,j,kc,QPRES ) =    p_ref + (apright + amright)*csq

          ! Transverse velocities
          du    = Im(i,j,kc,3,2,QU)
          dv    = Im(i,j,kc,3,2,QV)

          if (version_2 .eq. 1) then
              du  = du  - halfdt*srcQ(i,j,k3d,QU)/a_old
              dv  = dv  - halfdt*srcQ(i,j,k3d,QV)/a_old
          else if (version_2 .eq. 2) then
              du  = du  - halfdt*Im_g(i,j,kc,3,2,igx)/a_old
              dv  = dv  - halfdt*Im_g(i,j,kc,3,2,igy)/a_old
          end if

          if (w > ZERO) then 
             qzp(i,j,kc,QU    ) = u
             qzp(i,j,kc,QV    ) = v
          else ! wave moving toward the interface
             qzp(i,j,kc,QU    ) = du
             qzp(i,j,kc,QV    ) = dv
          endif

          ! We may have already dealt with flattening in the parabola
          if (ppm_flatten_before_integrals == 0) then
             xi  = flatn(i,j,k3d)
             xi1 = ONE - flatn(i,j,k3d)

             qzp(i,j,kc,QRHO  ) = xi1*rho  + xi*qzp(i,j,kc,QRHO  )
             qzp(i,j,kc,QW    ) = xi1*w    + xi*qzp(i,j,kc,QW    )
             qzp(i,j,kc,QU    ) = xi1*u    + xi*qzp(i,j,kc,QU    )
             qzp(i,j,kc,QV    ) = xi1*v    + xi*qzp(i,j,kc,QV    )
             qzp(i,j,kc,QREINT) = xi1*rhoe + xi*qzp(i,j,kc,QREINT)
             qzp(i,j,kc,QPRES ) = xi1*p    + xi*qzp(i,j,kc,QPRES )
          endif

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

          ! Set the reference state
          if (ppm_reference == 0 .or. &
             (ppm_reference == 1 .and. w + cc <= 0.0d0) ) then
              rho_ref = rho
                u_ref = u
                v_ref = v
                w_ref = w
                p_ref = p
             rhoe_ref = rhoe
          else
              ! This will be the fastest moving state to the right
              rho_ref = Ip(i,j,km,3,3,QRHO)
                u_ref = Ip(i,j,km,3,3,QU)
                v_ref = Ip(i,j,km,3,3,QV)
                w_ref = Ip(i,j,km,3,3,QW)
                p_ref = Ip(i,j,km,3,3,QPRES)
             rhoe_ref = Ip(i,j,km,3,3,QREINT)
          endif

          ! *m are the jumps carried by w-c
          ! *p are the jumps carried by w+c

          ! Note: for the transverse velocities, the jump is carried
          !       only by the w wave (the contact)

          dwm   = (   w_ref - Ip(i,j,km,3,1,QW))
          dpm   = (   p_ref - Ip(i,j,km,3,1,QPRES))

          drho  = ( rho_ref - Ip(i,j,km,3,2,QRHO))
          dp    = (   p_ref - Ip(i,j,km,3,2,QPRES))
          drhoe = (rhoe_ref - Ip(i,j,km,3,2,QREINT))

          dwp   = (   w_ref - Ip(i,j,km,3,3,QW))
          dpp   = (   p_ref - Ip(i,j,km,3,3,QPRES))

          if (version_2 .eq. 1) then
              dwm = dwm - halfdt*srcQ(i,j,k3d-1,QW)/a_old
              dwp = dwp - halfdt*srcQ(i,j,k3d-1,QW)/a_old
          else if (version_2 .eq. 2) then
              dwm = dwm - halfdt*Ip_g(i,j,km,3,1,igz)/a_old
              dwp = dwp - halfdt*Ip_g(i,j,km,3,3,igz)/a_old
          end if

          ! These are analogous to the beta's from the original PPM
          ! paper.  This is simply (l . dq), where dq = qref - I(q)

          alpham = HALF*(dpm/(rho*cc) - dwm)*rho/cc
          alphap = HALF*(dpp/(rho*cc) + dwp)*rho/cc
          alpha0r = drho - dp/csq
          alpha0e = drhoe - dp*enth
             
          if (w-cc .gt. ZERO) then
             amleft = -alpham
          else if (w-cc .lt. ZERO) then
             amleft = ZERO
          else
             amleft = -HALF*alpham
          endif
          if (w+cc .gt. ZERO) then
             apleft = -alphap
          else if (w+cc .lt. ZERO) then
             apleft = ZERO
          else
             apleft = -HALF*alphap
          endif
          if (w .gt. ZERO) then
             azrleft = -alpha0r
             azeleft = -alpha0e
          else if (w .lt. ZERO) then
             azrleft = ZERO
             azeleft = ZERO
          else
             azrleft = -HALF*alpha0r
             azeleft = -HALF*alpha0e
          endif
          
          ! The final interface states are just
          ! q_s = q_ref - sum (l . dq) r
          qzm(i,j,kc,QRHO  ) =  rho_ref +  apleft + amleft + azrleft
          qzm(i,j,kc,QW    ) =    w_ref + (apleft - amleft)*cc/rho
          qzm(i,j,kc,QREINT) = rhoe_ref + (apleft + amleft)*enth*csq + azeleft
          qzm(i,j,kc,QPRES ) =    p_ref + (apleft + amleft)*csq

          ! Transverse velocity
          du = Ip(i,j,km,3,2,QU)
          dv = Ip(i,j,km,3,2,QV)

          if (version_2 .eq. 1) then
              du  = du  - halfdt*srcQ(i,j,k3d-1,QU)/a_old
              dv  = dv  - halfdt*srcQ(i,j,k3d-1,QV)/a_old
          else if (version_2 .eq. 2) then
              du  = du  - halfdt*Ip_g(i,j,km,3,2,igx)/a_old
              dv  = dv  - halfdt*Ip_g(i,j,km,3,2,igy)/a_old
          end if
 
          if (w < ZERO) then
             qzm(i,j,kc,QU    ) = u
             qzm(i,j,kc,QV    ) = v
          else ! wave moving toward the interface
             qzm(i,j,kc,QU    ) = du
             qzm(i,j,kc,QV    ) = dv
          endif

          ! We may have already taken care of flattening in the parabola
          if (ppm_flatten_before_integrals == 0) then
             xi  = flatn(i,j,k3d-1)
             xi1 = ONE - flatn(i,j,k3d-1)
 
             qzm(i,j,kc,QRHO  ) = xi1*rho  + xi*qzm(i,j,kc,QRHO  )
             qzm(i,j,kc,QW    ) = xi1*w    + xi*qzm(i,j,kc,QW    )
             qzm(i,j,kc,QU    ) = xi1*u    + xi*qzm(i,j,kc,QU    )
             qzm(i,j,kc,QV    ) = xi1*v    + xi*qzm(i,j,kc,QV    )
             qzm(i,j,kc,QREINT) = xi1*rhoe + xi*qzm(i,j,kc,QREINT)
             qzm(i,j,kc,QPRES ) = xi1*p    + xi*qzm(i,j,kc,QPRES )
 
          endif

          ! If rho or p too small, set all the slopes to zero
          if (qzm(i,j,kc,QRHO ) .lt. small_dens .or. &
              qzm(i,j,kc,QPRES) .lt. small_pres) then
             qzm(i,j,kc,QRHO)  = rho
             qzm(i,j,kc,QPRES) = p
             qzm(i,j,kc,QW)    = w
          end if

          qzm(i,j,kc,QREINT) = qzm(i,j,kc,QPRES) / gamma_minus_1

       end do
    end do

    !--------------------------------------------------------------------------
    ! Passively advected quantities
    !--------------------------------------------------------------------------

    ! Do all of the passively advected quantities in one loop
    do ipassive = 1, npassive
         n = qpass_map(ipassive)
         do j = ilo2-1, ihi2+1
               do i = ilo1-1, ihi1+1

                  ! Plus state on face kc
                  w = q(i,j,k3d,QW)

                  if (ppm_flatten_before_integrals == 0) then
                     xi = flatn(i,j,k3d)
                  else
                     xi = ONE
                  endif

                  if (w .gt. ZERO) then
                     qzp(i,j,kc,n) = q(i,j,k3d,n)
                  else if (w .lt. ZERO) then
                     qzp(i,j,kc,n) = q(i,j,k3d,n) + xi*(Im(i,j,kc,3,2,n) - q(i,j,k3d,n))
                  else
                     qzp(i,j,kc,n) = q(i,j,k3d,n) + HALF*xi*(Im(i,j,kc,3,2,n) - q(i,j,k3d,n))
                  endif

                  ! Minus state on face k
                  w = q(i,j,k3d-1,QW)
                  
                  if (ppm_flatten_before_integrals == 0) then
                     xi = flatn(i,j,k3d-1)
                  else
                     xi = ONE
                  endif

                  if (w .gt. ZERO) then
                     qzm(i,j,kc,n) = q(i,j,k3d-1,n) + xi*(Ip(i,j,km,3,2,n) - q(i,j,k3d-1,n))
                  else if (w .lt. ZERO) then
                     qzm(i,j,kc,n) = q(i,j,k3d-1,n)
                  else
                     qzm(i,j,kc,n) = q(i,j,k3d-1,n) + HALF*xi*(Ip(i,j,km,3,2,n) - q(i,j,k3d-1,n))
                  endif

               enddo
         enddo
    enddo

    end subroutine tracez_ppm

end module trace_ppm_module
