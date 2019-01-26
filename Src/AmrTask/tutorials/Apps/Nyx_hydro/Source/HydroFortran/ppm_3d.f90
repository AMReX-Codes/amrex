module ppm_module

  implicit none

  private

  public ppm

contains
  !
  ! characteristics based on u
  !
  subroutine ppm(s,s_l1,s_l2,s_l3,s_h1,s_h2,s_h3, &
                 u,cspd,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                 flatn,f_l1,f_l2,f_l3,f_h1,f_h2,f_h3, &
                 Ip,Im, &
                 ilo1,ilo2,ihi1,ihi2,dx,dy,dz,dt,k3d,kc,a_old)

    use amrex_error_module
    use amrex_fort_module, only : rt => amrex_real
    use meth_params_module, only : ppm_type

    implicit none

    integer         , intent(in   ) ::   s_l1, s_l2, s_l3, s_h1, s_h2, s_h3
    integer         , intent(in   ) ::  qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
    integer         , intent(in   ) ::   f_l1, f_l2, f_l3, f_h1, f_h2, f_h3
    integer         , intent(in   ) ::  ilo1,ilo2,ihi1,ihi2
 
    real(rt), intent(in   ) ::      s( s_l1: s_h1, s_l2: s_h2, s_l3: s_h3)
    real(rt), intent(in   ) ::      u(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,3)
    real(rt), intent(in   ) ::   cspd(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
    real(rt), intent(in   ) ::  flatn( f_l1: f_h1, f_l2: f_h2, f_l3: f_h3)

    real(rt), intent(inout) :: Ip(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3)
    real(rt), intent(inout) :: Im(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3)

    real(rt), intent(in   ) ::  dx,dy,dz,dt,a_old

    real(rt) :: dt_over_a
    integer          :: k3d,kc

    dt_over_a = dt / a_old
   
    if (ppm_type .eq. 1) then

        call ppm_type1(s,s_l1,s_l2,s_l3,s_h1,s_h2,s_h3, &
                       u,cspd,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                       flatn,f_l1,f_l2,f_l3,f_h1,f_h2,f_h3, &
                       Ip,Im,ilo1,ilo2,ihi1,ihi2,dx,dy,dz,dt_over_a,k3d,kc)

    else if (ppm_type .eq. 2) then

        call ppm_type2(s,s_l1,s_l2,s_l3,s_h1,s_h2,s_h3, &
                       u,cspd,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                       flatn,f_l1,f_l2,f_l3,f_h1,f_h2,f_h3, &
                       Ip,Im,ilo1,ilo2,ihi1,ihi2,dx,dy,dz,dt_over_a,k3d,kc)

    end if

  end subroutine ppm

  ! :::
  ! ::: ----------------------------------------------------------------
  ! :::
    
  subroutine ppm_type1(s,s_l1,s_l2,s_l3,s_h1,s_h2,s_h3, &
                       u,cspd,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                       flatn,f_l1,f_l2,f_l3,f_h1,f_h2,f_h3, &
                       Ip,Im,ilo1,ilo2,ihi1,ihi2,dx,dy,dz,dt_over_a,k3d,kc)

    use amrex_error_module
    use amrex_mempool_module, only: amrex_allocate, amrex_deallocate
    use amrex_fort_module, only : rt => amrex_real
    use meth_params_module, only : ppm_type, ppm_flatten_before_integrals
    use amrex_constants_module

    implicit none

    integer           s_l1, s_l2, s_l3, s_h1, s_h2, s_h3
    integer          qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
    integer           f_l1, f_l2, f_l3, f_h1, f_h2, f_h3
    integer          ilo1,ilo2,ihi1,ihi2

    real(rt), intent(in) :: s( s_l1: s_h1, s_l2: s_h2, s_l3: s_h3)
    real(rt), intent(in) :: u(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,3)
    real(rt), intent(in) :: cspd(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
    real(rt), intent(in) :: flatn(f_l1: f_h1, f_l2: f_h2, f_l3: f_h3)

    real(rt), intent(out) :: Ip(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3)
    real(rt), intent(out) :: Im(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3)

    ! Note that dt_over_a = dt / a_old
    real(rt), intent(in) :: dx,dy,dz,dt_over_a
    integer, intent(in)          :: k3d,kc

    ! local
    integer i,j,k

    real(rt) :: dxinv,dyinv,dzinv

    real(rt), pointer :: dsl(:), dsr(:), dsc(:)
    real(rt), pointer :: sigma(:), s6(:)

    ! s_{\ib,+}, s_{\ib,-}
    real(rt), pointer :: sp(:)
    real(rt), pointer :: sm(:)

    ! \delta s_{\ib}^{vL}
    real(rt), pointer :: dsvl(:,:)
    real(rt), pointer :: dsvlm(:,:)
    real(rt), pointer :: dsvlp(:,:)

    ! s_{i+\half}^{H.O.}
    real(rt), pointer :: sedge(:,:)
    real(rt), pointer :: sedgez(:,:,:)

    dxinv = 1.0d0/dx
    dyinv = 1.0d0/dy
    dzinv = 1.0d0/dz

    ! cell-centered indexing
    call amrex_allocate(sp,ilo1-1,ihi1+1)
    call amrex_allocate(sm,ilo1-1,ihi1+1)

    call amrex_allocate(sigma,ilo1-1,ihi1+1)
    call amrex_allocate(s6,ilo1-1,ihi1+1)

    if (ppm_type .ne. 1) &
         call amrex_error("Should have ppm_type = 1 in ppm_type1")

    if (s_l1 .gt. ilo1-3 .or. s_l2 .gt. ilo2-3) then
         print *,'Low bounds of array: ',s_l1, s_l2
         print *,'Low bounds of  loop: ',ilo1 , ilo2
         call amrex_error("Need more ghost cells on array in ppm_type1")
    end if

    if (s_h1 .lt. ihi1+3 .or. s_h2 .lt. ihi2+3) then
         print *,'Hi  bounds of array: ',s_h1, s_h2
         print *,'Hi  bounds of  loop: ',ihi1 , ihi2
         call amrex_error("Need more ghost cells on array in ppm_type1")
    end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! x-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! cell-centered indexing w/extra x-ghost cell
    call amrex_allocate(dsvl,ilo1-2,ihi1+2,ilo2-1,ihi2+1)

    ! edge-centered indexing for x-faces -- ppm_type = 1 only
    call amrex_allocate(sedge,ilo1-1,ihi1+2,ilo2-1,ihi2+1)

    ! cell-centered indexing
    call amrex_allocate(dsc,ilo1-2,ihi1+2)
    call amrex_allocate(dsl,ilo1-2,ihi1+2)
    call amrex_allocate(dsr,ilo1-2,ihi1+2)

    ! compute s at x-edges

    ! compute van Leer slopes in x-direction
    dsvl = ZERO
    do j=ilo2-1,ihi2+1
       do i=ilo1-2,ihi1+2
          dsc(i) = HALF * (s(i+1,j,k3d) - s(i-1,j,k3d))
          dsl(i) = TWO  * (s(i  ,j,k3d) - s(i-1,j,k3d))
          dsr(i) = TWO  * (s(i+1,j,k3d) - s(i  ,j,k3d))
          if (dsl(i)*dsr(i) .gt. ZERO) &
               dsvl(i,j) = sign(ONE,dsc(i))*min(abs(dsc(i)),abs(dsl(i)),abs(dsr(i)))
       end do

       ! interpolate s to x-edges
       do i=ilo1-1,ihi1+2
          sedge(i,j) = HALF*(s(i,j,k3d)+s(i-1,j,k3d)) &
               - SIXTH*(dsvl(i,j)-dsvl(i-1,j))
          ! make sure sedge lies in between adjacent cell-centered values
          sedge(i,j) = max(sedge(i,j),min(s(i,j,k3d),s(i-1,j,k3d)))
          sedge(i,j) = min(sedge(i,j),max(s(i,j,k3d),s(i-1,j,k3d)))
       end do

       ! copy sedge into sp and sm
       sp(ilo1-1:ihi1+1) = sedge(ilo1:ihi1+2,j)
       sm(ilo1-1:ihi1+1) = sedge(ilo1-1:ihi1+1,j)

       if (ppm_flatten_before_integrals == 1) then
          ! Flatten the parabola BEFORE doing the other
          ! monotonization -- this is the method that Flash does
          sm = flatn(ilo1-1:ihi1+1,j,k3d)*sm + (ONE-flatn(ilo1-1:ihi1+1,j,k3d))*s(ilo1-1:ihi1+1,j,k3d)
          sp = flatn(ilo1-1:ihi1+1,j,k3d)*sp + (ONE-flatn(ilo1-1:ihi1+1,j,k3d))*s(ilo1-1:ihi1+1,j,k3d)
       endif

       ! Modify using quadratic limiters -- note this version of the limiting comes
       ! from Colella and Sekora (2008), not the original PPM paper.
       do i = ilo1-1, ihi1+1
          if ((sp(i)-s(i,j,k3d))*(s(i,j,k3d)-sm(i)) .le. ZERO) then
             sp(i) = s(i,j,k3d)
             sm(i) = s(i,j,k3d)
          else if (abs(sp(i)-s(i,j,k3d)) .ge. TWO*abs(sm(i)-s(i,j,k3d))) then
             sp(i) = THREE*s(i,j,k3d) - TWO*sm(i)
          else if (abs(sm(i)-s(i,j,k3d)) .ge. TWO*abs(sp(i)-s(i,j,k3d))) then
             sm(i) = THREE*s(i,j,k3d) - TWO*sp(i)
          end if
       end do

       if (ppm_flatten_before_integrals == 2) then
          ! Flatten the parabola AFTER doing the monotonization --
          ! this is the method that Miller & Colella do
          sm = flatn(ilo1-1:ihi1+1,j,k3d)*sm + (ONE-flatn(ilo1-1:ihi1+1,j,k3d))*s(ilo1-1:ihi1+1,j,k3d)
          sp = flatn(ilo1-1:ihi1+1,j,k3d)*sp + (ONE-flatn(ilo1-1:ihi1+1,j,k3d))*s(ilo1-1:ihi1+1,j,k3d)
       endif

       ! compute x-component of Ip and Im
       s6 = SIX*s(ilo1-1:ihi1+1,j,k3d) - THREE*(sm+sp)

       ! Ip/m is the integral under the parabola for the extent
       ! that a wave can travel over a timestep
       !
       ! Ip integrates to the right edge of a cell
       ! Im integrates to the left edge of a cell

       ! u-c wave
       sigma = abs(u(ilo1-1:ihi1+1,j,k3d,1)-cspd(ilo1-1:ihi1+1,j,k3d))*dt_over_a*dxinv

       do i = ilo1-1, ihi1+1
          if (u(i,j,k3d,1)-cspd(i,j,k3d) <= ZERO) then
             Ip(i,j,kc,1,1) = sp(i)
          else
             Ip(i,j,kc,1,1) = sp(i) - &
               HALF*sigma(i)*(sp(i)-sm(i)-(ONE-TWO3RD*sigma(i))*s6(i))
          endif
       end do

       do i = ilo1-1, ihi1+1
          if (u(i,j,k3d,1)-cspd(i,j,k3d) >= ZERO) then
             Im(i,j,kc,1,1) = sm(i)
          else
             Im(i,j,kc,1,1) = sm(i) + &
               HALF*sigma(i)*(sp(i)-sm(i)+(ONE-TWO3RD*sigma(i))*s6(i))
          endif
       end do

       ! u wave
       sigma = abs(u(ilo1-1:ihi1+1,j,k3d,1))*dt_over_a*dxinv

       do i = ilo1-1, ihi1+1
          if (u(i,j,k3d,1) <= ZERO) then
             Ip(i,j,kc,1,2) = sp(i)
          else
             Ip(i,j,kc,1,2) = sp(i) - &
               HALF*sigma(i)*(sp(i)-sm(i)-(ONE-TWO3RD*sigma(i))*s6(i))
          endif
       end do

       do i = ilo1-1, ihi1+1
          if (u(i,j,k3d,1) >= ZERO) then
             Im(i,j,kc,1,2) = sm(i)
          else
             Im(i,j,kc,1,2) = sm(i) + &
               HALF*sigma(i)*(sp(i)-sm(i)+(ONE-TWO3RD*sigma(i))*s6(i))
          endif
       end do

       ! u+c wave
       sigma = abs(u(ilo1-1:ihi1+1,j,k3d,1)+cspd(ilo1-1:ihi1+1,j,k3d))*dt_over_a*dxinv

       do i = ilo1-1, ihi1+1
          if (u(i,j,k3d,1)+cspd(i,j,k3d) <= ZERO) then
             Ip(i,j,kc,1,3) = sp(i)
          else
             Ip(i,j,kc,1,3) = sp(i) - &
               HALF*sigma(i)*(sp(i)-sm(i)-(ONE-TWO3RD*sigma(i))*s6(i))
          endif
       end do

       do i = ilo1-1, ihi1+1
          if (u(i,j,k3d,1)+cspd(i,j,k3d) >= ZERO) then
             Im(i,j,kc,1,3) = sm(i)
          else
             Im(i,j,kc,1,3) = sm(i) + &
               HALF*sigma(i)*(sp(i)-sm(i)+(ONE-TWO3RD*sigma(i))*s6(i))
          endif
       end do

    end do

    call amrex_deallocate(dsc)
    call amrex_deallocate(dsl)
    call amrex_deallocate(dsr)
    call amrex_deallocate(sedge)
    call amrex_deallocate(dsvl)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! y-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! cell-centered indexing w/extra y-ghost cell
    call amrex_allocate( dsvl,ilo1-1,ihi1+1,ilo2-2,ihi2+2)

    ! edge-centered indexing for y-faces
    call amrex_allocate(sedge,ilo1-1,ihi1+1,ilo2-1,ihi2+2)

    ! cell-centered indexing
    call amrex_allocate(dsc,ilo1-1,ihi1+1)
    call amrex_allocate(dsl,ilo1-1,ihi1+1)
    call amrex_allocate(dsr,ilo1-1,ihi1+1)

    ! compute s at y-edges

    ! compute van Leer slopes in y-direction
    dsvl = ZERO
    do j=ilo2-2,ihi2+2
       do i=ilo1-1,ihi1+1
          dsc(i) = HALF * (s(i,j+1,k3d) - s(i,j-1,k3d))
          dsl(i) = TWO  * (s(i,j  ,k3d) - s(i,j-1,k3d))
          dsr(i) = TWO  * (s(i,j+1,k3d) - s(i,j  ,k3d))
          if (dsl(i)*dsr(i) .gt. ZERO) &
               dsvl(i,j) = sign(ONE,dsc(i))*min(abs(dsc(i)),abs(dsl(i)),abs(dsr(i)))
       end do
    end do

    ! interpolate s to y-edges
    do j=ilo2-1,ihi2+2
       do i=ilo1-1,ihi1+1
          sedge(i,j) = HALF*(s(i,j,k3d)+s(i,j-1,k3d)) &
               - SIXTH*(dsvl(i,j)-dsvl(i,j-1))
          ! make sure sedge lies in between adjacent cell-centered values
          sedge(i,j) = max(sedge(i,j),min(s(i,j,k3d),s(i,j-1,k3d)))
          sedge(i,j) = min(sedge(i,j),max(s(i,j,k3d),s(i,j-1,k3d)))
       end do
    end do

    do j=ilo2-1,ihi2+1

       ! copy sedge into sp and sm
       sp = sedge(ilo1-1:ihi1+1,j+1)
       sm = sedge(ilo1-1:ihi1+1,j  )

       if (ppm_flatten_before_integrals == 1) then
          ! Flatten the parabola BEFORE doing the other
          ! monotonization -- this is the method that Flash does
          sm = flatn(ilo1-1:ihi1+1,j,k3d)*sm + (ONE-flatn(ilo1-1:ihi1+1,j,k3d))*s(ilo1-1:ihi1+1,j,k3d)
          sp = flatn(ilo1-1:ihi1+1,j,k3d)*sp + (ONE-flatn(ilo1-1:ihi1+1,j,k3d))*s(ilo1-1:ihi1+1,j,k3d)
       endif

       ! Modify using quadratic limiters
       do i = ilo1-1, ihi1+1
          if ((sp(i)-s(i,j,k3d))*(s(i,j,k3d)-sm(i)) .le. ZERO) then
             sp(i) = s(i,j,k3d)
             sm(i) = s(i,j,k3d)
          else if (abs(sp(i)-s(i,j,k3d)) .ge. TWO*abs(sm(i)-s(i,j,k3d))) then
             sp(i) = THREE*s(i,j,k3d) - TWO*sm(i)
          else if (abs(sm(i)-s(i,j,k3d)) .ge. TWO*abs(sp(i)-s(i,j,k3d))) then
             sm(i) = THREE*s(i,j,k3d) - TWO*sp(i)
          end if
       end do

       if (ppm_flatten_before_integrals == 2) then
          ! Flatten the parabola AFTER doing the monotonization --
          ! this is the method that Miller & Colella do
          sm = flatn(ilo1-1:ihi1+1,j,k3d)*sm + (ONE-flatn(ilo1-1:ihi1+1,j,k3d))*s(ilo1-1:ihi1+1,j,k3d)
          sp = flatn(ilo1-1:ihi1+1,j,k3d)*sp + (ONE-flatn(ilo1-1:ihi1+1,j,k3d))*s(ilo1-1:ihi1+1,j,k3d)
       endif

       ! compute y-component of Ip and Im
       s6 = SIX*s(ilo1-1:ihi1+1,j,k3d) - THREE*(sm+sp)

       ! v-c wave
       sigma = abs(u(ilo1-1:ihi1+1,j,k3d,2)-cspd(ilo1-1:ihi1+1,j,k3d))*dt_over_a*dyinv

       do i = ilo1-1, ihi1+1
          if (u(i,j,k3d,2)-cspd(i,j,k3d) <= ZERO) then
             Ip(i,j,kc,2,1) = sp(i)
          else
             Ip(i,j,kc,2,1) = sp(i) - &
               HALF*sigma(i)*(sp(i)-sm(i)-(ONE-TWO3RD*sigma(i))*s6(i))
          endif
       end do

       do i = ilo1-1, ihi1+1
          if (u(i,j,k3d,2)-cspd(i,j,k3d) >= ZERO) then
             Im(i,j,kc,2,1) = sm(i)
          else
             Im(i,j,kc,2,1) = sm(i) + &
               HALF*sigma(i)*(sp(i)-sm(i)+(ONE-TWO3RD*sigma(i))*s6(i))
          endif
       end do

       ! v wave
       sigma = abs(u(ilo1-1:ihi1+1,j,k3d,2))*dt_over_a*dyinv

       do i = ilo1-1, ihi1+1
          if (u(i,j,k3d,2) <= ZERO) then
             Ip(i,j,kc,2,2) = sp(i)
          else
             Ip(i,j,kc,2,2) = sp(i) - &
               HALF*sigma(i)*(sp(i)-sm(i)-(ONE-TWO3RD*sigma(i))*s6(i))
          endif
       end do

       do i = ilo1-1, ihi1+1
          if (u(i,j,k3d,2) >= ZERO) then
             Im(i,j,kc,2,2) = sm(i)
          else
             Im(i,j,kc,2,2) = sm(i) + &
               HALF*sigma(i)*(sp(i)-sm(i)+(ONE-TWO3RD*sigma(i))*s6(i))
          endif
       end do

       ! v+c wave
       sigma = abs(u(ilo1-1:ihi1+1,j,k3d,2)+cspd(ilo1-1:ihi1+1,j,k3d))*dt_over_a*dyinv

       do i = ilo1-1, ihi1+1
          if (u(i,j,k3d,2)+cspd(i,j,k3d) <= ZERO) then
             Ip(i,j,kc,2,3) = sp(i)
          else
             Ip(i,j,kc,2,3) = sp(i) - &
               HALF*sigma(i)*(sp(i)-sm(i)-(ONE-TWO3RD*sigma(i))*s6(i))
          endif
       end do

       do i = ilo1-1, ihi1+1
          if (u(i,j,k3d,2)+cspd(i,j,k3d) >= ZERO) then
             Im(i,j,kc,2,3) = sm(i)
          else
             Im(i,j,kc,2,3) = sm(i) + &
               HALF*sigma(i)*(sp(i)-sm(i)+(ONE-TWO3RD*sigma(i))*s6(i))
          endif
       end do

    end do

    call amrex_deallocate(dsc)
    call amrex_deallocate(dsl)
    call amrex_deallocate(dsr)
    call amrex_deallocate(dsvl)
    call amrex_deallocate(sedge)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! z-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! cell-centered indexing
    call amrex_allocate( dsvl,ilo1-1,ihi1+1,ilo2-1,ihi2+1)
    call amrex_allocate(dsvlm,ilo1-1,ihi1+1,ilo2-1,ihi2+1)
    call amrex_allocate(dsvlp,ilo1-1,ihi1+1,ilo2-1,ihi2+1)

    ! cell-centered indexing
    call amrex_allocate(dsc,ilo1-1,ihi1+1)
    call amrex_allocate(dsl,ilo1-1,ihi1+1)
    call amrex_allocate(dsr,ilo1-1,ihi1+1)

    call amrex_allocate(sedgez,ilo1-1,ihi1+1,ilo2-2,ihi2+3,k3d-1,k3d+2)

    ! compute s at z-edges

    ! compute van Leer slopes in z-direction
    dsvl  = ZERO
    dsvlm = ZERO
    dsvlp = ZERO

    do j=ilo2-1,ihi2+1

       ! compute on slab below
       k = k3d-1
       dsc = HALF * (s(ilo1-1:ihi1+1,j,k+1) - s(ilo1-1:ihi1+1,j,k-1))
       dsl = TWO  * (s(ilo1-1:ihi1+1,j,k  ) - s(ilo1-1:ihi1+1,j,k-1))
       dsr = TWO  * (s(ilo1-1:ihi1+1,j,k+1) - s(ilo1-1:ihi1+1,j,k  ))
       do i = ilo1-1, ihi1+1
          if (dsl(i)*dsr(i) .gt. ZERO) &
               dsvlm(i,j) = sign(ONE,dsc(i))*min(abs(dsc(i)),abs(dsl(i)),abs(dsr(i)))
       end do

       ! compute on slab above
       k = k3d+1
       dsc = HALF * (s(ilo1-1:ihi1+1,j,k+1) - s(ilo1-1:ihi1+1,j,k-1))
       dsl = TWO  * (s(ilo1-1:ihi1+1,j,k  ) - s(ilo1-1:ihi1+1,j,k-1))
       dsr = TWO  * (s(ilo1-1:ihi1+1,j,k+1) - s(ilo1-1:ihi1+1,j,k  ))
       do i = ilo1-1, ihi1+1
          if (dsl(i)*dsr(i) .gt. ZERO) &
               dsvlp(i,j) = sign(ONE,dsc(i))*min(abs(dsc(i)),abs(dsl(i)),abs(dsr(i)))
       end do

       ! compute on current slab
       k = k3d
       dsc = HALF * (s(ilo1-1:ihi1+1,j,k+1) - s(ilo1-1:ihi1+1,j,k-1))
       dsl = TWO  * (s(ilo1-1:ihi1+1,j,k  ) - s(ilo1-1:ihi1+1,j,k-1))
       dsr = TWO  * (s(ilo1-1:ihi1+1,j,k+1) - s(ilo1-1:ihi1+1,j,k  ))
       do i = ilo1-1, ihi1+1
          if (dsl(i)*dsr(i) .gt. ZERO) &
               dsvl(i,j) = sign(ONE,dsc(i))*min(abs(dsc(i)),abs(dsl(i)),abs(dsr(i)))
       end do

       ! interpolate to lo face
       k = k3d
       sm = HALF*(s(ilo1-1:ihi1+1,j,k)+s(ilo1-1:ihi1+1,j,k-1)) - SIXTH*(dsvl(ilo1-1:ihi1+1,j)-dsvlm(ilo1-1:ihi1+1,j))
       ! make sure sedge lies in between adjacent cell-centered values
       sm = max(sm,min(s(ilo1-1:ihi1+1,j,k),s(ilo1-1:ihi1+1,j,k-1)))
       sm = min(sm,max(s(ilo1-1:ihi1+1,j,k),s(ilo1-1:ihi1+1,j,k-1)))

       ! interpolate to hi face
       k = k3d+1
       sp = HALF*(s(ilo1-1:ihi1+1,j,k)+s(ilo1-1:ihi1+1,j,k-1)) - SIXTH*(dsvlp(ilo1-1:ihi1+1,j)-dsvl(ilo1-1:ihi1+1,j))
       ! make sure sedge lies in between adjacent cell-centered values
       sp = max(sp,min(s(ilo1-1:ihi1+1,j,k),s(ilo1-1:ihi1+1,j,k-1)))
       sp = min(sp,max(s(ilo1-1:ihi1+1,j,k),s(ilo1-1:ihi1+1,j,k-1)))

       if (ppm_flatten_before_integrals == 1) then
          ! flatten the parabola BEFORE doing the other
          ! monotonization -- this is the method that Flash does
          sm = flatn(ilo1-1:ihi1+1,j,k3d)*sm + (ONE-flatn(ilo1-1:ihi1+1,j,k3d))*s(ilo1-1:ihi1+1,j,k3d)
          sp = flatn(ilo1-1:ihi1+1,j,k3d)*sp + (ONE-flatn(ilo1-1:ihi1+1,j,k3d))*s(ilo1-1:ihi1+1,j,k3d)
       endif

       ! modify using quadratic limiters
       do i = ilo1-1, ihi1+1
          if ((sp(i)-s(i,j,k3d))*(s(i,j,k3d)-sm(i)) .le. ZERO) then
             sp(i) = s(i,j,k3d)
             sm(i) = s(i,j,k3d)
          else if (abs(sp(i)-s(i,j,k3d)) .ge. TWO*abs(sm(i)-s(i,j,k3d))) then
             sp(i) = THREE*s(i,j,k3d) - TWO*sm(i)
          else if (abs(sm(i)-s(i,j,k3d)) .ge. TWO*abs(sp(i)-s(i,j,k3d))) then
             sm(i) = THREE*s(i,j,k3d) - TWO*sp(i)
          end if
       end do

       if (ppm_flatten_before_integrals == 2) then
          ! flatten the parabola AFTER doing the monotonization --
          ! this is the method that Miller & Colella do
          sm = flatn(ilo1-1:ihi1+1,j,k3d)*sm + (ONE-flatn(ilo1-1:ihi1+1,j,k3d))*s(ilo1-1:ihi1+1,j,k3d)
          sp = flatn(ilo1-1:ihi1+1,j,k3d)*sp + (ONE-flatn(ilo1-1:ihi1+1,j,k3d))*s(ilo1-1:ihi1+1,j,k3d)
       endif

       ! compute z-component of Ip and Im
       s6 = SIX*s(ilo1-1:ihi1+1,j,k3d) - THREE*(sm+sp)

       ! w-c wave
       sigma = abs(u(ilo1-1:ihi1+1,j,k3d,3)-cspd(ilo1-1:ihi1+1,j,k3d))*dt_over_a*dzinv

       do i = ilo1-1, ihi1+1
          if (u(i,j,k3d,3)-cspd(i,j,k3d) <= ZERO) then
             Ip(i,j,kc,3,1) = sp(i)
          else
             Ip(i,j,kc,3,1) = sp(i) - &
               HALF*sigma(i)*(sp(i)-sm(i)-(ONE-TWO3RD*sigma(i))*s6(i))
          endif
       end do

       do i = ilo1-1, ihi1+1
          if (u(i,j,k3d,3)-cspd(i,j,k3d) >= ZERO) then
             Im(i,j,kc,3,1) = sm(i)
          else
             Im(i,j,kc,3,1) = sm(i) + &
               HALF*sigma(i)*(sp(i)-sm(i)+(ONE-TWO3RD*sigma(i))*s6(i))
          endif
       end do

       ! w wave
       sigma = abs(u(ilo1-1:ihi1+1,j,k3d,3))*dt_over_a*dzinv

       do i = ilo1-1, ihi1+1
          if (u(i,j,k3d,3) <= ZERO) then
             Ip(i,j,kc,3,2) = sp(i)
          else
             Ip(i,j,kc,3,2) = sp(i) - &
               HALF*sigma(i)*(sp(i)-sm(i)-(ONE-TWO3RD*sigma(i))*s6(i))
          endif
       end do

       do i = ilo1-1, ihi1+1
          if (u(i,j,k3d,3) >= ZERO) then
             Im(i,j,kc,3,2) = sm(i)
          else
             Im(i,j,kc,3,2) = sm(i) + &
               HALF*sigma(i)*(sp(i)-sm(i)+(ONE-TWO3RD*sigma(i))*s6(i))
          endif
       end do

       ! w+c wave
       sigma = abs(u(ilo1-1:ihi1+1,j,k3d,3)+cspd(ilo1-1:ihi1+1,j,k3d))*dt_over_a*dzinv

       do i = ilo1-1, ihi1+1
          if (u(i,j,k3d,3)+cspd(i,j,k3d) <= ZERO) then
             Ip(i,j,kc,3,3) = sp(i)
          else
             Ip(i,j,kc,3,3) = sp(i) - &
               HALF*sigma(i)*(sp(i)-sm(i)-(ONE-TWO3RD*sigma(i))*s6(i))
          endif
       end do

       do i = ilo1-1, ihi1+1
          if (u(i,j,k3d,3)+cspd(i,j,k3d) >= ZERO) then
             Im(i,j,kc,3,3) = sm(i)
          else
             Im(i,j,kc,3,3) = sm(i) + &
               HALF*sigma(i)*(sp(i)-sm(i)+(ONE-TWO3RD*sigma(i))*s6(i))
          endif
       end do

    end do

    call amrex_deallocate(dsc)
    call amrex_deallocate(dsl)
    call amrex_deallocate(dsr)
    call amrex_deallocate(dsvl)
    call amrex_deallocate(dsvlm)
    call amrex_deallocate(dsvlp)
    call amrex_deallocate(sp)
    call amrex_deallocate(sm)
    call amrex_deallocate(sedgez)
    call amrex_deallocate(sigma)
    call amrex_deallocate(s6)

  end subroutine ppm_type1

  ! :::
  ! ::: ----------------------------------------------------------------
  ! :::

  subroutine ppm_type2(s,s_l1,s_l2,s_l3,s_h1,s_h2,s_h3, &
                       u,cspd,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                       flatn,f_l1,f_l2,f_l3,f_h1,f_h2,f_h3, &
                       Ip,Im,ilo1,ilo2,ihi1,ihi2,dx,dy,dz,dt_over_a,k3d,kc)

    use amrex_error_module
    use amrex_fort_module, only : rt => amrex_real
    use meth_params_module, only : ppm_type, ppm_flatten_before_integrals
    use amrex_constants_module

    implicit none

    integer           s_l1, s_l2, s_l3, s_h1, s_h2, s_h3
    integer          qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
    integer           f_l1, f_l2, f_l3, f_h1, f_h2, f_h3
    integer          ilo1,ilo2,ihi1,ihi2

    real(rt)    s( s_l1: s_h1, s_l2: s_h2, s_l3: s_h3)
    real(rt)    u(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,3)
    real(rt) cspd(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
    real(rt) flatn( f_l1: f_h1, f_l2: f_h2, f_l3: f_h3)

    real(rt) Ip(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3)
    real(rt) Im(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3)

    ! Note that dt_over_a = dt / a_old
    real(rt) dx,dy,dz,dt_over_a
    real(rt) dxinv,dyinv,dzinv
    integer          k3d,kc

    ! local
    integer i,j,k
    logical extremum, bigp, bigm

    real(rt) D2, D2C, D2L, D2R, D2LIM, alphap, alpham
    real(rt) sgn, sigma, s6
    real(rt) dafacem, dafacep, dabarm, dabarp, dafacemin, dabarmin
    real(rt) dachkm, dachkp
    real(rt) amax, delam, delap

    ! s_{\ib,+}, s_{\ib,-}
    real(rt), allocatable :: sp(:,:)
    real(rt), allocatable :: sm(:,:)

    ! \delta s_{\ib}^{vL}
    real(rt), allocatable :: dsvl(:,:)
    real(rt), allocatable :: dsvlm(:,:)
    real(rt), allocatable :: dsvlp(:,:)

    ! s_{i+\half}^{H.O.}
    real(rt), allocatable :: sedge(:,:)
    real(rt), allocatable :: sedgez(:,:,:)

    ! constant used in Colella 2008
    real(rt), parameter :: C = 1.25d0

    dxinv = 1.0d0/dx
    dyinv = 1.0d0/dy
    dzinv = 1.0d0/dz

    ! cell-centered indexing
    allocate(sp(ilo1-1:ihi1+1,ilo2-1:ihi2+1))
    allocate(sm(ilo1-1:ihi1+1,ilo2-1:ihi2+1))

    if (ppm_type .ne. 2) &
         call amrex_error("Should have ppm_type = 2 in ppm_type2")

    if (s_l1 .gt. ilo1-3 .or. s_l2 .gt. ilo2-3) then
         print *,'Low bounds of array: ',s_l1, s_l2
         print *,'Low bounds of  loop: ',ilo1 , ilo2
         call amrex_error("Need more ghost cells on array in ppm_type2")
    end if

    if (s_h1 .lt. ihi1+3 .or. s_h2 .lt. ihi2+3) then
         print *,'Hi  bounds of array: ',s_h1, s_h2
         print *,'Hi  bounds of  loop: ',ihi1 , ihi2
         call amrex_error("Need more ghost cells on array in ppm_type2")
    end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! x-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! cell-centered indexing w/extra x-ghost cell
    allocate(dsvl(ilo1-2:ihi1+2,ilo2-1:ihi2+1))

    ! edge-centered indexing for x-faces
    allocate(sedge(ilo1-2:ihi1+3,ilo2-1:ihi2+1))

    ! compute s at x-edges

    ! interpolate s to x-edges
    do j=ilo2-1,ihi2+1
       do i=ilo1-2,ihi1+3
          sedge(i,j) = SEVEN12TH*(s(i-1,j,k3d)+s(i  ,j,k3d)) &
               - TWELFTH*(s(i-2,j,k3d)+s(i+1,j,k3d))
          !
          ! limit sedge
          !
          if ((sedge(i,j)-s(i-1,j,k3d))*(s(i,j,k3d)-sedge(i,j)) .lt. ZERO) then
             D2  = THREE*(s(i-1,j,k3d)-TWO*sedge(i,j)+s(i,j,k3d))
             D2L = s(i-2,j,k3d)-TWO*s(i-1,j,k3d)+s(i,j,k3d)
             D2R = s(i-1,j,k3d)-TWO*s(i,j,k3d)+s(i+1,j,k3d)
             sgn = sign(ONE,D2)
             D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),ZERO)
             sedge(i,j) = HALF*(s(i-1,j,k3d)+s(i,j,k3d)) - SIXTH*D2LIM
          end if
       end do
    end do
    !
    ! Use Colella 2008 limiters.
    !
    ! This is a new version of the algorithm to eliminate sensitivity to roundoff.
    !
    do j=ilo2-1,ihi2+1
       do i=ilo1-1,ihi1+1

          alphap   = sedge(i+1,j)-s(i,j,k3d)
          alpham   = sedge(i  ,j)-s(i,j,k3d)
          bigp     = abs(alphap).gt.TWO*abs(alpham)
          bigm     = abs(alpham).gt.TWO*abs(alphap)
          extremum = .false.

          if (alpham*alphap .ge. ZERO) then
             extremum = .true.
          else if (bigp .or. bigm) then
             !
             ! Possible extremum. We look at cell centered values and face
             ! centered values for a change in sign in the differences adjacent to
             ! the cell. We use the pair of differences whose minimum magnitude is the
             ! largest, and thus least susceptible to sensitivity to roundoff.
             !
             dafacem   = sedge(i,j) - sedge(i-1,j)
             dafacep   = sedge(i+2,j) - sedge(i+1,j)
             dabarm    = s(i,j,k3d) - s(i-1,j,k3d)
             dabarp    = s(i+1,j,k3d) - s(i,j,k3d)
             dafacemin = min(abs(dafacem),abs(dafacep))
             dabarmin  = min(abs(dabarm),abs(dabarp))
             if (dafacemin.ge.dabarmin) then
                dachkm = dafacem
                dachkp = dafacep
             else
                dachkm = dabarm
                dachkp = dabarp
             endif
             extremum = (dachkm*dachkp .le. ZERO)
          end if

          if (extremum) then
             D2     = SIX*(alpham + alphap)
             D2L    = s(i-2,j,k3d)-TWO*s(i-1,j,k3d)+s(i,j,k3d)
             D2R    = s(i,j,k3d)-TWO*s(i+1,j,k3d)+s(i+2,j,k3d)
             D2C    = s(i-1,j,k3d)-TWO*s(i,j,k3d)+s(i+1,j,k3d)
             sgn    = sign(ONE,D2)
             D2LIM  = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
             alpham = alpham*D2LIM/max(abs(D2),1.d-10)
             alphap = alphap*D2LIM/max(abs(D2),1.d-10)
          else
             if (bigp) then
                sgn   = sign(ONE,alpham)
                amax  = -alphap**2 / (4*(alpham + alphap))
                delam = s(i-1,j,k3d) - s(i,j,k3d)
                if (sgn*amax .ge. sgn*delam) then
                   if (sgn*(delam - alpham).ge.1.d-10) then
                      alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                   else 
                      alphap = -TWO*alpham
                   endif
                endif
             end if
             if (bigm) then
                sgn   = sign(ONE,alphap)
                amax  = -alpham**2 / (4*(alpham + alphap))
                delap = s(i+1,j,k3d) - s(i,j,k3d)
                if (sgn*amax .ge. sgn*delap) then
                   if (sgn*(delap - alphap).ge.1.d-10) then
                      alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                   else
                      alpham = -TWO*alphap
                   endif
                endif
             end if
          end if

          sm(i,j) = s(i,j,k3d) + alpham
          sp(i,j) = s(i,j,k3d) + alphap

          if (ppm_flatten_before_integrals > 0) then
             ! flatten the parabola AFTER doing the monotonization
             sm(i,j) = flatn(i,j,k3d)*sm(i,j) + (ONE-flatn(i,j,k3d))*s(i,j,k3d)
             sp(i,j) = flatn(i,j,k3d)*sp(i,j) + (ONE-flatn(i,j,k3d))*s(i,j,k3d)
          endif

          !
          ! Compute x-component of Ip and Im.
          !
          s6    = SIX*s(i,j,k3d) - THREE*(sm(i,j)+sp(i,j))

          ! u-c wave
          sigma = abs(u(i,j,k3d,1)-cspd(i,j,k3d))*dt_over_a*dxinv

          if (u(i,j,k3d,1)-cspd(i,j,k3d) <= ZERO) then
             Ip(i,j,kc,1,1) = sp(i,j)
          else
             Ip(i,j,kc,1,1) = sp(i,j) - &
                  HALF*sigma*(sp(i,j)-sm(i,j)-(ONE-TWO3RD*sigma)*s6)
          endif

          if (u(i,j,k3d,1)-cspd(i,j,k3d) >= ZERO) then
             Im(i,j,kc,1,1) = sm(i,j)
          else
             Im(i,j,kc,1,1) = sm(i,j) + &
                  HALF*sigma*(sp(i,j)-sm(i,j)+(ONE-TWO3RD*sigma)*s6)
          endif

          ! u wave
          sigma = abs(u(i,j,k3d,1))*dt_over_a*dxinv

          if (u(i,j,k3d,1) <= ZERO) then
             Ip(i,j,kc,1,2) = sp(i,j)
          else
             Ip(i,j,kc,1,2) = sp(i,j) - &
                  HALF*sigma*(sp(i,j)-sm(i,j)-(ONE-TWO3RD*sigma)*s6)
          endif

          if (u(i,j,k3d,1) >= ZERO) then
             Im(i,j,kc,1,2) = sm(i,j)
          else
             Im(i,j,kc,1,2) = sm(i,j) + &
                  HALF*sigma*(sp(i,j)-sm(i,j)+(ONE-TWO3RD*sigma)*s6)
          endif

          ! u+c wave
          sigma = abs(u(i,j,k3d,1)+cspd(i,j,k3d))*dt_over_a*dxinv

          if (u(i,j,k3d,1)+cspd(i,j,k3d) <= ZERO) then
             Ip(i,j,kc,1,3) = sp(i,j) 
          else
             Ip(i,j,kc,1,3) = sp(i,j) - &
                  HALF*sigma*(sp(i,j)-sm(i,j)-(ONE-TWO3RD*sigma)*s6)
          endif

          if (u(i,j,k3d,1)+cspd(i,j,k3d) >= ZERO) then
             Im(i,j,kc,1,3) = sm(i,j) 
          else
             Im(i,j,kc,1,3) = sm(i,j) + &
                  HALF*sigma*(sp(i,j)-sm(i,j)+(ONE-TWO3RD*sigma)*s6)
          endif

       end do
    end do

    deallocate(sedge,dsvl)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! y-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! cell-centered indexing w/extra y-ghost cell
    allocate( dsvl(ilo1-1:ihi1+1,ilo2-2:ihi2+2))

    ! edge-centered indexing for y-faces
    allocate(sedge(ilo1-1:ihi1+1,ilo2-2:ihi2+3))

    ! compute s at y-edges

    ! interpolate s to y-edges
    do j=ilo2-2,ihi2+3
       do i=ilo1-1,ihi1+1
          sedge(i,j) = SEVEN12TH*(s(i,j-1,k3d)+s(i,j,k3d)) &
               - TWELFTH*(s(i,j-2,k3d)+s(i,j+1,k3d))
          !
          ! limit sedge
          !
          if ((sedge(i,j)-s(i,j-1,k3d))*(s(i,j,k3d)-sedge(i,j)) .lt. ZERO) then
             D2  = THREE*(s(i,j-1,k3d)-TWO*sedge(i,j)+s(i,j,k3d))
             D2L = s(i,j-2,k3d)-TWO*s(i,j-1,k3d)+s(i,j,k3d)
             D2R = s(i,j-1,k3d)-TWO*s(i,j,k3d)+s(i,j+1,k3d)
             sgn = sign(ONE,D2)
             D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),ZERO)
             sedge(i,j) = HALF*(s(i,j-1,k3d)+s(i,j,k3d)) - SIXTH*D2LIM
          end if
       end do
    end do
    !
    ! Use Colella 2008 limiters.
    !
    ! This is a new version of the algorithm to eliminate sensitivity to roundoff.
    !
    do j=ilo2-1,ihi2+1
       do i=ilo1-1,ihi1+1

          alphap   = sedge(i,j+1)-s(i,j,k3d)
          alpham   = sedge(i,j  )-s(i,j,k3d)
          bigp     = abs(alphap).gt.TWO*abs(alpham)
          bigm     = abs(alpham).gt.TWO*abs(alphap)
          extremum = .false.

          if (alpham*alphap .ge. ZERO) then
             extremum = .true.
          else if (bigp .or. bigm) then
             !
             ! Possible extremum. We look at cell centered values and face
             ! centered values for a change in sign in the differences adjacent to
             ! the cell. We use the pair of differences whose minimum magnitude is the
             ! largest, and thus least susceptible to sensitivity to roundoff.
             !
             dafacem   = sedge(i,j) - sedge(i,j-1)
             dafacep   = sedge(i,j+2) - sedge(i,j+1)
             dabarm    = s(i,j,k3d) - s(i,j-1,k3d)
             dabarp    = s(i,j+1,k3d) - s(i,j,k3d)
             dafacemin = min(abs(dafacem),abs(dafacep))
             dabarmin  = min(abs(dabarm),abs(dabarp))
             if (dafacemin.ge.dabarmin) then
                dachkm = dafacem
                dachkp = dafacep
             else
                dachkm = dabarm
                dachkp = dabarp
             endif
             extremum = (dachkm*dachkp .le. ZERO)
          end if

          if (extremum) then
             D2     = SIX*(alpham + alphap)
             D2L    = s(i,j-2,k3d)-TWO*s(i,j-1,k3d)+s(i,j,k3d)
             D2R    = s(i,j,k3d)-TWO*s(i,j+1,k3d)+s(i,j+2,k3d)
             D2C    = s(i,j-1,k3d)-TWO*s(i,j,k3d)+s(i,j+1,k3d)
             sgn    = sign(ONE,D2)
             D2LIM  = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
             alpham = alpham*D2LIM/max(abs(D2),1.d-10)
             alphap = alphap*D2LIM/max(abs(D2),1.d-10)
          else
             if (bigp) then
                sgn   = sign(ONE,alpham)
                amax  = -alphap**2 / (4*(alpham + alphap))
                delam = s(i,j-1,k3d) - s(i,j,k3d)
                if (sgn*amax .ge. sgn*delam) then
                   if (sgn*(delam - alpham).ge.1.d-10) then
                      alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                   else 
                      alphap = -TWO*alpham
                   endif
                endif
             end if
             if (bigm) then
                sgn   = sign(ONE,alphap)
                amax  = -alpham**2 / (4*(alpham + alphap))
                delap = s(i,j+1,k3d) - s(i,j,k3d)
                if (sgn*amax .ge. sgn*delap) then
                   if (sgn*(delap - alphap).ge.1.d-10) then
                      alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                   else
                      alpham = -TWO*alphap
                   endif
                endif
             end if
          end if

          sm(i,j) = s(i,j,k3d) + alpham
          sp(i,j) = s(i,j,k3d) + alphap

          if (ppm_flatten_before_integrals > 0) then
             ! flatten the parabola AFTER doing the monotonization
             sm(i,j) = flatn(i,j,k3d)*sm(i,j) + (ONE-flatn(i,j,k3d))*s(i,j,k3d)
             sp(i,j) = flatn(i,j,k3d)*sp(i,j) + (ONE-flatn(i,j,k3d))*s(i,j,k3d)
          endif

          !
          ! Compute y-component of Ip and Im.
          !
          s6    = SIX*s(i,j,k3d) - THREE*(sm(i,j)+sp(i,j))

          ! v-c wave
          sigma = abs(u(i,j,k3d,2)-cspd(i,j,k3d))*dt_over_a*dyinv

          if (u(i,j,k3d,2)-cspd(i,j,k3d) <= ZERO) then
             Ip(i,j,kc,2,1) = sp(i,j) 
          else
             Ip(i,j,kc,2,1) = sp(i,j) - &
                  HALF*sigma*(sp(i,j)-sm(i,j)-(ONE-TWO3RD*sigma)*s6)
          endif

          if (u(i,j,k3d,2)-cspd(i,j,k3d) >= ZERO) then
             Im(i,j,kc,2,1) = sm(i,j) 
          else
             Im(i,j,kc,2,1) = sm(i,j) + &
                  HALF*sigma*(sp(i,j)-sm(i,j)+(ONE-TWO3RD*sigma)*s6)
          endif

          ! v wave
          sigma = abs(u(i,j,k3d,2))*dt_over_a*dyinv

          if (u(i,j,k3d,2) <= ZERO) then
             Ip(i,j,kc,2,2) = sp(i,j) 
          else
             Ip(i,j,kc,2,2) = sp(i,j) - &
                  HALF*sigma*(sp(i,j)-sm(i,j)-(ONE-TWO3RD*sigma)*s6)
          endif

          if (u(i,j,k3d,2) >= ZERO) then
             Im(i,j,kc,2,2) = sm(i,j) 
          else
             Im(i,j,kc,2,2) = sm(i,j) + &
                  HALF*sigma*(sp(i,j)-sm(i,j)+(ONE-TWO3RD*sigma)*s6)
          endif

          ! v+c wave
          sigma = abs(u(i,j,k3d,2)+cspd(i,j,k3d))*dt_over_a*dyinv

          if (u(i,j,k3d,2)+cspd(i,j,k3d) <= ZERO) then
             Ip(i,j,kc,2,3) = sp(i,j) 
          else
             Ip(i,j,kc,2,3) = sp(i,j) - &
                  HALF*sigma*(sp(i,j)-sm(i,j)-(ONE-TWO3RD*sigma)*s6)
          endif

          if (u(i,j,k3d,2)+cspd(i,j,k3d) >= ZERO) then
             Im(i,j,kc,2,3) = sm(i,j) 
          else
             Im(i,j,kc,2,3) = sm(i,j) + &
                  HALF*sigma*(sp(i,j)-sm(i,j)+(ONE-TWO3RD*sigma)*s6)
          endif
       end do
    end do

    deallocate(dsvl,sedge)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! z-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! cell-centered indexing
    allocate( dsvl(ilo1-1:ihi1+1,ilo2-1:ihi2+1))
    allocate(dsvlm(ilo1-1:ihi1+1,ilo2-1:ihi2+1))
    allocate(dsvlp(ilo1-1:ihi1+1,ilo2-1:ihi2+1))

    allocate(sedgez(ilo1-1:ihi1+1,ilo2-2:ihi2+3,k3d-1:k3d+2))

    ! compute s at z-edges

    ! interpolate s to z-edges
    do k=k3d-1,k3d+2
       do j=ilo2-1,ihi2+1
          do i=ilo1-1,ihi1+1
             sedgez(i,j,k) = SEVEN12TH*(s(i,j,k-1)+s(i,j,k)) &
                  - TWELFTH*(s(i,j,k-2)+s(i,j,k+1))
             !
             ! limit sedgez
             !
             if ((sedgez(i,j,k)-s(i,j,k-1))*(s(i,j,k)-sedgez(i,j,k)) .lt. ZERO) then
                D2  = THREE*(s(i,j,k-1)-TWO*sedgez(i,j,k)+s(i,j,k))
                D2L = s(i,j,k-2)-TWO*s(i,j,k-1)+s(i,j,k)
                D2R = s(i,j,k-1)-TWO*s(i,j,k)+s(i,j,k+1)
                sgn = sign(ONE,D2)
                D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),ZERO)
                sedgez(i,j,k) = HALF*(s(i,j,k-1)+s(i,j,k)) - SIXTH*D2LIM
             end if
          end do
       end do
    end do
    !
    ! Use Colella 2008 limiters.
    !
    ! This is a new version of the algorithm to eliminate sensitivity to roundoff.
    !
    k = k3d
    do j=ilo2-1,ihi2+1
       do i=ilo1-1,ihi1+1

          alphap   = sedgez(i,j,k+1)-s(i,j,k)
          alpham   = sedgez(i,j,k  )-s(i,j,k)
          bigp     = abs(alphap).gt.TWO*abs(alpham)
          bigm     = abs(alpham).gt.TWO*abs(alphap)
          extremum = .false.

          if (alpham*alphap .ge. ZERO) then
             extremum = .true.
          else if (bigp .or. bigm) then
             !
             ! Possible extremum. We look at cell centered values and face
             ! centered values for a change in sign in the differences adjacent to
             ! the cell. We use the pair of differences whose minimum magnitude is the
             ! largest, and thus least susceptible to sensitivity to roundoff.
             !
             dafacem   = sedgez(i,j,k) - sedgez(i,j,k-1)
             dafacep   = sedgez(i,j,k+2) - sedgez(i,j,k+1)
             dabarm    = s(i,j,k) - s(i,j,k-1)
             dabarp    = s(i,j,k+1) - s(i,j,k)
             dafacemin = min(abs(dafacem),abs(dafacep))
             dabarmin  = min(abs(dabarm),abs(dabarp))
             if (dafacemin.ge.dabarmin) then
                dachkm = dafacem
                dachkp = dafacep
             else
                dachkm = dabarm
                dachkp = dabarp
             endif
             extremum = (dachkm*dachkp .le. ZERO)
          end if

          if (extremum) then
             D2     = SIX*(alpham + alphap)
             D2L    = s(i,j,k-2)-TWO*s(i,j,k-1)+s(i,j,k)
             D2R    = s(i,j,k)-TWO*s(i,j,k+1)+s(i,j,k+2)
             D2C    = s(i,j,k-1)-TWO*s(i,j,k)+s(i,j,k+1)
             sgn    = sign(ONE,D2)
             D2LIM  = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
             alpham = alpham*D2LIM/max(abs(D2),1.d-10)
             alphap = alphap*D2LIM/max(abs(D2),1.d-10)
          else
             if (bigp) then
                sgn   = sign(ONE,alpham)
                amax  = -alphap**2 / (4*(alpham + alphap))
                delam = s(i,j,k-1) - s(i,j,k)
                if (sgn*amax .ge. sgn*delam) then
                   if (sgn*(delam - alpham).ge.1.d-10) then
                      alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                   else 
                      alphap = -TWO*alpham
                   endif
                endif
             end if
             if (bigm) then
                sgn   = sign(ONE,alphap)
                amax  = -alpham**2 / (4*(alpham + alphap))
                delap = s(i,j,k+1) - s(i,j,k)
                if (sgn*amax .ge. sgn*delap) then
                   if (sgn*(delap - alphap).ge.1.d-10) then
                      alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                   else
                      alpham = -TWO*alphap
                   endif
                endif
             end if
          end if

          sm(i,j) = s(i,j,k) + alpham
          sp(i,j) = s(i,j,k) + alphap

          if (ppm_flatten_before_integrals > 0) then
             ! flatten the parabola AFTER doing the monotonization (note k = k3d here)
             sm(i,j) = flatn(i,j,k3d)*sm(i,j) + (ONE-flatn(i,j,k3d))*s(i,j,k3d)
             sp(i,j) = flatn(i,j,k3d)*sp(i,j) + (ONE-flatn(i,j,k3d))*s(i,j,k3d)
          endif

          !
          ! Compute z-component of Ip and Im.
          !
          s6    = SIX*s(i,j,k3d) - THREE*(sm(i,j)+sp(i,j))
          
          ! w-c wave
          sigma = abs(u(i,j,k3d,3)-cspd(i,j,k3d))*dt_over_a*dzinv
          
          if (u(i,j,k3d,3)-cspd(i,j,k3d) <= ZERO) then
             Ip(i,j,kc,3,1) = sp(i,j) 
          else
             Ip(i,j,kc,3,1) = sp(i,j) - &
                  HALF*sigma*(sp(i,j)-sm(i,j)-(ONE-TWO3RD*sigma)*s6)
          endif

          if (u(i,j,k3d,3)-cspd(i,j,k3d) >= ZERO) then
             Im(i,j,kc,3,1) = sm(i,j) 
          else
             Im(i,j,kc,3,1) = sm(i,j) + &
                  HALF*sigma*(sp(i,j)-sm(i,j)+(ONE-TWO3RD*sigma)*s6)
          endif

          ! w wave
          sigma = abs(u(i,j,k3d,3))*dt_over_a*dzinv

          if (u(i,j,k3d,3) <= ZERO) then
             Ip(i,j,kc,3,2) = sp(i,j)
          else
             Ip(i,j,kc,3,2) = sp(i,j) - &
                  HALF*sigma*(sp(i,j)-sm(i,j)-(ONE-TWO3RD*sigma)*s6)
          endif

          if (u(i,j,k3d,3) >= ZERO) then
             Im(i,j,kc,3,2) = sm(i,j)
          else
             Im(i,j,kc,3,2) = sm(i,j) + &
                  HALF*sigma*(sp(i,j)-sm(i,j)+(ONE-TWO3RD*sigma)*s6)
          endif

          ! w+c wave
          sigma = abs(u(i,j,k3d,3)+cspd(i,j,k3d))*dt_over_a*dzinv

          if (u(i,j,k3d,3)+cspd(i,j,k3d) <= ZERO) then
             Ip(i,j,kc,3,3) = sp(i,j) 
          else
             Ip(i,j,kc,3,3) = sp(i,j) - &
                  HALF*sigma*(sp(i,j)-sm(i,j)-(ONE-TWO3RD*sigma)*s6)
          endif

          if (u(i,j,k3d,3)+cspd(i,j,k3d) >= ZERO) then
             Im(i,j,kc,3,3) = sm(i,j) 
          else
             Im(i,j,kc,3,3) = sm(i,j) + &
                  HALF*sigma*(sp(i,j)-sm(i,j)+(ONE-TWO3RD*sigma)*s6)
          endif

       end do
    end do

    deallocate(dsvl,dsvlm,dsvlp,sp,sm,sedgez)

  end subroutine ppm_type2

end module ppm_module

