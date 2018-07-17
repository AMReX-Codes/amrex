
subroutine amrex_interp_div_free_bfield (lo, hi, bx, bxlo, bxhi, by, bylo, byhi, &
     cx, cxlo, cxhi, cy, cylo, cyhi, dx, rr, use_limiter) bind(c)
  use amrex_fort_module, only : amrex_real, amrex_coarsen_intvect
  use amrex_mempool_module
  implicit none

  integer, intent(in) :: lo(2), hi(2), bxlo(2), bxhi(2), bylo(2), byhi(2), &
       cxlo(2), cxhi(2), cylo(2), cyhi(2), rr, use_limiter
  real(amrex_real), intent(in) :: dx(2)
  real(amrex_real), intent(inout) :: bx(bxlo(1):bxhi(1),bxlo(2):bxhi(2))
  real(amrex_real), intent(inout) :: by(bylo(1):byhi(1),bylo(2):byhi(2))
  real(amrex_real), intent(in   ) :: cx(cxlo(1):cxhi(1),cxlo(2):cxhi(2))
  real(amrex_real), intent(in   ) :: cy(cylo(1):cyhi(1),cylo(2):cyhi(2))

  integer :: i,j, ir, jr, ii, jj, clo(2), chi(2)
  
  real(amrex_real) :: dxinv(2), c1, c2, c3, x, y, rrinv
  real(amrex_real), dimension(:,:), contiguous, pointer :: dCxdy, dCydx
  real(amrex_real) :: coe_a0, coe_ax, coe_ay, coe_axy, coe_axx
  real(amrex_real) :: coe_b0, coe_bx, coe_by, coe_bxy, coe_byy
  real(amrex_real), parameter :: theta = 2.d0 ! 1: minmod, 2: MC
  real(amrex_real), parameter :: one = 1._amrex_real

  clo = amrex_coarsen_intvect(2,lo,rr);
  chi = amrex_coarsen_intvect(2,hi,rr);

  dxinv = 1.d0/dx
  rrinv = 1.d0/rr

  call bl_allocate(dCxdy, clo(1),chi(1)+1,clo(2),chi(2)  )
  call bl_allocate(dCydx, clo(1),chi(1)  ,clo(2),chi(2)+1)

  do    j = clo(2), chi(2)
     if (use_limiter .eq. 0) then
        do i = clo(1), chi(1)+1
           dCxdy(i,j) = 0.5d0*(cx(i,j+1)-cx(i,j-1))
        end do
     else
        do i = clo(1), chi(1)+1
           c1 = theta*(cx(i,j  )-cx(i,j-1))
           c2 = 0.5d0*(cx(i,j+1)-cx(i,j-1))
           c3 = theta*(cx(i,j+1)-cx(i,j  ))
           dCxdy(i,j) = 0.25d0*(sign(one,c1)+sign(one,c2))*(sign(one,c1)+sign(one,c3)) &
                *min(abs(c1),abs(c2),abs(c3))
        end do
     end if
  end do

  do    j = clo(2), chi(2)+1
     do i = clo(1), chi(1)
        if (use_limiter .eq. 0) then
           dCydx(i,j) = 0.5d0*(cy(i+1,j)-cy(i-1,j))
        else
           c1 = theta*(cy(i  ,j)-cy(i-1,j))
           c2 = 0.5d0*(cy(i+1,j)-cy(i-1,j))
           c3 = theta*(cy(i+1,j)-cy(i  ,j))
           dCydx(i,j) = 0.25d0*(sign(one,c1)+sign(one,c2))*(sign(one,c1)+sign(one,c3)) &
                *min(abs(c1),abs(c2),abs(c3))
        endif
     end do
  end do

  do    j = clo(2), chi(2)
     do i = clo(1), chi(1)
        coe_ax = (Cx(i+1,j)- Cx(i,j))*dxinv(1)
        coe_by = (Cy(i,j+1)- Cy(i,j))*dxinv(2)

        coe_ay = (dCxdy(i+1,j)+dCxdy(i,j))*0.5d0*dxinv(2)
        coe_bx = (dCydx(i,j+1)+dCydx(i,j))*0.5d0*dxinv(1)

        coe_axy = (dCxdy(i+1,j)-dCxdy(i,j))*dxinv(1)*dxinv(2)
        coe_bxy = (dCydx(i,j+1)-dCydx(i,j))*dxinv(1)*dxinv(2)

        coe_axx = -0.5d0*coe_bxy
        coe_byy = -0.5d0*coe_axy

        coe_a0 = 0.5d0*(Cx(i+1,j)+Cx(i,j)) - coe_axx*0.25d0*dx(1)*dx(1)
        coe_b0 = 0.5d0*(Cy(i,j+1)+Cy(i,j)) - coe_byy*0.25d0*dx(2)*dx(2)

        ir = i*rr
        jr = j*rr

        do jj = 0, rr-1
           y = (jj*rrinv-0.25d0) * dx(2)
           do ii = 0, rr
              x = (ii*rrinv-0.5d0) * dx(1)
              Bx(ir+ii,jr+jj) = coe_a0 + coe_ax*x + coe_ay*y &
                   + coe_axx*x*x + coe_axy*x*y
           end do
        end do

        do jj = 0, rr
           y = (jj*rrinv-0.5d0) * dx(2)
           do ii = 0, rr-1
              x = (ii*rrinv-0.25d0) * dx(1)
              By(ir+ii,jr+jj) = coe_b0 + coe_bx*x + coe_by*y &
                   + coe_bxy*x*y + coe_byy*y*y
           end do
        end do
     end do
  end do

  call bl_deallocate(dCxdy)
  call bl_deallocate(dCydx)

end subroutine amrex_interp_div_free_bfield

subroutine amrex_interp_efield (lo, hi, ex, exlo, exhi, ey, eylo, eyhi, &
     cx, cxlo, cxhi, cy, cylo, cyhi, rr, use_limiter) bind(c)
  use amrex_fort_module, only : amrex_real, amrex_coarsen_intvect
  use amrex_mempool_module
  implicit none

  integer, intent(in) :: lo(2), hi(2), exlo(2), exhi(2), eylo(2), eyhi(2), &
       cxlo(2), cxhi(2), cylo(2), cyhi(2), rr, use_limiter
  real(amrex_real), intent(inout) :: ex(exlo(1):exhi(1),exlo(2):exhi(2))
  real(amrex_real), intent(inout) :: ey(eylo(1):eyhi(1),eylo(2):eyhi(2))
  real(amrex_real), intent(in   ) :: cx(cxlo(1):cxhi(1),cxlo(2):cxhi(2))
  real(amrex_real), intent(in   ) :: cy(cylo(1):cyhi(1),cylo(2):cyhi(2))

  integer :: i,j,ic,jc
  integer :: clo(2), chi(2)
  real(amrex_real) :: c1, c2, c3, dc
  real(amrex_real), pointer, contiguous :: tmpx(:,:), tmpy(:,:)
  real(amrex_real), parameter :: theta = 2.d0 ! 1: minmod, 2: MC
  real(amrex_real), parameter :: one = 1._amrex_real

  clo = amrex_coarsen_intvect(2,lo,rr);
  chi = amrex_coarsen_intvect(2,hi,rr);

  call bl_allocate(tmpx,clo(1)-1,chi(1)+1, lo(2)  , hi(2)+1)
  call bl_allocate(tmpy, lo(1)  , hi(1)+1,clo(2)-1,chi(2)+1)

  do    jc = clo(2)  , chi(2)+1
     do ic = clo(1)-1, chi(1)+1
        tmpx(ic,jc*2) = cx(ic,jc)
     end do
     if (jc .lt. chi(2)+1) then
        do ic = clo(1)-1, chi(1)+1
           tmpx(ic,jc*2+1) = 0.5d0*(cx(ic,jc)+cx(ic,jc+1))
        end do
     end if
  end do

  do    j  =  lo(2), hi(2)+1
     if (use_limiter .eq. 0) then
        do ic = clo(1),chi(1)
           dc = 0.5d0*(tmpx(ic+1,j) - tmpx(ic-1,j))           
           ex(ic*2  ,j) = tmpx(ic,j) - 0.25d0*dc
           ex(ic*2+1,j) = tmpx(ic,j) + 0.25d0*dc
        end do
     else
        do ic = clo(1),chi(1)
           c1 = theta*(tmpx(ic  ,j) - tmpx(ic-1,j))
           c2 = 0.5d0*(tmpx(ic+1,j) - tmpx(ic-1,j))
           c3 = theta*(tmpx(ic+1,j) - tmpx(ic  ,j))
           dc = 0.25d0*(sign(one,c1)+sign(one,c2))*(sign(one,c1)+sign(one,c3)) &
                *min(abs(c1),abs(c2),abs(c3))
           
           ex(ic*2  ,j) = tmpx(ic,j) - 0.25d0*dc
           ex(ic*2+1,j) = tmpx(ic,j) + 0.25d0*dc
        end do
     end if
  end do

  ! ey
  do    jc = clo(2)-1, chi(2)+1
     do ic = clo(1)  , chi(1)
        tmpy(ic*2  ,jc) = cy(ic,jc)
        tmpy(ic*2+1,jc) = 0.5d0*(cy(ic,jc)+cy(ic+1,jc))
     end do
     ic = chi(1)+1
     tmpy(ic*2,jc) = cy(ic,jc)
  end do

  do    jc = clo(2), chi(2)
     if (use_limiter .eq. 0) then
        do i  =  lo(1),  hi(1)+1
           dc = 0.5d0*(tmpy(i,jc+1) - tmpy(i,jc-1))           
           ey(i,jc*2  ) = tmpy(i,jc) - 0.25d0*dc
           ey(i,jc*2+1) = tmpy(i,jc) + 0.25d0*dc
        end do
     else
        do i  =  lo(1),  hi(1)+1
           c1 = theta*(tmpy(i,jc  ) - tmpy(i,jc-1))
           c2 = 0.5d0*(tmpy(i,jc+1) - tmpy(i,jc-1))
           c3 = theta*(tmpy(i,jc+1) - tmpy(i,jc  ))
           dc = 0.25d0*(sign(one,c1)+sign(one,c2))*(sign(one,c1)+sign(one,c3)) &
                *min(abs(c1),abs(c2),abs(c3))
           
           ey(i,jc*2  ) = tmpy(i,jc) - 0.25d0*dc
           ey(i,jc*2+1) = tmpy(i,jc) + 0.25d0*dc
        end do
     end if
  end do

  call bl_deallocate(tmpx)
  call bl_deallocate(tmpy)

end subroutine amrex_interp_efield

subroutine amrex_interp_cc_bfield (lo, hi, by, bylo, byhi, cy, cylo, cyhi, rr, use_limiter) bind(c)
  use amrex_fort_module, only : amrex_real, amrex_coarsen_intvect
  use amrex_mempool_module
  implicit none

  integer, intent(in) :: lo(2), hi(2), bylo(2), byhi(2), cylo(2), cyhi(2), rr, use_limiter
  real(amrex_real), intent(inout) :: by(bylo(1):byhi(1),bylo(2):byhi(2))
  real(amrex_real), intent(in   ) :: cy(cylo(1):cyhi(1),cylo(2):cyhi(2))

  integer :: i, ic, jc
  integer :: clo(2), chi(2)
  real(amrex_real), pointer, contiguous :: tmp(:,:)
  real(amrex_real) :: c1, c2, c3, dc
  real(amrex_real), parameter :: theta = 2.d0 ! 1: minmod, 2: MC
  real(amrex_real), parameter :: one = 1._amrex_real
  
  clo = amrex_coarsen_intvect(2,lo,rr);
  chi = amrex_coarsen_intvect(2,hi,rr);

  call bl_allocate(tmp,lo(1),hi(1),clo(2)-1,chi(2)+1)

  do    jc = clo(2)-1, chi(2)+1
     if (use_limiter .eq. 0) then
        do ic = clo(1)  , chi(1)
           dc = 0.5d0*(cy(ic+1,jc) - cy(ic-1,jc))
           tmp(ic*2  ,jc) = cy(ic,jc) - 0.25d0*dc
           tmp(ic*2+1,jc) = cy(ic,jc) + 0.25d0*dc
        end do
     else
        do ic = clo(1)  , chi(1)
           c1 = theta*(cy(ic  ,jc) - cy(ic-1,jc))
           c2 = 0.5d0*(cy(ic+1,jc) - cy(ic-1,jc))
           c3 = theta*(cy(ic+1,jc) - cy(ic  ,jc))
           dc = 0.25d0*(sign(one,c1)+sign(one,c2))*(sign(one,c1)+sign(one,c3)) &
                *min(abs(c1),abs(c2),abs(c3))
           tmp(ic*2  ,jc) = cy(ic,jc) - 0.25d0*dc
           tmp(ic*2+1,jc) = cy(ic,jc) + 0.25d0*dc
        end do
     end if
  end do

  do    jc = clo(2), chi(2)
     if (use_limiter .eq. 0) then
        do i  =  lo(1),  hi(1)
           dc = 0.5d0*(tmp(i,jc+1) - tmp(i,jc-1))
           by(i,jc*2  ) = tmp(i,jc) - 0.25d0*dc
           by(i,jc*2+1) = tmp(i,jc) + 0.25d0*dc
        end do
     else
        do i  =  lo(1),  hi(1)
           c1 = theta*(tmp(i,jc  ) - tmp(i,jc-1))
           c2 = 0.5d0*(tmp(i,jc+1) - tmp(i,jc-1))
           c3 = theta*(tmp(i,jc+1) - tmp(i,jc  ))
           dc = 0.25d0*(sign(one,c1)+sign(one,c2))*(sign(one,c1)+sign(one,c3)) &
                *min(abs(c1),abs(c2),abs(c3))
           by(i,jc*2  ) = tmp(i,jc) - 0.25d0*dc
           by(i,jc*2+1) = tmp(i,jc) + 0.25d0*dc
        end do
     end if
  end do

  call bl_deallocate(tmp)

end subroutine amrex_interp_cc_bfield

subroutine amrex_interp_nd_efield (lo, hi, ey, eylo, eyhi, cy, cylo, cyhi, rr) bind(c)
  use amrex_fort_module, only : amrex_real, amrex_coarsen_intvect
  implicit none

  integer, intent(in) :: lo(2), hi(2), eylo(2), eyhi(2), cylo(2), cyhi(2), rr
  real(amrex_real), intent(inout) :: ey(eylo(1):eyhi(1),eylo(2):eyhi(2))
  real(amrex_real), intent(in   ) :: cy(cylo(1):cyhi(1),cylo(2):cyhi(2))

  integer :: i,j,clo(2),chi(2)
  
  clo = amrex_coarsen_intvect(2,lo,rr);
  chi = amrex_coarsen_intvect(2,hi,rr);

  do    j = clo(2), chi(2)+1
     do i = clo(1), chi(1)
        ey(i*2  ,j*2) = cy(i,j)
        ey(i*2+1,j*2) = 0.5d0*(cy(i,j)+cy(i+1,j))
     end do
     i = chi(1)+1
     ey(i*2,j*2) = cy(i,j)     
     if (j .le. chi(2)) then
        do i = clo(1), chi(1)
           ey(i*2  ,j*2+1) = 0.5d0*(cy(i,j)+cy(i,j+1))
           ey(i*2+1,j*2+1) = 0.25d0*(cy(i,j)+cy(i+1,j)+cy(i,j+1)+cy(i+1,j+1))
        end do
        i = chi(1)+1
        ey(i*2,j*2+1) = 0.5d0*(cy(i,j)+cy(i,j+1))
     end if
  end do

end subroutine amrex_interp_nd_efield
