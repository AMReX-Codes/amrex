
subroutine amrex_interp_div_free_bfield (lo, hi, bx, bxlo, bxhi, by, bylo, byhi, bz, bzlo, bzhi, &
     cx, cxlo, cxhi, cy, cylo, cyhi, cz, czlo, czhi, dx, rr, use_limiter) bind(c)
  use amrex_fort_module, only : amrex_real, amrex_coarsen_intvect
  use amrex_mempool_module
  implicit none

  integer, intent(in) :: lo(3), hi(3), bxlo(3), bxhi(3), bylo(3), byhi(3), bzlo(3), bzhi(3), &
       cxlo(3), cxhi(3), cylo(3), cyhi(3), czlo(3), czhi(3), rr, use_limiter
  real(amrex_real), intent(in) :: dx(3)
  real(amrex_real), intent(inout) :: bx(bxlo(1):bxhi(1),bxlo(2):bxhi(2),bxlo(3):bxhi(3))
  real(amrex_real), intent(inout) :: by(bylo(1):byhi(1),bylo(2):byhi(2),bylo(3):byhi(3))
  real(amrex_real), intent(inout) :: bz(bzlo(1):bzhi(1),bzlo(2):bzhi(2),bzlo(3):bzhi(3))
  real(amrex_real), intent(in   ) :: cx(cxlo(1):cxhi(1),cxlo(2):cxhi(2),cxlo(3):cxhi(3))
  real(amrex_real), intent(in   ) :: cy(cylo(1):cyhi(1),cylo(2):cyhi(2),cylo(3):cyhi(3))
  real(amrex_real), intent(in   ) :: cz(czlo(1):czhi(1),czlo(2):czhi(2),czlo(3):czhi(3))

  integer :: i,j,k, ir, jr, kr, ii, jj, kk, clo(3), chi(3)
  real(amrex_real) :: dxinv(3), c1, c2, c3, x, y, z, rrinv
  real(amrex_real), dimension(:,:,:), contiguous, pointer :: dCxdy, dCxdz, dCydx, dCydz, dCzdx, dCzdy
  real(amrex_real) :: coe_a0, coe_ax, coe_ay, coe_az, coe_axx, coe_axy, coe_axz
  real(amrex_real) :: coe_b0, coe_bx, coe_by, coe_bz, coe_bxy, coe_byy, coe_byz
  real(amrex_real) :: coe_c0, coe_cx, coe_cy, coe_cz, coe_cxz, coe_cyz, coe_czz
  real(amrex_real), parameter :: theta = 2.d0 ! 1: minmod, 2: MC
  real(amrex_real), parameter :: one = 1.0_amrex_real
  
  clo = amrex_coarsen_intvect(3,lo,rr);
  chi = amrex_coarsen_intvect(3,hi,rr);

  dxinv = 1.d0/dx
  rrinv = 1.d0/rr

  call bl_allocate(dCxdy, clo(1),chi(1)+1,clo(2),chi(2)  ,clo(3),chi(3)  )
  call bl_allocate(dCxdz, clo(1),chi(1)+1,clo(2),chi(2)  ,clo(3),chi(3)  )
  call bl_allocate(dCydx, clo(1),chi(1)  ,clo(2),chi(2)+1,clo(3),chi(3)  )
  call bl_allocate(dCydz, clo(1),chi(1)  ,clo(2),chi(2)+1,clo(3),chi(3)  )
  call bl_allocate(dCzdx, clo(1),chi(1)  ,clo(2),chi(2)  ,clo(3),chi(3)+1)
  call bl_allocate(dCzdy, clo(1),chi(1)  ,clo(2),chi(2)  ,clo(3),chi(3)+1)

  do       k = clo(3), chi(3)
     do    j = clo(2), chi(2)
        if (use_limiter .eq. 0) then
           do i = clo(1), chi(1)+1
              dCxdy(i,j,k) = 0.5d0*(cx(i,j+1,k)-cx(i,j-1,k))
              dCxdz(i,j,k) = 0.5d0*(cx(i,j,k+1)-cx(i,j,k-1))
           end do
        else
           do i = clo(1), chi(1)+1
              c1 = theta*(cx(i,j  ,k)-cx(i,j-1,k))
              c2 = 0.5d0*(cx(i,j+1,k)-cx(i,j-1,k))
              c3 = theta*(cx(i,j+1,k)-cx(i,j  ,k))
              dCxdy(i,j,k) = 0.25d0*(sign(one,c1)+sign(one,c2))*(sign(one,c1)+sign(one,c3)) &
                   *min(abs(c1),abs(c2),abs(c3))
              
              c1 = theta*(cx(i,j,k  )-cx(i,j,k-1))
              c2 = 0.5d0*(cx(i,j,k+1)-cx(i,j,k-1))
              c3 = theta*(cx(i,j,k+1)-cx(i,j,k  ))
              dCxdz(i,j,k) = 0.25d0*(sign(one,c1)+sign(one,c2))*(sign(one,c1)+sign(one,c3)) &
                   *min(abs(c1),abs(c2),abs(c3))
           end do
        end if
     end do
  end do

  do       k = clo(3), chi(3)
     do    j = clo(2), chi(2)+1
        if (use_limiter .eq. 0) then
           do i = clo(1), chi(1)
              dCydx(i,j,k) = 0.5d0*(cy(i+1,j,k)-cy(i-1,j,k))
              dCydz(i,j,k) = 0.5d0*(cy(i,j,k+1)-cy(i,j,k-1))
           end do
        else
           do i = clo(1), chi(1)
              c1 = theta*(cy(i  ,j,k)-cy(i-1,j,k))
              c2 = 0.5d0*(cy(i+1,j,k)-cy(i-1,j,k))
              c3 = theta*(cy(i+1,j,k)-cy(i  ,j,k))
              dCydx(i,j,k) = 0.25d0*(sign(one,c1)+sign(one,c2))*(sign(one,c1)+sign(one,c3)) &
                   *min(abs(c1),abs(c2),abs(c3))
              
              c1 = theta*(cy(i,j,k  )-cy(i,j,k-1))
              c2 = 0.5d0*(cy(i,j,k+1)-cy(i,j,k-1))
              c3 = theta*(cy(i,j,k+1)-cy(i,j,k  ))
              dCydz(i,j,k) = 0.25d0*(sign(one,c1)+sign(one,c2))*(sign(one,c1)+sign(one,c3)) &
                   *min(abs(c1),abs(c2),abs(c3))
           end do
        end if
     end do
  end do

  do       k = clo(3), chi(3)+1
     do    j = clo(2), chi(2)
        if (use_limiter .eq. 0) then
           do i = clo(1), chi(1)
              dCzdx(i,j,k) = 0.5d0*(cz(i+1,j,k)-cz(i-1,j,k))
              dCzdy(i,j,k) = 0.5d0*(cz(i,j+1,k)-cz(i,j-1,k))
           end do
        else
           do i = clo(1), chi(1)
              c1 = theta*(cz(i  ,j,k)-cz(i-1,j,k))
              c2 = 0.5d0*(cz(i+1,j,k)-cz(i-1,j,k))
              c3 = theta*(cz(i+1,j,k)-cz(i  ,j,k))
              dCzdx(i,j,k) = 0.25d0*(sign(one,c1)+sign(one,c2))*(sign(one,c1)+sign(one,c3)) &
                   *min(abs(c1),abs(c2),abs(c3))
              
              c1 = theta*(cz(i,j  ,k)-cz(i,j-1,k))
              c2 = 0.5d0*(cz(i,j+1,k)-cz(i,j-1,k))
              c3 = theta*(cz(i,j+1,k)-cz(i,j  ,k))
              dCzdy(i,j,k) = 0.25d0*(sign(one,c1)+sign(one,c2))*(sign(one,c1)+sign(one,c3)) &
                   *min(abs(c1),abs(c2),abs(c3))
           end do
        end if
     end do
  end do

  do       k = clo(3), chi(3)
     do    j = clo(2), chi(2)
        do i = clo(1), chi(1)
           coe_ax = ( Cx  (i+1,j,k)- Cx  (i,j,k))*      dxinv(1)
           coe_ay = (dCxdy(i+1,j,k)+dCxdy(i,j,k))*0.5d0*dxinv(2)
           coe_az = (dCxdz(i+1,j,k)+dCxdz(i,j,k))*0.5d0*dxinv(3)

           coe_bx = (dCydx(i,j+1,k)+dCydx(i,j,k))*0.5d0*dxinv(1)
           coe_by = ( Cy  (i,j+1,k)- Cy  (i,j,k))*      dxinv(2)
           coe_bz = (dCydz(i,j+1,k)+dCydz(i,j,k))*0.5d0*dxinv(3)

           coe_cx = (dCzdx(i,j,k+1)+dCzdx(i,j,k))*0.5d0*dxinv(1)
           coe_cy = (dCzdy(i,j,k+1)+dCzdy(i,j,k))*0.5d0*dxinv(2)
           coe_cz = ( Cz  (i,j,k+1)- Cz  (i,j,k))*      dxinv(3)

           coe_axy = (dCxdy(i+1,j,k)-dCxdy(i,j,k))*dxinv(1)*dxinv(2)
           coe_axz = (dCxdz(i+1,j,k)-dCxdz(i,j,k))*dxinv(1)*dxinv(3)

           coe_bxy = (dCydx(i,j+1,k)-dCydx(i,j,k))*dxinv(1)*dxinv(2)
           coe_byz = (dCydz(i,j+1,k)-dCydz(i,j,k))*dxinv(2)*dxinv(3)

           coe_cxz = (dCzdx(i,j,k+1)-dCzdx(i,j,k))*dxinv(1)*dxinv(3)
           coe_cyz = (dCzdy(i,j,k+1)-dCzdy(i,j,k))*dxinv(2)*dxinv(3)

           coe_axx = -0.5d0*(coe_bxy+coe_cxz)
           coe_byy = -0.5d0*(coe_axy+coe_cyz)
           coe_czz = -0.5d0*(coe_axz+coe_byz)

           coe_a0 = 0.5d0*(Cx(i+1,j,k)+Cx(i,j,k)) - coe_axx*0.25d0*dx(1)*dx(1)
           coe_b0 = 0.5d0*(Cy(i,j+1,k)+Cy(i,j,k)) - coe_byy*0.25d0*dx(2)*dx(2)
           coe_c0 = 0.5d0*(Cz(i,j,k+1)+Cz(i,j,k)) - coe_czz*0.25d0*dx(3)*dx(3)

           ir = i*rr
           jr = j*rr
           kr = k*rr

           do kk = 0, rr-1
              z = (kk*rrinv-0.25d0) * dx(3)
              do jj = 0, rr-1
                 y = (jj*rrinv-0.25d0) * dx(2)
                 do ii = 0, rr
                    x = (ii*rrinv-0.5d0) * dx(1)
                    Bx(ir+ii,jr+jj,kr+kk) = coe_a0 + coe_ax*x + coe_ay*y + coe_az*z &
                         + coe_axx*x*x + coe_axy*x*y + coe_axz*x*z
                 end do
              end do
           end do

           do kk = 0, rr-1
              z = (kk*rrinv-0.25d0) * dx(3)
              do jj = 0, rr
                 y = (jj*rrinv-0.5d0) * dx(2)
                 do ii = 0, rr-1
                    x = (ii*rrinv-0.25d0) * dx(1)
                    By(ir+ii,jr+jj,kr+kk) = coe_b0 + coe_bx*x + coe_by*y + coe_bz*z &
                         + coe_bxy*x*y + coe_byy*y*y + coe_byz*y*z
                 end do
              end do
           end do

           do kk = 0, rr
              z = (kk*rrinv-0.5d0) * dx(3)
              do jj = 0, rr-1
                 y = (jj*rrinv-0.25d0) * dx(2)
                 do ii = 0, rr-1
                    x = (ii*rrinv-0.25d0) * dx(1)
                    Bz(ir+ii,jr+jj,kr+kk) = coe_c0 + coe_cx*x + coe_cy*y + coe_cz*z &
                         + coe_cxz*x*z + coe_cyz*y*z + coe_czz*z*z
                 end do
              end do
           end do

        end do
     end do
  end do

  call bl_deallocate(dCxdy)
  call bl_deallocate(dCxdz)
  call bl_deallocate(dCydx)
  call bl_deallocate(dCydz)
  call bl_deallocate(dCzdx)
  call bl_deallocate(dCzdy)

end subroutine amrex_interp_div_free_bfield


subroutine amrex_interp_efield (lo, hi, ex, exlo, exhi, ey, eylo, eyhi, ez, ezlo, ezhi, &
     cx, cxlo, cxhi, cy, cylo, cyhi, cz, czlo, czhi, rr, use_limiter) bind(c)
  use amrex_fort_module, only : amrex_real, amrex_coarsen_intvect
  use amrex_mempool_module
  implicit none

  integer, intent(in) :: lo(3), hi(3), exlo(3), exhi(3), eylo(3), eyhi(3), ezlo(3), ezhi(3), &
       cxlo(3), cxhi(3), cylo(3), cyhi(3), czlo(3), czhi(3), rr, use_limiter

  real(amrex_real), intent(inout) :: ex(exlo(1):exhi(1),exlo(2):exhi(2),exlo(3):exhi(3))
  real(amrex_real), intent(inout) :: ey(eylo(1):eyhi(1),eylo(2):eyhi(2),eylo(3):eyhi(3))
  real(amrex_real), intent(inout) :: ez(ezlo(1):ezhi(1),ezlo(2):ezhi(2),ezlo(3):ezhi(3))
  real(amrex_real), intent(in   ) :: cx(cxlo(1):cxhi(1),cxlo(2):cxhi(2),cxlo(3):cxhi(3))
  real(amrex_real), intent(in   ) :: cy(cylo(1):cyhi(1),cylo(2):cyhi(2),cylo(3):cyhi(3))
  real(amrex_real), intent(in   ) :: cz(czlo(1):czhi(1),czlo(2):czhi(2),czlo(3):czhi(3))

  integer :: i, j, k, ic, jc, kc
  integer :: clo(3), chi(3)
  real(amrex_real) :: c1, c2, c3, dc
  real(amrex_real), pointer, contiguous :: tmpx(:,:,:), tmpy(:,:,:), tmpz(:,:,:)
  real(amrex_real), parameter :: theta = 2.d0 ! 1: minmod, 2: MC
  real(amrex_real), parameter :: one = 1.0_amrex_real

  clo = amrex_coarsen_intvect(3,lo,rr);
  chi = amrex_coarsen_intvect(3,hi,rr);

  call bl_allocate(tmpx,clo(1)-1,chi(1)+1, lo(2)  , hi(2)+1, lo(3)  , hi(3)+1)
  call bl_allocate(tmpy, lo(1)  , hi(1)+1,clo(2)-1,chi(2)+1, lo(3)  , hi(3)+1)
  call bl_allocate(tmpz, lo(1)  , hi(1)+1, lo(2)  , hi(2)+1,clo(3)-1,chi(3)+1)

  do       kc = clo(3)  , chi(3)+1
     do    jc = clo(2)  , chi(2)+1
        do ic = clo(1)-1, chi(1)+1
           tmpx(ic,jc*2,kc*2) = cx(ic,jc,kc)
        end do
        if (jc .lt. chi(2)+1) then
           do ic = clo(1)-1, chi(1)+1
              tmpx(ic,jc*2+1,kc*2) = 0.5d0*(cx(ic,jc,kc)+cx(ic,jc+1,kc))
           end do
        end if
     end do
  end do

  do       k  =  lo(3)+1,  hi(3)  , 2
     do    j  =  lo(2)  ,  hi(2)+1
        do ic = clo(1)-1, chi(1)+1
           tmpx(ic,j,k) = 0.5d0*(tmpx(ic,j,k-1)+tmpx(ic,j,k+1))
        end do
     end do
  end do
           
  do       k  =  lo(3), hi(3)+1
     do    j  =  lo(2), hi(2)+1
        if (use_limiter .eq. 0) then
           do ic = clo(1),chi(1)
              dc = 0.5d0*(tmpx(ic+1,j,k) - tmpx(ic-1,j,k))              
              ex(ic*2  ,j,k) = tmpx(ic,j,k) - 0.25d0*dc
              ex(ic*2+1,j,k) = tmpx(ic,j,k) + 0.25d0*dc
           end do
        else
           do ic = clo(1),chi(1)
              c1 = theta*(tmpx(ic  ,j,k) - tmpx(ic-1,j,k))
              c2 = 0.5d0*(tmpx(ic+1,j,k) - tmpx(ic-1,j,k))
              c3 = theta*(tmpx(ic+1,j,k) - tmpx(ic  ,j,k))
              dc = 0.25d0*(sign(one,c1)+sign(one,c2))*(sign(one,c1)+sign(one,c3)) &
                   *min(abs(c1),abs(c2),abs(c3))
              
              ex(ic*2  ,j,k) = tmpx(ic,j,k) - 0.25d0*dc
              ex(ic*2+1,j,k) = tmpx(ic,j,k) + 0.25d0*dc
           end do
        end if
     end do
  end do

  ! ey
  do       kc = clo(3)  , chi(3)+1
     do    jc = clo(2)-1, chi(2)+1
        do ic = clo(1)  , chi(1)
           tmpy(ic*2  ,jc,kc*2) = cy(ic,jc,kc)
           tmpy(ic*2+1,jc,kc*2) = 0.5d0*(cy(ic,jc,kc)+cy(ic+1,jc,kc))
        end do
        ic = chi(1)+1
        tmpy(ic*2,jc,kc*2) = cy(ic,jc,kc)        
     end do
  end do

  do       k  =  lo(3)+1,  hi(3)  , 2
     do    jc = clo(2)-1, chi(2)+1
        do i  =  lo(1)  ,  hi(1)+1
           tmpy(i,jc,k) = 0.5d0*(tmpy(i,jc,k-1)+tmpy(i,jc,k+1))
        end do
     end do
  end do
  
  do       k  =  lo(3),  hi(3)+1
     do    jc = clo(2), chi(2)
        if (use_limiter .eq. 0) then
           do i  =  lo(1),  hi(1)+1
              dc = 0.5d0*(tmpy(i,jc+1,k) - tmpy(i,jc-1,k))              
              ey(i,jc*2  ,k) = tmpy(i,jc,k) - 0.25d0*dc
              ey(i,jc*2+1,k) = tmpy(i,jc,k) + 0.25d0*dc
           end do
        else
           do i  =  lo(1),  hi(1)+1
              c1 = theta*(tmpy(i,jc  ,k) - tmpy(i,jc-1,k))
              c2 = 0.5d0*(tmpy(i,jc+1,k) - tmpy(i,jc-1,k))
              c3 = theta*(tmpy(i,jc+1,k) - tmpy(i,jc  ,k))
              dc = 0.25d0*(sign(one,c1)+sign(one,c2))*(sign(one,c1)+sign(one,c3)) &
                   *min(abs(c1),abs(c2),abs(c3))
              
              ey(i,jc*2  ,k) = tmpy(i,jc,k) - 0.25d0*dc
              ey(i,jc*2+1,k) = tmpy(i,jc,k) + 0.25d0*dc
           end do
        end if
     end do
  end do

  
  ! ez
  do       kc = clo(3)-1, chi(3)+1
     do    jc = clo(2)  , chi(2)+1
        do ic = clo(1)  , chi(1)
           tmpz(ic*2  ,jc*2,kc) = cz(ic,jc,kc)
           tmpz(ic*2+1,jc*2,kc) = 0.5d0*(cz(ic,jc,kc)+cz(ic+1,jc,kc))
        end do
        ic = chi(1)+1
        tmpz(ic*2,jc*2,kc) = cz(ic,jc,kc)
     end do

     do    j = lo(2)+1, hi(2)  , 2
        do i = lo(1)  , hi(1)+1
           tmpz(i,j,kc) = 0.5d0*(tmpz(i,j-1,kc)+tmpz(i,j+1,kc))
        end do
     end do
  end do

  do       kc = clo(3), chi(3)
     do    j  =  lo(2),  hi(2)+1
        if (use_limiter .eq. 0) then
           do i  =  lo(1),  hi(1)+1
              dc = 0.5d0*(tmpz(i,j,kc+1) - tmpz(i,j,kc-1))              
              ez(i,j,kc*2  ) = tmpz(i,j,kc) - 0.25d0*dc
              ez(i,j,kc*2+1) = tmpz(i,j,kc) + 0.25d0*dc
           end do
        else
           do i  =  lo(1),  hi(1)+1
              c1 = theta*(tmpz(i,j,kc  ) - tmpz(i,j,kc-1))
              c2 = 0.5d0*(tmpz(i,j,kc+1) - tmpz(i,j,kc-1))
              c3 = theta*(tmpz(i,j,kc+1) - tmpz(i,j,kc  ))
              dc = 0.25d0*(sign(one,c1)+sign(one,c2))*(sign(one,c1)+sign(one,c3)) &
                   *min(abs(c1),abs(c2),abs(c3))
              
              ez(i,j,kc*2  ) = tmpz(i,j,kc) - 0.25d0*dc
              ez(i,j,kc*2+1) = tmpz(i,j,kc) + 0.25d0*dc
           end do
        end if
     end do
  end do
  
  call bl_deallocate(tmpx)
  call bl_deallocate(tmpy)
  call bl_deallocate(tmpz)
  
end subroutine amrex_interp_efield
