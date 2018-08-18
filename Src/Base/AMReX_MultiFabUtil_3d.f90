subroutine amrex_fort_avg_nd_to_cc (lo, hi, ncomp, &
     cc, ccl1, ccl2, ccl3, cch1, cch2, cch3, &
     nd, ndl1, ndl2, ndl3, ndh1, ndh2, ndh3 ) bind(c)

  use amrex_fort_module, only : amrex_real
  implicit none
  integer          :: lo(3),hi(3), ncomp
  integer          :: ccl1, ccl2, ccl3, cch1, cch2, cch3
  integer          :: ndl1, ndl2, ndl3, ndh1, ndh2, ndh3
  real(amrex_real) :: cc(ccl1:cch1, ccl2:cch2, ccl3:cch3, ncomp)
  real(amrex_real) :: nd(ndl1:ndh1, ndl2:ndh2, ndl3:ndh3, ncomp)

  ! Local variables
  integer          :: i,j,k,n
  
  do n = 1, ncomp
     do k=lo(3),hi(3)
        do j=lo(2),hi(2)
           do i=lo(1),hi(1)
              cc(i,j,k,n) = 0.125_amrex_real * ( nd(i,j  ,k  ,n) + nd(i+1,j  ,k  ,n) &
                   +                             nd(i,j+1,k  ,n) + nd(i+1,j+1,k  ,n) &
                   +                             nd(i,j  ,k+1,n) + nd(i+1,j  ,k+1,n) &
                   +                             nd(i,j+1,k+1,n) + nd(i+1,j+1,k+1,n) )
           end do
        end do
     end do
  end do

end subroutine amrex_fort_avg_nd_to_cc


! ***************************************************************************************
! subroutine bl_avg_eg_to_cc 
! ***************************************************************************************

subroutine bl_avg_eg_to_cc (lo, hi, &
     cc, ccl1, ccl2, ccl3, cch1, cch2, cch3, &
     Ex, Exl1, Exl2, Exl3, Exh1, Exh2, Exh3, &
     Ey, Eyl1, Eyl2, Eyl3, Eyh1, Eyh2, Eyh3, &
     Ez, Ezl1, Ezl2, Ezl3, Ezh1, Ezh2, Ezh3)

  use amrex_fort_module, only : amrex_real
  implicit none
  integer          :: lo(3),hi(3)
  integer          :: ccl1, ccl2, ccl3, cch1, cch2, cch3
  integer          :: Exl1, Exl2, Exl3, Exh1, Exh2, Exh3
  integer          :: Eyl1, Eyl2, Eyl3, Eyh1, Eyh2, Eyh3
  integer          :: Ezl1, Ezl2, Ezl3, Ezh1, Ezh2, Ezh3
  real(amrex_real) :: cc(ccl1:cch1, ccl2:cch2, ccl3:cch3, 3)
  real(amrex_real) :: Ex(Exl1:Exh1, Exl2:Exh2, Exl3:Exh3)
  real(amrex_real) :: Ey(Eyl1:Eyh1, Eyl2:Eyh2, Eyl3:Eyh3)
  real(amrex_real) :: Ez(Ezl1:Ezh1, Ezl2:Ezh2, Ezl3:Ezh3)

  ! Local variables
  integer          :: i,j,k
  
  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)
           cc(i,j,k,1) = 0.25d0 * ( Ex(i,j,k) + Ex(i,j+1,k) + Ex(i,j,k+1) + Ex(i,j+1,k+1) )
           cc(i,j,k,2) = 0.25d0 * ( Ey(i,j,k) + Ey(i+1,j,k) + Ey(i,j,k+1) + Ey(i+1,j,k+1) )
           cc(i,j,k,3) = 0.25d0 * ( Ez(i,j,k) + Ez(i+1,j,k) + Ez(i,j+1,k) + Ez(i+1,j+1,k) )
        enddo
     enddo
  enddo

end subroutine bl_avg_eg_to_cc

! ***************************************************************************************
! subroutine bl_avg_fc_to_cc 
! ***************************************************************************************

subroutine bl_avg_fc_to_cc (lo, hi, &
     cc, ccl1, ccl2, ccl3, cch1, cch2, cch3, &
     fx, fxl1, fxl2, fxl3, fxh1, fxh2, fxh3, &
     fy, fyl1, fyl2, fyl3, fyh1, fyh2, fyh3, &
     fz, fzl1, fzl2, fzl3, fzh1, fzh2, fzh3, &
     dx, problo, coord_type)

  use amrex_fort_module, only : amrex_real
  implicit none
  integer          :: lo(3),hi(3), coord_type
  integer          :: ccl1, ccl2, ccl3, cch1, cch2, cch3
  integer          :: fxl1, fxl2, fxl3, fxh1, fxh2, fxh3
  integer          :: fyl1, fyl2, fyl3, fyh1, fyh2, fyh3
  integer          :: fzl1, fzl2, fzl3, fzh1, fzh2, fzh3
  real(amrex_real) :: cc(ccl1:cch1, ccl2:cch2, ccl3:cch3, 3)
  real(amrex_real) :: fx(fxl1:fxh1, fxl2:fxh2, fxl3:fxh3)
  real(amrex_real) :: fy(fyl1:fyh1, fyl2:fyh2, fyl3:fyh3)
  real(amrex_real) :: fz(fzl1:fzh1, fzl2:fzh2, fzl3:fzh3)
  real(amrex_real) :: dx(3), problo(3)

  ! Local variables
  integer          :: i,j,k
  
  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)
           cc(i,j,k,1) = 0.5d0 * ( fx(i,j,k) + fx(i+1,j,k) )
           cc(i,j,k,2) = 0.5d0 * ( fy(i,j,k) + fy(i,j+1,k) )
           cc(i,j,k,3) = 0.5d0 * ( fz(i,j,k) + fz(i,j,k+1) )
        enddo
     enddo
  enddo

end subroutine bl_avg_fc_to_cc

! ***************************************************************************************
! subroutine bl_avg_cc_to_fc 
! ***************************************************************************************

subroutine bl_avg_cc_to_fc (xlo, xhi, ylo, yhi, zlo, zhi, &
     fx, fxl1, fxl2, fxl3, fxh1, fxh2, fxh3, &
     fy, fyl1, fyl2, fyl3, fyh1, fyh2, fyh3, &
     fz, fzl1, fzl2, fzl3, fzh1, fzh2, fzh3, &
     cc, ccl1, ccl2, ccl3, cch1, cch2, cch3, &
     dx, problo, coord_type)

  use amrex_fort_module, only : amrex_real
  implicit none
  integer          :: xlo(3),xhi(3), ylo(3),yhi(3), zlo(3), zhi(3), coord_type
  integer          :: fxl1, fxl2, fxl3, fxh1, fxh2, fxh3
  integer          :: fyl1, fyl2, fyl3, fyh1, fyh2, fyh3
  integer          :: fzl1, fzl2, fzl3, fzh1, fzh2, fzh3
  integer          :: ccl1, ccl2, ccl3, cch1, cch2, cch3
  real(amrex_real) :: cc(ccl1:cch1, ccl2:cch2, ccl3:cch3)
  real(amrex_real) :: fx(fxl1:fxh1, fxl2:fxh2, fxl3:fxh3)
  real(amrex_real) :: fy(fyl1:fyh1, fyl2:fyh2, fyl3:fyh3)
  real(amrex_real) :: fz(fzl1:fzh1, fzl2:fzh2, fzl3:fzh3)
  real(amrex_real) :: dx(3), problo(3)

  ! Local variables
  integer          :: i,j,k
  
  do k=xlo(3),xhi(3)
     do j=xlo(2),xhi(2)
        do i=xlo(1),xhi(1)
           fx(i,j,k) = 0.5d0 * (cc(i-1,j,k) + cc(i,j,k))
        enddo
     enddo
  end do

  do k=ylo(3),yhi(3)
     do j=ylo(2),yhi(2)
        do i=ylo(1),yhi(1)
           fy(i,j,k) = 0.5d0 * (cc(i,j-1,k) + cc(i,j,k))
        enddo
     enddo
  end do

  do k=zlo(3),zhi(3)
     do j=zlo(2),zhi(2)
        do i=zlo(1),zhi(1)
           fz(i,j,k) = 0.5d0 * (cc(i,j,k-1) + cc(i,j,k))
        enddo
     enddo
  end do

end subroutine bl_avg_cc_to_fc

! ***************************************************************************************
! subroutine bl_avgdown_faces
! ***************************************************************************************

subroutine bl_avgdown_faces (lo, hi, &
     f, f_l1, f_l2, f_l3, f_h1, f_h2, f_h3, &
     c, c_l1, c_l2, c_l3, c_h1, c_h2, c_h3, &
     ratio,idir,nc)

  use amrex_fort_module, only : amrex_real
  implicit none
  integer          :: lo(3),hi(3)
  integer          :: f_l1, f_l2, f_l3, f_h1, f_h2, f_h3
  integer          :: c_l1, c_l2, c_l3, c_h1, c_h2, c_h3
  integer          :: ratio(3), idir, nc
  real(amrex_real) :: f(f_l1:f_h1, f_l2:f_h2, f_l3:f_h3, nc)
  real(amrex_real) :: c(c_l1:c_h1, c_l2:c_h2, c_l3:c_h3, nc)

  ! Local variables
  integer :: i, j, k, n, facx, facy, facz, iref, jref, kref, ii, jj, kk
  real(amrex_real) :: facInv

  facx = ratio(1)
  facy = ratio(2)
  facz = ratio(3)

  if (idir .eq. 0) then

     facInv = 1.d0 / (facy*facz)

     do n = 1, nc
        do k        = lo(3), hi(3)
           kk       = k * facz
           do j     = lo(2), hi(2)
              jj    = j * facy
              do i  = lo(1), hi(1)
                 ii = i * facx
                 c(i,j,k,n) = 0.d0
                 do    kref = 0, facz-1
                    do jref = 0, facy-1
                       c(i,j,k,n) = c(i,j,k,n) + f(ii,jj+jref,kk+kref,n)
                    end do
                 end do
                 c(i,j,k,n) = c(i,j,k,n) * facInv
              end do
           end do
        end do
     end do

  else if (idir .eq. 1) then

     facInv = 1.d0 / (facx*facz)

     do n = 1, nc
        do k        = lo(3), hi(3)
           kk       = k * facz
           do j     = lo(2), hi(2)
              jj    = j * facy
              do i  = lo(1), hi(1)
                 ii = i * facx
                 c(i,j,k,n) = 0.d0
                 do    kref = 0, facz-1
                    do iref = 0, facx-1
                       c(i,j,k,n) = c(i,j,k,n) + f(ii+iref,jj,kk+kref,n)
                    end do
                 end do
                 c(i,j,k,n) = c(i,j,k,n) * facInv
              end do
           end do
        end do
     end do

  else

     facInv = 1.d0 / (facx*facy)

     do n = 1, nc
        do k        = lo(3), hi(3)
           kk       = k * facz
           do j     = lo(2), hi(2)
              jj    = j * facy
              do i  = lo(1), hi(1)
                 ii = i * facx
                 c(i,j,k,n) = 0.d0
                 do    jref = 0, facy-1
                    do iref = 0, facx-1
                       c(i,j,k,n) = c(i,j,k,n) + f(ii+iref,jj+jref,kk,n)
                    end do
                 end do
                 c(i,j,k,n) = c(i,j,k,n) * facInv
              end do
           end do
        end do
     end do

  end if

end subroutine bl_avgdown_faces

! ***************************************************************************************
! subroutine bl_avgdown_faces
! ***************************************************************************************

subroutine bl_avgdown_edges (lo, hi, &
     f, f_l1, f_l2, f_l3, f_h1, f_h2, f_h3, &
     c, c_l1, c_l2, c_l3, c_h1, c_h2, c_h3, &
     ratio,idir,nc)

  use amrex_fort_module, only : amrex_real
  implicit none
  integer          :: lo(3),hi(3)
  integer          :: f_l1, f_l2, f_l3, f_h1, f_h2, f_h3
  integer          :: c_l1, c_l2, c_l3, c_h1, c_h2, c_h3
  integer          :: ratio(3), idir, nc
  real(amrex_real) :: f(f_l1:f_h1, f_l2:f_h2, f_l3:f_h3, nc)
  real(amrex_real) :: c(c_l1:c_h1, c_l2:c_h2, c_l3:c_h3, nc)

  ! Local variables
  integer :: i, j, k, n, facx, facy, facz, iref, jref, kref, ii, jj, kk
  real(amrex_real) :: facInv

  facx = ratio(1)
  facy = ratio(2)
  facz = ratio(3)

  if (idir .eq. 0) then

     facInv = 1.d0 / facx

     do n = 1, nc
        do k        = lo(3), hi(3)
           kk       = k * facz
           do j     = lo(2), hi(2)
              jj    = j * facy
              do i  = lo(1), hi(1)
                 ii = i * facx
                 c(i,j,k,n) = 0.d0
                 do iref = 0, facx-1
                    c(i,j,k,n) = c(i,j,k,n) + f(ii+iref,jj,kk,n)
                 end do
                 c(i,j,k,n) = c(i,j,k,n) * facInv
              end do
           end do
        end do
     end do

  else if (idir .eq. 1) then

     facInv = 1.d0 / facy

     do n = 1, nc
        do k        = lo(3), hi(3)
           kk       = k * facz
           do j     = lo(2), hi(2)
              jj    = j * facy
              do i  = lo(1), hi(1)
                 ii = i * facx
                 c(i,j,k,n) = 0.d0
                 do jref = 0, facy-1
                    c(i,j,k,n) = c(i,j,k,n) + f(ii,jj+jref,kk,n)
                 end do
                 c(i,j,k,n) = c(i,j,k,n) * facInv
              end do
           end do
        end do
     end do

  else

     facInv = 1.d0 / facz

     do n = 1, nc
        do k        = lo(3), hi(3)
           kk       = k * facz
           do j     = lo(2), hi(2)
              jj    = j * facy
              do i  = lo(1), hi(1)
                 ii = i * facx
                 c(i,j,k,n) = 0.d0
                 do kref = 0, facz-1
                    c(i,j,k,n) = c(i,j,k,n) + f(ii,jj,kk+kref,n)
                 end do
                 c(i,j,k,n) = c(i,j,k,n) * facInv
              end do
           end do
        end do
     end do

  end if

end subroutine bl_avgdown_edges

! ***************************************************************************************
! subroutine bl_avgdown
! ***************************************************************************************

subroutine bl_avgdown (lo,hi,&
     fine,f_l1,f_l2,f_l3,f_h1,f_h2,f_h3, &
     crse,c_l1,c_l2,c_l3,c_h1,c_h2,c_h3, &
     lrat,ncomp)

  use amrex_fort_module, only : amrex_real
  implicit none
  
  integer f_l1,f_l2,f_l3,f_h1,f_h2,f_h3
  integer c_l1,c_l2,c_l3,c_h1,c_h2,c_h3
  integer lo(3), hi(3)
  integer lrat(3), ncomp
  real(amrex_real) crse(c_l1:c_h1,c_l2:c_h2,c_l3:c_h3,ncomp)
  real(amrex_real) fine(f_l1:f_h1,f_l2:f_h2,f_l3:f_h3,ncomp)
  
  integer :: i, j, k, ii, jj, kk, n, iref, jref, kref
  real(amrex_real) :: volfrac

  volfrac = 1.d0 / dble(lrat(1)*lrat(2)*lrat(3))
  
  do n = 1, ncomp
     do k        = lo(3), hi(3)
        kk       = k * lrat(3)
        do j     = lo(2), hi(2)
           jj    = j * lrat(2)
           do i  = lo(1), hi(1)
              ii = i * lrat(1)
              crse(i,j,k,n) = 0.d0
              do       kref = 0, lrat(3)-1
                 do    jref = 0, lrat(2)-1
                    do iref = 0, lrat(1)-1
                       crse(i,j,k,n) = crse(i,j,k,n) + fine(ii+iref,jj+jref,kk+kref,n)
                    end do
                 end do
              end do
              crse(i,j,k,n) = volfrac * crse(i,j,k,n)
           end do
        end do
     end do
  end do
end subroutine bl_avgdown

subroutine bl_avgdown_nodes (lo,hi,&
     fine,f_l1,f_l2,f_l3,f_h1,f_h2,f_h3, &
     crse,c_l1,c_l2,c_l3,c_h1,c_h2,c_h3, &
     lrat,ncomp)

  use amrex_fort_module, only : amrex_real
  implicit none
  
  integer f_l1,f_l2,f_l3,f_h1,f_h2,f_h3
  integer c_l1,c_l2,c_l3,c_h1,c_h2,c_h3
  integer lo(3), hi(3)
  integer lrat(3), ncomp
  real(amrex_real) crse(c_l1:c_h1,c_l2:c_h2,c_l3:c_h3,ncomp)
  real(amrex_real) fine(f_l1:f_h1,f_l2:f_h2,f_l3:f_h3,ncomp)
  
  integer :: i, j, k, ii, jj, kk, n
  
  do n = 1, ncomp
     do k        = lo(3), hi(3)
        kk       = k * lrat(3)
        do j     = lo(2), hi(2)
           jj    = j * lrat(2)
           do i  = lo(1), hi(1)
              ii = i * lrat(1)
              crse(i,j,k,n) = fine(ii,jj,kk,n)
           end do
        end do
     end do
  end do
end subroutine bl_avgdown_nodes


subroutine amrex_compute_divergence (lo, hi, divu, dlo, dhi, u, ulo, uhi, &
     v, vlo, vhi, w, wlo, whi, dxinv) bind(c)
  use amrex_fort_module, only : amrex_real
  implicit none
  integer, dimension(3), intent(in) :: lo, hi, dlo, dhi, ulo, uhi, vlo, vhi, wlo, whi
  real(amrex_real), intent(inout) :: divu(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3))
  real(amrex_real), intent(in   ) ::    u(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
  real(amrex_real), intent(in   ) ::    v(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
  real(amrex_real), intent(in   ) ::    w(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))
  real(amrex_real), intent(in) :: dxinv(3)
  integer :: i,j,k
  do       k = lo(3), hi(3)
     do    j = lo(2), hi(2)
        do i = lo(1), hi(1)
           divu(i,j,k) = dxinv(1) * (u(i+1,j,k)-u(i,j,k)) &
                +        dxinv(2) * (v(i,j+1,k)-v(i,j,k)) &
                +        dxinv(3) * (w(i,j,k+1)-w(i,j,k))
        end do
     end do
  end do
end subroutine amrex_compute_divergence
