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
