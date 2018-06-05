#include "AMReX_LO_BCTYPES.H"


module amrex_habec_module

  ! habec is Hypre abec, where abec is the form of the linear equation
  ! we are solving:
  ! 
  ! alpha*phi - div(beta*grad phi) + div(\vec{c}*phi) 

  use amrex_fort_module, only : rt => amrex_real
  implicit none

contains

subroutine amrex_hbvec(vec, &
                 reg_l1,reg_l2,reg_l3,reg_h1,reg_h2,reg_h3, &
                 cdir, bct, bho, bcl, &
                 bcval, bcv_l1,bcv_l2,bcv_l3,bcv_h1,bcv_h2,bcv_h3, &
                 mask, msk_l1,msk_l2,msk_l3,msk_h1,msk_h2,msk_h3, &
                 b, bbox_l1,bbox_l2,bbox_l3,bbox_h1,bbox_h2,bbox_h3, &
                 beta, dx) bind(C, name="amrex_hbvec")

  use amrex_fort_module, only : rt => amrex_real
  integer :: reg_l1,reg_l2,reg_l3,reg_h1,reg_h2,reg_h3
  integer :: bcv_l1,bcv_l2,bcv_l3,bcv_h1,bcv_h2,bcv_h3
  integer :: msk_l1,msk_l2,msk_l3,msk_h1,msk_h2,msk_h3
  integer :: bbox_l1,bbox_l2,bbox_l3,bbox_h1,bbox_h2,bbox_h3
  integer :: cdir, bct, bho
  real(rt)         :: bcl, beta, dx(3)
  real(rt)         :: vec(reg_l1:reg_h1,reg_l2:reg_h2,reg_l3:reg_h3)
  real(rt)         :: bcval(bcv_l1:bcv_h1,bcv_l2:bcv_h2,bcv_l3:bcv_h3)
  integer :: mask(msk_l1:msk_h1,msk_l2:msk_h2,msk_l3:msk_h3)
  real(rt)         :: b(bbox_l1:bbox_h1,bbox_l2:bbox_h2,bbox_l3:bbox_h3)
  real(rt)         :: h, bfv
  real(rt)         :: h2, th2
  integer :: i, j, k
  if (cdir == 0 .OR. cdir == 3) then
     h = dx(1)
  elseif (cdir == 1 .OR. cdir == 4) then
     h = dx(2)
  else
     h = dx(3)
  endif
  if (bct == LO_DIRICHLET) then
     if (bho >= 1) then
        h2 = 0.5e0_rt * h
        th2 = 3.e0_rt * h2
        bfv = 2.e0_rt * beta / ((bcl + h2) * (bcl + th2))
     else
        bfv = (beta / h) / (0.5e0_rt * h + bcl)
     endif
  else if (bct == LO_NEUMANN) then
     bfv = beta / h
  else
     print *, "hbvec: unsupported boundary type"
     stop
  endif
  if (cdir == 0) then
     i = reg_l1
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           if (mask(i-1,j,k) > 0) then
              vec(i,j,k) = vec(i,j,k) + &
                   bfv * b(i,j,k) * bcval(i-1,j,k)
           endif
        enddo
     enddo
  else if (cdir == 3) then
     i = reg_h1
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           if (mask(i+1,j,k) > 0) then
              vec(i,j,k) = vec(i,j,k) + &
                   bfv * b(i+1,j,k) * bcval(i+1,j,k)
           endif
        enddo
     enddo
  else if (cdir == 1) then
     j = reg_l2
     do k = reg_l3, reg_h3
        do i = reg_l1, reg_h1
           if (mask(i,j-1,k) > 0) then
              vec(i,j,k) = vec(i,j,k) + &
                   bfv * b(i,j,k) * bcval(i,j-1,k)
           endif
        enddo
     enddo
  else if (cdir == 4) then
     j = reg_h2
     do k = reg_l3, reg_h3
        do i = reg_l1, reg_h1
           if (mask(i,j+1,k) > 0) then
              vec(i,j,k) = vec(i,j,k) + &
                   bfv * b(i,j+1,k) * bcval(i,j+1,k)
           endif
        enddo
     enddo
  else if (cdir == 2) then
     k = reg_l3
     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1
           if (mask(i,j,k-1) > 0) then
              vec(i,j,k) = vec(i,j,k) + &
                   bfv * b(i,j,k) * bcval(i,j,k-1)
           endif
        enddo
     enddo
  else if (cdir == 5) then
     k = reg_h3
     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1
           if (mask(i,j,k+1) > 0) then
              vec(i,j,k) = vec(i,j,k) + &
                   bfv * b(i,j,k+1) * bcval(i,j,k+1)
           endif
        enddo
     enddo
  else
     print *, "hbvec: impossible face orientation"
  endif
end subroutine amrex_hbvec

subroutine amrex_hbvec3(vec, &
                  reg_l1,reg_l2,reg_l3,reg_h1,reg_h2,reg_h3, &
                  cdir, bctype, bho, bcl, &
                  bcval, bcv_l1,bcv_l2,bcv_l3,bcv_h1,bcv_h2,bcv_h3, &
                  mask, msk_l1,msk_l2,msk_l3,msk_h1,msk_h2,msk_h3, &
                  b, bbox_l1,bbox_l2,bbox_l3,bbox_h1,bbox_h2,bbox_h3, &
                  beta, dx) bind(C, name="amrex_hbvec3")

  use amrex_fort_module, only : rt => amrex_real
  integer :: reg_l1,reg_l2,reg_l3,reg_h1,reg_h2,reg_h3
  integer :: bcv_l1,bcv_l2,bcv_l3,bcv_h1,bcv_h2,bcv_h3
  integer :: msk_l1,msk_l2,msk_l3,msk_h1,msk_h2,msk_h3
  integer :: bbox_l1,bbox_l2,bbox_l3,bbox_h1,bbox_h2,bbox_h3
  integer :: cdir, bctype, bho
  real(rt)         :: bcl, beta, dx(3)
  real(rt)         :: vec(reg_l1:reg_h1,reg_l2:reg_h2,reg_l3:reg_h3)
  real(rt)         :: bcval(bcv_l1:bcv_h1,bcv_l2:bcv_h2,bcv_l3:bcv_h3)
  integer :: mask(msk_l1:msk_h1,msk_l2:msk_h2,msk_l3:msk_h3)
  real(rt)         :: b(bbox_l1:bbox_h1,bbox_l2:bbox_h2,bbox_l3:bbox_h3)
  real(rt)         :: h, bfv
  real(rt)         :: h2, th2
  integer :: i, j, k, bct
  if (cdir == 0 .OR. cdir == 3) then
     h = dx(1)
  elseif (cdir == 1 .OR. cdir == 4) then
     h = dx(2)
  else
     h = dx(3)
  endif
  bct = bctype
  if (cdir == 0) then
     i = reg_l1
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           if (mask(i-1,j,k) > 0) then
              if (bct == LO_DIRICHLET) then
                 if (bho >= 1) then
                    h2 = 0.5e0_rt * h
                    th2 = 3.e0_rt * h2
                    bfv = 2.e0_rt * beta / ((bcl + h2) * (bcl + th2))
                 else
                    bfv = (beta / h) / (0.5e0_rt * h + bcl)
                 endif
                 bfv = bfv * b(i,j,k)
              else if (bct == LO_NEUMANN) then
                 bfv = beta / h
              else
                 print *, "hbvec3: unsupported boundary type"
                 stop
              endif
              vec(i,j,k) = vec(i,j,k) + bfv * bcval(i-1,j,k)
           endif
        enddo
     enddo
  else if (cdir == 3) then
     i = reg_h1
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           if (mask(i+1,j,k) > 0) then
              if (bct == LO_DIRICHLET) then
                 if (bho >= 1) then
                    h2 = 0.5e0_rt * h
                    th2 = 3.e0_rt * h2
                    bfv = 2.e0_rt * beta / ((bcl + h2) * (bcl + th2))
                 else
                    bfv = (beta / h) / (0.5e0_rt * h + bcl)
                 endif
                 bfv = bfv * b(i+1,j,k)
              else if (bct == LO_NEUMANN) then
                 bfv = beta / h
              else
                 print *, "hbvec3: unsupported boundary type"
                 stop
              endif
              vec(i,j,k) = vec(i,j,k) + bfv * bcval(i+1,j,k)
           endif
        enddo
     enddo
  else if (cdir == 1) then
     j = reg_l2
     do k = reg_l3, reg_h3
        do i = reg_l1, reg_h1
           if (mask(i,j-1,k) > 0) then
              if (bct == LO_DIRICHLET) then
                 if (bho >= 1) then
                    h2 = 0.5e0_rt * h
                    th2 = 3.e0_rt * h2
                    bfv = 2.e0_rt * beta / ((bcl + h2) * (bcl + th2))
                 else
                    bfv = (beta / h) / (0.5e0_rt * h + bcl)
                 endif
                 bfv = bfv * b(i,j,k)
              else if (bct == LO_NEUMANN) then
                 bfv = beta / h
              else
                 print *, "hbvec3: unsupported boundary type"
                 stop
              endif
              vec(i,j,k) = vec(i,j,k) + bfv * bcval(i,j-1,k)
           endif
        enddo
     enddo
  else if (cdir == 4) then
     j = reg_h2
     do k = reg_l3, reg_h3
        do i = reg_l1, reg_h1
           if (mask(i,j+1,k) > 0) then
              if (bct == LO_DIRICHLET) then
                 if (bho >= 1) then
                    h2 = 0.5e0_rt * h
                    th2 = 3.e0_rt * h2
                    bfv = 2.e0_rt * beta / ((bcl + h2) * (bcl + th2))
                 else
                    bfv = (beta / h) / (0.5e0_rt * h + bcl)
                 endif
                 bfv = bfv * b(i,j+1,k)
              else if (bct == LO_NEUMANN) then
                 bfv = beta / h
              else
                 print *, "hbvec3: unsupported boundary type"
                 stop
              endif
              vec(i,j,k) = vec(i,j,k) + bfv * bcval(i,j+1,k)
           endif
        enddo
     enddo
  else if (cdir == 2) then
     k = reg_l3
     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1
           if (mask(i,j,k-1) > 0) then
              if (bct == LO_DIRICHLET) then
                 if (bho >= 1) then
                    h2 = 0.5e0_rt * h
                    th2 = 3.e0_rt * h2
                    bfv = 2.e0_rt * beta / ((bcl + h2) * (bcl + th2))
                 else
                    bfv = (beta / h) / (0.5e0_rt * h + bcl)
                 endif
                 bfv = bfv * b(i,j,k)
              else if (bct == LO_NEUMANN) then
                 bfv = beta / h
              else
                 print *, "hbvec3: unsupported boundary type"
                 stop
              endif
              vec(i,j,k) = vec(i,j,k) + bfv * bcval(i,j,k-1)
           endif
        enddo
     enddo
  else if (cdir == 5) then
     k = reg_h3
     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1
           if (mask(i,j,k+1) > 0) then
              if (bct == LO_DIRICHLET) then
                 if (bho >= 1) then
                    h2 = 0.5e0_rt * h
                    th2 = 3.e0_rt * h2
                    bfv = 2.e0_rt * beta / ((bcl + h2) * (bcl + th2))
                 else
                    bfv = (beta / h) / (0.5e0_rt * h + bcl)
                 endif
                 bfv = bfv * b(i,j,k+1)
              else if (bct == LO_NEUMANN) then
                 bfv = beta / h
              else
                 print *, "hbvec3: unsupported boundary type"
                 stop
              endif
              vec(i,j,k) = vec(i,j,k) + bfv * bcval(i,j,k+1)
           endif
        enddo
     enddo
  else
     print *, "hbvec3: impossible face orientation"
  endif
end subroutine amrex_hbvec3


subroutine amrex_hmac(mat, a, &
                abox_l1,abox_l2,abox_l3,abox_h1,abox_h2,abox_h3, &
                reg_l1,reg_l2,reg_l3,reg_h1,reg_h2,reg_h3, &
                alpha) bind(C, name="amrex_hmac")

  use amrex_fort_module, only : rt => amrex_real
  integer :: abox_l1,abox_l2,abox_l3,abox_h1,abox_h2,abox_h3
  integer :: reg_l1,reg_l2,reg_l3,reg_h1,reg_h2,reg_h3
  real(rt)         :: a(abox_l1:abox_h1,abox_l2:abox_h2,abox_l3:abox_h3)
  real(rt)         :: mat(0:6, reg_l1:reg_h1,reg_l2:reg_h2,reg_l3:reg_h3)
  real(rt)         :: alpha
  integer :: i, j, k
  if (alpha == 0.e0_rt) then
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           do i = reg_l1, reg_h1
              mat(0,i,j,k) = 0.e0_rt
           enddo
        enddo
     enddo
  else
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           do i = reg_l1, reg_h1
              mat(0,i,j,k) = alpha * a(i,j,k)
           enddo
        enddo
     enddo
  endif
end subroutine amrex_hmac

subroutine amrex_hmbc(mat, b, &
                bbox_l1,bbox_l2,bbox_l3,bbox_h1,bbox_h2,bbox_h3, &
                reg_l1,reg_l2,reg_l3,reg_h1,reg_h2,reg_h3, &
                beta, dx, n) bind(C, name="amrex_hmbc")

  use amrex_fort_module, only : rt => amrex_real
  integer :: bbox_l1,bbox_l2,bbox_l3,bbox_h1,bbox_h2,bbox_h3
  integer :: reg_l1,reg_l2,reg_l3,reg_h1,reg_h2,reg_h3
  integer :: n
  real(rt)         :: b(bbox_l1:bbox_h1,bbox_l2:bbox_h2,bbox_l3:bbox_h3)
  real(rt)         :: mat(0:6, reg_l1:reg_h1,reg_l2:reg_h2,reg_l3:reg_h3)
  real(rt)         :: beta, dx(3)
  real(rt)         :: fac
  integer :: i, j, k
  if (n == 0) then
     fac = beta / (dx(1)**2)
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           do i = reg_l1, reg_h1
              mat(0,i,j,k) = mat(0,i,j,k) + fac * (b(i,j,k) + b(i+1,j,k))
              mat(1,i,j,k) = - fac * b(i,j,k)
              mat(2,i,j,k) = - fac * b(i+1,j,k)
           enddo
        enddo
     enddo
  elseif (n == 1) then
     fac = beta / (dx(2)**2)
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           do i = reg_l1, reg_h1
              mat(0,i,j,k) = mat(0,i,j,k) + fac * (b(i,j,k) + b(i,j+1,k))
              mat(3,i,j,k) = - fac * b(i,j,k)
              mat(4,i,j,k) = - fac * b(i,j+1,k)
           enddo
        enddo
     enddo
  else
     fac = beta / (dx(3)**2)
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           do i = reg_l1, reg_h1
              mat(0,i,j,k) = mat(0,i,j,k) + fac * (b(i,j,k) + b(i,j,k+1))
              mat(5,i,j,k) = - fac * b(i,j,k)
              mat(6,i,j,k) = - fac * b(i,j,k+1)
           enddo
        enddo
     enddo
  endif
end subroutine amrex_hmbc


subroutine amrex_hmmat(mat, &
                 reg_l1,reg_l2,reg_l3,reg_h1,reg_h2,reg_h3, &
                 cdir, bct, bho, bcl, &
                 mask, msk_l1,msk_l2,msk_l3,msk_h1,msk_h2,msk_h3, &
                 b, bbox_l1,bbox_l2,bbox_l3,bbox_h1,bbox_h2,bbox_h3, &
                 beta, dx) bind(C, name="amrex_hmmat")

  use amrex_fort_module, only : rt => amrex_real
  integer :: reg_l1,reg_l2,reg_l3,reg_h1,reg_h2,reg_h3
  integer :: msk_l1,msk_l2,msk_l3,msk_h1,msk_h2,msk_h3
  integer :: bbox_l1,bbox_l2,bbox_l3,bbox_h1,bbox_h2,bbox_h3
  integer :: cdir, bct, bho
  real(rt)         :: bcl, beta, dx(3)
  real(rt)         :: mat(0:6, reg_l1:reg_h1,reg_l2:reg_h2,reg_l3:reg_h3)
  integer :: mask(msk_l1:msk_h1,msk_l2:msk_h2,msk_l3:msk_h3)
  real(rt)         :: b(bbox_l1:bbox_h1,bbox_l2:bbox_h2,bbox_l3:bbox_h3)
  real(rt)         :: h, fac, bfm, bfv
  real(rt)         :: bfm2, h2, th2
  integer :: i, j, k
  if (cdir == 0 .OR. cdir == 3) then
     h = dx(1)
  elseif (cdir == 1 .OR. cdir == 4) then
     h = dx(2)
  else
     h = dx(3)
  endif
  fac = beta / (h**2)
  if (bct == LO_DIRICHLET) then
     if (bho >= 1) then
        h2 = 0.5e0_rt * h
        th2 = 3.e0_rt * h2
        bfm = fac * (th2 - bcl) / (bcl + h2) - fac
        bfm2 = fac * (bcl - h2) / (bcl + th2)
     else
        bfv = (beta / h) / (0.5e0_rt * h + bcl)
        bfm = bfv - fac
     endif
  else if (bct == LO_NEUMANN) then
     bfm = -fac
     bfm2 = 0.e0_rt
  else
     print *, "hmmat: unsupported boundary type"
     stop
  endif
  if (cdir == 0) then
     i = reg_l1
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           if (mask(i-1,j,k) > 0) then
              mat(0,i,j,k) = mat(0,i,j,k) + bfm * b(i,j,k)
              mat(1,i,j,k) = 0.e0_rt
              if (bho >= 1) then
                 mat(2,i,j,k) = mat(2,i,j,k) + bfm2 * b(i,j,k)
              endif
           endif
        enddo
     enddo
  else if (cdir == 3) then
     i = reg_h1
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           if (mask(i+1,j,k) > 0) then
              mat(0,i,j,k) = mat(0,i,j,k) + bfm * b(i+1,j,k)
              mat(2,i,j,k) = 0.e0_rt
              if (bho >= 1) then
                 mat(1,i,j,k) = mat(1,i,j,k) + bfm2 * b(i+1,j,k)
              endif
           endif
        enddo
     enddo
  else if (cdir == 1) then
     j = reg_l2
     do k = reg_l3, reg_h3
        do i = reg_l1, reg_h1
           if (mask(i,j-1,k) > 0) then
              mat(0,i,j,k) = mat(0,i,j,k) + bfm * b(i,j,k)
              mat(3,i,j,k) = 0.e0_rt
              if (bho >= 1) then
                 mat(4,i,j,k) = mat(4,i,j,k) + bfm2 * b(i,j,k)
              endif
           endif
        enddo
     enddo
  else if (cdir == 4) then
     j = reg_h2
     do k = reg_l3, reg_h3
        do i = reg_l1, reg_h1
           if (mask(i,j+1,k) > 0) then
              mat(0,i,j,k) = mat(0,i,j,k) + bfm * b(i,j+1,k)
              mat(4,i,j,k) = 0.e0_rt
              if (bho >= 1) then
                 mat(3,i,j,k) = mat(3,i,j,k) + bfm2 * b(i,j+1,k)
              endif
           endif
        enddo
     enddo
  else if (cdir == 2) then
     k = reg_l3
     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1
           if (mask(i,j,k-1) > 0) then
              mat(0,i,j,k) = mat(0,i,j,k) + bfm * b(i,j,k)
              mat(5,i,j,k) = 0.e0_rt
              if (bho >= 1) then
                 mat(6,i,j,k) = mat(6,i,j,k) + bfm2 * b(i,j,k)
              endif
           endif
        enddo
     enddo
  else if (cdir == 5) then
     k = reg_h3
     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1
           if (mask(i,j,k+1) > 0) then
              mat(0,i,j,k) = mat(0,i,j,k) + bfm * b(i,j,k+1)
              mat(6,i,j,k) = 0.e0_rt
              if (bho >= 1) then
                 mat(5,i,j,k) = mat(5,i,j,k) + bfm2 * b(i,j,k+1)
              endif
           endif
        enddo
     enddo
  else
     print *, "hmmat: impossible face orientation"
  endif
end subroutine amrex_hmmat

subroutine amrex_hmmat3(mat, &
                  reg_l1,reg_l2,reg_l3,reg_h1,reg_h2,reg_h3, &
                  cdir, bctype, bho, bcl, &
                  mask, msk_l1,msk_l2,msk_l3,msk_h1,msk_h2,msk_h3, &
                  b, bbox_l1,bbox_l2,bbox_l3,bbox_h1,bbox_h2,bbox_h3, &
                  beta, dx) bind(C, name="amrex_hmmat3")

  use amrex_fort_module, only : rt => amrex_real
  integer :: reg_l1,reg_l2,reg_l3,reg_h1,reg_h2,reg_h3
  integer :: msk_l1,msk_l2,msk_l3,msk_h1,msk_h2,msk_h3
  integer :: bbox_l1,bbox_l2,bbox_l3,bbox_h1,bbox_h2,bbox_h3
  integer :: cdir, bctype, bho
  real(rt)         :: bcl, beta, dx(3)
  real(rt)         :: mat(0:6, reg_l1:reg_h1,reg_l2:reg_h2,reg_l3:reg_h3)
  integer :: mask(msk_l1:msk_h1,msk_l2:msk_h2,msk_l3:msk_h3)
  real(rt)         :: b(bbox_l1:bbox_h1,bbox_l2:bbox_h2,bbox_l3:bbox_h3)
  real(rt)         :: h, fac, bfm, bfv
  real(rt)         :: bfm2, h2, th2
  integer :: i, j, k, bct
  ! The -fac * b(i,j,k) term applied to the matrix diagonal is the contribution
  ! from the interior stencil which must be removed at the boundary.
  if (cdir == 0 .OR. cdir == 3) then
     h = dx(1)
  elseif (cdir == 1 .OR. cdir == 4) then
     h = dx(2)
  else
     h = dx(3)
  endif
  bct = bctype
  fac = beta / (h**2)
  if (cdir == 0) then
     i = reg_l1
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           if (mask(i-1,j,k) > 0) then
              if (bct == LO_DIRICHLET) then
                 if (bho >= 1) then
                    h2 = 0.5e0_rt * h
                    th2 = 3.e0_rt * h2
                    bfm = fac * (th2 - bcl) / (bcl + h2)  * b(i,j,k)
                    bfm2 = fac * (bcl - h2) / (bcl + th2) * b(i,j,k)
                 else
                    bfv = (beta / h) / (0.5e0_rt * h + bcl)
                    bfm = bfv * b(i,j,k)
                 endif
              else if (bct == LO_NEUMANN) then
                 bfm  = 0.e0_rt
                 bfm2 = 0.e0_rt
              else
                 print *, "hmmat3: unsupported boundary type"
                 stop
              endif
              mat(0,i,j,k) = mat(0,i,j,k) + bfm - fac * b(i,j,k)
              mat(1,i,j,k) = 0.e0_rt
              if (bho >= 1) then
                 mat(2,i,j,k) = mat(2,i,j,k) + bfm2
              endif
           endif
        enddo
     enddo
  else if (cdir == 3) then
     i = reg_h1
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           if (mask(i+1,j,k) > 0) then
              if (bct == LO_DIRICHLET) then
                 if (bho >= 1) then
                    h2 = 0.5e0_rt * h
                    th2 = 3.e0_rt * h2
                    bfm = fac * (th2 - bcl) / (bcl + h2)  * b(i+1,j,k)
                    bfm2 = fac * (bcl - h2) / (bcl + th2) * b(i+1,j,k)
                 else
                    bfv = (beta / h) / (0.5e0_rt * h + bcl)
                    bfm = bfv * b(i+1,j,k)
                 endif
              else if (bct == LO_NEUMANN) then
                 bfm  = 0.e0_rt
                 bfm2 = 0.e0_rt
              else
                 print *, "hmmat3: unsupported boundary type"
                 stop
              endif
              mat(0,i,j,k) = mat(0,i,j,k) + bfm - fac * b(i+1,j,k)
              mat(2,i,j,k) = 0.e0_rt
              if (bho >= 1) then
                 mat(1,i,j,k) = mat(1,i,j,k) + bfm2
              endif
           endif
        enddo
     enddo
  else if (cdir == 1) then
     j = reg_l2
     do k = reg_l3, reg_h3
        do i = reg_l1, reg_h1
           if (mask(i,j-1,k) > 0) then
              if (bct == LO_DIRICHLET) then
                 if (bho >= 1) then
                    h2 = 0.5e0_rt * h
                    th2 = 3.e0_rt * h2
                    bfm = fac * (th2 - bcl) / (bcl + h2)  * b(i,j,k)
                    bfm2 = fac * (bcl - h2) / (bcl + th2) * b(i,j,k)
                 else
                    bfv = (beta / h) / (0.5e0_rt * h + bcl)
                    bfm = bfv * b(i,j,k)
                 endif
              else if (bct == LO_NEUMANN) then
                 bfm  = 0.e0_rt
                 bfm2 = 0.e0_rt
              else
                 print *, "hmmat3: unsupported boundary type"
                 stop
              endif
              mat(0,i,j,k) = mat(0,i,j,k) + bfm - fac * b(i,j,k)
              mat(3,i,j,k) = 0.e0_rt
              if (bho >= 1) then
                 mat(4,i,j,k) = mat(4,i,j,k) + bfm2
              endif
           endif
        enddo
     enddo
  else if (cdir == 4) then
     j = reg_h2
     do k = reg_l3, reg_h3
        do i = reg_l1, reg_h1
           if (mask(i,j+1,k) > 0) then
              if (bct == LO_DIRICHLET) then
                 if (bho >= 1) then
                    h2 = 0.5e0_rt * h
                    th2 = 3.e0_rt * h2
                    bfm = fac * (th2 - bcl) / (bcl + h2)  * b(i,j+1,k)
                    bfm2 = fac * (bcl - h2) / (bcl + th2) * b(i,j+1,k)
                 else
                    bfv = (beta / h) / (0.5e0_rt * h + bcl)
                    bfm = bfv * b(i,j+1,k)
                 endif
              else if (bct == LO_NEUMANN) then
                 bfm  = 0.e0_rt
                 bfm2 = 0.e0_rt
              else
                 print *, "hmmat3: unsupported boundary type"
                 stop
              endif
              mat(0,i,j,k) = mat(0,i,j,k) + bfm - fac * b(i,j+1,k)
              mat(4,i,j,k) = 0.e0_rt
              if (bho >= 1) then
                 mat(3,i,j,k) = mat(3,i,j,k) + bfm2
              endif
           endif
        enddo
     enddo
  else if (cdir == 2) then
     k = reg_l3
     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1
           if (mask(i,j,k-1) > 0) then
              if (bct == LO_DIRICHLET) then
                 if (bho >= 1) then
                    h2 = 0.5e0_rt * h
                    th2 = 3.e0_rt * h2
                    bfm = fac * (th2 - bcl) / (bcl + h2)  * b(i,j,k)
                    bfm2 = fac * (bcl - h2) / (bcl + th2) * b(i,j,k)
                 else
                    bfv = (beta / h) / (0.5e0_rt * h + bcl)
                    bfm = bfv * b(i,j,k)
                 endif
              else if (bct == LO_NEUMANN) then
                 bfm  = 0.e0_rt
                 bfm2 = 0.e0_rt
              else
                 print *, "hmmat3: unsupported boundary type"
                 stop
              endif
              mat(0,i,j,k) = mat(0,i,j,k) + bfm - fac * b(i,j,k)
              mat(5,i,j,k) = 0.e0_rt
              if (bho >= 1) then
                 mat(6,i,j,k) = mat(6,i,j,k) + bfm2
              endif
           endif
        enddo
     enddo
  else if (cdir == 5) then
     k = reg_h3
     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1
           if (mask(i,j,k+1) > 0) then
              if (bct == LO_DIRICHLET) then
                 if (bho >= 1) then
                    h2 = 0.5e0_rt * h
                    th2 = 3.e0_rt * h2
                    bfm = fac * (th2 - bcl) / (bcl + h2)  * b(i,j,k+1)
                    bfm2 = fac * (bcl - h2) / (bcl + th2) * b(i,j,k+1)
                 else
                    bfv = (beta / h) / (0.5e0_rt * h + bcl)
                    bfm = bfv * b(i,j,k+1)
                 endif
              else if (bct == LO_NEUMANN) then
                 bfm  = 0.e0_rt
                 bfm2 = 0.e0_rt
              else
                 print *, "hmmat3: unsupported boundary type"
                 stop
              endif
              mat(0,i,j,k) = mat(0,i,j,k) + bfm - fac * b(i,j,k+1)
              mat(6,i,j,k) = 0.e0_rt
              if (bho >= 1) then
                 mat(5,i,j,k) = mat(5,i,j,k) + bfm2
              endif
           endif
        enddo
     enddo
  else
     print *, "hmmat3: impossible face orientation"
  endif
end subroutine amrex_hmmat3

end module amrex_habec_module
