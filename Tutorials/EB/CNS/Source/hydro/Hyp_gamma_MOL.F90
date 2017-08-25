module advection_module
  use amrex_fort_module, only : rt=>amrex_real
  implicit none
  private
  public :: hyp_mol_gam_3d

  logical, parameter :: debug = .false.

contains

  subroutine hyp_mol_gam_3d(q, qd_lo, qd_hi, &
                     lo, hi, dx, flux1, flux2, flux3)

    use mempool_module, only : amrex_allocate, amrex_deallocate
    use cns_module, only : urho, umx, umy, umz, ueden, ueint, utemp, nvar, &
         qrho,qu,qv,qw,qp,qc,qeint,qtemp,qvar, smallp, smallr
    use cns_physics_module, only : gamma
    use slope_module, only : slopex , slopey, slopez
    use riemann_module, only : analriem

    integer, intent(in) :: qd_lo(3), qd_hi(3)
    integer, intent(in) :: lo(3), hi(3)
    real(rt), intent(in) :: dx(3)
    real(rt), intent(in   ) ::     q( qd_lo(1): qd_hi(1), qd_lo(2): qd_hi(2), qd_lo(3): qd_hi(3),QVAR)
    real(rt), intent(  out) :: flux1(lo(1):hi(1)+1,lo(2):hi(2)  ,lo(3):hi(3)  ,5)
    real(rt), intent(  out) :: flux2(lo(1):hi(1)  ,lo(2):hi(2)+1,lo(3):hi(3)  ,5)
    real(rt), intent(  out) :: flux3(lo(1):hi(1)  ,lo(2):hi(2)  ,lo(3):hi(3)+1,5)

    real(rt) :: qtempl(5), qtempr(5), fluxtemp(5)
    real(rt) :: dxinv, dyinv, dzinv, cspeed
    integer :: i, j, k
    integer :: qtlo(3), qthi(3)
    real(rt), pointer, contiguous :: dq(:,:,:,:)

    integer, parameter :: nextra = 0

    qtlo = lo - nextra - 1
    qthi = hi + nextra + 1

    call amrex_allocate ( dq, qtlo(1), qthi(1), qtlo(2), qthi(2), qtlo(3), qthi(3), 1, 5)

    ! Local constants
    dxinv = 1.d0/dx(1)
    dyinv = 1.d0/dx(2)
    dzinv = 1.d0/dx(3)

    call slopex(q,qd_lo,qd_hi, &
         dq,qtlo,qthi, &
         lo(1)-nextra,lo(2)-nextra,lo(3)-nextra,  &
         hi(1)+nextra,hi(2)+nextra,hi(3)+nextra,QVAR)

    do       k = lo(3)-nextra, hi(3)+nextra
       do    j = lo(2)-nextra, hi(2)+nextra
          do i = lo(1)-nextra, hi(1)+nextra+1

!  alphas   1,2,3 correspond to u-c, u, u+c repsectively, 4 and 5 are transverse velocities

! right eigenvectors are rho, u, p, v, w

!   in qtemp, 1 is rho, 2 is u, 3 is p , 4 is v and 5 is w

             cspeed = q(i-1,j,k,QC)
             qtempl(1) = q(i-1,j,k,QRHO) + 0.5d0 * ( (dq(i-1,j,k,1)+dq(i-1,j,k,3))/cspeed + dq(i-1,j,k,2))
             qtempl(2) = q(i-1,j,k,QU) + 0.5d0 * ( (dq(i-1,j,k,3)-dq(i-1,j,k,1))/q(i-1,j,k,QRHO))
             qtempl(3)=  q(i-1,j,k,QP) + 0.5d0 *  (dq(i-1,j,k,1)+dq(i-1,j,k,3))*cspeed 

             qtempl(4) = q(i-1,j,k,QV) + 0.5d0 * dq(i-1,j,k,4)
             qtempl(5) = q(i-1,j,k,Qw) + 0.5d0 * dq(i-1,j,k,5)
             
             cspeed = q(i,j,k,QC)
             qtempr(1) = q(i,j,k,QRHO) - 0.5d0 * ( (dq(i,j,k,1)+dq(i,j,k,3))/cspeed + dq(i,j,k,2))
             qtempr(2) = q(i,j,k,QU) - 0.5d0 * ( (dq(i,j,k,3)-dq(i,j,k,1))/q(i,j,k,QRHO))
             qtempr(3)=  q(i,j,k,QP) - 0.5d0 *  (dq(i,j,k,1)+dq(i,j,k,3))*cspeed 
             qtempr(4) = q(i,j,k,QV) - 0.5d0 * dq(i,j,k,4)
             qtempr(5) = q(i,j,k,Qw) - 0.5d0 *  dq(i,j,k,5)
             
             call analriem(gamma,qtempl, qtempr,  smallp,smallr,fluxtemp,  debug)

             flux1(i,j,k,URHO) = fluxtemp(1)
             flux1(i,j,k,UMX) = fluxtemp(2)
             flux1(i,j,k,UMY) = fluxtemp(4)
             flux1(i,j,k,UMZ) = fluxtemp(5)
             flux1(i,j,k,UEDEN) = fluxtemp(3)
          enddo
       enddo
    enddo

    call slopey(q,qd_lo,qd_hi, &
         dq,qtlo,qthi, &
         lo(1)-nextra,lo(2)-nextra,lo(3)-nextra,  &
         hi(1)+nextra,hi(2)+nextra,hi(3)+nextra,QVAR)

    do       k = lo(3)-nextra, hi(3)+nextra
       do    j = lo(2)-nextra, hi(2)+nextra+1
          do i = lo(1)-nextra, hi(1)+nextra

!     1,2,3 correspond to u-c, u, u+c repsectively

! right eigenvectors are rho, v, p, u, w

!   in qtemp, 1 is rho, 2 is v, 3 is q , 4 is u and 5 is w

             cspeed = q(i,j-1,k,QC)
             qtempl(1) = q(i,j-1,k,QRHO) + 0.5d0 * ( (dq(i,j-1,k,1)+dq(i,j-1,k,3))/cspeed + dq(i,j-1,k,2))
             qtempl(2) = q(i,j-1,k,QV) + 0.5d0 * ( (dq(i,j-1,k,3)-dq(i,j-1,k,1))/q(i,j-1,k,QRHO))
             qtempl(3) = q(i,j-1,k,QP) + 0.5d0 *  (dq(i,j-1,k,1)+dq(i,j-1,k,3))*cspeed 
             qtempl(4) = q(i,j-1,k,QU) + 0.5d0 * dq(i,j-1,k,4)
             qtempl(5) = q(i,j-1,k,Qw) + 0.5d0 * dq(i,j-1,k,5)
             
             cspeed = q(i,j,k,QC)
             qtempr(1) = q(i,j,k,QRHO) - 0.5d0 * ( (dq(i,j,k,1)+dq(i,j,k,3))/cspeed + dq(i,j,k,2))
             qtempr(2) = q(i,j,k,QV) - 0.5d0 * ( (dq(i,j,k,3)-dq(i,j,k,1))/q(i,j,k,QRHO))
             qtempr(3) = q(i,j,k,QP) - 0.5d0 *  (dq(i,j,k,1)+dq(i,j,k,3))*cspeed 
             qtempr(4) = q(i,j,k,QU) - 0.5d0 * dq(i,j,k,4)
             qtempr(5) = q(i,j,k,Qw) - 0.5d0 *  dq(i,j,k,5)

             call analriem(gamma,qtempl, qtempr,  smallp,smallr,fluxtemp,  debug)

             flux2(i,j,k,URHO) = fluxtemp(1)
             flux2(i,j,k,UMX) = fluxtemp(4)
             flux2(i,j,k,UMY) = fluxtemp(2)
             flux2(i,j,k,UMZ) = fluxtemp(5)
             flux2(i,j,k,UEDEN) = fluxtemp(3)
          enddo
       enddo
    enddo
     
    call slopez(q,qd_lo,qd_hi, &
         dq,qtlo,qthi, &
         lo(1)-nextra,lo(2)-nextra,lo(3)-nextra,   &
         hi(1)+nextra,hi(2)+nextra,hi(3)+nextra,QVAR)
    
    do       k = lo(3)-nextra, hi(3)+nextra+1
       do    j = lo(2)-nextra, hi(2)+nextra
          do i = lo(1)-nextra, hi(1)+nextra
             
!     1,2,3 correspond to u-c, u, u+c repsectively

! right eigenvectors are rho, v, p, u, w

!   in qtemp, 1 is rho, 2 is v, 3 is q , 4 is u and 5 is w

             cspeed = q(i,j,k-1,QC)
             qtempl(1) = q(i,j,k-1,QRHO) + 0.5d0 * ( (dq(i,j,k-1,1)+dq(i,j,k-1,3))/cspeed + dq(i,j,k-1,2))
             qtempl(2) = q(i,j,k-1,QW) + 0.5d0 * ( (dq(i,j,k-1,3)-dq(i,j,k-1,1))/q(i,j,k-1,QRHO))
             qtempl(3) = q(i,j,k-1,QP) + 0.5d0 *  (dq(i,j,k-1,1)+dq(i,j,k-1,3))*cspeed 
             qtempl(4) = q(i,j,k-1,QU) + 0.5d0 * dq(i,j,k-1,4)
             qtempl(5) = q(i,j,k-1,QV) + 0.5d0 * dq(i,j,k-1,5)

             cspeed = q(i,j,k,QC)
             qtempr(1) = q(i,j,k,QRHO) - 0.5d0 * ( (dq(i,j,k,1)+dq(i,j,k,3))/cspeed + dq(i,j,k,2))
             qtempr(2) = q(i,j,k,QW) - 0.5d0 * ( (dq(i,j,k,3)-dq(i,j,k,1))/q(i,j,k,QRHO))
             qtempr(3)=  q(i,j,k,QP) - 0.5d0 *  (dq(i,j,k,1)+dq(i,j,k,3))*cspeed 
             qtempr(4) = q(i,j,k,QU) - 0.5d0 * dq(i,j,k,4)
             qtempr(5) = q(i,j,k,QV) - 0.5d0 *  dq(i,j,k,5)

             call analriem(gamma,qtempl, qtempr,  smallp,smallr,fluxtemp,  debug)

             flux3(i,j,k,URHO) = fluxtemp(1)
             flux3(i,j,k,UMX) = fluxtemp(4)
             flux3(i,j,k,UMY) = fluxtemp(5)
             flux3(i,j,k,UMZ) = fluxtemp(2)
             flux3(i,j,k,UEDEN) = fluxtemp(3)
          enddo
       enddo
    enddo

    call amrex_deallocate(dq)

  end subroutine hyp_mol_gam_3d

end module advection_module
