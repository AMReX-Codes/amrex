module eb_advection_module

  use amrex_fort_module, only : rt=>amrex_real
  implicit none
  private
  public :: hyp_mol_gam_eb_3d

  integer, parameter, public :: nextra_eb = 2

contains

  subroutine hyp_mol_gam_eb_3d(q, qd_lo, qd_hi, &
                     lo, hi, dx, &
                     flux1, fd1_lo, fd1_hi, &
                     flux2, fd2_lo, fd2_hi, &
                     flux3, fd3_lo, fd3_hi, &
                     flag, fg_lo, fg_hi)

    use amrex_mempool_module, only : amrex_allocate, amrex_deallocate
    use cns_module, only : urho, umx, umy, umz, ueden, ueint, utemp, nvar, &
         qrho,qu,qv,qw,qp,qc,qeint,qtemp,qvar, smallp, smallr
    use cns_physics_module, only : gamma
    use ebslope_module, only : ebslopex , ebslopey, ebslopez
    use riemann_module, only : analriem

    integer, intent(in) :: qd_lo(3), qd_hi(3)
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: fd1_lo(3), fd1_hi(3)
    integer, intent(in) :: fd2_lo(3), fd2_hi(3)
    integer, intent(in) :: fd3_lo(3), fd3_hi(3)
    integer, intent(in) :: fg_lo(3), fg_hi(3)
    real(rt), intent(in) :: dx(3)
    real(rt), intent(in   ) ::     q( qd_lo(1): qd_hi(1), qd_lo(2): qd_hi(2), qd_lo(3): qd_hi(3),QVAR)
    real(rt), intent(inout) :: flux1(fd1_lo(1):fd1_hi(1),fd1_lo(2):fd1_hi(2),fd1_lo(3):fd1_hi(3),5)
    real(rt), intent(inout) :: flux2(fd2_lo(1):fd2_hi(1),fd2_lo(2):fd2_hi(2),fd2_lo(3):fd2_hi(3),5)
    real(rt), intent(inout) :: flux3(fd3_lo(1):fd3_hi(1),fd3_lo(2):fd3_hi(2),fd3_lo(3):fd3_hi(3),5)
    integer,  intent(in   ) ::  flag( fg_lo(1): fg_hi(1), fg_lo(2): fg_hi(2), fg_lo(3): fg_hi(3))

    real(rt) :: cspeed
    integer :: i, j, k
    integer :: qtlo(3), qthi(3)
    real(rt), pointer, contiguous :: dq(:,:,:,:)
    real(rt), dimension(lo(1)-nextra_eb-1:hi(1)+nextra_eb+1) :: rl, ul, pl, ut1l, ut2l
    real(rt), dimension(lo(1)-nextra_eb-1:hi(1)+nextra_eb+1) :: rr, ur, pr, ut1r, ut2r
    integer, parameter :: nextra = nextra_eb

    qtlo = lo - nextra - 1
    qthi = hi + nextra + 1

    call amrex_allocate ( dq, qtlo(1), qthi(1), qtlo(2), qthi(2), qtlo(3), qthi(3), 1, 5)

    call ebslopex(q,qd_lo,qd_hi, &
         dq,qtlo,qthi, flag, fg_lo, fg_hi, &
         lo(1)-nextra, lo(2)-nextra-1, lo(3)-nextra-1,  &
         hi(1)+nextra, hi(2)+nextra+1, hi(3)+nextra+1, QVAR)

    do       k = lo(3)-nextra-1, hi(3)+nextra+1
       do    j = lo(2)-nextra-1, hi(2)+nextra+1
          do i = lo(1)-nextra  , hi(1)+nextra+1

             cspeed = q(i-1,j,k,QC)
             rl(i) = q(i-1,j,k,QRHO) + 0.5d0 * ( (dq(i-1,j,k,1)+dq(i-1,j,k,3))/cspeed + dq(i-1,j,k,2))
             rl(i) = max(rl(i),smallr)
             ul(i) = q(i-1,j,k,QU) + 0.5d0 * ( (dq(i-1,j,k,3)-dq(i-1,j,k,1))/q(i-1,j,k,QRHO))
             pl(i)=  q(i-1,j,k,QP) + 0.5d0 *  (dq(i-1,j,k,1)+dq(i-1,j,k,3))*cspeed 
             pl(i) = max(pl(i),smallp)
             ut1l(i) = q(i-1,j,k,QV) + 0.5d0 * dq(i-1,j,k,4)
             ut2l(i) = q(i-1,j,k,Qw) + 0.5d0 * dq(i-1,j,k,5)
             
             cspeed = q(i,j,k,QC)
             rr(i) = q(i,j,k,QRHO) - 0.5d0 * ( (dq(i,j,k,1)+dq(i,j,k,3))/cspeed + dq(i,j,k,2))
             rr(i) = max(rr(i),smallr)
             ur(i) = q(i,j,k,QU) - 0.5d0 * ( (dq(i,j,k,3)-dq(i,j,k,1))/q(i,j,k,QRHO))
             pr(i)=  q(i,j,k,QP) - 0.5d0 *  (dq(i,j,k,1)+dq(i,j,k,3))*cspeed 
             pr(i) = max(pr(i),smallp)
             ut1r(i) = q(i,j,k,QV) - 0.5d0 * dq(i,j,k,4)
             ut2r(i) = q(i,j,k,Qw) - 0.5d0 *  dq(i,j,k,5)
          end do

          call analriem(gamma, smallp, smallr, lo(1)-nextra, hi(1)+nextra+1, j, k, &
               rl(lo(1)-nextra:), ul(lo(1)-nextra:), pl(lo(1)-nextra:), ut1l(lo(1)-nextra:), ut2l(lo(1)-nextra:), &
               rr(lo(1)-nextra:), ur(lo(1)-nextra:), pr(lo(1)-nextra:), ut1r(lo(1)-nextra:), ut2r(lo(1)-nextra:), &
               flux1, fd1_lo, fd1_hi, 2, 3, 4)
       enddo
    enddo

    call ebslopey(q,qd_lo,qd_hi, &
         dq,qtlo,qthi, flag, fg_lo, fg_hi, &
         lo(1)-nextra-1, lo(2)-nextra, lo(3)-nextra-1,  &
         hi(1)+nextra+1, hi(2)+nextra, hi(3)+nextra+1, QVAR)

    do       k = lo(3)-nextra-1, hi(3)+nextra+1
       do    j = lo(2)-nextra  , hi(2)+nextra+1
          do i = lo(1)-nextra-1, hi(1)+nextra+1

             cspeed = q(i,j-1,k,QC)
             rl(i) = q(i,j-1,k,QRHO) + 0.5d0 * ( (dq(i,j-1,k,1)+dq(i,j-1,k,3))/cspeed + dq(i,j-1,k,2))
             rl(i) = max(rl(i),smallr)
             ul(i) = q(i,j-1,k,QV) + 0.5d0 * ( (dq(i,j-1,k,3)-dq(i,j-1,k,1))/q(i,j-1,k,QRHO))
             pl(i) = q(i,j-1,k,QP) + 0.5d0 *  (dq(i,j-1,k,1)+dq(i,j-1,k,3))*cspeed 
             pl(i) = max(pl(i),smallp)
             ut1l(i) = q(i,j-1,k,QU) + 0.5d0 * dq(i,j-1,k,4)
             ut2l(i) = q(i,j-1,k,Qw) + 0.5d0 * dq(i,j-1,k,5)
             
             cspeed = q(i,j,k,QC)
             rr(i) = q(i,j,k,QRHO) - 0.5d0 * ( (dq(i,j,k,1)+dq(i,j,k,3))/cspeed + dq(i,j,k,2))
             rr(i) = max(rr(i),smallr)
             ur(i) = q(i,j,k,QV) - 0.5d0 * ( (dq(i,j,k,3)-dq(i,j,k,1))/q(i,j,k,QRHO))
             pr(i) = q(i,j,k,QP) - 0.5d0 *  (dq(i,j,k,1)+dq(i,j,k,3))*cspeed 
             pr(i) = max(pr(i),smallp)
             ut1r(i) = q(i,j,k,QU) - 0.5d0 * dq(i,j,k,4)
             ut2r(i) = q(i,j,k,Qw) - 0.5d0 *  dq(i,j,k,5)
          end do

          call analriem(gamma, smallp, smallr, lo(1)-nextra-1, hi(1)+nextra+1, j, k, &
               rl, ul, pl, ut1l, ut2l, &
               rr, ur, pr, ut1r, ut2r, &
               flux2, fd2_lo, fd2_hi, 3, 2, 4)
       enddo
    enddo
    
    call ebslopez(q,qd_lo,qd_hi, &
         dq,qtlo,qthi, flag, fg_lo, fg_hi, &
         lo(1)-nextra-1, lo(2)-nextra-1, lo(3)-nextra,   &
         hi(1)+nextra+1, hi(2)+nextra+1, hi(3)+nextra, QVAR)
    
    do       k = lo(3)-nextra  , hi(3)+nextra+1
       do    j = lo(2)-nextra-1, hi(2)+nextra+1
          do i = lo(1)-nextra-1, hi(1)+nextra+1

             cspeed = q(i,j,k-1,QC)
             rl(i) = q(i,j,k-1,QRHO) + 0.5d0 * ( (dq(i,j,k-1,1)+dq(i,j,k-1,3))/cspeed + dq(i,j,k-1,2))
             rl(i) = max(rl(i),smallr)
             ul(i) = q(i,j,k-1,QW) + 0.5d0 * ( (dq(i,j,k-1,3)-dq(i,j,k-1,1))/q(i,j,k-1,QRHO))
             pl(i)=  q(i,j,k-1,QP) + 0.5d0 *  (dq(i,j,k-1,1)+dq(i,j,k-1,3))*cspeed 
             pl(i) = max(pl(i),smallp)
             ut1l(i) = q(i,j,k-1,QU) + 0.5d0 * dq(i,j,k-1,4)
             ut2l(i) = q(i,j,k-1,QV) + 0.5d0 * dq(i,j,k-1,5)

             cspeed = q(i,j,k,QC)
             rr(i) = q(i,j,k,QRHO) - 0.5d0 * ( (dq(i,j,k,1)+dq(i,j,k,3))/cspeed + dq(i,j,k,2))
             rr(i) = max(rr(i),smallr)
             ur(i) = q(i,j,k,QW) - 0.5d0 * ( (dq(i,j,k,3)-dq(i,j,k,1))/q(i,j,k,QRHO))
             pr(i)=  q(i,j,k,QP) - 0.5d0 *  (dq(i,j,k,1)+dq(i,j,k,3))*cspeed 
             pr(i) = max(pr(i),smallp)
             ut1r(i) = q(i,j,k,QU) - 0.5d0 * dq(i,j,k,4)
             ut2r(i) = q(i,j,k,QV) - 0.5d0 *  dq(i,j,k,5)
          end do

          call analriem(gamma, smallp, smallr, lo(1)-nextra-1, hi(1)+nextra+1, j, k, &
               rl, ul, pl, ut1l, ut2l, &
               rr, ur, pr, ut1r, ut2r, &
               flux3, fd3_lo, fd3_hi, 4, 2, 3)
       enddo
    enddo

    call amrex_deallocate(dq)

  end subroutine hyp_mol_gam_eb_3d

end module eb_advection_module
