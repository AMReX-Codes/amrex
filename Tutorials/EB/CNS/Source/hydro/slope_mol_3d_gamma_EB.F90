module ebslope_module

  use amrex_ebcellflag_module, only : get_neighbor_cells, is_regular_cell
  use amrex_fort_module, only : rt=>amrex_real
  use cns_module, only : qrho,qu,qv,qw,qp,qc,qeint,qtemp,qvar

  implicit none

  private

  public ebslopex, ebslopey, ebslopez

  integer, parameter :: plm_iorder = 2
  real(rt), parameter :: plm_theta = 2.0   ! [1,2]; 1: minmod; 2: van Leer's MC

contains

  subroutine ebslopex(q,qd_lo,qd_hi, &
       dqxal,qpd_lo,qpd_hi, &
       flag, fglo, fghi, &
       ilo1,ilo2,ilo3,ihi1,ihi2,ihi3,nv)

    integer, intent(in) :: qd_lo(3), qd_hi(3)
    integer, intent(in) :: qpd_lo(3),qpd_hi(3)
    integer, intent(in) :: fglo(3), fghi(3)
    integer, intent(in) :: ilo1, ilo2, ihi1, ihi2, ilo3, ihi3, nv
    real(rt), intent(in) :: q(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),nv)
    real(rt), intent(out) :: dqxal(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),5)
    integer, intent(in) :: flag(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3))

    integer i, j, k, n
    real(rt) :: slop, dsgn,dlim,dcen
    real(rt) :: dlft(ilo1-1:ihi1+1,5), drgt(ilo1-1:ihi1+1,5)
    integer :: nbr(-1:1,-1:1,-1:1)

    if(plm_iorder.eq.1) then

       do n = 1, 5
          do k = ilo3, ihi3
            do j = ilo2, ihi2
               do i = ilo1-1, ihi1+1
                  dqxal(i,j,k,n) = 0.d0
               enddo
            enddo
         enddo
      enddo
      
   else

      ! Compute slopes in first coordinate direction
       do k = ilo3, ihi3
          do j = ilo2, ihi2
             do i = ilo1-1, ihi1+1

                call get_neighbor_cells(flag(i,j,k), nbr)

                if (is_regular_cell(flag(i,j,k)) .or. nbr(-1,0,0).eq.1) then
                   dlft(i,1) = 0.5d0*(q(i,j,k,QP)-q(i-1,j,k,QP))/q(i,j,k,QC) - 0.5d0*q(i,j,k,QRHO)*(q(i,j,k,QU) - q(i-1,j,k,QU))
                   dlft(i,2) = (q(i,j,k,QRHO)-q(i-1,j,k,QRHO))- (q(i,j,k,QP) - q(i-1,j,k,QP))/q(i,j,k,QC)**2
                   dlft(i,3) = 0.5d0*(q(i,j,k,QP)-q(i-1,j,k,QP))/q(i,j,k,QC) + 0.5d0*q(i,j,k,QRHO)*(q(i,j,k,QU) - q(i-1,j,k,QU))
                   dlft(i,4) = q(i,j,k,QV) - q(i-1,j,k,QV) 
                   dlft(i,5) = q(i,j,k,QW) - q(i-1,j,k,QW) 
                else
                   dlft(i,:) = 0.d0
                end if

                if (is_regular_cell(flag(i,j,k)) .or. nbr(1,0,0).eq.1) then
                   drgt(i,1) = 0.5d0*(q(i+1,j,k,QP)-q(i,j,k,QP))/q(i,j,k,QC) - 0.5d0*q(i,j,k,QRHO)*(q(i+1,j,k,QU) - q(i,j,k,QU))
                   drgt(i,2) = (q(i+1,j,k,QRHO)-q(i,j,k,QRHO))- (q(i+1,j,k,QP) - q(i,j,k,QP))/q(i,j,k,QC)**2
                   drgt(i,3) = 0.5d0*(q(i+1,j,k,QP)-q(i,j,k,QP))/q(i,j,k,QC) + 0.5d0*q(i,j,k,QRHO)*(q(i+1,j,k,QU) - q(i,j,k,QU))
                   drgt(i,4) = q(i+1,j,k,QV) - q(i,j,k,QV) 
                   drgt(i,5) = q(i+1,j,k,QW) - q(i,j,k,QW) 
                else
                   drgt(i,:) = 0.d0
                end if

             enddo

             do n=1,5
                ! First compute Fromm slopes
                do i = ilo1-1, ihi1+1
                   dcen = 0.5d0 * (dlft(i,n)+drgt(i,n))
                   dsgn = sign(1.d0, dcen)
                   slop = plm_theta * min( abs(dlft(i,n)), abs(drgt(i,n)) )
                   if (dlft(i,n)*drgt(i,n) .ge. 0.d0) then
                      dlim = slop
                   else
                      dlim = 0.d0
                   endif
                   dqxal(i,j,k,n) = dsgn*min( dlim, abs(dcen) )
                enddo
             end do
          end do
       end do
    end if

  end subroutine ebslopex

  subroutine ebslopey(q,qd_lo,qd_hi, &
       dqyal,qpd_lo,qpd_hi, &
       flag, fglo, fghi, &
       ilo1,ilo2,ilo3,ihi1,ihi2,ihi3,nv)

    integer, intent(in) :: qd_lo(3), qd_hi(3)
    integer, intent(in) :: qpd_lo(3),qpd_hi(3)
    integer, intent(in) :: fglo(3), fghi(3)
    integer, intent(in) :: ilo1, ilo2, ihi1, ihi2, ilo3, ihi3, nv
    real(rt), intent(in) :: q(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),nv)
    real(rt), intent(out) :: dqyal(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),5)
    integer, intent(in) :: flag(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3))
    
    integer i, j, k, n
    real(rt) slop, dsgn,dlim,dcen
    real(rt) :: dlft(ilo1:ihi1,5), drgt(ilo1:ihi1,5)
    integer :: nbr(-1:1,-1:1,-1:1)
    
    if(plm_iorder.eq.1) then
       
       do n = 1, 5
          do k = ilo3, ihi3
             do j = ilo2-1, ihi2+1
                do i = ilo1, ihi1
                   dqyal(i,j,k,n) = 0.d0
                enddo
             enddo
          enddo
       enddo
       
    else
       
       ! Compute slopes in first coordinate direction
       do k = ilo3, ihi3
          do j = ilo2-1, ihi2+1
             do i = ilo1, ihi1
                call get_neighbor_cells( flag(i,j,k), nbr )

                if (is_regular_cell(flag(i,j,k)) .or. nbr(0,-1,0).eq.1) then
                   dlft(i,1) = 0.5d0*(q(i,j,k,QP)-q(i,j-1,k,QP))/q(i,j,k,QC) - 0.5d0*q(i,j,k,QRHO)*(q(i,j,k,QV) - q(i,j-1,k,QV))
                   dlft(i,2) = (q(i,j,k,QRHO)-q(i,j-1,k,QRHO))- (q(i,j,k,QP) - q(i,j-1,k,QP))/q(i,j,k,QC)**2
                   dlft(i,3) = 0.5d0*(q(i,j,k,QP)-q(i,j-1,k,QP))/q(i,j,k,QC) + 0.5d0*q(i,j,k,QRHO)*(q(i,j,k,QV) - q(i,j-1,k,QV))
                   dlft(i,4) = q(i,j,k,QU) - q(i,j-1,k,QU) 
                   dlft(i,5) = q(i,j,k,QW) - q(i,j-1,k,QW) 
                else
                   dlft(i,:) = 0.d0
                end if

                if (is_regular_cell(flag(i,j,k)) .or. nbr(0,1,0).eq.1) then
                   drgt(i,1) = 0.5d0*(q(i,j+1,k,QP)-q(i,j,k,QP))/q(i,j,k,QC) - 0.5d0*q(i,j,k,QRHO)*(q(i,j+1,k,QV) - q(i,j,k,QV))
                   drgt(i,2) = (q(i,j+1,k,QRHO)-q(i,j,k,QRHO))- (q(i,j+1,k,QP) - q(i,j,k,QP))/q(i,j,k,QC)**2
                   drgt(i,3) = 0.5d0*(q(i,j+1,k,QP)-q(i,j,k,QP))/q(i,j,k,QC) + 0.5d0*q(i,j,k,QRHO)*(q(i,j+1,k,QV) - q(i,j,k,QV))
                   drgt(i,4) = q(i,j+1,k,QU) - q(i,j,k,QU) 
                   drgt(i,5) = q(i,j+1,k,QW) - q(i,j,k,QW)
                else
                   drgt(i,:) = 0.d0
                end if
             enddo
             
             do n=1,5
                ! First compute Fromm slopes
                do i = ilo1, ihi1
                   dcen = 0.5d0 * (dlft(i,n)+drgt(i,n))
                   dsgn = sign(1.d0, dcen)
                   slop = plm_theta * min( abs(dlft(i,n)), abs(drgt(i,n)) )
                   if (dlft(i,n)*drgt(i,n) .ge. 0.d0) then
                      dlim = slop
                   else
                      dlim = 0.d0
                   endif
                   dqyal(i,j,k,n) = dsgn*min( dlim, abs(dcen) )
                enddo
             enddo
          enddo
       enddo
    end if
    
  end subroutine ebslopey
  
  subroutine ebslopez(q,qd_lo,qd_hi, &
       dqzal,qpd_lo,qpd_hi, &
       flag, fglo, fghi, &
       ilo1,ilo2,ilo3,ihi1,ihi2,ihi3,nv)
    
    integer, intent(in) :: qd_lo(3), qd_hi(3)
    integer, intent(in) :: qpd_lo(3),qpd_hi(3)
    integer, intent(in) :: fglo(3), fghi(3)
    integer, intent(in) :: ilo1, ilo2, ihi1, ihi2, ilo3, ihi3, nv
    real(rt), intent(in) :: q(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),nv)
    real(rt), intent(out) :: dqzal(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),5)
    integer, intent(in) :: flag(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3))
    
    integer i, j, k, n
    real(rt) slop,dsgn,dlim,dcen
    real(rt) :: dlft(ilo1:ihi1,nv), drgt(ilo1:ihi1,nv)
    integer :: nbr(-1:1,-1:1,-1:1)
    
    if(plm_iorder.eq.1) then
       
       do n = 1, 5
          do k = ilo3-1, ihi3+1
             do j = ilo2, ihi2
                do i = ilo1, ihi1
                   dqzal(i,j,k,n) = 0.d0
                enddo
             enddo
          enddo
       enddo

    else

       do k = ilo3-1, ihi3+1
            do j = ilo2, ihi2
               do i = ilo1, ihi1
                  call get_neighbor_cells( flag(i,j,k), nbr )

                  if (is_regular_cell(flag(i,j,k)) .or. nbr(0,0,-1).eq.1) then
                     dlft(i,1) = 0.5d0*(q(i,j,k,QP)-q(i,j,k-1,QP))/q(i,j,k,QC) - 0.5d0*q(i,j,k,QRHO)*(q(i,j,k,QW) - q(i,j,k-1,QW))
                     dlft(i,2) = (q(i,j,k,QRHO)-q(i,j,k-1,QRHO))- (q(i,j,k,QP) - q(i,j,k-1,QP))/q(i,j,k,QC)**2
                     dlft(i,3) = 0.5d0*(q(i,j,k,QP)-q(i,j,k-1,QP))/q(i,j,k,QC) + 0.5d0*q(i,j,k,QRHO)*(q(i,j,k,QW) - q(i,j,k-1,QW))
                     dlft(i,4) = q(i,j,k,QU) - q(i,j,k-1,QU) 
                     dlft(i,5) = q(i,j,k,QV) - q(i,j,k-1,QV) 
                  else
                     dlft(i,:) = 0.d0
                  end if

                  if (is_regular_cell(flag(i,j,k)) .or. nbr(0,0,1).eq.1) then
                     drgt(i,1) = 0.5d0*(q(i,j,k+1,QP)-q(i,j,k,QP))/q(i,j,k,QC) - 0.5d0*q(i,j,k,QRHO)*(q(i,j,k+1,QW) - q(i,j,k,QW))
                     drgt(i,2) = (q(i,j,k+1,QRHO)-q(i,j,k,QRHO))- (q(i,j,k+1,QP) - q(i,j,k,QP))/q(i,j,k,QC)**2
                     drgt(i,3) = 0.5d0*(q(i,j,k+1,QP)-q(i,j,k,QP))/q(i,j,k,QC) + 0.5d0*q(i,j,k,QRHO)*(q(i,j,k+1,QW) - q(i,j,k,QW))
                     drgt(i,4) = q(i,j,k+1,QU) - q(i,j,k,QU) 
                     drgt(i,5) = q(i,j,k+1,QV) - q(i,j,k,QV) 
                  else
                     drgt(i,:) = 0.d0
                  end if
               enddo

               do n=1,5
                  do i = ilo1, ihi1
                     dcen = 0.5d0 * (dlft(i,n)+drgt(i,n))
                     dsgn = sign(1.d0, dcen)
                     slop = plm_theta * min( abs(dlft(i,n)), abs(drgt(i,n)) )
                     if (dlft(i,n)*drgt(i,n) .ge. 0.d0) then
                        dlim = slop
                     else
                        dlim = 0.d0
                     endif
                     dqzal(i,j,k,n) = dsgn*min( dlim, abs(dcen) )
                  enddo
               enddo
            enddo
         enddo
      endif
      
    end subroutine ebslopez

end module ebslope_module
