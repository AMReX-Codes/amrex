
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_BC_TYPES.H"
#include "AMReX_ArrayLim.H"

#define SDIM 3

c :: ----------------------------------------------------------
c :: MACDIV:  compute the MAC divergence in generalized coordinates
c ::
c :: INPUTS / OUTPUTS:
c ::  dmac        <=  MAC divergence (cell centered)
c ::  DIMS(dmac)   => index limits for dmac
c ::  lo,hi        => index limits of grid interior
c ::  ux           => X edge velocity
c ::  DIMS(ux)     => index limits for ux
c ::  uy           => Y edge velocity
c ::  DIMS(uy)     => index limits for uy
c ::  uz           => Y edge velocity
c ::  DIMS(uz)     => index limits for uz
c ::  xarea        => area of cell faces in X dircetion
c ::  DIMS(ax)     => index limits for xarea
c ::  yarea        => area of cell faces in Y dircetion
c ::  DIMS(ay)     => index limits for yarea
c ::  zarea        => area of cell faces in Z dircetion
c ::  DIMS(az)     => index limits for zarea
c ::  vol          => cell volume
c ::  DIMS(vol)    => index limits for vol
c :: ----------------------------------------------------------
c ::
       subroutine FORT_MACDIV (dmac,DIMS(dmac),lo,hi,
     &                         ux,DIMS(ux),uy,DIMS(uy),uz,DIMS(uz),
     &                         xarea,DIMS(ax),yarea,DIMS(ay),
     &                         zarea,DIMS(az),vol,DIMS(vol))

       implicit none
       integer DIMDEC(dmac)
       integer lo(SDIM), hi(SDIM)
       integer DIMDEC(ux)
       integer DIMDEC(uy)
       integer DIMDEC(uz)
       integer DIMDEC(ax)
       integer DIMDEC(ay)
       integer DIMDEC(az)
       integer DIMDEC(vol)
       REAL_T  dmac(DIMV(dmac))
       REAL_T    ux(DIMV(ux))
       REAL_T    uy(DIMV(uy))
       REAL_T    uz(DIMV(uz))
       REAL_T xarea(DIMV(ax))
       REAL_T yarea(DIMV(ay))
       REAL_T zarea(DIMV(az))
       REAL_T   vol(DIMV(vol))

       integer i, j, k

!$omp parallel do private(i,j,k)
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                dmac(i,j,k) = (xarea(i+1,j,k)*ux(i+1,j,k) - xarea(i,j,k)*ux(i,j,k)
     &               +     yarea(i,j+1,k)*uy(i,j+1,k) - yarea(i,j,k)*uy(i,j,k)
     &               +     zarea(i,j,k+1)*uz(i,j,k+1) - zarea(i,j,k)*uz(i,j,k)
     &               )/vol(i,j,k)
             end do
          end do
       end do
!$omp end parallel do

       end

c :: ----------------------------------------------------------
c :: SCALEAREA
c ::          area = area * anel_coeff
c ::                 OR 
c ::          area = area / anel_coeff
c :: ----------------------------------------------------------

       subroutine FORT_SCALEAREA (xarea,DIMS(ax),yarea,DIMS(ay),
     &                            zarea,DIMS(az),
     &                            anel_coeff,anel_coeff_lo,anel_coeff_hi,lo,hi,mult)

       implicit none
       integer lo(SDIM), hi(SDIM)
       integer anel_coeff_lo,anel_coeff_hi
       integer DIMDEC(ax)
       integer DIMDEC(ay)
       integer DIMDEC(az)
       REAL_T xarea(DIMV(ax))
       REAL_T yarea(DIMV(ay))
       REAL_T zarea(DIMV(az))
       REAL_T anel_coeff(anel_coeff_lo:anel_coeff_hi)
       integer mult

       integer i,j,k

       if (mult .eq. 1) then

          do k = lo(3), hi(3)
          do j = lo(2), hi(2)
          do i = lo(1), hi(1)+1
             xarea(i,j,k) =  xarea(i,j,k) * anel_coeff(k)
          end do
          end do
          end do

          do k = lo(3), hi(3)
          do j = lo(2), hi(2)+1
          do i = lo(1), hi(1)
             yarea(i,j,k) =  yarea(i,j,k) * anel_coeff(k)
          end do
          end do
          end do

          do k = lo(3), hi(3)+1
          do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             zarea(i,j,k) =  zarea(i,j,k) * 0.5d0 * (anel_coeff(k)+anel_coeff(k-1))
          end do
          end do
          end do

          if (lo(3) .eq. 0) then
             k = lo(3)

             do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                zarea(i,j,k) =  zarea(i,j,k) * anel_coeff(k-1) /
     $                          (0.5d0 * (anel_coeff(k)+anel_coeff(k-1)))
             end do
             end do

          end if

       else if (mult .eq. -1) then

          do k = lo(3), hi(3)
          do j = lo(2), hi(2)
          do i = lo(1), hi(1)+1
             xarea(i,j,k) =  xarea(i,j,k) / anel_coeff(k)
          end do
          end do
          end do

          do k = lo(3), hi(3)
          do j = lo(2), hi(2)+1
          do i = lo(1), hi(1)
             yarea(i,j,k) =  yarea(i,j,k) / anel_coeff(k)
          end do
          end do
          end do

          do k = lo(3), hi(3)+1
          do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             zarea(i,j,k) =  zarea(i,j,k) / (0.5d0 * (anel_coeff(k)+anel_coeff(k-1)))
          end do
          end do
          end do

          if (lo(3) .eq. 0) then
             k = lo(3)

             do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                zarea(i,j,k) =  zarea(i,j,k) / anel_coeff(k-1) *
     $                         (0.5d0 * (anel_coeff(k)+anel_coeff(k-1)))
             end do
             end do

          end if

       else 
 
          print *,'BOGUS MULT IN SCALEAREA '
          stop

       end if

       end
