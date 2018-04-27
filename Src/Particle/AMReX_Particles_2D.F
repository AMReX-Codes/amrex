#if 0
! ::
! :: ----------------------------------------------------------
! :: Sum the fine grid mass into the coarse grid.
! :: Overlap is given in coarse grid coordinates.
! :: ----------------------------------------------------------
! ::
      subroutine part_summassdown(crse,c_l1,c_l2,c_h1,c_h2,
     &                            fine,f_l1,f_l2,f_h1,f_h2,
     &                            lo,hi,lrat)
      implicit none

      integer          c_l1,c_l2,c_h1,c_h2
      integer          f_l1,f_l2,f_h1,f_h2
      integer          lo(2), hi(2), lrat(2)
      double precision crse(c_l1:c_h1,c_l2:c_h2)
      double precision fine(f_l1:f_h1,f_l2:f_h2)

      integer          i, j, ic, jc, ioff, joff, lratx, lraty

      lratx = lrat(1)
      lraty = lrat(2)

      do joff = 0, lraty-1
         do jc = lo(2), hi(2)
            j = jc*lraty + joff
            do ioff = 0, lratx-1
               do ic = lo(1), hi(1)
                  i = ic*lratx + ioff
                  crse(ic,jc) = crse(ic,jc) + fine(i,j)
               enddo
            enddo
         enddo
      enddo

      end subroutine part_summassdown
! ::
! :: ----------------------------------------------------------
! :: Sum the crse grid mass into the fine grid.
! :: Overlap is given in coarse grid coordinates.
! :: ----------------------------------------------------------
! ::
      subroutine part_summassup(crse,c_l1,c_l2,c_h1,c_h2,
     &                          fine,f_l1,f_l2,f_h1,f_h2,
     &                          lo,hi,lrat)
      implicit none

      integer          c_l1,c_l2,c_h1,c_h2
      integer          f_l1,f_l2,f_h1,f_h2
      integer          lo(2), hi(2), lrat(2)
      double precision crse(c_l1:c_h1,c_l2:c_h2)
      double precision fine(f_l1:f_h1,f_l2:f_h2)

      integer          i, j, ic, jc, ioff, joff, lratx, lraty
      double precision volfrac

      lratx = lrat(1)
      lraty = lrat(2)

      volfrac = 1.d0/float(lratx*lraty)

      do joff = 0, lraty-1
         do jc = lo(2), hi(2)
            j = jc*lraty + joff
            do ioff = 0, lratx-1
               do ic = lo(1), hi(1)
                  i = ic*lratx + ioff
                  fine(i,j) = fine(i,j) + volfrac*crse(ic,jc)
               enddo
            enddo
         enddo
      enddo

      end subroutine part_summassup

! ::
! :: ----------------------------------------------------------
! :: Sum the crse grid density into the fine grid.
! :: Overlap is given in coarse grid coordinates.
! :: This differs from summassup in that the crse contribution
! :: is *not* weighted by the volume fraction, since this is
! :: density, not mass
! :: ----------------------------------------------------------
! ::
      subroutine part_sumdensup(crse,c_l1,c_l2,c_h1,c_h2,
     &                          fine,f_l1,f_l2,f_h1,f_h2,
     &                          lo,hi,lrat)
      implicit none

      integer          c_l1,c_l2,c_h1,c_h2
      integer          f_l1,f_l2,f_h1,f_h2
      integer          lo(2), hi(2), lrat(2)
      double precision crse(c_l1:c_h1,c_l2:c_h2)
      double precision fine(f_l1:f_h1,f_l2:f_h2)

      integer          i, j, ic, jc, ioff, joff, lratx, lraty

      lratx = lrat(1)
      lraty = lrat(2)

      do joff = 0, lraty-1
         do jc = lo(2), hi(2)
            j = jc*lraty + joff
            do ioff = 0, lratx-1
               do ic = lo(1), hi(1)
                  i = ic*lratx + ioff
                  fine(i,j) = fine(i,j) + crse(ic,jc)
               enddo
            enddo
         enddo
      enddo

      end subroutine part_sumdensup
! ::
! :: ----------------------------------------------------------
! ::
!
!     Bicubic interpolation of cell-centered data.
!
!     Needs two grow cells to work properly.
!
      SUBROUTINE BICUBICINTERP(X,Y,GRID,ILO,IHI,PLO,PHI,DX,Z)

      IMPLICIT DOUBLE PRECISION (A-H,K-L,O-Z)

      INTEGER          ILO(2),IHI(2)
      DOUBLE PRECISION GRID(ILO(1):IHI(1),ILO(2):IHI(2))
      DOUBLE PRECISION DX(2),PLO(2),PHI(2),F(4),F1(4),F2(4),F12(4)
      DOUBLE PRECISION DUMMY1,DUMMY2,DUMMY3,DUMMY4,DUMMY5

      XMIN = PLO(1)
      XMAX = PHI(1)
      YMIN = PLO(2)
      YMAX = PHI(2)

      M  = IHI(1)-ILO(1)
      N  = IHI(2)-ILO(2)
      H  = DX(1)
      XX = ( X - XMIN ) / H
      MM = INT ( XX )
      K  = DX(2)
      YY = ( Y - YMIN ) / K
      NN = INT ( YY )

C     Next line insures all elements of GRID used are in range lo-hi

      MM = MAX ( ILO(1) + 1 , MIN ( IHI(1) - 2, MM ) )
      XR = XX - MM
      NN = MAX ( ILO(2) + 1 , MIN ( IHI(2) - 2, NN ) )
      YR = YY - NN

#if 0
      if ( MM .ne. INT ( XX ) ) then
         print*, 'BICUBIC: MM shifted, extrapolating ...'
      end if
      if ( NN .ne. INT ( YY ) ) then
         print*, 'BICUBIC: NN shifted, extrapolating ...'
      end if
#endif

C     Assemble information needed by BCUINT()
 
      XL = MM * H + XMIN
      XU = XL + H
      YL = NN * K + YMIN
      YU = YL + K
 
      DO J = NN, NN+1
         DO  I = MM, MM+1
            II = 1 + (I-MM) + 3 * (J-NN)
            IF ( II .EQ. 5 ) II = 3
            F(II)   =  GRID(I,J)
            F1(II)  = (GRID(I+1,J) - GRID(I-1,J)) / (2D0 * H)
            F2(II)  = (GRID(I,J+1) - GRID(I,J-1)) / (2D0 * K)
            F12(II) = (GRID(I+1,J+1) - GRID(I+1,J-1) - GRID(I-1,J+1)
     &           +  GRID(I-1,J-1)) / ( 4D0 * H * K )
         END DO
      END DO
 
      CALL BCUINT( F, F1, F2, F12, XL, XU, YL, YU, X, Y, Z,
     &             DUMMY1, DUMMY2, DUMMY3, DUMMY4, DUMMY5)

      END
 
      SUBROUTINE BCUINT(Y,Y1,Y2,Y12,X1L,X1U,X2L,X2U,X1,X2,A,A1,A2,
     &                  A11,A22,A12)
 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 
      DIMENSION Y(4),Y1(4),Y2(4),Y12(4),C(4,4)
 
      CALL BCUCOF(Y,Y1,Y2,Y12,X1U-X1L,X2U-X2L,C)
 
      T   = (X1 - X1L) / (X1U - X1L)
      U   = (X2 - X2L) / (X2U - X2L)
      A   = 0D0
      A1  = 0D0
      A2  = 0D0
      A11 = 0D0
      A22 = 0D0
 
      DO  I = 4, 1, -1
         A   = T * A  + ((C(I,4)*U + C(I,3))*U + C(I,2))*U + C(I,1)
         A2  = T * A2 + (3D0*C(I,4)*U + 2D0*C(I,3))*U + C(I,2)
         A22 = T * A22 + 6D0*C(I,4)*U + 2D0*C(I,3)
         A1  = U * A1 + (3D0*C(4,I)*T + 2D0*C(3,I))*T + C(2,I)
         A11 = U * A11 + 6D0*C(4,I)*T + 2D0*C(3,I)
      END DO
 
      A12 =          C(2,2) + U*(2D0*C(2,3) + 3D0*C(2,4)*U) +
     &      T * (2D0*C(3,2) + U*(4D0*C(3,3) + 6D0*C(3,4)*U) +
     &      T * (3D0*C(4,2) + U*(6D0*C(4,3) + 9D0*C(4,4)*U)))
 
      A1  = A1  /  (X1U - X1L)
      A11 = A11 / ((X1U - X1L) * (X1U - X1L))
      A2  = A2  /  (X2U - X2L)
      A22 = A22 / ((X2U - X2L) * (X2U - X2L))
      A12 = A12 / ((X1U - X1L) * (X2U - X2L))

      END
 
      SUBROUTINE BCUCOF(Y,Y1,Y2,Y12,D1,D2,C)
 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 
      DIMENSION Y(4),Y1(4),Y2(4),Y12(4),C(4,4),CL(16),X(16),WT(16,16)
 
      DATA WT/1,0,-3,2,4*0,-3,0,9,-6,2,0,-6,4,8*0,3,0,-9,6,-2,0,6,-4,
     &  10*0,9,-6,2*0,-6,4,2*0,3,-2,6*0,-9,6,2*0,6,-4,
     &  4*0,1,0,-3,2,-2,0,6,-4,1,0,-3,2,8*0,-1,0,3,-2,1,0,-3,2,
     &  10*0,-3,2,2*0,3,-2,6*0,3,-2,2*0,-6,4,2*0,3,-2,
     &  0,1,-2,1,5*0,-3,6,-3,0,2,-4,2,9*0,3,-6,3,0,-2,4,-2,
     &  10*0,-3,3,2*0,2,-2,2*0,-1,1,6*0,3,-3,2*0,-2,2,
     &  5*0,1,-2,1,0,-2,4,-2,0,1,-2,1,9*0,-1,2,-1,0,1,-2,1,
     &  10*0,1,-1,2*0,-1,1,6*0,-1,1,2*0,2,-2,2*0,-1,1/
 
      D1D2 = D1 * D2
 
      DO I = 1, 4
         X(I)    = Y(I)
         X(I+4)  = Y1(I)  * D1
         X(I+8)  = Y2(I)  * D2
         X(I+12) = Y12(I) * D1D2
      END DO
 
      DO I = 1, 16
         XX = 0D0
         DO K = 1, 16
            XX = XX + WT(I,K) * X(K)
          END DO
         CL(I) = XX
      END DO
 
      L = 0
 
      DO I = 1, 4
         DO J = 1, 4
            L = L + 1
            C(I,J) = CL(L)
         END DO
      END DO

      END

#endif
