/* -------------------------------------------------------------------- */
/*  HyptermOnly_F.F                                                     */
/* -------------------------------------------------------------------- */
#include <REAL.H>
#include <ArrayLim.H>
#include <CONSTANTS.H>
#include <HyptermOnly_F.H>
#include <BC_TYPES.H>
      
#define SDIM 3


      REAL_T function parallel_wtime()
        REAL_T tws
	call WALL_SECOND(tws)
        parallel_wtime = tws
      return
      end


      subroutine FORT_INITDATA(cons, DIMS(cons), q, dx, nvars)

      implicit none

      integer DIMDEC(cons)
      integer DIMDEC(grid)
      integer nvars
      REAL_T  dx(SDIM)
      REAL_T  cons(DIMV(cons), nvars)
      REAL_T  q(DIMV(cons), nvars+1)

      REAL_T  scale, xloc, yloc, zloc, rholoc, uvel, vvel, wvel, eloc, rhoinv
      integer irho, imx, imy, imz, iene, qu, qv, qw, qpres
      integer lo(SDIM),hi(SDIM)
      integer i, j, k

      REAL_T GAMMA, CV, CVinv

      irho = 1
      imx  = 2
      imy  = 3
      imz  = 4
      iene = 5

      lo(1) = cons_l1
      lo(2) = cons_l2
      lo(3) = cons_l3
      hi(1) = cons_h1
      hi(2) = cons_h2
      hi(3) = cons_h3

      scale = 1.0d0


!$OMP PARALLEL DO PRIVATE(i,j,k,xloc,yloc,zloc,
!$OMP&uvel,vvel,wvel,rholoc,eloc)
      do k=lo(3),hi(3)
         zloc = dfloat(k)*dx(3)/scale
         do j=lo(2),hi(2)
            yloc = dfloat(j)*dx(2)/scale
            do i=lo(1),hi(1)
               xloc = dfloat(i)*dx(1)/scale

               uvel   = 1.1d4*sin(1*xloc)*sin(2*yloc)*sin(3*zloc)
               vvel   = 1.0d4*sin(2*xloc)*sin(4*yloc)*sin(1*zloc)
               wvel   = 1.2d4*sin(3*xloc)*cos(2*yloc)*sin(2*zloc)
               rholoc = 1.0d-3 + 1.0d-5*sin(1*xloc)*cos(2*yloc)*cos(3*zloc)
               eloc   = 2.5d9  + 1.0d-3*sin(2*xloc)*cos(2*yloc)*sin(2*zloc)

               cons(i,j,k,irho) = rholoc
               cons(i,j,k,imx)  = rholoc*uvel
               cons(i,j,k,imy)  = rholoc*vvel
               cons(i,j,k,imz)  = rholoc*wvel
               cons(i,j,k,iene) = rholoc*(eloc + (uvel**2+vvel**2+wvel**2)/2)

            enddo
         enddo
      enddo
!$OMP END PARALLEL DO


      GAMMA = 1.4d0
      CV    = 8.3333333333d6
      CVinv = 1.0d0 / CV

!$OMP PARALLEL DO PRIVATE(i,j,k,eloc,rhoinv)
      do k = lo(3),hi(3)
        do j = lo(2),hi(2)
          do i = lo(1),hi(1)

             rhoinv     = 1.0d0/cons(i,j,k,1)
             q(i,j,k,1) = cons(i,j,k,1)
             q(i,j,k,2) = cons(i,j,k,2)*rhoinv
             q(i,j,k,3) = cons(i,j,k,3)*rhoinv
             q(i,j,k,4) = cons(i,j,k,4)*rhoinv

             eloc = cons(i,j,k,5)*rhoinv-0.5d0*(q(i,j,k,2)**2+q(i,j,k,3)**2+q(i,j,k,4)**2)

             q(i,j,k,5) = (GAMMA-1.d0)*eloc*cons(i,j,k,1)
             q(i,j,k,6) = eloc * CVinv

          enddo
        enddo
      enddo
!$OMP END PARALLEL DO


      end




      subroutine FORT_HYPTERM(cons, DIMS(cons), DIMS(grid), q, flux, dx, nvars)

      implicit none

      integer DIMDEC(cons)
      integer DIMDEC(grid)
      integer nvars
      REAL_T  dx(SDIM)
      REAL_T  cons(DIMV(cons), nvars)
      REAL_T  q(DIMV(cons), nvars+1)
      REAL_T  flux(DIMV(grid), nvars)

      REAL_T unp1,unp2,unp3,unp4,unm1,unm2,unm3,unm4
      REAL_T dxinvX, dxinvY, dxinvZ

      REAL_T L1_start, L1_end
      REAL_T L2_start, L2_end
      REAL_T L3_start, L3_end

      integer irho, imx, imy, imz, iene, qu, qv, qw, qpres
      REAL_T ALP, BET, GAM, DEL

      REAL_T t, parallel_wtime

      integer lo(SDIM),hi(SDIM)
      integer JBlocks, JBlockSize, jb
      integer L1iters, L2iters, L3iters
      integer i, j, k



      logical bIOP

!        write (6,*) grid_l1,grid_h1
!        write (6,*) grid_l2,grid_h2
!        write (6,*) grid_l3,grid_h3
!        write (6,*) cons_l1,cons_h1
!        write (6,*) cons_l2,cons_h2
!        write (6,*) cons_l3,cons_h3


      lo(1) = grid_l1
      hi(1) = grid_h1
      lo(2) = grid_l2
      hi(2) = grid_h2
      lo(3) = grid_l3
      hi(3) = grid_h3


      irho = 1
      imx  = 2
      imy  = 3
      imz  = 4
      iene = 5

      qu    = 2
      qv    = 3
      qw    = 4
      qpres = 5

      ALP =  0.8d0
      BET = -0.2d0
      GAM =  4.d0/105.d0
      DEL = -1.d0/280.d0


      L1iters = 0
      L2iters = 0
      L3iters = 0
      bIOP = .true.


      dxinvX = 1.0d0 / dx(1)
      dxinvY = 1.0d0 / dx(2)
      dxinvZ = 1.0d0 / dx(3)


!     if ( bIOP ) then
!      print*, "lo(1) hi(1) = ", lo(1), hi(1)
!      print*, "lo(2) hi(2) = ", lo(2), hi(2)
!      print*, "lo(3) hi(3) = ", lo(3), hi(3)
!     end if


      L1_start = parallel_wtime()

! OMP PARALLEL DO PRIVATE(i,j,k,unp1,unp2,unp3,unp4,
! OMP&unm1,unm2,unm3,unm4) reduction(+ : L1iters)
!$OMP PARALLEL DO PRIVATE(i,j,k,unp1,unp2,unp3,unp4,
!$OMP&unm1,unm2,unm3,unm4)  reduction(+ : L1iters)
      do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             L1iters = L1iters + 1

             unm4 = q(i-4,j,k,qu)
             unm3 = q(i-3,j,k,qu)
             unm2 = q(i-2,j,k,qu)
             unm1 = q(i-1,j,k,qu)

             unp1 = q(i+1,j,k,qu)
             unp2 = q(i+2,j,k,qu)
             unp3 = q(i+3,j,k,qu)
             unp4 = q(i+4,j,k,qu)

             flux(i,j,k,irho)= -dxinvX *
     $             (ALP*(cons(i+1,j,k,imx)-cons(i-1,j,k,imx)) 
     $            + BET*(cons(i+2,j,k,imx)-cons(i-2,j,k,imx)) 
     $            + GAM*(cons(i+3,j,k,imx)-cons(i-3,j,k,imx)) 
     $            + DEL*(cons(i+4,j,k,imx)-cons(i-4,j,k,imx)))

             flux(i,j,k,imx)= -dxinvX *
     $             (ALP*(cons(i+1,j,k,imx)*unp1-cons(i-1,j,k,imx)*unm1 
     $            + (q(i+1,j,k,qpres)-q(i-1,j,k,qpres)))               
     $            + BET*(cons(i+2,j,k,imx)*unp2-cons(i-2,j,k,imx)*unm2 
     $            + (q(i+2,j,k,qpres)-q(i-2,j,k,qpres)))               
     $            + GAM*(cons(i+3,j,k,imx)*unp3-cons(i-3,j,k,imx)*unm3 
     $            + (q(i+3,j,k,qpres)-q(i-3,j,k,qpres)))               
     $            + DEL*(cons(i+4,j,k,imx)*unp4-cons(i-4,j,k,imx)*unm4 
     $            + (q(i+4,j,k,qpres)-q(i-4,j,k,qpres))))

             flux(i,j,k,imy)= -dxinvX *
     $             (ALP*(cons(i+1,j,k,imy)*unp1-cons(i-1,j,k,imy)*unm1) 
     $            + BET*(cons(i+2,j,k,imy)*unp2-cons(i-2,j,k,imy)*unm2) 
     $            + GAM*(cons(i+3,j,k,imy)*unp3-cons(i-3,j,k,imy)*unm3) 
     $            + DEL*(cons(i+4,j,k,imy)*unp4-cons(i-4,j,k,imy)*unm4))

             flux(i,j,k,imz)= -dxinvX *
     $             (ALP*(cons(i+1,j,k,imz)*unp1-cons(i-1,j,k,imz)*unm1) 
     $            + BET*(cons(i+2,j,k,imz)*unp2-cons(i-2,j,k,imz)*unm2) 
     $            + GAM*(cons(i+3,j,k,imz)*unp3-cons(i-3,j,k,imz)*unm3) 
     $            + DEL*(cons(i+4,j,k,imz)*unp4-cons(i-4,j,k,imz)*unm4))

             flux(i,j,k,iene)= -dxinvX *
     $             (ALP*(cons(i+1,j,k,iene)*unp1-cons(i-1,j,k,iene)*unm1 
     $            + (q(i+1,j,k,qpres)*unp1-q(i-1,j,k,qpres)*unm1))       
     $            + BET*(cons(i+2,j,k,iene)*unp2-cons(i-2,j,k,iene)*unm2 
     $            + (q(i+2,j,k,qpres)*unp2-q(i-2,j,k,qpres)*unm2))       
     $            + GAM*(cons(i+3,j,k,iene)*unp3-cons(i-3,j,k,iene)*unm3 
     $            + (q(i+3,j,k,qpres)*unp3-q(i-3,j,k,qpres)*unm3))       
     $            + DEL*(cons(i+4,j,k,iene)*unp4-cons(i-4,j,k,iene)*unm4 
     $            + (q(i+4,j,k,qpres)*unp4-q(i-4,j,k,qpres)*unm4)))
          enddo
       enddo
      enddo
!$OMP END PARALLEL DO
      L1_end = parallel_wtime()

      L2_start = parallel_wtime()
!$OMP PARALLEL DO PRIVATE(i,j,k,unp1,unp2,unp3,unp4,
!$OMP&unm1,unm2,unm3,unm4) reduction(+ : L2iters)
      do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             L2iters = L2iters + 1

             unp1 = q(i,j+1,k,qv)
             unp2 = q(i,j+2,k,qv)
             unp3 = q(i,j+3,k,qv)
             unp4 = q(i,j+4,k,qv)

             unm1 = q(i,j-1,k,qv)
             unm2 = q(i,j-2,k,qv)
             unm3 = q(i,j-3,k,qv)
             unm4 = q(i,j-4,k,qv)

             flux(i,j,k,irho)=flux(i,j,k,irho) - 
     $             (ALP*(cons(i,j+1,k,imy)-cons(i,j-1,k,imy)) 
     $            + BET*(cons(i,j+2,k,imy)-cons(i,j-2,k,imy)) 
     $            + GAM*(cons(i,j+3,k,imy)-cons(i,j-3,k,imy)) 
     $            + DEL*(cons(i,j+4,k,imy)-cons(i,j-4,k,imy)))*dxinvY
          enddo

          do i=lo(1),hi(1)

             unp1 = q(i,j+1,k,qv)
             unp2 = q(i,j+2,k,qv)
             unp3 = q(i,j+3,k,qv)
             unp4 = q(i,j+4,k,qv)

             unm1 = q(i,j-1,k,qv)
             unm2 = q(i,j-2,k,qv)
             unm3 = q(i,j-3,k,qv)
             unm4 = q(i,j-4,k,qv)

             flux(i,j,k,imx)=flux(i,j,k,imx) - 
     $             (ALP*(cons(i,j+1,k,imx)*unp1-cons(i,j-1,k,imx)*unm1) 
     $            + BET*(cons(i,j+2,k,imx)*unp2-cons(i,j-2,k,imx)*unm2) 
     $            + GAM*(cons(i,j+3,k,imx)*unp3-cons(i,j-3,k,imx)*unm3) 
     $            + DEL*(cons(i,j+4,k,imx)*unp4-cons(i,j-4,k,imx)*unm4))*dxinvY
          enddo

          do i=lo(1),hi(1)

             unp1 = q(i,j+1,k,qv)
             unp2 = q(i,j+2,k,qv)
             unp3 = q(i,j+3,k,qv)
             unp4 = q(i,j+4,k,qv)

             unm1 = q(i,j-1,k,qv)
             unm2 = q(i,j-2,k,qv)
             unm3 = q(i,j-3,k,qv)
             unm4 = q(i,j-4,k,qv)

             flux(i,j,k,imy)=flux(i,j,k,imy) - 
     $             (ALP*(cons(i,j+1,k,imy)*unp1-cons(i,j-1,k,imy)*unm1 
     $            + (q(i,j+1,k,qpres)-q(i,j-1,k,qpres)))               
     $            + BET*(cons(i,j+2,k,imy)*unp2-cons(i,j-2,k,imy)*unm2 
     $            + (q(i,j+2,k,qpres)-q(i,j-2,k,qpres)))               
     $            + GAM*(cons(i,j+3,k,imy)*unp3-cons(i,j-3,k,imy)*unm3 
     $            + (q(i,j+3,k,qpres)-q(i,j-3,k,qpres)))               
     $            + DEL*(cons(i,j+4,k,imy)*unp4-cons(i,j-4,k,imy)*unm4 
     $            + (q(i,j+4,k,qpres)-q(i,j-4,k,qpres))))*dxinvY
          enddo

          do i=lo(1),hi(1)

             unp1 = q(i,j+1,k,qv)
             unp2 = q(i,j+2,k,qv)
             unp3 = q(i,j+3,k,qv)
             unp4 = q(i,j+4,k,qv)

             unm1 = q(i,j-1,k,qv)
             unm2 = q(i,j-2,k,qv)
             unm3 = q(i,j-3,k,qv)
             unm4 = q(i,j-4,k,qv)

             flux(i,j,k,imz)=flux(i,j,k,imz) - 
     $             (ALP*(cons(i,j+1,k,imz)*unp1-cons(i,j-1,k,imz)*unm1) 
     $            + BET*(cons(i,j+2,k,imz)*unp2-cons(i,j-2,k,imz)*unm2) 
     $            + GAM*(cons(i,j+3,k,imz)*unp3-cons(i,j-3,k,imz)*unm3) 
     $            + DEL*(cons(i,j+4,k,imz)*unp4-cons(i,j-4,k,imz)*unm4))*dxinvY
          enddo

          do i=lo(1),hi(1)

             unp1 = q(i,j+1,k,qv)
             unp2 = q(i,j+2,k,qv)
             unp3 = q(i,j+3,k,qv)
             unp4 = q(i,j+4,k,qv)

             unm1 = q(i,j-1,k,qv)
             unm2 = q(i,j-2,k,qv)
             unm3 = q(i,j-3,k,qv)
             unm4 = q(i,j-4,k,qv)

             flux(i,j,k,iene)=flux(i,j,k,iene) - 
     $             (ALP*(cons(i,j+1,k,iene)*unp1-cons(i,j-1,k,iene)*unm1 
     $            + (q(i,j+1,k,qpres)*unp1-q(i,j-1,k,qpres)*unm1))       
     $            + BET*(cons(i,j+2,k,iene)*unp2-cons(i,j-2,k,iene)*unm2 
     $            + (q(i,j+2,k,qpres)*unp2-q(i,j-2,k,qpres)*unm2))       
     $            + GAM*(cons(i,j+3,k,iene)*unp3-cons(i,j-3,k,iene)*unm3 
     $            + (q(i,j+3,k,qpres)*unp3-q(i,j-3,k,qpres)*unm3))       
     $            + DEL*(cons(i,j+4,k,iene)*unp4-cons(i,j-4,k,iene)*unm4 
     $            + (q(i,j+4,k,qpres)*unp4-q(i,j-4,k,qpres)*unm4)))*dxinvY
          enddo
       enddo
      enddo
!$OMP END PARALLEL DO
      L2_end = parallel_wtime()

      JBlockSize=11
      JBlocks = ((hi(2)-lo(2) + JBlockSize-1)/JBlockSize)

      L3_start = parallel_wtime()

!$OMP PARALLEL DO PRIVATE(jb,i,j,k,unp1,unp2,unp3,unp4,
!$OMP&unm1,unm2,unm3,unm4)  schedule(static,1)  reduction(+ : L3iters)
      do jb=0,JBlocks-1
      do k=lo(3),hi(3)
       !!!do j=lo(2),hi(2)
       do j=jb*JBlockSize, ((jb+1)*JBlockSize)-1
         if(j.le.hi(2)) then
          do i=lo(1),hi(1)

             L3iters = L3iters + 1

             unp1 = q(i,j,k+1,qw)
             unp2 = q(i,j,k+2,qw)
             unp3 = q(i,j,k+3,qw)
             unp4 = q(i,j,k+4,qw)

             unm1 = q(i,j,k-1,qw)
             unm2 = q(i,j,k-2,qw)
             unm3 = q(i,j,k-3,qw)
             unm4 = q(i,j,k-4,qw)

             flux(i,j,k,irho)=flux(i,j,k,irho) - 
     $             (ALP*(cons(i,j,k+1,imz)-cons(i,j,k-1,imz)) 
     $            + BET*(cons(i,j,k+2,imz)-cons(i,j,k-2,imz)) 
     $            + GAM*(cons(i,j,k+3,imz)-cons(i,j,k-3,imz)) 
     $            + DEL*(cons(i,j,k+4,imz)-cons(i,j,k-4,imz)))*dxinvZ
          enddo

          do i=lo(1),hi(1)

             unp1 = q(i,j,k+1,qw)
             unp2 = q(i,j,k+2,qw)
             unp3 = q(i,j,k+3,qw)
             unp4 = q(i,j,k+4,qw)

             unm1 = q(i,j,k-1,qw)
             unm2 = q(i,j,k-2,qw)
             unm3 = q(i,j,k-3,qw)
             unm4 = q(i,j,k-4,qw)

             flux(i,j,k,imx)=flux(i,j,k,imx) - 
     $             (ALP*(cons(i,j,k+1,imx)*unp1-cons(i,j,k-1,imx)*unm1) 
     $            + BET*(cons(i,j,k+2,imx)*unp2-cons(i,j,k-2,imx)*unm2) 
     $            + GAM*(cons(i,j,k+3,imx)*unp3-cons(i,j,k-3,imx)*unm3) 
     $            + DEL*(cons(i,j,k+4,imx)*unp4-cons(i,j,k-4,imx)*unm4))*dxinvZ
          enddo

          do i=lo(1),hi(1)

             unp1 = q(i,j,k+1,qw)
             unp2 = q(i,j,k+2,qw)
             unp3 = q(i,j,k+3,qw)
             unp4 = q(i,j,k+4,qw)

             unm1 = q(i,j,k-1,qw)
             unm2 = q(i,j,k-2,qw)
             unm3 = q(i,j,k-3,qw)
             unm4 = q(i,j,k-4,qw)

             flux(i,j,k,imy)=flux(i,j,k,imy) - 
     $             (ALP*(cons(i,j,k+1,imy)*unp1-cons(i,j,k-1,imy)*unm1) 
     $            + BET*(cons(i,j,k+2,imy)*unp2-cons(i,j,k-2,imy)*unm2) 
     $            + GAM*(cons(i,j,k+3,imy)*unp3-cons(i,j,k-3,imy)*unm3) 
     $            + DEL*(cons(i,j,k+4,imy)*unp4-cons(i,j,k-4,imy)*unm4))*dxinvZ
          enddo

          do i=lo(1),hi(1)

             unp1 = q(i,j,k+1,qw)
             unp2 = q(i,j,k+2,qw)
             unp3 = q(i,j,k+3,qw)
             unp4 = q(i,j,k+4,qw)

             unm1 = q(i,j,k-1,qw)
             unm2 = q(i,j,k-2,qw)
             unm3 = q(i,j,k-3,qw)
             unm4 = q(i,j,k-4,qw)

             flux(i,j,k,imz)=flux(i,j,k,imz) - 
     $             (ALP*(cons(i,j,k+1,imz)*unp1-cons(i,j,k-1,imz)*unm1 
     $            + (q(i,j,k+1,qpres)-q(i,j,k-1,qpres)))               
     $            + BET*(cons(i,j,k+2,imz)*unp2-cons(i,j,k-2,imz)*unm2 
     $            + (q(i,j,k+2,qpres)-q(i,j,k-2,qpres)))               
     $            + GAM*(cons(i,j,k+3,imz)*unp3-cons(i,j,k-3,imz)*unm3 
     $            + (q(i,j,k+3,qpres)-q(i,j,k-3,qpres)))               
     $            + DEL*(cons(i,j,k+4,imz)*unp4-cons(i,j,k-4,imz)*unm4 
     $            + (q(i,j,k+4,qpres)-q(i,j,k-4,qpres))))*dxinvZ
          enddo

          do i=lo(1),hi(1)

             unp1 = q(i,j,k+1,qw)
             unp2 = q(i,j,k+2,qw)
             unp3 = q(i,j,k+3,qw)
             unp4 = q(i,j,k+4,qw)

             unm1 = q(i,j,k-1,qw)
             unm2 = q(i,j,k-2,qw)
             unm3 = q(i,j,k-3,qw)
             unm4 = q(i,j,k-4,qw)

             flux(i,j,k,iene)=flux(i,j,k,iene) - 
     $             (ALP*(cons(i,j,k+1,iene)*unp1-cons(i,j,k-1,iene)*unm1 
     $            + (q(i,j,k+1,qpres)*unp1-q(i,j,k-1,qpres)*unm1))       
     $            + BET*(cons(i,j,k+2,iene)*unp2-cons(i,j,k-2,iene)*unm2 
     $            + (q(i,j,k+2,qpres)*unp2-q(i,j,k-2,qpres)*unm2))       
     $            + GAM*(cons(i,j,k+3,iene)*unp3-cons(i,j,k-3,iene)*unm3 
     $            + (q(i,j,k+3,qpres)*unp3-q(i,j,k-3,qpres)*unm3))       
     $            + DEL*(cons(i,j,k+4,iene)*unp4-cons(i,j,k-4,iene)*unm4 
     $            + (q(i,j,k+4,qpres)*unp4-q(i,j,k-4,qpres)*unm4)))*dxinvZ
          enddo
	 endif
       enddo
      enddo
      enddo
!$OMP END PARALLEL DO
      L3_end = parallel_wtime()

      if ( bIOP ) then
       print*, "L1iters = ", L1iters
       print*, "L2iters = ", L2iters
       print*, "L3iters = ", L3iters
      end if

      if ( bIOP ) then
       print *,"-----------------"
       write(6,42),"L1 time     (s) =",L1_end-L1_start
       write(6,42),"L2 time     (s) =",L2_end-L2_start
       write(6,42),"L3 time     (s) =",L3_end-L3_start
       print *,"-----------------"
       write(6,42),"Total time  (s) =",L3_end-L1_start
       print *,"-----------------"
      end if

42    format(a,f12.8)



       end


      subroutine FORT_HYPTERM_UNOPT(cons, DIMS(cons), DIMS(grid), q, flux, dx, nvars)

      implicit none

      integer DIMDEC(cons)
      integer DIMDEC(grid)
      integer nvars
      REAL_T  dx(SDIM)
      REAL_T  cons(DIMV(cons), nvars)
      REAL_T  q(DIMV(cons), nvars+1)
      REAL_T  flux(DIMV(grid), nvars)

      REAL_T unp1,unp2,unp3,unp4,unm1,unm2,unm3,unm4
      REAL_T dxinv(SDIM)

      REAL_T L1_start, L1_end
      REAL_T L2_start, L2_end
      REAL_T L3_start, L3_end

      integer irho, imx, imy, imz, iene, qu, qv, qw, qpres
      REAL_T ALP, BET, GAM, DEL

      REAL_T t, parallel_wtime

      integer lo(SDIM),hi(SDIM)
      integer JBlocks, JBlockSize, jb
      integer L1iters, L2iters, L3iters
      integer i, j, k



      logical bIOP

         write (6,*) grid_l1,grid_h1
         write (6,*) grid_l2,grid_h2
         write (6,*) grid_l3,grid_h3
         write (6,*) cons_l1,cons_h1
         write (6,*) cons_l2,cons_h2
         write (6,*) cons_l3,cons_h3


      lo(1) = grid_l1
      hi(1) = grid_h1
      lo(2) = grid_l2
      hi(2) = grid_h2
      lo(3) = grid_l3
      hi(3) = grid_h3


      irho = 1
      imx  = 2
      imy  = 3
      imz  = 4
      iene = 5

      qu    = 2
      qv    = 3
      qw    = 4
      qpres = 5

      ALP =  0.8d0
      BET = -0.2d0
      GAM =  4.d0/105.d0
      DEL = -1.d0/280.d0


      L1iters = 0
      L2iters = 0
      L3iters = 0
      bIOP = .true.


      do i=1,3
         dxinv(i) = 1.0d0 / dx(i)
      end do


      if ( bIOP ) then
       print*, "lo(1) hi(1) = ", lo(1), hi(1)
       print*, "lo(2) hi(2) = ", lo(2), hi(2)
       print*, "lo(3) hi(3) = ", lo(3), hi(3)
      end if


      L1_start = parallel_wtime()

!$OMP PARALLEL DO PRIVATE(i,j,k,unp1,unp2,unp3,unp4,unm1,unm2,unm3,unm4) reduction(+ : L1iters)
      do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             L1iters = L1iters + 1

             unm4 = q(i-4,j,k,qu)
             unm3 = q(i-3,j,k,qu)
             unm2 = q(i-2,j,k,qu)
             unm1 = q(i-1,j,k,qu)

             unp1 = q(i+1,j,k,qu)
             unp2 = q(i+2,j,k,qu)
             unp3 = q(i+3,j,k,qu)
             unp4 = q(i+4,j,k,qu)

             flux(i,j,k,irho)= - 
     $             (ALP*(cons(i+1,j,k,imx)-cons(i-1,j,k,imx)) 
     $            + BET*(cons(i+2,j,k,imx)-cons(i-2,j,k,imx)) 
     $            + GAM*(cons(i+3,j,k,imx)-cons(i-3,j,k,imx)) 
     $            + DEL*(cons(i+4,j,k,imx)-cons(i-4,j,k,imx)))*dxinv(1)

             flux(i,j,k,imx)= - 
     $             (ALP*(cons(i+1,j,k,imx)*unp1-cons(i-1,j,k,imx)*unm1 
     $            + (q(i+1,j,k,qpres)-q(i-1,j,k,qpres)))               
     $            + BET*(cons(i+2,j,k,imx)*unp2-cons(i-2,j,k,imx)*unm2 
     $            + (q(i+2,j,k,qpres)-q(i-2,j,k,qpres)))               
     $            + GAM*(cons(i+3,j,k,imx)*unp3-cons(i-3,j,k,imx)*unm3 
     $            + (q(i+3,j,k,qpres)-q(i-3,j,k,qpres)))               
     $            + DEL*(cons(i+4,j,k,imx)*unp4-cons(i-4,j,k,imx)*unm4 
     $            + (q(i+4,j,k,qpres)-q(i-4,j,k,qpres))))*dxinv(1)

             flux(i,j,k,imy)= - 
     $             (ALP*(cons(i+1,j,k,imy)*unp1-cons(i-1,j,k,imy)*unm1) 
     $            + BET*(cons(i+2,j,k,imy)*unp2-cons(i-2,j,k,imy)*unm2) 
     $            + GAM*(cons(i+3,j,k,imy)*unp3-cons(i-3,j,k,imy)*unm3) 
     $            + DEL*(cons(i+4,j,k,imy)*unp4-cons(i-4,j,k,imy)*unm4))*dxinv(1)

             flux(i,j,k,imz)= - 
     $             (ALP*(cons(i+1,j,k,imz)*unp1-cons(i-1,j,k,imz)*unm1) 
     $            + BET*(cons(i+2,j,k,imz)*unp2-cons(i-2,j,k,imz)*unm2) 
     $            + GAM*(cons(i+3,j,k,imz)*unp3-cons(i-3,j,k,imz)*unm3) 
     $            + DEL*(cons(i+4,j,k,imz)*unp4-cons(i-4,j,k,imz)*unm4))*dxinv(1)

             flux(i,j,k,iene)= - 
     $             (ALP*(cons(i+1,j,k,iene)*unp1-cons(i-1,j,k,iene)*unm1 
     $            + (q(i+1,j,k,qpres)*unp1-q(i-1,j,k,qpres)*unm1))       
     $            + BET*(cons(i+2,j,k,iene)*unp2-cons(i-2,j,k,iene)*unm2 
     $            + (q(i+2,j,k,qpres)*unp2-q(i-2,j,k,qpres)*unm2))       
     $            + GAM*(cons(i+3,j,k,iene)*unp3-cons(i-3,j,k,iene)*unm3 
     $            + (q(i+3,j,k,qpres)*unp3-q(i-3,j,k,qpres)*unm3))       
     $            + DEL*(cons(i+4,j,k,iene)*unp4-cons(i-4,j,k,iene)*unm4 
     $            + (q(i+4,j,k,qpres)*unp4-q(i-4,j,k,qpres)*unm4)))*dxinv(1) 
          enddo
       enddo
      enddo
!$OMP END PARALLEL DO
      L1_end = parallel_wtime()

      L2_start = parallel_wtime()
!$OMP PARALLEL DO PRIVATE(i,j,k,unp1,unp2,unp3,unp4,unm1,unm2,unm3,unm4) reduction(+ : L2iters)
      do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             L2iters = L2iters + 1

             unp1 = q(i,j+1,k,qv)
             unp2 = q(i,j+2,k,qv)
             unp3 = q(i,j+3,k,qv)
             unp4 = q(i,j+4,k,qv)

             unm1 = q(i,j-1,k,qv)
             unm2 = q(i,j-2,k,qv)
             unm3 = q(i,j-3,k,qv)
             unm4 = q(i,j-4,k,qv)

             flux(i,j,k,irho)=flux(i,j,k,irho) - 
     $             (ALP*(cons(i,j+1,k,imy)-cons(i,j-1,k,imy)) 
     $            + BET*(cons(i,j+2,k,imy)-cons(i,j-2,k,imy)) 
     $            + GAM*(cons(i,j+3,k,imy)-cons(i,j-3,k,imy)) 
     $            + DEL*(cons(i,j+4,k,imy)-cons(i,j-4,k,imy)))*dxinv(2)

             flux(i,j,k,imx)=flux(i,j,k,imx) - 
     $             (ALP*(cons(i,j+1,k,imx)*unp1-cons(i,j-1,k,imx)*unm1) 
     $            + BET*(cons(i,j+2,k,imx)*unp2-cons(i,j-2,k,imx)*unm2) 
     $            + GAM*(cons(i,j+3,k,imx)*unp3-cons(i,j-3,k,imx)*unm3) 
     $            + DEL*(cons(i,j+4,k,imx)*unp4-cons(i,j-4,k,imx)*unm4))*dxinv(2)

             flux(i,j,k,imy)=flux(i,j,k,imy) - 
     $             (ALP*(cons(i,j+1,k,imy)*unp1-cons(i,j-1,k,imy)*unm1 
     $            + (q(i,j+1,k,qpres)-q(i,j-1,k,qpres)))               
     $            + BET*(cons(i,j+2,k,imy)*unp2-cons(i,j-2,k,imy)*unm2 
     $            + (q(i,j+2,k,qpres)-q(i,j-2,k,qpres)))               
     $            + GAM*(cons(i,j+3,k,imy)*unp3-cons(i,j-3,k,imy)*unm3 
     $            + (q(i,j+3,k,qpres)-q(i,j-3,k,qpres)))               
     $            + DEL*(cons(i,j+4,k,imy)*unp4-cons(i,j-4,k,imy)*unm4 
     $            + (q(i,j+4,k,qpres)-q(i,j-4,k,qpres))))*dxinv(2)

             flux(i,j,k,imz)=flux(i,j,k,imz) - 
     $             (ALP*(cons(i,j+1,k,imz)*unp1-cons(i,j-1,k,imz)*unm1) 
     $            + BET*(cons(i,j+2,k,imz)*unp2-cons(i,j-2,k,imz)*unm2) 
     $            + GAM*(cons(i,j+3,k,imz)*unp3-cons(i,j-3,k,imz)*unm3) 
     $            + DEL*(cons(i,j+4,k,imz)*unp4-cons(i,j-4,k,imz)*unm4))*dxinv(2)

             flux(i,j,k,iene)=flux(i,j,k,iene) - 
     $             (ALP*(cons(i,j+1,k,iene)*unp1-cons(i,j-1,k,iene)*unm1 
     $            + (q(i,j+1,k,qpres)*unp1-q(i,j-1,k,qpres)*unm1))       
     $            + BET*(cons(i,j+2,k,iene)*unp2-cons(i,j-2,k,iene)*unm2 
     $            + (q(i,j+2,k,qpres)*unp2-q(i,j-2,k,qpres)*unm2))       
     $            + GAM*(cons(i,j+3,k,iene)*unp3-cons(i,j-3,k,iene)*unm3 
     $            + (q(i,j+3,k,qpres)*unp3-q(i,j-3,k,qpres)*unm3))       
     $            + DEL*(cons(i,j+4,k,iene)*unp4-cons(i,j-4,k,iene)*unm4 
     $            + (q(i,j+4,k,qpres)*unp4-q(i,j-4,k,qpres)*unm4)))*dxinv(2)
          enddo
       enddo
      enddo
!$OMP END PARALLEL DO
      L2_end = parallel_wtime()


      L3_start = parallel_wtime()

!$OMP PARALLEL DO PRIVATE(jb,i,j,k,unp1,unp2,unp3,unp4,
!$OMP&unm1,unm2,unm3,unm4)  schedule(static,1) reduction(+ : L3iters)
      do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             L3iters = L3iters + 1

             unp1 = q(i,j,k+1,qw)
             unp2 = q(i,j,k+2,qw)
             unp3 = q(i,j,k+3,qw)
             unp4 = q(i,j,k+4,qw)

             unm1 = q(i,j,k-1,qw)
             unm2 = q(i,j,k-2,qw)
             unm3 = q(i,j,k-3,qw)
             unm4 = q(i,j,k-4,qw)

             flux(i,j,k,irho)=flux(i,j,k,irho) - 
     $             (ALP*(cons(i,j,k+1,imz)-cons(i,j,k-1,imz)) 
     $            + BET*(cons(i,j,k+2,imz)-cons(i,j,k-2,imz)) 
     $            + GAM*(cons(i,j,k+3,imz)-cons(i,j,k-3,imz)) 
     $            + DEL*(cons(i,j,k+4,imz)-cons(i,j,k-4,imz)))*dxinv(3)

             flux(i,j,k,imx)=flux(i,j,k,imx) - 
     $             (ALP*(cons(i,j,k+1,imx)*unp1-cons(i,j,k-1,imx)*unm1) 
     $            + BET*(cons(i,j,k+2,imx)*unp2-cons(i,j,k-2,imx)*unm2) 
     $            + GAM*(cons(i,j,k+3,imx)*unp3-cons(i,j,k-3,imx)*unm3) 
     $            + DEL*(cons(i,j,k+4,imx)*unp4-cons(i,j,k-4,imx)*unm4))*dxinv(3)

             flux(i,j,k,imy)=flux(i,j,k,imy) - 
     $             (ALP*(cons(i,j,k+1,imy)*unp1-cons(i,j,k-1,imy)*unm1) 
     $            + BET*(cons(i,j,k+2,imy)*unp2-cons(i,j,k-2,imy)*unm2) 
     $            + GAM*(cons(i,j,k+3,imy)*unp3-cons(i,j,k-3,imy)*unm3) 
     $            + DEL*(cons(i,j,k+4,imy)*unp4-cons(i,j,k-4,imy)*unm4))*dxinv(3)

             flux(i,j,k,imz)=flux(i,j,k,imz) - 
     $             (ALP*(cons(i,j,k+1,imz)*unp1-cons(i,j,k-1,imz)*unm1 
     $            + (q(i,j,k+1,qpres)-q(i,j,k-1,qpres)))               
     $            + BET*(cons(i,j,k+2,imz)*unp2-cons(i,j,k-2,imz)*unm2 
     $            + (q(i,j,k+2,qpres)-q(i,j,k-2,qpres)))               
     $            + GAM*(cons(i,j,k+3,imz)*unp3-cons(i,j,k-3,imz)*unm3 
     $            + (q(i,j,k+3,qpres)-q(i,j,k-3,qpres)))               
     $            + DEL*(cons(i,j,k+4,imz)*unp4-cons(i,j,k-4,imz)*unm4 
     $            + (q(i,j,k+4,qpres)-q(i,j,k-4,qpres))))*dxinv(3)

             flux(i,j,k,iene)=flux(i,j,k,iene) - 
     $             (ALP*(cons(i,j,k+1,iene)*unp1-cons(i,j,k-1,iene)*unm1 
     $            + (q(i,j,k+1,qpres)*unp1-q(i,j,k-1,qpres)*unm1))       
     $            + BET*(cons(i,j,k+2,iene)*unp2-cons(i,j,k-2,iene)*unm2 
     $            + (q(i,j,k+2,qpres)*unp2-q(i,j,k-2,qpres)*unm2))       
     $            + GAM*(cons(i,j,k+3,iene)*unp3-cons(i,j,k-3,iene)*unm3 
     $            + (q(i,j,k+3,qpres)*unp3-q(i,j,k-3,qpres)*unm3))       
     $            + DEL*(cons(i,j,k+4,iene)*unp4-cons(i,j,k-4,iene)*unm4 
     $            + (q(i,j,k+4,qpres)*unp4-q(i,j,k-4,qpres)*unm4)))*dxinv(3)
          enddo
       enddo
      enddo
!$OMP END PARALLEL DO
      L3_end = parallel_wtime()

      if ( bIOP ) then
       print*, "L1iters = ", L1iters
       print*, "L2iters = ", L2iters
       print*, "L3iters = ", L3iters
      end if

      if ( bIOP ) then
       write(6,42),"L1 time     (s) =",L1_end-L1_start
       write(6,42),"L2 time     (s) =",L2_end-L2_start
       write(6,42),"L3 time     (s) =",L3_end-L3_start
       write(6,42),"Total time  (s) =",L3_end-L1_start
      end if

42    format(a,f12.8)



       end


