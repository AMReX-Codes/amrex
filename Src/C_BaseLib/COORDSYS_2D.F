c
c $Id: COORDSYS_2D.F,v 1.6 1999-02-03 00:16:31 lijewski Exp $
c

#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "COORDSYS_F.H"
#include "ArrayLim.H"

#define SDIM 2

c :: ----------------------------------------------------------
c :: SETVOL
c ::             Compute the volume of each cell
c ::
c :: INPUTS / OUTPUTS:
c ::  vol         <=  volume array
c ::  vlo,vhi      => index limits of vol array
c ::  offset       => shift to origin of computational domain
c ::  dx           => cell size
c ::  coord        => coordinate flag (0 = cartesian, 1 = RZ)
c :: ----------------------------------------------------------
c ::
       subroutine FORT_SETVOL(vol,DIMS(vol),offset,dx,coord)
       integer    DIMDEC(vol)
       integer    coord
       REAL_T     dx(SDIM), offset(SDIM)
       REAL_T     vol(DIMV(vol))
       
       integer    i, j
       REAL_T     ri, ro, v
       REAL_T     RZFACTOR
       parameter (RZFACTOR = 2*Pi)
       
       if (coord .eq. 0) then
c
c         ::::: cartesian
c
          v = dx(1)*dx(2)
          do j = ARG_L2(vol), ARG_H2(vol)
             do i = ARG_L1(vol), ARG_H1(vol)
                vol(i,j) = v
             end do
          end do
       else
c
c         ::::: R-Z
c
          do i = ARG_L1(vol), ARG_H1(vol)
             ri = offset(1) + dx(1)*i
             ro = ri + dx(1)
             v = (half*RZFACTOR)*dx(2)*dx(1)*(ro + ri)
             do j = ARG_L2(vol), ARG_H2(vol)
                vol(i,j) = v
             end do
          end do
       end if
       
       end

c========================================================
c========================================================
        subroutine FORT_SETVOLPT(vol,
     $          ro, ri, po, pi,  dx, coord)

        integer coord
        REAL_T dx(SDIM)
        REAL_T     vol
        REAL_T     ro, po, pi
        REAL_T     ri

        integer i,j
        REAL_T     RZFACTOR
        parameter (RZFACTOR = 2*Pi)

        if(coord .eq. 0) then
           vol = (ro-ri)*dx(2)
        elseif(coord .eq. 1) then
           vol = half*RZFACTOR*dx(2)*(ro**2 - ri**2)
           vol = abs(vol)
        elseif(coord .eq. 2) then
           vol = RZFACTOR*(ro**3-ri**3)*(cos(pi)-cos(po))/three
        else
            write(6,*) ' bogus value of coord ... bndrylib::SETVOLPT'
            stop
        endif

        return
        end



c :: ----------------------------------------------------------
c :: SETDLOGA
c ::             Compute  d(log(A))/dr in each cell
c ::
c :: INPUTS / OUTPUTS:
c ::  dloga        <=  dloga array
c ::  dlo,dhi      => index limits of dloga array
c ::  offset       => shift to origin of computational domain
c ::  dx           => cell size
c ::  coord        => coordinate flag (0 = cartesian, 1 = RZ)
c :: ----------------------------------------------------------
c ::
       subroutine FORT_SETDLOGA(dloga,DIMS(dloga),offset,dx,dir,coord)

       integer    DIMDEC(dloga)
       integer    coord
       REAL_T     dx(SDIM), offset(SDIM)
       REAL_T     dloga(DIMV(dloga))
       integer dir
       
       integer    i, j
       REAL_T     ri, ro, dlga
       
       if (coord .eq. 0) then
c
c         ::::: cartesian
c
          do j = ARG_L2(dloga), ARG_H2(dloga)
             do i = ARG_L1(dloga), ARG_H1(dloga)
                dloga(i,j) = zero
             end do
          end do

       else if( coord .eq. 1 ) then
c
c         ::::: R-Z
c
          if( dir .eq. 0 ) then
             do i = ARG_L1(dloga), ARG_H1(dloga)
                ri = offset(1) + dx(1)*i
                ro = ri + dx(1)
                dlga = two/(ro+ri)                                        
                do j = ARG_L2(dloga), ARG_H2(dloga)
                   dloga(i,j) = dlga
                end do
             end do
          else if( dir .eq. 1 ) then
             do i = ARG_L1(dloga), ARG_H1(dloga)
                ri = offset(1) + dx(1)*i
                ro = ri + dx(1)
                dlga = two/(ro+ri)                                        
                do j = ARG_L2(dloga), ARG_H2(dloga)
                   dloga(i,j) = zero
                end do
             end do
          else
             write(6,*)'setdloga: illega direction ',dir
             stop
          end if
       else
          write(6,*)'setdloga: illegal coordinate system ',coord
          stop
       end if
       
       end

c :: ----------------------------------------------------------
c :: SETAREA
c ::             Compute the area of given cell face
c ::
c :: INPUTS / OUTPUTS:
c ::  area        <=  area array
c ::  alo,ahi      => index limits of area array
c ::  offset       => shift to origin of computational domain
c ::  dx           => cell size
c ::  coord        => coordinate flag (0 =cartesian, 1 = RZ)
c :: ----------------------------------------------------------
c ::
       subroutine FORT_SETAREA(area,DIMS(area),offset,dx,dir,coord)

       integer    DIMDEC(area)
       integer    coord, dir
       REAL_T     dx(SDIM), offset(SDIM)
       REAL_T     area(DIMV(area))
       
       integer    i, j
       REAL_T     ri, ro, a
       REAL_T     RZFACTOR
       parameter (RZFACTOR = 2*Pi)
       
       if (coord .eq. 0) then
c
c         ::::: cartesian
c
          if (dir .eq. 0) then
             do j = ARG_L2(area), ARG_H2(area)
                do i = ARG_L1(area), ARG_H1(area)
                   area(i,j) = dx(2)
                end do
             end do
          else
             do j = ARG_L2(area), ARG_H2(area)
                do i = ARG_L1(area), ARG_H1(area)
                   area(i,j) = dx(1)
                end do
             end do
          end if
       else
c
c         ::::: R-Z
c
          if (dir .eq. 0) then
             do i = ARG_L1(area), ARG_H1(area)
                ri = offset(1) + dx(1)*i
                a = RZFACTOR*ri*dx(2)
                do j = ARG_L2(area), ARG_H2(area)
                   area(i,j) = a
                end do
             end do
          else
             do i = ARG_L1(area), ARG_H1(area)
                ri = offset(1) + dx(1)*i
                ro = ri + dx(1)
                a = dx(1)*(half*RZFACTOR)*(ri + ro)
                do j = ARG_L2(area), ARG_H2(area)
                   area(i,j) = a
                end do
             end do
          end if
       end if
       
       end
