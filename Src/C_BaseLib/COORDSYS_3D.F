c
c $Id: COORDSYS_3D.F,v 1.7 2001-08-15 19:23:55 lijewski Exp $
c

#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "COORDSYS_F.H"
#include "ArrayLim.H"

#define SDIM 3

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
       
       integer    i, j, k
       REAL_T     ri, ro, v
       
       if (coord .eq. 0) then
c
c         ::::: cartesian
c
          v = dx(1)*dx(2)*dx(3)
          do k = ARG_L3(vol), ARG_H3(vol)
             do j = ARG_L2(vol), ARG_H2(vol)
                do i = ARG_L1(vol), ARG_H1(vol)
                   vol(i,j,k) = v
                end do
             end do
          end do
       else
          write(6,*) "FORT_SETVOLUME not define for coord = ",coord
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
       
       integer    i, j, k
       REAL_T     fa
       
       if (coord .eq. 0) then
          if (dir .eq. 0) then
             fa = dx(2)*dx(3)
          else if (dir .eq. 1) then
             fa = dx(1)*dx(3)
          else if (dir .eq. 2) then
             fa = dx(1)*dx(2)
          else
             write(6,*) "FORT_SETAREA: invalid dir = ",dir
             stop
          end if
          do k = ARG_L3(area), ARG_H3(area)
             do j = ARG_L2(area), ARG_H2(area)
                do i = ARG_L1(area), ARG_H1(area)
                   area(i,j,k) = fa
                end do
             end do
          end do
       else
          write(6,*) "FORT_SETAREA not define for coord = ",coord
          stop
       end if
       
       end

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
       
       integer    i, j, k
       REAL_T     ri, ro, dlga, po, pi
       
       if (coord .eq. 0) then
c
c         ::::: cartesian
c
          do k = ARG_L3(dloga), ARG_H3(dloga)
          do j = ARG_L2(dloga), ARG_H2(dloga)
             do i = ARG_L1(dloga), ARG_H1(dloga)
                dloga(i,j,k) = zero
             end do
          end do
	  enddo

       else 
	 write(6,*)' non-cartesian not allowed in 3D yet'
	stop

	endif

	return
	end
