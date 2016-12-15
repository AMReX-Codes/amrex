 
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif
 
#include "REAL.H"
#include "CONSTANTS.H"
#include "COORDSYS_F.H"
#include "ArrayLim.H"
 
#define SDIM 1

c :: ----------------------------------------------------------
c :: SETVOL
c ::             Compute the volume of each cell
c ::
c :: INPUTS / OUTPUTS:
c ::  vol         <=  volume array
c ::  vlo,vhi      => index limits of vol array
c ::  offset       => shift to origin of computational domain
c ::  dx           => cell size
c ::  coord        => coordinate flag (0 = cartesian, 1 = RZ, 2 = spherical)
c :: ----------------------------------------------------------
c ::
       subroutine FORT_SETVOL(DIMS(reg),vol,DIMS(vol),offset,dx,coord)
       implicit none
       integer    DIMDEC(reg)
       integer    DIMDEC(vol)
       integer    coord
       REAL_T     dx(SDIM), offset(SDIM)
       REAL_T     vol(DIMV(vol))
       
       integer    i
       REAL_T     ri, ro, v
       REAL_T     RZFACTOR
       parameter (RZFACTOR = 2.d0*3.14159265358979323846d0)
       
       if (coord .eq. 0) then
          v = dx(1)
             do i = ARG_L1(reg), ARG_H1(reg)
                vol(i) = v
             enddo

       else if (coord .eq. 1) then

          do i = ARG_L1(reg), ARG_H1(reg)
             ri = offset(1) + dx(1)*i
             ro = ri + dx(1)
             v = half*(RZFACTOR)*(ro**2 - ri**2)
             vol(i) = abs(v)
          enddo

       else if (coord .eq. 2) then

          do i = ARG_L1(reg), ARG_H1(reg)
             ri = offset(1) + dx(1)*i
             ro = ri + dx(1)
             v = (two3rd*RZFACTOR)*(ro**3 - ri**3)
             vol(i) = abs(v)
          enddo
       else 
          call bl_abort('bogus value of coord... bndrylib::SETVOL')
       endif
       
       return
       end
	
	subroutine FORT_SETVOLPT(vol, volloi1, volhii1, 
     $		ro, roloi1, rohii1, ri, riloi1, rihii1, dx, coord)

	integer volloi1, volhii1
	integer roloi1, rohii1, riloi1, rihii1
	integer coord
	REAL_T dx(SDIM)
	REAL_T vol(volloi1:volhii1)
	REAL_T ro(roloi1:rohii1)
	REAL_T ri(riloi1:rihii1)
	
	integer i
        REAL_T     RZFACTOR
        parameter (RZFACTOR = 2.d0*3.14159265358979323846d0)

c
c  note that dx is usually unity.  dx not unity is used by the nfluid
c  slic reconstruction
c

	if (coord .eq. 0) then
	   do i = roloi1, rohii1
	      vol(i) = (ro(i)-ri(i))
	   enddo

	else if (coord .eq. 1) then
           do i = roloi1, rohii1
              vol(i) = half*RZFACTOR*(ro(i)**2 - ri(i)**2)
              vol(i) = abs(vol(i))
           enddo

        else if (coord .eq. 2) then
           do i = roloi1, rohii1
	      vol(i) = two3rd*RZFACTOR*(ro(i)**3-ri(i)**3)
	      vol(i) = abs(vol(i))
           enddo

        else
           call bl_abort('bogus value of coord ... bndrylib::SETVOLPT')
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
       
       integer    i
       REAL_T     rc
       
       if (coord .eq. 0) then

          do i = ARG_L1(dloga), ARG_H1(dloga)
             dloga(i) = zero
          enddo

       else if (coord .eq. 1) then

          do i = ARG_L1(dloga), ARG_H1(dloga)
             rc = offset(1) + dx(1)*(dble(i) + 0.5d0)
             dloga(i) = 1.d0 / rc
          enddo

       else if (coord .eq. 2) then

          do i = ARG_L1(dloga), ARG_H1(dloga)
             rc = offset(1) + dx(1)*(dfloat(i) + 0.5d0)
             dloga(i) = 2.d0/rc
          enddo

       else
          call bl_abort('SETDLOGA: illegal coordinate system')
       endif
       return
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
c ::  coord        => coordinate flag (0 =cartesian, 1 = RZ, 2 = spherical)
c :: ----------------------------------------------------------
c ::
       subroutine FORT_SETAREA(DIMS(reg),area,DIMS(area),offset,dx,dir,coord)

       integer    DIMDEC(reg)
       integer    DIMDEC(area)
       integer    coord, dir
       REAL_T     dx(SDIM), offset(SDIM)
       REAL_T     area(DIMV(area))
       
       integer    i
       REAL_T     ri, a
       REAL_T     RZFACTOR
       parameter (RZFACTOR = 2.d0*3.14159265358979323846d0)
       
       if (coord .eq. 0) then
         do i = ARG_L1(reg), ARG_H1(reg)
            area(i) = one
         enddo

       else if (coord .eq. 1) then
          do i = ARG_L1(reg), ARG_H1(reg)
             ri = offset(1) + dx(1)*dble(i)
             a = RZFACTOR*ri
             area(i) = abs(a)
          enddo

       else if( coord .eq. 2) then
          do i = ARG_L1(reg), ARG_H1(reg)
             ri = offset(1) + dx(1)*dble(i)
             a = two*RZFACTOR*ri*ri
             area(i) = abs(a)
          enddo

       else
          call bl_abort('bogus value for coord... SETAREA')
       endif
       
       return
       end
