#undef  BL_LANG_CC
#define BL_LANG_FORT

#include "REAL.H"
#include "CONSTANTS.H"
#include "COORDSYS_F.H"
#include "DIMS.H"

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
       subroutine FORT_SETVOL(vol,DIMS(vol),offset,dx,coord)
       integer    DIMDEC(vol)
       integer    coord
       REAL_T     dx(SDIM), offset(SDIM)
       REAL_T     vol(DIMV(vol))
       
       integer    i, j
       REAL_T     ri, ro, v
       REAL_T     RZFACTOR
       parameter (RZFACTOR = 2.d0*3.14159265358979323846d0)
       
       if (coord .eq. 0) then
c         ::::: cartesian
          v = dx(1)
             do i = ARG_L1(vol), ARG_H1(vol)
                vol(i) = v
             enddo
       else
	 if(coord .eq. 1) then
c  ::: 1D cylindrical
           do i = ARG_L1(vol), ARG_H1(vol)
              ri = offset(1) + dx(1)*i
              ro = ri + dx(1)
              v = half*(RZFACTOR)*(ro**2 - ri**2)
              vol(i) = abs(v)
           enddo
	 else
	    if(coord .eq. 2) then
c         ::::: 1D spherical
               do i = ARG_L1(vol), ARG_H1(vol)
                  ri = offset(1) + dx(1)*i
                  ro = ri + dx(1)
                  v = (two3rd*RZFACTOR)*(ro**3 - ri**3)
                  vol(i) = abs(v)
               enddo
	     else
	        write(6,*)' bogus value of coord... bndrylib::SETVOL '
	        stop
	     endif
	 endif
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
	
	integer i,j
        REAL_T     RZFACTOR
        parameter (RZFACTOR = 2.d0*3.14159265358979323846d0)

c
c  note that dx is usually unity.  dx not unity is used by the nfluid
c  slic reconstruction
c

	if(coord .eq. 0) then
	   do i = roloi1, rohii1
	      vol(i) = (ro(i)-ri(i))
	   enddo

	else
	   if(coord .eq. 1) then
	      do i = roloi1, rohii1
	         vol(i) = half*RZFACTOR*(ro(i)**2 - ri(i)**2)
	 	 vol(i) = abs(vol(i))
	      enddo
	   else
	      if (coord .eq. 2) then
	         do i = roloi1, rohii1
	            vol(i) = two3rd*RZFACTOR*(ro(i)**3-ri(i)**3)
	 	 vol(i) = abs(vol(i))
	         enddo
	      else
	         write(6,*) ' bogus value of coord ... bndrylib::SETVOLPT'
	         stop
	      endif
	   endif
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

c         ::::: cartesian
             do i = ARG_L1(dloga), ARG_H1(dloga)
                dloga(i) = zero
             enddo

       else 
c  :::  1D cylindrical
	  if(coord .eq. 1) then
            do i = ARG_L1(dloga), ARG_H1(dloga)
               ri = offset(1) + dx(1)*i
	       ri = abs(ri)
               ro = ri + dx(1)
	       ro = abs(ro)
               dlga = two/(ro+ri)                                        
               dloga(i) = dlga
            enddo
	  else
	     if( coord .eq. 2 ) then
c         ::::: 1d spherical
                do i = ARG_L1(dloga), ARG_H1(dloga)
                   ri = offset(1) + dx(1)*i
	           ri = abs(ri)
                   ro = ri + dx(1)
	           ro = abs(ro)
                   dlga = four/(ro+ri)                                        
                   dloga(i) = dlga
                enddo
             else
                write(6,*)'setdloga: illegal coordinate system ',coord
                stop
             endif
	  endif
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
       subroutine FORT_SETAREA(area,DIMS(area),offset,dx,dir,coord)

       integer    DIMDEC(area)
       integer    coord, dir
       REAL_T     dx(SDIM), offset(SDIM)
       REAL_T     area(DIMV(area))
       
       integer    i, j
       REAL_T     ri, ro, a
       REAL_T     RZFACTOR
       parameter (RZFACTOR = 2.d0*3.14159265358979323846d0)
       
       if (coord .eq. 0) then
c         ::::: cartesian
            do i = ARG_L1(area), ARG_H1(area)
               area(i) = one
            enddo
       else
c         ::::: 1d spherical

	  if(coord .eq. 1) then
            do i = ARG_L1(area), ARG_H1(area)
               ri = offset(1) + dx(1)*i
	       a = RZFACTOR*ri
	       area(i) = abs(a)
            enddo
	  else
	      if( coord .eq. 2) then
                do i = ARG_L1(area), ARG_H1(area)
                   ri = offset(1) + dx(1)*i
	           a = two*RZFACTOR*ri*ri
	           area(i) = abs(a)
                enddo
	      else
	        write(6,*)' bogus value for coord... SETAREA '
	        stop
	      endif
	  endif
       endif
       
       return
       end
