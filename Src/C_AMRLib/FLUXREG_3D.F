#undef  BL_LANG_CC
#define BL_LANG_FORT

#include "REAL.H"
#include "CONSTANTS.H"
#include "FLUXREG_F.H"
#include "DIMS.H"

#define SDIM 3

c ::: -----------------------------------------------------------
c ::: Init coarse grid flux into flux register
c :::
c ::: INPUTS/OUTPUTS:
c ::: flx        => flux array
c ::: DIMS(flx)  => index limits for flx
c ::: reg       <=  flux register
c ::: DIMS(reg)  => index limits for reg
c ::: lo,hi      => region of update
c ::: mult       => scalar multiplicative factor
c ::: -----------------------------------------------------------

      subroutine FORT_FRCRSEINIT (reg,DIMS(reg),flx,DIMS(flx),lo,hi,
     &                            numcomp,dir,mult)
      integer    DIMDEC(reg)
      integer    DIMDEC(flx)
      integer    lo(SDIM), hi(SDIM)
      integer    numcomp, dir
      REAL_T     mult
      REAL_T     reg(DIMV(reg),numcomp)
      REAL_T     flx(DIMV(flx),numcomp)

      integer    n, i, j, k

      do n = 1, numcomp
      do k = lo(3), hi(3)
      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            reg(i,j,k,n) =  mult*flx(i,j,k,n)
         enddo
      enddo
      enddo
      enddo

      return
      end

c ::: -----------------------------------------------------------
c ::: Init coarse grid flux (times area) into flux register
c :::
c ::: INPUTS/OUTPUTS:
c ::: flx        => flux array
c ::: DIMS(flx)  => index limits for flx
c ::: area       => aread of cell edge
c ::: DIMS(area) => index limits for area
c ::: reg       <=  flux register
c ::: DIMS(reg)  => index limits for reg
c ::: lo,hi      => region of update
c ::: mult       => scalar multiplicative factor
c ::: -----------------------------------------------------------

      subroutine FORT_FRCAINIT (reg,DIMS(reg),flx,DIMS(flx),
     &                          area,DIMS(area),lo,hi,
     &                          numcomp,dir,mult)
      integer    DIMDEC(reg)
      integer    DIMDEC(flx)
      integer    DIMDEC(area)
      integer    lo(SDIM), hi(SDIM)
      integer    numcomp, dir
      REAL_T     mult
      REAL_T     reg(DIMV(reg),numcomp)
      REAL_T     flx(DIMV(flx),numcomp)
      REAL_T    area(DIMV(area))

      integer    n, i, j, k

      do n = 1, numcomp
      do k = lo(3), hi(3)
      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            reg(i,j,k,n) = mult*area(i,j,k)*flx(i,j,k,n)
         enddo
      enddo
      enddo
      enddo

      return
      end

c ::: -----------------------------------------------------------
c ::: Add fine grid flux to flux register.  Flux array is a fine grid
c ::: edge based object, Register is a coarse grid edge based object.      
c ::: It is assumed that the coarsened flux region contains the register
c ::: region.
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: reg       <=> edge centered coarse grid flux register
c ::: DIMS(reg)  => index limits for reg
c ::: flx        => edge centered fine grid flux array
c ::: DIMS(flx)  => index limits for flx
c ::: numcomp    => number of components to update
c ::: dir        => direction normal to flux register
c ::: ratio(3)   => refinement ratios between coarse and fine
c ::: mult       => scalar multiplicative factor      
c ::: -----------------------------------------------------------

      subroutine FORT_FRFINEADD(reg,DIMS(reg),flx,DIMS(flx),
     &                          numcomp,dir,ratio,mult)
      integer    DIMDEC(reg)
      integer    DIMDEC(flx)
      integer    ratio(3), dir, numcomp
      REAL_T     mult
      REAL_T     reg(DIMV(reg),numcomp)
      REAL_T     flx(DIMV(flx),numcomp)
      
      integer    n, i, j, k, ic, jc, kc, ioff, joff, koff
      integer    ratiox, ratioy, ratioz

      ratiox = ratio(1)
      ratioy = ratio(2)
      ratioz = ratio(3)
      
      if (dir .eq. 0) then
c        ::::: flux normal to X direction
         ic = ARG_L1(reg)
         i = ic*ratiox
         if (ARG_L1(reg) .ne. ARG_H1(reg)) then
            write(6,*) "FORT_FRFINEADD: bad register direction"
            stop
         endif
         if (i .lt. ARG_L1(flx) .or. i .gt. ARG_H1(flx)) then
            write(6,*) "FORT_FRFINEADD: index outside flux range"
            stop
         endif
         do n = 1, numcomp
         do koff = 0, ratioz-1
         do kc = ARG_L3(reg), ARG_H3(reg)
            k = ratioz*kc + koff
            do joff = 0, ratioy-1            
            do jc = ARG_L2(reg), ARG_H2(reg)
               j = ratioy*jc + joff
               reg(ic,jc,kc,n) = reg(ic,jc,kc,n) + mult*flx(i,j,k,n)
            enddo
            enddo
         enddo
         enddo
         enddo
      else if (dir .eq. 1) then
c        ::::: flux normal to Y direction
         jc = ARG_L2(reg)
         j = jc*ratioy
         if (ARG_L2(reg) .ne. ARG_H2(reg)) then
            write(6,*) "FORT_FRFINEADD: bad register direction"
            stop
         endif
         if (j .lt. ARG_L2(flx) .or. j .gt. ARG_H2(flx)) then
            write(6,*) "FORT_FRFINEADD: index outside flux range"
            stop
         endif
         do n = 1, numcomp
         do koff = 0, ratioz-1
         do kc = ARG_L3(reg), ARG_H3(reg)
            k = ratioz*kc + koff
            do ioff = 0, ratiox-1            
            do ic = ARG_L1(reg), ARG_H1(reg)
               i = ratiox*ic + ioff
               reg(ic,jc,kc,n) = reg(ic,jc,kc,n) + mult*flx(i,j,k,n)
            enddo
            enddo
         enddo
         enddo
         enddo
      else
c        ::::: flux normal to Z direction
         kc = ARG_L3(reg)
         k = kc*ratioz
         if (ARG_L3(reg) .ne. ARG_H3(reg)) then
            write(6,*) "FORT_FRFINEADD: bad register direction"
            stop
         endif
         if (k .lt. ARG_L3(flx) .or. k .gt. ARG_H3(flx)) then
            write(6,*) "FORT_FRFINEADD: index outside flux range"
            stop
         endif
         do n = 1, numcomp
         do joff = 0, ratioy-1
         do jc = ARG_L2(reg), ARG_H2(reg)
            j = ratioy*jc + joff
            do ioff = 0, ratiox-1            
            do ic = ARG_L1(reg), ARG_H1(reg)
               i = ratiox*ic + ioff
               reg(ic,jc,kc,n) = reg(ic,jc,kc,n) + mult*flx(i,j,k,n)
            enddo
            enddo
         enddo
         enddo
         enddo
      endif
      
      return
      end

c ::: -----------------------------------------------------------
c ::: Add fine grid flux times area to flux register.  
c ::: Flux array is a fine grid edge based object, Register is a 
c ::: coarse grid edge based object.      
c ::: It is assumed that the coarsened flux region contains the register
c ::: region.
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: reg       <=> edge centered coarse grid flux register
c ::: DIMS(reg)  => index limits for reg
c ::: flx        => edge centered fine grid flux array
c ::: DIMS(flx)  => index limits for flx
c ::: area       => edge centered area array
c ::: DIMS(area) => index limits for area
c ::: numcomp    => number of components to update
c ::: dir        => direction normal to flux register
c ::: ratio(3)   => refinement ratios between coarse and fine
c ::: mult       => scalar multiplicative factor      
c ::: -----------------------------------------------------------

      subroutine FORT_FRFAADD(reg,DIMS(reg),flx,DIMS(flx),area,DIMS(area),
     &                        numcomp,dir,ratio,mult)
      integer    DIMDEC(reg)
      integer    DIMDEC(flx)
      integer    DIMDEC(area)
      integer    ratio(3), dir, numcomp
      REAL_T     mult
      REAL_T     reg(DIMV(reg),numcomp)
      REAL_T     flx(DIMV(flx),numcomp)
      REAL_T     area(DIMV(area))
      
      integer    n, i, j, k, ic, jc, kc, ioff, joff, koff
      integer    ratiox, ratioy, ratioz

      ratiox = ratio(1)
      ratioy = ratio(2)
      ratioz = ratio(3)
      
      if (dir .eq. 0) then
c        ::::: flux normal to X direction
         ic = ARG_L1(reg)
         i = ic*ratiox
         if (ARG_L1(reg) .ne. ARG_H1(reg)) then
            write(6,*) "FORT_FRFINEADD: bad register direction"
            stop
         endif
         if (i .lt. ARG_L1(flx) .or. i .gt. ARG_H1(flx)) then
            write(6,*) "FORT_FRFINEADD: index outside flux range"
            stop
         endif
         do n = 1, numcomp
         do koff = 0, ratioz-1
         do kc = ARG_L3(reg), ARG_H3(reg)
            k = ratioz*kc + koff
            do joff = 0, ratioy-1            
            do jc = ARG_L2(reg), ARG_H2(reg)
               j = ratioy*jc + joff
               reg(ic,jc,kc,n) = reg(ic,jc,kc,n) + 
     &                           mult*area(i,j,k)*flx(i,j,k,n)
            enddo
            enddo
         enddo
         enddo
         enddo
      else if (dir .eq. 1) then
c        ::::: flux normal to Y direction
         jc = ARG_L2(reg)
         j = jc*ratioy
         if (ARG_L2(reg) .ne. ARG_H2(reg)) then
            write(6,*) "FORT_FRFINEADD: bad register direction"
            stop
         endif
         if (j .lt. ARG_L2(flx) .or. j .gt. ARG_H2(flx)) then
            write(6,*) "FORT_FRFINEADD: index outside flux range"
            stop
         endif
         do n = 1, numcomp
         do koff = 0, ratioz-1
         do kc = ARG_L3(reg), ARG_H3(reg)
            k = ratioz*kc + koff
            do ioff = 0, ratiox-1            
            do ic = ARG_L1(reg), ARG_H1(reg)
               i = ratiox*ic + ioff
               reg(ic,jc,kc,n) = reg(ic,jc,kc,n) + 
     &                           mult*area(i,j,k)*flx(i,j,k,n)
            enddo
            enddo
         enddo
         enddo
         enddo
      else
c        ::::: flux normal to Z direction
         kc = ARG_L3(reg)
         k = kc*ratioz
         if (ARG_L3(reg) .ne. ARG_H3(reg)) then
            write(6,*) "FORT_FRFINEADD: bad register direction"
            stop
         endif
         if (k .lt. ARG_L3(flx) .or. k .gt. ARG_H3(flx)) then
            write(6,*) "FORT_FRFINEADD: index outside flux range"
            stop
         endif
         do n = 1, numcomp
         do joff = 0, ratioy-1
         do jc = ARG_L2(reg), ARG_H2(reg)
            j = ratioy*jc + joff
            do ioff = 0, ratiox-1            
            do ic = ARG_L1(reg), ARG_H1(reg)
               i = ratiox*ic + ioff
               reg(ic,jc,kc,n) = reg(ic,jc,kc,n) + 
     &                           mult*area(i,j,k)*flx(i,j,k,n)
            enddo
            enddo
         enddo
         enddo
         enddo
      endif
      
      return
      end

c :: 
c :: --------------------------------------------------------------
c :: reflux:   reflux the data on the outer boundary of
c ::           a fine grid.  
c :: 
c :: Inputs/Outputs
c :: s           <=>  state data array
c :: DIMS(s)      =>  index limits of s array
c :: vol          =>  volume array 
c :: DIMS(vol)    =>  index limits of vol array
c :: reg          =>  flux register
c :: DIMS(reg)    =>  index limits of reg array
c :: lo,hi        =>  subregion of s array to be updated
c :: numcomp      =>  number of components to update
c :: mult         =>  multiplative factor (+1 or -1 depending on nomal)
c :: --------------------------------------------------------------
c :: 
      subroutine FORT_FRREFLUX (s,DIMS(s),vol,DIMS(vol),reg,DIMS(reg),
     &                          lo,hi,numcomp,mult)

      integer    DIMDEC(s)
      integer    DIMDEC(vol)
      integer    DIMDEC(reg)
      integer    lo(SDIM), hi(SDIM)
      integer    numcomp
      REAL_T     mult
      REAL_T       s(DIMV(s),numcomp)
      REAL_T     vol(DIMV(vol))
      REAL_T     reg(DIMV(reg),numcomp)
      
      integer n, i, j, k, lenx, leny, lenz, maxlen
      lenx = hi(1) - lo(1) + 1
      leny = hi(2) - lo(2) + 1
      lenz = hi(3) - lo(3) + 1
      maxlen = max(lenx, leny, lenz)
      
      if (lenx .eq. maxlen) then
         do n = 1, numcomp
            do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  s(i,j,k,n) = s(i,j,k,n) + mult*reg(i,j,k,n)/vol(i,j,k)
               enddo
            enddo
            enddo
         enddo
      else if (leny .eq. maxlen) then
         do n = 1, numcomp
            do i = lo(1), hi(1)
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  s(i,j,k,n) = s(i,j,k,n) + mult*reg(i,j,k,n)/vol(i,j,k)
               enddo
            enddo
            enddo
         enddo
      else
         do n = 1, numcomp
            do i = lo(1), hi(1)
            do j = lo(2), hi(2)
               do k = lo(3), hi(3)
                  s(i,j,k,n) = s(i,j,k,n) + mult*reg(i,j,k,n)/vol(i,j,k)
               enddo
            enddo
            enddo
         enddo
      endif
      
      return
      end

c :: 
c :: --------------------------------------------------------------
c :: cvreflux:   constant volume version of reflux
c :: 
c :: Inputs/Outputs
c :: s           <=>  state data array
c :: DIMS(s)      =>  index limits of s array
c :: dx           =>  cell size
c :: reg          =>  flux register
c :: DIMS(reg)    =>  index limits of reg array
c :: lo,hi        =>  subregion of s array to be updated
c :: numcomp      =>  number of components to update
c :: mult         =>  multiplative factor (+1 or -1 depending on nomal)
c :: --------------------------------------------------------------
c :: 
      subroutine FORT_FRCVREFLUX (s,DIMS(s),dx,reg,DIMS(reg),lo,hi,
     &                            numcomp,mult)

      integer    DIMDEC(s)
      integer    DIMDEC(reg)
      integer    lo(SDIM), hi(SDIM)
      integer    numcomp
      REAL_T     mult, dx(SDIM)
      REAL_T     reg(DIMV(reg),numcomp)
      REAL_T       s(DIMV(s),numcomp)
      
      integer n, i, j, k, lenx, leny, lenz, maxlen
      REAL_T vol

      lenx = hi(1) - lo(1) + 1
      leny = hi(2) - lo(2) + 1
      lenz = hi(3) - lo(3) + 1
      maxlen = max(lenx, leny, lenz)
      vol = dx(1)*dx(2)*dx(3)
      
      if (lenx .eq. maxlen) then
         do n = 1, numcomp
            do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  s(i,j,k,n) = s(i,j,k,n) + mult*reg(i,j,k,n)/vol
               enddo
            enddo
            enddo
         enddo
      else if (leny .eq. maxlen) then
         do n = 1, numcomp
            do i = lo(1), hi(1)
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  s(i,j,k,n) = s(i,j,k,n) + mult*reg(i,j,k,n)/vol
               enddo
            enddo
            enddo
         enddo
      else
         do n = 1, numcomp
            do i = lo(1), hi(1)
            do j = lo(2), hi(2)
               do k = lo(3), hi(3)
                  s(i,j,k,n) = s(i,j,k,n) + mult*reg(i,j,k,n)/vol
               enddo
            enddo
            enddo
         enddo
      endif
      
      return
      end

c ::: -----------------------------------------------------------
c ::: Add flux register into coarse grid flux (opposite of standard reflux ops)
c :::
c ::: INPUTS/OUTPUTS:
c ::: flx        => flux array
c ::: DIMS(flx)  => index limits for flx
c ::: reg       <=  flux register
c ::: DIMS(reg)  => index limits for reg
c ::: lo,hi      => region of update
c ::: mult       => scalar multiplicative factor
c ::: -----------------------------------------------------------

      subroutine FORT_SCALADDTO (flx,DIMS(flx),area,DIMS(area),reg,DIMS(reg),
     $                           lo,hi,numcomp,mult)

      integer    DIMDEC(flx)
      integer    DIMDEC(area)
      integer    DIMDEC(reg)
      integer    lo(SDIM), hi(SDIM)
      integer    numcomp
      REAL_T     mult
      REAL_T      flx(DIMV(flx),numcomp)
      REAL_T     area(DIMV(area))
      REAL_T      reg(DIMV(reg),numcomp)

      integer    n, i, j, k, istart

      if (area(lo(1),lo(2),lo(3)) .lt. 1.d-8) then
        do n = 1, numcomp
        do k = lo(3), hi(3)
        do j = lo(2), hi(2)
          flx(lo(1),j,k,n) = flx(lo(1),j,k,n)
        enddo
        enddo
        enddo
        istart = lo(1)+1
      else
        istart = lo(1)
      endif

      do n = 1, numcomp
       do k = lo(3), hi(3)
       do j = lo(2), hi(2)
        do i = istart, hi(1)
          flx(i,j,k,n) = flx(i,j,k,n) + mult*reg(i,j,k,numcomp)/area(i,j,k)
        enddo
       enddo
       enddo
      enddo

      return
      end

