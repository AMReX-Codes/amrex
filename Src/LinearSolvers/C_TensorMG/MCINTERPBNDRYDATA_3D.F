#include "BC_TYPES.H"
#include "REAL.H"
#include "CONSTANTS.H"
#include "MCINTERPBNDRYDATA_F.H"
#include "ArrayLim.H"

#define SDIM 3
#define NUMDERIV 9
#define DX    1
#define DY    2
#define DZ    3
#define D2X   4
#define D2Y   5
#define D2Z   6
#define DXY   7
#define DXZ   8
#define DYZ   9

c ::: CX,CY,CZ are the "component" number for storing derivative information
c ::: see the assignments to bdry at the bottom of each of the subroutines.
c ::: This should be more suggestive of how info is being stored than
c ::: writing plain numbers.
#define CX   1
#define CY   2
#define CZ   3




c ---------------------------------------------------------------
c ::  FORT_BDINTERPXLO : Interpolation on Xlo Face
c ::       Quadratic Interpolation from crse data
c ::       in directions transverse to face of grid
c ::
c ::  Inputs/Outputs:
c ::  bdry       <=  fine grid bndry data strip
c ::  DIMS(bdry)  => index limits of bdry
c ::  lo,hi       => index limits of grd interior
c ::  DIMS(cb)    => index limits of coarsened grid interior
c ::  nvar        => number of variables in solution, half number in bdry
c ::  ratio       => refinement ratio
c ::  not_covered => mask is set to this value if cell is not
c ::                 covered by another fine grid and not outside the domain.
c ::  mask        => fine grid mask bndry strip
c ::  DIMS(mask)  => index limits of mask array
c ::  crse        => crse grid bndry data strip
c ::  DIMS(crse)  => index limits of crse array
c ::  derives     => crse grid tmp array for derivatives
c ---------------------------------------------------------------

      subroutine FORT_BDINTERPXLO (bdry,DIMS(bdry),
     &           lo,hi,DIMS(cb),nvar,ratio,not_covered,
     &           mask,DIMS(mask),crse,DIMS(crse),derives,hfine)

      integer  nvar, ratio, not_covered
      integer  lo(SDIM), hi(SDIM)
      integer  DIMDEC(bdry)
      integer  DIMDEC(cb)
      integer  DIMDEC(mask)
      integer  DIMDEC(crse)
      REAL_T   bdry(DIMV(bdry),(1+3)*nvar)
      REAL_T   derives(DIM23(cb),NUMDERIV)
      integer  mask(DIMV(mask))
      REAL_T   crse(DIMV(crse),nvar)
      REAL_T   hfine(BL_SPACEDIM-1)

      REAL_T   yy, zz
      integer  i, j, k, ic, jc, kc, joff, koff, n
      integer  jclo, jchi, kclo, kchi
      REAL_T   hc(3),hf(3)

c ::: ------------ functions ------------------------------
      logical notUsable,Usable
      REAL_T ddxcen,ddx2cen,ddxplus,ddxminus
      REAL_T ddycen,ddy2cen,ddyplus,ddyminus
      REAL_T ddzcen,ddz2cen,ddzplus,ddzminus
      REAL_T ddxdycen,ddxdypp,ddxdypm,ddxdymp,ddxdymm
      REAL_T ddxdzcen,ddxdzpp,ddxdzpm,ddxdzmp,ddxdzmm
      REAL_T ddydzcen,ddydzpp,ddydzpm,ddydzmp,ddydzmm
      logical okxypp,okxypm,okxymp,okxymm,okxycen
      logical okxzpp,okxzpm,okxzmp,okxzmm,okxzcen
      logical okyzpp,okyzpm,okyzmp,okyzmm,okyzcen
      REAL_T ddxlowp,ddxlowm
      REAL_T ddylowp,ddylowm
      REAL_T ddzlowp,ddzlowm

      notUsable(i,j,k) = (mask(i,j,k).ne.not_covered)
      Usable(i,j,k) = (mask(i,j,k).eq.not_covered)

      ddxcen(i,j,k,n) = (crse(i+1,j,k,n)-crse(i-1,j,k,n))/(2*hc(1))
      ddycen(i,j,k,n) = (crse(i,j+1,k,n)-crse(i,j-1,k,n))/(2*hc(2))
      ddzcen(i,j,k,n) = (crse(i,j,k+1,n)-crse(i,j,k-1,n))/(2*hc(3))

      ddx2cen(i,j,k,n)=(crse(i+1,j,k,n)-2*crse(i,j,k,n)+crse(i-1,j,k,n))
     &     /(hc(1)**2)
      ddy2cen(i,j,k,n)=(crse(i,j+1,k,n)-2*crse(i,j,k,n)+crse(i,j-1,k,n))
     &     /(hc(2)**2)
      ddz2cen(i,j,k,n)=(crse(i,j,k+1,n)-2*crse(i,j,k,n)+crse(i,j,k-1,n))
     &     /(hc(3)**2)

      ddxplus(i,j,k,n)=(-crse(i+2,j,k,n)+4*crse(i+1,j,k,n)-3*crse(i,j,k,n))
     &     /(2*hc(1))
      ddyplus(i,j,k,n)=(-crse(i,j+2,k,n)+4*crse(i,j+1,k,n)-3*crse(i,j,k,n))
     &     /(2*hc(2))
      ddzplus(i,j,k,n)=(-crse(i,j,k+2,n)+4*crse(i,j,k+1,n)-3*crse(i,j,k,n))
     &     /(2*hc(3))

      ddxminus(i,j,k,n)=(+crse(i-2,j,k,n)-4*crse(i-1,j,k,n)+3*crse(i,j,k,n))
     &     /(2*hc(1))
      ddyminus(i,j,k,n)=(+crse(i,j-2,k,n)-4*crse(i,j-1,k,n)+3*crse(i,j,k,n))
     &     /(2*hc(2))
      ddzminus(i,j,k,n)=(+crse(i,j,k-2,n)-4*crse(i,j,k-1,n)+3*crse(i,j,k,n))
     &     /(2*hc(3))

      ddxdycen(i,j,k,n)=(crse(i+1,j+1,k,n)-crse(i-1,j+1,k,n)
     &     -crse(i+1,j-1,k,n)+crse(i-1,j-1,k,n))/(4*hc(1)*hc(2))
      okxycen(i,j,k)=(Usable(i+1,j+1,k).and.Usable(i-1,j+1,k)
     &     .and.Usable(i+1,j-1,k).and.Usable(i-1,j-1,k))
      ddxdypp(i,j,k,n)=(crse(i+1,j+1,k,n)-crse(i,j+1,k,n)
     &     -crse(i+1,j,k,n)+crse(i,j,k,n))/(hc(1)*hc(2))
      okxypp(i,j,k)=(Usable(i+1,j+1,k).and.Usable(i,j+1,k)
     &     .and.Usable(i+1,j,k))
      ddxdypm(i,j,k,n)=(crse(i+1,j,k,n)-crse(i,j,k,n)
     &     -crse(i+1,j-1,k,n)+crse(i,j-1,k,n))/(4*hc(1)*hc(2))
      okxypm(i,j,k)=(Usable(i+1,j,k)
     &     .and.Usable(i+1,j-1,k).and.Usable(i,j-1,k))
      ddxdymp(i,j,k,n)=(crse(i,j+1,k,n)-crse(i-1,j+1,k,n)
     &     -crse(i,j,k,n)+crse(i-1,j,k,n))/(hc(1)*hc(2))
      okxymp(i,j,k)=(Usable(i,j+1,k).and.Usable(i-1,j+1,k)
     &     .and.Usable(i-1,j,k))
      ddxdymm(i,j,k,n)=(crse(i,j,k,n)-crse(i-1,j,k,n)
     &     -crse(i,j-1,k,n)+crse(i-1,j-1,k,n))/(hc(1)*hc(2))
      okxymm(i,j,k)=(Usable(i-1,j,k)
     &     .and.Usable(i,j-1,k).and.Usable(i-1,j-1,k))

      ddxdzcen(i,j,k,n)=(crse(i+1,j,k+1,n)-crse(i-1,j,k+1,n)
     &     -crse(i+1,j,k-1,n)+crse(i-1,j,k-1,n))/(4*hc(1)*hc(3))
      okxzcen(i,j,k)=(Usable(i+1,j,k+1).and.Usable(i-1,j,k+1)
     &     .and.Usable(i+1,j,k-1).and.Usable(i-1,j,k-1))
      ddxdzpp(i,j,k,n)=(crse(i+1,j,k+1,n)-crse(i,j,k+1,n)
     &     -crse(i+1,j,k,n)+crse(i,j,k,n))/(hc(1)*hc(3))
      okxzpp(i,j,k)=(Usable(i+1,j,k+1).and.Usable(i,j,k+1)
     &     .and.Usable(i+1,j,k))
      ddxdzpm(i,j,k,n)=(crse(i+1,j,k,n)-crse(i,j,k,n)
     &     -crse(i+1,j,k-1,n)+crse(i,j,k-1,n))/(4*hc(1)*hc(3))
      okxzpm(i,j,k)=(Usable(i+1,j,k)
     &     .and.Usable(i+1,j,k-1).and.Usable(i,j,k-1))
      ddxdzmp(i,j,k,n)=(crse(i,j,k+1,n)-crse(i-1,j,k+1,n)
     &     -crse(i,j,k,n)+crse(i-1,j,k,n))/(hc(1)*hc(3))
      okxzmp(i,j,k)=(Usable(i,j,k+1).and.Usable(i-1,j,k+1)
     &     .and.Usable(i-1,j,k))
      ddxdzmm(i,j,k,n)=(crse(i,j,k,n)-crse(i-1,j,k,n)
     &     -crse(i,j,k-1,n)+crse(i-1,j,k-1,n))/(hc(1)*hc(3))
      okxzmm(i,j,k)=(Usable(i-1,j,k)
     &     .and.Usable(i,j,k-1).and.Usable(i-1,j,k-1))

      ddydzcen(i,j,k,n)=(crse(i,j+1,k+1,n)-crse(i,j-1,k+1,n)
     &     -crse(i,j+1,k-1,n)+crse(i,j-1,k-1,n))/(4*hc(2)*hc(3))
      okyzcen(i,j,k)=(Usable(i,j+1,k+1).and.Usable(i,j-1,k+1)
     &     .and.Usable(i,j+1,k-1).and.Usable(i,j-1,k-1))
      ddydzpp(i,j,k,n)=(crse(i,j+1,k+1,n)-crse(i,j,k+1,n)
     &     -crse(i,j+1,k,n)+crse(i,j,k,n))/(hc(2)*hc(3))
      okyzpp(i,j,k)=(Usable(i,j+1,k+1).and.Usable(i,j,k+1)
     &     .and.Usable(i,j+1,k))
      ddydzpm(i,j,k,n)=(crse(i,j+1,k,n)-crse(i,j,k,n)
     &     -crse(i,j+1,k-1,n)+crse(i,j,k-1,n))/(4*hc(2)*hc(3))
      okyzpm(i,j,k)=(Usable(i,j+1,k)
     &     .and.Usable(i,j+1,k-1).and.Usable(i,j,k-1))
      ddydzmp(i,j,k,n)=(crse(i,j,k+1,n)-crse(i,j-1,k+1,n)
     &     -crse(i,j,k,n)+crse(i,j-1,k,n))/(hc(2)*hc(3))
      okyzmp(i,j,k)=(Usable(i,j,k+1).and.Usable(i,j-1,k+1)
     &     .and.Usable(i,j-1,k))
      ddydzmm(i,j,k,n)=(crse(i,j,k,n)-crse(i,j-1,k,n)
     &     -crse(i,j,k-1,n)+crse(i,j-1,k-1,n))/(hc(2)*hc(3))
      okyzmm(i,j,k)=(Usable(i,j-1,k)
     &     .and.Usable(i,j,k-1).and.Usable(i,j-1,k-1))

      ddxlowp(i,j,k,n)=(crse(i+1,j,k,n)-crse(i,j,k,n))/hc(1)
      ddylowp(i,j,k,n)=(crse(i,j+1,k,n)-crse(i,j,k,n))/hc(2)
      ddzlowp(i,j,k,n)=(crse(i,j,k+1,n)-crse(i,j,k,n))/hc(3)

      ddxlowm(i,j,k,n)=(-crse(i-1,j,k,n)+crse(i,j,k,n))/hc(1)
      ddylowm(i,j,k,n)=(-crse(i,j-1,k,n)+crse(i,j,k,n))/hc(2)
      ddzlowm(i,j,k,n)=(-crse(i,j,k-1,n)+crse(i,j,k,n))/hc(3)
c ::: ------------ functions ------------------------------

      kclo = ARG_L3(cb)
      kchi = ARG_H3(cb)
      jclo = ARG_L2(cb)
      jchi = ARG_H2(cb)
      ic   = ARG_L1(cb)-1
      i    = lo(1)-1
c ::: :::
      hf(1) = 0
      hf(2) = hfine(1)
      hf(3) = hfine(2)
      hc(1) = 0
      hc(2) = hfine(1)*ratio
      hc(3) = hfine(2)*ratio
c ::: :::

      do n = 1, nvar
c     ::::: calculate derivatives for interpolation
         do kc = kclo, kchi
            k = ratio*kc
            do jc = jclo, jchi
               j = ratio*jc

               derives(jc,kc,DY) = ddycen(ic,jc,kc,n)
               derives(jc,kc,D2Y) = ddy2cen(ic,jc,kc,n)

               if( notUsable(i,j-1,k) .and. notUsable(i,j+ratio,k)  ) then
                  derives(jc,kc,DY)  = zero
                  derives(jc,kc,D2Y) = zero
               else if( notUsable(i,j-1,k) .and. Usable(i,j+ratio,k) ) then
c     :::: :::: : positive side usable
                  if( jc+2.le.jchi+1 ) then
                     if( Usable(i,j+2*ratio,k) ) then
                        derives(jc,kc,DY)=ddyplus(ic,jc,kc,n)
                        derives(jc,kc,D2Y)=ddy2cen(ic,jc+1,kc,n)
                     else
                        derives(jc,kc,DY)=ddylowp(ic,jc,kc,n)
                        derives(jc,kc,D2Y)=zero
                     endif
                  else
                     derives(jc,kc,DY)=ddylowp(ic,jc,kc,n)
                     derives(jc,kc,D2Y)=zero
                  endif
               else if( Usable(i,j-1,k) .and. notUsable(i,j+ratio,k) ) then
c     :::: :::: : negative side usable
                  if( jc-2.ge.jclo-1 ) then
                     if( Usable(i,j-1,k) ) then
                        derives(jc,kc,DY)=ddyminus(ic,jc,kc,n)
                        derives(jc,kc,D2Y)=ddy2cen(ic,jc-1,kc,n)
                     else
                        derives(jc,kc,DY)=ddylowm(ic,jc,kc,n)
                        derives(jc,kc,D2Y)=zero
                     endif
                  else
                     derives(jc,kc,DY)=ddylowm(ic,jc,kc,n)
                     derives(jc,kc,D2Y)=zero
                  endif
               endif

               derives(jc,kc,DZ)  = ddzcen(ic,jc,kc,n)
               derives(jc,kc,D2Z) = ddz2cen(ic,jc,kc,n)

               if( notUsable(i,j,k-1) .and. notUsable(i,j,k+ratio)) then
                  derives(jc,kc,DZ)  = zero
                  derives(jc,kc,D2Z) = zero
               else if( notUsable(i,j,k-1) .and. Usable(i,j,k+ratio)) then
c     :::: :::: : positive size usable
                  if( kc+2.le.kchi+1 ) then
                     if( Usable(i,j,k+2*ratio) ) then
                        derives(jc,kc,DZ)  = ddzplus(ic,jc,kc,n)
                        derives(jc,kc,D2Z) =  ddz2cen(ic,jc,kc+1,n)
                     else
                        derives(jc,kc,DZ)  = ddzlowp(ic,jc,kc,n)
                        derives(jc,kc,D2Z) =  zero
                     endif
                  else
                     derives(jc,kc,DZ)  = ddzlowp(ic,jc,kc,n)
                     derives(jc,kc,D2Z) =  zero
                  endif
               else if( Usable(i,j,k-1) .and. notUsable(i,j,k+ratio)) then
c     :::: :::: : negative size usable
                  if( kc-2.ge.kclo-1 ) then
                     if( Usable(i,j,k-1) ) then
                        derives(jc,kc,DZ)  = ddzminus(ic,jc,kc,n)
                        derives(jc,kc,D2Z) =  ddz2cen(ic,jc,kc-1,n)
                     else
                        derives(jc,kc,DZ)  = ddzlowm(ic,jc,kc,n)
                        derives(jc,kc,D2Z) =  zero
                     endif
                  else
                     derives(jc,kc,DZ)  = ddzlowm(ic,jc,kc,n)
                     derives(jc,kc,D2Z) =  zero
                  endif
               endif

               if( okyzcen(i,j,k) ) then
                  derives(jc,kc,DYZ) = ddydzcen(ic,jc,kc,n)
               else if( okyzpp(i,j,k) ) then
                  derives(jc,kc,DYZ) = ddydzpp(ic,jc,kc,n)
               else if( okyzpm(i,j,k) ) then
                  derives(jc,kc,DYZ) = ddydzpm(ic,jc,kc,n)
               else if( okyzmp(i,j,k) ) then
                  derives(jc,kc,DYZ) = ddydzmp(ic,jc,kc,n)
               else if( okyzmm(i,j,k) ) then
                  derives(jc,kc,DYZ) = ddydzmm(ic,jc,kc,n)
               else
                  derives(jc,kc,DYZ) = zero
               endif

            enddo
         enddo

c     ::::: interpolate to fine grid
         do koff = 0, ratio - 1
            zz = (koff+half)*hf(3)-hc(3)/2
            do kc = kclo,kchi
               k = ratio*kc + koff
               do joff = 0, ratio - 1
                  yy = (joff+half)*hf(2)-hc(2)/2
                  do jc = jclo, jchi
                     j = ratio*jc + joff
                     bdry(i,j,k,n) = crse(ic,jc,kc,n) + 
     &                    yy*derives(jc,kc,DY) +
     &                    half*yy**2*derives(jc,kc,D2Y) + 
     &                    zz*derives(jc,kc,DZ) +
     &                    half*zz**2*derives(jc,kc,D2Z) + 
     &                    yy*zz*derives(jc,kc,DYZ) 
                     bdry(i,j,k,n+CY*nvar) =
     &                    derives(jc,kc,DY)+
     &                    yy*derives(jc,kc,D2Y)+
     &                    zz*derives(jc,kc,DYZ)
                     bdry(i,j,k,n+CZ*nvar) =
     &                    derives(jc,kc,DZ)+
     &                    zz*derives(jc,kc,D2Z)+
     &                    yy*derives(jc,kc,DYZ)
                  enddo
               enddo
            enddo
         enddo

      enddo
            
      return
      end


c ---------------------------------------------------------------
c ::  FORT_BDINTERPXHI : Interpolation on Xhi Face
c ::       Quadratic Interpolation from crse data
c ::       in directions transverse to face of grid
c ::
c ::  Inputs/Outputs:
c ::  bdry       <=  fine grid bndry data strip
c ::  DIMS(bdry)  => index limits of bdry
c ::  lo,hi       => index limits of grd interior
c ::  DIMS(cb)    => index limits of coarsened grid interior
c ::  nvar        => number of variables in solution, half number in bdry
c ::  ratio       => refinement ratio
c ::  not_covered => mask is set to this value if cell is not
c ::                 covered by another fine grid and not outside the domain.
c ::  mask        => fine grid mask bndry strip
c ::  DIMS(mask)  => index limits of mask array
c ::  crse        => crse grid bndry data strip
c ::  DIMS(crse)  => index limits of crse array
c ::  derives     => crse grid tmp array for derivatives
c ---------------------------------------------------------------

      subroutine FORT_BDINTERPXHI (bdry,DIMS(bdry),
     &           lo,hi,DIMS(cb),nvar,ratio,not_covered,
     &           mask,DIMS(mask),crse,DIMS(crse),derives,hfine)

      integer  nvar, ratio, not_covered
      integer  lo(SDIM), hi(SDIM)
      integer  DIMDEC(bdry)
      integer  DIMDEC(cb)
      integer  DIMDEC(mask)
      integer  DIMDEC(crse)
      REAL_T   bdry(DIMV(bdry),(1+3)*nvar)
      REAL_T   derives(DIM23(cb),NUMDERIV)
      integer  mask(DIMV(mask))
      REAL_T   crse(DIMV(crse),nvar)
      REAL_T   hfine(BL_SPACEDIM-1)

      REAL_T   yy, zz
      integer  i, j, k, ic, jc, kc, joff, koff, n
      integer  jclo, jchi, kclo, kchi
      REAL_T   hc(3),hf(3)

c ::: ------------ functions ------------------------------
      logical notUsable,Usable
      REAL_T ddxcen,ddx2cen,ddxplus,ddxminus
      REAL_T ddycen,ddy2cen,ddyplus,ddyminus
      REAL_T ddzcen,ddz2cen,ddzplus,ddzminus
      REAL_T ddxdycen,ddxdypp,ddxdypm,ddxdymp,ddxdymm
      REAL_T ddxdzcen,ddxdzpp,ddxdzpm,ddxdzmp,ddxdzmm
      REAL_T ddydzcen,ddydzpp,ddydzpm,ddydzmp,ddydzmm
      logical okxypp,okxypm,okxymp,okxymm,okxycen
      logical okxzpp,okxzpm,okxzmp,okxzmm,okxzcen
      logical okyzpp,okyzpm,okyzmp,okyzmm,okyzcen
      REAL_T ddxlowp,ddxlowm
      REAL_T ddylowp,ddylowm
      REAL_T ddzlowp,ddzlowm

      notUsable(i,j,k) = (mask(i,j,k).ne.not_covered)
      Usable(i,j,k) = (mask(i,j,k).eq.not_covered)

      ddxcen(i,j,k,n) = (crse(i+1,j,k,n)-crse(i-1,j,k,n))/(2*hc(1))
      ddycen(i,j,k,n) = (crse(i,j+1,k,n)-crse(i,j-1,k,n))/(2*hc(2))
      ddzcen(i,j,k,n) = (crse(i,j,k+1,n)-crse(i,j,k-1,n))/(2*hc(3))

      ddx2cen(i,j,k,n)=(crse(i+1,j,k,n)-2*crse(i,j,k,n)+crse(i-1,j,k,n))
     &     /(hc(1)**2)
      ddy2cen(i,j,k,n)=(crse(i,j+1,k,n)-2*crse(i,j,k,n)+crse(i,j-1,k,n))
     &     /(hc(2)**2)
      ddz2cen(i,j,k,n)=(crse(i,j,k+1,n)-2*crse(i,j,k,n)+crse(i,j,k-1,n))
     &     /(hc(3)**2)

      ddxplus(i,j,k,n)=(-crse(i+2,j,k,n)+4*crse(i+1,j,k,n)-3*crse(i,j,k,n))
     &     /(2*hc(1))
      ddyplus(i,j,k,n)=(-crse(i,j+2,k,n)+4*crse(i,j+1,k,n)-3*crse(i,j,k,n))
     &     /(2*hc(2))
      ddzplus(i,j,k,n)=(-crse(i,j,k+2,n)+4*crse(i,j,k+1,n)-3*crse(i,j,k,n))
     &     /(2*hc(3))

      ddxminus(i,j,k,n)=(+crse(i-2,j,k,n)-4*crse(i-1,j,k,n)+3*crse(i,j,k,n))
     &     /(2*hc(1))
      ddyminus(i,j,k,n)=(+crse(i,j-2,k,n)-4*crse(i,j-1,k,n)+3*crse(i,j,k,n))
     &     /(2*hc(2))
      ddzminus(i,j,k,n)=(+crse(i,j,k-2,n)-4*crse(i,j,k-1,n)+3*crse(i,j,k,n))
     &     /(2*hc(3))

      ddxdycen(i,j,k,n)=(crse(i+1,j+1,k,n)-crse(i-1,j+1,k,n)
     &     -crse(i+1,j-1,k,n)+crse(i-1,j-1,k,n))/(4*hc(1)*hc(2))
      okxycen(i,j,k)=(Usable(i+1,j+1,k).and.Usable(i-1,j+1,k)
     &     .and.Usable(i+1,j-1,k).and.Usable(i-1,j-1,k))
      ddxdypp(i,j,k,n)=(crse(i+1,j+1,k,n)-crse(i,j+1,k,n)
     &     -crse(i+1,j,k,n)+crse(i,j,k,n))/(hc(1)*hc(2))
      okxypp(i,j,k)=(Usable(i+1,j+1,k).and.Usable(i,j+1,k)
     &     .and.Usable(i+1,j,k))
      ddxdypm(i,j,k,n)=(crse(i+1,j,k,n)-crse(i,j,k,n)
     &     -crse(i+1,j-1,k,n)+crse(i,j-1,k,n))/(4*hc(1)*hc(2))
      okxypm(i,j,k)=(Usable(i+1,j,k)
     &     .and.Usable(i+1,j-1,k).and.Usable(i,j-1,k))
      ddxdymp(i,j,k,n)=(crse(i,j+1,k,n)-crse(i-1,j+1,k,n)
     &     -crse(i,j,k,n)+crse(i-1,j,k,n))/(hc(1)*hc(2))
      okxymp(i,j,k)=(Usable(i,j+1,k).and.Usable(i-1,j+1,k)
     &     .and.Usable(i-1,j,k))
      ddxdymm(i,j,k,n)=(crse(i,j,k,n)-crse(i-1,j,k,n)
     &     -crse(i,j-1,k,n)+crse(i-1,j-1,k,n))/(hc(1)*hc(2))
      okxymm(i,j,k)=(Usable(i-1,j,k)
     &     .and.Usable(i,j-1,k).and.Usable(i-1,j-1,k))

      ddxdzcen(i,j,k,n)=(crse(i+1,j,k+1,n)-crse(i-1,j,k+1,n)
     &     -crse(i+1,j,k-1,n)+crse(i-1,j,k-1,n))/(4*hc(1)*hc(3))
      okxzcen(i,j,k)=(Usable(i+1,j,k+1).and.Usable(i-1,j,k+1)
     &     .and.Usable(i+1,j,k-1).and.Usable(i-1,j,k-1))
      ddxdzpp(i,j,k,n)=(crse(i+1,j,k+1,n)-crse(i,j,k+1,n)
     &     -crse(i+1,j,k,n)+crse(i,j,k,n))/(hc(1)*hc(3))
      okxzpp(i,j,k)=(Usable(i+1,j,k+1).and.Usable(i,j,k+1)
     &     .and.Usable(i+1,j,k))
      ddxdzpm(i,j,k,n)=(crse(i+1,j,k,n)-crse(i,j,k,n)
     &     -crse(i+1,j,k-1,n)+crse(i,j,k-1,n))/(4*hc(1)*hc(3))
      okxzpm(i,j,k)=(Usable(i+1,j,k)
     &     .and.Usable(i+1,j,k-1).and.Usable(i,j,k-1))
      ddxdzmp(i,j,k,n)=(crse(i,j,k+1,n)-crse(i-1,j,k+1,n)
     &     -crse(i,j,k,n)+crse(i-1,j,k,n))/(hc(1)*hc(3))
      okxzmp(i,j,k)=(Usable(i,j,k+1).and.Usable(i-1,j,k+1)
     &     .and.Usable(i-1,j,k))
      ddxdzmm(i,j,k,n)=(crse(i,j,k,n)-crse(i-1,j,k,n)
     &     -crse(i,j,k-1,n)+crse(i-1,j,k-1,n))/(hc(1)*hc(3))
      okxzmm(i,j,k)=(Usable(i-1,j,k)
     &     .and.Usable(i,j,k-1).and.Usable(i-1,j,k-1))

      ddydzcen(i,j,k,n)=(crse(i,j+1,k+1,n)-crse(i,j-1,k+1,n)
     &     -crse(i,j+1,k-1,n)+crse(i,j-1,k-1,n))/(4*hc(2)*hc(3))
      okyzcen(i,j,k)=(Usable(i,j+1,k+1).and.Usable(i,j-1,k+1)
     &     .and.Usable(i,j+1,k-1).and.Usable(i,j-1,k-1))
      ddydzpp(i,j,k,n)=(crse(i,j+1,k+1,n)-crse(i,j,k+1,n)
     &     -crse(i,j+1,k,n)+crse(i,j,k,n))/(hc(2)*hc(3))
      okyzpp(i,j,k)=(Usable(i,j+1,k+1).and.Usable(i,j,k+1)
     &     .and.Usable(i,j+1,k))
      ddydzpm(i,j,k,n)=(crse(i,j+1,k,n)-crse(i,j,k,n)
     &     -crse(i,j+1,k-1,n)+crse(i,j,k-1,n))/(4*hc(2)*hc(3))
      okyzpm(i,j,k)=(Usable(i,j+1,k)
     &     .and.Usable(i,j+1,k-1).and.Usable(i,j,k-1))
      ddydzmp(i,j,k,n)=(crse(i,j,k+1,n)-crse(i,j-1,k+1,n)
     &     -crse(i,j,k,n)+crse(i,j-1,k,n))/(hc(2)*hc(3))
      okyzmp(i,j,k)=(Usable(i,j,k+1).and.Usable(i,j-1,k+1)
     &     .and.Usable(i,j-1,k))
      ddydzmm(i,j,k,n)=(crse(i,j,k,n)-crse(i,j-1,k,n)
     &     -crse(i,j,k-1,n)+crse(i,j-1,k-1,n))/(hc(2)*hc(3))
      okyzmm(i,j,k)=(Usable(i,j-1,k)
     &     .and.Usable(i,j,k-1).and.Usable(i,j-1,k-1))

      ddxlowp(i,j,k,n)=(crse(i+1,j,k,n)-crse(i,j,k,n))/hc(1)
      ddylowp(i,j,k,n)=(crse(i,j+1,k,n)-crse(i,j,k,n))/hc(2)
      ddzlowp(i,j,k,n)=(crse(i,j,k+1,n)-crse(i,j,k,n))/hc(3)

      ddxlowm(i,j,k,n)=(-crse(i-1,j,k,n)+crse(i,j,k,n))/hc(1)
      ddylowm(i,j,k,n)=(-crse(i,j-1,k,n)+crse(i,j,k,n))/hc(2)
      ddzlowm(i,j,k,n)=(-crse(i,j,k-1,n)+crse(i,j,k,n))/hc(3)
c ::: ------------ functions ------------------------------

      kclo = ARG_L3(cb)
      kchi = ARG_H3(cb)
      jclo = ARG_L2(cb)
      jchi = ARG_H2(cb)
      ic   = ARG_H1(cb)+1
      i    = hi(1)+1
c ::: :::
      hf(1) = 0
      hf(2) = hfine(1)
      hf(3) = hfine(2)
      hc(1) = 0
      hc(2) = hfine(1)*ratio
      hc(3) = hfine(2)*ratio
c ::: :::

      do n = 1, nvar
c     ::::: calculate derivatives for interpolation
         do kc = kclo, kchi
            k = ratio*kc
            do jc = jclo, jchi
               j = ratio*jc

               derives(jc,kc,DY) = ddycen(ic,jc,kc,n)
               derives(jc,kc,D2Y) = ddy2cen(ic,jc,kc,n)

               if( notUsable(i,j-1,k) .and. notUsable(i,j+ratio,k)  ) then
                  derives(jc,kc,DY)  = zero
                  derives(jc,kc,D2Y) = zero
               else if( notUsable(i,j-1,k) .and. Usable(i,j+ratio,k) ) then
c     :::: :::: : positive side usable
                  if( jc+2.le.jchi+1 ) then
                     if( Usable(i,j+2*ratio,k) ) then
                        derives(jc,kc,DY)=ddyplus(ic,jc,kc,n)
                        derives(jc,kc,D2Y)=ddy2cen(ic,jc+1,kc,n)
                     else
                        derives(jc,kc,DY)=ddylowp(ic,jc,kc,n)
                        derives(jc,kc,D2Y)=zero
                     endif
                  else
                     derives(jc,kc,DY)=ddylowp(ic,jc,kc,n)
                     derives(jc,kc,D2Y)=zero
                  endif
               else if( Usable(i,j-1,k) .and. notUsable(i,j+ratio,k) ) then
c     :::: :::: : negative side usable
                  if( jc-2.ge.jclo-1 ) then
                     if( Usable(i,j-1,k) ) then
                        derives(jc,kc,DY)=ddyminus(ic,jc,kc,n)
                        derives(jc,kc,D2Y)=ddy2cen(ic,jc-1,kc,n)
                     else
                        derives(jc,kc,DY)=ddylowm(ic,jc,kc,n)
                        derives(jc,kc,D2Y)=zero
                     endif
                  else
                     derives(jc,kc,DY)=ddylowm(ic,jc,kc,n)
                     derives(jc,kc,D2Y)=zero
                  endif
               endif

               derives(jc,kc,DZ)  = ddzcen(ic,jc,kc,n)
               derives(jc,kc,D2Z) = ddz2cen(ic,jc,kc,n)

               if( notUsable(i,j,k-1) .and. notUsable(i,j,k+ratio)) then
                  derives(jc,kc,DZ)  = zero
                  derives(jc,kc,D2Z) = zero
               else if( notUsable(i,j,k-1) .and. Usable(i,j,k+ratio)) then
c     :::: :::: : positive size usable
                  if( kc+2.le.kchi+1 ) then
                     if( Usable(i,j,k+2*ratio) ) then
                        derives(jc,kc,DZ)  = ddzplus(ic,jc,kc,n)
                        derives(jc,kc,D2Z) =  ddz2cen(ic,jc,kc+1,n)
                     else
                        derives(jc,kc,DZ)  = ddzlowp(ic,jc,kc,n)
                        derives(jc,kc,D2Z) =  zero
                     endif
                  else
                     derives(jc,kc,DZ)  = ddzlowp(ic,jc,kc,n)
                     derives(jc,kc,D2Z) =  zero
                  endif
               else if( Usable(i,j,k-1) .and. notUsable(i,j,k+ratio)) then
c     :::: :::: : negative size usable
                  if( kc-2.ge.kclo-1 ) then
                     if( Usable(i,j,k-1) ) then
                        derives(jc,kc,DZ)  = ddzminus(ic,jc,kc,n)
                        derives(jc,kc,D2Z) =  ddz2cen(ic,jc,kc-1,n)
                     else
                        derives(jc,kc,DZ)  = ddzlowm(ic,jc,kc,n)
                        derives(jc,kc,D2Z) =  zero
                     endif
                  else
                     derives(jc,kc,DZ)  = ddzlowm(ic,jc,kc,n)
                     derives(jc,kc,D2Z) =  zero
                  endif
               endif

               if( okyzcen(i,j,k) ) then
                  derives(jc,kc,DYZ) = ddydzcen(ic,jc,kc,n)
               else if( okyzpp(i,j,k) ) then
                  derives(jc,kc,DYZ) = ddydzpp(ic,jc,kc,n)
               else if( okyzpm(i,j,k) ) then
                  derives(jc,kc,DYZ) = ddydzpm(ic,jc,kc,n)
               else if( okyzmp(i,j,k) ) then
                  derives(jc,kc,DYZ) = ddydzmp(ic,jc,kc,n)
               else if( okyzmm(i,j,k) ) then
                  derives(jc,kc,DYZ) = ddydzmm(ic,jc,kc,n)
               else
                  derives(jc,kc,DYZ) = zero
               endif

            enddo
         enddo

c     ::::: interpolate to fine grid
         do koff = 0, ratio - 1
            zz = (koff+half)*hf(3)-hc(3)/2
            do kc = kclo,kchi
               k = ratio*kc + koff
               do joff = 0, ratio - 1
                  yy = (joff+half)*hf(2)-hc(2)/2
                  do jc = jclo, jchi
                     j = ratio*jc + joff
                     bdry(i,j,k,n) = crse(ic,jc,kc,n) + 
     &                    yy*derives(jc,kc,DY) +
     &                    half*yy**2*derives(jc,kc,D2Y) + 
     &                    zz*derives(jc,kc,DZ) +
     &                    half*zz**2*derives(jc,kc,D2Z) + 
     &                    yy*zz*derives(jc,kc,DYZ) 
                     bdry(i,j,k,n+CY*nvar) =
     &                    derives(jc,kc,DY)+
     &                    yy*derives(jc,kc,D2Y)+
     &                    zz*derives(jc,kc,DYZ)
                     bdry(i,j,k,n+CZ*nvar) =
     &                    derives(jc,kc,DZ)+
     &                    zz*derives(jc,kc,D2Z)+
     &                    yy*derives(jc,kc,DYZ)
                  enddo
               enddo
            enddo
         enddo

      enddo
            
      return
      end





c ---------------------------------------------------------------
c ::  FORT_BDINTERPYLO : Interpolation on Ylo Face
c ::       Quadratic Interpolation from crse data
c ::       in directions transverse to face of grid
c ::
c ::  Inputs/Outputs:
c ::  bdry       <=  fine grid bndry data strip
c ::  DIMS(bdry)  => index limits of bdry
c ::  lo,hi       => index limits of grd interior
c ::  DIMS(cb)    => index limits of coarsened grid interior
c ::  nvar        => number of variables in solution, half number in bdry
c ::  ratio       => refinement ratio
c ::  not_covered => mask is set to this value if cell is not
c ::                 covered by another fine grid and not outside the domain.
c ::  mask        => fine grid mask bndry strip
c ::  DIMS(mask)  => index limits of mask array
c ::  crse        => crse grid bndry data strip
c ::  DIMS(crse)  => index limits of crse array
c ::  derives     => crse grid tmp array for derivatives
c ---------------------------------------------------------------

      subroutine FORT_BDINTERPYLO (bdry,DIMS(bdry),
     &           lo,hi,DIMS(cb),nvar,ratio,not_covered,
     &           mask,DIMS(mask),crse,DIMS(crse),derives,hfine)

      integer  nvar, ratio, not_covered
      integer  lo(SDIM), hi(SDIM)
      integer  DIMDEC(bdry)
      integer  DIMDEC(cb)
      integer  DIMDEC(mask)
      integer  DIMDEC(crse)
      REAL_T   bdry(DIMV(bdry),(1+3)*nvar)
      REAL_T   derives(DIM13(cb),NUMDERIV)
      integer  mask(DIMV(mask))
      REAL_T   crse(DIMV(crse),nvar)
      REAL_T   hfine(BL_SPACEDIM-1)

      REAL_T   xx, zz
      integer  i, j, k, ic, jc, kc, ioff, koff, n
      integer  iclo, ichi, kclo, kchi
      REAL_T   hc(3),hf(3)

c ::: ------------ functions ------------------------------
      logical notUsable,Usable
      REAL_T ddxcen,ddx2cen,ddxplus,ddxminus
      REAL_T ddycen,ddy2cen,ddyplus,ddyminus
      REAL_T ddzcen,ddz2cen,ddzplus,ddzminus
      REAL_T ddxdycen,ddxdypp,ddxdypm,ddxdymp,ddxdymm
      REAL_T ddxdzcen,ddxdzpp,ddxdzpm,ddxdzmp,ddxdzmm
      REAL_T ddydzcen,ddydzpp,ddydzpm,ddydzmp,ddydzmm
      logical okxypp,okxypm,okxymp,okxymm,okxycen
      logical okxzpp,okxzpm,okxzmp,okxzmm,okxzcen
      logical okyzpp,okyzpm,okyzmp,okyzmm,okyzcen
      REAL_T ddxlowp,ddxlowm
      REAL_T ddylowp,ddylowm
      REAL_T ddzlowp,ddzlowm

      notUsable(i,j,k) = (mask(i,j,k).ne.not_covered)
      Usable(i,j,k) = (mask(i,j,k).eq.not_covered)

      ddxcen(i,j,k,n) = (crse(i+1,j,k,n)-crse(i-1,j,k,n))/(2*hc(1))
      ddycen(i,j,k,n) = (crse(i,j+1,k,n)-crse(i,j-1,k,n))/(2*hc(2))
      ddzcen(i,j,k,n) = (crse(i,j,k+1,n)-crse(i,j,k-1,n))/(2*hc(3))

      ddx2cen(i,j,k,n)=(crse(i+1,j,k,n)-2*crse(i,j,k,n)+crse(i-1,j,k,n))
     &     /(hc(1)**2)
      ddy2cen(i,j,k,n)=(crse(i,j+1,k,n)-2*crse(i,j,k,n)+crse(i,j-1,k,n))
     &     /(hc(2)**2)
      ddz2cen(i,j,k,n)=(crse(i,j,k+1,n)-2*crse(i,j,k,n)+crse(i,j,k-1,n))
     &     /(hc(3)**2)

      ddxplus(i,j,k,n)=(-crse(i+2,j,k,n)+4*crse(i+1,j,k,n)-3*crse(i,j,k,n))
     &     /(2*hc(1))
      ddyplus(i,j,k,n)=(-crse(i,j+2,k,n)+4*crse(i,j+1,k,n)-3*crse(i,j,k,n))
     &     /(2*hc(2))
      ddzplus(i,j,k,n)=(-crse(i,j,k+2,n)+4*crse(i,j,k+1,n)-3*crse(i,j,k,n))
     &     /(2*hc(3))

      ddxminus(i,j,k,n)=(+crse(i-2,j,k,n)-4*crse(i-1,j,k,n)+3*crse(i,j,k,n))
     &     /(2*hc(1))
      ddyminus(i,j,k,n)=(+crse(i,j-2,k,n)-4*crse(i,j-1,k,n)+3*crse(i,j,k,n))
     &     /(2*hc(2))
      ddzminus(i,j,k,n)=(+crse(i,j,k-2,n)-4*crse(i,j,k-1,n)+3*crse(i,j,k,n))
     &     /(2*hc(3))

      ddxdycen(i,j,k,n)=(crse(i+1,j+1,k,n)-crse(i-1,j+1,k,n)
     &     -crse(i+1,j-1,k,n)+crse(i-1,j-1,k,n))/(4*hc(1)*hc(2))
      okxycen(i,j,k)=(Usable(i+1,j+1,k).and.Usable(i-1,j+1,k)
     &     .and.Usable(i+1,j-1,k).and.Usable(i-1,j-1,k))
      ddxdypp(i,j,k,n)=(crse(i+1,j+1,k,n)-crse(i,j+1,k,n)
     &     -crse(i+1,j,k,n)+crse(i,j,k,n))/(hc(1)*hc(2))
      okxypp(i,j,k)=(Usable(i+1,j+1,k).and.Usable(i,j+1,k)
     &     .and.Usable(i+1,j,k))
      ddxdypm(i,j,k,n)=(crse(i+1,j,k,n)-crse(i,j,k,n)
     &     -crse(i+1,j-1,k,n)+crse(i,j-1,k,n))/(4*hc(1)*hc(2))
      okxypm(i,j,k)=(Usable(i+1,j,k)
     &     .and.Usable(i+1,j-1,k).and.Usable(i,j-1,k))
      ddxdymp(i,j,k,n)=(crse(i,j+1,k,n)-crse(i-1,j+1,k,n)
     &     -crse(i,j,k,n)+crse(i-1,j,k,n))/(hc(1)*hc(2))
      okxymp(i,j,k)=(Usable(i,j+1,k).and.Usable(i-1,j+1,k)
     &     .and.Usable(i-1,j,k))
      ddxdymm(i,j,k,n)=(crse(i,j,k,n)-crse(i-1,j,k,n)
     &     -crse(i,j-1,k,n)+crse(i-1,j-1,k,n))/(hc(1)*hc(2))
      okxymm(i,j,k)=(Usable(i-1,j,k)
     &     .and.Usable(i,j-1,k).and.Usable(i-1,j-1,k))

      ddxdzcen(i,j,k,n)=(crse(i+1,j,k+1,n)-crse(i-1,j,k+1,n)
     &     -crse(i+1,j,k-1,n)+crse(i-1,j,k-1,n))/(4*hc(1)*hc(3))
      okxzcen(i,j,k)=(Usable(i+1,j,k+1).and.Usable(i-1,j,k+1)
     &     .and.Usable(i+1,j,k-1).and.Usable(i-1,j,k-1))
      ddxdzpp(i,j,k,n)=(crse(i+1,j,k+1,n)-crse(i,j,k+1,n)
     &     -crse(i+1,j,k,n)+crse(i,j,k,n))/(hc(1)*hc(3))
      okxzpp(i,j,k)=(Usable(i+1,j,k+1).and.Usable(i,j,k+1)
     &     .and.Usable(i+1,j,k))
      ddxdzpm(i,j,k,n)=(crse(i+1,j,k,n)-crse(i,j,k,n)
     &     -crse(i+1,j,k-1,n)+crse(i,j,k-1,n))/(4*hc(1)*hc(3))
      okxzpm(i,j,k)=(Usable(i+1,j,k)
     &     .and.Usable(i+1,j,k-1).and.Usable(i,j,k-1))
      ddxdzmp(i,j,k,n)=(crse(i,j,k+1,n)-crse(i-1,j,k+1,n)
     &     -crse(i,j,k,n)+crse(i-1,j,k,n))/(hc(1)*hc(3))
      okxzmp(i,j,k)=(Usable(i,j,k+1).and.Usable(i-1,j,k+1)
     &     .and.Usable(i-1,j,k))
      ddxdzmm(i,j,k,n)=(crse(i,j,k,n)-crse(i-1,j,k,n)
     &     -crse(i,j,k-1,n)+crse(i-1,j,k-1,n))/(hc(1)*hc(3))
      okxzmm(i,j,k)=(Usable(i-1,j,k)
     &     .and.Usable(i,j,k-1).and.Usable(i-1,j,k-1))

      ddydzcen(i,j,k,n)=(crse(i,j+1,k+1,n)-crse(i,j-1,k+1,n)
     &     -crse(i,j+1,k-1,n)+crse(i,j-1,k-1,n))/(4*hc(2)*hc(3))
      okyzcen(i,j,k)=(Usable(i,j+1,k+1).and.Usable(i,j-1,k+1)
     &     .and.Usable(i,j+1,k-1).and.Usable(i,j-1,k-1))
      ddydzpp(i,j,k,n)=(crse(i,j+1,k+1,n)-crse(i,j,k+1,n)
     &     -crse(i,j+1,k,n)+crse(i,j,k,n))/(hc(2)*hc(3))
      okyzpp(i,j,k)=(Usable(i,j+1,k+1).and.Usable(i,j,k+1)
     &     .and.Usable(i,j+1,k))
      ddydzpm(i,j,k,n)=(crse(i,j+1,k,n)-crse(i,j,k,n)
     &     -crse(i,j+1,k-1,n)+crse(i,j,k-1,n))/(4*hc(2)*hc(3))
      okyzpm(i,j,k)=(Usable(i,j+1,k)
     &     .and.Usable(i,j+1,k-1).and.Usable(i,j,k-1))
      ddydzmp(i,j,k,n)=(crse(i,j,k+1,n)-crse(i,j-1,k+1,n)
     &     -crse(i,j,k,n)+crse(i,j-1,k,n))/(hc(2)*hc(3))
      okyzmp(i,j,k)=(Usable(i,j,k+1).and.Usable(i,j-1,k+1)
     &     .and.Usable(i,j-1,k))
      ddydzmm(i,j,k,n)=(crse(i,j,k,n)-crse(i,j-1,k,n)
     &     -crse(i,j,k-1,n)+crse(i,j-1,k-1,n))/(hc(2)*hc(3))
      okyzmm(i,j,k)=(Usable(i,j-1,k)
     &     .and.Usable(i,j,k-1).and.Usable(i,j-1,k-1))

      ddxlowp(i,j,k,n)=(crse(i+1,j,k,n)-crse(i,j,k,n))/hc(1)
      ddylowp(i,j,k,n)=(crse(i,j+1,k,n)-crse(i,j,k,n))/hc(2)
      ddzlowp(i,j,k,n)=(crse(i,j,k+1,n)-crse(i,j,k,n))/hc(3)

      ddxlowm(i,j,k,n)=(-crse(i-1,j,k,n)+crse(i,j,k,n))/hc(1)
      ddylowm(i,j,k,n)=(-crse(i,j-1,k,n)+crse(i,j,k,n))/hc(2)
      ddzlowm(i,j,k,n)=(-crse(i,j,k-1,n)+crse(i,j,k,n))/hc(3)
c ::: ------------ functions ------------------------------

      kclo = ARG_L3(cb)
      kchi = ARG_H3(cb)
      iclo = ARG_L1(cb)
      ichi = ARG_H1(cb)
      jc   = ARG_L2(cb)-1
      j    = lo(2)-1
c ::: :::
      hf(1) = hfine(1)
      hf(2) = 0
      hf(3) = hfine(2)
      hc(1) = hfine(1)*ratio
      hc(2) = 0
      hc(3) = hfine(2)*ratio
c ::: :::

      do n = 1, nvar
c     ::::: calculate derivatives for interpolation
         do kc = kclo, kchi
            k = ratio*kc
            do ic = iclo, ichi
               i = ratio*ic

               derives(ic,kc,DX) = ddxcen(ic,jc,kc,n)
               derives(ic,kc,D2X) = ddx2cen(ic,jc,kc,n)

               if( notUsable(i-1,j,k) .and. notUsable(i+ratio,j,k)  ) then
                  derives(ic,kc,DX)  = zero
                  derives(ic,kc,D2X) = zero
               else if( notUsable(i-1,j,k) .and. Usable(i+ratio,j,k) ) then
c     :::: :::: : positive side usable
                  if( ic+2.le.ichi+1 ) then
                     if( Usable(i+2*ratio,j,k) ) then
                        derives(ic,kc,DX)=ddxplus(ic,jc,kc,n)
                        derives(ic,kc,D2X)=ddx2cen(ic+1,jc,kc,n)
                     else
                        derives(ic,kc,DX)=ddxlowp(ic,jc,kc,n)
                        derives(ic,kc,D2X)=zero
                     endif
                  else
                     derives(ic,kc,DX)=ddxlowp(ic,jc,kc,n)
                     derives(ic,kc,D2X)=zero
                  endif
               else if( Usable(i-1,j,k) .and. notUsable(i+ratio,j,k) ) then
c     :::: :::: : negative side usable
                  if( ic-2.ge.iclo-1 ) then
                     if( Usable(i-1,j,k) ) then
                        derives(ic,kc,DX)=ddxminus(ic,jc,kc,n)
                        derives(ic,kc,D2X)=ddx2cen(ic-1,jc,kc,n)
                     else
                        derives(ic,kc,DX)=ddxlowm(ic,jc,kc,n)
                        derives(ic,kc,D2X)=zero
                     endif
                  else
                     derives(ic,kc,DX)=ddxlowm(ic,jc,kc,n)
                     derives(ic,kc,D2X)=zero
                  endif
               endif

               derives(ic,kc,DZ)  = ddzcen(ic,jc,kc,n)
               derives(ic,kc,D2Z) = ddz2cen(ic,jc,kc,n)

               if( notUsable(i,j,k-1) .and. notUsable(i,j,k+ratio)) then
                  derives(ic,kc,DZ)  = zero
                  derives(ic,kc,D2Z) = zero
               else if( notUsable(i,j,k-1) .and. Usable(i,j,k+ratio)) then
c     :::: :::: : positive size usable
                  if( kc+2.le.kchi+1 ) then
                     if( Usable(i,j,k+2*ratio) ) then
                        derives(ic,kc,DZ)  = ddzplus(ic,jc,kc,n)
                        derives(ic,kc,D2Z) =  ddz2cen(ic,jc,kc+1,n)
                     else
                        derives(ic,kc,DZ)  = ddzlowp(ic,jc,kc,n)
                        derives(ic,kc,D2Z) =  zero
                     endif
                  else
                     derives(ic,kc,DZ)  = ddzlowp(ic,jc,kc,n)
                     derives(ic,kc,D2Z) =  zero
                  endif
               else if( Usable(i,j,k-1) .and. notUsable(i,j,k+ratio)) then
c     :::: :::: : negative size usable
                  if( kc-2.ge.kclo-1 ) then
                     if( Usable(i,j,k-1) ) then
                        derives(ic,kc,DZ)  = ddzminus(ic,jc,kc,n)
                        derives(ic,kc,D2Z) =  ddz2cen(ic,jc,kc-1,n)
                     else
                        derives(ic,kc,DZ)  = ddzlowm(ic,jc,kc,n)
                        derives(ic,kc,D2Z) =  zero
                     endif
                  else
                     derives(ic,kc,DZ)  = ddzlowm(ic,jc,kc,n)
                     derives(ic,kc,D2Z) =  zero
                  endif
               endif

               if( okxzcen(i,j,k) ) then
                  derives(ic,kc,DXZ) = ddxdzcen(ic,jc,kc,n)
               else if( okxzpp(i,j,k) ) then
                  derives(ic,kc,DXZ) = ddxdzpp(ic,jc,kc,n)
               else if( okxzpm(i,j,k) ) then
                  derives(ic,kc,DXZ) = ddxdzpm(ic,jc,kc,n)
               else if( okxzmp(i,j,k) ) then
                  derives(ic,kc,DXZ) = ddxdzmp(ic,jc,kc,n)
               else if( okxzmm(i,j,k) ) then
                  derives(ic,kc,DXZ) = ddxdzmm(ic,jc,kc,n)
               else
                  derives(ic,kc,DXZ) = zero
               endif

            enddo
         enddo

c     ::::: interpolate to fine grid
         do koff = 0, ratio - 1
            zz = (koff+half)*hf(3)-hc(3)/2
            do kc = kclo,kchi
               k = ratio*kc + koff
               do ioff = 0, ratio - 1
                  xx = (ioff+half)*hf(1)-hc(1)/2
                  do ic = iclo, ichi
                     i = ratio*ic + ioff
                     bdry(i,j,k,n) = crse(ic,jc,kc,n) + 
     &                    xx*derives(ic,kc,DX) +
     &                    half*xx**2*derives(ic,kc,D2X) + 
     &                    zz*derives(ic,kc,DZ) +
     &                    half*zz**2*derives(ic,kc,D2Z) + 
     &                    xx*zz*derives(ic,kc,DXZ) 
                     bdry(i,j,k,n+CX*nvar) =
     &                    derives(ic,kc,DX)+
     &                    xx*derives(ic,kc,D2X)+
     &                    zz*derives(ic,kc,DXZ)
                     bdry(i,j,k,n+CZ*nvar) =
     &                    derives(ic,kc,DZ)+
     &                    zz*derives(ic,kc,D2Z)+
     &                    xx*derives(ic,kc,DXZ)
                  enddo
               enddo
            enddo
         enddo

      enddo
            
      return
      end


c ---------------------------------------------------------------
c ::  FORT_BDINTERPYHI : Interpolation on Yhi Face
c ::       Quadratic Interpolation from crse data
c ::       in directions transverse to face of grid
c ::
c ::  Inputs/Outputs:
c ::  bdry       <=  fine grid bndry data strip
c ::  DIMS(bdry)  => index limits of bdry
c ::  lo,hi       => index limits of grd interior
c ::  DIMS(cb)    => index limits of coarsened grid interior
c ::  nvar        => number of variables in solution, half number in bdry
c ::  ratio       => refinement ratio
c ::  not_covered => mask is set to this value if cell is not
c ::                 covered by another fine grid and not outside the domain.
c ::  mask        => fine grid mask bndry strip
c ::  DIMS(mask)  => index limits of mask array
c ::  crse        => crse grid bndry data strip
c ::  DIMS(crse)  => index limits of crse array
c ::  derives     => crse grid tmp array for derivatives
c ---------------------------------------------------------------

      subroutine FORT_BDINTERPYHI (bdry,DIMS(bdry),
     &           lo,hi,DIMS(cb),nvar,ratio,not_covered,
     &           mask,DIMS(mask),crse,DIMS(crse),derives,hfine)

      integer  nvar, ratio, not_covered
      integer  lo(SDIM), hi(SDIM)
      integer  DIMDEC(bdry)
      integer  DIMDEC(cb)
      integer  DIMDEC(mask)
      integer  DIMDEC(crse)
      REAL_T   bdry(DIMV(bdry),(1+3)*nvar)
      REAL_T   derives(DIM13(cb),NUMDERIV)
      integer  mask(DIMV(mask))
      REAL_T   crse(DIMV(crse),nvar)
      REAL_T   hfine(BL_SPACEDIM-1)

      REAL_T   xx, zz
      integer  i, j, k, ic, jc, kc, ioff, koff, n
      integer  iclo, ichi, kclo, kchi
      REAL_T   hc(3),hf(3)

c ::: ------------ functions ------------------------------
      logical notUsable,Usable
      REAL_T ddxcen,ddx2cen,ddxplus,ddxminus
      REAL_T ddycen,ddy2cen,ddyplus,ddyminus
      REAL_T ddzcen,ddz2cen,ddzplus,ddzminus
      REAL_T ddxdycen,ddxdypp,ddxdypm,ddxdymp,ddxdymm
      REAL_T ddxdzcen,ddxdzpp,ddxdzpm,ddxdzmp,ddxdzmm
      REAL_T ddydzcen,ddydzpp,ddydzpm,ddydzmp,ddydzmm
      logical okxypp,okxypm,okxymp,okxymm,okxycen
      logical okxzpp,okxzpm,okxzmp,okxzmm,okxzcen
      logical okyzpp,okyzpm,okyzmp,okyzmm,okyzcen
      REAL_T ddxlowp,ddxlowm
      REAL_T ddylowp,ddylowm
      REAL_T ddzlowp,ddzlowm

      notUsable(i,j,k) = (mask(i,j,k).ne.not_covered)
      Usable(i,j,k) = (mask(i,j,k).eq.not_covered)

      ddxcen(i,j,k,n) = (crse(i+1,j,k,n)-crse(i-1,j,k,n))/(2*hc(1))
      ddycen(i,j,k,n) = (crse(i,j+1,k,n)-crse(i,j-1,k,n))/(2*hc(2))
      ddzcen(i,j,k,n) = (crse(i,j,k+1,n)-crse(i,j,k-1,n))/(2*hc(3))

      ddx2cen(i,j,k,n)=(crse(i+1,j,k,n)-2*crse(i,j,k,n)+crse(i-1,j,k,n))
     &     /(hc(1)**2)
      ddy2cen(i,j,k,n)=(crse(i,j+1,k,n)-2*crse(i,j,k,n)+crse(i,j-1,k,n))
     &     /(hc(2)**2)
      ddz2cen(i,j,k,n)=(crse(i,j,k+1,n)-2*crse(i,j,k,n)+crse(i,j,k-1,n))
     &     /(hc(3)**2)

      ddxplus(i,j,k,n)=(-crse(i+2,j,k,n)+4*crse(i+1,j,k,n)-3*crse(i,j,k,n))
     &     /(2*hc(1))
      ddyplus(i,j,k,n)=(-crse(i,j+2,k,n)+4*crse(i,j+1,k,n)-3*crse(i,j,k,n))
     &     /(2*hc(2))
      ddzplus(i,j,k,n)=(-crse(i,j,k+2,n)+4*crse(i,j,k+1,n)-3*crse(i,j,k,n))
     &     /(2*hc(3))

      ddxminus(i,j,k,n)=(+crse(i-2,j,k,n)-4*crse(i-1,j,k,n)+3*crse(i,j,k,n))
     &     /(2*hc(1))
      ddyminus(i,j,k,n)=(+crse(i,j-2,k,n)-4*crse(i,j-1,k,n)+3*crse(i,j,k,n))
     &     /(2*hc(2))
      ddzminus(i,j,k,n)=(+crse(i,j,k-2,n)-4*crse(i,j,k-1,n)+3*crse(i,j,k,n))
     &     /(2*hc(3))

      ddxdycen(i,j,k,n)=(crse(i+1,j+1,k,n)-crse(i-1,j+1,k,n)
     &     -crse(i+1,j-1,k,n)+crse(i-1,j-1,k,n))/(4*hc(1)*hc(2))
      okxycen(i,j,k)=(Usable(i+1,j+1,k).and.Usable(i-1,j+1,k)
     &     .and.Usable(i+1,j-1,k).and.Usable(i-1,j-1,k))
      ddxdypp(i,j,k,n)=(crse(i+1,j+1,k,n)-crse(i,j+1,k,n)
     &     -crse(i+1,j,k,n)+crse(i,j,k,n))/(hc(1)*hc(2))
      okxypp(i,j,k)=(Usable(i+1,j+1,k).and.Usable(i,j+1,k)
     &     .and.Usable(i+1,j,k))
      ddxdypm(i,j,k,n)=(crse(i+1,j,k,n)-crse(i,j,k,n)
     &     -crse(i+1,j-1,k,n)+crse(i,j-1,k,n))/(4*hc(1)*hc(2))
      okxypm(i,j,k)=(Usable(i+1,j,k)
     &     .and.Usable(i+1,j-1,k).and.Usable(i,j-1,k))
      ddxdymp(i,j,k,n)=(crse(i,j+1,k,n)-crse(i-1,j+1,k,n)
     &     -crse(i,j,k,n)+crse(i-1,j,k,n))/(hc(1)*hc(2))
      okxymp(i,j,k)=(Usable(i,j+1,k).and.Usable(i-1,j+1,k)
     &     .and.Usable(i-1,j,k))
      ddxdymm(i,j,k,n)=(crse(i,j,k,n)-crse(i-1,j,k,n)
     &     -crse(i,j-1,k,n)+crse(i-1,j-1,k,n))/(hc(1)*hc(2))
      okxymm(i,j,k)=(Usable(i-1,j,k)
     &     .and.Usable(i,j-1,k).and.Usable(i-1,j-1,k))

      ddxdzcen(i,j,k,n)=(crse(i+1,j,k+1,n)-crse(i-1,j,k+1,n)
     &     -crse(i+1,j,k-1,n)+crse(i-1,j,k-1,n))/(4*hc(1)*hc(3))
      okxzcen(i,j,k)=(Usable(i+1,j,k+1).and.Usable(i-1,j,k+1)
     &     .and.Usable(i+1,j,k-1).and.Usable(i-1,j,k-1))
      ddxdzpp(i,j,k,n)=(crse(i+1,j,k+1,n)-crse(i,j,k+1,n)
     &     -crse(i+1,j,k,n)+crse(i,j,k,n))/(hc(1)*hc(3))
      okxzpp(i,j,k)=(Usable(i+1,j,k+1).and.Usable(i,j,k+1)
     &     .and.Usable(i+1,j,k))
      ddxdzpm(i,j,k,n)=(crse(i+1,j,k,n)-crse(i,j,k,n)
     &     -crse(i+1,j,k-1,n)+crse(i,j,k-1,n))/(4*hc(1)*hc(3))
      okxzpm(i,j,k)=(Usable(i+1,j,k)
     &     .and.Usable(i+1,j,k-1).and.Usable(i,j,k-1))
      ddxdzmp(i,j,k,n)=(crse(i,j,k+1,n)-crse(i-1,j,k+1,n)
     &     -crse(i,j,k,n)+crse(i-1,j,k,n))/(hc(1)*hc(3))
      okxzmp(i,j,k)=(Usable(i,j,k+1).and.Usable(i-1,j,k+1)
     &     .and.Usable(i-1,j,k))
      ddxdzmm(i,j,k,n)=(crse(i,j,k,n)-crse(i-1,j,k,n)
     &     -crse(i,j,k-1,n)+crse(i-1,j,k-1,n))/(hc(1)*hc(3))
      okxzmm(i,j,k)=(Usable(i-1,j,k)
     &     .and.Usable(i,j,k-1).and.Usable(i-1,j,k-1))

      ddydzcen(i,j,k,n)=(crse(i,j+1,k+1,n)-crse(i,j-1,k+1,n)
     &     -crse(i,j+1,k-1,n)+crse(i,j-1,k-1,n))/(4*hc(2)*hc(3))
      okyzcen(i,j,k)=(Usable(i,j+1,k+1).and.Usable(i,j-1,k+1)
     &     .and.Usable(i,j+1,k-1).and.Usable(i,j-1,k-1))
      ddydzpp(i,j,k,n)=(crse(i,j+1,k+1,n)-crse(i,j,k+1,n)
     &     -crse(i,j+1,k,n)+crse(i,j,k,n))/(hc(2)*hc(3))
      okyzpp(i,j,k)=(Usable(i,j+1,k+1).and.Usable(i,j,k+1)
     &     .and.Usable(i,j+1,k))
      ddydzpm(i,j,k,n)=(crse(i,j+1,k,n)-crse(i,j,k,n)
     &     -crse(i,j+1,k-1,n)+crse(i,j,k-1,n))/(4*hc(2)*hc(3))
      okyzpm(i,j,k)=(Usable(i,j+1,k)
     &     .and.Usable(i,j+1,k-1).and.Usable(i,j,k-1))
      ddydzmp(i,j,k,n)=(crse(i,j,k+1,n)-crse(i,j-1,k+1,n)
     &     -crse(i,j,k,n)+crse(i,j-1,k,n))/(hc(2)*hc(3))
      okyzmp(i,j,k)=(Usable(i,j,k+1).and.Usable(i,j-1,k+1)
     &     .and.Usable(i,j-1,k))
      ddydzmm(i,j,k,n)=(crse(i,j,k,n)-crse(i,j-1,k,n)
     &     -crse(i,j,k-1,n)+crse(i,j-1,k-1,n))/(hc(2)*hc(3))
      okyzmm(i,j,k)=(Usable(i,j-1,k)
     &     .and.Usable(i,j,k-1).and.Usable(i,j-1,k-1))

      ddxlowp(i,j,k,n)=(crse(i+1,j,k,n)-crse(i,j,k,n))/hc(1)
      ddylowp(i,j,k,n)=(crse(i,j+1,k,n)-crse(i,j,k,n))/hc(2)
      ddzlowp(i,j,k,n)=(crse(i,j,k+1,n)-crse(i,j,k,n))/hc(3)

      ddxlowm(i,j,k,n)=(-crse(i-1,j,k,n)+crse(i,j,k,n))/hc(1)
      ddylowm(i,j,k,n)=(-crse(i,j-1,k,n)+crse(i,j,k,n))/hc(2)
      ddzlowm(i,j,k,n)=(-crse(i,j,k-1,n)+crse(i,j,k,n))/hc(3)
c ::: ------------ functions ------------------------------

      kclo = ARG_L3(cb)
      kchi = ARG_H3(cb)
      iclo = ARG_L1(cb)
      ichi = ARG_H1(cb)
      jc   = ARG_H2(cb)+1
      j    = hi(2)+1
c ::: :::
      hf(1) = hfine(1)
      hf(2) = 0
      hf(3) = hfine(2)
      hc(1) = hfine(1)*ratio
      hc(2) = 0
      hc(3) = hfine(2)*ratio
c ::: :::

      do n = 1, nvar
c     ::::: calculate derivatives for interpolation
         do kc = kclo, kchi
            k = ratio*kc
            do ic = iclo, ichi
               i = ratio*ic

               derives(ic,kc,DX) = ddxcen(ic,jc,kc,n)
               derives(ic,kc,D2X) = ddx2cen(ic,jc,kc,n)

               if( notUsable(i-1,j,k) .and. notUsable(i+ratio,j,k)  ) then
                  derives(ic,kc,DX)  = zero
                  derives(ic,kc,D2X) = zero
               else if( notUsable(i-1,j,k) .and. Usable(i+ratio,j,k) ) then
c     :::: :::: : positive side usable
                  if( ic+2.le.ichi+1 ) then
                     if( Usable(i+2*ratio,j,k) ) then
                        derives(ic,kc,DX)=ddxplus(ic,jc,kc,n)
                        derives(ic,kc,D2X)=ddx2cen(ic+1,jc,kc,n)
                     else
                        derives(ic,kc,DX)=ddxlowp(ic,jc,kc,n)
                        derives(ic,kc,D2X)=zero
                     endif
                  else
                     derives(ic,kc,DX)=ddxlowp(ic,jc,kc,n)
                     derives(ic,kc,D2X)=zero
                  endif
               else if( Usable(i-1,j,k) .and. notUsable(i+ratio,j,k) ) then
c     :::: :::: : negative side usable
                  if( ic-2.ge.iclo-1 ) then
                     if( Usable(i-1,j,k) ) then
                        derives(ic,kc,DX)=ddxminus(ic,jc,kc,n)
                        derives(ic,kc,D2X)=ddx2cen(ic-1,jc,kc,n)
                     else
                        derives(ic,kc,DX)=ddxlowm(ic,jc,kc,n)
                        derives(ic,kc,D2X)=zero
                     endif
                  else
                     derives(ic,kc,DX)=ddxlowm(ic,jc,kc,n)
                     derives(ic,kc,D2X)=zero
                  endif
               endif

               derives(ic,kc,DZ)  = ddzcen(ic,jc,kc,n)
               derives(ic,kc,D2Z) = ddz2cen(ic,jc,kc,n)

               if( notUsable(i,j,k-1) .and. notUsable(i,j,k+ratio)) then
                  derives(ic,kc,DZ)  = zero
                  derives(ic,kc,D2Z) = zero
               else if( notUsable(i,j,k-1) .and. Usable(i,j,k+ratio)) then
c     :::: :::: : positive size usable
                  if( kc+2.le.kchi+1 ) then
                     if( Usable(i,j,k+2*ratio) ) then
                        derives(ic,kc,DZ)  = ddzplus(ic,jc,kc,n)
                        derives(ic,kc,D2Z) =  ddz2cen(ic,jc,kc+1,n)
                     else
                        derives(ic,kc,DZ)  = ddzlowp(ic,jc,kc,n)
                        derives(ic,kc,D2Z) =  zero
                     endif
                  else
                     derives(ic,kc,DZ)  = ddzlowp(ic,jc,kc,n)
                     derives(ic,kc,D2Z) =  zero
                  endif
               else if( Usable(i,j,k-1) .and. notUsable(i,j,k+ratio)) then
c     :::: :::: : negative size usable
                  if( kc-2.ge.kclo-1 ) then
                     if( Usable(i,j,k-1) ) then
                        derives(ic,kc,DZ)  = ddzminus(ic,jc,kc,n)
                        derives(ic,kc,D2Z) =  ddz2cen(ic,jc,kc-1,n)
                     else
                        derives(ic,kc,DZ)  = ddzlowm(ic,jc,kc,n)
                        derives(ic,kc,D2Z) =  zero
                     endif
                  else
                     derives(ic,kc,DZ)  = ddzlowm(ic,jc,kc,n)
                     derives(ic,kc,D2Z) =  zero
                  endif
               endif

               if( okxzcen(i,j,k) ) then
                  derives(ic,kc,DXZ) = ddxdzcen(ic,jc,kc,n)
               else if( okxzpp(i,j,k) ) then
                  derives(ic,kc,DXZ) = ddxdzpp(ic,jc,kc,n)
               else if( okxzpm(i,j,k) ) then
                  derives(ic,kc,DXZ) = ddxdzpm(ic,jc,kc,n)
               else if( okxzmp(i,j,k) ) then
                  derives(ic,kc,DXZ) = ddxdzmp(ic,jc,kc,n)
               else if( okxzmm(i,j,k) ) then
                  derives(ic,kc,DXZ) = ddxdzmm(ic,jc,kc,n)
               else
                  derives(ic,kc,DXZ) = zero
               endif

            enddo
         enddo

c     ::::: interpolate to fine grid
         do koff = 0, ratio - 1
            zz = (koff+half)*hf(3)-hc(3)/2
            do kc = kclo,kchi
               k = ratio*kc + koff
               do ioff = 0, ratio - 1
                  xx = (ioff+half)*hf(1)-hc(1)/2
                  do ic = iclo, ichi
                     i = ratio*ic + ioff
                     bdry(i,j,k,n) = crse(ic,jc,kc,n) + 
     &                    xx*derives(ic,kc,DX) +
     &                    half*xx**2*derives(ic,kc,D2X) + 
     &                    zz*derives(ic,kc,DZ) +
     &                    half*zz**2*derives(ic,kc,D2Z) + 
     &                    xx*zz*derives(ic,kc,DXZ) 
                     bdry(i,j,k,n+CX*nvar) =
     &                    derives(ic,kc,DX)+
     &                    xx*derives(ic,kc,D2X)+
     &                    zz*derives(ic,kc,DXZ)
                     bdry(i,j,k,n+CZ*nvar) =
     &                    derives(ic,kc,DZ)+
     &                    zz*derives(ic,kc,D2Z)+
     &                    xx*derives(ic,kc,DXZ)
                  enddo
               enddo
            enddo
         enddo

      enddo
            
      return
      end






c ---------------------------------------------------------------
c ::  FORT_BDINTERPZLO : Interpolation on Zlo Face
c ::       Quadratic Interpolation from crse data
c ::       in directions transverse to face of grid
c ::
c ::  Inputs/Outputs:
c ::  bdry       <=  fine grid bndry data strip
c ::  DIMS(bdry)  => index limits of bdry
c ::  lo,hi       => index limits of grd interior
c ::  DIMS(cb)    => index limits of coarsened grid interior
c ::  nvar        => number of variables in solution, half number in bdry
c ::  ratio       => refinement ratio
c ::  not_covered => mask is set to this value if cell is not
c ::                 covered by another fine grid and not outside the domain.
c ::  mask        => fine grid mask bndry strip
c ::  DIMS(mask)  => index limits of mask array
c ::  crse        => crse grid bndry data strip
c ::  DIMS(crse)  => index limits of crse array
c ::  derives     => crse grid tmp array for derivatives
c ---------------------------------------------------------------

      subroutine FORT_BDINTERPZLO (bdry,DIMS(bdry),
     &           lo,hi,DIMS(cb),nvar,ratio,not_covered,
     &           mask,DIMS(mask),crse,DIMS(crse),derives,hfine)

      integer  nvar, ratio, not_covered
      integer  lo(SDIM), hi(SDIM)
      integer  DIMDEC(bdry)
      integer  DIMDEC(cb)
      integer  DIMDEC(mask)
      integer  DIMDEC(crse)
      REAL_T   bdry(DIMV(bdry),(1+3)*nvar)
      REAL_T   derives(DIM12(cb),NUMDERIV)
      integer  mask(DIMV(mask))
      REAL_T   crse(DIMV(crse),nvar)
      REAL_T   hfine(BL_SPACEDIM-1)

      REAL_T   yy, xx 
      integer  i, j, k, ic, jc, kc, ioff, joff, n
      integer  iclo, ichi, jclo, jchi
      REAL_T   hc(3),hf(3)

c ::: ------------ functions ------------------------------
      logical notUsable,Usable
      REAL_T ddxcen,ddx2cen,ddxplus,ddxminus
      REAL_T ddycen,ddy2cen,ddyplus,ddyminus
      REAL_T ddzcen,ddz2cen,ddzplus,ddzminus
      REAL_T ddxdycen,ddxdypp,ddxdypm,ddxdymp,ddxdymm
      REAL_T ddxdzcen,ddxdzpp,ddxdzpm,ddxdzmp,ddxdzmm
      REAL_T ddydzcen,ddydzpp,ddydzpm,ddydzmp,ddydzmm
      logical okxypp,okxypm,okxymp,okxymm,okxycen
      logical okxzpp,okxzpm,okxzmp,okxzmm,okxzcen
      logical okyzpp,okyzpm,okyzmp,okyzmm,okyzcen
      REAL_T ddxlowp,ddxlowm
      REAL_T ddylowp,ddylowm
      REAL_T ddzlowp,ddzlowm

      notUsable(i,j,k) = (mask(i,j,k).ne.not_covered)
      Usable(i,j,k) = (mask(i,j,k).eq.not_covered)

      ddxcen(i,j,k,n) = (crse(i+1,j,k,n)-crse(i-1,j,k,n))/(2*hc(1))
      ddycen(i,j,k,n) = (crse(i,j+1,k,n)-crse(i,j-1,k,n))/(2*hc(2))
      ddzcen(i,j,k,n) = (crse(i,j,k+1,n)-crse(i,j,k-1,n))/(2*hc(3))

      ddx2cen(i,j,k,n)=(crse(i+1,j,k,n)-2*crse(i,j,k,n)+crse(i-1,j,k,n))
     &     /(hc(1)**2)
      ddy2cen(i,j,k,n)=(crse(i,j+1,k,n)-2*crse(i,j,k,n)+crse(i,j-1,k,n))
     &     /(hc(2)**2)
      ddz2cen(i,j,k,n)=(crse(i,j,k+1,n)-2*crse(i,j,k,n)+crse(i,j,k-1,n))
     &     /(hc(3)**2)

      ddxplus(i,j,k,n)=(-crse(i+2,j,k,n)+4*crse(i+1,j,k,n)-3*crse(i,j,k,n))
     &     /(2*hc(1))
      ddyplus(i,j,k,n)=(-crse(i,j+2,k,n)+4*crse(i,j+1,k,n)-3*crse(i,j,k,n))
     &     /(2*hc(2))
      ddzplus(i,j,k,n)=(-crse(i,j,k+2,n)+4*crse(i,j,k+1,n)-3*crse(i,j,k,n))
     &     /(2*hc(3))

      ddxminus(i,j,k,n)=(+crse(i-2,j,k,n)-4*crse(i-1,j,k,n)+3*crse(i,j,k,n))
     &     /(2*hc(1))
      ddyminus(i,j,k,n)=(+crse(i,j-2,k,n)-4*crse(i,j-1,k,n)+3*crse(i,j,k,n))
     &     /(2*hc(2))
      ddzminus(i,j,k,n)=(+crse(i,j,k-2,n)-4*crse(i,j,k-1,n)+3*crse(i,j,k,n))
     &     /(2*hc(3))

      ddxdycen(i,j,k,n)=(crse(i+1,j+1,k,n)-crse(i-1,j+1,k,n)
     &     -crse(i+1,j-1,k,n)+crse(i-1,j-1,k,n))/(4*hc(1)*hc(2))
      okxycen(i,j,k)=(Usable(i+1,j+1,k).and.Usable(i-1,j+1,k)
     &     .and.Usable(i+1,j-1,k).and.Usable(i-1,j-1,k))
      ddxdypp(i,j,k,n)=(crse(i+1,j+1,k,n)-crse(i,j+1,k,n)
     &     -crse(i+1,j,k,n)+crse(i,j,k,n))/(hc(1)*hc(2))
      okxypp(i,j,k)=(Usable(i+1,j+1,k).and.Usable(i,j+1,k)
     &     .and.Usable(i+1,j,k))
      ddxdypm(i,j,k,n)=(crse(i+1,j,k,n)-crse(i,j,k,n)
     &     -crse(i+1,j-1,k,n)+crse(i,j-1,k,n))/(4*hc(1)*hc(2))
      okxypm(i,j,k)=(Usable(i+1,j,k)
     &     .and.Usable(i+1,j-1,k).and.Usable(i,j-1,k))
      ddxdymp(i,j,k,n)=(crse(i,j+1,k,n)-crse(i-1,j+1,k,n)
     &     -crse(i,j,k,n)+crse(i-1,j,k,n))/(hc(1)*hc(2))
      okxymp(i,j,k)=(Usable(i,j+1,k).and.Usable(i-1,j+1,k)
     &     .and.Usable(i-1,j,k))
      ddxdymm(i,j,k,n)=(crse(i,j,k,n)-crse(i-1,j,k,n)
     &     -crse(i,j-1,k,n)+crse(i-1,j-1,k,n))/(hc(1)*hc(2))
      okxymm(i,j,k)=(Usable(i-1,j,k)
     &     .and.Usable(i,j-1,k).and.Usable(i-1,j-1,k))

      ddxdzcen(i,j,k,n)=(crse(i+1,j,k+1,n)-crse(i-1,j,k+1,n)
     &     -crse(i+1,j,k-1,n)+crse(i-1,j,k-1,n))/(4*hc(1)*hc(3))
      okxzcen(i,j,k)=(Usable(i+1,j,k+1).and.Usable(i-1,j,k+1)
     &     .and.Usable(i+1,j,k-1).and.Usable(i-1,j,k-1))
      ddxdzpp(i,j,k,n)=(crse(i+1,j,k+1,n)-crse(i,j,k+1,n)
     &     -crse(i+1,j,k,n)+crse(i,j,k,n))/(hc(1)*hc(3))
      okxzpp(i,j,k)=(Usable(i+1,j,k+1).and.Usable(i,j,k+1)
     &     .and.Usable(i+1,j,k))
      ddxdzpm(i,j,k,n)=(crse(i+1,j,k,n)-crse(i,j,k,n)
     &     -crse(i+1,j,k-1,n)+crse(i,j,k-1,n))/(4*hc(1)*hc(3))
      okxzpm(i,j,k)=(Usable(i+1,j,k)
     &     .and.Usable(i+1,j,k-1).and.Usable(i,j,k-1))
      ddxdzmp(i,j,k,n)=(crse(i,j,k+1,n)-crse(i-1,j,k+1,n)
     &     -crse(i,j,k,n)+crse(i-1,j,k,n))/(hc(1)*hc(3))
      okxzmp(i,j,k)=(Usable(i,j,k+1).and.Usable(i-1,j,k+1)
     &     .and.Usable(i-1,j,k))
      ddxdzmm(i,j,k,n)=(crse(i,j,k,n)-crse(i-1,j,k,n)
     &     -crse(i,j,k-1,n)+crse(i-1,j,k-1,n))/(hc(1)*hc(3))
      okxzmm(i,j,k)=(Usable(i-1,j,k)
     &     .and.Usable(i,j,k-1).and.Usable(i-1,j,k-1))

      ddydzcen(i,j,k,n)=(crse(i,j+1,k+1,n)-crse(i,j-1,k+1,n)
     &     -crse(i,j+1,k-1,n)+crse(i,j-1,k-1,n))/(4*hc(2)*hc(3))
      okyzcen(i,j,k)=(Usable(i,j+1,k+1).and.Usable(i,j-1,k+1)
     &     .and.Usable(i,j+1,k-1).and.Usable(i,j-1,k-1))
      ddydzpp(i,j,k,n)=(crse(i,j+1,k+1,n)-crse(i,j,k+1,n)
     &     -crse(i,j+1,k,n)+crse(i,j,k,n))/(hc(2)*hc(3))
      okyzpp(i,j,k)=(Usable(i,j+1,k+1).and.Usable(i,j,k+1)
     &     .and.Usable(i,j+1,k))
      ddydzpm(i,j,k,n)=(crse(i,j+1,k,n)-crse(i,j,k,n)
     &     -crse(i,j+1,k-1,n)+crse(i,j,k-1,n))/(4*hc(2)*hc(3))
      okyzpm(i,j,k)=(Usable(i,j+1,k)
     &     .and.Usable(i,j+1,k-1).and.Usable(i,j,k-1))
      ddydzmp(i,j,k,n)=(crse(i,j,k+1,n)-crse(i,j-1,k+1,n)
     &     -crse(i,j,k,n)+crse(i,j-1,k,n))/(hc(2)*hc(3))
      okyzmp(i,j,k)=(Usable(i,j,k+1).and.Usable(i,j-1,k+1)
     &     .and.Usable(i,j-1,k))
      ddydzmm(i,j,k,n)=(crse(i,j,k,n)-crse(i,j-1,k,n)
     &     -crse(i,j,k-1,n)+crse(i,j-1,k-1,n))/(hc(2)*hc(3))
      okyzmm(i,j,k)=(Usable(i,j-1,k)
     &     .and.Usable(i,j,k-1).and.Usable(i,j-1,k-1))

      ddxlowp(i,j,k,n)=(crse(i+1,j,k,n)-crse(i,j,k,n))/hc(1)
      ddylowp(i,j,k,n)=(crse(i,j+1,k,n)-crse(i,j,k,n))/hc(2)
      ddzlowp(i,j,k,n)=(crse(i,j,k+1,n)-crse(i,j,k,n))/hc(3)

      ddxlowm(i,j,k,n)=(-crse(i-1,j,k,n)+crse(i,j,k,n))/hc(1)
      ddylowm(i,j,k,n)=(-crse(i,j-1,k,n)+crse(i,j,k,n))/hc(2)
      ddzlowm(i,j,k,n)=(-crse(i,j,k-1,n)+crse(i,j,k,n))/hc(3)
c ::: ------------ functions ------------------------------

      jclo = ARG_L2(cb)
      jchi = ARG_H2(cb)
      iclo = ARG_L1(cb)
      ichi = ARG_H1(cb)
      kc   = ARG_L3(cb)-1
      k    = lo(3)-1
c ::: :::
      hf(1) = hfine(1)
      hf(2) = hfine(2)
      hf(3) = 0
      hc(1) = hfine(1)*ratio
      hc(2) = hfine(2)*ratio
      hc(3) = 0
c ::: :::

      do n = 1, nvar
c     ::::: calculate derivatives for interpolation
         do jc = jclo, jchi
            j = ratio*jc
            do ic = iclo, ichi
               i = ratio*ic

               derives(ic,jc,DX) = ddxcen(ic,jc,kc,n)
               derives(ic,jc,D2X) = ddx2cen(ic,jc,kc,n)

               if( notUsable(i-1,j,k) .and. notUsable(i+ratio,j,k)  ) then
                  derives(ic,jc,DX)  = zero
                  derives(ic,jc,D2X) = zero
               else if( notUsable(i-1,j,k) .and. Usable(i+ratio,j,k) ) then
c     :::: :::: : positive side usable
                  if( ic+2.le.ichi+1 ) then
                     if( Usable(i+2*ratio,j,k) ) then
                        derives(ic,jc,DX)=ddxplus(ic,jc,kc,n)
                        derives(ic,jc,D2X)=ddx2cen(ic+1,jc,kc,n)
                     else
                        derives(ic,jc,DX)=ddxlowp(ic,jc,kc,n)
                        derives(ic,jc,D2X)=zero
                     endif
                  else
                     derives(ic,jc,DX)=ddxlowp(ic,jc,kc,n)
                     derives(ic,jc,D2X)=zero
                  endif
               else if( Usable(i-1,j,k) .and. notUsable(i+ratio,j,k) ) then
c     :::: :::: : negative side usable
                  if( ic-2.ge.iclo-1 ) then
                     if( Usable(i-1,j,k) ) then
                        derives(ic,jc,DX)=ddxminus(ic,jc,kc,n)
                        derives(ic,jc,D2X)=ddx2cen(ic-1,jc,kc,n)
                     else
                        derives(ic,jc,DX)=ddxlowm(ic,jc,kc,n)
                        derives(ic,jc,D2X)=zero
                     endif
                  else
                     derives(ic,jc,DX)=ddxlowm(ic,jc,kc,n)
                     derives(ic,jc,D2X)=zero
                  endif
               endif

               derives(ic,jc,DY) = ddycen(ic,jc,kc,n)
               derives(ic,jc,D2Y) = ddy2cen(ic,jc,kc,n)

               if( notUsable(i,j-1,k) .and. notUsable(i,j+ratio,k)  ) then
                  derives(ic,jc,DY)  = zero
                  derives(ic,jc,D2Y) = zero
               else if( notUsable(i,j-1,k) .and. Usable(i,j+ratio,k) ) then
c     :::: :::: : positive side usable
                  if( jc+2.le.jchi+1 ) then
                     if( Usable(i,j+2*ratio,k) ) then
                        derives(ic,jc,DY)=ddyplus(ic,jc,kc,n)
                        derives(ic,jc,D2Y)=ddy2cen(ic,jc+1,kc,n)
                     else
                        derives(ic,jc,DY)=ddylowp(ic,jc,kc,n)
                        derives(ic,jc,D2Y)=zero
                     endif
                  else
                     derives(ic,jc,DY)=ddylowp(ic,jc,kc,n)
                     derives(ic,jc,D2Y)=zero
                  endif
               else if( Usable(i,j-1,k) .and. notUsable(i,j+ratio,k) ) then
c     :::: :::: : negative side usable
                  if( jc-2.ge.jclo-1 ) then
                     if( Usable(i,j-1,k) ) then
                        derives(ic,jc,DY)=ddyminus(ic,jc,kc,n)
                        derives(ic,jc,D2Y)=ddy2cen(ic,jc-1,kc,n)
                     else
                        derives(ic,jc,DY)=ddylowm(ic,jc,kc,n)
                        derives(ic,jc,D2Y)=zero
                     endif
                  else
                     derives(ic,jc,DY)=ddylowm(ic,jc,kc,n)
                     derives(ic,jc,D2Y)=zero
                  endif
               endif


               if( okxycen(i,j,k) ) then
                  derives(ic,jc,DXY) = ddxdycen(ic,jc,kc,n)
               else if( okxypp(i,j,k) ) then
                  derives(ic,jc,DXY) = ddxdypp(ic,jc,kc,n)
               else if( okxypm(i,j,k) ) then
                  derives(ic,jc,DXY) = ddxdypm(ic,jc,kc,n)
               else if( okxymp(i,j,k) ) then
                  derives(ic,jc,DXY) = ddxdymp(ic,jc,kc,n)
               else if( okxymm(i,j,k) ) then
                  derives(ic,jc,DXY) = ddxdymm(ic,jc,kc,n)
               else
                  derives(ic,jc,DXY) = zero
               endif

            enddo
         enddo

c     ::::: interpolate to fine grid
         do joff = 0, ratio - 1
            yy = (joff+half)*hf(2)-hc(2)/2
            do jc = jclo,jchi
               j = ratio*jc + joff
               do ioff = 0, ratio - 1
                  xx = (ioff+half)*hf(1)-hc(1)/2
                  do ic = iclo, ichi
                     i = ratio*ic + ioff
                     bdry(i,j,k,n) = crse(ic,jc,kc,n) + 
     &                    xx*derives(ic,jc,DX) +
     &                    half*xx**2*derives(ic,jc,D2X) + 
     &                    yy*derives(ic,jc,DY) +
     &                    half*yy**2*derives(ic,jc,D2Y) + 
     &                    xx*yy*derives(ic,jc,DXY) 
                     bdry(i,j,k,n+CX*nvar) =
     &                    derives(ic,jc,DX)+
     &                    xx*derives(ic,jc,D2X)+
     &                    yy*derives(ic,jc,DXY)
                     bdry(i,j,k,n+CY*nvar) =
     &                    derives(ic,jc,DY)+
     &                    yy*derives(ic,jc,D2Y)+
     &                    xx*derives(ic,jc,DXY)
                  enddo
               enddo
            enddo
         enddo

      enddo
            
      return
      end


c ---------------------------------------------------------------
c ::  FORT_BDINTERPZHI : Interpolation on Zhi Face
c ::       Quadratic Interpolation from crse data
c ::       in directions transverse to face of grid
c ::
c ::  Inputs/Outputs:
c ::  bdry       <=  fine grid bndry data strip
c ::  DIMS(bdry)  => index limits of bdry
c ::  lo,hi       => index limits of grd interior
c ::  DIMS(cb)    => index limits of coarsened grid interior
c ::  nvar        => number of variables in solution, half number in bdry
c ::  ratio       => refinement ratio
c ::  not_covered => mask is set to this value if cell is not
c ::                 covered by another fine grid and not outside the domain.
c ::  mask        => fine grid mask bndry strip
c ::  DIMS(mask)  => index limits of mask array
c ::  crse        => crse grid bndry data strip
c ::  DIMS(crse)  => index limits of crse array
c ::  derives     => crse grid tmp array for derivatives
c ---------------------------------------------------------------

      subroutine FORT_BDINTERPZHI (bdry,DIMS(bdry),
     &           lo,hi,DIMS(cb),nvar,ratio,not_covered,
     &           mask,DIMS(mask),crse,DIMS(crse),derives,hfine)

      integer  nvar, ratio, not_covered
      integer  lo(SDIM), hi(SDIM)
      integer  DIMDEC(bdry)
      integer  DIMDEC(cb)
      integer  DIMDEC(mask)
      integer  DIMDEC(crse)
      REAL_T   bdry(DIMV(bdry),(1+3)*nvar)
      REAL_T   derives(DIM12(cb),NUMDERIV)
      integer  mask(DIMV(mask))
      REAL_T   crse(DIMV(crse),nvar)
      REAL_T   hfine(BL_SPACEDIM-1)

      REAL_T   yy, xx 
      integer  i, j, k, ic, jc, kc, ioff, joff, n
      integer  iclo, ichi, jclo, jchi
      REAL_T   hc(3),hf(3)

c ::: ------------ functions ------------------------------
      logical notUsable,Usable
      REAL_T ddxcen,ddx2cen,ddxplus,ddxminus
      REAL_T ddycen,ddy2cen,ddyplus,ddyminus
      REAL_T ddzcen,ddz2cen,ddzplus,ddzminus
      REAL_T ddxdycen,ddxdypp,ddxdypm,ddxdymp,ddxdymm
      REAL_T ddxdzcen,ddxdzpp,ddxdzpm,ddxdzmp,ddxdzmm
      REAL_T ddydzcen,ddydzpp,ddydzpm,ddydzmp,ddydzmm
      logical okxypp,okxypm,okxymp,okxymm,okxycen
      logical okxzpp,okxzpm,okxzmp,okxzmm,okxzcen
      logical okyzpp,okyzpm,okyzmp,okyzmm,okyzcen
      REAL_T ddxlowp,ddxlowm
      REAL_T ddylowp,ddylowm
      REAL_T ddzlowp,ddzlowm

      notUsable(i,j,k) = (mask(i,j,k).ne.not_covered)
      Usable(i,j,k) = (mask(i,j,k).eq.not_covered)

      ddxcen(i,j,k,n) = (crse(i+1,j,k,n)-crse(i-1,j,k,n))/(2*hc(1))
      ddycen(i,j,k,n) = (crse(i,j+1,k,n)-crse(i,j-1,k,n))/(2*hc(2))
      ddzcen(i,j,k,n) = (crse(i,j,k+1,n)-crse(i,j,k-1,n))/(2*hc(3))

      ddx2cen(i,j,k,n)=(crse(i+1,j,k,n)-2*crse(i,j,k,n)+crse(i-1,j,k,n))
     &     /(hc(1)**2)
      ddy2cen(i,j,k,n)=(crse(i,j+1,k,n)-2*crse(i,j,k,n)+crse(i,j-1,k,n))
     &     /(hc(2)**2)
      ddz2cen(i,j,k,n)=(crse(i,j,k+1,n)-2*crse(i,j,k,n)+crse(i,j,k-1,n))
     &     /(hc(3)**2)

      ddxplus(i,j,k,n)=(-crse(i+2,j,k,n)+4*crse(i+1,j,k,n)-3*crse(i,j,k,n))
     &     /(2*hc(1))
      ddyplus(i,j,k,n)=(-crse(i,j+2,k,n)+4*crse(i,j+1,k,n)-3*crse(i,j,k,n))
     &     /(2*hc(2))
      ddzplus(i,j,k,n)=(-crse(i,j,k+2,n)+4*crse(i,j,k+1,n)-3*crse(i,j,k,n))
     &     /(2*hc(3))

      ddxminus(i,j,k,n)=(+crse(i-2,j,k,n)-4*crse(i-1,j,k,n)+3*crse(i,j,k,n))
     &     /(2*hc(1))
      ddyminus(i,j,k,n)=(+crse(i,j-2,k,n)-4*crse(i,j-1,k,n)+3*crse(i,j,k,n))
     &     /(2*hc(2))
      ddzminus(i,j,k,n)=(+crse(i,j,k-2,n)-4*crse(i,j,k-1,n)+3*crse(i,j,k,n))
     &     /(2*hc(3))

      ddxdycen(i,j,k,n)=(crse(i+1,j+1,k,n)-crse(i-1,j+1,k,n)
     &     -crse(i+1,j-1,k,n)+crse(i-1,j-1,k,n))/(4*hc(1)*hc(2))
      okxycen(i,j,k)=(Usable(i+1,j+1,k).and.Usable(i-1,j+1,k)
     &     .and.Usable(i+1,j-1,k).and.Usable(i-1,j-1,k))
      ddxdypp(i,j,k,n)=(crse(i+1,j+1,k,n)-crse(i,j+1,k,n)
     &     -crse(i+1,j,k,n)+crse(i,j,k,n))/(hc(1)*hc(2))
      okxypp(i,j,k)=(Usable(i+1,j+1,k).and.Usable(i,j+1,k)
     &     .and.Usable(i+1,j,k))
      ddxdypm(i,j,k,n)=(crse(i+1,j,k,n)-crse(i,j,k,n)
     &     -crse(i+1,j-1,k,n)+crse(i,j-1,k,n))/(4*hc(1)*hc(2))
      okxypm(i,j,k)=(Usable(i+1,j,k)
     &     .and.Usable(i+1,j-1,k).and.Usable(i,j-1,k))
      ddxdymp(i,j,k,n)=(crse(i,j+1,k,n)-crse(i-1,j+1,k,n)
     &     -crse(i,j,k,n)+crse(i-1,j,k,n))/(hc(1)*hc(2))
      okxymp(i,j,k)=(Usable(i,j+1,k).and.Usable(i-1,j+1,k)
     &     .and.Usable(i-1,j,k))
      ddxdymm(i,j,k,n)=(crse(i,j,k,n)-crse(i-1,j,k,n)
     &     -crse(i,j-1,k,n)+crse(i-1,j-1,k,n))/(hc(1)*hc(2))
      okxymm(i,j,k)=(Usable(i-1,j,k)
     &     .and.Usable(i,j-1,k).and.Usable(i-1,j-1,k))

      ddxdzcen(i,j,k,n)=(crse(i+1,j,k+1,n)-crse(i-1,j,k+1,n)
     &     -crse(i+1,j,k-1,n)+crse(i-1,j,k-1,n))/(4*hc(1)*hc(3))
      okxzcen(i,j,k)=(Usable(i+1,j,k+1).and.Usable(i-1,j,k+1)
     &     .and.Usable(i+1,j,k-1).and.Usable(i-1,j,k-1))
      ddxdzpp(i,j,k,n)=(crse(i+1,j,k+1,n)-crse(i,j,k+1,n)
     &     -crse(i+1,j,k,n)+crse(i,j,k,n))/(hc(1)*hc(3))
      okxzpp(i,j,k)=(Usable(i+1,j,k+1).and.Usable(i,j,k+1)
     &     .and.Usable(i+1,j,k))
      ddxdzpm(i,j,k,n)=(crse(i+1,j,k,n)-crse(i,j,k,n)
     &     -crse(i+1,j,k-1,n)+crse(i,j,k-1,n))/(4*hc(1)*hc(3))
      okxzpm(i,j,k)=(Usable(i+1,j,k)
     &     .and.Usable(i+1,j,k-1).and.Usable(i,j,k-1))
      ddxdzmp(i,j,k,n)=(crse(i,j,k+1,n)-crse(i-1,j,k+1,n)
     &     -crse(i,j,k,n)+crse(i-1,j,k,n))/(hc(1)*hc(3))
      okxzmp(i,j,k)=(Usable(i,j,k+1).and.Usable(i-1,j,k+1)
     &     .and.Usable(i-1,j,k))
      ddxdzmm(i,j,k,n)=(crse(i,j,k,n)-crse(i-1,j,k,n)
     &     -crse(i,j,k-1,n)+crse(i-1,j,k-1,n))/(hc(1)*hc(3))
      okxzmm(i,j,k)=(Usable(i-1,j,k)
     &     .and.Usable(i,j,k-1).and.Usable(i-1,j,k-1))

      ddydzcen(i,j,k,n)=(crse(i,j+1,k+1,n)-crse(i,j-1,k+1,n)
     &     -crse(i,j+1,k-1,n)+crse(i,j-1,k-1,n))/(4*hc(2)*hc(3))
      okyzcen(i,j,k)=(Usable(i,j+1,k+1).and.Usable(i,j-1,k+1)
     &     .and.Usable(i,j+1,k-1).and.Usable(i,j-1,k-1))
      ddydzpp(i,j,k,n)=(crse(i,j+1,k+1,n)-crse(i,j,k+1,n)
     &     -crse(i,j+1,k,n)+crse(i,j,k,n))/(hc(2)*hc(3))
      okyzpp(i,j,k)=(Usable(i,j+1,k+1).and.Usable(i,j,k+1)
     &     .and.Usable(i,j+1,k))
      ddydzpm(i,j,k,n)=(crse(i,j+1,k,n)-crse(i,j,k,n)
     &     -crse(i,j+1,k-1,n)+crse(i,j,k-1,n))/(4*hc(2)*hc(3))
      okyzpm(i,j,k)=(Usable(i,j+1,k)
     &     .and.Usable(i,j+1,k-1).and.Usable(i,j,k-1))
      ddydzmp(i,j,k,n)=(crse(i,j,k+1,n)-crse(i,j-1,k+1,n)
     &     -crse(i,j,k,n)+crse(i,j-1,k,n))/(hc(2)*hc(3))
      okyzmp(i,j,k)=(Usable(i,j,k+1).and.Usable(i,j-1,k+1)
     &     .and.Usable(i,j-1,k))
      ddydzmm(i,j,k,n)=(crse(i,j,k,n)-crse(i,j-1,k,n)
     &     -crse(i,j,k-1,n)+crse(i,j-1,k-1,n))/(hc(2)*hc(3))
      okyzmm(i,j,k)=(Usable(i,j-1,k)
     &     .and.Usable(i,j,k-1).and.Usable(i,j-1,k-1))

      ddxlowp(i,j,k,n)=(crse(i+1,j,k,n)-crse(i,j,k,n))/hc(1)
      ddylowp(i,j,k,n)=(crse(i,j+1,k,n)-crse(i,j,k,n))/hc(2)
      ddzlowp(i,j,k,n)=(crse(i,j,k+1,n)-crse(i,j,k,n))/hc(3)

      ddxlowm(i,j,k,n)=(-crse(i-1,j,k,n)+crse(i,j,k,n))/hc(1)
      ddylowm(i,j,k,n)=(-crse(i,j-1,k,n)+crse(i,j,k,n))/hc(2)
      ddzlowm(i,j,k,n)=(-crse(i,j,k-1,n)+crse(i,j,k,n))/hc(3)
c ::: ------------ functions ------------------------------

      jclo = ARG_L2(cb)
      jchi = ARG_H2(cb)
      iclo = ARG_L1(cb)
      ichi = ARG_H1(cb)
      kc   = ARG_H3(cb)+1
      k    = hi(3)+1
c ::: :::
      hf(1) = hfine(1)
      hf(2) = hfine(2)
      hf(3) = 0
      hc(1) = hfine(1)*ratio
      hc(2) = hfine(2)*ratio
      hc(3) = 0
c ::: :::

      do n = 1, nvar
c     ::::: calculate derivatives for interpolation
         do jc = jclo, jchi
            j = ratio*jc
            do ic = iclo, ichi
               i = ratio*ic

               derives(ic,jc,DX) = ddxcen(ic,jc,kc,n)
               derives(ic,jc,D2X) = ddx2cen(ic,jc,kc,n)

               if( notUsable(i-1,j,k) .and. notUsable(i+ratio,j,k)  ) then
                  derives(ic,jc,DX)  = zero
                  derives(ic,jc,D2X) = zero
               else if( notUsable(i-1,j,k) .and. Usable(i+ratio,j,k) ) then
c     :::: :::: : positive side usable
                  if( ic+2.le.ichi+1 ) then
                     if( Usable(i+2*ratio,j,k) ) then
                        derives(ic,jc,DX)=ddxplus(ic,jc,kc,n)
                        derives(ic,jc,D2X)=ddx2cen(ic+1,jc,kc,n)
                     else
                        derives(ic,jc,DX)=ddxlowp(ic,jc,kc,n)
                        derives(ic,jc,D2X)=zero
                     endif
                  else
                     derives(ic,jc,DX)=ddxlowp(ic,jc,kc,n)
                     derives(ic,jc,D2X)=zero
                  endif
               else if( Usable(i-1,j,k) .and. notUsable(i+ratio,j,k) ) then
c     :::: :::: : negative side usable
                  if( ic-2.ge.iclo-1 ) then
                     if( Usable(i-1,j,k) ) then
                        derives(ic,jc,DX)=ddxminus(ic,jc,kc,n)
                        derives(ic,jc,D2X)=ddx2cen(ic-1,jc,kc,n)
                     else
                        derives(ic,jc,DX)=ddxlowm(ic,jc,kc,n)
                        derives(ic,jc,D2X)=zero
                     endif
                  else
                     derives(ic,jc,DX)=ddxlowm(ic,jc,kc,n)
                     derives(ic,jc,D2X)=zero
                  endif
               endif

               derives(ic,jc,DY) = ddycen(ic,jc,kc,n)
               derives(ic,jc,D2Y) = ddy2cen(ic,jc,kc,n)

               if( notUsable(i,j-1,k) .and. notUsable(i,j+ratio,k)  ) then
                  derives(ic,jc,DY)  = zero
                  derives(ic,jc,D2Y) = zero
               else if( notUsable(i,j-1,k) .and. Usable(i,j+ratio,k) ) then
c     :::: :::: : positive side usable
                  if( jc+2.le.jchi+1 ) then
                     if( Usable(i,j+2*ratio,k) ) then
                        derives(ic,jc,DY)=ddyplus(ic,jc,kc,n)
                        derives(ic,jc,D2Y)=ddy2cen(ic,jc+1,kc,n)
                     else
                        derives(ic,jc,DY)=ddylowp(ic,jc,kc,n)
                        derives(ic,jc,D2Y)=zero
                     endif
                  else
                     derives(ic,jc,DY)=ddylowp(ic,jc,kc,n)
                     derives(ic,jc,D2Y)=zero
                  endif
               else if( Usable(i,j-1,k) .and. notUsable(i,j+ratio,k) ) then
c     :::: :::: : negative side usable
                  if( jc-2.ge.jclo-1 ) then
                     if( Usable(i,j-1,k) ) then
                        derives(ic,jc,DY)=ddyminus(ic,jc,kc,n)
                        derives(ic,jc,D2Y)=ddy2cen(ic,jc-1,kc,n)
                     else
                        derives(ic,jc,DY)=ddylowm(ic,jc,kc,n)
                        derives(ic,jc,D2Y)=zero
                     endif
                  else
                     derives(ic,jc,DY)=ddylowm(ic,jc,kc,n)
                     derives(ic,jc,D2Y)=zero
                  endif
               endif


               if( okxycen(i,j,k) ) then
                  derives(ic,jc,DXY) = ddxdycen(ic,jc,kc,n)
               else if( okxypp(i,j,k) ) then
                  derives(ic,jc,DXY) = ddxdypp(ic,jc,kc,n)
               else if( okxypm(i,j,k) ) then
                  derives(ic,jc,DXY) = ddxdypm(ic,jc,kc,n)
               else if( okxymp(i,j,k) ) then
                  derives(ic,jc,DXY) = ddxdymp(ic,jc,kc,n)
               else if( okxymm(i,j,k) ) then
                  derives(ic,jc,DXY) = ddxdymm(ic,jc,kc,n)
               else
                  derives(ic,jc,DXY) = zero
               endif

            enddo
         enddo

c     ::::: interpolate to fine grid
         do joff = 0, ratio - 1
            yy = (joff+half)*hf(2)-hc(2)/2
            do jc = jclo,jchi
               j = ratio*jc + joff
               do ioff = 0, ratio - 1
                  xx = (ioff+half)*hf(1)-hc(1)/2
                  do ic = iclo, ichi
                     i = ratio*ic + ioff
                     bdry(i,j,k,n) = crse(ic,jc,kc,n) + 
     &                    xx*derives(ic,jc,DX) +
     &                    half*xx**2*derives(ic,jc,D2X) + 
     &                    yy*derives(ic,jc,DY) +
     &                    half*yy**2*derives(ic,jc,D2Y) + 
     &                    xx*yy*derives(ic,jc,DXY) 
                     bdry(i,j,k,n+CX*nvar) =
     &                    derives(ic,jc,DX)+
     &                    xx*derives(ic,jc,D2X)+
     &                    yy*derives(ic,jc,DXY)
                     bdry(i,j,k,n+CY*nvar) =
     &                    derives(ic,jc,DY)+
     &                    yy*derives(ic,jc,D2Y)+
     &                    xx*derives(ic,jc,DXY)
                  enddo
               enddo
            enddo
         enddo

      enddo
            
      return
      end

c ------------------------------------------------------------------------
c :: FORT_BDIDERIVXLO : copy from Xlo face of fine fab (getting values
c ::      which were stuck there by the physical boundary condition routines)
c ::      to get the first nvar values in bdry.  The following 3*nvar values
c ::      are gotten by computing the transverse derivative.  Two formulas
c ::      are used: 1) on the interior a centered difference is used; 
c ::      2) at the end points a one-sided difference is used.  No check
c ::      is made to check that the one-sided difference does not poke out
c ::      the other side of the grid.  But if grids are that small, a lot
c ::      of things will be blowing up.
c ::      NOTE: only the transverse derivs are computed, so one of the
c ::      deriv's if filled with error values

c :: Inputs/Outputs:
c ::  bdry <= fine grid bndry data strip
c ::  lo,hi => index limits of "grid". So all action takes place in cells
c ::           outside the grid
c ::  fine  => fine grid data which has phys-bc data stuffed in it
c ::  nvar  => number of components in fine.  But you write (3+1) times as many
c ::           components of data in bdry because you compute transverse
c ::           derivatives as well
c ::  hfine => delta x in transverse direction; needed for derivs
c ------------------------------------------------------------------------
      subroutine FORT_BDIDERIVXLO( bdry, DIMS(bdry),
     &     lo, hi,
     &     fine, DIMS(fine),
     &     nvar, hfine)
      integer nvar
      REAL_T  hfine(3)
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(bdry)
      integer DIMDEC(fine)
      REAL_T  bdry(DIMV(bdry),(3+1)*nvar)
      REAL_T  fine(DIMV(fine),  nvar)

      integer i,j,k,n
      REAL_T h(3)
c ::: ------ begin functions -------
      REAL_T ddxcen
      REAL_T ddycen
      REAL_T ddzcen
      REAL_T ddxplus,ddxminus,ddyplus,ddyminus,ddzplus,ddzminus

      ddxcen(i,j,k,n) = (fine(i+1,j,k,n)-fine(i-1,j,k,n))/(2*h(1))
      ddycen(i,j,k,n) = (fine(i,j+1,k,n)-fine(i,j-1,k,n))/(2*h(2))
      ddzcen(i,j,k,n) = (fine(i,j,k+1,n)-fine(i,j,k-1,n))/(2*h(3))

      ddxplus(i,j,k,n)=(-fine(i+2,j,k,n)+4*fine(i+1,j,k,n)-3*fine(i,j,k,n))
     &     /(2*h(1))
      ddyplus(i,j,k,n)=(-fine(i,j+2,k,n)+4*fine(i,j+1,k,n)-3*fine(i,j,k,n))
     &     /(2*h(2))
      ddzplus(i,j,k,n)=(-fine(i,j,k+2,n)+4*fine(i,j,k+1,n)-3*fine(i,j,k,n))
     &     /(2*h(3))

      ddxminus(i,j,k,n)=(+fine(i-2,j,k,n)-4*fine(i-1,j,k,n)+3*fine(i,j,k,n))
     &     /(2*h(1))
      ddyminus(i,j,k,n)=(+fine(i,j-2,k,n)-4*fine(i,j-1,k,n)+3*fine(i,j,k,n))
     &     /(2*h(2))
      ddzminus(i,j,k,n)=(+fine(i,j,k-2,n)-4*fine(i,j,k-1,n)+3*fine(i,j,k,n))
     &     /(2*h(3))
c ::: ------ end   functions -------

      i = lo(1)-1
      h(1) = 0
      h(2) = hfine(1)
      h(3) = hfine(2)

      do n=1,nvar
c ::: ::: copy boundary values over
         do k=lo(3),hi(3)
            do j=lo(2),hi(2)
               bdry(i,j,k,n) = fine(i,j,k,n)
            enddo
         enddo

c ::: ::: do stinking calculation of derivativ using centered difference
c ::: ::: this ONLY makes sense for DIRICHLET BC.  For other BC,
c ::: ::: we will not use these values.

c ::: ::: First, do y deriv
         do k=lo(3),hi(3)
            j=lo(2)
            bdry(i,j,k,n+CY*nvar) = ddyplus(i,j,k,n)
            do j=lo(2)+1,hi(2)-1
               bdry(i,j,k,n+CY*nvar) = ddycen(i,j,k,n)
            enddo
            j=hi(2)
            bdry(i,j,k,n+CY*nvar) = ddyminus(i,j,k,n)
         enddo
c ::: ::: next, do z deriv
         do j=lo(2),hi(2)
            k=lo(3)
            bdry(i,j,k,n+CZ*nvar) = ddzplus(i,j,k,n)
            do k=lo(3)+1,hi(3)-1
               bdry(i,j,k,n+CZ*nvar) = ddzcen(i,j,k,n)
            enddo
            k=hi(3)
            bdry(i,j,k,n+CZ*nvar) = ddzminus(i,j,k,n)
         enddo

      enddo

      return
      end
      

c ------------------------------------------------------------------------
c :: FORT_BDIDERIVXHI : copy from Xhi face of fine fab (getting values
c ::      which were stuck there by the physical boundary condition routines)
c ::      to get the first nvar values in bdry.  The following 3*nvar values
c ::      are gotten by computing the transverse derivative.  Two formulas
c ::      are used: 1) on the interior a centered difference is used; 
c ::      2) at the end points a one-sided difference is used.  No check
c ::      is made to check that the one-sided difference does not poke out
c ::      the other side of the grid.  But if grids are that small, a lot
c ::      of things will be blowing up.
c ::      NOTE: only the transverse derivs are computed, so one of the
c ::      deriv's if filled with error values

c :: Inputs/Outputs:
c ::  bdry <= fine grid bndry data strip
c ::  lo,hi => index limits of "grid". So all action takes place in cells
c ::           outside the grid
c ::  fine  => fine grid data which has phys-bc data stuffed in it
c ::  nvar  => number of components in fine.  But you write (3+1) times as many
c ::           components of data in bdry because you compute transverse
c ::           derivatives as well
c ::  hfine => delta x in transverse direction; needed for derivs
c ------------------------------------------------------------------------
      subroutine FORT_BDIDERIVXHI( bdry, DIMS(bdry),
     &     lo, hi,
     &     fine, DIMS(fine),
     &     nvar, hfine)
      integer nvar
      REAL_T  hfine(3)
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(bdry)
      integer DIMDEC(fine)
      REAL_T  bdry(DIMV(bdry),(3+1)*nvar)
      REAL_T  fine(DIMV(fine),  nvar)

      integer i,j,k,n
      REAL_T h(3)
c ::: ------ begin functions -------
      REAL_T ddxcen
      REAL_T ddycen
      REAL_T ddzcen
      REAL_T ddxplus,ddxminus,ddyplus,ddyminus,ddzplus,ddzminus

      ddxcen(i,j,k,n) = (fine(i+1,j,k,n)-fine(i-1,j,k,n))/(2*h(1))
      ddycen(i,j,k,n) = (fine(i,j+1,k,n)-fine(i,j-1,k,n))/(2*h(2))
      ddzcen(i,j,k,n) = (fine(i,j,k+1,n)-fine(i,j,k-1,n))/(2*h(3))

      ddxplus(i,j,k,n)=(-fine(i+2,j,k,n)+4*fine(i+1,j,k,n)-3*fine(i,j,k,n))
     &     /(2*h(1))
      ddyplus(i,j,k,n)=(-fine(i,j+2,k,n)+4*fine(i,j+1,k,n)-3*fine(i,j,k,n))
     &     /(2*h(2))
      ddzplus(i,j,k,n)=(-fine(i,j,k+2,n)+4*fine(i,j,k+1,n)-3*fine(i,j,k,n))
     &     /(2*h(3))

      ddxminus(i,j,k,n)=(+fine(i-2,j,k,n)-4*fine(i-1,j,k,n)+3*fine(i,j,k,n))
     &     /(2*h(1))
      ddyminus(i,j,k,n)=(+fine(i,j-2,k,n)-4*fine(i,j-1,k,n)+3*fine(i,j,k,n))
     &     /(2*h(2))
      ddzminus(i,j,k,n)=(+fine(i,j,k-2,n)-4*fine(i,j,k-1,n)+3*fine(i,j,k,n))
     &     /(2*h(3))
c ::: ------ end   functions -------

      i = hi(1)+1
      h(1) = 0
      h(2) = hfine(1)
      h(3) = hfine(2)
      do n=1,nvar
c ::: ::: copy boundary values over
         do k=lo(3),hi(3)
            do j=lo(2),hi(2)
               bdry(i,j,k,n) = fine(i,j,k,n)
            enddo
         enddo

c ::: ::: do stinking calculation of derivativ using centered difference
c ::: ::: this ONLY makes sense for DIRICHLET BC.  For other BC,
c ::: ::: we will not use these values.

c ::: ::: First, do y deriv
         do k=lo(3),hi(3)
            j=lo(2)
            bdry(i,j,k,n+CY*nvar) = ddyplus(i,j,k,n)
            do j=lo(2)+1,hi(2)-1
               bdry(i,j,k,n+CY*nvar) = ddycen(i,j,k,n)
            enddo
            j=hi(2)
            bdry(i,j,k,n+CY*nvar) = ddyminus(i,j,k,n)
         enddo
c ::: ::: next, do z deriv
         do j=lo(2),hi(2)
            k=lo(3)
            bdry(i,j,k,n+CZ*nvar) = ddzplus(i,j,k,n)
            do k=lo(3)+1,hi(3)-1
               bdry(i,j,k,n+CZ*nvar) = ddzcen(i,j,k,n)
            enddo
            k=hi(3)
            bdry(i,j,k,n+CZ*nvar) = ddzminus(i,j,k,n)
         enddo

      enddo

      return
      end
      

c ------------------------------------------------------------------------
c :: FORT_BDIDERIVYLO : copy from Ylo face of fine fab (getting values
c ::      which were stuck there by the physical boundary condition routines)
c ::      to get the first nvar values in bdry.  The following 3*nvar values
c ::      are gotten by computing the transverse derivative.  Two formulas
c ::      are used: 1) on the interior a centered difference is used; 
c ::      2) at the end points a one-sided difference is used.  No check
c ::      is made to check that the one-sided difference does not poke out
c ::      the other side of the grid.  But if grids are that small, a lot
c ::      of things will be blowing up.
c ::      NOTE: only the transverse derivs are computed, so one of the
c ::      deriv's if filled with error values

c :: Inputs/Outputs:
c ::  bdry <= fine grid bndry data strip
c ::  lo,hi => index limits of "grid". So all action takes place in cells
c ::           outside the grid
c ::  fine  => fine grid data which has phys-bc data stuffed in it
c ::  nvar  => number of components in fine.  But you write (3+1) times as many
c ::           components of data in bdry because you compute transverse
c ::           derivatives as well
c ::  hfine => delta x in transverse direction; needed for derivs
c ------------------------------------------------------------------------
      subroutine FORT_BDIDERIVYLO( bdry, DIMS(bdry),
     &     lo, hi,
     &     fine, DIMS(fine),
     &     nvar, hfine)
      integer nvar
      REAL_T  hfine(3)
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(bdry)
      integer DIMDEC(fine)
      REAL_T  bdry(DIMV(bdry),(3+1)*nvar)
      REAL_T  fine(DIMV(fine),  nvar)

      integer i,j,k,n
      REAL_T h(3)
c ::: ------ begin functions -------
      REAL_T ddxcen
      REAL_T ddycen
      REAL_T ddzcen
      REAL_T ddxplus,ddxminus,ddyplus,ddyminus,ddzplus,ddzminus

      ddxcen(i,j,k,n) = (fine(i+1,j,k,n)-fine(i-1,j,k,n))/(2*h(1))
      ddycen(i,j,k,n) = (fine(i,j+1,k,n)-fine(i,j-1,k,n))/(2*h(2))
      ddzcen(i,j,k,n) = (fine(i,j,k+1,n)-fine(i,j,k-1,n))/(2*h(3))

      ddxplus(i,j,k,n)=(-fine(i+2,j,k,n)+4*fine(i+1,j,k,n)-3*fine(i,j,k,n))
     &     /(2*h(1))
      ddyplus(i,j,k,n)=(-fine(i,j+2,k,n)+4*fine(i,j+1,k,n)-3*fine(i,j,k,n))
     &     /(2*h(2))
      ddzplus(i,j,k,n)=(-fine(i,j,k+2,n)+4*fine(i,j,k+1,n)-3*fine(i,j,k,n))
     &     /(2*h(3))

      ddxminus(i,j,k,n)=(+fine(i-2,j,k,n)-4*fine(i-1,j,k,n)+3*fine(i,j,k,n))
     &     /(2*h(1))
      ddyminus(i,j,k,n)=(+fine(i,j-2,k,n)-4*fine(i,j-1,k,n)+3*fine(i,j,k,n))
     &     /(2*h(2))
      ddzminus(i,j,k,n)=(+fine(i,j,k-2,n)-4*fine(i,j,k-1,n)+3*fine(i,j,k,n))
     &     /(2*h(3))
c ::: ------ end   functions -------

      j = lo(2)-1
      h(1) = hfine(1)
      h(2) = 0
      h(3) = hfine(2)
      do n=1,nvar
c ::: ::: copy boundary values over
         do k=lo(3),hi(3)
            do i=lo(1),hi(1)
               bdry(i,j,k,n) = fine(i,j,k,n)
            enddo
         enddo

c ::: ::: do stinking calculation of derivativ using centered difference
c ::: ::: this ONLY makes sense for DIRICHLET BC.  For other BC,
c ::: ::: we will not use these values.

c ::: ::: First, do x deriv
         do k=lo(3),hi(3)
            i=lo(1)
            bdry(i,j,k,n+CX*nvar) = ddxplus(i,j,k,n)
            do i=lo(1)+1,hi(1)-1
               bdry(i,j,k,n+CX*nvar) = ddxcen(i,j,k,n)
            enddo
            i=hi(1)
            bdry(i,j,k,n+CX*nvar) = ddxminus(i,j,k,n)
         enddo
c ::: ::: next, do z deriv
         do i=lo(1),hi(1)
            k=lo(3)
            bdry(i,j,k,n+CZ*nvar) = ddzplus(i,j,k,n)
            do k=lo(3)+1,hi(3)-1
               bdry(i,j,k,n+CZ*nvar) = ddzcen(i,j,k,n)
            enddo
            k=hi(3)
            bdry(i,j,k,n+CZ*nvar) = ddzminus(i,j,k,n)
         enddo

      enddo

      return
      end
      


c ------------------------------------------------------------------------
c :: FORT_BDIDERIVYHI : copy from Yhi face of fine fab (getting values
c ::      which were stuck there by the physical boundary condition routines)
c ::      to get the first nvar values in bdry.  The following 3*nvar values
c ::      are gotten by computing the transverse derivative.  Two formulas
c ::      are used: 1) on the interior a centered difference is used; 
c ::      2) at the end points a one-sided difference is used.  No check
c ::      is made to check that the one-sided difference does not poke out
c ::      the other side of the grid.  But if grids are that small, a lot
c ::      of things will be blowing up.
c ::      NOTE: only the transverse derivs are computed, so one of the
c ::      deriv's if filled with error values

c :: Inputs/Outputs:
c ::  bdry <= fine grid bndry data strip
c ::  lo,hi => index limits of "grid". So all action takes place in cells
c ::           outside the grid
c ::  fine  => fine grid data which has phys-bc data stuffed in it
c ::  nvar  => number of components in fine.  But you write (3+1) times as many
c ::           components of data in bdry because you compute transverse
c ::           derivatives as well
c ::  hfine => delta x in transverse direction; needed for derivs
c ------------------------------------------------------------------------
      subroutine FORT_BDIDERIVYHI( bdry, DIMS(bdry),
     &     lo, hi,
     &     fine, DIMS(fine),
     &     nvar, hfine)
      integer nvar
      REAL_T  hfine(3)
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(bdry)
      integer DIMDEC(fine)
      REAL_T  bdry(DIMV(bdry),(3+1)*nvar)
      REAL_T  fine(DIMV(fine),  nvar)

      integer i,j,k,n
      REAL_T h(3)
c ::: ------ begin functions -------
      REAL_T ddxcen
      REAL_T ddycen
      REAL_T ddzcen
      REAL_T ddxplus,ddxminus,ddyplus,ddyminus,ddzplus,ddzminus

      ddxcen(i,j,k,n) = (fine(i+1,j,k,n)-fine(i-1,j,k,n))/(2*h(1))
      ddycen(i,j,k,n) = (fine(i,j+1,k,n)-fine(i,j-1,k,n))/(2*h(2))
      ddzcen(i,j,k,n) = (fine(i,j,k+1,n)-fine(i,j,k-1,n))/(2*h(3))

      ddxplus(i,j,k,n)=(-fine(i+2,j,k,n)+4*fine(i+1,j,k,n)-3*fine(i,j,k,n))
     &     /(2*h(1))
      ddyplus(i,j,k,n)=(-fine(i,j+2,k,n)+4*fine(i,j+1,k,n)-3*fine(i,j,k,n))
     &     /(2*h(2))
      ddzplus(i,j,k,n)=(-fine(i,j,k+2,n)+4*fine(i,j,k+1,n)-3*fine(i,j,k,n))
     &     /(2*h(3))

      ddxminus(i,j,k,n)=(+fine(i-2,j,k,n)-4*fine(i-1,j,k,n)+3*fine(i,j,k,n))
     &     /(2*h(1))
      ddyminus(i,j,k,n)=(+fine(i,j-2,k,n)-4*fine(i,j-1,k,n)+3*fine(i,j,k,n))
     &     /(2*h(2))
      ddzminus(i,j,k,n)=(+fine(i,j,k-2,n)-4*fine(i,j,k-1,n)+3*fine(i,j,k,n))
     &     /(2*h(3))
c ::: ------ end   functions -------

      j = hi(2)+1
      h(1) = hfine(1)
      h(2) = 0
      h(3) = hfine(2)

      do n=1,nvar
c ::: ::: copy boundary values over
         do k=lo(3),hi(3)
            do i=lo(1),hi(1)
               bdry(i,j,k,n) = fine(i,j,k,n)
            enddo
         enddo

c ::: ::: do stinking calculation of derivativ using centered difference
c ::: ::: this ONLY makes sense for DIRICHLET BC.  For other BC,
c ::: ::: we will not use these values.

c ::: ::: First, do x deriv
         do k=lo(3),hi(3)
            i=lo(1)
            bdry(i,j,k,n+CX*nvar) = ddxplus(i,j,k,n)
            do i=lo(1)+1,hi(1)-1
               bdry(i,j,k,n+CX*nvar) = ddxcen(i,j,k,n)
            enddo
            i=hi(1)
            bdry(i,j,k,n+CX*nvar) = ddxminus(i,j,k,n)
         enddo
c ::: ::: next, do z deriv
         do i=lo(1),hi(1)
            k=lo(3)
            bdry(i,j,k,n+CZ*nvar) = ddzplus(i,j,k,n)
            do k=lo(3)+1,hi(3)-1
               bdry(i,j,k,n+CZ*nvar) = ddzcen(i,j,k,n)
            enddo
            k=hi(3)
            bdry(i,j,k,n+CZ*nvar) = ddzminus(i,j,k,n)
         enddo

      enddo

      return
      end
      

c ------------------------------------------------------------------------
c :: FORT_BDIDERIVZLO : copy from Zlo face of fine fab (getting values
c ::      which were stuck there by the physical boundary condition routines)
c ::      to get the first nvar values in bdry.  The following 3*nvar values
c ::      are gotten by computing the transverse derivative.  Two formulas
c ::      are used: 1) on the interior a centered difference is used; 
c ::      2) at the end points a one-sided difference is used.  No check
c ::      is made to check that the one-sided difference does not poke out
c ::      the other side of the grid.  But if grids are that small, a lot
c ::      of things will be blowing up.
c ::      NOTE: only the transverse derivs are computed, so one of the
c ::      deriv's if filled with error values

c :: Inputs/Outputs:
c ::  bdry <= fine grid bndry data strip
c ::  lo,hi => index limits of "grid". So all action takes place in cells
c ::           outside the grid
c ::  fine  => fine grid data which has phys-bc data stuffed in it
c ::  nvar  => number of components in fine.  But you write (3+1) times as many
c ::           components of data in bdry because you compute transverse
c ::           derivatives as well
c ::  hfine => delta x in transverse direction; needed for derivs
c ------------------------------------------------------------------------
      subroutine FORT_BDIDERIVZLO( bdry, DIMS(bdry),
     &     lo, hi,
     &     fine, DIMS(fine),
     &     nvar, hfine)
      integer nvar
      REAL_T  hfine(3)
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(bdry)
      integer DIMDEC(fine)
      REAL_T  bdry(DIMV(bdry),(3+1)*nvar)
      REAL_T  fine(DIMV(fine),  nvar)

      integer i,j,k,n
      REAL_T h(3)
c ::: ------ begin functions -------
      REAL_T ddxcen
      REAL_T ddycen
      REAL_T ddzcen
      REAL_T ddxplus,ddxminus,ddyplus,ddyminus,ddzplus,ddzminus

      ddxcen(i,j,k,n) = (fine(i+1,j,k,n)-fine(i-1,j,k,n))/(2*h(1))
      ddycen(i,j,k,n) = (fine(i,j+1,k,n)-fine(i,j-1,k,n))/(2*h(2))
      ddzcen(i,j,k,n) = (fine(i,j,k+1,n)-fine(i,j,k-1,n))/(2*h(3))

      ddxplus(i,j,k,n)=(-fine(i+2,j,k,n)+4*fine(i+1,j,k,n)-3*fine(i,j,k,n))
     &     /(2*h(1))
      ddyplus(i,j,k,n)=(-fine(i,j+2,k,n)+4*fine(i,j+1,k,n)-3*fine(i,j,k,n))
     &     /(2*h(2))
      ddzplus(i,j,k,n)=(-fine(i,j,k+2,n)+4*fine(i,j,k+1,n)-3*fine(i,j,k,n))
     &     /(2*h(3))

      ddxminus(i,j,k,n)=(+fine(i-2,j,k,n)-4*fine(i-1,j,k,n)+3*fine(i,j,k,n))
     &     /(2*h(1))
      ddyminus(i,j,k,n)=(+fine(i,j-2,k,n)-4*fine(i,j-1,k,n)+3*fine(i,j,k,n))
     &     /(2*h(2))
      ddzminus(i,j,k,n)=(+fine(i,j,k-2,n)-4*fine(i,j,k-1,n)+3*fine(i,j,k,n))
     &     /(2*h(3))
c ::: ------ end   functions -------

      k = lo(3)-1
      h(1) = hfine(1)
      h(2) = hfine(2)
      h(3) = 0

      do n=1,nvar
c ::: ::: copy boundary values over
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               bdry(i,j,k,n) = fine(i,j,k,n)
            enddo
         enddo

c ::: ::: do stinking calculation of derivativ using centered difference
c ::: ::: this ONLY makes sense for DIRICHLET BC.  For other BC,
c ::: ::: we will not use these values.

c ::: ::: First, do y deriv
         do i=lo(1),hi(1)
            j=lo(2)
            bdry(i,j,k,n+CY*nvar) = ddyplus(i,j,k,n)
            do j=lo(2)+1,hi(2)-1
               bdry(i,j,k,n+CY*nvar) = ddycen(i,j,k,n)
            enddo
            j=hi(2)
            bdry(i,j,k,n+CY*nvar) = ddyminus(i,j,k,n)
         enddo
c ::: ::: next, do x deriv
         do j=lo(2),hi(2)
            i=lo(1)
            bdry(i,j,k,n+CX*nvar) = ddxplus(i,j,k,n)
            do i=lo(1)+1,hi(1)-1
               bdry(i,j,k,n+CX*nvar) = ddxcen(i,j,k,n)
            enddo
            i=hi(1)
            bdry(i,j,k,n+CX*nvar) = ddxminus(i,j,k,n)
         enddo

      enddo

      return
      end
      
      

c ------------------------------------------------------------------------
c :: FORT_BDIDERIVZHI : copy from Zhi face of fine fab (getting values
c ::      which were stuck there by the physical boundary condition routines)
c ::      to get the first nvar values in bdry.  The following 3*nvar values
c ::      are gotten by computing the transverse derivative.  Two formulas
c ::      are used: 1) on the interior a centered difference is used; 
c ::      2) at the end points a one-sided difference is used.  No check
c ::      is made to check that the one-sided difference does not poke out
c ::      the other side of the grid.  But if grids are that small, a lot
c ::      of things will be blowing up.
c ::      NOTE: only the transverse derivs are computed, so one of the
c ::      deriv's if filled with error values

c :: Inputs/Outputs:
c ::  bdry <= fine grid bndry data strip
c ::  lo,hi => index limits of "grid". So all action takes place in cells
c ::           outside the grid
c ::  fine  => fine grid data which has phys-bc data stuffed in it
c ::  nvar  => number of components in fine.  But you write (3+1) times as many
c ::           components of data in bdry because you compute transverse
c ::           derivatives as well
c ::  hfine => delta x in transverse direction; needed for derivs
c ------------------------------------------------------------------------
      subroutine FORT_BDIDERIVZHI( bdry, DIMS(bdry),
     &     lo, hi,
     &     fine, DIMS(fine),
     &     nvar, hfine)
      integer nvar
      REAL_T  hfine(3)
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(bdry)
      integer DIMDEC(fine)
      REAL_T  bdry(DIMV(bdry),(3+1)*nvar)
      REAL_T  fine(DIMV(fine),  nvar)

      integer i,j,k,n
      REAL_T h(3)
c ::: ------ begin functions -------
      REAL_T ddxcen
      REAL_T ddycen
      REAL_T ddzcen
      REAL_T ddxplus,ddxminus,ddyplus,ddyminus,ddzplus,ddzminus

      ddxcen(i,j,k,n) = (fine(i+1,j,k,n)-fine(i-1,j,k,n))/(2*h(1))
      ddycen(i,j,k,n) = (fine(i,j+1,k,n)-fine(i,j-1,k,n))/(2*h(2))
      ddzcen(i,j,k,n) = (fine(i,j,k+1,n)-fine(i,j,k-1,n))/(2*h(3))

      ddxplus(i,j,k,n)=(-fine(i+2,j,k,n)+4*fine(i+1,j,k,n)-3*fine(i,j,k,n))
     &     /(2*h(1))
      ddyplus(i,j,k,n)=(-fine(i,j+2,k,n)+4*fine(i,j+1,k,n)-3*fine(i,j,k,n))
     &     /(2*h(2))
      ddzplus(i,j,k,n)=(-fine(i,j,k+2,n)+4*fine(i,j,k+1,n)-3*fine(i,j,k,n))
     &     /(2*h(3))

      ddxminus(i,j,k,n)=(+fine(i-2,j,k,n)-4*fine(i-1,j,k,n)+3*fine(i,j,k,n))
     &     /(2*h(1))
      ddyminus(i,j,k,n)=(+fine(i,j-2,k,n)-4*fine(i,j-1,k,n)+3*fine(i,j,k,n))
     &     /(2*h(2))
      ddzminus(i,j,k,n)=(+fine(i,j,k-2,n)-4*fine(i,j,k-1,n)+3*fine(i,j,k,n))
     &     /(2*h(3))
c ::: ------ end   functions -------

      k = hi(3)+1
      h(1) = hfine(1)
      h(2) = hfine(2)
      h(3) = 0

      do n=1,nvar
c ::: ::: copy boundary values over
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               bdry(i,j,k,n) = fine(i,j,k,n)
            enddo
         enddo

c ::: ::: do stinking calculation of derivativ using centered difference
c ::: ::: this ONLY makes sense for DIRICHLET BC.  For other BC,
c ::: ::: we will not use these values.

c ::: ::: First, do y deriv
         do i=lo(1),hi(1)
            j=lo(2)
            bdry(i,j,k,n+CY*nvar) = ddyplus(i,j,k,n)
            do j=lo(2)+1,hi(2)-1
               bdry(i,j,k,n+CY*nvar) = ddycen(i,j,k,n)
            enddo
            j=hi(2)
            bdry(i,j,k,n+CY*nvar) = ddyminus(i,j,k,n)
         enddo
c ::: ::: next, do x deriv
         do j=lo(2),hi(2)
            i=lo(1)
            bdry(i,j,k,n+CX*nvar) = ddxplus(i,j,k,n)
            do i=lo(1)+1,hi(1)-1
               bdry(i,j,k,n+CX*nvar) = ddxcen(i,j,k,n)
            enddo
            i=hi(1)
            bdry(i,j,k,n+CX*nvar) = ddxminus(i,j,k,n)
         enddo

      enddo

      return
      end
      




#undef SDIM
#undef NUMDERIV
#undef DX
#undef DY
#undef DZ
#undef D2X
#undef D2Y
#undef D2Z
#undef DXY
#undef DXZ
#undef DYZ
