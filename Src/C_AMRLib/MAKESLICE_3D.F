c
c $Id: MAKESLICE_3D.F,v 1.6 1998-03-31 00:24:00 lijewski Exp $
c

#undef  BL_LANG_CC         
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif
        
#include "REAL.H"
#include "CONSTANTS.H"
#include "MAKESLICE_F.H"

#define DIMENSION(lo,hi) lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)
#define SDIM 3         
#define npart 20


      subroutine FORT_MAKEMASKSLICE(data1,dlo1,dhi1,
     &  data2,dlo2,dhi2,data_min,data_max,slice_val,
     &  slice,dx)
      integer dlo1(SDIM),dhi1(SDIM)
      integer dlo2(SDIM),dhi2(SDIM)
      integer data1(DIMENSION(dlo1,dhi1))
      integer data2(DIMENSION(dlo2,dhi2))
      REAL_T dx(SDIM)
      REAL_T data_min,data_max,slice_val
      integer slice

      integer index,idx1,idx2,found,i,j,k,kval,kval1,kval2
      integer leftval,rightval
      REAL_T zval,tval,floatdata

      if ((dlo2(3).ne.0).or.(dhi2(3).ne.0)) then
         write (6,*) "second data must be thin!"
      end if

      if (slice.eq.0) then
         idx1=1
         idx2=2
         index=3
      else if (slice.eq.1) then
         idx1=1
         idx2=3
         index=2
      else
         idx1=2
         idx2=3
         index=1
      end if

      if ( (dlo1(idx1).ne.dlo2(1)).or.(dlo1(idx2).ne.dlo2(2)).or.
     &     (dhi1(idx1).ne.dhi2(1)).or.(dhi1(idx2).ne.dhi2(2)) ) then
         write (6,*) "thin data must jive with 3d data!"
      end if

      found=0
      do k=dlo1(index)-1,dhi1(index)
         zval=(k+half)*dx(index)
         if ((slice_val.ge.zval).and.(slice_val.le.zval+dx(index)).and.
     &        (found.eq.0)) then
            found=1
            kval=k 
         end if
      end do
      if (found.eq.0) then
         write (6,*) "slice data out of range!"
      end if
      kval1=kval
      kval2=kval+1
      if (kval1.lt.dlo1(index)) then
         kval1=kval2
      end if
      if (kval2.gt.dhi1(index)) then
         kval2=kval1
      end if

      do i=dlo1(idx1),dhi1(idx1)
         do j=dlo1(idx2),dhi1(idx2)
            zval=(kval+half)*dx(index)
            if (slice.eq.0) then
               leftval=data1(i,j,kval1)
               rightval=data1(i,j,kval2)
            else if (slice.eq.1) then
               leftval=data1(i,kval1,j)
               rightval=data1(i,kval2,j)
            else if (slice.eq.2) then
               leftval=data1(kval1,i,j)
               rightval=data1(kval2,i,j)
            end if
            tval=(slice_val-zval)/dx(index)
            floatdata=leftval*(one-tval)+rightval*tval
            
            if (floatdata.gt.zero) then
               data2(i,j,0)=1
            else
               data2(i,j,0)=0
            end if
         end do
      end do

      end

      subroutine FORT_MAKESLICE(data1,dlo1,dhi1,
     &  data2,dlo2,dhi2,data_min,data_max,slice_val,
     &  slice,dx)

      integer dlo1(SDIM),dhi1(SDIM)
      integer dlo2(SDIM),dhi2(SDIM)
      REAL_T data1(DIMENSION(dlo1,dhi1))
      REAL_T data2(DIMENSION(dlo2,dhi2))
      REAL_T dx(SDIM)
      REAL_T data_min,data_max,slice_val
      integer slice

      integer index,idx1,idx2,found,i,j,k,kval,kval1,kval2
      REAL_T zval,leftval,rightval,tval

      if ((dlo2(3).ne.0).or.(dhi2(3).ne.0)) then
         write (6,*) "second data must be thin!"
      end if

      if (slice.eq.0) then
         idx1=1
         idx2=2
         index=3
      else if (slice.eq.1) then
         idx1=1
         idx2=3
         index=2
      else
         idx1=2
         idx2=3
         index=1
      end if

      if ( (dlo1(idx1).ne.dlo2(1)).or.(dlo1(idx2).ne.dlo2(2)).or.
     &     (dhi1(idx1).ne.dhi2(1)).or.(dhi1(idx2).ne.dhi2(2)) ) then
         write (6,*) "thin data must jive with 3d data!"
      end if

      found=0
      do k=dlo1(index)-1,dhi1(index)
         zval=(k+half)*dx(index)
         if ((slice_val.ge.zval).and.(slice_val.le.zval+dx(index)).and.
     &        (found.eq.0)) then
            found=1
            kval=k 
         end if
      end do
      if (found.eq.0) then
         write (6,*) "slice data out of range!"
      end if
      kval1=kval
      kval2=kval+1
      if (kval1.lt.dlo1(index)) then
         kval1=kval2
      end if
      if (kval2.gt.dhi1(index)) then
         kval2=kval1
      end if

      do i=dlo1(idx1),dhi1(idx1)
         do j=dlo1(idx2),dhi1(idx2)
            zval=(kval+half)*dx(index)
            if (slice.eq.0) then
               leftval=data1(i,j,kval1)
               rightval=data1(i,j,kval2)
            else if (slice.eq.1) then
               leftval=data1(i,kval1,j)
               rightval=data1(i,kval2,j)
            else if (slice.eq.2) then
               leftval=data1(kval1,i,j)
               rightval=data1(kval2,i,j)
            end if
            tval=(slice_val-zval)/dx(index)
            data2(i,j,0)=leftval*(one-tval)+rightval*tval
         end do
      end do

      end

