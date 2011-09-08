      subroutine setsupr(nx,xmn,xsiz,ny,ymn,ysiz,nz,zmn,zsiz,nd,x,y,z,
     +                   vr,tmp,nsec,sec1,sec2,sec3,MAXSBX,MAXSBY,
     +                   MAXSBZ,nisb,nxsup,xmnsup,xsizsup,nysup,ymnsup,
     +                   ysizsup,nzsup,zmnsup,zsizsup)
c-----------------------------------------------------------------------
c
c           Establish Super Block Search Limits and Sort Data
c           *************************************************
c
c This subroutine sets up a 3-D "super block" model and orders the data
c by super block number.  The limits of the super block is set to the
c minimum and maximum limits of the grid; data outside are assigned to
c the nearest edge block.
c
c The idea is to establish a 3-D block network that contains all the
c relevant data.  The data are then sorted by their index location in
c the search network, i.e., the index location is given after knowing
c the block index in each coordinate direction (ix,iy,iz):
c          ii = (iz-1)*nxsup*nysup + (iy-1)*nxsup + ix
c An array, the same size as the number of super blocks, is constructed
c that contains the cumulative number of data in the model.  With this
c array it is easy to quickly check what data are located near any given
c location.
c
c
c
c INPUT VARIABLES:
c
c   nx,xmn,xsiz      Definition of the X grid being considered
c   ny,ymn,ysiz      Definition of the Y grid being considered
c   nz,zmn,zsiz      Definition of the Z grid being considered
c   nd               Number of data
c   x(nd)            X coordinates of the data
c   y(nd)            Y coordinates of the data
c   z(nd)            Z coordinates of the data
c   vr(nd)           Variable at each location.
c   tmp(nd)          Temporary storage to keep track of the super block
c                      index associated to each data (uses the same
c                      storage already allocated for the simulation)
c   nsec             Number of secondary variables to carry with vr
c   sec1(nd)         First secondary variable (if nsec >= 1)
c   sec2(nd)         Second secondary variable (if nsec >= 2)
c   sec3(nd)         Third secondary variable (if nsec = 3)
c   MAXSB[X,Y,Z]     Maximum size of super block network
c
c
c
c OUTPUT VARIABLES:
c
c   nisb()                Array with cumulative number of data in each
c                           super block.
c   nxsup,xmnsup,xsizsup  Definition of the X super block grid
c   nysup,ymnsup,ysizsup  Definition of the Y super block grid
c   nzsup,zmnsup,zsizsup  Definition of the Z super block grid
c
c
c
c EXTERNAL REFERENCES:
c
c   sortem           Sorting routine to sort the data
c
c
c
c-----------------------------------------------------------------------
      implicit integer (i-n)
      implicit double precision (a-h,o-z)

      real*8  x(*),y(*),z(*),vr(*),tmp(*),sec1(*),sec2(*),sec3(*)
      integer nisb(*)
      logical inflag
c
c Establish the number and size of the super blocks:
c
      nxsup   = min(nx,MAXSBX)
      nysup   = min(ny,MAXSBY)
      nzsup   = min(nz,MAXSBZ)
      xsizsup = real(nx)*xsiz/real(nxsup)
      ysizsup = real(ny)*ysiz/real(nysup)
      zsizsup = real(nz)*zsiz/real(nzsup)
      xmnsup  = (xmn-0.5*xsiz)+0.5*xsizsup
      ymnsup  = (ymn-0.5*ysiz)+0.5*ysizsup
      zmnsup  = (zmn-0.5*zsiz)+0.5*zsizsup
c
c Initialize the extra super block array to zeros:
c
      do i=1,nxsup*nysup*nzsup
            nisb(i) = 0
      end do
c
c Loop over all the data assigning the data to a super block and
c accumulating how many data are in each super block:
c
      do i=1,nd
            call getindx(nxsup,xmnsup,xsizsup,x(i),ix,inflag)
            call getindx(nysup,ymnsup,ysizsup,y(i),iy,inflag)
            call getindx(nzsup,zmnsup,zsizsup,z(i),iz,inflag)
            ii = ix + (iy-1)*nxsup + (iz-1)*nxsup*nysup
            tmp(i)   = ii
            nisb(ii) = nisb(ii) + 1
      end do
c
c Sort the data by ascending super block number:
c
      nsort = 4 + nsec
      call sortem(1,nd,tmp,nsort,x,y,z,vr,sec1,sec2,sec3)
c
c Set up array nisb with the starting address of the block data:
c
      do i=1,(nxsup*nysup*nzsup-1)
            nisb(i+1) = nisb(i) + nisb(i+1)
      end do
c
c Finished:
c
      return
      end
