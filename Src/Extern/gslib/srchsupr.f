      subroutine srchsupr(xloc,yloc,zloc,radsqd,irot,MAXROT,rotmat,
     +                    nsbtosr,ixsbtosr,iysbtosr,izsbtosr,noct,nd,
     +                    x,y,z,tmp,nisb,nxsup,xmnsup,xsizsup,
     +                    nysup,ymnsup,ysizsup,nzsup,zmnsup,zsizsup,
     +                    nclose,close,infoct)
c-----------------------------------------------------------------------
c
c              Search Within Super Block Search Limits
c              ***************************************
c
c
c This subroutine searches through all the data that have been tagged in
c the super block subroutine.  The close data are passed back in the
c index array "close".  An octant search is allowed.
c
c
c
c INPUT VARIABLES:
c
c   xloc,yloc,zloc   location of point being estimated/simulated
c   radsqd           squared search radius
c   irot             index of the rotation matrix for searching
c   MAXROT           size of rotation matrix arrays
c   rotmat           rotation matrices
c   nsbtosr          Number of super blocks to search
c   ixsbtosr         X offsets for super blocks to search
c   iysbtosr         Y offsets for super blocks to search
c   izsbtosr         Z offsets for super blocks to search
c   noct             If >0 then data will be partitioned into octants
c   nd               Number of data
c   x(nd)            X coordinates of the data
c   y(nd)            Y coordinates of the data
c   z(nd)            Z coordinates of the data
c   tmp(nd)          Temporary storage to keep track of the squared
c                      distance associated with each data
c   nisb()                Array with cumulative number of data in each
c                           super block.
c   nxsup,xmnsup,xsizsup  Definition of the X super block grid
c   nysup,ymnsup,ysizsup  Definition of the X super block grid
c   nzsup,zmnsup,zsizsup  Definition of the X super block grid
c
c
c
c OUTPUT VARIABLES:
c
c   nclose           Number of close data
c   close()          Index of close data
c   infoct           Number of informed octants (only computes if
c                      performing an octant search)
c
c
c
c EXTERNAL REFERENCES:
c
c   sqdist           Computes anisotropic squared distance
c   sortem           Sorts multiple arrays in ascending order
c
c
c
c-----------------------------------------------------------------------
      implicit integer (i-n)
      implicit double precision (a-h,o-z)

      real*8  x(*),y(*),z(*),tmp(*),close(*)
      real*8  rotmat(MAXROT,3,3),hsqd,sqdist
      integer nisb(*),inoct(8)
      integer ixsbtosr(*),iysbtosr(*),izsbtosr(*)
      logical inflag
c
c Determine the super block location of point being estimated:
c
      call getindx(nxsup,xmnsup,xsizsup,xloc,ix,inflag)
      call getindx(nysup,ymnsup,ysizsup,yloc,iy,inflag)
      call getindx(nzsup,zmnsup,zsizsup,zloc,iz,inflag)
c
c Loop over all the possible Super Blocks:
c
      nclose = 0
      do 1 isup=1,nsbtosr
c
c Is this super block within the grid system:
c
            ixsup = ix + ixsbtosr(isup)
            iysup = iy + iysbtosr(isup)
            izsup = iz + izsbtosr(isup)
            if(ixsup.le.0.or.ixsup.gt.nxsup.or.
     +         iysup.le.0.or.iysup.gt.nysup.or.
     +         izsup.le.0.or.izsup.gt.nzsup) go to 1
c
c Figure out how many samples in this super block:
c
            ii = ixsup + (iysup-1)*nxsup + (izsup-1)*nxsup*nysup
            if(ii.eq.1) then
                  nums = nisb(ii)
                  i    = 0
            else
                  nums = nisb(ii) - nisb(ii-1)
                  i    = nisb(ii-1)
            endif
c
c Loop over all the data in this super block:
c
            do 2 ii=1,nums
                  i = i + 1
c
c Check squared distance:
c
                  hsqd = sqdist(xloc,yloc,zloc,x(i),y(i),z(i),irot,
     +                          MAXROT,rotmat)
                  if(real(hsqd).gt.radsqd) go to 2
c
c Accept this sample:
c
                  nclose = nclose + 1
                  close(nclose) = real(i)
                  tmp(nclose)  = real(hsqd)
 2          continue
 1    continue
c
c Sort the nearby samples by distance to point being estimated:
c
      call sortem(1,nclose,tmp,1,close,c,d,e,f,g,h)
c
c If we aren't doing an octant search then just return:
c
      if(noct.le.0) return
c
c PARTITION THE DATA INTO OCTANTS:
c
      do i=1,8
            inoct(i) = 0
      end do
c
c Now pick up the closest samples in each octant:
c
      nt = 8*noct
      na = 0
      do j=1,nclose
            i  = int(close(j))
            h  = tmp(j)
            dx = x(i) - xloc
            dy = y(i) - yloc
            dz = z(i) - zloc
            if(dz.lt.0.) go to 5
            iq=4
            if(dx.le.0.d0 .and. dy.gt.0.d0) iq=1
            if(dx.gt.0.d0 .and. dy.ge.0.d0) iq=2
            if(dx.lt.0.d0 .and. dy.le.0.d0) iq=3
            go to 6
 5          iq=8
            if(dx.le.0.d0 .and. dy.gt.0.d0) iq=5
            if(dx.gt.0.d0 .and. dy.ge.0.d0) iq=6
            if(dx.lt.0.d0 .and. dy.le.0.d0) iq=7
 6          continue
            inoct(iq) = inoct(iq) + 1
c
c Keep this sample if the maximum has not been exceeded:
c
            if(inoct(iq).le.noct) then
                  na = na + 1
                  close(na) = i
                  tmp(na)   = h
                  if(na.eq.nt) go to 7
            endif
      end do
c
c End of data selection. Compute number of informed octants and return:
c
 7    nclose = na
      infoct = 0
      do i=1,8
            if(inoct(i).gt.0) infoct = infoct + 1
      end do
c
c Finished:
c
      return
      end
