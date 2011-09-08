      subroutine picksup(nxsup,xsizsup,nysup,ysizsup,nzsup,zsizsup,
     +                   irot,MAXROT,rotmat,radsqd,nsbtosr,ixsbtosr,
     +                   iysbtosr,izsbtosr)
c-----------------------------------------------------------------------
c
c             Establish Which Super Blocks to Search
c             **************************************
c
c This subroutine establishes which super blocks must be searched given
c that a point being estimated/simulated falls within a super block
c centered at 0,0,0.
c
c
c
c INPUT VARIABLES:
c
c   nxsup,xsizsup    Definition of the X super block grid
c   nysup,ysizsup    Definition of the Y super block grid
c   nzsup,zsizsup    Definition of the Z super block grid
c   irot             index of the rotation matrix for searching
c   MAXROT           size of rotation matrix arrays
c   rotmat           rotation matrices
c   radsqd           squared search radius
c
c
c
c OUTPUT VARIABLES:
c
c   nsbtosr          Number of super blocks to search
c   ixsbtosr         X offsets for super blocks to search
c   iysbtosr         Y offsets for super blocks to search
c   izsbtosr         Z offsets for super blocks to search
c
c
c
c EXTERNAL REFERENCES:
c
c   sqdist           Computes anisotropic squared distance
c
c
c
c-----------------------------------------------------------------------

      implicit integer (i-n)
      implicit double precision (a-h,o-z)

      real*8  rotmat(MAXROT,3,3),hsqd,sqdist,shortest
      integer ixsbtosr(*),iysbtosr(*),izsbtosr(*)
c
c MAIN Loop over all possible super blocks:
c
      nsbtosr = 0
      do i=-(nxsup-1),(nxsup-1)
      do j=-(nysup-1),(nysup-1)
      do k=-(nzsup-1),(nzsup-1)
            xo = real(i)*xsizsup
            yo = real(j)*ysizsup
            zo = real(k)*zsizsup
c
c Find the closest distance between the corners of the super blocks:
c
            shortest = 1.0e21
            do i1=-1,1
            do j1=-1,1
            do k1=-1,1
                  do i2=-1,1
                  do j2=-1,1
                  do k2=-1,1
                        if(i1.ne.0.and.j1.ne.0.and.k1.ne.0.and.
     +                     i2.ne.0.and.j2.ne.0.and.k2.ne.0) then
                              xdis = real(i1-i2)*0.5*xsizsup + xo
                              ydis = real(j1-j2)*0.5*ysizsup + yo
                              zdis = real(k1-k2)*0.5*zsizsup + zo
                              hsqd = sqdist(0.d0,0.d0,0.d0,xdis,ydis,
     +                                      ydis,irot,MAXROT,rotmat)
                              if(hsqd.lt.shortest) shortest = hsqd
                        end if
                  end do
                  end do
                  end do
            end do
            end do
            end do
c
c Keep this super block if it is close enoutgh:
c
            if(real(shortest).le.radsqd) then
                  nsbtosr = nsbtosr + 1
                  ixsbtosr(nsbtosr) = i
                  iysbtosr(nsbtosr) = j
                  izsbtosr(nsbtosr) = k
            end if
      end do
      end do
      end do
c
c Finished:
c
      return
      end
