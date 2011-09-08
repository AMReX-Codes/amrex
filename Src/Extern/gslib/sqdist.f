      real*8 function sqdist(x1,y1,z1,x2,y2,z2,ind,MAXROT,rotmat)
c-----------------------------------------------------------------------
c
c    Squared Anisotropic Distance Calculation Given Matrix Indicator
c    ***************************************************************
c
c This routine calculates the anisotropic distance between two points
c  given the coordinates of each point and a definition of the
c  anisotropy.
c
c
c INPUT VARIABLES:
c
c   x1,y1,z1         Coordinates of first point
c   x2,y2,z2         Coordinates of second point
c   ind              The rotation matrix to use
c   MAXROT           The maximum number of rotation matrices dimensioned
c   rotmat           The rotation matrices
c
c
c
c OUTPUT VARIABLES:
c
c   sqdist           The squared distance accounting for the anisotropy
c                      and the rotation of coordinates (if any).
c
c
c NO EXTERNAL REFERENCES
c
c
c-----------------------------------------------------------------------

      implicit integer (i-n)
      implicit double precision (a-h,o-z)

      real*8 rotmat(MAXROT,3,3),cont,dx,dy,dz
c
c Compute component distance vectors and the squared distance:
c
      dx = dble(x1 - x2)
      dy = dble(y1 - y2)
      dz = dble(z1 - z2)
      sqdist = 0.d0
      do i=1,3
            cont   = rotmat(ind,i,1) * dx
     +             + rotmat(ind,i,2) * dy
     +             + rotmat(ind,i,3) * dz
            sqdist = sqdist + cont * cont
      end do
      return
      end
