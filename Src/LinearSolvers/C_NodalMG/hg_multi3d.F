

#define dims(a) a##l0,a##h0,a##l1,a##h1,a##l2,a##h2
#define dimdec(a) a##l0:a##h0,a##l1:a##h1,a##l2:a##h2
c#define dimdecu(a) (a##h0-a##l0+1)*(a##h1-a##l1+1)*(a##h2-a##l2+1)
#define dimdecu(a) *

#ifndef HG_CONSTANT

#  ifndef HG_TERRAIN
#    define HG_SIGMA_NODE
#  endif

#ifdef HG_TERRAIN 

c twenty-seven point terrain version:

#define MXXC(ii,jj,kk) (4.d0 *  cdst(ii,j,k)   + \
                        2.d0 * (cdst(ii,jj,k)  - cdst(i,jj,k) + \
                                cdst(ii,j,kk)  - cdst(i,j,kk)) + \
                               (cdst(ii,jj,kk) - cdst(i,jj,kk)))

#define MYYC(ii,jj,kk) (4.d0 *  cdst(i,jj,k)   + \
                        2.d0 * (cdst(ii,jj,k)  - cdst(ii,j,k) + \
                                cdst(i,jj,kk)  - cdst(i,j,kk)) + \
                               (cdst(ii,jj,kk) - cdst(ii,j,kk)))

#define MZZC(ii,jj,kk) (4.d0 *  cdst(i,j,kk)   + \
                        2.d0 * (cdst(ii,j,kk)  - cdst(ii,j,k) + \
                                cdst(i,jj,kk)  - cdst(i,jj,k)) + \
                               (cdst(ii,jj,kk) - cdst(ii,jj,k)))

#define MXZC(ii,jj,kk) (6.d0 *  cdst(ii,j,kk)  + \
                        3.d0 * (cdst(ii,jj,kk) - cdst(i,jj,k)))

#define MYZC(ii,jj,kk) (6.d0 *  cdst(i,jj,kk)  + \
                        3.d0 * (cdst(ii,jj,kk) - cdst(ii,j,k)))

#define MXXF(i1,j1,k1,ii,jj,kk) (4.d0 *  fdst(ii,j1,k1)   + \
                                 2.d0 * (fdst(ii,jj,k1) - fdst(i1,jj,k1) + \
                                         fdst(ii,j1,kk) - fdst(i1,j1,kk)) + \
                                        (fdst(ii,jj,kk) - fdst(i1,jj,kk)))

#define MYYF(i1,j1,k1,ii,jj,kk) (4.d0 *  fdst(i1,jj,k1)   + \
                                 2.d0 * (fdst(ii,jj,k1) - fdst(ii,j1,k1) + \
                                         fdst(i1,jj,kk) - fdst(i1,j1,kk)) + \
                                        (fdst(ii,jj,kk) - fdst(ii,j1,kk)))

#define MZZF(i1,j1,k1,ii,jj,kk) (4.d0 *  fdst(i1,j1,kk)   + \
                                 2.d0 * (fdst(ii,j1,kk) - fdst(ii,j1,k1) + \
                                         fdst(i1,jj,kk) - fdst(i1,jj,k1)) + \
                                        (fdst(ii,jj,kk) - fdst(ii,jj,k1)))

#define MXZF(i1,j1,k1,ii,jj,kk) (6.d0 *  fdst(ii,j1,kk)  + \
                                 3.d0 * (fdst(ii,jj,kk) - fdst(i1,jj,k1)))

#define MYZF(i1,j1,k1,ii,jj,kk) (6.d0 *  fdst(i1,jj,kk)  + \
                                 3.d0 * (fdst(ii,jj,kk) - fdst(ii,j1,k1)))


c Note---assumes fdst linearly interpolated from cdst along face
      subroutine hgfres(res, dims(res),
     @ src,  dims(src),
     @ fdst, dims(fdst),
     @ cdst, dims(cdst),
     @ sigmaf, dims(sf),
     @ sigmac, dims(sc),
     @ dims(reg),
     @ ir, jr, kr, idim, idir)
      integer dims(res)
      integer dims(src)
      integer dims(fdst)
      integer dims(cdst)
      integer dims(sf)
      integer dims(sc)
      integer dims(reg)
      real*8 res(dimdec(res))
      real*8 src(dimdec(src))
      real*8 fdst(dimdec(fdst))
      real*8 cdst(dimdec(cdst))
      real*8 sigmaf(dimdec(sf), 5)
      real*8 sigmac(dimdec(sc), 5)
      integer ir, jr, kr, idim, idir
      real*8 fac00, fac0, fac1, fac2, tmp, cen
      integer i, j, k, is, js, ks, l, m, n

      fac00 = 1.d0 / 36.d0

      if (idim .eq. 0) then
         i = regl0
         if (idir .eq. 1) then
            is = i - 1
         else
            is = i
         endif
         do 10 k = regl2, regh2
            do 10 j = regl1, regh1

                  tmp = sigmac(is ,j-1,k-1,1) * MXXC(i-idir,j-1,k-1)
     @                 +sigmac(is ,j-1,k  ,1) * MXXC(i-idir,j-1,k+1)
     @                 +sigmac(is ,j  ,k-1,1) * MXXC(i-idir,j+1,k-1)
     @                 +sigmac(is ,j  ,k  ,1) * MXXC(i-idir,j+1,k+1)
                  tmp = tmp +
     @                 +sigmac(is ,j-1,k-1,2) * MYYC(i-idir,j-1,k-1)
     @                 +sigmac(is ,j-1,k  ,2) * MYYC(i-idir,j-1,k+1)
     @                 +sigmac(is ,j  ,k-1,2) * MYYC(i-idir,j+1,k-1)
     @                 +sigmac(is ,j  ,k  ,2) * MYYC(i-idir,j+1,k+1)
                  tmp = tmp +
     @                 +sigmac(is ,j-1,k-1,3) * MZZC(i-idir,j-1,k-1)
     @                 +sigmac(is ,j-1,k  ,3) * MZZC(i-idir,j-1,k+1)
     @                 +sigmac(is ,j  ,k-1,3) * MZZC(i-idir,j+1,k-1)
     @                 +sigmac(is ,j  ,k  ,3) * MZZC(i-idir,j+1,k+1)
                  tmp = tmp + idir * (
     @                 -sigmac(is ,j-1,k-1,4) * MXZC(i-idir,j-1,k-1)
     @                 +sigmac(is ,j-1,k  ,4) * MXZC(i-idir,j-1,k+1)
     @                 -sigmac(is ,j  ,k-1,4) * MXZC(i-idir,j+1,k-1)
     @                 +sigmac(is ,j  ,k  ,4) * MXZC(i-idir,j+1,k+1) )
                  tmp = tmp
     @                 -sigmac(is ,j-1,k-1,5) * MYZC(i-idir,j-1,k-1)
     @                 +sigmac(is ,j-1,k  ,5) * MYZC(i-idir,j-1,k+1)
     @                 +sigmac(is ,j  ,k-1,5) * MYZC(i-idir,j+1,k-1)
     @                 -sigmac(is ,j  ,k  ,5) * MYZC(i-idir,j+1,k+1)

              cen        = 
     @            4.d0 * ( sigmac(is,j-1,k-1,1) + sigmac(is,j-1,k,1)
     @                    +sigmac(is,j  ,k-1,1) + sigmac(is,j  ,k,1)
     @                    +sigmac(is,j-1,k-1,2) + sigmac(is,j-1,k,2)
     @                    +sigmac(is,j  ,k-1,2) + sigmac(is,j  ,k,2)
     @                    +sigmac(is,j-1,k-1,3) + sigmac(is,j-1,k,3)
     @                    +sigmac(is,j  ,k-1,3) + sigmac(is,j  ,k,3)) 
     @   + idir * 6.d0 * (-sigmac(is,j-1,k-1,4) + sigmac(is,j-1,k,4)
     @                    -sigmac(is,j  ,k-1,4) + sigmac(is,j  ,k,4))
     @          + 6.d0 * ( sigmac(is,j-1,k  ,5) - sigmac(is,j-1,k-1,5)
     @                    +sigmac(is,j  ,k-1,5) - sigmac(is,j  ,k,5)) 

               res(i*ir,j*jr,k*kr) =
     @            src(i*ir,j*jr,k*kr) - fac00 * (tmp - cen * cdst(i,j,k))
 10      continue

         fac0 = fac00 / (jr * kr)
         i = i * ir
         if (idir .eq. 1) then
            is = i
         else
            is = i - 1
         endif
         do 20 l = -(kr-1), kr-1
            fac2 = (kr-abs(l)) * fac0
            do 20 n = -(jr-1), jr-1
               fac1 = (jr-abs(n)) * fac2
               do 20 k = kr*regl2, kr*regh2, kr
                  do 20 j = jr*regl1, jr*regh1, jr

            tmp = sigmaf(is,j+n-1,k+l-1,1) * MXXF(i,j+n,k+l,i+idir,j+n-1,k+l-1)
     @           +sigmaf(is,j+n-1,k+l  ,1) * MXXF(i,j+n,k+l,i+idir,j+n-1,k+l+1)
     @           +sigmaf(is,j+n  ,k+l-1,1) * MXXF(i,j+n,k+l,i+idir,j+n+1,k+l-1)
     @           +sigmaf(is,j+n  ,k+l  ,1) * MXXF(i,j+n,k+l,i+idir,j+n+1,k+l+1)
            tmp = tmp
     @           +sigmaf(is,j+n-1,k+l-1,2) * MYYF(i,j+n,k+l,i+idir,j+n-1,k+l-1)
     @           +sigmaf(is,j+n-1,k+l  ,2) * MYYF(i,j+n,k+l,i+idir,j+n-1,k+l+1)
     @           +sigmaf(is,j+n  ,k+l-1,2) * MYYF(i,j+n,k+l,i+idir,j+n+1,k+l-1)
     @           +sigmaf(is,j+n  ,k+l  ,2) * MYYF(i,j+n,k+l,i+idir,j+n+1,k+l+1)
            tmp = tmp
     @           +sigmaf(is,j+n-1,k+l-1,3) * MZZF(i,j+n,k+l,i+idir,j+n-1,k+l-1)
     @           +sigmaf(is,j+n-1,k+l  ,3) * MZZF(i,j+n,k+l,i+idir,j+n-1,k+l+1)
     @           +sigmaf(is,j+n  ,k+l-1,3) * MZZF(i,j+n,k+l,i+idir,j+n+1,k+l-1)
     @           +sigmaf(is,j+n  ,k+l  ,3) * MZZF(i,j+n,k+l,i+idir,j+n+1,k+l+1)
            tmp = tmp + idir * (
     @            sigmaf(is,j+n-1,k+l-1,4) * MXZF(i,j+n,k+l,i+idir,j+n-1,k+l-1)
     @           -sigmaf(is,j+n-1,k+l  ,4) * MXZF(i,j+n,k+l,i+idir,j+n-1,k+l+1)
     @           +sigmaf(is,j+n  ,k+l-1,4) * MXZF(i,j+n,k+l,i+idir,j+n+1,k+l-1)
     @           -sigmaf(is,j+n  ,k+l  ,4) * MXZF(i,j+n,k+l,i+idir,j+n+1,k+l+1))
            tmp = tmp
     @           -sigmaf(is,j+n-1,k+l-1,5) * MYZF(i,j+n,k+l,i+idir,j+n-1,k+l-1)
     @           +sigmaf(is,j+n-1,k+l  ,5) * MYZF(i,j+n,k+l,i+idir,j+n-1,k+l+1)
     @           +sigmaf(is,j+n  ,k+l-1,5) * MYZF(i,j+n,k+l,i+idir,j+n+1,k+l-1)
     @           -sigmaf(is,j+n  ,k+l  ,5) * MYZF(i,j+n,k+l,i+idir,j+n+1,k+l+1)

            cen =
     @       4.d0 * (sigmaf(is,j+n-1,k+l-1,1) + sigmaf(is,j+n-1,k+l  ,1)
     @              +sigmaf(is,j+n  ,k+l-1,1) + sigmaf(is,j+n  ,k+l  ,1)
     @              +sigmaf(is,j+n-1,k+l-1,2) + sigmaf(is,j+n-1,k+l  ,2)
     @              +sigmaf(is,j+n  ,k+l-1,2) + sigmaf(is,j+n  ,k+l  ,2)
     @              +sigmaf(is,j+n-1,k+l-1,3) + sigmaf(is,j+n-1,k+l  ,3)
     @              +sigmaf(is,j+n  ,k+l-1,3) + sigmaf(is,j+n  ,k+l  ,3))
     @     + idir * 
     @       6.d0 * (-sigmaf(is,j+n-1,k+l  ,4) + sigmaf(is,j+n-1,k+l-1,4)
     @               -sigmaf(is,j+n  ,k+l  ,4) + sigmaf(is,j+n  ,k+l-1,4))
     @     + 6.d0 * ( sigmaf(is,j+n-1,k+l  ,5) - sigmaf(is,j+n-1,k+l-1,5)
     @               +sigmaf(is,j+n  ,k+l-1,5) - sigmaf(is,j+n  ,k+l  ,5))

             res(i,j,k) = res(i,j,k) - fac1 * (tmp - cen * fdst(i,j+n,k+l))
 20      continue
              
      else if (idim .eq. 1) then
         j = regl1
         if (idir .eq. 1) then
            js = j - 1
         else
            js = j
         endif
         do 30 k = regl2, regh2
            do 30 i = regl0, regh0

                  tmp = sigmac(i-1,js,k-1,1) * MXXC(i-1,j-idir,k-1)
     @                 +sigmac(i-1,js,k  ,1) * MXXC(i-1,j-idir,k+1)
     @                 +sigmac(i  ,js,k-1,1) * MXXC(i+1,j-idir,k-1)
     @                 +sigmac(i  ,js,k  ,1) * MXXC(i+1,j-idir,k+1)
                  tmp = tmp
     @                 +sigmac(i-1,js,k-1,2) * MYYC(i-1,j-idir,k-1)
     @                 +sigmac(i-1,js,k  ,2) * MYYC(i-1,j-idir,k+1)
     @                 +sigmac(i  ,js,k-1,2) * MYYC(i+1,j-idir,k-1)
     @                 +sigmac(i  ,js,k  ,2) * MYYC(i+1,j-idir,k+1)
                  tmp = tmp
     @                 +sigmac(i-1,js,k-1,3) * MZZC(i-1,j-idir,k-1)
     @                 +sigmac(i-1,js,k  ,3) * MZZC(i-1,j-idir,k+1)
     @                 +sigmac(i  ,js,k-1,3) * MZZC(i+1,j-idir,k-1)
     @                 +sigmac(i  ,js,k  ,3) * MZZC(i+1,j-idir,k+1)
                  tmp = tmp
     @                 -sigmac(i-1,js,k-1,4) * MXZC(i-1,j-idir,k-1)
     @                 +sigmac(i-1,js,k  ,4) * MXZC(i-1,j-idir,k+1)
     @                 +sigmac(i  ,js,k-1,4) * MXZC(i+1,j-idir,k-1)
     @                 -sigmac(i  ,js,k  ,4) * MXZC(i+1,j-idir,k+1)
                  tmp = tmp + idir * ( 
     @                 -sigmac(i-1,js,k-1,5) * MYZC(i-1,j-idir,k-1)
     @                 +sigmac(i-1,js,k  ,5) * MYZC(i-1,j-idir,k+1)
     @                 -sigmac(i  ,js,k-1,5) * MYZC(i+1,j-idir,k-1)
     @                 +sigmac(i  ,js,k  ,5) * MYZC(i+1,j-idir,k+1) )

               cen =
     @            4.d0 * (sigmac(i-1,js,k-1,1) + sigmac(i-1,js,k,1)
     @                   +sigmac(i  ,js,k-1,1) + sigmac(i  ,js,k,1)
     @                   +sigmac(i-1,js,k-1,2) + sigmac(i-1,js,k,2)
     @                   +sigmac(i  ,js,k-1,2) + sigmac(i  ,js,k,2)
     @                   +sigmac(i-1,js,k-1,3) + sigmac(i-1,js,k,3)
     @                   +sigmac(i  ,js,k-1,3) + sigmac(i  ,js,k,3))
     @          + 6.d0 * (sigmac(i-1,js,k  ,4) - sigmac(i-1,js,k-1,4)
     @                   +sigmac(i  ,js,k-1,4) - sigmac(i  ,js,k  ,4))
     @   + idir * 6.d0 * (-sigmac(i-1,js,k-1,5) + sigmac(i-1,js,k,5)
     @                    -sigmac(i  ,js,k-1,5) + sigmac(i  ,js,k,5))

               res(i*ir,j*jr,k*kr) =
     @            src(i*ir,j*jr,k*kr) - fac00 * (tmp - cen * cdst(i,j,k))
 30      continue

         fac0 = fac00 / (ir * kr)
         j = j * jr
         if (idir .eq. 1) then
            js = j
         else
            js = j - 1
         endif
         do 40 l = -(kr-1), kr-1
            fac2 = (kr-abs(l)) * fac0
            do 40 m = -(ir-1), ir-1
               fac1 = (ir-abs(m)) * fac2
               do 40 k = kr*regl2, kr*regh2, kr
                  do 40 i = ir*regl0, ir*regh0, ir

            tmp = sigmaf(i+m-1,js,k+l-1,1) * MXXF(i+m,j,k+l,i+m-1,j+idir,k+l-1)
     @           +sigmaf(i+m-1,js,k+l  ,1) * MXXF(i+m,j,k+l,i+m-1,j+idir,k+l+1)
     @           +sigmaf(i+m  ,js,k+l-1,1) * MXXF(i+m,j,k+l,i+m+1,j+idir,k+l-1)
     @           +sigmaf(i+m  ,js,k+l  ,1) * MXXF(i+m,j,k+l,i+m+1,j+idir,k+l+1)
            tmp = tmp
     @           +sigmaf(i+m-1,js,k+l-1,2) * MYYF(i+m,j,k+l,i+m-1,j+idir,k+l-1)
     @           +sigmaf(i+m-1,js,k+l  ,2) * MYYF(i+m,j,k+l,i+m-1,j+idir,k+l+1)
     @           +sigmaf(i+m  ,js,k+l-1,2) * MYYF(i+m,j,k+l,i+m+1,j+idir,k+l-1)
     @           +sigmaf(i+m  ,js,k+l  ,2) * MYYF(i+m,j,k+l,i+m+1,j+idir,k+l+1)
            tmp = tmp
     @           +sigmaf(i+m-1,js,k+l-1,3) * MZZF(i+m,j,k+l,i+m-1,j+idir,k+l-1)
     @           +sigmaf(i+m-1,js,k+l  ,3) * MZZF(i+m,j,k+l,i+m-1,j+idir,k+l+1)
     @           +sigmaf(i+m  ,js,k+l-1,3) * MZZF(i+m,j,k+l,i+m+1,j+idir,k+l-1)
     @           +sigmaf(i+m  ,js,k+l  ,3) * MZZF(i+m,j,k+l,i+m+1,j+idir,k+l+1)
            tmp = tmp
     @           -sigmaf(i+m-1,js,k+l-1,4) * MXZF(i+m,j,k+l,i+m-1,j+idir,k+l-1)
     @           +sigmaf(i+m-1,js,k+l  ,4) * MXZF(i+m,j,k+l,i+m-1,j+idir,k+l+1)
     @           +sigmaf(i+m  ,js,k+l-1,4) * MXZF(i+m,j,k+l,i+m+1,j+idir,k+l-1)
     @           -sigmaf(i+m  ,js,k+l  ,4) * MXZF(i+m,j,k+l,i+m+1,j+idir,k+l+1)
            tmp = tmp + idir * (
     @            sigmaf(i+m-1,js,k+l-1,5) * MYZF(i+m,j,k+l,i+m-1,j+idir,k+l-1)
     @           -sigmaf(i+m-1,js,k+l  ,5) * MYZF(i+m,j,k+l,i+m-1,j+idir,k+l+1)
     @           +sigmaf(i+m  ,js,k+l-1,5) * MYZF(i+m,j,k+l,i+m+1,j+idir,k+l-1)
     @           -sigmaf(i+m  ,js,k+l  ,5) * MYZF(i+m,j,k+l,i+m+1,j+idir,k+l+1))

         cen =
     @    4.d0 * (sigmaf(i+m-1,js,k+l-1,1) + sigmaf(i+m-1,js,k+l  ,1)
     @           +sigmaf(i+m  ,js,k+l-1,1) + sigmaf(i+m  ,js,k+l  ,1)
     @           +sigmaf(i+m-1,js,k+l-1,2) + sigmaf(i+m-1,js,k+l  ,2)
     @           +sigmaf(i+m  ,js,k+l-1,2) + sigmaf(i+m  ,js,k+l  ,2)
     @           +sigmaf(i+m-1,js,k+l-1,3) + sigmaf(i+m-1,js,k+l  ,3)
     @           +sigmaf(i+m  ,js,k+l-1,3) + sigmaf(i+m  ,js,k+l  ,3))
     @   +6.d0 * (sigmaf(i+m-1,js,k+l  ,4) - sigmaf(i+m-1,js,k+l-1,4)
     @           +sigmaf(i+m  ,js,k+l-1,4) - sigmaf(i+m  ,js,k+l  ,4))
     @   + idir * 
     @    6.d0 * (-sigmaf(i+m-1,js,k+l  ,5) + sigmaf(i+m-1,js,k+l-1,5)
     @            -sigmaf(i+m  ,js,k+l  ,5) + sigmaf(i+m  ,js,k+l-1,5))

          res(i,j,k) = res(i,j,k) - fac1 * (tmp - cen * fdst(i+m,j,k+l))

 40      continue

      else
         k = regl2
         if (idir .eq. 1) then
            ks = k - 1
         else
            ks = k
         endif
         do 50 j = regl1, regh1
            do 50 i = regl0, regh0

                  tmp = sigmac(i-1,j-1,ks,1) * MXXC(i-1,j-1,k-idir)
     @                 +sigmac(i-1,j  ,ks,1) * MXXC(i-1,j+1,k-idir)
     @                 +sigmac(i  ,j-1,ks,1) * MXXC(i+1,j-1,k-idir)
     @                 +sigmac(i  ,j  ,ks,1) * MXXC(i+1,j+1,k-idir)
                  tmp = tmp
     @                 +sigmac(i-1,j-1,ks,2) * MYYC(i-1,j-1,k-idir)
     @                 +sigmac(i-1,j  ,ks,2) * MYYC(i-1,j+1,k-idir)
     @                 +sigmac(i  ,j-1,ks,2) * MYYC(i+1,j-1,k-idir)
     @                 +sigmac(i  ,j  ,ks,2) * MYYC(i+1,j+1,k-idir)
                  tmp = tmp
     @                 +sigmac(i-1,j-1,ks,3) * MZZC(i-1,j-1,k-idir)
     @                 +sigmac(i-1,j  ,ks,3) * MZZC(i-1,j+1,k-idir)
     @                 +sigmac(i  ,j-1,ks,3) * MZZC(i+1,j-1,k-idir)
     @                 +sigmac(i  ,j  ,ks,3) * MZZC(i+1,j+1,k-idir)
                  tmp = tmp + idir * (
     @                 -sigmac(i-1,j-1,ks,4) * MXZC(i-1,j-1,k-idir)
     @                 -sigmac(i-1,j  ,ks,4) * MXZC(i-1,j+1,k-idir)
     @                 +sigmac(i  ,j-1,ks,4) * MXZC(i+1,j-1,k-idir)
     @                 +sigmac(i  ,j  ,ks,4) * MXZC(i+1,j+1,k-idir) )
                  tmp = tmp + idir * (
     @                 -sigmac(i-1,j-1,ks,5) * MYZC(i-1,j-1,k-idir)
     @                 +sigmac(i-1,j  ,ks,5) * MYZC(i-1,j+1,k-idir)
     @                 -sigmac(i  ,j-1,ks,5) * MYZC(i+1,j-1,k-idir)
     @                 +sigmac(i  ,j  ,ks,5) * MYZC(i+1,j+1,k-idir) )

               cen = 
     @            4.d0 * ( sigmac(i-1,j-1,ks,1)
     @                    +sigmac(i-1,j  ,ks,1) 
     @                    +sigmac(i  ,j-1,ks,1)
     @                    +sigmac(i  ,j  ,ks,1)
     @                    +sigmac(i-1,j-1,ks,2)
     @                    +sigmac(i-1,j  ,ks,2)
     @                    +sigmac(i  ,j-1,ks,2)
     @                    +sigmac(i  ,j  ,ks,2)
     @                    +sigmac(i-1,j-1,ks,3)
     @                    +sigmac(i-1,j  ,ks,3)
     @                    +sigmac(i  ,j-1,ks,3)
     @                    +sigmac(i  ,j  ,ks,3))
     @   + idir * 6.d0 * (-sigmac(i-1,j-1,ks,4)   
     @                    -sigmac(i-1,j  ,ks,4)
     @                    +sigmac(i  ,j-1,ks,4) 
     @                    +sigmac(i  ,j  ,ks,4) )
     @   + idir * 6.d0 * (-sigmac(i-1,j-1,ks,5)  
     @                    +sigmac(i-1,j,ks,5)
     @                    -sigmac(i  ,j-1,ks,5)     
     @                    +sigmac(i  ,j,ks,5))

               res(i*ir,j*jr,k*kr) =
     @            src(i*ir,j*jr,k*kr) - fac00 * (tmp - cen * cdst(i,j,k))
 50      continue

         fac0 = fac00 / (ir * jr)
         k = k * kr
         if (idir .eq. 1) then
            ks = k
         else
            ks = k - 1
         endif
         do 60 n = -(jr-1), jr-1
            fac2 = (jr-abs(n)) * fac0
            do 60 m = -(ir-1), ir-1
               fac1 = (ir-abs(m)) * fac2
               do 60 j = jr*regl1, jr*regh1, jr
                  do 60 i = ir*regl0, ir*regh0, ir

            tmp = sigmaf(i+m-1,j+n-1,ks,1) * MXXF(i+m,j+n,k,i+m-1,j+n-1,k+idir)
     @           +sigmaf(i+m-1,j+n  ,ks,1) * MXXF(i+m,j+n,k,i+m-1,j+n+1,k+idir)
     @           +sigmaf(i+m  ,j+n-1,ks,1) * MXXF(i+m,j+n,k,i+m+1,j+n-1,k+idir)
     @           +sigmaf(i+m  ,j+n  ,ks,1) * MXXF(i+m,j+n,k,i+m+1,j+n+1,k+idir)
            tmp = tmp
     @           +sigmaf(i+m-1,j+n-1,ks,2) * MYYF(i+m,j+n,k,i+m-1,j+n-1,k+idir)
     @           +sigmaf(i+m-1,j+n  ,ks,2) * MYYF(i+m,j+n,k,i+m-1,j+n+1,k+idir)
     @           +sigmaf(i+m  ,j+n-1,ks,2) * MYYF(i+m,j+n,k,i+m+1,j+n-1,k+idir)
     @           +sigmaf(i+m  ,j+n  ,ks,2) * MYYF(i+m,j+n,k,i+m+1,j+n+1,k+idir)
            tmp = tmp
     @           +sigmaf(i+m-1,j+n-1,ks,3) * MZZF(i+m,j+n,k,i+m-1,j+n-1,k+idir)
     @           +sigmaf(i+m-1,j+n  ,ks,3) * MZZF(i+m,j+n,k,i+m-1,j+n+1,k+idir)
     @           +sigmaf(i+m  ,j+n-1,ks,3) * MZZF(i+m,j+n,k,i+m+1,j+n-1,k+idir)
     @           +sigmaf(i+m  ,j+n  ,ks,3) * MZZF(i+m,j+n,k,i+m+1,j+n+1,k+idir)
            tmp = tmp + idir * (
     @            sigmaf(i+m-1,j+n-1,ks,4) * MXZF(i+m,j+n,k,i+m-1,j+n-1,k+idir)
     @           +sigmaf(i+m-1,j+n  ,ks,4) * MXZF(i+m,j+n,k,i+m-1,j+n+1,k+idir)
     @           -sigmaf(i+m  ,j+n-1,ks,4) * MXZF(i+m,j+n,k,i+m+1,j+n-1,k+idir)
     @           -sigmaf(i+m  ,j+n  ,ks,4) * MXZF(i+m,j+n,k,i+m+1,j+n+1,k+idir))
            tmp = tmp + idir * (
     @            sigmaf(i+m-1,j+n-1,ks,5) * MYZF(i+m,j+n,k,i+m-1,j+n-1,k+idir)
     @           -sigmaf(i+m-1,j+n  ,ks,5) * MYZF(i+m,j+n,k,i+m-1,j+n+1,k+idir)
     @           +sigmaf(i+m  ,j+n-1,ks,5) * MYZF(i+m,j+n,k,i+m+1,j+n-1,k+idir)
     @           -sigmaf(i+m  ,j+n  ,ks,5) * MYZF(i+m,j+n,k,i+m+1,j+n+1,k+idir))

         cen =
     @      4.d0 * (sigmaf(i+m-1,j+n-1,ks,1)
     @             +sigmaf(i+m-1,j+n  ,ks,1)
     @             +sigmaf(i+m  ,j+n-1,ks,1)
     @             +sigmaf(i+m  ,j+n  ,ks,1)
     @             +sigmaf(i+m-1,j+n-1,ks,2)
     @             +sigmaf(i+m-1,j+n  ,ks,2)
     @             +sigmaf(i+m  ,j+n-1,ks,2)
     @             +sigmaf(i+m  ,j+n  ,ks,2)
     @             +sigmaf(i+m-1,j+n-1,ks,3)
     @             +sigmaf(i+m-1,j+n  ,ks,3)
     @             +sigmaf(i+m  ,j+n-1,ks,3)
     @             +sigmaf(i+m  ,j+n  ,ks,3))
     @    + idir * 
     @      6.d0 * (sigmaf(i+m-1,j+n-1,ks,4)
     @             +sigmaf(i+m-1,j+n  ,ks,4)
     @             -sigmaf(i+m  ,j+n-1,ks,4)
     @             -sigmaf(i+m  ,j+n  ,ks,4))
     @    + idir * 
     @    + 6.d0 * (sigmaf(i+m-1,j+n-1,ks,5)
     @             -sigmaf(i+m-1,j+n  ,ks,5)
     @             +sigmaf(i+m  ,j+n-1,ks,5)
     @             -sigmaf(i+m  ,j+n  ,ks,5))

          res(i,j,k) = res(i,j,k) - fac1 * (tmp - cen * fdst(i+m,j+n,k))

 60      continue

      endif
      end

#undef MXXC
#undef MYYC
#undef MZZC
#undef MXZC
#undef MYZC

#define MXXC(ii,jj,kk) (4.d0 *  cdst(ii,jc,kc) + \
                        2.d0 * (cdst(ii,jj,kc) - cdst(ic,jj,kc) + \
                                cdst(ii,jc,kk) - cdst(ic,jc,kk)) + \
                               (cdst(ii,jj,kk) - cdst(ic,jj,kk)))

#define MYYC(ii,jj,kk) (4.d0 *  cdst(ic,jj,kc) + \
                        2.d0 * (cdst(ii,jj,kc) - cdst(ii,jc,kc) + \
                                cdst(ic,jj,kk) - cdst(ic,jc,kk)) + \
                               (cdst(ii,jj,kk) - cdst(ii,jc,kk)))

#define MZZC(ii,jj,kk) (4.d0 *  cdst(ic,jc,kk) + \
                        2.d0 * (cdst(ii,jc,kk) - cdst(ii,jc,kc) + \
                                cdst(ic,jj,kk) - cdst(ic,jj,kc)) + \
                               (cdst(ii,jj,kk) - cdst(ii,jj,kc)))

#define MXZC(ii,jj,kk) (6.d0 *  cdst(ii,jc,kk) + \
                        3.d0 * (cdst(ii,jj,kk) - cdst(ic,jj,kc)))

#define MYZC(ii,jj,kk) (6.d0 *  cdst(ic,jj,kk) + \
                        3.d0 * (cdst(ii,jj,kk) - cdst(ii,jc,kc)))

c Note---assumes fdst linearly interpolated from cdst along face
      subroutine hgeres(res, dims(res),
     @ src,  dims(src),
     @ fdst, dims(fdst),
     @ cdst, dims(cdst),
     @ sigmaf, dims(sf),
     @ sigmac, dims(sc),
     @ dims(reg),
     @ ir, jr, kr, ivect, ga)
      integer dims(res)
      integer dims(src)
      integer dims(fdst)
      integer dims(cdst)
      integer dims(sf)
      integer dims(sc)
      integer dims(reg)
      real*8 res(dimdec(res))
      real*8 src(dimdec(src))
      real*8 fdst(dimdec(fdst))
      real*8 cdst(dimdec(cdst))
      real*8 sigmaf(dimdec(sf),5)
      real*8 sigmac(dimdec(sc),5)
      integer ir, jr, kr, ivect(0:2), ga(0:1,0:1,0:1)
      real*8 fac00, fac0, fac1, fac, tmp, cen
      integer ic, jc, kc, if, jf, kf, ii, ji, ki, idir, jdir, kdir
      integer l, m, n

      ic = regl0
      jc = regl1
      kc = regl2
      if = ic * ir
      jf = jc * jr
      kf = kc * kr

      fac00 = 1.d0 / 36.d0

      if (ivect(0) .eq. 0) then
         do 10 if = ir*regl0, ir*regh0, ir
 10         res(if,jf,kf) = 0.d0
c quadrants
c each quadrant is two octants and their share of the two central edges
         fac0 = 1.d0 / ir
         do 40 ki = 0, 1
            kdir = 2 * ki - 1
            do 40 ji = 0, 1
               jdir = 2 * ji - 1
               if (ga(0,ji,ki) .eq. 1) then
                  do 20 m = -(ir-1), ir-1
                     fac = (ir-abs(m)) * fac0
                     do 20 if = ir*regl0, ir*regh0, ir

                  tmp = sigmaf(if+m-1,jf+ji-1,kf+ki-1,1) * 
     @                    MXXF(if+m,jf,kf,if+m-1,jf+jdir,kf+kdir)
     @                 +sigmaf(if+m  ,jf+ji-1,kf+ki-1,1) * 
     @                    MXXF(if+m,jf,kf,if+m+1,jf+jdir,kf+kdir)
                  tmp = tmp
     @                 +sigmaf(if+m-1,jf+ji-1,kf+ki-1,2) * 
     @                    MYYF(if+m,jf,kf,if+m-1,jf+jdir,kf+kdir)
     @                 +sigmaf(if+m  ,jf+ji-1,kf+ki-1,2) * 
     @                    MYYF(if+m,jf,kf,if+m+1,jf+jdir,kf+kdir)
                  tmp = tmp 
     @                 +sigmaf(if+m-1,jf+ji-1,kf+ki-1,3) * 
     @                    MZZF(if+m,jf,kf,if+m-1,jf+jdir,kf+kdir)
     @                 +sigmaf(if+m  ,jf+ji-1,kf+ki-1,3) * 
     @                    MZZF(if+m,jf,kf,if+m+1,jf+jdir,kf+kdir)
                  tmp = tmp + kdir * (
     @                  sigmaf(if+m-1,jf+ji-1,kf+ki-1,4) * 
     @                    MXZF(if+m,jf,kf,if+m-1,jf+jdir,kf+kdir)
     @                 -sigmaf(if+m  ,jf+ji-1,kf+ki-1,4) * 
     @                    MXZF(if+m,jf,kf,if+m+1,jf+jdir,kf+kdir) )
                  tmp = tmp + jdir * kdir * (
     @                 -sigmaf(if+m-1,jf+ji-1,kf+ki-1,5) * 
     @                    MYZF(if+m,jf,kf,if+m-1,jf+jdir,kf+kdir)
     @                 -sigmaf(if+m  ,jf+ji-1,kf+ki-1,5) * 
     @                    MYZF(if+m,jf,kf,if+m+1,jf+jdir,kf+kdir) )

                  cen =
     @               4.d0 * ( sigmaf(if+m-1,jf+ji-1,kf+ki-1,1)
     @                       +sigmaf(if+m  ,jf+ji-1,kf+ki-1,1)
     @                       +sigmaf(if+m-1,jf+ji-1,kf+ki-1,2)
     @                       +sigmaf(if+m  ,jf+ji-1,kf+ki-1,2)
     @                       +sigmaf(if+m-1,jf+ji-1,kf+ki-1,3)
     @                       +sigmaf(if+m  ,jf+ji-1,kf+ki-1,3))
     @             + 6.d0 * ( sigmaf(if+m-1,jf+ji-1,kf+ki-1,4)
     @                       -sigmaf(if+m  ,jf+ji-1,kf+ki-1,4)) * kdir 
     @             + 6.d0 * (-sigmaf(if+m-1,jf+ji-1,kf+ki-1,5)
     @                       -sigmaf(if+m  ,jf+ji-1,kf+ki-1,5)) * jdir * kdir 

                  tmp = tmp - cen * fdst(if+m,jf,kf)

                        res(if,jf,kf) = res(if,jf,kf) + fac * tmp
 20                  continue

               else
                  do 30 ic = regl0, regh0
                     if = ic * ir

                  tmp = sigmac(ic-1,jc+ji-1,kc+ki-1,1) * 
     @                    MXXC(ic-1,jc+jdir,kc+kdir)
     @                 +sigmac(ic  ,jc+ji-1,kc+ki-1,1) * 
     @                    MXXC(ic+1,jc+jdir,kc+kdir)
                  tmp = tmp
     @                 +sigmac(ic-1,jc+ji-1,kc+ki-1,2) * 
     @                    MYYC(ic-1,jc+jdir,kc+kdir)
     @                 +sigmac(ic  ,jc+ji-1,kc+ki-1,2) * 
     @                    MYYC(ic+1,jc+jdir,kc+kdir)
                  tmp = tmp
     @                 +sigmac(ic-1,jc+ji-1,kc+ki-1,3) * 
     @                    MZZC(ic-1,jc+jdir,kc+kdir)
     @                 +sigmac(ic  ,jc+ji-1,kc+ki-1,3) * 
     @                    MZZC(ic+1,jc+jdir,kc+kdir)
                  tmp = tmp + kdir * (
     @                  sigmac(ic-1,jc+ji-1,kc+ki-1,4) * 
     @                    MXZC(ic-1,jc+jdir,kc+kdir)
     @                 -sigmac(ic  ,jc+ji-1,kc+ki-1,4) * 
     @                    MXZC(ic+1,jc+jdir,kc+kdir) )
                  tmp = tmp + jdir * kdir * ( 
     @                 -sigmac(ic-1,jc+ji-1,kc+ki-1,5) * 
     @                    MYZC(ic-1,jc+jdir,kc+kdir)
     @                 -sigmac(ic  ,jc+ji-1,kc+ki-1,5) * 
     @                    MYZC(ic+1,jc+jdir,kc+kdir) )
                   cen =
     @               4.d0 * ( sigmac(ic-1,jc+ji-1,kc+ki-1,1)
     @                       +sigmac(ic  ,jc+ji-1,kc+ki-1,1)
     @                       +sigmac(ic-1,jc+ji-1,kc+ki-1,2)
     @                       +sigmac(ic  ,jc+ji-1,kc+ki-1,2)
     @                       +sigmac(ic-1,jc+ji-1,kc+ki-1,3)
     @                       +sigmac(ic  ,jc+ji-1,kc+ki-1,3))
     @             + 6.d0 * ( sigmac(ic-1,jc+ji-1,kc+ki-1,4)
     @                       -sigmac(ic  ,jc+ji-1,kc+ki-1,4)) * kdir
     @             + 6.d0 * (-sigmac(ic-1,jc+ji-1,kc+ki-1,5)
     @                       -sigmac(ic  ,jc+ji-1,kc+ki-1,5)) * jdir * kdir

                     tmp = tmp - cen * cdst(ic,jc,kc)
                     res(if,jf,kf) = res(if,jf,kf) + tmp
 30                  continue
               endif
 40            continue
c faces
c each face is two faces and two sides of an edge
         do 70 ki = 0, 1
            kdir = 2 * ki - 1
            do 70 ji = 0, 1
               jdir = 2 * ji - 1

               if (ga(0,ji,ki) - ga(0,ji,1-ki) .eq. 1) then
                  fac0 = 1.d0 / (ir * jr)
                  do 50 n = jdir, jdir*(jr-1), jdir
                     fac1 = (jr-abs(n)) * fac0
                     do 50 m = -(ir-1), ir-1
                        fac = (ir-abs(m)) * fac1
                        do 50 if = ir*regl0, ir*regh0, ir

                  tmp = sigmaf(if+m-1,jf+n-1,kf+ki-1,1) * 
     @                    MXXF(if+m,jf+n,kf,if+m-1,jf+n-1,kf+kdir)
     @                 +sigmaf(if+m-1,jf+n  ,kf+ki-1,1) * 
     @                    MXXF(if+m,jf+n,kf,if+m-1,jf+n+1,kf+kdir)
     @                 +sigmaf(if+m  ,jf+n-1,kf+ki-1,1) * 
     @                    MXXF(if+m,jf+n,kf,if+m+1,jf+n-1,kf+kdir)
     @                 +sigmaf(if+m  ,jf+n  ,kf+ki-1,1) * 
     @                    MXXF(if+m,jf+n,kf,if+m+1,jf+n+1,kf+kdir)
                  tmp = tmp
     @                 +sigmaf(if+m-1,jf+n-1,kf+ki-1,2) * 
     @                    MYYF(if+m,jf+n,kf,if+m-1,jf+n-1,kf+kdir)
     @                 +sigmaf(if+m-1,jf+n  ,kf+ki-1,2) * 
     @                    MYYF(if+m,jf+n,kf,if+m-1,jf+n+1,kf+kdir)
     @                 +sigmaf(if+m  ,jf+n-1,kf+ki-1,2) * 
     @                    MYYF(if+m,jf+n,kf,if+m+1,jf+n-1,kf+kdir)
     @                 +sigmaf(if+m  ,jf+n  ,kf+ki-1,2) * 
     @                    MYYF(if+m,jf+n,kf,if+m+1,jf+n+1,kf+kdir)
                  tmp = tmp
     @                 +sigmaf(if+m-1,jf+n-1,kf+ki-1,3) * 
     @                    MZZF(if+m,jf+n,kf,if+m-1,jf+n-1,kf+kdir)
     @                 +sigmaf(if+m-1,jf+n  ,kf+ki-1,3) * 
     @                    MZZF(if+m,jf+n,kf,if+m-1,jf+n+1,kf+kdir)
     @                 +sigmaf(if+m  ,jf+n-1,kf+ki-1,3) * 
     @                    MZZF(if+m,jf+n,kf,if+m+1,jf+n-1,kf+kdir)
     @                 +sigmaf(if+m  ,jf+n  ,kf+ki-1,3) * 
     @                    MZZF(if+m,jf+n,kf,if+m+1,jf+n+1,kf+kdir)
                  tmp = tmp + kdir * (
     @                 +sigmaf(if+m-1,jf+n-1,kf+ki-1,4) * 
     @                    MXZF(if+m,jf+n,kf,if+m-1,jf+n-1,kf+kdir)
     @                 +sigmaf(if+m-1,jf+n  ,kf+ki-1,4) * 
     @                    MXZF(if+m,jf+n,kf,if+m-1,jf+n+1,kf+kdir)
     @                 -sigmaf(if+m  ,jf+n-1,kf+ki-1,4) * 
     @                    MXZF(if+m,jf+n,kf,if+m+1,jf+n-1,kf+kdir)
     @                 -sigmaf(if+m  ,jf+n  ,kf+ki-1,4) * 
     @                    MXZF(if+m,jf+n,kf,if+m+1,jf+n+1,kf+kdir) )
                  tmp = tmp + kdir * ( 
     @                 +sigmaf(if+m-1,jf+n-1,kf+ki-1,5) * 
     @                    MYZF(if+m,jf+n,kf,if+m-1,jf+n-1,kf+kdir)
     @                 -sigmaf(if+m-1,jf+n  ,kf+ki-1,5) * 
     @                    MYZF(if+m,jf+n,kf,if+m-1,jf+n+1,kf+kdir)
     @                 +sigmaf(if+m  ,jf+n-1,kf+ki-1,5) * 
     @                    MYZF(if+m,jf+n,kf,if+m+1,jf+n-1,kf+kdir)
     @                 -sigmaf(if+m  ,jf+n  ,kf+ki-1,5) * 
     @                    MYZF(if+m,jf+n,kf,if+m+1,jf+n+1,kf+kdir) )

         cen =
     @      4.d0 * (sigmaf(if+m-1,jf+n-1,kf+ki-1,1)
     @             +sigmaf(if+m-1,jf+n  ,kf+ki-1,1)
     @             +sigmaf(if+m  ,jf+n-1,kf+ki-1,1)
     @             +sigmaf(if+m  ,jf+n  ,kf+ki-1,1)
     @             +sigmaf(if+m-1,jf+n-1,kf+ki-1,2)
     @             +sigmaf(if+m-1,jf+n  ,kf+ki-1,2)
     @             +sigmaf(if+m  ,jf+n-1,kf+ki-1,2)
     @             +sigmaf(if+m  ,jf+n  ,kf+ki-1,2)
     @             +sigmaf(if+m-1,jf+n-1,kf+ki-1,3)
     @             +sigmaf(if+m-1,jf+n  ,kf+ki-1,3)
     @             +sigmaf(if+m  ,jf+n-1,kf+ki-1,3)
     @             +sigmaf(if+m  ,jf+n  ,kf+ki-1,3))
     @    + 6.d0 * (sigmaf(if+m-1,jf+n-1,kf+ki-1,4)
     @             +sigmaf(if+m-1,jf+n  ,kf+ki-1,4)
     @             -sigmaf(if+m  ,jf+n-1,kf+ki-1,4)
     @             -sigmaf(if+m  ,jf+n  ,kf+ki-1,4)) * kdir
     @    + 6.d0 * (sigmaf(if+m-1,jf+n-1,kf+ki-1,5)
     @             -sigmaf(if+m-1,jf+n  ,kf+ki-1,5)
     @             +sigmaf(if+m  ,jf+n-1,kf+ki-1,5)
     @             -sigmaf(if+m  ,jf+n  ,kf+ki-1,5)) * kdir

                  tmp = tmp - cen * fdst(if+m,jf+n,kf)
                  res(if,jf,kf) = res(if,jf,kf) + fac * tmp

 50            continue
               endif

               if (ga(0,ji,ki) - ga(0,1-ji,ki) .eq. 1) then
                  fac0 = 1.d0 / (ir * kr)
                  do 60 l = kdir, kdir*(kr-1), kdir
                     fac1 = (kr-abs(l)) * fac0
                     do 60 m = -(ir-1), ir-1
                        fac = (ir-abs(m)) * fac1
                        do 60 if = ir*regl0, ir*regh0, ir

                  tmp = sigmaf(if+m-1,jf+ji-1,kf+l-1,1) * 
     @                    MXXF(if+m,jf,kf+l,if+m-1,jf+jdir,kf+l-1)
     @                 +sigmaf(if+m-1,jf+ji-1,kf+l  ,1) * 
     @                    MXXF(if+m,jf,kf+l,if+m-1,jf+jdir,kf+l+1)
     @                 +sigmaf(if+m  ,jf+ji-1,kf+l-1,1) * 
     @                    MXXF(if+m,jf,kf+l,if+m+1,jf+jdir,kf+l-1)
     @                 +sigmaf(if+m  ,jf+ji-1,kf+l  ,1) * 
     @                    MXXF(if+m,jf,kf+l,if+m+1,jf+jdir,kf+l+1)
                  tmp = tmp
     @                 +sigmaf(if+m-1,jf+ji-1,kf+l-1,2) * 
     @                    MYYF(if+m,jf,kf+l,if+m-1,jf+jdir,kf+l-1)
     @                 +sigmaf(if+m-1,jf+ji-1,kf+l  ,2) * 
     @                    MYYF(if+m,jf,kf+l,if+m-1,jf+jdir,kf+l+1)
     @                 +sigmaf(if+m  ,jf+ji-1,kf+l-1,2) * 
     @                    MYYF(if+m,jf,kf+l,if+m+1,jf+jdir,kf+l-1)
     @                 +sigmaf(if+m  ,jf+ji-1,kf+l  ,2) * 
     @                    MYYF(if+m,jf,kf+l,if+m+1,jf+jdir,kf+l+1)
                  tmp = tmp
     @                 +sigmaf(if+m-1,jf+ji-1,kf+l-1,3) * 
     @                    MZZF(if+m,jf,kf+l,if+m-1,jf+jdir,kf+l-1)
     @                 +sigmaf(if+m-1,jf+ji-1,kf+l  ,3) * 
     @                    MZZF(if+m,jf,kf+l,if+m-1,jf+jdir,kf+l+1)
     @                 +sigmaf(if+m  ,jf+ji-1,kf+l-1,3) * 
     @                    MZZF(if+m,jf,kf+l,if+m+1,jf+jdir,kf+l-1)
     @                 +sigmaf(if+m  ,jf+ji-1,kf+l  ,3) * 
     @                    MZZF(if+m,jf,kf+l,if+m+1,jf+jdir,kf+l+1)
                  tmp = tmp + ( 
     @                 -sigmaf(if+m-1,jf+ji-1,kf+l-1,4) * 
     @                    MXZF(if+m,jf,kf+l,if+m-1,jf+jdir,kf+l-1)
     @                 +sigmaf(if+m-1,jf+ji-1,kf+l  ,4) * 
     @                    MXZF(if+m,jf,kf+l,if+m-1,jf+jdir,kf+l+1)
     @                 +sigmaf(if+m  ,jf+ji-1,kf+l-1,4) * 
     @                    MXZF(if+m,jf,kf+l,if+m+1,jf+jdir,kf+l-1)
     @                 -sigmaf(if+m  ,jf+ji-1,kf+l  ,4) * 
     @                    MXZF(if+m,jf,kf+l,if+m+1,jf+jdir,kf+l+1) )
                  tmp = tmp + jdir * (
     @                 +sigmaf(if+m-1,jf+ji-1,kf+l-1,5) * 
     @                    MYZF(if+m,jf,kf+l,if+m-1,jf+jdir,kf+l-1)
     @                 -sigmaf(if+m-1,jf+ji-1,kf+l  ,5) * 
     @                    MYZF(if+m,jf,kf+l,if+m-1,jf+jdir,kf+l+1)
     @                 +sigmaf(if+m  ,jf+ji-1,kf+l-1,5) * 
     @                    MYZF(if+m,jf,kf+l,if+m+1,jf+jdir,kf+l-1)
     @                 -sigmaf(if+m  ,jf+ji-1,kf+l  ,5) * 
     @                    MYZF(if+m,jf,kf+l,if+m+1,jf+jdir,kf+l+1) )
         cen =
     @      4.d0 * (
     @             +sigmaf(if+m-1,jf+ji-1,kf+l-1,1)  
     @             +sigmaf(if+m-1,jf+ji-1,kf+l  ,1)
     @             +sigmaf(if+m  ,jf+ji-1,kf+l-1,1)  
     @             +sigmaf(if+m  ,jf+ji-1,kf+l  ,1)
     @             +sigmaf(if+m-1,jf+ji-1,kf+l-1,2)  
     @             +sigmaf(if+m-1,jf+ji-1,kf+l  ,2)
     @             +sigmaf(if+m  ,jf+ji-1,kf+l-1,2)  
     @             +sigmaf(if+m  ,jf+ji-1,kf+l  ,2)
     @             +sigmaf(if+m-1,jf+ji-1,kf+l-1,3)  
     @             +sigmaf(if+m-1,jf+ji-1,kf+l  ,3)
     @             +sigmaf(if+m  ,jf+ji-1,kf+l-1,3)  
     @             +sigmaf(if+m  ,jf+ji-1,kf+l  ,3))
     @    + 6.d0 * (sigmaf(if+m-1,jf+ji-1,kf+l  ,4)  
     @             -sigmaf(if+m-1,jf+ji-1,kf+l-1,4)
     @             +sigmaf(if+m  ,jf+ji-1,kf+l-1,4) 
     @             -sigmaf(if+m  ,jf+ji-1,kf+l  ,4))
     @    + 6.d0 * (sigmaf(if+m-1,jf+ji-1,kf+l-1,5)  
     @             -sigmaf(if+m-1,jf+ji-1,kf+l  ,5)
     @             +sigmaf(if+m  ,jf+ji-1,kf+l-1,5)  
     @             -sigmaf(if+m  ,jf+ji-1,kf+l  ,5)) * jdir

                  tmp = tmp - cen * fdst(if+m,jf,kf+l)
                  res(if,jf,kf) = res(if,jf,kf) + fac * tmp

 60            continue

               endif
 70            continue
c weighting
         do 80 if = ir*regl0, ir*regh0, ir
            res(if,jf,kf) = src(if,jf,kf) - res(if,jf,kf) * fac00
 80      continue

      else if (ivect(1) .eq. 0) then
         do 110 jf = jr*regl1, jr*regh1, jr
 110        res(if,jf,kf) = 0.d0
c quadrants
c each quadrant is two octants and their share of the two central edges
         fac0 = 1.d0 / jr
         do 140 ki = 0, 1
            kdir = 2 * ki - 1
            do 140 ii = 0, 1
               idir = 2 * ii - 1
               if (ga(ii,0,ki) .eq. 1) then
                  do 120 n = -(jr-1), jr-1
                     fac = (jr-abs(n)) * fac0
                     do 120 jf = jr*regl1, jr*regh1, jr

                  tmp = sigmaf(if+ii-1,jf+n-1,kf+ki-1,1) * 
     @                    MXXF(if,jf+n,kf,if+idir,jf+n-1,kf+kdir)
     @                 +sigmaf(if+ii-1,jf+n  ,kf+ki-1,1) * 
     @                    MXXF(if,jf+n,kf,if+idir,jf+n+1,kf+kdir)
                  tmp = tmp
     @                 +sigmaf(if+ii-1,jf+n-1,kf+ki-1,2) * 
     @                    MYYF(if,jf+n,kf,if+idir,jf+n-1,kf+kdir)
     @                 +sigmaf(if+ii-1,jf+n  ,kf+ki-1,2) * 
     @                    MYYF(if,jf+n,kf,if+idir,jf+n+1,kf+kdir) 
                  tmp = tmp 
     @                 +sigmaf(if+ii-1,jf+n-1,kf+ki-1,3) * 
     @                    MZZF(if,jf+n,kf,if+idir,jf+n-1,kf+kdir)
     @                 +sigmaf(if+ii-1,jf+n  ,kf+ki-1,3) * 
     @                    MZZF(if,jf+n,kf,if+idir,jf+n+1,kf+kdir)
                  tmp = tmp + idir * kdir * (
     @                 -sigmaf(if+ii-1,jf+n-1,kf+ki-1,4) * 
     @                    MXZF(if,jf+n,kf,if+idir,jf+n-1,kf+kdir)
     @                 -sigmaf(if+ii-1,jf+n  ,kf+ki-1,4) * 
     @                    MXZF(if,jf+n,kf,if+idir,jf+n+1,kf+kdir) )
                  tmp = tmp + kdir * (
     @                 +sigmaf(if+ii-1,jf+n-1,kf+ki-1,5) * 
     @                    MYZF(if,jf+n,kf,if+idir,jf+n-1,kf+kdir)
     @                 -sigmaf(if+ii-1,jf+n  ,kf+ki-1,5) * 
     @                    MYZF(if,jf+n,kf,if+idir,jf+n+1,kf+kdir) )
                  cen =
     @                  4.d0 * ( sigmaf(if+ii-1,jf+n-1,kf+ki-1,1)
     @                          +sigmaf(if+ii-1,jf+n  ,kf+ki-1,1)
     @                          +sigmaf(if+ii-1,jf+n-1,kf+ki-1,2)
     @                          +sigmaf(if+ii-1,jf+n  ,kf+ki-1,2)
     @                          +sigmaf(if+ii-1,jf+n-1,kf+ki-1,3)
     @                          +sigmaf(if+ii-1,jf+n  ,kf+ki-1,3))
     @                + 6.d0 * (-sigmaf(if+ii-1,jf+n-1,kf+ki-1,4)
     @                          -sigmaf(if+ii-1,jf+n  ,kf+ki-1,4)) * idir * kdir
     @                + 6.d0 * ( sigmaf(if+ii-1,jf+n-1,kf+ki-1,5) 
     @                          -sigmaf(if+ii-1,jf+n  ,kf+ki-1,5)) * kdir

                  tmp = tmp - cen * fdst(if,jf+n,kf)

                  res(if,jf,kf) = res(if,jf,kf) + fac * tmp
 120              continue

               else
                  do 130 jc = regl1, regh1
                     jf = jc * jr

                  tmp = sigmac(ic+ii-1,jc-1,kc+ki-1,1) * 
     @                    MXXC(ic+idir,jc-1,kc+kdir)
     @                 +sigmac(ic+ii-1,jc  ,kc+ki-1,1) * 
     @                    MXXC(ic+idir,jc+1,kc+kdir)
                  tmp = tmp
     @                 +sigmac(ic+ii-1,jc-1,kc+ki-1,2) * 
     @                    MYYC(ic+idir,jc-1,kc+kdir)
     @                 +sigmac(ic+ii-1,jc  ,kc+ki-1,2) * 
     @                    MYYC(ic+idir,jc+1,kc+kdir)
                  tmp = tmp
     @                 +sigmac(ic+ii-1,jc-1,kc+ki-1,3) * 
     @                    MZZC(ic+idir,jc-1,kc+kdir)
     @                 +sigmac(ic+ii-1,jc  ,kc+ki-1,3) * 
     @                    MZZC(ic+idir,jc+1,kc+kdir)
                  tmp = tmp + idir * kdir * (
     @                 -sigmac(ic+ii-1,jc-1,kc+ki-1,4) * 
     @                    MXZC(ic+idir,jc-1,kc+kdir)
     @                 -sigmac(ic+ii-1,jc  ,kc+ki-1,4) * 
     @                    MXZC(ic+idir,jc+1,kc+kdir) )
                  tmp = tmp + kdir * (
     @                 +sigmac(ic+ii-1,jc-1,kc+ki-1,5) * 
     @                    MYZC(ic+idir,jc-1,kc+kdir)
     @                 -sigmac(ic+ii-1,jc  ,kc+ki-1,5) * 
     @                    MYZC(ic+idir,jc+1,kc+kdir) )
                  cen =
     @               4.d0 * (
     @                       +sigmac(ic+ii-1,jc-1,kc+ki-1,1)
     @                       +sigmac(ic+ii-1,jc  ,kc+ki-1,1)
     @                       +sigmac(ic+ii-1,jc-1,kc+ki-1,2)
     @                       +sigmac(ic+ii-1,jc  ,kc+ki-1,2)
     @                       +sigmac(ic+ii-1,jc-1,kc+ki-1,3)
     @                       +sigmac(ic+ii-1,jc  ,kc+ki-1,3))
     @             + 6.d0 * (-sigmac(ic+ii-1,jc-1,kc+ki-1,4)
     @                       -sigmac(ic+ii-1,jc  ,kc+ki-1,4)) * idir * kdir
     @             + 6.d0 * ( sigmac(ic+ii-1,jc-1,kc+ki-1,5) 
     @                       -sigmac(ic+ii-1,jc  ,kc+ki-1,5)) * kdir

                     tmp = tmp - cen * cdst(ic,jc,kc)
                     res(if,jf,kf) = res(if,jf,kf) + tmp
 130                 continue

               endif
 140           continue
c faces
c each face is two faces and two sides of an edge
         do 170 ki = 0, 1
            kdir = 2 * ki - 1
            do 170 ii = 0, 1
               idir = 2 * ii - 1
               if (ga(ii,0,ki) - ga(ii,0,1-ki) .eq. 1) then
                  fac0 = 1.d0 / (ir * jr)
                  do 150 n = -(jr-1), jr-1
                     fac1 = (jr-abs(n)) * fac0
                     do 150 m = idir, idir*(ir-1), idir
                        fac = (ir-abs(m)) * fac1
                        do 150 jf = jr*regl1, jr*regh1, jr

                  tmp = sigmaf(if+m-1,jf+n-1,kf+ki-1,1) * 
     @                    MXXF(if+m,jf+n,kf,if+m-1,jf+n-1,kf+kdir)
     @                 +sigmaf(if+m-1,jf+n  ,kf+ki-1,1) * 
     @                    MXXF(if+m,jf+n,kf,if+m-1,jf+n+1,kf+kdir)
     @                 +sigmaf(if+m  ,jf+n-1,kf+ki-1,1) * 
     @                    MXXF(if+m,jf+n,kf,if+m+1,jf+n-1,kf+kdir)
     @                 +sigmaf(if+m  ,jf+n  ,kf+ki-1,1) * 
     @                    MXXF(if+m,jf+n,kf,if+m+1,jf+n+1,kf+kdir)
                  tmp = tmp
     @                 +sigmaf(if+m-1,jf+n-1,kf+ki-1,2) * 
     @                    MYYF(if+m,jf+n,kf,if+m-1,jf+n-1,kf+kdir)
     @                 +sigmaf(if+m-1,jf+n  ,kf+ki-1,2) * 
     @                    MYYF(if+m,jf+n,kf,if+m-1,jf+n+1,kf+kdir)
     @                 +sigmaf(if+m  ,jf+n-1,kf+ki-1,2) * 
     @                    MYYF(if+m,jf+n,kf,if+m+1,jf+n-1,kf+kdir)
     @                 +sigmaf(if+m  ,jf+n  ,kf+ki-1,2) * 
     @                    MYYF(if+m,jf+n,kf,if+m+1,jf+n+1,kf+kdir)
                  tmp = tmp
     @                 +sigmaf(if+m-1,jf+n-1,kf+ki-1,3) * 
     @                    MZZF(if+m,jf+n,kf,if+m-1,jf+n-1,kf+kdir)
     @                 +sigmaf(if+m-1,jf+n  ,kf+ki-1,3) * 
     @                    MZZF(if+m,jf+n,kf,if+m-1,jf+n+1,kf+kdir)
     @                 +sigmaf(if+m  ,jf+n-1,kf+ki-1,3) * 
     @                    MZZF(if+m,jf+n,kf,if+m+1,jf+n-1,kf+kdir)
     @                 +sigmaf(if+m  ,jf+n  ,kf+ki-1,3) * 
     @                    MZZF(if+m,jf+n,kf,if+m+1,jf+n+1,kf+kdir)
                  tmp = tmp + kdir * (
     @                  sigmaf(if+m-1,jf+n-1,kf+ki-1,4) * 
     @                    MXZF(if+m,jf+n,kf,if+m-1,jf+n-1,kf+kdir)
     @                 +sigmaf(if+m-1,jf+n  ,kf+ki-1,4) * 
     @                    MXZF(if+m,jf+n,kf,if+m-1,jf+n+1,kf+kdir)
     @                 -sigmaf(if+m  ,jf+n-1,kf+ki-1,4) * 
     @                    MXZF(if+m,jf+n,kf,if+m+1,jf+n-1,kf+kdir)
     @                 -sigmaf(if+m  ,jf+n  ,kf+ki-1,4) * 
     @                    MXZF(if+m,jf+n,kf,if+m+1,jf+n+1,kf+kdir) )
                  tmp = tmp + kdir * (
     @                  sigmaf(if+m-1,jf+n-1,kf+ki-1,5) * 
     @                    MYZF(if+m,jf+n,kf,if+m-1,jf+n-1,kf+kdir)
     @                 -sigmaf(if+m-1,jf+n  ,kf+ki-1,5) * 
     @                    MYZF(if+m,jf+n,kf,if+m-1,jf+n+1,kf+kdir)
     @                 +sigmaf(if+m  ,jf+n-1,kf+ki-1,5) * 
     @                    MYZF(if+m,jf+n,kf,if+m+1,jf+n-1,kf+kdir)
     @                 -sigmaf(if+m  ,jf+n  ,kf+ki-1,5) * 
     @                    MYZF(if+m,jf+n,kf,if+m+1,jf+n+1,kf+kdir) )

         cen =
     @      4.d0 * (sigmaf(if+m-1,jf+n-1,kf+ki-1,1)
     @             +sigmaf(if+m-1,jf+n  ,kf+ki-1,1)
     @             +sigmaf(if+m  ,jf+n-1,kf+ki-1,1)
     @             +sigmaf(if+m  ,jf+n  ,kf+ki-1,1)
     @             +sigmaf(if+m-1,jf+n-1,kf+ki-1,2)
     @             +sigmaf(if+m-1,jf+n  ,kf+ki-1,2)
     @             +sigmaf(if+m  ,jf+n-1,kf+ki-1,2)
     @             +sigmaf(if+m  ,jf+n  ,kf+ki-1,2)
     @             +sigmaf(if+m-1,jf+n-1,kf+ki-1,3)
     @             +sigmaf(if+m-1,jf+n  ,kf+ki-1,3)
     @             +sigmaf(if+m  ,jf+n-1,kf+ki-1,3)
     @             +sigmaf(if+m  ,jf+n  ,kf+ki-1,3))
     @    + 6.d0 * (sigmaf(if+m-1,jf+n-1,kf+ki-1,4)
     @             +sigmaf(if+m-1,jf+n  ,kf+ki-1,4)
     @             -sigmaf(if+m  ,jf+n-1,kf+ki-1,4)
     @             -sigmaf(if+m  ,jf+n  ,kf+ki-1,4)) * kdir
     @    + 6.d0 * (sigmaf(if+m-1,jf+n-1,kf+ki-1,5)
     @             -sigmaf(if+m-1,jf+n  ,kf+ki-1,5)
     @             +sigmaf(if+m  ,jf+n-1,kf+ki-1,5)
     @             -sigmaf(if+m  ,jf+n  ,kf+ki-1,5)) * kdir

                  tmp = tmp - cen * fdst(if+m,jf+n,kf)
                  res(if,jf,kf) = res(if,jf,kf) + fac * tmp

 150           continue

               endif
               if (ga(ii,0,ki) - ga(1-ii,0,ki) .eq. 1) then
                  fac0 = 1.d0 / (jr * kr)
                  do 160 l = kdir, kdir*(kr-1), kdir
                     fac1 = (kr-abs(l)) * fac0
                     do 160 n = -(jr-1), jr-1
                        fac = (jr-abs(n)) * fac1
                        do 160 jf = jr*regl1, jr*regh1, jr

                  tmp = sigmaf(if+ii-1,jf+n-1,kf+l-1,1) * 
     @                    MXXF(if,jf+n,kf+l,if+idir,jf+n-1,kf+l-1)
     @                 +sigmaf(if+ii-1,jf+n-1,kf+l  ,1) * 
     @                    MXXF(if,jf+n,kf+l,if+idir,jf+n-1,kf+l+1)
     @                 +sigmaf(if+ii-1,jf+n  ,kf+l-1,1) * 
     @                    MXXF(if,jf+n,kf+l,if+idir,jf+n+1,kf+l-1)
     @                 +sigmaf(if+ii-1,jf+n  ,kf+l  ,1) * 
     @                    MXXF(if,jf+n,kf+l,if+idir,jf+n+1,kf+l+1)
                  tmp = tmp
     @                 +sigmaf(if+ii-1,jf+n-1,kf+l-1,2) * 
     @                    MYYF(if,jf+n,kf+l,if+idir,jf+n-1,kf+l-1)
     @                 +sigmaf(if+ii-1,jf+n-1,kf+l  ,2) * 
     @                    MYYF(if,jf+n,kf+l,if+idir,jf+n-1,kf+l+1)
     @                 +sigmaf(if+ii-1,jf+n  ,kf+l-1,2) * 
     @                    MYYF(if,jf+n,kf+l,if+idir,jf+n+1,kf+l-1)
     @                 +sigmaf(if+ii-1,jf+n  ,kf+l  ,2) * 
     @                    MYYF(if,jf+n,kf+l,if+idir,jf+n+1,kf+l+1)
                  tmp = tmp
     @                 +sigmaf(if+ii-1,jf+n-1,kf+l-1,3) * 
     @                    MZZF(if,jf+n,kf+l,if+idir,jf+n-1,kf+l-1)
     @                 +sigmaf(if+ii-1,jf+n-1,kf+l  ,3) * 
     @                    MZZF(if,jf+n,kf+l,if+idir,jf+n-1,kf+l+1)
     @                 +sigmaf(if+ii-1,jf+n  ,kf+l-1,3) * 
     @                    MZZF(if,jf+n,kf+l,if+idir,jf+n+1,kf+l-1)
     @                 +sigmaf(if+ii-1,jf+n  ,kf+l  ,3) * 
     @                    MZZF(if,jf+n,kf+l,if+idir,jf+n+1,kf+l+1)
                  tmp = tmp + idir * (
     @                 +sigmaf(if+ii-1,jf+n-1,kf+l-1,4) * 
     @                    MXZF(if,jf+n,kf+l,if+idir,jf+n-1,kf+l-1)
     @                 -sigmaf(if+ii-1,jf+n-1,kf+l  ,4) * 
     @                    MXZF(if,jf+n,kf+l,if+idir,jf+n-1,kf+l+1)
     @                 +sigmaf(if+ii-1,jf+n  ,kf+l-1,4) * 
     @                    MXZF(if,jf+n,kf+l,if+idir,jf+n+1,kf+l-1)
     @                 -sigmaf(if+ii-1,jf+n  ,kf+l  ,4) * 
     @                    MXZF(if,jf+n,kf+l,if+idir,jf+n+1,kf+l+1) )
                  tmp = tmp + (
     @                 -sigmaf(if+ii-1,jf+n-1,kf+l-1,5) * 
     @                    MYZF(if,jf+n,kf+l,if+idir,jf+n-1,kf+l-1)
     @                 +sigmaf(if+ii-1,jf+n-1,kf+l  ,5) * 
     @                    MYZF(if,jf+n,kf+l,if+idir,jf+n-1,kf+l+1)
     @                 +sigmaf(if+ii-1,jf+n  ,kf+l-1,5) * 
     @                    MYZF(if,jf+n,kf+l,if+idir,jf+n+1,kf+l-1)
     @                 -sigmaf(if+ii-1,jf+n  ,kf+l  ,5) * 
     @                    MYZF(if,jf+n,kf+l,if+idir,jf+n+1,kf+l+1) )

                  cen =
     @               4.d0 * (sigmaf(if+ii-1,jf+n-1,kf+l-1,1) 
     @                      +sigmaf(if+ii-1,jf+n-1,kf+l  ,1)
     @                      +sigmaf(if+ii-1,jf+n  ,kf+l-1,1) 
     @                      +sigmaf(if+ii-1,jf+n  ,kf+l  ,1)
     @                      +sigmaf(if+ii-1,jf+n-1,kf+l-1,2) 
     @                      +sigmaf(if+ii-1,jf+n-1,kf+l  ,2)
     @                      +sigmaf(if+ii-1,jf+n  ,kf+l-1,2) 
     @                      +sigmaf(if+ii-1,jf+n  ,kf+l  ,2)
     @                      +sigmaf(if+ii-1,jf+n-1,kf+l-1,3) 
     @                      +sigmaf(if+ii-1,jf+n-1,kf+l  ,3)
     @                      +sigmaf(if+ii-1,jf+n  ,kf+l-1,3) 
     @                      +sigmaf(if+ii-1,jf+n  ,kf+l  ,3))
     @             + 6.d0 * (sigmaf(if+ii-1,jf+n-1,kf+l-1,4) 
     @                      -sigmaf(if+ii-1,jf+n-1,kf+l  ,4)
     @                      +sigmaf(if+ii-1,jf+n  ,kf+l-1,4) 
     @                      -sigmaf(if+ii-1,jf+n  ,kf+l  ,4)) * idir
     @             + 6.d0 * (sigmaf(if+ii-1,jf+n-1,kf+l  ,5) 
     @                      -sigmaf(if+ii-1,jf+n-1,kf+l-1,5)
     @                      +sigmaf(if+ii-1,jf+n  ,kf+l-1,5) 
     @                      -sigmaf(if+ii-1,jf+n  ,kf+l  ,5))

                  tmp = tmp - cen * fdst(if,jf+n,kf+l)
                  res(if,jf,kf) = res(if,jf,kf) + fac * tmp

 160           continue
               endif
 170           continue
c weighting
         do 180 jf = jr*regl1, jr*regh1, jr
            res(if,jf,kf) = src(if,jf,kf) - res(if,jf,kf) * fac00
 180     continue

      else
         do 210 kf = kr*regl2, kr*regh2, kr
 210        res(if,jf,kf) = 0.d0
c quadrants
c each quadrant is two octants and their share of the two central edges
         fac0 = 1.d0 / kr
         do 240 ji = 0, 1
            jdir = 2 * ji - 1
            do 240 ii = 0, 1
               idir = 2 * ii - 1
               if (ga(ii,ji,0) .eq. 1) then
                  do 220 l = -(kr-1), kr-1
                     fac = (kr-abs(l)) * fac0
                     do 220 kf = kr*regl2, kr*regh2, kr

                  tmp = sigmaf(if+ii-1,jf+ji-1,kf+l-1,1) * 
     @                    MXXF(if,jf,kf+l,if+idir,jf+jdir,kf+l-1)
     @                 +sigmaf(if+ii-1,jf+ji-1,kf+l  ,1) * 
     @                    MXXF(if,jf,kf+l,if+idir,jf+jdir,kf+l+1)
                  tmp = tmp
     @                 +sigmaf(if+ii-1,jf+ji-1,kf+l-1,2) * 
     @                    MYYF(if,jf,kf+l,if+idir,jf+jdir,kf+l-1)
     @                 +sigmaf(if+ii-1,jf+ji-1,kf+l  ,2) * 
     @                    MYYF(if,jf,kf+l,if+idir,jf+jdir,kf+l+1)
                  tmp = tmp
     @                 +sigmaf(if+ii-1,jf+ji-1,kf+l-1,3) * 
     @                    MZZF(if,jf,kf+l,if+idir,jf+jdir,kf+l-1)
     @                 +sigmaf(if+ii-1,jf+ji-1,kf+l  ,3) * 
     @                    MZZF(if,jf,kf+l,if+idir,jf+jdir,kf+l+1)
                  tmp = tmp + idir * (
     @                 +sigmaf(if+ii-1,jf+ji-1,kf+l-1,4) * 
     @                    MXZF(if,jf,kf+l,if+idir,jf+jdir,kf+l-1)
     @                 -sigmaf(if+ii-1,jf+ji-1,kf+l  ,4) * 
     @                    MXZF(if,jf,kf+l,if+idir,jf+jdir,kf+l+1) )
                  tmp = tmp + jdir * (
     @                 +sigmaf(if+ii-1,jf+ji-1,kf+l-1,5) * 
     @                    MYZF(if,jf,kf+l,if+idir,jf+jdir,kf+l-1)
     @                 -sigmaf(if+ii-1,jf+ji-1,kf+l  ,5) * 
     @                    MYZF(if,jf,kf+l,if+idir,jf+jdir,kf+l+1) )
                  cen =
     @               4.d0 * (sigmaf(if+ii-1,jf+ji-1,kf+l-1,1) 
     @                      +sigmaf(if+ii-1,jf+ji-1,kf+l  ,1)
     @                      +sigmaf(if+ii-1,jf+ji-1,kf+l-1,2) 
     @                      +sigmaf(if+ii-1,jf+ji-1,kf+l  ,2)
     @                      +sigmaf(if+ii-1,jf+ji-1,kf+l-1,3) 
     @                      +sigmaf(if+ii-1,jf+ji-1,kf+l  ,3))
     @             + 6.d0 * (sigmaf(if+ii-1,jf+ji-1,kf+l-1,4) 
     @                      -sigmaf(if+ii-1,jf+ji-1,kf+l  ,4)) * idir
     @             + 6.d0 * (sigmaf(if+ii-1,jf+ji-1,kf+l-1,5) 
     @                      -sigmaf(if+ii-1,jf+ji-1,kf+l  ,5)) * jdir

                  tmp = tmp - cen * fdst(if,jf,kf+l)

                        res(if,jf,kf) = res(if,jf,kf) + fac * tmp
 220                 continue

               else
                  do 230 kc = regl2, regh2
                     kf = kc * kr

                  tmp = sigmac(ic+ii-1,jc+ji-1,kc-1,1) * 
     @                    MXXC(ic+idir,jc+jdir,kc-1)
     @                 +sigmac(ic+ii-1,jc+ji-1,kc  ,1) * 
     @                    MXXC(ic+idir,jc+jdir,kc+1)
                  tmp = tmp
     @                 +sigmac(ic+ii-1,jc+ji-1,kc-1,2) * 
     @                    MYYC(ic+idir,jc+jdir,kc-1)
     @                 +sigmac(ic+ii-1,jc+ji-1,kc  ,2) * 
     @                    MYYC(ic+idir,jc+jdir,kc+1)
                  tmp = tmp
     @                 +sigmac(ic+ii-1,jc+ji-1,kc-1,3) * 
     @                    MZZC(ic+idir,jc+jdir,kc-1)
     @                 +sigmac(ic+ii-1,jc+ji-1,kc  ,3) * 
     @                    MZZC(ic+idir,jc+jdir,kc+1)
                  tmp = tmp + idir * (
     @                 +sigmac(ic+ii-1,jc+ji-1,kc-1,4) * 
     @                    MXZC(ic+idir,jc+jdir,kc-1)
     @                 -sigmac(ic+ii-1,jc+ji-1,kc  ,4) * 
     @                    MXZC(ic+idir,jc+jdir,kc+1) )
                  tmp = tmp + jdir * (
     @                 +sigmac(ic+ii-1,jc+ji-1,kc-1,5) * 
     @                    MYZC(ic+idir,jc+jdir,kc-1)
     @                 -sigmac(ic+ii-1,jc+ji-1,kc  ,5) * 
     @                    MYZC(ic+idir,jc+jdir,kc+1) )
                  cen =
     @               4.d0 * (sigmac(ic+ii-1,jc+ji-1,kc-1,1)
     @                      +sigmac(ic+ii-1,jc+ji-1,kc  ,1)
     @                      +sigmac(ic+ii-1,jc+ji-1,kc-1,2)
     @                      +sigmac(ic+ii-1,jc+ji-1,kc  ,2)
     @                      +sigmac(ic+ii-1,jc+ji-1,kc-1,3)
     @                      +sigmac(ic+ii-1,jc+ji-1,kc  ,3) )
     @             + 6.d0 * (sigmac(ic+ii-1,jc+ji-1,kc-1,4) 
     @                      -sigmac(ic+ii-1,jc+ji-1,kc  ,4) ) * idir
     @             + 6.d0 * (sigmac(ic+ii-1,jc+ji-1,kc-1,5)
     @                      -sigmac(ic+ii-1,jc+ji-1,kc  ,5) ) * jdir

                     tmp = tmp - cen * cdst(ic,jc,kc)
                     res(if,jf,kf) = res(if,jf,kf) + tmp
 230                 continue
               endif
 240           continue
c faces
c each face is two faces and two sides of an edge
         do 270 ji = 0, 1
            jdir = 2 * ji - 1
            do 270 ii = 0, 1
               idir = 2 * ii - 1
               if (ga(ii,ji,0) - ga(ii,1-ji,0) .eq. 1) then
                  fac0 = 1.d0 / (ir * kr)
                  do 250 l = -(kr-1), kr-1
                     fac1 = (kr-abs(l)) * fac0
                     do 250 m = idir, idir*(ir-1), idir
                        fac = (ir-abs(m)) * fac1
                        do 250 kf = kr*regl2, kr*regh2, kr

                  tmp = sigmaf(if+m-1,jf+ji-1,kf+l-1,1) * 
     @                    MXXF(if+m,jf,kf+l,if+m-1,jf+jdir,kf+l-1)
     @                 +sigmaf(if+m-1,jf+ji-1,kf+l  ,1) * 
     @                    MXXF(if+m,jf,kf+l,if+m-1,jf+jdir,kf+l+1)
     @                 +sigmaf(if+m  ,jf+ji-1,kf+l-1,1) * 
     @                    MXXF(if+m,jf,kf+l,if+m+1,jf+jdir,kf+l-1)
     @                 +sigmaf(if+m  ,jf+ji-1,kf+l  ,1) * 
     @                    MXXF(if+m,jf,kf+l,if+m+1,jf+jdir,kf+l+1)
                  tmp = tmp
     @                 +sigmaf(if+m-1,jf+ji-1,kf+l-1,2) * 
     @                    MYYF(if+m,jf,kf+l,if+m-1,jf+jdir,kf+l-1)
     @                 +sigmaf(if+m-1,jf+ji-1,kf+l  ,2) * 
     @                    MYYF(if+m,jf,kf+l,if+m-1,jf+jdir,kf+l+1)
     @                 +sigmaf(if+m  ,jf+ji-1,kf+l-1,2) * 
     @                    MYYF(if+m,jf,kf+l,if+m+1,jf+jdir,kf+l-1)
     @                 +sigmaf(if+m  ,jf+ji-1,kf+l  ,2) * 
     @                    MYYF(if+m,jf,kf+l,if+m+1,jf+jdir,kf+l+1)
                  tmp = tmp
     @                 +sigmaf(if+m-1,jf+ji-1,kf+l-1,3) * 
     @                    MZZF(if+m,jf,kf+l,if+m-1,jf+jdir,kf+l-1)
     @                 +sigmaf(if+m-1,jf+ji-1,kf+l  ,3) * 
     @                    MZZF(if+m,jf,kf+l,if+m-1,jf+jdir,kf+l+1)
     @                 +sigmaf(if+m  ,jf+ji-1,kf+l-1,3) * 
     @                    MZZF(if+m,jf,kf+l,if+m+1,jf+jdir,kf+l-1)
     @                 +sigmaf(if+m  ,jf+ji-1,kf+l  ,3) * 
     @                    MZZF(if+m,jf,kf+l,if+m+1,jf+jdir,kf+l+1)
                  tmp = tmp + (
     @                 -sigmaf(if+m-1,jf+ji-1,kf+l-1,4) * 
     @                    MXZF(if+m,jf,kf+l,if+m-1,jf+jdir,kf+l-1)
     @                 +sigmaf(if+m-1,jf+ji-1,kf+l  ,4) * 
     @                    MXZF(if+m,jf,kf+l,if+m-1,jf+jdir,kf+l+1)
     @                 +sigmaf(if+m  ,jf+ji-1,kf+l-1,4) * 
     @                    MXZF(if+m,jf,kf+l,if+m+1,jf+jdir,kf+l-1)
     @                 -sigmaf(if+m  ,jf+ji-1,kf+l  ,4) * 
     @                    MXZF(if+m,jf,kf+l,if+m+1,jf+jdir,kf+l+1) )
                  tmp = tmp + jdir * (
     @                 +sigmaf(if+m-1,jf+ji-1,kf+l-1,5) * 
     @                    MYZF(if+m,jf,kf+l,if+m-1,jf+jdir,kf+l-1)
     @                 -sigmaf(if+m-1,jf+ji-1,kf+l  ,5) * 
     @                    MYZF(if+m,jf,kf+l,if+m-1,jf+jdir,kf+l+1)
     @                 +sigmaf(if+m  ,jf+ji-1,kf+l-1,5) * 
     @                    MYZF(if+m,jf,kf+l,if+m+1,jf+jdir,kf+l-1)
     @                 -sigmaf(if+m  ,jf+ji-1,kf+l  ,5) * 
     @                    MYZF(if+m,jf,kf+l,if+m+1,jf+jdir,kf+l+1) )
                  cen =
     @               4.d0 * (sigmaf(if+m-1,jf+ji-1,kf+l-1,1)  
     @                      +sigmaf(if+m-1,jf+ji-1,kf+l  ,1)
     @                      +sigmaf(if+m  ,jf+ji-1,kf+l-1,1)  
     @                      +sigmaf(if+m  ,jf+ji-1,kf+l  ,1)
     @                      +sigmaf(if+m-1,jf+ji-1,kf+l-1,2)  
     @                      +sigmaf(if+m-1,jf+ji-1,kf+l  ,2)
     @                      +sigmaf(if+m  ,jf+ji-1,kf+l-1,2)  
     @                      +sigmaf(if+m  ,jf+ji-1,kf+l  ,2)
     @                      +sigmaf(if+m-1,jf+ji-1,kf+l-1,3)  
     @                      +sigmaf(if+m-1,jf+ji-1,kf+l  ,3)
     @                      +sigmaf(if+m  ,jf+ji-1,kf+l-1,3)  
     @                      +sigmaf(if+m  ,jf+ji-1,kf+l  ,3))
     @             + 6.d0 * (sigmaf(if+m-1,jf+ji-1,kf+l  ,4)  
     @                      -sigmaf(if+m-1,jf+ji-1,kf+l-1,4)
     @                      +sigmaf(if+m  ,jf+ji-1,kf+l-1,4) 
     @                      -sigmaf(if+m  ,jf+ji-1,kf+l  ,4))
     @             + 6.d0 * (sigmaf(if+m-1,jf+ji-1,kf+l-1,5)  
     @                      -sigmaf(if+m-1,jf+ji-1,kf+l  ,5)
     @                      +sigmaf(if+m  ,jf+ji-1,kf+l-1,5)  
     @                      -sigmaf(if+m  ,jf+ji-1,kf+l  ,5)) * jdir
                  tmp = tmp - cen * fdst(if+m,jf,kf+l)
                  res(if,jf,kf) = res(if,jf,kf) + fac * tmp

 250           continue
               endif
               if (ga(ii,ji,0) - ga(1-ii,ji,0) .eq. 1) then
                  fac0 = 1.d0 / (jr * kr)
                  do 260 l = -(kr-1), kr-1
                     fac1 = (kr-abs(l)) * fac0
                     do 260 n = jdir, jdir*(jr-1), jdir
                        fac = (jr-abs(n)) * fac1
                        do 260 kf = kr*regl2, kr*regh2, kr
                  tmp = sigmaf(if+ii-1,jf+n-1,kf+l-1,1) * 
     @                    MXXF(if,jf+n,kf+l,if+idir,jf+n-1,kf+l-1)
     @                 +sigmaf(if+ii-1,jf+n-1,kf+l  ,1) * 
     @                    MXXF(if,jf+n,kf+l,if+idir,jf+n-1,kf+l+1)
     @                 +sigmaf(if+ii-1,jf+n  ,kf+l-1,1) * 
     @                    MXXF(if,jf+n,kf+l,if+idir,jf+n+1,kf+l-1)
     @                 +sigmaf(if+ii-1,jf+n  ,kf+l  ,1) * 
     @                    MXXF(if,jf+n,kf+l,if+idir,jf+n+1,kf+l+1)
                  tmp = tmp
     @                 +sigmaf(if+ii-1,jf+n-1,kf+l-1,2) * 
     @                    MYYF(if,jf+n,kf+l,if+idir,jf+n-1,kf+l-1)
     @                 +sigmaf(if+ii-1,jf+n-1,kf+l  ,2) * 
     @                    MYYF(if,jf+n,kf+l,if+idir,jf+n-1,kf+l+1)
     @                 +sigmaf(if+ii-1,jf+n  ,kf+l-1,2) * 
     @                    MYYF(if,jf+n,kf+l,if+idir,jf+n+1,kf+l-1)
     @                 +sigmaf(if+ii-1,jf+n  ,kf+l  ,2) * 
     @                    MYYF(if,jf+n,kf+l,if+idir,jf+n+1,kf+l+1)
                  tmp = tmp
     @                 +sigmaf(if+ii-1,jf+n-1,kf+l-1,3) * 
     @                    MZZF(if,jf+n,kf+l,if+idir,jf+n-1,kf+l-1)
     @                 +sigmaf(if+ii-1,jf+n-1,kf+l  ,3) * 
     @                    MZZF(if,jf+n,kf+l,if+idir,jf+n-1,kf+l+1)
     @                 +sigmaf(if+ii-1,jf+n  ,kf+l-1,3) * 
     @                    MZZF(if,jf+n,kf+l,if+idir,jf+n+1,kf+l-1)
     @                 +sigmaf(if+ii-1,jf+n  ,kf+l  ,3) * 
     @                    MZZF(if,jf+n,kf+l,if+idir,jf+n+1,kf+l+1)
                  tmp = tmp + idir * (
     @                 +sigmaf(if+ii-1,jf+n-1,kf+l-1,4) * 
     @                    MXZF(if,jf+n,kf+l,if+idir,jf+n-1,kf+l-1)
     @                 -sigmaf(if+ii-1,jf+n-1,kf+l  ,4) * 
     @                    MXZF(if,jf+n,kf+l,if+idir,jf+n-1,kf+l+1)
     @                 +sigmaf(if+ii-1,jf+n  ,kf+l-1,4) * 
     @                    MXZF(if,jf+n,kf+l,if+idir,jf+n+1,kf+l-1)
     @                 -sigmaf(if+ii-1,jf+n  ,kf+l  ,4) * 
     @                    MXZF(if,jf+n,kf+l,if+idir,jf+n+1,kf+l+1) )
                  tmp = tmp + (
     @                 -sigmaf(if+ii-1,jf+n-1,kf+l-1,5) * 
     @                    MYZF(if,jf+n,kf+l,if+idir,jf+n-1,kf+l-1)
     @                 +sigmaf(if+ii-1,jf+n-1,kf+l  ,5) * 
     @                    MYZF(if,jf+n,kf+l,if+idir,jf+n-1,kf+l+1)
     @                 +sigmaf(if+ii-1,jf+n  ,kf+l-1,5) * 
     @                    MYZF(if,jf+n,kf+l,if+idir,jf+n+1,kf+l-1)
     @                 -sigmaf(if+ii-1,jf+n  ,kf+l  ,5) * 
     @                    MYZF(if,jf+n,kf+l,if+idir,jf+n+1,kf+l+1) )

                  cen =
     @               4.d0 * (sigmaf(if+ii-1,jf+n-1,kf+l-1,1) 
     @                      +sigmaf(if+ii-1,jf+n-1,kf+l  ,1)
     @                      +sigmaf(if+ii-1,jf+n  ,kf+l-1,1) 
     @                      +sigmaf(if+ii-1,jf+n  ,kf+l  ,1)
     @                      +sigmaf(if+ii-1,jf+n-1,kf+l-1,2) 
     @                      +sigmaf(if+ii-1,jf+n-1,kf+l  ,2)
     @                      +sigmaf(if+ii-1,jf+n  ,kf+l-1,2) 
     @                      +sigmaf(if+ii-1,jf+n  ,kf+l  ,2)
     @                      +sigmaf(if+ii-1,jf+n-1,kf+l-1,3) 
     @                      +sigmaf(if+ii-1,jf+n-1,kf+l  ,3)
     @                      +sigmaf(if+ii-1,jf+n  ,kf+l-1,3) 
     @                      +sigmaf(if+ii-1,jf+n  ,kf+l  ,3))
     @             + 6.d0 * (sigmaf(if+ii-1,jf+n-1,kf+l-1,4) 
     @                      -sigmaf(if+ii-1,jf+n-1,kf+l  ,4)
     @                      +sigmaf(if+ii-1,jf+n  ,kf+l-1,4) 
     @                      -sigmaf(if+ii-1,jf+n  ,kf+l  ,4)) * idir
     @             + 6.d0 * (sigmaf(if+ii-1,jf+n-1,kf+l  ,5) 
     @                      -sigmaf(if+ii-1,jf+n-1,kf+l-1,5)
     @                      +sigmaf(if+ii-1,jf+n  ,kf+l-1,5) 
     @                      -sigmaf(if+ii-1,jf+n  ,kf+l  ,5))
                  tmp = tmp - cen * fdst(if,jf+n,kf+l)
                  res(if,jf,kf) = res(if,jf,kf) + fac * tmp

 260           continue
               endif

 270           continue
c weighting
         do 280 kf = kr*regl2, kr*regh2, kr
 280        res(if,jf,kf) = src(if,jf,kf) - res(if,jf,kf) * fac00
      endif
      end

c Note---assumes fdst linearly interpolated from cdst along face
      subroutine hgcres(res, dims(res),
     @ src,  dims(src),
     @ fdst, dims(fdst),
     @ cdst, dims(cdst),
     @ sigmaf, dims(sf),
     @ sigmac, dims(sc),
     @ dims(reg),
     @ ir, jr, kr, ga)
      integer dims(res)
      integer dims(src)
      integer dims(fdst)
      integer dims(cdst)
      integer dims(sf)
      integer dims(sc)
      integer dims(reg)
      real*8 res(dimdec(res))
      real*8 src(dimdec(src))
      real*8 fdst(dimdec(fdst))
      real*8 cdst(dimdec(cdst))
      real*8 sigmaf(dimdec(sf),5)
      real*8 sigmac(dimdec(sc),5)
      integer ir, jr, kr, ga(0:1,0:1,0:1)
      real*8 cen, sum, fac00, fac0, fac1, fac2, fac
      integer ic, jc, kc, if, jf, kf, ii, ji, ki, idir, jdir, kdir
      integer l, m, n
      ic = regl0
      jc = regl1
      kc = regl2
      if = ic * ir
      jf = jc * jr
      kf = kc * kr
      sum = 0.d0

      fac00 = 1.d0 / 36.d0

c octants
      do 10 ki = 0, 1
         kdir = 2 * ki - 1
         do 10 ji = 0, 1
            jdir = 2 * ji - 1
            do 10 ii = 0, 1
               idir = 2 * ii - 1
               if (ga(ii,ji,ki) .eq. 1) then

                  sum = sum 
     @             +sigmaf(if+ii-1,jf+ji-1,kf+ki-1,1) * 
     @                MXXF(if,jf,kf,if+idir,jf+jdir,kf+kdir)
     @             +sigmaf(if+ii-1,jf+ji-1,kf+ki-1,2) * 
     @                MYYF(if,jf,kf,if+idir,jf+jdir,kf+kdir)
     @             +sigmaf(if+ii-1,jf+ji-1,kf+ki-1,3) * 
     @                MZZF(if,jf,kf,if+idir,jf+jdir,kf+kdir)
     @             -sigmaf(if+ii-1,jf+ji-1,kf+ki-1,4) * 
     @                MXZF(if,jf,kf,if+idir,jf+jdir,kf+kdir) * idir * kdir
     @             -sigmaf(if+ii-1,jf+ji-1,kf+ki-1,5) * 
     @                MYZF(if,jf,kf,if+idir,jf+jdir,kf+kdir) * jdir * kdir 

         cen =
     @      4.d0 * ( sigmaf(if+ii-1,jf+ji-1,kf+ki-1,1)
     @              +sigmaf(if+ii-1,jf+ji-1,kf+ki-1,2)
     @              +sigmaf(if+ii-1,jf+ji-1,kf+ki-1,3))
     @    + 6.d0 * (-sigmaf(if+ii-1,jf+ji-1,kf+ki-1,4)) * idir * kdir
     @    + 6.d0 * (-sigmaf(if+ii-1,jf+ji-1,kf+ki-1,5)) * jdir * kdir

                   sum = sum - cen * fdst(if,jf,kf)

               else

                  sum = sum 
     @             +sigmac(ic+ii-1,jc+ji-1,kc+ki-1,1) * 
     @                MXXC(ic+idir,jc+jdir,kc+kdir)
     @             +sigmac(ic+ii-1,jc+ji-1,kc+ki-1,2) * 
     @                MYYC(ic+idir,jc+jdir,kc+kdir)
     @             +sigmac(ic+ii-1,jc+ji-1,kc+ki-1,3) * 
     @                MZZC(ic+idir,jc+jdir,kc+kdir)
     @             -sigmac(ic+ii-1,jc+ji-1,kc+ki-1,4) * 
     @                MXZC(ic+idir,jc+jdir,kc+kdir) * idir * kdir
     @             -sigmac(ic+ii-1,jc+ji-1,kc+ki-1,5) * 
     @                MYZC(ic+idir,jc+jdir,kc+kdir) * jdir * kdir

         cen =
     @      4.d0 * ( sigmac(ic+ii-1,jc+ji-1,kc+ki-1,1)
     @              +sigmac(ic+ii-1,jc+ji-1,kc+ki-1,2)
     @              +sigmac(ic+ii-1,jc+ji-1,kc+ki-1,3))
     @    + 6.d0 * (-sigmac(ic+ii-1,jc+ji-1,kc+ki-1,4)) * idir * kdir
     @    + 6.d0 * (-sigmac(ic+ii-1,jc+ji-1,kc+ki-1,5)) * jdir * kdir

                   sum = sum - cen * cdst(ic,jc,kc)

               endif
 10            continue
c faces
      do 50 ki = 0, 1
         kdir = 2 * ki - 1
         do 50 ji = 0, 1
            jdir = 2 * ji - 1
            do 50 ii = 0, 1
               idir = 2 * ii - 1
               if (ga(ii,ji,ki) - ga(ii,ji,1-ki) .eq. 1) then
                  fac2 = 1.d0 / (ir * jr)
                  do 20 n = jdir, jdir*(jr-1), jdir
                     fac1 = (jr-abs(n)) * fac2
                     do 20 m = idir, idir*(ir-1), idir
                        fac = (ir-abs(m)) * fac1

         sum = sum + fac * (
     @     +sigmaf(if+m-1,jf+n-1,kf+ki-1,1) * 
     @        MXXF(if+m,jf+n,kf,if+m-1,jf+n-1,kf+kdir)
     @     +sigmaf(if+m-1,jf+n  ,kf+ki-1,1) * 
     @        MXXF(if+m,jf+n,kf,if+m-1,jf+n+1,kf+kdir)
     @     +sigmaf(if+m  ,jf+n-1,kf+ki-1,1) * 
     @        MXXF(if+m,jf+n,kf,if+m+1,jf+n-1,kf+kdir)
     @     +sigmaf(if+m  ,jf+n  ,kf+ki-1,1) * 
     @        MXXF(if+m,jf+n,kf,if+m+1,jf+n+1,kf+kdir))
         sum = sum + fac * (
     @     +sigmaf(if+m-1,jf+n-1,kf+ki-1,2) * 
     @        MYYF(if+m,jf+n,kf,if+m-1,jf+n-1,kf+kdir)
     @     +sigmaf(if+m-1,jf+n  ,kf+ki-1,2) * 
     @        MYYF(if+m,jf+n,kf,if+m-1,jf+n+1,kf+kdir)
     @     +sigmaf(if+m  ,jf+n-1,kf+ki-1,2) * 
     @        MYYF(if+m,jf+n,kf,if+m+1,jf+n-1,kf+kdir)
     @     +sigmaf(if+m  ,jf+n  ,kf+ki-1,2) * 
     @        MYYF(if+m,jf+n,kf,if+m+1,jf+n+1,kf+kdir))
         sum = sum + fac * (
     @     +sigmaf(if+m-1,jf+n-1,kf+ki-1,3) * 
     @        MZZF(if+m,jf+n,kf,if+m-1,jf+n-1,kf+kdir)
     @     +sigmaf(if+m-1,jf+n  ,kf+ki-1,3) * 
     @        MZZF(if+m,jf+n,kf,if+m-1,jf+n+1,kf+kdir)
     @     +sigmaf(if+m  ,jf+n-1,kf+ki-1,3) * 
     @        MZZF(if+m,jf+n,kf,if+m+1,jf+n-1,kf+kdir)
     @     +sigmaf(if+m  ,jf+n  ,kf+ki-1,3) * 
     @        MZZF(if+m,jf+n,kf,if+m+1,jf+n+1,kf+kdir))
         sum = sum + fac * kdir * (
     @     +sigmaf(if+m-1,jf+n-1,kf+ki-1,4) * 
     @        MXZF(if+m,jf+n,kf,if+m-1,jf+n-1,kf+kdir)
     @     +sigmaf(if+m-1,jf+n  ,kf+ki-1,4) * 
     @        MXZF(if+m,jf+n,kf,if+m-1,jf+n+1,kf+kdir)
     @     -sigmaf(if+m  ,jf+n-1,kf+ki-1,4) * 
     @        MXZF(if+m,jf+n,kf,if+m+1,jf+n-1,kf+kdir)
     @     -sigmaf(if+m  ,jf+n  ,kf+ki-1,4) * 
     @        MXZF(if+m,jf+n,kf,if+m+1,jf+n+1,kf+kdir))
         sum = sum + fac * kdir * (
     @     +sigmaf(if+m-1,jf+n-1,kf+ki-1,5) * 
     @        MYZF(if+m,jf+n,kf,if+m-1,jf+n-1,kf+kdir)
     @     -sigmaf(if+m-1,jf+n  ,kf+ki-1,5) * 
     @        MYZF(if+m,jf+n,kf,if+m-1,jf+n+1,kf+kdir)
     @     +sigmaf(if+m  ,jf+n-1,kf+ki-1,5) * 
     @        MYZF(if+m,jf+n,kf,if+m+1,jf+n-1,kf+kdir)
     @     -sigmaf(if+m  ,jf+n  ,kf+ki-1,5) * 
     @        MYZF(if+m,jf+n,kf,if+m+1,jf+n+1,kf+kdir))

         cen =
     @      4.d0 * ( sigmaf(if+m-1,jf+n-1,kf+ki-1,1)
     @              +sigmaf(if+m-1,jf+n  ,kf+ki-1,1)
     @              +sigmaf(if+m  ,jf+n-1,kf+ki-1,1)
     @              +sigmaf(if+m  ,jf+n  ,kf+ki-1,1)
     @              +sigmaf(if+m-1,jf+n-1,kf+ki-1,2)
     @              +sigmaf(if+m-1,jf+n  ,kf+ki-1,2)
     @              +sigmaf(if+m  ,jf+n-1,kf+ki-1,2)
     @              +sigmaf(if+m  ,jf+n  ,kf+ki-1,2)
     @              +sigmaf(if+m-1,jf+n-1,kf+ki-1,3)
     @              +sigmaf(if+m-1,jf+n  ,kf+ki-1,3)
     @              +sigmaf(if+m  ,jf+n-1,kf+ki-1,3)
     @              +sigmaf(if+m  ,jf+n  ,kf+ki-1,3))
     @    + 6.d0 * ( sigmaf(if+m-1,jf+n-1,kf+ki-1,4)
     @              +sigmaf(if+m-1,jf+n  ,kf+ki-1,4)
     @              -sigmaf(if+m  ,jf+n-1,kf+ki-1,4)
     @              -sigmaf(if+m  ,jf+n  ,kf+ki-1,4)) * kdir
     @    + 6.d0 * ( sigmaf(if+m-1,jf+n-1,kf+ki-1,5)
     @              -sigmaf(if+m-1,jf+n  ,kf+ki-1,5)
     @              +sigmaf(if+m  ,jf+n-1,kf+ki-1,5)
     @              -sigmaf(if+m  ,jf+n  ,kf+ki-1,5)) * kdir

          sum = sum - fac * cen * fdst(if+m,jf+n,kf)

 20               continue
               endif

               if (ga(ii,ji,ki) - ga(ii,1-ji,ki) .eq. 1) then
                  fac2 = 1.d0 / (ir * kr)
                  do 30 l = kdir, kdir*(kr-1), kdir
                     fac1 = (kr-abs(l)) * fac2
                     do 30 m = idir, idir*(ir-1), idir
                        fac = (ir-abs(m)) * fac1

         sum = sum + fac * (
     @                 +sigmaf(if+m-1,jf+ji-1,kf+l-1,1) * 
     @                    MXXF(if+m,jf,kf+l,if+m-1,jf+jdir,kf+l-1)
     @                 +sigmaf(if+m-1,jf+ji-1,kf+l  ,1) * 
     @                    MXXF(if+m,jf,kf+l,if+m-1,jf+jdir,kf+l+1)
     @                 +sigmaf(if+m  ,jf+ji-1,kf+l-1,1) * 
     @                    MXXF(if+m,jf,kf+l,if+m+1,jf+jdir,kf+l-1)
     @                 +sigmaf(if+m  ,jf+ji-1,kf+l  ,1) * 
     @                    MXXF(if+m,jf,kf+l,if+m+1,jf+jdir,kf+l+1) )
         sum = sum + fac * (
     @                 +sigmaf(if+m-1,jf+ji-1,kf+l-1,2) * 
     @                    MYYF(if+m,jf,kf+l,if+m-1,jf+jdir,kf+l-1)
     @                 +sigmaf(if+m-1,jf+ji-1,kf+l  ,2) * 
     @                    MYYF(if+m,jf,kf+l,if+m-1,jf+jdir,kf+l+1)
     @                 +sigmaf(if+m  ,jf+ji-1,kf+l-1,2) * 
     @                    MYYF(if+m,jf,kf+l,if+m+1,jf+jdir,kf+l-1)
     @                 +sigmaf(if+m  ,jf+ji-1,kf+l  ,2) * 
     @                    MYYF(if+m,jf,kf+l,if+m+1,jf+jdir,kf+l+1) )
         sum = sum + fac * (
     @                 +sigmaf(if+m-1,jf+ji-1,kf+l-1,3) * 
     @                    MZZF(if+m,jf,kf+l,if+m-1,jf+jdir,kf+l-1)
     @                 +sigmaf(if+m-1,jf+ji-1,kf+l  ,3) * 
     @                    MZZF(if+m,jf,kf+l,if+m-1,jf+jdir,kf+l+1)
     @                 +sigmaf(if+m  ,jf+ji-1,kf+l-1,3) * 
     @                    MZZF(if+m,jf,kf+l,if+m+1,jf+jdir,kf+l-1)
     @                 +sigmaf(if+m  ,jf+ji-1,kf+l  ,3) * 
     @                    MZZF(if+m,jf,kf+l,if+m+1,jf+jdir,kf+l+1) )
         sum = sum + fac * (
     @                 -sigmaf(if+m-1,jf+ji-1,kf+l-1,4) * 
     @                    MXZF(if+m,jf,kf+l,if+m-1,jf+jdir,kf+l-1)
     @                 +sigmaf(if+m-1,jf+ji-1,kf+l  ,4) * 
     @                    MXZF(if+m,jf,kf+l,if+m-1,jf+jdir,kf+l+1)
     @                 +sigmaf(if+m  ,jf+ji-1,kf+l-1,4) * 
     @                    MXZF(if+m,jf,kf+l,if+m+1,jf+jdir,kf+l-1)
     @                 -sigmaf(if+m  ,jf+ji-1,kf+l  ,4) * 
     @                    MXZF(if+m,jf,kf+l,if+m+1,jf+jdir,kf+l+1) )
         sum = sum + fac * jdir * (
     @                 +sigmaf(if+m-1,jf+ji-1,kf+l-1,5) * 
     @                    MYZF(if+m,jf,kf+l,if+m-1,jf+jdir,kf+l-1)
     @                 -sigmaf(if+m-1,jf+ji-1,kf+l  ,5) * 
     @                    MYZF(if+m,jf,kf+l,if+m-1,jf+jdir,kf+l+1)
     @                 +sigmaf(if+m  ,jf+ji-1,kf+l-1,5) * 
     @                    MYZF(if+m,jf,kf+l,if+m+1,jf+jdir,kf+l-1)
     @                 -sigmaf(if+m  ,jf+ji-1,kf+l  ,5) * 
     @                    MYZF(if+m,jf,kf+l,if+m+1,jf+jdir,kf+l+1) )

         cen =
     @      4.d0 * (
     @        +sigmaf(if+m-1,jf+ji-1,kf+l-1,1) + sigmaf(if+m-1,jf+ji-1,kf+l,1)
     @        +sigmaf(if+m  ,jf+ji-1,kf+l-1,1) + sigmaf(if+m  ,jf+ji-1,kf+l,1)
     @        +sigmaf(if+m-1,jf+ji-1,kf+l-1,2) + sigmaf(if+m-1,jf+ji-1,kf+l,2)
     @        +sigmaf(if+m  ,jf+ji-1,kf+l-1,2) + sigmaf(if+m  ,jf+ji-1,kf+l,2)
     @        +sigmaf(if+m-1,jf+ji-1,kf+l-1,3) + sigmaf(if+m-1,jf+ji-1,kf+l,3)
     @        +sigmaf(if+m  ,jf+ji-1,kf+l-1,3) + sigmaf(if+m  ,jf+ji-1,kf+l,3))
     @     +6.d0 * (
     @         sigmaf(if+m-1,jf+ji-1,kf+l  ,4) - sigmaf(if+m-1,jf+ji-1,kf+l-1,4)
     @        +sigmaf(if+m  ,jf+ji-1,kf+l-1,4) - sigmaf(if+m  ,jf+ji-1,kf+l,4))
     @     +6.d0 * jdir * (
     @         sigmaf(if+m-1,jf+ji-1,kf+l-1,5) - sigmaf(if+m-1,jf+ji-1,kf+l,5)
     @        +sigmaf(if+m  ,jf+ji-1,kf+l-1,5) - sigmaf(if+m  ,jf+ji-1,kf+l,5))

          sum = sum - fac * cen * fdst(if+m,jf,kf+l)

 30               continue
               endif
               if (ga(ii,ji,ki) - ga(1-ii,ji,ki) .eq. 1) then
                  fac2 = 1.d0 / (jr * kr)
                  do 40 l = kdir, kdir*(kr-1), kdir
                     fac1 = (kr-abs(l)) * fac2
                     do 40 n = jdir, jdir*(jr-1), jdir
                        fac = (jr-abs(n)) * fac1

         sum = sum + fac * (
     @                 +sigmaf(if+ii-1,jf+n-1,kf+l-1,1) * 
     @                    MXXF(if,jf+n,kf+l,if+idir,jf+n-1,kf+l-1)
     @                 +sigmaf(if+ii-1,jf+n-1,kf+l  ,1) * 
     @                    MXXF(if,jf+n,kf+l,if+idir,jf+n-1,kf+l+1)
     @                 +sigmaf(if+ii-1,jf+n  ,kf+l-1,1) * 
     @                    MXXF(if,jf+n,kf+l,if+idir,jf+n+1,kf+l-1)
     @                 +sigmaf(if+ii-1,jf+n  ,kf+l  ,1) * 
     @                    MXXF(if,jf+n,kf+l,if+idir,jf+n+1,kf+l+1) )
         sum = sum + fac * (
     @                 +sigmaf(if+ii-1,jf+n-1,kf+l-1,2) * 
     @                    MYYF(if,jf+n,kf+l,if+idir,jf+n-1,kf+l-1)
     @                 +sigmaf(if+ii-1,jf+n-1,kf+l  ,2) * 
     @                    MYYF(if,jf+n,kf+l,if+idir,jf+n-1,kf+l+1)
     @                 +sigmaf(if+ii-1,jf+n  ,kf+l-1,2) * 
     @                    MYYF(if,jf+n,kf+l,if+idir,jf+n+1,kf+l-1)
     @                 +sigmaf(if+ii-1,jf+n  ,kf+l  ,2) * 
     @                    MYYF(if,jf+n,kf+l,if+idir,jf+n+1,kf+l+1) )
         sum = sum + fac * (
     @                 +sigmaf(if+ii-1,jf+n-1,kf+l-1,3) * 
     @                    MZZF(if,jf+n,kf+l,if+idir,jf+n-1,kf+l-1)
     @                 +sigmaf(if+ii-1,jf+n-1,kf+l  ,3) * 
     @                    MZZF(if,jf+n,kf+l,if+idir,jf+n-1,kf+l+1)
     @                 +sigmaf(if+ii-1,jf+n  ,kf+l-1,3) * 
     @                    MZZF(if,jf+n,kf+l,if+idir,jf+n+1,kf+l-1)
     @                 +sigmaf(if+ii-1,jf+n  ,kf+l  ,3) * 
     @                    MZZF(if,jf+n,kf+l,if+idir,jf+n+1,kf+l+1) )
         sum = sum + fac * idir * (
     @                 +sigmaf(if+ii-1,jf+n-1,kf+l-1,4) * 
     @                    MXZF(if,jf+n,kf+l,if+idir,jf+n-1,kf+l-1)
     @                 -sigmaf(if+ii-1,jf+n-1,kf+l  ,4) * 
     @                    MXZF(if,jf+n,kf+l,if+idir,jf+n-1,kf+l+1)
     @                 +sigmaf(if+ii-1,jf+n  ,kf+l-1,4) * 
     @                    MXZF(if,jf+n,kf+l,if+idir,jf+n+1,kf+l-1)
     @                 -sigmaf(if+ii-1,jf+n  ,kf+l  ,4) * 
     @                    MXZF(if,jf+n,kf+l,if+idir,jf+n+1,kf+l+1) )
         sum = sum + fac * (
     @                 -sigmaf(if+ii-1,jf+n-1,kf+l-1,5) * 
     @                    MYZF(if,jf+n,kf+l,if+idir,jf+n-1,kf+l-1)
     @                 +sigmaf(if+ii-1,jf+n-1,kf+l  ,5) * 
     @                    MYZF(if,jf+n,kf+l,if+idir,jf+n-1,kf+l+1)
     @                 +sigmaf(if+ii-1,jf+n  ,kf+l-1,5) * 
     @                    MYZF(if,jf+n,kf+l,if+idir,jf+n+1,kf+l-1)
     @                 -sigmaf(if+ii-1,jf+n  ,kf+l  ,5) * 
     @                    MYZF(if,jf+n,kf+l,if+idir,jf+n+1,kf+l+1) )

         cen =
     @      4.d0 * (
     @        sigmaf(if+ii-1,jf+n-1,kf+l-1,1) + sigmaf(if+ii-1,jf+n-1,kf+l  ,1)
     @       +sigmaf(if+ii-1,jf+n  ,kf+l-1,1) + sigmaf(if+ii-1,jf+n  ,kf+l  ,1)
     @       +sigmaf(if+ii-1,jf+n-1,kf+l-1,2) + sigmaf(if+ii-1,jf+n-1,kf+l  ,2)
     @       +sigmaf(if+ii-1,jf+n  ,kf+l-1,2) + sigmaf(if+ii-1,jf+n  ,kf+l  ,2)
     @       +sigmaf(if+ii-1,jf+n-1,kf+l-1,3) + sigmaf(if+ii-1,jf+n-1,kf+l  ,3)
     @       +sigmaf(if+ii-1,jf+n  ,kf+l-1,3) + sigmaf(if+ii-1,jf+n  ,kf+l  ,3))
     @    + 6.d0 * idir * (
     @        sigmaf(if+ii-1,jf+n-1,kf+l-1,4) - sigmaf(if+ii-1,jf+n-1,kf+l  ,4)
     @       +sigmaf(if+ii-1,jf+n  ,kf+l-1,4) - sigmaf(if+ii-1,jf+n  ,kf+l  ,4))
     @    + 6.d0 * (
     @        sigmaf(if+ii-1,jf+n-1,kf+l  ,5) - sigmaf(if+ii-1,jf+n-1,kf+l-1,5)
     @       +sigmaf(if+ii-1,jf+n  ,kf+l-1,5) - sigmaf(if+ii-1,jf+n  ,kf+l  ,5))

          sum = sum - fac * cen * fdst(if,jf+n,kf+l)

 40               continue

            endif
 50         continue
c edges
      do 90 ki = 0, 1
         kdir = 2 * ki - 1
         do 90 ji = 0, 1
            jdir = 2 * ji - 1
            do 90 ii = 0, 1
               idir = 2 * ii - 1

               if (ga(ii,ji,ki) -
     @             min(ga(ii,ji,1-ki), ga(ii,1-ji,ki), ga(ii,1-ji,1-ki))
     @             .eq. 1) then
                  fac1 = 1.d0 / ir
                  do 60 m = idir, idir*(ir-1), idir
                     fac = (ir-abs(m)) * fac1

                  sum = sum + fac * (
     @                 +sigmaf(if+m-1,jf+ji-1,kf+ki-1,1) * 
     @                    MXXF(if+m,jf,kf,if+m-1,jf+jdir,kf+kdir)
     @                 +sigmaf(if+m  ,jf+ji-1,kf+ki-1,1) * 
     @                    MXXF(if+m,jf,kf,if+m+1,jf+jdir,kf+kdir) )
                  sum = sum + fac * (
     @                 +sigmaf(if+m-1,jf+ji-1,kf+ki-1,2) * 
     @                    MYYF(if+m,jf,kf,if+m-1,jf+jdir,kf+kdir)
     @                 +sigmaf(if+m  ,jf+ji-1,kf+ki-1,2) * 
     @                    MYYF(if+m,jf,kf,if+m+1,jf+jdir,kf+kdir) )
                  sum = sum + fac * (
     @                 +sigmaf(if+m-1,jf+ji-1,kf+ki-1,3) * 
     @                    MZZF(if+m,jf,kf,if+m-1,jf+jdir,kf+kdir)
     @                 +sigmaf(if+m  ,jf+ji-1,kf+ki-1,3) * 
     @                    MZZF(if+m,jf,kf,if+m+1,jf+jdir,kf+kdir) )
                  sum = sum + fac * kdir * (
     @                 +sigmaf(if+m-1,jf+ji-1,kf+ki-1,4) * 
     @                    MXZF(if+m,jf,kf,if+m-1,jf+jdir,kf+kdir)
     @                 -sigmaf(if+m  ,jf+ji-1,kf+ki-1,4) * 
     @                    MXZF(if+m,jf,kf,if+m+1,jf+jdir,kf+kdir) )
                  sum = sum  + fac * jdir * kdir * (
     @                 -sigmaf(if+m-1,jf+ji-1,kf+ki-1,5) * 
     @                    MYZF(if+m,jf,kf,if+m-1,jf+jdir,kf+kdir)
     @                 -sigmaf(if+m  ,jf+ji-1,kf+ki-1,5) * 
     @                    MYZF(if+m,jf,kf,if+m+1,jf+jdir,kf+kdir) )

                  cen =
     @               4.d0 * ( sigmaf(if+m-1,jf+ji-1,kf+ki-1,1)
     @                       +sigmaf(if+m  ,jf+ji-1,kf+ki-1,1)
     @                       +sigmaf(if+m-1,jf+ji-1,kf+ki-1,2)
     @                       +sigmaf(if+m  ,jf+ji-1,kf+ki-1,2)
     @                       +sigmaf(if+m-1,jf+ji-1,kf+ki-1,3)
     @                       +sigmaf(if+m  ,jf+ji-1,kf+ki-1,3))
     @             + 6.d0 * ( sigmaf(if+m-1,jf+ji-1,kf+ki-1,4)
     @                       -sigmaf(if+m  ,jf+ji-1,kf+ki-1,4)) * kdir 
     @             + 6.d0 * (-sigmaf(if+m-1,jf+ji-1,kf+ki-1,5)
     @                       -sigmaf(if+m  ,jf+ji-1,kf+ki-1,5)) * jdir * kdir 

                  sum = sum - fac * cen * fdst(if+m,jf,kf)
 60               continue
               endif

               if (ga(ii,ji,ki) -
     @             min(ga(ii,ji,1-ki), ga(1-ii,ji,ki), ga(1-ii,ji,1-ki))
     @             .eq. 1) then
                  fac1 = 1.d0 / jr
                  do 70 n = jdir, jdir*(jr-1), jdir
                     fac = (jr-abs(n)) * fac1

                  sum = sum + fac * (
     @                 +sigmaf(if+ii-1,jf+n-1,kf+ki-1,1) * 
     @                    MXXF(if,jf+n,kf,if+idir,jf+n-1,kf+kdir)
     @                 +sigmaf(if+ii-1,jf+n  ,kf+ki-1,1) * 
     @                    MXXF(if,jf+n,kf,if+idir,jf+n+1,kf+kdir) )
                  sum = sum + fac * (
     @                 +sigmaf(if+ii-1,jf+n-1,kf+ki-1,2) * 
     @                    MYYF(if,jf+n,kf,if+idir,jf+n-1,kf+kdir)
     @                 +sigmaf(if+ii-1,jf+n  ,kf+ki-1,2) * 
     @                    MYYF(if,jf+n,kf,if+idir,jf+n+1,kf+kdir) )
                  sum = sum + fac * (
     @                 +sigmaf(if+ii-1,jf+n-1,kf+ki-1,3) * 
     @                    MZZF(if,jf+n,kf,if+idir,jf+n-1,kf+kdir)
     @                 +sigmaf(if+ii-1,jf+n  ,kf+ki-1,3) * 
     @                    MZZF(if,jf+n,kf,if+idir,jf+n+1,kf+kdir) )
                  sum = sum + fac * idir * kdir * (
     @                 -sigmaf(if+ii-1,jf+n-1,kf+ki-1,4) * 
     @                    MXZF(if,jf+n,kf,if+idir,jf+n-1,kf+kdir)
     @                 -sigmaf(if+ii-1,jf+n  ,kf+ki-1,4) * 
     @                    MXZF(if,jf+n,kf,if+idir,jf+n+1,kf+kdir) )
                  sum = sum  + fac * kdir * (
     @                 +sigmaf(if+ii-1,jf+n-1,kf+ki-1,5) * 
     @                    MYZF(if,jf+n,kf,if+idir,jf+n-1,kf+kdir)
     @                 -sigmaf(if+ii-1,jf+n  ,kf+ki-1,5) * 
     @                    MYZF(if,jf+n,kf,if+idir,jf+n+1,kf+kdir) )
                  cen =
     @                  4.d0 * ( sigmaf(if+ii-1,jf+n-1,kf+ki-1,1)
     @                          +sigmaf(if+ii-1,jf+n  ,kf+ki-1,1)
     @                          +sigmaf(if+ii-1,jf+n-1,kf+ki-1,2)
     @                          +sigmaf(if+ii-1,jf+n  ,kf+ki-1,2)
     @                          +sigmaf(if+ii-1,jf+n-1,kf+ki-1,3)
     @                          +sigmaf(if+ii-1,jf+n  ,kf+ki-1,3))
     @                + 6.d0 * (-sigmaf(if+ii-1,jf+n-1,kf+ki-1,4)
     @                          -sigmaf(if+ii-1,jf+n  ,kf+ki-1,4)) * idir * kdir
     @                + 6.d0 * ( sigmaf(if+ii-1,jf+n-1,kf+ki-1,5) 
     @                          -sigmaf(if+ii-1,jf+n  ,kf+ki-1,5)) * kdir
                  sum = sum - fac * cen * fdst(if,jf+n,kf)

 70               continue

               endif

               if (ga(ii,ji,ki) -
     @             min(ga(ii,1-ji,ki), ga(1-ii,ji,ki), ga(1-ii,1-ji,ki))
     @             .eq. 1) then
                  fac1 = 1.d0 / kr
                  do 80 l = kdir, kdir*(kr-1), kdir
                     fac = (kr-abs(l)) * fac1

                  sum = sum + fac * (
     @                 +sigmaf(if+ii-1,jf+ji-1,kf+l-1,1) * 
     @                    MXXF(if,jf,kf+l,if+idir,jf+jdir,kf+l-1)
     @                 +sigmaf(if+ii-1,jf+ji-1,kf+l  ,1) * 
     @                    MXXF(if,jf,kf+l,if+idir,jf+jdir,kf+l+1) )
                  sum = sum + fac * (
     @                 +sigmaf(if+ii-1,jf+ji-1,kf+l-1,2) * 
     @                    MYYF(if,jf,kf+l,if+idir,jf+jdir,kf+l-1)
     @                 +sigmaf(if+ii-1,jf+ji-1,kf+l  ,2) * 
     @                    MYYF(if,jf,kf+l,if+idir,jf+jdir,kf+l+1) )
                  sum = sum + fac * (
     @                 +sigmaf(if+ii-1,jf+ji-1,kf+l-1,3) * 
     @                    MZZF(if,jf,kf+l,if+idir,jf+jdir,kf+l-1)
     @                 +sigmaf(if+ii-1,jf+ji-1,kf+l  ,3) * 
     @                    MZZF(if,jf,kf+l,if+idir,jf+jdir,kf+l+1) )
                  sum = sum + fac * idir * (
     @                 +sigmaf(if+ii-1,jf+ji-1,kf+l-1,4) * 
     @                    MXZF(if,jf,kf+l,if+idir,jf+jdir,kf+l-1)
     @                 -sigmaf(if+ii-1,jf+ji-1,kf+l  ,4) * 
     @                    MXZF(if,jf,kf+l,if+idir,jf+jdir,kf+l+1) )
                  sum = sum + fac * jdir * (
     @                 +sigmaf(if+ii-1,jf+ji-1,kf+l-1,5) * 
     @                    MYZF(if,jf,kf+l,if+idir,jf+jdir,kf+l-1)
     @                 -sigmaf(if+ii-1,jf+ji-1,kf+l  ,5) * 
     @                    MYZF(if,jf,kf+l,if+idir,jf+jdir,kf+l+1) )
                  cen =
     @               4.d0 * (sigmaf(if+ii-1,jf+ji-1,kf+l-1,1) 
     @                      +sigmaf(if+ii-1,jf+ji-1,kf+l  ,1)
     @                      +sigmaf(if+ii-1,jf+ji-1,kf+l-1,2) 
     @                      +sigmaf(if+ii-1,jf+ji-1,kf+l  ,2)
     @                      +sigmaf(if+ii-1,jf+ji-1,kf+l-1,3) 
     @                      +sigmaf(if+ii-1,jf+ji-1,kf+l  ,3))
     @             + 6.d0 * (sigmaf(if+ii-1,jf+ji-1,kf+l-1,4) 
     @                      -sigmaf(if+ii-1,jf+ji-1,kf+l  ,4)) * idir
     @             + 6.d0 * (sigmaf(if+ii-1,jf+ji-1,kf+l-1,5) 
     @                      -sigmaf(if+ii-1,jf+ji-1,kf+l  ,5)) * jdir

                  sum = sum - fac * cen * fdst(if,jf,kf+l)

 80               continue

               endif
 90            continue

c weighting
      res(if,jf,kf) = src(if,jf,kf) - sum * fac00
      end

c end of terrain stencils
#else 

c variable density versions:

c Note---assumes fdst linearly interpolated from cdst along face
      subroutine hgfres(res, dims(res),
     @ src,  dims(src),
     @ fdst, dims(fdst),
     @ cdst, dims(cdst),
     @ sigmaf, dims(sf),
     @ sigmac, dims(sc),
     @ dims(reg),
     @ hx, hy, hz, ir, jr, kr, idim, idir)
      integer dims(res)
      integer dims(src)
      integer dims(fdst)
      integer dims(cdst)
      integer dims(sf)
      integer dims(sc)
      integer dims(reg)
      real*8 res(dimdec(res))
      real*8 src(dimdec(src))
      real*8 fdst(dimdec(fdst))
      real*8 cdst(dimdec(cdst))
      real*8 sigmaf(dimdec(sf))
      real*8 sigmac(dimdec(sc))
      real*8 hx, hy, hz
      integer ir, jr, kr, idim, idir
      real*8 hxm2, hym2, hzm2, fac0, fac1, fac2, tmp
      integer i, j, k, is, js, ks, l, m, n
      if (idim .eq. 0) then
         i = regl0
         if (idir .eq. 1) then
            is = i - 1
         else
            is = i
         endif
         fac0 = 0.5d0 * ir / (ir + 1.d0)
         hxm2 = 1.d0 / (ir * ir * hx * hx)
         hym2 = 1.d0 / (jr * jr * hy * hy)
         hzm2 = 1.d0 / (kr * kr * hz * hz)
         do 10 k = regl2, regh2
            do 10 j = regl1, regh1
 10            res(i*ir,j*jr,k*kr) =
     @           src(i*ir,j*jr,k*kr) - fac0 *
     @             (hxm2 *
     @               ((sigmac(is,j-1,k-1) + sigmac(is,j-1,k) +
     @                 sigmac(is,j,k-1)   + sigmac(is,j,k)) *
     @                 (cdst(i-idir,j,k) - cdst(i,j,k))) +
     @              hym2 *
     @               ((sigmac(is,j-1,k-1) + sigmac(is,j-1,k)) *
     @                 (cdst(i,j-1,k) - cdst(i,j,k)) +
     @                (sigmac(is,j,k-1) + sigmac(is,j,k)) *
     @                 (cdst(i,j+1,k) - cdst(i,j,k))) +
     @              hzm2 *
     @               ((sigmac(is,j-1,k-1) + sigmac(is,j,k-1)) *
     @                 (cdst(i,j,k-1) - cdst(i,j,k)) +
     @                (sigmac(is,j-1,k) + sigmac(is,j,k)) *
     @                 (cdst(i,j,k+1) - cdst(i,j,k))))
         fac0 = fac0 / (ir * jr * kr * jr * kr)
         hxm2 = ir * ir * hxm2
         hym2 = jr * jr * hym2
         hzm2 = kr * kr * hzm2
         i = i * ir
         if (idir .eq. 1) then
            is = i
         else
            is = i - 1
         endif
         do 20 l = 0, kr-1
            fac2 = (kr-l) * fac0
            if (l .eq. 0) fac2 = 0.5d0 * fac2
            do 20 n = 0, jr-1
               fac1 = (jr-n) * fac2
               if (n .eq. 0) fac1 = 0.5d0 * fac1
               do 20 k = kr*regl2, kr*regh2, kr
                  do 20 j = jr*regl1, jr*regh1, jr
                     tmp = hxm2 *
     @                 ((sigmaf(is,j-n-1,k-l-1) + sigmaf(is,j-n-1,k-l) +
     @                   sigmaf(is,j-n,k-l-1) + sigmaf(is,j-n,k-l)) *
     @                   (fdst(i+idir,j-n,k-l) - fdst(i,j-n,k-l)) +
     @                  (sigmaf(is,j-n-1,k+l-1) + sigmaf(is,j-n-1,k+l) +
     @                   sigmaf(is,j-n,k+l-1) + sigmaf(is,j-n,k+l)) *
     @                   (fdst(i+idir,j-n,k+l) - fdst(i,j-n,k+l)) +
     @                  (sigmaf(is,j+n-1,k-l-1) + sigmaf(is,j+n-1,k-l) +
     @                   sigmaf(is,j+n,k-l-1) + sigmaf(is,j+n,k-l)) *
     @                   (fdst(i+idir,j+n,k-l) - fdst(i,j+n,k-l)) +
     @                  (sigmaf(is,j+n-1,k+l-1) + sigmaf(is,j+n-1,k+l) +
     @                   sigmaf(is,j+n,k+l-1) + sigmaf(is,j+n,k+l)) *
     @                   (fdst(i+idir,j+n,k+l) - fdst(i,j+n,k+l)))
                     tmp = tmp + hym2 *
     @                 ((sigmaf(is,j-n-1,k-l-1) + sigmaf(is,j-n-1,k-l)) *
     @                   (fdst(i,j-n-1,k-l) - fdst(i,j-n,k-l)) +
     @                  (sigmaf(is,j-n,k-l-1) + sigmaf(is,j-n,k-l)) *
     @                   (fdst(i,j-n+1,k-l) - fdst(i,j-n,k-l)) +
     @                  (sigmaf(is,j-n-1,k+l-1) + sigmaf(is,j-n-1,k+l)) *
     @                   (fdst(i,j-n-1,k+l) - fdst(i,j-n,k+l)) +
     @                  (sigmaf(is,j-n,k+l-1) + sigmaf(is,j-n,k+l)) *
     @                   (fdst(i,j-n+1,k+l) - fdst(i,j-n,k+l)) +
     @                  (sigmaf(is,j+n-1,k-l-1) + sigmaf(is,j+n-1,k-l)) *
     @                   (fdst(i,j+n-1,k-l) - fdst(i,j+n,k-l)) +
     @                  (sigmaf(is,j+n,k-l-1) + sigmaf(is,j+n,k-l)) *
     @                   (fdst(i,j+n+1,k-l) - fdst(i,j+n,k-l)) +
     @                  (sigmaf(is,j+n-1,k+l-1) + sigmaf(is,j+n-1,k+l)) *
     @                   (fdst(i,j+n-1,k+l) - fdst(i,j+n,k+l)) +
     @                  (sigmaf(is,j+n,k+l-1) + sigmaf(is,j+n,k+l)) *
     @                   (fdst(i,j+n+1,k+l) - fdst(i,j+n,k+l)))
 20            res(i,j,k) = res(i,j,k) - fac1 * (tmp + hzm2 *
     @                 ((sigmaf(is,j-n-1,k-l-1) + sigmaf(is,j-n,k-l-1)) *
     @                   (fdst(i,j-n,k-l-1) - fdst(i,j-n,k-l)) +
     @                  (sigmaf(is,j-n-1,k-l) + sigmaf(is,j-n,k-l)) *
     @                   (fdst(i,j-n,k-l+1) - fdst(i,j-n,k-l)) +
     @                  (sigmaf(is,j-n-1,k+l-1) + sigmaf(is,j-n,k+l-1)) *
     @                   (fdst(i,j-n,k+l-1) - fdst(i,j-n,k+l)) +
     @                  (sigmaf(is,j-n-1,k+l) + sigmaf(is,j-n,k+l)) *
     @                   (fdst(i,j-n,k+l+1) - fdst(i,j-n,k+l)) +
     @                  (sigmaf(is,j+n-1,k-l-1) + sigmaf(is,j+n,k-l-1)) *
     @                   (fdst(i,j+n,k-l-1) - fdst(i,j+n,k-l)) +
     @                  (sigmaf(is,j+n-1,k-l) + sigmaf(is,j+n,k-l)) *
     @                   (fdst(i,j+n,k-l+1) - fdst(i,j+n,k-l)) +
     @                  (sigmaf(is,j+n-1,k+l-1) + sigmaf(is,j+n,k+l-1)) *
     @                   (fdst(i,j+n,k+l-1) - fdst(i,j+n,k+l)) +
     @                  (sigmaf(is,j+n-1,k+l) + sigmaf(is,j+n,k+l)) *
     @                   (fdst(i,j+n,k+l+1) - fdst(i,j+n,k+l))))
      else if (idim .eq. 1) then
         j = regl1
         if (idir .eq. 1) then
            js = j - 1
         else
            js = j
         endif
         fac0 = 0.5d0 * jr / (jr + 1.d0)
         hxm2 = 1.d0 / (ir * ir * hx * hx)
         hym2 = 1.d0 / (jr * jr * hy * hy)
         hzm2 = 1.d0 / (kr * kr * hz * hz)
         do 30 k = regl2, regh2
            do 30 i = regl0, regh0
 30            res(i*ir,j*jr,k*kr) =
     @           src(i*ir,j*jr,k*kr) - fac0 *
     @             (hxm2 *
     @               ((sigmac(i-1,js,k-1) + sigmac(i-1,js,k)) *
     @                 (cdst(i-1,j,k) - cdst(i,j,k)) +
     @                (sigmac(i,js,k-1) + sigmac(i,js,k)) *
     @                 (cdst(i+1,j,k) - cdst(i,j,k))) +
     @              hym2 *
     @               ((sigmac(i-1,js,k-1) + sigmac(i-1,js,k) +
     @                 sigmac(i,js,k-1)   + sigmac(i,js,k)) *
     @                 (cdst(i,j-idir,k) - cdst(i,j,k))) +
     @              hzm2 *
     @               ((sigmac(i-1,js,k-1) + sigmac(i,js,k-1)) *
     @                 (cdst(i,j,k-1) - cdst(i,j,k)) +
     @                (sigmac(i-1,js,k) + sigmac(i,js,k)) *
     @                 (cdst(i,j,k+1) - cdst(i,j,k))))
         fac0 = fac0 / (ir * jr * kr * ir * kr)
         hxm2 = ir * ir * hxm2
         hym2 = jr * jr * hym2
         hzm2 = kr * kr * hzm2
         j = j * jr
         if (idir .eq. 1) then
            js = j
         else
            js = j - 1
         endif
         do 40 l = 0, kr-1
            fac2 = (kr-l) * fac0
            if (l .eq. 0) fac2 = 0.5d0 * fac2
            do 40 m = 0, ir-1
               fac1 = (ir-m) * fac2
               if (m .eq. 0) fac1 = 0.5d0 * fac1
               do 40 k = kr*regl2, kr*regh2, kr
                  do 40 i = ir*regl0, ir*regh0, ir
                     tmp = hxm2 *
     @                 ((sigmaf(i-m-1,js,k-l-1) + sigmaf(i-m-1,js,k-l)) *
     @                   (fdst(i-m-1,j,k-l) - fdst(i-m,j,k-l)) +
     @                  (sigmaf(i-m,js,k-l-1) + sigmaf(i-m,js,k-l)) *
     @                   (fdst(i-m+1,j,k-l) - fdst(i-m,j,k-l)) +
     @                  (sigmaf(i-m-1,js,k+l-1) + sigmaf(i-m-1,js,k+l)) *
     @                   (fdst(i-m-1,j,k+l) - fdst(i-m,j,k+l)) +
     @                  (sigmaf(i-m,js,k+l-1) + sigmaf(i-m,js,k+l)) *
     @                   (fdst(i-m+1,j,k+l) - fdst(i-m,j,k+l)) +
     @                  (sigmaf(i+m-1,js,k-l-1) + sigmaf(i+m-1,js,k-l)) *
     @                   (fdst(i+m-1,j,k-l) - fdst(i+m,j,k-l)) +
     @                  (sigmaf(i+m,js,k-l-1) + sigmaf(i+m,js,k-l)) *
     @                   (fdst(i+m+1,j,k-l) - fdst(i+m,j,k-l)) +
     @                  (sigmaf(i+m-1,js,k+l-1) + sigmaf(i+m-1,js,k+l)) *
     @                   (fdst(i+m-1,j,k+l) - fdst(i+m,j,k+l)) +
     @                  (sigmaf(i+m,js,k+l-1) + sigmaf(i+m,js,k+l)) *
     @                   (fdst(i+m+1,j,k+l) - fdst(i+m,j,k+l)))
                     tmp = tmp + hym2 *
     @                 ((sigmaf(i-m-1,js,k-l-1) + sigmaf(i-m-1,js,k-l) +
     @                   sigmaf(i-m,js,k-l-1) + sigmaf(i-m,js,k-l)) *
     @                   (fdst(i-m,j+idir,k-l) - fdst(i-m,j,k-l)) +
     @                  (sigmaf(i-m-1,js,k+l-1) + sigmaf(i-m-1,js,k+l) +
     @                   sigmaf(i-m,js,k+l-1) + sigmaf(i-m,js,k+l)) *
     @                   (fdst(i-m,j+idir,k+l) - fdst(i-m,j,k+l)) +
     @                  (sigmaf(i+m-1,js,k-l-1) + sigmaf(i+m-1,js,k-l) +
     @                   sigmaf(i+m,js,k-l-1) + sigmaf(i+m,js,k-l)) *
     @                   (fdst(i+m,j+idir,k-l) - fdst(i+m,j,k-l)) +
     @                  (sigmaf(i+m-1,js,k+l-1) + sigmaf(i+m-1,js,k+l) +
     @                   sigmaf(i+m,js,k+l-1) + sigmaf(i+m,js,k+l)) *
     @                   (fdst(i+m,j+idir,k+l) - fdst(i+m,j,k+l)))
 40            res(i,j,k) = res(i,j,k) - fac1 * (tmp + hzm2 *
     @                 ((sigmaf(i-m-1,js,k-l-1) + sigmaf(i-m,js,k-l-1)) *
     @                   (fdst(i-m,j,k-l-1) - fdst(i-m,j,k-l)) +
     @                  (sigmaf(i-m-1,js,k-l) + sigmaf(i-m,js,k-l)) *
     @                   (fdst(i-m,j,k-l+1) - fdst(i-m,j,k-l)) +
     @                  (sigmaf(i-m-1,js,k+l-1) + sigmaf(i-m,js,k+l-1)) *
     @                   (fdst(i-m,j,k+l-1) - fdst(i-m,j,k+l)) +
     @                  (sigmaf(i-m-1,js,k+l) + sigmaf(i-m,js,k+l)) *
     @                   (fdst(i-m,j,k+l+1) - fdst(i-m,j,k+l)) +
     @                  (sigmaf(i+m-1,js,k-l-1) + sigmaf(i+m,js,k-l-1)) *
     @                   (fdst(i+m,j,k-l-1) - fdst(i+m,j,k-l)) +
     @                  (sigmaf(i+m-1,js,k-l) + sigmaf(i+m,js,k-l)) *
     @                   (fdst(i+m,j,k-l+1) - fdst(i+m,j,k-l)) +
     @                  (sigmaf(i+m-1,js,k+l-1) + sigmaf(i+m,js,k+l-1)) *
     @                   (fdst(i+m,j,k+l-1) - fdst(i+m,j,k+l)) +
     @                  (sigmaf(i+m-1,js,k+l) + sigmaf(i+m,js,k+l)) *
     @                   (fdst(i+m,j,k+l+1) - fdst(i+m,j,k+l))))
      else
         k = regl2
         if (idir .eq. 1) then
            ks = k - 1
         else
            ks = k
         endif
         fac0 = 0.5d0 * kr / (kr + 1.d0)
         hxm2 = 1.d0 / (ir * ir * hx * hx)
         hym2 = 1.d0 / (jr * jr * hy * hy)
         hzm2 = 1.d0 / (kr * kr * hz * hz)
         do 50 j = regl1, regh1
            do 50 i = regl0, regh0
 50            res(i*ir,j*jr,k*kr) =
     @           src(i*ir,j*jr,k*kr) - fac0 *
     @             (hxm2 *
     @               ((sigmac(i-1,j-1,ks) + sigmac(i-1,j,ks)) *
     @                 (cdst(i-1,j,k) - cdst(i,j,k)) +
     @                (sigmac(i,j-1,ks) + sigmac(i,j,ks)) *
     @                 (cdst(i+1,j,k) - cdst(i,j,k))) +
     @              hym2 *
     @               ((sigmac(i-1,j-1,ks) + sigmac(i,j-1,ks)) *
     @                 (cdst(i,j-1,k) - cdst(i,j,k)) +
     @                (sigmac(i-1,j,ks) + sigmac(i,j,ks)) *
     @                 (cdst(i,j+1,k) - cdst(i,j,k))) +
     @              hzm2 *
     @               ((sigmac(i-1,j-1,ks) + sigmac(i-1,j,ks) +
     @                 sigmac(i,j-1,ks)   + sigmac(i,j,ks)) *
     @                 (cdst(i,j,k-idir) - cdst(i,j,k))))
         fac0 = fac0 / (ir * jr * kr * ir * jr)
         hxm2 = ir * ir * hxm2
         hym2 = jr * jr * hym2
         hzm2 = kr * kr * hzm2
         k = k * kr
         if (idir .eq. 1) then
            ks = k
         else
            ks = k - 1
         endif
         do 60 n = 0, jr-1
            fac2 = (jr-n) * fac0
            if (n .eq. 0) fac2 = 0.5d0 * fac2
            do 60 m = 0, ir-1
               fac1 = (ir-m) * fac2
               if (m .eq. 0) fac1 = 0.5d0 * fac1
               do 60 j = jr*regl1, jr*regh1, jr
                  do 60 i = ir*regl0, ir*regh0, ir
                     tmp = hxm2 *
     @                 ((sigmaf(i-m-1,j-n-1,ks) + sigmaf(i-m-1,j-n,ks)) *
     @                   (fdst(i-m-1,j-n,k) - fdst(i-m,j-n,k)) +
     @                  (sigmaf(i-m,j-n-1,ks) + sigmaf(i-m,j-n,ks)) *
     @                   (fdst(i-m+1,j-n,k) - fdst(i-m,j-n,k)) +
     @                  (sigmaf(i-m-1,j+n-1,ks) + sigmaf(i-m-1,j+n,ks)) *
     @                   (fdst(i-m-1,j+n,k) - fdst(i-m,j+n,k)) +
     @                  (sigmaf(i-m,j+n-1,ks) + sigmaf(i-m,j+n,ks)) *
     @                   (fdst(i-m+1,j+n,k) - fdst(i-m,j+n,k)) +
     @                  (sigmaf(i+m-1,j-n-1,ks) + sigmaf(i+m-1,j-n,ks)) *
     @                   (fdst(i+m-1,j-n,k) - fdst(i+m,j-n,k)) +
     @                  (sigmaf(i+m,j-n-1,ks) + sigmaf(i+m,j-n,ks)) *
     @                   (fdst(i+m+1,j-n,k) - fdst(i+m,j-n,k)) +
     @                  (sigmaf(i+m-1,j+n-1,ks) + sigmaf(i+m-1,j+n,ks)) *
     @                   (fdst(i+m-1,j+n,k) - fdst(i+m,j+n,k)) +
     @                  (sigmaf(i+m,j+n-1,ks) + sigmaf(i+m,j+n,ks)) *
     @                   (fdst(i+m+1,j+n,k) - fdst(i+m,j+n,k)))
                     tmp = tmp + hym2 *
     @                 ((sigmaf(i-m-1,j-n-1,ks) + sigmaf(i-m,j-n-1,ks)) *
     @                   (fdst(i-m,j-n-1,k) - fdst(i-m,j-n,k)) +
     @                  (sigmaf(i-m-1,j-n,ks) + sigmaf(i-m,j-n,ks)) *
     @                   (fdst(i-m,j-n+1,k) - fdst(i-m,j-n,k)) +
     @                  (sigmaf(i-m-1,j+n-1,ks) + sigmaf(i-m,j+n-1,ks)) *
     @                   (fdst(i-m,j+n-1,k) - fdst(i-m,j+n,k)) +
     @                  (sigmaf(i-m-1,j+n,ks) + sigmaf(i-m,j+n,ks)) *
     @                   (fdst(i-m,j+n+1,k) - fdst(i-m,j+n,k)) +
     @                  (sigmaf(i+m-1,j-n-1,ks) + sigmaf(i+m,j-n-1,ks)) *
     @                   (fdst(i+m,j-n-1,k) - fdst(i+m,j-n,k)) +
     @                  (sigmaf(i+m-1,j-n,ks) + sigmaf(i+m,j-n,ks)) *
     @                   (fdst(i+m,j-n+1,k) - fdst(i+m,j-n,k)) +
     @                  (sigmaf(i+m-1,j+n-1,ks) + sigmaf(i+m,j+n-1,ks)) *
     @                   (fdst(i+m,j+n-1,k) - fdst(i+m,j+n,k)) +
     @                  (sigmaf(i+m-1,j+n,ks) + sigmaf(i+m,j+n,ks)) *
     @                   (fdst(i+m,j+n+1,k) - fdst(i+m,j+n,k)))
 60            res(i,j,k) = res(i,j,k) - fac1 * (tmp + hzm2 *
     @                 ((sigmaf(i-m-1,j-n-1,ks) + sigmaf(i-m-1,j-n,ks) +
     @                   sigmaf(i-m,j-n-1,ks) + sigmaf(i-m,j-n,ks)) *
     @                   (fdst(i-m,j-n,k+idir) - fdst(i-m,j-n,k)) +
     @                  (sigmaf(i-m-1,j+n-1,ks) + sigmaf(i-m-1,j+n,ks) +
     @                   sigmaf(i-m,j+n-1,ks) + sigmaf(i-m,j+n,ks)) *
     @                   (fdst(i-m,j+n,k+idir) - fdst(i-m,j+n,k)) +
     @                  (sigmaf(i+m-1,j-n-1,ks) + sigmaf(i+m-1,j-n,ks) +
     @                   sigmaf(i+m,j-n-1,ks) + sigmaf(i+m,j-n,ks)) *
     @                   (fdst(i+m,j-n,k+idir) - fdst(i+m,j-n,k)) +
     @                  (sigmaf(i+m-1,j+n-1,ks) + sigmaf(i+m-1,j+n,ks) +
     @                   sigmaf(i+m,j+n-1,ks) + sigmaf(i+m,j+n,ks)) *
     @                   (fdst(i+m,j+n,k+idir) - fdst(i+m,j+n,k))))
      endif
      end

c end of variable density stencils

c Note---assumes fdst linearly interpolated from cdst along face
      subroutine hgeres(res, dims(res),
     @ src,  dims(src),
     @ fdst, dims(fdst),
     @ cdst, dims(cdst),
     @ sigmaf, dims(sf),
     @ sigmac, dims(sc),
     @ dims(reg),
     @ hx, hy, hz, ir, jr, kr, ivect, ga)
      integer dims(res)
      integer dims(src)
      integer dims(fdst)
      integer dims(cdst)
      integer dims(sf)
      integer dims(sc)
      integer dims(reg)
      real*8 res(dimdec(res))
      real*8 src(dimdec(src))
      real*8 fdst(dimdec(fdst))
      real*8 cdst(dimdec(cdst))
      real*8 sigmaf(dimdec(sf))
      real*8 sigmac(dimdec(sc))
      real*8 hx, hy, hz
      integer ir, jr, kr, ivect(0:2), ga(0:1,0:1,0:1)
      real*8 r3, hxm2, hym2, hzm2, hxm2c, hym2c, hzm2c
      real*8 center, cfac, ffac, fac0, fac1, fac, tmp
      integer ic, jc, kc, if, jf, kf, ii, ji, ki, idir, jdir, kdir
      integer l, m, n
      r3 = ir * jr * kr
      hxm2c = 1.d0 / (ir * ir * hx * hx)
      hym2c = 1.d0 / (jr * jr * hy * hy)
      hzm2c = 1.d0 / (kr * kr * hz * hz)
      hxm2 = ir * ir * hxm2c
      hym2 = jr * jr * hym2c
      hzm2 = kr * kr * hzm2c
      ic = regl0
      jc = regl1
      kc = regl2
      if = ic * ir
      jf = jc * jr
      kf = kc * kr
      center = 0.d0
      if (ivect(0) .eq. 0) then
         do 10 if = ir*regl0, ir*regh0, ir
 10         res(if,jf,kf) = 0.d0
c quadrants
c each quadrant is two octants and their share of the two central edges
         fac0 = 1.d0 / ir
         ffac = ir
         cfac = r3
         do 40 ki = 0, 1
            kdir = 2 * ki - 1
            do 40 ji = 0, 1
               jdir = 2 * ji - 1
               if (ga(0,ji,ki) .eq. 1) then
                  center = center + ffac
                  do 20 m = 0, ir-1
                     fac = (ir-m) * fac0
                     if (m .eq. 0) fac = 0.5d0 * fac
                     do 20 if = ir*regl0, ir*regh0, ir
                        tmp = hxm2 *
     @                    (sigmaf(if-m-1,jf+ji-1,kf+ki-1) *
     @                      (fdst(if-m-1,jf,kf) - fdst(if-m,jf,kf)) +
     @                     sigmaf(if-m,jf+ji-1,kf+ki-1) *
     @                      (fdst(if-m+1,jf,kf) - fdst(if-m,jf,kf)) +
     @                     sigmaf(if+m-1,jf+ji-1,kf+ki-1) *
     @                      (fdst(if+m-1,jf,kf) - fdst(if+m,jf,kf)) +
     @                     sigmaf(if+m,jf+ji-1,kf+ki-1) *
     @                      (fdst(if+m+1,jf,kf) - fdst(if+m,jf,kf)))
                        tmp = tmp + hym2 *
     @                    ((sigmaf(if-m-1,jf+ji-1,kf+ki-1) +
     @                      sigmaf(if-m,jf+ji-1,kf+ki-1)) *
     @                      (fdst(if-m,jf+jdir,kf) - fdst(if-m,jf,kf)) +
     @                     (sigmaf(if+m-1,jf+ji-1,kf+ki-1) +
     @                      sigmaf(if+m,jf+ji-1,kf+ki-1)) *
     @                      (fdst(if+m,jf+jdir,kf) - fdst(if+m,jf,kf)))
 20                     res(if,jf,kf) = res(if,jf,kf) + fac * (tmp + hzm2 *
     @                    ((sigmaf(if-m-1,jf+ji-1,kf+ki-1) +
     @                      sigmaf(if-m,jf+ji-1,kf+ki-1)) *
     @                      (fdst(if-m,jf,kf+kdir) - fdst(if-m,jf,kf)) +
     @                     (sigmaf(if+m-1,jf+ji-1,kf+ki-1) +
     @                      sigmaf(if+m,jf+ji-1,kf+ki-1)) *
     @                      (fdst(if+m,jf,kf+kdir) - fdst(if+m,jf,kf))))
               else
                  center = center + cfac
                  do 30 ic = regl0, regh0
                     if = ic * ir
 30                  res(if,jf,kf) = res(if,jf,kf) + r3 *
     @                (sigmac(ic-1,jc+ji-1,kc+ki-1) *
     @                   (hxm2c * (cdst(ic-1,jc,kc) - cdst(ic,jc,kc)) +
     @                    hym2c * (cdst(ic,jc+jdir,kc) - cdst(ic,jc,kc)) +
     @                    hzm2c * (cdst(ic,jc,kc+kdir) - cdst(ic,jc,kc))) +
     @                 sigmac(ic,jc+ji-1,kc+ki-1) *
     @                   (hxm2c * (cdst(ic+1,jc,kc) - cdst(ic,jc,kc)) +
     @                    hym2c * (cdst(ic,jc+jdir,kc) - cdst(ic,jc,kc)) +
     @                    hzm2c * (cdst(ic,jc,kc+kdir) - cdst(ic,jc,kc))))
               endif
 40            continue
c faces
c each face is two faces and two sides of an edge
         do 70 ki = 0, 1
            kdir = 2 * ki - 1
            do 70 ji = 0, 1
               jdir = 2 * ji - 1
               if (ga(0,ji,ki) - ga(0,ji,1-ki) .eq. 1) then
                  fac0 = 1.d0 / (ir * jr)
                  ffac = ir * (jr - 1)
                  center = center + ffac
                  do 50 n = jdir, jdir*(jr-1), jdir
                     fac1 = (jr-abs(n)) * fac0
                     do 50 m = 0, ir-1
                        fac = (ir-m) * fac1
                        if (m .eq. 0) fac = 0.5d0 * fac
                        do 50 if = ir*regl0, ir*regh0, ir
                           tmp = hxm2 *
     @       ((sigmaf(if-m-1,jf+n-1,kf+ki-1) + sigmaf(if-m-1,jf+n,kf+ki-1)) *
     @         (fdst(if-m-1,jf+n,kf) - fdst(if-m,jf+n,kf)) +
     @        (sigmaf(if-m,jf+n-1,kf+ki-1) + sigmaf(if-m,jf+n,kf+ki-1)) *
     @         (fdst(if-m+1,jf+n,kf) - fdst(if-m,jf+n,kf)) +
     @        (sigmaf(if+m-1,jf+n-1,kf+ki-1) + sigmaf(if+m-1,jf+n,kf+ki-1)) *
     @         (fdst(if+m-1,jf+n,kf) - fdst(if+m,jf+n,kf)) +
     @        (sigmaf(if+m,jf+n-1,kf+ki-1) + sigmaf(if+m,jf+n,kf+ki-1)) *
     @         (fdst(if+m+1,jf+n,kf) - fdst(if+m,jf+n,kf)))
                           tmp = tmp + hym2 *
     @       ((sigmaf(if-m-1,jf+n-1,kf+ki-1) + sigmaf(if-m,jf+n-1,kf+ki-1)) *
     @         (fdst(if-m,jf+n-1,kf) - fdst(if-m,jf+n,kf)) +
     @        (sigmaf(if-m-1,jf+n,kf+ki-1) + sigmaf(if-m,jf+n,kf+ki-1)) *
     @         (fdst(if-m,jf+n+1,kf) - fdst(if-m,jf+n,kf)) +
     @        (sigmaf(if+m-1,jf+n-1,kf+ki-1) + sigmaf(if+m,jf+n-1,kf+ki-1)) *
     @         (fdst(if+m,jf+n-1,kf) - fdst(if+m,jf+n,kf)) +
     @        (sigmaf(if+m-1,jf+n,kf+ki-1) + sigmaf(if+m,jf+n,kf+ki-1)) *
     @         (fdst(if+m,jf+n+1,kf) - fdst(if+m,jf+n,kf)))
 50                        res(if,jf,kf) = res(if,jf,kf) + fac * (tmp + hzm2 *
     @       ((sigmaf(if-m-1,jf+n-1,kf+ki-1) + sigmaf(if-m-1,jf+n,kf+ki-1) +
     @         sigmaf(if-m,jf+n-1,kf+ki-1) + sigmaf(if-m,jf+n,kf+ki-1)) *
     @         (fdst(if-m,jf+n,kf+kdir) - fdst(if-m,jf+n,kf)) +
     @        (sigmaf(if+m-1,jf+n-1,kf+ki-1) + sigmaf(if+m-1,jf+n,kf+ki-1) +
     @         sigmaf(if+m,jf+n-1,kf+ki-1) + sigmaf(if+m,jf+n,kf+ki-1)) *
     @         (fdst(if+m,jf+n,kf+kdir) - fdst(if+m,jf+n,kf))))
               endif
               if (ga(0,ji,ki) - ga(0,1-ji,ki) .eq. 1) then
                  fac0 = 1.d0 / (ir * kr)
                  ffac = ir * (kr - 1)
                  center = center + ffac
                  do 60 l = kdir, kdir*(kr-1), kdir
                     fac1 = (kr-abs(l)) * fac0
                     do 60 m = 0, ir-1
                        fac = (ir-m) * fac1
                        if (m .eq. 0) fac = 0.5d0 * fac
                        do 60 if = ir*regl0, ir*regh0, ir
                           tmp = hxm2 *
     @       ((sigmaf(if-m-1,jf+ji-1,kf+l-1) + sigmaf(if-m-1,jf+ji-1,kf+l)) *
     @         (fdst(if-m-1,jf,kf+l) - fdst(if-m,jf,kf+l)) +
     @        (sigmaf(if-m,jf+ji-1,kf+l-1) + sigmaf(if-m,jf+ji-1,kf+l)) *
     @         (fdst(if-m+1,jf,kf+l) - fdst(if-m,jf,kf+l)) +
     @        (sigmaf(if+m-1,jf+ji-1,kf+l-1) + sigmaf(if+m-1,jf+ji-1,kf+l)) *
     @         (fdst(if+m-1,jf,kf+l) - fdst(if+m,jf,kf+l)) +
     @        (sigmaf(if+m,jf+ji-1,kf+l-1) + sigmaf(if+m,jf+ji-1,kf+l)) *
     @         (fdst(if+m+1,jf,kf+l) - fdst(if+m,jf,kf+l)))
                           tmp = tmp + hym2 *
     @       ((sigmaf(if-m-1,jf+ji-1,kf+l-1) + sigmaf(if-m-1,jf+ji-1,kf+l) +
     @         sigmaf(if-m,jf+ji-1,kf+l-1) + sigmaf(if-m,jf+ji-1,kf+l)) *
     @         (fdst(if-m,jf+jdir,kf+l) - fdst(if-m,jf,kf+l)) +
     @        (sigmaf(if+m-1,jf+ji-1,kf+l-1) + sigmaf(if+m-1,jf+ji-1,kf+l) +
     @         sigmaf(if+m,jf+ji-1,kf+l-1) + sigmaf(if+m,jf+ji-1,kf+l)) *
     @         (fdst(if+m,jf+jdir,kf+l) - fdst(if+m,jf,kf+l)))
 60                        res(if,jf,kf) = res(if,jf,kf) + fac * (tmp + hzm2 *
     @       ((sigmaf(if-m-1,jf+ji-1,kf+l-1) + sigmaf(if-m,jf+ji-1,kf+l-1)) *
     @         (fdst(if-m,jf,kf+l-1) - fdst(if-m,jf,kf+l)) +
     @        (sigmaf(if-m-1,jf+ji-1,kf+l) + sigmaf(if-m,jf+ji-1,kf+l)) *
     @         (fdst(if-m,jf,kf+l+1) - fdst(if-m,jf,kf+l)) +
     @        (sigmaf(if+m-1,jf+ji-1,kf+l-1) + sigmaf(if+m,jf+ji-1,kf+l-1)) *
     @         (fdst(if+m,jf,kf+l-1) - fdst(if+m,jf,kf+l)) +
     @        (sigmaf(if+m-1,jf+ji-1,kf+l) + sigmaf(if+m,jf+ji-1,kf+l)) *
     @         (fdst(if+m,jf,kf+l+1) - fdst(if+m,jf,kf+l))))
               endif
 70            continue
c weighting
         do 80 if = ir*regl0, ir*regh0, ir
 80         res(if,jf,kf) = src(if,jf,kf) - res(if,jf,kf) / center
      else if (ivect(1) .eq. 0) then
         do 110 jf = jr*regl1, jr*regh1, jr
 110        res(if,jf,kf) = 0.d0
c quadrants
c each quadrant is two octants and their share of the two central edges
         fac0 = 1.d0 / jr
         ffac = jr
         cfac = r3
         do 140 ki = 0, 1
            kdir = 2 * ki - 1
            do 140 ii = 0, 1
               idir = 2 * ii - 1
               if (ga(ii,0,ki) .eq. 1) then
                  center = center + ffac
                  do 120 n = 0, jr-1
                     fac = (jr-n) * fac0
                     if (n .eq. 0) fac = 0.5d0 * fac
                     do 120 jf = jr*regl1, jr*regh1, jr
                        tmp = hxm2 *
     @                    ((sigmaf(if+ii-1,jf-n-1,kf+ki-1) +
     @                      sigmaf(if+ii-1,jf-n,kf+ki-1)) *
     @                      (fdst(if+idir,jf-n,kf) - fdst(if,jf-n,kf)) +
     @                     (sigmaf(if+ii-1,jf+n-1,kf+ki-1) +
     @                      sigmaf(if+ii-1,jf+n,kf+ki-1)) *
     @                      (fdst(if+idir,jf+n,kf) - fdst(if,jf+n,kf)))
                        tmp = tmp + hym2 *
     @                    (sigmaf(if+ii-1,jf-n-1,kf+ki-1) *
     @                      (fdst(if,jf-n-1,kf) - fdst(if,jf-n,kf)) +
     @                     sigmaf(if+ii-1,jf-n,kf+ki-1) *
     @                      (fdst(if,jf-n+1,kf) - fdst(if,jf-n,kf)) +
     @                     sigmaf(if+ii-1,jf+n-1,kf+ki-1) *
     @                      (fdst(if,jf+n-1,kf) - fdst(if,jf+n,kf)) +
     @                     sigmaf(if+ii-1,jf+n,kf+ki-1) *
     @                      (fdst(if,jf+n+1,kf) - fdst(if,jf+n,kf)))
 120                    res(if,jf,kf) = res(if,jf,kf) + fac * (tmp + hzm2 *
     @                    ((sigmaf(if+ii-1,jf-n-1,kf+ki-1) +
     @                      sigmaf(if+ii-1,jf-n,kf+ki-1)) *
     @                      (fdst(if,jf-n,kf+kdir) - fdst(if,jf-n,kf)) +
     @                     (sigmaf(if+ii-1,jf+n-1,kf+ki-1) +
     @                      sigmaf(if+ii-1,jf+n,kf+ki-1)) *
     @                      (fdst(if,jf+n,kf+kdir) - fdst(if,jf+n,kf))))
               else
                  center = center + cfac
                  do 130 jc = regl1, regh1
                     jf = jc * jr
 130                 res(if,jf,kf) = res(if,jf,kf) + r3 *
     @                (sigmac(ic+ii-1,jc-1,kc+ki-1) *
     @                   (hxm2c * (cdst(ic+idir,jc,kc) - cdst(ic,jc,kc)) +
     @                    hym2c * (cdst(ic,jc-1,kc) - cdst(ic,jc,kc)) +
     @                    hzm2c * (cdst(ic,jc,kc+kdir) - cdst(ic,jc,kc))) +
     @                 sigmac(ic+ii-1,jc,kc+ki-1) *
     @                   (hxm2c * (cdst(ic+idir,jc,kc) - cdst(ic,jc,kc)) +
     @                    hym2c * (cdst(ic,jc+1,kc) - cdst(ic,jc,kc)) +
     @                    hzm2c * (cdst(ic,jc,kc+kdir) - cdst(ic,jc,kc))))
               endif
 140           continue
c faces
c each face is two faces and two sides of an edge
         do 170 ki = 0, 1
            kdir = 2 * ki - 1
            do 170 ii = 0, 1
               idir = 2 * ii - 1
               if (ga(ii,0,ki) - ga(ii,0,1-ki) .eq. 1) then
                  fac0 = 1.d0 / (ir * jr)
                  ffac = jr * (ir - 1)
                  center = center + ffac
                  do 150 n = 0, jr-1
                     fac1 = (jr-n) * fac0
                     if (n .eq. 0) fac1 = 0.5d0 * fac1
                     do 150 m = idir, idir*(ir-1), idir
                        fac = (ir-abs(m)) * fac1
                        do 150 jf = jr*regl1, jr*regh1, jr
                           tmp = hxm2 *
     @       ((sigmaf(if+m-1,jf-n-1,kf+ki-1) + sigmaf(if+m-1,jf-n,kf+ki-1)) *
     @         (fdst(if+m-1,jf-n,kf) - fdst(if+m,jf-n,kf)) +
     @        (sigmaf(if+m,jf-n-1,kf+ki-1) + sigmaf(if+m,jf-n,kf+ki-1)) *
     @         (fdst(if+m+1,jf-n,kf) - fdst(if+m,jf-n,kf)) +
     @        (sigmaf(if+m-1,jf+n-1,kf+ki-1) + sigmaf(if+m-1,jf+n,kf+ki-1)) *
     @         (fdst(if+m-1,jf+n,kf) - fdst(if+m,jf+n,kf)) +
     @        (sigmaf(if+m,jf+n-1,kf+ki-1) + sigmaf(if+m,jf+n,kf+ki-1)) *
     @         (fdst(if+m+1,jf+n,kf) - fdst(if+m,jf+n,kf)))
                           tmp = tmp + hym2 *
     @       ((sigmaf(if+m-1,jf-n-1,kf+ki-1) + sigmaf(if+m,jf-n-1,kf+ki-1)) *
     @         (fdst(if+m,jf-n-1,kf) - fdst(if+m,jf-n,kf)) +
     @        (sigmaf(if+m-1,jf-n,kf+ki-1) + sigmaf(if+m,jf-n,kf+ki-1)) *
     @         (fdst(if+m,jf-n+1,kf) - fdst(if+m,jf-n,kf)) +
     @        (sigmaf(if+m-1,jf+n-1,kf+ki-1) + sigmaf(if+m,jf+n-1,kf+ki-1)) *
     @         (fdst(if+m,jf+n-1,kf) - fdst(if+m,jf+n,kf)) +
     @        (sigmaf(if+m-1,jf+n,kf+ki-1) + sigmaf(if+m,jf+n,kf+ki-1)) *
     @         (fdst(if+m,jf+n+1,kf) - fdst(if+m,jf+n,kf)))
 150                       res(if,jf,kf) = res(if,jf,kf) + fac * (tmp + hzm2 *
     @       ((sigmaf(if+m-1,jf-n-1,kf+ki-1) + sigmaf(if+m-1,jf-n,kf+ki-1) +
     @         sigmaf(if+m,jf-n-1,kf+ki-1) + sigmaf(if+m,jf-n,kf+ki-1)) *
     @         (fdst(if+m,jf-n,kf+kdir) - fdst(if+m,jf-n,kf)) +
     @        (sigmaf(if+m-1,jf+n-1,kf+ki-1) + sigmaf(if+m-1,jf+n,kf+ki-1) +
     @         sigmaf(if+m,jf+n-1,kf+ki-1) + sigmaf(if+m,jf+n,kf+ki-1)) *
     @         (fdst(if+m,jf+n,kf+kdir) - fdst(if+m,jf+n,kf))))
               endif
               if (ga(ii,0,ki) - ga(1-ii,0,ki) .eq. 1) then
                  fac0 = 1.d0 / (jr * kr)
                  ffac = jr * (kr - 1)
                  center = center + ffac
                  do 160 l = kdir, kdir*(kr-1), kdir
                     fac1 = (kr-abs(l)) * fac0
                     do 160 n = 0, jr-1
                        fac = (jr-n) * fac1
                        if (n .eq. 0) fac = 0.5d0 * fac
                        do 160 jf = jr*regl1, jr*regh1, jr
                           tmp = hxm2 *
     @       ((sigmaf(if+ii-1,jf-n-1,kf+l-1) + sigmaf(if+ii-1,jf-n-1,kf+l) +
     @         sigmaf(if+ii-1,jf-n,kf+l-1) + sigmaf(if+ii-1,jf-n,kf+l)) *
     @         (fdst(if+idir,jf-n,kf+l) - fdst(if,jf-n,kf+l)) +
     @        (sigmaf(if+ii-1,jf+n-1,kf+l-1) + sigmaf(if+ii-1,jf+n-1,kf+l) +
     @         sigmaf(if+ii-1,jf+n,kf+l-1) + sigmaf(if+ii-1,jf+n,kf+l)) *
     @         (fdst(if+idir,jf+n,kf+l) - fdst(if,jf+n,kf+l)))
                           tmp = tmp + hym2 *
     @       ((sigmaf(if+ii-1,jf-n-1,kf+l-1) + sigmaf(if+ii-1,jf-n-1,kf+l)) *
     @         (fdst(if,jf-n-1,kf+l) - fdst(if,jf-n,kf+l)) +
     @        (sigmaf(if+ii-1,jf-n,kf+l-1) + sigmaf(if+ii-1,jf-n,kf+l)) *
     @         (fdst(if,jf-n+1,kf+l) - fdst(if,jf-n,kf+l)) +
     @        (sigmaf(if+ii-1,jf+n-1,kf+l-1) + sigmaf(if+ii-1,jf+n-1,kf+l)) *
     @         (fdst(if,jf+n-1,kf+l) - fdst(if,jf+n,kf+l)) +
     @        (sigmaf(if+ii-1,jf+n,kf+l-1) + sigmaf(if+ii-1,jf+n,kf+l)) *
     @         (fdst(if,jf+n+1,kf+l) - fdst(if,jf+n,kf+l)))
 160                       res(if,jf,kf) = res(if,jf,kf) + fac * (tmp + hzm2 *
     @       ((sigmaf(if+ii-1,jf-n-1,kf+l-1) + sigmaf(if+ii-1,jf-n,kf+l-1)) *
     @         (fdst(if,jf-n,kf+l-1) - fdst(if,jf-n,kf+l)) +
     @        (sigmaf(if+ii-1,jf-n-1,kf+l) + sigmaf(if+ii-1,jf-n,kf+l)) *
     @         (fdst(if,jf-n,kf+l+1) - fdst(if,jf-n,kf+l)) +
     @        (sigmaf(if+ii-1,jf+n-1,kf+l-1) + sigmaf(if+ii-1,jf+n,kf+l-1)) *
     @         (fdst(if,jf+n,kf+l-1) - fdst(if,jf+n,kf+l)) +
     @        (sigmaf(if+ii-1,jf+n-1,kf+l) + sigmaf(if+ii-1,jf+n,kf+l)) *
     @         (fdst(if,jf+n,kf+l+1) - fdst(if,jf+n,kf+l))))
               endif
 170           continue
c weighting
         do 180 jf = jr*regl1, jr*regh1, jr
 180        res(if,jf,kf) = src(if,jf,kf) - res(if,jf,kf) / center
      else
         do 210 kf = kr*regl2, kr*regh2, kr
 210        res(if,jf,kf) = 0.d0
c quadrants
c each quadrant is two octants and their share of the two central edges
         fac0 = 1.d0 / kr
         ffac = kr
         cfac = r3
         do 240 ji = 0, 1
            jdir = 2 * ji - 1
            do 240 ii = 0, 1
               idir = 2 * ii - 1
               if (ga(ii,ji,0) .eq. 1) then
                  center = center + ffac
                  do 220 l = 0, kr-1
                     fac = (kr-l) * fac0
                     if (l .eq. 0) fac = 0.5d0 * fac
                     do 220 kf = kr*regl2, kr*regh2, kr
                        tmp = hxm2 *
     @                    ((sigmaf(if+ii-1,jf+ji-1,kf-l-1) +
     @                      sigmaf(if+ii-1,jf+ji-1,kf-l)) *
     @                      (fdst(if+idir,jf,kf-l) - fdst(if,jf,kf-l)) +
     @                     (sigmaf(if+ii-1,jf+ji-1,kf+l-1) +
     @                      sigmaf(if+ii-1,jf+ji-1,kf+l)) *
     @                      (fdst(if+idir,jf,kf+l) - fdst(if,jf,kf+l)))
                        tmp = tmp + hym2 *
     @                    ((sigmaf(if+ii-1,jf+ji-1,kf-l-1) +
     @                      sigmaf(if+ii-1,jf+ji-1,kf-l)) *
     @                      (fdst(if,jf+jdir,kf-l) - fdst(if,jf,kf-l)) +
     @                     (sigmaf(if+ii-1,jf+ji-1,kf+l-1) +
     @                      sigmaf(if+ii-1,jf+ji-1,kf+l)) *
     @                      (fdst(if,jf+jdir,kf+l) - fdst(if,jf,kf+l)))
 220                    res(if,jf,kf) = res(if,jf,kf) + fac * (tmp + hzm2 *
     @                    (sigmaf(if+ii-1,jf+ji-1,kf-l-1) *
     @                      (fdst(if,jf,kf-l-1) - fdst(if,jf,kf-l)) +
     @                     sigmaf(if+ii-1,jf+ji-1,kf-l) *
     @                      (fdst(if,jf,kf-l+1) - fdst(if,jf,kf-l)) +
     @                     sigmaf(if+ii-1,jf+ji-1,kf+l-1) *
     @                      (fdst(if,jf,kf+l-1) - fdst(if,jf,kf+l)) +
     @                     sigmaf(if+ii-1,jf+ji-1,kf+l) *
     @                      (fdst(if,jf,kf+l+1) - fdst(if,jf,kf+l))))
               else
                  center = center + cfac
                  do 230 kc = regl2, regh2
                     kf = kc * kr
 230                 res(if,jf,kf) = res(if,jf,kf) + r3 *
     @                (sigmac(ic+ii-1,jc+ji-1,kc-1) *
     @                   (hxm2c * (cdst(ic+idir,jc,kc) - cdst(ic,jc,kc)) +
     @                    hym2c * (cdst(ic,jc+jdir,kc) - cdst(ic,jc,kc)) +
     @                    hzm2c * (cdst(ic,jc,kc-1) - cdst(ic,jc,kc))) +
     @                 sigmac(ic+ii-1,jc+ji-1,kc) *
     @                   (hxm2c * (cdst(ic+idir,jc,kc) - cdst(ic,jc,kc)) +
     @                    hym2c * (cdst(ic,jc+jdir,kc) - cdst(ic,jc,kc)) +
     @                    hzm2c * (cdst(ic,jc,kc+1) - cdst(ic,jc,kc))))
               endif
 240           continue
c faces
c each face is two faces and two sides of an edge
         do 270 ji = 0, 1
            jdir = 2 * ji - 1
            do 270 ii = 0, 1
               idir = 2 * ii - 1
               if (ga(ii,ji,0) - ga(ii,1-ji,0) .eq. 1) then
                  fac0 = 1.d0 / (ir * kr)
                  ffac = kr * (ir - 1)
                  center = center + ffac
                  do 250 l = 0, kr-1
                     fac1 = (kr-l) * fac0
                     if (l .eq. 0) fac1 = 0.5d0 * fac1
                     do 250 m = idir, idir*(ir-1), idir
                        fac = (ir-abs(m)) * fac1
                        do 250 kf = kr*regl2, kr*regh2, kr
                           tmp = hxm2 *
     @       ((sigmaf(if+m-1,jf+ji-1,kf-l-1) + sigmaf(if+m-1,jf+ji-1,kf-l)) *
     @         (fdst(if+m-1,jf,kf-l) - fdst(if+m,jf,kf-l)) +
     @        (sigmaf(if+m,jf+ji-1,kf-l-1) + sigmaf(if+m,jf+ji-1,kf-l)) *
     @         (fdst(if+m+1,jf,kf-l) - fdst(if+m,jf,kf-l)) +
     @        (sigmaf(if+m-1,jf+ji-1,kf+l-1) + sigmaf(if+m-1,jf+ji-1,kf+l)) *
     @         (fdst(if+m-1,jf,kf+l) - fdst(if+m,jf,kf+l)) +
     @        (sigmaf(if+m,jf+ji-1,kf+l-1) + sigmaf(if+m,jf+ji-1,kf+l)) *
     @         (fdst(if+m+1,jf,kf+l) - fdst(if+m,jf,kf+l)))
                           tmp = tmp + hym2 *
     @       ((sigmaf(if+m-1,jf+ji-1,kf-l-1) + sigmaf(if+m-1,jf+ji-1,kf-l) +
     @         sigmaf(if+m,jf+ji-1,kf-l-1) + sigmaf(if+m,jf+ji-1,kf-l)) *
     @         (fdst(if+m,jf+jdir,kf-l) - fdst(if+m,jf,kf-l)) +
     @        (sigmaf(if+m-1,jf+ji-1,kf+l-1) + sigmaf(if+m-1,jf+ji-1,kf+l) +
     @         sigmaf(if+m,jf+ji-1,kf+l-1) + sigmaf(if+m,jf+ji-1,kf+l)) *
     @         (fdst(if+m,jf+jdir,kf+l) - fdst(if+m,jf,kf+l)))
 250                       res(if,jf,kf) = res(if,jf,kf) + fac * (tmp + hzm2 *
     @       ((sigmaf(if+m-1,jf+ji-1,kf-l-1) + sigmaf(if+m,jf+ji-1,kf-l-1)) *
     @         (fdst(if+m,jf,kf-l-1) - fdst(if+m,jf,kf-l)) +
     @        (sigmaf(if+m-1,jf+ji-1,kf-l) + sigmaf(if+m,jf+ji-1,kf-l)) *
     @         (fdst(if+m,jf,kf-l+1) - fdst(if+m,jf,kf-l)) +
     @        (sigmaf(if+m-1,jf+ji-1,kf+l-1) + sigmaf(if+m,jf+ji-1,kf+l-1)) *
     @         (fdst(if+m,jf,kf+l-1) - fdst(if+m,jf,kf+l)) +
     @        (sigmaf(if+m-1,jf+ji-1,kf+l) + sigmaf(if+m,jf+ji-1,kf+l)) *
     @         (fdst(if+m,jf,kf+l+1) - fdst(if+m,jf,kf+l))))
               endif
               if (ga(ii,ji,0) - ga(1-ii,ji,0) .eq. 1) then
                  fac0 = 1.d0 / (jr * kr)
                  ffac = kr * (jr - 1)
                  center = center + ffac
                  do 260 l = 0, kr-1
                     fac1 = (kr-l) * fac0
                     if (l .eq. 0) fac1 = 0.5d0 * fac1
                     do 260 n = jdir, jdir*(jr-1), jdir
                        fac = (jr-abs(n)) * fac1
                        do 260 kf = kr*regl2, kr*regh2, kr
                           tmp = hxm2 *
     @       ((sigmaf(if+ii-1,jf+n-1,kf-l-1) + sigmaf(if+ii-1,jf+n-1,kf-l) +
     @         sigmaf(if+ii-1,jf+n,kf-l-1) + sigmaf(if+ii-1,jf+n,kf-l)) *
     @         (fdst(if+idir,jf+n,kf-l) - fdst(if,jf+n,kf-l)) +
     @        (sigmaf(if+ii-1,jf+n-1,kf+l-1) + sigmaf(if+ii-1,jf+n-1,kf+l) +
     @         sigmaf(if+ii-1,jf+n,kf+l-1) + sigmaf(if+ii-1,jf+n,kf+l)) *
     @         (fdst(if+idir,jf+n,kf+l) - fdst(if,jf+n,kf+l)))
                           tmp = tmp + hym2 *
     @       ((sigmaf(if+ii-1,jf+n-1,kf-l-1) + sigmaf(if+ii-1,jf+n-1,kf-l)) *
     @         (fdst(if,jf+n-1,kf-l) - fdst(if,jf+n,kf-l)) +
     @        (sigmaf(if+ii-1,jf+n,kf-l-1) + sigmaf(if+ii-1,jf+n,kf-l)) *
     @         (fdst(if,jf+n+1,kf-l) - fdst(if,jf+n,kf-l)) +
     @        (sigmaf(if+ii-1,jf+n-1,kf+l-1) + sigmaf(if+ii-1,jf+n-1,kf+l)) *
     @         (fdst(if,jf+n-1,kf+l) - fdst(if,jf+n,kf+l)) +
     @        (sigmaf(if+ii-1,jf+n,kf+l-1) + sigmaf(if+ii-1,jf+n,kf+l)) *
     @         (fdst(if,jf+n+1,kf+l) - fdst(if,jf+n,kf+l)))
 260                       res(if,jf,kf) = res(if,jf,kf) + fac * (tmp + hzm2 *
     @       ((sigmaf(if+ii-1,jf+n-1,kf-l-1) + sigmaf(if+ii-1,jf+n,kf-l-1)) *
     @         (fdst(if,jf+n,kf-l-1) - fdst(if,jf+n,kf-l)) +
     @        (sigmaf(if+ii-1,jf+n-1,kf-l) + sigmaf(if+ii-1,jf+n,kf-l)) *
     @         (fdst(if,jf+n,kf-l+1) - fdst(if,jf+n,kf-l)) +
     @        (sigmaf(if+ii-1,jf+n-1,kf+l-1) + sigmaf(if+ii-1,jf+n,kf+l-1)) *
     @         (fdst(if,jf+n,kf+l-1) - fdst(if,jf+n,kf+l)) +
     @        (sigmaf(if+ii-1,jf+n-1,kf+l) + sigmaf(if+ii-1,jf+n,kf+l)) *
     @         (fdst(if,jf+n,kf+l+1) - fdst(if,jf+n,kf+l))))
               endif
 270           continue
c weighting
         do 280 kf = kr*regl2, kr*regh2, kr
 280        res(if,jf,kf) = src(if,jf,kf) - res(if,jf,kf) / center
      endif
      end

c Note---assumes fdst linearly interpolated from cdst along face
      subroutine hgcres(res, dims(res),
     @ src,  dims(src),
     @ fdst, dims(fdst),
     @ cdst, dims(cdst),
     @ sigmaf, dims(sf),
     @ sigmac, dims(sc),
     @ dims(reg),
     @ hx, hy, hz, ir, jr, kr, ga)
      integer dims(res)
      integer dims(src)
      integer dims(fdst)
      integer dims(cdst)
      integer dims(sf)
      integer dims(sc)
      integer dims(reg)
      real*8 res(dimdec(res))
      real*8 src(dimdec(src))
      real*8 fdst(dimdec(fdst))
      real*8 cdst(dimdec(cdst))
      real*8 sigmaf(dimdec(sf))
      real*8 sigmac(dimdec(sc))
      real*8 hx, hy, hz
      integer ir, jr, kr, ga(0:1,0:1,0:1)
      real*8 r3, hxm2, hym2, hzm2, hxm2c, hym2c, hzm2c
      real*8 sum, center, cfac, ffac, fac1, fac2, fac
      integer ic, jc, kc, if, jf, kf, ii, ji, ki, idir, jdir, kdir
      integer l, m, n
      r3 = ir * jr * kr
      hxm2c = 1.d0 / (ir * ir * hx * hx)
      hym2c = 1.d0 / (jr * jr * hy * hy)
      hzm2c = 1.d0 / (kr * kr * hz * hz)
      hxm2 = ir * ir * hxm2c
      hym2 = jr * jr * hym2c
      hzm2 = kr * kr * hzm2c
      ic = regl0
      jc = regl1
      kc = regl2
      if = ic * ir
      jf = jc * jr
      kf = kc * kr
      sum = 0.d0
      center = 0.d0
c octants
      fac = 1.d0
      ffac = 0.5d0
      cfac = 0.5d0 * r3
      do 10 ki = 0, 1
         kdir = 2 * ki - 1
         do 10 ji = 0, 1
            jdir = 2 * ji - 1
            do 10 ii = 0, 1
               idir = 2 * ii - 1
               if (ga(ii,ji,ki) .eq. 1) then
                  center = center + ffac
                  sum = sum + fac * sigmaf(if+ii-1,jf+ji-1,kf+ki-1) *
     @              (hxm2 * (fdst(if+idir,jf,kf) - fdst(if,jf,kf)) +
     @               hym2 * (fdst(if,jf+jdir,kf) - fdst(if,jf,kf)) +
     @               hzm2 * (fdst(if,jf,kf+kdir) - fdst(if,jf,kf)))
               else
                  center = center + cfac
                  sum = sum + r3 * sigmac(ic+ii-1,jc+ji-1,kc+ki-1) *
     @              (hxm2c * (cdst(ic+idir,jc,kc) - cdst(ic,jc,kc)) +
     @               hym2c * (cdst(ic,jc+jdir,kc) - cdst(ic,jc,kc)) +
     @               hzm2c * (cdst(ic,jc,kc+kdir) - cdst(ic,jc,kc)))
               endif
 10            continue
c faces
      do 50 ki = 0, 1
         kdir = 2 * ki - 1
         do 50 ji = 0, 1
            jdir = 2 * ji - 1
            do 50 ii = 0, 1
               idir = 2 * ii - 1
               if (ga(ii,ji,ki) - ga(ii,ji,1-ki) .eq. 1) then
                  fac2 = 1.d0 / (ir * jr)
                  ffac = 0.5d0 * (ir-1) * (jr-1)
                  center = center + ffac
                  do 20 n = jdir, jdir*(jr-1), jdir
                     fac1 = (jr-abs(n)) * fac2
                     do 20 m = idir, idir*(ir-1), idir
                        fac = (ir-abs(m)) * fac1
 20                     sum = sum + fac *
     @                    (hxm2 *
     @       ((sigmaf(if+m-1,jf+n-1,kf+ki-1) + sigmaf(if+m-1,jf+n,kf+ki-1)) *
     @         (fdst(if+m-1,jf+n,kf) - fdst(if+m,jf+n,kf)) +
     @        (sigmaf(if+m,jf+n-1,kf+ki-1) + sigmaf(if+m,jf+n,kf+ki-1)) *
     @         (fdst(if+m+1,jf+n,kf) - fdst(if+m,jf+n,kf))) +
     @                     hym2 *
     @       ((sigmaf(if+m-1,jf+n-1,kf+ki-1) + sigmaf(if+m,jf+n-1,kf+ki-1)) *
     @         (fdst(if+m,jf+n-1,kf) - fdst(if+m,jf+n,kf)) +
     @        (sigmaf(if+m-1,jf+n,kf+ki-1) + sigmaf(if+m,jf+n,kf+ki-1)) *
     @         (fdst(if+m,jf+n+1,kf) - fdst(if+m,jf+n,kf))) +
     @                     hzm2 *
     @       ((sigmaf(if+m-1,jf+n-1,kf+ki-1) + sigmaf(if+m-1,jf+n,kf+ki-1) +
     @         sigmaf(if+m,jf+n-1,kf+ki-1) + sigmaf(if+m,jf+n,kf+ki-1)) *
     @         (fdst(if+m,jf+n,kf+kdir) - fdst(if+m,jf+n,kf))))
               endif
               if (ga(ii,ji,ki) - ga(ii,1-ji,ki) .eq. 1) then
                  fac2 = 1.d0 / (ir * kr)
                  ffac = 0.5d0 * (ir-1) * (kr-1)
                  center = center + ffac
                  do 30 l = kdir, kdir*(kr-1), kdir
                     fac1 = (kr-abs(l)) * fac2
                     do 30 m = idir, idir*(ir-1), idir
                        fac = (ir-abs(m)) * fac1
 30                     sum = sum + fac *
     @                    (hxm2 *
     @       ((sigmaf(if+m-1,jf+ji-1,kf+l-1) + sigmaf(if+m-1,jf+ji-1,kf+l)) *
     @         (fdst(if+m-1,jf,kf+l) - fdst(if+m,jf,kf+l)) +
     @        (sigmaf(if+m,jf+ji-1,kf+l-1) + sigmaf(if+m,jf+ji-1,kf+l)) *
     @         (fdst(if+m+1,jf,kf+l) - fdst(if+m,jf,kf+l))) +
     @                     hym2 *
     @       ((sigmaf(if+m-1,jf+ji-1,kf+l-1) + sigmaf(if+m-1,jf+ji-1,kf+l) +
     @         sigmaf(if+m,jf+ji-1,kf+l-1) + sigmaf(if+m,jf+ji-1,kf+l)) *
     @         (fdst(if+m,jf+jdir,kf+l) - fdst(if+m,jf,kf+l))) +
     @                     hzm2 *
     @       ((sigmaf(if+m-1,jf+ji-1,kf+l-1) + sigmaf(if+m,jf+ji-1,kf+l-1)) *
     @         (fdst(if+m,jf,kf+l-1) - fdst(if+m,jf,kf+l)) +
     @        (sigmaf(if+m-1,jf+ji-1,kf+l) + sigmaf(if+m,jf+ji-1,kf+l)) *
     @         (fdst(if+m,jf,kf+l+1) - fdst(if+m,jf,kf+l))))
               endif
               if (ga(ii,ji,ki) - ga(1-ii,ji,ki) .eq. 1) then
                  fac2 = 1.d0 / (jr * kr)
                  ffac = 0.5d0 * (jr-1) * (kr-1)
                  center = center + ffac
                  do 40 l = kdir, kdir*(kr-1), kdir
                     fac1 = (kr-abs(l)) * fac2
                     do 40 n = jdir, jdir*(jr-1), jdir
                        fac = (jr-abs(n)) * fac1
 40                     sum = sum + fac *
     @                    (hxm2 *
     @       ((sigmaf(if+ii-1,jf+n-1,kf+l-1) + sigmaf(if+ii-1,jf+n-1,kf+l) +
     @         sigmaf(if+ii-1,jf+n,kf+l-1) + sigmaf(if+ii-1,jf+n,kf+l)) *
     @         (fdst(if+idir,jf+n,kf+l) - fdst(if,jf+n,kf+l))) +
     @                     hym2 *
     @       ((sigmaf(if+ii-1,jf+n-1,kf+l-1) + sigmaf(if+ii-1,jf+n-1,kf+l)) *
     @         (fdst(if,jf+n-1,kf+l) - fdst(if,jf+n,kf+l)) +
     @        (sigmaf(if+ii-1,jf+n,kf+l-1) + sigmaf(if+ii-1,jf+n,kf+l)) *
     @         (fdst(if,jf+n+1,kf+l) - fdst(if,jf+n,kf+l))) +
     @                     hzm2 *
     @       ((sigmaf(if+ii-1,jf+n-1,kf+l-1) + sigmaf(if+ii-1,jf+n,kf+l-1)) *
     @         (fdst(if,jf+n,kf+l-1) - fdst(if,jf+n,kf+l)) +
     @        (sigmaf(if+ii-1,jf+n-1,kf+l) + sigmaf(if+ii-1,jf+n,kf+l)) *
     @         (fdst(if,jf+n,kf+l+1) - fdst(if,jf+n,kf+l))))
            endif
 50         continue
c edges
      do 90 ki = 0, 1
         kdir = 2 * ki - 1
         do 90 ji = 0, 1
            jdir = 2 * ji - 1
            do 90 ii = 0, 1
               idir = 2 * ii - 1
               if (ga(ii,ji,ki) -
     @             min(ga(ii,ji,1-ki), ga(ii,1-ji,ki), ga(ii,1-ji,1-ki))
     @             .eq. 1) then
                  fac1 = 1.d0 / ir
                  ffac = 0.5d0 * (ir-1)
                  center = center + ffac
                  do 60 m = idir, idir*(ir-1), idir
                     fac = (ir-abs(m)) * fac1
 60                  sum = sum + fac *
     @                    (hxm2 *
     @       (sigmaf(if+m-1,jf+ji-1,kf+ki-1) *
     @         (fdst(if+m-1,jf,kf) - fdst(if+m,jf,kf)) +
     @        sigmaf(if+m,jf+ji-1,kf+ki-1) *
     @         (fdst(if+m+1,jf,kf) - fdst(if+m,jf,kf))) +
     @                     hym2 *
     @       (sigmaf(if+m-1,jf+ji-1,kf+ki-1) + sigmaf(if+m,jf+ji-1,kf+ki-1)) *
     @         (fdst(if+m,jf+jdir,kf) - fdst(if+m,jf,kf)) +
     @                     hzm2 *
     @       (sigmaf(if+m-1,jf+ji-1,kf+ki-1) + sigmaf(if+m,jf+ji-1,kf+ki-1)) *
     @         (fdst(if+m,jf,kf+kdir) - fdst(if+m,jf,kf)))
               endif
               if (ga(ii,ji,ki) -
     @             min(ga(ii,ji,1-ki), ga(1-ii,ji,ki), ga(1-ii,ji,1-ki))
     @             .eq. 1) then
                  fac1 = 1.d0 / jr
                  ffac = 0.5d0 * (jr-1)
                  center = center + ffac
                  do 70 n = jdir, jdir*(jr-1), jdir
                     fac = (jr-abs(n)) * fac1
 70                  sum = sum + fac *
     @                    (hxm2 *
     @       (sigmaf(if+ii-1,jf+n-1,kf+ki-1) + sigmaf(if+ii-1,jf+n,kf+ki-1)) *
     @         (fdst(if+idir,jf+n,kf) - fdst(if,jf+n,kf)) +
     @                     hym2 *
     @       (sigmaf(if+ii-1,jf+n-1,kf+ki-1) *
     @         (fdst(if,jf+n-1,kf) - fdst(if,jf+n,kf)) +
     @        sigmaf(if+ii-1,jf+n,kf+ki-1) *
     @         (fdst(if,jf+n+1,kf) - fdst(if,jf+n,kf))) +
     @                     hzm2 *
     @       (sigmaf(if+ii-1,jf+n-1,kf+ki-1) + sigmaf(if+ii-1,jf+n,kf+ki-1)) *
     @         (fdst(if,jf+n,kf+kdir) - fdst(if,jf+n,kf)))
               endif
               if (ga(ii,ji,ki) -
     @             min(ga(ii,1-ji,ki), ga(1-ii,ji,ki), ga(1-ii,1-ji,ki))
     @             .eq. 1) then
                  fac1 = 1.d0 / kr
                  ffac = 0.5d0 * (kr-1)
                  center = center + ffac
                  do 80 l = kdir, kdir*(kr-1), kdir
                     fac = (kr-abs(l)) * fac1
 80                  sum = sum + fac *
     @                    (hxm2 *
     @       (sigmaf(if+ii-1,jf+ji-1,kf+l-1) + sigmaf(if+ii-1,jf+ji-1,kf+l)) *
     @         (fdst(if+idir,jf,kf+l) - fdst(if,jf,kf+l)) +
     @                     hym2 *
     @       (sigmaf(if+ii-1,jf+ji-1,kf+l-1) + sigmaf(if+ii-1,jf+ji-1,kf+l)) *
     @         (fdst(if,jf+jdir,kf+l) - fdst(if,jf,kf+l)) +
     @                     hzm2 *
     @       (sigmaf(if+ii-1,jf+ji-1,kf+l-1) *
     @         (fdst(if,jf,kf+l-1) - fdst(if,jf,kf+l)) +
     @        sigmaf(if+ii-1,jf+ji-1,kf+l) *
     @         (fdst(if,jf,kf+l+1) - fdst(if,jf,kf+l))))
               endif
 90            continue
c weighting
      res(if,jf,kf) = src(if,jf,kf) - sum / center
      end

c end of variable density non-terrain version
#endif

c constant density versions:
#else

c Note---assumes fdst linearly interpolated from cdst along face
      subroutine hgfres(res, dims(res),
     @ src,  dims(src),
     @ fdst, dims(fdst),
     @ cdst, dims(cdst),
     @ dims(reg),
     @ hx, irat, idim, idir)
      implicit real*8(a-h,o-z), integer(i-n)
      integer dims(res)
      integer dims(src)
      integer dims(fdst)
      integer dims(cdst)
      integer dims(reg)
      real*8 res(dimdec(res))
      real*8 src(dimdec(src))
      real*8 fdst(dimdec(fdst))
      real*8 cdst(dimdec(cdst))
      real*8 hx
      integer irat, idim, idir
      hm2 = 1.d0 / (hx*hx)
      fac0 = 0.5d0 / (irat*irat+irat)
      if (idim .eq. 0) then
         i = regl0
         do 10 k = regl2, regh2
            do 10 j = regl1, regh1
 10            res(i*irat,j*irat,k*irat) =
     @           src(i*irat,j*irat,k*irat) - 2.d0 * fac0 * hm2 *
     @               (2.d0 * cdst(i-idir,j,k) +
     @                cdst(i,j-1,k) + cdst(i,j+1,k) +
     @                cdst(i,j,k-1) + cdst(i,j,k+1) -
     @                6.d0 * cdst(i,j,k))
         fac0 = fac0 / (irat*irat*irat)
         i = i * irat
         do 20 l = 0, irat-1
            fac2 = (irat-l) * fac0
            if (l .eq. 0) fac2 = 0.5d0 * fac2
            do 20 n = 0, irat-1
               fac1 = (irat-n) * fac2
               if (n .eq. 0) fac1 = 0.5d0 * fac1
               do 20 k = irat*regl2, irat*regh2, irat
                  do 20 j = irat*regl1, irat*regh1, irat
 20            res(i,j,k) = res(i,j,k) - 2.d0 * fac1 * hm2 *
     @                (2.d0 *
     @                 (fdst(i+idir,j-n,k-l) + fdst(i+idir,j-n,k+l) +
     @                  fdst(i+idir,j+n,k-l) + fdst(i+idir,j+n,k+l)) +
     @                 (fdst(i,j-n-1,k-l) + fdst(i,j-n+1,k-l) +
     @                  fdst(i,j-n-1,k+l) + fdst(i,j-n+1,k+l) +
     @                  fdst(i,j+n-1,k-l) + fdst(i,j+n+1,k-l) +
     @                  fdst(i,j+n-1,k+l) + fdst(i,j+n+1,k+l)) +
     @                 (fdst(i,j-n,k-l-1) + fdst(i,j-n,k-l+1) +
     @                  fdst(i,j-n,k+l-1) + fdst(i,j-n,k+l+1) +
     @                  fdst(i,j+n,k-l-1) + fdst(i,j+n,k-l+1) +
     @                  fdst(i,j+n,k+l-1) + fdst(i,j+n,k+l+1)) - 
     @                 6.d0 *
     @                 (fdst(i,j-n,k-l) + fdst(i,j-n,k+l) +
     @                  fdst(i,j+n,k-l) + fdst(i,j+n,k+l)))
      else if (idim .eq. 1) then
         j = regl1
         do 30 k = regl2, regh2
            do 30 i = regl0, regh0
 30            res(i*irat,j*irat,k*irat) =
     @           src(i*irat,j*irat,k*irat) - 2.d0 * fac0 * hm2 *
     @               (cdst(i-1,j,k) + cdst(i+1,j,k) +
     @                2.d0 * cdst(i,j-idir,k) +
     @                cdst(i,j,k-1) + cdst(i,j,k+1) -
     @                6.d0 * cdst(i,j,k))
         fac0 = fac0 / (irat*irat*irat)
         j = j * irat
         do 40 l = 0, irat-1
            fac2 = (irat-l) * fac0
            if (l .eq. 0) fac2 = 0.5d0 * fac2
            do 40 m = 0, irat-1
               fac1 = (irat-m) * fac2
               if (m .eq. 0) fac1 = 0.5d0 * fac1
               do 40 k = irat*regl2, irat*regh2, irat
                  do 40 i = irat*regl0, irat*regh0, irat
 40            res(i,j,k) = res(i,j,k) - 2.d0 * fac1 * hm2 *
     @                ((fdst(i-m-1,j,k-l) + fdst(i-m+1,j,k-l) +
     @                  fdst(i-m-1,j,k+l) + fdst(i-m+1,j,k+l) +
     @                  fdst(i+m-1,j,k-l) + fdst(i+m+1,j,k-l) +
     @                  fdst(i+m-1,j,k+l) + fdst(i+m+1,j,k+l)) +
     @                 2.d0 *
     @                 (fdst(i-m,j+idir,k-l) + fdst(i-m,j+idir,k+l) +
     @                  fdst(i+m,j+idir,k-l) + fdst(i+m,j+idir,k+l)) +
     @                 (fdst(i-m,j,k-l-1) + fdst(i-m,j,k-l+1) +
     @                  fdst(i-m,j,k+l-1) + fdst(i-m,j,k+l+1) +
     @                  fdst(i+m,j,k-l-1) + fdst(i+m,j,k-l+1) +
     @                  fdst(i+m,j,k+l-1) + fdst(i+m,j,k+l+1)) -
     @                 6.d0 *
     @                 (fdst(i-m,j,k-l) + fdst(i-m,j,k+l) +
     @                  fdst(i+m,j,k-l) + fdst(i+m,j,k+l)))
      else
         k = regl2
         do 50 j = regl1, regh1
            do 50 i = regl0, regh0
 50            res(i*irat,j*irat,k*irat) =
     @           src(i*irat,j*irat,k*irat) - 2.d0 * fac0 * hm2 *
     @               (cdst(i-1,j,k) + cdst(i+1,j,k) +
     @                cdst(i,j-1,k) + cdst(i,j+1,k) +
     @                2.d0 * cdst(i,j,k-idir) -
     @                6.d0 * cdst(i,j,k))
         fac0 = fac0 / (irat*irat*irat)
         k = k * irat
         do 60 n = 0, irat-1
            fac2 = (irat-n) * fac0
            if (n .eq. 0) fac2 = 0.5d0 * fac2
            do 60 m = 0, irat-1
               fac1 = (irat-m) * fac2
               if (m .eq. 0) fac1 = 0.5d0 * fac1
               do 60 j = irat*regl1, irat*regh1, irat
                  do 60 i = irat*regl0, irat*regh0, irat
 60            res(i,j,k) = res(i,j,k) - 2.d0 * fac1 * hm2 *
     @                ((fdst(i-m-1,j-n,k) + fdst(i-m+1,j-n,k) +
     @                  fdst(i-m-1,j+n,k) + fdst(i-m+1,j+n,k) +
     @                  fdst(i+m-1,j-n,k) + fdst(i+m+1,j-n,k) +
     @                  fdst(i+m-1,j+n,k) + fdst(i+m+1,j+n,k)) +
     @                 (fdst(i-m,j-n-1,k) + fdst(i-m,j-n+1,k) +
     @                  fdst(i-m,j+n-1,k) + fdst(i-m,j+n+1,k) +
     @                  fdst(i+m,j-n-1,k) + fdst(i+m,j-n+1,k) +
     @                  fdst(i+m,j+n-1,k) + fdst(i+m,j+n+1,k)) +
     @                 2.d0 *
     @                 (fdst(i-m,j-n,k+idir) + fdst(i-m,j+n,k+idir) +
     @                  fdst(i+m,j-n,k+idir) + fdst(i+m,j+n,k+idir)) -
     @                 6.d0 *
     @                 (fdst(i-m,j-n,k) + fdst(i-m,j+n,k) +
     @                  fdst(i+m,j-n,k) + fdst(i+m,j+n,k)))
      endif
      end

c Note---assumes fdst linearly interpolated from cdst along face
      subroutine hgeres(res, dims(res),
     @ src,  dims(src),
     @ fdst, dims(fdst),
     @ cdst, dims(cdst),
     @ dims(reg),
     @ hx, irat, ivect, ga)
      implicit real*8(a-h,o-z), integer(i-n)
      integer dims(res)
      integer dims(src)
      integer dims(fdst)
      integer dims(cdst)
      integer dims(reg)
      real*8 res(dimdec(res))
      real*8 src(dimdec(src))
      real*8 fdst(dimdec(fdst))
      real*8 cdst(dimdec(cdst))
      real*8 hx
      integer irat, ivect(0:2), ga(0:1,0:1,0:1)
      r2 = irat * irat
      rm1 = 1.d0 / irat
      rm2 = rm1 * rm1
      hm2 = 1.d0 / (hx*hx)
      ic = regl0
      jc = regl1
      kc = regl2
      if = ic * irat
      jf = jc * irat
      kf = kc * irat
      center = 0.d0
      if (ivect(0) .eq. 0) then
         do 10 if = irat*regl0, irat*regh0, irat
 10         res(if,jf,kf) = 0.d0
c quadrants
c each quadrant is two octants and their share of the two central edges
         fac0 = rm2
         ffac = 1.d0
         cfac = r2
         do 40 ki = 0, 1
            kdir = 2 * ki - 1
            do 40 ji = 0, 1
               jdir = 2 * ji - 1
               if (ga(0,ji,ki) .eq. 1) then
                  center = center + ffac
                  do 20 m = 0, irat-1
                     fac = (irat-m) * fac0
                     if (m .eq. 0) fac = 0.5d0 * fac
                     do 20 if = irat*regl0, irat*regh0, irat
 20                     res(if,jf,kf) = res(if,jf,kf) + hm2 * fac *
     @                    ((fdst(if-m-1,jf,kf) + fdst(if-m+1,jf,kf) +
     @                      fdst(if+m-1,jf,kf) + fdst(if+m+1,jf,kf)) +
     @                     2.d0 *
     @                     (fdst(if-m,jf+jdir,kf) + fdst(if+m,jf+jdir,kf) +
     @                      fdst(if-m,jf,kf+kdir) + fdst(if+m,jf,kf+kdir)) -
     @                     6.d0 *
     @                     (fdst(if-m,jf,kf) + fdst(if+m,jf,kf)))
               else
                  center = center + cfac
                  do 30 ic = regl0, regh0
                     if = ic * irat
 30                  res(if,jf,kf) = res(if,jf,kf) + hm2 *
     @                   (cdst(ic-1,jc,kc) + cdst(ic,jc+jdir,kc) +
     @                    cdst(ic,jc,kc+kdir) + cdst(ic+1,jc,kc) +
     @                    cdst(ic,jc+jdir,kc) + cdst(ic,jc,kc+kdir) -
     @                    6.d0 * cdst(ic,jc,kc))
               endif
 40            continue
c faces
c each face is two faces and two sides of an edge
         fac0 = rm1 * rm2
         ffac = irat - 1
         do 70 ki = 0, 1
            kdir = 2 * ki - 1
            do 70 ji = 0, 1
               jdir = 2 * ji - 1
               if (ga(0,ji,ki) - ga(0,ji,1-ki) .eq. 1) then
                  center = center + ffac
                  do 50 n = jdir, jdir*(irat-1), jdir
                     fac1 = (irat-abs(n)) * fac0
                     do 50 m = 0, irat-1
                        fac = (irat-m) * fac1
                        if (m .eq. 0) fac = 0.5d0 * fac
                        do 50 if = irat*regl0, irat*regh0, irat
 50                        res(if,jf,kf) = res(if,jf,kf) + 2.d0 * fac * hm2 *
     @       ((fdst(if-m-1,jf+n,kf) + fdst(if-m+1,jf+n,kf) +
     @         fdst(if+m-1,jf+n,kf) + fdst(if+m+1,jf+n,kf)) +
     @        (fdst(if-m,jf+n-1,kf) + fdst(if-m,jf+n+1,kf) +
     @         fdst(if+m,jf+n-1,kf) + fdst(if+m,jf+n+1,kf)) +
     @        2.d0 *
     @        (fdst(if-m,jf+n,kf+kdir) + fdst(if+m,jf+n,kf+kdir)) -
     @        6.d0 *
     @        (fdst(if-m,jf+n,kf) + fdst(if+m,jf+n,kf)))
               endif
               if (ga(0,ji,ki) - ga(0,1-ji,ki) .eq. 1) then
                  center = center + ffac
                  do 60 l = kdir, kdir*(irat-1), kdir
                     fac1 = (irat-abs(l)) * fac0
                     do 60 m = 0, irat-1
                        fac = (irat-m) * fac1
                        if (m .eq. 0) fac = 0.5d0 * fac
                        do 60 if = irat*regl0, irat*regh0, irat
 60                        res(if,jf,kf) = res(if,jf,kf) + 2.d0 * fac * hm2 *
     @       ((fdst(if-m-1,jf,kf+l) + fdst(if-m+1,jf,kf+l) +
     @         fdst(if+m-1,jf,kf+l) + fdst(if+m+1,jf,kf+l)) +
     @        (fdst(if-m,jf,kf+l-1) + fdst(if-m,jf,kf+l+1) +
     @         fdst(if+m,jf,kf+l-1) + fdst(if+m,jf,kf+l+1)) +
     @        2.d0 *
     @        (fdst(if-m,jf+jdir,kf+l) + fdst(if+m,jf+jdir,kf+l)) -
     @        6.d0 *
     @        (fdst(if-m,jf,kf+l) + fdst(if+m,jf,kf+l)))
               endif
 70            continue
c weighting
         do 80 if = irat*regl0, irat*regh0, irat
 80         res(if,jf,kf) = src(if,jf,kf) - res(if,jf,kf) / center
      else if (ivect(1) .eq. 0) then
         do 110 jf = irat*regl1, irat*regh1, irat
 110        res(if,jf,kf) = 0.d0
c quadrants
c each quadrant is two octants and their share of the two central edges
         fac0 = rm2
         ffac = 1.d0
         cfac = r2
         do 140 ki = 0, 1
            kdir = 2 * ki - 1
            do 140 ii = 0, 1
               idir = 2 * ii - 1
               if (ga(ii,0,ki) .eq. 1) then
                  center = center + ffac
                  do 120 n = 0, irat-1
                     fac = (irat-n) * fac0
                     if (n .eq. 0) fac = 0.5d0 * fac
                     do 120 jf = irat*regl1, irat*regh1, irat
 120                    res(if,jf,kf) = res(if,jf,kf) + hm2 * fac *
     @                    ((fdst(if,jf-n-1,kf) + fdst(if,jf-n+1,kf) +
     @                      fdst(if,jf+n-1,kf) + fdst(if,jf+n+1,kf)) +
     @                     2.d0 *
     @                     (fdst(if+idir,jf-n,kf) + fdst(if+idir,jf+n,kf) +
     @                      fdst(if,jf-n,kf+kdir) + fdst(if,jf+n,kf+kdir)) -
     @                     6.d0 *
     @                     (fdst(if,jf-n,kf) + fdst(if,jf+n,kf)))
               else
                  center = center + cfac
                  do 130 jc = regl1, regh1
                     jf = jc * irat
 130                 res(if,jf,kf) = res(if,jf,kf) + hm2 *
     @                   (cdst(ic+idir,jc,kc) + cdst(ic,jc-1,kc) +
     @                    cdst(ic,jc,kc+kdir) + cdst(ic+idir,jc,kc) +
     @                    cdst(ic,jc+1,kc) + cdst(ic,jc,kc+kdir) -
     @                    6.d0 * cdst(ic,jc,kc))
               endif
 140           continue
c faces
c each face is two faces and two sides of an edge
         fac0 = rm1 * rm2
         ffac = irat - 1
         do 170 ki = 0, 1
            kdir = 2 * ki - 1
            do 170 ii = 0, 1
               idir = 2 * ii - 1
               if (ga(ii,0,ki) - ga(ii,0,1-ki) .eq. 1) then
                  center = center + ffac
                  do 150 n = 0, irat-1
                     fac1 = (irat-n) * fac0
                     if (n .eq. 0) fac1 = 0.5d0 * fac1
                     do 150 m = idir, idir*(irat-1), idir
                        fac = (irat-abs(m)) * fac1
                        do 150 jf = irat*regl1, irat*regh1, irat
 150                       res(if,jf,kf) = res(if,jf,kf) + 2.d0 * fac * hm2 *
     @       ((fdst(if+m-1,jf-n,kf) + fdst(if+m+1,jf-n,kf) +
     @         fdst(if+m-1,jf+n,kf) + fdst(if+m+1,jf+n,kf)) +
     @        (fdst(if+m,jf-n-1,kf) + fdst(if+m,jf-n+1,kf) +
     @         fdst(if+m,jf+n-1,kf) + fdst(if+m,jf+n+1,kf)) +
     @        2.d0 *
     @        (fdst(if+m,jf-n,kf+kdir) + fdst(if+m,jf+n,kf+kdir)) -
     @        6.d0 *
     @        (fdst(if+m,jf-n,kf) + fdst(if+m,jf+n,kf)))
               endif
               if (ga(ii,0,ki) - ga(1-ii,0,ki) .eq. 1) then
                  center = center + ffac
                  do 160 l = kdir, kdir*(irat-1), kdir
                     fac1 = (irat-abs(l)) * fac0
                     do 160 n = 0, irat-1
                        fac = (irat-n) * fac1
                        if (n .eq. 0) fac = 0.5d0 * fac
                        do 160 jf = irat*regl1, irat*regh1, irat
 160                       res(if,jf,kf) = res(if,jf,kf) + 2.d0 * fac * hm2 *
     @       (2.d0 *
     @        (fdst(if+idir,jf-n,kf+l) + fdst(if+idir,jf+n,kf+l)) +
     @        (fdst(if,jf-n-1,kf+l) + fdst(if,jf-n+1,kf+l) +
     @         fdst(if,jf+n-1,kf+l) + fdst(if,jf+n+1,kf+l)) +
     @        (fdst(if,jf-n,kf+l-1) + fdst(if,jf-n,kf+l+1) +
     @         fdst(if,jf+n,kf+l-1) + fdst(if,jf+n,kf+l+1)) -
     @        6.d0 *
     @        (fdst(if,jf-n,kf+l) + fdst(if,jf+n,kf+l)))
               endif
 170           continue
c weighting
         do 180 jf = irat*regl1, irat*regh1, irat
 180        res(if,jf,kf) = src(if,jf,kf) - res(if,jf,kf) / center
      else
         do 210 kf = irat*regl2, irat*regh2, irat
 210        res(if,jf,kf) = 0.d0
c quadrants
c each quadrant is two octants and their share of the two central edges
         fac0 = rm2
         ffac = 1.d0
         cfac = r2
         do 240 ji = 0, 1
            jdir = 2 * ji - 1
            do 240 ii = 0, 1
               idir = 2 * ii - 1
               if (ga(ii,ji,0) .eq. 1) then
                  center = center + ffac
                  do 220 l = 0, irat-1
                     fac = (irat-l) * fac0
                     if (l .eq. 0) fac = 0.5d0 * fac
                     do 220 kf = irat*regl2, irat*regh2, irat
 220                    res(if,jf,kf) = res(if,jf,kf) + hm2 * fac *
     @                    ((fdst(if,jf,kf-l-1) + fdst(if,jf,kf-l+1) +
     @                      fdst(if,jf,kf+l-1) + fdst(if,jf,kf+l+1)) +
     @                     2.d0 *
     @                     (fdst(if+idir,jf,kf-l) + fdst(if+idir,jf,kf+l) +
     @                      fdst(if,jf+jdir,kf-l) + fdst(if,jf+jdir,kf+l)) -
     @                     6.d0 *
     @                     (fdst(if,jf,kf-l) + fdst(if,jf,kf+l)))
               else
                  center = center + cfac
                  do 230 kc = regl2, regh2
                     kf = kc * irat
 230                 res(if,jf,kf) = res(if,jf,kf) + hm2 *
     @                   (cdst(ic+idir,jc,kc) + cdst(ic,jc+jdir,kc) +
     @                    cdst(ic,jc,kc-1) + cdst(ic+idir,jc,kc) +
     @                    cdst(ic,jc+jdir,kc) + cdst(ic,jc,kc+1) -
     @                    6.d0 * cdst(ic,jc,kc))
               endif
 240           continue
c faces
c each face is two faces and two sides of an edge
         fac0 = rm1 * rm2
         ffac = irat - 1
         do 270 ji = 0, 1
            jdir = 2 * ji - 1
            do 270 ii = 0, 1
               idir = 2 * ii - 1
               if (ga(ii,ji,0) - ga(ii,1-ji,0) .eq. 1) then
                  center = center + ffac
                  do 250 l = 0, irat-1
                     fac1 = (irat-l) * fac0
                     if (l .eq. 0) fac1 = 0.5d0 * fac1
                     do 250 m = idir, idir*(irat-1), idir
                        fac = (irat-abs(m)) * fac1
                        do 250 kf = irat*regl2, irat*regh2, irat
 250                       res(if,jf,kf) = res(if,jf,kf) + 2.d0 * fac * hm2 *
     @       ((fdst(if+m-1,jf,kf-l) + fdst(if+m+1,jf,kf-l) +
     @         fdst(if+m-1,jf,kf+l) + fdst(if+m+1,jf,kf+l)) +
     @        (fdst(if+m,jf,kf-l-1) + fdst(if+m,jf,kf-l+1) +
     @         fdst(if+m,jf,kf+l-1) + fdst(if+m,jf,kf+l+1)) +
     @        2.d0 *
     @        (fdst(if+m,jf+jdir,kf-l) + fdst(if+m,jf+jdir,kf+l)) -
     @        6.d0 *
     @        (fdst(if+m,jf,kf-l) + fdst(if+m,jf,kf+l)))
               endif
               if (ga(ii,ji,0) - ga(1-ii,ji,0) .eq. 1) then
                  center = center + ffac
                  do 260 l = 0, irat-1
                     fac1 = (irat-l) * fac0
                     if (l .eq. 0) fac1 = 0.5d0 * fac1
                     do 260 n = jdir, jdir*(irat-1), jdir
                        fac = (irat-abs(n)) * fac1
                        do 260 kf = irat*regl2, irat*regh2, irat
 260                       res(if,jf,kf) = res(if,jf,kf) + 2.d0 * fac * hm2 *
     @       (2.d0 *
     @        (fdst(if+idir,jf+n,kf-l) + fdst(if+idir,jf+n,kf+l)) +
     @        (fdst(if,jf+n-1,kf-l) + fdst(if,jf+n+1,kf-l) +
     @         fdst(if,jf+n-1,kf+l) + fdst(if,jf+n+1,kf+l)) +
     @        (fdst(if,jf+n,kf-l-1) + fdst(if,jf+n,kf-l+1) +
     @         fdst(if,jf+n,kf+l-1) + fdst(if,jf+n,kf+l+1)) -
     @        6.d0 *
     @        (fdst(if,jf+n,kf-l) + fdst(if,jf+n,kf+l)))
               endif
 270           continue
c weighting
         do 280 kf = irat*regl2, irat*regh2, irat
 280        res(if,jf,kf) = src(if,jf,kf) - res(if,jf,kf) / center
      endif
      end

c Note---assumes fdst linearly interpolated from cdst along face
      subroutine hgcres(res, dims(res),
     @ src,  dims(src),
     @ fdst, dims(fdst),
     @ cdst, dims(cdst),
     @ dims(reg),
     @ hx, irat, ga)
      implicit real*8(a-h,o-z), integer(i-n)
      integer dims(res)
      integer dims(src)
      integer dims(fdst)
      integer dims(cdst)
      integer dims(reg)
      real*8 res(dimdec(res))
      real*8 src(dimdec(src))
      real*8 fdst(dimdec(fdst))
      real*8 cdst(dimdec(cdst))
      real*8 hx
      integer irat, ga(0:1,0:1,0:1)
      r2 = irat * irat
      rm1 = 1.d0 / irat
      rm2 = rm1 * rm1
      hm2 = 1.d0 / (hx*hx)
      ic = regl0
      jc = regl1
      kc = regl2
      if = ic * irat
      jf = jc * irat
      kf = kc * irat
      sum = 0.d0
      center = 0.d0
c octants
      fac = rm1
      ffac = 0.5d0 * rm1
      cfac = 0.5d0 * r2
      do 10 ki = 0, 1
         kdir = 2 * ki - 1
         do 10 ji = 0, 1
            jdir = 2 * ji - 1
            do 10 ii = 0, 1
               idir = 2 * ii - 1
               if (ga(ii,ji,ki) .eq. 1) then
                  center = center + ffac
                  sum = sum + fac * hm2 *
     @              (fdst(if+idir,jf,kf) + fdst(if,jf+jdir,kf) +
     @               fdst(if,jf,kf+kdir) - 3.d0 * fdst(if,jf,kf))
               else
                  center = center + cfac
                  sum = sum + hm2 *
     @              (cdst(ic+idir,jc,kc) + cdst(ic,jc+jdir,kc) +
     @               cdst(ic,jc,kc+kdir) - 3.d0 * cdst(ic,jc,kc))
               endif
 10            continue
c faces
      fac2 = rm1 * rm2
      ffac = 0.5d0 * (irat-1) * (irat-1) * rm1
      do 50 ki = 0, 1
         kdir = 2 * ki - 1
         do 50 ji = 0, 1
            jdir = 2 * ji - 1
            do 50 ii = 0, 1
               idir = 2 * ii - 1
               if (ga(ii,ji,ki) - ga(ii,ji,1-ki) .eq. 1) then
                  center = center + ffac
                  do 20 n = jdir, jdir*(irat-1), jdir
                     fac1 = (irat-abs(n)) * fac2
                     do 20 m = idir, idir*(irat-1), idir
                        fac = (irat-abs(m)) * fac1
 20                     sum = sum + 2.d0 * fac * hm2 *
     @       (fdst(if+m-1,jf+n,kf) + fdst(if+m+1,jf+n,kf) +
     @        fdst(if+m,jf+n-1,kf) + fdst(if+m,jf+n+1,kf) +
     @        2.d0 * fdst(if+m,jf+n,kf+kdir) -
     @        6.d0 * fdst(if+m,jf+n,kf))
               endif
               if (ga(ii,ji,ki) - ga(ii,1-ji,ki) .eq. 1) then
                  center = center + ffac
                  do 30 l = kdir, kdir*(irat-1), kdir
                     fac1 = (irat-abs(l)) * fac2
                     do 30 m = idir, idir*(irat-1), idir
                        fac = (irat-abs(m)) * fac1
 30                     sum = sum + 2.d0 * fac * hm2 *
     @       (fdst(if+m-1,jf,kf+l) + fdst(if+m+1,jf,kf+l) +
     @        fdst(if+m,jf,kf+l-1) + fdst(if+m,jf,kf+l+1) +
     @        2.d0 * fdst(if+m,jf+jdir,kf+l) -
     @        6.d0 * fdst(if+m,jf,kf+l))
               endif
               if (ga(ii,ji,ki) - ga(1-ii,ji,ki) .eq. 1) then
                  center = center + ffac
                  do 40 l = kdir, kdir*(irat-1), kdir
                     fac1 = (irat-abs(l)) * fac2
                     do 40 n = jdir, jdir*(irat-1), jdir
                        fac = (irat-abs(n)) * fac1
 40                     sum = sum + 2.d0 * fac * hm2 *
     @       (fdst(if,jf+n-1,kf+l) + fdst(if,jf+n+1,kf+l) +
     @        fdst(if,jf+n,kf+l-1) + fdst(if,jf+n,kf+l+1) +
     @        2.d0 * fdst(if+idir,jf+n,kf+l) -
     @        6.d0 * fdst(if,jf+n,kf+l))
            endif
 50         continue
c edges
      fac1 = rm2
      ffac = 0.5d0 * (irat-1) * rm1
      do 90 ki = 0, 1
         kdir = 2 * ki - 1
         do 90 ji = 0, 1
            jdir = 2 * ji - 1
            do 90 ii = 0, 1
               idir = 2 * ii - 1
               if (ga(ii,ji,ki) -
     @             min(ga(ii,ji,1-ki), ga(ii,1-ji,ki), ga(ii,1-ji,1-ki))
     @             .eq. 1) then
                  center = center + ffac
                  do 60 m = idir, idir*(irat-1), idir
                     fac = (irat-abs(m)) * fac1
 60                  sum = sum + fac * hm2 *
     @       (fdst(if+m-1,jf,kf) + fdst(if+m+1,jf,kf) +
     @        2.d0 * (fdst(if+m,jf+jdir,kf) +
     @                fdst(if+m,jf,kf+kdir)) -
     @        6.d0 * fdst(if+m,jf,kf))
               endif
               if (ga(ii,ji,ki) -
     @             min(ga(ii,ji,1-ki), ga(1-ii,ji,ki), ga(1-ii,ji,1-ki))
     @             .eq. 1) then
                  center = center + ffac
                  do 70 n = jdir, jdir*(irat-1), jdir
                     fac = (irat-abs(n)) * fac1
 70                  sum = sum + fac * hm2 *
     @       (2.d0 * (fdst(if+idir,jf+n,kf) +
     @                fdst(if,jf+n,kf+kdir)) +
     @        fdst(if,jf+n-1,kf) + fdst(if,jf+n+1,kf) -
     @        6.d0 * fdst(if,jf+n,kf))
               endif
               if (ga(ii,ji,ki) -
     @             min(ga(ii,1-ji,ki), ga(1-ii,ji,ki), ga(1-ii,1-ji,ki))
     @             .eq. 1) then
                  center = center + ffac
                  do 80 l = kdir, kdir*(irat-1), kdir
                     fac = (irat-abs(l)) * fac1
 80                  sum = sum + fac * hm2 *
     @       (2.d0 * (fdst(if+idir,jf,kf+l) +
     @                fdst(if,jf+jdir,kf+l)) +
     @        fdst(if,jf,kf+l-1) + fdst(if,jf,kf+l+1) -
     @        6.d0 * fdst(if,jf,kf+l))
               endif
 90            continue
c weighting
      res(if,jf,kf) = src(if,jf,kf) - sum / center
      end

#endif

#ifndef HG_CONSTANT

#  ifndef HG_SIGMA_NODE

c NODE-based data, factor of 2 only.
      subroutine hgints(dest,
     @ dims(dest),
     @ dims(reg),
     @ sigx, sigy, sigz,
     @ dims(sb),
     @ src,
     @ dims(src),
     @ dims(bb), ir, jr, kr)
      integer dims(dest)
      integer dims(reg)
      integer dims(sb)
      integer dims(src)
      integer dims(bb)
      integer ir, jr, kr
      real*8 dest(dimdec(dest))
      real*8 sigx(dimdec(sb))
      real*8 sigy(dimdec(sb))
      real*8 sigz(dimdec(sb))
      real*8 src(dimdec(src))
      integer i, j, k, ic, jc, kc
      do 10 kc = bbl2, bbh2
         do 10 jc = bbl1, bbh1
            do 10 ic = bbl0, bbh0
 10            dest(ir*ic,jr*jc,kr*kc) = src(ic,jc,kc)
      if (ir .eq. 2) then
         do 20 kc = bbl2, bbh2
            do 20 jc = bbl1, bbh1
               do 20 ic = bbl0, bbh0-1
                  i = ir * ic
                  j = jr * jc
                  k = kr * kc
 20               dest(i+1,j,k) = ((sigx(i,j-1,k-1) + sigx(i,j-1,k) +
     @                              sigx(i,j,k-1) + sigx(i,j,k)) *
     @                             src(ic,jc,kc) +
     @                             (sigx(i+1,j-1,k-1) + sigx(i+1,j-1,k) +
     @                              sigx(i+1,j,k-1) + sigx(i+1,j,k)) *
     @                             src(ic+1,jc,kc)) /
     @                             (sigx(i,j-1,k-1) + sigx(i,j-1,k) +
     @                              sigx(i,j,k-1) + sigx(i,j,k) +
     @                              sigx(i+1,j-1,k-1) + sigx(i+1,j-1,k) +
     @                              sigx(i+1,j,k-1) + sigx(i+1,j,k))
      endif
      if (jr .eq. 2) then
         do 30 kc = bbl2, bbh2
            do 30 jc = bbl1, bbh1-1
               do 30 ic = bbl0, bbh0
                  i = ir * ic
                  j = jr * jc
                  k = kr * kc
 30               dest(i,j+1,k) = ((sigy(i-1,j,k-1) + sigy(i-1,j,k) +
     @                              sigy(i,j,k-1) + sigy(i,j,k)) *
     @                             src(ic,jc,kc) +
     @                             (sigy(i-1,j+1,k-1) + sigy(i-1,j+1,k) +
     @                              sigy(i,j+1,k-1) + sigy(i,j+1,k)) *
     @                             src(ic,jc+1,kc)) /
     @                             (sigy(i-1,j,k-1) + sigy(i-1,j,k) +
     @                              sigy(i,j,k-1) + sigy(i,j,k) +
     @                              sigy(i-1,j+1,k-1) + sigy(i-1,j+1,k) +
     @                              sigy(i,j+1,k-1) + sigy(i,j+1,k))
      endif
      if (kr .eq. 2) then
         do 40 kc = bbl2, bbh2-1
            do 40 jc = bbl1, bbh1
               do 40 ic = bbl0, bbh0
                  i = ir * ic
                  j = jr * jc
                  k = kr * kc
 40               dest(i,j,k+1) = ((sigz(i-1,j-1,k) + sigz(i-1,j,k) +
     @                              sigz(i,j-1,k) + sigz(i,j,k)) *
     @                             src(ic,jc,kc) +
     @                             (sigz(i-1,j-1,k+1) + sigz(i-1,j,k+1) +
     @                              sigz(i,j-1,k+1) + sigz(i,j,k+1)) *
     @                             src(ic,jc,kc+1)) /
     @                             (sigz(i-1,j-1,k) + sigz(i-1,j,k) +
     @                              sigz(i,j-1,k) + sigz(i,j,k) +
     @                              sigz(i-1,j-1,k+1) + sigz(i-1,j,k+1) +
     @                              sigz(i,j-1,k+1) + sigz(i,j,k+1))
      endif
      if (ir .eq. 2 .and. jr .eq. 2) then
         do 50 kc = bbl2, bbh2
            do 50 jc = bbl1, bbh1-1
cdir$ ivdep
               do 50 ic = bbl0, bbh0-1
                  i = ir * ic
                  j = jr * jc
                  k = kr * kc
 50            dest(i+1,j+1,k) = ((sigx(i,j,k-1) + sigx(i,j,k) +
     @                             sigx(i,j+1,k-1) + sigx(i,j+1,k)) *
     @                            dest(i,j+1,k) +
     @                            (sigx(i+1,j,k-1) + sigx(i+1,j,k) +
     @                             sigx(i+1,j+1,k-1) + sigx(i+1,j+1,k)) *
     @                            dest(i+2,j+1,k) +
     @                            (sigy(i,j,k-1) + sigy(i,j,k) +
     @                             sigy(i+1,j,k-1) + sigy(i+1,j,k)) *
     @                            dest(i+1,j,k) +
     @                            (sigy(i,j+1,k-1) + sigy(i,j+1,k) +
     @                             sigy(i+1,j+1,k-1) + sigy(i+1,j+1,k)) *
     @                            dest(i+1,j+2,k)) /
     @                           (sigx(i,j,k-1) + sigx(i,j,k) +
     @                            sigx(i,j+1,k-1) + sigx(i,j+1,k) +
     @                            sigx(i+1,j,k-1) + sigx(i+1,j,k) +
     @                            sigx(i+1,j+1,k-1) + sigx(i+1,j+1,k) +
     @                            sigy(i,j,k-1) + sigy(i,j,k) +
     @                            sigy(i+1,j,k-1) + sigy(i+1,j,k) +
     @                            sigy(i,j+1,k-1) + sigy(i,j+1,k) +
     @                            sigy(i+1,j+1,k-1) + sigy(i+1,j+1,k))
      endif
      if (ir .eq. 2 .and. kr .eq. 2) then
         do 60 kc = bbl2, bbh2-1
            do 60 jc = bbl1, bbh1
cdir$ ivdep
               do 60 ic = bbl0, bbh0-1
                  i = ir * ic
                  j = jr * jc
                  k = kr * kc
 60            dest(i+1,j,k+1) = ((sigx(i,j-1,k) + sigx(i,j-1,k+1) +
     @                             sigx(i,j,k) + sigx(i,j,k+1)) *
     @                            dest(i,j,k+1) +
     @                            (sigx(i+1,j-1,k) + sigx(i+1,j-1,k+1) +
     @                             sigx(i+1,j,k) + sigx(i+1,j,k+1)) *
     @                            dest(i+2,j,k+1) +
     @                            (sigz(i,j-1,k) + sigz(i,j,k) +
     @                             sigz(i+1,j-1,k) + sigz(i+1,j,k)) *
     @                            dest(i+1,j,k) +
     @                            (sigz(i,j-1,k+1) + sigz(i,j,k+1) +
     @                             sigz(i+1,j-1,k+1) + sigz(i+1,j,k+1)) *
     @                            dest(i+1,j,k+2)) /
     @                           (sigx(i,j-1,k) + sigx(i,j-1,k+1) +
     @                            sigx(i,j,k) + sigx(i,j,k+1) +
     @                            sigx(i+1,j-1,k) + sigx(i+1,j-1,k+1) +
     @                            sigx(i+1,j,k) + sigx(i+1,j,k+1) +
     @                            sigz(i,j-1,k) + sigz(i,j,k) +
     @                            sigz(i+1,j-1,k) + sigz(i+1,j,k) +
     @                            sigz(i,j-1,k+1) + sigz(i,j,k+1) +
     @                            sigz(i+1,j-1,k+1) + sigz(i+1,j,k+1))
      endif
      if (jr .eq. 2 .and. kr .eq. 2) then
         do 70 kc = bbl2, bbh2-1
            do 70 jc = bbl1, bbh1-1
cdir$ ivdep
               do 70 ic = bbl0, bbh0
                  i = ir * ic
                  j = jr * jc
                  k = kr * kc
 70            dest(i,j+1,k+1) = ((sigy(i-1,j,k) + sigy(i-1,j,k+1) +
     @                             sigy(i,j,k) + sigy(i,j,k+1)) *
     @                            dest(i,j,k+1) +
     @                            (sigy(i-1,j+1,k) + sigy(i-1,j+1,k+1) +
     @                             sigy(i,j+1,k) + sigy(i,j+1,k+1)) *
     @                            dest(i,j+2,k+1) +
     @                            (sigz(i-1,j,k) + sigz(i-1,j+1,k) +
     @                             sigz(i,j,k) + sigz(i,j+1,k)) *
     @                            dest(i,j+1,k) +
     @                            (sigz(i-1,j,k+1) + sigz(i-1,j+1,k+1) +
     @                             sigz(i,j,k+1) + sigz(i,j+1,k+1)) *
     @                            dest(i,j+1,k+2)) /
     @                           (sigy(i-1,j,k) + sigy(i-1,j,k+1) +
     @                            sigy(i,j,k) + sigy(i,j,k+1) +
     @                            sigy(i-1,j+1,k) + sigy(i-1,j+1,k+1) +
     @                            sigy(i,j+1,k) + sigy(i,j+1,k+1) +
     @                            sigz(i-1,j,k) + sigz(i-1,j+1,k) +
     @                            sigz(i,j,k) + sigz(i,j+1,k) +
     @                            sigz(i-1,j,k+1) + sigz(i-1,j+1,k+1) +
     @                            sigz(i,j,k+1) + sigz(i,j+1,k+1))
      endif
      if (ir .eq. 2 .and. jr .eq. 2 .and. kr .eq. 2) then
         do 80 kc = bbl2, bbh2-1
            do 80 jc = bbl1, bbh1-1
cdir$ ivdep
               do 80 ic = bbl0, bbh0-1
                  i = ir * ic
                  j = jr * jc
                  k = kr * kc
               dest(i+1,j+1,k+1) = ((sigx(i,j,k) + sigx(i,j,k+1) +
     @                               sigx(i,j+1,k) + sigx(i,j+1,k+1)) *
     @                              dest(i,j+1,k+1) +
     @                              (sigx(i+1,j,k) + sigx(i+1,j,k+1) +
     @                               sigx(i+1,j+1,k) + sigx(i+1,j+1,k+1)) *
     @                              dest(i+2,j+1,k+1) +
     @                              (sigy(i,j,k) + sigy(i,j,k+1) +
     @                               sigy(i+1,j,k) + sigy(i+1,j,k+1)) *
     @                              dest(i+1,j,k+1) +
     @                              (sigy(i,j+1,k) + sigy(i,j+1,k+1) +
     @                               sigy(i+1,j+1,k) + sigy(i+1,j+1,k+1)) *
     @                              dest(i+1,j+2,k+1) +
     @                              (sigz(i,j,k) + sigz(i,j+1,k) +
     @                               sigz(i+1,j,k) + sigz(i+1,j+1,k)) *
     @                              dest(i+1,j+1,k) +
     @                              (sigz(i,j,k+1) + sigz(i,j+1,k+1) +
     @                               sigz(i+1,j,k+1) + sigz(i+1,j+1,k+1)) *
     @                              dest(i+1,j+1,k+2))
 80            dest(i+1,j+1,k+1) = dest(i+1,j+1,k+1) /
     @                             (sigx(i,j,k) + sigx(i,j,k+1) +
     @                              sigx(i,j+1,k) + sigx(i,j+1,k+1) +
     @                              sigx(i+1,j,k) + sigx(i+1,j,k+1) +
     @                              sigx(i+1,j+1,k) + sigx(i+1,j+1,k+1) +
     @                              sigy(i,j,k) + sigy(i,j,k+1) +
     @                              sigy(i+1,j,k) + sigy(i+1,j,k+1) +
     @                              sigy(i,j+1,k) + sigy(i,j+1,k+1) +
     @                              sigy(i+1,j+1,k) + sigy(i+1,j+1,k+1) +
     @                              sigz(i,j,k) + sigz(i,j+1,k) +
     @                              sigz(i+1,j,k) + sigz(i+1,j+1,k) +
     @                              sigz(i,j,k+1) + sigz(i,j+1,k+1) +
     @                              sigz(i+1,j,k+1) + sigz(i+1,j+1,k+1))
      endif
      end

#  else

c NODE-based data, factor of 2 only.
      subroutine hgints(dest,
     @ dims(dest),
     @ dims(reg),
     @ signd,
     @ dims(sn),
     @ src,
     @ dims(src),
     @ dims(bb), ir, jr, kr)
      integer dims(dest)
      integer dims(reg)
      integer dims(sn)
      integer dims(src)
      integer dims(bb)
      integer ir, jr, kr
      real*8 dest(dimdec(dest))
      real*8 signd(dimdec(sn), 3)
      real*8 src(dimdec(src))
      integer i, j, k, ic, jc, kc
      do 10 kc = bbl2, bbh2
         do 10 jc = bbl1, bbh1
            do 10 ic = bbl0, bbh0
 10            dest(ir*ic,jr*jc,kr*kc) = src(ic,jc,kc)
      if (ir .eq. 2) then
         do 20 kc = bbl2, bbh2
            do 20 jc = bbl1, bbh1
               do 20 ic = bbl0, bbh0-1
                  i = ir * ic
                  j = jr * jc
                  k = kr * kc
 20            dest(i+1,j,k) = (signd(i,j,k,1)   * src(ic,jc,kc) +
     @                          signd(i+1,j,k,1) * src(ic+1,jc,kc)) /
     @                         (signd(i,j,k,1) + signd(i+1,j,k,1))
      endif
      if (jr .eq. 2) then
         do 30 kc = bbl2, bbh2
            do 30 jc = bbl1, bbh1-1
               do 30 ic = bbl0, bbh0
                  i = ir * ic
                  j = jr * jc
                  k = kr * kc
 30            dest(i,j+1,k) = (signd(i,j,k,2)   * src(ic,jc,kc) +
     @                          signd(i,j+1,k,2) * src(ic,jc+1,kc)) /
     @                         (signd(i,j,k,2) + signd(i,j+1,k,2))
      endif
      if (kr .eq. 2) then
         do 40 kc = bbl2, bbh2-1
            do 40 jc = bbl1, bbh1
               do 40 ic = bbl0, bbh0
                  i = ir * ic
                  j = jr * jc
                  k = kr * kc
 40            dest(i,j,k+1) = (signd(i,j,k,3)   * src(ic,jc,kc) +
     @                          signd(i,j,k+1,3) * src(ic,jc,kc+1)) /
     @                         (signd(i,j,k,3) + signd(i,j,k+1,3))
      endif
      if (ir .eq. 2 .and. jr .eq. 2) then
         do 50 kc = bbl2, bbh2
            do 50 jc = bbl1, bbh1-1
cdir$ ivdep
               do 50 ic = bbl0, bbh0-1
                  i = ir * ic
                  j = jr * jc
                  k = kr * kc
 50            dest(i+1,j+1,k) = (signd(i,j+1,k,1)   * dest(i,j+1,k) +
     @                            signd(i+1,j+1,k,1) * dest(i+2,j+1,k) +
     @                            signd(i+1,j,k,2)   * dest(i+1,j,k) +
     @                            signd(i+1,j+1,k,2) * dest(i+1,j+2,k)) /
     @                           (signd(i,j+1,k,1) + signd(i+1,j+1,k,1) +
     @                            signd(i+1,j,k,2) + signd(i+1,j+1,k,2))
      endif
      if (ir .eq. 2 .and. kr .eq. 2) then
         do 60 kc = bbl2, bbh2-1
            do 60 jc = bbl1, bbh1
cdir$ ivdep
               do 60 ic = bbl0, bbh0-1
                  i = ir * ic
                  j = jr * jc
                  k = kr * kc
 60            dest(i+1,j,k+1) = (signd(i,j,k+1,1)   * dest(i,j,k+1) +
     @                            signd(i+1,j,k+1,1) * dest(i+2,j,k+1) +
     @                            signd(i+1,j,k,3)   * dest(i+1,j,k) +
     @                            signd(i+1,j,k+1,3) * dest(i+1,j,k+2)) /
     @                           (signd(i,j,k+1,1) + signd(i+1,j,k+1,1) +
     @                            signd(i+1,j,k,3) + signd(i+1,j,k+1,3))
      endif
      if (jr .eq. 2 .and. kr .eq. 2) then
         do 70 kc = bbl2, bbh2-1
            do 70 jc = bbl1, bbh1-1
cdir$ ivdep
               do 70 ic = bbl0, bbh0
                  i = ir * ic
                  j = jr * jc
                  k = kr * kc
 70            dest(i,j+1,k+1) = (signd(i,j,k+1,2)   * dest(i,j,k+1) +
     @                            signd(i,j+1,k+1,2) * dest(i,j+2,k+1) +
     @                            signd(i,j+1,k,3)   * dest(i,j+1,k) +
     @                            signd(i,j+1,k+1,3) * dest(i,j+1,k+2)) /
     @                           (signd(i,j,k+1,2) + signd(i,j+1,k+1,2) +
     @                            signd(i,j+1,k,3) + signd(i,j+1,k+1,3))
      endif
      if (ir .eq. 2 .and. jr .eq. 2 .and. kr .eq. 2) then
         do 80 kc = bbl2, bbh2-1
            do 80 jc = bbl1, bbh1-1
cdir$ ivdep
               do 80 ic = bbl0, bbh0-1
                  i = ir * ic
                  j = jr * jc
                  k = kr * kc
 80            dest(i+1,j+1,k+1) = (signd(i,j+1,k+1,1)   * dest(i,j+1,k+1) +
     @                              signd(i+1,j+1,k+1,1) * dest(i+2,j+1,k+1) +
     @                              signd(i+1,j,k+1,2)   * dest(i+1,j,k+1) +
     @                              signd(i+1,j+1,k+1,2) * dest(i+1,j+2,k+1) +
     @                              signd(i+1,j+1,k,3)   * dest(i+1,j+1,k) +
     @                              signd(i+1,j+1,k+1,3) * dest(i+1,j+1,k+2)) /
     @                             (signd(i,j+1,k+1,1) + signd(i+1,j+1,k+1,1) +
     @                              signd(i+1,j,k+1,2) + signd(i+1,j+1,k+1,2) +
     @                              signd(i+1,j+1,k,3) + signd(i+1,j+1,k+1,3))
      endif
      end

#  endif

#  ifndef HG_SIGMA_NODE_ALL

c CELL-based data only.
      subroutine hgsrst(destx, desty, destz,
     @ dims(dest),
     @ dims(reg),
     @ srcx, srcy, srcz,
     @ dims(src), ir, jr, kr)
      integer dims(dest)
      integer dims(reg)
      integer dims(src)
      integer ir, jr, kr
      real*8 destx(dimdec(dest))
      real*8 desty(dimdec(dest))
      real*8 destz(dimdec(dest))
      real*8 srcx(dimdec(src))
      real*8 srcy(dimdec(src))
      real*8 srcz(dimdec(src))
      integer i, j, k, i2, j2, k2
      if (ir .eq. 2 .and. jr .eq. 2 .and. kr .eq. 2) then
         do 10 k = regl2, regh2
            do 10 j = regl1, regh1
               do 10 i = regl0, regh0
                  i2 = 2 * i
                  j2 = 2 * j
                  k2 = 2 * k
                  destx(i,j,k) = 0.5d0 / (1.d0 / (srcx(i2,j2,k2) +
     @                                            srcx(i2,j2,k2+1) +
     @                                            srcx(i2,j2+1,k2) +
     @                                            srcx(i2,j2+1,k2+1)) +
     @                                    1.d0 / (srcx(i2+1,j2,k2) +
     @                                            srcx(i2+1,j2,k2+1) +
     @                                            srcx(i2+1,j2+1,k2) +
     @                                            srcx(i2+1,j2+1,k2+1)))
                  desty(i,j,k) = 0.5d0 / (1.d0 / (srcy(i2,j2,k2) +
     @                                            srcy(i2,j2,k2+1) +
     @                                            srcy(i2+1,j2,k2) +
     @                                            srcy(i2+1,j2,k2+1)) +
     @                                    1.d0 / (srcy(i2,j2+1,k2) +
     @                                            srcy(i2,j2+1,k2+1) +
     @                                            srcy(i2+1,j2+1,k2) +
     @                                            srcy(i2+1,j2+1,k2+1)))
                  destz(i,j,k) = 0.5d0 / (1.d0 / (srcz(i2,j2,k2) +
     @                                            srcz(i2,j2+1,k2) +
     @                                            srcz(i2+1,j2,k2) +
     @                                            srcz(i2+1,j2+1,k2)) +
     @                                    1.d0 / (srcz(i2,j2,k2+1) +
     @                                            srcz(i2,j2+1,k2+1) +
     @                                            srcz(i2+1,j2,k2+1) +
     @                                            srcz(i2+1,j2+1,k2+1)))
 10               continue
      else if (ir .eq. 2 .and. jr .eq. 2) then
         do 20 k = regl2, regh2
            do 20 j = regl1, regh1
               do 20 i = regl0, regh0
                  i2 = 2 * i
                  j2 = 2 * j
                  destx(i,j,k) = 1.d0 / (1.d0 / (srcx(i2,j2,k) +
     @                                           srcx(i2,j2+1,k)) +
     @                                   1.d0 / (srcx(i2+1,j2,k) +
     @                                           srcx(i2+1,j2+1,k)))
                  desty(i,j,k) = 1.d0 / (1.d0 / (srcy(i2,j2,k) +
     @                                           srcy(i2+1,j2,k)) +
     @                                   1.d0 / (srcy(i2,j2+1,k) +
     @                                           srcy(i2+1,j2+1,k)))
                  destz(i,j,k) = 0.25d0 * (srcz(i2,j2,k) +
     @                                     srcz(i2,j2+1,k) +
     @                                     srcz(i2+1,j2,k) +
     @                                     srcz(i2+1,j2+1,k))
 20               continue
      else if (ir .eq. 2 .and. kr .eq. 2) then
         do 30 k = regl2, regh2
            do 30 j = regl1, regh1
               do 30 i = regl0, regh0
                  i2 = 2 * i
                  k2 = 2 * k
                  destx(i,j,k) = 1.d0 / (1.d0 / (srcx(i2,j,k2) +
     @                                           srcx(i2,j,k2+1)) +
     @                                   1.d0 / (srcx(i2+1,j,k2) +
     @                                           srcx(i2+1,j,k2+1)))
                  desty(i,j,k) = 0.25d0 * (srcy(i2,j,k2) +
     @                                     srcy(i2,j,k2+1) +
     @                                     srcy(i2+1,j,k2) +
     @                                     srcy(i2+1,j,k2+1))
                  destz(i,j,k) = 1.d0 / (1.d0 / (srcz(i2,j,k2) +
     @                                           srcz(i2+1,j,k2)) +
     @                                   1.d0 / (srcz(i2,j,k2+1) +
     @                                           srcz(i2+1,j,k2+1)))
 30               continue
      else if (jr .eq. 2 .and. kr .eq. 2) then
         do 40 k = regl2, regh2
            do 40 j = regl1, regh1
               do 40 i = regl0, regh0
                  j2 = 2 * j
                  k2 = 2 * k
                  destx(i,j,k) = 0.25d0 * (srcx(i,j2,k2) +
     @                                     srcx(i,j2,k2+1) +
     @                                     srcx(i,j2+1,k2) +
     @                                     srcx(i,j2+1,k2+1))
                  desty(i,j,k) = 1.d0 / (1.d0 / (srcy(i,j2,k2) +
     @                                           srcy(i,j2,k2+1)) +
     @                                   1.d0 / (srcy(i,j2+1,k2) +
     @                                           srcy(i,j2+1,k2+1)))
                  destz(i,j,k) = 1.d0 / (1.d0 / (srcz(i,j2,k2) +
     @                                           srcz(i,j2+1,k2)) +
     @                                   1.d0 / (srcz(i,j2,k2+1) +
     @                                           srcz(i,j2+1,k2+1)))
 40               continue
      else if (ir .eq. 2) then
         do 50 k = regl2, regh2
            do 50 j = regl1, regh1
               do 50 i = regl0, regh0
                  i2 = 2 * i
                  destx(i,j,k) = 2.d0 / (1.d0 / srcx(i2,j,k) +
     @                                   1.d0 / srcx(i2+1,j,k))
                  desty(i,j,k) = 0.5d0 * (srcy(i2,j,k) +
     @                                    srcy(i2+1,j,k))
                  destz(i,j,k) = 0.5d0 * (srcz(i2,j,k) +
     @                                    srcz(i2+1,j,k))
 50               continue
      else if (jr .eq. 2) then
         do 60 k = regl2, regh2
            do 60 j = regl1, regh1
               do 60 i = regl0, regh0
                  j2 = 2 * j
                  destx(i,j,k) = 0.5d0 * (srcx(i,j2,k) +
     @                                    srcx(i,j2+1,k))
                  desty(i,j,k) = 2.d0 / (1.d0 / srcy(i,j2,k) +
     @                                   1.d0 / srcy(i,j2+1,k))
                  destz(i,j,k) = 0.5d0 * (srcz(i,j2,k) +
     @                                    srcz(i,j2+1,k))
 60               continue
      else if (kr .eq. 2) then
         do 70 k = regl2, regh2
            do 70 j = regl1, regh1
               do 70 i = regl0, regh0
                  k2 = 2 * k
                  destx(i,j,k) = 0.5d0 * (srcx(i,j,k2) +
     @                                    srcx(i,j,k2+1))
                  desty(i,j,k) = 0.5d0 * (srcy(i,j,k2) +
     @                                    srcy(i,j,k2+1))
                  destz(i,j,k) = 2.d0 / (1.d0 / srcz(i,j,k2) +
     @                                   1.d0 / srcz(i,j,k2+1))
 70               continue
      endif
      end

#  else

c NODE(EDGE)-based data only.
      subroutine hgsrst(destx, desty, destz,
     @ dims(dest),
     @ dims(reg),
     @ srcx, srcy, srcz,
     @ dims(src))
      implicit real*8(a-h,o-z), integer(i-n)
      integer dims(dest)
      integer dims(reg)
      integer dims(src)
      real*8 destx(dimdec(dest))
      real*8 desty(dimdec(dest))
      real*8 destz(dimdec(dest))
      real*8 srcx(dimdec(src))
      real*8 srcy(dimdec(src))
      real*8 srcz(dimdec(src))
      print *, "Warning:  Don't use this routine.  Not modified for "
      print *, "          IntVect refinement ratios"
      do 10 k = regl2, regh2
         do 10 j = regl1, regh1
            do 10 i = regl0, regh0
               i2 = 2 * i
               j2 = 2 * j
               k2 = 2 * k
c               destx(i,j,k) = 0.5d0 / (1.d0 / srcx(i2,j2,k2) +
c     @                                 1.d0 / srcx(i2+1,j2,k2)) +
c     @                        0.5d0 / (1.d0 / srcx(i2,j2,k2+1) +
c     @                                 1.d0 / srcx(i2+1,j2,k2+1)) +
c     @                        0.5d0 / (1.d0 / srcx(i2,j2+1,k2) +
c     @                                 1.d0 / srcx(i2+1,j2+1,k2)) +
c     @                        0.5d0 / (1.d0 / srcx(i2,j2+1,k2+1) +
c     @                                 1.d0 / srcx(i2+1,j2+1,k2+1))
c               desty(i,j,k) = 0.5d0 / (1.d0 / srcy(i2,j2,k2) +
c     @                                 1.d0 / srcy(i2,j2+1,k2)) +
c     @                        0.5d0 / (1.d0 / srcy(i2,j2,k2+1) +
c     @                                 1.d0 / srcy(i2,j2+1,k2+1)) +
c     @                        0.5d0 / (1.d0 / srcy(i2+1,j2,k2) +
c     @                                 1.d0 / srcy(i2+1,j2+1,k2)) +
c     @                        0.5d0 / (1.d0 / srcy(i2+1,j2,k2+1) +
c     @                                 1.d0 / srcy(i2+1,j2+1,k2+1))
c               destz(i,j,k) = 0.5d0 / (1.d0 / srcz(i2,j2,k2) +
c     @                                 1.d0 / srcz(i2,j2,k2+1)) +
c     @                        0.5d0 / (1.d0 / srcz(i2+1,j2,k2) +
c     @                                 1.d0 / srcz(i2+1,j2,k2+1)) +
c     @                        0.5d0 / (1.d0 / srcz(i2,j2+1,k2) +
c     @                                 1.d0 / srcz(i2,j2+1,k2+1)) +
c     @                        0.5d0 / (1.d0 / srcz(i2+1,j2+1,k2) +
c     @                                 1.d0 / srcz(i2+1,j2+1,k2+1))
               destx(i,j,k) = 0.5d0 / (1.d0 / (srcx(i2,j2,k2) +
     @                                         srcx(i2,j2,k2+1) +
     @                                         srcx(i2,j2+1,k2) +
     @                                         srcx(i2,j2+1,k2+1)) +
     @                                 1.d0 / (srcx(i2+1,j2,k2) +
     @                                         srcx(i2+1,j2,k2+1) +
     @                                         srcx(i2+1,j2+1,k2) +
     @                                         srcx(i2+1,j2+1,k2+1)))
               desty(i,j,k) = 0.5d0 / (1.d0 / (srcy(i2,j2,k2) +
     @                                         srcy(i2,j2,k2+1) +
     @                                         srcy(i2+1,j2,k2) +
     @                                         srcy(i2+1,j2,k2+1)) +
     @                                 1.d0 / (srcy(i2,j2+1,k2) +
     @                                         srcy(i2,j2+1,k2+1) +
     @                                         srcy(i2+1,j2+1,k2) +
     @                                         srcy(i2+1,j2+1,k2+1)))
               destz(i,j,k) = 0.5d0 / (1.d0 / (srcz(i2,j2,k2) +
     @                                         srcz(i2,j2+1,k2) +
     @                                         srcz(i2+1,j2,k2) +
     @                                         srcz(i2+1,j2+1,k2)) +
     @                                 1.d0 / (srcz(i2,j2,k2+1) +
     @                                         srcz(i2,j2+1,k2+1) +
     @                                         srcz(i2+1,j2,k2+1) +
     @                                         srcz(i2+1,j2+1,k2+1)))
 10         continue
      end

#  endif

#endif

#ifdef HG_CROSS_STENCIL

#ifdef HG_CONSTANT

#define AVG fac * (cor(i-1) + cor(i+1) + \
                   cor(i-jdiff) + cor(i+jdiff) + \
                   cor(i-kdiff) + cor(i+kdiff))

      subroutine hgrlxu(cor,
     @ res,
     @ dims(res),
     @ mask,
     @ dims(reg),
     @ hx)
      implicit real*8(a-h,o-z), integer(i-n)
      integer dims(res)
      integer dims(reg)
      real*8 cor(dimdecu(res))
      real*8 res(dimdecu(res))
      real*8 mask(dimdecu(res))
      real*8 hx
#if 0
      avg() = fac * (cor(i-1) + cor(i+1) +
     @               cor(i-jdiff) + cor(i+jdiff) +
     @               cor(i-kdiff) + cor(i+kdiff))
#endif
      hxm2 = 1.d0 / (hx*hx)
      fac = hxm2
      cen = 1.d0 / (6.d0 * hxm2)
      jdiff = resh0 - resl0 + 1
      kdiff = (resh1 - resl1 + 1) * jdiff
cdir$ ivdep
      do 10 i = (regl2 - resl2) * kdiff + (regl1 - resl1) * jdiff +
     @          (regl0 - resl0) + 1,
     @          (regh2 - resl2) * kdiff + (regh1 - resl1) * jdiff +
     @          (regh0 - resl0) + 1, 2
 10      cor(i) = cor(i) + mask(i) * ((AVG - res(i)) * cen - cor(i))
cdir$ ivdep
      do 20 i = (regl2 - resl2) * kdiff + (regl1 - resl1) * jdiff +
     @          (regl0 - resl0) + 2,
     @          (regh2 - resl2) * kdiff + (regh1 - resl1) * jdiff +
     @          (regh0 - resl0) + 1, 2
 20      cor(i) = cor(i) + mask(i) * ((AVG - res(i)) * cen - cor(i))
      end

      subroutine hgres(res,
     @ dims(res),
     @ src, dest,
     @ dims(reg),
     @ hx)
      implicit real*8(a-h,o-z), integer(i-n)
      integer dims(res)
      integer dims(reg)
      real*8 res(dimdec(res))
      real*8 src(dimdec(res))
      real*8 dest(dimdec(res))
      real*8 hx
      hxm2 = 1.d0 / (hx*hx)
      fac = hxm2
      if (regh2 - regl2 .gt. regh1 - regl1 .and.
     @    regh2 - regl2 .gt. regh0 - regl0) then
         do 10 j = regl1, regh1
            do 10 i = regl0, regh0
               do 10 k = regl2, regh2
 10               res(i,j,k) = src(i,j,k) - fac *
     @           (dest(i-1,j,k) + dest(i+1,j,k) +
     @            dest(i,j-1,k) + dest(i,j+1,k) +
     @            dest(i,j,k-1) + dest(i,j,k+1) -
     @            6.d0 * dest(i,j,k))
      else if (regh1 - regl1 .gt. regh0 - regl0) then
         do 20 k = regl2, regh2
            do 20 i = regl0, regh0
               do 20 j = regl1, regh1
 20               res(i,j,k) = src(i,j,k) - fac *
     @           (dest(i-1,j,k) + dest(i+1,j,k) +
     @            dest(i,j-1,k) + dest(i,j+1,k) +
     @            dest(i,j,k-1) + dest(i,j,k+1) -
     @            6.d0 * dest(i,j,k))
      else
         do 30 k = regl2, regh2
            do 30 j = regl1, regh1
               do 30 i = regl0, regh0
 30               res(i,j,k) = src(i,j,k) - fac *
     @           (dest(i-1,j,k) + dest(i+1,j,k) +
     @            dest(i,j-1,k) + dest(i,j+1,k) +
     @            dest(i,j,k-1) + dest(i,j,k+1) -
     @            6.d0 * dest(i,j,k))
      endif
      end

      subroutine hgresu(res,
     @ dims(res),
     @ src, dest,
     @ dims(reg),
     @ hx)
      implicit real*8(a-h,o-z), integer(i-n)
      integer dims(res)
      integer dims(reg)
      real*8 res(dimdecu(res))
      real*8 src(dimdecu(res))
      real*8 dest(dimdecu(res))
      real*8 hx
      hxm2 = 1.d0 / (hx*hx)
      fac = hxm2
      jdiff = resh0 - resl0 + 1
      kdiff = (resh1 - resl1 + 1) * jdiff
      do 10 i = (regl2 - resl2) * kdiff + (regl1 - resl1) * jdiff +
     @          (regl0 - resl0) + 1,
     @          (regh2 - resl2) * kdiff + (regh1 - resl1) * jdiff +
     @          (regh0 - resl0) + 1
 10      res(i) = src(i) - fac *
     @     (dest(i-1) + dest(i+1) +
     @      dest(i-jdiff) + dest(i+jdiff) +
     @      dest(i-kdiff) + dest(i+kdiff) -
     @      6.d0 * dest(i))
      end

#else

c seven-point variable stencils

#ifdef HG_SIGMA_NODE

      subroutine hgcen(cen,
     @ dims(cen),
     @ signd,
     @ dims(sn),
     @ dims(reg))
      integer dims(cen)
      integer dims(sn)
      integer dims(reg)
      real*8 cen(dimdec(cen))
      real*8 signd(dimdec(sn), 3)
      integer i, j, k
      do 10 k = regl2, regh2
         do 10 j = regl1, regh1
            do 10 i = regl0, regh0
 10            cen(i,j,k) = 1.d0 / (signd(i-1,j,k,1) + signd(i,j,k,1) +
     @                              signd(i,j-1,k,2) + signd(i,j,k,2) +
     @                              signd(i,j,k-1,3) + signd(i,j,k,3))
      end

#define AVGU (sig(i-1)        * cor(i-1) + \
              sig(i)          * cor(i+1) + \
              sig(i+ly-jdiff) * cor(i-jdiff) + \
              sig(i+ly)       * cor(i+jdiff) + \
              sig(i+lz-kdiff) * cor(i-kdiff) + \
              sig(i+lz)       * cor(i+kdiff))

c sig here contains three different directions all stored on "nodes"
      subroutine hgrlxu(cor, res, sig, cen,
     @ dims(res),
     @ mask,
     @ dims(reg))
      integer dims(res)
      integer dims(reg)
      real*8 cor(dimdecu(res))
      real*8 res(dimdecu(res))
      real*8 sig(dimdecu(res))
      real*8 cen(dimdecu(res))
      real*8 mask(dimdecu(res))
      integer i, jdiff, kdiff, ly, lz
#if 0
      avg() = (sig(i-1)        * cor(i-1) +
     @         sig(i)          * cor(i+1) +
     @         sig(i+ly-jdiff) * cor(i-jdiff) +
     @         sig(i+ly)       * cor(i+jdiff) +
     @         sig(i+lz-kdiff) * cor(i-kdiff) +
     @         sig(i+lz)       * cor(i+kdiff))
#endif
      jdiff =  resh0 - resl0 + 1
      kdiff = (resh1 - resl1 + 1) * jdiff
      ly    = (resh2 - resl2 + 1) * kdiff
      lz    = 2 * ly
cdir$ ivdep
      do 10 i = (regl2 - resl2) * kdiff + (regl1 - resl1) * jdiff +
     @          (regl0 - resl0) + 1,
     @          (regh2 - resl2) * kdiff + (regh1 - resl1) * jdiff +
     @          (regh0 - resl0) + 1, 2
 10   cor(i) = cor(i) + mask(i) * ((AVGU - res(i)) * cen(i) - cor(i))
cdir$ ivdep
      do 20 i = (regl2 - resl2) * kdiff + (regl1 - resl1) * jdiff +
     @          (regl0 - resl0) + 2,
     @          (regh2 - resl2) * kdiff + (regh1 - resl1) * jdiff +
     @          (regh0 - resl0) + 1, 2
 20   cor(i) = cor(i) + mask(i) * ((AVGU - res(i)) * cen(i) - cor(i))
      end

#define AVG (signd(i-1,j,k,1) * cor(i-1,j,k) + \
             signd(i,j,k,1)   * cor(i+1,j,k) + \
             signd(i,j-1,k,2) * cor(i,j-1,k) + \
             signd(i,j,k,2)   * cor(i,j+1,k) + \
             signd(i,j,k-1,3) * cor(i,j,k-1) + \
             signd(i,j,k,3)   * cor(i,j,k+1))

      subroutine hgrlx(cor,
     @ dims(cor),
     @ res,
     @ dims(res),
     @ signd,
     @ dims(sn),
     @ cen,
     @ dims(cen),
     @ dims(reg))
      integer dims(cor)
      integer dims(res)
      integer dims(sn)
      integer dims(cen)
      integer dims(reg)
      real*8 cor(dimdec(cor))
      real*8 res(dimdec(res))
      real*8 signd(dimdec(sn), 3)
      real*8 cen(dimdec(cen))
      integer i, j, k, ipar0, ipar, ipass
#if 0
      avg() = (signd(i-1,j,k,1) * cor(i-1,j,k) +
     @         signd(i,j,k,1)   * cor(i+1,j,k) +
     @         signd(i,j-1,k,2) * cor(i,j-1,k) +
     @         signd(i,j,k,2)   * cor(i,j+1,k) +
     @         signd(i,j,k-1,3) * cor(i,j,k-1) +
     @         signd(i,j,k,3)   * cor(i,j,k+1))
#endif
      if (regh2 - regl2 .gt. regh1 - regl1 .and.
     @    regh2 - regl2 .gt. regh0 - regl0) then
         do 10 ipass = 0, 1
            ipar0 = ipass
            do 10 j = regl1, regh1
               ipar0 = 1 - ipar0
               ipar = ipar0
               do 10 i = regl0, regh0
                  ipar = 1 - ipar
cdir$ ivdep
                  do 10 k = regl2 + ipar, regh2, 2
                     cor(i,j,k) = (AVG - res(i,j,k)) * cen(i,j,k)
 10                  continue
      else if (regh1 - regl1 .gt. regh0 - regl0) then
         do 20 ipass = 0, 1
            ipar0 = ipass
            do 20 k = regl2, regh2
               ipar0 = 1 - ipar0
               ipar = ipar0
               do 20 i = regl0, regh0
                  ipar = 1 - ipar
cdir$ ivdep
                  do 20 j = regl1 + ipar, regh1, 2
                     cor(i,j,k) = (AVG - res(i,j,k)) * cen(i,j,k)
 20                  continue
      else
         do 30 ipass = 0, 1
            ipar0 = ipass
            do 30 k = regl2, regh2
               ipar0 = 1 - ipar0
               ipar = ipar0
               do 30 j = regl1, regh1
                  ipar = 1 - ipar
cdir$ ivdep
                  do 30 i = regl0 + ipar, regh0, 2
                     cor(i,j,k) = (AVG - res(i,j,k)) * cen(i,j,k)
 30                  continue
      endif
      end

#define RHSL0 (res(i,j,k) - (signd(i,j-1,k,2) * cor(i,j-1,k) + \
                             signd(i,j,k,2)   * cor(i,j+1,k) + \
                             signd(i,j,k-1,3) * cor(i,j,k-1) + \
                             signd(i,j,k,3)   * cor(i,j,k+1)))

#define RHSL1 (res(i,j,k) - (signd(i-1,j,k,1) * cor(i-1,j,k) + \
                             signd(i,j,k,1)   * cor(i+1,j,k) + \
                             signd(i,j,k-1,3) * cor(i,j,k-1) + \
                             signd(i,j,k,3)   * cor(i,j,k+1)))

#define RHSL2 (res(i,j,k) - (signd(i-1,j,k,1) * cor(i-1,j,k) + \
                             signd(i,j,k,1)   * cor(i+1,j,k) + \
                             signd(i,j-1,k,2) * cor(i,j-1,k) + \
                             signd(i,j,k,2)   * cor(i,j+1,k)))

      subroutine hgrlxl(cor,
     @ dims(cor),
     @ res,
     @ dims(res),
     @ signd,
     @ dims(sn),
     @ cen,
     @ dims(cen),
     @ dims(reg),
     @ dims(dom),
     @ lsd)
      integer dims(cor)
      integer dims(res)
      integer dims(sn)
      integer dims(cen)
      integer dims(reg)
      integer dims(dom)
      real*8 cor(dimdec(cor))
      real*8 res(dimdec(res))
      real*8 signd(dimdec(sn), 3)
      real*8 cen(dimdec(cen))
      real*8 wrk(256)
      integer lsd
      real*8 betm, aj
      integer i, j, k, kw, ipass, ipar
      if (lsd .eq. 2) then
         do 20 ipass = 1, 0, -1
         ipar = ipass
         do 20 i = regl0, regh0
         ipar = 1 - ipar
         do 20 j = regl1 + ipar, regh1, 2
            k = regl2
            betm = -cen(i,j,k)
            if (betm .eq. 0.d0) then
c dirichlet bdy:
               cor(i,j,k) = 0.d0
               wrk(1) = 0.d0
            else if (regl2 .eq. doml2) then
c neumann bdy:
               cor(i,j,k) = RHSL2 * betm
               aj = signd(i,j,k,3)
               wrk(1) = 2.d0 * aj * betm
            else
c interface to grid at same level:
               aj = signd(i,j,k-1,3)
               cor(i,j,k) = (RHSL2 - aj * cor(i,j,k-1)) * betm
               aj = signd(i,j,k,3)
               wrk(1) = aj * betm
            endif
c forward solve loop:
            do 10 k = regl2 + 1, regh2 - 1
               kw = k - regl2
               if (cen(i,j,k) .eq. 0.d0) then
                  betm = 0.d0
               else
                  betm = 1.d0 / (-1.d0 / cen(i,j,k) - aj * wrk(kw))
               endif
               cor(i,j,k) = (RHSL2 - aj * cor(i,j,k-1)) * betm
               aj = signd(i,j,k,3)
               wrk(kw + 1) = aj * betm
 10            continue
            k = regh2
            kw = k - regl2
            if (cen(i,j,k) .eq. 0.d0) then
c dirichlet bdy:
               cor(i,j,k) = 0.d0
            else if (regh2 .eq. domh2) then
c neumann bdy:
               aj = 2.d0 * aj
               betm = 1.d0 / (-1.d0 / cen(i,j,k) - aj * wrk(kw))
               cor(i,j,k) = (RHSL2 - aj * cor(i,j,k-1)) * betm
            else if (kw .gt. 0) then
c interface to grid at same level:
               betm = 1.d0 / (-1.d0 / cen(i,j,k) - aj * wrk(kw))
               cor(i,j,k) = RHSL2 - aj * cor(i,j,k-1)
               aj = signd(i,j,k,3)
               cor(i,j,k) = (cor(i,j,k) - aj * cor(i,j,k+1)) * betm
            endif
c back substitution loop:
            do 20 k = regh2 - 1, regl2, -1
               kw = k - regl2
               cor(i,j,k) = cor(i,j,k) - wrk(kw + 1) * cor(i,j,k+1)
 20            continue
      else
         print *, "Line solve not implemented in dimension", lsd
      endif
      end

      subroutine hgrlnf(cor,
     @ dims(cor),
     @ res, dims(res),
     @ wrk, dims(wrk),
     @ signd, dims(sn),
     @ cen, dims(cen),
     @ dims(reg),
     @ dims(dom),
     @ lsd, ipass)
      integer dims(cor)
      integer dims(res)
      integer dims(wrk)
      integer dims(sn)
      integer dims(cen)
      integer dims(reg)
      integer dims(dom)
      real*8 cor(dimdec(cor))
      real*8 res(dimdec(res))
      real*8 wrk(dimdec(wrk))
      real*8 signd(dimdec(sn), 3)
      real*8 cen(dimdec(cen))
      integer lsd, ipass
      real*8 betm, aj
      integer i, j, k, ioff, ipar
      if (lsd .eq. 0) then
         if (mod(regl1 + regl2, 2) .eq. 0) then
            ioff = 1 - ipass
         else
            ioff = ipass
         endif
         i = regl0
         ipar = ioff
         do 10 k = regl2, regh2
            ipar = 1 - ipar
            do 10 j = regl1 + ipar, regh1, 2
               if (cen(i,j,k) .eq. 0.d0) then
c dirichlet bdy:
                  wrk(i,j,k) = 0.d0
               else if (regl0 .eq. doml0) then
c neumann bdy:
                  betm = -cen(i,j,k)
                  cor(i,j,k) = RHSL0 * betm
                  aj = signd(i,j,k,1)
                  wrk(i,j,k) = 2.d0 * aj * betm
               endif
 10            continue
c forward solve loop:
         do 20 i = regl0 + 1, regh0 - 1
            ipar = ioff
            do 20 k = regl2, regh2
               ipar = 1 - ipar
               do 20 j = regl1 + ipar, regh1, 2
                  aj = signd(i-1,j,k,1)
                  if (cen(i,j,k) .eq. 0.d0) then
                     betm = 0.d0
                  else
                     betm = 1.d0 / (-1.d0 / cen(i,j,k) - aj * wrk(i-1,j,k))
                  endif
                  cor(i,j,k) = (RHSL0 - aj * cor(i-1,j,k)) * betm
                  aj = signd(i,j,k,1)
 20               wrk(i,j,k) = aj * betm
         i = regh0
         ipar = ioff
         do 30 k = regl2, regh2
            ipar = 1 - ipar
            do 30 j = regl1 + ipar, regh1, 2
               if (cen(i,j,k) .eq. 0.d0) then
c dirichlet bdy:
               else if (regh0 .eq. domh0) then
c neumann bdy:
                  aj = 2.d0 * signd(i-1,j,k,1)
                  betm = 1.d0 / (-1.d0 / cen(i,j,k) - aj * wrk(i-1,j,k))
                  cor(i,j,k) = (RHSL0 - aj * cor(i-1,j,k)) * betm
               else if (i .gt. regl0) then
c interface to grid at same level:
                  aj = signd(i-1,j,k,1)
                  betm = 1.d0 / (-1.d0 / cen(i,j,k) - aj * wrk(i-1,j,k))
                  cor(i,j,k) = (RHSL0 - aj * cor(i-1,j,k)) * betm
                  aj = signd(i,j,k,1)
                  wrk(i,j,k) = aj * betm
               endif
 30            continue
      else if (lsd .eq. 1) then
         if (mod(regl0 + regl2, 2) .eq. 0) then
            ioff = 1 - ipass
         else
            ioff = ipass
         endif
         j = regl1
         ipar = ioff
         do 110 k = regl2, regh2
            ipar = 1 - ipar
            do 110 i = regl0 + ipar, regh0, 2
               if (cen(i,j,k) .eq. 0.d0) then
c dirichlet bdy:
                  wrk(i,j,k) = 0.d0
               else if (regl1 .eq. doml1) then
c neumann bdy:
                  betm = -cen(i,j,k)
                  cor(i,j,k) = RHSL1 * betm
                  aj = signd(i,j,k,2)
                  wrk(i,j,k) = 2.d0 * aj * betm
               endif
 110           continue
c forward solve loop:
         do 120 j = regl1 + 1, regh1 - 1
            ipar = ioff
            do 120 k = regl2, regh2
               ipar = 1 - ipar
               do 120 i = regl0 + ipar, regh0, 2
                  aj = signd(i,j-1,k,2)
                  if (cen(i,j,k) .eq. 0.d0) then
                     betm = 0.d0
                  else
                     betm = 1.d0 / (-1.d0 / cen(i,j,k) - aj * wrk(i,j-1,k))
                  endif
                  cor(i,j,k) = (RHSL1 - aj * cor(i,j-1,k)) * betm
                  aj = signd(i,j,k,2)
 120              wrk(i,j,k) = aj * betm
         j = regh1
         ipar = ioff
         do 130 k = regl2, regh2
            ipar = 1 - ipar
            do 130 i = regl0 + ipar, regh0, 2
               if (cen(i,j,k) .eq. 0.d0) then
c dirichlet bdy:
               else if (regh1 .eq. domh1) then
c neumann bdy:
                  aj = 2.d0 * signd(i,j-1,k,2)
                  betm = 1.d0 / (-1.d0 / cen(i,j,k) - aj * wrk(i,j-1,k))
                  cor(i,j,k) = (RHSL1 - aj * cor(i,j-1,k)) * betm
               else if (j .gt. regl1) then
c interface to grid at same level:
                  aj = signd(i,j-1,k,2)
                  betm = 1.d0 / (-1.d0 / cen(i,j,k) - aj * wrk(i,j-1,k))
                  cor(i,j,k) = (RHSL1 - aj * cor(i,j-1,k)) * betm
                  aj = signd(i,j,k,2)
                  wrk(i,j,k) = aj * betm
               endif
 130           continue
      else
         if (mod(regl0 + regl1, 2) .eq. 0) then
            ioff = 1 - ipass
         else
            ioff = ipass
         endif
         k = regl2
         ipar = ioff
         do 210 j = regl1, regh1
            ipar = 1 - ipar
            do 210 i = regl0 + ipar, regh0, 2
               if (cen(i,j,k) .eq. 0.d0) then
c dirichlet bdy:
                  wrk(i,j,k) = 0.d0
               else if (regl2 .eq. doml2) then
c neumann bdy:
                  betm = -cen(i,j,k)
                  cor(i,j,k) = RHSL2 * betm
                  aj = signd(i,j,k,3)
                  wrk(i,j,k) = 2.d0 * aj * betm
               endif
 210           continue
c forward solve loop:
         do 220 k = regl2 + 1, regh2 - 1
            ipar = ioff
            do 220 j = regl1, regh1
               ipar = 1 - ipar
               do 220 i = regl0 + ipar, regh0, 2
                  aj = signd(i,j,k-1,3)
                  if (cen(i,j,k) .eq. 0.d0) then
                     betm = 0.d0
                  else
                     betm = 1.d0 / (-1.d0 / cen(i,j,k) - aj * wrk(i,j,k-1))
                  endif
                  cor(i,j,k) = (RHSL2 - aj * cor(i,j,k-1)) * betm
                  aj = signd(i,j,k,3)
 220              wrk(i,j,k) = aj * betm
         k = regh2
         ipar = ioff
         do 230 j = regl1, regh1
            ipar = 1 - ipar
            do 230 i = regl0 + ipar, regh0, 2
               if (cen(i,j,k) .eq. 0.d0) then
c dirichlet bdy:
               else if (regh2 .eq. domh2) then
c neumann bdy:
                  aj = 2.d0 * signd(i,j,k-1,3)
                  betm = 1.d0 / (-1.d0 / cen(i,j,k) - aj * wrk(i,j,k-1))
                  cor(i,j,k) = (RHSL2 - aj * cor(i,j,k-1)) * betm
               else if (k .gt. regl2) then
c interface to grid at same level:
                  aj = signd(i,j,k-1,3)
                  betm = 1.d0 / (-1.d0 / cen(i,j,k) - aj * wrk(i,j,k-1))
                  cor(i,j,k) = (RHSL2 - aj * cor(i,j,k-1)) * betm
                  aj = signd(i,j,k,3)
                  wrk(i,j,k) = aj * betm
               endif
 230           continue
      endif
      end

      subroutine hgres(res,
     @ dims(res),
     @ src,
     @ dims(src),
     @ dest,
     @ dims(dest),
     @ signd,
     @ dims(sn),
     @ dims(reg))
      integer dims(res)
      integer dims(src)
      integer dims(dest)
      integer dims(sn)
      integer dims(reg)
      real*8 res(dimdec(res))
      real*8 src(dimdec(src))
      real*8 dest(dimdec(dest))
      real*8 signd(dimdec(sn), 3)
      integer i, j, k
      if (regh2 - regl2 .gt. regh1 - regl1 .and.
     @    regh2 - regl2 .gt. regh0 - regl0) then
         do 10 j = regl1, regh1
            do 10 i = regl0, regh0
               do 10 k = regl2, regh2
 10               res(i,j,k) = src(i,j,k) -
     @              (signd(i-1,j,k,1) * (dest(i-1,j,k) - dest(i,j,k)) +
     @               signd(i,j,k,1)   * (dest(i+1,j,k) - dest(i,j,k)) +
     @               signd(i,j-1,k,2) * (dest(i,j-1,k) - dest(i,j,k)) +
     @               signd(i,j,k,2)   * (dest(i,j+1,k) - dest(i,j,k)) +
     @               signd(i,j,k-1,3) * (dest(i,j,k-1) - dest(i,j,k)) +
     @               signd(i,j,k,3)   * (dest(i,j,k+1) - dest(i,j,k)))
      else if (regh1 - regl1 .gt. regh0 - regl0) then
         do 20 k = regl2, regh2
            do 20 i = regl0, regh0
               do 20 j = regl1, regh1
 20               res(i,j,k) = src(i,j,k) -
     @              (signd(i-1,j,k,1) * (dest(i-1,j,k) - dest(i,j,k)) +
     @               signd(i,j,k,1)   * (dest(i+1,j,k) - dest(i,j,k)) +
     @               signd(i,j-1,k,2) * (dest(i,j-1,k) - dest(i,j,k)) +
     @               signd(i,j,k,2)   * (dest(i,j+1,k) - dest(i,j,k)) +
     @               signd(i,j,k-1,3) * (dest(i,j,k-1) - dest(i,j,k)) +
     @               signd(i,j,k,3)   * (dest(i,j,k+1) - dest(i,j,k)))
      else
         do 30 k = regl2, regh2
            do 30 j = regl1, regh1
               do 30 i = regl0, regh0
 30               res(i,j,k) = src(i,j,k) -
     @              (signd(i-1,j,k,1) * (dest(i-1,j,k) - dest(i,j,k)) +
     @               signd(i,j,k,1)   * (dest(i+1,j,k) - dest(i,j,k)) +
     @               signd(i,j-1,k,2) * (dest(i,j-1,k) - dest(i,j,k)) +
     @               signd(i,j,k,2)   * (dest(i,j+1,k) - dest(i,j,k)) +
     @               signd(i,j,k-1,3) * (dest(i,j,k-1) - dest(i,j,k)) +
     @               signd(i,j,k,3)   * (dest(i,j,k+1) - dest(i,j,k)))
      endif
      end

      subroutine hgresu(res,
     @ dims(res),
     @ src, dest, signd, mask,
     @ dims(reg))
      integer dims(res)
      integer dims(reg)
      real*8 res(dimdecu(res))
      real*8 src(dimdecu(res))
      real*8 dest(dimdecu(res))
      real*8 signd(dimdecu(res))
      real*8 mask(dimdecu(res))
      integer i, jdiff, kdiff, ly, lz
      jdiff = resh0 - resl0 + 1
      kdiff = (resh1 - resl1 + 1) * jdiff
      ly    = (resh2 - resl2 + 1) * kdiff
      lz    = 2 * ly
      do 10 i = (regl2 - resl2) * kdiff + (regl1 - resl1) * jdiff +
     @          (regl0 - resl0) + 1,
     @          (regh2 - resl2) * kdiff + (regh1 - resl1) * jdiff +
     @          (regh0 - resl0) + 1
 10      res(i) = mask(i) * (src(i) -
     @     (signd(i-1)        * (dest(i-1) - dest(i)) +
     @      signd(i)          * (dest(i+1) - dest(i)) +
     @      signd(i+ly-jdiff) * (dest(i-jdiff) - dest(i)) +
     @      signd(i+ly)       * (dest(i+jdiff) - dest(i)) +
     @      signd(i+lz-kdiff) * (dest(i-kdiff) - dest(i)) +
     @      signd(i+lz)       * (dest(i+kdiff) - dest(i))))
      end

      subroutine hgscon(signd,
     @ dims(sn),
     @ sigx, sigy, sigz,
     @ dims(sc),
     @ dims(reg),
     @ hx, hy, hz)
      integer dims(sn)
      integer dims(sc)
      integer dims(reg)
      real*8 signd(dimdec(sn), 3)
      real*8 sigx(dimdec(sc))
      real*8 sigy(dimdec(sc))
      real*8 sigz(dimdec(sc))
      real*8 hx, hy, hz
      real*8 facx, facy, facz
      integer i, j, k
      facx = 0.25d0 / (hx*hx)
      facy = 0.25d0 / (hy*hy)
      facz = 0.25d0 / (hz*hz)
      do 10 k = regl2, regh2
         do 10 j = regl1, regh1
            do 10 i = regl0-1, regh0
 10            signd(i,j,k,1) = facx * (sigx(i,j-1,k-1) + sigx(i,j-1,k) +
     @                                  sigx(i,j,k-1)   + sigx(i,j,k))
      do 20 k = regl2, regh2
         do 20 j = regl1-1, regh1
            do 20 i = regl0, regh0
 20            signd(i,j,k,2) = facy * (sigy(i-1,j,k-1) + sigy(i-1,j,k) +
     @                                  sigy(i,j,k-1)   + sigy(i,j,k))
      do 30 k = regl2-1, regh2
         do 30 j = regl1, regh1
            do 30 i = regl0, regh0
 30            signd(i,j,k,3) = facz * (sigz(i-1,j-1,k) + sigz(i-1,j,k) +
     @                                  sigz(i,j-1,k)   + sigz(i,j,k))
      end

#else

      subroutine hgcen(cen,
     @ dims(cen),
     @ sigx, sigy, sigz,
     @ dims(sb),
     @ dims(reg),
     @ hx, hy, hz)
      implicit real*8(a-h,o-z), integer(i-n)
      integer dims(cen)
      integer dims(sb)
      integer dims(reg)
      real*8 cen(dimdec(cen))
      real*8 sigx(dimdec(sb))
      real*8 sigy(dimdec(sb))
      real*8 sigz(dimdec(sb))
      real*8 hx, hy, hz
      hxm2 = 1.d0 / (hx*hx)
      hym2 = 1.d0 / (hy*hy)
      hzm2 = 1.d0 / (hz*hz)
      do 10 k = regl2, regh2
         do 10 j = regl1, regh1
            do 10 i = regl0, regh0
 10            cen(i,j,k) = 4.d0 /
     @                     (hxm2 * (sigx(i-1,j-1,k-1) + sigx(i-1,j-1,k) +
     @                              sigx(i-1,j,k-1)   + sigx(i-1,j,k) +
     @                              sigx(i,j-1,k-1)   + sigx(i,j-1,k) +
     @                              sigx(i,j,k-1)     + sigx(i,j,k)) +
     @                      hym2 * (sigy(i-1,j-1,k-1) + sigy(i-1,j-1,k) +
     @                              sigy(i-1,j,k-1)   + sigy(i-1,j,k) +
     @                              sigy(i,j-1,k-1)   + sigy(i,j-1,k) +
     @                              sigy(i,j,k-1)     + sigy(i,j,k)) +
     @                      hzm2 * (sigz(i-1,j-1,k-1) + sigz(i-1,j-1,k) +
     @                              sigz(i-1,j,k-1)   + sigz(i-1,j,k) +
     @                              sigz(i,j-1,k-1)   + sigz(i,j-1,k) +
     @                              sigz(i,j,k-1)     + sigz(i,j,k)))
      end

#define AVG 0.25d0 * (hxm2 * \
               ((sigx(i-1,j-1,k-1) + sigx(i-1,j-1,k) + \
                 sigx(i-1,j,k-1)   + sigx(i-1,j,k)) * cor(i-1,j,k) + \
                (sigx(i,j-1,k-1)   + sigx(i,j-1,k) + \
                 sigx(i,j,k-1)     + sigx(i,j,k)) * cor(i+1,j,k)) + \
                      hym2 * \
               ((sigy(i-1,j-1,k-1) + sigy(i-1,j-1,k) + \
                 sigy(i,j-1,k-1)   + sigy(i,j-1,k)) * cor(i,j-1,k) + \
                (sigy(i-1,j,k-1)   + sigy(i-1,j,k) + \
                 sigy(i,j,k-1)     + sigy(i,j,k)) * cor(i,j+1,k)) + \
                      hzm2 * \
               ((sigz(i-1,j-1,k-1) + sigz(i-1,j,k-1) + \
                 sigz(i,j-1,k-1)   + sigz(i,j,k-1)) * cor(i,j,k-1) + \
                (sigz(i-1,j-1,k)   + sigz(i-1,j,k) + \
                 sigz(i,j-1,k)     + sigz(i,j,k)) * cor(i,j,k+1)))

      subroutine hgrlx(cor,
     @ dims(cor),
     @ res,
     @ dims(res),
     @ sigx, sigy, sigz,
     @ dims(sf),
     @ cen,
     @ dims(cen),
     @ dims(reg),
     @ hx, hy, hz)
      implicit real*8(a-h,o-z), integer(i-n)
      integer dims(cor)
      integer dims(res)
      integer dims(sf)
      integer dims(cen)
      integer dims(reg)
      real*8 cor(dimdec(cor))
      real*8 res(dimdec(res))
      real*8 sigx(dimdec(sf))
      real*8 sigy(dimdec(sf))
      real*8 sigz(dimdec(sf))
      real*8 cen(dimdec(cen))
      real*8 hx, hy, hz
#if 0
      avg() = 0.25d0 * (hxm2 *
     @          ((sigx(i-1,j-1,k-1) + sigx(i-1,j-1,k) +
     @            sigx(i-1,j,k-1)   + sigx(i-1,j,k)) * cor(i-1,j,k) +
     @           (sigx(i,j-1,k-1)   + sigx(i,j-1,k) +
     @            sigx(i,j,k-1)     + sigx(i,j,k)) * cor(i+1,j,k)) +
     @                 hym2 *
     @          ((sigy(i-1,j-1,k-1) + sigy(i-1,j-1,k) +
     @            sigy(i,j-1,k-1)   + sigy(i,j-1,k)) * cor(i,j-1,k) +
     @           (sigy(i-1,j,k-1)   + sigy(i-1,j,k) +
     @            sigy(i,j,k-1)     + sigy(i,j,k)) * cor(i,j+1,k)) +
     @                 hzm2 *
     @          ((sigz(i-1,j-1,k-1) + sigz(i-1,j,k-1) +
     @            sigz(i,j-1,k-1)   + sigz(i,j,k-1)) * cor(i,j,k-1) +
     @           (sigz(i-1,j-1,k)   + sigz(i-1,j,k) +
     @            sigz(i,j-1,k)     + sigz(i,j,k)) * cor(i,j,k+1)))
#endif
      hxm2 = 1.d0 / (hx*hx)
      hym2 = 1.d0 / (hy*hy)
      hzm2 = 1.d0 / (hz*hz)
      if (regh2 - regl2 .gt. regh1 - regl1 .and.
     @    regh2 - regl2 .gt. regh0 - regl0) then
         do 10 ipass = 0, 1
            ipar0 = ipass
            do 10 j = regl1, regh1
               ipar0 = 1 - ipar0
               ipar = ipar0
               do 10 i = regl0, regh0
                  ipar = 1 - ipar
cdir$ ivdep
                  do 10 k = regl2 + ipar, regh2, 2
                     cor(i,j,k) = (AVG - res(i,j,k)) * cen(i,j,k)
 10                  continue
      else if (regh1 - regl1 .gt. regh0 - regl0) then
         do 20 ipass = 0, 1
            ipar0 = ipass
            do 20 k = regl2, regh2
               ipar0 = 1 - ipar0
               ipar = ipar0
               do 20 i = regl0, regh0
                  ipar = 1 - ipar
cdir$ ivdep
                  do 20 j = regl1 + ipar, regh1, 2
                     cor(i,j,k) = (AVG - res(i,j,k)) * cen(i,j,k)
 20                  continue
      else
         do 30 ipass = 0, 1
            ipar0 = ipass
            do 30 k = regl2, regh2
               ipar0 = 1 - ipar0
               ipar = ipar0
               do 30 j = regl1, regh1
                  ipar = 1 - ipar
cdir$ ivdep
                  do 30 i = regl0 + ipar, regh0, 2
                     cor(i,j,k) = (AVG - res(i,j,k)) * cen(i,j,k)
 30                  continue
      endif
      end

      subroutine hgres(res,
     @ dims(res),
     @ src,
     @ dims(src),
     @ dest,
     @ dims(dest),
     @ sigx, sigy, sigz,
     @ dims(sf),
     @ dims(reg),
     @ hx, hy, hz)
      implicit real*8(a-h,o-z), integer(i-n)
      integer dims(res)
      integer dims(src)
      integer dims(dest)
      integer dims(sf)
      integer dims(reg)
      real*8 res(dimdec(res))
      real*8 src(dimdec(src))
      real*8 dest(dimdec(dest))
      real*8 sigx(dimdec(sf))
      real*8 sigy(dimdec(sf))
      real*8 sigz(dimdec(sf))
      real*8 hx, hy, hz
      hxm2 = 1.d0 / (hx*hx)
      hym2 = 1.d0 / (hy*hy)
      hzm2 = 1.d0 / (hz*hz)
      if (regh2 - regl2 .gt. regh1 - regl1 .and.
     @    regh2 - regl2 .gt. regh0 - regl0) then
         do 10 j = regl1, regh1
            do 10 i = regl0, regh0
               do 10 k = regl2, regh2
 10               res(i,j,k) = src(i,j,k) - 0.25d0 *        (hxm2 *
     @              ((sigx(i-1,j-1,k-1) + sigx(i-1,j-1,k) +
     @                sigx(i-1,j,k-1)   + sigx(i-1,j,k)) *
     @                (dest(i-1,j,k) - dest(i,j,k)) +
     @               (sigx(i,j-1,k-1)   + sigx(i,j-1,k) +
     @                sigx(i,j,k-1)     + sigx(i,j,k)) *
     @                (dest(i+1,j,k) - dest(i,j,k))) +       hym2 *
     @              ((sigy(i-1,j-1,k-1) + sigy(i-1,j-1,k) +
     @                sigy(i,j-1,k-1)   + sigy(i,j-1,k)) *
     @                (dest(i,j-1,k) - dest(i,j,k)) +
     @               (sigy(i-1,j,k-1)   + sigy(i-1,j,k) +
     @                sigy(i,j,k-1)     + sigy(i,j,k)) *
     @                (dest(i,j+1,k) - dest(i,j,k))) +       hzm2 *
     @              ((sigz(i-1,j-1,k-1) + sigz(i-1,j,k-1) +
     @                sigz(i,j-1,k-1)   + sigz(i,j,k-1)) *
     @                (dest(i,j,k-1) - dest(i,j,k)) +
     @               (sigz(i-1,j-1,k)   + sigz(i-1,j,k) +
     @                sigz(i,j-1,k)     + sigz(i,j,k)) *
     @                (dest(i,j,k+1) - dest(i,j,k)))             )
      else if (regh1 - regl1 .gt. regh0 - regl0) then
         do 20 k = regl2, regh2
            do 20 i = regl0, regh0
               do 20 j = regl1, regh1
 20               res(i,j,k) = src(i,j,k) - 0.25d0 *        (hxm2 *
     @              ((sigx(i-1,j-1,k-1) + sigx(i-1,j-1,k) +
     @                sigx(i-1,j,k-1)   + sigx(i-1,j,k)) *
     @                (dest(i-1,j,k) - dest(i,j,k)) +
     @               (sigx(i,j-1,k-1)   + sigx(i,j-1,k) +
     @                sigx(i,j,k-1)     + sigx(i,j,k)) *
     @                (dest(i+1,j,k) - dest(i,j,k))) +       hym2 *
     @              ((sigy(i-1,j-1,k-1) + sigy(i-1,j-1,k) +
     @                sigy(i,j-1,k-1)   + sigy(i,j-1,k)) *
     @                (dest(i,j-1,k) - dest(i,j,k)) +
     @               (sigy(i-1,j,k-1)   + sigy(i-1,j,k) +
     @                sigy(i,j,k-1)     + sigy(i,j,k)) *
     @                (dest(i,j+1,k) - dest(i,j,k))) +       hzm2 *
     @              ((sigz(i-1,j-1,k-1) + sigz(i-1,j,k-1) +
     @                sigz(i,j-1,k-1)   + sigz(i,j,k-1)) *
     @                (dest(i,j,k-1) - dest(i,j,k)) +
     @               (sigz(i-1,j-1,k)   + sigz(i-1,j,k) +
     @                sigz(i,j-1,k)     + sigz(i,j,k)) *
     @                (dest(i,j,k+1) - dest(i,j,k)))             )
      else
         do 30 k = regl2, regh2
            do 30 j = regl1, regh1
               do 30 i = regl0, regh0
 30               res(i,j,k) = src(i,j,k) - 0.25d0 *        (hxm2 *
     @              ((sigx(i-1,j-1,k-1) + sigx(i-1,j-1,k) +
     @                sigx(i-1,j,k-1)   + sigx(i-1,j,k)) *
     @                (dest(i-1,j,k) - dest(i,j,k)) +
     @               (sigx(i,j-1,k-1)   + sigx(i,j-1,k) +
     @                sigx(i,j,k-1)     + sigx(i,j,k)) *
     @                (dest(i+1,j,k) - dest(i,j,k))) +       hym2 *
     @              ((sigy(i-1,j-1,k-1) + sigy(i-1,j-1,k) +
     @                sigy(i,j-1,k-1)   + sigy(i,j-1,k)) *
     @                (dest(i,j-1,k) - dest(i,j,k)) +
     @               (sigy(i-1,j,k-1)   + sigy(i-1,j,k) +
     @                sigy(i,j,k-1)     + sigy(i,j,k)) *
     @                (dest(i,j+1,k) - dest(i,j,k))) +       hzm2 *
     @              ((sigz(i-1,j-1,k-1) + sigz(i-1,j,k-1) +
     @                sigz(i,j-1,k-1)   + sigz(i,j,k-1)) *
     @                (dest(i,j,k-1) - dest(i,j,k)) +
     @               (sigz(i-1,j-1,k)   + sigz(i-1,j,k) +
     @                sigz(i,j-1,k)     + sigz(i,j,k)) *
     @                (dest(i,j,k+1) - dest(i,j,k)))             )
      endif
      end

#endif

#endif

#elif (defined HG_TERRAIN)

c twenty-seven-point terrain stencils

      subroutine hgcen(cen, dims(cen),
     @ sig, dims(sb),
     @ dims(reg))
      integer dims(cen)
      integer dims(sb)
      integer dims(reg)
      real*8 cen(dimdec(cen))
      real*8 sig(dimdec(sb), 5)
      integer i, j, k
      do 10 k = regl2, regh2
         do 10 j = regl1, regh1
            do 10 i = regl0, regh0
 10            cen(i,j,k) = 36.d0 /
     @           (4.d0 * (sig(i-1,j-1,k-1,1) + sig(i-1,j-1,k,1) +
     @                    sig(i-1,j,k-1,1)   + sig(i-1,j,k,1) +
     @                    sig(i,j-1,k-1,1)   + sig(i,j-1,k,1) +
     @                    sig(i,j,k-1,1)     + sig(i,j,k,1) +
     @                    sig(i-1,j-1,k-1,2) + sig(i-1,j-1,k,2) +
     @                    sig(i-1,j,k-1,2)   + sig(i-1,j,k,2) +
     @                    sig(i,j-1,k-1,2)   + sig(i,j-1,k,2) +
     @                    sig(i,j,k-1,2)     + sig(i,j,k,2) +
     @                    sig(i-1,j-1,k-1,3) + sig(i-1,j-1,k,3) +
     @                    sig(i-1,j,k-1,3)   + sig(i-1,j,k,3) +
     @                    sig(i,j-1,k-1,3)   + sig(i,j-1,k,3) +
     @                    sig(i,j,k-1,3)     + sig(i,j,k,3)) +
     @            6.d0 * (sig(i-1,j-1,k,4)   - sig(i-1,j-1,k-1,4) +
     @                    sig(i-1,j,k,4)     - sig(i-1,j,k-1,4) +
     @                    sig(i,j-1,k-1,4)   - sig(i,j-1,k,4) +
     @                    sig(i,j,k-1,4)     - sig(i,j,k,4) +
     @                    sig(i-1,j-1,k,5)   - sig(i-1,j-1,k-1,5) +
     @                    sig(i-1,j,k-1,5)   - sig(i-1,j,k,5) +
     @                    sig(i,j-1,k,5)     - sig(i,j-1,k-1,5) +
     @                    sig(i,j,k-1,5)     - sig(i,j,k,5)))
      end

#define MXX(ii,jj,kk) (4.d0 *  cor(ii,j,k)   + \
                       2.d0 * (cor(ii,jj,k)  - cor(i,jj,k) + \
                               cor(ii,j,kk)  - cor(i,j,kk)) + \
                              (cor(ii,jj,kk) - cor(i,jj,kk)))

#define MYY(ii,jj,kk) (4.d0 *  cor(i,jj,k)   + \
                       2.d0 * (cor(ii,jj,k)  - cor(ii,j,k) + \
                               cor(i,jj,kk)  - cor(i,j,kk)) + \
                              (cor(ii,jj,kk) - cor(ii,j,kk)))

#define MZZ(ii,jj,kk) (4.d0 *  cor(i,j,kk)   + \
                       2.d0 * (cor(ii,j,kk)  - cor(ii,j,k) + \
                               cor(i,jj,kk)  - cor(i,jj,k)) + \
                              (cor(ii,jj,kk) - cor(ii,jj,k)))

#define MXZ(ii,jj,kk) (6.d0 *  cor(ii,j,kk)  + \
                       3.d0 * (cor(ii,jj,kk) - cor(i,jj,k)))

#define MYZ(ii,jj,kk) (6.d0 *  cor(i,jj,kk)  + \
                       3.d0 * (cor(ii,jj,kk) - cor(ii,j,k)))

      subroutine hgrlx(cor, dims(cor),
     @ res, dims(res),
     @ sig, dims(sf),
     @ cen, dims(cen),
     @ dims(reg))
      integer dims(cor)
      integer dims(res)
      integer dims(sf)
      integer dims(cen)
      integer dims(reg)
      real*8 cor(dimdec(cor))
      real*8 res(dimdec(res))
      real*8 sig(dimdec(sf), 5)
      real*8 cen(dimdec(cen))
      real*8 fac, tmp
      integer i, j, k
      fac = 1.d0 / 36.d0
      if (regh2 - regl2 .gt. regh1 - regl1 .and.
     @    regh2 - regl2 .gt. regh0 - regl0) then
         do 10 j = regl1, regh1
            do 10 i = regl0, regh0
cdir$ ivdep
               do 10 k = regl2, regh2
                  tmp = sig(i-1,j-1,k-1,1) * MXX(i-1,j-1,k-1) +
     @                  sig(i-1,j-1,k  ,1) * MXX(i-1,j-1,k+1) +
     @                  sig(i-1,j  ,k-1,1) * MXX(i-1,j+1,k-1) +
     @                  sig(i-1,j  ,k  ,1) * MXX(i-1,j+1,k+1)
                  tmp = tmp +
     @                  sig(i  ,j-1,k-1,1) * MXX(i+1,j-1,k-1) +
     @                  sig(i  ,j-1,k  ,1) * MXX(i+1,j-1,k+1) +
     @                  sig(i  ,j  ,k-1,1) * MXX(i+1,j+1,k-1) +
     @                  sig(i  ,j  ,k  ,1) * MXX(i+1,j+1,k+1)
                  tmp = tmp +
     @                  sig(i-1,j-1,k-1,2) * MYY(i-1,j-1,k-1) +
     @                  sig(i-1,j-1,k  ,2) * MYY(i-1,j-1,k+1) +
     @                  sig(i-1,j  ,k-1,2) * MYY(i-1,j+1,k-1) +
     @                  sig(i-1,j  ,k  ,2) * MYY(i-1,j+1,k+1)
                  tmp = tmp +
     @                  sig(i  ,j-1,k-1,2) * MYY(i+1,j-1,k-1) +
     @                  sig(i  ,j-1,k  ,2) * MYY(i+1,j-1,k+1) +
     @                  sig(i  ,j  ,k-1,2) * MYY(i+1,j+1,k-1) +
     @                  sig(i  ,j  ,k  ,2) * MYY(i+1,j+1,k+1)
                  tmp = tmp +
     @                  sig(i-1,j-1,k-1,3) * MZZ(i-1,j-1,k-1) +
     @                  sig(i-1,j-1,k  ,3) * MZZ(i-1,j-1,k+1) +
     @                  sig(i-1,j  ,k-1,3) * MZZ(i-1,j+1,k-1) +
     @                  sig(i-1,j  ,k  ,3) * MZZ(i-1,j+1,k+1)
                  tmp = tmp +
     @                  sig(i  ,j-1,k-1,3) * MZZ(i+1,j-1,k-1) +
     @                  sig(i  ,j-1,k  ,3) * MZZ(i+1,j-1,k+1) +
     @                  sig(i  ,j  ,k-1,3) * MZZ(i+1,j+1,k-1) +
     @                  sig(i  ,j  ,k  ,3) * MZZ(i+1,j+1,k+1)
                  tmp = tmp -
     @                  sig(i-1,j-1,k-1,4) * MXZ(i-1,j-1,k-1) +
     @                  sig(i-1,j-1,k  ,4) * MXZ(i-1,j-1,k+1) -
     @                  sig(i-1,j  ,k-1,4) * MXZ(i-1,j+1,k-1) +
     @                  sig(i-1,j  ,k  ,4) * MXZ(i-1,j+1,k+1) +
     @                  sig(i  ,j-1,k-1,4) * MXZ(i+1,j-1,k-1) -
     @                  sig(i  ,j-1,k  ,4) * MXZ(i+1,j-1,k+1) +
     @                  sig(i  ,j  ,k-1,4) * MXZ(i+1,j+1,k-1) -
     @                  sig(i  ,j  ,k  ,4) * MXZ(i+1,j+1,k+1)
                  tmp = tmp -
     @                  sig(i-1,j-1,k-1,5) * MYZ(i-1,j-1,k-1) +
     @                  sig(i-1,j-1,k  ,5) * MYZ(i-1,j-1,k+1) +
     @                  sig(i-1,j  ,k-1,5) * MYZ(i-1,j+1,k-1) -
     @                  sig(i-1,j  ,k  ,5) * MYZ(i-1,j+1,k+1) -
     @                  sig(i  ,j-1,k-1,5) * MYZ(i+1,j-1,k-1) +
     @                  sig(i  ,j-1,k  ,5) * MYZ(i+1,j-1,k+1) +
     @                  sig(i  ,j  ,k-1,5) * MYZ(i+1,j+1,k-1) -
     @                  sig(i  ,j  ,k  ,5) * MYZ(i+1,j+1,k+1)
                  cor(i,j,k) = (fac * tmp - res(i,j,k)) * cen(i,j,k)
 10               continue
      else if (regh1 - regl1 .gt. regh0 - regl0) then
         do 20 k = regl2, regh2
            do 20 i = regl0, regh0
cdir$ ivdep
               do 20 j = regl1, regh1
                  tmp = sig(i-1,j-1,k-1,1) * MXX(i-1,j-1,k-1) +
     @                  sig(i-1,j-1,k  ,1) * MXX(i-1,j-1,k+1) +
     @                  sig(i-1,j  ,k-1,1) * MXX(i-1,j+1,k-1) +
     @                  sig(i-1,j  ,k  ,1) * MXX(i-1,j+1,k+1)
                  tmp = tmp +
     @                  sig(i  ,j-1,k-1,1) * MXX(i+1,j-1,k-1) +
     @                  sig(i  ,j-1,k  ,1) * MXX(i+1,j-1,k+1) +
     @                  sig(i  ,j  ,k-1,1) * MXX(i+1,j+1,k-1) +
     @                  sig(i  ,j  ,k  ,1) * MXX(i+1,j+1,k+1)
                  tmp = tmp +
     @                  sig(i-1,j-1,k-1,2) * MYY(i-1,j-1,k-1) +
     @                  sig(i-1,j-1,k  ,2) * MYY(i-1,j-1,k+1) +
     @                  sig(i-1,j  ,k-1,2) * MYY(i-1,j+1,k-1) +
     @                  sig(i-1,j  ,k  ,2) * MYY(i-1,j+1,k+1)
                  tmp = tmp +
     @                  sig(i  ,j-1,k-1,2) * MYY(i+1,j-1,k-1) +
     @                  sig(i  ,j-1,k  ,2) * MYY(i+1,j-1,k+1) +
     @                  sig(i  ,j  ,k-1,2) * MYY(i+1,j+1,k-1) +
     @                  sig(i  ,j  ,k  ,2) * MYY(i+1,j+1,k+1)
                  tmp = tmp +
     @                  sig(i-1,j-1,k-1,3) * MZZ(i-1,j-1,k-1) +
     @                  sig(i-1,j-1,k  ,3) * MZZ(i-1,j-1,k+1) +
     @                  sig(i-1,j  ,k-1,3) * MZZ(i-1,j+1,k-1) +
     @                  sig(i-1,j  ,k  ,3) * MZZ(i-1,j+1,k+1)
                  tmp = tmp +
     @                  sig(i  ,j-1,k-1,3) * MZZ(i+1,j-1,k-1) +
     @                  sig(i  ,j-1,k  ,3) * MZZ(i+1,j-1,k+1) +
     @                  sig(i  ,j  ,k-1,3) * MZZ(i+1,j+1,k-1) +
     @                  sig(i  ,j  ,k  ,3) * MZZ(i+1,j+1,k+1)
                  tmp = tmp -
     @                  sig(i-1,j-1,k-1,4) * MXZ(i-1,j-1,k-1) +
     @                  sig(i-1,j-1,k  ,4) * MXZ(i-1,j-1,k+1) -
     @                  sig(i-1,j  ,k-1,4) * MXZ(i-1,j+1,k-1) +
     @                  sig(i-1,j  ,k  ,4) * MXZ(i-1,j+1,k+1) +
     @                  sig(i  ,j-1,k-1,4) * MXZ(i+1,j-1,k-1) -
     @                  sig(i  ,j-1,k  ,4) * MXZ(i+1,j-1,k+1) +
     @                  sig(i  ,j  ,k-1,4) * MXZ(i+1,j+1,k-1) -
     @                  sig(i  ,j  ,k  ,4) * MXZ(i+1,j+1,k+1)
                  tmp = tmp -
     @                  sig(i-1,j-1,k-1,5) * MYZ(i-1,j-1,k-1) +
     @                  sig(i-1,j-1,k  ,5) * MYZ(i-1,j-1,k+1) +
     @                  sig(i-1,j  ,k-1,5) * MYZ(i-1,j+1,k-1) -
     @                  sig(i-1,j  ,k  ,5) * MYZ(i-1,j+1,k+1) -
     @                  sig(i  ,j-1,k-1,5) * MYZ(i+1,j-1,k-1) +
     @                  sig(i  ,j-1,k  ,5) * MYZ(i+1,j-1,k+1) +
     @                  sig(i  ,j  ,k-1,5) * MYZ(i+1,j+1,k-1) -
     @                  sig(i  ,j  ,k  ,5) * MYZ(i+1,j+1,k+1)
                  cor(i,j,k) = (fac * tmp - res(i,j,k)) * cen(i,j,k)
 20               continue
      else
         do 30 k = regl2, regh2
            do 30 j = regl1, regh1
cdir$ ivdep
               do 30 i = regl0, regh0
                  tmp = sig(i-1,j-1,k-1,1) * MXX(i-1,j-1,k-1) +
     @                  sig(i-1,j-1,k  ,1) * MXX(i-1,j-1,k+1) +
     @                  sig(i-1,j  ,k-1,1) * MXX(i-1,j+1,k-1) +
     @                  sig(i-1,j  ,k  ,1) * MXX(i-1,j+1,k+1)
                  tmp = tmp +
     @                  sig(i  ,j-1,k-1,1) * MXX(i+1,j-1,k-1) +
     @                  sig(i  ,j-1,k  ,1) * MXX(i+1,j-1,k+1) +
     @                  sig(i  ,j  ,k-1,1) * MXX(i+1,j+1,k-1) +
     @                  sig(i  ,j  ,k  ,1) * MXX(i+1,j+1,k+1)
                  tmp = tmp +
     @                  sig(i-1,j-1,k-1,2) * MYY(i-1,j-1,k-1) +
     @                  sig(i-1,j-1,k  ,2) * MYY(i-1,j-1,k+1) +
     @                  sig(i-1,j  ,k-1,2) * MYY(i-1,j+1,k-1) +
     @                  sig(i-1,j  ,k  ,2) * MYY(i-1,j+1,k+1)
                  tmp = tmp +
     @                  sig(i  ,j-1,k-1,2) * MYY(i+1,j-1,k-1) +
     @                  sig(i  ,j-1,k  ,2) * MYY(i+1,j-1,k+1) +
     @                  sig(i  ,j  ,k-1,2) * MYY(i+1,j+1,k-1) +
     @                  sig(i  ,j  ,k  ,2) * MYY(i+1,j+1,k+1)
                  tmp = tmp +
     @                  sig(i-1,j-1,k-1,3) * MZZ(i-1,j-1,k-1) +
     @                  sig(i-1,j-1,k  ,3) * MZZ(i-1,j-1,k+1) +
     @                  sig(i-1,j  ,k-1,3) * MZZ(i-1,j+1,k-1) +
     @                  sig(i-1,j  ,k  ,3) * MZZ(i-1,j+1,k+1)
                  tmp = tmp +
     @                  sig(i  ,j-1,k-1,3) * MZZ(i+1,j-1,k-1) +
     @                  sig(i  ,j-1,k  ,3) * MZZ(i+1,j-1,k+1) +
     @                  sig(i  ,j  ,k-1,3) * MZZ(i+1,j+1,k-1) +
     @                  sig(i  ,j  ,k  ,3) * MZZ(i+1,j+1,k+1)
                  tmp = tmp -
     @                  sig(i-1,j-1,k-1,4) * MXZ(i-1,j-1,k-1) +
     @                  sig(i-1,j-1,k  ,4) * MXZ(i-1,j-1,k+1) -
     @                  sig(i-1,j  ,k-1,4) * MXZ(i-1,j+1,k-1) +
     @                  sig(i-1,j  ,k  ,4) * MXZ(i-1,j+1,k+1) +
     @                  sig(i  ,j-1,k-1,4) * MXZ(i+1,j-1,k-1) -
     @                  sig(i  ,j-1,k  ,4) * MXZ(i+1,j-1,k+1) +
     @                  sig(i  ,j  ,k-1,4) * MXZ(i+1,j+1,k-1) -
     @                  sig(i  ,j  ,k  ,4) * MXZ(i+1,j+1,k+1)
                  tmp = tmp -
     @                  sig(i-1,j-1,k-1,5) * MYZ(i-1,j-1,k-1) +
     @                  sig(i-1,j-1,k  ,5) * MYZ(i-1,j-1,k+1) +
     @                  sig(i-1,j  ,k-1,5) * MYZ(i-1,j+1,k-1) -
     @                  sig(i-1,j  ,k  ,5) * MYZ(i-1,j+1,k+1) -
     @                  sig(i  ,j-1,k-1,5) * MYZ(i+1,j-1,k-1) +
     @                  sig(i  ,j-1,k  ,5) * MYZ(i+1,j-1,k+1) +
     @                  sig(i  ,j  ,k-1,5) * MYZ(i+1,j+1,k-1) -
     @                  sig(i  ,j  ,k  ,5) * MYZ(i+1,j+1,k+1)
                  cor(i,j,k) = (fac * tmp - res(i,j,k)) * cen(i,j,k)
 30               continue
      endif
      end

c cor substituted for dest to use same macro:

      subroutine hgres(res, dims(res),
     @ src,  dims(src),
     @ cor, dims(dest),
     @ sig, dims(sf),
     @ cen, dims(cen),
     @ dims(reg))
      integer dims(res)
      integer dims(src)
      integer dims(dest)
      integer dims(sf)
      integer dims(cen)
      integer dims(reg)
      real*8 res(dimdec(res))
      real*8 src(dimdec(src))
      real*8 cor(dimdec(dest))
      real*8 sig(dimdec(sf), 5)
      real*8 cen(dimdec(cen))
      real*8 fac, tmp
      integer i, j, k
      fac = 1.d0 / 36.d0
      if (regh2 - regl2 .gt. regh1 - regl1 .and.
     @    regh2 - regl2 .gt. regh0 - regl0) then
         do 10 j = regl1, regh1
            do 10 i = regl0, regh0
               do 10 k = regl2, regh2
                  tmp = sig(i-1,j-1,k-1,1) * MXX(i-1,j-1,k-1) +
     @                  sig(i-1,j-1,k  ,1) * MXX(i-1,j-1,k+1) +
     @                  sig(i-1,j  ,k-1,1) * MXX(i-1,j+1,k-1) +
     @                  sig(i-1,j  ,k  ,1) * MXX(i-1,j+1,k+1)
                  tmp = tmp +
     @                  sig(i  ,j-1,k-1,1) * MXX(i+1,j-1,k-1) +
     @                  sig(i  ,j-1,k  ,1) * MXX(i+1,j-1,k+1) +
     @                  sig(i  ,j  ,k-1,1) * MXX(i+1,j+1,k-1) +
     @                  sig(i  ,j  ,k  ,1) * MXX(i+1,j+1,k+1)
                  tmp = tmp +
     @                  sig(i-1,j-1,k-1,2) * MYY(i-1,j-1,k-1) +
     @                  sig(i-1,j-1,k  ,2) * MYY(i-1,j-1,k+1) +
     @                  sig(i-1,j  ,k-1,2) * MYY(i-1,j+1,k-1) +
     @                  sig(i-1,j  ,k  ,2) * MYY(i-1,j+1,k+1)
                  tmp = tmp +
     @                  sig(i  ,j-1,k-1,2) * MYY(i+1,j-1,k-1) +
     @                  sig(i  ,j-1,k  ,2) * MYY(i+1,j-1,k+1) +
     @                  sig(i  ,j  ,k-1,2) * MYY(i+1,j+1,k-1) +
     @                  sig(i  ,j  ,k  ,2) * MYY(i+1,j+1,k+1)
                  tmp = tmp +
     @                  sig(i-1,j-1,k-1,3) * MZZ(i-1,j-1,k-1) +
     @                  sig(i-1,j-1,k  ,3) * MZZ(i-1,j-1,k+1) +
     @                  sig(i-1,j  ,k-1,3) * MZZ(i-1,j+1,k-1) +
     @                  sig(i-1,j  ,k  ,3) * MZZ(i-1,j+1,k+1)
                  tmp = tmp +
     @                  sig(i  ,j-1,k-1,3) * MZZ(i+1,j-1,k-1) +
     @                  sig(i  ,j-1,k  ,3) * MZZ(i+1,j-1,k+1) +
     @                  sig(i  ,j  ,k-1,3) * MZZ(i+1,j+1,k-1) +
     @                  sig(i  ,j  ,k  ,3) * MZZ(i+1,j+1,k+1)
                  tmp = tmp -
     @                  sig(i-1,j-1,k-1,4) * MXZ(i-1,j-1,k-1) +
     @                  sig(i-1,j-1,k  ,4) * MXZ(i-1,j-1,k+1) -
     @                  sig(i-1,j  ,k-1,4) * MXZ(i-1,j+1,k-1) +
     @                  sig(i-1,j  ,k  ,4) * MXZ(i-1,j+1,k+1) +
     @                  sig(i  ,j-1,k-1,4) * MXZ(i+1,j-1,k-1) -
     @                  sig(i  ,j-1,k  ,4) * MXZ(i+1,j-1,k+1) +
     @                  sig(i  ,j  ,k-1,4) * MXZ(i+1,j+1,k-1) -
     @                  sig(i  ,j  ,k  ,4) * MXZ(i+1,j+1,k+1)
                  tmp = tmp -
     @                  sig(i-1,j-1,k-1,5) * MYZ(i-1,j-1,k-1) +
     @                  sig(i-1,j-1,k  ,5) * MYZ(i-1,j-1,k+1) +
     @                  sig(i-1,j  ,k-1,5) * MYZ(i-1,j+1,k-1) -
     @                  sig(i-1,j  ,k  ,5) * MYZ(i-1,j+1,k+1) -
     @                  sig(i  ,j-1,k-1,5) * MYZ(i+1,j-1,k-1) +
     @                  sig(i  ,j-1,k  ,5) * MYZ(i+1,j-1,k+1) +
     @                  sig(i  ,j  ,k-1,5) * MYZ(i+1,j+1,k-1) -
     @                  sig(i  ,j  ,k  ,5) * MYZ(i+1,j+1,k+1)
                  if (cen(i,j,k) .ne. 0.d0) then
                     res(i,j,k) = src(i,j,k)
     @                  - fac * tmp + cor(i,j,k) / cen(i,j,k)
                  else
                     res(i,j,k) = 0.d0
                  endif
 10               continue
      else if (regh1 - regl1 .gt. regh0 - regl0) then
         do 20 k = regl2, regh2
            do 20 i = regl0, regh0
               do 20 j = regl1, regh1
                  tmp = sig(i-1,j-1,k-1,1) * MXX(i-1,j-1,k-1) +
     @                  sig(i-1,j-1,k  ,1) * MXX(i-1,j-1,k+1) +
     @                  sig(i-1,j  ,k-1,1) * MXX(i-1,j+1,k-1) +
     @                  sig(i-1,j  ,k  ,1) * MXX(i-1,j+1,k+1)
                  tmp = tmp +
     @                  sig(i  ,j-1,k-1,1) * MXX(i+1,j-1,k-1) +
     @                  sig(i  ,j-1,k  ,1) * MXX(i+1,j-1,k+1) +
     @                  sig(i  ,j  ,k-1,1) * MXX(i+1,j+1,k-1) +
     @                  sig(i  ,j  ,k  ,1) * MXX(i+1,j+1,k+1)
                  tmp = tmp +
     @                  sig(i-1,j-1,k-1,2) * MYY(i-1,j-1,k-1) +
     @                  sig(i-1,j-1,k  ,2) * MYY(i-1,j-1,k+1) +
     @                  sig(i-1,j  ,k-1,2) * MYY(i-1,j+1,k-1) +
     @                  sig(i-1,j  ,k  ,2) * MYY(i-1,j+1,k+1)
                  tmp = tmp +
     @                  sig(i  ,j-1,k-1,2) * MYY(i+1,j-1,k-1) +
     @                  sig(i  ,j-1,k  ,2) * MYY(i+1,j-1,k+1) +
     @                  sig(i  ,j  ,k-1,2) * MYY(i+1,j+1,k-1) +
     @                  sig(i  ,j  ,k  ,2) * MYY(i+1,j+1,k+1)
                  tmp = tmp +
     @                  sig(i-1,j-1,k-1,3) * MZZ(i-1,j-1,k-1) +
     @                  sig(i-1,j-1,k  ,3) * MZZ(i-1,j-1,k+1) +
     @                  sig(i-1,j  ,k-1,3) * MZZ(i-1,j+1,k-1) +
     @                  sig(i-1,j  ,k  ,3) * MZZ(i-1,j+1,k+1)
                  tmp = tmp +
     @                  sig(i  ,j-1,k-1,3) * MZZ(i+1,j-1,k-1) +
     @                  sig(i  ,j-1,k  ,3) * MZZ(i+1,j-1,k+1) +
     @                  sig(i  ,j  ,k-1,3) * MZZ(i+1,j+1,k-1) +
     @                  sig(i  ,j  ,k  ,3) * MZZ(i+1,j+1,k+1)
                  tmp = tmp -
     @                  sig(i-1,j-1,k-1,4) * MXZ(i-1,j-1,k-1) +
     @                  sig(i-1,j-1,k  ,4) * MXZ(i-1,j-1,k+1) -
     @                  sig(i-1,j  ,k-1,4) * MXZ(i-1,j+1,k-1) +
     @                  sig(i-1,j  ,k  ,4) * MXZ(i-1,j+1,k+1) +
     @                  sig(i  ,j-1,k-1,4) * MXZ(i+1,j-1,k-1) -
     @                  sig(i  ,j-1,k  ,4) * MXZ(i+1,j-1,k+1) +
     @                  sig(i  ,j  ,k-1,4) * MXZ(i+1,j+1,k-1) -
     @                  sig(i  ,j  ,k  ,4) * MXZ(i+1,j+1,k+1)
                  tmp = tmp -
     @                  sig(i-1,j-1,k-1,5) * MYZ(i-1,j-1,k-1) +
     @                  sig(i-1,j-1,k  ,5) * MYZ(i-1,j-1,k+1) +
     @                  sig(i-1,j  ,k-1,5) * MYZ(i-1,j+1,k-1) -
     @                  sig(i-1,j  ,k  ,5) * MYZ(i-1,j+1,k+1) -
     @                  sig(i  ,j-1,k-1,5) * MYZ(i+1,j-1,k-1) +
     @                  sig(i  ,j-1,k  ,5) * MYZ(i+1,j-1,k+1) +
     @                  sig(i  ,j  ,k-1,5) * MYZ(i+1,j+1,k-1) -
     @                  sig(i  ,j  ,k  ,5) * MYZ(i+1,j+1,k+1)
                  if (cen(i,j,k) .ne. 0.d0) then
                     res(i,j,k) = src(i,j,k)
     @                  - fac * tmp + cor(i,j,k) / cen(i,j,k)
                  else
                     res(i,j,k) = 0.d0
                  endif
 20               continue
      else
         do 30 k = regl2, regh2
            do 30 j = regl1, regh1
               do 30 i = regl0, regh0
                  tmp = sig(i-1,j-1,k-1,1) * MXX(i-1,j-1,k-1) +
     @                  sig(i-1,j-1,k  ,1) * MXX(i-1,j-1,k+1) +
     @                  sig(i-1,j  ,k-1,1) * MXX(i-1,j+1,k-1) +
     @                  sig(i-1,j  ,k  ,1) * MXX(i-1,j+1,k+1)
                  tmp = tmp +
     @                  sig(i  ,j-1,k-1,1) * MXX(i+1,j-1,k-1) +
     @                  sig(i  ,j-1,k  ,1) * MXX(i+1,j-1,k+1) +
     @                  sig(i  ,j  ,k-1,1) * MXX(i+1,j+1,k-1) +
     @                  sig(i  ,j  ,k  ,1) * MXX(i+1,j+1,k+1)
                  tmp = tmp +
     @                  sig(i-1,j-1,k-1,2) * MYY(i-1,j-1,k-1) +
     @                  sig(i-1,j-1,k  ,2) * MYY(i-1,j-1,k+1) +
     @                  sig(i-1,j  ,k-1,2) * MYY(i-1,j+1,k-1) +
     @                  sig(i-1,j  ,k  ,2) * MYY(i-1,j+1,k+1)
                  tmp = tmp +
     @                  sig(i  ,j-1,k-1,2) * MYY(i+1,j-1,k-1) +
     @                  sig(i  ,j-1,k  ,2) * MYY(i+1,j-1,k+1) +
     @                  sig(i  ,j  ,k-1,2) * MYY(i+1,j+1,k-1) +
     @                  sig(i  ,j  ,k  ,2) * MYY(i+1,j+1,k+1)
                  tmp = tmp +
     @                  sig(i-1,j-1,k-1,3) * MZZ(i-1,j-1,k-1) +
     @                  sig(i-1,j-1,k  ,3) * MZZ(i-1,j-1,k+1) +
     @                  sig(i-1,j  ,k-1,3) * MZZ(i-1,j+1,k-1) +
     @                  sig(i-1,j  ,k  ,3) * MZZ(i-1,j+1,k+1)
                  tmp = tmp +
     @                  sig(i  ,j-1,k-1,3) * MZZ(i+1,j-1,k-1) +
     @                  sig(i  ,j-1,k  ,3) * MZZ(i+1,j-1,k+1) +
     @                  sig(i  ,j  ,k-1,3) * MZZ(i+1,j+1,k-1) +
     @                  sig(i  ,j  ,k  ,3) * MZZ(i+1,j+1,k+1)
                  tmp = tmp -
     @                  sig(i-1,j-1,k-1,4) * MXZ(i-1,j-1,k-1) +
     @                  sig(i-1,j-1,k  ,4) * MXZ(i-1,j-1,k+1) -
     @                  sig(i-1,j  ,k-1,4) * MXZ(i-1,j+1,k-1) +
     @                  sig(i-1,j  ,k  ,4) * MXZ(i-1,j+1,k+1) +
     @                  sig(i  ,j-1,k-1,4) * MXZ(i+1,j-1,k-1) -
     @                  sig(i  ,j-1,k  ,4) * MXZ(i+1,j-1,k+1) +
     @                  sig(i  ,j  ,k-1,4) * MXZ(i+1,j+1,k-1) -
     @                  sig(i  ,j  ,k  ,4) * MXZ(i+1,j+1,k+1)
                  tmp = tmp -
     @                  sig(i-1,j-1,k-1,5) * MYZ(i-1,j-1,k-1) +
     @                  sig(i-1,j-1,k  ,5) * MYZ(i-1,j-1,k+1) +
     @                  sig(i-1,j  ,k-1,5) * MYZ(i-1,j+1,k-1) -
     @                  sig(i-1,j  ,k  ,5) * MYZ(i-1,j+1,k+1) -
     @                  sig(i  ,j-1,k-1,5) * MYZ(i+1,j-1,k-1) +
     @                  sig(i  ,j-1,k  ,5) * MYZ(i+1,j-1,k+1) +
     @                  sig(i  ,j  ,k-1,5) * MYZ(i+1,j+1,k-1) -
     @                  sig(i  ,j  ,k  ,5) * MYZ(i+1,j+1,k+1)
                  if (cen(i,j,k) .ne. 0.d0) then
                     res(i,j,k) = src(i,j,k)
     @                  - fac * tmp + cor(i,j,k) / cen(i,j,k)
                  else
                     res(i,j,k) = 0.d0
                  endif
 30               continue
      endif
      end

      subroutine hgrlnf(cor,
     @ dims(cor),
     @ res, dims(res),
     @ wrk, dims(wrk),
     @ signd, dims(sn),
     @ cen, dims(cen),
     @ dims(reg),
     @ dims(dom),
     @ lsd, ipass)
      integer dims(cor)
      integer dims(res)
      integer dims(wrk)
      integer dims(sn)
      integer dims(cen)
      integer dims(reg)
      integer dims(dom)
      real*8 cor(dimdec(cor))
      real*8 res(dimdec(res))
      real*8 wrk(dimdec(wrk))
      real*8 signd(dimdec(sn), 3)
      real*8 cen(dimdec(cen))
      integer lsd, ipass
      real*8 betm, aj
      integer i, j, k, ioff, ipar
      end

#endif

      subroutine hgrlnb(cor, dims(cor),
     @ wrk, dims(wrk),
     @ dims(reg),
     @ lsd, ipass)
      integer dims(cor)
      integer dims(wrk)
      integer dims(reg)
      real*8 cor(dimdec(cor))
      real*8 wrk(dimdec(wrk))
      integer lsd, ipass
      integer ioff, ipar, i, j, k
      if (lsd .eq. 0) then
         if (mod(regl1 + regl2, 2) .eq. 0) then
            ioff = 1 - ipass
         else
            ioff = ipass
         endif
c back substitution loop:
         do 10 i = regh0 - 1, regl0, -1
            ipar = ioff
            do 10 k = regl2, regh2
               ipar = 1 - ipar
               do 10 j = regl1 + ipar, regh1, 2
                  cor(i,j,k) = cor(i,j,k) - wrk(i,j,k) * cor(i+1,j,k)
 10               continue
      else if (lsd .eq. 1) then
         if (mod(regl0 + regl2, 2) .eq. 0) then
            ioff = 1 - ipass
         else
            ioff = ipass
         endif
c back substitution loop:
         do 20 j = regh1 - 1, regl1, -1
            ipar = ioff
            do 20 k = regl2, regh2
               ipar = 1 - ipar
               do 20 i = regl0 + ipar, regh0, 2
                  cor(i,j,k) = cor(i,j,k) - wrk(i,j,k) * cor(i,j+1,k)
 20               continue
      else
         if (mod(regl0 + regl1, 2) .eq. 0) then
            ioff = 1 - ipass
         else
            ioff = ipass
         endif
c back substitution loop:
         do 30 k = regh2 - 1, regl2, -1
            ipar = ioff
            do 30 j = regl1, regh1
               ipar = 1 - ipar
               do 30 i = regl0 + ipar, regh0, 2
                  cor(i,j,k) = cor(i,j,k) - wrk(i,j,k) * cor(i,j,k+1)
 30               continue
      endif
      end

#ifdef HG_CONSTANT
#  define CGOPT 2
#else
#  define CGOPT 1
#endif

#if (CGOPT == 1)

c Unrolled indexing in these 3 routines uses the fact that each array
c has a border of width 1

c Works for NODE-based data.
      subroutine hgip(v0, v1, mask,
     @ dims(reg),
     @ sum)
      integer dims(reg)
      real*8 v0(dimdecu(reg))
      real*8 v1(dimdecu(reg))
      real*8 mask(dimdecu(reg))
      real*8 sum
      integer i, jdiff, kdiff
      jdiff = regh0 - regl0 + 1
      kdiff = (regh1 - regl1 + 1) * jdiff
      do 10 i = kdiff + jdiff + 2, kdiff * (regh2 - regl2) - jdiff - 1
 10      sum = sum + mask(i) * v0(i) * v1(i)
      end

      subroutine hgcg1(r, p, z, x, w, c, mask,
     @ dims(reg),
     @ alpha, rho)
      integer dims(reg)
      real*8 r(dimdecu(reg))
      real*8 p(dimdecu(reg))
      real*8 z(dimdecu(reg))
      real*8 x(dimdecu(reg))
      real*8 w(dimdecu(reg))
      real*8 c(dimdecu(reg))
      real*8 mask(dimdecu(reg))
      real*8 alpha, rho
      integer i, jdiff, kdiff
      jdiff = regh0 - regl0 + 1
      kdiff = (regh1 - regl1 + 1) * jdiff
      do 10 i = kdiff + jdiff + 2, kdiff * (regh2 - regl2) - jdiff - 1
         r(i) = r(i) - alpha * w(i)
         x(i) = x(i) + alpha * p(i)
         z(i) = r(i) * c(i)
 10      rho = rho + mask(i) * z(i) * r(i)
      end

      subroutine hgcg2(p, z,
     @ dims(reg),
     @ alpha)
      integer dims(reg)
      real*8 p(dimdecu(reg))
      real*8 z(dimdecu(reg))
      real*8 alpha
      integer i, jdiff, kdiff
      jdiff = regh0 - regl0 + 1
      kdiff = (regh1 - regl1 + 1) * jdiff
      do 10 i = kdiff + jdiff + 2, kdiff * (regh2 - regl2) - jdiff - 1
 10      p(i) = alpha * p(i) + z(i)
      end

#endif
#if (CGOPT == 2)

      subroutine hgcg(r, p, z, x, w, c, mask,
     @ ngrids, strid1, strid2, nvals,
     @ rs, ps, zs, xs, ws, cs, ms,
     @ nsets, edge, nvals1, nvals2, dstart, sstart,
     @ dstrd1, dstrd2, sstrd1, sstrd2,
     @ hx, alpha, rho, it, pcode)
      implicit real*8(a-h,o-z), integer(i-n)
      real*8 r(0:1), p(0:1), z(0:1)
      real*8 x(0:1), w(0:1), c(0:1), mask(0:1)
      integer ngrids
      integer strid1(ngrids), strid2(ngrids), nvals(ngrids)
      integer rs(ngrids), ps(ngrids), zs(ngrids)
      integer xs(ngrids), ws(ngrids), cs(ngrids), ms(ngrids)
      real*8 edge(0:1)
      integer nsets
      integer nvals1(nsets), nvals2(nsets)
      integer dstart(nsets), sstart(nsets)
      integer dstrd1(nsets), dstrd2(nsets)
      integer sstrd1(nsets), sstrd2(nsets)
      real*8 hx, alpha, rho
      integer it, pcode
      hxm2 = 1.d0 / (hx*hx)
#ifdef HG_CROSS_STENCIL
      facm = -hxm2
#else
      facm = -hxm2 / 3.d0
#endif
      rho = 0.d0
      do 10 k = 1, ngrids
c must do this loop over entire array to initialize borders of p
         do 10 i = 0, nvals(k) - 1
            z(zs(k)+i) = r(rs(k)+i) * c(cs(k)+i)
            p(ps(k)+i) = z(zs(k)+i)
 10         rho = rho + mask(ms(k)+i) * z(zs(k)+i) * r(rs(k)+i)
      tol = 1.d-3 * rho
      if (tol .le. 0.d0) then
         return
      endif
 100  continue
      it = it + 1
      if (it .gt. 250) then
         print *, "FORTRAN hgcg---conjugate-gradient iteration failed"
      endif
      rhoold = rho
      alpha = 0.d0
      do 110 k = 1, nsets
         do 110 j = 0, nvals2(k) - 1
            do 110 i = 0, nvals1(k) - 1
 110           edge(dstart(k) + j * dstrd2(k) + i * dstrd1(k)) =
     @         edge(sstart(k) + j * sstrd2(k) + i * sstrd1(k))
      do 120 k = 1, ngrids
         do 120 i = strid1(k) + strid2(k) + 1,
     @              nvals(k) - strid1(k) - strid2(k) - 2
#ifdef HG_CROSS_STENCIL
            w(ws(k)+i) = facm *
     @        (p(ps(k)+i-1) + p(ps(k)+i+1) +
     @         p(ps(k)+i-strid1(k)) + p(ps(k)+i+strid1(k)) +
     @         p(ps(k)+i-strid2(k)) + p(ps(k)+i+strid2(k)) -
     @         6.d0 * p(ps(k)+i))
#else
            print *, "FORTRAN hgcg---full box stencils not implemented"
#endif
 120        alpha = alpha + mask(ms(k)+i) * p(ps(k)+i) * w(ws(k)+i)
      alpha = rho / alpha
      rho = 0.d0
      do 130 k = 1, ngrids
         do 130 i = strid1(k) + strid2(k) + 1,
     @              nvals(k) - strid1(k) - strid2(k) - 2
            r(rs(k)+i) = r(rs(k)+i) - alpha * w(ws(k)+i)
            x(xs(k)+i) = x(xs(k)+i) + alpha * p(ps(k)+i)
            z(zs(k)+i) = r(rs(k)+i) * c(cs(k)+i)
 130        rho = rho + mask(ms(k)+i) * z(zs(k)+i) * r(rs(k)+i)
      if (pcode .ge. 3) then
         print *, it, rho
      endif
      if (rho .le. tol .or. it .gt. 250) then
         return
      endif
      alpha = rho / rhoold
      do 140 k = 1, ngrids
         do 140 i = strid1(k) + strid2(k) + 1,
     @              nvals(k) - strid1(k) - strid2(k) - 2
 140        p(ps(k)+i) = alpha * p(ps(k)+i) + z(zs(k)+i)
      goto 100
      end

#endif
