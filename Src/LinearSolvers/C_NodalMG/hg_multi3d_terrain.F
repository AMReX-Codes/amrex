c twenty-seven point terrain version:

c-----------------------------------------------------------------------
c Note---assumes fdst linearly interpolated from cdst along face
      subroutine hgfres_terrain(
     & res,    resl0,resh0,resl1,resh1,resl2,resh2,
     & src,    srcl0,srch0,srcl1,srch1,srcl2,srch2,
     & fdst,   fdstl0,fdsth0,fdstl1,fdsth1,fdstl2,fdsth2,
     & cdst,   cdstl0,cdsth0,cdstl1,cdsth1,cdstl2,cdsth2,
     & sigmaf, sfl0,sfh0,sfl1,sfh1,sfl2,sfh2,
     & sigmac, scl0,sch0,scl1,sch1,scl2,sch2,
     &         regl0,regh0,regl1,regh1,regl2,regh2,
     & hx, hy, hz, ir, jr, kr, idim, idir)
      integer resl0,resh0,resl1,resh1,resl2,resh2
      integer srcl0,srch0,srcl1,srch1,srcl2,srch2
      integer fdstl0,fdsth0,fdstl1,fdsth1,fdstl2,fdsth2
      integer cdstl0,cdsth0,cdstl1,cdsth1,cdstl2,cdsth2
      integer sfl0,sfh0,sfl1,sfh1,sfl2,sfh2
      integer scl0,sch0,scl1,sch1,scl2,sch2
      integer regl0,regh0,regl1,regh1,regl2,regh2
      double precision hx, hy, hz
      double precision res(resl0:resh0,resl1:resh1,resl2:resh2)
      double precision src(srcl0:srch0,srcl1:srch1,srcl2:srch2)
      double precision fdst(fdstl0:fdsth0,fdstl1:fdsth1,fdstl2:fdsth2)
      double precision cdst(cdstl0:cdsth0,cdstl1:cdsth1,cdstl2:cdsth2)
      double precision sigmaf(sfl0:sfh0,sfl1:sfh1,sfl2:sfh2, 5)
      double precision sigmac(scl0:sch0,scl1:sch1,scl2:sch2, 5)
      integer ir, jr, kr, idim, idir
      double precision fac00, fac0, fac1, fac2, tmp, cen
      integer i, j, k, is, js, ks, l, m, n
      double precision MXXC, MYYC, MZZC, MXZC, MYZC
      double precision MXXF, MYYF, MZZF, MXZF, MYZF

      MXXC(ii,jj,kk) = (4.d0 *  cdst(ii,j,k)   + 
     &                   2.d0 * (cdst(ii,jj,k)  - cdst(i,jj,k) + 
     &                           cdst(ii,j,kk)  - cdst(i,j,kk)) + 
     &                          (cdst(ii,jj,kk) - cdst(i,jj,kk)))

      MYYC(ii,jj,kk) = (4.d0 *  cdst(i,jj,k)   + 
     &                   2.d0 * (cdst(ii,jj,k)  - cdst(ii,j,k) + 
     &                           cdst(i,jj,kk)  - cdst(i,j,kk)) + 
     &                          (cdst(ii,jj,kk) - cdst(ii,j,kk)))

      MZZC(ii,jj,kk) = (4.d0 *  cdst(i,j,kk)   + 
     &                   2.d0 * (cdst(ii,j,kk)  - cdst(ii,j,k) + 
     &                           cdst(i,jj,kk)  - cdst(i,jj,k)) + 
     &                          (cdst(ii,jj,kk) - cdst(ii,jj,k)))

      MXZC(ii,jj,kk) = (6.d0 *  cdst(ii,j,kk)  + 
     &                   3.d0 * (cdst(ii,jj,kk) - cdst(i,jj,k)))

      MYZC(ii,jj,kk) = (6.d0 *  cdst(i,jj,kk)  + 
     &                   3.d0 * (cdst(ii,jj,kk) - cdst(ii,j,k)))

      MXXF(i1,j1,k1,ii,jj,kk) = (4.d0 *  fdst(ii,j1,k1)   + 
     &                         2.d0 * (fdst(ii,jj,k1) - fdst(i1,jj,k1) +
     &                                fdst(ii,j1,kk) - fdst(i1,j1,kk)) +
     &                                (fdst(ii,jj,kk) - fdst(i1,jj,kk)))

      MYYF(i1,j1,k1,ii,jj,kk) = (4.d0 *  fdst(i1,jj,k1)   + 
     &                         2.d0 * (fdst(ii,jj,k1) - fdst(ii,j1,k1) +
     &                                fdst(i1,jj,kk) - fdst(i1,j1,kk)) +
     &                                (fdst(ii,jj,kk) - fdst(ii,j1,kk)))

      MZZF(i1,j1,k1,ii,jj,kk) = (4.d0 *  fdst(i1,j1,kk)   + 
     &                         2.d0 * (fdst(ii,j1,kk) - fdst(ii,j1,k1) +
     &                                fdst(i1,jj,kk) - fdst(i1,jj,k1)) +
     &                                (fdst(ii,jj,kk) - fdst(ii,jj,k1)))

      MXZF(i1,j1,k1,ii,jj,kk) = (6.d0 *  fdst(ii,j1,kk)  +
     &                         3.d0 * (fdst(ii,jj,kk) - fdst(i1,jj,k1)))

      MYZF(i1,j1,k1,ii,jj,kk) = (6.d0 *  fdst(i1,jj,kk)  +
     &                         3.d0 * (fdst(ii,jj,kk) - fdst(ii,j1,k1)))


      fac00 = 1.d0 / 36.d0

      if (idim .eq. 0) then
         i = regl0
         if (idir .eq. 1) then
            is = i - 1
         else
            is = i
         end if
         do k = regl2, regh2
            do j = regl1, regh1

                  tmp = sigmac(is ,j-1,k-1,1) * MXXC(i-idir,j-1,k-1)
     &                 +sigmac(is ,j-1,k  ,1) * MXXC(i-idir,j-1,k+1)
     &                 +sigmac(is ,j  ,k-1,1) * MXXC(i-idir,j+1,k-1)
     &                 +sigmac(is ,j  ,k  ,1) * MXXC(i-idir,j+1,k+1)
                  tmp = tmp
     &                 +sigmac(is ,j-1,k-1,2) * MYYC(i-idir,j-1,k-1)
     &                 +sigmac(is ,j-1,k  ,2) * MYYC(i-idir,j-1,k+1)
     &                 +sigmac(is ,j  ,k-1,2) * MYYC(i-idir,j+1,k-1)
     &                 +sigmac(is ,j  ,k  ,2) * MYYC(i-idir,j+1,k+1)
                  tmp = tmp
     &                 +sigmac(is ,j-1,k-1,3) * MZZC(i-idir,j-1,k-1)
     &                 +sigmac(is ,j-1,k  ,3) * MZZC(i-idir,j-1,k+1)
     &                 +sigmac(is ,j  ,k-1,3) * MZZC(i-idir,j+1,k-1)
     &                 +sigmac(is ,j  ,k  ,3) * MZZC(i-idir,j+1,k+1)
                  tmp = tmp + idir * (
     &                 -sigmac(is ,j-1,k-1,4) * MXZC(i-idir,j-1,k-1)
     &                 +sigmac(is ,j-1,k  ,4) * MXZC(i-idir,j-1,k+1)
     &                 -sigmac(is ,j  ,k-1,4) * MXZC(i-idir,j+1,k-1)
     &                 +sigmac(is ,j  ,k  ,4) * MXZC(i-idir,j+1,k+1) )
                  tmp = tmp
     &                 -sigmac(is ,j-1,k-1,5) * MYZC(i-idir,j-1,k-1)
     &                 +sigmac(is ,j-1,k  ,5) * MYZC(i-idir,j-1,k+1)
     &                 +sigmac(is ,j  ,k-1,5) * MYZC(i-idir,j+1,k-1)
     &                 -sigmac(is ,j  ,k  ,5) * MYZC(i-idir,j+1,k+1)

              cen        = 
     &            4.d0 * ( sigmac(is,j-1,k-1,1) + sigmac(is,j-1,k,1)
     &                    +sigmac(is,j  ,k-1,1) + sigmac(is,j  ,k,1)
     &                    +sigmac(is,j-1,k-1,2) + sigmac(is,j-1,k,2)
     &                    +sigmac(is,j  ,k-1,2) + sigmac(is,j  ,k,2)
     &                    +sigmac(is,j-1,k-1,3) + sigmac(is,j-1,k,3)
     &                    +sigmac(is,j  ,k-1,3) + sigmac(is,j  ,k,3)) 
     &   + idir * 6.d0 * (-sigmac(is,j-1,k-1,4) + sigmac(is,j-1,k,4)
     &                    -sigmac(is,j  ,k-1,4) + sigmac(is,j  ,k,4))
     &          + 6.d0 * ( sigmac(is,j-1,k  ,5) - sigmac(is,j-1,k-1,5)
     &                    +sigmac(is,j  ,k-1,5) - sigmac(is,j  ,k,5)) 

               res(i*ir,j*jr,k*kr) =
     &          src(i*ir,j*jr,k*kr) - fac00 * (tmp - cen * cdst(i,j,k))
            end do
         end do

         fac0 = fac00 / (jr * kr)
         i = i * ir
         if (idir .eq. 1) then
            is = i
         else
            is = i - 1
         end if
         do l = -(kr-1), kr-1
            fac2 = (kr-abs(l)) * fac0
            do n = -(jr-1), jr-1
               fac1 = (jr-abs(n)) * fac2
               do k = kr*regl2, kr*regh2, kr
                  do j = jr*regl1, jr*regh1, jr

            tmp = sigmaf(is,j+n-1,k+l-1,1) 
     &                * MXXF(i,j+n,k+l,i+idir,j+n-1,k+l-1)
     &           +sigmaf(is,j+n-1,k+l  ,1) 
     &                * MXXF(i,j+n,k+l,i+idir,j+n-1,k+l+1)
     &           +sigmaf(is,j+n  ,k+l-1,1) 
     &                * MXXF(i,j+n,k+l,i+idir,j+n+1,k+l-1)
     &           +sigmaf(is,j+n  ,k+l  ,1) 
     &                * MXXF(i,j+n,k+l,i+idir,j+n+1,k+l+1)
            tmp = tmp
     &           +sigmaf(is,j+n-1,k+l-1,2) 
     &                * MYYF(i,j+n,k+l,i+idir,j+n-1,k+l-1)
     &           +sigmaf(is,j+n-1,k+l  ,2) 
     &                * MYYF(i,j+n,k+l,i+idir,j+n-1,k+l+1)
     &           +sigmaf(is,j+n  ,k+l-1,2) 
     &                * MYYF(i,j+n,k+l,i+idir,j+n+1,k+l-1)
     &           +sigmaf(is,j+n  ,k+l  ,2) 
     &                * MYYF(i,j+n,k+l,i+idir,j+n+1,k+l+1)
            tmp = tmp
     &           +sigmaf(is,j+n-1,k+l-1,3) 
     &                * MZZF(i,j+n,k+l,i+idir,j+n-1,k+l-1)
     &           +sigmaf(is,j+n-1,k+l  ,3) 
     &                * MZZF(i,j+n,k+l,i+idir,j+n-1,k+l+1)
     &           +sigmaf(is,j+n  ,k+l-1,3) 
     &                * MZZF(i,j+n,k+l,i+idir,j+n+1,k+l-1)
     &           +sigmaf(is,j+n  ,k+l  ,3) 
     &                * MZZF(i,j+n,k+l,i+idir,j+n+1,k+l+1)
            tmp = tmp + idir * (
     &            sigmaf(is,j+n-1,k+l-1,4) 
     &                * MXZF(i,j+n,k+l,i+idir,j+n-1,k+l-1)
     &           -sigmaf(is,j+n-1,k+l  ,4) 
     &                * MXZF(i,j+n,k+l,i+idir,j+n-1,k+l+1)
     &           +sigmaf(is,j+n  ,k+l-1,4) 
     &                * MXZF(i,j+n,k+l,i+idir,j+n+1,k+l-1)
     &           -sigmaf(is,j+n  ,k+l  ,4) 
     &                * MXZF(i,j+n,k+l,i+idir,j+n+1,k+l+1))
            tmp = tmp
     &           -sigmaf(is,j+n-1,k+l-1,5) 
     &                * MYZF(i,j+n,k+l,i+idir,j+n-1,k+l-1)
     &           +sigmaf(is,j+n-1,k+l  ,5) 
     &                * MYZF(i,j+n,k+l,i+idir,j+n-1,k+l+1)
     &           +sigmaf(is,j+n  ,k+l-1,5) 
     &                * MYZF(i,j+n,k+l,i+idir,j+n+1,k+l-1)
     &           -sigmaf(is,j+n  ,k+l  ,5) 
     &                * MYZF(i,j+n,k+l,i+idir,j+n+1,k+l+1)

            cen =
     &       4.d0 *(sigmaf(is,j+n-1,k+l-1,1) + sigmaf(is,j+n-1,k+l  ,1)
     &             +sigmaf(is,j+n  ,k+l-1,1) + sigmaf(is,j+n  ,k+l  ,1)
     &             +sigmaf(is,j+n-1,k+l-1,2) + sigmaf(is,j+n-1,k+l  ,2)
     &             +sigmaf(is,j+n  ,k+l-1,2) + sigmaf(is,j+n  ,k+l  ,2)
     &             +sigmaf(is,j+n-1,k+l-1,3) + sigmaf(is,j+n-1,k+l  ,3)
     &             +sigmaf(is,j+n  ,k+l-1,3) + sigmaf(is,j+n  ,k+l  ,3))
     &     + idir * 
     &       6.d0 *(-sigmaf(is,j+n-1,k+l  ,4) + sigmaf(is,j+n-1,k+l-1,4)
     &             -sigmaf(is,j+n  ,k+l  ,4) + sigmaf(is,j+n  ,k+l-1,4))
     &     + 6.d0 *( sigmaf(is,j+n-1,k+l  ,5) - sigmaf(is,j+n-1,k+l-1,5)
     &             +sigmaf(is,j+n  ,k+l-1,5) - sigmaf(is,j+n  ,k+l  ,5))

             res(i,j,k) = res(i,j,k) 
     &            - fac1 * (tmp - cen * fdst(i,j+n,k+l))
                  end do
               end do
            end do
         end do
              
      else if (idim .eq. 1) then
         j = regl1
         if (idir .eq. 1) then
            js = j - 1
         else
            js = j
         end if
         do k = regl2, regh2
            do i = regl0, regh0

                  tmp = sigmac(i-1,js,k-1,1) * MXXC(i-1,j-idir,k-1)
     &                 +sigmac(i-1,js,k  ,1) * MXXC(i-1,j-idir,k+1)
     &                 +sigmac(i  ,js,k-1,1) * MXXC(i+1,j-idir,k-1)
     &                 +sigmac(i  ,js,k  ,1) * MXXC(i+1,j-idir,k+1)
                  tmp = tmp
     &                 +sigmac(i-1,js,k-1,2) * MYYC(i-1,j-idir,k-1)
     &                 +sigmac(i-1,js,k  ,2) * MYYC(i-1,j-idir,k+1)
     &                 +sigmac(i  ,js,k-1,2) * MYYC(i+1,j-idir,k-1)
     &                 +sigmac(i  ,js,k  ,2) * MYYC(i+1,j-idir,k+1)
                  tmp = tmp
     &                 +sigmac(i-1,js,k-1,3) * MZZC(i-1,j-idir,k-1)
     &                 +sigmac(i-1,js,k  ,3) * MZZC(i-1,j-idir,k+1)
     &                 +sigmac(i  ,js,k-1,3) * MZZC(i+1,j-idir,k-1)
     &                 +sigmac(i  ,js,k  ,3) * MZZC(i+1,j-idir,k+1)
                  tmp = tmp
     &                 -sigmac(i-1,js,k-1,4) * MXZC(i-1,j-idir,k-1)
     &                 +sigmac(i-1,js,k  ,4) * MXZC(i-1,j-idir,k+1)
     &                 +sigmac(i  ,js,k-1,4) * MXZC(i+1,j-idir,k-1)
     &                 -sigmac(i  ,js,k  ,4) * MXZC(i+1,j-idir,k+1)
                  tmp = tmp + idir * ( 
     &                 -sigmac(i-1,js,k-1,5) * MYZC(i-1,j-idir,k-1)
     &                 +sigmac(i-1,js,k  ,5) * MYZC(i-1,j-idir,k+1)
     &                 -sigmac(i  ,js,k-1,5) * MYZC(i+1,j-idir,k-1)
     &                 +sigmac(i  ,js,k  ,5) * MYZC(i+1,j-idir,k+1) )

               cen =
     &            4.d0 * (sigmac(i-1,js,k-1,1) + sigmac(i-1,js,k,1)
     &                   +sigmac(i  ,js,k-1,1) + sigmac(i  ,js,k,1)
     &                   +sigmac(i-1,js,k-1,2) + sigmac(i-1,js,k,2)
     &                   +sigmac(i  ,js,k-1,2) + sigmac(i  ,js,k,2)
     &                   +sigmac(i-1,js,k-1,3) + sigmac(i-1,js,k,3)
     &                   +sigmac(i  ,js,k-1,3) + sigmac(i  ,js,k,3))
     &          + 6.d0 * (sigmac(i-1,js,k  ,4) - sigmac(i-1,js,k-1,4)
     &                   +sigmac(i  ,js,k-1,4) - sigmac(i  ,js,k  ,4))
     &   + idir * 6.d0 * (-sigmac(i-1,js,k-1,5) + sigmac(i-1,js,k,5)
     &                    -sigmac(i  ,js,k-1,5) + sigmac(i  ,js,k,5))

               res(i*ir,j*jr,k*kr) =
     &          src(i*ir,j*jr,k*kr) - fac00 * (tmp - cen * cdst(i,j,k))
            end do
         end do

         fac0 = fac00 / (ir * kr)
         j = j * jr
         if (idir .eq. 1) then
            js = j
         else
            js = j - 1
         end if
         do l = -(kr-1), kr-1
            fac2 = (kr-abs(l)) * fac0
            do m = -(ir-1), ir-1
               fac1 = (ir-abs(m)) * fac2
               do k = kr*regl2, kr*regh2, kr
                  do i = ir*regl0, ir*regh0, ir

            tmp = sigmaf(i+m-1,js,k+l-1,1) 
     &                    * MXXF(i+m,j,k+l,i+m-1,j+idir,k+l-1)
     &           +sigmaf(i+m-1,js,k+l  ,1) 
     &                    * MXXF(i+m,j,k+l,i+m-1,j+idir,k+l+1)
     &           +sigmaf(i+m  ,js,k+l-1,1) 
     &                    * MXXF(i+m,j,k+l,i+m+1,j+idir,k+l-1)
     &           +sigmaf(i+m  ,js,k+l  ,1) 
     &                    * MXXF(i+m,j,k+l,i+m+1,j+idir,k+l+1)
            tmp = tmp
     &           +sigmaf(i+m-1,js,k+l-1,2) 
     &                    * MYYF(i+m,j,k+l,i+m-1,j+idir,k+l-1)
     &           +sigmaf(i+m-1,js,k+l  ,2) 
     &                    * MYYF(i+m,j,k+l,i+m-1,j+idir,k+l+1)
     &           +sigmaf(i+m  ,js,k+l-1,2) 
     &                    * MYYF(i+m,j,k+l,i+m+1,j+idir,k+l-1)
     &           +sigmaf(i+m  ,js,k+l  ,2) 
     &                    * MYYF(i+m,j,k+l,i+m+1,j+idir,k+l+1)
            tmp = tmp
     &           +sigmaf(i+m-1,js,k+l-1,3) 
     &                    * MZZF(i+m,j,k+l,i+m-1,j+idir,k+l-1)
     &           +sigmaf(i+m-1,js,k+l  ,3) 
     &                    * MZZF(i+m,j,k+l,i+m-1,j+idir,k+l+1)
     &           +sigmaf(i+m  ,js,k+l-1,3) 
     &                    * MZZF(i+m,j,k+l,i+m+1,j+idir,k+l-1)
     &           +sigmaf(i+m  ,js,k+l  ,3) 
     &                    * MZZF(i+m,j,k+l,i+m+1,j+idir,k+l+1)
            tmp = tmp
     &           -sigmaf(i+m-1,js,k+l-1,4) 
     &                    * MXZF(i+m,j,k+l,i+m-1,j+idir,k+l-1)
     &           +sigmaf(i+m-1,js,k+l  ,4) 
     &                    * MXZF(i+m,j,k+l,i+m-1,j+idir,k+l+1)
     &           +sigmaf(i+m  ,js,k+l-1,4) 
     &                    * MXZF(i+m,j,k+l,i+m+1,j+idir,k+l-1)
     &           -sigmaf(i+m  ,js,k+l  ,4) 
     &                    * MXZF(i+m,j,k+l,i+m+1,j+idir,k+l+1)
            tmp = tmp + idir * (
     &            sigmaf(i+m-1,js,k+l-1,5) 
     &                    * MYZF(i+m,j,k+l,i+m-1,j+idir,k+l-1)
     &           -sigmaf(i+m-1,js,k+l  ,5) 
     &                    * MYZF(i+m,j,k+l,i+m-1,j+idir,k+l+1)
     &           +sigmaf(i+m  ,js,k+l-1,5) 
     &                    * MYZF(i+m,j,k+l,i+m+1,j+idir,k+l-1)
     &           -sigmaf(i+m  ,js,k+l  ,5) 
     &                    * MYZF(i+m,j,k+l,i+m+1,j+idir,k+l+1))

         cen =
     &    4.d0 * (sigmaf(i+m-1,js,k+l-1,1) + sigmaf(i+m-1,js,k+l  ,1)
     &           +sigmaf(i+m  ,js,k+l-1,1) + sigmaf(i+m  ,js,k+l  ,1)
     &           +sigmaf(i+m-1,js,k+l-1,2) + sigmaf(i+m-1,js,k+l  ,2)
     &           +sigmaf(i+m  ,js,k+l-1,2) + sigmaf(i+m  ,js,k+l  ,2)
     &           +sigmaf(i+m-1,js,k+l-1,3) + sigmaf(i+m-1,js,k+l  ,3)
     &           +sigmaf(i+m  ,js,k+l-1,3) + sigmaf(i+m  ,js,k+l  ,3))
     &   +6.d0 * (sigmaf(i+m-1,js,k+l  ,4) - sigmaf(i+m-1,js,k+l-1,4)
     &           +sigmaf(i+m  ,js,k+l-1,4) - sigmaf(i+m  ,js,k+l  ,4))
     &   + idir * 
     &    6.d0 * (-sigmaf(i+m-1,js,k+l  ,5) + sigmaf(i+m-1,js,k+l-1,5)
     &            -sigmaf(i+m  ,js,k+l  ,5) + sigmaf(i+m  ,js,k+l-1,5))

          res(i,j,k) = res(i,j,k) - fac1 * (tmp - cen * fdst(i+m,j,k+l))
                  end do
               end do
            end do
         end do
      else
         k = regl2
         if (idir .eq. 1) then
            ks = k - 1
         else
            ks = k
         end if
         do j = regl1, regh1
            do i = regl0, regh0

                  tmp = sigmac(i-1,j-1,ks,1) * MXXC(i-1,j-1,k-idir)
     &                 +sigmac(i-1,j  ,ks,1) * MXXC(i-1,j+1,k-idir)
     &                 +sigmac(i  ,j-1,ks,1) * MXXC(i+1,j-1,k-idir)
     &                 +sigmac(i  ,j  ,ks,1) * MXXC(i+1,j+1,k-idir)
                  tmp = tmp
     &                 +sigmac(i-1,j-1,ks,2) * MYYC(i-1,j-1,k-idir)
     &                 +sigmac(i-1,j  ,ks,2) * MYYC(i-1,j+1,k-idir)
     &                 +sigmac(i  ,j-1,ks,2) * MYYC(i+1,j-1,k-idir)
     &                 +sigmac(i  ,j  ,ks,2) * MYYC(i+1,j+1,k-idir)
                  tmp = tmp
     &                 +sigmac(i-1,j-1,ks,3) * MZZC(i-1,j-1,k-idir)
     &                 +sigmac(i-1,j  ,ks,3) * MZZC(i-1,j+1,k-idir)
     &                 +sigmac(i  ,j-1,ks,3) * MZZC(i+1,j-1,k-idir)
     &                 +sigmac(i  ,j  ,ks,3) * MZZC(i+1,j+1,k-idir)
                  tmp = tmp + idir * (
     &                 -sigmac(i-1,j-1,ks,4) * MXZC(i-1,j-1,k-idir)
     &                 -sigmac(i-1,j  ,ks,4) * MXZC(i-1,j+1,k-idir)
     &                 +sigmac(i  ,j-1,ks,4) * MXZC(i+1,j-1,k-idir)
     &                 +sigmac(i  ,j  ,ks,4) * MXZC(i+1,j+1,k-idir) )
                  tmp = tmp + idir * (
     &                 -sigmac(i-1,j-1,ks,5) * MYZC(i-1,j-1,k-idir)
     &                 +sigmac(i-1,j  ,ks,5) * MYZC(i-1,j+1,k-idir)
     &                 -sigmac(i  ,j-1,ks,5) * MYZC(i+1,j-1,k-idir)
     &                 +sigmac(i  ,j  ,ks,5) * MYZC(i+1,j+1,k-idir) )

               cen = 
     &            4.d0 * ( sigmac(i-1,j-1,ks,1)
     &                    +sigmac(i-1,j  ,ks,1) 
     &                    +sigmac(i  ,j-1,ks,1)
     &                    +sigmac(i  ,j  ,ks,1)
     &                    +sigmac(i-1,j-1,ks,2)
     &                    +sigmac(i-1,j  ,ks,2)
     &                    +sigmac(i  ,j-1,ks,2)
     &                    +sigmac(i  ,j  ,ks,2)
     &                    +sigmac(i-1,j-1,ks,3)
     &                    +sigmac(i-1,j  ,ks,3)
     &                    +sigmac(i  ,j-1,ks,3)
     &                    +sigmac(i  ,j  ,ks,3))
     &   + idir * 6.d0 * (-sigmac(i-1,j-1,ks,4)   
     &                    -sigmac(i-1,j  ,ks,4)
     &                    +sigmac(i  ,j-1,ks,4) 
     &                    +sigmac(i  ,j  ,ks,4) )
     &   + idir * 6.d0 * (-sigmac(i-1,j-1,ks,5)  
     &                    +sigmac(i-1,j,ks,5)
     &                    -sigmac(i  ,j-1,ks,5)     
     &                    +sigmac(i  ,j,ks,5))

               res(i*ir,j*jr,k*kr) =
     &          src(i*ir,j*jr,k*kr) - fac00 * (tmp - cen * cdst(i,j,k))
            end do
         end do

         fac0 = fac00 / (ir * jr)
         k = k * kr
         if (idir .eq. 1) then
            ks = k
         else
            ks = k - 1
         end if
         do n = -(jr-1), jr-1
            fac2 = (jr-abs(n)) * fac0
            do m = -(ir-1), ir-1
               fac1 = (ir-abs(m)) * fac2
               do j = jr*regl1, jr*regh1, jr
                  do i = ir*regl0, ir*regh0, ir

            tmp = sigmaf(i+m-1,j+n-1,ks,1) 
     &                * MXXF(i+m,j+n,k,i+m-1,j+n-1,k+idir)
     &           +sigmaf(i+m-1,j+n  ,ks,1) 
     &                * MXXF(i+m,j+n,k,i+m-1,j+n+1,k+idir)
     &           +sigmaf(i+m  ,j+n-1,ks,1) 
     &                * MXXF(i+m,j+n,k,i+m+1,j+n-1,k+idir)
     &           +sigmaf(i+m  ,j+n  ,ks,1) 
     &                * MXXF(i+m,j+n,k,i+m+1,j+n+1,k+idir)
            tmp = tmp
     &           +sigmaf(i+m-1,j+n-1,ks,2) 
     &                * MYYF(i+m,j+n,k,i+m-1,j+n-1,k+idir)
     &           +sigmaf(i+m-1,j+n  ,ks,2) 
     &                * MYYF(i+m,j+n,k,i+m-1,j+n+1,k+idir)
     &           +sigmaf(i+m  ,j+n-1,ks,2) 
     &                * MYYF(i+m,j+n,k,i+m+1,j+n-1,k+idir)
     &           +sigmaf(i+m  ,j+n  ,ks,2) 
     &                * MYYF(i+m,j+n,k,i+m+1,j+n+1,k+idir)
            tmp = tmp
     &           +sigmaf(i+m-1,j+n-1,ks,3) 
     &                * MZZF(i+m,j+n,k,i+m-1,j+n-1,k+idir)
     &           +sigmaf(i+m-1,j+n  ,ks,3) 
     &                * MZZF(i+m,j+n,k,i+m-1,j+n+1,k+idir)
     &           +sigmaf(i+m  ,j+n-1,ks,3) 
     &                * MZZF(i+m,j+n,k,i+m+1,j+n-1,k+idir)
     &           +sigmaf(i+m  ,j+n  ,ks,3) 
     &                * MZZF(i+m,j+n,k,i+m+1,j+n+1,k+idir)
            tmp = tmp + idir * (
     &            sigmaf(i+m-1,j+n-1,ks,4) 
     &                * MXZF(i+m,j+n,k,i+m-1,j+n-1,k+idir)
     &           +sigmaf(i+m-1,j+n  ,ks,4) 
     &                * MXZF(i+m,j+n,k,i+m-1,j+n+1,k+idir)
     &           -sigmaf(i+m  ,j+n-1,ks,4) 
     &                * MXZF(i+m,j+n,k,i+m+1,j+n-1,k+idir)
     &           -sigmaf(i+m  ,j+n  ,ks,4) 
     &                * MXZF(i+m,j+n,k,i+m+1,j+n+1,k+idir))
            tmp = tmp + idir * (
     &            sigmaf(i+m-1,j+n-1,ks,5) 
     &                * MYZF(i+m,j+n,k,i+m-1,j+n-1,k+idir)
     &           -sigmaf(i+m-1,j+n  ,ks,5) 
     &                * MYZF(i+m,j+n,k,i+m-1,j+n+1,k+idir)
     &           +sigmaf(i+m  ,j+n-1,ks,5) 
     &                * MYZF(i+m,j+n,k,i+m+1,j+n-1,k+idir)
     &           -sigmaf(i+m  ,j+n  ,ks,5) 
     &                * MYZF(i+m,j+n,k,i+m+1,j+n+1,k+idir))

         cen =
     &      4.d0 * (sigmaf(i+m-1,j+n-1,ks,1)
     &             +sigmaf(i+m-1,j+n  ,ks,1)
     &             +sigmaf(i+m  ,j+n-1,ks,1)
     &             +sigmaf(i+m  ,j+n  ,ks,1)
     &             +sigmaf(i+m-1,j+n-1,ks,2)
     &             +sigmaf(i+m-1,j+n  ,ks,2)
     &             +sigmaf(i+m  ,j+n-1,ks,2)
     &             +sigmaf(i+m  ,j+n  ,ks,2)
     &             +sigmaf(i+m-1,j+n-1,ks,3)
     &             +sigmaf(i+m-1,j+n  ,ks,3)
     &             +sigmaf(i+m  ,j+n-1,ks,3)
     &             +sigmaf(i+m  ,j+n  ,ks,3))
     &    + idir * 
     &      6.d0 * (sigmaf(i+m-1,j+n-1,ks,4)
     &             +sigmaf(i+m-1,j+n  ,ks,4)
     &             -sigmaf(i+m  ,j+n-1,ks,4)
     &             -sigmaf(i+m  ,j+n  ,ks,4))
     &    + idir * 
     &    + 6.d0 * (sigmaf(i+m-1,j+n-1,ks,5)
     &             -sigmaf(i+m-1,j+n  ,ks,5)
     &             +sigmaf(i+m  ,j+n-1,ks,5)
     &             -sigmaf(i+m  ,j+n  ,ks,5))

          res(i,j,k) = res(i,j,k) - fac1 * (tmp - cen * fdst(i+m,j+n,k))
                  end do
               end do
            end do
         end do

      end if
      end
c-----------------------------------------------------------------------
c Note---assumes fdst linearly interpolated from cdst along face
      subroutine hgeres_terrain(
     & res,    resl0,resh0,resl1,resh1,resl2,resh2,
     & src,    srcl0,srch0,srcl1,srch1,srcl2,srch2,
     & fdst,   fdstl0,fdsth0,fdstl1,fdsth1,fdstl2,fdsth2,
     & cdst,   cdstl0,cdsth0,cdstl1,cdsth1,cdstl2,cdsth2,
     & sigmaf, sfl0,sfh0,sfl1,sfh1,sfl2,sfh2,
     & sigmac, scl0,sch0,scl1,sch1,scl2,sch2,
     &         regl0,regh0,regl1,regh1,regl2,regh2,
     & hx, hy, hz, ir, jr, kr, ivect, ga)
      integer resl0,resh0,resl1,resh1,resl2,resh2
      integer srcl0,srch0,srcl1,srch1,srcl2,srch2
      integer fdstl0,fdsth0,fdstl1,fdsth1,fdstl2,fdsth2
      integer cdstl0,cdsth0,cdstl1,cdsth1,cdstl2,cdsth2
      integer sfl0,sfh0,sfl1,sfh1,sfl2,sfh2
      integer scl0,sch0,scl1,sch1,scl2,sch2
      integer regl0,regh0,regl1,regh1,regl2,regh2
      double precision hx, hy, hz
      double precision res(resl0:resh0,resl1:resh1,resl2:resh2)
      double precision src(srcl0:srch0,srcl1:srch1,srcl2:srch2)
      double precision fdst(fdstl0:fdsth0,fdstl1:fdsth1,fdstl2:fdsth2)
      double precision cdst(cdstl0:cdsth0,cdstl1:cdsth1,cdstl2:cdsth2)
      double precision sigmaf(sfl0:sfh0,sfl1:sfh1,sfl2:sfh2,5)
      double precision sigmac(scl0:sch0,scl1:sch1,scl2:sch2,5)
      integer ir, jr, kr, ivect(0:2), ga(0:1,0:1,0:1)
      double precision fac00, fac0, fac1, fac, tmp, cen
      integer ic, jc, kc, if, jf, kf, ii, ji, ki, idir, jdir, kdir
      integer l, m, n
      double precision MXXC, MYYC, MZZC, MXZC, MYZC
      double precision MXXF, MYYF, MZZF, MXZF, MYZF
      MXXC(ii,jj,kk) = (4.d0 *  cdst(ii,jc,kc) + 
     &                   2.d0 * (cdst(ii,jj,kc) - cdst(ic,jj,kc) + 
     &                           cdst(ii,jc,kk) - cdst(ic,jc,kk)) +
     &                          (cdst(ii,jj,kk) - cdst(ic,jj,kk)))

      MYYC(ii,jj,kk) = (4.d0 *  cdst(ic,jj,kc) + 
     &                   2.d0 * (cdst(ii,jj,kc) - cdst(ii,jc,kc) + 
     &                           cdst(ic,jj,kk) - cdst(ic,jc,kk)) +
     &                          (cdst(ii,jj,kk) - cdst(ii,jc,kk)))

      MZZC(ii,jj,kk) = (4.d0 *  cdst(ic,jc,kk) + 
     &                   2.d0 * (cdst(ii,jc,kk) - cdst(ii,jc,kc) + 
     &                           cdst(ic,jj,kk) - cdst(ic,jj,kc)) +
     &                          (cdst(ii,jj,kk) - cdst(ii,jj,kc)))

      MXZC(ii,jj,kk) = (6.d0 *  cdst(ii,jc,kk) +
     &                   3.d0 * (cdst(ii,jj,kk) - cdst(ic,jj,kc)))

      MYZC(ii,jj,kk) = (6.d0 *  cdst(ic,jj,kk) +
     &                   3.d0 * (cdst(ii,jj,kk) - cdst(ii,jc,kc)))

      MXXF(i1,j1,k1,ii,jj,kk) = (4.d0 *  fdst(ii,j1,k1)   + 
     &                         2.d0 * (fdst(ii,jj,k1) - fdst(i1,jj,k1) +
     &                                fdst(ii,j1,kk) - fdst(i1,j1,kk)) +
     &                                (fdst(ii,jj,kk) - fdst(i1,jj,kk)))

      MYYF(i1,j1,k1,ii,jj,kk) = (4.d0 *  fdst(i1,jj,k1)   + 
     &                         2.d0 * (fdst(ii,jj,k1) - fdst(ii,j1,k1) +
     &                                fdst(i1,jj,kk) - fdst(i1,j1,kk)) +
     &                                (fdst(ii,jj,kk) - fdst(ii,j1,kk)))

      MZZF(i1,j1,k1,ii,jj,kk) = (4.d0 *  fdst(i1,j1,kk)   + 
     &                         2.d0 * (fdst(ii,j1,kk) - fdst(ii,j1,k1) +
     &                                fdst(i1,jj,kk) - fdst(i1,jj,k1)) +
     &                                (fdst(ii,jj,kk) - fdst(ii,jj,k1)))

      MXZF(i1,j1,k1,ii,jj,kk) = (6.d0 *  fdst(ii,j1,kk)  +
     &                         3.d0 * (fdst(ii,jj,kk) - fdst(i1,jj,k1)))

      MYZF(i1,j1,k1,ii,jj,kk) = (6.d0 *  fdst(i1,jj,kk)  +
     &                         3.d0 * (fdst(ii,jj,kk) - fdst(ii,j1,k1)))


      ic = regl0
      jc = regl1
      kc = regl2
      if = ic * ir
      jf = jc * jr
      kf = kc * kr

      fac00 = 1.d0 / 36.d0

      if (ivect(0) .eq. 0) then
         do if = ir*regl0, ir*regh0, ir
            res(if,jf,kf) = 0.d0
         end do
c quadrants
c each quadrant is two octants and their share of the two central edges
         fac0 = 1.d0 / ir
         do ki = 0, 1
            kdir = 2 * ki - 1
            do ji = 0, 1
               jdir = 2 * ji - 1
               if (ga(0,ji,ki) .eq. 1) then
                  do m = -(ir-1), ir-1
                     fac = (ir-abs(m)) * fac0
                     do if = ir*regl0, ir*regh0, ir

                  tmp = sigmaf(if+m-1,jf+ji-1,kf+ki-1,1) * 
     &                    MXXF(if+m,jf,kf,if+m-1,jf+jdir,kf+kdir)
     &                 +sigmaf(if+m  ,jf+ji-1,kf+ki-1,1) * 
     &                    MXXF(if+m,jf,kf,if+m+1,jf+jdir,kf+kdir)
                  tmp = tmp
     &                 +sigmaf(if+m-1,jf+ji-1,kf+ki-1,2) * 
     &                    MYYF(if+m,jf,kf,if+m-1,jf+jdir,kf+kdir)
     &                 +sigmaf(if+m  ,jf+ji-1,kf+ki-1,2) * 
     &                    MYYF(if+m,jf,kf,if+m+1,jf+jdir,kf+kdir)
                  tmp = tmp 
     &                 +sigmaf(if+m-1,jf+ji-1,kf+ki-1,3) * 
     &                    MZZF(if+m,jf,kf,if+m-1,jf+jdir,kf+kdir)
     &                 +sigmaf(if+m  ,jf+ji-1,kf+ki-1,3) * 
     &                    MZZF(if+m,jf,kf,if+m+1,jf+jdir,kf+kdir)
                  tmp = tmp + kdir * (
     &                  sigmaf(if+m-1,jf+ji-1,kf+ki-1,4) * 
     &                    MXZF(if+m,jf,kf,if+m-1,jf+jdir,kf+kdir)
     &                 -sigmaf(if+m  ,jf+ji-1,kf+ki-1,4) * 
     &                    MXZF(if+m,jf,kf,if+m+1,jf+jdir,kf+kdir) )
                  tmp = tmp + jdir * kdir * (
     &                 -sigmaf(if+m-1,jf+ji-1,kf+ki-1,5) * 
     &                    MYZF(if+m,jf,kf,if+m-1,jf+jdir,kf+kdir)
     &                 -sigmaf(if+m  ,jf+ji-1,kf+ki-1,5) * 
     &                    MYZF(if+m,jf,kf,if+m+1,jf+jdir,kf+kdir) )

                  cen =
     &               4.d0 * ( sigmaf(if+m-1,jf+ji-1,kf+ki-1,1)
     &                       +sigmaf(if+m  ,jf+ji-1,kf+ki-1,1)
     &                       +sigmaf(if+m-1,jf+ji-1,kf+ki-1,2)
     &                       +sigmaf(if+m  ,jf+ji-1,kf+ki-1,2)
     &                       +sigmaf(if+m-1,jf+ji-1,kf+ki-1,3)
     &                       +sigmaf(if+m  ,jf+ji-1,kf+ki-1,3))
     &             + 6.d0 * ( sigmaf(if+m-1,jf+ji-1,kf+ki-1,4)
     &                       -sigmaf(if+m  ,jf+ji-1,kf+ki-1,4)) * kdir 
     &             + 6.d0 * (-sigmaf(if+m-1,jf+ji-1,kf+ki-1,5)
     &                  -sigmaf(if+m  ,jf+ji-1,kf+ki-1,5)) * jdir * kdir

                  tmp = tmp - cen * fdst(if+m,jf,kf)

                        res(if,jf,kf) = res(if,jf,kf) + fac * tmp
                     end do
                  end do

               else
                  do ic = regl0, regh0
                     if = ic * ir

                  tmp = sigmac(ic-1,jc+ji-1,kc+ki-1,1) * 
     &                    MXXC(ic-1,jc+jdir,kc+kdir)
     &                 +sigmac(ic  ,jc+ji-1,kc+ki-1,1) * 
     &                    MXXC(ic+1,jc+jdir,kc+kdir)
                  tmp = tmp
     &                 +sigmac(ic-1,jc+ji-1,kc+ki-1,2) * 
     &                    MYYC(ic-1,jc+jdir,kc+kdir)
     &                 +sigmac(ic  ,jc+ji-1,kc+ki-1,2) * 
     &                    MYYC(ic+1,jc+jdir,kc+kdir)
                  tmp = tmp
     &                 +sigmac(ic-1,jc+ji-1,kc+ki-1,3) * 
     &                    MZZC(ic-1,jc+jdir,kc+kdir)
     &                 +sigmac(ic  ,jc+ji-1,kc+ki-1,3) * 
     &                    MZZC(ic+1,jc+jdir,kc+kdir)
                  tmp = tmp + kdir * (
     &                  sigmac(ic-1,jc+ji-1,kc+ki-1,4) * 
     &                    MXZC(ic-1,jc+jdir,kc+kdir)
     &                 -sigmac(ic  ,jc+ji-1,kc+ki-1,4) * 
     &                    MXZC(ic+1,jc+jdir,kc+kdir) )
                  tmp = tmp + jdir * kdir * ( 
     &                 -sigmac(ic-1,jc+ji-1,kc+ki-1,5) * 
     &                    MYZC(ic-1,jc+jdir,kc+kdir)
     &                 -sigmac(ic  ,jc+ji-1,kc+ki-1,5) * 
     &                    MYZC(ic+1,jc+jdir,kc+kdir) )
                   cen =
     &               4.d0 * ( sigmac(ic-1,jc+ji-1,kc+ki-1,1)
     &                       +sigmac(ic  ,jc+ji-1,kc+ki-1,1)
     &                       +sigmac(ic-1,jc+ji-1,kc+ki-1,2)
     &                       +sigmac(ic  ,jc+ji-1,kc+ki-1,2)
     &                       +sigmac(ic-1,jc+ji-1,kc+ki-1,3)
     &                       +sigmac(ic  ,jc+ji-1,kc+ki-1,3))
     &             + 6.d0 * ( sigmac(ic-1,jc+ji-1,kc+ki-1,4)
     &                       -sigmac(ic  ,jc+ji-1,kc+ki-1,4)) * kdir
     &             + 6.d0 * (-sigmac(ic-1,jc+ji-1,kc+ki-1,5)
     &                    -sigmac(ic  ,jc+ji-1,kc+ki-1,5)) * jdir * kdir

                     tmp = tmp - cen * cdst(ic,jc,kc)
                     res(if,jf,kf) = res(if,jf,kf) + tmp
                  end do
               end if
            end do
         end do
c faces
c each face is two faces and two sides of an edge
         do ki = 0, 1
            kdir = 2 * ki - 1
            do ji = 0, 1
               jdir = 2 * ji - 1

               if (ga(0,ji,ki) - ga(0,ji,1-ki) .eq. 1) then
                  fac0 = 1.d0 / (ir * jr)
                  do  n = jdir, jdir*(jr-1), jdir
                     fac1 = (jr-abs(n)) * fac0
                     do m = -(ir-1), ir-1
                        fac = (ir-abs(m)) * fac1
                        do if = ir*regl0, ir*regh0, ir

                  tmp = sigmaf(if+m-1,jf+n-1,kf+ki-1,1) * 
     &                    MXXF(if+m,jf+n,kf,if+m-1,jf+n-1,kf+kdir)
     &                 +sigmaf(if+m-1,jf+n  ,kf+ki-1,1) * 
     &                    MXXF(if+m,jf+n,kf,if+m-1,jf+n+1,kf+kdir)
     &                 +sigmaf(if+m  ,jf+n-1,kf+ki-1,1) * 
     &                    MXXF(if+m,jf+n,kf,if+m+1,jf+n-1,kf+kdir)
     &                 +sigmaf(if+m  ,jf+n  ,kf+ki-1,1) * 
     &                    MXXF(if+m,jf+n,kf,if+m+1,jf+n+1,kf+kdir)
                  tmp = tmp
     &                 +sigmaf(if+m-1,jf+n-1,kf+ki-1,2) * 
     &                    MYYF(if+m,jf+n,kf,if+m-1,jf+n-1,kf+kdir)
     &                 +sigmaf(if+m-1,jf+n  ,kf+ki-1,2) * 
     &                    MYYF(if+m,jf+n,kf,if+m-1,jf+n+1,kf+kdir)
     &                 +sigmaf(if+m  ,jf+n-1,kf+ki-1,2) * 
     &                    MYYF(if+m,jf+n,kf,if+m+1,jf+n-1,kf+kdir)
     &                 +sigmaf(if+m  ,jf+n  ,kf+ki-1,2) * 
     &                    MYYF(if+m,jf+n,kf,if+m+1,jf+n+1,kf+kdir)
                  tmp = tmp
     &                 +sigmaf(if+m-1,jf+n-1,kf+ki-1,3) * 
     &                    MZZF(if+m,jf+n,kf,if+m-1,jf+n-1,kf+kdir)
     &                 +sigmaf(if+m-1,jf+n  ,kf+ki-1,3) * 
     &                    MZZF(if+m,jf+n,kf,if+m-1,jf+n+1,kf+kdir)
     &                 +sigmaf(if+m  ,jf+n-1,kf+ki-1,3) * 
     &                    MZZF(if+m,jf+n,kf,if+m+1,jf+n-1,kf+kdir)
     &                 +sigmaf(if+m  ,jf+n  ,kf+ki-1,3) * 
     &                    MZZF(if+m,jf+n,kf,if+m+1,jf+n+1,kf+kdir)
                  tmp = tmp + kdir * (
     &                 +sigmaf(if+m-1,jf+n-1,kf+ki-1,4) * 
     &                    MXZF(if+m,jf+n,kf,if+m-1,jf+n-1,kf+kdir)
     &                 +sigmaf(if+m-1,jf+n  ,kf+ki-1,4) * 
     &                    MXZF(if+m,jf+n,kf,if+m-1,jf+n+1,kf+kdir)
     &                 -sigmaf(if+m  ,jf+n-1,kf+ki-1,4) * 
     &                    MXZF(if+m,jf+n,kf,if+m+1,jf+n-1,kf+kdir)
     &                 -sigmaf(if+m  ,jf+n  ,kf+ki-1,4) * 
     &                    MXZF(if+m,jf+n,kf,if+m+1,jf+n+1,kf+kdir) )
                  tmp = tmp + kdir * ( 
     &                 +sigmaf(if+m-1,jf+n-1,kf+ki-1,5) * 
     &                    MYZF(if+m,jf+n,kf,if+m-1,jf+n-1,kf+kdir)
     &                 -sigmaf(if+m-1,jf+n  ,kf+ki-1,5) * 
     &                    MYZF(if+m,jf+n,kf,if+m-1,jf+n+1,kf+kdir)
     &                 +sigmaf(if+m  ,jf+n-1,kf+ki-1,5) * 
     &                    MYZF(if+m,jf+n,kf,if+m+1,jf+n-1,kf+kdir)
     &                 -sigmaf(if+m  ,jf+n  ,kf+ki-1,5) * 
     &                    MYZF(if+m,jf+n,kf,if+m+1,jf+n+1,kf+kdir) )

         cen =
     &      4.d0 * (sigmaf(if+m-1,jf+n-1,kf+ki-1,1)
     &             +sigmaf(if+m-1,jf+n  ,kf+ki-1,1)
     &             +sigmaf(if+m  ,jf+n-1,kf+ki-1,1)
     &             +sigmaf(if+m  ,jf+n  ,kf+ki-1,1)
     &             +sigmaf(if+m-1,jf+n-1,kf+ki-1,2)
     &             +sigmaf(if+m-1,jf+n  ,kf+ki-1,2)
     &             +sigmaf(if+m  ,jf+n-1,kf+ki-1,2)
     &             +sigmaf(if+m  ,jf+n  ,kf+ki-1,2)
     &             +sigmaf(if+m-1,jf+n-1,kf+ki-1,3)
     &             +sigmaf(if+m-1,jf+n  ,kf+ki-1,3)
     &             +sigmaf(if+m  ,jf+n-1,kf+ki-1,3)
     &             +sigmaf(if+m  ,jf+n  ,kf+ki-1,3))
     &    + 6.d0 * (sigmaf(if+m-1,jf+n-1,kf+ki-1,4)
     &             +sigmaf(if+m-1,jf+n  ,kf+ki-1,4)
     &             -sigmaf(if+m  ,jf+n-1,kf+ki-1,4)
     &             -sigmaf(if+m  ,jf+n  ,kf+ki-1,4)) * kdir
     &    + 6.d0 * (sigmaf(if+m-1,jf+n-1,kf+ki-1,5)
     &             -sigmaf(if+m-1,jf+n  ,kf+ki-1,5)
     &             +sigmaf(if+m  ,jf+n-1,kf+ki-1,5)
     &             -sigmaf(if+m  ,jf+n  ,kf+ki-1,5)) * kdir

                  tmp = tmp - cen * fdst(if+m,jf+n,kf)
                  res(if,jf,kf) = res(if,jf,kf) + fac * tmp
                        end do
                     end do
                  end do
               end if

               if (ga(0,ji,ki) - ga(0,1-ji,ki) .eq. 1) then
                  fac0 = 1.d0 / (ir * kr)
                  do l = kdir, kdir*(kr-1), kdir
                     fac1 = (kr-abs(l)) * fac0
                     do m = -(ir-1), ir-1
                        fac = (ir-abs(m)) * fac1
                        do if = ir*regl0, ir*regh0, ir

                  tmp = sigmaf(if+m-1,jf+ji-1,kf+l-1,1) * 
     &                    MXXF(if+m,jf,kf+l,if+m-1,jf+jdir,kf+l-1)
     &                 +sigmaf(if+m-1,jf+ji-1,kf+l  ,1) * 
     &                    MXXF(if+m,jf,kf+l,if+m-1,jf+jdir,kf+l+1)
     &                 +sigmaf(if+m  ,jf+ji-1,kf+l-1,1) * 
     &                    MXXF(if+m,jf,kf+l,if+m+1,jf+jdir,kf+l-1)
     &                 +sigmaf(if+m  ,jf+ji-1,kf+l  ,1) * 
     &                    MXXF(if+m,jf,kf+l,if+m+1,jf+jdir,kf+l+1)
                  tmp = tmp
     &                 +sigmaf(if+m-1,jf+ji-1,kf+l-1,2) * 
     &                    MYYF(if+m,jf,kf+l,if+m-1,jf+jdir,kf+l-1)
     &                 +sigmaf(if+m-1,jf+ji-1,kf+l  ,2) * 
     &                    MYYF(if+m,jf,kf+l,if+m-1,jf+jdir,kf+l+1)
     &                 +sigmaf(if+m  ,jf+ji-1,kf+l-1,2) * 
     &                    MYYF(if+m,jf,kf+l,if+m+1,jf+jdir,kf+l-1)
     &                 +sigmaf(if+m  ,jf+ji-1,kf+l  ,2) * 
     &                    MYYF(if+m,jf,kf+l,if+m+1,jf+jdir,kf+l+1)
                  tmp = tmp
     &                 +sigmaf(if+m-1,jf+ji-1,kf+l-1,3) * 
     &                    MZZF(if+m,jf,kf+l,if+m-1,jf+jdir,kf+l-1)
     &                 +sigmaf(if+m-1,jf+ji-1,kf+l  ,3) * 
     &                    MZZF(if+m,jf,kf+l,if+m-1,jf+jdir,kf+l+1)
     &                 +sigmaf(if+m  ,jf+ji-1,kf+l-1,3) * 
     &                    MZZF(if+m,jf,kf+l,if+m+1,jf+jdir,kf+l-1)
     &                 +sigmaf(if+m  ,jf+ji-1,kf+l  ,3) * 
     &                    MZZF(if+m,jf,kf+l,if+m+1,jf+jdir,kf+l+1)
                  tmp = tmp + ( 
     &                 -sigmaf(if+m-1,jf+ji-1,kf+l-1,4) * 
     &                    MXZF(if+m,jf,kf+l,if+m-1,jf+jdir,kf+l-1)
     &                 +sigmaf(if+m-1,jf+ji-1,kf+l  ,4) * 
     &                    MXZF(if+m,jf,kf+l,if+m-1,jf+jdir,kf+l+1)
     &                 +sigmaf(if+m  ,jf+ji-1,kf+l-1,4) * 
     &                    MXZF(if+m,jf,kf+l,if+m+1,jf+jdir,kf+l-1)
     &                 -sigmaf(if+m  ,jf+ji-1,kf+l  ,4) * 
     &                    MXZF(if+m,jf,kf+l,if+m+1,jf+jdir,kf+l+1) )
                  tmp = tmp + jdir * (
     &                 +sigmaf(if+m-1,jf+ji-1,kf+l-1,5) * 
     &                    MYZF(if+m,jf,kf+l,if+m-1,jf+jdir,kf+l-1)
     &                 -sigmaf(if+m-1,jf+ji-1,kf+l  ,5) * 
     &                    MYZF(if+m,jf,kf+l,if+m-1,jf+jdir,kf+l+1)
     &                 +sigmaf(if+m  ,jf+ji-1,kf+l-1,5) * 
     &                    MYZF(if+m,jf,kf+l,if+m+1,jf+jdir,kf+l-1)
     &                 -sigmaf(if+m  ,jf+ji-1,kf+l  ,5) * 
     &                    MYZF(if+m,jf,kf+l,if+m+1,jf+jdir,kf+l+1) )
         cen =
     &      4.d0 * (
     &             +sigmaf(if+m-1,jf+ji-1,kf+l-1,1)  
     &             +sigmaf(if+m-1,jf+ji-1,kf+l  ,1)
     &             +sigmaf(if+m  ,jf+ji-1,kf+l-1,1)  
     &             +sigmaf(if+m  ,jf+ji-1,kf+l  ,1)
     &             +sigmaf(if+m-1,jf+ji-1,kf+l-1,2)  
     &             +sigmaf(if+m-1,jf+ji-1,kf+l  ,2)
     &             +sigmaf(if+m  ,jf+ji-1,kf+l-1,2)  
     &             +sigmaf(if+m  ,jf+ji-1,kf+l  ,2)
     &             +sigmaf(if+m-1,jf+ji-1,kf+l-1,3)  
     &             +sigmaf(if+m-1,jf+ji-1,kf+l  ,3)
     &             +sigmaf(if+m  ,jf+ji-1,kf+l-1,3)  
     &             +sigmaf(if+m  ,jf+ji-1,kf+l  ,3))
     &    + 6.d0 * (sigmaf(if+m-1,jf+ji-1,kf+l  ,4)  
     &             -sigmaf(if+m-1,jf+ji-1,kf+l-1,4)
     &             +sigmaf(if+m  ,jf+ji-1,kf+l-1,4) 
     &             -sigmaf(if+m  ,jf+ji-1,kf+l  ,4))
     &    + 6.d0 * (sigmaf(if+m-1,jf+ji-1,kf+l-1,5)  
     &             -sigmaf(if+m-1,jf+ji-1,kf+l  ,5)
     &             +sigmaf(if+m  ,jf+ji-1,kf+l-1,5)  
     &             -sigmaf(if+m  ,jf+ji-1,kf+l  ,5)) * jdir

                  tmp = tmp - cen * fdst(if+m,jf,kf+l)
                  res(if,jf,kf) = res(if,jf,kf) + fac * tmp
                        end do
                     end do
                  end do
               end if
            end do
         end do
c weighting
         do if = ir*regl0, ir*regh0, ir
            res(if,jf,kf) = src(if,jf,kf) - res(if,jf,kf) * fac00
         end do

      else if (ivect(1) .eq. 0) then
         do jf = jr*regl1, jr*regh1, jr
            res(if,jf,kf) = 0.d0
         end do
c quadrants
c each quadrant is two octants and their share of the two central edges
         fac0 = 1.d0 / jr
         do  ki = 0, 1
            kdir = 2 * ki - 1
            do ii = 0, 1
               idir = 2 * ii - 1
               if (ga(ii,0,ki) .eq. 1) then
                  do n = -(jr-1), jr-1
                     fac = (jr-abs(n)) * fac0
                     do jf = jr*regl1, jr*regh1, jr

                  tmp = sigmaf(if+ii-1,jf+n-1,kf+ki-1,1) * 
     &                    MXXF(if,jf+n,kf,if+idir,jf+n-1,kf+kdir)
     &                 +sigmaf(if+ii-1,jf+n  ,kf+ki-1,1) * 
     &                    MXXF(if,jf+n,kf,if+idir,jf+n+1,kf+kdir)
                  tmp = tmp
     &                 +sigmaf(if+ii-1,jf+n-1,kf+ki-1,2) * 
     &                    MYYF(if,jf+n,kf,if+idir,jf+n-1,kf+kdir)
     &                 +sigmaf(if+ii-1,jf+n  ,kf+ki-1,2) * 
     &                    MYYF(if,jf+n,kf,if+idir,jf+n+1,kf+kdir) 
                  tmp = tmp 
     &                 +sigmaf(if+ii-1,jf+n-1,kf+ki-1,3) * 
     &                    MZZF(if,jf+n,kf,if+idir,jf+n-1,kf+kdir)
     &                 +sigmaf(if+ii-1,jf+n  ,kf+ki-1,3) * 
     &                    MZZF(if,jf+n,kf,if+idir,jf+n+1,kf+kdir)
                  tmp = tmp + idir * kdir * (
     &                 -sigmaf(if+ii-1,jf+n-1,kf+ki-1,4) * 
     &                    MXZF(if,jf+n,kf,if+idir,jf+n-1,kf+kdir)
     &                 -sigmaf(if+ii-1,jf+n  ,kf+ki-1,4) * 
     &                    MXZF(if,jf+n,kf,if+idir,jf+n+1,kf+kdir) )
                  tmp = tmp + kdir * (
     &                 +sigmaf(if+ii-1,jf+n-1,kf+ki-1,5) * 
     &                    MYZF(if,jf+n,kf,if+idir,jf+n-1,kf+kdir)
     &                 -sigmaf(if+ii-1,jf+n  ,kf+ki-1,5) * 
     &                    MYZF(if,jf+n,kf,if+idir,jf+n+1,kf+kdir) )
                  cen =
     &                  4.d0 * ( sigmaf(if+ii-1,jf+n-1,kf+ki-1,1)
     &                          +sigmaf(if+ii-1,jf+n  ,kf+ki-1,1)
     &                          +sigmaf(if+ii-1,jf+n-1,kf+ki-1,2)
     &                          +sigmaf(if+ii-1,jf+n  ,kf+ki-1,2)
     &                          +sigmaf(if+ii-1,jf+n-1,kf+ki-1,3)
     &                          +sigmaf(if+ii-1,jf+n  ,kf+ki-1,3))
     &                + 6.d0 * (-sigmaf(if+ii-1,jf+n-1,kf+ki-1,4)
     &                  -sigmaf(if+ii-1,jf+n  ,kf+ki-1,4)) * idir * kdir
     &                + 6.d0 * ( sigmaf(if+ii-1,jf+n-1,kf+ki-1,5) 
     &                         -sigmaf(if+ii-1,jf+n  ,kf+ki-1,5)) * kdir

                  tmp = tmp - cen * fdst(if,jf+n,kf)

                  res(if,jf,kf) = res(if,jf,kf) + fac * tmp
                     end do
                  end do
               else
                  do jc = regl1, regh1
                     jf = jc * jr

                  tmp = sigmac(ic+ii-1,jc-1,kc+ki-1,1) * 
     &                    MXXC(ic+idir,jc-1,kc+kdir)
     &                 +sigmac(ic+ii-1,jc  ,kc+ki-1,1) * 
     &                    MXXC(ic+idir,jc+1,kc+kdir)
                  tmp = tmp
     &                 +sigmac(ic+ii-1,jc-1,kc+ki-1,2) * 
     &                    MYYC(ic+idir,jc-1,kc+kdir)
     &                 +sigmac(ic+ii-1,jc  ,kc+ki-1,2) * 
     &                    MYYC(ic+idir,jc+1,kc+kdir)
                  tmp = tmp
     &                 +sigmac(ic+ii-1,jc-1,kc+ki-1,3) * 
     &                    MZZC(ic+idir,jc-1,kc+kdir)
     &                 +sigmac(ic+ii-1,jc  ,kc+ki-1,3) * 
     &                    MZZC(ic+idir,jc+1,kc+kdir)
                  tmp = tmp + idir * kdir * (
     &                 -sigmac(ic+ii-1,jc-1,kc+ki-1,4) * 
     &                    MXZC(ic+idir,jc-1,kc+kdir)
     &                 -sigmac(ic+ii-1,jc  ,kc+ki-1,4) * 
     &                    MXZC(ic+idir,jc+1,kc+kdir) )
                  tmp = tmp + kdir * (
     &                 +sigmac(ic+ii-1,jc-1,kc+ki-1,5) * 
     &                    MYZC(ic+idir,jc-1,kc+kdir)
     &                 -sigmac(ic+ii-1,jc  ,kc+ki-1,5) * 
     &                    MYZC(ic+idir,jc+1,kc+kdir) )
                  cen =
     &               4.d0 * (
     &                       +sigmac(ic+ii-1,jc-1,kc+ki-1,1)
     &                       +sigmac(ic+ii-1,jc  ,kc+ki-1,1)
     &                       +sigmac(ic+ii-1,jc-1,kc+ki-1,2)
     &                       +sigmac(ic+ii-1,jc  ,kc+ki-1,2)
     &                       +sigmac(ic+ii-1,jc-1,kc+ki-1,3)
     &                       +sigmac(ic+ii-1,jc  ,kc+ki-1,3))
     &             + 6.d0 * (-sigmac(ic+ii-1,jc-1,kc+ki-1,4)
     &                    -sigmac(ic+ii-1,jc  ,kc+ki-1,4)) * idir * kdir
     &             + 6.d0 * ( sigmac(ic+ii-1,jc-1,kc+ki-1,5) 
     &                       -sigmac(ic+ii-1,jc  ,kc+ki-1,5)) * kdir

                     tmp = tmp - cen * cdst(ic,jc,kc)
                     res(if,jf,kf) = res(if,jf,kf) + tmp
                  end do
               end if
            end do
         end do
c faces
c each face is two faces and two sides of an edge
         do ki = 0, 1
            kdir = 2 * ki - 1
            do ii = 0, 1
               idir = 2 * ii - 1
               if (ga(ii,0,ki) - ga(ii,0,1-ki) .eq. 1) then
                  fac0 = 1.d0 / (ir * jr)
                  do n = -(jr-1), jr-1
                     fac1 = (jr-abs(n)) * fac0
                     do m = idir, idir*(ir-1), idir
                        fac = (ir-abs(m)) * fac1
                        do jf = jr*regl1, jr*regh1, jr

                  tmp = sigmaf(if+m-1,jf+n-1,kf+ki-1,1) * 
     &                    MXXF(if+m,jf+n,kf,if+m-1,jf+n-1,kf+kdir)
     &                 +sigmaf(if+m-1,jf+n  ,kf+ki-1,1) * 
     &                    MXXF(if+m,jf+n,kf,if+m-1,jf+n+1,kf+kdir)
     &                 +sigmaf(if+m  ,jf+n-1,kf+ki-1,1) * 
     &                    MXXF(if+m,jf+n,kf,if+m+1,jf+n-1,kf+kdir)
     &                 +sigmaf(if+m  ,jf+n  ,kf+ki-1,1) * 
     &                    MXXF(if+m,jf+n,kf,if+m+1,jf+n+1,kf+kdir)
                  tmp = tmp
     &                 +sigmaf(if+m-1,jf+n-1,kf+ki-1,2) * 
     &                    MYYF(if+m,jf+n,kf,if+m-1,jf+n-1,kf+kdir)
     &                 +sigmaf(if+m-1,jf+n  ,kf+ki-1,2) * 
     &                    MYYF(if+m,jf+n,kf,if+m-1,jf+n+1,kf+kdir)
     &                 +sigmaf(if+m  ,jf+n-1,kf+ki-1,2) * 
     &                    MYYF(if+m,jf+n,kf,if+m+1,jf+n-1,kf+kdir)
     &                 +sigmaf(if+m  ,jf+n  ,kf+ki-1,2) * 
     &                    MYYF(if+m,jf+n,kf,if+m+1,jf+n+1,kf+kdir)
                  tmp = tmp
     &                 +sigmaf(if+m-1,jf+n-1,kf+ki-1,3) * 
     &                    MZZF(if+m,jf+n,kf,if+m-1,jf+n-1,kf+kdir)
     &                 +sigmaf(if+m-1,jf+n  ,kf+ki-1,3) * 
     &                    MZZF(if+m,jf+n,kf,if+m-1,jf+n+1,kf+kdir)
     &                 +sigmaf(if+m  ,jf+n-1,kf+ki-1,3) * 
     &                    MZZF(if+m,jf+n,kf,if+m+1,jf+n-1,kf+kdir)
     &                 +sigmaf(if+m  ,jf+n  ,kf+ki-1,3) * 
     &                    MZZF(if+m,jf+n,kf,if+m+1,jf+n+1,kf+kdir)
                  tmp = tmp + kdir * (
     &                  sigmaf(if+m-1,jf+n-1,kf+ki-1,4) * 
     &                    MXZF(if+m,jf+n,kf,if+m-1,jf+n-1,kf+kdir)
     &                 +sigmaf(if+m-1,jf+n  ,kf+ki-1,4) * 
     &                    MXZF(if+m,jf+n,kf,if+m-1,jf+n+1,kf+kdir)
     &                 -sigmaf(if+m  ,jf+n-1,kf+ki-1,4) * 
     &                    MXZF(if+m,jf+n,kf,if+m+1,jf+n-1,kf+kdir)
     &                 -sigmaf(if+m  ,jf+n  ,kf+ki-1,4) * 
     &                    MXZF(if+m,jf+n,kf,if+m+1,jf+n+1,kf+kdir) )
                  tmp = tmp + kdir * (
     &                  sigmaf(if+m-1,jf+n-1,kf+ki-1,5) * 
     &                    MYZF(if+m,jf+n,kf,if+m-1,jf+n-1,kf+kdir)
     &                 -sigmaf(if+m-1,jf+n  ,kf+ki-1,5) * 
     &                    MYZF(if+m,jf+n,kf,if+m-1,jf+n+1,kf+kdir)
     &                 +sigmaf(if+m  ,jf+n-1,kf+ki-1,5) * 
     &                    MYZF(if+m,jf+n,kf,if+m+1,jf+n-1,kf+kdir)
     &                 -sigmaf(if+m  ,jf+n  ,kf+ki-1,5) * 
     &                    MYZF(if+m,jf+n,kf,if+m+1,jf+n+1,kf+kdir) )

         cen =
     &      4.d0 * (sigmaf(if+m-1,jf+n-1,kf+ki-1,1)
     &             +sigmaf(if+m-1,jf+n  ,kf+ki-1,1)
     &             +sigmaf(if+m  ,jf+n-1,kf+ki-1,1)
     &             +sigmaf(if+m  ,jf+n  ,kf+ki-1,1)
     &             +sigmaf(if+m-1,jf+n-1,kf+ki-1,2)
     &             +sigmaf(if+m-1,jf+n  ,kf+ki-1,2)
     &             +sigmaf(if+m  ,jf+n-1,kf+ki-1,2)
     &             +sigmaf(if+m  ,jf+n  ,kf+ki-1,2)
     &             +sigmaf(if+m-1,jf+n-1,kf+ki-1,3)
     &             +sigmaf(if+m-1,jf+n  ,kf+ki-1,3)
     &             +sigmaf(if+m  ,jf+n-1,kf+ki-1,3)
     &             +sigmaf(if+m  ,jf+n  ,kf+ki-1,3))
     &    + 6.d0 * (sigmaf(if+m-1,jf+n-1,kf+ki-1,4)
     &             +sigmaf(if+m-1,jf+n  ,kf+ki-1,4)
     &             -sigmaf(if+m  ,jf+n-1,kf+ki-1,4)
     &             -sigmaf(if+m  ,jf+n  ,kf+ki-1,4)) * kdir
     &    + 6.d0 * (sigmaf(if+m-1,jf+n-1,kf+ki-1,5)
     &             -sigmaf(if+m-1,jf+n  ,kf+ki-1,5)
     &             +sigmaf(if+m  ,jf+n-1,kf+ki-1,5)
     &             -sigmaf(if+m  ,jf+n  ,kf+ki-1,5)) * kdir

                  tmp = tmp - cen * fdst(if+m,jf+n,kf)
                  res(if,jf,kf) = res(if,jf,kf) + fac * tmp
                        end do
                     end do
                  end do
               end if
               if (ga(ii,0,ki) - ga(1-ii,0,ki) .eq. 1) then
                  fac0 = 1.d0 / (jr * kr)
                  do l = kdir, kdir*(kr-1), kdir
                     fac1 = (kr-abs(l)) * fac0
                     do n = -(jr-1), jr-1
                        fac = (jr-abs(n)) * fac1
                        do jf = jr*regl1, jr*regh1, jr

                  tmp = sigmaf(if+ii-1,jf+n-1,kf+l-1,1) * 
     &                    MXXF(if,jf+n,kf+l,if+idir,jf+n-1,kf+l-1)
     &                 +sigmaf(if+ii-1,jf+n-1,kf+l  ,1) * 
     &                    MXXF(if,jf+n,kf+l,if+idir,jf+n-1,kf+l+1)
     &                 +sigmaf(if+ii-1,jf+n  ,kf+l-1,1) * 
     &                    MXXF(if,jf+n,kf+l,if+idir,jf+n+1,kf+l-1)
     &                 +sigmaf(if+ii-1,jf+n  ,kf+l  ,1) * 
     &                    MXXF(if,jf+n,kf+l,if+idir,jf+n+1,kf+l+1)
                  tmp = tmp
     &                 +sigmaf(if+ii-1,jf+n-1,kf+l-1,2) * 
     &                    MYYF(if,jf+n,kf+l,if+idir,jf+n-1,kf+l-1)
     &                 +sigmaf(if+ii-1,jf+n-1,kf+l  ,2) * 
     &                    MYYF(if,jf+n,kf+l,if+idir,jf+n-1,kf+l+1)
     &                 +sigmaf(if+ii-1,jf+n  ,kf+l-1,2) * 
     &                    MYYF(if,jf+n,kf+l,if+idir,jf+n+1,kf+l-1)
     &                 +sigmaf(if+ii-1,jf+n  ,kf+l  ,2) * 
     &                    MYYF(if,jf+n,kf+l,if+idir,jf+n+1,kf+l+1)
                  tmp = tmp
     &                 +sigmaf(if+ii-1,jf+n-1,kf+l-1,3) * 
     &                    MZZF(if,jf+n,kf+l,if+idir,jf+n-1,kf+l-1)
     &                 +sigmaf(if+ii-1,jf+n-1,kf+l  ,3) * 
     &                    MZZF(if,jf+n,kf+l,if+idir,jf+n-1,kf+l+1)
     &                 +sigmaf(if+ii-1,jf+n  ,kf+l-1,3) * 
     &                    MZZF(if,jf+n,kf+l,if+idir,jf+n+1,kf+l-1)
     &                 +sigmaf(if+ii-1,jf+n  ,kf+l  ,3) * 
     &                    MZZF(if,jf+n,kf+l,if+idir,jf+n+1,kf+l+1)
                  tmp = tmp + idir * (
     &                 +sigmaf(if+ii-1,jf+n-1,kf+l-1,4) * 
     &                    MXZF(if,jf+n,kf+l,if+idir,jf+n-1,kf+l-1)
     &                 -sigmaf(if+ii-1,jf+n-1,kf+l  ,4) * 
     &                    MXZF(if,jf+n,kf+l,if+idir,jf+n-1,kf+l+1)
     &                 +sigmaf(if+ii-1,jf+n  ,kf+l-1,4) * 
     &                    MXZF(if,jf+n,kf+l,if+idir,jf+n+1,kf+l-1)
     &                 -sigmaf(if+ii-1,jf+n  ,kf+l  ,4) * 
     &                    MXZF(if,jf+n,kf+l,if+idir,jf+n+1,kf+l+1) )
                  tmp = tmp + (
     &                 -sigmaf(if+ii-1,jf+n-1,kf+l-1,5) * 
     &                    MYZF(if,jf+n,kf+l,if+idir,jf+n-1,kf+l-1)
     &                 +sigmaf(if+ii-1,jf+n-1,kf+l  ,5) * 
     &                    MYZF(if,jf+n,kf+l,if+idir,jf+n-1,kf+l+1)
     &                 +sigmaf(if+ii-1,jf+n  ,kf+l-1,5) * 
     &                    MYZF(if,jf+n,kf+l,if+idir,jf+n+1,kf+l-1)
     &                 -sigmaf(if+ii-1,jf+n  ,kf+l  ,5) * 
     &                    MYZF(if,jf+n,kf+l,if+idir,jf+n+1,kf+l+1) )

                  cen =
     &               4.d0 * (sigmaf(if+ii-1,jf+n-1,kf+l-1,1) 
     &                      +sigmaf(if+ii-1,jf+n-1,kf+l  ,1)
     &                      +sigmaf(if+ii-1,jf+n  ,kf+l-1,1) 
     &                      +sigmaf(if+ii-1,jf+n  ,kf+l  ,1)
     &                      +sigmaf(if+ii-1,jf+n-1,kf+l-1,2) 
     &                      +sigmaf(if+ii-1,jf+n-1,kf+l  ,2)
     &                      +sigmaf(if+ii-1,jf+n  ,kf+l-1,2) 
     &                      +sigmaf(if+ii-1,jf+n  ,kf+l  ,2)
     &                      +sigmaf(if+ii-1,jf+n-1,kf+l-1,3) 
     &                      +sigmaf(if+ii-1,jf+n-1,kf+l  ,3)
     &                      +sigmaf(if+ii-1,jf+n  ,kf+l-1,3) 
     &                      +sigmaf(if+ii-1,jf+n  ,kf+l  ,3))
     &             + 6.d0 * (sigmaf(if+ii-1,jf+n-1,kf+l-1,4) 
     &                      -sigmaf(if+ii-1,jf+n-1,kf+l  ,4)
     &                      +sigmaf(if+ii-1,jf+n  ,kf+l-1,4) 
     &                      -sigmaf(if+ii-1,jf+n  ,kf+l  ,4)) * idir
     &             + 6.d0 * (sigmaf(if+ii-1,jf+n-1,kf+l  ,5) 
     &                      -sigmaf(if+ii-1,jf+n-1,kf+l-1,5)
     &                      +sigmaf(if+ii-1,jf+n  ,kf+l-1,5) 
     &                      -sigmaf(if+ii-1,jf+n  ,kf+l  ,5))

                  tmp = tmp - cen * fdst(if,jf+n,kf+l)
                  res(if,jf,kf) = res(if,jf,kf) + fac * tmp
                        end do
                     end do
                  end do
               end if
            end do
         end do
c weighting
         do  jf = jr*regl1, jr*regh1, jr
            res(if,jf,kf) = src(if,jf,kf) - res(if,jf,kf) * fac00
         end do
      else
         do kf = kr*regl2, kr*regh2, kr
            res(if,jf,kf) = 0.d0
         end do
c quadrants
c each quadrant is two octants and their share of the two central edges
         fac0 = 1.d0 / kr
         do  ji = 0, 1
            jdir = 2 * ji - 1
            do ii = 0, 1
               idir = 2 * ii - 1
               if (ga(ii,ji,0) .eq. 1) then
                  do l = -(kr-1), kr-1
                     fac = (kr-abs(l)) * fac0
                     do kf = kr*regl2, kr*regh2, kr

                  tmp = sigmaf(if+ii-1,jf+ji-1,kf+l-1,1) * 
     &                    MXXF(if,jf,kf+l,if+idir,jf+jdir,kf+l-1)
     &                 +sigmaf(if+ii-1,jf+ji-1,kf+l  ,1) * 
     &                    MXXF(if,jf,kf+l,if+idir,jf+jdir,kf+l+1)
                  tmp = tmp
     &                 +sigmaf(if+ii-1,jf+ji-1,kf+l-1,2) * 
     &                    MYYF(if,jf,kf+l,if+idir,jf+jdir,kf+l-1)
     &                 +sigmaf(if+ii-1,jf+ji-1,kf+l  ,2) * 
     &                    MYYF(if,jf,kf+l,if+idir,jf+jdir,kf+l+1)
                  tmp = tmp
     &                 +sigmaf(if+ii-1,jf+ji-1,kf+l-1,3) * 
     &                    MZZF(if,jf,kf+l,if+idir,jf+jdir,kf+l-1)
     &                 +sigmaf(if+ii-1,jf+ji-1,kf+l  ,3) * 
     &                    MZZF(if,jf,kf+l,if+idir,jf+jdir,kf+l+1)
                  tmp = tmp + idir * (
     &                 +sigmaf(if+ii-1,jf+ji-1,kf+l-1,4) * 
     &                    MXZF(if,jf,kf+l,if+idir,jf+jdir,kf+l-1)
     &                 -sigmaf(if+ii-1,jf+ji-1,kf+l  ,4) * 
     &                    MXZF(if,jf,kf+l,if+idir,jf+jdir,kf+l+1) )
                  tmp = tmp + jdir * (
     &                 +sigmaf(if+ii-1,jf+ji-1,kf+l-1,5) * 
     &                    MYZF(if,jf,kf+l,if+idir,jf+jdir,kf+l-1)
     &                 -sigmaf(if+ii-1,jf+ji-1,kf+l  ,5) * 
     &                    MYZF(if,jf,kf+l,if+idir,jf+jdir,kf+l+1) )
                  cen =
     &               4.d0 * (sigmaf(if+ii-1,jf+ji-1,kf+l-1,1) 
     &                      +sigmaf(if+ii-1,jf+ji-1,kf+l  ,1)
     &                      +sigmaf(if+ii-1,jf+ji-1,kf+l-1,2) 
     &                      +sigmaf(if+ii-1,jf+ji-1,kf+l  ,2)
     &                      +sigmaf(if+ii-1,jf+ji-1,kf+l-1,3) 
     &                      +sigmaf(if+ii-1,jf+ji-1,kf+l  ,3))
     &             + 6.d0 * (sigmaf(if+ii-1,jf+ji-1,kf+l-1,4) 
     &                      -sigmaf(if+ii-1,jf+ji-1,kf+l  ,4)) * idir
     &             + 6.d0 * (sigmaf(if+ii-1,jf+ji-1,kf+l-1,5) 
     &                      -sigmaf(if+ii-1,jf+ji-1,kf+l  ,5)) * jdir

                  tmp = tmp - cen * fdst(if,jf,kf+l)

                        res(if,jf,kf) = res(if,jf,kf) + fac * tmp
                     end do
                  end do
               else
                  do kc = regl2, regh2
                     kf = kc * kr

                  tmp = sigmac(ic+ii-1,jc+ji-1,kc-1,1) * 
     &                    MXXC(ic+idir,jc+jdir,kc-1)
     &                 +sigmac(ic+ii-1,jc+ji-1,kc  ,1) * 
     &                    MXXC(ic+idir,jc+jdir,kc+1)
                  tmp = tmp
     &                 +sigmac(ic+ii-1,jc+ji-1,kc-1,2) * 
     &                    MYYC(ic+idir,jc+jdir,kc-1)
     &                 +sigmac(ic+ii-1,jc+ji-1,kc  ,2) * 
     &                    MYYC(ic+idir,jc+jdir,kc+1)
                  tmp = tmp
     &                 +sigmac(ic+ii-1,jc+ji-1,kc-1,3) * 
     &                    MZZC(ic+idir,jc+jdir,kc-1)
     &                 +sigmac(ic+ii-1,jc+ji-1,kc  ,3) * 
     &                    MZZC(ic+idir,jc+jdir,kc+1)
                  tmp = tmp + idir * (
     &                 +sigmac(ic+ii-1,jc+ji-1,kc-1,4) * 
     &                    MXZC(ic+idir,jc+jdir,kc-1)
     &                 -sigmac(ic+ii-1,jc+ji-1,kc  ,4) * 
     &                    MXZC(ic+idir,jc+jdir,kc+1) )
                  tmp = tmp + jdir * (
     &                 +sigmac(ic+ii-1,jc+ji-1,kc-1,5) * 
     &                    MYZC(ic+idir,jc+jdir,kc-1)
     &                 -sigmac(ic+ii-1,jc+ji-1,kc  ,5) * 
     &                    MYZC(ic+idir,jc+jdir,kc+1) )
                  cen =
     &               4.d0 * (sigmac(ic+ii-1,jc+ji-1,kc-1,1)
     &                      +sigmac(ic+ii-1,jc+ji-1,kc  ,1)
     &                      +sigmac(ic+ii-1,jc+ji-1,kc-1,2)
     &                      +sigmac(ic+ii-1,jc+ji-1,kc  ,2)
     &                      +sigmac(ic+ii-1,jc+ji-1,kc-1,3)
     &                      +sigmac(ic+ii-1,jc+ji-1,kc  ,3) )
     &             + 6.d0 * (sigmac(ic+ii-1,jc+ji-1,kc-1,4) 
     &                      -sigmac(ic+ii-1,jc+ji-1,kc  ,4) ) * idir
     &             + 6.d0 * (sigmac(ic+ii-1,jc+ji-1,kc-1,5)
     &                      -sigmac(ic+ii-1,jc+ji-1,kc  ,5) ) * jdir

                     tmp = tmp - cen * cdst(ic,jc,kc)
                     res(if,jf,kf) = res(if,jf,kf) + tmp
                  end do
               end if
            end do
         end do
c faces
c each face is two faces and two sides of an edge
         do ji = 0, 1
            jdir = 2 * ji - 1
            do ii = 0, 1
               idir = 2 * ii - 1
               if (ga(ii,ji,0) - ga(ii,1-ji,0) .eq. 1) then
                  fac0 = 1.d0 / (ir * kr)
                  do l = -(kr-1), kr-1
                     fac1 = (kr-abs(l)) * fac0
                     do m = idir, idir*(ir-1), idir
                        fac = (ir-abs(m)) * fac1
                        do kf = kr*regl2, kr*regh2, kr

                  tmp = sigmaf(if+m-1,jf+ji-1,kf+l-1,1) * 
     &                    MXXF(if+m,jf,kf+l,if+m-1,jf+jdir,kf+l-1)
     &                 +sigmaf(if+m-1,jf+ji-1,kf+l  ,1) * 
     &                    MXXF(if+m,jf,kf+l,if+m-1,jf+jdir,kf+l+1)
     &                 +sigmaf(if+m  ,jf+ji-1,kf+l-1,1) * 
     &                    MXXF(if+m,jf,kf+l,if+m+1,jf+jdir,kf+l-1)
     &                 +sigmaf(if+m  ,jf+ji-1,kf+l  ,1) * 
     &                    MXXF(if+m,jf,kf+l,if+m+1,jf+jdir,kf+l+1)
                  tmp = tmp
     &                 +sigmaf(if+m-1,jf+ji-1,kf+l-1,2) * 
     &                    MYYF(if+m,jf,kf+l,if+m-1,jf+jdir,kf+l-1)
     &                 +sigmaf(if+m-1,jf+ji-1,kf+l  ,2) * 
     &                    MYYF(if+m,jf,kf+l,if+m-1,jf+jdir,kf+l+1)
     &                 +sigmaf(if+m  ,jf+ji-1,kf+l-1,2) * 
     &                    MYYF(if+m,jf,kf+l,if+m+1,jf+jdir,kf+l-1)
     &                 +sigmaf(if+m  ,jf+ji-1,kf+l  ,2) * 
     &                    MYYF(if+m,jf,kf+l,if+m+1,jf+jdir,kf+l+1)
                  tmp = tmp
     &                 +sigmaf(if+m-1,jf+ji-1,kf+l-1,3) * 
     &                    MZZF(if+m,jf,kf+l,if+m-1,jf+jdir,kf+l-1)
     &                 +sigmaf(if+m-1,jf+ji-1,kf+l  ,3) * 
     &                    MZZF(if+m,jf,kf+l,if+m-1,jf+jdir,kf+l+1)
     &                 +sigmaf(if+m  ,jf+ji-1,kf+l-1,3) * 
     &                    MZZF(if+m,jf,kf+l,if+m+1,jf+jdir,kf+l-1)
     &                 +sigmaf(if+m  ,jf+ji-1,kf+l  ,3) * 
     &                    MZZF(if+m,jf,kf+l,if+m+1,jf+jdir,kf+l+1)
                  tmp = tmp + (
     &                 -sigmaf(if+m-1,jf+ji-1,kf+l-1,4) * 
     &                    MXZF(if+m,jf,kf+l,if+m-1,jf+jdir,kf+l-1)
     &                 +sigmaf(if+m-1,jf+ji-1,kf+l  ,4) * 
     &                    MXZF(if+m,jf,kf+l,if+m-1,jf+jdir,kf+l+1)
     &                 +sigmaf(if+m  ,jf+ji-1,kf+l-1,4) * 
     &                    MXZF(if+m,jf,kf+l,if+m+1,jf+jdir,kf+l-1)
     &                 -sigmaf(if+m  ,jf+ji-1,kf+l  ,4) * 
     &                    MXZF(if+m,jf,kf+l,if+m+1,jf+jdir,kf+l+1) )
                  tmp = tmp + jdir * (
     &                 +sigmaf(if+m-1,jf+ji-1,kf+l-1,5) * 
     &                    MYZF(if+m,jf,kf+l,if+m-1,jf+jdir,kf+l-1)
     &                 -sigmaf(if+m-1,jf+ji-1,kf+l  ,5) * 
     &                    MYZF(if+m,jf,kf+l,if+m-1,jf+jdir,kf+l+1)
     &                 +sigmaf(if+m  ,jf+ji-1,kf+l-1,5) * 
     &                    MYZF(if+m,jf,kf+l,if+m+1,jf+jdir,kf+l-1)
     &                 -sigmaf(if+m  ,jf+ji-1,kf+l  ,5) * 
     &                    MYZF(if+m,jf,kf+l,if+m+1,jf+jdir,kf+l+1) )
                  cen =
     &               4.d0 * (sigmaf(if+m-1,jf+ji-1,kf+l-1,1)  
     &                      +sigmaf(if+m-1,jf+ji-1,kf+l  ,1)
     &                      +sigmaf(if+m  ,jf+ji-1,kf+l-1,1)  
     &                      +sigmaf(if+m  ,jf+ji-1,kf+l  ,1)
     &                      +sigmaf(if+m-1,jf+ji-1,kf+l-1,2)  
     &                      +sigmaf(if+m-1,jf+ji-1,kf+l  ,2)
     &                      +sigmaf(if+m  ,jf+ji-1,kf+l-1,2)  
     &                      +sigmaf(if+m  ,jf+ji-1,kf+l  ,2)
     &                      +sigmaf(if+m-1,jf+ji-1,kf+l-1,3)  
     &                      +sigmaf(if+m-1,jf+ji-1,kf+l  ,3)
     &                      +sigmaf(if+m  ,jf+ji-1,kf+l-1,3)  
     &                      +sigmaf(if+m  ,jf+ji-1,kf+l  ,3))
     &             + 6.d0 * (sigmaf(if+m-1,jf+ji-1,kf+l  ,4)  
     &                      -sigmaf(if+m-1,jf+ji-1,kf+l-1,4)
     &                      +sigmaf(if+m  ,jf+ji-1,kf+l-1,4) 
     &                      -sigmaf(if+m  ,jf+ji-1,kf+l  ,4))
     &             + 6.d0 * (sigmaf(if+m-1,jf+ji-1,kf+l-1,5)  
     &                      -sigmaf(if+m-1,jf+ji-1,kf+l  ,5)
     &                      +sigmaf(if+m  ,jf+ji-1,kf+l-1,5)  
     &                      -sigmaf(if+m  ,jf+ji-1,kf+l  ,5)) * jdir
                  tmp = tmp - cen * fdst(if+m,jf,kf+l)
                  res(if,jf,kf) = res(if,jf,kf) + fac * tmp
                        end do
                     end do
                  end do
               end if
               if (ga(ii,ji,0) - ga(1-ii,ji,0) .eq. 1) then
                  fac0 = 1.d0 / (jr * kr)
                  do l = -(kr-1), kr-1
                     fac1 = (kr-abs(l)) * fac0
                     do n = jdir, jdir*(jr-1), jdir
                        fac = (jr-abs(n)) * fac1
                        do kf = kr*regl2, kr*regh2, kr
                  tmp = sigmaf(if+ii-1,jf+n-1,kf+l-1,1) * 
     &                    MXXF(if,jf+n,kf+l,if+idir,jf+n-1,kf+l-1)
     &                 +sigmaf(if+ii-1,jf+n-1,kf+l  ,1) * 
     &                    MXXF(if,jf+n,kf+l,if+idir,jf+n-1,kf+l+1)
     &                 +sigmaf(if+ii-1,jf+n  ,kf+l-1,1) * 
     &                    MXXF(if,jf+n,kf+l,if+idir,jf+n+1,kf+l-1)
     &                 +sigmaf(if+ii-1,jf+n  ,kf+l  ,1) * 
     &                    MXXF(if,jf+n,kf+l,if+idir,jf+n+1,kf+l+1)
                  tmp = tmp
     &                 +sigmaf(if+ii-1,jf+n-1,kf+l-1,2) * 
     &                    MYYF(if,jf+n,kf+l,if+idir,jf+n-1,kf+l-1)
     &                 +sigmaf(if+ii-1,jf+n-1,kf+l  ,2) * 
     &                    MYYF(if,jf+n,kf+l,if+idir,jf+n-1,kf+l+1)
     &                 +sigmaf(if+ii-1,jf+n  ,kf+l-1,2) * 
     &                    MYYF(if,jf+n,kf+l,if+idir,jf+n+1,kf+l-1)
     &                 +sigmaf(if+ii-1,jf+n  ,kf+l  ,2) * 
     &                    MYYF(if,jf+n,kf+l,if+idir,jf+n+1,kf+l+1)
                  tmp = tmp
     &                 +sigmaf(if+ii-1,jf+n-1,kf+l-1,3) * 
     &                    MZZF(if,jf+n,kf+l,if+idir,jf+n-1,kf+l-1)
     &                 +sigmaf(if+ii-1,jf+n-1,kf+l  ,3) * 
     &                    MZZF(if,jf+n,kf+l,if+idir,jf+n-1,kf+l+1)
     &                 +sigmaf(if+ii-1,jf+n  ,kf+l-1,3) * 
     &                    MZZF(if,jf+n,kf+l,if+idir,jf+n+1,kf+l-1)
     &                 +sigmaf(if+ii-1,jf+n  ,kf+l  ,3) * 
     &                    MZZF(if,jf+n,kf+l,if+idir,jf+n+1,kf+l+1)
                  tmp = tmp + idir * (
     &                 +sigmaf(if+ii-1,jf+n-1,kf+l-1,4) * 
     &                    MXZF(if,jf+n,kf+l,if+idir,jf+n-1,kf+l-1)
     &                 -sigmaf(if+ii-1,jf+n-1,kf+l  ,4) * 
     &                    MXZF(if,jf+n,kf+l,if+idir,jf+n-1,kf+l+1)
     &                 +sigmaf(if+ii-1,jf+n  ,kf+l-1,4) * 
     &                    MXZF(if,jf+n,kf+l,if+idir,jf+n+1,kf+l-1)
     &                 -sigmaf(if+ii-1,jf+n  ,kf+l  ,4) * 
     &                    MXZF(if,jf+n,kf+l,if+idir,jf+n+1,kf+l+1) )
                  tmp = tmp + (
     &                 -sigmaf(if+ii-1,jf+n-1,kf+l-1,5) * 
     &                    MYZF(if,jf+n,kf+l,if+idir,jf+n-1,kf+l-1)
     &                 +sigmaf(if+ii-1,jf+n-1,kf+l  ,5) * 
     &                    MYZF(if,jf+n,kf+l,if+idir,jf+n-1,kf+l+1)
     &                 +sigmaf(if+ii-1,jf+n  ,kf+l-1,5) * 
     &                    MYZF(if,jf+n,kf+l,if+idir,jf+n+1,kf+l-1)
     &                 -sigmaf(if+ii-1,jf+n  ,kf+l  ,5) * 
     &                    MYZF(if,jf+n,kf+l,if+idir,jf+n+1,kf+l+1) )

                  cen =
     &               4.d0 * (sigmaf(if+ii-1,jf+n-1,kf+l-1,1) 
     &                      +sigmaf(if+ii-1,jf+n-1,kf+l  ,1)
     &                      +sigmaf(if+ii-1,jf+n  ,kf+l-1,1) 
     &                      +sigmaf(if+ii-1,jf+n  ,kf+l  ,1)
     &                      +sigmaf(if+ii-1,jf+n-1,kf+l-1,2) 
     &                      +sigmaf(if+ii-1,jf+n-1,kf+l  ,2)
     &                      +sigmaf(if+ii-1,jf+n  ,kf+l-1,2) 
     &                      +sigmaf(if+ii-1,jf+n  ,kf+l  ,2)
     &                      +sigmaf(if+ii-1,jf+n-1,kf+l-1,3) 
     &                      +sigmaf(if+ii-1,jf+n-1,kf+l  ,3)
     &                      +sigmaf(if+ii-1,jf+n  ,kf+l-1,3) 
     &                      +sigmaf(if+ii-1,jf+n  ,kf+l  ,3))
     &             + 6.d0 * (sigmaf(if+ii-1,jf+n-1,kf+l-1,4) 
     &                      -sigmaf(if+ii-1,jf+n-1,kf+l  ,4)
     &                      +sigmaf(if+ii-1,jf+n  ,kf+l-1,4) 
     &                      -sigmaf(if+ii-1,jf+n  ,kf+l  ,4)) * idir
     &             + 6.d0 * (sigmaf(if+ii-1,jf+n-1,kf+l  ,5) 
     &                      -sigmaf(if+ii-1,jf+n-1,kf+l-1,5)
     &                      +sigmaf(if+ii-1,jf+n  ,kf+l-1,5) 
     &                      -sigmaf(if+ii-1,jf+n  ,kf+l  ,5))
                  tmp = tmp - cen * fdst(if,jf+n,kf+l)
                  res(if,jf,kf) = res(if,jf,kf) + fac * tmp
                        end do
                     end do
                  end do
               end if
            end do
         end do
c weighting
         do kf = kr*regl2, kr*regh2, kr
            res(if,jf,kf) = src(if,jf,kf) - res(if,jf,kf) * fac00
         end do 
      end if
      end
c-----------------------------------------------------------------------
c Note---assumes fdst linearly interpolated from cdst along face
      subroutine hgcres_terrain(
     & res,    resl0,resh0,resl1,resh1,resl2,resh2,
     & src,    srcl0,srch0,srcl1,srch1,srcl2,srch2,
     & fdst,   fdstl0,fdsth0,fdstl1,fdsth1,fdstl2,fdsth2,
     & cdst,   cdstl0,cdsth0,cdstl1,cdsth1,cdstl2,cdsth2,
     & sigmaf, sfl0,sfh0,sfl1,sfh1,sfl2,sfh2,
     & sigmac, scl0,sch0,scl1,sch1,scl2,sch2,
     &         regl0,regh0,regl1,regh1,regl2,regh2,
     & hx, hy, hz, ir, jr, kr, ga, idd)
      integer resl0,resh0,resl1,resh1,resl2,resh2
      integer srcl0,srch0,srcl1,srch1,srcl2,srch2
      integer fdstl0,fdsth0,fdstl1,fdsth1,fdstl2,fdsth2
      integer cdstl0,cdsth0,cdstl1,cdsth1,cdstl2,cdsth2
      integer sfl0,sfh0,sfl1,sfh1,sfl2,sfh2
      integer scl0,sch0,scl1,sch1,scl2,sch2
      integer regl0,regh0,regl1,regh1,regl2,regh2
      double precision hx, hy, hz
      double precision res(resl0:resh0,resl1:resh1,resl2:resh2)
      double precision src(srcl0:srch0,srcl1:srch1,srcl2:srch2)
      double precision fdst(fdstl0:fdsth0,fdstl1:fdsth1,fdstl2:fdsth2)
      double precision cdst(cdstl0:cdsth0,cdstl1:cdsth1,cdstl2:cdsth2)
      double precision sigmaf(sfl0:sfh0,sfl1:sfh1,sfl2:sfh2,5)
      double precision sigmac(scl0:sch0,scl1:sch1,scl2:sch2,5)
      integer ir, jr, kr, ga(0:1,0:1,0:1), idd
      double precision cen, sum, fac00, fac1, fac2, fac
      integer ic, jc, kc, if, jf, kf, ii, ji, ki, idir, jdir, kdir
      integer l, m, n
      double precision MXXC, MYYC, MZZC, MXZC, MYZC
      double precision MXXF, MYYF, MZZF, MXZF, MYZF
      MXXC(ii,jj,kk) = (4.d0 *  cdst(ii,jc,kc) + 
     &                   2.d0 * (cdst(ii,jj,kc) - cdst(ic,jj,kc) + 
     &                           cdst(ii,jc,kk) - cdst(ic,jc,kk)) +
     &                          (cdst(ii,jj,kk) - cdst(ic,jj,kk)))

      MYYC(ii,jj,kk) = (4.d0 *  cdst(ic,jj,kc) + 
     &                   2.d0 * (cdst(ii,jj,kc) - cdst(ii,jc,kc) + 
     &                           cdst(ic,jj,kk) - cdst(ic,jc,kk)) +
     &                          (cdst(ii,jj,kk) - cdst(ii,jc,kk)))

      MZZC(ii,jj,kk) = (4.d0 *  cdst(ic,jc,kk) + 
     &                   2.d0 * (cdst(ii,jc,kk) - cdst(ii,jc,kc) + 
     &                           cdst(ic,jj,kk) - cdst(ic,jj,kc)) +
     &                          (cdst(ii,jj,kk) - cdst(ii,jj,kc)))

      MXZC(ii,jj,kk) = (6.d0 *  cdst(ii,jc,kk) +
     &                   3.d0 * (cdst(ii,jj,kk) - cdst(ic,jj,kc)))

      MYZC(ii,jj,kk) = (6.d0 *  cdst(ic,jj,kk) +
     &                   3.d0 * (cdst(ii,jj,kk) - cdst(ii,jc,kc)))

      MXXF(i1,j1,k1,ii,jj,kk) = (4.d0 *  fdst(ii,j1,k1)   + 
     &                         2.d0 * (fdst(ii,jj,k1) - fdst(i1,jj,k1) +
     &                               fdst(ii,j1,kk) - fdst(i1,j1,kk)) +
     &                                (fdst(ii,jj,kk) - fdst(i1,jj,kk)))

      MYYF(i1,j1,k1,ii,jj,kk) = (4.d0 *  fdst(i1,jj,k1)   + 
     &                         2.d0 * (fdst(ii,jj,k1) - fdst(ii,j1,k1) +
     &                                fdst(i1,jj,kk) - fdst(i1,j1,kk)) +
     &                                (fdst(ii,jj,kk) - fdst(ii,j1,kk)))

      MZZF(i1,j1,k1,ii,jj,kk) = (4.d0 *  fdst(i1,j1,kk)   + 
     &                         2.d0 * (fdst(ii,j1,kk) - fdst(ii,j1,k1) +
     &                                fdst(i1,jj,kk) - fdst(i1,jj,k1)) +
     &                                (fdst(ii,jj,kk) - fdst(ii,jj,k1)))

      MXZF(i1,j1,k1,ii,jj,kk) = (6.d0 *  fdst(ii,j1,kk)  +
     &                         3.d0 * (fdst(ii,jj,kk) - fdst(i1,jj,k1)))

      MYZF(i1,j1,k1,ii,jj,kk) = (6.d0 *  fdst(i1,jj,kk)  +
     &                         3.d0 * (fdst(ii,jj,kk) - fdst(ii,j1,k1)))
      ic = regl0
      jc = regl1
      kc = regl2
      if = ic * ir
      jf = jc * jr
      kf = kc * kr
      sum = 0.d0

      fac00 = 1.d0 / 36.d0

c octants
      do ki = 0, 1
         kdir = 2 * ki - 1
         do ji = 0, 1
            jdir = 2 * ji - 1
            do ii = 0, 1
               idir = 2 * ii - 1
               if (ga(ii,ji,ki) .eq. 1) then

                  sum = sum 
     &             +sigmaf(if+ii-1,jf+ji-1,kf+ki-1,1) * 
     &                MXXF(if,jf,kf,if+idir,jf+jdir,kf+kdir)
     &             +sigmaf(if+ii-1,jf+ji-1,kf+ki-1,2) * 
     &                MYYF(if,jf,kf,if+idir,jf+jdir,kf+kdir)
     &             +sigmaf(if+ii-1,jf+ji-1,kf+ki-1,3) * 
     &                MZZF(if,jf,kf,if+idir,jf+jdir,kf+kdir)
     &             -sigmaf(if+ii-1,jf+ji-1,kf+ki-1,4) * 
     &              MXZF(if,jf,kf,if+idir,jf+jdir,kf+kdir) * idir * kdir
     &             -sigmaf(if+ii-1,jf+ji-1,kf+ki-1,5) * 
     &              MYZF(if,jf,kf,if+idir,jf+jdir,kf+kdir) * jdir * kdir

         cen =
     &      4.d0 * ( sigmaf(if+ii-1,jf+ji-1,kf+ki-1,1)
     &              +sigmaf(if+ii-1,jf+ji-1,kf+ki-1,2)
     &              +sigmaf(if+ii-1,jf+ji-1,kf+ki-1,3))
     &    + 6.d0 * (-sigmaf(if+ii-1,jf+ji-1,kf+ki-1,4)) * idir * kdir
     &    + 6.d0 * (-sigmaf(if+ii-1,jf+ji-1,kf+ki-1,5)) * jdir * kdir

                   sum = sum - cen * fdst(if,jf,kf)

               else

                  sum = sum 
     &             +sigmac(ic+ii-1,jc+ji-1,kc+ki-1,1) * 
     &                MXXC(ic+idir,jc+jdir,kc+kdir)
     &             +sigmac(ic+ii-1,jc+ji-1,kc+ki-1,2) * 
     &                MYYC(ic+idir,jc+jdir,kc+kdir)
     &             +sigmac(ic+ii-1,jc+ji-1,kc+ki-1,3) * 
     &                MZZC(ic+idir,jc+jdir,kc+kdir)
     &             -sigmac(ic+ii-1,jc+ji-1,kc+ki-1,4) * 
     &                MXZC(ic+idir,jc+jdir,kc+kdir) * idir * kdir
     &             -sigmac(ic+ii-1,jc+ji-1,kc+ki-1,5) * 
     &                MYZC(ic+idir,jc+jdir,kc+kdir) * jdir * kdir

         cen =
     &      4.d0 * ( sigmac(ic+ii-1,jc+ji-1,kc+ki-1,1)
     &              +sigmac(ic+ii-1,jc+ji-1,kc+ki-1,2)
     &              +sigmac(ic+ii-1,jc+ji-1,kc+ki-1,3))
     &    + 6.d0 * (-sigmac(ic+ii-1,jc+ji-1,kc+ki-1,4)) * idir * kdir
     &    + 6.d0 * (-sigmac(ic+ii-1,jc+ji-1,kc+ki-1,5)) * jdir * kdir

                   sum = sum - cen * cdst(ic,jc,kc)

               end if
            end do
         end do
      end do
c faces
      do ki = 0, 1
         kdir = 2 * ki - 1
         do ji = 0, 1
            jdir = 2 * ji - 1
            do ii = 0, 1
               idir = 2 * ii - 1
               if (ga(ii,ji,ki) - ga(ii,ji,1-ki) .eq. 1) then
                  fac2 = 1.d0 / (ir * jr)
                  do n = jdir, jdir*(jr-1), jdir
                     fac1 = (jr-abs(n)) * fac2
                     do m = idir, idir*(ir-1), idir
                        fac = (ir-abs(m)) * fac1

         sum = sum + fac * (
     &     +sigmaf(if+m-1,jf+n-1,kf+ki-1,1) * 
     &        MXXF(if+m,jf+n,kf,if+m-1,jf+n-1,kf+kdir)
     &     +sigmaf(if+m-1,jf+n  ,kf+ki-1,1) * 
     &        MXXF(if+m,jf+n,kf,if+m-1,jf+n+1,kf+kdir)
     &     +sigmaf(if+m  ,jf+n-1,kf+ki-1,1) * 
     &        MXXF(if+m,jf+n,kf,if+m+1,jf+n-1,kf+kdir)
     &     +sigmaf(if+m  ,jf+n  ,kf+ki-1,1) * 
     &        MXXF(if+m,jf+n,kf,if+m+1,jf+n+1,kf+kdir))
         sum = sum + fac * (
     &     +sigmaf(if+m-1,jf+n-1,kf+ki-1,2) * 
     &        MYYF(if+m,jf+n,kf,if+m-1,jf+n-1,kf+kdir)
     &     +sigmaf(if+m-1,jf+n  ,kf+ki-1,2) * 
     &        MYYF(if+m,jf+n,kf,if+m-1,jf+n+1,kf+kdir)
     &     +sigmaf(if+m  ,jf+n-1,kf+ki-1,2) * 
     &        MYYF(if+m,jf+n,kf,if+m+1,jf+n-1,kf+kdir)
     &     +sigmaf(if+m  ,jf+n  ,kf+ki-1,2) * 
     &        MYYF(if+m,jf+n,kf,if+m+1,jf+n+1,kf+kdir))
         sum = sum + fac * (
     &     +sigmaf(if+m-1,jf+n-1,kf+ki-1,3) * 
     &        MZZF(if+m,jf+n,kf,if+m-1,jf+n-1,kf+kdir)
     &     +sigmaf(if+m-1,jf+n  ,kf+ki-1,3) * 
     &        MZZF(if+m,jf+n,kf,if+m-1,jf+n+1,kf+kdir)
     &     +sigmaf(if+m  ,jf+n-1,kf+ki-1,3) * 
     &        MZZF(if+m,jf+n,kf,if+m+1,jf+n-1,kf+kdir)
     &     +sigmaf(if+m  ,jf+n  ,kf+ki-1,3) * 
     &        MZZF(if+m,jf+n,kf,if+m+1,jf+n+1,kf+kdir))
         sum = sum + fac * kdir * (
     &     +sigmaf(if+m-1,jf+n-1,kf+ki-1,4) * 
     &        MXZF(if+m,jf+n,kf,if+m-1,jf+n-1,kf+kdir)
     &     +sigmaf(if+m-1,jf+n  ,kf+ki-1,4) * 
     &        MXZF(if+m,jf+n,kf,if+m-1,jf+n+1,kf+kdir)
     &     -sigmaf(if+m  ,jf+n-1,kf+ki-1,4) * 
     &        MXZF(if+m,jf+n,kf,if+m+1,jf+n-1,kf+kdir)
     &     -sigmaf(if+m  ,jf+n  ,kf+ki-1,4) * 
     &        MXZF(if+m,jf+n,kf,if+m+1,jf+n+1,kf+kdir))
         sum = sum + fac * kdir * (
     &     +sigmaf(if+m-1,jf+n-1,kf+ki-1,5) * 
     &        MYZF(if+m,jf+n,kf,if+m-1,jf+n-1,kf+kdir)
     &     -sigmaf(if+m-1,jf+n  ,kf+ki-1,5) * 
     &        MYZF(if+m,jf+n,kf,if+m-1,jf+n+1,kf+kdir)
     &     +sigmaf(if+m  ,jf+n-1,kf+ki-1,5) * 
     &        MYZF(if+m,jf+n,kf,if+m+1,jf+n-1,kf+kdir)
     &     -sigmaf(if+m  ,jf+n  ,kf+ki-1,5) * 
     &        MYZF(if+m,jf+n,kf,if+m+1,jf+n+1,kf+kdir))

         cen =
     &      4.d0 * ( sigmaf(if+m-1,jf+n-1,kf+ki-1,1)
     &              +sigmaf(if+m-1,jf+n  ,kf+ki-1,1)
     &              +sigmaf(if+m  ,jf+n-1,kf+ki-1,1)
     &              +sigmaf(if+m  ,jf+n  ,kf+ki-1,1)
     &              +sigmaf(if+m-1,jf+n-1,kf+ki-1,2)
     &              +sigmaf(if+m-1,jf+n  ,kf+ki-1,2)
     &              +sigmaf(if+m  ,jf+n-1,kf+ki-1,2)
     &              +sigmaf(if+m  ,jf+n  ,kf+ki-1,2)
     &              +sigmaf(if+m-1,jf+n-1,kf+ki-1,3)
     &              +sigmaf(if+m-1,jf+n  ,kf+ki-1,3)
     &              +sigmaf(if+m  ,jf+n-1,kf+ki-1,3)
     &              +sigmaf(if+m  ,jf+n  ,kf+ki-1,3))
     &    + 6.d0 * ( sigmaf(if+m-1,jf+n-1,kf+ki-1,4)
     &              +sigmaf(if+m-1,jf+n  ,kf+ki-1,4)
     &              -sigmaf(if+m  ,jf+n-1,kf+ki-1,4)
     &              -sigmaf(if+m  ,jf+n  ,kf+ki-1,4)) * kdir
     &    + 6.d0 * ( sigmaf(if+m-1,jf+n-1,kf+ki-1,5)
     &              -sigmaf(if+m-1,jf+n  ,kf+ki-1,5)
     &              +sigmaf(if+m  ,jf+n-1,kf+ki-1,5)
     &              -sigmaf(if+m  ,jf+n  ,kf+ki-1,5)) * kdir

          sum = sum - fac * cen * fdst(if+m,jf+n,kf)
                     end do
                  end do
               end if

               if (ga(ii,ji,ki) - ga(ii,1-ji,ki) .eq. 1) then
                  fac2 = 1.d0 / (ir * kr)
                  do l = kdir, kdir*(kr-1), kdir
                     fac1 = (kr-abs(l)) * fac2
                     do m = idir, idir*(ir-1), idir
                        fac = (ir-abs(m)) * fac1

         sum = sum + fac * (
     &                 +sigmaf(if+m-1,jf+ji-1,kf+l-1,1) * 
     &                    MXXF(if+m,jf,kf+l,if+m-1,jf+jdir,kf+l-1)
     &                 +sigmaf(if+m-1,jf+ji-1,kf+l  ,1) * 
     &                    MXXF(if+m,jf,kf+l,if+m-1,jf+jdir,kf+l+1)
     &                 +sigmaf(if+m  ,jf+ji-1,kf+l-1,1) * 
     &                    MXXF(if+m,jf,kf+l,if+m+1,jf+jdir,kf+l-1)
     &                 +sigmaf(if+m  ,jf+ji-1,kf+l  ,1) * 
     &                    MXXF(if+m,jf,kf+l,if+m+1,jf+jdir,kf+l+1) )
         sum = sum + fac * (
     &                 +sigmaf(if+m-1,jf+ji-1,kf+l-1,2) * 
     &                    MYYF(if+m,jf,kf+l,if+m-1,jf+jdir,kf+l-1)
     &                 +sigmaf(if+m-1,jf+ji-1,kf+l  ,2) * 
     &                    MYYF(if+m,jf,kf+l,if+m-1,jf+jdir,kf+l+1)
     &                 +sigmaf(if+m  ,jf+ji-1,kf+l-1,2) * 
     &                    MYYF(if+m,jf,kf+l,if+m+1,jf+jdir,kf+l-1)
     &                 +sigmaf(if+m  ,jf+ji-1,kf+l  ,2) * 
     &                    MYYF(if+m,jf,kf+l,if+m+1,jf+jdir,kf+l+1) )
         sum = sum + fac * (
     &                 +sigmaf(if+m-1,jf+ji-1,kf+l-1,3) * 
     &                    MZZF(if+m,jf,kf+l,if+m-1,jf+jdir,kf+l-1)
     &                 +sigmaf(if+m-1,jf+ji-1,kf+l  ,3) * 
     &                    MZZF(if+m,jf,kf+l,if+m-1,jf+jdir,kf+l+1)
     &                 +sigmaf(if+m  ,jf+ji-1,kf+l-1,3) * 
     &                    MZZF(if+m,jf,kf+l,if+m+1,jf+jdir,kf+l-1)
     &                 +sigmaf(if+m  ,jf+ji-1,kf+l  ,3) * 
     &                    MZZF(if+m,jf,kf+l,if+m+1,jf+jdir,kf+l+1) )
         sum = sum + fac * (
     &                 -sigmaf(if+m-1,jf+ji-1,kf+l-1,4) * 
     &                    MXZF(if+m,jf,kf+l,if+m-1,jf+jdir,kf+l-1)
     &                 +sigmaf(if+m-1,jf+ji-1,kf+l  ,4) * 
     &                    MXZF(if+m,jf,kf+l,if+m-1,jf+jdir,kf+l+1)
     &                 +sigmaf(if+m  ,jf+ji-1,kf+l-1,4) * 
     &                    MXZF(if+m,jf,kf+l,if+m+1,jf+jdir,kf+l-1)
     &                 -sigmaf(if+m  ,jf+ji-1,kf+l  ,4) * 
     &                    MXZF(if+m,jf,kf+l,if+m+1,jf+jdir,kf+l+1) )
         sum = sum + fac * jdir * (
     &                 +sigmaf(if+m-1,jf+ji-1,kf+l-1,5) * 
     &                    MYZF(if+m,jf,kf+l,if+m-1,jf+jdir,kf+l-1)
     &                 -sigmaf(if+m-1,jf+ji-1,kf+l  ,5) * 
     &                    MYZF(if+m,jf,kf+l,if+m-1,jf+jdir,kf+l+1)
     &                 +sigmaf(if+m  ,jf+ji-1,kf+l-1,5) * 
     &                    MYZF(if+m,jf,kf+l,if+m+1,jf+jdir,kf+l-1)
     &                 -sigmaf(if+m  ,jf+ji-1,kf+l  ,5) * 
     &                    MYZF(if+m,jf,kf+l,if+m+1,jf+jdir,kf+l+1) )

         cen =
     &      4.d0 * (
     &            + sigmaf(if+m-1,jf+ji-1,kf+l-1,1) 
     &            + sigmaf(if+m-1,jf+ji-1,kf+l,1)
     &            + sigmaf(if+m  ,jf+ji-1,kf+l-1,1) 
     &            + sigmaf(if+m  ,jf+ji-1,kf+l,1)
     &        +sigmaf(if+m-1,jf+ji-1,kf+l-1,2) 
     &            + sigmaf(if+m-1,jf+ji-1,kf+l,2)
     &        +sigmaf(if+m  ,jf+ji-1,kf+l-1,2) 
     &            + sigmaf(if+m  ,jf+ji-1,kf+l,2)
     &        +sigmaf(if+m-1,jf+ji-1,kf+l-1,3) 
     &            + sigmaf(if+m-1,jf+ji-1,kf+l,3)
     &        +sigmaf(if+m  ,jf+ji-1,kf+l-1,3) 
     &            + sigmaf(if+m  ,jf+ji-1,kf+l,3))
     &     +6.d0 * (
     &         sigmaf(if+m-1,jf+ji-1,kf+l  ,4) 
     &            - sigmaf(if+m-1,jf+ji-1,kf+l-1,4)
     &        +sigmaf(if+m  ,jf+ji-1,kf+l-1,4) 
     &            - sigmaf(if+m  ,jf+ji-1,kf+l,4))
     &     +6.d0 * jdir * (
     &         sigmaf(if+m-1,jf+ji-1,kf+l-1,5) 
     &            - sigmaf(if+m-1,jf+ji-1,kf+l,5)
     &        +sigmaf(if+m  ,jf+ji-1,kf+l-1,5) 
     &            - sigmaf(if+m  ,jf+ji-1,kf+l,5))

          sum = sum - fac * cen * fdst(if+m,jf,kf+l)
                     end do
                  end do
               end if
               if (ga(ii,ji,ki) - ga(1-ii,ji,ki) .eq. 1) then
                  fac2 = 1.d0 / (jr * kr)
                  do l = kdir, kdir*(kr-1), kdir
                     fac1 = (kr-abs(l)) * fac2
                     do n = jdir, jdir*(jr-1), jdir
                        fac = (jr-abs(n)) * fac1

         sum = sum + fac * (
     &                 +sigmaf(if+ii-1,jf+n-1,kf+l-1,1) * 
     &                    MXXF(if,jf+n,kf+l,if+idir,jf+n-1,kf+l-1)
     &                 +sigmaf(if+ii-1,jf+n-1,kf+l  ,1) * 
     &                    MXXF(if,jf+n,kf+l,if+idir,jf+n-1,kf+l+1)
     &                 +sigmaf(if+ii-1,jf+n  ,kf+l-1,1) * 
     &                    MXXF(if,jf+n,kf+l,if+idir,jf+n+1,kf+l-1)
     &                 +sigmaf(if+ii-1,jf+n  ,kf+l  ,1) * 
     &                    MXXF(if,jf+n,kf+l,if+idir,jf+n+1,kf+l+1) )
         sum = sum + fac * (
     &                 +sigmaf(if+ii-1,jf+n-1,kf+l-1,2) * 
     &                    MYYF(if,jf+n,kf+l,if+idir,jf+n-1,kf+l-1)
     &                 +sigmaf(if+ii-1,jf+n-1,kf+l  ,2) * 
     &                    MYYF(if,jf+n,kf+l,if+idir,jf+n-1,kf+l+1)
     &                 +sigmaf(if+ii-1,jf+n  ,kf+l-1,2) * 
     &                    MYYF(if,jf+n,kf+l,if+idir,jf+n+1,kf+l-1)
     &                 +sigmaf(if+ii-1,jf+n  ,kf+l  ,2) * 
     &                    MYYF(if,jf+n,kf+l,if+idir,jf+n+1,kf+l+1) )
         sum = sum + fac * (
     &                 +sigmaf(if+ii-1,jf+n-1,kf+l-1,3) * 
     &                    MZZF(if,jf+n,kf+l,if+idir,jf+n-1,kf+l-1)
     &                 +sigmaf(if+ii-1,jf+n-1,kf+l  ,3) * 
     &                    MZZF(if,jf+n,kf+l,if+idir,jf+n-1,kf+l+1)
     &                 +sigmaf(if+ii-1,jf+n  ,kf+l-1,3) * 
     &                    MZZF(if,jf+n,kf+l,if+idir,jf+n+1,kf+l-1)
     &                 +sigmaf(if+ii-1,jf+n  ,kf+l  ,3) * 
     &                    MZZF(if,jf+n,kf+l,if+idir,jf+n+1,kf+l+1) )
         sum = sum + fac * idir * (
     &                 +sigmaf(if+ii-1,jf+n-1,kf+l-1,4) * 
     &                    MXZF(if,jf+n,kf+l,if+idir,jf+n-1,kf+l-1)
     &                 -sigmaf(if+ii-1,jf+n-1,kf+l  ,4) * 
     &                    MXZF(if,jf+n,kf+l,if+idir,jf+n-1,kf+l+1)
     &                 +sigmaf(if+ii-1,jf+n  ,kf+l-1,4) * 
     &                    MXZF(if,jf+n,kf+l,if+idir,jf+n+1,kf+l-1)
     &                 -sigmaf(if+ii-1,jf+n  ,kf+l  ,4) * 
     &                    MXZF(if,jf+n,kf+l,if+idir,jf+n+1,kf+l+1) )
         sum = sum + fac * (
     &                 -sigmaf(if+ii-1,jf+n-1,kf+l-1,5) * 
     &                    MYZF(if,jf+n,kf+l,if+idir,jf+n-1,kf+l-1)
     &                 +sigmaf(if+ii-1,jf+n-1,kf+l  ,5) * 
     &                    MYZF(if,jf+n,kf+l,if+idir,jf+n-1,kf+l+1)
     &                 +sigmaf(if+ii-1,jf+n  ,kf+l-1,5) * 
     &                    MYZF(if,jf+n,kf+l,if+idir,jf+n+1,kf+l-1)
     &                 -sigmaf(if+ii-1,jf+n  ,kf+l  ,5) * 
     &                    MYZF(if,jf+n,kf+l,if+idir,jf+n+1,kf+l+1) )

         cen =
     &      4.d0 * (
     &        sigmaf(if+ii-1,jf+n-1,kf+l-1,1) 
     &        + sigmaf(if+ii-1,jf+n-1,kf+l  ,1)
     &       +sigmaf(if+ii-1,jf+n  ,kf+l-1,1) 
     &        + sigmaf(if+ii-1,jf+n  ,kf+l  ,1)
     &       +sigmaf(if+ii-1,jf+n-1,kf+l-1,2) 
     &        + sigmaf(if+ii-1,jf+n-1,kf+l  ,2)
     &       +sigmaf(if+ii-1,jf+n  ,kf+l-1,2) 
     &        + sigmaf(if+ii-1,jf+n  ,kf+l  ,2)
     &       +sigmaf(if+ii-1,jf+n-1,kf+l-1,3) 
     &        + sigmaf(if+ii-1,jf+n-1,kf+l  ,3)
     &       +sigmaf(if+ii-1,jf+n  ,kf+l-1,3) 
     &        + sigmaf(if+ii-1,jf+n  ,kf+l  ,3))
     &    + 6.d0 * idir * (
     &        sigmaf(if+ii-1,jf+n-1,kf+l-1,4) 
     &        - sigmaf(if+ii-1,jf+n-1,kf+l  ,4)
     &       +sigmaf(if+ii-1,jf+n  ,kf+l-1,4) 
     &        - sigmaf(if+ii-1,jf+n  ,kf+l  ,4))
     &    + 6.d0 * (
     &        sigmaf(if+ii-1,jf+n-1,kf+l  ,5) 
     &        - sigmaf(if+ii-1,jf+n-1,kf+l-1,5)
     &       +sigmaf(if+ii-1,jf+n  ,kf+l-1,5) 
     &        - sigmaf(if+ii-1,jf+n  ,kf+l  ,5))

          sum = sum - fac * cen * fdst(if,jf+n,kf+l)
                     end do
                  end do
               end if
            end do
         end do
      end do
c edges
      do ki = 0, 1
         kdir = 2 * ki - 1
         do ji = 0, 1
            jdir = 2 * ji - 1
            do ii = 0, 1
               idir = 2 * ii - 1

               if (ga(ii,ji,ki) -
     &             min(ga(ii,ji,1-ki), ga(ii,1-ji,ki), ga(ii,1-ji,1-ki))
     &             .eq. 1) then
                  fac1 = 1.d0 / ir
                  do 60 m = idir, idir*(ir-1), idir
                     fac = (ir-abs(m)) * fac1

                  sum = sum + fac * (
     &                 +sigmaf(if+m-1,jf+ji-1,kf+ki-1,1) * 
     &                    MXXF(if+m,jf,kf,if+m-1,jf+jdir,kf+kdir)
     &                 +sigmaf(if+m  ,jf+ji-1,kf+ki-1,1) * 
     &                    MXXF(if+m,jf,kf,if+m+1,jf+jdir,kf+kdir) )
                  sum = sum + fac * (
     &                 +sigmaf(if+m-1,jf+ji-1,kf+ki-1,2) * 
     &                    MYYF(if+m,jf,kf,if+m-1,jf+jdir,kf+kdir)
     &                 +sigmaf(if+m  ,jf+ji-1,kf+ki-1,2) * 
     &                    MYYF(if+m,jf,kf,if+m+1,jf+jdir,kf+kdir) )
                  sum = sum + fac * (
     &                 +sigmaf(if+m-1,jf+ji-1,kf+ki-1,3) * 
     &                    MZZF(if+m,jf,kf,if+m-1,jf+jdir,kf+kdir)
     &                 +sigmaf(if+m  ,jf+ji-1,kf+ki-1,3) * 
     &                    MZZF(if+m,jf,kf,if+m+1,jf+jdir,kf+kdir) )
                  sum = sum + fac * kdir * (
     &                 +sigmaf(if+m-1,jf+ji-1,kf+ki-1,4) * 
     &                    MXZF(if+m,jf,kf,if+m-1,jf+jdir,kf+kdir)
     &                 -sigmaf(if+m  ,jf+ji-1,kf+ki-1,4) * 
     &                    MXZF(if+m,jf,kf,if+m+1,jf+jdir,kf+kdir) )
                  sum = sum  + fac * jdir * kdir * (
     &                 -sigmaf(if+m-1,jf+ji-1,kf+ki-1,5) * 
     &                    MYZF(if+m,jf,kf,if+m-1,jf+jdir,kf+kdir)
     &                 -sigmaf(if+m  ,jf+ji-1,kf+ki-1,5) * 
     &                    MYZF(if+m,jf,kf,if+m+1,jf+jdir,kf+kdir) )

                  cen =
     &               4.d0 * ( sigmaf(if+m-1,jf+ji-1,kf+ki-1,1)
     &                       +sigmaf(if+m  ,jf+ji-1,kf+ki-1,1)
     &                       +sigmaf(if+m-1,jf+ji-1,kf+ki-1,2)
     &                       +sigmaf(if+m  ,jf+ji-1,kf+ki-1,2)
     &                       +sigmaf(if+m-1,jf+ji-1,kf+ki-1,3)
     &                       +sigmaf(if+m  ,jf+ji-1,kf+ki-1,3))
     &             + 6.d0 * ( sigmaf(if+m-1,jf+ji-1,kf+ki-1,4)
     &                       -sigmaf(if+m  ,jf+ji-1,kf+ki-1,4)) * kdir 
     &             + 6.d0 * (-sigmaf(if+m-1,jf+ji-1,kf+ki-1,5)
     &                  -sigmaf(if+m  ,jf+ji-1,kf+ki-1,5)) * jdir * kdir

                  sum = sum - fac * cen * fdst(if+m,jf,kf)
 60               continue
               end if

               if (ga(ii,ji,ki) -
     &             min(ga(ii,ji,1-ki), ga(1-ii,ji,ki), ga(1-ii,ji,1-ki))
     &             .eq. 1) then
                  fac1 = 1.d0 / jr
                  do 70 n = jdir, jdir*(jr-1), jdir
                     fac = (jr-abs(n)) * fac1

                  sum = sum + fac * (
     &                 +sigmaf(if+ii-1,jf+n-1,kf+ki-1,1) * 
     &                    MXXF(if,jf+n,kf,if+idir,jf+n-1,kf+kdir)
     &                 +sigmaf(if+ii-1,jf+n  ,kf+ki-1,1) * 
     &                    MXXF(if,jf+n,kf,if+idir,jf+n+1,kf+kdir) )
                  sum = sum + fac * (
     &                 +sigmaf(if+ii-1,jf+n-1,kf+ki-1,2) * 
     &                    MYYF(if,jf+n,kf,if+idir,jf+n-1,kf+kdir)
     &                 +sigmaf(if+ii-1,jf+n  ,kf+ki-1,2) * 
     &                    MYYF(if,jf+n,kf,if+idir,jf+n+1,kf+kdir) )
                  sum = sum + fac * (
     &                 +sigmaf(if+ii-1,jf+n-1,kf+ki-1,3) * 
     &                    MZZF(if,jf+n,kf,if+idir,jf+n-1,kf+kdir)
     &                 +sigmaf(if+ii-1,jf+n  ,kf+ki-1,3) * 
     &                    MZZF(if,jf+n,kf,if+idir,jf+n+1,kf+kdir) )
                  sum = sum + fac * idir * kdir * (
     &                 -sigmaf(if+ii-1,jf+n-1,kf+ki-1,4) * 
     &                    MXZF(if,jf+n,kf,if+idir,jf+n-1,kf+kdir)
     &                 -sigmaf(if+ii-1,jf+n  ,kf+ki-1,4) * 
     &                    MXZF(if,jf+n,kf,if+idir,jf+n+1,kf+kdir) )
                  sum = sum  + fac * kdir * (
     &                 +sigmaf(if+ii-1,jf+n-1,kf+ki-1,5) * 
     &                    MYZF(if,jf+n,kf,if+idir,jf+n-1,kf+kdir)
     &                 -sigmaf(if+ii-1,jf+n  ,kf+ki-1,5) * 
     &                    MYZF(if,jf+n,kf,if+idir,jf+n+1,kf+kdir) )
                  cen =
     &                  4.d0 * ( sigmaf(if+ii-1,jf+n-1,kf+ki-1,1)
     &                          +sigmaf(if+ii-1,jf+n  ,kf+ki-1,1)
     &                          +sigmaf(if+ii-1,jf+n-1,kf+ki-1,2)
     &                          +sigmaf(if+ii-1,jf+n  ,kf+ki-1,2)
     &                          +sigmaf(if+ii-1,jf+n-1,kf+ki-1,3)
     &                          +sigmaf(if+ii-1,jf+n  ,kf+ki-1,3))
     &                + 6.d0 * (-sigmaf(if+ii-1,jf+n-1,kf+ki-1,4)
     &                  -sigmaf(if+ii-1,jf+n  ,kf+ki-1,4)) * idir * kdir
     &                + 6.d0 * ( sigmaf(if+ii-1,jf+n-1,kf+ki-1,5) 
     &                         -sigmaf(if+ii-1,jf+n  ,kf+ki-1,5)) * kdir
                  sum = sum - fac * cen * fdst(if,jf+n,kf)

 70               continue

               end if

               if (ga(ii,ji,ki) -
     &             min(ga(ii,1-ji,ki), ga(1-ii,ji,ki), ga(1-ii,1-ji,ki))
     &             .eq. 1) then
                  fac1 = 1.d0 / kr
                  do 80 l = kdir, kdir*(kr-1), kdir
                     fac = (kr-abs(l)) * fac1

                  sum = sum + fac * (
     &                 +sigmaf(if+ii-1,jf+ji-1,kf+l-1,1) * 
     &                    MXXF(if,jf,kf+l,if+idir,jf+jdir,kf+l-1)
     &                 +sigmaf(if+ii-1,jf+ji-1,kf+l  ,1) * 
     &                    MXXF(if,jf,kf+l,if+idir,jf+jdir,kf+l+1) )
                  sum = sum + fac * (
     &                 +sigmaf(if+ii-1,jf+ji-1,kf+l-1,2) * 
     &                    MYYF(if,jf,kf+l,if+idir,jf+jdir,kf+l-1)
     &                 +sigmaf(if+ii-1,jf+ji-1,kf+l  ,2) * 
     &                    MYYF(if,jf,kf+l,if+idir,jf+jdir,kf+l+1) )
                  sum = sum + fac * (
     &                 +sigmaf(if+ii-1,jf+ji-1,kf+l-1,3) * 
     &                    MZZF(if,jf,kf+l,if+idir,jf+jdir,kf+l-1)
     &                 +sigmaf(if+ii-1,jf+ji-1,kf+l  ,3) * 
     &                    MZZF(if,jf,kf+l,if+idir,jf+jdir,kf+l+1) )
                  sum = sum + fac * idir * (
     &                 +sigmaf(if+ii-1,jf+ji-1,kf+l-1,4) * 
     &                    MXZF(if,jf,kf+l,if+idir,jf+jdir,kf+l-1)
     &                 -sigmaf(if+ii-1,jf+ji-1,kf+l  ,4) * 
     &                    MXZF(if,jf,kf+l,if+idir,jf+jdir,kf+l+1) )
                  sum = sum + fac * jdir * (
     &                 +sigmaf(if+ii-1,jf+ji-1,kf+l-1,5) * 
     &                    MYZF(if,jf,kf+l,if+idir,jf+jdir,kf+l-1)
     &                 -sigmaf(if+ii-1,jf+ji-1,kf+l  ,5) * 
     &                    MYZF(if,jf,kf+l,if+idir,jf+jdir,kf+l+1) )
                  cen =
     &               4.d0 * (sigmaf(if+ii-1,jf+ji-1,kf+l-1,1) 
     &                      +sigmaf(if+ii-1,jf+ji-1,kf+l  ,1)
     &                      +sigmaf(if+ii-1,jf+ji-1,kf+l-1,2) 
     &                      +sigmaf(if+ii-1,jf+ji-1,kf+l  ,2)
     &                      +sigmaf(if+ii-1,jf+ji-1,kf+l-1,3) 
     &                      +sigmaf(if+ii-1,jf+ji-1,kf+l  ,3))
     &             + 6.d0 * (sigmaf(if+ii-1,jf+ji-1,kf+l-1,4) 
     &                      -sigmaf(if+ii-1,jf+ji-1,kf+l  ,4)) * idir
     &             + 6.d0 * (sigmaf(if+ii-1,jf+ji-1,kf+l-1,5) 
     &                      -sigmaf(if+ii-1,jf+ji-1,kf+l  ,5)) * jdir

                  sum = sum - fac * cen * fdst(if,jf,kf+l)

 80               continue

               end if
            end do
         end do
      end do
c weighting
      res(if,jf,kf) = src(if,jf,kf) - sum * fac00
      end
c twenty-seven-point terrain stencils
c-----------------------------------------------------------------------
      subroutine hgcen_terrain(
     & cen, cenl0,cenh0,cenl1,cenh1,cenl2,cenh2,
     & sig, sbl0,sbh0,sbl1,sbh1,sbl2,sbh2,
     &      regl0,regh0,regl1,regh1,regl2,regh2)
      integer cenl0,cenh0,cenl1,cenh1,cenl2,cenh2
      integer sbl0,sbh0,sbl1,sbh1,sbl2,sbh2
      integer regl0,regh0,regl1,regh1,regl2,regh2
      double precision cen(cenl0:cenh0,cenl1:cenh1,cenl2:cenh2)
      double precision sig(sbl0:sbh0,sbl1:sbh1,sbl2:sbh2, 5)
      integer i, j, k
      do k = regl2, regh2
         do j = regl1, regh1
            do i = regl0, regh0
               cen(i,j,k) = 36.d0 /
     &           (4.d0 * (sig(i-1,j-1,k-1,1) + sig(i-1,j-1,k,1) +
     &                    sig(i-1,j,k-1,1)   + sig(i-1,j,k,1) +
     &                    sig(i,j-1,k-1,1)   + sig(i,j-1,k,1) +
     &                    sig(i,j,k-1,1)     + sig(i,j,k,1) +
     &                    sig(i-1,j-1,k-1,2) + sig(i-1,j-1,k,2) +
     &                    sig(i-1,j,k-1,2)   + sig(i-1,j,k,2) +
     &                    sig(i,j-1,k-1,2)   + sig(i,j-1,k,2) +
     &                    sig(i,j,k-1,2)     + sig(i,j,k,2) +
     &                    sig(i-1,j-1,k-1,3) + sig(i-1,j-1,k,3) +
     &                    sig(i-1,j,k-1,3)   + sig(i-1,j,k,3) +
     &                    sig(i,j-1,k-1,3)   + sig(i,j-1,k,3) +
     &                    sig(i,j,k-1,3)     + sig(i,j,k,3)) +
     &            6.d0 * (sig(i-1,j-1,k,4)   - sig(i-1,j-1,k-1,4) +
     &                    sig(i-1,j,k,4)     - sig(i-1,j,k-1,4) +
     &                    sig(i,j-1,k-1,4)   - sig(i,j-1,k,4) +
     &                    sig(i,j,k-1,4)     - sig(i,j,k,4) +
     &                    sig(i-1,j-1,k,5)   - sig(i-1,j-1,k-1,5) +
     &                    sig(i-1,j,k-1,5)   - sig(i-1,j,k,5) +
     &                    sig(i,j-1,k,5)     - sig(i,j-1,k-1,5) +
     &                    sig(i,j,k-1,5)     - sig(i,j,k,5)))
            end do
         end do
      end do
      end
c-----------------------------------------------------------------------
      subroutine hgrlx_terrain(
     & cor, corl0,corh0,corl1,corh1,corl2,corh2,
     & res, resl0,resh0,resl1,resh1,resl2,resh2,
     & sig, sfl0,sfh0,sfl1,sfh1,sfl2,sfh2,
     & cen, cenl0,cenh0,cenl1,cenh1,cenl2,cenh2,
     &      regl0,regh0,regl1,regh1,regl2,regh2)
      integer corl0,corh0,corl1,corh1,corl2,corh2
      integer resl0,resh0,resl1,resh1,resl2,resh2
      integer sfl0,sfh0,sfl1,sfh1,sfl2,sfh2
      integer cenl0,cenh0,cenl1,cenh1,cenl2,cenh2
      integer regl0,regh0,regl1,regh1,regl2,regh2
      double precision cor(corl0:corh0,corl1:corh1,corl2:corh2)
      double precision res(resl0:resh0,resl1:resh1,resl2:resh2)
      double precision sig(sfl0:sfh0,sfl1:sfh1,sfl2:sfh2, 5)
      double precision cen(cenl0:cenh0,cenl1:cenh1,cenl2:cenh2)
      double precision fac, tmp
      integer i, j, k
      double precision MXX, MYY, MZZ, MXZ, MYZ
      MXX(ii,jj,kk) = (4.d0 *  cor(ii,j,k)   + 
     &                  2.d0 * (cor(ii,jj,k)  - cor(i,jj,k) + 
     &                          cor(ii,j,kk)  - cor(i,j,kk)) + 
     &                         (cor(ii,jj,kk) - cor(i,jj,kk)))

      MYY(ii,jj,kk) = (4.d0 *  cor(i,jj,k)   + 
     &                  2.d0 * (cor(ii,jj,k)  - cor(ii,j,k) + 
     &                          cor(i,jj,kk)  - cor(i,j,kk)) + 
     &                         (cor(ii,jj,kk) - cor(ii,j,kk)))

      MZZ(ii,jj,kk) = (4.d0 *  cor(i,j,kk)   + 
     &                  2.d0 * (cor(ii,j,kk)  - cor(ii,j,k) + 
     &                          cor(i,jj,kk)  - cor(i,jj,k)) + 
     &                         (cor(ii,jj,kk) - cor(ii,jj,k)))

      MXZ(ii,jj,kk) = (6.d0 *  cor(ii,j,kk)  + 
     &                  3.d0 * (cor(ii,jj,kk) - cor(i,jj,k)))

      MYZ(ii,jj,kk) = (6.d0 *  cor(i,jj,kk)  + 
     &                  3.d0 * (cor(ii,jj,kk) - cor(ii,j,k)))
      fac = 1.d0 / 36.d0
         do k = regl2, regh2
            do j = regl1, regh1
cdir$ ivdep
               do i = regl0, regh0
                  tmp = sig(i-1,j-1,k-1,1) * MXX(i-1,j-1,k-1) +
     &                  sig(i-1,j-1,k  ,1) * MXX(i-1,j-1,k+1) +
     &                  sig(i-1,j  ,k-1,1) * MXX(i-1,j+1,k-1) +
     &                  sig(i-1,j  ,k  ,1) * MXX(i-1,j+1,k+1)
                  tmp = tmp +
     &                  sig(i  ,j-1,k-1,1) * MXX(i+1,j-1,k-1) +
     &                  sig(i  ,j-1,k  ,1) * MXX(i+1,j-1,k+1) +
     &                  sig(i  ,j  ,k-1,1) * MXX(i+1,j+1,k-1) +
     &                  sig(i  ,j  ,k  ,1) * MXX(i+1,j+1,k+1)
                  tmp = tmp +
     &                  sig(i-1,j-1,k-1,2) * MYY(i-1,j-1,k-1) +
     &                  sig(i-1,j-1,k  ,2) * MYY(i-1,j-1,k+1) +
     &                  sig(i-1,j  ,k-1,2) * MYY(i-1,j+1,k-1) +
     &                  sig(i-1,j  ,k  ,2) * MYY(i-1,j+1,k+1)
                  tmp = tmp +
     &                  sig(i  ,j-1,k-1,2) * MYY(i+1,j-1,k-1) +
     &                  sig(i  ,j-1,k  ,2) * MYY(i+1,j-1,k+1) +
     &                  sig(i  ,j  ,k-1,2) * MYY(i+1,j+1,k-1) +
     &                  sig(i  ,j  ,k  ,2) * MYY(i+1,j+1,k+1)
                  tmp = tmp +
     &                  sig(i-1,j-1,k-1,3) * MZZ(i-1,j-1,k-1) +
     &                  sig(i-1,j-1,k  ,3) * MZZ(i-1,j-1,k+1) +
     &                  sig(i-1,j  ,k-1,3) * MZZ(i-1,j+1,k-1) +
     &                  sig(i-1,j  ,k  ,3) * MZZ(i-1,j+1,k+1)
                  tmp = tmp +
     &                  sig(i  ,j-1,k-1,3) * MZZ(i+1,j-1,k-1) +
     &                  sig(i  ,j-1,k  ,3) * MZZ(i+1,j-1,k+1) +
     &                  sig(i  ,j  ,k-1,3) * MZZ(i+1,j+1,k-1) +
     &                  sig(i  ,j  ,k  ,3) * MZZ(i+1,j+1,k+1)
                  tmp = tmp -
     &                  sig(i-1,j-1,k-1,4) * MXZ(i-1,j-1,k-1) +
     &                  sig(i-1,j-1,k  ,4) * MXZ(i-1,j-1,k+1) -
     &                  sig(i-1,j  ,k-1,4) * MXZ(i-1,j+1,k-1) +
     &                  sig(i-1,j  ,k  ,4) * MXZ(i-1,j+1,k+1) +
     &                  sig(i  ,j-1,k-1,4) * MXZ(i+1,j-1,k-1) -
     &                  sig(i  ,j-1,k  ,4) * MXZ(i+1,j-1,k+1) +
     &                  sig(i  ,j  ,k-1,4) * MXZ(i+1,j+1,k-1) -
     &                  sig(i  ,j  ,k  ,4) * MXZ(i+1,j+1,k+1)
                  tmp = tmp -
     &                  sig(i-1,j-1,k-1,5) * MYZ(i-1,j-1,k-1) +
     &                  sig(i-1,j-1,k  ,5) * MYZ(i-1,j-1,k+1) +
     &                  sig(i-1,j  ,k-1,5) * MYZ(i-1,j+1,k-1) -
     &                  sig(i-1,j  ,k  ,5) * MYZ(i-1,j+1,k+1) -
     &                  sig(i  ,j-1,k-1,5) * MYZ(i+1,j-1,k-1) +
     &                  sig(i  ,j-1,k  ,5) * MYZ(i+1,j-1,k+1) +
     &                  sig(i  ,j  ,k-1,5) * MYZ(i+1,j+1,k-1) -
     &                  sig(i  ,j  ,k  ,5) * MYZ(i+1,j+1,k+1)
                  cor(i,j,k) = (fac * tmp - res(i,j,k)) * cen(i,j,k)
               end do
            end do
         end do
      end

c cor substituted for dest to use same macro:
c-----------------------------------------------------------------------
      subroutine hgres_terrain(
     & res, resl0,resh0,resl1,resh1,resl2,resh2,
     & src, srcl0,srch0,srcl1,srch1,srcl2,srch2,
     & cor, destl0,desth0,destl1,desth1,destl2,desth2,
     & sig, sfl0,sfh0,sfl1,sfh1,sfl2,sfh2,
     & cen, cenl0,cenh0,cenl1,cenh1,cenl2,cenh2,
     &      regl0,regh0,regl1,regh1,regl2,regh2,
     & hx, hy, hz)
      integer resl0,resh0,resl1,resh1,resl2,resh2
      integer srcl0,srch0,srcl1,srch1,srcl2,srch2
      integer destl0,desth0,destl1,desth1,destl2,desth2
      integer sfl0,sfh0,sfl1,sfh1,sfl2,sfh2
      integer cenl0,cenh0,cenl1,cenh1,cenl2,cenh2
      integer regl0,regh0,regl1,regh1,regl2,regh2
      double precision hx, hy, hz
      double precision res(resl0:resh0,resl1:resh1,resl2:resh2)
      double precision src(srcl0:srch0,srcl1:srch1,srcl2:srch2)
      double precision cor(destl0:desth0,destl1:desth1,destl2:desth2)
      double precision sig(sfl0:sfh0,sfl1:sfh1,sfl2:sfh2, 5)
      double precision cen(cenl0:cenh0,cenl1:cenh1,cenl2:cenh2)
      double precision fac, tmp
      integer i, j, k
      double precision MXX, MYY, MZZ, MXZ, MYZ
      MXX(ii,jj,kk) = (4.d0 *  cor(ii,j,k)   + 
     &                  2.d0 * (cor(ii,jj,k)  - cor(i,jj,k) + 
     &                          cor(ii,j,kk)  - cor(i,j,kk)) + 
     &                         (cor(ii,jj,kk) - cor(i,jj,kk)))

      MYY(ii,jj,kk) = (4.d0 *  cor(i,jj,k)   + 
     &                  2.d0 * (cor(ii,jj,k)  - cor(ii,j,k) + 
     &                          cor(i,jj,kk)  - cor(i,j,kk)) + 
     &                         (cor(ii,jj,kk) - cor(ii,j,kk)))

      MZZ(ii,jj,kk) = (4.d0 *  cor(i,j,kk)   + 
     &                  2.d0 * (cor(ii,j,kk)  - cor(ii,j,k) + 
     &                          cor(i,jj,kk)  - cor(i,jj,k)) + 
     &                         (cor(ii,jj,kk) - cor(ii,jj,k)))

      MXZ(ii,jj,kk) = (6.d0 *  cor(ii,j,kk)  + 
     &                  3.d0 * (cor(ii,jj,kk) - cor(i,jj,k)))

      MYZ(ii,jj,kk) = (6.d0 *  cor(i,jj,kk)  + 
     &                  3.d0 * (cor(ii,jj,kk) - cor(ii,j,k)))
      fac = 1.d0 / 36.d0
         do k = regl2, regh2
            do j = regl1, regh1
               do i = regl0, regh0
                  tmp = sig(i-1,j-1,k-1,1) * MXX(i-1,j-1,k-1) +
     &                  sig(i-1,j-1,k  ,1) * MXX(i-1,j-1,k+1) +
     &                  sig(i-1,j  ,k-1,1) * MXX(i-1,j+1,k-1) +
     &                  sig(i-1,j  ,k  ,1) * MXX(i-1,j+1,k+1)
                  tmp = tmp +
     &                  sig(i  ,j-1,k-1,1) * MXX(i+1,j-1,k-1) +
     &                  sig(i  ,j-1,k  ,1) * MXX(i+1,j-1,k+1) +
     &                  sig(i  ,j  ,k-1,1) * MXX(i+1,j+1,k-1) +
     &                  sig(i  ,j  ,k  ,1) * MXX(i+1,j+1,k+1)
                  tmp = tmp +
     &                  sig(i-1,j-1,k-1,2) * MYY(i-1,j-1,k-1) +
     &                  sig(i-1,j-1,k  ,2) * MYY(i-1,j-1,k+1) +
     &                  sig(i-1,j  ,k-1,2) * MYY(i-1,j+1,k-1) +
     &                  sig(i-1,j  ,k  ,2) * MYY(i-1,j+1,k+1)
                  tmp = tmp +
     &                  sig(i  ,j-1,k-1,2) * MYY(i+1,j-1,k-1) +
     &                  sig(i  ,j-1,k  ,2) * MYY(i+1,j-1,k+1) +
     &                  sig(i  ,j  ,k-1,2) * MYY(i+1,j+1,k-1) +
     &                  sig(i  ,j  ,k  ,2) * MYY(i+1,j+1,k+1)
                  tmp = tmp +
     &                  sig(i-1,j-1,k-1,3) * MZZ(i-1,j-1,k-1) +
     &                  sig(i-1,j-1,k  ,3) * MZZ(i-1,j-1,k+1) +
     &                  sig(i-1,j  ,k-1,3) * MZZ(i-1,j+1,k-1) +
     &                  sig(i-1,j  ,k  ,3) * MZZ(i-1,j+1,k+1)
                  tmp = tmp +
     &                  sig(i  ,j-1,k-1,3) * MZZ(i+1,j-1,k-1) +
     &                  sig(i  ,j-1,k  ,3) * MZZ(i+1,j-1,k+1) +
     &                  sig(i  ,j  ,k-1,3) * MZZ(i+1,j+1,k-1) +
     &                  sig(i  ,j  ,k  ,3) * MZZ(i+1,j+1,k+1)
                  tmp = tmp -
     &                  sig(i-1,j-1,k-1,4) * MXZ(i-1,j-1,k-1) +
     &                  sig(i-1,j-1,k  ,4) * MXZ(i-1,j-1,k+1) -
     &                  sig(i-1,j  ,k-1,4) * MXZ(i-1,j+1,k-1) +
     &                  sig(i-1,j  ,k  ,4) * MXZ(i-1,j+1,k+1) +
     &                  sig(i  ,j-1,k-1,4) * MXZ(i+1,j-1,k-1) -
     &                  sig(i  ,j-1,k  ,4) * MXZ(i+1,j-1,k+1) +
     &                  sig(i  ,j  ,k-1,4) * MXZ(i+1,j+1,k-1) -
     &                  sig(i  ,j  ,k  ,4) * MXZ(i+1,j+1,k+1)
                  tmp = tmp -
     &                  sig(i-1,j-1,k-1,5) * MYZ(i-1,j-1,k-1) +
     &                  sig(i-1,j-1,k  ,5) * MYZ(i-1,j-1,k+1) +
     &                  sig(i-1,j  ,k-1,5) * MYZ(i-1,j+1,k-1) -
     &                  sig(i-1,j  ,k  ,5) * MYZ(i-1,j+1,k+1) -
     &                  sig(i  ,j-1,k-1,5) * MYZ(i+1,j-1,k-1) +
     &                  sig(i  ,j-1,k  ,5) * MYZ(i+1,j-1,k+1) +
     &                  sig(i  ,j  ,k-1,5) * MYZ(i+1,j+1,k-1) -
     &                  sig(i  ,j  ,k  ,5) * MYZ(i+1,j+1,k+1)
                  if (cen(i,j,k) .ne. 0.d0) then
                     res(i,j,k) = src(i,j,k)
     &                  - fac * tmp + cor(i,j,k) / cen(i,j,k)
                  else
                     res(i,j,k) = 0.d0
                  end if
               end do
            end do
         end do
      end
c-----------------------------------------------------------------------
      subroutine hgrlnf_terrain(
     & cor,   corl0,corh0,corl1,corh1,corl2,corh2,
     & res,   resl0,resh0,resl1,resh1,resl2,resh2,
     & wrk,   wrkl0,wrkh0,wrkl1,wrkh1,wrkl2,wrkh2,
     & signd, snl0,snh0,snl1,snh1,snl2,snh2,
     & cen,   cenl0,cenh0,cenl1,cenh1,cenl2,cenh2,
     &        regl0,regh0,regl1,regh1,regl2,regh2,
     &        doml0,domh0,doml1,domh1,doml2,domh2,
     & lsd, ipass)
      integer corl0,corh0,corl1,corh1,corl2,corh2
      integer resl0,resh0,resl1,resh1,resl2,resh2
      integer wrkl0,wrkh0,wrkl1,wrkh1,wrkl2,wrkh2
      integer snl0,snh0,snl1,snh1,snl2,snh2
      integer cenl0,cenh0,cenl1,cenh1,cenl2,cenh2
      integer regl0,regh0,regl1,regh1,regl2,regh2
      integer doml0,domh0,doml1,domh1,doml2,domh2
      double precision cor(corl0:corh0,corl1:corh1,corl2:corh2)
      double precision res(resl0:resh0,resl1:resh1,resl2:resh2)
      double precision wrk(wrkl0:wrkh0,wrkl1:wrkh1,wrkl2:wrkh2)
      double precision signd(snl0:snh0,snl1:snh1,snl2:snh2, 3)
      double precision cen(cenl0:cenh0,cenl1:cenh1,cenl2:cenh2)
      integer lsd, ipass
      stop 'hgrlnf_terrain: no code'
      end
