       subroutine comp(imax,jmax,dx,dy,fluxx,fluxy,ebdiffop,
     1  delm,divc,vfrac,apx,apy,centx,centy,itype,ndim)

       implicit none

       double precision dx,dy
       integer imax,jmax,ndim

       double precision fluxy(0:ndim,0:ndim)
       double precision fluxx(0:ndim,0:ndim)
       double precision ebdiffop(0:ndim,0:ndim)
       double precision delm(0:ndim,0:ndim)
       double precision divc(0:ndim,0:ndim)
       double precision vfrac(0:ndim,0:ndim)
       double precision apx(0:ndim,0:ndim)
       double precision apy(0:ndim,0:ndim)
       double precision centx(0:ndim,0:ndim)
       double precision centy(0:ndim,0:ndim)

       integer itype(0:ndim,0:ndim)

       double precision frac,fxp,fxm,fyp,fym,divnc,vtot
       double precision testint

       integer i,j,ii,jj

c  OLD VERSION
c  assumed update was U_t = (\nabla \cdot F) ; not  with minus sign
c  assumed centroids wrt bottom of cell


       do j=1,jmax
       do i=1,imax

           if(itype(i,j) .eq. -1)then
              divc(i,j) = 0.d0
           elseif(itype(i,j).eq. 0)then
              divc(i,j) = (fluxx(i+1,j) - fluxx(i,j))/dx
     1             + (fluxy(i,j+1)-fluxy(i,j))/dy
           else

              if(apx(i,j).lt.1.d0)then
                 if(centx(i,j).le. 0.5d0)then
                     frac = 0.5d0 - centx(i,j)
                     fxm = frac*fluxx(i,j-1)+(1.d0-frac)*fluxx(i,j)
                 else
                     frac = centx(i,j) - 0.5d0
                     fxm = frac*fluxx(i,j+1)+(1.d0-frac)*fluxx(i,j)
                 endif
              else
                 fxm = fluxx(i,j)
              endif

              if(apx(i+1,j).lt.1.d0)then
                 if(centx(i+1,j).le. 0.5d0)then
                     frac = 0.5d0 - centx(i+1,j)
                     fxp = frac*fluxx(i+1,j-1)+(1.d0-frac)*fluxx(i+1,j)
                 else
                     frac = centx(i+1,j) - 0.5d0
                     fxp = frac*fluxx(i+1,j+1)+(1.d0-frac)*fluxx(i+1,j)
                 endif
              else
                 fxp = fluxx(i+1,j)
              endif

              if(apy(i,j).lt.1.d0)then
                 if(centy(i,j).le. 0.5d0)then
                     frac = 0.5d0 - centy(i,j)
                     fym = frac*fluxy(i-1,j)+(1.d0-frac)*fluxy(i,j)
                 else
                     frac = centy(i,j) - 0.5d0
                     fym = frac*fluxy(i+1,j)+(1.d0-frac)*fluxy(i,j)
                 endif
              else
                 fym = fluxy(i,j)
              endif

              if(apy(i,j+1).lt.1.d0)then
                 if(centy(i,j+1).le. 0.5d0)then
                     frac = 0.5d0 - centy(i,j+1)
                     fyp = frac*fluxy(i-1,j+1)+(1.d0-frac)*fluxy(i,j+1)
                 else
                     frac = centy(i,j+1) - 0.5d0
                     fyp = frac*fluxy(i+1,j+1)+(1.d0-frac)*fluxy(i,j+1)
                 endif
              else
                 fyp = fluxy(i,j+1)
              endif

              divc(i,j) =((apx(i+1,j)*fxp - apx(i,j)*fxm) / dx
     1                  +(apy(i,j+1)*fyp - apy(i,j)*fym) / dy) /
     2                    vfrac(i,j)


           endif

        enddo
        enddo

        testint = 0.d0
        do j=1,jmax
        do i=1,imax
             testint = testint + divc(i,j)*vfrac(i,j)*dx*dy
        enddo  
        enddo
        write(6,*)" integral of divc", testint

c       write(6,*)" before definito"
c       do j=13,18
c       do i=9,12
c          write(6,*)i,j,itype(i,j),vfrac(i,j),ebdiffop(i,j),divc(i,j)
c       enddo
c       enddo
              
        do j=1,jmax
        do i=1,imax
          if(itype(i,j).ne.1)then
            ebdiffop(i,j) = divc(i,j)
          else
            vtot = vfrac(i-1,j-1)+vfrac(i-1,j)+vfrac(i-1,j+1)
     1            + vfrac(i,j-1)+vfrac(i,j)+vfrac(i,j+1)
     2             + vfrac(i+1,j-1)+vfrac(i+1,j)+vfrac(i+1,j+1)
            divnc = 0.d0
            do jj = -1,1
            do ii = -1,1
               divnc = divnc + vfrac(i+ii,j+jj)*divc(i+ii,j+jj)
            enddo
            enddo
            divnc = divnc / vtot
            ebdiffop(i,j) = vfrac(i,j)*divc(i,j)+(1.d0-vfrac(i,j))*divnc
            delm(i,j) = vfrac(i,j)*(1.d0-vfrac(i,j))*(divc(i,j)-divnc)
          endif
        enddo
        enddo

c       write(6,*)" before redist"
c       do j=13,18
c       do i=9,12
c          write(6,*)i,j,itype(i,j),vfrac(i,j),ebdiffop(i,j),divc(i,j)
c       enddo
c       enddo

        do j=1,jmax
        do i=1,imax

           if(itype(i,j) .eq. 1)then
            vtot = vfrac(i-1,j-1)+vfrac(i-1,j)+vfrac(i-1,j+1)
     1            + vfrac(i,j-1)+vfrac(i,j+1)
     2            + vfrac(i+1,j-1)+vfrac(i+1,j)+vfrac(i+1,j+1)
            do jj = -1,1
            do ii = -1,1
            if((ii.ne. 0 .or. jj.ne.0).and.itype(i+ii,j+jj).ne. -1)then
                 ebdiffop(i+ii,j+jj) = ebdiffop(i+ii,j+jj)+
     1           delm(i,j)/vtot
c    1           vfrac(i+ii,j+jj)*delm(i,j)/vtot
               endif
            enddo
            enddo


           endif
      enddo
      enddo

c       write(6,*)" after redist"
c       do j=13,18
c       do i=9,12
c          write(6,*)i,j,itype(i,j),vfrac(i,j),ebdiffop(i,j)
c       enddo
c       enddo


      return
      end
