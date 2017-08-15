       subroutine comp3d(imax,jmax,kmax,dx,dy,dz,fluxx,fluxy,fluxz,
     1   ebdiffop, delm,divc,vfrac,apx,apy,apz,centx_y,centx_z_,centy_x,
     1   centy_z,centz_x,centz_y,itype,ndim)

       implicit none

       double precision dx,dy,dz
       integer imax,jmax,kmax,ndim

       double precision fluxy(0:ndim,0:ndim,0:ndim)
       double precision fluxx(0:ndim,0:ndim,0:ndim)
       double precision ebdiffop(0:ndim,0:ndim,0:ndim)
       double precision delm(0:ndim,0:ndim,0:ndim)
       double precision divc(0:ndim,0:ndim,0:ndim)
       double precision vfrac(0:ndim,0:ndim,0:ndim)
       double precision apx(0:ndim,0:ndim,0:ndim)
       double precision apy(0:ndim,0:ndim,0:ndim)
       double precision apz(0:ndim,0:ndim,0:ndim)
       double precision centx_y(0:ndim,0:ndim,0:ndim)
       double precision centx_z(0:ndim,0:ndim,0:ndim)
       double precision centy_x(0:ndim,0:ndim,0:ndim)
       double precision centy_z(0:ndim,0:ndim,0:ndim)
       double precision centz_x(0:ndim,0:ndim,0:ndim)
       double precision centz_y(0:ndim,0:ndim,0:ndim)

       integer itype(0:ndim,0:ndim,0:ndim)

       double precision fxp,fxm,fyp,fym,fzp,fzm,divnc,vtot
       double precision fracx,fracy,fracz
       double precision testint

       integer i,j,k,ii,jj,kk

c     itype -1 for body , 0 for regular, 1 for mixed


       do k=1,kmax
       do j=1,jmax
       do i=1,imax

           if(itype(i,j,k) .eq. -1)then
              divc(i,j,k) = 0.d0
           elseif(itype(i,j,k).eq. 0)then
              divc(i,j,k) = - ( (fluxx(i+1,j,k) - fluxx(i,j,k))/dx
     1             + (fluxy(i,j+1,k)-fluxy(i,j,k))/dy  
     1             + (fluxz(i,j,k+1)-fluxz(i,j,k))/dz )
           else

              if(apx(i,j,k).lt.1.d0)then
c  centroid nondimensional  and zero at face center
                 if(centx_y(i,j,k).le. 0.0d0)then
                     fracy = - centx_y(i,j,k)
                     if(centx_z(i,j,k).le. 0.0d0)then
                       fracz = - centx_z(i,j,k)
                       fxm = (1.d0-fracz)*(fracy*fluxx(i,j-1,k)+
     1                       (1.d0-fracy)*fluxx(i,j,k)) + 
     1                  fracz*( fracy*fluxx(i,j-1,k-1)+
     1                   (1.d0-fracy)*fluxx(i,j,k-1))
                     else
                       fracz =  centx_z(i,j,k)
                       fxm = (1.d0-fracz)*(fracy*fluxx(i,j-1,k)+
     1                       (1.d0-fracy)*fluxx(i,j,k)) + 
     1                  fracz*( fracy*fluxx(i,j-1,k+1)+
     1                   (1.d0-fracy)*fluxx(i,j,k+1))
                     endif
                 else
                     fracy =  centx_y(i,j,k)
                     if(centx_z(i,j,k).le. 0.0d0)then
                        fracz = -centx_z(i,j,k) 
                        fxm = (1.d0-fracz)*(fracy*fluxx(i,j+1,k)+
     1                       (1.d0-fracy)*fluxx(i,j,k)) + 
     1                  fracz*( fracy*fluxx(i,j+1,k-1)+
     1                   (1.d0-fracy)*fluxx(i,j,k-1))
                     else
                        fracz = centx_z(i,j,k) 
                        fxm = (1.d0-fracz)*(fracy*fluxx(i,j+1,k)+
     1                       (1.d0-fracy)*fluxx(i,j,k)) + 
     1                  fracz*( fracy*fluxx(i,j+1,k+1)+
     1                   (1.d0-fracy)*fluxx(i,j,k+1))

                     endif
                 endif
              else
                 fxm = fluxx(i,j,k)
              endif

              if(apx(i+1,j,k).lt.1.d0)then
c  centroid nondimensional  and zero at face center
                 if(centx_y(i+1,j,k).le. 0.0d0)then
                     fracy = - centx_y(i+1,j,k)
                     if(centx_z(i+1,j,k).le. 0.0d0)then
                       fracz = - centx_z(i+1,j,k)
                       fxp = (1.d0-fracz)*(fracy*fluxx(i+1,j-1,k)+
     1                       (1.d0-fracy)*fluxx(i+1,j,k)) + 
     1                  fracz*( fracy*fluxx(i+1,j-1,k-1)+
     1                   (1.d0-fracy)*fluxx(i+1,j,k-1))
                     else
                       fracz =  centx_z(i+1,j,k)
                       fxp = (1.d0-fracz)*(fracy*fluxx(i+1,j-1,k)+
     1                       (1.d0-fracy)*fluxx(i+1,j,k)) + 
     1                  fracz*( fracy*fluxx(i+1,j-1,k+1)+
     1                   (1.d0-fracy)*fluxx(i+1,j,k+1))
                     endif
                 else
                     fracy =  centx_y(i+1,j,k)
                     if(centx_z(i+1,j,k).le. 0.0d0)then
                        fracz = -centx_z(i+1,j,k) 
                        fxp = (1.d0-fracz)*(fracy*fluxx(i+1,j+1,k)+
     1                       (1.d0-fracy)*fluxx(i+1,j,k)) + 
     1                  fracz*( fracy*fluxx(i+1,j+1,k-1)+
     1                   (1.d0-fracy)*fluxx(i+1,j,k-1))
                     else
                        frpcz = centx_z(i+1,j,k) 
                        fxp = (1.d0-fracz)*(fracy*fluxx(i+1,j+1,k)+
     1                       (1.d0-fracy)*fluxx(i+1,j,k)) + 
     1                  fracz*( fracy*fluxx(i+1,j+1,k+1)+
     1                   (1.d0-fracy)*fluxx(i+1,j,k+1))

                     endif
                 endif
              else
                 fxp = fluxx(i+1,j,k)
              endif

              if(apy(i,j,k).lt.1.d0)then
c  centroid nondimensional  and zero at face center
                 if(centy_x(i,j,k).le. 0.0d0)then
                     fracx = - centy_x(i,j,k)
                     if(centy_z(i,j,k).le. 0.0d0)then
                       fracz = - centy_z(i,j,k)
                       fym = (1.d0-fracz)*(fracx*fluxy(i-1,j,k)+
     1                       (1.d0-fracx)*fluxy(i,j,k)) + 
     1                  fracz*( fracx*fluxy(i-1,j,k-1)+
     1                   (1.d0-fracx)*fluxy(i,j,k-1))
                     else
                       fracz =  centy_z(i,j,k)
                       fym = (1.d0-fracz)*(fracx*fluxy(i-1,j,k)+
     1                       (1.d0-fracx)*fluxy(i,j,k)) + 
     1                  fracz*( fracx*fluxy(i-1,j,k+1)+
     1                   (1.d0-fracyx*fluxy(i,j,k+1))
                     endif
                 else
                     fracx =  centy_x(i,j,k)
                     if(centy_z(i,j,k).le. 0.0d0)then
                        fracz = -centx_z(i,j,k) 
                        fym = (1.d0-fracz)*(fracy*fluxy(i+1,j,k)+
     1                       (1.d0-fracy)*fluxy(i,j,k)) + 
     1                  fracz*( fracy*fluxy(i+1,j,k-1)+
     1                   (1.d0-fracy)*fluxy(i,j,k-1))
                     else
                        fracz = centx_z(i,j,k) 
                        fym = (1.d0-fracz)*(fracy*fluxy(i+1,j,k)+
     1                       (1.d0-fracy)*fluxy(i,j,k)) + 
     1                  fracz*( fracy*fluxy(i+1,j,k+1)+
     1                   (1.d0-fracy)*fluxy(i,j,k+1))

                     endif
                 endif
              else
                 fym = fluxy(i,j,k)
              endif

              if(apy(i,j+1,k).lt.1.d0)then
c  centroid nondimensional  and zero at face center
                 if(centy_x(i,j+1,k).le. 0.0d0)then
                     fracx = - centy_x(i,j+1,k)
                     if(centy_z(i,j+1,k).le. 0.0d0)then
                       fracz = - centy_z(i,j+1,k)
                       fyp = (1.d0-fracz)*(fracx*fluxy(i-1,j+1,k)+
     1                       (1.d0-fracx)*fluxy(i,j+1,k)) + 
     1                  fracz*( fracx*fluxy(i-1,j+1,k-1)+
     1                   (1.d0-fracx)*fluxy(i,j+1,k-1))
                     else
                       fracz =  centy_z(i,j+1,k)
                       fyp = (1.d0-fracz)*(fracx*fluxy(i-1,j+1,k)+
     1                       (1.d0-fracx)*fluxy(i,j+1,k)) + 
     1                  fracz*( fracx*fluxy(i-1,j+1,k+1)+
     1                   (1.d0-fracyx*fluxy(i,j+1,k+1))
                     endif
                 else
                     fracx =  centy_x(i,j+1,k)
                     if(centy_z(i,j+1,k).le. 0.0d0)then
                        fracz = -centx_z(i,j+1,k) 
                        fyp = (1.d0-fracz)*(fracy*fluxy(i+1,j+1,k)+
     1                       (1.d0-fracy)*fluxy(i,j+1,k)) + 
     1                  fracz*( fracy*fluxy(i+1,j+1,k-1)+
     1                   (1.d0-fracy)*fluxy(i,j+1,k-1))
                     else
                        fracz = centx_z(i,j+1,k) 
                        fyp = (1.d0-fracz)*(fracy*fluxy(i+1,j+1,k)+
     1                       (1.d0-fracy)*fluxy(i,j+1,k)) + 
     1                  fracz*( fracy*fluxy(i,j+1,k+1)+
     1                   (1.d0-fracy)*fluxy(i,j+1,k+1))

                     endif
                 endif
              else
                 fyp = fluxy(i,j+1,k)
              endif

              if(apz(i,j,k).lt.1.d0)then
c  centroid nondimensional  and zero at face center
                 if(centz_x(i,j,k).le. 0.0d0)then
                     fracx = - centz_x(i,j,k)
                     if(centz_y(i,j,k).le. 0.0d0)then
                       fracy = - centz_y(i,j,k)
                       fzm = (1.d0-fracy)*(fracx*fluxy(i-1,j,k)+
     1                       (1.d0-fracx)*fluxy(i,j,k)) + 
     1                  fracy*( fracx*fluxy(i-1,j-1,k)+
     1                   (1.d0-fracx)*fluxy(i,j-1,k))
                     else
                       fracy =  centz_y(i,j,k)
                       fzm = (1.d0-fracy)*(fracx*fluxy(i-1,j,k)+
     1                       (1.d0-fracx)*fluxy(i,j,k)) + 
     1                  fracy*( fracx*fluxy(i-1,j+1,k)+
     1                   (1.d0-fracyx*fluxy(i,j+1,k))
                     endif
                 else
                     fracx =  centz_x(i,j,k)
                     if(centz_y(i,j,k).le. 0.0d0)then
                        fracy = -centz_y(i,j,k) 
                        fzm = (1.d0-fracy)*(fracx*fluxz(i+1,j,k)+
     1                       (1.d0-fracx)*fluxz(i,j,k)) + 
     1                  fracy*( fracx*fluxz(i+1,j-1,k)+
     1                   (1.d0-fracx)*fluxz(i,j-1,k))
                     else
                        fracy = centz_y(i,j,k) 
                        fzm = (1.d0-fracy)*(fracx*fluxz(i+1,j,k)+
     1                       (1.d0-fracx)*fluxz(i,j,k)) + 
     1                  fracy*( fracx*fluxz(i+1,j+1,k)+
     1                   (1.d0-fracx)*fluxz(i,j+1,k))

                     endif
                 endif
              else
                 fzm = fluxz(i,j,k)
              endif

              if(apz(i,j,k+1).lt.1.d0)then
c  centroid nondimensional  and zero at face center
                 if(centz_x(i,j,k+1).le. 0.0d0)then
                     fracx = - centz_x(i,j,k+1)
                     if(centz_y(i,j,k+1).le. 0.0d0)then
                       fracy = - centz_y(i,j,k+1)
                       fzp = (1.d0-fracy)*(fracx*fluxy(i-1,j,k+1)+
     1                       (1.d0-fracx)*fluxy(i,j,k+1)) + 
     1                  fracy*( fracx*fluxy(i-1,j-1,k+1)+
     1                   (1.d0-fracx)*fluxy(i,j-1,k+1))
                     else
                       fracy =  centz_y(i,j,k+1)
                       fzp = (1.d0-fracy)*(fracx*fluxy(i-1,j,k+1)+
     1                       (1.d0-fracx)*fluxy(i,j,k+1)) + 
     1                  fracy*( fracx*fluxy(i-1,j+1,k+1)+
     1                   (1.d0-fracyx*fluxy(i,j+1,k+1))
                     endif
                 else
                     fracx =  centz_x(i,j,k+1)
                     if(centz_y(i,j,k+1).le. 0.0d0)then
                        fracy = -centz_y(i,j,k+1) 
                        fzp = (1.d0-fracy)*(fracx*fluxz(i+1,j,k+1)+
     1                       (1.d0-fracx)*fluxz(i,j,k+1)) + 
     1                  fracy*( fracx*fluxz(i+1,j-1,k+1)+
     1                   (1.d0-fracx)*fluxz(i,j-1,k+1))
                     else
                        fracy = centz_y(i,j,k+1) 
                        fzp = (1.d0-fracy)*(fracx*fluxz(i+1,j,k+1)+
     1                       (1.d0-fracx)*fluxz(i,j,k+1)) + 
     1                  fracy*( fracx*fluxz(i+1,j+1,k+1)+
     1                   (1.d0-fracx)*fluxz(i,j+1,k+1))

                     endif
                 endif
              else
                 fzp = fluxz(i,j,k+1)
              endif


              divc(i,j,k) = -((apx(i+1,j,k)*fxp - apx(i,j,k)*fxm) / dx
     1                  +(apy(i,j+1,k)*fyp - apy(i,j,k)*fym) / dy /
     1                  +(apz(i,j,k+1)*fzp - apz(i,j,k)*fzm) / dz) /
     2                    vfrac(i,j)


           endif

        enddo
        enddo
              
        do k=1,kmax
        do j=1,jmax
        do i=1,imax
          if(itype(i,j).ne.1)then
C  for noneb cells
            ebdiffop(i,j) = divc(i,j)
          else
            vtot = vfrac(i-1,j-1,k)+vfrac(i-1,j,k)+vfrac(i-1,j+1,k)
     1         + vfrac(i,j-1,k)+vfrac(i,j,k)+vfrac(i,j+1,k)
     2         + vfrac(i+1,j-1,k)+vfrac(i+1,j,k)+vfrac(i+1,j+1,k)
               + vfrac(i-1,j-1,k-1)+vfrac(i-1,j,k-1)+vfrac(i-1,j+1,k-1)
     1         + vfrac(i,j-1,k-1)+vfrac(i,j,k-1)+vfrac(i,j+1,k-1)
     2         + vfrac(i+1,j-1,k-1)+vfrac(i+1,j,k-1)+vfrac(i+1,j+1,k-1)
               + vfrac(i-1,j-1,k+1)+vfrac(i-1,j,k+1)+vfrac(i-1,j+1,k+1)
     1         + vfrac(i,j-1,k+1)+vfrac(i,j,k+1)+vfrac(i,j+1,k+1)
     2         + vfrac(i+1,j-1,k+1)+vfrac(i+1,j,k+1)+vfrac(i+1,j+1,k+1)
            divnc = 0.d0
            do jj = -1,1
            do ii = -1,1
               divnc = divnc + vfrac(i+ii,j+jj)*divc(i+ii,j+jj)
            enddo
            enddo
            divnc = divnc / vtot
            ebdiffop(i,j,k) = vfrac(i,j,k)*divc(i,j,k)+
     1                    (1.d0-vfrac(i,j,k))*divnc
            delm(i,j,k) = vfrac(i,j,k)*(1.d0-vfrac(i,j,k))*
     1                    (divc(i,j,k)-divnc)
          endif
        enddo
        enddo

        do k=1,kmax
        do j=1,jmax
        do i=1,imax

           if(itype(i,j.k) .eq. 1)then
            vtot = vfrac(i-1,j-1,k)+vfrac(i-1,j,k)+vfrac(i-1,j+1,k)
     1         + vfrac(i,j-1,k)+vfrac(i,j+1,k)
     2         + vfrac(i+1,j-1,k)+vfrac(i+1,j,k)+vfrac(i+1,j+1,k)
               + vfrac(i-1,j-1,k-1)+vfrac(i-1,j,k-1)+vfrac(i-1,j+1,k-1)
     1         + vfrac(i,j-1,k-1)+vfrac(i,j,k-1)+vfrac(i,j+1,k-1)
     2         + vfrac(i+1,j-1,k-1)+vfrac(i+1,j,k-1)+vfrac(i+1,j+1,k-1)
               + vfrac(i-1,j-1,k+1)+vfrac(i-1,j,k+1)+vfrac(i-1,j+1,k+1)
     1         + vfrac(i,j-1,k+1)+vfrac(i,j,k+1)+vfrac(i,j+1,k+1)
     2         + vfrac(i+1,j-1,k+1)+vfrac(i+1,j,k+1)+vfrac(i+1,j+1,k+1)
            do kk = -1,1
            do jj = -1,1
            do ii = -1,1
            if((ii.ne. 0 .or. jj.ne.0 .or. kk.ne. 0)
     1           .and.itype(i+ii,j+jj).ne. -1)then
                 ebdiffop(i+ii,j+jj,k+kk) = ebdiffop(i+ii,j+jj,k+kk)+
     1           delm(i,j,k)/vtot
               endif
            enddo
            enddo


           endif
      enddo
      enddo

      return
      end
