
      subroutine tag_denerror(tag,tagl1,tagl2,tagl3,tagh1,tagh2,tagh3, &
                              set,clear, &
                              den,denl1,denl2,denl3,denh1,denh2,denh3, &
                              lo,hi,nd,domlo,domhi, &
                              delta,xlo,problo,time,level)
      use amrex_fort_module, only : rt => amrex_real
      use probdata_module
      implicit none

      integer set, clear, nd, level
      integer tagl1,tagl2,tagl3,tagh1,tagh2,tagh3
      integer denl1,denl2,denl3,denh1,denh2,denh3
      integer lo(3), hi(3), domlo(3), domhi(3)
      integer tag(tagl1:tagh1,tagl2:tagh2,tagl3:tagh3)
      real(rt) den(denl1:denh1,denl2:denh2,denl3:denh3,nd)
      real(rt) delta(3), xlo(3), problo(3), time

      real(rt) ax,ay,az
      integer i, j, k

!     Tag on regions of high density
      if (level .lt. max_denerr_lev) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  if (den(i,j,k,1) .ge. denerr) then
                     tag(i,j,k) = set
                  endif
               enddo
            enddo
         enddo
      endif

!     Tag on regions of high density gradient
      if (level .lt. max_dengrad_lev) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  ax = ABS(den(i+1,j,k,1) - den(i,j,k,1))
                  ay = ABS(den(i,j+1,k,1) - den(i,j,k,1))
                  az = ABS(den(i,j,k+1,1) - den(i,j,k,1))
                  ax = MAX(ax,ABS(den(i,j,k,1) - den(i-1,j,k,1)))
                  ay = MAX(ay,ABS(den(i,j,k,1) - den(i,j-1,k,1)))
                  az = MAX(az,ABS(den(i,j,k,1) - den(i,j,k-1,1)))
                  if ( MAX(ax,ay,az) .ge. dengrad) then
                     tag(i,j,k) = set
                  endif
               enddo
            enddo
         enddo
      endif

      end subroutine tag_denerror

! ::: -----------------------------------------------------------

      subroutine tag_presserror(tag,tagl1,tagl2,tagl3,tagh1,tagh2,tagh3, &
                                set,clear, &
                                press,pressl1,pressl2,pressl3,pressh1, &
                                pressh2,pressh3, &
                                lo,hi,np,domlo,domhi, &
                                delta,xlo,problo,time,level)

      use probdata_module
      implicit none

      integer set, clear, np, level
      integer tagl1,tagl2,tagl3,tagh1,tagh2,tagh3
      integer pressl1,pressl2,pressl3,pressh1,pressh2,pressh3
      integer lo(3), hi(3), domlo(3), domhi(3)
      integer tag(tagl1:tagh1,tagl2:tagh2,tagl3:tagh3)
      real(rt) press(pressl1:pressh1,pressl2:pressh2, &
                             pressl3:pressh3,np)
      real(rt) delta(3), xlo(3), problo(3), time

      real(rt) ax,ay,az
      integer i, j, k

!     Tag on regions of high pressure
      if (level .lt. max_presserr_lev) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  if (press(i,j,k,1) .ge. presserr) then
                     tag(i,j,k) = set
                  endif
               enddo
            enddo
         enddo
      endif

!     Tag on regions of high pressure gradient
      if (level .lt. max_pressgrad_lev) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  ax = ABS(press(i+1,j,k,1) - press(i,j,k,1))
                  ay = ABS(press(i,j+1,k,1) - press(i,j,k,1))
                  az = ABS(press(i,j,k+1,1) - press(i,j,k,1))
                  ax = MAX(ax,ABS(press(i,j,k,1) - press(i-1,j,k,1)))
                  ay = MAX(ay,ABS(press(i,j,k,1) - press(i,j-1,k,1)))
                  az = MAX(az,ABS(press(i,j,k,1) - press(i,j,k-1,1)))
                  if ( MAX(ax,ay,az) .ge. pressgrad) then
                     tag(i,j,k) = set
                  endif
               enddo
            enddo
         enddo
      endif

      end subroutine tag_presserror

! ::: -----------------------------------------------------------

      subroutine tag_velerror(tag,tagl1,tagl2,tagl3,tagh1,tagh2,tagh3, &
                              set,clear, &
                              vel,vell1,vell2,vell3,velh1,velh2,velh3, &
                              lo,hi,nv,domlo,domhi, &
                              delta,xlo,problo,time,level)

      use probdata_module
      implicit none

      integer set, clear, nv, level
      integer tagl1,tagl2,tagl3,tagh1,tagh2,tagh3
      integer vell1,vell2,vell3,velh1,velh2,velh3
      integer lo(3), hi(3), domlo(3), domhi(3)
      integer tag(tagl1:tagh1,tagl2:tagh2,tagl3:tagh3)
      real(rt) vel(vell1:velh1,vell2:velh2,vell3:velh3,nv)
      real(rt) delta(3), xlo(3), problo(3), time

      real(rt) ax,ay,az
      integer i, j, k

!     Tag on regions of high velocity gradient
      if (level .lt. max_velgrad_lev) then
         !$OMP PARALLEL DO PRIVATE(i,j,k,ax,ay,az)
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  ax = ABS(vel(i+1,j,k,1) - vel(i,j,k,1))
                  ay = ABS(vel(i,j+1,k,1) - vel(i,j,k,1))
                  az = ABS(vel(i,j,k+1,1) - vel(i,j,k,1))
                  ax = MAX(ax,ABS(vel(i,j,k,1) - vel(i-1,j,k,1)))
                  ay = MAX(ay,ABS(vel(i,j,k,1) - vel(i,j-1,k,1)))
                  az = MAX(az,ABS(vel(i,j,k,1) - vel(i,j,k-1,1)))
                  if ( MAX(ax,ay,az) .ge. velgrad) then
                     tag(i,j,k) = set
                  endif
               enddo
            enddo
         enddo
         !$OMP END PARALLEL DO
      endif

      end subroutine tag_velerror
