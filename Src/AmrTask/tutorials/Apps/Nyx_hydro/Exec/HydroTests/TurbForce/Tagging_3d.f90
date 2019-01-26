      subroutine tag_denerror(tag,tagl1,tagl2,tagl3,tagh1,tagh2,tagh3, &
                              set,clear, &
                              den,denl1,denl2,denl3,denh1,denh2,denh3, &
                              lo,hi,nd,domlo,domhi, &
                              delta,xlo,problo,time,level)
      use amrex_fort_module, only : rt => amrex_real
      use probdata_module
      use turbforce_module, only : forcing_time_scale_max, stop_forcing
      implicit none

      integer set, clear, nd, level
      integer tagl1,tagl2,tagl3,tagh1,tagh2,tagh3
      integer denl1,denl2,denl3,denh1,denh2,denh3
      integer lo(3), hi(3), domlo(3), domhi(3)
      integer tag(tagl1:tagh1,tagl2:tagh2,tagl3:tagh3)
      real(rt) den(denl1:denh1,denl2:denh2,denl3:denh3,nd)
      real(rt) delta(3), xlo(3), problo(3), time

      real(rt) :: ax,ay,az
      integer          :: n_lo, n_hi
      integer          :: i, j, k, domhalf
      domhalf = (domhi(1)-domlo(1))/2
      n_lo = domlo(1) + domhalf - domhalf/(2**(level+1))
      n_hi = domlo(1) + domhalf + domhalf/(2**(level+1))

      ! use this to start at t = 0.25
      !if (time .ge. 0.125*(1 + 0.5*level)*stop_forcing*forcing_time_scale_max) then
      ! use this to test with a = 1.0
      if (time .ge. (1 + 0.5*level)*stop_forcing*forcing_time_scale_max) then
      ! use this to test with a = 0.5
      !if (time .ge. 0.5*(1 + 0.5*level)*stop_forcing*forcing_time_scale_max) then
         !$OMP parallel do private(k,j,i)
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  if (i.ge.n_lo .and. i.lt.n_hi .and. &
                      j.ge.n_lo .and. j.lt.n_hi .and. &
                      k.ge.n_lo .and. k.lt.n_hi) then
                     tag(i,j,k) = set
                  endif
               enddo
            enddo
         enddo
      endif

!     Tag on regions of high density
!     if (level .lt. max_denerr_lev) then
!        do k = lo(3), hi(3)
!           do j = lo(2), hi(2)
!              do i = lo(1), hi(1)
!                if (den(i,j,k,1) .ge. denerr) then
!                    tag(i,j,k) = set
!                 endif
!              enddo
!           enddo
!        enddo
!     endif

!     Tag on regions of high density gradient
!     if (level .lt. max_dengrad_lev) then
!        do k = lo(3), hi(3)
!           do j = lo(2), hi(2)
!              do i = lo(1), hi(1)
!                 ax = ABS(den(i+1,j,k,1) - den(i,j,k,1))
!                 ay = ABS(den(i,j+1,k,1) - den(i,j,k,1))
!                 az = ABS(den(i,j,k+1,1) - den(i,j,k,1))
!                 ax = MAX(ax,ABS(den(i,j,k,1) - den(i-1,j,k,1)))
!                 ay = MAX(ay,ABS(den(i,j,k,1) - den(i,j-1,k,1)))
!                 az = MAX(az,ABS(den(i,j,k,1) - den(i,j,k-1,1)))
!                 if ( MAX(ax,ay,az) .ge. dengrad) then
!                    tag(i,j,k) = set
!                 endif
!              enddo
!           enddo
!        enddo
!     endif

      end subroutine tag_denerror
