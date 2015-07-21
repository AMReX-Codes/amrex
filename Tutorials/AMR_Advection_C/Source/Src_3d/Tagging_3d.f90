! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on the density
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: tag      <=  integer tag array
! ::: lo,hi     => index extent of tag array
! ::: set       => integer value to tag cell for refinement
! ::: clear     => integer value to untag cell
! ::: spec      => species array
! ::: nd        => number of components in den array (should be 1)
! ::: domlo,hi  => index extent of problem domain
! ::: delta     => cell spacing
! ::: xlo       => physical location of lower left hand
! :::              corner of tag array
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: level     => refinement level of this array
! ::: -----------------------------------------------------------
      subroutine state_error(tag,tagl1,tagl2,tagl3,tagh1,tagh2,tagh3, &
                             set,clear, &
                             var,varl1,varl2,varl3,varh1,varh2,varh3, &
                             lo,hi,nd,domlo,domhi, &
                             delta,xlo,problo,time,level)
      use probdata_module
      implicit none

      integer          :: set, clear, nd, level
      integer          :: tagl1,tagl2,tagl3,tagh1,tagh2,tagh3
      integer          :: varl1,varl2,varl3,varh1,varh2,varh3
      integer          :: lo(3), hi(3), domlo(3), domhi(3)
      integer          :: tag(tagl1:tagh1,tagl2:tagh2,tagl3:tagh3)
      double precision :: var(varl1:varh1,varl2:varh2,varl3:varh3,nd)
      double precision :: delta(3), xlo(3), problo(3), time
      double precision :: ax,ay,az
      integer          :: i,j,k

!     Tag on regions of high concentration
!     In this example we only test on whether the 2nd species has conc > specerr
      if (level .lt. max_specerr_lev) then
         do k = lo(3), hi(3)
         do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            if (var(i,j,k,2) .ge. specerr) then
               tag(i,j,k) = set
            endif
         enddo
         enddo
         enddo
      endif

!     Tag on regions of high gradient
      if (level .lt. max_specgrad_lev) then
         do k = lo(3), hi(3)
         do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            ax = abs(var(i+1,j,k,1) - var(i,j,k,1))
            ay = abs(var(i,j+1,k,1) - var(i,j,k,1))
            az = abs(var(i,j,k+1,1) - var(i,j,k,1))
            ax = max(ax,abs(var(i,j,k,1) - var(i-1,j,k,1)))
            ay = max(ay,abs(var(i,j,k,1) - var(i,j-1,k,1)))
            az = max(az,abs(var(i,j,k,1) - var(i,j,k-1,1)))
            if ( max(max(ax,ay),az) .ge. specgrad) then
               tag(i,j,k) = set
            endif
            enddo
         enddo
         enddo
      endif
      
      end subroutine state_error
