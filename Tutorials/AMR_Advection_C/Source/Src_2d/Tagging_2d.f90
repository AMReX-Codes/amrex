
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
      subroutine state_error(tag,tagl1,tagl2,tagh1,tagh2, &
                             set,clear, &
                             spec,specl1,specl2,spech1,spech2, &
                             lo,hi,nd,domlo,domhi, &
                             delta,xlo,problo,time,level)
      use probdata_module
      implicit none

      integer set, clear, nd, level
      integer tagl1,tagl2,tagh1,tagh2
      integer specl1,specl2,spech1,spech2
      integer lo(2), hi(2), domlo(2), domhi(2)
      integer tag(tagl1:tagh1,tagl2:tagh2)
      double precision spec(specl1:spech1,specl2:spech2,nd)
      double precision delta(2), xlo(2), problo(2), time

      double precision ax,ay
      integer i, j

!     Tag on regions of high concentration
!     In this example we only test on whether the 2nd species has conc > specerr
      if (level .lt. max_specerr_lev) then
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               if (spec(i,j,2) .ge. specerr) then
                  tag(i,j) = set
               endif
            enddo
         enddo
      endif

!     Tag on regions of high gradient
      if (level .lt. max_specgrad_lev) then
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               ax = ABS(spec(i+1,j,1) - spec(i,j,1))
               ay = ABS(spec(i,j+1,1) - spec(i,j,1))
               ax = MAX(ax,ABS(spec(i,j,1) - spec(i-1,j,1)))
               ay = MAX(ay,ABS(spec(i,j,1) - spec(i,j-1,1)))
               if ( MAX(ax,ay) .ge. specgrad) then
                  tag(i,j) = set
               endif
            enddo
         enddo
      endif
      
      end subroutine state_error
