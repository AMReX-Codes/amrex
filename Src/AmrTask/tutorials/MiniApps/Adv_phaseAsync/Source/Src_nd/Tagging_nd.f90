
! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on the state
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: tag        <=  integer tag array
! ::: tag_lo,hi   => index extent of tag array
! ::: state       => state array
! ::: state_lo,hi => index extent of state array
! ::: set         => integer value to tag cell for refinement
! ::: clear       => integer value to untag cell
! ::: lo,hi       => work region we are allowed to change
! ::: dx          => cell size
! ::: problo      => phys loc of lower left corner of prob domain
! ::: time        => problem evolution time
! ::: level       => refinement level of this array
! ::: -----------------------------------------------------------

subroutine state_error(tag,tag_lo,tag_hi, &
                       state,state_lo,state_hi, &
                       set,clear,&
                       lo,hi,&
                       dx,problo,time,level) bind(C, name="state_error")

  use tagging_params_module, only : phierr, phigrad, max_phierr_lev, max_phigrad_lev
  implicit none
  
  integer          :: lo(3),hi(3)
  integer          :: state_lo(3),state_hi(3)
  integer          :: tag_lo(3),tag_hi(3)
  double precision :: state(state_lo(1):state_hi(1), &
                            state_lo(2):state_hi(2), &
                            state_lo(3):state_hi(3))
  integer          :: tag(tag_lo(1):tag_hi(1),tag_lo(2):tag_hi(2),tag_lo(3):tag_hi(3))
  double precision :: problo(3),dx(3),time
  integer          :: level,set,clear

  double precision :: ax, ay, az
  integer          :: i, j, k, dim

  if (state_lo(3) .eq. state_hi(3)) then
     dim = 2
  else
     dim = 3
  end if

  ! Tag on regions of high phi
  if (level .lt. max_phierr_lev) then
     do       k = lo(3), hi(3)
        do    j = lo(2), hi(2)
           do i = lo(1), hi(1)
              if (state(i,j,k) .ge. phierr(level)) then
                 tag(i,j,k) = set
              endif
           enddo
        enddo
     enddo
  endif

  ! Tag on regions of high phisity gradient
  if (level .lt. max_phigrad_lev) then
     do       k = lo(3), hi(3)
        do    j = lo(2), hi(2)
           do i = lo(1), hi(1)
              ax = abs(state(i-1,j,k)-state(i,j,k))
              ax = max(ax, abs(state(i,j,k)-state(i+1,j,k)))
              ay = abs(state(i,j-1,k)-state(i,j,k))
              ay = max(ay, abs(state(i,j,k)-state(i,j+1,k)))
              if (dim .eq. 2) then
                 az = 0.d0
              else
                 az = abs(state(i,j,k-1)-state(i,j,k))
                 az = max(az, abs(state(i,j,k)-state(i,j,k+1)))
              end if
              if (max(ax,ay,az) .ge. phigrad(level)) then
                 tag(i,j,k) = set
              end if
            enddo
         enddo
      end do
   endif
  
end subroutine state_error

