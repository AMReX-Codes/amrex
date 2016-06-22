program fdata
  implicit none

  call flamelength

end program fdata

! Estimates the length of the flame.
! Note that we don't handle periodic boundaries yet, so we
! restrict the domain to accrete(domain,-1) of the original domain.
subroutine flamelength
  use plotfile_module
  use filler_module
  use bl_IO_module
  implicit none
  integer :: i, j, n, ii, jj, f, nc
  real(kind=dp_t), allocatable :: c_fab(:,:,:)
  real(kind=dp_t) :: dt, dy, y1, t1
  integer :: lo(2), hi(2), plo(2), phi(2), comp
  type(plotfile) pf
  integer :: unit
  integer :: level
  character(len=64) :: fname
  integer :: narg, farg

  integer :: nx

  real(kind=dp_t) :: dx

  integer, dimension(2) :: tlo(2), thi(2)

  logical :: found

  integer :: divide
  integer :: cnt
  real(kind=dp_t) :: cut
  real(kind=dp_t) :: fblo2m1, fblo1m1, fbhi2p1, fbhi1p1
  integer :: max_level
  
  narg = command_argument_count()
  farg = 1
  cut = .25_dp_t
  divide = 1
  max_level = -huge(max_level)

  do while ( farg <= narg )
     call get_command_argument(farg, value = fname)
     select case (fname)

     case ('-c', '--cut')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname, *) cut

     case ('--max_level')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname,*) max_level

     case ('--divide')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname,*) divide

     case default
        exit
     end select
     farg = farg + 1
  end do

  if ( farg > narg ) return

  unit = unit_new()

  call get_command_argument(farg, value = fname)
  call build(pf, fname, unit)
  level = max(pf%flevel, max_level)
  plo = lwb(plotfile_get_pd_box(pf, level))
  phi = upb(plotfile_get_pd_box(pf, level))

! figure out which component is the carbon mass fraction
  nc = 1
  
  found = .false.
  do i = 1, pf%nvars
     if (pf%names(i) == "Y(C12)") then
        comp = i
        found = .true.
        exit
     endif
  enddo

  if (.NOT. found) then
     print *, 'ERROR: carbon mass fraction not found'
     stop
  endif

  nx = phi(1) - plo(1) + 1
  dx = nx/divide


  do f = farg, narg
     cnt = 0
     call get_command_argument(f, value = fname)
     if ( f /= farg ) call build(pf, fname, unit)


     ! we are going to blow this out to a fab in pieces, so we can
     ! handle large domains.  We cut in the x-direction
     do n = 1, divide

        ! allocate the storage for the current subdomain
        tlo(1) = (n-1)*dx
        thi(1) = n*dx -1

        if (n == divide) then
           thi(1) = nx - 1
        endif

        if (n > 1) then
           tlo(1) = tlo(1) - 1
        endif
        
        if (n < divide) then
           thi(1) = thi(1) + 1
        endif

        tlo(2) = plo(2)
        thi(2) = phi(2)

        allocate(c_fab(tlo(1):thi(1), tlo(2):thi(2),nc))
        call blow_out_to_sub_fab(c_fab, tlo, thi, pf, (/comp/), level)

        
        do jj = tlo(2)+1, thi(2)-1
           do ii = tlo(1)+1, thi(1)-1
              fblo2m1 = c_fab(ii,jj-1,1)
              fbhi2p1 = c_fab(ii,jj+1,1)
              fblo1m1 = c_fab(ii-1,jj,1)
              fbhi1p1 = c_fab(ii+1,jj,1)
              if ( .TRUE. ) then
                 if ( c_fab(ii,jj,1) > CUT .and. fblo2m1 < CUT ) cnt = cnt + 1
                 if ( c_fab(ii,jj,1) < CUT .and. fblo2m1 > CUT ) cnt = cnt + 1
                 if ( c_fab(ii,jj,1) > CUT .and. fbhi2p1 < CUT ) cnt = cnt + 1
                 if ( c_fab(ii,jj,1) < CUT .and. fbhi2p1 > CUT ) cnt = cnt + 1
              
                 if ( c_fab(ii,jj,1) > CUT .and. fblo1m1 < CUT ) cnt = cnt + 1
                 if ( c_fab(ii,jj,1) < CUT .and. fblo1m1 > CUT ) cnt = cnt + 1
                 if ( c_fab(ii,jj,1) > CUT .and. fbhi1p1 < CUT ) cnt = cnt + 1
                 if ( c_fab(ii,jj,1) < CUT .and. fbhi1p1 > CUT ) cnt = cnt + 1
              else
                 fblo2m1 = fblo2m1 - CUT
                 fbhi2p1 = fbhi2p1 - CUT
                 fblo1m1 = fblo1m1 - CUT
                 fbhi1p1 = fbhi1p1 - CUT
                 if ( (c_fab(ii,jj,1)-CUT)*fblo2m1 < 0 ) cnt = cnt + 1
                 if ( (c_fab(ii,jj,1)-CUT)*fbhi2p1 < 0 ) cnt = cnt + 1
                 if ( (c_fab(ii,jj,1)-CUT)*fblo1m1 < 0 ) cnt = cnt + 1
                 if ( (c_fab(ii,jj,1)-CUT)*fbhi1p1 < 0 ) cnt = cnt + 1
              end if
           end do
        end do

        deallocate (c_fab)

     enddo

     ! divide by 2 to avoid double counting
     print '(3(g20.10,1x))', pf%tm, cnt/2
     call destroy(pf)
  end do

end subroutine flamelength
