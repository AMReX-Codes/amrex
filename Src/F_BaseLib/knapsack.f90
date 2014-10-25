
module sfc_module
  use box_module
  !
  ! This is more than a bit of a hack.  We're going to store a few values here
  ! while inside the sfc_i() routine, so that they're accessible via
  ! sfc_greater_i().
  !
  integer              :: dm, mpower
  integer, allocatable :: boxkeys(:,:)
  type(box), pointer   :: pbxs(:)
end module sfc_module

module knapsack_module

  use bl_types

  implicit none

  logical, private :: knapsack_verbose = .false.

  real(kind=dp_t), private :: knapsack_threshold = 0.9_dp_t

  private :: greater_d, sfc_greater_i

contains

  subroutine knapsack_set_verbose(yesorno)
    logical, intent(in) :: yesorno
    knapsack_verbose = yesorno
  end subroutine knapsack_set_verbose

  pure function greater_d(a,b) result(r)
    logical :: r
    real(kind=dp_t), intent(in) :: a, b
    r = a > b
  end function greater_d

  subroutine knapsack_i(prc, ibxs, np, verbose, threshold)

    use parallel
    use box_module
    use vector_i_module
    use sort_d_module
    use bl_error_module
    use bl_prof_module

    integer,         intent(out), dimension(:) :: prc
    integer,         intent(in ), dimension(:) :: ibxs
    integer,         intent(in )               :: np
    real(kind=dp_t), intent(in ), optional     :: threshold
    logical,         intent(in ), optional     :: verbose

    logical         :: lverb
    integer         :: i, j, k, isizes(size(ibxs))
    real(kind=dp_t) :: w1o, wjo, w1n, wjn, dif, efficiency, lthresh, t1, t2
    real(kind=dp_t) :: total_weight, xmean, stddev, sizes(size(ibxs))

    integer,         allocatable :: iprocs(:)
    type(vector_i),  allocatable :: procs(:)
    real(kind=dp_t), allocatable :: pweights(:)

    type(bl_prof_timer), save :: bpt

    if ( np < 1 ) call bl_error('knapsack_i(): np < 1')

    call build(bpt, 'knapsack')

    call cpu_time(t1)

    lverb   = knapsack_verbose  ; if ( present(verbose  ) ) lverb   = verbose
    lthresh = knapsack_threshold; if ( present(threshold) ) lthresh = threshold

    allocate(procs(np), iprocs(np), pweights(np))
    !
    ! Each processor maintains a list of its boxes in a dynamic array,
    ! which is addressed indirectly through the array iprocs.
    !
    do i = 1, size(procs)
       call build(procs(i))
       iprocs(i) = i
    end do
    sizes = ibxs
    total_weight = sum(sizes)
    call sort(sizes, isizes, greater_d)
    !
    ! Place the box in the least loaded processor.
    ! Before we start, the procs array satisfies a heap property.
    !
    pweights = weights(np)
    do i = 1, size(ibxs)
       call push_back(procs(iprocs(1)), isizes(i))
       call reheap_procs()
    end do
    call sort_procs() ! now, they need to be sorted.
    outer: do
       efficiency = total_weight/np/weight(1)
       if ( efficiency > lthresh ) exit outer
       !
       ! For each box in the most loaded processor.
       !
       do i = 1, size(procs(iprocs(1)))
          !
          ! Check each less lightly loaded processor...
          !
          do j = 2, size(procs)
             !
             ! By looking at each box in that processor.
             !
             do k = 1, size(procs(iprocs(j)))
                w1o = weight(1)
                wjo = weight(j)
                dif = weight_ball(1,i) - weight_ball(j,k)
                w1n = w1o - dif
                wjn = wjo + dif
                if ( w1n < w1o .AND. wjn < w1o ) then
                   !
                   ! If we decrease the overall load in both the processors by
                   ! comparing to the most loaded processor, and we decrease 
                   ! the load in the less loaded processor, move the box.
                   !
                   call swap_balls(1,i,j,k)
                   !
                   ! Reorder the processors.
                   !
                   call sort_procs()
                   !
                   ! Start over again.
                   !
                   cycle outer
                end if
             end do
          end do
       end do
       exit outer
    end do outer

    do i = 1, size(procs)
       do j = 1, size(procs(iprocs(i)))
          prc(ball(i,j)) = i
       end do
    end do

    prc = prc - 1   ! Processor numbers start at 0.

    call cpu_time(t2)

     if ( .false. .and. lverb .and. np > 1 .and. parallel_ioprocessor() ) then
        xmean  = sum(ibxs)/size(ibxs)
        stddev = sqrt(sum((ibxs-xmean)**2)/(size(ibxs)-1))
        print *, 'np               = ', np
        print *, 'n                = ', size(ibxs)
        print *, 'max box weight   = ', maxval(ibxs)
        print *, 'min box weight   = ', minval(ibxs)
        print *, 'mean bx weight   = ', xmean
        print *, 'stdd bx weight   = ', stddev
        print *, 'max weight       = ', weight(1)
        print *, 'max weight       = ', weight(np)
        print *, 'total weight     = ', sum(weights(np))
     end if

    efficiency = sum(weights(np))/(np*weight(1))

    do i = 1, size(procs)
       call destroy(procs(i))
    end do

    call cpu_time(t2)

    if ( lverb .and. np > 1 .and. parallel_ioprocessor() ) then
       print *, 'KNAPSACK effi = ', efficiency
       print *, 'KNAPSACK time = ', t2-t1
    end if

    call destroy(bpt)

  contains

    subroutine swap_balls(m,i,j,k)
      integer, intent(in) :: m, i, j, k
      integer :: dmi, djk, ii, kk
      dmi = ball(m,i)
      ii = erase(procs(iprocs(m)), i)
      djk = ball(j,k)
      kk = erase(procs(iprocs(j)), k)
      call push_back(procs(iprocs(j)), dmi)
      call push_back(procs(iprocs(m)), djk)
    end subroutine swap_balls

    function weights(np) result(r)
      integer, intent(in) :: np
      real(kind=dp_t), dimension(np) :: r
      integer :: j
      r = 0
      do j = 1, np
         r(j) = weight(j)
      end do
    end function weights

    function weight(i) result(r)
      integer, intent(in) :: i
      real(kind=dp_t) :: r
      integer :: j
      r = 0
      do j = 1, size(procs(iprocs(i)))
         r = r + sizes(at(procs(iprocs(i)),j))
      end do
    end function weight

    function weight_ball(i,j) result(r)
      integer, intent(in) :: i, j
      real(kind=dp_t) :: r
      r = sizes(ball(i,j))
    end function weight_ball

    function ball(i,j) result(r)
      integer, intent(in) :: i, j
      integer :: r
      r = at(procs(iprocs(i)),j)
    end function ball

    subroutine sort_procs()
      integer :: iii(size(procs))
      real(kind=dp_t) :: iss(size(procs))
      iss = weights(np)
      call sort(iss, iii, greater_d)
      iprocs = iprocs(iii)
    end subroutine sort_procs

    subroutine reheap_procs()
      integer :: i, dd, sz, tmp
      real(kind=dp_t) :: wt
      sz = size(procs)
      wt = weight(1)
      dd = iprocs(1)
      ! first remove the last
      tmp = iprocs(sz)
      sz = sz - 1
      i = 2
      do while ( i <= sz+1/2 )
         if ( i < sz .AND. &
              (pweights(iprocs(i)) > pweights(iprocs(i+1))) ) i = i + 1
         if ( .not. ( pweights(tmp) >= pweights(iprocs(i))) ) exit
         iprocs(i/2) = iprocs(i)
         i = 2*i
      end do
      iprocs(i/2) = tmp
      ! now, put it back!
      sz = sz + 1
      i = sz
      iprocs(i) = dd
      pweights(dd) = wt
      do while ( i > 1 )
         if ( .not. (pweights(iprocs(i/2)) > pweights(dd)) ) exit
         iprocs(i) = iprocs(i/2)
         i = i/2
      end do
      iprocs(i) = dd
    end subroutine reheap_procs

  end subroutine knapsack_i

  subroutine make_sfc(iorder, bxs)
    use sfc_module
    use box_module
    use sort_i_module
    use bl_error_module
    use bl_prof_module
    
    integer,   intent(out)           :: iorder(:)
    type(box), intent(in ), target   :: bxs(:)
    
    integer :: i
    type(bl_prof_timer), save :: bpt

    call build(bpt, 'sfc')
    
    call init_sfc_module()

    do i = 1, size(bxs)
       iorder(i) = i
    end do

    call mergesort_i(iorder, sfc_greater_i)
    
    call close_sfc_module()

    call destroy(bpt)

  contains

    subroutine init_sfc_module
      integer :: i
      dm     =  get_dim(bxs(1))
      allocate(boxkeys(dm,size(bxs)))
      do i = 1, size(bxs)
         boxkeys(:,i) = make_box_key(i)
      end do
      mpower =  maxpower()
      pbxs   => bxs
    end subroutine init_sfc_module

    subroutine close_sfc_module()
      deallocate(boxkeys)
      nullify(pbxs)
    end subroutine close_sfc_module

    function make_box_key(i) result(r)
      integer :: r(dm)
      integer, intent(in) :: i
!      r = lwb(bxs(i))
      r = (lwb(bxs(i)) + upb(bxs(i))) / 2
    end function make_box_key

    function maxpower() result(r)
      
      integer :: r, maxijk, i
      
      maxijk = -Huge(1)
      
      do i = 1,size(bxs)
         maxijk = max(maxijk, maxval(boxkeys(:,i)))
      end do
      
      r = 0
      do while ( ishft(1,r) <= maxijk )
         r = r + 1
      end do
      
    end function maxpower
    
  end subroutine make_sfc

  subroutine distribute_sfc(prc, iorder, ibxs, np, verbose)

    use parallel

    integer,   intent(out)           :: prc(:)
    integer,   intent(in )           :: iorder(:)
    integer,   intent(in )           :: ibxs(:)
    integer,   intent(in )           :: np
    logical,   intent(in ), optional :: verbose

    integer :: i, k, cnt, sz, nbxs
    integer, allocatable :: whichcpu(:)
    real(kind=dp_t) :: totalvol, volpercpu, vol, maxvol, efficiency
    logical :: lverb

    lverb = knapsack_verbose ; if ( present(verbose) ) lverb = verbose

    nbxs = size(prc)

    allocate(whichcpu(nbxs))

    maxvol = -Huge(1.0_dp_t)
    
    volpercpu = 0.0_dp_t
    do i = 1, size(ibxs)
       volpercpu = volpercpu + ibxs(i)
    end do
    volpercpu = volpercpu / np
    
    k        = 1
    sz       = size(ibxs)
    totalvol = 0.0_dp_t

    do i = 1, np
       
       cnt = 0
       vol = 0.0_dp_t
       
       do while ( (k <= sz) .and. ((i == np) .or. (vol < volpercpu)) )
          whichcpu(k) = i
          vol = vol + ibxs(iorder(k))
          k   = k   + 1
          cnt = cnt + 1
       end do
       
       totalvol = totalvol + vol
       
       if ( (totalvol/i) > volpercpu .and. (cnt > 1) .and. (k <= sz) ) then
          k        = k - 1
          vol      = vol - ibxs(iorder(k))
          totalvol = totalvol - ibxs(iorder(k))
       endif
       
       maxvol = max(maxvol,vol)
       
    end do
    !
    ! Force "whichcpu" values to be zero-based instead of one-based.
    !
    whichcpu = whichcpu - 1

    do i = 1, nbxs
       prc(iorder(i)) = whichcpu(i)
    end do
    
    efficiency = volpercpu / maxvol

    if ( lverb .and. np > 1 .and. parallel_ioprocessor() ) then
       print *, 'SFC effi = ', efficiency
    end if
    
  end subroutine distribute_sfc

  function sfc_greater_i(ilhs,irhs) result(r)

    use sfc_module

    logical             :: r
    integer, intent(in) :: ilhs,irhs

    integer m,d,NNN,lkey,rkey

    do m = mpower-1,0,-1

       NNN = ishft(1,m)

       do d = 1, dm

          lkey = boxkeys(d,ilhs) / NNN
          rkey = boxkeys(d,irhs) / NNN

          if ( lkey < rkey ) then
             r = .true.
             return
          else if ( lkey > rkey ) then
             r = .false.
             return
          end if
       end do
    end do

    r = .false.

  end function sfc_greater_i

end module knapsack_module
