!! Knapsack performs a form of load balancing where the 
module knapsack_module

  use bl_types

  implicit none

  integer, private :: knapsack_verbose = 0
  real(kind=dp_t), private :: knapsack_threshold = 1_dp_t

  private :: greater_d

contains

  function greater_d(a,b) result(r)
    logical :: r
    real(kind=dp_t), intent(in) :: a, b
    r = a > b
  end function greater_d

  subroutine knapsack_i(prc, ibxs, np, verbose, threshold)
    use vector_i_module
    use sort_d_module
    use bl_error_module
    !    use named_comparisons_module, only : greater_d
    integer, intent(out), dimension(:) :: prc
    integer, intent(in), dimension(:) :: ibxs
    integer, intent(in) :: np
    real(kind=dp_t), intent(in), optional :: threshold
    integer, intent(in), optional :: verbose

    real(kind=dp_t), dimension(size(ibxs)) :: sizes
    integer, dimension(size(ibxs)) :: isizes
    integer :: i, j, k
    real(kind=dp_t) :: w1o, wjo, w1n, wjn, dif
    type(vector_i), allocatable, dimension(:) :: procs
    integer, allocatable, dimension(:) :: iprocs
    real(kind=dp_t), allocatable, dimension(:) :: pweights
    real(kind=dp_t) :: efficiency, lthresh
    real(kind=dp_t) :: total_weight
    integer :: ierr
    integer :: lverb

    lverb   = knapsack_verbose  ; if ( present(verbose  ) ) lverb   = verbose
    lthresh = knapsack_threshold; if ( present(threshold) ) lthresh = threshold
    if ( np < 1 ) then
       call bl_error('LAYOUT_KNAPSACK: np < 1')
    end if

    allocate(procs(np), iprocs(np), pweights(np), stat=ierr)

    if ( ierr /= 0 ) then
       call bl_error('layout_knapsack: failed to allocate memory')
    end if
    ! each processor maintains a list of its boxes in a
    ! dynamic array, which is addressed indirectly through the
    ! array iprocs
    do i = 1, size(procs)
       call build(procs(i))
       iprocs(i) = i
    end do
    sizes = ibxs
    total_weight = sum(sizes)
    call sort(sizes, isizes, greater_d)
    if ( lverb > 1 ) then
       print *, 'total weight = ', total_weight
       print *, 'sizes = ', sizes
       print *, 'sizes(isizes) =  ', sizes(isizes)
       print *, 'isizes = ', isizes
    end if
    ! place the box in the least loaded processor.
    ! Note: before we start, the procs array satisfies a heap property.
    pweights = weights()
    do i = 1, size(ibxs)
       call push_back(procs(iprocs(1)), isizes(i))
       call reheap_procs()
    end do
    call sort_procs()           ! now, they need to be sorted.
    outer: do
       efficiency = total_weight/np/weight(1)
       if ( efficiency > lthresh ) exit outer
       if ( lverb > 1) then
          do i = 1, size(procs)
             print *, 'dp=', dataptr(procs(i))
             print *, 'ss=', sizes(dataptr(procs(i)))
          end do
          print *, 'weights = ', weights()
       end if
       ! For each box in the most loaded processor
       do i = 1, size(procs(iprocs(1)))
          ! check each less lightly loaded processor... 
          do j = 2, size(procs)
             ! by looking at each box in that processor
             do k = 1, size(procs(iprocs(j)))
                w1o = weight(1)
                wjo = weight(j)
                dif = weight_ball(1,i) - weight_ball(j,k)
                w1n = w1o - dif
                wjn = wjo + dif
                if ( w1n < w1o .AND. wjn < w1o ) then
                   ! if we decrease the overall load in both the processors by
                   ! comparing to the most loaded processor, and we decrease 
                   ! the load in the less loaded processor, move the box.
                   call swap_balls(1,i,j,k)
                   ! reorder the processors.
                   call sort_procs()
                   ! start over again.
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
    prc = prc - 1               ! processor numbers start at 0.
    if ( lverb > 0 ) then
       print *, 'np           = ', np
       print *, 'n            = ', size(ibxs)
       print *, 'max weight   = ', weight(1)
       print *, 'total weight = ', sum(weights())
       print *, 'efficiency   = ', sum(weights())/(np*weight(1))
    end if
    do i = 1, size(procs)
       call destroy(procs(i))
    end do

    !! INTEL FIXME !!
    deallocate(procs, iprocs, pweights)

  contains

    subroutine swap_balls(m,i,j,k)
      integer, intent(in) :: m, i, j, k
      integer dmi, djk, ii, kk

      dmi = ball(m,i)
      ii = erase(procs(iprocs(m)), i)
      djk = ball(j,k)
      kk = erase(procs(iprocs(j)), k)
      call push_back(procs(iprocs(j)), dmi)
      call push_back(procs(iprocs(m)), djk)

    end subroutine swap_balls

    function weights() result(r)
      real(kind=dp_t), dimension(size(procs)) :: r
      integer j

      r = 0
      do j = 1, size(procs)
         r(j) = weight(j)
      end do

    end function weights

    function weight(i) result(r)
      integer, intent(in) :: i
      real(kind=dp_t) :: r
      integer j

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

      iss = weights()
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

end module knapsack_module
