subroutine t_knapsack
  use bl_IO_module
  use f2kcli
  use knapsack_module
  integer :: un, np, n, idm
  real(kind=dp_t) :: t1, t2
  real(kind=dp_t) :: thresh
  integer :: verbose
  integer, allocatable :: iweights(:), prc(:)
  real(kind=dp_t) :: maxprc, minprc, xmean, stddev
  real(kind=dp_t), allocatable :: weights(:)
  integer :: i
  character(len=128) fname

  if ( command_argument_count() < 1 ) then
     np = 128
  else
     call get_command_argument(1, value = fname)
     read(fname,*) np
  end if

  un = unit_new()
  open(un, file="conn_defs", status = 'old', action = 'read')
  read(un,fmt=*) n
  allocate(iweights(n))
  do i = 1, n
     read(un,fmt=*) idm, iweights(i)
  end do
  close(un)
  allocate(prc(n))
  verbose = 0; thresh = 1.0_dp_t
  call cpu_time(t1)
  call knapsack_i(prc, iweights, np, verbose, thresh)
  call cpu_time(t2)
  allocate(weights(n))
  weights = iweights
  maxprc = -Huge(maxprc)
  minprc =  Huge(minprc)
  do i = 0, np-1
     maxprc = max(maxprc, sum(weights, mask = prc==i))
     minprc = min(minprc, sum(weights, mask = prc==i))
  end do
  print *, 'np               = ', np
  print *, 'n                = ', size(iweights)
  print *, 'max box weight   = ', maxval(weights)
  print *, 'min box weight   = ', minval(weights)
  xmean = sum(weights)/size(weights)
  print *, 'mean bx weight   = ', xmean
  stddev = sqrt(sum((weights-xmean)**2)/(size(weights)-1))
  print *, 'stdd bx weight   = ', stddev
  print *, 'max prc weight   = ', maxprc
  print *, 'min prc weight   = ', minprc
  print *, 'weight prc(0)    = ', sum(weights, mask = prc==0)
  print *, 'total weight     = ', sum(weights)
  print *, 'efficiency       = ', sum(weights)/np/sum(weights,mask = prc==0)
  print *, 'knapsack time    = ', t2-t1

  un = unit_new()
  open(un,file="knapsack.out", status='replace', action = 'write')
  write(un,fmt='(i5)') n
  do i = 1, n
     write(un,fmt='(i5,1x,i5,1x,i10)') i, prc(i), iweights(i)
  end do
  close(un)

end subroutine t_knapsack
