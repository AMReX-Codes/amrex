
subroutine t_smoother
  use mg_smoother_module
  implicit none
  integer  n, nn
  real(kind=dp_t),allocatable :: ss(:,:,:)
  real(kind=dp_t),allocatable :: uu(:,:)
  real(kind=dp_t),allocatable :: ff(:,:)
  integer        ,allocatable :: mm(:,:)
  logical :: skwd
  integer :: ng, j, i
  integer :: lo(2), ll, hh
  real(kind=dp_t) :: omega

  n = 8
  allocate(ss(2*n,2*n,0:4))
  allocate(uu(0:2*n+1,0:2*n+1))
  allocate(ff(2*n,2*n))
  allocate(mm(n,n))

  ss(:,:,0) = -4
  ss(:,:,1) =  1
  ss(:,:,2) =  1
  ss(:,:,3) =  1
  ss(:,:,4) =  1

  uu        =  0
  ff        =  0
  ll = n/2-1
  hh = n/2
  ff(ll:hh    ,ll:hh    ) = 1
  ff(ll+n:hh+n,ll:hh    ) = 1
  ff(ll:hh    ,ll+n:hh+n) = 1
  ff(ll+n:hh+n,ll+n:hh+n) = 1
  omega = 1
  lo = 1
  ng = 1
  skwd = .false.
  
  do nn = 0, 1
     call gs_rb_smoother_2d(omega, ss, uu(0:n+1  ,0:n+1  ), ff(1:n    ,1:n    ), mm, &
          lo + (/0,0/), ng, nn, skwd)
     call gs_rb_smoother_2d(omega, ss, uu(n:2*n+1,0:n-1  ), ff(n+1:2*n,1:n    ), mm, &
          lo + (/n,0/), ng, nn, skwd)
     call gs_rb_smoother_2d(omega, ss, uu(0:n+1  ,n:2*n+1), ff(1:n    ,n+1:2*n), mm, &
          lo + (/0,n/), ng, nn, skwd)
     call gs_rb_smoother_2d(omega, ss, uu(n:2*n+1,n:2*n+1), ff(n+1:2*n,n+1:2*n), mm, &
          lo + (/n,n/), ng, nn, skwd)
  end do

  do j = 1, n
     do i = 1, n
        write(1,fmt=*) i,j,uu(i,j)
     end do
  end do

  uu = 0
  do nn = 0, 1
     call gs_rb_smoother_2d(omega, ss, uu, ff, mm, lo, ng, nn, skwd)
  end do
  do j = 1, n
     do i = 1, n
        write(2,fmt=*) i,j,uu(i,j)
     end do
  end do
end subroutine t_smoother
