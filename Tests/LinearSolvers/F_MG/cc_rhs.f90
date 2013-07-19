module cc_rhs_module 
    use BoxLib
    use ml_layout_module
    use multifab_module
    use mt19937_module

    implicit none

contains

  subroutine cc_rhs(mla, pd, rh, rhs_type)

    type(ml_layout), intent(inout) :: mla
    type(box      ), intent(in   ) :: pd
    type( multifab), intent(inout) :: rh(:)
    integer        , intent(in   ) :: rhs_type

    integer                        :: n, dm, nlevs

    dm = mla%dim

    nlevs = mla%nlevel

    do n = nlevs, 1, -1
       call setval(rh(n), val = 0.d0, all=.true.)
    end do

    if (rhs_type .eq. 1) then
       call mf_init_1(rh(nlevs),pd)
    else if (rhs_type .eq. 2) then
       call mf_init_2(rh(nlevs),pd)
    else if (rhs_type .eq. 3) then
       call mf_init_3(rh(nlevs),pd)
    else if (rhs_type .eq. 4) then
       call mf_init_sins(rh(nlevs),pd)
    else if (rhs_type .eq. 5) then
       call mf_init_exact(rh(nlevs),pd)
    else if (rhs_type .eq. 6) then
       call mf_init_rand(rh(nlevs),pd)
    end if

  end subroutine cc_rhs

  subroutine mf_init_1(mf,pd)

    type(multifab), intent(inout) :: mf
    type(box)     , intent(in   ) :: pd
    type(box) :: bx

    bx = pd
    bx%lo(1:bx%dim) = (bx%hi(1:bx%dim) + bx%lo(1:bx%dim))/4
    bx%hi(1:bx%dim) = bx%lo(1:bx%dim)
    call setval(mf, 1.d0, bx)

    call print(bx,'Setting to 1 here ')

    bx = pd
    bx%lo(1:bx%dim) = 3*(bx%hi(1:bx%dim) + bx%lo(1:bx%dim))/4
    bx%hi(1:bx%dim) = bx%lo(1:bx%dim)
    call setval(mf, -1.d0, bx)

    call print(bx,'Setting to -1 here ')

  end subroutine mf_init_1

  subroutine mf_init_2(mf,pd)

    type(multifab), intent(inout) :: mf
    type(box     ), intent(in   ) :: pd

    integer   :: i
    type(box) :: bx

    do i = 1, nfabs(mf)

       bx = get_box(mf,i)
       bx%lo(1:bx%dim) = (bx%hi(1:bx%dim) + bx%lo(1:bx%dim))/2
       bx%hi(1:bx%dim) = bx%lo(1:bx%dim)
       call setval(mf%fbs(i), 1.0_dp_t, bx)

!      Single point of non-zero RHS: use this to make system solvable
       bx = get_box(mf,i)
       bx%lo(1       ) = (bx%hi(1       ) + bx%lo(1       ))/2 + 1
       bx%lo(2:bx%dim) = (bx%hi(2:bx%dim) + bx%lo(2:bx%dim))/2
       bx%hi(1:bx%dim) = bx%lo(1:bx%dim)
       call setval(mf%fbs(i), -1.0_dp_t, bx)

!      1-d Strip: Variation in x-direction
!      bx%lo(1) = (bx%hi(1) + bx%lo(1))/2
!      bx%hi(1) = bx%lo(1)+1

!      1-d Strip: Variation in y-direction
!      bx%lo(2) = (bx%hi(2) + bx%lo(2))/2
!      bx%hi(2) = bx%lo(2)+1

!      1-d Strip: Variation in z-direction
!      bx%lo(3) = (bx%hi(3) + bx%lo(3))/2
!      bx%hi(3) = bx%lo(3)+1

    end do
  end subroutine mf_init_2

  subroutine mf_init_3(mf,pd)

    type(multifab), intent(inout) :: mf
    type(box     )  , intent(in   ) :: pd

    integer   :: i
    type(box) :: bx, rhs_box, rhs_intersect_box

    rhs_box%dim = mf%dim
    rhs_box%lo(1:rhs_box%dim) = 7
    rhs_box%hi(1:rhs_box%dim) = 8

    do i = 1, nfabs(mf)
       bx = get_ibox(mf,i)
       rhs_intersect_box = box_intersection(bx,rhs_box)
       if (.not. empty(rhs_intersect_box)) then
         bx%lo(1:bx%dim) = lwb(rhs_intersect_box)
         bx%hi(1:bx%dim) = upb(rhs_intersect_box)
!         print *,'SETTING RHS IN BOX ',i,' : ', bx%lo(1:bx%dim),bx%hi(1:bx%dim)
         call setval(mf%fbs(i), 1.d0, bx)
       end if
    end do

  end subroutine mf_init_3

  subroutine mf_init_sins(mf,pd)

    type(multifab)  , intent(inout) :: mf
    type(box     )  , intent(in   ) :: pd

    type(box)                :: bx
    real(kind=dp_t)          :: dx
    real(kind=dp_t), pointer :: rp(:,:,:,:)

    integer   :: i,dm,nx,ny

    nx = pd%hi(1) - pd%lo(1) + 1
    ny = pd%hi(2) - pd%lo(2) + 1

    if (nx .ne. ny) then
       print *,'Not sure what to do with nx .neq. ny in mf_init_sins'
       stop
    end if

    dx = 1.d0 / dble(nx)
 
    dm = mf%dim

    print *,'Setting rhs to a sum of sins '
    do i = 1, nfabs(mf)
       bx = get_ibox(mf, i)
       rp => dataptr(mf,i,bx)
       if (dm.eq.2) then
          call init_sin_2d(rp(:,:,1,1),bx,dx)
       else if (dm.eq.3) then 
          call init_sin_3d(rp(:,:,:,1),bx,dx)
       end if
    end do

  end subroutine mf_init_sins

  subroutine init_sin_2d(rhs,bx,dx)

      type(box)       , intent(in   ) :: bx
      double precision, intent(inout) :: rhs(bx%lo(1):,bx%lo(2):)
      double precision, intent(in   ) :: dx

      integer :: i,j
      double precision :: tpi, epi, spi
      double precision :: x, y

      tpi = 8.d0 * datan(1.d0)

      do j = bx%lo(2),bx%hi(2)
      do i = bx%lo(1),bx%hi(1)
         x = (dble(i)+0.5d0) * dx
         y = (dble(j)+0.5d0) * dx
         rhs(i,j) = sin(tpi*x) + sin(tpi*y)
      end do
      end do

  end subroutine init_sin_2d

  subroutine init_sin_3d(rhs,bx,dx)

      type(box)       , intent(in   ) :: bx
      double precision, intent(inout) :: rhs(bx%lo(1):,bx%lo(2):,bx%lo(3):)
      double precision, intent(in   ) :: dx

      integer :: i,j,k
      double precision :: tpi, epi, spi
      double precision :: x, y, z

      tpi = 8.d0 * datan(1.d0)
      epi =  4.d0 * tpi
      spi = 16.d0 * tpi

      do k = bx%lo(3),bx%hi(3)
      do j = bx%lo(2),bx%hi(2)
      do i = bx%lo(1),bx%hi(1)
         x = (dble(i)+0.5d0) * dx
         y = (dble(j)+0.5d0) * dx
         z = (dble(k)+0.5d0) * dx
         rhs(i,j,k) = sin(tpi*x) + sin(tpi*y) + sin(tpi*z)
      end do
      end do
      end do

  end subroutine init_sin_3d

  subroutine mf_init_exact(mf,pd)

    type(multifab)  , intent(inout) :: mf
    type(box     )  , intent(in   ) :: pd

    type(multifab)           :: phi
    type(box)                :: bx
    real(kind=dp_t)          :: dx
    real(kind=dp_t), pointer :: rp(:,:,:,:)
    real(kind=dp_t), pointer :: pp(:,:,:,:)

    integer   :: i,dm,nx,ny

    nx = pd%hi(1) - pd%lo(1) + 1
    ny = pd%hi(2) - pd%lo(2) + 1

    call multifab_build(phi,mf%la,1,1)

    if (nx .ne. ny) then
       print *,'Not sure what to do with nx .neq. ny in mf_init_exact'
       stop
    end if

    dx = 1.d0 / dble(nx)
 
    dm = mf%dim

    print *,'Setting rhs to Lap(phi) where phi = (x^4 - x^2) * (y^4 - y^2) * (z^4 - z^2)'
    do i = 1, nfabs(mf)
       bx = get_ibox(mf, i)
       rp => dataptr(mf,i)
       pp => dataptr(phi,i)
       if (dm.eq.2) then
          call init_exact_2d(rp(:,:,1,1),pp(:,:,1,1),bx,dx)
       else if (dm.eq.3) then 
          call init_exact_3d(rp(:,:,:,1),pp(:,:,:,1),bx,dx)
       end if
    end do

  end subroutine mf_init_exact

  subroutine init_exact_2d(rhs,phi,bx,dx)

      type(box)       , intent(in   ) :: bx
      double precision, intent(inout) :: rhs(bx%lo(1)  :,bx%lo(2)  :)
      double precision, intent(inout) :: phi(bx%lo(1)-1:,bx%lo(2)-1:)
      double precision, intent(in   ) :: dx

      integer :: i,j
      double precision :: tpi, epi, spi
      double precision :: x, y

      tpi = 8.d0 * datan(1.d0)

      do j = bx%lo(2)-1,bx%hi(2)+1
      do i = bx%lo(1)-1,bx%hi(1)+1
         x = (dble(i)+0.5d0) * dx
         y = (dble(j)+0.5d0) * dx
         phi(i,j) = (x**4 - x**2) * (y**4 - y**2)
      end do
      end do

      do j = bx%lo(2),bx%hi(2)
      do i = bx%lo(1),bx%hi(1)
         rhs(i,j) = (phi(i+1,j) + phi(i-1,j) + phi(i,j+1) + phi(i,j-1) - 4.d0 * phi(i,j)) / (dx*dx)
      end do
      end do

  end subroutine init_exact_2d

  subroutine init_exact_3d(rhs,phi,bx,dx)

      type(box)       , intent(in   ) :: bx
      double precision, intent(inout) :: rhs(bx%lo(1)  :,bx%lo(2)  :,bx%lo(3)  :)
      double precision, intent(inout) :: phi(bx%lo(1)-1:,bx%lo(2)-1:,bx%lo(3)-1:)
      double precision, intent(in   ) :: dx

      integer :: i,j,k
      double precision :: tpi, epi, spi
      double precision :: x, y, z

      tpi = 8.d0 * datan(1.d0)
      epi =  4.d0 * tpi
      spi = 16.d0 * tpi

      do k = bx%lo(3)-1,bx%hi(3)+1
      do j = bx%lo(2)-1,bx%hi(2)+1
      do i = bx%lo(1)-1,bx%hi(1)+1
         x = (dble(i)+0.5d0) * dx
         y = (dble(j)+0.5d0) * dx
         z = (dble(k)+0.5d0) * dx
         phi(i,j,k) = (x**4 - x**2) * (y**4 - y**2) * (z**4 - z**2)
      end do
      end do
      end do

      do k = bx%lo(3),bx%hi(3)
      do j = bx%lo(2),bx%hi(2)
      do i = bx%lo(1),bx%hi(1)
         rhs(i,j,k) = (phi(i+1,j,k) + phi(i-1,j,k) + phi(i,j+1,k) + phi(i,j-1,k) &
                      +phi(i,j,k+1) + phi(i,j,k-1) - 6.d0 * phi(i,j,k)) / (dx*dx)
      end do
      end do
      end do

  end subroutine init_exact_3d

  subroutine mf_init_rand(mf,pd)

    type(multifab)  , intent(inout) :: mf
    type(box     )  , intent(in   ) :: pd

    type(box)                :: bx
    real(kind=dp_t)          :: dx
    real(kind=dp_t), pointer :: rp(:,:,:,:)

    integer   :: i,dm,nx,ny

    nx = pd%hi(1) - pd%lo(1) + 1
    ny = pd%hi(2) - pd%lo(2) + 1

    dx = 1.d0 / dble(nx)
 
    dm = mf%dim

    print *,'Setting rhs to random numbers on the interval [-0.5,0.5], dm = ', dm
    do i = 1, nfabs(mf)
       bx = get_ibox(mf, i)
       rp => dataptr(mf,i,bx)
       if (dm.eq.2) then
          call init_rand_2d(rp(:,:,1,1),bx,dx)
       else if (dm.eq.3) then 
          call init_rand_3d(rp(:,:,:,1),bx,dx)
       end if
    end do

  end subroutine mf_init_rand

  subroutine init_rand_2d(rhs,bx,dx)

      type(box)       , intent(in   ) :: bx
      double precision, intent(inout) :: rhs(bx%lo(1):,bx%lo(2):)
      double precision, intent(in   ) :: dx

      integer :: i,j
      double precision :: rhsAvg

      call init_genrand(1) 

      do j = bx%lo(2),bx%hi(2)
      do i = bx%lo(1),bx%hi(1)
         rhs(i,j) = genrand_real3() 
      end do
      end do

      rhsAvg = sum(rhs)/((bx%hi(1)-bx%lo(1)+1)*(bx%hi(2)-bx%lo(2)+1))
      rhs = rhs - rhsAvg

  end subroutine init_rand_2d

  subroutine init_rand_3d(rhs,bx,dx)

      type(box)       , intent(in   ) :: bx
      double precision, intent(inout) :: rhs(bx%lo(1):,bx%lo(2):,bx%lo(3):)
      double precision, intent(in   ) :: dx

      integer :: i,j,k
      double precision :: rhsAvg

      call init_genrand(1)

      do k = bx%lo(3),bx%hi(3)
      do j = bx%lo(2),bx%hi(2)
      do i = bx%lo(1),bx%hi(1)
         rhs(i,j,k) = genrand_real3()
      end do
      end do
      end do

      rhsAvg = sum(rhs)/((bx%hi(1)-bx%lo(1)+1)*(bx%hi(2)-bx%lo(2)+1)*(bx%hi(3)-bx%lo(3)+1))
      rhs = rhs - rhsAvg

  end subroutine init_rand_3d

end module cc_rhs_module
