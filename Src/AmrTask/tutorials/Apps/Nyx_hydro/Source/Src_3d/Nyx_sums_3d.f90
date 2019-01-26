
! ::
! :: ----------------------------------------------------------
! :: INPUTS / OUTPUTS:
! ::  dat        => field to be summed
! ::  rlo,rhi    => index limits of dat array
! ::  lo,hi      => index limits of tilebox
! ::  dx         => cell size
! ::  s         <=  volume-weighted sum
! :: ----------------------------------------------------------
! ::

      subroutine sum_over_level(dat,r_l1,r_l2,r_l3,r_h1,r_h2,r_h3,lo,hi,dx,s) &
        bind(C, name="sum_over_level")

        use amrex_fort_module, only : rt => amrex_real
        implicit none
        integer          :: r_l1,r_l2,r_l3,r_h1,r_h2,r_h3
        integer          :: lo(3), hi(3)
        real(rt) :: s, dx(3)
        real(rt) :: dat(r_l1:r_h1,r_l2:r_h2,r_l3:r_h3)
        
        integer          :: i, j, k

        s    = 0.d0

        do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)
                 s = s + dat(i,j,k)
              enddo
           enddo
        enddo

        s = s * (dx(1)*dx(2)*dx(3))

      end subroutine sum_over_level

! ::
! :: ----------------------------------------------------------
! :: INPUTS / OUTPUTS:
! ::  dat1       =>  first field
! ::  dat2       => second field
! ::  rlo,rhi    => index limits of dat1 array
! ::  slo,shi    => index limits of dat2 array
! ::  lo,hi      => index limits of tilebox
! ::  dx         => cell size
! ::  s         <=  volume-weighted sum
! :: ----------------------------------------------------------
! ::

      subroutine sum_product(dat1,r_l1,r_l2,r_l3,r_h1,r_h2,r_h3,&
                             dat2,s_l1,s_l2,s_l3,s_h1,s_h2,s_h3,lo,hi,dx,s) &
        bind(C, name="sum_product")

        use amrex_fort_module, only : rt => amrex_real
        implicit none
        integer          :: r_l1,r_l2,r_l3,r_h1,r_h2,r_h3
        integer          :: s_l1,s_l2,s_l3,s_h1,s_h2,s_h3
        integer          :: lo(3), hi(3)
        real(rt) :: s, dx(3)
        real(rt) :: dat1(r_l1:r_h1,r_l2:r_h2,r_l3:r_h3)
        real(rt) :: dat2(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3)

        integer          :: i, j, k

        s    = 0.d0

        do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)
                 s = s + dat1(i,j,k)*dat2(i,j,k)
              enddo
           enddo
        enddo

        s = s * (dx(1)*dx(2)*dx(3))

      end subroutine sum_product

! ::
! :: ----------------------------------------------------------
! :: INPUTS / OUTPUTS:
! ::  dat1       =>  first field
! ::  dat2       => second field
! ::  dat3       =>  third field
! ::  rlo,rhi    => index limits of rho array
! ::  lo,hi      => index limits of tilebox
! ::  dx         => cell size
! ::  s         <=  volume-weighted sum
! :: ----------------------------------------------------------
! ::

      subroutine sum_prod_prod(dat1,r_l1,r_l2,r_l3,r_h1,r_h2,r_h3,&
                               dat2,s_l1,s_l2,s_l3,s_h1,s_h2,s_h3,&
                               dat3,t_l1,t_l2,t_l3,t_h1,t_h2,t_h3,lo,hi,dx,s) &
        bind(C, name="sum_prod_prod")

        use amrex_fort_module, only : rt => amrex_real
        implicit none

        integer          :: r_l1,r_l2,r_l3,r_h1,r_h2,r_h3
        integer          :: s_l1,s_l2,s_l3,s_h1,s_h2,s_h3
        integer          :: t_l1,t_l2,t_l3,t_h1,t_h2,t_h3
        integer          :: lo(3), hi(3)
        real(rt) :: s, dx(3)
        real(rt) :: dat1(r_l1:r_h1,r_l2:r_h2,r_l3:r_h3)
        real(rt) :: dat2(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3)
        real(rt) :: dat3(t_l1:t_h1,t_l2:t_h2,t_l3:t_h3)

        integer          :: i, j, k

        s    = 0.d0

        do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)
                 s = s + dat1(i,j,k)*dat2(i,j,k)*dat3(i,j,k)
              enddo
           enddo
        enddo

        s = s * (dx(1)*dx(2)*dx(3))

      end subroutine sum_prod_prod
