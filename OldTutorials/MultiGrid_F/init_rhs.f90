module init_rhs_mod

contains


    subroutine init_rhs(rhs,rhs_type,prob_lo,dx)

        use constants
        use multifab_module

        integer         , intent(in)     :: rhs_type
        double precision, intent(in)     :: dx
        double precision, intent(in)     :: prob_lo(:)
        type(multifab)  , intent(inout)  :: rhs

        integer :: i, ng
        integer :: lo(rhs%dim), hi(rhs%dim)

        ! 1st three indices are position; 4th index tells us which component
        double precision, pointer :: dp(:,:,:,:) 

        ! When constructing RHS, # ghost cells = 0
        ng = rhs%ng

        ! Loop through all fabs in this multifab
        do i = 1,nfabs(rhs)
            dp => dataptr(rhs,i)
            lo = lwb(get_box(rhs,i))
            hi = upb(get_box(rhs,i))
            select case(rhs%dim)
                case(2)
                    select case(rhs_type)
                        case(SIN_RHS)
                            call init_rhs_sin_2d(dp(:,:,1,1),ng,lo,hi,prob_lo,dx)
                        case(RAND_RHS)
                            call init_rhs_rand_2d(dp(:,:,1,1),ng,lo,hi)
                        case(DELTA_RHS)
                            call init_rhs_delta_2d(dp(:,:,1,1),ng,lo,hi)
                    end select
                case(3)
                    select case(rhs_type)
                        case(SIN_RHS)
                            call init_rhs_sin_3d(dp(:,:,:,1),ng,lo,hi,prob_lo,dx)
                        case(RAND_RHS)
                            call init_rhs_rand_3d(dp(:,:,:,1),ng,lo,hi)
                        case(DELTA_RHS)
                            call init_rhs_delta_3d(dp(:,:,:,1),ng,lo,hi)
                    end select
            end select
        end do

    end subroutine init_rhs
                

    subroutine init_rhs_sin_2d(rhs,ng,lo,hi,prob_lo,dx)

        use constants

        implicit none

        integer         , intent(in)  :: ng
        integer         , intent(in)  :: lo(2), hi(2)
        double precision, intent(in)  :: dx
        double precision, intent(in)  :: prob_lo(2)
        double precision, intent(out) :: rhs(lo(1)-ng:,lo(2)-ng:)

        integer          :: i, j
        double precision :: x, y

        do j = lo(2),hi(2)
            y = prob_lo(2) + (dble(j)+0.5d0) * dx
            do i = lo(1),hi(1)
                x = prob_lo(1) + (dble(i)+0.5d0) * dx
                rhs(i,j) = sin(2.d0*pi*x) + sin(2.d0*pi*y)
            end do
        end do

    end subroutine init_rhs_sin_2d


    subroutine init_rhs_sin_3d(rhs,ng,lo,hi,prob_lo,dx)

        use constants

        implicit none

        integer         , intent(in)  :: ng
        integer         , intent(in)  :: lo(3), hi(3)
        double precision, intent(in)  :: dx
        double precision, intent(in)  :: prob_lo(3)
        double precision, intent(out) :: rhs(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)

        integer          :: i, j, k
        double precision :: x, y, z

        do k = lo(3),hi(3)
            z = prob_lo(3) + (dble(k)+0.5d0) * dx
            do j = lo(2),hi(2)
                y = prob_lo(2) + (dble(j)+0.5d0) * dx
                do i = lo(1),hi(1)
                    x = prob_lo(1) + (dble(i)+0.5d0) * dx
                    rhs(i,j,k) = sin(2.d0*pi*x) + sin(2.d0*pi*y) + sin(2.d0*pi*z)
                end do
            end do
        end do

    end subroutine init_rhs_sin_3d


    subroutine init_rhs_rand_2d(rhs,ng,lo,hi)

        use constants
        use mt19937_module ! in BoxLib/Src/F_BaseLib

        implicit none

        integer         , intent(in)  :: ng
        integer         , intent(in)  :: lo(2), hi(2)
        double precision, intent(out) :: rhs(lo(1)-ng:,lo(2)-ng:)

        integer          :: i, j
        double precision :: rhsAvg

        ! Just an arbitrary seed; random seed can also be used
        call init_genrand(lo(1)*hi(1))

        do j = lo(2),hi(2)
            do i = lo(1),hi(1)
                rhs(i,j) = genrand_real3()
            end do
        end do

        ! Offset RHS to ensure that RHS sums to zero
        rhsAvg = sum(rhs(:,:)) / (hi(1)-lo(1)+1)**2
        rhs(:,:) = rhs(:,:) - rhsAvg

    end subroutine init_rhs_rand_2d


    subroutine init_rhs_rand_3d(rhs,ng,lo,hi)

        use constants
        use mt19937_module ! in BoxLib/Src/F_BaseLib

        implicit none

        integer         , intent(in)  :: ng
        integer         , intent(in)  :: lo(3), hi(3)
        double precision, intent(out) :: rhs(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)

        integer          :: i, j, k
        double precision :: rhsAvg

        ! Just an arbitrary seed; random seed can also be used
        call init_genrand(lo(1)*hi(1))

        do k = lo(3),hi(3)
            do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                    rhs(i,j,k) = genrand_real3()
                end do
            end do
        end do

        ! Offset RHS to ensure that RHS sums to zero
        rhsAvg = sum(rhs(:,:,:)) / (hi(1)-lo(1)+1)**3
        rhs(:,:,:) = rhs(:,:,:) - rhsAvg

    end subroutine init_rhs_rand_3d


    subroutine init_rhs_delta_2d(rhs,ng,lo,hi)
    ! Constructs an RHS in which each fab is zero everywhere except for a -1 in the
    ! third quadrant and 1 in the first quadrant.

        use constants

        implicit none

        integer         , intent(in)  :: ng
        integer         , intent(in)  :: lo(2), hi(2)
        double precision, intent(out) :: rhs(lo(1)-ng:,lo(2)-ng:)

        integer :: n

        n = (hi(1)-lo(1) + 1)

        rhs(:,:) = 0
        rhs(lo(1)+n/4  ,lo(2)+n/4) = -1
        rhs(lo(1)+3*n/4,lo(2)+3*n/4) = 1


    end subroutine init_rhs_delta_2d


    subroutine init_rhs_delta_3d(rhs,ng,lo,hi)
    ! Constructs an RHS in which each fab is zero everywhere except for a -1 in the
    ! third quadrant and 1 in the first quadrant.

        use constants

        implicit none

        integer         , intent(in)  :: ng
        integer         , intent(in)  :: lo(3), hi(3)
        double precision, intent(out) :: rhs(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)

        integer :: n 

        n = (hi(1)-lo(1) + 1)

        rhs(:,:,:) = 0
        rhs(lo(1)+n/4  ,lo(2)+n/4  ,lo(3)+n/4) = -1
        rhs(lo(1)+3*n/4,lo(2)+3*n/4,lo(3)+3*n/4) = 1

    end subroutine init_rhs_delta_3d


end module init_rhs_mod
