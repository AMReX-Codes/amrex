module init_coeffs_mod

contains


    subroutine init_coeffs(coeffs,coeffs_type)

        use constants
        use multifab_module

        implicit none

        integer       , intent(in)    :: coeffs_type
        type(multifab), intent(inout) :: coeffs ! Cell-centered coeffs

        integer :: i, ng
        integer :: lo(coeffs%dim), hi(coeffs%dim)

        ! 1st three indices are position; 4th index tells us which component
        double precision, pointer :: dp(:,:,:,:) 

        ! When constructing cell-centered coefficients, # ghost cells = 1
        ng = coeffs%ng 

        do i = 1,nfabs(coeffs)
            dp => dataptr(coeffs,i)
            lo = lwb(get_box(coeffs,i))
            hi = upb(get_box(coeffs,i))
            select case(coeffs%dim)
                case(2)
                    select case(coeffs_type)
                        case(CONST_COEFFS)
                            call init_coeffs_const_2d(dp(:,:,1,1),ng,lo,hi)
                        case(RAND_COEFFS)
                            call init_coeffs_rand_2d(dp(:,:,1,1),ng,lo,hi)
                    end select
                case(3)
                    select case(coeffs_type)
                        case(CONST_COEFFS)
                            call init_coeffs_const_3d(dp(:,:,:,1),ng,lo,hi)
                        case(RAND_COEFFS)
                            call init_coeffs_rand_3d(dp(:,:,:,1),ng,lo,hi)
                    end select
            end select
        end do

        call multifab_fill_boundary(coeffs)

    end subroutine init_coeffs


    subroutine init_coeffs_const_2d(coeffs,ng,lo,hi)

        implicit none

        integer         , intent(in) :: ng
        integer         , intent(in) :: lo(2), hi(2)
        double precision, intent(out):: coeffs(lo(1)-ng:,lo(2)-ng:)

        integer :: i,j

        do j = lo(2),hi(2)
            do i  = lo(1),hi(1)
                coeffs(i,j) = 1.d0
            end do
        end do

    end subroutine init_coeffs_const_2d


    subroutine init_coeffs_const_3d(coeffs,ng,lo,hi)

        implicit none

        integer         , intent(in)  :: ng
        integer         , intent(in)  :: lo(3), hi(3)
        double precision, intent(out) :: coeffs(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)

        integer :: i,j,k

        do k = lo(3),hi(3)
            do j = lo(2),hi(2)
                do i  = lo(1),hi(1)
                    coeffs(i,j,k) = 1.d0
                end do
            end do
        end do

    end subroutine init_coeffs_const_3d


    subroutine init_coeffs_rand_2d(coeffs,ng,lo,hi)

        use mt19937_module

        implicit none

        integer         , intent(in)  :: ng
        integer         , intent(in)  :: lo(2), hi(2)
        double precision, intent(out) :: coeffs(lo(1)-ng:,lo(2)-ng:)

        integer :: i,j

        ! Just an arbitrary seed; random seed can also be used
        call init_genrand(lo(1)*hi(1))

        do j = lo(2),hi(2)
            do i  = lo(1),hi(1)
                coeffs(i,j) = 1 + genrand_real3()
            end do
        end do

    end subroutine init_coeffs_rand_2d


    subroutine init_coeffs_rand_3d(coeffs,ng,lo,hi) 

        use mt19937_module

        implicit none

        integer         , intent(in)  :: ng
        integer         , intent(in)  :: lo(3), hi(3)
        double precision, intent(out) :: coeffs(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)

        integer :: i,j,k

        ! Just an arbitrary seed; random seed can also be used
        call init_genrand(lo(1)*hi(1))

        do k = lo(3),hi(3)
            do j = lo(2),hi(2)
                do i  = lo(1),hi(1)
                    coeffs(i,j,k) = 1 + genrand_real3()
                end do
            end do
        end do

    end subroutine init_coeffs_rand_3d


    subroutine fill_nodal_coeffs(nd,cc,dir)

        use multifab_module

        implicit none

        integer       , intent(in)    :: dir ! Direction to nodalize
        type(multifab), intent(in)    :: cc ! Cell-centered coefficients
        type(multifab), intent(inout) :: nd ! Nodalized version of coefficients

        integer :: i
        integer :: lo(nd%dim), hi(nd%dim)
        
        ! Extract data. (1st three indices are position; 4th index tells us which component)
        double precision, pointer :: dp_nd(:,:,:,:)
        double precision, pointer :: dp_cc(:,:,:,:)

        do i = 1,nfabs(cc)
            dp_nd => dataptr(nd,i)
            dp_cc => dataptr(cc,i)
            lo = lwb(get_ibox(nd,i)) ! Lower (nodal) indices of fab being operated on
            hi = upb(get_ibox(nd,i)) ! Upper (nodal) indices of fab being operated on
            select case(cc%dim)
                case(2)
                    call fill_nodal_coeffs_2d(dp_nd(:,:,1,1),dp_cc(:,:,1,1),nd%ng,cc%ng,lo,hi,dir)
                case(3)
                    call fill_nodal_coeffs_3d(dp_nd(:,:,:,1),dp_cc(:,:,:,1),nd%ng,cc%ng,lo,hi,dir)
            end select

        end do

    end subroutine fill_nodal_coeffs
    

    subroutine fill_nodal_coeffs_2d(nd,cc,nd_ng,cc_ng,lo,hi,dir)

        implicit none

        integer         , intent(in)  :: dir ! Direction to nodalize
        integer         , intent(in)  :: cc_ng ! We'll always have one ghost cell for cell-centered coeffs
        integer         , intent(in)  :: nd_ng ! We'll always have zero ghost cells for nodal coeffs
        integer         , intent(in)  :: lo(2), hi(2) 
        double precision, intent(in)  :: cc(lo(1)-cc_ng:, lo(2)-cc_ng:)
        double precision, intent(out) :: nd(lo(1)-nd_ng:, lo(2)-nd_ng:)

        integer :: i, j

        select case(dir)
            case(1)
                do j = lo(2),hi(2)
                    do i = lo(1),hi(1)
                        nd(i,j) = (cc(i-1,j)+cc(i,j))/2.d0
                    end do
                end do
            case(2)
                do i = lo(1),hi(1)
                    do j = lo(2),hi(2)
                        nd(i,j) = (cc(i,j-1)+cc(i,j))/2.d0
                    end do
                end do
        end select

    end subroutine fill_nodal_coeffs_2d


    subroutine fill_nodal_coeffs_3d(nd,cc,nd_ng,cc_ng,lo,hi,dir)

        implicit none

        integer         , intent(in)  :: dir ! Direction to nodalize
        integer         , intent(in)  :: cc_ng ! We'll always have one ghost cell for cell-centered coeffs
        integer         , intent(in)  :: nd_ng ! We'll always have zero ghost cells for nodal coeffs
        integer         , intent(in)  :: lo(3), hi(3) 
        double precision, intent(in)  :: cc(lo(1)-cc_ng:, lo(2)-cc_ng:, lo(3)-cc_ng:)
        double precision, intent(out) :: nd(lo(1)-nd_ng:, lo(2)-nd_ng:, lo(3)-nd_ng:)

        integer :: i, j, k

        select case(dir)
            case(1) 
                do k = lo(3),hi(3)
                    do j = lo(2),hi(2)
                        do i = lo(1),hi(1)
                            nd(i,j,k) = (cc(i-1,j,k)+cc(i,j,k))/2.d0
                        end do
                    end do
                end do
            case(2)
                do k = lo(3),hi(3)
                    do i = lo(1),hi(1)
                        do j = lo(2),hi(2)
                            nd(i,j,k) = (cc(i,j-1,k)+cc(i,j,k))/2.d0
                        end do
                    end do
                end do
            case(3)
                do i = lo(1),hi(1)
                    do j = lo(2),hi(2)
                        do k = lo(3),hi(3)
                            nd(i,j,k) = (cc(i,j,k-1)+cc(i,j,k))/2.d0
                        end do
                    end do
                end do
        end select
            
    end subroutine fill_nodal_coeffs_3d


    subroutine restrict_coeffs(crse,fine,dir)

        use multifab_module

        implicit none

        integer       , intent(in)    :: dir ! Direction of restriction
        type(multifab), intent(in)    :: fine ! Nodal coefficients on fine grid
        type(multifab), intent(inout) :: crse ! Nodal coefficients on coarse grid

        integer :: i, ng
        integer :: fine_lo(fine%dim), fine_hi(fine%dim), crse_lo(crse%dim), crse_hi(crse%dim)

        ! 1st three indices are position; 4th index tells us which component
        double precision, pointer :: dp_fine(:,:,:,:)
        double precision, pointer :: dp_crse(:,:,:,:)

        ! We'll always have zero ghost cells for nodal coeffs
        ng = fine%ng

        do i = 1,nfabs(fine)
            dp_fine => dataptr(fine,i)
            dp_crse => dataptr(crse,i)
            fine_lo = lwb(get_ibox(fine,i)) ! Use ibox instead of box to retrieve nodal indexing
            fine_hi = upb(get_ibox(fine,i))
            crse_lo = lwb(get_ibox(crse,i))
            crse_hi = upb(get_ibox(crse,i))
            select case(fine%dim)
                case(2)
                    call restrict_coeffs_2d(dp_crse(:,:,1,1),dp_fine(:,:,1,1),ng,crse_lo,crse_hi,fine_lo,fine_hi,dir)
                case(3)
                    call restrict_coeffs_3d(dp_crse(:,:,:,1),dp_fine(:,:,:,1),ng,crse_lo,crse_hi,fine_lo,fine_hi,dir)
            end select
        end do

    end subroutine restrict_coeffs


    subroutine restrict_coeffs_2d(crse,fine,ng,crse_lo,crse_hi,fine_lo,fine_hi,dir)

        implicit none

        integer         , intent(in)  :: dir, ng
        integer         , intent(in)  :: crse_lo(2), crse_hi(2), fine_lo(2), fine_hi(2)
        double precision, intent(in)  :: fine(fine_lo(1)-ng:,fine_lo(2)-ng:)
        double precision, intent(out) :: crse(crse_lo(1)-ng:,crse_lo(2)-ng:)
        
        integer :: i, j

        select case(dir)
            case(1)
                do j = crse_lo(2),crse_hi(2)
                    do i = crse_lo(1),crse_hi(1)
                        crse(i,j) = (fine(2*i,2*j)+fine(2*i,2*j+1))/2.d0
                    end do
                end do
            case(2)
                do i = crse_lo(1),crse_hi(1)
                    do j = crse_lo(2),crse_hi(2)
                        crse(i,j) = (fine(2*i,2*j)+fine(2*i+1,2*j))/2.d0
                    end do
                end do
        end select

    end subroutine restrict_coeffs_2d


    subroutine restrict_coeffs_3d(crse,fine,ng,crse_lo,crse_hi,fine_lo,fine_hi,dir)

        implicit none

        integer         , intent(in)  :: dir, ng
        integer         , intent(in)  :: crse_lo(3), crse_hi(3), fine_lo(3), fine_hi(3)
        double precision, intent(in)  :: fine(fine_lo(1)-ng:,fine_lo(2)-ng:,fine_lo(3)-ng:)
        double precision, intent(out) :: crse(crse_lo(1)-ng:,crse_lo(2)-ng:,crse_lo(3)-ng:)

        integer :: i, j, k

        select case(dir)
            case(1)
                do k = crse_lo(3),crse_hi(3)
                    do j = crse_lo(2),crse_hi(2)
                        do i = crse_lo(1),crse_hi(1)
                            crse(i,j,k) = (fine(2*i,2*j+1,2*k)+fine(2*i,2*j,2*k+1)&
                                                 +fine(2*i,2*j,2*k)  +fine(2*i,2*j+1,2*k+1))/4.d0
                        end do
                    end do
                end do
            case(2)
                do k = crse_lo(3),crse_hi(3)
                    do i = crse_lo(1),crse_hi(1)
                        do j = crse_lo(2),crse_hi(2)
                            crse(i,j,k) = (fine(2*i+1,2*j,2*k)+fine(2*i,2*j,2*k+1)&
                                                 +fine(2*i,2*j,2*k)  +fine(2*i+1,2*j,2*k+1))/4.d0
                        end do
                    end do
                end do
            case(3)
                do j = crse_lo(2),crse_hi(2)
                    do i = crse_lo(1),crse_hi(1)
                        do k = crse_lo(3),crse_hi(3)
                            crse(i,j,k) = (fine(2*i+1,2*j,2*k)+fine(2*i,2*j+1,2*k)&
                                                 +fine(2*i,2*j,2*k)  +fine(2*i+1,2*j+1,2*k))/4.d0
                        end do
                    end do
                end do
        end select 

    end subroutine restrict_coeffs_3d


end module init_coeffs_mod
