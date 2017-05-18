module traverse_mod

    implicit none

contains

    subroutine gsrb(u,f,bb,dx)

        use multifab_module

        implicit none

        double precision, intent(in)  :: dx ! Grid spacing
        type(multifab), intent(in)    :: f ! RHS
        type(multifab), intent(in)    :: bb(:) ! Nodal coefficients
        type(multifab), intent(inout) :: u ! Solution

        integer :: i
        integer :: lo(f%dim), hi(f%dim)
        double precision, pointer :: dp_u(:,:,:,:)
        double precision, pointer :: dp_f(:,:,:,:)
        double precision, pointer :: dp_bx(:,:,:,:)
        double precision, pointer :: dp_by(:,:,:,:)
        double precision, pointer :: dp_bz(:,:,:,:)

        ! Loop through all fabs in multifab, smoothing red cells only
        do i = 1,nfabs(u)

            ! Extract data from fab
            dp_u => dataptr(u,i)
            dp_f => dataptr(f,i) 
            dp_bx => dataptr(bb(1),i)
            dp_by => dataptr(bb(2),i)

            ! Get lower and upper (cell-centered) indices of fab
            lo = lwb(get_box(f,i))
            hi = upb(get_box(f,i))

            ! Red sweep
            select case(u%dim)
                case(2)
                    call gsrb_r_2d(dp_u(:,:,1,1), dp_f(:,:,1,1), dp_bx(:,:,1,1), dp_by(:,:,1,1), u%ng, lo, hi, dx)
                case(3)
                    dp_bz => dataptr(bb(3),i)
                    call gsrb_r_3d(dp_u(:,:,:,1), dp_f(:,:,:,1), dp_bx(:,:,:,1), dp_by(:,:,:,1), dp_bz(:,:,:,1), u%ng, lo, hi, dx)
            end select

        end do

        ! Fill in red ghost cells
        call multifab_fill_boundary(u)

        ! Loop through all fabs in multifab, smoothing black cells only
        do i = 1,nfabs(u)

            ! Extract data from fab
            dp_u => dataptr(u,i)
            dp_f => dataptr(f,i) 
            dp_bx => dataptr(bb(1),i)
            dp_by => dataptr(bb(2),i)

            ! Get lower and upper (cell-centered) indices of fab
            lo = lwb(get_box(f,i))
            hi = upb(get_box(f,i))

            ! Black sweep
            select case(u%dim)
                case(2)
                    call gsrb_b_2d(dp_u(:,:,1,1), dp_f(:,:,1,1), dp_bx(:,:,1,1), dp_by(:,:,1,1), u%ng, lo, hi, dx)
                case(3)
                    dp_bz => dataptr(bb(3),i)
                    call gsrb_b_3d(dp_u(:,:,:,1), dp_f(:,:,:,1), dp_bx(:,:,:,1), dp_by(:,:,:,1), dp_bz(:,:,:,1), u%ng, lo, hi, dx)
            end select

        end do

        ! Fill in black ghost cells
        call multifab_fill_boundary(u)

    end subroutine gsrb


    subroutine gsrb_r_2d(u,f,bx,by,ng,lo,hi,dx)

        implicit none

        integer         , intent(in)  :: ng ! Will always be 1 for us
        integer         , intent(in)  :: lo(2), hi(2)
        double precision, intent(in)  :: dx
        double precision, intent(in)  :: f(lo(1):, lo(2):)
        double precision, intent(in)  :: bx(lo(1):, lo(2):)
        double precision, intent(in)  :: by(lo(1):, lo(2):)
        double precision, intent(out) :: u(lo(1)-ng:, lo(2)-ng:)

        integer :: i, j

        ! Red cells
        do j = lo(2),hi(2),2
            do i = lo(1)+1,hi(1),2 ! odd x, even y
                u(i,j) = ( bx(i+1,j  )*u(i+1,j  ) + bx(i,j)*u(i-1,j  ) &
                         + by(i  ,j+1)*u(i  ,j+1) + by(i,j)*u(i  ,j-1) &
                         - dx**2*f(i,j) ) / &
                         ( bx(i+1,j) + bx(i,j) &
                         + by(i,j+1) + by(i,j) )
            end do
        end do
        do j = lo(2)+1,hi(2),2
            do i = lo(1),hi(1),2 ! even x, odd y
                u(i,j) = ( bx(i+1,j  )*u(i+1,j  ) + bx(i,j)*u(i-1,j  ) &
                         + by(i  ,j+1)*u(i  ,j+1) + by(i,j)*u(i  ,j-1) &
                         - dx**2*f(i,j) ) / &
                         ( bx(i+1,j) + bx(i,j) &
                         + by(i,j+1) + by(i,j) )
            end do
        end do

    end subroutine gsrb_r_2d


    subroutine gsrb_b_2d(u,f,bx,by,ng,lo,hi,dx)

        implicit none

        integer         , intent(in)  :: ng ! Will always be 1 for us
        integer         , intent(in)  :: lo(2), hi(2)
        double precision, intent(in)  :: dx
        double precision, intent(in)  :: f(lo(1):, lo(2):)
        double precision, intent(in)  :: bx(lo(1):, lo(2):)
        double precision, intent(in)  :: by(lo(1):, lo(2):)
        double precision, intent(out) :: u(lo(1)-ng:, lo(2)-ng:)

        integer :: i, j

        ! Black cells
        do j = lo(2),hi(2),2
            do i = lo(1),hi(1),2 ! even x, even y
                u(i,j) = ( bx(i+1,j  )*u(i+1,j  ) + bx(i,j)*u(i-1,j  ) &
                         + by(i  ,j+1)*u(i  ,j+1) + by(i,j)*u(i  ,j-1) &
                         - dx**2*f(i,j) ) / &
                         ( bx(i+1,j) + bx(i,j) &
                         + by(i,j+1) + by(i,j) )
            end do
        end do
        do j = lo(2)+1,hi(2),2 
            do i = lo(1)+1,hi(1),2 ! odd x, odd y
                u(i,j) = ( bx(i+1,j  )*u(i+1,j  ) + bx(i,j)*u(i-1,j  ) &
                         + by(i  ,j+1)*u(i  ,j+1) + by(i,j)*u(i  ,j-1) &
                         - dx**2*f(i,j) ) / &
                         ( bx(i+1,j) + bx(i,j) &
                         + by(i,j+1) + by(i,j) )
            end do
        end do

    end subroutine gsrb_b_2d


    subroutine gsrb_r_3d(u,f,bx,by,bz,ng,lo,hi,dx)

        implicit none

        integer         , intent(in)  :: ng ! Will always be 1 for us
        integer         , intent(in)  :: lo(3), hi(3)
        double precision, intent(in)  :: dx
        double precision, intent(in)  :: f(lo(1):, lo(2):, lo(3):)
        double precision, intent(in)  :: bx(lo(1):, lo(2):, lo(3):)
        double precision, intent(in)  :: by(lo(1):, lo(2):, lo(3):)
        double precision, intent(in)  :: bz(lo(1):, lo(2):, lo(3):)
        double precision, intent(out) :: u(lo(1)-ng:, lo(2)-ng:, lo(3)-ng:)

        integer :: i, j, k
    
        ! RED
        do k = lo(3),hi(3)
            if (mod(k,2)==0) then
                do j = lo(2),hi(2)
                    if (mod(j,2)==0) then
                        do i = lo(1)+1,hi(1),2 ! odd x, even y, even z
                            u(i,j,k) = ( bx(i+1,j,  k  )*u(i+1,j,  k  ) + bx(i,j,k)*u(i-1,j  ,k  ) &
                                       + by(i  ,j+1,k  )*u(i  ,j+1,k  ) + by(i,j,k)*u(i  ,j-1,k  ) &
                                       + bz(i  ,j  ,k+1)*u(i  ,j  ,k+1) + bz(i,j,k)*u(i  ,j  ,k-1) &
                                       - dx**2*f(i,j,k)) / &
                                       ( bx(i+1,j,k) + by(i,j+1,k) + bz(i,j,k+1) &
                                       + bx(i  ,j,k) + by(i,j  ,k) + bz(i,j,k  ) )
                        end do
                    else
                        do i = lo(1),hi(1),2 ! even x, odd y, even z
                            u(i,j,k) = ( bx(i+1,j,  k  )*u(i+1,j,  k  ) + bx(i,j,k)*u(i-1,j  ,k  ) &
                                       + by(i  ,j+1,k  )*u(i  ,j+1,k  ) + by(i,j,k)*u(i  ,j-1,k  ) &
                                       + bz(i  ,j  ,k+1)*u(i  ,j  ,k+1) + bz(i,j,k)*u(i  ,j  ,k-1) &
                                       - dx**2*f(i,j,k)) / &
                                       ( bx(i+1,j,k) + by(i,j+1,k) + bz(i,j,k+1) &
                                       + bx(i  ,j,k) + by(i,j  ,k) + bz(i,j,k  ) )
                        end do
                    end if
                end do
            else
                do j = lo(2),hi(2)
                    if (mod(j,2)==0) then
                        do i = lo(1),hi(1),2 ! even x, even y, odd z
                            u(i,j,k) = ( bx(i+1,j,  k  )*u(i+1,j,  k  ) + bx(i,j,k)*u(i-1,j  ,k  ) &
                                       + by(i  ,j+1,k  )*u(i  ,j+1,k  ) + by(i,j,k)*u(i  ,j-1,k  ) &
                                       + bz(i  ,j  ,k+1)*u(i  ,j  ,k+1) + bz(i,j,k)*u(i  ,j  ,k-1) &
                                       - dx**2*f(i,j,k)) / &
                                       ( bx(i+1,j,k) + by(i,j+1,k) + bz(i,j,k+1) &
                                       + bx(i  ,j,k) + by(i,j  ,k) + bz(i,j,k  ) )
                        end do
                    else
                        do i = lo(1)+1,hi(1),2 ! odd x, odd y, odd z
                            u(i,j,k) = ( bx(i+1,j,  k  )*u(i+1,j,  k  ) + bx(i,j,k)*u(i-1,j  ,k  ) &
                                       + by(i  ,j+1,k  )*u(i  ,j+1,k  ) + by(i,j,k)*u(i  ,j-1,k  ) &
                                       + bz(i  ,j  ,k+1)*u(i  ,j  ,k+1) + bz(i,j,k)*u(i  ,j  ,k-1) &
                                       - dx**2*f(i,j,k)) / &
                                       ( bx(i+1,j,k) + by(i,j+1,k) + bz(i,j,k+1) &
                                       + bx(i  ,j,k) + by(i,j  ,k) + bz(i,j,k  ) )
                        end do
                    end if
                end do
            end if
        end do

    end subroutine gsrb_r_3d

    subroutine gsrb_b_3d(u,f,bx,by,bz,ng,lo,hi,dx)

        implicit none

        integer         , intent(in)  :: ng ! Will always be 1 for us
        integer         , intent(in)  :: lo(3), hi(3)
        double precision, intent(in)  :: dx
        double precision, intent(in)  :: f(lo(1):, lo(2):, lo(3):)
        double precision, intent(in)  :: bx(lo(1):, lo(2):, lo(3):)
        double precision, intent(in)  :: by(lo(1):, lo(2):, lo(3):)
        double precision, intent(in)  :: bz(lo(1):, lo(2):, lo(3):)
        double precision, intent(out) :: u(lo(1)-ng:, lo(2)-ng:, lo(3)-ng:)

        integer :: i, j, k
        
        ! BLACK
        do k = lo(3),hi(3)
            if (mod(k,2)==0) then
                do j = lo(2),hi(2)
                    if (mod(j,2)==0) then
                        do i = lo(1),hi(1),2 ! even x, even y, even z
                            u(i,j,k) = ( bx(i+1,j,  k  )*u(i+1,j,  k  ) + bx(i,j,k)*u(i-1,j  ,k  ) &
                                       + by(i  ,j+1,k  )*u(i  ,j+1,k  ) + by(i,j,k)*u(i  ,j-1,k  ) &
                                       + bz(i  ,j  ,k+1)*u(i  ,j  ,k+1) + bz(i,j,k)*u(i  ,j  ,k-1) &
                                       - dx**2*f(i,j,k)) / &
                                       ( bx(i+1,j,k) + by(i,j+1,k) + bz(i,j,k+1) &
                                       + bx(i  ,j,k) + by(i,j  ,k) + bz(i,j,k  ) )
                        end do
                    else
                        do i = lo(1)+1,hi(1),2 ! odd x, odd y, even z
                            u(i,j,k) = ( bx(i+1,j,  k  )*u(i+1,j,  k  ) + bx(i,j,k)*u(i-1,j  ,k  ) &
                                       + by(i  ,j+1,k  )*u(i  ,j+1,k  ) + by(i,j,k)*u(i  ,j-1,k  ) &
                                       + bz(i  ,j  ,k+1)*u(i  ,j  ,k+1) + bz(i,j,k)*u(i  ,j  ,k-1) &
                                       - dx**2*f(i,j,k)) / &
                                       ( bx(i+1,j,k) + by(i,j+1,k) + bz(i,j,k+1) &
                                       + bx(i  ,j,k) + by(i,j  ,k) + bz(i,j,k  ) )
                        end do
                    end if
                end do
            else
                do j = lo(2),hi(2)
                    if (mod(j,2)==0) then
                        do i = lo(1)+1,hi(1),2 ! odd x, even y, odd z
                            u(i,j,k) = ( bx(i+1,j,  k  )*u(i+1,j,  k  ) + bx(i,j,k)*u(i-1,j  ,k  ) &
                                       + by(i  ,j+1,k  )*u(i  ,j+1,k  ) + by(i,j,k)*u(i  ,j-1,k  ) &
                                       + bz(i  ,j  ,k+1)*u(i  ,j  ,k+1) + bz(i,j,k)*u(i  ,j  ,k-1) &
                                       - dx**2*f(i,j,k)) / &
                                       ( bx(i+1,j,k) + by(i,j+1,k) + bz(i,j,k+1) &
                                       + bx(i  ,j,k) + by(i,j  ,k) + bz(i,j,k  ) )
                        end do
                    else
                        do i = lo(1),hi(1),2 ! even x, odd y, odd z
                            u(i,j,k) = ( bx(i+1,j,  k  )*u(i+1,j,  k  ) + bx(i,j,k)*u(i-1,j  ,k  ) &
                                       + by(i  ,j+1,k  )*u(i  ,j+1,k  ) + by(i,j,k)*u(i  ,j-1,k  ) &
                                       + bz(i  ,j  ,k+1)*u(i  ,j  ,k+1) + bz(i,j,k)*u(i  ,j  ,k-1) &
                                       - dx**2*f(i,j,k)) / &
                                       ( bx(i+1,j,k) + by(i,j+1,k) + bz(i,j,k+1) &
                                       + bx(i  ,j,k) + by(i,j  ,k) + bz(i,j,k  ) )
                        end do
                    end if
                end do
            end if
        end do

    end subroutine gsrb_b_3d 


    subroutine restriction(crse,fine)

        use multifab_module

        implicit none

        type(multifab), intent(in)    :: fine
        type(multifab), intent(inout) :: crse

        double precision, pointer :: dp_crse(:,:,:,:), dp_fine(:,:,:,:)

        integer :: i
        integer :: fine_lo(fine%dim), fine_hi(fine%dim), crse_lo(crse%dim), crse_hi(crse%dim)

        do i = 1,nfabs(fine)
            dp_fine => dataptr(fine,i)
            dp_crse => dataptr(crse,i)
            crse_lo = lwb(get_box(crse,i))
            crse_hi = upb(get_box(crse,i))
            fine_lo = lwb(get_box(fine,i))
            fine_hi = upb(get_box(fine,i))
            select case(fine%dim)
                case(2)
                    call restriction_2d(dp_crse(:,:,1,1),dp_fine(:,:,1,1),crse_lo,crse_hi,fine_lo,fine_hi)
                case(3)
                    call restriction_3d(dp_crse(:,:,:,1),dp_fine(:,:,:,1),crse_lo,crse_hi,fine_lo,fine_hi)
            end select
        end do

    end subroutine restriction


    subroutine restriction_2d(crse,fine,crse_lo,crse_hi,fine_lo,fine_hi)

        implicit none

        integer, intent(in) :: crse_lo(2), crse_hi(2), fine_lo(2), fine_hi(2)
        double precision, intent(in)  :: fine(fine_lo(1):,fine_lo(2):)
        double precision, intent(out) :: crse(crse_lo(1):,crse_lo(2):)

        integer :: i,j

        do j = crse_lo(2),crse_hi(2)
            do i = crse_lo(1),crse_hi(1)
                crse(i,j) = (fine(2*i,2*j)+fine(2*i+1,2*j)+fine(2*i,2*j+1)+fine(2*i+1,2*j+1)) / 4.d0
            end do
        end do

    end subroutine restriction_2d
    
        
    subroutine restriction_3d(crse,fine,crse_lo,crse_hi,fine_lo,fine_hi)

        implicit none

        integer, intent(in) :: crse_lo(3), crse_hi(3), fine_lo(3), fine_hi(3)
        double precision, intent(in)  :: fine(fine_lo(1):,fine_lo(2):,fine_lo(3):)
        double precision, intent(out) :: crse(crse_lo(1):,crse_lo(2):,crse_lo(3):)

        integer :: i,j,k

        do k = crse_lo(3), crse_hi(3)
            do j = crse_lo(2), crse_hi(2)
                do i = crse_lo(1), crse_hi(1)
                    crse(i,j,k) = (fine(2*i  ,2*j  ,2*k) + fine(2*i  ,2*j  ,2*k+1) &
                                 + fine(2*i+1,2*j  ,2*k) + fine(2*i+1,2*j  ,2*k+1) &
                                 + fine(2*i  ,2*j+1,2*k) + fine(2*i  ,2*j+1,2*k+1) &
                                 + fine(2*i+1,2*j+1,2*k) + fine(2*i+1,2*j+1,2*k+1)) / 8.d0
                end do
            end do
        end do

    end subroutine restriction_3d


    subroutine prolongation(fine,crse,prolongation_type)

        use multifab_module

        implicit none

        integer, intent(in) :: prolongation_type ! Either piecewise constant or piecewise linear
        type(multifab), intent(in)    :: crse
        type(multifab), intent(inout) :: fine

        integer :: i
        integer :: crse_lo(crse%dim), crse_hi(crse%dim), fine_lo(fine%dim), fine_hi(fine%dim)
        double precision, pointer :: dp_crse(:,:,:,:), dp_fine(:,:,:,:)

        do i = 1,nfabs(fine)
            dp_fine => dataptr(fine,i)
            dp_crse => dataptr(crse,i)
            crse_lo = lwb(get_box(crse,i))
            crse_hi = upb(get_box(crse,i))
            fine_lo = lwb(get_box(fine,i))
            fine_hi = upb(get_box(fine,i))
            select case(fine%dim)
                case(2)
                    call prolongation_2d(dp_fine(:,:,1,1),dp_crse(:,:,1,1),crse%ng,&
                                         fine_lo,fine_hi,crse_lo,crse_hi,prolongation_type)
                case(3)
                    call prolongation_3d(dp_fine(:,:,:,1),dp_crse(:,:,:,1),crse%ng,&
                                         fine_lo,fine_hi,crse_lo,crse_hi,prolongation_type)
            end select
        end do

        call multifab_fill_boundary(fine)

    end subroutine prolongation


    subroutine prolongation_2d(fine,crse,ng,fine_lo,fine_hi,crse_lo,crse_hi,prolongation_type)

        use constants

        implicit none

        integer, intent(in) :: ng, prolongation_type
        integer, intent(in) :: crse_lo(2), crse_hi(2), fine_lo(2), fine_hi(2)
        double precision, intent(in)  :: crse(crse_lo(1)-ng:,crse_lo(2)-ng:)
        double precision, intent(out) :: fine(fine_lo(1)-ng:,fine_lo(2)-ng:)

        integer :: i, j, k, ii, jj, kk
        double precision, parameter :: c1 = 9.d0/16.d0
        double precision, parameter :: c2 = 3.d0/16.d0
        double precision, parameter :: c3 = 1.d0/16.d0

        if (prolongation_type == CONST_INTERP) then
            do j = fine_lo(2), fine_hi(2)
                do i = fine_lo(1), fine_hi(1)
                    fine(i,j) = fine(i,j) + crse(i/2,j/2)
                end do
            end do 
        else if (prolongation_type == LIN_INTERP) then
            do  j = fine_lo(2),fine_hi(2)
                jj = j/2
                do i = fine_lo(1),fine_hi(1)
                    ii = i/2
                    if (mod(i,2)==0 .and. mod(j,2)==0) then ! even x, even y
                        fine(i,j) = fine(i,j) + c1*crse(ii,jj) + c2*crse(ii-1,jj) + c2*crse(ii,jj-1) + c3*crse(ii-1,jj-1)
                    else if (mod(i,2)==0 .and. mod(j,2)==1) then ! even x, odd y
                        fine(i,j) = fine(i,j) + c1*crse(ii,jj) + c2*crse(ii-1,jj) + c2*crse(ii,jj+1) + c3*crse(ii-1,jj+1)
                    else if (mod(i,2)==1 .and. mod(j,2)==0) then ! odd x, even y
                        fine(i,j) = fine(i,j) + c1*crse(ii,jj) + c2*crse(ii+1,jj) + c2*crse(ii,jj-1) + c3*crse(ii+1,jj-1)
                    else if (mod(i,2)==1 .and. mod(j,2)==1) then ! odd x, odd y
                        fine(i,j) = fine(i,j) + c1*crse(ii,jj) + c2*crse(ii+1,jj) + c2*crse(ii,jj+1) + c3*crse(ii+1,jj+1)
                    end if
                end do
            end do
        end if

    end subroutine prolongation_2d


    subroutine prolongation_3d(fine,crse,ng,fine_lo,fine_hi,crse_lo,crse_hi,prolongation_type)

        use constants

        implicit none

        integer, intent(in) :: ng, prolongation_type
        integer, intent(in) :: crse_lo(3), crse_hi(3), fine_lo(3), fine_hi(3)
        double precision, intent(in)  :: crse(crse_lo(1)-ng:,crse_lo(2)-ng:,crse_lo(3)-ng:)
        double precision, intent(out) :: fine(fine_lo(1)-ng:,fine_lo(2)-ng:,fine_lo(3)-ng:)

        integer :: i, j, k, ii, jj, kk
        double precision, parameter :: c1 = 27.d0/64.d0
        double precision, parameter :: c2 = 9.d0/64.d0
        double precision, parameter :: c3 = 3.d0/64.d0
        double precision, parameter :: c4 = 1.d0/64.d0

        if (prolongation_type == CONST_INTERP) then
            do k = fine_lo(3),fine_hi(3)
                do j = fine_lo(2),fine_hi(2)
                    do i = fine_lo(1),fine_hi(1)
                       fine(i,j,k) = fine(i,j,k) + crse(i/2,j/2,k/2)
                    end do
                end do
            end do
        else if (prolongation_type == LIN_INTERP) then
            do k = fine_lo(3),fine_hi(3)
                do j = fine_lo(2),fine_hi(2)
                    do i = fine_lo(1),fine_hi(1)
                        ii = i/2
                        jj = j/2
                        kk = k/2
                        if (mod(i,2)==0 .and. mod(j,2)==0 .and. mod(k,2)==0) then !  even x, even y, even z

                            fine(i,j,k) = c1*crse(ii,jj,kk) + c4*crse(ii-1,jj-1,kk-1) + fine(i,j,k) &
                                        + c2*(crse(ii-1,jj  ,kk  ) + crse(ii  ,jj-1,kk  ) + crse(ii  ,jj  ,kk-1)) &
                                        + c3*(crse(ii  ,jj-1,kk-1) + crse(ii-1,jj  ,kk-1) + crse(ii-1,jj-1,kk  ))

                        else if (mod(i,2)==0 .and. mod(j,2)==0 .and. mod(k,2)==1) then !  even x, even y, odd z

                            fine(i,j,k) = c1*crse(ii,jj,kk) + c4*crse(ii-1,jj-1,kk+1) + fine(i,j,k) &
                                        + c2*(crse(ii-1,jj  ,kk  ) + crse(ii  ,jj-1,kk  ) + crse(ii  ,jj  ,kk+1)) &
                                        + c3*(crse(ii  ,jj-1,kk+1) + crse(ii-1,jj  ,kk+1) + crse(ii-1,jj-1,kk  ))

                        else if (mod(i,2)==0 .and. mod(j,2)==1 .and. mod(k,2)==0) then !  even x, odd y, even z

                            fine(i,j,k) = c1*crse(ii,jj,kk) + c4*crse(ii-1,jj+1,kk-1) + fine(i,j,k) &
                                        + c2*(crse(ii-1,jj  ,kk  ) + crse(ii  ,jj+1,kk  ) + crse(ii  ,jj  ,kk-1)) &
                                        + c3*(crse(ii  ,jj+1,kk-1) + crse(ii-1,jj  ,kk-1) + crse(ii-1,jj+1,kk  ))

                        else if (mod(i,2)==0 .and. mod(j,2)==1 .and. mod(k,2)==1) then !  even x, odd y, odd z

                            fine(i,j,k) = c1*crse(ii,jj,kk) + c4*crse(ii-1,jj+1,kk+1) + fine(i,j,k) &
                                        + c2*(crse(ii-1,jj  ,kk  ) + crse(ii  ,jj+1,kk  ) + crse(ii  ,jj  ,kk+1)) &
                                        + c3*(crse(ii  ,jj+1,kk+1) + crse(ii-1,jj  ,kk+1) + crse(ii-1,jj+1,kk  ))

                        else if (mod(i,2)==1 .and. mod(j,2)==0 .and. mod(k,2)==0) then ! odd x, even y, even z

                            fine(i,j,k) = c1*crse(ii,jj,kk) + c4*crse(ii+1,jj-1,kk-1) + fine(i,j,k) &
                                        + c2*(crse(ii+1,jj  ,kk  ) + crse(ii  ,jj-1,kk  ) + crse(ii  ,jj  ,kk-1)) &
                                        + c3*(crse(ii  ,jj-1,kk-1) + crse(ii+1,jj  ,kk-1) + crse(ii+1,jj-1,kk  ))

                        else if (mod(i,2)==1 .and. mod(j,2)==0 .and. mod(k,2)==1) then ! odd x, even y, odd z

                            fine(i,j,k) = c1*crse(ii,jj,kk) + c4*crse(ii+1,jj-1,kk+1) + fine(i,j,k) &
                                        + c2*(crse(ii+1,jj  ,kk  ) + crse(ii  ,jj-1,kk  ) + crse(ii  ,jj  ,kk+1)) &
                                        + c3*(crse(ii  ,jj-1,kk+1) + crse(ii+1,jj  ,kk+1) + crse(ii+1,jj-1,kk  ))

                        else if (mod(i,2)==1 .and. mod(j,2)==1 .and. mod(k,2)==0) then ! odd x, odd y, even z

                            fine(i,j,k) = c1*crse(ii,jj,kk) + c4*crse(ii+1,jj+1,kk-1) + fine(i,j,k) &
                                        + c2*(crse(ii+1,jj  ,kk  ) + crse(ii  ,jj+1,kk  ) + crse(ii  ,jj  ,kk-1)) &
                                        + c3*(crse(ii  ,jj+1,kk-1) + crse(ii+1,jj  ,kk-1) + crse(ii+1,jj+1,kk  )) 

                        else if (mod(i,2)==1 .and. mod(j,2)==1 .and. mod(k,2)==1) then ! odd x, odd y, odd z

                            fine(i,j,k) = c1*crse(ii,jj,kk) + c4*crse(ii+1,jj+1,kk+1) + fine(i,j,k) &
                                        + c2*(crse(ii+1,jj  ,kk  ) + crse(ii  ,jj+1,kk  ) + crse(ii  ,jj  ,kk+1)) &
                                        + c3*(crse(ii  ,jj+1,kk+1) + crse(ii+1,jj  ,kk+1) + crse(ii+1,jj+1,kk  )) 

                        end if
                    end do
                end do
            end do
        end if

    end subroutine prolongation_3d


    subroutine residual(res,u,f,bb,dx)

        use multifab_module

        implicit none

        double precision, intent(in) :: dx ! Grid spacing
        type(multifab), intent(in)   :: u ! Solution
        type(multifab), intent(in)   :: f ! RHS
        type(multifab), intent(in)   :: bb(:) ! Nodal coefficients
        type(multifab), intent(inout)  :: res ! Residual

        integer :: i
        integer :: lo(u%dim), hi(u%dim)
        double precision, pointer :: dp_u(:,:,:,:)
        double precision, pointer :: dp_f(:,:,:,:)
        double precision, pointer :: dp_res(:,:,:,:)
        double precision, pointer :: dp_bx(:,:,:,:)
        double precision, pointer :: dp_by(:,:,:,:)
        double precision, pointer :: dp_bz(:,:,:,:)

        ! Loop through all fabs in multifab
        do i = 1,nfabs(u)

            ! Extract data from each fab
            dp_res => dataptr(res,i)
            dp_u => dataptr(u,i)
            dp_f => dataptr(f,i)
            dp_bx => dataptr(bb(1),i)
            dp_by => dataptr(bb(2),i)

            ! Get lower and upper (cell-centered) indices of fab
            lo = lwb(get_box(u,i))
            hi = upb(get_box(u,i))

            ! Calculate residual
            select case(u%dim)
                case(2)
                    call residual_2d(dp_res(:,:,1,1), dp_u(:,:,1,1), dp_f(:,:,1,1), dp_bx(:,:,1,1), dp_by(:,:,1,1),&
                                     u%ng, lo, hi, dx)
                case(3)
                    dp_bz => dataptr(bb(3),i)
                    call residual_3d(dp_res(:,:,:,1), dp_u(:,:,:,1), dp_f(:,:,:,1), dp_bx(:,:,:,1), dp_by(:,:,:,1), dp_bz(:,:,:,1),&
                                     u%ng, lo, hi, dx)
            end select

        end do

    end subroutine residual


    subroutine residual_2d(res,u,f,bx,by,ng,lo,hi,dx)

        implicit none

        integer         , intent(in)  :: ng ! Will always be 1 for us
        integer         , intent(in)  :: lo(2), hi(2)
        double precision, intent(in)  :: dx
        double precision, intent(in)  :: f(lo(1):,lo(2):)
        double precision, intent(in)  :: bx(lo(1):,lo(2):)
        double precision, intent(in)  :: by(lo(1):,lo(2):)
        double precision, intent(in)  :: u(lo(1)-ng:,lo(2)-ng:)
        double precision, intent(out) :: res(lo(1):,lo(2):)

        integer :: i, j
        double precision :: lapl

        do j = lo(2),hi(2)
            do i = lo(1),hi(1)
                lapl = ( bx(i+1,j  )*u(i+1,j  ) + bx(i,j)*u(i-1,j  ) &
                       + by(i  ,j+1)*u(i  ,j+1) + by(i,j)*u(i  ,j-1) &
                       -(bx(i+1,j  ) + bx(i,j) &
                       + by(i  ,j+1) + by(i,j))*u(i,j) ) / dx**2
                res(i,j) = f(i,j) - lapl
            end do
        end do

    end subroutine residual_2d

    subroutine residual_3d(res,u,f,bx,by,bz,ng,lo,hi,dx)

        implicit none

        integer         , intent(in)  :: ng ! Will always be 1 for us
        integer         , intent(in)  :: lo(3), hi(3)
        double precision, intent(in)  :: dx
        double precision, intent(in)  :: f(lo(1):,lo(2):,lo(3):)
        double precision, intent(in)  :: bx(lo(1):,lo(2):,lo(3):)
        double precision, intent(in)  :: by(lo(1):,lo(2):,lo(3):)
        double precision, intent(in)  :: bz(lo(1):,lo(2):,lo(3):)
        double precision, intent(in)  :: u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
        double precision, intent(out) :: res(lo(1):,lo(2):,lo(3):)

        integer :: i, j, k
        double precision :: lapl

        do k = lo(3),hi(3)
            do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                    lapl = ( bx(i+1,j  ,k  )*u(i+1,j  ,k  ) + bx(i,j,k)*u(i-1,j  ,k  ) &
                           + by(i  ,j+1,k  )*u(i  ,j+1,k  ) + by(i,j,k)*u(i  ,j-1,k  ) &
                           + bz(i  ,j  ,k+1)*u(i  ,j  ,k+1) + bz(i,j,k)*u(i  ,j  ,k-1) &
                           -(bx(i+1,j,k) + by(i,j+1,k) + bz(i,j,k+1) &
                           + bx(i  ,j,k) + by(i,j  ,k) + bz(i,j,k  ))*u(i,j,k) ) / dx**2
                    res(i,j,k) = f(i,j,k) - lapl
                end do
            end do
        end do

    end subroutine residual_3d


end module traverse_mod



