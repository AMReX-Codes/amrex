module coarsen_coeffs_module

  use bl_types
  use multifab_module

  implicit none

  real(kind=dp_t), private, parameter :: HALF   = 0.5_dp_t
  real(kind=dp_t), private, parameter :: FOURTH = 0.25_dp_t
  real(kind=dp_t), private, parameter :: EIGHTH = 0.125_dp_t

  private :: crse_coeffs_1d, crse_coeffs_mc_1d
  private :: crse_coeffs_2d, crse_coeffs_mc_2d
  private :: crse_coeffs_3d, crse_coeffs_mc_3d

contains

  subroutine coarsen_coeffs(cf,cc)
    type(multifab), intent(in   ) :: cf
    type(multifab), intent(inout) :: cc

    real(kind=dp_t), pointer :: cfp(:,:,:,:)
    real(kind=dp_t), pointer :: ccp(:,:,:,:)
    integer :: i,ng,nc_edge,dm
    integer :: lof(cf%dim)
    integer :: loc(cc%dim),hic(cc%dim)
    logical :: single_comp

    dm = cf%dim

    if (cf%nc .ne. cc%nc) then
       print *,'ncomp_fine not equal to ncomp_crse in coarsen-coeffs'
       print *,'ncomp_fine = ',cf%nc
       print *,'ncomp_crse = ',cc%nc
       call bl_error("coarsen_coeffs.f90 :: coarsen_coeffs")
    end if

    ! In the case of the nodal solver:
       ! we expect cf%nc = cc%nc = 1, so there is only a single cell-centered component,
       ! and nc_edge = 0

    ! In the case of the cell-centered solver::
       ! If the stencil is (alpha - del dot beta grad) then nc_edge = 1, and there are
       !     (1+dm            ) coefficients stored in coeffs, e.g. alpha,betax,betay,betaz in 3d
       ! If the stencil is (alpha - sum_{i=1}^N beta0_i (del dot beta_i grad) ) then nc_edge = N, and there are
       !     (1+(dm+1)*nc_edge) coefficients stored in coeffs, e.g. alpha,beta0_ibetax_i,betay_i,betaz_i in 3d

    ! Nodal case
    if (cf%nc .eq. 1) then

       single_comp = .true.

    ! Standard cell-centered case
    else if (cf%nc .eq. (dm+1) ) then

       single_comp = .true.

    ! Special stencil with sum over nc_edge components *and* additional cell-centered values
    else

       nc_edge = (cf%nc-1) / (dm+1)
       single_comp = .false.

       if (cf%nc .ne. (nc_edge*dm+1)+1) then
          print *,'ncomp_fine = ',cf%nc
          print *,'dm         = ',dm
          print *,'nc_edge    = ',nc_edge
          print *,'ncomp_fine != (nc_edge*(dm+1)) + 1 '
          call bl_error("coarsen_coeffs.f90 :: coarsen_coeffs")
       end if

    end if

    ng = cc%ng

    ! This is either our coarsening for the nodal solve or
    ! our standard coarsening for the cell-centered solve, where coeffs contains (dm+1) components --
    !   in 3d, for example, there are 4 components, for alpha, betax, betay, betaz
    if (single_comp) then

       do i = 1, cf%nboxes
          if ( multifab_remote(cf,i) ) cycle
          cfp => dataptr(cf, i)
          ccp => dataptr(cc, i)
   
          loc =  lwb(get_box(cc, i))
          hic =  upb(get_box(cc, i))
   
          lof =  lwb(get_box(cf, i))

          select case (cf%dim)
          case (1)
             call crse_coeffs_1d(ccp(:,1,1,:), cfp(:,1,1,:), ng, loc, hic, lof)
          case (2)
             call crse_coeffs_2d(ccp(:,:,1,:), cfp(:,:,1,:), ng, loc, hic, lof)
          case (3)
             call crse_coeffs_3d(ccp(:,:,:,:), cfp(:,:,:,:), ng, loc, hic, lof)
          end select
       end do

    else
    ! This is our special coarsening for the porous media solve, 
    !   where coeffs contains (dm*nc_edge+1) components --
    !   in 3d, for example, 
    !   component       0       is filled with alpha, 
    !   components      1:  nc_edge are filled with cell-centered beta0(1:nc_edge),
    !   components   nc+1:2*nc_edge are filled with edge-centered betax(1:nc_edge),
    !   components 2*nc+1:3*nc_edge are filled with edge-centered betay(1:nc_edge),
    !   components 3*nc+1:4*nc_edge are filled with edge-centered betaz(1:nc_edge)

       do i = 1, cf%nboxes
          if ( multifab_remote(cf,i) ) cycle
          cfp => dataptr(cf, i)
          ccp => dataptr(cc, i)
   
          loc =  lwb(get_box(cc, i))
          hic =  upb(get_box(cc, i))
   
          lof =  lwb(get_box(cf, i))

          select case (cf%dim)
          case (1)
             call crse_coeffs_mc_1d(ccp(:,1,1,:), cfp(:,1,1,:), ng, loc, hic, lof, nc_edge)
          case (2)
             call crse_coeffs_mc_2d(ccp(:,:,1,:), cfp(:,:,1,:), ng, loc, hic, lof, nc_edge)
          case (3)
             call crse_coeffs_mc_3d(ccp(:,:,:,:), cfp(:,:,:,:), ng, loc, hic, lof, nc_edge)
          end select
       end do

    end if
    
  end subroutine coarsen_coeffs

  subroutine crse_coeffs_1d(cc, cf, ng, loc, hic, lof)

    integer        , intent(in   ) :: ng,lof(:),loc(:),hic(:)
    real(kind=dp_t), intent(inout) :: cc(loc(1)-ng:,0:)
    real(kind=dp_t), intent(in)    :: cf(lof(1)-ng:,0:)
    integer :: i, i2

    !   Cell-centered alpha array
    do i = loc(1),hic(1)
       i2 = 2*i
       cc(i,0) = HALF * (cf(i2,0) + cf(i2+1,0))
    end do

    !   Edge-centered betax array (if it exists)
    if (size(cc,dim=2) > 1) then
      do i = loc(1),hic(1)+1
         i2 = 2*i
         cc(i,1) = cf(i2,1)
      end do
    end if

  end subroutine crse_coeffs_1d

  subroutine crse_coeffs_2d(cc, cf, ng, loc, hic, lof)

    integer        , intent(in   ) :: ng,lof(:),loc(:),hic(:)
    real(kind=dp_t), intent(inout) :: cc(loc(1)-ng:,loc(2)-ng:,0:)
    real(kind=dp_t), intent(in)    :: cf(lof(1)-ng:,lof(2)-ng:,0:)
    integer :: i, i2, j, j2

    !   Cell-centered alpha array
    do j = loc(2),hic(2)
       do i = loc(1),hic(1)
          i2 = 2*i
          j2 = 2*j
          cc(i,j,0) = FOURTH * ( &
               + cf(i2,j2  ,0) + cf(i2+1,j2  ,0) & 
               + cf(i2,j2+1,0) + cf(i2+1,j2+1,0) &
               )
       end do
    end do

   !   Edge-centered betax array (if it exists)
    if (size(cc,dim=3) > 1) then
      do j = loc(2),hic(2)
         do i = loc(1),hic(1)+1
            i2 = 2*i
            j2 = 2*j
            cc(i,j,1) = HALF * (cf(i2,j2,1) + cf(i2,j2+1,1))
         end do
      end do
    end if

    !   Edge-centered betay array (if it exists)
    if (size(cc,dim=3) > 2) then
      do j = loc(2),hic(2)+1
         do i = loc(1),hic(1)
            i2 = 2*i
            j2 = 2*j
            cc(i,j,2) = HALF * (cf(i2,j2,2) + cf(i2+1,j2,2))
         end do
      end do
    end if

  end subroutine crse_coeffs_2d

  subroutine crse_coeffs_3d(cc, cf, ng, loc, hic, lof)

    integer        , intent(in   ) :: ng,lof(:),loc(:),hic(:)
    real(kind=dp_t), intent(inout) :: cc(loc(1)-ng:,loc(2)-ng:,loc(3)-ng:,0:)
    real(kind=dp_t), intent(in)    :: cf(lof(1)-ng:,lof(2)-ng:,lof(3)-ng:,0:)
    integer :: i, i2, j, j2, k, k2

    !   Cell-centered alpha array
    do k = loc(3),hic(3)
       do j = loc(2),hic(2)
          do i = loc(1),hic(1)
             i2 = 2*i
             j2 = 2*j
             k2 = 2*k
             cc(i,j,k,0) = EIGHTH * ( &
                  + cf(i2,j2  ,k2  ,0) + cf(i2+1,j2  ,k2  ,0) & 
                  + cf(i2,j2+1,k2  ,0) + cf(i2+1,j2+1,k2  ,0) &
                  + cf(i2,j2  ,k2+1,0) + cf(i2+1,j2  ,k2+1,0) &
                  + cf(i2,j2+1,k2+1,0) + cf(i2+1,j2+1,k2+1,0) &
                  )
          end do
       end do
    end do

    ! Edge-centered betax array (if it exists)
    if (size(cc,dim=4) > 1) then
      do k = loc(3),hic(3)
         do j = loc(2),hic(2)
            do i = loc(1),hic(1)+1
               i2 = 2*i
               j2 = 2*j
               k2 = 2*k
               cc(i,j,k,1) = FOURTH * ( &
                    + cf(i2,j2,k2  ,1) + cf(i2,j2+1,k2  ,1) &
                    + cf(i2,j2,k2+1,1) + cf(i2,j2+1,k2+1,1) &
                    )
            end do
         end do
      end do
    end if

    ! Edge-centered betay array (if it exists)
    if (size(cc,dim=4) > 2) then
      do k = loc(3),hic(3)
         do j = loc(2),hic(2)+1
            do i = loc(1),hic(1)
               i2 = 2*i
               j2 = 2*j
               k2 = 2*k
               cc(i,j,k,2) = FOURTH * ( &
                    + cf(i2,j2,k2  ,2) + cf(i2+1,j2,k2  ,2) &
                    + cf(i2,j2,k2+1,2) + cf(i2+1,j2,k2+1,2) &
                    )
            end do
         end do
      end do
    end if

    ! Edge-centered betaz array (if it exists)
    if (size(cc,dim=4) > 3) then
      do k = loc(3),hic(3)+1
         do j = loc(2),hic(2)
            do i = loc(1),hic(1)
               i2 = 2*i
               j2 = 2*j
               k2 = 2*k
               cc(i,j,k,3) = FOURTH * ( &
                    + cf(i2,j2  ,k2,3) + cf(i2+1,j2  ,k2,3) &
                    + cf(i2,j2+1,k2,3) + cf(i2+1,j2+1,k2,3) &
                    )
            end do
         end do
      end do
    end if

  end subroutine crse_coeffs_3d

  subroutine crse_coeffs_mc_1d(cc, cf, ng, loc, hic, lof, nc)

    integer        , intent(in   ) :: ng,lof(:),loc(:),hic(:),nc
    real(kind=dp_t), intent(inout) :: cc(loc(1)-ng:,0:)
    real(kind=dp_t), intent(in)    :: cf(lof(1)-ng:,0:)
    integer :: i, i2

    ! This is our special coarsening for the porous media solve, 
    !   where coeffs contains (dm*nc+1) components --
    !   in 1d
    !   component       0       is filled with alpha, 
    !   components      1:  nc are filled with cell-centered beta0(1:nc),
    !   components   nc+1:2*nc are filled with edge-centered betax(1:nc),

    !   Cell-centered alpha and beta0 arrays
    do i = loc(1),hic(1)
       i2 = 2*i
       cc(i,0:nc) = HALF * (cf(i2,0:nc) + cf(i2+1,0:nc))
    end do

    !   Edge-centered betax array
    do i = loc(1),hic(1)+1
       i2 = 2*i
       cc(i,nc+1:2*nc) = cf(i2,nc+1:2*nc)
    end do

  end subroutine crse_coeffs_mc_1d

  subroutine crse_coeffs_mc_2d(cc, cf, ng, loc, hic, lof, nc)

    integer        , intent(in   ) :: ng,lof(:),loc(:),hic(:),nc
    real(kind=dp_t), intent(inout) :: cc(loc(1)-ng:,loc(2)-ng:,0:)
    real(kind=dp_t), intent(in)    :: cf(lof(1)-ng:,lof(2)-ng:,0:)
    integer :: i, i2, j, j2

    ! This is our special coarsening for the porous media solve, 
    !   where coeffs contains (dm*nc+1) components --
    !   in 2d
    !   component       0       is filled with alpha, 
    !   components      1:  nc are filled with cell-centered beta0(1:nc),
    !   components   nc+1:2*nc are filled with edge-centered betax(1:nc),
    !   components 2*nc+1:3*nc are filled with edge-centered betay(1:nc),

    !   Cell-centered alpha and beta0 arrays
    do j = loc(2),hic(2)
       do i = loc(1),hic(1)
          i2 = 2*i
          j2 = 2*j
          cc(i,j,0:nc) = FOURTH * ( &
               + cf(i2,j2  ,0:nc) + cf(i2+1,j2  ,0:nc) & 
               + cf(i2,j2+1,0:nc) + cf(i2+1,j2+1,0:nc) &
               )
       end do
    end do

    !   Edge-centered betax array
    do j = loc(2),hic(2)
       do i = loc(1),hic(1)+1
          i2 = 2*i
          j2 = 2*j
          cc(i,j,nc+1:2*nc) = HALF * (cf(i2,j2,nc+1:2*nc) + cf(i2,j2+1,nc+1:2*nc))
       end do
    end do

    !   Edge-centered betay array
    do j = loc(2),hic(2)+1
       do i = loc(1),hic(1)
          i2 = 2*i
        j2 = 2*j
          cc(i,j,2*nc+1:3*nc) = HALF * (cf(i2,j2,2*nc+1:3*nc) + cf(i2+1,j2,2*nc+1:3*nc))
       end do
    end do

  end subroutine crse_coeffs_mc_2d

  subroutine crse_coeffs_mc_3d(cc, cf, ng, loc, hic, lof, nc)

    integer        , intent(in   ) :: ng,lof(:),loc(:),hic(:), nc
    real(kind=dp_t), intent(inout) :: cc(loc(1)-ng:,loc(2)-ng:,loc(3)-ng:,0:)
    real(kind=dp_t), intent(in)    :: cf(lof(1)-ng:,lof(2)-ng:,lof(3)-ng:,0:)
    integer :: i, i2, j, j2, k, k2

    ! This is our special coarsening for the porous media solve, 
    !   where coeffs contains (dm*nc+1) components --
    !   in 3d
    !   component       0       is filled with alpha, 
    !   components      1:  nc are filled with cell-centered beta0(1:nc),
    !   components   nc+1:2*nc are filled with edge-centered betax(1:nc),
    !   components 2*nc+1:3*nc are filled with edge-centered betay(1:nc),
    !   components 3*nc+1:4*nc are filled with edge-centered betaz(1:nc)

    !   Cell-centered alpha and beta0 arrays
    do k = loc(3),hic(3)
       do j = loc(2),hic(2)
          do i = loc(1),hic(1)
             i2 = 2*i
             j2 = 2*j
             k2 = 2*k
             cc(i,j,k,0:nc) = EIGHTH * ( &
                  + cf(i2,j2  ,k2  ,0:nc) + cf(i2+1,j2  ,k2  ,0:nc) & 
                  + cf(i2,j2+1,k2  ,0:nc) + cf(i2+1,j2+1,k2  ,0:nc) &
                  + cf(i2,j2  ,k2+1,0:nc) + cf(i2+1,j2  ,k2+1,0:nc) &
                  + cf(i2,j2+1,k2+1,0:nc) + cf(i2+1,j2+1,k2+1,0:nc) &
                  )
          end do
       end do
    end do

   !   Edge-centered betax array
   do k = loc(3),hic(3)
      do j = loc(2),hic(2)
         do i = loc(1),hic(1)+1
            i2 = 2*i
            j2 = 2*j
            k2 = 2*k
            cc(i,j,k,1:nc) = FOURTH * ( &
                 + cf(i2,j2,k2  ,nc+1:2*nc) + cf(i2,j2+1,k2  ,nc+1:2*nc) &
                 + cf(i2,j2,k2+1,nc+1:2*nc) + cf(i2,j2+1,k2+1,nc+1:2*nc) )
         end do
      end do
   end do

   !   Edge-centered betay array
   do k = loc(3),hic(3)
      do j = loc(2),hic(2)+1
         do i = loc(1),hic(1)
            i2 = 2*i
            j2 = 2*j
            k2 = 2*k
            cc(i,j,k,nc+1:2*nc) = FOURTH * ( &
                 + cf(i2,j2,k2  ,2*nc+1:3*nc) + cf(i2+1,j2,k2  ,2*nc+1:3*nc) &
                 + cf(i2,j2,k2+1,2*nc+1:3*nc) + cf(i2+1,j2,k2+1,2*nc+1:3*nc) )
         end do
      end do
   end do

   !   Edge-centered betaz array
   do k = loc(3),hic(3)+1
      do j = loc(2),hic(2)
         do i = loc(1),hic(1)
            i2 = 2*i
            j2 = 2*j
            k2 = 2*k
            cc(i,j,k,2*nc+1:3*nc) = FOURTH * ( &
                 + cf(i2,j2  ,k2,3*nc+1:4*nc) + cf(i2+1,j2  ,k2,3*nc+1:4*nc) &
                 + cf(i2,j2+1,k2,3*nc+1:4*nc) + cf(i2+1,j2+1,k2,3*nc+1:4*nc) )
         end do
      end do
   end do

  end subroutine crse_coeffs_mc_3d

end module coarsen_coeffs_module
