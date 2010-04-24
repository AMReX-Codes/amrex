module coarsen_coeffs_module

  use bl_types
  use multifab_module

  implicit none

  real(kind=dp_t), private, parameter :: HALF   = 0.5_dp_t
  real(kind=dp_t), private, parameter :: FOURTH = 0.25_dp_t
  real(kind=dp_t), private, parameter :: EIGHTH = 0.125_dp_t

  private :: crse_cell_coeffs_1d, crse_xedge_coeffs_1d
  private :: crse_cell_coeffs_2d, crse_xedge_coeffs_2d, crse_yedge_coeffs_2d
  private :: crse_cell_coeffs_3d, crse_xedge_coeffs_3d, crse_yedge_coeffs_3d, crse_zedge_coeffs_3d

contains

  subroutine coarsen_cell_coeffs(cf,cc)

    type(multifab), intent(in   ) :: cf
    type(multifab), intent(inout) :: cc

    real(kind=dp_t), pointer :: cfp(:,:,:,:)
    real(kind=dp_t), pointer :: ccp(:,:,:,:)
    integer :: i,ng,dm
    integer :: lof(cf%dim)
    integer :: loc(cc%dim),hic(cc%dim)

    dm = cf%dim

    if (cf%nc .ne. cc%nc) then
       print *,'ncomp_fine not equal to ncomp_crse in coarsen_cell_coeffs'
       print *,'ncomp_fine = ',cf%nc
       print *,'ncomp_crse = ',cc%nc
       call bl_error("coarsen_coeffs.f90 :: coarsen_cell_coeffs")
    end if

    ng = cc%ng

    do i = 1, cf%nboxes
       if ( multifab_remote(cf,i) ) cycle
       cfp => dataptr(cf, i)
       ccp => dataptr(cc, i)

       loc =  lwb(get_box(cc, i))
       hic =  upb(get_box(cc, i))

       lof =  lwb(get_box(cf, i))

       select case (cf%dim)
       case (1)
          call crse_cell_coeffs_1d(ccp(:,1,1,:), cfp(:,1,1,:), ng, loc, hic, lof)
       case (2)
          call crse_cell_coeffs_2d(ccp(:,:,1,:), cfp(:,:,1,:), ng, loc, hic, lof)
       case (3)
          call crse_cell_coeffs_3d(ccp(:,:,:,:), cfp(:,:,:,:), ng, loc, hic, lof)
       end select
    end do

  end subroutine coarsen_cell_coeffs

  subroutine coarsen_edge_coeffs(cf,cc)

    type(multifab), intent(in   ) :: cf(:)
    type(multifab), intent(inout) :: cc(:)

    real(kind=dp_t), pointer :: cxp(:,:,:,:)
    real(kind=dp_t), pointer :: cyp(:,:,:,:)
    real(kind=dp_t), pointer :: czp(:,:,:,:)
    real(kind=dp_t), pointer :: fxp(:,:,:,:)
    real(kind=dp_t), pointer :: fyp(:,:,:,:)
    real(kind=dp_t), pointer :: fzp(:,:,:,:)

    integer :: i,ng,dm
    integer :: lof(cf(1)%dim)
    integer :: loc(cc(1)%dim),hic(cc(1)%dim)

    dm = cf(1)%dim

    if (cf(1)%nc .ne. cc(1)%nc) then
       print *,'ncomp_fine not equal to ncomp_crse in coarsen_edge_coeffs'
       print *,'ncomp_fine = ',cf(1)%nc
       print *,'ncomp_crse = ',cc(1)%nc
       call bl_error("coarsen_coeffs.f90 :: coarsen_edge_coeffs")
    end if

    ng = cc(1)%ng

    do i = 1, cf(1)%nboxes
       if ( multifab_remote(cf(1),i) ) cycle

       select case (cf(1)%dim)
       case (1)
          loc =  lwb(get_box(cc(1), i))
          hic =  upb(get_box(cc(1), i))
          lof =  lwb(get_box(cf(1), i))
          cxp => dataptr(cc(1), i) 
          fxp => dataptr(cf(1), i) 
          call crse_xedge_coeffs_1d(cxp(:,1,1,:), fxp(:,1,1,:), ng, loc, hic, lof)
       case (2)
          loc =  lwb(get_box(cc(1), i))
          hic =  upb(get_box(cc(1), i))
          lof =  lwb(get_box(cf(1), i))
          cxp => dataptr(cc(1), i) 
          fxp => dataptr(cf(1), i) 
          call crse_xedge_coeffs_2d(cxp(:,:,1,:), fxp(:,:,1,:), ng, loc, hic, lof)
          loc =  lwb(get_box(cc(2), i))
          hic =  upb(get_box(cc(2), i))
          lof =  lwb(get_box(cf(2), i))
          cyp => dataptr(cc(2), i) 
          fyp => dataptr(cf(2), i) 
          call crse_yedge_coeffs_2d(cyp(:,:,1,:), fyp(:,:,1,:), ng, loc, hic, lof)
       case (3)
          loc =  lwb(get_box(cc(1), i))
          hic =  upb(get_box(cc(1), i))
          lof =  lwb(get_box(cf(1), i))
          cxp => dataptr(cc(1), i) 
          fxp => dataptr(cf(1), i) 
          call crse_xedge_coeffs_3d(cxp(:,:,:,:), fxp(:,:,:,:), ng, loc, hic, lof)
          loc =  lwb(get_box(cc(2), i))
          hic =  upb(get_box(cc(2), i))
          lof =  lwb(get_box(cf(2), i))
          cyp => dataptr(cc(2), i) 
          fyp => dataptr(cf(2), i) 
          call crse_yedge_coeffs_3d(cyp(:,:,:,:), fyp(:,:,:,:), ng, loc, hic, lof)
          loc =  lwb(get_box(cc(3), i))
          hic =  upb(get_box(cc(3), i))
          lof =  lwb(get_box(cf(3), i))
          czp => dataptr(cc(3), i) 
          fzp => dataptr(cf(3), i) 
          call crse_zedge_coeffs_3d(czp(:,:,:,:), fzp(:,:,:,:), ng, loc, hic, lof)
       end select
    end do 

  end subroutine coarsen_edge_coeffs

  subroutine crse_cell_coeffs_1d(cc, cf, ng, loc, hic, lof)

    integer        , intent(in   ) :: ng,lof(:),loc(:),hic(:)
    real(kind=dp_t), intent(inout) :: cc(loc(1)-ng:,:)
    real(kind=dp_t), intent(in)    :: cf(lof(1)-ng:,:)
    integer :: i, i2

    do i = loc(1),hic(1)
       i2 = 2*i
       cc(i,:) = HALF * (cf(i2,:) + cf(i2+1,:))
    end do

  end subroutine crse_cell_coeffs_1d

  subroutine crse_xedge_coeffs_1d(cc, cf, ng, loc, hic, lof)

    integer        , intent(in   ) :: ng,lof(:),loc(:),hic(:)
    real(kind=dp_t), intent(inout) :: cc(loc(1)-ng:,:)
    real(kind=dp_t), intent(in)    :: cf(lof(1)-ng:,:)
    integer :: i, i2

   do i = loc(1),hic(1)+1
      i2 = 2*i
      cc(i,:) = cf(i2,:)
   end do

  end subroutine crse_xedge_coeffs_1d

  subroutine crse_cell_coeffs_2d(cc, cf, ng, loc, hic, lof)

    integer        , intent(in   ) :: ng,lof(:),loc(:),hic(:)
    real(kind=dp_t), intent(inout) :: cc(loc(1)-ng:,loc(2)-ng:,:)
    real(kind=dp_t), intent(in)    :: cf(lof(1)-ng:,lof(2)-ng:,:)
    integer :: i, i2, j, j2

    do j = loc(2),hic(2)
       do i = loc(1),hic(1)
          i2 = 2*i
          j2 = 2*j
          cc(i,j,:) = FOURTH * ( &
               + cf(i2,j2  ,:) + cf(i2+1,j2  ,:) & 
               + cf(i2,j2+1,:) + cf(i2+1,j2+1,:) &
               )
       end do
    end do

  end subroutine crse_cell_coeffs_2d

  subroutine crse_xedge_coeffs_2d(cc, cf, ng, loc, hic, lof)

    integer        , intent(in   ) :: ng,lof(:),loc(:),hic(:)
    real(kind=dp_t), intent(inout) :: cc(loc(1)-ng:,loc(2)-ng:,:)
    real(kind=dp_t), intent(in)    :: cf(lof(1)-ng:,lof(2)-ng:,:)
    integer :: i, i2, j, j2

    do j = loc(2),hic(2)
       do i = loc(1),hic(1)+1
          i2 = 2*i
          j2 = 2*j
          cc(i,j,:) = HALF * (cf(i2,j2,:) + cf(i2,j2+1,:))
        end do
    end do

  end subroutine crse_xedge_coeffs_2d

  subroutine crse_yedge_coeffs_2d(cc, cf, ng, loc, hic, lof)

    integer        , intent(in   ) :: ng,lof(:),loc(:),hic(:)
    real(kind=dp_t), intent(inout) :: cc(loc(1)-ng:,loc(2)-ng:,:)
    real(kind=dp_t), intent(in)    :: cf(lof(1)-ng:,lof(2)-ng:,:)
    integer :: i, i2, j, j2

    do j = loc(2),hic(2)+1
       do i = loc(1),hic(1)
          i2 = 2*i
          j2 = 2*j
          cc(i,j,:) = HALF * (cf(i2,j2,:) + cf(i2+1,j2,:))
       end do
    end do

  end subroutine crse_yedge_coeffs_2d

  subroutine crse_cell_coeffs_3d(cc, cf, ng, loc, hic, lof)

    integer        , intent(in   ) :: ng,lof(:),loc(:),hic(:)
    real(kind=dp_t), intent(inout) :: cc(loc(1)-ng:,loc(2)-ng:,loc(3)-ng:,:)
    real(kind=dp_t), intent(in)    :: cf(lof(1)-ng:,lof(2)-ng:,lof(3)-ng:,:)
    integer :: i, i2, j, j2, k, k2

    do k = loc(3),hic(3)
       do j = loc(2),hic(2)
          do i = loc(1),hic(1)
             i2 = 2*i
             j2 = 2*j
             k2 = 2*k
             cc(i,j,k,:) = EIGHTH * ( &
                  + cf(i2,j2  ,k2  ,:) + cf(i2+1,j2  ,k2  ,:) & 
                  + cf(i2,j2+1,k2  ,:) + cf(i2+1,j2+1,k2  ,:) &
                  + cf(i2,j2  ,k2+1,:) + cf(i2+1,j2  ,k2+1,:) &
                  + cf(i2,j2+1,k2+1,:) + cf(i2+1,j2+1,k2+1,:) &
                  )
          end do
       end do
    end do

  end subroutine crse_cell_coeffs_3d

  subroutine crse_xedge_coeffs_3d(cc, cf, ng, loc, hic, lof)

    integer        , intent(in   ) :: ng,lof(:),loc(:),hic(:)
    real(kind=dp_t), intent(inout) :: cc(loc(1)-ng:,loc(2)-ng:,loc(3)-ng:,:)
    real(kind=dp_t), intent(in)    :: cf(lof(1)-ng:,lof(2)-ng:,lof(3)-ng:,:)
    integer :: i, i2, j, j2, k, k2

    do k = loc(3),hic(3)
       do j = loc(2),hic(2)
          do i = loc(1),hic(1)+1
             i2 = 2*i
             j2 = 2*j
             k2 = 2*k
             cc(i,j,k,:) = FOURTH * ( &
                  + cf(i2,j2,k2  ,:) + cf(i2,j2+1,k2  ,:) &
                  + cf(i2,j2,k2+1,:) + cf(i2,j2+1,k2+1,:) &
                  )
          end do
       end do
    end do

  end subroutine crse_xedge_coeffs_3d

  subroutine crse_yedge_coeffs_3d(cc, cf, ng, loc, hic, lof)

    integer        , intent(in   ) :: ng,lof(:),loc(:),hic(:)
    real(kind=dp_t), intent(inout) :: cc(loc(1)-ng:,loc(2)-ng:,loc(3)-ng:,:)
    real(kind=dp_t), intent(in)    :: cf(lof(1)-ng:,lof(2)-ng:,lof(3)-ng:,:)
    integer :: i, i2, j, j2, k, k2

    do k = loc(3),hic(3)
       do j = loc(2),hic(2)+1
          do i = loc(1),hic(1)
             i2 = 2*i
             j2 = 2*j
             k2 = 2*k
             cc(i,j,k,:) = FOURTH * ( &
                  + cf(i2,j2,k2  ,:) + cf(i2+1,j2,k2  ,:) &
                  + cf(i2,j2,k2+1,:) + cf(i2+1,j2,k2+1,:) &
                  )
          end do
       end do
    end do

  end subroutine crse_yedge_coeffs_3d

  subroutine crse_zedge_coeffs_3d(cc, cf, ng, loc, hic, lof)

    integer        , intent(in   ) :: ng,lof(:),loc(:),hic(:)
    real(kind=dp_t), intent(inout) :: cc(loc(1)-ng:,loc(2)-ng:,loc(3)-ng:,:)
    real(kind=dp_t), intent(in)    :: cf(lof(1)-ng:,lof(2)-ng:,lof(3)-ng:,:)
    integer :: i, i2, j, j2, k, k2

    do k = loc(3),hic(3)+1
       do j = loc(2),hic(2)
          do i = loc(1),hic(1)
             i2 = 2*i
             j2 = 2*j
             k2 = 2*k
             cc(i,j,k,:) = FOURTH * ( &
                  + cf(i2,j2  ,k2,:) + cf(i2+1,j2  ,k2,:) &
                  + cf(i2,j2+1,k2,:) + cf(i2+1,j2+1,k2,:) &
                  )
          end do
       end do
    end do

  end subroutine crse_zedge_coeffs_3d

end module coarsen_coeffs_module
