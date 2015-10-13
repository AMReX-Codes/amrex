module coarsen_coeffs_module

  use bl_constants_module
  use bl_types
  use multifab_module

  implicit none

  private
  public :: coarsen_cell_coeffs, coarsen_edge_coeffs

contains

  subroutine coarsen_cell_coeffs(cf,cc)

    type(multifab), intent(in   ) :: cf
    type(multifab), intent(inout) :: cc

    real(kind=dp_t), pointer :: cfp(:,:,:,:)
    real(kind=dp_t), pointer :: ccp(:,:,:,:)
    integer :: i,ng,dm
    integer :: lof(get_dim(cf))
    integer :: loc(get_dim(cc)),hic(get_dim(cc))

    dm = get_dim(cf)

    if ( ncomp(cf) .ne. ncomp(cc) ) then
       print *,'ncomp_fine not equal to ncomp_crse in coarsen_cell_coeffs'
       print *,'ncomp_fine = ',ncomp(cf)
       print *,'ncomp_crse = ',ncomp(cc)
       call bl_error("coarsen_coeffs.f90 :: coarsen_cell_coeffs")
    end if

    ng = nghost(cc)

    !$OMP PARALLEL DO PRIVATE(i,cfp,ccp,loc,hic,lof)
    do i = 1, nfabs(cf)
       cfp => dataptr(cf, i)
       ccp => dataptr(cc, i)

       loc =  lwb(get_box(cc, i))
       hic =  upb(get_box(cc, i))

       lof =  lwb(get_box(cf, i))

       select case (dm)
       case (1)
          call crse_cell_coeffs_1d(ccp(:,1,1,:), cfp(:,1,1,:), ng, loc, hic, lof)
       case (2)
          call crse_cell_coeffs_2d(ccp(:,:,1,:), cfp(:,:,1,:), ng, loc, hic, lof)
       case (3)
          call crse_cell_coeffs_3d(ccp(:,:,:,:), cfp(:,:,:,:), ng, loc, hic, lof)
       end select
    end do
    !$OMP END PARALLEL DO

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
    integer :: lof(get_dim(cf(1)))
    integer :: loc(get_dim(cc(1))), hic(get_dim(cc(1)))

    dm = get_dim(cf(1))

    if ( ncomp(cf(1)) .ne. ncomp(cc(1)) ) then
       print *,'ncomp_fine not equal to ncomp_crse in coarsen_edge_coeffs'
       print *,'ncomp_fine = ',ncomp(cf(1))
       print *,'ncomp_crse = ',ncomp(cc(1))
       call bl_error("coarsen_coeffs.f90 :: coarsen_edge_coeffs")
    end if

    ng = nghost(cc(1))

    !$OMP PARALLEL DO PRIVATE(i,loc,hic,lof,cxp,cyp,czp,fxp,fyp,fzp)
    do i = 1, nfabs(cf(1))

       select case (dm)
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
          call crse_xedge_coeffs_3d(cxp, fxp, ng, loc, hic, lof)
          loc =  lwb(get_box(cc(2), i))
          hic =  upb(get_box(cc(2), i))
          lof =  lwb(get_box(cf(2), i))
          cyp => dataptr(cc(2), i) 
          fyp => dataptr(cf(2), i) 
          call crse_yedge_coeffs_3d(cyp, fyp, ng, loc, hic, lof)
          loc =  lwb(get_box(cc(3), i))
          hic =  upb(get_box(cc(3), i))
          lof =  lwb(get_box(cf(3), i))
          czp => dataptr(cc(3), i) 
          fzp => dataptr(cf(3), i) 
          call crse_zedge_coeffs_3d(czp, fzp, ng, loc, hic, lof)
       end select
    end do 
    !$OMP END PARALLEL DO

  end subroutine coarsen_edge_coeffs

  subroutine crse_cell_coeffs_1d(cc, cf, ng, loc, hic, lof)

    integer        , intent(in   ) :: ng,lof(:),loc(:),hic(:)
    real(kind=dp_t), intent(inout) :: cc(loc(1)-ng:,:)
    real(kind=dp_t), intent(in)    :: cf(lof(1)-ng:,:)
    integer :: i, i2, n

    do n = 1, size(cc,2)
       do i = loc(1),hic(1)
          i2 = 2*i
          cc(i,n) = HALF * (cf(i2,n) + cf(i2+1,n))
       end do
    end do

  end subroutine crse_cell_coeffs_1d

  subroutine crse_xedge_coeffs_1d(cc, cf, ng, loc, hic, lof)

    integer        , intent(in   ) :: ng,lof(:),loc(:),hic(:)
    real(kind=dp_t), intent(inout) :: cc(loc(1)-ng:,:)
    real(kind=dp_t), intent(in)    :: cf(lof(1)-ng:,:)
    integer :: i, i2, n

    do n = 1, size(cc,2)
       do i = loc(1),hic(1)+1
          i2 = 2*i
          cc(i,n) = cf(i2,n)
       end do
    end do

  end subroutine crse_xedge_coeffs_1d

  subroutine crse_cell_coeffs_2d(cc, cf, ng, loc, hic, lof)

    integer        , intent(in   ) :: ng,lof(:),loc(:),hic(:)
    real(kind=dp_t), intent(inout) :: cc(loc(1)-ng:,loc(2)-ng:,:)
    real(kind=dp_t), intent(in)    :: cf(lof(1)-ng:,lof(2)-ng:,:)
    integer :: i, i2, j, j2, n

    do n = 1, size(cc,3)
       do j = loc(2),hic(2)
          do i = loc(1),hic(1)
             i2 = 2*i
             j2 = 2*j
             cc(i,j,n) = FOURTH * (  &
                  cf(i2,j2  ,n) + cf(i2+1,j2  ,n) + & 
                  cf(i2,j2+1,n) + cf(i2+1,j2+1,n) )
          end do
       end do
    end do

  end subroutine crse_cell_coeffs_2d

  subroutine crse_xedge_coeffs_2d(cc, cf, ng, loc, hic, lof)

    integer        , intent(in   ) :: ng,lof(:),loc(:),hic(:)
    real(kind=dp_t), intent(inout) :: cc(loc(1)-ng:,loc(2)-ng:,:)
    real(kind=dp_t), intent(in)    :: cf(lof(1)-ng:,lof(2)-ng:,:)
    integer :: i, i2, j, j2, n

    do n = 1, size(cc,3)
       do j = loc(2),hic(2)
          do i = loc(1),hic(1)+1
             i2 = 2*i
             j2 = 2*j
             cc(i,j,n) = HALF * (cf(i2,j2,n) + cf(i2,j2+1,n))
          end do
       end do
    end do

  end subroutine crse_xedge_coeffs_2d

  subroutine crse_yedge_coeffs_2d(cc, cf, ng, loc, hic, lof)

    integer        , intent(in   ) :: ng,lof(:),loc(:),hic(:)
    real(kind=dp_t), intent(inout) :: cc(loc(1)-ng:,loc(2)-ng:,:)
    real(kind=dp_t), intent(in)    :: cf(lof(1)-ng:,lof(2)-ng:,:)
    integer :: i, i2, j, j2, n

    do n = 1, size(cc,3)
       do j = loc(2),hic(2)+1
          do i = loc(1),hic(1)
             i2 = 2*i
             j2 = 2*j
             cc(i,j,n) = HALF * (cf(i2,j2,n) + cf(i2+1,j2,n))
          end do
       end do
    end do

  end subroutine crse_yedge_coeffs_2d

  subroutine crse_cell_coeffs_3d(cc, cf, ng, loc, hic, lof)

    integer        , intent(in   ) :: ng,lof(:),loc(:),hic(:)
    real(kind=dp_t), intent(inout) :: cc(loc(1)-ng:,loc(2)-ng:,loc(3)-ng:,1:)
    real(kind=dp_t), intent(in)    :: cf(lof(1)-ng:,lof(2)-ng:,lof(3)-ng:,1:)
    integer :: i, i2, j, j2, k, k2, n

    do n = 1, size(cc,4)
       do k = loc(3),hic(3)
          k2 = 2*k
          do j = loc(2),hic(2)
             j2 = 2*j
             do i = loc(1),hic(1)
                i2 = 2*i
                cc(i,j,k,n) = EIGHTH * ( &
                     + cf(i2,j2  ,k2  ,n) + cf(i2+1,j2  ,k2  ,n) & 
                     + cf(i2,j2+1,k2  ,n) + cf(i2+1,j2+1,k2  ,n) &
                     + cf(i2,j2  ,k2+1,n) + cf(i2+1,j2  ,k2+1,n) &
                     + cf(i2,j2+1,k2+1,n) + cf(i2+1,j2+1,k2+1,n) &
                     )
             end do
          end do
       end do
    end do

  end subroutine crse_cell_coeffs_3d

  subroutine crse_xedge_coeffs_3d(cc, cf, ng, loc, hic, lof)

    integer        , intent(in   ) :: ng,lof(:),loc(:),hic(:)
    real(kind=dp_t), intent(inout) :: cc(loc(1)-ng:,loc(2)-ng:,loc(3)-ng:,1:)
    real(kind=dp_t), intent(in)    :: cf(lof(1)-ng:,lof(2)-ng:,lof(3)-ng:,1:)
    integer :: i, i2, j, j2, k, k2, n

    do n = 1, size(cc,4)
       do k = loc(3),hic(3)
          k2 = 2*k
          do j = loc(2),hic(2)
             j2 = 2*j
             do i = loc(1),hic(1)+1
                i2 = 2*i
                cc(i,j,k,n) = FOURTH * ( &
                     + cf(i2,j2,k2  ,n) + cf(i2,j2+1,k2  ,n) &
                     + cf(i2,j2,k2+1,n) + cf(i2,j2+1,k2+1,n) &
                     )
             end do
          end do
       end do
    end do

  end subroutine crse_xedge_coeffs_3d

  subroutine crse_yedge_coeffs_3d(cc, cf, ng, loc, hic, lof)

    integer        , intent(in   ) :: ng,lof(:),loc(:),hic(:)
    real(kind=dp_t), intent(inout) :: cc(loc(1)-ng:,loc(2)-ng:,loc(3)-ng:,1:)
    real(kind=dp_t), intent(in)    :: cf(lof(1)-ng:,lof(2)-ng:,lof(3)-ng:,1:)
    integer :: i, i2, j, j2, k, k2, n

    do n = 1, size(cc,4)
       do k = loc(3),hic(3)
          k2 = 2*k
          do j = loc(2),hic(2)+1
             j2 = 2*j
             do i = loc(1),hic(1)
                i2 = 2*i
                cc(i,j,k,n) = FOURTH * ( &
                     + cf(i2,j2,k2  ,n) + cf(i2+1,j2,k2  ,n) &
                     + cf(i2,j2,k2+1,n) + cf(i2+1,j2,k2+1,n) &
                     )
             end do
          end do
       end do
    end do

  end subroutine crse_yedge_coeffs_3d

  subroutine crse_zedge_coeffs_3d(cc, cf, ng, loc, hic, lof)

    integer        , intent(in   ) :: ng,lof(:),loc(:),hic(:)
    real(kind=dp_t), intent(inout) :: cc(loc(1)-ng:,loc(2)-ng:,loc(3)-ng:,1:)
    real(kind=dp_t), intent(in)    :: cf(lof(1)-ng:,lof(2)-ng:,lof(3)-ng:,1:)
    integer :: i, i2, j, j2, k, k2, n

    do n = 1, size(cc,4)
       do k = loc(3),hic(3)+1
          k2 = 2*k
          do j = loc(2),hic(2)
             j2 = 2*j
             do i = loc(1),hic(1)
                i2 = 2*i
                cc(i,j,k,n) = FOURTH * ( &
                     + cf(i2,j2  ,k2,n) + cf(i2+1,j2  ,k2,n) &
                     + cf(i2,j2+1,k2,n) + cf(i2+1,j2+1,k2,n) &
                     )
             end do
          end do
       end do
    end do

  end subroutine crse_zedge_coeffs_3d

end module coarsen_coeffs_module
