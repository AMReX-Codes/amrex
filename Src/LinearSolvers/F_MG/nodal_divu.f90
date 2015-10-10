module nodal_divu_module

  use bl_constants_module
  use bc_functions_module
  use bndry_reg_module
  use mg_tower_module
  use ml_nodal_restriction_module
  use define_bc_module

  implicit none

  private
  public :: divu, divucc, subtract_divu_from_rh

contains

    subroutine divu(nlevs,mgt,unew,rh,ref_ratio,nodal, lo_inflow, hi_inflow)

      integer        , intent(in   ) :: nlevs
      type(mg_tower) , intent(inout) :: mgt(:)
      type(multifab) , intent(inout) :: unew(:)
      type(multifab) , intent(inout) :: rh(:)
      integer        , intent(in   ) :: ref_ratio(:,:)
      logical        , intent(in   ) :: nodal(:)
      integer        , intent(in   ) :: lo_inflow(:), hi_inflow(:)

      real(kind=dp_t), pointer :: unp(:,:,:,:) 
      real(kind=dp_t), pointer :: rhp(:,:,:,:) 
      integer        , pointer ::  mp(:,:,:,:) 

      integer         :: i,n,dm,ng_u,ng_r,mglev_fine
      integer         :: lo(get_dim(unew(nlevs))),hi(get_dim(unew(nlevs)))
      type(      box) :: pdc
      type(   layout) :: la_crse,la_fine
      type(bndry_reg) :: brs_flx

      dm = get_dim(unew(nlevs))
      ng_u = nghost(unew(nlevs))
      ng_r = nghost(rh(nlevs))

!     Create the regular single-level divergence.
      do n = 1, nlevs
         mglev_fine = mgt(n)%nlevels
         call multifab_fill_boundary(unew(n))
         do i = 1, nfabs(unew(n))
            unp => dataptr(unew(n), i)
            rhp => dataptr(rh(n)  , i)
            mp  => dataptr(mgt(n)%mm(mglev_fine),i)
            lo = lwb(get_box(unew(n),i))
            hi = upb(get_box(unew(n),i))
            select case (dm)
               case (1)
                 call divu_1d(unp(:,1,1,1), rhp(:,1,1,1), &
                               mp(:,1,1,1), mgt(n)%dh(:,mglev_fine), &
                              mgt(n)%face_type(i,:,:), lo, hi, ng_u, ng_r)
               case (2)
                 call divu_2d(unp(:,:,1,:), rhp(:,:,1,1), &
                               mp(:,:,1,1), mgt(n)%dh(:,mglev_fine), &
                              mgt(n)%face_type(i,:,:), lo, hi, ng_u, ng_r)
               case (3)
                 call divu_3d(unp(:,:,:,:), rhp(:,:,:,1), &
                               mp(:,:,:,1), mgt(n)%dh(:,mglev_fine), &
                              mgt(n)%face_type(i,:,:), lo, hi, ng_u, ng_r)
            end select
         end do
      end do

!     Modify the divu above at coarse-fine interfaces.
      do n = nlevs,2,-1

         la_fine = get_layout(unew(n))
         la_crse = get_layout(unew(n-1))
         pdc = get_pd(la_crse)

         call bndry_reg_rr_build_nd(brs_flx,la_fine, ref_ratio(n-1,:), pdc, nodal)
         call crse_fine_divu(n,nlevs,rh(n-1),unew,brs_flx,ref_ratio(n-1,:),mgt, &
              lo_inflow, hi_inflow)
         call bndry_reg_destroy(brs_flx)
      end do

    end subroutine divu

    subroutine divu_1d(u,rh,mm,dx,face_type,lo,hi,ng_u,ng_r)

      integer        , intent(in   ) :: ng_u,ng_r,lo(:),hi(:)
      real(kind=dp_t), intent(inout) ::  u(lo(1)-ng_u:)
      real(kind=dp_t), intent(inout) :: rh(lo(1)-ng_r:)
      integer        , intent(inout) :: mm(lo(1):)
      real(kind=dp_t), intent(in   ) :: dx(:)
      integer        , intent(in   ) :: face_type(:,:)

      integer         :: i

      rh = ZERO

      do i = lo(1),hi(1)+1
         if (.not.bc_dirichlet(mm(i),1,0)) then 
           rh(i) = (u(i)-u(i-1)) / dx(1) 
         end if
      end do

      if (face_type(1,1) == BC_NEU) rh( 0) = TWO*rh( 0)
      if (face_type(1,2) == BC_NEU) rh(hi(1)+1) = TWO*rh(hi(1)+1)

    end subroutine divu_1d

    subroutine divu_2d(u,rh,mm,dx,face_type,lo,hi,ng_u,ng_r)

      integer        , intent(in   ) :: ng_u,ng_r,lo(:),hi(:)
      real(kind=dp_t), intent(inout) ::  u(lo(1)-ng_u:,lo(2)-ng_u:,:)
      real(kind=dp_t), intent(inout) :: rh(lo(1)-ng_r:,lo(2)-ng_r:)
      integer        , intent(inout) :: mm(lo(1):,lo(2):)
      real(kind=dp_t), intent(in   ) :: dx(:)
      integer        , intent(in   ) :: face_type(:,:)

      integer         :: i,j

      rh = ZERO

      do j = lo(2),hi(2)+1
      do i = lo(1),hi(1)+1
         if (.not.bc_dirichlet(mm(i,j),1,0)) then 
           rh(i,j) = (u(i  ,j,1) + u(i  ,j-1,1) &
                     -u(i-1,j,1) - u(i-1,j-1,1)) / dx(1) + &
                     (u(i,j  ,2) + u(i-1,j  ,2) &
                     -u(i,j-1,2) - u(i-1,j-1,2)) / dx(2)
           rh(i,j) = HALF * rh(i,j)
         end if
      end do
      end do

      if (face_type(1,1) == BC_NEU) rh( 0,:) = TWO*rh( 0,:)
      if (face_type(1,2) == BC_NEU) rh(hi(1)+1,:) = TWO*rh(hi(1)+1,:)
      if (face_type(2,1) == BC_NEU) rh(:, 0) = TWO*rh(:, 0)
      if (face_type(2,2) == BC_NEU) rh(:,hi(2)+1) = TWO*rh(:,hi(2)+1)

    end subroutine divu_2d

    subroutine divu_3d(u,rh,mm,dx,face_type,lo,hi,ng_u,ng_r)

      integer        , intent(in   ) :: ng_u,ng_r,lo(:),hi(:)
      real(kind=dp_t), intent(inout) ::  u(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:,:)
      real(kind=dp_t), intent(inout) :: rh(lo(1)-ng_r:,lo(2)-ng_r:,lo(3)-ng_r:)
      integer        , intent(inout) :: mm(lo(1):,lo(2):,lo(3):)
      real(kind=dp_t), intent(in   ) :: dx(:)
      integer        , intent(in   ) :: face_type(:,:)

      integer         :: i,j,k
      real(kind=dp_t) :: ivdx, ivdy, ivdz
      
      rh = ZERO

      ivdx = ONE / dx(1)
      ivdy = ONE / dx(2)
      ivdz = ONE / dx(3)

      !$OMP PARALLEL DO PRIVATE(i,j,k)
      do k = lo(3),hi(3)+1
      do j = lo(2),hi(2)+1
      do i = lo(1),hi(1)+1
         if (.not. bc_dirichlet(mm(i,j,k),1,0)) then
           rh(i,j,k) = (u(i  ,j,k  ,1) + u(i  ,j-1,k  ,1) &
                       +u(i  ,j,k-1,1) + u(i  ,j-1,k-1,1) &
                       -u(i-1,j,k  ,1) - u(i-1,j-1,k  ,1) &
                       -u(i-1,j,k-1,1) - u(i-1,j-1,k-1,1)) * ivdx + &
                       (u(i,j  ,k  ,2) + u(i-1,j  ,k  ,2) &
                       +u(i,j  ,k-1,2) + u(i-1,j  ,k-1,2) &
                       -u(i,j-1,k  ,2) - u(i-1,j-1,k  ,2) &
                       -u(i,j-1,k-1,2) - u(i-1,j-1,k-1,2)) * ivdy + &
                       (u(i,j  ,k  ,3) + u(i-1,j  ,k  ,3) &
                       +u(i,j-1,k  ,3) + u(i-1,j-1,k  ,3) &
                       -u(i,j  ,k-1,3) - u(i-1,j  ,k-1,3) &
                       -u(i,j-1,k-1,3) - u(i-1,j-1,k-1,3)) * ivdz
           rh(i,j,k) = FOURTH * rh(i,j,k)
         end if
      end do
      end do
      end do
      !$OMP END PARALLEL DO

      if (face_type(1,1) == BC_NEU) rh( 0,:,:) = TWO*rh( 0,:,:)
      if (face_type(1,2) == BC_NEU) rh(hi(1)+1,:,:) = TWO*rh(hi(1)+1,:,:)
      if (face_type(2,1) == BC_NEU) rh(:, 0,:) = TWO*rh(:, 0,:)
      if (face_type(2,2) == BC_NEU) rh(:,hi(2)+1,:) = TWO*rh(:,hi(2)+1,:)
      if (face_type(3,1) == BC_NEU) rh(:,:, 0) = TWO*rh(:,:, 0)
      if (face_type(3,2) == BC_NEU) rh(:,:,hi(3)+1) = TWO*rh(:,:,hi(3)+1)

    end subroutine divu_3d

    subroutine crse_fine_divu(n_fine,nlevs,rh_crse,u,brs_flx,ref_ratio,mgt, &
         lo_inflow, hi_inflow)

      use nodal_interface_stencil_module, only : ml_fine_contrib

      integer        , intent(in   ) :: n_fine,nlevs
      type(multifab) , intent(inout) :: rh_crse
      type(multifab) , intent(inout) :: u(:)
      type(bndry_reg), intent(inout) :: brs_flx
      integer        , intent(in   ) :: ref_ratio(:)
      type(mg_tower) , intent(in   ) :: mgt(:)
      integer        , intent(in   ) :: lo_inflow(:), hi_inflow(:)

      real(kind=dp_t), pointer :: unp(:,:,:,:) 
      real(kind=dp_t), pointer :: rhp(:,:,:,:) 

      type(multifab) :: temp_rhs, temp_rhs_crse
      type(  layout) :: la_crse,la_fine
      type(     box) :: pdc
      integer        :: i,dm,n_crse,ng,mglev_fine, mglev_crse
      logical        :: nodal(get_dim(u(n_fine)))

      dm     = get_dim(u(n_fine))
      n_crse = n_fine-1

      ng    = nghost(u(nlevs))
      nodal = .true.

      la_crse = get_layout(u(n_crse))
      la_fine = get_layout(u(n_fine))

      mglev_crse = mgt(n_crse)%nlevels
      mglev_fine = mgt(n_fine)%nlevels

      call multifab_build(temp_rhs, la_fine, 1, 1, nodal)
      call setval(temp_rhs, ZERO, 1, all=.true.)

!     Zero out the flux registers which will hold the fine contributions
      call bndry_reg_setval(brs_flx, ZERO, all = .true.)

!     Compute the fine contributions at faces, edges and corners.

!     First compute a residual which only takes contributions from the
!        grid on which it is calculated.
      do i = 1, nfabs(u(n_fine))
          unp => dataptr(u(n_fine), i)
          rhp => dataptr( temp_rhs, i)
          select case (dm)
             case (1)
               call grid_divu_1d(unp(:,1,1,1), rhp(:,1,1,1), mgt(n_fine)%dh(:,mglev_fine), &
                                 mgt(n_fine)%face_type(i,:,:), ng, lo_inflow, hi_inflow)
             case (2)
               call grid_divu_2d(unp(:,:,1,:), rhp(:,:,1,1), mgt(n_fine)%dh(:,mglev_fine), &
                                 mgt(n_fine)%face_type(i,:,:), ng, lo_inflow, hi_inflow)
             case (3)
               call grid_divu_3d(unp(:,:,:,:), rhp(:,:,:,1), mgt(n_fine)%dh(:,mglev_fine), &
                                 mgt(n_fine)%face_type(i,:,:), ng, lo_inflow, hi_inflow)
          end select
      end do

      pdc = get_pd(la_crse)

      do i = 1,dm
         call ml_fine_contrib(brs_flx%bmf(i,0), &
                              temp_rhs,mgt(n_fine)%mm(mglev_fine),ref_ratio,pdc,-i)
         call ml_fine_contrib(brs_flx%bmf(i,1), &
                              temp_rhs,mgt(n_fine)%mm(mglev_fine),ref_ratio,pdc,+i)
      end do

!     Compute the crse contributions at edges and corners and add to rh(n-1).

      call multifab_build(temp_rhs_crse, la_crse, 1, 1, nodal)
      call setval(temp_rhs_crse, ZERO, 1, all=.true.)

      do i = 1,dm
         call ml_crse_divu_contrib(temp_rhs_crse, brs_flx%bmf(i,0), u(n_crse), &
                                   mgt(n_fine)%mm(mglev_fine), mgt(n_crse)%dh(:,mglev_crse), &
                                   pdc,ref_ratio, -i, lo_inflow, hi_inflow)
         call ml_crse_divu_contrib(temp_rhs_crse, brs_flx%bmf(i,1), u(n_crse), &
                                   mgt(n_fine)%mm(mglev_fine), mgt(n_crse)%dh(:,mglev_crse), &
                                   pdc,ref_ratio, +i, lo_inflow, hi_inflow)
      end do

      call multifab_plus_plus(rh_crse,temp_rhs_crse)
      call periodic_add_copy(rh_crse,temp_rhs_crse,synced=.true.)

      call multifab_destroy(temp_rhs)
      call multifab_destroy(temp_rhs_crse)

    end subroutine crse_fine_divu

    subroutine grid_divu_1d(u,rh,dx,face_type,ng, lo_inflow, hi_inflow)

      integer        , intent(in   ) :: ng
      real(kind=dp_t), intent(inout) ::  u(-ng:)
      real(kind=dp_t), intent(inout) :: rh(-1:)
      real(kind=dp_t), intent(in   ) :: dx(:)
      integer        , intent(in   ) :: face_type(:,:)
      integer        , intent(in   ) :: lo_inflow(:), hi_inflow(:)

      integer :: i,nx
      nx = size(rh,dim=1) - 3

      if (face_type(1,1) .ne. BC_NEU .or. lo_inflow(1) .ne. 1) then
         u(-1) = ZERO
      end if

      if (face_type(1,2) .ne. BC_NEU .or. hi_inflow(1) .ne. 1) then
         u(nx) = ZERO
      end if

      do i = 0,nx
         rh(i) = (u(i)-u(i-1)) / dx(1)
      end do

      if (face_type(1,1) == BC_NEU) rh( 0) = TWO*rh( 0)
      if (face_type(1,2) == BC_NEU) rh(nx) = TWO*rh(nx)

    end subroutine grid_divu_1d

    subroutine grid_divu_2d(u,rh,dx,face_type,ng, lo_inflow, hi_inflow)

      integer        , intent(in   ) :: ng
      real(kind=dp_t), intent(inout) ::  u(-ng:,-ng:,1:)
      real(kind=dp_t), intent(inout) :: rh(-1:,-1:)
      real(kind=dp_t), intent(in   ) :: dx(:)
      integer        , intent(in   ) :: face_type(:,:)
      integer        , intent(in   ) :: lo_inflow(:), hi_inflow(:)

      integer :: i,j,nx,ny
      nx = size(rh,dim=1) - 3
      ny = size(rh,dim=2) - 3

      ! x-veclocity
      u(:,-1,1) = ZERO
      u(:,ny,1) = ZERO
      if (lo_inflow(1) .ne. 1 .or. face_type(1,1) .ne. BC_NEU) then
         u(-1,0:ny-1,1) = ZERO
      end if
      if (hi_inflow(1) .ne. 1 .or. face_type(1,2) .ne. BC_NEU) then
         u(nx,0:ny-1,1) = ZERO
      end if

      ! y-velocity
      u(-1,:,2) = ZERO
      u(nx,:,2) = ZERO
      if (lo_inflow(2) .ne. 1 .or. face_type(2,1) .ne. BC_NEU) then
         u(0:nx-1,-1,2) = ZERO
      end if
      if (hi_inflow(2) .ne. 1 .or. face_type(2,2) .ne. BC_NEU) then
         u(0:nx-1,ny,2) = ZERO
      end if
      
      do j = 0,ny
      do i = 0,nx
         rh(i,j) = HALF * (u(i  ,j,1) + u(i  ,j-1,1) &
                          -u(i-1,j,1) - u(i-1,j-1,1)) / dx(1) + &
                   HALF * (u(i,j  ,2) + u(i-1,j  ,2) &
                          -u(i,j-1,2) - u(i-1,j-1,2)) / dx(2)
      end do
      end do

      if (face_type(1,1) == BC_NEU) rh( 0,:) = TWO*rh( 0,:)
      if (face_type(1,2) == BC_NEU) rh(nx,:) = TWO*rh(nx,:)
      if (face_type(2,1) == BC_NEU) rh(:, 0) = TWO*rh(:, 0)
      if (face_type(2,2) == BC_NEU) rh(:,ny) = TWO*rh(:,ny)

    end subroutine grid_divu_2d

    subroutine grid_divu_3d(u,rh,dx,face_type,ng, lo_inflow, hi_inflow)

      integer        , intent(in   ) :: ng
      real(kind=dp_t), intent(inout) ::  u(-ng:,-ng:,-ng:,1:)
      real(kind=dp_t), intent(inout) :: rh(-1:,-1:,-1:)
      real(kind=dp_t), intent(in   ) :: dx(:)
      integer        , intent(in   ) :: face_type(:,:)
      integer        , intent(in   ) :: lo_inflow(:), hi_inflow(:)

      integer         :: i,j,k,nx,ny,nz
      real(kind=dp_t) :: ivdx, ivdy, ivdz

      nx = size(rh,dim=1) - 3
      ny = size(rh,dim=2) - 3
      nz = size(rh,dim=3) - 3

      ! x-velocity
      u(:,-1,:,1) = ZERO
      u(:,ny,:,1) = ZERO
      u(:,:,-1,1) = ZERO
      u(:,:,nz,1) = ZERO
      if (lo_inflow(1) .ne. 1 .or. face_type(1,1) .ne. BC_NEU) then
         u(-1,0:ny-1,0:nz-1,1) = ZERO
      end if
      if (hi_inflow(1) .ne. 1 .or. face_type(1,2) .ne. BC_NEU) then
         u(nx,0:ny-1,0:nz-1,1) = ZERO
      end if

      ! y-velocity
      u(-1,:,:,2) = ZERO
      u(nx,:,:,2) = ZERO
      u(:,:,-1,2) = ZERO
      u(:,:,nz,2) = ZERO
      if (lo_inflow(2) .ne. 1 .or. face_type(2,1) .ne. BC_NEU) then
         u(0:nx-1,-1,0:nz-1,2) = ZERO
      end if
      if (hi_inflow(2) .ne. 1 .or. face_type(2,2) .ne. BC_NEU) then
         u(0:nx-1,ny,0:nz-1,2) = ZERO
      end if

      ! z-velocity
      u(-1,:,:,3) = ZERO
      u(nx,:,:,3) = ZERO
      u(:,-1,:,3) = ZERO
      u(:,ny,:,3) = ZERO
      if (lo_inflow(3) .ne. 1 .or. face_type(3,1) .ne. BC_NEU) then
         u(0:nx-1,0:ny-1,-1,3) = ZERO
      end if
      if (hi_inflow(3) .ne. 1 .or. face_type(3,2) .ne. BC_NEU) then
         u(0:nx-1,0:ny-1,nz,3) = ZERO
      end if

      ivdx = ONE / dx(1)
      ivdy = ONE / dx(2)
      ivdz = ONE / dx(3)

      !$OMP PARALLEL DO PRIVATE(i,j,k)
      do k = 0,nz
      do j = 0,ny
      do i = 0,nx
         rh(i,j,k) = (u(i  ,j,k  ,1) + u(i  ,j-1,k  ,1) &
                     +u(i  ,j,k-1,1) + u(i  ,j-1,k-1,1) &
                     -u(i-1,j,k  ,1) - u(i-1,j-1,k  ,1) &
                     -u(i-1,j,k-1,1) - u(i-1,j-1,k-1,1)) * ivdx + &
                     (u(i,j  ,k  ,2) + u(i-1,j  ,k  ,2) &
                     +u(i,j  ,k-1,2) + u(i-1,j  ,k-1,2) &
                     -u(i,j-1,k  ,2) - u(i-1,j-1,k  ,2) &
                     -u(i,j-1,k-1,2) - u(i-1,j-1,k-1,2)) * ivdy + &
                     (u(i,j  ,k  ,3) + u(i-1,j  ,k  ,3) &
                     +u(i,j-1,k  ,3) + u(i-1,j-1,k  ,3) &
                     -u(i,j  ,k-1,3) - u(i-1,j  ,k-1,3) &
                     -u(i,j-1,k-1,3) - u(i-1,j-1,k-1,3)) * ivdz
         rh(i,j,k) = FOURTH*rh(i,j,k)
      end do
      end do
      end do
      !$OMP END PARALLEL DO

      if (face_type(1,1) == BC_NEU) rh( 0,:,:) = TWO*rh( 0,:,:)
      if (face_type(1,2) == BC_NEU) rh(nx,:,:) = TWO*rh(nx,:,:)
      if (face_type(2,1) == BC_NEU) rh(:, 0,:) = TWO*rh(:, 0,:)
      if (face_type(2,2) == BC_NEU) rh(:,ny,:) = TWO*rh(:,ny,:)
      if (face_type(3,1) == BC_NEU) rh(:,:, 0) = TWO*rh(:,:, 0)
      if (face_type(3,2) == BC_NEU) rh(:,:,nz) = TWO*rh(:,:,nz)

    end subroutine grid_divu_3d

    subroutine ml_crse_divu_contrib(rh, flux, u, mm, dx, crse_domain, ir, side, &
         lo_inflow, hi_inflow)
     type(multifab), intent(inout) :: rh
     type(multifab), intent(inout) :: flux
     type(multifab), intent(in   ) :: u
     type(imultifab),intent(in   ) :: mm
     real(kind=dp_t),intent(in   ) :: dx(:)
     type(box)      ,intent(in   ) :: crse_domain
     integer        ,intent(in   ) :: ir(:)
     integer        ,intent(in   ) :: side
     integer        ,intent(in   ) :: lo_inflow(:), hi_inflow(:)

     type(box) :: fbox, ubox, mbox, isect

     integer   :: lo (get_dim(rh)), hi (get_dim(rh)), lou(get_dim(rh)), dims(4), dm
     integer   :: lof(get_dim(rh)), hif(get_dim(rh)), lor(get_dim(rh)), lom(get_dim(rh))
     integer   :: lodom(get_dim(rh)), hidom(get_dim(rh)), dir, i, j, k, proc, ii, jj
     logical   :: nodal(get_dim(rh))
     logical   :: pmask(get_dim(rh))

     type(layout) :: flux_la

     integer,               parameter :: tag = 1371
     real(kind=dp_t),       pointer   :: rp(:,:,:,:), fp(:,:,:,:), up(:,:,:,:)
     integer,               pointer   :: mp(:,:,:,:)
     type(box_intersector), pointer   :: bi(:)
     type(bl_prof_timer),   save      :: bpt

     call build(bpt, "ml_crse_divu_contrib")

     if ( .not. cell_centered_q(flux) ) call bl_error('ml_crse_divu_contrib(): flux NOT cell centered')

     dims    = 1;
     nodal   = .true.
     dir     = iabs(side)
     lodom   = lwb(crse_domain)
     hidom   = upb(crse_domain)+1
     flux_la = get_layout(flux)
     dm      = get_dim(rh)
     pmask   = get_pmask(get_layout(rh))

     do j = 1, nboxes(u%la)

       ubox = box_nodalize(get_box(u%la,j),nodal)
       lou  = lwb(grow(box_nodalize(get_box( u%la,j), u%nodal), u%ng))
       lor  = lwb(grow(box_nodalize(get_box(rh%la,j),rh%nodal),rh%ng))

       bi => layout_get_box_intersector(flux_la, ubox)

       do k = 1, size(bi)

          i = bi(k)%i

          if ( remote(flux%la,i) .and. remote(u%la,j) ) cycle

          fbox  = box_nodalize(get_box(flux%la,i),flux%nodal)
          isect = bi(k)%bx
          lof   = lwb(fbox)
          hif   = upb(fbox)

          if ( (lof(dir) == lodom(dir) .or. lof(dir) == hidom(dir)) .and. (.not. pmask(dir)) ) cycle

          lo = lwb(isect)
          hi = upb(isect)

          if ( local(flux%la,i) .and. local(u%la,j) ) then
             ii   =  local_index(mm,i)
             jj   =  local_index(u, j)
             lom  =  lwb(get_pbox(mm,ii))
             fp   => dataptr(flux,ii)
             mp   => dataptr(mm  ,ii)
             up   => dataptr(u   ,jj)
             rp   => dataptr(rh  ,jj)
             select case (dm)
             case (1)
                call ml_interface_1d_divu(rp(:,1,1,1), lor, &
                     fp(:,1,1,1), lof, &
                     up(:,1,1,1), lou, mp(:,1,1,1), lom, lo, ir, side, dx)
             case (2)
                call ml_interface_2d_divu(rp(:,:,1,1), lor, &
                     fp(:,:,1,1), lof, lof, hif, &
                     up(:,:,1,:), lou, mp(:,:,1,1), lom, lo, hi, ir, side, dx, &
                     lo_inflow, hi_inflow)
             case (3)
                call ml_interface_3d_divu(rp(:,:,:,1), lor, &
                     fp(:,:,:,1), lof, lof, hif,  &
                     up(:,:,:,:), lou, mp(:,:,:,1), lom, lo, hi, ir, side, dx, &
                     lo_inflow, hi_inflow)
             end select

          else if ( local(flux%la,i) ) then
             !
             ! Must send flux & mm.
             !
             ii   = local_index(mm,i)
             mbox =  intersection(refine(isect,ir), get_pbox(mm,ii))
             fp   => dataptr(flux, ii, isect, 1, ncomp(flux))
             mp   => dataptr(mm,   ii, mbox,  1, ncomp(mm))
             proc =  get_proc(get_layout(u), j)
             call parallel_send(fp, proc, tag)
             call parallel_send(mp, proc, tag)

          else if ( local(u%la,j) ) then
             !
             ! Must receive flux & mm.
             !
             jj   = local_index(u,j)
             proc = get_proc(flux_la, i)
             mbox = intersection(refine(isect,ir), grow(box_nodalize(get_box(mm%la,i),mm%nodal),mm%ng))
             lom  = lwb(mbox)
             dims(1:dm) = extent(isect)
             allocate(fp(dims(1),dims(2),dims(3),ncomp(flux)))
             dims(1:dm) = extent(mbox)
             allocate(mp(dims(1),dims(2),dims(3),ncomp(mm)))
             call parallel_recv(fp, proc, tag)
             call parallel_recv(mp, proc, tag)
             up => dataptr(u, jj)
             rp => dataptr(rh,jj)
             select case (dm)
             case (1)
                call ml_interface_1d_divu(rp(:,1,1,1), lor, &
                     fp(:,1,1,1), lo, &
                     up(:,1,1,1), lou, mp(:,1,1,1), lom, lo, ir, side, dx)
             case (2)
                call ml_interface_2d_divu(rp(:,:,1,1), lor, &
                     fp(:,:,1,1), lo, lof, hif, &
                     up(:,:,1,:), lou, mp(:,:,1,1), lom, lo, hi, ir, side, dx, &
                     lo_inflow, hi_inflow)
             case (3)
                call ml_interface_3d_divu(rp(:,:,:,1), lor, &
                     fp(:,:,:,1), lo, lof, hif, &
                     up(:,:,:,:), lou, mp(:,:,:,1), lom, lo, hi, ir, side, dx, &
                     lo_inflow, hi_inflow)
             end select
             deallocate(fp,mp)
          end if
       end do
       deallocate(bi)
    end do
    call destroy(bpt)

   end subroutine ml_crse_divu_contrib

    subroutine ml_interface_1d_divu(rh, lor, fine_flux, lof, uc, loc, &
                                    mm, lom, lo, ir, side, dx)
    integer, intent(in) :: lor(:)
    integer, intent(in) :: loc(:)
    integer, intent(in) :: lom(:)
    integer, intent(in) :: lof(:)
    integer, intent(in) :: lo(:)
    real (kind = dp_t), intent(inout) ::        rh(lor(1):)
    real (kind = dp_t), intent(in   ) :: fine_flux(lof(1):)
    real (kind = dp_t), intent(in   ) ::        uc(loc(1):)
    integer           , intent(in   ) ::        mm(lom(1):)
    integer           , intent(in   ) :: ir(:)
    integer           , intent(in   ) :: side
    real(kind=dp_t)   , intent(in   ) :: dx(:)

    integer :: i
    real (kind = dp_t) :: crse_flux,fac

    i = lo(1)

!   NOTE: MM IS ON THE FINE GRID, NOT THE CRSE

    fac = ONE / ir(1)

!   Lo i side
    if (side == -1) then

          if (bc_dirichlet(mm(ir(1)*i),1,0)) then

             crse_flux = uc(i) / dx(1)

             rh(i) = rh(i) - crse_flux + fac*fine_flux(i)

          end if

!   Hi i side
    else if (side ==  1) then

          if (bc_dirichlet(mm(ir(1)*i),1,0)) then

             crse_flux = -uc(i-1) / dx(1)
             
             rh(i) = rh(i) - crse_flux + fac*fine_flux(i)

          end if

    end if

  end subroutine ml_interface_1d_divu

    subroutine ml_interface_2d_divu(rh, lor, fine_flux, lof, loflx, hiflx, uc, loc, &
                                    mm, lom, lo, hi, ir, side, dx, lo_inflow, hi_inflow)
    integer, intent(in) :: lor(:)
    integer, intent(in) :: loc(:)
    integer, intent(in) :: lom(:)
    integer, intent(in) :: lof(:)
    integer, intent(in) :: loflx(:), hiflx(:)
    integer, intent(in) :: lo(:), hi(:)
    real (kind = dp_t), intent(inout) ::        rh(lor(1):,lor(2):)
    real (kind = dp_t), intent(in   ) :: fine_flux(lof(1):,lof(2):)
    real (kind = dp_t), intent(in   ) ::        uc(loc(1):,loc(2):,:)
    integer           , intent(in   ) ::        mm(lom(1):,lom(2):)
    integer           , intent(in   ) :: ir(:)
    integer           , intent(in   ) :: side
    real(kind=dp_t)   , intent(in   ) :: dx(:)
    integer, intent(in) :: lo_inflow(:), hi_inflow(:)

    integer :: i, j
    real (kind = dp_t) :: crse_flux,fac,ivdx,ivdy

    i = lo(1)
    j = lo(2)

    ivdx = ONE / dx(1)
    ivdy = ONE / dx(2)

!   NOTE: THESE STENCILS ONLY WORK FOR DX == DY.

!   NOTE: MM IS ON THE FINE GRID, NOT THE CRSE

    fac = ONE / (ir(1)*ir(2))

!   Lo i side
    if (side == -1) then

       do j = lo(2),hi(2)

          if (bc_dirichlet(mm(ir(1)*i,ir(2)*j),1,0)) then
             if (j == loflx(2)) then
                if (bc_neumann(mm(ir(1)*i,ir(2)*j),2,-1)) then
                   if (lo_inflow(2) .eq. 1) then
                      crse_flux = uc(i,j,1)*ivdx + (uc(i,j,2) - uc(i,j-1,2))*ivdy
                   else
                      crse_flux =        (uc(i,j,1)*ivdx + uc(i,j,2)*ivdy)
                   end if
                else 
                   ! We have FOURTH rather than HALF here because
                   ! point (i,j) will be touched again when side == -2.
                   ! So in the end, the total crse_flux subtracted from rh(i,j)
                   ! will be HALF*(uc(i,j,1)*ivdx + uc(i,j,2)*ivdy)
                   crse_flux = FOURTH*(uc(i,j,1)*ivdx + uc(i,j,2)*ivdy)
                end if
             else if (j == hiflx(2)) then
                if (bc_neumann(mm(ir(1)*i,ir(2)*j),2,+1)) then
                   if (hi_inflow(2) .eq. 1) then
                      crse_flux = uc(i,j-1,1)*ivdx + (uc(i,j,2) - uc(i,j-1,2))*ivdy
                   else
                      crse_flux =        (uc(i,j-1,1)*ivdx - uc(i,j-1,2)*ivdy)
                   end if
                else 
                   crse_flux = FOURTH*(uc(i,j-1,1)*ivdx - uc(i,j-1,2)*ivdy)
                end if
             else
                crse_flux = (HALF*(uc(i,j,1) + uc(i,j-1,1))*ivdx &
                            +HALF*(uc(i,j,2) - uc(i,j-1,2))*ivdy)
             end if

             rh(i,j) = rh(i,j) - crse_flux + fac*fine_flux(i,j)

          end if

       end do

!   Hi i side
    else if (side ==  1) then

       do j = lo(2),hi(2)

          if (bc_dirichlet(mm(ir(1)*i,ir(2)*j),1,0)) then
             if (j == loflx(2)) then
                if (bc_neumann(mm(ir(1)*i,ir(2)*j),2,-1)) then
                   if (lo_inflow(2) .eq. 1) then
                      crse_flux = -uc(i-1,j,1)*ivdx + (uc(i-1,j,2)-uc(i-1,j-1,2))*ivdy
                   else
                      crse_flux =        (-uc(i-1,j,1)*ivdx + uc(i-1,j,2)*ivdy)
                   end if
                else
                   crse_flux = FOURTH*(-uc(i-1,j,1)*ivdx + uc(i-1,j,2)*ivdy)
                end if 
             else if (j == hiflx(2)) then
                if (bc_neumann(mm(ir(1)*i,ir(2)*j),2,+1)) then
                   if (hi_inflow(2) .eq. 1) then
                      crse_flux = -uc(i-1,j-1,1)*ivdx + (uc(i-1,j,2)-uc(i-1,j-1,2))*ivdy
                   else
                      crse_flux =        (-uc(i-1,j-1,1)*ivdx - uc(i-1,j-1,2)*ivdy)
                   end if
                else 
                   crse_flux = FOURTH*(-uc(i-1,j-1,1)*ivdx - uc(i-1,j-1,2)*ivdy) 
                end if
             else
                crse_flux = (HALF*(-uc(i-1,j,1)-uc(i-1,j-1,1))*ivdx  &
                            +HALF*( uc(i-1,j,2)-uc(i-1,j-1,2))*ivdy)
             end if
             
             rh(i,j) = rh(i,j) - crse_flux + fac*fine_flux(i,j)

          end if

       end do

    ! Lo j side
    else if (side == -2) then

       do i = lo(1),hi(1)

          if (bc_dirichlet(mm(ir(1)*i,ir(2)*j),1,0)) then
             if (i == loflx(1)) then
                if (bc_neumann(mm(ir(1)*i,ir(2)*j),1,-1)) then
                   if (lo_inflow(1) .eq. 1) then
                      crse_flux = (uc(i,j,1)-uc(i-1,j,1))*ivdx + uc(i,j,2)*ivdy
                   else
                      crse_flux =        (uc(i,j,1)*ivdx + uc(i,j,2)*ivdy)
                   end if
                else 
                   crse_flux = FOURTH*(uc(i,j,1)*ivdx + uc(i,j,2)*ivdy)
                end if
             else if (i == hiflx(1)) then
                if (bc_neumann(mm(ir(1)*i,ir(2)*j),1,+1)) then
                   if (hi_inflow(1) .eq. 1) then
                      crse_flux = (uc(i,j,1)-uc(i-1,j,1))*ivdx + uc(i-1,j,2)*ivdy
                   else
                      crse_flux =        (-uc(i-1,j,1)*ivdx + uc(i-1,j,2)*ivdy)
                   end if
                else 
                   crse_flux = FOURTH*(-uc(i-1,j,1)*ivdx + uc(i-1,j,2)*ivdy)
                end if
             else
                crse_flux = (HALF*(uc(i,j,1)-uc(i-1,j,1))*ivdx  &
                            +HALF*(uc(i,j,2)+uc(i-1,j,2))*ivdy)
             end if
             rh(i,j) = rh(i,j) - crse_flux + fac*fine_flux(i,j)

          end if

       end do

    ! Hi j side
    else if (side ==  2) then

       do i = lo(1),hi(1)

          if (bc_dirichlet(mm(ir(1)*i,ir(2)*j),1,0)) then

             if (i == loflx(1)) then
                if (bc_neumann(mm(ir(1)*i,ir(2)*j),1,-1)) then
                   if (lo_inflow(1) .eq. 1) then
                      crse_flux = (uc(i,j-1,1)-uc(i-1,j-1,1))*ivdx - uc(i,j-1,2)*ivdy
                   else
                      crse_flux =        (uc(i,j-1,1)*ivdx - uc(i,j-1,2)*ivdy)
                   end if
                else
                   crse_flux = FOURTH*(uc(i,j-1,1)*ivdx - uc(i,j-1,2)*ivdy)
                end if

             else if (i == hiflx(1)) then
                if (bc_neumann(mm(ir(1)*i,ir(2)*j),1,+1)) then
                   if (hi_inflow(1) .eq. 1) then
                      crse_flux = (uc(i,j-1,1)-uc(i-1,j-1,1))*ivdx - uc(i-1,j-1,2)*ivdy
                   else
                      crse_flux =        (-uc(i-1,j-1,1)*ivdx - uc(i-1,j-1,2)*ivdy)
                   end if
                else 
                   crse_flux = FOURTH*(-uc(i-1,j-1,1)*ivdx - uc(i-1,j-1,2)*ivdy)
                end if
             else
                crse_flux = (HALF*( uc(i,j-1,1)-uc(i-1,j-1,1))*ivdx &
                            +HALF*(-uc(i,j-1,2)-uc(i-1,j-1,2))*ivdy)
             end if

             rh(i,j) = rh(i,j) - crse_flux + fac*fine_flux(i,j)

          end if
       end do

    end if

  end subroutine ml_interface_2d_divu

    subroutine ml_interface_3d_divu(rh, lor, fine_flux, lof, loflx, hiflx, uc, loc, &
                                    mm, lom, lo, hi, ir, side, dx, lo_inflow, hi_inflow)
    integer, intent(in) :: lor(:)
    integer, intent(in) :: loc(:)
    integer, intent(in) :: lom(:)
    integer, intent(in) :: lof(:)
    integer, intent(in) :: loflx(:), hiflx(:)
    integer, intent(in) :: lo(:), hi(:)
    real (kind = dp_t), intent(inout) ::        rh(lor(1):,lor(2):,lor(3):)
    real (kind = dp_t), intent(in   ) :: fine_flux(lof(1):,lof(2):,lof(3):)
    real (kind = dp_t), intent(in   ) ::        uc(loc(1):,loc(2):,loc(3):,:)
    integer           , intent(in   ) ::        mm(lom(1):,lom(2):,lom(3):)
    integer           , intent(in   ) :: ir(:)
    integer           , intent(in   ) :: side
    real(kind=dp_t)   , intent(in   ) :: dx(:)
    integer, intent(in) :: lo_inflow(:), hi_inflow(:)

    integer :: i, j, k, ii, jj, kk, ifac
    logical :: lo_i_neu,lo_j_neu,lo_k_neu,hi_i_neu,hi_j_neu,hi_k_neu
    logical :: lo_i_not,lo_j_not,lo_k_not,hi_i_not,hi_j_not,hi_k_not
    real (kind = dp_t) :: cell_pp,cell_mp,cell_pm,cell_mm
    real (kind = dp_t) :: crse_flux,fac,ivdx,ivdy,ivdz

    ii = lo(1)
    jj = lo(2)
    kk = lo(3)

    ivdx = ONE / dx(1)
    ivdy = ONE / dx(2)
    ivdz = ONE / dx(3)

!   NOTE: THESE STENCILS ONLY WORK FOR DX == DY == DZ.

!   NOTE: MM IS ON THE FINE GRID, NOT THE CRSE

    fac = ONE / (ir(1)*ir(2)*ir(3))

!   Lo/Hi i side
    if (( side == -1) .or. (side == 1) ) then
 
      if (side == -1) then
        i    = ii
        ifac = 1
      else
        i    = ii-1
        ifac = -1
      end if

      !$OMP PARALLEL DO PRIVATE(j,k,lo_j_not,hi_j_not,lo_k_not,hi_k_not) &
      !$OMP PRIVATE(lo_j_neu,hi_j_neu,lo_k_neu,hi_k_neu) &
      !$OMP PRIVATE(cell_mm,cell_pm,cell_mp,cell_pp,crse_flux)
      do k = lo(3),hi(3)
      do j = lo(2),hi(2)

        if (bc_dirichlet(mm(ir(1)*ii,ir(2)*j,ir(3)*k),1,0)) then

          lo_j_not = .false.
          hi_j_not = .false.
          lo_j_neu = .false.
          hi_j_neu = .false.
          lo_k_not = .false.
          hi_k_not = .false.
          lo_k_neu = .false.
          hi_k_neu = .false.

          cell_pp =  (uc(i,j  ,k  ,1) )*ivdx * ifac &
                   + (uc(i,j  ,k  ,2) )*ivdy &
                   + (uc(i,j  ,k  ,3) )*ivdz
          cell_pm =  (uc(i,j  ,k-1,1) )*ivdx * ifac &
                   + (uc(i,j  ,k-1,2) )*ivdy &
                   - (uc(i,j  ,k-1,3) )*ivdz
          cell_mp =  (uc(i,j-1,k  ,1) )*ivdx * ifac &
                   - (uc(i,j-1,k  ,2) )*ivdy &
                   + (uc(i,j-1,k  ,3) )*ivdz
          cell_mm =  (uc(i,j-1,k-1,1) )*ivdx * ifac &
                   - (uc(i,j-1,k-1,2) )*ivdy &
                   - (uc(i,j-1,k-1,3) )*ivdz

          if (j == loflx(2)) then
             if (.not. bc_neumann(mm(ir(1)*ii,ir(2)*j,ir(3)*k),2,-1)) lo_j_not = .true.
             if (bc_neumann(mm(ir(1)*ii,ir(2)*j,ir(3)*k),2,-1))       lo_j_neu = .true.
          end if

          if (j == hiflx(2)) then
             if (.not. bc_neumann(mm(ir(1)*ii,ir(2)*j,ir(3)*k),2,+1)) hi_j_not = .true.
             if (bc_neumann(mm(ir(1)*ii,ir(2)*j,ir(3)*k),2,+1))       hi_j_neu = .true.
          end if

          if (k == loflx(3)) then
             if (.not. bc_neumann(mm(ir(1)*ii,ir(2)*j,ir(3)*k),3,-1)) lo_k_not = .true.
             if (bc_neumann(mm(ir(1)*ii,ir(2)*j,ir(3)*k),3,-1))       lo_k_neu = .true.
          end if

          if (k == hiflx(3)) then
             if (.not. bc_neumann(mm(ir(1)*ii,ir(2)*j,ir(3)*k),3,+1)) hi_k_not = .true.
             if (bc_neumann(mm(ir(1)*ii,ir(2)*j,ir(3)*k),3,+1))       hi_k_neu = .true.
          end if

          if (lo_k_not) then ! lo-k: coarse-fine
             if (lo_j_not) then  ! lo-j: coarse-fine
                crse_flux = THIRD*cell_pp
             else if (lo_j_neu) then 
                crse_flux = cell_pp
                if (lo_inflow(2) .eq. 1) then
                   crse_flux = crse_flux - uc(i,j-1,k,2)*ivdy
                end if
             else if (hi_j_not) then
                crse_flux = THIRD*cell_mp
             else if (hi_j_neu) then
                crse_flux = cell_mp
                if (hi_inflow(2) .eq. 1) then
                   crse_flux = crse_flux + uc(i,j,k,2)*ivdy
                end if
             else
                crse_flux = HALF*(cell_pp + cell_mp)
             end if
          else if (lo_k_neu) then
             if (lo_j_not) then
                crse_flux = cell_pp
                if (lo_inflow(3) .eq. 1) then
                   crse_flux = crse_flux - uc(i,j,k-1,3)*ivdz
                end if
             else if (lo_j_neu) then
                crse_flux = FOUR*cell_pp
                if (lo_inflow(2) .eq. 1) then
                   crse_flux = crse_flux - FOUR*uc(i,j-1,k,2)*ivdy
                end if
                if (lo_inflow(3) .eq. 1) then
                   crse_flux = crse_flux - FOUR*uc(i,j,k-1,3)*ivdz
                end if
             else if (hi_j_not) then
                crse_flux = cell_mp
                if (lo_inflow(3) .eq. 1) then
                   crse_flux = crse_flux - uc(i,j-1,k-1,3)*ivdz
                end if
             else if (hi_j_neu) then
                crse_flux = FOUR*cell_mp
                if (hi_inflow(2) .eq. 1) then
                   crse_flux = crse_flux + FOUR*uc(i,j,k,2)*ivdy
                end if
                if (lo_inflow(3) .eq. 1) then
                   crse_flux = crse_flux - FOUR*uc(i,j-1,k-1,3)*ivdz
                end if
             else
                crse_flux = TWO*(cell_pp + cell_mp)
                if (lo_inflow(3) .eq. 1) then
                   crse_flux = crse_flux - TWO*(uc(i,j,k-1,3) + uc(i,j-1,k-1,3))*ivdz
                end if
             end if
          else if (hi_k_not) then
             if (lo_j_not) then
                crse_flux = THIRD*cell_pm
             else if (lo_j_neu) then
                crse_flux = cell_pm
                if (lo_inflow(2) .eq. 1) then
                   crse_flux = crse_flux - uc(i,j-1,k-1,2)*ivdy
                end if
             else if (hi_j_not) then
                crse_flux = THIRD*cell_mm
             else if (hi_j_neu) then
                crse_flux = cell_mm
                if (hi_inflow(2) .eq. 1) then
                   crse_flux = crse_flux + uc(i,j,k-1,2)*ivdy
                end if
             else
                crse_flux = HALF*(cell_pm  + cell_mm)
             end if
          else if (hi_k_neu) then
             if (lo_j_not) then
                crse_flux = cell_pm
                if (hi_inflow(3) .eq. 1) then
                   crse_flux = crse_flux + uc(i,j,k,3)*ivdz
                end if
             else if (lo_j_neu) then
                crse_flux = FOUR*cell_pm
                if (lo_inflow(2) .eq. 1) then
                   crse_flux = crse_flux - FOUR*uc(i,j-1,k-1,2)*ivdy
                end if
                if (hi_inflow(3) .eq. 1) then
                   crse_flux = crse_flux + FOUR*uc(i,j,k,3)*ivdz
                end if
             else if (hi_j_not) then
                crse_flux = cell_mm
                if (hi_inflow(3) .eq. 1) then
                   crse_flux = crse_flux + uc(i,j-1,k,3)*ivdz
                end if
             else if (hi_j_neu) then
                crse_flux = FOUR*cell_mm
                if (hi_inflow(2) .eq. 1) then
                   crse_flux = crse_flux + FOUR*uc(i,j,k-1,2)*ivdy
                end if
                if (hi_inflow(3) .eq. 1) then
                   crse_flux = crse_flux + FOUR*uc(i,j-1,k,3)*ivdz
                end if
             else
                crse_flux = TWO*(cell_pm  + cell_mm)
                if (hi_inflow(3) .eq. 1) then
                   crse_flux = crse_flux + TWO*(uc(i,j,k,3)+uc(i,j-1,k,3))*ivdz
                end if
             end if
          else
             if (lo_j_not) then
                crse_flux = HALF*(cell_pm  + cell_pp)
             else if (lo_j_neu) then
                crse_flux = TWO*(cell_pm  + cell_pp)
                if (lo_inflow(2) .eq. 1) then
                   crse_flux = crse_flux - TWO*(uc(i,j-1,k-1,2)+uc(i,j-1,k,2))*ivdy
                end if
             else if (hi_j_not) then
                crse_flux = HALF*(cell_mm  + cell_mp)
             else if (hi_j_neu) then
                crse_flux = TWO*(cell_mm  + cell_mp)
                if (hi_inflow(2) .eq. 1) then
                   crse_flux = crse_flux + TWO*(uc(i,j,k-1,2)+uc(i,j,k,2))*ivdy
                end if
             else
                crse_flux = cell_mm  + cell_mp + cell_pm + cell_pp
             end if
          end if

          rh(ii,j,k) = rh(ii,j,k) - FOURTH*crse_flux + fac*fine_flux(ii,j,k)
        end if

      end do
      end do
      !$OMP END PARALLEL DO

!   Lo/Hi j side
    else if (( side == -2) .or. (side == 2) ) then
 
      if (side == -2) then
        j    = jj
        ifac = 1
      else
        j    = jj-1
        ifac = -1
      end if

      !$OMP PARALLEL DO PRIVATE(i,k,lo_i_not,hi_i_not,lo_k_not,hi_k_not) &
      !$OMP PRIVATE(lo_i_neu,hi_i_neu,lo_k_neu,hi_k_neu) &
      !$OMP PRIVATE(cell_mm,cell_pm,cell_mp,cell_pp,crse_flux)
      do k = lo(3),hi(3)
      do i = lo(1),hi(1)

        if (bc_dirichlet(mm(ir(1)*i,ir(2)*jj,ir(3)*k),1,0)) then

          lo_i_not = .false.
          hi_i_not = .false.
          lo_i_neu = .false.
          hi_i_neu = .false.
          lo_k_not = .false.
          hi_k_not = .false.
          lo_k_neu = .false.
          hi_k_neu = .false.

          if (i == loflx(1)) then
             if (.not. bc_neumann(mm(ir(1)*i,ir(2)*jj,ir(3)*k),1,-1)) lo_i_not = .true.
             if (bc_neumann(mm(ir(1)*i,ir(2)*jj,ir(3)*k),1,-1))       lo_i_neu = .true.
          end if

          if (i == hiflx(1)) then
             if (.not. bc_neumann(mm(ir(1)*i,ir(2)*jj,ir(3)*k),1,+1)) hi_i_not = .true.
             if (bc_neumann(mm(ir(1)*i,ir(2)*jj,ir(3)*k),1,+1))       hi_i_neu = .true.
          end if

          if (k == loflx(3)) then
             if (.not. bc_neumann(mm(ir(1)*i,ir(2)*jj,ir(3)*k),3,-1)) lo_k_not = .true.
             if (bc_neumann(mm(ir(1)*i,ir(2)*jj,ir(3)*k),3,-1))       lo_k_neu = .true.
          end if
          if (k == hiflx(3)) then
             if (.not. bc_neumann(mm(ir(1)*i,ir(2)*jj,ir(3)*k),3,+1)) hi_k_not = .true.
             if (bc_neumann(mm(ir(1)*i,ir(2)*jj,ir(3)*k),3,+1))       hi_k_neu = .true.
          end if

          cell_pp =  (uc(i  ,j,k  ,1) )*ivdx &
                    +(uc(i  ,j,k  ,2) )*ivdy * ifac &
                    +(uc(i  ,j,k  ,3) )*ivdz
          cell_pm =  (uc(i  ,j,k-1,1) )*ivdx &
                    +(uc(i  ,j,k-1,2) )*ivdy * ifac &
                    -(uc(i  ,j,k-1,3) )*ivdz
          cell_mp = -(uc(i-1,j,k  ,1) )*ivdx &
                    +(uc(i-1,j,k  ,2) )*ivdy * ifac &
                    +(uc(i-1,j,k  ,3) )*ivdz
          cell_mm = -(uc(i-1,j,k-1,1) )*ivdx &
                    +(uc(i-1,j,k-1,2) )*ivdy * ifac &
                    -(uc(i-1,j,k-1,3) )*ivdz

          if (lo_k_not) then
             if (lo_i_not) then
                crse_flux = THIRD*cell_pp
             else if (lo_i_neu) then
                crse_flux = cell_pp
                if (lo_inflow(1) .eq. 1) then
                   crse_flux = crse_flux - uc(i-1,j,k,1)*ivdx
                end if
             else if (hi_i_not) then
                crse_flux = THIRD*cell_mp
             else if (hi_i_neu) then
                crse_flux = cell_mp
                if (hi_inflow(1) .eq. 1) then
                   crse_flux = crse_flux + uc(i,j,k,1)*ivdx
                end if
             else
                crse_flux = HALF*(cell_pp + cell_mp)
             end if
          else if (lo_k_neu) then
             if (lo_i_not) then
                crse_flux = cell_pp
                if (lo_inflow(3) .eq. 1) then
                   crse_flux = crse_flux - uc(i,j,k-1,3)*ivdz
                end if
             else if (lo_i_neu) then
                crse_flux = FOUR*cell_pp
                if (lo_inflow(1) .eq. 1) then
                   crse_flux = crse_flux - FOUR*uc(i-1,j,k,1)*ivdx
                end if
                if (lo_inflow(3) .eq. 1) then
                   crse_flux = crse_flux - FOUR*uc(i,j,k-1,3)*ivdz
                end if
             else if (hi_i_not) then
                crse_flux = cell_mp
                if (lo_inflow(3) .eq. 1) then
                   crse_flux = crse_flux - uc(i-1,j,k-1,3)*ivdz
                end if
             else if (hi_i_neu) then
                crse_flux = FOUR*cell_mp
                if (hi_inflow(1) .eq. 1) then
                   crse_flux = crse_flux + FOUR*uc(i,j,k,1)*ivdx
                end if
                if (lo_inflow(3) .eq. 1) then
                   crse_flux = crse_flux - FOUR*uc(i-1,j,k-1,3)*ivdz
                end if
             else
                crse_flux = TWO*(cell_pp + cell_mp)
                if (lo_inflow(3) .eq. 1) then
                   crse_flux = crse_flux - TWO*(uc(i,j,k-1,3)+uc(i-1,j,k-1,3))*ivdz
                end if
             end if
          else if (hi_k_not) then
             if (lo_i_not) then
                crse_flux = THIRD*cell_pm
             else if (lo_i_neu) then
                crse_flux = cell_pm
                if (lo_inflow(1) .eq. 1) then
                   crse_flux = crse_flux - uc(i-1,j,k-1,1)*ivdx
                end if
             else if (hi_i_not) then
                crse_flux = THIRD*cell_mm
             else if (hi_i_neu) then
                crse_flux = cell_mm
                if (hi_inflow(1) .eq. 1) then
                   crse_flux = crse_flux + uc(i,j,k-1,1)*ivdx
                end if
             else
                crse_flux = HALF*(cell_pm  + cell_mm)
             end if
          else if (hi_k_neu) then
             if (lo_i_not) then
                crse_flux = cell_pm
                if (hi_inflow(3) .eq. 1) then
                   crse_flux = crse_flux + uc(i,j,k,3)*ivdz
                end if
             else if (lo_i_neu) then
                crse_flux = FOUR*cell_pm
                if (lo_inflow(1) .eq. 1) then
                   crse_flux = crse_flux - FOUR*uc(i-1,j,k-1,1)*ivdx
                end if
                if (hi_inflow(3) .eq. 1) then
                   crse_flux = crse_flux + FOUR*uc(i,j,k,3)*ivdz
                end if
             else if (hi_i_not) then
                crse_flux = cell_mm
                if (hi_inflow(3) .eq. 1) then
                   crse_flux = crse_flux + uc(i-1,j,k,3)*ivdz
                end if
             else if (hi_i_neu) then
                crse_flux = FOUR*cell_mm
                if (hi_inflow(1) .eq. 1) then
                   crse_flux = crse_flux + FOUR*uc(i,j,k-1,1)*ivdx
                end if
                if (hi_inflow(3) .eq. 1) then
                   crse_flux = crse_flux + FOUR*uc(i-1,j,k,3)*ivdz
                end if
             else
                crse_flux = TWO*(cell_pm  + cell_mm)
                if (hi_inflow(3) .eq. 1) then
                   crse_flux = crse_flux + TWO*(uc(i,j,k,3)+uc(i-1,j,k,3))*ivdz
                end if
             end if
          else
             if (lo_i_not) then
                crse_flux = HALF*(cell_pm  + cell_pp)
             else if (lo_i_neu) then
                crse_flux = TWO*(cell_pm  + cell_pp)
                if (lo_inflow(1) .eq. 1) then
                   crse_flux = crse_flux -TWO*(uc(i-1,j,k-1,1)+uc(i-1,j,k,1))*ivdx
                end if
             else if (hi_i_not) then
                crse_flux = HALF*(cell_mm  + cell_mp)
             else if (hi_i_neu) then
                crse_flux = TWO*(cell_mm  + cell_mp)
                if (hi_inflow(1) .eq. 1) then
                   crse_flux = crse_flux + TWO*(uc(i,j,k-1,1)+uc(i,j,k,1))*ivdx
                end if
             else
                crse_flux = cell_mm  + cell_mp + cell_pm + cell_pp
             end if
          end if

          rh(i,jj,k) = rh(i,jj,k) - FOURTH*crse_flux + fac*fine_flux(i,jj,k)
        end if

      end do
      end do
      !$OMP END PARALLEL DO

!   Lo/Hi k side
    else if (( side == -3) .or. (side == 3) ) then
 
      if (side == -3) then
        k    = kk
        ifac = 1
      else
        k    = kk-1
        ifac = -1
      end if

      !$OMP PARALLEL DO PRIVATE(i,j,lo_i_not,hi_i_not,lo_j_not,hi_j_not) &
      !$OMP PRIVATE(lo_i_neu,hi_i_neu,lo_j_neu,hi_j_neu) &
      !$OMP PRIVATE(cell_mm,cell_pm,cell_mp,cell_pp,crse_flux)
      do j = lo(2),hi(2)
      do i = lo(1),hi(1)

        if (bc_dirichlet(mm(ir(1)*i,ir(2)*j,ir(3)*kk),1,0)) then

          lo_i_not = .false.
          hi_i_not = .false.
          lo_i_neu = .false.
          hi_i_neu = .false.
          lo_j_not = .false.
          hi_j_not = .false.
          lo_j_neu = .false.
          hi_j_neu = .false.

          if (i == loflx(1)) then
             if (.not. bc_neumann(mm(ir(1)*i,ir(2)*j,ir(3)*kk),1,-1)) lo_i_not = .true.
             if (bc_neumann(mm(ir(1)*i,ir(2)*j,ir(3)*kk),1,-1))       lo_i_neu = .true.
          end if

          if (i == hiflx(1)) then
             if (.not. bc_neumann(mm(ir(1)*i,ir(2)*j,ir(3)*kk),1,+1)) hi_i_not = .true.
             if (bc_neumann(mm(ir(1)*i,ir(2)*j,ir(3)*kk),1,+1))       hi_i_neu = .true.
          end if

          if (j == loflx(2)) then
             if (.not. bc_neumann(mm(ir(1)*i,ir(2)*j,ir(3)*kk),2,-1)) lo_j_not = .true.
             if (bc_neumann(mm(ir(1)*i,ir(2)*j,ir(3)*kk),2,-1))       lo_j_neu = .true.
          end if

          if (j == hiflx(2)) then
             if (.not. bc_neumann(mm(ir(1)*i,ir(2)*j,ir(3)*kk),2,+1)) hi_j_not = .true.
             if (bc_neumann(mm(ir(1)*i,ir(2)*j,ir(3)*kk),2,+1))       hi_j_neu = .true.
          end if

          cell_pp =  (uc(i  ,j  ,k,1) )*ivdx &
                    +(uc(i  ,j  ,k,2) )*ivdy &
                    +(uc(i  ,j  ,k,3) )*ivdz * ifac
          cell_pm =  (uc(i  ,j-1,k,1) )*ivdx &
                    -(uc(i  ,j-1,k,2) )*ivdy &
                    +(uc(i  ,j-1,k,3) )*ivdz * ifac
          cell_mp = -(uc(i-1,j  ,k,1) )*ivdx &
                    +(uc(i-1,j  ,k,2) )*ivdy &
                    +(uc(i-1,j  ,k,3) )*ivdz * ifac
          cell_mm = -(uc(i-1,j-1,k,1) )*ivdx &
                    -(uc(i-1,j-1,k,2) )*ivdy &
                    +(uc(i-1,j-1,k,3) )*ivdz * ifac

          if (lo_j_not) then
             if (lo_i_not) then
                crse_flux = THIRD*cell_pp
             else if (lo_i_neu) then
                crse_flux = cell_pp
                if (lo_inflow(1) .eq. 1) then
                   crse_flux = crse_flux - uc(i-1,j,k,1)*ivdx
                end if
             else if (hi_i_not) then
                crse_flux = THIRD*cell_mp
             else if (hi_i_neu) then
                crse_flux = cell_mp
                if (hi_inflow(1) .eq. 1) then
                   crse_flux = crse_flux + uc(i,j,k,1)*ivdx
                end if
             else
                crse_flux = HALF*(cell_pp + cell_mp)
             end if
          else if (lo_j_neu) then
             if (lo_i_not) then
                crse_flux = cell_pp
                if (lo_inflow(2) .eq. 1) then
                   crse_flux = crse_flux - uc(i,j-1,k,2)*ivdy
                end if
             else if (lo_i_neu) then
                crse_flux = FOUR*cell_pp
                if (lo_inflow(1) .eq. 1) then
                   crse_flux = crse_flux - FOUR*uc(i-1,j,k,1)*ivdx
                end if
                if (lo_inflow(2) .eq. 1) then
                   crse_flux = crse_flux - FOUR*uc(i,j-1,k,2)*ivdy
                end if
             else if (hi_i_not) then
                crse_flux = cell_mp
                if (lo_inflow(2) .eq. 1) then
                   crse_flux = crse_flux - uc(i-1,j-1,k,2)*ivdy
                end if
             else if (hi_i_neu) then
                crse_flux = FOUR*cell_mp
                if (hi_inflow(1) .eq. 1) then
                   crse_flux = crse_flux + FOUR*uc(i,j,k,1)*ivdx
                end if
                if (lo_inflow(2) .eq. 1) then
                   crse_flux = crse_flux - FOUR*uc(i-1,j-1,k,2)*ivdy
                end if
             else
                crse_flux = TWO*(cell_pp + cell_mp)
                if (lo_inflow(2) .eq. 1) then
                   crse_flux = crse_flux - TWO*(uc(i,j-1,k,2)+uc(i-1,j-1,k,2))*ivdy
                end if
             end if
          else if (hi_j_not) then
             if (lo_i_not) then
                crse_flux = THIRD*cell_pm
             else if (lo_i_neu) then
                crse_flux = cell_pm
                if (lo_inflow(1) .eq. 1) then
                   crse_flux = crse_flux - uc(i-1,j-1,k,1)*ivdx
                end if
             else if (hi_i_not) then
                crse_flux = THIRD*cell_mm
             else if (hi_i_neu) then
                crse_flux = cell_mm
                if (hi_inflow(1) .eq. 1) then
                   crse_flux = crse_flux + uc(i,j-1,k,1)*ivdx
                end if
             else
                crse_flux = HALF*(cell_pm  + cell_mm)
             end if
          else if (hi_j_neu) then
             if (lo_i_not) then
                crse_flux = cell_pm
                if (hi_inflow(2) .eq. 1) then
                   crse_flux = crse_flux + uc(i,j,k,2)*ivdy
                end if
             else if (lo_i_neu) then
                crse_flux = FOUR*cell_pm
                if (lo_inflow(1) .eq. 1) then
                   crse_flux = crse_flux - FOUR*uc(i-1,j-1,k,1)*ivdx
                end if
                if (hi_inflow(2) .eq. 1) then
                   crse_flux = crse_flux + FOUR*uc(i,j,k,2)*ivdy
                end if
             else if (hi_i_not) then
                crse_flux = cell_mm
                if (hi_inflow(2) .eq. 1) then
                   crse_flux = crse_flux + uc(i-1,j,k,2)*ivdy
                end if
             else if (hi_i_neu) then
                crse_flux = FOUR*cell_mm
                if (hi_inflow(1) .eq. 1) then
                   crse_flux = crse_flux + FOUR*uc(i,j-1,k,1)*ivdx
                end if
                if (hi_inflow(2) .eq. 1) then
                   crse_flux = crse_flux + FOUR*uc(i-1,j,k,2)*ivdy
                end if
             else
                crse_flux = TWO*(cell_pm  + cell_mm)
                if (hi_inflow(2) .eq. 1) then
                   crse_flux = crse_flux + TWO*(uc(i,j,k,2)+uc(i-1,j,k,2))*ivdy
                end if
             end if
          else
             if (lo_i_not) then
                crse_flux = HALF*(cell_pm  + cell_pp)
             else if (lo_i_neu) then
                crse_flux = TWO*(cell_pm  + cell_pp)
                if (lo_inflow(1) .eq. 1) then
                   crse_flux = crse_flux - TWO*(uc(i-1,j-1,k,1)+uc(i-1,j,k,1))*ivdx
                end if
             else if (hi_i_not) then
                crse_flux = HALF*(cell_mm  + cell_mp)
             else if (hi_i_neu) then
                crse_flux = TWO*(cell_mm  + cell_mp)
                if (hi_inflow(1) .eq. 1) then
                   crse_flux = crse_flux + TWO*(uc(i,j-1,k,1)+uc(i,j,k,1))*ivdx
                end if
             else
                crse_flux = cell_mm  + cell_mp + cell_pm + cell_pp
             end if
          end if
  
          rh(i,j,kk) = rh(i,j,kk) - FOURTH*crse_flux + fac*fine_flux(i,j,kk)
        end if

      end do
      end do
      !$OMP END PARALLEL DO

    end if

  end subroutine ml_interface_3d_divu


    subroutine subtract_divu_from_rh(nlevs,mgt,rh,divu_rhs,subtract_flag_in)

      integer        , intent(in   ) :: nlevs
      type(mg_tower) , intent(inout) :: mgt(:)
      type(multifab) , intent(inout) :: rh(:)
      type(multifab) , intent(in   ) :: divu_rhs(:)
      logical, intent(in), optional  :: subtract_flag_in

      real(kind=dp_t), pointer :: dp(:,:,:,:) 
      real(kind=dp_t), pointer :: rp(:,:,:,:) 
      integer        , pointer :: mp(:,:,:,:) 

      integer :: i,n,dm,ng_r,ng_d
      integer :: mglev_fine
      logical :: subtract_flag

      subtract_flag = .true.
      if (present(subtract_flag_in)) subtract_flag = subtract_flag_in

      dm   = get_dim(rh(nlevs))
      ng_r = nghost(rh(nlevs))
      ng_d = nghost(divu_rhs(nlevs))

!     Create the regular single-level divergence.
      do n = 1, nlevs
         mglev_fine = mgt(n)%nlevels
         do i = 1, nfabs(rh(n))
            rp => dataptr(rh(n)      , i)
            dp => dataptr(divu_rhs(n), i)
            mp => dataptr(mgt(n)%mm(mglev_fine),i)
            select case (dm)
               case (1)
                 call subtract_divu_from_rh_1d(rp(:,1,1,1), ng_r, &
                                               dp(:,1,1,1), ng_d, &
                                               mp(:,1,1,1), subtract_flag )
               case (2)
                 call subtract_divu_from_rh_2d(rp(:,:,1,1), ng_r, &
                                               dp(:,:,1,1), ng_d, &
                                               mp(:,:,1,1), subtract_flag )
               case (3)
                 call subtract_divu_from_rh_3d(rp(:,:,:,1), ng_r, &
                                               dp(:,:,:,1), ng_d, &
                                               mp(:,:,:,1), subtract_flag )
            end select
         end do
      end do

    end subroutine subtract_divu_from_rh

    subroutine subtract_divu_from_rh_1d(rh,ng_rh,divu_rhs,ng_divu,mm,subtract_flag)

      integer        , intent(in   ) :: ng_rh,ng_divu
      real(kind=dp_t), intent(inout) ::       rh(  -ng_rh:)
      real(kind=dp_t), intent(inout) :: divu_rhs(-ng_divu:)
      integer        , intent(inout) :: mm(0:)
      logical        , intent(in   ) :: subtract_flag

      integer         :: i,nx

      nx = size(mm,dim=1) - 1

      if (subtract_flag) then
         do i = 0,nx
            if (.not. bc_dirichlet(mm(i),1,0)) &
                 rh(i) = rh(i) - divu_rhs(i)
         end do
      else
         do i = 0,nx
            if (.not. bc_dirichlet(mm(i),1,0)) &
                 rh(i) = rh(i) + divu_rhs(i)
         end do
      end if

    end subroutine subtract_divu_from_rh_1d

    subroutine subtract_divu_from_rh_2d(rh,ng_rh,divu_rhs,ng_divu,mm,subtract_flag)

      integer        , intent(in   ) :: ng_rh,ng_divu
      real(kind=dp_t), intent(inout) ::       rh(  -ng_rh:,  -ng_rh:)
      real(kind=dp_t), intent(inout) :: divu_rhs(-ng_divu:,-ng_divu:)
      integer        , intent(inout) :: mm(0:,0:)
      logical        , intent(in   ) :: subtract_flag

      integer         :: i,j,nx,ny

      nx = size(mm,dim=1) - 1
      ny = size(mm,dim=2) - 1

      if (subtract_flag) then
         do j = 0,ny
         do i = 0,nx
            if (.not. bc_dirichlet(mm(i,j),1,0)) &
                 rh(i,j) = rh(i,j) - divu_rhs(i,j)
         end do
         end do
      else
         do j = 0,ny
         do i = 0,nx
            if (.not. bc_dirichlet(mm(i,j),1,0)) &
                 rh(i,j) = rh(i,j) + divu_rhs(i,j)
         end do
         end do
      end if

    end subroutine subtract_divu_from_rh_2d

    subroutine subtract_divu_from_rh_3d(rh,ng_rh,divu_rhs,ng_divu,mm,subtract_flag)

      integer        , intent(in   ) :: ng_rh,ng_divu
      real(kind=dp_t), intent(inout) ::       rh(  -ng_rh:,  -ng_rh:,   -ng_rh:)
      real(kind=dp_t), intent(inout) :: divu_rhs(-ng_divu:,-ng_divu:, -ng_divu:)
      integer        , intent(inout) :: mm(0:,0:,0:)
      logical        , intent(in   ) :: subtract_flag

      integer         :: i,j,k,nx,ny,nz

      nx = size(mm,dim=1) - 1
      ny = size(mm,dim=2) - 1
      nz = size(mm,dim=3) - 1

      if (subtract_flag) then
         !$OMP PARALLEL DO PRIVATE(i,j,k)
         do k = 0,nz
         do j = 0,ny
         do i = 0,nx
            if (.not. bc_dirichlet(mm(i,j,k),1,0)) &
                 rh(i,j,k) = rh(i,j,k) - divu_rhs(i,j,k)
         end do
         end do
         end do
         !$OMP END PARALLEL DO
      else
         !$OMP PARALLEL DO PRIVATE(i,j,k)
         do k = 0,nz
         do j = 0,ny
         do i = 0,nx
            if (.not. bc_dirichlet(mm(i,j,k),1,0)) &
                 rh(i,j,k) = rh(i,j,k) + divu_rhs(i,j,k)
         end do
         end do
         end do
         !$OMP END PARALLEL DO
      end if

    end subroutine subtract_divu_from_rh_3d

  ! have_divu stuff

  subroutine divucc(nlevs,mgt,rhcc,rh,ref_ratio,nodal)

    integer        , intent(in   ) :: nlevs
    type(mg_tower) , intent(inout) :: mgt(:)
    type(multifab) , intent(inout) :: rhcc(:)
    type(multifab) , intent(inout) :: rh(:)
    integer        , intent(in   ) :: ref_ratio(:,:)
    logical        , intent(in   ) :: nodal(:)
    
    real(kind=dp_t), pointer :: rcp(:,:,:,:) 
    real(kind=dp_t), pointer :: rhp(:,:,:,:) 
    integer        , pointer ::  mp(:,:,:,:) 
    
    integer         :: i,n,dm,ng,mglev_fine
    type(      box) :: pdc
    type(   layout) :: la_crse,la_fine
    type(bndry_reg) :: brs_flx

    dm = get_dim(rhcc(nlevs))
    ng = nghost(rhcc(nlevs))
    
    ! regular single-level divergence
    do n = 1, nlevs
       mglev_fine = mgt(n)%nlevels
       call multifab_fill_boundary(rhcc(n))
       do i = 1, nfabs(rhcc(n))
          rcp => dataptr(rhcc(n), i)
          rhp => dataptr(rh(n)  , i)
          mp  => dataptr(mgt(n)%mm(mglev_fine),i)
          select case (dm)
          case (1)
             call bl_error('divucc_1d not implemented')
          case (2)
             call divucc_2d(rcp(:,:,1,1), rhp(:,:,1,1), &
                  &          mp(:,:,1,1), mgt(n)%face_type(i,:,:), ng)
          case (3)
             call divucc_3d(rcp(:,:,:,1), rhp(:,:,:,1), &
                  &          mp(:,:,:,1), mgt(n)%face_type(i,:,:), ng)
          end select
       end do
    end do

    ! Modify rh at coarse-fine interfaces.
    do n = nlevs,2,-1

       la_fine = get_layout(rhcc(n))
       la_crse = get_layout(rhcc(n-1))
       pdc = get_pd(la_crse)
       
       call bndry_reg_rr_build_nd(brs_flx,la_fine,ref_ratio(n-1,:), pdc, nodal=nodal)
       call crse_fine_divucc(n,rh(n-1),rhcc,brs_flx,ref_ratio(n-1,:),mgt)
       call bndry_reg_destroy(brs_flx)
    end do
    
  end subroutine divucc

  
  subroutine divucc_2d(rc,rh,mm,face_type,ng)

    integer        , intent(in   ) :: ng
    real(kind=dp_t), intent(inout) :: rc(-ng:,-ng:)
    real(kind=dp_t), intent(inout) :: rh(-1:,-1:)
    integer        , intent(inout) :: mm(0:,0:)
    integer        , intent(in   ) :: face_type(:,:)
    
    integer         :: i,j,nx,ny
    real(kind=dp_t), pointer   :: rhtmp(:,:)
    
    nx = size(rh,dim=1) - 3
    ny = size(rh,dim=2) - 3
    
    allocate(rhtmp(0:nx,0:ny))

    do j = 0,ny
    do i = 0,nx
       if (.not.bc_dirichlet(mm(i,j),1,0)) then 
          rhtmp(i,j) = (rc(i-1,j-1)+rc(i,j-1)+rc(i-1,j)+rc(i,j))*FOURTH
       else
          rhtmp(i,j) = ZERO 
       end if
    end do
    end do

    if (face_type(1,1) == BC_NEU) rhtmp( 0,:) = TWO*rhtmp( 0,:)
    if (face_type(1,2) == BC_NEU) rhtmp(nx,:) = TWO*rhtmp(nx,:)
    if (face_type(2,1) == BC_NEU) rhtmp(:, 0) = TWO*rhtmp(:, 0)
    if (face_type(2,2) == BC_NEU) rhtmp(:,ny) = TWO*rhtmp(:,ny)

    do j = 0,ny
    do i = 0,nx
       rh(i,j) = rh(i,j) + rhtmp(i,j)
    end do
    end do

    deallocate(rhtmp)
    
  end subroutine divucc_2d
  
  subroutine divucc_3d(rc,rh,mm,face_type,ng)

    integer        , intent(in   ) :: ng
    real(kind=dp_t), intent(inout) :: rc(-ng:,-ng:,-ng:)
    real(kind=dp_t), intent(inout) :: rh(-1:,-1:,-1:)
    integer        , intent(inout) :: mm(0:,0:,0:)
    integer        , intent(in   ) :: face_type(:,:)
    
    integer         :: i,j,k,nx,ny,nz
    real(kind=dp_t), pointer   :: rhtmp(:,:,:)
    
    nx = size(rh,dim=1) - 3
    ny = size(rh,dim=2) - 3
    nz = size(rh,dim=3) - 3
    
    allocate(rhtmp(0:nx,0:ny,0:nz))

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = 0,nz
    do j = 0,ny
    do i = 0,nx
       if (.not.bc_dirichlet(mm(i,j,k),1,0)) then 
          rhtmp(i,j,k) = EIGHTH *  &
               (rc(i-1,j-1,k-1)+rc(i,j-1,k-1)+rc(i-1,j,k-1)+rc(i,j,k-1) &
               +rc(i-1,j-1,k  )+rc(i,j-1,k  )+rc(i-1,j,k  )+rc(i,j,k  ) )
       else
          rhtmp(i,j,k) = ZERO 
       end if
    end do
    end do
    end do
    !$OMP END PARALLEL DO

    if (face_type(1,1) == BC_NEU) rhtmp( 0,:,:) = TWO*rhtmp( 0,:,:)
    if (face_type(1,2) == BC_NEU) rhtmp(nx,:,:) = TWO*rhtmp(nx,:,:)
    if (face_type(2,1) == BC_NEU) rhtmp(:, 0,:) = TWO*rhtmp(:, 0,:)
    if (face_type(2,2) == BC_NEU) rhtmp(:,ny,:) = TWO*rhtmp(:,ny,:)
    if (face_type(3,1) == BC_NEU) rhtmp(:,:, 0) = TWO*rhtmp(:,:, 0)
    if (face_type(3,2) == BC_NEU) rhtmp(:,:,nz) = TWO*rhtmp(:,:,nz)

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = 0,nz
    do j = 0,ny
    do i = 0,nx
       rh(i,j,k) = rh(i,j,k) + rhtmp(i,j,k)
    end do
    end do
    end do
    !$OMP END PARALLEL DO

    deallocate(rhtmp)
    
  end subroutine divucc_3d
  
  subroutine crse_fine_divucc(n_fine,rh_crse,rhcc,brs_flx,ref_ratio,mgt)

    integer        , intent(in   ) :: n_fine
    type(multifab) , intent(inout) :: rh_crse
    type(multifab) , intent(in   ) :: rhcc(:)
    type(bndry_reg), intent(inout) :: brs_flx
    integer        , intent(in   ) :: ref_ratio(:)
    type(mg_tower) , intent(in   ) :: mgt(:)

    type(multifab) :: temp_rhs_crse
    type(  layout) :: la_crse,la_fine
    type(     box) :: pdc
    integer :: i,dm,n_crse, ng
    integer :: mglev_fine, mglev_crse
    logical :: nodal(get_dim(rhcc(n_fine)))

    dm = get_dim(rhcc(n_fine))
    ng = nghost(rhcc(n_fine))
    n_crse = n_fine-1
    
    nodal = .true.
    
    la_crse = get_layout(rhcc(n_crse))
    la_fine = get_layout(rhcc(n_fine))
    
    mglev_crse = mgt(n_crse)%nlevels
    mglev_fine = mgt(n_fine)%nlevels
  
    call multifab_build(temp_rhs_crse, la_crse, 1, 1, nodal)
    call setval(temp_rhs_crse, ZERO, 1, all=.true.)

    ! Zero out the flux registers which will hold the fine contributions
    call bndry_reg_setval(brs_flx, ZERO, all = .true.)

    ! compute the fine contributions

    pdc = get_pd(la_crse)

    do i = 1, dm
       call ml_fine_rhcc_contrib(brs_flx%bmf(i,0), &
            rhcc(n_fine),mgt(n_fine)%mm(mglev_fine),ref_ratio,pdc,-i)
       call ml_fine_rhcc_contrib(brs_flx%bmf(i,1), &
            rhcc(n_fine),mgt(n_fine)%mm(mglev_fine),ref_ratio,pdc,+i)
    end do

    ! compute the crse contributions
    do i = 1, dm
       call ml_crse_rhcc_contrib(temp_rhs_crse, brs_flx%bmf(i,0), &
            rhcc(n_crse), mgt(n_fine)%mm(mglev_fine), pdc,ref_ratio, -i) 
       call ml_crse_rhcc_contrib(temp_rhs_crse, brs_flx%bmf(i,1), &
            rhcc(n_crse), mgt(n_fine)%mm(mglev_fine), pdc,ref_ratio, +i) 
    end do

    call multifab_plus_plus(rh_crse,temp_rhs_crse)
    call periodic_add_copy(rh_crse,temp_rhs_crse,synced=.true.)

    call multifab_destroy(temp_rhs_crse)

  end subroutine crse_fine_divucc

  subroutine ml_fine_rhcc_contrib(flux, rhcc, mm, ratio, crse_domain, side)

    use bl_prof_module
    type(multifab), intent(inout) :: flux
    type(multifab), intent(in   ) :: rhcc
    type(imultifab), intent(in) :: mm
    type(box) :: crse_domain
    type(box) :: fbox
    integer :: side
    integer :: ratio(:)
    integer :: lof(get_dim(flux)), dm
    integer :: lo_dom(get_dim(flux)), hi_dom(get_dim(flux))
    integer :: i, dir
    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: rp(:,:,:,:)
    integer        , pointer :: mp(:,:,:,:)
    logical :: pmask(get_dim(rhcc))
    type(bl_prof_timer), save :: bpt

    call build(bpt, "ml_fine_rhcc_contrib")

    lo_dom = lwb(crse_domain)
    hi_dom = upb(crse_domain) + 1
    dir    = iabs(side)
    pmask  = get_pmask(get_layout(rhcc))
    dm     = get_dim(flux)

    do i = 1, nfabs(flux)
       fbox =  get_ibox(flux,i)
       lof  =  lwb(fbox)
       fp   => dataptr(flux, i)
       rp   => dataptr(rhcc, i)
       mp   => dataptr(mm, i)
       if ( pmask(dir) .or. &
            (lof(dir) /= lo_dom(dir) .and. lof(dir) /= hi_dom(dir)) ) then
          select case(dm)
          case (1)
             call bl_error("ml_fine_rhcc_contrib_1d not implemented");
          case (2)
             call ml_fine_rhcc_contrib_2d(fp(:,:,1,1), rp(:,:,1,1), mp(:,:,1,1), ratio, side)
          case (3)
             call ml_fine_rhcc_contrib_3d(fp(:,:,:,1), rp(:,:,:,1), mp(:,:,:,1), ratio, side)
          end select
       end if
    enddo

    call destroy(bpt)

  end subroutine ml_fine_rhcc_contrib

  subroutine ml_fine_rhcc_contrib_2d(flx, rhcc, mm, ratio, side)

    real (kind = dp_t), intent(inout) :: flx( 0:, 0:)
    real (kind = dp_t), intent(in   ) :: rhcc(-1:,-1:)
    integer           , intent(in   ) ::   mm( 0:, 0:)
    integer, intent(in) :: ratio(:), side
    integer :: nxcc, nycc, nxf, nyf
    integer :: icc, jcc, iif, jjf
    integer :: istart, iend, jstart, jend
    real (kind = dp_t) :: xc, yc, rrfac, foo, fx, fy, freflect

    nxf = size(flx,dim=1)
    nyf = size(flx,dim=2)

    nxcc = size(rhcc,dim=1)-2
    nycc = size(rhcc,dim=2)-2

    rrfac = ONE / (ratio(1) * ratio(2))

    if (side == -1 .or. side == 1) then ! Lo/Hi i side

       iif = 0

       if (side == -1) then
          istart = 0
          iend = ratio(1)-1
          xc = -HALF
       else
          istart = nxcc-ratio(1)
          iend = nxcc-1
          xc = iend + HALF
       endif

       do jjf = 0, nyf-1
          if (jjf == 0) then
             jstart = 0
             jend = ratio(2)-1
             yc = -HALF
             if (bc_neumann(mm(istart,jstart),2,-1)) then
                freflect = TWO
             else
                freflect = ONE
             end if
          else if (jjf == nyf-1) then
             jstart = nycc - ratio(2)
             jend = nycc-1
             yc = jend + HALF
             if (bc_neumann(mm(istart,jend+1),2,+1)) then
                freflect = TWO
             else
                freflect = ONE
             end if
          else
             jstart = jjf * ratio(2) - ratio(2)
             jend = jjf * ratio(2) + ratio(2) - 1
             yc = jjf * ratio(2) - HALF
             freflect = ONE
          end if

          foo = ZERO
          do jcc = jstart, jend
             fy = ONE - abs((jcc-yc)/ratio(2))
             do icc = istart, iend
                fx = ONE - abs((icc-xc)/ratio(1))
                foo = foo + rhcc(icc,jcc) * fx * fy
             end do
          enddo
          flx(iif,jjf) = flx(iif,jjf) + foo * rrfac * freflect

       enddo

    else if (side == -2 .or. side == 2) then ! Lo/Hi j side

       jjf = 0

       if (side == -2) then
          jstart = 0
          jend = ratio(2)-1
          yc = -HALF
       else
          jstart = nycc-ratio(2)
          jend = nycc-1
          yc = jend + HALF
       endif

       do iif = 0, nxf-1

          if (iif == 0) then
             istart = 0
             iend = ratio(1)-1
             xc = -HALF
             if (bc_neumann(mm(istart,jstart),1,-1)) then
                freflect = TWO
             else
                freflect = ONE
             end if
          else if (iif == nxf-1) then
             istart = nxcc - ratio(1)
             iend = nxcc-1
             xc = iend + HALF
             if (bc_neumann(mm(iend+1,jstart),1,+1)) then
                freflect = TWO
             else
                freflect = ONE
             end if
          else
             istart = iif * ratio(1) - ratio(1)
             iend = iif * ratio(1) + ratio(1) - 1
             xc = iif * ratio(1) - HALF
             freflect = ONE
          end if

          foo = ZERO
          do jcc = jstart, jend
             fy = ONE - abs((jcc-yc)/ratio(2))
             do icc = istart, iend
                fx = ONE - abs((icc-xc)/ratio(1))
                foo = foo + rhcc(icc,jcc) * fx * fy
             end do
          enddo
          flx(iif,jjf) = flx(iif,jjf) + foo * rrfac * freflect

       enddo

    end if

  end subroutine ml_fine_rhcc_contrib_2d


  subroutine ml_fine_rhcc_contrib_3d(flx, rhcc, mm, ratio, side)

    real (kind = dp_t), intent(inout) :: flx( 0:, 0:, 0:)
    real (kind = dp_t), intent(in   ) :: rhcc(-1:,-1:,-1:)
    integer           , intent(in   ) ::   mm( 0:, 0:, 0:)
    integer, intent(in) :: ratio(:), side
    integer :: nxcc, nycc, nzcc, nxf, nyf, nzf
    integer :: icc, jcc, kcc, iif, jjf, kkf
    integer :: istart, iend, jstart, jend, kstart, kend
    real (kind = dp_t) :: xc, yc, zc, rrfac, foo, fx, fy, fz, freflect

    nxf = size(flx,dim=1)
    nyf = size(flx,dim=2)
    nzf = size(flx,dim=3)

    nxcc = size(rhcc,dim=1)-2
    nycc = size(rhcc,dim=2)-2
    nzcc = size(rhcc,dim=3)-2

    rrfac = ONE / (ratio(1) * ratio(2) * ratio(3))

    if (side == -1 .or. side == 1) then 

       iif = 0

       if (side == -1) then
          istart = 0
          iend = ratio(1)-1
          xc = -HALF
       else
          istart = nxcc-ratio(1)
          iend = nxcc-1
          xc = iend + HALF
       endif

       do kkf = 0, nzf-1

          if (kkf == 0) then
             kstart = 0
             kend = ratio(3)-1
             zc = -HALF
          else if (kkf == nzf-1) then
             kstart = nzcc - ratio(3)
             kend = nzcc-1
             zc = kend + HALF
          else
             kstart = kkf * ratio(3) - ratio(3)
             kend = kkf * ratio(3) + ratio(3) - 1
             zc = kkf * ratio(3) - HALF
          end if

          do jjf = 0, nyf-1

             if (jjf == 0) then
                jstart = 0
                jend = ratio(2)-1
                yc = -HALF
             else if (jjf == nyf-1) then
                jstart = nycc - ratio(2)
                jend = nycc-1
                yc = jend + HALF
             else
                jstart = jjf * ratio(2) - ratio(2)
                jend = jjf * ratio(2) + ratio(2) - 1
                yc = jjf * ratio(2) - HALF
             end if

             freflect = ONE

             if (bc_neumann(mm(istart,jstart,kstart),2,-1)) then
                freflect = freflect * TWO
             else if (bc_neumann(mm(istart,jend+1,kstart),2,+1)) then
                freflect = freflect * TWO
             end if

             if (bc_neumann(mm(istart,jstart,kstart),3,-1)) then
                freflect = freflect * TWO
             else if (bc_neumann(mm(istart,jstart,kend+1),3,+1)) then
                freflect = freflect * TWO
             end if

             foo = ZERO
             do kcc = kstart, kend
                fz = ONE - abs((kcc-zc)/ratio(3))
                do jcc = jstart, jend
                   fy = ONE - abs((jcc-yc)/ratio(2))
                   do icc = istart, iend
                      fx = ONE - abs((icc-xc)/ratio(1))
                      foo = foo + rhcc(icc,jcc,kcc) * fx * fy * fz
                   end do
                end do
             enddo

             flx(iif,jjf,kkf) = flx(iif,jjf,kkf) + foo * rrfac * freflect
             
          enddo
       enddo

    else if (side == -2 .or. side == 2) then 

       jjf = 0

       if (side == -2) then
          jstart = 0
          jend = ratio(2)-1
          yc = -HALF
       else
          jstart = nycc-ratio(2)
          jend = nycc-1
          yc = jend + HALF
       endif

       do kkf = 0, nzf-1

          if (kkf == 0) then
             kstart = 0
             kend = ratio(3)-1
             zc = -HALF
          else if (kkf == nzf-1) then
             kstart = nzcc - ratio(3)
             kend = nzcc-1
             zc = kend + HALF
          else
             kstart = kkf * ratio(3) - ratio(3)
             kend = kkf * ratio(3) + ratio(3) - 1
             zc = kkf * ratio(3) - HALF
          end if

          do iif = 0, nxf-1

             if (iif == 0) then
                istart = 0
                iend = ratio(1)-1
                xc = -HALF
             else if (iif == nxf-1) then
                istart = nxcc - ratio(1)
                iend = nxcc-1
                xc = iend + HALF
             else
                istart = iif * ratio(1) - ratio(1)
                iend = iif * ratio(1) + ratio(1) - 1
                xc = iif * ratio(1) - HALF
             end if

             freflect = ONE

             if (bc_neumann(mm(istart,jstart,kstart),1,-1)) then
                freflect = freflect * TWO
             else if (bc_neumann(mm(iend+1,jstart,kstart),1,+1)) then
                freflect = freflect * TWO
             end if

             if (bc_neumann(mm(istart,jstart,kstart),3,-1)) then
                freflect = freflect * TWO
             else if (bc_neumann(mm(istart,jstart,kend+1),3,+1)) then
                freflect = freflect * TWO
             end if

             foo = ZERO
             do kcc = kstart, kend
                fz = ONE - abs((kcc-zc)/ratio(3))
                do jcc = jstart, jend
                   fy = ONE - abs((jcc-yc)/ratio(2))
                   do icc = istart, iend
                      fx = ONE - abs((icc-xc)/ratio(1))
                      foo = foo + rhcc(icc,jcc,kcc) * fx * fy * fz
                   end do
                end do
             enddo

             flx(iif,jjf,kkf) = flx(iif,jjf,kkf) + foo * rrfac * freflect
             
          enddo
       enddo

    else if (side == -3 .or. side == 3) then 

       kkf = 0

       if (side == -3) then
          kstart = 0
          kend = ratio(3)-1
          zc = -HALF
       else
          kstart = nzcc-ratio(3)
          kend = nzcc-1
          zc = kend + HALF
       endif

       do jjf = 0, nyf-1

          if (jjf == 0) then
             jstart = 0
             jend = ratio(2)-1
             yc = -HALF
          else if (jjf == nyf-1) then
             jstart = nycc - ratio(2)
             jend = nycc-1
             yc = jend + HALF
          else
             jstart = jjf * ratio(2) - ratio(2)
             jend = jjf * ratio(2) + ratio(2) - 1
             yc = jjf * ratio(2) - HALF
          end if

          do iif = 0, nxf-1

             if (iif == 0) then
                istart = 0
                iend = ratio(1)-1
                xc = -HALF
             else if (iif == nxf-1) then
                istart = nxcc - ratio(1)
                iend = nxcc-1
                xc = iend + HALF
             else
                istart = iif * ratio(1) - ratio(1)
                iend = iif * ratio(1) + ratio(1) - 1
                xc = iif * ratio(1) - HALF
             end if

             freflect = ONE

             if (bc_neumann(mm(istart,jstart,kstart),1,-1)) then
                freflect = freflect * TWO
             else if (bc_neumann(mm(iend+1,jstart,kstart),1,+1)) then
                freflect = freflect * TWO
             end if

             if (bc_neumann(mm(istart,jstart,kstart),2,-1)) then
                freflect = freflect * TWO
             else if (bc_neumann(mm(istart,jend+1,kstart),2,+1)) then
                freflect = freflect * TWO
             end if

             foo = ZERO
             do kcc = kstart, kend
                fz = ONE - abs((kcc-zc)/ratio(3))
                do jcc = jstart, jend
                   fy = ONE - abs((jcc-yc)/ratio(2))
                   do icc = istart, iend
                      fx = ONE - abs((icc-xc)/ratio(1))
                      foo = foo + rhcc(icc,jcc,kcc) * fx * fy * fz
                   end do
                end do
             enddo

             flx(iif,jjf,kkf) = flx(iif,jjf,kkf) + foo * rrfac * freflect
             
          enddo
       enddo
       
    end if

  end subroutine ml_fine_rhcc_contrib_3d

  subroutine ml_crse_rhcc_contrib(rh, flux, rhcc, mm, crse_domain, ir, side)
    type(multifab), intent(inout) :: rh
    type(multifab), intent(inout) :: flux
    type(multifab), intent(in   ) :: rhcc
    type(imultifab),intent(in   ) :: mm
    type(box)      ,intent(in   ) :: crse_domain
    integer        ,intent(in   ) :: ir(:)
    integer        ,intent(in   ) :: side
    
    type(box) :: fbox, rcbox, mbox, isect

    integer   :: lo (get_dim(rh)), hi (get_dim(rh)), lorc(get_dim(rh)), dims(4), dm
    integer   :: lof(get_dim(rh)), hif(get_dim(rh)), lorh(get_dim(rh)), lom(get_dim(rh))
    integer   :: lodom(get_dim(rh)), hidom(get_dim(rh)), dir, i, j, k, proc, ii, jj
    logical   :: nodal(get_dim(rh))
    logical   :: pmask(get_dim(rh))

    type(layout) :: flux_la
    
    integer,               parameter :: tag = 1731
    real(kind=dp_t),       pointer   :: rhp(:,:,:,:), fp(:,:,:,:), rcp(:,:,:,:)
    integer,               pointer   :: mp(:,:,:,:)
    type(box_intersector), pointer   :: bi(:)
    type(bl_prof_timer),   save      :: bpt
    
    call build(bpt, "ml_crse_rhcc_contrib")

    if ( .not. cell_centered_q(flux) ) call bl_error('ml_crse_rhcc_contrib(): flux NOT cell centered')

    dims    = 1;
    nodal   = .true.
    dir     = iabs(side)
    lodom   = lwb(crse_domain)
    hidom   = upb(crse_domain)+1
    flux_la = get_layout(flux)
    dm      = get_dim(rh)
    pmask   = get_pmask(get_layout(rh))

    do j = 1, nboxes(rhcc%la)

       rcbox = box_nodalize(get_box(rhcc%la,j),nodal)
       lorc  = lwb(grow(rcbox,rhcc%ng))
       lorh  = lwb(grow(box_nodalize(get_box(rh%la,j),rh%nodal),rh%ng))
       
       bi => layout_get_box_intersector(flux_la, rcbox)

       do k = 1, size(bi)
          
          i = bi(k)%i

          if ( remote(flux%la,i) .and. remote(rhcc%la,j) ) cycle
          
          fbox  = get_box(flux%la,i)
          isect = bi(k)%bx
          lof   = lwb(fbox)
          hif   = upb(fbox)

          if ( (lof(dir) == lodom(dir) .or. lof(dir) == hidom(dir)) .and. (.not. pmask(dir)) ) cycle

          lo = lwb(isect)
          hi = upb(isect)

          if ( local(flux%la,i) .and. local(rhcc%la,j) ) then
             ii   =  local_index(mm,i)
             jj   =  local_index(rh,j)
             lom  =  lwb(get_pbox(mm,ii))
             fp   => dataptr(flux,ii)
             mp   => dataptr(mm  ,ii)
             rcp  => dataptr(rhcc,jj)
             rhp  => dataptr(rh  ,jj)
             select case (dm)
             case (1)
                call bl_error("ml_interface_rhcc_1d not done")
             case (2)
                call ml_interface_rhcc_2d(rhp(:,:,1,1), lorh, &
                     fp(:,:,1,1), lof, lof, hif, &
                     rcp(:,:,1,1), lorh, mp(:,:,1,1), lom, lo, hi, ir, side)
             case (3)
                call ml_interface_rhcc_3d(rhp(:,:,:,1), lorh, &
                     fp(:,:,:,1), lof, lof, hif, &
                     rcp(:,:,:,1), lorh, mp(:,:,:,1), lom, lo, hi, ir, side)
             end select
          else if ( local(flux%la,i) ) then
             !
             ! Must send flux & mm.
             !
             ii   =  local_index(mm,i)
             mbox =  intersection(refine(isect,ir), get_pbox(mm,ii))
             fp   => dataptr(flux, ii, isect, 1, ncomp(flux))
             mp   => dataptr(mm,   ii, mbox,  1, ncomp(mm))
             proc =  get_proc(get_layout(rhcc), j)
             call parallel_send(fp, proc, tag)
             call parallel_send(mp, proc, tag)
          else
             !
             ! Must receive flux & mm.
             !
             jj   = local_index(rh,j)
             proc = get_proc(flux_la, i)
             mbox = intersection(refine(isect,ir), grow(box_nodalize(get_box(mm%la,i),mm%nodal),mm%ng))
             lom  = lwb(mbox)
             dims(1:dm) = extent(isect)
             allocate(fp(dims(1),dims(2),dims(3),ncomp(flux)))
             dims(1:dm) = extent(mbox)
             allocate(mp(dims(1),dims(2),dims(3),ncomp(mm)))
             call parallel_recv(fp, proc, tag)
             call parallel_recv(mp, proc, tag)

             rcp => dataptr(rhcc,jj)
             rhp => dataptr(rh  ,jj)
             select case (dm)
             case (1)
                call bl_error("ml_interface_rhcc_1d not done")
             case (2)
                call ml_interface_rhcc_2d(rhp(:,:,1,1), lorh, &
                     fp(:,:,1,1), lo, lof, hif, &
                     rcp(:,:,1,1), lorh, mp(:,:,1,1), lom, lo, hi, ir, side)
             case (3)
                call ml_interface_rhcc_3d(rhp(:,:,:,1), lorh, &
                     fp(:,:,:,1), lo, lof, hif, &
                     rcp(:,:,:,1), lorh, mp(:,:,:,1), lom, lo, hi, ir, side)
             end select
             deallocate(fp,mp)
          end if
       end do
       deallocate(bi)
    end do

    call destroy(bpt)
  end subroutine ml_crse_rhcc_contrib

  subroutine ml_interface_rhcc_2d(rh, lor, fine_flux, lof, loflx, hiflx, rc, loc, &
       &                          mm, lom, lo, hi, ir, side)
    integer, intent(in) :: lor(:)
    integer, intent(in) :: loc(:)
    integer, intent(in) :: lom(:)
    integer, intent(in) :: lof(:)
    integer, intent(in) :: loflx(:), hiflx(:)
    integer, intent(in) :: lo(:), hi(:)
    real (kind = dp_t), intent(inout) ::        rh(lor(1):,lor(2):)
    real (kind = dp_t), intent(in   ) :: fine_flux(lof(1):,lof(2):)
    real (kind = dp_t), intent(in   ) ::        rc(loc(1):,loc(2):)
    integer           , intent(in   ) ::        mm(lom(1):,lom(2):)
    integer           , intent(in   ) :: ir(:)
    integer           , intent(in   ) :: side

    integer :: i, j
    real (kind = dp_t) :: crse_flux, fac

    i = lo(1)
    j = lo(2)

!   NOTE: MM IS ON THE FINE GRID, NOT THE CRSE

!   Lo i side
    if (side == -1) then

       do j = lo(2), hi(2)

          if (bc_dirichlet(mm(ir(1)*i,ir(2)*j),1,0)) then

             fac = ONE

             if (j == loflx(2)) then
                if (bc_neumann(mm(ir(1)*i,ir(2)*j),2,-1)) then
                   crse_flux = rc(i,j)*HALF                   
                else
                   crse_flux = rc(i,j)*FOURTH
                   fac = HALF  ! because the corner will touch again when side==-2
                end if
             else if (j == hiflx(2)) then
                if (bc_neumann(mm(ir(1)*i,ir(2)*j),2,+1)) then
                   crse_flux = rc(i,j-1)*HALF
                else
                   crse_flux = rc(i,j-1)*FOURTH
                   fac = HALF
                end if
             else
                crse_flux = (rc(i,j-1)+rc(i,j))*FOURTH
             end if

             rh(i,j) = rh(i,j) + fac * (fine_flux(i,j) - crse_flux)

          end if

       end do

!   Hi i side
    else if (side ==  1) then

       do j = lo(2),hi(2)

          if (bc_dirichlet(mm(ir(1)*i,ir(2)*j),1,0)) then

             fac = ONE

             if (j == loflx(2)) then
                if (bc_neumann(mm(ir(1)*i,ir(2)*j),2,-1)) then
                   crse_flux = rc(i-1,j)*HALF
                else
                   crse_flux = rc(i-1,j)*FOURTH
                   fac = HALF
                end if
             else if (j == hiflx(2)) then
                if (bc_neumann(mm(ir(1)*i,ir(2)*j),2,+1)) then
                   crse_flux = rc(i-1,j-1)*HALF
                else 
                   crse_flux = rc(i-1,j-1)*FOURTH
                   fac = HALF
                end if
             else
                crse_flux = (rc(i-1,j-1)+rc(i-1,j))*FOURTH
             end if

             rh(i,j) = rh(i,j) + fac * (fine_flux(i,j) - crse_flux)

          end if

       end do

! Lo j side
    else if (side == -2) then

       do i = lo(1),hi(1)
          
          if (bc_dirichlet(mm(ir(1)*i,ir(2)*j),1,0)) then

             fac = ONE

             if (i == loflx(1)) then
                if (bc_neumann(mm(ir(1)*i,ir(2)*j),1,-1)) then
                   crse_flux = rc(i,j)*HALF
                else 
                   crse_flux = rc(i,j)*FOURTH
                   fac = HALF
                end if
             else if (i == hiflx(1)) then
                if (bc_neumann(mm(ir(1)*i,ir(2)*j),1,+1)) then
                   crse_flux = rc(i-1,j)*HALF
                else 
                   crse_flux = rc(i-1,j)*FOURTH
                   fac = HALF
                end if
             else
                crse_flux = (rc(i-1,j)+rc(i,j))*FOURTH
             end if

             rh(i,j) = rh(i,j) + fac * (fine_flux(i,j) - crse_flux)
             
          end if
       
       end do

! Hi j side
    else if (side ==  2) then

       do i = lo(1),hi(1)

          if (bc_dirichlet(mm(ir(1)*i,ir(2)*j),1,0)) then

             fac = ONE

             if (i == loflx(1)) then
                if (bc_neumann(mm(ir(1)*i,ir(2)*j),1,-1)) then
                   crse_flux = rc(i,j-1)*HALF
                else
                   crse_flux = rc(i,j-1)*FOURTH
                   fac = HALF
                end if

             else if (i == hiflx(1)) then
                if (bc_neumann(mm(ir(1)*i,ir(2)*j),1,+1)) then
                   crse_flux = rc(i,j-1)*HALF
                else 
                   crse_flux = rc(i-1,j-1)*FOURTH
                   fac = HALF
                end if
             else
                crse_flux = (rc(i-1,j-1)+rc(i,j-1))*FOURTH
             end if

             rh(i,j) = rh(i,j) + fac * (fine_flux(i,j) - crse_flux)

          end if

       end do

    end if

  end subroutine ml_interface_rhcc_2d

  subroutine ml_interface_rhcc_3d(rh, lor, fine_flux, lof, loflx, hiflx, rc, loc, &
       &                          mm, lom, lo, hi, ir, side)
    integer, intent(in) :: lor(:)
    integer, intent(in) :: loc(:)
    integer, intent(in) :: lom(:)
    integer, intent(in) :: lof(:)
    integer, intent(in) :: loflx(:), hiflx(:)
    integer, intent(in) :: lo(:), hi(:)
    real (kind = dp_t), intent(inout) ::        rh(lor(1):,lor(2):,lor(3):)
    real (kind = dp_t), intent(in   ) :: fine_flux(lof(1):,lof(2):,lof(3):)
    real (kind = dp_t), intent(in   ) ::        rc(loc(1):,loc(2):,loc(3):)
    integer           , intent(in   ) ::        mm(lom(1):,lom(2):,lom(3):)
    integer           , intent(in   ) :: ir(:)
    integer           , intent(in   ) :: side

    integer :: i, j, k, ii, jj, kk
    logical :: lo_i_neu,lo_j_neu,lo_k_neu,hi_i_neu,hi_j_neu,hi_k_neu
    logical :: lo_i_not,lo_j_not,lo_k_not,hi_i_not,hi_j_not,hi_k_not
    real (kind = dp_t) :: cell_pp,cell_mp,cell_pm,cell_mm
    real (kind = dp_t) :: crse_flux,fac

    ii = lo(1)
    jj = lo(2)
    kk = lo(3)

!   NOTE: MM IS ON THE FINE GRID, NOT THE CRSE

!   Lo/Hi i side
    if (( side == -1) .or. (side == 1) ) then
 
      if (side == -1) then
        i    = ii
      else
        i    = ii-1
      end if

      !$OMP PARALLEL DO PRIVATE(j,k,lo_j_not,hi_j_not,lo_k_not,hi_k_not) &
      !$OMP PRIVATE(lo_j_neu,hi_j_neu,lo_k_neu,hi_k_neu) &
      !$OMP PRIVATE(cell_mm,cell_pm,cell_mp,cell_pp,crse_flux, fac)
      do k = lo(3),hi(3)
      do j = lo(2),hi(2)

        if (bc_dirichlet(mm(ir(1)*ii,ir(2)*j,ir(3)*k),1,0)) then

          cell_pp = rc(i,j  ,k  ) * EIGHTH
          cell_pm = rc(i,j  ,k-1) * EIGHTH
          cell_mp = rc(i,j-1,k  ) * EIGHTH
          cell_mm = rc(i,j-1,k-1) * EIGHTH

          lo_j_not = .false.
          hi_j_not = .false.
          lo_j_neu = .false.
          hi_j_neu = .false.
          lo_k_not = .false.
          hi_k_not = .false.
          lo_k_neu = .false.
          hi_k_neu = .false.

          if (j == loflx(2)) then
             if (.not. bc_neumann(mm(ir(1)*ii,ir(2)*j,ir(3)*k),2,-1)) lo_j_not = .true.
             if (bc_neumann(mm(ir(1)*ii,ir(2)*j,ir(3)*k),2,-1))       lo_j_neu = .true.
          end if

          if (j == hiflx(2)) then
             if (.not. bc_neumann(mm(ir(1)*ii,ir(2)*j,ir(3)*k),2,+1)) hi_j_not = .true.
             if (bc_neumann(mm(ir(1)*ii,ir(2)*j,ir(3)*k),2,+1))       hi_j_neu = .true.
          end if

          if (k == loflx(3)) then
             if (.not. bc_neumann(mm(ir(1)*ii,ir(2)*j,ir(3)*k),3,-1)) lo_k_not = .true.
             if (bc_neumann(mm(ir(1)*ii,ir(2)*j,ir(3)*k),3,-1))       lo_k_neu = .true.
          end if

          if (k == hiflx(3)) then
             if (.not. bc_neumann(mm(ir(1)*ii,ir(2)*j,ir(3)*k),3,+1)) hi_k_not = .true.
             if (bc_neumann(mm(ir(1)*ii,ir(2)*j,ir(3)*k),3,+1))       hi_k_neu = .true.
          end if

          if (lo_k_not) then
             if (lo_j_not) then
                crse_flux = THIRD*cell_pp
                fac = THIRD
             else if (lo_j_neu) then
                crse_flux = cell_pp
                fac = HALF
             else if (hi_j_not) then
                crse_flux = THIRD*cell_mp
                fac = THIRD
             else if (hi_j_neu) then
                crse_flux = cell_mp
                fac = HALF
             else
                crse_flux = HALF*(cell_pp + cell_mp)
                fac = HALF
             end if
          else if (lo_k_neu) then
             if (lo_j_not) then
                crse_flux = cell_pp
                fac = HALF
             else if (lo_j_neu) then
                crse_flux = FOUR*cell_pp
                fac = ONE
             else if (hi_j_not) then
                crse_flux = cell_mp
                fac = HALF
             else if (hi_j_neu) then
                crse_flux = FOUR*cell_mp
                fac = ONE
             else
                crse_flux = TWO*(cell_pp + cell_mp)
                fac = one
             end if
          else if (hi_k_not) then
             if (lo_j_not) then
                crse_flux = THIRD*cell_pm
                fac = THIRD
             else if (lo_j_neu) then
                crse_flux = cell_pm
                fac = HALF
             else if (hi_j_not) then
                crse_flux = THIRD*cell_mm
                fac = THIRD
             else if (hi_j_neu) then
                crse_flux = cell_mm
                fac = HALF
             else
                crse_flux = HALF*(cell_pm  + cell_mm)
                fac = HALF
             end if
          else if (hi_k_neu) then
             if (lo_j_not) then
                crse_flux = cell_pm
                fac = HALF
             else if (lo_j_neu) then
                crse_flux = FOUR*cell_pm
                fac = ONE
             else if (hi_j_not) then
                crse_flux = cell_mm
                fac = HALF
             else if (hi_j_neu) then
                crse_flux = FOUR*cell_mm
                fac = ONE
             else
                crse_flux = TWO*(cell_pm  + cell_mm)
                fac = ONE
             end if
          else
             if (lo_j_not) then
                crse_flux = HALF*(cell_pm  + cell_pp)
                fac = HALF
             else if (lo_j_neu) then
                crse_flux = TWO*(cell_pm  + cell_pp)
                fac = ONE
             else if (hi_j_not) then
                crse_flux = HALF*(cell_mm  + cell_mp)
                fac = HALF
             else if (hi_j_neu) then
                crse_flux = TWO*(cell_mm  + cell_mp)
                fac = ONE
             else
                crse_flux = cell_mm  + cell_mp + cell_pm + cell_pp
                fac = ONE
             end if
          end if

          rh(ii,j,k) = rh(ii,j,k) - crse_flux + fac*fine_flux(ii,j,k)
        end if

      end do
      end do
      !$OMP END PARALLEL DO

!   Lo/Hi j side
    else if (( side == -2) .or. (side == 2) ) then
 
      if (side == -2) then
        j    = jj
      else
        j    = jj-1
      end if

      !$OMP PARALLEL DO PRIVATE(i,k,lo_i_not,hi_i_not,lo_k_not,hi_k_not) &
      !$OMP PRIVATE(lo_i_neu,hi_i_neu,lo_k_neu,hi_k_neu) &
      !$OMP PRIVATE(cell_mm,cell_pm,cell_mp,cell_pp,crse_flux, fac)
      do k = lo(3),hi(3)
      do i = lo(1),hi(1)

        if (bc_dirichlet(mm(ir(1)*i,ir(2)*jj,ir(3)*k),1,0)) then

          lo_i_not = .false.
          hi_i_not = .false.
          lo_i_neu = .false.
          hi_i_neu = .false.
          lo_k_not = .false.
          hi_k_not = .false.
          lo_k_neu = .false.
          hi_k_neu = .false.

          if (i == loflx(1)) then
             if (.not. bc_neumann(mm(ir(1)*i,ir(2)*jj,ir(3)*k),1,-1)) lo_i_not = .true.
             if (bc_neumann(mm(ir(1)*i,ir(2)*jj,ir(3)*k),1,-1))       lo_i_neu = .true.
          end if

          if (i == hiflx(1)) then
             if (.not. bc_neumann(mm(ir(1)*i,ir(2)*jj,ir(3)*k),1,+1)) hi_i_not = .true.
             if (bc_neumann(mm(ir(1)*i,ir(2)*jj,ir(3)*k),1,+1))       hi_i_neu = .true.
          end if

          if (k == loflx(3)) then
             if (.not. bc_neumann(mm(ir(1)*i,ir(2)*jj,ir(3)*k),3,-1)) lo_k_not = .true.
             if (bc_neumann(mm(ir(1)*i,ir(2)*jj,ir(3)*k),3,-1))       lo_k_neu = .true.
          end if
          if (k == hiflx(3)) then
             if (.not. bc_neumann(mm(ir(1)*i,ir(2)*jj,ir(3)*k),3,+1)) hi_k_not = .true.
             if (bc_neumann(mm(ir(1)*i,ir(2)*jj,ir(3)*k),3,+1))       hi_k_neu = .true.
          end if

          cell_pp = rc(i  ,j,k  ) * EIGHTH
          cell_pm = rc(i  ,j,k-1) * EIGHTH
          cell_mp = rc(i-1,j,k  ) * EIGHTH
          cell_mm = rc(i-1,j,k-1) * EIGHTH

          if (lo_k_not) then
             if (lo_i_not) then
                crse_flux = THIRD*cell_pp
                fac = THIRD
             else if (lo_i_neu) then
                crse_flux = cell_pp
                fac = HALF
             else if (hi_i_not) then
                crse_flux = THIRD*cell_mp
                fac = THIRD
             else if (hi_i_neu) then
                crse_flux = cell_mp
                fac = HALF
             else
                crse_flux = HALF*(cell_pp + cell_mp)
                fac = HALF
             end if
          else if (lo_k_neu) then
             if (lo_i_not) then
                crse_flux = cell_pp
                fac = HALF
             else if (lo_i_neu) then
                crse_flux = FOUR*cell_pp
                fac = ONE
             else if (hi_i_not) then
                crse_flux = cell_mp
                fac = HALF
             else if (hi_i_neu) then
                crse_flux = FOUR*cell_mp
                fac = ONE
             else
                crse_flux = TWO*(cell_pp + cell_mp)
                fac = ONE
             end if
          else if (hi_k_not) then
             if (lo_i_not) then
                crse_flux = THIRD*cell_pm
                fac = THIRD
             else if (lo_i_neu) then
                crse_flux = cell_pm
                fac = HALF
             else if (hi_i_not) then
                crse_flux = THIRD*cell_mm
                fac = THIRD
             else if (hi_i_neu) then
                crse_flux = cell_mm
                fac = HALF
             else
                crse_flux = HALF*(cell_pm  + cell_mm)
                fac = HALF
             end if
          else if (hi_k_neu) then
             if (lo_i_not) then
                crse_flux = cell_pm
                fac = HALF
             else if (lo_i_neu) then
                crse_flux = FOUR*cell_pm
                fac = ONE
             else if (hi_i_not) then
                crse_flux = cell_mm
                fac = HALF
             else if (hi_i_neu) then
                crse_flux = FOUR*cell_mm
                fac = ONE
             else
                crse_flux = TWO*(cell_pm  + cell_mm)
                fac = ONE
             end if
          else
             if (lo_i_not) then
                crse_flux = HALF*(cell_pm  + cell_pp)
                fac = HALF
             else if (lo_i_neu) then
                crse_flux = TWO*(cell_pm  + cell_pp)
                fac = ONE
             else if (hi_i_not) then
                crse_flux = HALF*(cell_mm  + cell_mp)
                fac = HALF
             else if (hi_i_neu) then
                crse_flux = TWO*(cell_mm  + cell_mp)
                fac = ONE
             else
                crse_flux = cell_mm  + cell_mp + cell_pm + cell_pp
                fac = ONE
             end if
          end if

          rh(i,jj,k) = rh(i,jj,k) - crse_flux + fac*fine_flux(i,jj,k)
        end if

      end do
      end do
      !$OMP END PARALLEL DO

!   Lo/Hi k side
    else if (( side == -3) .or. (side == 3) ) then
 
      if (side == -3) then
        k    = kk
      else
        k    = kk-1
      end if

      !$OMP PARALLEL DO PRIVATE(i,j,lo_i_not,hi_i_not,lo_j_not,hi_j_not) &
      !$OMP PRIVATE(lo_i_neu,hi_i_neu,lo_j_neu,hi_j_neu) &
      !$OMP PRIVATE(cell_mm,cell_pm,cell_mp,cell_pp,crse_flux,fac)
      do j = lo(2),hi(2)
      do i = lo(1),hi(1)

        if (bc_dirichlet(mm(ir(1)*i,ir(2)*j,ir(3)*kk),1,0)) then

          lo_i_not = .false.
          hi_i_not = .false.
          lo_i_neu = .false.
          hi_i_neu = .false.
          lo_j_not = .false.
          hi_j_not = .false.
          lo_j_neu = .false.
          hi_j_neu = .false.

          if (i == loflx(1)) then
             if (.not. bc_neumann(mm(ir(1)*i,ir(2)*j,ir(3)*kk),1,-1)) lo_i_not = .true.
             if (bc_neumann(mm(ir(1)*i,ir(2)*j,ir(3)*kk),1,-1))       lo_i_neu = .true.
          end if

          if (i == hiflx(1)) then
             if (.not. bc_neumann(mm(ir(1)*i,ir(2)*j,ir(3)*kk),1,+1)) hi_i_not = .true.
             if (bc_neumann(mm(ir(1)*i,ir(2)*j,ir(3)*kk),1,+1))       hi_i_neu = .true.
          end if

          if (j == loflx(2)) then
             if (.not. bc_neumann(mm(ir(1)*i,ir(2)*j,ir(3)*kk),2,-1)) lo_j_not = .true.
             if (bc_neumann(mm(ir(1)*i,ir(2)*j,ir(3)*kk),2,-1))       lo_j_neu = .true.
          end if

          if (j == hiflx(2)) then
             if (.not. bc_neumann(mm(ir(1)*i,ir(2)*j,ir(3)*kk),2,+1)) hi_j_not = .true.
             if (bc_neumann(mm(ir(1)*i,ir(2)*j,ir(3)*kk),2,+1))       hi_j_neu = .true.
          end if

          cell_pp = rc(i  ,j  ,k) * EIGHTH
          cell_pm = rc(i  ,j-1,k) * EIGHTH
          cell_mp = rc(i-1,j  ,k) * EIGHTH
          cell_mm = rc(i-1,j-1,k) * EIGHTH

          if (lo_j_not) then
             if (lo_i_not) then
                crse_flux = THIRD*cell_pp
                fac = THIRD
             else if (lo_i_neu) then
                crse_flux = cell_pp
                fac = HALF
             else if (hi_i_not) then
                crse_flux = THIRD*cell_mp
                fac = THIRD
             else if (hi_i_neu) then
                crse_flux = cell_mp
                fac = HALF
             else
                crse_flux = HALF*(cell_pp + cell_mp)
                fac = HALF
             end if
          else if (lo_j_neu) then
             if (lo_i_not) then
                crse_flux = cell_pp
                fac = HALF
             else if (lo_i_neu) then
                crse_flux = FOUR*cell_pp
                fac = ONE
             else if (hi_i_not) then
                crse_flux = cell_mp
                fac = HALF
             else if (hi_i_neu) then
                crse_flux = FOUR*cell_mp
                fac = ONE
             else
                crse_flux = TWO*(cell_pp + cell_mp)
                fac = ONE
             end if
          else if (hi_j_not) then
             if (lo_i_not) then
                crse_flux = THIRD*cell_pm
                fac = THIRD
             else if (lo_i_neu) then
                crse_flux = cell_pm
                fac = HALF
             else if (hi_i_not) then
                crse_flux = THIRD*cell_mm
                fac = THIRD
             else if (hi_i_neu) then
                crse_flux = cell_mm
                fac = HALF
             else
                crse_flux = HALF*(cell_pm  + cell_mm)
                fac = HALF
             end if
          else if (hi_j_neu) then
             if (lo_i_not) then
                crse_flux = cell_pm
                fac = HALF
             else if (lo_i_neu) then
                crse_flux = FOUR*cell_pm
                fac = ONE
             else if (hi_i_not) then
                crse_flux = cell_mm
                fac = HALF
             else if (hi_i_neu) then
                crse_flux = FOUR*cell_mm
                fac = ONE
             else
                crse_flux = TWO*(cell_pm  + cell_mm)
                fac = ONE
             end if
          else
             if (lo_i_not) then
                crse_flux = HALF*(cell_pm  + cell_pp)
                fac = HALF
             else if (lo_i_neu) then
                crse_flux = TWO*(cell_pm  + cell_pp)
                fac = ONE
             else if (hi_i_not) then
                crse_flux = HALF*(cell_mm  + cell_mp)
                fac = HALF
             else if (hi_i_neu) then
                crse_flux = TWO*(cell_mm  + cell_mp)
                fac = ONE
             else
                crse_flux = cell_mm  + cell_mp + cell_pm + cell_pp
                fac = ONE
             end if
          end if
  
          rh(i,j,kk) = rh(i,j,kk) - crse_flux + fac*fine_flux(i,j,kk)
        end if

      end do
      end do
      !$OMP END PARALLEL DO

    end if

  end subroutine ml_interface_rhcc_3d

end module nodal_divu_module


