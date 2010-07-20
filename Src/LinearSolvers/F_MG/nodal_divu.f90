module nodal_divu_module

  use bl_constants_module
  use bndry_reg_module
  use mg_tower_module
  use ml_restriction_module

  implicit none

contains

!   ********************************************************************************************* !

    subroutine divu(nlevs,mgt,unew,rh,ref_ratio,nodal)

      integer        , intent(in   ) :: nlevs
      type(mg_tower) , intent(inout) :: mgt(:)
      type(multifab) , intent(inout) :: unew(:)
      type(multifab) , intent(inout) :: rh(:)
      integer        , intent(in   ) :: ref_ratio(:,:)
      logical        , intent(in   ) :: nodal(:)

      real(kind=dp_t), pointer :: unp(:,:,:,:) 
      real(kind=dp_t), pointer :: rhp(:,:,:,:) 
      integer        , pointer ::  mp(:,:,:,:) 

      integer :: i,n,dm,ng
      integer :: mglev_fine
      type(      box) :: pdc
      type(   layout) :: la_crse,la_fine
      type(bndry_reg) :: brs_flx

      dm = get_dim(unew(nlevs))
      ng = nghost(unew(nlevs))

!     Create the regular single-level divergence.
      do n = 1, nlevs
         mglev_fine = mgt(n)%nlevels
         call multifab_fill_boundary(unew(n))
         do i = 1, nboxes(unew(n))
            if ( multifab_remote(unew(n), i) ) cycle
            unp => dataptr(unew(n), i)
            rhp => dataptr(rh(n)  , i)
            mp   => dataptr(mgt(n)%mm(mglev_fine),i)
            select case (dm)
               case (1)
                 call divu_1d(unp(:,1,1,1), rhp(:,1,1,1), &
                               mp(:,1,1,1), mgt(n)%dh(:,mglev_fine), &
                              mgt(n)%face_type(i,:,:), ng)
               case (2)
                 call divu_2d(unp(:,:,1,:), rhp(:,:,1,1), &
                               mp(:,:,1,1), mgt(n)%dh(:,mglev_fine), &
                              mgt(n)%face_type(i,:,:), ng)
               case (3)
                 call divu_3d(unp(:,:,:,:), rhp(:,:,:,1), &
                               mp(:,:,:,1), mgt(n)%dh(:,mglev_fine), &
                              mgt(n)%face_type(i,:,:), ng)
            end select
         end do
      end do

!     Modify the divu above at coarse-fine interfaces.
      do n = nlevs,2,-1

         la_fine = get_layout(unew(n))
         mglev_fine = mgt(n)%nlevels

         la_crse = get_layout(unew(n-1))
         pdc = get_pd(la_crse)

         call bndry_reg_rr_build(brs_flx,la_fine,la_crse, ref_ratio(n-1,:), &
                                 pdc, nodal = nodal, other = .false.)
         call crse_fine_divu(n,nlevs,rh(n-1),unew,brs_flx,ref_ratio(n-1,:),mgt)
         call bndry_reg_destroy(brs_flx)
      end do

    end subroutine divu

!   ********************************************************************************************* !

    subroutine divu_1d(u,rh,mm,dx,face_type,ng)

      integer        , intent(in   ) :: ng
      real(kind=dp_t), intent(inout) ::  u(-ng:)
      real(kind=dp_t), intent(inout) :: rh(-1:)
      integer        , intent(inout) :: mm(0:)
      real(kind=dp_t), intent(in   ) :: dx(:)
      integer        , intent(in   ) :: face_type(:,:)

      integer         :: i,nx

      nx = size(rh,dim=1) - 3

      rh = ZERO

      do i = 0,nx
         if (.not.bc_dirichlet(mm(i),1,0)) then 
           rh(i) = (u(i)-u(i-1)) / dx(1) 
         end if
      end do

      if (face_type(1,1) == BC_NEU) rh( 0) = TWO*rh( 0)
      if (face_type(1,2) == BC_NEU) rh(nx) = TWO*rh(nx)

    end subroutine divu_1d

!   ********************************************************************************************* !

    subroutine divu_2d(u,rh,mm,dx,face_type,ng)

      integer        , intent(in   ) :: ng
      real(kind=dp_t), intent(inout) ::  u(-ng:,-ng:,1:)
      real(kind=dp_t), intent(inout) :: rh(-1:,-1:)
      integer        , intent(inout) :: mm(0:,0:)
      real(kind=dp_t), intent(in   ) :: dx(:)
      integer        , intent(in   ) :: face_type(:,:)

      integer         :: i,j,nx,ny

      nx = size(rh,dim=1) - 3
      ny = size(rh,dim=2) - 3

      rh = ZERO

      do j = 0,ny
      do i = 0,nx
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
      if (face_type(1,2) == BC_NEU) rh(nx,:) = TWO*rh(nx,:)
      if (face_type(2,1) == BC_NEU) rh(:, 0) = TWO*rh(:, 0)
      if (face_type(2,2) == BC_NEU) rh(:,ny) = TWO*rh(:,ny)

    end subroutine divu_2d

!   ********************************************************************************************* !

    subroutine divu_3d(u,rh,mm,dx,face_type,ng)

      integer        , intent(in   ) :: ng
      real(kind=dp_t), intent(inout) ::  u(-ng:,-ng:,-ng:,1:)
      real(kind=dp_t), intent(inout) :: rh(-1:,-1:,-1:)
      integer        , intent(inout) :: mm(0:,0:,0:)
      real(kind=dp_t), intent(in   ) :: dx(:)
      integer        , intent(in   ) :: face_type(:,:)

      integer         :: i,j,k,nx,ny,nz

      nx = size(rh,dim=1) - 3
      ny = size(rh,dim=2) - 3
      nz = size(rh,dim=3) - 3

      rh = ZERO

      !$OMP PARALLEL DO PRIVATE(i,j,k)
      do k = 0,nz
      do j = 0,ny
      do i = 0,nx
         if (.not. bc_dirichlet(mm(i,j,k),1,0)) then
           rh(i,j,k) = (u(i  ,j,k  ,1) + u(i  ,j-1,k  ,1) &
                       +u(i  ,j,k-1,1) + u(i  ,j-1,k-1,1) &
                       -u(i-1,j,k  ,1) - u(i-1,j-1,k  ,1) &
                       -u(i-1,j,k-1,1) - u(i-1,j-1,k-1,1)) / dx(1) + &
                       (u(i,j  ,k  ,2) + u(i-1,j  ,k  ,2) &
                       +u(i,j  ,k-1,2) + u(i-1,j  ,k-1,2) &
                       -u(i,j-1,k  ,2) - u(i-1,j-1,k  ,2) &
                       -u(i,j-1,k-1,2) - u(i-1,j-1,k-1,2)) / dx(2) + &
                       (u(i,j  ,k  ,3) + u(i-1,j  ,k  ,3) &
                       +u(i,j-1,k  ,3) + u(i-1,j-1,k  ,3) &
                       -u(i,j  ,k-1,3) - u(i-1,j  ,k-1,3) &
                       -u(i,j-1,k-1,3) - u(i-1,j-1,k-1,3)) / dx(3)
           rh(i,j,k) = FOURTH * rh(i,j,k)
         end if
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

    end subroutine divu_3d

!   ********************************************************************************************* !

    subroutine crse_fine_divu(n_fine,nlevs,rh_crse,u,brs_flx,ref_ratio,mgt)

      use ml_util_module, only : ml_fine_contrib

      integer        , intent(in   ) :: n_fine,nlevs
      type(multifab) , intent(inout) :: rh_crse
      type(multifab) , intent(inout) :: u(:)
      type(bndry_reg), intent(inout) :: brs_flx
      integer        , intent(in   ) :: ref_ratio(:)
      type(mg_tower) , intent(in   ) :: mgt(:)

      real(kind=dp_t), pointer :: unp(:,:,:,:) 
      real(kind=dp_t), pointer :: rhp(:,:,:,:) 

      type(multifab) :: temp_rhs, temp_rhs_crse
      type(  layout) :: la_crse,la_fine
      type(     box) :: pdc
      integer :: i,dm,n_crse,ng
      integer :: mglev_fine, mglev_crse
      logical :: nodal(get_dim(u(n_fine)))

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
       do i = 1, nboxes(u(n_fine))
          if ( multifab_remote(u(n_fine), i) ) cycle
          unp => dataptr(u(n_fine), i)
          rhp => dataptr( temp_rhs, i)
          select case (dm)
             case (1)
               call grid_divu_1d(unp(:,1,1,1), rhp(:,1,1,1), mgt(n_fine)%dh(:,mglev_fine), &
                                 mgt(n_fine)%face_type(i,:,:), ng)
             case (2)
               call grid_divu_2d(unp(:,:,1,:), rhp(:,:,1,1), mgt(n_fine)%dh(:,mglev_fine), &
                                 mgt(n_fine)%face_type(i,:,:), ng)
             case (3)
               call grid_divu_3d(unp(:,:,:,:), rhp(:,:,:,1), mgt(n_fine)%dh(:,mglev_fine), &
                                 mgt(n_fine)%face_type(i,:,:), ng)
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
                                   pdc,ref_ratio, -i)
         call ml_crse_divu_contrib(temp_rhs_crse, brs_flx%bmf(i,1), u(n_crse), &
                                   mgt(n_fine)%mm(mglev_fine), mgt(n_crse)%dh(:,mglev_crse), &
                                   pdc,ref_ratio, +i)
      end do

      call multifab_plus_plus(rh_crse,temp_rhs_crse)
      call periodic_add_copy(rh_crse,temp_rhs_crse,synced=.true.)

      call multifab_destroy(temp_rhs)
      call multifab_destroy(temp_rhs_crse)

    end subroutine crse_fine_divu

!   ********************************************************************************************* !

    subroutine grid_divu_1d(u,rh,dx,face_type,ng)

      integer        , intent(in   ) :: ng
      real(kind=dp_t), intent(inout) ::  u(-ng:)
      real(kind=dp_t), intent(inout) :: rh(-1:)
      real(kind=dp_t), intent(in   ) :: dx(:)
      integer        , intent(in   ) :: face_type(:,:)

      integer :: i,nx
      nx = size(rh,dim=1) - 3

      u(-1) = ZERO
      u(nx) = ZERO

      do i = 0,nx
         rh(i) = (u(i)-u(i-1)) / dx(1)
      end do

      if (face_type(1,1) == BC_NEU) rh( 0) = TWO*rh( 0)
      if (face_type(1,2) == BC_NEU) rh(nx) = TWO*rh(nx)

    end subroutine grid_divu_1d

!   ********************************************************************************************* !

    subroutine grid_divu_2d(u,rh,dx,face_type,ng)

      integer        , intent(in   ) :: ng
      real(kind=dp_t), intent(inout) ::  u(-ng:,-ng:,1:)
      real(kind=dp_t), intent(inout) :: rh(-1:,-1:)
      real(kind=dp_t), intent(in   ) :: dx(:)
      integer        , intent(in   ) :: face_type(:,:)

      integer :: i,j,nx,ny
      nx = size(rh,dim=1) - 3
      ny = size(rh,dim=2) - 3

      u(-1,:,:) = ZERO
      u(nx,:,:) = ZERO
      u(:,-1,:) = ZERO
      u(:,ny,:) = ZERO

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

!   ********************************************************************************************* !

    subroutine grid_divu_3d(u,rh,dx,face_type,ng)

      integer        , intent(in   ) :: ng
      real(kind=dp_t), intent(inout) ::  u(-ng:,-ng:,-ng:,1:)
      real(kind=dp_t), intent(inout) :: rh(-1:,-1:,-1:)
      real(kind=dp_t), intent(in   ) :: dx(:)
      integer        , intent(in   ) :: face_type(:,:)

      integer :: i,j,k,nx,ny,nz

      nx = size(rh,dim=1) - 3
      ny = size(rh,dim=2) - 3
      nz = size(rh,dim=3) - 3

      u(-1,:,:,:) = ZERO
      u(nx,:,:,:) = ZERO
      u(:,-1,:,:) = ZERO
      u(:,ny,:,:) = ZERO
      u(:,:,-1,:) = ZERO
      u(:,:,nz,:) = ZERO

      !$OMP PARALLEL DO PRIVATE(i,j,k)
      do k = 0,nz
      do j = 0,ny
      do i = 0,nx
         rh(i,j,k) = (u(i  ,j,k  ,1) + u(i  ,j-1,k  ,1) &
                     +u(i  ,j,k-1,1) + u(i  ,j-1,k-1,1) &
                     -u(i-1,j,k  ,1) - u(i-1,j-1,k  ,1) &
                     -u(i-1,j,k-1,1) - u(i-1,j-1,k-1,1)) / dx(1) + &
                     (u(i,j  ,k  ,2) + u(i-1,j  ,k  ,2) &
                     +u(i,j  ,k-1,2) + u(i-1,j  ,k-1,2) &
                     -u(i,j-1,k  ,2) - u(i-1,j-1,k  ,2) &
                     -u(i,j-1,k-1,2) - u(i-1,j-1,k-1,2)) / dx(2) + &
                     (u(i,j  ,k  ,3) + u(i-1,j  ,k  ,3) &
                     +u(i,j-1,k  ,3) + u(i-1,j-1,k  ,3) &
                     -u(i,j  ,k-1,3) - u(i-1,j  ,k-1,3) &
                     -u(i,j-1,k-1,3) - u(i-1,j-1,k-1,3)) / dx(3)
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

!   ********************************************************************************************* !

    !
    ! TODO - cache the communication pattern here in the flux layout?
    !
    subroutine ml_crse_divu_contrib(rh, flux, u, mm, dx, crse_domain, ir, side)
     type(multifab), intent(inout) :: rh
     type(multifab), intent(inout) :: flux
     type(multifab), intent(in   ) :: u
     type(imultifab),intent(in   ) :: mm
     real(kind=dp_t),intent(in   ) :: dx(:)
     type(box)      ,intent(in   ) :: crse_domain
     integer        ,intent(in   ) :: ir(:)
     integer        ,intent(in   ) :: side

     type(box) :: fbox, ubox, mbox, isect

     integer   :: lo (get_dim(rh)), hi (get_dim(rh)), lou(get_dim(rh)), dims(4), dm
     integer   :: lof(get_dim(rh)), hif(get_dim(rh)), lor(get_dim(rh)), lom(get_dim(rh))
     integer   :: lodom(get_dim(rh)), hidom(get_dim(rh)), dir, i, j, k, proc
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

     do j = 1, nboxes(u)
       ubox = box_nodalize(get_ibox(u,j),nodal)
       lou  = lwb(get_pbox(u,j))
       lor  = lwb(get_pbox(rh,j))

       bi => layout_get_box_intersector(flux_la, ubox)

       do k = 1, size(bi)

          i = bi(k)%i

          if ( remote(flux,i) .and. remote(u,j) ) cycle
          
          fbox  = get_ibox(flux,i)
          isect = bi(k)%bx
          lof   = lwb(fbox)
          hif   = upb(fbox)

          if ( (lof(dir) == lodom(dir) .or. lof(dir) == hidom(dir)) .and. &
               .not. pmask(dir) ) cycle

          lo = lwb(isect)
          hi = upb(isect)

          if ( local(flux,i) .and. local(u,j) ) then
             lom  =  lwb(get_pbox(mm,i))
             fp   => dataptr(flux,i)
             mp   => dataptr(mm  ,i)
             up   => dataptr(u   ,j)
             rp   => dataptr(rh  ,j)
             select case (dm)
             case (1)
                call ml_interface_1d_divu(rp(:,1,1,1), lor, &
                     fp(:,1,1,1), lof, &
                     up(:,1,1,1), lou, mp(:,1,1,1), lom, lo, ir, side, dx)
             case (2)
                call ml_interface_2d_divu(rp(:,:,1,1), lor, &
                     fp(:,:,1,1), lof, lof, hif, &
                     up(:,:,1,:), lou, mp(:,:,1,1), lom, lo, hi, ir, side, dx)
             case (3)
                call ml_interface_3d_divu(rp(:,:,:,1), lor, &
                     fp(:,:,:,1), lof, lof, hif,  &
                     up(:,:,:,:), lou, mp(:,:,:,1), lom, lo, hi, ir, side, dx)
             end select

          else if ( local(flux,i) ) then
             !
             ! Must send flux & mm.
             !
             mbox =  intersection(refine(isect,ir), get_pbox(mm,i))
             fp   => dataptr(flux, i, isect, 1, ncomp(flux))
             mp   => dataptr(mm,   i, mbox,  1, ncomp(mm))
             proc =  get_proc(get_layout(u), j)
             call parallel_send(fp, proc, tag)
             call parallel_send(mp, proc, tag)

          else if ( local(u,j) ) then
             !
             ! Must receive flux & mm.
             !
             proc = get_proc(flux_la, i)
             mbox = intersection(refine(isect,ir), get_pbox(mm,i))
             lom  = lwb(mbox)
             dims(1:dm) = extent(isect)
             allocate(fp(dims(1),dims(2),dims(3),ncomp(flux)))
             dims(1:dm) = extent(mbox)
             allocate(mp(dims(1),dims(2),dims(3),ncomp(mm)))
             call parallel_recv(fp, proc, tag)
             call parallel_recv(mp, proc, tag)
             up => dataptr(u  ,j)
             rp => dataptr(rh ,j)
             select case (dm)
             case (1)
                call ml_interface_1d_divu(rp(:,1,1,1), lor, &
                     fp(:,1,1,1), lo, &
                     up(:,1,1,1), lou, mp(:,1,1,1), lom, lo, ir, side, dx)
             case (2)
                call ml_interface_2d_divu(rp(:,:,1,1), lor, &
                     fp(:,:,1,1), lo, lof, hif, &
                     up(:,:,1,:), lou, mp(:,:,1,1), lom, lo, hi, ir, side, dx)
             case (3)
                call ml_interface_3d_divu(rp(:,:,:,1), lor, &
                     fp(:,:,:,1), lo, lof, hif, &
                     up(:,:,:,:), lou, mp(:,:,:,1), lom, lo, hi, ir, side, dx)
             end select
             deallocate(fp,mp)
          end if
       end do
       deallocate(bi)
    end do
    call destroy(bpt)

   end subroutine ml_crse_divu_contrib

!   ********************************************************************************************* !

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

    fac = 1.0_dp_t / dble(ir(1))

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

!   ********************************************************************************************* !

    subroutine ml_interface_2d_divu(rh, lor, fine_flux, lof, loflx, hiflx, uc, loc, &
                                    mm, lom, lo, hi, ir, side, dx)
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

    integer :: i, j
    real (kind = dp_t) :: crse_flux,fac

    i = lo(1)
    j = lo(2)

!   NOTE: THESE STENCILS ONLY WORK FOR DX == DY.

!   NOTE: MM IS ON THE FINE GRID, NOT THE CRSE

    fac = 1.0_dp_t / dble(ir(1)*ir(2))

!   Lo i side
    if (side == -1) then

       do j = lo(2),hi(2)

          if (bc_dirichlet(mm(ir(1)*i,ir(2)*j),1,0)) then
             if (j == loflx(2)) then
                if (bc_neumann(mm(ir(1)*i,ir(2)*j),2,-1)) then
                   crse_flux =        (uc(i,j,1)/dx(1) + uc(i,j,2)/dx(2))
                else 
                   crse_flux = FOURTH*(uc(i,j,1)/dx(1) + uc(i,j,2)/dx(2))
                end if
             else if (j == hiflx(2)) then
                if (bc_neumann(mm(ir(1)*i,ir(2)*j),2,+1)) then
                   crse_flux =        (uc(i,j-1,1)/dx(1) - uc(i,j-1,2)/dx(2))
                else 
                   crse_flux = FOURTH*(uc(i,j-1,1)/dx(1) - uc(i,j-1,2)/dx(2))
                end if
             else
                crse_flux = (HALF*(uc(i,j,1) + uc(i,j-1,1))/dx(1) &
                            +HALF*(uc(i,j,2) - uc(i,j-1,2))/dx(2))
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
                   crse_flux =        (-uc(i-1,j,1)/dx(1) + uc(i-1,j,2)/dx(2))
                else
                   crse_flux = FOURTH*(-uc(i-1,j,1)/dx(1) + uc(i-1,j,2)/dx(2))
                end if 
             else if (j == hiflx(2)) then
                if (bc_neumann(mm(ir(1)*i,ir(2)*j),2,+1)) then
                   crse_flux =        (-uc(i-1,j-1,1)/dx(1) - uc(i-1,j-1,2)/dx(2))
                else 
                   crse_flux = FOURTH*(-uc(i-1,j-1,1)/dx(1) - uc(i-1,j-1,2)/dx(2)) 
                end if
             else
                crse_flux = (HALF*(-uc(i-1,j,1)-uc(i-1,j-1,1))/dx(1)  &
                            +HALF*( uc(i-1,j,2)-uc(i-1,j-1,2))/dx(2))
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
                   crse_flux =        (uc(i,j,1)/dx(1) + uc(i,j,2)/dx(2))
                else 
                   crse_flux = FOURTH*(uc(i,j,1)/dx(1) + uc(i,j,2)/dx(2))
                end if
             else if (i == hiflx(1)) then
                if (bc_neumann(mm(ir(1)*i,ir(2)*j),1,+1)) then
                   crse_flux =        (-uc(i-1,j,1)/dx(1) + uc(i-1,j,2)/dx(2))
                else 
                   crse_flux = FOURTH*(-uc(i-1,j,1)/dx(1) + uc(i-1,j,2)/dx(2))
                end if
             else
                crse_flux = (HALF*(uc(i,j,1)-uc(i-1,j,1))/dx(1)  &
                            +HALF*(uc(i,j,2)+uc(i-1,j,2))/dx(2))
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
                   crse_flux =        (uc(i,j-1,1)/dx(1) - uc(i,j-1,2)/dx(2))
                else
                   crse_flux = FOURTH*(uc(i,j-1,1)/dx(1) - uc(i,j-1,2)/dx(2))
                end if

             else if (i == hiflx(1)) then
                if (bc_neumann(mm(ir(1)*i,ir(2)*j),1,+1)) then
                   crse_flux =        (-uc(i-1,j-1,1)/dx(1) - uc(i-1,j-1,2)/dx(2))
                else 
                   crse_flux = FOURTH*(-uc(i-1,j-1,1)/dx(1) - uc(i-1,j-1,2)/dx(2))
                end if
             else
                crse_flux = (HALF*( uc(i,j-1,1)-uc(i-1,j-1,1))/dx(1) &
                            +HALF*(-uc(i,j-1,2)-uc(i-1,j-1,2))/dx(2))
             end if

             rh(i,j) = rh(i,j) - crse_flux + fac*fine_flux(i,j)

          end if
       end do

    end if

  end subroutine ml_interface_2d_divu

!   ********************************************************************************************* !

    subroutine ml_interface_3d_divu(rh, lor, fine_flux, lof, loflx, hiflx, uc, loc, &
                                    mm, lom, lo, hi, ir, side, dx)
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

    integer :: i, j, k, ii, jj, kk, ifac
    logical :: lo_i_neu,lo_j_neu,lo_k_neu,hi_i_neu,hi_j_neu,hi_k_neu
    logical :: lo_i_not,lo_j_not,lo_k_not,hi_i_not,hi_j_not,hi_k_not
    real (kind = dp_t) :: cell_pp,cell_mp,cell_pm,cell_mm
    real (kind = dp_t) :: crse_flux,fac

    ii = lo(1)
    jj = lo(2)
    kk = lo(3)

!   NOTE: THESE STENCILS ONLY WORK FOR DX == DY == DZ.

!   NOTE: MM IS ON THE FINE GRID, NOT THE CRSE

    fac = 1.0_dp_t / (ir(1)*ir(2)*ir(3))

!   Lo/Hi i side
    if (( side == -1) .or. (side == 1) ) then
 
      if (side == -1) then
        i    = ii
        ifac = 1
      else
        i    = ii-1
        ifac = -1
      end if

      do k = lo(3),hi(3)
      do j = lo(2),hi(2)

        if (bc_dirichlet(mm(ir(1)*ii,ir(2)*j,ir(3)*k),1,0)) then

          cell_pp =  (uc(i,j  ,k  ,1) )/dx(1) * ifac &
                   + (uc(i,j  ,k  ,2) )/dx(2) &
                   + (uc(i,j  ,k  ,3) )/dx(3)
          cell_pm =  (uc(i,j  ,k-1,1) )/dx(1) * ifac &
                   + (uc(i,j  ,k-1,2) )/dx(2) &
                   - (uc(i,j  ,k-1,3) )/dx(3)
          cell_mp =  (uc(i,j-1,k  ,1) )/dx(1) * ifac &
                   - (uc(i,j-1,k  ,2) )/dx(2) &
                   + (uc(i,j-1,k  ,3) )/dx(3)
          cell_mm =  (uc(i,j-1,k-1,1) )/dx(1) * ifac &
                   - (uc(i,j-1,k-1,2) )/dx(2) &
                   - (uc(i,j-1,k-1,3) )/dx(3)

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
          end if
          if (j == hiflx(2)) then
             if (.not. bc_neumann(mm(ir(1)*ii,ir(2)*j,ir(3)*k),2,+1)) hi_j_not = .true.
          end if
          if (j == loflx(2)) then
             if (bc_neumann(mm(ir(1)*ii,ir(2)*j,ir(3)*k),2,-1)) lo_j_neu = .true.
          end if
          if (j == hiflx(2)) then
             if (bc_neumann(mm(ir(1)*ii,ir(2)*j,ir(3)*k),2,+1)) hi_j_neu = .true.
          end if
          if (k == loflx(3)) then
             if (.not. bc_neumann(mm(ir(1)*ii,ir(2)*j,ir(3)*k),3,-1)) lo_k_not = .true.
          end if
          if (k == hiflx(3)) then
             if (.not. bc_neumann(mm(ir(1)*ii,ir(2)*j,ir(3)*k),3,+1)) hi_k_not = .true.
          end if
          if (k == loflx(3)) then
             if (bc_neumann(mm(ir(1)*ii,ir(2)*j,ir(3)*k),3,-1)) lo_k_neu = .true.
          end if
          if (k == hiflx(3)) then
             if (bc_neumann(mm(ir(1)*ii,ir(2)*j,ir(3)*k),3,+1)) hi_k_neu = .true.
          end if

          if (lo_k_not) then
             if (lo_j_not) then
                crse_flux = THIRD*cell_pp
             else if (lo_j_neu) then
                crse_flux = cell_pp
             else if (hi_j_not) then
                crse_flux = THIRD*cell_mp
             else if (hi_j_neu) then
                crse_flux = cell_mp
             else
                crse_flux = HALF*(cell_pp + cell_mp)
             end if
          else if (lo_k_neu) then
             if (lo_j_not) then
                crse_flux = cell_pp
             else if (lo_j_neu) then
                crse_flux = FOUR*cell_pp
             else if (hi_j_not) then
                crse_flux = cell_mp
             else if (hi_j_neu) then
                crse_flux = FOUR*cell_mp
             else
                crse_flux = TWO*(cell_pp + cell_mp)
             end if
          else if (hi_k_not) then
             if (lo_j_not) then
                crse_flux = THIRD*cell_pm
             else if (lo_j_neu) then
                crse_flux = cell_pm
             else if (hi_j_not) then
                crse_flux = THIRD*cell_mm
             else if (hi_j_neu) then
                crse_flux = cell_mm
             else
                crse_flux = HALF*(cell_pm  + cell_mm)
             end if
          else if (hi_k_neu) then
             if (lo_j_not) then
                crse_flux = cell_pm
             else if (lo_j_neu) then
                crse_flux = FOUR*cell_pm
             else if (hi_j_not) then
                crse_flux = cell_mm
             else if (hi_j_neu) then
                crse_flux = FOUR*cell_mm
             else
                crse_flux = TWO*(cell_pm  + cell_mm)
             end if
          else
             if (lo_j_not) then
                crse_flux = HALF*(cell_pm  + cell_pp)
             else if (lo_j_neu) then
                crse_flux = TWO*(cell_pm  + cell_pp)
             else if (hi_j_not) then
                crse_flux = HALF*(cell_mm  + cell_mp)
             else if (hi_j_neu) then
                crse_flux = TWO*(cell_mm  + cell_mp)
             else
                crse_flux = cell_mm  + cell_mp + cell_pm + cell_pp
             end if
          end if

          rh(ii,j,k) = rh(ii,j,k) - FOURTH*crse_flux + fac*fine_flux(ii,j,k)
        end if

      end do
      end do

!   Lo/Hi j side
    else if (( side == -2) .or. (side == 2) ) then
 
      if (side == -2) then
        j    = jj
        ifac = 1
      else
        j    = jj-1
        ifac = -1
      end if

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
          end if
          if (i == hiflx(1)) then
             if (.not. bc_neumann(mm(ir(1)*i,ir(2)*jj,ir(3)*k),1,+1)) hi_i_not = .true.
          end if
          if (i == loflx(1)) then
             if (bc_neumann(mm(ir(1)*i,ir(2)*jj,ir(3)*k),1,-1)) lo_i_neu = .true.
          end if
          if (i == hiflx(1)) then
             if (bc_neumann(mm(ir(1)*i,ir(2)*jj,ir(3)*k),1,+1)) hi_i_neu = .true.
          end if
          if (k == loflx(3)) then
             if (.not. bc_neumann(mm(ir(1)*i,ir(2)*jj,ir(3)*k),3,-1)) lo_k_not = .true.
          end if
          if (k == hiflx(3)) then
             if (.not. bc_neumann(mm(ir(1)*i,ir(2)*jj,ir(3)*k),3,+1)) hi_k_not = .true.
          end if
          if (k == loflx(3)) then
             if (bc_neumann(mm(ir(1)*i,ir(2)*jj,ir(3)*k),3,-1)) lo_k_neu = .true.
          end if
          if (k == hiflx(3)) then
             if (bc_neumann(mm(ir(1)*i,ir(2)*jj,ir(3)*k),3,+1)) hi_k_neu = .true.
          end if

          cell_pp =  (uc(i  ,j,k  ,1) )/dx(1) &
                    +(uc(i  ,j,k  ,2) )/dx(2) * ifac &
                    +(uc(i  ,j,k  ,3) )/dx(3)
          cell_pm =  (uc(i  ,j,k-1,1) )/dx(1) &
                    +(uc(i  ,j,k-1,2) )/dx(2) * ifac &
                    -(uc(i  ,j,k-1,3) )/dx(3)
          cell_mp = -(uc(i-1,j,k  ,1) )/dx(1) &
                    +(uc(i-1,j,k  ,2) )/dx(2) * ifac &
                    +(uc(i-1,j,k  ,3) )/dx(3)
          cell_mm = -(uc(i-1,j,k-1,1) )/dx(1) &
                    +(uc(i-1,j,k-1,2) )/dx(2) * ifac &
                    -(uc(i-1,j,k-1,3) )/dx(3)

          if (lo_k_not) then
             if (lo_i_not) then
                crse_flux = THIRD*cell_pp
             else if (lo_i_neu) then
                crse_flux = cell_pp
             else if (hi_i_not) then
                crse_flux = THIRD*cell_mp
             else if (hi_i_neu) then
                crse_flux = cell_mp
             else
                crse_flux = HALF*(cell_pp + cell_mp)
             end if
          else if (lo_k_neu) then
             if (lo_i_not) then
                crse_flux = cell_pp
             else if (lo_i_neu) then
                crse_flux = FOUR*cell_pp
             else if (hi_i_not) then
                crse_flux = cell_mp
             else if (hi_i_neu) then
                crse_flux = FOUR*cell_mp
             else
                crse_flux = TWO*(cell_pp + cell_mp)
             end if
          else if (hi_k_not) then
             if (lo_i_not) then
                crse_flux = THIRD*cell_pm
             else if (lo_i_neu) then
                crse_flux = cell_pm
             else if (hi_i_not) then
                crse_flux = THIRD*cell_mm
             else if (hi_i_neu) then
                crse_flux = cell_mm
             else
                crse_flux = HALF*(cell_pm  + cell_mm)
             end if
          else if (hi_k_neu) then
             if (lo_i_not) then
                crse_flux = cell_pm
             else if (lo_i_neu) then
                crse_flux = FOUR*cell_pm
             else if (hi_i_not) then
                crse_flux = cell_mm
             else if (hi_i_neu) then
                crse_flux = FOUR*cell_mm
             else
                crse_flux = TWO*(cell_pm  + cell_mm)
             end if
          else
             if (lo_i_not) then
                crse_flux = HALF*(cell_pm  + cell_pp)
             else if (lo_i_neu) then
                crse_flux = TWO*(cell_pm  + cell_pp)
             else if (hi_i_not) then
                crse_flux = HALF*(cell_mm  + cell_mp)
             else if (hi_i_neu) then
                crse_flux = TWO*(cell_mm  + cell_mp)
             else
                crse_flux = cell_mm  + cell_mp + cell_pm + cell_pp
             end if
          end if

          rh(i,jj,k) = rh(i,jj,k) - FOURTH*crse_flux + fac*fine_flux(i,jj,k)
        end if

      end do
      end do

!   Lo/Hi k side
    else if (( side == -3) .or. (side == 3) ) then
 
      if (side == -3) then
        k    = kk
        ifac = 1
      else
        k    = kk-1
        ifac = -1
      end if

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
          end if
          if (i == hiflx(1)) then
             if (.not. bc_neumann(mm(ir(1)*i,ir(2)*j,ir(3)*kk),1,+1)) hi_i_not = .true.
          end if
          if (i == loflx(1)) then
             if (bc_neumann(mm(ir(1)*i,ir(2)*j,ir(3)*kk),1,-1)) lo_i_neu = .true.
          end if
          if (i == hiflx(1)) then
             if (bc_neumann(mm(ir(1)*i,ir(2)*j,ir(3)*kk),1,+1)) hi_i_neu = .true.
          end if
          if (j == loflx(2)) then
             if (.not. bc_neumann(mm(ir(1)*i,ir(2)*j,ir(3)*kk),2,-1)) lo_j_not = .true.
          end if
          if (j == hiflx(2)) then
             if (.not. bc_neumann(mm(ir(1)*i,ir(2)*j,ir(3)*kk),2,+1)) hi_j_not = .true.
          end if
          if (j == loflx(2)) then
             if (bc_neumann(mm(ir(1)*i,ir(2)*j,ir(3)*kk),2,-1)) lo_j_neu = .true.
          end if
          if (j == hiflx(2)) then
             if (bc_neumann(mm(ir(1)*i,ir(2)*j,ir(3)*kk),2,+1)) hi_j_neu = .true.
          end if

          cell_pp =  (uc(i  ,j  ,k,1) )/dx(1) &
                    +(uc(i  ,j  ,k,2) )/dx(2) &
                    +(uc(i  ,j  ,k,3) )/dx(3) * ifac
          cell_pm =  (uc(i  ,j-1,k,1) )/dx(1) &
                    -(uc(i  ,j-1,k,2) )/dx(2) &
                    +(uc(i  ,j-1,k,3) )/dx(3) * ifac
          cell_mp = -(uc(i-1,j  ,k,1) )/dx(1) &
                    +(uc(i-1,j  ,k,2) )/dx(2) &
                    +(uc(i-1,j  ,k,3) )/dx(3) * ifac
          cell_mm = -(uc(i-1,j-1,k,1) )/dx(1) &
                    -(uc(i-1,j-1,k,2) )/dx(2) &
                    +(uc(i-1,j-1,k,3) )/dx(3) * ifac

          if (lo_j_not) then
             if (lo_i_not) then
                crse_flux = THIRD*cell_pp
             else if (lo_i_neu) then
                crse_flux = cell_pp
             else if (hi_i_not) then
                crse_flux = THIRD*cell_mp
             else if (hi_i_neu) then
                crse_flux = cell_mp
             else
                crse_flux = HALF*(cell_pp + cell_mp)
             end if
          else if (lo_j_neu) then
             if (lo_i_not) then
                crse_flux = cell_pp
             else if (lo_i_neu) then
                crse_flux = FOUR*cell_pp
             else if (hi_i_not) then
                crse_flux = cell_mp
             else if (hi_i_neu) then
                crse_flux = FOUR*cell_mp
             else
                crse_flux = TWO*(cell_pp + cell_mp)
             end if
          else if (hi_j_not) then
             if (lo_i_not) then
                crse_flux = THIRD*cell_pm
             else if (lo_i_neu) then
                crse_flux = cell_pm
             else if (hi_i_not) then
                crse_flux = THIRD*cell_mm
             else if (hi_i_neu) then
                crse_flux = cell_mm
             else
                crse_flux = HALF*(cell_pm  + cell_mm)
             end if
          else if (hi_j_neu) then
             if (lo_i_not) then
                crse_flux = cell_pm
             else if (lo_i_neu) then
                crse_flux = FOUR*cell_pm
             else if (hi_i_not) then
                crse_flux = cell_mm
             else if (hi_i_neu) then
                crse_flux = FOUR*cell_mm
             else
                crse_flux = TWO*(cell_pm  + cell_mm)
             end if
          else
             if (lo_i_not) then
                crse_flux = HALF*(cell_pm  + cell_pp)
             else if (lo_i_neu) then
                crse_flux = TWO*(cell_pm  + cell_pp)
             else if (hi_i_not) then
                crse_flux = HALF*(cell_mm  + cell_mp)
             else if (hi_i_neu) then
                crse_flux = TWO*(cell_mm  + cell_mp)
             else
                crse_flux = cell_mm  + cell_mp + cell_pm + cell_pp
             end if
          end if
  
          rh(i,j,kk) = rh(i,j,kk) - FOURTH*crse_flux + fac*fine_flux(i,j,kk)
        end if

      end do
      end do

    end if

  end subroutine ml_interface_3d_divu

    subroutine subtract_divu_from_rh(nlevs,mgt,rh,divu_rhs)

      integer        , intent(in   ) :: nlevs
      type(mg_tower) , intent(inout) :: mgt(:)
      type(multifab) , intent(inout) :: rh(:)
      type(multifab) , intent(in   ) :: divu_rhs(:)

      real(kind=dp_t), pointer :: dp(:,:,:,:) 
      real(kind=dp_t), pointer :: rp(:,:,:,:) 
      integer        , pointer :: mp(:,:,:,:) 

      integer :: i,n,dm,ng_r,ng_d
      integer :: mglev_fine

      dm   = get_dim(rh(nlevs))
      ng_r = nghost(rh(nlevs))
      ng_d = nghost(divu_rhs(nlevs))

!     Create the regular single-level divergence.
      do n = 1, nlevs
         mglev_fine = mgt(n)%nlevels
         do i = 1, nboxes(rh(n))
            if ( multifab_remote(rh(n), i) ) cycle
            rp => dataptr(rh(n)      , i)
            dp => dataptr(divu_rhs(n), i)
            mp => dataptr(mgt(n)%mm(mglev_fine),i)
            select case (dm)
               case (1)
                 call subtract_divu_from_rh_1d(rp(:,1,1,1), ng_r, &
                                               dp(:,1,1,1), ng_d, &
                                               mp(:,1,1,1) )
               case (2)
                 call subtract_divu_from_rh_2d(rp(:,:,1,1), ng_r, &
                                               dp(:,:,1,1), ng_d, &
                                               mp(:,:,1,1) )
               case (3)
                 call subtract_divu_from_rh_3d(rp(:,:,:,1), ng_r, &
                                               dp(:,:,:,1), ng_d, &
                                               mp(:,:,:,1) )
            end select
         end do
      end do

    end subroutine subtract_divu_from_rh

!   ********************************************************************************************* !

    subroutine subtract_divu_from_rh_1d(rh,ng_rh,divu_rhs,ng_divu,mm)

      integer        , intent(in   ) :: ng_rh,ng_divu
      real(kind=dp_t), intent(inout) ::       rh(  -ng_rh:)
      real(kind=dp_t), intent(inout) :: divu_rhs(-ng_divu:)
      integer        , intent(inout) :: mm(0:)

      integer         :: i,nx

      nx = size(mm,dim=1) - 1

      do i = 0,nx
         if (.not. bc_dirichlet(mm(i),1,0)) &
           rh(i) = rh(i) - divu_rhs(i)
      end do

    end subroutine subtract_divu_from_rh_1d

!   ********************************************************************************************* !

    subroutine subtract_divu_from_rh_2d(rh,ng_rh,divu_rhs,ng_divu,mm)

      integer        , intent(in   ) :: ng_rh,ng_divu
      real(kind=dp_t), intent(inout) ::       rh(  -ng_rh:,  -ng_rh:)
      real(kind=dp_t), intent(inout) :: divu_rhs(-ng_divu:,-ng_divu:)
      integer        , intent(inout) :: mm(0:,0:)

      integer         :: i,j,nx,ny

      nx = size(mm,dim=1) - 1
      ny = size(mm,dim=2) - 1

      do j = 0,ny
      do i = 0,nx
         if (.not. bc_dirichlet(mm(i,j),1,0)) &
           rh(i,j) = rh(i,j) - divu_rhs(i,j)
      end do
      end do

    end subroutine subtract_divu_from_rh_2d

!   ********************************************************************************************* !

    subroutine subtract_divu_from_rh_3d(rh,ng_rh,divu_rhs,ng_divu,mm)

      integer        , intent(in   ) :: ng_rh,ng_divu
      real(kind=dp_t), intent(inout) ::       rh(  -ng_rh:,  -ng_rh:,   -ng_rh:)
      real(kind=dp_t), intent(inout) :: divu_rhs(-ng_divu:,-ng_divu:, -ng_divu:)
      integer        , intent(inout) :: mm(0:,0:,0:)

      integer         :: i,j,k,nx,ny,nz

      nx = size(mm,dim=1) - 1
      ny = size(mm,dim=2) - 1
      nz = size(mm,dim=3) - 1

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

    end subroutine subtract_divu_from_rh_3d

end module nodal_divu_module
