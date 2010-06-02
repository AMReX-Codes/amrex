module stencil_nodal_module

  use bl_types
  use multifab_module
  use stencil_module
  use impose_neumann_bcs_module

  implicit none

  real (kind = dp_t), private, parameter :: ZERO  = 0.0_dp_t
  real (kind = dp_t), private, parameter :: ONE   = 1.0_dp_t
  real (kind = dp_t), private, parameter :: TWO   = 2.0_dp_t
  real (kind = dp_t), private, parameter :: THREE = 3.0_dp_t
  real (kind = dp_t), private, parameter :: FOUR  = 4.0_dp_t
  real (kind = dp_t), private, parameter :: FIVE  = 5.0_dp_t
  real (kind = dp_t), private, parameter :: SIX   = 6.0_dp_t
  real (kind = dp_t), private, parameter :: SEVEN = 7.0_dp_t
  real (kind = dp_t), private, parameter :: EIGHT = 8.0_dp_t
  real (kind = dp_t), private, parameter :: TEN   = 10.0_dp_t
  real (kind = dp_t), private, parameter :: HALF  = 0.5_dp_t
  real (kind = dp_t), private, parameter :: FOURTH= 0.25_dp_t
  real (kind = dp_t), private, parameter :: THIRD = 1.0_dp_t/3.0_dp_t
  real (kind = dp_t), private, parameter :: SIXTH = 1.0_dp_t/6.0_dp_t
  real (kind = dp_t), private, parameter :: FOUR_THIRD = 4.0_dp_t/3.0_dp_t

  private :: set_nodal_faces_edges_corners_3d

contains

  subroutine stencil_set_bc_nodal(sdim, bx, nbx, idx, mask, face_type, pd_periodic, bxa_periodic)
    integer,         intent(in   ) :: sdim
    type(box),       intent(in   ) :: bx, nbx
    type(imultifab), intent(inout) :: mask
    integer,         intent(in   ) :: idx
    integer,         intent(in   ) :: face_type(:,:,:)
    type(box),       intent(in   ) :: pd_periodic
    type(boxarray),  intent(in   ) :: bxa_periodic

    integer, pointer :: mp(:,:,:,:)
    type(box)        :: bx1
    type(boxarray)   :: ba
    integer          :: ii, dm, ib, jb, kb, jb_lo, kb_lo
    logical          :: nodal(sdim)

    nodal = .true.
    !
    ! Set the mask to BC_DIR or BC_NEU based on face_type at a physical boundary.
    !
    do dm = 1, bx%dim
       !
       ! Lo side
       !
       bx1 = nbx
       bx1%hi(dm) = bx1%lo(dm)
       mp => dataptr(mask%fbs(idx), bx1)
       if (face_type(idx,dm,1) == BC_NEU) mp = ibset(mp, BC_BIT(BC_NEU, dm, -1))
       if (face_type(idx,dm,1) == BC_DIR) mp = ibset(mp, BC_BIT(BC_DIR,  1,  0))
       !
       ! Hi side
       !
       bx1 = nbx
       bx1%lo(dm) = bx1%hi(dm)
       mp => dataptr(mask%fbs(idx), bx1)
       if (face_type(idx,dm,2) == BC_NEU) mp = ibset(mp, BC_BIT(BC_NEU, dm, +1))
       if (face_type(idx,dm,2) == BC_DIR) mp = ibset(mp, BC_BIT(BC_DIR,  1,  0))
    end do
    !
    ! Set the mask to BC_DIR at coarse-fine boundaries.
    !
    jb_lo = -1; if (dm .lt. 2) jb_lo = 1
    kb_lo = -1; if (dm .lt. 3) kb_lo = 1

    do kb = kb_lo, 1
       do jb = jb_lo, 1
          do ib = -1, 1
             bx1 = shift(bx,ib,1)
             if (dm > 1) bx1 = shift(bx1,jb,2)
             if (dm > 2) bx1 = shift(bx1,kb,3)
             bx1 = intersection(bx1, pd_periodic)
             if ( empty(bx1) ) cycle
             call boxarray_boxarray_diff(ba, bx1, bxa_periodic)
             do ii = 1, ba%nboxes
                bx1 = intersection(box_nodalize(ba%bxs(ii),nodal), nbx)
                if ( empty(bx1) ) cycle
                mp => dataptr(mask%fbs(idx), bx1)
                mp = ibset(mp, BC_BIT(BC_DIR,1,0))
             end do
             call destroy(ba)
          end do
       end do
    end do

  end subroutine stencil_set_bc_nodal

  subroutine s_simple_1d_nodal(ss, sg, mm, dh)
    real (kind = dp_t), intent(inout) :: ss(:,0:)
    real (kind = dp_t), intent(inout) :: sg(0:)
    integer           , intent(in   ) :: mm(:)
    real (kind = dp_t), intent(in   ) :: dh(:)

    integer i,nx
    real (kind = dp_t) f1

    f1 = ONE/dh(1)**2

    ! x derivatives
    nx = size(ss,dim=1)

    if (bc_neumann(mm( 1),1,-1)) sg( 0) = sg(   1)
    if (bc_neumann(mm(nx),1,+1)) sg(nx) = sg(nx-1)

    do i = 1,nx
       ss(i,1) = sg(i  )*f1
       ss(i,2) = sg(i-1)*f1
       ss(i,0) = -(sg(i)+sg(i-1))*f1
    end do

  end subroutine s_simple_1d_nodal

  subroutine s_cross_2d_nodal(ss, sg, mm, face_type, dh)
    real (kind = dp_t), intent(inout) :: ss(:,:,0:)
    real (kind = dp_t), intent(inout) :: sg(0:,0:)
    integer           , intent(in   ) :: mm(:,:)
    integer           , intent(in   ) :: face_type(:,:)
    real (kind = dp_t), intent(in   ) :: dh(:)

    integer            :: i, j, nx, ny

    nx = size(ss,dim=1)
    ny = size(ss,dim=2)

    ! Set sg on edges at a Neumann boundary.
    do i = 1,nx-1
       if (bc_neumann(mm(i, 1),2,-1)) sg(i, 0) = sg(i,1)
       if (bc_neumann(mm(i,ny),2,+1)) sg(i,ny) = sg(i,ny-1)
    end do

    do j = 1,ny-1
       if (bc_neumann(mm( 1,j),1,-1)) sg( 0,j) = sg(   1,j)
       if (bc_neumann(mm(nx,j),1,+1)) sg(nx,j) = sg(nx-1,j)
    end do

!   Note: we do the corners *after* each of the edge has been done.
    if (face_type(1,1) == BC_NEU) then
       sg(0, 0) = sg(1, 0)
       sg(0,ny) = sg(1,ny)
    end if
    if (face_type(1,2) == BC_NEU) then
       sg(nx, 0) = sg(nx-1,0)
       sg(nx,ny) = sg(nx-1,ny)
    end if
    if (face_type(2,1) == BC_NEU) then
       sg( 0,0) = sg( 0,1)
       sg(nx,0) = sg(nx,1)
    end if
    if (face_type(2,2) == BC_NEU) then
       sg( 0,ny) = sg( 0,ny-1)
       sg(nx,ny) = sg(nx,ny-1)
    end if

    do j = 1,ny
      do i = 1,nx
          ss(i,j,1) = sg(i  ,j-1) + sg(i  ,j  )
          ss(i,j,2) = sg(i-1,j-1) + sg(i-1,j  )
          ss(i,j,3) = sg(i-1,j  ) + sg(i  ,j  )
          ss(i,j,4) = sg(i-1,j-1) + sg(i  ,j-1)
          ss(i,j,0) = -(ss(i,j,1) + ss(i,j,2) + ss(i,j,3) + ss(i,j,4))
      end do
    end do

    ss = ss * (HALF / (dh(1))**2)

  end subroutine s_cross_2d_nodal

  subroutine s_simple_2d_one_sided(ss, sg, mm, face_type, dh)
    real (kind = dp_t), intent(inout) :: ss(:,:,0:)
    real (kind = dp_t), intent(inout) :: sg(0:,0:)
    integer           , intent(in   ) :: mm(:,:)
    integer           , intent(in   ) :: face_type(:,:)
    real (kind = dp_t), intent(in   ) :: dh(:)

    real (kind = dp_t), allocatable :: sg_int(:,:)

    integer            :: i, j, nx, ny

    nx = size(ss,dim=1)
    ny = size(ss,dim=2)

    allocate(sg_int(0:size(sg,dim=1)-1,0:size(sg,dim=2)-1))

    sg_int = ZERO
    do j = 1, ny-1
      do i = 1, nx-1
         sg_int(i,j) = sg(i,j)
      end do
    end do

    ! Set sg on edges at a Neumann boundary.
    do i = 1,nx-1
       if (bc_neumann(mm(i, 1),2,-1)) sg_int(i, 0) = sg_int(i,1)
       if (bc_neumann(mm(i,ny),2,+1)) sg_int(i,ny) = sg_int(i,ny-1)
    end do

    do j = 1,ny-1
       if (bc_neumann(mm( 1,j),1,-1)) sg_int( 0,j) = sg_int(   1,j)
       if (bc_neumann(mm(nx,j),1,+1)) sg_int(nx,j) = sg_int(nx-1,j)
    end do

!   Note: we do the corners *after* each of the edge has been done.
    if (face_type(1,1) == BC_NEU) then
       sg_int(0, 0) = sg_int(1, 0)
       sg_int(0,ny) = sg_int(1,ny)
    end if
    if (face_type(1,2) == BC_NEU) then
       sg_int(nx, 0) = sg_int(nx-1,0)
       sg_int(nx,ny) = sg_int(nx-1,ny)
    end if
    if (face_type(2,1) == BC_NEU) then
       sg_int( 0,0) = sg_int( 0,1)
       sg_int(nx,0) = sg_int(nx,1)
    end if
    if (face_type(2,2) == BC_NEU) then
       sg_int( 0,ny) = sg_int( 0,ny-1)
       sg_int(nx,ny) = sg_int(nx,ny-1)
    end if

    do j = 1,ny
      do i = 1,nx
          ss(i,j,1) = sg_int(i  ,j-1) + sg_int(i  ,j  )
          ss(i,j,2) = sg_int(i-1,j-1) + sg_int(i-1,j  )
          ss(i,j,3) = sg_int(i-1,j  ) + sg_int(i  ,j  )
          ss(i,j,4) = sg_int(i-1,j-1) + sg_int(i  ,j-1)
          ss(i,j,0) = -(ss(i,j,1) + ss(i,j,2) + ss(i,j,3) + ss(i,j,4))
      end do
    end do

    ss = ss * (HALF / (dh(1))**2)

    deallocate(sg_int)

  end subroutine s_simple_2d_one_sided

  subroutine s_dense_2d_nodal(ss, sg, mm, face_type, dh)
    real (kind = dp_t), intent(inout) :: ss(:,:,0:)
    real (kind = dp_t), intent(inout) :: sg(0:,0:)
    integer           , intent(in   ) :: mm(:,:)
    integer           , intent(in   ) :: face_type(:,:)
    real (kind = dp_t), intent(in   ) :: dh(:)

    integer            :: i, j, nx, ny

    nx = size(ss,dim=1)
    ny = size(ss,dim=2)

    ! Set sg on edges at a Neumann boundary.
    do i = 1,nx-1
       if (bc_neumann(mm(i, 1),2,-1)) sg(i, 0) = sg(i,1)
       if (bc_neumann(mm(i,ny),2,+1)) sg(i,ny) = sg(i,ny-1)
    end do

    do j = 1,ny-1
       if (bc_neumann(mm( 1,j),1,-1)) sg( 0,j) = sg(   1,j)
       if (bc_neumann(mm(nx,j),1,+1)) sg(nx,j) = sg(nx-1,j)
    end do

!   Note: we do the corners *after* each of the edge has been done.
    if (face_type(1,1) == BC_NEU) then
       sg(0, 0) = sg(1, 0)
       sg(0,ny) = sg(1,ny)
    end if
    if (face_type(1,2) == BC_NEU) then
       sg(nx, 0) = sg(nx-1,0)
       sg(nx,ny) = sg(nx-1,ny)
    end if
    if (face_type(2,1) == BC_NEU) then
       sg( 0,0) = sg( 0,1)
       sg(nx,0) = sg(nx,1)
    end if
    if (face_type(2,2) == BC_NEU) then
       sg( 0,ny) = sg( 0,ny-1)
       sg(nx,ny) = sg(nx,ny-1)
    end if

    do j = 1,ny
      do i = 1,nx

          ! Faces
          ss(i,j,2) = HALF * (sg(i  ,j-1) + sg(i-1,j-1))
          ss(i,j,4) = HALF * (sg(i-1,j  ) + sg(i-1,j-1))
          ss(i,j,5) = HALF * (sg(i  ,j  ) + sg(i  ,j-1))
          ss(i,j,7) = HALF * (sg(i  ,j  ) + sg(i-1,j  ))

          ! Corners
          ss(i,j,1) = sg(i-1,j-1)
          ss(i,j,3) = sg(i  ,j-1)
          ss(i,j,6) = sg(i-1,j  ) 
          ss(i,j,8) = sg(i,  j  )
 
          ss(i,j,0) = -(ss(i,j,1) + ss(i,j,2) + ss(i,j,3) + ss(i,j,4) &
                       +ss(i,j,5) + ss(i,j,6) + ss(i,j,7) + ss(i,j,8) ) 

      end do
    end do

    ss = ss * (THIRD / (dh(1))**2)

  end subroutine s_dense_2d_nodal

  subroutine set_nodal_faces_edges_corners_3d(nx, ny, nz, sg, mm)
    integer           , intent(in   ) :: nx, ny, nz
    real (kind = dp_t), intent(inout) :: sg(0:,0:,0:)
    integer           , intent(in   ) :: mm(:,:,:)

    integer :: i, j, k
    !
    ! Set sg on faces at a Neumann boundary.
    !
    !$OMP PARALLEL DO PRIVATE(i,j) IF(ny.ge.4)
    do j = 1,ny-1
       do i = 1,nx-1
          if (bc_neumann(mm(i,j, 1),3,-1)) sg(i,j, 0) = sg(i,j,1)
          if (bc_neumann(mm(i,j,nz),3,+1)) sg(i,j,nz) = sg(i,j,nz-1)
       end do
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i,k) IF(nz.ge.4)
    do k = 1,nz-1
       do i = 1,nx-1
          if (bc_neumann(mm(i, 1,k),2,-1)) sg(i, 0,k) = sg(i,1,k)
          if (bc_neumann(mm(i,ny,k),2,+1)) sg(i,ny,k) = sg(i,ny-1,k)
       end do
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(j,k) IF(nz.ge.4)
    do k = 1,nz-1
       do j = 1,ny-1
          if (bc_neumann(mm( 1,j,k),1,-1)) sg( 0,j,k) = sg(   1,j,k)
          if (bc_neumann(mm(nx,j,k),1,+1)) sg(nx,j,k) = sg(nx-1,j,k)
       end do
    end do
    !$OMP END PARALLEL DO
    !
    ! Set sg on edges at a Neumann boundary.
    !
    do i = 1,nx-1
       if (bc_neumann(mm(i, 1, 1),2,-1)) sg(i, 0, 0) = sg(i,1, 0) 
       if (bc_neumann(mm(i, 1,nz),2,-1)) sg(i, 0,nz) = sg(i,1,nz) 

       if (bc_neumann(mm(i,ny, 1),2,+1)) sg(i,ny, 0) = sg(i,ny-1, 0)
       if (bc_neumann(mm(i,ny,nz),2,+1)) sg(i,ny,nz) = sg(i,ny-1,nz)

       if (bc_neumann(mm(i, 1, 1),3,-1)) sg(i, 0, 0) = sg(i, 0,1)
       if (bc_neumann(mm(i,ny, 1),3,-1)) sg(i,ny, 0) = sg(i,ny,1)

       if (bc_neumann(mm(i, 1,nz),3,+1)) sg(i, 0,nz) = sg(i, 0,nz-1)
       if (bc_neumann(mm(i,ny,nz),3,+1)) sg(i,ny,nz) = sg(i,ny,nz-1)
    end do

    do j = 1,ny-1
       if (bc_neumann(mm( 1,j, 1),1,-1)) sg( 0,j, 0) = sg(1,j, 0)
       if (bc_neumann(mm( 1,j,nz),1,-1)) sg( 0,j,nz) = sg(1,j,nz)

       if (bc_neumann(mm(nx,j, 1),1,+1)) sg(nx,j, 0) = sg(nx-1,j, 0)
       if (bc_neumann(mm(nx,j,nz),1,+1)) sg(nx,j,nz) = sg(nx-1,j,nz)

       if (bc_neumann(mm( 1,j, 1),3,-1)) sg( 0,j, 0) = sg( 0,j,1)
       if (bc_neumann(mm(nx,j, 1),3,-1)) sg(nx,j, 0) = sg(nx,j,1)

       if (bc_neumann(mm( 1,j,nz),3,+1)) sg( 0,j,nz) = sg( 0,j,nz-1)
       if (bc_neumann(mm(nx,j,nz),3,+1)) sg(nx,j,nz) = sg(nx,j,nz-1)
    end do

    do k = 1,nz-1
       if (bc_neumann(mm( 1, 1,k),1,-1)) sg( 0, 0,k) = sg(1, 0,k)
       if (bc_neumann(mm( 1,ny,k),1,-1)) sg( 0,ny,k) = sg(1,ny,k)

       if (bc_neumann(mm(nx, 1,k),1,+1)) sg(nx, 0,k) = sg(nx-1, 0,k)
       if (bc_neumann(mm(nx,ny,k),1,+1)) sg(nx,ny,k) = sg(nx-1,ny,k)

       if (bc_neumann(mm( 1, 1,k),2,-1)) sg( 0, 0,k) = sg( 0,1,k)
       if (bc_neumann(mm(nx, 1,k),2,-1)) sg(nx, 0,k) = sg(nx,1,k)

       if (bc_neumann(mm( 1,ny,k),2,+1)) sg( 0,ny,k) = sg( 0,ny-1,k)
       if (bc_neumann(mm(nx,ny,k),2,+1)) sg(nx,ny,k) = sg(nx,ny-1,k)
    end do

    if (bc_neumann(mm( 1, 1, 1),1,-1)) sg( 0, 0, 0) = sg( 1, 0, 0) 
    if (bc_neumann(mm( 1, 1, 1),2,-1)) sg( 0, 0, 0) = sg( 0, 1, 0) 
    if (bc_neumann(mm( 1, 1, 1),3,-1)) sg( 0, 0, 0) = sg( 0, 0, 1) 

    if (bc_neumann(mm(nx, 1, 1),1,+1)) sg(nx, 0, 0) = sg(nx-1, 0, 0) 
    if (bc_neumann(mm(nx, 1, 1),2,-1)) sg(nx, 0, 0) = sg(nx  , 1, 0) 
    if (bc_neumann(mm(nx, 1, 1),3,-1)) sg(nx, 0, 0) = sg(nx  , 0, 1) 

    if (bc_neumann(mm( 1,ny, 1),1,-1)) sg( 0,ny, 0) = sg( 1,ny  , 0) 
    if (bc_neumann(mm( 1,ny, 1),2,+1)) sg( 0,ny, 0) = sg( 0,ny-1, 0) 
    if (bc_neumann(mm( 1,ny, 1),3,-1)) sg( 0,ny, 0) = sg( 0,ny  , 1) 

    if (bc_neumann(mm( 1, 1,nz),1,-1)) sg( 0, 0,nz) = sg( 1, 0,nz  ) 
    if (bc_neumann(mm( 1, 1,nz),2,-1)) sg( 0, 0,nz) = sg( 0, 1,nz  ) 
    if (bc_neumann(mm( 1, 1,nz),3,+1)) sg( 0, 0,nz) = sg( 0, 0,nz-1) 

    if (bc_neumann(mm(nx,ny, 1),1,+1)) sg(nx,ny, 0) = sg(nx-1,ny  , 0) 
    if (bc_neumann(mm(nx,ny, 1),2,+1)) sg(nx,ny, 0) = sg(nx  ,ny-1, 0) 
    if (bc_neumann(mm(nx,ny, 1),3,-1)) sg(nx,ny, 0) = sg(nx  ,ny  , 1) 

    if (bc_neumann(mm(nx, 1,nz),1,+1)) sg(nx, 0,nz) = sg(nx-1, 0,nz  ) 
    if (bc_neumann(mm(nx, 1,nz),2,-1)) sg(nx, 0,nz) = sg(nx  , 1,nz  ) 
    if (bc_neumann(mm(nx, 1,nz),3,+1)) sg(nx, 0,nz) = sg(nx  , 0,nz-1) 

    if (bc_neumann(mm( 1,ny,nz),1,-1)) sg( 0,ny,nz) = sg( 1,ny  ,nz  ) 
    if (bc_neumann(mm( 1,ny,nz),2,+1)) sg( 0,ny,nz) = sg( 0,ny-1,nz  ) 
    if (bc_neumann(mm( 1,ny,nz),3,+1)) sg( 0,ny,nz) = sg( 0,ny  ,nz-1) 

    if (bc_neumann(mm(nx,ny,nz),1,+1)) sg(nx,ny,nz) = sg(nx-1,ny  ,nz  ) 
    if (bc_neumann(mm(nx,ny,nz),2,+1)) sg(nx,ny,nz) = sg(nx  ,ny-1,nz  ) 
    if (bc_neumann(mm(nx,ny,nz),3,+1)) sg(nx,ny,nz) = sg(nx  ,ny  ,nz-1) 

  end subroutine set_nodal_faces_edges_corners_3d

  subroutine s_cross_3d_nodal(ss, sg, mm, dh)
    real (kind = dp_t), intent(inout) :: ss(:,:,:,0:)
    real (kind = dp_t), intent(inout) :: sg(0:,0:,0:)
    integer           , intent(inout) :: mm(:,:,:)
    real (kind = dp_t), intent(in   ) :: dh(:)

    integer :: i, j, k, nx, ny, nz

    nx = size(ss,dim=1)
    ny = size(ss,dim=2)
    nz = size(ss,dim=3)
    !
    !   BEGIN STENCIL
    !
    !   Stencil applies as follows :  (i is left to right, j is down to up)
    !
    !        at k-1        at k         at k+1
    !
    !                        3                 
    !          6          2  0  1          5   
    !                        4                  
    !
    !   END STENCIL
    !
    ss = ZERO

    call set_nodal_faces_edges_corners_3d(nx, ny, nz, sg, mm)

    !$OMP PARALLEL DO PRIVATE(i,j,k) IF(nz.ge.4)
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             !
             ! Faces in x-direction
             !
             ss(i,j,k,2) = (sg(i-1,j-1,k-1) + sg(i-1,j-1,k  ) &
                  +sg(i-1,j  ,k-1) + sg(i-1,j  ,k  ))
             ss(i,j,k,1) = (sg(i  ,j-1,k-1) + sg(i  ,j-1,k  ) &
                  +sg(i  ,j  ,k-1) + sg(i  ,j  ,k  ))
             !
             ! Faces in y-direction
             !
             ss(i,j,k,4) = (sg(i-1,j-1,k-1) + sg(i-1,j-1,k  ) &
                  +sg(i  ,j-1,k-1) + sg(i  ,j-1,k  ))
             ss(i,j,k,3) = (sg(i-1,j  ,k-1) + sg(i-1,j  ,k  ) &
                  +sg(i  ,j  ,k-1) + sg(i  ,j  ,k  ))
             !
             ! Faces in z-direction
             !
             ss(i,j,k,6) = (sg(i-1,j-1,k-1) + sg(i-1,j  ,k-1) &
                  +sg(i  ,j-1,k-1) + sg(i  ,j  ,k-1))
             ss(i,j,k,5) = (sg(i-1,j-1,k  ) + sg(i-1,j  ,k  ) &
                  +sg(i  ,j-1,k  ) + sg(i  ,j  ,k  ))

             ss(i,j,k,0) = -( ss(i,j,k,1) + ss(i,j,k,2) &
                  +ss(i,j,k,3) + ss(i,j,k,4) &
                  +ss(i,j,k,5) + ss(i,j,k,6) )
          end do
       end do
    end do
    !$OMP END PARALLEL DO

    ss = ss * (FOURTH / (dh(1))**2)

  end subroutine s_cross_3d_nodal

  subroutine s_simple_3d_one_sided(ss, sg, mm, dh)

    real (kind = dp_t), intent(inout) :: ss(:,:,:,0:)
    real (kind = dp_t), intent(inout) :: sg(0:,0:,0:)
    integer           , intent(inout) :: mm(:,:,:)
    real (kind = dp_t), intent(in   ) :: dh(:)

    integer :: i, j, k, nx, ny, nz
    real (kind = dp_t), allocatable :: sg_int(:,:,:)

    nx = size(ss,dim=1)
    ny = size(ss,dim=2)
    nz = size(ss,dim=3)
    !
    !   BEGIN STENCIL
    !
    !   Stencil applies as follows :  (i is left to right, j is down to up)
    !
    !        at k-1        at k         at k+1
    !
    !                        3                 
    !          6          2  0  1          5   
    !                        4                  
    !
    !   END STENCIL
    !
    ss = ZERO

    allocate(sg_int(0:size(sg,dim=1)-1,0:size(sg,dim=2)-1,0:size(sg,dim=3)-1))

    sg_int = ZERO

    do k = 1, nz-1
       do j = 1, ny-1
          do i = 1, nx-1
             sg_int(i,j,k) = sg(i,j,k)
          end do
       end do
    end do

    call set_nodal_faces_edges_corners_3d(nx, ny, nz, sg_int, mm)

    !$OMP PARALLEL DO PRIVATE(i,j,k) IF(nz.ge.4)
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             !
             ! Faces in x-direction
             !
             ss(i,j,k,2) = (sg_int(i-1,j-1,k-1) + sg_int(i-1,j-1,k  ) &
                  +sg_int(i-1,j  ,k-1) + sg_int(i-1,j  ,k  ))
             ss(i,j,k,1) = (sg_int(i  ,j-1,k-1) + sg_int(i  ,j-1,k  ) &
                  +sg_int(i  ,j  ,k-1) + sg_int(i  ,j  ,k  ))
             !
             ! Faces in y-direction
             !
             ss(i,j,k,4) = (sg_int(i-1,j-1,k-1) + sg_int(i-1,j-1,k  ) &
                  +sg_int(i  ,j-1,k-1) + sg_int(i  ,j-1,k  ))
             ss(i,j,k,3) = (sg_int(i-1,j  ,k-1) + sg_int(i-1,j  ,k  ) &
                  +sg_int(i  ,j  ,k-1) + sg_int(i  ,j  ,k  ))
             !
             ! Faces in z-direction
             !
             ss(i,j,k,6) = (sg_int(i-1,j-1,k-1) + sg_int(i-1,j  ,k-1) &
                  +sg_int(i  ,j-1,k-1) + sg_int(i  ,j  ,k-1))
             ss(i,j,k,5) = (sg_int(i-1,j-1,k  ) + sg_int(i-1,j  ,k  ) &
                  +sg_int(i  ,j-1,k  ) + sg_int(i  ,j  ,k  ))

             ss(i,j,k,0) = -( ss(i,j,k,1) + ss(i,j,k,2) &
                  +ss(i,j,k,3) + ss(i,j,k,4) &
                  +ss(i,j,k,5) + ss(i,j,k,6) )
          end do
       end do
    end do
    !$OMP END PARALLEL DO

    ss = ss * (FOURTH / (dh(1))**2)

    deallocate(sg_int)

  end subroutine s_simple_3d_one_sided

  subroutine s_dense_3d_nodal(ss, sg, mm, dh)
    real (kind = dp_t), intent(inout) :: ss(:,:,:,0:)
    real (kind = dp_t), intent(inout) :: sg(0:,0:,0:)
    integer           , intent(inout) :: mm(:,:,:)
    real (kind = dp_t), intent(in   ) :: dh(:)

    integer            :: i, j, k, nx, ny, nz
    real (kind = dp_t) :: fx,fy,fz,f0

    nx = size(ss,dim=1)
    ny = size(ss,dim=2)
    nz = size(ss,dim=3)
    !
    !   BEGIN STENCIL
    !
    !   Stencil applies as follows :  (i is left to right, j is down to up)
    !
    !        at k-1        at k         at k+1
    !
    !       6  7  8      11 24  12     18 19 20
    !       4 25  5      21  0  22     16 26 17
    !       1  2  3       9 23  10     13 14 15
    !
    !   END STENCIL
    !
    call set_nodal_faces_edges_corners_3d(nx, ny, nz, sg, mm)

    fx = ONE / (36.0_dp_t)
    fy = ONE / (36.0_dp_t)
    fz = ONE / (36.0_dp_t)
    f0 = FOUR * (fx + fy + fz)

    !$OMP PARALLEL DO PRIVATE(i,j,k) IF(nz.ge.4)
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             !
             ! Corners
             !
             ss(i,j,k, 1) = (fx+fy+fz)* sg(i-1,j-1,k-1)
             ss(i,j,k, 3) = (fx+fy+fz)* sg(i  ,j-1,k-1)
             ss(i,j,k, 6) = (fx+fy+fz)* sg(i-1,j  ,k-1)
             ss(i,j,k, 8) = (fx+fy+fz)* sg(i  ,j  ,k-1)
             ss(i,j,k,13) = (fx+fy+fz)* sg(i-1,j-1,k  )
             ss(i,j,k,15) = (fx+fy+fz)* sg(i  ,j-1,k  )
             ss(i,j,k,18) = (fx+fy+fz)* sg(i-1,j  ,k  )
             ss(i,j,k,20) = (fx+fy+fz)* sg(i  ,j  ,k  )
             !
             ! Edges in x-direction
             !
             ss(i,j,k, 2) = (TWO*fy+TWO*fz-fx)*(sg(i  ,j-1,k-1) + sg(i-1,j-1,k-1))
             ss(i,j,k, 7) = (TWO*fy+TWO*fz-fx)*(sg(i  ,j  ,k-1) + sg(i-1,j  ,k-1))
             ss(i,j,k,14) = (TWO*fy+TWO*fz-fx)*(sg(i  ,j-1,k  ) + sg(i-1,j-1,k  ))
             ss(i,j,k,19) = (TWO*fy+TWO*fz-fx)*(sg(i  ,j  ,k  ) + sg(i-1,j  ,k  ))
             !
             ! Edges in y-direction
             !
             ss(i,j,k, 4) = (TWO*fx+TWO*fz-fy)*(sg(i-1,j-1,k-1) + sg(i-1,j  ,k-1))
             ss(i,j,k, 5) = (TWO*fx+TWO*fz-fy)*(sg(i  ,j-1,k-1) + sg(i  ,j  ,k-1))
             ss(i,j,k,16) = (TWO*fx+TWO*fz-fy)*(sg(i-1,j-1,k  ) + sg(i-1,j  ,k  ))
             ss(i,j,k,17) = (TWO*fx+TWO*fz-fy)*(sg(i  ,j-1,k  ) + sg(i  ,j  ,k  ))
             !
             ! Edges in z-direction
             !
             ss(i,j,k, 9) = (TWO*fx+TWO*fy-fz)*(sg(i-1,j-1,k-1) + sg(i-1,j-1,k  ))
             ss(i,j,k,10) = (TWO*fx+TWO*fy-fz)*(sg(i  ,j-1,k-1) + sg(i  ,j-1,k  ))
             ss(i,j,k,11) = (TWO*fx+TWO*fy-fz)*(sg(i-1,j  ,k-1) + sg(i-1,j  ,k  ))
             ss(i,j,k,12) = (TWO*fx+TWO*fy-fz)*(sg(i  ,j  ,k-1) + sg(i  ,j  ,k  ))

             if (size(ss,dim=4) .eq. 27) then
                !
                ! Faces in x-direction (only non-zero for non-uniform dx)
                !
                ss(i,j,k,21) = (FOUR*fx-TWO*fy-TWO*fz)*(sg(i-1,j-1,k-1) + sg(i-1,j-1,k  ) &
                     +sg(i-1,j  ,k-1) + sg(i-1,j  ,k  ))
                ss(i,j,k,22) = (FOUR*fx-TWO*fy-TWO*fz)*(sg(i  ,j-1,k-1) + sg(i  ,j-1,k  ) &
                     +sg(i  ,j  ,k-1) + sg(i  ,j  ,k  ))
                !
                ! Faces in y-direction (only non-zero for non-uniform dx)
                !
                ss(i,j,k,23) = (FOUR*fy-TWO*fx-TWO*fz)*(sg(i-1,j-1,k-1) + sg(i-1,j-1,k  ) &
                     +sg(i  ,j-1,k-1) + sg(i  ,j-1,k  ))
                ss(i,j,k,24) = (FOUR*fy-TWO*fx-TWO*fz)*(sg(i-1,j  ,k-1) + sg(i-1,j  ,k  ) &
                     +sg(i  ,j  ,k-1) + sg(i  ,j  ,k  ))
                !
                ! Faces in z-direction (only non-zero for non-uniform dx)
                !
                ss(i,j,k,25) = (FOUR*fz-TWO*fx-TWO*fy)*(sg(i-1,j-1,k-1) + sg(i-1,j  ,k-1) &
                     +sg(i  ,j-1,k-1) + sg(i  ,j  ,k-1))
                ss(i,j,k,26) = (FOUR*fz-TWO*fx-TWO*fy)*(sg(i-1,j-1,k  ) + sg(i-1,j  ,k  ) &
                     +sg(i  ,j-1,k  ) + sg(i  ,j  ,k  ))
             end if

             ss(i,j,k,0) = -(sg(i-1,j-1,k-1) + sg(i,j-1,k-1) &
                  +sg(i-1,j  ,k-1) + sg(i,j  ,k-1) &
                  +sg(i-1,j-1,k  ) + sg(i,j-1,k  ) &
                  +sg(i-1,j  ,k  ) + sg(i,j  ,k  ) ) * f0

          end do
       end do
    end do
    !$OMP END PARALLEL DO

    ss = ss * (ONE / ((dh(1))**2))

  end subroutine s_dense_3d_nodal

  subroutine stencil_apply_1d_nodal(ss, dd, uu, mm, ng)
    integer, intent(in) :: ng
    real (kind = dp_t), intent(in   ) :: ss(:,0:)
    real (kind = dp_t), intent(inout) :: dd(0:)
    real (kind = dp_t), intent(inout) :: uu(1-ng:)
    integer           , intent(in   ) :: mm(:)
    integer i,lo(1)
 
    dd(:) = ZERO

    lo(:) = 1
    call impose_neumann_bcs_1d(uu,mm,lo,ng)
   
    i = 1
    if (.not. bc_dirichlet(mm(i),1,0)) then
      dd(i) = ss(i,0)*uu(i) + ss(i,1)*uu(i+1) + ss(i,2)*uu(i-1)
    end if

    do i = 2,size(ss,dim=1)-1
      dd(i) = ss(i,0)*uu(i) + ss(i,1) * uu(i+1) + ss(i,2) * uu(i-1)
    end do

    i = size(ss,dim=1)
    if (.not. bc_dirichlet(mm(i),1,0)) then
      dd(i) = ss(i,0)*uu(i) + ss(i,1)*uu(i+1) + ss(i,2)*uu(i-1)
    end if

  end subroutine stencil_apply_1d_nodal

  subroutine stencil_apply_2d_nodal(ss, dd, uu, mm, ng)
    integer, intent(in) :: ng
    real (kind = dp_t), intent(inout) :: uu(1-ng:,1-ng:)
    real (kind = dp_t), intent(inout) :: dd(0:,0:)
    real (kind = dp_t), intent(in   ) :: ss(:,:,0:)
    integer           , intent(in   ) :: mm(:,:)

    integer i,j,lo(2),nx,ny
    logical zeroit,iface,jface

    lo(:) = 1
    call impose_neumann_bcs_2d(uu,mm,lo,ng)
 
    nx = size(ss,dim=1)
    ny = size(ss,dim=2)

    if (size(ss,dim=3) .eq. 5) then

       do j = 1,ny
          jface = .false. ; if ( (j.eq.1).or.(j.eq.ny) ) jface = .true.
          do i = 1,nx
             iface = .false. ; if ( (i.eq.1).or.(i.eq.nx) ) iface = .true.

             zeroit = .false.

             if ( iface .or. jface ) then
                if (bc_dirichlet(mm(i,j),1,0)) zeroit = .true.
             end if

             if (zeroit) then
                dd(i,j) = ZERO
             else
                dd(i,j) = ss(i,j,0)*uu(i,j) + ss(i,j,1) * uu(i+1,j  ) &
                     + ss(i,j,2) * uu(i-1,j  ) &
                     + ss(i,j,3) * uu(i  ,j+1) &
                     + ss(i,j,4) * uu(i  ,j-1) 
             end if
          end do
       end do

    else if (size(ss,dim=3) .eq. 9) then

       do j = 1,ny
          jface = .false. ; if ( (j.eq.1).or.(j.eq.ny) ) jface = .true.
          do i = 1,nx
             iface = .false. ; if ( (i.eq.1).or.(i.eq.nx) ) iface = .true.

             zeroit = .false.

             if ( iface .or. jface ) then
                if (bc_dirichlet(mm(i,j),1,0)) zeroit = .true.
             end if

             if (zeroit) then
                dd(i,j) = ZERO
             else
                dd(i,j) = ss(i,j,0)*uu(i,j) + ss(i,j,1) * uu(i-1,j-1) &
                     + ss(i,j,2) * uu(i  ,j-1) &
                     + ss(i,j,3) * uu(i+1,j-1) &
                     + ss(i,j,4) * uu(i-1,j  ) &
                     + ss(i,j,5) * uu(i+1,j  ) &
                     + ss(i,j,6) * uu(i-1,j+1) &
                     + ss(i,j,7) * uu(i  ,j+1) &
                     + ss(i,j,8) * uu(i+1,j+1)
             end if
          end do
       end do

    end if

  end subroutine stencil_apply_2d_nodal

  subroutine stencil_apply_3d_nodal(ss, dd, uu, mm, ng, uniform_dh)
    integer           , intent(in   ) :: ng
    real (kind = dp_t), intent(inout) :: uu(1-ng:,1-ng:,1-ng:)
    real (kind = dp_t), intent(inout) :: dd(0:,0:,0:)
    real (kind = dp_t), intent(in   ) :: ss(:,:,:,0:)
    integer           , intent(in   ) :: mm(:,:,:)
    logical           , intent(in   ) :: uniform_dh

    integer i,j,k,lo(3),nx,ny,nz

    logical zeroit,jface,kface

    lo = 1

    call impose_neumann_bcs_3d(uu,mm,lo,ng)

    nz = size(ss,dim=3)
    ny = size(ss,dim=2)
    nx = size(ss,dim=1)

    if (size(ss,dim=4) .eq. 7) then

       !$OMP PARALLEL DO PRIVATE(i,j,k,zeroit,jface,kface) IF(nz.ge.4)
       do k = 1,nz
          kface = .false. ; if ( (k.eq.1) .or. (k.eq.nz) ) kface = .true.

          do j = 1,ny
             jface = .false. ; if ( (j.eq.1) .or. (j.eq.ny) ) jface = .true.

             do i = 1,nx

                zeroit = .false.

                if ( jface .or. kface .or. (i.eq.1) .or. (i.eq.nx) ) then
                   if (bc_dirichlet(mm(i,j,k),1,0)) zeroit = .true.
                end if

                if (zeroit) then
                   dd(i,j,k) = ZERO
                else
                   dd(i,j,k) = &
                        ss(i,j,k,0) * uu(i,j,k)       + &
                        ss(i,j,k,1) * uu(i+1,j  ,k  ) + &
                        ss(i,j,k,2) * uu(i-1,j  ,k  ) + &
                        ss(i,j,k,3) * uu(i  ,j+1,k  ) + &
                        ss(i,j,k,4) * uu(i  ,j-1,k  ) + &
                        ss(i,j,k,5) * uu(i  ,j  ,k+1) + &
                        ss(i,j,k,6) * uu(i  ,j  ,k-1)
                end if

             end do
          end do
       end do
       !$OMP END PARALLEL DO

    else if ((size(ss,dim=4) .eq. 21) .or. (size(ss,dim=4) .eq. 27)) then

       !$OMP PARALLEL DO PRIVATE(i,j,k,zeroit,jface,kface) IF(nz.ge.4)
       do k = 1,nz
          kface = .false. ; if ( (k.eq.1) .or. (k.eq.nz) ) kface = .true.

          do j = 1,ny
             jface = .false. ; if ( (j.eq.1) .or. (j.eq.ny) ) jface = .true.

             do i = 1,nx

                zeroit = .false.

                if ( jface .or. kface .or. (i.eq.1) .or. (i.eq.nx) ) then
                   if (bc_dirichlet(mm(i,j,k),1,0)) zeroit = .true.
                end if

                if (zeroit) then
                   dd(i,j,k) = ZERO
                else
                   dd(i,j,k) = ss(i,j,k,0)*uu(i,j,k) &
                        + ss(i,j,k, 1) * uu(i-1,j-1,k-1) + ss(i,j,k, 2) * uu(i  ,j-1,k-1) &
                        + ss(i,j,k, 3) * uu(i+1,j-1,k-1) + ss(i,j,k, 4) * uu(i-1,j  ,k-1) &
                        + ss(i,j,k, 5) * uu(i+1,j  ,k-1) + ss(i,j,k, 6) * uu(i-1,j+1,k-1) &
                        + ss(i,j,k, 7) * uu(i  ,j+1,k-1) + ss(i,j,k, 8) * uu(i+1,j+1,k-1) &
                        + ss(i,j,k, 9) * uu(i-1,j-1,k  ) + ss(i,j,k,10) * uu(i+1,j-1,k  ) &
                        + ss(i,j,k,11) * uu(i-1,j+1,k  ) + ss(i,j,k,12) * uu(i+1,j+1,k  ) &
                        + ss(i,j,k,13) * uu(i-1,j-1,k+1) + ss(i,j,k,14) * uu(i  ,j-1,k+1) &
                        + ss(i,j,k,15) * uu(i+1,j-1,k+1) + ss(i,j,k,16) * uu(i-1,j  ,k+1) &
                        + ss(i,j,k,17) * uu(i+1,j  ,k+1) + ss(i,j,k,18) * uu(i-1,j+1,k+1) &
                        + ss(i,j,k,19) * uu(i  ,j+1,k+1) + ss(i,j,k,20) * uu(i+1,j+1,k+1)

                   if ((size(ss,dim=4) .eq. 27) .and. (.not. uniform_dh)) then
                      !
                      ! Add faces (only non-zero for non-uniform dx)
                      !
                      dd(i,j,k) = dd(i,j,k) + &
                           ss(i,j,k,21) * uu(i-1,j  ,k  ) + ss(i,j,k,22) * uu(i+1,j  ,k  ) &
                           + ss(i,j,k,23) * uu(i  ,j-1,k  ) + ss(i,j,k,24) * uu(i  ,j+1,k  ) &
                           + ss(i,j,k,25) * uu(i  ,j  ,k-1) + ss(i,j,k,26) * uu(i  ,j  ,k+1)
                   end if
                end if

             end do
          end do
       end do
       !$OMP END PARALLEL DO
    else 
      print*,'BAD STENCIL SIZE IN APPLY_3D_NODAL ',size(ss,dim=4)
      call bl_error(' ')

    end if

  end subroutine stencil_apply_3d_nodal

end module stencil_nodal_module

