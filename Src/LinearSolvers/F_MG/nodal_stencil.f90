module nodal_stencil_module

  use bl_types
  use multifab_module
  use bc_functions_module
  use impose_neumann_bcs_module
  use stencil_types_module

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

  public :: set_faces_edges_corners_2d,set_faces_edges_corners_3d

contains

  subroutine stencil_set_bc_nodal(sdim, bx, nbx, idx, mask, face_type, pd_periodic, la_periodic)
    integer,         intent(in   ) :: sdim
    type(box),       intent(in   ) :: bx, nbx
    type(imultifab), intent(inout) :: mask
    integer,         intent(in   ) :: idx
    integer,         intent(in   ) :: face_type(:,:,:)
    type(box),       intent(in   ) :: pd_periodic
    type(layout),    intent(inout) :: la_periodic

    integer, pointer :: mp(:,:,:,:)
    type(box)        :: bx1
    type(boxarray)   :: ba
    integer          :: ii, dm, ib, jb, kb, jb_lo, kb_lo
    logical          :: nodal(sdim)

    nodal = .true.
    !
    ! Set the mask to BC_DIR or BC_NEU based on face_type at a physical boundary.
    !
    do dm = 1, sdim
       !
       ! Lo side
       !
       bx1 = nbx
       bx1%hi(dm) = bx1%lo(dm)
       mp => dataptr(mask, idx, bx1)
       if (face_type(idx,dm,1) == BC_NEU) then
          mp = ibset(mp, BC_BIT(BC_NEU, dm, -1))
       else if (face_type(idx,dm,1) == BC_DIR) then
          mp = ibset(mp, BC_BIT(BC_DIR, 1, 0))
       end if
       !
       ! Hi side
       !
       bx1 = nbx
       bx1%lo(dm) = bx1%hi(dm)
       mp => dataptr(mask, idx, bx1)
       if (face_type(idx,dm,2) == BC_NEU) then
          mp = ibset(mp, BC_BIT(BC_NEU, dm, +1))
       else if (face_type(idx,dm,2) == BC_DIR) then
          mp = ibset(mp, BC_BIT(BC_DIR, 1, 0))
       end if
    end do
    !
    ! Set the mask to BC_DIR at coarse-fine boundaries.
    !
    jb_lo = -1; if (sdim .lt. 2) jb_lo = 1
    kb_lo = -1; if (sdim .lt. 3) kb_lo = 1

    do kb = kb_lo, 1
       do jb = jb_lo, 1
          do ib = -1, 1
             bx1 = shift(bx,ib,1)
             if (sdim > 1) bx1 = shift(bx1,jb,2)
             if (sdim > 2) bx1 = shift(bx1,kb,3)
             bx1 = intersection(bx1, pd_periodic)
             if ( empty(bx1) ) cycle
             call layout_boxarray_diff(ba, bx1, la_periodic)
             do ii = 1, nboxes(ba)
                bx1 = intersection(box_nodalize(get_box(ba,ii),nodal), nbx)
                if ( empty(bx1) ) cycle
                mp => dataptr(mask, idx, bx1)
                mp = ibset(mp, BC_BIT(BC_DIR,1,0))
             end do
             call destroy(ba)
          end do
       end do
    end do

  end subroutine stencil_set_bc_nodal

  subroutine s_1d_nodal(ss, sg, mm, dh)
    real (kind = dp_t), intent(inout) :: ss(0:)
    real (kind = dp_t), intent(inout) :: sg(0:)
    integer           , intent(in   ) :: mm(:)
    real (kind = dp_t), intent(in   ) :: dh(:)

    integer            :: nx
    real (kind = dp_t) :: f1

    nx = size(sg,dim=1) - 2

    if (bc_neumann(mm( 1),1,-1)) sg( 0) = sg(   1)
    if (bc_neumann(mm(nx),1,+1)) sg(nx) = sg(nx-1)

    f1 = ONE/dh(1)**2

    ss(:) = f1 * sg(:)

  end subroutine s_1d_nodal

  subroutine s_2d_nodal(ss, sg, mm, dh)
    real (kind = dp_t), intent(inout) :: ss(0:,0:)
    real (kind = dp_t), intent(inout) :: sg(0:,0:)
    integer           , intent(in   ) :: mm(:,:)
    real (kind = dp_t), intent(in   ) :: dh(:)

    integer            :: nx, ny
    real (kind = dp_t) :: fac

    nx = size(sg,dim=1) - 2
    ny = size(sg,dim=2) - 2

    call set_faces_edges_corners_2d(nx, ny, sg, mm)

    ! This is not the full scaling for either stencil
    !   -- the cross stencil will need a factor of (1/2)
    !   -- the dense stencil will need a factor of (1/3)
    fac = 1.d0 / (dh(1))**2

    ss(:,:) = fac * sg(:,:)

  end subroutine s_2d_nodal

  subroutine s_3d_nodal(ss, sg, mm, dh)
    real (kind = dp_t), intent(inout) :: ss(0:,0:,0:)
    real (kind = dp_t), intent(inout) :: sg(0:,0:,0:)
    integer           , intent(inout) :: mm(:,:,:)
    real (kind = dp_t), intent(in   ) :: dh(:)

    integer            :: nx, ny, nz
    real (kind = dp_t) :: fac

    nx = size(sg,dim=1) - 2
    ny = size(sg,dim=2) - 2
    nz = size(sg,dim=3) - 2

    call set_faces_edges_corners_3d(nx, ny, nz, sg, mm)

    ! This is the right scaling for the cross stencil
    ! We end up needing to multiply by four later when we
    !    use this in the dense stencil.
    fac = (FOURTH / (dh(1))**2)

    ss(:,:,:) = fac * sg(:,:,:)

  end subroutine s_3d_nodal

  subroutine set_faces_edges_corners_2d(nx, ny, sg, mm)
    integer           , intent(in   ) :: nx, ny
    real (kind = dp_t), intent(inout) :: sg(0:,0:)
    integer           , intent(in   ) :: mm(:,:)

    integer :: i, j
    !
    ! Set sg on edges at a Neumann boundary.
    !
    do i = 1,nx-1
       if (bc_neumann(mm(i, 1),2,-1)) sg(i, 0) = sg(i,1)
       if (bc_neumann(mm(i,ny),2,+1)) sg(i,ny) = sg(i,ny-1)
    end do

    do j = 1,ny-1
       if (bc_neumann(mm( 1,j),1,-1)) sg( 0,j) = sg(   1,j)
       if (bc_neumann(mm(nx,j),1,+1)) sg(nx,j) = sg(nx-1,j)
    end do

    !
    ! Note: we do the corners *after* each of the edges has been done.
    !
    if (bc_neumann(mm( 1, 1),1,-1)) sg( 0, 0) = sg(1, 0)
    if (bc_neumann(mm( 1,ny),1,-1)) sg( 0,ny) = sg(1,ny)

    if (bc_neumann(mm(nx, 1),1,+1)) sg(nx, 0) = sg(nx-1, 0)
    if (bc_neumann(mm(nx,ny),1,+1)) sg(nx,ny) = sg(nx-1,ny)

    if (bc_neumann(mm( 1, 1),2,-1)) sg( 0, 0) = sg( 0,1)
    if (bc_neumann(mm(nx, 1),2,-1)) sg(nx, 0) = sg(nx,1)

    if (bc_neumann(mm( 1,ny),2,+1)) sg( 0,ny) = sg( 0,ny-1)
    if (bc_neumann(mm(nx,ny),2,+1)) sg(nx,ny) = sg(nx,ny-1)

  end subroutine set_faces_edges_corners_2d

  subroutine set_faces_edges_corners_3d(nx, ny, nz, sg, mm)
    integer           , intent(in   ) :: nx, ny, nz
    real (kind = dp_t), intent(inout) :: sg(0:,0:,0:)
    integer           , intent(in   ) :: mm(:,:,:)

    integer :: i, j, k
    !
    ! Set sg on faces at a Neumann boundary.
    !
    do j = 1,ny-1
       do i = 1,nx-1
          if (bc_neumann(mm(i,j, 1),3,-1)) sg(i,j, 0) = sg(i,j,1)
          if (bc_neumann(mm(i,j,nz),3,+1)) sg(i,j,nz) = sg(i,j,nz-1)
       end do
    end do

    do k = 1,nz-1
       do i = 1,nx-1
          if (bc_neumann(mm(i, 1,k),2,-1)) sg(i, 0,k) = sg(i,1,k)
          if (bc_neumann(mm(i,ny,k),2,+1)) sg(i,ny,k) = sg(i,ny-1,k)
       end do
    end do

    do k = 1,nz-1
       do j = 1,ny-1
          if (bc_neumann(mm( 1,j,k),1,-1)) sg( 0,j,k) = sg(   1,j,k)
          if (bc_neumann(mm(nx,j,k),1,+1)) sg(nx,j,k) = sg(nx-1,j,k)
       end do
    end do
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

  end subroutine set_faces_edges_corners_3d

end module nodal_stencil_module
