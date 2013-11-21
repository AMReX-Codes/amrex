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

  private :: set_faces_edges_corners_2d,set_faces_edges_corners_3d

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

  subroutine s_simple_1d_nodal(ss, sg, mm, dh)
    real (kind = dp_t), intent(inout) :: ss(0:,:)
    real (kind = dp_t), intent(inout) :: sg(0:)
    integer           , intent(in   ) :: mm(:)
    real (kind = dp_t), intent(in   ) :: dh(:)

    integer            :: i,nx
    real (kind = dp_t) :: f1

    f1 = ONE/dh(1)**2

    nx = size(ss,dim=2)

    if (bc_neumann(mm( 1),1,-1)) sg( 0) = sg(   1)
    if (bc_neumann(mm(nx),1,+1)) sg(nx) = sg(nx-1)

    do i = 1,nx
       ss(1,i) = sg(i  )*f1
       ss(2,i) = sg(i-1)*f1
       ss(0,i) = -(sg(i)+sg(i-1))*f1
    end do

  end subroutine s_simple_1d_nodal

  subroutine s_cross_2d_nodal(ss, sg, mm, face_type, dh)
    real (kind = dp_t), intent(inout) :: ss(0:,:,:)
    real (kind = dp_t), intent(inout) :: sg(0:,0:)
    integer           , intent(in   ) :: mm(:,:)
    integer           , intent(in   ) :: face_type(:,:)
    real (kind = dp_t), intent(in   ) :: dh(:)

    integer            :: i, j, nx, ny
    real (kind = dp_t) :: fac

    nx = size(ss,dim=2)
    ny = size(ss,dim=3)

    call set_faces_edges_corners_2d(nx, ny, sg, mm, face_type)

    fac = (HALF / (dh(1))**2)

    do j = 1,ny
      do i = 1,nx
          ss(1,i,j) = fac*(sg(i  ,j-1) + sg(i  ,j  ))
          ss(2,i,j) = fac*(sg(i-1,j-1) + sg(i-1,j  ))
          ss(3,i,j) = fac*(sg(i-1,j  ) + sg(i  ,j  ))
          ss(4,i,j) = fac*(sg(i-1,j-1) + sg(i  ,j-1))
          ss(0,i,j) = -(ss(1,i,j) + ss(2,i,j) + ss(3,i,j) + ss(4,i,j))
      end do
    end do

  end subroutine s_cross_2d_nodal

  subroutine s_simple_2d_one_sided(ss, sg, mm, face_type, dh)
    real (kind = dp_t), intent(inout) :: ss(0:,:,:)
    real (kind = dp_t), intent(inout) :: sg(0:,0:)
    integer           , intent(in   ) :: mm(:,:)
    integer           , intent(in   ) :: face_type(:,:)
    real (kind = dp_t), intent(in   ) :: dh(:)

    integer            :: i, j, nx, ny
    real (kind = dp_t) :: fac
    real (kind = dp_t) :: sg_int(0:size(sg,dim=1)-1,0:size(sg,dim=2)-1)

    nx = size(ss,dim=2)
    ny = size(ss,dim=3)

    sg_int = ZERO

    do j = 1, ny-1
      do i = 1, nx-1
         sg_int(i,j) = sg(i,j)
      end do
    end do

    call set_faces_edges_corners_2d(nx, ny, sg_int, mm, face_type)

    fac = (HALF / (dh(1))**2)

    do j = 1,ny
      do i = 1,nx
          ss(1,i,j) = fac*(sg_int(i  ,j-1) + sg_int(i  ,j  ))
          ss(2,i,j) = fac*(sg_int(i-1,j-1) + sg_int(i-1,j  ))
          ss(3,i,j) = fac*(sg_int(i-1,j  ) + sg_int(i  ,j  ))
          ss(4,i,j) = fac*(sg_int(i-1,j-1) + sg_int(i  ,j-1))
          ss(0,i,j) = -(ss(1,i,j) + ss(2,i,j) + ss(3,i,j) + ss(4,i,j))
      end do
    end do

  end subroutine s_simple_2d_one_sided

  subroutine s_dense_2d_nodal(ss, sg, ng_s, mm, face_type, dh, lo, hi)

    integer           , intent(in   ) :: lo(:),hi(:),ng_s 
    real (kind = dp_t), intent(inout) :: ss(0:,lo(1):,lo(2):)
    real (kind = dp_t), intent(inout) :: sg(lo(1)-ng_s:,lo(2)-ng_s:,:)
    integer           , intent(in   ) :: mm(lo(1):,lo(2):)
    integer           , intent(in   ) :: face_type(:,:)
    real (kind = dp_t), intent(in   ) :: dh(:)

    integer            :: i, j
    real (kind = dp_t) :: fac

    ! Set sg on edges at a Neumann boundary.
    do i = lo(1),hi(1)
       if (bc_neumann(mm(i,lo(2)  ),2,-1)) sg(i,lo(2)-1,:) = sg(i,lo(2),:) 
       if (bc_neumann(mm(i,hi(2)+1),2,+1)) sg(i,hi(2)+1,:) = sg(i,hi(2),:) 
    end do

    do j = lo(2),hi(2)
       if (bc_neumann(mm(lo(1)  ,j),1,-1)) sg(lo(1)-1,j,:) = sg(lo(1),j,:)
       if (bc_neumann(mm(hi(1)+1,j),1,+1)) sg(hi(1)+1,j,:) = sg(hi(1),j,:)
    end do

    ! Note: we do the corners *after* each of the edge has been done.
    ! Lo i
    if (face_type(1,1) == BC_NEU) then
       sg(lo(1)-1,lo(2)-1,:) = sg(lo(1),lo(2)-1,:)
       sg(lo(1)-1,hi(2)+1,:) = sg(lo(1),hi(2)+1,:)
    end if

    ! Hi i
    if (face_type(1,2) == BC_NEU) then
       sg(hi(1)+1,lo(2)-1,:) = sg(hi(1),lo(2)-1,:)
       sg(hi(1)+1,hi(2)+1,:) = sg(hi(1),hi(2)+1,:)
    end if

    ! Lo j
    if (face_type(2,1) == BC_NEU) then
       sg(lo(1)-1,lo(2)-1,:) = sg(lo(1)-1,lo(2),:)
       sg(hi(1)+1,lo(2)-1,:) = sg(hi(1)+1,lo(2),:)
    end if

    ! Hi j
    if (face_type(2,2) == BC_NEU) then
       sg(lo(1)-1,hi(2)+1,:) = sg(lo(1)-1,hi(2),:)
       sg(hi(1)+1,hi(2)+1,:) = sg(hi(1)+1,hi(2),:)
    end if

    fac = (THIRD / (dh(1))**2)

    do j = lo(2),hi(2)+1
      do i = lo(1),hi(1)+1

          ! Faces
          ss(2,i,j) = fac*(HALF * (sg(i  ,j-1,1) + sg(i-1,j-1,1)))
          ss(4,i,j) = fac*(HALF * (sg(i-1,j  ,1) + sg(i-1,j-1,1)))
          ss(5,i,j) = fac*(HALF * (sg(i  ,j  ,1) + sg(i  ,j-1,1)))
          ss(7,i,j) = fac*(HALF * (sg(i  ,j  ,1) + sg(i-1,j  ,1)))

          ! Corners
          ss(1,i,j) = fac*sg(i-1,j-1,1)
          ss(3,i,j) = fac*sg(i  ,j-1,1)
          ss(6,i,j) = fac*sg(i-1,j  ,1)
          ss(8,i,j) = fac*sg(i,  j  ,1)
 
          ss(0,i,j) = -(ss(1,i,j) + ss(2,i,j) + ss(3,i,j) + ss(4,i,j) &
                      + ss(5,i,j) + ss(6,i,j) + ss(7,i,j) + ss(8,i,j) )
      end do
    end do

  end subroutine s_dense_2d_nodal

  subroutine set_faces_edges_corners_2d(nx, ny, sg, mm, face_type)
    integer           , intent(in   ) :: nx, ny
    real (kind = dp_t), intent(inout) :: sg(0:,0:)
    integer           , intent(in   ) :: mm(:,:)
    integer           , intent(in   ) :: face_type(:,:)

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

  subroutine s_cross_3d_nodal(ss, sg, mm, dh)
    real (kind = dp_t), intent(inout) :: ss(0:,:,:,:)
    real (kind = dp_t), intent(inout) :: sg(0:,0:,0:)
    integer           , intent(inout) :: mm(:,:,:)
    real (kind = dp_t), intent(in   ) :: dh(:)

    integer            :: i, j, k, nx, ny, nz
    real (kind = dp_t) :: fac

    nx = size(ss,dim=2)
    ny = size(ss,dim=3)
    nz = size(ss,dim=4)
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
    call set_faces_edges_corners_3d(nx, ny, nz, sg, mm)

    fac = (FOURTH / (dh(1))**2)

    !$OMP PARALLEL DO PRIVATE(i,j,k) IF(nz.ge.4)
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             !
             ! Faces in x-direction
             !
             ss(2,i,j,k) = fac*((sg(i-1,j-1,k-1) + sg(i-1,j-1,k  ) &
                                +sg(i-1,j  ,k-1) + sg(i-1,j  ,k  )))
             ss(1,i,j,k) = fac*((sg(i  ,j-1,k-1) + sg(i  ,j-1,k  ) &
                                +sg(i  ,j  ,k-1) + sg(i  ,j  ,k  )))
             !
             ! Faces in y-direction
             !
             ss(4,i,j,k) = fac*((sg(i-1,j-1,k-1) + sg(i-1,j-1,k  ) &
                                +sg(i  ,j-1,k-1) + sg(i  ,j-1,k  )))
             ss(3,i,j,k) = fac*((sg(i-1,j  ,k-1) + sg(i-1,j  ,k  ) &
                                +sg(i  ,j  ,k-1) + sg(i  ,j  ,k  )))
             !
             ! Faces in z-direction
             !
             ss(6,i,j,k) = fac*((sg(i-1,j-1,k-1) + sg(i-1,j  ,k-1) &
                                +sg(i  ,j-1,k-1) + sg(i  ,j  ,k-1)))
             ss(5,i,j,k) = fac*((sg(i-1,j-1,k  ) + sg(i-1,j  ,k  ) &
                                +sg(i  ,j-1,k  ) + sg(i  ,j  ,k  )))

             ss(0,i,j,k) = -( ss(1,i,j,k) + ss(2,i,j,k) &
                  +ss(3,i,j,k) + ss(4,i,j,k) &
                  +ss(5,i,j,k) + ss(6,i,j,k) )
          end do
       end do
    end do
    !$OMP END PARALLEL DO

  end subroutine s_cross_3d_nodal

  subroutine s_simple_3d_one_sided(ss, sg, mm, dh)

    real (kind = dp_t), intent(inout) :: ss(0:,:,:,:)
    real (kind = dp_t), intent(inout) :: sg(0:,0:,0:)
    integer           , intent(inout) :: mm(:,:,:)
    real (kind = dp_t), intent(in   ) :: dh(:)

    integer :: i, j, k, nx, ny, nz
    real (kind = dp_t) :: fac
    real (kind = dp_t), allocatable :: sg_int(:,:,:)

    nx = size(ss,dim=2)
    ny = size(ss,dim=3)
    nz = size(ss,dim=4)
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
    allocate(sg_int(0:size(sg,dim=1)-1,0:size(sg,dim=2)-1,0:size(sg,dim=3)-1))

    sg_int = ZERO

    do k = 1, nz-1
       do j = 1, ny-1
          do i = 1, nx-1
             sg_int(i,j,k) = sg(i,j,k)
          end do
       end do
    end do

    call set_faces_edges_corners_3d(nx, ny, nz, sg_int, mm)

    fac = (FOURTH / (dh(1))**2)

    !$OMP PARALLEL DO PRIVATE(i,j,k) IF(nz.ge.4)
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             !
             ! Faces in x-direction
             !
             ss(2,i,j,k) = fac*((sg_int(i-1,j-1,k-1) + sg_int(i-1,j-1,k  ) &
                  +sg_int(i-1,j  ,k-1) + sg_int(i-1,j  ,k  )))
             ss(1,i,j,k) = fac*((sg_int(i  ,j-1,k-1) + sg_int(i  ,j-1,k  ) &
                  +sg_int(i  ,j  ,k-1) + sg_int(i  ,j  ,k  )))
             !
             ! Faces in y-direction
             !
             ss(4,i,j,k) = fac*((sg_int(i-1,j-1,k-1) + sg_int(i-1,j-1,k  ) &
                  +sg_int(i  ,j-1,k-1) + sg_int(i  ,j-1,k  )))
             ss(3,i,j,k) = fac*((sg_int(i-1,j  ,k-1) + sg_int(i-1,j  ,k  ) &
                  +sg_int(i  ,j  ,k-1) + sg_int(i  ,j  ,k  )))
             !
             ! Faces in z-direction
             !
             ss(6,i,j,k) = fac*((sg_int(i-1,j-1,k-1) + sg_int(i-1,j  ,k-1) &
                  +sg_int(i  ,j-1,k-1) + sg_int(i  ,j  ,k-1)))
             ss(5,i,j,k) = fac*((sg_int(i-1,j-1,k  ) + sg_int(i-1,j  ,k  ) &
                  +sg_int(i  ,j-1,k  ) + sg_int(i  ,j  ,k  )))

             ss(0,i,j,k) = -( ss(1,i,j,k) + ss(2,i,j,k) &
                  +ss(3,i,j,k) + ss(4,i,j,k) &
                  +ss(5,i,j,k) + ss(6,i,j,k) )
          end do
       end do
    end do
    !$OMP END PARALLEL DO

    deallocate(sg_int)

  end subroutine s_simple_3d_one_sided

  subroutine s_dense_3d_nodal(ss, sg, mm, dh)
    real (kind = dp_t), intent(inout) :: ss(0:,:,:,:)
    real (kind = dp_t), intent(inout) :: sg(0:,0:,0:)
    integer           , intent(inout) :: mm(:,:,:)
    real (kind = dp_t), intent(in   ) :: dh(:)

    integer            :: i, j, k, nx, ny, nz
    real (kind = dp_t) :: fx,fy,fz,f0,fac,fxyz,f2y2zx,f2x2zy,f2x2yz

    real (kind = dp_t), parameter :: ONETHIRTYSIXTH = ONE / (36.0_dp_t)

    nx = size(ss,dim=2)
    ny = size(ss,dim=3)
    nz = size(ss,dim=4)
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
    call set_faces_edges_corners_3d(nx, ny, nz, sg, mm)

    fx     = ONETHIRTYSIXTH
    fy     = ONETHIRTYSIXTH
    fz     = ONETHIRTYSIXTH
    f0     = FOUR * (fx + fy + fz)
    fac    = (ONE / ((dh(1))**2))
    fxyz   = (fx+fy+fz)
    f2y2zx = (TWO*fy+TWO*fz-fx)
    f2x2zy = (TWO*fx+TWO*fz-fy)
    f2x2yz = (TWO*fx+TWO*fy-fz)

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             !
             ! Corners
             !
             ss( 1,i,j,k) = fac*fxyz*sg(i-1,j-1,k-1)
             ss( 3,i,j,k) = fac*fxyz*sg(i  ,j-1,k-1)
             ss( 6,i,j,k) = fac*fxyz*sg(i-1,j  ,k-1)
             ss( 8,i,j,k) = fac*fxyz*sg(i  ,j  ,k-1)
             ss(13,i,j,k) = fac*fxyz*sg(i-1,j-1,k  )
             ss(15,i,j,k) = fac*fxyz*sg(i  ,j-1,k  )
             ss(18,i,j,k) = fac*fxyz*sg(i-1,j  ,k  )
             ss(20,i,j,k) = fac*fxyz*sg(i  ,j  ,k  )
             !
             ! Edges in x-direction
             !
             ss( 2,i,j,k) = fac*(f2y2zx*(sg(i  ,j-1,k-1) + sg(i-1,j-1,k-1)))
             ss( 7,i,j,k) = fac*(f2y2zx*(sg(i  ,j  ,k-1) + sg(i-1,j  ,k-1)))
             ss(14,i,j,k) = fac*(f2y2zx*(sg(i  ,j-1,k  ) + sg(i-1,j-1,k  )))
             ss(19,i,j,k) = fac*(f2y2zx*(sg(i  ,j  ,k  ) + sg(i-1,j  ,k  )))
             !
             ! Edges in y-direction
             !
             ss( 4,i,j,k) = fac*(f2x2zy*(sg(i-1,j-1,k-1) + sg(i-1,j  ,k-1)))
             ss( 5,i,j,k) = fac*(f2x2zy*(sg(i  ,j-1,k-1) + sg(i  ,j  ,k-1)))
             ss(16,i,j,k) = fac*(f2x2zy*(sg(i-1,j-1,k  ) + sg(i-1,j  ,k  )))
             ss(17,i,j,k) = fac*(f2x2zy*(sg(i  ,j-1,k  ) + sg(i  ,j  ,k  )))
             !
             ! Edges in z-direction
             !
             ss( 9,i,j,k) = fac*(f2x2yz*(sg(i-1,j-1,k-1) + sg(i-1,j-1,k  )))
             ss(10,i,j,k) = fac*(f2x2yz*(sg(i  ,j-1,k-1) + sg(i  ,j-1,k  )))
             ss(11,i,j,k) = fac*(f2x2yz*(sg(i-1,j  ,k-1) + sg(i-1,j  ,k  )))
             ss(12,i,j,k) = fac*(f2x2yz*(sg(i  ,j  ,k-1) + sg(i  ,j  ,k  )))

             if (size(ss,dim=1) .eq. 27) then
                !
                ! Faces in x-direction (only non-zero for non-uniform dx)
                !
                ss(21,i,j,k) = fac*((FOUR*fx-TWO*fy-TWO*fz)*(sg(i-1,j-1,k-1) + sg(i-1,j-1,k  ) &
                     +sg(i-1,j  ,k-1) + sg(i-1,j  ,k  )))
                ss(22,i,j,k) = fac*((FOUR*fx-TWO*fy-TWO*fz)*(sg(i  ,j-1,k-1) + sg(i  ,j-1,k  ) &
                     +sg(i  ,j  ,k-1) + sg(i  ,j  ,k  )))
                !
                ! Faces in y-direction (only non-zero for non-uniform dx)
                !
                ss(23,i,j,k) = fac*((FOUR*fy-TWO*fx-TWO*fz)*(sg(i-1,j-1,k-1) + sg(i-1,j-1,k  ) &
                     +sg(i  ,j-1,k-1) + sg(i  ,j-1,k  )))
                ss(24,i,j,k) = fac*((FOUR*fy-TWO*fx-TWO*fz)*(sg(i-1,j  ,k-1) + sg(i-1,j  ,k  ) &
                     +sg(i  ,j  ,k-1) + sg(i  ,j  ,k  )))
                !
                ! Faces in z-direction (only non-zero for non-uniform dx)
                !
                ss(25,i,j,k) = fac*((FOUR*fz-TWO*fx-TWO*fy)*(sg(i-1,j-1,k-1) + sg(i-1,j  ,k-1) &
                     +sg(i  ,j-1,k-1) + sg(i  ,j  ,k-1)))
                ss(26,i,j,k) = fac*((FOUR*fz-TWO*fx-TWO*fy)*(sg(i-1,j-1,k  ) + sg(i-1,j  ,k  ) &
                     +sg(i  ,j-1,k  ) + sg(i  ,j  ,k  )))
             end if

             ss(0,i,j,k) = -fac*f0*( sg(i-1,j-1,k-1) + sg(i,j-1,k-1) &
                                    +sg(i-1,j  ,k-1) + sg(i,j  ,k-1) &
                                    +sg(i-1,j-1,k  ) + sg(i,j-1,k  ) &
                                    +sg(i-1,j  ,k  ) + sg(i,j  ,k  ) )

          end do
       end do
    end do
    !$OMP END PARALLEL DO

  end subroutine s_dense_3d_nodal

end module nodal_stencil_module
