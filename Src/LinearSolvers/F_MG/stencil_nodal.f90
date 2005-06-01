module stencil_nodal_module

  use bl_types
  use multifab_module
  use stencil_module

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
  real (kind = dp_t), private, parameter :: FOUR_THIRD = 4.0_dp_t/3.0_dp_t

contains

  subroutine stencil_fill_nodal(ss, sg, dh_local, dh, mask, face_type)
    type(multifab ), intent(inout) :: ss
    type(multifab ), intent(inout) :: sg
    real(kind=dp_t), intent(in   ) :: dh_local(:), dh(:)
    type(imultifab), intent(inout) :: mask
    integer        , intent(in   ) :: face_type(:,:,:)

    real(kind=dp_t), pointer :: sp(:,:,:,:)
    real(kind=dp_t), pointer :: cp(:,:,:,:)
    integer        , pointer :: mp(:,:,:,:)
    integer                  :: i
    !
    ! Do this just to set everything in the mask to zero.
    !
    call setval(mask,BC_INT)

    do i = 1, ss%nboxes
       if ( multifab_remote(ss,i) ) cycle

       sp => dataptr(ss,   i)
       cp => dataptr(sg,   i)
       mp => dataptr(mask, i)

       call stencil_set_bc_nodal(ss, i, mask, face_type)

       select case (ss%dim)
       case (1)
          call s_simple_1d_nodal(sp(:,1,1,:), cp(:,1,1,1), mp(:,1,1,1), face_type(i,1,:), dh)
       case (2)
          call s_simple_2d_nodal(sp(:,:,1,:), cp(:,:,1,1), mp(:,:,1,1), face_type(i,:,:), dh)
       case (3)
          call s_simple_3d_nodal(sp(:,:,:,:), cp(:,:,:,1), mp(:,:,:,1), face_type(i,:,:), dh, dh_local)
       end select
    end do

!   call mask_pretty_print(mask, "mask", nodal = .true.)

  end subroutine stencil_fill_nodal

  subroutine stencil_set_bc_nodal(ss, idx, mask, face_type)
    type(multifab),  intent(in)    :: ss
    type(imultifab), intent(inout) :: mask
    integer,         intent(in)    :: idx
    integer,         intent(in)    :: face_type(:,:,:)

    integer, pointer :: mp(:,:,:,:)
    type(box)        :: bx, bx1, nbx
    type(boxarray)   :: ba, dba, pba
    integer          :: ii, dm, ib, jb, kb, jb_lo, kb_lo

    type(box)        :: bxs(3**ss%dim), sbx
    integer          :: shft(3**ss%dim,ss%dim), cnt, i, j

    bx  = get_box(ss,idx)
    nbx = get_ibox(ss, idx)
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

       if ((face_type(idx,dm,1) == BC_NEU) .or. (face_type(idx,dm,1) == BC_DIR)) then
          call boxarray_add_clean(dba, bx1)
       end if
       !
       ! Hi side
       !
       bx1 = nbx
       bx1%lo(dm) = bx1%hi(dm)
       mp => dataptr(mask%fbs(idx), bx1)
       if (face_type(idx,dm,2) == BC_NEU) mp = ibset(mp, BC_BIT(BC_NEU, dm, +1))
       if (face_type(idx,dm,2) == BC_DIR) mp = ibset(mp, BC_BIT(BC_DIR,  1,  0))

       if ((face_type(idx,dm,2) == BC_NEU) .or. (face_type(idx,dm,2) == BC_DIR)) then
          call boxarray_add_clean(dba, bx1)
       end if
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
             bx1 = box_intersection(bx1, ss%la%lap%pd)
             call boxarray_boxarray_diff(ba, bx1, ss%la%lap%bxa)
             do ii = 1, ba%nboxes
                bx1 = box_intersection(box_nodalize(ba%bxs(ii),ss%nodal), nbx)
                if (.not. empty(bx1)) then
                   mp => dataptr(mask%fbs(idx), bx1)
                   mp = ibset(mp, BC_BIT(BC_DIR,1,0))
                end if
             end do
             call destroy(ba)
          end do
       end do
    end do

    !
    ! Reset any Fine-Fine boundaries due to periodicity.
    !
    call multifab_internal_sync_shift(ss%la%lap%pd, nbx, ss%la%lap%pmask, ss%nodal, shft, cnt)
    !
    ! The "i == 1" component here is the box "nbx" itself which we ignore.
    !
    do i = 2, cnt
       sbx = shift(nbx, shft(i,:))
       do j = 1, ss%nboxes
          bx1 = intersection(get_ibox(ss,j), sbx)
          if ( .not. empty(bx1) ) then
             bx1 = shift(bx1, -shft(i,:))
             call boxarray_boxarray_diff(pba, bx1, dba)
             do ii = 1, pba%nboxes
                mp => dataptr(mask%fbs(idx), pba%bxs(ii))
                mp = BC_INT
             end do
             call destroy(pba)
          end if
       end do
    end do

    if ( built_q(dba) ) call destroy(dba)

  end subroutine stencil_set_bc_nodal

  subroutine s_simple_1d_nodal(ss, sg, mm, face_type, dh)
    real (kind = dp_t), intent(  out) :: ss(:,0:)
    real (kind = dp_t), intent(inout) :: sg(0:)
    integer           , intent(in   ) :: mm(:)
    real (kind = dp_t), intent(in   ) :: dh(:)
    integer           , intent(in   ) :: face_type(:)

    integer i,nx
    real (kind = dp_t) f1

    f1 = ONE/dh(1)**2

    ! x derivatives
    nx = size(ss,dim=1)

    do i = 1,nx
       ss(i,0) = (sg(i)+sg(i-1))*f1
       ss(i,1) = -sg(i  )*f1
       ss(i,2) = -sg(i-1)*f1
    end do

    if (bc_neumann(mm( 1),1,-1)) sg( 0) = sg(   1)
    if (bc_neumann(mm(nx),1,+1)) sg(nx) = sg(nx-1)

  end subroutine s_simple_1d_nodal

  subroutine s_simple_2d_nodal(ss, sg, mm, face_type, dh)
    real (kind = dp_t), intent(  out) :: ss(:,:,0:)
    real (kind = dp_t), intent(inout) :: sg(0:,0:)
    integer           , intent(in   ) :: mm(:,:)
    integer           , intent(in   ) :: face_type(:,:)
    real (kind = dp_t), intent(in   ) :: dh(:)

    integer            :: i, j, nx, ny
    real (kind = dp_t) :: fx, fy

!   fx = HALF*THIRD/dh(1)**2
!   fy = HALF*THIRD/dh(2)**2

    fx = HALF*THIRD
    fy = HALF*THIRD

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
          ss(i,j,1) = -(fx+fy)* sg(i-1,j-1)
          ss(i,j,2) = -(TWO*fy-fx)*(sg(i  ,j-1) + sg(i-1,j-1))
          ss(i,j,3) = -(fx+fy)* sg(i  ,j-1)

          ss(i,j,4) = -(TWO*fx-fy)*(sg(i-1,j  ) + sg(i-1,j-1))
          ss(i,j,5) = -(TWO*fx-fy)*(sg(i  ,j  ) + sg(i  ,j-1))

          ss(i,j,6) = -(fx+fy)* sg(i-1,j  ) 
          ss(i,j,7) = -(TWO*fy-fx)*(sg(i  ,j) + sg(i-1,j))
          ss(i,j,8) = -(fx+fy)    * sg(i,  j  )
 
          ss(i,j,0) = -ss(i,j,1) - ss(i,j,2) - ss(i,j,3) - ss(i,j,4) &
                      -ss(i,j,5) - ss(i,j,6) - ss(i,j,7) - ss(i,j,8)
      end do
    end do

  end subroutine s_simple_2d_nodal

  subroutine s_simple_3d_nodal(ss, sg, mm, face_type, dh, dh_local)
    real (kind = dp_t), intent(  out) :: ss(:,:,:,0:)
    real (kind = dp_t), intent(inout) :: sg(0:,0:,0:)
    integer           , intent(inout) :: mm(:,:,:)
    integer, intent(in)               :: face_type(:,:)
    real (kind = dp_t), intent(in   ) :: dh(:), dh_local(:)

    type(box     ) :: bx1
    type(boxarray) :: ba
    integer :: i, j, k, ilo, jlo, klo, ib, jb, ii, jj, nx, ny, nz
    real (kind = dp_t) :: fx,fy,fz,f0
    real (kind = dp_t) :: ratio

    nx = size(ss,dim=1)
    ny = size(ss,dim=2)
    nz = size(ss,dim=3)

    fx = ONE / (36.0_dp_t*dh(1)**2)
    fy = ONE / (36.0_dp_t*dh(2)**2)
    fz = ONE / (36.0_dp_t*dh(3)**2)
    f0 = FOUR * (fx + fy + fz)
!   f0 = THIRD * THIRD * ( ONE/dh(1)**2 + ONE/dh(2)**2 + ONE/dh(3)**2)

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

    ! Set sg on faces at a Neumann boundary.
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

    ! Set sg on edges at a Neumann boundary.
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

    do k = 1, nz
    do j = 1, ny
      do i = 1, nx

!         Corners
          ss(i,j,k, 1) = -(fx+fy+fz)* sg(i-1,j-1,k-1)
          ss(i,j,k, 3) = -(fx+fy+fz)* sg(i  ,j-1,k-1)
          ss(i,j,k, 6) = -(fx+fy+fz)* sg(i-1,j  ,k-1)
          ss(i,j,k, 8) = -(fx+fy+fz)* sg(i  ,j  ,k-1)
          ss(i,j,k,13) = -(fx+fy+fz)* sg(i-1,j-1,k  )
          ss(i,j,k,15) = -(fx+fy+fz)* sg(i  ,j-1,k  )
          ss(i,j,k,18) = -(fx+fy+fz)* sg(i-1,j  ,k  )
          ss(i,j,k,20) = -(fx+fy+fz)* sg(i  ,j  ,k  )

!         Edges in x-direction
          ss(i,j,k, 2) = -(TWO*fy+TWO*fz-fx)*(sg(i  ,j-1,k-1) + sg(i-1,j-1,k-1))
          ss(i,j,k, 7) = -(TWO*fy+TWO*fz-fx)*(sg(i  ,j  ,k-1) + sg(i-1,j  ,k-1))
          ss(i,j,k,14) = -(TWO*fy+TWO*fz-fx)*(sg(i  ,j-1,k  ) + sg(i-1,j-1,k  ))
          ss(i,j,k,19) = -(TWO*fy+TWO*fz-fx)*(sg(i  ,j  ,k  ) + sg(i-1,j  ,k  ))

!         Edges in y-direction
          ss(i,j,k, 4) = -(TWO*fx+TWO*fz-fy)*(sg(i-1,j-1,k-1) + sg(i-1,j  ,k-1))
          ss(i,j,k, 5) = -(TWO*fx+TWO*fz-fy)*(sg(i  ,j-1,k-1) + sg(i  ,j  ,k-1))
          ss(i,j,k,16) = -(TWO*fx+TWO*fz-fy)*(sg(i-1,j-1,k  ) + sg(i-1,j  ,k  ))
          ss(i,j,k,17) = -(TWO*fx+TWO*fz-fy)*(sg(i  ,j-1,k  ) + sg(i  ,j  ,k  ))

!         Edges in z-direction
          ss(i,j,k, 9) = -(TWO*fx+TWO*fy-fz)*(sg(i-1,j-1,k-1) + sg(i-1,j-1,k  ))
          ss(i,j,k,10) = -(TWO*fx+TWO*fy-fz)*(sg(i  ,j-1,k-1) + sg(i  ,j-1,k  ))
          ss(i,j,k,11) = -(TWO*fx+TWO*fy-fz)*(sg(i-1,j  ,k-1) + sg(i-1,j  ,k  ))
          ss(i,j,k,12) = -(TWO*fx+TWO*fy-fz)*(sg(i  ,j  ,k-1) + sg(i  ,j  ,k  ))

!         Faces in x-direction (only non-zero for non-uniform dx)
          ss(i,j,k,21) = -(FOUR*fx-TWO*fy-TWO*fz)*(sg(i-1,j-1,k-1) + sg(i-1,j-1,k  ) &
                                                  +sg(i-1,j  ,k-1) + sg(i-1,j  ,k  ))
          ss(i,j,k,22) = -(FOUR*fx-TWO*fy-TWO*fz)*(sg(i  ,j-1,k-1) + sg(i  ,j-1,k  ) &
                                                  +sg(i  ,j  ,k-1) + sg(i  ,j  ,k  ))

!         Faces in y-direction (only non-zero for non-uniform dx)
          ss(i,j,k,23) = -(FOUR*fy-TWO*fx-TWO*fz)*(sg(i-1,j-1,k-1) + sg(i-1,j-1,k  ) &
                                                  +sg(i  ,j-1,k-1) + sg(i  ,j-1,k  ))
          ss(i,j,k,24) = -(FOUR*fy-TWO*fx-TWO*fz)*(sg(i-1,j  ,k-1) + sg(i-1,j  ,k  ) &
                                                  +sg(i  ,j  ,k-1) + sg(i  ,j  ,k  ))

!         Faces in z-direction (only non-zero for non-uniform dx)
          ss(i,j,k,25) = -(FOUR*fz-TWO*fx-TWO*fy)*(sg(i-1,j-1,k-1) + sg(i-1,j  ,k-1) &
                                                  +sg(i  ,j-1,k-1) + sg(i  ,j  ,k-1))
          ss(i,j,k,26) = -(FOUR*fz-TWO*fx-TWO*fy)*(sg(i-1,j-1,k  ) + sg(i-1,j  ,k  ) &
                                                  +sg(i  ,j-1,k  ) + sg(i  ,j  ,k  ))

          ss(i,j,k,0) = (sg(i-1,j-1,k-1) + sg(i,j-1,k-1) &
                        +sg(i-1,j  ,k-1) + sg(i,j  ,k-1) &
                        +sg(i-1,j-1,k  ) + sg(i,j-1,k  ) &
                        +sg(i-1,j  ,k  ) + sg(i,j  ,k  ) ) * f0
                  
      end do
    end do
    end do

    ratio = dh_local(1) / dh(1)
    ss = ss*ratio

  end subroutine s_simple_3d_nodal

  subroutine stencil_apply_1d_nodal(ss, dd, uu, mm, ng)
    integer, intent(in) :: ng
    real (kind = dp_t), intent(in)  :: ss(:,0:)
    real (kind = dp_t), intent(out) :: dd(0:)
    real (kind = dp_t), intent(in)  :: uu(1-ng:)
    integer           , intent(in)  :: mm(:)
    integer i
   
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
    real (kind = dp_t), intent(  out) :: dd(0:,0:)
    real (kind = dp_t), intent(in   ) :: ss(:,:,0:)
    integer           , intent(in   ) :: mm(:,:)
    integer i,j,lo(2)

    lo(:) = 1
    call impose_neumann_bcs_2d(uu,mm,lo,ng)

    do j = 1,size(ss,dim=2)
    do i = 1,size(ss,dim=1)
       if (bc_dirichlet(mm(i,j),1,0)) then
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

  end subroutine stencil_apply_2d_nodal

  subroutine stencil_apply_3d_nodal(ss, dd, uu, mm, ng)
    integer, intent(in) :: ng
    real (kind = dp_t), intent(inout) :: uu(1-ng:,1-ng:,1-ng:)
    real (kind = dp_t), intent(  out) :: dd(0:,0:,0:)
    real (kind = dp_t), intent(in   ) :: ss(:,:,:,0:)
    integer           , intent(in   ) :: mm(:,:,:)
    integer i,j,k,lo(3)

    lo(:) = 1
    call impose_neumann_bcs_3d(uu,mm,lo,ng)

    do k = 1,size(ss,dim=3)
    do j = 1,size(ss,dim=2)
    do i = 1,size(ss,dim=1)

       if (bc_dirichlet(mm(i,j,k),1,0)) then
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

!       Add faces (only non-zero for non-uniform dx)
        dd(i,j,k) = dd(i,j,k) + &
            ss(i,j,k,21) * uu(i-1,j  ,k  ) + ss(i,j,k,22) * uu(i+1,j  ,k  ) &
          + ss(i,j,k,23) * uu(i  ,j-1,k  ) + ss(i,j,k,24) * uu(i  ,j+1,k  ) &
          + ss(i,j,k,25) * uu(i  ,j  ,k-1) + ss(i,j,k,26) * uu(i  ,j  ,k+1) 

       end if

    end do
    end do
    end do

  end subroutine stencil_apply_3d_nodal

  subroutine fine_edge_resid_1d(ss, dd, uu, res, mm, ratio, side)
    real (kind = dp_t), intent(in)  ::  ss(0:,0:)
    real (kind = dp_t), intent(out) ::  dd(:)
    real (kind = dp_t), intent(in)  ::  uu(-1:)
    real (kind = dp_t), intent(in)  :: res(-1:)
    integer           , intent(in)  ::  mm(:)
    integer, intent(in) :: ratio(:), side

    integer i,nx
    real (kind = dp_t) :: fac

    nx = size(ss,dim=1)-1

    fac = ONE / float(ratio(1))

    if (side == -1) then
      i = 0
      dd(1) = fac*ss(i,1)*(uu(i+1)-uu(i))
    else if (side == 1) then
      i = nx
      dd(1) = fac*ss(i,2)*(uu(i-1)-uu(i))
    end if

  end subroutine fine_edge_resid_1d

  subroutine fine_edge_resid_2d(dd, res, mm, ratio, side, lod)
    integer           , intent(in   ) :: lod(2)
    real (kind = dp_t), intent(  out) ::  dd(lod(1):,lod(2):)
    real (kind = dp_t), intent(inout) :: res(-1:,-1:)
    integer           , intent(in   ) ::  mm(0:,0:)
    integer, intent(in) :: ratio(:), side
    integer :: nx, ny, nxc, nyc
    integer :: hid(2)
    integer :: lo(2),ng_res
    integer :: i,j,ic,jc,m,n,isign,ioff,joff
    integer :: ileft,irght,jbot,jtop
    real (kind = dp_t) :: fac, fac0, fac1

    nx = size(mm,dim=1)-1
    ny = size(mm,dim=2)-1

    nxc = size(dd,dim=1)
    nyc = size(dd,dim=2)

    hid(1) = lod(1) + nxc-1
    hid(2) = lod(2) + nyc-1

    lo(:) = 0
    ng_res = 1
    call impose_neumann_bcs_2d(res,mm,lo,ng_res)

!   Lo/Hi i side
    if (side == -1 .or. side == 1) then

      if (side == -1) then
         i  = 0
         isign =  1
      else
         i  = nx
         isign = -1
      end if

      ic = lod(1)
      fac0 = ONE / ratio(2)

!     First average along the coarse-fine edge
      do jc = lod(2),hid(2)
         n = 0
         fac = HALF*ratio(2)*fac0
         j = (jc-lod(2))*ratio(2)
         if (j >  0) dd(ic,jc) = dd(ic,jc) + fac * res(i,j)
         if (j < ny) dd(ic,jc) = dd(ic,jc) + fac * res(i,j)

         do n = 1,ratio(2)-1
            fac = (ratio(2)-n)*fac0

            j = (jc-lod(2))*ratio(2) + n
            if (j < ny) then
               if (jc==lod(2).and..not.bc_dirichlet(mm(i,j),1,0)) &
                 fac = HALF * fac
               dd(ic,jc) = dd(ic,jc) + fac * res(i,j)
            end if

            j = (jc-lod(2))*ratio(2) - n
            if (j > 0) then
               if (jc==hid(2).and..not.bc_dirichlet(mm(i,j),1,0)) &
                 fac = HALF * fac
               dd(ic,jc) = dd(ic,jc) + fac * res(i,j)
            end if
         end do

      end do

      j = 0
      if (bc_neumann(mm(i,j),2,-1)) dd(ic,lod(2)) = TWO*dd(ic,lod(2))

      j = (hid(2)-lod(2))*ratio(2)
      if (bc_neumann(mm(i,j),2, 1)) dd(ic,hid(2)) = TWO*dd(ic,hid(2))

!     Now average towards the interior of the fine grid
      fac0 = fac0 / ratio(1)
      do n = 0,ratio(2)-1
         fac1 = (ratio(2)-n) * fac0
         if (n == 0) fac1 = HALF * fac1
         do m = 1,ratio(1)-1
            fac = (ratio(1)-m) * fac1
            ioff = i+isign*m
            do jc = lod(2),hid(2)
               j = (jc-lod(2))*ratio(2)
               jbot = j-n
               jtop = j+n
               if (j==0  .and. bc_neumann(mm(i,j),2,-1)) jbot = jtop 
               if (j==ny .and. bc_neumann(mm(i,j),2,+1)) jtop = jbot 
     
               if (j==0 .and. .not. bc_neumann(mm(i,j),2,-1)) then

                  if (n==0 .and. .not.bc_dirichlet(mm(ioff,j),1,0)) then

                    dd(ic,jc) = dd(ic,jc) + fac * res(ioff,j)

                  else if (n==0 .and. bc_dirichlet(mm(ioff,j),1,0)) then
                  else
                    dd(ic,jc) = dd(ic,jc) + HALF * fac * res(ioff,jtop)
                  end if
               else if (j==ny .and. .not. bc_neumann(mm(i,j),2,+1)) then
                  if (n==0 .and. .not.bc_dirichlet(mm(ioff,j),1,0)) then

                    dd(ic,jc) = dd(ic,jc) + fac * res(ioff,j)

                  else if (n==0 .and. bc_dirichlet(mm(ioff,j),1,0)) then
                  else
                    dd(ic,jc) = dd(ic,jc) + HALF * fac * res(ioff,jbot)
                  end if
               else
                  dd(ic,jc) = dd(ic,jc) + fac * ( res(ioff,jtop) + &
                                                  res(ioff,jbot) )
               end if 
            end do
         end do
      end do

      do jc = lod(2),hid(2)
         if (.not.bc_dirichlet(mm(i,(jc-lod(2))*ratio(2)),1,0)) dd(ic,jc) = ZERO
      end do

!   Lo/Hi j side
    else if (side == -2 .or. side == 2) then

      if (side == -2) then
         j  = 0
         isign =  1
      else
         j  = ny
         isign = -1
      end if

      jc = lod(2)
      fac0 = ONE / ratio(1) 

!     First average along the coarse-fine edge
      do ic = lod(1),hid(1)
         do n = 0,ratio(1)-1
            fac = (ratio(1)-n)*fac0
            if (n == 0) fac = HALF * fac

            i = (ic-lod(1))*ratio(1) + n
            if (i == 0) then

               dd(ic,jc) = dd(ic,jc) + fac * res(i,j)

            else if (i < nx) then

               if (ic==lod(1).and.n>0.and..not.bc_dirichlet(mm(i,j),1,0)) &
                 fac = HALF * fac

               dd(ic,jc) = dd(ic,jc) + fac * res(i,j)

            end if

            i = (ic-lod(1))*ratio(1) - n
            if (i == nx) then

              dd(ic,jc) = dd(ic,jc) + fac * res(i,j)

            else if (i > 0) then

               if (ic==hid(1).and.n>0.and..not.bc_dirichlet(mm(i,j),1,0)) &
                 fac = HALF * fac

               dd(ic,jc) = dd(ic,jc) + fac * res(i,j)

            end if
         end do
      end do

      i = 0
      if (bc_neumann(mm(i,j),1,-1)) dd(lod(1),jc) = TWO*dd(lod(1),jc)

      i = (hid(1)-lod(1))*ratio(1)
      if (bc_neumann(mm(i,j),1, 1)) dd(hid(1),jc) = TWO*dd(hid(1),jc)

!     Now average towards the interior of the fine grid
      fac0 = fac0 / ratio(2)
      do n = 0,ratio(1)-1
         fac1 = (ratio(1)-n) * fac0
         if (n == 0) fac1 = HALF * fac1
         do m = 1,ratio(2)-1
            joff = j + isign*m
            fac = (ratio(2)-m) * fac1
            do ic = lod(1),hid(1)
               i = (ic-lod(1))*ratio(1)
               ileft = i-n
               irght = i+n
               if (i==0  .and. bc_neumann(mm(i,j),1,-1)) ileft = irght 
               if (i==nx .and. bc_neumann(mm(i,j),1,+1)) irght = ileft 

               if (i==0 .and. .not. bc_neumann(mm(i,j),1,-1)) then
                  if (n==0 .and. .not.bc_dirichlet(mm(i,joff),1,0)) then

                    dd(ic,jc) = dd(ic,jc) + fac * res(i,joff)

                  else if (n==0 .and. bc_dirichlet(mm(i,joff),1,0)) then
                  else
                    dd(ic,jc) = dd(ic,jc) + HALF * fac * res(irght,joff)
                  end if
               else if (i==nx .and. .not. bc_neumann(mm(i,j),1,+1)) then
                  if (n==0 .and. .not.bc_dirichlet(mm(i,joff),1,0)) then

                    dd(ic,jc) = dd(ic,jc) + fac * res(i,joff)

                  else if (n==0 .and. bc_dirichlet(mm(i,joff),1,0)) then
                  else
                    dd(ic,jc) = dd(ic,jc) + HALF * fac * res(ileft,joff)
                  end if
               else
                  dd(ic,jc) = dd(ic,jc) + fac * ( res(irght,joff) + &
                                                  res(ileft,joff) )
               end if
            end do
         end do
      end do

      do ic = lod(1),hid(1)
         if (.not.bc_dirichlet(mm((ic-lod(1))*ratio(1),j),1,0)) dd(ic,jc) = ZERO
      end do

    end if

  end subroutine fine_edge_resid_2d

  subroutine fine_edge_resid_3d(dd, res, mm, ratio, side, lod)
    integer, intent(in) :: lod(:)
    real (kind = dp_t), intent(  out) ::  dd(lod(1):,lod(2):,lod(3):)
    real (kind = dp_t), intent(inout) :: res(-1:,-1:,-1:)
    integer           , intent(in   ) ::  mm(0:,0:,0:)
    integer, intent(in) :: ratio(:),side
    integer :: nx, ny, nz, nxc, nyc, nzc
    integer :: hid(3),lo(3),ng_res
    integer :: i,j,k,l,ic,jc,kc,m,n
    integer :: isign,ioff,joff,koff
    integer :: ileft,irght,jbot,jtop,kdwn,kup
    real (kind = dp_t) :: fac, fac0, fac1, fac2
    real (kind = dp_t) :: corner_fac

    nx = size(mm,dim=1)-1
    ny = size(mm,dim=2)-1
    nz = size(mm,dim=3)-1

    nxc = size(dd,dim=1)
    nyc = size(dd,dim=2)
    nzc = size(dd,dim=3)

    hid(1) = lod(1) + nxc-1
    hid(2) = lod(2) + nyc-1
    hid(3) = lod(3) + nzc-1

    lo(:) = 0
    ng_res = 1
    call impose_neumann_bcs_3d(res,mm,lo,ng_res)

    if (side == -1 .or. side == 1) then

      if (side == -1) then
         i  = 0
         isign =  1
      else
         i  = nx
         isign = -1
      end if

      ic = lod(1)
      fac0 = 1.0_dp_t / (ratio(2)*ratio(3))

!     First average along the coarse-fine face.
      do kc = lod(3),hid(3)
      do jc = lod(2),hid(2)
         do n = 0,ratio(2)-1
           fac2 = (ratio(2)-n)*fac0
           if (n == 0) fac2 = HALF * fac2

           j = (jc-lod(2))*ratio(2) + n
           if (j < ny) then
            do l = 0,ratio(3)-1
               fac = (ratio(3)-l)*fac2
               if (l == 0) fac = HALF * fac

               k = (kc-lod(3))*ratio(3) + l
               if (k < nz) then
                  dd(ic,jc,kc) = dd(ic,jc,kc) + fac * res(i,j,k)
               end if

               k = (kc-lod(3))*ratio(3) - l
               if (k >  0) then
                  dd(ic,jc,kc) = dd(ic,jc,kc) + fac * res(i,j,k)
               end if
            end do
           end if

           j = (jc-lod(2))*ratio(2) - n
           if (j > 0) then
            do l = 0,ratio(3)-1
               fac = (ratio(3)-l)*fac2
               if (l == 0) fac = HALF * fac

               k = (kc-lod(3))*ratio(3) + l
               if (k < nz) then
                  dd(ic,jc,kc) = dd(ic,jc,kc) + fac * res(i,j,k)
               end if

               k = (kc-lod(3))*ratio(3) - l
               if (k >  0) then
                  dd(ic,jc,kc) = dd(ic,jc,kc) + fac * res(i,j,k)
               end if
            end do
           end if

         end do

      end do
      end do

      jc = lod(2)
      kc = lod(3)
      j = 0
      k = 0
      if (.not. bc_neumann(mm(i,j,k),2,-1) .and. &
          .not. bc_neumann(mm(i,j,k),3,-1) ) &
         dd(ic,jc,kc) = dd(ic,jc,kc) + 0.25_dp_t * res(i,j,k) / 3.0_dp_t

      jc = hid(2)
      kc = lod(3)
      j = ny
      k = 0
      if (.not. bc_neumann(mm(i,j,k),2,+1) .and. &
          .not. bc_neumann(mm(i,j,k),3,-1) ) &
         dd(ic,jc,kc) = dd(ic,jc,kc) + 0.25_dp_t * res(i,j,k) / 3.0_dp_t

      jc = lod(2)
      kc = hid(3)
      j = 0
      k = nz
      if (.not. bc_neumann(mm(i,j,k),2,-1) .and. &
          .not. bc_neumann(mm(i,j,k),3,+1) ) &
         dd(ic,jc,kc) = dd(ic,jc,kc) + 0.25_dp_t * res(i,j,k) / 3.0_dp_t

      jc = hid(2)
      kc = hid(3)
      j = ny
      k = nz
      if (.not. bc_neumann(mm(i,j,k),2,+1) .and. &
          .not. bc_neumann(mm(i,j,k),3,+1) ) &
         dd(ic,jc,kc) = dd(ic,jc,kc) + 0.25_dp_t * res(i,j,k) / 3.0_dp_t

      j = 0
      do kc = lod(3),hid(3)
         k = (kc-lod(3))*ratio(3)
         if (bc_neumann(mm(i,j,k),2,-1)) dd(ic,lod(2),kc) = TWO*dd(ic,lod(2),kc)
      end do

      j = (hid(2)-lod(2))*ratio(2)
      do kc = lod(3),hid(3)
         k = (kc-lod(3))*ratio(3)
         if (bc_neumann(mm(i,j,k),2, 1)) dd(ic,hid(2),kc) = TWO*dd(ic,hid(2),kc)
      end do

      k = 0
      do jc = lod(2),hid(2)
         j = (jc-lod(2))*ratio(2)
         if (bc_neumann(mm(i,j,k),3,-1)) dd(ic,jc,lod(3)) = TWO*dd(ic,jc,lod(3))
      end do

      k = (hid(3)-lod(3))*ratio(3)
      do jc = lod(2),hid(2)
         j = (jc-lod(2))*ratio(2)
         if (bc_neumann(mm(i,j,k),3, 1)) dd(ic,jc,hid(3)) = TWO*dd(ic,jc,hid(3))
      end do

!     Now average towards the interior of the grid.
      fac0 = fac0 / ratio(1)
      ic = lod(1)
      do l = 0, ratio(3)-1
        fac2 = (ratio(3)-l) * fac0
        if (l == 0) fac2 = HALF * fac2
        do n = 0, ratio(2)-1
          fac1 = (ratio(2)-n) * fac2
          if (n == 0) fac1 = HALF * fac1
          do m = 1, ratio(1)-1
            ioff = i+isign*m
            fac = (ratio(1)-m) * fac1
            if (m == 0) fac = HALF * fac
            do kc = lod(3),hid(3)
              k = (kc-lod(3))*ratio(3)
              do jc = lod(2),hid(2)
                j = (jc-lod(2))*ratio(2)
                jtop = j+n
                jbot = j-n
                kup  = k+l
                kdwn = k-l
                if (j==0  .and. bc_neumann(mm(i,j,k),2,-1)) jbot = jtop
                if (j==ny .and. bc_neumann(mm(i,j,k),2,+1)) jtop = jbot
                if (k==0  .and. bc_neumann(mm(i,j,k),3,-1)) kdwn = kup 
                if (k==nz .and. bc_neumann(mm(i,j,k),3,+1)) kup  = kdwn

                if ( ( (jc==lod(2) .and. .not. bc_neumann(mm(i,j,k),2,-1)) .or. &
                       (jc==hid(2) .and. .not. bc_neumann(mm(i,j,k),2,+1)) ) .and. &
                     ( (kc==lod(3) .and. .not. bc_neumann(mm(i,j,k),3,-1)) .or. &
                       (kc==hid(3) .and. .not. bc_neumann(mm(i,j,k),3,+1)) ) ) then
                   corner_fac = 1.0_dp_t / 3.0_dp_t
                else if ( ( (jc==lod(2) .and. .not. bc_neumann(mm(i,j,k),2,-1)) .or. &
                            (jc==hid(2) .and. .not. bc_neumann(mm(i,j,k),2,+1)) ) .and. &
                          .not. &
                          ( (kc==lod(3) .and. .not. bc_neumann(mm(i,j,k),3,-1)) .or. &
                            (kc==hid(3) .and. .not. bc_neumann(mm(i,j,k),3,+1)) ) ) then
                   corner_fac = 1.0_dp_t  / 2.0_dp_t
                else if ( .not. &
                          ( (jc==lod(2) .and. .not. bc_neumann(mm(i,j,k),2,-1)) .or. &
                            (jc==hid(2) .and. .not. bc_neumann(mm(i,j,k),2,+1)) ) .and. &
                          ( (kc==lod(3) .and. .not. bc_neumann(mm(i,j,k),3,-1)) .or. &
                            (kc==hid(3) .and. .not. bc_neumann(mm(i,j,k),3,+1)) ) ) then
                   corner_fac = 1.0_dp_t / 2.0_dp_t
                else
                   corner_fac = 1.0_dp_t
                end if

                if ( (j-n >  0 .or. bc_neumann(mm(i,j,k),2,-1)) .and. &
                     (j-n < ny .or. bc_neumann(mm(i,j,k),2,+1)) ) then
                   if ( (k-l >  0 .or. bc_neumann(mm(i,j,k),3,-1)) .and. &
                        (k-l < nz .or. bc_neumann(mm(i,j,k),3,+1)) ) then
                      dd(ic,jc,kc) = dd(ic,jc,kc) + corner_fac * &
                           fac*res(ioff,jbot,kdwn) 
                   end if
                   if ( (k+l >  0 .or. bc_neumann(mm(i,j,k),3,-1)) .and. &
                        (k+l < nz .or. bc_neumann(mm(i,j,k),3,+1)) ) then
                      dd(ic,jc,kc) = dd(ic,jc,kc) + corner_fac * &
                           fac*res(ioff,jbot,kup) 
                   end if
                end if

                if ( (j+n >  0 .or. bc_neumann(mm(i,j,k),2,-1)) .and. &
                     (j+n < ny .or. bc_neumann(mm(i,j,k),2,+1)) ) then
                   if ( (k-l >  0 .or. bc_neumann(mm(i,j,k),3,-1)) .and. &
                        (k-l < nz .or. bc_neumann(mm(i,j,k),3,+1)) ) then
                      dd(ic,jc,kc) = dd(ic,jc,kc) + corner_fac * &
                           fac*res(ioff,jtop,kdwn) 
                   end if
                   if ( (k+l >  0 .or. bc_neumann(mm(i,j,k),3,-1)) .and. &
                        (k+l < nz .or. bc_neumann(mm(i,j,k),3,+1)) ) then
                      dd(ic,jc,kc) = dd(ic,jc,kc) + corner_fac * &
                           fac*res(ioff,jtop,kup) 
                   end if
                end if

              end do
            end do
          end do
        end do
      end do

      do kc = lod(3),hid(3)
      do jc = lod(2),hid(2)
         j = (jc-lod(2))*ratio(2)
         k = (kc-lod(3))*ratio(3)
         if (.not.bc_dirichlet(mm(i,j,k),1,0)) dd(ic,jc,kc) = ZERO
      end do
      end do

    else if (side == -2 .or. side == 2) then

      if (side == -2) then
         j  = 0
         isign =  1
      else
         j  = ny
         isign = -1
      end if
      jc = lod(2)
      fac0 = 1.0_dp_t / (ratio(1)*ratio(3))

!     First average along the coarse-fine face.
      do kc = lod(3),hid(3)
      do ic = lod(1),hid(1)
         do n = 0,ratio(1)-1
           fac2 = (ratio(1)-n)*fac0
           if (n == 0) fac2 = HALF * fac2

           i = (ic-lod(1))*ratio(1) + n
           if (i < nx) then
            do l = 0,ratio(3)-1
               fac = (ratio(3)-l)*fac2
               if (l == 0) fac = HALF * fac

               k = (kc-lod(3))*ratio(3) + l
               if (k < nz) then
                  dd(ic,jc,kc) = dd(ic,jc,kc) + fac * res(i,j,k)
               end if

               k = (kc-lod(3))*ratio(3) - l
               if (k >  0) then
                  dd(ic,jc,kc) = dd(ic,jc,kc) + fac * res(i,j,k)
               end if
            end do
           end if

           i = (ic-lod(1))*ratio(1) - n
           if (i > 0) then
            do l = 0,ratio(3)-1
               fac = (ratio(3)-l)*fac2
               if (l == 0) fac = HALF * fac

               k = (kc-lod(3))*ratio(3) + l
               if (k < nz) then
                  dd(ic,jc,kc) = dd(ic,jc,kc) + fac * res(i,j,k)
               end if

               k = (kc-lod(3))*ratio(3) - l
               if (k >  0) then
                  dd(ic,jc,kc) = dd(ic,jc,kc) + fac * res(i,j,k)
               end if
            end do
           end if

         end do

      end do
      end do

      ic = lod(1)
      kc = lod(3)
      i = 0
      k = 0
      if (.not. bc_neumann(mm(i,j,k),1,-1) .and. &
          .not. bc_neumann(mm(i,j,k),3,-1) ) &
         dd(ic,jc,kc) = dd(ic,jc,kc) + 0.25_dp_t * res(i,j,k) / 3.0_dp_t

      ic = hid(1)
      kc = lod(3)
      i = nx
      k = 0
      if (.not. bc_neumann(mm(i,j,k),1,+1) .and. &
          .not. bc_neumann(mm(i,j,k),3,-1) ) &
         dd(ic,jc,kc) = dd(ic,jc,kc) + 0.25_dp_t * res(i,j,k) / 3.0_dp_t

      ic = lod(1)
      kc = hid(3)
      i = 0
      k = nz
      if (.not. bc_neumann(mm(i,j,k),1,-1) .and. &
          .not. bc_neumann(mm(i,j,k),3,+1) ) &
         dd(ic,jc,kc) = dd(ic,jc,kc) + 0.25_dp_t * res(i,j,k) / 3.0_dp_t

      ic = hid(1)
      kc = hid(3)
      i = nx
      k = nz
      if (.not. bc_neumann(mm(i,j,k),1,+1) .and. &
          .not. bc_neumann(mm(i,j,k),3,+1) ) &
         dd(ic,jc,kc) = dd(ic,jc,kc) + 0.25_dp_t * res(i,j,k) / 3.0_dp_t

      i = 0
      do kc = lod(3),hid(3)
         k = (kc-lod(3))*ratio(3)
         if (bc_neumann(mm(i,j,k),1,-1)) dd(lod(1),jc,kc) = TWO*dd(lod(1),jc,kc)
      end do

      i = (hid(1)-lod(1))*ratio(1)
      do kc = lod(3),hid(3)
         k = (kc-lod(3))*ratio(3)
         if (bc_neumann(mm(i,j,k),1, 1)) dd(hid(1),jc,kc) = TWO*dd(hid(1),jc,kc)
      end do

      k = 0
      do ic = lod(1),hid(1)
         i = (ic-lod(1))*ratio(1)
         if (bc_neumann(mm(i,j,k),3,-1)) dd(ic,jc,lod(3)) = TWO*dd(ic,jc,lod(3))
      end do

      k = (hid(3)-lod(3))*ratio(3)
      do ic = lod(1),hid(1)
         i = (ic-lod(1))*ratio(1)
         if (bc_neumann(mm(i,j,k),3, 1)) dd(ic,jc,hid(3)) = TWO*dd(ic,jc,hid(3))
      end do

!     Now average towards the interior of the grid.
      fac0 = fac0 / ratio(2)
      jc = lod(2)
      do l = 0, ratio(3)-1
        fac2 = (ratio(3)-l) * fac0
        if (l == 0) fac2 = HALF * fac2
        do n = 0, ratio(1)-1
          fac1 = (ratio(1)-n) * fac2
          if (n == 0) fac1 = HALF * fac1
          do m = 1, ratio(2)-1
            joff = j+isign*m
            fac = (ratio(2)-m) * fac1
            if (m == 0) fac = HALF * fac
            do kc = lod(3),hid(3)
              k = (kc-lod(3))*ratio(3)
              do ic = lod(1),hid(1)
                i = (ic-lod(1))*ratio(1)
                irght = i+n
                ileft = i-n
                kup  = k+l
                kdwn = k-l
                if (i==0  .and. bc_neumann(mm(i,j,k),1,-1)) ileft = irght
                if (i==nx .and. bc_neumann(mm(i,j,k),1,+1)) irght = ileft
                if (k==0  .and. bc_neumann(mm(i,j,k),3,-1)) kdwn = kup 
                if (k==nz .and. bc_neumann(mm(i,j,k),3,+1)) kup  = kdwn

                if ( ( (ic==lod(1) .and. .not. bc_neumann(mm(i,j,k),1,-1)) .or. &
                       (ic==hid(1) .and. .not. bc_neumann(mm(i,j,k),1,+1)) ) .and. &
                     ( (kc==lod(3) .and. .not. bc_neumann(mm(i,j,k),3,-1)) .or. &
                       (kc==hid(3) .and. .not. bc_neumann(mm(i,j,k),3,+1)) ) ) then
                   corner_fac = 1.0_dp_t / 3.0_dp_t
                else if ( ( (ic==lod(1) .and. .not. bc_neumann(mm(i,j,k),1,-1)) .or. &
                            (ic==hid(1) .and. .not. bc_neumann(mm(i,j,k),1,+1)) ) .and. &
                          .not. &
                          ( (kc==lod(3) .and. .not. bc_neumann(mm(i,j,k),3,-1)) .or. &
                            (kc==hid(3) .and. .not. bc_neumann(mm(i,j,k),3,+1)) ) ) then
                   corner_fac = 1.0_dp_t  / 2.0_dp_t
                else if ( .not. &
                          ( (ic==lod(1) .and. .not. bc_neumann(mm(i,j,k),1,-1)) .or. &
                            (ic==hid(1) .and. .not. bc_neumann(mm(i,j,k),1,+1)) ) .and. &
                          ( (kc==lod(3) .and. .not. bc_neumann(mm(i,j,k),3,-1)) .or. &
                            (kc==hid(3) .and. .not. bc_neumann(mm(i,j,k),3,+1)) ) ) then
                   corner_fac = 1.0_dp_t / 2.0_dp_t
                else
                   corner_fac = 1.0_dp_t
                end if

                if ( (i-n >  0 .or. bc_neumann(mm(i,j,k),1,-1)) .and. &
                     (i-n < nx .or. bc_neumann(mm(i,j,k),1,+1)) ) then
                   if ( (k-l >  0 .or. bc_neumann(mm(i,j,k),3,-1)) .and. &
                        (k-l < nz .or. bc_neumann(mm(i,j,k),3,+1)) ) then
                      dd(ic,jc,kc) = dd(ic,jc,kc) + corner_fac * &
                           fac*res(ileft,joff,kdwn) 
                   end if
                   if ( (k+l >  0 .or. bc_neumann(mm(i,j,k),3,-1)) .and. &
                        (k+l < nz .or. bc_neumann(mm(i,j,k),3,+1)) ) then
                      dd(ic,jc,kc) = dd(ic,jc,kc) + corner_fac * &
                           fac*res(ileft,joff,kup) 
                   end if
                end if

                if ( (i+n >  0 .or. bc_neumann(mm(i,j,k),1,-1)) .and. &
                     (i+n < nx .or. bc_neumann(mm(i,j,k),1,+1)) ) then
                   if ( (k-l >  0 .or. bc_neumann(mm(i,j,k),3,-1)) .and. &
                        (k-l < nz .or. bc_neumann(mm(i,j,k),3,+1)) ) then
                      dd(ic,jc,kc) = dd(ic,jc,kc) + corner_fac * &
                           fac*res(irght,joff,kdwn) 
                   end if
                   if ( (k+l >  0 .or. bc_neumann(mm(i,j,k),3,-1)) .and. &
                        (k+l < nz .or. bc_neumann(mm(i,j,k),3,+1)) ) then
                      dd(ic,jc,kc) = dd(ic,jc,kc) + corner_fac * &
                           fac*res(irght,joff,kup) 
                   end if
                end if

              end do
            end do
          end do
        end do
      end do

      do kc = lod(3),hid(3)
      do ic = lod(1),hid(1)
         i = (ic-lod(1))*ratio(1)
         k = (kc-lod(3))*ratio(3)
         if (.not.bc_dirichlet(mm(i,j,k),1,0)) dd(ic,jc,kc) = ZERO
      end do
      end do

    else 

      if (side == -3) then
         k  = 0
         isign =  1
      else
         k  = nz
         isign = -1
      end if
      kc = lod(3)
      fac0 = 1.0_dp_t / (ratio(1)*ratio(2))

!     First average along the coarse-fine face.
      do jc = lod(2),hid(2)
      do ic = lod(1),hid(1)
         do n = 0,ratio(1)-1
           fac2 = (ratio(1)-n)*fac0
           if (n == 0) fac2 = HALF * fac2

           i = (ic-lod(1))*ratio(1) + n
           if (i < nx) then
            do l = 0,ratio(2)-1
               fac = (ratio(2)-l)*fac2
               if (l == 0) fac = HALF * fac

               j = (jc-lod(2))*ratio(2) + l
               if (j < ny) then
                  dd(ic,jc,kc) = dd(ic,jc,kc) + fac * res(i,j,k)
               end if

               j = (jc-lod(2))*ratio(2) - l
               if (j >  0) then
                  dd(ic,jc,kc) = dd(ic,jc,kc) + fac * res(i,j,k)
               end if
            end do
           end if

           i = (ic-lod(1))*ratio(1) - n
           if (i > 0) then
            do l = 0,ratio(2)-1
               fac = (ratio(2)-l)*fac2
               if (l == 0) fac = HALF * fac

               j = (jc-lod(2))*ratio(2) + l
               if (j < ny) then
                  dd(ic,jc,kc) = dd(ic,jc,kc) + fac * res(i,j,k)
               end if

               j = (jc-lod(2))*ratio(2) - l
               if (j >  0) then
                  dd(ic,jc,kc) = dd(ic,jc,kc) + fac * res(i,j,k)
               end if
            end do
           end if

         end do

      end do
      end do

      ic = lod(1)
      jc = lod(2)
      i = 0
      j = 0
      if (.not. bc_neumann(mm(i,j,k),1,-1) .and. &
          .not. bc_neumann(mm(i,j,k),2,-1) ) &
         dd(ic,jc,kc) = dd(ic,jc,kc) + 0.25_dp_t * res(i,j,k) / 3.0_dp_t

      ic = hid(1)
      jc = lod(2)
      i = nx
      j = 0
      if (.not. bc_neumann(mm(i,j,k),1,+1) .and. &
          .not. bc_neumann(mm(i,j,k),2,-1) ) &
         dd(ic,jc,kc) = dd(ic,jc,kc) + 0.25_dp_t * res(i,j,k) / 3.0_dp_t

      ic = lod(1)
      jc = hid(2)
      i = 0
      j = ny
      if (.not. bc_neumann(mm(i,j,k),1,-1) .and. &
          .not. bc_neumann(mm(i,j,k),2,+1) ) &
         dd(ic,jc,kc) = dd(ic,jc,kc) + 0.25_dp_t * res(i,j,k) / 3.0_dp_t

      ic = hid(1)
      jc = hid(2)
      i = nx
      j = ny
      if (.not. bc_neumann(mm(i,j,k),1,+1) .and. &
          .not. bc_neumann(mm(i,j,k),2,+1) ) &
         dd(ic,jc,kc) = dd(ic,jc,kc) + 0.25_dp_t * res(i,j,k) / 3.0_dp_t

      i = 0
      do jc = lod(2),hid(2)
         j = (jc-lod(2))*ratio(2)
         if (bc_neumann(mm(i,j,k),1,-1)) dd(lod(1),jc,kc) = TWO*dd(lod(1),jc,kc)
      end do

      i = (hid(1)-lod(1))*ratio(1)
      do jc = lod(2),hid(2)
         j = (jc-lod(2))*ratio(2)
         if (bc_neumann(mm(i,j,k),1,+1)) dd(hid(1),jc,kc) = TWO*dd(hid(1),jc,kc)
      end do

      j = 0
      do ic = lod(1),hid(1)
         i = (ic-lod(1))*ratio(1)
         if (bc_neumann(mm(i,j,k),2,-1)) dd(ic,lod(2),kc) = TWO*dd(ic,lod(2),kc)
      end do

      j = (hid(2)-lod(2))*ratio(2)
      do ic = lod(1),hid(1)
         i = (ic-lod(1))*ratio(1)
         if (bc_neumann(mm(i,j,k),2,+1)) dd(ic,hid(2),kc) = TWO*dd(ic,hid(2),kc)
      end do

!     Now average towards the interior of the grid.
      fac0 = fac0 / ratio(3)
      kc = lod(3)
      do l = 0, ratio(2)-1
        fac2 = (ratio(2)-l) * fac0
        if (l == 0) fac2 = HALF * fac2
        do n = 0, ratio(1)-1
          fac1 = (ratio(1)-n) * fac2
          if (n == 0) fac1 = HALF * fac1
          do m = 1, ratio(3)-1
            koff = k+isign*m
            fac = (ratio(3)-m) * fac1
            if (m == 0) fac = HALF * fac
            do jc = lod(2),hid(2)
              j = (jc-lod(2))*ratio(2)
              do ic = lod(1),hid(1)
                i = (ic-lod(1))*ratio(1)
                irght = i+n
                ileft = i-n
                jtop  = j+l
                jbot  = j-l
                if (i==0  .and. bc_neumann(mm(i,j,k),1,-1)) ileft = irght
                if (i==nx .and. bc_neumann(mm(i,j,k),1,+1)) irght = ileft
                if (j==0  .and. bc_neumann(mm(i,j,k),2,-1)) jbot  = jtop
                if (j==ny .and. bc_neumann(mm(i,j,k),2,+1)) jtop  = jbot
                if ( ( (ic==lod(1) .and. .not. bc_neumann(mm(i,j,k),1,-1)) .or. &
                       (ic==hid(1) .and. .not. bc_neumann(mm(i,j,k),1,+1)) ) .and. &
                     ( (jc==lod(2) .and. .not. bc_neumann(mm(i,j,k),2,-1)) .or. &
                       (jc==hid(2) .and. .not. bc_neumann(mm(i,j,k),2,+1)) ) ) then
                   corner_fac = 1.0_dp_t / 3.0_dp_t
                else if ( ( (ic==lod(1) .and. .not. bc_neumann(mm(i,j,k),1,-1)) .or. &
                            (ic==hid(1) .and. .not. bc_neumann(mm(i,j,k),1,+1)) ) .and. &
                          .not. &
                          ( (jc==lod(2) .and. .not. bc_neumann(mm(i,j,k),2,-1)) .or. &
                            (jc==hid(2) .and. .not. bc_neumann(mm(i,j,k),2,+1)) ) ) then
                   corner_fac = 1.0_dp_t  / 2.0_dp_t
                else if ( .not. &
                          ( (ic==lod(1) .and. .not. bc_neumann(mm(i,j,k),1,-1)) .or. &
                            (ic==hid(1) .and. .not. bc_neumann(mm(i,j,k),1,+1)) ) .and. &
                          ( (jc==lod(2) .and. .not. bc_neumann(mm(i,j,k),2,-1)) .or. &
                            (jc==hid(2) .and. .not. bc_neumann(mm(i,j,k),2,+1)) ) ) then
                   corner_fac = 1.0_dp_t / 2.0_dp_t
                else
                   corner_fac = 1.0_dp_t
                end if

                if ( (i-n >  0 .or. bc_neumann(mm(i,j,k),1,-1)) .and. &
                     (i-n < nx .or. bc_neumann(mm(i,j,k),1,+1)) ) then
                   if ( (j-l >  0 .or. bc_neumann(mm(i,j,k),2,-1)) .and. &
                        (j-l < ny .or. bc_neumann(mm(i,j,k),2,+1)) ) then
                      dd(ic,jc,kc) = dd(ic,jc,kc) + corner_fac * &
                           fac*res(ileft,jbot,koff) 
                   end if
                   if ( (j+l >  0 .or. bc_neumann(mm(i,j,k),2,-1)) .and. &
                        (j+l < ny .or. bc_neumann(mm(i,j,k),2,+1)) ) then
                      dd(ic,jc,kc) = dd(ic,jc,kc) + corner_fac * &
                           fac*res(ileft,jtop,koff) 
                   end if
                end if

                if ( (i+n >  0 .or. bc_neumann(mm(i,j,k),1,-1)) .and. &
                     (i+n < nx .or. bc_neumann(mm(i,j,k),1,+1)) ) then
                   if ( (j-l >  0 .or. bc_neumann(mm(i,j,k),2,-1)) .and. &
                        (j-l < ny .or. bc_neumann(mm(i,j,k),2,+1)) ) then
                      dd(ic,jc,kc) = dd(ic,jc,kc) + corner_fac * &
                           fac*res(irght,jbot,koff) 
                   end if
                   if ( (j+l >  0 .or. bc_neumann(mm(i,j,k),2,-1)) .and. &
                        (j+l < ny .or. bc_neumann(mm(i,j,k),2,+1)) ) then
                      dd(ic,jc,kc) = dd(ic,jc,kc) + corner_fac * &
                           fac*res(irght,jtop,koff) 
                   end if
                end if

              end do
            end do
          end do
        end do
      end do

      do jc = lod(2),hid(2)
      do ic = lod(1),hid(1)
         i = (ic-lod(1))*ratio(1)
         j = (jc-lod(2))*ratio(2)
         if (.not.bc_dirichlet(mm(i,j,k),1,0)) dd(ic,jc,kc) = ZERO
      end do
      end do

    end if

  end subroutine fine_edge_resid_3d

  subroutine impose_neumann_bcs_2d(uu,mm,lo,ng)

    integer, intent(in) :: ng,lo(:)
    real (kind = dp_t), intent(inout) :: uu(lo(1)-ng:,lo(2)-ng:)
    integer           , intent(in   ) :: mm(lo(1)   :,lo(2)   :)
    integer :: i, j, hi(2)

    hi(1) = lo(1) + size(mm,dim=1)-1
    hi(2) = lo(2) + size(mm,dim=2)-1

    do i = lo(1),hi(1)
       if (bc_neumann(mm(i,lo(2)),2,-1)) uu(i,lo(2)-1) = uu(i,lo(2)+1)
       if (bc_neumann(mm(i,hi(2)),2,+1)) uu(i,hi(2)+1) = uu(i,hi(2)-1)
    end do

    do j = lo(2),hi(2)
       if (bc_neumann(mm(lo(1),j),1,-1)) uu(lo(1)-1,j) = uu(lo(1)+1,j)
       if (bc_neumann(mm(hi(1),j),1,+1)) uu(hi(1)+1,j) = uu(hi(1)-1,j)
    end do

    if (bc_neumann(mm(lo(1),lo(2)),1,-1)) uu(lo(1)-1,lo(2)-1) = uu(lo(1)+1,lo(2)-1) 
    if (bc_neumann(mm(lo(1),lo(2)),2,-1)) uu(lo(1)-1,lo(2)-1) = uu(lo(1)-1,lo(2)+1) 

    if (bc_neumann(mm(hi(1),lo(2)),1,+1)) uu(hi(1)+1,lo(2)-1) = uu(hi(1)-1,lo(2)-1) 
    if (bc_neumann(mm(hi(1),lo(2)),2,-1)) uu(hi(1)+1,lo(2)-1) = uu(hi(1)+1,lo(2)+1) 

    if (bc_neumann(mm(lo(1),hi(2)),1,-1)) uu(lo(1)-1,hi(2)+1) = uu(lo(1)+1,hi(2)+1) 
    if (bc_neumann(mm(lo(1),hi(2)),2,+1)) uu(lo(1)-1,hi(2)+1) = uu(lo(1)-1,hi(2)-1) 

    if (bc_neumann(mm(hi(1),hi(2)),1,+1)) uu(hi(1)+1,hi(2)+1) = uu(hi(1)-1,hi(2)+1) 
    if (bc_neumann(mm(hi(1),hi(2)),2,+1)) uu(hi(1)+1,hi(2)+1) = uu(hi(1)+1,hi(2)-1) 

  end subroutine impose_neumann_bcs_2d

  subroutine impose_neumann_bcs_3d(uu,mm,lo,ng)

    integer, intent(in) :: ng,lo(:)
    real (kind = dp_t), intent(inout) :: uu(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    integer           , intent(in   ) :: mm(lo(1)   :,lo(2)   :,lo(3)   :)
    integer :: i, j, k, hi(3)

    hi(1) = lo(1) + size(mm,dim=1)-1
    hi(2) = lo(2) + size(mm,dim=2)-1
    hi(3) = lo(3) + size(mm,dim=3)-1

!   Faces
    do k = lo(3),hi(3)
    do j = lo(2),hi(2)
      i = lo(1)
      if (bc_neumann(mm(i,j,k),1,-1)) uu(i-1,j,k) = uu(i+1,j,k)
      i = hi(1)
      if (bc_neumann(mm(i,j,k),1,+1)) uu(i+1,j,k) = uu(i-1,j,k)
    end do
    end do
   
    do k = lo(3),hi(3)
    do i = lo(1),hi(1)
      j = lo(2)
      if (bc_neumann(mm(i,j,k),2,-1)) uu(i,j-1,k) = uu(i,j+1,k)
      j = hi(2)
      if (bc_neumann(mm(i,j,k),2,+1)) uu(i,j+1,k) = uu(i,j-1,k)
    end do
    end do
   
    do j = lo(2),hi(2)
    do i = lo(1),hi(1)
      k = lo(3)
      if (bc_neumann(mm(i,j,k),3,-1)) uu(i,j,k-1) = uu(i,j,k+1)
      k = hi(3)
      if (bc_neumann(mm(i,j,k),3,+1)) uu(i,j,k+1) = uu(i,j,k-1)
    end do
    end do

!   Edges
    do i = lo(1),hi(1)
      if (bc_neumann(mm(i,lo(2),lo(3)),2,-1)) uu(i,lo(2)-1,lo(3)-1) = uu(i,lo(2)+1,lo(3)-1)
      if (bc_neumann(mm(i,lo(2),hi(3)),2,-1)) uu(i,lo(2)-1,hi(3)+1) = uu(i,lo(2)+1,hi(3)+1)

      if (bc_neumann(mm(i,hi(2),lo(3)),2,+1)) uu(i,hi(2)+1,lo(3)-1) = uu(i,hi(2)-1,lo(3)-1)
      if (bc_neumann(mm(i,hi(2),hi(3)),2,+1)) uu(i,hi(2)+1,hi(3)+1) = uu(i,hi(2)-1,hi(3)+1)

      if (bc_neumann(mm(i,lo(2),lo(3)),3,-1)) uu(i,lo(2)-1,lo(3)-1) = uu(i,lo(2)-1,lo(3)+1)
      if (bc_neumann(mm(i,hi(2),lo(3)),3,-1)) uu(i,hi(2)+1,lo(3)-1) = uu(i,hi(2)+1,lo(3)+1)

      if (bc_neumann(mm(i,lo(2),hi(3)),3,+1)) uu(i,lo(2)-1,hi(3)+1) = uu(i,lo(2)-1,hi(3)-1)
      if (bc_neumann(mm(i,hi(2),hi(3)),3,+1)) uu(i,hi(2)+1,hi(3)+1) = uu(i,hi(2)+1,hi(3)-1)
    end do

    do j = lo(2),hi(2)
      if (bc_neumann(mm(lo(1),j,lo(3)),1,-1)) uu(lo(1)-1,j,lo(3)-1) = uu(lo(1)+1,j,lo(3)-1)
      if (bc_neumann(mm(lo(1),j,hi(3)),1,-1)) uu(lo(1)-1,j,hi(3)+1) = uu(lo(1)+1,j,hi(3)+1)

      if (bc_neumann(mm(hi(1),j,lo(3)),1,+1)) uu(hi(1)+1,j,lo(3)-1) = uu(hi(1)-1,j,lo(3)-1)
      if (bc_neumann(mm(hi(1),j,hi(3)),1,+1)) uu(hi(1)+1,j,hi(3)+1) = uu(hi(1)-1,j,hi(3)+1)

      if (bc_neumann(mm(lo(1),j,lo(3)),3,-1)) uu(lo(1)-1,j,lo(3)-1) = uu(lo(1)-1,j,lo(3)+1)
      if (bc_neumann(mm(hi(1),j,lo(3)),3,-1)) uu(hi(1)+1,j,lo(3)-1) = uu(hi(1)+1,j,lo(3)+1)

      if (bc_neumann(mm(lo(1),j,hi(3)),3,+1)) uu(lo(1)-1,j,hi(3)+1) = uu(lo(1)-1,j,hi(3)-1)
      if (bc_neumann(mm(hi(1),j,hi(3)),3,+1)) uu(hi(1)+1,j,hi(3)+1) = uu(hi(1)+1,j,hi(3)-1)
    end do

    do k = lo(3),hi(3)
      if (bc_neumann(mm(lo(1),lo(2),k),1,-1)) uu(lo(1)-1,lo(2)-1,k) = uu(lo(1)+1,lo(2)-1,k)
      if (bc_neumann(mm(lo(1),hi(2),k),1,-1)) uu(lo(1)-1,hi(2)+1,k) = uu(lo(1)+1,hi(2)+1,k)

      if (bc_neumann(mm(hi(1),lo(2),k),1,+1)) uu(hi(1)+1,lo(2)-1,k) = uu(hi(1)-1,lo(2)-1,k)
      if (bc_neumann(mm(hi(1),hi(2),k),1,+1)) uu(hi(1)+1,hi(2)+1,k) = uu(hi(1)-1,hi(2)+1,k)

      if (bc_neumann(mm(lo(1),lo(2),k),2,-1)) uu(lo(1)-1,lo(2)-1,k) = uu(lo(1)-1,lo(2)+1,k)
      if (bc_neumann(mm(hi(1),lo(2),k),2,-1)) uu(hi(1)+1,lo(2)-1,k) = uu(hi(1)+1,lo(2)+1,k)

      if (bc_neumann(mm(lo(1),hi(2),k),2,+1)) uu(lo(1)-1,hi(2)+1,k) = uu(lo(1)-1,hi(2)-1,k)
      if (bc_neumann(mm(hi(1),hi(2),k),2,+1)) uu(hi(1)+1,hi(2)+1,k) = uu(hi(1)+1,hi(2)-1,k)

    end do

!   Corners
    i = lo(1)
    j = lo(2)
    k = lo(3)
    if (bc_neumann(mm(i,j,k),1,-1))  uu(i-1,j-1,k-1) = uu(i+1,j-1,k-1)
    if (bc_neumann(mm(i,j,k),2,-1))  uu(i-1,j-1,k-1) = uu(i-1,j+1,k-1)
    if (bc_neumann(mm(i,j,k),3,-1))  uu(i-1,j-1,k-1) = uu(i-1,j-1,k+1)

    i = hi(1)
    j = lo(2)
    k = lo(3)
    if (bc_neumann(mm(i,j,k),1,+1))  uu(i+1,j-1,k-1) = uu(i-1,j-1,k-1)
    if (bc_neumann(mm(i,j,k),2,-1))  uu(i+1,j-1,k-1) = uu(i+1,j+1,k-1)
    if (bc_neumann(mm(i,j,k),3,-1))  uu(i+1,j-1,k-1) = uu(i+1,j-1,k+1)

    i = lo(1)
    j = hi(2)
    k = lo(3)
    if (bc_neumann(mm(i,j,k),1,-1))  uu(i-1,j+1,k-1) = uu(i+1,j+1,k-1)
    if (bc_neumann(mm(i,j,k),2,+1))  uu(i-1,j+1,k-1) = uu(i-1,j-1,k-1)
    if (bc_neumann(mm(i,j,k),3,-1))  uu(i-1,j+1,k-1) = uu(i-1,j+1,k+1)

    i = lo(1)
    j = lo(2)
    k = hi(3)
    if (bc_neumann(mm(i,j,k),1,-1))  uu(i-1,j-1,k+1) = uu(i+1,j-1,k+1)
    if (bc_neumann(mm(i,j,k),2,-1))  uu(i-1,j-1,k+1) = uu(i-1,j+1,k+1)
    if (bc_neumann(mm(i,j,k),3,+1))  uu(i-1,j-1,k+1) = uu(i-1,j-1,k-1)

    i = hi(1)
    j = hi(2)
    k = lo(3)
    if (bc_neumann(mm(i,j,k),1,+1))  uu(i+1,j+1,k-1) = uu(i-1,j+1,k-1)
    if (bc_neumann(mm(i,j,k),2,+1))  uu(i+1,j+1,k-1) = uu(i+1,j-1,k-1)
    if (bc_neumann(mm(i,j,k),3,-1))  uu(i+1,j+1,k-1) = uu(i+1,j+1,k+1)

    i = hi(1)
    j = lo(2)
    k = hi(3)
    if (bc_neumann(mm(i,j,k),1,+1))  uu(i+1,j-1,k+1) = uu(i-1,j-1,k+1)
    if (bc_neumann(mm(i,j,k),2,-1))  uu(i+1,j-1,k+1) = uu(i+1,j+1,k+1)
    if (bc_neumann(mm(i,j,k),3,+1))  uu(i+1,j-1,k+1) = uu(i+1,j-1,k-1)

    i = lo(1)
    j = hi(2)
    k = hi(3)
    if (bc_neumann(mm(i,j,k),1,-1))  uu(i-1,j+1,k+1) = uu(i+1,j+1,k+1)
    if (bc_neumann(mm(i,j,k),2,+1))  uu(i-1,j+1,k+1) = uu(i-1,j-1,k+1)
    if (bc_neumann(mm(i,j,k),3,+1))  uu(i-1,j+1,k+1) = uu(i-1,j+1,k-1)

    i = hi(1)
    j = hi(2)
    k = hi(3)
    if (bc_neumann(mm(i,j,k),1,+1))  uu(i+1,j+1,k+1) = uu(i-1,j+1,k+1)
    if (bc_neumann(mm(i,j,k),2,+1))  uu(i+1,j+1,k+1) = uu(i+1,j-1,k+1)
    if (bc_neumann(mm(i,j,k),3,+1))  uu(i+1,j+1,k+1) = uu(i+1,j+1,k-1)

  end subroutine impose_neumann_bcs_3d

end module stencil_nodal_module

