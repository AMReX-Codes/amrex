module create_umac_grown_module

  use multifab_module
  use bl_constants_module

  private
  public :: create_umac_grown

contains

  ! The goal of create_umac_grown is to take input face-centered normal velocities
  ! and fill the fine ghost cells which are not filled with multifab_fill_boundary.

  ! Here is a 2D illustration.
  ! x and y are input face-centered normal velocities at the fine level.
  ! X and Y are input face-centered normal velocities at the coarse level.
  ! We need to compute x-velocities at points "a" and y-velocities at points "b"

  !      |--y--|--y--|--y--|--y--|--y--|--y--|--y--|--y--|--b--Y-----|
  !      |     |     |     |     |     |     |     |     |           |
  !      |     |     |     |     |     |     |     |     |           |
  !      x     x     x     x     x     x     x     x     x     a     |
  !      |     |     |     |     |     |     |     |     |           |
  !      |     |     |     |     |     |     |     |     |           |
  !      |--y--|--y--|--y--|--y--|--y--|--y--|--y--|--y--X  b        X
  !      |     |     |     |     |     |     |     |     |           |
  !      |     |     |     |     |     |     |     |     |           |
  !      x     x     x     x     x     x     x     x     x     a     |
  !      |     |     |     |     |     |     |     |     |           |
  !      |     |     |     |     |     |     |     |     |           |
  !      |--y--|--y--|--y--|--y--|--y--|--y--|--y--|--y--|--b--Y-----|
  !      |     |     |     |     |     |     |     |     |           |
  !      |     |     |     |     |     |     |     |     |           |
  !      x     x     x     x     x     x     x     x     x     a     |
  !      |     |     |     |     |     |     |     |     |           |
  !      |     |     |     |     |     |     |     |     |           |
  !      |--y--|--y--|--y--|--y--|--y--|--y--|--y--|--y--X  b        X
  !      |     |     |     |     |     |     |     |     |           |
  !      |     |     |     |     |     |     |     |     |           |
  !      x     x     x     x     x     x     x     x     x     a     |
  !      |     |     |     |     |     |     |     |     |           |
  !      |     |     |     |     |     |     |     |     |           |
  !      |--y--Y--y--|--y--Y--y--|--y--Y--y--|--y--Y--y--|--b--Y-----|
  !      |           |           |           |           |           |
  !      |           |           |           |           |           |
  !      a     a     a     a     a     a     a     a     a     a     |
  !      |           |           |           |           |           |
  !      |           |           |           |           |           |
  !      X  b     b  X  b     b  X  b     b  X  b     b  X  b        X
  !      |           |           |           |           |           |
  !      |           |           |           |           |           |
  !      |           |           |           |           |           |
  !      |           |           |           |           |           |
  !      |           |           |           |           |           |
  !      |-----Y-----|-----Y-----|-----Y-----|-----Y-----|-----Y-----|

  subroutine create_umac_grown(finelev,fine,crse,bc_crse,bc_fine)

    use define_bc_module
    use multifab_physbc_edgevel_module, only: multifab_physbc_edgevel

    integer       , intent(in   ) :: finelev
    type(multifab), intent(inout) :: fine(:)
    type(multifab), intent(inout) :: crse(:)
    type(bc_level), intent(in   ) :: bc_crse,bc_fine

    ! local
    integer        :: i,j,k,ng_f,ng_c,dm
    integer        :: c_lo(get_dim(crse(1))),c_hi(get_dim(crse(1)))
    integer        :: f_lo(get_dim(crse(1))),f_hi(get_dim(crse(1)))

    type(fgassoc)  :: fgasc
    type(boxarray) :: f_ba,c_ba,tba
    type(multifab) :: f_mf,c_mf,tcrse,tfine
    type(layout)   :: f_la,c_la,tla,fine_la

    real(kind=dp_t), pointer :: fp(:,:,:,:), cp(:,:,:,:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "create_umac_grown")

    dm = get_dim(crse(1))

    ! Call multifab_fill_boundary and impose_phys_bcs_on_edges on the coarse level.
    ! We will call multifab_fill_boundary and impose_phys_bcs_on_edges on the fine level
    ! later in this subroutine
    ! this can be optimized by putting if (finelev .eq. nlevs) around these lines, but
    ! nlevs isn't available in this subroutine yet
    do i=1,dm
       call multifab_fill_boundary(crse(i))
    end do
    call multifab_physbc_edgevel(crse,bc_crse)

    ! Grab the cached boxarray of all ghost cells not covered by valid region.
    fine_la = get_layout(fine(1))
    fgasc   = layout_fgassoc(fine_la, 1)

    call boxarray_build_copy(f_ba,fgasc%ba)
    call boxarray_build_copy(c_ba,fgasc%ba)

    do i=1,nboxes(f_ba)
       call set_box(c_ba,i,coarsen(get_box(f_ba,i),2))
       call set_box(f_ba,i,refine(get_box(c_ba,i),2))
    end do

    do i=1,dm

       call build(f_la,f_ba,get_pd(get_layout(fine(i))),get_pmask(get_layout(fine(i))))
       call build(c_la,c_ba,get_pd(get_layout(crse(i))),get_pmask(get_layout(crse(i))), &
                  explicit_mapping=get_proc(f_la))

       ! Create c_mf and f_mf on the same proc.
       call build(f_mf,f_la,1,0,nodal_flags(fine(i)))
       call build(c_mf,c_la,1,0,nodal_flags(crse(i)))

       ! Update c_mf with valid and ghost regions from crse.
       call boxarray_build_copy(tba, get_boxarray(crse(i)))
       do j=1,nboxes(tba)
          call set_box(tba,j,grow(get_box(crse(i)%la,j),1))
       end do
       call build(tla, tba, get_pd(get_layout(crse(i))), get_pmask(get_layout(crse(i))), &
                  explicit_mapping = get_proc(get_layout(crse(i))))
       call destroy(tba)
       call build(tcrse, tla, 1, 0, nodal_flags(crse(i)))

       if (nghost(crse(i)) .lt. 1) &
           call bl_error("Need at least one ghost cell in crse in create_umac_grown")

       !$OMP PARALLEL DO PRIVATE(j,fp,cp)
       do j=1,nfabs(tcrse)
          fp => dataptr(tcrse,  j)
          cp => dataptr(crse(i),j, get_ibox(tcrse,j))
          call cpy_d(fp,cp)
       end do
       !$OMP END PARALLEL DO

       call copy(c_mf, tcrse)

       call destroy(tcrse)
       call destroy(tla)
       
       ng_f = nghost(f_mf)
       ng_c = nghost(c_mf)
       !
       ! Fill in some of the fine ghost cells from crse.
       !
       !$OMP PARALLEL DO PRIVATE(j,fp,cp,f_lo,c_lo,c_hi)
       do j=1,nfabs(f_mf)
          fp   => dataptr(f_mf,j)
          cp   => dataptr(c_mf,j)
          f_lo =  lwb(get_box(f_mf,j))
          c_lo =  lwb(get_box(c_mf,j))
          c_hi =  upb(get_box(c_mf,j))
          select case(dm)
          case (2)
             call pc_edge_interp_2d(i,f_lo,c_lo,c_hi,fp(:,:,1,1),ng_f,cp(:,:,1,1),ng_c)
          case (3)
             call pc_edge_interp_3d(i,f_lo,c_lo,c_hi,fp(:,:,:,1),ng_f,cp(:,:,:,1),ng_c)
          end select
       end do
       !$OMP END PARALLEL DO

       call copy(f_mf,fine(i))
       !
       ! Fill in the rest of the fine ghost cells.
       !
       !$OMP PARALLEL DO PRIVATE(j,fp,f_lo,c_lo,c_hi)
       do j=1,nfabs(f_mf)
          fp   => dataptr(f_mf,j)
          f_lo =  lwb(get_box(f_mf,j))
          c_lo =  lwb(get_box(c_mf,j))
          c_hi =  upb(get_box(c_mf,j))
          select case(dm)
          case (2)
             call lin_edge_interp_2d(i,f_lo,c_lo,c_hi,fp(:,:,1,1),ng_f)
          case (3)
             call lin_edge_interp_3d(i,f_lo,c_lo,c_hi,fp(:,:,:,1),ng_f)
          end select
       end do
       !$OMP END PARALLEL DO

       call destroy(c_mf)
       call destroy(c_la)

       ! Update ghost regions of fine where they overlap with f_mf.
       call boxarray_build_copy(tba, get_boxarray(fine(i)))
       do j = 1, nboxes(tba)
          call set_box(tba,j,grow(get_box(fine(i)%la,j),nghost(fine(i))))
       end do
       call build(tla, tba, get_pd(get_layout(fine(i))), get_pmask(get_layout(fine(i))), &
                  explicit_mapping=get_proc(get_layout(fine(i))))
       call destroy(tba)
       call build(tfine, tla, 1, 0, nodal_flags(fine(i)))

       call copy(tfine, f_mf)

       call destroy(f_mf)
       call destroy(f_la)

       !$OMP PARALLEL DO PRIVATE(j,k,fp,cp,tba)
       do j=1,nfabs(tfine)
          call boxarray_box_diff(tba, get_ibox(tfine,j), get_ibox(fine(i),j))
          do k = 1, nboxes(tba)
             fp => dataptr(fine(i), j, get_box(tba,k))
             cp => dataptr(tfine,   j, get_box(tba,k))
             call cpy_d(fp,cp)
          end do
          call destroy(tba)
       end do
       !$OMP END PARALLEL DO

       call destroy(tfine)
       call destroy(tla)

       ng_f = nghost(fine(1))

       call multifab_fill_boundary(fine(i))
       !
       ! Now fix up umac grown due to the low order interpolation we used.
       !
       !$OMP PARALLEL DO PRIVATE(j,fp,f_lo,f_hi)
       do j=1, nfabs(fine(i))
          fp   => dataptr(fine(i), j)
          f_lo =  lwb(get_box(fine(i), j))
          f_hi =  upb(get_box(fine(i), j))
          select case(dm)
          case (2)
             call correct_umac_grown_2d(fp(:,:,1,1),ng_f,f_lo,f_hi,i)
          case (3)
             call correct_umac_grown_3d(fp(:,:,:,1),ng_f,f_lo,f_hi,i)
          end select
       end do
       !$OMP END PARALLEL DO

       call multifab_fill_boundary(fine(i))
    end do

    call multifab_physbc_edgevel(fine,bc_fine)

    call destroy(f_ba)
    call destroy(c_ba)
    call destroy(bpt)

  end subroutine create_umac_grown

  subroutine pc_edge_interp_2d(dir,f_lo,c_lo,c_hi,fine,ng_f,crse,ng_c)

    integer,         intent(in   ) :: dir,f_lo(:),c_lo(:),c_hi(:)
    integer,         intent(in   ) :: ng_f,ng_c
    real(kind=dp_t), intent(inout) :: fine(f_lo(1)-ng_f:,f_lo(2)-ng_f:)
    real(kind=dp_t), intent(inout) :: crse(c_lo(1)-ng_c:,c_lo(2)-ng_c:)

    ! local
    integer :: i,ii,j,jj

    if (dir .eq. 1) then

       do j=c_lo(2),c_hi(2)
          do i=c_lo(1),c_hi(1)+1
             do jj=0,1
                fine(2*i,2*j+jj) = crse(i,j)
             end do
          end do
       end do

    else

       do j=c_lo(2),c_hi(2)+1
          do i=c_lo(1),c_hi(1)
             do ii=0,1
                fine(2*i+ii,2*j) = crse(i,j)
             end do
          end do
       end do

    end if    

  end subroutine pc_edge_interp_2d

  subroutine pc_edge_interp_3d(dir,f_lo,c_lo,c_hi,fine,ng_f,crse,ng_c)

    integer,         intent(in   ) :: dir,f_lo(:),c_lo(:),c_hi(:)
    integer,         intent(in   ) :: ng_f,ng_c
    real(kind=dp_t), intent(inout) :: fine(f_lo(1)-ng_f:,f_lo(2)-ng_f:,f_lo(3)-ng_f:)
    real(kind=dp_t), intent(inout) :: crse(c_lo(1)-ng_c:,c_lo(2)-ng_c:,c_lo(3)-ng_c:)

    ! local
    integer :: i,ii,j,jj,k,kk

    if (dir .eq. 1) then

       do k=c_lo(3),c_hi(3)
          do j=c_lo(2),c_hi(2)
             do i=c_lo(1),c_hi(1)+1
                do kk=0,1
                   do jj=0,1
                      fine(2*i,2*j+jj,2*k+kk) = crse(i,j,k)
                   end do
                end do
             end do
          end do
       end do

    else if (dir .eq. 2) then

       do k=c_lo(3),c_hi(3)
          do j=c_lo(2),c_hi(2)+1
             do i=c_lo(1),c_hi(1)
                do kk=0,1
                   do ii=0,1
                      fine(2*i+ii,2*j,2*k+kk) = crse(i,j,k)
                   end do
                end do
             end do
          end do
       end do

    else

       do k=c_lo(3),c_hi(3)+1
          do j=c_lo(2),c_hi(2)
             do i=c_lo(1),c_hi(1)
                do jj=0,1
                   do ii=0,1
                      fine(2*i+ii,2*j+jj,2*k) = crse(i,j,k)
                   end do
                end do
             end do
          end do
       end do

    end if    

  end subroutine pc_edge_interp_3d

  subroutine lin_edge_interp_2d(dir,f_lo,c_lo,c_hi,fine,ng_f)

    integer,         intent(in   ) :: dir,f_lo(:),c_lo(:),c_hi(:),ng_f
    real(kind=dp_t), intent(inout) :: fine(f_lo(1)-ng_f:,f_lo(2)-ng_f:)

    ! local
    integer :: i,ii,j,jj

    if (dir .eq. 1) then

       do j=c_lo(2),c_hi(2)
          do i=c_lo(1),c_hi(1)
             do jj=0,1
                fine(2*i+1,2*j+jj) = HALF*(fine(2*i,2*j+jj)+fine(2*i+2,2*j+jj))
             end do
          end do
       end do

    else

       do j=c_lo(2),c_hi(2)
          do i=c_lo(1),c_hi(1)
             do ii=0,1
                fine(2*i+ii,2*j+1) = HALF*(fine(2*i+ii,2*j)+fine(2*i+ii,2*j+2))
             end do
          end do
       end do

    end if    

  end subroutine lin_edge_interp_2d

  subroutine lin_edge_interp_3d(dir,f_lo,c_lo,c_hi,fine,ng_f)

    integer,         intent(in   ) :: dir,f_lo(:),c_lo(:),c_hi(:),ng_f
    real(kind=dp_t), intent(inout) :: fine(f_lo(1)-ng_f:,f_lo(2)-ng_f:,f_lo(3)-ng_f:)

    ! local
    integer :: i,ii,j,jj,k,kk

    if (dir .eq. 1) then

       do k=c_lo(3),c_hi(3)
          do j=c_lo(2),c_hi(2)
             do i=c_lo(1),c_hi(1)
                do kk=0,1
                   do jj=0,1
                      fine(2*i+1,2*j+jj,2*k+kk) = &
                           HALF*(fine(2*i,2*j+jj,2*k+kk)+fine(2*i+2,2*j+jj,2*k+kk))
                   end do
                end do
             end do
          end do
       end do

    else if (dir .eq. 2) then

       do k=c_lo(3),c_hi(3)
          do j=c_lo(2),c_hi(2)
             do i=c_lo(1),c_hi(1)
                do kk=0,1
                   do ii=0,1
                      fine(2*i+ii,2*j+1,2*k+kk) = &
                           HALF*(fine(2*i+ii,2*j,2*k+kk)+fine(2*i+ii,2*j+2,2*k+kk))
                   end do
                end do
             end do
          end do
       end do

    else

       do k=c_lo(3),c_hi(3)
          do j=c_lo(2),c_hi(2)
             do i=c_lo(1),c_hi(1)
                do jj=0,1
                   do ii=0,1
                      fine(2*i+ii,2*j+jj,2*k+1) = &
                           HALF*(fine(2*i+ii,2*j+jj,2*k)+fine(2*i+ii,2*j+jj,2*k+2))
                   end do
                end do
             end do
          end do
       end do

    end if    

  end subroutine lin_edge_interp_3d

  subroutine correct_umac_grown_2d(vel,ng_f,lo,hi,dir)
    
    integer        , intent(in   ) :: lo(:),hi(:),dir,ng_f
    real(kind=dp_t), intent(inout) :: vel(lo(1)-ng_f:,lo(2)-ng_f:)

    ! local
    integer         :: i,j,signx,signy
    real(kind=dp_t) :: temp_velx_lo(lo(2)-1:hi(2)+1)
    real(kind=dp_t) :: temp_velx_hi(lo(2)-1:hi(2)+1)
    real(kind=dp_t) :: temp_vely_lo(lo(1)-1:hi(1)+1)
    real(kind=dp_t) :: temp_vely_hi(lo(1)-1:hi(1)+1)

    if (dir .eq. 1) then

       ! for each normal velocity in the first fine ghost cell in the normal direction
       ! compute what the coarse velocity was that came from the first coarse ghost cell
       ! store this value in the first fine ghost cell
       do j=lo(2)-1,hi(2)+1
          vel(lo(1)-1,j) = TWO*vel(lo(1)-1,j) - vel(lo(1),j)
          vel(hi(1)+2,j) = TWO*vel(hi(1)+2,j) - vel(hi(1)+1,j)
       end do

       ! store the coarse velocity in a temporary array
       temp_velx_lo(lo(2)-1:hi(2)+1) = vel(lo(1)-1,lo(2)-1:hi(2)+1)
       temp_velx_hi(lo(2)-1:hi(2)+1) = vel(hi(1)+2,lo(2)-1:hi(2)+1)

       ! linearly interpolate to obtain a better estimate of the velocity from
       ! the first coarse ghost cell
       ! store this value in the first fine ghost cell
       do j=lo(2)-1,hi(2)+1
          if (abs(mod(j,2)) .eq. 1) then
             signy = 1
          else
             signy = -1
          end if
          vel(lo(1)-1,j) = THREE4TH*vel(lo(1)-1,j) + FOURTH*temp_velx_lo(j+signy)
          vel(hi(1)+2,j) = THREE4TH*vel(hi(1)+2,j) + FOURTH*temp_velx_hi(j+signy)
       end do

       ! average the grid edge value with the velocity from the first coarse ghost cell
       ! (which is currently stored in the first fine ghost cell)
       ! to get a better estimate of the first fine ghost cell
       do j=lo(2)-1,hi(2)+1
          vel(lo(1)-1,j) = HALF*(vel(lo(1)-1,j)+vel(lo(1),j))
          vel(hi(1)+2,j) = HALF*(vel(hi(1)+2,j)+vel(hi(1)+1,j))
       end do

       ! at transverse faces, the first fine ghost value was set to the 
       ! first coarse ghost cell value
       ! we linearly interpolate to get a better estimate of this value
       do i=lo(1),hi(1)+1
          vel(i,lo(2)-1) = THREE4TH*vel(i,lo(2)-1) + EIGHTH*(vel(i,lo(2))+vel(i,lo(2)+1))
          vel(i,hi(2)+1) = THREE4TH*vel(i,hi(2)+1) + EIGHTH*(vel(i,hi(2))+vel(i,hi(2)-1))
       end do

    else

       ! for each normal velocity in the first fine ghost cell in the normal direction
       ! compute what the coarse velocity was that came from the first coarse ghost cell
       ! store this value in the first fine ghost cell
       do i=lo(1)-1,hi(1)+1
          vel(i,lo(2)-1) = TWO*vel(i,lo(2)-1) - vel(i,lo(2))
          vel(i,hi(2)+2) = TWO*vel(i,hi(2)+2) - vel(i,hi(2)+1)
       end do

       ! store the coarse velocity in a temporary array
       temp_vely_lo(lo(1)-1:hi(1)+1) = vel(lo(1)-1:hi(1)+1,lo(2)-1)
       temp_vely_hi(lo(1)-1:hi(1)+1) = vel(lo(1)-1:hi(1)+1,hi(2)+2)

       ! linearly interpolate to obtain a better estimate of the velocity from
       ! the first coarse ghost cell
       ! store this value in the first fine ghost cell
       do i=lo(1)-1,hi(1)+1
          if (abs(mod(i,2)) .eq. 1) then
             signx = 1
          else
             signx = -1
          end if
          vel(i,lo(2)-1) = THREE4TH*vel(i,lo(2)-1) + FOURTH*temp_vely_lo(i+signx)
          vel(i,hi(2)+2) = THREE4TH*vel(i,hi(2)+2) + FOURTH*temp_vely_hi(i+signx)
       end do

       ! average the grid edge value with the velocity from the first coarse ghost cell
       ! (which is currently stored in the first fine ghost cell)
       ! to get a better estimate of the first fine ghost cell
       do i=lo(1)-1,hi(1)+1
          vel(i,lo(2)-1) = HALF*(vel(i,lo(2)-1)+vel(i,lo(2)))
          vel(i,hi(2)+2) = HALF*(vel(i,hi(2)+2)+vel(i,hi(2)+1))
       end do

       ! at transverse faces, the first fine ghost value was set to the 
       ! first coarse ghost cell value
       ! we linearly interpolate to get a better estimate of this value
       do j=lo(2),hi(2)+1
          vel(lo(1)-1,j) = THREE4TH*vel(lo(1)-1,j) + EIGHTH*(vel(lo(1),j)+vel(lo(1)+1,j))
          vel(hi(1)+1,j) = THREE4TH*vel(hi(1)+1,j) + EIGHTH*(vel(hi(1),j)+vel(hi(1)-1,j))
       end do

    end if

  end subroutine correct_umac_grown_2d

  subroutine correct_umac_grown_3d(vel,ng_f,lo,hi,dir)
    
    integer        , intent(in   ) :: lo(:),hi(:),dir,ng_f
    real(kind=dp_t), intent(inout) :: vel(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:)

    integer         :: i,j,k,signx,signy,signz
    ! real(kind=dp_t) :: temp_velx_lo(lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)
    ! real(kind=dp_t) :: temp_velx_hi(lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)
    ! real(kind=dp_t) :: temp_vely_lo(lo(1)-1:hi(1)+1,lo(3)-1:hi(3)+1)
    ! real(kind=dp_t) :: temp_vely_hi(lo(1)-1:hi(1)+1,lo(3)-1:hi(3)+1)
    ! real(kind=dp_t) :: temp_velz_lo(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1)
    ! real(kind=dp_t) :: temp_velz_hi(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1)

     real(kind=dp_t), allocatable :: temp_velx_lo(:,:)
     real(kind=dp_t), allocatable :: temp_velx_hi(:,:)
     real(kind=dp_t), allocatable :: temp_vely_lo(:,:)
     real(kind=dp_t), allocatable :: temp_vely_hi(:,:)
     real(kind=dp_t), allocatable :: temp_velz_lo(:,:)
     real(kind=dp_t), allocatable :: temp_velz_hi(:,:)

     real(kind=dp_t), parameter :: onesixteenth    = (1.d0/16.d0)
     real(kind=dp_t), parameter :: threesixteenths = (3.d0/16.d0)
     real(kind=dp_t), parameter :: ninesixteenths  = (9.d0/16.d0)

    if (dir .eq. 1) then

       allocate(temp_velx_lo(lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))
       allocate(temp_velx_hi(lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))

       ! for each normal velocity in the first fine ghost cell in the normal direction
       ! compute what the coarse velocity was that came from the first coarse ghost cell
       ! store this value in the first fine ghost cell
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             vel(lo(1)-1,j,k) = TWO*vel(lo(1)-1,j,k) - vel(lo(1),j,k)
             vel(hi(1)+2,j,k) = TWO*vel(hi(1)+2,j,k) - vel(hi(1)+1,j,k)
          end do
       end do

       ! store the coarse velocity in a temporary array
       temp_velx_lo(lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1) = &
            vel(lo(1)-1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)
       temp_velx_hi(lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1) = &
            vel(hi(1)+2,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)

       ! linearly interpolate to obtain a better estimate of the velocity from
       ! the first coarse ghost cell
       ! store this value in the first fine ghost cell
       do k=lo(3)-1,hi(3)+1
          if (abs(mod(k,2)) .eq. 1) then
             signz = 1
          else
             signz = -1
          end if
          do j=lo(2)-1,hi(2)+1
             if (abs(mod(j,2)) .eq. 1) then
                signy = 1
             else
                signy = -1
             end if
             vel(lo(1)-1,j,k) = ninesixteenths*vel(lo(1)-1,j,k) &
                  + threesixteenths*temp_velx_lo(j+signy,k) &
                  + threesixteenths*temp_velx_lo(j,k+signz) &
                  + onesixteenth*temp_velx_lo(j+signy,k+signz)
             vel(hi(1)+2,j,k) = ninesixteenths*vel(hi(1)+2,j,k) &
                  + threesixteenths*temp_velx_hi(j+signy,k) &
                  + threesixteenths*temp_velx_hi(j,k+signz) &
                  + onesixteenth*temp_velx_hi(j+signy,k+signz)
          end do
       end do

       ! average the grid edge value with the velocity from the first coarse ghost cell
       ! (which is currently stored in the first fine ghost cell)
       ! to get a better estimate of the first fine ghost cell
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             vel(lo(1)-1,j,k) = HALF*(vel(lo(1)-1,j,k)+vel(lo(1),j,k))
             vel(hi(1)+2,j,k) = HALF*(vel(hi(1)+2,j,k)+vel(hi(1)+1,j,k))
          end do
       end do

       ! at transverse faces, the first fine ghost value was set to the 
       ! first coarse ghost cell value
       ! we linearly interpolate to get a better estimate of this value
       do j=lo(2)-1,hi(2)+1
          do i=lo(1),hi(1)+1
             vel(i,j,lo(3)-1) = THREE4TH*vel(i,j,lo(3)-1) &
                  + EIGHTH*(vel(i,j,lo(3))+vel(i,j,lo(3)+1))
             vel(i,j,hi(3)+1) = THREE4TH*vel(i,j,hi(3)+1) &
                  + EIGHTH*(vel(i,j,hi(3))+vel(i,j,hi(3)-1))
          end do
       end do

       do k=lo(3)-1,hi(3)+1
          do i=lo(1),hi(1)+1
             vel(i,lo(2)-1,k) = THREE4TH*vel(i,lo(2)-1,k) &
                  + EIGHTH*(vel(i,lo(2),k)+vel(i,lo(2)+1,k))
             vel(i,hi(2)+1,k) = THREE4TH*vel(i,hi(2)+1,k) &
                  + EIGHTH*(vel(i,hi(2),k)+vel(i,hi(2)-1,k))
          end do
       end do

       deallocate(temp_velx_lo, temp_velx_hi)

    else if (dir .eq. 2) then

       allocate(temp_vely_lo(lo(1)-1:hi(1)+1,lo(3)-1:hi(3)+1))
       allocate(temp_vely_hi(lo(1)-1:hi(1)+1,lo(3)-1:hi(3)+1))

       ! for each normal velocity in the first fine ghost cell in the normal direction
       ! compute what the coarse velocity was that came from the first coarse ghost cell
       ! store this value in the first fine ghost cell
       do k=lo(3)-1,hi(3)+1
          do i=lo(1)-1,hi(1)+1
             vel(i,lo(2)-1,k) = TWO*vel(i,lo(2)-1,k) - vel(i,lo(2),k)
             vel(i,hi(2)+2,k) = TWO*vel(i,hi(2)+2,k) - vel(i,hi(2)+1,k)
          end do
       end do

       ! store the coarse velocity in a temporary array
       temp_vely_lo(lo(1)-1:hi(1)+1,lo(3)-1:hi(3)+1) = &
            vel(lo(1)-1:hi(1)+1,lo(2)-1,lo(3)-1:hi(3)+1)
       temp_vely_hi(lo(1)-1:hi(1)+1,lo(3)-1:hi(3)+1) = &
            vel(lo(1)-1:hi(1)+1,hi(2)+2,lo(3)-1:hi(3)+1)

       ! linearly interpolate to obtain a better estimate of the velocity from
       ! the first coarse ghost cell
       ! store this value in the first fine ghost cell
       do k=lo(3)-1,hi(3)+1
          if (abs(mod(k,2)) .eq. 1) then
             signz = 1
          else
             signz = -1
          end if
          do i=lo(1)-1,hi(1)+1
             if (abs(mod(i,2)) .eq. 1) then
                signx = 1
             else
                signx = -1
             end if
             vel(i,lo(2)-1,k) = ninesixteenths*vel(i,lo(2)-1,k) &
                  + threesixteenths*temp_vely_lo(i+signx,k) &
                  + threesixteenths*temp_vely_lo(i,k+signz) &
                  + onesixteenth*temp_vely_lo(i+signx,k+signz)
             vel(i,hi(2)+2,k) = ninesixteenths*vel(i,hi(2)+2,k) &
                  + threesixteenths*temp_vely_hi(i+signx,k) &
                  + threesixteenths*temp_vely_hi(i,k+signz) &
                  + onesixteenth*temp_vely_hi(i+signx,k+signz)
          end do
       end do

       ! average the grid edge value with the velocity from the first coarse ghost cell
       ! (which is currently stored in the first fine ghost cell)
       ! to get a better estimate of the first fine ghost cell
       do k=lo(3)-1,hi(3)+1
          do i=lo(1)-1,hi(1)+1
             vel(i,lo(2)-1,k) = HALF*(vel(i,lo(2)-1,k)+vel(i,lo(2),k))
             vel(i,hi(2)+2,k) = HALF*(vel(i,hi(2)+2,k)+vel(i,hi(2)+1,k))
          end do
       end do

       ! at transverse faces, the first fine ghost value was set to the 
       ! first coarse ghost cell value
       ! we linearly interpolate to get a better estimate of this value
       do j=lo(2),hi(2)+1
          do i=lo(1)-1,hi(1)+1
             vel(i,j,lo(3)-1) = THREE4TH*vel(i,j,lo(3)-1) &
                  + EIGHTH*(vel(i,j,lo(3))+vel(i,j,lo(3)+1))
             vel(i,j,hi(3)+1) = THREE4TH*vel(i,j,hi(3)+1) &
                  + EIGHTH*(vel(i,j,hi(3))+vel(i,j,hi(3)-1))
          end do
       end do

       do k=lo(3)-1,hi(3)+1
          do j=lo(2),hi(2)+1
             vel(lo(1)-1,j,k) = THREE4TH*vel(lo(1)-1,j,k) &
                  + EIGHTH*(vel(lo(1),j,k)+vel(lo(1)+1,j,k))
             vel(hi(1)+1,j,k) = THREE4TH*vel(hi(1)+1,j,k) &
                  + EIGHTH*(vel(hi(1),j,k)+vel(hi(1)-1,j,k))
          end do
       end do

       deallocate(temp_vely_lo, temp_vely_hi)

    else

       allocate(temp_velz_lo(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1))
       allocate(temp_velz_hi(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1))

       ! for each normal velocity in the first fine ghost cell in the normal direction
       ! compute what the coarse velocity was that came from the first coarse ghost cell
       ! store this value in the first fine ghost cell
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1
             vel(i,j,lo(3)-1) = TWO*vel(i,j,lo(3)-1) - vel(i,j,lo(3))
             vel(i,j,hi(3)+2) = TWO*vel(i,j,hi(3)+2) - vel(i,j,hi(3)+1)
          end do
       end do

       ! store the coarse velocity in a temporary array
       temp_velz_lo(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1) = &
            vel(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1)
       temp_velz_hi(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1) = &
            vel(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,hi(3)+2)

       ! linearly interpolate to obtain a better estimate of the velocity from
       ! the first coarse ghost cell
       ! store this value in the first fine ghost cell
       do j=lo(2)-1,hi(2)+1
          if (abs(mod(j,2)) .eq. 1) then
             signy = 1
          else
             signy = -1
          end if
          do i=lo(1)-1,hi(1)+1
             if (abs(mod(i,2)) .eq. 1) then
                signx = 1
             else
                signx = -1
             end if
             vel(i,j,lo(3)-1) = ninesixteenths*vel(i,j,lo(3)-1) &
                  + threesixteenths*temp_velz_lo(i+signx,j) &
                  + threesixteenths*temp_velz_lo(i,j+signy) &
                  + onesixteenth*temp_velz_lo(i+signx,j+signy)
             vel(i,j,hi(3)+2) = ninesixteenths*vel(i,j,hi(3)+2) &
                  + threesixteenths*temp_velz_hi(i+signx,j) &
                  + threesixteenths*temp_velz_hi(i,j+signy) &
                  + onesixteenth*temp_velz_hi(i+signx,j+signy)
          end do
       end do

       ! average the grid edge value with the velocity from the first coarse ghost cell
       ! (which is currently stored in the first fine ghost cell)
       ! to get a better estimate of the first fine ghost cell
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1
             vel(i,j,lo(3)-1) = HALF*(vel(i,j,lo(3)-1)+vel(i,j,lo(3)))
             vel(i,j,hi(3)+2) = HALF*(vel(i,j,hi(3)+2)+vel(i,j,hi(3)+1))
          end do
       end do

       ! at transverse faces, the first fine ghost value was set to the 
       ! first coarse ghost cell value
       ! we linearly interpolate to get a better estimate of this value
       do k=lo(3),hi(3)+1
          do i=lo(1)-1,hi(1)+1
             vel(i,lo(2)-1,k) = THREE4TH*vel(i,lo(2)-1,k) &
                  + EIGHTH*(vel(i,lo(2),k)+vel(i,lo(2)+1,k))
             vel(i,hi(2)+1,k) = THREE4TH*vel(i,hi(2)+1,k) &
                  + EIGHTH*(vel(i,hi(2),k)+vel(i,hi(2)-1,k))
          end do
       end do

       do k=lo(3),hi(3)+1
          do j=lo(2)-1,hi(2)+1
             vel(lo(1)-1,j,k) = THREE4TH*vel(lo(1)-1,j,k) &
                  + EIGHTH*(vel(lo(1),j,k)+vel(lo(1)+1,j,k))
             vel(hi(1)+1,j,k) = THREE4TH*vel(hi(1)+1,j,k) &
                  + EIGHTH*(vel(hi(1),j,k)+vel(hi(1)-1,j,k))
          end do
       end do

       deallocate(temp_velz_lo, temp_velz_hi)

    end if

  end subroutine correct_umac_grown_3d

end module create_umac_grown_module
