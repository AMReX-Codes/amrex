module sparse_solve_module

  use bl_types
  use multifab_module
  use stencil_module
  use sort_box_module

  implicit none

  type sparse_matrix
     integer :: nrow = -1
     integer :: numa = -1
     integer   , pointer:: ia(:) => Null()
     integer   , pointer:: ja(:) => Null()
     real(dp_t), pointer:: aa(:) => Null()
  end type sparse_matrix

  type sparse
     real(dp_t) :: Anorm = 0.0_dp_t
     type(imultifab) index_into_aa 
     type(sparse_matrix) :: smt
     type(sparse_matrix) :: sil
  end type sparse

  interface built_q
     module procedure sparse_built_q
  end interface

  interface destroy
     module procedure sparse_destroy
  end interface

  real(kind=dp_t), private, parameter :: zero = 0.0_dp_t
  real(kind=dp_t), private, parameter :: TWO = 2.0_dp_t

contains

  subroutine sparse_matrix_build(sm, numa, nrow)
    type(sparse_matrix), intent(out) :: sm
    integer, intent(in) :: numa, nrow
    sm%nrow = nrow
    sm%numa = numa
    allocate(sm%aa(numa), sm%ja(numa), sm%ia(nrow+1))
  end subroutine sparse_matrix_build

  subroutine sparse_matrix_destroy(sm)
    type(sparse_matrix), intent(inout) :: sm
    sm%nrow = 0
    sm%numa = 0
    deallocate(sm%aa, sm%ja, sm%ia)
  end subroutine sparse_matrix_destroy

  function sparse_built_q(spo) result(r)
    logical :: r
    type(sparse), intent(in) :: spo
    r = associated(spo%smt%ia)
  end function sparse_built_q

  subroutine sparse_destroy(spo)
    type(sparse), intent(inout) :: spo
    call sparse_matrix_destroy(spo%smt)
    call sparse_matrix_destroy(spo%sil)
    spo%Anorm = 0.0_dp_t
    call destroy(spo%index_into_aa)
  end subroutine sparse_destroy

  subroutine sparse_build_st(spo, st, order, verbose)
    type(sparse), intent(out) :: spo
    type(stencil), intent(in) :: st
    integer, intent(in) :: order, verbose
    type(layout) :: la
    la = get_layout(st%ss)
    call sparse_build(spo, st%ss, st%mm, la, order, verbose)
  end subroutine sparse_build_st

  subroutine sparse_build(spo, ss, mm, la, order, verbose)
    type(   layout), intent(in) :: la
    type(imultifab), intent(in) :: mm
    type(multifab) , intent(in) :: ss
    type(sparse)                :: spo
    integer, intent(in)         :: verbose

    integer        , pointer :: mp(:,:,:,:)
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    real(kind=dp_t), allocatable :: temp_aa(:)
    type(box) :: b,sbx,mbx,ibx,jbx,shifted_box
    type(box) :: bbox

    type(box), allocatable :: barray(:)
    integer  , allocatable :: iarray(:)
    integer  , allocatable :: tmp_iarray(:)
    integer  , allocatable :: new_iarray_i(:)
    integer  , allocatable :: new_iarray_ij(:,:)
    integer  , allocatable :: new_iarray_ijk(:,:,:)

    integer :: order

    integer :: numpts, num_aa, num_bdy_pts
    integer :: dm
    integer :: ig,igrid,jgrid,ngrids
    integer :: i,j,k,n,ns,nx,ny,nz
    integer :: dir
    integer :: jlo, jhi, klo, khi
    integer :: inode, iedge

    integer i_lbnd(ss%dim)
    integer i_ubnd(ss%dim)

    integer max_neighbors
    integer, allocatable :: neighbors(:,:)
    integer, allocatable :: sides(:,:)
    integer, allocatable :: num_nbors(:)

    integer, allocatable :: num_grids_for_j(:)
    integer, allocatable :: num_grids_for_jk(:,:)

    integer    , pointer :: ind(:,:,:,:)

    ns = ss%nc
    dm = ss%dim

    ngrids = ss%nboxes

    spo%Anorm = stencil_norm(ss)

    call imultifab_build(spo%index_into_aa,la,1,2)
    call setval(spo%index_into_aa,-1,all=.true.)


    select case(dm)
    case(1)
       max_neighbors = get_max_neighbors_1d(ss)
    case(2)
       max_neighbors = get_max_neighbors_2d(ss)
    case(3)
       max_neighbors = get_max_neighbors_3d(ss)
    end select

    if ( verbose > 0 ) then
       print *, 'MAX_NEIGHBORS = ', max_neighbors
    end if

    allocate(neighbors(ngrids,max_neighbors))
    allocate(    sides(ngrids,max_neighbors))
    allocate(num_nbors(ngrids))

    num_nbors = 0
    neighbors = -1

    select case(dm)
    case(1)
       call make_neighbors_1d(ss, neighbors, sides, num_nbors)
    case(2)
       call make_neighbors_2d(ss, neighbors, sides, num_nbors)
    case(3)
       call make_neighbors_3d(ss, neighbors, sides, num_nbors)
    end select

    numpts = 0
    num_aa = 0
    num_bdy_pts = 0

    do igrid = 1, ngrids
       ibx = get_box(ss, igrid)
       i_lbnd = ibx%lo(1:dm)
       i_ubnd = ibx%hi(1:dm)

       nx = box_extent_d(ibx,1)
       select case(dm)
       case(1)
          numpts = numpts +    nx
          num_aa = num_aa + ns*nx - 2
          num_bdy_pts = num_bdy_pts + 2
       case(2)
          ny = box_extent_d(ibx,2)
          numpts = numpts +    nx*ny
          num_aa = num_aa + ns*nx*ny - 2*(nx+ny)
          num_bdy_pts = num_bdy_pts  + 2*(nx+ny)
       case(3)
          ny = box_extent_d(ibx,2)
          nz = box_extent_d(ibx,3)
          numpts = numpts +    nx*ny*nz
          num_aa = num_aa + ns*nx*ny*nz - 2*(nx*ny + ny*nz + nx*nz)
          num_bdy_pts = num_bdy_pts     + 2*(nx*ny + ny*nz + nx*nz)
       end select

       !      Add back in any fine-fine edges
       do j = 1,num_nbors(igrid)
          jgrid = neighbors(igrid,j)
          jbx = get_box(ss, jgrid)
          dir = sides(igrid,j)
          shifted_box = shift(ibx,sign(1,dir),abs(dir))
          num_aa = num_aa + volume(intersection(shifted_box, jbx))
          num_bdy_pts = num_bdy_pts - volume(intersection(shifted_box, jbx))
       end do
    end do
    if ( verbose > 0 ) then
       print *,'NUMPTS NUM_AA NUM_BDY_PTS ',numpts, num_aa, num_bdy_pts
    end if

    allocate(temp_aa(num_aa+num_bdy_pts*(ns+order-2)))
    allocate(spo%smt%ia(numpts+1))

    !   Create the minimum box containing all the grids
    bbox = get_box(ss, 1)
    do igrid = 2,ngrids
       bbox = box_bbox(bbox,get_box(ss,igrid))
    end do

    !   Convert the ss into CSR format.

    select case(dm)
    case(1)
       iedge = 1
       inode = 1

       allocate(barray(ngrids))
       allocate(iarray(ngrids))
       allocate(tmp_iarray(ngrids))
       allocate(new_iarray_i(ngrids))

       do igrid = 1,ngrids
          barray(igrid) = get_box(ss,igrid)
          iarray(igrid) = igrid
       end do
       call heapsort_indirect_box(barray,tmp_iarray)
       do igrid = 1,ngrids
          new_iarray_i(igrid) = iarray(tmp_iarray(igrid))
       end do

       deallocate(barray)
       deallocate(iarray)
       deallocate(tmp_iarray)

       do igrid = 1, ngrids
          ig = new_iarray_i(igrid)
          sp => dataptr(ss, ig)
          mp => dataptr(mm     , ig)
          ind => dataptr(spo%index_into_aa, ig)
          call create_aa_1d(sp(:,1,1,:), &
               mp(:,1,1,1),              &
               temp_aa,                  &
               spo%smt%ia,     & 
               ind(:,1,1,1),             &
               inode, iedge)
       end do
       spo%smt%ia(inode) = iedge
       if ( verbose > 0 ) then
          print *,'FINAL INODE IEDGE ',inode,iedge
       end if

       num_aa = iedge - 1
       allocate(spo%smt%aa(num_aa))
       do i = 1, num_aa
          spo%smt%aa(i) = temp_aa(i)
       end do

       if ( verbose > 0 ) then
          print *,'REAL NUMPTS NUM_AA ',numpts, num_aa
       end if

       call imultifab_fill_boundary(spo%index_into_aa)

       allocate(spo%smt%ja(num_aa))

       iedge = 1
       do igrid = 1, ngrids
          ig = new_iarray_i(igrid)
          mp => dataptr(mm     , ig)
          sbx = get_pbox(spo%index_into_aa,ig)
          ind => dataptr(spo%index_into_aa,ig)
          call create_ja_1d(spo%smt%ja,mp(:,1,1,1),ind(:,1,1,1),iedge)
       end do

    case(2)

       jlo = bbox%lo(2)
       jhi = bbox%hi(2)

       allocate(num_grids_for_j(jlo:jhi))
       allocate(new_iarray_ij(jlo:jhi,ngrids))
       do j = jlo,jhi
          n = 0
          do igrid = 1,ngrids
             b = get_box(ss,igrid)
             if (j >= b%lo(2) .and. j <= b%hi(2)) n = n+1
          end do

          allocate(barray(n))
          allocate(iarray(n))
          allocate(tmp_iarray(n))

          n = 0
          do igrid = 1,ngrids
             b = get_box(ss,igrid)
             if (j >= b%lo(2) .and. j <= b%hi(2)) then
                n = n+1
                barray(n) = get_box(ss,igrid)
                iarray(n) = igrid
             end if
          end do
          num_grids_for_j(j) = n

          call heapsort_indirect_box(barray,tmp_iarray)

          do igrid = 1, num_grids_for_j(j)
             new_iarray_ij(j,igrid) = iarray(tmp_iarray(igrid))
          end do

          deallocate(barray)
          deallocate(iarray)
          deallocate(tmp_iarray)

       end do

       iedge = 1
       inode = 1

       do j = jlo,jhi
          do igrid = 1, num_grids_for_j(j)
             ig = new_iarray_ij(j,igrid)
             sbx = get_box(ss,ig)
             ibx = sbx
             call set_lwb(ibx,2,j)
             call set_upb(ibx,2,j)

             sp => dataptr(ss, ig, ibx)
             mp => dataptr(mm     , ig, ibx)
             ind => dataptr(spo%index_into_aa, ig, ibx)

             call create_aa_2d(sp(:,:,1,:), & 
                  mp(:,:,1,1),              &
                  temp_aa,                  &
                  spo%smt%ia,     & 
                  ind(:,:,1,1),             &
                  inode,iedge)

          end do
       end do

       if ( verbose > 0 ) then
          print *,'FINAL INODE IEDGE ',inode,iedge
       end if

       num_aa = iedge - 1
       allocate(spo%smt%aa(num_aa))
       do i = 1, num_aa
          spo%smt%aa(i) = temp_aa(i)
       end do

       if ( verbose > 0 ) then
          print *,'REAL NUMPTS NUM_AA ',numpts, num_aa
       end if

       spo%smt%ia(inode) = iedge

       call imultifab_fill_boundary(spo%index_into_aa)

       allocate(spo%smt%ja(num_aa))

       iedge = 1
       do j = jlo,jhi
          do igrid = 1, num_grids_for_j(j)
             ig = new_iarray_ij(j,igrid)

             sbx = get_pbox(spo%index_into_aa,ig)
             call set_lwb(sbx,2,j-2)
             call set_upb(sbx,2,j+2)
             ind => dataptr(spo%index_into_aa, ig, sbx)

             mbx = get_box(mm,ig)
             call set_lwb(mbx,2,j)
             call set_upb(mbx,2,j)
             mp => dataptr(mm,ig,mbx)

             call create_ja_2d(spo%smt%ja,mp(:,:,1,1),ind(:,:,1,1),iedge)
          end do
       end do

    case(3)

       jlo = bbox%lo(2)
       jhi = bbox%hi(2)
       klo = bbox%lo(3)
       khi = bbox%hi(3)

       allocate(num_grids_for_jk(jlo:jhi,klo:khi))
       allocate(new_iarray_ijk(jlo:jhi,klo:khi,ngrids))

       do k = klo,khi
          do j = jlo,jhi
             n = 0
             do igrid = 1,ngrids
                b = get_box(ss,igrid)
                if (j >= b%lo(2) .and. j <= b%hi(2) .and. &
                     k >= b%lo(3) .and. k <= b%hi(3)) n = n+1
             end do

             allocate(barray(n))
             allocate(iarray(n))
             allocate(tmp_iarray(n))

             n = 0
             do igrid = 1,ngrids
                b = get_box(ss,igrid)
                if (j >= b%lo(2) .and. j <= b%hi(2) .and. &
                     k >= b%lo(3) .and. k <= b%hi(3)) then
                   n = n+1
                   barray(n) = get_box(ss,igrid)
                   iarray(n) = igrid
                end if
             end do
             num_grids_for_jk(j,k) = n

             call heapsort_indirect_box(barray,tmp_iarray)

             do igrid = 1,num_grids_for_jk(j,k)
                new_iarray_ijk(j,k,igrid) = iarray(tmp_iarray(igrid))
             end do

             deallocate(barray)
             deallocate(iarray)
             deallocate(tmp_iarray)

          end do
       end do

       iedge = 1
       inode = 1

       do k = klo,khi
          do j = jlo,jhi
             do igrid = 1, num_grids_for_jk(j,k)
                ig = new_iarray_ijk(j,k,igrid)
                sbx = get_box(ss,ig)
                ibx = sbx
                call set_lwb(ibx,2,j)
                call set_upb(ibx,2,j)
                call set_lwb(ibx,3,k)
                call set_upb(ibx,3,k)
                sp => dataptr(ss, ig, ibx)
                ind => dataptr(spo%index_into_aa, ig, ibx)
                mp => dataptr(mm,ig,ibx)
                call create_aa_3d(sp(:,:,:,:), &
                     mp(:,:,:,1),              &
                     temp_aa,                  &
                     spo%smt%ia,     & 
                     ind(:,:,:,1),             &
                     inode,iedge)
             end do
          end do
       end do
       spo%smt%ia(inode) = iedge
       if ( verbose > 0 ) then
          print *,'FINAL INODE IEDGE ',inode,iedge
       end if

       num_aa = iedge - 1
       allocate(spo%smt%aa(num_aa))
       do i = 1, num_aa
          spo%smt%aa(i) = temp_aa(i)
       end do

       if ( verbose > 0 ) then
          print *,'REAL NUMPTS NUM_AA ',numpts, num_aa
       end if

       call imultifab_fill_boundary(spo%index_into_aa)

       allocate(spo%smt%ja(num_aa))

       iedge = 1
       inode = 1

       do k = klo,khi
          do j = jlo,jhi
             do igrid = 1, num_grids_for_jk(j,k)
                ig = new_iarray_ijk(j,k,igrid)

                sbx = get_pbox(spo%index_into_aa,ig)
                call set_lwb(sbx,2,j-2)
                call set_upb(sbx,2,j+2)
                call set_lwb(sbx,3,k-2)
                call set_upb(sbx,3,k+2)
                ind => dataptr(spo%index_into_aa, ig, sbx)

                mbx = get_box(mm,ig)
                call set_lwb(mbx,2,j)
                call set_upb(mbx,2,j)
                call set_lwb(mbx,3,k)
                call set_upb(mbx,3,k)
                mp => dataptr(mm,ig,mbx)

                call create_ja_3d(spo%smt%ja,mp(:,:,:,1),ind(:,:,:,1),iedge)
             end do
          end do
       end do
    end select

    ! complete the definition of smt
    spo%smt%numa = num_aa
    spo%smt%nrow = numpts

    call ilut_build(spo%smt, spo%sil)

  contains

    function get_max_neighbors_1d(rh) result(r)
      integer :: r
      type(multifab), intent(in) :: rh
      type(box) :: ibx,jbx
      integer igrid,jgrid,ngrids
      integer i_lbnd
      integer i_ubnd
      integer j_lbnd
      integer j_ubnd
      integer ir
      ngrids = rh%nboxes

      r = 0
      do igrid = 1,ngrids
         ibx = get_box(rh, igrid)
         i_lbnd = ibx%lo(1)
         i_ubnd = ibx%hi(1)

         ir = 0
         do jgrid = 1, ngrids
            if (jgrid == igrid) cycle
            jbx = get_box(rh, jgrid)
            j_lbnd = jbx%lo(1)
            j_ubnd = jbx%hi(1)
            if (j_ubnd == i_lbnd-1) then
               ir = ir + 1
            end if
            if (j_lbnd == i_ubnd+1) then
               ir = ir + 1
            end if
         end do
         r = max(r, ir)
      end do

    end function get_max_neighbors_1d

    function get_max_neighbors_2d(rh) result(r)
      integer :: r
      type(multifab), intent(in) :: rh
      type(box) :: ibx,jbx
      integer igrid,jgrid,ngrids
      integer i_lbnd(rh%dim)
      integer i_ubnd(rh%dim)
      integer j_lbnd(rh%dim)
      integer j_ubnd(rh%dim)
      integer ir
      ngrids = rh%nboxes

      r = 0
      do igrid = 1, ngrids
         ibx = get_box(rh, igrid)
         i_lbnd = ibx%lo(1:rh%dim)
         i_ubnd = ibx%hi(1:rh%dim)

         ir = 0
         do jgrid = 1,ngrids
            if (jgrid == igrid) cycle
            jbx = get_box(rh, jgrid)
            j_lbnd = jbx%lo(1:rh%dim)
            j_ubnd = jbx%hi(1:rh%dim)
            if ( j_ubnd(1) == i_lbnd(1)-1 .and. &
                 j_ubnd(2) >= i_lbnd(2)   .and. &
                 j_lbnd(2) <= i_ubnd(2)   ) then
               ir = ir + 1
            end if
            if ( j_lbnd(1) == i_ubnd(1)+1 .and. &
                 j_ubnd(2) >= i_lbnd(2)   .and. &
                 j_lbnd(2) <= i_ubnd(2)   ) then
               ir = ir + 1
            end if
            if ( j_ubnd(2) == i_lbnd(2)-1 .and. &
                 j_ubnd(1) >= i_lbnd(1)   .and. &
                 j_lbnd(1) <= i_ubnd(1)   ) then
               ir = ir + 1
            end if
            if ( j_lbnd(2) == i_ubnd(2)+1 .and. &
                 j_ubnd(1) >= i_lbnd(1)   .and. &
                 j_lbnd(1) <= i_ubnd(1)   ) then
               ir = ir + 1
            end if
         end do
         r = max(r, ir)
      end do

    end function get_max_neighbors_2d

    function get_max_neighbors_3d(rh) result(r)
      integer :: r
      type(multifab), intent(in) :: rh
      type(box) :: ibx,jbx
      integer igrid,jgrid,ngrids
      integer i_lbnd(rh%dim)
      integer i_ubnd(rh%dim)
      integer j_lbnd(rh%dim)
      integer j_ubnd(rh%dim)
      integer ir
      ngrids = rh%nboxes

      r = 0
      do igrid = 1,ngrids
         ibx = get_box(rh, igrid)
         i_lbnd = ibx%lo(1:rh%dim)
         i_ubnd = ibx%hi(1:rh%dim)

         ir = 0
         do jgrid = 1,ngrids
            if (jgrid == igrid) cycle
            jbx = get_box(rh, jgrid)
            j_lbnd = jbx%lo(1:rh%dim)
            j_ubnd = jbx%hi(1:rh%dim)
            if ( j_ubnd(1) == i_lbnd(1)-1 .and. &
                 j_ubnd(2) >= i_lbnd(2)   .and. &
                 j_lbnd(2) <= i_ubnd(2)   .and. &
                 j_ubnd(3) >= i_lbnd(3)   .and. &
                 j_lbnd(3) <= i_ubnd(3)   ) then
               ir = ir + 1
            end if
            if ( j_lbnd(1) == i_ubnd(1)+1 .and. &
                 j_ubnd(2) >= i_lbnd(2)   .and. &
                 j_lbnd(2) <= i_ubnd(2)   .and. &
                 j_ubnd(3) >= i_lbnd(3)   .and. &
                 j_lbnd(3) <= i_ubnd(3)   ) then
               ir = ir + 1
            end if
            if ( j_ubnd(2) == i_lbnd(2)-1 .and. &
                 j_ubnd(1) >= i_lbnd(1)   .and. &
                 j_lbnd(1) <= i_ubnd(1)   .and. &
                 j_ubnd(3) >= i_lbnd(3)   .and. &
                 j_lbnd(3) <= i_ubnd(3)   ) then
               ir = ir + 1
            end if
            if ( j_lbnd(2) == i_ubnd(2)+1 .and. &
                 j_ubnd(1) >= i_lbnd(1)   .and. &
                 j_lbnd(1) <= i_ubnd(1)   .and. &
                 j_ubnd(3) >= i_lbnd(3)   .and. &
                 j_lbnd(3) <= i_ubnd(3)   ) then
               ir = ir + 1
            end if
            if ( j_ubnd(3) == i_lbnd(3)-1 .and. &
                 j_ubnd(1) >= i_lbnd(1)   .and. &
                 j_lbnd(1) <= i_ubnd(1)   .and. &
                 j_ubnd(2) >= i_lbnd(2)   .and. &
                 j_lbnd(2) <= i_ubnd(2)   ) then
               ir = ir + 1
            end if
            if ( j_lbnd(3) == i_ubnd(3)+1 .and. &
                 j_ubnd(1) >= i_lbnd(1)   .and. &
                 j_lbnd(1) <= i_ubnd(1)   .and. &
                 j_ubnd(2) >= i_lbnd(2)   .and. &
                 j_lbnd(2) <= i_ubnd(2)   ) then
               ir = ir + 1
            end if
         end do
         r = max(r, ir)
      end do

    end function get_max_neighbors_3d

    subroutine make_neighbors_1d(rh, neighbors, sides, num_nbors)

      type(multifab), intent(in) :: rh
      type(box) :: ibx,jbx
      integer igrid,jgrid,ngrids
      integer i_lbnd
      integer i_ubnd
      integer j_lbnd
      integer j_ubnd

      integer :: neighbors(:,:)
      integer :: sides(:,:)
      integer :: num_nbors(:)

      integer inbr 

      ngrids = rh%nboxes

      do igrid = 1,ngrids
         ibx = get_box(rh, igrid)
         i_lbnd = ibx%lo(1)
         i_ubnd = ibx%hi(1)

         inbr = 1
         do jgrid = 1, ngrids
            if (jgrid == igrid) cycle
            jbx = get_box(rh, jgrid)
            j_lbnd = jbx%lo(1)
            j_ubnd = jbx%hi(1)
            if (j_ubnd == i_lbnd-1) then
               neighbors(igrid,inbr) = jgrid
               sides(igrid,inbr) = -1
               inbr = inbr + 1
            end if
            if (j_lbnd == i_ubnd+1) then
               neighbors(igrid,inbr) = jgrid
               sides(igrid,inbr) = 1
               inbr = inbr + 1
            end if
         end do
         num_nbors(igrid) = inbr - 1
      end do

    end subroutine make_neighbors_1d

    subroutine make_neighbors_2d(rh,neighbors,sides,num_nbors)

      type(multifab), intent(in) :: rh
      type(box) :: ibx,jbx
      integer igrid,jgrid,ngrids
      integer i_lbnd(rh%dim)
      integer i_ubnd(rh%dim)
      integer j_lbnd(rh%dim)
      integer j_ubnd(rh%dim)

      integer :: neighbors(:,:)
      integer :: sides(:,:)
      integer :: num_nbors(:)

      integer inbr 

      ngrids = rh%nboxes

      do igrid = 1,ngrids
         ibx = get_box(rh, igrid)
         i_lbnd = ibx%lo(1:rh%dim)
         i_ubnd = ibx%hi(1:rh%dim)

         inbr = 1
         do jgrid = 1,ngrids
            if (jgrid == igrid) cycle
            jbx = get_box(rh, jgrid)
            j_lbnd = jbx%lo(1:rh%dim)
            j_ubnd = jbx%hi(1:rh%dim)
            if ( j_ubnd(1) == i_lbnd(1)-1 .and. &
                 j_ubnd(2) >= i_lbnd(2)   .and. &
                 j_lbnd(2) <= i_ubnd(2)   ) then
               neighbors(igrid,inbr) = jgrid
               sides(igrid,inbr) = -1
               inbr = inbr + 1
            end if
            if ( j_lbnd(1) == i_ubnd(1)+1 .and. &
                 j_ubnd(2) >= i_lbnd(2)   .and. &
                 j_lbnd(2) <= i_ubnd(2)   ) then
               neighbors(igrid,inbr) = jgrid
               sides(igrid,inbr) = 1
               inbr = inbr + 1
            end if
            if ( j_ubnd(2) == i_lbnd(2)-1 .and. &
                 j_ubnd(1) >= i_lbnd(1)   .and. &
                 j_lbnd(1) <= i_ubnd(1)   ) then
               neighbors(igrid,inbr) = jgrid
               sides(igrid,inbr) = -2
               inbr = inbr + 1
            end if
            if ( j_lbnd(2) == i_ubnd(2)+1 .and. &
                 j_ubnd(1) >= i_lbnd(1)   .and. &
                 j_lbnd(1) <= i_ubnd(1)   ) then
               neighbors(igrid,inbr) = jgrid
               sides(igrid,inbr) = 2
               inbr = inbr + 1
            end if
         end do
         num_nbors(igrid) = inbr - 1
      end do

    end subroutine make_neighbors_2d

    subroutine make_neighbors_3d(rh,neighbors,sides,num_nbors)

      type(multifab), intent(in) :: rh
      type(box) :: ibx,jbx
      integer igrid,jgrid,ngrids
      integer i_lbnd(rh%dim)
      integer i_ubnd(rh%dim)
      integer j_lbnd(rh%dim)
      integer j_ubnd(rh%dim)

      integer :: neighbors(:,:)
      integer :: sides(:,:)
      integer :: num_nbors(:)

      integer inbr 

      ngrids = rh%nboxes

      do igrid = 1,ngrids
         ibx = get_box(rh, igrid)
         i_lbnd = ibx%lo(1:rh%dim)
         i_ubnd = ibx%hi(1:rh%dim)

         inbr = 1
         do jgrid = 1,ngrids
            if (jgrid == igrid) cycle
            jbx = get_box(rh, jgrid)
            j_lbnd = jbx%lo(1:rh%dim)
            j_ubnd = jbx%hi(1:rh%dim)
            if ( j_ubnd(1) == i_lbnd(1)-1 .and. &
                 j_ubnd(2) >= i_lbnd(2)   .and. &
                 j_lbnd(2) <= i_ubnd(2)   .and. &
                 j_ubnd(3) >= i_lbnd(3)   .and. &
                 j_lbnd(3) <= i_ubnd(3)   ) then
               neighbors(igrid,inbr) = jgrid
               sides(igrid,inbr) = -1
               inbr = inbr + 1
            end if
            if ( j_lbnd(1) == i_ubnd(1)+1 .and. &
                 j_ubnd(2) >= i_lbnd(2)   .and. &
                 j_lbnd(2) <= i_ubnd(2)   .and. &
                 j_ubnd(3) >= i_lbnd(3)   .and. &
                 j_lbnd(3) <= i_ubnd(3)   ) then
               neighbors(igrid,inbr) = jgrid
               sides(igrid,inbr) = 1
               inbr = inbr + 1
            end if
            if ( j_ubnd(2) == i_lbnd(2)-1 .and. &
                 j_ubnd(1) >= i_lbnd(1)   .and. &
                 j_lbnd(1) <= i_ubnd(1)   .and. &
                 j_ubnd(3) >= i_lbnd(3)   .and. &
                 j_lbnd(3) <= i_ubnd(3)   ) then
               neighbors(igrid,inbr) = jgrid
               sides(igrid,inbr) = -2
               inbr = inbr + 1
            end if
            if ( j_lbnd(2) == i_ubnd(2)+1 .and. &
                 j_ubnd(1) >= i_lbnd(1)   .and. &
                 j_lbnd(1) <= i_ubnd(1)   .and. &
                 j_ubnd(3) >= i_lbnd(3)   .and. &
                 j_lbnd(3) <= i_ubnd(3)   ) then
               neighbors(igrid,inbr) = jgrid
               sides(igrid,inbr) = 2
               inbr = inbr + 1
            end if
            if ( j_ubnd(3) == i_lbnd(3)-1 .and. &
                 j_ubnd(1) >= i_lbnd(1)   .and. &
                 j_lbnd(1) <= i_ubnd(1)   .and. &
                 j_ubnd(2) >= i_lbnd(2)   .and. &
                 j_lbnd(2) <= i_ubnd(2)   ) then
               neighbors(igrid,inbr) = jgrid
               sides(igrid,inbr) = -3
               inbr = inbr + 1
            end if
            if ( j_lbnd(3) == i_ubnd(3)+1 .and. &
                 j_ubnd(1) >= i_lbnd(1)   .and. &
                 j_lbnd(1) <= i_ubnd(1)   .and. &
                 j_ubnd(2) >= i_lbnd(2)   .and. &
                 j_lbnd(2) <= i_ubnd(2)   ) then
               neighbors(igrid,inbr) = jgrid
               sides(igrid,inbr) = 3
               inbr = inbr + 1
            end if
         end do
         num_nbors(igrid) = inbr - 1
      end do

    end subroutine make_neighbors_3d

    subroutine create_aa_1d(sp, mp, aa, ia, ind, inode, iedge)
      integer           , intent(in) :: mp(:)
      real (kind = dp_t), intent(in) :: sp(:,0:)
      real (kind = dp_t), intent(inout) :: aa(:)
      integer, intent(inout) :: ia(:)
      integer, intent(inout) :: ind(-1:)
      integer, intent(inout) :: iedge,inode
      integer nx
      integer i
      integer, parameter :: XBC = 3

      nx = size(sp,dim=1)

      i = 1
      ia(inode) = iedge
      ind(i) = inode

      if (bc_interior(mp(i),1,-1)) then
         aa(iedge) = sp(i,2) ; iedge = iedge + 1           ! Left
      end if

      aa(iedge)  = sp(i,0) ; iedge = iedge + 1            ! Center

      if (bc_interior(mp(i),1,+1)) then
         aa(iedge) = sp(i,1) ; iedge = iedge + 1         ! Right by 1
      end if

      if (bc_skewed(mp(i),1,+1)) then
         aa(iedge) = sp(i,XBC) ; iedge = iedge + 1         ! Right by 2
      end if

      inode = inode + 1

      do i = 2, nx-1
         ia(inode) = iedge
         ind(i) = inode

         aa(iedge) = sp(i,2) ; iedge = iedge + 1           ! Left
         aa(iedge) = sp(i,0) ; iedge = iedge + 1           ! Center
         aa(iedge) = sp(i,1) ; iedge = iedge + 1           ! Right

         inode = inode + 1
      end do

      if (nx /= 1) then
         i = nx
         ia(inode) = iedge
         ind(i) = inode

         if (bc_skewed(mp(i),1,-1)) then
            aa(iedge) = sp(i,XBC) ; iedge = iedge + 1         ! Left by 2
            aa(iedge) = sp(i,2) ; iedge = iedge + 1         ! Left
         else 
            aa(iedge) = sp(i,2) ; iedge = iedge + 1         ! Left
         end if

         aa(iedge) = sp(i,0) ; iedge = iedge + 1           ! Center

         if (bc_interior(mp(i),1,+1)) then
            aa(iedge) = sp(i,1) ; iedge = iedge + 1         ! Right
         end if

         inode = inode + 1
      end if

    end subroutine create_aa_1d

    subroutine create_ja_1d(ja, mp, ind, iedge)

      integer, intent(inout) :: ja(:)
      integer, intent(in   ) :: mp(:)
      integer, intent(inout) :: ind(-1:)
      integer iedge

      integer nx
      integer i

      nx = size(mp,dim=1)

      i = 1
      if (bc_interior(mp(i),1,-1)) then
         ja(iedge) = ind(i-1) ; iedge = iedge + 1          ! Left
      end if

      ja(iedge) = ind(i  ) ; iedge = iedge + 1          ! Center

      if (bc_interior(mp(i),1,+1)) then
         ja(iedge) = ind(i+1) ; iedge = iedge + 1          ! Right
      end if

      if (bc_skewed(mp(i),1,+1)) then
         ja(iedge) = ind(i+2) ; iedge = iedge + 1          ! Right by 2
      end if

      do i = 2, nx-1
         ja(iedge) = ind(i-1) ; iedge = iedge + 1          ! Left
         ja(iedge) = ind(i  ) ; iedge = iedge + 1          ! Center
         ja(iedge) = ind(i+1) ; iedge = iedge + 1          ! Right
      end do

      if (nx /= 1) then
         i = nx
         if (bc_skewed(mp(i),1,-1)) then
            ja(iedge) = ind(i-2) ; iedge = iedge + 1        ! Left by 2
         end if

         ja(iedge) = ind(i-1) ; iedge = iedge + 1        ! Left
         ja(iedge) = ind(i  ) ; iedge = iedge + 1        ! Center

         if (bc_interior(mp(i),1,+1)) then
            ja(iedge) = ind(i+1) ; iedge = iedge + 1        ! Right
         end if
      end if


    end subroutine create_ja_1d

    subroutine create_aa_2d(sp, mp, aa, ia, ind, inode, iedge)

      integer           , intent(in)    :: mp(:,:)
      real (kind = dp_t), intent(in)    :: sp(:,:,0:)
      real (kind = dp_t), intent(inout) :: aa(:)
      integer, intent(inout) :: ind(:,:)
      integer, intent(inout) :: ia(:)
      integer nx,ny
      integer i,j
      integer iedge,inode
      integer, parameter :: XBC = 5, YBC = 6

      nx = size(sp,dim=1)
      ny = size(sp,dim=2)

      j = 1

      !   **************************************************************************

      i = 1
      ia(inode) = iedge
      ind(i,j) = inode

      if (bc_skewed(mp(i,j),2,-1)) then
         aa(iedge) = sp(i,j,YBC) ; iedge = iedge + 1       ! Down by 2
         aa(iedge) = sp(i,j,4) ; iedge = iedge + 1       ! Down
      else if (bc_interior(mp(i,j),2,-1)) then
         aa(iedge) = sp(i,j,4) ; iedge = iedge + 1       ! Down 
      end if

      if (bc_interior(mp(i,j),1,-1)) then
         aa(iedge) = sp(i,j,2) ; iedge = iedge + 1       ! Left
      end if

      aa(iedge)  = sp(i,j,0) ; iedge = iedge + 1        ! Center

      if (bc_interior(mp(i,j),1,+1)) then
         aa(iedge) = sp(i,j,1) ; iedge = iedge + 1     ! Right
      end if

      if (bc_skewed(mp(i,j),1,+1)) then
         aa(iedge) = sp(i,j,XBC) ; iedge = iedge + 1     ! Right by 2
      end if

      if (bc_interior(mp(i,j),2,+1)) then
         aa(iedge) = sp(i,j,3) ; iedge = iedge + 1       ! Up
      end if
      if (bc_skewed(mp(i,j),2,+1)) then
         aa(iedge) = sp(i,j,YBC) ; iedge = iedge + 1       ! Up by 2
      end if

      inode = inode + 1

      !   **************************************************************************

      do i = 2, nx-1
         ia(inode) = iedge
         ind(i,j) = inode

         if (bc_skewed(mp(i,j),2,-1)) then
            aa(iedge) = sp(i,j,YBC) ; iedge = iedge + 1     ! Down by 2
            aa(iedge) = sp(i,j,4) ; iedge = iedge + 1     ! Down
         else if (bc_interior(mp(i,j),2,-1)) then
            aa(iedge) = sp(i,j,4) ; iedge = iedge + 1     ! Down 
         end if

         aa(iedge) = sp(i,j,2) ; iedge = iedge + 1       ! Left
         aa(iedge) = sp(i,j,0) ; iedge = iedge + 1       ! Center
         aa(iedge) = sp(i,j,1) ; iedge = iedge + 1       ! Right

         if (bc_interior(mp(i,j),2,+1)) then
            aa(iedge) = sp(i,j,3) ; iedge = iedge + 1     ! Up
         end if
         if (bc_skewed(mp(i,j),2,+1)) then
            aa(iedge) = sp(i,j,YBC) ; iedge = iedge + 1     ! Up by 2
         end if

         inode = inode + 1
      end do

      !   **************************************************************************

      if (nx /= 1) then
         i = nx
         ia(inode) = iedge
         ind(i,j) = inode

         if (bc_skewed(mp(i,j),2,-1)) then
            aa(iedge) = sp(i,j,YBC) ; iedge = iedge + 1     ! Down by 2
            aa(iedge) = sp(i,j,4) ; iedge = iedge + 1     ! Down
         else if (bc_interior(mp(i,j),2,-1)) then
            aa(iedge) = sp(i,j,4) ; iedge = iedge + 1     ! Down 
         end if

         if (bc_skewed(mp(i,j),1,-1)) then
            aa(iedge) = sp(i,j,XBC) ; iedge = iedge + 1     ! Left by 2
            aa(iedge) = sp(i,j,2) ; iedge = iedge + 1     ! Left
         else
            aa(iedge) = sp(i,j,2) ; iedge = iedge + 1     ! Left
         end if

         aa(iedge) = sp(i,j,0) ; iedge = iedge + 1       ! Center

         if (bc_interior(mp(i,j),1,+1)) then
            aa(iedge) = sp(i,j,1) ; iedge = iedge + 1     ! Right
         end if

         if (bc_interior(mp(i,j),2,+1)) then
            aa(iedge) = sp(i,j,3) ; iedge = iedge + 1     ! Up
         end if
         if (bc_skewed(mp(i,j),2,+1)) then
            aa(iedge) = sp(i,j,YBC) ; iedge = iedge + 1     ! Up by 2
         end if

         inode = inode + 1
      end if

    end subroutine create_aa_2d

    subroutine create_ja_2d(ja, mp, ind, iedge)

      integer, intent(inout) :: ja(:)
      integer, intent(in   ) :: mp(:,:)
      integer, intent(inout) :: ind(-1:,-1:)
      integer iedge

      integer nx
      integer i,j

      nx = size(mp,dim=1)

      j = 1

      !**************************************************************************

      i = 1

      if (bc_skewed(mp(i,j),2,-1)) then
         ja(iedge) = ind(i,j-2) ; iedge = iedge + 1      ! Down by 2
      end if

      if (bc_interior(mp(i,j),2,-1)) then
         ja(iedge) = ind(i,j-1) ; iedge = iedge + 1       ! Down 
      end if

      if (bc_interior(mp(i,j),1,-1)) then
         ja(iedge) = ind(i-1,j) ; iedge = iedge + 1       ! Left
      end if

      ja(iedge)  = ind(i,j) ; iedge = iedge + 1          ! Center

      if (bc_interior(mp(i,j),1,+1)) then
         ja(iedge) = ind(i+1,j) ; iedge = iedge + 1       ! Right
      end if

      if (bc_skewed(mp(i,j),1,+1)) then
         ja(iedge) = ind(i+2,j) ; iedge = iedge + 1       ! Right by 2
      end if

      if (bc_interior(mp(i,j),2,+1)) then
         ja(iedge) = ind(i,j+1) ; iedge = iedge + 1       ! Up
      end if
      if (bc_skewed(mp(i,j),2,+1)) then
         ja(iedge) = ind(i,j+2) ; iedge = iedge + 1       ! Up by 2
      end if

      !**************************************************************************

      do i = 2, nx-1

         if (bc_skewed(mp(i,j),2,-1)) then
            ja(iedge) = ind(i,j-2) ; iedge = iedge + 1     ! Down by 2
         end if
         if (bc_interior(mp(i,j),2,-1)) then
            ja(iedge) = ind(i,j-1) ; iedge = iedge + 1     ! Down 
         end if

         ja(iedge) = ind(i-1,j) ; iedge = iedge + 1       ! Left
         ja(iedge) = ind(i  ,j) ; iedge = iedge + 1       ! Center
         ja(iedge) = ind(i+1,j) ; iedge = iedge + 1       ! Right

         if (bc_interior(mp(i,j),2,+1)) then
            ja(iedge) = ind(i,j+1) ; iedge = iedge + 1     ! Up
         end if
         if (bc_skewed(mp(i,j),2,+1)) then
            ja(iedge) = ind(i,j+2) ; iedge = iedge + 1     ! Up by 2
         end if

      end do

      !**************************************************************************

      if (nx /= 1) then
         i = nx

         if (bc_skewed(mp(i,j),2,-1)) then
            ja(iedge) = ind(i,j-2) ; iedge = iedge + 1     ! Down by 2
         end if
         if (bc_interior(mp(i,j),2,-1)) then
            ja(iedge) = ind(i,j-1) ; iedge = iedge + 1     ! Down 
         end if

         if (bc_skewed(mp(i,j),1,-1)) then
            ja(iedge) = ind(i-2,j) ; iedge = iedge + 1     ! Left by 2
         end if

         ja(iedge) = ind(i-1,j) ; iedge = iedge + 1       ! Left
         ja(iedge) = ind(i  ,j) ; iedge = iedge + 1       ! Center

         if (bc_interior(mp(i,j),1,+1)) then
            ja(iedge) = ind(i+1,j) ; iedge = iedge + 1     ! Right
         end if

         if (bc_interior(mp(i,j),2,+1)) then
            ja(iedge) = ind(i,j+1) ; iedge = iedge + 1       ! Up
         end if
         if (bc_skewed(mp(i,j),2,+1)) then
            ja(iedge) = ind(i,j+2) ; iedge = iedge + 1       ! Up by 2
         end if

      end if

    end subroutine create_ja_2d

    subroutine create_aa_3d(sp, mp, aa, ia, ind, inode, iedge)
      integer           , intent(in) :: mp(:,:,:)
      real (kind = dp_t), intent(in) :: sp(:,:,:,0:)
      real (kind = dp_t), intent(inout) :: aa(:)
      integer, intent(inout) :: ia(:)
      integer, intent(inout) :: ind(:,:,:)
      integer iedge,inode

      integer nx
      integer i,j,k
      integer, parameter :: XBC = 7, YBC = 8, ZBC = 9

      nx = size(sp,dim=1)

      k = 1
      j = 1

      !**************************************************************************

      i = 1
      ia(inode) = iedge
      ind(i,j,k) = inode

      if (bc_skewed(mp(i,j,k),3,-1)) then
         aa(iedge) = sp(i,j,k,ZBC) ; iedge = iedge + 1       ! Back by 2
         aa(iedge) = sp(i,j,k,6) ; iedge = iedge + 1       ! Back
      else if (bc_interior(mp(i,j,k),3,-1)) then
         aa(iedge) = sp(i,j,k,6) ; iedge = iedge + 1       ! Back 
      end if

      if (bc_skewed(mp(i,j,k),2,-1)) then
         aa(iedge) = sp(i,j,k,YBC) ; iedge = iedge + 1       ! Down by 2
         aa(iedge) = sp(i,j,k,4) ; iedge = iedge + 1       ! Down
      else if (bc_interior(mp(i,j,k),2,-1)) then
         aa(iedge) = sp(i,j,k,4) ; iedge = iedge + 1       ! Down 
      end if

      if (bc_interior(mp(i,j,k),1,-1)) then
         aa(iedge) = sp(i,j,k,2) ; iedge = iedge + 1       ! Left
      end if

      aa(iedge)  = sp(i,j,k,0) ; iedge = iedge + 1        ! Center

      if (bc_interior(mp(i,j,k),1,+1)) then
         aa(iedge) = sp(i,j,k,1) ; iedge = iedge + 1     ! Right
      end if

      if (bc_skewed(mp(i,j,k),1,+1)) then
         aa(iedge) = sp(i,j,k,XBC) ; iedge = iedge + 1     ! Right by 2
      end if

      if (bc_interior(mp(i,j,k),2,+1)) then
         aa(iedge) = sp(i,j,k,3) ; iedge = iedge + 1       ! Up
      end if
      if (bc_skewed(mp(i,j,k),2,+1)) then
         aa(iedge) = sp(i,j,k,YBC) ; iedge = iedge + 1       ! Up by 2
      end if

      if (bc_interior(mp(i,j,k),3,+1)) then
         aa(iedge) = sp(i,j,k,5) ; iedge = iedge + 1       ! Fwd
      end if
      if (bc_skewed(mp(i,j,k),3,+1)) then
         aa(iedge) = sp(i,j,k,ZBC) ; iedge = iedge + 1       ! Fwd by 2
      end if

      inode = inode + 1

      !**************************************************************************

      do i = 2, nx-1
         ia(inode) = iedge
         ind(i,j,k) = inode

         if (bc_skewed(mp(i,j,k),3,-1)) then
            aa(iedge) = sp(i,j,k,ZBC) ; iedge = iedge + 1     ! Back by 2
            aa(iedge) = sp(i,j,k,6) ; iedge = iedge + 1     ! Back
         else if (bc_interior(mp(i,j,k),3,-1)) then
            aa(iedge) = sp(i,j,k,6) ; iedge = iedge + 1     ! Back 
         end if

         if (bc_skewed(mp(i,j,k),2,-1)) then
            aa(iedge) = sp(i,j,k,YBC) ; iedge = iedge + 1     ! Down by 2
            aa(iedge) = sp(i,j,k,4) ; iedge = iedge + 1     ! Down
         else if (bc_interior(mp(i,j,k),2,-1)) then
            aa(iedge) = sp(i,j,k,4) ; iedge = iedge + 1     ! Down 
         end if

         aa(iedge) = sp(i,j,k,2) ; iedge = iedge + 1       ! Left
         aa(iedge) = sp(i,j,k,0) ; iedge = iedge + 1       ! Center
         aa(iedge) = sp(i,j,k,1) ; iedge = iedge + 1       ! Right

         if (bc_interior(mp(i,j,k),2,+1)) then
            aa(iedge) = sp(i,j,k,3) ; iedge = iedge + 1     ! Up
         end if
         if (bc_skewed(mp(i,j,k),2,+1)) then
            aa(iedge) = sp(i,j,k,YBC) ; iedge = iedge + 1     ! Up by 2
         end if

         if (bc_interior(mp(i,j,k),3,+1)) then
            aa(iedge) = sp(i,j,k,5) ; iedge = iedge + 1     ! Fwd
         end if
         if (bc_skewed(mp(i,j,k),3,+1)) then
            aa(iedge) = sp(i,j,k,ZBC) ; iedge = iedge + 1     ! Fwd by 2
         end if

         inode = inode + 1
      end do

      !**************************************************************************

      if (nx /= 1) then
         i = nx
         ia(inode) = iedge
         ind(i,j,k) = inode

         if (bc_skewed(mp(i,j,k),3,-1)) then
            aa(iedge) = sp(i,j,k,ZBC) ; iedge = iedge + 1     ! Back by 2
            aa(iedge) = sp(i,j,k,6) ; iedge = iedge + 1     ! Back
         else if (bc_interior(mp(i,j,k),3,-1)) then
            aa(iedge) = sp(i,j,k,6) ; iedge = iedge + 1     ! Back 
         end if

         if (bc_skewed(mp(i,j,k),2,-1)) then
            aa(iedge) = sp(i,j,k,YBC) ; iedge = iedge + 1     ! Down by 2
            aa(iedge) = sp(i,j,k,4) ; iedge = iedge + 1     ! Down
         else if (bc_interior(mp(i,j,k),2,-1)) then
            aa(iedge) = sp(i,j,k,4) ; iedge = iedge + 1     ! Down 
         end if

         if (bc_skewed(mp(i,j,k),1,-1)) then
            aa(iedge) = sp(i,j,k,XBC) ; iedge = iedge + 1     ! Left by 2
            aa(iedge) = sp(i,j,k,2) ; iedge = iedge + 1     ! Left
         else
            aa(iedge) = sp(i,j,k,2) ; iedge = iedge + 1     ! Left
         end if

         aa(iedge) = sp(i,j,k,0) ; iedge = iedge + 1       ! Center

         if (bc_interior(mp(i,j,k),1,+1)) then
            aa(iedge) = sp(i,j,k,1) ; iedge = iedge + 1     ! Right
         end if

         if (bc_interior(mp(i,j,k),2,+1)) then
            aa(iedge) = sp(i,j,k,3) ; iedge = iedge + 1     ! Up
         end if
         if (bc_skewed(mp(i,j,k),2,+1)) then
            aa(iedge) = sp(i,j,k,YBC) ; iedge = iedge + 1     ! Up by 2
         end if

         if (bc_interior(mp(i,j,k),3,+1)) then
            aa(iedge) = sp(i,j,k,5) ; iedge = iedge + 1     ! Fwd
         end if
         if (bc_skewed(mp(i,j,k),3,+1)) then
            aa(iedge) = sp(i,j,k,ZBC) ; iedge = iedge + 1     ! Fwd by 2
         end if

         inode = inode + 1
      end if

    end subroutine create_aa_3d

    subroutine create_ja_3d(ja, mp, ind, iedge)
      integer, intent(inout) :: ja(:)
      integer, intent(in   ) :: mp(:,:,:)
      integer, intent(inout) :: ind(-1:,-1:,-1:)
      integer iedge

      integer nx
      integer i,j,k

      nx = size(mp,dim=1)

      k = 1
      j = 1

      !**************************************************************************

      i = 1

      if (bc_skewed(mp(i,j,k),3,-1)) then
         ja(iedge) = ind(i,j,k-2) ; iedge = iedge + 1     ! Back by 2
      end if

      if (bc_interior(mp(i,j,k),3,-1)) then
         ja(iedge) = ind(i,j,k-1) ; iedge = iedge + 1     ! Back 
      end if

      if (bc_skewed(mp(i,j,k),2,-1)) then
         ja(iedge) = ind(i,j-2,k) ; iedge = iedge + 1     ! Down by 2
      end if

      if (bc_interior(mp(i,j,k),2,-1)) then
         ja(iedge) = ind(i,j-1,k) ; iedge = iedge + 1     ! Down 
      end if

      if (bc_interior(mp(i,j,k),1,-1)) then
         ja(iedge) = ind(i-1,j,k) ; iedge = iedge + 1     ! Left
      end if

      ja(iedge)  = ind(i,j,k) ; iedge = iedge + 1        ! Center

      if (bc_interior(mp(i,j,k),1,+1)) then
         ja(iedge) = ind(i+1,j,k) ; iedge = iedge + 1     ! Right
      end if

      if (bc_skewed(mp(i,j,k),1,+1)) then
         ja(iedge) = ind(i+2,j,k) ; iedge = iedge + 1     ! Right by 2
      end if

      if (bc_interior(mp(i,j,k),2,+1)) then
         ja(iedge) = ind(i,j+1,k) ; iedge = iedge + 1     ! Up
      end if
      if (bc_skewed(mp(i,j,k),2,+1)) then
         ja(iedge) = ind(i,j+2,k) ; iedge = iedge + 1     ! Up by 2
      end if

      if (bc_interior(mp(i,j,k),3,+1)) then
         ja(iedge) = ind(i,j,k+1) ; iedge = iedge + 1     ! Fwd
      end if
      if (bc_skewed(mp(i,j,k),3,+1)) then
         ja(iedge) = ind(i,j,k+2) ; iedge = iedge + 1     ! Fwd by 2
      end if

      !**************************************************************************

      do i = 2, nx-1

         if (bc_skewed(mp(i,j,k),3,-1)) then
            ja(iedge) = ind(i,j,k-2) ; iedge = iedge + 1   ! Back by 2
         end if
         if (bc_interior(mp(i,j,k),3,-1)) then
            ja(iedge) = ind(i,j,k-1) ; iedge = iedge + 1   ! Back 
         end if
         if (bc_skewed(mp(i,j,k),2,-1)) then
            ja(iedge) = ind(i,j-2,k) ; iedge = iedge + 1   ! Down by 2
         end if
         if (bc_interior(mp(i,j,k),2,-1)) then
            ja(iedge) = ind(i,j-1,k) ; iedge = iedge + 1   ! Down 
         end if

         ja(iedge) = ind(i-1,j,k) ; iedge = iedge + 1     ! Left
         ja(iedge) = ind(i  ,j,k) ; iedge = iedge + 1     ! Center
         ja(iedge) = ind(i+1,j,k) ; iedge = iedge + 1     ! Right

         if (bc_interior(mp(i,j,k),2,+1)) then
            ja(iedge) = ind(i,j+1,k) ; iedge = iedge + 1   ! Up
         end if
         if (bc_skewed(mp(i,j,k),2,+1)) then
            ja(iedge) = ind(i,j+2,k) ; iedge = iedge + 1   ! Up by 2
         end if

         if (bc_interior(mp(i,j,k),3,+1)) then
            ja(iedge) = ind(i,j,k+1) ; iedge = iedge + 1   ! Fwd
         end if
         if (bc_skewed(mp(i,j,k),3,+1)) then
            ja(iedge) = ind(i,j,k+2) ; iedge = iedge + 1   ! Fwd by 2
         end if

      end do

      !**************************************************************************

      if (nx /= 1) then
         i = nx

         if (bc_skewed(mp(i,j,k),3,-1)) then
            ja(iedge) = ind(i,j,k-2) ; iedge = iedge + 1   ! Back by 2
         end if
         if (bc_interior(mp(i,j,k),3,-1)) then
            ja(iedge) = ind(i,j,k-1) ; iedge = iedge + 1   ! Back 
         end if

         if (bc_skewed(mp(i,j,k),2,-1)) then
            ja(iedge) = ind(i,j-2,k) ; iedge = iedge + 1   ! Down by 2
         end if
         if (bc_interior(mp(i,j,k),2,-1)) then
            ja(iedge) = ind(i,j-1,k) ; iedge = iedge + 1   ! Down 
         end if

         if (bc_skewed(mp(i,j,k),1,-1)) then
            ja(iedge) = ind(i-2,j,k) ; iedge = iedge + 1   ! Left by 2
         end if

         ja(iedge) = ind(i-1,j,k) ; iedge = iedge + 1     ! Left
         ja(iedge) = ind(i  ,j,k) ; iedge = iedge + 1     ! Center

         if (bc_interior(mp(i,j,k),1,+1)) then
            ja(iedge) = ind(i+1,j,k) ; iedge = iedge + 1   ! Right
         end if

         if (bc_interior(mp(i,j,k),2,+1)) then
            ja(iedge) = ind(i,j+1,k) ; iedge = iedge + 1   ! Up
         end if
         if (bc_skewed(mp(i,j,k),2,+1)) then
            ja(iedge) = ind(i,j+2,k) ; iedge = iedge + 1   ! Up by 2
         end if

         if (bc_interior(mp(i,j,k),3,+1)) then
            ja(iedge) = ind(i,j,k+1) ; iedge = iedge + 1   ! Fwd
         end if
         if (bc_skewed(mp(i,j,k),3,+1)) then
            ja(iedge) = ind(i,j,k+2) ; iedge = iedge + 1   ! Fwd by 2
         end if

      end if

    end subroutine create_ja_3d

  subroutine ilut_build(amat, sil)
    type(sparse_matrix), intent(in) :: amat
    type(sparse_matrix), intent(out) :: sil
    integer :: ierr
    integer, allocatable :: iw(:)
    real(kind=dp_t), allocatable ::  wk(:)
    real(kind=dp_t) :: tol
    integer :: nwk
    integer :: lfil, nrow
    external ilut

    nrow = amat%nrow

    lfil = 4
    tol = 1.e-5_dp_t

    ! FIXME: this is a guess as to the maximum size of these arrays.
    allocate(wk(nrow+1))
    allocate(iw(2*nrow))

    call sparse_matrix_build(sil, (2*lfil+1)*nrow, nrow)

    nwk = size(sil%aa)
    call ilut (nrow, amat%aa, amat%ja, amat%ia, lfil, tol, sil%aa, sil%ja, sil%ia, nwk, wk, iw, ierr)
    if ( ierr /= 0 ) then
       write(*,*) 'ILUT: IERR = ', ierr
       call bl_error("ILUT: failed")
    end if

  end subroutine ilut_build

  end subroutine sparse_build

  subroutine sparse_nodal_build(spo, ss, mm, la, face_type, verbose)
    type(   layout), intent(in) :: la
    type(imultifab), intent(in) :: mm
    type(multifab) , intent(in) :: ss
    type(sparse)                :: spo
    integer, intent(in)         :: face_type(:,:,:)
    integer, intent(in)         :: verbose

    integer        , pointer :: mp(:,:,:,:)
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    real(kind=dp_t), allocatable :: temp_aa(:)
    type(box) :: b,sbx,mbx,ibx,jbx,shifted_box
    type(box) :: bbox,jbox

    type(imultifab) :: mm_grown

    type(box), allocatable :: barray(:)
    integer  , allocatable :: iarray(:)
    integer  , allocatable :: tmp_iarray(:)
    integer  , allocatable :: new_iarray_i(:)
    integer  , allocatable :: new_iarray_ij(:,:)
    integer  , allocatable :: new_iarray_ijk(:,:,:)

    integer :: pts_to_add,numpts, num_aa
    integer :: dm
    integer :: ig,igrid,jgrid,ngrids
    integer :: i,j,k,n,ns,nx,ny,nz
    integer :: dir
    integer :: jlo, jhi, klo, khi
    integer :: ii,jj,kk
    integer :: inode, iedge
    logical :: at_jhi,at_khi

    integer i_lbnd(ss%dim)
    integer i_ubnd(ss%dim)

    integer max_neighbors
    integer, allocatable :: neighbors(:,:)
    integer, allocatable :: sides(:,:)
    integer, allocatable :: num_nbors(:)

    integer, allocatable :: num_grids_for_j(:)
    integer, allocatable :: num_grids_for_jk(:,:)

    integer    , pointer :: ind(:,:,:,:)

    ns = ss%nc
    dm = ss%dim

    ngrids = ss%nboxes

    spo%Anorm = stencil_norm(ss)

!   Build a mask with one ghost cell for use in creating aa
    call imultifab_build(mm_grown,la,1,1,mm%nodal)
    do igrid = 1, ngrids
       mp => dataptr(mm_grown,igrid)
       mp(:,:,:,:) = ibset(mp(:,:,:,:), BC_BIT(BC_DIR,1,0))
    end do
    call copy(mm_grown,mm,all=.false.)
    call imultifab_fill_boundary(mm_grown)

!   Build the ia array.
    call imultifab_build(spo%index_into_aa,la,1,1,ss%nodal)
    call setval(spo%index_into_aa,-1,all=.true.)

    numpts = 0
    do igrid = 1, ngrids
       numpts = numpts + volume(get_ibox(ss,igrid))
    end do
    num_aa = numpts * ns

    if ( verbose > 0 ) then
       print *,'PROJECTED NUMPTS NUM_AA ',numpts, num_aa
    end if

    allocate(temp_aa(num_aa))
    allocate(spo%smt%ia(numpts+1))

    !   Create the minimum box containing all the grids
    bbox = get_ibox(ss, 1)
    do igrid = 2,ngrids
       bbox = box_bbox(bbox,get_ibox(ss,igrid))
    end do

    !   Convert the ss into CSR format.

    select case(dm)
    case(1)
       iedge = 1
       inode = 1

       allocate(barray(ngrids))
       allocate(iarray(ngrids))
       allocate(tmp_iarray(ngrids))
       allocate(new_iarray_i(ngrids))

       do igrid = 1,ngrids
          barray(igrid) = get_box(ss,igrid)
          iarray(igrid) = igrid
       end do
       call heapsort_indirect_box(barray,tmp_iarray)
       do igrid = 1,ngrids
          new_iarray_i(igrid) = iarray(tmp_iarray(igrid))
       end do

       deallocate(barray)
       deallocate(iarray)
       deallocate(tmp_iarray)

       do igrid = 1, ngrids
          ig = new_iarray_i(igrid)
          sp => dataptr(ss, ig)
          mp => dataptr(mm_grown,ig)
          ind => dataptr(spo%index_into_aa, ig)
          call create_nodal_aa_1d(sp(:,1,1,:), &
               mp(:,1,1,1),              &
               temp_aa,                  &
               spo%smt%ia,     & 
               ind(:,1,1,1),             &
               inode, iedge)
       end do
       spo%smt%ia(inode) = iedge

       num_aa = iedge - 1
       numpts = inode - 1
       allocate(spo%smt%aa(num_aa))
       do i = 1, num_aa
          spo%smt%aa(i) = temp_aa(i)
       end do

       if ( verbose > 0 ) then
          print *,'ACTUAL NUMPTS NUM_AA ',numpts, num_aa
       end if

       call imultifab_fill_boundary(spo%index_into_aa)
       call copy_nodal_ind_on_intersect(spo%index_into_aa)

       allocate(spo%smt%ja(num_aa))

       iedge = 1
       do igrid = 1, ngrids
          ig = new_iarray_i(igrid)
          mp => dataptr(mm_grown,ig)
          sbx = get_pbox(spo%index_into_aa,ig)
          ind => dataptr(spo%index_into_aa,ig)
          call create_nodal_ja_1d(spo%smt%ja,mp(:,1,1,1),ind(:,1,1,1),iedge)
       end do

    case(2)

       jlo = bbox%lo(2)
       jhi = bbox%hi(2)

       allocate(num_grids_for_j(jlo:jhi))
       allocate(new_iarray_ij(jlo:jhi,ngrids))
       do j = jlo,jhi
          n = 0
          do igrid = 1,ngrids
             b = get_ibox(ss,igrid)
             if (j >= b%lo(2) .and. j <= b%hi(2)) n = n+1
          end do

          allocate(barray(n))
          allocate(iarray(n))
          allocate(tmp_iarray(n))

          n = 0
          do igrid = 1,ngrids
             b = get_ibox(ss,igrid)
             if (j >= b%lo(2) .and. j <= b%hi(2)) then
                n = n+1
                barray(n) = get_ibox(ss,igrid)
                iarray(n) = igrid
             end if
          end do
          num_grids_for_j(j) = n

          call heapsort_indirect_box(barray,tmp_iarray)

          do igrid = 1, num_grids_for_j(j)
             new_iarray_ij(j,igrid) = iarray(tmp_iarray(igrid))
          end do

          deallocate(barray)
          deallocate(iarray)
          deallocate(tmp_iarray)

       end do

       iedge = 1
       inode = 1

       do j = jlo,jhi
          do igrid = 1, num_grids_for_j(j)
             ig = new_iarray_ij(j,igrid)

             ibx = get_ibox(ss,ig)
             at_jhi = (j == box_upb_d(ibx,2))

             call set_lwb(ibx,2,j)
             call set_upb(ibx,2,j)
             sp  => dataptr(ss               , ig, ibx)
             ind => dataptr(spo%index_into_aa, ig, ibx)

             ibx = get_pbox(mm_grown,ig)
             call set_lwb(ibx,2,j-1)
             call set_upb(ibx,2,j+1)
             mp => dataptr(mm_grown, ig, ibx)

             call create_nodal_aa_2d(sp(:,:,1,:), & 
                  mp(:,:,1,1),              &
                  temp_aa,                  &
                  spo%smt%ia,     & 
                  ind(:,:,1,1),             &
                  inode,iedge,at_jhi)

          end do
       end do

       num_aa = iedge - 1
       numpts = inode - 1
       allocate(spo%smt%aa(num_aa))
       do i = 1, num_aa
          spo%smt%aa(i) = temp_aa(i)
       end do

       if ( verbose > 0 ) then
          print *,'ACTUAL NUMPTS NUM_AA ',numpts, num_aa
       end if

       spo%smt%ia(inode) = iedge

       call imultifab_fill_boundary(spo%index_into_aa)
       call copy_nodal_ind_on_intersect(spo%index_into_aa)

       allocate(spo%smt%ja(num_aa))

       iedge = 1
       do j = jlo,jhi
          do igrid = 1, num_grids_for_j(j)
             ig = new_iarray_ij(j,igrid)

             jbox = get_ibox(spo%index_into_aa,ig)
             at_jhi = (j == box_upb_d(jbox,2))

             ibx = get_pbox(spo%index_into_aa,ig)
             call set_lwb(ibx,2,j-1)
             call set_upb(ibx,2,j+1)
             ind => dataptr(spo%index_into_aa, ig, ibx)

             ibx = get_pbox(mm_grown,ig)
             call set_lwb(ibx,2,j-1)
             call set_upb(ibx,2,j+1)
             mp => dataptr(mm_grown, ig, ibx)

             call create_nodal_ja_2d(spo%smt%ja,mp(:,:,1,1),ind(:,:,1,1),iedge,at_jhi)
          end do
       end do

    case(3)

       jlo = bbox%lo(2)
       jhi = bbox%hi(2)
       klo = bbox%lo(3)
       khi = bbox%hi(3)

       allocate(num_grids_for_jk(jlo:jhi,klo:khi))
       allocate(new_iarray_ijk(jlo:jhi,klo:khi,ngrids))

       do k = klo,khi
          do j = jlo,jhi
             n = 0
             do igrid = 1,ngrids
                b = get_box(ss,igrid)
                if (j >= b%lo(2) .and. j <= b%hi(2) .and. &
                     k >= b%lo(3) .and. k <= b%hi(3)) n = n+1
             end do

             allocate(barray(n))
             allocate(iarray(n))
             allocate(tmp_iarray(n))

             n = 0
             do igrid = 1,ngrids
                b = get_box(ss,igrid)
                if (j >= b%lo(2) .and. j <= b%hi(2) .and. &
                     k >= b%lo(3) .and. k <= b%hi(3)) then
                   n = n+1
                   barray(n) = get_box(ss,igrid)
                   iarray(n) = igrid
                end if
             end do
             num_grids_for_jk(j,k) = n

             call heapsort_indirect_box(barray,tmp_iarray)

             do igrid = 1,num_grids_for_jk(j,k)
                new_iarray_ijk(j,k,igrid) = iarray(tmp_iarray(igrid))
             end do

             deallocate(barray)
             deallocate(iarray)
             deallocate(tmp_iarray)

          end do
       end do

       iedge = 1
       inode = 1

       do k = klo,khi
          do j = jlo,jhi
             do igrid = 1, num_grids_for_jk(j,k)
                ig = new_iarray_ijk(j,k,igrid)
                sbx = get_ibox(ss,ig)
                ibx = sbx

                at_jhi = (j == box_upb_d(ibx,2))
                at_khi = (k == box_upb_d(ibx,3))

                call set_lwb(ibx,2,j)
                call set_upb(ibx,2,j)
                call set_lwb(ibx,3,k)
                call set_upb(ibx,3,k)
                sp => dataptr(ss, ig, ibx)
                ind => dataptr(spo%index_into_aa, ig, ibx)
                mp => dataptr(mm,ig,ibx)
                call create_nodal_aa_3d(sp(:,:,:,:), &
                     mp(:,:,:,1),              &
                     temp_aa,                  &
                     spo%smt%ia,     & 
                     ind(:,:,:,1),             &
                     inode,iedge,at_jhi,at_khi)
             end do
          end do
       end do
       spo%smt%ia(inode) = iedge

       num_aa = iedge - 1
       numpts = inode - 1
       allocate(spo%smt%aa(num_aa))
       do i = 1, num_aa
          spo%smt%aa(i) = temp_aa(i)
       end do

       if ( verbose > 0 ) then
          print *,'ACTUAL NUMPTS NUM_AA ',numpts, num_aa
       end if

       call imultifab_fill_boundary(spo%index_into_aa)
       call copy_nodal_ind_on_intersect(spo%index_into_aa)

       allocate(spo%smt%ja(num_aa))

       iedge = 1
       inode = 1

       do k = klo,khi
          do j = jlo,jhi
             do igrid = 1, num_grids_for_jk(j,k)
                ig = new_iarray_ijk(j,k,igrid)

                sbx = get_pbox(spo%index_into_aa,ig)
                call set_lwb(sbx,2,j-2)
                call set_upb(sbx,2,j+2)
                call set_lwb(sbx,3,k-2)
                call set_upb(sbx,3,k+2)
                ind => dataptr(spo%index_into_aa, ig, sbx)

                mbx = get_box(mm,ig)
                call set_lwb(mbx,2,j)
                call set_upb(mbx,2,j)
                call set_lwb(mbx,3,k)
                call set_upb(mbx,3,k)
                mp => dataptr(mm,ig,mbx)

                at_jhi = (j == box_upb_d(ibx,2))
                at_khi = (k == box_upb_d(ibx,3))

                call create_nodal_ja_3d(spo%smt%ja,mp(:,:,:,1),ind(:,:,:,1),iedge, &
                                        at_jhi,at_khi)
             end do
          end do
       end do

    end select

    ! complete the definition of smt
    spo%smt%numa = num_aa
    spo%smt%nrow = numpts

    call ilut_build(spo%smt, spo%sil)

  contains

    subroutine create_nodal_aa_1d(sp, mp, aa, ia, ind, inode, iedge)
      real (kind = dp_t), intent(in   ) :: sp(:,0:)
      integer           , intent(in   ) :: mp(0:)
      real (kind = dp_t), intent(inout) :: aa(:)
      integer           , intent(inout) :: ia(:)
      integer           , intent(inout) :: ind(0:)
      integer           , intent(inout) :: inode,iedge
      integer nx
      integer i
      integer, parameter :: XBC = 3

      nx = size(sp,dim=1)

      i = 1
      ia(inode) = iedge
      ind(i) = inode

      if (.not. bc_dirichlet(mp(i),1,0)) then
        if (.not. bc_dirichlet(mp(i-1),1,0)) then
          aa(iedge) = sp(i,2) ; iedge = iedge + 1        ! Left
        end if
          aa(iedge) = sp(i,0) ; iedge = iedge + 1        ! Center
        if (.not. bc_dirichlet(mp(i+1),1,0)) then
          aa(iedge) = sp(i,1) ; iedge = iedge + 1        ! Right
        end if
      end if
      inode = inode + 1

      do i = 2, nx-1
         ia(inode) = iedge
         ind(i) = inode
         aa(iedge) = sp(i,2) ; iedge = iedge + 1        ! Left
         aa(iedge) = sp(i,0) ; iedge = iedge + 1        ! Center
         aa(iedge) = sp(i,1) ; iedge = iedge + 1        ! Right
         inode = inode + 1
      end do

!     Only do high side if Neumann; otherwise is Dirichlet (in which case 
!       we don't include it) or it has been taken care of by the grid on the
!       other side.
      i = nx
      ia(inode) = iedge
      ind(i) = inode
      if (bc_neumann(mp(i),1,1)) then
         aa(iedge) = TWO*sp(i,2) ; iedge = iedge + 1    ! Left + Reflection
         aa(iedge) = sp(i,0) ; iedge = iedge + 1        ! Center
      end if
      inode = inode + 1

    end subroutine create_nodal_aa_1d

    subroutine create_nodal_ja_1d(ja, mp, ind, iedge)

      integer, intent(inout) :: ja(:)
      integer, intent(in   ) :: mp(0:)
      integer, intent(inout) :: ind(0:)
      integer, intent(inout) :: iedge

      integer :: i,nx

      nx = size(mp,dim=1)

      i = 1
      if (.not. bc_dirichlet(mp(i),1,0)) then
        if (.not. bc_dirichlet(mp(i-1),1,0)) then
         ja(iedge) = ind(i-1) ; iedge = iedge + 1          ! Left
        end if
         ja(iedge) = ind(i  ) ; iedge = iedge + 1          ! Center
        if (.not. bc_dirichlet(mp(i+1),1,0)) then
         ja(iedge) = ind(i+1) ; iedge = iedge + 1          ! Right
        end if
      end if

      do i = 2, nx-1
         ja(iedge) = ind(i-1) ; iedge = iedge + 1          ! Left
         ja(iedge) = ind(i  ) ; iedge = iedge + 1          ! Center
         ja(iedge) = ind(i+1) ; iedge = iedge + 1          ! Right
      end do

!     Only do high side if Neumann; otherwise is Dirichlet (in which case 
!       we don't include it) or it has been taken care of by the grid on the
!       other side.
      i = nx
      if (bc_neumann(mp(i),1,1)) then
         ja(iedge) = ind(i-1) ; iedge = iedge + 1       ! Left + Reflection
         ja(iedge) = ind(i  ) ; iedge = iedge + 1       ! Center
      end if


    end subroutine create_nodal_ja_1d

    subroutine create_nodal_aa_2d(sp, mp, aa, ia, ind, inode, iedge, at_jhi)

      real (kind = dp_t), intent(in)    :: sp(:,:,0:)
      integer           , intent(in)    :: mp(0:,0:)
      real (kind = dp_t), intent(inout) :: aa(:)
      integer           , intent(inout) :: ia(:)
      integer           , intent(inout) :: ind(:,:)
      integer           , intent(inout) :: inode,iedge
      logical           , intent(in   ) :: at_jhi

      integer :: i,j,nx,ny
      integer :: iedge0

      nx = size(sp,dim=1)
      ny = size(sp,dim=2)

      j = 1
      iedge0 = iedge

      !   **************************************************************************

      do i = 1, nx-1
         if (.not. bc_dirichlet(mp(i,j),1,0)) then
           if (.not. at_jhi .or. bc_neumann(mp(i,j),2,1)) then
            ia(inode) = iedge
            ind(i,j) = inode
            if (.not. bc_dirichlet(mp(i-1,j-1),1,0)) then
              aa(iedge) = sp(i,j,1)                           ! SW
              if (bc_neumann(mp(i,j),2,1)) aa(iedge) = TWO*aa(iedge)
              iedge = iedge + 1
            end if
            if (.not. bc_dirichlet(mp(i  ,j-1),1,0)) then
              aa(iedge) = sp(i,j,2)                           ! S
              if (bc_neumann(mp(i,j),2,1)) aa(iedge) = TWO*aa(iedge)
              iedge = iedge + 1
            end if
            if (.not. bc_dirichlet(mp(i+1,j-1),1,0)) then
              aa(iedge) = sp(i,j,3)                           ! SE
              if (bc_neumann(mp(i,j),1,-1)) aa(iedge) = TWO*aa(iedge)
              if (bc_neumann(mp(i,j),2, 1)) aa(iedge) = TWO*aa(iedge)
              iedge = iedge + 1
            end if
            if (.not. bc_dirichlet(mp(i-1,j  ),1,0)) then
              aa(iedge) = sp(i,j,4)                           ! W
              if (bc_neumann(mp(i,j),1,-1)) aa(iedge) = TWO*aa(iedge)
              iedge = iedge + 1
            end if
              aa(iedge) = sp(i,j,0)                           ! Center
              iedge = iedge + 1
            if (.not. bc_dirichlet(mp(i+1,j  ),1,0)) then
              aa(iedge) = sp(i,j,5)                           ! E
              if (bc_neumann(mp(i,j),1,-1)) aa(iedge) = TWO*aa(iedge)
              iedge = iedge + 1
            end if
            if (.not. bc_dirichlet(mp(i-1,j+1),1,0)) then
              aa(iedge) = sp(i,j,6)                           ! NW
              if (bc_neumann(mp(i,j),2,-1)) aa(iedge) = TWO*aa(iedge)
              iedge = iedge + 1
            end if
            if (.not. bc_dirichlet(mp(i  ,j+1),1,0)) then
              aa(iedge) = sp(i,j,7)                           ! N
              if (bc_neumann(mp(i,j),2,-1)) aa(iedge) = TWO*aa(iedge)
              iedge = iedge + 1
            end if
            if (.not. bc_dirichlet(mp(i+1,j+1),1,0)) then
              aa(iedge) = sp(i,j,8)                           ! NE
              if (bc_neumann(mp(i,j),1,-1)) aa(iedge) = TWO*aa(iedge)
              if (bc_neumann(mp(i,j),2,-1)) aa(iedge) = TWO*aa(iedge)
              iedge = iedge + 1
            end if
            inode = inode + 1
           end if
         end if
      end do

      !   **************************************************************************

      i = nx
!     Only do high side if Neumann; otherwise is Dirichlet (in which case 
!     we don't include it) or it has been taken care of by the grid on the
!     other side.
      if (.not. bc_dirichlet(mp(i,j),1,0) .and. &
                bc_neumann(  mp(i,j),1,1)) then
        if (.not. at_jhi .or. bc_neumann(mp(i,j),2,1)) then
            ia(inode) = iedge
            ind(i,j) = inode
            if (.not. bc_dirichlet(mp(i-1,j-1),1,0)) then
              aa(iedge) = sp(i,j,1)                           ! SW
              if (bc_neumann(mp(i,j),2, 1)) aa(iedge) = TWO*aa(iedge)
              if (bc_neumann(mp(i,j),1, 1)) aa(iedge) = TWO*aa(iedge)
              iedge = iedge + 1
            end if
            if (.not. bc_dirichlet(mp(i  ,j-1),1,0)) then
              aa(iedge) = sp(i,j,2)                           ! S
              if (bc_neumann(mp(i,j),2, 1)) aa(iedge) = TWO*aa(iedge)
              iedge = iedge + 1
            end if
            if (.not. bc_dirichlet(mp(i-1,j  ),1,0)) then
              aa(iedge) = sp(i,j,4)                           ! W
              if (bc_neumann(mp(i,j),1, 1)) aa(iedge) = TWO*aa(iedge)
              iedge = iedge + 1
            end if
              aa(iedge) = sp(i,j,0) ; iedge = iedge + 1       ! Center
            if (.not. bc_dirichlet(mp(i-1,j+1),1,0)) then
              aa(iedge) = sp(i,j,6)                           ! NW
              if (bc_neumann(mp(i,j),2,-1)) aa(iedge) = TWO*aa(iedge)
              if (bc_neumann(mp(i,j),1, 1)) aa(iedge) = TWO*aa(iedge)
              iedge = iedge + 1
            end if
            if (.not. bc_dirichlet(mp(i  ,j+1),1,0)) then
              aa(iedge) = sp(i,j,7)                           ! N
              if (bc_neumann(mp(i,j),2,-1)) aa(iedge) = TWO*aa(iedge)
              iedge = iedge + 1
            end if
            inode = inode + 1
        end if
      end if

    end subroutine create_nodal_aa_2d

    subroutine create_nodal_ja_2d(ja, mp, ind, iedge, at_jhi)

      integer, intent(inout) ::  ja(:)
      integer, intent(in   ) ::  mp(0:,0:)
      integer, intent(inout) :: ind(0:,0:)
      integer, intent(inout) :: iedge
      logical, intent(in   ) :: at_jhi

      integer nx
      integer i,j,ie
      integer iedge0

      nx = size(mp,dim=1)-2
      j = 1

      !**************************************************************************

      iedge0 = iedge
      do i = 1, nx-1
         if (.not. bc_dirichlet(mp(i,j),1,0)) then
           if (.not. at_jhi .or. bc_neumann(mp(i,j),2,1)) then
            if (.not. bc_dirichlet(mp(i-1,j-1),1,0)) then
              ja(iedge) = ind(i-1,j-1) ; iedge = iedge + 1       ! SW
            end if
            if (.not. bc_dirichlet(mp(i  ,j-1),1,0)) then
              ja(iedge) = ind(i  ,j-1) ; iedge = iedge + 1       ! S
            end if
            if (.not. bc_dirichlet(mp(i+1,j-1),1,0)) then
              ja(iedge) = ind(i+1,j-1) ; iedge = iedge + 1       ! SE
            end if
            if (.not. bc_dirichlet(mp(i-1,j  ),1,0)) then
              ja(iedge) = ind(i-1,j  ) ; iedge = iedge + 1       ! W
            end if
              ja(iedge) = ind(i  ,j  ) ; iedge = iedge + 1       ! Center
            if (.not. bc_dirichlet(mp(i+1,j  ),1,0)) then
              ja(iedge) = ind(i+1,j  ) ; iedge = iedge + 1       ! E
            end if
            if (.not. bc_dirichlet(mp(i-1,j+1),1,0)) then
              ja(iedge) = ind(i-1,j+1) ; iedge = iedge + 1       ! NW
            end if
            if (.not. bc_dirichlet(mp(i  ,j+1),1,0)) then
              ja(iedge) = ind(i  ,j+1) ; iedge = iedge + 1       ! N
            end if
            if (.not. bc_dirichlet(mp(i+1,j+1),1,0)) then
              ja(iedge) = ind(i+1,j+1) ; iedge = iedge + 1       ! NE
            end if
           end if
         end if
      end do

      !**************************************************************************

      i = nx
!     Only do high side if Neumann; otherwise is Dirichlet (in which case 
!     we don't include it) or it is taken care of by the grid on the
!     other side.
      if (.not. bc_dirichlet(mp(i,j),1,0) .and. &
                   bc_neumann(  mp(i,j),1,1)) then
        if (.not. at_jhi .or. bc_neumann(mp(i,j),2,1)) then
          if (.not. bc_dirichlet(mp(i-1,j-1),1,0)) then
            ja(iedge) = ind(i-1,j-1) ; iedge = iedge + 1       ! SW
          end if
          if (.not. bc_dirichlet(mp(i  ,j-1),1,0)) then
            ja(iedge) = ind(i  ,j-1) ; iedge = iedge + 1       ! S
          end if
          if (.not. bc_dirichlet(mp(i-1,j  ),1,0)) then
            ja(iedge) = ind(i-1,j  ) ; iedge = iedge + 1       ! W
          end if
            ja(iedge) = ind(i  ,j  ) ; iedge = iedge + 1       ! Center
          if (.not. bc_dirichlet(mp(i-1,j+1),1,0)) then
            ja(iedge) = ind(i-1,j+1) ; iedge = iedge + 1       ! NW
          end if
          if (.not. bc_dirichlet(mp(i  ,j+1),1,0)) then
            ja(iedge) = ind(i  ,j+1) ; iedge = iedge + 1       ! N
          end if
        end if
      end if

    end subroutine create_nodal_ja_2d

    subroutine create_nodal_aa_3d(sp, mp, aa, ia, ind, inode, iedge, at_jhi, at_khi)
      integer           , intent(in) :: mp(:,:,:)
      real (kind = dp_t), intent(in) :: sp(:,:,:,0:)
      real (kind = dp_t), intent(inout) :: aa(:)
      integer, intent(inout) :: ia(:)
      integer, intent(inout) :: ind(:,:,:)
      logical, intent(in   ) :: at_jhi, at_khi
      integer iedge,inode

      integer nx
      integer i,j,k

      nx = size(sp,dim=1)

      k = 1
      j = 1

      !**************************************************************************

      do i = 1, nx-1
       if (.not. bc_dirichlet(mp(i,j,k),1,0)) then
        if ( (.not. at_jhi .or. bc_neumann(mp(i,j,k),2,1)) .and. &
             (.not. at_khi .or. bc_neumann(mp(i,j,k),3,1)) ) then
         ia(inode) = iedge
         ind(i,j,k) = inode
         aa(iedge) = sp(i,j,k,1 ) ; iedge = iedge + 1       ! SW : KLO
         aa(iedge) = sp(i,j,k,2 ) ; iedge = iedge + 1       ! S  : KLO
         aa(iedge) = sp(i,j,k,3 ) ; iedge = iedge + 1       ! SE : KLO
         aa(iedge) = sp(i,j,k,4 ) ; iedge = iedge + 1       ! W  : KLO
         aa(iedge) = sp(i,j,k,25) ; iedge = iedge + 1       ! C  : KLO
         aa(iedge) = sp(i,j,k,5 ) ; iedge = iedge + 1       ! E  : KLO
         aa(iedge) = sp(i,j,k,6 ) ; iedge = iedge + 1       ! NW : KLO
         aa(iedge) = sp(i,j,k,7 ) ; iedge = iedge + 1       ! N  : KLO
         aa(iedge) = sp(i,j,k,8 ) ; iedge = iedge + 1       ! NE : KLO
         aa(iedge) = sp(i,j,k,9 ) ; iedge = iedge + 1       ! SW : KMED
         aa(iedge) = sp(i,j,k,23) ; iedge = iedge + 1       ! S  : KMED
         aa(iedge) = sp(i,j,k,10) ; iedge = iedge + 1       ! SE : KMED
         aa(iedge) = sp(i,j,k,21) ; iedge = iedge + 1       ! W  : KMED
         aa(iedge) = sp(i,j,k,0 ) ; iedge = iedge + 1       ! C  : KMED
         aa(iedge) = sp(i,j,k,22) ; iedge = iedge + 1       ! E  : KMED
         aa(iedge) = sp(i,j,k,11) ; iedge = iedge + 1       ! NW : KMED
         aa(iedge) = sp(i,j,k,24) ; iedge = iedge + 1       ! N  : KMED
         aa(iedge) = sp(i,j,k,12) ; iedge = iedge + 1       ! NE : KMED
         aa(iedge) = sp(i,j,k,13) ; iedge = iedge + 1       ! SW : KHI
         aa(iedge) = sp(i,j,k,14) ; iedge = iedge + 1       ! S  : KHI
         aa(iedge) = sp(i,j,k,15) ; iedge = iedge + 1       ! SE : KHI
         aa(iedge) = sp(i,j,k,16) ; iedge = iedge + 1       ! W  : KHI
         aa(iedge) = sp(i,j,k,26) ; iedge = iedge + 1       ! C  : KHI
         aa(iedge) = sp(i,j,k,17) ; iedge = iedge + 1       ! E  : KHI
         aa(iedge) = sp(i,j,k,18) ; iedge = iedge + 1       ! NW : KHI
         aa(iedge) = sp(i,j,k,19) ; iedge = iedge + 1       ! N  : KHI
         aa(iedge) = sp(i,j,k,20) ; iedge = iedge + 1       ! NE : KHI
         inode = inode + 1
        end if
       end if
      end do

      !**************************************************************************

      if (nx /= 1) then
         i = nx
         ia(inode) = iedge
         ind(i,j,k) = inode
!        Only do high side if Neumann; otherwise is Dirichlet (in which case 
!        we don't include it) or it has been taken care of by the grid on the
!        other side.
         if (.not. bc_dirichlet(mp(i,j,k),1,0) .and. &
                   bc_neumann(  mp(i,j,k),1,1)) then
          if ( (.not. at_jhi .or. bc_neumann(mp(i,j,k),2,1)) .and. &
               (.not. at_khi .or. bc_neumann(mp(i,j,k),3,1)) ) then
            aa(iedge) = sp(i,j,k,1 ) ; iedge = iedge + 1       ! SW : KLO
            aa(iedge) = sp(i,j,k,2 ) ; iedge = iedge + 1       ! S  : KLO
            aa(iedge) = sp(i,j,k,3 ) ; iedge = iedge + 1       ! SE : KLO
            aa(iedge) = sp(i,j,k,4 ) ; iedge = iedge + 1       ! W  : KLO
            aa(iedge) = sp(i,j,k,25) ; iedge = iedge + 1       ! C  : KLO
            aa(iedge) = sp(i,j,k,5 ) ; iedge = iedge + 1       ! E  : KLO
            aa(iedge) = sp(i,j,k,6 ) ; iedge = iedge + 1       ! NW : KLO
            aa(iedge) = sp(i,j,k,7 ) ; iedge = iedge + 1       ! N  : KLO
            aa(iedge) = sp(i,j,k,8 ) ; iedge = iedge + 1       ! NE : KLO
            aa(iedge) = sp(i,j,k,9 ) ; iedge = iedge + 1       ! SW : KMED
            aa(iedge) = sp(i,j,k,23) ; iedge = iedge + 1       ! S  : KMED
            aa(iedge) = sp(i,j,k,10) ; iedge = iedge + 1       ! SE : KMED
            aa(iedge) = sp(i,j,k,21) ; iedge = iedge + 1       ! W  : KMED
            aa(iedge) = sp(i,j,k,0 ) ; iedge = iedge + 1       ! C  : KMED
            aa(iedge) = sp(i,j,k,22) ; iedge = iedge + 1       ! E  : KMED
            aa(iedge) = sp(i,j,k,11) ; iedge = iedge + 1       ! NW : KMED
            aa(iedge) = sp(i,j,k,24) ; iedge = iedge + 1       ! N  : KMED
            aa(iedge) = sp(i,j,k,12) ; iedge = iedge + 1       ! NE : KMED
            aa(iedge) = sp(i,j,k,13) ; iedge = iedge + 1       ! SW : KHI
            aa(iedge) = sp(i,j,k,14) ; iedge = iedge + 1       ! S  : KHI
            aa(iedge) = sp(i,j,k,15) ; iedge = iedge + 1       ! SE : KHI
            aa(iedge) = sp(i,j,k,16) ; iedge = iedge + 1       ! W  : KHI
            aa(iedge) = sp(i,j,k,26) ; iedge = iedge + 1       ! C  : KHI
            aa(iedge) = sp(i,j,k,17) ; iedge = iedge + 1       ! E  : KHI
            aa(iedge) = sp(i,j,k,18) ; iedge = iedge + 1       ! NW : KHI
            aa(iedge) = sp(i,j,k,19) ; iedge = iedge + 1       ! N  : KHI
            aa(iedge) = sp(i,j,k,20) ; iedge = iedge + 1       ! NE : KHI
            inode = inode + 1
          end if
         end if
      end if

    end subroutine create_nodal_aa_3d

    subroutine create_nodal_ja_3d(ja, mp, ind, iedge, at_jhi, at_khi)
      integer, intent(inout) :: ja(:)
      integer, intent(in   ) :: mp(:,:,:)
      integer, intent(inout) :: ind(-1:,-1:,-1:)
      logical, intent(in   ) :: at_jhi, at_khi
      integer iedge

      integer nx
      integer i,j,k

      nx = size(mp,dim=1)

      k = 1
      j = 1

      !**************************************************************************

      i = 1
      if (.not. bc_dirichlet(mp(i,j,k),1,0)) then
       if ( (.not. at_jhi .or. bc_neumann(mp(i,j,k),2,1)) .and. &
            (.not. at_khi .or. bc_neumann(mp(i,j,k),3,1)) ) then
         ja(iedge) = ind(i-1,j-1,k-1) ; iedge = iedge + 1       ! SW : KLO
         ja(iedge) = ind(i  ,j-1,k-1) ; iedge = iedge + 1       ! S  : KLO
         ja(iedge) = ind(i+1,j-1,k-1) ; iedge = iedge + 1       ! SE : KLO
         ja(iedge) = ind(i-1,j  ,k-1) ; iedge = iedge + 1       ! W  : KLO
         ja(iedge) = ind(i  ,j  ,k-1) ; iedge = iedge + 1       ! C  : KLO
         ja(iedge) = ind(i+1,j  ,k-1) ; iedge = iedge + 1       ! E  : KLO
         ja(iedge) = ind(i-1,j+1,k-1) ; iedge = iedge + 1       ! NW : KLO
         ja(iedge) = ind(i  ,j+1,k-1) ; iedge = iedge + 1       ! N  : KLO
         ja(iedge) = ind(i+1,j+1,k-1) ; iedge = iedge + 1       ! NE : KLO
         ja(iedge) = ind(i-1,j-1,k  ) ; iedge = iedge + 1       ! SW : KMED
         ja(iedge) = ind(i  ,j-1,k  ) ; iedge = iedge + 1       ! S  : KMED
         ja(iedge) = ind(i+1,j-1,k  ) ; iedge = iedge + 1       ! SE : KMED
         ja(iedge) = ind(i-1,j  ,k  ) ; iedge = iedge + 1       ! W  : KMED
         ja(iedge) = ind(i  ,j  ,k  ) ; iedge = iedge + 1       ! C  : KMED
         ja(iedge) = ind(i+1,j  ,k  ) ; iedge = iedge + 1       ! E  : KMED
         ja(iedge) = ind(i-1,j+1,k  ) ; iedge = iedge + 1       ! NW : KMED
         ja(iedge) = ind(i  ,j+1,k  ) ; iedge = iedge + 1       ! N  : KMED
         ja(iedge) = ind(i+1,j+1,k  ) ; iedge = iedge + 1       ! NE : KMED
         ja(iedge) = ind(i-1,j-1,k+1) ; iedge = iedge + 1       ! SW : KHI
         ja(iedge) = ind(i  ,j-1,k+1) ; iedge = iedge + 1       ! S  : KHI
         ja(iedge) = ind(i+1,j-1,k+1) ; iedge = iedge + 1       ! SE : KHI
         ja(iedge) = ind(i-1,j  ,k+1) ; iedge = iedge + 1       ! W  : KHI
         ja(iedge) = ind(i  ,j  ,k+1) ; iedge = iedge + 1       ! C  : KHI
         ja(iedge) = ind(i+1,j  ,k+1) ; iedge = iedge + 1       ! E  : KHI
         ja(iedge) = ind(i-1,j+1,k+1) ; iedge = iedge + 1       ! NW : KHI
         ja(iedge) = ind(i  ,j+1,k+1) ; iedge = iedge + 1       ! N  : KHI
         ja(iedge) = ind(i+1,j+1,k+1) ; iedge = iedge + 1       ! NE : KHI
       end if
      end if

      !**************************************************************************

      do i = 2, nx-1
       if ( (.not. at_jhi .or. bc_neumann(mp(i,j,k),2,1)) .and. &
            (.not. at_khi .or. bc_neumann(mp(i,j,k),3,1)) ) then
         ja(iedge) = ind(i-1,j-1,k-1) ; iedge = iedge + 1       ! SW : KLO
         ja(iedge) = ind(i  ,j-1,k-1) ; iedge = iedge + 1       ! S  : KLO
         ja(iedge) = ind(i+1,j-1,k-1) ; iedge = iedge + 1       ! SE : KLO
         ja(iedge) = ind(i-1,j  ,k-1) ; iedge = iedge + 1       ! W  : KLO
         ja(iedge) = ind(i  ,j  ,k-1) ; iedge = iedge + 1       ! C  : KLO
         ja(iedge) = ind(i+1,j  ,k-1) ; iedge = iedge + 1       ! E  : KLO
         ja(iedge) = ind(i-1,j+1,k-1) ; iedge = iedge + 1       ! NW : KLO
         ja(iedge) = ind(i  ,j+1,k-1) ; iedge = iedge + 1       ! N  : KLO
         ja(iedge) = ind(i+1,j+1,k-1) ; iedge = iedge + 1       ! NE : KLO
         ja(iedge) = ind(i-1,j-1,k  ) ; iedge = iedge + 1       ! SW : KMED
         ja(iedge) = ind(i  ,j-1,k  ) ; iedge = iedge + 1       ! S  : KMED
         ja(iedge) = ind(i+1,j-1,k  ) ; iedge = iedge + 1       ! SE : KMED
         ja(iedge) = ind(i-1,j  ,k  ) ; iedge = iedge + 1       ! W  : KMED
         ja(iedge) = ind(i  ,j  ,k  ) ; iedge = iedge + 1       ! C  : KMED
         ja(iedge) = ind(i+1,j  ,k  ) ; iedge = iedge + 1       ! E  : KMED
         ja(iedge) = ind(i-1,j+1,k  ) ; iedge = iedge + 1       ! NW : KMED
         ja(iedge) = ind(i  ,j+1,k  ) ; iedge = iedge + 1       ! N  : KMED
         ja(iedge) = ind(i+1,j+1,k  ) ; iedge = iedge + 1       ! NE : KMED
         ja(iedge) = ind(i-1,j-1,k+1) ; iedge = iedge + 1       ! SW : KHI
         ja(iedge) = ind(i  ,j-1,k+1) ; iedge = iedge + 1       ! S  : KHI
         ja(iedge) = ind(i+1,j-1,k+1) ; iedge = iedge + 1       ! SE : KHI
         ja(iedge) = ind(i-1,j  ,k+1) ; iedge = iedge + 1       ! W  : KHI
         ja(iedge) = ind(i  ,j  ,k+1) ; iedge = iedge + 1       ! C  : KHI
         ja(iedge) = ind(i+1,j  ,k+1) ; iedge = iedge + 1       ! E  : KHI
         ja(iedge) = ind(i-1,j+1,k+1) ; iedge = iedge + 1       ! NW : KHI
         ja(iedge) = ind(i  ,j+1,k+1) ; iedge = iedge + 1       ! N  : KHI
         ja(iedge) = ind(i+1,j+1,k+1) ; iedge = iedge + 1       ! NE : KHI
       end if
      end do

      !**************************************************************************

      if (nx /= 1) then
         i = nx
!        Only do high side if Neumann; otherwise is Dirichlet (in which case 
!        we don't include it) or it has been taken care of by the grid on the
!        other side.
         if (.not. bc_dirichlet(mp(i,j,k),1,0) .and. &
                   bc_neumann(  mp(i,j,k),1,1)) then
          if ( (.not. at_jhi .or. bc_neumann(mp(i,j,k),2,1)) .and. &
               (.not. at_khi .or. bc_neumann(mp(i,j,k),3,1)) ) then
            ja(iedge) = ind(i-1,j-1,k-1) ; iedge = iedge + 1       ! SW : KLO
            ja(iedge) = ind(i  ,j-1,k-1) ; iedge = iedge + 1       ! S  : KLO
            ja(iedge) = ind(i+1,j-1,k-1) ; iedge = iedge + 1       ! SE : KLO
            ja(iedge) = ind(i-1,j  ,k-1) ; iedge = iedge + 1       ! W  : KLO
            ja(iedge) = ind(i  ,j  ,k-1) ; iedge = iedge + 1       ! C  : KLO
            ja(iedge) = ind(i+1,j  ,k-1) ; iedge = iedge + 1       ! E  : KLO
            ja(iedge) = ind(i-1,j+1,k-1) ; iedge = iedge + 1       ! NW : KLO
            ja(iedge) = ind(i  ,j+1,k-1) ; iedge = iedge + 1       ! N  : KLO
            ja(iedge) = ind(i+1,j+1,k-1) ; iedge = iedge + 1       ! NE : KLO
            ja(iedge) = ind(i-1,j-1,k  ) ; iedge = iedge + 1       ! SW : KMED
            ja(iedge) = ind(i  ,j-1,k  ) ; iedge = iedge + 1       ! S  : KMED
            ja(iedge) = ind(i+1,j-1,k  ) ; iedge = iedge + 1       ! SE : KMED
            ja(iedge) = ind(i-1,j  ,k  ) ; iedge = iedge + 1       ! W  : KMED
            ja(iedge) = ind(i  ,j  ,k  ) ; iedge = iedge + 1       ! C  : KMED
            ja(iedge) = ind(i+1,j  ,k  ) ; iedge = iedge + 1       ! E  : KMED
            ja(iedge) = ind(i-1,j+1,k  ) ; iedge = iedge + 1       ! NW : KMED
            ja(iedge) = ind(i  ,j+1,k  ) ; iedge = iedge + 1       ! N  : KMED
            ja(iedge) = ind(i+1,j+1,k  ) ; iedge = iedge + 1       ! NE : KMED
            ja(iedge) = ind(i-1,j-1,k+1) ; iedge = iedge + 1       ! SW : KHI
            ja(iedge) = ind(i  ,j-1,k+1) ; iedge = iedge + 1       ! S  : KHI
            ja(iedge) = ind(i+1,j-1,k+1) ; iedge = iedge + 1       ! SE : KHI
            ja(iedge) = ind(i-1,j  ,k+1) ; iedge = iedge + 1       ! W  : KHI
            ja(iedge) = ind(i  ,j  ,k+1) ; iedge = iedge + 1       ! C  : KHI
            ja(iedge) = ind(i+1,j  ,k+1) ; iedge = iedge + 1       ! E  : KHI
            ja(iedge) = ind(i-1,j+1,k+1) ; iedge = iedge + 1       ! NW : KHI
            ja(iedge) = ind(i  ,j+1,k+1) ; iedge = iedge + 1       ! N  : KHI
            ja(iedge) = ind(i+1,j+1,k+1) ; iedge = iedge + 1       ! NE : KHI
          end if
         end if
      end if

    end subroutine create_nodal_ja_3d

  subroutine ilut_build(amat, sil)
    type(sparse_matrix), intent(in) :: amat
    type(sparse_matrix), intent(out) :: sil
    integer :: ierr
    integer, allocatable :: iw(:)
    real(kind=dp_t), allocatable ::  wk(:)
    real(kind=dp_t) :: tol
    integer :: nwk
    integer :: lfil, nrow
    external ilut

    nrow = amat%nrow

    lfil = 4
    tol = 1.e-5_dp_t

    ! FIXME: this is a guess as to the maximum size of these arrays.
    allocate(wk(nrow+1))
    allocate(iw(2*nrow))

    call sparse_matrix_build(sil, (2*lfil+1)*nrow, nrow)

    nwk = size(sil%aa)
    call ilut (nrow, amat%aa, amat%ja, amat%ia, lfil, tol, sil%aa, sil%ja, sil%ia, nwk, wk, iw, ierr)
    if ( ierr /= 0 ) then
       write(*,*) 'ILUT: IERR = ', ierr
       call bl_error("ILUT: failed")
    end if

  end subroutine ilut_build

  subroutine copy_nodal_ind_on_intersect(ind)

    type(imultifab) , intent(inout) :: ind

    integer        , pointer :: ind_i(:,:,:,:)
    integer        , pointer :: ind_j(:,:,:,:)

    integer :: lo_i(ind%dim),hi_i(ind%dim)
    integer :: lo_j(ind%dim),hi_j(ind%dim)
    integer :: lo_a(ind%dim),hi_a(ind%dim)
    integer :: i,j,igrid,jgrid,ngrids,idm
    logical :: do_copy

    type(box) :: ibx, jbx, abx

    ngrids = ind%nboxes

!   Copy on intersect from the grid which owns ind() to the grid which doesn't.
    if (ngrids > 1) then
      do igrid = 1,ngrids
        ibx = get_ibox(ind,igrid)

        do jgrid = 1,ngrids
         jbx = get_ibox(ind,jgrid)

         abx = intersection(ibx,jbx)
         if (.not. empty(abx)) then
           ind_i => dataptr(ind, igrid)
           ind_j => dataptr(ind, jgrid)
           lo_i = lwb(ibx)
           hi_i = upb(ibx)
           lo_j = lwb(jbx)
           hi_j = upb(jbx)
           lo_a = lwb(abx)
           hi_a = upb(abx)
           do_copy = .false.
           do idm = 1,ind%dim
             if (lo_i(idm) == hi_j(idm)) do_copy = .true.
           end do
           if (do_copy) then
             select case(ind%dim)
             case(1)
               if (ind_i(lo_a(1),1,1,1) .gt. -1) ind_j(i,1,1,1) = ind_i(i,1,1,1)
             case(2)
               do j = lo_a(2),hi_a(2)
               do i = lo_a(1),hi_a(1)
                 if (ind_i(i,j,1,1) .gt. -1) then
                  ind_j(i,j,1,1) = ind_i(i,j,1,1)
                 end if
               end do
               end do
             case(3)
               do k = lo_a(3),hi_a(3)
               do j = lo_a(2),hi_a(2)
               do i = lo_a(1),hi_a(1)
                 if (ind_i(i,j,k,1) .gt. -1) &
                  ind_j(i,j,k,1) = ind_i(i,j,k,1)
               end do
               end do
               end do
             end select
           end if
         end if
        end do
       end do
    end if

  end subroutine copy_nodal_ind_on_intersect

  end subroutine sparse_nodal_build

  subroutine sparse_solve(spo, uu, rh, eps, max_iter, verbose, stat)

    type(sparse), intent(in) :: spo
    type(multifab), intent(inout) :: uu
    type(multifab), intent(in)    :: rh

    real(dp_t), intent(in) :: eps
    integer, intent(in) :: max_iter, verbose
    integer, intent(out) :: stat

    real(kind=dp_t), pointer :: rp(:,:,:,:)

    integer        , pointer :: ind(:,:,:,:)

    real (kind = dp_t), allocatable ::  rhs(:)
    real (kind = dp_t), allocatable ::  soln(:)

    type(box) ibx

    integer lo(rh%dim)
    integer hi(rh%dim)
    integer i,j,k
    integer igrid,ngrids
    integer numpts, num_aa

    ngrids = rh%nboxes

    numpts = spo%smt%nrow
    num_aa = spo%smt%numa

    allocate(soln(numpts))
    allocate(rhs(numpts))
    soln = zero

    !   Put the RHS into the 1-d RHS vector
    do igrid = 1, ngrids
       rp => dataptr(rh,igrid)
       ind => dataptr(spo%index_into_aa, igrid)

       ibx = get_box(rh,igrid)
       lo = lwb(ibx)
       hi = upb(ibx)

       select case(rh%dim)
       case(1)
          do i = lo(1),hi(1)
             rhs(ind(i,1,1,1)) = rp(i,1,1,1)
          end do
       case (2)
          !$OMP PARALLEL DO PRIVATE(i,j)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                rhs(ind(i,j,1,1)) = rp(i,j,1,1)
             end do
          end do
          !$OMP END PARALLEL DO
       case (3)
          !$OMP PARALLEL DO PRIVATE(i,j,k)
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                   rhs(ind(i,j,k,1)) = rp(i,j,k,1)
                end do
             end do
          end do
          !$OMP END PARALLEL DO
       end select
    end do

    call sparse_solver(spo%smt, spo%sil, eps, max_iter, rhs, soln, spo%Anorm, verbose, stat)

    !   Put the 1-d solution back into the N-d array.
    do igrid = 1, ngrids
       rp => dataptr(uu,igrid)
       ind => dataptr(spo%index_into_aa, igrid)

       ibx = get_box(uu,igrid)
       lo = lwb(ibx)
       hi = upb(ibx)

       select case(rh%dim)
       case(1)
          do i = lo(1),hi(1)
             rp(i,1,1,1) = soln(ind(i,1,1,1))
          end do
       case(2)
          !$OMP PARALLEL DO PRIVATE(i,j)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                rp(i,j,1,1) = soln(ind(i,j,1,1))
             end do
          end do
          !$OMP END PARALLEL DO
       case(3)
          !$OMP PARALLEL DO PRIVATE(i,j,k)
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                   rp(i,j,k,1) = soln(ind(i,j,k,1))
                end do
             end do
          end do
          !$OMP END PARALLEL DO
       end select
    end do

  end subroutine sparse_solve

  subroutine sparse_nodal_solve(spo, uu, rh, eps, max_iter, verbose, stat)

    type(sparse), intent(in) :: spo
    type(multifab), intent(inout) :: uu
    type(multifab), intent(in)    :: rh

    real(dp_t), intent(in) :: eps
    integer, intent(in) :: max_iter, verbose
    integer, intent(out) :: stat

    real(kind=dp_t), pointer :: rp(:,:,:,:)
    real(kind=dp_t), pointer :: rp_j(:,:,:,:)

    integer        , pointer :: ind(:,:,:,:)

    real (kind = dp_t), allocatable ::  rhs(:)
    real (kind = dp_t), allocatable ::  soln(:)

    type(box) :: ibx,jbx,abx
    logical   :: do_copy

    integer lo(rh%dim),lo_j(rh%dim),lo_a(rh%dim)
    integer hi(rh%dim),hi_j(rh%dim),hi_a(rh%dim)
    integer i,j,k,idm
    integer igrid,jgrid,ngrids
    integer numpts, num_aa

    ngrids = rh%nboxes

    numpts = spo%smt%nrow
    num_aa = spo%smt%numa

    allocate(soln(numpts))
    allocate(rhs(numpts))
    soln = zero

    !   Put the RHS into the 1-d RHS vector
    do igrid = 1, ngrids
       rp => dataptr(rh,igrid)
       ind => dataptr(spo%index_into_aa, igrid)

       ibx = get_ibox(rh,igrid)
       lo = lwb(ibx)
       hi = upb(ibx)

       select case(rh%dim)
       case(1)
          do i = lo(1),hi(1)
             if (ind(i,1,1,1) .gt. -1) &
                rhs(ind(i,1,1,1)) = rp(i,1,1,1)
          end do
       case (2)
          !$OMP PARALLEL DO PRIVATE(i,j)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
               if (ind(i,j,1,1) .gt. -1) then
                 rhs(ind(i,j,1,1)) = rp(i,j,1,1)
               end if
             end do
          end do
          !$OMP END PARALLEL DO
       case (3)
          !$OMP PARALLEL DO PRIVATE(i,j,k)
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                   if (ind(i,j,k,1) .gt. -1) &
                      rhs(ind(i,j,k,1)) = rp(i,j,k,1)
                end do
             end do
          end do
          !$OMP END PARALLEL DO
       end select
    end do

    call sparse_solver(spo%smt, spo%sil, eps, max_iter, rhs, soln, spo%Anorm, verbose, stat)

    !   Put the 1-d solution back into the N-d array.
    do igrid = 1, ngrids
       rp => dataptr(uu,igrid)
       ind => dataptr(spo%index_into_aa, igrid)

       ibx = get_ibox(uu,igrid)
       lo = lwb(ibx)
       hi = upb(ibx)

       select case(rh%dim)
       case(1)
          do i = lo(1),hi(1)
             if (ind(i,1,1,1) .gt. -1) &
                rp(i,1,1,1) = soln(ind(i,1,1,1))
          end do
       case(2)
          !$OMP PARALLEL DO PRIVATE(i,j)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
               if (ind(i,j,1,1) > -1) then
                 rp(i,j,1,1) = soln(ind(i,j,1,1))
               end if
             end do
          end do
          !$OMP END PARALLEL DO
       case(3)
          !$OMP PARALLEL DO PRIVATE(i,j,k)
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                   if (ind(i,j,k,1) .gt. -1) &
                   rp(i,j,k,1) = soln(ind(i,j,k,1))
                end do
             end do
          end do
          !$OMP END PARALLEL DO
       end select

!      Copy on intersect from the grid which owns the solution to the grid which doesn't.
       if (ngrids > 1) then
        do jgrid = 1,ngrids
         jbx = get_ibox(uu,jgrid)
         abx = intersection(ibx,jbx)
         if (.not. empty(abx)) then
           rp_j => dataptr(uu,jgrid)
           lo_j = lwb(jbx)
           hi_j = upb(jbx)
           lo_a = lwb(abx)
           hi_a = upb(abx)
           do_copy = .false.
           do idm = 1,rh%dim
             if (lo(idm) == hi_j(idm)) do_copy = .true.
           end do
           if (do_copy) then
             select case(rh%dim)
             case(1)
               if (ind(lo_a(1),1,1,1) .gt. -1) rp_j(i,1,1,1) = rp(i,1,1,1)
             case(2)
               do j = lo_a(2),hi_a(2)
               do i = lo_a(1),hi_a(1)
                 if (ind(i,j,1,1) .gt. -1) &
                  rp_j(i,j,1,1) = rp(i,j,1,1)
               end do
               end do
             case(3)
               do k = lo_a(3),hi_a(3)
               do j = lo_a(2),hi_a(2)
               do i = lo_a(1),hi_a(1)
                 if (ind(i,j,k,1) .gt. -1) &
                  rp_j(i,j,k,1) = rp(i,j,k,1)
               end do
               end do
               end do
             end select
           end if
         end if
        end do
       end if
    end do

  end subroutine sparse_nodal_solve

  subroutine sparse_solver( &
       smt, sil, &
       eps, maxits, &
       rhs, sol, &
       Anorm, verbose, stat)
    type(sparse_matrix), intent(in) :: smt
    type(sparse_matrix), intent(in) :: sil
    integer, intent(in) :: maxits
    integer, intent(in) :: verbose
    real(kind=dp_t), intent(inout)  :: sol(:), rhs(:)
    real(kind=dp_t), intent(in) :: eps, Anorm
    integer, intent(out) :: stat

    integer :: numa
    integer :: nrow
    integer :: ipar(16)
    integer :: i
    real(kind=dp_t), allocatable ::  wk(:)
    real(kind=dp_t) fpar(16)
    real(kind=dp_t) bnorm
    integer :: nwrk
    integer :: kryord, slvr

    integer, parameter :: SOLVER_CG       = 1
    integer, parameter :: SOLVER_CGNR     = 2
    integer, parameter :: SOLVER_BCG      = 3
    integer, parameter :: SOLVER_DBCG     = 4
    integer, parameter :: SOLVER_BCGSTAB  = 5
    integer, parameter :: SOLVER_TFQMR    = 6
    integer, parameter :: SOLVER_FOM      = 7
    integer, parameter :: SOLVER_GMRES    = 8
    integer, parameter :: SOLVER_FGMRES   = 9
    integer, parameter :: SOLVER_DQGMRES  = 10

    external cg, cgnr, bcg, dbcg, bcgstab, tfqmr, fom, gmres, fgmres, dqgmres

    stat = 0
    slvr = SOLVER_BCGSTAB

    nrow = smt%nrow
    numa = smt%numa
    if (nrow == 1 .and. numa == 1) then
       sol(1) = rhs(1) / smt%aa(1)
       return
    end if

    bnorm = maxval(abs(rhs))
    kryord = 15

    select case (slvr)
    case (SOLVER_CG,SOLVER_CGNR)
       nwrk = 5*nrow
    case (SOLVER_BCG)
       nwrk = 7*nrow
    case (SOLVER_DBCG)
       nwrk = 11*nrow
    case (SOLVER_BCGSTAB)
       nwrk = 8*nrow
    case (SOLVER_TFQMR)
       nwrk = 11*nrow
    case (SOLVER_FOM,SOLVER_GMRES)
       nwrk = (nrow+3)*(kryord+2) + (kryord+1)*kryord/2
    case (SOLVER_FGMRES)
       nwrk = 2*nrow*(kryord+1) + (kryord+1)*kryord/2 + 3*kryord + 2
    case (SOLVER_DQGMRES)
       nwrk = nrow + (kryord + 1)*(2*nrow+4)
    case default
       call bl_error("SPARSE_SOLVER: unknown SOLVER")
    end select

    !       CG      == 5 * n
    !       CGNR    == 5 * n
    !       BCG     == 7 * n
    !       DBCG    == 11 * n
    !       BCGSTAB == 8 * n
    !       TFQMR   == 11 * n
    !       FOM     == (n+3)*(m+2) + (m+1)*m/2 (m = ipar(5), default m=15)
    !       GMRES   == (n+3)*(m+2) + (m+1)*m/2 (m = ipar(5), default m=15)
    !       FGMRES  == 2*n*(m+1) + (m+1)*m/2 + 3*m + 2 (m = ipar(5), default m=15)
    !       DQGMRES == n + lb * (2*n+4) (lb=ipar(5)+1, default lb = 16)

    !
    !     set the parameters for the iterative solvers
    !

    ipar(2) = 2                 ! right preconditioning only
    ipar(3) = 1 ! 999           ! stopping criteria
    ipar(4) = nwrk              ! length of wk array
    ipar(5) = kryord            ! size of KRYLOV subspace 
    ipar(6) = maxits            ! max number of Matrix Vector Multiplies
    fpar(1) = eps               ! The relative tolerance
    fpar(2) = 1.0e-10_dp_t      ! The absoulte tolerance
    fpar(2) = 0.0_dp_t
    fpar(11) = 0.0_dp_t         ! Initialize Flop Count

    allocate(wk(nwrk))

    do i = 1, nrow
       sol(i) = 0.0_dp_t
    end do

    select case (slvr)
    case (SOLVER_CG     )
       call runrc(nrow, rhs, sol, ipar, fpar, wk, smt, sil, eps, Anorm, bnorm, verbose, cg)
    case (SOLVER_CGNR   )
       call runrc(nrow, rhs, sol, ipar, fpar, wk, smt, sil, eps, Anorm, bnorm, verbose, cgnr)
    case (SOLVER_BCG    )
       call runrc(nrow, rhs, sol, ipar, fpar, wk, smt, sil, eps, Anorm, bnorm, verbose, bcg)
    case (SOLVER_DBCG   )
       call runrc(nrow, rhs, sol, ipar, fpar, wk, smt, sil, eps, Anorm, bnorm, verbose, dbcg)
    case (SOLVER_BCGSTAB)
       call runrc(nrow, rhs, sol, ipar, fpar, wk, smt, sil, eps, Anorm, bnorm, verbose, bcgstab)
    case (SOLVER_TFQMR  )
       call runrc(nrow, rhs, sol, ipar, fpar, wk, smt, sil, eps, Anorm, bnorm, verbose, tfqmr)
    case (SOLVER_FOM    )
       call runrc(nrow, rhs, sol, ipar, fpar, wk, smt, sil, eps, Anorm, bnorm, verbose, fom)
    case (SOLVER_GMRES  )
       call runrc(nrow, rhs, sol, ipar, fpar, wk, smt, sil, eps, Anorm, bnorm, verbose, gmres)
    case (SOLVER_FGMRES )
       call runrc(nrow, rhs, sol, ipar, fpar, wk, smt, sil, eps, Anorm, bnorm, verbose, fgmres)
    case (SOLVER_DQGMRES)
       call runrc(nrow, rhs, sol, ipar, fpar, wk, smt, sil, eps, Anorm, bnorm, verbose, dqgmres)
    case default
       call bl_error("SPARSE_SOLVER: unknown solver")
    end select
    if ( ipar(1) == 999 ) then
       stat = -1
       return
    end if
    if ( ipar(1) /= 0 ) then
       write(*,*) 'RUNRC: ipar(1) ', ipar(1)
       call bl_error("RUNRC: failed")
    end if

  contains

    subroutine runrc(n, rhs, sol, ipar, fpar, wk, &
         smt, sil, &
         eps, Anorm, bnorm, verbose, &
         solver)
      type(sparse_matrix) :: smt, sil
      integer, intent(in) :: n
      real(kind = dp_t) :: rhs(:), sol(:), fpar(16), wk(*)
      real(kind = dp_t), intent(in) :: Anorm, bnorm, eps
      integer :: ipar(16)
      external solver
      integer, intent(in) ::  verbose

      !     the actual tester. It starts the iterative linear system solvers
      !     with a initial guess suppied by the user.
      !
      !     The structure {au, jau, ju} is assumed to have the output from
      !     the ILU* routines in ilut.f.
      !
      integer i, its
      real(kind = dp_t) res, dnrm2
      real(kind = dp_t), allocatable :: rr(:)
      external dnrm2
      !
      !     ipar(2) can be 0, 1, 2, please don't use 3
      !
      allocate(rr(n))
      if (ipar(2).gt.2) then
         print *, "I can not do both left and right preconditioning."
         ipar(1) = 1
         return
      end if
      !
      !     normal execution
      !
      its = 0
      res = 0.0_dp_t
      !
      do i = 1, n
         sol(i) = 0.0_dp_t
      end do
      !
      ipar(1) = 0

      do
         call solver(n, rhs, sol, ipar, fpar, wk)
         !
         !     output the residuals
         !

         if (ipar(7).ne.its) then
            if ( verbose .gt. 1 ) then
               write (*, *) its, real(res)
            end if
            its = ipar(7)
         end if

         res = fpar(5)
         !
         !      print *, ipar(1), ipar(3), ipar(8), ipar(9)
         if ( ipar(1) == 0 ) exit

         !      call amux(n, sol, rr, a, ja, ia)
         !      print *, 'bnorm = ', bnorm
         !      print *, 'Anorm = ', Anorm
         !     print *, 'newtest ', maxval(abs(rr-rhs)), eps*(Anorm*maxval(abs(sol)) + bnorm )
         select case (ipar(1))
         case (1)
            call amux(n, wk(ipar(8)), wk(ipar(9)), smt%aa, smt%ja, smt%ia)
            cycle
         case(2)
            call atmux(n, wk(ipar(8)), wk(ipar(9)), smt%aa, smt%ja, smt%ia)
            cycle
         case (3,5)
            call lusol(n,wk(ipar(8)),wk(ipar(9)), sil%aa, sil%ja, sil%ia)
            cycle
         case (4,6)
            call lutsol(n,wk(ipar(8)),wk(ipar(9)),sil%aa, sil%ja, sil%ia)
            cycle
         case (10)
            ipar(11) = 0
            ! convergence test goes here
            ! res < eps*(Anorm*norm_inf(uu) + bnorm)
            if ( ipar(8) < 0 ) then
               call bl_error("can't use this convergence test")
            end if
            rr = wk(ipar(8):ipar(8)+n-1)/sil%aa(1:n)
            print  *, 'amux 2 = ', maxval(abs(rr))
            call lusol(n, wk(ipar(8)), rr, sil%aa, sil%ja, sil%ia)
            print  *, 'amux 1 = ', maxval(abs(rr))
            call amux(n, rr, wk(ipar(9)), smt%aa, smt%ja, smt%ia)
            if ( maxval(abs(wk(ipar(9):ipar(9)+n-1) - rhs)) &
                 < eps*(Anorm*maxval(abs(wk(ipar(8):ipar(8)+n-1) + bnorm))) ) then
               ipar(11) = 1
            end if
            cycle
         case (-1)
            print *, 'Iterative solver has iterated too many times.'
         case (-2)
            print *, 'Iterative solver was not given enough work space.'
            print *, 'The work space should at least have ', ipar(4), ' elements.'
         case (-3)
            !! FIXME! this is a bad hack!
            !! unforuntately the test in case(10) above can't be used because the 
            !! current solution is not in wk(ipar(8))
            sol = 0.0_dp_t
            ipar(1) = 999
            exit
            print *, 'Iterative solver is facing a break-down.'
         case default
            print *, 'Iterative solver terminated. code =', ipar(1)
         end select
         exit
      end do

      if ( verbose .gt. 0 ) then
         write (*, fmt = '("SPARSE: iterations: ", i5)' ) ipar(7)
      end if
      if ( verbose > 1 ) then
         write (*, *) '# return code =', ipar(1), '        convergence rate =', fpar(7)

         ! check the error

         call amux(n, sol, wk, smt%aa,smt%ja,smt%ia)
         do i = 1, n
            wk(i) = wk(i) - rhs(i)
         end do
         write (*, *) '# the actual residual norm is', dnrm2(n,wk,1)

      end if

    end subroutine runrc

  end subroutine sparse_solver

end module sparse_solve_module

function distdot(n,x,ix,y,iy)
  use bl_types
  integer n, ix, iy
  real(kind=dp_t) distdot, x(*), y(*), ddot
  external ddot
  distdot = ddot(n,x,ix,y,iy)
  return
end function distdot
