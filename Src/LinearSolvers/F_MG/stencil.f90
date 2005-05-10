module stencil_module

  use bl_types
  use bc_module
  use multifab_module

  implicit none

  integer, parameter :: ST_FILL_LAPLACE_2        = 1
  integer, parameter :: ST_FILL_LAPLACE_4        = 2
  integer, parameter :: ST_FILL_WEIGHTED_LAPLACE = 3
  integer, parameter :: ST_FILL_ABEC_LAPLACE     = 4
  integer, parameter :: ST_FILL_TENSOR_ABEC      = 5
  integer, parameter :: ST_FILL_USER_DEFINED     = 100

  integer, parameter :: ST_CROSS = 1
  integer, parameter :: ST_DENSE = 2
  integer, parameter :: ST_DIAG  = 3
  integer, parameter :: ST_TENSOR = 4

  type stencil_extrap_fab
     integer :: dim = 0
  end type stencil_extrap_fab

  type stencil_extrap
     integer :: dim = 0
     type(stencil_extrap_fab), pointer :: ext(:) => Null()
  end type stencil_extrap

  type new_stencil
     integer :: dim = 0
     integer :: type = 0
     type(multifab) :: ss
     type(imultifab) :: mm
  end type new_stencil

  type stencil
     integer :: dim = 0
     integer :: ns  = 0
     integer :: type = 0
     type(multifab)  :: ss
     type(imultifab) :: mm
     logical, pointer :: skewed(:) => Null()
     logical, pointer :: diag_0(:) => NUll()
     logical :: extrap_bc = .false.
     real(kind=dp_t), pointer :: xa(:) => Null()
     real(kind=dp_t), pointer :: xb(:) => Null()
     real(kind=dp_t), pointer :: pxa(:) => Null()
     real(kind=dp_t), pointer :: pxb(:) => Null()
     real(kind=dp_t), pointer :: dh(:) => Null()
     integer :: extrap_max_order = 0
  end type stencil

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
  real (kind = dp_t), private, parameter :: THIRD = 1.0_dp_t/3.0_dp_t
  real (kind = dp_t), private, parameter :: FOUR_THIRD = 4.0_dp_t/3.0_dp_t

  integer, parameter :: BC_GEOM = 3

  integer, dimension(3,MAX_SPACEDIM,-1:1) , parameter :: BC_BIT = &
       reshape((/ &
        0,   9,  18, &          ! BC_BIT(:,1,-1)
        1,  10,  19, &          ! BC_BIT(:,2,-1)
        2,  11,  20, &          ! BC_BIT(:,3,-1)
        3,  12,  21, &          ! BC_BIT(:,1, 0)
        4,  13,  22, &          ! BC_BIT(:,2, 0)
        5,  14,  23, &          ! BC_BIT(:,3, 0)
        6,  15,  24, &          ! BC_BIT(:,1, 1)
        7,  16,  25, &          ! BC_BIT(:,2, 1)
        8,  17,  26  &          ! BC_BIT(:,3, 1)
       /), (/3,MAX_SPACEDIM,3/))

  integer, parameter :: ST_DENSE_MAP_1D(-1:1) = &
       reshape((/ &
       1, 0, 2  &
       /), (/3/))

  integer, parameter :: ST_DENSE_MAP_2D(-1:1,-1:1) = &
       reshape((/ &
       1, 2, 3, &
       4, 0, 5, &
       6, 7, 8  &
       /), (/3,3/))

  integer, parameter :: ST_DENSE_MAP_3D(-1:1,-1:1,-1:1) = &
       reshape((/ &
       1, 2, 3, &
       4, 5, 6, &
       7, 8, 9, &
       10, 11, 12, &
       13,  0, 14, &
       15, 16, 17, &
       18, 19, 20, &
       21, 22, 23, &
       24, 25, 26  &
       /), (/3,3,3/))

! private f_bc_bit

  interface destroy
     module procedure stencil_destroy
  end interface

contains

  subroutine bc_set_bc_face(bc_face, dim, face, bc_type)
    integer, intent(inout) :: bc_face(:,:)
    integer, intent(in)    :: dim
    integer, intent(in)    :: face
    integer, intent(in)    :: bc_type
    integer :: lface

    if ( dim < 1 .or. dim > size(bc_face,dim=1)) then
       call bl_error("BC_SET_BC_FACE: out of bounds dim: ", dim)
    end if
    if ( face .ne. -1 .or. face .ne. 1 ) then
       call bl_error("BC_SET_BC_FACE: out of bounds face: ", face)
    end if
    lface = face; if ( lface == -1 ) lface = 0
    select case ( bc_type )
    case (BC_PER, BC_INT, BC_DIR, BC_NEU)
       bc_face(dim,lface) = bc_type
    case default
       call bl_error("BC_SET_BC_FACE: unknown bc_type: ", bc_type)
    end select
  end subroutine bc_set_bc_face

  function skewed_q(mm) result(r)
    logical :: r
    integer, intent(in) :: mm(:,:,:,:)
    integer :: i, j
    r = .false.
    do i = 1, 3
       do j = -1, 1
          r = r .or. any( ibits(mm(:,:,:,:),BC_BIT(BC_GEOM,i,j),1) /= 0 )
          if ( r ) return
       end do
    end do
  end function skewed_q

  function diag_0_q(ss) result(r)
    logical :: r
    real(kind=dp_t), intent(in) :: ss(:,:,:,0:)
    r = any(ss(:,:,:,0) == ZERO)
  end function diag_0_q

  elemental function bc_skewed(m, dim, face) result(r)
    logical :: r
    integer, intent(in) :: m, dim, face
    r = btest(m,BC_BIT(BC_GEOM,dim,face))
  end function bc_skewed

  elemental function bc_dirichlet(m, dim, face) result(r)
    logical :: r
    integer, intent(in) :: m, dim, face
    r = btest(m,BC_BIT(BC_DIR,dim,face))
  end function bc_dirichlet

  elemental function bc_neumann(m, dim, face) result(r)
    logical :: r
    integer, intent(in) :: m, dim, face
    r = btest(m,BC_BIT(BC_NEU,dim,face))
  end function bc_neumann

!   elemental function bc_robin(m,dim,face) result(r)
!     logical :: r
!     integer, intent(in) :: m, dim, face
!     r = btest(m,BC_BIT(BC_ROB,dim,face))
!   end function bc_robin

  elemental function bc_interior(m, dim, face) result(r)
    logical :: r
    integer, intent(in) :: m, dim, face
    r = ( ibits(m,BC_BIT(BC_DIR,dim,face),1) == 0 .and. &
          ibits(m,BC_BIT(BC_NEU,dim,face),1) == 0 )
  end function bc_interior

  elemental function bc_pretty_bit(m, dim, face) result(r)
    character :: r
    integer, intent(in) :: m, dim, face
    if ( ibits(m,BC_BIT(BC_DIR,dim,face),1) == 1 ) then
       r = 'D'
    else if ( ibits(m,BC_BIT(BC_NEU,dim,face),1) == 1 ) then
       r = 'N'
    else
       r = 'I'
    end if
  end function bc_pretty_bit

  subroutine mask_pretty_print(mf, str, unit, all, data, skip, nodal)
    use bl_IO_module
    type(imultifab), intent(in) :: mf
    character (len=*), intent(in), optional :: str
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: all, data, nodal
    integer, intent(in), optional :: skip
    integer :: i, ii
    integer :: un
    character(len=5) :: fn
    un = unit_stdout(unit)
    call unit_skip(un, skip)
    write(unit=un, fmt='("IMULTIFAB")', advance = 'no')
    if ( present(str) ) then
       write(unit=un, fmt='(": ",A)') str
    else
       write(unit=un, fmt='()')
    end if
    call unit_skip(un, skip)
    write(unit=un, fmt='(" DIM     = ",i2)') mf%dim
    call unit_skip(un, skip)
    write(unit=un, fmt='(" NC      = ",i2)') mf%nc
    call unit_skip(un, skip)
    write(unit=un, fmt='(" NG      = ",i2)') mf%ng
    call unit_skip(un, skip)
    write(unit=un, fmt='(" NODAL   = ",3(L2,1X))') mf%nodal
    call unit_skip(un, skip)
    write(unit=un, fmt='(" NBOXES  = ",i2)') mf%nboxes
    do ii = 0, parallel_nprocs()
       if ( ii == parallel_myproc() ) then
          do i = 1, mf%nboxes; if ( remote(mf,i) ) cycle
             write(unit=fn, fmt='(i5)') i
             call mask_pretty_print_fab(mf%fbs(i), str = fn, unit = unit, &
                  all = all, data = data, skip = unit_get_skip(skip) + 2, &
                  nodal = nodal)
          end do
       end if
       call parallel_barrier()
    end do
  end subroutine mask_pretty_print

  subroutine mask_pretty_print_fab(fb, str, unit, all, data, bx, skip, nodal)
    use bl_IO_module
    type(ifab), intent(in) :: fb
    character(len=*), intent(in), optional :: str
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: all, data, nodal
    integer, intent(in), optional :: skip
    type(box), intent(in), optional :: bx
    integer :: un
    logical :: lall, ldata, lnodal
    type(box) :: lbx
    lbx  = box_allbox(fb%dim); if ( present(bx) ) lbx  = bx
    lall = .TRUE.; if ( present(all) ) lall = all
    ldata = .TRUE.; if ( present(data) ) ldata = data
    lnodal = .FALSE.; if ( present(nodal) ) lnodal = nodal
    un = unit_stdout(unit)
    call unit_skip(un, skip)
    write(unit=un, fmt='("IFAB")', advance = 'no')
    if ( present(str) ) then
       write(unit=un, fmt='(": ",A)') str
    else
       write(unit=un, fmt='()')
    end if
    call unit_skip(un, skip)
    write(unit=un, fmt='(" DIM     = ",i2)') fb%dim
    call unit_skip(un, skip)
    write(unit=un, fmt='(" NG      = ",i2)') fb%ng
    call unit_skip(un, skip)
    write(unit=un, fmt='(" NC      = ",i2)') fb%nc
    call unit_skip(un, skip)
    write(unit=un, fmt='(" IBX     = ",i2)', advance = 'no')
    call print(fb%ibx, unit = un)
    call unit_skip(un, skip)
    write(unit=un, fmt='(" PBX     = ",i2)', advance = 'no')
    call print(fb%pbx, unit = un)
    call unit_skip(un, skip)
    write(unit=un, fmt='(" BX      = ",i2)', advance = 'no')
    call print(fb%bx, unit = un)
    if ( .not. associated(fb%p) ) then
       call unit_skip(un, skip)
       write(unit=un) 'NOT ASSOCIATED'
    else
       select case (fb%dim)
       case (1)
          call print_1d(fb%p(:,1,1,:), lbound(fb%p), intersection(fb%ibx,lbx), lnodal)
       case (2)
          call print_2d(fb%p(:,:,1,:), lbound(fb%p), intersection(fb%ibx,lbx), lnodal)
       case (3)
          call print_3d(fb%p(:,:,:,:), lbound(fb%p), intersection(fb%ibx,lbx), lnodal)
       end select
    end if

  contains

    subroutine print_1d(fb, lo, bx, lnodal)
      integer, intent(in) :: lo(:)
      type(box), intent(in) :: bx
      integer, intent(in) :: fb(lo(1):,:)
      logical, intent(in) :: lnodal
      integer n, i
      integer nc, hi(1)
      character(len=1) c
      nc = size(fb,dim=2)
      hi(1) = lo(1) + size(fb,dim=1) - 1
      if ( ldata ) then
         do n = 1, nc
            do i = lo(1), hi(1)
               if ( .not. ( lall .or. contains(bx, (/i/)) ) ) cycle
               c = ' '
               if ( .not. contains(bx, (/i/)) ) c = '*'
               call unit_skip(un, skip)

               if ( lnodal ) then
                  write(unit=un, fmt='(A1,1X,I3,1(1X,I5),1X,A2,A2,A2)') &
                       c, n, i, &
                       bc_pretty_bit(fb(i,n),1, 0), &
                       bc_pretty_bit(fb(i,n),1,-1), bc_pretty_bit(fb(i,n),1,1)
               else
                  write(unit=un, fmt='(A1,1X,I3,1(1X,I5),1X,A2,A2)') &
                       c, n, i, &
                       bc_pretty_bit(fb(i,n),1,-1), bc_pretty_bit(fb(i,n),1,1)
               end if


            end do
         end do
      end if
    end subroutine print_1d

    subroutine print_2d(fb, lo, bx, lnodal)
      integer, intent(in) :: lo(:)
      type(box), intent(in) :: bx
      integer, intent(in) :: fb(lo(1):,lo(2):,:)
      logical, intent(in) :: lnodal
      integer n, j, i
      integer nc, hi(2)
      character(len=1) c
      nc = size(fb,dim=3)
      do i = 1, 2
         hi(i) = lo(i) + size(fb,dim=i) - 1
      end do
      if ( ldata ) then
         do n = 1, nc
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  if ( .not. ( lall .or. contains(bx, (/i,j/)) ) ) cycle
                  c = ' '
                  if ( .not. contains(bx, (/i,j/)) ) c = '*'
                  call unit_skip(un, skip)

                  if ( lnodal) then
                     write(unit=un, fmt='(A1,1X,I3,2(1X,I5),1X,A2,A2,A2,A2,A2)') &
                          c, n, i, j, &
                          bc_pretty_bit(fb(i,j,n),1, 0), &
                          bc_pretty_bit(fb(i,j,n),1,-1), bc_pretty_bit(fb(i,j,n),1,1), &
                          bc_pretty_bit(fb(i,j,n),2,-1), bc_pretty_bit(fb(i,j,n),2,1)
                  else
                     write(unit=un, fmt='(A1,1X,I3,2(1X,I5),1X,A2,A2,A2,A2)') &
                          c, n, i, j, &
                          bc_pretty_bit(fb(i,j,n),1,-1), bc_pretty_bit(fb(i,j,n),1,1), &
                          bc_pretty_bit(fb(i,j,n),2,-1), bc_pretty_bit(fb(i,j,n),2,1)
                  end if

               end do
            end do
         end do
      end if
    end subroutine print_2d

    subroutine print_3d(fb, lo, bx, lnodal)
      integer, intent(in) :: lo(:)
      type(box), intent(in) :: bx
      integer, intent(in) :: fb(lo(1):,lo(2):,lo(3):,:)
      logical, intent(in) :: lnodal
      integer :: n, k, j, i
      integer :: nc, hi(3)
      character(len=1) :: c
      nc = size(fb,dim=4)
      do i = 1, 3
         hi(i) = lo(i) + size(fb,dim=i) - 1
      end do
      if ( ldata ) then
         do n = 1, nc
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     if ( .not. ( lall .or. contains(bx, (/i,j,k/)) ) ) cycle
                     c = ' '
                     if ( .not. contains(bx, (/i,j,k/)) ) c = '*'
                     call unit_skip(un, skip)

                     if ( lnodal) then
                        write(unit=un, fmt='(A1,1X,I3,3(1X,I5),1X,A2,A2,A2,A2,A2,A2,A2)') &
                             c, n, i, j, k, &
                             bc_pretty_bit(fb(i,j,k,n),1, 0),  &
                             bc_pretty_bit(fb(i,j,k,n),1,-1), bc_pretty_bit(fb(i,j,k,n),1,1), &
                             bc_pretty_bit(fb(i,j,k,n),2,-1), bc_pretty_bit(fb(i,j,k,n),2,1), &
                             bc_pretty_bit(fb(i,j,k,n),3,-1), bc_pretty_bit(fb(i,j,k,n),3, 1)
                     else
                        write(unit=un, fmt='(A1,1X,I3,3(1X,I5),1X,A2,A2,A2,A2,A2,A2)') &
                             c, n, i, j, k, &
                             bc_pretty_bit(fb(i,j,k,n),1,-1), bc_pretty_bit(fb(i,j,k,n),1,1), &
                             bc_pretty_bit(fb(i,j,k,n),2,-1), bc_pretty_bit(fb(i,j,k,n),2,1), &
                             bc_pretty_bit(fb(i,j,k,n),3,-1), bc_pretty_bit(fb(i,j,k,n),3, 1)
                     end if

                  end do
               end do
            end do
         end do
      end if
    end subroutine print_3d

  end subroutine mask_pretty_print_fab

! elemental function f_bc_bit(bc, dim, face) result(r)
!   integer :: r
!   integer, intent(in) :: bc
!   integer, intent(in) :: dim
!   integer, intent(in) :: face
!   r = 0 + (bc-1)*2*MAX_SPACEDIM + (dim-1) + ((face+3)/2-1)*MAX_SPACEDIM
! end function f_bc_bit

  subroutine stencil_destroy(st)
    type(stencil), intent(inout) :: st
    call destroy(st%ss)
    call destroy(st%mm)
    deallocate(st%skewed, st%xa, st%xb, st%dh, st%diag_0,st%pxa,st%pxb)
    st%dim = 0
    st%ns = -1
  end subroutine stencil_destroy

  subroutine stencil_build(st, la, dh, type, nc, nodal)
    type(stencil), intent(out) :: st
    type(layout), intent(in) :: la
    real(kind=dp_t), intent(in) :: dh(:)
    integer, intent(in) :: type
    integer, intent(in), optional :: nc
    logical, intent(in), optional :: nodal(:)
    st%dim = layout_dim(la)
    st%type = type
    if ( present(nc) ) then
       st%ns = nc
    else
       select case ( type )
       case (ST_CROSS)
          st%ns = 1 + 2*st%dim
       case (ST_DENSE)
          st%ns = 3**st%dim
       case default
          call bl_error("STENCIL_BUILD: TYPE not known: ", type)
       end select
    end if
    call multifab_build(st%ss, la, st%ns, 0, nodal = nodal)
    call imultifab_build(st%mm, la,     1, 0, nodal = nodal)
    allocate(st%skewed(nboxes(la)))
    allocate(st%diag_0(nboxes(la)))
    allocate(st%xa(st%dim), st%xb(st%dim), st%dh(st%dim),st%pxa(st%dim),st%pxb(st%dim))
    st%xa = ZERO
    st%xb = ZERO
    st%pxa = ZERO
    st%pxb = ZERO
    st%dh = dh
  end subroutine stencil_build

  function stencil_norm_st(st, mask) result(r)
    type(stencil), intent(in) :: st
    type(lmultifab), intent(in), optional :: mask
    real(kind=dp_t) :: r
    r = stencil_norm(st%ss, mask)
  end function stencil_norm_st

  function stencil_norm(ss, mask) result(r)
    real(kind=dp_t) :: r
    type(multifab), intent(in) :: ss
    type(lmultifab), intent(in), optional :: mask
    integer :: i
    real(kind=dp_t) :: r1
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    logical, pointer :: lp(:,:,:,:)
    r1 = -Huge(r1)
    if ( present(mask) ) then
       !$OMP PARALLEL DO PRIVATE(i,sp) REDUCTION(MAX:r1)
       do i = 1, ss%nboxes
          if ( remote(ss,i) ) cycle
          sp => dataptr(ss, i)
          lp => dataptr(mask, i)
          r1 = max(r1, maxval(sum(abs(sp),dim=4),mask=lp(:,:,:,1)))
       end do
       !$OMP END PARALLEL DO
    else
       !$OMP PARALLEL DO PRIVATE(i,sp) REDUCTION(MAX:r1)
       do i = 1, ss%nboxes
          if ( multifab_remote(ss,i) ) cycle
          sp => dataptr(ss, i)
          r1 = max(r1, maxval(sum(abs(sp(:,:,:,:)),dim=4)))
       end do
       !$OMP END PARALLEL DO
    end if
    call parallel_reduce(r,r1,MPI_MAX)
  end function stencil_norm

  subroutine stencil_defect_st(st, rr, ff, uu)
    type(stencil), intent(in) :: st
    type(multifab), intent(inout) :: uu, rr
    type(multifab), intent(in) :: ff
    call stencil_apply_st(st, rr, uu)
    call saxpy(rr, ff, -ONE, rr)
  end subroutine stencil_defect_st

  subroutine stencil_apply_st_c(st, rr, cr, uu, cu)
    type(stencil), intent(in) :: st
    type(multifab), intent(inout) :: rr
    type(multifab), intent(inout) :: uu
    integer, intent(in) :: cr, cu
    real(kind=dp_t), pointer :: rp(:,:,:,:)
    real(kind=dp_t), pointer :: up(:,:,:,:)
    real(kind=dp_t), pointer :: upn(:,:,:,:)
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    integer        , pointer :: mp(:,:,:,:)
    integer i
    !   integer :: nn = 21
    logical :: skwd

    call multifab_fill_boundary_c(uu, cu, 1)
    if ( st%extrap_bc) then
       ! call multifab_print(uu, unit=nn); nn = nn + 1
       ! call multifab_print(uu, 'before')
       call stencil_extrap_bc(st, uu, cu)
       ! call multifab_print(uu, 'after')
       ! call multifab_print(uu, unit=nn); nn = nn + 1
    end if

    do i = 1, rr%nboxes; if ( multifab_remote(rr, i) ) cycle
       rp => dataptr(rr, i, get_ibox(rr,i), cr)
       up => dataptr(uu, i, cu)
       upn => dataptr(uu, i, get_ibox(uu,i), cu)
       sp => dataptr(st%ss, i)
       mp => dataptr(st%mm, i)
       select case ( st%type )
       case (ST_DIAG)
          rp = sp*upn
       case (ST_CROSS)
          skwd = st%skewed(i)
          select case( st%dim )
          case (1)
             call stencil_apply_1d(sp(:,1,1,:), rp(:,1,1,1), up(:,1,1,1), mp(:,1,1,1), uu%ng, skwd)
          case (2)
             call stencil_apply_2d(sp(:,:,1,:), rp(:,:,1,1), up(:,:,1,1), mp(:,:,1,1), uu%ng, skwd)
          case (3)
             call stencil_apply_3d(sp(:,:,:,:), rp(:,:,:,1), up(:,:,:,1), mp(:,:,:,1), uu%ng, skwd)
          end select
       case (ST_DENSE)
          select case( st%dim )
          case (1)
             call stencil_dense_apply_1d(sp(:,1,1,:), rp(:,1,1,1), up(:,1,1,1), mp(:,1,1,1), uu%ng)
          case (2)
             call stencil_dense_apply_2d(sp(:,:,1,:), rp(:,:,1,1), up(:,:,1,1), mp(:,:,1,1), uu%ng)
          case (3)
             call stencil_dense_apply_3d(sp(:,:,:,:), rp(:,:,:,1), up(:,:,:,1), mp(:,:,:,1), uu%ng)
          end select
       end select
    end do
  end subroutine stencil_apply_st_c

  subroutine stencil_apply_st(st, rr, uu, c, mask)
    type(stencil), intent(in) :: st
    type(multifab), intent(inout) :: rr
    type(multifab), intent(inout) :: uu
    integer, intent(in), optional :: c
    logical, intent(in), optional :: mask(:)
    real(kind=dp_t), pointer :: rp(:,:,:,:)
    real(kind=dp_t), pointer :: up(:,:,:,:)
    real(kind=dp_t), pointer :: upn(:,:,:,:)
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    integer        , pointer :: mp(:,:,:,:)
    integer i, n
!   integer :: nn = 21
    logical :: skwd

    call multifab_fill_boundary(uu)
    if ( st%extrap_bc) then
!      call multifab_print(uu, unit=nn); nn = nn + 1
       call stencil_extrap_bc(st, uu)
!      call multifab_print(uu, unit=nn); nn = nn + 1
    end if

    if ( present(c) ) then
       n = c
       if ( present(mask) ) then
          if ( .not. mask(n) ) return
       end if
       do i = 1, rr%nboxes; if ( multifab_remote(rr, i) ) cycle
          rp => dataptr(rr, i, get_ibox(rr,i))
          up => dataptr(uu, i, n)
          upn => dataptr(uu, i, get_ibox(uu,i), n)
          sp => dataptr(st%ss, i)
          mp => dataptr(st%mm, i)
          select case ( st%type )
          case (ST_DIAG)
             rp = sp*upn
          case (ST_CROSS)
             skwd = st%skewed(i)
             select case( st%dim )
             case (1)
                call stencil_apply_1d(sp(:,1,1,:), rp(:,1,1,1), up(:,1,1,1), mp(:,1,1,1), uu%ng, skwd)
             case (2)
                call stencil_apply_2d(sp(:,:,1,:), rp(:,:,1,1), up(:,:,1,1), mp(:,:,1,1), uu%ng, skwd)
             case (3)
                call stencil_apply_3d(sp(:,:,:,:), rp(:,:,:,1), up(:,:,:,1), mp(:,:,:,1), uu%ng, skwd)
             end select
          case (ST_DENSE)
             select case( st%dim )
             case (1)
                call stencil_dense_apply_1d(sp(:,1,1,:), rp(:,1,1,1), up(:,1,1,1), mp(:,1,1,1), uu%ng)
             case (2)
                call stencil_dense_apply_2d(sp(:,:,1,:), rp(:,:,1,1), up(:,:,1,1), mp(:,:,1,1), uu%ng)
             case (3)
                call stencil_dense_apply_3d(sp(:,:,:,:), rp(:,:,:,1), up(:,:,:,1), mp(:,:,:,1), uu%ng)
             end select
          end select
       end do
    else
       do i = 1, rr%nboxes; if ( multifab_remote(rr, i) ) cycle
          sp => dataptr(st%ss, i)
          mp => dataptr(st%mm, i)
          do n = 1, rr%nc
             if ( present(mask) ) then
                if ( .not. mask(n) ) cycle
             end if
             rp => dataptr(rr, i, get_ibox(rr,i), n)
             up => dataptr(uu, i, n)
             upn => dataptr(uu, i, get_ibox(uu,i), n)
             select case ( st%type )
             case (ST_DIAG)
                rp = sp*upn
             case (ST_CROSS)
                skwd = st%skewed(i)
                select case( st%dim )
                case (1)
                   call stencil_apply_1d(sp(:,1,1,:), rp(:,1,1,1), up(:,1,1,1), mp(:,1,1,1), uu%ng, skwd)
                case (2)
                   call stencil_apply_2d(sp(:,:,1,:), rp(:,:,1,1), up(:,:,1,1), mp(:,:,1,1), uu%ng, skwd)
                case (3)
                   call stencil_apply_3d(sp(:,:,:,:), rp(:,:,:,1), up(:,:,:,1), mp(:,:,:,1), uu%ng, skwd)
                end select
             case (ST_DENSE)
                select case( st%dim )
                case (1)
                   call stencil_dense_apply_1d(sp(:,1,1,:), rp(:,1,1,1), up(:,1,1,1), mp(:,1,1,1), uu%ng)
                case (2)
                   call stencil_dense_apply_2d(sp(:,:,1,:), rp(:,:,1,1), up(:,:,1,1), mp(:,:,1,1), uu%ng)
                case (3)
                   call stencil_dense_apply_3d(sp(:,:,:,:), rp(:,:,:,1), up(:,:,:,1), mp(:,:,:,1), uu%ng)
                end select
             end select
          end do
       end do
    end if
  end subroutine stencil_apply_st

  subroutine stencil_flux_fill_st(st, br, uu, rr, c, mask)
    use bndry_reg_module
    type(stencil), intent(in) :: st
    type(bndry_reg), intent(inout) :: br
    type(multifab), intent(inout) :: uu
    integer, intent(in) :: rr(:)
    integer, intent(in), optional :: c
    logical, intent(in), optional :: mask(:)
    real(kind=dp_t), pointer :: rp(:,:,:,:)
    real(kind=dp_t), pointer :: up(:,:,:,:)
    real(kind=dp_t), pointer :: upn(:,:,:,:)
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    integer        , pointer :: mp(:,:,:,:)
    integer i, n, ff, dm, dd
    logical :: skwd
    integer :: lrr

    call multifab_fill_boundary(uu)
    if ( st%extrap_bc) then
       call stencil_extrap_bc(st, uu)
    end if

    dm = st%dim
    if ( present(c) ) then
       n = c
       if ( present(mask) ) then
          if ( .not. mask(n) ) return
       end if
       do i = 1, uu%nboxes; if ( multifab_remote(uu, i) ) cycle
          do ff = -1, 1, 2
             do dd = 1, dm
                up => dataptr(uu, i, n)
                upn => dataptr(uu, i, get_ibox(uu,i), n)
                sp => dataptr(st%ss, i)
                mp => dataptr(st%mm, i)
                select case ( st%type )
                case (ST_CROSS)
                   skwd = st%skewed(i)
                   select case( st%dim )
                   case (1)
                      call stencil_flux_1d( &
                           sp(:,1,1,:), rp(:,1,1,1), up(:,1,1,1), mp(:,1,1,1), &
                           uu%ng, lrr, ff, dd, skwd)
                   case (2)
                      call stencil_flux_2d( &
                           sp(:,:,1,:), rp(:,:,1,1), up(:,:,1,1), mp(:,:,1,1), &
                           uu%ng, lrr, ff, dd, skwd)
                   case (3)
                      call stencil_flux_3d( &
                           sp(:,:,:,:), rp(:,:,:,1), up(:,:,:,1), mp(:,:,:,1), &
                           uu%ng, lrr, ff, dd, skwd)
                   end select
                case default
                   call bl_error("STENCIL_FLUX_FILL_ST: not ready yet ", st%type)
                end select
             end do
          end do
       end do
    else
       do i = 1, uu%nboxes; if ( multifab_remote(uu, i) ) cycle
          sp => dataptr(st%ss, i)
          mp => dataptr(st%mm, i)
          do ff = -1, 1, 2
             do dd = 1, dm
                do n = 1, uu%nc
                   if ( present(mask) ) then
                      if ( .not. mask(n) ) cycle
                   end if
                   up => dataptr(uu, i, n)
                   upn => dataptr(uu, i, get_ibox(uu,i), n)
                   select case ( st%type )
                   case (ST_CROSS)
                      skwd = st%skewed(i)
                      select case( st%dim )
                      case (1)
                         call stencil_flux_1d( &
                           sp(:,1,1,:), rp(:,1,1,1), up(:,1,1,1), mp(:,1,1,1), &
                           uu%ng, lrr, ff, dd, skwd)
                      case (2)
                         call stencil_flux_2d( &
                           sp(:,:,1,:), rp(:,:,1,1), up(:,:,1,1), mp(:,:,1,1), &
                           uu%ng, lrr, ff, dd, skwd)
                      case (3)
                         call stencil_flux_3d( &
                           sp(:,:,:,:), rp(:,:,:,1), up(:,:,:,1), mp(:,:,:,1), &
                           uu%ng, lrr, ff, dd, skwd)
                      end select
                   case default
                      call bl_error("STENCIL_FLUX_FILL_ST: not ready yet ", st%type)
                   end select
                end do
             end do
          end do
       end do
    end if
  end subroutine stencil_flux_fill_st

  subroutine stencil_fill(st, pdv, bc_face, fill, coeffs, iparm, rparm, fill_fcn, alpha, beta)
    type(stencil), intent(inout)  :: st
    type(multifab), intent(inout), optional :: coeffs
    type(boxarray), intent(in) :: pdv
    integer, intent(in) :: bc_face(:,:)
    integer, intent(in) :: fill
    real(kind=dp_t), intent(in), optional :: alpha, beta
    integer, intent(in), optional :: iparm(:)
    real(kind=dp_t), intent(in), optional :: rparm(:)
    interface
       subroutine fill_fcn(dim, ss, cc, dh, mm, xa, xb, iparm, rparm)
         use bl_types
         integer,         intent(in   ) :: dim
         real(kind=dp_t), intent(inout) :: ss(:,:,:,0:)
         real(kind=dp_t), intent(in   ) :: cc(0:,0:,0:,:)
         integer,         intent(inout) :: mm(:,:,:,:)
         real(kind=dp_t), intent(in   ) :: dh(:), xa(:), xb(:)
         integer, intent(in), optional :: iparm(:)
         real(kind=dp_t), intent(in), optional :: rparm(:)
       end subroutine fill_fcn
    end interface
    optional :: fill_fcn
    integer :: i, id
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    real(kind=dp_t), pointer :: cp(:,:,:,:)
    integer, pointer :: mp(:,:,:,:)
    type(box) :: bx, pd
    real(kind=dp_t) lxa(st%dim), lxb(st%dim)
    integer :: lorder
    real(kind=dp_t) :: lal, lbe

    pd = get_pd(get_layout(st%ss))

    lal = ZERO; if ( present(alpha) ) lal = alpha
    lbe =  ONE; if ( present(beta)  ) lbe = beta

    if ( present(fill_fcn) ) then
       if ( fill /= ST_FILL_USER_DEFINED ) then
          call bl_error("STENCIL_FILL_ST: fill_fcn and fill /= ST_FILL_USER_DEFINED")
       end if
    else
       if ( fill == ST_FILL_USER_DEFINED ) then
          call bl_error("STENCIL_FILL_ST: no fill_fcn and fill == ST_FILL_USER_DEFINED")
       end if
    end if

    do i = 1, st%ss%nboxes; if ( multifab_remote(st%ss,i) ) cycle
       sp => dataptr(st%ss, i)
       select case ( fill )
       case ( ST_FILL_LAPLACE_2 )
          select case (st%dim)
          case (1)
             call s_laplacian_1d(sp(:,1,1,:), st%dh, lal, lbe)
          case (2)
             call s_laplacian_2d(sp(:,:,1,:), st%dh, lal, lbe)
          case (3)
             call s_laplacian_3d(sp(:,:,:,:), st%dh, lal, lbe)
          end select
       case ( ST_FILL_LAPLACE_4 )
          select case (st%dim)
          case (1)
             call s_mehrstellen_1d(sp(:,1,1,:), st%dh, lal, lbe)
          case (2)
             call s_mehrstellen_2d(sp(:,:,1,:), st%dh, lal, lbe)
          case (3)
             call s_mehrstellen_3d(sp(:,:,:,:), st%dh, lal, lbe)
          end select
       case ( ST_FILL_WEIGHTED_LAPLACE )
          if ( .not. present(coeffs) ) then
             call bl_error("STENCIL_FILL: ST_FILL_ABEC_LAPLACE requires coeffs")
          end if
          cp => dataptr(coeffs, i)
          select case (st%dim)
          case (1)
             call s_weighted_laplacian_1d(sp(:,1,1,:), cp(:,1,1,1), st%dh, lal, lbe)
          case (2)
             call s_weighted_laplacian_2d(sp(:,:,1,:), cp(:,:,1,1), st%dh, lal, lbe)
          case (3)
             call s_weighted_laplacian_3d(sp(:,:,:,:), cp(:,:,:,1), st%dh, lal, lbe)
          end select
       case ( ST_FILL_ABEC_LAPLACE )
          if ( .not. present(coeffs) ) then
             call bl_error("STENCIL_FILL: ST_FILL_ABEC_LAPLACE requires coeffs")
          end if
          cp => dataptr(coeffs, i)
          mp => dataptr(st%mm, i)
          lorder = 2; if ( present(iparm) ) lorder = iparm(1)
          bx = get_box(st%ss,i)
          lxa = st%xa
          lxb = st%xb
          do id = 1, pd%dim
             if ( .not. st%ss%la%lap%pmask(id) ) then
                if ( bx%lo(id) == pd%lo(id) ) then
                   lxa(id) = st%pxa(id)
                end if
                if ( bx%hi(id) == pd%hi(id) ) then
                   lxb(id) = st%pxb(id)
                end if
             end if
          end do
          select case (st%dim)
          case (1)
             call s_simple_1d_cc(sp(:,1,1,:), cp(:,1,1,:), st%dh, &
                  mp(:,1,1,1), lxa, lxb, lorder)
          case (2)
             call s_simple_2d_cc(sp(:,:,1,:), cp(:,:,1,:), st%dh, &
                  mp(:,:,1,1), lxa, lxb, lorder)
          case (3)
             call s_simple_3d_cc(sp(:,:,:,:), cp(:,:,:,:), st%dh, &
                  mp(:,:,:,1), lxa, lxb, lorder)
          end select
       case ( ST_FILL_USER_DEFINED )
          if ( .not. present(coeffs) ) then
             call bl_error("STENCIL_FILL: ST_FILL_USER_DEFINED requires coeffs")
          end if
          cp => dataptr(coeffs, i)
          mp => dataptr(st%mm, i)
          call fill_fcn(st%dim, sp, cp, st%dh, mp, st%xa, st%xb, iparm, rparm)
       end select
       st%diag_0(i) = diag_0_q(sp)
    end do
    if ( fill /= ST_FILL_USER_DEFINED .and. present(coeffs) ) then
       call plus_plus(st%ss, 1, coeffs, 1)
    end if
  end subroutine stencil_fill

  subroutine stencil_set_extrap_bc(st, max_order)
    type(stencil), intent(inout) :: st
    integer, intent(in) :: max_order
    st%extrap_bc = .true.
    st%extrap_max_order = max_order
  end subroutine stencil_set_extrap_bc

  subroutine stencil_set_bc_st(st, bc_face)
    type(stencil), intent(inout) :: st
    integer,       intent(in)    :: bc_face(:,:)
    integer,       pointer       :: mp(:,:,:,:)
    integer                      :: i
    do i = 1, st%ss%nboxes
       if ( multifab_remote(st%ss,i) ) cycle
       mp => dataptr(st%mm, i)
       call stencil_set_bc(st%ss, i, st%mm%fbs(i), bc_face)
       st%skewed(i) = skewed_q(mp)
    end do
  end subroutine stencil_set_bc_st

  subroutine stencil_print(st, str, unit, legacy, skip)
    use bl_IO_module
    type(stencil), intent(in) :: st
    character (len=*), intent(in), optional :: str
    integer, intent(in), optional :: unit
    integer, intent(in), optional :: skip
    logical, intent(in), optional :: legacy
    integer :: un
    un = unit_stdout(unit)
    if ( parallel_IOProcessor() ) then
       call unit_skip(un, skip)
       write(unit=un, fmt='("STENCIL ", i1)', advance = 'NO') 
       if ( present(str) ) then
          write(unit=un, fmt='(": ",A)') str
       else
          write(unit=un, fmt='()')
       end if
       call unit_skip(un, skip)
       write(unit=un, fmt='(" DIM     = ",i2)') st%dim
       call unit_skip(un, skip)
       write(unit=un, fmt='(" NS      = ",i2)') st%ns
       call unit_skip(un, skip)
       write(unit=un, fmt='(" TYPE    = ",i2)') st%type
       if ( st%extrap_bc) then
          call unit_skip(un, skip)
          write(unit=un, fmt='(" EXTRAP_BC")')
          call unit_skip(un, skip)
          write(unit=un, fmt='("   ORDER = ",i2)') st%extrap_max_order
       end if
       call unit_skip(un, skip)
       write(unit=un, fmt='(" SKWD    = ",i10,"/",i10  )') count(st%skewed), size(st%skewed)
       call unit_skip(un, skip)
       write(unit=un, fmt='(" XA      = ",3(ES20.10,1x))') st%xa
       call unit_skip(un, skip)
       write(unit=un, fmt='(" XB      = ",3(ES20.10,1x))') st%xb
       call unit_skip(un, skip)
       write(unit=un, fmt='(" PXA     = ",3(ES20.10,1x))') st%pxa
       call unit_skip(un, skip)
       write(unit=un, fmt='(" PXB     = ",3(ES20.10,1x))') st%pxb
       call unit_skip(un, skip)
       write(unit=un, fmt='(" DH      = ",3(ES20.10,1x))') st%dh
    end if
  end subroutine stencil_print

  subroutine stencil_extrap_bc(st, uu, c)
    type(stencil), intent(in) :: st
    type(multifab), intent(inout) :: uu
    integer, intent(in), optional :: c
    real(kind=dp_t), pointer :: up(:,:,:,:)
    integer, pointer :: mp(:,:,:,:)
    integer i, n

    if ( st%extrap_max_order < 1 ) then
       call bl_error("STENCIL_EXTRAP_BC: extrap_max_order < 1: ", st%extrap_max_order)
    end if

    if ( present(c) ) then
    do i = 1, uu%nboxes
       if ( multifab_remote(uu, i) ) cycle
       up => dataptr(uu, i)
       mp => dataptr(st%mm, i)
       n = c
          select case ( st%dim )
          case (1)
             call extrap_1d(up(:,1,1,n), mp(:,1,1,1), st%xa, st%xb, st%dh, &
                  st%extrap_max_order, st%type == ST_CROSS)
          case (2)
             call extrap_2d(up(:,:,1,n), mp(:,:,1,1), st%xa, st%xb, st%dh, &
                  st%extrap_max_order, st%type == ST_CROSS)
          case (3)
             call extrap_3d(up(:,:,:,n), mp(:,:,:,1), st%xa, st%xb, st%dh, &
                  st%extrap_max_order, st%type == ST_CROSS)
          end select
    end do
    else
    do i = 1, uu%nboxes
       if ( multifab_remote(uu, i) ) cycle
       up => dataptr(uu, i)
       mp => dataptr(st%mm, i)
       do n = 1, uu%nc
          select case ( st%dim )
          case (1)
             call extrap_1d(up(:,1,1,n), mp(:,1,1,1), st%xa, st%xb, st%dh, &
                  st%extrap_max_order, st%type == ST_CROSS)
          case (2)
             call extrap_2d(up(:,:,1,n), mp(:,:,1,1), st%xa, st%xb, st%dh, &
                  st%extrap_max_order, st%type == ST_CROSS)
          case (3)
             call extrap_3d(up(:,:,:,n), mp(:,:,:,1), st%xa, st%xb, st%dh, &
                  st%extrap_max_order, st%type == ST_CROSS)
          end select
       end do
    end do
    end if
  end subroutine stencil_extrap_bc

  subroutine extrap_1d(ph, mm, xa, xb, dh, max_order, cross)
    integer, intent(in) :: mm(:)
    real(kind=dp_t), intent(inout) :: ph(0:)
    real(kind=dp_t), intent(in) :: xa(:), xb(:), dh(:)
    integer, intent(in) :: max_order
    logical, intent(in) :: cross
    integer i, nn(1)
    integer :: norder(1)
    real(dp_t) :: xx(0:max_order,2)
    real(dp_t) :: cc(0:max_order,3)
    integer, save :: saved_max_order = -1

    nn = ubound(mm)
    if ( max_order < 1 ) then
       call bl_error("EXTRAP_1D: max_order < 1: ", max_order)
    end if

    xx(:,1) = (/ Huge(xx), (i+HALF, i=0, saved_max_order-2) /)
    xx(:,2) = (/ Huge(xx), (i+HALF, i=0, saved_max_order-2) /)

    norder = min(nn, max_order)

    xx(0,1) = -xa(1)/dh(1)
    xx(0,2) = -xb(1)/dh(1)
    call poly_interp_coeff(cc(0:norder(1),1), -HALF, xx(0:norder(1),1))
    call poly_interp_coeff(cc(0:norder(1),2), -HALF, xx(0:norder(1),2))
    if ( bc_dirichlet(mm( 1), 1, -1) ) then
       ph(0) = sum(ph(0:norder(1))*cc(0:norder(1),1))
    end if
    if ( bc_dirichlet(mm(nn(1)), 1, +1) ) then
       ph(nn(1)+1) = sum(ph(nn(1)+1:nn(1)-(norder(1)-1):-1)*cc(0:norder(1),2))
    end if

  end subroutine extrap_1d

  subroutine extrap_2d(ph, mm, xa, xb, dh, max_order, cross)
    integer, intent(in) :: mm(:,:)
    real(kind=dp_t), intent(inout) :: ph(0:,0:)
    real(kind=dp_t), intent(in) :: xa(:), xb(:), dh(:)
    integer, intent(in) :: max_order
    logical, intent(in) :: cross
    integer i, j, nn(2)
    integer :: norder(2)
    real(dp_t) :: xx(0:max_order, 2)
    real(dp_t) :: cc(0:max_order, 2)
    real(dp_t) :: xc(1:max_order+1)

    nn = ubound(mm)
    if ( max_order < 1 ) then
       call bl_error("EXTRAP_2D: max_order < 1: ", max_order)
    end if

    xx(:,1) = (/ Huge(xx), (i+HALF, i=0, max_order-1) /)
    xx(:,2) = (/ Huge(xx), (i+HALF, i=0, max_order-1) /)

    norder = min(nn, max_order)

    xx(0,1) = -xa(1)/dh(1)
    xx(0,2) = -xb(1)/dh(1)
    call poly_interp_coeff(cc(0:norder(1),1), -HALF, xx(0:norder(1),1))
    call poly_interp_coeff(cc(0:norder(1),2), -HALF, xx(0:norder(1),2))
    do j = 1, nn(2)
       if ( bc_dirichlet(mm( 1,j), 1, -1) ) then
          ph(0,j) = sum(ph(0:norder(1),j)*cc(0:norder(1),1))
       end if
       if ( bc_dirichlet(mm(nn(1),j), 1, +1) ) then
          ph(nn(1)+1,j) = sum(ph(nn(1)+1:nn(1)-(norder(1)-1):-1,j)*cc(0:norder(1),2))
       else
          ph(nn(1)+1,j) = ph(nn(1),j)
       end if
    end do

    xx(0,1) = -xa(2)/dh(2)
    xx(0,2) = -xb(2)/dh(2)
    call poly_interp_coeff(cc(0:norder(2),1), -HALF, xx(0:norder(2),1))
    call poly_interp_coeff(cc(0:norder(2),2), -HALF, xx(0:norder(2),2))
    do i = 1, nn(1)
       if ( bc_dirichlet(mm(i, 1), 2, -1) ) then
          ph(i,0) = sum(ph(i,0:norder(2))*cc(0:norder(2),1))
       end if
       if ( bc_dirichlet(mm(i,nn(2)), 2, +1) ) then
          ph(i,nn(2)+1) = sum(ph(i,nn(2)+1:nn(2)-(norder(2)-1):-1)*cc(0:norder(2),2))
       end if
    end do


    if ( .not. cross ) then
       norder = min(nn, max_order+1)
       xc(:)  = (/(i+HALF, i=0, max_order)/)
       !! Corners
       do i = 1, 2
          call poly_interp_coeff(cc(0:norder(i)-1,i), -HALF, xc(1:norder(i)))
       end do

       if ( bc_dirichlet(mm(1,1),1,-1) .and. bc_dirichlet(mm(1,1),2,-1) ) then
          ph(0,0) = ( &
               + sum(ph(1:norder(1),0)*cc(0:norder(1)-1,1)) &
               + sum(ph(0,1:norder(2))*cc(0:norder(2)-1,2)) &
               ) * HALF
       end if
       if ( bc_dirichlet(mm(1,nn(2)),1,-1) .and. bc_dirichlet(mm(1,nn(2)),2,+1) ) then
          ph(0,nn(2)+1) = ( &
               + sum(ph(1:norder(1),nn(2)+1)           *cc(0:norder(1)-1,1)) &
               + sum(ph(0,nn(2):nn(2)-(norder(2)-1):-1)*cc(0:norder(2)-1,2)) &
               ) * HALF
       end if
       if ( bc_dirichlet(mm(nn(1),1),1,+1) .and. bc_dirichlet(mm(nn(1),1),2,-1) ) then
          ph(nn(1)+1,0) = ( &
               + sum(ph(nn(1):nn(1)-(norder(1)-1):-1,0)*cc(0:norder(1)-1,1)) &
               + sum(ph(nn(1)+1,1:norder(2))           *cc(0:norder(2)-1,2)) &
               ) * HALF
       end if
       if ( bc_dirichlet(mm(nn(1),nn(2)),1,+1) .and. bc_dirichlet(mm(nn(1),nn(2)),2,+1) ) then
          ph(nn(1)+1,nn(2)+1) = ( &
               + sum(ph(nn(1):nn(1)-(norder(1)-1):-1,nn(2)+1)*cc(0:norder(1)-1,1)) &
               + sum(ph(nn(1)+1,nn(2):nn(2)-(norder(2)-1):-1)*cc(0:norder(2)-1,2)) &
               ) * HALF
       end if
    end if

  end subroutine extrap_2d

  subroutine extrap_3d(ph, mm, xa, xb, dh, max_order, cross)
    integer, intent(in) :: mm(:,:,:)
    real(kind=dp_t), intent(inout) :: ph(0:,0:,0:)
    real(kind=dp_t), intent(in) :: xa(:), xb(:), dh(:)
    integer, intent(in) :: max_order
    logical, intent(in) :: cross
    integer i, j, k, nn(3)
    integer :: norder(3)
    real(dp_t) :: xx(0:max_order,2)
    real(dp_t) :: cc(0:max_order,3)
    real(dp_t) :: xc(1:max_order+1)

    nn = ubound(mm)

    if ( max_order < 1 ) then
       call bl_error("EXTRAP_3D: max_order < 1: ", max_order)
    end if

    xx(:,1) = (/ Huge(xx), (i+HALF, i=0, max_order-1) /)
    xx(:,2) = (/ Huge(xx), (i+HALF, i=0, max_order-1) /)

    norder = min(nn, max_order)

    xx(0,1) = -xa(1)/dh(1)
    xx(0,2) = -xb(1)/dh(1)
    call poly_interp_coeff(cc(0:norder(1),1), -HALF, xx(0:norder(1),1))
    call poly_interp_coeff(cc(0:norder(1),2), -HALF, xx(0:norder(1),2))
    do k = 1, nn(3)
       do j = 1, nn(2)
          if ( bc_dirichlet(mm(1,j,k),1,-1) ) then
             ph(0,j,k)       = sum(ph(0:norder(1),j,k)                   *cc(0:norder(1),1))
          end if
          if ( bc_dirichlet(mm(nn(1),j,k),1,+1) ) then
             ph(nn(1)+1,j,k) = sum(ph(nn(1)+1:nn(1)-(norder(1)-1):-1,j,k)*cc(0:norder(1),2))
          end if
       end do
    end do

    xx(0,1) = -xa(2)/dh(2)
    xx(0,2) = -xb(2)/dh(2)
    call poly_interp_coeff(cc(0:norder(2),1), -HALF, xx(0:norder(2),1))
    call poly_interp_coeff(cc(0:norder(2),2), -HALF, xx(0:norder(2),2))
    do k = 1, nn(3)
       do i = 1, nn(1)
          if ( bc_dirichlet(mm(i, 1,k), 2, -1) ) then
             ph(i,0,k)       = sum(ph(i,0:norder(2),k)                   *cc(0:norder(2),1))
          end if
          if ( bc_dirichlet(mm(i,nn(2),k), 2, +1) ) then
             ph(i,nn(2)+1,k) = sum(ph(i,nn(2)+1:nn(2)-(norder(2)-1):-1,k)*cc(0:norder(2),2))
          end if
       end do
    end do

    xx(0,1) = -xa(3)/dh(3)
    xx(0,2) = -xb(3)/dh(3)
    call poly_interp_coeff(cc(0:norder(3),1), -HALF, xx(0:norder(3),1))
    call poly_interp_coeff(cc(0:norder(3),2), -HALF, xx(0:norder(3),2))
    do j = 1, nn(2)
       do i = 1, nn(1)
          if ( bc_dirichlet(mm(i,j, 1), 3, -1) ) then
             ph(i,j,0)       = sum(ph(i,j,0:norder(3))                   *cc(0:norder(3),1))
          end if
          if ( bc_dirichlet(mm(i,j,nn(3)), 2, +1) ) then
             ph(i,j,nn(3)+1) = sum(ph(i,j,nn(3)+1:nn(3)-(norder(3)-1):-1)*cc(0:norder(3),2))
          end if
       end do
    end do

    if ( .not. cross ) then
       norder = min(nn, max_order+1)
       xc(:)  = (/(i+HALF, i=0, max_order)/)
       do i = 1, 3
          call poly_interp_coeff(cc(0:norder(i)-1,1), -HALF, xc(1:norder(i)))
       end do

       !! Corners
       if ( bc_dirichlet(mm(1,1,1),1,-1) .and. &
            bc_dirichlet(mm(1,1,1),2,-1) .and. &
            bc_dirichlet(mm(1,1,1),3,-1) ) then
          ph(0,0,0) = ( &
               + sum(ph(1:norder(1),0,0)*cc(0:norder(1)-1,1)) &
               + sum(ph(0,1:norder(2),0)*cc(0:norder(2)-1,2)) &
               + sum(ph(0,0,1:norder(3))*cc(0:norder(3)-1,3)) &
               ) * THIRD
       end if
       if ( bc_dirichlet(mm(1,nn(2),1),1,-1) .and. &
            bc_dirichlet(mm(1,nn(2),1),2,+1) .and. &
            bc_dirichlet(mm(1,nn(2),1),3,-1) ) then
          ph(0,nn(2)+1,0) = ( &
               + sum(ph(1:norder(1),nn(2)+1,0)           *cc(0:norder(1)-1,1)) &
               + sum(ph(0,nn(2):nn(2)-(norder(2)-1):-1,0)*cc(0:norder(2)-1,2)) &
               + sum(ph(0,nn(2)+1,1:norder(3))           *cc(0:norder(3)-1,3)) &
               ) * THIRD
       end if
       if ( bc_dirichlet(mm(nn(1),1,1),1,+1) .and. &
            bc_dirichlet(mm(nn(1),1,1),2,-1) .and. &
            bc_dirichlet(mm(nn(1),1,1),3,-1) ) then
          ph(nn(1)+1,0,0) = ( &
               + sum(ph(nn(1):nn(1)-(norder(1)-1):-1,0,0)*cc(0:norder(1)-1,1)) &
               + sum(ph(nn(1)+1,1:norder(2),0)           *cc(0:norder(2)-1,2)) &
               + sum(ph(nn(1)+1,0,1:norder(3))           *cc(0:norder(3)-1,3)) &
               ) * THIRD
       end if
       if ( bc_dirichlet(mm(nn(1),nn(2),1),1,+1) .and. &
            bc_dirichlet(mm(nn(1),nn(2),1),2,+1) .and. &
            bc_dirichlet(mm(nn(1),nn(2),1),3,-1) ) then
          ph(nn(1)+1,nn(2)+1,0) = ( &
               + sum(ph(nn(1):nn(1)-(norder(1)-1):-1,nn(2)+1,0)*cc(0:norder(1)-1,1)) &
               + sum(ph(nn(1)+1,nn(2):nn(2)-(norder(2)-1):-1,0)*cc(0:norder(2)-1,2)) &
               + sum(ph(nn(1)+1,nn(2)+1,1:norder(3))           *cc(0:norder(3)-1,3)) &
               ) * THIRD
       end if
       if ( bc_dirichlet(mm(1,1,nn(3)),1,-1) .and. &
            bc_dirichlet(mm(1,1,nn(3)),2,-1) .and. &
            bc_dirichlet(mm(1,1,nn(3)),3,+1) ) then
          ph(0,0,nn(3)+1) = ( &
               + sum(ph(1:norder(1),0,0)                 *cc(0:norder(1)-1,1)) &
               + sum(ph(0,1:norder(2),0)                 *cc(0:norder(2)-1,2)) &
               + sum(ph(0,0,nn(3):nn(3)-(norder(3)-1):-1)*cc(0:norder(3)-1,3)) &
               ) * THIRD
       end if
       if ( bc_dirichlet(mm(1,nn(2),nn(3)),1,-1) .and. &
            bc_dirichlet(mm(1,nn(2),nn(3)),2,+1) .and. &
            bc_dirichlet(mm(1,nn(2),nn(3)),3,+1) ) then
          ph(0,nn(2)+1,nn(3)+1) = ( &
               + sum(ph(1:norder(1),nn(2)+1,0)                 *cc(0:norder(1)-1,1)) &
               + sum(ph(0,nn(2):nn(2)-(norder(2)-1):-1,0)      *cc(0:norder(2)-1,2)) &
               + sum(ph(0,nn(2)+1,nn(3):nn(3)-(norder(3)-1):-1)*cc(0:norder(3)-1,3)) &
               ) * THIRD
       end if
       if ( bc_dirichlet(mm(nn(1),1,nn(3)),1,+1) .and. &
            bc_dirichlet(mm(nn(1),1,nn(3)),2,-1) .and. &
            bc_dirichlet(mm(nn(1),1,nn(3)),3,+1) ) then
          ph(nn(1)+1,0,nn(3)+1) = ( &
               + sum(ph(nn(1):nn(1)-(norder(1)-1):-1,0,0)      *cc(0:norder(1)-1,1)) &
               + sum(ph(nn(1)+1,1:norder(2),0)                 *cc(0:norder(2)-1,2)) &
               + sum(ph(nn(1)+1,0,nn(3):nn(3)-(norder(3)-1):-1)*cc(0:norder(3)-1,3)) &
               ) * THIRD
       end if
       if ( bc_dirichlet(mm(nn(1),nn(2),nn(3)),1,+1) .and. &
            bc_dirichlet(mm(nn(1),nn(2),nn(3)),2,+1) .and. &
            bc_dirichlet(mm(nn(1),nn(2),nn(3)),3,+1) ) then
          ph(nn(1)+1,nn(2)+1,nn(3)+1) = ( &
               + sum(ph(nn(1):nn(1)-(norder(1)-1):-1,nn(2)+1,0)      *cc(0:norder(1)-1,1)) &
               + sum(ph(nn(1)+1,nn(2):nn(2)-(norder(2)-1):-1,0)      *cc(0:norder(2)-1,2)) &
               + sum(ph(nn(1)+1,nn(2)+1,nn(3):nn(3)-(norder(3)-1):-1)*cc(0:norder(3)-1,3)) &
               ) * THIRD
       end if
       !! X edges
       do i = 1, nn(1)
          if ( bc_dirichlet(mm(i,1,1),2,-1) .and. &
               bc_dirichlet(mm(i,1,1),3,-1) ) then
             ph(i,0,0) = ( &
                  + sum(ph(i,1:norder(2),0)*cc(0:norder(2)-1,2)) &
                  + sum(ph(i,0,1:norder(3))*cc(0:norder(3)-1,3)) &
                  ) * TWO
          end if
          if ( bc_dirichlet(mm(i,nn(2),1),2,+1) .and. &
               bc_dirichlet(mm(i,nn(2),1),3,-1) ) then
             ph(i,nn(2)+1,0) = ( &
                  + sum(ph(i,nn(2):nn(2)-(norder(2)-1):-1,0)*cc(0:norder(2)-1,2)) &
                  + sum(ph(i,nn(2)+1,1:norder(3))           *cc(0:norder(3)-1,3)) &
                  ) * TWO
          end if
          if ( bc_dirichlet(mm(i,1,nn(3)),2,-1) .and. &
               bc_dirichlet(mm(i,1,nn(3)),3,+1) ) then
             ph(i,0,nn(3)+1) = ( &
                  + sum(ph(i,1:norder(2),nn(3)+1)           *cc(0:norder(2)-1,2)) &
                  + sum(ph(i,0,nn(3):nn(3)-(norder(3)-1):-1)*cc(0:norder(3)-1,3)) &
                  ) * TWO
          end if
          if ( bc_dirichlet(mm(i,nn(2),nn(3)),2,+1) .and. &
               bc_dirichlet(mm(i,nn(2),nn(3)),3,+1) ) then
             ph(i,nn(2)+1,nn(3)+1) = ( &
                  + sum(ph(i,nn(2):nn(2)-(norder(2)-1):-1,nn(3)+1)*cc(0:norder(2)-1,2)) &
                  + sum(ph(i,nn(2)+1,nn(3):nn(3)-(norder(3)-1):-1)*cc(0:norder(3)-1,3)) &
                  ) * TWO
          end if
       end do
       !! Y edges
       do j = 1, nn(2)
          if ( bc_dirichlet(mm(1,j,1),1,-1) .and. &
               bc_dirichlet(mm(1,j,1),3,-1) ) then
             ph(0,j,0) = ( &
                  + sum(ph(1:norder(1),j,0)*cc(0:norder(1)-1,1)) &
                  + sum(ph(0,j,1:norder(3))*cc(0:norder(3)-1,3)) &
                  ) * TWO
          end if
          if ( bc_dirichlet(mm(nn(1),j,1),1,+1) .and. &
               bc_dirichlet(mm(nn(1),j,1),3,-1) ) then
             ph(nn(1)+1,j,0) = ( &
                  + sum(ph(nn(1):nn(1)-(norder(1)-1):-1,j,0)*cc(0:norder(1)-1,1)) &
                  + sum(ph(nn(1)+1,j,1:norder(3))           *cc(0:norder(3)-1,3)) &
                  ) * TWO
          end if
          if ( bc_dirichlet(mm(1,j,nn(3)),1,-1) .and. &
               bc_dirichlet(mm(1,j,nn(3)),3,+1) ) then
             ph(0,j,nn(3)+1) = ( &
                  + sum(ph(1:norder(1),j,nn(3)+1)           *cc(0:norder(1)-1,1)) &
                  + sum(ph(0,j,nn(3):nn(3)-(norder(3)-1):-1)*cc(0:norder(3)-1,3)) &
                  ) * TWO
          end if
          if ( bc_dirichlet(mm(nn(1),j,nn(3)),1,+1) .and. &
               bc_dirichlet(mm(nn(1),j,nn(3)),3,+1) ) then
             ph(nn(1)+1,j,nn(3)+1) = ( &
                  + sum(ph(nn(1):nn(1)-(norder(1)-1):-1,j,nn(3)+1)*cc(0:norder(1)-1,1)) &
                  + sum(ph(nn(1)+1,j,nn(3):nn(3)-(norder(3)-1):-1)*cc(0:norder(3)-1,3)) &
                  ) * TWO
          end if
       end do
       !! Z edges
       do k = 1, nn(3)
          if ( bc_dirichlet(mm(1,1,k),1,-1) .and. &
               bc_dirichlet(mm(1,1,k),2,-1) ) then
             ph(0,0,k) = ( &
                  + sum(ph(1:norder(1),0,k)*cc(0:norder(1)-1,1)) &
                  + sum(ph(0,1:norder(2),k)*cc(0:norder(2)-1,2)) &
                  ) * TWO
          end if
          if ( bc_dirichlet(mm(nn(1),1,k),1,+1) .and. &
               bc_dirichlet(mm(nn(1),1,k),2,-1) ) then
             ph(nn(1)+1,0,k) = ( &
                  + sum(ph(nn(1):nn(1)-(norder(1)-1):-1,0,k)*cc(0:norder(1)-1,1)) &
                  + sum(ph(nn(1)+1,1:norder(2),k)           *cc(0:norder(2)-1,2)) &
                  ) * TWO
          end if
          if ( bc_dirichlet(mm(1,nn(2),k),1,-1) .and. &
               bc_dirichlet(mm(1,nn(2),k),2,+1) ) then
             ph(0,nn(2)+1,k) = ( &
                  + sum(ph(1:norder(1),nn(2)+1,k)           *cc(0:norder(1)-1,1)) &
                  + sum(ph(0,nn(2):nn(2)-(norder(2)-1):-1,k)*cc(0:norder(2)-1,2)) &
                  ) * TWO
          end if
          if ( bc_dirichlet(mm(nn(1),nn(2),k),1,+1) .and. &
               bc_dirichlet(mm(nn(1),nn(2),k),2,+1) ) then
             ph(nn(1)+1,nn(2)+1,k) = ( &
                  + sum(ph(nn(1):nn(1)-(norder(1)-1):-1,nn(2)+1,k)*cc(0:norder(1)-1,1)) &
                  + sum(ph(nn(1)+1,nn(2):nn(2)-(norder(2)-1):-1,k)*cc(0:norder(2)-1,2)) &
                  ) * TWO
          end if
       end do
    end if

  end subroutine extrap_3d

  subroutine stencil_set_bc(st, idx, mask, bc_face, cf_face)
    type(multifab), intent(in)           :: st
    integer,        intent(in)           :: idx
    type(ifab),     intent(inout)        :: mask
    integer,        intent(in)           :: bc_face(:,:)
    integer,        intent(in), optional :: cf_face(:,:)

    type(box)        :: bx1, src
    type(boxarray)   :: ba, sba
    integer          :: i, j, ii, jj, k, ldom
    integer, pointer :: mp(:,:,:,:)
    integer          :: lcf_face(size(bc_face, 1), size(bc_face, 2))
    !
    ! The Coarse-Fine boundary is Dirichlet unless specified.
    !
    lcf_face = BC_DIR; if ( present(cf_face) ) lcf_face = cf_face
    !
    ! Initialize every border to Fine-Fine (Interior).
    !
    call setval(mask, BC_INT)

    do i = 1, st%dim
       if ( bc_face(i,1) == BC_PER .and. ( bc_face(i,1) /= bc_face(i,2) )) then
          call bl_error("STENCIL_SET_BC: confusion in bc_face")
       end if
       do j = -1, 1, 2
          bx1 = shift(get_box(st, idx), j, i)
          jj = (3 + j)/2
          if ( contains(st%la%lap%pd, bx1) ) then
             !
             ! We're not touching a physical boundary -- set any/all C-F bndrys.
             !
             call boxarray_boxarray_diff(ba, bx1, st%la%lap%bxa)
             do ii = 1, ba%nboxes
                bx1 = shift(ba%bxs(ii), -j, i)
                mp => dataptr(mask, bx1)
                mp = ibset(mp, BC_BIT(lcf_face(i, jj), i, j))
             end do
             call destroy(ba)
          else
             !
             ! We touch a physical boundary in that direction.
             !
             if ( .not. st%la%lap%pmask(i) ) then
                !
                ! We're not periodic in that direction -- use physical BCs.
                !
                call boxarray_box_diff(ba, bx1, st%la%lap%pd)
                do ii = 1, ba%nboxes
                   bx1 = shift(ba%bxs(ii), -j, i)
                   mp => dataptr(mask, bx1)
                   mp = ibset(mp, BC_BIT(bc_face(i, jj), i, j))
                end do
                call destroy(ba)
             else
                !
                ! Remove any/all Fine-Fine intersections.
                !
                ldom = extent(st%la%lap%pd, i)
                call boxarray_build_bx(ba, bx1)
                do k = 1, st%nboxes
                   src = shift(get_box(st, k), j*ldom, i)
                   if ( intersects(bx1, src) ) then
                      call boxarray_build_bx(sba, src)
                      call boxarray_diff(ba, sba)
                      call destroy(sba)
                   end if
                end do
                !
                ! Set any remaining boxes to C-F.
                !
                do ii = 1, ba%nboxes
                   bx1 = shift(ba%bxs(ii), -j, i)
                   mp => dataptr(mask, bx1)
                   mp = ibset(mp, BC_BIT(lcf_face(i, jj), i, j))
                end do
                call destroy(ba)
             end if
          end if
       end do
    end do

  end subroutine stencil_set_bc

  subroutine stencil_fill_cc(ss, coeffs, dh, pdv, mask, xa, xb, pxa, pxb, pd, order, bc_face, fnc)
    type(multifab), intent(inout) :: ss
    type(multifab), intent(in   ) :: coeffs
    real(kind=dp_t), intent(in) :: dh(:)
    type(boxarray), intent(in) :: pdv
    type(box), intent(in) :: pd
    type(imultifab), intent(inout) :: mask
    integer, intent(in) :: order
    integer, intent(in) :: bc_face(:,:)
    real(kind=dp_t), intent(in) :: xa(:), xb(:), pxa(:), pxb(:)
    interface
       function fnc(i, j, k, n) result(r)
         integer, intent(in) :: i, j, k, n
         integer :: r
       end function fnc
    end interface
    optional :: fnc
    type(box) :: bx
    real(kind=dp_t) lxa(ss%dim), lxb(ss%dim)

    real(kind=dp_t), pointer :: sp(:,:,:,:)
    real(kind=dp_t), pointer :: cp(:,:,:,:)
    integer, pointer :: mp(:,:,:,:)
    integer i,id

    ! Store the incoming values.

    do i = 1, ss%nboxes
       if ( multifab_remote(ss,i) ) cycle
       bx = get_box(ss,i)
       call stencil_set_bc(ss, i, mask%fbs(i), bc_face)
       lxa = xa
       lxb = xb
       do id = 1,pd%dim
          if ( .not. ss%la%lap%pmask(id) ) then
             if ( bx%lo(id) == pd%lo(id) ) then
                lxa(id) = pxa(id)
             end if
             if ( bx%hi(id) == pd%hi(id) ) then
                lxb(id) = pxb(id)
             end if
          end if
       end do

       sp => dataptr(ss, i)
       cp => dataptr(coeffs, i)
       mp => dataptr(mask, i)

       select case (ss%dim)
       case (1)
          call s_simple_1d_cc(sp(:,1,1,:), cp(:,1,1,:), dh, mp(:,1,1,1), lxa, lxb, order)
       case (2)
          call s_simple_2d_cc(sp(:,:,1,:), cp(:,:,1,:), dh, mp(:,:,1,1), lxa, lxb, order)
       case (3)
          call s_simple_3d_cc(sp(:,:,:,:), cp(:,:,:,:), dh, mp(:,:,:,1), lxa, lxb, order)
       end select
    end do

!    call print(ss,unit=22)
!    call mask_pretty_print(mask,unit=22)

  end subroutine stencil_fill_cc

  elemental function stencil_bc_type(mask, dir, face) result(r)
    integer, intent(in) :: mask, dir, face
    integer :: r
    r = BC_INT
    if      ( bc_dirichlet(mask,dir,face) ) then
       r = BC_DIR
    else if ( bc_neumann  (mask,dir,face) ) then
       r = bc_NEU
!     else if ( bc_robin    (mask,dir,face) ) then
!        r = bc_ROB
    end if
  end function stencil_bc_type
    
  subroutine stencil_bndry_aaa(maxo, nx, dir, face, mask, &
       d_s0, d_sp, d_sm, d_ss, &
       d_b0, d_b1, d_xa, d_xb, dh, d_bclo, d_bchi)
    integer, intent(in) :: maxo
    integer, intent(in) :: nx, face, dir
    integer, intent(inout) :: mask
    real(kind=dp_t), intent(inout) :: d_s0, d_sm, d_sp, d_ss
    real(kind=dp_t), intent(in) :: d_xa, d_xb, dh
    real(kind=dp_t), intent(in) :: d_b0, d_b1
    integer, intent(in) :: d_bclo, d_bchi
    real(kind=dp_t) :: f1 
    real(kind=dp_t) :: xa, xb, s0, sm, sp, ss, b0, b1
    integer :: bclo, bchi
    logical :: skewed
    integer :: imaxo
    logical :: old_old = .TRUE.

    f1 = ONE/dh**2
    skewed = .FALSE.

    if ( face == 1 ) then
       xa  = d_xb/dh
       xb  = d_xa/dh
       b0  = d_b1
       b1  = d_b0
       bclo = d_bchi
       bchi = d_bclo
    else if ( face == -1 ) then
       xa  = d_xa/dh
       xb  = d_xb/dh
       b0  = d_b0
       b1  = d_b1
       bclo = d_bclo
       bchi = d_bchi
    else 
       call bl_error("STENCIL_BNDRY_AAA: face not -1 or 1")
    end if

    !     if ( bclo == BC_ROB .and. (.not.present(aa1) .and. .not.present(bb1)) ) &
    !          call bl_error("ROBIN BC's not ready yet")
    !     if ( bchi == BC_ROB .and. (.not.present(aa2) .and. .not.present(bb2)) ) &
    !          call bl_error("ROBIN BC's not ready yet")
    if ( nx == 1 .and. face == 1 ) call bl_error("STENCIL_BNDRY_AAA: Shouldn't happen!")

    s0 = 0
    ss = 0
    sm = 0
    sp = 0
    !
    ! TODO -- this stuff is just not quite right.
    ! Some of this logic needs to be moved into the bc_?? routines themselves.
    !
    if ( nx > 1 ) bchi = BC_INT
    imaxo = maxo
    if ( nx == 1 ) imaxo = 1
    if ( nx == 2 ) imaxo = 2

    select case ( bclo ) 
    case ( BC_INT )
       select case (bchi)
       case (BC_INT)
          call bc_ii
       case (BC_DIR)
          call bc_id
       case (BC_NEU)
          call bc_in
       case default
          call bl_error("STENCIL_BNDRY_AAA: Strange BCHI ", bchi)
       end select
    case (BC_DIR)
       select case (bchi)
       case (BC_INT)
          call bc_di
       case (BC_DIR)
          call bc_dd
       case (BC_NEU)
          call bc_dn
       case default
          call bl_error("STENCIL_BNDRY_AAA: Strange BCHI ", bchi)
       end select
    case (BC_NEU)
       select case (bchi)
       case (BC_INT)
          call bc_ni
       case (BC_DIR)
          call bc_nd
       case (BC_NEU)
          call bc_nn
       case default
          call bl_error("STENCIL_BNDRY_AAA: Strange BCHI ", bchi)
       end select
    case default
       call bl_error("STENCIL_BNDRY_AAA: Strange BCLO ", bclo)
    end select
    ! One Sided Neumann

    d_s0 = d_s0 - s0*f1
    d_ss = - ss*f1
    if ( face == 1 ) then
       d_sm = - sp*f1
       d_sp = - sm*f1
       if ( skewed ) &
            mask = ibset(mask, BC_BIT(BC_GEOM,dir,-1))
    else if ( face == -1 ) then
       d_sm = - sm*f1
       d_sp = - sp*f1
       if ( skewed ) &
            mask = ibset(mask, BC_BIT(BC_GEOM,dir,+1))
    else 
       call bl_error("STENCIL_BNDRY_AAA: face not -1 or 1")
    end if
  contains

    subroutine bc_ii
      sm  = b0
      s0  = -(b0+b1)
      sp  = b1
      ss  = 0
      skewed = .false.
    end subroutine bc_ii

    subroutine bc_id
      if ( nx > 1 ) then
         call bc_ii
      else
         sm =  b0 + ( -1 + 4/(3 + 2*xb)) * b1
         s0 = -b0 + (( -3 + 2*xb )/(1 + 2*xb)) * b1
         sp =  8*b1/(3 + 4*xb*(2 + xb))
         ss = 0
         skewed = .false.
      end if
    end subroutine bc_id

    subroutine bc_in
      if ( nx > 1 ) then
         call bc_ii
      else
         sm =  b0 - xb*b1/(1 + xb)
         s0 = -b0 + xb*b1/(1 + xb)
         sp =  b1/(1 + xb)
         ss = 0
         skewed = .false.
      end if
    end subroutine bc_in

    subroutine bc_di
      select case (imaxo)
      case (1)
         sm = 2*b0/(1 + 2*xa)
         s0 = -2*b0/(1 + 2*xa) - b1
         sp = b1
         ss = 0
         skewed = .false.
      case (2)
         sm = 8*b0/(3 + 4*xa*(2 + xa))
         s0 = ((-3 + 2*xa)/(1 + 2*xa))*b0 - b1
         sp = ((1-2*xa)/(3 + 2*xa))*b0    + b1
         ss = 0
         skewed = .false.
      case(3)
         if ( old_old ) then
            sm = 48*b0/(15 + 46*xa + 36*xa**2 + 8*xa**3)
            s0 = 4*((-1 + xa)/(1 + 2*xa))*b0 -  b1
            sp = 3*((1-2*xa)/(3 + 2*xa))*b0 + b1
            ss = (-1 + 2*xa)*b0/(5 + 2*xa)
         else
            sm = 46*b0/((1 + 2*xa)*(3 + 2*xa)*(5+2*xa))
            s0 = -((15 - 16*xa)*b0 + (4 + 8*xa)*b1)/(4*(1 + 2*xa))
            sp = ((5 - 12*xa)*b0 + (6 + 4*xa)*b1)/(2*(3 + 2*xa))
            ss = (-3 + 8*xa)*b0/(4*( 5 + 2*xa))
         end if
         skewed = .true.
      end select
    end subroutine bc_di

    subroutine bc_dd
      select case ( imaxo )
      case (1)
         sm = ((3+2*xb)*b0 + (1-2*xb)*b1)/((1+2*xa)*(1+xa+xb))
         s0 = 4*((-1 + xa - xb)*b0 + (-1-xa+xb)*b1)/((1+2*xa)*(1+2*xb))
         sp = ((1-2*xa)*b0 + (3+2*xa)*b1)/((1+2*xb)*(1+xa+xb))
         ss = 0
         skewed = .false.
      case (2)
         sm = ((3+2*xb)*b0/((1+2*xa)*(1+xa+xb)))
         s0 = 4*(-1+xa-xb)*b0/((1+2*xa)*(1+2*xb)) - b1
         sp = b1
         ss = (1-2*xa)*b0/((1+xa*xb)*(1+2*xb))
         skewed = .true.
      case (3)
         if ( old_old ) then
            sm = 5*(5+2*xb)*b0/((3+4*xa*(2+xa))*(2+xa+xb))
            s0 = (-13-6*xb + 2*xa*(7+2*xb))*b0/((1+2*xa)*(3+2*xb)) - b1
            sp = - ((-1 + 2*xa)*(5+2*xb))*b0/((3+2*xa)*(1+2*xb))   + b1
            ss = 4*(-1 + 2*xa)*b0/((2+xa*xb)*(3+4*xb*(2+xb)))
         else 
            sm = (19 + 8*xb)*b0/((1 + 2*xa)*(3 + 2*xa)*(2 + xa + xb))
            s0 = ( &
                 + (-12 + 14*xa-6*xb+4*xa*xb)*b0 &
                 + (-3 - 6*xa - 2*xb - 4*xa*xb)*b1 &
                 ) /((1 + 2*xa)*(3 + 2*xb))
            sp = -( &
                 + (-4 + 10*xa - 2*xb + 4*xa*xb)*b0 &
                 + (-3 - 2*xa - 6*xb - 4*xa*xb)*b1 &
                 )/((1 + 2*xb)*(3 + 2*xa))
            ss = (-3 + 8*xa)*b0/((1 + 2*xb)*(3 + 2*xb)*(2 + xa + xb))
         end if
         skewed = .true.
      end select
    end subroutine bc_dd

    subroutine bc_dn
      select case ( imaxo )
      case (1)
         sm = 8*((1+xb)*b0 - xb*b1)/((1+2*xa)*(3+2*xa+4*xb))
         s0 = -8*((1+xb)*b0 + xb*b1)/((1+2*xa)*(4+2*xa+4*xb))
         sp = ((1-2*xa)*b0 + (3+2*xa)*b1)/(3+2*xa+4*xb)
         ss = 0
         skewed = .false.
      case (2)
         sm = 4*(3+2*xb)*b0/((1+2*xa)*(1+xa+xb))
         s0 = 4*(-1+xa-xb)*b0/((1+2*xa)*(1+2*xb)) - b1
         sp = b1
         ss = (1-2*xa)*b0/((1+xa+xb)*(1+2*xb))
         skewed = .true.
      case (3)
         sm = 4*(5+2*xb)*b0/((3+4*xa*(2+xa))*(2+xa+xb))
         s0 = (-13 + 6*xb +2*xa*(7+2*xb))*b0/((1+2*xa)*(3+2*xb)) - b1
         sp = -(-1+2*xa)*(5+2*xb)*b0/((3+2*xa)*(1+2*xb)) + b1
         ss = 4*(-1 + 2*xa)*b0/((2+xa+xb)*(3+4*xb*(2+xb)))
         skewed = .true.
      end select
    end subroutine bc_dn

    subroutine bc_ni
      select case ( imaxo )
      case (1)
         sm = -b0
         s0 = -b1
         sp =  b1
         ss = 0
         skewed = .false.
      case (2)
         sm = -b0/(1 + xa)
         s0 = xa*b0/(1 + xa) - b1
         sp = -xa*b0/(1 + xa) + b1
         ss = 0
         skewed = .false.
      case (3)
         sm = -24*b0/(23 + 12*xa*(3+xa))
         s0 = 2*((-1 + 12*xa*(2+xa))/(23 + 12*xa*(3+xa)))*b0 - b1
         sp = -3*((-1 + 4*xa*(5+3*xa))/(23 + 12*xa*(3+xa)))*b0 + b1
         ss = ((-1 + 12*xa*(1+xa))/(23 + 12*xa*(3+xa)))*b0
         skewed = .true.
      end select
    end subroutine bc_ni

    subroutine bc_nd
      select case ( imaxo )
      case (1)
         sm = - ((3+2*xb)*b0 + (1-2*xb)*b1)/(3+4*xa+2*xb)
         s0 = 8*(xa*b0 -(1+xa)*b1)/((1+2*xb)*(3+4*xa+2*xb))
         sp = 8*(-xa*b0 + (1+xb)*b1)/((1+2*xb)*(3+4*xa+2*xb))
         ss = 0
         skewed = .false.
      case (2)
         sm =  -(3+2*xb)*b0/(3+4*xa+2*xb)
         s0 = 8*xa*b0/((1+2*xb)*(3+4*xa+2*xb)) - b1
         sp = b1 
         ss = -8*xa*b0/((1+2*xb)*(3+4*xa+2*xb))
         skewed = .true.
      case (3)
         sm = (-4*(5 + 2*xb)*b0)/(19 + 8*xb + 4*xa*(8 + 3*xa + 2*xb))
         s0 = (-7 - 2*xb + 4*xa*(36 + 22*xb + 4*xb**2 + 3*xa*(7 + 2*xb)))*b0 - b1
         sp = -(((5 + 2*xb)* (-1 + 4*xa*(4 + 3*xa + 2*xb))* b0) &
              /((1 + 2*xb)*(19 + 8*xb + 4*xa*(8 + 3*xa + 2*xb)))) + b1
         ss = (8*(-1 + 12*xa*(1 + xa))*  b0) &
              /((3 + 4*xb*(2 + xb))*(19 + 8*xb +4*xa*(8 + 3*xa + 2*xb)))
         skewed = .true.
      end select
    end subroutine bc_nd

    subroutine bc_nn
      select case ( imaxo )
      case (1)
         sm = (-(1+xb)*b0 + xb*b1)/(1+xa+xb)
         s0 = 0
         sp = (-xa*b0 + (1+xa)*b1)/(1+xa+xb)
         ss = 0
         skewed = .false.
      case (2)
         sm = -(1+xb)*b0/(1+xa+xb)
         s0 = -b1
         sp =  b1
         ss = -xa*b0/(1+xa+xb)
         skewed = .true.
      case (3)
         sm = -(23+12*xb*(3+xb))*b0/((2+xa+xb)*(11+12*xb+12*xa*(1+xb)))
         s0 = (-1 + 12*xa*(2*xb))*b0/((11+12*xb+12*xa*(1+xb))) - b1
         sp = -(-1 + 12*xa*(2*xb))*b0/((11+12*xb+12*xa*(1+xb))) + b1
         ss = (-1 + 12*xa*(1+xa))*b0/((2+xa*xb)*(11+12*xb + 12*xa*(1+xb)))
         skewed = .true.
      end select
    end subroutine bc_nn

  end subroutine stencil_bndry_aaa
  
  !! -COMPUTES D^2(BETA*PHI)'s stencil
  subroutine s_laplacian_1d(ss, dh, alpha, beta)
    real (kind = dp_t), intent(out) :: ss(:,0:)
    real (kind = dp_t), intent(in) :: dh(:)
    real(kind=dp_t), intent(in) :: alpha, beta
    real (kind = dp_t) :: f1(size(dh))
    integer nx
    nx = size(ss,1)
    f1 = beta*ONE/dh**2
    ss(:,0) = -TWO*f1(1) + alpha
    ss(:,1) = +ONE*f1(1)
    ss(:,2) = +ONE*f1(1)
  end subroutine s_laplacian_1d

  subroutine s_laplacian_2d(ss, dh, alpha, beta)
    real (kind = dp_t), intent(out) :: ss(:,:,0:)
    real (kind = dp_t), intent(in) :: dh(:)
    real(kind=dp_t), intent(in) :: alpha, beta
    real (kind = dp_t) :: f1(size(dh))
    integer nx, ny
    nx = size(ss,1)
    ny = size(ss,2)
    f1 = beta*ONE/dh**2
    ss(:,:,0) = -TWO*sum(f1) + alpha
    ss(:,:,1) =  f1(1)
    ss(:,:,2) =  f1(1)
    ss(:,:,3) =  f1(2)
    ss(:,:,4) =  f1(2)
  end subroutine s_laplacian_2d

  subroutine s_laplacian_3d(ss, dh, alpha, beta)
    real (kind = dp_t), intent(out) :: ss(:,:,:,0:)
    real (kind = dp_t), intent(in) :: dh(:)
    real(kind=dp_t), intent(in) :: alpha, beta
    real (kind = dp_t) :: f1(size(dh))
    integer nx, ny, nz
    nx = size(ss,1)
    ny = size(ss,2)
    nz = size(ss,3)
    f1 = beta*ONE/dh**2
    ss(:,:,:,0) = -TWO*sum(f1) + alpha
    ss(:,:,:,1) =  f1(1)
    ss(:,:,:,2) =  f1(1)
    ss(:,:,:,3) =  f1(2)
    ss(:,:,:,4) =  f1(2)
    ss(:,:,:,5) =  f1(3)
    ss(:,:,:,6) =  f1(3)
  end subroutine s_laplacian_3d

  subroutine s_weighted_laplacian_1d(ss, bcoef, dh, alpha, beta)
    real (kind = dp_t), intent(out) :: ss(:,0:)
    real (kind = dp_t), intent(in)  :: bcoef(0:)
    real (kind = dp_t), intent(in) :: dh(:)
    real(kind=dp_t), intent(in) :: alpha, beta
    real (kind = dp_t) :: f1(size(dh))
    integer nx
    nx = size(ss,1)
    f1 = beta*ONE/dh**2
    ss(:,0) = -bcoef(1:nx  )*TWO*sum(f1) + alpha
    ss(:,1) =  bcoef(2:nx+1)
    ss(:,2) =  bcoef(0:nx-1)
  end subroutine s_weighted_laplacian_1d

  subroutine s_weighted_laplacian_2d(ss, bcoef, dh, alpha, beta)
    real (kind = dp_t), intent(out) :: ss(:,:,0:)
    real (kind = dp_t), intent(in)  :: bcoef(0:,0:)
    real (kind = dp_t), intent(in) :: dh(:)
    real(kind=dp_t), intent(in) :: alpha, beta
    real (kind = dp_t) :: f1(size(dh))
    integer nx, ny
    nx = size(ss,1)
    ny = size(ss,2)
    f1 = beta*ONE/dh**2
    ss(:,:,0) = -bcoef(1:nx  ,1:ny  )*TWO*sum(f1) + alpha
    ss(:,:,1) =  bcoef(2:nx+1,1:ny  )*f1(1)
    ss(:,:,2) =  bcoef(0:nx-1,1:ny  )*f1(1)
    ss(:,:,3) =  bcoef(1:nx  ,2:ny+1)*f1(2)
    ss(:,:,4) =  bcoef(1:nx  ,0:ny-1)*f1(2)
  end subroutine s_weighted_laplacian_2d

  subroutine s_weighted_laplacian_3d(ss, bcoef, dh, alpha, beta)
    real (kind = dp_t), intent(out) :: ss(:,:,:,0:)
    real (kind = dp_t), intent(in)  :: bcoef(0:,0:,0:)
    real (kind = dp_t), intent(in) :: dh(:)
    real(kind=dp_t), intent(in) :: alpha, beta
    real (kind = dp_t) :: f1(size(dh))
    integer nx, ny, nz
    nx = size(ss,1)
    ny = size(ss,2)
    nz = size(ss,3)
    f1 = beta*ONE/dh**2
    ss(:,:,:,0) = -bcoef(1:nx  ,1:ny  ,1:nz  )*TWO*sum(f1) + alpha
    ss(:,:,:,1) =  bcoef(2:nx+1,1:ny  ,1:nz  )*f1(1)
    ss(:,:,:,2) =  bcoef(0:nx-1,1:ny  ,1:nz  )*f1(1)
    ss(:,:,:,3) =  bcoef(1:nx  ,2:ny+1,1:nz  )*f1(2)
    ss(:,:,:,4) =  bcoef(1:nx  ,0:ny-1,1:nz  )*f1(2)
    ss(:,:,:,5) =  bcoef(1:nx  ,2:ny+1,2:nz+1)*f1(3)
    ss(:,:,:,6) =  bcoef(1:nx  ,0:ny-1,0:nz-1)*f1(3)
  end subroutine s_weighted_laplacian_3d

  subroutine s_mehrstellen_1d(ss, dh, alpha, beta)
    real (kind = dp_t), intent(out) :: ss(:,0:)
    real (kind = dp_t), intent(in) :: dh(:)
    real(kind=dp_t), intent(in) :: alpha, beta
    call bl_error("S_MEHRSTELLEN_1D: not ready")
    ss(:,1) = 1 ; ss(:,0) =  -2 ; ss(:,2) = 1
    ss = beta*ss/(dh(1)**2)
    ss(:,0)  = ss(:,0) + alpha
  end subroutine s_mehrstellen_1d

  subroutine s_mehrstellen_2d(ss, dh, alpha, beta)
    real (kind = dp_t), intent(out) :: ss(:,:,0:)
    real (kind = dp_t), intent(in) :: dh(:)
    real(kind=dp_t), intent(in) :: alpha, beta
    ss(:,:,1) =  1 ; ss(:,:,2) =   4 ; ss(:,:,3) =  1
    ss(:,:,4) =  4 ; ss(:,:,0) = -20 ; ss(:,:,5) =  4 
    ss(:,:,6) =  1 ; ss(:,:,7) =   4 ; ss(:,:,8) =  1
    ss = beta*ss/(6*dh(1)**2)
    ss(:,:,0) = ss(:,:,0) + alpha
  end subroutine s_mehrstellen_2d

  subroutine s_mehrstellen_3d(ss, dh, alpha, beta)
    real (kind = dp_t), intent(out) :: ss(:,:,:,0:)
    real (kind = dp_t), intent(in) :: dh(:)
    real(kind=dp_t), intent(in) :: alpha, beta
    ss(:,:,:, 1) =   0 ; ss(:,:,:, 2) =   1 ; ss(:,:,:, 3) =   0
    ss(:,:,:, 4) =   1 ; ss(:,:,:, 5) =   2 ; ss(:,:,:, 6) =   1
    ss(:,:,:, 7) =   0 ; ss(:,:,:, 8) =   1 ; ss(:,:,:, 9) =   0

    ss(:,:,:,10) =   1 ; ss(:,:,:,11) =   2 ; ss(:,:,:,12) =   1
    ss(:,:,:,13) =   2 ; ss(:,:,:, 0) = -24 ; ss(:,:,:,14) =   2 
    ss(:,:,:,15) =   1 ; ss(:,:,:,16) =   2 ; ss(:,:,:,17) =   1

    ss(:,:,:,18) =   0 ; ss(:,:,:,19) =   1 ; ss(:,:,:,20) =   0
    ss(:,:,:,21) =   1 ; ss(:,:,:,22) =   2 ; ss(:,:,:,23) =   1
    ss(:,:,:,24) =   0 ; ss(:,:,:,25) =   1 ; ss(:,:,:,26) =   0
    ss = beta*ss/(6*dh(1)**2)
    ss(:,:,:,0)  = ss(:,:,:,0) + alpha
  end subroutine s_mehrstellen_3d

  subroutine s_simple_1d_cc(ss, beta, dh, mask, xa, xb, order)
    integer, intent(inout) :: mask(:)
    real (kind = dp_t), intent(out) :: ss(:,0:)
    real (kind = dp_t), intent(in)  :: beta(0:,0:)
    real (kind = dp_t), intent(in) :: dh(:)
    real (kind = dp_t), intent(in) :: xa(:), xb(:)
    integer, intent(in) :: order
    integer nx
    integer i
    real (kind = dp_t) f1(1)

    integer, parameter :: XBC = 3
    integer :: bclo, bchi

    nx = size(ss,dim=1)

    f1 = ONE/dh**2

    mask = ibclr(mask, BC_BIT(BC_GEOM,1,-1))
    mask = ibclr(mask, BC_BIT(BC_GEOM,1,+1))

    do i = 1, nx
       ss(i,0) =   ZERO
       ss(i,1) = -beta(i+1,1)*f1(1)
       ss(i,2) = -beta(i  ,1)*f1(1)
       ss(i,XBC) = ZERO
    end do

    ! x derivatives

    do i = 2, nx-1
       ss(i,0) = ss(i,0) + (beta(i+1,1)+beta(i,1))*f1(1)
    end do

    bclo = stencil_bc_type(mask( 1),1,-1)
    bchi = stencil_bc_type(mask(nx),1,+1)

    call stencil_bndry_aaa(order, nx, 1, -1, mask(1), &
         ss(1,0), ss(1,1), ss(1,2), ss(1,XBC), &
         beta(1,1), beta(2,1), xa(1), xb(1), dh(1), bclo, bchi)
    if ( nx > 1 ) then
       call stencil_bndry_aaa(order, nx, 1, 1, mask(nx), &
            ss(nx,0), ss(nx,1), ss(nx,2), ss(nx,XBC), &
            beta(nx,1), beta(nx+1,1), xa(1), xb(1), dh(1), bclo, bchi)
    end if

    do i = 1, nx
       ss(i,0) = ss(i,0) + beta(i,0)
    end do

  end subroutine s_simple_1d_cc

  subroutine s_simple_2d_cc(ss, beta, dh, mask, xa, xb, order)
    integer, intent(inout) :: mask(:,:)
    real (kind = dp_t), intent(out) :: ss(:,:,0:)
    real (kind = dp_t), intent(in)  :: beta(0:,0:,0:)
    real (kind = dp_t), intent(in)  :: dh(:)
    real (kind = dp_t), intent(in)  :: xa(:), xb(:)
    integer, intent(in) :: order
    integer nx, ny
    integer i, j
    real (kind = dp_t) f1(2)

    integer :: bclo, bchi

    integer, parameter :: XBC = 5, YBC = 6

    nx = size(ss,dim=1)
    ny = size(ss,dim=2)

    f1 = ONE/dh**2

    mask = ibclr(mask, BC_BIT(BC_GEOM,1,-1))
    mask = ibclr(mask, BC_BIT(BC_GEOM,1,+1))
    mask = ibclr(mask, BC_BIT(BC_GEOM,2,-1))
    mask = ibclr(mask, BC_BIT(BC_GEOM,2,+1))

    do j = 1, ny
       do i = 1, nx
          ss(i,j,0) = ZERO
          ss(i,j,1) = -beta(i+1,j,1)*f1(1)
          ss(i,j,2) = -beta(i  ,j,1)*f1(1)
          ss(i,j,3) = -beta(i,j+1,2)*f1(2)
          ss(i,j,4) = -beta(i,j  ,2)*f1(2)
          ss(i,j,XBC) = ZERO
          ss(i,j,YBC) = ZERO
       end do
    end do

    ! x derivatives

    do j = 1, ny
       do i = 2, nx-1
          ss(i,j,0) = ss(i,j,0) + (beta(i,j,1)+beta(i+1,j,1))*f1(1)
       end do
    end do

    do j = 1,ny
       bclo = stencil_bc_type(mask( 1, j),1,-1)
       bchi = stencil_bc_type(mask(nx, j),1,+1)

       call stencil_bndry_aaa(order, nx, 1, -1, mask(1,j), &
            ss(1,j,0), ss(1,j,1), ss(1,j,2), ss(1,j,XBC), &
            beta(1,j,1), beta(2,j,1), xa(1), xb(1), dh(1), bclo, bchi)
       if ( nx > 1 ) then
          call stencil_bndry_aaa(order, nx, 1, 1, mask(nx,j), &
               ss(nx,j,0), ss(nx,j,1), ss(nx,j,2), ss(nx,j,XBC), &
               beta(nx,j,1), beta(nx+1,j,1), xa(1), xb(1), dh(1), bclo, bchi)
       end if
    end do

    ! y derivatives

    do i = 1, nx
       do j = 2, ny-1
          ss(i,j,0) = ss(i,j,0) + (beta(i,j,2)+beta(i,j+1,2))*f1(2)
       end do
    end do

    do i = 1, nx
       bclo = stencil_bc_type(mask( i, 1),2,-1)
       bchi = stencil_bc_type(mask( i,ny),2,+1)

       call stencil_bndry_aaa(order, ny, 2, -1, mask(i,1), &
            ss(i,1,0), ss(i,1,3), ss(i,1,4),ss(i,1,YBC), &
            beta(i,1,2), beta(i,2,2), xa(2), xb(2), dh(2), bclo, bchi)
       if ( ny > 1 ) then
          call stencil_bndry_aaa(order, ny, 2, 1, mask(i,ny), &
               ss(i,ny,0), ss(i,ny,3), ss(i,ny,4), ss(i,ny,YBC), &
               beta(i,ny,2), beta(i,ny+1,2), xa(2), xb(2), dh(2), bclo, bchi)
       end if
    end do

    if ( .false. ) then
       print *, 'xa = ', xa, ' xb = ', xb, ' dh = ', dh
       do i = 2, nx-1
          do j = 2, ny-1
             write(*,fmt='(i3," ", i3,":",7g20.10)') i,j, ss(i,j,0:)
          end do
       end do
       print *, 'x bcs'
       do j = 1, ny
          i = 1
          write(*,fmt='(i3," ", i3,":",7g20.10)') i,j, ss(i,j,0:)
          i = nx
          write(*,fmt='(i3," ", i3,":",7g20.10)') i,j, ss(i,j,0:)
       end do
       print *, 'y bcs'
       do i = 2, nx-1
          j = 1
          write(*,fmt='(i3," ", i3,":",7g20.10)') i,j, ss(i,j,0:)
          j = ny
          write(*,fmt='(i3," ", i3,":",7g20.10)') i,j, ss(i,j,0:)
       end do
       !   stop
    end if

    do j = 1, ny
       do i = 1, nx
          ss(i,j,0) = ss(i,j,0) + beta(i,j,0)
       end do
    end do

  end subroutine s_simple_2d_cc

  subroutine s_simple_3d_cc(ss, beta, dh, mask, xa, xb, order)
    integer nx,ny,nz
    integer, intent(inout) :: mask(:,:,:)
    real (kind = dp_t), intent(out) :: ss(:,:,:,0:)
    real (kind = dp_t), intent(in)  :: beta(0:,0:,0:,0:)
    real (kind = dp_t), intent(in) :: dh(:)
    real (kind = dp_t), intent(in) :: xa(:), xb(:)
    integer, intent(in) :: order
    integer i, j, k
    real (kind = dp_t) f1(3)

    integer :: bclo, bchi

    integer, parameter :: XBC = 7, YBC = 8, ZBC = 9

    nx = size(ss,dim=1)
    ny = size(ss,dim=2)
    nz = size(ss,dim=3)

    f1 = ONE/dh**2

    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             ss(i,j,k,0) = ZERO
             ss(i,j,k,1) = -beta(i+1,j,k,1)*f1(1)
             ss(i,j,k,2) = -beta(i  ,j,k,1)*f1(1)
             ss(i,j,k,3) = -beta(i,j+1,k,2)*f1(2)
             ss(i,j,k,4) = -beta(i,j  ,k,2)*f1(2)
             ss(i,j,k,5) = -beta(i,j,k+1,3)*f1(3)
             ss(i,j,k,6) = -beta(i,j,k  ,3)*f1(3)
             ss(i,j,k,XBC) = ZERO
             ss(i,j,k,YBC) = ZERO
             ss(i,j,k,ZBC) = ZERO
          end do
       end do
    end do

    mask = ibclr(mask, BC_BIT(BC_GEOM,1,-1))
    mask = ibclr(mask, BC_BIT(BC_GEOM,1,+1))
    mask = ibclr(mask, BC_BIT(BC_GEOM,2,-1))
    mask = ibclr(mask, BC_BIT(BC_GEOM,2,+1))
    mask = ibclr(mask, BC_BIT(BC_GEOM,3,-1))
    mask = ibclr(mask, BC_BIT(BC_GEOM,3,+1))

    ! x derivatives

    do k = 1, nz
       do j = 1, ny
          do i = 2, nx-1
             ss(i,j,k,0) = ss(i,j,k,0) + (beta(i,j,k,1)+beta(i+1,j,k,1))*f1(1)
          end do
       end do
    end do

    do k = 1,nz
       do j = 1,ny
          bclo = stencil_bc_type(mask( 1, j,k),1,-1)
          bchi = stencil_bc_type(mask(nx, j,k),1,+1)

          call stencil_bndry_aaa(order, nx, 1, -1, mask(1,j,k), &
               ss(1,j,k,0), ss(1,j,k,1), ss(1,j,k,2), ss(1,j,k,XBC), &
               beta(1,j,k,1), beta(2,j,k,1), xa(1), xb(1), dh(1), bclo, bchi)
          if ( nx > 1 ) then
             call stencil_bndry_aaa(order, nx, 1, 1, mask(nx,j,k), &
                  ss(nx,j,k,0), ss(nx,j,k,1), ss(nx,j,k,2), ss(nx,j,k,XBC), &
                  beta(nx,j,k,1), beta(nx+1,j,k,1), xa(1), xb(1), dh(1), bclo, bchi)
          end if
       end do
    end do

    ! y derivatives

    do k = 1,nz
       do i = 1,nx
          do j = 2,ny-1
             ss(i,j,k,0) = ss(i,j,k,0) + (beta(i,j,k,2)+beta(i,j+1,k,2))*f1(2)
          end do
       end do
    end do

    do k = 1, nz
       do i = 1, nx
          bclo = stencil_bc_type(mask(i,1,k) ,2,-1)
          bchi = stencil_bc_type(mask(i,ny,k),2,+1)

          call stencil_bndry_aaa(order, ny, 2, -1, mask(i,1,k), &
               ss(i,1,k,0), ss(i,1,k,3), ss(i,1,k,4),ss(i,1,k,YBC), &
               beta(i,1,k,2), beta(i,2,k,2), xa(2), xb(2), dh(2), bclo, bchi)
          if ( ny > 1 ) then
             call stencil_bndry_aaa(order, ny, 2, 1, mask(i,ny,k), &
                  ss(i,ny,k,0), ss(i,ny,k,3), ss(i,ny,k,4), ss(i,ny,k,YBC), &
                  beta(i,ny,k,2), beta(i,ny+1,k,2), xa(2), xb(2), dh(2), bclo, bchi)
          end if
       end do
    end do

    ! z derivatives

    do j = 1, ny
       do i = 1, nx
          do k = 2, nz-1
             ss(i,j,k,0) = ss(i,j,k,0) + (beta(i,j,k,3)+beta(i,j,k+1,3))*f1(3)
          end do
       end do
    end do

    do j = 1, ny
       do i = 1, nx
          bclo = stencil_bc_type(mask(i,j,1) ,3,-1)
          bchi = stencil_bc_type(mask(i,j,nz),3,+1)

          call stencil_bndry_aaa(order, nz, 3, -1, mask(i,j,1), &
               ss(i,j,1,0), ss(i,j,1,5), ss(i,j,1,6),ss(i,j,1,ZBC), &
               beta(i,j,1,3), beta(i,j,2,3), xa(3), xb(3), dh(3), bclo, bchi)
          if ( nz > 1 ) then
             call stencil_bndry_aaa(order, nz, 3, 1, mask(i,j,nz), &
                  ss(i,j,nz,0), ss(i,j,nz,5), ss(i,j,nz,6), ss(i,j,nz,ZBC), &
                  beta(i,j,nz,3), beta(i,j,nz+1,3), xa(3), xb(3), dh(3), bclo, bchi)
          end if
       end do
    end do

    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             ss(i,j,k,0) = ss(i,j,k,0) + beta(i,j,k,0)
          end do
       end do
    end do

  end subroutine s_simple_3d_cc

  subroutine stencil_apply_1d(ss, dd, uu, mm, ng, skwd)
    integer, intent(in) :: ng
    real (kind = dp_t), intent(in)  :: ss(:,0:)
    real (kind = dp_t), intent(out) :: dd(:)
    real (kind = dp_t), intent(in)  :: uu(1-ng:)
    integer           , intent(in)  :: mm(:)
    logical, intent(in), optional :: skwd
    integer nx
    integer i
    integer, parameter :: XBC = 3
    logical :: lskwd
   
    lskwd = .true.; if ( present(skwd) ) lskwd = skwd

    nx = size(ss,dim=1)

    do i = 1,nx
       dd(i) = ss(i,0)*uu(i) + ss(i,1)*uu(i+1) + ss(i,2)*uu(i-1)
    end do
    if ( lskwd ) then
       if (nx >= 2) then
          i = 1
          if (bc_skewed(mm(i),1,+1)) then
             dd(i) = dd(i) + ss(i,XBC)*uu(i+2)
          end if
  
          i = nx
          if (bc_skewed(mm(i),1,-1)) then
             dd(i) = dd(i) + ss(i,XBC)*uu(i-2)
          end if
       end if
    end if

  end subroutine stencil_apply_1d

  subroutine stencil_flux_1d(ss, flux, uu, mm, ng, ratio, face, dim, skwd)
    integer, intent(in) :: ng
    real (kind = dp_t), intent(in)  :: ss(:,0:)
    real (kind = dp_t), intent(out) :: flux(:)
    real (kind = dp_t), intent(in)  :: uu(1-ng:)
    integer           , intent(in)  :: mm(:)
    logical, intent(in), optional :: skwd
    integer, intent(in) :: ratio, face, dim
    integer nx
    integer i
    integer, parameter :: XBC = 3

    real (kind = dp_t) :: fac

    logical :: lskwd

    lskwd = .true. ; if ( present(skwd) ) lskwd = skwd

    nx = size(ss,dim=1)

    !   This factor is dx^fine / dx^crse
    fac = ONE / real(ratio, kind=dp_t)

    if ( dim == 1 ) then
       if ( face == -1 ) then
          i = 1
          if (bc_dirichlet(mm(1),1,-1)) then
             flux(1) = ss(i,1)*(uu(i+1)-uu(i)) + ss(i,2)*(uu(i-1)-uu(i)) &
                  - ss(i+1,2)*(uu(i+1)-uu(i))
             if (bc_skewed(mm(i),1,+1)) then
                flux(1) =  flux(1) + ss(i,XBC)*uu(i+2)
             end if
          else 
             flux(1) = Huge(flux)
          end if
          flux(1) = fac*flux(1)
       else if ( face == 1 ) then
          i = nx
          if (bc_dirichlet(mm(i),1,+1)) then
             flux(1) = ss(i,1)*(uu(i+1)-uu(i)) + ss(i,2)*(uu(i-1)-uu(i)) &
                  - ss(i-1,1)*(uu(i-1)-uu(i))
             if (bc_skewed(mm(i),1,-1)) then
                flux(1) =  flux(1) + ss(i,XBC)*uu(i-2)
             end if
          else 
             flux(1) = Huge(flux)
          end if
          flux(1) = fac*flux(1)
       end if
    end if

  end subroutine stencil_flux_1d

  subroutine stencil_apply_2d(ss, dd, uu, mm, ng, skwd)
    integer, intent(in) :: ng
    real (kind = dp_t), intent(in ) :: uu(1-ng:,1-ng:)
    real (kind = dp_t), intent(out) :: dd(:,:)
    real (kind = dp_t), intent(in ) :: ss(:,:,0:)
    logical, intent(in), optional :: skwd
    integer           , intent(in)  :: mm(:,:)
    integer nx,ny
    integer i,j

    integer, parameter :: XBC = 5, YBC = 6

    logical :: lskwd

    lskwd = .true.; if ( present(skwd) ) lskwd = skwd

    nx = size(ss,dim=1)
    ny = size(ss,dim=2)

    do j = 1,ny
       do i = 1,nx
          dd(i,j) = ss(i,j,0)*uu(i,j) &
               + ss(i,j,1)*uu(i+1,j  ) + ss(i,j,2)*uu(i-1,j  ) &
               + ss(i,j,3)*uu(i  ,j+1) + ss(i,j,4)*uu(i  ,j-1)
       end do
    end do

    if ( lskwd ) then
       ! Corrections for skewed stencils
       if (nx > 1) then
          do j = 1, ny

             i = 1
             if (bc_skewed(mm(i,j),1,+1)) then
                dd(i,j) = dd(i,j) + ss(i,j,XBC)*uu(i+2,j)
             end if

             i = nx
             if (bc_skewed(mm(i,j),1,-1)) then
                dd(i,j) = dd(i,j) + ss(i,j,XBC)*uu(i-2,j)
             end if
          end do
       end if

       if (ny > 1) then
          do i = 1, nx

             j = 1
             if (bc_skewed(mm(i,j),2,+1)) then
                dd(i,j) = dd(i,j) + ss(i,j,YBC)*uu(i,j+2)
             end if

             j = ny
             if (bc_skewed(mm(i,j),2,-1)) then
                dd(i,j) = dd(i,j) + ss(i,j,YBC)*uu(i,j-2)
             end if

          end do
       end if
    end if
  end subroutine stencil_apply_2d

  subroutine stencil_flux_2d(ss, flux, uu, mm, ng, ratio, face, dim, skwd)
    integer, intent(in) :: ng
    real (kind = dp_t), intent(in ) :: uu(1-ng:,1-ng:)
    real (kind = dp_t), intent(out) :: flux(:,:)
    real (kind = dp_t), intent(in ) :: ss(:,:,0:)
    integer           , intent(in)  :: mm(:,:)
    logical, intent(in), optional :: skwd
    integer, intent(in) :: ratio, face, dim
    integer nx,ny
    integer i,j,ic,jc
    real (kind = dp_t) :: fac
    integer, parameter :: XBC = 5, YBC = 6
    logical :: lskwd

    lskwd = .true. ; if ( present(skwd) ) lskwd = skwd

    nx = size(ss,dim=1)
    ny = size(ss,dim=2)

    !   Note that one factor of ratio is the tangential averaging, while the
    !     other is the normal factor
    fac = ONE/real(ratio*ratio, kind=dp_t)

!   Lo i face
    if ( dim == 1 ) then
       if (face == -1) then

          i = 1
          flux(1,:) = ZERO
          do j = 1,ny
             jc = (j-1)/ratio+1
             if (bc_dirichlet(mm(i,j),1,-1)) then
                flux(1,jc) = flux(1,jc)  &
                     + ss(i,j,1)*(uu(i+1,j)-uu(i,j)) &
                     + ss(i,j,2)*(uu(i-1,j)-uu(i,j)) - ss(i+1,j,2)*(uu(i+1,j)-uu(i,j))
                if (bc_skewed(mm(i,j),1,+1)) &
                     flux(1,jc) = flux(1,jc) + ss(i,j,XBC)*(uu(i+2,j)-uu(i,j)) 
             else   
                flux(1,jc) = Huge(flux)
             end if
          end do
          flux(1,:) = fac * flux(1,:)

!      Hi i face
       else if (face == 1) then

          i = nx
          flux(1,:) = ZERO
          do j = 1,ny
             jc = (j-1)/ratio+1
             if (bc_dirichlet(mm(i,j),1,+1)) then

                flux(1,jc) = flux(1,jc) &
                     + ss(i,j,1)*(uu(i+1,j)-uu(i,j)) &
                     + ss(i,j,2)*(uu(i-1,j)-uu(i,j)) - ss(i-1,j,1)*(uu(i-1,j)-uu(i,j))
                if (bc_skewed(mm(i,j),1,-1)) &
                     flux(1,jc) = flux(1,jc) + ss(i,j,XBC)*(uu(i-2,j)-uu(i,j))
             else 
                flux(1,jc) = Huge(flux)
             end if
          end do
          flux(1,:) = fac * flux(1,:)

       end if

!   Lo j face
    else if ( dim == 2 ) then
       if (face == -1) then

          j = 1
          flux(:,1) = ZERO
          do i = 1,nx
             ic = (i-1)/ratio+1
             if (bc_dirichlet(mm(i,j),2,-1)) then
                flux(ic,1) = flux(ic,1)  &
                     + ss(i,j,3)*(uu(i,j+1)-uu(i,j)) &
                     + ss(i,j,4)*(uu(i,j-1)-uu(i,j)) - ss(i,j+1,4)*(uu(i,j+1)-uu(i,j))
                if (bc_skewed(mm(i,j),2,+1)) &
                     flux(ic,1) =  flux(ic,1) + ss(i,j,YBC)*(uu(i,j+2)-uu(i,j))
             else 
                flux(ic,1) = Huge(flux)
             end if
          end do
          flux(:,1) = fac * flux(:,1)


!      Hi j face
       else if (face == 1) then

          j = ny
          flux(:,1) = ZERO
          do i = 1,nx
             ic = (i-1)/ratio+1
             if (bc_dirichlet(mm(i,j),2,+1)) then
                flux(ic,1) = flux(ic,1)  &
                     + ss(i,j,3)*(uu(i,j+1)-uu(i,j)) &
                     + ss(i,j,4)*(uu(i,j-1)-uu(i,j)) - ss(i,j-1,3)*(uu(i,j-1)-uu(i,j))
                if (bc_skewed(mm(i,j),2,-1)) &
                     flux(ic,1) = flux(ic,1) + ss(i,j,YBC)*(uu(i,j-2)-uu(i,j))
             else
                flux(ic,1) = Huge(flux)
             end if
          end do
          flux(:,1) = fac * flux(:,1)

       end if
    end if

  end subroutine stencil_flux_2d

  subroutine stencil_apply_3d(ss, dd, uu, mm, ng, skwd)
    integer, intent(in) :: ng
    real (kind = dp_t), intent(in) :: uu(1-ng:,1-ng:,1-ng:)
    real (kind = dp_t), intent(out) :: dd(:,:,:)
    real (kind = dp_t), intent(in) :: ss(:,:,:,0:)
    logical, intent(in), optional :: skwd
    integer           , intent(in) :: mm(:,:,:)
    integer nx,ny,nz
    integer i,j,k

    integer, parameter :: XBC = 7, YBC = 8, ZBC = 9
    logical :: lskwd

    lskwd = .true.; if ( present(skwd) ) lskwd = skwd

    nx = size(ss,dim=1)
    ny = size(ss,dim=2)
    nz = size(ss,dim=3)

    do k = 1,nz
       do j = 1,ny
          do i = 1,nx
             dd(i,j,k) = ss(i,j,k,0)*uu(i,j,k) &
                  + ss(i,j,k,1)*uu(i+1,j  ,k  ) + ss(i,j,k,2)*uu(i-1,j  ,k  ) &
                  + ss(i,j,k,3)*uu(i  ,j+1,k  ) + ss(i,j,k,4)*uu(i  ,j-1,k  ) &
                  + ss(i,j,k,5)*uu(i  ,j  ,k+1) + ss(i,j,k,6)*uu(i  ,j  ,k-1)
          end do
       end do
    end do

    if ( lskwd ) then
       ! Corrections for skewed stencils
       if (nx > 1) then
          do k = 1, nz
             do j = 1, ny
                i = 1
                if (bc_skewed(mm(i,j,k),1,+1)) then
                   dd(i,j,k) = dd(i,j,k) + ss(i,j,k,XBC)*uu(i+2,j,k)
                end if

                i = nx
                if (bc_skewed(mm(i,j,k),1,-1)) then
                   dd(i,j,k) = dd(i,j,k) + ss(i,j,k,XBC)*uu(i-2,j,k)
                end if
             end do
          end do
       end if

       if (ny > 1) then
          do k = 1,nz
             do i = 1,nx
                j = 1
                if (bc_skewed(mm(i,j,k),2,+1)) then
                   dd(i,j,k) = dd(i,j,k) + ss(i,j,k,YBC)*uu(i,j+2,k)
                end if

                j = ny
                if (bc_skewed(mm(i,j,k),2,-1)) then
                   dd(i,j,k) = dd(i,j,k) + ss(i,j,k,YBC)*uu(i,j-2,k)
                end if
             end do
          end do
       end if

       if (nz > 1) then
          do j = 1,ny
             do i = 1,nx
                k = 1
                if (bc_skewed(mm(i,j,k),3,+1)) then
                   dd(i,j,k) = dd(i,j,k) + ss(i,j,k,ZBC)*uu(i,j,k+2)
                end if

                k = nz
                if (bc_skewed(mm(i,j,k),3,-1)) then
                   dd(i,j,k) = dd(i,j,k) + ss(i,j,k,ZBC)*uu(i,j,k-2)
                end if
             end do
          end do
       end if
    end if
  end subroutine stencil_apply_3d

  subroutine stencil_flux_3d(ss, flux, uu, mm, ng, ratio, face, dim, skwd)
    integer, intent(in) :: ng
    real (kind = dp_t), intent(in ) :: uu(1-ng:,1-ng:,1-ng:)
    real (kind = dp_t), intent(out) :: flux(:,:,:)
    real (kind = dp_t), intent(in ) :: ss(:,:,:,0:)
    integer           , intent(in)  :: mm(:,:,:)
    logical, intent(in), optional :: skwd
    integer, intent(in) :: ratio, face, dim
    integer nx, ny, nz
    integer i,j,k,ic,jc,kc
    integer, parameter :: XBC = 7, YBC = 8, ZBC = 9
    real (kind = dp_t) :: fac
    logical :: lskwd

    lskwd = .true. ; if ( present(skwd) ) lskwd = skwd

    nx = size(ss,dim=1)
    ny = size(ss,dim=2)
    nz = size(ss,dim=3)

    ! Note that two factors of ratio is from the tangential averaging, while the
    ! other is the normal factor
    fac = ONE/real(ratio*ratio*ratio, kind=dp_t)

    !   Lo i face
    if ( dim ==  1 ) then
       if (face == -1) then

          i = 1
          flux(1,:,:) = ZERO
          do k = 1,nz
             do j = 1,ny
                jc = (j-1)/ratio + 1
                kc = (k-1)/ratio + 1
                if (bc_dirichlet(mm(i,j,k),1,-1)) then
                   flux(1,jc,kc) =  flux(1,jc,kc) &
                        + ss(i,j,k,1)*(uu(i+1,j,k)-uu(i,j,k)) &
                        + ss(i,j,k,2)*(uu(i-1,j,k)-uu(i,j,k)) &
                        - ss(i+1,j,k,2)*(uu(i+1,j,k)-uu(i,j,k))
                   if (bc_skewed(mm(i,j,k),1,+1)) &
                        flux(1,jc,kc) =  flux(1,jc,kc) + ss(i,j,k,XBC)*(uu(i+2,j,k)-uu(i,j,k))
                else 
                   flux(1,jc,kc) = Huge(flux)
                end if
             end do
          end do
          flux(1,:,:) = flux(1,:,:) * fac

          !   Hi i face
       else if (face ==  1) then

          i = nx
          flux(1,:,:) = ZERO
          do k = 1,nz
             do j = 1,ny
                jc = (j-1)/ratio + 1
                kc = (k-1)/ratio + 1
                if (bc_dirichlet(mm(i,j,k),1,+1)) then
                   flux(1,jc,kc) =  flux(1,jc,kc) &
                        + ss(i,j,k,1)*(uu(i+1,j,k)-uu(i,j,k)) &
                        + ss(i,j,k,2)*(uu(i-1,j,k)-uu(i,j,k)) &
                        - ss(i-1,j,k,1)*(uu(i-1,j,k)-uu(i,j,k))
                   if (bc_skewed(mm(i,j,k),1,-1)) &
                        flux(1,jc,kc) =  flux(1,jc,kc) + ss(i,j,k,XBC)*(uu(i-2,j,k)-uu(i,j,k))
                else 
                   flux(1,jc,kc) = Huge(flux)
                end if

             end do
          end do
          flux(1,:,:) = flux(1,:,:) * fac
       end if
       !   Lo j face
    else if ( dim == 2 ) then
       if (face == -1) then
          j = 1
          flux(:,1,:) = ZERO
          do k = 1,nz
             do i = 1,nx
                ic = (i-1)/ratio + 1
                kc = (k-1)/ratio + 1
                if (bc_dirichlet(mm(i,j,k),2,-1)) then
                   flux(ic,1,kc) =  flux(ic,1,kc) &
                        + ss(i,j,k,3)*(uu(i,j+1,k)-uu(i,j,k)) &
                        + ss(i,j,k,4)*(uu(i,j-1,k)-uu(i,j,k)) &
                        - ss(i,j+1,k,4)*(uu(i,j+1,k)-uu(i,j,k))
                   if (bc_skewed(mm(i,j,k),2,+1)) &
                        flux(ic,1,kc) =  flux(ic,1,kc) + ss(i,j,k,YBC)*(uu(i,j+2,k)-uu(i,j,k))
                else 
                   flux(ic,1,kc) = Huge(flux)
                end if
             end do
          end do
          flux(:,1,:) = flux(:,1,:) * fac
          !   Hi j face
       else if (face ==  1) then
          j = ny
          flux(:,1,:) = ZERO
          do k = 1,nz
             do i = 1,nx
                ic = (i-1)/ratio + 1
                kc = (k-1)/ratio + 1

                if (bc_dirichlet(mm(i,j,k),2,+1)) then
                   flux(ic,1,kc) =  flux(ic,1,kc) &
                        + ss(i,j,k,3)*(uu(i,j+1,k)-uu(i,j,k)) &
                        + ss(i,j,k,4)*(uu(i,j-1,k)-uu(i,j,k)) &
                        - ss(i,j-1,k,3)*(uu(i,j-1,k)-uu(i,j,k))
                   if (bc_skewed(mm(i,j,k),2,-1)) &
                        flux(ic,1,kc) =  flux(ic,1,kc) + ss(i,j,k,YBC)*(uu(i,j-2,k)-uu(i,j,k))
                else
                   flux(ic,1,kc) = Huge(flux)
                end if
             end do
          end do
          flux(:,1,:) = flux(:,1,:) * fac
          !   Lo k face
       end if
    else if ( dim == 3 ) then
       if (face == -1) then

          k = 1
          flux(:,:,1) = ZERO
          do j = 1,ny
             do i = 1,nx
                ic = (i-1)/ratio + 1
                jc = (j-1)/ratio + 1
                if (bc_dirichlet(mm(i,j,k),3,-1)) then
                   flux(ic,jc,1) =  flux(ic,jc,1) &
                        + ss(i,j,k,5)*(uu(i,j,k+1)-uu(i,j,k)) &
                        + ss(i,j,k,6)*(uu(i,j,k-1)-uu(i,j,k)) &
                        - ss(i,j,k+1,6)*(uu(i,j,k+1)-uu(i,j,k))
                   if (bc_skewed(mm(i,j,k),3,+1)) &
                        flux(ic,jc,1) =  flux(ic,jc,1) + ss(i,j,k,ZBC)*(uu(i,j,k+2)-uu(i,j,k)) 
                else 
                   flux(ic,jc,1) = Huge(flux)
                end if
             end do
          end do
          flux(:,:,1) = flux(:,:,1) * fac

          !   Hi k face
       else if (face ==  1) then

          k = nz
          flux(:,:,1) = ZERO
          do j = 1,ny
             do i = 1,nx
                ic = (i-1)/ratio + 1
                jc = (j-1)/ratio + 1
                if (bc_dirichlet(mm(i,j,k),3,+1)) then
                   flux(ic,jc,1) =  flux(ic,jc,1) &
                        + ss(i,j,k,5)*(uu(i,j,k+1)-uu(i,j,k)) &
                        + ss(i,j,k,6)*(uu(i,j,k-1)-uu(i,j,k)) &
                        - ss(i,j,k-1,5)*(uu(i,j,k-1)-uu(i,j,k))
                   if (bc_skewed(mm(i,j,k),3,-1)) &
                        flux(ic,jc,1) =  flux(ic,jc,1) + ss(i,j,k,ZBC)*(uu(i,j,k-2)-uu(i,j,k))
                else
                   flux(ic,jc,1) = Huge(flux)
                end if
             end do
          end do
          flux(:,:,1) = flux(:,:,1) * fac

       end if
    end if

  end subroutine stencil_flux_3d

  subroutine stencil_dense_apply_1d(ss, dd, uu, mm, ng)
    integer, intent(in) :: ng
    real (kind = dp_t), intent(in   ) :: ss(:,0:)
    real (kind = dp_t), intent(  out) :: dd(:)
    real (kind = dp_t), intent(in   ) :: uu(1-ng:)
    integer           , intent(in   ) :: mm(:)
    integer i, nx
   
    nx = size(ss,dim=1)
    do i = 1, nx
      dd(i) = ss(i,1)*uu(i-1) + ss(i,0)*uu(i) + ss(i,2)*uu(i+1)
    end do

  end subroutine stencil_dense_apply_1d

  subroutine stencil_dense_apply_2d(ss, dd, uu, mm, ng)
    integer, intent(in) :: ng
    real (kind = dp_t), intent(in   ) :: uu(1-ng:,1-ng:)
    real (kind = dp_t), intent(  out) :: dd(:,:)
    real (kind = dp_t), intent(in   ) :: ss(:,:,0:)
    integer           , intent(in   ) :: mm(:,:)
    integer i, j, nx, ny

    nx = size(ss,dim=1)
    ny = size(ss,dim=2)

    do j = 1, ny
       do i = 1, nx
          dd(i,j) = &
               + ss(i,j,1)*uu(i-1,j-1) + ss(i,j,2)*uu(i  ,j-1) + ss(i,j,3)*uu(i+1,j-1) &
               + ss(i,j,4)*uu(i-1,j  ) + ss(i,j,0)*uu(i  ,j  ) + ss(i,j,5)*uu(i+1,j  ) &
               + ss(i,j,6)*uu(i-1,j+1) + ss(i,j,7)*uu(i  ,j+1) + ss(i,j,8)*uu(i+1,j+1)
       end do
    end do

  end subroutine stencil_dense_apply_2d

  subroutine stencil_dense_apply_3d(ss, dd, uu, mm, ng)
    integer, intent(in) :: ng
    real (kind = dp_t), intent(in   ) :: uu(1-ng:,1-ng:,1-ng:)
    real (kind = dp_t), intent(  out) :: dd(:,:,:)
    real (kind = dp_t), intent(in   ) :: ss(:,:,:,0:)
    integer           , intent(in   ) :: mm(:,:,:)
    integer i, j, k, nx, ny, nz

    nx = size(ss,dim=1)
    ny = size(ss,dim=2)
    nz = size(ss,dim=3)

    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             dd(i,j,k) = &
                  + ss(i,j,k, 1)*uu(i-1,j-1,k-1) &
                  + ss(i,j,k, 2)*uu(i  ,j-1,k-1) &
                  + ss(i,j,k, 3)*uu(i+1,j-1,k-1) &
                  + ss(i,j,k, 4)*uu(i-1,j  ,k-1) &
                  + ss(i,j,k, 5)*uu(i  ,j  ,k-1) &
                  + ss(i,j,k, 6)*uu(i+1,j  ,k-1) &
                  + ss(i,j,k, 7)*uu(i-1,j+1,k-1) &
                  + ss(i,j,k, 8)*uu(i  ,j+1,k-1) &
                  + ss(i,j,k, 9)*uu(i+1,j+1,k-1) &

                  + ss(i,j,k,10)*uu(i-1,j-1,k  ) &
                  + ss(i,j,k,11)*uu(i  ,j-1,k  ) &
                  + ss(i,j,k,12)*uu(i+1,j-1,k  ) &
                  + ss(i,j,k,13)*uu(i-1,j  ,k  ) &
                  + ss(i,j,k, 0)*uu(i  ,j  ,k  ) &
                  + ss(i,j,k,14)*uu(i+1,j  ,k  ) &
                  + ss(i,j,k,15)*uu(i-1,j+1,k  ) &
                  + ss(i,j,k,16)*uu(i  ,j+1,k  ) &
                  + ss(i,j,k,17)*uu(i+1,j+1,k  ) &

                  + ss(i,j,k,18)*uu(i-1,j-1,k+1) &
                  + ss(i,j,k,19)*uu(i  ,j-1,k+1) &
                  + ss(i,j,k,20)*uu(i+1,j-1,k+1) &
                  + ss(i,j,k,21)*uu(i-1,j  ,k+1) &
                  + ss(i,j,k,22)*uu(i  ,j  ,k+1) &
                  + ss(i,j,k,23)*uu(i+1,j  ,k+1) &
                  + ss(i,j,k,24)*uu(i-1,j+1,k+1) &
                  + ss(i,j,k,25)*uu(i  ,j+1,k+1) &
                  + ss(i,j,k,26)*uu(i+1,j+1,k+1)
          end do
       end do
    end do

  end subroutine stencil_dense_apply_3d
  
  ! polyInterpCoeff:
  !  
  ! This routine returns the Lagrange interpolating coefficients for a
  ! polynomial through N points, evaluated at xInt (see Numerical Recipes,
  ! v2, p102, e.g.):
  !
  !          (x-x2)(x-x3)...(x-xN)              (x-x1)(x-x2)...(x-x(N-1))
  ! P(x) = ----------------------- y1  + ... + ------------------------  yN
  !         (x1-x2)(x1-x3)...(x1-xN)            (x1-x2)(x1-x3)...(x1-xN)
  !
  ! P(xInt) = sum_(i=1)^(N) y[i]*c[i]
  !

  subroutine poly_interp_coeff(c, xInt, x)
    real(kind=dp_t), intent(in) :: xInt, x(:)
    real(kind=dp_t), intent(out) :: c(:)
    real(kind=dp_t) num, den
    integer i, j, N
    N = size(x)
!print *, 'N = ', N
    do j = 1, N
       num = ONE
       den = ONE
       do i = 1, j - 1
          num = num*(xInt - x(i))
          den = den*(x(j) - x(i))
       end do
       do i = j + 1, N
          num = num*(xInt - x(i))
          den = den*(x(j) - x(i))
       end do
       if (den == ZERO) then
          print *, 'xInt = ', x
          print *, 'j    = ', j
          print *, 'x    = ', x
          print *, 'c    = ', c
          call bl_error('polyInterpCoeff::invalid data')
       end if
       c(j) = num/den
    end do
  end subroutine poly_interp_coeff

  !     
  !     This is a test driver for the routine polyInterpCoeff.  Sample data
  !     is created from the statement function, and the location of the 
  !     boundary node and internal nodes are set, as apporpriate.  The
  !     number of points created is equal to the test NORDER set at the
  !     top of this file through a define.  The coefficients are computed,
  !     and then the ghost cell value is constructed from the resulting
  !     coefficients and written out.
  !

  subroutine t_polyInterpCoeffTest(norder)
    integer, intent(in) :: NORDER
    integer j
    real(kind=dp_t) c(0:NORDER-1), ci(0:NORDER-1)
    real(kind=dp_t) y(0:NORDER-1)
    real(kind=dp_t) x(0:NORDER-1)
    real(kind=dp_t) xInt

    call random_number(ci)
    
    x = (/ ZERO, (j+HALF,j=0,NORDER-2) /)
    do j = 0, NORDER-2
       y(j) = horner(x(j), ci)
    end do

    xInt = -HALF

    call poly_interp_coeff(c, xInt, x)

    print *, 'x = ', x
    print *, 'y = ', y
    print *, 'c = ', c
    print *, 'Interpolated y = ', sum(c*y)

  contains

    function Horner(xx, cc) result(r)
      real(kind=dp_t) :: r
      real(kind=dp_t), intent(in) :: xx
      real(kind=dp_t), intent(in) :: cc(:)
      integer :: i

      r = cc(1)
      do i = 2, size(cc)
         r = xx*r + cc(i)
      end do

    end function Horner

  end subroutine t_polyInterpCoeffTest

  subroutine stencil_fine_flux_1d(ss, flux, uu, mm, ng, face, dim, skwd)
    integer, intent(in) :: ng
    real (kind = dp_t), intent(in)  :: ss(:,0:)
    real (kind = dp_t), intent(out) :: flux(:)
    real (kind = dp_t), intent(in)  :: uu(1-ng:)
    integer           , intent(in)  :: mm(:)
    logical, intent(in), optional :: skwd
    integer, intent(in) :: face, dim
    integer nx
    integer i
    integer, parameter :: XBC = 3

    real (kind = dp_t) :: fac

    logical :: lskwd

    lskwd = .true. ; if ( present(skwd) ) lskwd = skwd

    nx = size(ss,dim=1)

    if ( dim == 1 ) then
       if ( face == -1 ) then
!         Lo i face
          i = 1
          if (bc_dirichlet(mm(1),1,-1)) then
             flux(1) = ss(i,1)*(uu(i+1)-uu(i)) + ss(i,2)*(uu(i-1)-uu(i)) &
                  - ss(i+1,2)*(uu(i+1)-uu(i))
             if (bc_skewed(mm(i),1,+1)) then
                flux(1) =  flux(1) + ss(i,XBC)*uu(i+2)
             end if
          else 
             flux(1) = ss(i,2)*(uu(i-1)-uu(i))
          end if
       else if ( face == 1 ) then

!         Hi i face
          i = nx
          if (bc_dirichlet(mm(i),1,+1)) then
             flux(1) = ss(i,1)*(uu(i+1)-uu(i)) + ss(i,2)*(uu(i-1)-uu(i)) &
                  - ss(i-1,1)*(uu(i-1)-uu(i))
             if (bc_skewed(mm(i),1,-1)) then
                flux(1) =  flux(1) + ss(i,XBC)*uu(i-2)
             end if
          else 
             flux(1) = ss(i,1)*(uu(i+1)-uu(i))
          end if
       end if
    end if

  end subroutine stencil_fine_flux_1d

  subroutine stencil_fine_flux_2d(ss, flux, uu, mm, ng, face, dim, skwd)
    integer, intent(in) :: ng
    real (kind = dp_t), intent(in ) :: uu(1-ng:,1-ng:)
    real (kind = dp_t), intent(out) :: flux(:,:)
    real (kind = dp_t), intent(in ) :: ss(:,:,0:)
    integer           , intent(in)  :: mm(:,:)
    logical, intent(in), optional :: skwd
    integer, intent(in) :: face, dim
    integer nx,ny
    integer i,j
    real (kind = dp_t) :: fac
    integer, parameter :: XBC = 5, YBC = 6
    logical :: lskwd

    lskwd = .true. ; if ( present(skwd) ) lskwd = skwd

    nx = size(ss,dim=1)
    ny = size(ss,dim=2)

!   Lo i face
    if ( dim == 1 ) then
       if (face == -1) then

          i = 1
          flux(1,:) = ZERO
          do j = 1,ny
             if (bc_dirichlet(mm(i,j),1,-1)) then
                flux(1,j) = &
                       ss(i,j,1)*(uu(i+1,j)-uu(i,j)) &
                     + ss(i,j,2)*(uu(i-1,j)-uu(i,j)) - ss(i+1,j,2)*(uu(i+1,j)-uu(i,j))
                if (bc_skewed(mm(i,j),1,+1)) &
                     flux(1,j) = flux(1,j) + ss(i,j,XBC)*(uu(i+2,j)-uu(i,j)) 
             else if (bc_neumann(mm(i,j),1,-1)) then
                flux(1,j) = ss(i,j,2)*uu(i-1,j)
             else   
                flux(1,j) = ss(i,j,2)*(uu(i-1,j)-uu(i,j))
             end if
          end do

!      Hi i face
       else if (face == 1) then

          i = nx
          flux(1,:) = ZERO
          do j = 1,ny
             if (bc_dirichlet(mm(i,j),1,+1)) then
                flux(1,j) = &
                       ss(i,j,1)*(uu(i+1,j)-uu(i,j)) &
                     + ss(i,j,2)*(uu(i-1,j)-uu(i,j)) - ss(i-1,j,1)*(uu(i-1,j)-uu(i,j))
                if (bc_skewed(mm(i,j),1,-1)) &
                     flux(1,j) = flux(1,j) + ss(i,j,XBC)*(uu(i-2,j)-uu(i,j))
             else if (bc_neumann(mm(i,j),1,+1)) then
                flux(1,j) = ss(i,j,1)*uu(i+1,j)
             else 
                flux(1,j) = ss(i,j,1)*(uu(i+1,j)-uu(i,j))
             end if
          end do

       end if

!   Lo j face
    else if ( dim == 2 ) then
       if (face == -1) then

          j = 1
          flux(:,1) = ZERO
          do i = 1,nx
             if (bc_dirichlet(mm(i,j),2,-1)) then
                flux(i,1) = &
                       ss(i,j,3)*(uu(i,j+1)-uu(i,j)) &
                     + ss(i,j,4)*(uu(i,j-1)-uu(i,j)) - ss(i,j+1,4)*(uu(i,j+1)-uu(i,j))
                if (bc_skewed(mm(i,j),2,+1)) &
                     flux(i,1) =  flux(i,1) + ss(i,j,YBC)*(uu(i,j+2)-uu(i,j))
             else if (bc_neumann(mm(i,j),2,-1)) then
                flux(i,1) = ss(i,j,4)*uu(i,j-1)
             else 
                flux(i,1) = ss(i,j,4)*(uu(i,j-1)-uu(i,j))
             end if
          end do


!      Hi j face
       else if (face == 1) then

          j = ny
          flux(:,1) = ZERO
          do i = 1,nx
             if (bc_dirichlet(mm(i,j),2,+1)) then
                flux(i,1) = &
                       ss(i,j,3)*(uu(i,j+1)-uu(i,j)) &
                     + ss(i,j,4)*(uu(i,j-1)-uu(i,j)) - ss(i,j-1,3)*(uu(i,j-1)-uu(i,j))
                if (bc_skewed(mm(i,j),2,-1)) &
                     flux(i,1) = flux(i,1) + ss(i,j,YBC)*(uu(i,j-2)-uu(i,j))
             else if (bc_neumann(mm(i,j),2,+1)) then
                flux(i,1) = ss(i,j,3)*uu(i,j+1)
             else
                flux(i,1) = ss(i,j,3)*(uu(i,j+1)-uu(i,j))
             end if
          end do

       end if
    end if

  end subroutine stencil_fine_flux_2d

  subroutine stencil_fine_flux_3d(ss, flux, uu, mm, ng, face, dim, skwd)
    integer, intent(in) :: ng
    real (kind = dp_t), intent(in ) :: uu(1-ng:,1-ng:,1-ng:)
    real (kind = dp_t), intent(out) :: flux(:,:,:)
    real (kind = dp_t), intent(in ) :: ss(:,:,:,0:)
    integer           , intent(in)  :: mm(:,:,:)
    logical, intent(in), optional :: skwd
    integer, intent(in) :: face, dim
    integer nx, ny, nz
    integer i,j,k
    integer, parameter :: XBC = 7, YBC = 8, ZBC = 9
    logical :: lskwd

    lskwd = .true. ; if ( present(skwd) ) lskwd = skwd

    nx = size(ss,dim=1)
    ny = size(ss,dim=2)
    nz = size(ss,dim=3)

    if ( dim ==  1 ) then
       !   Lo i face
       if (face == -1) then

          i = 1
          flux(1,:,:) = ZERO
          do k = 1,nz
             do j = 1,ny
                if (bc_dirichlet(mm(i,j,k),1,-1)) then
                   flux(1,j,k) =  &
                          ss(i,j,k,1)*(uu(i+1,j,k)-uu(i,j,k)) &
                        + ss(i,j,k,2)*(uu(i-1,j,k)-uu(i,j,k)) &
                        - ss(i+1,j,k,2)*(uu(i+1,j,k)-uu(i,j,k))
                   if (bc_skewed(mm(i,j,k),1,+1)) &
                        flux(1,j,k) =  flux(1,j,k) + ss(i,j,k,XBC)*(uu(i+2,j,k)-uu(i,j,k))
                else if (bc_neumann(mm(i,j,k),1,-1)) then
                   flux(1,j,k) = ss(i,j,k,2)*uu(i-1,j,k)
                else 
                   flux(1,j,k) = ss(i,j,k,2)*(uu(i-1,j,k)-uu(i,j,k))
                end if
             end do
          end do

       !   Hi i face
       else if (face ==  1) then

          i = nx
          flux(1,:,:) = ZERO
          do k = 1,nz
             do j = 1,ny
                if (bc_dirichlet(mm(i,j,k),1,+1)) then
                   flux(1,j,k) = &
                          ss(i,j,k,1)*(uu(i+1,j,k)-uu(i,j,k)) &
                        + ss(i,j,k,2)*(uu(i-1,j,k)-uu(i,j,k)) &
                        - ss(i-1,j,k,1)*(uu(i-1,j,k)-uu(i,j,k))
                   if (bc_skewed(mm(i,j,k),1,-1)) &
                        flux(1,j,k) =  flux(1,j,k) + ss(i,j,k,XBC)*(uu(i-2,j,k)-uu(i,j,k))
                else if (bc_neumann(mm(i,j,k),1,+1)) then
                   flux(1,j,k) = ss(i,j,k,1)*uu(i+1,j,k)
                else 
                   flux(1,j,k) = ss(i,j,k,1)*(uu(i+1,j,k)-uu(i,j,k))
                end if
             end do
          end do
       end if

    else if ( dim == 2 ) then

       !   Lo j face
       if (face == -1) then
          j = 1
          flux(:,1,:) = ZERO
          do k = 1,nz
             do i = 1,nx
                if (bc_dirichlet(mm(i,j,k),2,-1)) then
                   flux(i,1,k) = &
                          ss(i,j,k,3)*(uu(i,j+1,k)-uu(i,j,k)) &
                        + ss(i,j,k,4)*(uu(i,j-1,k)-uu(i,j,k)) &
                        - ss(i,j+1,k,4)*(uu(i,j+1,k)-uu(i,j,k))
                   if (bc_skewed(mm(i,j,k),2,+1)) &
                        flux(i,1,k) =  flux(i,1,k) + ss(i,j,k,YBC)*(uu(i,j+2,k)-uu(i,j,k))
                else if (bc_neumann(mm(i,j,k),2,-1)) then
                   flux(i,1,k) = ss(i,j,k,4)*uu(i,j-1,k)
                else 
                   flux(i,1,k) = ss(i,j,k,4)*(uu(i,j-1,k)-uu(i,j,k))
                end if
             end do
          end do

       !   Hi j face
       else if (face ==  1) then

          j = ny
          flux(:,1,:) = ZERO
          do k = 1,nz
             do i = 1,nx
                if (bc_dirichlet(mm(i,j,k),2,+1)) then
                   flux(i,1,k) =  &
                          ss(i,j,k,3)*(uu(i,j+1,k)-uu(i,j,k)) &
                        + ss(i,j,k,4)*(uu(i,j-1,k)-uu(i,j,k)) &
                        - ss(i,j-1,k,3)*(uu(i,j-1,k)-uu(i,j,k))
                   if (bc_skewed(mm(i,j,k),2,-1)) &
                        flux(i,1,k) =  flux(i,1,k) + ss(i,j,k,YBC)*(uu(i,j-2,k)-uu(i,j,k))
                else if (bc_neumann(mm(i,j,k),2,+1)) then
                   flux(i,1,k) = ss(i,j,k,3)*uu(i,j+1,k)
                else
                   flux(i,1,k) = ss(i,j,k,3)*(uu(i,j+1,k)-uu(i,j,k))
                end if
             end do
          end do
       end if

    else if ( dim == 3 ) then

       !   Lo k face
       if (face == -1) then

          k = 1
          flux(:,:,1) = ZERO
          do j = 1,ny
             do i = 1,nx
                if (bc_dirichlet(mm(i,j,k),3,-1)) then
                   flux(i,j,1) =  &
                          ss(i,j,k,5)*(uu(i,j,k+1)-uu(i,j,k)) &
                        + ss(i,j,k,6)*(uu(i,j,k-1)-uu(i,j,k)) &
                        - ss(i,j,k+1,6)*(uu(i,j,k+1)-uu(i,j,k))
                   if (bc_skewed(mm(i,j,k),3,+1)) &
                        flux(i,j,1) =  flux(i,j,1) + ss(i,j,k,ZBC)*(uu(i,j,k+2)-uu(i,j,k)) 
                else if (bc_neumann(mm(i,j,k),3,-1)) then
                   flux(i,j,1) = ss(i,j,k,6)*uu(i,j,k-1)
                else 
                   flux(i,j,1) = ss(i,j,k,6)*(uu(i,j,k-1)-uu(i,j,k))
                end if
             end do
          end do

       !   Hi k face
       else if (face ==  1) then

          k = nz
          flux(:,:,1) = ZERO
          do j = 1,ny
             do i = 1,nx
                if (bc_dirichlet(mm(i,j,k),3,+1)) then
                   flux(i,j,1) =  &
                          ss(i,j,k,5)*(uu(i,j,k+1)-uu(i,j,k)) &
                        + ss(i,j,k,6)*(uu(i,j,k-1)-uu(i,j,k)) &
                        - ss(i,j,k-1,5)*(uu(i,j,k-1)-uu(i,j,k))
                   if (bc_skewed(mm(i,j,k),3,-1)) &
                        flux(i,j,1) =  flux(i,j,1) + ss(i,j,k,ZBC)*(uu(i,j,k-2)-uu(i,j,k))
                else if (bc_neumann(mm(i,j,k),3,+1)) then
                   flux(i,j,1) = ss(i,j,k,5)*uu(i,j,k+1)
                else
                   flux(i,j,1) = ss(i,j,k,5)*(uu(i,j,k+1)-uu(i,j,k))
                end if
             end do
          end do

       end if
    end if

  end subroutine stencil_fine_flux_3d

end module stencil_module
