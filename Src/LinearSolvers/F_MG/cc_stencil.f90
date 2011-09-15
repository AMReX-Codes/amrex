module cc_stencil_module

  use bl_types
  use bc_module
  use bc_functions_module
  use multifab_module

  implicit none

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

  real(kind=dp_t), parameter, private :: ZERO = 0.0_dp_t
  real(kind=dp_t), parameter, private :: HALF = 0.5_dp_t
  real(kind=dp_t), parameter, private :: ONE  = 1.0_dp_t

  private :: stencil_bc_type, stencil_bndry_aaa

contains

  function stencil_norm_st(st, mask) result(r)
    type(stencil), intent(in) :: st
    type(lmultifab), intent(in), optional :: mask
    real(kind=dp_t) :: r
    r = stencil_norm(st%ss, mask)
  end function stencil_norm_st

  function stencil_norm(ss, mask) result(r)
    use bl_prof_module
    real(kind=dp_t) :: r
    type(multifab), intent(in) :: ss
    type(lmultifab), intent(in), optional :: mask
    integer :: i,j,k,n,b
    real(kind=dp_t) :: r1, sum_comps
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    logical, pointer :: lp(:,:,:,:)
    type(bl_prof_timer), save :: bpt

    call build(bpt, "st_norm")

    r1 = -Huge(r1)

    if ( present(mask) ) then
       do b = 1, nboxes(ss)
          if ( remote(ss,b) ) cycle
          sp => dataptr(ss, b)
          lp => dataptr(mask, b)
          !$OMP PARALLEL DO PRIVATE(i,j,k,n,sum_comps) REDUCTION(max : r1)
          do k = lbound(sp,dim=3), ubound(sp,dim=3)
             do j = lbound(sp,dim=2), ubound(sp,dim=2)
                do i = lbound(sp,dim=1), ubound(sp,dim=1)
                   if ( lp(i,j,k,1) ) then
                      sum_comps = ZERO
                      do n = lbound(sp,dim=4), ubound(sp,dim=4)
                         sum_comps = sum_comps + abs(sp(i,j,k,n))
                      end do
                      r1 = max(r1,sum_comps)
                   end if
                end do
             end do
          end do
          !$OMP END PARALLEL DO
       end do
    else

       do b = 1, nboxes(ss)
          if ( multifab_remote(ss,b) ) cycle
          sp => dataptr(ss, b)
          !$OMP PARALLEL DO PRIVATE(i,j,k,n,sum_comps) REDUCTION(max : r1)
          do k = lbound(sp,dim=3), ubound(sp,dim=3)
             do j = lbound(sp,dim=2), ubound(sp,dim=2)
                do i = lbound(sp,dim=1), ubound(sp,dim=1)
                   sum_comps = ZERO
                   do n = lbound(sp,dim=4), ubound(sp,dim=4)
                      sum_comps = sum_comps + abs(sp(i,j,k,n))
                   end do
                   r1 = max(r1,sum_comps)
                end do
             end do
          end do
          !$OMP END PARALLEL DO
       end do
    end if

    call parallel_reduce(r,r1,MPI_MAX)
    call destroy(bpt)
  end function stencil_norm

  function max_of_stencil_sum(ss, mask) result(r)
    use bl_prof_module
    real(kind=dp_t) :: r
    type(multifab), intent(in) :: ss
    type(lmultifab), intent(in), optional :: mask
    integer :: i,j,k,n,b
    real(kind=dp_t) :: r1, sum_comps
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    logical, pointer :: lp(:,:,:,:)
    type(bl_prof_timer), save :: bpt

    ! NOTE: this is exactly the same as the stencil_norm function except that we sum the
    !       components of the stencil, not the absolute value of each component

    call build(bpt, "st_sum")
    r1 = -Huge(r1)
    if ( present(mask) ) then
       do b = 1, nboxes(ss)
          if ( remote(ss,b) ) cycle
          sp => dataptr(ss, b)
          lp => dataptr(mask, b)
          do k = lbound(sp,dim=3), ubound(sp,dim=3)
             do j = lbound(sp,dim=2), ubound(sp,dim=2)
                do i = lbound(sp,dim=1), ubound(sp,dim=1)
                   if ( lp(i,j,k,1) ) then
                      sum_comps = ZERO
                      do n = lbound(sp,dim=4), ubound(sp,dim=4)
                         sum_comps = sum_comps + sp(i,j,k,n)
                      end do
                      r1 = max(r1,sum_comps)
                   end if
                end do
             end do
          end do

       end do
    else
       do b = 1, nboxes(ss)
          if ( multifab_remote(ss,b) ) cycle
          sp => dataptr(ss, b)
          do k = lbound(sp,dim=3), ubound(sp,dim=3)
             do j = lbound(sp,dim=2), ubound(sp,dim=2)
                do i = lbound(sp,dim=1), ubound(sp,dim=1)
                   sum_comps = ZERO
                   do n = lbound(sp,dim=4), ubound(sp,dim=4)
                      sum_comps = sum_comps + sp(i,j,k,n)
                   end do
                   r1 = max(r1,sum_comps)
                end do
             end do
          end do
       end do
    end if

    call parallel_reduce(r,r1,MPI_MAX)
    call destroy(bpt)
  end function max_of_stencil_sum

  subroutine stencil_set_extrap_bc(st, max_order)
    type(stencil), intent(inout) :: st
    integer, intent(in) :: max_order
    st%extrap_bc = .true.
    st%extrap_max_order = max_order
  end subroutine stencil_set_extrap_bc

  subroutine stencil_print(st, str, unit, skip)
    use bl_IO_module
    type(stencil), intent(in) :: st
    character (len=*), intent(in), optional :: str
    integer, intent(in), optional :: unit
    integer, intent(in), optional :: skip
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

  subroutine stencil_set_bc(st, idx, mask, bc_face, cf_face)
    type(multifab),  intent(in)           :: st
    integer,         intent(in)           :: idx
    type(imultifab), intent(inout)        :: mask
    integer,         intent(in)           :: bc_face(:,:)
    integer,         intent(in), optional :: cf_face(:,:)

    type(box)        :: bx1, src, pd
    type(boxarray)   :: ba, sba
    integer          :: i, j, ii, jj, k, ldom
    integer, pointer :: mp(:,:,:,:)
    integer          :: lcf_face(size(bc_face, 1), size(bc_face, 2))
    logical          :: pmask(get_dim(st))
    !
    ! The Coarse-Fine boundary is Dirichlet unless specified.
    !
    lcf_face = BC_DIR; if ( present(cf_face) ) lcf_face = cf_face
    !
    ! Initialize every border to Fine-Fine (Interior).
    !
    mp => dataptr(mask, idx)

    mp = BC_INT

    pd    = get_pd(get_layout(st))
    pmask = get_pmask(get_layout(st))

    do i = 1, get_dim(st)
       if ( bc_face(i,1) == BC_PER .and. ( bc_face(i,1) /= bc_face(i,2) )) then
          call bl_error("STENCIL_SET_BC: confusion in bc_face")
       end if
       do j = -1, 1, 2
          bx1 = shift(get_box(st, idx), j, i)
          jj = (3 + j)/2
          if ( contains(pd, bx1) ) then
             !
             ! We're not touching a physical boundary -- set any/all C-F bndrys.
             !
             call boxarray_boxarray_diff(ba, bx1, get_boxarray(st))
             do ii = 1, nboxes(ba)
                bx1 = shift(get_box(ba,ii), -j, i)
                mp => dataptr(mask, idx, bx1)
                mp = ibset(mp, BC_BIT(lcf_face(i, jj), i, j))
             end do
             call destroy(ba)
          else
             !
             ! We touch a physical boundary in that direction.
             !
             if ( .not. pmask(i) ) then
                !
                ! We're not periodic in that direction -- use physical BCs.
                !
                call boxarray_box_diff(ba, bx1, pd)
                do ii = 1, nboxes(ba)
                   bx1 = shift(get_box(ba,ii), -j, i)
                   mp => dataptr(mask, idx, bx1)
                   mp = ibset(mp, BC_BIT(bc_face(i, jj), i, j))
                end do
                call destroy(ba)
             else
                !
                ! Remove any/all Fine-Fine intersections.
                !
                ldom = extent(pd, i)
                call boxarray_build_bx(ba, bx1)
                do k = 1, nboxes(st)
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
                do ii = 1, nboxes(ba)
                   bx1 = shift(get_box(ba,ii), -j, i)
                   mp => dataptr(mask, idx, bx1)
                   mp = ibset(mp, BC_BIT(lcf_face(i, jj), i, j))
                end do
                call destroy(ba)
             end if
          end if
       end do
    end do

  end subroutine stencil_set_bc

  elemental function stencil_bc_type(mask, dir, face) result(r)
    integer, intent(in) :: mask, dir, face
    integer :: r
    r = BC_INT
    if      ( bc_dirichlet(mask,dir,face) ) then
       r = BC_DIR
    else if ( bc_neumann  (mask,dir,face) ) then
       r = bc_NEU
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

    logical, parameter :: old_old = .TRUE.

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

    s0 = ZERO
    ss = ZERO
    sm = ZERO
    sp = ZERO
    !
    ! TODO -- this stuff is just not quite right.
    ! Some of this logic needs to be moved into the bc_?? routines themselves.
    !
    if ( nx > 1 ) bchi = BC_INT
    imaxo = maxo
    if ( nx == 1 ) imaxo = 1
    if ( nx == 2 ) imaxo = min(imaxo,2)

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
      call bl_error("STENCIL_BNDRY_AAA: should never reach bc_ii")
!     sm  = b0
!     s0  = -(b0+b1)
!     sp  = b1
!     ss  = ZERO
!     skewed = .false.
    end subroutine bc_ii

    subroutine bc_id
      if ( nx > 1 ) then
         call bc_ii
      else
         sm =  b0 + ( -1 + 4/(3 + 2*xb)) * b1
         s0 = -b0 + (( -3 + 2*xb )/(1 + 2*xb)) * b1
         sp =  8*b1/(3 + 4*xb*(2 + xb))
         ss = ZERO
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
         ss = ZERO
         skewed = .false.
      end if
    end subroutine bc_in

    subroutine bc_di
      select case (imaxo)
      case (1)
         sm = 2*b0/(1 + 2*xa)
         s0 = -2*b0/(1 + 2*xa) - b1
         sp = b1
         ss = ZERO
         skewed = .false.
      case (2)
         sm = 8*b0/(3 + 4*xa*(2 + xa))
         s0 = ((-3 + 2*xa)/(1 + 2*xa))*b0 - b1
         sp = ((1-2*xa)/(3 + 2*xa))*b0    + b1
         ss = ZERO
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
         ss = ZERO
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
         ss = ZERO
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
!        sm = -b0
         sm = ZERO
         s0 = -b1
         sp =  b1
         ss = ZERO
         skewed = .false.
      case (2)
!        sm = -b0/(1 + xa)
         sm = ZERO
         s0 = xa*b0/(1 + xa) - b1
         sp = -xa*b0/(1 + xa) + b1
         ss = ZERO
         skewed = .false.
      case (3)
!        sm = -24*b0/(23 + 12*xa*(3+xa))
!        s0 = 2*((-1 + 12*xa*(2+xa))/(23 + 12*xa*(3+xa)))*b0 - b1
!        sp = -3*((-1 + 4*xa*(5+3*xa))/(23 + 12*xa*(3+xa)))*b0 + b1
!        ss = ((-1 + 12*xa*(1+xa))/(23 + 12*xa*(3+xa)))*b0
!        skewed = .true.

         ! NOTE: we cant do anything higher-order for Neumann or we will lose solvability
         sm = -b0/(1 + xa)
         s0 = xa*b0/(1 + xa) - b1
         sp = -xa*b0/(1 + xa) + b1
         ss = ZERO
         skewed = .false.
      end select
    end subroutine bc_ni

    subroutine bc_nd
      select case ( imaxo )
      case (1)
         sm = - ((3+2*xb)*b0 + (1-2*xb)*b1)/(3+4*xa+2*xb)
         s0 = 8*(xa*b0 -(1+xa)*b1)/((1+2*xb)*(3+4*xa+2*xb))
         sp = 8*(-xa*b0 + (1+xb)*b1)/((1+2*xb)*(3+4*xa+2*xb))
         ss = ZERO
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
         s0 = ZERO
         sp = (-xa*b0 + (1+xa)*b1)/(1+xa+xb)
         ss = ZERO
         skewed = .false.
      case (2)
         sm = -(1+xb)*b0/(1+xa+xb)
         s0 = -b1
         sp =  b1
         ss = -xa*b0/(1+xa+xb)
         skewed = .true.
      case (3)
!        sm = -(23+12*xb*(3+xb))*b0/((2+xa+xb)*(11+12*xb+12*xa*(1+xb)))
!        s0 = (-1 + 12*xa*(2*xb))*b0/((11+12*xb+12*xa*(1+xb))) - b1
!        sp = -(-1 + 12*xa*(2*xb))*b0/((11+12*xb+12*xa*(1+xb))) + b1
!        ss = (-1 + 12*xa*(1+xa))*b0/((2+xa*xb)*(11+12*xb + 12*xa*(1+xb)))
!        skewed = .true.

         ! NOTE: we cant do anything higher-order for Neumann or we will lose solvability
         sm = -(1+xb)*b0/(1+xa+xb)
         s0 = -b1
         sp =  b1
         ss = -xa*b0/(1+xa+xb)
         skewed = .true.
      end select
    end subroutine bc_nn

  end subroutine stencil_bndry_aaa
  
  subroutine s_simple_1d_cc(ss, alpha, ng_a, betax, ng_b, dh, mask, lo, hi, xa, xb, order)

    integer           , intent(in   ) :: ng_a, ng_b, lo(:), hi(:), order
    integer           , intent(inout) :: mask(lo(1):)
    real (kind = dp_t), intent(  out) ::   ss(lo(1)  :,0:)
    real (kind = dp_t), intent(in   ) :: alpha(lo(1)-ng_a:)
    real (kind = dp_t), intent(in   ) :: betax(lo(1)-ng_b:)
    real (kind = dp_t), intent(in   ) :: dh(:)
    real (kind = dp_t), intent(in   ) :: xa(:), xb(:)

    real (kind = dp_t) :: f1(1)
    integer            :: i,bclo,bchi,nx
    integer, parameter :: XBC = 3

    nx = hi(1)-lo(1)+1 
    f1 = ONE/dh**2

    mask = ibclr(mask, BC_BIT(BC_GEOM,1,-1))
    mask = ibclr(mask, BC_BIT(BC_GEOM,1,+1))

    do i = lo(1),hi(1)
       ss(i,0) =   ZERO
       ss(i,1) = -betax(i+1)*f1(1)
       ss(i,2) = -betax(i  )*f1(1)
       ss(i,XBC) = ZERO
    end do

    ! x derivatives

    do i = lo(1)+1, hi(1)-1
       ss(i,0) = ss(i,0) + (betax(i+1)+betax(i))*f1(1)
    end do

    bclo = stencil_bc_type(mask(lo(1)),1,-1)
    bchi = stencil_bc_type(mask(hi(1)),1,+1)

    i = lo(1)
    if (bclo .eq. BC_INT) then
       ss(i,0) = ss(i,0) + (betax(i)+betax(i+1))*f1(1)
    else
       call stencil_bndry_aaa(order, nx, 1, -1, mask(i), &
            ss(i,0), ss(i,1), ss(i,2), ss(i,XBC), &
            betax(i), betax(i+1), xa(1), xb(1), dh(1), bclo, bchi)
    end if

    if ( hi(1) > lo(1) ) then
       i = hi(1)
       if (bchi .eq. BC_INT) then
          ss(i,0) = ss(i,0) + (betax(i)+betax(i+1))*f1(1)
       else
          call stencil_bndry_aaa(order, nx, 1, 1, mask(i), &
               ss(i,0), ss(i,1), ss(i,2), ss(i,XBC), &
               betax(i), betax(i+1), xa(1), xb(1), dh(1), bclo, bchi)
       end if
    end if

    do i = lo(1),hi(1)
       ss(i,0) = ss(i,0) + alpha(i)
    end do

  end subroutine s_simple_1d_cc

  subroutine s_simple_2d_cc(ss, alpha, ng_a, betax, betay, ng_b, dh, mask, lo, hi, xa, xb, order)

    integer           , intent(in   ) :: ng_a, ng_b, lo(:), hi(:), order
    integer           , intent(inout) :: mask(lo(1)  :,lo(2)  :)
    real (kind = dp_t), intent(  out) ::   ss(lo(1)  :,lo(2)  :,0:)
    real (kind = dp_t), intent(in   ) :: alpha(lo(1)-ng_a:,lo(2)-ng_a:)
    real (kind = dp_t), intent(in   ) :: betax(lo(1)-ng_b:,lo(2)-ng_b:)
    real (kind = dp_t), intent(in   ) :: betay(lo(1)-ng_b:,lo(2)-ng_b:)
    real (kind = dp_t), intent(in   ) :: dh(:)
    real (kind = dp_t), intent(in   ) :: xa(:), xb(:)

    real (kind = dp_t) :: f1(2)
    integer            :: i, j, bclo, bchi, nx, ny
    integer, parameter :: XBC = 5, YBC = 6

    nx = hi(1)-lo(1)+1
    ny = hi(2)-lo(2)+1

    f1 = ONE/dh**2

    mask = ibclr(mask, BC_BIT(BC_GEOM,1,-1))
    mask = ibclr(mask, BC_BIT(BC_GEOM,1,+1))
    mask = ibclr(mask, BC_BIT(BC_GEOM,2,-1))
    mask = ibclr(mask, BC_BIT(BC_GEOM,2,+1))

    do j = lo(2),hi(2)
       do i = lo(1),hi(1)
          ss(i,j,0) = ZERO
          ss(i,j,1) = -betax(i+1,j)*f1(1)
          ss(i,j,2) = -betax(i  ,j)*f1(1)
          ss(i,j,3) = -betay(i,j+1)*f1(2)
          ss(i,j,4) = -betay(i,j  )*f1(2)
          ss(i,j,XBC) = ZERO
          ss(i,j,YBC) = ZERO
       end do
    end do

    ! x derivatives

    do j = lo(2),hi(2)
       do i = lo(1)+1,hi(1)-1
          ss(i,j,0) = ss(i,j,0) + (betax(i,j)+betax(i+1,j))*f1(1)
       end do
    end do

    do j = lo(2),hi(2)
       bclo = stencil_bc_type(mask(lo(1),j),1,-1)
       bchi = stencil_bc_type(mask(hi(1),j),1,+1)
 
       i = lo(1)
       if (bclo .eq. BC_INT) then
          ss(i,j,0) = ss(i,j,0) + (betax(i,j)+betax(i+1,j))*f1(1)
       else
          call stencil_bndry_aaa(order, nx, 1, -1, mask(i,j), &
               ss(i,j,0), ss(i,j,1), ss(i,j,2), ss(i,j,XBC), &
               betax(i,j), betax(i+1,j), &
               xa(1), xb(1), dh(1), bclo, bchi)
       end if

       if ( hi(1) > lo(1) ) then
          i = hi(1)
          if (bchi .eq. BC_INT) then
             ss(i,j,0) = ss(i,j,0) + (betax(i,j)+betax(i+1,j))*f1(1)
          else
             call stencil_bndry_aaa(order, nx, 1, 1, mask(i,j), &
                  ss(i,j,0), ss(i,j,1), ss(i,j,2), ss(i,j,XBC), &
                  betax(i,j), betax(i+1,j), &
                  xa(1), xb(1), dh(1), bclo, bchi) 
          end if 
       end if
    end do

    ! y derivatives

    do i = lo(1),hi(1)
       do j = lo(2)+1,hi(2)-1
          ss(i,j,0) = ss(i,j,0) + (betay(i,j)+betay(i,j+1))*f1(2)
       end do
    end do

    do i = lo(1),hi(1)
       bclo = stencil_bc_type(mask( i,lo(2)),2,-1)
       bchi = stencil_bc_type(mask( i,hi(2)),2,+1)

       j = lo(2)
       if (bclo .eq. BC_INT) then
          ss(i,j,0) = ss(i,j,0) + (betay(i,j)+betay(i,j+1))*f1(2)
       else
          call stencil_bndry_aaa(order, ny, 2, -1, mask(i,j), &
               ss(i,j,0), ss(i,j,3), ss(i,j,4),ss(i,j,YBC), &
               betay(i,j), betay(i,j+1), &
               xa(2), xb(2), dh(2), bclo, bchi)
       end if

       if ( hi(2) > lo(2) ) then
          j = hi(2)
          if (bchi .eq. BC_INT) then
             ss(i,j,0) = ss(i,j,0) + (betay(i,j)+betay(i,j+1))*f1(2)
          else
             call stencil_bndry_aaa(order, ny, 2, 1, mask(i,j), &
                  ss(i,j,0), ss(i,j,3), ss(i,j,4), ss(i,j,YBC), &
                  betay(i,j), betay(i,j+1), &
                  xa(2), xb(2), dh(2), bclo, bchi)
          end if
       end if
    end do

    do j = lo(2),hi(2)
       do i = lo(1),hi(1)
          ss(i,j,0) = ss(i,j,0) + alpha(i,j)
       end do
    end do

  end subroutine s_simple_2d_cc

  subroutine s_simplen_2d_cc(ss, alpha, ng_a, betax, betay, ng_b, dh, mask, lo, hi, xa, xb, order)

    integer           , intent(in   ) :: ng_a, ng_b, lo(:), hi(:), order
    integer           , intent(inout) :: mask(lo(1)  :,lo(2)  :)
    real (kind = dp_t), intent(  out) ::   ss(lo(1)  :,lo(2)  :,0:)
    real (kind = dp_t), intent(in   ) :: alpha(lo(1)-ng_a:,lo(2)-ng_a:,0:)
    real (kind = dp_t), intent(in   ) :: betax(lo(1)-ng_b:,lo(2)-ng_b:,:)
    real (kind = dp_t), intent(in   ) :: betay(lo(1)-ng_b:,lo(2)-ng_b:,:)
    real (kind = dp_t), intent(in   ) :: dh(:)
    real (kind = dp_t), intent(in   ) :: xa(:), xb(:)

    real (kind = dp_t) :: f1(2), blo, bhi
    integer            :: i, j, dm, n, bclo, bchi, nx, ny, nc
    integer, parameter :: XBC = 5, YBC = 6

    nx = hi(1)-lo(1)+1
    ny = hi(2)-lo(2)+1

    dm = 2
    nc = size(betax,dim=3)
    f1 = ONE/dh**2

    mask = ibclr(mask, BC_BIT(BC_GEOM,1,-1))
    mask = ibclr(mask, BC_BIT(BC_GEOM,1,+1))
    mask = ibclr(mask, BC_BIT(BC_GEOM,2,-1))
    mask = ibclr(mask, BC_BIT(BC_GEOM,2,+1))
 
    ss(:,:,:) = 0.d0

    ! Consider the operator  ( alpha - sum_n (beta0_n del dot beta_n grad) )
    ! Components alpha(i,j,   0) = alpha
    ! Components alpha(i,j,1:nc) = beta0_n
    ! Components betax(i,j,1:nc) = betax_n
    ! Components betay(i,j,1:nc) = betay_n

    ! ss(i,j,1) is the coefficient of phi(i+1,j  )
    ! ss(i,j,2) is the coefficient of phi(i-1,j  )
    ! ss(i,j,3) is the coefficient of phi(i  ,j+1)
    ! ss(i,j,4) is the coefficient of phi(i  ,j-1)
    ! ss(i,j,0) is the coefficient of phi(i  ,j  )
    
    do j = lo(2),hi(2)
       do i = lo(1),hi(1)
          do n = 1,nc
             ss(i,j,1) = ss(i,j,1) - betax(i+1,j,n)*f1(1) / alpha(i,j,n)
             ss(i,j,2) = ss(i,j,2) - betax(i  ,j,n)*f1(1) / alpha(i,j,n)
             ss(i,j,3) = ss(i,j,3) - betay(i,j+1,n)*f1(2) / alpha(i,j,n)
             ss(i,j,4) = ss(i,j,4) - betay(i,j  ,n)*f1(2) / alpha(i,j,n)
          end do 
       end do
    end do

    ! x derivatives

    do j = lo(2),hi(2)
       do i = lo(1)+1,hi(1)-1
          do n = 1, nc
             ss(i,j,0) = ss(i,j,0) + (betax(i,j,n)+betax(i+1,j,n))*f1(1) / alpha(i,j,n)
          end do
       end do
    end do

    do j = lo(2),hi(2)
       bclo = stencil_bc_type(mask(lo(1),j),1,-1)
       bchi = stencil_bc_type(mask(hi(1),j),1,+1)
 
       i = lo(1)
       if (bclo .eq. BC_INT) then
          do n = 1, nc
             ss(i,j,0) = ss(i,j,0) + (betax(i,j,n)+betax(i+1,j,n))*f1(1) / alpha(i,j,n)
          end do
       else
          blo = 0.d0
          bhi = 0.d0
          do n = 1,nc
            blo = blo + betax(i  ,j,n) / alpha(i,j,n)
            bhi = bhi + betax(i+1,j,n) / alpha(i,j,n)
          end do
          call stencil_bndry_aaa(order, nx, 1, -1, mask(i,j), &
               ss(i,j,0), ss(i,j,1), ss(i,j,2), ss(i,j,XBC), &
               blo, bhi, xa(1), xb(1), dh(1), bclo, bchi)
       end if

       if ( hi(1) > lo(1) ) then
          i = hi(1)
          if (bchi .eq. BC_INT) then
             do n = 1, nc
                ss(i,j,0) = ss(i,j,0) + (betax(i,j,n)+betax(i+1,j,n))*f1(1) / alpha(i,j,n)
             end do
          else
             blo = 0.d0
             bhi = 0.d0
             do n = 1,nc
                blo = blo + betax(i  ,j,n) / alpha(i,j,n)
                bhi = bhi + betax(i+1,j,n) / alpha(i,j,n)
             end do
             call stencil_bndry_aaa(order, nx, 1, 1, mask(i,j), &
                  ss(i,j,0), ss(i,j,1), ss(i,j,2), ss(i,j,XBC), &
                  blo, bhi, xa(1), xb(1), dh(1), bclo, bchi)
          end if
       end if
    end do

    ! y derivatives

    do i = lo(1),hi(1)
       do j = lo(2)+1,hi(2)-1
          do n = 1,nc
             ss(i,j,0) = ss(i,j,0) + (betay(i,j,n)+betay(i,j+1,n))*f1(2) / alpha(i,j,n)
          end do
       end do
    end do

    do i = lo(1),hi(1)
       bclo = stencil_bc_type(mask( i,lo(2)),2,-1)
       bchi = stencil_bc_type(mask( i,hi(2)),2,+1)

       j = lo(2)
       if (bclo .eq. BC_INT) then
          do n = 1,nc
             ss(i,j,0) = ss(i,j,0) + (betay(i,j,n)+betay(i,j+1,n))*f1(2) / alpha(i,j,n)
          end do
       else
          blo = 0.d0
          bhi = 0.d0
          do n = 1,nc
             blo = blo + betay(i  ,j,n) / alpha(i,j,n) 
             bhi = bhi + betay(i,j+1,n) / alpha(i,j,n) 
          end do
          call stencil_bndry_aaa(order, ny, 2, -1, mask(i,j), &
               ss(i,j,0), ss(i,j,3), ss(i,j,4),ss(i,j,YBC), &
               blo, bhi, xa(2), xb(2), dh(2), bclo, bchi)
       end if

       if ( hi(2) > lo(2) ) then
          j = hi(2)
          if (bchi .eq. BC_INT) then
             do n = 1,nc
                ss(i,j,0) = ss(i,j,0) + (betay(i,j,n)+betay(i,j+1,n))*f1(2) / alpha(i,j,n)
             end do
          else
             blo = 0.d0
             bhi = 0.d0
             do n = 1,nc
                blo = blo + betay(i  ,j,n) / alpha(i,j,n) 
                bhi = bhi + betay(i,j+1,n) / alpha(i,j,n) 
             end do
             call stencil_bndry_aaa(order, ny, 2, 1, mask(i,j), &
                  ss(i,j,0), ss(i,j,3), ss(i,j,4), ss(i,j,YBC), &
                  blo, bhi, xa(2), xb(2), dh(2), bclo, bchi)
          end if
       end if
    end do

    do j = lo(2),hi(2)
       do i = lo(1),hi(1)
          ss(i,j,0) = ss(i,j,0) + alpha(i,j,0) 
       end do
    end do

  end subroutine s_simplen_2d_cc

 subroutine s_simplem_2d_cc(ss, alpha, ng_a, betax, betay, ng_b, dh, mask, lo, hi, xa, xb, order)

    integer           , intent(in   ) :: ng_a, ng_b, lo(:), hi(:), order
    integer           , intent(inout) :: mask(lo(1)  :,lo(2)  :)
    real (kind = dp_t), intent(  out) :: ss(lo(1)  :,lo(2)  :,0:)
    real (kind = dp_t), intent(in   ) :: alpha(lo(1)-ng_a:,lo(2)-ng_a:,0:)
    real (kind = dp_t), intent(in   ) :: betax(lo(1)-ng_b:,lo(2)-ng_b:,:)
    real (kind = dp_t), intent(in   ) :: betay(lo(1)-ng_b:,lo(2)-ng_b:,:)
    real (kind = dp_t), intent(in   ) :: dh(:)
    real (kind = dp_t), intent(in   ) :: xa(:), xb(:)

    real (kind = dp_t) :: f1(2), blo, bhi
    integer            :: i, j, dm, bclo, bchi, nx, ny, nc
    integer, parameter :: XBC = 5, YBC = 6

    nx = hi(1)-lo(1)+1
    ny = hi(2)-lo(2)+1

    dm = 2
    nc = size(betax,dim=3)
    f1 = ONE/dh**2

    mask = ibclr(mask, BC_BIT(BC_GEOM,1,-1))
    mask = ibclr(mask, BC_BIT(BC_GEOM,1,+1))
    mask = ibclr(mask, BC_BIT(BC_GEOM,2,-1))
    mask = ibclr(mask, BC_BIT(BC_GEOM,2,+1))
 
    ss(:,:,:) = 0.d0

    ! Consider the operator  ( alpha - sum_n (beta0_n del dot beta_n grad) )
    ! Components alpha(i,j,   0) = alpha
    ! Components alpha(i,j,1:nc) = beta0_n
    ! Components betax(i,j,1:nc) = betax_n
    ! Components betay(i,j,1:nc) = betay_n

    ! ss(i,j,1) is the coefficient of phi(i+1,j  )
    ! ss(i,j,2) is the coefficient of phi(i-1,j  )
    ! ss(i,j,3) is the coefficient of phi(i  ,j+1)
    ! ss(i,j,4) is the coefficient of phi(i  ,j-1)
    ! ss(i,j,0) is the coefficient of phi(i  ,j  )
    
    do j = lo(2),hi(2)
       do i = lo(1),hi(1)
          ss(i,j,1) = ss(i,j,1) - (betax(i+1,j,1)+betax(i+1,j,2))*f1(1) 
          ss(i,j,2) = ss(i,j,2) - (betax(i  ,j,1)-betax(i,  j,2))*f1(1) 
          ss(i,j,3) = ss(i,j,3) - (betay(i,j+1,1)+betay(i,j+1,2))*f1(2) 
          ss(i,j,4) = ss(i,j,4) - (betay(i,j  ,1)-betay(i,j  ,2))*f1(2) 
       end do
    end do

    ! x derivatives

    do j = lo(2),hi(2)
       do i = lo(1)+1,hi(1)-1
            ss(i,j,0) = ss(i,j,0) - ss(i,j,1) - ss(i,j,2) 
       end do
    end do

    do j = lo(2),hi(2)
       bclo = stencil_bc_type(mask(lo(1),j),1,-1)
       bchi = stencil_bc_type(mask(hi(1),j),1,+1)
 
       i = lo(1)
       if (bclo .eq. BC_INT) then
          ss(i,j,0) = ss(i,j,0) - ss(i,j,1) - ss(i,j,2)
       else
          blo = -ss(i,j,2)/f1(1)  
          bhi = -ss(i,j,1)/f1(1)
          call stencil_bndry_aaa(order, nx, 1, -1, mask(i,j), &
               ss(i,j,0), ss(i,j,1), ss(i,j,2), ss(i,j,XBC), &
               blo, bhi, xa(1), xb(1), dh(1), bclo, bchi)
       end if

       if ( hi(1) > lo(1) ) then
          i = hi(1)
          if (bchi .eq. BC_INT) then
             ss(i,j,0) = ss(i,j,0) - ss(i,j,1) - ss(i,j,2) 
          else
             blo = -ss(i,j,2)/f1(1)  
             bhi = -ss(i,j,1)/f1(1)
             call stencil_bndry_aaa(order, nx, 1, 1, mask(i,j), &
                  ss(i,j,0), ss(i,j,1), ss(i,j,2), ss(i,j,XBC), &
                  blo, bhi, xa(1), xb(1), dh(1), bclo, bchi)
          end if
       end if
    end do

    ! y derivatives

    do i = lo(1),hi(1)
       do j = lo(2)+1,hi(2)-1
             ss(i,j,0) = ss(i,j,0) - ss(i,j,3) - ss(i,j,4) 
       end do
    end do

    do i = lo(1),hi(1)
       bclo = stencil_bc_type(mask( i,lo(2)),2,-1)
       bchi = stencil_bc_type(mask( i,hi(2)),2,+1)

       j = lo(2)
       if (bclo .eq. BC_INT) then
          ss(i,j,0) = ss(i,j,0) - ss(i,j,3) - ss(i,j,4) 
       else
          blo = -ss(i,j,4)/f1(2)  
          bhi = -ss(i,j,3)/f1(2)         
          call stencil_bndry_aaa(order, ny, 2, -1, mask(i,j), &
               ss(i,j,0), ss(i,j,3), ss(i,j,4),ss(i,j,YBC), &
               blo, bhi, xa(2), xb(2), dh(2), bclo, bchi)
       end if

       if ( hi(2) > lo(2) ) then
          j = hi(2)
          if (bchi .eq. BC_INT) then
             ss(i,j,0) = ss(i,j,0) - ss(i,j,3) - ss(i,j,4) 
          else
             blo = -ss(i,j,4)/f1(2)  
             bhi = -ss(i,j,3)/f1(2)
             call stencil_bndry_aaa(order, ny, 2, 1, mask(i,j), &
                  ss(i,j,0), ss(i,j,3), ss(i,j,4), ss(i,j,YBC), &
                  blo, bhi, xa(2), xb(2), dh(2), bclo, bchi)
          end if
       end if
    end do

    do j = lo(2),hi(2)
       do i = lo(1),hi(1)
          ss(i,j,0) = ss(i,j,0) + alpha(i,j,0) 
       end do
    end do

  end subroutine s_simplem_2d_cc

  subroutine s_simpleg_2d_cc(ss, alpha, ng_a, betax, betay, ng_b, dh, mask, lo, hi)

    integer           , intent(in   ) :: ng_a, ng_b, lo(:), hi(:)
    integer           , intent(inout) :: mask(lo(1)  :,lo(2)  :)
    real (kind = dp_t), intent(  out) :: ss(lo(1)  :,lo(2)  :,0:)
    real (kind = dp_t), intent(in   ) :: alpha(lo(1)-ng_a:,lo(2)-ng_a:,0:)
    real (kind = dp_t), intent(in   ) :: betax(lo(1)-ng_b:,lo(2)-ng_b:,:)
    real (kind = dp_t), intent(in   ) :: betay(lo(1)-ng_b:,lo(2)-ng_b:,:)
    real (kind = dp_t), intent(in   ) :: dh(:)

    real (kind = dp_t) :: f1(2)
    integer            :: i, j, dm, bclo, bchi, nx, ny, nc
    integer, parameter :: XBC = 5, YBC = 6

    nx = hi(1)-lo(1)+1
    ny = hi(2)-lo(2)+1

    dm = 2
    nc = size(betax,dim=3)
    f1 = ONE/dh**2

    mask = ibclr(mask, BC_BIT(BC_GEOM,1,-1))
    mask = ibclr(mask, BC_BIT(BC_GEOM,1,+1))
    mask = ibclr(mask, BC_BIT(BC_GEOM,2,-1))
    mask = ibclr(mask, BC_BIT(BC_GEOM,2,+1))
 
    ss(:,:,:) = 0.d0

    ! Consider the operator  ( alpha - sum_n (beta0_n del dot beta_n grad) )
    ! Components alpha(i,j,   0) = alpha
    ! Components alpha(i,j,1:nc) = beta0_n
    ! Components betax(i,j,1:nc) = betax_n
    ! Components betay(i,j,1:nc) = betay_n

    ! ss(i,j,1) is the coefficient of phi(i+1,j  )
    ! ss(i,j,2) is the coefficient of phi(i-1,j  )
    ! ss(i,j,3) is the coefficient of phi(i  ,j+1)
    ! ss(i,j,4) is the coefficient of phi(i  ,j-1)
    ! ss(i,j,0) is the coefficient of phi(i  ,j  )
    
    do j = lo(2),hi(2)
       do i = lo(1),hi(1)
          ss(i,j,1) = ss(i,j,1) - betax(i+1,j,1)
          ss(i,j,2) = ss(i,j,2) - betax(i  ,j,2)
          ss(i,j,3) = ss(i,j,3) - betay(i,j+1,1) 
          ss(i,j,4) = ss(i,j,4) - betay(i,j  ,2) 
       end do
    end do

    ! x derivatives

    do j = lo(2),hi(2)
       do i = lo(1)+1,hi(1)-1
          ss(i,j,0) = ss(i,j,0) - betax(i,j,3)
       end do
    end do

    do j = lo(2),hi(2)
       bclo = stencil_bc_type(mask(lo(1),j),1,-1)
       bchi = stencil_bc_type(mask(hi(1),j),1,+1)
 
       i = lo(1)
       if (bclo .eq. BC_INT) then
          ss(i,j,0) = ss(i,j,0) - betax(i,j,3)
       elseif (bclo .eq. BC_NEU) then
          ss(i,j,0) = ss(i,j,0) - betax(i,j,3) - betax(i,j,2)
          ss(i,j,2) = 0.d0
          ss(i,j,XBC) = 0.d0
       elseif (bclo .eq. BC_DIR) then
          ss(i,j,0) = ss(i,j,0) - betax(i,j,3)
          ss(i,j,2) = 0.d0
          ss(i,j,XBC) = 0.d0
       end if

       if ( hi(1) > lo(1) ) then
          i = hi(1)
          if (bchi .eq. BC_INT) then
             ss(i,j,0) = ss(i,j,0) - betax(i,j,3)
          elseif (bchi .eq. BC_NEU) then
             ss(i,j,0) = ss(i,j,0) - betax(i,j,3) - betax(i+1,j,1)
             ss(i,j,1) = 0.d0
             ss(i,j,XBC) = 0.d0
          elseif (bchi .eq. BC_DIR) then
             ss(i,j,0) = ss(i,j,0) - betax(i,j,3)
             ss(i,j,1) = 0.d0
             ss(i,j,XBC) = 0.d0
          end if
       end if
    end do

    ! y derivatives
    do i = lo(1),hi(1)
       do j = lo(2)+1,hi(2)-1
          ss(i,j,0) = ss(i,j,0) - betay(i,j,3)
       end do
    end do

    do i = lo(1),hi(1)
       bclo = stencil_bc_type(mask( i,lo(2)),2,-1)
       bchi = stencil_bc_type(mask( i,hi(2)),2,+1)

       j = lo(2)
       if (bclo .eq. BC_INT) then
          ss(i,j,0) = ss(i,j,0) - betay(i,j,3)
       elseif (bclo .eq. BC_NEU) then
          ss(i,j,0)   = ss(i,j,0) - betay(i,j,3) - betay(i,j,2)
          ss(i,j,4)   = 0.d0
          ss(i,j,YBC) = 0.d0
       elseif (bclo .eq. BC_DIR) then
          ss(i,j,0) = ss(i,j,0) - betay(i,j,3) 
          ss(i,j,4) = 0.d0
          ss(i,j,YBC) = 0.d0
       end if

       if ( hi(2) > lo(2) ) then
          j = hi(2)
          if (bchi .eq. BC_INT) then
             ss(i,j,0) = ss(i,j,0) - betay(i,j,3)
          elseif (bchi .eq. BC_NEU) then
             ss(i,j,0) = ss(i,j,0) - betay(i,j,3) - betay(i,j+1,1)
             ss(i,j,3) = 0.d0
             ss(i,j,YBC) = 0.d0
          elseif (bchi .eq. BC_DIR) then
             ss(i,j,0) = ss(i,j,0) - betay(i,j,3) 
             ss(i,j,3) = 0.d0
             ss(i,j,YBC) = 0.d0
          end if
       end if
    end do

    do j = lo(2),hi(2)
       do i = lo(1),hi(1)
          ss(i,j,0) = ss(i,j,0) + alpha(i,j,0) 
       end do
    end do

  end subroutine s_simpleg_2d_cc

  subroutine s_simple_3d_cc(ss, alpha, ng_a, betax, betay, betaz, ng_b, dh, mask, lo, hi, xa, xb, order)


    integer           , intent(in   ) :: ng_a, ng_b, lo(:), hi(:), order
    integer           , intent(inout) :: mask(lo(1)  :,lo(2)  :,lo(3)  :)
    real (kind = dp_t), intent(  out) ::   ss(lo(1)  :,lo(2)  :,lo(3)  :,0:)
    real (kind = dp_t), intent(in   ) :: alpha(lo(1)-ng_a:,lo(2)-ng_a:,lo(3)-ng_a:)
    real (kind = dp_t), intent(in   ) :: betax(lo(1)-ng_b:,lo(2)-ng_b:,lo(3)-ng_b:)
    real (kind = dp_t), intent(in   ) :: betay(lo(1)-ng_b:,lo(2)-ng_b:,lo(3)-ng_b:)
    real (kind = dp_t), intent(in   ) :: betaz(lo(1)-ng_b:,lo(2)-ng_b:,lo(3)-ng_b:)
    real (kind = dp_t), intent(in   ) :: dh(:)
    real (kind = dp_t), intent(in   ) :: xa(:), xb(:)

    real (kind = dp_t) :: f1(3)
    integer            :: i, j, k, bclo, bchi, nx, ny, nz
    integer            :: lnx, lny, lnz, lorder
    integer, parameter :: XBC = 7, YBC = 8, ZBC = 9

    nx = hi(1)-lo(1)+1
    ny = hi(2)-lo(2)+1
    nz = hi(3)-lo(3)+1
    f1 = ONE/dh**2

    ss(:,:,:,0) = ZERO

    ss(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1) = -betax(lo(1)+1:hi(1)+1,lo(2)  :hi(2),  lo(3)  :hi(3)  )*f1(1)
    ss(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),2) = -betax(lo(1)  :hi(1),  lo(2)  :hi(2),  lo(3)  :hi(3)  )*f1(1)

    ss(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),3) = -betay(lo(1)  :hi(1),  lo(2)+1:hi(2)+1,lo(3)  :hi(3)  )*f1(2)
    ss(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),4) = -betay(lo(1)  :hi(1),  lo(2)  :hi(2),  lo(3)  :hi(3)  )*f1(2)

    ss(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),5) = -betaz(lo(1)  :hi(1),  lo(2)  :hi(2),  lo(3)+1:hi(3)+1)*f1(3)
    ss(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),6) = -betaz(lo(1)  :hi(1),  lo(2)  :hi(2),  lo(3)  :hi(3)  )*f1(3)

    ss(:,:,:,XBC) = ZERO
    ss(:,:,:,YBC) = ZERO
    ss(:,:,:,ZBC) = ZERO

    mask = ibclr(mask, BC_BIT(BC_GEOM,1,-1))
    mask = ibclr(mask, BC_BIT(BC_GEOM,1,+1))
    mask = ibclr(mask, BC_BIT(BC_GEOM,2,-1))
    mask = ibclr(mask, BC_BIT(BC_GEOM,2,+1))
    mask = ibclr(mask, BC_BIT(BC_GEOM,3,-1))
    mask = ibclr(mask, BC_BIT(BC_GEOM,3,+1))

    lnx = nx; lny = ny; lnz = nz; lorder = order

    ! x derivatives

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1)+1,hi(1)-1
             ss(i,j,k,0) = ss(i,j,k,0) + (betax(i,j,k)+betax(i+1,j,k))*f1(1)
          end do
       end do
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i,j,k,bclo,bchi) FIRSTPRIVATE(lorder,lnx)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          bclo = stencil_bc_type(mask(lo(1),j,k),1,-1)
          bchi = stencil_bc_type(mask(hi(1),j,k),1,+1)

          i = lo(1)
          if (bclo .eq. BC_INT) then
             ss(i,j,k,0) = ss(i,j,k,0) + (betax(i,j,k)+betax(i+1,j,k))*f1(1)
          else
             call stencil_bndry_aaa(lorder, lnx, 1, -1, mask(i,j,k), &
                  ss(i,j,k,0), ss(i,j,k,1), ss(i,j,k,2), ss(i,j,k,XBC), &
                  betax(i,j,k), betax(i+1,j,k), xa(1), xb(1), dh(1), bclo, bchi)
          end if

          if ( hi(1) > lo(1) ) then
             i = hi(1)
             if (bchi .eq. BC_INT) then
                ss(i,j,k,0) = ss(i,j,k,0) + (betax(i,j,k)+betax(i+1,j,k))*f1(1)
             else
                call stencil_bndry_aaa(lorder, lnx, 1, 1, mask(i,j,k), &
                     ss(i,j,k,0), ss(i,j,k,1), ss(i,j,k,2), ss(i,j,k,XBC), &
                     betax(i,j,k), betax(i+1,j,k), xa(1), xb(1), dh(1), bclo, bchi)
             end if
          end if
       end do
    end do
    !$OMP END PARALLEL DO

    ! y derivatives

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3),hi(3)
       do j = lo(2)+1,hi(2)-1
          do i = lo(1),hi(1)
             ss(i,j,k,0) = ss(i,j,k,0) + (betay(i,j,k)+betay(i,j+1,k))*f1(2)
          end do
       end do
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i,j,k,bclo,bchi) FIRSTPRIVATE(lorder,lny)
    do k = lo(3),hi(3)
       do i = lo(1),hi(1)
          bclo = stencil_bc_type(mask(i,lo(2),k),2,-1)
          bchi = stencil_bc_type(mask(i,hi(2),k),2,+1)

          j = lo(2)
          if (bclo .eq. BC_INT) then
             ss(i,j,k,0) = ss(i,j,k,0) + (betay(i,j,k)+betay(i,j+1,k))*f1(2)
          else
             call stencil_bndry_aaa(lorder, lny, 2, -1, mask(i,j,k), &
                  ss(i,j,k,0), ss(i,j,k,3), ss(i,j,k,4),ss(i,j,k,YBC), &
                  betay(i,j,k), betay(i,j+1,k), xa(2), xb(2), dh(2), bclo, bchi)
          end if
          if ( hi(2) > lo(2) ) then
             j = hi(2)
             if (bchi .eq. BC_INT) then
                ss(i,j,k,0) = ss(i,j,k,0) + (betay(i,j,k)+betay(i,j+1,k))*f1(2)
             else
                call stencil_bndry_aaa(lorder, lny, 2, 1, mask(i,j,k), &
                     ss(i,j,k,0), ss(i,j,k,3), ss(i,j,k,4), ss(i,j,k,YBC), &
                     betay(i,j,k), betay(i,j+1,k), xa(2), xb(2), dh(2), bclo, bchi)
             end if
          end if
       end do
    end do
    !$OMP END PARALLEL DO

    ! z derivatives

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3)+1,hi(3)-1
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             ss(i,j,k,0) = ss(i,j,k,0) + (betaz(i,j,k)+betaz(i,j,k+1))*f1(3)
          end do
       end do
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i,j,k,bclo,bchi) FIRSTPRIVATE(lorder,lnz)
    do j = lo(2),hi(2)
       do i = lo(1),hi(1)
          bclo = stencil_bc_type(mask(i,j,lo(3)),3,-1)
          bchi = stencil_bc_type(mask(i,j,hi(3)),3,+1)

          k = lo(3)
          if (bclo .eq. BC_INT) then
             ss(i,j,k,0) = ss(i,j,k,0) + (betaz(i,j,k)+betaz(i,j,k+1))*f1(3)
          else
             call stencil_bndry_aaa(lorder, lnz, 3, -1, mask(i,j,k), &
                  ss(i,j,k,0), ss(i,j,k,5), ss(i,j,k,6),ss(i,j,k,ZBC), &
                  betaz(i,j,k), betaz(i,j,k+1), xa(3), xb(3), dh(3), bclo, bchi)
          end if
          if ( hi(3) > lo(3) ) then
             k = hi(3)
             if (bchi .eq. BC_INT) then
                ss(i,j,k,0) = ss(i,j,k,0) + (betaz(i,j,k)+betaz(i,j,k+1))*f1(3)
             else
                call stencil_bndry_aaa(lorder, lnz, 3, 1, mask(i,j,k), &
                     ss(i,j,k,0), ss(i,j,k,5), ss(i,j,k,6), ss(i,j,k,ZBC), &
                     betaz(i,j,k), betaz(i,j,k+1), xa(3), xb(3), dh(3), bclo, bchi)
             end if
          end if
       end do
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             ss(i,j,k,0) = ss(i,j,k,0) + alpha(i,j,k)
          end do
       end do
    end do
    !$OMP END PARALLEL DO

  end subroutine s_simple_3d_cc

  subroutine s_simpleg_3d_cc(ss, alpha, ng_a, betax, betay, betaz, ng_b, dh, mask, lo, hi)


    integer           , intent(in   ) :: ng_a, ng_b, lo(:), hi(:)
    integer           , intent(inout) :: mask(lo(1)  :,lo(2)  :,lo(3)  :)
    real (kind = dp_t), intent(  out) ::   ss(lo(1)  :,lo(2)  :,lo(3)  :,0:)
    real (kind = dp_t), intent(in   ) :: alpha(lo(1)-ng_a:,lo(2)-ng_a:,lo(3)-ng_a:,0:)
    real (kind = dp_t), intent(in   ) :: betax(lo(1)-ng_b:,lo(2)-ng_b:,lo(3)-ng_b:,:)
    real (kind = dp_t), intent(in   ) :: betay(lo(1)-ng_b:,lo(2)-ng_b:,lo(3)-ng_b:,:)
    real (kind = dp_t), intent(in   ) :: betaz(lo(1)-ng_b:,lo(2)-ng_b:,lo(3)-ng_b:,:)
    real (kind = dp_t), intent(in   ) :: dh(:)

    real (kind = dp_t) :: f1(3)
    integer            :: i, j, k, bclo, bchi, nx, ny, nz
    integer, parameter :: XBC = 7, YBC = 8, ZBC = 9

    nx = hi(1)-lo(1)+1
    ny = hi(2)-lo(2)+1
    nz = hi(3)-lo(3)+1
    f1 = ONE/dh**2

    ss(:,:,:,:) = ZERO

    ss(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1) = -betax(lo(1)+1:hi(1)+1,lo(2)  :hi(2),  lo(3)  :hi(3)  ,1)
    ss(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),2) = -betax(lo(1)  :hi(1),  lo(2)  :hi(2),  lo(3)  :hi(3)  ,2)

    ss(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),3) = -betay(lo(1)  :hi(1),  lo(2)+1:hi(2)+1,lo(3)  :hi(3)  ,1)
    ss(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),4) = -betay(lo(1)  :hi(1),  lo(2)  :hi(2),  lo(3)  :hi(3)  ,2)

    ss(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),5) = -betaz(lo(1)  :hi(1),  lo(2)  :hi(2),  lo(3)+1:hi(3)+1,1)
    ss(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),6) = -betaz(lo(1)  :hi(1),  lo(2)  :hi(2),  lo(3)  :hi(3)  ,2)

    ss(:,:,:,XBC) = ZERO
    ss(:,:,:,YBC) = ZERO
    ss(:,:,:,ZBC) = ZERO

    mask = ibclr(mask, BC_BIT(BC_GEOM,1,-1))
    mask = ibclr(mask, BC_BIT(BC_GEOM,1,+1))
    mask = ibclr(mask, BC_BIT(BC_GEOM,2,-1))
    mask = ibclr(mask, BC_BIT(BC_GEOM,2,+1))
    mask = ibclr(mask, BC_BIT(BC_GEOM,3,-1))
    mask = ibclr(mask, BC_BIT(BC_GEOM,3,+1))

    ! x derivatives

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1)+1,hi(1)-1
             ss(i,j,k,0) = ss(i,j,k,0) - betax(i,j,k,3)
          end do
       end do
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i,j,k,bclo,bchi) 
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          bclo = stencil_bc_type(mask(lo(1),j,k),1,-1)
          bchi = stencil_bc_type(mask(hi(1),j,k),1,+1)

          i = lo(1)
          if (bclo .eq. BC_INT) then
             ss(i,j,k,0) = ss(i,j,k,0) - betax(i,j,k,3)
          elseif (bclo .eq. BC_NEU) then
             ss(i,j,k,0) = ss(i,j,k,0) - betax(i,j,k,3) - betax(i,j,k,2)
             ss(i,j,k,2) = 0.d0
             ss(i,j,k,XBC) = 0.d0
          elseif (bclo .eq. BC_DIR) then
             ss(i,j,k,0) = ss(i,j,k,0) - betax(i,j,k,3)
             ss(i,j,k,2) = 0.d0
             ss(i,j,k,XBC) = 0.d0
          end if

          if ( hi(1) > lo(1) ) then
             i = hi(1)
             if (bchi .eq. BC_INT) then
                ss(i,j,k,0) = ss(i,j,k,0) - betax(i,j,k,3)
             elseif (bchi .eq. BC_NEU) then
                ss(i,j,k,0) = ss(i,j,k,0) - betax(i,j,k,3) - betax(i+1,j,k,1)
                ss(i,j,k,1) = 0.d0
                ss(i,j,k,XBC) = 0.d0
             elseif (bchi .eq. BC_DIR) then
                ss(i,j,k,0) = ss(i,j,k,0) - betax(i,j,k,3)
                ss(i,j,k,1) = 0.d0
                ss(i,j,k,XBC) = 0.d0
             end if
          end if
       end do
    end do
    !$OMP END PARALLEL DO

    ! y derivatives

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3),hi(3)
       do j = lo(2)+1,hi(2)-1
          do i = lo(1),hi(1)
             ss(i,j,k,0) = ss(i,j,k,0) - betay(i,j,k,3)
          end do
       end do
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i,j,k,bclo,bchi) 
    do k = lo(3),hi(3)
       do i = lo(1),hi(1)
          bclo = stencil_bc_type(mask(i,lo(2),k),2,-1)
          bchi = stencil_bc_type(mask(i,hi(2),k),2,+1)

          j = lo(2)
          if (bclo .eq. BC_INT) then
             ss(i,j,k,0) = ss(i,j,k,0) - betay(i,j,k,3)
          elseif (bclo .eq. BC_NEU) then
             ss(i,j,k,0)   = ss(i,j,k,0) - betay(i,j,k,3) - betay(i,j,k,2)
             ss(i,j,k,4)   = 0.d0
             ss(i,j,k,YBC) = 0.d0
          elseif (bclo .eq. BC_DIR) then
             ss(i,j,k,0)   = ss(i,j,k,0) - betay(i,j,k,3) 
             ss(i,j,k,4)   = 0.d0
             ss(i,j,k,YBC) = 0.d0
          end if

          if ( hi(2) > lo(2) ) then
             j = hi(2)
             if (bchi .eq. BC_INT) then
                ss(i,j,k,0) = ss(i,j,k,0) - betay(i,j,k,3)
             elseif (bchi .eq. BC_NEU) then
                ss(i,j,k,0)   = ss(i,j,k,0) - betay(i,j,k,3) - betay(i,j+1,k,1)
                ss(i,j,k,3)   = 0.d0
                ss(i,j,k,YBC) = 0.d0
             elseif (bchi .eq. BC_DIR) then
                ss(i,j,k,0)   = ss(i,j,k,0) - betay(i,j,k,3) 
                ss(i,j,k,3)   = 0.d0
                ss(i,j,k,YBC) = 0.d0
             end if
          end if
       end do
    end do
    !$OMP END PARALLEL DO

    ! z derivatives

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3)+1,hi(3)-1
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             ss(i,j,k,0) = ss(i,j,k,0) - betaz(i,j,k,3)
          end do
       end do
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i,j,k,bclo,bchi) 
    do j = lo(2),hi(2)
       do i = lo(1),hi(1)
          bclo = stencil_bc_type(mask(i,j,lo(3)),3,-1)
          bchi = stencil_bc_type(mask(i,j,hi(3)),3,+1)

          k = lo(3)
          if (bclo .eq. BC_INT) then
             ss(i,j,k,0) = ss(i,j,k,0) - betaz(i,j,k,3)
          elseif (bclo .eq. BC_NEU) then
             ss(i,j,k,0)   = ss(i,j,k,0) - betaz(i,j,k,3) - betaz(i,j,k,2)
             ss(i,j,k,6)   = 0.d0
             ss(i,j,k,ZBC) = 0.d0
          elseif (bclo .eq. BC_DIR) then
             ss(i,j,k,0)   = ss(i,j,k,0) - betaz(i,j,k,3) 
             ss(i,j,k,6)   = 0.d0
             ss(i,j,k,ZBC) = 0.d0
          end if

          if ( hi(3) > lo(3) ) then
             k = hi(3)
             if (bchi .eq. BC_INT) then
                ss(i,j,k,0) = ss(i,j,k,0) - betaz(i,j,k,3)
             elseif (bchi .eq. BC_NEU) then
                ss(i,j,k,0)   = ss(i,j,k,0) - betaz(i,j,k,3) - betaz(i,j,k+1,1)
                ss(i,j,k,5)   = 0.d0
                ss(i,j,k,ZBC) = 0.d0
             elseif (bchi .eq. BC_DIR) then
                ss(i,j,k,0)   = ss(i,j,k,0) - betaz(i,j,k,3) 
                ss(i,j,k,5)   = 0.d0
                ss(i,j,k,ZBC) = 0.d0
             end if
          end if
       end do
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             ss(i,j,k,0) = ss(i,j,k,0) + alpha(i,j,k,0)
          end do
       end do
    end do
    !$OMP END PARALLEL DO

  end subroutine s_simpleg_3d_cc

  subroutine s_minion_cross_fill_2d(ss, alpha, ng_a, betax, betay, ng_b, dh, mask, lo, hi)

    integer           , intent(in   ) :: ng_a,ng_b
    integer           , intent(in   ) :: lo(:), hi(:)
    integer           , intent(inout) :: mask(lo(1):,lo(2):)
    real (kind = dp_t), intent(  out) :: ss(lo(1):,lo(2):,0:)
    real (kind = dp_t), intent(in   ) :: alpha(lo(1)-ng_a:,lo(2)-ng_a:)
    real (kind = dp_t), intent(in   ) :: betax(lo(1)-ng_b:,lo(2)-ng_b:)
    real (kind = dp_t), intent(in   ) :: betay(lo(1)-ng_b:,lo(2)-ng_b:)
    real (kind = dp_t), intent(in   ) :: dh(:)
    integer nx, ny
    integer i, j

    nx = hi(1)-lo(1)+1
    ny = hi(2)-lo(2)+1

    mask = ibclr(mask, BC_BIT(BC_GEOM,1,-1))
    mask = ibclr(mask, BC_BIT(BC_GEOM,1,+1))
    mask = ibclr(mask, BC_BIT(BC_GEOM,2,-1))
    mask = ibclr(mask, BC_BIT(BC_GEOM,2,+1))

    ss = 0.0d0

    ! We only include the beta's here to get the viscous coefficients in here for now.
    ! The projection has beta == 1.
    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          ss(i,j,1) =   1.d0 * betax(i  ,j)
          ss(i,j,2) = -16.d0 * betax(i  ,j)
          ss(i,j,3) = -16.d0 * betax(i+1,j)
          ss(i,j,4) =   1.d0 * betax(i+1,j)
          ss(i,j,5) =   1.d0 * betay(i,j  )
          ss(i,j,6) = -16.d0 * betay(i,j  )
          ss(i,j,7) = -16.d0 * betay(i,j+1)
          ss(i,j,8) =   1.d0 * betay(i,j+1)
          ss(i,j,0) = -(ss(i,j,1) + ss(i,j,2) + ss(i,j,3) + ss(i,j,4) &
                       +ss(i,j,5) + ss(i,j,6) + ss(i,j,7) + ss(i,j,8) )
       end do
    end do

    ss = ss * (ONE / (12.d0 * dh(1)**2))

    ! This adds the "alpha" term in (alpha - del dot beta grad) phi = RHS.
    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          ss(i,j,0) = ss(i,j,0) + alpha(i,j)
       end do
    end do

  end subroutine s_minion_cross_fill_2d

  subroutine s_minion_full_old_2d(ss, beta, ng_b, dh, mask, lo, hi)

    integer           , intent(in   ) :: ng_b
    integer           , intent(in   ) :: lo(:), hi(:)
    integer           , intent(inout) :: mask(:,:)
    real (kind = dp_t), intent(  out) :: ss(:,:,0:)
    real (kind = dp_t), intent(inout) :: beta(1-ng_b:,1-ng_b:,0:)
    real (kind = dp_t), intent(in   ) :: dh(:)

    integer            :: i, j, nx, ny
    real (kind = dp_t) :: fac

    nx = hi(1)-lo(1)+1
    ny = hi(2)-lo(2)+1

    mask = ibclr(mask, BC_BIT(BC_GEOM,1,-1))
    mask = ibclr(mask, BC_BIT(BC_GEOM,1,+1))
    mask = ibclr(mask, BC_BIT(BC_GEOM,2,-1))
    mask = ibclr(mask, BC_BIT(BC_GEOM,2,+1))

    ss = 0.0d0

    ! First use the betax coefficients
    do j = 1, ny
       do i = 1, nx

          ss(i,j,12) = 27648.d0*beta(i+1,j,1) + 414720.d0 * beta(i  ,j,1)
          ss(i,j,13) = 27648.d0*beta(i  ,j,1) + 414720.d0 * beta(i+1,j,1)

          ss(i,j,11) = -27648.d0 * beta(i  ,j,1)
          ss(i,j,14) = -27648.d0 * beta(i+1,j,1)

          ss(i,j, 8) = -2550.d0 * beta(i,j+2,1)  - 2550.d0 * beta(i+1,j+2,1) &
                      +17340.d0 * beta(i,j+1,1) + 17340.d0 * beta(i+1,j+1,1) &
                      -17340.d0 * beta(i,j-1,1) - 17340.d0 * beta(i+1,j-1,1) &
                      + 2550.d0 * beta(i,j-2,1)  + 2550.d0 * beta(i+1,j-2,1)

          ss(i,j,17) = -2550.d0 * beta(i,j-2,1)  - 2550.d0 * beta(i+1,j-2,1) &
                      +17340.d0 * beta(i,j-1,1) + 17340.d0 * beta(i+1,j-1,1) &
                      -17340.d0 * beta(i,j+1,1) - 17340.d0 * beta(i+1,j+1,1) &
                      + 2550.d0 * beta(i,j+2,1)  + 2550.d0 * beta(i+1,j+2,1)

          ss(i,j, 7) =  170.d0 * beta(i+1,j+2,1) +  2550.d0 * beta(i,j+2,1) &
                      -1156.d0 * beta(i+1,j+1,1) - 17340.d0 * beta(i,j+1,1) &
                      +1156.d0 * beta(i+1,j-1,1) + 17340.d0 * beta(i,j-1,1) &
                      - 170.d0 * beta(i+1,j-2,1) -  2550.d0 * beta(i,j-2,1) 

          ss(i,j,16) =  170.d0 * beta(i+1,j-2,1) +  2550.d0 * beta(i,j-2,1) &
                      -1156.d0 * beta(i+1,j-1,1) - 17340.d0 * beta(i,j-1,1) &
                      +1156.d0 * beta(i+1,j+1,1) + 17340.d0 * beta(i,j+1,1) &
                      - 170.d0 * beta(i+1,j+2,1) -  2550.d0 * beta(i,j+2,1) 

          ss(i,j, 9) =  170.d0 * beta(i,j+2,1) +  2550.d0 * beta(i+1,j+2,1) &
                      -1156.d0 * beta(i,j+1,1) - 17340.d0 * beta(i+1,j+1,1) &
                      +1156.d0 * beta(i,j-1,1) + 17340.d0 * beta(i+1,j-1,1) &
                      - 170.d0 * beta(i,j-2,1) -  2550.d0 * beta(i+1,j-2,1) 

          ss(i,j,18) =  170.d0 * beta(i,j-2,1) +  2550.d0 * beta(i+1,j-2,1) &
                      -1156.d0 * beta(i,j-1,1) - 17340.d0 * beta(i+1,j-1,1) &
                      +1156.d0 * beta(i,j+1,1) + 17340.d0 * beta(i+1,j+1,1) &
                      - 170.d0 * beta(i,j+2,1) -  2550.d0 * beta(i+1,j+2,1) 

          ss(i,j, 6) = -170.d0 * beta(i,j+2,1) +  1156.d0 * beta(i,j+1,1) &
                       +170.d0 * beta(i,j-2,1) -  1156.d0 * beta(i,j-1,1)
          ss(i,j,15) = -170.d0 * beta(i,j-2,1) +  1156.d0 * beta(i,j-1,1) &
                       +170.d0 * beta(i,j+2,1) -  1156.d0 * beta(i,j+1,1)
          ss(i,j,10) = -170.d0 * beta(i+1,j+2,1) +  1156.d0 * beta(i+1,j+1,1) &
                       +170.d0 * beta(i+1,j-2,1) -  1156.d0 * beta(i+1,j-1,1)
          ss(i,j,19) = -170.d0 * beta(i+1,j-2,1) +  1156.d0 * beta(i+1,j-1,1) &
                       +170.d0 * beta(i+1,j+2,1) -  1156.d0 * beta(i+1,j+1,1)

          ss(i,j, 3) =   375.d0 * beta(i,j+2,1) +  375.d0 * beta(i+1,j+2,1) &
                      - 2550.d0 * beta(i,j+1,1) - 2550.d0 * beta(i+1,j+1,1) &
                      + 2550.d0 * beta(i,j-1,1) + 2550.d0 * beta(i+1,j-1,1) &
                      -  375.d0 * beta(i,j-2,1) -  375.d0 * beta(i+1,j-2,1)
          ss(i,j,22) =  375.d0 * beta(i,j-2,1) +  375.d0 * beta(i+1,j-2,1) &
                      -2550.d0 * beta(i,j-1,1) - 2550.d0 * beta(i+1,j-1,1) &
                      +2550.d0 * beta(i,j+1,1) + 2550.d0 * beta(i+1,j+1,1) &
                      - 375.d0 * beta(i,j+2,1) -  375.d0 * beta(i+1,j+2,1)

          ss(i,j, 2) = - 25.d0 * beta(i+1,j+2,1) -  375.d0 * beta(i,j+2,1) &
                       +170.d0 * beta(i+1,j+1,1) + 2550.d0 * beta(i,j+1,1) &
                       -170.d0 * beta(i+1,j-1,1) - 2550.d0 * beta(i,j-1,1) &
                       + 25.d0 * beta(i+1,j-2,1) +  375.d0 * beta(i,j-2,1)
          ss(i,j,21) = - 25.d0 * beta(i+1,j-2,1) -  375.d0 * beta(i,j-2,1) &
                       +170.d0 * beta(i+1,j-1,1) + 2550.d0 * beta(i,j-1,1) &
                       -170.d0 * beta(i+1,j+1,1) - 2550.d0 * beta(i,j+1,1) &
                       + 25.d0 * beta(i+1,j+2,1) +  375.d0 * beta(i,j+2,1)
          ss(i,j, 4) = - 25.d0 * beta(i,j+2,1) -  375.d0 * beta(i+1,j+2,1) &
                       +170.d0 * beta(i,j+1,1) + 2550.d0 * beta(i+1,j+1,1) &
                       -170.d0 * beta(i,j-1,1) - 2550.d0 * beta(i+1,j-1,1) &
                       + 25.d0 * beta(i,j-2,1) +  375.d0 * beta(i+1,j-2,1)
          ss(i,j,23) = - 25.d0 * beta(i,j-2,1) -  375.d0 * beta(i+1,j-2,1) &
                       +170.d0 * beta(i,j-1,1) + 2550.d0 * beta(i+1,j-1,1) &
                       -170.d0 * beta(i,j+1,1) - 2550.d0 * beta(i+1,j+1,1) &
                       + 25.d0 * beta(i,j+2,1) +  375.d0 * beta(i+1,j+2,1)

          ss(i,j, 1) =   25.d0 * beta(i,j+2,1) -  170.d0 * beta(i,j+1,1) &
                        -25.d0 * beta(i,j-2,1) +  170.d0 * beta(i,j-1,1)
          ss(i,j, 5) =   25.d0 * beta(i+1,j+2,1) -  170.d0 * beta(i+1,j+1,1) &
                        -25.d0 * beta(i+1,j-2,1) +  170.d0 * beta(i+1,j-1,1)
          ss(i,j,20) =   25.d0 * beta(i,j-2,1) -  170.d0 * beta(i,j-1,1) &
                        -25.d0 * beta(i,j+2,1) +  170.d0 * beta(i,j+1,1)
          ss(i,j,24) =   25.d0 * beta(i+1,j-2,1) -  170.d0 * beta(i+1,j-1,1) &
                        -25.d0 * beta(i+1,j+2,1) +  170.d0 * beta(i+1,j+1,1)

          ss(i,j, 0) = -414720.d0 * (beta(i,j,1) + beta(i+1,j,1))

       end do
    end do

    ! Then use the betay coefficients
    do j = 1, ny
       do i = 1, nx

          ss(i,j, 8) = ss(i,j, 8) + 27648.d0*beta(i,j+1,2) + 414720.d0 * beta(i,j  ,2)
          ss(i,j,17) = ss(i,j,17) + 27648.d0*beta(i,j  ,2) + 414720.d0 * beta(i,j+1,2)

          ss(i,j, 3) = ss(i,j, 3) - 27648.d0 * beta(i,j  ,2)
          ss(i,j,22) = ss(i,j,22) - 27648.d0 * beta(i,j+1,2)

          ss(i,j,12) = ss(i,j,12) & 
                       -2550.d0 * beta(i+2,j,2)  - 2550.d0 * beta(i+2,j+1,2) &
                      +17340.d0 * beta(i+1,j,2) + 17340.d0 * beta(i+1,j+1,2) &
                      -17340.d0 * beta(i-1,j,2) - 17340.d0 * beta(i-1,j+1,2) &
                      + 2550.d0 * beta(i-2,j,2)  + 2550.d0 * beta(i-2,j+1,2)

          ss(i,j,13) = ss(i,j,13) & 
                       -2550.d0 * beta(i-2,j,2)  - 2550.d0 * beta(i-2,j+1,2) &
                      +17340.d0 * beta(i-1,j,2) + 17340.d0 * beta(i-1,j+1,2) &
                      -17340.d0 * beta(i+1,j,2) - 17340.d0 * beta(i+1,j+1,2) &
                      + 2550.d0 * beta(i+2,j,2)  + 2550.d0 * beta(i+2,j+1,2)

          ss(i,j, 7) = ss(i,j, 7) &
                      + 170.d0 * beta(i+2,j+1,2) +  2550.d0 * beta(i+2,j  ,2) &
                      -1156.d0 * beta(i+1,j+1,2) - 17340.d0 * beta(i+1,j  ,2) &
                      +1156.d0 * beta(i-1,j+1,2) + 17340.d0 * beta(i-1,j  ,2) &
                      - 170.d0 * beta(i-2,j+1,2) -  2550.d0 * beta(i-2,j  ,2) 

          ss(i,j,16) = ss(i,j,16) &  
                      + 170.d0 * beta(i+2,j  ,2) +  2550.d0 * beta(i+2,j+1,2) &
                      -1156.d0 * beta(i+1,j  ,2) - 17340.d0 * beta(i+1,j+1,2) &
                      +1156.d0 * beta(i-1,j  ,2) + 17340.d0 * beta(i-1,j+1,2) &
                      - 170.d0 * beta(i-2,j  ,2) -  2550.d0 * beta(i-2,j+1,2) 

          ss(i,j, 9) = ss(i,j, 9) &  
                     +  170.d0 * beta(i-2,j+1,2) +  2550.d0 * beta(i-2,j  ,2) &
                      -1156.d0 * beta(i-1,j+1,2) - 17340.d0 * beta(i-1,j  ,2) &
                      +1156.d0 * beta(i+1,j+1,2) + 17340.d0 * beta(i+1,j  ,2) &
                      - 170.d0 * beta(i+2,j+1,2) -  2550.d0 * beta(i+2,j  ,2) 

          ss(i,j,18) = ss(i,j,18) &  
                     +  170.d0 * beta(i-2,j  ,2) +  2550.d0 * beta(i-2,j+1,2) &
                      -1156.d0 * beta(i-1,j  ,2) - 17340.d0 * beta(i-1,j+1,2) &
                      +1156.d0 * beta(i+1,j  ,2) + 17340.d0 * beta(i+1,j+1,2) &
                      - 170.d0 * beta(i+2,j  ,2) -  2550.d0 * beta(i+2,j+1,2) 

          ss(i,j, 2) = ss(i,j, 2) &
                       -170.d0 * beta(i+2,j,2) +  1156.d0 * beta(i+1,j,2) &
                       +170.d0 * beta(i-2,j,2) -  1156.d0 * beta(i-1,j,2)

          ss(i,j,21) = ss(i,j,21) &
                       -170.d0 * beta(i+2,j+1,2) +  1156.d0 * beta(i+1,j+1,2) &
                       +170.d0 * beta(i-2,j+1,2) -  1156.d0 * beta(i-1,j+1,2)

          ss(i,j, 4) = ss(i,j, 4) &
                       -170.d0 * beta(i-2,j,2) +  1156.d0 * beta(i-1,j,2) &
                       +170.d0 * beta(i+2,j,2) -  1156.d0 * beta(i+1,j,2)

          ss(i,j,23) = ss(i,j,23) &
                       -170.d0 * beta(i-2,j+1,2) +  1156.d0 * beta(i-1,j+1,2) &
                       +170.d0 * beta(i+2,j+1,2) -  1156.d0 * beta(i+1,j+1,2)

          ss(i,j,11) = ss(i,j,11) &
                      +  375.d0 * beta(i+2,j,2) +  375.d0 * beta(i+2,j+1,2) &
                      - 2550.d0 * beta(i+1,j,2) - 2550.d0 * beta(i+1,j+1,2) &
                      + 2550.d0 * beta(i-1,j,2) + 2550.d0 * beta(i-1,j+1,2) &
                      -  375.d0 * beta(i-2,j,2) -  375.d0 * beta(i-2,j+1,2)

          ss(i,j,14) = ss(i,j,14) &
                     +  375.d0 * beta(i-2,j,2) +  375.d0 * beta(i-2,j+1,2) &
                      -2550.d0 * beta(i-1,j,2) - 2550.d0 * beta(i-1,j+1,2) &
                      +2550.d0 * beta(i+1,j,2) + 2550.d0 * beta(i+1,j+1,2) &
                      - 375.d0 * beta(i+2,j,2) -  375.d0 * beta(i+2,j+1,2)

          ss(i,j, 6) = ss(i,j, 6) &
                       - 25.d0 * beta(i+2,j+1,2) -  375.d0 * beta(i+2,j,2) &
                       +170.d0 * beta(i+1,j+1,2) + 2550.d0 * beta(i+1,j,2) &
                       -170.d0 * beta(i-1,j+1,2) - 2550.d0 * beta(i-1,j,2) &
                       + 25.d0 * beta(i-2,j+1,2) +  375.d0 * beta(i-2,j,2)

          ss(i,j,15) = ss(i,j,15) &
                       - 25.d0 * beta(i+2,j,2) -  375.d0 * beta(i+2,j+1,2) &
                       +170.d0 * beta(i+1,j,2) + 2550.d0 * beta(i+1,j+1,2) &
                       -170.d0 * beta(i-1,j,2) - 2550.d0 * beta(i-1,j+1,2) &
                       + 25.d0 * beta(i-2,j,2) +  375.d0 * beta(i-2,j+1,2)

          ss(i,j,10) = ss(i,j,10) &
                       - 25.d0 * beta(i-2,j+1,2) -  375.d0 * beta(i-2,j,2) &
                       +170.d0 * beta(i-1,j+1,2) + 2550.d0 * beta(i-1,j,2) &
                       -170.d0 * beta(i+1,j+1,2) - 2550.d0 * beta(i+1,j,2) &
                       + 25.d0 * beta(i+2,j+1,2) +  375.d0 * beta(i+2,j,2)

          ss(i,j,19) = ss(i,j,19) &
                       - 25.d0 * beta(i-2,j,2) -  375.d0 * beta(i-2,j+1,2) &
                       +170.d0 * beta(i-1,j,2) + 2550.d0 * beta(i-1,j+1,2) &
                       -170.d0 * beta(i+1,j,2) - 2550.d0 * beta(i+1,j+1,2) &
                       + 25.d0 * beta(i+2,j,2) +  375.d0 * beta(i+2,j+1,2)

          ss(i,j, 1) = ss(i,j, 1) &
                       + 25.d0 * beta(i+2,j,2) -  170.d0 * beta(i+1,j,2) &
                        -25.d0 * beta(i-2,j,2) +  170.d0 * beta(i-1,j,2)
          ss(i,j, 5) = ss(i,j, 5) &
                       + 25.d0 * beta(i-2,j,2) -  170.d0 * beta(i-1,j,2) &
                        -25.d0 * beta(i+2,j,2) +  170.d0 * beta(i+1,j,2)
          ss(i,j,20) = ss(i,j,20) &
                       + 25.d0 * beta(i+2,j+1,2) -  170.d0 * beta(i+1,j+1,2) &
                        -25.d0 * beta(i-2,j+1,2) +  170.d0 * beta(i-1,j+1,2)
          ss(i,j,24) = ss(i,j,24) &
                       + 25.d0 * beta(i-2,j+1,2) -  170.d0 * beta(i-1,j+1,2) &
                        -25.d0 * beta(i+2,j+1,2) +  170.d0 * beta(i+1,j+1,2)

          ss(i,j, 0) = ss(i,j,0) -414720.d0 * ( beta(i,j,2) + beta(i,j+1,2) )

       end do
    end do
  
    fac = -1.d0 / (12.d0**2 * 48.d0**2 * dh(1)**2)

    ss = fac * ss

    ! This adds the "alpha" term in (alpha - del dot beta grad) phi = RHS.
    do j = 1, ny
       do i = 1, nx
          ss(i,j,0) = ss(i,j,0) + beta(i,j,0)
       end do
    end do

  end subroutine s_minion_full_old_2d

  subroutine s_minion_full_fill_2d(ss, alpha, ng_a, betax, betay, ng_b, dh, mask, lo, hi)

    integer           , intent(in   ) :: ng_a,ng_b
    integer           , intent(in   ) :: lo(:), hi(:)
    integer           , intent(inout) :: mask(lo(1):,lo(2):)
    real (kind = dp_t), intent(  out) ::    ss(lo(1):     ,lo(2):     ,0:)
    real (kind = dp_t), intent(in   ) :: alpha(lo(1)-ng_a:,lo(2)-ng_a:)
    real (kind = dp_t), intent(in   ) :: betax(lo(1)-ng_b:,lo(2)-ng_b:)
    real (kind = dp_t), intent(in   ) :: betay(lo(1)-ng_b:,lo(2)-ng_b:)
    real (kind = dp_t), intent(in   ) :: dh(:)

    integer          :: i, j, nsten
    double precision :: hx2,hy2,hx22,hy22,ss_sum
    double precision :: rholy,rhory,rhotx,rhobx  !  Transverse derivatives
    double precision :: s1,s2,scale              ! Coefficients for slopes

    mask = ibclr(mask, BC_BIT(BC_GEOM,1,-1))
    mask = ibclr(mask, BC_BIT(BC_GEOM,1,+1))
    mask = ibclr(mask, BC_BIT(BC_GEOM,2,-1))
    mask = ibclr(mask, BC_BIT(BC_GEOM,2,+1))

    !  Remember we are doing -del beta grad here (the alpha is done very last)
    !
    !    The stencil ordering for phi is
    !     20 21 22 23 24   
    !     15 16 17 18 19   
    !     11 12 0  13 14   
    !     6  7  8  9  10
    !     1  2  3  4  5
    !    The points for beta are at i,j,1 and i+1,j,1 for left and right
    !                           and i,j,2 and i,j+1,2 for top  and bottom
    !
    !  The stencil has two parts: the regular stencil and the correction.

    !  start with zero
    ss=0.0d0

    ! These are the coefficients in the second order part
    hx22 = -1.d0 / (12.d0 * dh(1)**2 )
    hy22 = -1.d0 / (12.d0 * dh(2)**2 )

    !  These are the coefficients for the tangential slopes
    !  We don't divide by h because the product of slopes is multiplied by h^2/12
    s1 = 34.d0 / 48.d0
    s2 = -5.d0 / 48.d0

    !  In the output of make_stencil2d.f90, the coefficents of the correction are multiplied by
    !  all the denominators so that they are integers.  So we have to put the denominators back in
    scale = 1.0d0/(12.0d0*48.d0)

    !  The coefficients hx2 and hy2 are defined by  (the minus sign is because it is minus beta*Lap)
    hx2 = -1.0d0/(12.d0*dh(1)**2)*scale
    hy2 = -1.0d0/(12.d0*dh(2)**2)*scale

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          ss_sum=0.0d0

          rholy = s1*(betax(i  ,j+1)- betax(i  ,j-1)) + s2*(betax(i  ,j+2)- betax(i  ,j-2))
          rhory = s1*(betax(i+1,j+1)- betax(i+1,j-1)) + s2*(betax(i+1,j+2)- betax(i+1,j-2))
          rhotx = s1*(betay(i+1,j+1)- betay(i-1,j+1)) + s2*(betay(i+2,j+1)- betay(i-2,j+1))
          rhobx = s1*(betay(i+1,j  )- betay(i-1,j  )) + s2*(betay(i+2,j  )- betay(i-2,j  ))

          !   DOING CONTRIB AT           -2          -2
         nsten =   1
           ss(i,j,nsten) = (   &
                             -5.0d0*rholy*hx2  &
                             -5.0d0*rhobx*hy2  )
                                               
            !   DOING CONTRIB AT           -1          -2

         nsten =   2
           ss(i,j,nsten) = (   &
                         +     5.0d0*rhory*hx2  &
                         +    75.0d0*rholy*hx2  &
                         +    34.0d0*rhobx*hy2  )
                                               
            !   DOING CONTRIB AT            0          -2
         nsten =   3
           ss(i,j,nsten) = (   &
                            -75.0d0*rhory*hx2  &
                            -75.0d0*rholy*hx2  )
                                               
            !   DOING CONTRIB AT            1          -2
         nsten =   4
           ss(i,j,nsten) = (   &
                         +    75.0d0*rhory*hx2  &
                         +     5.0d0*rholy*hx2  &
                            -34.0d0*rhobx*hy2   )
                                               
            !   DOING CONTRIB AT            2          -2
         nsten =   5
           ss(i,j,nsten) = (   &
                             -5.0d0*rhory*hx2  &
                         +     5.0d0*rhobx*hy2 )
                                               
            !   DOING CONTRIB AT           -2          -1
         nsten =   6
           ss(i,j,nsten) = (   &
                         +    34.0d0*rholy*hx2  &
                         +     5.0d0*rhotx*hy2  &
                         +    75.0d0*rhobx*hy2  )
                                               
            !   DOING CONTRIB AT           -1          -1
         nsten =   7
           ss(i,j,nsten) = (   &
                            -34.0d0*rhory*hx2  &
                           -510.0d0*rholy*hx2  &
                            -34.0d0*rhotx*hy2  &
                           -510.0d0*rhobx*hy2  )
                                               
            !   DOING CONTRIB AT            0          -1
         nsten =   8
           ss(i,j,nsten) = (   &
                         +   510.0d0*rhory*hx2  &
                         +   510.0d0*rholy*hx2  )
                                               
            !   DOING CONTRIB AT            1          -1
         nsten =   9
           ss(i,j,nsten) = (   &
                           -510.0d0*rhory*hx2  &
                            -34.0d0*rholy*hx2  &
                         +    34.0d0*rhotx*hy2  &
                         +   510.0d0*rhobx*hy2  )
                                               
            !   DOING CONTRIB AT            2          -1
         nsten =  10
           ss(i,j,nsten) = (   &
                         +    34.0d0*rhory*hx2  &
                             -5.0d0*rhotx*hy2  &
                            -75.0d0*rhobx*hy2  )
                                               
            !   DOING CONTRIB AT           -2           0
         nsten =  11
           ss(i,j,nsten) = (   &
                            -75.0d0*rhotx*hy2  &
                            -75.0d0*rhobx*hy2  )
                                               
            !   DOING CONTRIB AT           -1           0
         nsten =  12
           ss(i,j,nsten) = (   &
                         +   510.0d0*rhotx*hy2  &
                         +   510.0d0*rhobx*hy2  )
                                               
            !   DOING CONTRIB AT            1           0
         nsten =  13
           ss(i,j,nsten) = (   &
                           -510.0d0*rhotx*hy2  &
                           -510.0d0*rhobx*hy2  )
                                               
            !   DOING CONTRIB AT            2           0
         nsten =  14
           ss(i,j,nsten) = (   &
                         +    75.0d0*rhotx*hy2  &
                         +    75.0d0*rhobx*hy2  )
                                               
            !   DOING CONTRIB AT           -2           1
         nsten =  15
           ss(i,j,nsten) = (   &
                            -34.0d0*rholy*hx2  &
                         +    75.0d0*rhotx*hy2  &
                         +     5.0d0*rhobx*hy2  )
                                               
            !   DOING CONTRIB AT           -1           1
         nsten =  16
           ss(i,j,nsten) = (   &
                         +    34.0d0*rhory*hx2  &
                         +   510.0d0*rholy*hx2  &
                           -510.0d0*rhotx*hy2  &
                            -34.0d0*rhobx*hy2  )
                                               
            !   DOING CONTRIB AT            0           1
         nsten =  17
           ss(i,j,nsten) = (   &
                           -510.0d0*rhory*hx2  &
                           -510.0d0*rholy*hx2  )
                                               
            !   DOING CONTRIB AT            1           1
         nsten =  18
           ss(i,j,nsten) = (   &
                         +   510.0d0*rhory*hx2  &
                         +    34.0d0*rholy*hx2  &
                         +   510.0d0*rhotx*hy2  &
                         +    34.0d0*rhobx*hy2  )
                                               
            !   DOING CONTRIB AT            2           1
         nsten =  19
           ss(i,j,nsten) = (   &
                            -34.0d0*rhory*hx2  &
                            -75.0d0*rhotx*hy2  &
                             -5.0d0*rhobx*hy2  )
                                               
            !   DOING CONTRIB AT           -2           2
         nsten =  20
           ss(i,j,nsten) = (   &
                         +     5.0d0*rholy*hx2 &
                             -5.0d0*rhotx*hy2  )
                                              
            !   DOING CONTRIB AT           -1           2
         nsten =  21
           ss(i,j,nsten) = (   &
                             -5.0d0*rhory*hx2  &
                            -75.0d0*rholy*hx2  &
                         +    34.0d0*rhotx*hy2 )
                                              
            !   DOING CONTRIB AT            0           2
         nsten =  22
           ss(i,j,nsten) = (   &
                         +    75.0d0*rhory*hx2  &
                         +    75.0d0*rholy*hx2  )
                                              
            !   DOING CONTRIB AT            1           2
         nsten =  23
           ss(i,j,nsten) = (   &
                            -75.0d0*rhory*hx2  &
                             -5.0d0*rholy*hx2  &
                            -34.0d0*rhotx*hy2  )

            !   DOING CONTRIB AT            2           2
         nsten =  24
           ss(i,j,nsten) = (   &
                         +     5.0d0*rhory*hx2  &
                         +     5.0d0*rhotx*hy2  )

          !  Now we add in the 2nd order stencil
          ss(i,j,11) = ss(i,j,11) + (                            - betax(i,j))*hx22
          ss(i,j,12) = ss(i,j,12) + (        betax(i+1,j) + 15.0d0*betax(i,j))*hx22
          ss(i,j,0) =  ss(i,j,0 ) + (-15.0d0*betax(i+1,j) - 15.0d0*betax(i,j))*hx22
          ss(i,j,13) = ss(i,j,13) + ( 15.0d0*betax(i+1,j) +        betax(i,j))*hx22
          ss(i,j,14) = ss(i,j,14) + (       -betax(i+1,j)                    )*hx22

          ss(i,j,3) = ss(i,j,3)   + (                            - betay(i,j))*hy22
          ss(i,j,8) = ss(i,j,8)   + (        betay(i,j+1) + 15.0d0*betay(i,j))*hy22
          ss(i,j,0) =  ss(i,j,0 ) + (-15.0d0*betay(i,j+1) - 15.0d0*betay(i,j))*hy22
          ss(i,j,17) = ss(i,j,17) + ( 15.0d0*betay(i,j+1) +        betay(i,j))*hy22
          ss(i,j,22) = ss(i,j,22) + (       -betay(i,j+1)                    )*hy22
       end do
    end do

    ! This adds the "alpha" term in (alpha - del dot beta grad) phi = RHS.
    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          ss(i,j,0) = ss(i,j,0) + alpha(i,j)
       end do
    end do

  end subroutine s_minion_full_fill_2d

  subroutine s_minion_cross_fill_3d(ss, alpha, ng_a, betax, betay, betaz, ng_b, dh, mask, lo, hi)

    integer           , intent(in   ) :: ng_a,ng_b
    integer           , intent(in   ) :: lo(:), hi(:)
    integer           , intent(inout) :: mask(lo(1)      :,lo(2)     :,lo(3)     :)
    real (kind = dp_t), intent(  out) ::    ss(lo(1)     :,lo(2)     :,lo(3)     :,0:)
    real (kind = dp_t), intent(in   ) :: alpha(lo(1)-ng_a:,lo(2)-ng_a:,lo(3)-ng_a:)
    real (kind = dp_t), intent(in   ) :: betax(lo(1)-ng_b:,lo(2)-ng_b:,lo(3)-ng_b:)
    real (kind = dp_t), intent(in   ) :: betay(lo(1)-ng_b:,lo(2)-ng_b:,lo(3)-ng_b:)
    real (kind = dp_t), intent(in   ) :: betaz(lo(1)-ng_b:,lo(2)-ng_b:,lo(3)-ng_b:)
    real (kind = dp_t), intent(in   ) :: dh(:)
    integer i, j, k

    mask = ibclr(mask, BC_BIT(BC_GEOM,1,-1))
    mask = ibclr(mask, BC_BIT(BC_GEOM,1,+1))
    mask = ibclr(mask, BC_BIT(BC_GEOM,2,-1))
    mask = ibclr(mask, BC_BIT(BC_GEOM,2,+1))
    mask = ibclr(mask, BC_BIT(BC_GEOM,3,-1))
    mask = ibclr(mask, BC_BIT(BC_GEOM,3,+1))

    ss = 0.d0
    !
    ! We only include the beta's here to get the viscous coefficients in here for now.
    ! The projection has beta == 1.
    !
    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             ss(i,j,k, 1) =   1.d0 * betax(i  ,j,k)
             ss(i,j,k, 2) = -16.d0 * betax(i  ,j,k)
             ss(i,j,k, 3) = -16.d0 * betax(i+1,j,k)
             ss(i,j,k, 4) =   1.d0 * betax(i+1,j,k)
             ss(i,j,k, 5) =   1.d0 * betay(i,j  ,k)
             ss(i,j,k, 6) = -16.d0 * betay(i,j  ,k)
             ss(i,j,k, 7) = -16.d0 * betay(i,j+1,k)
             ss(i,j,k, 8) =   1.d0 * betay(i,j+1,k)
             ss(i,j,k, 9) =   1.d0 * betaz(i,j,k  )
             ss(i,j,k,10) = -16.d0 * betaz(i,j,k  )
             ss(i,j,k,11) = -16.d0 * betaz(i,j,k+1)
             ss(i,j,k,12) =   1.d0 * betaz(i,j,k+1)
             ss(i,j,k,0)  = -(ss(i,j,k,1) + ss(i,j,k, 2) + ss(i,j,k, 3) + ss(i,j,k, 4) &
                             +ss(i,j,k,5) + ss(i,j,k, 6) + ss(i,j,k, 7) + ss(i,j,k, 8) &
                             +ss(i,j,k,9) + ss(i,j,k,10) + ss(i,j,k,11) + ss(i,j,k,12) )
          end do
       end do
    end do
    !$OMP END PARALLEL DO

    ss = ss * (ONE / (12.d0 * dh(1)**2))
    !
    ! This adds the "alpha" term in (alpha - del dot beta grad) phi = RHS.
    !
    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             ss(i,j,k,0) = ss(i,j,k,0) + alpha(i,j,k)
          end do
       end do
    end do

  end subroutine s_minion_cross_fill_3d

  subroutine s_minion_full_fill_3d(ss, alpha, ng_a, betax, betay, betaz, ng_b, dh, mask, lo, hi)

    integer           , intent(in   ) :: ng_a,ng_b
    integer           , intent(in   ) :: lo(:), hi(:)
    integer           , intent(inout) :: mask(lo(1)      :,lo(2)     :,lo(3)     :)
    real (kind = dp_t), intent(  out) ::    ss(lo(1)     :,lo(2)     :,lo(3)     :,0:)
    real (kind = dp_t), intent(in   ) :: alpha(lo(1)-ng_a:,lo(2)-ng_a:,lo(3)-ng_a:)
    real (kind = dp_t), intent(inout) :: betax(lo(1)-ng_b:,lo(2)-ng_b:,lo(3)-ng_b:)
    real (kind = dp_t), intent(inout) :: betay(lo(1)-ng_b:,lo(2)-ng_b:,lo(3)-ng_b:)
    real (kind = dp_t), intent(inout) :: betaz(lo(1)-ng_b:,lo(2)-ng_b:,lo(3)-ng_b:)
    real (kind = dp_t), intent(in   ) :: dh(:)

    integer          :: i, j, k, nsten
    double precision :: hx2,hy2,hz2,hx22,hy22,hz22

   ! Transverse derivatives: left,right,top,bottom,fore,aft
    double precision :: rhoax,rhobx,rhofx,rhotx
    double precision :: rhoby,rholy,rhory,rhoty
    double precision :: rhofz,rholz,rhorz,rhoaz
    double precision :: s1,s2,scale
    double precision :: sum

    mask = ibclr(mask, BC_BIT(BC_GEOM,1,-1))
    mask = ibclr(mask, BC_BIT(BC_GEOM,1,+1))
    mask = ibclr(mask, BC_BIT(BC_GEOM,2,-1))
    mask = ibclr(mask, BC_BIT(BC_GEOM,2,+1))
    mask = ibclr(mask, BC_BIT(BC_GEOM,3,-1))
    mask = ibclr(mask, BC_BIT(BC_GEOM,3,+1))

    !  These are the coefficients in the second order part
    hx22 = -1.d0 / (12.d0 * dh(1)**2 )
    hy22 = -1.d0 / (12.d0 * dh(2)**2 )
    hz22 = -1.d0 / (12.d0 * dh(3)**2 )

    !  These are the coefficients for the tangential slopes
    !  We don't divide by h because the product of slopes is multiplied by h^2/12
    s1 = 34.d0 / 48.d0
    s2 = -5.d0 / 48.d0

    !  In the output of make_stencil3d.f90, the coefficents of the correction are multiplied by
    !  all the denominators so that they are integers.  So we have to put the denominators back in
    scale = 1.0d0/(12.0d0*48.d0)

    !  The coefficients hx2 and hy2 are defined by  (the minus sign is because it is minus beta*Lap)
    hx2 = -1.0d0/(12.d0*dh(1)**2)*scale
    hy2 = -1.0d0/(12.d0*dh(2)**2)*scale
    hz2 = -1.0d0/(12.d0*dh(3)**2)*scale

    !$OMP PARALLEL DO PRIVATE(i,j,k)

    !  Initialize to zero.
    ss = 0.d0

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)

          rholy = s1*(betax(i  ,j+1,k)- betax(i  ,j-1,k)) + s2*(betax(i  ,j+2,k)- betax(i  ,j-2,k))
          rhory = s1*(betax(i+1,j+1,k)- betax(i+1,j-1,k)) + s2*(betax(i+1,j+2,k)- betax(i+1,j-2,k))
          rholz = s1*(betax(i  ,j,k+1)- betax(i  ,j,k-1)) + s2*(betax(i  ,j,k+2)- betax(i  ,j,k-2))
          rhorz = s1*(betax(i+1,j,k+1)- betax(i+1,j,k-1)) + s2*(betax(i+1,j,k+2)- betax(i+1,j,k-2))

          rhofx = s1*(betay(i+1,j+1,k)- betay(i-1,j+1,k)) + s2*(betay(i+2,j+1,k)- betay(i-2,j+1,k))
          rhoax = s1*(betay(i+1,j  ,k)- betay(i-1,j  ,k)) + s2*(betay(i+2,j  ,k)- betay(i-2,j  ,k))
          rhofz = s1*(betay(i,j+1,k+1)- betay(i,j+1,k-1)) + s2*(betay(i,j+1,k+2)- betay(i,j+1,k-2))
          rhoaz = s1*(betay(i,j  ,k+1)- betay(i,j  ,k-1)) + s2*(betay(i,j  ,k+2)- betay(i,j  ,k-2))

          rhotx = s1*(betaz(i+1,j,k+1)- betaz(i-1,j,k+1)) + s2*(betaz(i+2,j,k+1)- betaz(i-2,j,k+1))
          rhobx = s1*(betaz(i+1,j  ,k)- betaz(i-1,j  ,k)) + s2*(betaz(i+2,j  ,k)- betaz(i-2,j  ,k))
          rhoty = s1*(betaz(i,j+1,k+1)- betaz(i,j-1,k+1)) + s2*(betaz(i,j+2,k+1)- betaz(i,j-2,k+1))
          rhoby = s1*(betaz(i,j+1,k  )- betaz(i,j-1,k  )) + s2*(betaz(i,j+2  ,k)- betaz(i,j-2  ,k))

 ! DOING CONTRIB AT            0          -2          -2 nsten =            1
        nsten =   1
          ss(i,j,k,nsten) = ( &
                        -5.0d0*rhoby*hz2 &
                        -5.0d0*rhoaz*hy2 )
                                              
 ! DOING CONTRIB AT            0          -1          -2 nsten =            2
        nsten =   2
          ss(i,j,k,nsten) = ( &
                   +    34.0d0*rhoby*hz2 &
                   +     5.0d0*rhofz*hy2 &
                   +    75.0d0*rhoaz*hy2 )
                                              
 ! DOING CONTRIB AT           -2           0          -2 nsten =            3
        nsten =   3
          ss(i,j,k,nsten) = ( &
                        -5.0d0*rholz*hx2 &
                        -5.0d0*rhobx*hz2 )
                                              
 ! DOING CONTRIB AT           -1           0          -2 nsten =            4
        nsten =   4
          ss(i,j,k,nsten) = ( &
                   +    75.0d0*rholz*hx2 &
                   +    34.0d0*rhobx*hz2 )
                                              
 ! DOING CONTRIB AT            0           0          -2 nsten =            5
        nsten =   5
          ss(i,j,k,nsten) = ( &
                       -75.0d0*rholz*hx2 &
                       -75.0d0*rhofz*hy2 &
                       -75.0d0*rhoaz*hy2 )
                                              
 ! DOING CONTRIB AT            1           0          -2 nsten =            6
        nsten =   6
          ss(i,j,k,nsten) = ( &
                   +     5.0d0*rholz*hx2 &
                       -34.0d0*rhobx*hz2 )
                                              
 ! DOING CONTRIB AT            2           0          -2 nsten =            7
        nsten =   7
          ss(i,j,k,nsten) = ( &
                   +     5.0d0*rhobx*hz2 )
                                              
 ! DOING CONTRIB AT            0           1          -2 nsten =            8
        nsten =   8
          ss(i,j,k,nsten) = ( &
                       -34.0d0*rhoby*hz2 &
                   +    75.0d0*rhofz*hy2 &
                   +     5.0d0*rhoaz*hy2 )
                                              
 ! DOING CONTRIB AT            0           2          -2 nsten =            9
        nsten =   9
          ss(i,j,k,nsten) = ( &
                   +     5.0d0*rhoby*hz2 &
                        -5.0d0*rhofz*hy2 )
                                              
 ! DOING CONTRIB AT            0          -2          -1 nsten =           10
        nsten =  10
          ss(i,j,k,nsten) = ( &
                   +     5.0d0*rhoty*hz2 &
                   +    75.0d0*rhoby*hz2 &
                   +    34.0d0*rhoaz*hy2 )
                                              
 ! DOING CONTRIB AT            0          -1          -1 nsten =           11
        nsten =  11
          ss(i,j,k,nsten) = ( &
                       -34.0d0*rhoty*hz2 &
                      -510.0d0*rhoby*hz2 &
                       -34.0d0*rhofz*hy2 &
                      -510.0d0*rhoaz*hy2 )
                                              
 ! DOING CONTRIB AT           -2           0          -1 nsten =           12
        nsten =  12
          ss(i,j,k,nsten) = ( &
                   +    34.0d0*rholz*hx2 &
                   +     5.0d0*rhotx*hz2 &
                   +    75.0d0*rhobx*hz2 )
                                              
 ! DOING CONTRIB AT           -1           0          -1 nsten =           13
        nsten =  13
          ss(i,j,k,nsten) = ( &
                      -510.0d0*rholz*hx2 &
                       -34.0d0*rhotx*hz2 &
                      -510.0d0*rhobx*hz2 )
                                              
 ! DOING CONTRIB AT            0           0          -1 nsten =           14
        nsten =  14
          ss(i,j,k,nsten) = ( &
                   +   510.0d0*rholz*hx2 &
                   +   510.0d0*rhofz*hy2 &
                   +   510.0d0*rhoaz*hy2 )
                                              
 ! DOING CONTRIB AT            1           0          -1 nsten =           15
        nsten =  15
          ss(i,j,k,nsten) = ( &
                       -34.0d0*rholz*hx2 &
                   +    34.0d0*rhotx*hz2 &
                   +   510.0d0*rhobx*hz2 )
                                              
 ! DOING CONTRIB AT            2           0          -1 nsten =           16
        nsten =  16
          ss(i,j,k,nsten) = ( &
                        -5.0d0*rhotx*hz2 &
                       -75.0d0*rhobx*hz2 )
                                              
 ! DOING CONTRIB AT            0           1          -1 nsten =           17
        nsten =  17
          ss(i,j,k,nsten) = ( &
                   +    34.0d0*rhoty*hz2 &
                   +   510.0d0*rhoby*hz2 &
                      -510.0d0*rhofz*hy2 &
                       -34.0d0*rhoaz*hy2 )
                                              
 ! DOING CONTRIB AT            0           2          -1 nsten =           18
        nsten =  18
          ss(i,j,k,nsten) = ( &
                        -5.0d0*rhoty*hz2 &
                       -75.0d0*rhoby*hz2 &
                   +    34.0d0*rhofz*hy2 )
                                              
 ! DOING CONTRIB AT           -2          -2           0 nsten =           19
        nsten =  19
          ss(i,j,k,nsten) = ( &
                        -5.0d0*rholy*hx2 &
                        -5.0d0*rhoax*hy2 )
                                              
 ! DOING CONTRIB AT           -1          -2           0 nsten =           20
        nsten =  20
          ss(i,j,k,nsten) = ( &
                   +     5.0d0*rhory*hx2 &
                   +     5.0d0*rhorz*hx2 &
                   +    75.0d0*rholy*hx2 &
                   +    34.0d0*rhoax*hy2 )
                                              
 ! DOING CONTRIB AT            0          -2           0 nsten =           21
        nsten =  21
          ss(i,j,k,nsten) = ( &
                       -75.0d0*rhory*hx2 &
                       -75.0d0*rhorz*hx2 &
                       -75.0d0*rholy*hx2 &
                       -75.0d0*rhoty*hz2 &
                       -75.0d0*rhoby*hz2 )
                                              
 ! DOING CONTRIB AT            1          -2           0 nsten =           22
        nsten =  22
          ss(i,j,k,nsten) = ( &
                   +    75.0d0*rhory*hx2 &
                   +    75.0d0*rhorz*hx2 &
                   +     5.0d0*rholy*hx2 &
                       -34.0d0*rhoax*hy2 )
                                              
 ! DOING CONTRIB AT            2          -2           0 nsten =           23
        nsten =  23
          ss(i,j,k,nsten) = ( &
                        -5.0d0*rhory*hx2 &
                        -5.0d0*rhorz*hx2 &
                   +     5.0d0*rhoax*hy2 )
                                              
 ! DOING CONTRIB AT           -2          -1           0 nsten =           24
        nsten =  24
          ss(i,j,k,nsten) = ( &
                   +    34.0d0*rholy*hx2 &
                   +     5.0d0*rhofx*hy2 &
                   +    75.0d0*rhoax*hy2 )
                                              
 ! DOING CONTRIB AT           -1          -1           0 nsten =           25
        nsten =  25
          ss(i,j,k,nsten) = ( &
                       -34.0d0*rhory*hx2 &
                       -34.0d0*rhorz*hx2 &
                      -510.0d0*rholy*hx2 &
                       -34.0d0*rhofx*hy2 &
                      -510.0d0*rhoax*hy2 )
                                              
 ! DOING CONTRIB AT            0          -1           0 nsten =           26
        nsten =  26
          ss(i,j,k,nsten) = ( &
                   +   510.0d0*rhory*hx2 &
                   +   510.0d0*rhorz*hx2 &
                   +   510.0d0*rholy*hx2 &
                   +   510.0d0*rhoty*hz2 &
                   +   510.0d0*rhoby*hz2 )
                                              
 ! DOING CONTRIB AT            1          -1           0 nsten =           27
        nsten =  27
          ss(i,j,k,nsten) = ( &
                      -510.0d0*rhory*hx2 &
                      -510.0d0*rhorz*hx2 &
                       -34.0d0*rholy*hx2 &
                   +    34.0d0*rhofx*hy2 &
                   +   510.0d0*rhoax*hy2 )
                                              
 ! DOING CONTRIB AT            2          -1           0 nsten =           28
        nsten =  28
          ss(i,j,k,nsten) = ( &
                   +    34.0d0*rhory*hx2 &
                   +    34.0d0*rhorz*hx2 &
                        -5.0d0*rhofx*hy2 &
                       -75.0d0*rhoax*hy2 )
                                              
 ! DOING CONTRIB AT           -2           0           0 nsten =           29
        nsten =  29
          ss(i,j,k,nsten) = ( &
                       -75.0d0*rhotx*hz2 &
                       -75.0d0*rhobx*hz2 &
                       -75.0d0*rhofx*hy2 &
                       -75.0d0*rhoax*hy2 )
                                              
 ! DOING CONTRIB AT           -1           0           0 nsten =           30
        nsten =  30
          ss(i,j,k,nsten) = ( &
                   +   510.0d0*rhotx*hz2 &
                   +   510.0d0*rhobx*hz2 &
                   +   510.0d0*rhofx*hy2 &
                   +   510.0d0*rhoax*hy2 )
                                              
 ! DOING CONTRIB AT            1           0           0 nsten =           31
        nsten =  31
          ss(i,j,k,nsten) = ( &
                      -510.0d0*rhotx*hz2 &
                      -510.0d0*rhobx*hz2 &
                      -510.0d0*rhofx*hy2 &
                      -510.0d0*rhoax*hy2 )
                                              
 ! DOING CONTRIB AT            2           0           0 nsten =           32
        nsten =  32
          ss(i,j,k,nsten) = ( &
                   +    75.0d0*rhotx*hz2 &
                   +    75.0d0*rhobx*hz2 &
                   +    75.0d0*rhofx*hy2 &
                   +    75.0d0*rhoax*hy2 )
                                              
 ! DOING CONTRIB AT           -2           1           0 nsten =           33
        nsten =  33
          ss(i,j,k,nsten) = ( &
                       -34.0d0*rholy*hx2 &
                   +    75.0d0*rhofx*hy2 &
                   +     5.0d0*rhoax*hy2 )
                                              
 ! DOING CONTRIB AT           -1           1           0 nsten =           34
        nsten =  34
          ss(i,j,k,nsten) = ( &
                   +    34.0d0*rhory*hx2 &
                   +    34.0d0*rhorz*hx2 &
                   +   510.0d0*rholy*hx2 &
                      -510.0d0*rhofx*hy2 &
                       -34.0d0*rhoax*hy2 )
                                              
 ! DOING CONTRIB AT            0           1           0 nsten =           35
        nsten =  35
          ss(i,j,k,nsten) = ( &
                      -510.0d0*rhory*hx2 &
                      -510.0d0*rhorz*hx2 &
                      -510.0d0*rholy*hx2 &
                      -510.0d0*rhoty*hz2 &
                      -510.0d0*rhoby*hz2 )
                                              
 ! DOING CONTRIB AT            1           1           0 nsten =           36
        nsten =  36
          ss(i,j,k,nsten) = ( &
                   +   510.0d0*rhory*hx2 &
                   +   510.0d0*rhorz*hx2 &
                   +    34.0d0*rholy*hx2 &
                   +   510.0d0*rhofx*hy2 &
                   +    34.0d0*rhoax*hy2 )
                                              
 ! DOING CONTRIB AT            2           1           0 nsten =           37
        nsten =  37
          ss(i,j,k,nsten) = ( &
                       -34.0d0*rhory*hx2 &
                       -34.0d0*rhorz*hx2 &
                       -75.0d0*rhofx*hy2 &
                        -5.0d0*rhoax*hy2 )
                                              
 ! DOING CONTRIB AT           -2           2           0 nsten =           38
        nsten =  38
          ss(i,j,k,nsten) = ( &
                   +     5.0d0*rholy*hx2 &
                        -5.0d0*rhofx*hy2 )
                                              
 ! DOING CONTRIB AT           -1           2           0 nsten =           39
        nsten =  39
          ss(i,j,k,nsten) = ( &
                        -5.0d0*rhory*hx2 &
                        -5.0d0*rhorz*hx2 &
                       -75.0d0*rholy*hx2 &
                   +    34.0d0*rhofx*hy2 )
                                              
 ! DOING CONTRIB AT            0           2           0 nsten =           40
        nsten =  40
          ss(i,j,k,nsten) = ( &
                   +    75.0d0*rhory*hx2 &
                   +    75.0d0*rhorz*hx2 &
                   +    75.0d0*rholy*hx2 &
                   +    75.0d0*rhoty*hz2 &
                   +    75.0d0*rhoby*hz2 )
                                              
 ! DOING CONTRIB AT            1           2           0 nsten =           41
        nsten =  41
          ss(i,j,k,nsten) = ( &
                       -75.0d0*rhory*hx2 &
                       -75.0d0*rhorz*hx2 &
                        -5.0d0*rholy*hx2 &
                       -34.0d0*rhofx*hy2 )
                                              
 ! DOING CONTRIB AT            2           2           0 nsten =           42
        nsten =  42
          ss(i,j,k,nsten) = ( &
                   +     5.0d0*rhory*hx2 &
                   +     5.0d0*rhorz*hx2 &
                   +     5.0d0*rhofx*hy2 )
                                              
 ! DOING CONTRIB AT            0          -2           1 nsten =           43
        nsten =  43
          ss(i,j,k,nsten) = ( &
                   +    75.0d0*rhoty*hz2 &
                   +     5.0d0*rhoby*hz2 &
                       -34.0d0*rhoaz*hy2 )
                                              
 ! DOING CONTRIB AT            0          -1           1 nsten =           44
        nsten =  44
          ss(i,j,k,nsten) = ( &
                      -510.0d0*rhoty*hz2 &
                       -34.0d0*rhoby*hz2 &
                   +    34.0d0*rhofz*hy2 &
                   +   510.0d0*rhoaz*hy2 )
                                              
 ! DOING CONTRIB AT           -2           0           1 nsten =           45
        nsten =  45
          ss(i,j,k,nsten) = ( &
                       -34.0d0*rholz*hx2 &
                   +    75.0d0*rhotx*hz2 &
                   +     5.0d0*rhobx*hz2 )
                                              
 ! DOING CONTRIB AT           -1           0           1 nsten =           46
        nsten =  46
          ss(i,j,k,nsten) = ( &
                   +   510.0d0*rholz*hx2 &
                      -510.0d0*rhotx*hz2 &
                       -34.0d0*rhobx*hz2 )
                                              
 ! DOING CONTRIB AT            0           0           1 nsten =           47
        nsten =  47
          ss(i,j,k,nsten) = ( &
                      -510.0d0*rholz*hx2 &
                      -510.0d0*rhofz*hy2 &
                      -510.0d0*rhoaz*hy2 )
                                              
 ! DOING CONTRIB AT            1           0           1 nsten =           48
        nsten =  48
          ss(i,j,k,nsten) = ( &
                   +    34.0d0*rholz*hx2 &
                   +   510.0d0*rhotx*hz2 &
                   +    34.0d0*rhobx*hz2 )
                                              
 ! DOING CONTRIB AT            2           0           1 nsten =           49
        nsten =  49
          ss(i,j,k,nsten) = ( &
                       -75.0d0*rhotx*hz2 &
                        -5.0d0*rhobx*hz2 )
                                              
 ! DOING CONTRIB AT            0           1           1 nsten =           50
        nsten =  50
          ss(i,j,k,nsten) = ( &
                   +   510.0d0*rhoty*hz2 &
                   +    34.0d0*rhoby*hz2 &
                   +   510.0d0*rhofz*hy2 &
                   +    34.0d0*rhoaz*hy2 )
                                              
 ! DOING CONTRIB AT            0           2           1 nsten =           51
        nsten =  51
          ss(i,j,k,nsten) = ( &
                       -75.0d0*rhoty*hz2 &
                        -5.0d0*rhoby*hz2 &
                       -34.0d0*rhofz*hy2 )
                                              
 ! DOING CONTRIB AT            0          -2           2 nsten =           52
        nsten =  52
          ss(i,j,k,nsten) = ( &
                        -5.0d0*rhoty*hz2 &
                   +     5.0d0*rhoaz*hy2 )
                                              
 ! DOING CONTRIB AT            0          -1           2 nsten =           53
        nsten =  53
          ss(i,j,k,nsten) = ( &
                   +    34.0d0*rhoty*hz2 &
                        -5.0d0*rhofz*hy2 &
                       -75.0d0*rhoaz*hy2 )
                                              
 ! DOING CONTRIB AT           -2           0           2 nsten =           54
        nsten =  54
          ss(i,j,k,nsten) = ( &
                   +     5.0d0*rholz*hx2 &
                        -5.0d0*rhotx*hz2 )
                                              
 ! DOING CONTRIB AT           -1           0           2 nsten =           55
        nsten =  55
          ss(i,j,k,nsten) = ( &
                       -75.0d0*rholz*hx2 &
                   +    34.0d0*rhotx*hz2 )
                                              
 ! DOING CONTRIB AT            0           0           2 nsten =           56
        nsten =  56
          ss(i,j,k,nsten) = ( &
                   +    75.0d0*rholz*hx2 &
                   +    75.0d0*rhofz*hy2 &
                   +    75.0d0*rhoaz*hy2 )
                                              
 ! DOING CONTRIB AT            1           0           2 nsten =           57
        nsten =  57
          ss(i,j,k,nsten) = ( &
                        -5.0d0*rholz*hx2 &
                       -34.0d0*rhotx*hz2 )
                                              
 ! DOING CONTRIB AT            2           0           2 nsten =           58
        nsten =  58
          ss(i,j,k,nsten) = ( &
                   +     5.0d0*rhotx*hz2 )
                                              
 ! DOING CONTRIB AT            0           1           2 nsten =           59
        nsten =  59
          ss(i,j,k,nsten) = ( &
                       -34.0d0*rhoty*hz2 &
                       -75.0d0*rhofz*hy2 &
                        -5.0d0*rhoaz*hy2 )
                                              
 ! DOING CONTRIB AT            0           2           2 nsten =           60
        nsten =  60
          ss(i,j,k,nsten) = ( &
                   +     5.0d0*rhoty*hz2 &
                   +     5.0d0*rhofz*hy2 )
                                              

          !  Now we add in the 2nd order stencil
          ss(i,j,k,29) = ss(i,j,k,29) + (                       -        betax(i,j,k))*hx22
          ss(i,j,k,30) = ss(i,j,k,30) + (        betax(i+1,j,k) + 15.0d0*betax(i,j,k))*hx22
          ! ss(i,j,k, 0) = ss(i,j,k, 0) + (-15.0d0*betax(i+1,j,k) - 15.0d0*betax(i,j,k))*hx22
          ss(i,j,k,31) = ss(i,j,k,31) + ( 15.0d0*betax(i+1,j,k) +        betax(i,j,k))*hx22
          ss(i,j,k,32) = ss(i,j,k,32) + (       -betax(i+1,j,k)                      )*hx22
 
          ss(i,j,k,21) = ss(i,j,k,21) + (                      -        betay(i,j,k))*hy22
          ss(i,j,k,26) = ss(i,j,k,26) + (        betay(i,j+1,k) + 15.0d0*betay(i,j,k))*hy22
          ! ss(i,j,k, 0) = ss(i,j,k, 0) + (-15.0d0*betay(i,j+1,k) - 15.0d0*betay(i,j,k))*hy22
          ss(i,j,k,35) = ss(i,j,k,35) + ( 15.0d0*betay(i,j+1,k) +        betay(i,j,k))*hy22
          ss(i,j,k,40) = ss(i,j,k,40) + (       -betay(i,j+1,k)                      )*hy22
 
          ss(i,j,k, 5) = ss(i,j,k, 5) + (                       -        betaz(i,j,k))*hz22
          ss(i,j,k,14) = ss(i,j,k,14) + (        betaz(i,j,k+1) + 15.0d0*betaz(i,j,k))*hz22
          ! ss(i,j,k, 0) = ss(i,j,k, 0) + (-15.0d0*betaz(i,j,k+1) - 15.0d0*betaz(i,j,k))*hz22
          ss(i,j,k,47) = ss(i,j,k,47) + ( 15.0d0*betaz(i,j,k+1) +        betaz(i,j,k))*hz22
          ss(i,j,k,56) = ss(i,j,k,56) + (       -betaz(i,j,k+1)                      )*hz22

          end do
       end do
    end do

    !$OMP END PARALLEL DO

    !
    ! Define the center stencil.
    !
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
          sum = 0.d0
          do nsten = 1,60
             sum = sum + ss(i,j,k,nsten)
          end do
          ss(i,j,k,0) = -sum
          end do
       end do
    end do

    !
    ! This adds the "alpha" term in (alpha - del dot beta grad) phi = RHS.
    !
    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             ss(i,j,k,0) = ss(i,j,k,0) + alpha(i,j,k)
          end do
       end do
    end do

  end subroutine s_minion_full_fill_3d

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

    j = 0
    
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

end module cc_stencil_module
