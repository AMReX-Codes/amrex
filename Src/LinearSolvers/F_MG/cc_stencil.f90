module cc_stencil_module

  use bl_types
  use bl_constants_module
  use bc_module
  use bc_functions_module
  use multifab_module
  use stencil_util_module, only : is_ibc_stencil

  implicit none

  private
  public :: stencil_norm, max_of_stencil_sum, stencil_set_bc, &
       simple_2d_const, simple_3d_const, s_simple_1d_cc, s_simple_2d_cc, s_simple_3d_cc, &
       s_minion_cross_fill_2d, s_minion_cross_fill_3d, s_minion_full_fill_2d, s_minion_full_fill_3d, &
       s_simplen_2d_cc, s_simplem_2d_cc, s_simpleg_2d_cc, s_simpleg_3d_cc

contains

  function stencil_norm(ss, mask, local) result(r)
    use bl_prof_module
    real(kind=dp_t) :: r
    type(multifab), intent(in) :: ss
    type(lmultifab), intent(in), optional :: mask
    logical, intent(in), optional :: local
    integer :: i,j,k,n,b,lo(4),hi(4)
    real(kind=dp_t) :: r1, sum_comps
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    logical, pointer :: lp(:,:,:,:)
    logical :: llocal
    type(bl_prof_timer), save :: bpt

    call build(bpt, "st_norm_cc")

    llocal = .false.; if ( present(local) ) llocal = local

    r1 = -Huge(r1)

    if ( present(mask) ) then
       !$OMP PARALLEL DO PRIVATE(i,j,k,n,sum_comps,sp,lp,lo,hi) REDUCTION(max:r1)
       do b = 1, nfabs(ss)
          sp => dataptr(ss, b)
          lp => dataptr(mask, b)

          if (is_ibc_stencil(ss,b)) then
             if (any(lp .eqv. .true.)) then
                r1 = max(r1, abs(sp(1,1,1,1))+2.d0*sum(abs(sp(2:,1,1,1))))
             end if
          else
             lo =  lbound(sp)
             hi =  ubound(sp)
             do k = lo(4), hi(4)
                do j = lo(3), hi(3)
                   do i = lo(2), hi(2)
                      if ( lp(i,j,k,1) ) then
                         sum_comps = ZERO
                         do n = lo(1), hi(1)
                            sum_comps = sum_comps + abs(sp(n,i,j,k))
                         end do
                         r1 = max(r1,sum_comps)
                      end if
                   end do
                end do
             end do
          end if
       end do
       !$OMP END PARALLEL DO
    else
       !$OMP PARALLEL DO PRIVATE(i,j,k,n,sum_comps,sp,lo,hi) REDUCTION(max:r1)
       do b = 1, nfabs(ss)
          sp => dataptr(ss, b)
          
          if (is_ibc_stencil(ss,b)) then
             r1 = max(r1, abs(sp(1,1,1,1))+2.d0*sum(abs(sp(2:,1,1,1))))
          else
             lo =  lbound(sp)
             hi =  ubound(sp)
             do k = lo(4), hi(4)
                do j = lo(3), hi(3)
                   do i = lo(2), hi(2)
                      sum_comps = ZERO
                      do n = lo(1), hi(1)
                         sum_comps = sum_comps + abs(sp(n,i,j,k))
                      end do
                      r1 = max(r1,sum_comps)
                   end do
                end do
             end do
          end if
       end do
       !$OMP END PARALLEL DO
    end if

    r = r1

    if ( .not. llocal ) call parallel_reduce(r, r1, MPI_MAX)

    call destroy(bpt)
  end function stencil_norm

  function max_of_stencil_sum(ss, mask, local) result(r)
    use bl_prof_module
    real(kind=dp_t) :: r
    type(multifab), intent(in) :: ss
    type(lmultifab), intent(in), optional :: mask
    logical, intent(in), optional :: local
    integer :: i,j,k,n,b,lo(4),hi(4)
    real(kind=dp_t) :: r1, sum_comps
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    logical, pointer :: lp(:,:,:,:)
    logical :: llocal
    type(bl_prof_timer), save :: bpt

    llocal = .false.; if ( present(local) ) llocal = local
    !
    ! NOTE: this is exactly the same as the stencil_norm function except that we sum the
    ! components of the stencil, not the absolute value of each component.
    !
    call build(bpt, "st_sum")

    r1 = -Huge(r1)

    if ( present(mask) ) then
       !$OMP PARALLEL DO PRIVATE(i,j,k,n,sum_comps,sp,lp,lo,hi) REDUCTION(max:r1)
       do b = 1, nfabs(ss)
          sp => dataptr(ss, b)
          lp => dataptr(mask, b)

          if (is_ibc_stencil(ss,b)) then
             if (any(lp .eqv. .true.)) then
                r1 = max(r1, sp(1,1,1,1)+2.d0*sum(sp(2:,1,1,1)))
             end if
          else
             lo =  lbound(sp)
             hi =  ubound(sp)
             do k = lo(4), hi(4)
                do j = lo(3), hi(3)
                   do i = lo(2), hi(2)
                      if ( lp(i,j,k,1) ) then
                         sum_comps = ZERO
                         do n = lo(1), hi(1)
                            sum_comps = sum_comps + sp(n,i,j,k)
                         end do
                         r1 = max(r1,sum_comps)
                      end if
                   end do
                end do
             end do
          end if
       end do
       !$OMP END PARALLEL DO
    else
       !$OMP PARALLEL DO PRIVATE(i,j,k,n,sum_comps,sp,lo,hi) REDUCTION(max:r1)
       do b = 1, nfabs(ss)
          sp => dataptr(ss, b)

          if (is_ibc_stencil(ss,b)) then
             r1 = max(r1, sp(1,1,1,1)+2.d0*sum(sp(2:,1,1,1)))
          else
             lo =  lbound(sp)
             hi =  ubound(sp)
             do k = lo(4), hi(4)
                do j = lo(3), hi(3)
                   do i = lo(2), hi(2)
                      sum_comps = ZERO
                      do n = lo(1), hi(1)
                         sum_comps = sum_comps + sp(n,i,j,k)
                      end do
                      r1 = max(r1,sum_comps)
                   end do
                end do
             end do
          end if
       end do
       !$OMP END PARALLEL DO
    end if

    r = r1

    if ( .not. llocal ) call parallel_reduce(r, r1, MPI_MAX)

    call destroy(bpt)
  end function max_of_stencil_sum

  subroutine stencil_set_bc(st, idx, mask, bc_face, cf_face, intbox)
    type(multifab),  intent(in)           :: st
    integer,         intent(in)           :: idx
    type(imultifab), intent(inout)        :: mask
    integer,         intent(in)           :: bc_face(:,:)
    integer,         intent(in), optional :: cf_face(:,:)
    logical,         intent(out), optional:: intbox

    type(box)        :: bx1, pd
    type(boxarray)   :: ba
    type(layout)     :: la
    integer          :: i, j, ii, jj, ldom, blo, bhi
    integer, pointer :: mp(:,:,:,:)
    integer          :: lcf_face(size(bc_face, 1), size(bc_face, 2))
    logical          :: pmask(get_dim(st)), intflag
    !
    ! The Coarse-Fine boundary is Dirichlet unless specified.
    !
    lcf_face = BC_DIR; if ( present(cf_face) ) lcf_face = cf_face
    !
    ! Initialize every border to Fine-Fine (Interior).
    !
    mp => dataptr(mask, idx)
    mp = BC_INT

    la    = get_layout(st)
    pd    = get_pd(get_layout(st))
    pmask = get_pmask(get_layout(st))

    intflag = .true.

    do i = 1, get_dim(st)
       do j = -1, 1, 2
          bx1 = get_box(st,idx)
          if (j .eq. -1) then
             jj = 1
             blo = box_lwb_d(bx1, i)
             call box_set_lwb_d(bx1, i, blo-1)
             call box_set_upb_d(bx1, i, blo-1)
          else
             jj = 2
             bhi = box_upb_d(bx1, i)
             call box_set_lwb_d(bx1, i, bhi+1)
             call box_set_upb_d(bx1, i, bhi+1)
          end if
          if ( contains(pd, bx1) ) then
             !
             ! We're not touching a physical boundary -- set any/all C-F bndrys.
             !
             call layout_boxarray_diff(ba, bx1, la)
             do ii = 1, nboxes(ba)
                bx1 = shift(get_box(ba,ii), -j, i)
                mp => dataptr(mask, idx, bx1)
                mp = ibset(mp, BC_BIT(lcf_face(i, jj), i, j))
                intflag = .false.
             end do
             call boxarray_destroy(ba)
          else
             !
             ! We touch a physical boundary in that direction.
             !
             if ( .not. pmask(i) ) then
                !
                ! We're not periodic in that direction -- use physical BCs.
                !
                bx1 = shift(bx1, -j, i)
                mp => dataptr(mask, idx, bx1)
                mp = ibset(mp, BC_BIT(bc_face(i, jj), i, j))
                intflag = .false.
             else
                !
                ! Remove any/all Fine-Fine intersections.
                !
                ldom = extent(pd, i)
                bx1 = shift(bx1, -j*ldom, i)
                call layout_boxarray_diff(ba, bx1, la)
                do ii = 1, nboxes(ba)
                   bx1 = shift(get_box(ba,ii), j*ldom-j, i)
                   mp => dataptr(mask, idx, bx1)
                   mp = ibset(mp, BC_BIT(lcf_face(i, jj), i, j))
                   intflag = .false.
                end do
                call boxarray_destroy(ba)
             end if
          end if
       end do
    end do

    if (present(intbox)) intbox = intflag

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
    real (kind = dp_t), intent(  out) ::   ss(0:,lo(1):)
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
       ss(0,i) =   ZERO
       ss(1,i) = -betax(i+1)*f1(1)
       ss(2,i) = -betax(i  )*f1(1)
       ss(XBC,i) = ZERO
    end do

    ! x derivatives

    do i = lo(1)+1, hi(1)-1
       ss(0,i) = ss(0,i) + (betax(i+1)+betax(i))*f1(1)
    end do

    bclo = stencil_bc_type(mask(lo(1)),1,-1)
    bchi = stencil_bc_type(mask(hi(1)),1,+1)

    i = lo(1)
    if (bclo .eq. BC_INT) then
       ss(0,i) = ss(0,i) + (betax(i)+betax(i+1))*f1(1)
    else
       call stencil_bndry_aaa(order, nx, 1, -1, mask(i), &
            ss(0,i), ss(1,i), ss(2,i), ss(XBC,i), &
            betax(i), betax(i+1), xa(1), xb(1), dh(1), bclo, bchi)
    end if

    if ( hi(1) > lo(1) ) then
       i = hi(1)
       if (bchi .eq. BC_INT) then
          ss(0,i) = ss(0,i) + (betax(i)+betax(i+1))*f1(1)
       else
          call stencil_bndry_aaa(order, nx, 1, 1, mask(i), &
               ss(0,i), ss(1,i), ss(2,i), ss(XBC,i), &
               betax(i), betax(i+1), xa(1), xb(1), dh(1), bclo, bchi)
       end if
    end if

    do i = lo(1),hi(1)
       ss(0,i) = ss(0,i) + alpha(i)
    end do

  end subroutine s_simple_1d_cc

  subroutine s_simple_2d_cc(ss, alpha, ng_a, betax, betay, ng_b, dh, mask, lo, hi, xa, xb, order)

    integer           , intent(in   ) :: ng_a, ng_b, lo(:), hi(:), order
    integer           , intent(inout) :: mask(lo(1)  :,lo(2)  :)
    real (kind = dp_t), intent(  out) ::   ss(0:, lo(1)  :,lo(2)  :)
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
          ss(0,i,j) = ZERO
          ss(1,i,j) = -betax(i+1,j)*f1(1)
          ss(2,i,j) = -betax(i  ,j)*f1(1)
          ss(3,i,j) = -betay(i,j+1)*f1(2)
          ss(4,i,j) = -betay(i,j  )*f1(2)
          ss(XBC,i,j) = ZERO
          ss(YBC,i,j) = ZERO
       end do
    end do

    ! x derivatives

    do j = lo(2),hi(2)
       do i = lo(1)+1,hi(1)-1
          ss(0,i,j) = ss(0,i,j) + (betax(i,j)+betax(i+1,j))*f1(1)
       end do
    end do

    do j = lo(2),hi(2)
       bclo = stencil_bc_type(mask(lo(1),j),1,-1)
       bchi = stencil_bc_type(mask(hi(1),j),1,+1)
 
       i = lo(1)
       if (bclo .eq. BC_INT) then
          ss(0,i,j) = ss(0,i,j) + (betax(i,j)+betax(i+1,j))*f1(1)
       else
          call stencil_bndry_aaa(order, nx, 1, -1, mask(i,j), &
               ss(0,i,j), ss(1,i,j), ss(2,i,j), ss(XBC,i,j), &
               betax(i,j), betax(i+1,j), &
               xa(1), xb(1), dh(1), bclo, bchi)
       end if

       if ( hi(1) > lo(1) ) then
          i = hi(1)
          if (bchi .eq. BC_INT) then
             ss(0,i,j) = ss(0,i,j) + (betax(i,j)+betax(i+1,j))*f1(1)
          else
             call stencil_bndry_aaa(order, nx, 1, 1, mask(i,j), &
                  ss(0,i,j), ss(1,i,j), ss(2,i,j), ss(XBC,i,j), &
                  betax(i,j), betax(i+1,j), &
                  xa(1), xb(1), dh(1), bclo, bchi) 
          end if 
       end if
    end do

    ! y derivatives

    do i = lo(1),hi(1)
       do j = lo(2)+1,hi(2)-1
          ss(0,i,j) = ss(0,i,j) + (betay(i,j)+betay(i,j+1))*f1(2)
       end do
    end do

    do i = lo(1),hi(1)
       bclo = stencil_bc_type(mask( i,lo(2)),2,-1)
       bchi = stencil_bc_type(mask( i,hi(2)),2,+1)

       j = lo(2)
       if (bclo .eq. BC_INT) then
          ss(0,i,j) = ss(0,i,j) + (betay(i,j)+betay(i,j+1))*f1(2)
       else
          call stencil_bndry_aaa(order, ny, 2, -1, mask(i,j), &
               ss(0,i,j), ss(3,i,j), ss(4,i,j),ss(YBC,i,j), &
               betay(i,j), betay(i,j+1), &
               xa(2), xb(2), dh(2), bclo, bchi)
       end if

       if ( hi(2) > lo(2) ) then
          j = hi(2)
          if (bchi .eq. BC_INT) then
             ss(0,i,j) = ss(0,i,j) + (betay(i,j)+betay(i,j+1))*f1(2)
          else
             call stencil_bndry_aaa(order, ny, 2, 1, mask(i,j), &
                  ss(0,i,j), ss(3,i,j), ss(4,i,j), ss(YBC,i,j), &
                  betay(i,j), betay(i,j+1), &
                  xa(2), xb(2), dh(2), bclo, bchi)
          end if
       end if
    end do

    do j = lo(2),hi(2)
       do i = lo(1),hi(1)
          ss(0,i,j) = ss(0,i,j) + alpha(i,j)
       end do
    end do

  end subroutine s_simple_2d_cc

  subroutine s_simplen_2d_cc(ss, alpha, ng_a, betax, betay, ng_b, dh, mask, lo, hi, xa, xb, order)

    integer           , intent(in   ) :: ng_a, ng_b, lo(:), hi(:), order
    integer           , intent(inout) :: mask(lo(1)  :,lo(2)  :)
    real (kind = dp_t), intent(  out) ::   ss(0:, lo(1)  :,lo(2)  :)
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
 
    ss = zero

    ! Consider the operator  ( alpha - sum_n (beta0_n del dot beta_n grad) )
    ! Components alpha(i,j,   0) = alpha
    ! Components alpha(i,j,1:nc) = beta0_n
    ! Components betax(i,j,1:nc) = betax_n
    ! Components betay(i,j,1:nc) = betay_n

    ! ss(1,i,j) is the coefficient of phi(i+1,j  )
    ! ss(2,i,j) is the coefficient of phi(i-1,j  )
    ! ss(3,i,j) is the coefficient of phi(i  ,j+1)
    ! ss(4,i,j) is the coefficient of phi(i  ,j-1)
    ! ss(0,i,j) is the coefficient of phi(i  ,j  )
    
    do j = lo(2),hi(2)
       do i = lo(1),hi(1)
          do n = 1,nc
             ss(1,i,j) = ss(1,i,j) - betax(i+1,j,n)*f1(1) / alpha(i,j,n)
             ss(2,i,j) = ss(2,i,j) - betax(i  ,j,n)*f1(1) / alpha(i,j,n)
             ss(3,i,j) = ss(3,i,j) - betay(i,j+1,n)*f1(2) / alpha(i,j,n)
             ss(4,i,j) = ss(4,i,j) - betay(i,j  ,n)*f1(2) / alpha(i,j,n)
          end do 
       end do
    end do

    ! x derivatives

    do j = lo(2),hi(2)
       do i = lo(1)+1,hi(1)-1
          do n = 1, nc
             ss(0,i,j) = ss(0,i,j) + (betax(i,j,n)+betax(i+1,j,n))*f1(1) / alpha(i,j,n)
          end do
       end do
    end do

    do j = lo(2),hi(2)
       bclo = stencil_bc_type(mask(lo(1),j),1,-1)
       bchi = stencil_bc_type(mask(hi(1),j),1,+1)
 
       i = lo(1)
       if (bclo .eq. BC_INT) then
          do n = 1, nc
             ss(0,i,j) = ss(0,i,j) + (betax(i,j,n)+betax(i+1,j,n))*f1(1) / alpha(i,j,n)
          end do
       else
          blo = zero
          bhi = zero
          do n = 1,nc
            blo = blo + betax(i  ,j,n) / alpha(i,j,n)
            bhi = bhi + betax(i+1,j,n) / alpha(i,j,n)
          end do
          call stencil_bndry_aaa(order, nx, 1, -1, mask(i,j), &
               ss(0,i,j), ss(1,i,j), ss(2,i,j), ss(XBC,i,j), &
               blo, bhi, xa(1), xb(1), dh(1), bclo, bchi)
       end if

       if ( hi(1) > lo(1) ) then
          i = hi(1)
          if (bchi .eq. BC_INT) then
             do n = 1, nc
                ss(0,i,j) = ss(0,i,j) + (betax(i,j,n)+betax(i+1,j,n))*f1(1) / alpha(i,j,n)
             end do
          else
             blo = zero
             bhi = zero
             do n = 1,nc
                blo = blo + betax(i  ,j,n) / alpha(i,j,n)
                bhi = bhi + betax(i+1,j,n) / alpha(i,j,n)
             end do
             call stencil_bndry_aaa(order, nx, 1, 1, mask(i,j), &
                  ss(0,i,j), ss(1,i,j), ss(2,i,j), ss(XBC,i,j), &
                  blo, bhi, xa(1), xb(1), dh(1), bclo, bchi)
          end if
       end if
    end do

    ! y derivatives

    do i = lo(1),hi(1)
       do j = lo(2)+1,hi(2)-1
          do n = 1,nc
             ss(0,i,j) = ss(0,i,j) + (betay(i,j,n)+betay(i,j+1,n))*f1(2) / alpha(i,j,n)
          end do
       end do
    end do

    do i = lo(1),hi(1)
       bclo = stencil_bc_type(mask( i,lo(2)),2,-1)
       bchi = stencil_bc_type(mask( i,hi(2)),2,+1)

       j = lo(2)
       if (bclo .eq. BC_INT) then
          do n = 1,nc
             ss(0,i,j) = ss(0,i,j) + (betay(i,j,n)+betay(i,j+1,n))*f1(2) / alpha(i,j,n)
          end do
       else
          blo = zero
          bhi = zero
          do n = 1,nc
             blo = blo + betay(i  ,j,n) / alpha(i,j,n) 
             bhi = bhi + betay(i,j+1,n) / alpha(i,j,n) 
          end do
          call stencil_bndry_aaa(order, ny, 2, -1, mask(i,j), &
               ss(0,i,j), ss(3,i,j), ss(4,i,j),ss(YBC,i,j), &
               blo, bhi, xa(2), xb(2), dh(2), bclo, bchi)
       end if

       if ( hi(2) > lo(2) ) then
          j = hi(2)
          if (bchi .eq. BC_INT) then
             do n = 1,nc
                ss(0,i,j) = ss(0,i,j) + (betay(i,j,n)+betay(i,j+1,n))*f1(2) / alpha(i,j,n)
             end do
          else
             blo = zero
             bhi = zero
             do n = 1,nc
                blo = blo + betay(i  ,j,n) / alpha(i,j,n) 
                bhi = bhi + betay(i,j+1,n) / alpha(i,j,n) 
             end do
             call stencil_bndry_aaa(order, ny, 2, 1, mask(i,j), &
                  ss(0,i,j), ss(3,i,j), ss(4,i,j), ss(YBC,i,j), &
                  blo, bhi, xa(2), xb(2), dh(2), bclo, bchi)
          end if
       end if
    end do

    do j = lo(2),hi(2)
       do i = lo(1),hi(1)
          ss(0,i,j) = ss(0,i,j) + alpha(i,j,0) 
       end do
    end do

  end subroutine s_simplen_2d_cc

 subroutine s_simplem_2d_cc(ss, alpha, ng_a, betax, betay, ng_b, dh, mask, lo, hi, xa, xb, order)

    integer           , intent(in   ) :: ng_a, ng_b, lo(:), hi(:), order
    integer           , intent(inout) :: mask(lo(1)  :,lo(2)  :)
    real (kind = dp_t), intent(  out) :: ss(0:, lo(1)  :,lo(2)  :)
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
 
    ss = zero

    ! Consider the operator  ( alpha - sum_n (beta0_n del dot beta_n grad) )
    ! Components alpha(i,j,   0) = alpha
    ! Components alpha(i,j,1:nc) = beta0_n
    ! Components betax(i,j,1:nc) = betax_n
    ! Components betay(i,j,1:nc) = betay_n

    ! ss(1,i,j) is the coefficient of phi(i+1,j  )
    ! ss(2,i,j) is the coefficient of phi(i-1,j  )
    ! ss(3,i,j) is the coefficient of phi(i  ,j+1)
    ! ss(4,i,j) is the coefficient of phi(i  ,j-1)
    ! ss(0,i,j) is the coefficient of phi(i  ,j  )
    
    do j = lo(2),hi(2)
       do i = lo(1),hi(1)
          ss(1,i,j) = ss(1,i,j) - (betax(i+1,j,1)+betax(i+1,j,2))*f1(1) 
          ss(2,i,j) = ss(2,i,j) - (betax(i  ,j,1)-betax(i,  j,2))*f1(1) 
          ss(3,i,j) = ss(3,i,j) - (betay(i,j+1,1)+betay(i,j+1,2))*f1(2) 
          ss(4,i,j) = ss(4,i,j) - (betay(i,j  ,1)-betay(i,j  ,2))*f1(2) 
       end do
    end do

    ! x derivatives

    do j = lo(2),hi(2)
       do i = lo(1)+1,hi(1)-1
            ss(0,i,j) = ss(0,i,j) - ss(1,i,j) - ss(2,i,j) 
       end do
    end do

    do j = lo(2),hi(2)
       bclo = stencil_bc_type(mask(lo(1),j),1,-1)
       bchi = stencil_bc_type(mask(hi(1),j),1,+1)
 
       i = lo(1)
       if (bclo .eq. BC_INT) then
          ss(0,i,j) = ss(0,i,j) - ss(1,i,j) - ss(2,i,j)
       else
          blo = -ss(2,i,j)/f1(1)  
          bhi = -ss(1,i,j)/f1(1)
          call stencil_bndry_aaa(order, nx, 1, -1, mask(i,j), &
               ss(0,i,j), ss(1,i,j), ss(2,i,j), ss(XBC,i,j), &
               blo, bhi, xa(1), xb(1), dh(1), bclo, bchi)
       end if

       if ( hi(1) > lo(1) ) then
          i = hi(1)
          if (bchi .eq. BC_INT) then
             ss(0,i,j) = ss(0,i,j) - ss(1,i,j) - ss(2,i,j) 
          else
             blo = -ss(2,i,j)/f1(1)  
             bhi = -ss(1,i,j)/f1(1)
             call stencil_bndry_aaa(order, nx, 1, 1, mask(i,j), &
                  ss(0,i,j), ss(1,i,j), ss(2,i,j), ss(XBC,i,j), &
                  blo, bhi, xa(1), xb(1), dh(1), bclo, bchi)
          end if
       end if
    end do

    ! y derivatives

    do i = lo(1),hi(1)
       do j = lo(2)+1,hi(2)-1
             ss(0,i,j) = ss(0,i,j) - ss(3,i,j) - ss(4,i,j) 
       end do
    end do

    do i = lo(1),hi(1)
       bclo = stencil_bc_type(mask( i,lo(2)),2,-1)
       bchi = stencil_bc_type(mask( i,hi(2)),2,+1)

       j = lo(2)
       if (bclo .eq. BC_INT) then
          ss(0,i,j) = ss(0,i,j) - ss(3,i,j) - ss(4,i,j) 
       else
          blo = -ss(4,i,j)/f1(2)  
          bhi = -ss(3,i,j)/f1(2)         
          call stencil_bndry_aaa(order, ny, 2, -1, mask(i,j), &
               ss(0,i,j), ss(3,i,j), ss(4,i,j),ss(YBC,i,j), &
               blo, bhi, xa(2), xb(2), dh(2), bclo, bchi)
       end if

       if ( hi(2) > lo(2) ) then
          j = hi(2)
          if (bchi .eq. BC_INT) then
             ss(0,i,j) = ss(0,i,j) - ss(3,i,j) - ss(4,i,j) 
          else
             blo = -ss(4,i,j)/f1(2)  
             bhi = -ss(3,i,j)/f1(2)
             call stencil_bndry_aaa(order, ny, 2, 1, mask(i,j), &
                  ss(0,i,j), ss(3,i,j), ss(4,i,j), ss(YBC,i,j), &
                  blo, bhi, xa(2), xb(2), dh(2), bclo, bchi)
          end if
       end if
    end do

    do j = lo(2),hi(2)
       do i = lo(1),hi(1)
          ss(0,i,j) = ss(0,i,j) + alpha(i,j,0) 
       end do
    end do

  end subroutine s_simplem_2d_cc

  subroutine s_simpleg_2d_cc(ss, alpha, ng_a, betax, betay, ng_b, dh, mask, lo, hi)

    integer           , intent(in   ) :: ng_a, ng_b, lo(:), hi(:)
    integer           , intent(inout) :: mask(lo(1)  :,lo(2)  :)
    real (kind = dp_t), intent(  out) :: ss(0:, lo(1)  :,lo(2)  :)
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
 
    ss = zero

    ! Consider the operator  ( alpha - sum_n (beta0_n del dot beta_n grad) )
    ! Components alpha(i,j,   0) = alpha
    ! Components alpha(i,j,1:nc) = beta0_n
    ! Components betax(i,j,1:nc) = betax_n
    ! Components betay(i,j,1:nc) = betay_n

    ! ss(1,i,j) is the coefficient of phi(i+1,j  )
    ! ss(2,i,j) is the coefficient of phi(i-1,j  )
    ! ss(3,i,j) is the coefficient of phi(i  ,j+1)
    ! ss(4,i,j) is the coefficient of phi(i  ,j-1)
    ! ss(0,i,j) is the coefficient of phi(i  ,j  )
    
    do j = lo(2),hi(2)
       do i = lo(1),hi(1)
          ss(1,i,j) = ss(1,i,j) - betax(i,j,1)
          ss(2,i,j) = ss(2,i,j) - betax(i,j,2)
          ss(3,i,j) = ss(3,i,j) - betay(i,j,1) 
          ss(4,i,j) = ss(4,i,j) - betay(i,j,2) 
       end do
    end do

    ! x derivatives

    do j = lo(2),hi(2)
       do i = lo(1)+1,hi(1)-1
          ss(0,i,j) = ss(0,i,j) - betax(i,j,3)
       end do
    end do

    do j = lo(2),hi(2)
       bclo = stencil_bc_type(mask(lo(1),j),1,-1)
       bchi = stencil_bc_type(mask(hi(1),j),1,+1)
 
       i = lo(1)
       if (bclo .eq. BC_INT) then
          ss(0,i,j) = ss(0,i,j) - betax(i,j,3)
       elseif (bclo .eq. BC_NEU) then
          ss(0,i,j) = ss(0,i,j) - betax(i,j,3) - betax(i,j,2)
          ss(2,i,j) = zero
          ss(XBC,i,j) = zero
       elseif (bclo .eq. BC_DIR) then
          ss(0,i,j) = ss(0,i,j) - betax(i,j,3)
          ss(2,i,j) = zero
          ss(XBC,i,j) = zero
       end if

       if ( hi(1) > lo(1) ) then
          i = hi(1)
          if (bchi .eq. BC_INT) then
             ss(0,i,j) = ss(0,i,j) - betax(i,j,3)
          elseif (bchi .eq. BC_NEU) then
             ss(0,i,j) = ss(0,i,j) - betax(i,j,3) - betax(i,j,1)
             ss(1,i,j) = zero
             ss(XBC,i,j) = zero
          elseif (bchi .eq. BC_DIR) then
             ss(0,i,j) = ss(0,i,j) - betax(i,j,3)
             ss(1,i,j) = zero
             ss(XBC,i,j) = zero
          end if
       end if
    end do

    ! y derivatives
    do i = lo(1),hi(1)
       do j = lo(2)+1,hi(2)-1
          ss(0,i,j) = ss(0,i,j) - betay(i,j,3)
       end do
    end do

    do i = lo(1),hi(1)
       bclo = stencil_bc_type(mask( i,lo(2)),2,-1)
       bchi = stencil_bc_type(mask( i,hi(2)),2,+1)

       j = lo(2)
       if (bclo .eq. BC_INT) then
          ss(0,i,j) = ss(0,i,j) - betay(i,j,3)
       elseif (bclo .eq. BC_NEU) then
          ss(0,i,j)   = ss(0,i,j) - betay(i,j,3) - betay(i,j,2)
          ss(4,i,j)   = zero
          ss(YBC,i,j) = zero
       elseif (bclo .eq. BC_DIR) then
          ss(0,i,j) = ss(0,i,j) - betay(i,j,3) 
          ss(4,i,j) = zero
          ss(YBC,i,j) = zero
       end if

       if ( hi(2) > lo(2) ) then
          j = hi(2)
          if (bchi .eq. BC_INT) then
             ss(0,i,j) = ss(0,i,j) - betay(i,j,3)
          elseif (bchi .eq. BC_NEU) then
             ss(0,i,j) = ss(0,i,j) - betay(i,j,3) - betay(i,j,1)
             ss(3,i,j) = zero
             ss(YBC,i,j) = zero
          elseif (bchi .eq. BC_DIR) then
             ss(0,i,j) = ss(0,i,j) - betay(i,j,3) 
             ss(3,i,j) = zero
             ss(YBC,i,j) = zero
          end if
       end if
    end do

    do j = lo(2),hi(2)
       do i = lo(1),hi(1)
          ss(0,i,j) = ss(0,i,j) + alpha(i,j,0) 
       end do
    end do

  end subroutine s_simpleg_2d_cc

  subroutine s_simple_3d_cc(ss, alpha, ng_a, betax, betay, betaz, ng_b, dh, mask, lo, hi, xa, xb, order)


    integer           , intent(in   ) :: ng_a, ng_b, lo(:), hi(:), order
    integer           , intent(inout) :: mask(lo(1)  :,lo(2)  :,lo(3)  :)
    real (kind = dp_t), intent(  out) ::   ss(0:, lo(1)  :,lo(2)  :,lo(3)  :)
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

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             ss(0,i,j,k)   = ZERO
             ss(1,i,j,k)   = -betax(i+1, j  , k  )*f1(1)
             ss(2,i,j,k)   = -betax(i  , j  , k  )*f1(1)
             ss(3,i,j,k)   = -betay(i  , j+1, k  )*f1(2)
             ss(4,i,j,k)   = -betay(i  , j  , k  )*f1(2)
             ss(5,i,j,k)   = -betaz(i  , j  , k+1)*f1(3)
             ss(6,i,j,k)   = -betaz(i  , j  , k  )*f1(3)
             ss(XBC,i,j,k) = ZERO
             ss(YBC,i,j,k) = ZERO
             ss(ZBC,i,j,k) = ZERO
          end do
       end do
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             mask(i,j,k) = ibclr(mask(i,j,k), BC_BIT(BC_GEOM,1,-1))
             mask(i,j,k) = ibclr(mask(i,j,k), BC_BIT(BC_GEOM,1,+1))
             mask(i,j,k) = ibclr(mask(i,j,k), BC_BIT(BC_GEOM,2,-1))
             mask(i,j,k) = ibclr(mask(i,j,k), BC_BIT(BC_GEOM,2,+1))
             mask(i,j,k) = ibclr(mask(i,j,k), BC_BIT(BC_GEOM,3,-1))
             mask(i,j,k) = ibclr(mask(i,j,k), BC_BIT(BC_GEOM,3,+1))
          end do
       end do
    end do
    !$OMP END PARALLEL DO

    lnx = nx; lny = ny; lnz = nz; lorder = order

    ! x derivatives

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1)+1,hi(1)-1
             ss(0,i,j,k) = ss(0,i,j,k) + (betax(i,j,k)+betax(i+1,j,k))*f1(1)
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
             ss(0,i,j,k) = ss(0,i,j,k) + (betax(i,j,k)+betax(i+1,j,k))*f1(1)
          else
             call stencil_bndry_aaa(lorder, lnx, 1, -1, mask(i,j,k), &
                  ss(0,i,j,k), ss(1,i,j,k), ss(2,i,j,k), ss(XBC,i,j,k), &
                  betax(i,j,k), betax(i+1,j,k), xa(1), xb(1), dh(1), bclo, bchi)
          end if

          if ( hi(1) > lo(1) ) then
             i = hi(1)
             if (bchi .eq. BC_INT) then
                ss(0,i,j,k) = ss(0,i,j,k) + (betax(i,j,k)+betax(i+1,j,k))*f1(1)
             else
                call stencil_bndry_aaa(lorder, lnx, 1, 1, mask(i,j,k), &
                     ss(0,i,j,k), ss(1,i,j,k), ss(2,i,j,k), ss(XBC,i,j,k), &
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
             ss(0,i,j,k) = ss(0,i,j,k) + (betay(i,j,k)+betay(i,j+1,k))*f1(2)
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
             ss(0,i,j,k) = ss(0,i,j,k) + (betay(i,j,k)+betay(i,j+1,k))*f1(2)
          else
             call stencil_bndry_aaa(lorder, lny, 2, -1, mask(i,j,k), &
                  ss(0,i,j,k), ss(3,i,j,k), ss(4,i,j,k),ss(YBC,i,j,k), &
                  betay(i,j,k), betay(i,j+1,k), xa(2), xb(2), dh(2), bclo, bchi)
          end if
          if ( hi(2) > lo(2) ) then
             j = hi(2)
             if (bchi .eq. BC_INT) then
                ss(0,i,j,k) = ss(0,i,j,k) + (betay(i,j,k)+betay(i,j+1,k))*f1(2)
             else
                call stencil_bndry_aaa(lorder, lny, 2, 1, mask(i,j,k), &
                     ss(0,i,j,k), ss(3,i,j,k), ss(4,i,j,k), ss(YBC,i,j,k), &
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
             ss(0,i,j,k) = ss(0,i,j,k) + (betaz(i,j,k)+betaz(i,j,k+1))*f1(3)
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
             ss(0,i,j,k) = ss(0,i,j,k) + (betaz(i,j,k)+betaz(i,j,k+1))*f1(3)
          else
             call stencil_bndry_aaa(lorder, lnz, 3, -1, mask(i,j,k), &
                  ss(0,i,j,k), ss(5,i,j,k), ss(6,i,j,k),ss(ZBC,i,j,k), &
                  betaz(i,j,k), betaz(i,j,k+1), xa(3), xb(3), dh(3), bclo, bchi)
          end if
          if ( hi(3) > lo(3) ) then
             k = hi(3)
             if (bchi .eq. BC_INT) then
                ss(0,i,j,k) = ss(0,i,j,k) + (betaz(i,j,k)+betaz(i,j,k+1))*f1(3)
             else
                call stencil_bndry_aaa(lorder, lnz, 3, 1, mask(i,j,k), &
                     ss(0,i,j,k), ss(5,i,j,k), ss(6,i,j,k), ss(ZBC,i,j,k), &
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
             ss(0,i,j,k) = ss(0,i,j,k) + alpha(i,j,k)
          end do
       end do
    end do
    !$OMP END PARALLEL DO

  end subroutine s_simple_3d_cc

  subroutine s_simpleg_3d_cc(ss, alpha, ng_a, betax, betay, betaz, ng_b, dh, mask, lo, hi)

    integer           , intent(in   ) :: ng_a, ng_b, lo(:), hi(:)
    integer           , intent(inout) :: mask(lo(1)  :,lo(2)  :,lo(3)  :)
    real (kind = dp_t), intent(  out) ::   ss(0:, lo(1)  :,lo(2)  :,lo(3)  :)
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

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             ss(0,i,j,k)   = ZERO
             ss(1,i,j,k)   = -betax(i,j,k,1)
             ss(2,i,j,k)   = -betax(i,j,k,2)
             ss(3,i,j,k)   = -betay(i,j,k,1)
             ss(4,i,j,k)   = -betay(i,j,k,2)
             ss(5,i,j,k)   = -betaz(i,j,k,1)
             ss(6,i,j,k)   = -betaz(i,j,k,2)
             ss(XBC,i,j,k) = ZERO
             ss(YBC,i,j,k) = ZERO
             ss(ZBC,i,j,k) = ZERO
          end do
       end do
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             mask(i,j,k) = ibclr(mask(i,j,k), BC_BIT(BC_GEOM,1,-1))
             mask(i,j,k) = ibclr(mask(i,j,k), BC_BIT(BC_GEOM,1,+1))
             mask(i,j,k) = ibclr(mask(i,j,k), BC_BIT(BC_GEOM,2,-1))
             mask(i,j,k) = ibclr(mask(i,j,k), BC_BIT(BC_GEOM,2,+1))
             mask(i,j,k) = ibclr(mask(i,j,k), BC_BIT(BC_GEOM,3,-1))
             mask(i,j,k) = ibclr(mask(i,j,k), BC_BIT(BC_GEOM,3,+1))
          end do
       end do
    end do
    !$OMP END PARALLEL DO

    ! x derivatives

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1)+1,hi(1)-1
             ss(0,i,j,k) = ss(0,i,j,k) - betax(i,j,k,3)
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
             ss(0,i,j,k) = ss(0,i,j,k) - betax(i,j,k,3)
          elseif (bclo .eq. BC_NEU) then
             ss(0,i,j,k) = ss(0,i,j,k) - betax(i,j,k,3) - betax(i,j,k,2)
             ss(2,i,j,k) = zero
             ss(XBC,i,j,k) = zero
          elseif (bclo .eq. BC_DIR) then
             ss(0,i,j,k) = ss(0,i,j,k) - betax(i,j,k,3)
             ss(2,i,j,k) = zero
             ss(XBC,i,j,k) = zero
          end if

          if ( hi(1) > lo(1) ) then
             i = hi(1)
             if (bchi .eq. BC_INT) then
                ss(0,i,j,k) = ss(0,i,j,k) - betax(i,j,k,3)
             elseif (bchi .eq. BC_NEU) then
                ss(0,i,j,k) = ss(0,i,j,k) - betax(i,j,k,3) - betax(i,j,k,1)
                ss(1,i,j,k) = zero
                ss(XBC,i,j,k) = zero
             elseif (bchi .eq. BC_DIR) then
                ss(0,i,j,k) = ss(0,i,j,k) - betax(i,j,k,3)
                ss(1,i,j,k) = zero
                ss(XBC,i,j,k) = zero
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
             ss(0,i,j,k) = ss(0,i,j,k) - betay(i,j,k,3)
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
             ss(0,i,j,k) = ss(0,i,j,k) - betay(i,j,k,3)
          elseif (bclo .eq. BC_NEU) then
             ss(0,i,j,k)   = ss(0,i,j,k) - betay(i,j,k,3) - betay(i,j,k,2)
             ss(4,i,j,k)   = zero
             ss(YBC,i,j,k) = zero
          elseif (bclo .eq. BC_DIR) then
             ss(0,i,j,k)   = ss(0,i,j,k) - betay(i,j,k,3) 
             ss(4,i,j,k)   = zero
             ss(YBC,i,j,k) = zero
          end if

          if ( hi(2) > lo(2) ) then
             j = hi(2)
             if (bchi .eq. BC_INT) then
                ss(0,i,j,k) = ss(0,i,j,k) - betay(i,j,k,3)
             elseif (bchi .eq. BC_NEU) then
                ss(0,i,j,k)   = ss(0,i,j,k) - betay(i,j,k,3) - betay(i,j,k,1)
                ss(3,i,j,k)   = zero
                ss(YBC,i,j,k) = zero
             elseif (bchi .eq. BC_DIR) then
                ss(0,i,j,k)   = ss(0,i,j,k) - betay(i,j,k,3) 
                ss(3,i,j,k)   = zero
                ss(YBC,i,j,k) = zero
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
             ss(0,i,j,k) = ss(0,i,j,k) - betaz(i,j,k,3)
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
             ss(0,i,j,k) = ss(0,i,j,k) - betaz(i,j,k,3)
          elseif (bclo .eq. BC_NEU) then
             ss(0,i,j,k)   = ss(0,i,j,k) - betaz(i,j,k,3) - betaz(i,j,k,2)
             ss(6,i,j,k)   = zero
             ss(ZBC,i,j,k) = zero
          elseif (bclo .eq. BC_DIR) then
             ss(0,i,j,k)   = ss(0,i,j,k) - betaz(i,j,k,3) 
             ss(6,i,j,k)   = zero
             ss(ZBC,i,j,k) = zero
          end if

          if ( hi(3) > lo(3) ) then
             k = hi(3)
             if (bchi .eq. BC_INT) then
                ss(0,i,j,k) = ss(0,i,j,k) - betaz(i,j,k,3)
             elseif (bchi .eq. BC_NEU) then
                ss(0,i,j,k)   = ss(0,i,j,k) - betaz(i,j,k,3) - betaz(i,j,k,1)
                ss(5,i,j,k)   = zero
                ss(ZBC,i,j,k) = zero
             elseif (bchi .eq. BC_DIR) then
                ss(0,i,j,k)   = ss(0,i,j,k) - betaz(i,j,k,3) 
                ss(5,i,j,k)   = zero
                ss(ZBC,i,j,k) = zero
             end if
          end if
       end do
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             ss(0,i,j,k) = ss(0,i,j,k) + alpha(i,j,k,0)
          end do
       end do
    end do
    !$OMP END PARALLEL DO

  end subroutine s_simpleg_3d_cc

  subroutine s_minion_cross_fill_2d(ss, alpha, ng_a, betax, betay, ng_b, dh, mask, lo, hi)

    integer           , intent(in   ) :: ng_a,ng_b
    integer           , intent(in   ) :: lo(:), hi(:)
    integer           , intent(inout) :: mask(lo(1):,lo(2):)
    real (kind = dp_t), intent(  out) :: ss(0:,lo(1):,lo(2):)
    real (kind = dp_t), intent(in   ) :: alpha(lo(1)-ng_a:,lo(2)-ng_a:)
    real (kind = dp_t), intent(in   ) :: betax(lo(1)-ng_b:,lo(2)-ng_b:)
    real (kind = dp_t), intent(in   ) :: betay(lo(1)-ng_b:,lo(2)-ng_b:)
    real (kind = dp_t), intent(in   ) :: dh(:)

    integer            :: i, j
    real (kind = dp_t) :: fac

    fac = (ONE / ( TWELVE * dh(1)**2))

    mask = ibclr(mask, BC_BIT(BC_GEOM,1,-1))
    mask = ibclr(mask, BC_BIT(BC_GEOM,1,+1))
    mask = ibclr(mask, BC_BIT(BC_GEOM,2,-1))
    mask = ibclr(mask, BC_BIT(BC_GEOM,2,+1))
    !
    ! We only include the beta's here to get the viscous coefficients in here for now.
    ! The projection has beta == 1.
    !
    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          ss(1,i,j) =      ONE * betax(i  ,j) * fac
          ss(2,i,j) = -SIXTEEN * betax(i  ,j) * fac
          ss(3,i,j) = -SIXTEEN * betax(i+1,j) * fac
          ss(4,i,j) =      ONE * betax(i+1,j) * fac
          ss(5,i,j) =      ONE * betay(i,j  ) * fac
          ss(6,i,j) = -SIXTEEN * betay(i,j  ) * fac
          ss(7,i,j) = -SIXTEEN * betay(i,j+1) * fac
          ss(8,i,j) =      ONE * betay(i,j+1) * fac

          ss(0,i,j) = -( ss(1,i,j) + ss(2,i,j) + ss(3,i,j) + ss(4,i,j)   &
                        +ss(5,i,j) + ss(6,i,j) + ss(7,i,j) + ss(8,i,j) ) + alpha(i,j)
       end do
    end do

  end subroutine s_minion_cross_fill_2d

  ! subroutine s_minion_full_old_2d(ss, beta, ng_b, dh, mask, lo, hi)

  !   integer           , intent(in   ) :: ng_b
  !   integer           , intent(in   ) :: lo(:), hi(:)
  !   integer           , intent(inout) :: mask(:,:)
  !   real (kind = dp_t), intent(  out) :: ss(0:,:,:)
  !   real (kind = dp_t), intent(inout) :: beta(1-ng_b:,1-ng_b:,0:)
  !   real (kind = dp_t), intent(in   ) :: dh(:)

  !   integer            :: i, j, nx, ny
  !   real (kind = dp_t) :: fac

  !   nx = hi(1)-lo(1)+1
  !   ny = hi(2)-lo(2)+1

  !   mask = ibclr(mask, BC_BIT(BC_GEOM,1,-1))
  !   mask = ibclr(mask, BC_BIT(BC_GEOM,1,+1))
  !   mask = ibclr(mask, BC_BIT(BC_GEOM,2,-1))
  !   mask = ibclr(mask, BC_BIT(BC_GEOM,2,+1))

  !   ! First use the betax coefficients
  !   do j = 1, ny
  !      do i = 1, nx

  !         ss(12,i,j) = 27648._dp_t*beta(i+1,j,1) + 414720._dp_t * beta(i  ,j,1)
  !         ss(13,i,j) = 27648._dp_t*beta(i  ,j,1) + 414720._dp_t * beta(i+1,j,1)

  !         ss(11,i,j) = -27648._dp_t * beta(i  ,j,1)
  !         ss(14,i,j) = -27648._dp_t * beta(i+1,j,1)

  !         ss( 8,i,j) = -2550._dp_t * beta(i,j+2,1)  - 2550._dp_t * beta(i+1,j+2,1) &
  !                     +17340._dp_t * beta(i,j+1,1) + 17340._dp_t * beta(i+1,j+1,1) &
  !                     -17340._dp_t * beta(i,j-1,1) - 17340._dp_t * beta(i+1,j-1,1) &
  !                     + 2550._dp_t * beta(i,j-2,1)  + 2550._dp_t * beta(i+1,j-2,1)

  !         ss(17,i,j) = -2550._dp_t * beta(i,j-2,1)  - 2550._dp_t * beta(i+1,j-2,1) &
  !                     +17340._dp_t * beta(i,j-1,1) + 17340._dp_t * beta(i+1,j-1,1) &
  !                     -17340._dp_t * beta(i,j+1,1) - 17340._dp_t * beta(i+1,j+1,1) &
  !                     + 2550._dp_t * beta(i,j+2,1)  + 2550._dp_t * beta(i+1,j+2,1)

  !         ss( 7,i,j) =  170._dp_t * beta(i+1,j+2,1) +  2550._dp_t * beta(i,j+2,1) &
  !                     -1156._dp_t * beta(i+1,j+1,1) - 17340._dp_t * beta(i,j+1,1) &
  !                     +1156._dp_t * beta(i+1,j-1,1) + 17340._dp_t * beta(i,j-1,1) &
  !                     - 170._dp_t * beta(i+1,j-2,1) -  2550._dp_t * beta(i,j-2,1) 

  !         ss(16,i,j) =  170._dp_t * beta(i+1,j-2,1) +  2550._dp_t * beta(i,j-2,1) &
  !                     -1156._dp_t * beta(i+1,j-1,1) - 17340._dp_t * beta(i,j-1,1) &
  !                     +1156._dp_t * beta(i+1,j+1,1) + 17340._dp_t * beta(i,j+1,1) &
  !                     - 170._dp_t * beta(i+1,j+2,1) -  2550._dp_t * beta(i,j+2,1) 

  !         ss( 9,i,j) =  170._dp_t * beta(i,j+2,1) +  2550._dp_t * beta(i+1,j+2,1) &
  !                     -1156._dp_t * beta(i,j+1,1) - 17340._dp_t * beta(i+1,j+1,1) &
  !                     +1156._dp_t * beta(i,j-1,1) + 17340._dp_t * beta(i+1,j-1,1) &
  !                     - 170._dp_t * beta(i,j-2,1) -  2550._dp_t * beta(i+1,j-2,1) 

  !         ss(18,i,j) =  170._dp_t * beta(i,j-2,1) +  2550._dp_t * beta(i+1,j-2,1) &
  !                     -1156._dp_t * beta(i,j-1,1) - 17340._dp_t * beta(i+1,j-1,1) &
  !                     +1156._dp_t * beta(i,j+1,1) + 17340._dp_t * beta(i+1,j+1,1) &
  !                     - 170._dp_t * beta(i,j+2,1) -  2550._dp_t * beta(i+1,j+2,1) 

  !         ss( 6,i,j) = -170._dp_t * beta(i,j+2,1) +  1156._dp_t * beta(i,j+1,1) &
  !                      +170._dp_t * beta(i,j-2,1) -  1156._dp_t * beta(i,j-1,1)
  !         ss(15,i,j) = -170._dp_t * beta(i,j-2,1) +  1156._dp_t * beta(i,j-1,1) &
  !                      +170._dp_t * beta(i,j+2,1) -  1156._dp_t * beta(i,j+1,1)
  !         ss(10,i,j) = -170._dp_t * beta(i+1,j+2,1) +  1156._dp_t * beta(i+1,j+1,1) &
  !                      +170._dp_t * beta(i+1,j-2,1) -  1156._dp_t * beta(i+1,j-1,1)
  !         ss(19,i,j) = -170._dp_t * beta(i+1,j-2,1) +  1156._dp_t * beta(i+1,j-1,1) &
  !                      +170._dp_t * beta(i+1,j+2,1) -  1156._dp_t * beta(i+1,j+1,1)

  !         ss( 3,i,j) =   375._dp_t * beta(i,j+2,1) +  375._dp_t * beta(i+1,j+2,1) &
  !                     - 2550._dp_t * beta(i,j+1,1) - 2550._dp_t * beta(i+1,j+1,1) &
  !                     + 2550._dp_t * beta(i,j-1,1) + 2550._dp_t * beta(i+1,j-1,1) &
  !                     -  375._dp_t * beta(i,j-2,1) -  375._dp_t * beta(i+1,j-2,1)
  !         ss(22,i,j) =  375._dp_t * beta(i,j-2,1) +  375._dp_t * beta(i+1,j-2,1) &
  !                     -2550._dp_t * beta(i,j-1,1) - 2550._dp_t * beta(i+1,j-1,1) &
  !                     +2550._dp_t * beta(i,j+1,1) + 2550._dp_t * beta(i+1,j+1,1) &
  !                     - 375._dp_t * beta(i,j+2,1) -  375._dp_t * beta(i+1,j+2,1)

  !         ss( 2,i,j) = - 25._dp_t * beta(i+1,j+2,1) -  375._dp_t * beta(i,j+2,1) &
  !                      +170._dp_t * beta(i+1,j+1,1) + 2550._dp_t * beta(i,j+1,1) &
  !                      -170._dp_t * beta(i+1,j-1,1) - 2550._dp_t * beta(i,j-1,1) &
  !                      + 25._dp_t * beta(i+1,j-2,1) +  375._dp_t * beta(i,j-2,1)
  !         ss(21,i,j) = - 25._dp_t * beta(i+1,j-2,1) -  375._dp_t * beta(i,j-2,1) &
  !                      +170._dp_t * beta(i+1,j-1,1) + 2550._dp_t * beta(i,j-1,1) &
  !                      -170._dp_t * beta(i+1,j+1,1) - 2550._dp_t * beta(i,j+1,1) &
  !                      + 25._dp_t * beta(i+1,j+2,1) +  375._dp_t * beta(i,j+2,1)
  !         ss( 4,i,j) = - 25._dp_t * beta(i,j+2,1) -  375._dp_t * beta(i+1,j+2,1) &
  !                      +170._dp_t * beta(i,j+1,1) + 2550._dp_t * beta(i+1,j+1,1) &
  !                      -170._dp_t * beta(i,j-1,1) - 2550._dp_t * beta(i+1,j-1,1) &
  !                      + 25._dp_t * beta(i,j-2,1) +  375._dp_t * beta(i+1,j-2,1)
  !         ss(23,i,j) = - 25._dp_t * beta(i,j-2,1) -  375._dp_t * beta(i+1,j-2,1) &
  !                      +170._dp_t * beta(i,j-1,1) + 2550._dp_t * beta(i+1,j-1,1) &
  !                      -170._dp_t * beta(i,j+1,1) - 2550._dp_t * beta(i+1,j+1,1) &
  !                      + 25._dp_t * beta(i,j+2,1) +  375._dp_t * beta(i+1,j+2,1)

  !         ss( 1,i,j) =   25._dp_t * beta(i,j+2,1) -  170._dp_t * beta(i,j+1,1) &
  !                       -25._dp_t * beta(i,j-2,1) +  170._dp_t * beta(i,j-1,1)
  !         ss( 5,i,j) =   25._dp_t * beta(i+1,j+2,1) -  170._dp_t * beta(i+1,j+1,1) &
  !                       -25._dp_t * beta(i+1,j-2,1) +  170._dp_t * beta(i+1,j-1,1)
  !         ss(20,i,j) =   25._dp_t * beta(i,j-2,1) -  170._dp_t * beta(i,j-1,1) &
  !                       -25._dp_t * beta(i,j+2,1) +  170._dp_t * beta(i,j+1,1)
  !         ss(24,i,j) =   25._dp_t * beta(i+1,j-2,1) -  170._dp_t * beta(i+1,j-1,1) &
  !                       -25._dp_t * beta(i+1,j+2,1) +  170._dp_t * beta(i+1,j+1,1)

  !         ss( 0,i,j) = -414720._dp_t * (beta(i,j,1) + beta(i+1,j,1))

  !      end do
  !   end do

  !   fac = -ONE / (TWELVE**2 * 48._dp_t**2 * dh(1)**2)

  !   ! Then use the betay coefficients
  !   do j = 1, ny
  !      do i = 1, nx

  !         ss( 8,i,j) = ( ss( 8,i,j) + 27648._dp_t*beta(i,j+1,2) + 414720._dp_t * beta(i,j  ,2) ) * fac
  !         ss(17,i,j) = ( ss(17,i,j) + 27648._dp_t*beta(i,j  ,2) + 414720._dp_t * beta(i,j+1,2) ) * fac

  !         ss( 3,i,j) = ( ss( 3,i,j) - 27648._dp_t * beta(i,j  ,2) ) * fac
  !         ss(22,i,j) = ( ss(22,i,j) - 27648._dp_t * beta(i,j+1,2) ) * fac

  !         ss(12,i,j) = ( ss(12,i,j) & 
  !                      -2550._dp_t * beta(i+2,j,2)  - 2550._dp_t * beta(i+2,j+1,2) &
  !                     +17340._dp_t * beta(i+1,j,2) + 17340._dp_t * beta(i+1,j+1,2) &
  !                     -17340._dp_t * beta(i-1,j,2) - 17340._dp_t * beta(i-1,j+1,2) &
  !                     + 2550._dp_t * beta(i-2,j,2)  + 2550._dp_t * beta(i-2,j+1,2) ) * fac

  !         ss(13,i,j) = ( ss(13,i,j) & 
  !                      -2550._dp_t * beta(i-2,j,2)  - 2550._dp_t * beta(i-2,j+1,2) &
  !                     +17340._dp_t * beta(i-1,j,2) + 17340._dp_t * beta(i-1,j+1,2) &
  !                     -17340._dp_t * beta(i+1,j,2) - 17340._dp_t * beta(i+1,j+1,2) &
  !                     + 2550._dp_t * beta(i+2,j,2)  + 2550._dp_t * beta(i+2,j+1,2) ) * fac

  !         ss( 7,i,j) = ( ss( 7,i,j) &
  !                     + 170._dp_t * beta(i+2,j+1,2) +  2550._dp_t * beta(i+2,j  ,2) &
  !                     -1156._dp_t * beta(i+1,j+1,2) - 17340._dp_t * beta(i+1,j  ,2) &
  !                     +1156._dp_t * beta(i-1,j+1,2) + 17340._dp_t * beta(i-1,j  ,2) &
  !                     - 170._dp_t * beta(i-2,j+1,2) -  2550._dp_t * beta(i-2,j  ,2) ) * fac

  !         ss(16,i,j) = ( ss(16,i,j) &  
  !                     + 170._dp_t * beta(i+2,j  ,2) +  2550._dp_t * beta(i+2,j+1,2) &
  !                     -1156._dp_t * beta(i+1,j  ,2) - 17340._dp_t * beta(i+1,j+1,2) &
  !                     +1156._dp_t * beta(i-1,j  ,2) + 17340._dp_t * beta(i-1,j+1,2) &
  !                     - 170._dp_t * beta(i-2,j  ,2) -  2550._dp_t * beta(i-2,j+1,2) ) * fac

  !         ss( 9,i,j) = ( ss( 9,i,j) &  
  !                    +  170._dp_t * beta(i-2,j+1,2) +  2550._dp_t * beta(i-2,j  ,2) &
  !                     -1156._dp_t * beta(i-1,j+1,2) - 17340._dp_t * beta(i-1,j  ,2) &
  !                     +1156._dp_t * beta(i+1,j+1,2) + 17340._dp_t * beta(i+1,j  ,2) &
  !                     - 170._dp_t * beta(i+2,j+1,2) -  2550._dp_t * beta(i+2,j  ,2) ) * fac

  !         ss(18,i,j) = ( ss(18,i,j) &  
  !                    +  170._dp_t * beta(i-2,j  ,2) +  2550._dp_t * beta(i-2,j+1,2) &
  !                     -1156._dp_t * beta(i-1,j  ,2) - 17340._dp_t * beta(i-1,j+1,2) &
  !                     +1156._dp_t * beta(i+1,j  ,2) + 17340._dp_t * beta(i+1,j+1,2) &
  !                     - 170._dp_t * beta(i+2,j  ,2) -  2550._dp_t * beta(i+2,j+1,2) ) * fac

  !         ss( 2,i,j) = ( ss( 2,i,j) &
  !                      -170._dp_t * beta(i+2,j,2) +  1156._dp_t * beta(i+1,j,2) &
  !                      +170._dp_t * beta(i-2,j,2) -  1156._dp_t * beta(i-1,j,2) ) * fac

  !         ss(21,i,j) = ( ss(21,i,j) &
  !                      -170._dp_t * beta(i+2,j+1,2) +  1156._dp_t * beta(i+1,j+1,2) &
  !                      +170._dp_t * beta(i-2,j+1,2) -  1156._dp_t * beta(i-1,j+1,2) ) * fac

  !         ss( 4,i,j) = ( ss( 4,i,j) &
  !                      -170._dp_t * beta(i-2,j,2) +  1156._dp_t * beta(i-1,j,2) &
  !                      +170._dp_t * beta(i+2,j,2) -  1156._dp_t * beta(i+1,j,2) ) * fac

  !         ss(23,i,j) = ( ss(23,i,j) &
  !                      -170._dp_t * beta(i-2,j+1,2) +  1156._dp_t * beta(i-1,j+1,2) &
  !                      +170._dp_t * beta(i+2,j+1,2) -  1156._dp_t * beta(i+1,j+1,2) ) * fac

  !         ss(11,i,j) = ( ss(11,i,j) &
  !                     +  375._dp_t * beta(i+2,j,2) +  375._dp_t * beta(i+2,j+1,2) &
  !                     - 2550._dp_t * beta(i+1,j,2) - 2550._dp_t * beta(i+1,j+1,2) &
  !                     + 2550._dp_t * beta(i-1,j,2) + 2550._dp_t * beta(i-1,j+1,2) &
  !                     -  375._dp_t * beta(i-2,j,2) -  375._dp_t * beta(i-2,j+1,2) ) * fac

  !         ss(14,i,j) = ( ss(14,i,j) &
  !                    +  375._dp_t * beta(i-2,j,2) +  375._dp_t * beta(i-2,j+1,2) &
  !                     -2550._dp_t * beta(i-1,j,2) - 2550._dp_t * beta(i-1,j+1,2) &
  !                     +2550._dp_t * beta(i+1,j,2) + 2550._dp_t * beta(i+1,j+1,2) &
  !                     - 375._dp_t * beta(i+2,j,2) -  375._dp_t * beta(i+2,j+1,2) ) * fac

  !         ss( 6,i,j) = ( ss( 6,i,j) &
  !                      - 25._dp_t * beta(i+2,j+1,2) -  375._dp_t * beta(i+2,j,2) &
  !                      +170._dp_t * beta(i+1,j+1,2) + 2550._dp_t * beta(i+1,j,2) &
  !                      -170._dp_t * beta(i-1,j+1,2) - 2550._dp_t * beta(i-1,j,2) &
  !                      + 25._dp_t * beta(i-2,j+1,2) +  375._dp_t * beta(i-2,j,2) ) * fac

  !         ss(15,i,j) = ( ss(15,i,j) &
  !                      - 25._dp_t * beta(i+2,j,2) -  375._dp_t * beta(i+2,j+1,2) &
  !                      +170._dp_t * beta(i+1,j,2) + 2550._dp_t * beta(i+1,j+1,2) &
  !                      -170._dp_t * beta(i-1,j,2) - 2550._dp_t * beta(i-1,j+1,2) &
  !                      + 25._dp_t * beta(i-2,j,2) +  375._dp_t * beta(i-2,j+1,2) ) * fac

  !         ss(10,i,j) = ( ss(10,i,j) &
  !                      - 25._dp_t * beta(i-2,j+1,2) -  375._dp_t * beta(i-2,j,2) &
  !                      +170._dp_t * beta(i-1,j+1,2) + 2550._dp_t * beta(i-1,j,2) &
  !                      -170._dp_t * beta(i+1,j+1,2) - 2550._dp_t * beta(i+1,j,2) &
  !                      + 25._dp_t * beta(i+2,j+1,2) +  375._dp_t * beta(i+2,j,2) ) * fac

  !         ss(19,i,j) = ( ss(19,i,j) &
  !                      - 25._dp_t * beta(i-2,j,2) -  375._dp_t * beta(i-2,j+1,2) &
  !                      +170._dp_t * beta(i-1,j,2) + 2550._dp_t * beta(i-1,j+1,2) &
  !                      -170._dp_t * beta(i+1,j,2) - 2550._dp_t * beta(i+1,j+1,2) &
  !                      + 25._dp_t * beta(i+2,j,2) +  375._dp_t * beta(i+2,j+1,2) ) * fac

  !         ss( 1,i,j) = ( ss( 1,i,j) &
  !                      + 25._dp_t * beta(i+2,j,2) -  170._dp_t * beta(i+1,j,2) &
  !                       -25._dp_t * beta(i-2,j,2) +  170._dp_t * beta(i-1,j,2) ) * fac
  !         ss( 5,i,j) = ( ss( 5,i,j) &
  !                      + 25._dp_t * beta(i-2,j,2) -  170._dp_t * beta(i-1,j,2) &
  !                       -25._dp_t * beta(i+2,j,2) +  170._dp_t * beta(i+1,j,2) ) * fac
  !         ss(20,i,j) = ( ss(20,i,j) &
  !                      + 25._dp_t * beta(i+2,j+1,2) -  170._dp_t * beta(i+1,j+1,2) &
  !                       -25._dp_t * beta(i-2,j+1,2) +  170._dp_t * beta(i-1,j+1,2) ) * fac
  !         ss(24,i,j) = ( ss(24,i,j) &
  !                      + 25._dp_t * beta(i-2,j+1,2) -  170._dp_t * beta(i-1,j+1,2) &
  !                       -25._dp_t * beta(i+2,j+1,2) +  170._dp_t * beta(i+1,j+1,2) ) * fac

  !         ss( 0,i,j) = ( ss(0,i,j) -414720._dp_t * ( beta(i,j,2) + beta(i,j+1,2) ) ) * fac + beta(i,j,0)

  !      end do
  !   end do

  ! end subroutine s_minion_full_old_2d

  subroutine s_minion_full_fill_2d(ss, alpha, ng_a, betax, betay, ng_b, dh, mask, lo, hi)

    integer           , intent(in   ) :: ng_a,ng_b
    integer           , intent(in   ) :: lo(:), hi(:)
    integer           , intent(inout) :: mask(lo(1):,lo(2):)
    real (kind = dp_t), intent(  out) :: ss(0:,lo(1):     ,lo(2):     )
    real (kind = dp_t), intent(in   ) :: alpha(lo(1)-ng_a:,lo(2)-ng_a:)
    real (kind = dp_t), intent(in   ) :: betax(lo(1)-ng_b:,lo(2)-ng_b:)
    real (kind = dp_t), intent(in   ) :: betay(lo(1)-ng_b:,lo(2)-ng_b:)
    real (kind = dp_t), intent(in   ) :: dh(:)

    integer          :: i, j, nsten
    real (kind=dp_t) :: hx2,hy2,hx22,hy22
    real (kind=dp_t) :: rholy,rhory,rhotx,rhobx  ! Transverse derivatives
    real (kind=dp_t) :: s1,s2,scale              ! Coefficients for slopes

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
    ss = ZERO

    ! These are the coefficients in the second order part
    hx22 = -ONE / (TWELVE * dh(1)**2 )
    hy22 = -ONE / (TWELVE * dh(2)**2 )

    !  These are the coefficients for the tangential slopes
    !  We don't divide by h because the product of slopes is multiplied by h^2/12
    s1 = 34._dp_t / 48._dp_t
    s2 = -5._dp_t / 48._dp_t

    !  In the output of make_stencil2d.f90, the coefficents of the correction are multiplied by
    !  all the denominators so that they are integers.  So we have to put the denominators back in
    scale = ONE/(TWELVE*48._dp_t)

    !  The coefficients hx2 and hy2 are defined by  (the minus sign is because it is minus beta*Lap)
    hx2 = -ONE/(TWELVE*dh(1)**2)*scale
    hy2 = -ONE/(TWELVE*dh(2)**2)*scale

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)

          rholy = s1*(betax(i  ,j+1)- betax(i  ,j-1)) + s2*(betax(i  ,j+2)- betax(i  ,j-2))
          rhory = s1*(betax(i+1,j+1)- betax(i+1,j-1)) + s2*(betax(i+1,j+2)- betax(i+1,j-2))
          rhotx = s1*(betay(i+1,j+1)- betay(i-1,j+1)) + s2*(betay(i+2,j+1)- betay(i-2,j+1))
          rhobx = s1*(betay(i+1,j  )- betay(i-1,j  )) + s2*(betay(i+2,j  )- betay(i-2,j  ))

          !   DOING CONTRIB AT           -2          -2
         nsten =   1
           ss(nsten,i,j) = (   &
                             -5.0_dp_t*rholy*hx2  &
                             -5.0_dp_t*rhobx*hy2  )
                                               
            !   DOING CONTRIB AT           -1          -2

         nsten =   2
           ss(nsten,i,j) = (   &
                         +     5.0_dp_t*rhory*hx2  &
                         +    75.0_dp_t*rholy*hx2  &
                         +    34.0_dp_t*rhobx*hy2  )
                                               
            !   DOING CONTRIB AT            0          -2
         nsten =   3
           ss(nsten,i,j) = (   &
                            -75.0_dp_t*rhory*hx2  &
                            -75.0_dp_t*rholy*hx2  )
                                               
            !   DOING CONTRIB AT            1          -2
         nsten =   4
           ss(nsten,i,j) = (   &
                         +    75.0_dp_t*rhory*hx2  &
                         +     5.0_dp_t*rholy*hx2  &
                            -34.0_dp_t*rhobx*hy2   )
                                               
            !   DOING CONTRIB AT            2          -2
         nsten =   5
           ss(nsten,i,j) = (   &
                             -5.0_dp_t*rhory*hx2  &
                         +    5.0_dp_t*rhobx*hy2 )
                                               
            !   DOING CONTRIB AT           -2          -1
         nsten =   6
           ss(nsten,i,j) = (   &
                         +    34.0_dp_t*rholy*hx2  &
                         +     5.0_dp_t*rhotx*hy2  &
                         +    75.0_dp_t*rhobx*hy2  )
                                               
            !   DOING CONTRIB AT           -1          -1
         nsten =   7
           ss(nsten,i,j) = (   &
                            -34.0_dp_t*rhory*hx2  &
                           -510.0_dp_t*rholy*hx2  &
                            -34.0_dp_t*rhotx*hy2  &
                           -510.0_dp_t*rhobx*hy2  )
                                               
            !   DOING CONTRIB AT            0          -1
         nsten =   8
           ss(nsten,i,j) = (   &
                         +   510.0_dp_t*rhory*hx2  &
                         +   510.0_dp_t*rholy*hx2  )
                                               
            !   DOING CONTRIB AT            1          -1
         nsten =   9
           ss(nsten,i,j) = (   &
                           -510.0_dp_t*rhory*hx2  &
                            -34.0_dp_t*rholy*hx2  &
                         +   34.0_dp_t*rhotx*hy2  &
                         +  510.0_dp_t*rhobx*hy2  )
                                               
            !   DOING CONTRIB AT            2          -1
         nsten =  10
           ss(nsten,i,j) = (   &
                         +    34.0_dp_t*rhory*hx2  &
                              -5.0_dp_t*rhotx*hy2  &
                             -75.0_dp_t*rhobx*hy2  )
                                               
            !   DOING CONTRIB AT           -2           0
         nsten =  11
           ss(nsten,i,j) = (   &
                            -75.0_dp_t*rhotx*hy2  &
                            -75.0_dp_t*rhobx*hy2  )
                                               
            !   DOING CONTRIB AT           -1           0
         nsten =  12
           ss(nsten,i,j) = (   &
                         +   510.0_dp_t*rhotx*hy2  &
                         +   510.0_dp_t*rhobx*hy2  )
                                               
            !   DOING CONTRIB AT            1           0
         nsten =  13
           ss(nsten,i,j) = (   &
                           -510.0_dp_t*rhotx*hy2  &
                           -510.0_dp_t*rhobx*hy2  )
                                               
            !   DOING CONTRIB AT            2           0
         nsten =  14
           ss(nsten,i,j) = (   &
                         +    75.0_dp_t*rhotx*hy2  &
                         +    75.0_dp_t*rhobx*hy2  )
                                               
            !   DOING CONTRIB AT           -2           1
         nsten =  15
           ss(nsten,i,j) = (   &
                            -34.0_dp_t*rholy*hx2  &
                         +    75.0_dp_t*rhotx*hy2  &
                         +     5.0_dp_t*rhobx*hy2  )
                                               
            !   DOING CONTRIB AT           -1           1
         nsten =  16
           ss(nsten,i,j) = (   &
                         +    34.0_dp_t*rhory*hx2  &
                         +   510.0_dp_t*rholy*hx2  &
                           -510.0_dp_t*rhotx*hy2  &
                            -34.0_dp_t*rhobx*hy2  )
                                               
            !   DOING CONTRIB AT            0           1
         nsten =  17
           ss(nsten,i,j) = (   &
                           -510.0_dp_t*rhory*hx2  &
                           -510.0_dp_t*rholy*hx2  )
                                               
            !   DOING CONTRIB AT            1           1
         nsten =  18
           ss(nsten,i,j) = (   &
                         +   510.0_dp_t*rhory*hx2  &
                         +    34.0_dp_t*rholy*hx2  &
                         +   510.0_dp_t*rhotx*hy2  &
                         +    34.0_dp_t*rhobx*hy2  )
                                               
            !   DOING CONTRIB AT            2           1
         nsten =  19
           ss(nsten,i,j) = (   &
                            -34.0_dp_t*rhory*hx2  &
                            -75.0_dp_t*rhotx*hy2  &
                             -5.0_dp_t*rhobx*hy2  )
                                               
            !   DOING CONTRIB AT           -2           2
         nsten =  20
           ss(nsten,i,j) = (   &
                         +     5.0_dp_t*rholy*hx2 &
                             -5.0_dp_t*rhotx*hy2  )
                                              
            !   DOING CONTRIB AT           -1           2
         nsten =  21
           ss(nsten,i,j) = (   &
                             -5.0_dp_t*rhory*hx2  &
                            -75.0_dp_t*rholy*hx2  &
                         +    34.0_dp_t*rhotx*hy2 )
                                              
            !   DOING CONTRIB AT            0           2
         nsten =  22
           ss(nsten,i,j) = (   &
                         +    75.0_dp_t*rhory*hx2  &
                         +    75.0_dp_t*rholy*hx2  )
                                              
            !   DOING CONTRIB AT            1           2
         nsten =  23
           ss(nsten,i,j) = (   &
                            -75.0_dp_t*rhory*hx2  &
                             -5.0_dp_t*rholy*hx2  &
                            -34.0_dp_t*rhotx*hy2  )

            !   DOING CONTRIB AT            2           2
         nsten =  24
           ss(nsten,i,j) = (   &
                         +     5.0_dp_t*rhory*hx2  &
                         +     5.0_dp_t*rhotx*hy2  )

          !  Now we add in the 2nd order stencil
          ss(11,i,j) = ss(11,i,j) + (                                  - betax(i,j))*hx22
          ss(12,i,j) = ss(12,i,j) + (           betax(i+1,j) + 15.0_dp_t*betax(i,j))*hx22
          ss( 0,i,j) = ss( 0,i,j) + (-15.0_dp_t*betax(i+1,j) - 15.0_dp_t*betax(i,j))*hx22
          ss(13,i,j) = ss(13,i,j) + ( 15.0_dp_t*betax(i+1,j) +           betax(i,j))*hx22
          ss(14,i,j) = ss(14,i,j) + (          -betax(i+1,j)                       )*hx22

          ss( 3,i,j) = ss( 3,i,j) + (                                  - betay(i,j))*hy22
          ss( 8,i,j) = ss( 8,i,j) + (           betay(i,j+1) + 15.0_dp_t*betay(i,j))*hy22
          ss( 0,i,j) = ss( 0,i,j) + (-15.0_dp_t*betay(i,j+1) - 15.0_dp_t*betay(i,j))*hy22
          ss(17,i,j) = ss(17,i,j) + ( 15.0_dp_t*betay(i,j+1) +           betay(i,j))*hy22
          ss(22,i,j) = ss(22,i,j) + (          -betay(i,j+1)                       )*hy22
       end do
    end do

    ! This adds the "alpha" term in (alpha - del dot beta grad) phi = RHS.
    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          ss(0,i,j) = ss(0,i,j) + alpha(i,j)
       end do
    end do

  end subroutine s_minion_full_fill_2d

  subroutine s_minion_cross_fill_3d(ss, alpha, ng_a, betax, betay, betaz, ng_b, dh, mask, lo, hi)

    integer           , intent(in   ) :: ng_a,ng_b
    integer           , intent(in   ) :: lo(:), hi(:)
    integer           , intent(inout) :: mask(lo(1) :,lo(2) :,lo(3) :)
    real (kind = dp_t), intent(  out) :: ss(0:,lo(1) :,lo(2) :,lo(3) :)
    real (kind = dp_t), intent(in   ) :: alpha(lo(1)-ng_a:,lo(2)-ng_a:,lo(3)-ng_a:)
    real (kind = dp_t), intent(in   ) :: betax(lo(1)-ng_b:,lo(2)-ng_b:,lo(3)-ng_b:)
    real (kind = dp_t), intent(in   ) :: betay(lo(1)-ng_b:,lo(2)-ng_b:,lo(3)-ng_b:)
    real (kind = dp_t), intent(in   ) :: betaz(lo(1)-ng_b:,lo(2)-ng_b:,lo(3)-ng_b:)
    real (kind = dp_t), intent(in   ) :: dh(:)

    integer            :: i, j, k
    real (kind = dp_t) :: fac

    fac = (ONE / (TWELVE * dh(1)**2))

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             mask(i,j,k) = ibclr(mask(i,j,k), BC_BIT(BC_GEOM,1,-1))
             mask(i,j,k) = ibclr(mask(i,j,k), BC_BIT(BC_GEOM,1,+1))
             mask(i,j,k) = ibclr(mask(i,j,k), BC_BIT(BC_GEOM,2,-1))
             mask(i,j,k) = ibclr(mask(i,j,k), BC_BIT(BC_GEOM,2,+1))
             mask(i,j,k) = ibclr(mask(i,j,k), BC_BIT(BC_GEOM,3,-1))
             mask(i,j,k) = ibclr(mask(i,j,k), BC_BIT(BC_GEOM,3,+1))
          end do
       end do
    end do
    !$OMP END PARALLEL DO
    !
    ! We only include the beta's here to get the viscous coefficients in here for now.
    ! The projection has beta == 1.
    !
    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             ss( 1,i,j,k) =      ONE * betax(i  ,j,k) * fac
             ss( 2,i,j,k) = -SIXTEEN * betax(i  ,j,k) * fac
             ss( 3,i,j,k) = -SIXTEEN * betax(i+1,j,k) * fac
             ss( 4,i,j,k) =      ONE * betax(i+1,j,k) * fac
             ss( 5,i,j,k) =      ONE * betay(i,j  ,k) * fac
             ss( 6,i,j,k) = -SIXTEEN * betay(i,j  ,k) * fac
             ss( 7,i,j,k) = -SIXTEEN * betay(i,j+1,k) * fac
             ss( 8,i,j,k) =      ONE * betay(i,j+1,k) * fac
             ss( 9,i,j,k) =      ONE * betaz(i,j,k  ) * fac
             ss(10,i,j,k) = -SIXTEEN * betaz(i,j,k  ) * fac
             ss(11,i,j,k) = -SIXTEEN * betaz(i,j,k+1) * fac
             ss(12,i,j,k) =      ONE * betaz(i,j,k+1) * fac
             ss(0,i,j,k)  = -( ss(1,i,j,k) + ss( 2,i,j,k) + ss( 3,i,j,k) + ss( 4,i,j,k) &
                              +ss(5,i,j,k) + ss( 6,i,j,k) + ss( 7,i,j,k) + ss( 8,i,j,k) &
                              +ss(9,i,j,k) + ss(10,i,j,k) + ss(11,i,j,k) + ss(12,i,j,k) ) + alpha(i,j,k)
          end do
       end do
    end do
    !$OMP END PARALLEL DO

  end subroutine s_minion_cross_fill_3d

  subroutine s_minion_full_fill_3d(ss, alpha, ng_a, betax, betay, betaz, ng_b, dh, mask, lo, hi)

    integer           , intent(in   ) :: ng_a,ng_b
    integer           , intent(in   ) :: lo(:), hi(:)
    integer           , intent(inout) :: mask(lo(1)      :,lo(2)     :,lo(3)     :)
    real (kind = dp_t), intent(  out) :: ss(0:,lo(1)     :,lo(2)     :,lo(3)     :)
    real (kind = dp_t), intent(in   ) :: alpha(lo(1)-ng_a:,lo(2)-ng_a:,lo(3)-ng_a:)
    real (kind = dp_t), intent(inout) :: betax(lo(1)-ng_b:,lo(2)-ng_b:,lo(3)-ng_b:)
    real (kind = dp_t), intent(inout) :: betay(lo(1)-ng_b:,lo(2)-ng_b:,lo(3)-ng_b:)
    real (kind = dp_t), intent(inout) :: betaz(lo(1)-ng_b:,lo(2)-ng_b:,lo(3)-ng_b:)
    real (kind = dp_t), intent(in   ) :: dh(:)

    integer          :: i, j, k, nsten
    real (kind=dp_t) :: hx2,hy2,hz2,hx22,hy22,hz22

   ! Transverse derivatives: left,right,top,bottom,fore,aft
    real (kind=dp_t) :: rhoax,rhobx,rhofx,rhotx
    real (kind=dp_t) :: rhoby,rholy,rhory,rhoty
    real (kind=dp_t) :: rhofz,rholz,rhorz,rhoaz
    real (kind=dp_t) :: s1,s2,scale
    real (kind=dp_t) :: sum

    !  These are the coefficients in the second order part
    hx22 = -ONE / (TWELVE * dh(1)**2 )
    hy22 = -ONE / (TWELVE * dh(2)**2 )
    hz22 = -ONE / (TWELVE * dh(3)**2 )

    !  These are the coefficients for the tangential slopes
    !  We don't divide by h because the product of slopes is multiplied by h^2/12
    s1 = 34._dp_t / 48._dp_t
    s2 = -5._dp_t / 48._dp_t

    !  In the output of make_stencil3d.f90, the coefficents of the correction are multiplied by
    !  all the denominators so that they are integers.  So we have to put the denominators back in
    scale = ONE/(TWELVE*48._dp_t)

    !  The coefficients hx2 and hy2 are defined by  (the minus sign is because it is minus beta*Lap)
    hx2 = -ONE/(TWELVE*dh(1)**2)*scale
    hy2 = -ONE/(TWELVE*dh(2)**2)*scale
    hz2 = -ONE/(TWELVE*dh(3)**2)*scale

    !  Initialize to zero.
    ss = 0._dp_t

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             mask(i,j,k) = ibclr(mask(i,j,k), BC_BIT(BC_GEOM,1,-1))
             mask(i,j,k) = ibclr(mask(i,j,k), BC_BIT(BC_GEOM,1,+1))
             mask(i,j,k) = ibclr(mask(i,j,k), BC_BIT(BC_GEOM,2,-1))
             mask(i,j,k) = ibclr(mask(i,j,k), BC_BIT(BC_GEOM,2,+1))
             mask(i,j,k) = ibclr(mask(i,j,k), BC_BIT(BC_GEOM,3,-1))
             mask(i,j,k) = ibclr(mask(i,j,k), BC_BIT(BC_GEOM,3,+1))
          end do
       end do
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i,j,k,nsten,rhoax,rhobx,rhofx,rhotx,rhoby,rholy,rhory,rhoty,rhofz,rholz,rhorz,rhoaz)
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
          ss(nsten,i,j,k) = ( &
                        -5.0_dp_t*rhoby*hz2 &
                        -5.0_dp_t*rhoaz*hy2 )
                                              
 ! DOING CONTRIB AT            0          -1          -2 nsten =            2
        nsten =   2
          ss(nsten,i,j,k) = ( &
                   +    34.0_dp_t*rhoby*hz2 &
                   +     5.0_dp_t*rhofz*hy2 &
                   +    75.0_dp_t*rhoaz*hy2 )
                                              
 ! DOING CONTRIB AT           -2           0          -2 nsten =            3
        nsten =   3
          ss(nsten,i,j,k) = ( &
                        -5.0_dp_t*rholz*hx2 &
                        -5.0_dp_t*rhobx*hz2 )
                                              
 ! DOING CONTRIB AT           -1           0          -2 nsten =            4
        nsten =   4
          ss(nsten,i,j,k) = ( &
                   +    75.0_dp_t*rholz*hx2 &
                   +    34.0_dp_t*rhobx*hz2 )
                                              
 ! DOING CONTRIB AT            0           0          -2 nsten =            5
        nsten =   5
          ss(nsten,i,j,k) = ( &
                       -75.0_dp_t*rholz*hx2 &
                       -75.0_dp_t*rhofz*hy2 &
                       -75.0_dp_t*rhoaz*hy2 )
                                              
 ! DOING CONTRIB AT            1           0          -2 nsten =            6
        nsten =   6
          ss(nsten,i,j,k) = ( &
                   +     5.0_dp_t*rholz*hx2 &
                       -34.0_dp_t*rhobx*hz2 )
                                              
 ! DOING CONTRIB AT            2           0          -2 nsten =            7
        nsten =   7
          ss(nsten,i,j,k) = ( &
                   +     5.0_dp_t*rhobx*hz2 )
                                              
 ! DOING CONTRIB AT            0           1          -2 nsten =            8
        nsten =   8
          ss(nsten,i,j,k) = ( &
                       -34.0_dp_t*rhoby*hz2 &
                   +    75.0_dp_t*rhofz*hy2 &
                   +     5.0_dp_t*rhoaz*hy2 )
                                              
 ! DOING CONTRIB AT            0           2          -2 nsten =            9
        nsten =   9
          ss(nsten,i,j,k) = ( &
                   +     5.0_dp_t*rhoby*hz2 &
                        -5.0_dp_t*rhofz*hy2 )
                                              
 ! DOING CONTRIB AT            0          -2          -1 nsten =           10
        nsten =  10
          ss(nsten,i,j,k) = ( &
                   +     5.0_dp_t*rhoty*hz2 &
                   +    75.0_dp_t*rhoby*hz2 &
                   +    34.0_dp_t*rhoaz*hy2 )
                                              
 ! DOING CONTRIB AT            0          -1          -1 nsten =           11
        nsten =  11
          ss(nsten,i,j,k) = ( &
                       -34.0_dp_t*rhoty*hz2 &
                      -510.0_dp_t*rhoby*hz2 &
                       -34.0_dp_t*rhofz*hy2 &
                      -510.0_dp_t*rhoaz*hy2 )
                                              
 ! DOING CONTRIB AT           -2           0          -1 nsten =           12
        nsten =  12
          ss(nsten,i,j,k) = ( &
                   +    34.0_dp_t*rholz*hx2 &
                   +     5.0_dp_t*rhotx*hz2 &
                   +    75.0_dp_t*rhobx*hz2 )
                                              
 ! DOING CONTRIB AT           -1           0          -1 nsten =           13
        nsten =  13
          ss(nsten,i,j,k) = ( &
                      -510.0_dp_t*rholz*hx2 &
                       -34.0_dp_t*rhotx*hz2 &
                      -510.0_dp_t*rhobx*hz2 )
                                              
 ! DOING CONTRIB AT            0           0          -1 nsten =           14
        nsten =  14
          ss(nsten,i,j,k) = ( &
                   +   510.0_dp_t*rholz*hx2 &
                   +   510.0_dp_t*rhofz*hy2 &
                   +   510.0_dp_t*rhoaz*hy2 )
                                              
 ! DOING CONTRIB AT            1           0          -1 nsten =           15
        nsten =  15
          ss(nsten,i,j,k) = ( &
                       -34.0_dp_t*rholz*hx2 &
                   +    34.0_dp_t*rhotx*hz2 &
                   +   510.0_dp_t*rhobx*hz2 )
                                              
 ! DOING CONTRIB AT            2           0          -1 nsten =           16
        nsten =  16
          ss(nsten,i,j,k) = ( &
                        -5.0_dp_t*rhotx*hz2 &
                       -75.0_dp_t*rhobx*hz2 )
                                              
 ! DOING CONTRIB AT            0           1          -1 nsten =           17
        nsten =  17
          ss(nsten,i,j,k) = ( &
                   +    34.0_dp_t*rhoty*hz2 &
                   +   510.0_dp_t*rhoby*hz2 &
                      -510.0_dp_t*rhofz*hy2 &
                       -34.0_dp_t*rhoaz*hy2 )
                                              
 ! DOING CONTRIB AT            0           2          -1 nsten =           18
        nsten =  18
          ss(nsten,i,j,k) = ( &
                        -5.0_dp_t*rhoty*hz2 &
                       -75.0_dp_t*rhoby*hz2 &
                   +    34.0_dp_t*rhofz*hy2 )
                                              
 ! DOING CONTRIB AT           -2          -2           0 nsten =           19
        nsten =  19
          ss(nsten,i,j,k) = ( &
                        -5.0_dp_t*rholy*hx2 &
                        -5.0_dp_t*rhoax*hy2 )
                                              
 ! DOING CONTRIB AT           -1          -2           0 nsten =           20
        nsten =  20
          ss(nsten,i,j,k) = ( &
                   +     5.0_dp_t*rhory*hx2 &
                   +     5.0_dp_t*rhorz*hx2 &
                   +    75.0_dp_t*rholy*hx2 &
                   +    34.0_dp_t*rhoax*hy2 )
                                              
 ! DOING CONTRIB AT            0          -2           0 nsten =           21
        nsten =  21
          ss(nsten,i,j,k) = ( &
                       -75.0_dp_t*rhory*hx2 &
                       -75.0_dp_t*rhorz*hx2 &
                       -75.0_dp_t*rholy*hx2 &
                       -75.0_dp_t*rhoty*hz2 &
                       -75.0_dp_t*rhoby*hz2 )
                                              
 ! DOING CONTRIB AT            1          -2           0 nsten =           22
        nsten =  22
          ss(nsten,i,j,k) = ( &
                   +    75.0_dp_t*rhory*hx2 &
                   +    75.0_dp_t*rhorz*hx2 &
                   +     5.0_dp_t*rholy*hx2 &
                       -34.0_dp_t*rhoax*hy2 )
                                              
 ! DOING CONTRIB AT            2          -2           0 nsten =           23
        nsten =  23
          ss(nsten,i,j,k) = ( &
                        -5.0_dp_t*rhory*hx2 &
                        -5.0_dp_t*rhorz*hx2 &
                   +     5.0_dp_t*rhoax*hy2 )
                                              
 ! DOING CONTRIB AT           -2          -1           0 nsten =           24
        nsten =  24
          ss(nsten,i,j,k) = ( &
                   +    34.0_dp_t*rholy*hx2 &
                   +     5.0_dp_t*rhofx*hy2 &
                   +    75.0_dp_t*rhoax*hy2 )
                                              
 ! DOING CONTRIB AT           -1          -1           0 nsten =           25
        nsten =  25
          ss(nsten,i,j,k) = ( &
                       -34.0_dp_t*rhory*hx2 &
                       -34.0_dp_t*rhorz*hx2 &
                      -510.0_dp_t*rholy*hx2 &
                       -34.0_dp_t*rhofx*hy2 &
                      -510.0_dp_t*rhoax*hy2 )
                                              
 ! DOING CONTRIB AT            0          -1           0 nsten =           26
        nsten =  26
          ss(nsten,i,j,k) = ( &
                   +   510.0_dp_t*rhory*hx2 &
                   +   510.0_dp_t*rhorz*hx2 &
                   +   510.0_dp_t*rholy*hx2 &
                   +   510.0_dp_t*rhoty*hz2 &
                   +   510.0_dp_t*rhoby*hz2 )
                                              
 ! DOING CONTRIB AT            1          -1           0 nsten =           27
        nsten =  27
          ss(nsten,i,j,k) = ( &
                      -510.0_dp_t*rhory*hx2 &
                      -510.0_dp_t*rhorz*hx2 &
                       -34.0_dp_t*rholy*hx2 &
                   +    34.0_dp_t*rhofx*hy2 &
                   +   510.0_dp_t*rhoax*hy2 )
                                              
 ! DOING CONTRIB AT            2          -1           0 nsten =           28
        nsten =  28
          ss(nsten,i,j,k) = ( &
                   +    34.0_dp_t*rhory*hx2 &
                   +    34.0_dp_t*rhorz*hx2 &
                        -5.0_dp_t*rhofx*hy2 &
                       -75.0_dp_t*rhoax*hy2 )
                                              
 ! DOING CONTRIB AT           -2           0           0 nsten =           29
        nsten =  29
          ss(nsten,i,j,k) = ( &
                       -75.0_dp_t*rhotx*hz2 &
                       -75.0_dp_t*rhobx*hz2 &
                       -75.0_dp_t*rhofx*hy2 &
                       -75.0_dp_t*rhoax*hy2 )
                                              
 ! DOING CONTRIB AT           -1           0           0 nsten =           30
        nsten =  30
          ss(nsten,i,j,k) = ( &
                   +   510.0_dp_t*rhotx*hz2 &
                   +   510.0_dp_t*rhobx*hz2 &
                   +   510.0_dp_t*rhofx*hy2 &
                   +   510.0_dp_t*rhoax*hy2 )
                                              
 ! DOING CONTRIB AT            1           0           0 nsten =           31
        nsten =  31
          ss(nsten,i,j,k) = ( &
                      -510.0_dp_t*rhotx*hz2 &
                      -510.0_dp_t*rhobx*hz2 &
                      -510.0_dp_t*rhofx*hy2 &
                      -510.0_dp_t*rhoax*hy2 )
                                              
 ! DOING CONTRIB AT            2           0           0 nsten =           32
        nsten =  32
          ss(nsten,i,j,k) = ( &
                   +    75.0_dp_t*rhotx*hz2 &
                   +    75.0_dp_t*rhobx*hz2 &
                   +    75.0_dp_t*rhofx*hy2 &
                   +    75.0_dp_t*rhoax*hy2 )
                                              
 ! DOING CONTRIB AT           -2           1           0 nsten =           33
        nsten =  33
          ss(nsten,i,j,k) = ( &
                       -34.0_dp_t*rholy*hx2 &
                   +    75.0_dp_t*rhofx*hy2 &
                   +     5.0_dp_t*rhoax*hy2 )
                                              
 ! DOING CONTRIB AT           -1           1           0 nsten =           34
        nsten =  34
          ss(nsten,i,j,k) = ( &
                   +    34.0_dp_t*rhory*hx2 &
                   +    34.0_dp_t*rhorz*hx2 &
                   +   510.0_dp_t*rholy*hx2 &
                      -510.0_dp_t*rhofx*hy2 &
                       -34.0_dp_t*rhoax*hy2 )
                                              
 ! DOING CONTRIB AT            0           1           0 nsten =           35
        nsten =  35
          ss(nsten,i,j,k) = ( &
                      -510.0_dp_t*rhory*hx2 &
                      -510.0_dp_t*rhorz*hx2 &
                      -510.0_dp_t*rholy*hx2 &
                      -510.0_dp_t*rhoty*hz2 &
                      -510.0_dp_t*rhoby*hz2 )
                                              
 ! DOING CONTRIB AT            1           1           0 nsten =           36
        nsten =  36
          ss(nsten,i,j,k) = ( &
                   +   510.0_dp_t*rhory*hx2 &
                   +   510.0_dp_t*rhorz*hx2 &
                   +    34.0_dp_t*rholy*hx2 &
                   +   510.0_dp_t*rhofx*hy2 &
                   +    34.0_dp_t*rhoax*hy2 )
                                              
 ! DOING CONTRIB AT            2           1           0 nsten =           37
        nsten =  37
          ss(nsten,i,j,k) = ( &
                       -34.0_dp_t*rhory*hx2 &
                       -34.0_dp_t*rhorz*hx2 &
                       -75.0_dp_t*rhofx*hy2 &
                        -5.0_dp_t*rhoax*hy2 )
                                              
 ! DOING CONTRIB AT           -2           2           0 nsten =           38
        nsten =  38
          ss(nsten,i,j,k) = ( &
                   +     5.0_dp_t*rholy*hx2 &
                        -5.0_dp_t*rhofx*hy2 )
                                              
 ! DOING CONTRIB AT           -1           2           0 nsten =           39
        nsten =  39
          ss(nsten,i,j,k) = ( &
                        -5.0_dp_t*rhory*hx2 &
                        -5.0_dp_t*rhorz*hx2 &
                       -75.0_dp_t*rholy*hx2 &
                   +    34.0_dp_t*rhofx*hy2 )
                                              
 ! DOING CONTRIB AT            0           2           0 nsten =           40
        nsten =  40
          ss(nsten,i,j,k) = ( &
                   +    75.0_dp_t*rhory*hx2 &
                   +    75.0_dp_t*rhorz*hx2 &
                   +    75.0_dp_t*rholy*hx2 &
                   +    75.0_dp_t*rhoty*hz2 &
                   +    75.0_dp_t*rhoby*hz2 )
                                              
 ! DOING CONTRIB AT            1           2           0 nsten =           41
        nsten =  41
          ss(nsten,i,j,k) = ( &
                       -75.0_dp_t*rhory*hx2 &
                       -75.0_dp_t*rhorz*hx2 &
                        -5.0_dp_t*rholy*hx2 &
                       -34.0_dp_t*rhofx*hy2 )
                                              
 ! DOING CONTRIB AT            2           2           0 nsten =           42
        nsten =  42
          ss(nsten,i,j,k) = ( &
                   +     5.0_dp_t*rhory*hx2 &
                   +     5.0_dp_t*rhorz*hx2 &
                   +     5.0_dp_t*rhofx*hy2 )
                                              
 ! DOING CONTRIB AT            0          -2           1 nsten =           43
        nsten =  43
          ss(nsten,i,j,k) = ( &
                   +    75.0_dp_t*rhoty*hz2 &
                   +     5.0_dp_t*rhoby*hz2 &
                       -34.0_dp_t*rhoaz*hy2 )
                                              
 ! DOING CONTRIB AT            0          -1           1 nsten =           44
        nsten =  44
          ss(nsten,i,j,k) = ( &
                      -510.0_dp_t*rhoty*hz2 &
                       -34.0_dp_t*rhoby*hz2 &
                   +    34.0_dp_t*rhofz*hy2 &
                   +   510.0_dp_t*rhoaz*hy2 )
                                              
 ! DOING CONTRIB AT           -2           0           1 nsten =           45
        nsten =  45
          ss(nsten,i,j,k) = ( &
                       -34.0_dp_t*rholz*hx2 &
                   +    75.0_dp_t*rhotx*hz2 &
                   +     5.0_dp_t*rhobx*hz2 )
                                              
 ! DOING CONTRIB AT           -1           0           1 nsten =           46
        nsten =  46
          ss(nsten,i,j,k) = ( &
                   +   510.0_dp_t*rholz*hx2 &
                      -510.0_dp_t*rhotx*hz2 &
                       -34.0_dp_t*rhobx*hz2 )
                                              
 ! DOING CONTRIB AT            0           0           1 nsten =           47
        nsten =  47
          ss(nsten,i,j,k) = ( &
                      -510.0_dp_t*rholz*hx2 &
                      -510.0_dp_t*rhofz*hy2 &
                      -510.0_dp_t*rhoaz*hy2 )
                                              
 ! DOING CONTRIB AT            1           0           1 nsten =           48
        nsten =  48
          ss(nsten,i,j,k) = ( &
                   +    34.0_dp_t*rholz*hx2 &
                   +   510.0_dp_t*rhotx*hz2 &
                   +    34.0_dp_t*rhobx*hz2 )
                                              
 ! DOING CONTRIB AT            2           0           1 nsten =           49
        nsten =  49
          ss(nsten,i,j,k) = ( &
                       -75.0_dp_t*rhotx*hz2 &
                        -5.0_dp_t*rhobx*hz2 )
                                              
 ! DOING CONTRIB AT            0           1           1 nsten =           50
        nsten =  50
          ss(nsten,i,j,k) = ( &
                   +   510.0_dp_t*rhoty*hz2 &
                   +    34.0_dp_t*rhoby*hz2 &
                   +   510.0_dp_t*rhofz*hy2 &
                   +    34.0_dp_t*rhoaz*hy2 )
                                              
 ! DOING CONTRIB AT            0           2           1 nsten =           51
        nsten =  51
          ss(nsten,i,j,k) = ( &
                       -75.0_dp_t*rhoty*hz2 &
                        -5.0_dp_t*rhoby*hz2 &
                       -34.0_dp_t*rhofz*hy2 )
                                              
 ! DOING CONTRIB AT            0          -2           2 nsten =           52
        nsten =  52
          ss(nsten,i,j,k) = ( &
                        -5.0_dp_t*rhoty*hz2 &
                   +     5.0_dp_t*rhoaz*hy2 )
                                              
 ! DOING CONTRIB AT            0          -1           2 nsten =           53
        nsten =  53
          ss(nsten,i,j,k) = ( &
                   +    34.0_dp_t*rhoty*hz2 &
                        -5.0_dp_t*rhofz*hy2 &
                       -75.0_dp_t*rhoaz*hy2 )
                                              
 ! DOING CONTRIB AT           -2           0           2 nsten =           54
        nsten =  54
          ss(nsten,i,j,k) = ( &
                   +     5.0_dp_t*rholz*hx2 &
                        -5.0_dp_t*rhotx*hz2 )
                                              
 ! DOING CONTRIB AT           -1           0           2 nsten =           55
        nsten =  55
          ss(nsten,i,j,k) = ( &
                       -75.0_dp_t*rholz*hx2 &
                   +    34.0_dp_t*rhotx*hz2 )
                                              
 ! DOING CONTRIB AT            0           0           2 nsten =           56
        nsten =  56
          ss(nsten,i,j,k) = ( &
                   +    75.0_dp_t*rholz*hx2 &
                   +    75.0_dp_t*rhofz*hy2 &
                   +    75.0_dp_t*rhoaz*hy2 )
                                              
 ! DOING CONTRIB AT            1           0           2 nsten =           57
        nsten =  57
          ss(nsten,i,j,k) = ( &
                        -5.0_dp_t*rholz*hx2 &
                       -34.0_dp_t*rhotx*hz2 )
                                              
 ! DOING CONTRIB AT            2           0           2 nsten =           58
        nsten =  58
          ss(nsten,i,j,k) = ( &
                   +     5.0_dp_t*rhotx*hz2 )
                                              
 ! DOING CONTRIB AT            0           1           2 nsten =           59
        nsten =  59
          ss(nsten,i,j,k) = ( &
                       -34.0_dp_t*rhoty*hz2 &
                       -75.0_dp_t*rhofz*hy2 &
                        -5.0_dp_t*rhoaz*hy2 )
                                              
 ! DOING CONTRIB AT            0           2           2 nsten =           60
        nsten =  60
          ss(nsten,i,j,k) = ( &
                   +     5.0_dp_t*rhoty*hz2 &
                   +     5.0_dp_t*rhofz*hy2 )
                                              

          !  Now we add in the 2nd order stencil
          ss(29,i,j,k) = ss(29,i,j,k) + (                          -           betax(i,j,k))*hx22
          ss(30,i,j,k) = ss(30,i,j,k) + (           betax(i+1,j,k) + 15.0_dp_t*betax(i,j,k))*hx22
          !ss(0,i,j,k) = ss( 0,i,j,k) + (-15.0_dp_t*betax(i+1,j,k) - 15.0_dp_t*betax(i,j,k))*hx22
          ss(31,i,j,k) = ss(31,i,j,k) + ( 15.0_dp_t*betax(i+1,j,k) +           betax(i,j,k))*hx22
          ss(32,i,j,k) = ss(32,i,j,k) + (          -betax(i+1,j,k)                         )*hx22
 
          ss(21,i,j,k) = ss(21,i,j,k) + (                          -           betay(i,j,k))*hy22
          ss(26,i,j,k) = ss(26,i,j,k) + (           betay(i,j+1,k) + 15.0_dp_t*betay(i,j,k))*hy22
          !ss(0,i,j,k) = ss( 0,i,j,k) + (-15.0_dp_t*betay(i,j+1,k) - 15.0_dp_t*betay(i,j,k))*hy22
          ss(35,i,j,k) = ss(35,i,j,k) + ( 15.0_dp_t*betay(i,j+1,k) +           betay(i,j,k))*hy22
          ss(40,i,j,k) = ss(40,i,j,k) + (          -betay(i,j+1,k)                         )*hy22
 
          ss( 5,i,j,k) = ss( 5,i,j,k) + (                          -           betaz(i,j,k))*hz22
          ss(14,i,j,k) = ss(14,i,j,k) + (           betaz(i,j,k+1) + 15.0_dp_t*betaz(i,j,k))*hz22
          !ss(0,i,j,k) = ss( 0,i,j,k) + (-15.0_dp_t*betaz(i,j,k+1) - 15.0_dp_t*betaz(i,j,k))*hz22
          ss(47,i,j,k) = ss(47,i,j,k) + ( 15.0_dp_t*betaz(i,j,k+1) +           betaz(i,j,k))*hz22
          ss(56,i,j,k) = ss(56,i,j,k) + (          -betaz(i,j,k+1)                         )*hz22

          end do
       end do
    end do
    !$OMP END PARALLEL DO

    !
    ! Define the center stencil and add the "alpha" term in 
    !     (alpha - del dot beta grad) phi = RHS.
    !$OMP PARALLEL DO PRIVATE(i,j,k,sum,nsten)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
          sum = ZERO
          do nsten = 1,60
             sum = sum + ss(nsten,i,j,k)
          end do
          ss(0,i,j,k) = alpha(i,j,k) - sum
          end do
       end do
    end do
    !$OMP END PARALLEL DO

  end subroutine s_minion_full_fill_3d

  subroutine simple_2d_const(ss, alpha_const, beta_const, dh, mask, lo, hi, xa, xb, order)


    integer           , intent(in   ) :: lo(:), hi(:), order
    integer           , intent(inout) :: mask(lo(1)  :,lo(2)  :)
    real (kind = dp_t), intent(  out) ::   ss(0:, lo(1)  :,lo(2)  :)
    real (kind = dp_t), intent(in   ) :: alpha_const, beta_const
    real (kind = dp_t), intent(in   ) :: dh(:), xa(:), xb(:)

    real (kind = dp_t) :: f1(2)
    integer            :: i, j, bclo, bchi, nx, ny
    integer            :: lnx, lny, lorder
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
          ss(0,i,j)   = ZERO
          ss(1,i,j)   = -beta_const*f1(1)
          ss(2,i,j)   = -beta_const*f1(1)
          ss(3,i,j)   = -beta_const*f1(2)
          ss(4,i,j)   = -beta_const*f1(2)
          ss(XBC,i,j) = ZERO
          ss(YBC,i,j) = ZERO
       end do
    end do

    lnx = nx; lny = ny; lorder = order

    ! x derivatives

    !$OMP PARALLEL DO PRIVATE(i,j)
    do j = lo(2),hi(2)
       do i = lo(1)+1,hi(1)-1
          ss(0,i,j) = ss(0,i,j) + TWO*beta_const*f1(1)
       end do
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i,j,bclo,bchi) FIRSTPRIVATE(lorder,lnx)
    do j = lo(2),hi(2)
       bclo = stencil_bc_type(mask(lo(1),j),1,-1)
       bchi = stencil_bc_type(mask(hi(1),j),1,+1)

       i = lo(1)
       if (bclo .eq. BC_INT) then
          ss(0,i,j) = ss(0,i,j) + TWO*beta_const*f1(1)
       else
          call stencil_bndry_aaa(lorder, lnx, 1, -1, mask(i,j), &
               ss(0,i,j), ss(1,i,j), ss(2,i,j), ss(XBC,i,j), &
               beta_const, beta_const, xa(1), xb(1), dh(1), bclo, bchi)
       end if

       if ( hi(1) > lo(1) ) then
          i = hi(1)
          if (bchi .eq. BC_INT) then
             ss(0,i,j) = ss(0,i,j) + TWO*beta_const*f1(1)
          else
             call stencil_bndry_aaa(lorder, lnx, 1, 1, mask(i,j), &
                  ss(0,i,j), ss(1,i,j), ss(2,i,j), ss(XBC,i,j), &
                  beta_const, beta_const, xa(1), xb(1), dh(1), bclo, bchi)
          end if
       end if
    end do
    !$OMP END PARALLEL DO

    ! y derivatives

    !$OMP PARALLEL DO PRIVATE(i,j)
    do j = lo(2)+1,hi(2)-1
       do i = lo(1),hi(1)
          ss(0,i,j) = ss(0,i,j) + TWO*beta_const*f1(2)
       end do
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i,j,bclo,bchi) FIRSTPRIVATE(lorder,lny)
    do i = lo(1),hi(1)
       bclo = stencil_bc_type(mask(i,lo(2)),2,-1)
       bchi = stencil_bc_type(mask(i,hi(2)),2,+1)

       j = lo(2)
       if (bclo .eq. BC_INT) then
          ss(0,i,j) = ss(0,i,j) + TWO*beta_const*f1(2)
       else
          call stencil_bndry_aaa(lorder, lny, 2, -1, mask(i,j), &
               ss(0,i,j), ss(3,i,j), ss(4,i,j),ss(YBC,i,j), &
               beta_const, beta_const, xa(2), xb(2), dh(2), bclo, bchi)
       end if
       if ( hi(2) > lo(2) ) then
          j = hi(2)
          if (bchi .eq. BC_INT) then
             ss(0,i,j) = ss(0,i,j) + TWO*beta_const*f1(2)
          else
             call stencil_bndry_aaa(lorder, lny, 2, 1, mask(i,j), &
               ss(0,i,j), ss(3,i,j), ss(4,i,j), ss(YBC,i,j), &
                  beta_const, beta_const, xa(2), xb(2), dh(2), bclo, bchi)
          end if
       end if
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i,j)
    do j = lo(2),hi(2)
       do i = lo(1),hi(1)
          ss(0,i,j) = ss(0,i,j) + alpha_const
       end do
    end do
    !$OMP END PARALLEL DO

  end subroutine simple_2d_const

  subroutine simple_3d_const(ss, alpha_const, beta_const, dh, mask, lo, hi, xa, xb, order)


    integer           , intent(in   ) :: lo(:), hi(:), order
    integer           , intent(inout) :: mask(lo(1)  :,lo(2)  :,lo(3)  :)
    real (kind = dp_t), intent(  out) ::   ss(0:, lo(1)  :,lo(2)  :,lo(3)  :)
    real (kind = dp_t), intent(in   ) :: alpha_const, beta_const
    real (kind = dp_t), intent(in   ) :: dh(:), xa(:), xb(:)

    real (kind = dp_t) :: f1(3)
    integer            :: i, j, k, bclo, bchi, nx, ny, nz
    integer            :: lnx, lny, lnz, lorder
    integer, parameter :: XBC = 7, YBC = 8, ZBC = 9

    nx = hi(1)-lo(1)+1
    ny = hi(2)-lo(2)+1
    nz = hi(3)-lo(3)+1
    f1 = ONE/dh**2

    lnx = nx; lny = ny; lnz = nz; lorder = order

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             ss(0,i,j,k)   = ZERO
             ss(1,i,j,k)   = -beta_const*f1(1)
             ss(2,i,j,k)   = -beta_const*f1(1)
             ss(3,i,j,k)   = -beta_const*f1(2)
             ss(4,i,j,k)   = -beta_const*f1(2)
             ss(5,i,j,k)   = -beta_const*f1(3)
             ss(6,i,j,k)   = -beta_const*f1(3)
             ss(XBC,i,j,k) = ZERO
             ss(YBC,i,j,k) = ZERO
             ss(ZBC,i,j,k) = ZERO
          end do
       end do
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             mask(i,j,k) = ibclr(mask(i,j,k), BC_BIT(BC_GEOM,1,-1))
             mask(i,j,k) = ibclr(mask(i,j,k), BC_BIT(BC_GEOM,1,+1))
             mask(i,j,k) = ibclr(mask(i,j,k), BC_BIT(BC_GEOM,2,-1))
             mask(i,j,k) = ibclr(mask(i,j,k), BC_BIT(BC_GEOM,2,+1))
             mask(i,j,k) = ibclr(mask(i,j,k), BC_BIT(BC_GEOM,3,-1))
             mask(i,j,k) = ibclr(mask(i,j,k), BC_BIT(BC_GEOM,3,+1))
          end do
       end do
    end do
    !$OMP END PARALLEL DO

    ! x derivatives

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1)+1,hi(1)-1
             ss(0,i,j,k) = ss(0,i,j,k) + TWO*beta_const*f1(1)
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
             ss(0,i,j,k) = ss(0,i,j,k) + TWO*beta_const*f1(1)
          else
             call stencil_bndry_aaa(lorder, lnx, 1, -1, mask(i,j,k), &
                  ss(0,i,j,k), ss(1,i,j,k), ss(2,i,j,k), ss(XBC,i,j,k), &
                  beta_const, beta_const, xa(1), xb(1), dh(1), bclo, bchi)
          end if

          if ( hi(1) > lo(1) ) then
             i = hi(1)
             if (bchi .eq. BC_INT) then
                ss(0,i,j,k) = ss(0,i,j,k) + TWO*beta_const*f1(1)
             else
                call stencil_bndry_aaa(lorder, lnx, 1, 1, mask(i,j,k), &
                     ss(0,i,j,k), ss(1,i,j,k), ss(2,i,j,k), ss(XBC,i,j,k), &
                     beta_const, beta_const, xa(1), xb(1), dh(1), bclo, bchi)
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
             ss(0,i,j,k) = ss(0,i,j,k) + TWO*beta_const*f1(2)
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
             ss(0,i,j,k) = ss(0,i,j,k) + TWO*beta_const*f1(2)
          else
             call stencil_bndry_aaa(lorder, lny, 2, -1, mask(i,j,k), &
                  ss(0,i,j,k), ss(3,i,j,k), ss(4,i,j,k),ss(YBC,i,j,k), &
                  beta_const, beta_const, xa(2), xb(2), dh(2), bclo, bchi)
          end if
          if ( hi(2) > lo(2) ) then
             j = hi(2)
             if (bchi .eq. BC_INT) then
                ss(0,i,j,k) = ss(0,i,j,k) + TWO*beta_const*f1(2)
             else
                call stencil_bndry_aaa(lorder, lny, 2, 1, mask(i,j,k), &
                     ss(0,i,j,k), ss(3,i,j,k), ss(4,i,j,k), ss(YBC,i,j,k), &
                     beta_const, beta_const, xa(2), xb(2), dh(2), bclo, bchi)
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
             ss(0,i,j,k) = ss(0,i,j,k) + TWO*beta_const*f1(3)
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
             ss(0,i,j,k) = ss(0,i,j,k) + TWO*beta_const*f1(3)
          else
             call stencil_bndry_aaa(lorder, lnz, 3, -1, mask(i,j,k), &
                  ss(0,i,j,k), ss(5,i,j,k), ss(6,i,j,k),ss(ZBC,i,j,k), &
                  beta_const, beta_const, xa(3), xb(3), dh(3), bclo, bchi)
          end if
          if ( hi(3) > lo(3) ) then
             k = hi(3)
             if (bchi .eq. BC_INT) then
                ss(0,i,j,k) = ss(0,i,j,k) + TWO*beta_const*f1(3)
             else
                call stencil_bndry_aaa(lorder, lnz, 3, 1, mask(i,j,k), &
                     ss(0,i,j,k), ss(5,i,j,k), ss(6,i,j,k), ss(ZBC,i,j,k), &
                     beta_const, beta_const, xa(3), xb(3), dh(3), bclo, bchi)
             end if
          end if
       end do
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             ss(0,i,j,k) = ss(0,i,j,k) + alpha_const
          end do
       end do
    end do
    !$OMP END PARALLEL DO

  end subroutine simple_3d_const

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

  ! subroutine t_polyInterpCoeffTest(norder)
  !   integer, intent(in) :: NORDER
  !   integer j
  !   real(kind=dp_t) c(0:NORDER-1), ci(0:NORDER-1)
  !   real(kind=dp_t) y(0:NORDER-1)
  !   real(kind=dp_t) x(0:NORDER-1)
  !   real(kind=dp_t) xInt

  !   call random_number(ci)

  !   j = 0
    
  !   x = (/ ZERO, (j+HALF,j=0,NORDER-2) /)
  !   do j = 0, NORDER-2
  !      y(j) = horner(x(j), ci)
  !   end do

  !   xInt = -HALF

  !   call poly_interp_coeff(c, xInt, x)

  !   print *, 'x = ', x
  !   print *, 'y = ', y
  !   print *, 'c = ', c
  !   print *, 'Interpolated y = ', sum(c*y)

  ! contains

  !   function Horner(xx, cc) result(r)
  !     real(kind=dp_t) :: r
  !     real(kind=dp_t), intent(in) :: xx
  !     real(kind=dp_t), intent(in) :: cc(:)
  !     integer :: i

  !     r = cc(1)
  !     do i = 2, size(cc)
  !        r = xx*r + cc(i)
  !     end do

  !   end function Horner

  ! end subroutine t_polyInterpCoeffTest

end module cc_stencil_module
