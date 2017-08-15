module multifab_c_module

  use iso_c_binding
  use layout_module
  use multifab_module

  implicit none

  type(layout), save :: la
  integer, save :: ndim

  integer, parameter :: len_name = 16

  type multifab_c
     character(len=len_name) :: name
     type(multifab) :: mf
  end type multifab_c

  integer, parameter :: max_mfs = 8
  type(multifab_c), save :: mfs(max_mfs)
  integer, save :: nmf = 0

  private

  public :: init_multifab_c, build_layout_from_c, destroy_multifab_c, &
       share_multifab_with_f, share_fab_with_f, &
       get_mf_c

contains

  subroutine init_multifab_c (comm) bind(c, name='init_multifab_c')
    integer(c_int), intent(in), value :: comm
    call parallel_initialize(comm)
  end subroutine init_multifab_c

  subroutine build_layout_from_c (nb, dm, lo, hi, pd_lo, pd_hi, pm, pmap) &
       bind(c, name='build_layout_from_c')
    integer(c_int), intent(in), value :: nb, dm
    integer(c_int), intent(in)        :: lo(dm,nb), hi(dm,nb), &
         &                               pd_lo(dm), pd_hi(dm),&
         &                               pm(dm), pmap(nb)

    integer :: i
    logical :: pmask(dm)
    type(box) :: pd
    type(box), allocatable :: bxs(:)
    type(boxarray) :: ba

    ndim = dm

    pmask = (pm /= 0)

    call build(pd, pd_lo, pd_hi)

    allocate(bxs(nb))
    do i = 1, nb
       call build(bxs(i), lo(:,i), hi(:,i)) 
    end do

    call boxarray_build_v(ba, bxs)

    call layout_build_ba(la, ba, pd, pmask, mapping = LA_EXPLICIT, &
         explicit_mapping = pmap)

    call boxarray_destroy(ba)

  end subroutine build_layout_from_c

  subroutine destroy_multifab_c () bind(c, name='destroy_multifab_c')
    integer :: imf, i, n
    n = nlocal(la)
    do imf = 1, nmf
       do i = 1, n
          nullify(mfs(imf)%mf%fbs(i)%p)
       end do
       mfs(nmf)%name = repeat(" ", len_name)
       call multifab_destroy(mfs(imf)%mf)
    end do
    nmf = 0
    call layout_destroy(la)
  end subroutine destroy_multifab_c

  subroutine share_multifab_with_f (mf_name, nc, ng, nd) &
       bind(c, name='share_multifab_with_f')
    character(kind=c_char), intent(in) :: mf_name(*)
    integer(c_int) , intent(in), value :: nc, ng
    integer(c_int), intent(in) :: nd(ndim)

    integer :: i, n
    logical :: nodal(ndim)

    nodal = (nd == 1)

    nmf = nmf + 1
    
    mfs(nmf)%name = repeat(" ", len_name)
    do i = 1, len_name
       if (mf_name(i) == C_NULL_CHAR) exit
       mfs(nmf)%name(i:i) = mf_name(i)
    end do
    
    call multifab_build(mfs(nmf)%mf, la, nc, ng, nodal, fab_alloc=.false.)
    
    n = nlocal(la)
    
    allocate(mfs(nmf)%mf%fbs(n))
    
    do i = 1, n
       call fab_build(mfs(nmf)%mf%fbs(i), &
            get_box(la, global_index(la,i)), &
            nc, ng, nodal, alloc = .false.)
    end do

  end subroutine share_multifab_with_f

  subroutine share_fab_with_f (i, cp) bind(c, name='share_fab_with_f')
    integer(c_int), intent(in), value :: i
    type(c_ptr), intent(in), value :: cp
    
    type(box) :: bx
    integer :: lo(4), hi(4), sz(4), fidx
    double precision, pointer :: fp(:,:,:,:)

    fidx = i + 1 ! Fortran index starts with 1

    bx = get_pbox(mfs(nmf)%mf, fidx)

    lo = 1
    hi = 1
    lo(1:ndim) = lwb(bx)
    hi(1:ndim) = upb(bx)
    hi(4) = ncomp(mfs(nmf)%mf)
    sz = hi - lo + 1

    call c_f_pointer(cp, fp, shape=sz)

    call shift_bound_d4(fp, lo, mfs(nmf)%mf%fbs(fidx)%p)
    ! mfs(nmf)%mf%fbs(i)(lo(1):,lo(2):,lo(3):,lo(4):) => fp  ! some compilers do not support it yet
    
  end subroutine share_fab_with_f

  subroutine shift_bound_d4 (fp, lo, a)
      integer, intent(in) :: lo(4)
      double precision, target, intent(in) :: fp(lo(1):,lo(2):,lo(3):,lo(4):)
      double precision, pointer, intent(inout) :: a(:,:,:,:)
      a => fp
    end subroutine shift_bound_d4

  subroutine get_mf_c (mf, name, ierr)
    type(multifab), intent(inout) :: mf
    character(len=*), intent(in) :: name
    integer, intent(out) :: ierr

    integer :: i

    ierr = 1
    
    do i = 1, nmf
       if (trim(name) .eq. trim(mfs(i)%name)) then
          mf = mfs(i)%mf
          ierr = 0
          exit
       end if
    end do

  end subroutine get_mf_c
  
end module multifab_c_module
