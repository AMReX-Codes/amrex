module ml_layout_module

  use layout_module
  use multifab_module
  use ml_boxarray_module
  use bl_prof_module

  implicit none

  type ml_layout
     integer                  :: dim    = 0
     integer                  :: nlevel = 0
     type(ml_boxarray)        :: mba
     type(layout)   , pointer ::    la(:) => Null()
     type(lmultifab), pointer ::  mask(:) => Null() ! cell-centered mask
     logical        , pointer :: pmask(:) => Null() ! periodic mask
  end type ml_layout

  ! wz: this is my note for myself. Will make this clearer when finishing implementation

  ! 0: do nothing
  ! 1: sfc on each level;  ignore fine when distribute;  keep sfc order
  ! 2: do sfc on the finest level first and keep its order; work our way
  !    down; mark coarse grids with fine proc id; honor them if not over
  !    volpercpu; do rest with sfc
  ! 3: work our way up; mark fine grids with coarse proc id.  Do sfc and
  !    cut into chunks.  Let the ones that can benefit most pick first.  Then
  !    let the ones with most works pick.  Try to think how to minimize mpi
  !    gather.
  integer, private, save :: ml_layout_strategy = 1

  interface build
     module procedure ml_layout_build
     module procedure ml_layout_build_n
     module procedure ml_layout_build_mla
  end interface

  interface destroy
     module procedure ml_layout_destroy
  end interface

  interface operator(.eq.)
     module procedure ml_layout_equal
  end interface
  interface operator(.ne.)
     module procedure ml_layout_not_equal
  end interface

  interface print
     module procedure ml_layout_print
  end interface

  interface nlevels
     module procedure ml_layout_nlevels
  end interface

  interface nboxes
     module procedure ml_layout_nboxes
  end interface

  interface get_box
     module procedure ml_layout_get_box
  end interface

  interface built_q
     module procedure ml_layout_built_q
  end interface

contains

  subroutine ml_layout_set_strategy(i)
    integer, intent(in) :: i
    ml_layout_strategy = i
  end subroutine ml_layout_set_strategy

  function ml_layout_built_q(mla) result(r)
    logical :: r
    type(ml_layout), intent(in) :: mla
    r = associated(mla%la)
  end function ml_layout_built_q

  function ml_layout_nlevels(mla) result(r)
    integer :: r
    type(ml_layout), intent(in) :: mla
    r = mla%nlevel
  end function ml_layout_nlevels

  function ml_layout_nboxes(mla, lev) result(r)
    integer :: r
    type(ml_layout), intent(in) :: mla
    integer, intent(in) :: lev
    r = nboxes(mla%mba, lev)
  end function ml_layout_nboxes

  function ml_layout_equal(mla1, mla2) result(r)
    logical :: r
    type(ml_layout), intent(in) :: mla1, mla2
    r = associated(mla1%la, mla2%la)
  end function ml_layout_equal
  
  function ml_layout_not_equal(mla1, mla2) result(r)
    logical :: r
    type(ml_layout), intent(in) :: mla1, mla2
    r = .not. associated(mla1%la, mla2%la)
  end function ml_layout_not_equal

  function ml_layout_get_layout(mla, n) result(r)
    type(layout) :: r
    type(ml_layout), intent(in) :: mla
    integer, intent(in) :: n
    r = mla%la(n)
  end function ml_layout_get_layout

  function ml_layout_get_pd(mla, n) result(r)
    type(box) :: r
    type(ml_layout), intent(in) :: mla
    integer, intent(in) :: n
    r = ml_boxarray_get_pd(mla%mba, n)
  end function ml_layout_get_pd

  function ml_layout_get_box(mla, lev, n) result(r)
    type(box) :: r
    type(ml_layout), intent(in) :: mla
    integer, intent(in) :: n, lev
    r = get_box(mla%la(lev), n)
  end function ml_layout_get_box

  subroutine ml_layout_build_n(mla, nlevel, dm)
    type(ml_layout), intent(out) :: mla
    integer, intent(in) :: nlevel, dm

    mla%nlevel              = nlevel
    mla%dim                 = dm
    allocate(mla%pmask(mla%dim))
    allocate(mla%la(mla%nlevel))
    allocate(mla%mask(mla%nlevel-1))
    call build(mla%mba, nlevel, dm)
  end subroutine ml_layout_build_n

  ! The behavior of this subroutine has changed!!!
  ! The layouts in mla are no longer simple copies of la_array.
  ! They might be a simple copy, or a different one built on the same boxarray.
  ! It is now caller's responsibility to check and delete the layouts in la_array
  ! that are not used in mla.  We cannot do it for the caller because there
  ! could be multifabs that are still using those layouts.  
  ! See Src/F_BaseLib/regrid.f90 for examples.
  subroutine ml_layout_build_la_array(mla, la_array, mba, pmask, nlevel)

    type(ml_layout  ), intent(  out) :: mla
    type(   layout  ), intent(inout) :: la_array(:)
    type(ml_boxarray), intent(in   ) :: mba
    integer,           intent(in   ) :: nlevel
    logical                        :: pmask(:)

    type(boxarray) :: bac
    integer        :: n

    mla%nlevel = nlevel
    mla%dim    = get_dim(mba)

    ! Copy only nlevel levels of the mba
    call build(mla%mba,nlevel,mla%dim)

    mla%mba%pd(1:nlevel) = mba%pd(1:nlevel)
    do n = 1, mla%nlevel-1
      mla%mba%rr(n,:) = mba%rr(n,:)
    end do

    do n = 1, mla%nlevel
      call copy(mla%mba%bas(n),mba%bas(n))
    end do

    ! Build the pmask
    allocate(mla%pmask(mla%dim))
    mla%pmask  = pmask

    ! Point to the existing la_array(:)
    allocate(mla%la(mla%nlevel))
    
    call optimize_layouts(mla%la, la_array, mla%nlevel, mla%mba%rr)

    allocate(mla%mask(mla%nlevel-1))

    do n = mla%nlevel-1,  1, -1
       call lmultifab_build(mla%mask(n), mla%la(n), nc = 1, ng = 0)
       call setval(mla%mask(n), val = .TRUE.)
       call copy(bac, mba%bas(n+1))
       call boxarray_coarsen(bac, mba%rr(n,:))
       call setval(mla%mask(n), .false., bac)
       call destroy(bac)
    end do

  end subroutine ml_layout_build_la_array

  subroutine ml_layout_build_mla(mla, mla_in)
    type(ml_layout), intent(inout) :: mla
    type(ml_layout), intent(in   ) :: mla_in

    integer :: n

    mla%dim    = mla_in%dim
    mla%nlevel = mla_in%nlevel

    allocate(mla%pmask(mla%dim))
    mla%pmask = mla_in%pmask

    call copy(mla%mba, mla_in%mba)

    allocate(mla%la(mla%nlevel))
    allocate(mla%mask(mla%nlevel-1))
    do n = 1, mla%nlevel
       call build(mla%la(n), mla%mba%bas(n), mla%mba%pd(n), pmask=mla%pmask, &
            explicit_mapping=get_proc(mla_in%la(n)))
    end do
    do n = 1, mla%nlevel-1
       call lmultifab_build(mla%mask(n), mla%la(n), nc = 1, ng = 0)
       call lmultifab_copy(mla%mask(n), mla_in%mask(n))
    end do

  end subroutine ml_layout_build_mla

  subroutine ml_layout_build(mla, mba, pmask)
    type(ml_layout)  , intent(inout) :: mla
    type(ml_boxarray), intent(in   ) :: mba
    logical, optional                :: pmask(:)
    call ml_layout_restricted_build(mla, mba, mba%nlevel, pmask)
  end subroutine ml_layout_build

  subroutine ml_layout_restricted_build(mla, mba, nlevs, pmask)

    ! this subroutine is the same thing as ml_layout_build except that
    ! the mla will only have nlevs instead of mba%nlevel

    type(ml_layout)  , intent(inout) :: mla
    type(ml_boxarray), intent(in   ) :: mba
    integer          , intent(in   ) :: nlevs
    logical, optional                :: pmask(:)

    type(boxarray) :: bac
    type(layout), allocatable :: la_array(:)
    integer :: n
    logical :: lpmask(mba%dim)

    lpmask = .false.; if (present(pmask)) lpmask = pmask
    allocate(mla%pmask(mba%dim))
    mla%pmask  = lpmask

    mla%nlevel = nlevs
    mla%dim    = mba%dim

!   Have to copy only nlevs of the mba
!   Replace 
!   call copy(mla%mba, mba)
!   by these lines
    call build(mla%mba,nlevs,mla%dim)
    mla%mba%pd(1:nlevs) = mba%pd(1:nlevs)
    do n = 1, mla%nlevel-1
      mla%mba%rr(n,:) = mba%rr(n,:)
    end do
    do n = 1, mla%nlevel
      call copy(mla%mba%bas(n),mba%bas(n))
    end do

    allocate(mla%la(mla%nlevel), la_array(mla%nlevel))

    do n = 1, mla%nlevel
       call build(la_array(n), mba%bas(n), mba%pd(n), pmask=lpmask)
    end do

    call optimize_layouts(mla%la, la_array, mla%nlevel, mba%rr)

    do n = 1, mla%nlevel
       if (mla%la(n) .ne. la_array(n)) then
          call destroy(la_array(n))
       end if
    end do

    allocate(mla%mask(mla%nlevel-1))

    do n = mla%nlevel-1,  1, -1
       call lmultifab_build(mla%mask(n), mla%la(n), nc = 1, ng = 0)
       call setval(mla%mask(n), val = .TRUE.)
       call copy(bac, mba%bas(n+1))
       call boxarray_coarsen(bac, mba%rr(n,:))
       call setval(mla%mask(n), .false., bac)
       call destroy(bac)
    end do

  end subroutine ml_layout_restricted_build

  subroutine ml_layout_destroy(mla, keep_coarse_layout)
    type(ml_layout), intent(inout) :: mla
    logical, intent(in), optional :: keep_coarse_layout
    integer :: n, n0
    logical :: lkeepcoarse

    lkeepcoarse = .false.;  if (present(keep_coarse_layout)) lkeepcoarse = keep_coarse_layout

    do n = 1, mla%nlevel-1
       if (built_q(mla%mask(n))) call destroy(mla%mask(n))
    end do
    call destroy(mla%mba)

    if (lkeepcoarse) then
       n0 = 2
    else
       n0 = 1
    end if
    do n = n0, mla%nlevel
       call destroy(mla%la(n))
    end do

    deallocate(mla%la, mla%mask)
    mla%dim = 0
    mla%nlevel = 0
    deallocate(mla%pmask)
  end subroutine ml_layout_destroy

  subroutine ml_layout_print(mla, str, unit, skip)
    use bl_IO_module
    type(ml_layout), intent(in) :: mla
    character (len=*), intent(in), optional :: str
    integer, intent(in), optional :: unit
    integer, intent(in), optional :: skip
    integer :: i, j
    integer :: un
    un = unit_stdout(unit)
    call unit_skip(un, skip)
    write(unit=un, fmt = '("MLLAYOUT[(*")', advance = 'no')
    if ( present(str) ) then
       write(unit=un, fmt='(" ",A)') str
    else
       write(unit=un, fmt='()')
    end if
    call unit_skip(un, skip)
    write(unit=un, fmt='(" DIM     = ",i2)') mla%dim
    call unit_skip(un, skip)
    write(unit=un, fmt='(" NLEVEL  = ",i2)') mla%nlevel
    call unit_skip(un, skip)
    write(unit=un, fmt='(" *) {")')
    do i = 1, mla%nlevel
       call unit_skip(un, unit_get_skip(skip)+1)
       write(unit=un, fmt = '("(* LEVEL ", i2)') i
       call unit_skip(un, unit_get_skip(skip)+1)
       write(unit=un, fmt = '(" PD = ")', advance = 'no')
       call print(mla%mba%pd(i), unit=un, advance = 'NO')
       write(unit=un, fmt = '(" *) {")')
       do j = 1, nboxes(mla%mba%bas(i))
           call unit_skip(un, unit_get_skip(skip)+2)
           write(unit=un, fmt = '("{")', advance = 'no')
           call print(get_box(mla%mba%bas(i),j), unit = unit, advance = 'NO')
           write(unit=un, fmt = '(", ", I0, "}")', advance = 'no') get_proc(mla%la(i), j)
           if ( j == nboxes(mla%mba%bas(i)) ) then
              call unit_skip(un, unit_get_skip(skip)+1)
              write(unit=un, fmt = '("}")')
           else
              write(unit=un, fmt = '(",")')
           end if
       end do
       if ( i == mla%nlevel ) then
          call unit_skip(un, skip)
          write(unit=un, fmt = '("}]")')
       else
          write(unit=un, fmt = '(",")')
       end if
    end do
  end subroutine ml_layout_print

  subroutine optimize_layouts(la_new, la_old, nlevs, rr)
    use fab_module
    use knapsack_module
    type(layout), intent(out  ) :: la_new(:)
    type(layout), intent(inout) :: la_old(:)
    integer, intent(in) :: nlevs, rr(:,:)

    integer :: n
    logical :: mc_flag
    integer(kind=ll_t), allocatable :: lucvol(:)
    type(bl_prof_timer), save :: bpt

    if (ml_layout_strategy .eq. 0 .or. nlevs.eq.1) then
       la_new(1:nlevs) = la_old(1:nlevs)
       return
    end if

    call build(bpt, "optimize_layouts")       

    mc_flag = get_manual_control_least_used_cpus_flag()
    call set_manual_control_least_used_cpus(.true.)

    allocate(lucvol(nlevs))
    do n = 1, nlevs
       lucvol(n) = layout_local_volume(la_old(n))         
    end do

    select case (ml_layout_strategy)
    case (1)
       call layout_opt_ignore_fine()
    case (2)
       call layout_opt_from_fine_to_coarse()
    case default
       call bl_error("unknown ml_layout_strategy", ml_layout_strategy)
    end select

    deallocate(lucvol)
    
    call set_manual_control_least_used_cpus(mc_flag)

    call destroy(bpt)

  contains

    subroutine build_new_layout(lao, lai)
      type(layout), intent(out) :: lao
      type(layout), intent(in ) :: lai

      logical :: natural_order
      integer :: nbxs, nprocs, i
      integer, pointer     :: luc(:)
      integer, allocatable :: sfc_order(:), prc(:), ibxs(:)
      type(boxarray) :: ba
      type(box), pointer :: pbxs(:)

      nbxs = nboxes(lai)
      nprocs = parallel_nprocs()
      
      ba = get_boxarray(lai)

      allocate(sfc_order(nbxs), prc(nbxs), ibxs(nbxs))

      if (sfc_order_built_q(lai)) then
         sfc_order = get_sfc_order(lai)
      else
         pbxs => dataptr(ba)
         call make_sfc(sfc_order, pbxs)
         nullify(pbxs)
      end if

      luc => least_used_cpus()  ! 0 based

      natural_order = .true.
      do i = 0, nprocs-1
         if (luc(i) .ne. i) then
            natural_order = .false.
            exit
         end if
      end do

      do i = 1, nbxs
         ibxs(i) = volume(get_box(lai,i))
      end do
      
      ! If the processors are not in natural order, let's sort sfc chunks in reverse order
      call distribute_sfc(prc, sfc_order, ibxs, nprocs, sort_chunks = .not.natural_order)
      
      do i = 1, nbxs
         prc(i) = luc(prc(i))
      end do
      
      deallocate(luc)

      if (all(prc .eq. get_proc(lai))) then
         lao = lai
      else
         call layout_build_ba(lao, ba, get_pd(lai), get_pmask(lai), explicit_mapping=prc)
      end if

      call set_sfc_order(lao, sfc_order)

    end subroutine build_new_layout


    subroutine layout_opt_ignore_fine() ! aka strategy 1
      ! build the coarest level first
      call set_luc_vol(0_ll_t) 
      call build_new_layout(la_new(1), la_old(1))
      do n = 2, nlevs
         ! ignore the finer levels in least_used_cpus() in fab.f90
         ! Note that layout is affected by the results of least_used_cpus()
         call set_luc_vol(sum(lucvol(1:n-1)))
         call build_new_layout(la_new(n), la_old(n))
      end do
    end subroutine layout_opt_ignore_fine


    subroutine layout_opt_from_fine_to_coarse() ! aka strategy 2      

      use vector_i_module
      use sort_i_module

      logical :: overlap
      integer :: i, j, k, m, lo(4), hi(4), iproc, ntmp
      integer :: myproc, nprocs, nbxs, nlft
      integer, allocatable :: prc(:), prc2(:), cnt(:), vbx(:), itmp(:), novlp(:), novlp2(:),&
           novlptmp(:), vbxlft(:), prclft(:)
      integer, pointer :: p(:,:,:,:)
      integer (kind=ll_t) :: vtot, maxvol
      integer (kind=ll_t), allocatable :: vprc(:)
      real(kind=dp_t) :: volpercpu
      type(vector_i), allocatable :: bxid(:)
      type(boxarray) :: bac
      type(layout) :: lacfine
      type(imultifab) :: mff, mfc

      myproc = parallel_myproc()
      nprocs = parallel_nprocs()

      allocate(cnt(nprocs), vprc(nprocs))
      allocate(bxid(nprocs))

      do iproc = 1,nprocs
         call build(bxid(iproc))
      end do

      ! First, let build the finest layout
      call set_luc_vol(0_ll_t)
      call build_new_layout(la_new(nlevs), la_old(nlevs))

      do n = nlevs-1, 1, -1

         call set_luc_vol(sum(lucvol(n+1:nlevs)))

         nbxs = nboxes(la_old(n))
         allocate(prc(nbxs), prc2(nbxs), vbx(nbxs))
         allocate(novlp(nbxs), novlp2(nbxs))

         call layout_build_coarse(lacfine, la_new(n+1), rr(n,:))

         call imultifab_build(mff, lacfine  , 1, 0)
         call imultifab_build(mfc, la_old(n), 1, 0)

         do m = 1, nfabs(mff)
            call setval(mff%fbs(m), myproc+1)  ! make it 1-based
         end do

         call setval(mfc, 0)
         call imultifab_copy(mfc, mff)

         prc2   = 0
         novlp2 = 0

         do m = 1, nfabs(mfc)

            p => dataptr(mfc, m)
            lo = lbound(p)
            hi = ubound(p)

            cnt = 0
            overlap = .false.
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     if (p(i,j,k,1) .gt. 0) then
                        overlap = .true.
                        cnt(p(i,j,k,1)) = cnt(p(i,j,k,1)) + 1
                     end if
                  end do
               end do
            end do

            if (overlap) then
               i = layout_global_index(la_old(n),m)
               prc2  (i) = maxloc(cnt,1)
               novlp2(i) = maxval(cnt)
            end if

         end do

         call destroy(mff)
         call destroy(mfc)

         call parallel_reduce(prc  ,   prc2, MPI_SUM)
         call parallel_reduce(novlp, novlp2, MPI_SUM)
         ! prc contains the preferred proc. no. for each box
         ! novlp contains the number of points overlaping with fine level on the preferred proc. 
         deallocate(prc2,novlp2)

         do iproc=1,nprocs
            call clear(bxid(iproc))
         end do

         bac = get_boxarray(la_old(n))

         vtot = 0_ll_t
         vprc = 0_ll_t
         do m = 1, nbxs
            vbx(m) = volume(get_box(bac,m))
            vtot = vtot + vbx(m)
            iproc = prc(m)
            if (iproc .gt. 0) then
               vprc(iproc) = vprc(iproc) + vbx(m)
               call push_back(bxid(iproc), m)
            end if
         end do
         volpercpu = real(vtot,dp_t) / real(nprocs,dp_t)

         do iproc=1,nprocs
            if (vprc(iproc) .gt. volpercpu) then
               ntmp = size(bxid(iproc))
               if (ntmp .gt. 1) then
                  allocate(novlptmp(ntmp), itmp(ntmp))
                  do j = 1, ntmp
                     novlptmp(j) = novlp(at(bxid(iproc),j))
                  end do
                  call sort(novlptmp,itmp)
                  do j = 1, ntmp-1
                     m = at(bxid(iproc),itmp(j))
                     prc(m) = 0
                     vprc(iproc) = vprc(iproc) - vbx(m)
                     if (vprc(iproc) .le. volpercpu) exit
                  end do
               end if
               deallocate(novlptmp,itmp)
            end if
            call clear(bxid(iproc))
         end do

         nlft = count(prc .eq. 0)

         if (nlft .gt. 0) then

            allocate(vbxlft(nlft),prclft(nlft))
            vbxlft = 0
            
            j = 1
            do m = 1, nbxs
               if (prc(m) .eq. 0) then
                  vbxlft(j) = vbx(m)
                  j = j + 1
               end if
            end do

            call knapsack_pf(prclft, vbxlft, nprocs, vprc, verbose=.false.)

            j = 1
            do m = 1, nbxs
               if (prc(m) .eq. 0) then
                  prc(m) = prclft(j) + 1
                  j = j+1
               end if
            end do

            if (.true.) then
               if (parallel_ioprocessor()) then
                  vprc = 0_ll_t
                  do m = 1, nbxs
                     vprc(prc(m)) =  vprc(prc(m)) + vbx(m)
                  end do
                  maxvol = maxval(vprc)
               end if
            end if

            deallocate(vbxlft,prclft)
         end if

         prc = prc - 1
         call layout_build_ba(la_new(n), bac, get_pd(la_old(n)), get_pmask(la_old(n)), &
              explicit_mapping=prc)

         deallocate(prc,vbx,novlp)
      end do

      deallocate(cnt,vprc)
      do iproc = 1, nprocs
         call destroy(bxid(iproc))
      end do
      deallocate(bxid)

    end subroutine layout_opt_from_fine_to_coarse

  end subroutine optimize_layouts

end module ml_layout_module
