! Take 2 single-level plotfiles as input and compare them in finite-difference way.
! 
! The two plotfiles must have the same physical domain.  They do not need to
! have the same number of grids.  But they must have some kind of coarse-fine relation.
!

program fcoarsen

  use BoxLib
  use parallel
  use bl_space
  use bl_error_module
  use bl_IO_module
  use fabio_module
  use layout_module
  use multifab_module
  use plotfile_module

  implicit none

  character (len=256) :: plotfile_f, plotfile_c, input_field
  integer :: narg, farg
  character (len=256) :: fname

  type(multifab), pointer :: pltdata_f(:)
  type(multifab) :: mf_c, mf_d
  type(layout) :: la_f, la_c

  type(box) :: pb_f
  type(boxarray) :: ba_f, ba_c
  integer :: rr(3)

  integer :: i,j,k,m,n, lo(3), hi(3)
  double precision, pointer, dimension(:,:,:,:) :: fp, cp

  type(plotfile) :: pf_f
  character (len=20), pointer :: vnames(:) => Null()
  character (len=30) :: afoo

  double precision :: prob_lo(3), prob_hi(3), time, dx(3)
  type(multifab) :: err(1)
  integer :: rr0(0)

  integer :: dm

  !---------------------------------------------------------------------------
  ! process the command line arguments

  narg = command_argument_count()

  ! defaults
  plotfile_f = ""
  plotfile_c = ""
  input_field = ""

  ! use coarsening factor of 2 by default
  rr(:) = 2

  farg = 1
  do while (farg <= narg)
     call get_command_argument(farg, value = fname)
     
     select case (fname)

     case ('-i', '--infile')
        farg = farg + 1
        call get_command_argument(farg, value = plotfile_f)

     case ('-e', '--outfile')
        farg = farg + 1
        call get_command_argument(farg, value = plotfile_c)

     case ('-c', '--coarsen')
        farg = farg + 1
        call get_command_argument(farg, value = input_field)
        read(input_field,*) rr(1)
        rr(2:)=rr(1)

     case default
        exit

     end select
     farg = farg + 1
  enddo

  if (len_trim(plotfile_f) == 0) then
     print *, " "
     print *, "Take a plotfile and write out the coarsened plotfile."
     print *, "usage (all arguments are REQUIRED):"
     print *, "   fcoarsen --infile file1 --outfile file2"
     print *, " "
     stop
  endif

  call boxlib_initialize()

  ! get variable names, prob_lo/hi, time
  call build(pf_f, plotfile_f, unit_new())
  allocate(vnames(pf_f%nvars))
  do m=1,pf_f%nvars
     vnames(m) = trim(pf_f%names(m))
  end do
  prob_lo = pf_f%plo
  prob_hi = pf_f%phi
  time = pf_f%tm
  call destroy(pf_f)

  ! read plotfile
  call fabio_ml_multifab_read_d(pltdata_f, plotfile_f)

  if (size(pltdata_f) .ne. 1) then
     call bl_error("ERROR: plotfile has more than one level")
  end if

  la_f = get_layout(pltdata_f(1))
  pb_f = layout_get_pd(la_f)

  dm = layout_dim(la_f)

  ba_f = layout_boxarray(la_f)
  call boxarray_build_copy(ba_c, ba_f)
  call boxarray_coarsen(ba_c, rr)

  call print(ba_c)

  call layout_build_ba(la_c,ba_c,boxarray_bbox(ba_c),pmask=get_pmask(la_f))
  call multifab_build(mf_c, la_c, pltdata_f(1)%nc, 0)

  do n=1,nfabs(mf_c)

     fp => dataptr(pltdata_f(1), n)
     cp => dataptr(mf_c, n)

     lo = lwb(get_box(mf_c,n))
     hi = upb(get_box(mf_c,n))

     if (dm .eq. 2) then

        do m=1,mf_c%nc
           do j = lo(2),hi(2)
              do i = lo(1),hi(1)
                 cp(i,j,1,m) = sum( fp(i*rr(1):i*rr(1)+rr(1)-1, &
                                       j*rr(2):j*rr(2)+rr(2)-1, 1,m) ) / product(rr(1:dm))                 
              end do
           end do
        end do

     else if (dm .eq. 3) then

        do m=1,mf_c%nc
           do k = lo(3),hi(3)
              do j = lo(2),hi(2)
                 do i = lo(1),hi(1)
                    cp(i,j,k,m) = sum( fp(i*rr(1):i*rr(1)+rr(1)-1, &
                                          j*rr(2):j*rr(2)+rr(2)-1, &
                                          k*rr(3):k*rr(3)+rr(3)-1, m) ) / product(rr(1:dm))  
                 end do
              end do
           end do
        end do

     end if

  end do

  if (len_trim(plotfile_c) .ne. 0) then
     err(1) = mf_c
     dx = (prob_hi - prob_lo) / box_extent(boxarray_bbox(layout_boxarray(la_f)))

     call fabio_ml_multifab_write_d(err, rr0, plotfile_c, vnames, &
                                    la_c%lap%pd, prob_lo, prob_hi, time, dx(1:dm)*rr(1:dm))
  end if

  call destroy(ba_c)
  call destroy(mf_c)
  call destroy(pltdata_f(1))
  call destroy(la_c)
  call destroy(la_f)

  deallocate(pltdata_f)
  deallocate(vnames)

end program fcoarsen
