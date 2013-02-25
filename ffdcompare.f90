! Take 2 single-level plotfiles as input and compare them in finite-difference way.
! 
! The two plotfiles must have the same physical domain.  They do not need to
! have the same number of grids.  But they must have some kind of coarse-fine relation.
!

program ffdcompare

  use f2kcli
  use bl_space
  use bl_error_module
  use bl_IO_module
  use fabio_module
  use layout_module
  use multifab_module
  use plotfile_module

  implicit none

  character (len=256) :: plotfile_a, plotfile_b
  integer :: narg, farg
  character (len=256) :: fname

  type(multifab), pointer :: pltdata_a(:), pltdata_b(:)
  type(multifab) :: mf_c, mf_d
  type(layout) :: la_a, la_b, la_c

  type(box) :: pb_a, pb_b
  type(boxarray) :: ba_b, ba_c
  integer :: rr

  integer :: i,j,k,m,n, lo(3), hi(3), ncells, pblo(3)
  double precision, pointer, dimension(:,:,:,:) :: bp, cp

  type(plotfile) :: pf_a
  character (len=256), pointer :: vnames(:) => Null()
  character (len=30) :: afoo

  !---------------------------------------------------------------------------
  ! process the command line arguments

  narg = command_argument_count()

  ! defaults
  plotfile_a = ""
  plotfile_b = ""

  farg = 1
  do while (farg <= narg)
     call get_command_argument(farg, value = fname)
     
     select case (fname)

     case ('--infile1')
        farg = farg + 1
        call get_command_argument(farg, value = plotfile_a)

     case ('--infile2')
        farg = farg + 1
        call get_command_argument(farg, value = plotfile_b)

     case default
        exit

     end select
     farg = farg + 1
  enddo

  if (len_trim(plotfile_a) == 0 .OR. len_trim(plotfile_b) == 0) then
     print *, " "
     print *, "Take 2 single-level plotfiles as input and compare them in finite-difference way."
     print *, "The two plotfiles must have the same physical domain.  They do not need to"
     print *, "have the same number of grids.  But they must have some kind of coarse-fine relation. "
     print *, " "
     print *, "usage:"
     print *, "   ffdcompare --infile1 file1 --infile2 file2 [--norm 0]"
     print *, " "
     stop
  endif

  ! get variable names
  call build(pf_a, plotfile_a, unit_new())
  allocate(vnames(pf_a%nvars))
  do m=1,pf_a%nvars
     vnames(m) = trim(pf_a%names(m))
  end do
  call destroy(pf_a)

  ! read plotfiles
  call fabio_ml_multifab_read_d(pltdata_a, plotfile_a)
  call fabio_ml_multifab_read_d(pltdata_b, plotfile_b)

  if (size(pltdata_a) .ne. 1 .or. size(pltdata_b) .ne. 1) then
     call bl_error("ERROR: plotfiles have more than one level")
  end if

  la_a = get_layout(pltdata_a(1))
  la_b = get_layout(pltdata_b(1))

  pb_a = layout_get_pd(la_a)
  pb_b = layout_get_pd(la_b)

  pblo = pb_a%lo

  rr = (pb_b%hi(1) - pb_b%lo(1) + 1) / (pb_a%hi(1) - pb_a%lo(1) + 1)
  
  ba_b = layout_boxarray(la_b)
  call boxarray_build_copy(ba_c, ba_b)
  call boxarray_coarsen(ba_c, rr)

  call layout_build_ba(la_c,ba_c,boxarray_bbox(ba_c),pmask=get_pmask(la_b))
  call multifab_build(mf_c, la_c, pltdata_b(1)%nc, 0)
  call destroy(ba_c)

  do n=1,nfabs(mf_c)

     bp => dataptr(pltdata_b(1), n)
     cp => dataptr(mf_c, n)

     lo = lwb(get_box(mf_c,n))
     hi = upb(get_box(mf_c,n))

     do m=1,mf_c%nc
        do k = lo(3),hi(3)
           do j = lo(2),hi(2)
              do i = lo(1),hi(1)
                 cp(i,j,k,m) = bp((i-pblo(1))*rr+pblo(1), &
                      &           (j-pblo(2))*rr+pblo(2), &
                      &           (k-pblo(3))*rr+pblo(3), m)
                end do
             end do
          end do
     end do
  end do

  call destroy(pltdata_b(1))
  call destroy(la_b)

  call multifab_build(mf_d, la_a, mf_c%nc, 0)
  call multifab_copy_c(mf_d,1,mf_c,1,mf_c%nc)

  call destroy(mf_c)
  call destroy(la_c)

  call multifab_sub_sub(mf_d, pltdata_a(1))

  ncells = volume(get_boxarray(mf_d))

1000 format(1x,a24,2x,g20.10,2x,g20.10,2x,g20.10)
1001 format(1x,a24,2x,a20,2x,a20,2x,a20)

  print *, ""

  afoo = "variable" 
  write (*,1001) afoo, "0-norm", "1-norm", "2-norm"
  write (*,*) "------------------------------------------------------------------------------------------"
  do m=1, mf_d%nc
     write (*,1000) vnames(m), multifab_norm_inf_c(mf_d, m), &
          multifab_norm_l1_c(mf_d, m)/dble(ncells), &
          multifab_norm_l2_c(mf_d, m)/sqrt(dble(ncells)) 
  end do

  call destroy(mf_d)
  call destroy(pltdata_a(1))
  call destroy(la_a)

  deallocate(pltdata_a)
  deallocate(pltdata_b)
  deallocate(vnames)

end program ffdcompare
