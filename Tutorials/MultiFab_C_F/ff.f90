module ff_module

  implicit none

  public :: ff

contains

  subroutine ff () bind(c, name='ff')
    use multifab_module
    use multifab_c_module
    use bl_types
    
    integer :: ierr, nc, ng
    double precision :: norm
    type(multifab) :: mf

    if (parallel_IOProcessor()) then
       print *, ''
       print *, "In Fortran"
    end if

    ! Assuming there is a mf "data" in multifab_c_module,
    ! which has been shared by C++.
    call get_mf_c (mf, "data", ierr)
    if (ierr .eq. 0) then
       if (parallel_IOProcessor()) then
          print *, '    Successfully get data mf!'
       end if

       nc = multifab_ncomp(mf)
       ng = multifab_nghost(mf)

       if (parallel_IOProcessor()) then
          print *, '    ncomp  = ', nc
          print *, '    nghost = ', ng
       end if

       !-----------------------------------------

       call multifab_setval(mf, 1.d0, all=.true.)

       norm = norm_l1(mf, nc, 1, all=.true.)

       if (parallel_IOProcessor()) then
          print *, "    1-norm of last component", norm
       end if

       !-----------------------------------------

       call multifab_setval(mf, 3.d0, all=.false.)

       norm = norm_l1(mf, nc, 1, all=.true.)

       if (parallel_IOProcessor()) then
          print *, "    1-norm of last component", norm
       end if

       !-----------------------------------------

       call multifab_fill_boundary(mf)

       norm = norm_l1(mf, nc, 1, all=.true.)

       if (parallel_IOProcessor()) then
          print *, "    1-norm of last component", norm
       end if

    end if

    call get_mf_c (mf, "foo", ierr)
    if (ierr .ne. 0) then
       if (parallel_IOProcessor()) then
          print *, '    Cannot get foo mf!'
       end if
    end if

    ! no need to destroy mf!
  end subroutine ff

end module ff_module
