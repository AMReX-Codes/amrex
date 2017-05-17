module cycle_tower_mod

    use constants
    use multifab_module

    implicit none

    type multi_array
        double precision, dimension(:,:,:), allocatable :: arr
    end type multi_array

    type cycle_tower
        integer                       :: dim = 3 ! # dimension of problem
        integer                       :: s1 = 2 ! # of smooths on each grid level, going down
        integer                       :: s2 = 2 ! # of smooths on each grid level, going up
        integer                       :: v_bot = 8 ! # of smooths at the bottom of the V-Cycle
        integer                       :: fmg_bot = 20 ! # of smooths at the bottom of the FMG-Cycle
        integer                       :: nlvl = 0 ! # of multigrid levels
        integer                       :: max_iter = 20 ! # of cycles before automatic program termination
        integer                       :: cycle_type = 1 ! V-Cycle
        integer                       :: interp_type = 1 ! Constant Interpolation
        integer                       :: rhs_type = 1 ! Sine wave
        integer                       :: coeffs_type = 1 ! Constant coefficients
        integer                       :: verbose = 1 ! How much to print
        integer         , allocatable :: n(:) ! Resolution at each grid level
        double precision              :: eps = 1E-10 ! Convergence threshold
        double precision, allocatable :: dx(:) ! Grid spacing, uniform in all directions, at each grid level
        type(multifab)  , allocatable :: uu(:) ! Solution at all grid levels, used for multigrid cycle
        type(multifab)  , allocatable :: ff(:) ! RHS at all grid levels, used for multigrid cycle
        type(multifab)  , allocatable :: bb(:,:) ! Beta coefficients, nodal in direction of second index
        type(multifab)  , allocatable :: res(:) ! Residual at all grid levels
        type(multifab)                :: sol ! Solution to our original problem
        type(multifab)                :: rhs ! RHS of our original problem
        type(multifab)                :: coeffs ! Cell-centered coefficients of our original problem
    end type cycle_tower

contains

   subroutine build_cycle_tower(ct,rhs,coeffs,la,dim,n_cells,dx,s1,s2,v_bot,fmg_bot,nlvl,max_iter,&
                                eps,cycle_type,interp_type,rhs_type,coeffs_type,verbose)

        use multifab_module
        use init_rhs_mod
        use init_coeffs_mod

        implicit none

        integer          , intent(in)    :: dim
        integer          , intent(in)    :: n_cells
        integer          , intent(in)    :: s1
        integer          , intent(in)    :: s2
        integer          , intent(in)    :: v_bot
        integer          , intent(in)    :: fmg_bot
        integer          , intent(in)    :: nlvl
        integer          , intent(in)    :: max_iter
        integer          , intent(in)    :: cycle_type
        integer          , intent(in)    :: interp_type
        integer          , intent(in)    :: rhs_type
        integer          , intent(in)    :: coeffs_type
        integer          , intent(in)    :: verbose
        double precision , intent(in)    :: eps
        double precision , intent(in)    :: dx
        type(multifab)   , intent(in)    :: rhs
        type(multifab)   , intent(in)    :: coeffs
        type(layout)     , intent(inout) :: la
        type(cycle_tower), intent(out)   :: ct
    
        integer                   :: lvl, d, i
        type(layout)              :: la_crse, la_temp

        integer :: lo(dim), hi(dim)
        double precision, pointer :: dp(:,:,:,:)

        ct%dim = dim
        ct%s1 = s1
        ct%s2 = s2
        ct%v_bot = v_bot
        ct%fmg_bot = fmg_bot
        ct%nlvl = nlvl
        ct%max_iter = max_iter
        ct%eps = eps
        ct%cycle_type = cycle_type
        ct%interp_type = interp_type
        ct%rhs_type = rhs_type
        ct%coeffs_type = coeffs_type
        ct%verbose = verbose

        allocate(ct%n(nlvl))
        allocate(ct%dx(nlvl))
        allocate(ct%uu(nlvl))
        allocate(ct%ff(nlvl))
        allocate(ct%res(nlvl))
        allocate(ct%bb(nlvl,dim))

        ct%rhs = rhs
        ct%coeffs = coeffs

        ! Build multifab for final solution, beginning with initial guess of zero
        call multifab_build(ct%sol,la,nc=1,ng=1)
        call setval(ct%sol,0.d0,all=.true.)

        ! Make a copy to be coarsened, so the actual layout isn't lost
        la_temp = la
        do lvl = nlvl,1,-1

            ct%n(lvl) = n_cells/2**(nlvl-lvl)
            ct%dx(lvl) = dx*2**(nlvl-lvl)

            ! Build multifabs at this grid level, to be used later in multigrid cycle
            call multifab_build(ct%uu(lvl),la_temp,nc=1,ng=1) 
            call multifab_build(ct%ff(lvl),la_temp,nc=1,ng=0)
            call multifab_build(ct%res(lvl),la_temp,nc=1,ng=0)

            ! For each direction, build a beta-coefficient multifab that is nodal in only that direction
            do d = 1,dim
                call multifab_build_edge(ct%bb(lvl,d),la_temp,nc=1,ng=0,dir=d)
            end do
           
            ! Build coarsened version of layout, to be used for next multigrid level 
            if (lvl > 1) then
                call layout_build_coarse(la_crse,la_temp,(/(2,lvl=1,ct%dim)/))
            end if

            la_temp = la_crse

        end do
        
        do d = 1,dim
            ! Nodalize coefficients in one direction and fill in appropriate multifab
            call fill_nodal_coeffs(ct%bb(nlvl,d),coeffs,d)
            do lvl = nlvl-1,1,-1
                ! Restrict nodalized coefficients at each multigrid level
                call restrict_coeffs(ct%bb(lvl,d),ct%bb(lvl+1,d),d)
            end do
        end do

        ! Set RHS of top multigrid level to be equal to RHS of original problem
        call copy(ct%ff(nlvl),ct%rhs)

        ! Set solution of top multigrid level to be zero
        call setval(ct%uu(nlvl),0.d0,all=.true.)

    end subroutine build_cycle_tower


    subroutine destroy_cycle_tower(ct)

        use multifab_module

        implicit none

        type(cycle_tower), intent(inout) :: ct

        integer :: lvl, d

        do lvl = 1,ct%nlvl

            call destroy(ct%uu(lvl))
            call destroy(ct%ff(lvl))
            call destroy(ct%res(lvl))

            do d = 1,ct%dim
                call destroy(ct%bb(lvl,d))
            end do

        end do

        call destroy(ct%sol)
        call destroy(ct%rhs)
        call destroy(ct%coeffs)

        deallocate(ct%dx)
        deallocate(ct%n)
        deallocate(ct%uu)
        deallocate(ct%ff)
        deallocate(ct%bb)
        deallocate(ct%res)

    end subroutine destroy_cycle_tower


    subroutine print_cycle_tower(ct,lvl,string)

        use traverse_mod

        implicit none

        integer          , intent(in) :: lvl
        character(len=*) , intent(in) :: string
        type(cycle_tower), intent(inout) :: ct

        logical :: p1,p2,p3,p4,p5
        integer :: n

        call residual(ct%res(lvl),ct%uu(lvl),ct%ff(lvl),ct%bb(lvl,:),ct%dx(lvl))

        if (ct%verbose <= 0) then
            p1 = .false.; p2 = .false.; p3 = .false.; p4 = .false.; p5 = .false.
        else if (ct%verbose == 1) then
            p1 = .false.; p2 = .false.; p3 = .false.; p4 = .false.; p5 = .false.
        else if (ct%verbose == 2) then
            p1 = .true.; p2 = .false.; p3 = .false.; p4 = .false.; p5 = .false.
        else if (ct%verbose == 3) then
            p1 = .true.; p2 = .true.; p3 = .false.; p4 = .false.; p5 = .false.
        else if (ct%verbose == 4) then
            p1 = .true.; p2 = .true.; p3 = .true.; p4 = .false.; p5 = .false.
        else if (ct%verbose == 5) then
            p1 = .true.; p2 = .true.; p3 = .true.; p4 = .true.; p5 = .false.
        else 
            p1 = .true.; p2 = .true.; p3 = .true.; p4 = .true.; p5 = .true.
        end if 

        if (p1 .or. p2 .or. p3 .or. p4) then
            write(*,'(a,a,$)') '   ',string
            if (p1) then
                write(*,'(a,i1,a,i3,a,es26.18)')'InfNorm residual(lvl=',lvl,' n=',ct%n(lvl),') = ',norm_inf(ct%res(lvl))
            end if 
            if (p2) then
                write(*,'(a,i1,a)')'uu(lvl=',lvl,')'
                call print_multifab(ct%uu(lvl),'uu')
            end if
            if (p3) then
                write(*,'(a,i1,a)')'ff(lvl=',lvl,')'
                call print_multifab(ct%ff(lvl),'ff')
            end if
            if (p4) then
                write(*,'(a,i1,a)')'residual(lvl=',lvl,')'
                call print_multifab(ct%res(lvl),'res')
            end if
            if (p5) then
                write(*,'(a,i1,a)')'bbx(lvl=',lvl,')'
                call print_multifab(ct%bb(lvl,1),'bx')
                write(*,'(a,i1,a)')'bby(lvl=',lvl,')'
                call print_multifab(ct%bb(lvl,2),'by')
                if (ct%dim==3) then
                    write(*,'(a,i1,a)')'bbz(lvl=',lvl,')'
                    call print_multifab(ct%bb(lvl,3),'bz')
                end if
            end if
        end if

    end subroutine print_cycle_tower

 
    subroutine print_multifab(mf,string)
 
        use multifab_module

        implicit none

        type(multifab), intent(in) :: mf

        character(len=*)          :: string
        integer                   :: i, j, k, f
        integer                   :: lo(mf%dim),  hi(mf%dim)
        double precision, pointer :: dp(:,:,:,:)

        do f = 1,nfabs(mf)
            lo = lwb(get_pbox(mf,f))
            hi = upb(get_pbox(mf,f))
            write(*,'(a,a,i2,a)') string, ' fab ',f, ':'
            if (mf%dim==3) then
                dp => dataptr(mf,f)
                do k = lo(3),hi(3)
                    write(*,'(a,i2)') 'k =', k
                    do j = hi(2),lo(2),-1
                        print *, j, dp(:,j,k,1)
                    end do 
                    write(*,'(a,$)') '           ' 
                    do i = lo(1),hi(1)
                        write(*,'(a,i0,a,$)') '            ', i, '            '
                    end do
                    print *, ''
                end do
            else if (mf%dim==2) then
                    dp => dataptr(mf,f)
                    do j = hi(2),lo(2),-1
                        write(*,'(i2,$)') j
                        do i = lo(1),hi(1)
                            write(*,'(es26.16,$)') dp(i,j,1,1)
                        end do
                        print *, ''
                    end do 
                    write(*,'(a,$)') '  ' 
                    do i = lo(1),hi(1)
                        write(*,'(a,i0,a,$)') '             ', i, '            '
                    end do
                    print *, ''
            end if
        end do
        print *, ''
     
    end subroutine print_multifab


end module cycle_tower_mod

