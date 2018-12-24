! print out information regarding the AMR levels and boxes

program fboxinfo

  use bl_space
  use bl_error_module
  use bl_constants_module
  use bl_IO_module
  use plotfile_module

  implicit none

  type(plotfile) pf
  integer :: unit
  integer :: i, j, f

  integer :: ncells, ncells_domain

  integer :: lo(MAX_SPACEDIM), hi(MAX_SPACEDIM)
  integer :: dim

  integer :: narg, farg
  character(len=256) :: fname

  logical :: full, gridfile, castro, levels

  unit = unit_new()

  full = .false.
  gridfile = .false.
  castro = .false.
  levels = .false.

  narg = command_argument_count()

  farg = 1
  do while ( farg <= narg )
     call get_command_argument(farg, value = fname)

     select case (fname)

     case ('-f', '--full')
        full = .true.

     case ('-g', '--gridfile')
         gridfile = .true.

     case ('-l', '--levels')
         levels = .true.

     case ('-c', '--castro')
         castro = .true.

     case default
        exit

     end select

     farg = farg + 1
  end do

  if (gridfile .and. full) then
     print *, "ERROR: cannot specify both full and gridfile modes"
     stop
  endif
  if (castro .and. full) then
     print *, "ERROR: cannot specify both full and castro modes"
     stop
  endif

  if (farg > narg) then
     print *, " "
     print *, "Dump out information about the AMR levels and boxes"
     print *, "Works with 1-, 2-, or 3-d datasets."
     print *, " "
     print *, "usage:"
     print *, "   fboxinfo [-f|--full] plotfile"
     print *, " "
     print *, "args:"
     print *, "   [-f|--full]     output detailed information about the boxes"
     print *, "   [-g|--gridfile] output a gridfile for use with test_average"
     print *, "   [-c|--castro]   output a gridfile for use with castro restart"
     print *, "   [-l|--levels]   just output the number of levels"
     print *, " "
     stop
  end if


  do f = farg, narg

     call get_command_argument(f, value = fname)
     call build(pf, fname, unit)

     if (.not. gridfile .and. .not. castro .and. .not. levels) then
        print *, "plotfile: ", trim(fname)
     end if

     dim = pf%dim


     ! output the number of boxes at each level
1000 format(1x,"level ", i3, ": number of boxes = ", i6, ", volume = ", f6.2, "%")
1001 format(1x,"      ", 1x, "  maximum zones =   ", i7)
1002 format(1x,"      ", 1x, "  maximum zones =   ", i7, " x ", i7)
1003 format(1x,"      ", 1x, "  maximum zones =   ", i7, " x ", i7, " x " i7)

     if (levels) then
        print *, pf%flevel
        cycle
     endif


     if (.not. gridfile .and. .not. castro) then
        do i = 1, pf%flevel
           ncells = 0
           do j = 1, nboxes(pf, i)
              ncells = ncells + volume(get_pbox(pf, i, j))
           enddo
           ncells_domain = volume(plotfile_get_pd_box(pf, i)) 

           write (*,1000) i, nboxes(pf, i), 100.0*dble(ncells)/ncells_domain

           lo = lwb(plotfile_get_pd_box(pf, i))
           hi = upb(plotfile_get_pd_box(pf, i))

           select case (dim)
           case (1)
              write (*,1001) hi(1)-lo(1)+1
              
           case (2)
              write (*,1002) hi(1)-lo(1)+1, hi(2)-lo(2)+1
              
           case (3)
              write (*,1003) hi(1)-lo(1)+1, hi(2)-lo(2)+1, hi(3)-lo(3)+1
              
           end select
           
           print *, " "
        enddo
     endif

2001 format(1x,"  box ", i5, ":  ("i5,")   ("i5,")")
2002 format(1x,"  box ", i5, ":  ("i5,",",i5,")   ("i5,",",i5,")")
2003 format(1x,"  box ", i5, ":  ("i5,",",i5,",",i5")   ("i5,",",i5,",",i5")")

     ! if desired, output detailed box infomation
     if (full) then

        do i = 1, pf%flevel
           
           print *, " "
           print *, " level ", i
              
           do j = 1, nboxes(pf, i)
              
              lo = lwb(get_box(pf, i, j))
              hi = upb(get_box(pf, i, j))
              
              select case (dim)
              case (1)
                 write (*,2001) j, lo(1), hi(1)
                 
              case (2)
                 write (*,2002) j, lo(1), lo(2), hi(1), hi(2)
                 
              case (3)
                 write (*,2003) j, lo(1), lo(2), lo(3), hi(1), hi(2), hi(3)
                 
              end select
              
           end do
        end do
        
     endif
     
2990 format(1x,i2)
     
2997 format(3x,"((",i5,              ") ("i5,") ( 0 ))  ", i5 )
2998 format(3x,"((",i5,",",i5,       ") ("i5,",",i5,") ( 0, 0 ))  ", i5)
2999 format(3x,"((",i5,",",i5,",",i5,") ("i5,",",i5,",",i5,") ( 0, 0, 0 ))  ", i5)
     
3001 format(6x,"((",i5,              ") ("i5,") ( 0 ))  ")
3002 format(6x,"((",i5,",",i5,       ") ("i5,",",i5,") ( 0, 0 ))  ")
3003 format(6x,"((",i5,",",i5,",",i5,") ("i5,",",i5,",",i5,") ( 0, 0, 0 ))  ")
     
     
     if (gridfile) then
        
        write (*,2990) pf%flevel
        
        do i = 1, pf%flevel
           
           lo = lwb(plotfile_get_pd_box(pf, i))
           hi = upb(plotfile_get_pd_box(pf, i))
           
           select case (dim)
           case (1)
              write (*,2997) lo(1), hi(1), nboxes(pf, i)
              
           case (2)
              write (*,2998) lo(1), lo(2), hi(1), hi(2), nboxes(pf, i)
              
           case (3)
              write (*,2999) lo(1), lo(2), lo(3), hi(1), hi(2), hi(3), nboxes(pf, i)
           end select
           
           do j = 1, nboxes(pf, i)
              
              lo = lwb(get_box(pf, i, j))
              hi = upb(get_box(pf, i, j))
              
              select case (dim)
              case (1)
                 write (*,3001) lo(1), hi(1)
                 
              case (2)
                 write (*,3002) lo(1), lo(2), hi(1), hi(2)
                 
              case (3)
                 write (*,3003) lo(1), lo(2), lo(3), hi(1), hi(2), hi(3)
                 
              end select
              
           end do
        end do
        
     endif
     
     if (castro) then
        
        if (pf%flevel .eq. 1) then
           print*, "ERROR: only use --castro option for multilevel"
           stop
        end if
        
        print*,pf%flevel-1
        
        do i = 2, pf%flevel
           
           print*, nboxes(pf, i)
           
           do j = 1, nboxes(pf, i)
              
              lo = lwb(get_box(pf, i, j))
              hi = upb(get_box(pf, i, j))
              
              select case (dim)
              case (1)
                 write (*,3001) lo(1)/2, (hi(1)-1)/2
                 
              case (2)
                 write (*,3002) lo(1)/2, lo(2)/2, (hi(1)-1)/2, (hi(2)-1)/2
                 
              case (3)
                 write (*,3003) lo(1)/2, lo(2)/2, lo(3)/2, (hi(1)-1)/2, (hi(2)-1)/2, (hi(3)-1)/2
                 
              end select
              
           end do
        end do
        
     endif
     
     
     print *, " "
     call destroy(pf)
     
  enddo

end program fboxinfo
