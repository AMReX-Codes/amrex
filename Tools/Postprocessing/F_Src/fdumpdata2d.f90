!! fdumpdata2d
!!
!! usage: fdumpdata2d --pf (pfs) [--var (vars)] [--lev (levs)] [--separate]
!!        where (pfs) = list of plotfiles,
!!              (vars) = list of variables,
!!              (levs) = list of levels,
!!              --var, --lev, --seprate are optional.
!!
!! This program generates plain-text-formatted data files from plotfiles.
!!
!! A resulting dump data file has columns (or blocks of columns separated 
!! by two empty lines). Values in columns correspond to:
!!   (x-coord) (y-coord) (var1_val) (var2_val) ... (varN_val)
!! If there are several blocks, each block contains data for a level.
!!
!! If plotfile (actually directory) is abc/def/ghi/,
!! resulting dump data file(s) is abc/def/dat.ghi (without --separate option)
!! or abc/def/dat.ghi.lev* (with --separate option)
!!
!! Using --pf option, you can specify plotfiles:
!!   --pf plotfile1 plotfile2 ... plotfileM
!!
!! Using --var option, you can specify variables:
!!   --var var1 var2 ... varN 
!! If not given, all variables are used.
!!
!! Using --lev option, you can specify levels:
!!   --lev lev1 lev2 ... levK
!! If not given, all levels are used. 
!!
!! If --separate option is used, a seperate file for each level is generated.
!! That is, a total of M*K files are generated:
!!   plotfileI.dat.levJ (I=1..M, J=1..K) 
!! Otherwise, for each plotfile, a single data file is generated.
!! That is, a total of M files are generated:
!!   plotfileI.dat (I=1..M) 
!! In this case, there are K blocks of data seprated by two empty lines.
!!

!! example commands for gnuplot
!!
!! 1. plot cell points of kth level
!! gnuplot> plot "data.plotfile.levk" u 1:2 w p
!! gnuplot> plot "data.plotfile" index k-1 u 1:2 w p
!!
!! 2. plot values of nth variable at kth level
!! gnuplot> splot "dat.plotfile.levk" using 1:2:3 with impulse
!! gnuplot> splot "dat.plotfile" index k-1 using 1:2:3 with impulse


program fdumpdata2d 

  use bl_space, only: MAX_SPACEDIM
  use bl_constants_module, only: HALF
  use bl_IO_module
  use plotfile_module

  implicit none

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! variables for handling command-line arguments !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! maximum length of an argument in command line
  integer, parameter :: MAX_CMD_ARG_LEN = 256
  ! maximum length of a long string (see below, FMT, str1, str2) 
  integer, parameter :: MAX_LONG_STR_LEN = 1024
  ! maximum length of a filename
  integer, parameter :: MAX_FILENAME_LEN = 256
  ! separator for file path / or \
  character, parameter :: FILE_PATH_SEPARATOR = '/'

  ! for --pf option (pf=plotfile)
  logical :: is_defined_pf = .false.   ! whether it is defined in cmd line
  integer :: cmd_arg_npf               ! number of pfs
  character(len=MAX_CMD_ARG_LEN), allocatable :: cmd_arg_pf_list(:)
                                       ! list of pfs in cmd line

  ! for --var option (var=variable)
  logical :: is_defined_var = .false.  ! whether it is defined in cmd line
  integer :: cmd_arg_nvar              ! number of vars in cmd line
  character(len=MAX_CMD_ARG_LEN), allocatable :: cmd_arg_var_list(:)
                                       ! list of vars in cmd line 

  ! for --lev option (lev=level)
  logical :: is_defined_lev = .false.  ! whether it is defined in cmd line
  integer :: cmd_arg_nlev              ! number of levs in cmd line
  integer, allocatable :: cmd_arg_lev_list(:)
                                       ! list of levs in cmd line 

  ! for --separate option
  ! if used, a data file is generated for each lev data
  logical :: is_defined_separate = .false.  
                                       ! whether it is defined in cmd line

  ! nvar and var_list will be used in main part.
  ! if --var option is given, nvar and var_list are the same as
  ! cmd_arg_nvar and cmd_arg_var_list.
  ! otherwise, for each plot file, nvar and var_list will be obtained.
  integer :: nvar
  character(len=MAX_CMD_ARG_LEN), allocatable :: var_list(:)

  ! nlev and lev_list will be used in main part.
  ! if --lev option is given, nlev and lev_list are the same as
  ! cmd_arg_nlev and cmd_arg_lev_list.
  ! otherwise, for each plot file, nlev and lev_list will be obtained.
  integer :: nlev             
  integer, allocatable :: lev_list(:)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! variables for main part !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

  ! pf is a plotfile object that will contain both the meta-data
  ! describing the layout of the boxes that make up the AMR levels and
  ! the actual data itself (although, not necessarily all of the data
  ! at any one time, to save memory).
  type(plotfile) pf

  ! the pointer p will point to the data for a single box in the
  ! computational domain. Regardless of the dimensionality of the
  ! simulation, p will always have 4 indices (3 space + 1 component).
  real(kind=dp_t), pointer :: p(:,:,:,:)

  integer :: ipf   ! index for cmd_arg_pf_list 
  integer :: ntvar ! number of total variables in a plot file  cf. nvar
  integer :: ivar  ! index for variables 
  integer :: ntlev ! number of total levles in a plot file  cf. nlev
  integer :: ilev  ! index for levels 
  integer :: lev   ! current level (i.e. lev = lev_list(ilev))
  integer :: nbox  ! number of boxes at a level
  integer :: ibox  ! index for boxes 
  integer :: ix,iy ! indices for cells

  real(kind=dp_t) :: dx(MAX_SPACEDIM)  ! grid size of a level
  integer :: lo(MAX_SPACEDIM)          ! lower bound for a box
  integer :: hi(MAX_SPACEDIM)          ! upper bound for a box

  real(kind=dp_t) :: xx,yy             ! center position of a cell

  ! list of variable indices for variables in var_list
  integer, allocatable :: var_index_list(:)  

  integer :: unit1 ! for file input (i.e. plotfiles)
  integer :: unit2 ! for file output (i.e. dump files)
  character(len=MAX_FILENAME_LEN) :: dump_file ! filename for output

  !! temporary long strings 
  character(len=MAX_LONG_STR_LEN) :: str1, str2  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! handling command-line arguments !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! read arguments in command line
  call read_cmd_arg

  ! check arguments
  call check_cmd_arg

  ! if --var option has been given, copy nvar and var_list
  if (is_defined_var) then
    nvar = cmd_arg_nvar
    allocate(var_list(nvar))
    var_list = cmd_arg_var_list
  end if

  ! if --lev option has been given, copy nlev and lev_list
  if (is_defined_lev) then
    nlev = cmd_arg_nlev
    allocate(lev_list(nlev))
    lev_list = cmd_arg_lev_list
  end if

  !!!!!!!!!!!!!!!
  !! main part !!
  !!!!!!!!!!!!!!!

  unit1 = unit_new()  !! will be used for file input
  unit2 = unit_new()  !! will be used for file output

  do ipf = 1, cmd_arg_npf

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    print '(A,A,A)', "** reading ", trim(cmd_arg_pf_list(ipf)), "..."

    ! build the plotfile object that contains the data describing
    ! the plotfile
    call build(pf, cmd_arg_pf_list(ipf), unit1)

    ! make sure we are 2-d
    if (pf%dim /= 2) then
      print *, "ERROR: not a 2-d dataset"
      stop
    end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! get ntvar and ntlev 
    ntvar = pf%nvars
    ntlev = pf%flevel

    ! if --var option has not been given, assign values for nvar and var_list
    if (.not. is_defined_var) then
      nvar = ntvar 
      allocate(var_list(nvar))
      do ivar = 1, nvar
        var_list(ivar) = plotfile_var_name(pf,ivar) 
      end do
    end if

    ! if --lev option has not been given, assign values for nlev and lev_list
    if (.not. is_defined_lev) then
      nlev = ntlev
      allocate(lev_list(nlev))
      do ilev = 1, nlev
        lev_list(ilev) = ilev
      end do
    end if

    ! construct var_index_list
    allocate(var_index_list(nvar))
    do ivar = 1, nvar
      var_index_list(ivar) = plotfile_var_index(pf,var_list(ivar))
      if (var_index_list(ivar)<0) then
        print *, "ERROR: invalid variable ", trim(var_list(ivar))
        stop
      end if
    end do

    ! check lev_list
    do ilev = 1, nlev
      lev = lev_list(ilev)
      if (lev < 1 .or. lev > ntlev) then
        print *, "ERROR: level is out of range"
        stop
      end if
    end do

    ! print whole list of vars defined in plot file and var_list
    print '(A,A)', " - variables defined in plotfile: ", &
          trim(conv_str_array_to_str_with_comma( &
          (/(plotfile_var_name(pf,ivar), ivar = 1, ntvar)/) ))

    print '(A,A)', " - requested variables: ", &
          trim(conv_str_array_to_str_with_comma(var_list))

    ! print ntlev and lev_list
    print '(A,I0)', " - number of levels in plotfile = ", ntlev 
    print '(A,A)', " - requested levels: ", &
          trim(conv_int_array_to_str_with_comma(lev_list))

    ! print refinement ratios (only for multi-level case)
    if (ntlev>1) then
      print '(A,"x=(",A,"), y=(",A,")")', " - refinement ratios: ", &
            trim(conv_int_array_to_str_with_comma(pf%refrat(:,1))), &
            trim(conv_int_array_to_str_with_comma(pf%refrat(:,2)))
    end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! now explore each level !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! if --separate option has not been invoked, open out_file at this point 
    if (.not. is_defined_separate) then
      ! file open for output
      ! note that plotfile is actually a directory.
      ! if plotfile is abc/def/ghi/,
      ! then output file name is abc/def/dat.ghi
      call separate_last_item(trim(cmd_arg_pf_list(ipf)),str1,str2)
      write (dump_file,'(A,A,A)') trim(str1),"dat.",trim(str2)
      open(unit=unit2,file=dump_file,action='write')
    end if

    do ilev = 1, nlev

      ! get info of level 
      lev = lev_list(ilev)
      dx = plotfile_get_dx(pf,lev)
      nbox = nboxes(pf,lev)

      ! if --separate option has been invoked, open out_file for each level
      if (is_defined_separate) then
        ! file open for output
        ! note that plotfile is actually a directory.
        ! if plotfile is abc/def/ghi/,
        ! then output file name is abc/def/dat.ghi.levN
        call separate_last_item(trim(cmd_arg_pf_list(ipf)),str1,str2)
        write (dump_file,'(A,A,A,A,I0)') &
              trim(str1),"dat.",trim(str2),".lev",ilev
        open(unit=unit2,file=dump_file,action='write')
      end if

      ! print info 
      print '(A,I0,": dx=",E10.3,", dy=",E10.3,", nbox=",I0,", dump_file=",A)', &
            " * level ", lev, dx(1), dx(2), nbox, trim(dump_file)

      ! for no --separate option, indicate the level with gnuplot comment '#'
      if (.not. is_defined_separate) write (unit2,'(A,I0)') "# level ", ilev

      ! go into each box
      do ibox = 1, nbox

        ! read in the data 1 patch at a time
        call fab_bind(pf,ilev,ibox)

        ! get the integer bounds of the current box, in terms of
        ! this level's index space
        lo = lwb(get_box(pf,ilev,ibox))
        hi = upb(get_box(pf,ilev,ibox))

        ! get a pointer to the current patch
        p => dataptr(pf,ilev,ibox)

        ! print info 
        print '(A,I0,": lo=(",I0,",",I0,"), hi=(",I0,",",I0,")")', &
              "  + box ",ibox,lo(1),lo(2),hi(1),hi(2)

        ! loop over cells
        do iy = lo(2), hi(2)
          yy = pf%plo(2) + (iy+HALF)*dx(2)

          do ix = lo(1), hi(1)
            xx = pf%plo(1) + (ix+HALF)*dx(1)

            write (unit2,*) xx, yy, (p(ix,iy,1,var_index_list(ivar)),ivar=1,nvar)
          end do

        end do

      end do

      ! for no --separate option, put two empty lines to make a block 
      ! for each level data.
      ! this is for gnuplot index feature.
      ! e.g. splot "dat.plot000" index 0 u 1:2:3 w i  (for first level)
      !      splot "dat.plot000" index 1 u 1:2:3 w i  (for second level)
      if (.not. is_defined_separate) then
        write (unit2,*) " "
        write (unit2,*) " "
      end if

      if (is_defined_separate) close(unit=unit2)  ! close at this point

    end do

    if (.not. is_defined_separate) close(unit=unit2)  ! close at this point

    ! if --var option has not been given, deallocate var_list for next use
    if (.not. is_defined_var) deallocate(var_list)
   
    ! if --lev option has not been given, deallocate lev_list for next use
    if (.not. is_defined_lev) deallocate(lev_list)

    ! whether --var option has been given or not, deallocate var_index
    deallocate(var_index_list)
   
  end do 

  contains

  !! print usage of this program
  !! usually when error occurs
  subroutine print_usage

    print *, "************************************************************************"
    print *, "usage: fdumpdata2d --pf (pfs) [--var (vars)] [--lev (levs)] [--separate]"
    print *, "where (pfs) = list of plotfiles,"
    print *, "      (vars) = list of variables,"
    print *, "      (levs) = list of levels,"
    print *, "      --var, --lev, --seprate are optional."
    print *, "************************************************************************"

  end subroutine print_usage

  !! read command line and get items for each -- options
  subroutine read_cmd_arg

    integer :: i, pos
    character(len=MAX_CMD_ARG_LEN) :: arg
    integer :: narg

    ! get the number of arguments in command line 

    narg = command_argument_count()

    ! analyze cmd args

    pos = 1  ! position to be read at

    do i = 1, narg

      if (i<pos) cycle  ! skip if current position is not the one for reading 
                        ! actually, this position was already read by
                        ! cmd_arg_get_items_str (or cmd_arg_get_items_int)

      ! otherwise, read

      call get_command_argument(i, value = arg)  ! i and pos are the same

      ! if it contains --, then it is an option. 
      if (does_start_with(arg,"--")) then  
                                           
        ! --pf option
        if (arg=="--pf") then

          if (is_defined_pf) then  ! already defined
            print *, "ERROR: --pf option is repeated"
            call print_usage
            stop
          else
            is_defined_pf = .true.
            call cmd_arg_get_items_str(i,cmd_arg_npf,cmd_arg_pf_list) 
            pos = i + cmd_arg_npf + 1
          end if

        ! --var option
        else if (arg=="--var") then

          if (is_defined_var) then  ! already defined
            print *, "ERROR: --var option is repeated"
            call print_usage
            stop
          else
            is_defined_var = .true.
            call cmd_arg_get_items_str(i,cmd_arg_nvar,cmd_arg_var_list) 
            pos = i + cmd_arg_nvar + 1
          end if

        ! --lev option
        else if (arg=="--lev") then

          if (is_defined_lev) then  ! already defined
            print *, "ERROR: --lev option is repeated"
            call print_usage
            stop
          else
            is_defined_lev = .true.
            call cmd_arg_get_items_int(i,cmd_arg_nlev,cmd_arg_lev_list) 
            pos = i + cmd_arg_nlev + 1
          end if

        ! --separate option
        else if (arg=="--separate") then
     
          if (is_defined_separate) then  ! already defined
            print *, "ERROR: --separate option is repeated"
            call print_usage
            stop
          else
            is_defined_separate = .true.
            pos = i+1
          end if 

        ! otherwise, it must be a strange option
        else  

          print *, "ERROR: not supported option ", trim(arg)
          call print_usage
          stop

        end if
    
      ! unexpected item which does not follow any -- option 
      else  

        print *, "ERROR: unexpected item ", trim(arg)
        call print_usage
        stop

      endif
    end do

  end subroutine read_cmd_arg

  !! after getting options by using read_cmd_arg,
  !! check options are correct
  subroutine check_cmd_arg

    integer :: i

    ! check --pf option

    if (.not. is_defined_pf) then
      print *, "ERROR: no --pf option"
      call print_usage
      stop
    else if (cmd_arg_npf==0) then
      print *, "ERROR: no plotfiles are provided after --pf option"
      call print_usage
      stop
    end if

    ! check --var option

    if (is_defined_var .and. cmd_arg_nvar==0) then
      print *, "ERROR: no variables are provided after --var option"
      call print_usage
      stop
    endif

    ! check --lev option

    if (is_defined_lev .and. cmd_arg_nlev==0) then
      print *, "ERROR: no levels are provided after --lev option"
      call print_usage
      stop
    endif

  end subroutine check_cmd_arg

  function does_start_with(str,pattern) result (r)

    character(len=*), intent(in) :: str
    character(len=*), intent(in) :: pattern
    logical :: r

    if (len(str)<len(pattern)) then       ! if str is shorter than pattern
      r = .false.  
    else    
      r = (str(1:len(pattern))==pattern)  ! otherwise, compare 
    end if
 
  end function does_start_with

  !! get the list of items for the -- option located at nth argument position
  !! for example, exec --pf a b c --lev 1 3 --separate
  !! in order to get items for --pf, use pos=1
  !! then will get n=3 and list=(a,b,c)
  subroutine cmd_arg_get_items_str(pos,n,list)

    integer, intent(in)  :: pos  ! position of --option
    integer, intent(out) :: n    ! number of items btw pos+1 & next --option
    character(len=MAX_CMD_ARG_LEN), intent(out), allocatable :: list(:)
                                 ! list of items
    integer :: i, j
    integer :: narg
    character(len=MAX_CMD_ARG_LEN) :: arg

    ! find the position of next --option

    narg = command_argument_count()

    do i = pos+1, narg   ! otherwise, find next --option
      call get_command_argument(i, value = arg)
      if (does_start_with(arg,"--")) exit
    end do

    ! get the number of items

    n = i-pos-1

    ! allocate list and get the list

    if (n>0) allocate(list(n))

    do i = 1, n
      call get_command_argument(pos+i, value = arg)
      list(i) = arg
    end do 

  end subroutine cmd_arg_get_items_str

  !! get the list of items for the -- option located at nth argument position
  !! for example, exec --pf a b c --lev 1 3 --separate
  !! in order to get items for --lev, use pos=5
  !! then will get n=2 and list=(1,3)
  subroutine cmd_arg_get_items_int(pos,n,list)

    integer, intent(in)  :: pos  ! position of --option
    integer, intent(out) :: n    ! number of items btw pos+1 & next --option
    integer, intent(out), allocatable :: list(:)
                                 ! list of items
    integer :: i, j
    integer :: narg
    character(len=MAX_CMD_ARG_LEN) :: arg

    ! find the position of next --option

    narg = command_argument_count()

    do i = pos+1, narg   ! otherwise, find next --option
      call get_command_argument(i, value = arg)
      if (does_start_with(arg,"--")) exit
    end do

    ! get the number of items

    n = i-pos-1

    ! allocate list and get the list

    if (n>0) allocate(list(n))

    do i = 1, n
      call get_command_argument(pos+i, value = arg)
      read (arg,'(I20)') list(i)
    end do 

  end subroutine cmd_arg_get_items_int

  !! For a given integer array (e.g., array = (/1,2,3/)),
  !! this function generates a string for this array,
  !! where components are separated by commas (e.g., str = "1,2,3")
  function conv_int_array_to_str_with_comma(array) result (str)

    integer :: array(:)
    character(len=MAX_LONG_STR_LEN) :: str
    integer :: i
    integer :: n
    character(len=MAX_LONG_STR_lEN) :: FMT

    n = size(array)

    ! if there is one integer
    if (n==1) then
      write (str,'(I0)') array(1)
      return
    end if

    ! otherwise
    write (fmt,'("(",I0,"(I0,A),I0)")') n-1
    write (str,fmt) (array(i), ",", i = 1, n-1), array(n)
    return

  end function conv_int_array_to_str_with_comma

  !! For a given array of strings (e.g., array = (/"abc","de ","f  "/)),
  !! this function generates a string for this array,
  !! where components are separated by commas (e.g., str = "abc,de,f")
  function conv_str_array_to_str_with_comma(array) result (str)

    character(len=*) :: array(:)
    character(len=MAX_LONG_STR_LEN) :: str
    integer :: i
    integer :: n
    character(len=MAX_LONG_STR_lEN) :: FMT

    n = size(array)

    ! if there is one string 
    if (n==1) then
      write (str,'(A)') trim(array(1))
      return
    end if

    ! otherwise
    write (fmt,'("(",I0,"(A,A),A)")') n-1
    write (str,fmt) (trim(array(i)), "," , i = 1, n-1), trim(array(n))
    return

  end function conv_str_array_to_str_with_comma 

  !! separate the last item from path
  !! for example, for path1=abc/def/ghi, get path2=abc/def/ and last=ghi
  !! additionally, if the last character is the file path separate,
  !! it is deleted.
  !! for example, for path1=abc/def/ghi/, then path2=abc/def/ and last=ghi
  subroutine separate_last_item(path1,path2,last)

    character(len=*), intent(in)  :: path1  ! path1=path2/last  
    character(len=*), intent(out) :: path2
    character(len=*), intent(out) :: last
    integer :: n, i
    logical :: found = .false.

    n = len(path1)
    if (path1(n:n)==FILE_PATH_SEPARATOR) n = n-1  ! if any, it will be ignored 

    do i = n, 2, -1
      if (path1(i:i)==FILE_PATH_SEPARATOR) then
        found = .true.
        exit
      end if
    end do

    if (found) then
      path2 = path1(1:i)  ! separator included
      last = path1(i+1:n) 
    else
      path2 = ''
      last = path1(1:n)
    end if

  end subroutine

end program fdumpdata2d 
