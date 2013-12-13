!
! this program calculates the buoyancy as a function of height, r:
!
!    b(r) = \int_{r0}^{r} N^2 dr'
!
! where N^2 is the Brunt-Vaisala frequency
!    
!    N^2 = - g * ( dln(rho)/dr - (dln(rho)/dr)_ad )
!        = (g/rho) * ((rho/(Gamma1*p))*dp/dr - drho/dr)
!
!
! ********** NOTE **********
!   The above quantities may need to be re-evaluated
! **************************
!
!
! another way of looking at this is to calculate the adiabatic excess,
! del nabla, using the standard adiabatic gradients, nabla =
! dln(T)/dln(p):
!
!    del nabla = nabla - nabla_ad
!
! the adiabatic nabla is related to the second gamma exponent (CG 9.88)
!  
!    (gamma2 - 1 ) / gamma2 = ( dlnT / dlnP )_ad = nabla_ad
!
! our EOS returns gamma1 only - it can be shown that gamma2 is related
! to gamma1 by:
!
!    nabla_ad = (gamma1 - chi_rho ) / (chi_T * gamma1)
!
! where
!
!    chi_rho = ( dlnP / dlnrho )_T         chi_T = ( dlnP / dlnT )_rho
!
! the input for this routine is an outputfile from the faverage
! routine.  this file should contain the following columns
!
!  #  y density_avg density_rms temp_avg temp_rms (species_avg species_rms)
!
! where temp_* can be either tfromp_* or tfromh_* and there is a
! (species_avg, species_rms) pair for each species in the network
!

program fbuoyancy

  use f2kcli
  use bl_error_module
  use bl_constants_module
  use bl_types
  use eos_module
  use eos_type_module
  use network

  implicit none

  character(len=256) :: input, output
  integer :: nlines
  real(kind=dp_t) :: gconst

  integer :: narg, farg
  character(len=256) :: fname, format_string

  integer :: i, j, ipos
  character(len=1024) :: line

  real(kind=dp_t) :: time, rms
  integer :: npoints
  real(kind=dp_t), allocatable :: height(:), dens(:), temp(:), comp(:,:)
  real(kind=dp_t), allocatable :: pres(:), gamma1(:), nabla_ad(:), n2(:)
  real(kind=dp_t) :: b, chi_t, chi_rho, adiabatic_excess

  real(kind=dp_t) :: dpdr, drhodr, dtdr, dtdP

  logical :: do_diag = .false.

  real(kind=dp_t), parameter :: small = 1.e-14

  type (eos_t) :: eos_state


  ! defaults
  input = ''
  output = 'out.dat'
  nlines = 1
  gconst = 2.450d14
  format_string = '("# ",g15.10)'

  ! get command line stuff
  narg = command_argument_count()
  farg = 1
  do while (farg <= narg)
     call get_command_argument(farg,value=fname)

     select case(fname)

     case('-i','--input')
        farg = farg + 1
        call get_command_argument(farg, value=input)

     case('-o', '--output')
        farg = farg + 1
        call get_command_argument(farg, value=output)

     case('-n', '--nlines')
        farg = farg + 1
        call get_command_argument(farg, value=fname)

        read(fname,*) nlines

     case('-g', '--gconst')
        farg = farg + 1
        call get_command_argument(farg, value=fname)

        read(fname,*) gconst

!      case('-f', '--format')
!         farg = farg + 1
!         call get_command_argument(farg, value=fname)

!         ichomp=index(trim(format_string), " ", back=.true.)
!         format_string = format_string(:ichomp) // trim(fname) // ')'

     case default
        exit

     end select
     farg = farg + 1
  enddo

  ! sanity checks
  if (input == '' .or. nlines == 1) then
     call print_usage()
     stop
  endif

  ! initialize
  call network_init()
  call eos_init()

  print *, ''
  print *, 'input file', trim(input)
  print *, 'output file', trim(output)
  print *, 'gconst', gconst
  print *, 'nlines', nlines
  print *, ''

  ! let's do this
  open(unit=99,file=trim(input))

  ! read in the time
  read(99,'(a1024)') line
  ipos = index(line,':')+1

  read(line(ipos:),*) time

  ! make sure we have all the species by parsing the second header line
  read(99,'(a1024)') line

  if (.not. all_species_included(line)) &
       call bl_error("not all species in network were included in inputfile")

  ! we have (nlines-2) number of points in the model
  npoints = nlines -2
  allocate(height(npoints), dens(npoints), temp(npoints), comp(npoints,nspec))
  allocate(pres(npoints), gamma1(npoints), nabla_ad(npoints), n2(npoints))

  ! read the rest of the file
  do i = 1, npoints
     read(99,*) height(i), dens(i), rms, temp(i), rms, &
          (comp(i,j), rms, j = 1, nspec)

     eos_state%rho = dens(i)
     eos_state%T = temp(i)
     eos_state%xn = comp(i,:)

     call eos(eos_input_rt, eos_state, do_diag)

     pres(i) = eos_state%p
     gamma1(i) = eos_state%gam1
     
     chi_t = temp(i) * eos_state%dpdt / pres(i)
     chi_rho = dens(i) * eos_state%dpdr / pres(i)

     nabla_ad(i) = (gamma1(i) - chi_rho) / (gamma1(i) * chi_t)
     
  enddo

  close(99)

  open(unit=88,file=trim(output))
  write(88,format_string) time
  write(88,'(a)') "# height n2 b adiabatic_excess"

  b = ZERO
  ! calculate b
  do i = 1, npoints
     
     ! build the gradients
     if (i == 1) then
        drhodr = (dens(i+1) - dens(i)) / (height(i+1) - height(i))
        dpdr   = (pres(i+1) - pres(i)) / (height(i+1) - height(i))
        dtdr   = (temp(i+1) - temp(i)) / (height(i+1) - height(i))

     else if (i == npoints) then
        drhodr = (dens(i) - dens(i-1)) / (height(i) - height(i-1))
        dpdr   = (pres(i) - pres(i-1)) / (height(i) - height(i-1))
        dtdr   = (temp(i) - temp(i-1)) / (height(i) - height(i-1))

     else
        ! center difference
        drhodr = (dens(i+1) - dens(i-1)) / (height(i+1) - height(i-1))
        dpdr   = (pres(i+1) - pres(i-1)) / (height(i+1) - height(i-1))
        dtdr   = (temp(i+1) - temp(i-1)) / (height(i+1) - height(i-1))

     endif

     if (abs(dpdr) .le. small ) then
        dtdP = ZERO
     else
        dtdP = dtdr / dpdr
     endif

     n2(i) = (gconst/dens(i)) * (dens(i)*dpdr/(gamma1(i)*pres(i)) - drhodr)

     ! build the integral for b using simple trapezoid rule
     if (i > 1) b = b + HALF*(n2(i)+n2(i-1))*(height(i)-height(i-1))

     adiabatic_excess = pres(i) * dtdP / temp(i) - nabla_ad(i)

     ! limit the adiabatic excess to lie between -10 and 10; in the
     ! fluff region, there are some zero gradients which can cause bad
     ! things for the finite differences
!     adiabatic_excess = max(-1.e1,min(1.e1,adiabatic_excess))

     write(88,*) height(i), n2(i), b, adiabatic_excess

  enddo

  close(88)

  contains 

    subroutine print_usage()
      
      implicit none

      print *, ''
      print *, 'Description: '
      print *, '  This program takes as input an outputfile from the '
      print *, '  faverage.f90 routine and calculates the Brunt-Vaisala'
      print *, '  frequency, buoyancy and adiabatic excess as a function of'
      print *, '  height.'
      print *, '  The input file must have the following format:'
      print *, '    #  y density_avg density_rms temp_avg temp_rms \\'
      print *, '    (species_avg species_rms) <whatever else>'
      print *, '  where temp can be either tfromp or tfromh and there is a'
      print *, '  (avg, rms) pair for each species in the network used to'
      print *, '  compile the fbuoyancy routine.'
      print *, ''
      print *, 'Usage: '
      print *, '  fbuoyancy <args>'
      print *, ''
      print *, 'Arguments: '
      print *, '  [-i|--input]   <filename>: '
      print *, '      specifies the input file; (required) '
      print *, '  [-o|--output]  <filename>:'
      print *, '      specifies the output file; default is "out.dat"'
      print *, '  [-n|--nlines]  <integer>:'
      print *, '      specifies total number of lines in inputfile; (required)'
      print *, '      one can use: -n `wc -l <inputfile> | awk "{print $1}"`'
      print *, '  [-g|--gconst]  <float>:'
      print *, '      specifies the gravitational constant; default is 2.45e14'
!       print *, '  [-f|--format] <format_string>:'
!       print *, '      This allows the user to set a FORTRAN format statement'
!       print *, '      for the output of the plotfile time at the beginning of.'
!       print *, '      the outputfile.'
!       print *, ' '
!       print *, '      <format_string> usually contains both a comment '
!       print *, '      character string and floating point format specifier, in'
!       print *, '      a standard FORTRAN format.  The default is set to '
!      print *, '      ("# ",g15.10).'
      print *, ''

    end subroutine print_usage


    ! parse the header string to see how many species are listed we
    ! will count the number of "X"'s in the string and divide this
    ! number by 2 (b/c of avg and rms quantities)
    ! if this matches nspec from network, then we say they are all included
    function all_species_included(string) result(r)

      character(len=*) :: string
      logical :: r

      integer :: spec_count

      ipos = index(string, 'X')
      spec_count = 0

      r = .false.

      do i = ipos, len(trim(string))
         if (string(i:i) == 'X') spec_count = spec_count+1
      enddo

      if (spec_count/2 == nspec) r = .true.

    end function all_species_included
         

end program fbuoyancy
