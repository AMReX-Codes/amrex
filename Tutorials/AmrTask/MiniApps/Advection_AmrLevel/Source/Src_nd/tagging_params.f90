module tagging_params_module

  double precision, save :: phierr(0:15), phigrad(0:15)

  integer, save :: max_phierr_lev, max_phigrad_lev

contains
  
  subroutine get_tagging_params(name, namlen) bind(C, name="get_tagging_params")

    ! Initialize the tagging parameters

    integer, intent(in) :: namlen
    integer, intent(in) :: name(namlen)

    integer :: un, i, status

    integer, parameter :: maxlen = 256
    character (len=maxlen) :: probin

    namelist /tagging/ phierr, phigrad, max_phierr_lev, max_phigrad_lev

    ! Set namelist defaults
    phierr(:) = 1.d20
    phigrad(:) = 1.d20
    max_phierr_lev = -1
    max_phigrad_lev = -1

    ! create the filename
    if (namlen > maxlen) then
       print *, 'probin file name too long'
       stop
    endif

    do i = 1, namlen
       probin(i:i) = char(name(i))
    end do

    ! read in the namelist
    un = 9
    open (unit=un, file=probin(1:namlen), form='formatted', status='old')
    read (unit=un, nml=tagging, iostat=status)

    if (status < 0) then
       ! the namelist does not exist, so we just go with the defaults
       continue

    else if (status > 0) then
       ! some problem in the namelist
       print *, 'ERROR: problem in the tagging namelist'
       stop
    endif

    close (unit=un)

  end subroutine get_tagging_params

end module tagging_params_module
