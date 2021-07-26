
subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

  use amrex_fort_module
  use probdata_module
  use, intrinsic :: iso_fortran_env, only : unterr => error_unit ! Fortran 2003

  implicit none

  integer(c_int), intent(in) :: init, namlen
  integer(c_int), intent(in) :: name(namlen)
  real(amrex_real), intent(in) :: problo(*), probhi(*)

  integer :: untin, i

  namelist /fortin/ adv_vel

  !
  ! Build "probin" filename -- the name of file containing fortin namelist.
  !
  integer, parameter :: maxlen = 256
  character(len=maxlen) :: probin

  if (namlen > maxlen) then
     write(unterr,*) 'probin file name too long'
     stop
  end if

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do

  ! set the namelist default
  adv_vel(:) = 1.d0

  ! Read namelists
  untin = 9
  open(untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin,fortin)
  close(unit=untin)

end subroutine amrex_probinit


subroutine initdata(level, time, lo, hi, &
     phi, phi_lo, phi_hi, &
     dx, prob_lo) bind(C, name="initdata")

  use amrex_fort_module
  implicit none
  integer(c_int), intent(in) :: level, lo(3), hi(3), phi_lo(3), phi_hi(3)
  real(amrex_real), intent(in) :: time
  real(amrex_real), intent(inout) :: phi(phi_lo(1):phi_hi(1), &
       &                                 phi_lo(2):phi_hi(2), &
       &                                 phi_lo(3):phi_hi(3))
  real(amrex_real), intent(in) :: dx(3), prob_lo(3)

  integer          :: dm
  integer          :: i,j,k
  real(amrex_real) :: x,y,z,r2

  if (phi_lo(3) == 0 .and. phi_hi(3) == 0) then
     dm = 2
  else
     dm = 3
  end if

  !$omp parallel do private(i,j,k,x,y,z,r2) collapse(2)
  do k = lo(3), hi(3)
     do j = lo(2), hi(2)

        z = prob_lo(3) + (real(k,kind=amrex_real)+0.5d0) * dx(3)
        y = prob_lo(2) + (real(j,kind=amrex_real)+0.5d0) * dx(2)

        ! The compiler should automatically convert this innermost loop to SIMD
        do i = lo(1), hi(1)

           x = prob_lo(1) + (real(i,kindo=amrex_real)+0.5d0) * dx(1)

           if ( dm == 2 ) then
              r2 = ((x-0.0d0)**2 + (y-0.0d0)**2) / 0.01d0
              phi(i,j,k) = 1.d0 + exp(-r2)
           else
              r2 = ((x-0.0d0)**2 + (y-0.0d0)**2 + (z-0.0d0)**2) / 0.01d0
              phi(i,j,k) = 1.d0 + exp(-r2)
           end if ! dm

        end do ! i
     end do ! j
  end do ! k
  !$omp end parallel do

end subroutine initdata
