
subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

  use probdata_module

  implicit none

  integer, intent(in) :: init, namlen
  integer, intent(in) :: name(namlen)
  double precision, intent(in) :: problo(*), probhi(*)

  integer untin,i

  namelist /fortin/ adv_vel
  
  !
  ! Build "probin" filename -- the name of file containing fortin namelist.
  !     
  integer maxlen
  parameter (maxlen=256)
  character probin*(maxlen)
  
  if (namlen .gt. maxlen) then
     write(6,*) 'probin file name too long'
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

  implicit none
  integer, intent(in) :: level, lo(3), hi(3), phi_lo(3), phi_hi(3)
  double precision, intent(in) :: time
  double precision, intent(inout) :: phi(phi_lo(1):phi_hi(1), &
       &                                 phi_lo(2):phi_hi(2), &
       &                                 phi_lo(3):phi_hi(3))
  double precision, intent(in) :: dx(3), prob_lo(3)

  integer          :: dm
  integer          :: i,j,k
  double precision :: x,y,z,r2
  
  if (phi_lo(3) .eq. 0 .and. phi_hi(3) .eq. 0) then
     dm = 2
  else
     dm = 3
  end if

  !$omp parallel do private(i,j,k,x,y,z,r2) collapse(2)
  do k=lo(3),hi(3)
     z = prob_lo(3) + (dble(k)+0.5d0) * dx(3)
     do j=lo(2),hi(2)
        y = prob_lo(2) + (dble(j)+0.5d0) * dx(2)
        do i=lo(1),hi(1)
           x = prob_lo(1) + (dble(i)+0.5d0) * dx(1)
           
           if ( dm.eq. 2) then
              r2 = ((x-0.0d0)**2 + (y-0.0d0)**2) / 0.01d0
              phi(i,j,k) = 1.d0 + exp(-r2)
           else
              r2 = ((x-0.0d0)**2 + (y-0.0d0)**2 + (z-0.0d0)**2) / 0.01d0
              phi(i,j,k) = 1.d0 + exp(-r2)
           end if
        end do
     end do
  end do
  !$omp end parallel do

end subroutine initdata
