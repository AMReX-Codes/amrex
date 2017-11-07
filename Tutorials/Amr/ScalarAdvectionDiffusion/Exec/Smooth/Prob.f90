
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

  use amrex_fort_module, only : amrex_real, dim=>bl_spacedim
  implicit none
  integer, intent(in) :: level, lo(3), hi(3), phi_lo(3), phi_hi(3)
  double precision, intent(in) :: time
  double precision, intent(inout) :: phi(phi_lo(1):phi_hi(1), &
       &                                 phi_lo(2):phi_hi(2), &
       &                                 phi_lo(3):phi_hi(3))
  double precision, intent(in) :: dx(3), prob_lo(3)

  integer          :: i,j,k
  double precision :: xlo,ylo,xhi,yhi,zlo,zhi, pi, denom, integralval, freq
  
  freq = 1.0d0
  pi = 3.14159265358979323846264338327950288d0
  if(dim.eq.2) then
     denom = (dx(1)*dx(2)*pi*pi*freq*freq)
  else
     denom = (dx(1)*dx(2)*dx(3)*pi*pi*pi*freq*freq*freq)
  endif
  !$omp parallel do private(i,j,k,x,y,z,r2) collapse(2)
  do k=lo(3),hi(3)
     zlo = dx(3)*(k)
     zhi = dx(3)*(k+1)

     do j=lo(2),hi(2)
        ylo = dx(2)*(j)
        yhi = dx(2)*(j+1)

        do i=lo(1),hi(1)
           xlo = dx(1)*(i)
           xhi = dx(1)*(i+1)

           if(dim .eq. 2) then
              integralval = (cos(freq*pi*xhi) - cos(freq*pi*xlo))*(cos(freq*pi*yhi) - cos(freq*pi*ylo))
           else
              integralval = -(cos(freq*pi*xhi) - cos(freq*pi*xlo))*(cos(freq*pi*yhi) - cos(freq*pi*ylo))*(cos(freq*pi*zhi) - cos(freq*pi*zlo))
           end if

!!debug set to 1d           
!           integralval = (cos(freq*pi*xhi) - cos(freq*pi*xlo))
!           denom = dx(1)*pi*freq
!           if(j.eq.0) then
!!              if((i.eq.0).or.(i.eq.63)) then
!              if(i.eq.63) then
!!                 print*, "i xlo xhi sinxhi sinxlo integralval, denom", i, xlo, xhi, cos(freq*pi*xhi), cos(freq*pi*xlo), integralval, denom
!                 print*, "i xlo xhi sinxhi sinxlo integralval, denom", i, xlo, xhi, cos(freq*pi*1.0d0), cos(freq*pi*xlo), integralval, denom
!              endif
!           endif
!!!end debug
           phi(i,j,k) = integralval/denom

        end do
     end do
  end do
  !$omp end parallel do

!  j = 0
!  k = 0
!  print*, "*** init phi 0  1  62 63= ", j, phi(0,j,k),phi(1,j,k),phi(62,j,k),phi(63,j,k), "****"

end subroutine initdata
