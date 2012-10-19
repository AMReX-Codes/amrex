
      subroutine PROBINIT (init,name,namlen,problo,probhi)

      use probdata_module
      implicit none

      integer init, namlen
      integer name(namlen)
      double precision problo(2), probhi(2)

      integer untin,i

      namelist /fortin/ probtype,specerr,specgrad,max_specerr_lev,max_specgrad_lev,&
                        diff_coeff

!
!     Build "probin" filename -- the name of file containing fortin namelist.
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
         
!     Read namelists
      untin = 9
      open(untin,file=probin(1:namlen),form='formatted',status='old')
      read(untin,fortin)
      close(unit=untin)

      end

      subroutine setregionnum(id,id_l1,id_l2,id_h1,id_h2,lo,hi,dx,plo)

      use probdata_module
      implicit none

      integer lo(2), hi(2)
      integer id_l1,id_l2,id_h1,id_h2
      double precision dx(2), plo(2)
      integer id(id_l1:id_h1,id_l2:id_h2)

      integer          :: i,j
      double precision :: x,y

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)

             x = plo(1) + (i-lo(1)+0.5d0)*dx(1) 
             y = plo(2) + (j-lo(2)+0.5d0)*dx(2) 

             if (x.ge.0.25d0 .and. x.le.0.75d0 .and. y.ge.0.25d0 .and. y.le.0.75d0) then
                id(i,j) = 1
             else
                id(i,j) = 2
             endif
         enddo
      enddo

      end subroutine setregionnum
