      subroutine fillfab(d,nx,ny)
      implicit none

      double precision d(nx,ny,2)
      integer nx, ny
      character :: filename*100
      real delx
      integer i, j

      filename = './field_full_native_double.dat'

      open(unit=2345, file=trim(filename), status='old',
     &     form='unformatted')

      read(2345) i, j
      read(2345) delx

      if (nx .ne. i) stop 'nx != i'
      if (ny .ne. j) stop 'ny != j'

      read(2345) d(1:nx,1:ny,1:2)
      close(2345)

      print*, 'i: ', i, 'j: ', j
      print*, 'delx: ', delx

      end subroutine fillfab
