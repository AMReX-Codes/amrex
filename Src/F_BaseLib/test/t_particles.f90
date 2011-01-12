subroutine t_particle

  use particle_module

  integer i, id
  type(particle) p, q
  type(particle_vector) v

  type(box)         :: bx
  type(ml_boxarray) :: mba
  type(ml_layout)   :: mlla
  type(boxarray)    :: ba
  double precision  :: dx(1,3)
  double precision  :: problo(3)
  double precision  :: probhi(3)

  id = 1

  p%pos = (/ 1.d0, 2.d0, 3.d0 /)
  q%pos = (/ 1.d1, 2.d1, 3.d1 /)

  p%id = id; id = id + 1
  q%id = id; id = id + 1

  call print(p)
  call print(q)

  call build(v)

  call print(v)

  call print(v, 'PV')

  do i = 1, 10
     p%id = id; id = id + 1
     q%id = id; id = id + 1
     call add(v,p)
     call add(v,q)
  end do

!  call print(v, 'PV')

!  call print(v, 'PV', .true.)

  print*, 'size = ', size(v)

  !
  ! Now fill the vector to current capacity.
  !
  do while (size(v) < capacity(v))
     call add(v,p)
  end do

  p%id = 12345

  call remove(v, 1); call add(v,p); print*, 'size = ', size(v)

  call remove(v, size(v)); call add(v,p); print*, 'size = ', size(v)

  call remove(v, size(v)/2); call add(v,p); print*, 'size = ', size(v)

  call remove(v, 2)
  call remove(v, 3)
  call remove(v, 4)
  call add(v,p); print*, 'size = ', size(v)
  call add(v,p); print*, 'size = ', size(v)
  call add(v,p); print*, 'size = ', size(v)

  call print(v, 'PV')

  call clear(v)

  call build(v)

  problo = -1.0d0
  probhi = +1.0d0

  bx = make_box( (/0,0,0/), (/31,31,31/) )

  call build(ba,bx)

  call boxarray_maxsize(ba,16)

  call print(ba)

  do i = 1, 3
     dx(1,i) = (probhi(i) - problo(i)) / extent(bx,i)
  end do

  call build(mba, ba, bx)

  call destroy(ba)

  call build(mlla, mba)

  call destroy(mba)

  call init_random(v,100,17971,mlla,dx,problo,probhi)

  call destroy(mlla)

  print*, '**************************************************'

  print*, 'size(v): ', size(v)

  call print(v, 'after init_random')

  print*, '**************************************************'

  call destroy(v)

end subroutine t_particle
