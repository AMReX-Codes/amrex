subroutine t_particle

  use particle_module

  integer i, id
  type(particle) p, q
  type(particle_vector) v

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

  call destroy(v)

end subroutine t_particle
