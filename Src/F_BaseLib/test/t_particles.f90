subroutine t_particle

  use particle_module

  type(particle) p, q

  p%pos = (/ 1.d0, 2.d0, 3.d0 /)
  q%pos = (/ 1.d1, 2.d1, 3.d1 /)

  call print(p)
  call print(q)

end subroutine t_particle
