subroutine t_particle

  use particle_module

  integer i, id
  type(particle) p, q
  type(particle_vector) v

  type(box)         :: bx,bx2
  type(ml_boxarray) :: mba
  type(ml_layout)   :: mla
  type(boxarray)    :: ba
  double precision  :: dx(1,MAX_SPACEDIM), dx2(2,MAX_SPACEDIM)
  double precision  :: problo(3)
  double precision  :: probhi(3)
  logical           :: pmask(3)

  id = 1

  pmask = .true.

  p%pos = (/ 1.d0, 2.d0, 3.d0 /)
  q%pos = (/ 1.d1, 2.d1, 3.d1 /)

  p%id = id; id = id + 1
  q%id = id; id = id + 1

  call print(p)
  call print(q)

  call build(v)

  if (parallel_IOProcessor()) then
     call print(v)
     call print(v, 'PV')
  end if

  do i = 1, 10
     p%id = id; id = id + 1
     q%id = id; id = id + 1
     call add(v,p)
     call add(v,q)
  end do

  if (parallel_IOProcessor()) then
     call print(v, 'PV')
     call print(v, 'PV', .true.)
     print*, 'size = ', size(v)
  end if
  !
  ! Now fill the vector to current capacity.
  !
  do while (size(v) < capacity(v))
     call add(v,p)
  end do

  p%id = 12345

  call remove(v, 1); call add(v,p)

  call remove(v, size(v)); call add(v,p)

  call remove(v, size(v)/2); call add(v,p)

  call remove(v, 2)
  call remove(v, 3)
  call remove(v, 4)
  call add(v,p)
  call add(v,p)
  call add(v,p)

  if (parallel_IOProcessor()) then
     call print(v, 'PV')
     print*, 'size = ', size(v)
  end if

  call clear(v)
  !
  ! Let's build a single level mla
  !
  call build(v)

  problo = -1.0d0
  probhi = +1.0d0

  bx = make_box( (/0,0,0/), (/31,31,31/) )

  call build(ba,bx)

  call boxarray_maxsize(ba,16)

  if (parallel_IOProcessor()) call print(ba)

  do i = 1, 3
     dx(1,i) = (probhi(i) - problo(i)) / extent(bx,i)
  end do

  call build(mba, ba, bx)

  call destroy(ba)

  if (parallel_IOProcessor()) print*, 'pmask: ', pmask

  call build(mla, mba, pmask)

  call destroy(mba)

  call init_random(v,1000,17971,mla,dx,problo,probhi)

  if (parallel_IOProcessor()) then
     print*, ''
     print*, 'size(v): ', size(v)
     print*, ''
     call print(v, 'after init_random')
     print*, ''
  end if

  call parallel_barrier()
  !
  ! Let's move the particles a bit.
  !
  do i = 1,100
     if (parallel_IOProcessor()) print*, i, 'Calling move_random(one-level mla) ...'
     call move_random(v,mla,dx,problo,probhi)
  end do

  if (parallel_IOProcessor()) then
     print*, ''
     print*, 'size(v): ', size(v)
     print*, ''
     call print(v, 'after move_random')
     print*, ''
     call flush(6)
  end if

  call parallel_barrier()
  call bl_error('Got Here')

  call destroy(mla)
  !
  ! Now let's try for a multi-level mla
  !
  call clear(v)

  call build(v)

  call destroy(ba)

  call build(mba,2,MAX_SPACEDIM)

  mba%rr(1,:) = 2

  call build(ba,bx)

  call boxarray_maxsize(ba,16)

  if (parallel_IOProcessor()) call print(ba, 'level 1')

  do i = 1, MAX_SPACEDIM
     dx2(1,i) = (probhi(i) - problo(i)) / extent(bx,i)
     dx2(2,i) = dx2(1,i) / 2.0d0
  end do

  call copy(mba%bas(1),ba)
  mba%pd(1) = bx

  if (parallel_IOProcessor()) call print(mba%pd(1), 'pd(1)')

  call destroy(ba)

  bx2 = make_box( (/16,16,16/), (/47,47,47/) )

  call build(ba,bx2)

  call boxarray_maxsize(ba,16)

  if (parallel_IOProcessor()) call print(ba, 'level 2')

  call copy(mba%bas(2),ba)
  mba%pd(2) = refine(mba%pd(1),2)

  if (parallel_IOProcessor()) call print(mba%pd(2), 'pd(2)')

  call destroy(ba)

  call build(mla, mba, pmask)

  call destroy(mba)

  call init_random(v,1000,171717171,mla,dx2,problo,probhi)

  if (parallel_IOProcessor()) then
     print*, ''
     print*, 'size(v): ', size(v)
     print*, ''
     call print(v, 'after init_random using 2-level mla')
     print*, ''
  end if
  !
  ! Let's move the particles a bit.
  !
  do i = 1,100
     if (parallel_IOProcessor()) print*, i, 'Calling move_random(two-level mla) ...'
     call move_random(v,mla,dx2,problo,probhi)
  end do

  call destroy(mla)

  call destroy(v)

end subroutine t_particle
