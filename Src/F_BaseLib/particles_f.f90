
!
! The data structures in this code aren't particularly sophisticated.
! It can handle some tens or hundreds of thousands of particles OK.
! It's when you get into the millions of particles range when there
! could start to be issues.  I'm thinking of how to beef up the data
! structures.
!

module particle_module

  use bl_space
  use bl_constants_module
  use vector_i_module
  use multifab_module
  use ml_layout_module

  implicit none

  type particle

     integer :: id   = -1
     integer :: cpu  = -1
     integer :: lev  = -1
     integer :: grd  = -1

     integer :: cell(MAX_SPACEDIM)

     double precision :: pos(MAX_SPACEDIM)
     double precision :: origpos(MAX_SPACEDIM)

  end type particle

  interface index
     module procedure particle_index
  end interface index

  interface where
     module procedure particle_where
  end interface where

  interface periodic_shift
     module procedure particle_periodic_shift
  end interface periodic_shift

  interface print
     module procedure particle_print
  end interface print

  type particle_container
     private
     integer                 :: size = 0
     type(particle), pointer :: d(:) => NULL()
     type(vector_i)          :: invalid
     logical                 :: built = .false.
  end type particle_container

  interface build
     module procedure particle_container_build
  end interface build

  interface destroy
     module procedure particle_container_destroy
  end interface destroy
  !
  ! This is just the number of "valid" particles.
  !
  ! Valid particles are those for which "id" > 0.
  !
  ! Invalid are those for which "id" <= 0.
  !
  interface size
     module procedure particle_container_size
  end interface size

  interface empty
     module procedure particle_container_empty
  end interface empty

  interface ok
     module procedure particle_container_okay
  end interface ok

  interface swap
     module procedure particle_container_swap
  end interface swap
  !
  ! This symbolically removes particles from the vector
  ! by negating the "id" of the particle.  This will be
  ! used by add() to cut down on memory allocation and copying.
  !
  interface remove
     module procedure particle_container_remove
  end interface remove

  !
  ! This "usually" does a push_back on the space holding
  ! the particles.  If however that space is at capacity,
  ! it will try to add the particle by overwriting one that was
  ! previously removed, if possible, before allocating more space.
  ! 
  interface add
     module procedure particle_container_add
     module procedure particle_container_add_point
  end interface add

  !
  ! This is the total number of "valid" + "invalid" particles.
  !
  interface capacity
     module procedure particle_container_capacity
  end interface capacity

  interface clear
     module procedure particle_container_clear
  end interface clear

  interface print
     module procedure particle_container_print
  end interface print
  !
  ! Gives direct access to the underlying vector of particles.
  !
  ! This includes both "valid" and "invalid" particles.
  !
  ! Do NOT try to "add" or "remove" particles using this interface.
  !
  ! You can change/update valid particles using this interface, but if that
  ! causes any particles to be "owned" by another CPU you'll have to call
  ! redistribute() as well.
  !
  interface dataptr
    module procedure particle_container_dataptr
  end interface dataptr

  interface redistribute
     module procedure particle_container_redistribute
  end interface redistribute
  !
  ! These two are useful for testing purposes and/or as examples of how
  ! to write init routines and move routines.
  !
  interface init_random
     module procedure particle_container_init_random
  end interface init_random
  interface move_random
     module procedure particle_container_move_random
  end interface move_random

  interface move_advect
     module procedure particle_container_move_advect
  end interface move_advect

  interface timestamp
     module procedure particle_container_timestamp
  end interface timestamp

  private :: particle_container_reserve, particle_container_swap
  
  logical, save, private :: pVerboseFlag   = .false.
  logical, save, private :: pDebuggingFlag = .false.

  integer, save :: particle_id = 1

  !
  ! Used by checkpoint/restart.
  !
  character(len=*), parameter, private :: Hdr         = 'HDR'
  character(len=*), parameter, private :: TheData     = 'DATA'
  character(len=*), parameter, private :: ParticleDir = 'Particles'
  character(len=*), parameter, private :: Version     = 'Version_One_Dot_Zero'

contains

  subroutine particle_index(p,lev,mla,dx,prob_lo,iv)

    type(particle),   intent(in)    :: p
    integer,          intent(in)    :: lev
    type(ml_layout),  intent(in)    :: mla
    double precision, intent(in)    :: dx(:,:)
    double precision, intent(in)    :: prob_lo(:)
    integer,          intent(inout) :: iv(:)

    integer i

    do i = 1, mla%dim
       iv(i) = floor((p%pos(i)-prob_lo(i))/dx(lev,i)) + lwb(mla%la(lev)%lap%pd,i)
    end do

  end subroutine particle_index
  !
  ! A return value of true means we found the particle in our grid hierarchy.
  !
  ! A value of false means it's outside our domain.
  !
  function particle_where(p,mla,dx,prob_lo,update) result(r)

    use bl_error_module

    type(particle),   intent(inout) :: p
    type(ml_layout),  intent(inout) :: mla
    double precision, intent(in   ) :: dx(:,:)
    double precision, intent(in   ) :: prob_lo(:)

    logical, intent(in), optional :: update

    logical                         :: r, lupdate
    integer                         :: lev, iv(MAX_SPACEDIM), dm
    type(box_intersector), pointer  :: bi(:)
    type(box)                       :: bx

    lupdate = .false. ; if ( present(update) ) lupdate = update

    dm = mla%dim

    if ( lupdate ) then
       !
       ! We have a valid particle whose position has changed slightly.
       !
       ! Try to update it smartly; i.e. less costly than the whole enchilada.
       !
       call bl_assert(p%id  > 0, 'particle_where: p%id must be  > 0')
       call bl_assert(p%grd > 0, 'particle_where: p%grd must be > 0')
       call bl_assert(p%lev > 0, 'particle_where: p%lev must be > 0')

       call bl_assert(p%lev <= size(mla%la), 'particle_where: lev out of bounds')

       call bl_assert(p%grd <= nboxes(mla%la(p%lev)%lap%bxa), 'particle_where: p%grd out of bounds')

       call particle_index(p,p%lev,mla,dx,prob_lo,iv)

       if ( all(p%cell(1:dm) == iv(1:dm)) ) then
          !
          ! The particle hasn't left its cell.
          !
          r = .true.

          return
       end if

       if ( p%lev == size(mla%la) ) then

          p%cell(1:dm) = iv(1:dm)
          
          if ( contains(get_box(mla%la(p%lev)%lap%bxa,p%grd),iv(1:dm)) ) then
             !
             ! It's left its cell but's still in the same grid.
             !
             r = .true.

             return
          end if
       end if
    end if

    do lev = size(mla%la), 1, -1

       call particle_index(p,lev,mla,dx,prob_lo,iv)

       call build(bx,iv(1:dm))

       bi => layout_get_box_intersector(mla%la(lev),bx)

       if ( size(bi) > 0 ) then

          call bl_assert(size(bi) == 1, 'particle_where: should only be one box intersector')

          p%lev        = lev
          p%grd        = bi(1)%i
          p%cell(1:dm) = iv(1:dm)

          deallocate(bi)
          !
          ! Found where the particle belongs !!!
          !
          r = .true.

          return
       end if

       deallocate(bi)

    end do

    r = .false.

  end function particle_where

  subroutine particle_periodic_shift(p,mla,dx,prob_lo,prob_hi)

    use bl_error_module

    type(particle),   intent(inout) :: p
    type(ml_layout),  intent(in   ) :: mla
    double precision, intent(in)    :: dx(:,:)
    double precision, intent(in   ) :: prob_lo(:)
    double precision, intent(in   ) :: prob_hi(:)

    integer                     :: i, dm, iv(MAX_SPACEDIM)
    type(box)                   :: pd
    double precision            :: plen
    double precision, parameter :: eps = 1.0d-13

    call bl_assert(p%id > 0, 'periodic_shift: not a valid particle')

    dm = mla%dim

    pd = mla%la(p%lev)%lap%pd
    
    call particle_index(p,p%lev,mla,dx,prob_lo,iv)

    do i = 1,dm
       if ( .not. mla%pmask(i) ) cycle

       plen = (prob_hi(i) - prob_lo(i))

       if ( iv(i) > upb(pd,i) ) then

          if ( p%pos(i) == prob_hi(i) ) then
             !
             ! Don't let particles lie exactly on the domain face.
             ! Force the particle to be outside the domain so the
             ! periodic shift will bring it back inside.
             !
             p%pos(i) = p%pos(i) + eps;
          end if

          p%pos(i) = p%pos(i) - plen

       else if ( iv(i) < lwb(pd,i) ) then

          if ( p%pos(i) == prob_lo(i) ) then
             !
             ! Don't let particles lie exactly on the domain face.
             ! Force the particle to be outside the domain so the
             ! periodic shift will bring it back inside.
             !
             p%pos(i) = p%pos(i) - eps;
          end if

          p%pos(i) = p%pos(i) + plen

       end if
    end do

  end subroutine particle_periodic_shift

  subroutine particle_print(p)

    type(particle), intent(in) :: p

    print*, 'id   = ', p%id
    print*, 'cpu  = ', p%cpu
    print*, 'lev  = ', p%lev
    print*, 'grd  = ', p%grd
    print*, 'cell = ', p%cell
    print*, 'pos  = ', p%pos

  end subroutine particle_print

  subroutine particle_container_build(d)
    type(particle_container), intent(out) :: d
    if ( .not. d%built ) then
       allocate(d%d(d%size))
       call build(d%invalid)
       d%built = .true.
    end if
  end subroutine particle_container_build

  subroutine particle_container_destroy(d)
    type(particle_container), intent(inout) :: d
    if ( d%built ) then
       call particle_container_clear(d)
       call vector_clear_i(d%invalid)
       d%size  = 0
       d%d     => Null()
       d%built = .false.
    end if
  end subroutine particle_container_destroy

  pure function particle_container_empty(d) result(r)
    logical :: r
    type(particle_container), intent(in) :: d
    r = (d%size == 0)
  end function particle_container_empty

  pure function particle_container_okay(d) result(r)
    logical :: r
    type(particle_container), intent(in) :: d
    integer i, cnt
    cnt = 0
    do i = 1, capacity(d)
       if ( d%d(i)%id <= 0 ) cycle
       cnt = cnt + 1
    end do
    r = (d%size == cnt)
  end function particle_container_okay

  pure function particle_container_size(d) result(r)
    integer :: r
    type(particle_container), intent(in) :: d
    r = d%size
  end function particle_container_size

  subroutine particle_container_swap(a,b)
    type(particle_container), intent(inout) :: a,b
    type(particle_container) :: t
    t = a; a = b; b = t
  end subroutine particle_container_swap

  pure function particle_container_capacity(d) result(r)
    integer :: r
    type(particle_container), intent(in) :: d
    r = size(d%d)
  end function particle_container_capacity

  subroutine particle_container_reserve(d, size)
    type(particle_container), intent(inout) :: d
    integer,                  intent(in   ) :: size
    type(particle), pointer :: np(:)
    if ( size <= particle_container_capacity(d) ) return
    allocate(np(size))
    if ( associated(d%d) ) then
       np(1:d%size) = d%d(1:d%size)
       deallocate(d%d)
    end if
    d%d => np
  end subroutine particle_container_reserve

  subroutine particle_container_add_point(d,loc,mla,dx,prob_lo,conditional)
    type(particle_container), intent(inout) :: d
    real(kind=dp_t),          intent(in   ) :: loc(:)
    type(ml_layout),          intent(inout) :: mla
    double precision,         intent(in   ) :: dx(:,:)
    double precision,         intent(in   ) :: prob_lo(:)
    logical, optional,        intent(in   ) :: conditional

    type(particle)          :: p
    type(particle), pointer :: p_stored
    integer                 :: i,dm
    logical                 :: lconditional
    logical                 :: same_cell

    lconditional = .false. ; if ( present(conditional) ) lconditional = conditional

    dm = mla%dim

    do i=1,dm
       p%pos(i) = loc(i)
    end do

     if ( .not. particle_where(p,mla,dx,prob_lo) ) then
        call bl_error('particle_container_add_point: invalid particle')
     end if

     if ( local(mla%la(p%lev),p%grd) ) then
        !
        ! We own it.
        !
        if ( lconditional ) then
           !
           ! Make sure that we don't already have a particle in the same grid cell.
           !
           do i = 1, capacity(d)
              if ( d%d(i)%id <= 0 ) cycle

              p_stored => d%d(i)

              same_cell = all(p_stored%cell(1:dm) == p%cell(1:dm))

              if ( same_cell .and. p_stored%lev == p%lev ) then
                 !
                 ! Already have a particle here.
                 !
                 return                 
              endif
           enddo
        endif

        p%id  = get_particle_id()
        p%cpu = parallel_myproc()
        
        call add(d,p)
     end if

  end subroutine particle_container_add_point

  subroutine particle_container_add(d,v)
    type(particle_container), intent(inout) :: d
    type(particle),           intent(in   ) :: v
    if ( d%size >= particle_container_capacity(d) ) then
       call particle_container_reserve(d,max(d%size+1,2*d%size))
       d%size      = d%size + 1
       d%d(d%size) = v
    else
       if ( empty(d%invalid) ) then
          d%size      = d%size + 1
          d%d(d%size) = v
       else
          d%size               = d%size + 1
          d%d(back(d%invalid)) = v
          call pop_back(d%invalid)
       end if
    end if
  end subroutine particle_container_add

  subroutine particle_container_remove(d,i)
    type(particle_container), intent(inout) :: d
    integer,                  intent(in   ) :: i
    d%size    =  d%size - 1
    d%d(i)%id = -d%d(i)%id
    call push_back(d%invalid,i)
  end subroutine particle_container_remove

  subroutine particle_container_clear(d)
    type(particle_container), intent(inout) :: d
    d%size = 0
    if ( associated(d%d) ) then
       deallocate(d%d)
       d%d => Null()
    end if
    call clear(d%invalid)
  end subroutine particle_container_clear

  subroutine particle_container_print(d, str, valid)
    type(particle_container), intent(in)           :: d
    character (len=*),        intent(in), optional :: str
    logical,                  intent(in), optional :: valid

    logical :: lvalid
    integer :: i
    !
    ! If "valid" .eq. .true. only print valid particles.
    !
    lvalid = .false. ; if ( present(valid) ) lvalid = valid

    if ( present(str) ) print*, str
    if ( empty(d) ) then
       print*, '"Empty"'
    else
       do i = 1, d%size
          if ( lvalid ) then
             if ( d%d(i)%id > 0 ) then
                call print(d%d(i))
             end if
          else
             call print(d%d(i))
          end if
       end do
    end if
  end subroutine particle_container_print

  function particle_container_dataptr(d) result(r)
    type(particle), pointer :: r(:)
    type(particle_container), intent(in) :: d
    r => d%d
  end function particle_container_dataptr

  pure function particle_verbose() result(r)
    logical :: r
    r = pVerboseFlag
  end function particle_verbose
  subroutine particle_setverbose(val)
    logical :: val
    pVerboseFlag = val
  end subroutine particle_setverbose

  pure function particle_container_debugging() result(r)
    logical :: r
    r = pDebuggingFlag
  end function particle_container_debugging
  subroutine particle_container_setdebugging(val)
    logical :: val
    pDebuggingFlag = val
  end subroutine particle_container_setdebugging

  function get_particle_id() result(r)
    integer :: r
    r = particle_id
    particle_id = particle_id + 1
  end function get_particle_id

  subroutine set_particle_id(r)
    integer :: r
    particle_id = r
  end subroutine set_particle_id

  subroutine particle_container_init_random(particles,icnt,iseed,mla,dx,prob_lo,prob_hi)

    use mt19937_module
    use bl_error_module

    type(particle_container), intent(inout) :: particles
    integer,                  intent(in   ) :: icnt
    integer,                  intent(in   ) :: iseed
    type(ml_layout),          intent(inout) :: mla
    double precision,         intent(in   ) :: dx(:,:)
    double precision,         intent(in   ) :: prob_lo(:)
    double precision,         intent(in   ) :: prob_hi(:)

    integer          :: i, j, dm, nparticles, nparticles_tot
    double precision :: rnd, len(MAX_SPACEDIM)
    type(particle)   :: p

    dm = mla%dim

    call bl_assert(icnt  > 0, 'particle_container_init_random: icnt must be > 0')
    call bl_assert(iseed > 0, 'particle_container_init_random: iseed must be > 0')

    call bl_assert(empty(particles), 'particle_container_init_random: particle vector should be empty')

    do i = 1,dm
       len(i) = prob_hi(i) - prob_lo(i)
    end do
    !
    ! All CPUs get the same random numbers.
    !
    ! Hence they all generate the same particles.
    !
    ! But only the CPU that "owns" the particle keeps it.
    !
    call init_genrand(iseed)

    do i = 1,icnt

       do j = 1,dm
          !
          ! A random number in (0,1).
          !
          rnd = genrand_real3()

          p%pos(j) = prob_lo(j) + (rnd * len(j))

          call bl_assert(p%pos(j) < prob_hi(j), 'particle_container_init_random: particle out of bounds')
       end do

       if ( .not. particle_where(p,mla,dx,prob_lo) ) then
          call bl_error('particle_container_init_random: invalid particle')
       end if

       if ( local(mla%la(p%lev),p%grd) ) then
          !
          ! We own it.
          !
          p%id  = particle_id
          p%cpu = parallel_myproc()
          particle_id = particle_id + 1

          call add(particles,p)
       end if

    end do

    nparticles = size(particles)

    call parallel_reduce(nparticles_tot, nparticles, MPI_SUM, proc = parallel_IOProcessorNode())

    if ( parallel_IOProcessor() ) then
       print*, 'particle_container_init_random(): nparticles_tot: ', nparticles_tot
    end if
    !
    ! Note that there is no need to call redistribute() here as we've
    ! guaranteed the particles are owned by the correct CPU.
    !
  end subroutine particle_container_init_random

  subroutine particle_container_move_random(particles,mla,dx,prob_lo,prob_hi)

    use mt19937_module
    use bl_error_module

    type(particle_container), intent(inout) :: particles
    type(ml_layout),          intent(inout) :: mla
    double precision,         intent(in   ) :: dx(:,:)
    double precision,         intent(in   ) :: prob_lo(:)
    double precision,         intent(in   ) :: prob_hi(:)

    integer                 :: i, j, dm, sgn
    double precision        :: pos
    type(particle), pointer :: p

    dm = mla%dim

    if ( particle_container_debugging() ) then
       call bl_assert(ok(particles), 'particle_container_move_random: not OK on entry')
    end if

    do i = 1, capacity(particles)
       !
       ! Make sure to ignore invalid particles.
       !
       if ( particles%d(i)%id <= 0 ) cycle

       p => particles%d(i)

       do j = 1, dm
          sgn = 1

          if ( genrand_real3() >= 0.5d0 ) sgn = -1

          pos = sgn * 0.25d0 * dx(p%lev,j) * genrand_real3()

          p%pos(j) = p%pos(j) + pos
       end do
       !
       ! The particle has moved.
       !
       ! Try to put it in the right place in the hierarchy.
       !
       if ( .not. particle_where(p,mla,dx,prob_lo,update=.true.) ) then
          !
          ! Try to shift particle back across any periodic boundary.
          !
          call particle_periodic_shift(p,mla,dx,prob_lo,prob_hi)

          if ( .not. particle_where(p,mla,dx,prob_lo) ) then
             !
             ! TODO - deal with non-periodic boundary conditions !!!
             !
             print*, 'particle leaving the domain:'

             call print(p)

             call particle_container_remove(particles,i)
          end if
       end if
    end do
    !
    ! Call redistribute() to give particles to the CPU that owns'm.
    !
    call particle_container_redistribute(particles,mla,dx,prob_lo,.true.)

    if ( particle_container_debugging() ) then
       call bl_assert(ok(particles), 'particle_container_move_random: not OK on exit')
    end if

  end subroutine particle_container_move_random

  subroutine particle_container_move_advect(particles,mla,umac,dx,dt,prob_lo,prob_hi)

    use bl_error_module

    type(particle_container), intent(inout) :: particles
    type(ml_layout),          intent(inout) :: mla
    type(multifab),           intent(in   ) :: umac(:,:)
    double precision,         intent(in   ) :: dx(:,:)
    double precision,         intent(in   ) :: dt
    double precision,         intent(in   ) :: prob_lo(:)
    double precision,         intent(in   ) :: prob_hi(:)

    integer                  :: dm, d, p_id
    type(particle), pointer  :: p
    real(kind=dp_t), pointer :: ump(:,:,:,:)
    real(kind=dp_t), pointer :: vmp(:,:,:,:)
    real(kind=dp_t), pointer :: wmp(:,:,:,:)
    double precision         :: umac_lo(mla%dim), umac_hi(mla%dim)
    double precision         :: vel, delta
    integer                  :: ipass

    dm = mla%dim

    if ( particle_container_debugging() ) then
       call bl_assert(ok(particles), 'particle_container_move_advect: not OK on entry')
    end if
    !
    !-------------------------------------------------------------------------
    ! do the advection using a predictor-corrector scheme.  ipass = 1 is the
    ! predictor and updates the particle to the 1/2-time.  ipass = 2 is the 
    ! corrector and updates the particle to the final time.
    ! -------------------------------------------------------------------------
    !
    do ipass = 1,2
       do p_id = 1, capacity(particles)
          !
          ! Make sure to ignore invalid particles.
          !
          if ( particles%d(p_id)%id <= 0 ) cycle

          p => particles%d(p_id)
          
          ump => dataptr(umac(p%lev,1),local_index(umac(p%lev,1),p%grd))

          select case (dm)
          case (1)
             umac_lo(1) = ump(p%cell(1)  ,1,1,1)
             umac_hi(1) = ump(p%cell(1)+1,1,1,1)
          case (2)
             vmp => dataptr(umac(p%lev,2),local_index(umac(p%lev,2),p%grd))
             umac_lo(1) = ump(p%cell(1)  ,p%cell(2),1,1)
             umac_hi(1) = ump(p%cell(1)+1,p%cell(2),1,1)
             umac_lo(2) = vmp(p%cell(1),p%cell(2)  ,1,1)
             umac_hi(2) = vmp(p%cell(1),p%cell(2)+1,1,1)
          case (3)
             vmp => dataptr(umac(p%lev,2),local_index(umac(p%lev,2),p%grd))
             wmp => dataptr(umac(p%lev,3),local_index(umac(p%lev,3),p%grd))
             umac_lo(1) = ump(p%cell(1)  ,p%cell(2),p%cell(3),1)
             umac_hi(1) = ump(p%cell(1)+1,p%cell(2),p%cell(3),1)
             umac_lo(2) = vmp(p%cell(1),p%cell(2)  ,p%cell(3),1)
             umac_hi(2) = vmp(p%cell(1),p%cell(2)+1,p%cell(3),1)
             umac_lo(3) = wmp(p%cell(1),p%cell(2),p%cell(3)  ,1)
             umac_hi(3) = wmp(p%cell(1),p%cell(2),p%cell(3)+1,1)
          end select

          do d=1,dm

             delta = (p%pos(d) - prob_lo(d)) / dx(p%lev,d) - p%cell(d)
             vel   = umac_lo(d) + delta * (umac_hi(d) - umac_lo(d))

             if (ipass == 1) then
                !
                ! predictor stage
                !
                ! save the original particle position and update the
                ! particle to the 1/2-time position
                !
                p%origpos(d) = p%pos(d)
                
                p%pos(d) = p%pos(d) + HALF*dt*vel
             else
                ! corrector stage
                !
                ! update to the final time, using the original position
                ! and the velocity at the 1/2-time position
                !
                p%pos(d) = p%origpos(d) + dt*vel
             endif

          end do
          !
          ! The particle has moved.
          !
          ! Try to put it in the right place in the hierarchy.
          !
          if ( .not. particle_where(p,mla,dx,prob_lo,update=.true.) ) then
             !
             ! Try to shift particle back across any periodic boundary.
             !
             call particle_periodic_shift(p,mla,dx,prob_lo,prob_hi)
             
             if ( .not. particle_where(p,mla,dx,prob_lo) ) then
                !
                ! TODO - deal with non-periodic boundary conditions !!!
                !
                print*, 'particle leaving the domain:'
                
                call print(p)
                
                call particle_container_remove(particles,p_id)
             end if
          end if
       end do
       !
       ! Call redistribute() to give particles to the CPU that owns'm.
       !
       call particle_container_redistribute(particles,mla,dx,prob_lo,.true.)

    enddo

    if ( particle_container_debugging() ) then
       call bl_assert(ok(particles), 'particle_container_move_advect: not OK on exit')
    end if

  end subroutine particle_container_move_advect

  subroutine particle_container_redistribute(particles,mla,dx,prob_lo,where)

    use parallel
    use bl_error_module

    type(particle_container), intent(inout) :: particles
    type(ml_layout),          intent(inout) :: mla
    double precision,         intent(in   ) :: dx(:,:)
    double precision,         intent(in   ) :: prob_lo(:)
    !
    ! Has particle_where() been called on all the particles?
    !
    logical, intent(in), optional :: where

    type(particle)   :: p
    integer          :: i, myproc, nprocs, proc, sCnt, rCnt, iN, rN, ioff, roff, dm
    logical          :: lwhere
    double precision :: rbeg, rend, rtime

    integer, allocatable, save :: nSnd(:), nRcv(:), nSndOff(:), nRcvOff(:)
    integer, allocatable, save :: indx(:), nSnd2(:), nRcv2(:)

    double precision, allocatable, save :: SndDataR(:), RcvDataR(:)
    integer,          allocatable, save :: SndDataI(:), RcvDataI(:)
    integer,                       save :: sCntMax = 0, rCntMax = 0

    type(particle), pointer :: pp

    if ( particle_container_debugging() ) then
       call bl_assert(ok(particles), 'particle_container_redistribute: not OK on entry')
    end if

    rbeg = parallel_wtime()

    lwhere = .false. ; if ( present(where) ) lwhere = where

    if ( .not. lwhere ) then
       do i = 1, capacity(particles)
          if ( particles%d(i)%id <= 0 ) cycle
          if ( .not. particle_where(particles%d(i),mla,dx,prob_lo) ) then
             call bl_error('particle_container_redistribute: invalid particle in original vector')
          end if
       end do
    end if

    myproc = parallel_myproc()
    nprocs = parallel_nprocs()

    if ( nprocs == 1 ) return

    dm = mla%dim
    iN = 2        ! The count of integers in each particle sent/received.
    rN = 2*dm     ! The count of reals    in each particle sent/received.
                  ! The "2" here refers to "pos" and "origpos"

    if ( .not. allocated(nSnd) ) then
       allocate(indx   (0:nprocs-1)                    )
       allocate(nSnd   (0:nprocs-1),nRcv   (0:nprocs-1))
       allocate(nSnd2  (0:nprocs-1),nRcv2  (0:nprocs-1))
       allocate(nSndOff(0:nprocs-1),nRcvOff(0:nprocs-1))
    end if

    nSnd = 0
    nRcv = 0

    do i = 1, capacity(particles)
       if ( particles%d(i)%id <= 0 ) cycle

       pp => particles%d(i)

       proc = get_proc(mla%la(pp%lev),pp%grd)

       if ( proc == myproc ) cycle

       nSnd(proc) = nSnd(proc) + 1
    end do

    call bl_assert(nSnd(myproc) == 0, 'particle_container_redistribute: no sending to oneself')

    call parallel_alltoall(nRcv, nSnd, 1)

    call bl_assert(nRcv(myproc) == 0, 'particle_container_redistribute: no receiving from oneself')
    !
    ! Save off copies of nSnd and nRcv
    !
    nSnd2 = nSnd
    nRcv2 = nRcv

    sCnt = SUM(nSnd)
    rCnt = SUM(nRcv)

    if ( ( sCnt > sCntMax ) .or. ( .not. allocated(SndDataI) ) ) then
       if ( allocated(SndDataI) ) then
          deallocate(SndDataI,SndDataR)
       end if
       allocate( SndDataI (0:iN*sCnt-1), SndDataR (0:rN*sCnt-1) )
       sCntMax = sCnt
    end if

    if ( ( rCnt > rCntMax ) .or. ( .not. allocated(RcvDataI) ) ) then
       if ( allocated(RcvDataI) ) then
          deallocate(RcvDataI,RcvDataR)
       end if
       allocate( RcvDataI (0:iN*rCnt-1), RcvDataR (0:rN*rCnt-1) )
       rCntMax = rCnt
    end if
    !
    ! Now populate SndDataI and SndDataR with particle data.
    !
    nSndOff(0) = 0
    nRcvOff(0) = 0
    do i = 1, nprocs-1
       nSndOff(i) = nSndOff(i-1) + nSnd(i-1)
       nRcvOff(i) = nRcvOff(i-1) + nRcv(i-1)
    end do

    indx = nSndOff

    do i = 1, capacity(particles)
       if ( particles%d(i)%id <= 0 ) cycle

       pp => particles%d(i)

       proc = get_proc(mla%la(pp%lev),pp%grd)

       if ( proc == myproc ) cycle

       ioff = iN * indx(proc)
       roff = rN * indx(proc)

       SndDataI(ioff          )      = pp%id
       SndDataI(ioff+1        )      = pp%cpu
       SndDataR(roff:roff+dm-1)      = pp%pos(1:dm)
       SndDataR(roff+dm:roff+2*dm-1) = pp%origpos(1:dm)
       
       indx(proc) = indx(proc) + 1
    end do
    !
    ! First get the integer data in chunks of "iN".
    !
    nSnd    = iN * nSnd2
    nRcv    = iN * nRcv2
    nSndOff = iN * nSndOff
    nRcvOff = iN * nRcvOff

    call parallel_alltoall(RcvDataI, nRcv, nRcvOff, SndDataI, nSnd, nSndOff)
    !
    ! Now the real data in chunks of "rN".
    !
    nSnd = rN * nSnd2
    nRcv = rN * nRcv2

    nSndOff(0) = 0
    nRcvOff(0) = 0
    do i = 1, NProcs-1
       nSndOff(i) = nSndOff(i-1) + nSnd(i-1)
       nRcvOff(i) = nRcvOff(i-1) + nRcv(i-1)
    end do

    call parallel_alltoall(RcvDataR, nRcv, nRcvOff, SndDataR, nSnd, nSndOff)
    !
    ! Let's remove() sent particles to make space for received ones.
    !
    do i = 1, capacity(particles)
       if ( particles%d(i)%id <= 0 ) cycle

       pp => particles%d(i)

       if ( local(mla%la(pp%lev),pp%grd) ) cycle

       call remove(particles, i)
    end do

    do i = 0, rCnt-1

       ioff = iN * i
       roff = rN * i

       p%id            = RcvDataI(ioff          )
       p%cpu           = RcvDataI(ioff+1        )
       p%pos(1:dm)     = RcvDataR(roff:roff+dm-1)
       p%origpos(1:dm) = RcvDataR(roff+dm:roff+2*dm-1)
       !
       ! Got to set the members of the particle that we didn't transfer.
       !
       if ( .not. particle_where(p,mla,dx,prob_lo) ) then
          call bl_error('particle_container_redistribute: invalid particle after particle transfer')
       end if

       call add(particles,p)
    end do

    if ( particle_verbose() ) then
       rend = parallel_wtime() - rbeg

       call parallel_reduce(rtime, rend, MPI_MAX, proc = parallel_IOProcessorNode())

       if ( parallel_IOProcessor() ) then
          print*, '    particle_container_redistribute(): time: ', rtime
       end if
    end if

    if ( particle_container_debugging() ) then
       call bl_assert(ok(particles), 'particle_container_redistribute: not OK on exit')
    end if

  end subroutine particle_container_redistribute

  subroutine particle_container_checkpoint(particles,dir,mla)

    use parallel
    use fabio_module, only: fabio_mkdir, fabio_open, fabio_close, FABIO_WRONLY, &
                            fabio_write_raw_array_i, fabio_write_raw_array_d
    use bl_IO_module, only: unit_new
    use bl_error_module

    type(particle_container), intent(in   ) :: particles
    character(len=*),         intent(in   ) :: dir
    type(ml_layout),          intent(inout) :: mla

    character(len=256) :: pdir
    integer            :: i, j, k, nparticles, nparticles_tot, ioproc
    integer            :: un, nprocs, iN, dN, fd, cnt, pid, maxpid
    double precision   :: rbeg, rend, rtime

    integer,          allocatable :: isnd(:), ircv(:)
    integer, save,    allocatable :: rcvc(:), rcvd(:)
    double precision, allocatable :: dsnd(:), drcv(:)

    rbeg = parallel_wtime()

    if ( particle_container_debugging() ) then
       call bl_assert(ok(particles), 'particle_container_checkpoint: not OK on entry')
    end if

    call bl_assert(trim(dir) .ne. '' , 'particle_container_checkpoint: dir must be non-empty')

    pdir = trim(dir) // '/' // ParticleDir

    if ( parallel_IOProcessor() ) then
       call fabio_mkdir(pdir)
    end if
    !
    ! Gotta wait till the directory gets built.
    !
    call parallel_barrier()

    iN = 2       ! # of integers to snd/rcv for each particle.
    dN = mla%dim ! # of double precisions to snd/rcv for each particle.
    !
    ! Note that we do NOT checkpoint/restart the origpos member.
    !
    nprocs     = parallel_nprocs()
    ioproc     = parallel_IOProcessorNode()
    nparticles = size(particles)

    call parallel_reduce(nparticles_tot, nparticles, MPI_SUM)
    !
    ! We have to save the maximum value of particle ID.
    !
    pid = get_particle_id()
    !
    ! Got to set it back as get_particle_id() increments it.
    !
    call set_particle_id(pid)

    call parallel_reduce(maxpid, pid, MPI_MAX)

    if ( parallel_IOProcessor() ) then

       un = unit_new()

       open(unit   = un,                             &
            file   = trim(pdir) // '/' // trim(Hdr), &
            form   = 'formatted',                    &
            access = 'sequential',                   &
            status = 'replace',                      &
            action = 'write')
       !
       ! First thing written is our Checkpoint/Restart version string.
       !
       write(unit = un, fmt = '(A)') Version
       !
       ! Then dim for sanity checking.
       !
       write(unit = un, fmt = '(I1)') mla%dim
       !
       ! Then the total number of particles.
       !
       write(unit = un, fmt = '(I9)') nparticles_tot
       !
       ! Then the maximum particle ID that needs to be restored on restart.
       !
       write(unit = un, fmt = '(I9)') maxpid

       close(un)
    end if

    if ( nparticles_tot == 0 ) return

    if ( .not. allocated(rcvc) ) then
       allocate(rcvc(0:nprocs-1), rcvd(0:nprocs-1))
    end if
    !
    ! Since the number of particles is expected to be small relative to the
    ! amount of other data in the problem we'll write all the particle data
    ! into a single file.  I'll send all the data to the IO processor. We'll
    ! only write out the id and cpu of the integer data.  The rest can be
    ! rebuilt on restart() via redistribute().
    !
    if ( parallel_IOProcessor() ) then

       allocate(ircv(iN * nparticles_tot))
       !
       ! Let's open the file into which we'll stuff the particle data.
       !
       call fabio_open(fd, trim(pdir) // '/' // trim(TheData), FABIO_WRONLY)
    end if
    !
    ! We add one to the allocation here to guarantee we have at least one element.
    !
    allocate(isnd(iN * nparticles + 1))
    !
    ! Got to first send the counts of number of things to expect from each CPU.
    !
    isnd(1) = nparticles

    call parallel_gather(isnd, rcvc, 1, root = ioproc)
    !
    ! First the integers.
    !
    if ( parallel_IOProcessor() ) then
       rcvc = iN * rcvc
       rcvd = 0
       do i = 1,nprocs-1
          rcvd(i) = rcvd(i-1) + rcvc(i-1)
       end do
    end if

    cnt = 0

    do i = 1, capacity(particles)
       if ( particles%d(i)%id <= 0 ) cycle

       j           = iN * cnt + 1
       isnd(j    ) = particles%d(i)%id
       isnd(j + 1) = particles%d(i)%cpu

       cnt = cnt + 1

       call bl_assert(particles%d(i)%id > 0, 'particle_container_checkpoint: got an invalid particle')
    end do

    call bl_assert(cnt == nparticles, 'particle_container_checkpoint: particle count wrong 1')

    call parallel_gather(isnd, iN*nparticles, ircv, rcvc, rcvd, root = ioproc)

    deallocate(isnd)

    if ( parallel_IOProcessor() ) then
       !
       ! Append ircv to a file.
       !
       call fabio_write_raw_array_i(fd, ircv, iN*nparticles_tot)

       deallocate(ircv)
       !
       ! Now the real data.
       !
       allocate(drcv(dN * nparticles_tot))
    end if
    !
    ! Add one to the allocation to guarantee we have at least one element.
    !
    allocate(dsnd(dN * nparticles + 1))

    if ( parallel_IOProcessor() ) then
       rcvc = dN * rcvc
       rcvc =      rcvc / iN
       rcvd = 0
       do i = 1,nprocs-1
          rcvd(i) = rcvd(i-1) + rcvc(i-1)
       end do
    end if

    cnt = 0

    do i = 1, capacity(particles)
       if ( particles%d(i)%id <= 0 ) cycle

       j = dN * cnt
       do k = 1, dN
          dsnd(j + k) = particles%d(i)%pos(k)
       end do

       cnt = cnt + 1
    end do

    call bl_assert(cnt == nparticles, 'particle_container_checkpoint: particle count wrong 2')

    call parallel_gather(dsnd, dN*nparticles, drcv, rcvc, rcvd, root = ioproc)

    deallocate(dsnd)

    if ( parallel_IOProcessor() ) then
       !
       ! Append drcv to file.
       !
       call fabio_write_raw_array_d(fd, drcv, dN*nparticles_tot)

       call fabio_close(fd)

       deallocate(drcv)
    end if

    if ( particle_verbose() ) then
       rend = parallel_wtime() - rbeg

       call parallel_reduce(rtime, rend, MPI_MAX, proc = ioproc)

       if ( parallel_IOProcessor() ) then
          print*, '    particle_container_checkpoint(): time: ', rtime
       end if
    end if

    if ( particle_container_debugging() ) then
       call bl_assert(ok(particles), 'particle_container_checkpoint: not OK on exit')
    end if

  end subroutine particle_container_checkpoint

  subroutine particle_container_restart(particles,dir,mla,dx,prob_lo)

    use parallel
    use fabio_module, only: fabio_open, fabio_close, FABIO_RDONLY, &
                            fabio_write_raw_array_i, fabio_write_raw_array_d
    use bl_IO_module, only: unit_new
    use bl_error_module

    type(particle_container), intent(inout) :: particles
    character(len=*),         intent(in   ) :: dir
    type(ml_layout),          intent(inout) :: mla
    double precision,         intent(in   ) :: dx(:,:)
    double precision,         intent(in   ) :: prob_lo(:)

    type(particle)           :: p
    type(particle_container) :: tparticles
    character(len=256)       :: pdir, the_version_string
    integer                  :: i, j, k, dm, nparticles, ioproc
    integer                  :: un, nprocs, iN, dN, fd, maxpid
    double precision         :: rbeg, rend, rtime

    integer,          allocatable :: ircv(:)
    double precision, allocatable :: drcv(:)

    rbeg = parallel_wtime()

    if ( particle_container_debugging() ) then
       call bl_assert(ok(particles), 'particle_container_restart: not OK on entry')
    end if

    call bl_assert(empty(particles), 'particle_container_restart: particle vector must be empty')

    call bl_assert(trim(dir) .ne. '' , 'particle_container_restart: dir must be non-empty')

    pdir = trim(dir) // '/' // ParticleDir

    iN = 2       ! # of integers to snd/rcv for each particle.
    dN = mla%dim ! # of double precisions to snd/rcv for each particle.
    !
    ! Note that we do NOT checkpoint/restart the origpos member.
    !
    nprocs = parallel_nprocs()
    ioproc = parallel_IOProcessorNode()

    call build(particles)

    if ( parallel_IOProcessor() ) then

       un = unit_new()

       open(unit   = un,                             &
            file   = trim(pdir) // '/' // trim(Hdr), &
            form   = 'formatted',                    &
            access = 'sequential',                   &
            status = 'old',                          &
            action = 'read')
       !
       ! First thing read is our Checkpoint/Restart version string.
       !
       read(unit = un, fmt = '(A)') the_version_string
       !
       ! Then dim for sanity checking.
       !
       read(unit = un, fmt = '(I1)') dm
       !
       ! Then the total number of particles.
       !
       read(unit = un, fmt = '(I9)') nparticles
       !
       ! Then the maximum particle ID that needs to be restored on restart.
       !
       read(unit = un, fmt = '(I9)') maxpid

       close(un)

       call bl_assert(dm == dN, 'particle_container_restart: dimension mis-match')

       call bl_assert(nparticles >= 0, 'particle_container_restart: nparticles must be >= 0')
    end if

    call parallel_bcast(nparticles, root = ioproc)

    call parallel_bcast(maxpid, root = ioproc)
    !
    ! Don't forget to restore particle_id.
    !
    call set_particle_id(maxpid)

    if ( nparticles == 0 ) return

    if ( parallel_IOProcessor() ) then

       allocate(ircv(iN * nparticles), drcv(dN * nparticles))
       !
       ! Let's open the file into which we stuffed the particle data.
       !
       call fabio_open(fd, trim(pdir) // '/' // trim(TheData), FABIO_RDONLY)
       !
       ! And read the data.
       !
       call fabio_read_raw_array_i(fd, ircv, iN*nparticles)

       call fabio_read_raw_array_d(fd, drcv, dN*nparticles)

       call fabio_close(fd)
       !
       ! Now assemble the particles.
       !
       call particle_container_reserve(particles, nparticles)

       do i = 1, nparticles
          j     = iN * (i-1) + 1
          p%id  = ircv(j    )
          p%cpu = ircv(j + 1)

          call bl_assert(p%id > 0, 'particle_container_restart: read an invalid particle')

          j = dN * (i-1)
          do k = 1, dN
             p%pos(k) = drcv(j + k)
          end do

          call add(particles,p)
       end do

       deallocate(ircv,drcv)
    end if
    !
    ! Let redistribute() do its thing.
    !
    call redistribute(particles, mla, dx, prob_lo)

    if ( parallel_IOProcessor() ) then
       !
       ! Trim out all the unused space in particles.
       !
       call build(tparticles)

       call particle_container_reserve(tparticles, size(particles))

       do i = 1, capacity(particles)
          if ( particles%d(i)%id <= 0 ) cycle
          call add(tparticles, particles%d(i))
       end do

       call bl_assert(size(tparticles) == size(particles), 'particle_container_restart: Hmmmm')

       call swap(particles,tparticles)

       call destroy(tparticles)
    end if

    if ( particle_verbose() ) then
       rend = parallel_wtime() - rbeg

       call parallel_reduce(rtime, rend, MPI_MAX, proc = ioproc)

       if ( parallel_IOProcessor() ) then
          print*, '    particle_container_restart(): time: ', rtime
       end if
    end if

    if ( particle_container_debugging() ) then
       call bl_assert(ok(particles), 'particle_container_restart: not OK on exit')
    end if

  end subroutine particle_container_restart

  subroutine particle_container_timestamp(particles,basename,mf,idx,names,time,vel_mf)

    use parallel
    use bl_IO_module, only: unit_new
    use fabio_module, only: fabio_unlink_if_empty
    use bl_error_module

    type(particle_container), intent(in) :: particles
    character(len=*),         intent(in) :: basename
    type(multifab),           intent(in) :: mf(:)
    integer,                  intent(in) :: idx(:)
    character(len=*),         intent(in) :: names(:)
    double precision,         intent(in) :: time
    type(multifab), optional, intent(in) :: vel_mf(:)

    integer             :: un, iSet, nOutFiles, nSets, mySet, MyProc, NProcs
    integer             :: i, j, n, dm, iBuff(1), wakeUpPID, waitForPID, tag
    double precision    :: rbeg, rend, rtime
    character(len=64)   :: filename
    character(len=2)    :: index
    logical             :: itexists

    logical, save :: need_header = .true.

    type(particle), pointer :: p
    real(dp_t),     pointer :: r(:,:,:,:)
    real(dp_t),     pointer :: v(:,:,:,:)

    logical :: print_vels

    logical,                       save :: first = .true.
    character(len=64),             save :: dfmtstr, hfmtstr
    real(dp_t), allocatable, save :: values(:)
    real(dp_t) :: vels(MAX_SPACEDIM)

    integer, parameter :: MaxOpenFiles = 32

    ! if the optional vel_mf is present, then we will also store velocities
    print_vels = .false. ; if (present(vel_mf)) print_vels = .true.


    dm        = mf(1)%dim
    NProcs    = parallel_nprocs()
    MyProc    = parallel_myproc()
    nOutFiles = min(MaxOpenFiles,parallel_nprocs())
    nSets     = (NProcs + (nOutFiles - 1)) / nOutFiles
    mySet     = MyProc / nOutFiles

    if ( particle_container_debugging() ) then
       call bl_assert(ok(particles), 'particle_container_timestamp: not OK on entry')
    end if

    rbeg = parallel_wtime()

    if ( first) then

       first = .false.
       !
       ! We assume that for any given run of a program that what goes into
       ! the timestamp file doesn't change.  For this reason we'll build the
       ! required format string for the timestamp file only once on the first
       ! call to this routine.
       !
       ! We have to build up the format string we want to use.
       !
       ! We always output these per each particle we own:
       !
       !   id cpu pos(1:dm) time
       !
       ! These are then followed by size(idx) double precision values pulled
       ! from index idx(i) from the mf for each particle.
       !
       call bl_assert(size(idx) <= 99, 'timestamp: size(idx) too big')

       if (print_vels) then
          write(unit=index, fmt='(I2)') (size(idx) + 2*dm + 1)
       else
          write(unit=index, fmt='(I2)') (size(idx) + dm + 1)
       endif

       index = adjustl(index)

       dfmtstr = '(2x, I10, 1X, I10, ' // trim(index) // '(1X, G17.10))'
       hfmtstr = '(a2, a10, 1x, a10, ' // trim(index) // '(1x, a17))'

       allocate(values(size(idx)))

       do i = 1, size(idx)
          call bl_assert(idx(i) >= 1 .and. idx(i) <= mf(1)%nc, 'timestamp: invalid indices into mf')
       end do

    end if
    !
    ! A simple sanity check.
    !
    call bl_assert(size(idx) == size(values), 'timestamp: size(idx) inconsistent across calls')

    do iSet = 0, nSets

       if (mySet == iSet) then

          if ( size(particles) > 0 ) then
             !
             ! Gotta build up the file name ...
             !
             write(unit=index, fmt='(i2.2)') mod(MyProc,nOutFiles)

             filename = trim(basename) // '_' // index
             !
             ! Fortran is somewhat braindead in the ways you're allowed to open
             ! files.  There doesn't appear to be any way to open a file for
             ! appending that will also create the file if it doesn't exist.
             !
             inquire(file=filename, exist=itexists)

             un = unit_new()

             if ( itexists ) then
                open(unit     = un,             &
                     file     = trim(filename), &
                     form     = 'formatted',    &
                     access   = 'sequential',   &
                     status   = 'old',          &
                     position = 'append',       &
                     action   = 'write')
             else
                open(unit     = un,             &
                     file     = trim(filename), &
                     form     = 'formatted',    &
                     access   = 'sequential',   &
                     status   = 'new',          &
                     action   = 'write')
             end if

             ! write the header listing the column names
             if (need_header) then
                write (unit=un, fmt="(a1)") " "

                select case (dm)
                case (1)
                   if (print_vels) then
                      write(unit=un, fmt=trim(hfmtstr)) "# ", "part-ID", "origin-CPU", "x", "time", &
                           "vx", &
                           (names(n), n=1, size(idx))
                   else
                      write(unit=un, fmt=trim(hfmtstr)) "# ", "part-ID", "origin-CPU", "x", "time", &
                           (names(n), n=1, size(idx))   
                   endif
                case (2)
                   if (print_vels) then
                      write(unit=un, fmt=trim(hfmtstr)) "# ", "part-ID", "origin-CPU", "x", "y", "time", &
                           "vx", "vy", &
                           (names(n), n=1, size(idx))
                   else
                      write(unit=un, fmt=trim(hfmtstr)) "# ", "part-ID", "origin-CPU", "x", "y", "time", &
                           (names(n), n=1, size(idx))
                   endif
                case (3)
                   if (print_vels) then
                      write(unit=un, fmt=trim(hfmtstr)) "# ", "part-ID", "origin-CPU", "x", "y", "z", "time", &
                           "vx", "vy", "vz", &
                           (names(n), n=1, size(idx))                
                   else
                      write(unit=un, fmt=trim(hfmtstr)) "# ", "part-ID", "origin-CPU", "x", "y", "z", "time", &
                        (names(n), n=1, size(idx))                
                   endif
                end select
                
                need_header = .false.
             endif


             do i = 1, capacity(particles)
                if ( particles%d(i)%id <= 0 ) cycle

                p => particles%d(i)
                !
                ! Populate values
                !
                do j = 1, size(idx)
                   r => dataptr(mf(p%lev),local_index(mf(p%lev),p%grd))

                   select case (dm)
                   case (1)
                      values(j) = r(p%cell(1),1,1,idx(j))
                   case (2)
                      values(j) = r(p%cell(1),p%cell(2),1,idx(j))
                   case (3)
                      values(j) = r(p%cell(1),p%cell(2),p%cell(3),idx(j))
                   end select
                end do

                if (print_vels) then
                   v => dataptr(vel_mf(p%lev),local_index(vel_mf(p%lev),p%grd))
                   select case (dm)
                   case (1)
                      vels(1:dm) = v(p%cell(1),1,1,1:dm)
                   case (2)
                      vels(1:dm) = v(p%cell(1),p%cell(2),1,1:dm)
                   case (3)
                      vels(1:dm) = v(p%cell(1),p%cell(2),p%cell(3),1:dm)
                   end select
                endif

                !
                ! Write to the file.
                !
                if (print_vels) then
                   write(unit=un, fmt=trim(dfmtstr)) p%id, p%cpu, p%pos(1:dm), time, vels(1:dm), values
                else
                   write(unit=un, fmt=trim(dfmtstr)) p%id, p%cpu, p%pos(1:dm), time, values
                endif
             end do

             close(un)
          end if
          !
          ! Wake up next writer.
          !
          iBuff(1)  = 0
          wakeUpPID = (MyProc + nOutFiles)
          tag       = mod(MyProc,nOutFiles)

          if (wakeUpPID < NProcs) then
             call parallel_send(iBuff, wakeUpPID, tag)
          end if
       end if

       if (mySet == (iSet + 1)) then
          !
          ! Next set waits.
          !
          waitForPID = (MyProc - nOutFiles)
          tag        = mod(MyProc,nOutFiles)

          call parallel_recv(iBuff, waitForPID, tag)
       end if
    end do

    if ( particle_verbose() ) then
       rend = parallel_wtime() - rbeg

       call parallel_reduce(rtime, rend, MPI_MAX, proc = parallel_IOProcessorNode())

       if ( parallel_IOProcessor() ) then
          print*, '    particle_container_timestamp(): time: ', rtime
       end if
    end if

    if ( particle_container_debugging() ) then
       call bl_assert(ok(particles), 'particle_container_timestamp: not OK on exit')
    end if

  end subroutine particle_container_timestamp

  subroutine make_particle_count(mla,plotdata,comp,particles)

    use bl_error_module

    type(ml_layout        ),  intent(in)    :: mla
    type(particle_container), intent(inout) :: particles
    type(multifab),           intent(inout) :: plotdata(:)
    integer                 , intent(in   ) :: comp

    type(particle), pointer  :: p
    real(kind=dp_t), pointer :: ptr(:,:,:,:)
    integer                  :: i,dm

    dm = mla%dim

    do i = 1, capacity(particles)

       if ( particles%d(i)%id <= 0 ) cycle

       p => particles%d(i)

       ptr => dataptr(plotdata(p%lev),local_index(plotdata(p%lev),p%grd))

       select case (dm)
       case (1)
          ptr(p%cell(1),1,1,comp) = ptr(p%cell(1),1,1,comp)  + 1
       case (2)
          ptr(p%cell(1),p%cell(2),1,comp) = ptr(p%cell(1),p%cell(2),1,comp)  + 1
       case (3)
          ptr(p%cell(1),p%cell(2),p%cell(3),comp) = ptr(p%cell(1),p%cell(2),p%cell(3),comp) + 1
       end select

       call bl_assert(particles%d(i)%id > 0, 'make_particle_count: got an invalid particle')
    end do

  end subroutine make_particle_count


  function particle_global_numparticles(particles) result(num_tot)

    type(particle_container), intent(inout) :: particles

    integer  :: num_tot

    integer :: nparticles, nparticles_tot

    nparticles = size(particles)
    
    call parallel_reduce(nparticles_tot, nparticles, MPI_SUM)
     
    num_tot = nparticles_tot

    return
  end function particle_global_numparticles        

end module particle_module
