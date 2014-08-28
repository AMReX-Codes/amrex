module bc_functions_module

  use bl_types
  use bl_space
  use bc_module

  implicit none

  private :: bc_pretty_bit

  integer, parameter :: BC_GEOM = 3

  integer, dimension(3,MAX_SPACEDIM,-1:1) , parameter :: BC_BIT = &
       reshape((/ &
        0,   9,  18, &          ! BC_BIT(:,1,-1)
        1,  10,  19, &          ! BC_BIT(:,2,-1)
        2,  11,  20, &          ! BC_BIT(:,3,-1)
        3,  12,  21, &          ! BC_BIT(:,1, 0)
        4,  13,  22, &          ! BC_BIT(:,2, 0)
        5,  14,  23, &          ! BC_BIT(:,3, 0)
        6,  15,  24, &          ! BC_BIT(:,1, 1)
        7,  16,  25, &          ! BC_BIT(:,2, 1)
        8,  17,  26  &          ! BC_BIT(:,3, 1)
       /), (/3,MAX_SPACEDIM,3/))

contains

  pure function skewed_q(mm) result(r)
    logical :: r
    integer, intent(in) :: mm(:,:,:,:)
    integer :: dim, face
    r = .false.
    do dim = 1, 3
       do face = -1, 1, 2
          r = r .or. any( ibits(mm(:,:,:,:),BC_BIT(BC_GEOM,dim,face),1) /= 0 )
          if ( r ) return
       end do
    end do
  end function skewed_q

  elemental function bc_skewed(m, dim, face) result(r)
    logical :: r
    integer, intent(in) :: m, dim, face
    r = btest(m,BC_BIT(BC_GEOM,dim,face))
  end function bc_skewed

  elemental function bc_dirichlet(m, dim, face) result(r)
    logical :: r
    integer, intent(in) :: m, dim, face
    r = btest(m,BC_BIT(BC_DIR,dim,face))
  end function bc_dirichlet

  elemental function bc_neumann(m, dim, face) result(r)
    logical :: r
    integer, intent(in) :: m, dim, face
    r = btest(m,BC_BIT(BC_NEU,dim,face))
  end function bc_neumann

  elemental function bc_interior(m, dim, face) result(r)
    logical :: r
    integer, intent(in) :: m, dim, face
    r = ( ibits(m,BC_BIT(BC_DIR,dim,face),1) == 0 .and. &
          ibits(m,BC_BIT(BC_NEU,dim,face),1) == 0 )
  end function bc_interior

  elemental function bc_pretty_bit(m, dim, face) result(r)
    character :: r
    integer, intent(in) :: m, dim, face
    if ( ibits(m,BC_BIT(BC_DIR,dim,face),1) == 1 ) then
       r = 'D'
    else if ( ibits(m,BC_BIT(BC_NEU,dim,face),1) == 1 ) then
       r = 'N'
    else
       r = 'I'
    end if
  end function bc_pretty_bit

end module bc_functions_module
