module multifab_physbc_edgevel_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use define_bc_module

  implicit none

  private

  public :: multifab_physbc_edgevel

contains

  subroutine multifab_physbc_edgevel(uedge,the_bc_level)

    use bl_prof_module

    type(multifab) , intent(inout) :: uedge(:)
    type(bc_level) , intent(in   ) :: the_bc_level

    ! Local variables
    real(kind=dp_t), pointer :: utp(:,:,:,:)
    real(kind=dp_t), pointer :: vtp(:,:,:,:)
    real(kind=dp_t), pointer :: wtp(:,:,:,:)
    integer                  :: lo(get_dim(uedge(1))),hi(get_dim(uedge(1)))
    integer                  :: i,ng_ut,dm

    type(bl_prof_timer), save :: bpt

    call build(bpt, "multifab_physbc_edgevel")

    dm = get_dim(uedge(1))
    ng_ut = nghost(uedge(1))

    do i=1,nfabs(uedge(1))
       utp => dataptr(uedge(1),i)
       lo =  lwb(get_box(uedge(1),i))
       hi =  upb(get_box(uedge(1),i))
       select case (dm)
       case (1)
          call physbc_edgevel_1d(utp(:,1,1,1), ng_ut, lo,hi, &
                                the_bc_level%phys_bc_level_array(i,:,:))
       case (2)
          vtp => dataptr(uedge(2),i)
          call physbc_edgevel_2d(utp(:,:,1,1), vtp(:,:,1,1), ng_ut, lo,hi, &
                                the_bc_level%phys_bc_level_array(i,:,:))
       case (3)
          vtp => dataptr(uedge(2),i)
          wtp => dataptr(uedge(3),i)
          call physbc_edgevel_3d(utp(:,:,:,1), vtp(:,:,:,1), wtp(:,:,:,1), ng_ut, &
                                lo, hi, the_bc_level%phys_bc_level_array(i,:,:))
       end select
    end do

    call destroy(bpt)

  end subroutine multifab_physbc_edgevel

  subroutine physbc_edgevel_1d(uedge,ng_ut,lo,hi,phys_bc)

    use bc_module

    integer,         intent(in   ) :: lo(:),hi(:),ng_ut
    real(kind=dp_t), intent(inout) :: uedge(lo(1)-ng_ut:)
    integer        , intent(in   ) :: phys_bc(:,:)
    
    ! Local variables
    integer :: is,ie

    is = lo(1)
    ie = hi(1)
    
    ! impose lo i side bc's
    select case(phys_bc(1,1))
    case (OUTLET)
       uedge(is-1) = uedge(is)
    case (SYMMETRY)
       uedge(is-1) = -uedge(is+1)
    case (INTERIOR, PERIODIC, INLET, SLIP_WALL, NO_SLIP_WALL)
    case  default 
       call bl_error("physbc_edgevel_1d: invalid boundary type phys_bc(1,1)")
    end select

    ! impose hi i side bc's
    select case(phys_bc(1,2))
    case (OUTLET)
       uedge(ie+2) = uedge(ie+1)
    case (SYMMETRY)
       uedge(ie+2) = -uedge(ie)
    case (INTERIOR, PERIODIC, INLET, SLIP_WALL, NO_SLIP_WALL)
    case  default
       call bl_error("physbc_edgevel_1d: invalid boundary type phys_bc(1,2)")
    end select

  end subroutine physbc_edgevel_1d


  subroutine physbc_edgevel_2d(uedge,vedge,ng_ut,lo,hi,phys_bc)

    use bc_module

    integer,         intent(in   ) :: lo(:),hi(:),ng_ut
    real(kind=dp_t), intent(inout) :: uedge(lo(1)-ng_ut:,lo(2)-ng_ut:)
    real(kind=dp_t), intent(inout) :: vedge(lo(1)-ng_ut:,lo(2)-ng_ut:)
    integer        , intent(in   ) :: phys_bc(:,:)
    
    ! Local variables
    integer :: is,ie,js,je

    is = lo(1)
    js = lo(2)
    ie = hi(1)
    je = hi(2)
    
    ! impose lo i side bc's
    select case(phys_bc(1,1))
    case (INLET, SLIP_WALL, NO_SLIP_WALL)
       vedge(is-1,:) = ZERO
    case (OUTLET)
       uedge(is-1,:) = uedge(is,:)
       vedge(is-1,:) = vedge(is,:)
    case (SYMMETRY)
       uedge(is-1,:) = -uedge(is+1,:)
       vedge(is-1,:) = vedge(is  ,:)
    case (INTERIOR, PERIODIC)
    case  default 
       call bl_error("physbc_edgevel_2d: invalid boundary type phys_bc(1,1)")
    end select

    ! impose hi i side bc's
    select case(phys_bc(1,2))
    case (INLET, SLIP_WALL, NO_SLIP_WALL)
       vedge(ie+1,:) = ZERO
    case (OUTLET)
       uedge(ie+2,:) = uedge(ie+1,:)
       vedge(ie+1,:) = vedge(ie  ,:)
    case (SYMMETRY)
       uedge(ie+2,:) = -uedge(ie,:)
       vedge(ie+1,:) = vedge(ie,:)
    case (INTERIOR, PERIODIC)
    case  default
       call bl_error("physbc_edgevel_2d: invalid boundary type phys_bc(1,2)")
    end select

    ! impose lo j side bc's
    select case(phys_bc(2,1))
    case (INLET, SLIP_WALL, NO_SLIP_WALL)
       uedge(:,js-1) = ZERO
    case (OUTLET)
       uedge(:,js-1) = uedge(:,js)
       vedge(:,js-1) = vedge(:,js)
    case (SYMMETRY)
       uedge(:,js-1) = uedge(:,js  )
       vedge(:,js-1) = -vedge(:,js+1)
    case (INTERIOR, PERIODIC)
    case  default 
       call bl_error("physbc_edgevel_2d: invalid boundary type phys_bc(2,1)")
    end select

    ! impose hi j side bc's
    select case(phys_bc(2,2))
    case (INLET, SLIP_WALL, NO_SLIP_WALL)
       uedge(:,je+1) = ZERO
    case (OUTLET)
       uedge(:,je+1) = uedge(:,je)
       vedge(:,je+2) = vedge(:,je+1)
    case (SYMMETRY)
       uedge(:,je+1) = uedge(:,je)
       vedge(:,je+2) = -vedge(:,je)
    case (INTERIOR, PERIODIC)
    case  default 
       call bl_error("physbc_edgevel_2d: invalid boundary type phys_bc(2,2)")
    end select

  end subroutine physbc_edgevel_2d


  subroutine physbc_edgevel_3d(uedge,vedge,wedge,ng_ut,lo,hi,phys_bc)

    use bc_module

    integer,         intent(in   ) :: lo(:),hi(:),ng_ut
    real(kind=dp_t), intent(inout) :: uedge(lo(1)-ng_ut:,lo(2)-ng_ut:,lo(3)-ng_ut:)
    real(kind=dp_t), intent(inout) :: vedge(lo(1)-ng_ut:,lo(2)-ng_ut:,lo(3)-ng_ut:)
    real(kind=dp_t), intent(inout) :: wedge(lo(1)-ng_ut:,lo(2)-ng_ut:,lo(3)-ng_ut:)
    integer        , intent(in   ) :: phys_bc(:,:)
    
    ! Local variables
    integer :: is,ie,js,je,ks,ke

    is = lo(1)
    js = lo(2)
    ks = lo(3)
    ie = hi(1)
    je = hi(2)
    ke = hi(3)


    ! impose lo i side bc's
    select case(phys_bc(1,1))
    case (INLET, SLIP_WALL, NO_SLIP_WALL)
       vedge(is-1,:,:) = ZERO
       wedge(is-1,:,:) = ZERO
    case (OUTLET)
       uedge(is-1,:,:) = uedge(is,:,:)
       vedge(is-1,:,:) = vedge(is,:,:)
       wedge(is-1,:,:) = wedge(is,:,:)
    case (SYMMETRY)
       uedge(is-1,:,:) = -uedge(is+1,:,:)
       vedge(is-1,:,:) = vedge(is  ,:,:)
       wedge(is-1,:,:) = wedge(is  ,:,:)
    case (INTERIOR, PERIODIC)
    case  default
       call bl_error("physbc_edgevel_3d: invalid boundary type phys_bc(1,1)")
    end select

    ! impose hi i side bc's
    select case(phys_bc(1,2))
    case (INLET, SLIP_WALL, NO_SLIP_WALL)
       vedge(ie+1,:,:) = ZERO
       wedge(ie+1,:,:) = ZERO
    case (OUTLET)
       uedge(ie+2,:,:) = uedge(ie+1,:,:)
       vedge(ie+1,:,:) = vedge(ie  ,:,:)
       wedge(ie+1,:,:) = wedge(ie  ,:,:)
    case (SYMMETRY)
       uedge(ie+2,:,:) = -uedge(ie,:,:)
       vedge(ie+1,:,:) = vedge(ie,:,:)
       wedge(ie+1,:,:) = wedge(ie,:,:)
    case (INTERIOR, PERIODIC)
    case  default
       call bl_error("physbc_edgevel_3d: invalid boundary type phys_bc(1,2)")
    end select

    ! impose lo j side bc's
    select case(phys_bc(2,1))
    case (INLET, SLIP_WALL, NO_SLIP_WALL)
       uedge(:,js-1,:) = ZERO
       wedge(:,js-1,:) = ZERO
    case (OUTLET)
       uedge(:,js-1,:) = uedge(:,js,:)
       vedge(:,js-1,:) = vedge(:,js,:)
       wedge(:,js-1,:) = wedge(:,js,:)
    case (SYMMETRY)
       uedge(:,js-1,:) = uedge(:,js  ,:)
       vedge(:,js-1,:) = -vedge(:,js+1,:)
       wedge(:,js-1,:) = wedge(:,js  ,:)
    case (INTERIOR, PERIODIC)
    case  default
       call bl_error("physbc_edgevel_3d: invalid boundary type phys_bc(2,1)")
    end select

    ! impose hi j side bc's
    select case(phys_bc(2,2))
    case (INLET, SLIP_WALL, NO_SLIP_WALL)
       uedge(:,je+1,:) = ZERO
       wedge(:,je+1,:) = ZERO
    case (OUTLET)
       uedge(:,je+1,:) = uedge(:,je  ,:)
       vedge(:,je+2,:) = vedge(:,je+1,:)
       wedge(:,je+1,:) = wedge(:,je  ,:)
    case (SYMMETRY)
       uedge(:,je+1,:) = uedge(:,je,:)
       vedge(:,je+2,:) = -vedge(:,je,:)
       wedge(:,je+1,:) = wedge(:,je,:)
    case (INTERIOR, PERIODIC)
    case  default
       call bl_error("physbc_edgevel_3d: invalid boundary type phys_bc(2,2)")
    end select

    ! impose lo k side bc's
    select case(phys_bc(3,1))
    case (INLET, SLIP_WALL, NO_SLIP_WALL)
       uedge(:,:,ks-1) = ZERO
       vedge(:,:,ks-1) = ZERO
    case (OUTLET)
       uedge(:,:,ks-1) = uedge(:,:,ks)
       vedge(:,:,ks-1) = vedge(:,:,ks)
       wedge(:,:,ks-1) = wedge(:,:,ks)
    case (SYMMETRY)
       uedge(:,:,ks-1) = uedge(:,:,ks)
       vedge(:,:,ks-1) = vedge(:,:,ks)
       wedge(:,:,ks-1) = -wedge(:,:,ks+1)
    case (INTERIOR, PERIODIC)
    case  default 
       call bl_error("physbc_edgevel_3d: invalid boundary type phys_bc(3,1)")
    end select

    ! impose hi k side bc's
    select case(phys_bc(3,2))
    case (INLET, SLIP_WALL, NO_SLIP_WALL)
       uedge(:,:,ke+1) = ZERO
       vedge(:,:,ke+1) = ZERO
    case (OUTLET)
       uedge(:,:,ke+1) = uedge(:,:,ke)
       vedge(:,:,ke+1) = vedge(:,:,ke)
       wedge(:,:,ke+2) = wedge(:,:,ke+1)
    case (SYMMETRY)
       uedge(:,:,ke+1) = uedge(:,:,ke)
       vedge(:,:,ke+1) = vedge(:,:,ke)
       wedge(:,:,ke+2) = -wedge(:,:,ke)
    case (INTERIOR, PERIODIC)
    case  default
       call bl_error("physbc_edgevel_3d: invalid boundary type phys_bc(3,2)")
    end select

  end subroutine physbc_edgevel_3d
  
end module multifab_physbc_edgevel_module
