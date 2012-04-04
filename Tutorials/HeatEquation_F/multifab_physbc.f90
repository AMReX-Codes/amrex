module multifab_physbc_module

  use multifab_module
  use define_bc_module

  implicit none

  private

  public :: multifab_physbc

contains

  subroutine multifab_physbc(s,icomp,bccomp,ncomp,the_bc_level)

    use setbc_module

    integer        , intent(in   ) :: icomp,bccomp,ncomp
    type(multifab) , intent(inout) :: s
    type(bc_level) , intent(in   ) :: the_bc_level

    ! Local
    integer                  :: lo(get_dim(s))
    integer                  :: i,ng,dm
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    
    ng = nghost(s)
    dm = get_dim(s)
    
    do i=1,nboxes(s)
       if ( multifab_remote(s,i) ) cycle
       sp => dataptr(s,i)
       lo = lwb(get_box(s,i))
       select case (dm)
       case (2)
          call setbc_2d(sp(:,:,1,icomp), lo, ng, &
                        the_bc_level%adv_bc_level_array(i,:,:,bccomp))
       case (3)
          call setbc_3d(sp(:,:,:,icomp), lo, ng, &
                        the_bc_level%adv_bc_level_array(i,:,:,bccomp))
       end select
    end do
 
end subroutine multifab_physbc

end module multifab_physbc_module
