module averagedown_module

  use amrex_amr_module

  use amr_data_module, only : phi_new
  use my_amr_module

  implicit none
  private
  
  public :: averagedown, averagedownto

contains

  subroutine averagedown ()
    integer :: lev, finest_level, icomp
    finest_level = amrex_get_finest_level()

   do icomp=1,ncomp
      do lev = finest_level-1, 0, -1
       call amrex_average_down(phi_new(lev+1), phi_new(lev), amrex_geom(lev+1), amrex_geom(lev), &
            icomp, 1, amrex_ref_ratio(lev))
       enddo
    end do
  end subroutine averagedown

  subroutine averagedownto (clev)
    integer, intent(in) :: clev
    integer :: icomp
   
    do icomp=1,ncomp
	    call amrex_average_down(phi_new(clev+1), phi_new(clev), amrex_geom(clev+1), amrex_geom(clev), &
            icomp, 1, amrex_ref_ratio(clev))    
    enddo	

  end subroutine averagedownto

end module averagedown_module
