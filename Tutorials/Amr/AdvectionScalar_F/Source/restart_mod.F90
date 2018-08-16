module restart_module

  use iso_c_binding
  use amrex_fort_module
  use amrex_base_module 

  implicit none

 interface
     subroutine writecheckpointfile (stepno, finest_level,dt,time,pba,phi) bind(c)
       import
       implicit none
       integer(c_int), intent(in) :: stepno(*)
       integer(c_int), value :: finest_level
       real(amrex_real), intent (in) :: dt(*),time(*)
       type(c_ptr), intent(in) :: pba(*), phi(*)
     end subroutine writecheckpointfile

     subroutine readheaderandboxarraydata (finest_level,stepno,dt,time,pba,pdm) bind(c)
       import
       implicit none
       integer(c_int), intent(out) :: finest_level(*)
       integer(c_int), intent(out) :: stepno(*)
       real(amrex_real), intent(out) :: dt(*),time(*)
       type(c_ptr), intent(out) :: pba(*), pdm(*)
     end subroutine readheaderandboxarraydata

     subroutine readmultifabdata (finest_level,phi) bind(c)
       import
       implicit none
       integer (c_int), value :: finest_level 
       type(c_ptr), intent(out) :: phi(*)
     end subroutine readmultifabdata

     subroutine amrex_fi_set_boxarray (lev, pba, amrcore) bind(c)
       import
       implicit none
       type(c_ptr), value :: pba
       integer(c_int), value :: lev
       type(c_ptr), value :: amrcore
     end subroutine amrex_fi_set_boxarray

     subroutine amrex_fi_set_distromap  (lev, pdm, amrcore) bind(c)
       import
       implicit none
       type(c_ptr), value :: pdm
       integer(c_int), value :: lev
       type(c_ptr), value :: amrcore
     end subroutine amrex_fi_set_distromap

     subroutine amrex_fi_clone_boxarray (bao, bai) bind(c)
       import
       implicit none
       type(c_ptr) :: bao
       type(c_ptr), value :: bai
     end subroutine amrex_fi_clone_boxarray

     subroutine amrex_fi_set_finest_level (lev, amrcore) bind(c)
       import
       implicit none
       integer(c_int), value :: lev
       type(c_ptr), value :: amrcore
     end subroutine amrex_fi_set_finest_level


 end interface

 contains

  subroutine readcheckpointfile ()

     use amrex_amr_module
     use amr_data_module, only : t_new, phi_new, phi_old, flux_reg, stepno_vec, do_reflux, dt_vec

     implicit none
     integer :: lev, finest_level(1), ncomp
     type(c_ptr) :: pba(0:amrex_max_level), pdm(0:amrex_max_level), phi(0:amrex_max_level)
     type(amrex_boxarray) :: ba(0:amrex_max_level)
     type(amrex_distromap) :: dm(0:amrex_max_level)
     type(c_ptr) :: amrcore
     type(amrex_box) :: domain, dom
     type(amrex_geometry) :: geom(0:amrex_max_level)
     
     amrcore = amrex_get_amrcore()

     ! Dummy variables 
     domain = amrex_box([0,0,0], [1,1,1])

     do lev = 0, amrex_max_level
       call amrex_geometry_build(geom(lev), domain)
     end do

     dom=geom(0)%domain

     do lev = 0, amrex_max_level
        call amrex_boxarray_build(ba(lev), dom)
        call amrex_distromap_build(dm(lev), ba(lev))
     enddo

     pba(0:amrex_max_level)=ba(0:amrex_max_level)%p
     pdm(0:amrex_max_level)=dm(0:amrex_max_level)%p
    
     call readheaderandboxarraydata (finest_level,stepno_vec,dt_vec,t_new,pba(0:amrex_max_level),pdm(0:amrex_max_level))

     do lev=0, finest_level(1)
        ba(lev)=pba(lev)
        dm(lev)=pdm(lev)
     end do

     do lev=0, finest_level(1)
        call amrex_fi_set_boxarray(lev, ba(lev)%p, amrcore)
        call amrex_fi_set_distromap(lev, dm(lev)%p, amrcore)
     enddo

     ncomp=1

     do lev=0, finest_level(1)
        call amrex_multifab_build(phi_new(lev), ba(lev), dm(lev), ncomp, 0)
        call amrex_multifab_build(phi_old(lev), ba(lev), dm(lev), ncomp, 0)
        if (lev > 0 .and. do_reflux) then
           call amrex_fluxregister_build(flux_reg(lev), ba(lev), dm(lev), amrex_ref_ratio(lev-1), lev, ncomp)
        end if
      enddo

      phi(0:amrex_max_level) = phi_new(0:amrex_max_level)%p
      call readmultifabdata(finest_level(1),phi(0:amrex_max_level))

      call amrex_fi_set_finest_level(finest_level(1),amrcore)
	
 end subroutine readcheckpointfile 

end module restart_module
