module amrex_multifabutil_module

  use iso_c_binding
  use amrex_fort_module
  use amrex_multifab_module
  use amrex_geometry_module

  implicit none
  private

  public :: amrex_average_down, & ! volume weighted average down of cell data
       &    amrex_average_down_cell,  & ! average down of cell data
       &    amrex_average_down_node,  & ! average down of nodal data
       &    amrex_average_down_faces, & ! average down of face data
       &    amrex_average_cellcenter_to_face ! average from cell centers to faces

  interface
     subroutine amrex_fi_average_down (fmf, cmf, fgeom, cgeom, scomp, ncomp, rr) bind(c)
       import
       implicit none
       type(c_ptr), value :: fmf, cmf, fgeom, cgeom
       integer(c_int), value :: scomp, ncomp, rr
     end subroutine amrex_fi_average_down

     subroutine amrex_fi_average_down_cell_node (fmf, cmf, scomp, ncomp, rr) bind(c)
       import
       implicit none
       type(c_ptr), value :: fmf, cmf
       integer(c_int), value :: scomp, ncomp, rr
     end subroutine amrex_fi_average_down_cell_node

     subroutine amrex_fi_average_down_faces (fmf, cmf, cgeom, scomp, ncomp, rr) bind(c)
       import
       implicit none
       type(c_ptr), intent(in) :: fmf(*), cmf(*)
       type(c_ptr), value :: cgeom
       integer(c_int), value :: scomp, ncomp, rr
     end subroutine amrex_fi_average_down_faces

     subroutine amrex_fi_average_cellcenter_to_face (facemf, ccmf, geom) bind(c)
       import
       implicit none
       type(c_ptr), intent(in) :: facemf(*)
       type(c_ptr), value :: ccmf, geom
     end subroutine amrex_fi_average_cellcenter_to_face

  end interface

contains

  subroutine amrex_average_down (fmf, cmf, fgeom, cgeom, scomp, ncomp, rr)
    type(amrex_multifab), intent(in   ) :: fmf
    type(amrex_multifab), intent(inout) :: cmf
    type(amrex_geometry), intent(in) :: fgeom, cgeom
    integer, intent(in) :: scomp, ncomp, rr
    call amrex_fi_average_down(fmf%p, cmf%p, fgeom%p, cgeom%p, scomp-1, ncomp, rr)
  end subroutine amrex_average_down

  subroutine amrex_average_down_cell (fmf, cmf, scomp, ncomp, rr)
    type(amrex_multifab), intent(in   ) :: fmf
    type(amrex_multifab), intent(inout) :: cmf
    integer, intent(in) :: scomp, ncomp, rr
    call amrex_fi_average_down_cell_node(fmf%p, cmf%p, scomp-1, ncomp, rr)
  end subroutine amrex_average_down_cell

  subroutine amrex_average_down_node (fmf, cmf, scomp, ncomp, rr)
    type(amrex_multifab), intent(in   ) :: fmf
    type(amrex_multifab), intent(inout) :: cmf
    integer, intent(in) :: scomp, ncomp, rr
    call amrex_fi_average_down_cell_node(fmf%p, cmf%p, scomp-1, ncomp, rr)
  end subroutine amrex_average_down_node

  subroutine amrex_average_down_faces (fmf, cmf, cgeom, scomp, ncomp, rr)
    type(amrex_multifab), intent(in   ) :: fmf(amrex_spacedim)
    type(amrex_multifab), intent(inout) :: cmf(amrex_spacedim)
    type(amrex_geometry), intent(in) :: cgeom ! coarse level geometry
    integer, intent(in) :: scomp, ncomp, rr
    type(c_ptr) :: cp(amrex_spacedim), fp(amrex_spacedim)
    integer :: idim
    do idim = 1, amrex_spacedim
       fp(idim) = fmf(idim)%p
       cp(idim) = cmf(idim)%p
    end do
    call amrex_fi_average_down_faces(fp, cp, cgeom%p, scomp-1, ncomp, rr)
  end subroutine amrex_average_down_faces


  subroutine amrex_average_cellcenter_to_face (facemf, ccmf, geom)
    type(amrex_multifab), intent(inout) :: facemf(amrex_spacedim)
    type(amrex_multifab), intent(in) :: ccmf
    type(amrex_geometry), intent(in) :: geom
    type(c_ptr) :: cp(amrex_spacedim)
    integer :: idim
    do idim = 1, amrex_spacedim
       cp(idim) = facemf(idim)%p
    end do
    call amrex_fi_average_cellcenter_to_face(cp, ccmf%p, geom%p)
  end subroutine amrex_average_cellcenter_to_face

end module amrex_multifabutil_module
