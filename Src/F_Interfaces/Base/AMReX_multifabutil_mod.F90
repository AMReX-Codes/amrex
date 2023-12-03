module amrex_multifabutil_module

  use iso_c_binding
  use amrex_fort_module
  use amrex_multifab_module
  use amrex_geometry_module

  implicit none
  private

  public :: amrex_average_down, amrex_average_down_faces, amrex_average_cellcenter_to_face

  public :: amrex_average_down_dg, amrex_average_down_cg

  interface
     subroutine amrex_fi_average_down (fmf, cmf, fgeom, cgeom, scomp, ncomp, rr) bind(c)
       import
       implicit none
       type(c_ptr), value :: fmf, cmf, fgeom, cgeom
       integer(c_int), value :: scomp, ncomp, rr
     end subroutine amrex_fi_average_down

     subroutine amrex_fi_average_down_dg &
       ( fmf, cmf, fgeom, cgeom, scomp, ncomp, rr, nDOFX, nFine, VolumeRatio, &
         ProjectionMatrix_T, WeightsX_q ) bind(c)
       import
       implicit none
       type(c_ptr)     , value :: &
         fmf, cmf, fgeom, cgeom, ProjectionMatrix_T, WeightsX_q
       integer(c_int)  , value :: scomp, ncomp, rr, nDOFX, nFine
       real(amrex_real), value :: VolumeRatio
     end subroutine amrex_fi_average_down_dg

     subroutine amrex_fi_average_down_cg &
       ( fmf, cmf, fgeom, cgeom, scomp, ncomp, rr, nDOFX, nFine, &
         G2L, L2G, F2C ) bind(c)
       import
       implicit none
       type(c_ptr)     , value :: &
         fmf, cmf, fgeom, cgeom, G2L, L2G, F2C
       integer(c_int)  , value :: scomp, ncomp, rr, nDOFX, nFine
     end subroutine amrex_fi_average_down_cg

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

  subroutine amrex_average_down_dg &
    ( fmf, cmf, fgeom, cgeom, scomp, ncomp, rr, nDOFX, nFine, VolumeRatio, &
      ProjectionMatrix_T, WeightsX_q)
    type(amrex_multifab), intent(in   ) :: fmf
    type(amrex_multifab), intent(inout) :: cmf
    type(amrex_geometry), intent(in) :: fgeom, cgeom
    integer             , intent(in) :: scomp, ncomp, rr, nDOFX, nFine
    real(amrex_real)    , intent(in) :: VolumeRatio
    type(c_ptr)         , intent(in) :: ProjectionMatrix_T, WeightsX_q
    call amrex_fi_average_down_dg &
           ( fmf%p, cmf%p, fgeom%p, cgeom%p, scomp-1, ncomp, rr, &
             nDOFX, nFine, VolumeRatio, ProjectionMatrix_T, WeightsX_q )
  end subroutine amrex_average_down_dg

  subroutine amrex_average_down_cg &
    ( fmf, cmf, fgeom, cgeom, scomp, ncomp, rr, nDOFX, nFine, &
      G2L, L2G, F2C )
    type(amrex_multifab), intent(in   ) :: fmf
    type(amrex_multifab), intent(inout) :: cmf
    type(amrex_geometry), intent(in) :: fgeom, cgeom
    integer             , intent(in) :: scomp, ncomp, rr, nDOFX, nFine
    type(c_ptr)         , intent(in) :: G2L, L2G, F2C
    call amrex_fi_average_down_cg &
           ( fmf%p, cmf%p, fgeom%p, cgeom%p, scomp-1, ncomp, rr, &
             nDOFX, nFine, G2L, L2G, F2C )
  end subroutine amrex_average_down_cg

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
