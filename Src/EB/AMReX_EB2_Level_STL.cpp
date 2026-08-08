#include <AMReX_EB2_Level_STL.H>

namespace amrex::EB2 {

STLLevel::STLLevel (IndexSpace const* is, STLtools const& stl_tools, const Geometry& geom,
                    int max_grid_size, int ngrow, bool extend_domain_face, int num_crse_opt,
                    bool support_mvmc)
    : GShopLevel<STLtools>(is, geom)
{
    BL_PROFILE("EB2::STLLevel()-fine");

    if (stl_tools.usesMarchingCubes()) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            !support_mvmc,
            "The supported marching-cubes STL generator is single-valued; "
            "use ordinary EB2::Build instead of BuildMultiValuedMultiCut");
#if (AMREX_SPACEDIM == 3)
        define_fine_marching_cubes(
            stl_tools,geom,max_grid_size,ngrow,extend_domain_face,num_crse_opt);
        return;
#else
        amrex::Abort(
            "eb2.stl_geometry_method=marching_cubes is available only in 3D");
#endif
    }

    if (amrex::Verbose() && support_mvmc) {
        amrex::Warning("STLlevel: support_mvmc = true is not supported yet");
    }
    define_fine(stl_tools,geom,max_grid_size,ngrow,extend_domain_face,num_crse_opt);
}

STLLevel::STLLevel (IndexSpace const* is, int ilev, int max_grid_size, int ngrow,
                    const Geometry& geom, STLLevel& fineLevel)
    : GShopLevel<STLtools>(is, ilev, max_grid_size, ngrow, geom, fineLevel)
{}

}
