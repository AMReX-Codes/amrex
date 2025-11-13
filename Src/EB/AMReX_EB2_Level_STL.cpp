#include <AMReX_EB2_Level_STL.H>

namespace amrex::EB2 {

STLLevel::STLLevel (IndexSpace const* is, STLtools const& stl_tools, const Geometry& geom,
                    int max_grid_size, int ngrow, bool extend_domain_face, int num_crse_opt,
                    bool support_mvmc)
    : GShopLevel<STLtools>(is, geom)
{
    BL_PROFILE("EB2::STLLevel()-fine");

#if (AMREX_SPACEDIM == 3)
    bool test_marching_cubes = false;
    {
        ParmParse pp("eb2");
        pp.query("test_marching_cubes", test_marching_cubes);
    }
    if (support_mvmc && test_marching_cubes) {
        define_fine_mvmc(stl_tools, geom, max_grid_size, ngrow, extend_domain_face, num_crse_opt);
    } else
#endif
    {
        if (amrex::Verbose() && support_mvmc) {
            amrex::Warning("STLlevel: support_mvmc = true is not supported yet");
        }
        define_fine(stl_tools, geom, max_grid_size, ngrow, extend_domain_face, num_crse_opt);
    }
}

STLLevel::STLLevel (IndexSpace const* is, int ilev, int max_grid_size, int ngrow,
                    const Geometry& geom, STLLevel& fineLevel)
    : GShopLevel<STLtools>(is, ilev, max_grid_size, ngrow, geom, fineLevel)
{}

}
