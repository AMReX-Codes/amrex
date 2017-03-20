
#include <AMReX.H>
#include <AMReX_Array.H>
#include <AMReX_ParmParse.H>

using namespace amrex;

extern "C" {

    void amrex_fi_init_octree ()
    {
        ParmParse pp("amr");
        int cnt = pp.countval("max_grid_size");
        int max_grid_size;
        if (cnt == 0) {
            max_grid_size = 8;
            pp.add("max_grid_size", max_grid_size);
        } else if (cnt == 1) {
            pp.get("max_grid_size", max_grid_size);
        } else {
            amrex::Abort("amrex_fi_init_octree: must use the same max_grid_size for all levels");
        }

        int blocking_factor = 2*max_grid_size;
        pp.add("blocking_factor", blocking_factor);

        pp.add("grid_eff", 1.0);

        int max_level;
        pp.get("max_level", max_level);

        Array<int> ref_ratio(max_level, 2);
        pp.addarr("ref_ratio", ref_ratio);
    }
}
