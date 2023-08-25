#ifndef AMREX_EB2_LEVEL_STL_H_
#define AMREX_EB2_LEVEL_STL_H_
#include <AMReX_Config.H>

#include <AMReX_EB2_Level.H>
#include <AMReX_EB_STL_utils.H>

namespace amrex::EB2 {

class STLLevel
    : public GShopLevel<STLtools>
{
public:

    STLLevel (IndexSpace const* is, STLtools const& stl_tools, const Geometry& geom,
              int max_grid_size, int ngrow, bool extend_domain_face, int num_crse_opt);

    STLLevel (IndexSpace const* is, int ilev, int max_grid_size, int ngrow,
              const Geometry& geom, STLLevel& fineLevel);

};

}

#endif
