#ifndef AMREX_EB2_LEVEL_STL_H_
#define AMREX_EB2_LEVEL_STL_H_
#include <AMReX_Config.H>

#include <AMReX_EB2_Level.H>
#include <AMReX_EB_STL_utils.H>

/**
 * \file AMReX_EB2_Level_STL.H
 *
 * EB2 level builder backed by triangulated STL geometry.
 */

namespace amrex::EB2 {

//! EB2 level that samples geometry from STL files.
class STLLevel
    : public GShopLevel<STLtools>
{
public:

    /**
     * \brief Construct a level from the triangulation stored in \p stl_tools.
     */
    STLLevel (IndexSpace const* is, STLtools const& stl_tools, const Geometry& geom,
              int max_grid_size, int ngrow, bool extend_domain_face,
              int num_crse_opt, bool support_mvmc, bool finest_level);

    //! Coarsen from an existing fine-level STL description.
    STLLevel (IndexSpace const* is, int ilev, int max_grid_size, int ngrow,
              const Geometry& geom, STLLevel& fineLevel);
};

}

#endif
