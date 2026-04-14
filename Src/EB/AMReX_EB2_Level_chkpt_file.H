#ifndef AMREX_EB2_LEVEL_CHKPT_FILE_H_
#define AMREX_EB2_LEVEL_CHKPT_FILE_H_
#include <AMReX_Config.H>

#include <AMReX_EB2_Level.H>
#include <AMReX_EB_chkpt_file.H>

/**
 * \file AMReX_EB2_Level_chkpt_file.H
 *
 * EB2 level implementation that sources its data from checkpoint archives.
 */

namespace amrex::EB2 {

//! EB level specialized for reading checkpoint EB data.
class ChkptFileLevel
    : public GShopLevel<ChkptFile>
{
public:

    //! Construct a level from the provided checkpoint file data.
    ChkptFileLevel (IndexSpace const* is, ChkptFile const& chkpt_file, const Geometry& geom,
                    int max_grid_size, int ngrow, bool extend_domain_face);

    //! Coarsening constructor used to build levels below an already loaded one.
    ChkptFileLevel (IndexSpace const* is, int ilev, int max_grid_size, int ngrow,
                    const Geometry& geom, ChkptFileLevel& fineLevel);

    //! Fill fine-level checkpoint buffers for CUDA execution.
    void define_fine_chkpt_file (ChkptFile const& chkpt_file,
                                 Geometry const& geom, int max_grid_size, int ngrow,
                                 bool extend_domain_face);

    //! Finalize EBCellFlag connectivity adjustments after loading from disk.
    void finalize_cell_flags ();
};

}

#endif
