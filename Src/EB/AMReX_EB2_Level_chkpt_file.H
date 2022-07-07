#ifndef AMREX_EB2_LEVEL_CHKPT_FILE_H_
#define AMREX_EB2_LEVEL_CHKPT_FILE_H_
#include <AMReX_Config.H>

#include <AMReX_EB2_Level.H>
#include <AMReX_EB_chkpt_file.H>

namespace amrex { namespace EB2 {

class ChkptFileLevel
    : public GShopLevel<ChkptFile>
{
public:

    ChkptFileLevel (IndexSpace const* is, ChkptFile const& chkpt_file, const Geometry& geom,
            int max_grid_size, int ngrow, bool extend_domain_face);

    ChkptFileLevel (IndexSpace const* is, int ilev, int max_grid_size, int ngrow,
            const Geometry& geom, ChkptFileLevel& fineLevel);

};

}}

#endif
