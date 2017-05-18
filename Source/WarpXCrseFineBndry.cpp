
#include <WarpXCrseFineBndry.H>

using namespace amrex;

WarpXCrseFineBndry::WarpXCrseFineBndry (const amrex::BoxArray& crse_ba,
                                        const amrex::BoxArray& fine_ba,
                                        const amrex::DistributionMapping& crse_dm,
                                        const amrex::DistributionMapping& fine_dm,
                                        const amrex::Geometry& crse_gm,
                                        const amrex::Geometry& fine_gm,
                                        int ngrow)
{
    MultiFab fmf(fine_ba, fine_dm, 1, 0, MFInfo().SetAlloc(false));
    // ngrow+1 is needed so the can include particles that could deposit current on the face.
    const FabArrayBase::CFinfo& cfinfo = FabArrayBase::TheCFinfo(fmf, fine_gm, ngrow+1);

    ba_fine_layout = cfinfo.ba_cfb;
    dm_fine_layout = cfinfo.dm_cfb;

    
}
