
#include <WarpXCrseFineBndry.H>
#include <WarpX.H>

using namespace amrex;

WarpXCrseFineBndry::WarpXCrseFineBndry (const amrex::BoxArray& crse_ba,
                                        const amrex::BoxArray& fine_ba,
                                        const amrex::DistributionMapping& crse_dm,
                                        const amrex::DistributionMapping& fine_dm,
                                        const amrex::Geometry& crse_gm,
                                        const amrex::Geometry& fine_gm,
                                        const amrex::IntVect& ref_ratio,
                                        int ngrow)
{
    MultiFab fmf(fine_ba, fine_dm, 1, 0, MFInfo().SetAlloc(false));
    // ngrow+1 is needed so the can include particles that could deposit current on the face.
    const FabArrayBase::CFinfo& cfinfo = FabArrayBase::TheCFinfo(fmf, fine_gm, ngrow+1);

    ba_fine_layout = cfinfo.ba_cfb;
    dm_fine_layout = cfinfo.dm_cfb;

    BoxList bl(fine_ba.ixType());
    Array<int> iproc;
    for (int i = 0, N = ba_fine_layout.size(); i < N; ++i)
    {
        const Box& fbx = ba_fine_layout[i];
        const Box& cbx = amrex::coarsen(fbx, ref_ratio);
        const auto& isects = crse_ba.intersections(cbx);
        for (const auto& isec : isects)
        {
            const Box& b = amrex::refine(isec.second, ref_ratio);
            bl.push_back(b);
            iproc.push_back(crse_dm[isec.first]);

            
        }
    }

    ba_crse_layout.define(bl);
    dm_crse_layout.define(iproc);
    
    jx_crse_layout.define(amrex::convert(ba_crse_layout,WarpX::jx_nodal_flag),
                          dm_crse_layout, 1, ngrow);
    jy_crse_layout.define(amrex::convert(ba_crse_layout,WarpX::jy_nodal_flag),
                          dm_crse_layout, 1, ngrow);
    jz_crse_layout.define(amrex::convert(ba_crse_layout,WarpX::jz_nodal_flag),
                          dm_crse_layout, 1, ngrow);
}
