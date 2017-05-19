
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
: m_fine_geom(fine_gm)
{
    MultiFab fmf(fine_ba, fine_dm, 1, 0, MFInfo().SetAlloc(false));
    // ngrow+1 is needed so the can include particles that could deposit current on the face.
    const FabArrayBase::CFinfo& cfinfo = FabArrayBase::TheCFinfo(fmf, fine_gm, ngrow+1);

    BoxArray ba_fine_layout = cfinfo.ba_cfb;
    DistributionMapping dm_fine_layout = cfinfo.dm_cfb;

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

    BoxArray ba_crse_layout(bl);
    DistributionMapping dm_crse_layout(iproc);
    
    jx_crse_layout.define(amrex::convert(ba_crse_layout,WarpX::jx_nodal_flag),
                          dm_crse_layout, 1, ngrow);
    jy_crse_layout.define(amrex::convert(ba_crse_layout,WarpX::jy_nodal_flag),
                          dm_crse_layout, 1, ngrow);
    jz_crse_layout.define(amrex::convert(ba_crse_layout,WarpX::jz_nodal_flag),
                          dm_crse_layout, 1, ngrow);
}


Array<std::pair<Box,std::array<FArrayBox*,3> > >
WarpXCrseFineBndry::intersections(const Box& bx)
{
    BL_ASSERT(bx.cellCentered());

    Array<std::pair<Box,std::array<FArrayBox*,3> > > isects;

    unsigned char flag = MFIter::AllBoxes;
    for (MFIter mfi(jx_crse_layout, flag); mfi.isValid(); ++mfi)
    {
        Box cc = amrex::enclosedCells(mfi.validbox());
        cc &= bx;
        if (cc.ok())
        {
            isects.push_back({cc, { &jx_crse_layout[mfi],
                                    &jy_crse_layout[mfi],
                                    &jz_crse_layout[mfi] } });
        }
    }
    return isects;
}

void
WarpXCrseFineBndry::AddCurrent (MultiFab& jx, MultiFab& jy, MultiFab& jz) const
{
    const auto& period = m_fine_geom.periodicity();
    const int ng = jx_crse_layout.nGrow();
    jx.copy(jx_crse_layout, 0, 0, 1, ng, 0, period, FabArrayBase::ADD);
    jy.copy(jy_crse_layout, 0, 0, 1, ng, 0, period, FabArrayBase::ADD);
    jz.copy(jz_crse_layout, 0, 0, 1, ng, 0, period, FabArrayBase::ADD);
}

void
WarpXCrseFineBndry::setVal (Real v)
{
    jx_crse_layout.setVal(v);
    jy_crse_layout.setVal(v);
    jz_crse_layout.setVal(v);
}
