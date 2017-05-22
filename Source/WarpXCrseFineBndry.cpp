
#include <WarpXCrseFineBndry.H>
#include <WarpX.H>

using namespace amrex;

WarpXBndryForFine::WarpXBndryForFine (const amrex::BoxArray& crse_ba,
                                      const amrex::BoxArray& fine_ba,
                                      const amrex::DistributionMapping& crse_dm,
                                      const amrex::DistributionMapping& fine_dm,
                                      const amrex::Geometry& crse_gm,
                                      const amrex::Geometry& fine_gm,
                                      const amrex::IntVect& ref_ratio,
                                      int ngrow)
  : m_fine_geom(fine_gm),
    m_ref_ratio(ref_ratio),
    m_ngrow(ngrow)
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

    BoxArray ba(bl);
    DistributionMapping dm(iproc);
    
    jx.define(amrex::convert(ba,WarpX::jx_nodal_flag), dm, 1, ngrow);
    jy.define(amrex::convert(ba,WarpX::jy_nodal_flag), dm, 1, ngrow);
    jz.define(amrex::convert(ba,WarpX::jz_nodal_flag), dm, 1, ngrow);
}

bool
WarpXBndryForFine::intersects (const Box& bx) const
{
    BL_ASSERT(bx.cellCentered());
    unsigned char flag = MFIter::AllBoxes;
    for (MFIter mfi(jx, flag); mfi.isValid(); ++mfi)
    {
        const Box& cc = amrex::enclosedCells(mfi.validbox());
        if (cc.intersects(bx)) return true;
    }
    return false;
}

void
WarpXBndryForFine::addFrom (amrex::FArrayBox& jx_in, amrex::FArrayBox& jy_in, amrex::FArrayBox& jz_in)
{
    unsigned char flag = MFIter::AllBoxes;
    for (MFIter mfi(jx, flag); mfi.isValid(); ++mfi)
    {
        FArrayBox& dst = jx[mfi];
        FArrayBox& src = jx_in;
        Box bx = dst.box();
        bx &= src.box();
        if (bx.ok()) {
            dst.plus(src,bx,bx,0,0);
            src.setVal(0.0,bx,0);  // to avoid being added to another fab
        }
    }

    for (MFIter mfi(jy, flag); mfi.isValid(); ++mfi)
    {
        FArrayBox& dst = jy[mfi];
        FArrayBox& src = jy_in;
        Box bx = dst.box();
        bx &= src.box();
        if (bx.ok()) {
            dst.plus(src,bx,bx,0,0);
            src.setVal(0.0,bx,0);  // to avoid being added to another fab
        }
    }

    for (MFIter mfi(jz, flag); mfi.isValid(); ++mfi)
    {
        FArrayBox& dst = jz[mfi];
        FArrayBox& src = jz_in;
        Box bx = dst.box();
        bx &= src.box();
        if (bx.ok()) {
            dst.plus(src,bx,bx,0,0);
            src.setVal(0.0,bx,0);  // to avoid being added to another fab
        }
    }
}

void
WarpXBndryForFine::addTo (MultiFab& jx_out, MultiFab& jy_out, MultiFab& jz_out) const
{
    const auto& period = m_fine_geom.periodicity();
    const int ng = jx.nGrow();
    jx_out.copy(jx, 0, 0, 1, ng, 0, period, FabArrayBase::ADD);
    jy_out.copy(jy, 0, 0, 1, ng, 0, period, FabArrayBase::ADD);
    jz_out.copy(jz, 0, 0, 1, ng, 0, period, FabArrayBase::ADD);
}

void
WarpXBndryForFine::setVal (Real v)
{
    jx.setVal(v);
    jy.setVal(v);
    jz.setVal(v);
}

WarpXBndryForCrse::WarpXBndryForCrse (const amrex::BoxArray& crse_ba,
                                      const amrex::BoxArray& fine_ba,
                                      const amrex::DistributionMapping& crse_dm,
                                      const amrex::DistributionMapping& fine_dm,
                                      const amrex::Geometry& crse_gm,
                                      const amrex::Geometry& fine_gm,
                                      const amrex::IntVect& ref_ratio,
                                      int ngrow)
: m_crse_geom(crse_gm),
    m_ref_ratio(ref_ratio),
    m_ngrow(ngrow)
{
    MultiFab fmf(fine_ba, fine_dm, 1, 0, MFInfo().SetAlloc(false));
    const FabArrayBase::CFinfo& cfinfo = FabArrayBase::TheCFinfo(fmf, fine_gm, ngrow*ref_ratio[0]);

    BoxArray ba = cfinfo.ba_cfb;
    ba.coarsen(ref_ratio);

    jx.define(amrex::convert(ba,WarpX::jx_nodal_flag), cfinfo.dm_cfb, 1, 0);
    jy.define(amrex::convert(ba,WarpX::jy_nodal_flag), cfinfo.dm_cfb, 1, 0);
    jz.define(amrex::convert(ba,WarpX::jz_nodal_flag), cfinfo.dm_cfb, 1, 0);
}

void
WarpXBndryForCrse::setVal (Real v)
{
    jx.setVal(v);
    jy.setVal(v);
    jz.setVal(v);
}

bool
WarpXBndryForCrse::intersects (const Box& bx_in) const
{
    const Box& bx = amrex::grow(bx_in, m_ngrow);
    BL_ASSERT(bx.cellCentered());
    unsigned char flag = MFIter::AllBoxes;
    for (MFIter mfi(jx, flag); mfi.isValid(); ++mfi)
    {
        const Box& cc = amrex::enclosedCells(mfi.validbox());
        if (cc.intersects(bx)) return true;
    }
    return false;
}

void
WarpXBndryForCrse::addFrom (amrex::FArrayBox& jx_in, amrex::FArrayBox& jy_in, amrex::FArrayBox& jz_in)
{
    unsigned char flag = MFIter::AllBoxes;
    for (MFIter mfi(jx, flag); mfi.isValid(); ++mfi)
    {
        FArrayBox& dst = jx[mfi];
        FArrayBox& src = jx_in;
        Box bx = dst.box();
        bx &= src.box();
        if (bx.ok()) {
            dst.plus(src,bx,bx,0,0);
            src.setVal(0.0,bx,0);  // to avoid being added to another fab
        }
    }

    for (MFIter mfi(jy, flag); mfi.isValid(); ++mfi)
    {
        FArrayBox& dst = jy[mfi];
        FArrayBox& src = jy_in;
        Box bx = dst.box();
        bx &= src.box();
        if (bx.ok()) {
            dst.plus(src,bx,bx,0,0);
            src.setVal(0.0,bx,0);  // to avoid being added to another fab
        }
    }

    for (MFIter mfi(jz, flag); mfi.isValid(); ++mfi)
    {
        FArrayBox& dst = jz[mfi];
        FArrayBox& src = jz_in;
        Box bx = dst.box();
        bx &= src.box();
        if (bx.ok()) {
            dst.plus(src,bx,bx,0,0);
            src.setVal(0.0,bx,0);  // to avoid being added to another fab
        }
    }
}

void
WarpXBndryForCrse::addTo (MultiFab& jx_out, MultiFab& jy_out, MultiFab& jz_out) const
{
    const auto& period = m_crse_geom.periodicity();
    const int ng = jx.nGrow();
    jx_out.copy(jx, 0, 0, 1, ng, 0, period, FabArrayBase::ADD);
    jy_out.copy(jy, 0, 0, 1, ng, 0, period, FabArrayBase::ADD);
    jz_out.copy(jz, 0, 0, 1, ng, 0, period, FabArrayBase::ADD);
}
