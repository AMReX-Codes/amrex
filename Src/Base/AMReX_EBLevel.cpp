
#include <AMReX_EBLevel.H>
#include <type_traits>

namespace amrex {

EBLevel::EBLevel ()
{
}

EBLevel::~EBLevel ()
{
}

EBLevel::EBLevel (const BoxArray& ba, const DistributionMapping& dm, const Box& domain, const int ng)
    : EBLevelGrid(ba,dm,domain,ng),
      m_flags(ba, dm, 1, ng)
{
    static_assert(sizeof(EBCellFlag) == 4, "sizeof EBCellFlag != 4");
    static_assert(std::is_standard_layout<EBCellFlag>::value == true, "EBCellFlag is not pod");

    for (MFIter mfi(m_flags); mfi.isValid(); ++mfi)
    {
        auto& fab = m_flags[mfi];
        const Box& bx = fab.box();
        const EBISBox& ebis = m_ebisl[mfi];
        for (BoxIterator bi(bx); bi.ok(); ++bi)
        {
            const IntVect& iv = bi();
            auto& cellflag = fab(iv);
            if (ebis.isRegular(iv)) {
                cellflag.setRegular();
            } else if (ebis.isCovered(iv)) {
                cellflag.setCovered();
            } else if (ebis.isMultiValued(iv)) {
                cellflag.setMultiValued();
                amrex::Abort("EBLevel: multi-value not supported yet");
            } else {
                cellflag.setSingleValued();
            }
        }
        const Box& ibx = amrex::grow(bx,-1);
        for (BoxIterator bi(bx); bi.ok(); ++bi)
        {
            const IntVect& iv = bi();
            auto& cellflag = fab(iv);
            cellflag.setDisconnected();
            cellflag.setConnected(IntVect::TheZeroVector());
            const auto& vofs = ebis.getVoFs(iv);
            for (const auto& vi : vofs)
            {
                for (int idir = 0; idir < AMREX_SPACEDIM; ++idir)
                {
                    for (SideIterator sit; sit.ok(); ++sit)
                    {
//                        const auto& 
                    }
                }
            }
        }        
    }
}

}
