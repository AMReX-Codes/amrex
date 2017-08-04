
#include <AMReX_EBLevel.H>
#include <type_traits>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace amrex {

EBLevel::EBLevel ()
{
}

EBLevel::~EBLevel ()
{
}

EBLevel::EBLevel (const BoxArray& ba, const DistributionMapping& dm, const Box& domain, const int ng)
    : EBLevelGrid(ba,dm,domain,ng),
      m_flags(std::make_shared<FabArray<EBFlagFab> >(ba, dm, 1, ng))
{
    static_assert(sizeof(EBCellFlag) == 4, "sizeof EBCellFlag != 4");
    static_assert(std::is_standard_layout<EBCellFlag>::value == true, "EBCellFlag is not pod");

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*m_flags); mfi.isValid(); ++mfi)
    {
        auto& fab = (*m_flags)[mfi];
        const Box& bx = fab.box() & domain;
        const EBISBox& ebis = m_ebisl[mfi];
        int nregular=0, nsingle=0, nmulti=0, ncovered=0;
        int ncells = bx.numPts();
        for (BoxIterator bi(bx); bi.ok(); ++bi)
        {
            const IntVect& iv = bi();
            auto& cellflag = fab(iv);
            if (ebis.isRegular(iv)) {
                cellflag.setRegular();
                ++nregular;
            } else if (ebis.isCovered(iv)) {
                cellflag.setCovered();
                ++ncovered;
            } else if (ebis.isMultiValued(iv)) {
                cellflag.setMultiValued(ebis.numVoFs(iv));
                ++nmulti;
                amrex::Abort("EBLevel: multi-value not supported yet");
            } else {
                cellflag.setSingleValued();
                ++nsingle;
            }
        }
        
        if (nregular == ncells) {
            fab.setType(FabType::regular);
        } else if (ncovered == ncells) {
            fab.setType(FabType::covered);
        } else if (nmulti > 0) {
            fab.setType(FabType::multivalued);
        } else {
            fab.setType(FabType::singlevalued);
        }

        const Box& ibx = amrex::grow(fab.box(),-1) & domain;
        for (BoxIterator bi(ibx); bi.ok(); ++bi)
        {
            const IntVect& iv = bi();
            auto& cellflag = fab(iv);
            cellflag.setDisconnected();
            cellflag.setConnected(IntVect::TheZeroVector());
            const auto& vofs = ebis.getVoFs(iv);
            for (const auto& vi : vofs)
            {
                for (int dir_0 = 0; dir_0 < AMREX_SPACEDIM; ++dir_0)
                {
                    int dir_1 = (dir_0+1) % AMREX_SPACEDIM;
                    int dir_2 = (dir_0+2) % AMREX_SPACEDIM;
                    for (SideIterator sit_0; sit_0.ok(); ++sit_0)
                    {
                        const auto& vofs_0 = ebis.getVoFs(vi, dir_0, sit_0(), 1);
                        for (const auto& vi_0 : vofs_0)
                        {
                            for (SideIterator sit_1; sit_1.ok(); ++sit_1)
                            {
                                const auto& vofs_1 = ebis.getVoFs(vi_0, dir_1, sit_1(), 1);
                                for (const auto& vi_1 : vofs_1)
                                {
#if (AMREX_SPACEDIM == 2)
                                    IntVect iv;
                                    iv[dir_0] = amrex::sign(sit_0());
                                    iv[dir_1] = amrex::sign(sit_1());
                                    cellflag.setConnected(iv);
#else
                                    for (SideIterator sit_2; sit_2.ok(); ++sit_2)
                                    {
                                        const auto& vofs_2 = ebis.getVoFs(vi_1, dir_2, sit_2(), 1);
                                        for (const auto& vi_2 : vofs_2)
                                        {
                                            IntVect iv;
                                            iv[dir_0] = amrex::sign(sit_0());
                                            iv[dir_1] = amrex::sign(sit_1());
                                            iv[dir_2] = amrex::sign(sit_2());
                                            cellflag.setConnected(iv);
                                        }
                                    }
#endif
                                }
                            }
                        }
                    }
                }
            }
        }        
    }
}

}
