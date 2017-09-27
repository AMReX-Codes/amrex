
#include <AMReX_EBTower.H>
#include <AMReX_EBISLevel.H>
#include <AMReX_EBLevel.H>

namespace amrex {

EBTower* EBTower::m_instance = nullptr;

void
EBTower::Build ()
{
    if (!m_instance) {
        m_instance = new EBTower();
    }
}

void
EBTower::Destroy ()
{
    delete m_instance;
}

EBTower::EBTower ()
{
    const EBIndexSpace* ebis = AMReX_EBIS::instance();

    m_domains = ebis->getDomains();

    const int nlevels = m_domains.size();

    m_irregular_ba.resize(nlevels);
    m_covered_ba.resize(nlevels);

    for (int lev = 0; lev < nlevels; ++lev)
    {
        const auto& ebisl = ebis->getEBISLevel(lev);
        const auto& graph = ebisl.getGraph();

        Array<Box> covb;
        Array<Box> irrb;

        for (MFIter mfi(graph); mfi.isValid(); ++mfi)
        {
            const auto& g = graph[mfi];
            if (g.hasIrregular()) {
                irrb.push_back(mfi.validbox());
            } else if (g.isAllCovered()) {
                covb.push_back(mfi.validbox());
            }
        }

        amrex::AllGatherBoxes(covb);
        amrex::AllGatherBoxes(irrb);

        if (covb.size() > 0) {
            m_covered_ba[lev].define(BoxList{std::move(covb)});
        }

        if (irrb.size() > 0)
        {
            m_irregular_ba[lev].define(BoxList{std::move(irrb)});
        }
    }

    m_cellflags.resize(nlevels);

    for (int lev = 0; lev < nlevels; ++lev)
    {
        if (!m_irregular_ba[lev].empty()) 
        {
            const BoxArray& ba = m_irregular_ba[lev];
            DistributionMapping dm {ba};
            EBLevel eblevel(ba, dm, m_domains[lev], 1);

            m_cellflags[lev].define(ba, dm, 1, 0);
            m_cellflags[lev].copy(eblevel.getMultiEBCellFlagFab());
        }
    }
}

EBTower::~EBTower ()
{
}

}

