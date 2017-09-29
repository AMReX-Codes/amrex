
#include <AMReX_EBTower.H>
#include <AMReX_EBISLevel.H>
#include <AMReX_EBLevelGrid.H>

#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

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
    BL_PROFILE("EBTower::EBTower()");

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
            EBLevelGrid eblg(ba, dm, m_domains[lev], 1);

            m_cellflags[lev].define(ba, dm, 1, 0);
            initCellFlags(lev, eblg);
        }
    }
}

EBTower::~EBTower ()
{
}

void
EBTower::initCellFlags (int lev, const EBLevelGrid& eblg)
{
    FabArray<EBCellFlagFab>& ebcf = m_cellflags[lev];
    const auto& ebisl = eblg.getEBISL();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(ebcf); mfi.isValid(); ++mfi)
    {
        auto& fab = ebcf[mfi];
        const EBISBox& ebis = ebisl[mfi];
        fab.copy(ebis.getEBGraph().getEBCellFlagFab());
        fab.setType(ebis.getEBGraph().getEBCellFlagFab().getType());
    }
}

int
EBTower::getIndex (const Box& domain) const
{
    auto bx_it = std::find(m_domains.begin(), m_domains.end(), domain);
    AMREX_ALWAYS_ASSERT(bx_it != m_domains.end());
    return std::distance(m_domains.begin(), bx_it);
}

void
EBTower::fillEBCellFlag (FabArray<EBCellFlagFab>& a_flag, const Geometry& a_geom)
{
    BL_PROFILE("EBTower::fillEBCellFlag()");

    const Box& domain = a_geom.Domain();

    int lev = m_instance->getIndex(domain);

    const auto& src_flag = m_instance->m_cellflags[lev];
    a_flag.ParallelCopy(src_flag, 0, 0, 1, 0, a_flag.nGrow(), a_geom.periodicity());

    const BoxArray& cov_ba = m_instance->m_covered_ba[lev];
    auto cov_val = EBCellFlag::TheCoveredCell();
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        std::vector<std::pair<int,Box> > isects;
        for (MFIter mfi(a_flag); mfi.isValid(); ++mfi)
        {
            auto& fab = a_flag[mfi];
            const Box& bx = fab.box() & domain;

            // covered cells
            cov_ba.intersections(bx, isects);
            for (const auto& is: isects) {
                fab.setVal(cov_val, is.second);
            }

            // fix type and region for each fab
            fab.setRegion(bx);
            fab.setType(FabType::undefined);
            auto typ = fab.getType(bx);
            fab.setType(typ);
        }
    }
}

}
