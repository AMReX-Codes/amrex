
#include <AMReX_EBTower.H>
#include <AMReX_EBISLevel.H>

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
    const Array<Box>& domains = ebis->getDomains();

    const int nlevels = domains.size();

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
        if (irrb.size() > 0) {
            m_irregular_ba[lev].define(BoxList{std::move(irrb)});
        }
    }
}

EBTower::~EBTower ()
{
}

}

