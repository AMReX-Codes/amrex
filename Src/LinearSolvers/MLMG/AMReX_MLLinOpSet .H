#include <AMReX_FabArraySet.H>
#include <AMReX_MLLinOp.H>

using namespace amrex;

class MLLinOpSet : public MLLinOpT<MultiFabSet>
{
public:

    using MFSet = amrex::MultiFabSet;

    MFSet make (int amrlev, int mglev, IntVect const& ng) const;

    MFSet makeCoarseMG (int amrlev, int mglev, IntVect const& ng) const;

    MFSet makeCoarseAmr (int famrlev, IntVect const& ng) const;

    int getNSet () const { return m_ixtype_set.size(); }

    int getNComp () const { return m_NComp; }

protected:

    Vector<IndexType> m_ixtype_set;

    int m_NComp;
};

using MFSet = amrex::MultiFabSet;

MFSet
MLLinOpSet::make (int amrlev, int mglev, IntVect const& ng) const
{
    auto bxa = m_grids[amrlev][mglev];
    const DistributionMapping& dm = m_dmap[amrlev][mglev];
    auto fact = *m_factory[amrlev][mglev];
    return MFSet(bxa, m_ixtype_set, dm, getNComp(), ng, MFInfo(), fact);
}

// template <typename MFSet>
MFSet
MLLinOpSet::makeCoarseMG (int amrlev, int mglev, IntVect const& ng) const
{
    BoxArray cba = m_grids[amrlev][mglev];
    IntVect ratio = (amrlev > 0) ? IntVect(2) : mg_coarsen_ratio_vec[mglev];
    cba.coarsen(ratio);
    return MFSet(cba, m_ixtype_set, m_dmap[amrlev][mglev], getNComp(), ng);
}

// template <typename MFSet>
MFSet
MLLinOpSet::makeCoarseAmr (int famrlev, IntVect const& ng) const
{
    BoxArray cba = m_grids[famrlev][0];
    IntVect ratio(AMRRefRatio(famrlev-1));
    cba.coarsen(ratio);
    return MFSet(cba, m_ixtype_set, m_dmap[famrlev][0], getNComp(), ng);
}