#include <AMReX_MultiFabSet.H>
#include <AMReX_MLLinOp.H>

using namespace amrex;

class MLLinOpSet : public MLLinOpT<MultiFabSet> 
{
public:

    using MFSet = typename amrex::MultiFabSet;

    MultiFabSet make (int amrlev, int mglev, IntVect const& ng) const;

    MFSet makeCoarseMG (int amrlev, int mglev, IntVect const& ng) const;

protected:
    Vector<MLLinOp> m_linop_set;
    Vector<IntVect> m_ixtype_set;
};

MultiFabSet
MLLinOpSet::make (int amrlev, int mglev, IntVect const& ng) const
{
    
    auto dm = m_dmap[amrlev][mglev];
    auto bxa = m_grids[amrlev][mglev];
    auto fact = *m_factory[amrlev][mglev];
    return MultiFabSet(bxa, dm, getNComp(), ng, MFInfo(), fact);
    // MultiFabSet<N> result;
    //                       getNComp(), ng, MFInfo(),
    



    // for (int i = 0; i < N; ++i) {
    //     result[i] = FAB(amrex::convert(bxa, m_ixtype_set[i]),
    //                    m_dmap[amrlev][mglev], getNComp(), ng, MFInfo(),
    //                    *m_factory[amrlev][mglev]);
    // }
    // return result;
}

MultiFabSet
MLLinOpSet::makeCoarseMG (int amrlev, int mglev, IntVect const& ng) const
{
    MultiFabSet mfs;
    for (int i = 0; i < getNComp(); i++) {
        mfs.setElem(i, &m_linop_set[i].makeCoarseMG(amrlev, mglev, ng));
    }
}