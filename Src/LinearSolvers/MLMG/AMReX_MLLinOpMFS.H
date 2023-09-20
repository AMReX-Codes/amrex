#include <AMReX_MultiFabSet.H>
#include <AMReX_MLLinOp.H>

using namespace amrex;

template <std::size_t N>
class MLLinOPMFSetT : public MLLinOpT<MultiFabSet<N>> 
{
public:

    MultiFabSet<N> make (int amrlev, int mglev, IntVect const& ng) const;

private:
    Array<IntVect,N> m_ixtype_array;
};

template <std::size_t N>
MultiFabSet<N>
MLLinOPMFSetT<N>::make (int amrlev, int mglev, IntVect const& ng) const
{
    
    auto dm = m_dmap[amrlev][mglev];
    auto bxa = m_grids[amrlev][mglev];
    auto fact = *m_factory[amrlev][mglev];
    return MultiFabSet<N>(bxa, dm, getNComp(), ng, MFInfo(), fact);
    // MultiFabSet<N> result;
    //                       getNComp(), ng, MFInfo(),
    



    // for (int i = 0; i < N; ++i) {
    //     result[i] = FAB(amrex::convert(bxa, m_ixtype_array[i]),
    //                    m_dmap[amrlev][mglev], getNComp(), ng, MFInfo(),
    //                    *m_factory[amrlev][mglev]);
    // }
    // return result;
}