
#include <AMReX_EBDataCollection.H>
#include <AMReX_EBTower.H>
#include <AMReX_MultiFab.H>

namespace amrex {

EBDataCollection::EBDataCollection (const Geometry& a_geom,
                                    const BoxArray& a_ba,
                                    const DistributionMapping& a_dm,
                                    int a_ngrow, EBSupport a_support)
    : m_ngrow(a_ngrow),
      m_support(a_support),
      m_geom(a_geom),
      m_cellflags(new FabArray<EBCellFlagFab>(a_ba, a_dm, 1, a_ngrow))
{
    AMREX_ALWAYS_ASSERT(EBTower::get() != nullptr);

    EBTower::fillEBCellFlag(*m_cellflags, m_geom);
}

EBDataCollection::~EBDataCollection ()
{
    delete m_cellflags;
}

const FabArray<EBCellFlagFab>&
EBDataCollection::getMultiEBCellFlagFab () const
{
    return *m_cellflags;
}


}

