
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
      m_geom(a_geom)
{
    if (m_support >= EBSupport::basic)
    {
        AMREX_ALWAYS_ASSERT(EBTower::get() != nullptr);

        m_cellflags = new FabArray<EBCellFlagFab>(a_ba, a_dm, 1, m_ngrow);
        EBTower::fillEBCellFlag(*m_cellflags, m_geom);

        if (m_support >= EBSupport::volume)
        {
            m_volfrac = new MultiFab(a_ba, a_dm, 1, m_ngrow);
        }
    }
}

EBDataCollection::~EBDataCollection ()
{
    delete m_cellflags;
    delete m_volfrac;
    delete m_bndrycent;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        delete m_areafrac[idim];
        delete m_facecent[idim];
    }
}

const FabArray<EBCellFlagFab>&
EBDataCollection::getMultiEBCellFlagFab () const
{
    BL_ASSERT(m_cellflags != nullptr);
    return *m_cellflags;
}

}

