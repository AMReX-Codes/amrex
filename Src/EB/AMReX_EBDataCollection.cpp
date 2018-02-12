
#include <AMReX_EBDataCollection.H>
#include <AMReX_EBTower.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiCutFab.H>

namespace amrex {

EBDataCollection::EBDataCollection (const Geometry& a_geom,
                                    const BoxArray& a_ba_in,
                                    const DistributionMapping& a_dm,
                                    const Vector<int>& a_ngrow, EBSupport a_support)
    : m_ngrow(a_ngrow),
      m_support(a_support),
      m_geom(a_geom)
{
    if (m_support >= EBSupport::basic)
    {
        AMREX_ALWAYS_ASSERT(EBTower::get() != nullptr);

        // The BoxArray argument may not be cell-centered BoxArray.
        const BoxArray& a_ba = amrex::convert(a_ba_in, IntVect::TheZeroVector());

        m_cellflags = new FabArray<EBCellFlagFab>(a_ba, a_dm, 1, m_ngrow[0],
                                                  MFInfo(), DefaultFabFactory<EBCellFlagFab>());
        EBTower::fillEBCellFlag(*m_cellflags, m_geom);

        if (m_support >= EBSupport::volume)
        {
            m_volfrac = new MultiFab(a_ba, a_dm, 1, m_ngrow[1], MFInfo(), FArrayBoxFactory());
            EBTower::fillVolFrac(*m_volfrac, m_geom);
        }

        if (m_support == EBSupport::full)
        {
            const int ng = m_ngrow[2];
            m_bndrycent = new MultiCutFab(a_ba, a_dm, 3, ng, *m_cellflags);
            EBTower::fillBndryCent(*m_bndrycent, m_geom);
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                const BoxArray& faceba = amrex::convert(a_ba, IntVect::TheDimensionVector(idim));
                m_areafrac[idim] = new MultiCutFab(faceba, a_dm, 1, ng, *m_cellflags);
                m_facecent[idim] = new MultiCutFab(faceba, a_dm, AMREX_SPACEDIM-1, ng, *m_cellflags);
            }
            EBTower::fillFaceGeometry(m_areafrac, m_facecent, m_geom);
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
    AMREX_ASSERT(m_cellflags != nullptr);
    return *m_cellflags;
}

const MultiFab&
EBDataCollection::getVolFrac () const
{
    AMREX_ASSERT(m_volfrac != nullptr);
    return *m_volfrac;
}

const MultiCutFab&
EBDataCollection::getBndryCent () const
{
    AMREX_ASSERT(m_bndrycent != nullptr);
    return *m_bndrycent;
}

std::array<const MultiCutFab*, AMREX_SPACEDIM>
EBDataCollection::getAreaFrac () const
{
    AMREX_ASSERT(m_areafrac[0] != nullptr);
    return {AMREX_D_DECL(m_areafrac[0], m_areafrac[1], m_areafrac[2])};
}

std::array<const MultiCutFab*, AMREX_SPACEDIM>
EBDataCollection::getFaceCent () const
{
    AMREX_ASSERT(m_facecent[0] != nullptr);
    return {AMREX_D_DECL(m_facecent[0], m_facecent[1], m_facecent[2])};
}

}

