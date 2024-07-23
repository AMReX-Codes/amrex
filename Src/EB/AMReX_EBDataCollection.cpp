
#include <AMReX_EBDataCollection.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_MultiCutFab.H>

#include <AMReX_EB2_Level.H>
#include <algorithm>
#include <utility>

namespace amrex {

EBDataCollection::EBDataCollection (const EB2::Level& a_level,
                                    const Geometry& a_geom,
                                    const BoxArray& a_ba_in,
                                    const DistributionMapping& a_dm,
                                    Vector<int> a_ngrow, EBSupport a_support)
    : m_ngrow(std::move(a_ngrow)),
      m_support(a_support),
      m_geom(a_geom)
{
    // The BoxArray argument may not be cell-centered BoxArray.
    const BoxArray& a_ba = amrex::convert(a_ba_in, IntVect::TheZeroVector());

    if (m_support >= EBSupport::basic)
    {
        m_cellflags = new FabArray<EBCellFlagFab>(a_ba, a_dm, 1, m_ngrow[0], MFInfo(),
                                                  DefaultFabFactory<EBCellFlagFab>());
        a_level.fillEBCellFlag(*m_cellflags, m_geom);
        m_levelset = new MultiFab(amrex::convert(a_ba,IntVect::TheUnitVector()), a_dm,
                                  1, m_ngrow[0], MFInfo(), FArrayBoxFactory());
        a_level.fillLevelSet(*m_levelset, m_geom);
    }

    if (m_support >= EBSupport::volume)
    {
        AMREX_ALWAYS_ASSERT(m_ngrow[1] <= m_ngrow[0]);

        m_volfrac = new MultiFab(a_ba, a_dm, 1, m_ngrow[1], MFInfo(), FArrayBoxFactory());
        a_level.fillVolFrac(*m_volfrac, m_geom);

        m_centroid = new MultiCutFab(a_ba, a_dm, AMREX_SPACEDIM, m_ngrow[1], *m_cellflags);
        a_level.fillCentroid(*m_centroid, m_geom);
    }

    if (m_support == EBSupport::full)
    {
        AMREX_ALWAYS_ASSERT(m_ngrow[2] <= m_ngrow[0]);

        const int ng = m_ngrow[2];

        m_bndrycent = new MultiCutFab(a_ba, a_dm, AMREX_SPACEDIM, ng, *m_cellflags);
        a_level.fillBndryCent(*m_bndrycent, m_geom);

        m_bndryarea = new MultiCutFab(a_ba, a_dm, 1, ng, *m_cellflags);
        a_level.fillBndryArea(*m_bndryarea, m_geom);

        m_bndrynorm = new MultiCutFab(a_ba, a_dm, AMREX_SPACEDIM, ng, *m_cellflags);
        a_level.fillBndryNorm(*m_bndrynorm, m_geom);

        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            const BoxArray& faceba = amrex::convert(a_ba, IntVect::TheDimensionVector(idim));
            m_areafrac[idim] = new MultiCutFab(faceba, a_dm, 1, ng, *m_cellflags);
            m_facecent[idim] = new MultiCutFab(faceba, a_dm, AMREX_SPACEDIM-1, ng, *m_cellflags);
            IntVect edge_type{1}; edge_type[idim] = 0;
            m_edgecent[idim] = new MultiCutFab(amrex::convert(a_ba, edge_type), a_dm,
                                               1, ng, *m_cellflags);
        }

        a_level.fillAreaFrac(m_areafrac, m_geom);
        a_level.fillFaceCent(m_facecent, m_geom);
        a_level.fillEdgeCent(m_edgecent, m_geom);
    }

    if (! a_level.hasEBInfo()) {
        m_cutcellmask = new iMultiFab(a_ba, a_dm, 1, 0, MFInfo(),
                                      DefaultFabFactory<IArrayBox>());
        a_level.fillCutCellMask(*m_cutcellmask, m_geom);
    }

    if (! a_level.isAllRegular()) {
        extendDataOutsideDomain(a_level.nGrowVect());
    }
}

void EBDataCollection::extendDataOutsideDomain (IntVect const& level_ng)
{
    if (m_cellflags == nullptr) { return; }

    int const ngrow = m_ngrow[0];
    Box const& data_domain = amrex::grow(m_geom.Domain(), ngrow);

    Box level_domain = m_geom.Domain();
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (m_geom.isPeriodic(idim)) {
            level_domain.grow(idim, ngrow);
        } else {
            level_domain.grow(idim, level_ng[idim]);
        }
    }

    if (level_domain.contains(data_domain)) { return; }

    Box const& level_nodal_domain = amrex::surroundingNodes(level_domain);
    Array<Box,AMREX_SPACEDIM> lev_ap_domain
        {AMREX_D_DECL(amrex::surroundingNodes(level_domain,0),
                      amrex::surroundingNodes(level_domain,1),
                      amrex::surroundingNodes(level_domain,2))};

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*m_cellflags); mfi.isValid(); ++mfi) {
        auto const current_type = (*m_cellflags)[mfi].getType();
        // Because the default value for EBCellFlagFab is regular cell, we
        // can skip a regular fab.
        if (current_type == FabType::regular) { continue; }

        Box const& bx = mfi.fabbox();
        if (! level_domain.contains(bx)) {
            Box const& nbx = amrex::surroundingNodes(bx);
            auto const& ls_a = m_levelset->array(mfi);
            auto const& flag_a = m_cellflags->array(mfi);
            Array4<Real> vfrc_a;
            if (m_volfrac) {
                vfrc_a = m_volfrac->array(mfi);
            }
            amrex::ParallelFor(nbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                if (! level_nodal_domain.contains(i,j,k)) {
                    int ii = std::clamp(i, level_nodal_domain.smallEnd(0),
                                           level_nodal_domain.bigEnd(0));
                    int jj = std::clamp(j, level_nodal_domain.smallEnd(1),
                                           level_nodal_domain.bigEnd(1));
#if (AMREX_SPACEDIM > 2)
                    int kk = std::clamp(k, level_nodal_domain.smallEnd(2),
                                           level_nodal_domain.bigEnd(2));
#else
                    int kk = 0;
#endif
                    ls_a(i,j,k) = ls_a(ii,jj,kk);
                }
            });
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                if (! level_domain.contains(i,j,k)) {
                    EBCellFlag flag;
                    int ii = std::clamp(i, level_domain.smallEnd(0),
                                           level_domain.bigEnd(0));
                    int jj = std::clamp(j, level_domain.smallEnd(1),
                                           level_domain.bigEnd(1));
#if (AMREX_SPACEDIM > 2)
                    int kk = std::clamp(k, level_domain.smallEnd(2),
                                           level_domain.bigEnd(2));
#else
                    int kk = 0;
#endif
                    if (flag_a(ii,jj,kk).isCovered()) {
                        flag.setCovered();
                    } else if (flag_a(ii,jj,kk).isRegular()) {
                        flag.setRegular();
                    } else {
                        int ncov = 0;
#if (AMREX_SPACEDIM > 2)
                        for (int ko = k; ko <= k+1; ++ko) {
#else
                        for (int ko = k; ko <= k  ; ++ko) {
#endif
                        for (int jo = j; jo <= j+1; ++jo) {
                        for (int io = i; io <= i+1; ++io) {
                            if (ls_a(io,jo,ko) >= Real(0.0)) { ++ncov; }
                        }}}
                        constexpr int nnodes = (AMREX_SPACEDIM == 2) ? 4 : 8;
                        if (ncov == nnodes) {
                            flag.setCovered();
                        } else if (ncov == 0) {
                            flag.setRegular();
                        } else {
                            flag.setSingleValued();
                        }
                    }
                    if (flag.isCovered()) {
                        flag_a(i,j,k).setCovered();
                        if (vfrc_a && vfrc_a.contains(i,j,k)) {
                            vfrc_a(i,j,k) = Real(0.0);
                        }
                    } else if (flag.isRegular()) {
                        flag_a(i,j,k).setRegular();
                        if (vfrc_a && vfrc_a.contains(i,j,k)) {
                            vfrc_a(i,j,k) = Real(1.0);
                        }
                    } else {
                        flag_a(i,j,k).setSingleValued();
                        if (vfrc_a && vfrc_a.contains(i,j,k)) {
                            vfrc_a(i,j,k) = vfrc_a(ii,jj,kk);
                        }
                    }
                }
            });

            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                if (m_areafrac[idim] && m_areafrac[idim]->ok(mfi)) {
                    auto const& ap = m_areafrac[idim]->array(mfi);
                    Box apbx = Box(ap);
                    if (apbx.smallEnd(idim) == nbx.smallEnd(idim)) {
                        apbx.growLo(idim,-1);
                    }
                    if (apbx.bigEnd(idim) == nbx.bigEnd(idim)) {
                        apbx.growHi(idim,-1);
                    }
                    auto const& lev_apidim_domain = lev_ap_domain[idim];
                    Dim3 const& off = IntVect::TheDimensionVector(idim).dim3();
                    amrex::ParallelFor(apbx,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        if (! lev_apidim_domain.contains(i,j,k)) {
                            if (flag_a(i-off.x,j-off.y,k-off.z).isCovered() ||
                                flag_a(i,j,k).isCovered())
                            {
                                ap(i,j,k) = Real(0.0);
                            }
                        }
                    });
                }
            }

            if ((*m_cellflags)[mfi].getType(mfi.validbox()) != FabType::singlevalued) {
                // If it's already a cut fab in the valid region, the change in
                // ghost cells will not change the fab types.
                (*m_cellflags)[mfi].resetType(ngrow);
            }
        }
    }
}

EBDataCollection::~EBDataCollection ()
{
    delete m_cellflags;
    delete m_levelset;
    delete m_volfrac;
    delete m_centroid;
    delete m_bndrycent;
    delete m_bndrynorm;
    delete m_bndryarea;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        delete m_areafrac[idim];
        delete m_facecent[idim];
        delete m_edgecent[idim];
    }
    delete m_cutcellmask;
}

const FabArray<EBCellFlagFab>&
EBDataCollection::getMultiEBCellFlagFab () const
{
    AMREX_ASSERT(m_cellflags != nullptr);
    return *m_cellflags;
}

const MultiFab&
EBDataCollection::getLevelSet () const
{
    AMREX_ASSERT(m_levelset != nullptr);
    return *m_levelset;
}

const MultiFab&
EBDataCollection::getVolFrac () const
{
    AMREX_ASSERT(m_volfrac != nullptr);
    return *m_volfrac;
}

const MultiCutFab&
EBDataCollection::getCentroid () const
{
    AMREX_ASSERT(m_centroid != nullptr);
    return *m_centroid;
}

const MultiCutFab&
EBDataCollection::getBndryCent () const
{
    AMREX_ASSERT(m_bndrycent != nullptr);
    return *m_bndrycent;
}

const MultiCutFab&
EBDataCollection::getBndryArea () const
{
    AMREX_ASSERT(m_bndryarea != nullptr);
    return *m_bndryarea;
}

Array<const MultiCutFab*, AMREX_SPACEDIM>
EBDataCollection::getAreaFrac () const
{
    AMREX_ASSERT(m_areafrac[0] != nullptr);
    return {AMREX_D_DECL(m_areafrac[0], m_areafrac[1], m_areafrac[2])};
}

Array<const MultiCutFab*, AMREX_SPACEDIM>
EBDataCollection::getFaceCent () const
{
    AMREX_ASSERT(m_facecent[0] != nullptr);
    return {AMREX_D_DECL(m_facecent[0], m_facecent[1], m_facecent[2])};
}

Array<const MultiCutFab*, AMREX_SPACEDIM>
EBDataCollection::getEdgeCent () const
{
    AMREX_ASSERT(m_edgecent[0] != nullptr);
    return {AMREX_D_DECL(m_edgecent[0], m_edgecent[1], m_edgecent[2])};
}

const MultiCutFab&
EBDataCollection::getBndryNormal () const
{
    AMREX_ASSERT(m_bndrynorm != nullptr);
    return *m_bndrynorm;
}

const iMultiFab*
EBDataCollection::getCutCellMask () const
{
    return m_cutcellmask;
}

}
