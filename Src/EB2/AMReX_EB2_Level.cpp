
#include <AMReX_EB2_Level.H>
#include <AMReX_IArrayBox.H>
#include <algorithm>

namespace amrex { namespace EB2 {

int
Level::coarsenFromFine (Level& fineLevel)
{
    const BoxArray& fine_grids = fineLevel.m_grids;
    const BoxArray& fine_covered_grids = fineLevel.m_covered_grids;
    const DistributionMapping& fine_dmap = fineLevel.m_dmap;
    m_grids = amrex::coarsen(fine_grids,2);
    m_covered_grids = amrex::coarsen(fine_covered_grids, 2);
    m_dmap = fine_dmap;

    auto const& f_levelset = fineLevel.m_levelset;
    m_levelset.define(amrex::convert(m_grids,IntVect::TheNodeVector()), m_dmap, 1, 0);
    int mvmc_error = 0;
#ifdef _OPENMP
#pragma omp parallel reduction(max:mvmc_error)
#endif
    {
#if (AMREX_SPACEDIM == 3)
        IArrayBox ncuts;
#endif
        for (MFIter mfi(m_levelset,true); mfi.isValid(); ++mfi)
            {
                const Box& ccbx = mfi.tilebox(IntVect::TheCellVector());
                const Box& ndbx = mfi.tilebox();
#if (AMREX_SPACEDIM == 3)
                ncuts.resize(amrex::surroundingNodes(ccbx),6);
#endif
                int tile_error = 0;
                amrex_eb2_check_mvmc(BL_TO_FORTRAN_BOX(ccbx),
                                     BL_TO_FORTRAN_BOX(ndbx),
                                     BL_TO_FORTRAN_ANYD(m_levelset[mfi]),
                                     BL_TO_FORTRAN_ANYD(f_levelset[mfi]),
#if (AMREX_SPACEDIM == 3)
                                     BL_TO_FORTRAN_ANYD(ncuts),
#endif
                                     &tile_error);
                mvmc_error = std::max(mvmc_error, tile_error);
            }
    }
    if (mvmc_error) return mvmc_error;
    
    const int ng = 2;
    m_cellflag.define(m_grids, m_dmap, 1, ng);
    m_volfrac.define(m_grids, m_dmap, 1, ng);
    m_centroid.define(m_grids, m_dmap, AMREX_SPACEDIM, ng);
    m_bndryarea.define(m_grids, m_dmap, 1, ng);
    m_bndrycent.define(m_grids, m_dmap, AMREX_SPACEDIM, ng);
    m_bndrynorm.define(m_grids, m_dmap, AMREX_SPACEDIM, ng);
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        m_areafrac[idim].define(amrex::convert(m_grids, IntVect::TheDimensionVector(idim)),
                                m_dmap, 1, ng);
        m_facecent[idim].define(amrex::convert(m_grids, IntVect::TheDimensionVector(idim)),
                                m_dmap, AMREX_SPACEDIM-1, ng);
    }

    const Geometry& fine_geom = fineLevel.m_geom;
    const auto& fine_period = fine_geom.periodicity();
    auto& f_cellflag = fineLevel.m_cellflag;
    f_cellflag.FillBoundary(fine_period);
    MultiFab& f_volfrac = fineLevel.m_volfrac;
    f_volfrac.FillBoundary(fine_period);
    MultiFab& f_centroid = fineLevel.m_centroid;
    f_centroid.FillBoundary(fine_period);
    MultiFab& f_bndryarea = fineLevel.m_bndryarea;
    f_bndryarea.FillBoundary(fine_period);
    MultiFab& f_bndrycent = fineLevel.m_bndrycent;
    f_bndrycent.FillBoundary(fine_period);
    MultiFab& f_bndrynorm = fineLevel.m_bndrynorm;
    f_bndrynorm.FillBoundary(fine_period);
    auto& f_areafrac = fineLevel.m_areafrac;
    auto& f_facecent = fineLevel.m_facecent;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        f_areafrac[idim].FillBoundary(fine_period);
        f_facecent[idim].FillBoundary(fine_period);
    }

    if (!fine_covered_grids.empty())
    {
        const std::vector<IntVect>& pshifts = fine_period.shiftIntVect();

#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            std::vector<std::pair<int,Box> > isects;
            for (MFIter mfi(f_volfrac); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.fabbox();
                for (const auto& iv : pshifts)
                {
                    fine_covered_grids.intersections(bx+iv, isects);
                    for (const auto& is : isects)
                    {
                        Box ibox = is.second - iv;
                        f_volfrac[mfi].setVal(0.0, ibox, 0, 1);
                        f_cellflag[mfi].setVal(EBCellFlag::TheCoveredCell(), ibox, 0, 1);
                        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                            const Box& fbx = amrex::surroundingNodes(ibox,idim);
                            f_areafrac[idim].setVal(0.0, fbx, 0, 1);
                        }
                    }
                }
            }
        }
    }

    int error = 0;
#ifdef _OPENMP
#pragma omp parallel reduction(max:error)
#endif
    for (MFIter mfi(m_volfrac,true); mfi.isValid(); ++mfi)
    {
        // default
        {
            const Box& gbx = mfi.growntilebox(2);
            m_volfrac[mfi].setVal(1.0, gbx, 0, 1);
            m_centroid[mfi].setVal(0.0, gbx, 0, AMREX_SPACEDIM);
            m_bndryarea[mfi].setVal(0.0, gbx, 0, 1);
            m_bndrycent[mfi].setVal(-1.0, gbx, 0, AMREX_SPACEDIM);
            m_bndrynorm[mfi].setVal(0.0, gbx, 0, AMREX_SPACEDIM);
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                const Box& fbx = mfi.grownnodaltilebox(idim,2);
                m_areafrac[idim][mfi].setVal(1.0, fbx, 0, 1);
                m_facecent[idim][mfi].setVal(0.0, fbx, 0, AMREX_SPACEDIM-1);
            }
        }

        const Box&  bx = mfi.tilebox();
        const Box& xbx = mfi.nodaltilebox(0);
        const Box& ybx = mfi.nodaltilebox(1);
#if (AMREX_SPACEDIM == 3)
        const Box& zbx = mfi.nodaltilebox(2);
#endif

        int tile_error = 0;
        amrex_eb2_coarsen_from_fine(BL_TO_FORTRAN_BOX( bx),
                                    BL_TO_FORTRAN_BOX(xbx),
                                    BL_TO_FORTRAN_BOX(ybx),
#if (AMREX_SPACEDIM == 3)
                                    BL_TO_FORTRAN_BOX(zbx),
#endif
                                    BL_TO_FORTRAN_ANYD(m_volfrac[mfi]),
                                    BL_TO_FORTRAN_ANYD(f_volfrac[mfi]),
                                    BL_TO_FORTRAN_ANYD(m_centroid[mfi]),
                                    BL_TO_FORTRAN_ANYD(f_centroid[mfi]),
                                    BL_TO_FORTRAN_ANYD(m_bndryarea[mfi]),
                                    BL_TO_FORTRAN_ANYD(f_bndryarea[mfi]),
                                    BL_TO_FORTRAN_ANYD(m_bndrycent[mfi]),
                                    BL_TO_FORTRAN_ANYD(f_bndrycent[mfi]),
                                    BL_TO_FORTRAN_ANYD(m_bndrynorm[mfi]),
                                    BL_TO_FORTRAN_ANYD(f_bndrynorm[mfi]),
                                    BL_TO_FORTRAN_ANYD(m_areafrac[0][mfi]),
                                    BL_TO_FORTRAN_ANYD(f_areafrac[0][mfi]),
                                    BL_TO_FORTRAN_ANYD(m_areafrac[1][mfi]),
                                    BL_TO_FORTRAN_ANYD(f_areafrac[1][mfi]),
#if (AMREX_SPACEDIM == 3)
                                    BL_TO_FORTRAN_ANYD(m_areafrac[2][mfi]),
                                    BL_TO_FORTRAN_ANYD(f_areafrac[2][mfi]),
#endif
                                    BL_TO_FORTRAN_ANYD(m_facecent[0][mfi]),
                                    BL_TO_FORTRAN_ANYD(f_facecent[0][mfi]),
                                    BL_TO_FORTRAN_ANYD(m_facecent[1][mfi]),
                                    BL_TO_FORTRAN_ANYD(f_facecent[1][mfi]),
#if (AMREX_SPACEDIM == 3)
                                    BL_TO_FORTRAN_ANYD(m_facecent[2][mfi]),
                                    BL_TO_FORTRAN_ANYD(f_facecent[2][mfi]),
#endif
                                    BL_TO_FORTRAN_ANYD(m_cellflag[mfi]),
                                    BL_TO_FORTRAN_ANYD(f_cellflag[mfi]),
                                    &tile_error);
        error = std::max(error,tile_error);
    }

    if (!error) {
        buildCellFlag();
    }

    return error;
}

void
Level::buildCellFlag ()
{
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        m_areafrac[idim].FillBoundary(0,1,{AMREX_D_DECL(1,1,1)},m_geom.periodicity());
    }

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(m_cellflag,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        amrex_eb2_build_cellflag_from_ap
            (BL_TO_FORTRAN_BOX(bx),
             BL_TO_FORTRAN_ANYD(m_areafrac[0][mfi]),
             BL_TO_FORTRAN_ANYD(m_areafrac[1][mfi]),
#if (AMREX_SPACEDIM == 3)
             BL_TO_FORTRAN_ANYD(m_areafrac[2][mfi]),
#endif
             BL_TO_FORTRAN_ANYD(m_cellflag[mfi]));
    }
}

void
Level::fillEBCellFlag (FabArray<EBCellFlagFab>& cellflag, const Geometry& geom) const
{
    if (isAllRegular()) {
        cellflag.setVal(EBCellFlag::TheDefaultCell());
        for (MFIter mfi(cellflag); mfi.isValid(); ++mfi)
        {
            auto& fab = cellflag[mfi];
            const Box& bx = fab.box();
            fab.setRegion(bx);
            fab.setType(FabType::regular);
        }
        return;
    }

    cellflag.ParallelCopy(m_cellflag,0,0,1,0,cellflag.nGrow(),geom.periodicity());

    const std::vector<IntVect>& pshifts = geom.periodicity().shiftIntVect();

    Box gdomain = geom.Domain();
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (geom.isPeriodic(idim)) {
            gdomain.grow(idim, cellflag.nGrow());
        }
    }

    auto cov_val = EBCellFlag::TheCoveredCell();
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        std::vector<std::pair<int,Box> > isects;
        for (MFIter mfi(cellflag); mfi.isValid(); ++mfi)
        {
            auto& fab = cellflag[mfi];
            const Box& bx = fab.box();
            if (!m_covered_grids.empty())
            {
                for (const auto& iv : pshifts)
                {
                    m_covered_grids.intersections(bx+iv, isects);
                    for (const auto& is : isects) {
                        fab.setVal(cov_val, is.second-iv, 0, 1);
                    }
                }
            }

            // fix type and region for each fab
            const Box& regbx = bx & gdomain;
            fab.setRegion(regbx);
            fab.setType(FabType::undefined);
            auto typ = fab.getType(regbx);
            fab.setType(typ);
        }
    }
}

void
Level::fillVolFrac (MultiFab& vfrac, const Geometry& geom) const
{
    vfrac.setVal(1.0);
    if (isAllRegular()) return;

    vfrac.ParallelCopy(m_volfrac,0,0,1,0,vfrac.nGrow(),geom.periodicity());

    const std::vector<IntVect>& pshifts = geom.periodicity().shiftIntVect();

    Real cov_val = 0.0; // for covered cells

#ifdef _OPENMP
#pragma omp parallel
#endif
    if (!m_covered_grids.empty())
    {
        std::vector<std::pair<int,Box> > isects;
        for (MFIter mfi(vfrac); mfi.isValid(); ++mfi)
        {
            auto& fab = vfrac[mfi];
            const Box& bx = fab.box();
            for (const auto& iv : pshifts)
            {
                m_covered_grids.intersections(bx+iv, isects);
                for (const auto& is : isects) {
                    fab.setVal(cov_val, is.second-iv, 0, 1);
                }
            }
        }
    }
}

void
Level::fillCentroid (MultiCutFab& centroid, const Geometry& geom) const
{
    if (isAllRegular()) {
        centroid.setVal(0.0);
        return;
    }

    MultiFab tmp(centroid.boxArray(), centroid.DistributionMap(),
                 AMREX_SPACEDIM, centroid.nGrow());
    tmp.setVal(0.0);
    tmp.ParallelCopy(m_centroid,0,0,AMREX_SPACEDIM,0,centroid.nGrow(),geom.periodicity());
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(centroid.data()); mfi.isValid(); ++mfi)
    {
        if (centroid.ok(mfi)) {
            CutFab& dst = centroid[mfi];
            const void* src = static_cast<const void*>(tmp[mfi].dataPtr());
            dst.copyFromMem(src);
        }
    }
}

void
Level::fillBndryCent (MultiCutFab& bndrycent, const Geometry& geom) const
{
    if (isAllRegular()) {
        bndrycent.setVal(-1.0);
        return;
    }

    MultiFab tmp(bndrycent.boxArray(), bndrycent.DistributionMap(),
                 bndrycent.nComp(), bndrycent.nGrow());
    tmp.setVal(-1.0);
    tmp.ParallelCopy(m_bndrycent,0,0,bndrycent.nComp(),0,bndrycent.nGrow(),geom.periodicity());
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(bndrycent.data()); mfi.isValid(); ++mfi)
    {
        if (bndrycent.ok(mfi)) {
            CutFab& dst = bndrycent[mfi];
            const void* src = static_cast<const void*>(tmp[mfi].dataPtr());
            dst.copyFromMem(src);
        }
    }
}

void
Level::fillAreaFrac (Array<MultiCutFab*,AMREX_SPACEDIM>& a_areafrac, const Geometry& geom) const
{
    if (isAllRegular()) {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            a_areafrac[idim]->setVal(1.0);
        }
        return;
    }

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        MultiCutFab& areafrac = *a_areafrac[idim];
        MultiFab tmp(areafrac.boxArray(), areafrac.DistributionMap(),
                     areafrac.nComp(), areafrac.nGrow());
        tmp.setVal(1.0);
        tmp.ParallelCopy(m_areafrac[idim],0,0,areafrac.nComp(),0,areafrac.nGrow(),geom.periodicity());

        for (MFIter mfi(areafrac.data()); mfi.isValid(); ++mfi)
        {
            if (areafrac.ok(mfi)) {
                CutFab& dst = areafrac[mfi];
                const void* src = static_cast<const void*>(tmp[mfi].dataPtr());
                dst.copyFromMem(src);
            }
        }
    }

    const std::vector<IntVect>& pshifts = geom.periodicity().shiftIntVect();

    Real cov_val = 0.0; // for covered cells
        
#ifdef _OPENMP
#pragma omp parallel
#endif
    if (!m_covered_grids.empty())
    {
        std::vector<std::pair<int,Box> > isects;
        for (MFIter mfi(a_areafrac[0]->data()); mfi.isValid(); ++mfi)
        {
            if (a_areafrac[0]->ok(mfi))
            {
                const Box& ccbx = amrex::enclosedCells((*a_areafrac[0])[mfi].box());
                for (const auto& iv : pshifts)
                {
                    m_covered_grids.intersections(ccbx+iv, isects);
                    for (const auto& is : isects) {
                        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                            const Box& fbx = amrex::surroundingNodes(is.second-iv,idim);
                            (*a_areafrac[idim])[mfi].setVal(cov_val, fbx, 0, 1);
                        }
                    }
                }
            }
        }
    }
}

void
Level::fillFaceCent (Array<MultiCutFab*,AMREX_SPACEDIM>& a_facecent, const Geometry& geom) const
{
    if (isAllRegular()) {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            a_facecent[idim]->setVal(0.0);
        }        
        return;
    }

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        MultiCutFab& facecent = *a_facecent[idim];
        MultiFab tmp(facecent.boxArray(), facecent.DistributionMap(),
                     facecent.nComp(), facecent.nGrow());
        tmp.setVal(0.0);
        tmp.ParallelCopy(m_facecent[idim],0,0,facecent.nComp(),0,facecent.nGrow(),geom.periodicity());
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(facecent.data()); mfi.isValid(); ++mfi)
        {
            if (facecent.ok(mfi)) {
                CutFab& dst = facecent[mfi];
                const void* src = static_cast<const void*>(tmp[mfi].dataPtr());
                dst.copyFromMem(src);
            }
        }
    }
}

}}
