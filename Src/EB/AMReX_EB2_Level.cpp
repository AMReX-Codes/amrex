
#include <AMReX_EB2_Level.H>
#include <AMReX_IArrayBox.H>
#include <algorithm>

namespace amrex { namespace EB2 {

void
Level::prepareForCoarsening (const Level& rhs, int max_grid_size, IntVect ngrow)
{
    BoxArray all_grids(amrex::grow(m_geom.Domain(),ngrow));
    all_grids.maxSize(max_grid_size);
    FabArray<EBCellFlagFab> cflag(all_grids, DistributionMapping{all_grids}, 1, 1);
    rhs.fillEBCellFlag(cflag, m_geom);
    
    Vector<Box> cut_boxes;
    Vector<Box> covered_boxes;

    for (MFIter mfi(cflag); mfi.isValid(); ++mfi)
    {
        FabType t = cflag[mfi].getType();
        AMREX_ASSERT(t != FabType::undefined);
        const Box& vbx = mfi.validbox();
        if (t == FabType::covered) {
            covered_boxes.push_back(vbx);
        } else if (t != FabType::regular) {
            cut_boxes.push_back(vbx);
        }
    }

    amrex::AllGatherBoxes(cut_boxes);
    amrex::AllGatherBoxes(covered_boxes);
    
    if (!covered_boxes.empty()) {
        m_covered_grids = BoxArray(BoxList(std::move(covered_boxes)));
    }

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!cut_boxes.empty(),
                                     "EB2::Level: how come there are no cut boxes?");

    m_grids = BoxArray(BoxList(std::move(cut_boxes)));
    m_dmap = DistributionMapping(m_grids);

    m_levelset.define(amrex::convert(m_grids,IntVect::TheNodeVector()), m_dmap, 1, 0);
    rhs.fillLevelSet(m_levelset, m_geom);
    
//    m_mgf.define(m_grids, m_dmap);
    const int ng = 2;
    m_cellflag.define(m_grids, m_dmap, 1, ng);
    rhs.fillEBCellFlag(m_cellflag, m_geom);

    m_volfrac.define(m_grids, m_dmap, 1, ng);
    rhs.fillVolFrac(m_volfrac, m_geom);

    m_centroid.define(m_grids, m_dmap, AMREX_SPACEDIM, ng);
    rhs.fillCentroid(m_centroid, m_geom);

    m_bndryarea.define(m_grids, m_dmap, 1, ng);
    rhs.fillBndryArea(m_bndryarea, m_geom);

    m_bndrycent.define(m_grids, m_dmap, AMREX_SPACEDIM, ng);
    rhs.fillBndryCent(m_bndrycent, m_geom);

    m_bndrynorm.define(m_grids, m_dmap, AMREX_SPACEDIM, ng);
    rhs.fillBndryNorm(m_bndrynorm, m_geom);

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        m_areafrac[idim].define(amrex::convert(m_grids, IntVect::TheDimensionVector(idim)),
                                m_dmap, 1, ng);
        m_facecent[idim].define(amrex::convert(m_grids, IntVect::TheDimensionVector(idim)),
                                m_dmap, AMREX_SPACEDIM-1, ng);
    }
    rhs.fillAreaFrac(amrex::GetArrOfPtrs(m_areafrac), m_geom);
    rhs.fillFaceCent(amrex::GetArrOfPtrs(m_facecent), m_geom);

    m_ok = true;
}

int
Level::coarsenFromFine (Level& fineLevel, bool fill_boundary)
{
    const BoxArray& fine_grids = fineLevel.m_grids;
    const BoxArray& fine_covered_grids = fineLevel.m_covered_grids;
    const DistributionMapping& fine_dmap = fineLevel.m_dmap;
    m_grids = amrex::coarsen(fine_grids,2);
    m_covered_grids = amrex::coarsen(fine_covered_grids, 2);
    m_dmap = fine_dmap;

    if (! (fine_grids.coarsenable(2,2) &&
           (fine_covered_grids.empty() || fine_covered_grids.coarsenable(2,2)))) {
        return 1;
    }

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
    {
        bool b = mvmc_error;
        ParallelDescriptor::ReduceBoolOr(b);
        mvmc_error = b;
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

    auto& f_cellflag = fineLevel.m_cellflag;
    MultiFab& f_volfrac = fineLevel.m_volfrac;
    MultiFab& f_centroid = fineLevel.m_centroid;
    MultiFab& f_bndryarea = fineLevel.m_bndryarea;
    MultiFab& f_bndrycent = fineLevel.m_bndrycent;
    MultiFab& f_bndrynorm = fineLevel.m_bndrynorm;
    auto& f_areafrac = fineLevel.m_areafrac;
    auto& f_facecent = fineLevel.m_facecent;

    if (fill_boundary)
    {
        const Geometry& fine_geom = fineLevel.m_geom;
        const auto& fine_period = fine_geom.periodicity();
        f_cellflag.FillBoundary(fine_period);
        f_volfrac.FillBoundary(fine_period);
        f_centroid.FillBoundary(fine_period);
        f_bndryarea.FillBoundary(fine_period);
        f_bndrycent.FillBoundary(fine_period);
        f_bndrynorm.FillBoundary(fine_period);
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
    {
        bool b = error;
        ParallelDescriptor::ReduceBoolOr(b);
        error = b;
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

namespace {
    void copyMultiFabToMultiCutFab (MultiCutFab& dstmf, const MultiFab& srcmf)
    {
        const int ncomp = srcmf.nComp();
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(dstmf.data()); mfi.isValid(); ++mfi)
        {
            if (dstmf.ok(mfi)) {
                const auto dstfab = dstmf.array(mfi);
                const auto srcfab = srcmf.array(mfi);
                const Box& box = mfi.fabbox();
                AMREX_HOST_DEVICE_FOR_4D (box,ncomp,i,j,k,n,
                {
                    dstfab(i,j,k,n) = srcfab(i,j,k,n);
                });
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
    fillCentroid(tmp, geom);
    copyMultiFabToMultiCutFab(centroid, tmp);
}

void
Level::fillCentroid (MultiFab& centroid, const Geometry& geom) const
{
    centroid.setVal(0.0);
    if (!isAllRegular()) {
        centroid.ParallelCopy(m_centroid,0,0,AMREX_SPACEDIM,0,centroid.nGrow(),
                              geom.periodicity());
    }
}

void
Level::fillBndryArea (MultiCutFab& bndryarea, const Geometry& geom) const
{
    if (isAllRegular()) {
        bndryarea.setVal(0.0);
        return;
    }

    MultiFab tmp(bndryarea.boxArray(), bndryarea.DistributionMap(),
                 1, bndryarea.nGrow());
    fillBndryArea(tmp, geom);
    copyMultiFabToMultiCutFab(bndryarea, tmp);
}
        
void
Level::fillBndryArea (   MultiFab& bndryarea, const Geometry& geom) const
{
    bndryarea.setVal(0.0);
    if (!isAllRegular()) {
        bndryarea.ParallelCopy(m_bndryarea,0,0,1,0,bndryarea.nGrow(),geom.periodicity());
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
    fillBndryCent(tmp, geom);
    copyMultiFabToMultiCutFab(bndrycent, tmp);
}

void
Level::fillBndryCent (MultiFab& bndrycent, const Geometry& geom) const
{
    bndrycent.setVal(-1.0);
    if (!isAllRegular()) {
        bndrycent.ParallelCopy(m_bndrycent,0,0,bndrycent.nComp(),0,bndrycent.nGrow(),
                               geom.periodicity());
    }
}

void
Level::fillBndryNorm (MultiCutFab& bndrynorm, const Geometry& geom) const
{
    if (isAllRegular()) {
        bndrynorm.setVal(0.0);
        return;
    }

    MultiFab tmp(bndrynorm.boxArray(), bndrynorm.DistributionMap(),
                 bndrynorm.nComp(), bndrynorm.nGrow());
    fillBndryNorm(tmp, geom);
    copyMultiFabToMultiCutFab(bndrynorm, tmp);
}

void
Level::fillBndryNorm (   MultiFab& bndrynorm, const Geometry& geom) const
{
    bndrynorm.setVal(0.0);
    if (!isAllRegular()) {
        bndrynorm.ParallelCopy(m_bndrynorm,0,0,bndrynorm.nComp(),0,bndrynorm.nGrow(),
                               geom.periodicity());
    }
}
        
void
Level::fillAreaFrac (Array<MultiCutFab*,AMREX_SPACEDIM> const& a_areafrac, const Geometry& geom) const
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
        tmp.ParallelCopy(m_areafrac[idim],0,0,areafrac.nComp(),
                         0,areafrac.nGrow(),geom.periodicity());
        copyMultiFabToMultiCutFab(areafrac, tmp);
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
Level::fillAreaFrac (Array<MultiFab*,AMREX_SPACEDIM> const& a_areafrac, const Geometry& geom) const
{
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        a_areafrac[idim]->setVal(1.0);
    }

    if (isAllRegular()) return;

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        auto& areafrac = *a_areafrac[idim];
        areafrac.ParallelCopy(m_areafrac[idim],0,0,areafrac.nComp(),
                              0,areafrac.nGrow(),geom.periodicity());
    }

    const std::vector<IntVect>& pshifts = geom.periodicity().shiftIntVect();

    Real cov_val = 0.0; // for covered cells
        
#ifdef _OPENMP
#pragma omp parallel
#endif
    if (!m_covered_grids.empty())
    {
        std::vector<std::pair<int,Box> > isects;
        for (MFIter mfi(*a_areafrac[0]); mfi.isValid(); ++mfi)
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

void
Level::fillFaceCent (Array<MultiCutFab*,AMREX_SPACEDIM> const& a_facecent, const Geometry& geom) const
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
        tmp.ParallelCopy(m_facecent[idim],0,0,facecent.nComp(),
                         0,facecent.nGrow(),geom.periodicity());
        copyMultiFabToMultiCutFab(facecent,tmp);
    }
}

void
Level::fillFaceCent (Array<MultiFab*,AMREX_SPACEDIM> const& a_facecent, const Geometry& geom) const
{
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        a_facecent[idim]->setVal(0.0);
    }
    if (!isAllRegular()) {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            auto& facecent = *a_facecent[idim];
            a_facecent[idim]->ParallelCopy(m_facecent[idim],0,0,facecent.nComp(),
                                           0,facecent.nGrow(),geom.periodicity());
        }
    }
}

void
Level::fillLevelSet (MultiFab& levelset, const Geometry& geom) const
{
    levelset.setVal(-1.0);
    levelset.ParallelCopy(m_levelset,0,0,1,0,0);

    const std::vector<IntVect>& pshifts = geom.periodicity().shiftIntVect();

    Real cov_val = 1.0; // for covered cells

#ifdef _OPENMP
#pragma omp parallel
#endif
    if (!m_covered_grids.empty())
    {
        std::vector<std::pair<int,Box> > isects;
        for (MFIter mfi(levelset); mfi.isValid(); ++mfi)
        {
            FArrayBox& lsfab = levelset[mfi];
            const Box& ccbx = amrex::enclosedCells(lsfab.box());
            for (const auto& iv : pshifts)
            {
                m_covered_grids.intersections(ccbx+iv, isects);
                for (const auto& is : isects) {
                    const Box& fbx = amrex::surroundingNodes(is.second-iv);
                    lsfab.setVal(cov_val, fbx, 0, 1);
                }
            }
        }
    }
}
        
}}
