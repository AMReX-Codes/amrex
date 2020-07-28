
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

    if (Gpu::notInLaunchRegion())
    {
#ifdef _OPENMP
#pragma omp parallel reduction(max:mvmc_error)
#endif
        for (MFIter mfi(m_levelset,true); mfi.isValid(); ++mfi)
        {
            const Box& ccbx = mfi.tilebox(IntVect::TheCellVector());
            const Box& ndbx = mfi.tilebox();
            auto const& crse = m_levelset.array(mfi);
            auto const& fine = f_levelset.const_array(mfi);

            amrex::LoopConcurrentOnCpu(ndbx,
            [=] (int i, int j, int k) noexcept
            {
                crse(i,j,k) = fine(2*i,2*j,2*k);
            });

            int tile_error = 0;
            amrex::LoopOnCpu(ccbx,
            [=,&tile_error] (int i, int j, int k) noexcept
            {
                int ierror = check_mvmc(i,j,k,fine);
                tile_error = std::max(tile_error,ierror);
            });

            mvmc_error = std::max(mvmc_error, tile_error);
        }
    }
    else
    {
        ReduceOps<ReduceOpMax> reduce_op;
        ReduceData<int> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        for (MFIter mfi(m_levelset); mfi.isValid(); ++mfi)
        {
            const Box& ndbx = mfi.validbox();
            const Box& ccbx = amrex::enclosedCells(ndbx);
            auto const& crse = m_levelset.array(mfi);
            auto const& fine = f_levelset.const_array(mfi);
            reduce_op.eval(ndbx, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                crse(i,j,k) = fine(2*i,2*j,2*k);
                int ierror;
                if (ccbx.contains(IntVect{AMREX_D_DECL(i,j,k)})) {
                    ierror = check_mvmc(i,j,k,fine);
                } else {
                    ierror = 0;
                }
                return {ierror};
            });
        }
        ReduceTuple rv = reduce_data.value();
        mvmc_error = amrex::max(0, amrex::get<0>(rv));
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
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            {
                std::vector<std::pair<int,Box> > isects;
                for (MFIter mfi(f_volfrac); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.fabbox();
                    auto const& vfrac = f_volfrac.array(mfi);
                    auto const& cflag = f_cellflag.array(mfi);
                    AMREX_D_TERM(auto const& apx = f_areafrac[0].array(mfi);,
                                 auto const& apy = f_areafrac[1].array(mfi);,
                                 auto const& apz = f_areafrac[2].array(mfi););

                    for (const auto& iv : pshifts)
                    {
                        fine_covered_grids.intersections(bx+iv, isects);
                        for (const auto& is : isects)
                        {
                            Box const& ibox = is.second - iv;
                            Box const& indbox = amrex::surroundingNodes(ibox);
                            AMREX_D_TERM(auto const& xbx = amrex::surroundingNodes(ibox,0);,
                                         auto const& ybx = amrex::surroundingNodes(ibox,1);,
                                         auto const& zbx = amrex::surroundingNodes(ibox,2););
                            AMREX_HOST_DEVICE_FOR_3D(indbox, i,j,k,
                            {
                                IntVect cell(AMREX_D_DECL(i,j,k));
                                if (ibox.contains(cell)) {
                                    vfrac(i,j,k) = 0.0;
                                    cflag(i,j,k) = EBCellFlag::TheCoveredCell();
                                }
                                AMREX_D_TERM(if (xbx.contains(cell)) apx(i,j,k) = 0.0;,
                                             if (ybx.contains(cell)) apy(i,j,k) = 0.0;,
                                             if (zbx.contains(cell)) apz(i,j,k) = 0.0;);
                            });
                        }
                    }
                }
            }
        }
    }

    int error = 0;

    if (Gpu::notInLaunchRegion())
    {
#ifdef _OPENMP
#pragma omp parallel reduction(max:error)
#endif
        for (MFIter mfi(m_volfrac,true); mfi.isValid(); ++mfi)
        {
            auto const& cvol = m_volfrac.array(mfi);
            auto const& ccent = m_centroid.array(mfi);
            auto const& cba = m_bndryarea.array(mfi);
            auto const& cbc = m_bndrycent.array(mfi);
            auto const& cbn = m_bndrynorm.array(mfi);
            AMREX_D_TERM(auto const& capx = m_areafrac[0].array(mfi);,
                         auto const& capy = m_areafrac[1].array(mfi);,
                         auto const& capz = m_areafrac[2].array(mfi););
            AMREX_D_TERM(auto const& cfcx = m_facecent[0].array(mfi);,
                         auto const& cfcy = m_facecent[1].array(mfi);,
                         auto const& cfcz = m_facecent[2].array(mfi););
            auto const& cflag = m_cellflag.array(mfi);

            auto const& fvol = f_volfrac.const_array(mfi);
            auto const& fcent = f_centroid.const_array(mfi);
            auto const& fba = f_bndryarea.const_array(mfi);
            auto const& fbc = f_bndrycent.const_array(mfi);
            auto const& fbn = f_bndrynorm.const_array(mfi);
            AMREX_D_TERM(auto const& fapx = f_areafrac[0].const_array(mfi);,
                         auto const& fapy = f_areafrac[1].const_array(mfi);,
                         auto const& fapz = f_areafrac[2].const_array(mfi););
            AMREX_D_TERM(auto const& ffcx = f_facecent[0].const_array(mfi);,
                         auto const& ffcy = f_facecent[1].const_array(mfi);,
                         auto const& ffcz = f_facecent[2].const_array(mfi););
            auto const& fflag = f_cellflag.const_array(mfi);

            Box const& bx = mfi.validbox();
            AMREX_D_TERM(Box const& xbx = amrex::surroundingNodes(bx,0);,
                         Box const& ybx = amrex::surroundingNodes(bx,1);,
                         Box const& zbx = amrex::surroundingNodes(bx,2););
            Box const& gbx = amrex::grow(bx,2);
            AMREX_D_TERM(Box const& xgbx = amrex::surroundingNodes(gbx,0);,
                         Box const& ygbx = amrex::surroundingNodes(gbx,1);,
                         Box const& zgbx = amrex::surroundingNodes(gbx,2););

            Box const& ndgbx = mfi.grownnodaltilebox(-1,2);

            int tile_error = 0;
            amrex::LoopOnCpu(ndgbx,
            [=,&tile_error] (int i, int j, int k) noexcept
            {
                amrex::ignore_unused(j,k);
                int ierr = coarsen_from_fine(AMREX_D_DECL(i,j,k),
                                             bx, gbx,
                                             AMREX_D_DECL(xbx,ybx,zbx),
                                             AMREX_D_DECL(xgbx,ygbx,zgbx),
                                             cvol,ccent,cba,cbc,cbn,
                                             AMREX_D_DECL(capx,capy,capz),
                                             AMREX_D_DECL(cfcx,cfcy,cfcz),
                                             cflag,
                                             fvol,fcent,fba,fbc,fbn,
                                             AMREX_D_DECL(fapx,fapy,fapz),
                                             AMREX_D_DECL(ffcx,ffcy,ffcz),
                                             fflag);
                tile_error = std::max(tile_error,ierr);
            });

            error = std::max(error,tile_error);
        }
    }
    else
    {
        ReduceOps<ReduceOpMax> reduce_op;
        ReduceData<int> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        for (MFIter mfi(m_volfrac); mfi.isValid(); ++mfi)
        {
            auto const& cvol = m_volfrac.array(mfi);
            auto const& ccent = m_centroid.array(mfi);
            auto const& cba = m_bndryarea.array(mfi);
            auto const& cbc = m_bndrycent.array(mfi);
            auto const& cbn = m_bndrynorm.array(mfi);
            AMREX_D_TERM(auto const& capx = m_areafrac[0].array(mfi);,
                         auto const& capy = m_areafrac[1].array(mfi);,
                         auto const& capz = m_areafrac[2].array(mfi););
            AMREX_D_TERM(auto const& cfcx = m_facecent[0].array(mfi);,
                         auto const& cfcy = m_facecent[1].array(mfi);,
                         auto const& cfcz = m_facecent[2].array(mfi););
            auto const& cflag = m_cellflag.array(mfi);

            auto const& fvol = f_volfrac.const_array(mfi);
            auto const& fcent = f_centroid.const_array(mfi);
            auto const& fba = f_bndryarea.const_array(mfi);
            auto const& fbc = f_bndrycent.const_array(mfi);
            auto const& fbn = f_bndrynorm.const_array(mfi);
            AMREX_D_TERM(auto const& fapx = f_areafrac[0].const_array(mfi);,
                         auto const& fapy = f_areafrac[1].const_array(mfi);,
                         auto const& fapz = f_areafrac[2].const_array(mfi););
            AMREX_D_TERM(auto const& ffcx = f_facecent[0].const_array(mfi);,
                         auto const& ffcy = f_facecent[1].const_array(mfi);,
                         auto const& ffcz = f_facecent[2].const_array(mfi););
            auto const& fflag = f_cellflag.const_array(mfi);

            Box const& bx = mfi.validbox();
            AMREX_D_TERM(Box const& xbx = amrex::surroundingNodes(bx,0);,
                         Box const& ybx = amrex::surroundingNodes(bx,1);,
                         Box const& zbx = amrex::surroundingNodes(bx,2););
            Box const& gbx = amrex::grow(bx,2);
            AMREX_D_TERM(Box const& xgbx = amrex::surroundingNodes(gbx,0);,
                         Box const& ygbx = amrex::surroundingNodes(gbx,1);,
                         Box const& zgbx = amrex::surroundingNodes(gbx,2););
            Box const& ndgbx = amrex::surroundingNodes(gbx);

#ifdef AMREX_USE_DPCPP
            // xxxxx DPCPP todo: kernel size
            Vector<Array4<Real const> > htmp = {fvol, fcent, fba, fbc, fbn,
                                                AMREX_D_DECL(fapx,fapy,fapz),
                                                AMREX_D_DECL(ffcx,ffcy,ffcz)};
            Vector<Array4<Real> > htmp2 = {AMREX_D_DECL(capx,capy,capz)};
            Gpu::AsyncArray<Array4<Real const> > dtmp(htmp.data(), 5+2*AMREX_SPACEDIM);
            Gpu::AsyncArray<Array4<Real> > dtmp2(htmp2.data(), AMREX_SPACEDIM);
            auto dp = dtmp.data();
            auto dp2 = dtmp2.data();
#endif

            reduce_op.eval(ndgbx, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                amrex::ignore_unused(j,k);
                int ierr = coarsen_from_fine(AMREX_D_DECL(i,j,k),
                                             bx, gbx,
                                             AMREX_D_DECL(xbx,ybx,zbx),
                                             AMREX_D_DECL(xgbx,ygbx,zgbx),
                                             cvol,ccent,cba,cbc,cbn,
#ifdef AMREX_USE_DPCPP
                                             AMREX_D_DECL(dp2[0],dp2[1],dp2[2]),
#else
                                             AMREX_D_DECL(capx,capy,capz),
#endif
                                             AMREX_D_DECL(cfcx,cfcy,cfcz),
                                             cflag,
#ifdef AMREX_USE_DPCPP
                                             dp[0], dp[1], dp[2], dp[3], dp[4],
                                             dp[5], dp[6], dp[7], dp[8],
#if (AMREX_SPACEDIM == 3)
                                             dp[9], dp[10],
#endif
#else
                                             fvol,fcent,fba,fbc,fbn,
                                             AMREX_D_DECL(fapx,fapy,fapz),
                                             AMREX_D_DECL(ffcx,ffcy,ffcz),
#endif
                                             fflag);
                return {ierr};
            });
        }

        ReduceTuple rv = reduce_data.value();
        error = amrex::max(0, amrex::get<0>(rv));
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
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(m_cellflag,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        auto const& cflag = m_cellflag.array(mfi);
        AMREX_D_TERM(auto const& apx = m_areafrac[0].const_array(mfi);,
                     auto const& apy = m_areafrac[1].const_array(mfi);,
                     auto const& apz = m_areafrac[2].const_array(mfi););
        AMREX_HOST_DEVICE_FOR_3D ( bx, i, j, k,
        {
            build_cellflag_from_ap(AMREX_D_DECL(i,j,k),
                                   cflag, AMREX_D_DECL(apx,apy,apz));
        });
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
            fab.setType(FabType::regular);
        }
        return;
    }

    const int ng = cellflag.nGrow();

    cellflag.ParallelCopy(m_cellflag,0,0,1,0,ng,geom.periodicity());

    const std::vector<IntVect>& pshifts = geom.periodicity().shiftIntVect();

    auto cov_val = EBCellFlag::TheCoveredCell();
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        std::vector<std::pair<int,Box> > isects;
        for (MFIter mfi(cellflag,MFItInfo().UseDefaultStream()); mfi.isValid(); ++mfi)
        {
            auto& fab = cellflag[mfi];
            auto const& a = fab.array();
            const Box& bx = fab.box();
            if (!m_covered_grids.empty())
            {
                for (const auto& iv : pshifts)
                {
                    m_covered_grids.intersections(bx+iv, isects);
                    for (const auto& is : isects) {
                        Box const& ibox = is.second-iv;
                        AMREX_HOST_DEVICE_PARALLEL_FOR_3D(ibox, i, j, k,
                        {
                            a(i,j,k) = cov_val;
                        });
                    }
                }
            }

            // fix type for each fab
            fab.setType(FabType::undefined);
            auto typ = fab.getType(bx);
            fab.setType(typ);
            for (int nshrink = 1; nshrink < ng; ++nshrink) {
                const Box& b = amrex::grow(bx,-nshrink);
                fab.getType(b);
            }
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

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    if (!m_covered_grids.empty())
    {
        std::vector<std::pair<int,Box> > isects;
        for (MFIter mfi(vfrac); mfi.isValid(); ++mfi)
        {
            auto const& fab = vfrac.array(mfi);
            const Box& bx = mfi.fabbox();
            for (const auto& iv : pshifts)
            {
                m_covered_grids.intersections(bx+iv, isects);
                for (const auto& is : isects) {
                    Box const& ibox = is.second-iv;
                    AMREX_HOST_DEVICE_PARALLEL_FOR_3D(ibox, i, j, k,
                    {
                        fab(i,j,k) = 0.0;  // covered cells
                    });
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
                const auto srcfab = srcmf.const_array(mfi);
                const Box& box = mfi.fabbox();
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D (box,ncomp,i,j,k,n,
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

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    if (!m_covered_grids.empty())
    {
        std::vector<std::pair<int,Box> > isects;
        for (MFIter mfi(a_areafrac[0]->data()); mfi.isValid(); ++mfi)
        {
            if (a_areafrac[0]->ok(mfi))
            {
                const Box& ccbx = amrex::enclosedCells((*a_areafrac[0])[mfi].box());
                AMREX_D_TERM(auto const& apx = a_areafrac[0]->array(mfi);,
                             auto const& apy = a_areafrac[1]->array(mfi);,
                             auto const& apz = a_areafrac[2]->array(mfi););
                for (const auto& iv : pshifts)
                {
                    m_covered_grids.intersections(ccbx+iv, isects);
                    for (const auto& is : isects) {
                        if (Gpu::inLaunchRegion()) {
                            const Box& bx = is.second-iv;
                            AMREX_D_TERM(const Box& xbx = amrex::surroundingNodes(bx,0);,
                                         const Box& ybx = amrex::surroundingNodes(bx,1);,
                                         const Box& zbx = amrex::surroundingNodes(bx,2););
                            amrex::ParallelFor(AMREX_D_DECL(xbx,ybx,zbx),
                              [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                              {
                                  apx(i,j,k) = 0.0;
                              }
#if (AMREX_SPACEDIM >= 2)
                            , [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                              {
                                  apy(i,j,k) = 0.0;
                              }
#if (AMREX_SPACEDIM == 3)
                            , [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                              {
                                  apz(i,j,k) = 0.0;
                              }
#endif
#endif
                            );
                        } else {
                            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                                const Box& fbx = amrex::surroundingNodes(is.second-iv,idim);
                                (*a_areafrac[idim])[mfi].setVal<RunOn::Host>(0.0, fbx, 0, 1);
                            }
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

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    if (!m_covered_grids.empty())
    {
        std::vector<std::pair<int,Box> > isects;
        for (MFIter mfi(*a_areafrac[0]); mfi.isValid(); ++mfi)
        {
            const Box& ccbx = amrex::enclosedCells((*a_areafrac[0])[mfi].box());
            AMREX_D_TERM(auto const& apx = a_areafrac[0]->array(mfi);,
                         auto const& apy = a_areafrac[1]->array(mfi);,
                         auto const& apz = a_areafrac[2]->array(mfi););
            for (const auto& iv : pshifts)
            {
                m_covered_grids.intersections(ccbx+iv, isects);
                for (const auto& is : isects) {
                    if (Gpu::inLaunchRegion()) {
                        const Box& bx = is.second-iv;
                        AMREX_D_TERM(const Box& xbx = amrex::surroundingNodes(bx,0);,
                                     const Box& ybx = amrex::surroundingNodes(bx,1);,
                                     const Box& zbx = amrex::surroundingNodes(bx,2););
                        amrex::ParallelFor(AMREX_D_DECL(xbx,ybx,zbx),
                          [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                          {
                              apx(i,j,k) = 0.0;
                          }
#if (AMREX_SPACEDIM >= 2)
                        , [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                          {
                              apy(i,j,k) = 0.0;
                          }
#if (AMREX_SPACEDIM == 3)
                        , [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                          {
                              apz(i,j,k) = 0.0;
                          }
#endif
#endif
                        );
                    } else {
                        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                            const Box& fbx = amrex::surroundingNodes(is.second-iv,idim);
                            (*a_areafrac[idim])[mfi].setVal<RunOn::Host>(0.0, fbx, 0, 1);
                        }
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
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    if (!m_covered_grids.empty())
    {
        std::vector<std::pair<int,Box> > isects;
        for (MFIter mfi(levelset); mfi.isValid(); ++mfi)
        {
            const auto& lsfab = levelset.array(mfi);
            const Box& ccbx = amrex::enclosedCells(mfi.fabbox());
            for (const auto& iv : pshifts)
            {
                m_covered_grids.intersections(ccbx+iv, isects);
                for (const auto& is : isects) {
                    const Box& fbx = amrex::surroundingNodes(is.second-iv);
                    AMREX_HOST_DEVICE_PARALLEL_FOR_3D(fbx, i, j, k,
                    {
                        lsfab(i,j,k) = cov_val;
                    });
                }
            }
        }
    }
}
        
}}
