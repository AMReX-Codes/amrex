
#include <AMReX_EB2_Level.H>

namespace amrex { namespace EB2 {

void
Level::coarsenFromFine (Level& fineLevel)
{
    const BoxArray& fine_grids = fineLevel.m_grids;
    const BoxArray& fine_covered_grids = fineLevel.m_covered_grids;
    const DistributionMapping& fine_dmap = fineLevel.m_dmap;
    m_grids = amrex::coarsen(fine_grids,2);
    m_covered_grids = amrex::coarsen(fine_covered_grids, 2);
    m_dmap = fine_dmap;

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
                                m_dmap, AMREX_SPACEDIM, ng);
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

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(m_volfrac,true); mfi.isValid(); ++mfi)
    {
        // default
        {
            const Box& gbx = mfi.growntilebox(2);
            m_volfrac[mfi].setVal(1.0, gbx, 0, 1);
            m_centroid[mfi].setVal(0.0, gbx, 0, AMREX_SPACEDIM);
            m_bndryarea[mfi].setVal(0.0, gbx, 0, 1);
            m_bndrycent[mfi].setVal(0.0, gbx, 0, AMREX_SPACEDIM);
            m_bndrynorm[mfi].setVal(0.0, gbx, 0, AMREX_SPACEDIM);
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                const Box& fbx = mfi.grownnodaltilebox(idim,2);
                m_areafrac[idim][mfi].setVal(1.0, fbx, 0, 1);
                m_facecent[idim][mfi].setVal(0.0, fbx, 0, AMREX_SPACEDIM);
            }
        }

        const Box&  bx = mfi.growntilebox(1);
        const Box& xbx = mfi.grownnodaltilebox(0,1);
        const Box& ybx = mfi.grownnodaltilebox(1,1);
#if (AMREX_SPACEDIM == 3)
        const Box& zbx = mfi.grownnodaltilebox(2,1);
#endif

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
                                    BL_TO_FORTRAN_ANYD(f_cellflag[mfi]));
    }

    // xxxxx todo
    // multivalue and multicut detection
    // celflag neighbors

    m_valid = true;
}

void
Level::fillEBCellFlag (FabArray<EBCellFlagFab>& cellflag, const Geometry& geom) const
{
    cellflag.ParallelCopy(m_cellflag,0,0,1,0,cellflag.nGrow(),geom.periodicity());

    const std::vector<IntVect>& pshifts = geom.periodicity().shiftIntVect();

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
            for (const auto& iv : pshifts)
            {
                m_covered_grids.intersections(bx+iv, isects);
                for (const auto& is : isects) {
                    fab.setVal(cov_val, is.second-iv, 0, 1);
                }
            }

            // fix type and region for each fab
            fab.setRegion(bx);
            fab.setType(FabType::undefined);
            auto typ = fab.getType(bx);
            fab.setType(typ);
        }
    }
}

void
Level::fillVolFrac (MultiFab& vfrac, const Geometry& geom) const
{
    vfrac.setVal(1.0);
    vfrac.ParallelCopy(m_volfrac,0,0,1,0,vfrac.nGrow(),geom.periodicity());

    const std::vector<IntVect>& pshifts = geom.periodicity().shiftIntVect();

    Real cov_val = 0.0; // for covered cells
#ifdef _OPENMP
#pragma omp parallel
#endif
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

}}
