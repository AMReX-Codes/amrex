
#include <AMReX_MultiFabUtil.H>
#include <sstream>
#include <iostream>

namespace {

    using namespace amrex;

    Box
    getIndexBox(const RealBox& real_box, const Geometry& geom) {
        IntVect slice_lo, slice_hi;

        AMREX_D_TERM(slice_lo[0]=static_cast<int>(std::floor((real_box.lo(0) - geom.ProbLo(0))/geom.CellSize(0)));,
                     slice_lo[1]=static_cast<int>(std::floor((real_box.lo(1) - geom.ProbLo(1))/geom.CellSize(1)));,
                     slice_lo[2]=static_cast<int>(std::floor((real_box.lo(2) - geom.ProbLo(2))/geom.CellSize(2))););

        AMREX_D_TERM(slice_hi[0]=static_cast<int>(std::floor((real_box.hi(0) - geom.ProbLo(0))/geom.CellSize(0)));,
                     slice_hi[1]=static_cast<int>(std::floor((real_box.hi(1) - geom.ProbLo(1))/geom.CellSize(1)));,
                     slice_hi[2]=static_cast<int>(std::floor((real_box.hi(2) - geom.ProbLo(2))/geom.CellSize(2))););

        return Box(slice_lo, slice_hi) & geom.Domain();
    }


    std::unique_ptr<MultiFab> allocateSlice(int dir, const MultiFab& cell_centered_data,
                                            int ncomp, const Geometry& geom, Real dir_coord,
                                            Vector<int>& slice_to_full_ba_map) {

        // Get our slice and convert to index space
        RealBox real_slice = geom.ProbDomain();
        real_slice.setLo(dir, dir_coord);
        real_slice.setHi(dir, dir_coord);
        Box slice_box = getIndexBox(real_slice, geom);

        // define the multifab that stores slice
        BoxArray const& ba = cell_centered_data.boxArray();
        const DistributionMapping& dm = cell_centered_data.DistributionMap();
        std::vector< std::pair<int, Box> > isects;
        ba.intersections(slice_box, isects, false, 0);
        Vector<Box> boxes;
        Vector<int> procs;
        for (auto const& is : isects) {
            procs.push_back(dm[is.first]);
            boxes.push_back(is.second);
            slice_to_full_ba_map.push_back(is.first);
        }
        BoxArray slice_ba(&boxes[0], static_cast<int>(boxes.size()));
        DistributionMapping slice_dmap(std::move(procs));
        std::unique_ptr<MultiFab> slice(new MultiFab(slice_ba, slice_dmap, ncomp, 0,
                                                     MFInfo(), cell_centered_data.Factory()));
        return slice;
    }
}

namespace amrex
{
    void average_node_to_cellcenter (MultiFab& cc, int dcomp,
         const MultiFab& nd, int scomp, int ncomp, int ngrow)
    {
#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion() && cc.isFusingCandidate()) {
            auto const& ccma = cc.arrays();
            auto const& ndma = nd.const_arrays();
            ParallelFor(cc, IntVect(ngrow), ncomp,
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
            {
                amrex_avg_nd_to_cc(i, j, k, n, ccma[box_no], ndma[box_no], dcomp, scomp);
            });
            if (!Gpu::inNoSyncRegion()) {
                Gpu::streamSynchronize();
            }
        } else
#endif
        {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(cc,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box bx = mfi.growntilebox(ngrow);
                Array4<Real> const& ccarr = cc.array(mfi);
                Array4<Real const> const& ndarr = nd.const_array(mfi);

                AMREX_HOST_DEVICE_PARALLEL_FOR_4D( bx, ncomp, i, j, k, n,
                {
                    amrex_avg_nd_to_cc(i, j, k, n, ccarr, ndarr, dcomp, scomp); // NOLINT(readability-suspicious-call-argument)
                });
            }
        }
    }

    void average_edge_to_cellcenter (MultiFab& cc, int dcomp,
        const Vector<const MultiFab*>& edge, int ngrow)
    {
        AMREX_ASSERT(cc.nComp() >= dcomp + AMREX_SPACEDIM);
        AMREX_ASSERT(edge.size() == AMREX_SPACEDIM);
        AMREX_ASSERT(edge[0]->nComp() == 1);
#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion() && cc.isFusingCandidate()) {
            auto const& ccma = cc.arrays();
            AMREX_D_TERM(auto const& exma = edge[0]->const_arrays();,
                         auto const& eyma = edge[1]->const_arrays();,
                         auto const& ezma = edge[2]->const_arrays(););
            ParallelFor(cc, IntVect(ngrow),
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
            {
                amrex_avg_eg_to_cc(i, j, k, ccma[box_no],
                                   AMREX_D_DECL(exma[box_no], eyma[box_no], ezma[box_no]),
                                   dcomp);
            });
            if (!Gpu::inNoSyncRegion()) {
                Gpu::streamSynchronize();
            }
        } else
#endif
        {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(cc,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box bx = mfi.growntilebox(ngrow);
                Array4<Real> const& ccarr = cc.array(mfi);
                AMREX_D_TERM(Array4<Real const> const& exarr = edge[0]->const_array(mfi);,
                             Array4<Real const> const& eyarr = edge[1]->const_array(mfi);,
                             Array4<Real const> const& ezarr = edge[2]->const_array(mfi););

                AMREX_HOST_DEVICE_PARALLEL_FOR_3D( bx, i, j, k,
                {
                    amrex_avg_eg_to_cc(i, j, k, ccarr, AMREX_D_DECL(exarr,eyarr,ezarr), dcomp);
                });
            }
        }
    }

    void average_face_to_cellcenter (MultiFab& cc, int dcomp,
        const Vector<const MultiFab*>& fc, int ngrow)
    {
        average_face_to_cellcenter(cc, dcomp,
            Array<MultiFab const*,AMREX_SPACEDIM>{{AMREX_D_DECL(fc[0],fc[1],fc[2])}},
            ngrow);
    }

    void average_face_to_cellcenter (MultiFab& cc, const Vector<const MultiFab*>& fc,
                                     const Geometry& geom)
    {
        average_face_to_cellcenter(cc,
                                   Array<MultiFab const*,AMREX_SPACEDIM>
                                   {{AMREX_D_DECL(fc[0],fc[1],fc[2])}},
                                   geom);
    }

    void average_face_to_cellcenter (MultiFab& cc,
                                     const Array<const MultiFab*,AMREX_SPACEDIM>& fc,
                                     const Geometry& geom)
    {
        AMREX_ASSERT(cc.nComp() >= AMREX_SPACEDIM);
        AMREX_ASSERT(fc[0]->nComp() == 1); // We only expect fc to have the gradient perpendicular to the face

        const GeometryData gd = geom.data();
        amrex::ignore_unused(gd);

#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion() && cc.isFusingCandidate()) {
            auto const& ccma = cc.arrays();
            AMREX_D_TERM(auto const& fxma = fc[0]->const_arrays();,
                         auto const& fyma = fc[1]->const_arrays();,
                         auto const& fzma = fc[2]->const_arrays(););
            ParallelFor(cc,
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
            {
                amrex_avg_fc_to_cc(i,j,k, ccma[box_no], AMREX_D_DECL(fxma[box_no],
                                                                     fyma[box_no],
                                                                     fzma[box_no]),
                                   0
#if (AMREX_SPACEDIM == 1)
                                   , gd
#endif
                                   );
            });
            if (!Gpu::inNoSyncRegion()) {
                Gpu::streamSynchronize();
            }
        } else
#endif
        {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(cc,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box bx = mfi.tilebox();
                Array4<Real> const& ccarr = cc.array(mfi);
                AMREX_D_TERM(Array4<Real const> const& fxarr = fc[0]->const_array(mfi);,
                             Array4<Real const> const& fyarr = fc[1]->const_array(mfi);,
                             Array4<Real const> const& fzarr = fc[2]->const_array(mfi););

#if (AMREX_SPACEDIM == 1)
                AMREX_HOST_DEVICE_PARALLEL_FOR_3D( bx, i, j, k,
                {
                    amrex_avg_fc_to_cc(i,j,k, ccarr, fxarr, 0, gd);
                });
#else
                AMREX_HOST_DEVICE_PARALLEL_FOR_3D( bx, i, j, k,
                {
                    amrex_avg_fc_to_cc(i,j,k, ccarr, AMREX_D_DECL(fxarr,fyarr,fzarr), 0);
                });
#endif
            }
        }
    }

    void average_cellcenter_to_face (const Vector<MultiFab*>& fc, const MultiFab& cc,
                                     const Geometry& geom, int ncomp, bool use_harmonic_averaging)
    {
        average_cellcenter_to_face(Array<MultiFab*,AMREX_SPACEDIM>{{AMREX_D_DECL(fc[0],fc[1],fc[2])}},
                                   cc, geom, ncomp, use_harmonic_averaging);
    }


    void average_cellcenter_to_face (const Array<MultiFab*,AMREX_SPACEDIM>& fc, const MultiFab& cc,
                                    const Geometry& geom, int ncomp, bool use_harmonic_averaging)
    {
        AMREX_ASSERT(cc.nComp() == ncomp);
        AMREX_ASSERT(cc.nGrow() >= 1);
        AMREX_ASSERT(fc[0]->nComp() == ncomp); // We only expect fc to have the gradient perpendicular to the face
#if (AMREX_SPACEDIM >= 2)
        AMREX_ASSERT(fc[1]->nComp() == ncomp); // We only expect fc to have the gradient perpendicular to the face
#endif
#if (AMREX_SPACEDIM == 3)
        AMREX_ASSERT(fc[2]->nComp() == ncomp); // We only expect fc to have the gradient perpendicular to the face
#endif


#if (AMREX_SPACEDIM == 1)
        const GeometryData& gd = geom.data();
        if (use_harmonic_averaging) {
            AMREX_ASSERT(gd.Coord() == 0);
        }
#else
        amrex::ignore_unused(geom);
#endif

#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion() && cc.isFusingCandidate()) {
            auto const& ccma = cc.const_arrays();
            AMREX_D_TERM(auto const& fxma = fc[0]->arrays();,
                         auto const& fyma = fc[1]->arrays();,
                         auto const& fzma = fc[2]->arrays(););
            MultiFab foo(amrex::convert(cc.boxArray(),IntVect(1)), cc.DistributionMap(), 1, 0,
                         MFInfo().SetAlloc(false));
            IntVect ng = -cc.nGrowVect();
            ParallelFor(foo, IntVect(0), ncomp,
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
            {
                Box ccbx(ccma[box_no]);
                ccbx.grow(ng);
                AMREX_D_TERM(Box const& xbx = amrex::surroundingNodes(ccbx,0);,
                             Box const& ybx = amrex::surroundingNodes(ccbx,1);,
                             Box const& zbx = amrex::surroundingNodes(ccbx,2););
#if (AMREX_SPACEDIM == 1)
                amrex_avg_cc_to_fc(i,j,k,n, xbx, fxma[box_no], ccma[box_no], gd, use_harmonic_averaging);
#else
                amrex_avg_cc_to_fc(i,j,k,n, AMREX_D_DECL(xbx,ybx,zbx),
                                   AMREX_D_DECL(fxma[box_no],fyma[box_no],fzma[box_no]),
                                   ccma[box_no], use_harmonic_averaging);
#endif
            });
            if (!Gpu::inNoSyncRegion()) {
                Gpu::streamSynchronize();
            }
        } else
#endif
        {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(cc,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                AMREX_D_TERM(const Box& xbx = mfi.nodaltilebox(0);,
                             const Box& ybx = mfi.nodaltilebox(1);,
                             const Box& zbx = mfi.nodaltilebox(2););
                const auto& index_bounds = amrex::getIndexBounds(AMREX_D_DECL(xbx,ybx,zbx));

                AMREX_D_TERM(Array4<Real> const& fxarr = fc[0]->array(mfi);,
                             Array4<Real> const& fyarr = fc[1]->array(mfi);,
                             Array4<Real> const& fzarr = fc[2]->array(mfi););
                Array4<Real const> const& ccarr = cc.const_array(mfi);

#if (AMREX_SPACEDIM == 1)
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D(index_bounds, ncomp, i, j, k, n,
                {
                    amrex_avg_cc_to_fc(i,j,k,n, xbx, fxarr, ccarr, gd, use_harmonic_averaging);
                });
#else
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D(index_bounds, ncomp, i, j, k, n,
                {
                    amrex_avg_cc_to_fc(i,j,k,n, AMREX_D_DECL(xbx,ybx,zbx),
                                       AMREX_D_DECL(fxarr,fyarr,fzarr), ccarr,
                                       use_harmonic_averaging);
                });
#endif
            }
        }
    }

// *************************************************************************************************************

    // Average fine cell-based MultiFab onto crse cell-centered MultiFab.
    // We do NOT assume that the coarse layout is a coarsened version of the fine layout.
    // This version DOES use volume-weighting.
    void average_down (const MultiFab& S_fine, MultiFab& S_crse,
                       const Geometry& fgeom, const Geometry& cgeom,
                       int scomp, int ncomp, int rr)
    {
         average_down(S_fine,S_crse,fgeom,cgeom,scomp,ncomp,rr*IntVect::TheUnitVector());
    }

    void average_down (const MultiFab& S_fine, MultiFab& S_crse,
                       const Geometry& fgeom, const Geometry& cgeom,
                       int scomp, int ncomp, const IntVect& ratio)
    {
        amrex::ignore_unused(fgeom,cgeom);

        BL_PROFILE("amrex::average_down_w_geom");

        if (S_fine.is_nodal() || S_crse.is_nodal())
        {
            amrex::Error("Can't use amrex::average_down for nodal MultiFab!");
        }

#if (AMREX_SPACEDIM == 3)
        amrex::average_down(S_fine, S_crse, scomp, ncomp, ratio);
#else

        AMREX_ASSERT(S_crse.nComp() == S_fine.nComp());

        //
        // Coarsen() the fine stuff on processors owning the fine data.
        //
        const BoxArray& fine_BA = S_fine.boxArray();
        const DistributionMapping& fine_dm = S_fine.DistributionMap();
        BoxArray crse_S_fine_BA = fine_BA;
        crse_S_fine_BA.coarsen(ratio);

        MultiFab crse_S_fine(crse_S_fine_BA,fine_dm,ncomp,0,MFInfo(),FArrayBoxFactory());

        MultiFab fvolume;
        fgeom.GetVolume(fvolume, fine_BA, fine_dm, 0);

#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion() && crse_S_fine.isFusingCandidate()) {
            auto const& crsema = crse_S_fine.arrays();
            auto const& finema = S_fine.const_arrays();
            auto const& finevolma = fvolume.const_arrays();
            ParallelFor(crse_S_fine, IntVect(0), ncomp,
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
            {
                amrex_avgdown_with_vol(i,j,k,n,crsema[box_no],finema[box_no],finevolma[box_no],
                                       0,scomp,ratio);
            });
            if (!Gpu::inNoSyncRegion()) {
                Gpu::streamSynchronize();
            }
        } else
#endif
        {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(crse_S_fine,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                //  NOTE: The tilebox is defined at the coarse level.
                const Box& bx = mfi.tilebox();
                Array4<Real> const& crsearr = crse_S_fine.array(mfi);
                Array4<Real const> const& finearr = S_fine.const_array(mfi);
                Array4<Real const> const& finevolarr = fvolume.const_array(mfi);

                AMREX_HOST_DEVICE_PARALLEL_FOR_4D(bx, ncomp, i, j, k, n,
                {
                    amrex_avgdown_with_vol(i,j,k,n,crsearr,finearr,finevolarr,0,scomp,ratio);
                });
            }
        }

        S_crse.ParallelCopy(crse_S_fine,0,scomp,ncomp);
#endif
   }

    void sum_fine_to_coarse(const MultiFab& S_fine, MultiFab& S_crse,
                            int scomp, int ncomp, const IntVect& ratio,
                            const Geometry& cgeom, const Geometry& /*fgeom*/)
    {
        AMREX_ASSERT(S_crse.nComp() == S_fine.nComp());
        AMREX_ASSERT(ratio == ratio[0]);
        AMREX_ASSERT(S_fine.nGrow() % ratio[0] == 0);

        const int nGrow = S_fine.nGrow() / ratio[0];

        //
        // Coarsen() the fine stuff on processors owning the fine data.
        //
        BoxArray crse_S_fine_BA = S_fine.boxArray(); crse_S_fine_BA.coarsen(ratio);

        MultiFab crse_S_fine(crse_S_fine_BA, S_fine.DistributionMap(), ncomp, nGrow, MFInfo(), FArrayBoxFactory());

#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion() && crse_S_fine.isFusingCandidate()) {
            auto const& crsema = crse_S_fine.arrays();
            auto const& finema = S_fine.const_arrays();
            ParallelFor(crse_S_fine, IntVect(nGrow), ncomp,
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
            {
                amrex_avgdown(i,j,k,n,crsema[box_no],finema[box_no],0,scomp,ratio);
            });
            if (!Gpu::inNoSyncRegion()) {
                Gpu::streamSynchronize();
            }
        } else
#endif
        {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(crse_S_fine, TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                //  NOTE: The tilebox is defined at the coarse level.
                const Box& bx = mfi.growntilebox(nGrow);
                Array4<Real> const& crsearr = crse_S_fine.array(mfi);
                Array4<Real const> const& finearr = S_fine.const_array(mfi);

                AMREX_HOST_DEVICE_PARALLEL_FOR_4D(bx, ncomp, i, j, k, n,
                {
                    amrex_avgdown(i,j,k,n,crsearr,finearr,0,scomp,ratio);
                });
            }
        }

        S_crse.ParallelCopy(crse_S_fine, 0, scomp, ncomp, nGrow, 0,
                            cgeom.periodicity(), FabArrayBase::ADD);
    }

    //! Average fine edge-based MultiFab onto crse edge-based MultiFab.
    //! This routine assumes that the crse BoxArray is a coarsened version of the fine BoxArray.
    void average_down_edges (const Vector<const MultiFab*>& fine, const Vector<MultiFab*>& crse,
                             const IntVect& ratio, int ngcrse)
    {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            average_down_edges(*fine[idim], *crse[idim], ratio, ngcrse);
        }
    }

    void average_down_edges (const Array<const MultiFab*,AMREX_SPACEDIM>& fine,
                             const Array<MultiFab*,AMREX_SPACEDIM>& crse,
                             const IntVect& ratio, int ngcrse)
    {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            average_down_edges(*fine[idim], *crse[idim], ratio, ngcrse);
        }
    }

    void average_down_edges (const MultiFab& fine, MultiFab& crse,
                             const IntVect& ratio, int ngcrse)
    {
        AMREX_ASSERT(crse.nComp() == fine.nComp());
        AMREX_ASSERT(fine.ixType() == crse.ixType());
        const auto type = fine.ixType();
        int dir;
        for (dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            if (type.cellCentered(dir)) break;
        }
        auto tmptype = type;
        tmptype.set(dir);
        if (dir >= AMREX_SPACEDIM || !tmptype.nodeCentered()) {
            amrex::Abort("average_down_edges: not face index type");
        }
        const int ncomp = crse.nComp();
        if (isMFIterSafe(fine, crse))
        {
#ifdef AMREX_USE_GPU
            if (Gpu::inLaunchRegion() && crse.isFusingCandidate()) {
                auto const& crsema = crse.arrays();
                auto const& finema = fine.const_arrays();
                ParallelFor(crse, IntVect(ngcrse), ncomp,
                [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
                {
                    amrex_avgdown_edges(i,j,k,n, crsema[box_no], finema[box_no], 0, 0, ratio, dir);
                });
                if (!Gpu::inNoSyncRegion()) {
                    Gpu::streamSynchronize();
                }
            } else
#endif
            {
#ifdef AMREX_USE_OMP
#pragma omp parallel if(Gpu::notInLaunchRegion())
#endif
                for (MFIter mfi(crse,TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.growntilebox(ngcrse);
                    Array4<Real> const& crsearr = crse.array(mfi);
                    Array4<Real const> const& finearr = fine.const_array(mfi);

                    AMREX_HOST_DEVICE_PARALLEL_FOR_4D(bx, ncomp, i, j, k, n,
                    {
                        amrex_avgdown_edges(i,j,k,n, crsearr, finearr, 0, 0, ratio, dir);
                    });
                }
            }
        }
        else
        {
            MultiFab ctmp(amrex::coarsen(fine.boxArray(),ratio), fine.DistributionMap(),
                          ncomp, ngcrse, MFInfo(), FArrayBoxFactory());
            average_down_edges(fine, ctmp, ratio, ngcrse);
            crse.ParallelCopy(ctmp,0,0,ncomp,ngcrse,ngcrse);
        }
    }

    void print_state(const MultiFab& mf, const IntVect& cell, const int n, const IntVect& ng)
    {
        printCell(mf, cell, n, ng);
    }

    void writeFabs (const MultiFab& mf, const std::string& name)
    {
        writeFabs (mf, 0, mf.nComp(), name);
    }

    void writeFabs (const MultiFab& mf, int comp, int ncomp, const std::string& name)
    {
        for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
            std::ofstream ofs(name+"-fab-"+std::to_string(mfi.index()));
            mf[mfi].writeOn(ofs, comp, ncomp);
        }
    }

    MultiFab ToMultiFab (const iMultiFab& imf)
    {
        return amrex::cast<MultiFab>(imf);
    }

    FabArray<BaseFab<Long> > ToLongMultiFab (const iMultiFab& imf)
    {
        return amrex::cast<FabArray<BaseFab<Long> > > (imf);
    }

    std::unique_ptr<MultiFab> get_slice_data(int dir, Real coord, const MultiFab& cc, const Geometry& geom, int start_comp, int ncomp, bool interpolate) {

        BL_PROFILE("amrex::get_slice_data");

        if (interpolate) {
            AMREX_ASSERT(cc.nGrow() >= 1);
        }

        const auto geomdata = geom.data();

        Vector<int> slice_to_full_ba_map;
        std::unique_ptr<MultiFab> slice = allocateSlice(dir, cc, ncomp, geom, coord, slice_to_full_ba_map);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*slice, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            int slice_gid = mfi.index();
            int full_gid = slice_to_full_ba_map[slice_gid];
            auto& slice_fab = (*slice)[mfi];
            auto const& full_fab = cc[full_gid];
            Array4<Real> const& slice_arr = slice_fab.array();
            Array4<Real const> const& full_arr = full_fab.const_array();

            const Box& tile_box  = mfi.tilebox();

            if (interpolate)
            {
                AMREX_LAUNCH_HOST_DEVICE_LAMBDA( tile_box, thread_box,
                {
                    amrex_fill_slice_interp(thread_box, slice_arr, full_arr,
                                            0, start_comp, ncomp,
                                            dir, coord, geomdata);
                });
            }
            else
            {
                slice_fab.copy<RunOn::Device>(full_fab, tile_box, start_comp, tile_box, 0, ncomp);
            }
        }

        return slice;
    }

    iMultiFab makeFineMask (const BoxArray& cba, const DistributionMapping& cdm,
                            const BoxArray& fba, const IntVect& ratio,
                            int crse_value, int fine_value)
    {
        return makeFineMask(cba, cdm, IntVect{0}, fba, ratio, Periodicity::NonPeriodic(),
                            crse_value, fine_value);
    }

    template <typename FAB>
    void makeFineMask_doit (FabArray<FAB>& mask, const BoxArray& fba,
                            const IntVect& ratio, Periodicity const& period,
                            typename FAB::value_type crse_value,
                            typename FAB::value_type fine_value)
    {
        using value_type = typename FAB::value_type;

        Vector<Array4BoxTag<value_type> > tags;

        bool run_on_gpu = Gpu::inLaunchRegion();

        const BoxArray& cfba = amrex::coarsen(fba,ratio);
        const std::vector<IntVect>& pshifts = period.shiftIntVect();
#ifdef AMREX_USE_OMP
#pragma omp parallel if (!run_on_gpu)
#endif
        {
            std::vector <std::pair<int,Box> > isects;
            for (MFIter mfi(mask); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.fabbox();
                Array4<value_type> const& arr = mask.array(mfi);
                auto& fab = mask[mfi];

                AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
                {
                    arr(i,j,k) = crse_value;
                });

                for (const auto& iv : pshifts) {
                    cfba.intersections(bx+iv, isects);
                    for (const auto& is : isects) {
                        Box const& b = is.second-iv;
                        if (run_on_gpu) {
                            tags.push_back({arr,b});
                        } else {
                            fab.template setVal<RunOn::Host>(fine_value, b);
                        }
                    }
                }
            }
        }

#ifdef AMREX_USE_GPU
        amrex::ParallelFor(tags, 1,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n,
                              Array4BoxTag<value_type> const& tag) noexcept
        {
            tag.dfab(i,j,k,n) = fine_value;
        });
#endif
    }

    iMultiFab makeFineMask (const BoxArray& cba, const DistributionMapping& cdm,
                            const IntVect& cnghost, const BoxArray& fba, const IntVect& ratio,
                            Periodicity const& period, int crse_value, int fine_value)
    {
        iMultiFab mask(cba, cdm, 1, cnghost);
        makeFineMask_doit(mask, fba, ratio, period, crse_value, fine_value);
        return mask;
    }

    MultiFab makeFineMask (const BoxArray& cba, const DistributionMapping& cdm,
                           const BoxArray& fba, const IntVect& ratio,
                           Real crse_value, Real fine_value)
    {
        MultiFab mask(cba, cdm, 1, 0);
        makeFineMask_doit(mask, fba, ratio, Periodicity::NonPeriodic(), crse_value, fine_value);
        return mask;
    }

    void computeDivergence (MultiFab& divu, const Array<MultiFab const*,AMREX_SPACEDIM>& umac,
                            const Geometry& geom)
    {
        AMREX_D_TERM(AMREX_ASSERT(divu.nComp()==umac[0]->nComp());,
                     AMREX_ASSERT(divu.nComp()==umac[1]->nComp());,
                     AMREX_ASSERT(divu.nComp()==umac[2]->nComp()));

#if (AMREX_SPACEDIM==2)
        const auto& ba = divu.boxArray();
        const auto& dm = divu.DistributionMap();
        MultiFab volume, areax, areay;
        if (geom.IsRZ()) {
            geom.GetVolume(volume, ba, dm, 0);
            geom.GetFaceArea(areax, ba, dm, 0, 0);
            geom.GetFaceArea(areay, ba, dm, 1, 0);
        }
#endif

        const GpuArray<Real,AMREX_SPACEDIM> dxinv = geom.InvCellSizeArray();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(divu,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Array4<Real> const& divuarr = divu.array(mfi);
            AMREX_D_TERM(Array4<Real const> const& uarr = umac[0]->const_array(mfi);,
                         Array4<Real const> const& varr = umac[1]->const_array(mfi);,
                         Array4<Real const> const& warr = umac[2]->const_array(mfi););

#if (AMREX_SPACEDIM==2)
            if (geom.IsRZ()) {
                Array4<Real const> const&  ax =  areax.array(mfi);
                Array4<Real const> const&  ay =  areay.array(mfi);
                Array4<Real const> const& vol = volume.array(mfi);

                AMREX_LAUNCH_HOST_DEVICE_LAMBDA (bx, tbx,
                {
                    amrex_compute_divergence_rz(tbx,divuarr,AMREX_D_DECL(uarr,varr,warr),ax,ay,vol);
                });
            } else
#endif
            {
                AMREX_LAUNCH_HOST_DEVICE_LAMBDA (bx, tbx,
                {
                    amrex_compute_divergence(tbx,divuarr,AMREX_D_DECL(uarr,varr,warr),dxinv);
                });
            }
        }
    }

    void computeGradient (MultiFab& grad,  const Array<MultiFab const*,AMREX_SPACEDIM>& umac,
                          const Geometry& geom)
    {
        AMREX_ASSERT(grad.nComp() >= AMREX_SPACEDIM);

#if (AMREX_SPACEDIM==2)
        const auto& ba = grad.boxArray();
        const auto& dm = grad.DistributionMap();
        MultiFab volume, areax, areay;
        if (geom.IsRZ()) {
            geom.GetVolume(volume, ba, dm, 0);
            geom.GetFaceArea(areax, ba, dm, 0, 0);
            geom.GetFaceArea(areay, ba, dm, 1, 0);
        }
#endif

        const GpuArray<Real,AMREX_SPACEDIM> dxinv = geom.InvCellSizeArray();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(grad,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            const auto& gradfab = grad.array(mfi);
            AMREX_D_TERM(const auto& ufab = umac[0]->const_array(mfi);,
                         const auto& vfab = umac[1]->const_array(mfi);,
                         const auto& wfab = umac[2]->const_array(mfi););
#if (AMREX_SPACEDIM==2)
            if (geom.IsRZ()) {
                Array4<Real const> const&  ax =  areax.array(mfi);
                Array4<Real const> const&  ay =  areay.array(mfi);
                Array4<Real const> const& vol = volume.array(mfi);

                AMREX_LAUNCH_HOST_DEVICE_LAMBDA (bx, tbx,
                {
                    amrex_compute_gradient_rz(tbx,gradfab,AMREX_D_DECL(ufab,vfab,wfab),ax,ay,vol);
                });
            } else
#endif
            {
                AMREX_LAUNCH_HOST_DEVICE_LAMBDA (bx, tbx,
                {
                    amrex_compute_gradient(tbx,gradfab,AMREX_D_DECL(ufab,vfab,wfab),dxinv);
                });
            }
        }
    }

    MultiFab periodicShift (MultiFab const& mf, IntVect const& offset,
                            Periodicity const& period)
    {
        MultiFab r(mf.boxArray(), mf.DistributionMap(), mf.nComp(), 0);

        BoxList bl = mf.boxArray().boxList();
        for (auto& b : bl) {
            b.shift(offset);
        }
        BoxArray nba(std::move(bl));
        MultiFab nmf(nba, mf.DistributionMap(), mf.nComp(), 0, MFInfo().SetAlloc(false));

        for (MFIter mfi(r); mfi.isValid(); ++mfi) {
            auto const& rfab = r[mfi];
            nmf.setFab(mfi, FArrayBox(amrex::shift(rfab.box(),offset), rfab.nComp(),
                                      rfab.dataPtr()));
        }

        nmf.ParallelCopy(mf, period);

        return r;
    }


    Gpu::HostVector<Real> sumToLine (MultiFab const& mf, int icomp, int ncomp,
                                     Box const& domain, int direction, bool local)
    {
        int n1d = domain.length(direction) * ncomp;
        Gpu::HostVector<Real> hv(n1d);

#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion())
        {
            Gpu::DeviceVector<Real> dv(domain.length(direction), Real(0.0));
            Real* p = dv.data();

            for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
                Box const& b = mfi.validbox();
                const auto lo = amrex::lbound(b);
                const auto len = amrex::length(b);
                auto const& fab = mf.const_array(mfi);

                int n2d, n2dx;
                if (direction == 0) {
                    n2d = len.y * len.z;
                    n2dx = len.y;
                } else if (direction == 1) {
                    n2d = len.x * len.z;
                    n2dx = len.x;
                } else {
                    n2d = len.x * len.y;
                    n2dx = len.x;
                }
                int n2dblocks = (n2d+AMREX_GPU_MAX_THREADS-1)/AMREX_GPU_MAX_THREADS;
                int nblocks = n2dblocks * b.length(direction);
#ifdef AMREX_USE_SYCL
                std::size_t shared_mem_byte = sizeof(Real)*Gpu::Device::warp_size;
                amrex::launch(nblocks, AMREX_GPU_MAX_THREADS, shared_mem_byte, Gpu::gpuStream(),
                              [=] AMREX_GPU_DEVICE (Gpu::Handler const& h) noexcept
#else
                amrex::launch(nblocks, AMREX_GPU_MAX_THREADS, Gpu::gpuStream(),
                              [=] AMREX_GPU_DEVICE () noexcept
#endif
                {
#ifdef AMREX_USE_SYCL
                    int i1d = h.blockIdx() / n2dblocks;
                    int i2d = h.threadIdx() + h.blockDim()*(h.blockIdx()-i1d*n2dblocks);
#else
                    int i1d = blockIdx.x / n2dblocks;
                    int i2d = threadIdx.x + blockDim.x*(blockIdx.x-i1d*n2dblocks);
#endif
                    int i2dy = i2d / n2dx;
                    int i2dx = i2d - i2dy*n2dx;
                    int i, j, k, idir;
                    if (direction == 0) {
                        i = i1d  + lo.x;
                        j = i2dx + lo.y;
                        k = i2dy + lo.z;
                        idir = i;
                    } else if (direction == 1) {
                        i = i2dx + lo.x;
                        j = i1d  + lo.y;
                        k = i2dy + lo.z;
                        idir = j;
                    } else {
                        i = i2dx + lo.x;
                        j = i2dy + lo.y;
                        k = i1d  + lo.z;
                        idir = k;
                    }
                    for (int n = 0; n < ncomp; ++n) {
                        Real r = (i2d < n2d) ? fab(i,j,k,n+icomp) : Real(0.0);
#ifdef AMREX_USE_SYCL
                        Gpu::deviceReduceSum_full(p+n+ncomp*idir, r, h);
#else
                        Gpu::deviceReduceSum_full(p+n+ncomp*idir, r);
#endif
                    }
                });
            }

            Gpu::copyAsync(Gpu::deviceToHost, dv.begin(), dv.end(), hv.begin());
            Gpu::streamSynchronize();
        }
        else
#endif
        {
            for (auto& x : hv) {
                x = Real(0.0);
            }

            Vector<Gpu::HostVector<Real> > other_hv(OpenMP::get_max_threads()-1);

            Vector<Real*> pp(OpenMP::get_max_threads());
            if (!pp.empty()) {
                pp[0] = hv.data();
            }
            for (int i = 1; i < OpenMP::get_max_threads(); ++i) {
                other_hv[i-1].resize(n1d, Real(0.0));
                pp[i] = other_hv[i-1].data();
            }

#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
            for (MFIter mfi(mf,true); mfi.isValid(); ++mfi) {
                Box const& b = mfi.tilebox();
                auto const& fab = mf.const_array(mfi);
                Real * AMREX_RESTRICT p = pp[OpenMP::get_thread_num()];
                if (direction == 0) {
                    amrex::LoopOnCpu(b, ncomp, [&] (int i, int j, int k, int n) noexcept
                    {
                        p[n+ncomp*i] += fab(i,j,k,n+icomp);
                    });
                } else if (direction == 1) {
                    amrex::LoopOnCpu(b, ncomp, [&] (int i, int j, int k, int n) noexcept
                    {
                        p[n+ncomp*j] += fab(i,j,k,n+icomp);
                    });
                } else {
                    amrex::LoopOnCpu(b, ncomp, [&] (int i, int j, int k, int n) noexcept
                    {
                        p[n+ncomp*k] += fab(i,j,k,n+icomp);
                    });
                }
            }

            if (! other_hv.empty()) {
#ifdef AMREX_USE_OMP
#pragma omp parallel for
#endif
                for (int i = 0; i < n1d; ++i) {
                    for (auto const& v : other_hv) {
                        hv[i] += v[i];
                    }
                }
            }
        }

        if (!local) {
            ParallelAllReduce::Sum(hv.data(), static_cast<int>(hv.size()),
                                   ParallelContext::CommunicatorSub());
        }
        return hv;
    }

    Real volumeWeightedSum (Vector<MultiFab const*> const& mf, int icomp,
                            Vector<Geometry> const& geom,
                            Vector<IntVect> const& ratio,
                            bool local)
    {
        ReduceOps<ReduceOpSum> reduce_op;
        ReduceData<Real> reduce_data(reduce_op);

#ifdef AMREX_USE_EB
        bool has_eb = !(mf[0]->isAllRegular());
#endif

        int nlevels = static_cast<int>(mf.size());
        for (int ilev = 0; ilev < nlevels-1; ++ilev) {
            iMultiFab mask = makeFineMask(*mf[ilev], *mf[ilev+1], IntVect(0),
                                          ratio[ilev],Periodicity::NonPeriodic(),
                                          0, 1);
            auto const& m = mask.const_arrays();
            auto const& a = mf[ilev]->const_arrays();
            auto const dx = geom[ilev].CellSizeArray();
            Real dv = AMREX_D_TERM(dx[0],*dx[1],*dx[2]);
#ifdef AMREX_USE_EB
            if (has_eb) {
                AMREX_ASSERT(mf[ilev]->hasEBFabFactory());
                auto const& f = dynamic_cast<EBFArrayBoxFactory const&>
                    (mf[ilev]->Factory());
                auto const& vfrac = f.getVolFrac();
                auto const& va = vfrac.const_arrays();
                reduce_op.eval(*mf[ilev], IntVect(0), reduce_data,
                [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
                               -> Real
                {
                    return m[box_no](i,j,k) ? Real(0.)
                        : dv*a[box_no](i,j,k,icomp)*va[box_no](i,j,k);
                });
            } else
#endif
            {
#if (AMREX_SPACEDIM == 1)
                if (geom[ilev].IsSPHERICAL()) {
                    const auto rlo = geom[ilev].CellSize(0);
                    reduce_op.eval(*mf[ilev], IntVect(0), reduce_data,
                    [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k)
                                   noexcept -> Real
                    {
                        if (m[box_no](i,j,k)) {
                            return Real(0.);
                        } else {
                            constexpr Real pi = Real(3.1415926535897932);
                            Real ri = rlo + dx[0]*i;
                            Real ro = ri + dx[0];
                            return Real(4./3.)*pi*(ro-ri)*(ro*ro+ro*ri+ri*ri)
                                * a[box_no](i,j,k,icomp);
                        }
                    });
                } else
#elif (AMREX_SPACEDIM == 2)
                if (geom[ilev].IsRZ()) {
                    const auto rlo = geom[ilev].CellSize(0);
                    reduce_op.eval(*mf[ilev], IntVect(0), reduce_data,
                    [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k)
                                   noexcept -> Real
                    {
                        if (m[box_no](i,j,k)) {
                            return Real(0.);
                        } else {
                            Real ri = rlo + dx[0]*i;
                            Real ro = ri + dx[0];
                            constexpr Real pi = Real(3.1415926535897932);
                            return pi*dx[1]*dx[0]*(ro+ri)
                                * a[box_no](i,j,k,icomp);
                        }
                    });
                } else
#endif
                {
                    reduce_op.eval(*mf[ilev], IntVect(0), reduce_data,
                    [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k)
                                   noexcept -> Real
                    {
                        return m[box_no](i,j,k) ? Real(0.)
                            : dv*a[box_no](i,j,k,icomp);
                    });
                }
            }
            Gpu::streamSynchronize();
        }

        auto const& a = mf.back()->const_arrays();
        auto const dx = geom[nlevels-1].CellSizeArray();
        Real dv = AMREX_D_TERM(dx[0],*dx[1],*dx[2]);
#ifdef AMREX_USE_EB
        if (has_eb) {
            AMREX_ASSERT(mf.back()->hasEBFabFactory());
            auto const& f = dynamic_cast<EBFArrayBoxFactory const&>
                (mf.back()->Factory());
            auto const& vfrac = f.getVolFrac();
            auto const& va = vfrac.const_arrays();
            reduce_op.eval(*mf.back(), IntVect(0), reduce_data,
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
                           -> Real
            {
                return dv*a[box_no](i,j,k,icomp)*va[box_no](i,j,k);
            });
        } else
#endif
        {
#if (AMREX_SPACEDIM == 1)
            if (geom[nlevels-1].IsSPHERICAL()) {
                const auto rlo = geom[nlevels-1].CellSize(0);
                reduce_op.eval(*mf.back(), IntVect(0), reduce_data,
                [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k)
                               noexcept -> Real
                {
                    constexpr Real pi = Real(3.1415926535897932);
                    Real ri = rlo + dx[0]*i;
                    Real ro = ri + dx[0];
                    return Real(4./3.)*pi*(ro-ri)*(ro*ro+ro*ri+ri*ri)
                        * a[box_no](i,j,k,icomp);
                });
            } else
#elif (AMREX_SPACEDIM == 2)
            if (geom[nlevels-1].IsRZ()) {
                const auto rlo = geom[nlevels-1].CellSize(0);
                reduce_op.eval(*mf.back(), IntVect(0), reduce_data,
                [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k)
                               noexcept -> Real
                {
                    Real ri = rlo + dx[0]*i;
                    Real ro = ri + dx[0];
                    constexpr Real pi = Real(3.1415926535897932);
                    return pi*dx[1]*dx[0]*(ro+ri)
                        * a[box_no](i,j,k,icomp);
                });
            } else
#endif
            {
                reduce_op.eval(*mf.back(), IntVect(0), reduce_data,
                [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
                {
                    return dv*a[box_no](i,j,k,icomp);
                });
            }
        }

        auto const& hv = reduce_data.value(reduce_op);
        Real r = amrex::get<0>(hv);

        if (!local) {
            ParallelAllReduce::Sum(r, ParallelContext::CommunicatorSub());
        }
        return r;
    }

    void FourthOrderInterpFromFineToCoarse (MultiFab& cmf, int scomp, int ncomp,
                                            MultiFab const& fmf,
                                            IntVect const& ratio)
    {
        AMREX_ASSERT(AMREX_D_TERM(   (ratio[0] == 2 || ratio[0] == 4),
                                  && (ratio[1] == 2 || ratio[1] == 4),
                                  && (ratio[2] == 2 || ratio[2] == 4)));

        MultiFab tmp(amrex::coarsen(fmf.boxArray(), ratio), fmf.DistributionMap(),
                     ncomp, 0);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        {
#if (AMREX_SPACEDIM > 1)
            FArrayBox xtmp;
#if (AMREX_SPACEDIM > 2)
            FArrayBox ytmp;
#endif
#endif
            for (MFIter mfi(tmp,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                Box const& bx = mfi.tilebox();
                auto const& fa = fmf.const_array(mfi,scomp);

                Box xbx = bx;
#if (AMREX_SPACEDIM == 1)
                auto const& xa = tmp.array(mfi);
#else
                xbx.refine(IntVect(AMREX_D_DECL(1,ratio[1],ratio[2])));
                if (ratio[1] == 2) { xbx.grow(1,1); }
#if (AMREX_SPACEDIM == 3)
                if (ratio[2] == 2) { xbx.grow(2,1); }
#endif
                xtmp.resize(xbx,ncomp);
                Elixir eli = xtmp.elixir();
                auto const& xa = xtmp.array();
#endif
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D(xbx, ncomp, i, j, k, n,
                {
                    int ii = 2*i;
                    xa(i,j,k,n) = Real(1./16)*(Real(9.)*(fa(ii  ,j,k,n) +
                                                         fa(ii+1,j,k,n))
                                               -         fa(ii-1,j,k,n)
                                               -         fa(ii+2,j,k,n));
                });

#if (AMREX_SPACEDIM > 1)
                Box ybx = bx;
                auto const& xca = xtmp.const_array();
#if (AMREX_SPACEDIM == 2)
                auto const& ya = tmp.array(mfi);
#else
                ybx.refine(IntVect(AMREX_D_DECL(1,1,ratio[2])));
                if (ratio[2] == 2) { ybx.grow(2,1); }
                ytmp.resize(ybx,ncomp);
                eli.append(ytmp.elixir());
                auto const& ya = ytmp.array();
#endif
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D(ybx, ncomp, i, j, k, n,
                {
                    int jj = 2*j;
                    ya(i,j,k,n) = Real(1./16)*(Real(9.)*(xca(i,jj  ,k,n) +
                                                         xca(i,jj+1,k,n))
                                               -         xca(i,jj-1,k,n)
                                               -         xca(i,jj+2,k,n));
                });

#if (AMREX_SPACEDIM == 3)
                auto const& yca = ytmp.const_array();
                auto const& ca = tmp.array(mfi);
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D(bx, ncomp, i, j, k, n,
                {
                    int kk = 2*k;
                    ca(i,j,k,n) = Real(1./16)*(Real(9.)*(yca(i,j,kk  ,n) +
                                                         yca(i,j,kk+1,n))
                                               -         yca(i,j,kk-1,n)
                                               -         yca(i,j,kk+2,n));
                });
#endif
#endif
            }
        }

        cmf.ParallelCopy(tmp, 0, scomp, ncomp);
    }
}
