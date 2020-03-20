
#include <AMReX_MultiFabUtil.H>
#include <sstream>
#include <iostream>

namespace {

    using namespace amrex;

    static Box
    getIndexBox(const RealBox& real_box, const Geometry& geom) {
        IntVect slice_lo, slice_hi;

        AMREX_D_TERM(slice_lo[0]=std::floor((real_box.lo(0) - geom.ProbLo(0))/geom.CellSize(0));,
                     slice_lo[1]=std::floor((real_box.lo(1) - geom.ProbLo(1))/geom.CellSize(1));,
                     slice_lo[2]=std::floor((real_box.lo(2) - geom.ProbLo(2))/geom.CellSize(2)););

        AMREX_D_TERM(slice_hi[0]=std::floor((real_box.hi(0) - geom.ProbLo(0))/geom.CellSize(0));,
                     slice_hi[1]=std::floor((real_box.hi(1) - geom.ProbLo(1))/geom.CellSize(1));,
                     slice_hi[2]=std::floor((real_box.hi(2) - geom.ProbLo(2))/geom.CellSize(2)););

        return Box(slice_lo, slice_hi) & geom.Domain();
    }

    static
    std::unique_ptr<MultiFab> allocateSlice(int dir, const MultiFab& cell_centered_data,
                                            int ncomp, const Geometry& geom, Real dir_coord,
                                            Vector<int>& slice_to_full_ba_map) {

        // Get our slice and convert to index space
        RealBox real_slice = geom.ProbDomain();
        real_slice.setLo(dir, dir_coord);
        real_slice.setHi(dir, dir_coord);
        Box slice_box = getIndexBox(real_slice, geom);

        // define the multifab that stores slice
        BoxArray ba = cell_centered_data.boxArray();
        const DistributionMapping& dm = cell_centered_data.DistributionMap();
        std::vector< std::pair<int, Box> > isects;
        ba.intersections(slice_box, isects, false, 0);
        Vector<Box> boxes;
        Vector<int> procs;
        for (int i = 0; i < isects.size(); ++i) {
            procs.push_back(dm[isects[i].first]);
            boxes.push_back(isects[i].second);
            slice_to_full_ba_map.push_back(isects[i].first);
        }
        BoxArray slice_ba(&boxes[0], boxes.size());
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
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(cc,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box bx = mfi.growntilebox(ngrow);
            Array4<Real> const& ccarr = cc.array(mfi);
            Array4<Real const> const& ndarr = nd.const_array(mfi);

            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
            {
                amrex_avg_nd_to_cc(tbx, ccarr, ndarr, dcomp, scomp, ncomp);
            });
        }
    }

    void average_edge_to_cellcenter (MultiFab& cc, int dcomp,
        const Vector<const MultiFab*>& edge, int ngrow)
    {
        AMREX_ASSERT(cc.nComp() >= dcomp + AMREX_SPACEDIM);
        AMREX_ASSERT(edge.size() == AMREX_SPACEDIM);
        AMREX_ASSERT(edge[0]->nComp() == 1);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(cc,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box bx = mfi.growntilebox(ngrow);
            Array4<Real> const& ccarr = cc.array(mfi);
            AMREX_D_TERM(Array4<Real const> const& exarr = edge[0]->const_array(mfi);,
                         Array4<Real const> const& eyarr = edge[1]->const_array(mfi);,
                         Array4<Real const> const& ezarr = edge[2]->const_array(mfi););

            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
            {
                amrex_avg_eg_to_cc(tbx, ccarr, AMREX_D_DECL(exarr,eyarr,ezarr), dcomp);
            });
        }
    }

    void average_face_to_cellcenter (MultiFab& cc, int dcomp,
        const Vector<const MultiFab*>& fc, int ngrow)
    {
        average_face_to_cellcenter(cc, dcomp,
            Array<MultiFab const*,AMREX_SPACEDIM>{AMREX_D_DECL(fc[0],fc[1],fc[2])},
            ngrow);
    }

    void average_face_to_cellcenter (MultiFab& cc, const Vector<const MultiFab*>& fc,
				     const Geometry& geom)
    {
        average_face_to_cellcenter(cc,
                                   Array<MultiFab const*,AMREX_SPACEDIM>
                                                  {AMREX_D_DECL(fc[0],fc[1],fc[2])},
                                   geom);
    }

    void average_face_to_cellcenter (MultiFab& cc, int dcomp,
                                     const Array<const MultiFab*,AMREX_SPACEDIM>& fc, int ngrow)
    {
        AMREX_ASSERT(cc.nComp() >= dcomp + AMREX_SPACEDIM);
        AMREX_ASSERT(fc[0]->nComp() == 1);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(cc,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box bx = mfi.growntilebox(ngrow);
            Array4<Real> const& ccarr = cc.array(mfi);
            AMREX_D_TERM(Array4<Real const> const& fxarr = fc[0]->const_array(mfi);,
                         Array4<Real const> const& fyarr = fc[1]->const_array(mfi);,
                         Array4<Real const> const& fzarr = fc[2]->const_array(mfi););

            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
            {
#if (AMREX_SPACEDIM == 1)
                amrex_avg_fc_to_cc(tbx, ccarr, AMREX_D_DECL(fxarr,fyarr,fzarr), dcomp, GeometryData());
#else
                amrex_avg_fc_to_cc(tbx, ccarr, AMREX_D_DECL(fxarr,fyarr,fzarr), dcomp);
#endif
            });
        }
    }

    void average_face_to_cellcenter (MultiFab& cc,
                                     const Array<const MultiFab*,AMREX_SPACEDIM>& fc,
				     const Geometry& geom)
    {
	AMREX_ASSERT(cc.nComp() >= AMREX_SPACEDIM);
	AMREX_ASSERT(fc[0]->nComp() == 1); // We only expect fc to have the gradient perpendicular to the face

        const GeometryData gd = geom.data();
        amrex::ignore_unused(gd);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(cc,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box bx = mfi.tilebox();
            Array4<Real> const& ccarr = cc.array(mfi);
            AMREX_D_TERM(Array4<Real const> const& fxarr = fc[0]->const_array(mfi);,
                         Array4<Real const> const& fyarr = fc[1]->const_array(mfi);,
                         Array4<Real const> const& fzarr = fc[2]->const_array(mfi););

            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
            {
#if (AMREX_SPACEDIM == 1)
                amrex_avg_fc_to_cc(tbx, ccarr, AMREX_D_DECL(fxarr,fyarr,fzarr), 0, gd);
#else
                amrex_avg_fc_to_cc(tbx, ccarr, AMREX_D_DECL(fxarr,fyarr,fzarr), 0);
#endif
            });
        }
    }

    void average_cellcenter_to_face (const Vector<MultiFab*>& fc, const MultiFab& cc,
				     const Geometry& geom)
    {
        average_cellcenter_to_face(Array<MultiFab*,AMREX_SPACEDIM>{AMREX_D_DECL(fc[0],fc[1],fc[2])},
                                   cc, geom);
    }


    void average_cellcenter_to_face (const Array<MultiFab*,AMREX_SPACEDIM>& fc, const MultiFab& cc,
                                    const Geometry& geom)
    {
	AMREX_ASSERT(cc.nComp() == 1);
	AMREX_ASSERT(cc.nGrow() >= 1);
	AMREX_ASSERT(fc[0]->nComp() == 1); // We only expect fc to have the gradient perpendicular to the face

#if (AMREX_SPACEDIM == 1)
        const GeometryData& gd = geom.data();
#endif

#ifdef _OPENMP
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

            AMREX_LAUNCH_HOST_DEVICE_LAMBDA (index_bounds, tbx,
            {
#if (AMREX_SPACEDIM == 1)
                amrex_avg_cc_to_fc(tbx, AMREX_D_DECL(xbx,ybx,zbx),
                                   AMREX_D_DECL(fxarr,fyarr,fzarr), ccarr, gd);
#else
                amrex_avg_cc_to_fc(tbx, AMREX_D_DECL(xbx,ybx,zbx),
                                   AMREX_D_DECL(fxarr,fyarr,fzarr), ccarr);
#endif
            });
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
        BL_PROFILE("amrex::average_down_w_geom");

        if (S_fine.is_nodal() || S_crse.is_nodal())
        {
            amrex::Error("Can't use amrex::average_down for nodal MultiFab!");
        }

#if (AMREX_SPACEDIM == 3)
	amrex::average_down(S_fine, S_crse, scomp, ncomp, ratio);
	return;
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

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(crse_S_fine,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            //  NOTE: The tilebox is defined at the coarse level.
            const Box& bx = mfi.tilebox();
            Array4<Real> const& crsearr = crse_S_fine.array(mfi);
            Array4<Real const> const& finearr = S_fine.const_array(mfi);
            Array4<Real const> const& finevolarr = fvolume.const_array(mfi);

            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
            {
                amrex_avgdown_with_vol(tbx,crsearr,finearr,finevolarr,
                                       0,scomp,ncomp,ratio);
            });
	}

        S_crse.copy(crse_S_fine,0,scomp,ncomp);
#endif
   }

// *************************************************************************************************************

    // Average fine cell-based MultiFab onto crse cell-centered MultiFab.
    // We do NOT assume that the coarse layout is a coarsened version of the fine layout.
    // This version does NOT use volume-weighting
    void average_down (const MultiFab& S_fine, MultiFab& S_crse, int scomp, int ncomp, int rr)
    {
         average_down(S_fine,S_crse,scomp,ncomp,rr*IntVect::TheUnitVector());
    }


    void sum_fine_to_coarse(const MultiFab& S_fine, MultiFab& S_crse,
                            int scomp, int ncomp, const IntVect& ratio,
                            const Geometry& cgeom, const Geometry& fgeom)
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

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(crse_S_fine, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            //  NOTE: The tilebox is defined at the coarse level.
            const Box& bx = mfi.growntilebox(nGrow);
            Array4<Real> const& crsearr = crse_S_fine.array(mfi);
            Array4<Real const> const& finearr = S_fine.const_array(mfi);

            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
            {
                amrex_avgdown(tbx,crsearr,finearr,0,scomp,ncomp,ratio);
            });
        }

        S_crse.copy(crse_S_fine, 0, scomp, ncomp, nGrow, 0,
                    cgeom.periodicity(), FabArrayBase::ADD);
    }

    void average_down (const MultiFab& S_fine, MultiFab& S_crse,
                       int scomp, int ncomp, const IntVect& ratio)
    {
        BL_PROFILE("amrex::average_down");
        AMREX_ASSERT(S_crse.nComp() == S_fine.nComp());
        AMREX_ASSERT((S_crse.is_cell_centered() && S_fine.is_cell_centered()) ||
                     (S_crse.is_nodal()         && S_fine.is_nodal()));

        bool is_cell_centered = S_crse.is_cell_centered();

        //
        // Coarsen() the fine stuff on processors owning the fine data.
        //
        BoxArray crse_S_fine_BA = S_fine.boxArray(); crse_S_fine_BA.coarsen(ratio);

        if (crse_S_fine_BA == S_crse.boxArray() and S_fine.DistributionMap() == S_crse.DistributionMap())
        {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(S_crse,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                //  NOTE: The tilebox is defined at the coarse level.
                const Box& bx = mfi.tilebox();
                Array4<Real> const& crsearr = S_crse.array(mfi);
                Array4<Real const> const& finearr = S_fine.const_array(mfi);

                if (is_cell_centered) {
                    AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
                    {
                        amrex_avgdown(tbx,crsearr,finearr,scomp,scomp,ncomp,ratio);
                    });
                } else {
                    AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
                    {
                        amrex_avgdown_nodes(tbx,crsearr,finearr,scomp,scomp,ncomp,ratio);
                    });
                }
            }
        }
        else
        {
            MultiFab crse_S_fine(crse_S_fine_BA, S_fine.DistributionMap(), ncomp, 0, MFInfo(), FArrayBoxFactory());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(crse_S_fine,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                //  NOTE: The tilebox is defined at the coarse level.
                const Box& bx = mfi.tilebox();
                Array4<Real> const& crsearr = crse_S_fine.array(mfi);
                Array4<Real const> const& finearr = S_fine.const_array(mfi);

                //  NOTE: We copy from component scomp of the fine fab into component 0 of the crse fab
                //        because the crse fab is a temporary which was made starting at comp 0, it is
                //        not part of the actual crse multifab which came in.

                if (is_cell_centered) {
                    AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
                    {
                        amrex_avgdown(tbx,crsearr,finearr,0,scomp,ncomp,ratio);
                    });
                } else {
                    AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
                    {
                        amrex_avgdown_nodes(tbx,crsearr,finearr,0,scomp,ncomp,ratio);
                    });
                }
            }

            S_crse.copy(crse_S_fine,0,scomp,ncomp);
        }
   }

// *************************************************************************************************************

    void average_down_faces (const Vector<const MultiFab*>& fine,
                             const Vector<MultiFab*>& crse,
                             const IntVect& ratio, int ngcrse)
    {
        average_down_faces(Array<const MultiFab*,AMREX_SPACEDIM>
                                   {AMREX_D_DECL(fine[0],fine[1],fine[2])},
                           Array<MultiFab*,AMREX_SPACEDIM>
                                   {AMREX_D_DECL(crse[0],crse[1],crse[2])},
                           ratio, ngcrse);
    }

    void average_down_faces (const Vector<const MultiFab*>& fine,
                             const Vector<MultiFab*>& crse, int ratio, int ngcrse)
    {
        average_down_faces(fine,crse,IntVect{ratio},ngcrse);
    }

    void average_down_faces (const Array<const MultiFab*,AMREX_SPACEDIM>& fine,
                             const Array<MultiFab*,AMREX_SPACEDIM>& crse,
                             int ratio, int ngcrse)
    {
        average_down_faces(fine,crse,IntVect{ratio},ngcrse);
    }

    // Average fine face-based MultiFab onto crse face-based MultiFab.
    void average_down_faces (const Array<const MultiFab*,AMREX_SPACEDIM>& fine,
                             const Array<MultiFab*,AMREX_SPACEDIM>& crse,
			     const IntVect& ratio, int ngcrse)
    {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            average_down_faces(*fine[idim], *crse[idim], ratio, ngcrse);
        }
    }

    void average_down_faces (const MultiFab& fine, MultiFab& crse,
                             const IntVect& ratio, int ngcrse)
    {
	AMREX_ASSERT(crse.nComp() == fine.nComp());
        AMREX_ASSERT(fine.ixType() == crse.ixType());
        const auto type = fine.ixType();
        int dir;
        for (dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            if (type.nodeCentered(dir)) break;
        }
        auto tmptype = type;
        tmptype.unset(dir);
        if (dir >= AMREX_SPACEDIM or !tmptype.cellCentered()) {
            amrex::Abort("average_down_faces: not face index type");
        }
        const int ncomp = crse.nComp();
        if (isMFIterSafe(fine, crse))
        {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(crse,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.growntilebox(ngcrse);
                Array4<Real> const& crsearr = crse.array(mfi);
                Array4<Real const> const& finearr = fine.const_array(mfi);

                AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
                {
                    amrex_avgdown_faces(tbx, crsearr, finearr, 0, 0, ncomp, ratio, dir);
                });
            }
        }
        else
        {
            MultiFab ctmp(amrex::coarsen(fine.boxArray(),ratio), fine.DistributionMap(),
                          ncomp, ngcrse, MFInfo(), FArrayBoxFactory());
            average_down_faces(fine, ctmp, ratio, ngcrse);
            crse.ParallelCopy(ctmp,0,0,ncomp,ngcrse,ngcrse);
        }
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
        if (dir >= AMREX_SPACEDIM or !tmptype.nodeCentered()) {
            amrex::Abort("average_down_edges: not face index type");
        }
        const int ncomp = crse.nComp();
        if (isMFIterSafe(fine, crse))
        {
#ifdef _OPENMP
#pragma omp parallel if(Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(crse,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.growntilebox(ngcrse);
                Array4<Real> const& crsearr = crse.array(mfi);
                Array4<Real const> const& finearr = fine.const_array(mfi);

                AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
                {
                    amrex_avgdown_edges(tbx, crsearr, finearr, 0, 0, ncomp, ratio, dir);
                });
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

    FabArray<BaseFab<long> > ToLongMultiFab (const iMultiFab& imf)
    {
        return amrex::cast<FabArray<BaseFab<long> > > (imf);
    }

    std::unique_ptr<MultiFab> get_slice_data(int dir, Real coord, const MultiFab& cc, const Geometry& geom, int start_comp, int ncomp, bool interpolate) {

        BL_PROFILE("amrex::get_slice_data");

        if (interpolate) {
            AMREX_ASSERT(cc.nGrow() >= 1);
        }

        const auto geomdata = geom.data();

        Vector<int> slice_to_full_ba_map;
        std::unique_ptr<MultiFab> slice = allocateSlice(dir, cc, ncomp, geom, coord, slice_to_full_ba_map);

#ifdef _OPENMP
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
                AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( tile_box, thread_box,
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
#ifdef _OPENMP
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
                    for (const auto is : isects) {
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
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n, Array4<value_type> const& a) noexcept
        {
            a(i,j,k,n) = fine_value;
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

#ifdef _OPENMP
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

#ifdef _OPENMP
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
        MultiFab nmf(nba, mf.DistributionMap(), mf.nComp(), 0,
                     MFInfo().SetAlloc(false));

        for (MFIter mfi(r); mfi.isValid(); ++mfi) {
            auto const& rfab = r[mfi];
            nmf.setFab(mfi, new FArrayBox(nba[mfi.index()], rfab.nComp(),
                                          rfab.dataPtr()), false);
        }

        nmf.ParallelCopy(mf, period);

        return r;
    }

#ifdef AMREX_USE_HDF5

    SortedGrids::SortedGrids()
    {
      return;
    }

    SortedGrids::SortedGrids(const BoxArray &grids, const DistributionMapping &dmap)
    {
      update(grids, dmap);
      return;
    }

    void SortedGrids::update(const BoxArray &grids, const DistributionMapping &dmap)
    {
      const Vector<int>& pMap = dmap.ProcessorMap();

      gridMap.clear();
      for (int i = 0; i < grids.size(); ++i) {
        int gridProc = pMap[i];
        Vector<std::pair<int, Box> >& boxesAtProc = gridMap[gridProc];
        boxesAtProc.push_back(std::pair<int, Box>(i, grids[i]));
      }

      // sorted grids and processor ids
      sortedGrids.resize(grids.size());
      sortedProcs.resize(grids.size());
      int bIndex = 0;
      for (auto it = gridMap.begin(); it != gridMap.end(); ++it) {
        int proc = it->first;
        Vector<std::pair<int, Box> >& boxesAtProc = it->second;
        for (int ii = 0; ii < boxesAtProc.size(); ++ii) {
          sortedGrids.set(bIndex, boxesAtProc[ii].second);
          sortedProcs[bIndex] = proc;
          ++bIndex;
        }
      }
      return;
    }

    void SortedGrids::update(MultiFab* data_mf)
    {
      int nComp = data_mf->nComp();

      // offset of data for each grid
      offsets.resize(sortedGrids.size() + 1);
      unsigned long long currentOffset = 0L;
      for (int b = 0; b < sortedGrids.size(); ++b) {
        offsets[b] = currentOffset;
        currentOffset += sortedGrids[b].numPts() * nComp;
      }
      offsets[sortedGrids.size()] = currentOffset;

      int nProcs(ParallelDescriptor::NProcs());

      // processor offsets
      procOffsets.resize(nProcs);
      procBufferSize.resize(nProcs);
      unsigned long long totalOffset = 0L;
      for (auto it = gridMap.begin(); it != gridMap.end(); ++it) {
        int proc = it->first;
        Vector<std::pair<int, Box> >& boxesAtProc = it->second;
        procOffsets[proc] = totalOffset;
        procBufferSize[proc] = 0L;
        for (int b = 0; b < boxesAtProc.size(); ++b) {
          procBufferSize[proc] += boxesAtProc[b].second.numPts() * nComp;
        }
        totalOffset += procBufferSize[proc];
      }
      return;
    }

    void SortedGrids::update(const Vector<long> count, const int size)
    {
      // offset of data for each grid
      offsets.resize(sortedGrids.size() + 1);
      unsigned long long currentOffset = 0L;
      for (int b = 0; b < sortedGrids.size(); ++b) {
        offsets[b] = currentOffset;
        currentOffset += count[b] * size;
      }
      offsets[sortedGrids.size()] = currentOffset;

      int nProcs(ParallelDescriptor::NProcs());

      // processor offsets
      procOffsets.resize(nProcs);
      procBufferSize.resize(nProcs);
      unsigned long long totalOffset = 0L;
      for (auto it = gridMap.begin(); it != gridMap.end(); ++it) {
        int proc = it->first;
        Vector<std::pair<int, Box> >& boxesAtProc = it->second;
        procOffsets[proc] = totalOffset;
        procBufferSize[proc] = 0L;
        for (int b = 0; b < boxesAtProc.size(); ++b) {
          procBufferSize[proc] += count[boxesAtProc[b].first] * size;
        }
        totalOffset += procBufferSize[proc];
      }
      return;
    }

    void writeMultiFab(H5& h5, MultiFab* data_mf, const Real time,
                           const std::vector<std::string> data_names)
    {

      int myProc(ParallelDescriptor::MyProc());
      int nComp = data_mf->nComp();

      const BoxArray& grids = data_mf->boxArray();
      const DistributionMapping& dmap = data_mf->DistributionMap();

      // make sure our multifab doesn't have any ghost cells
      MultiFab tmp_fab(data_mf->boxArray(), data_mf->DistributionMap(), nComp, 0, MFInfo(), data_mf->Factory());
      MultiFab::Copy(tmp_fab, *data_mf, 0, 0, nComp, 0);

      //-----------------------------------------------------------------------------------
      // multifabs data attributes

      std::map<std::string, int> vMInt;
      std::map<std::string, Real> vMReal;
      std::map<std::string, std::string> vMString;

      vMInt["comps"] = nComp;
      vMInt["num_names"] = data_names.size();
      for (size_t i=0; i<data_names.size(); ++i) {
          vMString["component_"+num2str(i)] = data_names[i];
        }
      vMReal["time"] = time;

      h5.writeAttribute(vMInt, vMReal, vMString);

      // the number of ghost cells saved and output
      IntVect nGrow = tmp_fab.nGrowVect();
      hid_t intvect_id = makeH5IntVec();
      int_h5_t gint = writeH5IntVec(nGrow.getVect());
      h5.writeAttribute("ghost", gint, intvect_id);
      h5.writeAttribute("outputGhost", gint, intvect_id);
      H5Tclose(intvect_id);

      // box layout etc
      SortedGrids sg(grids, dmap);
      sg.update(&tmp_fab);
      sg.sortedGrids.writeOnHDF5(h5, "boxes");

      // offsets
      h5.writeType("data:offsets=0", {(hsize_t)sg.offsets.size()}, sg.offsets,
                   H5T_NATIVE_LLONG);

      // first gather
      std::vector<Real> a_buffer(sg.procBufferSize[myProc], -1.0);
      long dataCount(0);
      for (MFIter mfi(tmp_fab); mfi.isValid(); ++mfi) {
          const Box& vbox = mfi.validbox();

          const Real* dataPtr = (tmp_fab)[mfi].dataPtr();
          for (int i(0); i < vbox.numPts() * nComp; ++i) {
              a_buffer[dataCount++] = dataPtr[i];
            }
        }

      // define how big and where
      std::vector<hsize_t> full_dims = {(hsize_t)sg.offsets.back()};
      std::vector<hsize_t> local_dims = {(hsize_t)sg.procBufferSize[myProc]};
      std::vector<hsize_t> offset = {(hsize_t)sg.procOffsets[myProc]};

      // then write
      h5.writeSlab("data:datatype=0", full_dims, local_dims, offset, a_buffer,
                   H5T_NATIVE_DOUBLE);

      return;

    }

    void readMultiFab(H5& h5, MultiFab* data_mf) {

      // get the boxes
      BoxArray boxes;
      boxes.readFromHDF5(h5, "boxes");

      // get the data offsets
      Vector<unsigned long long> offsets;
      h5.readType("data:offsets=0", offsets);

      int nComp = data_mf->nComp();

      std::vector<Real> buffer;
      std::vector<hsize_t> local_dims(1);
      std::vector<hsize_t> offset(1);

      // make a temporary multifab
      MultiFab tmp_fab(data_mf->boxArray(), data_mf->DistributionMap(), nComp, 0, MFInfo(), data_mf->Factory());

      // retrieve data and load into multifab
      for (MFIter mfi(tmp_fab); mfi.isValid(); ++mfi) {
          const Box& vbox = mfi.validbox();

          // get offset in the hdf5 file for the current box
          for (int i=0; i<boxes.size(); ++i) {
            if (boxes.get(i) == vbox) {
              offset[0] = offsets[i];
              local_dims[0] = offsets[i+1] - offsets[i];
              break;
            }
            if (i == boxes.size()-1) {
              amrex::Abort("Error: readMultiFab() unable to find correct box in hdf5 file");
            }
          }

          // get the data from the hdf5 file
          h5.readSlab("data:datatype=0", local_dims, offset, buffer);

          // load the data into the temporary multifab
          Real* dataPtr = tmp_fab[mfi].dataPtr();
          for (int i(0); i < vbox.numPts() * nComp; ++i) {
              dataPtr[i] = buffer[i];
            }
        }

        // now copy to the real fab
        MultiFab::Copy(*data_mf, tmp_fab, 0, 0, nComp, 0);
    }

#endif

}
