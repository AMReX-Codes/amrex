
#include <AMReX_MultiFabUtil.H>
#include <AMReX_MultiFabUtil_F.H>
#include <sstream>
#include <iostream>

namespace {

    using namespace amrex;
   
    static Box
    getIndexBox(const RealBox& real_box, const Geometry& geom) {
        IntVect slice_lo, slice_hi;
        
        D_TERM(slice_lo[0]=floor((real_box.lo(0) - geom.ProbLo(0))/geom.CellSize(0));,
               slice_lo[1]=floor((real_box.lo(1) - geom.ProbLo(1))/geom.CellSize(1));,
               slice_lo[2]=floor((real_box.lo(2) - geom.ProbLo(2))/geom.CellSize(2)););
        
        D_TERM(slice_hi[0]=floor((real_box.hi(0) - geom.ProbLo(0))/geom.CellSize(0));,
               slice_hi[1]=floor((real_box.hi(1) - geom.ProbLo(1))/geom.CellSize(1));,
               slice_hi[2]=floor((real_box.hi(2) - geom.ProbLo(2))/geom.CellSize(2)););
        
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
        DistributionMapping slice_dmap(procs);
        std::unique_ptr<MultiFab> slice(new MultiFab(slice_ba, slice_dmap, ncomp, 0,
                                                     MFInfo(), cell_centered_data.Factory()));
        return slice;
    }
}

namespace amrex
{
    void average_node_to_cellcenter (MultiFab& cc, int dcomp, const MultiFab& nd, int scomp, int ncomp)
    {
#ifdef _OPENMP
#pragma omp parallel
#endif
	for (MFIter mfi(cc,true); mfi.isValid(); ++mfi) 
	{
            const Box& bx = mfi.tilebox();
            amrex_fort_avg_nd_to_cc(bx.loVect(), bx.hiVect(), &ncomp,
                                    BL_TO_FORTRAN_N(cc[mfi],dcomp),
                                    BL_TO_FORTRAN_N(nd[mfi],scomp));
        }
    }

    void average_edge_to_cellcenter (MultiFab& cc, int dcomp, const Vector<const MultiFab*>& edge)
    {
	BL_ASSERT(cc.nComp() >= dcomp + AMREX_SPACEDIM);
	BL_ASSERT(edge.size() == AMREX_SPACEDIM);
	BL_ASSERT(edge[0]->nComp() == 1);
#ifdef _OPENMP
#pragma omp parallel
#endif
	for (MFIter mfi(cc,true); mfi.isValid(); ++mfi) 
	{
	    const Box& bx = mfi.tilebox();

	    BL_FORT_PROC_CALL(BL_AVG_EG_TO_CC,bl_avg_eg_to_cc)
		(bx.loVect(), bx.hiVect(),
		 BL_TO_FORTRAN_N(cc[mfi],dcomp),
		 AMREX_D_DECL(BL_TO_FORTRAN((*edge[0])[mfi]),
			BL_TO_FORTRAN((*edge[1])[mfi]),
			BL_TO_FORTRAN((*edge[2])[mfi])));
	}	
    }

    void average_face_to_cellcenter (MultiFab& cc, int dcomp, const Vector<const MultiFab*>& fc)
    {
	BL_ASSERT(cc.nComp() >= dcomp + AMREX_SPACEDIM);
	BL_ASSERT(fc.size() == AMREX_SPACEDIM);
	BL_ASSERT(fc[0]->nComp() == 1);

	Real dx[3] = {1.0,1.0,1.0};
	Real problo[3] = {0.,0.,0.};
	int coord_type = 0;

#ifdef _OPENMP
#pragma omp parallel
#endif
	for (MFIter mfi(cc,true); mfi.isValid(); ++mfi) 
	{
	    const Box& bx = mfi.tilebox();

	    BL_FORT_PROC_CALL(BL_AVG_FC_TO_CC,bl_avg_fc_to_cc)
		(bx.loVect(), bx.hiVect(),
		 BL_TO_FORTRAN_N(cc[mfi],dcomp),
		 AMREX_D_DECL(BL_TO_FORTRAN((*fc[0])[mfi]),
			BL_TO_FORTRAN((*fc[1])[mfi]),
			BL_TO_FORTRAN((*fc[2])[mfi])),
		 dx, problo, coord_type);
	}
    }

    void average_face_to_cellcenter (MultiFab& cc, const Vector<const MultiFab*>& fc,
				     const Geometry& geom)
    {
	BL_ASSERT(cc.nComp() >= AMREX_SPACEDIM);
	BL_ASSERT(fc.size() == AMREX_SPACEDIM);
	BL_ASSERT(fc[0]->nComp() == 1); // We only expect fc to have the gradient perpendicular to the face

	const Real* dx     = geom.CellSize();
	const Real* problo = geom.ProbLo();
	int coord_type = Geometry::Coord();

#ifdef _OPENMP
#pragma omp parallel
#endif
	for (MFIter mfi(cc,true); mfi.isValid(); ++mfi) 
	{
	    const Box& bx = mfi.tilebox();

	    BL_FORT_PROC_CALL(BL_AVG_FC_TO_CC,bl_avg_fc_to_cc)
		(bx.loVect(), bx.hiVect(),
		 BL_TO_FORTRAN(cc[mfi]),
		 AMREX_D_DECL(BL_TO_FORTRAN((*fc[0])[mfi]),
			BL_TO_FORTRAN((*fc[1])[mfi]),
			BL_TO_FORTRAN((*fc[2])[mfi])),
		 dx, problo, coord_type);
	}
    }

    void average_cellcenter_to_face (const Vector<MultiFab*>& fc, const MultiFab& cc,
				     const Geometry& geom)
    {
	BL_ASSERT(cc.nComp() == 1);
	BL_ASSERT(cc.nGrow() >= 1);
	BL_ASSERT(fc.size() == AMREX_SPACEDIM);
	BL_ASSERT(fc[0]->nComp() == 1); // We only expect fc to have the gradient perpendicular to the face

	const Real* dx     = geom.CellSize();
	const Real* problo = geom.ProbLo();
	int coord_type = Geometry::Coord();

#ifdef _OPENMP
#pragma omp parallel
#endif
	for (MFIter mfi(cc,true); mfi.isValid(); ++mfi) 
	{
	    const Box& xbx = mfi.nodaltilebox(0);
#if (AMREX_SPACEDIM > 1)
	    const Box& ybx = mfi.nodaltilebox(1);
#endif
#if (AMREX_SPACEDIM == 3)
	    const Box& zbx = mfi.nodaltilebox(2);
#endif
	    
	    BL_FORT_PROC_CALL(BL_AVG_CC_TO_FC,bl_avg_cc_to_fc)
		(xbx.loVect(), xbx.hiVect(),
#if (AMREX_SPACEDIM > 1)
		 ybx.loVect(), ybx.hiVect(),
#endif
#if (AMREX_SPACEDIM == 3)
		 zbx.loVect(), zbx.hiVect(),
#endif
		 AMREX_D_DECL(BL_TO_FORTRAN((*fc[0])[mfi]),
			BL_TO_FORTRAN((*fc[1])[mfi]),
			BL_TO_FORTRAN((*fc[2])[mfi])),
		 BL_TO_FORTRAN(cc[mfi]),
		 dx, problo, coord_type);
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
  
        if (S_fine.is_nodal() || S_crse.is_nodal())
        {
            amrex::Error("Can't use amrex::average_down for nodal MultiFab!");
        }

#if (AMREX_SPACEDIM == 3)
	amrex::average_down(S_fine, S_crse, scomp, ncomp, ratio);
	return;
#else

        BL_ASSERT(S_crse.nComp() == S_fine.nComp());

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
#pragma omp parallel
#endif
        for (MFIter mfi(crse_S_fine,true); mfi.isValid(); ++mfi)
        {
            //  NOTE: The tilebox is defined at the coarse level.
            const Box& tbx = mfi.tilebox();

            //  NOTE: We copy from component scomp of the fine fab into component 0 of the crse fab
            //        because the crse fab is a temporary which was made starting at comp 0, it is
            //        not part of the actual crse multifab which came in.
            BL_FORT_PROC_CALL(BL_AVGDOWN_WITH_VOL,bl_avgdown_with_vol)
                (tbx.loVect(), tbx.hiVect(),
                 BL_TO_FORTRAN_N(S_fine[mfi],scomp),
                 BL_TO_FORTRAN_N(crse_S_fine[mfi],0),
                 BL_TO_FORTRAN(fvolume[mfi]),
                 ratio.getVect(),&ncomp);
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
        BL_ASSERT(S_crse.nComp() == S_fine.nComp());
        BL_ASSERT(ratio == ratio[0]);
        BL_ASSERT(S_fine.nGrow() % ratio[0] == 0);

        const int nGrow = S_fine.nGrow() / ratio[0];

        //
        // Coarsen() the fine stuff on processors owning the fine data.
        //
        BoxArray crse_S_fine_BA = S_fine.boxArray(); crse_S_fine_BA.coarsen(ratio);

        MultiFab crse_S_fine(crse_S_fine_BA, S_fine.DistributionMap(), ncomp, nGrow, MFInfo(), FArrayBoxFactory());

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(crse_S_fine, true); mfi.isValid(); ++mfi)
        {
            //  NOTE: The tilebox is defined at the coarse level.
            const Box& tbx = mfi.growntilebox(nGrow);
            
            BL_FORT_PROC_CALL(BL_AVGDOWN, bl_avgdown)
                (tbx.loVect(), tbx.hiVect(),
                 BL_TO_FORTRAN_N(S_fine[mfi] , scomp),
                 BL_TO_FORTRAN_N(crse_S_fine[mfi], 0),
                 ratio.getVect(),&ncomp);
        }
        
        S_crse.copy(crse_S_fine, 0, scomp, ncomp, nGrow, 0,
                    cgeom.periodicity(), FabArrayBase::ADD);
    }

    void average_down (const MultiFab& S_fine, MultiFab& S_crse, 
                       int scomp, int ncomp, const IntVect& ratio)
    {
        BL_ASSERT(S_crse.nComp() == S_fine.nComp());
        BL_ASSERT((S_crse.is_cell_centered() && S_fine.is_cell_centered()) ||
                  (S_crse.is_nodal()         && S_fine.is_nodal()));

        bool is_cell_centered = S_crse.is_cell_centered();
        
        //
        // Coarsen() the fine stuff on processors owning the fine data.
        //
        BoxArray crse_S_fine_BA = S_fine.boxArray(); crse_S_fine_BA.coarsen(ratio);
        
        if (crse_S_fine_BA == S_crse.boxArray() and S_fine.DistributionMap() == S_crse.DistributionMap())
        {
#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi(S_crse,true); mfi.isValid(); ++mfi)
            {
                //  NOTE: The tilebox is defined at the coarse level.
                const Box& tbx = mfi.tilebox();

                if (is_cell_centered) {
                    BL_FORT_PROC_CALL(BL_AVGDOWN,bl_avgdown)
                        (tbx.loVect(), tbx.hiVect(),
                         BL_TO_FORTRAN_N(S_fine[mfi],scomp),
                         BL_TO_FORTRAN_N(S_crse[mfi],scomp),
                         ratio.getVect(),&ncomp);
                } else {
                    BL_FORT_PROC_CALL(BL_AVGDOWN_NODES,bl_avgdown_nodes)
                        (tbx.loVect(),tbx.hiVect(),
                         BL_TO_FORTRAN_N(S_fine[mfi],scomp),
                         BL_TO_FORTRAN_N(S_crse[mfi],scomp),
                         ratio.getVect(),&ncomp);
                }
            }
        }
        else
        {
            MultiFab crse_S_fine(crse_S_fine_BA, S_fine.DistributionMap(), ncomp, 0, MFInfo(), FArrayBoxFactory());

#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi(crse_S_fine,true); mfi.isValid(); ++mfi)
            {
                //  NOTE: The tilebox is defined at the coarse level.
                const Box& tbx = mfi.tilebox();
                
                //  NOTE: We copy from component scomp of the fine fab into component 0 of the crse fab
                //        because the crse fab is a temporary which was made starting at comp 0, it is
                //        not part of the actual crse multifab which came in.

                if (is_cell_centered) {
                    BL_FORT_PROC_CALL(BL_AVGDOWN,bl_avgdown)
                        (tbx.loVect(), tbx.hiVect(),
                         BL_TO_FORTRAN_N(S_fine[mfi],scomp),
                         BL_TO_FORTRAN_N(crse_S_fine[mfi],0),
                         ratio.getVect(),&ncomp);
                } else {
                    BL_FORT_PROC_CALL(BL_AVGDOWN_NODES,bl_avgdown_nodes)
                        (tbx.loVect(), tbx.hiVect(),
                         BL_TO_FORTRAN_N(S_fine[mfi],scomp),
                         BL_TO_FORTRAN_N(crse_S_fine[mfi],0),
                         ratio.getVect(),&ncomp);
                }
            }
            
            S_crse.copy(crse_S_fine,0,scomp,ncomp);
        }
   }

// *************************************************************************************************************

    // Average fine face-based MultiFab onto crse face-based MultiFab.
    void average_down_faces (const Vector<const MultiFab*>& fine, const Vector<MultiFab*>& crse,
			     const IntVect& ratio, int ngcrse)
    {
	BL_ASSERT(crse.size()  == AMREX_SPACEDIM);
	BL_ASSERT(fine.size()  == AMREX_SPACEDIM);
	BL_ASSERT(crse[0]->nComp() == fine[0]->nComp());

	int ncomp = crse[0]->nComp();

        if (isMFIterSafe(*fine[0], *crse[0]))
        {
#ifdef _OPENMP
#pragma omp parallel
#endif
            for (int n=0; n<AMREX_SPACEDIM; ++n) {
                for (MFIter mfi(*crse[n],true); mfi.isValid(); ++mfi)
                {
                    const Box& tbx = mfi.growntilebox(ngcrse);
                    
                    BL_FORT_PROC_CALL(BL_AVGDOWN_FACES,bl_avgdown_faces)
                        (tbx.loVect(),tbx.hiVect(),
                         BL_TO_FORTRAN((*fine[n])[mfi]),
                         BL_TO_FORTRAN((*crse[n])[mfi]),
                         ratio.getVect(),n,ncomp);
                }
            }
        }
        else
        {
            std::array<MultiFab,AMREX_SPACEDIM> ctmp;
            Vector<MultiFab*> vctmp(AMREX_SPACEDIM);
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
            {
                BoxArray cba = fine[idim]->boxArray();
                cba.coarsen(ratio);
                ctmp[idim].define(cba, fine[idim]->DistributionMap(), ncomp, ngcrse, MFInfo(), FArrayBoxFactory());
                vctmp[idim] = &ctmp[idim];
            }
            average_down_faces(fine, vctmp, ratio, ngcrse);
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
            {
                crse[idim]->ParallelCopy(ctmp[idim],0,0,ncomp,ngcrse,ngcrse);
            }
        }
    }

    //! Average fine edge-based MultiFab onto crse edge-based MultiFab.
    //! This routine assumes that the crse BoxArray is a coarsened version of the fine BoxArray.
    void average_down_edges (const Vector<const MultiFab*>& fine, const Vector<MultiFab*>& crse,
                             const IntVect& ratio, int ngcrse)
    {
	BL_ASSERT(crse.size()  == AMREX_SPACEDIM);
	BL_ASSERT(fine.size()  == AMREX_SPACEDIM);
	BL_ASSERT(crse[0]->nComp() == fine[0]->nComp());

	int ncomp = crse[0]->nComp();

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (int n=0; n<AMREX_SPACEDIM; ++n) {
            for (MFIter mfi(*crse[n],true); mfi.isValid(); ++mfi)
            {
                const Box& tbx = mfi.growntilebox(ngcrse);

                BL_FORT_PROC_CALL(BL_AVGDOWN_EDGES,bl_avgdown_edges)
                    (tbx.loVect(),tbx.hiVect(),
                     BL_TO_FORTRAN((*fine[n])[mfi]),
                     BL_TO_FORTRAN((*crse[n])[mfi]),
                     ratio.getVect(),n,ncomp);
            }
        }
    }

    //! Average fine node-based MultiFab onto crse node-based MultiFab.
    //! This routine assumes that the crse BoxArray is a coarsened version of the fine BoxArray.
    void average_down_nodal (const MultiFab& fine, MultiFab& crse, const IntVect& ratio, int ngcrse)
    {
        BL_ASSERT(fine.is_nodal());
        BL_ASSERT(crse.is_nodal());
	BL_ASSERT(crse.nComp() == fine.nComp());

	int ncomp = crse.nComp();

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(crse,true); mfi.isValid(); ++mfi)
            {
                const Box& tbx = mfi.growntilebox(ngcrse);
                
                BL_FORT_PROC_CALL(BL_AVGDOWN_NODES,bl_avgdown_nodes)
                    (tbx.loVect(),tbx.hiVect(),
                     BL_TO_FORTRAN(fine[mfi]),
                     BL_TO_FORTRAN(crse[mfi]),
                     ratio.getVect(),&ncomp);
            }
    }

// *************************************************************************************************************

    void fill_boundary(MultiFab& mf, const Geometry& geom, bool cross)
    {
	amrex::fill_boundary(mf, 0, mf.nComp(), geom, cross);
    }

// *************************************************************************************************************

    void fill_boundary(MultiFab& mf, int scomp, int ncomp, const Geometry& geom, bool cross)
    {
	mf.FillBoundary(scomp,ncomp,geom.periodicity(),cross);
    }


    void print_state(const MultiFab& mf, const IntVect& cell, const int n)
    {


#ifdef _OPENMP
#pragma omp parallel
#endif
      for (MFIter mfi(mf); mfi.isValid(); ++mfi)
	{
	  if (mfi.validbox().contains(cell)) {
	    if (n >= 0) {
	      amrex::AllPrint().SetPrecision(17) << mf[mfi](cell, n) << std::endl;
	    } else {
	      std::ostringstream ss;
	      ss.precision(17);
	      const int ncomp = mf.nComp();
	      for (int i = 0; i < ncomp-1; ++i)
		{
		  ss << mf[mfi](cell,i) << ", ";
		}
	      ss << mf[mfi](cell,ncomp-1);
	      amrex::AllPrint() << ss.str() << std::endl;	    
	    }
	  }
	}
      
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
        MultiFab mf(imf.boxArray(), imf.DistributionMap(), imf.nComp(), imf.nGrow(),
                    MFInfo(), FArrayBoxFactory());

        const int ncomp = imf.nComp();
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(mf,true); mfi.isValid(); ++mfi)
        {
            const Box& tbx = mfi.growntilebox();
            amrex_fort_int_to_real(BL_TO_FORTRAN_BOX(tbx), &ncomp,
                                   BL_TO_FORTRAN_ANYD(mf[mfi]),
                                   BL_TO_FORTRAN_ANYD(imf[mfi]));
        }

        return mf;
    }

    std::unique_ptr<MultiFab> get_slice_data(int dir, Real coord, const MultiFab& cc, const Geometry& geom, int start_comp, int ncomp) {

        BL_PROFILE("amrex::get_slice_data");
        
        Vector<int> slice_to_full_ba_map;
        std::unique_ptr<MultiFab> slice = allocateSlice(dir, cc, ncomp, geom, coord, slice_to_full_ba_map);
        
        int nf = cc.nComp();        
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(*slice, true); mfi.isValid(); ++mfi) {
            int slice_gid = mfi.index();
            int full_gid = slice_to_full_ba_map[slice_gid];
        
            const Box& tile_box  = mfi.tilebox();
            amrex_fill_slice(BL_TO_FORTRAN_BOX(tile_box),
                             BL_TO_FORTRAN_ANYD(cc[full_gid]),
                             BL_TO_FORTRAN_ANYD((*slice)[slice_gid]),
                             &start_comp, &nf, &ncomp);
        }
        
        return slice;
    }
}
