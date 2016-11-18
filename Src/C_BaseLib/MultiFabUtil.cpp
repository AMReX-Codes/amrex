
#include <MultiFabUtil.H>
#include <MultiFabUtil_F.H>

namespace BoxLib
{
    void average_edge_to_cellcenter (MultiFab& cc, int dcomp, const std::vector<MultiFab*>& edge)
    {
	BL_ASSERT(cc.nComp() >= dcomp + BL_SPACEDIM);
	BL_ASSERT(edge.size() == BL_SPACEDIM);
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
		 D_DECL(BL_TO_FORTRAN((*edge[0])[mfi]),
			BL_TO_FORTRAN((*edge[1])[mfi]),
			BL_TO_FORTRAN((*edge[2])[mfi])));
	}	
    }

    void average_face_to_cellcenter (MultiFab& cc, int dcomp, const std::vector<MultiFab*>& fc)
    {
	BL_ASSERT(cc.nComp() >= dcomp + BL_SPACEDIM);
	BL_ASSERT(fc.size() == BL_SPACEDIM);
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
		 D_DECL(BL_TO_FORTRAN((*fc[0])[mfi]),
			BL_TO_FORTRAN((*fc[1])[mfi]),
			BL_TO_FORTRAN((*fc[2])[mfi])),
		 dx, problo, coord_type);
	}
    }

    void average_face_to_cellcenter (MultiFab& cc, const PArray<MultiFab>& fc, const Geometry& geom)
    {
	BL_ASSERT(cc.nComp() >= BL_SPACEDIM);
	BL_ASSERT(fc.size() == BL_SPACEDIM);
	BL_ASSERT(fc[0].nComp() == 1); // We only expect fc to have the gradient perpendicular to the face

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
		 D_DECL(BL_TO_FORTRAN(fc[0][mfi]),
			BL_TO_FORTRAN(fc[1][mfi]),
			BL_TO_FORTRAN(fc[2][mfi])),
		 dx, problo, coord_type);
	}
    }

    void average_cellcenter_to_face (PArray<MultiFab>& fc, const MultiFab& cc, const Geometry& geom)
    {
	BL_ASSERT(cc.nComp() == 1);
	BL_ASSERT(cc.nGrow() >= 1);
	BL_ASSERT(fc.size() == BL_SPACEDIM);
	BL_ASSERT(fc[0].nComp() == 1); // We only expect fc to have the gradient perpendicular to the face

	const Real* dx     = geom.CellSize();
	const Real* problo = geom.ProbLo();
	int coord_type = Geometry::Coord();

#ifdef _OPENMP
#pragma omp parallel
#endif
	for (MFIter mfi(cc,true); mfi.isValid(); ++mfi) 
	{
	    const Box& xbx = mfi.nodaltilebox(0);
#if (BL_SPACEDIM > 1)
	    const Box& ybx = mfi.nodaltilebox(1);
#endif
#if (BL_SPACEDIM == 3)
	    const Box& zbx = mfi.nodaltilebox(2);
#endif
	    
	    BL_FORT_PROC_CALL(BL_AVG_CC_TO_FC,bl_avg_cc_to_fc)
		(xbx.loVect(), xbx.hiVect(),
#if (BL_SPACEDIM > 1)
		 ybx.loVect(), ybx.hiVect(),
#endif
#if (BL_SPACEDIM == 3)
		 zbx.loVect(), zbx.hiVect(),
#endif
		 D_DECL(BL_TO_FORTRAN(fc[0][mfi]),
			BL_TO_FORTRAN(fc[1][mfi]),
			BL_TO_FORTRAN(fc[2][mfi])),
		 BL_TO_FORTRAN(cc[mfi]),
		 dx, problo, coord_type);
	}
    }

// *************************************************************************************************************

    // Average fine cell-based MultiFab onto crse cell-centered MultiFab.
    // We do NOT assume that the coarse layout is a coarsened version of the fine layout.
    // This version DOES use volume-weighting.
    void average_down (MultiFab& S_fine, MultiFab& S_crse, const Geometry& fgeom, const Geometry& cgeom, 
                       int scomp, int ncomp, int rr)
    {
         average_down(S_fine,S_crse,fgeom,cgeom,scomp,ncomp,rr*IntVect::TheUnitVector());
    }

    void average_down (MultiFab& S_fine, MultiFab& S_crse, const Geometry& fgeom, const Geometry& cgeom, 
                       int scomp, int ncomp, const IntVect& ratio)
    {
  
        if (S_fine.is_nodal() || S_crse.is_nodal())
        {
            BoxLib::Error("Can't use BoxLib::average_down for nodal MultiFab!");
        }

#if (BL_SPACEDIM == 3)
	BoxLib::average_down(S_fine, S_crse, scomp, ncomp, ratio);
	return;
#else

        BL_ASSERT(S_crse.nComp() == S_fine.nComp());

        //
        // Coarsen() the fine stuff on processors owning the fine data.
        //
	const BoxArray& fine_BA = S_fine.boxArray();
        BoxArray crse_S_fine_BA = fine_BA; 
	crse_S_fine_BA.coarsen(ratio);

        MultiFab crse_S_fine(crse_S_fine_BA,ncomp,0);

	MultiFab fvolume;
	fgeom.GetVolume(fvolume, fine_BA, 0);

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
    void average_down (MultiFab& S_fine, MultiFab& S_crse, int scomp, int ncomp, int rr)
    {
         average_down(S_fine,S_crse,scomp,ncomp,rr*IntVect::TheUnitVector());
    }

    void average_down (MultiFab& S_fine, MultiFab& S_crse, 
                       int scomp, int ncomp, const IntVect& ratio)
    {
        BL_ASSERT(S_crse.nComp() == S_fine.nComp());

        //
        // Coarsen() the fine stuff on processors owning the fine data.
        //
        BoxArray crse_S_fine_BA = S_fine.boxArray(); crse_S_fine_BA.coarsen(ratio);

        MultiFab crse_S_fine(crse_S_fine_BA,ncomp,0);

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

            BL_FORT_PROC_CALL(BL_AVGDOWN,bl_avgdown)
                (tbx.loVect(), tbx.hiVect(),
                 BL_TO_FORTRAN_N(S_fine[mfi],scomp),
                 BL_TO_FORTRAN_N(crse_S_fine[mfi],0),
                 ratio.getVect(),&ncomp);
        }

        S_crse.copy(crse_S_fine,0,scomp,ncomp);
   }

// *************************************************************************************************************

    // Average fine face-based MultiFab onto crse fine-centered MultiFab.
    // This routine assumes that the crse BoxArray is a coarsened version of the fine BoxArray.
    void average_down_faces (PArray<MultiFab>& fine, PArray<MultiFab>& crse, IntVect& ratio)
    {
	BL_ASSERT(crse.size()  == BL_SPACEDIM);
	BL_ASSERT(fine.size()  == BL_SPACEDIM);
	BL_ASSERT(crse[0].nComp() == fine[0].nComp());

	int ncomp = crse[0].nComp();

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (int n=0; n<BL_SPACEDIM; ++n) {
            for (MFIter mfi(crse[n],true); mfi.isValid(); ++mfi)
            {
                const Box& tbx = mfi.tilebox();

                BL_FORT_PROC_CALL(BL_AVGDOWN_FACES,bl_avgdown_faces)
                    (tbx.loVect(),tbx.hiVect(),
                     BL_TO_FORTRAN(fine[n][mfi]),
                     BL_TO_FORTRAN(crse[n][mfi]),
                     ratio.getVect(),n,ncomp);
            }
        }
    }

// *************************************************************************************************************

    void fill_boundary(MultiFab& mf, const Geometry& geom, bool cross)
    {
	BoxLib::fill_boundary(mf, 0, mf.nComp(), geom, cross);
    }

// *************************************************************************************************************

    void fill_boundary(MultiFab& mf, int scomp, int ncomp, const Geometry& geom, bool cross)
    {
	mf.FillBoundary(scomp,ncomp,geom.periodicity(),cross);
    }
}
