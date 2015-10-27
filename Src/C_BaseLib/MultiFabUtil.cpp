
#include <MultiFabUtil.H>
#include <MultiFabUtil_F.H>

namespace BoxLib
{
    void average_face_to_cellcenter (MultiFab& cc, const PArray<MultiFab>& fc, const Geometry& geom)
    {
	BL_ASSERT(cc.nComp() == BL_SPACEDIM);
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

    // Average fine face-based MultiFab onto crse fine-centered MultiFab.
    // This routine assumes that the crse BoxArray is a coarsened version of the fine BoxArray.
    void average_down_faces (PArray<MultiFab>& fine, PArray<MultiFab>& crse, IntVect& ratio)
    {
	BL_ASSERT(crse.size()  == BL_SPACEDIM);
	BL_ASSERT(fine.size()  == BL_SPACEDIM);
	BL_ASSERT(crse[0].nComp() == 1);
	BL_ASSERT(fine[0].nComp() == 1);

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
                     ratio.getVect(),n);
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
	if (mf.nGrow() <= 0) return;
	
	bool local = false;  // Don't think we ever want it to be true.
	mf.FillBoundary(scomp, ncomp, local, cross);

	bool do_corners = !cross;
	geom.FillPeriodicBoundary(mf, scomp, ncomp, do_corners, local);
    }
}
