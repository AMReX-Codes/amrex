
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

    void fill_boundary(MultiFab& mf, const Geometry& geom, bool cross)
    {
	BoxLib::fill_boundary(mf, 0, mf.nComp(), geom, cross);
    }

    void fill_boundary(MultiFab& mf, int scomp, int ncomp, const Geometry& geom, bool cross)
    {
	if (mf.nGrow() <= 0) return;
	
	bool local = false;  // Don't think we ever want it to be true.
	mf.FillBoundary(scomp, ncomp, local, cross);

	bool do_corners = !cross;
	geom.FillPeriodicBoundary(mf, scomp, ncomp, do_corners, local);
    }
}
