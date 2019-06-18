
#include <AMReX_Extrapolater.H>
#include <AMReX_iMultiFab.H>

#ifdef _OPENMP
#include <omp.h>
#endif

extern "C"
{
    void amrex_first_order_extrap(amrex::Real* u, const int* ulo, const int* uhi, const int& nu,
			    const int* msk, const int* mlo, const int* mhi,
			    const int* lo, const int* hi, 
			    const int& scomp, const int& ncomp);
}

namespace amrex {

namespace Extrapolater
{
    void FirstOrderExtrap (MultiFab& mf, const Geometry& geom, int scomp, int ncomp)
    {
        // TODO gpu

	BL_ASSERT(mf.nGrow() == 1);
	BL_ASSERT(scomp >= 0);
	BL_ASSERT(ncomp <= mf.nComp());

	iMultiFab mask(mf.boxArray(), mf.DistributionMap(), 1, 1, MFInfo(),
                       DefaultFabFactory<IArrayBox>());
	mask.BuildMask(geom.Domain(), geom.periodicity(),
		       finebnd, crsebnd, physbnd, interior);

	int N = mf.nComp();

#ifdef _OPENMP
#pragma omp parallel
#endif
	for (MFIter mfi(mf); mfi.isValid(); ++mfi)
	{
	    const Box& bx = mfi.validbox();
	    const IArrayBox& maskfab = mask[mfi];
	    const Box& maskbox = maskfab.box();
	    FArrayBox& datafab = mf[mfi];
	    const Box& databox = datafab.box();

	    amrex_first_order_extrap(datafab.dataPtr(), databox.loVect(), databox.hiVect(), N, 
			       maskfab.dataPtr(), maskbox.loVect(), maskbox.hiVect(),
			       bx.loVect(), bx.hiVect(), scomp, ncomp);
	}
    }
}

}
