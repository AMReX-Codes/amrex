#include <MultiFab_C_F.H>

namespace amrex {

int  MultiFab_C_to_F::count = 0;

MultiFab_C_to_F::MultiFab_C_to_F (const Geometry& geom,
				  const DistributionMapping& dmap,
				  const BoxArray& ba)
	  
{
    BL_ASSERT(count == 0);
    count++;

    int nb = ba.size();
    int dm = BL_SPACEDIM;

    std::vector<int> lo(nb*dm);
    std::vector<int> hi(nb*dm);

    for ( int i = 0; i < nb; ++i ) {
	const Box& bx = amrex::enclosedCells(ba[i]);
        for ( int j = 0; j < dm; ++j ) {
	    lo[j + i*dm] = bx.smallEnd(j);
	    hi[j + i*dm] = bx.bigEnd(j);
	}
    }

    const Box& domain = geom.Domain();

    int pm[dm];
    for ( int i = 0; i < dm; ++i ) {
	pm[i] = geom.isPeriodic(i)? 1 : 0;
    }

    const Vector<int>& pmap = dmap.ProcessorMap();

    build_layout_from_c(nb, dm, &lo[0], &hi[0], 
			domain.loVect(), domain.hiVect(), 
			pm, pmap.dataPtr());
}

MultiFab_C_to_F::~MultiFab_C_to_F ()
{
    count--;
    destroy_multifab_c();
}

void 
MultiFab_C_to_F::share (MultiFab& cmf, const std::string& fmf_name)
{
    const Box& bx = cmf.boxArray()[0];
    int nodal[BL_SPACEDIM];
    for ( int i = 0; i < BL_SPACEDIM; ++i ) {
	nodal[i] = (bx.type(i) == IndexType::NODE) ? 1 : 0;
    }

    share_multifab_with_f (fmf_name.c_str(), cmf.nComp(), cmf.nGrow(), nodal);

    for (MFIter mfi(cmf); mfi.isValid(); ++mfi)
    {
	int li = mfi.LocalIndex();
	const FArrayBox& fab = cmf[mfi];
	share_fab_with_f (li, fab.dataPtr());
    }
}

}
