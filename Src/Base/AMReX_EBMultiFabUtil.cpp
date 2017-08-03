
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_EBFabFactory.H>

namespace amrex
{

void
EB_set_covered (MultiFab& mf)
{
    EB_set_covered(mf, 0, mf.nComp());
}

void
EB_set_covered (MultiFab& mf, int icomp, int ncomp)
{
    Array<Real> minvals(ncomp);
    for (int i = icomp; i < icomp+ncomp; ++i) {
        minvals[i] = mf.min(i,0,true);
    }
    ParallelDescriptor::ReduceRealMin(minvals.data(), ncomp);

    const auto& factory = dynamic_cast<EBFArrayBoxFactory const&>(mf.Factory());
    const auto& layout = factory.getEBISLayout();

    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        FArrayBox& fab = mf[mfi];
        const EBISBox& ebis = layout[mfi];
        for (BoxIterator bi(mfi.validbox()); bi.ok(); ++bi)
        {
            const IntVect& iv = bi();
            if (ebis.isCovered(iv))
            {
                for (int i = icomp; i < icomp+ncomp; ++i) {
                    fab(iv, i) = minvals[i-icomp];
                }
            }
        }
    }
}

}
