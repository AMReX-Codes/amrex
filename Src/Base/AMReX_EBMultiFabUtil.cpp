
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_EBMultiFabUtil_F.H>

#ifdef _OPENMP
#include <omp.h>
#endif

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
    const auto& eblevel = factory.getEBLevel();
    const auto& flags = eblevel.Flags();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(mf,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        FArrayBox& fab = mf[mfi];
        const auto& flagfab = flags[mfi];
        amrex_eb_set_covered(BL_TO_FORTRAN_BOX(bx),
                             BL_TO_FORTRAN_N_ANYD(fab,icomp),
                             BL_TO_FORTRAN_ANYD(flagfab),
                             minvals.data(),&ncomp);
    }
}

void
EB_set_volume_fraction (MultiFab& mf)
{
    BL_ASSERT(mf.nComp() == 1);

    const auto& factory = dynamic_cast<EBFArrayBoxFactory const&>(mf.Factory());
    const auto& ebisl   = factory.getEBISLayout();
    const auto& eblevel = factory.getEBLevel();
    const auto& flags = eblevel.Flags();
    const Box& domain = eblevel.getDomain();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(mf,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox();
        FArrayBox& fab = mf[mfi];
        fab.setVal(1.0, bx);
        FabType typ = flags[mfi].getType();
        if (typ != FabType::regular)
        {
            if (typ == FabType::covered) {
                fab.setVal(0.0, bx);
            }
            else
            {
                const auto& ebisbox = ebisl[mfi];
                const Box& bx_sect = bx & domain;
                for (BoxIterator bi(bx_sect); bi.ok(); ++bi)
                {
                    const IntVect& iv = bi();
                    const auto& vofs = ebisbox.getVoFs(iv);
                    Real vtot = 0.0;
                    for (const auto& vi : vofs)
                    {
                        vtot += ebisbox.volFrac(vi);
                    }
                    fab(iv) = vtot;
                }
            }
        }
    }
}

}
