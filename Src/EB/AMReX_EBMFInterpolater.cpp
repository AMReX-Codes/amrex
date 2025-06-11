#include <AMReX_EBMFInterpolater.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

namespace amrex {

EBMFCellConsLinInterp eb_mf_cell_cons_interp(false);
EBMFCellConsLinInterp eb_mf_lincc_interp(true);

void
EBMFCellConsLinInterp::interp (MultiFab const& crsemf, int ccomp, MultiFab& finemf, int fcomp, int nc,
                               IntVect const& ng, Geometry const& cgeom, Geometry const& fgeom,
                               Box const& dest_domain, IntVect const& ratio,
                               Vector<BCRec> const& bcs, int bcomp)
{
    MFCellConsLinInterp::interp(crsemf, ccomp, finemf, fcomp, nc, ng, cgeom, fgeom, dest_domain,
                                ratio, bcs, bcomp);

    const auto *const ffact = dynamic_cast<EBFArrayBoxFactory const*>(&finemf.Factory());
    const auto *const cfact = dynamic_cast<EBFArrayBoxFactory const*>(&crsemf.Factory());
    if (!ffact || !cfact || ffact->isAllRegular()) {
        return;
    }

    auto const& cflags = cfact->getMultiEBCellFlagFab();

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion()) {
        auto const& crse = crsemf.const_arrays();
        auto const& fine = finemf.arrays();
        auto const& cflg = cflags.const_arrays();
        ParallelFor(finemf, ng, nc,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
        {
            if (dest_domain.contains(i,j,k)) {
                int ic = amrex::coarsen(i,ratio[0]);
                int jc = amrex::coarsen(j,ratio[1]);
#if (AMREX_SPACEDIM == 2)
                int kc = k;
#else
                int kc = amrex::coarsen(k,ratio[2]);
#endif
                if (cflg[box_no](ic,jc,kc).numNeighbors() < AMREX_D_TERM(3,*3,*3)) {
                    fine[box_no](i,j,k,n+fcomp) = crse[box_no](ic,jc,kc,n+ccomp);
                }
            }
        });
        Gpu::streamSynchronize();
    } else
#endif
    {
        auto const& fflags = ffact->getMultiEBCellFlagFab();
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
        for (MFIter mfi(finemf); mfi.isValid(); ++mfi) {
            Box const target_fine_region = amrex::grow(mfi.validbox(),ng) & dest_domain;
            Box const& crse_bx = CoarseBox(target_fine_region, ratio);
            const EBCellFlagFab& crse_flag_fab = cflags[mfi];
            if (crsemf[mfi].getType() != FabType::regular) {
                const EBCellFlagFab& fine_flag_fab = fflags[mfi];
                const FabType ftype = fine_flag_fab.getType(target_fine_region);
                const FabType ctype = crse_flag_fab.getType(crse_bx);
                if (ftype == FabType::multivalued || ctype == FabType::multivalued)
                {
                    amrex::Abort("EBCellConservativeLinear::interp: multivalued not implemented");
                }
                else if (ftype == FabType::covered)
                {
                    ; // don't need to do anything special
                }
                else
                {
                    const auto& crse_flag = cflags.const_array(mfi);
                    const auto& fine = finemf.array(mfi);
                    const auto& crse = crsemf.const_array(mfi);
                    amrex::LoopConcurrentOnCpu(target_fine_region, nc,
                    [&] (int i, int j, int k, int n) noexcept
                    {
                        int ic = amrex::coarsen(i,ratio[0]);
                        int jc = amrex::coarsen(j,ratio[1]);
#if (AMREX_SPACEDIM == 2)
                        int kc = k;
#else
                        int kc = amrex::coarsen(k,ratio[2]);
#endif
                        if (crse_flag(ic,jc,kc).numNeighbors() < AMREX_D_TERM(3,*3,*3)) {
                            fine(i,j,k,n+fcomp) = crse(ic,jc,kc,n+ccomp);
                        }
                    });
                }
            }
        }
    }
}

}
