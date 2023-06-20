#include <AMReX_Interp_C.H>
#include <AMReX_MFInterp_C.H>
#include <AMReX_MFInterpolater.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

namespace amrex {

// Cell centered
MFPCInterp          mf_pc_interp;
MFCellConsLinInterp mf_cell_cons_interp(false);
MFCellConsLinInterp mf_lincc_interp(true);
MFCellConsLinMinmaxLimitInterp mf_linear_slope_minmax_interp;
MFCellBilinear      mf_cell_bilinear_interp;

// Nodal
MFNodeBilinear      mf_node_bilinear_interp;

Box
MFPCInterp::CoarseBox (const Box& fine, const IntVect& ratio)
{
    return amrex::coarsen(fine,ratio);
}

Box
MFPCInterp::CoarseBox (const Box& fine, int ratio)
{
    return amrex::coarsen(fine,ratio);
}

void
MFPCInterp::interp (MultiFab const& crsemf, int ccomp, MultiFab& finemf, int fcomp, int nc,
                    IntVect const& ng, Geometry const&, Geometry const&,
                    Box const& dest_domain, IntVect const& ratio,
                    Vector<BCRec> const&, int)
{
    AMREX_ASSERT(crsemf.nGrowVect() == 0);

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion()) {
        auto const& crse = crsemf.const_arrays();
        auto const& fine = finemf.arrays();
        ParallelFor(finemf, ng, nc,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
        {
            if (dest_domain.contains(i,j,k)) {
                AMREX_D_TERM(int ic = amrex::coarsen(i,ratio[0]);,
                             int jc = amrex::coarsen(j,ratio[1]);,
                             int kc = amrex::coarsen(k,ratio[2]);)
                AMREX_D_PICK(fine[box_no](i,0,0,n+fcomp) = crse[box_no](ic, 0, 0,n+ccomp);,
                             fine[box_no](i,j,0,n+fcomp) = crse[box_no](ic,jc, 0,n+ccomp);,
                             fine[box_no](i,j,k,n+fcomp) = crse[box_no](ic,jc,kc,n+ccomp);)
            }
        });
        Gpu::streamSynchronize();
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
        for (MFIter mfi(finemf); mfi.isValid(); ++mfi) {
            auto const& fine = finemf.array(mfi);
            auto const& crse = crsemf.const_array(mfi);
            Box const& fbox = amrex::grow(mfi.validbox(), ng) & dest_domain;
            pcinterp_interp(fbox, fine, fcomp, nc, crse, ccomp, ratio);
        }
    }
}

Box
MFCellConsLinInterp::CoarseBox (const Box& fine, const IntVect& ratio)
{
    Box crse = amrex::coarsen(fine,ratio);
    crse.grow(1);
    return crse;
}

Box
MFCellConsLinInterp::CoarseBox (const Box& fine, int ratio)
{
    Box crse = amrex::coarsen(fine,ratio);
    crse.grow(1);
    return crse;
}

void
MFCellConsLinInterp::interp (MultiFab const& crsemf, int ccomp, MultiFab& finemf, int fcomp, int nc,
                             IntVect const& ng, Geometry const& cgeom, Geometry const& fgeom,
                             Box const& dest_domain, IntVect const& ratio,
                             Vector<BCRec> const& bcs, int bcomp)
{
    AMREX_ASSERT(crsemf.nGrowVect() == 0);
    amrex::ignore_unused(fgeom);

    Box const& cdomain = cgeom.Domain();

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion()) {
        MultiFab crse_tmp(crsemf.boxArray(), crsemf.DistributionMap(), AMREX_SPACEDIM*nc, 0);
        auto const& crse = crsemf.const_arrays();
        auto const& tmp = crse_tmp.arrays();
        auto const& ctmp = crse_tmp.const_arrays();
        auto const& fine = finemf.arrays();

        Gpu::DeviceVector<BCRec> d_bc(nc);
        BCRec const* pbc = d_bc.data();
        Gpu::copyAsync(Gpu::hostToDevice, bcs.begin()+bcomp, bcs.begin()+bcomp+nc, d_bc.begin());

#if (AMREX_SPACEDIM == 1)
        if (cgeom.IsSPHERICAL()) {
            Real drf = fgeom.CellSize(0);
            Real rlo = fgeom.Offset(0);
            if (do_linear_limiting) {
                ParallelFor(crsemf, IntVect(-1),
                [=] AMREX_GPU_DEVICE (int box_no, int i, int, int) noexcept
                {
                    mf_cell_cons_lin_interp_llslope(i,0,0, tmp[box_no], crse[box_no], ccomp, nc,
                                                    cdomain, pbc);
                });
            } else {
                ParallelFor(crsemf, IntVect(-1), nc,
                [=] AMREX_GPU_DEVICE (int box_no, int i, int, int, int n) noexcept
                {
                    mf_cell_cons_lin_interp_mcslope_sph(i, n, tmp[box_no], crse[box_no], ccomp, nc,
                                                        cdomain, ratio, pbc, drf, rlo);
                });
            }

            ParallelFor(finemf, ng, nc,
            [=] AMREX_GPU_DEVICE (int box_no, int i, int, int, int n) noexcept
            {
                if (dest_domain.contains(i,0,0)) {
                    mf_cell_cons_lin_interp_sph(i, n, fine[box_no], fcomp, ctmp[box_no],
                                                crse[box_no], ccomp, nc, ratio, drf, rlo);
                }
            });
        } else
#elif (AMREX_SPACEDIM == 2)
        if (cgeom.IsRZ()) {
            Real drf = fgeom.CellSize(0);
            Real rlo = fgeom.Offset(0);
            if (do_linear_limiting) {
                ParallelFor(crsemf, IntVect(-1),
                [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int) noexcept
                {
                    mf_cell_cons_lin_interp_llslope(i,j,0, tmp[box_no], crse[box_no], ccomp, nc,
                                                    cdomain, pbc);
                });
            } else {
                ParallelFor(crsemf, IntVect(-1), nc,
                [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int, int n) noexcept
                {
                    mf_cell_cons_lin_interp_mcslope_rz(i,j,n, tmp[box_no], crse[box_no], ccomp, nc,
                                                       cdomain, ratio, pbc, drf, rlo);
                });
            }

            ParallelFor(finemf, ng, nc,
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int, int n) noexcept
            {
                if (dest_domain.contains(i,j,0)) {
                    mf_cell_cons_lin_interp_rz(i, j, n, fine[box_no], fcomp, ctmp[box_no],
                                               crse[box_no], ccomp, nc, ratio, drf, rlo);
                }
            });
        } else
#endif
        {
            if (do_linear_limiting) {
                ParallelFor(crsemf, IntVect(-1),
                [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
                {
                    mf_cell_cons_lin_interp_llslope(i,j,k, tmp[box_no], crse[box_no], ccomp, nc,
                                                    cdomain, pbc);
                });
            } else {
                ParallelFor(crsemf, IntVect(-1), nc,
                [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
                {
                    mf_cell_cons_lin_interp_mcslope(i,j,k,n, tmp[box_no], crse[box_no], ccomp, nc,
                                                    cdomain, ratio, pbc);
                });
            }

            ParallelFor(finemf, ng, nc,
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
            {
                if (dest_domain.contains(i,j,k)) {
                    mf_cell_cons_lin_interp(i,j,k,n, fine[box_no], fcomp, ctmp[box_no],
                                            crse[box_no], ccomp, nc, ratio);
                }
            });
        }

        Gpu::streamSynchronize();
    } else
#endif
    {
        BCRec const* pbc = bcs.data() + bcomp;

#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
        {
            FArrayBox tmpfab;
            for (MFIter mfi(finemf); mfi.isValid(); ++mfi) {
                auto const& fine = finemf.array(mfi);
                auto const& crse = crsemf.const_array(mfi);

                Box const& cbox = amrex::grow(crsemf[mfi].box(), -1);
                tmpfab.resize(cbox, AMREX_SPACEDIM*nc);
                auto const& tmp = tmpfab.array();
                auto const& ctmp = tmpfab.const_array();

                Box const& fbox = amrex::grow(mfi.validbox(), ng) & dest_domain;

#if (AMREX_SPACEDIM == 1)
                if (cgeom.IsSPHERICAL()) {
                    Real drf = fgeom.CellSize(0);
                    Real rlo = fgeom.Offset(0);
                    if (do_linear_limiting) {
                        amrex::LoopConcurrentOnCpu(cbox,
                        [&] (int i, int, int) noexcept
                        {
                            mf_cell_cons_lin_interp_llslope(i,0,0, tmp, crse, ccomp, nc,
                                                            cdomain, pbc);
                        });
                    } else {
                        amrex::LoopConcurrentOnCpu(cbox, nc,
                        [&] (int i, int, int, int n) noexcept
                        {
                            mf_cell_cons_lin_interp_mcslope_sph(i, n, tmp, crse, ccomp, nc,
                                                                cdomain, ratio, pbc, drf, rlo);
                        });
                    }

                    amrex::LoopConcurrentOnCpu(fbox, nc,
                    [&] (int i, int, int, int n) noexcept
                    {
                        mf_cell_cons_lin_interp_sph(i, n, fine, fcomp, ctmp,
                                                    crse, ccomp, nc, ratio, drf, rlo);
                    });
                } else
#elif (AMREX_SPACEDIM == 2)
                if (cgeom.IsRZ()) {
                    Real drf = fgeom.CellSize(0);
                    Real rlo = fgeom.Offset(0);
                    if (do_linear_limiting) {
                        amrex::LoopConcurrentOnCpu(cbox,
                        [&] (int i, int j, int) noexcept
                        {
                            mf_cell_cons_lin_interp_llslope(i,j,0, tmp, crse, ccomp, nc,
                                                            cdomain, pbc);
                        });
                    } else {
                        amrex::LoopConcurrentOnCpu(cbox, nc,
                        [&] (int i, int j, int, int n) noexcept
                        {
                            mf_cell_cons_lin_interp_mcslope_rz(i, j, n, tmp, crse, ccomp, nc,
                                                               cdomain, ratio, pbc, drf, rlo);
                        });
                    }

                    amrex::LoopConcurrentOnCpu(fbox, nc,
                    [&] (int i, int j, int, int n) noexcept
                    {
                        mf_cell_cons_lin_interp_rz(i, j, n, fine, fcomp, ctmp,
                                                   crse, ccomp, nc, ratio, drf, rlo);
                    });
                } else
#endif
                {
                    if (do_linear_limiting) {
                        amrex::LoopConcurrentOnCpu(cbox,
                        [&] (int i, int j, int k) noexcept
                        {
                            mf_cell_cons_lin_interp_llslope(i,j,k, tmp, crse, ccomp, nc,
                                                            cdomain, pbc);
                        });
                    } else {
                        amrex::LoopConcurrentOnCpu(cbox, nc,
                        [&] (int i, int j, int k, int n) noexcept
                        {
                            mf_cell_cons_lin_interp_mcslope(i,j,k,n, tmp, crse, ccomp, nc,
                                                            cdomain, ratio, pbc);
                        });
                    }

                    amrex::LoopConcurrentOnCpu(fbox, nc,
                    [&] (int i, int j, int k, int n) noexcept
                    {
                        mf_cell_cons_lin_interp(i,j,k,n, fine, fcomp, ctmp,
                                                crse, ccomp, nc, ratio);
                    });
                }
            }
        }
    }
}


Box
MFCellConsLinMinmaxLimitInterp::CoarseBox (const Box& fine, const IntVect& ratio)
{
    Box crse = amrex::coarsen(fine,ratio);
    crse.grow(1);
    return crse;
}

Box
MFCellConsLinMinmaxLimitInterp::CoarseBox (const Box& fine, int ratio)
{
    Box crse = amrex::coarsen(fine,ratio);
    crse.grow(1);
    return crse;
}

void
MFCellConsLinMinmaxLimitInterp::interp (MultiFab const& crsemf, int ccomp, MultiFab& finemf, int fcomp, int nc,
                             IntVect const& ng, Geometry const& cgeom, Geometry const& fgeom,
                             Box const& dest_domain, IntVect const& ratio,
                             Vector<BCRec> const& bcs, int bcomp)
{
    AMREX_ASSERT(crsemf.nGrowVect() == 0);
    amrex::ignore_unused(fgeom);

    Box const& cdomain = cgeom.Domain();

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion()) {
        MultiFab crse_tmp(crsemf.boxArray(), crsemf.DistributionMap(), AMREX_SPACEDIM*nc, 0);
        auto const& crse = crsemf.const_arrays();
        auto const& tmp = crse_tmp.arrays();
        auto const& ctmp = crse_tmp.const_arrays();
        auto const& fine = finemf.arrays();

        Gpu::DeviceVector<BCRec> d_bc(nc);
        BCRec const* pbc = d_bc.data();
        Gpu::copyAsync(Gpu::hostToDevice, bcs.begin()+bcomp, bcs.begin()+bcomp+nc, d_bc.begin());

        ParallelFor(crsemf, IntVect(-1),
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
        {
            mf_cell_cons_lin_interp_limit_minmax_llslope(i,j,k, tmp[box_no], crse[box_no], ccomp, nc,
                                     cdomain, ratio, pbc);
        });

        ParallelFor(finemf, ng, nc,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
        {
            if (dest_domain.contains(i,j,k)) {
                mf_cell_cons_lin_interp(i,j,k,n, fine[box_no], fcomp, ctmp[box_no],
                                        crse[box_no], ccomp, nc, ratio);
            }
        });

        Gpu::streamSynchronize();
    } else
#endif
    {
        BCRec const* pbc = bcs.data() + bcomp;

#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
        {
            FArrayBox tmpfab;
            for (MFIter mfi(finemf); mfi.isValid(); ++mfi) {
                auto const& fine = finemf.array(mfi);
                auto const& crse = crsemf.const_array(mfi);

                Box const& cbox = amrex::grow(crsemf[mfi].box(), -1);
                tmpfab.resize(cbox, AMREX_SPACEDIM*nc);
                auto const& tmp = tmpfab.array();
                auto const& ctmp = tmpfab.const_array();

                Box const& fbox = amrex::grow(mfi.validbox(), ng) & dest_domain;

                amrex::LoopConcurrentOnCpu(cbox,
                [&] (int i, int j, int k) noexcept
                {
                    mf_cell_cons_lin_interp_limit_minmax_llslope(i,j,k, tmp, crse, ccomp, nc,
                                             cdomain, ratio, pbc);
                });

                amrex::LoopConcurrentOnCpu(fbox, nc,
                [&] (int i, int j, int k, int n) noexcept
                {
                    mf_cell_cons_lin_interp(i,j,k,n, fine, fcomp, ctmp,
                                            crse, ccomp, nc, ratio);
                });
            }
        }
    }
}

Box
MFCellBilinear::CoarseBox (const Box& fine, const IntVect& ratio)
{
    const int* lo = fine.loVect();
    const int* hi = fine.hiVect();

    Box crse(amrex::coarsen(fine,ratio));
    const int* clo = crse.loVect();
    const int* chi = crse.hiVect();

    for (int i = 0; i < AMREX_SPACEDIM; i++) {
        if ((lo[i]-clo[i]*ratio[i])*2 < ratio[i]) {
            crse.growLo(i,1);
        }
        if ((hi[i]-chi[i]*ratio[i])*2 >= ratio[i]) {
            crse.growHi(i,1);
        }
    }
    return crse;
}

Box
MFCellBilinear::CoarseBox (const Box& fine, int ratio)
{
    return CoarseBox(fine,IntVect(ratio));
}

void
MFCellBilinear::interp (MultiFab const& crsemf, int ccomp, MultiFab& finemf, int fcomp, int nc,
                        IntVect const& ng, Geometry const&, Geometry const&,
                        Box const& dest_domain, IntVect const& ratio,
                        Vector<BCRec> const&, int)
{
#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion()) {
        auto const& crse = crsemf.const_arrays();
        auto const& fine = finemf.arrays();
        ParallelFor(finemf, ng, nc,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
        {
            if (dest_domain.contains(i,j,k)) {
                mf_cell_bilin_interp(i,j,k,n, fine[box_no], fcomp, crse[box_no], ccomp, ratio);
            }
        });
        Gpu::streamSynchronize();
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
        for (MFIter mfi(finemf); mfi.isValid(); ++mfi) {
            auto const& fine = finemf.array(mfi);
            auto const& crse = crsemf.const_array(mfi);
            Box const& fbox = amrex::grow(mfi.validbox(),ng) & dest_domain;
            amrex::LoopConcurrentOnCpu(fbox, nc,
            [=] (int i, int j, int k, int n) noexcept
            {
                mf_cell_bilin_interp(i,j,k,n, fine, fcomp, crse, ccomp, ratio);
            });
        }
    }
}

Box
MFNodeBilinear::CoarseBox (const Box& fine, const IntVect& ratio)
{
    Box crse = amrex::coarsen(fine,ratio);
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        if (crse.length(i) < 2) {
            crse.growHi(i,1); // Don't want degenerate boxes
        }
    }
    return crse;
}

Box
MFNodeBilinear::CoarseBox (const Box& fine, int ratio)
{
    return CoarseBox(fine, IntVect(ratio));
}

void
MFNodeBilinear::interp (MultiFab const& crsemf, int ccomp, MultiFab& finemf, int fcomp, int nc,
                        IntVect const& ng, Geometry const&, Geometry const&,
                        Box const& dest_domain, IntVect const& ratio,
                        Vector<BCRec> const&, int)
{
#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion()) {
        auto const& crse = crsemf.const_arrays();
        auto const& fine = finemf.arrays();
        ParallelFor(finemf, ng, nc,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
        {
            if (dest_domain.contains(i,j,k)) {
                mf_nodebilin_interp(i,j,k,n, fine[box_no], fcomp, crse[box_no], ccomp, ratio);
            }
        });
        Gpu::streamSynchronize();
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
        for (MFIter mfi(finemf); mfi.isValid(); ++mfi) {
            auto const& fine = finemf.array(mfi);
            auto const& crse = crsemf.const_array(mfi);
            Box const& fbox = amrex::grow(mfi.validbox(), ng) & dest_domain;
            amrex::LoopConcurrentOnCpu(fbox, nc,
            [=] (int i, int j, int k, int n) noexcept
            {
                mf_nodebilin_interp(i,j,k,n, fine, fcomp, crse, ccomp, ratio);
            });
        }
    }
}

}
