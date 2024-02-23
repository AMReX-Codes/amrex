#ifndef AMREX_MF_PARALLEL_FOR_C_H_
#define AMREX_MF_PARALLEL_FOR_C_H_
#include <AMReX_Config.H>

#ifndef AMREX_USE_GPU

#include <AMReX_MFIter.H>

namespace amrex::experimental::detail {

template <typename MF, typename F>
std::enable_if_t<IsFabArray<MF>::value>
ParallelFor (MF const& mf, IntVect const& nghost, IntVect const& ts, bool dynamic, F&& f)
{
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
    for (MFIter mfi(mf,MFItInfo().EnableTiling(ts).SetDynamic(dynamic)); mfi.isValid(); ++mfi) {
        Box const& bx = mfi.growntilebox(nghost);
        int const lidx = mfi.LocalIndex();
        const auto lo = amrex::lbound(bx);
        const auto hi = amrex::ubound(bx);
        for (        int k = lo.z; k <= hi.z; ++k) {
            for (    int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    f(lidx,i,j,k);
                }
            }
        }
    }
}

template <typename MF, typename F>
std::enable_if_t<IsFabArray<MF>::value>
ParallelFor (MF const& mf, IntVect const& nghost, int ncomp, IntVect const& ts, bool dynamic, F&& f)
{
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
    for (MFIter mfi(mf,MFItInfo().EnableTiling(ts).SetDynamic(dynamic)); mfi.isValid(); ++mfi) {
        Box const& bx = mfi.growntilebox(nghost);
        int const lidx = mfi.LocalIndex();
        const auto lo = amrex::lbound(bx);
        const auto hi = amrex::ubound(bx);
        for (int n = 0; n < ncomp; ++n) {
            for (        int k = lo.z; k <= hi.z; ++k) {
                for (    int j = lo.y; j <= hi.y; ++j) {
                    AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x; ++i) {
                        f(lidx,i,j,k,n);
                    }
                }
            }
        }
    }
}

}

#endif
#endif
