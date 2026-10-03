#ifndef AMREX_MF_PARALLEL_FOR_C_H_
#define AMREX_MF_PARALLEL_FOR_C_H_
#include <AMReX_Config.H>

#ifndef AMREX_USE_GPU

#include <AMReX_Concepts.H>
#include <AMReX_MFIter.H>

/// \cond DOXYGEN_IGNORE
namespace amrex::detail {

template <FabArrayType MF, typename F>
void ParallelFor_doit (MF const& mf, IntVect const& nghost, IntVect const& ts, bool dynamic, F const& f)
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

template <FabArrayType MF, typename F>
void ParallelFor_doit (MF const& mf, IntVect const& nghost, int ncomp, IntVect const& ts, bool dynamic, F const& f)
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

template <FabArrayType MF, typename F>
void ParallelForStrided_doit (MF const& mf, IntVect const& stride, IntVect const& offset, F const& f)
{
    const Dim3 s = stride.dim3(1);
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
    for (MFIter mfi(mf,MFItInfo().EnableTiling(FabArrayBase::mfiter_tile_size)); mfi.isValid(); ++mfi) {
        Box const& bx = mfi.tilebox();
        int const lidx = mfi.LocalIndex();
        const auto lo = amrex::lbound(strided_box(bx, stride, offset));
        const auto hi = amrex::ubound(bx);
        for (        int k = lo.z; k <= hi.z; k += s.z) {
            for (    int j = lo.y; j <= hi.y; j += s.y) {
                for (int i = lo.x; i <= hi.x; i += s.x) {
                    f(lidx,i,j,k);
                }
            }
        }
    }
}

template <FabArrayType MF, typename F>
void ParallelForRedBlack_doit (MF const& mf, int redblack, F const& f)
{
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
    for (MFIter mfi(mf,MFItInfo().EnableTiling(FabArrayBase::mfiter_tile_size)); mfi.isValid(); ++mfi) {
        Box const& bx = mfi.tilebox();
        int const lidx = mfi.LocalIndex();
        const auto lo = amrex::lbound(bx);
        const auto hi = amrex::ubound(bx);
        for (        int k = lo.z; k <= hi.z; ++k) {
            for (    int j = lo.y; j <= hi.y; ++j) {
                for (int i = lo.x + ((lo.x+j+k+redblack) & 1); i <= hi.x; i += 2) {
                    f(lidx,i,j,k);
                }
            }
        }
    }
}

}
/// \endcond

#endif
#endif
