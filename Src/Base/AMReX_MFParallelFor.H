#ifndef AMREX_MF_PARALLEL_FOR_H_
#define AMREX_MF_PARALLEL_FOR_H_
#include <AMReX_Config.H>

#include <AMReX_FabArrayBase.H>
#include <AMReX_TypeTraits.H>

#ifdef AMREX_USE_GPU
#include <AMReX_MFParallelForG.H>
#else
#include <AMReX_MFParallelForC.H>
#endif

namespace amrex {
namespace experimental {

struct ParForMFTileSize {
    IntVect tile_size;
    ParForMFTileSize (int sx, int sy, int sz)
        : tile_size(AMREX_D_DECL(sx,sy,sz))
        { amrex::ignore_unused(sy,sz); }
};

/**
 * \brief ParallelFor for MultiFab/FabArray.
 *
 * This version launch a kernel to work on the valid region.  If built for
 * CPU, tiling will be enabled.  For GPU build, this funciton is
 * NON-BLOCKING on the host. Conceptually, this is a 4D loop.
 *
 * \tparam MF the MultiFab/FabArray type
 * \tparam F a callable type like lambda
 *
 * \param mf the MultiFab/FabArray object used to specify the iteration space
 * \param f a calable object void(int,int,int,int), where the first argument
 *           is the local box index, and the following three are spatial indices
 *           for x, y, and z-directions.
 */
template <typename MF, typename F>
std::enable_if_t<IsFabArray<MF>::value>
ParallelFor (MF const& mf, F&& f)
{
    detail::ParallelFor(mf, IntVect(0), FabArrayBase::mfiter_tile_size, std::forward<F>(f));
}

/**
 * \brief ParallelFor for MultiFab/FabArray.
 *
 * This version launch a kernel to work on the valid and ghost region.  If
 * built for CPU, tiling will be enabled.  For GPU build, this funciton is
 * NON-BLOCKING on the host. Conceptually, this is a 4D loop.
 *
 * \tparam MF the MultiFab/FabArray type
 * \tparam F a callable type like lambda
 *
 * \param mf the MultiFab/FabArray object used to specify the iteration space
 * \param ng the number of ghost cells around the valid region
 * \param f a calable object void(int,int,int,int), where the first argument
 *           is the local box index, and the following three are spatial indices
 *           for x, y, and z-directions.
 */
template <typename MF, typename F>
std::enable_if_t<IsFabArray<MF>::value>
ParallelFor (MF const& mf, IntVect const& ng, F&& f)
{
    detail::ParallelFor(mf, ng, FabArrayBase::mfiter_tile_size, std::forward<F>(f));
}

/**
 * \brief ParallelFor for MultiFab/FabArray.
 *
 * This version launch a kernel to work on the valid region.  If built for
 * CPU, tiling will be enabled.  For GPU build, this funciton is
 * NON-BLOCKING on the host. Conceptually, this is a 5D loop.
 *
 * \tparam MF the MultiFab/FabArray type
 * \tparam F a callable type like lambda
 *
 * \param mf the MultiFab/FabArray object used to specify the iteration space
 * \param ncomp the number of component
 * \param f a calable object void(int,int,int,int,int), where the first argument
 *           is the local box index, the following three are spatial indices
 *           for x, y, and z-directions, and the last is for component.
 */
template <typename MF, typename F>
std::enable_if_t<IsFabArray<MF>::value>
ParallelFor (MF const& mf, int ncomp, F&& f)
{
    detail::ParallelFor(mf, IntVect(0), ncomp, FabArrayBase::mfiter_tile_size, std::forward<F>(f));
}

/**
 * \brief ParallelFor for MultiFab/FabArray.
 *
 * This version launch a kernel to work on the valid and ghost region.  If
 * built for CPU, tiling will be enabled.  For GPU build, this funciton is
 * NON-BLOCKING on the host. Conceptually, this is a 5D loop.
 *
 * \tparam MF the MultiFab/FabArray type
 * \tparam F a callable type like lambda
 *
 * \param mf the MultiFab/FabArray object used to specify the iteration space
 * \param ng the number of ghost cells around the valid region
 * \param ncomp the number of component
 * \param f a calable object void(int,int,int,int,int), where the first argument
 *           is the local box index, the following three are spatial indices
 *           for x, y, and z-directions, and the last is for component.
 */
template <typename MF, typename F>
std::enable_if_t<IsFabArray<MF>::value>
ParallelFor (MF const& mf, IntVect const& ng, int ncomp, F&& f)
{
    detail::ParallelFor(mf, ng, ncomp, FabArrayBase::mfiter_tile_size, std::forward<F>(f));
}

/**
 * \brief ParallelFor for MultiFab/FabArray.
 *
 * This version launch a kernel to work on the valid region.  If built for
 * CPU, tiling will be enabled.  However, one could specify a huge tile size
 * to effectively disable tiling.  For GPU build, this funciton is
 * NON-BLOCKING on the host. Conceptually, this is a 4D loop.
 *
 * \tparam MF the MultiFab/FabArray type
 * \tparam F a callable type like lambda
 *
 * \param mf the MultiFab/FabArray object used to specify the iteration space
 * \param ts tile size, ignored by GPU build
 * \param f a calable object void(int,int,int,int), where the first argument
 *           is the local box index, and the following three are spatial indices
 *           for x, y, and z-directions.
 */
template <typename MF, typename F>
std::enable_if_t<IsFabArray<MF>::value>
ParallelFor (MF const& mf, ParForMFTileSize const& ts, F&& f)
{
    detail::ParallelFor(mf, IntVect(0), ts.tile_size, std::forward<F>(f));
}

/**
 * \brief ParallelFor for MultiFab/FabArray.
 *
 * This version launch a kernel to work on the valid and ghost region.  If
 * built for CPU, tiling will be enabled.  However, one could specify a huge
 * tile size to effectively disable tiling.  For GPU build, this funciton is
 * NON-BLOCKING on the host. Conceptually, this is a 4D loop.
 *
 * \tparam MF the MultiFab/FabArray type
 * \tparam F a callable type like lambda
 *
 * \param mf the MultiFab/FabArray object used to specify the iteration space
 * \param ng the number of ghost cells around the valid region
 * \param ts tile size, ignored by GPU build
 * \param f a calable object void(int,int,int,int), where the first argument
 *           is the local box index, and the following three are spatial indices
 *           for x, y, and z-directions.
 */
template <typename MF, typename F>
std::enable_if_t<IsFabArray<MF>::value>
ParallelFor (MF const& mf, IntVect const& ng, ParForMFTileSize const& ts, F&& f)
{
    detail::ParallelFor(mf, ng, ts.tile_size, std::forward<F>(f));
}

/**
 * \brief ParallelFor for MultiFab/FabArray.
 *
 * This version launch a kernel to work on the valid region.  If built for
 * CPU, tiling will be enabled.  However, one could specify a huge tile size
 * to effectively disable tiling.  For GPU build, this funciton is
 * NON-BLOCKING on the host. Conceptually, this is a 5D loop.
 *
 * \tparam MF the MultiFab/FabArray type
 * \tparam F a callable type like lambda
 *
 * \param mf the MultiFab/FabArray object used to specify the iteration space
 * \param ncomp the number of component
 * \param ts tile size, ignored by GPU build
 * \param f a calable object void(int,int,int,int,int), where the first argument
 *           is the local box index, the following three are spatial indices
 *           for x, y, and z-directions, and the last is for component.
 */
template <typename MF, typename F>
std::enable_if_t<IsFabArray<MF>::value>
ParallelFor (MF const& mf, int ncomp, ParForMFTileSize const& ts, F&& f)
{
    detail::ParallelFor(mf, IntVect(0), ncomp, ts.tile_size, std::forward<F>(f));
}

/**
 * \brief ParallelFor for MultiFab/FabArray.
 *
 * This version launch a kernel to work on the valid and ghost region.  If
 * built for CPU, tiling will be enabled.  However, one could specify a huge
 * tile size to effectively disable tiling.  For GPU build, this funciton is
 * NON-BLOCKING on the host. Conceptually, this is a 5D loop.
 *
 * \tparam MF the MultiFab/FabArray type
 * \tparam F a callable type like lambda
 *
 * \param mf the MultiFab/FabArray object used to specify the iteration space
 * \param ng the number of ghost cells around the valid region
 * \param ncomp the number of component
 * \param ts tile size, ignored by GPU build
 * \param f a calable object void(int,int,int,int,int), where the first argument
 *           is the local box index, the following three are spatial indices
 *           for x, y, and z-directions, and the last is for component.
 */
template <typename MF, typename F>
std::enable_if_t<IsFabArray<MF>::value>
ParallelFor (MF const& mf, IntVect const& ng, int ncomp, ParForMFTileSize const& ts, F&& f)
{
    detail::ParallelFor(mf, ng, ncomp, ts.tile_size, std::forward<F>(f));
}

}
}

#endif
