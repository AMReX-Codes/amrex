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

struct TileSize {
    IntVect tile_size;
    explicit TileSize (IntVect const& ts) noexcept : tile_size(ts) {}
};

struct DynamicTiling {
    bool dynamic;
    explicit DynamicTiling (bool f) noexcept : dynamic(f) {}
};

namespace experimental {

/**
 * \brief ParallelFor for MultiFab/FabArray.
 *
 * This version launches a kernel to work on the valid region.  If built for
 * CPU, tiling will be enabled.  For GPU builds, this function is
 * NON-BLOCKING on the host. Conceptually, this is a 4D loop.
 *
 * \tparam MF the MultiFab/FabArray type
 * \tparam F a callable type like lambda
 *
 * \param mf the MultiFab/FabArray object used to specify the iteration space
 * \param f a callable object void(int,int,int,int), where the first argument
 *           is the local box index, and the following three are spatial indices
 *           for x, y, and z-directions.
 */
template <typename MF, typename F>
std::enable_if_t<IsFabArray<MF>::value>
ParallelFor (MF const& mf, F&& f)
{
    detail::ParallelFor(mf, IntVect(0), FabArrayBase::mfiter_tile_size, false, std::forward<F>(f));
}

/**
 * \brief ParallelFor for MultiFab/FabArray.
 *
 * This version launches a kernel to work on the valid region.  If built for
 * CPU, tiling will be enabled.  For GPU builds, this function is
 * NON-BLOCKING on the host. Conceptually, this is a 4D loop.
 *
 * \tparam MT max threads in GPU blocks (Only relevant for GPU builds)
 * \tparam MF the MultiFab/FabArray type
 * \tparam F a callable type like lambda
 *
 * \param mf the MultiFab/FabArray object used to specify the iteration space
 * \param f a callable object void(int,int,int,int), where the first argument
 *           is the local box index, and the following three are spatial indices
 *           for x, y, and z-directions.
 */
template <int MT, typename MF, typename F>
std::enable_if_t<IsFabArray<MF>::value>
ParallelFor (MF const& mf, F&& f)
{
#ifdef AMREX_USE_GPU
    detail::ParallelFor<MT>(mf, IntVect(0), FabArrayBase::mfiter_tile_size, false, std::forward<F>(f));
#else
    detail::ParallelFor(mf, IntVect(0), FabArrayBase::mfiter_tile_size, false, std::forward<F>(f));
#endif
}

/**
 * \brief ParallelFor for MultiFab/FabArray.
 *
 * This version launches a kernel to work on the valid and ghost regions.  If
 * built for CPU, tiling will be enabled.  For GPU builds, this function is
 * NON-BLOCKING on the host. Conceptually, this is a 4D loop.
 *
 * \tparam MF the MultiFab/FabArray type
 * \tparam F a callable type like lambda
 *
 * \param mf the MultiFab/FabArray object used to specify the iteration space
 * \param ng the number of ghost cells around the valid region
 * \param f a callable object void(int,int,int,int), where the first argument
 *           is the local box index, and the following three are spatial indices
 *           for x, y, and z-directions.
 */
template <typename MF, typename F>
std::enable_if_t<IsFabArray<MF>::value>
ParallelFor (MF const& mf, IntVect const& ng, F&& f)
{
    detail::ParallelFor(mf, ng, FabArrayBase::mfiter_tile_size, false, std::forward<F>(f));
}

/**
 * \brief ParallelFor for MultiFab/FabArray.
 *
 * This version launches a kernel to work on the valid and ghost regions.  If
 * built for CPU, tiling will be enabled.  For GPU builds, this function is
 * NON-BLOCKING on the host. Conceptually, this is a 4D loop.
 *
 * \tparam MT max threads in GPU blocks (Only relevant for GPU builds)
 * \tparam MF the MultiFab/FabArray type
 * \tparam F a callable type like lambda
 *
 * \param mf the MultiFab/FabArray object used to specify the iteration space
 * \param ng the number of ghost cells around the valid region
 * \param f a callable object void(int,int,int,int), where the first argument
 *           is the local box index, and the following three are spatial indices
 *           for x, y, and z-directions.
 */
template <int MT, typename MF, typename F>
std::enable_if_t<IsFabArray<MF>::value>
ParallelFor (MF const& mf, IntVect const& ng, F&& f)
{
#ifdef AMREX_USE_GPU
    detail::ParallelFor<MT>(mf, ng, FabArrayBase::mfiter_tile_size, false, std::forward<F>(f));
#else
    detail::ParallelFor(mf, ng, FabArrayBase::mfiter_tile_size, false, std::forward<F>(f));
#endif
}

/**
 * \brief ParallelFor for MultiFab/FabArray.
 *
 * This version launches a kernel to work on the valid and ghost regions.  If
 * built for CPU, tiling will be enabled.  For GPU builds, this function is
 * NON-BLOCKING on the host. Conceptually, this is a 5D loop.
 *
 * \tparam MF the MultiFab/FabArray type
 * \tparam F a callable type like lambda
 *
 * \param mf the MultiFab/FabArray object used to specify the iteration space
 * \param ng the number of ghost cells around the valid region
 * \param ncomp the number of component
 * \param f a callable object void(int,int,int,int,int), where the first argument
 *           is the local box index, the following three are spatial indices
 *           for x, y, and z-directions, and the last is for component.
 */
template <typename MF, typename F>
std::enable_if_t<IsFabArray<MF>::value>
ParallelFor (MF const& mf, IntVect const& ng, int ncomp, F&& f)
{
    detail::ParallelFor(mf, ng, ncomp, FabArrayBase::mfiter_tile_size, false, std::forward<F>(f));
}

/**
 * \brief ParallelFor for MultiFab/FabArray.
 *
 * This version launches a kernel to work on the valid and ghost regions.  If
 * built for CPU, tiling will be enabled.  For gpu builds, this function is
 * NON-BLOCKING on the host. Conceptually, this is a 5D loop.
 *
 * \tparam MT max threads in GPU blocks (Only relevant for GPU builds)
 * \tparam MF the MultiFab/FabArray type
 * \tparam F a callable type like lambda
 *
 * \param mf the MultiFab/FabArray object used to specify the iteration space
 * \param ng the number of ghost cells around the valid region
 * \param ncomp the number of component
 * \param f a callable object void(int,int,int,int,int), where the first argument
 *           is the local box index, the following three are spatial indices
 *           for x, y, and z-directions, and the last is for component.
 */
template <int MT, typename MF, typename F>
std::enable_if_t<IsFabArray<MF>::value>
ParallelFor (MF const& mf, IntVect const& ng, int ncomp, F&& f)
{
#ifdef AMREX_USE_GPU
    detail::ParallelFor<MT>(mf, ng, ncomp, FabArrayBase::mfiter_tile_size, false, std::forward<F>(f));
#else
    detail::ParallelFor(mf, ng, ncomp, FabArrayBase::mfiter_tile_size, false, std::forward<F>(f));
#endif
}

/**
 * \brief ParallelFor for MultiFab/FabArray.
 *
 * This version launches a kernel to work on the valid region.  If built for
 * CPU, tiling will be enabled.  However, one could specify a huge tile size
 * to effectively disable tiling.  For gpu builds, this function is
 * NON-BLOCKING on the host. Conceptually, this is a 4D loop.
 *
 * \tparam MF the MultiFab/FabArray type
 * \tparam F a callable type like lambda
 *
 * \param mf the MultiFab/FabArray object used to specify the iteration space
 * \param ts tile size, ignored by GPU build
 * \param f a callable object void(int,int,int,int), where the first argument
 *           is the local box index, and the following three are spatial indices
 *           for x, y, and z-directions.
 */
template <typename MF, typename F>
std::enable_if_t<IsFabArray<MF>::value>
ParallelFor (MF const& mf, TileSize const& ts, F&& f)
{
    detail::ParallelFor(mf, IntVect(0), ts.tile_size, false, std::forward<F>(f));
}

/**
 * \brief ParallelFor for MultiFab/FabArray.
 *
 * This version launches a kernel to work on the valid region.  If built for
 * CPU, tiling will be enabled.  However, one could specify a huge tile size
 * to effectively disable tiling.  For gpu builds, this function is
 * NON-BLOCKING on the host. Conceptually, this is a 4D loop.
 *
 * \tparam MT max threads in GPU blocks (Only relevant for GPU builds)
 * \tparam MF the MultiFab/FabArray type
 * \tparam F a callable type like lambda
 *
 * \param mf the MultiFab/FabArray object used to specify the iteration space
 * \param ts tile size, ignored by GPU build
 * \param f a callable object void(int,int,int,int), where the first argument
 *           is the local box index, and the following three are spatial indices
 *           for x, y, and z-directions.
 */
template <int MT, typename MF, typename F>
std::enable_if_t<IsFabArray<MF>::value>
ParallelFor (MF const& mf, TileSize const& ts, F&& f)
{
#ifdef AMREX_USE_GPU
    detail::ParallelFor<MT>(mf, IntVect(0), ts.tile_size, false, std::forward<F>(f));
#else
    detail::ParallelFor(mf, IntVect(0), ts.tile_size, false, std::forward<F>(f));
#endif
}

/**
 * \brief ParallelFor for MultiFab/FabArray.
 *
 * This version launches a kernel to work on the valid and ghost regions.  If
 * built for CPU, tiling will be enabled.  However, one could specify a huge
 * tile size to effectively disable tiling.  For gpu builds, this function is
 * NON-BLOCKING on the host. Conceptually, this is a 4D loop.
 *
 * \tparam MF the MultiFab/FabArray type
 * \tparam F a callable type like lambda
 *
 * \param mf the MultiFab/FabArray object used to specify the iteration space
 * \param ng the number of ghost cells around the valid region
 * \param ts tile size, ignored by GPU build
 * \param f a callable object void(int,int,int,int), where the first argument
 *           is the local box index, and the following three are spatial indices
 *           for x, y, and z-directions.
 */
template <typename MF, typename F>
std::enable_if_t<IsFabArray<MF>::value>
ParallelFor (MF const& mf, IntVect const& ng, TileSize const& ts, F&& f)
{
    detail::ParallelFor(mf, ng, ts.tile_size, false, std::forward<F>(f));
}

/**
 * \brief ParallelFor for MultiFab/FabArray.
 *
 * This version launches a kernel to work on the valid and ghost regions.  If
 * built for CPU, tiling will be enabled.  However, one could specify a huge
 * tile size to effectively disable tiling.  For gpu builds, this function is
 * NON-BLOCKING on the host. Conceptually, this is a 4D loop.
 *
 * \tparam MT max threads in GPU blocks (Only relevant for GPU builds)
 * \tparam MF the MultiFab/FabArray type
 * \tparam F a callable type like lambda
 *
 * \param mf the MultiFab/FabArray object used to specify the iteration space
 * \param ng the number of ghost cells around the valid region
 * \param ts tile size, ignored by GPU build
 * \param f a callable object void(int,int,int,int), where the first argument
 *           is the local box index, and the following three are spatial indices
 *           for x, y, and z-directions.
 */
template <int MT, typename MF, typename F>
std::enable_if_t<IsFabArray<MF>::value>
ParallelFor (MF const& mf, IntVect const& ng, TileSize const& ts, F&& f)
{
#ifdef AMREX_USE_GPU
    detail::ParallelFor<MT>(mf, ng, ts.tile_size, false, std::forward<F>(f));
#else
    detail::ParallelFor(mf, ng, ts.tile_size, false, std::forward<F>(f));
#endif
}

/**
 * \brief ParallelFor for MultiFab/FabArray.
 *
 * This version launches a kernel to work on the valid and ghost regions.  If
 * built for CPU, tiling will be enabled.  However, one could specify a huge
 * tile size to effectively disable tiling.  For gpu builds, this function is
 * NON-BLOCKING on the host. Conceptually, this is a 5D loop.
 *
 * \tparam MF the MultiFab/FabArray type
 * \tparam F a callable type like lambda
 *
 * \param mf the MultiFab/FabArray object used to specify the iteration space
 * \param ng the number of ghost cells around the valid region
 * \param ncomp the number of component
 * \param ts tile size, ignored by GPU build
 * \param f a callable object void(int,int,int,int,int), where the first argument
 *           is the local box index, the following three are spatial indices
 *           for x, y, and z-directions, and the last is for component.
 */
template <typename MF, typename F>
std::enable_if_t<IsFabArray<MF>::value>
ParallelFor (MF const& mf, IntVect const& ng, int ncomp, TileSize const& ts, F&& f)
{
    detail::ParallelFor(mf, ng, ncomp, ts.tile_size, false, std::forward<F>(f));
}

/**
 * \brief ParallelFor for MultiFab/FabArray.
 *
 * This version launches a kernel to work on the valid and ghost regions.  If
 * built for CPU, tiling will be enabled.  However, one could specify a huge
 * tile size to effectively disable tiling.  For gpu builds, this function is
 * NON-BLOCKING on the host. Conceptually, this is a 5D loop.
 *
 * \tparam MT max threads in GPU blocks (Only relevant for GPU builds)
 * \tparam MF the MultiFab/FabArray type
 * \tparam F a callable type like lambda
 *
 * \param mf the MultiFab/FabArray object used to specify the iteration space
 * \param ng the number of ghost cells around the valid region
 * \param ncomp the number of component
 * \param ts tile size, ignored by GPU build
 * \param f a callable object void(int,int,int,int,int), where the first argument
 *           is the local box index, the following three are spatial indices
 *           for x, y, and z-directions, and the last is for component.
 */
template <int MT, typename MF, typename F>
std::enable_if_t<IsFabArray<MF>::value>
ParallelFor (MF const& mf, IntVect const& ng, int ncomp, TileSize const& ts, F&& f)
{
#ifdef AMREX_USE_GPU
    detail::ParallelFor<MT>(mf, ng, ncomp, ts.tile_size, false, std::forward<F>(f));
#else
    detail::ParallelFor(mf, ng, ncomp, ts.tile_size, false, std::forward<F>(f));
#endif
}

/**
 * \brief ParallelFor for MultiFab/FabArray.
 *
 * This version launches a kernel to work on the valid and ghost regions.  If
 * built for CPU, tiling will be enabled.  However, one could specify a huge
 * tile size to effectively disable tiling.  For GPU builds, this function is
 * NON-BLOCKING on the host. Conceptually, this is a 4D loop.
 *
 * \tparam MF the MultiFab/FabArray type
 * \tparam F a callable type like lambda
 *
 * \param mf the MultiFab/FabArray object used to specify the iteration space
 * \param ng the number of ghost cells around the valid region
 * \param ts tile size, ignored by GPU build
 * \param dt controls dynamic tiling for the cpu build
 * \param f a callable object void(int,int,int,int), where the first argument
 *           is the local box index, and the following three are spatial indices
 *           for x, y, and z-directions.
 */
template <typename MF, typename F>
std::enable_if_t<IsFabArray<MF>::value>
ParallelFor (MF const& mf, IntVect const& ng, TileSize const& ts, DynamicTiling dt, F&& f)
{
    detail::ParallelFor(mf, ng, ts.tile_size, dt.dynamic, std::forward<F>(f));
}

/**
 * \brief ParallelFor for MultiFab/FabArray.
 *
 * This version launches a kernel to work on the valid and ghost regions.  If
 * built for CPU, tiling will be enabled.  However, one could specify a huge
 * tile size to effectively disable tiling.  For gpu builds, this function is
 * NON-BLOCKING on the host. Conceptually, this is a 5D loop.
 *
 * \tparam MF the MultiFab/FabArray type
 * \tparam F a callable type like lambda
 *
 * \param mf the MultiFab/FabArray object used to specify the iteration space
 * \param ng the number of ghost cells around the valid region
 * \param ncomp the number of component
 * \param ts tile size, ignored by GPU build
 * \param dt controls dynamic tiling for the cpu build
 * \param f a callable object void(int,int,int,int,int), where the first argument
 *           is the local box index, the following three are spatial indices
 *           for x, y, and z-directions, and the last is for component.
 */
template <typename MF, typename F>
std::enable_if_t<IsFabArray<MF>::value>
ParallelFor (MF const& mf, IntVect const& ng, int ncomp, TileSize const& ts,
             DynamicTiling dt, F&& f)
{
    detail::ParallelFor(mf, ng, ncomp, ts.tile_size, dt.dynamic, std::forward<F>(f));
}

/**
 * \brief ParallelFor for MultiFab/FabArray.
 *
 * This version launches a kernel to work on the valid and ghost regions.  If
 * built for CPU, tiling will be enabled.  However, one could specify a huge
 * tile size to effectively disable tiling.  For GPU builds, this function is
 * NON-BLOCKING on the host. Conceptually, this is a 4D loop.
 *
 * \tparam MT max threads in GPU blocks (Only relevant for GPU builds)
 * \tparam MF the MultiFab/FabArray type
 * \tparam F a callable type like lambda
 *
 * \param mf the MultiFab/FabArray object used to specify the iteration space
 * \param ng the number of ghost cells around the valid region
 * \param ts tile size, ignored by GPU build
 * \param dt controls dynamic tiling for the cpu build
 * \param f a callable object void(int,int,int,int), where the first argument
 *           is the local box index, and the following three are spatial indices
 *           for x, y, and z-directions.
 */
template <int MT, typename MF, typename F>
std::enable_if_t<IsFabArray<MF>::value>
ParallelFor (MF const& mf, IntVect const& ng, TileSize const& ts,
             DynamicTiling dt, F&& f)
{
#ifdef AMREX_USE_GPU
    detail::ParallelFor<MT>(mf, ng, ts.tile_size, dt.dynamic, std::forward<F>(f));
#else
    detail::ParallelFor(mf, ng, ts.tile_size, dt.dynamic, std::forward<F>(f));
#endif
}

/**
 * \brief ParallelFor for MultiFab/FabArray.
 *
 * This version launches a kernel to work on the valid and ghost regions.  If
 * built for CPU, tiling will be enabled.  However, one could specify a huge
 * tile size to effectively disable tiling.  For gpu builds, this function is
 * NON-BLOCKING on the host. Conceptually, this is a 5D loop.
 *
 * \tparam MT max threads in GPU blocks (Only relevant for GPU builds)
 * \tparam MF the MultiFab/FabArray type
 * \tparam F a callable type like lambda
 *
 * \param mf the MultiFab/FabArray object used to specify the iteration space
 * \param ng the number of ghost cells around the valid region
 * \param ncomp the number of component
 * \param ts tile size, ignored by GPU build
 * \param dt controls dynamic tiling for the cpu build
 * \param f a callable object void(int,int,int,int,int), where the first argument
 *           is the local box index, the following three are spatial indices
 *           for x, y, and z-directions, and the last is for component.
 */
template <int MT, typename MF, typename F>
std::enable_if_t<IsFabArray<MF>::value>
ParallelFor (MF const& mf, IntVect const& ng, int ncomp, TileSize const& ts,
             DynamicTiling dt, F&& f)
{
#ifdef AMREX_USE_GPU
    detail::ParallelFor<MT>(mf, ng, ncomp, ts.tile_size, dt.dynamic, std::forward<F>(f));
#else
    detail::ParallelFor(mf, ng, ncomp, ts.tile_size, dt.dynamic, std::forward<F>(f));
#endif
}

}

using experimental::ParallelFor;

}

#endif
