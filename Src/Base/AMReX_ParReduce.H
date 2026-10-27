#ifndef AMREX_PAR_REDUCE_H_
#define AMREX_PAR_REDUCE_H_
#include <AMReX_Config.H>

#include <AMReX_Concepts.H>
#include <AMReX_Reduce.H>

namespace amrex {

/**
 * \brief Parallel reduce for MultiFab/FabArray. The reduce result is local
 * and it's the user's responsibility if MPI communication is needed.
 *
 * This performs reduction over a MultiFab's valid and specified ghost regions. For
 * example, the code below computes the minimum of the first MultiFab and the maximum of
 * the second MultiFab.
 * \code{.cpp}
 *     auto const& ma1 = mf1.const_arrays();
 *     auto const& ma2 = mf2.const_arrays();
 *     GpuTuple<Real,Real> mm = ParReduce(TypeList<ReduceOpMin,ReduceOpMax>{},
 *                                        TypeList<Real,Real>{},
 *                                        mf1, mf1.nGrowVect(),
 *     [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
 *         -> GpuTuple<Real,Real>
 *     {
 *         return { ma1[box_no](i,j,k), ma2[box_no](i,j,k) };
 *     });
 * \endcode
 *
 * \tparam Ops... reduce operators (e.g., ReduceOpSum, ReduceOpMin, ReduceOpMax,
 *                ReduceOpLogicalAnd, and ReduceOpLogicalOr)
 * \tparam Ts...  data types (e.g., Real, int, etc.)
 * \tparam FAB    MultiFab/FabArray type
 * \tparam F      callable type like a lambda function
 *
 * \param operation_list list of reduce operators
 * \param type_list      list of data types
 * \param fa     a MultiFab/FabArray object used to specify the iteration space
 * \param nghost the number of ghost cells included in the iteration space
 * \param f      a callable object returning GpuTuple<Ts...>.  It takes four ints,
 *               where the first int is the local box index and the others are
 *               spatial indices for x, y, and z-directions.
 *
 * \return reduction result (GpuTuple<Ts...>)
 */
template <typename... Ops, typename... Ts, BaseFabType FAB, typename F>
typename ReduceData<Ts...>::Type
ParReduce (TypeList<Ops...> operation_list, TypeList<Ts...> type_list,
           FabArray<FAB> const& fa, IntVect const& nghost, F&& f)
{
    amrex::ignore_unused(operation_list,type_list);
    ReduceOps<Ops...> reduce_op;
    ReduceData<Ts...> reduce_data(reduce_op);
    reduce_op.eval(fa, nghost, reduce_data, std::forward<F>(f));
    return reduce_data.value(reduce_op);
}

/**
 * \brief Parallel reduce for MultiFab/FabArray. The reduce result is local
 * and it's the user's responsibility if MPI communication is needed.
 *
 * This performs reduction over a MultiFab's valid and specified ghost regions. For
 * example, the code below computes the sum of the processed data in a MultiFab.
 * \code{.cpp}
 *     auto const& ma = mf.const_arrays();
 *     Real ektot = ParReduce(TypeList<ReduceOpSum>{}, TypeList<Real>{},
 *                            mf, IntVect(0),
 *     [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
 *         -> GpuTuple<Real>
 *     {
 *         auto rho = ma[box_no](i,j,k,0);
 *         auto mx = ma[box_no](i,j,k,1);
 *         auto my = ma[box_no](i,j,k,2);
 *         auto mz = ma[box_no](i,j,k,3);
 *         auto ek = (mx*mx+my*my+mz*mz)/(2.*rho);
 *         return { ek };
 *     });
 * \endcode
 *
 * \tparam Op  Reduce operator (e.g., ReduceOpSum, ReduceOpMin, ReduceOpMax,
 *             ReduceOpLogicalAnd, and ReduceOpLogicalOr)
 * \tparam T   data type (e.g., Real, int, etc.)
 * \tparam FAB MultiFab/FabArray type
 * \tparam F   callable type like a lambda function
 *
 * \param operation_list a reduce operator stored in TypeList
 * \param type_list      a data type stored in TypeList
 * \param fa     a MultiFab/FabArray object used to specify the iteration space
 * \param nghost the number of ghost cells included in the iteration space
 * \param f      a callable object returning GpuTuple<T>.  It takes four ints,
 *               where the first int is the local box index and the others are
 *               spatial indices for x, y, and z-directions.
 *
 * \return reduction result (T)
 */
template <typename Op, typename T, BaseFabType FAB, typename F>
T
ParReduce (TypeList<Op> operation_list, TypeList<T> type_list,
           FabArray<FAB> const& fa, IntVect const& nghost, F&& f)
{
    amrex::ignore_unused(operation_list,type_list);
    ReduceOps<Op> reduce_op;
    ReduceData<T> reduce_data(reduce_op);
    reduce_op.eval(fa, nghost, reduce_data, std::forward<F>(f));
    auto const& hv = reduce_data.value(reduce_op);
    return amrex::get<0>(hv);
}

/**
 * \brief Parallel reduce for MultiFab/FabArray. The reduce result is local
 * and it's the user's responsibility if MPI communication is needed.
 *
 * This performs reduction over a MultiFab's valid and specified ghost regions and
 * components.  For example, the code below computes the minimum of the first MultiFab
 * and the maximum of the second MultiFab.
 * \code{.cpp}
 *     auto const& ma1 = mf1.const_arrays();
 *     auto const& ma2 = mf2.const_arrays();
 *     GpuTuple<Real,Real> mm = ParReduce(TypeList<ReduceOpMin,ReduceOpMax>{},
 *                                        TypeList<Real,Real>{},
 *                                        mf1, mf1.nGrowVect(), mf1.nComp(),
 *     [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
 *         -> GpuTuple<Real,Real>
 *     {
 *         return { ma1[box_no](i,j,k,n), ma2[box_no](i,j,k,n) };
 *     });
 * \endcode
 *
 * \tparam Ops... reduce operators (e.g., ReduceOpSum, ReduceOpMin, ReduceOpMax,
 *                ReduceOpLogicalAnd, and ReduceOpLogicalOr)
 * \tparam Ts...  data types (e.g., Real, int, etc.)
 * \tparam FAB    MultiFab/FabArray type
 * \tparam F      callable type like a lambda function
 *
 * \param operation_list list of reduce operators
 * \param type_list      list of data types
 * \param fa     a MultiFab/FabArray object used to specify the iteration space
 * \param nghost the number of ghost cells included in the iteration space
 * \param ncomp  the number of components in the iteration space
 * \param f      a callable object returning GpuTuple<Ts...>.  It takes five ints,
 *               where the first int is the local box index, the next three are
 *               spatial indices for x, y, and z-directions, and the last is for
 *               component.
 *
 * \return reduction result (GpuTuple<Ts...>)
 */
template <typename... Ops, typename... Ts, BaseFabType FAB, typename F>
typename ReduceData<Ts...>::Type
ParReduce (TypeList<Ops...> operation_list, TypeList<Ts...> type_list,
           FabArray<FAB> const& fa, IntVect const& nghost, int ncomp, F&& f)
{
    amrex::ignore_unused(operation_list,type_list);
    ReduceOps<Ops...> reduce_op;
    ReduceData<Ts...> reduce_data(reduce_op);
    reduce_op.eval(fa, nghost, ncomp, reduce_data, std::forward<F>(f));
    return reduce_data.value(reduce_op);
}

/**
 * \brief Parallel reduce for MultiFab/FabArray. The reduce result is local
 * and it's the user's responsibility if MPI communication is needed.
 *
 * This performs reduction over a MultiFab's valid and specified ghost regions. For
 * example, the code below computes the sum of the data in a MultiFab.
 * \code{.cpp}
 *     auto const& ma = mf.const_arrays();
 *     Real ektot = ParReduce(TypeList<ReduceOpSum>{}, TypeList<Real>{},
 *                            mf, mf.nGrowVect(), mf.nComp(),
 *     [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
 *         -> GpuTuple<Real>
 *     {
 *         return { ma[box_no](i,j,k,n) };
 *     });
 * \endcode
 *
 * \tparam Op  Reduce operator (e.g., ReduceOpSum, ReduceOpMin, ReduceOpMax,
 *             ReduceOpLogicalAnd, and ReduceOpLogicalOr)
 * \tparam T   data type (e.g., Real, int, etc.)
 * \tparam FAB MultiFab/FabArray type
 * \tparam F   callable type like a lambda function
 *
 * \param operation_list a reduce operator stored in TypeList
 * \param type_list      a data type stored in TypeList
 * \param fa     a MultiFab/FabArray object used to specify the iteration space
 * \param nghost the number of ghost cells included in the iteration space
 * \param ncomp  the number of components in the iteration space
 * \param f      a callable object returning GpuTuple<T>.  It takes five ints,
 *               where the first int is the local box index, the next three are
 *               spatial indices for x, y, and z-directions, and the last is for
 *               component.
 *
 * \return reduction result (T)
 */
template <typename Op, typename T, BaseFabType FAB, typename F>
T
ParReduce (TypeList<Op> operation_list, TypeList<T> type_list,
           FabArray<FAB> const& fa, IntVect const& nghost, int ncomp, F&& f)
{
    amrex::ignore_unused(operation_list,type_list);
    ReduceOps<Op> reduce_op;
    ReduceData<T> reduce_data(reduce_op);
    reduce_op.eval(fa, nghost, ncomp, reduce_data, std::forward<F>(f));
    auto const& hv = reduce_data.value(reduce_op);
    return amrex::get<0>(hv);
}

/**
 * \brief Parallel reduce for MultiFab/FabArray. The reduce result is local
 * and it's the user's responsibility if MPI communication is needed.
 *
 * This performs reduction over a MultiFab's valid region. For
 * example, the code below computes the minimum of the first MultiFab and the maximum of
 * the second MultiFab.
 * \code{.cpp}
 *     auto const& ma1 = mf1.const_arrays();
 *     auto const& ma2 = mf2.const_arrays();
 *     GpuTuple<Real,Real> mm = ParReduce(TypeList<ReduceOpMin,ReduceOpMax>{},
 *                                        TypeList<Real,Real>{}, mf1,
 *     [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
 *         -> GpuTuple<Real,Real>
 *     {
 *         return { ma1[box_no](i,j,k), ma2[box_no](i,j,k) };
 *     });
 * \endcode
 *
 * \tparam Ops... reduce operators (e.g., ReduceOpSum, ReduceOpMin, ReduceOpMax,
 *                ReduceOpLogicalAnd, and ReduceOpLogicalOr)
 * \tparam Ts...  data types (e.g., Real, int, etc.)
 * \tparam FAB    MultiFab/FabArray type
 * \tparam F      callable type like a lambda function
 *
 * \param operation_list list of reduce operators
 * \param type_list      list of data types
 * \param fa     a MultiFab/FabArray object used to specify the iteration space
 * \param f      a callable object returning GpuTuple<Ts...>.  It takes four ints,
 *               where the first int is the local box index and the others are
 *               spatial indices for x, y, and z-directions.
 *
 * \return reduction result (GpuTuple<Ts...>)
 */
template <typename... Ops, typename... Ts, BaseFabType FAB, typename F>
typename ReduceData<Ts...>::Type
ParReduce (TypeList<Ops...> operation_list, TypeList<Ts...> type_list,
           FabArray<FAB> const& fa, F&& f)
{
    return ParReduce(operation_list, type_list, fa, IntVect(0), std::forward<F>(f));
}

/**
 * \brief Parallel reduce for MultiFab/FabArray. The reduce result is local
 * and it's the user's responsibility if MPI communication is needed.
 *
 * This performs reduction over a MultiFab's valid region. For
 * example, the code below computes the sum of the processed data in a MultiFab.
 * \code{.cpp}
 *     auto const& ma = mf.const_arrays();
 *     Real ektot = ParReduce(TypeList<ReduceOpSum>{}, TypeList<Real>{}, mf,
 *     [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
 *         -> GpuTuple<Real>
 *     {
 *         auto rho = ma[box_no](i,j,k,0);
 *         auto mx = ma[box_no](i,j,k,1);
 *         auto my = ma[box_no](i,j,k,2);
 *         auto mz = ma[box_no](i,j,k,3);
 *         auto ek = (mx*mx+my*my+mz*mz)/(2.*rho);
 *         return { ek };
 *     });
 * \endcode
 *
 * \tparam Op  Reduce operator (e.g., ReduceOpSum, ReduceOpMin, ReduceOpMax,
 *             ReduceOpLogicalAnd, and ReduceOpLogicalOr)
 * \tparam T   data type (e.g., Real, int, etc.)
 * \tparam FAB MultiFab/FabArray type
 * \tparam F   callable type like a lambda function
 *
 * \param operation_list a reduce operator stored in TypeList
 * \param type_list      a data type stored in TypeList
 * \param fa     a MultiFab/FabArray object used to specify the iteration space
 * \param f      a callable object returning GpuTuple<T>.  It takes four ints,
 *               where the first int is the local box index and the others are
 *               spatial indices for x, y, and z-directions.
 *
 * \return reduction result (T)
 */
template <typename Op, typename T, BaseFabType FAB, typename F>
T
ParReduce (TypeList<Op> operation_list, TypeList<T> type_list,
           FabArray<FAB> const& fa, F&& f)
{
    return ParReduce(operation_list, type_list, fa, IntVect(0), std::forward<F>(f));
}

/// \cond DOXYGEN_IGNORE
namespace Reduce::detail {

// Reduce f(box_no,i,j,k,n) separately for each component n (PerBox false) or
// for each local box and component (PerBox true, result index box_no*ncomp+n).
// A class template, so that nvcc can name the function of the device lambdas.
template <bool PerBox, typename... Ops>
struct ParReduceComp
{
    template <typename RT, BaseFabType FAB, typename F>
    static Vector<RT>
    eval (FabArray<FAB> const& fa, IntVect const& nghost, int ncomp, F const& f)
    {
        const int nlocal = fa.local_size();
        const int nkeys = PerBox ? nlocal*ncomp : ncomp;
        Vector<RT> r(std::max(nkeys,0));
        for (auto& x : r) { for_each_init<0, RT, Ops...>(x); }
        if (nlocal == 0 || ncomp <= 0) { return r; }

#ifdef AMREX_USE_GPU
        // Block b works on component n = b/nb1 of local box L and stores its
        // partial result in slot b, so each key's slots are contiguous.  Box L
        // owns blocks [hoff[L],hoff[L+1]) of every component.
        constexpr int MT = AMREX_GPU_MAX_THREADS;
        const int maxblocks = int(Gpu::Device::maxBlocksPerLaunch());
        Long ntot = 0;
        for (int L = 0; L < nlocal; ++L) {
            Box b = fa.box(fa.IndexArray()[L]);
            ntot += b.ok() ? amrex::grow(b,nghost).numPts() : Long(0);
        }
        // About one launch worth of blocks, with at least 4 cells per thread.
        const Long cells_per_block = std::max(Long(MT)*4, (ntot*ncomp + maxblocks - 1) / maxblocks);
        Gpu::PinnedVector<int> hoff(nlocal+1, 0);
        for (int L = 0; L < nlocal; ++L) {
            Box b = fa.box(fa.IndexArray()[L]);
            Long const npts = b.ok() ? amrex::grow(b,nghost).numPts() : Long(0);
            Long const nb = (npts + cells_per_block - 1) / cells_per_block;
            AMREX_ALWAYS_ASSERT(Long(hoff[L]) + nb < Long(std::numeric_limits<int>::max()));
            hoff[L+1] = hoff[L] + int(nb);
        }
        const int nb1 = hoff[nlocal];
        AMREX_ALWAYS_ASSERT(Long(nb1)*ncomp < Long(std::numeric_limits<int>::max()));
        const int nblocks = nb1*ncomp;
        if (nblocks == 0) { return r; }

        auto const& stream = Gpu::gpuStream();
        Gpu::DeviceVector<int> doff(nlocal+1);
        Gpu::copyAsync(Gpu::hostToDevice, hoff.begin(), hoff.end(), doff.begin());
        Gpu::DeviceVector<RT> dpart(nblocks);
        Gpu::DeviceVector<RT> dr(nkeys);
        int const* poff = doff.data();
        RT* ppart = dpart.data();
        RT* pr = dr.data();
        BoxIndexer const* pbx = fa.getParForInfo(nghost).getBoxes();

#ifdef AMREX_USE_SYCL
        constexpr std::size_t shared_mem_bytes = sizeof(RT)
            * std::max(Gpu::Device::warp_size, MT/Gpu::Device::warp_size);
        amrex::launch<MT>(std::min(nblocks,maxblocks), shared_mem_bytes, stream,
        [=] AMREX_GPU_DEVICE (Gpu::Handler const& gh) noexcept
        {
            Dim1 blockIdx {gh.blockIdx()};
            Dim1 threadIdx{gh.threadIdx()};
            Dim1 gridDim  {gh.gridDim()};
#else
        amrex::launch<MT>(std::min(nblocks,maxblocks), 0, stream,
        [=] AMREX_GPU_DEVICE () noexcept
        {
#endif
            for (int b = int(blockIdx.x); b < nblocks; b += int(gridDim.x)) {
                const int n = b / nb1;
                const int bb = b - n*nb1;
                const int L = amrex::bisect(poff, 0, nlocal, bb);
                const int nb = poff[L+1] - poff[L];
                BoxIndexer const& indexer = pbx[L];
                RT t;
                for_each_init<0, RT, Ops...>(t);
                for (auto icell = std::uint64_t(bb-poff[L])*MT + threadIdx.x;
                     icell < indexer.numPts(); icell += std::uint64_t(nb)*MT)
                {
                    auto [i, j, k] = indexer(icell);
                    for_each_local<0, RT, Ops...>(t, f(L,i,j,k,n));
                }
                RT& dst = ppart[b];
                if (threadIdx.x == 0) { for_each_init<0, RT, Ops...>(dst); }
#ifdef AMREX_USE_SYCL
                for_each_parallel<0, RT, Ops...>(dst, t, gh);
#else
                for_each_parallel<0, RT, Ops...>(dst, t);
#endif
            }
        });

        // Combine each key's partial results in a fixed order.
#ifdef AMREX_USE_SYCL
        amrex::launch<MT>(std::min(nkeys,maxblocks), shared_mem_bytes, stream,
        [=] AMREX_GPU_DEVICE (Gpu::Handler const& gh) noexcept
        {
            Dim1 blockIdx {gh.blockIdx()};
            Dim1 threadIdx{gh.threadIdx()};
            Dim1 gridDim  {gh.gridDim()};
#else
        amrex::launch<MT>(std::min(nkeys,maxblocks), 0, stream,
        [=] AMREX_GPU_DEVICE () noexcept
        {
#endif
            for (int key = int(blockIdx.x); key < nkeys; key += int(gridDim.x)) {
                int begin = 0;
                int end = 0;
                if (PerBox) { // not constexpr: nvcc cannot first-capture ncomp there
                    const int L = key / ncomp;
                    const int n = key - L*ncomp;
                    begin = n*nb1 + poff[L];
                    end   = n*nb1 + poff[L+1];
                } else {
                    begin = key*nb1;
                    end   = begin + nb1;
                }
                RT t;
                for_each_init<0, RT, Ops...>(t);
                for (int ib = begin + int(threadIdx.x); ib < end; ib += MT) {
                    for_each_local<0, RT, Ops...>(t, ppart[ib]);
                }
                RT& dst = pr[key];
                if (threadIdx.x == 0) { for_each_init<0, RT, Ops...>(dst); }
#ifdef AMREX_USE_SYCL
                for_each_parallel<0, RT, Ops...>(dst, t, gh);
#else
                for_each_parallel<0, RT, Ops...>(dst, t);
#endif
            }
        });

        Gpu::dtoh_memcpy_async(r.data(), pr, sizeof(RT)*nkeys);
        Gpu::streamSynchronize();
#else
        if constexpr (PerBox) {
            // With too few (box, component) pairs for the threads, split each
            // box along its last dimension into contiguous chunks of about 16K
            // cells. Items are ordered by box, component and chunk, which is
            // the memory order, and their partial results are combined in order.
            constexpr int d = AMREX_SPACEDIM-1;
            const int nthreads = OpenMP::in_parallel() ? 1 : OpenMP::get_max_threads();
            const bool split = nthreads > 1 && nkeys < 4*nthreads;
            Vector<Box> bx(nlocal);
            Vector<int> thick(nlocal), nchunk(nlocal), ioff(nlocal+1, 0);
            for (int L = 0; L < nlocal; ++L) {
                bx[L] = fa.box(fa.IndexArray()[L]);
                if (bx[L].ok()) { bx[L].grow(nghost); }
                // An empty box is one empty chunk.
                const int len = bx[L].ok() ? bx[L].length(d) : 1;
                thick[L] = (split && bx[L].ok()) ? int(std::max(Long(1), Long(16384)/(bx[L].numPts()/len))) : len;
                nchunk[L] = (len + thick[L] - 1) / thick[L];
                ioff[L+1] = ioff[L] + nchunk[L]*ncomp;
            }
            // Without splitting, item is the result index.
            const int nitems = ioff[nlocal];
            Vector<RT> part(split ? nitems : 0);
            RT* pr = split ? part.data() : r.data();
#ifdef AMREX_USE_OMP
#pragma omp parallel for schedule(dynamic) if (!system::regtest_reduction)
#endif
            for (int item = 0; item < nitems; ++item) {
                const int L = split ? amrex::bisect(ioff.data(), 0, nlocal, item) : item/ncomp;
                const int m = item - ioff[L];
                const int n = m / nchunk[L];
                const int c = m - n*nchunk[L];
                Box b = bx[L];
                b.setSmall(d, bx[L].smallEnd(d) + c*thick[L]);
                b.setBig(d, std::min(b.smallEnd(d)+thick[L]-1, bx[L].bigEnd(d)));
                const auto lo = amrex::lbound(b);
                const auto hi = amrex::ubound(b);
                RT rem;
                for_each_init<0, RT, Ops...>(rem);
                Lanes<RT> acc;
                lanes_set<0>(acc, rem);
                for (int k = lo.z; k <= hi.z; ++k) {
                for (int j = lo.y; j <= hi.y; ++j) {
                    reduce_row<Ops...>(lo.x, hi.x, acc, rem,
                        [&] (int i) { return f(L,i,j,k,n); });
                }}
                RT& x = pr[item];
                for_each_init<0, RT, Ops...>(x);
                lanes_fold<0, RT, Ops...>(acc, x);
                for_each_local<0, RT, Ops...>(x, rem);
            }
            if (split) {
                for (int L = 0; L < nlocal; ++L) {
                    for (int m = 0; m < nchunk[L]*ncomp; ++m) {
                        for_each_local<0, RT, Ops...>(r[L*ncomp+m/nchunk[L]], part[ioff[L]+m]);
                    }
                }
            }
        } else {
            const int nthreads = OpenMP::in_parallel() ? 1 : OpenMP::get_max_threads();
            Vector<RT> rt(std::size_t(nthreads)*ncomp);
            for (auto& x : rt) { for_each_init<0, RT, Ops...>(x); }
#ifdef AMREX_USE_OMP
#pragma omp parallel if (nthreads > 1 && !system::regtest_reduction)
#endif
            {
                RT* rr = rt.data() + std::size_t(nthreads > 1 ? OpenMP::get_thread_num() : 0)*ncomp;
                // Component loop outermost, so that each pass streams one component.
                for (int n = 0; n < ncomp; ++n) {
                    RT rem;
                    for_each_init<0, RT, Ops...>(rem);
                    Lanes<RT> acc;
                    lanes_set<0>(acc, rem);
                    for (MFIter mfi(fa,true); mfi.isValid(); ++mfi) {
                        Box const& b = mfi.growntilebox(nghost);
                        const int L = mfi.LocalIndex();
                        const auto lo = amrex::lbound(b);
                        const auto hi = amrex::ubound(b);
                        for (int k = lo.z; k <= hi.z; ++k) {
                        for (int j = lo.y; j <= hi.y; ++j) {
                            reduce_row<Ops...>(lo.x, hi.x, acc, rem,
                                [&] (int i) { return f(L,i,j,k,n); });
                        }}
                    }
                    lanes_fold<0, RT, Ops...>(acc, rr[n]);
                    for_each_local<0, RT, Ops...>(rr[n], rem);
                }
            }
            for (int it = 0; it < nthreads; ++it) {
                for (int n = 0; n < ncomp; ++n) {
                    for_each_local<0, RT, Ops...>(r[n], rt[std::size_t(it)*ncomp+n]);
                }
            }
        }
#endif
        return r;
    }
};

}
/// \endcond

/**
 * \brief Parallel reduce for MultiFab/FabArray with a separate result for
 * each component. The reduce result is local and it's the user's
 * responsibility if MPI communication is needed.
 *
 * Unlike ParReduce with ncomp, which also reduces over the components, this
 * returns ncomp results. For example, the code below computes the minimum
 * and maximum of every component of a MultiFab.
 * \code{.cpp}
 *     auto const& ma = mf.const_arrays();
 *     Vector<GpuTuple<Real,Real>> mm =
 *         ParReducePerComp(TypeList<ReduceOpMin,ReduceOpMax>{},
 *                          TypeList<Real,Real>{},
 *                          mf, IntVect(0), mf.nComp(),
 *         [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
 *             -> GpuTuple<Real,Real>
 *         {
 *             return { ma[box_no](i,j,k,n), ma[box_no](i,j,k,n) };
 *         });
 * \endcode
 *
 * \param operation_list list of reduce operators
 * \param type_list      list of data types
 * \param fa     a MultiFab/FabArray object used to specify the iteration space
 * \param nghost the number of ghost cells included in the iteration space
 * \param ncomp  the number of components
 * \param f      a callable object returning GpuTuple<Ts...>.  It takes five ints:
 *               the local box index, the spatial indices, and the component.
 *
 * \return reduction results, one per component (Vector<GpuTuple<Ts...>>)
 */
template <typename... Ops, typename... Ts, BaseFabType FAB, typename F>
Vector<typename ReduceData<Ts...>::Type>
ParReducePerComp (TypeList<Ops...> operation_list, TypeList<Ts...> type_list,
                  FabArray<FAB> const& fa, IntVect const& nghost, int ncomp, F&& f)
{
    amrex::ignore_unused(operation_list,type_list);
    return Reduce::detail::ParReduceComp<false,Ops...>::template
        eval<typename ReduceData<Ts...>::Type>(fa, nghost, ncomp, std::forward<F>(f));
}

/**
 * \brief ParReducePerComp with a single reduce operator.
 *
 * \return reduction results, one per component (Vector<T>)
 */
template <typename Op, typename T, BaseFabType FAB, typename F>
Vector<T>
ParReducePerComp (TypeList<Op> operation_list, TypeList<T> type_list,
                  FabArray<FAB> const& fa, IntVect const& nghost, int ncomp, F&& f)
{
    amrex::ignore_unused(operation_list,type_list);
    auto const& rt = Reduce::detail::ParReduceComp<false,Op>::template
        eval<typename ReduceData<T>::Type>(fa, nghost, ncomp, std::forward<F>(f));
    Vector<T> r(rt.size());
    for (int i = 0, N = int(rt.size()); i < N; ++i) { r[i] = amrex::get<0>(rt[i]); }
    return r;
}

/**
 * \brief Parallel reduce for MultiFab/FabArray with a separate result for
 * each local box and component. The reduce result is local.
 *
 * This takes the same arguments as ParReducePerComp.  The result for
 * component n of local box box_no is at index box_no*ncomp+n.
 *
 * \return reduction results (Vector<GpuTuple<Ts...>> of size local_size()*ncomp)
 */
template <typename... Ops, typename... Ts, BaseFabType FAB, typename F>
Vector<typename ReduceData<Ts...>::Type>
ParReducePerBoxComp (TypeList<Ops...> operation_list, TypeList<Ts...> type_list,
                     FabArray<FAB> const& fa, IntVect const& nghost, int ncomp, F&& f)
{
    amrex::ignore_unused(operation_list,type_list);
    return Reduce::detail::ParReduceComp<true,Ops...>::template
        eval<typename ReduceData<Ts...>::Type>(fa, nghost, ncomp, std::forward<F>(f));
}

/**
 * \brief ParReducePerBoxComp with a single reduce operator.
 *
 * \return reduction results (Vector<T> of size local_size()*ncomp)
 */
template <typename Op, typename T, BaseFabType FAB, typename F>
Vector<T>
ParReducePerBoxComp (TypeList<Op> operation_list, TypeList<T> type_list,
                     FabArray<FAB> const& fa, IntVect const& nghost, int ncomp, F&& f)
{
    amrex::ignore_unused(operation_list,type_list);
    auto const& rt = Reduce::detail::ParReduceComp<true,Op>::template
        eval<typename ReduceData<T>::Type>(fa, nghost, ncomp, std::forward<F>(f));
    Vector<T> r(rt.size());
    for (int i = 0, N = int(rt.size()); i < N; ++i) { r[i] = amrex::get<0>(rt[i]); }
    return r;
}

}
#endif
