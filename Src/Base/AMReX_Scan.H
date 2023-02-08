#ifndef AMREX_SCAN_H_
#define AMREX_SCAN_H_
#include <AMReX_Config.H>

#include <AMReX_Extension.H>
#include <AMReX_Gpu.H>
#include <AMReX_Arena.H>

#if defined(AMREX_USE_CUDA) && defined(__CUDACC__) && (__CUDACC_VER_MAJOR__ >= 11)
#  include <cub/cub.cuh>
#elif defined(AMREX_USE_HIP)
#  include <rocprim/rocprim.hpp>
#elif defined(AMREX_USE_SYCL) && defined(AMREX_USE_ONEDPL)
#  include <oneapi/dpl/execution>
#  include <oneapi/dpl/numeric>
#endif

#include <cstdint>
#include <numeric>
#include <type_traits>

namespace amrex {
namespace Scan {

struct RetSum {
    bool flag = true;
    explicit operator bool() const noexcept { return flag; }
};
static constexpr RetSum   retSum{true};
static constexpr RetSum noRetSum{false};

namespace Type {
    static constexpr struct Inclusive {} inclusive{};
    static constexpr struct Exclusive {} exclusive{};
}

#if defined(AMREX_USE_GPU)

namespace detail {

template <typename T>
struct STVA
{
    char status;
    T value;
};

template <typename T, bool SINGLE_WORD> struct BlockStatus {};

template <typename T>
struct BlockStatus<T, true>
{
    template<typename U>
    union Data {
        STVA<U> s;
        uint64_t i;
        void operator=(Data<U> const&) = delete;
        void operator=(Data<U> &&) = delete;
    };
    Data<T> d;

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void write (char a_status, T a_value) {
#if defined(AMREX_USE_CUDA)
        volatile uint64_t tmp;
        reinterpret_cast<STVA<T> volatile&>(tmp).status = a_status;
        reinterpret_cast<STVA<T> volatile&>(tmp).value  = a_value;
        reinterpret_cast<uint64_t&>(d.s) = tmp;
#else
        Data<T> tmp;
        tmp.s = {a_status, a_value};
        static_assert(sizeof(unsigned long long) == sizeof(uint64_t),
                      "HIP/SYCL: unsigned long long must be 64 bits");
        Gpu::Atomic::Exch(reinterpret_cast<unsigned long long*>(&d),
                          reinterpret_cast<unsigned long long&>(tmp));
#endif
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    T get_aggregate() const { return d.s.value; }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    STVA<T> read () volatile {
#if defined(AMREX_USE_CUDA)
        volatile uint64_t tmp = reinterpret_cast<uint64_t volatile&>(d);
        return {reinterpret_cast<STVA<T> volatile&>(tmp).status,
                reinterpret_cast<STVA<T> volatile&>(tmp).value };
#else
        static_assert(sizeof(unsigned long long) == sizeof(uint64_t),
                      "HIP/SYCL: unsigned long long must be 64 bits");
        unsigned long long tmp = Gpu::Atomic::Add
            (reinterpret_cast<unsigned long long*>(const_cast<Data<T>*>(&d)), 0ull);
        return (*reinterpret_cast<Data<T>*>(&tmp)).s;
#endif
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void set_status (char a_status) { d.s.status = a_status; }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    STVA<T> wait () volatile {
        STVA<T> r;
        do {
#if defined(AMREX_USE_SYCL)
            sycl::atomic_fence(sycl::memory_order::acq_rel, sycl::memory_scope::work_group);
#else
            __threadfence_block();
#endif
            r = read();
        } while (r.status == 'x');
        return r;
    }
};

template <typename T>
struct BlockStatus<T, false>
{
    T aggregate;
    T inclusive;
    char status;

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void write (char a_status, T a_value) {
        if (a_status == 'a') {
            aggregate = a_value;
        } else {
            inclusive = a_value;
        }
#if defined(AMREX_USE_SYCL)
        sycl::atomic_fence(sycl::memory_order::acq_rel, sycl::memory_scope::device);
#else
        __threadfence();
#endif
        status = a_status;
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    T get_aggregate() const { return aggregate; }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    STVA<T> read () volatile {
        if (status == 'x') {
            return {'x', 0};
        } else if (status == 'a') {
            return {'a', aggregate};
        } else {
            return {'p', inclusive};
        }
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void set_status (char a_status) { status = a_status; }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    STVA<T> wait () volatile {
        STVA<T> r;
        do {
            r = read();
#if defined(AMREX_USE_SYCL)
            sycl::atomic_fence(sycl::memory_order::acq_rel, sycl::memory_scope::device);
#else
            __threadfence();
#endif
        } while (r.status == 'x');
        return r;
    }
};

}

#if defined(AMREX_USE_SYCL)

#ifndef AMREX_SYCL_NO_MULTIPASS_SCAN
template <typename T, typename N, typename FIN, typename FOUT, typename TYPE>
T PrefixSum_mp (N n, FIN && fin, FOUT && fout, TYPE, RetSum a_ret_sum)
{
    if (n <= 0) return 0;
    constexpr int nwarps_per_block = 8;
    constexpr int nthreads = nwarps_per_block*Gpu::Device::warp_size;
    constexpr int nchunks = 12;
    constexpr int nelms_per_block = nthreads * nchunks;
    AMREX_ALWAYS_ASSERT(static_cast<Long>(n) < static_cast<Long>(std::numeric_limits<int>::max())*nelms_per_block);
    int nblocks = (static_cast<Long>(n) + nelms_per_block - 1) / nelms_per_block;
    std::size_t sm = sizeof(T) * (Gpu::Device::warp_size + nwarps_per_block);
    auto stream = Gpu::gpuStream();

    std::size_t nbytes_blockresult = Arena::align(sizeof(T)*n);
    std::size_t nbytes_blocksum = Arena::align(sizeof(T)*nblocks);
    std::size_t nbytes_totalsum = Arena::align(sizeof(T));
    auto dp = (char*)(The_Arena()->alloc(nbytes_blockresult
                                         + nbytes_blocksum
                                         + nbytes_totalsum));
    T* blockresult_p = (T*)dp;
    T* blocksum_p = (T*)(dp + nbytes_blockresult);
    T* totalsum_p = (T*)(dp + nbytes_blockresult + nbytes_blocksum);

    amrex::launch(nblocks, nthreads, sm, stream,
    [=] AMREX_GPU_DEVICE (Gpu::Handler const& gh) noexcept
    {
        sycl::sub_group const& sg = gh.item->get_sub_group();
        int lane = sg.get_local_id()[0];
        int warp = sg.get_group_id()[0];
        int nwarps = sg.get_group_range()[0];

        int threadIdxx = gh.item->get_local_id(0);
        int blockIdxx = gh.item->get_group_linear_id();
        int blockDimx = gh.item->get_local_range(0);

        T* shared = (T*)(gh.local);
        T* shared2 = shared + Gpu::Device::warp_size;

        // Each block processes [ibegin,iend).
        N ibegin = nelms_per_block * blockIdxx;
        N iend = amrex::min(static_cast<N>(ibegin+nelms_per_block), n);

        // Each block is responsible for nchunks chunks of data,
        // where each chunk has blockDim.x elements, one for each
        // thread in the block.
        T sum_prev_chunk = 0; // inclusive sum from previous chunks.
        for (int ichunk = 0; ichunk < nchunks; ++ichunk) {
            N offset = ibegin + ichunk*blockDimx;
            if (offset >= iend) break;

            offset += threadIdxx;
            T x0 = (offset < iend) ? fin(offset) : 0;
            T x = x0;
            // Scan within a warp
            for (int i = 1; i <= Gpu::Device::warp_size; i *= 2) {
                T s = sycl::shift_group_right(sg, x, i);
                if (lane >= i) x += s;
            }

            // x now holds the inclusive sum within the warp.  The
            // last thread in each warp holds the inclusive sum of
            // this warp.  We will store it in shared memory.
            if (lane == Gpu::Device::warp_size - 1) {
                shared[warp] = x;
            }

            gh.item->barrier(sycl::access::fence_space::local_space);

            // The first warp will do scan on the warp sums for the
            // whole block.
            if (warp == 0) {
                T y = (lane < nwarps) ? shared[lane] : 0;
                for (int i = 1; i <= Gpu::Device::warp_size; i *= 2) {
                    T s = sycl::shift_group_right(sg, y, i);
                    if (lane >= i) y += s;
                }

                if (lane < nwarps) shared2[lane] = y;
            }

            gh.item->barrier(sycl::access::fence_space::local_space);

            // shared[0:nwarps) holds the inclusive sum of warp sums.

            // Also note x still holds the inclusive sum within the
            // warp.  Given these two, we can compute the inclusive
            // sum within this chunk.
            T sum_prev_warp = (warp == 0) ? 0 : shared2[warp-1];
            T tmp_out = sum_prev_warp + sum_prev_chunk +
                (std::is_same<std::decay_t<TYPE>,Type::Inclusive>::value ? x : x-x0);
            sum_prev_chunk += shared2[nwarps-1];

            if (offset < iend) {
                blockresult_p[offset] = tmp_out;
            }
        }

        // sum_prev_chunk now holds the sum of the whole block.
        if (threadIdxx == 0) {
            blocksum_p[blockIdxx] = sum_prev_chunk;
        }
    });

    amrex::launch(1, nthreads, sm, stream,
    [=] AMREX_GPU_DEVICE (Gpu::Handler const& gh) noexcept
    {
        sycl::sub_group const& sg = gh.item->get_sub_group();
        int lane = sg.get_local_id()[0];
        int warp = sg.get_group_id()[0];
        int nwarps = sg.get_group_range()[0];

        int threadIdxx = gh.item->get_local_id(0);
        int blockDimx = gh.item->get_local_range(0);

        T* shared = (T*)(gh.local);
        T* shared2 = shared + Gpu::Device::warp_size;

        T sum_prev_chunk = 0;
        for (int offset = threadIdxx; offset - threadIdxx < nblocks; offset += blockDimx) {
            T x = (offset < nblocks) ? blocksum_p[offset] : 0;
            // Scan within a warp
            for (int i = 1; i <= Gpu::Device::warp_size; i *= 2) {
                T s = sycl::shift_group_right(sg, x, i);
                if (lane >= i) x += s;
            }

            // x now holds the inclusive sum within the warp.  The
            // last thread in each warp holds the inclusive sum of
            // this warp.  We will store it in shared memory.
            if (lane == Gpu::Device::warp_size - 1) {
                shared[warp] = x;
            }

            gh.item->barrier(sycl::access::fence_space::local_space);

            // The first warp will do scan on the warp sums for the
            // whole block.
            if (warp == 0) {
                T y = (lane < nwarps) ? shared[lane] : 0;
                for (int i = 1; i <= Gpu::Device::warp_size; i *= 2) {
                    T s = sycl::shift_group_right(sg, y, i);
                    if (lane >= i) y += s;
                }

                if (lane < nwarps) shared2[lane] = y;
            }

            gh.item->barrier(sycl::access::fence_space::local_space);

            // shared[0:nwarps) holds the inclusive sum of warp sums.

            // Also note x still holds the inclusive sum within the
            // warp.  Given these two, we can compute the inclusive
            // sum within this chunk.
            T sum_prev_warp = (warp == 0) ? 0 : shared2[warp-1];
            T tmp_out = sum_prev_warp + sum_prev_chunk + x;
            sum_prev_chunk += shared2[nwarps-1];

            if (offset < nblocks) {
                blocksum_p[offset] = tmp_out;
            }
        }

        // sum_prev_chunk now holds the total sum.
        if (threadIdxx == 0) {
            *totalsum_p = sum_prev_chunk;
        }
    });

    amrex::launch(nblocks, nthreads, 0, stream,
    [=] AMREX_GPU_DEVICE (Gpu::Handler const& gh) noexcept
    {
        int threadIdxx = gh.item->get_local_id(0);
        int blockIdxx = gh.item->get_group_linear_id();
        int blockDimx = gh.item->get_local_range(0);

        // Each block processes [ibegin,iend).
        N ibegin = nelms_per_block * blockIdxx;
        N iend = amrex::min(static_cast<N>(ibegin+nelms_per_block), n);
        T prev_sum = (blockIdxx == 0) ? 0 : blocksum_p[blockIdxx-1];
        for (N offset = ibegin + threadIdxx; offset < iend; offset += blockDimx) {
            fout(offset, prev_sum + blockresult_p[offset]);
        }
    });

    T totalsum = 0;
    if (a_ret_sum) {
        Gpu::dtoh_memcpy_async(&totalsum, totalsum_p, sizeof(T));
    }
    Gpu::streamSynchronize();
    The_Arena()->free(dp);

    AMREX_GPU_ERROR_CHECK();

    return totalsum;
}
#endif

template <typename T, typename N, typename FIN, typename FOUT, typename TYPE,
          typename M=std::enable_if_t<std::is_integral<N>::value &&
                                      (std::is_same<std::decay_t<TYPE>,Type::Inclusive>::value ||
                                       std::is_same<std::decay_t<TYPE>,Type::Exclusive>::value)> >
T PrefixSum (N n, FIN && fin, FOUT && fout, TYPE type, RetSum a_ret_sum = retSum)
{
    if (n <= 0) return 0;
    constexpr int nwarps_per_block = 8;
    constexpr int nthreads = nwarps_per_block*Gpu::Device::warp_size;
    constexpr int nchunks = 12;
    constexpr int nelms_per_block = nthreads * nchunks;
    AMREX_ALWAYS_ASSERT(static_cast<Long>(n) < static_cast<Long>(std::numeric_limits<int>::max())*nelms_per_block);
    int nblocks = (static_cast<Long>(n) + nelms_per_block - 1) / nelms_per_block;

#ifndef AMREX_SYCL_NO_MULTIPASS_SCAN
    if (nblocks > 1) {
        return PrefixSum_mp<T>(n, std::forward<FIN>(fin), std::forward<FOUT>(fout), type, retSum);
    }
#endif

    std::size_t sm = sizeof(T) * (Gpu::Device::warp_size + nwarps_per_block) + sizeof(int);
    auto stream = Gpu::gpuStream();

    typedef typename std::conditional<sizeof(detail::STVA<T>) <= 8,
        detail::BlockStatus<T,true>, detail::BlockStatus<T,false> >::type BlockStatusT;

    std::size_t nbytes_blockstatus = Arena::align(sizeof(BlockStatusT)*nblocks);
    std::size_t nbytes_blockid = Arena::align(sizeof(unsigned int));
    std::size_t nbytes_totalsum = Arena::align(sizeof(T));
    auto dp = (char*)(The_Arena()->alloc(  nbytes_blockstatus
                                                + nbytes_blockid
                                                + nbytes_totalsum));
    BlockStatusT* AMREX_RESTRICT block_status_p = (BlockStatusT*)dp;
    unsigned int* AMREX_RESTRICT virtual_block_id_p = (unsigned int*)(dp + nbytes_blockstatus);
    T* AMREX_RESTRICT totalsum_p = (T*)(dp + nbytes_blockstatus + nbytes_blockid);

    amrex::ParallelFor(nblocks, [=] AMREX_GPU_DEVICE (int i) noexcept {
        BlockStatusT& block_status = block_status_p[i];
        block_status.set_status('x');
        if (i == 0) {
            *virtual_block_id_p = 0;
            *totalsum_p = 0;
        }
    });

    amrex::launch(nblocks, nthreads, sm, stream,
    [=] AMREX_GPU_DEVICE (Gpu::Handler const& gh) noexcept
    {
        sycl::sub_group const& sg = gh.item->get_sub_group();
        int lane = sg.get_local_id()[0];
        int warp = sg.get_group_id()[0];
        int nwarps = sg.get_group_range()[0];

        int threadIdxx = gh.item->get_local_id(0);
        int blockDimx = gh.item->get_local_range(0);
        int gridDimx = gh.item->get_group_range(0);

        T* shared = (T*)(gh.local);
        T* shared2 = shared + Gpu::Device::warp_size;

        // First of all, get block virtual id.  We must do this to
        // avoid deadlock because blocks may be launched in any order.
        // Anywhere in this function, we should not use blockIdx.
        int virtual_block_id = 0;
        if (gridDimx > 1) {
            int& virtual_block_id_shared = *((int*)(shared2+nwarps));
            if (threadIdxx == 0) {
                unsigned int bid = Gpu::Atomic::Add(virtual_block_id_p, 1u);
                virtual_block_id_shared = bid;
            }
            gh.item->barrier(sycl::access::fence_space::local_space);
            virtual_block_id = virtual_block_id_shared;
        }

        // Each block processes [ibegin,iend).
        N ibegin = nelms_per_block * virtual_block_id;
        N iend = amrex::min(static_cast<N>(ibegin+nelms_per_block), n);
        BlockStatusT& block_status = block_status_p[virtual_block_id];

        //
        // The overall algorithm is based on "Single-pass Parallel
        // Prefix Scan with Decoupled Look-back" by D. Merrill &
        // M. Garland.
        //

        // Each block is responsible for nchunks chunks of data,
        // where each chunk has blockDim.x elements, one for each
        // thread in the block.
        T sum_prev_chunk = 0; // inclusive sum from previous chunks.
        T tmp_out[nchunks]; // block-wide inclusive sum for chunks
        for (int ichunk = 0; ichunk < nchunks; ++ichunk) {
            N offset = ibegin + ichunk*blockDimx;
            if (offset >= iend) break;

            offset += threadIdxx;
            T x0 = (offset < iend) ? fin(offset) : 0;
            if  (std::is_same<std::decay_t<TYPE>,Type::Exclusive>::value && offset == n-1) {
                *totalsum_p += x0;
            }
            T x = x0;
            // Scan within a warp
            for (int i = 1; i <= Gpu::Device::warp_size; i *= 2) {
                T s = sycl::shift_group_right(sg, x, i);
                if (lane >= i) x += s;
            }

            // x now holds the inclusive sum within the warp.  The
            // last thread in each warp holds the inclusive sum of
            // this warp.  We will store it in shared memory.
            if (lane == Gpu::Device::warp_size - 1) {
                shared[warp] = x;
            }

            gh.item->barrier(sycl::access::fence_space::local_space);

            // The first warp will do scan on the warp sums for the
            // whole block.
            if (warp == 0) {
                T y = (lane < nwarps) ? shared[lane] : 0;
                for (int i = 1; i <= Gpu::Device::warp_size; i *= 2) {
                    T s = sycl::shift_group_right(sg, y, i);
                    if (lane >= i) y += s;
                }

                if (lane < nwarps) shared2[lane] = y;
            }

            gh.item->barrier(sycl::access::fence_space::local_space);

            // shared[0:nwarps) holds the inclusive sum of warp sums.

            // Also note x still holds the inclusive sum within the
            // warp.  Given these two, we can compute the inclusive
            // sum within this chunk.
            T sum_prev_warp = (warp == 0) ? 0 : shared2[warp-1];
            tmp_out[ichunk] = sum_prev_warp + sum_prev_chunk +
                (std::is_same<std::decay_t<TYPE>,Type::Inclusive>::value ? x : x-x0);
            sum_prev_chunk += shared2[nwarps-1];
        }

        // sum_prev_chunk now holds the sum of the whole block.
        if (threadIdxx == 0 && gridDimx > 1) {
            block_status.write((virtual_block_id == 0) ? 'p' : 'a',
                               sum_prev_chunk);
        }

        if (virtual_block_id == 0) {
            for (int ichunk = 0; ichunk < nchunks; ++ichunk) {
                N offset = ibegin + ichunk*blockDimx + threadIdxx;
                if (offset >= iend) break;
                fout(offset, tmp_out[ichunk]);
                if (offset == n-1) {
                    *totalsum_p += tmp_out[ichunk];
                }
            }
        } else if (virtual_block_id > 0) {

            if (warp == 0) {
                T exclusive_prefix = 0;
                BlockStatusT volatile* pbs = block_status_p;
                for (int iblock0 = virtual_block_id-1; iblock0 >= 0; iblock0 -= Gpu::Device::warp_size)
                {
                    int iblock = iblock0-lane;
                    detail::STVA<T> stva{'p', 0};
                    if (iblock >= 0) {
                        stva = pbs[iblock].wait();
                    }

                    T x = stva.value;

                    // implement our own __ballot
                    unsigned status_bf = (stva.status == 'p') ? (0x1u << lane) : 0;
                    for (int i = 1; i < Gpu::Device::warp_size; i *= 2) {
                        status_bf |= sycl::permute_group_by_xor(sg, status_bf, i);
                    }

                    bool stop_lookback = status_bf & 0x1u;
                    if (stop_lookback == false) {
                        if (status_bf != 0) {
                            T y = x;
                            if (lane > 0) x = 0;
                            unsigned int bit_mask = 0x1u;
                            for (int i = 1; i < Gpu::Device::warp_size; ++i) {
                                bit_mask <<= 1;
                                if (i == lane) x = y;
                                if (status_bf & bit_mask) {
                                    stop_lookback = true;
                                    break;
                                }
                            }
                        }

                        for (int i = Gpu::Device::warp_size/2; i > 0; i /= 2) {
                            x += sycl::shift_group_left(sg, x,i);
                        }
                    }

                    if (lane == 0) { exclusive_prefix += x; }
                    if (stop_lookback) break;
                }

                if (lane == 0) {
                    block_status.write('p', block_status.get_aggregate() + exclusive_prefix);
                    shared[0] = exclusive_prefix;
                }
            }

            gh.item->barrier(sycl::access::fence_space::local_space);

            T exclusive_prefix = shared[0];

            for (int ichunk = 0; ichunk < nchunks; ++ichunk) {
                N offset = ibegin + ichunk*blockDimx + threadIdxx;
                if (offset >= iend) break;
                T t = tmp_out[ichunk] + exclusive_prefix;
                fout(offset, t);
                if (offset == n-1) {
                    *totalsum_p += t;
                }
            }
        }
    });

    T totalsum = 0;
    if (a_ret_sum) {
        // xxxxx SYCL todo: Should test if using pinned memory and thus
        // avoiding memcpy is faster.
        Gpu::dtoh_memcpy_async(&totalsum, totalsum_p, sizeof(T));
    }
    Gpu::streamSynchronize();
    The_Arena()->free(dp);

    AMREX_GPU_ERROR_CHECK();

    return totalsum;
}

#elif defined(AMREX_USE_HIP)

template <typename T, typename N, typename FIN, typename FOUT, typename TYPE,
          typename M=std::enable_if_t<std::is_integral<N>::value &&
                                      (std::is_same<std::decay_t<TYPE>,Type::Inclusive>::value ||
                                       std::is_same<std::decay_t<TYPE>,Type::Exclusive>::value)> >
T PrefixSum (N n, FIN && fin, FOUT && fout, TYPE, RetSum a_ret_sum = retSum)
{
    if (n <= 0) return 0;
    constexpr int nwarps_per_block = 4;
    constexpr int nthreads = nwarps_per_block*Gpu::Device::warp_size; // # of threads per block
    constexpr int nelms_per_thread = sizeof(T) >= 8 ? 8 : 16;
    constexpr int nelms_per_block = nthreads * nelms_per_thread;
    int nblocks = (n + nelms_per_block - 1) / nelms_per_block;
    std::size_t sm = 0;
    auto stream = Gpu::gpuStream();

    typedef rocprim::detail::lookback_scan_state<T> ScanTileState;
    typedef rocprim::detail::ordered_block_id<unsigned int> OrderedBlockId;

    std::size_t nbytes_tile_state = rocprim::detail::align_size
        (ScanTileState::get_storage_size(nblocks));
    std::size_t nbytes_block_id = OrderedBlockId::get_storage_size();

    auto dp = (char*)(The_Arena()->alloc(nbytes_tile_state+nbytes_block_id));

    ScanTileState tile_state = ScanTileState::create(dp, nblocks);
    auto ordered_block_id = OrderedBlockId::create
        (reinterpret_cast<OrderedBlockId::id_type*>(dp + nbytes_tile_state));

    // Init ScanTileState on device
    amrex::launch((nblocks+nthreads-1)/nthreads, nthreads, 0, stream, [=] AMREX_GPU_DEVICE ()
    {
        auto& scan_tile_state = const_cast<ScanTileState&>(tile_state);
        auto& scan_bid = const_cast<OrderedBlockId&>(ordered_block_id);
        const unsigned int gid = blockIdx.x*blockDim.x + threadIdx.x;
        if (gid == 0) { scan_bid.reset(); }
        scan_tile_state.initialize_prefix(gid, nblocks);
    });

    T* totalsum_p = (a_ret_sum) ? (T*)(The_Pinned_Arena()->alloc(sizeof(T))) : nullptr;

    amrex::launch_global<nthreads> <<<nblocks, nthreads, sm, stream>>> (
    [=] AMREX_GPU_DEVICE () noexcept
    {
        typedef rocprim::block_load<T, nthreads, nelms_per_thread,
                                    rocprim::block_load_method::block_load_transpose> BlockLoad;
        typedef rocprim::block_scan<T, nthreads,
                                    rocprim::block_scan_algorithm::using_warp_scan> BlockScan;
        typedef rocprim::block_exchange<T, nthreads, nelms_per_thread> BlockExchange;
        typedef rocprim::detail::lookback_scan_prefix_op
            <T, rocprim::plus<T>, ScanTileState> LookbackScanPrefixOp;

        __shared__ struct TempStorage {
            typename OrderedBlockId::storage_type ordered_bid;
            union {
                typename BlockLoad::storage_type     load;
                typename BlockExchange::storage_type exchange;
                typename BlockScan::storage_type     scan;
            };
        } temp_storage;

        // Lambda captured tile_state is const.  We have to cast the const away.
        auto& scan_tile_state = const_cast<ScanTileState&>(tile_state);
        auto& scan_bid = const_cast<OrderedBlockId&>(ordered_block_id);

        auto const virtual_block_id = scan_bid.get(threadIdx.x, temp_storage.ordered_bid);

        // Each block processes [ibegin,iend).
        N ibegin = nelms_per_block * virtual_block_id;
        N iend = amrex::min(static_cast<N>(ibegin+nelms_per_block), n);

        auto input_begin = rocprim::make_transform_iterator(
            rocprim::make_counting_iterator(N(0)),
            [&] (N i) -> T { return fin(i+ibegin); });

        T data[nelms_per_thread];
        if (static_cast<int>(iend-ibegin) == nelms_per_block) {
            BlockLoad().load(input_begin, data, temp_storage.load);
        } else {
            // padding with 0
            BlockLoad().load(input_begin, data, iend-ibegin, 0, temp_storage.load);
        }

        __syncthreads();

        constexpr bool is_exclusive = std::is_same<std::decay_t<TYPE>,Type::Exclusive>::value;

        if (virtual_block_id == 0) {
            T block_agg;
            AMREX_IF_CONSTEXPR(is_exclusive) {
                BlockScan().exclusive_scan(data, data, T{0}, block_agg, temp_storage.scan);
            } else {
                BlockScan().inclusive_scan(data, data, block_agg, temp_storage.scan);
            }
            if (threadIdx.x == 0) {
                if (nblocks > 1) {
                    scan_tile_state.set_complete(0, block_agg);
                } else if (nblocks == 1 && totalsum_p) {
                    *totalsum_p = block_agg;
                }
            }
        } else {
            T last = data[nelms_per_thread-1]; // Need this for the total sum in exclusive case

            LookbackScanPrefixOp prefix_op(virtual_block_id, rocprim::plus<T>(), scan_tile_state);
            AMREX_IF_CONSTEXPR(is_exclusive) {
                BlockScan().exclusive_scan(data, data, temp_storage.scan, prefix_op,
                                           rocprim::plus<T>());
            } else {
                BlockScan().inclusive_scan(data, data, temp_storage.scan, prefix_op,
                                           rocprim::plus<T>());
            }
            if (totalsum_p) {
                if (iend == n && threadIdx.x == blockDim.x-1) { // last thread of last block
                    T tsum = data[nelms_per_thread-1];
                    AMREX_IF_CONSTEXPR(is_exclusive) { tsum += last; }
                    *totalsum_p = tsum;
                }
            }
        }

        __syncthreads();

        BlockExchange().blocked_to_striped(data, data, temp_storage.exchange);

        for (int i = 0; i < nelms_per_thread; ++i) {
            N offset = ibegin + i*blockDim.x + threadIdx.x;
            if (offset < iend) { fout(offset, data[i]); }
        }
    });

    Gpu::streamSynchronize();
    AMREX_GPU_ERROR_CHECK();

    The_Arena()->free(dp);

    T ret = (a_ret_sum) ? *totalsum_p : T(0);
    if (totalsum_p) { The_Pinned_Arena()->free(totalsum_p); }

    return ret;
}

#elif defined(AMREX_USE_CUDA) && defined(__CUDACC__) && (__CUDACC_VER_MAJOR__ >= 11)

template <typename T, typename N, typename FIN, typename FOUT, typename TYPE,
          typename M=std::enable_if_t<std::is_integral<N>::value &&
                                      (std::is_same<std::decay_t<TYPE>,Type::Inclusive>::value ||
                                       std::is_same<std::decay_t<TYPE>,Type::Exclusive>::value)> >
T PrefixSum (N n, FIN && fin, FOUT && fout, TYPE, RetSum a_ret_sum = retSum)
{
    if (n <= 0) return 0;
    constexpr int nwarps_per_block = 8;
    constexpr int nthreads = nwarps_per_block*Gpu::Device::warp_size; // # of threads per block
    constexpr int nelms_per_thread = sizeof(T) >= 8 ? 4 : 8;
    constexpr int nelms_per_block = nthreads * nelms_per_thread;
    int nblocks = (n + nelms_per_block - 1) / nelms_per_block;
    std::size_t sm = 0;
    auto stream = Gpu::gpuStream();

    typedef cub::ScanTileState<T> ScanTileState;
    std::size_t tile_state_size = 0;
    ScanTileState::AllocationSize(nblocks, tile_state_size);

    std::size_t nbytes_tile_state = Arena::align(tile_state_size);
    auto tile_state_p = (char*)(The_Arena()->alloc(nbytes_tile_state));

    ScanTileState tile_state;
    tile_state.Init(nblocks, tile_state_p, tile_state_size); // Init ScanTileState on host

    if (nblocks > 1) {
        // Init ScanTileState on device
        amrex::launch((nblocks+nthreads-1)/nthreads, nthreads, 0, stream, [=] AMREX_GPU_DEVICE ()
        {
            const_cast<ScanTileState&>(tile_state).InitializeStatus(nblocks);
        });
    }

    T* totalsum_p = (a_ret_sum) ? (T*)(The_Pinned_Arena()->alloc(sizeof(T))) : nullptr;

    amrex::launch_global<nthreads> <<<nblocks, nthreads, sm, stream>>> (
    [=] AMREX_GPU_DEVICE () noexcept
    {
        typedef cub::BlockLoad<T, nthreads, nelms_per_thread, cub::BLOCK_LOAD_WARP_TRANSPOSE>
            BlockLoad;
        typedef cub::BlockScan<T, nthreads, cub::BLOCK_SCAN_WARP_SCANS> BlockScan;
        typedef cub::BlockExchange<T, nthreads, nelms_per_thread>       BlockExchange;
        typedef cub::TilePrefixCallbackOp<T, cub::Sum, ScanTileState>   TilePrefixCallbackOp;

        __shared__ union TempStorage
        {
            typename BlockLoad::TempStorage     load;
            typename BlockExchange::TempStorage exchange;
            struct ScanStorage {
                typename BlockScan::TempStorage            scan;
                typename TilePrefixCallbackOp::TempStorage prefix;
            } scan_storeage;
        } temp_storage;

        // Lambda captured tile_state is const.  We have to cast the const away.
        auto& scan_tile_state = const_cast<ScanTileState&>(tile_state);

        int virtual_block_id = blockIdx.x;

        // Each block processes [ibegin,iend).
        N ibegin = nelms_per_block * virtual_block_id;
        N iend = amrex::min(static_cast<N>(ibegin+nelms_per_block), n);

        auto input_lambda = [&] (N i) -> T { return fin(i+ibegin); };
        cub::TransformInputIterator<T,decltype(input_lambda),cub::CountingInputIterator<N> >
            input_begin(cub::CountingInputIterator<N>(0), input_lambda);

        T data[nelms_per_thread];
        if (static_cast<int>(iend-ibegin) == nelms_per_block) {
            BlockLoad(temp_storage.load).Load(input_begin, data);
        } else {
            BlockLoad(temp_storage.load).Load(input_begin, data, iend-ibegin, 0); // padding with 0
        }

        __syncthreads();

        constexpr bool is_exclusive = std::is_same<std::decay_t<TYPE>,Type::Exclusive>::value;

        if (virtual_block_id == 0) {
            T block_agg;
            AMREX_IF_CONSTEXPR(is_exclusive) {
                BlockScan(temp_storage.scan_storeage.scan).ExclusiveSum(data, data, block_agg);
            } else {
                BlockScan(temp_storage.scan_storeage.scan).InclusiveSum(data, data, block_agg);
            }
            if (threadIdx.x == 0) {
                if (nblocks > 1) {
                    scan_tile_state.SetInclusive(0, block_agg);
                } else if (nblocks == 1 && totalsum_p) {
                    *totalsum_p = block_agg;
                }
            }
        } else {
            T last = data[nelms_per_thread-1]; // Need this for the total sum in exclusive case

            TilePrefixCallbackOp prefix_op(scan_tile_state, temp_storage.scan_storeage.prefix,
                                           cub::Sum{}, virtual_block_id);
            AMREX_IF_CONSTEXPR(is_exclusive) {
                BlockScan(temp_storage.scan_storeage.scan).ExclusiveSum(data, data, prefix_op);
            } else {
                BlockScan(temp_storage.scan_storeage.scan).InclusiveSum(data, data, prefix_op);
            }
            if (totalsum_p) {
                if (iend == n && threadIdx.x == blockDim.x-1) { // last thread of last block
                    T tsum = data[nelms_per_thread-1];
                    AMREX_IF_CONSTEXPR(is_exclusive) { tsum += last; }
                    *totalsum_p = tsum;
                }
            }
        }

        __syncthreads();

        BlockExchange(temp_storage.exchange).BlockedToStriped(data);

        for (int i = 0; i < nelms_per_thread; ++i) {
            N offset = ibegin + i*blockDim.x + threadIdx.x;
            if (offset < iend) { fout(offset, data[i]); }
        }
    });

    Gpu::streamSynchronize();
    AMREX_GPU_ERROR_CHECK();

    The_Arena()->free(tile_state_p);

    T ret = (a_ret_sum) ? *totalsum_p : T(0);
    if (totalsum_p) { The_Pinned_Arena()->free(totalsum_p); }

    return ret;
}

#else

template <typename T, typename N, typename FIN, typename FOUT, typename TYPE,
          typename M=std::enable_if_t<std::is_integral<N>::value &&
                                      (std::is_same<std::decay_t<TYPE>,Type::Inclusive>::value ||
                                       std::is_same<std::decay_t<TYPE>,Type::Exclusive>::value)> >
T PrefixSum (N n, FIN && fin, FOUT && fout, TYPE, RetSum a_ret_sum = retSum)
{
    if (n <= 0) return 0;
    constexpr int nwarps_per_block = 4;
    constexpr int nthreads = nwarps_per_block*Gpu::Device::warp_size;
    constexpr int nchunks = 12;
    constexpr int nelms_per_block = nthreads * nchunks;
    AMREX_ALWAYS_ASSERT(static_cast<Long>(n) < static_cast<Long>(std::numeric_limits<int>::max())*nelms_per_block);
    int nblocks = (static_cast<Long>(n) + nelms_per_block - 1) / nelms_per_block;
    std::size_t sm = sizeof(T) * (Gpu::Device::warp_size + nwarps_per_block) + sizeof(int);
    auto stream = Gpu::gpuStream();

    typedef typename std::conditional<sizeof(detail::STVA<T>) <= 8,
        detail::BlockStatus<T,true>, detail::BlockStatus<T,false> >::type BlockStatusT;

    std::size_t nbytes_blockstatus = Arena::align(sizeof(BlockStatusT)*nblocks);
    std::size_t nbytes_blockid = Arena::align(sizeof(unsigned int));
    std::size_t nbytes_totalsum = Arena::align(sizeof(T));
    auto dp = (char*)(The_Arena()->alloc(  nbytes_blockstatus
                                         + nbytes_blockid
                                         + nbytes_totalsum));
    BlockStatusT* AMREX_RESTRICT block_status_p = (BlockStatusT*)dp;
    unsigned int* AMREX_RESTRICT virtual_block_id_p = (unsigned int*)(dp + nbytes_blockstatus);
    T* AMREX_RESTRICT totalsum_p = (T*)(dp + nbytes_blockstatus + nbytes_blockid);

    amrex::ParallelFor(nblocks, [=] AMREX_GPU_DEVICE (int i) noexcept {
        BlockStatusT& block_status = block_status_p[i];
        block_status.set_status('x');
        if (i == 0) {
            *virtual_block_id_p = 0;
            *totalsum_p = 0;
        }
    });

    amrex::launch(nblocks, nthreads, sm, stream,
    [=] AMREX_GPU_DEVICE () noexcept
    {
        int lane = threadIdx.x % Gpu::Device::warp_size;
        int warp = threadIdx.x / Gpu::Device::warp_size;
        int nwarps = blockDim.x / Gpu::Device::warp_size;

        amrex::Gpu::SharedMemory<T> gsm;
        T* shared = gsm.dataPtr();
        T* shared2 = shared + Gpu::Device::warp_size;

        // First of all, get block virtual id.  We must do this to
        // avoid deadlock because CUDA may launch blocks in any order.
        // Anywhere in this function, we should not use blockIdx.
        int virtual_block_id = 0;
        if (gridDim.x > 1) {
            int& virtual_block_id_shared = *((int*)(shared2+nwarps));
            if (threadIdx.x == 0) {
                unsigned int bid = Gpu::Atomic::Add(virtual_block_id_p, 1u);
                virtual_block_id_shared = bid;
            }
            __syncthreads();
            virtual_block_id = virtual_block_id_shared;
        }

        // Each block processes [ibegin,iend).
        N ibegin = nelms_per_block * virtual_block_id;
        N iend = amrex::min(static_cast<N>(ibegin+nelms_per_block), n);
        BlockStatusT& block_status = block_status_p[virtual_block_id];

        //
        // The overall algorithm is based on "Single-pass Parallel
        // Prefix Scan with Decoupled Look-back" by D. Merrill &
        // M. Garland.
        //

        // Each block is responsible for nchunks chunks of data,
        // where each chunk has blockDim.x elements, one for each
        // thread in the block.
        T sum_prev_chunk = 0; // inclusive sum from previous chunks.
        T tmp_out[nchunks]; // block-wide inclusive sum for chunks
        for (int ichunk = 0; ichunk < nchunks; ++ichunk) {
            N offset = ibegin + ichunk*blockDim.x;
            if (offset >= iend) break;

            offset += threadIdx.x;
            T x0 = (offset < iend) ? fin(offset) : 0;
            if  (std::is_same<std::decay_t<TYPE>,Type::Exclusive>::value && offset == n-1) {
                *totalsum_p += x0;
            }
            T x = x0;
            // Scan within a warp
            for (int i = 1; i <= Gpu::Device::warp_size; i *= 2) {
                AMREX_HIP_OR_CUDA( T s = __shfl_up(x,i);,
                                   T s = __shfl_up_sync(0xffffffff, x, i); )
                if (lane >= i) x += s;
            }

            // x now holds the inclusive sum within the warp.  The
            // last thread in each warp holds the inclusive sum of
            // this warp.  We will store it in shared memory.
            if (lane == Gpu::Device::warp_size - 1) {
                shared[warp] = x;
            }

            __syncthreads();

            // The first warp will do scan on the warp sums for the
            // whole block.  Not all threads in the warp need to
            // participate.
#ifdef AMREX_USE_CUDA
            if (warp == 0 && lane < nwarps) {
                T y = shared[lane];
                int mask = (1 << nwarps) - 1;
                for (int i = 1; i <= nwarps; i *= 2) {
                    T s = __shfl_up_sync(mask, y, i, nwarps);
                    if (lane >= i) y += s;
                }
                shared2[lane] = y;
            }
#else
            if (warp == 0) {
                T y = 0;
                if (lane < nwarps) {
                    y = shared[lane];
                }
                for (int i = 1; i <= nwarps; i *= 2) {
                    T s = __shfl_up(y, i, nwarps);
                    if (lane >= i) y += s;
                }
                if (lane < nwarps) {
                    shared2[lane] = y;
                }
            }
#endif

            __syncthreads();

            // shared[0:nwarps) holds the inclusive sum of warp sums.

            // Also note x still holds the inclusive sum within the
            // warp.  Given these two, we can compute the inclusive
            // sum within this chunk.
            T sum_prev_warp = (warp == 0) ? 0 : shared2[warp-1];
            tmp_out[ichunk] = sum_prev_warp + sum_prev_chunk +
                (std::is_same<std::decay_t<TYPE>,Type::Inclusive>::value ? x : x-x0);
            sum_prev_chunk += shared2[nwarps-1];
        }

        // sum_prev_chunk now holds the sum of the whole block.
        if (threadIdx.x == 0 && gridDim.x > 1) {
            block_status.write((virtual_block_id == 0) ? 'p' : 'a',
                               sum_prev_chunk);
        }

        if (virtual_block_id == 0) {
            for (int ichunk = 0; ichunk < nchunks; ++ichunk) {
                N offset = ibegin + ichunk*blockDim.x + threadIdx.x;
                if (offset >= iend) break;
                fout(offset, tmp_out[ichunk]);
                if (offset == n-1) {
                    *totalsum_p += tmp_out[ichunk];
                }
            }
        } else if (virtual_block_id > 0) {

            if (warp == 0) {
                T exclusive_prefix = 0;
                BlockStatusT volatile* pbs = block_status_p;
                for (int iblock0 = virtual_block_id-1; iblock0 >= 0; iblock0 -= Gpu::Device::warp_size)
                {
                    int iblock = iblock0-lane;
                    detail::STVA<T> stva{'p', 0};
                    if (iblock >= 0) {
                        stva = pbs[iblock].wait();
                    }

                    T x = stva.value;

                    AMREX_HIP_OR_CUDA( uint64_t const status_bf = __ballot(stva.status == 'p');,
                                       unsigned const status_bf = __ballot_sync(0xffffffff, stva.status == 'p'));
                    bool stop_lookback = status_bf & 0x1u;
                    if (stop_lookback == false) {
                        if (status_bf != 0) {
                            T y = x;
                            if (lane > 0) x = 0;
                            AMREX_HIP_OR_CUDA(uint64_t bit_mask = 0x1ull;,
                                              unsigned bit_mask = 0x1u);
                            for (int i = 1; i < Gpu::Device::warp_size; ++i) {
                                bit_mask <<= 1;
                                if (i == lane) x = y;
                                if (status_bf & bit_mask) {
                                    stop_lookback = true;
                                    break;
                                }
                            }
                        }

                        for (int i = Gpu::Device::warp_size/2; i > 0; i /= 2) {
                            AMREX_HIP_OR_CUDA( x += __shfl_down(x,i);,
                                               x += __shfl_down_sync(0xffffffff, x, i); )
                        }
                    }

                    if (lane == 0) { exclusive_prefix += x; }
                    if (stop_lookback) break;
                }

                if (lane == 0) {
                    block_status.write('p', block_status.get_aggregate() + exclusive_prefix);
                    shared[0] = exclusive_prefix;
                }
            }

            __syncthreads();

            T exclusive_prefix = shared[0];

            for (int ichunk = 0; ichunk < nchunks; ++ichunk) {
                N offset = ibegin + ichunk*blockDim.x + threadIdx.x;
                if (offset >= iend) break;
                T t = tmp_out[ichunk] + exclusive_prefix;
                fout(offset, t);
                if (offset == n-1) {
                    *totalsum_p += t;
                }
            }
        }
    });

    T totalsum = 0;
    if (a_ret_sum) {
        // xxxxx CUDA < 11 todo: Should test if using pinned memory and thus
        // avoiding memcpy is faster.
        Gpu::dtoh_memcpy_async(&totalsum, totalsum_p, sizeof(T));
    }
    Gpu::streamSynchronize();
    The_Arena()->free(dp);

    AMREX_GPU_ERROR_CHECK();

    return totalsum;
}

#endif

// The return value is the total sum if a_ret_sum is true.
template <typename N, typename T, typename M=std::enable_if_t<std::is_integral<N>::value> >
T InclusiveSum (N n, T const* in, T * out, RetSum a_ret_sum = retSum)
{
#if defined(AMREX_USE_CUDA) && defined(__CUDACC__) && (__CUDACC_VER_MAJOR__ >= 11)
    void* d_temp = nullptr;
    std::size_t temp_bytes = 0;
    AMREX_GPU_SAFE_CALL(cub::DeviceScan::InclusiveSum(d_temp, temp_bytes, in, out, n,
                                                      Gpu::gpuStream()));
    d_temp = The_Arena()->alloc(temp_bytes);
    AMREX_GPU_SAFE_CALL(cub::DeviceScan::InclusiveSum(d_temp, temp_bytes, in, out, n,
                                                      Gpu::gpuStream()));
    T totalsum = 0;
    if (a_ret_sum) {
        Gpu::dtoh_memcpy_async(&totalsum, out+(n-1), sizeof(T));
    }
    Gpu::streamSynchronize();
    The_Arena()->free(d_temp);
    AMREX_GPU_ERROR_CHECK();
    return totalsum;
#elif defined(AMREX_USE_HIP)
    void* d_temp = nullptr;
    std::size_t temp_bytes = 0;
    AMREX_GPU_SAFE_CALL(rocprim::inclusive_scan(d_temp, temp_bytes, in, out, n,
                                                rocprim::plus<T>(), Gpu::gpuStream()));
    d_temp = The_Arena()->alloc(temp_bytes);
    AMREX_GPU_SAFE_CALL(rocprim::inclusive_scan(d_temp, temp_bytes, in, out, n,
                                                rocprim::plus<T>(), Gpu::gpuStream()));
    T totalsum = 0;
    if (a_ret_sum) {
        Gpu::dtoh_memcpy_async(&totalsum, out+(n-1), sizeof(T));
    }
    Gpu::streamSynchronize();
    The_Arena()->free(d_temp);
    AMREX_GPU_ERROR_CHECK();
    return totalsum;
#elif defined(AMREX_USE_SYCL) && defined(AMREX_USE_ONEDPL)
    auto policy = oneapi::dpl::execution::make_device_policy(Gpu::Device::streamQueue());
    std::inclusive_scan(policy, in, in+n, out, std::plus<T>(), T(0));
    T totalsum = 0;
    if (a_ret_sum) {
        Gpu::dtoh_memcpy_async(&totalsum, out+(n-1), sizeof(T));
    }
    Gpu::streamSynchronize();
    AMREX_GPU_ERROR_CHECK();
    return totalsum;
#else
    if (static_cast<Long>(n) <= static_cast<Long>(std::numeric_limits<int>::max())) {
        return PrefixSum<T>(static_cast<int>(n),
                            [=] AMREX_GPU_DEVICE (N i) -> T { return in[i]; },
                            [=] AMREX_GPU_DEVICE (N i, T const& x) { out[i] = x; },
                            Type::inclusive, a_ret_sum);
    } else {
        return PrefixSum<T>(n,
                            [=] AMREX_GPU_DEVICE (N i) -> T { return in[i]; },
                            [=] AMREX_GPU_DEVICE (N i, T const& x) { out[i] = x; },
                            Type::inclusive, a_ret_sum);
    }
#endif
}

// The return value is the total sum if a_ret_sum is true.
template <typename N, typename T, typename M=std::enable_if_t<std::is_integral<N>::value> >
T ExclusiveSum (N n, T const* in, T * out, RetSum a_ret_sum = retSum)
{
    if (n <= 0) return 0;
#if defined(AMREX_USE_CUDA) && defined(__CUDACC__) && (__CUDACC_VER_MAJOR__ >= 11)
    T in_last = 0;
    if (a_ret_sum) {
        Gpu::dtoh_memcpy_async(&in_last, in+(n-1), sizeof(T));
    }
    void* d_temp = nullptr;
    std::size_t temp_bytes = 0;
    AMREX_GPU_SAFE_CALL(cub::DeviceScan::ExclusiveSum(d_temp, temp_bytes, in, out, n,
                                                      Gpu::gpuStream()));
    d_temp = The_Arena()->alloc(temp_bytes);
    AMREX_GPU_SAFE_CALL(cub::DeviceScan::ExclusiveSum(d_temp, temp_bytes, in, out, n,
                                                      Gpu::gpuStream()));
    T out_last = 0;
    if (a_ret_sum) {
        Gpu::dtoh_memcpy_async(&out_last, out+(n-1), sizeof(T));
    }
    Gpu::streamSynchronize();
    The_Arena()->free(d_temp);
    AMREX_GPU_ERROR_CHECK();
    return in_last+out_last;
#elif defined(AMREX_USE_HIP)
    T in_last = 0;
    if (a_ret_sum) {
        Gpu::dtoh_memcpy_async(&in_last, in+(n-1), sizeof(T));
    }
    void* d_temp = nullptr;
    std::size_t temp_bytes = 0;
    AMREX_GPU_SAFE_CALL(rocprim::exclusive_scan(d_temp, temp_bytes, in, out, T{0}, n,
                                                rocprim::plus<T>(), Gpu::gpuStream()));
    d_temp = The_Arena()->alloc(temp_bytes);
    AMREX_GPU_SAFE_CALL(rocprim::exclusive_scan(d_temp, temp_bytes, in, out, T{0}, n,
                                                rocprim::plus<T>(), Gpu::gpuStream()));
    T out_last = 0;
    if (a_ret_sum) {
        Gpu::dtoh_memcpy_async(&out_last, out+(n-1), sizeof(T));
    }
    Gpu::streamSynchronize();
    The_Arena()->free(d_temp);
    AMREX_GPU_ERROR_CHECK();
    return in_last+out_last;
#elif defined(AMREX_USE_SYCL) && defined(AMREX_USE_ONEDPL)
    T in_last = 0;
    if (a_ret_sum) {
        Gpu::dtoh_memcpy_async(&in_last, in+(n-1), sizeof(T));
    }
    auto policy = oneapi::dpl::execution::make_device_policy(Gpu::Device::streamQueue());
    std::exclusive_scan(policy, in, in+n, out, T(0), std::plus<T>());
    T out_last = 0;
    if (a_ret_sum) {
        Gpu::dtoh_memcpy_async(&out_last, out+(n-1), sizeof(T));
    }
    Gpu::streamSynchronize();
    AMREX_GPU_ERROR_CHECK();
    return in_last+out_last;
#else
    if (static_cast<Long>(n) <= static_cast<Long>(std::numeric_limits<int>::max())) {
        return PrefixSum<T>(static_cast<int>(n),
                            [=] AMREX_GPU_DEVICE (N i) -> T { return in[i]; },
                            [=] AMREX_GPU_DEVICE (N i, T const& x) { out[i] = x; },
                            Type::exclusive, a_ret_sum);
    } else {
        return PrefixSum<T>(n,
                            [=] AMREX_GPU_DEVICE (N i) -> T { return in[i]; },
                            [=] AMREX_GPU_DEVICE (N i, T const& x) { out[i] = x; },
                            Type::exclusive, a_ret_sum);
    }
#endif
}

#else
//  !defined(AMREX_USE_GPU)
template <typename T, typename N, typename FIN, typename FOUT, typename TYPE,
          typename M=std::enable_if_t<std::is_integral<N>::value &&
                                      (std::is_same<std::decay_t<TYPE>,Type::Inclusive>::value ||
                                       std::is_same<std::decay_t<TYPE>,Type::Exclusive>::value)> >
T PrefixSum (N n, FIN && fin, FOUT && fout, TYPE, RetSum = retSum)
{
    if (n <= 0) return 0;
    T totalsum = 0;
    for (N i = 0; i < n; ++i) {
        T x = fin(i);
        T y = totalsum;
        totalsum += x;
        AMREX_IF_CONSTEXPR (std::is_same<std::decay_t<TYPE>,Type::Inclusive>::value) {
            y += x;
        }
        fout(i, y);
    }
    return totalsum;
}

// The return value is the total sum.
template <typename N, typename T, typename M=std::enable_if_t<std::is_integral<N>::value> >
T InclusiveSum (N n, T const* in, T * out, RetSum /*a_ret_sum*/ = retSum)
{
#if (__cplusplus >= 201703L) && (!defined(_GLIBCXX_RELEASE) || _GLIBCXX_RELEASE >= 10)
    // GCC's __cplusplus is not a reliable indication for C++17 support
    std::inclusive_scan(in, in+n, out);
#else
    std::partial_sum(in, in+n, out);
#endif
    return (n > 0) ? out[n-1] : T(0);
}

// The return value is the total sum.
template <typename N, typename T, typename M=std::enable_if_t<std::is_integral<N>::value> >
T ExclusiveSum (N n, T const* in, T * out, RetSum /*a_ret_sum*/ = retSum)
{
    if (n <= 0) return 0;

    auto in_last = in[n-1];
#if (__cplusplus >= 201703L) && (!defined(_GLIBCXX_RELEASE) || _GLIBCXX_RELEASE >= 10)
    // GCC's __cplusplus is not a reliable indication for C++17 support
    std::exclusive_scan(in, in+n, out, 0);
#else
    out[0] = 0;
    std::partial_sum(in, in+n-1, out+1);
#endif
    return in_last + out[n-1];
}

#endif

}

namespace Gpu
{
    template<class InIter, class OutIter>
    OutIter inclusive_scan (InIter begin, InIter end, OutIter result)
    {
#if defined(AMREX_USE_GPU)
        auto N = std::distance(begin, end);
        Scan::InclusiveSum(N, &(*begin), &(*result), Scan::noRetSum);
        OutIter result_end = result;
        std::advance(result_end, N);
        return result_end;
#elif (__cplusplus >= 201703L) && (!defined(_GLIBCXX_RELEASE) || _GLIBCXX_RELEASE >= 10)
        // GCC's __cplusplus is not a reliable indication for C++17 support
        return std::inclusive_scan(begin, end, result);
#else
        return std::partial_sum(begin, end, result);
#endif
    }

    template<class InIter, class OutIter>
    OutIter exclusive_scan (InIter begin, InIter end, OutIter result)
    {
#if defined(AMREX_USE_GPU)
        auto N = std::distance(begin, end);
        Scan::ExclusiveSum(N, &(*begin), &(*result), Scan::noRetSum);
        OutIter result_end = result;
        std::advance(result_end, N);
        return result_end;
#elif (__cplusplus >= 201703L) && (!defined(_GLIBCXX_RELEASE) || _GLIBCXX_RELEASE >= 10)
        // GCC's __cplusplus is not a reliable indication for C++17 support
        return std::exclusive_scan(begin, end, result, 0);
#else
        if (begin == end) return result;

        typename std::iterator_traits<InIter>::value_type sum = *begin;
        *result++ = sum - *begin;

        while (++begin != end) {
            sum = std::move(sum) + *begin;
            *result++ = sum - *begin;
        }
        return ++result;
#endif
    }

}}

#endif
