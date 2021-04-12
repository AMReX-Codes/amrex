#ifndef AMREX_SCAN_H_
#define AMREX_SCAN_H_
#include <AMReX_Config.H>

#include <AMReX_Gpu.H>
#include <AMReX_Arena.H>

#if defined(AMREX_USE_CUDA) && defined(__CUDACC__) && (__CUDACC_VER_MAJOR__ >= 11)
#  include <cub/cub.cuh>
#elif defined(AMREX_USE_HIP)
#  include <rocprim/rocprim.hpp>
#endif

#include <cstdint>
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
                      "HIP/DPCPP: unsigned long long must be 64 bits");
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
                      "HIP/DPCPP: unsigned long long must be 64 bits");
        unsigned long long tmp = Gpu::Atomic::Add
            (reinterpret_cast<unsigned long long*>(const_cast<Data<T>*>(&d)), 0ull);
        return (*reinterpret_cast<Data<T>*>(&tmp)).s;
#endif
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void set_status (char a_status) { d.s.status = a_status; }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
#if defined(AMREX_USE_DPCPP)
    STVA<T> wait (sycl::nd_item<1> const& item) volatile {
#else
    STVA<T> wait () volatile {
#endif
        STVA<T> r;
        do {
#if defined(AMREX_USE_DPCPP)
            item.mem_fence();
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
#if defined(AMREX_USE_DPCPP)
    void write (char a_status, T a_value, sycl::nd_item<1> const& item) {
#else
    void write (char a_status, T a_value) {
#endif
        if (a_status == 'a') {
            aggregate = a_value;
        } else {
            inclusive = a_value;
        }
#if defined(AMREX_USE_DPCPP)
        item.mem_fence(); // xxxxx DPCPP todo: This is at block level, but needs to be device level fence, which is currently a PR in intel/llvm
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
#if defined(AMREX_USE_DPCPP)
    STVA<T> wait (sycl::nd_item<1> const& item) volatile {
#else
    STVA<T> wait () volatile {
#endif
        STVA<T> r;
        do {
            r = read();
#if defined(AMREX_USE_DPCPP)
            item.mem_fence(); // xxxxx DPCPP todo: This is at block level, but needs to be device level fence, which is currently a PR in intel/llvm
#else
            __threadfence();
#endif
        } while (r.status == 'x');
        return r;
    }
};

}

#if defined(AMREX_USE_DPCPP)

template <typename T, typename N, typename FIN, typename FOUT, typename TYPE,
          typename M=std::enable_if_t<std::is_integral<N>::value &&
                                      (std::is_same<std::decay_t<TYPE>,Type::Inclusive>::value ||
                                       std::is_same<std::decay_t<TYPE>,Type::Exclusive>::value)> >
T PrefixSum (N n, FIN && fin, FOUT && fout, TYPE, RetSum a_ret_sum = retSum)
{
    if (n <= 0) return 0;
    constexpr int nwarps_per_block = 8;
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
    [=] AMREX_GPU_DEVICE (Gpu::Handler const& gh) noexcept
    {
        sycl::ONEAPI::sub_group const& sg = gh.item->get_sub_group();
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
                unsigned int bid = Gpu::Atomic::Inc<sycl::access::address_space::global_space>
                    (virtual_block_id_p, gridDimx);
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
                T s = sg.shuffle_up(x, i);
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
                    T s = sg.shuffle_up(y, i);
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
                        stva = pbs[iblock].wait(*gh.item);
                    }

                    T x = stva.value;

                    // implement our own __ballot
                    unsigned status_bf = (stva.status == 'p') ? (0x1u << lane) : 0;
                    for (int i = 1; i < Gpu::Device::warp_size; i *= 2) {
                        status_bf |= sg.shuffle_xor(status_bf, i);
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
                            x += sg.shuffle_down(x,i);
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
        Gpu::dtoh_memcpy_async(&totalsum, totalsum_p, sizeof(T));
    }
    Gpu::streamSynchronize();
    The_Arena()->free(dp);

    AMREX_GPU_ERROR_CHECK();

    return totalsum;
}

#elif defined(AMREX_USE_CUDA) && defined(__CUDACC__) && (__CUDACC_VER_MAJOR__ >= 11)

template <typename T, typename N, typename FIN, typename FOUT, typename TYPE,
          typename M=std::enable_if_t<std::is_integral<N>::value &&
                                      (std::is_same<std::decay_t<TYPE>,Type::Inclusive>::value ||
                                       std::is_same<std::decay_t<TYPE>,Type::Exclusive>::value)> >
T PrefixSum (N n, FIN && fin, FOUT && fout, TYPE, RetSum a_ret_sum = retSum)
{
    if (n <= 0) return 0;
    constexpr int nwarps_per_block = AMREX_HIP_OR_CUDA_OR_DPCPP(4,4,8);
    constexpr int nthreads = nwarps_per_block*Gpu::Device::warp_size; // # of threads per block
    constexpr int nelms_per_thread = 15;
    constexpr int nelms_per_block = nthreads * nelms_per_thread;
    int nblocks = (n + nelms_per_block - 1) / nelms_per_block;
    std::size_t sm = 0;
    auto stream = Gpu::gpuStream();

    typedef cub::ScanTileState<T> ScanTileState;
    std::size_t tile_state_size = 0;
    ScanTileState::AllocationSize(nblocks, tile_state_size);

    std::size_t nbytes_tile_state = Arena::align(tile_state_size);
    std::size_t nbytes_block_id = Arena::align(sizeof(unsigned int));
    auto dp = (char*)(The_Arena()->alloc(nbytes_tile_state+nbytes_block_id));
    auto tile_state_p = dp;

    ScanTileState tile_state;
    tile_state.Init(nblocks, tile_state_p, tile_state_size); // Init ScanTileState on host

    // Init ScanTileState on device
    amrex::launch((nblocks+nthreads-1)/nthreads, nthreads, 0, stream, [=] AMREX_GPU_DEVICE ()
    {
        const_cast<ScanTileState&>(tile_state).InitializeStatus(nblocks);
    });

    struct InputIterator {
        InputIterator (FIN const& a_fin, N a_offset)
            : m_fin(a_fin), m_offset(a_offset) {}

        AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        T operator[] (N i) const { return m_fin(i+m_offset); }

    private:
        FIN m_fin;
        N m_offset;
    };

    T* totalsum_p = (a_ret_sum) ? (T*)(The_Pinned_Arena()->alloc(sizeof(T))) : nullptr;

    amrex::launch(nblocks, nthreads, sm, stream,
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

        T data[nelms_per_thread];
        if (static_cast<int>(iend-ibegin) == nelms_per_block) {
            BlockLoad(temp_storage.load).Load(InputIterator(fin, ibegin), data);
        } else {
            BlockLoad(temp_storage.load).Load(InputIterator(fin, ibegin), data,
                                              iend-ibegin, 0); // padding with 0
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
            if (offset >= iend) break;
            fout(offset, data[i]);
        }
    });

    Gpu::streamSynchronize();
    AMREX_GPU_ERROR_CHECK();

    The_Arena()->free(dp);

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
                unsigned int bid = Gpu::Atomic::Inc(virtual_block_id_p, gridDim.x);
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
#if defined(AMREX_USE_CUDA) && defined(__CUDACC__) && (__CUDACC_VER_MAJOR__ >= 11)
    void* d_temp = nullptr;
    std::size_t temp_bytes = 0;
    AMREX_GPU_SAFE_CALL(cub::DeviceScan::ExclusiveSum(d_temp, temp_bytes, in, out, n,
                                                      Gpu::gpuStream()));
    d_temp = The_Arena()->alloc(temp_bytes);
    AMREX_GPU_SAFE_CALL(cub::DeviceScan::ExclusiveSum(d_temp, temp_bytes, in, out, n,
                                                      Gpu::gpuStream()));
    T in_last = 0, out_last = 0;
    if (a_ret_sum) {
        Gpu::dtoh_memcpy_async(&in_last, in+(n-1), sizeof(T));
        Gpu::dtoh_memcpy_async(&out_last, out+(n-1), sizeof(T));
    }
    Gpu::streamSynchronize();
    The_Arena()->free(d_temp);
    AMREX_GPU_ERROR_CHECK();
    return in_last+out_last;
#elif defined(AMREX_USE_HIP)
    void* d_temp = nullptr;
    std::size_t temp_bytes = 0;
    AMREX_GPU_SAFE_CALL(rocprim::exclusive_scan(d_temp, temp_bytes, in, out, T{0}, n,
                                                rocprim::plus<T>(), Gpu::gpuStream()));
    d_temp = The_Arena()->alloc(temp_bytes);
    AMREX_GPU_SAFE_CALL(rocprim::exclusive_scan(d_temp, temp_bytes, in, out, T{0}, n,
                                                rocprim::plus<T>(), Gpu::gpuStream()));
    T in_last = 0, out_last = 0;
    if (a_ret_sum) {
        Gpu::dtoh_memcpy_async(&in_last, in+(n-1), sizeof(T));
        Gpu::dtoh_memcpy_async(&out_last, out+(n-1), sizeof(T));
    }
    Gpu::streamSynchronize();
    The_Arena()->free(d_temp);
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
T PrefixSum (N n, FIN && fin, FOUT && fout, TYPE)
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
    std::partial_sum(in, in+n, out);
    return (n > 0) ? out[n-1] : T(0);
}

// The return value is the total sum.
template <typename N, typename T, typename M=std::enable_if_t<std::is_integral<N>::value> >
T ExclusiveSum (N n, T const* in, T * out, RetSum /*a_ret_sum*/ = retSum)
{
    if (n > 0) {
        out[0] = 0;
        std::partial_sum(in, in+n-1, out+1);
        return in[n-1]+out[n-1];
    } else {
        return 0;
    }
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
#else
        return std::partial_sum(begin, end, result);
#endif
    }

    template<class InIter, class OutIter>
    OutIter exclusive_scan(InIter begin, InIter end, OutIter result)
    {
#if defined(AMREX_USE_GPU)
        auto N = std::distance(begin, end);
        Scan::ExclusiveSum(N, &(*begin), &(*result), Scan::noRetSum);
        OutIter result_end = result;
        std::advance(result_end, N);
        return result_end;
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
