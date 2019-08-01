#ifndef AMREX_SCAN_H_
#define AMREX_SCAN_H_

#include <AMReX_Gpu.H>
#include <AMReX_Arena.H>
#include <cstdint>
#include <type_traits>

namespace amrex {
namespace Scan {

namespace detail {

template <typename T>
struct alignas(8) STVA
{
    char status;
    T value;
};

template <typename T, bool SINGLE_WORD> struct BlockStatus {};

template <typename T>
struct alignas(8) BlockStatus<T, true>
{
    STVA<T> stva;

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void write (char a_status, T a_value) {
        static_assert(sizeof(BlockStatus<T,true>) == 8,
                      "BlockStatus: wrong size");
        int64_t tmp;
        auto ptmp = reinterpret_cast<BlockStatus<T,true>*>(&tmp);
        ptmp->stva.status = a_status;
        ptmp->stva.value = a_value;
        auto pthis = reinterpret_cast<int64_t*>(this);
        *pthis = tmp;
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    T get_aggregate() const { return stva.value; }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    STVA<T> read () volatile {
        int64_t tmp = *(reinterpret_cast<int64_t volatile*>(this));
        auto ptmp = reinterpret_cast<BlockStatus<T,true>*>(&tmp);
        return { ptmp->stva.status, ptmp->stva.value };
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void set_status (char a_status) { stva.status = a_status; }
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
        __threadfence();
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
};

enum class Type { inclusive, exclusive };

template <typename T>
void Sum (int n, T const* in, T * out, Type type)
{
    if (n <= 0) return;
    constexpr int nthreads = 256;
    constexpr int nchunks = 12;
    constexpr int nelms_per_block = nthreads * nchunks;
    int nblocks = (n + nelms_per_block - 1) / nelms_per_block;
    std::size_t sm = sizeof(T) * Gpu::Device::warp_size + sizeof(int);
    auto stream = Gpu::gpuStream();

    // temporary memory
    typedef typename std::conditional<sizeof(BlockStatus<T,true>) <= 8,
        BlockStatus<T,true>, BlockStatus<T,false> >::type BlockStatusT;
    amrex::AsyncArray<BlockStatusT > block_status_aa(nblocks);
    BlockStatusT* block_status_p = block_status_aa.data();

    Gpu::DeviceScalar<unsigned int> virtual_block_id_ds;
    unsigned int* virtual_block_id_p = virtual_block_id_ds.dataPtr();

    amrex::ParallelFor(nblocks, [=] AMREX_GPU_DEVICE (int i) noexcept {
        BlockStatusT& block_status = block_status_p[i];
        block_status.set_status('x');
        if (i == 0) *virtual_block_id_p = 0;
    });

    amrex::launch_global<<<nblocks, nthreads, sm, stream>>>(
    [=] AMREX_GPU_DEVICE () noexcept
    {
        int lane = threadIdx.x % Gpu::Device::warp_size;
        int warp = threadIdx.x / Gpu::Device::warp_size;
        int nwarps = blockDim.x / Gpu::Device::warp_size;

        amrex::Gpu::SharedMemory<T> gsm;
        T* shared = gsm.dataPtr();
        int& virtual_block_id = *((int*)(shared+Gpu::Device::warp_size));

        // First of all, get block virtual id.  We must do this to
        // avoid deadlock because CUDA may launch blocks in any order.
        // Anywhere in this function, we should not use blockIdx.
        if (gridDim.x > 1) {
            if (threadIdx.x == 0) {
                unsigned int bid = Gpu::Atomic::Inc(virtual_block_id_p, gridDim.x);
                virtual_block_id = bid;
            }
            __syncthreads();
        } else {
            virtual_block_id = 0;
        }

        // Each block processes [ibegin,iend).
        int ibegin = nelms_per_block * virtual_block_id;
        int iend = amrex::min(ibegin+nelms_per_block, n);
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
            int offset = ibegin + ichunk*blockDim.x;
            if (offset >= iend) break;

            offset += threadIdx.x;
            T x0 = (offset < iend) ? in[offset] : 0; 
            T x = x0;
            // Scan within a warp
            for (int i = 1; i <= Gpu::Device::warp_size; i *= 2) {
                T s = __shfl_up_sync(0xffffffff, x, i);
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
            if (warp == 0 && lane < nwarps) {
                T y = shared[lane];
                int mask = (1 << nwarps) - 1;
                for (int i = 1; i <= nwarps; i *= 2) {
                    T s = __shfl_up_sync(mask, y, i, nwarps);
                    if (lane >= i) y += s;
                }

                shared[lane] = y;
            }

            __syncthreads();

            // shared[0:nwarps) holds the inclusive sum of warp sums.
            
            // Also note x still holds the inclusive sum within the
            // warp.  Given these two, we can compute the inclusive
            // sum within this chunk.
            T sum_prev_warp = (warp == 0) ? 0 : shared[warp-1];
            tmp_out[ichunk] = sum_prev_warp + sum_prev_chunk +
                ((type == Type::inclusive) ? x : x-x0);
            sum_prev_chunk += shared[nwarps-1];
        }

        // sum_prev_chunk now holds the sum of the whole block.
        if (threadIdx.x == 0 && gridDim.x > 1) {
            block_status.write((virtual_block_id == 0) ? 'p' : 'a',
                               sum_prev_chunk);
        }

        if (virtual_block_id == 0) {
            for (int ichunk = 0; ichunk < nchunks; ++ichunk) {
                int offset = ibegin + ichunk*blockDim.x + threadIdx.x;
                if (offset >= iend) break;
                out[offset] = tmp_out[ichunk];
            }
        } else if (virtual_block_id > 0) {

            if (warp == 0) {
                T exclusive_prefix = 0;
                BlockStatusT volatile* pbs = block_status_p;
                for (int iblock0 = virtual_block_id-1; iblock0 >= 0; iblock0 -= Gpu::Device::warp_size)
                {
                    int iblock = iblock0-lane;
                    STVA<T> stva{'p', 0};
                    do {
                        if (iblock >= 0) stva = pbs[iblock].read();
                    } while (__any_sync(0xffffffff,stva.status=='x'));

                    T x = stva.value;

                    unsigned status_bf = __ballot_sync(0xffffffff, stva.status == 'p');
                    unsigned int bit_mask = 0x1;
                    bool stop_lookback = status_bf & bit_mask;
                    if (stop_lookback == false) {
                        for (int i = 1; i < Gpu::Device::warp_size; ++i) {
                            T s = __shfl_down_sync(0xffffffff, x, i);
                            if (lane == 0) { x += s; }
                            bit_mask <<= 1;
                            if (status_bf & bit_mask) {
                                stop_lookback = true;
                                break;
                            }
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
                int offset = ibegin + ichunk*blockDim.x + threadIdx.x;
                if (offset < iend) out[offset] = tmp_out[ichunk] + exclusive_prefix;
            }
        }
    });

    AMREX_GPU_ERROR_CHECK();
}

}

template <typename N, typename T, typename M=amrex::EnableIf_t<std::is_integral<N>::value> >
void InclusiveSum (N n, T const* in, T * out)
{
    AMREX_ALWAYS_ASSERT(static_cast<long>(n) < static_cast<long>(std::numeric_limits<int>::max()));
    detail::Sum(n, in, out, detail::Type::inclusive);
}

template <typename N, typename T, typename M=amrex::EnableIf_t<std::is_integral<N>::value> >
void ExclusiveSum (N n, T const* in, T * out)
{
    AMREX_ALWAYS_ASSERT(static_cast<long>(n) < static_cast<long>(std::numeric_limits<int>::max()));
    detail::Sum(n, in, out, detail::Type::exclusive);
}

}}

#endif
