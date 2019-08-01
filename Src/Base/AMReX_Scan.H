#ifndef AMREX_SCAN_H_
#define AMREX_SCAN_H_

#include <AMReX_Gpu.H>
#include <AMReX_Arena.H>

namespace amrex {
namespace Scan {

namespace detail {

template <typename T>
struct BlockStatus {
    T aggregate;
    T inclusive;
    char status;
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
    amrex::AsyncArray<BlockStatus<T> > block_status_aa(nblocks);
    BlockStatus<T>* block_status_p = block_status_aa.data();

    Gpu::DeviceScalar<unsigned int> virtual_block_id_ds;
    unsigned int* virtual_block_id_p = virtual_block_id_ds.dataPtr();

    amrex::ParallelFor(nblocks, [=] AMREX_GPU_DEVICE (int i) noexcept {
        BlockStatus<T>& block_status = block_status_p[i];
        block_status.status = 'x';
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
        BlockStatus<T>& block_status = block_status_p[virtual_block_id];

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
            block_status.aggregate = sum_prev_chunk;
            if (virtual_block_id == 0) block_status.inclusive = sum_prev_chunk;
            // TODO: optimize away the memory fence.
            __threadfence(); // This ensures the status flag is touched after the fence.
            block_status.status = (virtual_block_id == 0) ? 'p' : 'a';
        }

        if (virtual_block_id == 0) {
            for (int ichunk = 0; ichunk < nchunks; ++ichunk) {
                int offset = ibegin + ichunk*blockDim.x + threadIdx.x;
                if (offset >= iend) break;
                out[offset] = tmp_out[ichunk];
            }
        } else if (virtual_block_id > 0) {
            T& exclusive_prefix = shared[0];

            if (threadIdx.x == 0) {
                exclusive_prefix = 0;
                volatile BlockStatus<T>* pbs = block_status_p;
                for (int iblock = virtual_block_id-1; iblock >= 0; --iblock) {
                    char prev_status;
                    while ((prev_status = pbs[iblock].status) == 'x') {};
                    if (prev_status == 'a') {
                        exclusive_prefix += pbs[iblock].aggregate;
                    } else {
                        exclusive_prefix += pbs[iblock].inclusive;
                        break;
                    }
                }

                block_status.inclusive = block_status.aggregate + exclusive_prefix;
                __threadfence();
                block_status.status = 'p';
            }

            __syncthreads();

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
