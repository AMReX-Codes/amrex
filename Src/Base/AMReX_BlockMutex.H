#ifndef AMREX_BLOCK_MUTEX_H_
#define AMREX_BLOCK_MUTEX_H_

#include <AMReX_Gpu.H>

namespace amrex {

#ifdef AMREX_USE_GPU
    
struct BlockMutex
{
    union state_t
    {
        struct { int blockid; int count; };
        unsigned long long ull;
    };

    AMREX_GPU_HOST_DEVICE
    static constexpr state_t FreeState () noexcept {
        return state_t{{-1,0}};
    }

    static void init_states (state_t* state, int N) {
        launch_global<<<(N+255)/256,256>>>([=] AMREX_GPU_DEVICE () noexcept
        {
            int i = threadIdx.x + blockIdx.x*blockDim.x;
            if (i < N) state[i] = FreeState();
        });
    }

    explicit BlockMutex (int N)
        : m_nstates(N)
    {
        static_assert(sizeof(unsigned long long) == 2*sizeof(int) and
                      sizeof(unsigned long long) == sizeof(state_t),
                      "BlockMutex: wrong size");
        // The first 4 bytes of unsigned long stores blockIdx.
        // The second 4 bytes, count.
        // The initial values are -1 and 0.
        cudaMalloc(&m_state, sizeof(state_t)*m_nstates);
        // In amrex, we will use amrex::ParallelFor
        init_states(m_state, m_nstates);
    }

    ~BlockMutex () {
        cudaFree(m_state);
    }

    void operator= (BlockMutex const&) = delete;

    AMREX_GPU_DEVICE
    void lock (int i) {
        int blockid = blockIdx.z*blockDim.x*blockDim.y + blockIdx.y*blockDim.x + blockIdx.x;
        state_t old = m_state[i];
        state_t assumed;
        do {
            assumed = old;
            state_t val;
            val.blockid = blockid;
            if (assumed.blockid == blockid) {
                // Already locked by another thread in this block. Need to ++count.
                val.count = assumed.count + 1;
            } else {
                // Currently unlocked or locked by another block.  Need to lock.
                val.count = 1;
                assumed = FreeState();
            }
            old.ull = atomicCAS((unsigned long long*)(m_state+i), assumed.ull, val.ull);
        } while (assumed.ull != old.ull);
    }

    AMREX_GPU_DEVICE
    void unlock (int i) {
        state_t old = m_state[i];
        state_t assumed;
        do {
            assumed = old;
            state_t val;
            if (assumed.count == 1) {
                // Need to unlock
                val = FreeState();
            } else {
                // --count, but do NOT unlock
                val = assumed;
                --val.count;
            }
            old.ull = atomicCAS((unsigned long long*)(m_state+i), assumed.ull, val.ull);
        } while (assumed.ull != old.ull);
    }

private:

    int m_nstates;
    state_t* m_state;
};
#endif  
}
#endif
