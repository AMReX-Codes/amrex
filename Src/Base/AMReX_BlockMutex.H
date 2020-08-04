#ifndef AMREX_BLOCK_MUTEX_H_
#define AMREX_BLOCK_MUTEX_H_

#include <AMReX_Gpu.H>

namespace amrex {

#ifdef AMREX_USE_GPU
    
struct BlockMutex
{
    union state_t
    {
        struct II { int blockid; int count; } data;
        unsigned long long ull;
    };

    AMREX_GPU_HOST_DEVICE
    static constexpr state_t FreeState () noexcept {
        return state_t{{-1,0}};
    }

    static void init_states (state_t* state, int N) noexcept;

    explicit BlockMutex (int N) noexcept;
    
    ~BlockMutex ();

    void operator= (BlockMutex const&) = delete;

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void lock (int i) noexcept {
#ifdef AMREX_USE_DPCPP
// xxxxx DPCPP todo
        amrex::ignore_unused(i);
#else
        int blockid = blockIdx.z*blockDim.x*blockDim.y + blockIdx.y*blockDim.x + blockIdx.x;
        state_t old = m_state[i];
        state_t assumed;
        do {
            assumed = old;
            state_t val;
            val.data.blockid = blockid;
            if (assumed.data.blockid == blockid) {
                // Already locked by another thread in this block. Need to ++count.
                val.data.count = assumed.data.count + 1;
            } else {
                // Currently unlocked or locked by another block.  Need to lock.
                val.data.count = 1;
                assumed = FreeState();
            }
            old.ull = atomicCAS((unsigned long long*)(m_state+i), assumed.ull, val.ull);
        } while (assumed.ull != old.ull);
#endif
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void unlock (int i) noexcept {
#ifdef AMREX_USE_DPCPP
// xxxxx DPCPP todo
        amrex::ignore_unused(i);
#else
        state_t old = m_state[i];
        state_t assumed;
        do {
            assumed = old;
            state_t val;
            if (assumed.data.count == 1) {
                // Need to unlock
                val = FreeState();
            } else {
                // --count, but do NOT unlock
                val = assumed;
                --val.data.count;
            }
            old.ull = atomicCAS((unsigned long long*)(m_state+i), assumed.ull, val.ull);
        } while (assumed.ull != old.ull);
#endif
    }

private:

    int m_nstates;
    state_t* m_state;
};
#endif  
}
#endif
