#include <AMReX_BlockMutex.H>

namespace amrex {

#ifdef AMREX_USE_GPU

void BlockMutex::init_states (state_t* state, int N) noexcept {
#ifdef AMREX_USE_DPCPP
    amrex::ignore_unused(state,N);
    amrex::Abort("xxxxx DPCPP todo");
#else
    amrex::launch((N+255)/256, 256, Gpu::nullStream(),
    [=] AMREX_GPU_DEVICE () noexcept
    {
        int i = threadIdx.x + blockIdx.x*blockDim.x;
        if (i < N) state[i] = FreeState();
    });
#endif
}

BlockMutex::BlockMutex (int N) noexcept
    : m_nstates(N)
{
    static_assert(sizeof(unsigned long long) == 2*sizeof(int) and
                  sizeof(unsigned long long) == sizeof(state_t),
                  "BlockMutex: wrong size");
    // The first 4 bytes of unsigned long stores blockIdx.
    // The second 4 bytes, count.
    // The initial values are -1 and 0.
    m_state = static_cast<state_t*>(The_Device_Arena()->alloc(sizeof(state_t)*m_nstates));
    init_states(m_state, m_nstates);
}

BlockMutex::~BlockMutex () {
    The_Device_Arena()->free(m_state);
}

#endif
    
}

