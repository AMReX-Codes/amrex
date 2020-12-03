#ifndef KERNELTIMER_H_
#define KERNELTIMER_H_

#if defined(AMREX_USE_KERNEL_TIMER) && defined(AMREX_USE_GPU)
#include <AMReX_REAL.H>

namespace amrex {

/**
 * \brief Defines a timer object to be used on GPU; measures summed thread cycles.
 */
class KernelTimer
{
public:
    /** Constructor.
     * \param[in] do_timing Controls whether timer is active.
     * \param[in,out] cost Pointer to cost which holds summed thread cycles
     * (for performance, it is recommended to allocate pinned host memory).
     */
    AMREX_GPU_DEVICE
    KernelTimer (const bool do_timing, amrex::Real* cost)
        : m_do_timing(do_timing), m_cost(cost) {
	if (do_timing && cost) {
#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
            // Start the timer
            m_wt = clock64();

#elif defined(AMREX_USE_DPCPP)
            // To be updated
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE( do_timing == false,
                                              "KernelTimer not yet supported for this hardware." );
#endif
	}
    }

    //! Destructor.
    AMREX_GPU_DEVICE
    ~KernelTimer () {
	if (m_do_timing && m_cost) {
#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
            m_wt = clock64() - m_wt;
            amrex::Gpu::Atomic::Add( m_cost, amrex::Real(m_wt));
#elif defined(AMREX_USE_DPCPP)
            // To be updated
#endif
	}
    }

private:
    //! Stores whether kernel timer is active.
    bool m_do_timing;

    //! Location in which to accumulate costs from all threads.
    amrex::Real* m_cost;

    //! Store the time difference (cost) from a single thread.
    long long int m_wt;
};


// Experimental; not guaranteed thread safe, only for testing
class KernelTimerShared
{
public:
    /** Constructor.
     * \param[in] do_timing Controls whether timer is active.
     * \param[in,out] cost Pointer to cost which holds summed thread cycles.
     */
    AMREX_GPU_DEVICE
    KernelTimerShared (const bool do_timing, amrex::Real* cost)
	: m_do_timing(do_timing), m_cost(cost) {
	if (do_timing && cost) {
#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
	    // Store threads' accumulated cycles for present block
	    __shared__ amrex::Real cost_shared;
	    cost_shared = 0.;
	    m_cost_shared = &cost_shared;

	    // Start the timer
	    m_wt = clock64();

#elif defined(AMREX_USE_DPCPP)
	    // To be updated
	    AMREX_ALWAYS_ASSERT_WITH_MESSAGE( do_timing == false,
					      "KernelTimer not yet supported for this hardware." );
#endif
	}
    }

    //! Destructor.
    AMREX_GPU_DEVICE
    ~KernelTimerShared () {
	if (m_do_timing && m_cost) {
#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
	m_wt = clock64() - m_wt;
	amrex::Gpu::Atomic::Add( m_cost_shared, amrex::Real(m_wt));

	// Make sure everyone finished
	__syncthreads();

	// Copy accumulated cycles from each block's smem back to input cost location
	if (threadIdx.x == 0) {
	    amrex::Gpu::Atomic::Add( m_cost, *m_cost_shared);
	}
#elif defined(AMREX_USE_DPCPP)
	// To be updated
#endif
	}
    }

private:
    //! Stores whether kernel timer is active.
    bool m_do_timing;

    //! Location in which to accumulate costs from all threads.
    amrex::Real* m_cost;

    //! Store the time difference (cost) from a single thread.
    long long int m_wt;

#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
    //! Location in which which to accumulate from each thread.
    amrex::Real* m_cost_shared;
#endif
};

class KernelTimerSharedPTX
{
public:
    /** Constructor.
     * \param[in] do_timing Controls whether timer is active.
     * \param[in,out] cost Pointer to cost which holds summed thread cycles.
     */
    AMREX_GPU_DEVICE
    KernelTimerSharedPTX (const bool do_timing, amrex::Real* cost)
	: m_do_timing(do_timing), m_cost(cost) {
	if (do_timing && cost) {
#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
	    // Store threads' accumulated cycles for present block
	    __shared__ amrex::Real cost_shared;
	    cost_shared = 0.;
	    m_cost_shared = &cost_shared;

	    // Start the timer
        asm("mov.u32 %0, %%clock;" : "=r"(m_wt));

#elif defined(AMREX_USE_DPCPP)
	    // To be updated
	    AMREX_ALWAYS_ASSERT_WITH_MESSAGE( do_timing == false,
					      "KernelTimer not yet supported for this hardware." );
#endif
	}
    }

    //! Destructor.
    AMREX_GPU_DEVICE
    ~KernelTimerSharedPTX () {
	if (m_do_timing && m_cost) {
#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
    asm(".reg .u32 t1;\n\t"         // temp reg t1
        " mov.u32 t1, %%clock;\n\t" // t1 = clock
        " sub.u32 %0, t1, %1;"      // cost = t1 - cost
        : "=r"(m_wt) : "r" (m_wt));
	amrex::Gpu::Atomic::Add( m_cost_shared, amrex::Real(m_wt));

	// Make sure everyone finished
	__syncthreads();

	// Copy accumulated cycles from each block's smem back to input cost location
	if (threadIdx.x == 0) {
	    amrex::Gpu::Atomic::Add( m_cost, *m_cost_shared);
	}
#elif defined(AMREX_USE_DPCPP)
	// To be updated
#endif
	}
    }

private:
    //! Stores whether kernel timer is active.
    bool m_do_timing;

    //! Location in which to accumulate costs from all threads.
    amrex::Real* m_cost;

    //! Store the time difference (cost) from a single thread.
    unsigned m_wt;

#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
    //! Location in which which to accumulate from each thread.
    amrex::Real* m_cost_shared;
#endif
};
    
}

#endif /* AMREX_USE_KERNEL_TIMER && AMREX_USE_GPU*/
#endif /* KERNELTIMER_H_ */
