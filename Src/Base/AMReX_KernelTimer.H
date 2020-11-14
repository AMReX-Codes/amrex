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
     * \param[in,out] cost Pointer to cost of corresponding box.
     */
    AMREX_GPU_DEVICE
    KernelTimer (const bool do_timing, amrex::Real* cost)
	: m_do_timing(do_timing), m_cost(cost) {
	if (do_timing && cost) {
#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
	    // To store present thread's time for deposition
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
    ~KernelTimer () {
	if (m_do_timing && m_cost) {
#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
	m_wt = clock64() - m_wt;
	amrex::Gpu::Atomic::Add( m_cost_shared, amrex::Real(m_wt));

	// Make sure everyone finished current deposit
	__syncthreads();

	// Copy temp smem arrays back to global
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
    
}

#endif /* AMREX_USE_KERNEL_TIMER && AMREX_USE_GPU*/
#endif /* KERNELTIMER_H_ */
