/* Copyright 2019-2020 Michael Rowan, Axel Huebl, Kevin Gott
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef KERNELTIMER_H_
#define KERNELTIMER_H_

#ifdef AMREX_USE_GPU
/**
 * \brief Defines a timer object to be used on GPU; measured summed thread time.
 */
class
KernelTimer
{
public:
  /** Constructor.
   * \param[in] do_timing Controls whether timer is active.
   * \param[in,out] cost Pointer to cost of corresponding box.
   */
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    KernelTimer(const bool do_timing, amrex::Real* cost) :
      m_do_timing(do_timing), m_cost(cost) {
      if (do_timing && cost) {
#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
	// To store present thread's time for deposition
	__shared__ amrex::Real cost_shared;
	cost_shared = 0.;
	m_cost_shared = &cost_shared;

	// Start the timer
	m_wt = clock64();

#elif defined(AMREX_USE_DPCPP)
	// Todo: DPC++
	AMREX_ALWAYS_ASSERT_WITH_MESSAGE( do_timing == false,
					  "KernelTimer not supported for SYCL." );
#endif
      }
    }

    //** Destructor. */
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    ~KernelTimer() {
      if (m_do_timing && m_cost) {
#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
	m_wt = clock64() - m_wt;
	amrex::Gpu::Atomic::Add( m_cost_shared, amrex::Real(m_wt));

	// Make sure everyone finished current deposit
	__syncthreads();

	// Copy tmp smem arrays back to global
	if (threadIdx.x == 0) {
	  amrex::Gpu::Atomic::Add( m_cost, *m_cost_shared);
	}
#endif
      }
    }

private:
  /** Stores whether kernel timer is active. */
  bool m_do_timing;

  /** Location in which to accumulate costs from all blocks.
   *  (after threads accumulate costs in shared memory).
   */
  amrex::Real* m_cost;

  /** Store the time difference (cost) from a single thread. */
  long long int m_wt;

#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
  /** Location in which which to accumulate from each thread
   *  (this will be declared shared in constructor).
   */
  amrex::Real* m_cost_shared;
#endif
};

#endif
#endif
