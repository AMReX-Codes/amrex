#ifndef AMREX_SUNDIALS_CORE_H_
#define AMREX_SUNDIALS_CORE_H_

#include <AMReX.H>
#include <AMReX_OpenMP.H>
#include <sundials/sundials_context.h>

namespace amrex::sundials {

/**
 * \brief Initialize the AMReX-SUNDIALS interface.
 *
 * \param[in] nthreads The number of threads that will be used with SUNDIALS.
 *
 * This will create the nthreads SUNDIALS context objects that are needed by
 * the SUNDIALS solver and vector objects. Called by amrex::Initialize.
 */

void Initialize(int nthreads);

/**
 * \brief Cleanup everything allocated for the AMReX-SUNDIALS interface.
 *
 * Called by amrex::Finalize.
 */

void Finalize();

/**
 * \brief Get the i-th SUNDIALS context instance.
 *
 * \param[in] i The index of the SUNDIALS context to return.
 *
 * A SUNDIALS context should not be used concurrently from different threads.
 */

::sundials::Context* The_Sundials_Context(int i = amrex::OpenMP::get_thread_num());

}

#endif
