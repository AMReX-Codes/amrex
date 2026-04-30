#ifndef AMREX_SUNMEMORY_H_
#define AMREX_SUNMEMORY_H_

#include <AMReX_Arena.H>
#include <AMReX_Sundials_Core.H>
#include <sundials/sundials_context.h>
#include <sundials/sundials_memory.h>

/**
 * \file AMReX_SUNMemory.H
 *
 * SUNMemory helper implementations that let SUNDIALS use AMReX arenas.
 */

namespace amrex::sundials {

/**
 * \brief A SUNDIALS SUNMemoryHelper interface to the amrex::Arena.
 *
 * This class allows SUNDIALS to allocate memory using the amrex::Arena.
 */
class MemoryHelper {
public:
    //! Construct a helper bound to the supplied SUNDIALS context \p sunctx.
    MemoryHelper(::sundials::Context* sunctx);
    //! Destroy the associated SUNMemoryHelper and free resources.
    ~MemoryHelper();
    //! Copy-construct a helper (shares the SUNMemoryHelper handle).
    MemoryHelper(const MemoryHelper& rhs);
    //! Move-construct, transferring ownership of the helper handle.
    MemoryHelper(MemoryHelper&& rhs) noexcept;
    MemoryHelper& operator=(MemoryHelper const& rhs) = delete;
    MemoryHelper& operator=(MemoryHelper&& rhs) noexcept = delete;
    //! Implicit conversion to the underlying SUNMemoryHelper handle.
    operator SUNMemoryHelper() const { return helper; }

    /**
     * \brief Allocate one helper per OpenMP thread using the given count.
     *
     * \param nthreads Number of thread-local SUNMemoryHelper instances created.
     */
    static void Initialize(int nthreads);
    //! Destroy all helpers created by Initialize().
    static void Finalize();
private:
    SUNMemoryHelper helper;
    ::sundials::Context* sunctx;
};


/**
 * \brief Get the OpenMP-thread-local SUNMemoryHelper instance.
 *
 * \param[in] i Thread index (defaults to current thread).
 * \return Helper aligned with the i-th SUNContext.
 */
MemoryHelper* The_SUNMemory_Helper(int i = amrex::OpenMP::get_thread_num());

/**
 * \brief Return the AMReX arena used to satisfy allocations of \p type.
 *
 * \param[in] type SUNDIALS memory type requested (host, pinned, device, etc.).
 * \return Corresponding amrex::Arena pointer.
 */
amrex::Arena* getArena (SUNMemoryType type);

}//amrex::sundials

#endif
