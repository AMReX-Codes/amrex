#ifndef AMREX_SUNMEMORY_H_
#define AMREX_SUNMEMORY_H_

#include <sundials/sundials_context.h>
#include <sundials/sundials_memory.h>

namespace amrex {
namespace sundials {

/**
 * \brief A SUNDIALS SUNMemoryHelper interface to the amrex::Arena.
 *
 * This class allows SUNDIALS to allocate memory using the amrex::Arena.
 */

class MemoryHelper {
public:
    MemoryHelper(::sundials::Context* sunctx);
    ~MemoryHelper();
    MemoryHelper(const MemoryHelper& helper);
    MemoryHelper(MemoryHelper&& helper);
    MemoryHelper& operator=(MemoryHelper other);
    MemoryHelper& operator=(MemoryHelper&& other);
    operator SUNMemoryHelper() const { return helper; }

    static void Initialize(int nthreads);
    static void Finalize();
private:
    ::sundials::Context* sunctx;
    SUNMemoryHelper helper;
};


/**
 * \brief Get the i-th SUNMemoryHelper instance.
 *
 * \param[in] i The index of the SUNMemoryHelper to return.
 *
 * SUNMemoryHelper_i is paired with the i-th SUNContext.
 */

MemoryHelper* The_SUNMemory_Helper(int i = amrex::OpenMP::get_thread_num());

/**
 * \brief Get the proper amrex::Arena based on the SUNMemoryType.
 *
 * \param[in] type A SUNDIALS SUNMemoryType.
 */

amrex::Arena* getArena (SUNMemoryType type);

}//sundials
}//amrex

#endif
