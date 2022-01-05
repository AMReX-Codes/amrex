#ifndef AMREX_SUNMEMORY_H_
#define AMREX_SUNMEMORY_H_

#include <sundials/sundials_context.h>
#include <sundials/sundials_memory.h>

namespace amrex {
namespace sundials {

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

#ifdef AMREX_USE_OMP
MemoryHelper* The_SUNMemory_Helper(int i = omp_get_thread_num());
#else
MemoryHelper* The_SUNMemory_Helper(int i = 0);
#endif
amrex::Arena* getArena (SUNMemoryType type);

}//sundials
}//amrex

#endif
