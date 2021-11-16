#ifndef AMREX_SUNMEMORY_H_
#define AMREX_SUNMEMORY_H_

#include <AMReX.H>
#include <AMReX_Arena.H>
#include <sundials/sundials_context.h>
#include <sundials/sundials_memory.h>

namespace amrex {
namespace sundials {

class MemoryHelper {
public:
  MemoryHelper(::sundials::Context*);
  ~MemoryHelper();
  operator SUNMemoryHelper() const { return helper; }
  static void Initialize(::sundials::Context*);
  static void Finalize();
private:
  SUNMemoryHelper helper;
};

MemoryHelper* The_SUNMemory_Helper();
amrex::Arena* getArena (SUNMemoryType type);

}
}

#endif
