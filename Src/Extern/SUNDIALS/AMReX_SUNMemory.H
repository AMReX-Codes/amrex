#ifndef AMREX_SUNMEMORY_H_
#define AMREX_SUNMEMORY_H_

#include <AMReX.H>
#include <AMReX_Arena.H>
#include <sundials/sundials_memory.h>

namespace amrex {
namespace sundials {

class MemoryHelper {
public:
  MemoryHelper();
  ~MemoryHelper();
  operator SUNMemoryHelper() const { return helper; }
private:
  SUNMemoryHelper helper;
};

void Initialize();
void Finalize();

MemoryHelper* The_SUNMemory_Helper();
amrex::Arena* getArena (SUNMemoryType type);

}
}

#endif
