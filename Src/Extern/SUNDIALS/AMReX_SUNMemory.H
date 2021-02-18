#include <AMReX.H>
#include <sundials/sundials_memory.h>

namespace amrex {
namespace sundials {
  
void Initialize();
void Finalize();

class MemoryHelper {
public:
  MemoryHelper();
  ~MemoryHelper();
  operator SUNMemoryHelper() const { return helper; }
private:
  SUNMemoryHelper helper;
};

MemoryHelper* The_SUNMemory_Helper ();

}
}
