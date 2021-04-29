#include <AMReX.H>
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

}
}
