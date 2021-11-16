#include <AMReX_Sundials.H>

namespace amrex {
namespace sundials {

namespace {
  bool initialized = false;
  ::sundials::Context* the_sundials_context = nullptr;
}

void Initialize()
{
  if (initialized) return;
  initialized = true;

  BL_ASSERT(the_sundials_context == nullptr);
  the_sundials_context = new ::sundials::Context();
  MemoryHelper::Initialize(the_sundials_context);
}

void Finalize()
{
  BL_ASSERT(the_sundials_context != nullptr);
  MemoryHelper::Finalize();
  delete the_sundials_context;
  the_sundials_context = nullptr;
}

::sundials::Context* The_Sundials_Context()
{
  BL_ASSERT(the_sundials_context != nullptr);
  return the_sundials_context;
}

}//sundials
}//amrex
