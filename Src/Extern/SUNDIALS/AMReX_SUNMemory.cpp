#include <AMReX.H>
#include <AMReX_Gpu.H>
#include <AMReX_SUNMemory.H>
#if defined(AMREX_USE_HIP)
#include <sunmemory/sunmemory_hip.h>
#elif defined(AMREX_USE_CUDA)
#include <sunmemory/sunmemory_cuda.h>
#endif

namespace amrex {
namespace sundials {

namespace {
  int Alloc(SUNMemoryHelper, SUNMemory* memptr, size_t memsize, SUNMemoryType mem_type)
  {
    SUNMemory mem = SUNMemoryNewEmpty();

    mem->ptr = NULL;
    mem->own = SUNTRUE;
    mem->type = mem_type;

    if (mem_type == SUNMEMTYPE_HOST) {
      mem->ptr = The_Cpu_Arena()->alloc(memsize);
    } else if (mem_type == SUNMEMTYPE_UVM) {
      mem->ptr = The_Managed_Arena()->alloc(memsize);
    } else if (mem_type == SUNMEMTYPE_DEVICE) {
      mem->ptr = The_Device_Arena()->alloc(memsize);
    } else if (mem_type == SUNMEMTYPE_PINNED) {
      mem->ptr = The_Pinned_Arena()->alloc(memsize);
    } else {
      free(mem);
      return(-1);
    }

    *memptr = mem;
    return(0);
  }

  int Dealloc(SUNMemoryHelper, SUNMemory mem)
  {
    if (mem == NULL) return(0);

    if (mem->type == SUNMEMTYPE_HOST) {
      if (mem->own) The_Cpu_Arena()->free(mem->ptr);
    } else if (mem->type == SUNMEMTYPE_UVM) {
      if (mem->own) The_Managed_Arena()->free(mem->ptr);
    } else if (mem->type == SUNMEMTYPE_DEVICE) {
      if (mem->own) The_Device_Arena()->free(mem->ptr);
    } else if (mem->type == SUNMEMTYPE_PINNED) {
      if (mem->own) The_Pinned_Arena()->free(mem->ptr);
    } else {
      return(-1);
    }

    free(mem);
    return(0);
  }

  SUNMemoryHelper CloneMemoryHelper(SUNMemoryHelper)
  {
    return *The_SUNMemory_Helper();
  }

  int DestroyMemoryHelper(SUNMemoryHelper)
  {
    return(0);
  }

  SUNMemoryHelper CreateMemoryHelper()
  {
    SUNMemoryHelper helper;

    helper = SUNMemoryHelper_NewEmpty();
    
    helper->content        = NULL;
    helper->ops->clone     = CloneMemoryHelper;
    helper->ops->destroy   = DestroyMemoryHelper;
    helper->ops->alloc     = Alloc;
    helper->ops->dealloc   = Dealloc;
#if defined(AMREX_USE_HIP)
    helper->ops->copy      = SUNMemoryHelper_Copy_Hip;
    helper->ops->copyasync = SUNMemoryHelper_CopyAsync_Hip;
#elif defined(AMREX_USE_CUDA)
    helper->ops->copy      = SUNMemoryHelper_Copy_Cuda;
    helper->ops->copyasync = SUNMemoryHelper_CopyAsync_Cuda;
#endif

    return helper;
  }

  bool initialized = false;
  MemoryHelper* the_sunmemory_helper = nullptr;
}

MemoryHelper::MemoryHelper()
  : helper(CreateMemoryHelper())
{}

MemoryHelper::~MemoryHelper()
{
  if (helper->ops) free(helper->ops);
  free(helper);
}

void Initialize()
{
  if (initialized) return;
  initialized = true;

  BL_ASSERT(the_sunmemory_helper == nullptr);
  the_sunmemory_helper = new MemoryHelper();
}

void Finalize()
{
  initialized = false;
  delete the_sunmemory_helper;
  the_sunmemory_helper = nullptr;
}

MemoryHelper* The_SUNMemory_Helper()
{
  BL_ASSERT(the_sunmemory_helper != nullptr);
  return the_sunmemory_helper;
}

}
}
