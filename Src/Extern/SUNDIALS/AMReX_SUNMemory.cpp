#include <AMReX.H>
#include <AMReX_Gpu.H>
#include <AMReX_SUNMemory.H>
#if defined(AMREX_USE_HIP)
#include <sunmemory/sunmemory_hip.h>
#elif defined(AMREX_USE_CUDA)
#include <sunmemory/sunmemory_cuda.h>
#elif defined(AMREX_USE_DPCPP)
#include <sunmemory/sunmemory_sycl.h>
#endif

namespace amrex {
namespace sundials {

namespace {
  amrex::Arena* getArena (SUNMemoryType mem_type)
  {
    if (mem_type == SUNMEMTYPE_HOST) {
      return The_Cpu_Arena();
    } else if (mem_type == SUNMEMTYPE_UVM) {
        if (The_Arena()->isManaged()) {
          return The_Arena();
        } else if (The_Managed_Arena()->isManaged()) {
          return The_Managed_Arena();
        } else {
          return nullptr;
        }
    } else if (mem_type == SUNMEMTYPE_DEVICE) {
      return The_Device_Arena();
    } else if (mem_type == SUNMEMTYPE_PINNED) {
      return The_Pinned_Arena();
    } else {
      return nullptr;
    }
  }

  int Alloc(SUNMemoryHelper, SUNMemory* memptr, size_t memsize, SUNMemoryType mem_type)
  {
    SUNMemory mem = SUNMemoryNewEmpty();

    if (mem == nullptr) return -1;
    mem->ptr = NULL;
    mem->own = SUNTRUE;
    mem->type = mem_type;
    auto arena = getArena(mem->type);
    if (arena) {
      mem->ptr = arena->alloc(memsize);
      *memptr = mem;
      return 0;
    }
    else {
      free(mem);
      memptr = nullptr;
      return -1;
    }

    return -1;
  }

  int Dealloc(SUNMemoryHelper, SUNMemory mem)
  {

    if (mem == nullptr) return 0;
    auto arena = getArena(mem->type);
    if (arena) {
      if(mem->own)
      {
        arena->free(mem->ptr);
        free(mem);
        return 0;
      }
    }
    else {
      free(mem);
      return -1;
    }

    free(mem);
    return 0;
  }

  SUNMemoryHelper CloneMemoryHelper(SUNMemoryHelper)
  {
    return *The_SUNMemory_Helper();
  }

  int DestroyMemoryHelper(SUNMemoryHelper)
  {
    // We just return because we do not want our
    // single memory helper instance to be destroyed
    // except when Finalize is called.
    return 0;
  }

  SUNMemoryHelper CreateMemoryHelper()
  {
    SUNMemoryHelper helper;

    helper = SUNMemoryHelper_NewEmpty();

    helper->content        = NULL;
    helper->ops->clone     = CloneMemoryHelper;
    helper->ops->alloc     = Alloc;
    helper->ops->dealloc   = Dealloc;
    helper->ops->destroy   = DestroyMemoryHelper;
#if defined(AMREX_USE_HIP)
    helper->ops->copy      = SUNMemoryHelper_Copy_Hip;
    helper->ops->copyasync = SUNMemoryHelper_CopyAsync_Hip;
#elif defined(AMREX_USE_CUDA)
    helper->ops->copy      = SUNMemoryHelper_Copy_Cuda;
    helper->ops->copyasync = SUNMemoryHelper_CopyAsync_Cuda;
#elif defined(AMREX_USE_DPCPP)
    helper->ops->copy      = SUNMemoryHelper_Copy_Sycl;
    helper->ops->copyasync = SUNMemoryHelper_CopyAsync_Sycl;
    helper->ops->clone     = SUNMemoryHelper_Clone_Sycl;

    // Attach the queue pointer as the content
    helper->content = (void*) &amrex::Gpu::Device::streamQueue();
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
