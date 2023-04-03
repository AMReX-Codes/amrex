#include <AMReX.H>
#include <AMReX_Gpu.H>
#include <AMReX_SUNMemory.H>
#if defined(AMREX_USE_HIP)
#include <sunmemory/sunmemory_hip.h>
#elif defined(AMREX_USE_CUDA)
#include <sunmemory/sunmemory_cuda.h>
#elif defined(AMREX_USE_SYCL)
#include <sunmemory/sunmemory_sycl.h>
#endif

namespace amrex::sundials {

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

    int Alloc(SUNMemoryHelper, SUNMemory* memptr, size_t memsize, SUNMemoryType mem_type, void* /*queue*/)
    {
        SUNMemory mem = SUNMemoryNewEmpty();

        if (mem == nullptr) return -1;
        mem->ptr = nullptr;
        mem->own = SUNTRUE;
        mem->type = mem_type;
        auto* arena = getArena(mem->type);
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

    int Dealloc(SUNMemoryHelper, SUNMemory mem, void* /*queue*/)
    {

        if (mem == nullptr) return 0;
        auto* arena = getArena(mem->type);
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
        // memory helper instance to be destroyed
        // except when Finalize is called.
        return 0;
    }

    void ActuallyDestroySUNMemoryHelper(SUNMemoryHelper helper)
    {
        if (helper->ops) free(helper->ops);
        free(helper);
    }

    SUNMemoryHelper CreateMemoryHelper(::sundials::Context* sunctx)
    {
        SUNMemoryHelper helper;

        helper = SUNMemoryHelper_NewEmpty(*sunctx);

        helper->content          = nullptr;
        helper->ops->clone       = CloneMemoryHelper;
        helper->ops->alloc       = Alloc;
        helper->ops->dealloc     = Dealloc;
        helper->ops->destroy     = DestroyMemoryHelper;
#if defined(AMREX_USE_HIP)
        helper->ops->copy      = SUNMemoryHelper_Copy_Hip;
        helper->ops->copyasync = SUNMemoryHelper_CopyAsync_Hip;
#elif defined(AMREX_USE_CUDA)
        helper->ops->copy      = SUNMemoryHelper_Copy_Cuda;
        helper->ops->copyasync = SUNMemoryHelper_CopyAsync_Cuda;
#elif defined(AMREX_USE_SYCL)
        helper->ops->copy      = SUNMemoryHelper_Copy_Sycl;
        helper->ops->copyasync = SUNMemoryHelper_CopyAsync_Sycl;
#endif

        return helper;
    }

    Vector<int> initialized;
    Vector<MemoryHelper*> the_sunmemory_helper;
} //namespace

MemoryHelper::MemoryHelper(::sundials::Context* a_sunctx)
    : helper(CreateMemoryHelper(a_sunctx)),
        sunctx(a_sunctx)
{
}

MemoryHelper::~MemoryHelper()
{
    ActuallyDestroySUNMemoryHelper(helper);
}

MemoryHelper::MemoryHelper(const MemoryHelper& rhs)
    : helper(CreateMemoryHelper(rhs.sunctx)),
        sunctx(rhs.sunctx)
{}

MemoryHelper::MemoryHelper(MemoryHelper&& rhs) noexcept
    : helper(rhs.helper),
        sunctx(rhs.sunctx)
{
    rhs.helper = nullptr;
    rhs.sunctx = nullptr;
}

MemoryHelper& MemoryHelper::operator=(MemoryHelper&& rhs) noexcept
{
    if (this != &rhs)
    {
        ActuallyDestroySUNMemoryHelper(helper);
        helper = rhs.helper;
        rhs.helper = nullptr;
        delete sunctx;
        sunctx = rhs.sunctx;
        rhs.sunctx = nullptr;
    }
    return *this;
}

void MemoryHelper::Initialize(int nthreads)
{
    if (initialized.empty()) {
        initialized.resize(nthreads);
        std::fill(initialized.begin(), initialized.end(), 0);
        the_sunmemory_helper.resize(nthreads);
        std::fill(the_sunmemory_helper.begin(), the_sunmemory_helper.end(), nullptr);
    }
    for (int i = 0; i < nthreads; i++) {
        if (initialized[i]) continue;
        initialized[i] = 1;
        BL_ASSERT(the_sunmemory_helper[i] == nullptr);
        the_sunmemory_helper[i] = new MemoryHelper(The_Sundials_Context(i));
    }
}

void MemoryHelper::Finalize()
{
    for (int i = 0; i < initialized.size(); i++) {
        initialized[i] = 0;
        delete the_sunmemory_helper[i];
        the_sunmemory_helper[i] = nullptr;
    }
}

MemoryHelper* The_SUNMemory_Helper(int i)
{
    BL_ASSERT(the_sunmemory_helper[i] != nullptr);
    return the_sunmemory_helper[i];
}

}
