
#include <AMReX_FabAllocator.H>

#ifdef AMREX_USE_GPU

namespace amrex {

namespace {
    constexpr std::size_t nblocks = 682;
}

FabPoolAllocator* fab_pool_allocator = nullptr;

constexpr std::size_t FabPoolAllocator::block_size;

FabPoolAllocator::FabPoolAllocator ()
#ifdef AMREX_FAB_IS_PINNED
    : m_arena(The_Pinned_Arena())
#else
    : m_arena(The_Managed_Arena())
#endif    
{
    void* p = m_arena->alloc(block_size*nblocks);
    m_orig.push_back(p);
    for (std::size_t i = 0; i < nblocks; ++i) {
        m_pool.push_back(p);
        p = (void*)((char*)p + block_size);
    }
}

FabPoolAllocator::~FabPoolAllocator ()
{
    for (auto p : m_orig) {
        m_arena->free(p);
    }
}

void*
FabPoolAllocator::alloc (std::size_t nbytes)
{
    AMREX_ASSERT(nbytes <= block_size);
    void* p = nullptr;
    std::lock_guard<std::mutex> lock(m_mutex);
    if (!m_pool.empty()) {
        p = m_pool.back();
        m_pool.pop_back();
    } else {
        p = m_arena->alloc(block_size*nblocks);
        m_orig.push_back(p);
        for (std::size_t i = 1; i < nblocks; ++i) {
            m_pool.push_back(p);
            p = (void*)((char*)p + block_size);
        }
    }
    return p;
}

void
FabPoolAllocator::free (void* ptr)
{
    std::lock_guard<std::mutex> lock(m_mutex);
    m_pool.push_back(ptr);
}

void makeFabPoolAllocator ()
{
    fab_pool_allocator = new FabPoolAllocator();
}

void destroyFabPoolAllocator ()
{
    delete fab_pool_allocator;
    fab_pool_allocator = nullptr;
}

}

#endif
