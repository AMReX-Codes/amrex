#ifndef AMREX_FAB_ALLOCATOR_H_
#define AMREX_FAB_ALLOCATOR_H_

#include <vector>
#include <mutex>
#include <AMReX_Arena.H>
#include <AMReX_TypeTraits.H>

namespace amrex {

#ifndef AMREX_USE_GPU

struct FabAllocator {};

template <class T, class Enable = void> struct IsFabAllocatorSafe : std::true_type {};

#else

struct FabPoolAllocator
{
public:
    FabPoolAllocator ();
    ~FabPoolAllocator ();
    void* alloc (std::size_t nbytes);
    void free (void* ptr);
    static constexpr std::size_t block_size = 96;
private:
    std::vector<void*> m_pool;
    std::vector<void*> m_orig;
    std::mutex m_mutex;
    Arena* m_arena = nullptr;
};

extern FabPoolAllocator* fab_pool_allocator;
void makeFabPoolAllocator ();
void destroyFabPoolAllocator ();

struct FabAllocator {
    void* operator new (std::size_t nbytes) { return fab_pool_allocator->alloc(nbytes); }
    void operator delete (void* ptr) { fab_pool_allocator->free(ptr); }
};



template <class T, class Enable = void> struct IsFabAllocatorSafe : std::false_type {};

template <class T> struct IsFabAllocatorSafe<T, typename std::enable_if<
        !amrex::IsBaseFab<T>::value || sizeof(T) <= FabPoolAllocator::block_size
                                                                       >::type
                                            >
    : std::true_type {};

#endif

}

#endif
