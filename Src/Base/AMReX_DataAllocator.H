#ifndef AMREX_DATA_ALLOCATOR_H_
#define AMREX_DATA_ALLOCATOR_H_
#include <AMReX_Config.H>

#include <AMReX_Arena.H>

namespace amrex {

struct DataAllocator {
    Arena* m_arena = nullptr;

    DataAllocator ()  noexcept  = default;

    DataAllocator (Arena* ar) noexcept : m_arena(ar) {}

    [[nodiscard]] void* alloc (std::size_t sz) const noexcept {
        return arena()->alloc(sz);
    }

    void free (void* pt) const noexcept {
        arena()->free(pt);
    }

    [[nodiscard]] Arena* arena () const noexcept {
        return (m_arena) ? m_arena : The_Arena();
    }
};

struct DataDeleter {
    Arena* m_arena = nullptr;
    DataDeleter () noexcept  = default;
    DataDeleter (Arena* ar) noexcept : m_arena(ar) {}
    [[nodiscard]] Arena* arena () const noexcept { return (m_arena) ? m_arena : The_Arena(); }
    void operator() (void* pt) const noexcept { arena()->free(pt); }
};

}

#endif
